//! JPL SPK ephemeris reader (DAF container, segment types 2 and 3).
//!
//! Spec: NAIF "DAF Required Reading" and "SPK Required Reading" (public NAIF documentation).
//! Files are memory-mapped, so opening is cheap and lookups never block on I/O after the
//! pages are resident. [`Ephemeris`] is `Send + Sync` and is meant to be shared via `Arc`.

use std::collections::HashMap;
use std::fs::File;
use std::path::Path;

use memmap2::Mmap;

use crate::error::EphemError;
use crate::math::{Vec3, add, scale, sub};
use crate::models::DeltaTModel;

/// Seconds per day.
const DAY_S: f64 = 86_400.0;
/// J2000.0 as a Julian Day (TDB).
const J2000_JD: f64 = 2_451_545.0;
/// Astronomical unit in km (IAU 2012, as used by DE440).
pub const AU_KM: f64 = 149_597_870.7;

const RECORD_BYTES: usize = 1024;
/// Largest modified-difference-array table size supported (SPK types 1 and 21).
const MAX_MDA_DIM: usize = 25;

#[derive(Debug, Clone, Copy)]
struct Segment {
    file: usize,
    target: i32,
    center: i32,
    data_type: i32,
    start_et: f64,
    end_et: f64,
    /// Word (f64) index of the first data word, 0-based.
    start_word: usize,
    /// Directory values (type 2/3).
    init: f64,
    intlen: f64,
    rsize: usize,
    n_records: usize,
}

struct SpkFile {
    map: Mmap,
}

impl SpkFile {
    #[inline]
    fn word(&self, idx: usize) -> f64 {
        let b = &self.map[idx * 8..idx * 8 + 8];
        f64::from_le_bytes(b.try_into().expect("8 bytes"))
    }

    #[inline]
    fn i32_at(&self, byte: usize) -> i32 {
        i32::from_le_bytes(self.map[byte..byte + 4].try_into().expect("4 bytes"))
    }
}

/// One or more SPK files indexed by target body.
pub struct Ephemeris {
    files: Vec<SpkFile>,
    delta_t_model: DeltaTModel,
    /// Segments by target; later files take precedence (searched first).
    by_target: HashMap<i32, Vec<Segment>>,
}

impl std::fmt::Debug for Ephemeris {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Ephemeris")
            .field("files", &self.files.len())
            .field("targets", &self.by_target.len())
            .finish()
    }
}

impl Ephemeris {
    /// Opens a single SPK file.
    pub fn open(path: impl AsRef<Path>) -> Result<Self, EphemError> {
        let mut e = Ephemeris {
            files: Vec::new(),
            delta_t_model: DeltaTModel::default(),
            by_target: HashMap::new(),
        };
        e.add_file(path)?;
        Ok(e)
    }

    /// Opens every `.bsp` file in a directory: the planetary file(s) plus any small-body files
    /// (for example those fetched by `tools/fetch_horizons_spk.py`).
    pub fn open_dir(dir: impl AsRef<Path>) -> Result<Self, EphemError> {
        let dir = dir.as_ref();
        let io = |e: std::io::Error| EphemError::Io(dir.display().to_string(), e.to_string());
        let mut paths: Vec<_> = std::fs::read_dir(dir)
            .map_err(io)?
            .filter_map(|e| e.ok().map(|e| e.path()))
            .filter(|p| p.extension().is_some_and(|x| x.eq_ignore_ascii_case("bsp")))
            .collect();
        paths.sort();
        let mut e = Ephemeris { files: Vec::new(), delta_t_model: DeltaTModel::default(), by_target: HashMap::new() };
        for p in &paths {
            e.add_file(p)?;
        }
        if e.files.is_empty() {
            return Err(EphemError::Io(dir.display().to_string(), "no .bsp files found".into()));
        }
        Ok(e)
    }

    /// Adds another SPK file; its segments take precedence over earlier files.
    pub fn add_file(&mut self, path: impl AsRef<Path>) -> Result<(), EphemError> {
        let path = path.as_ref();
        let file = File::open(path).map_err(|e| EphemError::Io(path.display().to_string(), e.to_string()))?;
        // SAFETY: the file is opened read-only and treated as immutable for the lifetime of the
        // map. Modifying an SPK file while it is mapped is unsupported.
        let map = unsafe { Mmap::map(&file) }
            .map_err(|e| EphemError::Io(path.display().to_string(), e.to_string()))?;
        let spk = SpkFile { map };
        let file_idx = self.files.len();
        let segments = read_segments(&spk, file_idx, path)?;
        self.files.push(spk);
        for seg in segments {
            self.by_target.entry(seg.target).or_default().insert(0, seg);
        }
        Ok(())
    }

    /// Chooses the ΔT model used to convert UT to TT (default: observed values).
    pub fn with_delta_t_model(mut self, model: DeltaTModel) -> Self {
        self.delta_t_model = model;
        self
    }

    pub fn delta_t_model(&self) -> DeltaTModel {
        self.delta_t_model
    }

    /// NAIF ids of all targets available.
    pub fn targets(&self) -> Vec<i32> {
        let mut t: Vec<i32> = self.by_target.keys().copied().collect();
        t.sort_unstable();
        t
    }

    /// Time span (TDB Julian Days) covered for a target, merged across its segments.
    pub fn coverage(&self, target: i32) -> Option<(f64, f64)> {
        let segs = self.by_target.get(&target)?;
        let start = segs.iter().map(|s| s.start_et).fold(f64::INFINITY, f64::min);
        let end = segs.iter().map(|s| s.end_et).fold(f64::NEG_INFINITY, f64::max);
        Some((J2000_JD + start / DAY_S, J2000_JD + end / DAY_S))
    }

    /// Position (AU) and velocity (AU/day) of `target` relative to its segment centre.
    fn segment_state(&self, target: i32, tdb_jd: f64) -> Result<(i32, Vec3, Vec3), EphemError> {
        let et = (tdb_jd - J2000_JD) * DAY_S;
        let segs = self
            .by_target
            .get(&target)
            .ok_or(EphemError::NoData { target, jd: tdb_jd })?;
        let seg = segs
            .iter()
            .find(|s| et >= s.start_et && et <= s.end_et)
            .ok_or(EphemError::NoData { target, jd: tdb_jd })?;
        let (p, v) = match seg.data_type {
            1 | 21 => eval_mda(&self.files[seg.file], seg, et),
            _ => eval_segment(&self.files[seg.file], seg, et),
        };
        Ok((seg.center, scale(p, 1.0 / AU_KM), scale(v, DAY_S / AU_KM)))
    }

    /// Position (AU) and velocity (AU/day) of `target` relative to the solar-system
    /// barycentre (NAIF 0), in the ICRF, at a TDB Julian Day.
    pub fn barycentric(&self, target: i32, tdb_jd: f64) -> Result<(Vec3, Vec3), EphemError> {
        let mut pos = [0.0; 3];
        let mut vel = [0.0; 3];
        let mut body = target;
        let mut hops = 0;
        while body != 0 {
            let (center, p, v) = self.segment_state(body, tdb_jd)?;
            pos = add(pos, p);
            vel = add(vel, v);
            body = center;
            hops += 1;
            if hops > 8 {
                return Err(EphemError::NoData { target, jd: tdb_jd });
            }
        }
        Ok((pos, vel))
    }

    /// State of `target` relative to `center` (both NAIF ids), AU and AU/day.
    pub fn relative(&self, target: i32, center: i32, tdb_jd: f64) -> Result<(Vec3, Vec3), EphemError> {
        let (pt, vt) = self.barycentric(target, tdb_jd)?;
        let (pc, vc) = self.barycentric(center, tdb_jd)?;
        Ok((sub(pt, pc), sub(vt, vc)))
    }
}

fn read_segments(spk: &SpkFile, file_idx: usize, path: &Path) -> Result<Vec<Segment>, EphemError> {
    let bad = |why: &str| EphemError::Format(path.display().to_string(), why.to_string());
    if spk.map.len() < RECORD_BYTES {
        return Err(bad("file shorter than one record"));
    }
    if &spk.map[0..7] != b"DAF/SPK" {
        return Err(bad("not a DAF/SPK file"));
    }
    if &spk.map[88..96] != b"LTL-IEEE" {
        return Err(bad("only little-endian IEEE SPK files are supported"));
    }
    let nd = spk.i32_at(8) as usize;
    let ni = spk.i32_at(12) as usize;
    if nd != 2 || ni != 6 {
        return Err(bad("unexpected ND/NI for SPK"));
    }
    let summary_words = nd + ni.div_ceil(2);
    let mut segments = Vec::new();
    let mut record = spk.i32_at(76) as usize; // FWARD
    while record != 0 {
        let base_word = (record - 1) * RECORD_BYTES / 8;
        if (base_word + 128) * 8 > spk.map.len() {
            return Err(bad("summary record out of range"));
        }
        let next = spk.word(base_word) as usize;
        let nsum = spk.word(base_word + 2) as usize;
        for i in 0..nsum {
            let w = base_word + 3 + i * summary_words;
            let start_et = spk.word(w);
            let end_et = spk.word(w + 1);
            let ib = (w + 2) * 8;
            let target = spk.i32_at(ib);
            let center = spk.i32_at(ib + 4);
            let _frame = spk.i32_at(ib + 8);
            let data_type = spk.i32_at(ib + 12);
            let start_addr = spk.i32_at(ib + 16) as usize;
            let end_addr = spk.i32_at(ib + 20) as usize;
            let end_word = end_addr - 1;
            let (init, intlen, rsize, n_records) = match data_type {
                2 | 3 => {
                    let rsize = spk.word(end_word - 1) as usize;
                    let n_records = spk.word(end_word) as usize;
                    let sets = if data_type == 2 { 3 } else { 6 };
                    if rsize < 2 + sets || (rsize - 2) / sets > 32 || n_records == 0 {
                        return Err(bad("unsupported segment record layout"));
                    }
                    (spk.word(end_word - 3), spk.word(end_word - 2), rsize, n_records)
                }
                // Modified difference arrays: type 1 has a fixed table size of 15, type 21 stores it.
                1 | 21 => {
                    let n_records = spk.word(end_word) as usize;
                    let maxdim = if data_type == 21 { spk.word(end_word - 1) as usize } else { 15 };
                    if maxdim == 0 || maxdim > MAX_MDA_DIM || n_records == 0 {
                        return Err(bad("unsupported difference-array segment"));
                    }
                    (0.0, 0.0, 4 * maxdim + 11, n_records)
                }
                _ => continue, // other types are not needed
            };
            segments.push(Segment {
                file: file_idx,
                target,
                center,
                data_type,
                start_et,
                end_et,
                start_word: start_addr - 1,
                init,
                intlen,
                rsize,
                n_records,
            });
        }
        record = next;
    }
    Ok(segments)
}

/// Evaluates a type 2/3 segment: km and km/s.
#[allow(clippy::needless_range_loop)]
fn eval_segment(f: &SpkFile, seg: &Segment, et: f64) -> (Vec3, Vec3) {
    let idx = (((et - seg.init) / seg.intlen).floor() as isize).clamp(0, seg.n_records as isize - 1)
        as usize;
    let rec = seg.start_word + idx * seg.rsize;
    let mid = f.word(rec);
    let radius = f.word(rec + 1);
    let s = (et - mid) / radius;
    let sets = if seg.data_type == 2 { 3 } else { 6 };
    let ncoef = (seg.rsize - 2) / sets;

    // Chebyshev polynomials and their derivatives at s.
    let mut t = [0.0f64; 32];
    let mut dt = [0.0f64; 32];
    let n = ncoef.min(32);
    t[0] = 1.0;
    if n > 1 {
        t[1] = s;
        dt[1] = 1.0;
    }
    for k in 2..n {
        t[k] = 2.0 * s * t[k - 1] - t[k - 2];
        dt[k] = 2.0 * t[k - 1] + 2.0 * s * dt[k - 1] - dt[k - 2];
    }

    let mut pos = [0.0; 3];
    let mut vel = [0.0; 3];
    for c in 0..3 {
        let base = rec + 2 + c * ncoef;
        let mut p = 0.0;
        let mut v = 0.0;
        // sum smallest terms first for accuracy
        for k in (0..n).rev() {
            let coef = f.word(base + k);
            p += coef * t[k];
            v += coef * dt[k];
        }
        pos[c] = p;
        vel[c] = v / radius;
    }
    if seg.data_type == 3 {
        for c in 0..3 {
            let base = rec + 2 + (3 + c) * ncoef;
            let mut v = 0.0;
            for k in (0..n).rev() {
                v += f.word(base + k) * t[k];
            }
            vel[c] = v;
        }
    }
    (pos, vel)
}

/// Evaluates a type 1/21 segment (modified difference arrays): km and km/s.
///
/// Algorithm: NAIF SPK Required Reading, types 1 and 21 (the SPKE01/SPKE21 recurrences);
/// cross-checked against the MIT-licensed `spktype21` Python implementation.
fn eval_mda(f: &SpkFile, seg: &Segment, et: f64) -> (Vec3, Vec3) {
    let maxdim = (seg.rsize - 11) / 4;
    let n = seg.n_records;
    // Final epochs of the records follow the records; pick the first record whose epoch ≥ et.
    let epochs = seg.start_word + n * seg.rsize;
    let (mut lo, mut hi) = (0usize, n - 1);
    while lo < hi {
        let mid = (lo + hi) / 2;
        if f.word(epochs + mid) < et {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    let rec = seg.start_word + lo * seg.rsize;
    let w_at = |i: usize| f.word(rec + i);

    let tl = w_at(0);
    let g = |j: usize| w_at(1 + j);
    let refpos = [w_at(maxdim + 1), w_at(maxdim + 3), w_at(maxdim + 5)];
    let refvel = [w_at(maxdim + 2), w_at(maxdim + 4), w_at(maxdim + 6)];
    let dt = |j: usize, i: usize| w_at(maxdim + 7 + i * maxdim + j);
    let kqmax1 = w_at(4 * maxdim + 7) as usize;
    let kq = [w_at(4 * maxdim + 8) as usize, w_at(4 * maxdim + 9) as usize, w_at(4 * maxdim + 10) as usize];

    let delta = et - tl;
    let mut fc = [0.0f64; MAX_MDA_DIM + 2];
    let mut wc = [0.0f64; MAX_MDA_DIM + 2];
    let mut w = [0.0f64; MAX_MDA_DIM + 4];
    fc[0] = 1.0;
    let mut tp = delta;
    let mq2 = kqmax1.saturating_sub(2);
    let mut ks = kqmax1 - 1;
    for j in 1..=mq2 {
        fc[j] = tp / g(j - 1);
        wc[j - 1] = delta / g(j - 1);
        tp = delta + g(j - 1);
    }
    for (j, wj) in w.iter_mut().enumerate().take(kqmax1) {
        *wj = 1.0 / (j + 1) as f64;
    }
    let mut jx = 0usize;
    let mut ks1 = ks - 1;
    while ks >= 2 {
        jx += 1;
        for j in 1..=jx {
            w[j + ks - 1] = fc[j] * w[j + ks1 - 1] - wc[j - 1] * w[j + ks - 1];
        }
        ks = ks1;
        ks1 -= 1;
    }
    let mut pos = [0.0; 3];
    let mut vel = [0.0; 3];
    for i in 0..3 {
        let mut sum = 0.0;
        for j in (1..=kq[i]).rev() {
            sum += dt(j - 1, i) * w[j + ks - 1];
        }
        pos[i] = refpos[i] + delta * (refvel[i] + delta * sum);
    }
    for j in 1..=jx {
        w[j + ks - 1] = fc[j] * w[j + ks1 - 1] - wc[j - 1] * w[j + ks - 1];
    }
    ks -= 1;
    for i in 0..3 {
        let mut sum = 0.0;
        for j in (1..=kq[i]).rev() {
            sum += dt(j - 1, i) * w[j + ks - 1];
        }
        vel[i] = refvel[i] + delta * sum;
    }
    (pos, vel)
}
