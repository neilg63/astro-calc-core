//! Excerpts of SPK files: a span of a long ephemeris (DE441) as a smaller, valid `.bsp`,
//! fetched with HTTP range requests (only the bytes needed) or cut from a local copy.
//!
//! Type 2 and 3 segments hold fixed-length Chebyshev records, `INIT + i * INTLEN` apart, with
//! a four-word directory (INIT, INTLEN, RSIZE, N) at the end. An excerpt copies the records
//! covering the span byte for byte and writes a new directory, summaries and file record, so
//! positions inside the span are exactly those of the source file.
//!
//! DAF layout written: record 1 the file record (the source's, with new pointers), then a
//! summary record and a name record for every 25 segments, then the segments' data.

use std::fs::{self, File, OpenOptions};
use std::io::{BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};

use anyhow::{Context, Result, bail};
use astro_calc_core::{DateRange, Ephemeris};
use indicatif::{ProgressBar, ProgressStyle};
use ureq::Agent;

const REC: u64 = 1024;
const J2000: f64 = 2_451_545.0;
const DAY_S: f64 = 86_400.0;
/// Summaries per summary record for ND = 2, NI = 6: (128 - 3) / 5.
const PER_RECORD: usize = 25;

/// Where the full file is read from.
pub enum Source {
    Http { agent: Agent, url: String },
    File(PathBuf),
}

impl Source {
    pub fn describe(&self) -> String {
        match self {
            Source::Http { url, .. } => format!("range requests to {url}"),
            Source::File(p) => format!("cut from {}", p.display()),
        }
    }

    fn reader(&self, offset: u64, len: u64) -> Result<Box<dyn Read + '_>> {
        match self {
            Source::Http { agent, url } => {
                let resp = agent.get(url).header("Range", format!("bytes={offset}-{}", offset + len - 1)).call()?;
                if resp.status().as_u16() != 206 {
                    bail!("{url}: the server does not serve byte ranges (HTTP {}); download the whole file and use --from-file", resp.status().as_u16());
                }
                Ok(Box::new(resp.into_body().into_reader().take(len)))
            }
            Source::File(p) => {
                let mut f = File::open(p).with_context(|| format!("cannot open {}", p.display()))?;
                f.seek(SeekFrom::Start(offset))?;
                Ok(Box::new(BufReader::with_capacity(1 << 20, f).take(len)))
            }
        }
    }

    fn read(&self, offset: u64, len: u64) -> Result<Vec<u8>> {
        let mut buf = Vec::with_capacity(len as usize);
        self.reader(offset, len)?.read_to_end(&mut buf)?;
        if buf.len() as u64 != len {
            bail!("short read at byte {offset}: {} of {len} bytes", buf.len());
        }
        Ok(buf)
    }
}

fn f64_at(b: &[u8], o: usize) -> f64 {
    f64::from_le_bytes(b[o..o + 8].try_into().expect("8 bytes"))
}

fn i32_at(b: &[u8], o: usize) -> i32 {
    i32::from_le_bytes(b[o..o + 4].try_into().expect("4 bytes"))
}

/// One segment of the excerpt: where its records come from and what it becomes.
#[derive(Debug, Clone)]
pub struct Segment {
    /// target, centre, frame, type (begin and end addresses are recomputed).
    ints: [i32; 4],
    name: Vec<u8>,
    intlen: f64,
    rsize: u64,
    /// Byte offset in the source of the first record copied.
    src_offset: u64,
    /// New values.
    init: f64,
    n: u64,
    start_et: f64,
    end_et: f64,
}

impl Segment {
    fn data_bytes(&self) -> u64 {
        self.n * self.rsize * 8
    }

    fn out_bytes(&self) -> u64 {
        self.data_bytes() + 32
    }

    fn directory(&self) -> Vec<u8> {
        [self.init, self.intlen, self.rsize as f64, self.n as f64].iter().flat_map(|v| v.to_le_bytes()).collect()
    }
}

/// Everything needed to write an excerpt, worked out from the source's index alone.
pub struct Plan {
    pub header: Vec<u8>,
    pub segments: Vec<Segment>,
    /// Size of the excerpt.
    pub total: u64,
    /// The span shared by every body in it.
    pub range: DateRange,
}

/// Reads the source's index (a few small reads) and plans the excerpt for `from`..`to`
/// (Julian Days). Records are kept whole, so the excerpt reaches a little beyond the span.
pub fn plan(src: &Source, from_jd: f64, to_jd: f64, label: &str) -> Result<Plan> {
    let (a_et, b_et) = ((from_jd - J2000) * DAY_S, (to_jd - J2000) * DAY_S);
    let file_record = src.read(0, REC)?;
    if &file_record[0..7] != b"DAF/SPK" || &file_record[88..96] != b"LTL-IEEE" {
        bail!("not a little-endian SPK file");
    }
    if (i32_at(&file_record, 8), i32_at(&file_record, 12)) != (2, 6) {
        bail!("unexpected ND/NI for an SPK file");
    }
    // The source's summaries and names, following the record chain.
    let mut found = Vec::new();
    let mut record = i32_at(&file_record, 76) as u64;
    for _ in 0..1000 {
        if record == 0 {
            break;
        }
        let pair = src.read((record - 1) * REC, 2 * REC)?;
        let (summ, names) = pair.split_at(REC as usize);
        let next = f64_at(summ, 0) as u64;
        for i in 0..(f64_at(summ, 16) as usize).min(PER_RECORD) {
            let w = (3 + i * 5) * 8;
            let ints: Vec<i32> = (0..6).map(|k| i32_at(summ, w + 16 + 4 * k)).collect();
            found.push((f64_at(summ, w), f64_at(summ, w + 8), ints, names[i * 40..i * 40 + 40].to_vec()));
        }
        record = next;
    }
    let mut segments = Vec::new();
    for (start, end, ints, name) in found {
        let (target, data_type, begin, end_addr) = (ints[0], ints[3], ints[4] as u64, ints[5] as u64);
        let a = a_et.max(start);
        let b = b_et.min(end);
        if b < a {
            continue; // outside the span
        }
        if !matches!(data_type, 2 | 3) {
            bail!("segment for NAIF {target} is type {data_type}; only types 2 and 3 can be excerpted");
        }
        let dir = src.read((end_addr - 4) * 8, 32)?;
        let (init, intlen, rsize, n) = (f64_at(&dir, 0), f64_at(&dir, 8), f64_at(&dir, 16) as u64, f64_at(&dir, 24) as u64);
        if intlen <= 0.0 || rsize < 2 || n == 0 || (begin - 1) + n * rsize + 4 != end_addr {
            bail!("segment for NAIF {target} has an unexpected layout");
        }
        let last = n as f64 - 1.0;
        let i0 = ((a - init) / intlen).floor().clamp(0.0, last) as u64;
        let i1 = ((b - init) / intlen).floor().clamp(0.0, last) as u64;
        let count = i1 - i0 + 1;
        let new_init = init + i0 as f64 * intlen;
        segments.push(Segment {
            ints: [ints[0], ints[1], ints[2], ints[3]],
            name,
            intlen,
            rsize,
            src_offset: ((begin - 1) + i0 * rsize) * 8,
            init: new_init,
            n: count,
            start_et: start.max(new_init),
            end_et: end.min(new_init + count as f64 * intlen),
        });
    }
    if segments.is_empty() {
        bail!("the source holds nothing between the dates asked for");
    }

    // Shared span: per target the union of its segments, then the overlap of all targets.
    let mut targets: Vec<i32> = segments.iter().map(|s| s.ints[0]).collect();
    targets.sort_unstable();
    targets.dedup();
    let (mut start, mut end) = (f64::NEG_INFINITY, f64::INFINITY);
    for t in targets {
        let own = segments.iter().filter(|s| s.ints[0] == t);
        start = start.max(own.clone().map(|s| s.start_et).fold(f64::INFINITY, f64::min));
        end = end.min(own.map(|s| s.end_et).fold(f64::NEG_INFINITY, f64::max));
    }
    let range = DateRange { start: J2000 + start / DAY_S, end: J2000 + end / DAY_S };

    // Header: file record, then a summary and a name record per 25 segments.
    let chunks = segments.len().div_ceil(PER_RECORD);
    let header_records = 1 + 2 * chunks as u64;
    let mut word = header_records * REC / 8 + 1; // first data word (1-based)
    let mut addresses = Vec::new();
    for s in &segments {
        let words = s.out_bytes() / 8;
        addresses.push((word, word + words - 1));
        word += words;
    }
    let mut header = file_record.clone();
    let mut ifn = format!("astro-calc-build excerpt of {label}").into_bytes();
    ifn.resize(60, b' ');
    header[16..76].copy_from_slice(&ifn);
    header[76..80].copy_from_slice(&2i32.to_le_bytes()); // FWARD
    header[80..84].copy_from_slice(&((2 + 2 * (chunks as i32 - 1)).to_le_bytes())); // BWARD
    header[84..88].copy_from_slice(&(i32::try_from(word).context("excerpt too large for a DAF file")?.to_le_bytes())); // FREE
    for c in 0..chunks {
        let this = 2 + 2 * c;
        let mut summ = vec![0u8; REC as usize];
        let mut names = vec![b' '; REC as usize];
        let next = if c + 1 < chunks { this + 2 } else { 0 };
        let prev = if c > 0 { this - 2 } else { 0 };
        let part = &segments[c * PER_RECORD..segments.len().min((c + 1) * PER_RECORD)];
        summ[0..8].copy_from_slice(&(next as f64).to_le_bytes());
        summ[8..16].copy_from_slice(&(prev as f64).to_le_bytes());
        summ[16..24].copy_from_slice(&(part.len() as f64).to_le_bytes());
        for (i, s) in part.iter().enumerate() {
            let (b, e) = addresses[c * PER_RECORD + i];
            let w = (3 + i * 5) * 8;
            summ[w..w + 8].copy_from_slice(&s.start_et.to_le_bytes());
            summ[w + 8..w + 16].copy_from_slice(&s.end_et.to_le_bytes());
            for (k, v) in s.ints.iter().copied().chain([b as i32, e as i32]).enumerate() {
                summ[w + 16 + 4 * k..w + 20 + 4 * k].copy_from_slice(&v.to_le_bytes());
            }
            names[i * 40..i * 40 + 40].copy_from_slice(&s.name);
        }
        header.extend_from_slice(&summ);
        header.extend_from_slice(&names);
    }
    let total = (word - 1) * 8;
    Ok(Plan { header, segments, total, range })
}

/// Writes the excerpt to `dest`'s part file, resuming a previous partial run.
pub fn write(src: &Source, plan: &Plan, dest: &Path) -> Result<PathBuf> {
    let part = crate::net::part_path(dest);
    let mut have = fs::metadata(&part).map(|m| m.len()).unwrap_or(0);
    // A part from a different plan (other span or source) starts again.
    if have > 0 {
        let mut head = vec![0u8; plan.header.len().min(have as usize)];
        File::open(&part)?.read_exact(&mut head)?;
        if have < plan.header.len() as u64 || head != plan.header || have > plan.total {
            have = 0;
        }
    }
    let mut out = if have > 0 {
        eprintln!("  resuming {} at {}", dest.file_name().unwrap_or_default().to_string_lossy(), crate::catalog::fmt_bytes(have as f64));
        BufWriter::with_capacity(1 << 20, OpenOptions::new().append(true).open(&part)?)
    } else {
        let mut w = BufWriter::with_capacity(1 << 20, File::create(&part)?);
        w.write_all(&plan.header)?;
        w
    };
    let pb = ProgressBar::new(plan.total).with_style(
        ProgressStyle::with_template("  {msg:16} [{bar:30}] {bytes}/{total_bytes} {binary_bytes_per_sec} eta {eta}")
            .expect("template")
            .progress_chars("=> "),
    );
    pb.set_message(dest.file_name().unwrap_or_default().to_string_lossy().into_owned());
    let mut pos = have.max(plan.header.len() as u64);
    pb.set_position(pos);
    let mut skip = pos - plan.header.len() as u64;
    let result = (|| -> Result<()> {
        for s in &plan.segments {
            if skip >= s.out_bytes() {
                skip -= s.out_bytes();
                continue;
            }
            if skip < s.data_bytes() {
                let len = s.data_bytes() - skip;
                let mut r = pb.wrap_read(src.reader(s.src_offset + skip, len)?);
                let copied = std::io::copy(&mut r, &mut out)?;
                if copied != len {
                    bail!("short read: {copied} of {len} bytes");
                }
                pos += len;
                skip = s.data_bytes();
            }
            let dir = s.directory();
            out.write_all(&dir[(skip - s.data_bytes()) as usize..])?;
            pos += 32 - (skip - s.data_bytes());
            pb.set_position(pos);
            skip = 0;
        }
        Ok(())
    })();
    out.flush()?;
    pb.finish_and_clear();
    result.context("excerpt interrupted; run the same command again to resume")?;
    let got = fs::metadata(&part)?.len();
    if got != plan.total {
        bail!("excerpt has {got} bytes, expected {}", plan.total);
    }
    Ok(part)
}

/// Checks that every record sits where its directory says: each record's own midpoint and
/// radius must match INIT + (i + ½) INTLEN and INTLEN / 2. A byte offset wrong by even one
/// word would fail this. Then that the file opens and covers the span.
pub fn check(path: &Path, plan: &Plan, from_jd: f64, to_jd: f64) -> Result<()> {
    let mut f = BufReader::with_capacity(1 << 20, File::open(path)?);
    f.seek(SeekFrom::Start(plan.header.len() as u64))?;
    let mut rec = Vec::new();
    for s in &plan.segments {
        rec.resize(s.rsize as usize * 8, 0);
        for i in 0..s.n {
            f.read_exact(&mut rec)?;
            let (mid, radius) = (f64_at(&rec, 0), f64_at(&rec, 8));
            let want = s.init + (i as f64 + 0.5) * s.intlen;
            if (mid - want).abs() > 1e-6 * s.intlen || (radius - s.intlen / 2.0).abs() > 1e-6 * s.intlen {
                bail!("record {i} of NAIF {} is misplaced (midpoint {mid}, expected {want})", s.ints[0]);
            }
        }
        f.seek_relative(32)?;
    }
    let eph = Ephemeris::open(path).context("the excerpt does not open as an ephemeris")?;
    let r = eph.files().first().map(|f| f.range).context("empty excerpt")?;
    if r.start > from_jd.max(plan.range.start) + 1.0 || r.end < to_jd.min(plan.range.end) - 1.0 {
        bail!("the excerpt covers {r}, not the span asked for");
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::catalog::date_jd;

    /// An excerpt of a local DE440 file must give bit-identical positions to the full file.
    #[test]
    fn excerpt_matches_source() {
        let data = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../data");
        let Some(full) = ["de440s.bsp", "de440.bsp"].iter().map(|f| data.join(f)).find(|p| p.is_file()) else {
            eprintln!("skipped: no de440s.bsp or de440.bsp in {}", data.display());
            return;
        };
        let dir = std::env::temp_dir().join(format!("acb-excerpt-{}", std::process::id()));
        fs::create_dir_all(&dir).unwrap();
        let dest = dir.join("cut.bsp");
        let (from, to) = (date_jd(1900, 1, 1), date_jd(1950, 1, 1));
        let src = Source::File(full.clone());
        let p = plan(&src, from, to, "de440s").unwrap();
        let part = write(&src, &p, &dest).unwrap();
        check(&part, &p, from, to).unwrap();
        fs::rename(&part, &dest).unwrap();
        assert!(fs::metadata(&dest).unwrap().len() < fs::metadata(&full).unwrap().len() / 4);

        let (a, b) = (Ephemeris::open(&full).unwrap(), Ephemeris::open(&dest).unwrap());
        let r = b.date_range().unwrap();
        assert!(r.start <= from && r.end >= to, "{r}");
        for k in 0..400 {
            let jd = from + k as f64 * (to - from) / 400.0 + 0.123;
            for id in [10, 301, 399, 1, 2, 4, 5, 6, 7, 8, 9] {
                assert_eq!(a.barycentric(id, jd).unwrap(), b.barycentric(id, jd).unwrap(), "NAIF {id} at JD {jd}");
            }
        }
        // Outside the span the excerpt has no data.
        assert!(b.barycentric(399, date_jd(1870, 1, 1)).is_err());

        // Resuming a truncated part gives the same file.
        let bytes = fs::read(&dest).unwrap();
        fs::write(crate::net::part_path(&dest), &bytes[..bytes.len() * 3 / 5]).unwrap();
        let part = write(&src, &p, &dest).unwrap();
        assert_eq!(fs::read(&part).unwrap(), bytes);
        fs::remove_dir_all(&dir).unwrap();
    }
}
