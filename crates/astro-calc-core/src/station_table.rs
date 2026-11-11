//! A pre-computed table of planetary stations.
//!
//! Finding stations live costs milliseconds; looking them up in a table costs nanoseconds.
//! [`StationTable::build`] computes a table from whatever ephemeris files are loaded, and
//! the table converts to and from a compact binary form that can live in a file, be embedded
//! with `include_bytes!`, or be fetched by a browser.
//!
//! Binary layout (all little-endian):
//!
//! | bytes | content |
//! |---|---|
//! | 8 | magic `ACSTAT01` |
//! | 4 + 4 | format version, number of bodies |
//! | 8 + 8 | first and last UT Julian Day searched |
//! | 8 | ΔT model code (1 byte) and padding |
//! | 8 | ephemeris fingerprint |
//! | 32 per body | body key (24 bytes, zero-padded), first record index, record count |
//! | 16 per station | Julian Day (f64), speed (f32), kind (1 byte), padding |

use std::path::Path;

use crate::body::Body;
use crate::error::{CalcError, EphemError};
use crate::events::{JdSpeed, Station, StationKind, StationSet};
use crate::models::DeltaTModel;
use crate::position::naif_id;
use crate::spk::Ephemeris;
use crate::time::JdUt;

const MAGIC: &[u8; 8] = b"ACSTAT01";
const VERSION: u32 = 1;
const HEADER_BYTES: usize = 48;
const ENTRY_BYTES: usize = 32;
const KEY_BYTES: usize = 24;
const RECORD_BYTES: usize = 16;
/// Days kept clear of each end of the ephemeris, where positions cannot be computed.
const EDGE_DAYS: f64 = 2.0;
const NAIF_EARTH: i32 = 399;

const fn model_code(model: DeltaTModel) -> u8 {
    match model {
        DeltaTModel::Observed => 0,
        DeltaTModel::AstroCalcApi => 1,
    }
}

const fn kind_code(kind: StationKind) -> u8 {
    match kind {
        StationKind::RetroStart => 0,
        StationKind::RetroPeak => 1,
        StationKind::RetroEnd => 2,
        StationKind::Peak => 3,
    }
}

/// Identifies the ephemeris a table was built from: Earth's barycentric distance at J2000,
/// which differs between ephemeris versions.
fn fingerprint(eph: &Ephemeris) -> f64 {
    eph.barycentric(NAIF_EARTH, crate::time::J2000).map(|(p, _)| crate::math::norm(p)).unwrap_or(0.0)
}

/// Stations of a set of bodies over a span of time, sorted for fast lookup.
#[derive(Debug, Clone, PartialEq)]
pub struct StationTable {
    from: f64,
    to: f64,
    delta_t_model: DeltaTModel,
    fingerprint: f64,
    bodies: Vec<(Body, Vec<Station>)>,
}

impl StationTable {
    /// Computes the stations of `bodies` between two times from the loaded ephemeris.
    /// A body that never turns retrograde gets an empty list.
    pub fn build(eph: &Ephemeris, bodies: &[Body], from: impl Into<JdUt>, to: impl Into<JdUt>) -> Result<Self, CalcError> {
        let (from, to) = (from.into().value(), to.into().value());
        if !(to > from && from.is_finite() && to.is_finite()) {
            return Err(CalcError::InvalidInterval { days: to - from, interval: 0.0 });
        }
        let mut list = Vec::with_capacity(bodies.len());
        for &body in bodies {
            if list.iter().any(|(b, _)| *b == body) {
                continue;
            }
            list.push((body, eph.stations_in(body, from, to, true)?));
        }
        Ok(Self { from, to, delta_t_model: eph.delta_t_model(), fingerprint: fingerprint(eph), bodies: list })
    }

    /// Computes the stations of `bodies` over the whole span the ephemeris covers for all
    /// of them (less two days at each end, where positions cannot be computed).
    pub fn build_all(eph: &Ephemeris, bodies: &[Body]) -> Result<Self, CalcError> {
        let (mut from, mut to) = (f64::NEG_INFINITY, f64::INFINITY);
        for id in bodies.iter().filter_map(|&b| naif_id(b)).chain([NAIF_EARTH]) {
            let (a, b) = eph.coverage(id).ok_or(EphemError::NoData { target: id, jd: crate::time::J2000 })?;
            from = from.max(a);
            to = to.min(b);
        }
        Self::build(eph, bodies, from + EDGE_DAYS, to - EDGE_DAYS)
    }

    /// The span searched, as UT Julian Days.
    pub fn range(&self) -> (f64, f64) {
        (self.from, self.to)
    }

    /// The ΔT model in force when the table was built.
    pub fn delta_t_model(&self) -> DeltaTModel {
        self.delta_t_model
    }

    pub fn bodies(&self) -> impl Iterator<Item = Body> + '_ {
        self.bodies.iter().map(|(b, _)| *b)
    }

    /// Total number of stations held.
    pub fn len(&self) -> usize {
        self.bodies.iter().map(|(_, s)| s.len()).sum()
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Whether the table was built from this ephemeris with its current ΔT model. A table
    /// that does not match should be rebuilt.
    pub fn matches(&self, eph: &Ephemeris) -> bool {
        self.delta_t_model == eph.delta_t_model() && (self.fingerprint - fingerprint(eph)).abs() < 1e-12
    }

    /// Whether the table can answer for this body at this time.
    pub fn covers(&self, body: Body, jd: f64) -> bool {
        jd >= self.from && jd <= self.to && self.stations(body).is_some()
    }

    /// Every station of a body, in chronological order; `None` if the body is not in the table.
    pub fn stations(&self, body: Body) -> Option<&[Station]> {
        self.bodies.iter().find(|(b, _)| *b == body).map(|(_, s)| s.as_slice())
    }

    /// The last `before` stations at or before `jd` followed by the next `after` stations
    /// after it, as one chronological slice. Shorter where the table ends. `None` if the
    /// body is absent or `jd` lies outside the table.
    pub fn around(&self, body: Body, jd: f64, before: usize, after: usize) -> Option<&[Station]> {
        if !(jd >= self.from && jd <= self.to) {
            return None;
        }
        let all = self.stations(body)?;
        let at = all.partition_point(|s| s.jd() <= jd);
        Some(&all[at.saturating_sub(before)..(at + after).min(all.len())])
    }

    /// The next `count` stations after `jd`.
    pub fn next(&self, body: Body, jd: f64, count: usize) -> Option<&[Station]> {
        self.around(body, jd, 0, count)
    }

    /// Stations with `from <= jd <= to`.
    pub fn between(&self, body: Body, from: f64, to: f64) -> Option<&[Station]> {
        let all = self.stations(body)?;
        let (a, b) = (all.partition_point(|s| s.jd() < from), all.partition_point(|s| s.jd() <= to));
        Some(&all[a..b.max(a)])
    }

    /// Whether the body is retrograde at `jd`, from the table alone: it is if the last
    /// turning point was a retrograde start. `None` if the table cannot answer.
    pub fn is_retrograde(&self, body: Body, jd: f64) -> Option<bool> {
        if !(jd >= self.from && jd <= self.to) {
            return None;
        }
        let all = self.stations(body)?;
        let at = all.partition_point(|s| s.jd() <= jd);
        let turning = |s: &&Station| matches!(s.kind(), StationKind::RetroStart | StationKind::RetroEnd);
        if let Some(last) = all[..at].iter().rev().find(turning) {
            return Some(last.kind() == StationKind::RetroStart);
        }
        // Before the first turning point: retrograde if the next one ends a retrograde.
        Some(all[at..].iter().find(turning).is_some_and(|s| s.kind() == StationKind::RetroEnd))
    }

    /// The body's current speed with its next `count` stations and the `stations_before`
    /// preceding ones: the same result as [`Ephemeris::next_stations`], taken from the table
    /// where it covers the request and computed live otherwise. The speed is always live.
    pub fn station_set(&self, eph: &Ephemeris, body: Body, jd: impl Into<JdUt>, count: usize, stations_before: u8) -> Result<StationSet, CalcError> {
        let jd = jd.into().value();
        match self.around(body, jd, stations_before as usize, count) {
            Some(found) => Ok(StationSet { body, jd, speed: eph.lng_speed(body, jd)?.1, stations: found.to_vec() }),
            None => eph.next_stations(body, jd, count, stations_before),
        }
    }

    /// The table in its binary form.
    pub fn to_bytes(&self) -> Result<Vec<u8>, EphemError> {
        let total = self.len();
        let mut out = Vec::with_capacity(HEADER_BYTES + self.bodies.len() * ENTRY_BYTES + total * RECORD_BYTES);
        out.extend_from_slice(MAGIC);
        out.extend_from_slice(&VERSION.to_le_bytes());
        out.extend_from_slice(&(self.bodies.len() as u32).to_le_bytes());
        out.extend_from_slice(&self.from.to_le_bytes());
        out.extend_from_slice(&self.to.to_le_bytes());
        out.extend_from_slice(&[model_code(self.delta_t_model), 0, 0, 0, 0, 0, 0, 0]);
        out.extend_from_slice(&self.fingerprint.to_le_bytes());
        let mut first = 0u32;
        for (body, stations) in &self.bodies {
            let key = body.to_string();
            if key.len() > KEY_BYTES {
                return Err(EphemError::Table(format!("body key too long: {key}")));
            }
            let mut name = [0u8; KEY_BYTES];
            name[..key.len()].copy_from_slice(key.as_bytes());
            out.extend_from_slice(&name);
            out.extend_from_slice(&first.to_le_bytes());
            out.extend_from_slice(&(stations.len() as u32).to_le_bytes());
            first += stations.len() as u32;
        }
        for (_, stations) in &self.bodies {
            for s in stations {
                out.extend_from_slice(&s.jd().to_le_bytes());
                out.extend_from_slice(&(s.speed() as f32).to_le_bytes());
                out.extend_from_slice(&[kind_code(s.kind()), 0, 0, 0]);
            }
        }
        Ok(out)
    }

    /// Reads a table from its binary form.
    pub fn from_bytes(bytes: &[u8]) -> Result<Self, EphemError> {
        let bad = |why: &str| EphemError::Table(why.to_string());
        let u32_at = |i: usize| u32::from_le_bytes(bytes[i..i + 4].try_into().expect("4 bytes"));
        let f64_at = |i: usize| f64::from_le_bytes(bytes[i..i + 8].try_into().expect("8 bytes"));
        if bytes.len() < HEADER_BYTES || &bytes[..8] != MAGIC {
            return Err(bad("not a station table"));
        }
        if u32_at(8) != VERSION {
            return Err(bad("unsupported format version"));
        }
        let n_bodies = u32_at(12) as usize;
        let (from, to) = (f64_at(16), f64_at(24));
        let delta_t_model = match bytes[32] {
            0 => DeltaTModel::Observed,
            1 => DeltaTModel::AstroCalcApi,
            _ => return Err(bad("unknown delta-T model")),
        };
        let fingerprint = f64_at(40);
        let records_at = n_bodies.checked_mul(ENTRY_BYTES).and_then(|n| n.checked_add(HEADER_BYTES)).ok_or_else(|| bad("bad body count"))?;
        if bytes.len() < records_at || !(bytes.len() - records_at).is_multiple_of(RECORD_BYTES) {
            return Err(bad("truncated"));
        }
        let n_records = (bytes.len() - records_at) / RECORD_BYTES;
        let mut bodies = Vec::with_capacity(n_bodies);
        for i in 0..n_bodies {
            let e = HEADER_BYTES + i * ENTRY_BYTES;
            let name = &bytes[e..e + KEY_BYTES];
            let key = std::str::from_utf8(&name[..name.iter().position(|&c| c == 0).unwrap_or(KEY_BYTES)]).map_err(|_| bad("bad body key"))?;
            let body: Body = key.parse().map_err(|_| bad("unknown body key"))?;
            let (first, count) = (u32_at(e + KEY_BYTES) as usize, u32_at(e + KEY_BYTES + 4) as usize);
            if first.checked_add(count).is_none_or(|end| end > n_records) {
                return Err(bad("body records out of range"));
            }
            let mut stations = Vec::with_capacity(count);
            for r in first..first + count {
                let at = records_at + r * RECORD_BYTES;
                let jd = f64_at(at);
                let speed = f32::from_le_bytes(bytes[at + 8..at + 12].try_into().expect("4 bytes")) as f64;
                if stations.last().is_some_and(|prev: &Station| prev.jd() >= jd) || jd.is_nan() {
                    return Err(bad("stations out of order"));
                }
                stations.push(match bytes[at + 12] {
                    0 => Station::RetroStart(jd),
                    1 => Station::RetroPeak(JdSpeed(jd, speed)),
                    2 => Station::RetroEnd(jd),
                    3 => Station::Peak(JdSpeed(jd, speed)),
                    _ => return Err(bad("unknown station kind")),
                });
            }
            bodies.push((body, stations));
        }
        Ok(Self { from, to, delta_t_model, fingerprint, bodies })
    }

    /// Writes the binary form to a file.
    pub fn save(&self, path: impl AsRef<Path>) -> Result<(), EphemError> {
        let path = path.as_ref();
        std::fs::write(path, self.to_bytes()?).map_err(|e| EphemError::Io(path.display().to_string(), e.to_string()))
    }

    /// Reads a table from a file.
    pub fn open(path: impl AsRef<Path>) -> Result<Self, EphemError> {
        let path = path.as_ref();
        let bytes = std::fs::read(path).map_err(|e| EphemError::Io(path.display().to_string(), e.to_string()))?;
        Self::from_bytes(&bytes)
    }
}
