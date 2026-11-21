//! Planetary stations and lunar phases, found by root- and extremum-finding on apparent
//! geocentric ecliptic positions.

use crate::body::Body;
use crate::spk::{Calc, CalcEngine};
use crate::error::CalcError;
use crate::geo::GeoPos;
use crate::math::diff_deg;
use crate::position::{CoordSystem, Perspective};
use crate::spk::Ephemeris;
use crate::time::JdUt;

/// Points of a planet's retrograde cycle.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum StationKind {
    /// Direct motion stops; retrograde begins (speed falls through zero).
    RetroStart,
    /// Fastest retrograde motion (most negative speed).
    RetroPeak,
    /// Retrograde motion stops; direct resumes (speed rises through zero).
    RetroEnd,
    /// Fastest direct motion.
    Peak,
}

impl StationKind {
    /// Key used by the AstroCalc API.
    pub const fn key(self) -> &'static str {
        match self {
            StationKind::RetroStart => "retro-start",
            StationKind::RetroPeak => "retro-peak",
            StationKind::RetroEnd => "retro-end",
            StationKind::Peak => "peak",
        }
    }
}

/// A moment (UT Julian Day) and the longitude speed at it (degrees per day).
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct JdSpeed(pub f64, pub f64);

impl JdSpeed {
    pub const fn jd(&self) -> f64 {
        self.0
    }

    pub const fn speed(&self) -> f64 {
        self.1
    }
}

/// One station of a body's retrograde cycle. The two turning points carry only their time,
/// since the speed there is zero; the two peaks carry the speed reached.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum Station {
    /// Fastest direct motion.
    Peak(JdSpeed),
    /// Retrograde motion begins.
    RetroStart(f64),
    /// Fastest retrograde motion (a negative speed).
    RetroPeak(JdSpeed),
    /// Retrograde motion ends.
    RetroEnd(f64),
}

impl Station {
    /// When the station occurs (UT Julian Day).
    pub const fn jd(&self) -> f64 {
        match self {
            Station::Peak(p) | Station::RetroPeak(p) => p.0,
            Station::RetroStart(jd) | Station::RetroEnd(jd) => *jd,
        }
    }

    /// Longitude speed at the station in degrees per day: zero at the turning points.
    pub const fn speed(&self) -> f64 {
        match self {
            Station::Peak(p) | Station::RetroPeak(p) => p.1,
            Station::RetroStart(_) | Station::RetroEnd(_) => 0.0,
        }
    }

    pub const fn kind(&self) -> StationKind {
        match self {
            Station::Peak(_) => StationKind::Peak,
            Station::RetroStart(_) => StationKind::RetroStart,
            Station::RetroPeak(_) => StationKind::RetroPeak,
            Station::RetroEnd(_) => StationKind::RetroEnd,
        }
    }
}

/// A body's motion at a moment, with the stations around it.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct StationSet {
    pub body: Body,
    /// The moment asked about (UT Julian Day).
    pub jd: f64,
    /// Longitude speed at `jd` in degrees per day; negative when retrograde.
    pub speed: f64,
    /// Stations in chronological order: any requested ones before `jd`, then those after.
    pub stations: Vec<Station>,
    /// The engine that produced the speed and stations.
    #[cfg_attr(feature = "serde", serde(default))]
    pub engine: CalcEngine,
}

impl StationSet {
    /// Whether the body is moving backwards at `jd`.
    pub fn is_retrograde(&self) -> bool {
        self.speed < 0.0
    }

    /// Stations after `jd`.
    pub fn upcoming(&self) -> impl Iterator<Item = &Station> {
        self.stations.iter().filter(|s| s.jd() > self.jd)
    }

    /// Stations at or before `jd`.
    pub fn previous(&self) -> impl Iterator<Item = &Station> {
        self.stations.iter().filter(|s| s.jd() <= self.jd)
    }

    /// The next station of a given kind after `jd`, if within the set.
    pub fn next(&self, kind: StationKind) -> Option<&Station> {
        self.upcoming().find(|s| s.kind() == kind)
    }
}

/// The four principal lunar phases.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum MoonPhase {
    New,
    FirstQuarter,
    Full,
    LastQuarter,
}

impl MoonPhase {
    /// 1 = new, 2 = first quarter, 3 = full, 4 = last quarter.
    pub const fn num(self) -> u8 {
        self as u8 + 1
    }

    /// Moon − Sun elongation in degrees.
    pub const fn angle(self) -> f64 {
        self as u8 as f64 * 90.0
    }

    pub const fn is_waxing(self) -> bool {
        matches!(self, MoonPhase::New | MoonPhase::FirstQuarter)
    }

    const fn from_quadrant(q: i64) -> Self {
        match q.rem_euclid(4) {
            0 => MoonPhase::New,
            1 => MoonPhase::FirstQuarter,
            2 => MoonPhase::Full,
            _ => MoonPhase::LastQuarter,
        }
    }
}

/// A lunar phase and the moment it is exact.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct MoonPhaseEvent {
    pub phase: MoonPhase,
    pub jd: JdUt,
    /// The engine that produced the time.
    #[cfg_attr(feature = "serde", serde(default))]
    pub engine: CalcEngine,
}

/// The two halves of the lunar month.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum Paksha {
    /// The bright half: new moon to full moon, the Moon waxing.
    Shukla,
    /// The dark half: full moon to new moon, the Moon waning.
    Krishna,
}

impl Paksha {
    pub const fn name(self) -> &'static str {
        match self {
            Paksha::Shukla => "Shukla",
            Paksha::Krishna => "Krishna",
        }
    }
}

/// A tithi, the lunar day of the Indian calendar: one of thirty equal 12° steps of the
/// Moon − Sun elongation, numbered 1 to 30 from the new moon. It lasts from about 19 to 26
/// hours. The elongation is a difference of longitudes, so no ayanamsha is involved.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Tithi(u8);

impl Tithi {
    /// Degrees of elongation in one tithi.
    pub const SPAN: f64 = 12.0;

    const NAMES: [&'static str; 14] = [
        "Pratipada", "Dwitiya", "Tritiya", "Chaturthi", "Panchami", "Shashthi", "Saptami", "Ashtami", "Navami", "Dashami", "Ekadashi", "Dwadashi",
        "Trayodashi", "Chaturdashi",
    ];

    /// The tithi with this number, 1 to 30.
    pub const fn new(num: u8) -> Option<Self> {
        if num >= 1 && num <= 30 { Some(Self(num)) } else { None }
    }

    /// The tithi in force at a Moon − Sun elongation in degrees.
    pub fn from_elongation(elongation: f64) -> Self {
        Self(((elongation.rem_euclid(360.0) / Self::SPAN) as u8).min(29) + 1)
    }

    /// 1 to 30, counted from the new moon.
    pub const fn num(self) -> u8 {
        self.0
    }

    pub const fn paksha(self) -> Paksha {
        if self.0 <= 15 { Paksha::Shukla } else { Paksha::Krishna }
    }

    /// The day within its paksha, 1 to 15.
    pub const fn day(self) -> u8 {
        (self.0 - 1) % 15 + 1
    }

    /// The name without the paksha: Pratipada to Chaturdashi, then Purnima (the 15th, ending
    /// at the full moon) or Amavasya (the 30th, ending at the new moon).
    pub const fn name(self) -> &'static str {
        match self.0 {
            15 => "Purnima",
            30 => "Amavasya",
            _ => Self::NAMES[(self.day() - 1) as usize],
        }
    }

    /// The elongation at which it begins, in degrees.
    pub const fn start_angle(self) -> f64 {
        (self.0 - 1) as f64 * Self::SPAN
    }

    /// The elongation at which it ends, in degrees (360 for the 30th).
    pub const fn end_angle(self) -> f64 {
        self.0 as f64 * Self::SPAN
    }

    /// The principal phase at which it begins: the new moon for the 1st, the full moon for
    /// the 16th. The quarters fall in the middle of the 8th and the 23rd, not at a boundary.
    pub const fn phase_at_start(self) -> Option<MoonPhase> {
        match self.0 {
            1 => Some(MoonPhase::New),
            16 => Some(MoonPhase::Full),
            _ => None,
        }
    }
}

/// The beginning of a tithi.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct TithiEvent {
    /// The tithi that begins.
    pub tithi: Tithi,
    pub jd: JdUt,
    /// The engine that produced the time.
    #[cfg_attr(feature = "serde", serde(default))]
    pub engine: CalcEngine,
}

/// The eleven karana names. Seven of them recur eight times in a lunar month; four occur once.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum KaranaName {
    Bava,
    Balava,
    Kaulava,
    Taitila,
    Gara,
    Vanija,
    Vishti,
    Shakuni,
    Chatushpada,
    Naga,
    Kimstughna,
}

impl KaranaName {
    pub const fn name(self) -> &'static str {
        match self {
            KaranaName::Bava => "Bava",
            KaranaName::Balava => "Balava",
            KaranaName::Kaulava => "Kaulava",
            KaranaName::Taitila => "Taitila",
            KaranaName::Gara => "Gara",
            KaranaName::Vanija => "Vanija",
            KaranaName::Vishti => "Vishti",
            KaranaName::Shakuni => "Shakuni",
            KaranaName::Chatushpada => "Chatushpada",
            KaranaName::Naga => "Naga",
            KaranaName::Kimstughna => "Kimstughna",
        }
    }

    /// Whether it is one of the seven recurring (movable) karanas.
    pub const fn is_recurring(self) -> bool {
        (self as u8) < 7
    }
}

/// A karana, half a tithi: one of sixty equal 6° steps of the Moon − Sun elongation, numbered
/// 1 to 60 from the new moon. The four principal phases all fall on karana boundaries.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Karana(u8);

impl Karana {
    /// Degrees of elongation in one karana.
    pub const SPAN: f64 = 6.0;

    /// The karana with this number, 1 to 60.
    pub const fn new(num: u8) -> Option<Self> {
        if num >= 1 && num <= 60 { Some(Self(num)) } else { None }
    }

    /// The karana in force at a Moon − Sun elongation in degrees.
    pub fn from_elongation(elongation: f64) -> Self {
        Self(((elongation.rem_euclid(360.0) / Self::SPAN) as u8).min(59) + 1)
    }

    /// 1 to 60, counted from the new moon.
    pub const fn num(self) -> u8 {
        self.0
    }

    /// The tithi it belongs to.
    pub const fn tithi(self) -> Tithi {
        Tithi((self.0 - 1) / 2 + 1)
    }

    /// Whether it is the first half of its tithi.
    pub const fn is_first_half(self) -> bool {
        self.0 % 2 == 1
    }

    /// Its name: Kimstughna for the 1st, the seven recurring names in turn for the 2nd to the
    /// 57th, then Shakuni, Chatushpada and Naga.
    pub const fn name(self) -> KaranaName {
        match self.0 {
            1 => KaranaName::Kimstughna,
            58 => KaranaName::Shakuni,
            59 => KaranaName::Chatushpada,
            60 => KaranaName::Naga,
            n => match (n - 2) % 7 {
                0 => KaranaName::Bava,
                1 => KaranaName::Balava,
                2 => KaranaName::Kaulava,
                3 => KaranaName::Taitila,
                4 => KaranaName::Gara,
                5 => KaranaName::Vanija,
                _ => KaranaName::Vishti,
            },
        }
    }

    /// The elongation at which it begins, in degrees.
    pub const fn start_angle(self) -> f64 {
        (self.0 - 1) as f64 * Self::SPAN
    }

    /// The principal phase at which it begins, if any: the 1st at the new moon, the 16th at
    /// the first quarter, the 31st at the full moon and the 46th at the last quarter.
    pub const fn phase_at_start(self) -> Option<MoonPhase> {
        match self.0 {
            1 => Some(MoonPhase::New),
            16 => Some(MoonPhase::FirstQuarter),
            31 => Some(MoonPhase::Full),
            46 => Some(MoonPhase::LastQuarter),
            _ => None,
        }
    }
}

/// The beginning of a karana.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct KaranaEvent {
    /// The karana that begins.
    pub karana: Karana,
    pub jd: JdUt,
    /// The engine that produced the time.
    #[cfg_attr(feature = "serde", serde(default))]
    pub engine: CalcEngine,
}

/// Sampling step (days) for station searches: well under the shortest gap between a body's
/// consecutive stations.
fn station_step(body: Body) -> f64 {
    match body {
        Body::Mercury => 1.0,
        Body::Venus | Body::Mars => 2.0,
        _ => 4.0,
    }
}

/// Padding (days) longer than the body's longest retrograde cycle (four stations).
fn station_pad(body: Body) -> f64 {
    match body {
        Body::Mercury => 130.0,
        Body::Venus => 600.0,
        Body::Jupiter | Body::Saturn | Body::Uranus | Body::Neptune | Body::Pluto | Body::Chiron | Body::Pholus => 420.0,
        // Mars, the main-belt asteroids and anything unknown.
        _ => 800.0,
    }
}

const GOLD: f64 = 0.381_966_011_250_105;
/// Days in a sidereal month.
const SIDEREAL_MONTH: f64 = 27.321_66;

impl Ephemeris {
    pub(crate) fn lng_speed(&self, body: Body, t: f64) -> Result<(f64, f64), CalcError> {
        let p = self.calc_position(body, t, GeoPos::default(), CoordSystem::Ecliptic, Perspective::Geocentric)?;
        Ok((p.lng, p.lng_speed))
    }

    /// Stations of a body between two times, in chronological order.
    ///
    /// `RetroStart` / `RetroEnd` are the zero crossings of the apparent geocentric longitude
    /// speed. `RetroPeak` and `Peak` mark the fastest retrograde and direct motion of each
    /// period (one of each per cycle), placed at the top of the smooth speed curve: small
    /// monthly ripples can make the instantaneous speed a day or two away fractionally higher.
    pub fn stations(&self, body: Body, from: impl Into<JdUt>, to: impl Into<JdUt>) -> Result<Calc<Vec<Station>>, CalcError> {
        let from = from.into().value();
        Ok(self.calc(from, self.stations_in(body, from, to.into().value(), true)?))
    }

    /// A body's speed at `jd` with its stations over the following `days`, preceded by the
    /// last `stations_before` stations before `jd`.
    pub fn planetary_stations(&self, body: Body, jd: impl Into<JdUt>, days: f64, stations_before: u8) -> Result<StationSet, CalcError> {
        let jd = jd.into().value();
        if !(days >= 0.0 && days.is_finite()) {
            return Err(CalcError::InvalidInterval { days, interval: station_step(body) });
        }
        let speed = self.lng_speed(body, jd)?.1;
        let mut stations = self.stations_before(body, jd, stations_before as usize)?;
        stations.extend(self.stations_in(body, jd, jd + days, false)?.into_iter().filter(|s| s.jd() > jd));
        Ok(StationSet { body, jd, speed, stations, engine: self.engine_used(jd) })
    }

    /// A body's speed at `jd` with its next `count` stations, preceded by the last
    /// `stations_before` stations before `jd`. Fewer are returned for a body that never
    /// turns retrograde (the Sun and Moon) or where the ephemeris ends first.
    pub fn next_stations(&self, body: Body, jd: impl Into<JdUt>, count: usize, stations_before: u8) -> Result<StationSet, CalcError> {
        let jd = jd.into().value();
        let speed = self.lng_speed(body, jd)?.1;
        let mut stations = self.stations_before(body, jd, stations_before as usize)?;
        // A cycle has four stations and is shorter than the padding length.
        let span = (count / 4 + 1) as f64 * station_pad(body);
        let ahead = self.stations_in(body, jd, jd + span, false)?;
        stations.extend(ahead.into_iter().filter(|s| s.jd() > jd).take(count));
        Ok(StationSet { body, jd, speed, stations, engine: self.engine_used(jd) })
    }

    /// The last `count` stations at or before `jd`.
    fn stations_before(&self, body: Body, jd: f64, count: usize) -> Result<Vec<Station>, CalcError> {
        if count == 0 {
            return Ok(Vec::new());
        }
        let span = (count / 4 + 1) as f64 * station_pad(body);
        let mut before = self.stations_in(body, jd - span, jd, false)?;
        Ok(before.split_off(before.len().saturating_sub(count)))
    }

    /// `strict`: fail if any part of `from..to` lies outside the ephemeris; otherwise search
    /// only the part that is covered.
    pub(crate) fn stations_in(&self, body: Body, from: f64, to: f64, strict: bool) -> Result<Vec<Station>, CalcError> {
        let step = station_step(body);
        // Search a padded range so the direct/retrograde periods overlapping the edges are complete.
        let pad = station_pad(body);
        let speed = |t: f64| self.lng_speed(body, t).map(|x| x.1);
        let n = (((to - from) + 2.0 * pad) / step).ceil() as usize + 1;
        let t_at = |i: usize| from - pad + i as f64 * step;
        let mut samples = Vec::with_capacity(n);
        for i in 0..n {
            match speed(t_at(i)) {
                Ok(s) => samples.push(s),
                // the padding may run past the ephemeris; the requested range itself must not
                Err(e) if strict && t_at(i) >= from && t_at(i) <= to => return Err(e),
                Err(_) => samples.push(f64::NAN),
            }
        }

        // Zero crossings: (time, falling?) with the sample index just before each.
        let mut crossings: Vec<(f64, bool, usize)> = Vec::new();
        for i in 0..n - 1 {
            let (s1, s2) = (samples[i], samples[i + 1]);
            if s1.is_nan() || s2.is_nan() || (s1 > 0.0) == (s2 > 0.0) {
                continue;
            }
            let (mut a, mut b, falling) = (t_at(i), t_at(i + 1), s1 > 0.0);
            while b - a > 1e-7 {
                let m = 0.5 * (a + b);
                if (speed(m)? > 0.0) == falling {
                    a = m;
                } else {
                    b = m;
                }
            }
            crossings.push((0.5 * (a + b), falling, i));
        }

        let mut out = Vec::new();
        let mut push = |kind, t: f64| -> Result<(), CalcError> {
            if t >= from && t <= to {
                out.push(match kind {
                    StationKind::RetroStart => Station::RetroStart(t),
                    StationKind::RetroEnd => Station::RetroEnd(t),
                    StationKind::Peak => Station::Peak(JdSpeed(t, self.lng_speed(body, t)?.1)),
                    StationKind::RetroPeak => Station::RetroPeak(JdSpeed(t, self.lng_speed(body, t)?.1)),
                });
            }
            Ok(())
        };
        for (k, &(t, falling, _)) in crossings.iter().enumerate() {
            push(if falling { StationKind::RetroStart } else { StationKind::RetroEnd }, t)?;
            // Extreme speed of the period that follows this crossing.
            let Some(&(t_next, _, i_next)) = crossings.get(k + 1) else { continue };
            let sign = if falling { 1.0 } else { -1.0 }; // minimise speed when retrograde, −speed when direct
            let (i0, i1) = (crossings[k].2 + 1, i_next);
            let Some(best) = (i0..=i1).min_by(|&x, &y| (sign * samples[x]).total_cmp(&(sign * samples[y]))) else { continue };
            let (mut a, mut b) = ((t_at(best) - step).max(t), (t_at(best) + step).min(t_next));
            let mut x1 = a + GOLD * (b - a);
            let mut x2 = b - GOLD * (b - a);
            let (mut f1, mut f2) = (sign * speed(x1)?, sign * speed(x2)?);
            while b - a > 1e-5 {
                if f1 < f2 {
                    (b, x2, f2) = (x2, x1, f1);
                    x1 = a + GOLD * (b - a);
                    f1 = sign * speed(x1)?;
                } else {
                    (a, x1, f1) = (x1, x2, f2);
                    x2 = b - GOLD * (b - a);
                    f2 = sign * speed(x2)?;
                }
            }
            // Speed is almost flat at its extreme, and for all but Mercury it carries small
            // monthly ripples (the Earth circling the Earth–Moon barycentre, and nutation), so
            // the highest value is ambiguous by days. The peak is therefore defined as the
            // top of the smooth curve: the zero of the acceleration, taken as the second
            // difference of longitude over a sidereal month, which cancels those ripples.
            let mut peak = 0.5 * (a + b);
            let (h, reach) = if body == Body::Mercury { (step, step) } else { (SIDEREAL_MONTH.min((t_next - t) / 4.0), 10.0) };
            let slope = |x: f64| -> Result<f64, CalcError> {
                let (l0, l1, l2) = (self.lng_speed(body, x - h)?.0, self.lng_speed(body, x)?.0, self.lng_speed(body, x + h)?.0);
                Ok(diff_deg(l1, l2) - diff_deg(l0, l1))
            };
            let (mut lo, mut hi) = ((peak - reach).max(t), (peak + reach).min(t_next));
            // Near the ends of the ephemeris the interval may not fit; the estimate then stands.
            if hi > lo
                && let (Ok(mut s_lo), Ok(s_hi)) = (slope(lo), slope(hi))
                && (s_lo > 0.0) != (s_hi > 0.0)
            {
                let mut ok = true;
                while hi - lo > 1e-6 {
                    let m = 0.5 * (lo + hi);
                    let Ok(s_m) = slope(m) else {
                        ok = false;
                        break;
                    };
                    if (s_m > 0.0) == (s_lo > 0.0) {
                        (lo, s_lo) = (m, s_m);
                    } else {
                        hi = m;
                    }
                }
                if ok {
                    peak = 0.5 * (lo + hi);
                }
            }
            push(if falling { StationKind::RetroPeak } else { StationKind::Peak }, peak)?;
        }
        out.sort_by(|a, b| a.jd().total_cmp(&b.jd()));
        Ok(out)
    }

    fn elongation_at(&self, t: f64, observer: Option<GeoPos>) -> Result<(f64, f64), CalcError> {
        let (geo, persp) = match observer {
            Some(g) => (g, Perspective::Topocentric),
            None => (GeoPos::default(), Perspective::Geocentric),
        };
        let m = self.calc_position(Body::Moon, t, geo, CoordSystem::Ecliptic, persp)?;
        let s = self.calc_position(Body::Sun, t, geo, CoordSystem::Ecliptic, persp)?;
        Ok(((m.lng - s.lng).rem_euclid(360.0), m.lng_speed - s.lng_speed))
    }

    /// Geocentric Moon − Sun elongation in ecliptic longitude (degrees, 0–360) and its rate (°/day).
    pub fn moon_elongation(&self, jd: impl Into<JdUt>) -> Result<Calc<(f64, f64)>, CalcError> {
        let jd = jd.into().value();
        Ok(self.calc(jd, self.elongation_at(jd, None)?))
    }

    /// The first moment at or after `start` when the Moon − Sun elongation is a whole number
    /// of `span` degrees (which must divide 360). Returns that number, from 0 at the new moon
    /// up to `360 / span − 1`, and the time.
    fn next_elongation_step(&self, start: f64, span: f64, observer: Option<GeoPos>) -> Result<(u32, f64), CalcError> {
        let (e, _) = self.elongation_at(start, observer)?;
        let step = (e / span).ceil();
        let target = (step * span).rem_euclid(360.0);
        // first estimate from the mean synodic rate, then Newton steps with the actual rate
        let mut t = start + (step * span - e) / 12.190_749;
        for _ in 0..20 {
            let (e, rate) = self.elongation_at(t, observer)?;
            let dt = diff_deg(e, target) / rate;
            t += dt;
            if dt.abs() < 1e-9 {
                break;
            }
        }
        Ok(((step as u32) % (360.0 / span).round() as u32, t.max(start)))
    }

    /// Every such moment between two times, in order.
    fn elongation_steps_between(&self, from: f64, to: f64, span: f64, observer: Option<GeoPos>) -> Result<Vec<(u32, f64)>, CalcError> {
        // The elongation never gains more than about 15° a day, so the next step is at least
        // span / 15 days on: resume the search a little short of that.
        let skip = span / 20.0;
        let mut t = from;
        let mut out = Vec::new();
        loop {
            let (step, jd) = self.next_elongation_step(t, span, observer)?;
            if jd > to {
                return Ok(out);
            }
            out.push((step, jd));
            t = jd + skip;
        }
    }

    fn next_phase_at(&self, start: f64, observer: Option<GeoPos>) -> Result<MoonPhaseEvent, CalcError> {
        let (quadrant, jd) = self.next_elongation_step(start, 90.0, observer)?;
        Ok(MoonPhaseEvent { phase: MoonPhase::from_quadrant(quadrant as i64), jd: JdUt(jd), engine: self.engine_used(jd) })
    }

    fn phases_between(&self, from: f64, to: f64, observer: Option<GeoPos>) -> Result<Vec<MoonPhaseEvent>, CalcError> {
        let steps = self.elongation_steps_between(from, to, 90.0, observer)?;
        Ok(steps.into_iter().map(|(q, jd)| MoonPhaseEvent { phase: MoonPhase::from_quadrant(q as i64), jd: JdUt(jd), engine: self.engine_used(jd) }).collect())
    }

    fn tithi_event(&self, (step, jd): (u32, f64)) -> TithiEvent {
        TithiEvent { tithi: Tithi(step as u8 + 1), jd: JdUt(jd), engine: self.engine_used(jd) }
    }

    fn karana_event(&self, (step, jd): (u32, f64)) -> KaranaEvent {
        KaranaEvent { karana: Karana(step as u8 + 1), jd: JdUt(jd), engine: self.engine_used(jd) }
    }

    /// The tithi in force at a moment, from the geocentric elongation.
    pub fn tithi_at(&self, jd: impl Into<JdUt>) -> Result<Calc<Tithi>, CalcError> {
        let jd = jd.into().value();
        Ok(self.calc(jd, Tithi::from_elongation(self.elongation_at(jd, None)?.0)))
    }

    /// The first beginning of a tithi at or after `jd` (geocentric, as tithis always are).
    pub fn next_tithi(&self, jd: impl Into<JdUt>) -> Result<TithiEvent, CalcError> {
        Ok(self.tithi_event(self.next_elongation_step(jd.into().value(), Tithi::SPAN, None)?))
    }

    /// The beginnings of all tithis between two times. The 1st begins at the new moon and the
    /// 16th at the full moon, at the times [`moon_phases`](Self::moon_phases) gives.
    pub fn tithis(&self, from: impl Into<JdUt>, to: impl Into<JdUt>) -> Result<Vec<TithiEvent>, CalcError> {
        let steps = self.elongation_steps_between(from.into().value(), to.into().value(), Tithi::SPAN, None)?;
        Ok(steps.into_iter().map(|s| self.tithi_event(s)).collect())
    }

    /// The karana in force at a moment, from the geocentric elongation.
    pub fn karana_at(&self, jd: impl Into<JdUt>) -> Result<Calc<Karana>, CalcError> {
        let jd = jd.into().value();
        Ok(self.calc(jd, Karana::from_elongation(self.elongation_at(jd, None)?.0)))
    }

    /// The first beginning of a karana at or after `jd` (geocentric).
    pub fn next_karana(&self, jd: impl Into<JdUt>) -> Result<KaranaEvent, CalcError> {
        Ok(self.karana_event(self.next_elongation_step(jd.into().value(), Karana::SPAN, None)?))
    }

    /// The beginnings of all karanas between two times: every tithi boundary and its midpoint,
    /// the four principal phases among them.
    pub fn karanas(&self, from: impl Into<JdUt>, to: impl Into<JdUt>) -> Result<Vec<KaranaEvent>, CalcError> {
        let steps = self.elongation_steps_between(from.into().value(), to.into().value(), Karana::SPAN, None)?;
        Ok(steps.into_iter().map(|s| self.karana_event(s)).collect())
    }

    /// The first principal lunar phase at or after `jd` (geocentric, the astronomical convention).
    pub fn next_moon_phase(&self, jd: impl Into<JdUt>) -> Result<MoonPhaseEvent, CalcError> {
        self.next_phase_at(jd.into().value(), None)
    }

    /// All principal lunar phases between two times (geocentric).
    pub fn moon_phases(&self, from: impl Into<JdUt>, to: impl Into<JdUt>) -> Result<Vec<MoonPhaseEvent>, CalcError> {
        self.phases_between(from.into().value(), to.into().value(), None)
    }

    /// Lunar phases as seen from a place on Earth (topocentric elongation), which can differ
    /// from the geocentric times by up to about two hours. This is what the AstroCalc API returns.
    pub fn moon_phases_topocentric(&self, from: impl Into<JdUt>, to: impl Into<JdUt>, geo: GeoPos) -> Result<Vec<MoonPhaseEvent>, CalcError> {
        self.phases_between(from.into().value(), to.into().value(), Some(geo))
    }
}
