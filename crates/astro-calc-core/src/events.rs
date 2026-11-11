//! Planetary stations and lunar phases, found by root- and extremum-finding on apparent
//! geocentric ecliptic positions.

use crate::body::Body;
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
    pub fn stations(&self, body: Body, from: impl Into<JdUt>, to: impl Into<JdUt>) -> Result<Vec<Station>, CalcError> {
        self.stations_in(body, from.into().value(), to.into().value(), true)
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
        Ok(StationSet { body, jd, speed, stations })
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
        Ok(StationSet { body, jd, speed, stations })
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
    pub fn moon_elongation(&self, jd: impl Into<JdUt>) -> Result<(f64, f64), CalcError> {
        self.elongation_at(jd.into().value(), None)
    }

    fn next_phase_at(&self, start: f64, observer: Option<GeoPos>) -> Result<MoonPhaseEvent, CalcError> {
        let (e, _) = self.elongation_at(start, observer)?;
        let quadrant = (e / 90.0).ceil() as i64;
        let target = (quadrant as f64 * 90.0).rem_euclid(360.0);
        // first estimate from the mean synodic rate, then Newton steps with the actual rate
        let mut t = start + (quadrant as f64 * 90.0 - e) / 12.190_749;
        for _ in 0..20 {
            let (e, rate) = self.elongation_at(t, observer)?;
            let dt = diff_deg(e, target) / rate;
            t += dt;
            if dt.abs() < 1e-9 {
                break;
            }
        }
        Ok(MoonPhaseEvent { phase: MoonPhase::from_quadrant(quadrant), jd: JdUt(t.max(start)) })
    }

    fn phases_between(&self, from: f64, to: f64, observer: Option<GeoPos>) -> Result<Vec<MoonPhaseEvent>, CalcError> {
        let mut t = from;
        let mut out = Vec::new();
        loop {
            let ev = self.next_phase_at(t, observer)?;
            if ev.jd.value() > to {
                return Ok(out);
            }
            out.push(ev);
            t = ev.jd.value() + 1.0; // phases are at least ~6.5 days apart
        }
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
