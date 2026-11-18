//! Transposed transits: the rising, setting and meridian passages of a fixed point of the
//! zodiac, such as the place a body held at a birth, carried onto the sky of a later day.
//!
//! The point is given by its ecliptic longitude and latitude and is held at those
//! coordinates in the ecliptic of date. A fixed point turns with the sky at the sidereal
//! rate, so nothing has to be sampled: each event is the moment the point's hour angle
//! reaches a known value,
//!
//! * upper transit (`mc`): `H = 0`,
//! * lower transit (`ic`): `H = 180°`,
//! * rising and setting: `H = ∓H₀`, where `cos H₀ = (sin h₀ − sin φ sin δ) / (cos φ cos δ)`
//!
//! for latitude `φ`, declination `δ` and horizon altitude `h₀` (Meeus, *Astronomical
//! Algorithms*, ch. 15). The time follows from the sidereal time at the start, and is then
//! corrected with the sidereal time and obliquity at the event itself.

use std::f64::consts::{PI, TAU};

use crate::body::Body;
use crate::error::CalcError;
use crate::geo::GeoPos;
use crate::models::{DeltaTModel, EarthOrientation, ut_to_tt};
use crate::position::{BodyPos, BodyPosExt, CoordSystem, Perspective};
use crate::rise_set::{RiseSetOptions, local_day_start};
use crate::spk::Ephemeris;
use crate::time::JdUt;

/// Rotation of the Earth against the equinox, in radians per day of Universal Time.
const SIDEREAL_RATE: f64 = TAU * 1.002_737_909_35;
/// A correction smaller than this ends the refinement (days; about 0.1 ms).
const TOLERANCE: f64 = 1e-9;

/// The first rising, setting and meridian passages of a fixed point after a starting time.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct TransposedTransits {
    /// The body whose place was transposed.
    pub body: Body,
    /// The ecliptic longitude and latitude (degrees) at which the point is held.
    pub lng: f64,
    pub lat: f64,
    /// The time the search began. Every event falls within a sidereal day after it.
    pub start: JdUt,
    /// `None` when the point does not cross the horizon: it is always up or always down.
    pub rise: Option<JdUt>,
    /// Upper meridian transit (culmination).
    pub mc: JdUt,
    pub set: Option<JdUt>,
    /// Lower meridian transit.
    pub ic: JdUt,
    /// Geometric altitude at `ic` and `mc` (degrees).
    pub min_alt: f64,
    pub max_alt: f64,
}

#[derive(Clone, Copy)]
enum Event {
    Rise,
    Mc,
    Set,
    Ic,
}

/// The point and the observer, with the angles that do not change from moment to moment.
struct Transposed {
    /// Ecliptic longitude and latitude of the point, radians.
    lng: f64,
    lat: f64,
    /// Geographic longitude (east positive) and latitude, radians.
    geo_lng: f64,
    phi: f64,
    /// Sine of the altitude at which the point counts as rising or setting.
    sin_horizon: f64,
}

/// The point in the sky of one moment.
struct Moment {
    /// Right ascension, radians.
    ra: f64,
    /// Hour angle, radians.
    hour_angle: f64,
    sin_dec: f64,
    cos_dec: f64,
}

impl Transposed {
    fn at(&self, jd: f64) -> Moment {
        let eo = EarthOrientation::new(jd, ut_to_tt(DeltaTModel::Observed, jd));
        let (sin_eps, cos_eps) = eo.eps_true.sin_cos();
        let (sin_lng, cos_lng) = self.lng.sin_cos();
        let (sin_lat, cos_lat) = self.lat.sin_cos();
        let sin_dec = sin_lat * cos_eps + cos_lat * sin_eps * sin_lng;
        let ra = (sin_lng * cos_eps * cos_lat - sin_lat * sin_eps).atan2(cos_lng * cos_lat);
        Moment { ra, hour_angle: eo.gast + self.geo_lng - ra, sin_dec, cos_dec: (1.0 - sin_dec * sin_dec).max(0.0).sqrt() }
    }

    /// Altitude (degrees) on the meridian: the upper transit, or the lower one.
    fn meridian_altitude(&self, m: &Moment, upper: bool) -> f64 {
        let (sin_phi, cos_phi) = self.phi.sin_cos();
        let cos_part = cos_phi * m.cos_dec;
        (sin_phi * m.sin_dec + if upper { cos_part } else { -cos_part }).clamp(-1.0, 1.0).asin().to_degrees()
    }

    /// The hour angle of an event, or `None` if the point never reaches the horizon.
    fn target(&self, m: &Moment, event: Event) -> Option<f64> {
        let semi_arc = || {
            let (sin_phi, cos_phi) = self.phi.sin_cos();
            let cos_h0 = (self.sin_horizon - sin_phi * m.sin_dec) / (cos_phi * m.cos_dec);
            // Also false for the NaN of a point at the zenith of a pole.
            (cos_h0.abs() <= 1.0).then(|| cos_h0.acos())
        };
        match event {
            Event::Mc => Some(0.0),
            Event::Ic => Some(PI),
            Event::Rise => semi_arc().map(|h0| -h0),
            Event::Set => semi_arc(),
        }
    }

    /// The first time after `start` at which the event happens.
    fn first(&self, start: f64, begin: &Moment, event: Event) -> Option<f64> {
        let mut jd = start + (self.target(begin, event)? - begin.hour_angle).rem_euclid(TAU) / SIDEREAL_RATE;
        // The sidereal time is very nearly linear and the obliquity very nearly constant
        // over a day, so the first correction is already under a millisecond.
        for _ in 0..4 {
            let m = self.at(jd);
            let error = (self.target(&m, event)? - m.hour_angle + PI).rem_euclid(TAU) - PI;
            jd += error / SIDEREAL_RATE;
            if (error / SIDEREAL_RATE).abs() < TOLERANCE {
                break;
            }
        }
        Some(jd)
    }
}

impl Transposed {
    fn new(lng: f64, lat: f64, geo: &GeoPos, horizon_alt: f64) -> Self {
        Self {
            lng: lng.to_radians(),
            lat: lat.to_radians(),
            geo_lng: geo.lng.to_radians(),
            phi: geo.lat.to_radians(),
            sin_horizon: horizon_alt.to_radians().sin(),
        }
    }

    /// Azimuth (from south, westwards) and altitude, degrees.
    fn horizontal(&self, jd: f64) -> (f64, f64) {
        let m = self.at(jd);
        let (sin_phi, cos_phi) = self.phi.sin_cos();
        let (sin_h, cos_h) = m.hour_angle.sin_cos();
        let alt = (sin_phi * m.sin_dec + cos_phi * m.cos_dec * cos_h).clamp(-1.0, 1.0).asin();
        let azi = (m.cos_dec * sin_h).atan2(m.cos_dec * cos_h * sin_phi - m.sin_dec * cos_phi);
        (azi.to_degrees().rem_euclid(360.0), alt.to_degrees())
    }
}

/// The first rising, upper transit, setting and lower transit, after `start`, of the place
/// a body held in a source chart, seen from `geo`.
///
/// `pos` is the body's ecliptic position in the source chart: its longitude and latitude
/// are kept and its speeds are not used. `horizon_alt` is the geometric altitude, in
/// degrees, at which the point rises and sets: `0.0` for the level horizon without an
/// atmosphere, or the value of [`RiseSetOptions::point_altitude`] to allow for refraction,
/// a raised horizon or dip. Start from [`local_day_start`] for the events of a local day.
pub fn transposed_point_transits(pos: &BodyPos, start: impl Into<JdUt>, geo: GeoPos, horizon_alt: f64) -> TransposedTransits {
    let start = start.into();
    let point = Transposed::new(pos.lng, pos.lat, &geo, horizon_alt);
    let begin = point.at(start.value());
    let find = |event| point.first(start.value(), &begin, event);
    // The meridian is always reached: only the horizon can be missed.
    let (mc, ic) = (find(Event::Mc).unwrap_or(f64::NAN), find(Event::Ic).unwrap_or(f64::NAN));
    TransposedTransits {
        body: pos.body,
        lng: pos.lng,
        lat: pos.lat,
        start,
        rise: find(Event::Rise).map(JdUt),
        mc: JdUt(mc),
        set: find(Event::Set).map(JdUt),
        ic: JdUt(ic),
        // Taken at the transits themselves: nutation shifts the declination a little in a day.
        min_alt: point.meridian_altitude(&point.at(ic), false),
        max_alt: point.meridian_altitude(&point.at(mc), true),
    }
}

/// [`transposed_point_transits`] for each of `days` successive days, the first beginning at `start`.
pub fn transposed_point_transit_days(pos: &BodyPos, start: impl Into<JdUt>, geo: GeoPos, horizon_alt: f64, days: u32) -> Vec<TransposedTransits> {
    let start = start.into();
    (0..days).map(|day| transposed_point_transits(pos, start.add_days(f64::from(day)), geo, horizon_alt)).collect()
}

/// A body's place in a source chart as a position in the sky of another moment, in all
/// three coordinate systems at once.
///
/// `pos` is the body's ecliptic position in the source chart. Its longitude and latitude
/// are held in the ecliptic of date at `jd`, so `lng` and `lat` are unchanged; `ras` and
/// `dcl` follow from the obliquity at `jd`, and `azi` and `alt` from the sidereal time at
/// `jd` and the place `geo`.
///
/// The result is not the body itself: its `body` is [`Body::transposed`] of the source
/// body, written `tp:su`, and [`Body::source`] gives the source body back.
///
/// The point is fixed among the stars, so its ecliptic and equatorial speeds are zero. Its
/// horizontal speeds are those of the turning sky, in degrees per day. The distance is that
/// of the source position, with no speed.
pub fn transposed_position(pos: &BodyPos, jd: impl Into<JdUt>, geo: GeoPos) -> BodyPosExt {
    // Rates over one second either side.
    const HALF_STEP: f64 = 1.0 / 86_400.0;
    let jd = jd.into().value();
    let point = Transposed::new(pos.lng, pos.lat, &geo, 0.0);
    let m = point.at(jd);
    let (azi, alt) = point.horizontal(jd);
    let (before, after) = (point.horizontal(jd - HALF_STEP), point.horizontal(jd + HALF_STEP));
    let turn = (after.0 - before.0 + 180.0).rem_euclid(360.0) - 180.0;
    BodyPosExt {
        body: Body::transposed(pos.body),
        lng: pos.lng.rem_euclid(360.0),
        lat: pos.lat,
        ras: m.ra.to_degrees().rem_euclid(360.0),
        dcl: m.sin_dec.atan2(m.cos_dec).to_degrees(),
        azi,
        alt,
        lng_speed: 0.0,
        lat_speed: 0.0,
        ras_speed: 0.0,
        dcl_speed: 0.0,
        azi_speed: turn / (2.0 * HALF_STEP),
        alt_speed: (after.1 - before.1) / (2.0 * HALF_STEP),
        dist: pos.dist,
        dist_speed: 0.0,
        engine: pos.engine,
    }
}

impl Ephemeris {
    /// The ecliptic positions of `bodies` in a source chart, ready to be transposed: seen
    /// from `source_geo` if a place is given, otherwise from the centre of the Earth.
    pub fn transposed_sources(&self, bodies: &[Body], source_jd: impl Into<JdUt>, source_geo: Option<GeoPos>) -> Result<Vec<BodyPos>, CalcError> {
        let perspective = if source_geo.is_some() { Perspective::Topocentric } else { Perspective::Geocentric };
        self.calc_positions(bodies, source_jd, source_geo.unwrap_or_default(), CoordSystem::Ecliptic, perspective)
    }

    /// Transposed positions: where each of `bodies` stood at `h_jd`, seen from `h_geo` in
    /// the given perspective, held at that ecliptic longitude and latitude and placed in the
    /// sky of `jd` at `geo`, in all three coordinate systems. See [`transposed_position`].
    ///
    /// Each result names its source: see [`Body::source`]. `h_geo` matters only for a
    /// perspective on the Earth's surface; the Local Space perspective counts as
    /// topocentric. The result includes the horizon, so the perspective must be one on the
    /// Earth, as for [`calc_positions_ext`](Self::calc_positions_ext).
    pub fn calc_transposed_positions(
        &self,
        bodies: &[Body],
        jd: impl Into<JdUt>,
        h_jd: impl Into<JdUt>,
        geo: GeoPos,
        h_geo: GeoPos,
        perspective: Perspective,
    ) -> Result<Vec<BodyPosExt>, CalcError> {
        if !perspective.is_earth_based() {
            return Err(CalcError::NoHorizon(perspective));
        }
        // Local Space would give the source on its own horizon: take it as ecliptic.
        let perspective = if perspective == Perspective::LocalSpaceHorizon { Perspective::Topocentric } else { perspective };
        let jd = jd.into();
        let sources = self.calc_positions(bodies, h_jd, h_geo, CoordSystem::Ecliptic, perspective)?;
        Ok(sources.iter().map(|pos| transposed_position(pos, jd, geo)).collect())
    }

    /// Transposed transits: where each of `bodies` stood at `source_jd` (seen from
    /// `source_geo`, or from the centre of the Earth without one), held at that ecliptic
    /// longitude and latitude and carried onto the sky of the local day containing `jd` at
    /// `geo`. For each body it gives the first rise, upper transit, set and lower transit
    /// after local mean midnight.
    ///
    /// A fixed point has no disc, so of the options only the horizon, refraction and dip
    /// apply, with [`ignore_latitude`](RiseSetOptions::ignore_latitude) to place the point
    /// on the ecliptic and [`geocentric`](RiseSetOptions::geocentric) to take the source
    /// positions from the centre of the Earth even when a source place is given.
    pub fn transposed_transits(
        &self,
        bodies: &[Body],
        source_jd: impl Into<JdUt>,
        source_geo: Option<GeoPos>,
        jd: impl Into<JdUt>,
        geo: GeoPos,
        opts: RiseSetOptions,
    ) -> Result<Vec<TransposedTransits>, CalcError> {
        let start = local_day_start(jd.into().value(), geo.lng);
        let horizon = opts.point_altitude(geo.alt);
        let sources = self.transposed_sources(bodies, source_jd, source_geo.filter(|_| !opts.geocentric))?;
        Ok(sources
            .iter()
            .map(|pos| {
                let pos = BodyPos { lat: if opts.ignore_latitude { 0.0 } else { pos.lat }, ..*pos };
                transposed_point_transits(&pos, start, geo, horizon)
            })
            .collect())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::houses::angles;
    use crate::math::diff_deg;

    const JD: f64 = 2_460_270.5;

    fn point(lng: f64, lat: f64) -> BodyPos {
        BodyPos { body: Body::Sun, lng, lat, lng_speed: 0.0, lat_speed: 0.0, dist: 1.0, dist_speed: 0.0, engine: Default::default() }
    }

    fn transposed_transits(lng: f64, lat: f64, start: f64, geo: GeoPos, horizon_alt: f64) -> TransposedTransits {
        transposed_point_transits(&point(lng, lat), start, geo, horizon_alt)
    }

    /// A point on the ecliptic rises as the ascendant, sets as the descendant and
    /// culminates as the MC: an independent check through the house code.
    #[test]
    fn ecliptic_points_meet_the_chart_angles() {
        for (lat, lng) in [(51.5, -0.1), (-33.9, 151.2), (13.1, 80.3), (64.1, -21.9), (0.0, 30.0)] {
            let geo = GeoPos::new(lat, lng, 0.0).unwrap();
            for point in (0..360).step_by(17).map(f64::from) {
                let t = transposed_transits(point, 0.0, JD, geo, 0.0);
                for jd in [t.mc.value(), t.ic.value(), t.rise.unwrap().value(), t.set.unwrap().value()] {
                    assert!((JD - 1e-9..JD + 1.0).contains(&jd), "{jd} is outside the day");
                }
                let rise = angles(t.rise.unwrap(), geo);
                let set = angles(t.set.unwrap(), geo);
                assert!(diff_deg(rise.ascendant, point).abs() < 1e-6, "rise at {lat},{lng}: {} vs {point}", rise.ascendant);
                assert!(diff_deg(set.descendant(), point).abs() < 1e-6, "set at {lat},{lng}: {} vs {point}", set.descendant());
                assert!(diff_deg(angles(t.mc, geo).mc, point).abs() < 1e-6);
                assert!(diff_deg(angles(t.ic, geo).ic(), point).abs() < 1e-6);
            }
        }
    }

    /// Altitude of the point, by the direct formula, for checks away from the ecliptic.
    fn altitude(lng: f64, lat: f64, jd: f64, geo: &GeoPos) -> f64 {
        transposed_position(&point(lng, lat), jd, *geo).alt
    }

    #[test]
    fn events_agree_with_the_altitude() {
        let geo = GeoPos::new(56.1, -3.6, 0.0).unwrap();
        for (lng, lat, horizon) in [(123.4, 5.2, 0.0), (301.0, -4.9, -0.6123), (45.0, 1.0, 2.5)] {
            let t = transposed_transits(lng, lat, JD, geo, horizon);
            let (rise, set) = (t.rise.unwrap().value(), t.set.unwrap().value());
            assert!((altitude(lng, lat, rise, &geo) - horizon).abs() < 1e-7);
            assert!((altitude(lng, lat, set, &geo) - horizon).abs() < 1e-7);
            // Going up at the rise and down at the set.
            assert!(altitude(lng, lat, rise + 1e-4, &geo) > horizon && altitude(lng, lat, set + 1e-4, &geo) < horizon);
            // The meridian altitudes are the extremes.
            assert!((altitude(lng, lat, t.mc.value(), &geo) - t.max_alt).abs() < 1e-6);
            assert!((altitude(lng, lat, t.ic.value(), &geo) - t.min_alt).abs() < 1e-6);
            for dt in [-0.01, 0.01] {
                assert!(altitude(lng, lat, t.mc.value() + dt, &geo) < t.max_alt);
                assert!(altitude(lng, lat, t.ic.value() + dt, &geo) > t.min_alt);
            }
        }
    }

    #[test]
    fn points_that_never_cross_the_horizon() {
        let arctic = GeoPos::new(78.2, 15.6, 0.0).unwrap();
        // Near the summer solstice point: always up. Near the winter one: always down.
        let up = transposed_transits(90.0, 0.0, JD, arctic, 0.0);
        let down = transposed_transits(270.0, 0.0, JD, arctic, 0.0);
        assert!(up.rise.is_none() && up.set.is_none() && up.min_alt > 0.0);
        assert!(down.rise.is_none() && down.set.is_none() && down.max_alt < 0.0);
        assert!(up.mc.value().is_finite() && down.ic.value().is_finite());
    }

    #[test]
    fn successive_days_step_by_a_sidereal_day() {
        let geo = GeoPos::new(28.6, 77.2, 0.0).unwrap();
        let days = transposed_point_transit_days(&point(200.0, 2.0), JD, geo, 0.0, 5);
        assert_eq!(days.len(), 5);
        for pair in days.windows(2) {
            let gap = pair[1].mc.value() - pair[0].mc.value();
            // A sidereal day is 3 min 56 s short of a day; sometimes the next window catches a later turn.
            assert!((gap - 0.99727).abs() < 1e-4 || (gap - 2.0 * 0.99727).abs() < 1e-4, "{gap}");
        }
    }
}
