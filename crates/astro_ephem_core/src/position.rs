//! Apparent positions of solar-system bodies.
//!
//! Pipeline for each body (Explanatory Supplement 3rd ed. ch. 7; IERS Conventions 2010):
//! 1. UT → TT (ΔT) → TDB.
//! 2. Observer's barycentric position and velocity: Earth from the ephemeris, plus for
//!    topocentric mode the WGS84 site vector rotated by sidereal time and precession-nutation.
//! 3. Light-time iteration for the body's barycentric position.
//! 4. Gravitational light deflection by the Sun.
//! 5. Relativistic aberration from the observer's barycentric velocity (annual + diurnal).
//! 6. Rotation to the true equator and equinox of date (bias, precession, nutation), then to the
//!    true ecliptic of date if requested.
//!
//! Speeds are central differences of the full pipeline over ±[`SPEED_STEP`] days.

use crate::body::Body;
use crate::error::CalcError;
use crate::geo::GeoPos;
use crate::math::{Vec3, add, cross, dot, mat_t_vec, mat_vec, norm, scale, sub, to_spherical, unit};
use crate::models::{EARTH_ROT_RATE, EarthOrientation, geodetic_to_itrs, tdb_minus_tt, ut_to_tt};
use crate::spk::{AU_KM, Ephemeris};
use crate::time::JdUt;

/// Speed of light in AU per day.
const C_AU_DAY: f64 = 299_792.458 * 86_400.0 / AU_KM;
/// Schwarzschild radius of the Sun in AU (2GM☉/c²), as used in the deflection formula.
const SRS: f64 = 1.974_125_743_36e-8;
/// Half-width of the central difference used for speeds, in days.
pub const SPEED_STEP: f64 = 0.001;

const NAIF_SUN: i32 = 10;
const NAIF_EARTH: i32 = 399;

/// Output coordinate system.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum CoordSystem {
    /// Ecliptic longitude/latitude of date.
    #[default]
    Ecliptic,
    /// Right ascension/declination of date.
    Equatorial,
    /// Azimuth/altitude (not yet implemented).
    Horizontal,
}

/// Point of observation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum Perspective {
    /// Centre of the Earth.
    #[default]
    Geocentric,
    /// Observer on the WGS84 ellipsoid at the given location.
    Topocentric,
}

/// Position of one body. In ecliptic mode `lng`/`lat` are longitude and latitude; in
/// equatorial mode right ascension and declination. All angles in degrees, speeds in
/// degrees per day, distance in AU and AU per day.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct BodyPos {
    pub body: Body,
    pub lng: f64,
    pub lat: f64,
    pub lng_speed: f64,
    pub lat_speed: f64,
    pub dist: f64,
    pub dist_speed: f64,
}

/// A position at a given time.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct BodyJdPos {
    pub jd: JdUt,
    pub pos: BodyPos,
}

/// The positions of one body over a series of times.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct BodyPosSet {
    pub body: Body,
    pub items: Vec<BodyJdPos>,
}

/// NAIF id of the ephemeris target for a body, if it comes straight from the DE files.
fn naif_id(body: Body) -> Option<i32> {
    Some(match body {
        Body::Sun => NAIF_SUN,
        Body::Moon => 301,
        Body::Mercury => 199,
        Body::Venus => 299,
        Body::Mars => 4,
        Body::Jupiter => 5,
        Body::Saturn => 6,
        Body::Uranus => 7,
        Body::Neptune => 8,
        Body::Pluto => 9,
        _ => return None,
    })
}

/// Everything about one instant and observer that every body shares.
pub(crate) struct Frame {
    pub(crate) tdb: f64,
    pub(crate) eo: EarthOrientation,
    /// Observer barycentric position (AU) and velocity (AU/day).
    obs_pos: Vec3,
    obs_vel: Vec3,
    /// Sun barycentric position at `tdb` (AU).
    sun_pos: Vec3,
}

impl Frame {
    pub(crate) fn new(eph: &Ephemeris, jd_ut: f64, geo: &GeoPos, persp: Perspective) -> Result<Self, CalcError> {
        let tt = ut_to_tt(eph.delta_t_model(), jd_ut);
        let tdb = tt + tdb_minus_tt(tt) / 86_400.0;
        let eo = EarthOrientation::new(jd_ut, tt);
        let (mut obs_pos, mut obs_vel) = eph.barycentric(NAIF_EARTH, tdb)?;
        if persp == Perspective::Topocentric {
            let site = geodetic_to_itrs(geo.lat.to_radians(), geo.lng.to_radians(), geo.alt);
            let site_au = scale(site, 1e-3 / AU_KM);
            // Earth-fixed → true equator of date: rotate by −GAST about z.
            let (s, c) = eo.gast.sin_cos();
            let r_tod = [c * site_au[0] - s * site_au[1], s * site_au[0] + c * site_au[1], site_au[2]];
            let v_tod = cross([0.0, 0.0, EARTH_ROT_RATE], r_tod);
            obs_pos = add(obs_pos, mat_t_vec(&eo.npb, r_tod));
            obs_vel = add(obs_vel, mat_t_vec(&eo.npb, v_tod));
        }
        let (sun_pos, _) = eph.barycentric(NAIF_SUN, tdb)?;
        Ok(Self { tdb, eo, obs_pos, obs_vel, sun_pos })
    }

    /// Apparent direction (unit vector, true equator of date) and light-time-corrected distance.
    fn apparent(&self, eph: &Ephemeris, target: i32) -> Result<(Vec3, f64), CalcError> {
        // Light-time iteration (converges to < 1e-12 day within a few passes).
        let mut tau = 0.0;
        let mut p = [0.0; 3];
        let mut body_bary = [0.0; 3];
        for _ in 0..4 {
            let (b, _) = eph.barycentric(target, self.tdb - tau)?;
            body_bary = b;
            p = sub(b, self.obs_pos);
            let new_tau = norm(p) / C_AU_DAY;
            let done = (new_tau - tau).abs() < 1e-12;
            tau = new_tau;
            if done {
                break;
            }
        }
        let dist = norm(p);
        let dir = unit(p);

        let dir = if target != NAIF_SUN {
            self.deflect(dir, unit(sub(body_bary, self.sun_pos)))
        } else {
            dir
        };
        let dir = self.aberrate(dir);
        Ok((mat_vec(&self.eo.npb, dir), dist))
    }

    /// Gravitational light deflection by the Sun. `dir`: observer → source (unit);
    /// `q`: Sun → source (unit; equal to `dir` for a source at infinity).
    fn deflect(&self, dir: Vec3, q: Vec3) -> Vec3 {
        let sun_to_obs = sub(self.obs_pos, self.sun_pos);
        let em = norm(sun_to_obs);
        let e = scale(sun_to_obs, 1.0 / em);
        let qpe = add(q, e);
        let dlim = 1e-6 / (em * em).max(1.0);
        let w = SRS / em / dot(q, qpe).max(dlim);
        let peq = cross(dir, cross(e, q));
        unit(add(dir, scale(peq, w)))
    }

    /// Relativistic aberration from the observer's barycentric velocity.
    fn aberrate(&self, dir: Vec3) -> Vec3 {
        let v = scale(self.obs_vel, 1.0 / C_AU_DAY);
        let bm1 = (1.0 - dot(v, v)).sqrt();
        let s = norm(sub(self.obs_pos, self.sun_pos));
        let pdv = dot(dir, v);
        let w1 = 1.0 + pdv / (1.0 + bm1);
        let w2 = SRS / s;
        unit([0, 1, 2].map(|i| dir[i] * bm1 + w1 * v[i] + w2 * (v[i] - pdv * dir[i])))
    }

    /// Apparent true-ecliptic-of-date longitude (degrees) of a fixed ICRS direction.
    pub(crate) fn apparent_ecliptic_lng_of_star(&self, icrs_dir: Vec3) -> f64 {
        let dir = self.aberrate(self.deflect(icrs_dir, icrs_dir));
        let u = mat_vec(&self.eo.npb, dir);
        let (s, c) = self.eo.eps_true.sin_cos();
        let (lng, _) = to_spherical([u[0], c * u[1] + s * u[2], -s * u[1] + c * u[2]]);
        lng.to_degrees()
    }

    /// Spherical coordinates in degrees plus distance.
    fn coords(&self, eph: &Ephemeris, target: i32, cs: CoordSystem) -> Result<(f64, f64, f64), CalcError> {
        let (u, dist) = self.apparent(eph, target)?;
        let v = match cs {
            CoordSystem::Equatorial => u,
            CoordSystem::Ecliptic => {
                let (s, c) = self.eo.eps_true.sin_cos();
                [u[0], c * u[1] + s * u[2], -s * u[1] + c * u[2]]
            }
            CoordSystem::Horizontal => return Err(CalcError::UnsupportedCoordSystem(cs)),
        };
        let (lng, lat) = to_spherical(v);
        Ok((lng.to_degrees(), lat.to_degrees(), dist))
    }
}

/// Frames at t − h, t and t + h, shared by every body at one time.
struct FrameTriple([Frame; 3]);

impl FrameTriple {
    fn new(eph: &Ephemeris, jd: f64, geo: &GeoPos, persp: Perspective) -> Result<Self, CalcError> {
        Ok(Self([
            Frame::new(eph, jd - SPEED_STEP, geo, persp)?,
            Frame::new(eph, jd, geo, persp)?,
            Frame::new(eph, jd + SPEED_STEP, geo, persp)?,
        ]))
    }

    fn body_pos(&self, eph: &Ephemeris, body: Body, cs: CoordSystem) -> Result<BodyPos, CalcError> {
        let target = naif_id(body).ok_or(CalcError::UnsupportedBody(body))?;
        let [a, m, b] = &self.0;
        let (l0, b0, d0) = a.coords(eph, target, cs)?;
        let (l1, b1, d1) = m.coords(eph, target, cs)?;
        let (l2, b2, d2) = b.coords(eph, target, cs)?;
        let h2 = 2.0 * SPEED_STEP;
        Ok(BodyPos {
            body,
            lng: l1,
            lat: b1,
            lng_speed: crate::math::diff_deg(l0, l2) / h2,
            lat_speed: (b2 - b0) / h2,
            dist: d1,
            dist_speed: (d2 - d0) / h2,
        })
    }
}

fn check_cs(cs: CoordSystem) -> Result<(), CalcError> {
    if cs == CoordSystem::Horizontal {
        Err(CalcError::UnsupportedCoordSystem(cs))
    } else {
        Ok(())
    }
}

/// Times `start + i·interval` for `i` in `0..ceil(days / interval)` (at least one sample).
fn sample_times(start: f64, days: f64, interval: f64) -> Result<Vec<f64>, CalcError> {
    if !(interval.is_finite() && interval > 0.0 && days.is_finite() && days >= 0.0) {
        return Err(CalcError::InvalidInterval { days, interval });
    }
    let n = ((days / interval).ceil() as usize).max(1);
    Ok((0..n).map(|i| start + i as f64 * interval).collect())
}

impl Ephemeris {
    /// Apparent position of one body.
    pub fn calc_position(
        &self,
        body: Body,
        jd: impl Into<JdUt>,
        geo: GeoPos,
        cs: CoordSystem,
        perspective: Perspective,
    ) -> Result<BodyPos, CalcError> {
        check_cs(cs)?;
        let frames = FrameTriple::new(self, jd.into().value(), &geo, perspective)?;
        frames.body_pos(self, body, cs)
    }

    /// Apparent positions of several bodies at one time. Time-dependent quantities are computed
    /// once and shared by all bodies.
    pub fn calc_positions(
        &self,
        bodies: &[Body],
        jd: impl Into<JdUt>,
        geo: GeoPos,
        cs: CoordSystem,
        perspective: Perspective,
    ) -> Result<Vec<BodyPos>, CalcError> {
        check_cs(cs)?;
        let frames = FrameTriple::new(self, jd.into().value(), &geo, perspective)?;
        bodies.iter().map(|&b| frames.body_pos(self, b, cs)).collect()
    }

    /// Positions of one body at `start + i·interval` over `days` (end exclusive).
    #[allow(clippy::too_many_arguments)]
    pub fn calc_progress(
        &self,
        body: Body,
        start: impl Into<JdUt>,
        days: f64,
        interval: f64,
        geo: GeoPos,
        cs: CoordSystem,
        perspective: Perspective,
    ) -> Result<Vec<BodyJdPos>, CalcError> {
        let mut sets = self.calc_progress_positions(&[body], start, days, interval, geo, cs, perspective)?;
        Ok(sets.pop().map(|s| s.items).unwrap_or_default())
    }

    /// Positions of several bodies at `start + i·interval` over `days` (end exclusive),
    /// grouped by body.
    #[allow(clippy::too_many_arguments)]
    pub fn calc_progress_positions(
        &self,
        bodies: &[Body],
        start: impl Into<JdUt>,
        days: f64,
        interval: f64,
        geo: GeoPos,
        cs: CoordSystem,
        perspective: Perspective,
    ) -> Result<Vec<BodyPosSet>, CalcError> {
        check_cs(cs)?;
        let times = sample_times(start.into().value(), days, interval)?;
        let mut sets: Vec<BodyPosSet> = bodies
            .iter()
            .map(|&body| BodyPosSet { body, items: Vec::with_capacity(times.len()) })
            .collect();
        for &t in &times {
            let frames = FrameTriple::new(self, t, &geo, perspective)?;
            for set in sets.iter_mut() {
                let pos = frames.body_pos(self, set.body, cs)?;
                set.items.push(BodyJdPos { jd: JdUt(t), pos });
            }
        }
        Ok(sets)
    }
}
