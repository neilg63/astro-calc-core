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
/// GM of the Earth plus the Moon in AU³/day² (DE440: 398 600.435 507 + 4 902.800 118 km³/s²).
const GM_EARTH_MOON: f64 = 403_503.235_625 * 86_400.0 * 86_400.0 / (AU_KM * AU_KM * AU_KM);
/// Mean inclination of the lunar orbit to the ecliptic, degrees.
const MEAN_LUNAR_INCLINATION: f64 = 5.145_396_4;
/// Mean apogee distance of the Moon, a(1 + e) with a = 384 400 km and e = 0.0549, in AU.
const MEAN_LUNAR_APOGEE_AU: f64 = 405_504.0 / AU_KM;

/// Output coordinate system.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum CoordSystem {
    /// Ecliptic longitude/latitude of date.
    #[default]
    Ecliptic,
    /// Right ascension/declination of date.
    Equatorial,
    /// Azimuth (from south, westwards) and geometric altitude.
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
pub(crate) fn naif_id(body: Body) -> Option<i32> {
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
        // Small bodies use JPL Horizons SPK ids (20 000 000 + minor-planet number).
        Body::Ceres => 20_000_001,
        Body::Pallas => 20_000_002,
        Body::Juno => 20_000_003,
        Body::Vesta => 20_000_004,
        Body::Chiron => 20_002_060,
        Body::Pholus => 20_005_145,
        Body::Catalog { catalog: crate::body::Catalog::JplHorizons, id } => i32::try_from(id).ok()?,
        Body::Catalog { catalog: crate::body::Catalog::MinorPlanetCenter, id } => {
            i32::try_from(id.checked_add(20_000_000)?).ok()?
        }
        _ => return None,
    })
}

/// Everything about one instant and observer that every body shares.
pub(crate) struct Frame {
    pub(crate) tt: f64,
    pub(crate) tdb: f64,
    pub(crate) eo: EarthOrientation,
    /// Observer latitude and longitude (radians), used for horizontal coordinates.
    obs_lat: f64,
    obs_lng: f64,
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
        Ok(Self {
            tt,
            tdb,
            eo,
            obs_lat: geo.lat.to_radians(),
            obs_lng: geo.lng.to_radians(),
            obs_pos,
            obs_vel,
            sun_pos,
        })
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

    /// Apparent geocentric geometry of a body for phenomena.
    pub(crate) fn geometry(&self, eph: &Ephemeris, body: Body) -> Result<crate::pheno::Geometry, CalcError> {
        let target = naif_id(body).ok_or(CalcError::UnsupportedBody(body))?;
        let (dir_b, dist) = self.apparent(eph, target)?;
        let (dir_s, dist_s) = self.apparent(eph, NAIF_SUN)?;
        // apparent geocentric vectors (true equator of date)
        let (b, s) = (scale(dir_b, dist), scale(dir_s, dist_s));
        let body_to_sun = sub(s, b);
        let phase_angle = crate::pheno::angle_between(body_to_sun, scale(b, -1.0));
        let elongation = crate::pheno::angle_between(b, s);
        // Saturn: sine of the sub-Earth latitude from the IAU pole (α₀ = 40.589°, δ₀ = 83.537°, ICRF).
        let ring_tilt_sin = if body == Body::Saturn {
            let t = (self.tt - 2_451_545.0) / 36_525.0;
            let (ra, dec) = ((40.589 - 0.036 * t).to_radians(), (83.537 - 0.004 * t).to_radians());
            let pole = mat_vec(&self.eo.npb, [dec.cos() * ra.cos(), dec.cos() * ra.sin(), dec.sin()]);
            -dot(pole, dir_b)
        } else {
            0.0
        };
        Ok(crate::pheno::Geometry { phase_angle, elongation, dist, helio_dist: norm(body_to_sun), ring_tilt_sin })
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

    /// ICRS vector → true ecliptic and equinox of date.
    fn to_ecl(&self, x: Vec3) -> Vec3 {
        let u = mat_vec(&self.eo.npb, x);
        let (s, c) = self.eo.eps_true.sin_cos();
        [u[0], c * u[1] + s * u[2], -s * u[1] + c * u[2]]
    }

    /// True ecliptic of date → true equator of date.
    fn ecl_to_equ(&self, e: Vec3) -> Vec3 {
        let (s, c) = self.eo.eps_true.sin_cos();
        [e[0], c * e[1] - s * e[2], s * e[1] + c * e[2]]
    }

    /// Ecliptic longitude (degrees, true ecliptic and equinox of date) of a lunar node.
    fn node_lng(&self, eph: &Ephemeris, kind: Source) -> Result<f64, CalcError> {
        Ok(match kind {
            Source::MeanNode => crate::models::mean_lunar_node(self.tt) + self.eo.dpsi.to_degrees(),
            _ => {
                // Osculating node: direction of k × h for the geometric geocentric lunar orbit,
                // h = r × v, in the true ecliptic of date.
                let (r, v) = eph.relative(301, NAIF_EARTH, self.tdb)?;
                let h = cross(self.to_ecl(r), self.to_ecl(v));
                h[0].atan2(-h[1]).to_degrees()
            }
        })
    }

    /// A lunar apsis as a vector (AU) in the true ecliptic of date.
    fn apsis(&self, eph: &Ephemeris, kind: Apsis) -> Result<Vec3, CalcError> {
        match kind {
            Apsis::MeanApogee => {
                // The mean orbit: mean node, mean inclination, mean perigee + 180°.
                let node = crate::models::mean_lunar_node(self.tt).to_radians();
                let arg = (crate::models::mean_lunar_perigee(self.tt) + 180.0).to_radians() - node;
                let (si, ci) = MEAN_LUNAR_INCLINATION.to_radians().sin_cos();
                let lng = node + (ci * arg.sin()).atan2(arg.cos()) + self.eo.dpsi;
                let lat = (si * arg.sin()).asin();
                Ok(scale([lat.cos() * lng.cos(), lat.cos() * lng.sin(), lat.sin()], MEAN_LUNAR_APOGEE_AU))
            }
            Apsis::OscuApogee => {
                // Apogee of the two-body ellipse that fits the Moon's geometric geocentric
                // position and velocity at this instant: opposite the eccentricity vector.
                let (r, v) = eph.relative(301, NAIF_EARTH, self.tdb)?;
                let (r, v) = (self.to_ecl(r), self.to_ecl(v));
                let rn = norm(r);
                let ev = sub(scale(cross(v, cross(r, v)), 1.0 / GM_EARTH_MOON), scale(r, 1.0 / rn));
                let e = norm(ev);
                let a = 1.0 / (2.0 / rn - dot(v, v) / GM_EARTH_MOON);
                Ok(scale(ev, -a * (1.0 + e) / e))
            }
            Apsis::IntpApogee | Apsis::IntpPerigee => self.interpolated_apsis(eph, kind == Apsis::IntpApogee),
        }
    }

    /// The Moon's actual apogees (or perigees) joined into a continuous curve.
    ///
    /// The Moon is at its greatest (least) distance once per anomalistic month. Its position
    /// at each of those moments is a real apsis. Between them the apsis is taken from a cubic
    /// through the four nearest ones, interpolating the offset from the mean apogee
    /// (perigee) in longitude, the latitude and the distance.
    fn interpolated_apsis(&self, eph: &Ephemeris, apogee: bool) -> Result<Vec3, CalcError> {
        const ANOMALISTIC_MONTH: f64 = 27.554_549_9;
        let mean_lng = |tt: f64| crate::models::mean_lunar_perigee(tt) + if apogee { 180.0 } else { 0.0 };
        let radial = |t: f64| eph.relative(301, NAIF_EARTH, t).map(|(r, v)| dot(r, v));
        // Position of the mean cycle: whole numbers are mean apogees (perigees).
        let phase = (crate::models::mean_lunar_anomaly(self.tt) - if apogee { 180.0 } else { 0.0 }) / 360.0;
        let base = phase.floor();
        let mut events: Vec<(f64, f64, f64, f64)> = Vec::with_capacity(6);
        for j in -2..=3 {
            let guess = self.tdb - (phase - base - j as f64) * ANOMALISTIC_MONTH;
            // The radial velocity falls through zero at an apogee and rises through it at a perigee.
            let (mut a, mut b) = (guess - 5.0, guess + 5.0);
            let before = radial(a)? > 0.0;
            if before != apogee || (radial(b)? > 0.0) == before {
                continue;
            }
            while b - a > 1e-7 {
                let m = 0.5 * (a + b);
                if (radial(m)? > 0.0) == before {
                    a = m;
                } else {
                    b = m;
                }
            }
            let t = 0.5 * (a + b);
            let p = self.to_ecl(eph.relative(301, NAIF_EARTH, t)?.0);
            let (lng, lat) = to_spherical(p);
            events.push((t, crate::math::diff_deg(mean_lng(t), lng.to_degrees()), lat, norm(p)));
        }
        if events.len() < 4 {
            return Err(CalcError::Ephem(crate::error::EphemError::NoData { target: 301, jd: self.tdb }));
        }
        // The four events nearest in time, two either side where possible.
        let after = events.partition_point(|e| e.0 <= self.tdb);
        let first = after.saturating_sub(2).min(events.len() - 4);
        let pts = &events[first..first + 4];
        let lagrange = |pick: fn(&(f64, f64, f64, f64)) -> f64| {
            (0..4)
                .map(|i| {
                    let w: f64 = (0..4).filter(|&k| k != i).map(|k| (self.tdb - pts[k].0) / (pts[i].0 - pts[k].0)).product();
                    w * pick(&pts[i])
                })
                .sum::<f64>()
        };
        let lng = (mean_lng(self.tt) + lagrange(|e| e.1)).to_radians();
        let (lat, dist) = (lagrange(|e| e.2), lagrange(|e| e.3));
        Ok(scale([lat.cos() * lng.cos(), lat.cos() * lng.sin(), lat.sin()], dist))
    }

    /// Apparent direction in the true equator of date (unit vector) and distance (AU; 0 for
    /// the lunar nodes).
    fn direction(&self, eph: &Ephemeris, src: Source, opposite: bool) -> Result<(Vec3, f64), CalcError> {
        let (u, dist) = match src {
            Source::Naif(target) => self.apparent(eph, target)?,
            Source::MeanNode | Source::TrueNode => {
                let l = self.node_lng(eph, src)?.to_radians();
                let (s, c) = self.eo.eps_true.sin_cos();
                let (sl, cl) = l.sin_cos();
                ([cl, c * sl, s * sl], 0.0)
            }
            Source::Star(star) => {
                // Catalogue direction carried forward by proper motion, then shifted by
                // annual parallax, bent by the Sun and aberrated like any other source.
                let years = (self.tt - crate::stars::HIPPARCOS_EPOCH_JD) / 365.25;
                let (sa, ca) = star.ra.to_radians().sin_cos();
                let (sd, cd) = star.dec.to_radians().sin_cos();
                let mas = (1.0f64 / 3.6e6).to_radians();
                let (ma, md) = (star.pm_ra * mas * years, star.pm_dec * mas * years);
                let east = [-sa, ca, 0.0];
                let north = [-sd * ca, -sd * sa, cd];
                let mut p = [0, 1, 2].map(|i| [cd * ca, cd * sa, sd][i] + ma * east[i] + md * north[i]);
                if star.parallax > 0.0 {
                    p = sub(unit(p), scale(self.obs_pos, star.parallax * mas));
                }
                let dir = unit(p);
                let dir = self.aberrate(self.deflect(dir, dir));
                (mat_vec(&self.eo.npb, dir), star.distance_au())
            }
            Source::Kepler(el) => {
                // A heliocentric orbit, treated like a real body: light-time, deflection and
                // aberration all apply.
                let mut tau = 0.0;
                let mut p = [0.0; 3];
                let mut bary = [0.0; 3];
                for _ in 0..3 {
                    bary = add(el.heliocentric(self.tt - tau), eph.barycentric(NAIF_SUN, self.tdb - tau)?.0);
                    p = sub(bary, self.obs_pos);
                    tau = norm(p) / C_AU_DAY;
                }
                let dir = self.aberrate(self.deflect(unit(p), unit(sub(bary, self.sun_pos))));
                (mat_vec(&self.eo.npb, dir), norm(p))
            }
            Source::Apsis(kind) => {
                let p = self.ecl_to_equ(self.apsis(eph, kind)?);
                (unit(p), norm(p))
            }
        };
        Ok((if opposite { scale(u, -1.0) } else { u }, dist))
    }

    /// Geometric altitude and hour angle (radians) of a body's centre, plus distance (AU).
    pub(crate) fn alt_ha_dist(&self, eph: &Ephemeris, body: Body) -> Result<(f64, f64, f64), CalcError> {
        let (src, opp) = source(body).ok_or(CalcError::UnsupportedBody(body))?;
        let (u, dist) = self.direction(eph, src, opp)?;
        let (ra, dec) = to_spherical(u);
        let ha = (self.eo.gast + self.obs_lng - ra).rem_euclid(std::f64::consts::TAU);
        let (_, alt) = equatorial_to_horizontal(ha, dec, self.obs_lat);
        Ok((alt, ha, dist))
    }

    /// Spherical coordinates in degrees plus distance.
    fn coords(&self, eph: &Ephemeris, src: Source, opposite: bool, cs: CoordSystem) -> Result<(f64, f64, f64), CalcError> {
        let (u, dist) = self.direction(eph, src, opposite)?;
        let (lng, lat) = match cs {
            CoordSystem::Equatorial => to_spherical(u),
            CoordSystem::Ecliptic => {
                let (s, c) = self.eo.eps_true.sin_cos();
                to_spherical([u[0], c * u[1] + s * u[2], -s * u[1] + c * u[2]])
            }
            CoordSystem::Horizontal => {
                let (ra, dec) = to_spherical(u);
                let (az, alt) = equatorial_to_horizontal(self.eo.gast + self.obs_lng - ra, dec, self.obs_lat);
                (az, alt)
            }
        };
        Ok((lng.to_degrees(), lat.to_degrees(), dist))
    }
}

/// Hour angle and declination → (azimuth, altitude), radians. Azimuth is measured from the
/// south point westwards (0 = south, 90° = west), the convention of the reference API;
/// altitude is geometric (no refraction).
pub fn equatorial_to_horizontal(hour_angle: f64, dec: f64, lat: f64) -> (f64, f64) {
    let (sh, ch) = hour_angle.sin_cos();
    let (sd, cd) = dec.sin_cos();
    let (sp, cp) = lat.sin_cos();
    let az = sh.atan2(ch * sp - sd / cd * cp).rem_euclid(std::f64::consts::TAU);
    let alt = (sp * sd + cp * cd * ch).clamp(-1.0, 1.0).asin();
    (az, alt)
}

/// Where a body's position comes from.
#[derive(Debug, Clone, Copy, PartialEq)]
enum Source {
    Naif(i32),
    MeanNode,
    TrueNode,
    Apsis(Apsis),
    Kepler(crate::hypothetical::OrbitalElements),
    Star(crate::stars::Star),
}

/// The lunar apsides offered as bodies.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Apsis {
    MeanApogee,
    OscuApogee,
    IntpApogee,
    IntpPerigee,
}

/// Data source for a body and whether to take the diametrically opposite point.
fn source(body: Body) -> Option<(Source, bool)> {
    Some(match body {
        Body::MeanNode => (Source::MeanNode, false),
        Body::MeanSouthNode => (Source::MeanNode, true),
        Body::TrueNode => (Source::TrueNode, false),
        Body::SouthNode => (Source::TrueNode, true),
        Body::MeanApogee => (Source::Apsis(Apsis::MeanApogee), false),
        Body::OscuApogee => (Source::Apsis(Apsis::OscuApogee), false),
        Body::IntpApogee => (Source::Apsis(Apsis::IntpApogee), false),
        Body::IntpPerigee => (Source::Apsis(Apsis::IntpPerigee), false),
        b => match (crate::hypothetical::elements(b), crate::stars::star_of(b)) {
            (Some(el), _) => (Source::Kepler(el), false),
            (_, Some(star)) => (Source::Star(*star), false),
            _ => (Source::Naif(naif_id(b)?), false),
        },
    })
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
        let (src, opp) = source(body).ok_or(CalcError::UnsupportedBody(body))?;
        self.source_pos(eph, body, src, opp, cs)
    }

    fn source_pos(&self, eph: &Ephemeris, body: Body, src: Source, opp: bool, cs: CoordSystem) -> Result<BodyPos, CalcError> {
        let [a, m, b] = &self.0;
        let (l0, b0, d0) = a.coords(eph, src, opp, cs)?;
        let (l1, b1, d1) = m.coords(eph, src, opp, cs)?;
        let (l2, b2, d2) = b.coords(eph, src, opp, cs)?;
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
        let frames = FrameTriple::new(self, jd.into().value(), &geo, perspective)?;
        frames.body_pos(self, body, cs)
    }

    /// Apparent position of a star from catalogue data supplied by the caller, for stars
    /// that are not in the built-in table. The data must be ICRS at epoch J1991.25, as in
    /// the Hipparcos catalogue.
    pub fn calc_star_position(
        &self,
        star: &crate::stars::Star,
        jd: impl Into<JdUt>,
        geo: GeoPos,
        cs: CoordSystem,
        perspective: Perspective,
    ) -> Result<BodyPos, CalcError> {
        let frames = FrameTriple::new(self, jd.into().value(), &geo, perspective)?;
        frames.source_pos(self, star.body(), Source::Star(*star), false, cs)
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
