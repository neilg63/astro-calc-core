//! Hypothetical bodies defined by Keplerian orbital elements.
//!
//! The eight "Uranian" bodies of the Hamburg School follow the elements derived by
//! J. Neely (1988); Isis-Transpluto follows the elements of the 1946 ephemeris attributed to
//! Sevin and Landscheidt. None of these bodies exists; the elements are conventions.

use crate::body::Body;
use crate::math::{Mat3, Vec3, mat_mul, mat_t_vec, rot_x, rot_z};
use crate::models::precession_fw;

/// Gaussian gravitational constant: mean motion in radians per day is `GAUSS_K / a^1.5`.
const GAUSS_K: f64 = 0.017_202_098_95;
/// 1900 January 0.5.
const J1900: f64 = 2_415_020.0;

/// The reference direction the angular elements are measured from.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum Equinox {
    /// Mean ecliptic and equinox of a fixed Julian Day.
    Fixed(f64),
    /// Mean ecliptic and equinox of the date of calculation.
    OfDate,
}

/// A heliocentric Keplerian orbit. Angles in degrees, distance in AU.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct OrbitalElements {
    /// Julian Day (TT) at which `mean_anomaly` applies.
    pub epoch: f64,
    pub equinox: Equinox,
    pub mean_anomaly: f64,
    pub semi_major_axis: f64,
    pub eccentricity: f64,
    pub arg_perihelion: f64,
    pub node: f64,
    pub inclination: f64,
}

const fn uranian(m: f64, a: f64, e: f64, w: f64, node: f64, i: f64) -> OrbitalElements {
    OrbitalElements {
        epoch: J1900,
        equinox: Equinox::Fixed(J1900),
        mean_anomaly: m,
        semi_major_axis: a,
        eccentricity: e,
        arg_perihelion: w,
        node,
        inclination: i,
    }
}

/// Elements of the built-in hypothetical bodies.
pub const fn elements(body: Body) -> Option<OrbitalElements> {
    Some(match body {
        Body::Cupido => uranian(163.7409, 40.998_37, 0.004_60, 171.4333, 129.8325, 1.0833),
        Body::Hades => uranian(27.6496, 50.667_44, 0.002_45, 148.1796, 161.3339, 1.0500),
        Body::Zeus => uranian(165.1232, 59.214_36, 0.001_20, 299.0440, 0.0, 0.0),
        Body::Kronos => uranian(169.0193, 64.819_60, 0.003_05, 208.8801, 0.0, 0.0),
        Body::Apollon => uranian(138.0533, 70.299_49, 0.0, 0.0, 0.0, 0.0),
        Body::Admetos => uranian(351.3350, 73.627_65, 0.0, 0.0, 0.0, 0.0),
        Body::Vulkanus => uranian(55.8983, 77.255_68, 0.0, 0.0, 0.0, 0.0),
        Body::Poseidon => uranian(165.5163, 83.669_07, 0.0, 0.0, 0.0, 0.0),
        Body::Isis => OrbitalElements {
            epoch: 2_368_547.66,
            equinox: Equinox::Fixed(2_431_456.5),
            mean_anomaly: 0.0,
            semi_major_axis: 77.775,
            eccentricity: 0.3,
            arg_perihelion: 0.7,
            node: 0.0,
            inclination: 0.0,
        },
        _ => return None,
    })
}

/// Rotation from the GCRS to the mean ecliptic and equinox of a date (IAU 2006,
/// Fukushima–Williams angles).
pub(crate) fn gcrs_to_mean_ecliptic(jd_tt: f64) -> Mat3 {
    let (gamb, phib, psib, _) = precession_fw(jd_tt);
    mat_mul(&rot_z(-psib), &mat_mul(&rot_x(phib), &rot_z(gamb)))
}

impl OrbitalElements {
    /// Heliocentric position (AU) in the GCRS at a TT Julian Day.
    pub fn heliocentric(&self, jd_tt: f64) -> Vec3 {
        let a = self.semi_major_axis;
        let e = self.eccentricity;
        let m = (self.mean_anomaly.to_radians() + GAUSS_K / (a * a.sqrt()) * (jd_tt - self.epoch)).rem_euclid(std::f64::consts::TAU);
        // Kepler's equation by Newton's method.
        let mut ea = m + e * m.sin();
        for _ in 0..12 {
            let step = (ea - e * ea.sin() - m) / (1.0 - e * ea.cos());
            ea -= step;
            if step.abs() < 1e-14 {
                break;
            }
        }
        // Position in the orbital plane, x towards perihelion.
        let (x, y) = (a * (ea.cos() - e), a * (1.0 - e * e).sqrt() * ea.sin());
        let (sw, cw) = self.arg_perihelion.to_radians().sin_cos();
        let (sn, cn) = self.node.to_radians().sin_cos();
        let (si, ci) = self.inclination.to_radians().sin_cos();
        let (u, v) = (x * cw - y * sw, x * sw + y * cw);
        let ecl = [u * cn - v * ci * sn, u * sn + v * ci * cn, v * si];
        let frame = match self.equinox {
            Equinox::Fixed(jd) => jd,
            Equinox::OfDate => jd_tt,
        };
        mat_t_vec(&gcrs_to_mean_ecliptic(frame), ecl)
    }
}
