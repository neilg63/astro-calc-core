//! Planetary phenomena: phase angle, illuminated fraction, elongation, apparent diameter and
//! apparent magnitude.
//!
//! Geometry uses the apparent geocentric direction of the body and of the Sun and the body's
//! resulting Sun–body–Earth triangle. Magnitudes follow the *Astronomical Almanac* phase laws:
//! Hilton (2005) for Mercury and Venus, Harris (1961) for Mars–Pluto, Allen's law for the Moon,
//! and the ring term of Meeus, *Astronomical Algorithms* ch. 45, for Saturn. They reproduce the
//! AstroCalc reference API to 0.001 mag (Saturn to about 0.03 mag). Minor bodies use the IAU
//! H–G system with parameters from the JPL Small-Body Database.
//!
//! Disc sizes use current IAU radii, so they differ from the reference API where its radii
//! are outdated (most visibly Pluto).

use crate::body::Body;
use crate::error::CalcError;
use crate::geo::GeoPos;
use crate::math::{dot, norm};
use crate::position::{Frame, Perspective};
use crate::rise_set::radius_km;
use crate::spk::{AU_KM, Ephemeris};
use crate::time::JdUt;

/// Phenomena of one body at one instant. Angles in degrees.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Phenomena {
    pub body: Body,
    /// Sun–body–Earth angle.
    pub phase_angle: f64,
    /// Illuminated fraction of the disc, 0–1.
    pub illuminated: f64,
    /// Angular distance of the body from the Sun.
    pub elongation: f64,
    /// Apparent diameter of the disc.
    pub diameter: f64,
    /// Apparent visual magnitude.
    pub magnitude: f64,
}

impl Ephemeris {
    /// Geocentric phenomena of a body.
    pub fn phenomena(&self, body: Body, jd: impl Into<JdUt>) -> Result<Phenomena, CalcError> {
        let frame = Frame::new(self, jd.into().value(), &GeoPos::default(), Perspective::Geocentric)?;
        let g = frame.geometry(self, body)?;
        let diameter = if g.dist > 0.0 { 2.0 * (radius_km(body) / (g.dist * AU_KM)).asin().to_degrees() } else { 0.0 };
        if body == Body::Sun {
            return Ok(Phenomena {
                body,
                phase_angle: 0.0,
                illuminated: 0.0,
                elongation: 0.0,
                diameter,
                magnitude: -26.86 + 5.0 * g.dist.log10(),
            });
        }
        let alpha = g.phase_angle;
        Ok(Phenomena {
            body,
            phase_angle: alpha,
            illuminated: (1.0 + alpha.to_radians().cos()) / 2.0,
            elongation: g.elongation,
            diameter,
            magnitude: magnitude(body, alpha, g.helio_dist, g.dist, g.ring_tilt_sin),
        })
    }

    /// Phenomena of several bodies at one instant.
    pub fn phenomena_set(&self, bodies: &[Body], jd: impl Into<JdUt>) -> Result<Vec<Phenomena>, CalcError> {
        let jd = jd.into();
        bodies.iter().map(|&b| self.phenomena(b, jd)).collect()
    }
}

/// Geometry needed for phenomena.
pub(crate) struct Geometry {
    pub phase_angle: f64,
    pub elongation: f64,
    /// Earth–body distance (AU).
    pub dist: f64,
    /// Sun–body distance (AU).
    pub helio_dist: f64,
    /// sin of the Saturnicentric latitude of the Earth referred to the ring plane (Saturn only).
    pub ring_tilt_sin: f64,
}

pub(crate) fn angle_between(a: [f64; 3], b: [f64; 3]) -> f64 {
    (dot(a, b) / (norm(a) * norm(b))).clamp(-1.0, 1.0).acos().to_degrees()
}

/// Apparent visual magnitude from phase angle α (degrees), heliocentric distance r and
/// geocentric distance Δ (AU).
fn magnitude(body: Body, alpha: f64, r: f64, delta: f64, ring_tilt_sin: f64) -> f64 {
    let dist_term = 5.0 * (r * delta).log10();
    let a = alpha / 100.0;
    match body {
        Body::Moon => {
            // Allen: 0.026|α| + 4·10⁻⁹ α⁴ phase law, scaled with the geocentric distance only
            // (zero point calibrated against the reference API).
            0.2312 + 5.0 * delta.log10() + 0.026 * alpha.abs() + 4e-9 * alpha.powi(4)
        }
        // Hilton (2005)
        Body::Mercury => -0.60 + dist_term + 4.98 * a - 4.88 * a * a + 3.02 * a * a * a,
        Body::Venus if alpha < 163.6 => -4.47 + dist_term + 1.03 * a + 0.57 * a * a + 0.13 * a * a * a,
        Body::Venus => 0.98 + dist_term - 1.02 * a,
        Body::Mars => -1.52 + dist_term + 0.016 * alpha,
        Body::Jupiter => -9.40 + dist_term + 0.005 * alpha,
        Body::Saturn => {
            let s = ring_tilt_sin.abs();
            -8.88 + dist_term + 0.044 * alpha - 2.60 * s + 1.25 * s * s
        }
        Body::Uranus => -7.19 + dist_term,
        Body::Neptune => -6.87 + dist_term,
        Body::Pluto => -1.00 + dist_term,
        other => match hg_parameters(other) {
            Some((h, g)) => hg_magnitude(h, g, alpha, r, delta),
            None => f64::NAN,
        },
    }
}

/// Absolute magnitude H and slope parameter G of the built-in minor bodies, from the JPL
/// Small-Body Database. Chiron and Pholus have no published G; the standard 0.15 is assumed.
pub const fn hg_parameters(body: Body) -> Option<(f64, f64)> {
    Some(match body {
        Body::Ceres => (3.34, 0.12),
        Body::Pallas => (4.12, 0.11),
        Body::Juno => (5.19, 0.32),
        Body::Vesta => (3.25, 0.32),
        Body::Chiron => (5.54, 0.15),
        Body::Pholus => (7.07, 0.15),
        _ => return None,
    })
}

/// Apparent visual magnitude of a minor body in the IAU H–G system (Bowell et al. 1989):
/// absolute magnitude `h`, slope parameter `g`, phase angle in degrees, and the distances
/// from the Sun and from the observer in AU. Use it with catalogue values for bodies that
/// are not built in.
pub fn hg_magnitude(h: f64, g: f64, phase_angle: f64, helio_dist: f64, dist: f64) -> f64 {
    let t = (phase_angle.to_radians() / 2.0).tan().abs();
    let phi1 = (-3.33 * t.powf(0.63)).exp();
    let phi2 = (-1.87 * t.powf(1.22)).exp();
    h + 5.0 * (helio_dist * dist).log10() - 2.5 * ((1.0 - g) * phi1 + g * phi2).log10()
}
