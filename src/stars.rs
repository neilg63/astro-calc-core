//! Fixed stars.
//!
//! The built-in table holds the 410 stars that have a name approved by the IAU and a
//! Hipparcos number. Each is addressed as a body through its Hipparcos number
//! (`Body::Catalog { catalog: Catalog::Hipparcos, id }`, written `hip:65474`), so stars work
//! with the same position, rise/set and aspect methods as the planets. Other stars can be
//! supplied as a [`Star`] to [`Ephemeris::calc_star_position`](crate::Ephemeris::calc_star_position).
//!
//! Apparent places include proper motion from the catalogue epoch, annual parallax, light
//! deflection, aberration, precession and nutation. Radial velocity is not modelled.

use crate::body::{Body, Catalog};
use crate::star_table::STARS;

/// Julian Day of the Hipparcos catalogue epoch, J1991.25.
pub(crate) const HIPPARCOS_EPOCH_JD: f64 = 2_448_349.062_5;

/// Catalogue data for one star: ICRS position at epoch J1991.25.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Star {
    /// Hipparcos catalogue number.
    pub hip: u32,
    /// IAU-approved proper name; empty for a star supplied by the caller.
    pub name: &'static str,
    /// Right ascension, degrees.
    pub ra: f64,
    /// Declination, degrees.
    pub dec: f64,
    /// Proper motion in right ascension, μα·cos δ, milliarcseconds per year.
    pub pm_ra: f64,
    /// Proper motion in declination, milliarcseconds per year.
    pub pm_dec: f64,
    /// Parallax, milliarcseconds (zero or negative if unknown).
    pub parallax: f64,
    /// Visual magnitude (NaN if unknown).
    pub magnitude: f64,
}

impl Star {
    /// The star as a body, for the position, rise/set and aspect methods.
    pub const fn body(&self) -> Body {
        Body::Catalog { catalog: Catalog::Hipparcos, id: self.hip as u64 }
    }

    /// Distance in AU from the parallax, or 0 if the parallax is not positive.
    pub fn distance_au(&self) -> f64 {
        if self.parallax > 0.0 { 206_264.806_247 * 1_000.0 / self.parallax } else { 0.0 }
    }
}

/// The built-in named stars, ordered by Hipparcos number.
pub fn named_stars() -> &'static [Star] {
    &STARS
}

/// A built-in star by Hipparcos number.
pub fn star_by_hip(hip: u32) -> Option<&'static Star> {
    STARS.binary_search_by_key(&hip, |s| s.hip).ok().map(|i| &STARS[i])
}

/// A built-in star by its IAU name. Case, spaces, hyphens and apostrophes are ignored, so
/// `"barnards star"` finds Barnard's Star.
pub fn star_by_name(name: &str) -> Option<&'static Star> {
    let fold = |s: &str| s.chars().filter(|c| c.is_alphanumeric()).flat_map(char::to_lowercase).collect::<String>();
    let wanted = fold(name);
    if wanted.is_empty() {
        return None;
    }
    STARS.iter().find(|s| fold(s.name) == wanted)
}

/// The built-in star a body refers to, if it is a Hipparcos catalogue body in the table.
pub(crate) fn star_of(body: Body) -> Option<&'static Star> {
    match body {
        Body::Catalog { catalog: Catalog::Hipparcos, id } => star_by_hip(u32::try_from(id).ok()?),
        _ => None,
    }
}
