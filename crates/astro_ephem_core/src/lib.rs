//! # astro_ephem_core
//!
//! Pure-Rust ephemeris core, built from published standards and JPL data files.
//!
//! **Cleanroom rule:** nothing in this crate may be derived from Swiss Ephemeris source code.
//! Sources are public specifications only (NAIF SPK documentation, IAU SOFA/ERFA documentation,
//! the Explanatory Supplement, IERS conventions, published Swiss Ephemeris manuals used as a
//! behavioural spec). Cite the reference in each algorithm's doc comment.
//!
//! Design rules:
//! * All times are Julian Days as `f64`, typed as [`JdUt`] or [`JdTt`].
//! * No global mutable state: every option is passed per call, so any number of threads can
//!   calculate at once.
//! * Synchronous and CPU-bound. Callable from async code: short calls inline, bulk calls via
//!   `spawn_blocking` (optional `tokio` wrappers to follow behind a feature flag).

pub mod body;
pub mod derived;
pub mod error;
pub mod geo;
pub mod time;

pub use body::{Body, BodyKind, Catalog, parse_body_list};
pub use derived::{Angle, Derived, DerivedId, DerivedRegistry, Ref};
pub use error::{BodyParseError, DerivedError, GeoError, TimeError};
pub use geo::GeoPos;
pub use time::{J2000, JdTt, JdUt};

#[cfg(test)]
mod tests;
