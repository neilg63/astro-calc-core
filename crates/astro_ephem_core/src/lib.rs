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

pub mod aspects;
pub mod ayanamsha;
pub mod body;
pub mod derived;
pub mod error;
pub mod events;
pub mod geo;
pub mod houses;
pub mod indian_time;
pub mod math;
pub mod models;
pub mod pheno;
pub mod position;
pub mod rise_set;
pub mod spk;
pub mod time;

pub use aspects::{
    AspectMatch, AspectOrb, AspectPeriod, BodyAspect, ChartKind, ChartSubject, ProgressionOptions, TROPICAL_YEAR, aspect_orbs,
    aspects_between, best_aspect, match_aspects, separation,
};
pub use ayanamsha::{Ayanamsha, to_sidereal};
pub use body::{Body, BodyKind, Catalog, parse_body_list};
pub use derived::{Angle, Derived, DerivedId, DerivedRegistry, Ref};
pub use error::{BodyParseError, CalcError, DerivedError, EphemError, GeoError, TimeError};
pub use events::{MoonPhase, MoonPhaseEvent, Station, StationKind};
pub use geo::GeoPos;
pub use houses::{Angles, HouseSystem, Houses, angles, ascendant, house_sets, houses};
pub use indian_time::{IndianDay, SolarDayKind};
pub use models::DeltaTModel;
pub use pheno::Phenomena;
pub use position::{BodyJdPos, BodyPos, BodyPosSet, CoordSystem, Perspective};
pub use rise_set::{DayTransits, DiscAlignment, RiseSetOptions, TransitEvent, local_day_start};
pub use spk::Ephemeris;
pub use time::{J2000, JdTt, JdUt};

#[cfg(test)]
mod tests;
