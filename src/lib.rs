//! # astro_calc_core
//!
//! Pure-Rust ephemeris core, built from published standards and JPL data files.
//!
//! **Cleanroom rule:** this is an independent implementation. Nothing in this crate may be
//! derived from the source code of any other ephemeris library. Sources are public
//! specifications only (NAIF SPK documentation, IAU SOFA/ERFA documentation, the Explanatory
//! Supplement, IERS conventions) and published results used as a behavioural reference.
//! Cite the reference in each algorithm's doc comment.
//!
//! Design rules:
//! * All times are Julian Days as `f64`, typed as [`JdUt`] or [`JdTt`].
//! * No global mutable state: every option is passed per call, so any number of threads can
//!   calculate at once.
//! * Synchronous and CPU-bound. Callable from async code: short calls inline, bulk calls via
//!   `spawn_blocking` (optional `tokio` wrappers to follow behind a feature flag).
//! * Builds for WebAssembly: turn off the default `mmap` feature and load data with
//!   [`Ephemeris::from_bytes`].

pub mod aspects;
pub mod ayanamsha;
pub mod body;
pub mod derived;
pub mod error;
pub mod events;
pub mod geo;
pub mod houses;
pub mod hypothetical;
pub mod indian_time;
pub mod math;
pub mod models;
pub mod pheno;
pub mod position;
pub mod rise_set;
pub mod spk;
mod star_table;
pub mod stars;
pub mod station_table;
pub mod time;

pub use aspects::{
    AspectMatch, AspectOrb, AspectPeriod, BodyAspect, ChartKind, ChartSubject, ProgressionOptions, TROPICAL_YEAR, aspect_orbs,
    aspects_between, best_aspect, match_aspects, separation,
};
pub use ayanamsha::{Ayanamsha, to_sidereal};
pub use body::{Body, BodyKind, Catalog, parse_body_list};
pub use derived::{Angle, Derived, DerivedId, DerivedRegistry, Ref};
pub use error::{BodyParseError, CalcError, DerivedError, EphemError, GeoError, TimeError};
pub use events::{JdSpeed, MoonPhase, MoonPhaseEvent, Station, StationKind, StationSet};
pub use geo::GeoPos;
pub use hypothetical::{Equinox, OrbitalElements};
pub use houses::{Angles, HouseSystem, Houses, angles, ascendant, house_sets, houses};
pub use indian_time::{IndianDay, SolarDayKind};
pub use models::DeltaTModel;
pub use pheno::{Phenomena, hg_magnitude, hg_parameters};
pub use position::{BodyJdPos, BodyPos, BodyPosSet, CoordSystem, Perspective};
pub use rise_set::{DayTransits, DiscAlignment, RiseSetOptions, TransitEvent, local_day_start};
pub use spk::Ephemeris;
pub use stars::{Star, named_stars, star_by_hip, star_by_name};
pub use station_table::StationTable;
pub use time::{J2000, JdTt, JdUt};

#[cfg(test)]
mod tests;
