//! Error types.

use thiserror::Error;

/// Failure to parse a [`Body`](crate::Body) from a string.
#[derive(Debug, Clone, PartialEq, Eq, Error)]
pub enum BodyParseError {
    #[error("empty body key")]
    Empty,
    #[error("unknown body key: {0}")]
    Unknown(String),
    #[error("unknown catalogue prefix: {0}")]
    UnknownCatalog(String),
    #[error("invalid catalogue id in: {0}")]
    InvalidId(String),
}

/// Failure to parse or validate a [`GeoPos`](crate::GeoPos).
#[derive(Debug, Clone, PartialEq, Error)]
pub enum GeoError {
    #[error("expected 'lat,lng' or 'lat,lng,alt', got: {0}")]
    Format(String),
    #[error("latitude {0} out of range -90..=90")]
    Latitude(f64),
    #[error("longitude {0} out of range -180..=180")]
    Longitude(f64),
    #[error("altitude {0} is not finite")]
    Altitude(f64),
}

/// Failure to register a derived point.
#[derive(Debug, Clone, PartialEq, Eq, Error)]
pub enum DerivedError {
    #[error("name already registered: {0}")]
    DuplicateName(String),
    #[error("name clashes with a built-in body key: {0}")]
    ReservedName(String),
    #[error("reference to unregistered derived point {0}")]
    UnknownReference(u32),
}

/// Failure to convert between Julian Days and calendar date-times.
#[derive(Debug, Clone, PartialEq, Error)]
pub enum TimeError {
    #[error("Julian Day {0} is outside the supported calendar range")]
    OutOfRange(f64),
    #[error("cannot parse date-time: {0}")]
    Parse(String),
}

/// Failure to read or evaluate an SPK ephemeris.
#[derive(Debug, Clone, PartialEq, Error)]
pub enum EphemError {
    #[error("cannot open {0}: {1}")]
    Io(String, String),
    #[error("invalid SPK file {0}: {1}")]
    Format(String, String),
    #[error("no ephemeris data for NAIF id {target} at JD {jd}")]
    NoData { target: i32, jd: f64 },
}

/// Failure to calculate a position.
#[derive(Debug, Clone, PartialEq, Error)]
pub enum CalcError {
    #[error(transparent)]
    Ephem(#[from] EphemError),
    #[error("body not supported yet: {0}")]
    UnsupportedBody(crate::Body),
    #[error("invalid sampling: days {days}, interval {interval}")]
    InvalidInterval { days: f64, interval: f64 },
}
