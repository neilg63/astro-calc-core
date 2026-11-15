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
    /// The date is outside the data for a body, or no data for it is loaded. `covered` is
    /// the span the loaded files do hold for that body (Julian Days, TDB), if any.
    #[error("{}", no_data_message(*target, *jd, *covered))]
    NoData { target: i32, jd: f64, covered: Option<(f64, f64)> },
    /// The date is outside the JPL files and beyond the span of the analytical series too.
    /// `covered` is that span (Julian Days, TDB).
    #[error("{} (JD {jd}) is outside the JPL files and beyond the analytical series, which cover {} to {}", crate::spk::date_label(*jd), crate::spk::date_label(covered.0), crate::spk::date_label(covered.1))]
    BeyondAnalytic { jd: f64, covered: (f64, f64) },
    #[error("invalid station table: {0}")]
    Table(String),
}

fn no_data_message(target: i32, jd: f64, covered: Option<(f64, f64)>) -> String {
    let date = |jd: f64| crate::spk::date_label(jd);
    match covered {
        Some((start, end)) => format!(
            "{} (JD {jd}) is outside the JPL ephemeris data for NAIF id {target}, which covers {} to {}",
            date(jd),
            date(start),
            date(end)
        ),
        None => format!("no JPL ephemeris data is loaded for NAIF id {target}"),
    }
}

/// Failure to calculate a position.
#[derive(Debug, Clone, PartialEq, Error)]
pub enum CalcError {
    #[error(transparent)]
    Ephem(#[from] EphemError),
    #[error("body not supported yet: {0}")]
    UnsupportedBody(crate::Body),
    #[error("{body} cannot be observed from the {perspective:?} perspective")]
    NotVisibleFrom { body: crate::Body, perspective: crate::Perspective },
    #[error("horizontal coordinates need an observer on the Earth, not the {0:?} perspective")]
    NoHorizon(crate::Perspective),
    #[error("invalid sampling: days {days}, interval {interval}")]
    InvalidInterval { days: f64, interval: f64 },
}

impl CalcError {
    /// Whether the failure is a date outside the data the engine holds for a body. The span
    /// that is held is in [`EphemError::NoData`]'s `covered` field.
    pub fn is_out_of_range(&self) -> bool {
        matches!(self, CalcError::Ephem(EphemError::NoData { covered: Some(_), .. } | EphemError::BeyondAnalytic { .. }))
    }

    /// The engine that reported the failure, where one was involved.
    pub fn engine(&self) -> Option<crate::CalcEngine> {
        match self {
            CalcError::Ephem(EphemError::NoData { .. }) => Some(crate::CalcEngine::Jpl),
            CalcError::Ephem(EphemError::BeyondAnalytic { .. }) => Some(crate::CalcEngine::Analytic),
            _ => None,
        }
    }
}
