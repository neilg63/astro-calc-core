//! Julian Day time types.
//!
//! Internally the crate works only with Julian Days as `f64`. Two newtypes keep the time
//! scales apart at zero runtime cost:
//!
//! * [`JdUt`]: Universal Time. Every public calculation method takes this.
//! * [`JdTt`]: Terrestrial Time. Ephemeris lookups use this, via ΔT (not yet implemented).
//!
//! Calendar conversion happens only at the edges, through `julian_day_converter` and chrono,
//! at millisecond precision.

use chrono::{DateTime, NaiveDateTime, Utc};
use julian_day_converter::{
    JulianDay, julian_day_to_unix_millis, unix_millis_to_julian_day, unixtime_to_julian_day,
};

use crate::error::TimeError;

/// Julian Day in Universal Time.
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd, Default)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize), serde(transparent))]
pub struct JdUt(pub f64);

/// Julian Day in Terrestrial Time.
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd, Default)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize), serde(transparent))]
pub struct JdTt(pub f64);

/// J2000.0 epoch (2000-01-01T12:00 TT) as a Julian Day.
pub const J2000: f64 = 2_451_545.0;

impl JdUt {
    #[inline]
    pub const fn new(jd: f64) -> Self {
        Self(jd)
    }

    #[inline]
    pub const fn value(self) -> f64 {
        self.0
    }

    /// Adds a (possibly fractional or negative) number of days.
    #[inline]
    pub fn add_days(self, days: f64) -> Self {
        Self(self.0 + days)
    }

    /// Days elapsed from `earlier` to `self`.
    #[inline]
    pub fn days_since(self, earlier: JdUt) -> f64 {
        self.0 - earlier.0
    }

    pub fn from_naive(dt: &NaiveDateTime) -> Self {
        Self(dt.to_jd())
    }

    pub fn from_utc(dt: &DateTime<Utc>) -> Self {
        Self(dt.naive_utc().to_jd())
    }

    pub fn from_unix(secs: i64) -> Self {
        Self(unixtime_to_julian_day(secs))
    }

    pub fn from_unix_millis(millis: i64) -> Self {
        Self(unix_millis_to_julian_day(millis))
    }

    pub fn to_unix_millis(self) -> i64 {
        julian_day_to_unix_millis(self.0)
    }

    pub fn to_naive(self) -> Result<NaiveDateTime, TimeError> {
        NaiveDateTime::from_jd(self.0).ok_or(TimeError::OutOfRange(self.0))
    }

    pub fn to_utc(self) -> Result<DateTime<Utc>, TimeError> {
        self.to_naive().map(|n| n.and_utc())
    }
}

impl JdTt {
    #[inline]
    pub const fn new(jd: f64) -> Self {
        Self(jd)
    }

    #[inline]
    pub const fn value(self) -> f64 {
        self.0
    }

    /// Julian centuries of TT since J2000.0, the argument of most precession/nutation series.
    #[inline]
    pub fn centuries_since_j2000(self) -> f64 {
        (self.0 - J2000) / 36_525.0
    }
}

impl From<f64> for JdUt {
    fn from(jd: f64) -> Self {
        Self(jd)
    }
}

impl From<NaiveDateTime> for JdUt {
    fn from(dt: NaiveDateTime) -> Self {
        Self::from_naive(&dt)
    }
}

impl From<DateTime<Utc>> for JdUt {
    fn from(dt: DateTime<Utc>) -> Self {
        Self::from_utc(&dt)
    }
}
