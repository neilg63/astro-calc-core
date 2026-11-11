//! Observer location.

use std::fmt;
use std::str::FromStr;

use crate::error::GeoError;

/// Geographic position of an observer: degrees (north and east positive) and metres.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct GeoPos {
    pub lat: f64,
    pub lng: f64,
    pub alt: f64,
}

impl GeoPos {
    /// Creates a validated position.
    pub fn new(lat: f64, lng: f64, alt: f64) -> Result<Self, GeoError> {
        if !(-90.0..=90.0).contains(&lat) {
            return Err(GeoError::Latitude(lat));
        }
        if !(-180.0..=180.0).contains(&lng) {
            return Err(GeoError::Longitude(lng));
        }
        if !alt.is_finite() {
            return Err(GeoError::Altitude(alt));
        }
        Ok(Self { lat, lng, alt })
    }

    /// Local mean solar time offset from UT in seconds (1° = 240 s), as used for weekday
    /// calculations by `julian_day_converter`.
    pub fn solar_offset_secs(&self) -> i32 {
        (self.lng * 240.0).round() as i32
    }
}

/// Parses the API's `loc` format: `lat,lng` or `lat,lng,alt`.
impl FromStr for GeoPos {
    type Err = GeoError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let parts: Vec<f64> = s
            .split(',')
            .map(|p| p.trim().parse::<f64>())
            .collect::<Result<_, _>>()
            .map_err(|_| GeoError::Format(s.to_string()))?;
        match parts.as_slice() {
            [lat, lng] => Self::new(*lat, *lng, 0.0),
            [lat, lng, alt] => Self::new(*lat, *lng, *alt),
            _ => Err(GeoError::Format(s.to_string())),
        }
    }
}

impl fmt::Display for GeoPos {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{},{},{}", self.lat, self.lng, self.alt)
    }
}
