//! The Indian solar day and its traditional time units.
//!
//! An Indian day (*ahorātra*) runs from one sunrise to the next. It divides into
//! 30 *muhūrtas* or 60 *ghaṭīs* (≈ 24 minutes); a ghaṭī into 60 *vighaṭīs* (≈ 24 seconds);
//! a vighaṭī into 60 *liptās* (≈ 0.4 seconds). Because the day is measured between actual
//! sunrises, the units stretch and shrink slightly through the year.

use std::fmt;

use julian_day_converter::julian_day_to_weekday_index;

use crate::body::Body;
use crate::error::CalcError;
use crate::geo::GeoPos;
use crate::rise_set::{RiseSetOptions, TransitEvent};
use crate::spk::Ephemeris;
use crate::time::JdUt;

/// How the day's boundaries were determined.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum SolarDayKind {
    /// Ordinary day: the Sun rises and sets.
    #[default]
    Normal,
    /// Midnight sun: no sunrise or sunset. The day runs between lower meridian transits and
    /// counts entirely as daytime.
    PolarDay,
    /// Polar night: no sunrise or sunset. The day runs between upper meridian transits and
    /// counts entirely as night.
    PolarNight,
}

/// One Indian solar day around a moment in time. All times are UT Julian Days.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct IndianDay {
    /// The moment of interest.
    pub jd: f64,
    /// Sunrise that began the day (the last sunrise at or before `jd`).
    pub rise_jd: f64,
    /// Sunset of the day.
    pub set_jd: f64,
    /// Sunrise that ends the day.
    pub next_rise_jd: f64,
    pub kind: SolarDayKind,
}

impl IndianDay {
    /// Builds a day from known sunrise, sunset and following sunrise.
    pub const fn new(jd: f64, rise_jd: f64, set_jd: f64, next_rise_jd: f64) -> Self {
        Self { jd, rise_jd, set_jd, next_rise_jd, kind: SolarDayKind::Normal }
    }

    /// Length of the whole day, sunrise to sunrise, in (civil) days.
    pub fn day_length(&self) -> f64 {
        self.next_rise_jd - self.rise_jd
    }

    /// Length of daylight, sunrise to sunset, in days.
    pub fn daylight_length(&self) -> f64 {
        self.set_jd - self.rise_jd
    }

    /// Length of the night, sunset to next sunrise, in days.
    pub fn night_length(&self) -> f64 {
        self.next_rise_jd - self.set_jd
    }

    /// Fraction of the day elapsed at `jd`, 0 ≤ progress < 1.
    pub fn progress(&self) -> f64 {
        ((self.jd - self.rise_jd) / self.day_length()).clamp(0.0, 1.0 - f64::EPSILON)
    }

    pub fn is_day_time(&self) -> bool {
        self.jd < self.set_jd
    }

    /// Muhūrta of the day, 0–29 (1/30 of the day, ≈ 48 minutes).
    pub fn muhurta(&self) -> u8 {
        (self.progress() * 30.0).floor() as u8
    }

    /// Ghaṭī of the day, 0–59 (1/60 of the day, ≈ 24 minutes).
    pub fn ghati(&self) -> u8 {
        (self.progress() * 60.0).floor() as u8
    }

    /// Vighaṭī within the ghaṭī, 0–59 (≈ 24 seconds).
    pub fn vighati(&self) -> u8 {
        ((self.progress() * 3_600.0).floor() % 60.0) as u8
    }

    /// Liptā within the vighaṭī, 0 ≤ liptā < 60 (≈ 0.4 seconds), fractional.
    pub fn lipta(&self) -> f64 {
        (self.progress() * 216_000.0) % 60.0
    }

    /// Daytime muhūrta 0–14 counted from sunrise to sunset, or night muhūrta 0–14 counted
    /// from sunset to the next sunrise (each 1/15 of daylight or of the night).
    pub fn half_day_muhurta(&self) -> u8 {
        let frac = if self.is_day_time() {
            (self.jd - self.rise_jd) / self.daylight_length()
        } else {
            (self.jd - self.set_jd) / self.night_length()
        };
        ((frac * 15.0).floor() as u8).min(14)
    }

    /// Weekday of the Indian day (Sunday = 0 … Saturday = 6), taken from the local calendar
    /// date at sunrise. `offset_secs` is the local offset from UTC, e.g.
    /// [`GeoPos::solar_offset_secs`] for local mean time.
    pub fn weekday_index(&self, offset_secs: i32) -> u8 {
        julian_day_to_weekday_index(self.rise_jd, offset_secs)
    }

    /// `"GG:VV:LL.ll"`, e.g. `"41:49:57.11"`.
    pub fn to_ghati_string(&self) -> String {
        format!("{:02}:{:02}:{:05.2}", self.ghati(), self.vighati(), self.lipta())
    }

    /// `"41 gh 49 vi 57.11 li"`.
    pub fn to_unit_string(&self) -> String {
        format!("{} gh {} vi {:.2} li", self.ghati(), self.vighati(), self.lipta())
    }

    /// `"muhurta 20, ghati 41, vighati 49, lipta 57.11 (night)"`.
    pub fn to_long_string(&self) -> String {
        format!(
            "muhurta {}, ghati {}, vighati {}, lipta {:.2} ({})",
            self.muhurta(),
            self.ghati(),
            self.vighati(),
            self.lipta(),
            if self.is_day_time() { "day" } else { "night" }
        )
    }
}

/// Displays as [`IndianDay::to_ghati_string`].
impl fmt::Display for IndianDay {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.to_ghati_string())
    }
}

impl Ephemeris {
    /// Last occurrence of `event` at or before `jd` within `back` days, if any.
    fn last_sun_event(&self, jd: f64, geo: GeoPos, event: TransitEvent, opts: RiseSetOptions, back: f64) -> Result<Option<f64>, CalcError> {
        let mut found = None;
        let mut from = jd - back;
        while let Some(t) = self.next_transit(Body::Sun, from, geo, event, opts, jd - from)? {
            if t.value() > jd {
                break;
            }
            found = Some(t.value());
            from = t.value() + 0.01;
        }
        Ok(found)
    }

    /// The Indian solar day containing `jd` at a place: the sunrise before, the sunset and the
    /// next sunrise. Where the Sun does not rise or set, the day is bounded by meridian
    /// transits instead (see [`SolarDayKind`]).
    pub fn indian_day(&self, jd: impl Into<JdUt>, geo: GeoPos, opts: RiseSetOptions) -> Result<IndianDay, CalcError> {
        let jd = jd.into().value();
        let next = |from: f64, ev| self.next_transit(Body::Sun, from, geo, ev, opts, 1.6).map(|o| o.map(|t| t.value()));
        if let Some(rise) = self.last_sun_event(jd, geo, TransitEvent::Rise, opts, 1.2)?
            && let (Some(next_rise), Some(set)) = (next(rise + 0.01, TransitEvent::Rise)?, next(rise, TransitEvent::Set)?)
            && set < next_rise
        {
            return Ok(IndianDay { jd, rise_jd: rise, set_jd: set, next_rise_jd: next_rise, kind: SolarDayKind::Normal });
        }
        // No complete sunrise–sunset–sunrise cycle: the Sun stays above or below the horizon.
        let up = self.altitude(Body::Sun, jd, geo)? > 0.0;
        let event = if up { TransitEvent::Ic } else { TransitEvent::Mc };
        let start = self.last_sun_event(jd, geo, event, opts, 1.2)?.unwrap_or(jd);
        let end = next(start + 0.01, event)?.unwrap_or(start + 1.0);
        Ok(IndianDay {
            jd,
            rise_jd: start,
            set_jd: if up { end } else { start },
            next_rise_jd: end,
            kind: if up { SolarDayKind::PolarDay } else { SolarDayKind::PolarNight },
        })
    }

    /// Positions of bodies at the sunrise that began the Indian day containing `jd`.
    pub fn positions_at_sunrise(
        &self,
        bodies: &[Body],
        jd: impl Into<JdUt>,
        geo: GeoPos,
        cs: crate::position::CoordSystem,
        perspective: crate::position::Perspective,
        opts: RiseSetOptions,
    ) -> Result<(IndianDay, Vec<crate::position::BodyPos>), CalcError> {
        let day = self.indian_day(jd, geo, opts)?;
        let pos = self.calc_positions(bodies, day.rise_jd, geo, cs, perspective)?;
        Ok((day, pos))
    }

    /// Positions of bodies at the sunset of the Indian day containing `jd`.
    pub fn positions_at_sunset(
        &self,
        bodies: &[Body],
        jd: impl Into<JdUt>,
        geo: GeoPos,
        cs: crate::position::CoordSystem,
        perspective: crate::position::Perspective,
        opts: RiseSetOptions,
    ) -> Result<(IndianDay, Vec<crate::position::BodyPos>), CalcError> {
        let day = self.indian_day(jd, geo, opts)?;
        let pos = self.calc_positions(bodies, day.set_jd, geo, cs, perspective)?;
        Ok((day, pos))
    }
}
