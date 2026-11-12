//! Rising, setting and meridian transits.
//!
//! Events are found by sampling the body's topocentric geometric altitude (or hour angle) and
//! refining each bracketed crossing by bisection to about 1 ms. A body rises or sets when its
//! centre's geometric altitude equals
//!
//! `h₀ = −refraction ± semi-diameter` (upper limb −, lower limb +, centre 0),
//!
//! where refraction at the horizon is [`HORIZON_REFRACTION_ARCMIN`] when enabled.
//! Constants were calibrated against the AstroCalc reference API (all eight legacy modes agree
//! to ~0.03′ of altitude, i.e. a fraction of a second in time).

use crate::body::Body;
use crate::error::CalcError;
use crate::geo::GeoPos;
use crate::position::{Frame, Perspective};
use crate::spk::{AU_KM, Ephemeris};
use crate::time::JdUt;

/// Refraction at the horizon used for rise/set, in arcminutes.
pub const HORIZON_REFRACTION_ARCMIN: f64 = 36.739;
/// Semi-diameters used when [`RiseSetOptions::fixed_disc`] is set (arcminutes).
const FIXED_SD_SUN: f64 = 15.994;
const FIXED_SD_MOON: f64 = 15.542;
/// Sampling step for event search (days).
const STEP: f64 = 1.0 / 48.0;
/// Bisection tolerance (days).
const TOL: f64 = 1e-8;

/// Which point of the disc defines rising and setting.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum DiscAlignment {
    /// Upper limb touches the horizon (conventional sunrise/sunset).
    UpperLimb,
    /// Centre of the disc.
    #[default]
    Center,
    /// Lower limb.
    LowerLimb,
}

/// Rise/set options.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct RiseSetOptions {
    pub disc: DiscAlignment,
    /// Apply atmospheric refraction at the horizon.
    pub refraction: bool,
    /// Use a fixed mean semi-diameter for the Sun and Moon instead of the actual one.
    pub fixed_disc: bool,
}

impl Default for RiseSetOptions {
    /// Centre of the disc with refraction (legacy mode 3).
    fn default() -> Self {
        Self::from_legacy_mode(3)
    }
}

impl RiseSetOptions {
    /// Maps the AstroCalc API `mode` parameter (0–7), including its historical behaviour:
    /// mode 1 ("no refraction only") behaves like mode 3, and the "fixed disc" modes use the
    /// upper limb.
    pub const fn from_legacy_mode(mode: u8) -> Self {
        let (disc, refraction, fixed_disc) = match mode {
            1 | 3 => (DiscAlignment::Center, true, false),
            2 => (DiscAlignment::Center, false, false),
            4 => (DiscAlignment::LowerLimb, false, false),
            5 => (DiscAlignment::LowerLimb, true, false),
            6 => (DiscAlignment::UpperLimb, false, true),
            7 => (DiscAlignment::UpperLimb, true, true),
            _ => (DiscAlignment::UpperLimb, true, false),
        };
        Self { disc, refraction, fixed_disc }
    }
}

/// Kinds of horizon and meridian events.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum TransitEvent {
    Rise,
    Set,
    /// Upper meridian transit (culmination).
    Mc,
    /// Lower meridian transit.
    Ic,
}

/// Rise, set and transits for the local day containing a given time, matching the
/// `sunRiseSets` / `moonRiseSets` objects of the AstroCalc API. All searches start at local
/// mean midnight before the reference time; `None` means no such event within the search span
/// (e.g. circumpolar bodies).
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct DayTransits {
    pub body: Body,
    /// Start of the local day (local mean midnight).
    pub day_start: JdUt,
    /// First setting after the previous local midnight.
    pub prev_set: Option<JdUt>,
    pub rise: Option<JdUt>,
    pub mc: Option<JdUt>,
    pub set: Option<JdUt>,
    pub ic: Option<JdUt>,
    /// First rising after `set`.
    pub next_rise: Option<JdUt>,
    /// Geometric altitude at `ic` and `mc` (degrees).
    pub min_alt: Option<f64>,
    pub max_alt: Option<f64>,
}

/// Mean radius in km used for semi-diameters.
pub(crate) fn radius_km(body: Body) -> f64 {
    match body {
        Body::Sun => 696_000.0,
        Body::Moon => 1_738.0,
        Body::Mercury => 2_439.7,
        Body::Venus => 6_051.8,
        Body::Mars => 3_396.2,
        Body::Jupiter => 71_492.0,
        Body::Saturn => 60_268.0,
        Body::Uranus => 25_559.0,
        Body::Neptune => 24_764.0,
        Body::Pluto => 1_188.3,
        // Mean radii from the JPL Small-Body Database.
        Body::Ceres => 469.7,
        Body::Pallas => 256.5,
        Body::Juno => 123.3,
        Body::Vesta => 261.4,
        Body::Chiron => 83.0,
        Body::Pholus => 95.0,
        _ => 0.0,
    }
}

/// Local mean midnight starting the local day that contains `jd`.
pub fn local_day_start(jd: f64, lng: f64) -> f64 {
    let offset = lng / 360.0;
    (jd + 0.5 + offset).floor() - 0.5 - offset
}

impl Ephemeris {
    fn horizon_target(&self, body: Body, opts: &RiseSetOptions, dist_au: f64) -> f64 {
        let refr = if opts.refraction { HORIZON_REFRACTION_ARCMIN } else { 0.0 };
        let sd = match (opts.fixed_disc, body) {
            (true, Body::Sun) => FIXED_SD_SUN,
            (true, Body::Moon) => FIXED_SD_MOON,
            _ if dist_au > 0.0 => (radius_km(body) / (dist_au * AU_KM)).asin().to_degrees() * 60.0,
            _ => 0.0,
        };
        let limb = match opts.disc {
            DiscAlignment::UpperLimb => -sd,
            DiscAlignment::Center => 0.0,
            DiscAlignment::LowerLimb => sd,
        };
        (limb - refr) / 60.0
    }

    /// Signed function whose zero crossings in the right direction are the event.
    fn event_fn(&self, body: Body, geo: &GeoPos, event: TransitEvent, opts: &RiseSetOptions, t: f64) -> Result<f64, CalcError> {
        let frame = Frame::new(self, t, geo, Perspective::Topocentric)?;
        let (alt, ha, dist) = frame.alt_ha_dist(self, body)?;
        Ok(match event {
            TransitEvent::Rise | TransitEvent::Set => alt.to_degrees() - self.horizon_target(body, opts, dist),
            // sin(HA) goes − → + at upper transit (HA = 0) and + → − at lower transit (HA = 180°)
            TransitEvent::Mc | TransitEvent::Ic => ha.sin(),
        })
    }

    /// First occurrence of `event` strictly after `start`, searching up to `max_days` ahead.
    pub fn next_transit(
        &self,
        body: Body,
        start: impl Into<JdUt>,
        geo: GeoPos,
        event: TransitEvent,
        opts: RiseSetOptions,
        max_days: f64,
    ) -> Result<Option<JdUt>, CalcError> {
        let start = start.into().value();
        let rising = matches!(event, TransitEvent::Rise | TransitEvent::Mc);
        let crosses = |a: f64, b: f64| if rising { a < 0.0 && b >= 0.0 } else { a > 0.0 && b <= 0.0 };
        let mut t0 = start;
        let mut f0 = self.event_fn(body, &geo, event, &opts, t0)?;
        while t0 < start + max_days {
            let t1 = t0 + STEP;
            let f1 = self.event_fn(body, &geo, event, &opts, t1)?;
            // For transits, ignore the sin(HA) sign change at the opposite meridian.
            let valid = match event {
                TransitEvent::Mc | TransitEvent::Ic => (f1 - f0).abs() < 1.0,
                _ => true,
            };
            if valid && crosses(f0, f1) {
                let (mut a, mut b, mut fa) = (t0, t1, f0);
                while b - a > TOL {
                    let m = 0.5 * (a + b);
                    let fm = self.event_fn(body, &geo, event, &opts, m)?;
                    if crosses(fa, fm) {
                        b = m;
                    } else {
                        a = m;
                        fa = fm;
                    }
                }
                return Ok(Some(JdUt(0.5 * (a + b))));
            }
            t0 = t1;
            f0 = f1;
        }
        Ok(None)
    }

    /// Geometric topocentric altitude (degrees) of a body's centre.
    pub fn altitude(&self, body: Body, jd: impl Into<JdUt>, geo: GeoPos) -> Result<f64, CalcError> {
        let frame = Frame::new(self, jd.into().value(), &geo, Perspective::Topocentric)?;
        Ok(frame.alt_ha_dist(self, body)?.0.to_degrees())
    }

    /// Rise, set and transits for the local day containing `jd` (see [`DayTransits`]).
    pub fn day_transits(&self, body: Body, jd: impl Into<JdUt>, geo: GeoPos, opts: RiseSetOptions) -> Result<DayTransits, CalcError> {
        const SPAN: f64 = 2.0;
        let start = local_day_start(jd.into().value(), geo.lng);
        let next = |from: f64, ev| self.next_transit(body, from, geo, ev, opts, SPAN);
        let rise = next(start, TransitEvent::Rise)?;
        let set = next(start, TransitEvent::Set)?;
        let mc = next(start, TransitEvent::Mc)?;
        let ic = next(start, TransitEvent::Ic)?;
        let prev_set = next(start - 1.0, TransitEvent::Set)?;
        let next_rise = match set {
            Some(s) => next(s.value(), TransitEvent::Rise)?,
            None => None,
        };
        let alt = |t: Option<JdUt>| t.map(|t| self.altitude(body, t, geo)).transpose();
        Ok(DayTransits {
            body,
            day_start: JdUt(start),
            prev_set,
            rise,
            mc,
            set,
            ic,
            next_rise,
            min_alt: alt(ic)?,
            max_alt: alt(mc)?,
        })
    }
}
