//! Rising, setting and meridian transits.
//!
//! Events are found by sampling the body's topocentric geometric altitude (or hour angle) and
//! refining each bracketed crossing by bisection to about 1 ms. A body rises or sets when its
//! centre's geometric altitude equals
//!
//! `h₀ = horizon − refraction ± semi-diameter − dip` (upper limb −, lower limb +, centre 0),
//!
//! where refraction at the level horizon is [`HORIZON_REFRACTION_ARCMIN`] when enabled.
//! [`RiseSetOptions`] also covers the three twilights, a raised or lowered horizon,
//! refraction for a given pressure and temperature, the dip of the horizon for an elevated
//! observer, geocentric positions and bodies projected onto the ecliptic.
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

/// Atmospheric refraction at the horizon.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum Refraction {
    /// No atmosphere: the geometric horizon.
    Off,
    /// Standard conditions ([`HORIZON_REFRACTION_ARCMIN`] at the horizon).
    #[default]
    Standard,
    /// Scaled for the air at the observer: pressure in hPa (millibars) and temperature in °C.
    Atmosphere { pressure: f64, temperature: f64 },
}

/// The altitude a body must cross to count as rising or setting.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum Horizon {
    /// The level horizon, altitude 0°.
    #[default]
    Level,
    /// A horizon at an apparent altitude in degrees: a ridge above the level (positive) or
    /// a view downwards (negative). Disc alignment and refraction still apply.
    Altitude(f64),
    /// Civil twilight: the centre is 6° below the horizon. For the twilights the disc
    /// alignment, refraction and dip are not used.
    CivilTwilight,
    /// Nautical twilight: 12° below.
    NauticalTwilight,
    /// Astronomical twilight: 18° below.
    AstronomicalTwilight,
}

/// Rise/set options. `RiseSetOptions::default()` is the centre of the disc on the level
/// horizon with standard refraction; change fields directly or with the `with_` methods.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct RiseSetOptions {
    pub disc: DiscAlignment,
    pub refraction: Refraction,
    /// Use a fixed mean semi-diameter for the Sun and Moon instead of the actual one.
    pub fixed_disc: bool,
    pub horizon: Horizon,
    /// Lower the horizon by the dip due to the observer's height above the ground
    /// (the location's altitude, taken as height above the surrounding terrain or sea).
    pub dip: bool,
    /// Use the body's geocentric position, ignoring parallax.
    pub geocentric: bool,
    /// Place the body on the ecliptic at its longitude, ignoring its ecliptic latitude.
    pub ignore_latitude: bool,
}

impl Default for RiseSetOptions {
    /// Centre of the disc with refraction (legacy mode 3).
    fn default() -> Self {
        Self::from_legacy_mode(3)
    }
}

impl RiseSetOptions {
    /// The conventional almanac definition: upper limb on the level horizon with standard
    /// refraction.
    pub const fn upper_limb() -> Self {
        Self::from_legacy_mode(0)
    }

    /// The start or end of a twilight ([`Horizon::CivilTwilight`] and the others).
    pub const fn twilight(horizon: Horizon) -> Self {
        Self { horizon, ..Self::from_legacy_mode(3) }
    }

    /// Traditional Hindu rising: the centre of the disc on the geometric horizon, seen from
    /// the centre of the Earth, with the body placed on the ecliptic.
    pub const fn hindu_rising() -> Self {
        Self { refraction: Refraction::Off, geocentric: true, ignore_latitude: true, ..Self::from_legacy_mode(3) }
    }

    pub const fn with_disc(mut self, disc: DiscAlignment) -> Self {
        self.disc = disc;
        self
    }

    pub const fn with_refraction(mut self, refraction: Refraction) -> Self {
        self.refraction = refraction;
        self
    }

    pub const fn with_fixed_disc(mut self, fixed: bool) -> Self {
        self.fixed_disc = fixed;
        self
    }

    pub const fn with_horizon(mut self, horizon: Horizon) -> Self {
        self.horizon = horizon;
        self
    }

    pub const fn with_dip(mut self, dip: bool) -> Self {
        self.dip = dip;
        self
    }

    pub const fn with_geocentric(mut self, geocentric: bool) -> Self {
        self.geocentric = geocentric;
        self
    }

    pub const fn with_ignore_latitude(mut self, ignore: bool) -> Self {
        self.ignore_latitude = ignore;
        self
    }

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
        Self {
            disc,
            refraction: if refraction { Refraction::Standard } else { Refraction::Off },
            fixed_disc,
            horizon: Horizon::Level,
            dip: false,
            geocentric: false,
            ignore_latitude: false,
        }
    }
}

/// Refraction (arcminutes) at an apparent altitude (degrees), scaled so that the horizon
/// value is [`HORIZON_REFRACTION_ARCMIN`]. The variation with altitude is Bennett's formula.
fn refraction_at(apparent_alt: f64) -> f64 {
    let bennett = |h: f64| 1.0 / (h + 7.31 / (h + 4.4)).to_radians().tan();
    if apparent_alt == 0.0 {
        return HORIZON_REFRACTION_ARCMIN;
    }
    let h = apparent_alt.clamp(-1.0, 89.9);
    HORIZON_REFRACTION_ARCMIN * bennett(h) / bennett(0.0)
}

/// Dip of the horizon (arcminutes) for an observer `height` metres above it.
fn horizon_dip(height: f64, refracted: bool) -> f64 {
    const EARTH_RADIUS_M: f64 = 6_371_000.0;
    if height <= 0.0 {
        return 0.0;
    }
    if refracted {
        // The standard navigational allowance, which includes terrestrial refraction.
        1.76 * height.sqrt()
    } else {
        (EARTH_RADIUS_M / (EARTH_RADIUS_M + height)).acos().to_degrees() * 60.0
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
    /// Geometric altitude (degrees) of the body's centre at the moment of rising or setting.
    fn horizon_target(&self, body: Body, opts: &RiseSetOptions, dist_au: f64, height: f64) -> f64 {
        let base = match opts.horizon {
            Horizon::CivilTwilight => return -6.0,
            Horizon::NauticalTwilight => return -12.0,
            Horizon::AstronomicalTwilight => return -18.0,
            Horizon::Level => 0.0,
            Horizon::Altitude(alt) => alt,
        };
        let refr = match opts.refraction {
            Refraction::Off => 0.0,
            Refraction::Standard => refraction_at(base),
            // Meeus, Astronomical Algorithms ch. 16: proportional to pressure, inversely to
            // absolute temperature, relative to 1010 hPa and 10 °C.
            Refraction::Atmosphere { pressure, temperature } => refraction_at(base) * (pressure / 1010.0) * (283.0 / (273.0 + temperature)),
        };
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
        let dip = if opts.dip { horizon_dip(height, opts.refraction != Refraction::Off) } else { 0.0 };
        base + (limb - refr - dip) / 60.0
    }

    /// Signed function whose zero crossings in the right direction are the event.
    fn event_fn(&self, body: Body, geo: &GeoPos, event: TransitEvent, opts: &RiseSetOptions, t: f64) -> Result<f64, CalcError> {
        let perspective = if opts.geocentric { Perspective::Geocentric } else { Perspective::Topocentric };
        let frame = Frame::new(self, t, geo, perspective)?;
        let (alt, ha, dist) = frame.alt_ha_dist(self, body, opts.ignore_latitude)?;
        Ok(match event {
            TransitEvent::Rise | TransitEvent::Set => alt.to_degrees() - self.horizon_target(body, opts, dist, geo.alt),
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
        Ok(frame.alt_ha_dist(self, body, false)?.0.to_degrees())
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
