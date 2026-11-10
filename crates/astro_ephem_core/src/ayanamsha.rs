//! Ayanamshas: the offset between the tropical and a sidereal zodiac.
//!
//! Two families, matching the reference API to within the stated tolerances (1850–2150):
//!
//! * **Epoch-based**: a fixed value at J2000 carried forward by the IAU 2006 general precession
//!   in longitude, plus nutation in longitude (a "true" ayanamsha, consistent with apparent
//!   tropical positions). Residual vs reference ≤ 0.005″ (Lahiri, Krishnamurti, Raman,
//!   Yukteshwar, JN Bhasin, Ushashashi); ≤ 0.12″ for the others.
//! * **Star-based**: the apparent true-ecliptic longitude of a reference star (or the Galactic
//!   Centre) minus its defined sidereal longitude. Star data: Hipparcos new reduction
//!   (van Leeuwen 2007); Galactic Centre: Sgr A* (Reid & Brunthaler 2004).
//!
//! J2000 values were fitted to the AstroCalc reference API (tools/harvest_houses_aya.py,
//! examples/aya_diag.rs and examples/aya_star.rs).

use std::fmt;
use std::str::FromStr;

use crate::error::CalcError;
use crate::geo::GeoPos;
use crate::math::{Vec3, norm_deg};
use crate::models::{general_precession_longitude, nutation};
use crate::position::{Frame, Perspective};
use crate::spk::Ephemeris;
use crate::time::JdUt;

/// Sidereal zodiac definitions.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum Ayanamsha {
    Tropical,
    Lahiri,
    Krishnamurti,
    Raman,
    Yukteshwar,
    JnBhasin,
    Ushashashi,
    Sassanian,
    Hipparchos,
    Aldebaran15Tau,
    ValensMoon,
    TrueCitra,
    TrueRevati,
    TruePushya,
    TrueMula,
    TrueSheoran,
    GalcentCochrane,
    GalcentMulaWilhelm,
}

/// Reference object for star-based ayanamshas: ICRS position at `epoch` (Julian year) and
/// proper motion (μα·cosδ, μδ) in mas/yr.
#[derive(Debug, Clone, Copy)]
struct StarRef {
    ra_deg: f64,
    dec_deg: f64,
    pm_ra: f64,
    pm_dec: f64,
    epoch: f64,
}

const fn hms(h: f64, m: f64, s: f64) -> f64 {
    (h + m / 60.0 + s / 3600.0) * 15.0
}

const fn dms(d: f64, m: f64, s: f64) -> f64 {
    if d < 0.0 { d - m / 60.0 - s / 3600.0 } else { d + m / 60.0 + s / 3600.0 }
}

/// α Vir (Spica), HIP 65474.
const SPICA: StarRef = StarRef { ra_deg: hms(13.0, 25.0, 11.6045), dec_deg: dms(-11.0, 9.0, 40.482), pm_ra: -42.35, pm_dec: -30.67, epoch: 1991.25 };
/// ζ Psc, HIP 5737.
const ZETA_PSC: StarRef = StarRef { ra_deg: hms(1.0, 13.0, 43.8020), dec_deg: dms(7.0, 34.0, 31.762), pm_ra: 145.00, pm_dec: -55.69, epoch: 1991.25 };
/// δ Cnc (Asellus Australis), HIP 42911.
const DELTA_CNC: StarRef = StarRef { ra_deg: hms(8.0, 44.0, 41.1101), dec_deg: dms(18.0, 9.0, 17.509), pm_ra: -17.67, pm_dec: -229.26, epoch: 1991.25 };
/// λ Sco (Shaula), HIP 85927.
const LAMBDA_SCO: StarRef = StarRef { ra_deg: hms(17.0, 33.0, 36.5264), dec_deg: dms(-37.0, 6.0, 13.495), pm_ra: -8.53, pm_dec: -30.80, epoch: 1991.25 };
/// Galactic Centre (Sgr A*).
const GAL_CENTER: StarRef = StarRef { ra_deg: hms(17.0, 45.0, 40.0409), dec_deg: dms(-29.0, 0.0, 28.118), pm_ra: 0.0, pm_dec: 0.0, epoch: 2000.0 };

enum Definition {
    None,
    /// Ayanamsha at J2000 (degrees, including nutation).
    Epoch(f64),
    /// Star whose apparent longitude is fixed at the given sidereal longitude (degrees).
    Star(StarRef, f64),
}

impl Ayanamsha {
    pub const ALL: [Ayanamsha; 18] = [
        Ayanamsha::Tropical,
        Ayanamsha::Lahiri,
        Ayanamsha::Krishnamurti,
        Ayanamsha::Raman,
        Ayanamsha::Yukteshwar,
        Ayanamsha::JnBhasin,
        Ayanamsha::Ushashashi,
        Ayanamsha::Sassanian,
        Ayanamsha::Hipparchos,
        Ayanamsha::Aldebaran15Tau,
        Ayanamsha::ValensMoon,
        Ayanamsha::TrueCitra,
        Ayanamsha::TrueRevati,
        Ayanamsha::TruePushya,
        Ayanamsha::TrueMula,
        Ayanamsha::TrueSheoran,
        Ayanamsha::GalcentCochrane,
        Ayanamsha::GalcentMulaWilhelm,
    ];

    fn definition(self) -> Definition {
        use Ayanamsha::*;
        match self {
            Tropical => Definition::None,
            Lahiri => Definition::Epoch(23.857_054),
            Krishnamurti => Definition::Epoch(23.760_470),
            Raman => Definition::Epoch(22.411_021),
            Yukteshwar => Definition::Epoch(22.478_803),
            JnBhasin => Definition::Epoch(22.762_137),
            Ushashashi => Definition::Epoch(20.057_541),
            Sassanian => Definition::Epoch(19.992_959),
            Hipparchos => Definition::Epoch(20.247_787),
            Aldebaran15Tau => Definition::Epoch(24.759_425),
            ValensMoon => Definition::Epoch(22.795_609),
            TrueCitra => Definition::Star(SPICA, 180.0),
            TrueRevati => Definition::Star(ZETA_PSC, dms(359.0, 50.0, 0.0)),
            TruePushya => Definition::Star(DELTA_CNC, 106.0),
            TrueMula => Definition::Star(LAMBDA_SCO, 240.0),
            TrueSheoran => Definition::Star(DELTA_CNC, 103.492_644),
            GalcentCochrane => Definition::Star(GAL_CENTER, 270.0),
            // Galactic Centre in the middle of nakshatra Mula (246°40′). Not verifiable against
            // the reference API, which returns 0 for this system.
            GalcentMulaWilhelm => Definition::Star(GAL_CENTER, dms(246.0, 40.0, 0.0)),
        }
    }

    /// Two-letter key used by the AstroCalc API (`None` for tropical).
    pub const fn key(self) -> Option<&'static str> {
        use Ayanamsha::*;
        Some(match self {
            Tropical => return None,
            Lahiri => "lh",
            Krishnamurti => "kr",
            Raman => "ra",
            Yukteshwar => "yu",
            JnBhasin => "jb",
            Ushashashi => "us",
            Sassanian => "sa",
            Hipparchos => "hi",
            Aldebaran15Tau => "at",
            ValensMoon => "vm",
            TrueCitra => "tc",
            TrueRevati => "tr",
            TruePushya => "tp",
            TrueMula => "tm",
            TrueSheoran => "ts",
            GalcentCochrane => "gc",
            GalcentMulaWilhelm => "gm",
        })
    }

    /// Snake-case name used by the AstroCalc API.
    pub const fn name(self) -> &'static str {
        use Ayanamsha::*;
        match self {
            Tropical => "tropical",
            Lahiri => "lahiri",
            Krishnamurti => "krishnamurti",
            Raman => "raman",
            Yukteshwar => "yukteshwar",
            JnBhasin => "jnbhasin",
            Ushashashi => "ushashashi",
            Sassanian => "sassanian",
            Hipparchos => "hipparchos",
            Aldebaran15Tau => "aldebaran_15_tau",
            ValensMoon => "valensmoon",
            TrueCitra => "true_citra",
            TrueRevati => "true_revati",
            TruePushya => "true_pushya",
            TrueMula => "true_mula",
            TrueSheoran => "true_sheoran",
            GalcentCochrane => "galcent_cochrane",
            GalcentMulaWilhelm => "galcent_mula_wilhelm",
        }
    }

    /// Swiss Ephemeris `SE_SIDM_*` number, for compatibility with existing API parameters only.
    pub const fn legacy_num(self) -> u8 {
        use Ayanamsha::*;
        match self {
            Tropical => 0,
            Lahiri => 1,
            Raman => 3,
            Ushashashi => 4,
            Krishnamurti => 5,
            Yukteshwar => 7,
            JnBhasin => 8,
            Aldebaran15Tau => 14,
            Hipparchos => 15,
            Sassanian => 16,
            TrueCitra => 27,
            TrueRevati => 28,
            TruePushya => 29,
            TrueMula => 35,
            GalcentMulaWilhelm => 36,
            TrueSheoran => 39,
            GalcentCochrane => 40,
            ValensMoon => 42,
        }
    }

    pub fn from_legacy_num(n: u8) -> Option<Self> {
        Self::ALL.into_iter().find(|a| a.legacy_num() == n)
    }

    /// True for systems defined by a star's apparent position.
    pub fn is_star_based(self) -> bool {
        matches!(self.definition(), Definition::Star(..))
    }
}

impl fmt::Display for Ayanamsha {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.name())
    }
}

/// Parses a key (`"lh"`), name (`"true_citra"`, `"truecitra"`) or legacy number (`"27"`).
impl FromStr for Ayanamsha {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let k = s.trim().to_ascii_lowercase();
        if k.is_empty() || k == "0" || k == "tropical" {
            return Ok(Ayanamsha::Tropical);
        }
        if let Ok(n) = k.parse::<u8>() {
            return Self::from_legacy_num(n).ok_or_else(|| format!("unknown ayanamsha number: {s}"));
        }
        let compact = k.replace(['_', '-', ' '], "");
        Self::ALL
            .into_iter()
            .find(|a| a.key() == Some(k.as_str()) || a.name().replace('_', "") == compact)
            .ok_or_else(|| format!("unknown ayanamsha: {s}"))
    }
}

fn star_direction(star: &StarRef, jd_tt: f64) -> Vec3 {
    let years = (jd_tt - 2_451_545.0) / 365.25 + 2000.0 - star.epoch;
    let dec0 = star.dec_deg.to_radians();
    let dec = dec0 + (star.pm_dec / 3.6e6 * years).to_radians();
    let ra = star.ra_deg.to_radians() + (star.pm_ra / 3.6e6 * years / dec0.cos()).to_radians();
    [dec.cos() * ra.cos(), dec.cos() * ra.sin(), dec.sin()]
}

impl Ephemeris {
    /// Ayanamsha in degrees at a UT Julian Day. Subtract it from a tropical longitude to get
    /// the sidereal longitude. Star-based systems need the ephemeris for aberration.
    pub fn ayanamsha(&self, jd: impl Into<JdUt>, aya: Ayanamsha) -> Result<f64, CalcError> {
        let jd = jd.into().value();
        match aya.definition() {
            Definition::None => Ok(0.0),
            Definition::Epoch(at_j2000) => {
                let tt = crate::models::ut_to_tt(self.delta_t_model(), jd);
                let (dpsi, _) = nutation(tt);
                Ok(norm_deg(at_j2000 + (general_precession_longitude(tt) + dpsi).to_degrees()))
            }
            Definition::Star(star, sidereal_lng) => {
                let frame = Frame::new(self, jd, &GeoPos::default(), Perspective::Geocentric)?;
                let tt = crate::models::ut_to_tt(self.delta_t_model(), jd);
                let lng = frame.apparent_ecliptic_lng_of_star(star_direction(&star, tt));
                Ok(norm_deg(lng - sidereal_lng))
            }
        }
    }

    /// Several ayanamshas at once.
    pub fn ayanamshas(&self, jd: impl Into<JdUt>, ayas: &[Ayanamsha]) -> Result<Vec<f64>, CalcError> {
        let jd = jd.into();
        ayas.iter().map(|&a| self.ayanamsha(jd, a)).collect()
    }
}

/// Sidereal longitude from a tropical longitude and an ayanamsha value (degrees).
#[inline]
pub fn to_sidereal(tropical_lng: f64, ayanamsha: f64) -> f64 {
    norm_deg(tropical_lng - ayanamsha)
}
