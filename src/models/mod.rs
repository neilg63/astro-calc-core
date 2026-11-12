#![allow(clippy::excessive_precision)]
//! Time scales and Earth orientation models.
//!
//! References:
//! * IAU 2006 precession (Capitaine et al. 2003; Fukushima–Williams angles, Hilton et al. 2006)
//! * IAU 2000B nutation (McCarthy & Luzum 2003)
//! * Earth rotation angle and GMST (IERS Conventions 2010, ch. 5)
//! * Formulae and coefficients cross-checked against ERFA (BSD-3-Clause), see NOTICE.

#[allow(clippy::approx_constant)]
mod delta_t_table;
mod delta_t_api_table;
mod nut00b_table;

use std::f64::consts::TAU;

use crate::math::{Mat3, Vec3, mat_mul, rot_x, rot_z};

/// Arcseconds to radians.
pub const AS2R: f64 = 4.848_136_811_095_359_935_899_141e-6;
/// Julian Day of J2000.0.
pub const J2000: f64 = 2_451_545.0;
/// Days per Julian century.
pub const DJC: f64 = 36_525.0;

#[inline]
fn centuries(jd: f64) -> f64 {
    (jd - J2000) / DJC
}

// ---------------------------------------------------------------------------------------------
// Time scales
// ---------------------------------------------------------------------------------------------

/// Decimal year for a Julian Day (good enough for ΔT interpolation).
#[inline]
fn decimal_year(jd: f64) -> f64 {
    2000.0 + (jd - 2_451_544.5) / 365.2425
}

/// Long-term parabola of Morrison & Stephenson (2004): ΔT ≈ −20 + 32 u², u = (y − 1820)/100.
#[inline]
fn delta_t_parabola(year: f64) -> f64 {
    let u = (year - 1820.0) / 100.0;
    -20.0 + 32.0 * u * u
}

/// Half the long-term acceleration of ΔT implied by the parabola (s/yr²).
const DT_HALF_ACCEL: f64 = 32.0 / 10_000.0;

/// Which ΔT values to use.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum DeltaTModel {
    /// Observed USNO values (1657–present), with the Morrison–Stephenson long-term
    /// acceleration beyond the table. Best estimate of the real Earth rotation.
    #[default]
    Observed,
    /// Values implied by the AstroCalc reference API for 1850–2149,
    /// falling back to `Observed` outside that range. For regression tests against the API:
    /// they differ from observed values by up to ~3 s before 1950, and from 2020 on follow the
    /// 2019 forecast (about 3.4 s more than observed in 2026).
    AstroCalcApi,
}

fn interpolate(table: &[(f64, f64)], x: f64) -> f64 {
    let i = table.partition_point(|&(tx, _)| tx <= x).clamp(1, table.len() - 1);
    let (x0, v0) = table[i - 1];
    let (x1, v1) = table[i];
    v0 + (v1 - v0) * (x - x0) / (x1 - x0)
}

/// ΔT = TT − UT1 in seconds, observed model (see [`DeltaTModel::Observed`]).
pub fn delta_t(jd_ut: f64) -> f64 {
    use delta_t_table::DELTA_T_TABLE as T;
    let y = decimal_year(jd_ut);
    let (first_y, first_v) = T[0];
    let (last_y, last_v) = T[T.len() - 1];
    if y <= first_y {
        return delta_t_parabola(y) + (first_v - delta_t_parabola(first_y));
    }
    if y >= last_y {
        // Last year's mean rate plus the long-term acceleration.
        let (py, pv) = T[T.len() - 13];
        let rate = (last_v - pv) / (last_y - py);
        let dy = y - last_y;
        return last_v + rate * dy + DT_HALF_ACCEL * dy * dy;
    }
    interpolate(&T, y)
}

/// ΔT = TT − UT1 in seconds for the chosen model.
pub fn delta_t_with(model: DeltaTModel, jd_ut: f64) -> f64 {
    match model {
        DeltaTModel::Observed => delta_t(jd_ut),
        DeltaTModel::AstroCalcApi => {
            use delta_t_api_table::DELTA_T_API as A;
            if jd_ut < A[0].0 || jd_ut > A[A.len() - 1].0 {
                delta_t(jd_ut)
            } else {
                interpolate(&A, jd_ut)
            }
        }
    }
}

/// TT Julian Day for a UT Julian Day.
#[inline]
pub fn ut_to_tt(model: DeltaTModel, jd_ut: f64) -> f64 {
    jd_ut + delta_t_with(model, jd_ut) / 86_400.0
}

/// TDB − TT in seconds (main periodic terms, accurate to ~30 µs).
#[inline]
pub fn tdb_minus_tt(jd_tt: f64) -> f64 {
    let g = (357.53 + 0.985_600_28 * (jd_tt - J2000)).to_radians();
    0.001_657 * g.sin() + 0.000_014 * (2.0 * g).sin()
}

// ---------------------------------------------------------------------------------------------
// Precession, nutation, obliquity
// ---------------------------------------------------------------------------------------------

/// Mean obliquity of the ecliptic, IAU 2006 (radians).
pub fn obliquity_mean(jd_tt: f64) -> f64 {
    let t = centuries(jd_tt);
    (84_381.406
        + (-46.836_769 + (-0.000_183_1 + (0.002_003_40 + (-0.000_000_576 - 0.000_000_043_4 * t) * t) * t) * t)
            * t)
        * AS2R
}

/// Fukushima–Williams precession angles (γ̄, φ̄, ψ̄, εA), IAU 2006, radians.
pub fn precession_fw(jd_tt: f64) -> (f64, f64, f64, f64) {
    let t = centuries(jd_tt);
    let gamb = (-0.052_928
        + (10.556_378 + (0.493_204_4 + (-0.000_312_38 + (-0.000_002_788 + 0.000_000_026_0 * t) * t) * t) * t)
            * t)
        * AS2R;
    let phib = (84_381.412_819
        + (-46.811_016 + (0.051_126_8 + (0.000_532_89 + (-0.000_000_440 - 0.000_000_017_6 * t) * t) * t) * t)
            * t)
        * AS2R;
    let psib = (-0.041_775
        + (5_038.481_484 + (1.558_417_5 + (-0.000_185_22 + (-0.000_026_452 - 0.000_000_014_8 * t) * t) * t) * t)
            * t)
        * AS2R;
    (gamb, phib, psib, obliquity_mean(jd_tt))
}

/// General precession in longitude since J2000, IAU 2006 (radians).
pub fn general_precession_longitude(jd_tt: f64) -> f64 {
    let t = centuries(jd_tt);
    (5_028.796_195 + (1.105_434_8 + (0.000_079_64 + (-0.000_023_857 - 0.000_000_038_3 * t) * t) * t) * t)
        * t
        * AS2R
}

/// Mean longitude of the Moon's ascending node on the mean ecliptic of date
/// (Chapront et al. 1998, as given in Meeus, Astronomical Algorithms 2nd ed. ch. 47), degrees.
pub fn mean_lunar_node(jd_tt: f64) -> f64 {
    let t = centuries(jd_tt);
    125.044_547_9 - 1_934.136_289_1 * t + 0.002_075_4 * t * t + t * t * t / 467_441.0 - t.powi(4) / 60_616_000.0
}

/// Mean longitude of the Moon's perigee on the mean ecliptic of date, measured along the
/// ecliptic to the node and then along the orbit (Chapront et al., as given in Meeus,
/// Astronomical Algorithms 2nd ed. ch. 50), degrees.
pub fn mean_lunar_perigee(jd_tt: f64) -> f64 {
    let t = centuries(jd_tt);
    83.353_246_5 + 4_069.013_728_7 * t - 0.010_32 * t * t - t * t * t / 80_053.0 + t.powi(4) / 18_999_000.0
}

/// Mean anomaly of the Moon (Meeus ch. 47), degrees.
pub fn mean_lunar_anomaly(jd_tt: f64) -> f64 {
    let t = centuries(jd_tt);
    134.963_396_4 + 477_198.867_505_5 * t + 0.008_741_4 * t * t + t * t * t / 69_699.0 - t.powi(4) / 14_712_000.0
}

/// Nutation in longitude and obliquity, IAU 2000B (radians).
pub fn nutation(jd_tt: f64) -> (f64, f64) {
    use nut00b_table::NUT00B;
    const TURNAS: f64 = 1_296_000.0;
    const U2R: f64 = AS2R / 1e7;
    let t = centuries(jd_tt);
    let el = (485_868.249_036 + 1_717_915_923.217_8 * t).rem_euclid(TURNAS) * AS2R;
    let elp = (1_287_104.793_05 + 129_596_581.048_1 * t).rem_euclid(TURNAS) * AS2R;
    let f = (335_779.526_232 + 1_739_527_262.847_8 * t).rem_euclid(TURNAS) * AS2R;
    let d = (1_072_260.703_69 + 1_602_961_601.209_0 * t).rem_euclid(TURNAS) * AS2R;
    let om = (450_160.398_036 - 6_962_890.543_1 * t).rem_euclid(TURNAS) * AS2R;
    let (mut dp, mut de) = (0.0, 0.0);
    for (m, c) in NUT00B.iter().rev() {
        let arg = (m[0] as f64 * el + m[1] as f64 * elp + m[2] as f64 * f + m[3] as f64 * d + m[4] as f64 * om)
            .rem_euclid(TAU);
        let (s, co) = arg.sin_cos();
        dp += (c[0] + c[1] * t) * s + c[2] * co;
        de += (c[3] + c[4] * t) * co + c[5] * s;
    }
    // Fixed offsets in lieu of planetary nutation terms.
    const DPPLAN: f64 = -0.135e-3 * AS2R;
    const DEPLAN: f64 = 0.388e-3 * AS2R;
    (dp * U2R + DPPLAN, de * U2R + DEPLAN)
}

/// Earth orientation at one instant: everything that is shared by all bodies.
#[derive(Debug, Clone, Copy)]
pub struct EarthOrientation {
    /// Bias-precession-nutation matrix: GCRS → true equator and equinox of date.
    pub npb: Mat3,
    /// True obliquity of the ecliptic (radians).
    pub eps_true: f64,
    /// Mean obliquity (radians).
    pub eps_mean: f64,
    /// Nutation in longitude (radians).
    pub dpsi: f64,
    /// Greenwich apparent sidereal time (radians).
    pub gast: f64,
}

impl EarthOrientation {
    pub fn new(jd_ut: f64, jd_tt: f64) -> Self {
        let (gamb, phib, psib, epsa) = precession_fw(jd_tt);
        let (dpsi, deps) = nutation(jd_tt);
        // NPB = Rx(−ε) · Rz(−ψ) · Rx(φ̄) · Rz(γ̄)
        let npb = mat_mul(
            &rot_x(-(epsa + deps)),
            &mat_mul(&rot_z(-(psib + dpsi)), &mat_mul(&rot_x(phib), &rot_z(gamb))),
        );
        let gast = (gmst06(jd_ut, jd_tt) + equation_of_equinoxes(jd_tt, dpsi, epsa)).rem_euclid(TAU);
        Self {
            npb,
            eps_true: epsa + deps,
            eps_mean: epsa,
            dpsi,
            gast,
        }
    }
}

// ---------------------------------------------------------------------------------------------
// Earth rotation
// ---------------------------------------------------------------------------------------------

/// Earth rotation angle (radians), UT1 Julian Day.
pub fn era00(jd_ut1: f64) -> f64 {
    let t = jd_ut1 - J2000;
    let f = jd_ut1.rem_euclid(1.0);
    (TAU * (f + 0.779_057_273_264_0 + 0.002_737_811_911_354_48 * t)).rem_euclid(TAU)
}

/// Greenwich mean sidereal time, IAU 2006 (radians).
pub fn gmst06(jd_ut1: f64, jd_tt: f64) -> f64 {
    let t = centuries(jd_tt);
    (era00(jd_ut1)
        + (0.014_506
            + (4_612.156_534 + (1.391_581_7 + (-0.000_000_44 + (-0.000_029_956 - 0.000_000_036_8 * t) * t) * t) * t)
                * t)
            * AS2R)
        .rem_euclid(TAU)
}

/// Equation of the equinoxes: Δψ cos εA plus the two largest complementary terms.
pub fn equation_of_equinoxes(jd_tt: f64, dpsi: f64, epsa: f64) -> f64 {
    let t = centuries(jd_tt);
    let om = (450_160.398_036 - 6_962_890.543_1 * t).rem_euclid(1_296_000.0) * AS2R;
    dpsi * epsa.cos() + (2_640.96e-6 * om.sin() + 63.52e-6 * (2.0 * om).sin()) * AS2R
}

/// Earth's rotation rate in radians per UT day.
pub const EARTH_ROT_RATE: f64 = TAU * 1.002_737_811_911_354_48;

// ---------------------------------------------------------------------------------------------
// Observer
// ---------------------------------------------------------------------------------------------

/// WGS84 geodetic → Earth-fixed Cartesian (metres). Latitude/longitude in radians.
pub fn geodetic_to_itrs(lat: f64, lng: f64, height_m: f64) -> Vec3 {
    const A: f64 = 6_378_137.0;
    const F: f64 = 1.0 / 298.257_223_563;
    let (sp, cp) = lat.sin_cos();
    let w = (1.0 - F) * (1.0 - F);
    let ac = A / (cp * cp + w * sp * sp).sqrt();
    let r = (ac + height_m) * cp;
    [r * lng.cos(), r * lng.sin(), (w * ac + height_m) * sp]
}
