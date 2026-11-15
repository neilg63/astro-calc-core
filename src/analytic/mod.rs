//! Built-in analytical series, used where a date lies outside the JPL files (feature
//! `analytic`). No data files are needed: the coefficients are compiled in.
//!
//! Planets and Sun: VSOP87, version E (barycentric rectangular coordinates, dynamical
//! ecliptic and equinox J2000). Bretagnon & Francou, Astron. Astrophys. 202, 309 (1988).
//! The theory states a precision of about 1″ over 4000 years either side of J2000 for
//! Mercury to Mars, 2000 years for Jupiter and Saturn and 6000 years for Uranus and Neptune.
//!
//! Moon: ELP/MPP02, in [`elp`]. Chapront & Francou, Astron. Astrophys. 404, 735 (2003).
//!
//! Both tables are packed from the published files by `tools/build_analytic_tables.py`.
//! There is no Pluto: it has no analytical source here.

pub mod elp;

use crate::math::Vec3;
use crate::spk::DateRange;

/// The dates the analytical engine accepts: 4000 Julian years either side of J2000 (TDB
/// Julian Days), about 2001 BC to AD 6000.
pub const RANGE: DateRange = DateRange { start: 2_451_545.0 - 1_461_000.0, end: 2_451_545.0 + 1_461_000.0 };

static VSOP87E: &[u8] = include_bytes!("vsop87e.bin");

const BODIES: usize = 9;
const COUNTS: usize = BODIES * 3 * 6;
const HEADER: usize = 8 + 4 * COUNTS;
/// Days in a thousand Julian years, the theory's unit of time.
const MILLENNIUM: f64 = 365_250.0;

/// Index of a body in the table, by NAIF id. The outer planets are given for their system
/// barycentres as well, which the theory does not distinguish from the planet.
fn index(target: i32) -> Option<usize> {
    Some(match target {
        10 => 0,
        1 | 199 => 1,
        2 | 299 => 2,
        399 => 3,
        4 | 499 => 4,
        5 | 599 => 5,
        6 | 699 => 6,
        7 | 799 => 7,
        8 | 899 => 8,
        _ => return None,
    })
}

/// Whether the engine has series for a body (NAIF id).
pub fn has(target: i32) -> bool {
    matches!(target, 3 | 301) || index(target).is_some()
}

/// Earth/Moon mass ratio (DE405, to which the lunar theory is fitted).
const EARTH_MOON_MASS_RATIO: f64 = 81.300_56;

/// Position (AU) and velocity (AU/day) of a body relative to the solar-system barycentre in
/// the ICRF, at a TDB Julian Day. `None` for a body the engine has no series for.
pub fn barycentric(target: i32, tdb_jd: f64) -> Option<(Vec3, Vec3)> {
    // The Moon and the Earth-Moon barycentre: the Earth plus (a share of) the lunar theory.
    let share = match target {
        301 => 1.0,
        3 => 1.0 / (1.0 + EARTH_MOON_MASS_RATIO),
        _ => return planet(target, tdb_jd),
    };
    let (earth_pos, earth_vel) = planet(399, tdb_jd)?;
    let (moon_pos, moon_vel) = elp::geocentric_ecliptic_km(tdb_jd);
    let scale = share / crate::spk::AU_KM;
    let (p, v) = (to_icrf(moon_pos), to_icrf(moon_vel));
    Some((
        [earth_pos[0] + scale * p[0], earth_pos[1] + scale * p[1], earth_pos[2] + scale * p[2]],
        [earth_vel[0] + scale * v[0], earth_vel[1] + scale * v[1], earth_vel[2] + scale * v[2]],
    ))
}

#[inline]
fn count(i: usize) -> usize {
    let at = 8 + 4 * i;
    u32::from_le_bytes(VSOP87E[at..at + 4].try_into().expect("4 bytes")) as usize
}

#[inline]
fn f64_at(at: usize) -> f64 {
    f64::from_le_bytes(VSOP87E[at..at + 8].try_into().expect("8 bytes"))
}

/// The Sun or a planet from VSOP87.
fn planet(target: i32, tdb_jd: f64) -> Option<(Vec3, Vec3)> {
    debug_assert_eq!(&VSOP87E[..8], b"ACVSOP01");
    let body = index(target)?;
    let t = (tdb_jd - 2_451_545.0) / MILLENNIUM;
    // Start of this body's terms: everything before it in the table.
    let first = body * 18;
    let mut at = HEADER + (0..first).map(count).sum::<usize>() * 24;
    let mut pos = [0.0; 3];
    let mut vel = [0.0; 3];
    for c in 0..3 {
        let mut power = 1.0; // t^alpha
        let mut lower = 0.0; // alpha * t^(alpha - 1)
        for alpha in 0..6 {
            let (mut sum, mut rate) = (0.0, 0.0);
            // Smallest terms first, for accuracy: the series are in decreasing order.
            let n = count(first + c * 6 + alpha);
            for k in (0..n).rev() {
                let term = at + k * 24;
                let (a, b, freq) = (f64_at(term), f64_at(term + 8), f64_at(term + 16));
                let (s, cos) = (b + freq * t).sin_cos();
                sum += a * cos;
                rate -= a * freq * s;
            }
            at += n * 24;
            pos[c] += power * sum;
            vel[c] += lower * sum + power * rate;
            lower = (alpha + 1) as f64 * power;
            power *= t;
        }
        vel[c] /= MILLENNIUM;
    }
    Some((to_icrf(pos), to_icrf(vel)))
}

/// Dynamical ecliptic and equinox J2000 to the J2000 equatorial frame (FK5, taken as the
/// ICRF): the rotation given with the theory.
fn to_icrf(v: Vec3) -> Vec3 {
    [
        v[0] + 0.000_000_440_360 * v[1] - 0.000_000_190_919 * v[2],
        -0.000_000_479_966 * v[0] + 0.917_482_137_087 * v[1] - 0.397_776_982_902 * v[2],
        0.397_776_982_902 * v[1] + 0.917_482_137_087 * v[2],
    ]
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Undoes [`to_icrf`], to compare with the theory's own check values.
    fn to_ecliptic(v: Vec3) -> Vec3 {
        // The transpose of the rotation.
        [
            v[0] - 0.000_000_479_966 * v[1],
            0.000_000_440_360 * v[0] + 0.917_482_137_087 * v[1] + 0.397_776_982_902 * v[2],
            -0.000_000_190_919 * v[0] - 0.397_776_982_902 * v[1] + 0.917_482_137_087 * v[2],
        ]
    }

    /// Values from `vsop87.chk`, distributed with the series.
    #[test]
    fn matches_published_check_values() {
        let cases: [(i32, f64, Vec3, Vec3); 5] = [
            (10, 2_451_545.0, [-0.0071415279, -0.0027881715, 0.0002061418], [0.0000053774, -0.0000074073, -0.0000000943]),
            (10, 2_268_920.0, [-0.0071641326, 0.0015976455, 0.0001823330], [-0.0000015913, -0.0000080258, 0.0000000583]),
            (199, 2_451_545.0, [-0.1372349394, -0.4500758422, -0.0243922379], [0.0213717757, -0.0064553869, -0.0024879611]),
            (399, 2_451_545.0, [-0.1842769826, 0.9644534530, 0.0002022407], [-0.0172022466, -0.0031661954, 0.0000000125]),
            (6, 2_451_545.0, [6.3992653459, 6.5672047374, -0.3688706482], [-0.0042869773, 0.0038829091, 0.0001028562]),
        ];
        for (id, jd, p, v) in cases {
            let (pos, vel) = barycentric(id, jd).unwrap();
            let (pos, vel) = (to_ecliptic(pos), to_ecliptic(vel));
            // The embedded table is truncated (see the tool), so agreement is to a few parts
            // in ten million of the body's distance, not to the last printed digit.
            let scale = 1.0 + crate::math::norm(p);
            for i in 0..3 {
                assert!((pos[i] - p[i]).abs() < 5e-7 * scale, "{id} {jd} pos {i}: {} vs {}", pos[i], p[i]);
                assert!((vel[i] - v[i]).abs() < 5e-8 * scale, "{id} {jd} vel {i}: {} vs {}", vel[i], v[i]);
            }
        }
    }
}
