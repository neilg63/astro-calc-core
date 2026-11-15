//! The Moon from the lunar theory ELP/MPP02 (Chapront & Francou, Astron. Astrophys. 404,
//! 735, 2003), with the constants fitted to the JPL ephemeris DE405 and the corrections to
//! the secular terms that keep it close to DE406 over six millennia.
//!
//! The series give the Moon's geocentric longitude, latitude and distance on the mean
//! ecliptic of date; they are then turned to the fixed ecliptic and equinox of J2000 with
//! Laskar's precession quantities P and Q, as the authors' reference subroutine does. The
//! table of terms, with the constants already applied, is packed by
//! `tools/build_analytic_tables.py`.

use crate::math::Vec3;

static TABLE: &[u8] = include_bytes!("elpmpp02.bin");

/// Bytes before the first term: magic, the mean longitude polynomial, 15 counts, padding.
const HEADER: usize = 112;
const TERM: usize = 48;
/// Arcseconds per radian.
const RAD: f64 = 648_000.0 / std::f64::consts::PI;
/// Days in a Julian century, the theory's unit of time.
const CENTURY: f64 = 36_525.0;
/// Ratio of the semi-major axes of DE405 and the theory, applied to the distance.
const DISTANCE_SCALE: f64 = 384_747.961_370_172_5 / 384_747.980_674_318;

/// Laskar's (1986) precession of the ecliptic: polynomial coefficients of P and Q.
const P: [f64; 5] = [0.101_803_91e-4, 0.470_204_39e-6, -0.541_736_7e-9, -0.250_794_8e-11, 0.463_486e-14];
const Q: [f64; 5] = [-0.113_469_002e-3, 0.123_726_74e-6, 0.126_541_7e-8, -0.137_180_8e-11, -0.320_334e-14];

#[inline]
fn f64_at(at: usize) -> f64 {
    f64::from_le_bytes(TABLE[at..at + 8].try_into().expect("8 bytes"))
}

#[inline]
fn count(i: usize) -> usize {
    let at = 48 + 4 * i;
    u32::from_le_bytes(TABLE[at..at + 4].try_into().expect("4 bytes")) as usize
}

/// Sum of one series and its derivative with respect to time in centuries: terms
/// `amplitude · t^power · sin(argument)`, the argument a quartic in `t`.
fn series(at: &mut usize, n: usize, t: &[f64; 5], power: usize) -> (f64, f64) {
    let (mut sum, mut rate) = (0.0, 0.0);
    // Smallest terms first, for accuracy: the series are roughly in decreasing order.
    for k in (0..n).rev() {
        let term = *at + k * TERM;
        let x = f64_at(term);
        let (f0, f1, f2, f3, f4) = (f64_at(term + 8), f64_at(term + 16), f64_at(term + 24), f64_at(term + 32), f64_at(term + 40));
        let y = f0 + f1 * t[1] + f2 * t[2] + f3 * t[3] + f4 * t[4];
        let yp = f1 + 2.0 * f2 * t[1] + 3.0 * f3 * t[2] + 4.0 * f4 * t[3];
        let (s, c) = y.sin_cos();
        sum += x * t[power] * s;
        rate += x * t[power] * yp * c;
        if power > 0 {
            rate += power as f64 * x * t[power - 1] * s;
        }
    }
    *at += n * TERM;
    (sum, rate)
}

/// Geocentric position (km) and velocity (km/day) of the Moon in the mean ecliptic and
/// equinox of J2000, at a TDB Julian Day.
pub fn geocentric_ecliptic_km(tdb_jd: f64) -> (Vec3, Vec3) {
    debug_assert_eq!(&TABLE[..8], b"ACELP001");
    let t1 = (tdb_jd - 2_451_545.0) / CENTURY;
    let t = [1.0, t1, t1 * t1, t1 * t1 * t1, t1 * t1 * t1 * t1];
    // Longitude, latitude (arcseconds) and distance (km), each with its rate per century.
    let mut v = [(0.0, 0.0); 3];
    let mut at = HEADER;
    for (iv, value) in v.iter_mut().enumerate() {
        *value = series(&mut at, count(iv), &t, 0);
    }
    for (iv, value) in v.iter_mut().enumerate() {
        for power in 0..4 {
            let (sum, rate) = series(&mut at, count(3 + iv * 4 + power), &t, power);
            value.0 += sum;
            value.1 += rate;
        }
    }
    let w1 = |k: usize| f64_at(8 + 8 * k);
    let lng = v[0].0 / RAD + w1(0) + w1(1) * t[1] + w1(2) * t[2] + w1(3) * t[3] + w1(4) * t[4];
    let lng_rate = v[0].1 / RAD + w1(1) + 2.0 * w1(2) * t[1] + 3.0 * w1(3) * t[2] + 4.0 * w1(4) * t[3];
    let (lat, lat_rate) = (v[1].0 / RAD, v[1].1 / RAD);
    let (dist, dist_rate) = (v[2].0 * DISTANCE_SCALE, v[2].1);

    // Rectangular coordinates on the mean ecliptic of date.
    let (sl, cl) = lng.sin_cos();
    let (sb, cb) = lat.sin_cos();
    let (cw, sw) = (dist * cb, dist * sb);
    let x = [cw * cl, cw * sl, sw];
    let radial = dist_rate * cb - lat_rate * sw;
    let xp = [radial * cl - lng_rate * x[1], radial * sl + lng_rate * x[0], dist_rate * sb + lat_rate * cw];

    // Turn to the fixed ecliptic of J2000.
    let pw = (P[0] + P[1] * t[1] + P[2] * t[2] + P[3] * t[3] + P[4] * t[4]) * t[1];
    let qw = (Q[0] + Q[1] * t[1] + Q[2] * t[2] + Q[3] * t[3] + Q[4] * t[4]) * t[1];
    let ra = 2.0 * (1.0 - pw * pw - qw * qw).sqrt();
    let (pwqw, pw2, qw2) = (2.0 * pw * qw, 1.0 - 2.0 * pw * pw, 1.0 - 2.0 * qw * qw);
    let (pwra, qwra) = (pw * ra, qw * ra);
    let pos = [
        pw2 * x[0] + pwqw * x[1] + pwra * x[2],
        pwqw * x[0] + qw2 * x[1] - qwra * x[2],
        -pwra * x[0] + qwra * x[1] + (pw2 + qw2 - 1.0) * x[2],
    ];
    let ppw = P[0] + (2.0 * P[1] + 3.0 * P[2] * t[1] + 4.0 * P[3] * t[2] + 5.0 * P[4] * t[3]) * t[1];
    let qpw = Q[0] + (2.0 * Q[1] + 3.0 * Q[2] * t[1] + 4.0 * Q[3] * t[2] + 5.0 * Q[4] * t[3]) * t[1];
    let (ppw2, qpw2) = (-4.0 * pw * ppw, -4.0 * qw * qpw);
    let ppwqpw = 2.0 * (ppw * qw + pw * qpw);
    let rap = (ppw2 + qpw2) / ra;
    let (ppwra, qpwra) = (ppw * ra + pw * rap, qpw * ra + qw * rap);
    let vel = [
        (pw2 * xp[0] + pwqw * xp[1] + pwra * xp[2] + ppw2 * x[0] + ppwqpw * x[1] + ppwra * x[2]) / CENTURY,
        (pwqw * xp[0] + qw2 * xp[1] - qwra * xp[2] + ppwqpw * x[0] + qpw2 * x[1] - qpwra * x[2]) / CENTURY,
        (-pwra * xp[0] + qwra * xp[1] + (pw2 + qw2 - 1.0) * xp[2] - ppwra * x[0] + qpwra * x[1] + (ppw2 + qpw2) * x[2]) / CENTURY,
    ];
    (pos, vel)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The reference positions in the authors' note (table 8, fit to DE405): km and km/day.
    const CHECKS: [(f64, Vec3, Vec3); 5] = [
        (2_500_000.5, [274_034.591_03, 252_067.536_89, -18_998.755_19], [-62_463.613_38, 65_693.963_92, 6_595.328_90]),
        (2_300_000.5, [353_104.313_59, -195_254.118_08, 34_943.545_92], [39_543.136_78, 74_373.180_70, -700.653_51]),
        (2_100_000.5, [-19_851.276_74, -385_646.177_17, -27_597.661_34], [87_539.407_44, -7_599.684_84, -4_960.443_60]),
        (1_900_000.5, [-370_342.792_54, -37_574.255_33, -4_527.918_40], [12_255.287_46, -89_710.975_08, 7_649.442_85]),
        (1_700_000.5, [-164_673.047_20, 367_791.713_29, 31_603.980_27], [-75_884.688_15, -35_802.265_58, -4_239.598_95]),
    ];

    #[test]
    fn matches_published_reference_positions() {
        for (jd, p, v) in CHECKS {
            let (pos, vel) = geocentric_ecliptic_km(jd);
            for i in 0..3 {
                assert!((pos[i] - p[i]).abs() < TOLERANCE_KM, "{jd} pos {i}: {} vs {}", pos[i], p[i]);
                assert!((vel[i] - v[i]).abs() < 10.0 * TOLERANCE_KM, "{jd} vel {i}: {} vs {}", vel[i], v[i]);
            }
        }
    }

    /// The embedded table is truncated (see the tool), which costs a few hundred metres; the
    /// full series reproduce these values to the last printed digit.
    const TOLERANCE_KM: f64 = 0.5;
}
