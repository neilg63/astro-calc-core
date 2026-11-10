//! Minimal 3-vector and 3×3 matrix helpers.
//!
//! Rotation matrices follow the SOFA/ERFA convention: `rot_x(a)` and `rot_z(a)` rotate the
//! coordinate axes (not the vector) by `a` radians, anticlockwise looking from +axis towards
//! the origin.

pub type Vec3 = [f64; 3];
pub type Mat3 = [[f64; 3]; 3];

pub const IDENTITY: Mat3 = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];

#[inline]
pub fn add(a: Vec3, b: Vec3) -> Vec3 {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

#[inline]
pub fn sub(a: Vec3, b: Vec3) -> Vec3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
pub fn scale(a: Vec3, s: f64) -> Vec3 {
    [a[0] * s, a[1] * s, a[2] * s]
}

#[inline]
pub fn dot(a: Vec3, b: Vec3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
pub fn cross(a: Vec3, b: Vec3) -> Vec3 {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

#[inline]
pub fn norm(a: Vec3) -> f64 {
    dot(a, a).sqrt()
}

#[inline]
pub fn unit(a: Vec3) -> Vec3 {
    scale(a, 1.0 / norm(a))
}

#[inline]
pub fn mat_vec(m: &Mat3, v: Vec3) -> Vec3 {
    [dot(m[0], v), dot(m[1], v), dot(m[2], v)]
}

/// `mᵀ · v`
#[inline]
pub fn mat_t_vec(m: &Mat3, v: Vec3) -> Vec3 {
    [
        m[0][0] * v[0] + m[1][0] * v[1] + m[2][0] * v[2],
        m[0][1] * v[0] + m[1][1] * v[1] + m[2][1] * v[2],
        m[0][2] * v[0] + m[1][2] * v[1] + m[2][2] * v[2],
    ]
}

pub fn mat_mul(a: &Mat3, b: &Mat3) -> Mat3 {
    let mut r = [[0.0; 3]; 3];
    for (i, row) in r.iter_mut().enumerate() {
        for (j, cell) in row.iter_mut().enumerate() {
            *cell = a[i][0] * b[0][j] + a[i][1] * b[1][j] + a[i][2] * b[2][j];
        }
    }
    r
}

pub fn rot_x(a: f64) -> Mat3 {
    let (s, c) = a.sin_cos();
    [[1.0, 0.0, 0.0], [0.0, c, s], [0.0, -s, c]]
}

pub fn rot_z(a: f64) -> Mat3 {
    let (s, c) = a.sin_cos();
    [[c, s, 0.0], [-s, c, 0.0], [0.0, 0.0, 1.0]]
}

/// Spherical angles (longitude in `[0, 2π)`, latitude) of a vector.
#[inline]
pub fn to_spherical(v: Vec3) -> (f64, f64) {
    let lng = v[1].atan2(v[0]).rem_euclid(std::f64::consts::TAU);
    let lat = v[2].atan2(v[0].hypot(v[1]));
    (lng, lat)
}

/// Normalises degrees to `[0, 360)`.
#[inline]
pub fn norm_deg(d: f64) -> f64 {
    let r = d.rem_euclid(360.0);
    if r >= 360.0 { 0.0 } else { r }
}

/// Signed shortest difference `b − a` in degrees, in `(-180, 180]`.
#[inline]
pub fn diff_deg(a: f64, b: f64) -> f64 {
    let d = (b - a).rem_euclid(360.0);
    if d > 180.0 { d - 360.0 } else { d }
}
