//! Tests star-position models against the reference API's star-based ayanamshas.
use astro_calc_core::*;
use astro_calc_core::math::*;
fn hms(h: f64, m: f64, s: f64) -> f64 { (h + m / 60.0 + s / 3600.0) * 15.0 }
fn dms(d: f64, m: f64, s: f64) -> f64 { d.signum() * (d.abs() + m / 60.0 + s / 3600.0) }
// (key, ra°, dec°, pmra mas/yr (μα*), pmdec mas/yr, epoch year)
fn stars() -> Vec<(&'static str, f64, f64, f64, f64, f64)> { vec![
    ("true_citra", hms(13.,25.,11.6045), dms(-11.,9.,40.482), -42.35, -30.67, 1991.25),
    ("true_revati", hms(1.,13.,43.8020), dms(7.,34.,31.762), 145.00, -55.69, 1991.25),
    ("true_pushya", hms(8.,44.,41.1101), dms(18.,9.,17.509), -17.67, -229.26, 1991.25),
    ("true_sheoran", hms(8.,44.,41.1101), dms(18.,9.,17.509), -17.67, -229.26, 1991.25),
    ("true_mula", hms(17.,33.,36.5264), dms(-37.,6.,13.495), -8.53, -30.80, 1991.25),
    ("galcent_cochrane", hms(17.,45.,40.0409), dms(-29.,0.,28.118), 0.0, 0.0, 2000.0),
]}
fn aberrate(dir: Vec3, eph: &Ephemeris, tdb: f64) -> Vec3 {
    let c = 299_792.458 * 86_400.0 / astro_calc_core::spk::AU_KM;
    let (e, ev) = eph.barycentric(399, tdb).unwrap();
    let (sun, _) = eph.barycentric(10, tdb).unwrap();
    let sun_to_obs = sub(e, sun);
    let em = norm(sun_to_obs);
    let eu = scale(sun_to_obs, 1.0 / em);
    // deflection for a source at infinity (q = p)
    let srs = 1.97412574336e-8;
    let qpe = add(dir, eu);
    let w = srs / em / dot(dir, qpe).max(1e-6);
    let dir = unit(add(dir, scale(cross(dir, cross(eu, dir)), w)));
    let v = scale(ev, 1.0 / c);
    let bm1 = (1.0 - dot(v, v)).sqrt();
    let pdv = dot(dir, v);
    let w1 = 1.0 + pdv / (1.0 + bm1);
    let w2 = srs / em;
    unit([0, 1, 2].map(|i| dir[i] * bm1 + w1 * v[i] + w2 * (v[i] - pdv * dir[i])))
}
fn main() {
    let eph = Ephemeris::open("data/de440s.bsp").unwrap();
    let rows: Vec<serde_json::Value> = serde_json::from_str(&std::fs::read_to_string("tests/fixtures/api_ayanamshas.json").unwrap()).unwrap();
    for (key, ra0, de0, pma, pmd, ep) in stars() {
        let mut res: Vec<[f64; 3]> = vec![];
        for r in &rows {
            let jd = r["jd"].as_f64().unwrap();
            let tt = models::ut_to_tt(DeltaTModel::AstroCalcApi, jd);
            let yrs = (tt - 2451545.0) / 365.25 + 2000.0 - ep;
            let de = de0 + pmd / 3.6e6 * yrs;
            let ra = ra0 + pma / 3.6e6 * yrs / de0.to_radians().cos();
            let (ra, de) = (ra.to_radians(), de.to_radians());
            let v0 = [de.cos() * ra.cos(), de.cos() * ra.sin(), de.sin()];
            let tdb = tt + models::tdb_minus_tt(tt) / 86400.0;
            let v = aberrate(v0, &eph, tdb);
            let (g, p, s, e) = models::precession_fw(tt);
            let (dpsi, deps) = models::nutation(tt);
            let mean = mat_mul(&rot_x(-e), &mat_mul(&rot_z(-s), &mat_mul(&rot_x(p), &rot_z(g))));
            let truem = mat_mul(&rot_x(-(e + deps)), &mat_mul(&rot_z(-(s + dpsi)), &mat_mul(&rot_x(p), &rot_z(g))));
            let ecl = |m: &Mat3, eps: f64| { let u = mat_vec(m, v); let w = mat_vec(&rot_x(eps), u); to_spherical(w).0.to_degrees() };
            let lm = ecl(&mean, e);
            let lt = ecl(&truem, e + deps);
            let api = r["values"][key].as_f64().unwrap();
            res.push([diff_deg(lm, api) * 3600.0, diff_deg(lt, api) * 3600.0, lt]);
        }
        let stat = |i: usize| { let v: Vec<f64> = res.iter().map(|x| x[i]).collect(); let mx = v.iter().cloned().fold(f64::MIN, f64::max); let mn = v.iter().cloned().fold(f64::MAX, f64::min); ((mx + mn) / 2.0, mx - mn) };
        let (om, rm) = stat(0); let (ot, rt) = stat(1);
        println!("{key:18} mean: offset {:12.3}\" range {:7.3}\" | true: offset {:12.3}\" ({:.6}°) range {:7.3}\"", om, rm, ot, ot / 3600.0, rt);
    }
}
