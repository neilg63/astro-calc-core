//! Black-box analysis of the reference API's ayanamsha series.
use astro_calc_core::*;
fn p_a(t: f64) -> f64 { // IAU 2006 general precession in longitude, arcsec
    (5028.796195 + (1.1054348 + (0.00007964 + (-0.000023857 - 0.0000000383 * t) * t) * t) * t) * t
}
fn main() {
    let rows: Vec<serde_json::Value> = serde_json::from_str(&std::fs::read_to_string("tests/fixtures/api_ayanamshas.json").unwrap()).unwrap();
    let keys: Vec<String> = rows[0]["values"].as_object().unwrap().keys().cloned().collect();
    for k in keys {
        let mut a = vec![]; let mut b = vec![];
        for r in &rows {
            let jd = r["jd"].as_f64().unwrap();
            let tt = models::ut_to_tt(DeltaTModel::AstroCalcApi, jd);
            let t = (tt - 2451545.0) / 36525.0;
            let v = r["values"][&k].as_f64().unwrap() * 3600.0;
            let (dpsi, _) = models::nutation(tt);
            let p = p_a(t);
            a.push(v - p);
            b.push(v - p - dpsi.to_degrees() * 3600.0);
        }
        let range = |x: &Vec<f64>| x.iter().cloned().fold(f64::MIN, f64::max) - x.iter().cloned().fold(f64::MAX, f64::min);
        let mean = |x: &Vec<f64>| x.iter().sum::<f64>() / x.len() as f64;
        println!("{k:22} mean-model range {:9.3}\"  +nut range {:9.3}\"  A(J2000)={:.6}°  A+nut={:.6}°", range(&a), range(&b), mean(&a) / 3600.0, mean(&b) / 3600.0);
    }
}
