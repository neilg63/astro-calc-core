//! Altitude (arcmin) of each body's centre at the reference API's rise/set times, per mode.
use astro_ephem_core::*;
fn main() {
    let eph = Ephemeris::open("data/de440s.bsp").unwrap().with_delta_t_model(DeltaTModel::AstroCalcApi);
    let cases: Vec<serde_json::Value> = serde_json::from_str(&std::fs::read_to_string("crates/astro_ephem_core/tests/fixtures/api_rise_set.json").unwrap()).unwrap();
    for key in ["su", "mo", "ju"] {
        let body: Body = key.parse().unwrap();
        for mode in 0..8 {
            let mut v = vec![];
            for c in &cases {
                let l = c["loc"].as_array().unwrap();
                let geo = GeoPos::new(l[0].as_f64().unwrap(), l[1].as_f64().unwrap(), 0.0).unwrap();
                let set = &c["modes"][mode.to_string()][key];
                for ev in ["rise", "set"] {
                    let t = set[ev].as_f64().unwrap();
                    if t < 1000.0 { continue; }
                    let h = eph.calc_position(body, t, geo, CoordSystem::Horizontal, Perspective::Topocentric).unwrap();
                    let sd = if key == "su" { (696_000.0f64 / (h.dist * spk::AU_KM)).asin().to_degrees() * 60.0 } else if key == "mo" { (1737.4f64 / (h.dist * spk::AU_KM)).asin().to_degrees() * 60.0 } else { 0.0 };
                    v.push((h.lat * 60.0, sd));
                }
            }
            let n = v.len() as f64;
            let mean = v.iter().map(|x| x.0).sum::<f64>() / n;
            let sd = v.iter().map(|x| x.1).sum::<f64>() / n;
            let spread = v.iter().map(|x| (x.0 - mean).abs()).fold(0.0, f64::max);
            let rel: Vec<String> = v.iter().take(4).map(|x| format!("{:.3}/{:.3}", x.0, x.1)).collect();
            println!("{key} mode {mode}: mean alt {mean:8.3}' spread {spread:6.3}' mean sd {sd:6.3}'  e.g. {}", rel.join(" "));
        }
    }
}
