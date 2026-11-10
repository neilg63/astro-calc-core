//! Infers the reference API's ΔT from its geocentric Moon/Sun longitudes.
//! Prints: utc, jd, our ΔT, implied ΔT (Moon), implied ΔT (Sun).
use astro_ephem_core::*;
fn main() {
    let eph = Ephemeris::open("data/de440s.bsp").unwrap();
    let rows: Vec<serde_json::Value> = serde_json::from_str(&std::fs::read_to_string("crates/astro_ephem_core/tests/fixtures/api_moon_series.json").unwrap()).unwrap();
    let g = GeoPos::default();
    for r in rows {
        let jd = r["jd"].as_f64().unwrap();
        let ours = models::delta_t(jd);
        let mut out = vec![];
        for k in ["mo", "su"] {
            let p = eph.calc_position(k.parse().unwrap(), jd, g, CoordSystem::Ecliptic, Perspective::Geocentric).unwrap();
            let d = math::diff_deg(r[k].as_f64().unwrap(), p.lng);
            out.push(ours - d / p.lng_speed * 86400.0);
        }
        println!("{} {:.1} {:.3} {:.3} {:.3}", &r["utc"].as_str().unwrap()[..7], jd, ours, out[0], out[1]);
    }
}
