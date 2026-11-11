//! Compares calc_positions against reference values harvested from the AstroCalc API
//! (Swiss Ephemeris 2.08). Run with `--nocapture` to see the per-body report.
//! Needs `data/de440s.bsp`; skipped otherwise.

use std::collections::BTreeMap;
use std::path::PathBuf;

use astro_calc_core::{Body, CoordSystem, DeltaTModel, Ephemeris, GeoPos, Perspective};
use serde_json::Value;

const BODIES: [&str; 10] = ["su", "mo", "me", "ve", "ma", "ju", "sa", "ur", "ne", "pl"];

fn diff_deg(a: f64, b: f64) -> f64 {
    let d = (b - a).rem_euclid(360.0);
    if d > 180.0 { d - 360.0 } else { d }
}

#[derive(Default)]
struct Stat {
    lng: f64,
    lat: f64,
    lng_speed: f64,
    lat_speed: f64,
    ra: f64,
    dec: f64,
}

impl Stat {
    fn max(&mut self, o: &Stat) {
        self.lng = self.lng.max(o.lng);
        self.lat = self.lat.max(o.lat);
        self.lng_speed = self.lng_speed.max(o.lng_speed);
        self.lat_speed = self.lat_speed.max(o.lat_speed);
        self.ra = self.ra.max(o.ra);
        self.dec = self.dec.max(o.dec);
    }
}

#[test]
fn positions_match_reference_api() {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let bsp = root.join("../../data/de440s.bsp");
    if !bsp.exists() {
        eprintln!("skipping: {} not present", bsp.display());
        return;
    }
    let eph = Ephemeris::open(&bsp).unwrap().with_delta_t_model(DeltaTModel::AstroCalcApi);
    let cases: Vec<Value> =
        serde_json::from_str(&std::fs::read_to_string(root.join("tests/fixtures/api_chart_data.json")).unwrap())
            .unwrap();
    let bodies: Vec<Body> = BODIES.iter().map(|k| k.parse().unwrap()).collect();

    let mut worst: BTreeMap<(i64, &str), Stat> = BTreeMap::new();
    for case in &cases {
        let jd = case["jd"].as_f64().unwrap();
        let loc = case["loc"].as_array().unwrap();
        let geo = GeoPos::new(loc[0].as_f64().unwrap(), loc[1].as_f64().unwrap(), loc[2].as_f64().unwrap()).unwrap();
        let topo = case["topo"].as_i64().unwrap();
        let persp = if topo == 1 { Perspective::Topocentric } else { Perspective::Geocentric };
        let ecl = eph.calc_positions(&bodies, jd, geo, CoordSystem::Ecliptic, persp).unwrap();
        let eq = eph.calc_positions(&bodies, jd, geo, CoordSystem::Equatorial, persp).unwrap();
        for (i, key) in BODIES.iter().enumerate() {
            let r = case["bodies"].as_array().unwrap().iter().find(|b| b["key"] == *key).unwrap();
            let f = |k: &str| r[k].as_f64().unwrap();
            let s = Stat {
                lng: diff_deg(f("lng"), ecl[i].lng).abs() * 3600.0,
                lat: (ecl[i].lat - f("lat")).abs() * 3600.0,
                lng_speed: (ecl[i].lng_speed - f("lngSpeed")).abs(),
                lat_speed: (ecl[i].lat_speed - f("latSpeed")).abs(),
                ra: diff_deg(f("rectAscension"), eq[i].lng).abs() * 3600.0,
                dec: (eq[i].lat - f("declination")).abs() * 3600.0,
            };
            worst.entry((topo, key)).or_default().max(&s);
        }
    }

    eprintln!("max |ours − API| over {} cases (arcsec; speeds in °/day)", cases.len());
    eprintln!("{:>4} {:>3} {:>9} {:>9} {:>9} {:>9} {:>11} {:>11}", "topo", "key", "lng\"", "lat\"", "ra\"", "dec\"", "lngSpd", "latSpd");
    for ((topo, key), s) in &worst {
        eprintln!(
            "{:>4} {:>3} {:>9.4} {:>9.4} {:>9.4} {:>9.4} {:>11.2e} {:>11.2e}",
            topo, key, s.lng, s.lat, s.ra, s.dec, s.lng_speed, s.lat_speed
        );
    }
    for ((_, key), s) in &worst {
        // Uranus–Pluto differ by DE431 (reference) vs DE440 (ours), growing away from ~2000.
        let tol = match *key {
            "ur" | "ne" | "pl" => 0.5,
            "mo" => 0.2,
            _ => 0.05,
        };
        assert!(s.lng < tol && s.lat < tol && s.ra < tol && s.dec < tol, "{key} exceeds {tol}\"");
    }
}
