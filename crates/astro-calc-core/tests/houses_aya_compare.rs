//! Angles, houses and ayanamshas against reference values harvested from the AstroCalc API.
//! Run with `--nocapture` for the report.

use std::collections::BTreeMap;
use std::path::PathBuf;

use astro_calc_core::*;
use serde_json::Value;

fn root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

fn fixture(name: &str) -> Vec<Value> {
    serde_json::from_str(&std::fs::read_to_string(root().join("tests/fixtures").join(name)).unwrap()).unwrap()
}

fn year(c: &Value) -> i32 {
    c["utc"].as_str().unwrap()[..4].parse().unwrap()
}

#[test]
fn angles_and_houses_match_reference() {
    let mut worst: BTreeMap<String, f64> = BTreeMap::new();
    for c in fixture("api_houses.json") {
        // After 2050 the reference's sidereal time carries a constant 1.78″ offset (an artefact of
        // its long-term model), so later dates are held to a looser tolerance.
        let late = year(&c) >= 2050;
        let jd = c["jd"].as_f64().unwrap();
        let geo = GeoPos::new(c["lat"].as_f64().unwrap(), c["lng"].as_f64().unwrap(), 0.0).unwrap();
        let a = angles(jd, geo);
        assert!((ascendant(jd, geo) - a.ascendant).abs() < 1e-12);
        let p = &c["points"];
        for (k, v) in [
            ("ascendant", a.ascendant),
            ("mc", a.mc),
            ("armc", a.armc),
            ("vertex", a.vertex),
            ("equasc", a.equasc),
            ("coasc1", a.coasc1),
            ("coasc2", a.coasc2),
            ("polasc", a.polasc),
        ] {
            let d = math::diff_deg(p[k].as_f64().unwrap(), v).abs() * 3600.0;
            let tol = if late { 6.0 } else { 0.02 };
            assert!(d < tol, "{k} off by {d}\" at {} lat {}", c["utc"], geo.lat);
            let e = worst.entry(format!("point {k}{}", if late { " (≥2050)" } else { "" })).or_default();
            *e = e.max(d);
        }
        for s in c["sets"].as_array().unwrap() {
            let code = s["system"].as_str().unwrap();
            let sys: HouseSystem = code.parse().unwrap();
            let h = houses(jd, geo, sys);
            for (i, v) in s["houses"].as_array().unwrap().iter().enumerate() {
                let d = math::diff_deg(v.as_f64().unwrap(), h.cusps[i]).abs() * 3600.0;
                // Placidus at high latitude: the reference stops iterating earlier.
                let tol = match (late, sys) {
                    (true, _) => 6.0,
                    // Horizon houses are near-degenerate on the equator.
                    (false, HouseSystem::Horizontal) if geo.lat == 0.0 => 5.0,
                    (false, HouseSystem::Placidus) => 2.5,
                    _ => 0.02,
                };
                assert!(d < tol, "{code} cusp {} off by {d}\" at {} lat {}", i + 1, c["utc"], geo.lat);
                let e = worst.entry(format!("houses {code}{}", if late { " (≥2050)" } else { "" })).or_default();
                *e = e.max(d);
            }
        }
    }
    for (k, v) in worst {
        eprintln!("{k:26} max {v:9.4}\"");
    }
}

#[test]
fn polar_fallback_is_reported() {
    let geo = GeoPos::new(70.0, 20.0, 0.0).unwrap();
    let h = houses(2451545.0, geo, HouseSystem::Placidus);
    assert_eq!(h.fallback, Some(HouseSystem::Porphyry));
    assert_eq!(h.cusps, houses(2451545.0, geo, HouseSystem::Porphyry).cusps);
    let w = houses(2451545.0, geo, HouseSystem::WholeSign);
    assert_eq!(w.fallback, None);
    assert_eq!(w.house_of(w.cusps[3] + 1.0), 4);
}

#[test]
fn ayanamshas_match_reference() {
    let bsp = root().join("../../data/de440s.bsp");
    if !bsp.exists() {
        eprintln!("skipping: {} not present", bsp.display());
        return;
    }
    let eph = Ephemeris::open(&bsp).unwrap().with_delta_t_model(DeltaTModel::AstroCalcApi);
    let mut worst: BTreeMap<&str, f64> = BTreeMap::new();
    for r in fixture("api_ayanamshas.json") {
        let jd = r["jd"].as_f64().unwrap();
        for aya in Ayanamsha::ALL {
            // The reference API returns 0 for Galcent Mula Wilhelm (unsupported there).
            if matches!(aya, Ayanamsha::Tropical | Ayanamsha::GalcentMulaWilhelm) {
                continue;
            }
            let api = r["values"][aya.name()].as_f64().unwrap();
            let d = math::diff_deg(api, eph.ayanamsha(jd, aya).unwrap()).abs() * 3600.0;
            let tol = if aya.is_star_based() { 1.1 } else { 0.15 };
            assert!(d < tol, "{aya} off by {d}\" at {}", r["utc"]);
            let e = worst.entry(aya.name()).or_default();
            *e = e.max(d);
        }
    }
    for (k, v) in worst {
        eprintln!("{k:22} max {v:7.4}\"");
    }
}

#[test]
fn ayanamsha_parsing() {
    assert_eq!("lh".parse::<Ayanamsha>().unwrap(), Ayanamsha::Lahiri);
    assert_eq!("27".parse::<Ayanamsha>().unwrap(), Ayanamsha::TrueCitra);
    assert_eq!("true_citra".parse::<Ayanamsha>().unwrap(), Ayanamsha::TrueCitra);
    assert_eq!("TrueCitra".parse::<Ayanamsha>().unwrap(), Ayanamsha::TrueCitra);
    assert_eq!("".parse::<Ayanamsha>().unwrap(), Ayanamsha::Tropical);
    assert!("zz".parse::<Ayanamsha>().is_err());
    for a in Ayanamsha::ALL {
        assert_eq!(Ayanamsha::from_legacy_num(a.legacy_num()), Some(a));
        assert_eq!(a.name().parse::<Ayanamsha>().unwrap(), a);
    }
    assert!((to_sidereal(10.0, 23.5) - 346.5).abs() < 1e-12);
}
