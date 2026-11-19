//! Angles, houses and ayanamshas against reference values harvested from the AstroCalc API.
//! Run with `--nocapture` for the report.

mod common;

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
    let Some(eph) = common::planetary() else {
        eprintln!("skipping: no planetary file in data/");
        return;
    };
    let eph = eph.with_delta_t_model(DeltaTModel::AstroCalcApi);
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

/// Gauquelin sectors against the reference. Towards the polar circles the reference stops
/// iterating early (its cusps miss the semi-arc relation that ours meet to 1e-12°), so the
/// tolerance widens with latitude. Beyond the circles it returns no usable cusps: there the
/// fallback is checked for its own properties in `gauquelin_sectors_have_their_properties`.
#[test]
fn gauquelin_sectors_match_reference() {
    let mut worst: BTreeMap<String, f64> = BTreeMap::new();
    let mut compared = 0;
    for c in fixture("api_gauquelin.json") {
        let jd = c["jd"].as_f64().unwrap();
        let geo = GeoPos::new(c["lat"].as_f64().unwrap(), c["lng"].as_f64().unwrap(), 0.0).unwrap();
        let g = gauquelin_sectors(jd, geo);
        if geo.lat.abs() > 66.0 {
            assert_eq!(g.fallback, Some(HouseSystem::Porphyry), "lat {}", geo.lat);
            continue;
        }
        assert_eq!(g.fallback, None);
        let tol = match geo.lat.abs() {
            l if l > 62.0 => 60.0,
            l if l > 58.0 => 1.0,
            _ => 0.02,
        };
        let cusps = c["cusps"].as_array().unwrap();
        assert_eq!(cusps.len(), 36);
        for (i, v) in cusps.iter().enumerate() {
            let d = math::diff_deg(v.as_f64().unwrap(), g.cusps[i]).abs() * 3600.0;
            assert!(d < tol, "sector cusp {} off by {d}\" at {} lat {}", i + 1, c["utc"], geo.lat);
            let e = worst.entry(format!("lat {:6.2}", geo.lat)).or_default();
            *e = e.max(d);
            compared += 1;
        }
    }
    assert!(compared >= 36 * 30, "only {compared} cusps compared");
    for (k, v) in &worst {
        eprintln!("G {k:12} max {v:9.4}\"");
    }
}

#[test]
fn gauquelin_sectors_have_their_properties() {
    let jd = 2451545.0;
    for lat in [-72.5, -33.865, 0.0, 28.6, 51.5, 64.5, 70.0, 85.0] {
        let geo = GeoPos::new(lat, 77.2, 0.0).unwrap();
        let (a, g) = (angles(jd, geo), gauquelin_sectors(jd, geo));
        // The four angles begin sectors 1, 10, 19 and 28.
        for (i, v) in [(0, a.ascendant), (9, a.mc), (18, a.descendant()), (27, a.ic())] {
            assert!(math::diff_deg(g.cusps[i], v).abs() < 1e-9, "lat {lat} cusp {}", i + 1);
        }
        // Sectors run clockwise: every step goes back in longitude, and 36 of them go once round.
        let steps: Vec<f64> = (0..36).map(|i| (g.cusps[i] - g.cusps[(i + 1) % 36]).rem_euclid(360.0)).collect();
        assert!(steps.iter().all(|&s| s > 0.0 && s < 180.0), "lat {lat}: {steps:?}");
        assert!((steps.iter().sum::<f64>() - 360.0).abs() < 1e-9, "lat {lat}");
        // Opposite sectors begin at opposite degrees.
        for i in 0..18 {
            assert!(math::diff_deg(g.cusps[i] + 180.0, g.cusps[i + 18]).abs() < 1e-9, "lat {lat} cusp {}", i + 1);
        }
        // Every third cusp is a Placidus cusp, or a Porphyry one where Placidus falls back.
        let p = houses(jd, geo, HouseSystem::Placidus);
        assert_eq!(g.fallback, p.fallback, "lat {lat}");
        for (sector, house) in [(1, 1), (4, 12), (7, 11), (10, 10), (13, 9), (16, 8), (19, 7), (22, 6), (25, 5), (28, 4), (31, 3), (34, 2)] {
            assert!(math::diff_deg(g.cusps[sector - 1], p.cusps[house - 1]).abs() < 1e-9, "lat {lat} sector {sector}");
        }
        // A longitude just inside each sector, and the cusp itself, belong to that sector.
        for i in 0..36 {
            assert_eq!(g.sector_of(g.cusps[i]), i as u8 + 1, "lat {lat}");
            assert_eq!(g.sector_of(g.cusps[i] - steps[i] / 2.0), i as u8 + 1, "lat {lat}");
        }
    }
    // Where it falls back, each quadrant is in nine equal parts of longitude.
    let g = gauquelin_sectors(jd, GeoPos::new(70.0, 77.2, 0.0).unwrap());
    assert_eq!(g.fallback, Some(HouseSystem::Porphyry));
    for quadrant in 0..4 {
        let step = |i: usize| (g.cusps[quadrant * 9 + i] - g.cusps[(quadrant * 9 + i + 1) % 36]).rem_euclid(360.0);
        assert!((1..9).all(|i| (step(i) - step(0)).abs() < 1e-9));
    }
}

#[cfg(feature = "serde")]
#[test]
fn gauquelin_sectors_round_trip_through_json() {
    let g = gauquelin_sectors(2451545.0, GeoPos::new(51.5, -0.12, 0.0).unwrap());
    let json = serde_json::to_value(g).unwrap();
    assert_eq!(json["cusps"].as_array().unwrap().len(), 36);
    assert_eq!(serde_json::from_value::<GauquelinSectors>(json.clone()).unwrap(), g);
    let mut short = json;
    short["cusps"].as_array_mut().unwrap().pop();
    assert!(serde_json::from_value::<GauquelinSectors>(short).is_err());
}
