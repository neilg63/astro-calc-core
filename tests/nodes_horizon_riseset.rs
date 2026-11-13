//! Lunar nodes, horizontal coordinates and rise/set/transit times against fixtures harvested
//! from the AstroCalc API (/chart-data, /rise-set-times, /positions).

use std::path::PathBuf;

use astro_calc_core::*;
use serde_json::Value;

fn setup() -> Option<(Ephemeris, PathBuf)> {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let bsp = root.join("data/de440s.bsp");
    if !bsp.exists() {
        eprintln!("skipping: {} not present", bsp.display());
        return None;
    }
    Some((Ephemeris::open(&bsp).unwrap().with_delta_t_model(DeltaTModel::AstroCalcApi), root.join("tests/fixtures")))
}

fn load(dir: &std::path::Path, name: &str) -> Vec<Value> {
    serde_json::from_str(&std::fs::read_to_string(dir.join(name)).unwrap()).unwrap()
}

fn geo_of(c: &Value, with_alt: bool) -> GeoPos {
    let l = c["loc"].as_array().unwrap();
    let alt = if with_alt { l[2].as_f64().unwrap() } else { 0.0 };
    GeoPos::new(l[0].as_f64().unwrap(), l[1].as_f64().unwrap(), alt).unwrap()
}

fn late(c: &Value) -> bool {
    c["utc"].as_str().unwrap()[..4].parse::<i32>().unwrap() >= 2050
}

#[test]
fn nodes_and_horizontal_coordinates() {
    let Some((eph, dir)) = setup() else { return };
    let mut worst = [0.0f64; 4];
    for c in load(&dir, "api_nodes_azalt.json") {
        let jd = c["jd"].as_f64().unwrap();
        let geo = geo_of(&c, true);
        let persp = if c["topo"] == 1 { Perspective::Topocentric } else { Perspective::Geocentric };
        for b in c["bodies"].as_array().unwrap() {
            let key = b["key"].as_str().unwrap();
            let body: Body = key.parse().unwrap();
            let f = |k: &str| b[k].as_f64().unwrap();
            if matches!(body, Body::TrueNode | Body::MeanNode) {
                let p = eph.calc_position(body, jd, geo, CoordSystem::Ecliptic, persp).unwrap();
                let d = math::diff_deg(f("lng"), p.lng).abs() * 3600.0;
                let tol = if body == Body::TrueNode { 0.5 } else { 0.15 };
                assert!(d < tol, "{key} off by {d}\" at {}", c["utc"]);
                worst[usize::from(body == Body::MeanNode)] = worst[usize::from(body == Body::MeanNode)].max(d);
                assert!(p.lat.abs() < 1e-9);
            }
            // Horizontal: azimuth from south westwards, geometric altitude. After 2050 the
            // reference's sidereal time is offset by 1.78″.
            let h = eph.calc_position(body, jd, geo, CoordSystem::Horizontal, persp).unwrap();
            let daz = math::diff_deg(f("azimuth"), h.lng).abs() * 3600.0;
            let dalt = (f("altitude") - h.lat).abs() * 3600.0;
            // azimuth is ill-conditioned near the zenith; scale by cos(alt)
            let daz = daz * h.lat.to_radians().cos();
            let tol = if late(&c) { 4.0 } else if body == Body::TrueNode { 0.6 } else { 0.2 };
            assert!(daz < tol && dalt < tol, "{key} az/alt off by {daz}\"/{dalt}\" at {}", c["utc"]);
            worst[2] = worst[2].max(daz);
            worst[3] = worst[3].max(dalt);
        }
    }
    eprintln!("true node {:.3}\" mean node {:.3}\" az·cos(alt) {:.3}\" alt {:.3}\"", worst[0], worst[1], worst[2], worst[3]);
    let s = eph.calc_position(Body::SouthNode, 2451545.0, GeoPos::default(), CoordSystem::Ecliptic, Perspective::Geocentric).unwrap();
    let n = eph.calc_position(Body::TrueNode, 2451545.0, GeoPos::default(), CoordSystem::Ecliptic, Perspective::Geocentric).unwrap();
    assert!((math::diff_deg(n.lng, s.lng).abs() - 180.0).abs() < 1e-9);
}

fn secs(ours: Option<JdUt>, api: f64) -> f64 {
    match ours {
        Some(t) if api > 1000.0 => (t.value() - api).abs() * 86_400.0,
        None if api < 1000.0 => 0.0,
        _ => f64::INFINITY,
    }
}

#[test]
fn rise_set_and_transits_all_modes() {
    let Some((eph, dir)) = setup() else { return };
    let (mut worst_t, mut worst_alt) = (0.0f64, 0.0f64);
    for c in load(&dir, "api_rise_set.json") {
        let jd = c["jd"].as_f64().unwrap();
        // the reference ignores the observer's height for rise/set
        let geo = geo_of(&c, false);
        for (mode, sets) in c["modes"].as_object().unwrap() {
            let opts = RiseSetOptions::from_legacy_mode(mode.parse().unwrap());
            for (key, ev) in sets.as_object().unwrap() {
                let d = eph.day_transits(key.parse().unwrap(), jd, geo, opts).unwrap();
                for (name, ours) in [("rise", d.rise), ("set", d.set), ("mc", d.mc), ("ic", d.ic)] {
                    let s = secs(ours, ev[name].as_f64().unwrap());
                    assert!(s < 0.6, "{key} {name} mode {mode} off by {s}s at {}", c["utc"]);
                    worst_t = worst_t.max(s);
                }
                for (name, ours) in [("min", d.min_alt), ("max", d.max_alt)] {
                    let a = (ours.unwrap() - ev[name].as_f64().unwrap()).abs() * 3600.0;
                    assert!(a < 0.25, "{key} {name} mode {mode} off by {a}\"");
                    worst_alt = worst_alt.max(a);
                }
            }
        }
        for (key, obj) in [("su", "sunRiseSets"), ("mo", "moonRiseSets")] {
            let d = eph.day_transits(key.parse().unwrap(), jd, geo, RiseSetOptions::default()).unwrap();
            let p = &c["positions"][obj];
            for (name, ours) in [("prevSet", d.prev_set), ("rise", d.rise), ("set", d.set), ("mc", d.mc), ("ic", d.ic), ("nextRise", d.next_rise)] {
                let s = secs(ours, p[name].as_f64().unwrap());
                assert!(s < 0.6, "/positions {key} {name} off by {s}s at {}", c["utc"]);
            }
        }
    }
    eprintln!("rise/set/transit times max {worst_t:.3}s; min/max altitudes max {worst_alt:.3}\"");
}

#[test]
fn legacy_mode_mapping() {
    use DiscAlignment::*;
    let m = |n| {
        let o = RiseSetOptions::from_legacy_mode(n);
        (o.disc, o.refraction != Refraction::Off, o.fixed_disc)
    };
    assert_eq!(m(0), (UpperLimb, true, false));
    assert_eq!(m(1), m(3));
    assert_eq!(m(3), (Center, true, false));
    assert_eq!(m(4), (LowerLimb, false, false));
    assert_eq!(m(7), (UpperLimb, true, true));
    assert_eq!(RiseSetOptions::default(), RiseSetOptions::from_legacy_mode(3));
    assert!((local_day_start(2450094.9756944445, -3.6) - 2450094.51).abs() < 1e-9);
}

#[test]
fn circumpolar_sun_has_no_rise() {
    let Some((eph, _)) = setup() else { return };
    // Tromsø around the June solstice: midnight sun.
    let geo = GeoPos::new(69.65, 18.96, 0.0).unwrap();
    let d = eph.day_transits(Body::Sun, 2460482.0, geo, RiseSetOptions::default()).unwrap();
    assert_eq!((d.rise, d.set), (None, None));
    assert!(d.mc.is_some() && d.min_alt.unwrap() > 0.0);
}
