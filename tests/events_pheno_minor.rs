//! Stations, lunar phases, the Indian day, minor bodies and phenomena, against fixtures
//! harvested from the AstroCalc API and a few published almanac values.

use std::path::{Path, PathBuf};

use astro_calc_core::*;
use serde_json::Value;

fn root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

fn data() -> PathBuf {
    root().join("data")
}

/// Planetary ephemeris only, with the reference API's ΔT.
fn eph() -> Option<Ephemeris> {
    let bsp = data().join("de440s.bsp");
    bsp.exists().then(|| Ephemeris::open(bsp).unwrap().with_delta_t_model(DeltaTModel::AstroCalcApi))
}

fn fixture(name: &str) -> Vec<Value> {
    let p: &Path = &root().join("tests/fixtures").join(name);
    serde_json::from_str(&std::fs::read_to_string(p).unwrap()).unwrap()
}

#[test]
fn stations_match_reference_table() {
    let Some(eph) = eph() else { return };
    for set in fixture("api_stations.json") {
        let (from, to) = (set["from"].as_f64().unwrap(), set["to"].as_f64().unwrap());
        for item in set["items"].as_array().unwrap() {
            let key = item["key"].as_str().unwrap();
            let body: Body = key.parse().unwrap();
            let ours = eph.stations(body, from, to).unwrap();
            // The reference table is only as fine as its sampling: minutes for the inner
            // planets, hours for the slow ones whose speed crosses zero very gently.
            let tol_min = match key {
                "me" => 5.0,
                "ve" => 10.0,
                "ma" => 30.0,
                "ju" => 90.0,
                "sa" => 180.0,
                _ => 480.0,
            };
            let mut checked = 0;
            for a in item["stations"].as_array().unwrap() {
                let (jd, kind) = (a["jd"].as_f64().unwrap(), a["type"].as_str().unwrap());
                if jd < from + 1.0 || jd > to - 1.0 || kind.ends_with("peak") {
                    continue;
                }
                let nearest = ours
                    .iter()
                    .filter(|o| o.kind().key() == kind)
                    .map(|o| (o.jd() - jd).abs() * 1440.0)
                    .fold(f64::INFINITY, f64::min);
                assert!(nearest < tol_min, "{key} {kind} at {jd}: nearest ours {nearest:.1} min away");
                checked += 1;
            }
            assert!(checked >= 8, "{key}: only {checked} stations checked");
            // every station of ours is a genuine zero or extreme of speed, in cycle order
            for w in ours.windows(2) {
                assert!(w[0].jd() < w[1].jd());
                assert_ne!(w[0].kind(), w[1].kind());
            }
            for s in &ours {
                let p = eph.calc_position(body, s.jd(), GeoPos::default(), CoordSystem::Ecliptic, Perspective::Geocentric).unwrap();
                assert!((p.lng_speed - s.speed()).abs() < 1e-6, "{key} {s:?}: speed {}", p.lng_speed);
                match s {
                    Station::Peak(p) => assert!(p.speed() > 0.0),
                    Station::RetroPeak(p) => assert!(p.speed() < 0.0),
                    _ => {}
                }
            }
        }
    }
}

#[test]
fn mercury_station_matches_almanac() {
    let bsp = data().join("de440s.bsp");
    if !bsp.exists() {
        return;
    }
    // observed ΔT: Mercury stationed direct on 2024-01-02 at 03:07 UTC
    let eph = Ephemeris::open(bsp).unwrap();
    let st = eph.stations(Body::Mercury, JdUt::from_iso("2023-12-28").unwrap(), JdUt::from_iso("2024-01-05").unwrap()).unwrap();
    assert_eq!(st.len(), 1);
    assert_eq!(st[0].kind(), StationKind::RetroEnd);
    assert_eq!(&JdUt(st[0].jd()).to_iso().unwrap()[..16], "2024-01-02T03:07");
    // geocentric new moon of 2024-01-11 at 11:57 UTC
    let nm = eph.next_moon_phase(JdUt::from_iso("2024-01-10").unwrap()).unwrap();
    assert_eq!(nm.phase, MoonPhase::New);
    assert_eq!(&nm.jd.to_iso().unwrap()[..16], "2024-01-11T11:57");
    assert!(eph.stations(Body::Sun, 2460000.0, 2460400.0).unwrap().is_empty());
}

#[test]
fn moon_phases_match_reference() {
    let Some(eph) = eph() else { return };
    let mut worst = 0.0f64;
    for c in fixture("api_moon_phases.json") {
        let ph = c["phases"].as_array().unwrap();
        let (first, last) = (ph[0]["jd"].as_f64().unwrap(), ph[ph.len() - 1]["jd"].as_f64().unwrap());
        // the reference API's phases are topocentric (here for 0°N 0°E)
        let ours = eph.moon_phases_topocentric(first - 0.5, last + 0.5, GeoPos::default()).unwrap();
        assert_eq!(ours.len(), ph.len());
        for (i, (o, a)) in ours.iter().zip(ph).enumerate() {
            assert_eq!(o.phase.num() as i64, a["num"].as_i64().unwrap());
            assert_eq!(o.phase.angle(), a["angle"].as_f64().unwrap());
            assert_eq!(o.phase.is_waxing(), a["waxing"].as_bool().unwrap());
            // the first entry of each reference list can be hours off (its search start)
            if i > 0 {
                let s = (o.jd.value() - a["jd"].as_f64().unwrap()).abs() * 86_400.0;
                assert!(s < 45.0, "phase {} off by {s}s in list {}", a["num"], c["utc"]);
                worst = worst.max(s);
            }
        }
        // geocentric phases are the same events, within the lunar parallax (~2 h)
        let geo = eph.moon_phases(first - 0.5, last + 0.5).unwrap();
        for (g, t) in geo.iter().zip(&ours) {
            assert_eq!(g.phase, t.phase);
            assert!((g.jd.value() - t.jd.value()).abs() < 0.2);
        }
    }
    eprintln!("topocentric moon phases max {worst:.1}s");
}

#[test]
fn indian_day_units() {
    // Values from the original ITime unit test.
    let d = IndianDay::new(2450037.25, 2450036.552391257, 2450036.991859111, 2450037.5529648396);
    assert_eq!((d.muhurta(), d.ghati(), d.vighati()), (20, 41, 49));
    assert!(d.lipta() > 57.1 && d.lipta() < 57.2);
    assert!(d.progress() > 0.697 && d.progress() < 0.698);
    assert!(!d.is_day_time());
    assert_eq!(d.to_ghati_string(), "41:49:57.11");
    assert_eq!(d.to_string(), d.to_ghati_string());
    assert_eq!(d.to_unit_string(), "41 gh 49 vi 57.11 li");
    assert_eq!(d.to_long_string(), "muhurta 20, ghati 41, vighati 49, lipta 57.11 (night)");
    assert_eq!(d.half_day_muhurta(), 6);
    assert!((d.day_length() - (d.daylight_length() + d.night_length())).abs() < 1e-12);
    // 1995-11-15: sunrise fell on a Wednesday (Sunday = 0)
    assert_eq!(d.weekday_index(0), 3);
}

#[test]
fn indian_day_from_ephemeris() {
    let Some(eph) = eph() else { return };
    let opts = RiseSetOptions::default();
    let delhi = GeoPos::new(28.6139, 77.209, 216.0).unwrap();
    for (i, hours) in [0.0, 5.0, 11.0, 17.0, 23.0].into_iter().enumerate() {
        let jd = 2460325.5 + hours / 24.0; // 2024-01-16 UTC
        let d = eph.indian_day(jd, delhi, opts).unwrap();
        assert_eq!(d.kind, SolarDayKind::Normal);
        assert!(d.rise_jd <= jd && jd < d.next_rise_jd, "case {i}");
        assert!(d.rise_jd < d.set_jd && d.set_jd < d.next_rise_jd);
        assert!((d.day_length() - 1.0).abs() < 0.01);
        assert!(d.ghati() < 60 && d.muhurta() < 30);
        assert_eq!(d.is_day_time(), jd < d.set_jd);
        // boundaries are real sunrises/sunsets: the Sun's centre sits at the refracted horizon
        for t in [d.rise_jd, d.set_jd, d.next_rise_jd] {
            let alt = eph.altitude(Body::Sun, t, GeoPos { alt: 0.0, ..delhi }).unwrap();
            assert!((alt * 60.0 + 36.739).abs() < 0.2, "altitude {alt} at {t}");
        }
    }
    // Sunrise in Delhi on 2024-01-16 was at about 07:15 IST = 01:45 UTC.
    let d = eph.indian_day(2460325.5 + 0.25, delhi, opts).unwrap();
    assert_eq!(&JdUt(d.rise_jd).to_iso().unwrap()[..15], "2024-01-16T01:4");
    let (day, pos) = eph
        .positions_at_sunrise(&[Body::Sun, Body::Moon], 2460325.75, delhi, CoordSystem::Ecliptic, Perspective::Topocentric, opts)
        .unwrap();
    assert_eq!(day.rise_jd, d.rise_jd);
    assert_eq!(pos.len(), 2);

    // Midnight sun and polar night at Tromsø.
    let tromso = GeoPos::new(69.65, 18.96, 0.0).unwrap();
    let summer = eph.indian_day(2460482.3, tromso, opts).unwrap();
    assert_eq!(summer.kind, SolarDayKind::PolarDay);
    assert!(summer.is_day_time() && summer.rise_jd <= 2460482.3 && 2460482.3 < summer.next_rise_jd);
    let winter = eph.indian_day(2460665.3, tromso, opts).unwrap();
    assert_eq!(winter.kind, SolarDayKind::PolarNight);
    assert!(!winter.is_day_time() && (winter.day_length() - 1.0).abs() < 0.01);
}

#[test]
fn minor_bodies_match_reference() {
    let dir = data();
    if !dir.join("de440s.bsp").exists() || !dir.join("chiron.bsp").exists() {
        eprintln!("skipping: small-body SPK files not present (tools/fetch_horizons_spk.py)");
        return;
    }
    let eph = Ephemeris::open_dir(&dir).unwrap().with_delta_t_model(DeltaTModel::AstroCalcApi);
    let mut worst = 0.0f64;
    for c in fixture("api_minor_bodies.json") {
        let jd = c["jd"].as_f64().unwrap();
        for b in c["bodies"].as_array().unwrap() {
            let key = b["key"].as_str().unwrap();
            let p = eph.calc_position(key.parse().unwrap(), jd, GeoPos::default(), CoordSystem::Ecliptic, Perspective::Geocentric).unwrap();
            let d = (math::diff_deg(b["lng"].as_f64().unwrap(), p.lng).abs()).max((p.lat - b["lat"].as_f64().unwrap()).abs()) * 3600.0;
            // The reference uses older orbit solutions: sub-arcsecond near the present,
            // drifting to ~2″ at the ends of 1850–2150.
            assert!(d < 2.5, "{key} off by {d}\" at {}", c["utc"]);
            assert!((p.lng_speed - b["lngSpeed"].as_f64().unwrap()).abs() < 2e-5);
            worst = worst.max(d);
        }
    }
    eprintln!("Chiron/Ceres/Juno max {worst:.3}\"");
    for body in [Body::Pallas, Body::Vesta, Body::Pholus, Body::catalog(Catalog::MinorPlanetCenter, 2060)] {
        let p = eph.calc_position(body, 2451545.0, GeoPos::default(), CoordSystem::Ecliptic, Perspective::Geocentric).unwrap();
        assert!(p.dist > 1.0 && p.dist < 40.0);
    }
}

#[test]
fn phenomena_match_reference() {
    let Some(eph) = eph() else { return };
    // radii (km) behind the reference API's disc sizes, to compare diameters like for like
    let ref_radius = |b: Body| match b {
        Body::Sun => 696_000.0,
        Body::Moon => 1_738.15,
        Body::Mercury => 2_439.0,
        Body::Venus => 6_052.0,
        Body::Mars => 3_397.2,
        Body::Jupiter => 71_398.0,
        Body::Saturn => 60_000.0,
        Body::Uranus => 25_400.0,
        Body::Neptune => 24_300.0,
        _ => 2_500.0,
    };
    for c in fixture("api_pheno.json") {
        let jd = c["jd"].as_f64().unwrap();
        let bodies: Vec<Body> = c["result"].as_array().unwrap().iter().map(|r| r["key"].as_str().unwrap().parse().unwrap()).collect();
        let ours = eph.phenomena_set(&bodies, jd).unwrap();
        for (p, r) in ours.iter().zip(c["result"].as_array().unwrap()) {
            let f = |k: &str| r[k].as_f64().unwrap();
            let key = p.body;
            assert!((p.phase_angle - f("phaseAngle")).abs() * 3600.0 < 30.0, "{key} phase angle");
            assert!((p.illuminated - f("phaseIlluminated")).abs() < 1e-4, "{key} illuminated");
            assert!((p.elongation - f("elongationOfPlanet")).abs() * 3600.0 < 0.6, "{key} elongation");
            let mag_tol = if key == Body::Saturn { 0.04 } else { 0.001 };
            assert!((p.magnitude - f("apparentMagnitude")).abs() < mag_tol, "{key} magnitude {} vs {}", p.magnitude, f("apparentMagnitude"));
            let dist = eph.calc_position(key, jd, GeoPos::default(), CoordSystem::Ecliptic, Perspective::Geocentric).unwrap().dist;
            let ref_diam = 2.0 * (ref_radius(key) / (dist * spk::AU_KM)).asin().to_degrees();
            assert!((ref_diam / f("apparentDiameterOfDisc") - 1.0).abs() < 1e-6, "{key} disc");
            assert!(p.diameter > 0.0);
        }
    }
}

#[test]
fn station_sets_report_motion_and_surrounding_stations() {
    let bsp = data().join("de440s.bsp");
    if !bsp.exists() {
        return;
    }
    let eph = Ephemeris::open(bsp).unwrap();
    // Mercury was retrograde from 2023-12-13 until 2024-01-02.
    let jd = JdUt::from_iso("2023-12-25T12:00:00").unwrap().value();
    let set = eph.next_stations(Body::Mercury, jd, 8, 3).unwrap();
    assert!(set.is_retrograde() && set.speed < 0.0);
    assert_eq!(set.stations.len(), 11);
    assert_eq!(set.previous().count(), 3);
    assert_eq!(set.upcoming().count(), 8);
    let kinds: Vec<_> = set.stations.iter().map(|s| s.kind()).collect();
    use StationKind::*;
    assert_eq!(kinds[..5], [Peak, RetroStart, RetroPeak, RetroEnd, Peak]);
    for w in set.stations.windows(2) {
        assert!(w[0].jd() < w[1].jd());
    }
    assert_eq!(&JdUt(set.next(RetroEnd).unwrap().jd()).to_iso().unwrap()[..16], "2024-01-02T03:07");
    assert_eq!(&JdUt(set.stations[1].jd()).to_iso().unwrap()[..10], "2023-12-13");
    assert_eq!(set.next(RetroEnd).unwrap().speed(), 0.0);

    // The same stations whichever way they are asked for.
    let last = set.stations.last().unwrap().jd();
    let by_days = eph.planetary_stations(Body::Mercury, jd, last - jd + 0.01, 3).unwrap();
    assert_eq!(by_days.stations.len(), set.stations.len());
    for (a, b) in by_days.stations.iter().zip(&set.stations) {
        assert!(a.kind() == b.kind() && (a.jd() - b.jd()).abs() < 1e-4, "{a:?} vs {b:?}");
    }
    let ranged = eph.stations(Body::Mercury, jd, last + 0.01).unwrap();
    assert_eq!(ranged.len(), 8);

    // Direct motion, slow bodies and bodies that never turn back.
    let mars = eph.next_stations(Body::Mars, jd, 4, 0).unwrap();
    assert!(!mars.is_retrograde());
    assert_eq!(mars.stations.len(), 4);
    assert_eq!(mars.stations[0].kind(), Peak);
    assert_eq!(eph.next_stations(Body::Pluto, jd, 12, 4).unwrap().stations.len(), 16);
    let sun = eph.next_stations(Body::Sun, jd, 4, 2).unwrap();
    assert!(sun.stations.is_empty() && sun.speed > 0.9);
    // Close to the end of the data the search returns what there is; outside it, an error.
    let late = eph.next_stations(Body::Mercury, 2_506_300.0, 40, 1).unwrap();
    assert!(!late.stations.is_empty() && late.stations.len() < 41);
    assert!(eph.next_stations(Body::Mercury, 2_600_000.0, 4, 0).is_err());
    assert!(eph.planetary_stations(Body::Mercury, jd, -1.0, 0).is_err());
}
