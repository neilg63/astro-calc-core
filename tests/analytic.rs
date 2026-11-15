//! The built-in analytical series: agreement with the JPL files where both apply, and the
//! fallback outside the files.
#![cfg(feature = "analytic")]

use astro_calc_core::{Body, CalcEngine, CalcError, CoordSystem, EphemError, Ephemeris, GeoPos, Perspective, math};

fn jpl() -> Ephemeris {
    Ephemeris::open_dir(concat!(env!("CARGO_MANIFEST_DIR"), "/data")).expect("data files")
}

fn geo() -> GeoPos {
    GeoPos::new(51.5, -0.1, 0.0).unwrap()
}

/// Largest difference allowed from JPL in apparent geocentric direction, arcseconds.
const LIMITS: [(Body, f64); 10] = [
    (Body::Sun, 0.3),
    (Body::Moon, 0.5),
    (Body::TrueNode, 3.0),
    (Body::Mercury, 0.5),
    (Body::Venus, 0.5),
    (Body::Mars, 1.5),
    (Body::Jupiter, 1.5),
    (Body::Saturn, 1.5),
    (Body::Uranus, 3.0),
    (Body::Neptune, 4.0),
];

#[test]
fn series_agree_with_jpl_inside_the_files() {
    let (jpl, series) = (jpl(), Ephemeris::analytic());
    for i in 0..60 {
        let jd = 2_397_000.0 + 1_780.0 * i as f64; // 1850s to 2140s
        for (body, limit) in LIMITS {
            let a = jpl.calc_position(body, jd, geo(), CoordSystem::Ecliptic, Perspective::Geocentric).unwrap();
            let b = series.calc_position(body, jd, geo(), CoordSystem::Ecliptic, Perspective::Geocentric).unwrap();
            assert_eq!((a.engine, b.engine), (CalcEngine::Jpl, CalcEngine::Analytic));
            let (dl, db) = (math::diff_deg(a.lng, b.lng).abs() * 3600.0, (a.lat - b.lat).abs() * 3600.0);
            assert!(dl < limit && db < limit, "{body} at {jd}: Δlng {dl:.3}″ Δlat {db:.3}″");
            // The osculating node's speed swings widely and is sensitive to small differences.
            let speed_limit = if body == Body::TrueNode { 2e-3 } else { 1e-4 };
            assert!((a.lng_speed - b.lng_speed).abs() < speed_limit, "{body} speed: {} vs {}", a.lng_speed, b.lng_speed);
            assert!((a.dist - b.dist).abs() <= 2e-4 * a.dist, "{body} distance: {} vs {}", a.dist, b.dist);
        }
    }
}

#[test]
fn falls_back_outside_the_files() {
    let eph = jpl();
    let inside = 2_451_545.0;
    let outside = 2_300_000.5; // 1 February 1585
    assert_eq!(eph.engine_at(inside), Some(CalcEngine::Jpl));
    assert_eq!(eph.engine_at(outside), Some(CalcEngine::Analytic));
    let bodies = [Body::Ascendant, Body::Sun, Body::Moon, Body::TrueNode, Body::SouthNode, Body::MeanApogee, Body::Mercury, Body::Venus, Body::Mars, Body::Jupiter, Body::Saturn, Body::Uranus, Body::Neptune];
    // Inside the files the series are not used.
    assert!(eph.calc_positions(&bodies, inside, geo(), CoordSystem::Ecliptic, Perspective::Topocentric).unwrap().iter().all(|p| p.engine == CalcEngine::Jpl));
    // Outside them every result says so.
    let old = eph.calc_positions_ext(&bodies, outside, geo(), Perspective::Topocentric).unwrap();
    assert!(old.iter().all(|p| p.engine == CalcEngine::Analytic && p.lng.is_finite() && p.alt.is_finite()));
    let lngs = eph.calc_longitudes(&bodies, outside, geo(), Perspective::Geocentric, None).unwrap();
    assert!(lngs.iter().all(|l| l.engine == CalcEngine::Analytic));
    assert_eq!(eph.altitude(Body::Sun, outside, geo()).unwrap().engine, CalcEngine::Analytic);
    assert_eq!(eph.phenomena(Body::Venus, outside).unwrap().engine, CalcEngine::Analytic);
    // Heliocentric and planet-centred views work too.
    let helio = eph.calc_position(Body::Earth, outside, geo(), CoordSystem::Ecliptic, Perspective::Heliocentric).unwrap();
    let sun = eph.calc_position(Body::Sun, outside, geo(), CoordSystem::Ecliptic, Perspective::Geocentric).unwrap();
    assert!(math::diff_deg(helio.lng, sun.lng + 180.0).abs() < 0.02);
}

#[test]
fn the_join_at_the_edge_of_the_files_is_smooth() {
    let eph = jpl();
    let range = eph.date_range().unwrap();
    // (Not the Moon: its speed changes too much in a day for this simple check.)
    for body in [Body::Sun, Body::Mars, Body::Saturn] {
        let a = eph.calc_position(body, range.start - 0.6, geo(), CoordSystem::Ecliptic, Perspective::Geocentric).unwrap();
        let b = eph.calc_position(body, range.start + 0.6, geo(), CoordSystem::Ecliptic, Perspective::Geocentric).unwrap();
        assert_eq!((a.engine, b.engine), (CalcEngine::Analytic, CalcEngine::Jpl));
        // 1.2 days of motion at the mean of the two speeds, to within the series' error.
        let moved = math::diff_deg(a.lng, b.lng);
        let expected = 1.2 * 0.5 * (a.lng_speed + b.lng_speed);
        assert!((moved - expected).abs() * 3600.0 < 3.0, "{body}: moved {moved}, expected {expected}");
    }
}

#[test]
fn what_the_series_cannot_do_is_reported() {
    let eph = jpl();
    let outside = 2_300_000.5;
    // No Pluto or minor bodies in the series: the error names the JPL coverage.
    for body in [Body::Pluto, Body::Chiron] {
        let err = eph.calc_position(body, outside, geo(), CoordSystem::Ecliptic, Perspective::Geocentric).unwrap_err();
        assert!(matches!(err, CalcError::Ephem(EphemError::NoData { covered: Some(_), .. })), "{body}: {err}");
        assert!(err.is_out_of_range());
    }
    // Beyond the span of the series themselves.
    let err = eph.calc_position(Body::Mars, 500_000.0, geo(), CoordSystem::Ecliptic, Perspective::Geocentric).unwrap_err();
    assert!(matches!(err, CalcError::Ephem(EphemError::BeyondAnalytic { .. })), "{err}");
    assert_eq!(err.engine(), Some(CalcEngine::Analytic));
    assert!(err.is_out_of_range());
    assert!(err.to_string().contains("beyond the analytical series"), "{err}");
    assert_eq!(eph.engine_at(500_000.0), None);
}

#[test]
fn an_ephemeris_with_no_files() {
    let eph = Ephemeris::analytic();
    assert!(eph.files().is_empty() && eph.date_range().is_none());
    assert_eq!(eph.engine_at(2_451_545.0), Some(CalcEngine::Analytic));
    let p = eph.calc_position(Body::Sun, 2_451_545.0, geo(), CoordSystem::Ecliptic, Perspective::Geocentric).unwrap();
    assert!((p.lng - 280.37).abs() < 0.05, "{}", p.lng);
    let moon = eph.calc_position(Body::Moon, 2_451_545.0, geo(), CoordSystem::Ecliptic, Perspective::Geocentric).unwrap();
    assert!((moon.lng - 223.32).abs() < 0.05 && (moon.dist * 149_597_870.7 - 402_449.0).abs() < 100.0, "{} {}", moon.lng, moon.dist);
    // With no files there is nothing to say about coverage.
    let err = eph.calc_position(Body::Pluto, 2_451_545.0, geo(), CoordSystem::Ecliptic, Perspective::Geocentric).unwrap_err();
    assert!(matches!(err, CalcError::Ephem(EphemError::NoData { target: 9, covered: None, .. })), "{err}");
}
