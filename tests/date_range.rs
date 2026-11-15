//! The span of the loaded JPL files, what happens outside it, and the engine reported.
// Some tests here describe failures that the `analytic` feature turns into results.
#![cfg_attr(feature = "analytic", allow(unused_imports))]

//
// These run against whichever planetary file `data/` holds (any astro-calc-build tier), and check
// the exact dates of the files JPL publishes.

use astro_calc_core::{Body, BodyPosition, CalcEngine, CalcError, CoordSystem, DateRange, EphemError, Ephemeris, GeoPos, Perspective, RiseSetOptions, TransitEvent, scan_dir};

/// The span of each planetary file JPL publishes, as the library reports it.
const KNOWN: [(&str, &str); 3] = [
    ("de440s.bsp", "1849-12-26 to 2150-01-22"),
    ("de440.bsp", "1549-12-31 to 2650-01-25"),
    ("de442.bsp", "1549-12-31 to 2650-01-25"),
];

fn eph() -> Ephemeris {
    Ephemeris::open_dir(concat!(env!("CARGO_MANIFEST_DIR"), "/data")).expect("data files")
}

/// The installed planetary file's name.
fn planetary_name() -> String {
    let files = scan_dir(concat!(env!("CARGO_MANIFEST_DIR"), "/data")).expect("data files");
    files.into_iter().find(|f| f.is_planetary()).expect("a planetary file").name
}

/// A year before the planetary file starts.
fn before(range: DateRange) -> f64 {
    range.start - 365.0
}

#[cfg(not(feature = "analytic"))]
fn label(jd: f64) -> String {
    DateRange { start: jd, end: jd }.start_date()
}

fn geo() -> GeoPos {
    GeoPos::new(51.5, -0.1, 0.0).unwrap()
}

#[test]
fn date_range_matches_the_planetary_file() {
    let eph = eph();
    let range = eph.date_range().expect("planetary file loaded");
    assert_eq!(range.to_string(), format!("{} to {}", range.start_date(), range.end_date()));
    if let Some((_, span)) = KNOWN.iter().find(|(name, _)| *name == planetary_name()) {
        assert_eq!(range.to_string(), *span);
    }
    assert!(range.contains(2_451_545.0));
    assert!(!range.contains(before(range)));
    // Every core body shares it; nothing reaches further than the Earth's own data.
    let (a, b) = eph.coverage(399).unwrap();
    assert!(range.start >= a && range.end <= b);
}

#[test]
fn body_ranges() {
    let eph = eph();
    let all = eph.date_range().unwrap();
    assert_eq!(eph.body_range(Body::Sun), Some(all));
    assert_eq!(eph.body_range(Body::Moon), Some(all));
    assert_eq!(eph.body_range(Body::TrueNode), Some(all));
    assert_eq!(eph.body_range(Body::Ascendant), Some(all));
    // The asteroid files are shorter than the planetary one, and lie inside it.
    let chiron = eph.body_range(Body::Chiron).expect("chiron file loaded");
    assert!(chiron.start >= all.start && chiron.end <= all.end);
    assert!(chiron.days() < all.days());
}

#[test]
#[cfg(not(feature = "analytic"))]
fn positions_report_the_engine() {
    let eph = eph();
    let jd = 2_451_545.0;
    assert_eq!(eph.engine_at(jd), Some(CalcEngine::Jpl));
    let pos = eph.calc_position(Body::Mars, jd, geo(), CoordSystem::Ecliptic, Perspective::Geocentric).unwrap();
    assert_eq!(pos.engine, CalcEngine::Jpl);
    assert_eq!(pos.engine(), CalcEngine::Jpl);
    let ext = eph.calc_position_ext(Body::Moon, jd, geo(), Perspective::Topocentric).unwrap();
    assert_eq!(ext.engine, CalcEngine::Jpl);
    assert_eq!(ext.equatorial().engine, CalcEngine::Jpl);
    let lngs = eph.calc_longitudes(&[Body::Sun, Body::TrueNode, Body::SouthNode], jd, geo(), Perspective::Geocentric, None).unwrap();
    assert!(lngs.iter().all(|l| l.engine == CalcEngine::Jpl));
    assert_eq!(CalcEngine::Jpl.to_string(), "jpl");
}

#[test]
#[cfg(not(feature = "analytic"))]
fn out_of_range_says_what_is_covered() {
    let eph = eph();
    let range = eph.date_range().unwrap();
    let jd = before(range);
    assert_eq!(eph.engine_at(jd), None);
    let err = eph.calc_position(Body::Mars, jd, geo(), CoordSystem::Ecliptic, Perspective::Geocentric).unwrap_err();
    match &err {
        CalcError::Ephem(EphemError::NoData { covered: Some((a, b)), .. }) => {
            assert!(*a > jd && *b > *a);
        }
        other => panic!("unexpected error: {other:?}"),
    }
    let text = err.to_string();
    assert!(text.contains(&label(jd)), "{text}");
    assert!(text.contains(&range.to_string()), "{text}");
}

#[test]
fn unloaded_body_says_so() {
    let eph = eph();
    let err = eph.barycentric(20_099_999, 2_451_545.0).unwrap_err();
    assert_eq!(err, EphemError::NoData { target: 20_099_999, jd: 2_451_545.0, covered: None });
    assert!(err.to_string().contains("no JPL ephemeris data is loaded"));
}

#[test]
fn data_directory_listing() {
    let dir = concat!(env!("CARGO_MANIFEST_DIR"), "/data");
    let files = scan_dir(dir).expect("data files");
    // Same listing from a loaded ephemeris.
    assert_eq!(files, eph().files());
    let planetary: Vec<_> = files.iter().filter(|f| f.is_planetary()).collect();
    assert_eq!(planetary.len(), 1);
    assert_eq!(planetary[0].range, eph().date_range().unwrap());
    if let Some((_, span)) = KNOWN.iter().find(|(name, _)| *name == planetary[0].name) {
        assert_eq!(planetary[0].range.to_string(), *span);
    }
    assert!(planetary[0].targets.contains(&301) && planetary[0].size > 30_000_000);
    // The asteroid files hold one body each and are not planetary.
    let chiron = files.iter().find(|f| f.targets == [20_002_060]).expect("chiron file");
    assert!(!chiron.is_planetary());
    assert!(scan_dir(concat!(env!("CARGO_MANIFEST_DIR"), "/src")).is_err());
}

#[test]
fn every_kind_of_result_reports_the_engine() {
    let eph = eph();
    let jd = 2_451_545.0;
    let opts = RiseSetOptions::default();
    assert_eq!(eph.calc_longitude(Body::Sun, jd, geo(), Perspective::Geocentric, None).unwrap().engine, CalcEngine::Jpl);
    let alt = eph.altitude(Body::Sun, jd, geo()).unwrap();
    assert_eq!(alt.engine, CalcEngine::Jpl);
    assert!(alt.abs() <= 90.0 && *alt == alt.value); // dereferences to the number
    assert_eq!(eph.next_transit(Body::Sun, jd, geo(), TransitEvent::Rise, opts, 1.5).unwrap().engine, CalcEngine::Jpl);
    assert_eq!(eph.day_transits(Body::Moon, jd, geo(), opts).unwrap().engine, CalcEngine::Jpl);
    assert_eq!(eph.indian_day(jd, geo(), opts).unwrap().engine, CalcEngine::Jpl);
    assert_eq!(eph.phenomena(Body::Venus, jd).unwrap().engine, CalcEngine::Jpl);
    assert_eq!(eph.moon_elongation(jd).unwrap().engine, CalcEngine::Jpl);
    assert_eq!(eph.next_moon_phase(jd).unwrap().engine, CalcEngine::Jpl);
    assert_eq!(eph.stations(Body::Mercury, jd, jd + 200.0).unwrap().engine, CalcEngine::Jpl);
    assert_eq!(eph.next_stations(Body::Mars, jd, 2, 0).unwrap().engine, CalcEngine::Jpl);
}

#[test]
#[cfg(not(feature = "analytic"))]
fn errors_classify_themselves() {
    let eph = eph();
    let early = eph.altitude(Body::Sun, before(eph.date_range().unwrap()), geo()).unwrap_err();
    assert!(early.is_out_of_range());
    assert_eq!(early.engine(), Some(CalcEngine::Jpl));
    let unsupported = eph.calc_position(Body::Ascendant, 2_451_545.0, geo(), CoordSystem::Ecliptic, Perspective::Heliocentric).unwrap_err();
    assert!(!unsupported.is_out_of_range());
    assert_eq!(unsupported.engine(), None);
}
