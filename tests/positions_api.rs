//! Behaviour of the public position API. Needs `data/de440s.bsp`; skipped otherwise.

use std::path::PathBuf;

use astro_calc_core::*;

fn eph() -> Option<Ephemeris> {
    let p = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("data/de440s.bsp");
    p.exists().then(|| Ephemeris::open(p).unwrap())
}

const ROME: GeoPos = GeoPos { lat: 41.9, lng: 12.5, alt: 20.0 };

#[test]
fn single_and_bulk_agree() {
    let Some(eph) = eph() else { return };
    let bodies = Body::CORE[..10].to_vec();
    let all = eph.calc_positions(&bodies, 2451545.0, ROME, CoordSystem::Ecliptic, Perspective::Topocentric).unwrap();
    for (b, p) in bodies.iter().zip(&all) {
        let one = eph.calc_position(*b, 2451545.0, ROME, CoordSystem::Ecliptic, Perspective::Topocentric).unwrap();
        assert_eq!(&one, p);
    }
}

#[test]
fn progress_sampling() {
    let Some(eph) = eph() else { return };
    let items = eph
        .calc_progress(Body::Mars, 2460310.5, 28.0, 0.5, ROME, CoordSystem::Ecliptic, Perspective::Geocentric)
        .unwrap();
    assert_eq!(items.len(), 56);
    assert_eq!(items[1].jd.value(), 2460311.0);
    let sets = eph
        .calc_progress_positions(&[Body::Sun, Body::Moon], 2460310.5, 1.0, 0.25, ROME, CoordSystem::Equatorial, Perspective::Topocentric)
        .unwrap();
    assert_eq!(sets.len(), 2);
    assert!(sets.iter().all(|s| s.items.len() == 4));
    assert!(matches!(
        eph.calc_progress(Body::Sun, 2460310.5, 1.0, 0.0, ROME, CoordSystem::Ecliptic, Perspective::Geocentric),
        Err(CalcError::InvalidInterval { .. })
    ));
}

#[test]
fn unsupported_inputs_are_errors() {
    let Some(eph) = eph() else { return };
    assert!(matches!(
        eph.calc_position(Body::Io, 2451545.0, ROME, CoordSystem::Ecliptic, Perspective::Geocentric),
        Err(CalcError::UnsupportedBody(Body::Io))
    ));
    // a supported body whose data file is not loaded
    assert!(matches!(
        eph.calc_position(Body::Chiron, 2451545.0, ROME, CoordSystem::Ecliptic, Perspective::Geocentric),
        Err(CalcError::Ephem(EphemError::NoData { target: 20_002_060, .. }))
    ));
    assert!(matches!(
        eph.calc_position(Body::Sun, 2200000.0, ROME, CoordSystem::Ecliptic, Perspective::Geocentric),
        Err(CalcError::Ephem(EphemError::NoData { .. }))
    ));
}

#[test]
fn observed_delta_t_is_sane() {
    use astro_calc_core::models::{delta_t, delta_t_with};
    assert!((delta_t(2451544.5) - 63.83).abs() < 0.05); // 2000-01-01
    assert!((delta_t(2461041.5) - 69.1).abs() < 0.2); // 2026-01-01
    let api = delta_t_with(DeltaTModel::AstroCalcApi, 2461041.5);
    assert!((api - 72.5).abs() < 0.1);
    assert!(delta_t(2506147.5) > 100.0 && delta_t(2506147.5) < 160.0); // 2149
}
