//! Dates beyond DE440 (1550–2650). These run only when a longer planetary file is installed:
//! `astro-calc-build` tier 3 (1000 BCE–3000 CE) or tier 4 (all of DE441). With a DE440 file
//! they are skipped, and the same dates are covered, if at all, by the analytical series.

mod common;

use astro_calc_core::*;

/// 1000-01-01 and 2900-01-01, both outside DE440 and inside either long file.
const YEAR_1000: f64 = 2_086_308.0;
const YEAR_2900: f64 = 2_780_263.0;

fn geo() -> GeoPos {
    GeoPos::new(51.5, -0.1, 0.0).unwrap()
}

#[test]
fn the_long_file_is_used_beyond_de440() {
    let Some(eph) = common::planetary_covering(YEAR_1000, YEAR_2900) else { return };
    assert!(common::has_long_range());
    assert!(!common::DE440_RANGE.contains(YEAR_1000) && !common::DE440_RANGE.contains(YEAR_2900));
    let bodies = [Body::Sun, Body::Moon, Body::Mercury, Body::Venus, Body::Mars, Body::Jupiter, Body::Saturn, Body::Uranus, Body::Neptune, Body::Pluto];
    for jd in [YEAR_1000, YEAR_2900] {
        assert_eq!(eph.engine_at(jd), Some(CalcEngine::Jpl));
        let list = eph.calc_positions_ext(&bodies, jd, geo(), Perspective::Topocentric).unwrap();
        for p in &list {
            assert_eq!(p.engine, CalcEngine::Jpl);
            assert!([p.lng, p.lat, p.ras, p.dcl, p.azi, p.alt, p.lng_speed, p.dist].iter().all(|v| v.is_finite()), "{:?} at {jd}", p.body);
        }
        // The Sun stays on the ecliptic and moves about a degree a day.
        assert!(list[0].lat.abs() < 0.01 && (list[0].lng_speed - 0.9856).abs() < 0.05);
        // The Moon is within its orbit's tilt of the ecliptic.
        assert!(list[1].lat.abs() < 5.4 && (11.0..15.5).contains(&list[1].lng_speed));
    }
    // Rise and set work there too.
    let day = eph.day_transits(Body::Sun, YEAR_1000, geo(), RiseSetOptions::default()).unwrap();
    assert!(day.rise.is_some() && day.set.is_some());
}

/// Where both apply, the long file and the analytical series tell the same story.
#[cfg(feature = "analytic")]
#[test]
fn the_series_agree_with_the_long_file() {
    let Some(eph) = common::planetary_covering(YEAR_1000, YEAR_2900) else { return };
    let series = Ephemeris::analytic();
    for jd in [YEAR_1000, YEAR_2900] {
        for body in [Body::Sun, Body::Moon, Body::Mercury, Body::Venus, Body::Mars, Body::Jupiter, Body::Saturn] {
            let a = eph.calc_position(body, jd, geo(), CoordSystem::Ecliptic, Perspective::Geocentric).unwrap();
            let b = series.calc_position(body, jd, geo(), CoordSystem::Ecliptic, Perspective::Geocentric).unwrap();
            assert_eq!((a.engine, b.engine), (CalcEngine::Jpl, CalcEngine::Analytic));
            // A loose bound: the series drift with distance from the present.
            assert!(math::diff_deg(a.lng, b.lng).abs() < 0.05 && (a.lat - b.lat).abs() < 0.05, "{body} at {jd}: {} vs {}", a.lng, b.lng);
        }
    }
}
