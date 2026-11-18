//! Transposed positions and transits against the ephemeris itself.

mod common;

use astro_calc_core::math::diff_deg;
use astro_calc_core::*;

fn eph() -> Option<Ephemeris> {
    common::planetary()
}

const SOURCE: f64 = 2_439_302.100_694_444_5;
const TARGET: f64 = 2_461_363.0;

/// Transposed onto its own moment, a body is where the ephemeris puts it in every system.
#[test]
fn a_position_transposed_to_its_own_moment_is_unchanged() {
    let Some(eph) = eph() else { return };
    let geo = GeoPos::new(56.05, -3.67, 10.0).unwrap();
    let bodies = [Body::Sun, Body::Moon, Body::Venus, Body::Saturn];
    let sources = eph.transposed_sources(&bodies, SOURCE, None).unwrap();
    let real = eph.calc_positions_ext(&bodies, SOURCE, geo, Perspective::Geocentric).unwrap();
    for (source, real) in sources.iter().zip(&real) {
        let moved = transposed_position(source, SOURCE, geo);
        assert_eq!((moved.body.kind(), moved.body.source()), (BodyKind::Transposed, Some(real.body)));
        for (got, want) in [(moved.lng, real.lng), (moved.ras, real.ras), (moved.azi, real.azi)] {
            assert!(diff_deg(got, want).abs() < 1e-7, "{:?}: {moved:?} vs {real:?}", real.body);
        }
        for (got, want) in [(moved.lat, real.lat), (moved.dcl, real.dcl), (moved.alt, real.alt)] {
            assert!((got - want).abs() < 1e-7, "{:?}: {moved:?} vs {real:?}", real.body);
        }
        assert_eq!((moved.dist, moved.dist_speed, moved.lng_speed, moved.ras_speed), (real.dist, 0.0, 0.0, 0.0));
    }
}

#[test]
fn transits_of_a_set_of_bodies() {
    let Some(eph) = eph() else { return };
    let (source_geo, geo) = (GeoPos::new(56.05, -3.67, 10.0).unwrap(), GeoPos::new(45.35, -87.67, 100.0).unwrap());
    let bodies = [Body::Sun, Body::Moon, Body::Mercury, Body::Venus, Body::Mars, Body::Jupiter, Body::Saturn];
    let opts = RiseSetOptions::default().with_refraction(Refraction::Off);
    let list = eph.transposed_transits(&bodies, SOURCE, None, TARGET, geo, opts).unwrap();
    let sources = eph.transposed_sources(&bodies, SOURCE, None).unwrap();
    let start = local_day_start(TARGET, geo.lng);
    for ((t, source), body) in list.iter().zip(&sources).zip(bodies) {
        assert_eq!((t.body, t.lng, t.lat, t.start.value()), (body, source.lng, source.lat, start));
        // At each event the transposed point is on the horizon or at its highest or lowest.
        let alt = |jd: JdUt| transposed_position(source, jd, geo).alt;
        assert!(alt(t.rise.unwrap()).abs() < 1e-7 && alt(t.set.unwrap()).abs() < 1e-7);
        assert!((alt(t.mc) - t.max_alt).abs() < 1e-6 && (alt(t.ic) - t.min_alt).abs() < 1e-6);
        // Rising, the altitude is increasing.
        assert!(transposed_position(source, t.rise.unwrap(), geo).alt_speed > 0.0);
    }
    // A source place makes the source topocentric, which moves the Moon most; the
    // geocentric option sets the place aside again.
    let topo = eph.transposed_transits(&bodies, SOURCE, Some(source_geo), TARGET, geo, opts).unwrap();
    assert!((topo[1].lng - list[1].lng).abs() > 0.1 && (topo[6].lng - list[6].lng).abs() < 0.01);
    assert_eq!(eph.transposed_transits(&bodies, SOURCE, Some(source_geo), TARGET, geo, opts.with_geocentric(true)).unwrap(), list);
    // On the ecliptic, the point rises as the ascendant.
    let flat = eph.transposed_transits(&bodies, SOURCE, None, TARGET, geo, opts.with_ignore_latitude(true)).unwrap();
    assert!(flat[1].lat == 0.0 && diff_deg(ascendant(flat[1].rise.unwrap(), geo), flat[1].lng).abs() < 1e-6);
    // Refraction lifts the point into view earlier.
    let refracted = eph.transposed_transits(&bodies, SOURCE, None, TARGET, geo, RiseSetOptions::default()).unwrap();
    assert!(refracted[0].rise.unwrap().value() < list[0].rise.unwrap().value());
}

#[test]
fn positions_of_a_set_of_bodies() {
    let Some(eph) = eph() else { return };
    let (h_geo, geo) = (GeoPos::new(56.05, -3.67, 10.0).unwrap(), GeoPos::new(45.35, -87.67, 100.0).unwrap());
    let bodies = [Body::Sun, Body::Moon, Body::Mars];
    for ps in [Perspective::Geocentric, Perspective::Topocentric] {
        let sources = eph.calc_positions(&bodies, SOURCE, h_geo, CoordSystem::Ecliptic, ps).unwrap();
        let list = eph.calc_transposed_positions(&bodies, TARGET, SOURCE, geo, h_geo, ps).unwrap();
        for (got, source) in list.iter().zip(&sources) {
            assert_eq!(*got, transposed_position(source, TARGET, geo));
            assert_eq!((got.body.source(), got.body.to_string()), (Some(source.body), format!("tp:{}", source.body)));
            // The ecliptic place is the source's; the declination has barely moved in sixty years.
            assert_eq!((got.lng, got.lat), (source.lng, source.lat));
            let then = eph.calc_position(source.body, SOURCE, h_geo, CoordSystem::Equatorial, ps).unwrap();
            assert!((got.dcl - then.lat).abs() < 0.02 && diff_deg(got.ras, then.lng).abs() < 0.02);
        }
    }
    // Local Space takes the topocentric source.
    let local = eph.calc_transposed_positions(&bodies, TARGET, SOURCE, geo, h_geo, Perspective::LocalSpaceHorizon).unwrap();
    assert_eq!(local, eph.calc_transposed_positions(&bodies, TARGET, SOURCE, geo, h_geo, Perspective::Topocentric).unwrap());
    // Away from the Earth there is no horizon.
    assert!(matches!(
        eph.calc_transposed_positions(&[Body::Earth, Body::Mars], TARGET, SOURCE, geo, h_geo, Perspective::Heliocentric),
        Err(CalcError::NoHorizon(_))
    ));
}
