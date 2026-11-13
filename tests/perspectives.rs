//! Observers other than the Earth, and Local Space.

use std::path::PathBuf;

use astro_calc_core::*;

fn eph() -> Option<Ephemeris> {
    let bsp = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("data/de440s.bsp");
    bsp.exists().then(|| Ephemeris::open(bsp).unwrap())
}

#[test]
fn keys_and_names_round_trip() {
    assert_eq!(Perspective::ALL.len(), 11);
    for p in Perspective::ALL {
        assert_eq!(p.key().parse::<Perspective>(), Ok(p));
        assert_eq!(p.name().parse::<Perspective>(), Ok(p));
        assert_eq!(p.to_string().to_uppercase().parse::<Perspective>(), Ok(p));
    }
    assert_eq!("Local Space".parse::<Perspective>(), Ok(Perspective::LocalSpaceHorizon));
    assert_eq!("ls".parse::<Perspective>(), Ok(Perspective::LocalSpaceHorizon));
    assert!("martian".parse::<Perspective>().is_err());
    assert_eq!(Perspective::Areocentric.center(), Some(Body::Mars));
    assert_eq!(Perspective::Barycentric.center(), None);
    let on_earth: Vec<_> = Perspective::ALL.into_iter().filter(|p| p.needs_location()).collect();
    assert_eq!(on_earth, [Perspective::Topocentric, Perspective::LocalSpaceHorizon]);
}

#[test]
fn views_from_other_centres() {
    let Some(eph) = eph() else { return };
    let geo = GeoPos::default();
    let ecl = |body: Body, jd: f64, ps: Perspective| eph.calc_position(body, jd, geo, CoordSystem::Ecliptic, ps);
    let opposite = |a: f64, b: f64| math::diff_deg(a + 180.0, b).abs();
    use Perspective::*;
    for jd in [2_415_100.5, 2_451_545.0, 2_460_000.5, 2_488_000.25] {
        // From the Sun the Earth lies opposite where the Sun appears from the Earth, at the
        // same distance. The remaining difference is the aberration in the Earth's view (20″).
        let (sun, earth) = (ecl(Body::Sun, jd, Geocentric).unwrap(), ecl(Body::Earth, jd, Heliocentric).unwrap());
        assert!(opposite(sun.lng, earth.lng) < 0.01, "{jd}: sun {} earth {}", sun.lng, earth.lng);
        // Distances agree to the few hundred km either body moves along the line in the light-time.
        assert!((sun.lat + earth.lat).abs() < 1e-3 && (sun.dist - earth.dist).abs() < 5e-6);
        assert!(earth.lng_speed > 0.94 && earth.lng_speed < 1.03);

        // The same symmetry from the Moon and from each planet.
        let (moon, back) = (ecl(Body::Moon, jd, Geocentric).unwrap(), ecl(Body::Earth, jd, Selenocentric).unwrap());
        assert!(opposite(moon.lng, back.lng) < 0.02 && (moon.dist - back.dist).abs() < 1e-6, "{jd}: moon {} earth {}", moon.lng, back.lng);
        for (planet, view) in [(Body::Mercury, Hermeocentric), (Body::Venus, Cytherocentric), (Body::Mars, Areocentric), (Body::Jupiter, Jovicentric), (Body::Saturn, Cronocentric)] {
            let (from_sun, to_sun) = (ecl(planet, jd, Heliocentric).unwrap(), ecl(Body::Sun, jd, view).unwrap());
            assert!(opposite(from_sun.lng, to_sun.lng) < 0.02, "{planet} at {jd}: {} vs {}", from_sun.lng, to_sun.lng);
            assert!((from_sun.dist - to_sun.dist).abs() < 1e-4);
            // The Earth and the other planets can be seen from there; the planet itself cannot.
            assert!(ecl(Body::Earth, jd, view).is_ok() && ecl(Body::Pluto, jd, view).is_ok());
            assert!(matches!(ecl(planet, jd, view), Err(CalcError::NotVisibleFrom { body, perspective }) if body == planet && perspective == view));
        }
        // An inner planet seen from the Earth is never far from the Sun; seen from Jupiter
        // the Earth never strays more than 12° from it.
        let (sun_j, earth_j) = (ecl(Body::Sun, jd, Jovicentric).unwrap(), ecl(Body::Earth, jd, Jovicentric).unwrap());
        assert!(math::diff_deg(sun_j.lng, earth_j.lng).abs() < 12.5);

        // The barycentre is inside or just outside the Sun: the Sun is visible from it,
        // close by, and the planets are almost where the Sun sees them.
        let sun_b = ecl(Body::Sun, jd, Barycentric).unwrap();
        assert!(sun_b.dist < 0.012, "Sun {} AU from the barycentre", sun_b.dist);
        let (mars_b, mars_h) = (ecl(Body::Mars, jd, Barycentric).unwrap(), ecl(Body::Mars, jd, Heliocentric).unwrap());
        assert!(math::diff_deg(mars_b.lng, mars_h.lng).abs() < 0.6 && math::diff_deg(mars_b.lng, mars_h.lng).abs() > 1e-4);
    }

    // What cannot be seen: the observer's own body, and the Earth's lunar points from elsewhere.
    let jd = 2_460_000.5;
    for (body, ps) in [(Body::Sun, Heliocentric), (Body::Moon, Selenocentric), (Body::Earth, Geocentric), (Body::Earth, Topocentric), (Body::TrueNode, Heliocentric), (Body::MeanApogee, Areocentric), (Body::MeanNode, Barycentric)] {
        assert!(matches!(ecl(body, jd, ps), Err(CalcError::NotVisibleFrom { .. })), "{body} from {ps}");
    }
    assert!(matches!(eph.calc_position(Body::Mars, jd, geo, CoordSystem::Horizontal, Heliocentric), Err(CalcError::NoHorizon(Heliocentric))));
    // Stars and hypothetical bodies are visible from anywhere.
    let spica = star_by_name("Spica").unwrap().body();
    let (near, far) = (ecl(spica, jd, Geocentric).unwrap(), ecl(spica, jd, Cronocentric).unwrap());
    assert!(math::diff_deg(near.lng, far.lng).abs() < 0.02);
    assert!(ecl(Body::Kronos, jd, Heliocentric).is_ok());

    // Heliocentric positions carry light-time but no aberration or deflection: they are
    // the plain direction from the Sun in the ecliptic of date.
    let tt = jd + models::delta_t(jd) / 86_400.0;
    let eo = models::EarthOrientation::new(jd, tt);
    let mut tau = 0.0;
    let mut v = [0.0; 3];
    for _ in 0..3 {
        v = math::sub(eph.barycentric(4, tt - tau).unwrap().0, eph.barycentric(10, tt).unwrap().0);
        tau = math::norm(v) / 173.144_632_674_24;
    }
    let u = math::mat_vec(&eo.npb, v);
    let (s, c) = eo.eps_true.sin_cos();
    let expected = math::to_spherical([u[0], c * u[1] + s * u[2], -s * u[1] + c * u[2]]).0.to_degrees();
    let mars = ecl(Body::Mars, jd, Heliocentric).unwrap();
    assert!(math::diff_deg(expected, mars.lng).abs() < 1e-6, "{} vs {expected}", mars.lng);
}

#[test]
fn local_space() {
    let Some(eph) = eph() else { return };
    let geo = GeoPos::new(51.5, -0.1, 0.0).unwrap();
    let jd = JdUt::from_iso("2024-06-21T12:00:00").unwrap();
    use Perspective::*;
    for body in [Body::Sun, Body::Moon, Body::Jupiter] {
        let topo = eph.calc_position(body, jd, geo, CoordSystem::Horizontal, Topocentric).unwrap();
        // Whatever coordinate system is asked for, Local Space answers on the horizon,
        // with azimuth from north through east where the horizontal system counts from south.
        for cs in [CoordSystem::Ecliptic, CoordSystem::Equatorial, CoordSystem::Horizontal] {
            let ls = eph.calc_position(body, jd, geo, cs, LocalSpaceHorizon).unwrap();
            assert!(math::diff_deg(ls.lng, topo.lng + 180.0).abs() < 1e-9 && (ls.lat - topo.lat).abs() < 1e-12, "{body} {cs:?}");
            assert_eq!(ls.dist, topo.dist);
        }
    }
    // At midday in London in June the Sun is due south and high; it rose in the north-east.
    let noon = eph.calc_position(Body::Sun, jd, geo, CoordSystem::Ecliptic, LocalSpaceHorizon).unwrap();
    assert!((noon.lng - 180.0).abs() < 3.0 && (noon.lat - 61.9).abs() < 0.5, "noon Sun {noon:?}");
    let dawn = eph.calc_position(Body::Sun, JdUt::from_iso("2024-06-21T04:00:00").unwrap(), geo, CoordSystem::Ecliptic, LocalSpaceHorizon).unwrap();
    assert!(dawn.lng > 45.0 && dawn.lng < 60.0 && dawn.lat.abs() < 3.0, "dawn Sun {dawn:?}");
    // The Sun's azimuth increases through the day (it moves clockwise).
    assert!(noon.lng_speed > 0.0);
}
