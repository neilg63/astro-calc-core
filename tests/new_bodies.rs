//! Krusinski houses, asteroid magnitudes, lunar apsides and hypothetical bodies.

use std::path::PathBuf;

use astro_calc_core::*;
use serde_json::Value;

fn root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

fn eph() -> Option<Ephemeris> {
    let bsp = root().join("data/de440s.bsp");
    bsp.exists().then(|| Ephemeris::open(bsp).unwrap())
}

fn ecliptic(eph: &Ephemeris, body: Body, jd: f64) -> BodyPos {
    eph.calc_position(body, jd, GeoPos::default(), CoordSystem::Ecliptic, Perspective::Geocentric).unwrap()
}

#[test]
fn krusinski_houses_hold_their_defining_properties() {
    assert_eq!(HouseSystem::from_code('U'), Some(HouseSystem::Krusinski));
    assert_eq!(HouseSystem::Krusinski.code(), 'U');
    assert!(HouseSystem::ALL.contains(&HouseSystem::Krusinski));
    for (lat, lng) in [(51.5, -0.1), (-33.9, 151.2), (0.0, 30.0), (64.1, -21.9), (28.6, 77.2)] {
        let geo = GeoPos::new(lat, lng, 0.0).unwrap();
        for i in 0..24 {
            let jd = 2_451_545.0 + i as f64 * 37.39;
            let a = angles(jd, geo);
            let h = houses(jd, geo, HouseSystem::Krusinski);
            let near = |x: f64, y: f64| math::diff_deg(x, y).abs() < 1e-9;
            // The zenith and nadir give the MC and IC; the circle starts at the ascendant.
            assert!(near(h.cusps[0], a.ascendant) && near(h.cusps[9], a.mc), "{lat},{lng} {jd}: {h:?} {a:?}");
            assert!(near(h.cusps[3], a.mc + 180.0) && near(h.cusps[6], a.ascendant + 180.0));
            for k in 0..6 {
                assert!(near(h.cusps[k] + 180.0, h.cusps[k + 6]));
            }
            // Cusps advance through the zodiac in house order.
            for k in 0..12 {
                let step = math::diff_deg(h.cusps[k], h.cusps[(k + 1) % 12]);
                assert!(step > 0.0 && step < 180.0, "{lat},{lng} {jd}: cusp {} to next is {step}", k + 1);
            }
        }
    }
}

#[test]
fn asteroid_magnitudes_follow_the_hg_system() {
    // At zero phase and unit distances the magnitude is H itself.
    assert!((hg_magnitude(3.34, 0.12, 0.0, 1.0, 1.0) - 3.34).abs() < 1e-12);
    assert!(hg_magnitude(3.34, 0.12, 20.0, 1.0, 1.0) > 3.34);
    assert_eq!(hg_parameters(Body::Ceres), Some((3.34, 0.12)));
    assert_eq!(hg_parameters(Body::Mars), None);
    if !root().join("data/ceres.bsp").exists() {
        return;
    }
    let eph = Ephemeris::open_dir(root().join("data")).unwrap().with_delta_t_model(DeltaTModel::AstroCalcApi);
    let cases: Vec<Value> = serde_json::from_str(&std::fs::read_to_string(root().join("tests/fixtures/api_hypothetical.json")).unwrap()).unwrap();
    let mut worst = 0f64;
    for c in &cases {
        let p = eph.phenomena(Body::Ceres, c["jd"].as_f64().unwrap()).unwrap();
        worst = worst.max((p.magnitude - c["pheno"][0]["apparentMagnitude"].as_f64().unwrap()).abs());
        assert!(p.diameter > 0.0);
    }
    assert!(worst < 0.002, "Ceres magnitude differs by {worst}");
    for b in [Body::Pallas, Body::Juno, Body::Vesta, Body::Chiron, Body::Pholus] {
        assert!(eph.phenomena(b, 2_451_545.0).unwrap().magnitude.is_finite());
    }
}

#[test]
fn hypothetical_bodies_match_reference() {
    let Some(eph) = eph() else { return };
    let eph = eph.with_delta_t_model(DeltaTModel::AstroCalcApi);
    let cases: Vec<Value> = serde_json::from_str(&std::fs::read_to_string(root().join("tests/fixtures/api_hypothetical.json")).unwrap()).unwrap();
    assert!(cases.len() >= 20);
    for c in &cases {
        let jd = c["jd"].as_f64().unwrap();
        for b in c["bodies"].as_array().unwrap() {
            let body: Body = b["key"].as_str().unwrap().parse().unwrap();
            let p = ecliptic(&eph, body, jd);
            let f = |k: &str| b[k].as_f64().unwrap();
            assert!(math::diff_deg(f("lng"), p.lng).abs() * 3600.0 < 0.01, "{body} lng at {jd}");
            assert!((p.lat - f("lat")).abs() * 3600.0 < 0.01, "{body} lat at {jd}");
            assert!((p.lng_speed - f("lngSpeed")).abs() < 5e-6, "{body} speed at {jd}");
        }
    }
    // Every built-in hypothetical body computes and stays near its orbit.
    for body in [Body::Cupido, Body::Hades, Body::Zeus, Body::Kronos, Body::Apollon, Body::Admetos, Body::Vulkanus, Body::Poseidon, Body::Isis] {
        let el = hypothetical::elements(body).unwrap();
        let p = ecliptic(&eph, body, 2_460_000.5);
        let (near, far) = (el.semi_major_axis * (1.0 - el.eccentricity) - 1.1, el.semi_major_axis * (1.0 + el.eccentricity) + 1.1);
        assert!(p.dist > near && p.dist < far, "{body} at {} AU", p.dist);
        assert!(p.lat.abs() < 1.2 && p.lng_speed.abs() < 0.03);
        let r = math::norm(el.heliocentric(2_460_000.5));
        assert!(r >= el.semi_major_axis * (1.0 - el.eccentricity) - 1e-9 && r <= el.semi_major_axis * (1.0 + el.eccentricity) + 1e-9);
    }
    assert!(hypothetical::elements(Body::Mars).is_none());
}

#[test]
fn lunar_apsides() {
    let Some(eph) = eph() else { return };
    // Osculating apogee at J2000.0 TT against JPL Horizons osculating elements of the Moon
    // (ecliptic of J2000): longitude 252.99808°, latitude 4.07466°, apogee 405 988.838 km.
    // Ours is on the true equinox of date, so add the nutation in longitude (−13.93″).
    let jd = 2_451_545.0 - models::delta_t(2_451_545.0) / 86_400.0;
    let p = ecliptic(&eph, Body::OscuApogee, jd);
    assert!((p.lng - (252.998_08 - 13.93 / 3600.0)).abs() < 2e-4, "osculating apogee lng {}", p.lng);
    assert!((p.lat - 4.074_66).abs() < 2e-4, "osculating apogee lat {}", p.lat);
    assert!((p.dist * spk::AU_KM - 405_988.838).abs() < 0.5, "osculating apogee {} km", p.dist * spk::AU_KM);

    // Mean apogee: 263.46° at J2000 on the ecliptic, always within the orbit's inclination.
    let m = ecliptic(&eph, Body::MeanApogee, jd);
    assert!((m.lng - 263.464).abs() < 0.01 && m.lat.abs() < 5.15, "mean apogee {m:?}");
    assert!((m.lng_speed - 0.1114).abs() < 0.002);

    // The interpolated apsides coincide with the Moon when it is actually at an apsis:
    // where its geometric distance from the Earth stops growing (apogee) or shrinking.
    let start = 2_460_000.5;
    let radial = |t: f64| {
        let (r, v) = eph.relative(301, 399, t + models::delta_t(t) / 86_400.0).unwrap();
        math::dot(r, v)
    };
    let (mut apogees, mut perigees) = (0, 0);
    for i in 0..60 {
        let (mut a, mut b) = (start + i as f64, start + i as f64 + 1.0);
        let rising = radial(a) > 0.0;
        if (radial(b) > 0.0) == rising {
            continue;
        }
        while b - a > 1e-7 {
            let m = 0.5 * (a + b);
            if (radial(m) > 0.0) == rising { a = m } else { b = m }
        }
        let body = if rising { Body::IntpApogee } else { Body::IntpPerigee };
        if rising { apogees += 1 } else { perigees += 1 }
        let (moon, apsis) = (ecliptic(&eph, Body::Moon, a), ecliptic(&eph, body, a));
        // The Moon's position is apparent, the apsis geometric: they differ by aberration
        // and light-time, about 0.006° and 40 km.
        assert!(math::diff_deg(moon.lng, apsis.lng).abs() < 0.01, "{body}: moon {} apsis {}", moon.lng, apsis.lng);
        assert!((moon.lat - apsis.lat).abs() < 0.002 && (moon.dist - apsis.dist).abs() < 4e-7, "{body}: {moon:?} {apsis:?}");
    }
    assert!(apogees >= 2 && perigees >= 2);
    // Between apsides they stay in the neighbourhood of the mean apogee and perigee.
    for i in 0..400 {
        let t = start + i as f64;
        let mean = ecliptic(&eph, Body::MeanApogee, t).lng;
        let (ia, ip) = (ecliptic(&eph, Body::IntpApogee, t), ecliptic(&eph, Body::IntpPerigee, t));
        assert!(math::diff_deg(mean, ia.lng).abs() < 12.0, "apogee {} vs mean {mean} at {t}", ia.lng);
        assert!(math::diff_deg(mean + 180.0, ip.lng).abs() < 35.0, "perigee {} vs mean {} at {t}", ip.lng, mean + 180.0);
        assert!(ia.dist > 0.00268 && ia.dist < 0.00273 && ip.dist > 0.00237 && ip.dist < 0.00248);
    }
}

#[test]
fn fixed_stars() {
    let stars = named_stars();
    assert_eq!(stars.len(), 410);
    assert!(stars.windows(2).all(|w| w[0].hip < w[1].hip));
    let spica = star_by_name("spica").unwrap();
    assert_eq!((spica.hip, spica.name), (65474, "Spica"));
    assert_eq!(star_by_hip(32349).unwrap().name, "Sirius");
    assert_eq!(star_by_name("barnards star").unwrap().name, "Barnard's Star");
    assert!(star_by_name("not a star").is_none() && star_by_hip(1).is_none());
    assert_eq!("hip:65474".parse::<Body>().unwrap(), spica.body());

    let Some(eph) = eph() else { return };
    let geo = GeoPos::new(45.0, 10.0, 0.0).unwrap();
    // The star-based ayanamshas, already checked against the reference, are the apparent
    // longitudes of these stars less a fixed offset.
    for (name, aya, offset) in [("Spica", Ayanamsha::TrueCitra, 180.0), ("Revati", Ayanamsha::TrueRevati, 359.0 + 50.0 / 60.0)] {
        let star = star_by_name(name).unwrap();
        for jd in [2_415_020.5, 2_451_545.0, 2_460_000.5, 2_488_070.0] {
            let p = ecliptic(&eph, star.body(), jd);
            let expected = eph.ayanamsha(jd, aya).unwrap();
            assert!(math::diff_deg(expected, p.lng - offset).abs() < 2e-5, "{name} at {jd}: {} vs {expected}", p.lng - offset);
            assert!(p.lng_speed.abs() < 1e-3 && p.lat_speed.abs() < 1e-3);
        }
    }
    // Sirius: apparent place near its J2000 catalogue position, 8.6 light-years away, and
    // drifting south by its proper motion of 1.2″ a year.
    let sirius = star_by_name("Sirius").unwrap();
    let equ = eph.calc_position(sirius.body(), 2_451_545.0, geo, CoordSystem::Equatorial, Perspective::Geocentric).unwrap();
    assert!((equ.lng - 101.287).abs() < 0.02 && (equ.lat + 16.716).abs() < 0.02, "Sirius {equ:?}");
    assert!((equ.dist / 63_241.077 - 8.6).abs() < 0.05, "Sirius at {} light-years", equ.dist / 63_241.077);
    let drift = ecliptic(&eph, sirius.body(), 2_451_545.0 + 36_525.0).lat - ecliptic(&eph, sirius.body(), 2_451_545.0).lat;
    assert!((drift + 0.034).abs() < 0.015, "Sirius moved {drift} degrees in latitude in a century");
    // Polaris stays within a degree of the pole; topocentric and geocentric places agree.
    let polaris = star_by_name("Polaris").unwrap().body();
    let geoc = eph.calc_position(polaris, 2_460_000.5, geo, CoordSystem::Equatorial, Perspective::Geocentric).unwrap();
    let topo = eph.calc_position(polaris, 2_460_000.5, geo, CoordSystem::Equatorial, Perspective::Topocentric).unwrap();
    assert!(geoc.lat > 89.2 && (geoc.lat - topo.lat).abs() < 1e-4);
    let hor = eph.calc_position(polaris, 2_460_000.5, geo, CoordSystem::Horizontal, Perspective::Topocentric).unwrap();
    assert!((hor.lat - 45.0).abs() < 1.0, "Polaris altitude {}", hor.lat);

    // A star supplied by the caller gives the same result as the built-in one, and stars
    // work with the rise and set methods.
    let custom = eph.calc_star_position(spica, 2_460_000.5, geo, CoordSystem::Ecliptic, Perspective::Geocentric).unwrap();
    assert_eq!(custom, ecliptic(&eph, spica.body(), 2_460_000.5));
    let day = eph.day_transits(spica.body(), 2_460_000.5, geo, RiseSetOptions::default()).unwrap();
    assert!(day.rise.is_some() && day.set.is_some() && day.mc.is_some());
    assert!(matches!(
        eph.calc_position(Body::catalog(Catalog::Hipparcos, 1), 2_460_000.5, geo, CoordSystem::Ecliptic, Perspective::Geocentric),
        Err(CalcError::UnsupportedBody(_))
    ));
}
