//! Longitude-only methods, the ascendant as a body, and the `default` and `basic` sets.

use std::path::PathBuf;
use std::time::Instant;

use astro_calc_core::*;

fn eph() -> Option<Ephemeris> {
    let bsp = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("data/de440s.bsp");
    bsp.exists().then(|| Ephemeris::open(bsp).unwrap())
}

fn keys(bodies: &[Body]) -> String {
    bodies.iter().map(|b| b.to_string()).collect::<Vec<_>>().join(",")
}

#[test]
fn body_sets() {
    assert_eq!(keys(&Body::DEFAULT_BODIES), "as,su,mo,me,ve,ma,ju,sa,ur,ne,pl,ke,ra");
    assert_eq!(keys(&Body::BASIC), "as,su,mo,me,ve,ma,ju,sa,ke,ra");
    // The Sun and Moon are luminaries, not planets.
    assert_eq!(keys(&Body::LUMINARIES), "su,mo");
    assert_eq!(keys(&Body::PLANETS), "me,ve,ma,ju,sa,ur,ne,pl");
    assert_eq!(keys(&Body::LUMINARIES_AND_PLANETS), "su,mo,me,ve,ma,ju,sa,ur,ne,pl");
    assert!(Body::LUMINARIES.iter().all(|b| b.kind() == BodyKind::Luminary) && Body::PLANETS.iter().all(|b| b.kind() == BodyKind::Planet));
    assert_eq!(parse_body_list("default").unwrap(), Body::DEFAULT_BODIES);
    // `core`, the set's former name, is still accepted.
    assert_eq!(parse_body_list("core").unwrap(), Body::DEFAULT_BODIES);
    assert_eq!(parse_body_list_for("Core", Perspective::Heliocentric).unwrap(), parse_body_list_for("default", Perspective::Heliocentric).unwrap());
    assert_eq!(parse_body_list("basic").unwrap(), Body::BASIC);
    // Sets can be mixed with keys; nothing is repeated.
    assert_eq!(keys(&parse_body_list("basic,ch,su,ur").unwrap()), "as,su,mo,me,ve,ma,ju,sa,ke,ra,ch,ur");
    assert_eq!(keys(&parse_body_list("ce, default").unwrap()), "ce,as,su,mo,me,ve,ma,ju,sa,ur,ne,pl,ke,ra");
    assert!(parse_body_list("default,zz").is_err());
    assert_eq!("as".parse::<Body>(), Ok(Body::Ascendant));
    assert_eq!((Body::Ascendant.key(), Body::Ascendant.kind()), (Some("as"), BodyKind::Angle));

    // Away from the Earth: the Earth comes in after Venus; the ascendant, the nodes and
    // the centre itself go out.
    use Perspective::*;
    assert_eq!(keys(&parse_body_list_for("default", Heliocentric).unwrap()), "mo,me,ve,ea,ma,ju,sa,ur,ne,pl");
    assert_eq!(keys(&parse_body_list_for("default", Areocentric).unwrap()), "su,mo,me,ve,ea,ju,sa,ur,ne,pl");
    assert_eq!(keys(&parse_body_list_for("basic", Selenocentric).unwrap()), "su,me,ve,ea,ma,ju,sa");
    assert_eq!(keys(&parse_body_list_for("basic", Barycentric).unwrap()), "su,mo,me,ve,ea,ma,ju,sa");
    for ps in [Geocentric, Topocentric, LocalSpaceHorizon] {
        assert_eq!(parse_body_list_for("default", ps).unwrap(), Body::DEFAULT_BODIES);
    }
}

#[test]
fn longitudes_ascendant_and_nodes() {
    let Some(eph) = eph() else { return };
    let geo = GeoPos::new(56.1, -3.6, 0.0).unwrap();
    let jd = 2_450_094.975_694_444_5;
    for ps in [Perspective::Geocentric, Perspective::Topocentric, Perspective::LocalSpaceHorizon] {
        let full = eph.calc_positions(&Body::DEFAULT_BODIES, jd, geo, CoordSystem::Ecliptic, if ps == Perspective::LocalSpaceHorizon { Perspective::Topocentric } else { ps }).unwrap();
        let lngs = eph.calc_longitudes(&Body::DEFAULT_BODIES, jd, geo, ps, None).unwrap();
        assert_eq!(lngs.len(), 13);
        for (l, p) in lngs.iter().zip(&full) {
            // The same longitude as a full position, and as the single-body method.
            assert_eq!(l.body, p.body);
            assert!(math::diff_deg(l.lng, p.lng).abs() < 1e-10, "{} {ps}: {} vs {}", l.body, l.lng, p.lng);
            assert!(math::diff_deg(l.lng, eph.calc_longitude(l.body, jd, geo, ps, None).unwrap()).abs() < 1e-10);
        }
        // The ascendant listed as a body is the ascendant of the house functions.
        assert_eq!(lngs[0].body, Body::Ascendant);
        assert!(math::diff_deg(lngs[0].lng, ascendant(jd, geo)).abs() < 1e-6, "{} vs {}", lngs[0].lng, ascendant(jd, geo));
        assert!((full[0].lat, full[0].dist) == (0.0, 0.0) || full[0].lat.abs() < 1e-12);
        assert!(full[0].lng_speed > 150.0, "the ascendant turns about once a day: {}", full[0].lng_speed);
        // Ketu is exactly opposite Rahu, in the list and on its own.
        let (ke, ra) = (lngs[11].lng, lngs[12].lng);
        assert_eq!((lngs[11].body, lngs[12].body), (Body::SouthNode, Body::TrueNode));
        assert!(math::diff_deg(ke + 180.0, ra).abs() < 1e-9);
        assert!(math::diff_deg(ke, eph.calc_longitude(Body::SouthNode, jd, geo, ps, None).unwrap()).abs() < 1e-9);
        assert!((full[11].lng_speed - full[12].lng_speed).abs() < 1e-12 && full[11].lat == -full[12].lat);
    }
    // In either order, and in all three coordinate systems at once, the pair stays opposite.
    let both = eph.calc_positions_ext(&[Body::TrueNode, Body::SouthNode, Body::MeanSouthNode, Body::MeanNode], jd, geo, Perspective::Topocentric).unwrap();
    for (a, b) in [(&both[0], &both[1]), (&both[3], &both[2])] {
        let alone = eph.calc_position_ext(b.body, jd, geo, Perspective::Topocentric).unwrap();
        for (x, y) in [(b.lng, alone.lng), (b.ras, alone.ras), (b.azi, alone.azi)] {
            assert!(math::diff_deg(x, y).abs() < 1e-9, "{} derived {x} vs computed {y}", b.body);
        }
        assert!((b.dcl - alone.dcl).abs() < 1e-9 && (b.alt - alone.alt).abs() < 1e-9 && (b.alt + a.alt).abs() < 1e-12);
        assert!((b.azi_speed - alone.azi_speed).abs() < 1e-6 && (b.alt_speed - alone.alt_speed).abs() < 1e-6);
    }

    // With an ayanamsha every longitude is sidereal, the ascendant included.
    let tropical = eph.calc_longitudes(&Body::BASIC, jd, geo, Perspective::Topocentric, None).unwrap();
    let sidereal = eph.calc_longitudes(&Body::BASIC, jd, geo, Perspective::Topocentric, Some(Ayanamsha::Lahiri)).unwrap();
    let aya = eph.ayanamsha(jd, Ayanamsha::Lahiri).unwrap();
    for (t, s) in tropical.iter().zip(&sidereal) {
        assert!((math::diff_deg(s.lng, t.lng) - aya).abs() < 1e-10 && (0.0..360.0).contains(&s.lng), "{}", t.body);
    }
    assert!((eph.calc_longitude(Body::Moon, jd, geo, Perspective::Geocentric, Some(Ayanamsha::Lahiri)).unwrap() - to_sidereal(eph.calc_longitude(Body::Moon, jd, geo, Perspective::Geocentric, None).unwrap(), aya)).abs() < 1e-10);

    // From the Sun: the adapted set works, and the ascendant has no meaning there.
    let helio = parse_body_list_for("default", Perspective::Heliocentric).unwrap();
    assert_eq!(eph.calc_longitudes(&helio, jd, geo, Perspective::Heliocentric, None).unwrap().len(), 10);
    assert!(matches!(eph.calc_longitude(Body::Ascendant, jd, geo, Perspective::Heliocentric, None), Err(CalcError::NotVisibleFrom { .. })));

    // Longitudes alone cost well under a full set of positions.
    let time = |f: &dyn Fn(f64)| {
        let t = Instant::now();
        (0..3000).for_each(|i| f(jd + i as f64));
        t.elapsed().as_secs_f64() / 3000.0 * 1e6
    };
    time(&|t| drop(eph.calc_positions(&Body::DEFAULT_BODIES, t, geo, CoordSystem::Ecliptic, Perspective::Topocentric).unwrap())); // warm up
    let (lng_only, full) = (
        time(&|t| drop(eph.calc_longitudes(&Body::DEFAULT_BODIES, t, geo, Perspective::Topocentric, None).unwrap())),
        time(&|t| drop(eph.calc_positions(&Body::DEFAULT_BODIES, t, geo, CoordSystem::Ecliptic, Perspective::Topocentric).unwrap())),
    );
    eprintln!("default set: longitudes {lng_only:.1} µs, full positions {full:.1} µs");
    assert!(lng_only < full * 0.8);
}
