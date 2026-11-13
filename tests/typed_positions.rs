//! Positions with fields named after their coordinate system.

use std::path::PathBuf;
use std::time::Instant;

use astro_calc_core::*;

fn eph() -> Option<Ephemeris> {
    let bsp = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("data/de440s.bsp");
    bsp.exists().then(|| Ephemeris::open(bsp).unwrap())
}

#[test]
fn three_systems_agree_with_the_general_method() {
    let Some(eph) = eph() else { return };
    let geo = GeoPos::new(56.1, -3.6, 0.0).unwrap();
    let jd = 2_450_094.975_694_444_5;
    let bodies = parse_body_list("su,mo,me,ma,ju,ra,ml,kr,hip:65474").unwrap();
    for ps in [Perspective::Geocentric, Perspective::Topocentric] {
        let all = eph.calc_positions_ext(&bodies, jd, geo, ps).unwrap();
        let ecl = eph.calc_positions(&bodies, jd, geo, CoordSystem::Ecliptic, ps).unwrap();
        let equ = eph.calc_positions(&bodies, jd, geo, CoordSystem::Equatorial, ps).unwrap();
        let hor = eph.calc_positions(&bodies, jd, geo, CoordSystem::Horizontal, ps).unwrap();
        for (i, x) in all.iter().enumerate() {
            assert_eq!(x.body, bodies[i]);
            // Each named field is the matching angle of the general method.
            assert_eq!((x.lng, x.lat, x.lng_speed, x.lat_speed), (ecl[i].lng, ecl[i].lat, ecl[i].lng_speed, ecl[i].lat_speed));
            assert_eq!((x.ras, x.dcl, x.ras_speed, x.dcl_speed), (equ[i].lng, equ[i].lat, equ[i].lng_speed, equ[i].lat_speed));
            assert_eq!((x.azi, x.alt, x.azi_speed, x.alt_speed), (hor[i].lng, hor[i].lat, hor[i].lng_speed, hor[i].lat_speed));
            // One distance serves all three.
            assert_eq!((x.dist, x.dist_speed), (ecl[i].dist, ecl[i].dist_speed));
            assert_eq!((equ[i].dist, hor[i].dist), (ecl[i].dist, ecl[i].dist));

            // The single-system methods and the parts of the combined one are the same thing.
            let (e, q, h) = (x.ecliptic(), x.equatorial(), x.horizontal());
            assert_eq!(e, ecl[i]);
            assert_eq!(q, eph.calc_position_equ(bodies[i], jd, geo, ps).unwrap());
            assert_eq!(h, eph.calc_position_hor(bodies[i], jd, geo, ps).unwrap());
            assert_eq!(e, eph.calc_position_ecl(bodies[i], jd, geo, ps).unwrap());
            assert_eq!(*x, eph.calc_position_ext(bodies[i], jd, geo, ps).unwrap());
            assert_eq!(BodyPosExt::from_parts(&e, &q, &h), *x);
            // Right ascension is held in degrees; hours are a fifteenth of that.
            assert_eq!((q.ras_hours(), x.ras_hours()), (q.ras / 15.0, x.ras / 15.0));
            assert!((0.0..24.0).contains(&q.ras_hours()));

            // The shared trait reports each system's angles under the common names.
            fn angles(p: &impl BodyPosition) -> (CoordSystem, f64, f64, f64, f64, f64) {
                (p.coord_system(), p.lng(), p.lat(), p.lng_speed(), p.lat_speed(), p.dist())
            }
            assert_eq!(angles(&e), (CoordSystem::Ecliptic, x.lng, x.lat, x.lng_speed, x.lat_speed, x.dist));
            assert_eq!(angles(&q), (CoordSystem::Equatorial, x.ras, x.dcl, x.ras_speed, x.dcl_speed, x.dist));
            assert_eq!(angles(&h), (CoordSystem::Horizontal, x.azi, x.alt, x.azi_speed, x.alt_speed, x.dist));
            assert_eq!(angles(x), angles(&e));
            assert_eq!((q.body(), h.body(), h.dist_speed()), (x.body, x.body, x.dist_speed));
        }
        assert_eq!(eph.calc_positions_equ(&bodies, jd, geo, ps).unwrap(), all.iter().map(|x| x.equatorial()).collect::<Vec<_>>());
        assert_eq!(eph.calc_positions_hor(&bodies, jd, geo, ps).unwrap(), all.iter().map(|x| x.horizontal()).collect::<Vec<_>>());
    }

    // Local Space keeps true ecliptic and equatorial angles here; only the azimuth changes,
    // running from north through east.
    let (topo, local) = (eph.calc_position_ext(Body::Sun, jd, geo, Perspective::Topocentric).unwrap(), eph.calc_position_ext(Body::Sun, jd, geo, Perspective::LocalSpaceHorizon).unwrap());
    assert_eq!((local.lng, local.ras, local.alt), (topo.lng, topo.ras, topo.alt));
    assert!(math::diff_deg(topo.azi + 180.0, local.azi).abs() < 1e-9);

    // Away from the Earth there is no horizon: equatorial works, the others that need it do not.
    assert!(eph.calc_position_equ(Body::Mars, jd, geo, Perspective::Heliocentric).is_ok());
    assert!(eph.calc_position_ecl(Body::Mars, jd, geo, Perspective::Heliocentric).is_ok());
    assert!(matches!(eph.calc_position_hor(Body::Mars, jd, geo, Perspective::Heliocentric), Err(CalcError::NoHorizon(_))));
    assert!(matches!(eph.calc_positions_ext(&bodies, jd, geo, Perspective::Areocentric), Err(CalcError::NoHorizon(_))));

    // All three together cost well under three separate calls.
    let ten = parse_body_list("su,mo,me,ve,ma,ju,sa,ur,ne,pl").unwrap();
    let time = |f: &dyn Fn(f64)| {
        let t = Instant::now();
        (0..3000).for_each(|i| f(jd + i as f64));
        t.elapsed().as_secs_f64()
    };
    let together = time(&|t| drop(eph.calc_positions_ext(&ten, t, geo, Perspective::Topocentric).unwrap()));
    let apart = time(&|t| {
        for cs in [CoordSystem::Ecliptic, CoordSystem::Equatorial, CoordSystem::Horizontal] {
            drop(eph.calc_positions(&ten, t, geo, cs, Perspective::Topocentric).unwrap());
        }
    });
    eprintln!("all three together {:.1} µs, separately {:.1} µs per chart of ten bodies", together / 3000.0 * 1e6, apart / 3000.0 * 1e6);
    assert!(together < apart * 0.75, "together {together} s, apart {apart} s");
}
