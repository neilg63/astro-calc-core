//! Aspect matching and the progressed-chart period search.

use std::path::PathBuf;

use astro_ephem_core::*;

fn eph() -> Option<Ephemeris> {
    let bsp = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../data/de440s.bsp");
    bsp.exists().then(|| Ephemeris::open(bsp).unwrap())
}

fn targets() -> Vec<AspectOrb> {
    aspect_orbs(&[(0.0, 2.0), (60.0, 1.5), (90.0, 2.0), (120.0, 2.0), (180.0, 2.0)])
}

#[test]
fn instant_matching_follows_caller_definitions() {
    let targets = aspect_orbs(&[(30.0, 5.0), (45.0, 8.0), (360.0 / 7.0, 9.0), (60.0, 8.0), (90.0, 8.0), (120.0, 8.0), (150.0, 10.0)]);
    let found: Vec<_> = match_aspects(192.928202, 249.325729, &targets).collect();
    assert_eq!(found.len(), 2);
    assert_eq!(found[0].index, 2);
    assert_eq!(found[1].index, 3);
    assert_eq!(best_aspect(192.928202, 249.325729, &targets).unwrap().aspect.angle, 60.0);
    // Across the 0° boundary, and angles above 180 fold back.
    assert!((separation(359.0, 10.0) - 11.0).abs() < 1e-12);
    assert!((separation(10.0, 359.0) - 11.0).abs() < 1e-12);
    assert_eq!(AspectOrb::new(240.0, -3.0), AspectOrb { angle: 120.0, orb: 3.0 });
    assert!(best_aspect(0.0, 20.0, &targets).is_none());
}

#[test]
fn progressed_periods_agree_with_dense_sampling() {
    let Some(eph) = eph() else { return };
    let a = ChartSubject::progressed(JdUt::from_iso("1975-03-21T08:30:00").unwrap(), GeoPos::new(56.1, -3.6, 0.0).unwrap());
    let b = ChartSubject::progressed(JdUt::from_iso("1978-11-02T19:15:00").unwrap(), GeoPos::new(48.85, 2.35, 0.0).unwrap());
    let bodies: Vec<Body> = Body::CORE.to_vec();
    let targets = targets();
    let opts = ProgressionOptions::default();
    let (from, to) = (b.birth.value(), b.birth.value() + 80.0 * TROPICAL_YEAR);
    let periods = eph.progressed_aspects(&a, &bodies, &b, &bodies, &targets, from, to, &opts).unwrap();
    assert!(periods.len() > 50, "only {} periods", periods.len());

    let lng = |s: &ChartSubject, body: Body, t: f64| {
        eph.calc_position(body, s.chart_jd(t, opts.year_days), s.geo, CoordSystem::Ecliptic, opts.perspective).unwrap().lng
    };
    let dev = |p: &AspectPeriod, t: f64| (separation(lng(&a, p.body_a, t), lng(&b, p.body_b, t)) - p.aspect.angle).abs();
    let mut exact = 0;
    for p in &periods {
        for t in [p.enter_jd, p.leave_jd].into_iter().flatten() {
            assert!((dev(p, t) - p.aspect.orb).abs() < 2e-5, "{p:?}");
        }
        for &t in &p.exact_jds {
            assert!(dev(p, t) < 2e-5, "{p:?}");
            exact += 1;
        }
        assert!((dev(p, p.min_orb_jd) - p.min_orb).abs() < 2e-5, "{p:?}");
        assert!(p.min_orb <= p.aspect.orb + 1e-9);
        let (s, e) = (p.enter_jd.unwrap_or(from), p.leave_jd.unwrap_or(to));
        assert!(s <= p.min_orb_jd && p.min_orb_jd <= e && dev(p, 0.5 * (s + e)) <= p.aspect.orb);
    }
    assert!(exact > 20);

    // Independent check: every in-orb moment on a fine grid lies in a reported period, and
    // every moment well inside a period is in orb.
    let fine: Vec<f64> = (0..).map(|i| from + i as f64 * 7.0).take_while(|&t| t <= to).collect();
    let (mut hits, mut missed) = (0, 0);
    for (ia, &ba) in bodies.iter().enumerate() {
        let la: Vec<f64> = fine.iter().map(|&t| lng(&a, ba, t)).collect();
        for &bb in &bodies[..] {
            let _ = ia;
            let lb: Vec<f64> = fine.iter().map(|&t| lng(&b, bb, t)).collect();
            for (i, &t) in fine.iter().enumerate() {
                for m in match_aspects(la[i], lb[i], &targets) {
                    hits += 1;
                    let covered = periods.iter().any(|p| {
                        p.body_a == ba
                            && p.body_b == bb
                            && p.aspect == m.aspect
                            && p.enter_jd.unwrap_or(from) - 1e-3 <= t
                            && t <= p.leave_jd.unwrap_or(to) + 1e-3
                    });
                    if !covered {
                        missed += 1;
                    }
                }
            }
        }
    }
    assert!(hits > 1000);
    assert_eq!(missed, 0, "{missed} of {hits} in-orb samples not covered");
}

#[test]
fn natal_chart_stays_fixed_and_bad_ranges_are_rejected() {
    let Some(eph) = eph() else { return };
    let geo = GeoPos::new(56.1, -3.6, 0.0).unwrap();
    let a = ChartSubject::progressed(2_442_492.85, geo);
    let b = ChartSubject::natal(2_443_815.3, geo);
    let opts = ProgressionOptions::default();
    let from = b.birth.value();
    let periods = eph.progressed_aspects(&a, &[Body::Venus], &b, &[Body::Mars], &targets(), from, from + 60.0 * TROPICAL_YEAR, &opts).unwrap();
    let mars = eph.calc_position(Body::Mars, b.birth, geo, CoordSystem::Ecliptic, Perspective::Geocentric).unwrap().lng;
    for p in &periods {
        for &t in &p.exact_jds {
            let venus = eph.calc_position(Body::Venus, a.chart_jd(t, opts.year_days), geo, CoordSystem::Ecliptic, Perspective::Geocentric).unwrap().lng;
            assert!((separation(venus, mars) - p.aspect.angle).abs() < 2e-5);
        }
    }
    assert!(periods.iter().any(|p| !p.exact_jds.is_empty()));
    assert!(matches!(
        eph.progressed_aspects(&a, &[Body::Venus], &b, &[Body::Mars], &targets(), from, from - 1.0, &opts),
        Err(CalcError::InvalidInterval { .. })
    ));
}
