//! Bulk charts and aspect searches between charts.

mod common;

use astro_calc_core::*;

fn eph() -> Option<Ephemeris> {
    common::planetary()
}

/// `n` records at scattered dates between 1900 and 2100 and scattered places, unsorted.
fn records(n: usize) -> Vec<ChartInput<usize>> {
    let mut seed = 0x2545_F491_4F6C_DD1D_u64;
    let mut next = move || {
        seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        (seed >> 11) as f64 / (1u64 << 53) as f64
    };
    (0..n).map(|key| ChartInput::new(key, 2_415_020.5 + next() * 73_000.0, GeoPos::new(next() * 120.0 - 60.0, next() * 360.0 - 180.0, 0.0).unwrap())).collect()
}

fn quincunx() -> Vec<AspectOrb> {
    aspect_orbs(&[(150.0, 2.0 + 1.0 / 60.0)])
}

fn lngs(eph: &Ephemeris, bodies: &[Body], item: &ChartInput<usize>) -> Vec<f64> {
    eph.calc_longitudes(bodies, item.jd, item.geo, Perspective::Geocentric, None).unwrap().iter().map(|l| l.lng).collect()
}

/// Whether any body of `a` and any of `b` stand 150° apart within the orb.
fn in_quincunx(a: &[f64], b: &[f64]) -> bool {
    a.iter().any(|&x| b.iter().any(|&y| (separation(x, y) - 150.0).abs() <= 2.0 + 1.0 / 60.0))
}

#[test]
fn charts_keep_input_order_and_keys() {
    let Some(eph) = eph() else { return };
    let mut items = records(200);
    // A date no engine covers, in the middle of the list.
    items.insert(100, ChartInput::new(9999, 0.0, GeoPos::default()));
    let opts = ChartOptions::default();
    let found = eph.bulk_charts(&items, &opts);
    assert_eq!(found.charts.len(), 200);
    assert_eq!(found.errors.len(), 1);
    assert_eq!(found.errors[0].key, 9999);
    assert!(found.errors[0].error.is_out_of_range());
    let keys: Vec<usize> = found.charts.iter().map(|c| c.key).collect();
    assert_eq!(keys, (0..200).collect::<Vec<_>>());
    for chart in found.charts.iter().step_by(17) {
        let item = &items[if chart.key < 100 { chart.key } else { chart.key + 1 }];
        let single = eph.calc_positions(&opts.bodies, item.jd, item.geo, CoordSystem::Ecliptic, opts.perspective).unwrap();
        assert_eq!(chart.positions, single);
        assert_eq!(chart.positions[0].body, Body::Ascendant);
        assert_eq!(chart.engine, CalcEngine::Jpl);
    }
}

#[test]
fn sidereal_charts_are_offset_by_the_ayanamsha() {
    let Some(eph) = eph() else { return };
    let items = records(3);
    let tropical = eph.bulk_charts(&items, &ChartOptions::default());
    let sidereal = eph.bulk_charts(&items, &ChartOptions { ayanamsha: Some(Ayanamsha::Lahiri), ..Default::default() });
    for (t, s) in tropical.charts.iter().zip(&sidereal.charts) {
        let offset = eph.ayanamsha(t.jd, Ayanamsha::Lahiri).unwrap();
        for (pt, ps) in t.positions.iter().zip(&s.positions) {
            assert!(separation(pt.lng - offset, ps.lng) < 1e-9);
        }
    }
}

#[test]
fn one_to_many_returns_only_matching_records() {
    let Some(eph) = eph() else { return };
    let items = records(600);
    let (source, targets) = items.split_first().unwrap();
    let bodies = Body::LUMINARIES_AND_PLANETS;
    let criteria = AspectCriteria::new(&[Body::Mars], &bodies, &quincunx());
    let found = eph.aspects_one_to_many(source, targets, &criteria).unwrap();
    assert!(found.errors.is_empty());
    assert!(!found.matches.is_empty() && found.matches.len() < targets.len() / 2);

    let mars = lngs(&eph, &[Body::Mars], source);
    let expected: Vec<usize> = targets.iter().filter(|t| in_quincunx(&mars, &lngs(&eph, &bodies, t))).map(|t| t.key).collect();
    assert_eq!(found.matches.iter().map(|m| m.b).collect::<Vec<_>>(), expected);
    for m in &found.matches {
        assert_eq!(m.a, source.key);
        for hit in &m.hits {
            assert_eq!((hit.body_a, hit.angle), (Body::Mars, 150.0));
            assert!(hit.deviation.abs() <= hit.orb);
            assert!(((separation(hit.lng_a, hit.lng_b) - 150.0).abs() - hit.deviation.abs()).abs() < 1e-9);
        }
        assert!(m.closest() <= 2.0 + 1.0 / 60.0);
    }
}

#[test]
fn pairs_are_tested_one_against_the_other_only() {
    let Some(eph) = eph() else { return };
    let items = records(400);
    let pairs: Vec<_> = items.chunks(2).map(|p| (p[0], p[1])).collect();
    let bodies = Body::LUMINARIES_AND_PLANETS;
    let criteria = AspectCriteria::new(&bodies, &bodies, &quincunx());
    let found = eph.aspects_in_pairs(&pairs, &criteria);
    let expected: Vec<(usize, usize)> =
        pairs.iter().filter(|(a, b)| in_quincunx(&lngs(&eph, &bodies, a), &lngs(&eph, &bodies, b))).map(|(a, b)| (a.key, b.key)).collect();
    assert!(!expected.is_empty() && expected.len() < pairs.len());
    assert_eq!(found.matches.iter().map(|m| (m.a, m.b)).collect::<Vec<_>>(), expected);

    // A failing record removes its pair and is reported by key.
    let mut broken = pairs.clone();
    broken[0].1 = ChartInput::new(7777, 0.0, GeoPos::default());
    let found = eph.aspects_in_pairs(&broken, &criteria);
    assert_eq!(found.errors.len(), 1);
    assert_eq!(found.errors[0].key, 7777);
    assert!(found.matches.iter().all(|m| m.b != 7777));
}

#[test]
fn many_to_many_covers_every_pair_once() {
    let Some(eph) = eph() else { return };
    let items = records(120);
    let bodies = [Body::Sun, Body::Moon, Body::Venus, Body::Mars];
    let all: Vec<Vec<f64>> = items.iter().map(|i| lngs(&eph, &bodies, i)).collect();

    // Same bodies on both sides: each pair once, the earlier record as A.
    let found = eph.aspects_many_to_many(&items, &AspectCriteria::new(&bodies, &bodies, &quincunx()));
    let mut expected = Vec::new();
    for i in 0..items.len() {
        for j in i + 1..items.len() {
            if in_quincunx(&all[i], &all[j]) {
                expected.push((i, j));
            }
        }
    }
    assert!(!expected.is_empty());
    assert_eq!(found.matches.iter().map(|m| (m.a, m.b)).collect::<Vec<_>>(), expected);

    // Different bodies: both ways round, A being the record whose Mars was used.
    let found = eph.aspects_many_to_many(&items, &AspectCriteria::new(&[Body::Mars], &[Body::Venus], &quincunx()));
    let mut expected = Vec::new();
    for i in 0..items.len() {
        for j in 0..items.len() {
            if i != j && in_quincunx(&all[i][3..], &all[j][2..3]) {
                expected.push((i, j));
            }
        }
    }
    assert!(!expected.is_empty());
    assert_eq!(found.matches.iter().map(|m| (m.a, m.b)).collect::<Vec<_>>(), expected);
    assert!(found.matches.iter().flat_map(|m| &m.hits).all(|h| (h.body_a, h.body_b) == (Body::Mars, Body::Venus)));
}

#[test]
fn progressed_charts_use_a_day_for_a_year() {
    let Some(eph) = eph() else { return };
    let items = records(50);
    let (source, targets) = items.split_first().unwrap();
    let bodies = Body::LUMINARIES_AND_PLANETS;
    let at = 2_460_000.5;
    let orbs = aspect_orbs(&[(0.0, 3.0), (90.0, 3.0), (180.0, 3.0)]);
    let criteria = AspectCriteria::new(&bodies, &bodies, &orbs).progressed_to(at, ChartKind::Progressed, ChartKind::Natal);
    let found = eph.aspects_one_to_many(source, targets, &criteria).unwrap();

    // The same search with the source moved by hand to its progressed date.
    let moved = ChartInput::new(source.key, source.jd.value() + (at - source.jd.value()) / TROPICAL_YEAR, source.geo);
    let natal = eph.aspects_one_to_many(&moved, targets, &AspectCriteria::new(&bodies, &bodies, &orbs)).unwrap();
    assert!(!found.matches.is_empty());
    assert_eq!(found.matches, natal.matches);
    // And it differs from comparing the two birth charts.
    let births = eph.aspects_one_to_many(source, targets, &AspectCriteria::new(&bodies, &bodies, &orbs)).unwrap();
    assert_ne!(found.matches, births.matches);
}
