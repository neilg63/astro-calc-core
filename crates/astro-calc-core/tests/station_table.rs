//! The pre-computed station table against the live search.

use std::path::PathBuf;

use astro_calc_core::*;

fn eph() -> Option<Ephemeris> {
    let bsp = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../data/de440s.bsp");
    bsp.exists().then(|| Ephemeris::open(bsp).unwrap())
}

#[test]
fn table_round_trips_and_matches_live_search() {
    let Some(eph) = eph() else { return };
    let bodies = [Body::Mercury, Body::Mars, Body::Jupiter, Body::Sun];
    let (from, to) = (2_455_000.0, 2_462_000.0);
    let table = StationTable::build(&eph, &bodies, from, to).unwrap();
    assert!(table.matches(&eph));
    assert_eq!(table.range(), (from, to));
    assert_eq!(table.bodies().collect::<Vec<_>>(), bodies);
    assert!(table.stations(Body::Sun).unwrap().is_empty());
    assert!(table.stations(Body::Venus).is_none());

    // Binary form: 16 bytes a station, and the same table back (speeds to f32 precision).
    let bytes = table.to_bytes().unwrap();
    assert_eq!(bytes.len(), 48 + 32 * bodies.len() + 16 * table.len());
    let back = StationTable::from_bytes(&bytes).unwrap();
    assert_eq!(back.len(), table.len());
    assert!(back.matches(&eph));
    for body in bodies {
        for (a, b) in table.stations(body).unwrap().iter().zip(back.stations(body).unwrap()) {
            assert!(a.kind() == b.kind() && a.jd() == b.jd() && (a.speed() - b.speed()).abs() < 1e-6);
        }
    }
    let path = std::env::temp_dir().join(format!("astro-calc-stations-{}.bin", std::process::id()));
    table.save(&path).unwrap();
    assert_eq!(StationTable::open(&path).unwrap(), back);
    std::fs::remove_file(&path).unwrap();

    // Lookups agree with the live search, and retrograde state with the sign of the speed.
    let mut worst = (0f64, 0f64);
    for body in [Body::Mercury, Body::Mars, Body::Jupiter] {
        for i in 0..40 {
            let jd = from + 900.0 + i as f64 * 131.7;
            let live = eph.next_stations(body, jd, 5, 3).unwrap();
            let found = back.around(body, jd, 3, 5).unwrap();
            assert_eq!(found.len(), live.stations.len(), "{body} at {jd}");
            for (a, b) in found.iter().zip(&live.stations) {
                // Turning points are sharp; a peak is flat-topped and settles less exactly.
                let turning = matches!(a.kind(), StationKind::RetroStart | StationKind::RetroEnd);
                let diff = (a.jd() - b.jd()).abs();
                if turning { worst.0 = worst.0.max(diff) } else { worst.1 = worst.1.max(diff) }
                assert!(a.kind() == b.kind() && diff < if turning { 1e-4 } else { 1e-3 }, "{a:?} vs {b:?}");
            }
            assert_eq!(back.is_retrograde(body, jd), Some(live.is_retrograde()), "{body} at {jd}");
            let set = back.station_set(&eph, body, jd, 5, 3).unwrap();
            assert_eq!(set.speed, live.speed);
            assert_eq!(set.stations, found);
        }
    }
    eprintln!("largest table-vs-live time difference: turning points {:.2} s, peaks {:.1} s", worst.0 * 86_400.0, worst.1 * 86_400.0);
    let all = back.stations(Body::Mercury).unwrap();
    assert_eq!(back.between(Body::Mercury, all[3].jd(), all[9].jd()).unwrap().len(), 7);
    assert_eq!(back.next(Body::Mercury, to - 1.0, 50).unwrap().len(), all.iter().filter(|s| s.jd() > to - 1.0).count());

    // Outside the table: no answer from the table, a live answer from station_set.
    assert!(back.around(Body::Mercury, to + 10.0, 1, 1).is_none());
    assert!(back.is_retrograde(Body::Venus, from + 10.0).is_none());
    assert!(!back.covers(Body::Mercury, from - 1.0) && back.covers(Body::Mercury, from + 1.0));
    assert_eq!(back.station_set(&eph, Body::Mercury, to + 500.0, 4, 0).unwrap().stations.len(), 4);
    assert_eq!(back.station_set(&eph, Body::Venus, from + 500.0, 4, 0).unwrap().stations.len(), 4);

    // A table built under another ΔT model does not match this ephemeris.
    let other = Ephemeris::open(PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../data/de440s.bsp")).unwrap().with_delta_t_model(DeltaTModel::AstroCalcApi);
    assert!(!back.matches(&other));
}

#[test]
fn damaged_tables_are_rejected() {
    assert!(StationTable::from_bytes(&[]).is_err());
    assert!(StationTable::from_bytes(&[0u8; 64]).is_err());
    let Some(eph) = eph() else { return };
    let good = StationTable::build(&eph, &[Body::Mercury], 2_460_000.0, 2_460_800.0).unwrap().to_bytes().unwrap();
    assert!(StationTable::from_bytes(&good).is_ok());
    assert!(StationTable::from_bytes(&good[..good.len() - 5]).is_err());
    for (at, value) in [(8, 9u8), (12, 200), (32, 7), (48, b'?'), (48 + 24, 99), (48 + 32 + 12, 9)] {
        let mut bad = good.clone();
        bad[at] = value;
        assert!(StationTable::from_bytes(&bad).is_err(), "byte {at} = {value} accepted");
    }
    // Swapping two records breaks the time order.
    let mut bad = good.clone();
    let r = 48 + 32;
    let (a, b) = (good[r..r + 16].to_vec(), good[r + 16..r + 32].to_vec());
    bad[r..r + 16].copy_from_slice(&b);
    bad[r + 16..r + 32].copy_from_slice(&a);
    assert!(StationTable::from_bytes(&bad).is_err());
    assert!(StationTable::build(&eph, &[Body::Mercury], 2_460_000.0, 2_459_000.0).is_err());
    assert!(StationTable::build(&eph, &[Body::Mercury], 2_300_000.0, 2_460_000.0).is_err());
}
