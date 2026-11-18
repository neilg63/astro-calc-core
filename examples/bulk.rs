//! Bulk charts and aspect searches over many records, with timings.
//!
//! `cargo run --release --example bulk` (add `--features parallel` to use every core).

use astro_calc_core::*;
use std::time::Instant;

fn main() {
    let eph = Ephemeris::open_dir("data").unwrap();

    // 10,000 records at scattered dates between 1900 and 2100, keyed by a name.
    let mut seed = 42_u64;
    let mut next = move || {
        seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        (seed >> 11) as f64 / (1u64 << 53) as f64
    };
    let records: Vec<ChartInput<String>> = (0..10_000)
        .map(|i| ChartInput::new(format!("Person{i}"), 2_415_020.5 + next() * 73_000.0, GeoPos::new(next() * 120.0 - 60.0, next() * 360.0 - 180.0, 0.0).unwrap()))
        .collect();

    // 1. A chart for every record: the default bodies, ascendant included.
    let t = Instant::now();
    let found = eph.bulk_charts(&records, &ChartOptions::default());
    println!("{} charts of 13 bodies in {:?} ({} failed)", found.charts.len(), t.elapsed(), found.errors.len());

    // 2. One to many: whose Sun, Moon or planets are quincunx (150° ± 2°1') the first person's Mars?
    let quincunx = aspect_orbs(&[(150.0, 2.0 + 1.0 / 60.0)]);
    let criteria = AspectCriteria::new(&[Body::Mars], &Body::LUMINARIES_AND_PLANETS, &quincunx);
    let (source, targets) = records.split_first().unwrap();
    let t = Instant::now();
    let mut found = eph.aspects_one_to_many(source, targets, &criteria).unwrap();
    println!("one to many: {} of {} records match in {:?}", found.matches.len(), targets.len(), t.elapsed());
    found.matches.sort_by(|a, b| a.closest().total_cmp(&b.closest()));
    for m in found.matches.iter().take(3) {
        let hit = m.hits[0];
        println!("  {}: {:?} at {:.3}°, {:+.3}° from exact", m.b, hit.body_b, hit.lng_b, hit.deviation);
    }

    // 3. Within given pairs only.
    let pairs: Vec<_> = records.chunks(2).map(|p| (p[0].clone(), p[1].clone())).collect();
    let both = AspectCriteria::new(&Body::LUMINARIES_AND_PLANETS, &Body::LUMINARIES_AND_PLANETS, &quincunx);
    let t = Instant::now();
    let found = eph.aspects_in_pairs(&pairs, &both);
    println!("pairs: {} of {} pairs match in {:?}", found.matches.len(), pairs.len(), t.elapsed());

    // 4. Every two records of one list: Mars of one quincunx Venus of the other.
    let list = &records[..2_000];
    let criteria = AspectCriteria::new(&[Body::Mars], &[Body::Venus], &quincunx);
    let t = Instant::now();
    let found = eph.aspects_many_to_many(list, &criteria);
    println!("many to many: {} matches among {} ordered pairs in {:?}", found.matches.len(), list.len() * (list.len() - 1), t.elapsed());
}
