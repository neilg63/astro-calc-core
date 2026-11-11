//! Builds a station table from the ephemeris files in a folder.
//! Usage: build_station_table [data-dir] [output-file] [bodies]
//! Defaults: data, <data-dir>/stations.bin, me,ve,ma,ju,sa,ur,ne,pl
use astro_calc_core::*;
use std::time::Instant;
fn main() {
    let mut args = std::env::args().skip(1);
    let dir = args.next().unwrap_or_else(|| "data".into());
    let out = args.next().unwrap_or_else(|| format!("{dir}/stations.bin"));
    let bodies = parse_body_list(&args.next().unwrap_or_else(|| "me,ve,ma,ju,sa,ur,ne,pl".into())).expect("body list");
    let eph = Ephemeris::open_dir(&dir).expect("ephemeris files");
    let t = Instant::now();
    let table = StationTable::build_all(&eph, &bodies).expect("stations");
    let built = t.elapsed();
    table.save(&out).expect("write table");
    let (from, to) = table.range();
    println!(
        "{} stations for {} bodies, {} to {}, computed in {built:.2?}, written to {out} ({} bytes)",
        table.len(),
        table.bodies().count(),
        &JdUt(from).to_iso().unwrap()[..10],
        &JdUt(to).to_iso().unwrap()[..10],
        std::fs::metadata(&out).unwrap().len()
    );
    let table = StationTable::open(&out).expect("read back");
    let jd = JdUt::from_iso("2026-11-11T12:00:00").unwrap().value();
    let n = 1_000_000;
    let t = Instant::now();
    let mut seen = 0usize;
    for i in 0..n {
        let when = from + (i as f64 * 0.1037) % (to - from);
        seen += table.next(bodies[i % bodies.len()], when, 4).map_or(0, |s| s.len());
    }
    println!("table lookup of the next 4 stations: {:.1} ns each ({seen} returned)", t.elapsed().as_secs_f64() * 1e9 / n as f64);
    for &body in &bodies {
        let next: Vec<String> = table.next(body, jd, 2).unwrap().iter().map(|s| format!("{} {}", s.kind().key(), &JdUt(s.jd()).to_iso().unwrap()[..10])).collect();
        println!("{body}: retrograde now {:?}; next {}", table.is_retrograde(body, jd).unwrap(), next.join(", "));
    }
}
