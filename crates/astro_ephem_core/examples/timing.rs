use astro_ephem_core::*;
use std::time::Instant;
fn main() {
    let eph = Ephemeris::open("data/de440s.bsp").unwrap();
    let bodies: Vec<Body> = parse_body_list("su,mo,me,ve,ma,ju,sa,ur,ne,pl").unwrap();
    let geo: GeoPos = "45.22,11.228".parse().unwrap();
    let t = Instant::now();
    let sets = eph.calc_progress_positions(&bodies, 2451545.0, 10_000.0, 1.0, geo, CoordSystem::Ecliptic, Perspective::Topocentric).unwrap();
    let n = sets.len() * sets[0].items.len();
    let e = t.elapsed();
    println!("{n} topocentric positions with speeds in {:?} ({:.2} µs each, single thread)", e, e.as_secs_f64() * 1e6 / n as f64);
}
