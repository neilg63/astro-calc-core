use astro_calc_core::*;
use std::time::Instant;
fn main() {
    let eph = Ephemeris::open("data/de440s.bsp").unwrap();
    let jd = JdUt::from_iso("2026-11-11T12:00:00").unwrap();
    for body in [Body::Mercury, Body::Venus, Body::Mars, Body::Jupiter, Body::Pluto] {
        let t = Instant::now();
        let set = eph.next_stations(body, jd, 4, 0).unwrap();
        let e = t.elapsed();
        let list: Vec<String> = set.stations.iter().map(|s| format!("{} {}", s.kind().key(), &JdUt(s.jd()).to_iso().unwrap()[..10])).collect();
        println!("{:8} speed {:+.4} retro {} | {} | {:?}", body.to_string(), set.speed, set.is_retrograde(), list.join(", "), e);
    }
}
