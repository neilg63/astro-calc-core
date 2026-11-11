//! Compares computed stations with a table of `{key: [[jd, lng, speed, type], ...]}` where
//! type is 0 retro-start, 1 retro-peak, 2 retro-end, 3 peak. Usage: stations_vs_table <json>
use astro_calc_core::*;
fn main() {
    let path = std::env::args().nth(1).expect("path to table json");
    let table: serde_json::Value = serde_json::from_str(&std::fs::read_to_string(path).unwrap()).unwrap();
    let eph = Ephemeris::open("data/de440s.bsp").unwrap().with_delta_t_model(DeltaTModel::AstroCalcApi);
    let kinds = [StationKind::RetroStart, StationKind::RetroPeak, StationKind::RetroEnd, StationKind::Peak];
    println!("body kind          n  missing   median    p95     max (minutes) | max speed diff | max lng diff (\")");
    for (key, rows) in table.as_object().unwrap() {
        let rows: Vec<Vec<f64>> = rows.as_array().unwrap().iter().map(|r| r.as_array().unwrap().iter().map(|v| v.as_f64().unwrap()).collect()).collect();
        if rows.is_empty() {
            continue;
        }
        let body: Body = key.parse().unwrap();
        let (from, to) = (rows[0][0] - 5.0, rows[rows.len() - 1][0] + 5.0);
        let ours = eph.stations(body, from, to).unwrap();
        let theirs_n = rows.len();
        for (k, kind) in kinds.iter().enumerate() {
            let mut diffs = Vec::new();
            let (mut missing, mut dspeed, mut dlng) = (0, 0f64, 0f64);
            for r in rows.iter().filter(|r| r[3] as usize == k) {
                let best = ours.iter().filter(|s| s.kind() == *kind).min_by(|a, b| (a.jd() - r[0]).abs().total_cmp(&(b.jd() - r[0]).abs())).unwrap();
                let d = (best.jd() - r[0]).abs() * 1440.0;
                if d > 3.0 * 1440.0 {
                    missing += 1;
                    continue;
                }
                diffs.push(d);
                dspeed = dspeed.max((best.speed() - r[2]).abs());
                let p = eph.calc_position(body, best.jd(), GeoPos::default(), CoordSystem::Ecliptic, Perspective::Geocentric).unwrap();
                dlng = dlng.max(math::diff_deg(p.lng, r[1]).abs() * 3600.0);
            }
            diffs.sort_by(|a, b| a.total_cmp(b));
            let n = diffs.len();
            if n == 0 { println!("{key:4} {:12} none matched, {missing} missing", kind.key()); continue; }
            println!("{key:4} {:12} {n:4} {missing:6} {:8.1} {:8.1} {:8.1}               | {dspeed:.2e}       | {dlng:.2}", kind.key(), diffs[n / 2], diffs[n * 95 / 100], diffs[n - 1]);
        }
        let extra = ours.iter().filter(|s| s.jd() > from + 5.0 && s.jd() < to - 5.0).count() as i64 - theirs_n as i64;
        println!("{key:4} table {theirs_n} entries, ours in same span {} ({extra:+})", theirs_n as i64 + extra);
    }
}
