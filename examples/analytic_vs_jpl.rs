//! Compares the built-in analytical series with the JPL files over the span the files cover:
//! the largest difference in apparent geocentric longitude and latitude for each body.
//!
//! Run with: `cargo run --release --features analytic --example analytic_vs_jpl`

use std::time::Instant;

use astro_calc_core::{Body, CalcEngine, CoordSystem, Ephemeris, GeoPos, Perspective, math};

fn main() {
    let dir = std::env::args().nth(1).unwrap_or_else(|| concat!(env!("CARGO_MANIFEST_DIR"), "/data").to_string());
    let jpl = Ephemeris::open_dir(&dir).expect("data files");
    let series = Ephemeris::analytic();
    let range = jpl.date_range().expect("planetary file");
    println!("JPL files cover {range}");
    let bodies = [Body::Sun, Body::Moon, Body::TrueNode, Body::Mercury, Body::Venus, Body::Mars, Body::Jupiter, Body::Saturn, Body::Uranus, Body::Neptune];
    let (from, to) = (range.start + 5.0, range.end - 5.0);
    let samples = 4000;
    let geo = GeoPos::default();
    println!("{:<8} {:>12} {:>12} {:>14}", "body", "max Δlng ″", "max Δlat ″", "max Δdist au");
    for body in bodies {
        let (mut dl, mut db, mut dd) = (0.0f64, 0.0f64, 0.0f64);
        for i in 0..samples {
            let jd = from + (to - from) * i as f64 / (samples - 1) as f64;
            let a = jpl.calc_position(body, jd, geo, CoordSystem::Ecliptic, Perspective::Geocentric).unwrap();
            let b = series.calc_position(body, jd, geo, CoordSystem::Ecliptic, Perspective::Geocentric).unwrap();
            assert_eq!((a.engine, b.engine), (CalcEngine::Jpl, CalcEngine::Analytic));
            dl = dl.max(math::diff_deg(a.lng, b.lng).abs() * 3600.0);
            db = db.max((a.lat - b.lat).abs() * 3600.0);
            dd = dd.max((a.dist - b.dist).abs());
        }
        println!("{:<8} {:>12.3} {:>12.3} {:>14.2e}", body.to_string(), dl, db, dd);
    }
    let t = Instant::now();
    let n = 200;
    for i in 0..n {
        series.calc_positions(&bodies, from + i as f64, geo, CoordSystem::Ecliptic, Perspective::Geocentric).unwrap();
    }
    println!("analytic: {:?} per chart of {} bodies", t.elapsed() / n, bodies.len());
    let t = Instant::now();
    for i in 0..n {
        jpl.calc_positions(&bodies, from + i as f64, geo, CoordSystem::Ecliptic, Perspective::Geocentric).unwrap();
    }
    println!("jpl:      {:?} per chart of {} bodies", t.elapsed() / n, bodies.len());
}
