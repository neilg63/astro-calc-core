//! `info`: what a data folder holds and what astro-calc-core can do with it.

use std::path::Path;

use anyhow::Result;
use astro_calc_core::{Body, DateRange, Ephemeris, StationTable};

use crate::catalog::{Catalog, core_key, fmt_bytes};

pub fn summary(cat: &Catalog, dir: &Path) -> Result<()> {
    let has_bsp = std::fs::read_dir(dir)
        .map(|it| it.filter_map(|e| e.ok()).any(|e| e.path().extension().is_some_and(|x| x.eq_ignore_ascii_case("bsp"))))
        .unwrap_or(false);
    println!("\n{}", std::fs::canonicalize(dir).unwrap_or(dir.to_path_buf()).display());
    if !has_bsp {
        println!(
            "  No JPL files. astro-calc-core will use its analytic engine where enabled (the `analytic`\n  \
             feature); otherwise every calculation fails. Run `astro-calc-build install {}`.",
            dir.display()
        );
        return Ok(());
    }
    let eph = match Ephemeris::open_dir(dir) {
        Ok(e) => e,
        Err(e) => {
            println!("  Cannot load: {e}");
            return Ok(());
        }
    };
    let name_of = |id: i32| cat.bodies.iter().find(|b| b.naif_id == id).map_or_else(|| core_key(id), |b| b.name.clone());
    for f in eph.files() {
        let kind = if f.is_planetary() {
            "Sun, Moon, planets".to_string()
        } else {
            f.targets.iter().map(|&id| name_of(id)).collect::<Vec<_>>().join(", ")
        };
        println!("  {:<18} {:>10}  {}  {}", f.name, fmt_bytes(f.size as f64), f.range, kind);
    }
    match eph.date_range() {
        Some(r) => println!("\n  Sun, Moon and planets: {r}"),
        None => println!("\n  No planetary file: the Sun, Moon and planets cannot be computed from JPL data."),
    }
    // Every small body in the files, by the key astro-calc-core and the API accept.
    let mut ids: Vec<i32> = eph.targets().into_iter().filter(|&id| id > 1_000_000).collect();
    ids.sort_unstable();
    for id in ids {
        let key = core_key(id);
        if let Some(r) = key.parse::<Body>().ok().and_then(|b| eph.body_range(b)) {
            println!("  {:<22} {r}  (body key {key})", format!("{}:", name_of(id)));
        }
    }
    let table = dir.join("stations.bin");
    match StationTable::open(&table) {
        Ok(t) if t.matches(&eph) => {
            let (a, b) = t.range();
            println!("  Station table:         {} ({} stations)", DateRange { start: a, end: b }, t.len());
        }
        Ok(_) => println!("  Station table:         built from a different ephemeris; rebuild with `astro-calc-build stations`"),
        Err(_) => println!("  Station table:         none (`astro-calc-build stations` builds one)"),
    }
    Ok(())
}
