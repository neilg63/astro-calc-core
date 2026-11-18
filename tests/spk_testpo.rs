//! Verifies the SPK reader against JPL's own reference values (testpo.440).
//! Needs a DE440 planetary file (`de440s.bsp` or `de440.bsp`) and `data/testpo.440`; skipped
//! otherwise. The reference values are DE440's own, so other ephemerides are not checked.

mod common;

use std::path::PathBuf;

use astro_calc_core::Ephemeris;

fn data_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("data")
}

/// testpo body codes → NAIF ids.
fn naif(code: i32) -> Option<i32> {
    Some(match code {
        1..=2 | 4..=9 => code,
        3 => 399,
        10 => 301,
        11 => 10,
        12 => 0,
        13 => 3,
        _ => return None,
    })
}

#[test]
fn de440s_matches_testpo() {
    let dir = data_dir();
    let testpo = dir.join("testpo.440");
    let Some(bsp) = common::planetary_file().filter(|path| path.file_name().is_some_and(|name| name.to_string_lossy().starts_with("de440")) && testpo.exists()) else {
        eprintln!("skipping: no DE440 file with testpo.440 in {}", dir.display());
        return;
    };
    let eph = Ephemeris::open(&bsp).unwrap();
    let (start, end) = eph.coverage(399).unwrap();
    let text = std::fs::read_to_string(&testpo).unwrap();
    let mut checked = 0;
    let mut worst = 0.0f64;
    for line in text.lines().skip_while(|l| !l.starts_with("EOT")).skip(1) {
        let p: Vec<&str> = line.split_whitespace().collect();
        if p.len() != 7 {
            continue;
        }
        let jd: f64 = p[2].parse().unwrap();
        let (t, c, x): (i32, i32, usize) = (p[3].parse().unwrap(), p[4].parse().unwrap(), p[5].parse().unwrap());
        let expected: f64 = p[6].parse().unwrap();
        let (Some(tn), Some(cn)) = (naif(t), naif(c)) else { continue };
        if jd < start + 1.0 || jd > end - 1.0 {
            continue;
        }
        let (pos, vel) = eph.relative(tn, cn, jd).unwrap();
        let got = if x <= 3 { pos[x - 1] } else { vel[x - 4] };
        let err = (got - expected).abs();
        worst = worst.max(err);
        assert!(err < 1e-13, "jd {jd} t{t} c{c} x{x}: got {got}, expected {expected}, err {err:e}");
        checked += 1;
    }
    eprintln!("checked {checked} testpo values, worst error {worst:e} AU");
    assert!(checked > 1000);
}

/// Data passed in as bytes (as on WebAssembly) gives the same results as an opened file.
#[test]
fn from_bytes_matches_opened_file() {
    // The whole file is read into memory: not for the multi-gigabyte long-range files.
    let Some(path) = common::planetary_file().filter(|path| std::fs::metadata(path).is_ok_and(|m| m.len() < 500_000_000)) else {
        return;
    };
    let opened = Ephemeris::open(&path).unwrap();
    let loaded = Ephemeris::from_bytes(std::fs::read(&path).unwrap()).unwrap();
    assert_eq!(opened.targets(), loaded.targets());
    for target in opened.targets() {
        assert_eq!(opened.coverage(target), loaded.coverage(target));
        let jd = 2_451_545.0 + target as f64 * 0.37;
        assert_eq!(opened.barycentric(target, jd).unwrap(), loaded.barycentric(target, jd).unwrap());
    }
    assert!(Ephemeris::from_bytes(vec![0u8; 2048]).is_err());
    assert!(Ephemeris::from_bytes(Vec::new()).is_err());
}
