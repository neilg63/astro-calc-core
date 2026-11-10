//! Verifies the SPK reader against JPL's own reference values (testpo.440).
//! Needs `data/de440s.bsp` and `data/testpo.440` at the workspace root; skipped otherwise.

use std::path::PathBuf;

use astro_ephem_core::Ephemeris;

fn data_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../data")
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
    let (bsp, testpo) = (dir.join("de440s.bsp"), dir.join("testpo.440"));
    if !bsp.exists() || !testpo.exists() {
        eprintln!("skipping: data files not present in {}", dir.display());
        return;
    }
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
