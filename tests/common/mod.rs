//! Finding the installed data files for the integration tests.
//!
//! `astro-calc-build` installs one planetary file, whose name depends on the tier chosen:
//! `de440s.bsp` (1850–2150), `de440.bsp` (1550–2650), a cut of DE441 such as
//! `de441_-999_3000.bsp`, or `de441.bsp` itself. The tests take whichever is there. A test
//! of dates beyond the widest standard file, DE440, asks for the span it needs with
//! [`planetary_covering`] and is skipped unless a longer file is installed.
#![allow(dead_code)]

use std::path::PathBuf;

use astro_calc_core::{DateRange, EphemFile, Ephemeris, scan_dir};

/// What `de440.bsp` covers, as Julian Days: 1549-12-31 to 2650-01-25.
pub const DE440_RANGE: DateRange = DateRange { start: 2_287_184.5, end: 2_688_976.5 };

pub fn root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

pub fn data_dir() -> PathBuf {
    root().join("data")
}

/// The installed planetary file: the longest, should there be more than one. `None` when
/// no data is installed, in which case the tests that need it pass without running.
pub fn planetary_info() -> Option<EphemFile> {
    let files = scan_dir(data_dir()).ok()?;
    files.into_iter().filter(EphemFile::is_planetary).max_by(|a, b| a.range.days().total_cmp(&b.range.days()))
}

/// The path of the installed planetary file.
pub fn planetary_file() -> Option<PathBuf> {
    planetary_info().map(|file| data_dir().join(file.name))
}

/// The planetary file on its own, without the small bodies.
pub fn planetary() -> Option<Ephemeris> {
    planetary_file().map(|path| Ephemeris::open(path).unwrap())
}

/// Every installed file, when a planetary file is among them.
pub fn all_files() -> Option<Ephemeris> {
    planetary_info().map(|_| Ephemeris::open_dir(data_dir()).unwrap())
}

/// Whether the installed planetary file reaches beyond DE440 at either end.
pub fn has_long_range() -> bool {
    planetary_info().is_some_and(|file| file.range.start < DE440_RANGE.start || file.range.end > DE440_RANGE.end)
}

/// The planetary file, if it covers the whole span `from..=to` (Julian Days). Otherwise
/// says why the test is skipped and returns `None`.
pub fn planetary_covering(from: f64, to: f64) -> Option<Ephemeris> {
    let Some(file) = planetary_info() else {
        eprintln!("skipping: no planetary file in {}", data_dir().display());
        return None;
    };
    if file.range.contains(from) && file.range.contains(to) {
        return Some(Ephemeris::open(data_dir().join(&file.name)).unwrap());
    }
    eprintln!("skipping: {} covers {}, and this test needs a longer file (astro-calc-build tier 3 or 4)", file.name, file.range);
    None
}
