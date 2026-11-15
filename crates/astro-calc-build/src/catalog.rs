//! What can be installed: planetary files, small bodies, body sets and tiers. The built-in
//! catalogue is the default; a config file (see `config`) overrides or extends it.

use astro_calc_core::{Body, DateRange};

pub const DEFAULT_JPL_BSP_URL: &str = "https://ssd.jpl.nasa.gov/ftp/eph/planets/bsp/";
pub const DEFAULT_HORIZONS_API: &str = "https://ssd.jpl.nasa.gov/api/horizons.api";
/// The span JPL Horizons will generate small-body SPK files for. Requests outside it are
/// refused ("START time outside set SPK limits").
pub const DEFAULT_HORIZONS_FIRST: &str = "1600-01-01";
pub const DEFAULT_HORIZONS_LAST: &str = "2500-01-01";

/// Observed size of a Horizons small-body file per year of coverage (5.4 MB for 300 years,
/// 16.2 MB for 900). Distant bodies are usually smaller.
pub const SMALL_BODY_BYTES_PER_YEAR: f64 = 18_000.0;
/// Observed station table size and single-thread build time per year (eight planets).
pub const STATION_BYTES_PER_YEAR: f64 = 580.0;
pub const STATION_SECONDS_PER_YEAR: f64 = 0.013;
/// Planets in the station table by default.
pub const STATION_BODIES: &str = "me,ve,ma,ju,sa,ur,ne,pl";

/// A planetary ephemeris file: the Sun, Moon, planets and Pluto.
#[derive(Debug, Clone, PartialEq)]
pub struct PlanetFile {
    pub key: String,
    /// Name in the data folder.
    pub file: String,
    /// Where to download it; `None` for the JPL base URL plus `file`.
    pub url: Option<String>,
    /// Size in bytes, if known. Unknown sizes are asked of the server when planning.
    pub size: Option<u64>,
    /// Coverage (TDB Julian Days), if known. Unknown ranges are read from the remote
    /// file's index when planning.
    pub range: Option<DateRange>,
    pub note: String,
    pub sha256: Option<String>,
    /// Built by cutting a span out of another planetary file instead of downloading it whole.
    pub excerpt: Option<Excerpt>,
}

/// A span of years of another planetary file (astronomical numbering: -999 is 1000 BCE),
/// from the start of `from` to the start of `to`.
#[derive(Debug, Clone, PartialEq)]
pub struct Excerpt {
    pub of: String,
    pub from: i32,
    pub to: i32,
}

impl Excerpt {
    pub fn range(&self) -> DateRange {
        DateRange { start: year_start_jd(self.from), end: year_start_jd(self.to) }
    }
}

/// A planetary file cut from `of` over `from`..`to`, named after them.
pub fn excerpt_file(of: &PlanetFile, from: i32, to: i32) -> PlanetFile {
    let key = format!("{}_{}_{}", of.key, from, to);
    PlanetFile {
        file: format!("{key}.bsp"),
        key,
        url: None,
        size: None,
        range: Some(DateRange { start: year_start_jd(from), end: year_start_jd(to) }),
        note: format!("{} cut to {} - {}", of.key, year_label(from), year_label(to)),
        sha256: None,
        excerpt: Some(Excerpt { of: of.key.clone(), from, to }),
    }
}

/// An astronomical year for people: -999 is 1000 BCE.
pub fn year_label(y: i32) -> String {
    if y <= 0 { format!("{} BCE", 1 - y) } else { format!("{y} CE") }
}

#[derive(Debug, Clone, PartialEq)]
pub enum BodySource {
    /// Generated on request by JPL Horizons; `command` is its COMMAND, e.g. `2060;`.
    Horizons { command: String },
    /// A ready-made SPK file.
    Url(String),
}

/// A small body with its own SPK file.
#[derive(Debug, Clone, PartialEq)]
pub struct SmallBody {
    pub key: String,
    pub name: String,
    pub file: String,
    /// NAIF id in the file (Horizons: 20 000 000 + minor-planet number).
    pub naif_id: i32,
    pub source: BodySource,
    /// This body's own span (Horizons), in place of the tier's.
    pub from: Option<String>,
    pub to: Option<String>,
    pub size: Option<u64>,
    pub sha256: Option<String>,
}

impl SmallBody {
    /// How astro-calc-core and the API name the body: the two-letter key for the bodies it
    /// knows by name, otherwise `mpc:<number>` or `jpl:<naif id>`.
    pub fn core_key(&self) -> String {
        core_key(self.naif_id)
    }
}

/// The astro-calc-core body key for a NAIF id.
pub fn core_key(naif_id: i32) -> String {
    if let Some((_, name)) = NAMED_SMALL_BODIES.iter().find(|(id, _)| *id == naif_id)
        && let Ok(body) = name.parse::<Body>()
    {
        return body.to_string();
    }
    if naif_id > 20_000_000 { format!("mpc:{}", naif_id - 20_000_000) } else { format!("jpl:{naif_id}") }
}

/// Small bodies astro-calc-core has named variants for.
const NAMED_SMALL_BODIES: [(i32, &str); 6] = [
    (20_000_001, "ceres"),
    (20_000_002, "pallas"),
    (20_000_003, "juno"),
    (20_000_004, "vesta"),
    (20_002_060, "chiron"),
    (20_005_145, "pholus"),
];

/// The span of the station table.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum StationSpan {
    /// Everything the planetary file covers.
    Full,
    /// From the start of the first year to the start of the second (astronomical
    /// numbering: 0 is 1 BC, -2999 is 3000 BC), clamped to the planetary file.
    Years(i32, i32),
}

impl Tier {
    /// The title as a key: `Long range` is `long-range`.
    pub fn slug(&self) -> String {
        self.title.trim().to_ascii_lowercase().replace(' ', "-")
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Tier {
    pub key: String,
    pub title: String,
    /// Key of a planetary file.
    pub planets: String,
    /// Body selection: `all`, `none`, set names and body keys, comma-separated.
    pub bodies: String,
    pub bodies_from: String,
    pub bodies_to: String,
    pub stations: StationSpan,
    pub about: String,
}

/// Everything that can be installed, and from where.
#[derive(Debug, Clone, PartialEq)]
pub struct Catalog {
    pub jpl_bsp_url: String,
    pub horizons_api: String,
    pub horizons_first: String,
    pub horizons_last: String,
    pub planets: Vec<PlanetFile>,
    pub bodies: Vec<SmallBody>,
    pub body_sets: Vec<(String, Vec<String>)>,
    pub tiers: Vec<Tier>,
}

fn planet(key: &str, size: u64, start: f64, end: f64, note: &str) -> PlanetFile {
    PlanetFile {
        key: key.into(),
        file: format!("{key}.bsp"),
        url: None,
        size: Some(size),
        range: Some(DateRange { start, end }),
        note: note.into(),
        sha256: None,
        excerpt: None,
    }
}

/// A small body generated by Horizons from its minor-planet number.
pub fn horizons_body(key: &str, name: &str, mpc: u32) -> SmallBody {
    SmallBody {
        key: key.into(),
        name: name.into(),
        file: format!("{key}.bsp"),
        naif_id: 20_000_000 + mpc as i32,
        source: BodySource::Horizons { command: format!("{mpc};") },
        from: None,
        to: None,
        size: None,
        sha256: None,
    }
}

fn tier(key: &str, title: &str, planets: &str, from: &str, to: &str, stations: StationSpan, about: &str) -> Tier {
    Tier {
        key: key.into(),
        title: title.into(),
        planets: planets.into(),
        bodies: "default".into(),
        bodies_from: from.into(),
        bodies_to: to.into(),
        stations,
        about: about.into(),
    }
}

impl Catalog {
    pub fn builtin() -> Self {
        let set = |name: &str, keys: &[&str]| (name.to_string(), keys.iter().map(|k| k.to_string()).collect());
        Catalog {
            jpl_bsp_url: DEFAULT_JPL_BSP_URL.into(),
            horizons_api: DEFAULT_HORIZONS_API.into(),
            horizons_first: DEFAULT_HORIZONS_FIRST.into(),
            horizons_last: DEFAULT_HORIZONS_LAST.into(),
            // Coverage as read from each file's own index.
            planets: vec![
                planet("de440s", 32_726_016, 2_396_752.5, 2_506_352.5, "DE440, short span"),
                planet("de440", 119_799_808, 2_287_184.5, 2_688_976.5, "DE440, full span"),
                planet("de441", 3_307_878_400, -3_100_015.5, 8_000_016.5, "DE441, the long-span companion to DE440"),
                planet(
                    "de442",
                    119_771_136,
                    2_287_184.5,
                    2_688_976.5,
                    "DE442 (2026), DE440's successor; the test fixtures assume DE440",
                ),
                PlanetFile {
                    key: "de441-ancient".into(),
                    note: "DE441 cut to 1000 BCE - 3000 CE: only those records are fetched from de441.bsp".into(),
                    ..excerpt_file(&planet("de441", 0, 0.0, 0.0, ""), -999, 3000)
                },
            ],
            bodies: vec![
                horizons_body("ceres", "Ceres", 1),
                horizons_body("pallas", "Pallas", 2),
                horizons_body("juno", "Juno", 3),
                horizons_body("vesta", "Vesta", 4),
                horizons_body("chiron", "Chiron", 2060),
                horizons_body("pholus", "Pholus", 5145),
            ],
            body_sets: vec![
                set("default", &["ceres", "pallas", "juno", "vesta", "chiron", "pholus"]),
                set("main", &["ceres", "pallas", "juno", "vesta"]),
                set("centaurs", &["chiron", "pholus"]),
            ],
            tiers: vec![
                tier(
                    "1",
                    "Modern",
                    "de440s",
                    "1850-01-01",
                    "2150-01-01",
                    StationSpan::Full,
                    "Births and events from the mid-19th to mid-22nd century. Smallest download.",
                ),
                tier(
                    "2",
                    "Historical",
                    "de440",
                    DEFAULT_HORIZONS_FIRST,
                    DEFAULT_HORIZONS_LAST,
                    StationSpan::Full,
                    "The Renaissance to the 27th century. Small bodies are limited by Horizons to 1600-2500.",
                ),
                tier(
                    "3",
                    "Ancient",
                    "de441-ancient",
                    DEFAULT_HORIZONS_FIRST,
                    DEFAULT_HORIZONS_LAST,
                    StationSpan::Full,
                    "Classical antiquity to the 30th century, from DE441 without downloading all of it. \
                     About as far as the ΔT and precession models stay sound. Small bodies are limited \
                     by Horizons to 1600-2500.",
                ),
                tier(
                    "4",
                    "Long range",
                    "de441",
                    DEFAULT_HORIZONS_FIRST,
                    DEFAULT_HORIZONS_LAST,
                    StationSpan::Years(-2999, 3000),
                    "Antiquity and the far future. Beyond a few thousand years from today the ΔT and \
                     precession models dominate the error, not the ephemeris. The station table defaults \
                     to 3000 BC-AD 3000.",
                ),
            ],
        }
    }

    pub fn planet(&self, key: &str) -> Option<&PlanetFile> {
        let key = key.trim().trim_end_matches(".bsp").to_ascii_lowercase();
        self.planets.iter().find(|p| p.key == key)
    }

    /// Where a planetary file is downloaded from; for an excerpt, the file it is cut from.
    pub fn planet_url(&self, p: &PlanetFile) -> String {
        if let Some(src) = p.excerpt.as_ref().and_then(|e| self.planet(&e.of))
            && src.excerpt.is_none()
        {
            return self.planet_url(src);
        }
        p.url.clone().unwrap_or_else(|| format!("{}/{}", self.jpl_bsp_url.trim_end_matches('/'), p.file))
    }

    /// Size of a planetary file: declared, or for an excerpt its share of the source.
    pub fn planet_estimate(&self, p: &PlanetFile) -> Option<f64> {
        if let Some(s) = p.size {
            return Some(s as f64);
        }
        let e = p.excerpt.as_ref()?;
        let src = self.planet(&e.of)?;
        let (size, full) = (src.size? as f64, src.range?);
        let part = e.range().intersect(&full)?;
        Some(size * part.days() / full.days())
    }

    pub fn body(&self, key: &str) -> Option<&SmallBody> {
        let key = key.trim().to_ascii_lowercase();
        self.bodies.iter().find(|b| b.key == key)
    }

    pub fn body_set(&self, name: &str) -> Option<&[String]> {
        let name = name.trim().to_ascii_lowercase();
        self.body_sets.iter().find(|(n, _)| *n == name).map(|(_, v)| v.as_slice())
    }

    /// A tier by key (`3`) or by title (`ancient`, `long-range`).
    pub fn tier(&self, key: &str) -> Option<&Tier> {
        let key = key.trim().to_ascii_lowercase().replace(' ', "-");
        self.tiers.iter().find(|t| t.key == key).or_else(|| self.tiers.iter().find(|t| t.slug() == key))
    }

    /// Resolves `all`, `none`, set names and body keys (comma-separated) to bodies, in
    /// catalogue order without repeats.
    pub fn select_bodies(&self, selection: &str) -> Result<Vec<&SmallBody>, String> {
        let mut keys: Vec<String> = Vec::new();
        for token in selection.split(',').map(|t| t.trim().to_ascii_lowercase()).filter(|t| !t.is_empty()) {
            match token.as_str() {
                "none" => {}
                "all" => keys.extend(self.bodies.iter().map(|b| b.key.clone())),
                t => {
                    if let Some(set) = self.body_set(t) {
                        keys.extend(set.iter().cloned());
                    } else if self.body(t).is_some() {
                        keys.push(t.to_string());
                    } else {
                        return Err(format!(
                            "unknown body or set '{t}' (bodies: {}; sets: {})",
                            self.bodies.iter().map(|b| b.key.as_str()).collect::<Vec<_>>().join(", "),
                            self.body_sets.iter().map(|(n, _)| n.as_str()).collect::<Vec<_>>().join(", ")
                        ));
                    }
                }
            }
        }
        Ok(self.bodies.iter().filter(|b| keys.contains(&b.key)).collect())
    }

    /// The span to generate a Horizons body over: its own, else the tier's.
    pub fn body_span(&self, body: &SmallBody, tier_from: &str, tier_to: &str) -> (String, String) {
        (body.from.clone().unwrap_or_else(|| tier_from.into()), body.to.clone().unwrap_or_else(|| tier_to.into()))
    }

    /// Checks that everything refers to something that exists. Lists every problem.
    pub fn validate(&self) -> Result<(), String> {
        let mut errs = Vec::new();
        let horizons_span = |label: &str, from: &str, to: &str, errs: &mut Vec<String>| {
            for (what, d) in [("from", from), ("to", to)] {
                if parse_date(d).is_none() {
                    errs.push(format!("{label}: {what} = \"{d}\" is not a YYYY-MM-DD date"));
                }
            }
            if parse_date(from).is_some() && parse_date(to).is_some() {
                if to <= from {
                    errs.push(format!("{label}: the span ends before it starts ({from} to {to})"));
                }
                if from < self.horizons_first.as_str() || to > self.horizons_last.as_str() {
                    errs.push(format!(
                        "{label}: Horizons only generates small-body files from {} to {} (sources.horizons_first/last)",
                        self.horizons_first, self.horizons_last
                    ));
                }
            }
        };
        for d in [&self.horizons_first, &self.horizons_last] {
            if parse_date(d).is_none() {
                errs.push(format!("sources: \"{d}\" is not a YYYY-MM-DD date"));
            }
        }
        for p in &self.planets {
            check_sha(&format!("planets.{}", p.key), &p.sha256, &mut errs);
            if !p.file.to_ascii_lowercase().ends_with(".bsp") {
                errs.push(format!("planets.{}: file \"{}\" must end in .bsp to be loaded", p.key, p.file));
            }
            if let Some(e) = &p.excerpt {
                match self.planet(&e.of) {
                    None => errs.push(format!("planets.{}.excerpt: unknown planetary file '{}'", p.key, e.of)),
                    Some(src) if src.excerpt.is_some() => {
                        errs.push(format!("planets.{}.excerpt: '{}' is itself an excerpt; cut from the full file", p.key, e.of))
                    }
                    Some(src) if src.file.eq_ignore_ascii_case(&p.file) => {
                        errs.push(format!("planets.{}: an excerpt needs a file name of its own", p.key))
                    }
                    _ => {}
                }
                if e.to <= e.from {
                    errs.push(format!("planets.{}.excerpt: to ({}) must come after from ({})", p.key, e.to, e.from));
                }
            }
        }
        for b in &self.bodies {
            let label = format!("bodies.{}", b.key);
            check_sha(&label, &b.sha256, &mut errs);
            if !b.file.to_ascii_lowercase().ends_with(".bsp") {
                errs.push(format!("{label}: file \"{}\" must end in .bsp to be loaded", b.file));
            }
            if self.planets.iter().any(|p| p.file.eq_ignore_ascii_case(&b.file)) {
                errs.push(format!("{label}: file \"{}\" is also a planetary file's name", b.file));
            }
            if self.bodies.iter().filter(|o| o.file.eq_ignore_ascii_case(&b.file)).count() > 1 {
                errs.push(format!("{label}: file \"{}\" is used by another body", b.file));
            }
            if b.naif_id <= 0 {
                errs.push(format!("{label}: naif_id must be positive"));
            }
            if matches!(b.source, BodySource::Horizons { .. }) && (b.from.is_some() || b.to.is_some()) {
                let (f, t) = self.body_span(b, &self.horizons_first, &self.horizons_last);
                horizons_span(&label, &f, &t, &mut errs);
            }
        }
        for (name, keys) in &self.body_sets {
            if matches!(name.as_str(), "all" | "none") || self.body(name).is_some() {
                errs.push(format!("body_sets.{name}: the name is taken by a body or by all/none"));
            }
            for k in keys {
                if self.body(k).is_none() {
                    errs.push(format!("body_sets.{name}: unknown body '{k}'"));
                }
            }
        }
        for t in &self.tiers {
            let label = format!("tiers.{}", t.key);
            if self.planet(&t.planets).is_none() {
                errs.push(format!(
                    "{label}: unknown planetary file '{}' (known: {})",
                    t.planets,
                    self.planets.iter().map(|p| p.key.as_str()).collect::<Vec<_>>().join(", ")
                ));
            }
            if let Err(e) = self.select_bodies(&t.bodies) {
                errs.push(format!("{label}.bodies: {e}"));
            }
            horizons_span(&label, &t.bodies_from, &t.bodies_to, &mut errs);
            if let StationSpan::Years(a, b) = t.stations
                && b <= a
            {
                errs.push(format!("{label}.stations: to ({b}) must come after from ({a})"));
            }
        }
        if errs.is_empty() { Ok(()) } else { Err(errs.join("\n")) }
    }
}

fn check_sha(label: &str, sha: &Option<String>, errs: &mut Vec<String>) {
    if let Some(s) = sha
        && !(s.len() == 64 && s.chars().all(|c| c.is_ascii_hexdigit()))
    {
        errs.push(format!("{label}: sha256 must be 64 hexadecimal digits"));
    }
}

/// Julian Day of 1 January 00:00 of a proleptic Gregorian year (astronomical numbering).
pub fn year_start_jd(year: i32) -> f64 {
    date_jd(year as i64, 1, 1)
}

/// Julian Day at 00:00 of a proleptic Gregorian date (astronomical year numbering).
pub fn date_jd(y: i64, m: i64, d: i64) -> f64 {
    let a = (14 - m).div_euclid(12);
    let y = y + 4800 - a;
    let m = m + 12 * a - 3;
    let jdn = d + (153 * m + 2).div_euclid(5) + 365 * y + y.div_euclid(4) - y.div_euclid(100) + y.div_euclid(400) - 32045;
    jdn as f64 - 0.5
}

/// Parses `YYYY-MM-DD`, with a leading `-` for years before 1 BC (astronomical numbering),
/// to a Julian Day at 00:00.
pub fn parse_date(s: &str) -> Option<f64> {
    let (neg, rest) = match s.strip_prefix('-') {
        Some(r) => (true, r),
        None => (false, s),
    };
    let mut parts = rest.split('-');
    let (y, m, d) = (parts.next()?, parts.next()?, parts.next()?);
    if parts.next().is_some() || y.len() < 4 || m.len() != 2 || d.len() != 2 {
        return None;
    }
    let all_digits = |p: &str| p.chars().all(|c| c.is_ascii_digit());
    if !(all_digits(y) && all_digits(m) && all_digits(d)) {
        return None;
    }
    let (y, m, d): (i64, i64, i64) = (y.parse().ok()?, m.parse().ok()?, d.parse().ok()?);
    if !(1..=12).contains(&m) || !(1..=31).contains(&d) {
        return None;
    }
    Some(date_jd(if neg { -y } else { y }, m, d))
}

/// Years between two `YYYY-MM-DD` dates, for size estimates.
pub fn years_between(from: &str, to: &str) -> f64 {
    match (parse_date(from), parse_date(to)) {
        (Some(a), Some(b)) => ((b - a) / 365.25).max(0.0),
        _ => 0.0,
    }
}

/// Estimated size of a body's file: declared, or from the span for Horizons.
pub fn body_estimate(body: &SmallBody, span: &(String, String)) -> Option<f64> {
    body.size.map(|s| s as f64).or_else(|| match body.source {
        BodySource::Horizons { .. } => Some(years_between(&span.0, &span.1) * SMALL_BODY_BYTES_PER_YEAR),
        BodySource::Url(_) => None,
    })
}

/// Binary size with the unit JPL's file browser uses.
pub fn fmt_bytes(n: f64) -> String {
    const K: f64 = 1024.0;
    if n >= K * K * K {
        format!("{:.2} GiB", n / (K * K * K))
    } else if n >= K * K {
        format!("{:.1} MiB", n / (K * K))
    } else {
        format!("{:.0} KiB", n / K)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn julian_days() {
        assert_eq!(date_jd(2000, 1, 1), 2_451_544.5);
        assert_eq!(date_jd(1858, 11, 17), 2_400_000.5);
        // JPL file starts, as read from their indexes.
        assert_eq!(date_jd(1849, 12, 26), 2_396_752.5);
        assert_eq!(date_jd(1549, 12, 31), 2_287_184.5);
        assert_eq!(date_jd(-4713, 11, 24), -0.5); // JD 0 at noon, proleptic Gregorian
        assert_eq!(parse_date("-13200-05-06"), Some(-3_100_015.5));
        assert_eq!(parse_date("2000-01-01"), Some(2_451_544.5));
        assert_eq!(parse_date("2000-1-01"), None);
        assert_eq!(parse_date("2000-13-01"), None);
    }

    #[test]
    fn builtin_is_valid() {
        let c = Catalog::builtin();
        c.validate().unwrap();
        assert_eq!(c.planet("DE440.bsp").unwrap().key, "de440");
        assert_eq!(c.planet_url(c.planet("de440").unwrap()), "https://ssd.jpl.nasa.gov/ftp/eph/planets/bsp/de440.bsp");
        assert_eq!(c.select_bodies("default").unwrap().len(), 6);
        assert_eq!(c.select_bodies("centaurs,ceres,chiron").unwrap().len(), 3);
        assert!(c.select_bodies("none").unwrap().is_empty());
        assert!(c.select_bodies("eris").is_err());
        assert_eq!(c.tier("ancient").unwrap().key, "3");
        assert_eq!(c.tier("Long Range").unwrap().key, "4");
        assert_eq!(c.tier("2").unwrap().title, "Historical");
    }

    #[test]
    fn ancient_excerpt() {
        let c = Catalog::builtin();
        let p = c.planet("de441-ancient").unwrap();
        assert_eq!(p.file, "de441_-999_3000.bsp", "same name as --planets de441 --planets-from -999 --planets-to 3000");
        assert_eq!(c.planet_url(p), "https://ssd.jpl.nasa.gov/ftp/eph/planets/bsp/de441.bsp");
        assert_eq!(p.range.unwrap().start_date(), "-0999-01-01");
        // About 104 MiB per thousand years of DE441.
        let mib = c.planet_estimate(p).unwrap() / 1024.0 / 1024.0;
        assert!((400.0..430.0).contains(&mib), "{mib}");
        assert_eq!(year_label(-999), "1000 BCE");
        let cut = excerpt_file(c.planet("de441").unwrap(), -499, 2500);
        assert_eq!((cut.key.as_str(), cut.file.as_str()), ("de441_-499_2500", "de441_-499_2500.bsp"));
    }

    #[test]
    fn core_keys() {
        assert_eq!(core_key(20_000_001), "ce");
        assert_eq!(core_key(20_136_199), "mpc:136199");
        assert_eq!(core_key(-12), "jpl:-12");
    }
}
