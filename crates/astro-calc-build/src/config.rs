//! The config file: overrides and additions to the built-in catalogue.
//!
//! Looked for in this order: `--config`, `$ASTRO_CALC_BUILD_CONFIG`, `./astro-calc-build.toml`,
//! `$XDG_CONFIG_HOME/astro-calc-build/config.toml` (default `~/.config/...`). Entries are
//! merged by key: a key that exists is changed field by field, a new key is added.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use anyhow::{Context, Result, bail};
use astro_calc_core::DateRange;
use serde::{Deserialize, Serialize};

use crate::catalog::{BodySource, Catalog, Excerpt, PlanetFile, SmallBody, StationSpan, Tier, horizons_body, parse_date};

pub const FILE_NAME: &str = "astro-calc-build.toml";

#[derive(Debug, Default, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct RawConfig {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub sources: Option<RawSources>,
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub planets: BTreeMap<String, RawPlanet>,
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub bodies: BTreeMap<String, RawBody>,
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub body_sets: BTreeMap<String, Vec<String>>,
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub tiers: BTreeMap<String, RawTier>,
}

#[derive(Debug, Default, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct RawSources {
    #[serde(skip_serializing_if = "Option::is_none")]
    pub jpl_bsp_url: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub horizons_api: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub horizons_first: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub horizons_last: Option<String>,
}

#[derive(Debug, Default, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct RawPlanet {
    #[serde(skip_serializing_if = "Option::is_none")]
    pub file: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub url: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub size: Option<u64>,
    /// `YYYY-MM-DD`, astronomical years (`-13200-05-06`).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub start: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub end: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub note: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub sha256: Option<String>,
    /// Cut this file from another: `{ of = "de441", from = -999, to = 3000 }` (years).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub excerpt: Option<RawExcerpt>,
}

#[derive(Debug, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct RawExcerpt {
    pub of: String,
    pub from: i32,
    pub to: i32,
}

#[derive(Debug, Default, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct RawBody {
    #[serde(skip_serializing_if = "Option::is_none")]
    pub name: Option<String>,
    /// Minor Planet Center number: sets `command` and `naif_id` for Horizons.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub mpc: Option<u32>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub command: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub naif_id: Option<i32>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub url: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub file: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub from: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub to: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub size: Option<u64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub sha256: Option<String>,
}

#[derive(Debug, Default, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct RawTier {
    #[serde(skip_serializing_if = "Option::is_none")]
    pub title: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub planets: Option<String>,
    /// A set name, body keys, `all` or `none`: a string or a list.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub bodies: Option<StringOrList>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub bodies_from: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub bodies_to: Option<String>,
    /// `"full"` or `{ from = YEAR, to = YEAR }`.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub stations: Option<RawStations>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub about: Option<String>,
}

#[derive(Debug, Deserialize, Serialize)]
#[serde(untagged)]
pub enum StringOrList {
    One(String),
    Many(Vec<String>),
}

impl StringOrList {
    fn joined(self) -> String {
        match self {
            StringOrList::One(s) => s,
            StringOrList::Many(v) => v.join(","),
        }
    }
}

#[derive(Debug, Deserialize, Serialize)]
#[serde(untagged)]
pub enum RawStations {
    Keyword(String),
    Years(RawYears),
}

#[derive(Debug, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct RawYears {
    pub from: i32,
    pub to: i32,
}

/// The config file to use, if any: `explicit` (from `--config` or the environment) must
/// exist; the default locations are used only when present.
pub fn find(explicit: Option<&Path>) -> Result<Option<PathBuf>> {
    if let Some(p) = explicit {
        if !p.is_file() {
            bail!("config file {} not found", p.display());
        }
        return Ok(Some(p.to_path_buf()));
    }
    let local = PathBuf::from(FILE_NAME);
    if local.is_file() {
        return Ok(Some(local));
    }
    let base = std::env::var_os("XDG_CONFIG_HOME")
        .filter(|v| !v.is_empty())
        .map(PathBuf::from)
        .or_else(|| std::env::var_os("HOME").map(|h| PathBuf::from(h).join(".config")));
    Ok(base.map(|b| b.join("astro-calc-build").join("config.toml")).filter(|p| p.is_file()))
}

/// The built-in catalogue with the config file applied, and the file used.
pub fn load(explicit: Option<&Path>, disabled: bool) -> Result<(Catalog, Option<PathBuf>)> {
    let mut cat = Catalog::builtin();
    if disabled {
        return Ok((cat, None));
    }
    let Some(path) = find(explicit)? else { return Ok((cat, None)) };
    let text = std::fs::read_to_string(&path).with_context(|| format!("cannot read {}", path.display()))?;
    let raw: RawConfig = toml::from_str(&text).with_context(|| format!("in {}", path.display()))?;
    apply(&mut cat, raw).with_context(|| format!("in {}", path.display()))?;
    cat.validate().map_err(|e| anyhow::anyhow!("in {}:\n{e}", path.display()))?;
    Ok((cat, Some(path)))
}

fn date_range(label: &str, start: Option<String>, end: Option<String>) -> Result<Option<DateRange>> {
    match (start, end) {
        (None, None) => Ok(None),
        (Some(s), Some(e)) => {
            let jd = |d: &str| parse_date(d).with_context(|| format!("{label}: \"{d}\" is not a YYYY-MM-DD date"));
            Ok(Some(DateRange { start: jd(&s)?, end: jd(&e)? }))
        }
        _ => bail!("{label}: give both start and end, or neither"),
    }
}

/// Merges a config into the catalogue.
pub fn apply(cat: &mut Catalog, raw: RawConfig) -> Result<()> {
    if let Some(s) = raw.sources {
        if let Some(v) = s.jpl_bsp_url {
            cat.jpl_bsp_url = v;
        }
        if let Some(v) = s.horizons_api {
            cat.horizons_api = v;
        }
        if let Some(v) = s.horizons_first {
            cat.horizons_first = v;
        }
        if let Some(v) = s.horizons_last {
            cat.horizons_last = v;
        }
    }

    for (key, rp) in raw.planets {
        let key = key.to_ascii_lowercase();
        let label = format!("planets.{key}");
        let idx = match cat.planets.iter().position(|p| p.key == key) {
            Some(i) => i,
            None => {
                cat.planets.push(PlanetFile {
                    key: key.clone(),
                    file: format!("{key}.bsp"),
                    url: None,
                    size: None,
                    range: None,
                    note: "added in the config file".into(),
                    sha256: None,
                    excerpt: None,
                });
                cat.planets.len() - 1
            }
        };
        let p = &mut cat.planets[idx];
        // Another file (or a copy of it elsewhere): the built-in size and dates may not apply.
        if rp.url.is_some() || rp.file.is_some() {
            if rp.size.is_none() {
                p.size = None;
            }
            if rp.start.is_none() {
                p.range = None;
            }
        }
        if let Some(v) = rp.file {
            p.file = v;
        }
        if let Some(v) = rp.url {
            p.url = Some(v);
        }
        if let Some(v) = rp.size {
            p.size = Some(v);
        }
        if let Some(r) = date_range(&label, rp.start, rp.end)? {
            p.range = Some(r);
        }
        if let Some(v) = rp.note {
            p.note = v;
        }
        if let Some(v) = rp.sha256 {
            p.sha256 = Some(v.to_ascii_lowercase());
        }
        if let Some(e) = rp.excerpt {
            let excerpt = Excerpt { of: e.of.to_ascii_lowercase(), from: e.from, to: e.to };
            // The size follows from the source; the dates are the years asked for.
            p.size = rp.size;
            p.range = Some(excerpt.range());
            p.excerpt = Some(excerpt);
        }
    }

    for (key, rb) in raw.bodies {
        let key = key.to_ascii_lowercase();
        let label = format!("bodies.{key}");
        if rb.url.is_some() && rb.command.is_some() {
            bail!("{label}: give url (a ready-made file) or command (Horizons), not both");
        }
        let idx = match cat.bodies.iter().position(|b| b.key == key) {
            Some(i) => i,
            None => {
                let mut b = match rb.mpc {
                    Some(n) => horizons_body(&key, &key, n),
                    None => {
                        let Some(naif_id) = rb.naif_id else {
                            bail!("{label}: a new body needs mpc = <number>, or naif_id with command or url");
                        };
                        let source = match (&rb.url, &rb.command) {
                            (Some(u), _) => BodySource::Url(u.clone()),
                            (_, Some(c)) => BodySource::Horizons { command: c.clone() },
                            _ => bail!("{label}: a new body needs command (Horizons) or url (a ready-made file)"),
                        };
                        SmallBody {
                            key: key.clone(),
                            name: key.clone(),
                            file: format!("{key}.bsp"),
                            naif_id,
                            source,
                            from: None,
                            to: None,
                            size: None,
                            sha256: None,
                        }
                    }
                };
                let mut chars = key.chars();
                b.name = chars.next().map(|c| c.to_uppercase().chain(chars).collect()).unwrap_or_default();
                cat.bodies.push(b);
                cat.bodies.len() - 1
            }
        };
        let b = &mut cat.bodies[idx];
        if let Some(n) = rb.mpc {
            b.naif_id = 20_000_000 + n as i32;
            b.source = BodySource::Horizons { command: format!("{n};") };
        }
        if rb.url.is_some() || rb.command.is_some() || rb.mpc.is_some() {
            b.size = rb.size;
        }
        if let Some(v) = rb.name {
            b.name = v;
        }
        if let Some(v) = rb.command {
            b.source = BodySource::Horizons { command: v };
        }
        if let Some(v) = rb.url {
            b.source = BodySource::Url(v);
        }
        if let Some(v) = rb.naif_id {
            b.naif_id = v;
        }
        if let Some(v) = rb.file {
            b.file = v;
        }
        if let Some(v) = rb.from {
            b.from = Some(v);
        }
        if let Some(v) = rb.to {
            b.to = Some(v);
        }
        if let Some(v) = rb.size {
            b.size = Some(v);
        }
        if let Some(v) = rb.sha256 {
            b.sha256 = Some(v.to_ascii_lowercase());
        }
    }

    for (name, keys) in raw.body_sets {
        let name = name.to_ascii_lowercase();
        let keys: Vec<String> = keys.iter().map(|k| k.to_ascii_lowercase()).collect();
        match cat.body_sets.iter_mut().find(|(n, _)| *n == name) {
            Some((_, v)) => *v = keys,
            None => cat.body_sets.push((name, keys)),
        }
    }

    for (key, rt) in raw.tiers {
        let key = key.to_ascii_lowercase();
        let label = format!("tiers.{key}");
        let idx = match cat.tiers.iter().position(|t| t.key == key) {
            Some(i) => i,
            None => {
                let Some(planets) = rt.planets.clone() else { bail!("{label}: a new tier needs planets = \"<file key>\"") };
                cat.tiers.push(Tier {
                    key: key.clone(),
                    title: key.clone(),
                    planets,
                    bodies: "default".into(),
                    bodies_from: cat.horizons_first.clone(),
                    bodies_to: cat.horizons_last.clone(),
                    stations: StationSpan::Full,
                    about: String::new(),
                });
                cat.tiers.len() - 1
            }
        };
        let t = &mut cat.tiers[idx];
        if let Some(v) = rt.title {
            t.title = v;
        }
        if let Some(v) = rt.planets {
            t.planets = v.to_ascii_lowercase();
        }
        if let Some(v) = rt.bodies {
            t.bodies = v.joined();
        }
        if let Some(v) = rt.bodies_from {
            t.bodies_from = v;
        }
        if let Some(v) = rt.bodies_to {
            t.bodies_to = v;
        }
        if let Some(v) = rt.stations {
            t.stations = match v {
                RawStations::Keyword(k) if k.eq_ignore_ascii_case("full") => StationSpan::Full,
                RawStations::Keyword(k) => bail!("{label}.stations: expected \"full\" or {{ from = YEAR, to = YEAR }}, got \"{k}\""),
                RawStations::Years(y) => StationSpan::Years(y.from, y.to),
            };
        }
        if let Some(v) = rt.about {
            t.about = v;
        }
    }
    Ok(())
}

/// The effective catalogue as a complete config file.
pub fn to_toml(cat: &Catalog, origin: Option<&Path>) -> Result<String> {
    let mut raw = RawConfig {
        sources: Some(RawSources {
            jpl_bsp_url: Some(cat.jpl_bsp_url.clone()),
            horizons_api: Some(cat.horizons_api.clone()),
            horizons_first: Some(cat.horizons_first.clone()),
            horizons_last: Some(cat.horizons_last.clone()),
        }),
        ..Default::default()
    };
    for p in &cat.planets {
        raw.planets.insert(
            p.key.clone(),
            RawPlanet {
                file: Some(p.file.clone()),
                url: p.url.clone(),
                size: p.size,
                // An excerpt's dates follow from its years.
                start: p.range.filter(|_| p.excerpt.is_none()).map(|r| r.start_date()),
                end: p.range.filter(|_| p.excerpt.is_none()).map(|r| r.end_date()),
                note: Some(p.note.clone()),
                sha256: p.sha256.clone(),
                excerpt: p.excerpt.as_ref().map(|e| RawExcerpt { of: e.of.clone(), from: e.from, to: e.to }),
            },
        );
    }
    for b in &cat.bodies {
        let (command, url) = match &b.source {
            BodySource::Horizons { command } => (Some(command.clone()), None),
            BodySource::Url(u) => (None, Some(u.clone())),
        };
        raw.bodies.insert(
            b.key.clone(),
            RawBody {
                name: Some(b.name.clone()),
                mpc: None,
                command,
                naif_id: Some(b.naif_id),
                url,
                file: Some(b.file.clone()),
                from: b.from.clone(),
                to: b.to.clone(),
                size: b.size,
                sha256: b.sha256.clone(),
            },
        );
    }
    for (name, keys) in &cat.body_sets {
        raw.body_sets.insert(name.clone(), keys.clone());
    }
    for t in &cat.tiers {
        raw.tiers.insert(
            t.key.clone(),
            RawTier {
                title: Some(t.title.clone()),
                planets: Some(t.planets.clone()),
                bodies: Some(StringOrList::One(t.bodies.clone())),
                bodies_from: Some(t.bodies_from.clone()),
                bodies_to: Some(t.bodies_to.clone()),
                stations: Some(match t.stations {
                    StationSpan::Full => RawStations::Keyword("full".into()),
                    StationSpan::Years(from, to) => RawStations::Years(RawYears { from, to }),
                }),
                about: Some(t.about.clone()),
            },
        );
    }
    let origin = match origin {
        Some(p) => format!("the built-in catalogue with {} applied", p.display()),
        None => "the built-in catalogue (no config file found)".into(),
    };
    Ok(format!(
        "# astro-calc-build configuration: {origin}.\n\
         # Save as ./{FILE_NAME} or ~/.config/astro-calc-build/config.toml, or pass --config.\n\
         # Keep only what you change: entries merge with the built-ins by key.\n\
         # Dates are YYYY-MM-DD with astronomical years (-2999 is 3000 BC).\n\n{}",
        toml::to_string_pretty(&raw)?
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn with(text: &str) -> Result<Catalog> {
        let mut c = Catalog::builtin();
        apply(&mut c, toml::from_str(text)?)?;
        c.validate().map_err(|e| anyhow::anyhow!(e))?;
        Ok(c)
    }

    #[test]
    fn effective_config_round_trips() {
        let c = with(
            r#"
            [sources]
            jpl_bsp_url = "https://mirror.example.org/jpl"
            [planets.de441]
            url = "https://mirror.example.org/big/de441.bsp"
            [planets.classical]
            excerpt = { of = "de441", from = -499, to = 1600 }
            [bodies.eris]
            mpc = 136199
            [bodies.sedna]
            url = "https://example.org/sedna.bsp"
            naif_id = 20090377
            [body_sets]
            tno = ["eris", "sedna"]
            [tiers.medieval]
            planets = "de441"
            bodies = ["main", "tno"]
            bodies_to = "2100-01-01"
            stations = { from = 500, to = 1700 }
            "#,
        )
        .unwrap();
        let text = to_toml(&c, None).unwrap();
        let mut again = Catalog::builtin();
        apply(&mut again, toml::from_str(&text).unwrap()).unwrap();
        // Additions come back in key order; compare content.
        assert_eq!(again.planets.len(), c.planets.len());
        for p in &c.planets {
            assert_eq!(again.planet(&p.key), Some(p));
        }
        for b in &c.bodies {
            assert_eq!(again.body(&b.key), Some(b), "{}", b.key);
        }
        for t in &c.tiers {
            assert_eq!(again.tier(&t.key), Some(t));
        }
    }

    #[test]
    fn overrides_and_additions() {
        let c = with(
            r#"
            [planets.de441]
            url = "https://mirror.example.org/de441.bsp"
            [bodies.chiron]
            from = "1700-01-01"
            [bodies.eris]
            mpc = 136199
            name = "Eris"
            [tiers.medieval]
            planets = "de441"
            bodies = "main,eris"
            "#,
        )
        .unwrap();
        let p = c.planet("de441").unwrap();
        assert_eq!(c.planet_url(p), "https://mirror.example.org/de441.bsp");
        assert_eq!((p.size, p.range), (None, None), "a different source drops the built-in size and dates");
        assert_eq!(c.body("chiron").unwrap().from.as_deref(), Some("1700-01-01"));
        let eris = c.body("eris").unwrap();
        assert_eq!((eris.naif_id, eris.core_key()), (20_136_199, "mpc:136199".to_string()));
        assert_eq!(eris.source, BodySource::Horizons { command: "136199;".into() });
        let t = c.tier("medieval").unwrap();
        assert_eq!(c.select_bodies(&t.bodies).unwrap().len(), 5);
        assert_eq!(t.stations, StationSpan::Full);
    }

    #[test]
    fn mistakes_are_reported() {
        let err = |text: &str| format!("{:#}", with(text).unwrap_err());
        assert!(err("[planets.de440]\nsizee = 3").contains("unknown field"));
        assert!(err("[tiers.x]\nplanets = \"de999\"").contains("unknown planetary file 'de999'"));
        assert!(err("[tiers.x]\nplanets = \"de440\"\nbodies = \"ceres,eris\"").contains("unknown body or set 'eris'"));
        assert!(err("[bodies.x]\nname = \"X\"").contains("needs mpc"));
        assert!(err("[bodies.chiron]\nfrom = \"1500-01-01\"").contains("Horizons only generates"));
        assert!(err("[body_sets]\nmine = [\"ceres\", \"nope\"]").contains("unknown body 'nope'"));
        assert!(err("[planets.de440]\nsha256 = \"abc\"").contains("64 hexadecimal"));
        assert!(err("[tiers.1]\nstations = \"some\"").contains("expected \"full\""));
        assert!(err("[planets.x]\nexcerpt = { of = \"de999\", from = 1, to = 2 }").contains("unknown planetary file 'de999'"));
        assert!(err("[planets.x]\nexcerpt = { of = \"de441-ancient\", from = 1, to = 2 }").contains("itself an excerpt"));
        assert!(err("[planets.x]\nexcerpt = { of = \"de441\", from = 5, to = 2 }").contains("must come after"));
    }
}
