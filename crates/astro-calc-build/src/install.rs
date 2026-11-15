//! `install`: works out what a tier needs, shows the plan, then downloads, verifies,
//! removes conflicting planetary files and builds the station table.

use std::fs;
use std::io::{BufRead, IsTerminal, Write};
use std::path::{Path, PathBuf};
use std::time::Instant;

use anyhow::{Context, Result, bail};
use astro_calc_core::{DateRange, Ephemeris, JdUt, StationTable, parse_body_list};

use crate::catalog::{
    BodySource, Catalog, PlanetFile, STATION_BYTES_PER_YEAR, STATION_SECONDS_PER_YEAR, SmallBody, StationSpan, body_estimate,
    fmt_bytes, parse_date, year_start_jd,
};
use crate::excerpt;
use crate::net::{self, Expect};

/// Days kept clear of the ends of the planetary file for a station table built over chosen
/// years: positions there need ΔT, which reaches hours to days in the remote past.
const YEARS_EDGE_DAYS: f64 = 10.0;

pub struct Options<'a> {
    pub cat: &'a Catalog,
    pub dir: PathBuf,
    pub planets: PlanetFile,
    /// A local copy of the full file to cut an excerpt from, instead of range requests.
    pub from_file: Option<PathBuf>,
    pub bodies: Vec<&'a SmallBody>,
    /// The tier's span for Horizons bodies without their own.
    pub tier_span: (String, String),
    /// `--bodies-from` / `--bodies-to`: override every Horizons body's span.
    pub span_override: (Option<String>, Option<String>),
    pub stations: Option<StationSpan>,
    pub station_bodies: String,
    pub keep_other_planets: bool,
    pub force: bool,
    pub yes: bool,
    pub dry_run: bool,
}

impl Options<'_> {
    /// The span a body's file should cover: for Horizons the override, the body's own or
    /// the tier's; for a ready-made file only what the config declares.
    pub fn span(&self, b: &SmallBody) -> Option<(String, String)> {
        match b.source {
            BodySource::Horizons { .. } => {
                let (f, t) = self.cat.body_span(b, &self.tier_span.0, &self.tier_span.1);
                Some((self.span_override.0.clone().unwrap_or(f), self.span_override.1.clone().unwrap_or(t)))
            }
            BodySource::Url(_) => b.from.clone().zip(b.to.clone()),
        }
    }

    /// Where an excerpt is read from: the local copy, else the full file's URL.
    fn excerpt_source(&self, agent: &ureq::Agent) -> Option<excerpt::Source> {
        self.planets.excerpt.as_ref()?;
        Some(match &self.from_file {
            Some(p) => excerpt::Source::File(p.clone()),
            None => excerpt::Source::Http { agent: agent.clone(), url: self.cat.planet_url(&self.planets) },
        })
    }
}

#[derive(PartialEq)]
enum Step {
    Fetch,
    Keep,
}

struct BodyPlan<'a> {
    body: &'a SmallBody,
    step: Step,
    span: Option<(String, String)>,
    size: Option<f64>,
}

struct Plan<'a> {
    planets: Step,
    planets_url: String,
    planets_size: Option<u64>,
    planets_range: Option<DateRange>,
    resume_from: u64,
    /// For an excerpt: what to copy, worked out from the full file's index.
    excerpt: Option<excerpt::Plan>,
    bodies: Vec<BodyPlan<'a>>,
    stations: Step,
    remove: Vec<(PathBuf, String)>,
    warnings: Vec<String>,
}

fn span_range(span: &(String, String)) -> Option<DateRange> {
    Some(DateRange { start: parse_date(&span.0)?, end: parse_date(&span.1)? })
}

/// Whether a file already holds `target`, over `range` if one is given (give or take two days).
fn file_has(path: &Path, target: i32, range: Option<DateRange>) -> bool {
    let Ok(eph) = Ephemeris::open(path) else { return false };
    match (eph.coverage(target), range) {
        (Some((a, b)), Some(r)) => a <= r.start + 2.0 && b >= r.end - 2.0,
        (Some(_), None) => true,
        (None, _) => false,
    }
}

fn make_plan<'a>(o: &'a Options<'a>, agent: &ureq::Agent) -> Plan<'a> {
    let mut warnings = Vec::new();
    let p = &o.planets;
    let planets_dest = o.dir.join(&p.file);
    let planets_url = o.cat.planet_url(p);
    let (mut planets_size, mut planets_range) = (p.size, p.range);
    let planets = if !o.force && file_has(&planets_dest, 399, p.range) && (p.range.is_some() || is_planetary(&planets_dest)) {
        if planets_range.is_none() {
            planets_range = Ephemeris::open(&planets_dest).ok().and_then(|e| e.date_range());
        }
        Step::Keep
    } else {
        Step::Fetch
    };
    let mut excerpt_plan = None;
    if planets == Step::Fetch
        && let (Some(e), Some(src)) = (&p.excerpt, o.excerpt_source(agent))
    {
        // Read the full file's index to know exactly what the excerpt will hold.
        let r = e.range();
        match excerpt::plan(&src, r.start, r.end, &e.of) {
            Ok(ep) => {
                planets_size = Some(ep.total);
                planets_range = Some(ep.range);
                excerpt_plan = Some(ep);
            }
            Err(err) => warnings.push(format!("cannot plan the excerpt from {}: {err:#}", src.describe())),
        }
    } else if planets == Step::Fetch && (planets_size.is_none() || planets_range.is_none()) {
        match net::probe(agent, &planets_url) {
            Ok((size, range)) => {
                planets_size = planets_size.or(size);
                planets_range = planets_range.or(range);
            }
            Err(e) => warnings.push(format!("could not read the size and dates of {planets_url}: {e:#}")),
        }
    }
    let resume_from = if o.force { 0 } else { fs::metadata(net::part_path(&planets_dest)).map_or(0, |m| m.len()) };

    let bodies = o
        .bodies
        .iter()
        .map(|&b| {
            let span = o.span(b);
            let range = span.as_ref().and_then(span_range);
            let step = if !o.force && file_has(&o.dir.join(&b.file), b.naif_id, range) { Step::Keep } else { Step::Fetch };
            let mut size = span.as_ref().and_then(|s| body_estimate(b, s)).or(b.size.map(|s| s as f64));
            if let (Step::Fetch, None, BodySource::Url(url)) = (&step, size, &b.source) {
                size = net::probe(agent, url).ok().and_then(|(s, _)| s).map(|s| s as f64);
            }
            BodyPlan { body: b, step, span, size }
        })
        .collect();

    // Two files for one body: open_dir uses whichever it loads last.
    for b in &o.bodies {
        let twins = o.cat.bodies.iter().filter(|other| {
            other.key != b.key
                && other.naif_id == b.naif_id
                && (o.bodies.iter().any(|s| s.key == other.key) || o.dir.join(&other.file).is_file())
        });
        for other in twins {
            if b.key < other.key || !o.bodies.iter().any(|s| s.key == other.key) {
                warnings.push(format!(
                    "{} and {} both hold {} ({}); remove one, as Ephemeris::open_dir uses whichever file name sorts last",
                    b.file,
                    other.file,
                    b.name,
                    b.core_key()
                ));
            }
        }
    }

    let mut remove = Vec::new();
    if let Ok(entries) = fs::read_dir(&o.dir) {
        let body_files: Vec<String> = o.cat.bodies.iter().map(|b| b.file.to_ascii_lowercase()).collect();
        let mut paths: Vec<PathBuf> = entries
            .filter_map(|e| e.ok().map(|e| e.path()))
            .filter(|p| p.extension().is_some_and(|x| x.eq_ignore_ascii_case("bsp")) && *p != planets_dest)
            .collect();
        paths.sort();
        for path in paths {
            let name = path.file_name().unwrap_or_default().to_string_lossy().into_owned();
            match Ephemeris::open(&path) {
                Ok(e) if e.files().first().is_some_and(|f| f.is_planetary()) => {
                    if o.keep_other_planets {
                        warnings.push(format!(
                            "{name} is another planetary ephemeris; open_dir loads both and the later name wins where they overlap"
                        ));
                    } else {
                        remove.push((path, format!("another planetary ephemeris ({}); mixing them gives inconsistent results", e.files()[0].range)));
                    }
                }
                Ok(_) => {
                    if !body_files.contains(&name.to_ascii_lowercase()) {
                        warnings.push(format!("{name} is not in the catalogue; it stays and is loaded with the rest"));
                    }
                }
                Err(err) => warnings.push(format!("{name} cannot be read ({err}); Ephemeris::open_dir will fail until it is removed")),
            }
        }
    }
    let stations = if planets == Step::Keep && remove.is_empty() && !o.force && stations_current(o, planets_range) {
        Step::Keep
    } else {
        Step::Fetch
    };
    Plan { planets, planets_url, planets_size, planets_range, resume_from, excerpt: excerpt_plan, bodies, stations, remove, warnings }
}

fn is_planetary(path: &Path) -> bool {
    Ephemeris::open(path).is_ok_and(|e| e.files().first().is_some_and(|f| f.is_planetary()))
}

/// Whether stations.bin was built from these files, for these bodies, over the span asked for.
fn stations_current(o: &Options, planets_range: Option<DateRange>) -> bool {
    let (Some(span), Some(file)) = (o.stations, planets_range) else { return false };
    let (Ok(eph), Ok(table)) = (Ephemeris::open_dir(&o.dir), StationTable::open(o.dir.join("stations.bin"))) else { return false };
    let Ok(wanted) = parse_body_list(&o.station_bodies) else { return false };
    let r = station_range(file, span);
    let (a, b) = table.range();
    table.matches(&eph) && wanted.iter().all(|w| table.bodies().any(|b| b == *w)) && a <= r.start + 3.0 && b >= r.end - 3.0
}

fn station_range(file: DateRange, span: StationSpan) -> DateRange {
    match span {
        StationSpan::Full => file,
        StationSpan::Years(a, b) => DateRange {
            start: year_start_jd(a).max(file.start + YEARS_EDGE_DAYS),
            end: year_start_jd(b).min(file.end - YEARS_EDGE_DAYS),
        },
    }
}

fn print_plan(o: &Options, plan: &Plan) {
    let action = |s: &Step, what: &str| if *s == Step::Keep { "already installed".to_string() } else { what.to_string() };
    let unknown = || "?".to_string();
    let mut download = 0.0;
    let mut download_unknown = false;
    let dir = fs::canonicalize(&o.dir).unwrap_or_else(|_| std::env::current_dir().map(|c| c.join(&o.dir)).unwrap_or(o.dir.clone()));
    println!("Install into {}{}\n", dir.display(), if o.dir.exists() { "" } else { " (to be created)" });
    let p = &o.planets;
    // An excerpt cut from a local file downloads nothing.
    let local = p.excerpt.is_some() && o.from_file.is_some();
    if plan.planets == Step::Fetch && !local {
        match plan.planets_size {
            Some(s) => download += s.saturating_sub(plan.resume_from) as f64,
            None => download_unknown = true,
        }
    }
    let verb = match (&p.excerpt, local) {
        (Some(e), true) => format!("cut from {}", e.of),
        (Some(e), false) => format!("excerpt of {}", e.of),
        (None, _) => "download".into(),
    };
    println!(
        "  Sun, Moon, planets  {:<16} {}  {:>10}  {}",
        p.file,
        plan.planets_range.map_or("dates unknown".into(), |r| r.to_string()),
        plan.planets_size.map_or_else(unknown, |s| fmt_bytes(s as f64)),
        action(&plan.planets, &verb)
    );
    println!("                      {}", p.note);
    if plan.planets == Step::Fetch {
        match (&p.excerpt, &o.from_file) {
            (Some(_), Some(f)) => println!("                      from the local copy {}", f.display()),
            (Some(_), None) => println!(
                "                      only these records, by range requests to {} (the full file is not downloaded)",
                plan.planets_url
            ),
            (None, _) => println!("                      from {}", plan.planets_url),
        }
    }
    if plan.bodies.is_empty() {
        println!("  Small bodies        none");
    }
    for (i, bp) in plan.bodies.iter().enumerate() {
        let b = bp.body;
        let (source, verb) = match &b.source {
            BodySource::Horizons { command } => (format!("Horizons '{command}'"), "generate"),
            BodySource::Url(url) => (url.clone(), "download"),
        };
        if bp.step == Step::Fetch {
            match bp.size {
                Some(s) => download += s,
                None => download_unknown = true,
            }
        }
        println!(
            "  {:<19} {:<16} {:<24}  {:>10}  {}",
            if i == 0 { "Small bodies" } else { "" },
            b.file,
            bp.span.as_ref().map_or("as published".into(), |(f, t)| format!("{f} to {t}")),
            bp.size.map_or_else(unknown, |s| format!("~{}", fmt_bytes(s))),
            action(&bp.step, verb)
        );
        println!("                      {:<16} {} as {}, from {source}", "", b.name, b.core_key());
    }
    match o.stations {
        Some(span) => {
            let (range, size, secs) = match plan.planets_range {
                Some(file) => {
                    let r = station_range(file, span);
                    let years = (r.days() / 365.25).max(0.0);
                    (r.to_string(), fmt_bytes(years * STATION_BYTES_PER_YEAR), format!("about {:.0} s", years * STATION_SECONDS_PER_YEAR))
                }
                None => ("the planetary file's span".into(), unknown(), "time unknown".into()),
            };
            println!("  Station table       {:<16} {range}  {size:>10}  {}", "stations.bin", action(&plan.stations, &format!("compute locally, {secs}")));
            println!("                      bodies: {}", o.station_bodies);
        }
        None => println!("  Station table       skipped (stations are then computed live, in milliseconds instead of nanoseconds)"),
    }
    for (path, why) in &plan.remove {
        println!("  Remove              {:<16} {}", path.file_name().unwrap_or_default().to_string_lossy(), why);
    }
    for w in &plan.warnings {
        println!("  Note                {w}");
    }
    if download > 0.0 || download_unknown {
        println!("\n  Download about {}{}.", fmt_bytes(download), if download_unknown { " plus files of unknown size" } else { "" });
    }
}

fn confirm(question: &str) -> Result<bool> {
    print!("{question} [y/N] ");
    std::io::stdout().flush()?;
    let mut line = String::new();
    std::io::stdin().lock().read_line(&mut line)?;
    Ok(matches!(line.trim().to_ascii_lowercase().as_str(), "y" | "yes"))
}

/// Asks for a tier on a terminal.
pub fn choose_tier(cat: &Catalog) -> Result<String> {
    if !std::io::stdin().is_terminal() {
        bail!("choose a tier with --tier (see `astro-calc-build tiers`)");
    }
    crate::print_tiers(cat);
    let keys: Vec<&str> = cat.tiers.iter().map(|t| t.key.as_str()).collect();
    print!("Tier to install [{}]: ", keys.join(", "));
    std::io::stdout().flush()?;
    let mut line = String::new();
    std::io::stdin().lock().read_line(&mut line)?;
    let key = line.trim().to_ascii_lowercase();
    if cat.tier(&key).is_none() {
        bail!("no tier '{key}'; expected one of {}", keys.join(", "));
    }
    println!();
    Ok(key)
}

pub fn run(o: Options) -> Result<()> {
    let agent = net::agent();
    let plan = make_plan(&o, &agent);
    print_plan(&o, &plan);
    let nothing_to_do = plan.planets == Step::Keep
        && plan.bodies.iter().all(|b| b.step == Step::Keep)
        && plan.remove.is_empty()
        && (o.stations.is_none() || plan.stations == Step::Keep);
    if nothing_to_do {
        println!("\n  Everything is installed. Use --force to download and rebuild anyway.");
        return Ok(());
    }
    if o.dry_run {
        println!("  Dry run: nothing changed.");
        return Ok(());
    }
    if !o.yes {
        if std::io::stdin().is_terminal() {
            if !confirm("\nProceed?")? {
                println!("Nothing changed.");
                return Ok(());
            }
        } else if !plan.remove.is_empty() {
            bail!("this would remove files; pass --yes to confirm, or --keep-other-planets");
        }
    }
    println!();
    fs::create_dir_all(&o.dir).with_context(|| format!("cannot create {}", o.dir.display()))?;
    let started = Instant::now();

    if plan.planets == Step::Fetch {
        let p = &o.planets;
        let dest = o.dir.join(&p.file);
        if o.force {
            let _ = fs::remove_file(net::part_path(&dest));
        }
        let part = match (&p.excerpt, o.excerpt_source(&agent)) {
            (Some(e), Some(src)) => {
                let r = e.range();
                let ep = match &plan.excerpt {
                    Some(ep) => ep,
                    None => &excerpt::plan(&src, r.start, r.end, &e.of).with_context(|| format!("excerpt from {}", src.describe()))?,
                };
                let part = excerpt::write(&src, ep, &dest)?;
                if let Err(err) = excerpt::check(&part, ep, r.start, r.end) {
                    let _ = fs::remove_file(&part);
                    return Err(err.context(format!("{} rejected", p.file)));
                }
                part
            }
            _ => net::download(&agent, &plan.planets_url, &dest, plan.planets_size)?,
        };
        let n = net::install_part(&part, &dest, Expect::Planetary, p.sha256.as_deref())?;
        println!("  ok  {:<16} {}", p.file, fmt_bytes(n as f64));
    }
    let mut failed = Vec::new();
    for bp in plan.bodies.iter().filter(|b| b.step == Step::Fetch) {
        let b = bp.body;
        let dest = o.dir.join(&b.file);
        let result = match (&b.source, &bp.span) {
            (BodySource::Horizons { command }, Some((from, to))) => net::horizons(&agent, &o.cat.horizons_api, command, from, to, &dest),
            (BodySource::Horizons { .. }, None) => Err(anyhow::anyhow!("no span to generate")),
            (BodySource::Url(url), _) => {
                if o.force {
                    let _ = fs::remove_file(net::part_path(&dest));
                }
                net::download(&agent, url, &dest, b.size)
            }
        }
        .and_then(|part| net::install_part(&part, &dest, Expect::Target(b.naif_id), b.sha256.as_deref()));
        match result {
            Ok(n) => println!("  ok  {:<16} {}", b.file, fmt_bytes(n as f64)),
            Err(e) => {
                println!("  !!  {:<16} {e:#}", b.file);
                failed.push(b.key.as_str());
            }
        }
    }
    for (path, _) in &plan.remove {
        fs::remove_file(path).with_context(|| format!("cannot remove {}", path.display()))?;
        println!("  rm  {}", path.file_name().unwrap_or_default().to_string_lossy());
    }
    if let Some(span) = o.stations
        && plan.stations == Step::Fetch
    {
        let years = match span {
            StationSpan::Full => None,
            StationSpan::Years(a, b) => Some((a, b)),
        };
        build_stations(&o.dir, &o.dir.join("stations.bin"), &o.station_bodies, years)?;
    }
    println!("\nDone in {:.0?}.", started.elapsed());
    crate::info::summary(o.cat, &o.dir)?;
    if !failed.is_empty() {
        bail!("not installed: {}; run the same command again to retry", failed.join(", "));
    }
    Ok(())
}

/// Builds the station table from the files in `dir`, over the shared span of the bodies
/// (`years` = None) or over the given years clamped to the planetary file.
pub fn build_stations(dir: &Path, out: &Path, bodies: &str, years: Option<(i32, i32)>) -> Result<()> {
    let bodies = parse_body_list(bodies).map_err(|e| anyhow::anyhow!("station bodies: {e}"))?;
    let eph = Ephemeris::open_dir(dir).with_context(|| format!("cannot load the files in {}", dir.display()))?;
    let range = eph.date_range().context("the planetary file is missing")?;
    let t = Instant::now();
    print!("  ..  stations.bin: computing ");
    let table = match years {
        None => {
            println!("over {range}");
            StationTable::build_all(&eph, &bodies)?
        }
        Some((a, b)) => {
            let r = station_range(range, StationSpan::Years(a, b));
            if r.end <= r.start {
                bail!("years {a} to {b} are outside the planetary file ({range})");
            }
            println!("over {r}");
            StationTable::build(&eph, &bodies, JdUt(r.start), JdUt(r.end))?
        }
    };
    let part = net::part_path(out);
    table.save(&part)?;
    fs::rename(&part, out)?;
    println!("  ok  stations.bin     {} stations, {}, in {:.1?}", table.len(), fmt_bytes(fs::metadata(out)?.len() as f64), t.elapsed());
    Ok(())
}
