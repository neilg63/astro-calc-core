//! astro-calc-build: installs the JPL ephemeris files that astro-calc-core reads.

mod catalog;
mod config;
mod excerpt;
mod info;
mod install;
mod net;

use std::path::{Path, PathBuf};

use anyhow::{Result, anyhow, bail};
use clap::{Args, Parser, Subcommand};

use catalog::{Catalog, STATION_BODIES, StationSpan, body_estimate, fmt_bytes, parse_date};

const ABOUT: &str = "Installs the JPL ephemeris files that astro-calc-core reads.";

const LONG_ABOUT: &str = "\
Installs the JPL ephemeris files that astro-calc-core reads.

astro-calc-core is distributed without JPL data. Without it, calculations use the built-in
analytic series (with the `analytic` feature) or fail. This tool downloads a planetary file
from JPL, has JPL Horizons generate files for the small bodies, checks every file by opening
it, and builds the pre-computed station table, all in the data folder you name with
--data-dir (-d): any path on the system. Point the
apps that use the data (e.g. astro-calc-api's ASTRO_DATA_DIR) at the same folder.

Built-in tiers:

  1  modern       1849-2150 CE           de440s.bsp          about 62 MiB in all
  2  historical   1549-2650 CE           de440.bsp           about 207 MiB in all
  3  ancient      1000 BCE-3000 CE       de441_-999_3000.bsp about 510 MiB in all
  4  long-range   13201 BCE-17191 CE     de441.bsp           about 3.2 GiB in all

Tier 3 is cut from DE441: only the records for its span are fetched, by HTTP range requests.
Small bodies are limited by Horizons to 1600-2500 in tiers 2 to 4.

A config file can change the download URLs and file names, and add planetary files, small
bodies, body sets and tiers. It is read from --config, $ASTRO_CALC_BUILD_CONFIG,
./astro-calc-build.toml or ~/.config/astro-calc-build/config.toml, the first found.
`astro-calc-build config` prints the effective configuration as a starting point.";

const EXAMPLES: &str = "\
Examples:
  astro-calc-build                                       Overview
  astro-calc-build install -d /srv/astro/data            Choose a tier interactively
  astro-calc-build install -d /srv/astro/data --tier 2 --dry-run
  astro-calc-build install --data-dir ../astro-calc-core/data --tier 1 --yes
  astro-calc-build install -d ~/ephem --tier ancient                 1000 BCE-3000 CE, ~510 MiB
  astro-calc-build install -d ~/ephem --tier ancient --from-file ~/Downloads/de441.bsp
  astro-calc-build install -d ~/ephem --tier 2 --bodies-from 1800-01-01 --bodies-to 2200-01-01
  astro-calc-build install -d ~/ephem --planets de441 --planets-from -499 --planets-to 2500
  astro-calc-build install -d ~/ephem --tier 4 --bodies none --stations-from -999 --stations-to 2500
  astro-calc-build install -d ~/ephem --tier 2 --bodies main,chiron
  astro-calc-build info -d /srv/astro/data               What is installed, and the dates it covers
  astro-calc-build tiers                                 Tiers, files, bodies and sets in detail
  astro-calc-build config > astro-calc-build.toml        Start a config file

The data folder can also be given as a plain argument (install /srv/astro/data) or by
ASTRO_DATA_DIR.";

#[derive(Parser)]
#[command(name = "astro-calc-build", version, about = ABOUT, long_about = LONG_ABOUT, after_help = EXAMPLES)]
struct Cli {
    /// Config file (default: ./astro-calc-build.toml, then ~/.config/astro-calc-build/config.toml)
    #[arg(long, global = true, value_name = "FILE", env = "ASTRO_CALC_BUILD_CONFIG")]
    config: Option<PathBuf>,

    /// Ignore any config file and use the built-in catalogue
    #[arg(long, global = true, conflicts_with = "config")]
    no_config: bool,

    /// Without a command: an overview of the tool, the tiers and the data folder
    #[command(subcommand)]
    command: Option<Command>,
}

#[derive(Subcommand)]
enum Command {
    /// Download, verify and set up the files for a date-range tier
    #[command(after_help = EXAMPLES)]
    Install(InstallArgs),
    /// List the tiers, planetary files, small bodies and body sets
    Tiers,
    /// Show the files in a data folder and the dates they cover
    Info {
        #[command(flatten)]
        dir: DirArg,
    },
    /// Build (or rebuild) the station table from the files already installed
    Stations(StationsArgs),
    /// Print the effective configuration (built-ins plus any config file) as TOML
    Config,
}

#[derive(Args)]
struct DirArg {
    /// The data folder, any path on the system: where the files go, and where astro-calc-core's
    /// Ephemeris::open_dir and astro-calc-api's ASTRO_DATA_DIR read them. Created if missing.
    /// Default: $ASTRO_DATA_DIR
    #[arg(short = 'd', long = "data-dir", value_name = "PATH")]
    data_dir: Option<PathBuf>,

    /// The data folder, as a plain argument (same as --data-dir)
    #[arg(value_name = "DATA_DIR", conflicts_with = "data_dir")]
    dir: Option<PathBuf>,
}

impl DirArg {
    /// `--data-dir`, the plain argument or `$ASTRO_DATA_DIR`, with `~` expanded, as an
    /// absolute path.
    fn resolve(&self) -> Result<PathBuf> {
        let raw = self
            .data_dir
            .clone()
            .or_else(|| self.dir.clone())
            .or_else(|| std::env::var_os("ASTRO_DATA_DIR").filter(|v| !v.is_empty()).map(PathBuf::from))
            .ok_or_else(|| anyhow!("name the data folder with --data-dir <PATH> (or -d), or set ASTRO_DATA_DIR"))?;
        absolute(raw)
    }
}

/// A path with `~` expanded, made absolute.
fn absolute(raw: PathBuf) -> Result<PathBuf> {
    let path = match raw.to_str() {
        Some("~") => home()?,
        Some(s) if s.starts_with("~/") => home()?.join(&s[2..]),
        _ => raw,
    };
    Ok(if path.is_absolute() { path } else { std::env::current_dir()?.join(path) })
}

fn home() -> Result<PathBuf> {
    std::env::var_os("HOME").map(PathBuf::from).ok_or_else(|| anyhow!("HOME is not set, so ~ cannot be expanded"))
}

#[derive(Args)]
struct InstallArgs {
    #[command(flatten)]
    dir: DirArg,

    /// Tier, by number or name: 1 modern (1849-2150), 2 historical (1549-2650), 3 ancient
    /// (1000 BCE-3000 CE), 4 long-range (13201 BCE-17191 CE), or one in the config file.
    /// Small bodies cover 1600-2500 at most. Asked for when omitted on a terminal
    #[arg(long, short)]
    tier: Option<String>,

    /// Planetary file instead of the tier's: de440s, de440, de441, de442, de441-ancient
    #[arg(long, value_name = "KEY")]
    planets: Option<String>,

    /// Cut the planetary file to start at this year (astronomical: -999 is 1000 BCE). Only the
    /// records needed are fetched, by range requests
    #[arg(long, value_name = "YEAR", allow_negative_numbers = true)]
    planets_from: Option<i32>,

    /// Cut the planetary file to end at the start of this year
    #[arg(long, value_name = "YEAR", allow_negative_numbers = true)]
    planets_to: Option<i32>,

    /// Cut the excerpt from this local copy of the full file (e.g. de441.bsp) instead of using
    /// range requests
    #[arg(long, value_name = "PATH")]
    from_file: Option<PathBuf>,

    /// Small bodies instead of the tier's: `all`, `none`, set names and body keys,
    /// comma-separated (e.g. `main,chiron`)
    #[arg(long, value_name = "LIST")]
    bodies: Option<String>,

    /// First day of every Horizons small-body file, YYYY-MM-DD (Horizons allows 1600-01-01 at the earliest)
    #[arg(long, value_name = "DATE")]
    bodies_from: Option<String>,

    /// Last day of every Horizons small-body file, YYYY-MM-DD (Horizons allows 2500-01-01 at the latest)
    #[arg(long, value_name = "DATE")]
    bodies_to: Option<String>,

    /// Do not build the station table
    #[arg(long)]
    no_stations: bool,

    /// First year of the station table (astronomical numbering: -2999 is 3000 BC). Default: the
    /// start of the planetary file (tier 3: -2999)
    #[arg(long, value_name = "YEAR", allow_negative_numbers = true)]
    stations_from: Option<i32>,

    /// Year the station table ends (at its start). Default: the end of the planetary file (tier 3: 3000)
    #[arg(long, value_name = "YEAR", allow_negative_numbers = true)]
    stations_to: Option<i32>,

    /// Planets in the station table
    #[arg(long, value_name = "LIST", default_value = STATION_BODIES)]
    station_bodies: String,

    /// Leave other planetary files (for example de440s.bsp when installing de440) in the folder
    #[arg(long)]
    keep_other_planets: bool,

    /// Download again even when the files are already installed
    #[arg(long)]
    force: bool,

    /// Do not ask for confirmation
    #[arg(long, short)]
    yes: bool,

    /// Show the plan and stop
    #[arg(long)]
    dry_run: bool,
}

#[derive(Args)]
struct StationsArgs {
    #[command(flatten)]
    dir: DirArg,

    /// Output file. Default: stations.bin in the data folder, where astro-calc-api looks for it
    #[arg(long, short)]
    out: Option<PathBuf>,

    /// Planets in the table
    #[arg(long, value_name = "LIST", default_value = STATION_BODIES)]
    bodies: String,

    /// First year (astronomical numbering: -2999 is 3000 BC). Default: the start of the files
    #[arg(long, value_name = "YEAR", allow_negative_numbers = true)]
    from: Option<i32>,

    /// Year the table ends (at its start). Default: the end of the files
    #[arg(long, value_name = "YEAR", allow_negative_numbers = true)]
    to: Option<i32>,
}

pub fn print_tiers(cat: &Catalog) {
    println!("Tiers\n");
    let w = cat.tiers.iter().map(|t| t.key.len()).max().unwrap_or(1).max(1);
    let pad = " ".repeat(w + 3);
    for t in &cat.tiers {
        let Some(p) = cat.planet(&t.planets) else { continue };
        let bodies = cat.select_bodies(&t.bodies).unwrap_or_default();
        let mut est = 0.0;
        let mut unknown = false;
        for b in &bodies {
            match body_estimate(b, &cat.body_span(b, &t.bodies_from, &t.bodies_to)) {
                Some(s) => est += s,
                None => unknown = true,
            }
        }
        let range = p.range.map_or("dates read on download".into(), |r| r.to_string());
        let size = cat.planet_estimate(p).map_or("size unknown".into(), |s| format!("~{}", fmt_bytes(s)));
        println!("  {:<w$} {:<11} {range}  {} {size}", t.key, t.title, p.file);
        println!(
            "{pad}{:<11} small bodies: {} ({}), {} to {}, about {}{}",
            "",
            t.bodies,
            bodies.len(),
            t.bodies_from,
            t.bodies_to,
            fmt_bytes(est),
            if unknown { " plus files of unknown size" } else { "" }
        );
        if let Some(s) = cat.planet_estimate(p) {
            println!("{pad}{:<11} total about {}", "", fmt_bytes(s + est));
        }
        let stations = match t.stations {
            StationSpan::Full => "station table over the full span".to_string(),
            StationSpan::Years(a, b) => format!("station table {a} to {b}"),
        };
        println!("{pad}{:<11} {stations}", "");
        if !t.about.is_empty() {
            println!("{pad}{:<11} {}", "", t.about);
        }
        println!();
    }
}

fn print_catalogue(cat: &Catalog, origin: Option<&Path>) {
    match origin {
        Some(p) => println!("Config: {}\n", p.display()),
        None => println!("Config: built-in catalogue (no config file)\n"),
    }
    print_tiers(cat);
    println!("Planetary files (any tier can use another with --planets)\n");
    for p in &cat.planets {
        println!(
            "  {:<8} {:<30} {:>10}  {}",
            p.key,
            p.range.map_or("dates read on download".into(), |r| r.to_string()),
            cat.planet_estimate(p).map_or("?".into(), |s| format!("~{}", fmt_bytes(s))),
            p.note
        );
        println!("  {:<8} {}", "", cat.planet_url(p));
    }
    println!("\nSmall bodies (Horizons generates {} to {} at most)\n", cat.horizons_first, cat.horizons_last);
    for b in &cat.bodies {
        let source = match &b.source {
            catalog::BodySource::Horizons { command } => format!("Horizons '{command}'"),
            catalog::BodySource::Url(u) => u.clone(),
        };
        let span = match (&b.from, &b.to) {
            (None, None) => String::new(),
            (f, t) => format!(", {} to {}", f.as_deref().unwrap_or("tier start"), t.as_deref().unwrap_or("tier end")),
        };
        println!("  {:<8} {:<10} {:<12} {}{span}", b.key, b.name, b.core_key(), source);
    }
    println!("\nBody sets (for --bodies; also `all` and `none`)\n");
    for (name, keys) in &cat.body_sets {
        println!("  {:<10} {}", name, keys.join(", "));
    }
}

fn check_horizons_date(cat: &Catalog, label: &str, s: &str) -> Result<String> {
    if parse_date(s).is_none() {
        bail!("{label}: expected YYYY-MM-DD, got {s}");
    }
    if s < cat.horizons_first.as_str() || s > cat.horizons_last.as_str() {
        bail!("{label}: Horizons only generates small-body files from {} to {}", cat.horizons_first, cat.horizons_last);
    }
    Ok(s.to_string())
}

fn install(cat: &Catalog, a: InstallArgs) -> Result<()> {
    let dir = a.dir.resolve()?;
    let unknown_planets = |key: &str| anyhow!("unknown planetary file {key}; see `astro-calc-build tiers`");
    let tier_key = match (&a.tier, &a.planets) {
        (Some(k), _) => k.clone(),
        // A planetary file alone picks the first tier that uses it, else tier 2's settings.
        (None, Some(key)) => {
            let p = cat.planet(key).ok_or_else(|| unknown_planets(key))?;
            cat.tiers.iter().find(|t| t.planets == p.key).map_or("2".into(), |t| t.key.clone())
        }
        (None, None) => install::choose_tier(cat)?,
    };
    let t = cat.tier(&tier_key).ok_or_else(|| {
        anyhow!("no tier '{tier_key}'; choose from {}", cat.tiers.iter().map(|t| t.key.as_str()).collect::<Vec<_>>().join(", "))
    })?;
    let chosen = match &a.planets {
        Some(key) => cat.planet(key).ok_or_else(|| unknown_planets(key))?,
        None => cat.planet(&t.planets).expect("validated tier"),
    };
    // --planets-from / --planets-to: an excerpt of the full file over those years.
    let planets = match (a.planets_from, a.planets_to) {
        (None, None) => chosen.clone(),
        (from, to) => {
            let (base, years) = match &chosen.excerpt {
                Some(e) => (cat.planet(&e.of).expect("validated excerpt"), Some((e.from, e.to))),
                None => (chosen, None),
            };
            let need = || anyhow!("give both --planets-from and --planets-to to cut {}", base.key);
            let from = from.or(years.map(|y| y.0)).ok_or_else(need)?;
            let to = to.or(years.map(|y| y.1)).ok_or_else(need)?;
            if to <= from {
                bail!("--planets-to must come after --planets-from");
            }
            let cut = catalog::excerpt_file(base, from, to);
            if let (Some(full), Some(want)) = (base.range, cut.range)
                && want.intersect(&full).is_none()
            {
                bail!("{} covers {full}, which does not reach {} to {}", base.file, catalog::year_label(from), catalog::year_label(to));
            }
            cut
        }
    };
    let from_file = a.from_file.map(absolute).transpose()?;
    if let Some(f) = &from_file {
        if planets.excerpt.is_none() {
            bail!("--from-file is for excerpts: --tier ancient, or --planets-from / --planets-to");
        }
        if !f.is_file() {
            bail!("--from-file: {} not found", f.display());
        }
    }
    let bodies = cat.select_bodies(a.bodies.as_deref().unwrap_or(&t.bodies)).map_err(|e| anyhow!("--bodies: {e}"))?;
    let from = a.bodies_from.as_deref().map(|d| check_horizons_date(cat, "--bodies-from", d)).transpose()?;
    let to = a.bodies_to.as_deref().map(|d| check_horizons_date(cat, "--bodies-to", d)).transpose()?;
    let stations = if a.no_stations {
        None
    } else {
        Some(match (t.stations, a.stations_from, a.stations_to) {
            (span, None, None) => span,
            (StationSpan::Years(f, to), from, end) => StationSpan::Years(from.unwrap_or(f), end.unwrap_or(to)),
            (StationSpan::Full, from, end) => StationSpan::Years(from.unwrap_or(-100_000), end.unwrap_or(100_000)),
        })
    };
    if let Some(StationSpan::Years(f, to)) = stations
        && to <= f
    {
        bail!("--stations-to must come after --stations-from");
    }
    let opts = install::Options {
        cat,
        dir,
        planets,
        from_file,
        bodies,
        tier_span: (t.bodies_from.clone(), t.bodies_to.clone()),
        span_override: (from, to),
        stations,
        station_bodies: a.station_bodies,
        keep_other_planets: a.keep_other_planets,
        force: a.force,
        yes: a.yes,
        dry_run: a.dry_run,
    };
    for b in &opts.bodies {
        if let Some((f, t)) = opts.span(b)
            && t <= f
        {
            bail!("{}: the span ends before it starts ({f} to {t})", b.key);
        }
    }
    install::run(opts)
}

fn main() -> Result<()> {
    // Exit quietly when output is piped into something that stops reading (`| head`),
    // as other command-line tools do, instead of panicking on the broken pipe.
    #[cfg(unix)]
    // SAFETY: called once at start-up, before any other thread exists.
    unsafe {
        libc::signal(libc::SIGPIPE, libc::SIG_DFL);
    }
    let cli = Cli::parse();
    let (cat, origin) = config::load(cli.config.as_deref(), cli.no_config)?;
    let Some(command) = cli.command else {
        print_overview(&cat, origin.as_deref());
        return Ok(());
    };
    match command {
        Command::Install(a) => install(&cat, a),
        Command::Tiers => {
            print_catalogue(&cat, origin.as_deref());
            Ok(())
        }
        Command::Config => {
            print!("{}", config::to_toml(&cat, origin.as_deref())?);
            Ok(())
        }
        Command::Info { dir } => info::summary(&cat, &dir.resolve()?),
        Command::Stations(a) => {
            let dir = a.dir.resolve()?;
            let out = a.out.unwrap_or_else(|| dir.join("stations.bin"));
            let years = match (a.from, a.to) {
                (None, None) => None,
                (f, t) => Some((f.unwrap_or(-100_000), t.unwrap_or(100_000))),
            };
            install::build_stations(&dir, &out, &a.bodies, years)?;
            info::summary(&cat, &dir)
        }
    }
}

/// What the tool does, the tiers, the commands, and the state of `$ASTRO_DATA_DIR` if set.
fn print_overview(cat: &Catalog, origin: Option<&Path>) {
    println!("astro-calc-build {}: {ABOUT}\n", env!("CARGO_PKG_VERSION"));
    println!("astro-calc-core ships without JPL data. This tool downloads a planetary file from JPL,");
    println!("has JPL Horizons generate the small-body files, checks every file and builds the");
    println!("station table, in any data folder you name with --data-dir (-d).\n");
    println!("Tiers (`astro-calc-build tiers` for details)\n");
    let w = cat.tiers.iter().map(|t| t.key.len()).max().unwrap_or(1);
    for t in &cat.tiers {
        let Some(p) = cat.planet(&t.planets) else { continue };
        let bodies = cat.select_bodies(&t.bodies).unwrap_or_default();
        let est: f64 = bodies.iter().filter_map(|b| body_estimate(b, &cat.body_span(b, &t.bodies_from, &t.bodies_to))).sum();
        println!(
            "  {:<w$}  {:<11} {:<30} {:<20} + {} bodies  {}",
            t.key,
            t.title,
            p.range.map_or("dates read on download".into(), |r| r.to_string()),
            p.file,
            bodies.len(),
            cat.planet_estimate(p).map_or("size unknown".into(), |s| format!("~{}", fmt_bytes(s + est)))
        );
    }
    println!("\nCommands\n");
    for (cmd, what) in [
        ("install -d PATH [--tier T]", "Download, verify and set up a tier (--dry-run shows the plan)"),
        ("info -d PATH", "What a data folder holds and the dates it covers"),
        ("stations -d PATH", "Rebuild the station table from the installed files"),
        ("tiers", "Tiers, planetary files, small bodies and body sets in detail"),
        ("config", "The effective configuration as TOML, to start a config file"),
    ] {
        println!("  {cmd:<28} {what}");
    }
    println!("\n  `astro-calc-build <command> --help` lists every option; `--help` gives the full description.\n");
    match origin {
        Some(p) => println!("Config file:  {}", p.display()),
        None => println!("Config file:  none (built-in catalogue; see `astro-calc-build --help`)"),
    }
    match std::env::var_os("ASTRO_DATA_DIR").filter(|v| !v.is_empty()) {
        Some(d) => {
            let d = PathBuf::from(d);
            let state = match astro_calc_core::Ephemeris::open_dir(&d) {
                Ok(e) => e.date_range().map_or("no planetary file".into(), |r| format!("Sun, Moon and planets {r}")),
                Err(_) => "no usable .bsp files yet".into(),
            };
            println!("Data folder:  {} (ASTRO_DATA_DIR): {state}", d.display());
        }
        None => println!("Data folder:  not set; pass --data-dir <PATH> (-d), or set ASTRO_DATA_DIR"),
    }
}
