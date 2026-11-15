# astro-calc-core

[![mirror](https://img.shields.io/badge/mirror-github-blue)](https://github.com/neilg63/astro-calc-core)
[![crates.io](https://img.shields.io/crates/v/astro-calc-core.svg)](https://crates.io/crates/astro-calc-core)
[![docs.rs](https://docs.rs/astro-calc-core/badge.svg)](https://docs.rs/astro-calc-core)

A pure-Rust astronomical calculation library for astrology software. It computes
positions, houses, ayanamshas, rise and set times, stations and aspects directly from
NASA JPL ephemeris files, with no C dependencies.

This crate, `astro-calc-core`, is the library. Applications are planned on top of it as
separate projects: `astro-calc-api` (a REST API) and `astro-calc-ui` (a desktop GUI).

## Status: pre-alpha

This is a pre-alpha release: an early version published for review and experiment.

- The public API is still changing and may break between versions.
- Results have been checked against JPL's own test data, published almanac values and an
  independent reference implementation, but the library has not yet been used in production.
- Some features are incomplete; see [Not yet implemented](#not-yet-implemented).

Do not rely on it where an error would matter without checking the results yourself.

## Core features

All times are Julian Days in Universal Time (`f64`, wrapped as `JdUt`). ISO-8601-like
strings such as `1990-08-23T12` are accepted by wrapper methods and converted once.

**Positions**

- Sun, Moon, Mercury to Pluto and the Earth, plus Ceres, Pallas, Juno, Vesta, Chiron and
  Pholus when their data files are present.
- Mean and true lunar nodes. Ketu is derived from Rahu when both are asked for.
- The ascendant as a body (`as`), so it can be listed with the planets.
- Body sets: `default` (ascendant, Sun to Pluto, Ketu, Rahu) and `basic` (the same without
  Uranus, Neptune and Pluto), adapted automatically to the perspective.
- Longitude-only methods, tropical or sidereal, about twice as fast as full positions.
- Lunar apogee: mean and osculating (Black Moon Lilith), plus the interpolated apogee and
  perigee.
- Fixed stars: the 410 stars with an IAU-approved name and a Hipparcos number, addressed
  by name or as `hip:65474`, with proper motion and parallax. Other stars can be supplied
  from catalogue data.
- Hypothetical bodies from Keplerian elements: the eight Uranian points (Cupido to
  Poseidon) and Isis-Transpluto. `OrbitalElements` is public for defining others.
- Ecliptic, equatorial or horizontal coordinates, with fields named after the system:
  `lng`/`lat` (`BodyPos`), `ras`/`dcl` (`BodyPosEqu`) or `azi`/`alt` (`BodyPosHor`), and
  `BodyPosExt` for all three from a single calculation. A shared `BodyPosition` trait
  reads any of them through `lng()`, `lat()`, their speeds and `dist()`.
- Eleven perspectives: geocentric, topocentric, heliocentric, selenocentric (Moon),
  hermeocentric (Mercury), cytherocentric (Venus), areocentric (Mars), jovicentric
  (Jupiter), cronocentric (Saturn), barycentric, and Local Space (azimuth and altitude on
  the observer's horizon). From other centres the Earth becomes a body, and coordinates
  stay on the Earth's ecliptic or equator of date.
- Longitude, latitude, distance and their daily speeds.
- Bulk methods: many bodies at one time, or many bodies over a series of times.

**Angles and houses**

- Ascendant, MC, ARMC, vertex, equatorial ascendant, co-ascendants and polar ascendant.
- House systems: Whole Sign, Equal, Vehlow, Porphyry, Placidus, Koch, Alcabitius, Campanus,
  Regiomontanus, Topocentric (Polich/Page), Morinus, Meridian, Horizontal and
  Krusinski-Pisa-Goelzer.

**Ayanamshas**

- Seventeen sidereal reference systems selected by an enum, including Lahiri, Krishnamurti,
  Raman, Yukteshwar and several star-based definitions (Citra, Revati, Pushya, Mula,
  Galactic Centre).
- Calculated separately and applied as offsets, so any position can be made sidereal.

**Rise, set and transit**

- Rise, set, upper and lower meridian transit for any supported body.
- Disc alignment (upper limb, centre or lower limb), with the actual or a fixed disc size.
- Refraction on, off, or scaled for a given air pressure and temperature.
- Civil, nautical and astronomical twilight.
- A raised or lowered horizon, and the dip of the horizon for an elevated observer.
- Geocentric positions and bodies projected onto the ecliptic, including the traditional
  Hindu rising.

**Planetary stations**

- Retrograde start and end, and the fastest direct and retrograde motion of each cycle.
- Current speed and retrograde state.
- An optional pre-computed `StationTable`: a compact binary file (about 175 KB for eight
  planets over 300 years) that answers lookups in nanoseconds and can be rebuilt from any
  ephemeris file.

**Moon phases**

- New, first quarter, full and last quarter, geocentric or topocentric.

**Indian solar day**

- Sunrise-to-sunrise days with muhurta, ghati, vighati and lipta, and positions at sunrise
  or sunset.

**Aspects**

- Aspect angles and orbs are supplied by the caller; none are built in.
- Instant matching between two longitudes or two sets of positions.
- A period search between two people's natal or secondary-progressed charts, returning when
  each aspect entered its orb, became exact and left it.

**Phenomena**

- Phase angle, illuminated fraction, elongation, apparent diameter and magnitude.
- Minor bodies use the IAU H–G magnitude system; `hg_magnitude` accepts catalogue values
  for bodies that are not built in.

**Design**

- No global state: an `Ephemeris` is `Send + Sync` and can be shared between threads.
- Synchronous and CPU-bound, so it can be called from async code without a runtime dependency.
- Data files are memory-mapped by default. With the `mmap` feature off, data can be loaded
  from bytes, which allows a WebAssembly build.
- Optional `serde` support for the public types.

## Example

```rust
use astro_calc_core::*;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let eph = Ephemeris::open_dir("data")?;
    let jd = JdUt::from_iso("1996-01-12T11:25:00")?;
    let geo = GeoPos::new(56.1, -3.6, 0.0)?;

    let bodies = parse_body_list("su,mo,me,ve,ma,ju,sa")?;
    for p in eph.calc_positions(&bodies, jd, geo, CoordSystem::Ecliptic, Perspective::Topocentric)? {
        println!("{} {:.4} ({:+.4} deg/day)", p.body, p.lng, p.lng_speed);
    }

    let lahiri = eph.ayanamsha(jd, Ayanamsha::Lahiri)?;
    println!("Lahiri ayanamsha {lahiri:.4}");

    let mercury = eph.next_stations(Body::Mercury, jd, 4, 0)?;
    println!("Mercury retrograde: {}", mercury.is_retrograde());
    Ok(())
}
```

## Bulk calculation

For large sets of records (births, events) the `bulk` module calculates in one call. Every
record carries a key of the caller's choosing, which is handed back untouched, so records
can be filtered beforehand and results matched to them afterwards.

```rust
let records = vec![
    ChartInput::new("UniqueNameA", 2_443_121.236, GeoPos::new(44.8993, 12.39373, 0.0)?),
    ChartInput::new("UniqueNameB", 2_440_408.892, GeoPos::new(48.1462, 11.38975, 0.0)?),
];

// A chart for each record: the default bodies, ascendant included.
let found = eph.bulk_charts(&records, &ChartOptions::default());

// Records with a body quincunx (150° ± 2°1') the first record's Mars.
let criteria = AspectCriteria::new(&[Body::Mars], &Body::LUMINARIES_AND_PLANETS, &aspect_orbs(&[(150.0, 2.0 + 1.0 / 60.0)]));
let found = eph.aspects_one_to_many(&records[0], &records[1..], &criteria)?;
for m in &found.matches {
    println!("{}: {:?} is {:+.3}° from exact", m.b, m.hits[0].body_b, m.hits[0].deviation);
}
```

| Function | Compares |
|---|---|
| `bulk_charts` | nothing: positions of the chosen bodies for every record |
| `aspects_in_pairs` | the two records of each given pair |
| `aspects_one_to_many` | one record with each record of a list |
| `aspects_many_to_many` | every two records of one list |

The aspect searches return only the pairs that match, each with the bodies involved and the
deviation from exact. `AspectCriteria::progressed_to` compares secondary-progressed charts
(a day for a year) as they stand on a given date. Records need not be sorted: they are worked
through in date order internally and returned in input order. A record that cannot be
calculated is listed under `errors` with its key and does not fail the batch.

With the `parallel` feature the work is spread over all cores with `rayon`. On an 18-core
Apple M5 Pro, DE440s, release build:

| Task | One thread | `parallel` |
|---|---|---|
| 10,000 charts of 13 bodies with speeds | 145 ms | 16–20 ms |
| One record against 9,999 (1 × 10 bodies) | 54 ms | 5 ms |
| 5,000 pairs (10 × 10 bodies) | 65 ms | 8–19 ms |
| 2,000 records, every pair both ways (4 million) | 94 ms | 13 ms |

## Data files

The library reads JPL ephemeris files in the SPK (`.bsp`) format. They are not included in
the crate or this repository. Without them, calculations use the analytical series (with the
`analytic` feature, below) or fail with an error naming the dates covered.

### Installing the data: `astro-calc-build`

`astro-calc-build` is a companion command-line tool, a separate crate in this repository
(`crates/astro-calc-build`). It downloads the files into any folder you choose, checks each one,
and builds the station table:

```bash
cargo install astro-calc-build          # from this repository: cargo install --path crates/astro-calc-build
astro-calc-build                        # overview: tiers, commands, the state of ASTRO_DATA_DIR
astro-calc-build install -d /srv/astro/data --tier 1
```

Choose a tier by the dates you need:

| Tier | Sun, Moon and planets | Planetary file | Small bodies | Download |
|---|---|---|---|---|
| 1 `modern` | 1849-12-26 to 2150-01-22 | `de440s.bsp`, 31 MiB | 1850 to 2150 | ~62 MiB |
| 2 `historical` | 1549-12-31 to 2650-01-25 | `de440.bsp`, 114 MiB | 1600 to 2500 | ~207 MiB |
| 3 `ancient` | 1000 BCE to 3000 CE | `de441_-999_3000.bsp`, 415 MiB, cut from DE441 | 1600 to 2500 | ~508 MiB |
| 4 `long-range` | 13201 BCE to 17191 CE | `de441.bsp`, 3.08 GiB | 1600 to 2500 | ~3.2 GiB |

Tiers can be named by number or by name (`--tier ancient`). The planetary file comes from
JPL's server. JPL Horizons generates the small-body files (Ceres, Pallas, Juno, Vesta, Chiron,
Pholus) on request, but only between 1600 and 2500.

Tier 3 is an excerpt of DE441: the installer reads DE441's index, fetches only the records
covering 1000 BCE to 3000 CE with HTTP range requests, and writes them as a smaller SPK file
with a new index. The records are copied byte for byte, so positions are exactly DE441's. It
is about as far as the ΔT and precession models stay sound. In tier 4, dates thousands of
years from today are limited by those models rather than by the ephemeris, so its station
table defaults to 3000 BCE to 3000 CE. Any span can be cut the same way with
`--planets-from` / `--planets-to`, from JPL or from a local copy (`--from-file`).

**Commands**

| Command | What it does |
|---|---|
| `install -d PATH` | Shows a plan (files, dates, sizes, anything to be removed), asks, then downloads, verifies and builds `stations.bin`. Without `--tier` on a terminal it asks for one |
| `info -d PATH` | The files in a folder, the dates each covers, the span for the Sun, Moon and planets, every small body with its body key, and whether `stations.bin` matches |
| `stations -d PATH` | Rebuilds `stations.bin`; `--from` / `--to` years, `--bodies` |
| `tiers` | Tiers, planetary files, small bodies and body sets in detail |
| `config` | The effective configuration as TOML |

`-d` / `--data-dir` takes any path (relative, absolute or `~/…`); the folder can also be
given as a plain argument, and `ASTRO_DATA_DIR` supplies it when omitted. Point the apps that
read the data at the same folder: `Ephemeris::open_dir(path)` here, `ASTRO_DATA_DIR` for
astro-calc-api.

**`install` options**

| Option | Effect |
|---|---|
| `-t, --tier T` | `1`-`4`, a tier name (`modern`, `historical`, `ancient`, `long-range`) or one in the config file |
| `--planets KEY` | Another planetary file: `de440s`, `de440`, `de441`, `de442` (2026, DE440's successor), `de441-ancient` |
| `--planets-from YEAR`, `--planets-to YEAR` | Cut the planetary file to these years (astronomical: -999 is 1000 BCE), fetching only those records |
| `--from-file PATH` | Cut the excerpt from a local copy of the full file instead of using range requests |
| `--bodies LIST` | `all`, `none`, set names (`default`, `main`, `centaurs`) and body keys, comma-separated |
| `--bodies-from DATE`, `--bodies-to DATE` | Span of the Horizons small-body files, `YYYY-MM-DD` |
| `--stations-from YEAR`, `--stations-to YEAR` | Span of the station table, astronomical years (-2999 is 3000 BC) |
| `--station-bodies LIST` | Planets in the station table (default `me,ve,ma,ju,sa,ur,ne,pl`) |
| `--no-stations` | Skip the station table |
| `--keep-other-planets` | Keep other planetary files in the folder (removed by default, so ephemerides are not mixed) |
| `--dry-run` | Show the plan only |
| `--force` | Download and rebuild even what is already installed |
| `-y, --yes` | Do not ask for confirmation |

Each file is written as `<name>.part`, checked (an SPK index holding the expected bodies, and
SHA-256 when one is configured) and only then renamed, so an interrupted run never leaves a
broken file for `open_dir`. Re-running skips what is in place and resumes a partial download.

**Config file.** `astro-calc-build.toml` in the current folder or
`~/.config/astro-calc-build/config.toml` (or `--config FILE`) overrides the download URLs and
file names (a mirror, say), and adds planetary files, small bodies, body sets and named tiers.
A small body can be added by its Minor Planet Center number (`mpc = 136199` for Eris), which
this library then calculates as `mpc:136199`. `astro-calc-build config` prints the complete
effective configuration to start from; see `crates/astro-calc-build/README.md` for the format.

### What the library does with them

The date range the library supports is the range of the files loaded.
`scan_dir(dir)` lists the `.bsp` files in a directory with the dates and bodies each holds,
without loading them for calculation; `Ephemeris::files()` gives the same for a loaded
ephemeris. `Ephemeris::date_range()` reports the span shared by the Sun, Moon and planets,
and `body_range(body)` the span for one body (the minor-body files are often shorter).

Every calculation that uses the ephemeris reports the engine behind it, as a `CalcEngine`:
in an `engine` field on result structs, or in a `Calc<T>` wrapper (`engine` and `value`,
dereferencing to the value) where the result is a bare number or list. A date outside the
data fails with an error naming the body and the dates covered; `CalcError::is_out_of_range()`
picks those out. The engine is `CalcEngine::Jpl` unless the `analytic` feature is on (below). Functions that are pure
mathematics (houses, angles, aspects between given positions, formula-based ayanamshas) do
not report one.

## Analytical fallback (feature `analytic`)

With the optional `analytic` feature, a date outside the JPL files is calculated from
built-in series instead of failing, and the result reports `CalcEngine::Analytic`.
`Ephemeris::analytic()` builds an ephemeris with no files at all, which suits WebAssembly.

```toml
astro-calc-core = { version = "0.1", features = ["analytic"] }
```

- **Covered:** the Sun and Mercury to Neptune from the VSOP87 theory, and the Moon from the
  lunar theory ELP/MPP02, for 4000 years either side of J2000. Everything built on them
  works: the lunar nodes and apsides, all perspectives, rise and set, stations, moon
  phases, phenomena, stars, hypothetical bodies, houses and ayanamshas.
- **Not covered:** Pluto and the minor bodies. Outside the files these fail with an error
  naming the dates the files hold.
- **Accuracy inside 1850 to 2150**, against the JPL files, in apparent geocentric longitude:
  0.1″ for the Sun and Moon, 0.2″ for Mercury and Venus, 0.7″ for Mars and Jupiter, 0.4″
  for Saturn, 1.7″ for Uranus and 2.5″ for Neptune.
- **Accuracy before 1850**, from four spot checks against the AstroCalc reference API at
  AD 500, 1200, 1585 and 1700, after removing the difference in ΔT: Mercury and Venus within
  0.4″, Mars and Jupiter within 1″, Saturn within 3″, the Moon 1.5″ at 1700 growing to 6″ at
  AD 500. Uranus and Neptune drift much more, 3″ to 8″ around 1600 and about 40″ at AD 500,
  because the theory was fitted to an older JPL ephemeris in which their orbits were less
  well known. Nothing earlier than AD 500 has been checked.
- **ΔT matters more than the series at remote dates.** Published reconstructions of the
  Earth's rotation differ by tens of seconds before 1700 and by a minute or two in the
  Middle Ages, and each second moves the Moon by half an arcsecond.
- **Speed and size:** about two hundred times slower than the JPL files (about 2 ms for ten
  bodies) and 0.9 MB larger. The embedded tables are truncations of the published series;
  `tools/build_analytic_tables.py` rebuilds them at any level.
- **One engine per moment:** a date outside the span the files share for the Sun, Moon and
  planets is served entirely by the series, never pieced together from both.

## Data sources and references

| Source | Used for |
|---|---|
| [JPL Planetary and Lunar Ephemerides, DE440](https://ssd.jpl.nasa.gov/planets/eph_export.html) (Park et al. 2021) | Positions of the Sun, Moon and planets |
| [JPL Horizons](https://ssd.jpl.nasa.gov/horizons/) | Ephemeris files for asteroids and centaurs |
| [NAIF SPK and DAF Required Reading](https://naif.jpl.nasa.gov/naif/documentation.html) | The ephemeris file format |
| [ERFA](https://github.com/liberfa/erfa), derived from IAU SOFA (BSD 3-clause) | IAU 2006 precession, IAU 2000B nutation, Earth rotation, light deflection and aberration |
| IERS Conventions (2010) | Earth orientation models |
| US Naval Observatory ΔT tables | Difference between Terrestrial Time and Universal Time |
| Morrison & Stephenson (2004) | Long-term ΔT outside the observed range |
| Bretagnon & Francou (1988), VSOP87, as distributed by IMCCE | Analytical series for the Sun and planets (feature `analytic`) |
| Chapront & Francou (2003), ELP/MPP02, as distributed by Paris Observatory | Analytical series for the Moon (feature `analytic`) |
| Espenak & Meeus (2006), Five Millennium Canon of Solar Eclipses | ΔT before 1657 |
| Hipparcos, the New Reduction (van Leeuwen 2007), via NASA HEASARC | Star positions, proper motions and parallaxes |
| [IAU Catalog of Star Names](https://www.iau.org/public/themes/naming_stars/) (IAU Working Group on Star Names) | Star names and visual magnitudes |
| Reid & Brunthaler (2004) | Position of the Galactic Centre (Sgr A*) |
| Hilton (2005); Harris (1961); Meeus, *Astronomical Algorithms* ch. 45 | Planetary magnitudes and Saturn's rings |
| [JPL Small-Body Database](https://ssd.jpl.nasa.gov/tools/sbdb_lookup.html); Bowell et al. (1989) | Asteroid H and G values, radii and the H–G magnitude law |
| Chapront et al., as given in Meeus ch. 47 and 50 | Mean lunar node, perigee and anomaly |
| Neely (1988) | Orbital elements of the Uranian points |
| Krusinski (1995) | The Krusinski-Pisa-Goelzer house system |
| *Explanatory Supplement to the Astronomical Almanac*, 3rd ed. | The apparent-position pipeline |
| Holden, *The Elements of House Division*; Meeus, *Astronomical Algorithms* | House systems and general reference |

Licence terms for the transcribed tables and formulae are in [NOTICE](NOTICE).

## Known limits of the newer features

- The interpolated lunar apogee and perigee follow this library's own construction (a curve
  through the Moon's actual apsides). Other publishers construct theirs differently, so
  values between apsides will not agree exactly.
- Of the hypothetical bodies, Kronos and Isis have been checked against an independent
  reference. The other seven use the same method with published elements but have not been
  independently checked.
- The Krusinski house system has been tested for its defining properties only.

## Not yet implemented

- Async wrappers for the bulk functions.
- Pluto in the analytical fallback (see *Analytical fallback*).

## Examples

```bash
cargo run --release --example positions
```

`positions` is a quick start. `bulk` shows the bulk functions with timings, `timing` and `stations_timing` report speeds, and
`build_station_table` generates the station table file (as does `astro-calc-build stations`). The other examples and the Python
scripts in `tools/` are maintenance aids: they regenerate test fixtures and data, and are
not part of the published crate.

## Tests

```bash
cargo test --release
```

Tests that need ephemeris files are skipped when the files are absent.

## Licence

Licensed under either of

- the MIT licence ([LICENSE-MIT](LICENSE-MIT)), or
- the Apache License, Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE)),

at your option. Third-party attributions are in [NOTICE](NOTICE).

Unless you explicitly state otherwise, any contribution intentionally submitted for
inclusion in this work shall be dual licensed as above, without any additional terms or
conditions.
