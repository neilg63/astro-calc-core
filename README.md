# astro-calc-core

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

## Data files

The library reads JPL ephemeris files in the SPK (`.bsp`) format. They are not included in
this repository; fetch them into `data/`:

```bash
curl -o data/de440s.bsp https://ssd.jpl.nasa.gov/ftp/eph/planets/bsp/de440s.bsp
```

`de440s.bsp` (33 MB) covers the Sun, Moon and planets from 1850 to 2150 and is all that is
needed for a standard chart. The larger `de440.bsp` (1550 to 2650) uses the same format.

The minor bodies come from JPL Horizons, one file each:

```bash
python3 -I tools/fetch_horizons_spk.py data
```

The date range the library supports is the range of the files loaded.

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

- Parallel bulk calculation and async wrappers.

## Examples

```bash
cargo run --release --example positions
```

`positions` is a quick start. `timing` and `stations_timing` report speeds, and
`build_station_table` generates the station table file. The other examples and the Python
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
