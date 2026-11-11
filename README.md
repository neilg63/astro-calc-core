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
- Mean and true lunar nodes.
- Geocentric or topocentric; ecliptic, equatorial or horizontal coordinates.
- Longitude, latitude, distance and their daily speeds.
- Bulk methods: many bodies at one time, or many bodies over a series of times.

**Angles and houses**

- Ascendant, MC, ARMC, vertex, equatorial ascendant, co-ascendants and polar ascendant.
- House systems: Whole Sign, Equal, Vehlow, Porphyry, Placidus, Koch, Alcabitius, Campanus,
  Regiomontanus, Topocentric (Polich/Page), Morinus, Meridian and Horizontal.

**Ayanamshas**

- Seventeen sidereal reference systems selected by an enum, including Lahiri, Krishnamurti,
  Raman, Yukteshwar and several star-based definitions (Citra, Revati, Pushya, Mula,
  Galactic Centre).
- Calculated separately and applied as offsets, so any position can be made sidereal.

**Rise, set and transit**

- Rise, set, upper and lower meridian transit for any supported body.
- Options for disc centre or limb, refraction, and a fixed disc size.

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
| Hipparcos, the New Reduction (van Leeuwen 2007), via NASA HEASARC | Star positions for star-based ayanamshas |
| Reid & Brunthaler (2004) | Position of the Galactic Centre (Sgr A*) |
| Hilton (2005); Harris (1961); Meeus, *Astronomical Algorithms* ch. 45 | Planetary magnitudes and Saturn's rings |
| *Explanatory Supplement to the Astronomical Almanac*, 3rd ed. | The apparent-position pipeline |
| Holden, *The Elements of House Division*; Meeus, *Astronomical Algorithms* | House systems and general reference |

Licence terms for the transcribed tables and formulae are in [NOTICE](NOTICE).

## Not yet implemented

- Lunar apogee and perigee points.
- Hypothetical bodies and fixed stars as positions.
- The Krusinski house system.
- Magnitudes for asteroids.
- Parallel bulk calculation and async wrappers.
- Viewpoints other than the Earth (heliocentric or from another planet).

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
