# astro-calc-build

Installs the JPL ephemeris files that [astro-calc-core](https://crates.io/crates/astro-calc-core) reads.

astro-calc-core is distributed without JPL data. Without it, calculations use the built-in analytic
series (with the `analytic` feature) or fail. This tool, working in the data folder you name:

1. downloads a planetary file from JPL (resumable),
2. has JPL Horizons generate files for the small bodies (or downloads ready-made ones),
3. checks every file (an SPK index holding the expected bodies, and SHA-256 when configured)
   before moving it into place,
4. removes any other planetary file in the folder, so ephemerides are never mixed,
5. builds `stations.bin`, the pre-computed table of planetary stations.

```sh
cargo install astro-calc-build      # or: cargo install --path crates/astro-calc-build
astro-calc-build                                     # overview
astro-calc-build install -d /srv/astro/data         # choose a tier interactively
```

`install`, `info` and `stations` take the data folder with `-d` / `--data-dir`: any path,
relative, absolute or `~/…`. It can also be given as a plain argument, and `ASTRO_DATA_DIR`
supplies it when omitted. Point the apps that read the data at the same folder: for
astro-calc-api, set `ASTRO_DATA_DIR` to it.

## Tiers

| Tier | Sun, Moon, planets | File | Small bodies | Download |
|---|---|---|---|---|
| 1 `modern` | 1849-12-26 – 2150-01-22 | `de440s.bsp` 31 MiB | 1850 – 2150 | ~62 MiB |
| 2 `historical` | 1549-12-31 – 2650-01-25 | `de440.bsp` 114 MiB | 1600 – 2500 | ~207 MiB |
| 3 `ancient` | 1000 BCE – 3000 CE | `de441_-999_3000.bsp` 415 MiB | 1600 – 2500 | ~508 MiB |
| 4 `long-range` | 13201 BCE – 17191 CE | `de441.bsp` 3.08 GiB | 1600 – 2500 | ~3.2 GiB |

`--tier` takes the number or the name. JPL Horizons only generates small-body files between
1600-01-01 and 2500-01-01, so tiers 2 to 4 cover the small bodies over a shorter span than the
planets.

### Excerpts (tier 3)

JPL publishes nothing between `de440.bsp` (1550–2650) and the full `de441.bsp`, so tier 3 is cut
from DE441. Its segments hold fixed-length records, so the installer reads DE441's index with a
few small range requests, fetches only the records covering the span, and writes them with a
new index as `de441_-999_3000.bsp`. The records are copied byte for byte: positions are exactly
DE441's. Every record's own midpoint is checked against its place in the new file, so a single
misplaced word would be caught. An interrupted run resumes.

Any span of any planetary file can be cut the same way:

```sh
astro-calc-build install -d ~/ephem --planets de441 --planets-from -499 --planets-to 2500   # 500 BCE – 2500 CE
astro-calc-build install -d ~/ephem --tier ancient --from-file ~/Downloads/de441.bsp       # from a local copy
```

`--from-file` cuts from a full file you already have instead of using range requests. In
the config file, `[planets.mine] excerpt = { of = "de441", from = -499, to = 2500 }` adds one
to the catalogue.

In tier 4, positions thousands of years from today are limited by the ΔT and precession models,
not by the ephemeris. Its station table therefore defaults to 3000 BCE – 3000 CE
(`--stations-from` / `--stations-to` to change). Tier 3 already ends about where those models
stop being sound.

`--planets de442` installs DE442 (2026) in place of DE440. The library's test fixtures assume DE440.

## Commands

```sh
astro-calc-build install -d /srv/astro/data --tier 2 --dry-run    # show the plan only
astro-calc-build install --data-dir ../astro-calc-core/data --tier 1 --yes
astro-calc-build install -d ~/ephem --tier ancient                # 1000 BCE – 3000 CE, ~508 MiB
astro-calc-build install -d ~/ephem --tier 4 --bodies none --stations-from -999 --stations-to 2500
astro-calc-build install -d ~/ephem --tier 2 --bodies main,chiron
astro-calc-build info -d /srv/astro/data    # what is installed and the dates it covers
astro-calc-build stations -d /srv/astro/data --from 1000 --to 2200
astro-calc-build tiers                      # tiers, files, bodies and sets in detail
astro-calc-build config                     # the effective configuration, as TOML
```

Years are astronomical (0 is 1 BC, -2999 is 3000 BC), as in ISO 8601.

Re-running `install` skips files that are already in place, resumes an interrupted download, and
retries any small body that failed.

## Config file

A TOML file overrides the built-in URLs and file names and adds planetary files, small bodies,
body sets and tiers. The first found is used:

1. `--config <file>` or `$ASTRO_CALC_BUILD_CONFIG` (must exist)
2. `./astro-calc-build.toml`
3. `$XDG_CONFIG_HOME/astro-calc-build/config.toml` (default `~/.config/...`)

`--no-config` ignores them. Entries merge with the built-ins by key, so a file holds only what
it changes; `astro-calc-build config > astro-calc-build.toml` writes the complete effective
configuration as a starting point. Unknown keys and references to things that do not exist
are reported with their names before anything is downloaded.

```toml
[sources]
jpl_bsp_url    = "https://mirror.example.org/jpl/"   # base URL for planetary files
horizons_api   = "https://ssd.jpl.nasa.gov/api/horizons.api"
horizons_first = "1600-01-01"                        # Horizons' limits, in case JPL changes them
horizons_last  = "2500-01-01"

# Override a planetary file. A new url or file drops the built-in size and dates; the plan
# then reads them from the server (size, and dates from the file's index by range requests).
[planets.de441]
url = "https://mirror.example.org/big/de441.bsp"
sha256 = "…64 hex digits…"                           # optional check after download

# Add one (a future release, say)
[planets.de443]
file = "de443.bsp"

# Small bodies: change a built-in one's span, or add one by Minor Planet Center number
[bodies.chiron]
from = "1700-01-01"
to   = "2300-01-01"

[bodies.eris]
mpc = 136199                 # Horizons COMMAND "136199;" and NAIF id 20136199
name = "Eris"                # astro-calc-core and the API know it as mpc:136199

# …or point at a ready-made SPK file
[bodies.sedna]
url = "https://example.org/spk/sedna.bsp"
naif_id = 20090377

# Sets, for --bodies and tiers (built in: default, main, centaurs; and all / none)
[body_sets]
tno = ["eris", "sedna"]

# Tiers: change 1-3 or add named ones (--tier medieval)
[tiers.medieval]
title = "Medieval"
planets = "de441"
bodies = ["main", "tno"]     # or "main,tno", "all", "none"
bodies_from = "1600-01-01"
bodies_to = "2100-01-01"
stations = { from = 500, to = 1700 }   # or "full"
```

Every field is listed in the output of `astro-calc-build config`. Command-line options override
the config file, which overrides the built-ins.
