Not committed (see .gitignore). Install with the data CLI (crates/astro-calc-build):

    cargo install --path crates/astro-calc-build
    astro-calc-build install -d data --tier 1    # de440s + small bodies, 1850–2150
    astro-calc-build info -d data

Tiers 2 (1549–2650, de440), 3 (1000 BCE–3000 CE, cut from de441) and 4 (13201 BCE–17191 CE,
de441) cover longer spans; see `astro-calc-build tiers`.

The test reference file is fetched separately:

    curl -O https://ssd.jpl.nasa.gov/ftp/eph/planets/ascii/de440/testpo.440

`Ephemeris::open_dir("data")` loads every `.bsp` file here.
