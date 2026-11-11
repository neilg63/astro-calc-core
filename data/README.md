Not committed (see .gitignore). Fetch with:

    curl -O https://ssd.jpl.nasa.gov/ftp/eph/planets/bsp/de440s.bsp
    curl -O https://ssd.jpl.nasa.gov/ftp/eph/planets/ascii/de440/testpo.440
    python3 -I ../tools/fetch_horizons_spk.py .    # Ceres, Pallas, Juno, Vesta, Chiron, Pholus (1850–2150)

`Ephemeris::open_dir("data")` loads every `.bsp` file here.
