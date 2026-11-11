"""Download small-body SPK files (type 21) from JPL Horizons into data/.
Usage: python3 -I tools/fetch_horizons_spk.py data [start stop]
Bodies: Ceres, Pallas, Juno, Vesta, Chiron, Pholus."""
import base64, json, os, sys, urllib.parse, urllib.request

BODIES = {"ceres": "1;", "pallas": "2;", "juno": "3;", "vesta": "4;", "chiron": "2060;", "pholus": "5145;"}
out = sys.argv[1]
start, stop = (sys.argv[2], sys.argv[3]) if len(sys.argv) > 3 else ("1850-01-01", "2150-01-01")
for name, cmd in BODIES.items():
    q = urllib.parse.urlencode({"format": "json", "COMMAND": f"'{cmd}'", "OBJ_DATA": "NO", "MAKE_EPHEM": "YES",
                                "EPHEM_TYPE": "SPK", "START_TIME": f"'{start}'", "STOP_TIME": f"'{stop}'"})
    req = urllib.request.Request("https://ssd.jpl.nasa.gov/api/horizons.api?" + q, headers={"User-Agent": "astro-calc"})
    with urllib.request.urlopen(req, timeout=300) as r:
        d = json.load(r)
    if "spk" not in d:
        print(name, "FAILED:", (d.get("result") or d.get("error") or "")[:300])
        continue
    path = os.path.join(out, f"{name}.bsp")
    open(path, "wb").write(base64.b64decode(d["spk"]))
    print(name, d.get("spk_file_id"), os.path.getsize(path), "bytes")
