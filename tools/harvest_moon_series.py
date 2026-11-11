"""Harvest geocentric Moon/Sun longitudes from the AstroCalc API at 1 Jan and 1 Jul 00:00 UTC
each year, for inferring the ΔT model behind the reference values.
Usage: python3 -I tools/harvest_moon_series.py out.json [start_year end_year]"""
import json, sys, time, urllib.request
out = sys.argv[1]
y0, y1 = (int(sys.argv[2]), int(sys.argv[3])) if len(sys.argv) > 3 else (1850, 2149)
rows = []
for y in range(y0, y1 + 1):
    for m in ("01", "07"):
        q = f"https://astro.multifaceted.info/positions?dt={y:04d}-{m}-01T00:00:00&loc=0,0&bodies=su,mo&eq=0&topo=0"
        req = urllib.request.Request(q, headers={"User-Agent": "astro-calc-fixtures"})
        with urllib.request.urlopen(req, timeout=30) as r:
            d = json.load(r)
        rows.append({"utc": f"{y:04d}-{m}-01T00:00:00", "jd": d["date"]["jd"],
                     "su": d["longitudes"]["su"], "mo": d["longitudes"]["mo"]})
        time.sleep(0.1)
json.dump(rows, open(out, "w"))
print(len(rows), rows[0], rows[-1])
