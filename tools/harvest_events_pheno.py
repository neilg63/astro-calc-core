"""Harvest planetary stations, moon phases, phenomena and asteroid positions from the AstroCalc API.
Usage: python3 -I tools/harvest_events_pheno.py OUT_DIR"""
import datetime as dt, json, os, random, sys, time, urllib.request

API = "https://astro.multifaceted.info"
def get(path):
    req = urllib.request.Request(API + path, headers={"User-Agent": "astro_ephem-fixtures"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.load(r)

out = sys.argv[1]
rng = random.Random(5731)
start = dt.datetime(1852, 1, 1)
span = int((dt.datetime(2148, 12, 31) - start).total_seconds())
def rand_dt():
    return (start + dt.timedelta(seconds=rng.randrange(span))).strftime("%Y-%m-%dT%H:%M:%S")

stations = []
for y0 in (2001, 2030):  # the reference table covers 2000–2050
    d = get(f"/planet-stations?dt={y0}-01-01&dt2={y0 + 12}-01-01&bodies=me,ve,ma,ju,sa,ur,ne,pl")
    stations.append({"from": d["start"]["jd"], "to": d["end"]["jd"], "items": d["items"]})
    time.sleep(0.2)
json.dump(stations, open(os.path.join(out, "api_stations.json"), "w"), indent=1)

phases = []
for _ in range(12):
    when = rand_dt()
    d = get(f"/moon-phases?dt={when}&loc=0,0")
    phases.append({"utc": when, "jd": d["date"]["jd"], "phases": d["phases"]})
    time.sleep(0.1)
json.dump(phases, open(os.path.join(out, "api_moon_phases.json"), "w"), indent=1)

pheno, minor = [], []
for _ in range(30):
    when = rand_dt()
    d = get(f"/pheno?dt={when}&bodies=su,mo,me,ve,ma,ju,sa,ur,ne,pl")
    pheno.append({"utc": when, "jd": d["date"]["jd"], "result": d["result"]})
    d = get(f"/chart-data?dt={when}&loc=0,0&bodies=ch,ce,jn&eq=2&topo=0")
    minor.append({"utc": when, "jd": d["date"]["jd"], "bodies": d["bodies"]})
    time.sleep(0.1)
json.dump(pheno, open(os.path.join(out, "api_pheno.json"), "w"), indent=1)
json.dump(minor, open(os.path.join(out, "api_minor_bodies.json"), "w"), indent=1)
print(len(stations), sum(len(i["stations"]) for s in stations for i in s["items"]), len(phases), len(pheno), len(minor))
