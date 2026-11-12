"""Harvest Kronos and Isis positions and Ceres phenomena from the AstroCalc API.
Usage: python3 -I tools/harvest_hypothetical.py OUT_DIR"""
import datetime as dt, json, os, random, sys, time, urllib.request

API = "https://astro.multifaceted.info"
def get(path):
    req = urllib.request.Request(API + path, headers={"User-Agent": "astro-calc-fixtures"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.load(r)

out = sys.argv[1]
rng = random.Random(8841)
start = dt.datetime(1852, 1, 1)
span = int((dt.datetime(2148, 12, 31) - start).total_seconds())
cases = []
for _ in range(24):
    when = (start + dt.timedelta(seconds=rng.randrange(span))).strftime("%Y-%m-%dT%H:%M:%S")
    d = get(f"/chart-data?dt={when}&loc=0,0&bodies=kr,is&eq=2&topo=0")
    p = get(f"/pheno?dt={when}&bodies=ce")
    cases.append({"utc": when, "jd": d["date"]["jd"], "bodies": d["bodies"], "pheno": p["result"]})
    time.sleep(0.1)
json.dump(cases, open(os.path.join(out, "api_hypothetical.json"), "w"), indent=1)
print(len(cases))
