"""Harvest lunar nodes, azimuth/altitude and rise/set fixtures from the AstroCalc API.
Usage: python3 -I tools/harvest_nodes_riseset.py OUT_DIR"""
import datetime as dt, json, os, random, sys, time, urllib.request

API = "https://astro.multifaceted.info"
def get(path):
    req = urllib.request.Request(API + path, headers={"User-Agent": "astro-calc-fixtures"})
    with urllib.request.urlopen(req, timeout=30) as r:
        return json.load(r)

out = sys.argv[1]
rng = random.Random(2450094)
start = dt.datetime(1852, 1, 1)
span = int((dt.datetime(2148, 12, 31) - start).total_seconds())
def rand_case():
    when = (start + dt.timedelta(seconds=rng.randrange(span))).strftime("%Y-%m-%dT%H:%M:%S")
    return when, round(rng.uniform(-58, 62), 3), round(rng.uniform(-180, 180), 3), rng.choice([0, 0, 50, 1200])

charts = []
for _ in range(30):
    when, lat, lng, alt = rand_case()
    for topo in (0, 1):
        d = get(f"/chart-data?dt={when}&loc={lat},{lng},{alt}&bodies=su,mo,ma,ju,ra,mn&eq=2&topo={topo}")
        charts.append({"utc": when, "jd": d["date"]["jd"], "loc": [lat, lng, alt], "topo": topo, "bodies": d["bodies"]})
        time.sleep(0.1)
json.dump(charts, open(os.path.join(out, "api_nodes_azalt.json"), "w"), indent=1)

risesets = []
for _ in range(20):
    when, lat, lng, alt = rand_case()
    case = {"utc": when, "loc": [lat, lng, alt], "modes": {}}
    for mode in range(8):
        d = get(f"/rise-set-times?dt={when}&loc={lat},{lng},{alt}&bodies=su,mo,ma,ju,sa&days=1&mode={mode}")
        case["jd"] = d["date"]["jd"]
        case["modes"][str(mode)] = {s["key"]: {i["key"]: i["value"] for i in s["items"]} for s in d["sets"]}
        time.sleep(0.1)
    p = get(f"/positions?dt={when}&loc={lat},{lng},{alt}&topo=1")
    case["positions"] = {"longitudes": p["longitudes"], "sunRiseSets": p["sunRiseSets"], "moonRiseSets": p["moonRiseSets"]}
    risesets.append(case)
json.dump(risesets, open(os.path.join(out, "api_rise_set.json"), "w"), indent=1)
print(len(charts), len(risesets))
