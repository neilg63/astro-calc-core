"""Harvest house cusps/points (all systems) and ayanamsha values from the AstroCalc API.
Usage: python3 -I tools/harvest_houses_aya.py houses.json ayanamshas.json"""
import datetime as dt, json, random, sys, time, urllib.request

def get(url):
    req = urllib.request.Request(url, headers={"User-Agent": "astro_ephem-fixtures"})
    with urllib.request.urlopen(req, timeout=30) as r:
        return json.load(r)

API = "https://astro.multifaceted.info"
SYSTEMS = "W,E,A,V,O,P,K,B,C,M,R,T,X,H,U"
LATS = [0.0, 12.97, 23.4, 45.22, -33.865, 51.5, 60.17, 64.5, -55.0, 70.0]
rng = random.Random(1996)
start = dt.datetime(1851, 1, 1)
span = int((dt.datetime(2148, 12, 31) - start).total_seconds())
houses = []
for i in range(60):
    when = (start + dt.timedelta(seconds=rng.randrange(span))).strftime("%Y-%m-%dT%H:%M:%S")
    lat, lng = LATS[i % len(LATS)], round(rng.uniform(-180, 180), 4)
    d = get(f"{API}/houses?dt={when}&loc={lat},{lng}&hsys={SYSTEMS}&aya=0")
    houses.append({"utc": when, "jd": d["date"]["jd"], "lat": lat, "lng": lng,
                   "points": d["houseSets"]["points"], "sets": d["houseSets"]["sets"]})
    time.sleep(0.1)
json.dump(houses, open(sys.argv[1], "w"), indent=1)

ayas = []
for y in range(1850, 2150, 2):
    for m in ("01", "08"):
        d = get(f"{API}/chart-data?dt={y}-{m}-01T00:00:00&loc=0,0&bodies=su&aya=all")
        ayas.append({"utc": f"{y}-{m}-01T00:00:00", "jd": d["date"]["jd"],
                     "values": {a["key"]: a["value"] for a in d["ayanamshas"]}})
        time.sleep(0.1)
json.dump(ayas, open(sys.argv[2], "w"), indent=1)
print(len(houses), len(ayas), [s["system"] for s in houses[0]["sets"]])
