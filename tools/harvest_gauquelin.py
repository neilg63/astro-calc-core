"""Harvest Gauquelin sector cusps from the AstroCalc API, with the angles of each case.
Usage: python3 -I tools/harvest_gauquelin.py gauquelin.json"""
import datetime as dt, json, random, sys, time, urllib.request

def get(url):
    req = urllib.request.Request(url, headers={"User-Agent": "astro-calc-fixtures"})
    with urllib.request.urlopen(req, timeout=30) as r:
        return json.load(r)

API = "https://astro.multifaceted.info"
LATS = [0.0, 12.97, 23.4, 45.22, -33.865, 51.5, 60.17, 64.5, -55.0, 70.0, -72.5, 85.0]
rng = random.Random(1928)
start = dt.datetime(1851, 1, 1)
span = int((dt.datetime(2049, 12, 31) - start).total_seconds())
cases = []
for i in range(48):
    when = (start + dt.timedelta(seconds=rng.randrange(span))).strftime("%Y-%m-%dT%H:%M:%S")
    lat, lng = LATS[i % len(LATS)], round(rng.uniform(-180, 180), 4)
    d = get(f"{API}/houses?dt={when}&loc={lat},{lng}&hsys=G&aya=0")
    sets = d["houseSets"]["sets"]
    cases.append({"utc": when, "jd": d["date"]["jd"], "lat": lat, "lng": lng,
                  "ascendant": d["houseSets"]["points"]["ascendant"], "mc": d["houseSets"]["points"]["mc"],
                  "cusps": sets[0]["houses"]})
    time.sleep(0.15)
json.dump(cases, open(sys.argv[1], "w"), indent=1)
print(len(cases), sorted({len(c["cusps"]) for c in cases}))
