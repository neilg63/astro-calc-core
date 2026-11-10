"""Harvest reference positions from the AstroCalc API into a JSON fixture.

Usage: python3 -I tools/harvest_api_fixtures.py [out.json]
Sends whole-second ISO times (the API rounds to seconds) and records the JD it echoes.
"""
import datetime as dt, json, random, sys, time, urllib.request

BASE = "https://astro.multifaceted.info/chart-data"
BODIES = "su,mo,me,ve,ma,ju,sa,ur,ne,pl"
LOCS = [(45.22, 11.228, 0), (-33.865, 151.209, 58), (40.7128, -74.006, 10), (69.649, 18.955, 0), (-1.2921, 36.8219, 1795)]

def fetch(when, loc, topo):
    q = f"{BASE}?dt={when}&loc={loc[0]},{loc[1]},{loc[2]}&bodies={BODIES}&eq=2&topo={topo}"
    req = urllib.request.Request(q, headers={"User-Agent": "astro_ephem-fixtures"})
    with urllib.request.urlopen(req, timeout=30) as r:
        return json.load(r)

def main():
    out = sys.argv[1] if len(sys.argv) > 1 else "fixtures.json"
    rng = random.Random(440)
    start, end = dt.datetime(1851, 1, 1), dt.datetime(2148, 12, 31)
    span = int((end - start).total_seconds())
    cases = []
    for i in range(36):
        when = (start + dt.timedelta(seconds=rng.randrange(span))).strftime("%Y-%m-%dT%H:%M:%S")
        loc = LOCS[i % len(LOCS)]
        for topo in (0, 1):
            d = fetch(when, loc, topo)
            cases.append({"utc": when, "jd": d["date"]["jd"], "loc": list(loc), "topo": topo,
                          "bodies": d["bodies"]})
            time.sleep(0.15)
    json.dump(cases, open(out, "w"), indent=1)
    print(f"wrote {len(cases)} cases to {out}")

main()
