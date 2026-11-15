#!/usr/bin/env python3
"""Packs the published planetary and lunar series into the binary tables that the `analytic`
feature embeds.

Planets and Sun: VSOP87, version E. Bretagnon & Francou, "Planetary theories in rectangular
and spherical variables: VSOP87 solutions", Astron. Astrophys. 202, 309 (1988); files
VSOP87E.* as distributed by IMCCE (ftp://ftp.imcce.fr/pub/ephem/planets/vsop87/).

Moon: ELP/MPP02. Chapront & Francou, "The lunar theory ELP revisited. Introduction of new
planetary perturbations", Astron. Astrophys. 404, 735 (2003); files ELP_MAIN.S1-3 and
ELP_PERT.S1-3 as distributed by Paris Observatory
(ftp://polac.obspm.fr/pub/2_lunar_solutions/2_elpmpp02/). The constants and the way the
series are assembled follow the authors' note and reference subroutine in the same
directory, using the set fitted to the JPL ephemeris DE405 with the corrections to the
secular terms that bring it close to DE406 over six millennia.

Usage:
    python3 -I tools/build_analytic_tables.py --vsop87 VSOP87_DIR --elp ELP_DIR

VSOP87: a term `A cos(B + C T)` of the series in T^alpha is kept when its largest
contribution over the supported span, A * SPAN^alpha, is at least EPS times the body's mean
distance. The default, 1e-8, keeps about half of the 51,880 terms and adds at most about
0.15" to the theory's own error; `--eps 0` keeps everything.

Layout (little-endian): magic "ACVSOP01", then for each of 9 bodies (Sun, Mercury ... Neptune),
3 coordinates and 6 powers of T a u32 term count (162 counts), then the terms in the same
order, each three f64: A (au), B (rad), C (rad per thousand Julian years).

ELP/MPP02: a term is kept when its amplitude, times T^power at 30 centuries from J2000 for
the terms multiplied by a power of time, is at least `--elp-eps` arcseconds (longitude and
latitude) or the equivalent length at the Moon's distance (distance series). The default,
0.0005", keeps about 6,300 of the 35,901 terms and adds a few hundredths of an arcsecond;
`--elp-eps 0` keeps everything.

Layout (little-endian): magic "ACELP001"; the five polynomial coefficients of the Moon's
mean longitude W1 (f64, radians and radians per Julian century^k); 15 u32 term counts: the
main problem for longitude, latitude and distance, then the perturbations for each of those
three in powers 0-3 of time; 4 bytes of padding; then the terms in that order, each six
f64: amplitude (arcseconds, or km for distance) and the five polynomial coefficients of its
argument (radians), so that a term is `amplitude * t^power * sin(argument(t))`.
"""
import argparse
import math
import os
import struct
import sys

BODIES = ["sun", "mer", "ven", "ear", "mar", "jup", "sat", "ura", "nep"]
# Mean distances (au) used to scale the truncation; the Sun's is its offset from the barycentre.
SCALE = [0.005, 0.3871, 0.7233, 1.0, 1.5237, 5.2026, 9.5547, 19.2181, 30.1096]
SPAN = 4.0  # thousands of Julian years either side of J2000


def read(path):
    """Returns terms[coordinate][alpha] = [(A, B, C), ...] for one VSOP87E file."""
    terms = [[[] for _ in range(6)] for _ in range(3)]
    current = None
    with open(path, encoding="ascii") as f:
        for line in f:
            if line.startswith(" VSOP87"):
                current = terms[int(line[41]) - 1][int(line[59])]
            elif line.strip():
                current.append((float(line[79:97]), float(line[97:111]), float(line[111:131])))
    return terms


def dms(deg, minutes, sec):
    return math.radians(deg + minutes / 60.0 + sec / 3600.0)


def elp_constants():
    """The constants of ELP/MPP02 fitted to DE405: polynomial coefficients (radians, Julian
    centuries) of the lunar and planetary arguments, and the corrections to the amplitudes."""
    rad = 648000.0 / math.pi  # arcseconds per radian
    am, alpha = 0.074801329, 0.002571881
    dtasm, xa = (2.0 * alpha) / (3.0 * am), (2.0 * alpha) / 3.0
    dprec = -0.29965  # correction to the precession constant (IAU 2000A)
    # Corrections fitted to DE405 over 1950-2060.
    dw1_0, dw2_0, dw3_0, deart_0, dperi = -0.07008, 0.20794, -0.07215, -0.00033, -0.00749
    dw1_1, dgam, de, deart_1, dep = -0.35106, 0.00085, -0.00006, 0.00732, 0.00224
    dw2_1, dw3_1, dw1_2 = 0.08017, -0.04317, -0.03743
    w = [
        [dms(218, 18, 59.95571 + dw1_0), (1732559343.73604 + dw1_1) / rad, (-6.8084 + dw1_2) / rad, 0.66040e-2 / rad, -0.31690e-4 / rad],
        [dms(83, 21, 11.67475 + dw2_0), (14643420.3171 + dw2_1) / rad, -38.2631 / rad, -0.45047e-1 / rad, 0.21301e-3 / rad],
        [dms(125, 2, 40.39816 + dw3_0), (-6967919.5383 + dw3_1) / rad, 6.3590 / rad, 0.76250e-2 / rad, -0.35860e-4 / rad],
    ]
    eart = [dms(100, 27, 59.13885 + deart_0), (129597742.29300 + deart_1) / rad, -0.020200 / rad, 0.90000e-5 / rad, 0.15000e-6 / rad]
    peri = [dms(102, 56, 14.45766 + dperi), 1161.24342 / rad, 0.529265 / rad, -0.11814e-3 / rad, 0.11379e-4 / rad]
    # Corrections to the secular terms, for closeness to DE406 over six millennia.
    w[0][3] -= 0.00018865 / rad
    w[0][4] -= 0.00001024 / rad
    w[1][2] += 0.00470602 / rad
    w[1][3] -= 0.00025213 / rad
    w[2][2] -= 0.00261070 / rad
    w[2][3] -= 0.00010712 / rad
    # Corrections to the mean motions of perigee and node implied by the changed constants.
    bp = [[0.311079095, -0.4482398e-2, -0.110248500e-2, 0.1056062e-2, 0.50928e-4],
          [-0.103837907, 0.6682870e-3, -0.129807200e-2, -0.1780280e-3, -0.37342e-4]]
    for i in (1, 2):
        b = bp[i - 1]
        x = w[i][1] / w[0][1]
        y = am * b[0] + xa * b[4]
        w[i][1] += ((x - y) * dw1_1 + (y / am) * deart_1 + w[0][1] * (b[1] * dgam + b[2] * de + b[3] * dep)) / rad
    # Delaunay arguments D, F, l, l'.
    dela = [
        [w[0][k] - eart[k] for k in range(5)],
        [w[0][k] - w[2][k] for k in range(5)],
        [w[0][k] - w[1][k] for k in range(5)],
        [eart[k] - peri[k] for k in range(5)],
    ]
    dela[0][0] += math.pi
    # Mean longitudes of Mercury ... Neptune (the third is the Earth-Moon barycentre).
    lng0 = [(252, 15, 3.216919), (181, 58, 44.758419), (100, 27, 59.138850), (355, 26, 3.642778),
            (34, 21, 5.379392), (50, 4, 38.902495), (314, 3, 4.354234), (304, 20, 56.808371)]
    rate = [538101628.66888, 210664136.45777, 129597742.29300, 68905077.65936,
            10925660.57335, 4399609.33632, 1542482.57845, 786547.89700]
    planets = [[dms(*a), r / rad, 0.0, 0.0, 0.0] for a, r in zip(lng0, rate)]
    zeta = [w[0][0], w[0][1] + (5029.0966 + dprec) / rad, w[0][2], w[0][3], w[0][4]]
    corr = {
        "delnu": (0.55604 + dw1_1) / rad / w[0][1],
        "dele": (0.01789 + de) / rad,
        "delg": (-0.08066 + dgam) / rad,
        "delnp": (-0.06424 + deart_1) / rad / w[0][1],
        "delep": (-0.12879 + dep) / rad,
        "am": am,
        "dtasm": dtasm,
    }
    return w[0], dela, planets, zeta, corr


def build_elp(src, out, eps):
    w1, dela, planets, zeta, c = elp_constants()
    span = 30.0  # Julian centuries from J2000 at which time-dependent terms are weighed
    groups, total = [], 0
    # Main problem: one series per coordinate.
    for iv in (1, 2, 3):
        limit = eps if iv < 3 else eps / 206264.8 * 385000.0
        terms = []
        with open(os.path.join(src, "ELP_MAIN.S%d" % iv), encoding="ascii") as f:
            count = int(f.readline()[25:35])
            for _ in range(count):
                line = f.readline()
                ilu = [int(line[i * 3:i * 3 + 3]) for i in range(4)]
                a = float(line[14:27])
                b = [float(line[27 + i * 12:39 + i * 12]) for i in range(5)]
                tgv = b[0] + c["dtasm"] * b[4]
                if iv == 3:
                    a -= 2.0 * a * c["delnu"] / 3.0
                amp = a + tgv * (c["delnp"] - c["am"] * c["delnu"]) + b[1] * c["delg"] + b[2] * c["dele"] + b[3] * c["delep"]
                arg = [sum(ilu[i] * dela[i][k] for i in range(4)) for k in range(5)]
                if iv == 3:
                    arg[0] += math.pi / 2.0
                terms.append((amp, *arg))
        total += len(terms)
        groups.append([t for t in terms if abs(t[0]) >= limit])
    # Perturbations: per coordinate, one series for each power of time 0-3.
    for iv in (1, 2, 3):
        limit = eps if iv < 3 else eps / 206264.8 * 385000.0
        with open(os.path.join(src, "ELP_PERT.S%d" % iv), encoding="ascii") as f:
            for power in range(4):
                head = f.readline()
                count = int(head[25:35])
                terms = []
                for _ in range(count):
                    line = f.readline()
                    s = float(line[5:25].replace("D", "E"))
                    co = float(line[25:45].replace("D", "E"))
                    ifi = [int(line[45 + i * 3:48 + i * 3]) for i in range(16)]
                    arg = []
                    for k in range(5):
                        v = math.atan2(co, s) % (2.0 * math.pi) if k == 0 else 0.0
                        v += sum(ifi[i] * dela[i][k] for i in range(4))
                        v += sum(ifi[i] * planets[i - 4][k] for i in range(4, 12))
                        v += ifi[12] * zeta[k]
                        arg.append(v)
                    terms.append((math.hypot(co, s), *arg))
                total += len(terms)
                groups.append([t for t in terms if abs(t[0]) * span**power >= limit])
    kept = sum(len(g) for g in groups)
    with open(out, "wb") as f:
        f.write(b"ACELP001")
        f.write(struct.pack("<5d", *w1))
        f.write(struct.pack("<15I4x", *[len(g) for g in groups]))
        for g in groups:
            for t in g:
                f.write(struct.pack("<6d", *t))
    print(f"ELP/MPP02: {kept} of {total} terms kept, {112 + 48 * kept} bytes -> {out}", file=sys.stderr)


def build_vsop87(src, out, eps):
    counts, body_bytes, kept, total = [], bytearray(), 0, 0
    for name, scale in zip(BODIES, SCALE):
        terms = read(os.path.join(src, "VSOP87E." + name))
        for coord in terms:
            for alpha, series in enumerate(coord):
                total += len(series)
                series = [t for t in series if abs(t[0]) * SPAN**alpha >= eps * scale]
                kept += len(series)
                counts.append(len(series))
                for t in series:
                    body_bytes += struct.pack("<3d", *t)
    with open(out, "wb") as f:
        f.write(b"ACVSOP01")
        f.write(struct.pack("<%dI" % len(counts), *counts))
        f.write(body_bytes)
    print(f"VSOP87E: {kept} of {total} terms kept, {8 + 4 * len(counts) + len(body_bytes)} bytes -> {out}", file=sys.stderr)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vsop87", help="directory holding the VSOP87E.* files")
    ap.add_argument("--elp", help="directory holding the ELP_MAIN.S* and ELP_PERT.S* files")
    ap.add_argument("--eps", type=float, default=1e-8, help="VSOP87 truncation (relative)")
    ap.add_argument("--elp-eps", type=float, default=0.0005, help="ELP/MPP02 truncation (arcseconds)")
    ap.add_argument("--out-dir", default="src/analytic")
    args = ap.parse_args()
    if not (args.vsop87 or args.elp):
        ap.error("give --vsop87 and/or --elp")
    if args.vsop87:
        build_vsop87(args.vsop87, os.path.join(args.out_dir, "vsop87e.bin"), args.eps)
    if args.elp:
        build_elp(args.elp, os.path.join(args.out_dir, "elpmpp02.bin"), args.elp_eps)


if __name__ == "__main__":
    main()
