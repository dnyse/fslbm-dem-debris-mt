#!/usr/bin/env python3
"""Rod-rod penetration statistics from a floater-history folder (owner copies, sub-cycle 0), using the EXACT convex penetration depth
    PD(A,B) = min over unit n of [ h_A(n) + h_B(n) - d.n ],   h(n) = (L/2)|n.a| + R sqrt(1-(n.a)^2)
(see verify_cylinder_penetration_truth.py / CylinderClosedFormRangeTest.cpp).

    python3 analyze_rod_overlap.py <floater-history dir> [--steps 5900,5930,5960] [--R 3 --L 18]
"""
import argparse, csv, glob, os, collections
import numpy as np
ap = argparse.ArgumentParser(); ap.add_argument("dir"); ap.add_argument("--steps", default=""); ap.add_argument("--R", type=float, default=3.0); ap.add_argument("--L", type=float, default=18.0)
a = ap.parse_args(); R, H = a.R, a.L / 2
rows = collections.defaultdict(list)
want = set(int(x) for x in a.steps.split(",")) if a.steps else None
for f in sorted(glob.glob(os.path.join(a.dir, "rank_*_state.csv"))):
    for r in csv.DictReader(open(f)):
        if r["subCycle"] != "0" or r["isLocal"] != "1" or r["shapeType"] != "6": continue
        s = int(r["mainStep"])
        if want is None or s in want: rows[s].append(r)
N = 4000; i = np.arange(N) + 0.5; z = 1 - 2 * i / N; ph = np.pi * (3 - np.sqrt(5)) * np.arange(N); rr = np.sqrt(1 - z * z)
DIRS = np.stack([rr * np.cos(ph), rr * np.sin(ph), z], 1)
def axis(r):
    w, x, y, zq = [float(r[k]) for k in ("quatR", "quatI", "quatJ", "quatK")]; n = np.sqrt(w*w + x*x + y*y + zq*zq); w, x, y, zq = w/n, x/n, y/n, zq/n
    return np.array([2*(x*zq + w*y), 2*(y*zq - w*x), 1 - 2*(x*x + y*y)])
def h(a_):
    c = np.abs(DIRS @ a_); return H * c + R * np.sqrt(np.clip(1 - c * c, 0, None))
print(f"{'step':>6} {'rods':>5} {'pairs PD>0':>10} {'PD>0.1R':>8} {'PD>0.3R':>8} {'PD>1R':>7} {'max PD [LU]':>11} {'mean PD of overlapping':>22}")
for s in sorted(rows):
    P = np.array([[float(r[k]) for k in ("posX", "posY", "posZ")] for r in rows[s]]); A = np.array([axis(r) for r in rows[s]]); HH = [h(x) for x in A]
    pds = []
    for i_ in range(len(P)):
        for j_ in range(i_ + 1, len(P)):
            d = P[j_] - P[i_]
            if np.linalg.norm(d) > 2 * np.hypot(R, H): continue
            pd = np.min(HH[i_] + HH[j_] - DIRS @ d)
            if pd > 0: pds.append(pd)
    pds = np.array(pds) if pds else np.zeros(1)
    print(f"{s:6d} {len(P):5d} {int((pds>0).sum()):10d} {int((pds>0.1*R).sum()):8d} {int((pds>0.3*R).sum()):8d} {int((pds>R).sum()):7d} {pds.max():11.3f} {pds[pds>0].mean() if (pds>0).any() else 0:22.3f}")
