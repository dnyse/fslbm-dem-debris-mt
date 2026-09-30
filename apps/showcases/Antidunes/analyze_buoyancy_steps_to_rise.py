#!/usr/bin/env python3
"""Steps needed for the floater centre to rise the SAME RELATIVE distance at R1 / R2 / R4 (MEM-only vs analytic buoyancy).

[Claude, 2026-09-21, corrected 2026-09-22] The buoyancy_compare_* runs (R1/R2/R4 x sphere/cylinder/box x mem/analytic; scenario identical in
physical units, every LENGTH scaled by R -- box halfExtent 8/16/32, domain likewise) all start the floater at rest at the same relative height.

CORRECTED 2026-09-22 (was wrong before): the LU threshold D must scale WITH R, D_R = D_base * R, not stay fixed across resolutions. A fixed
absolute D (e.g. 10 LU at every R) is NOT the same journey at different R -- at R1 (box halfExtent 8) it is 1.25x the box's own half-size, at
R4 (box halfExtent 32) only 0.31x, so the R4 numbers understated how far it had really travelled relative to its own scale. Scaling D by R makes
'D_base * R lattice units' the same FRACTION of the floater's own size (and the same PHYSICAL distance in real units) at every resolution, so the
resulting step counts are a genuine resolution-convergence comparison -- if the method converges, they should scale close to R^2 (since dt ~
1/R^2 at fixed physical viscosity/gravity, matching the prms' own timesteps scaling 2600 -> 10400 -> 41600).

usage: analyze_buoyancy_steps_to_rise.py <root> [--dbase D] where <root> contains the three archive folders
       antidunes-buoyancy-r{1,2,4}-compare-all*/<shape>_<method>/eval/antidunes.sqlite (vault: $HPCVAULT/antidunes-archive/).
       --dbase sets the R1 threshold (default 3,6,10,14, i.e. D_R4 = 12,24,40,56).
"""
import argparse, glob, re, sqlite3

ap = argparse.ArgumentParser()
ap.add_argument("root")
ap.add_argument("--dbase", default="3,6,10,14", help="comma-separated D values AT R=1 (scaled by R for R=2,4)")
a = ap.parse_args()
D_BASE = [float(v) for v in a.dbase.split(",")]

S = {}
for f in glob.glob(f"{a.root}/antidunes-buoyancy-r*-compare-all*/*/eval/antidunes.sqlite"):
    R = int(re.search(r"-r(\d)-", f).group(1)); run = f.split("/")[-3]
    S[(R, run)] = sqlite3.connect(f"file:{f}?mode=ro", uri=True).execute(
        "select timestep,maxParticlePos_z from runs order by timestep").fetchall()

def steps_to(rows, D):
    z0 = rows[0][1]; prev = rows[0]
    for t, z in rows[1:]:
        if z - z0 >= D:
            t0, z_0 = prev
            return t0 + (D - (z_0 - z0)) / (z - z_0) * (t - t0)
        prev = (t, z)
    return None

for shape in ("sphere", "cylinder", "box"):
    print(f"\n{shape}: steps to rise D_base*R lattice units (same RELATIVE journey at every R)")
    print(" D_base | " + " | ".join(f"R{R} mem  R{R} ana (D_R{R})" for R in (1, 2, 4)))
    for Db in D_BASE:
        cells = []
        for R in (1, 2, 4):
            D = Db * R
            for m in ("mem", "analytic"):
                v = steps_to(S[(R, f"{shape}_{m}")], D) if (R, f"{shape}_{m}") in S else None
                cells.append("   -   " if v is None else f"{v:7.0f}")
        print(f" {Db:6.1f} | " + "   ".join(cells))

    # explicit R2/R1 and R4/R2 step-count ratios at the largest common D_base, for the expected-R^2 check
    Db = D_BASE[-1]
    print(f"  ratios at D_base={Db:g} (expect ~R^2 = 4 and 4 if the method converges cleanly):")
    for m in ("mem", "analytic"):
        vals = [steps_to(S[(R, f"{shape}_{m}")], Db * R) if (R, f"{shape}_{m}") in S else None for R in (1, 2, 4)]
        if all(v is not None for v in vals):
            print(f"    {m:9s}: steps={[round(v) for v in vals]}  R2/R1={vals[1]/vals[0]:.2f}  R4/R2={vals[2]/vals[1]:.2f}")
        else:
            print(f"    {m:9s}: missing data {vals}")
