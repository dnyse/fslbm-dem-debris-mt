#!/usr/bin/env python3
"""Analyse a fixed-heel righting-moment scan (run_65_bogner_righting_scan_s16.sbatch).

    python3 analyze_bogner_righting_scan.py <archive_dir> [--manifest prms/bogner_righting_scan_s16_manifest.json] [--plot out.png]

<archive_dir> contains one sub-folder per heel angle (a00, a05, ..., a26p565, a30), each with floater-history/rank_<N>_state.csv.
For every angle and box (identified by its X position: rho 0.75 / 0.5 / 0.25 left to right) it averages the raw MEM hydroTorqueY and hydroForceZ of the
OWNING rank's copy over the logged window, and compares with the ideal hydrostatic reference in the manifest.

Sign convention: a vertical buoyancy force F at horizontal offset d = x_B - x_G gives a torque about +Y of  tau_y = -F*d, so the ideal MEM torque is
-M_ideal_signed_xz (the manifest stores +F*d). The rho=0.75 vs rho=0.25 comparison does not depend on this convention.
"""
import argparse, csv, glob, json, os, sys
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("archive")
ap.add_argument("--manifest", default=os.path.join(os.path.dirname(__file__), "prms", "bogner_righting_scan_s16_manifest.json"))
ap.add_argument("--plot", default="")
a = ap.parse_args()
man = json.load(open(a.manifest)); DENS = man["densities"]; W0, W1 = man["window"]

def read_state(folder):
    rows = []
    for f in sorted(glob.glob(os.path.join(folder, "floater-history", "rank_*_state.csv"))):
        with open(f) as fh:
            for r in csv.DictReader(fh):
                if r["isLocal"] != "1": continue                   # hydro force/torque are only valid on the owning rank's local copy
                rows.append((int(r["mainStep"]), int(r["subCycle"]), int(r["uid"]), float(r["posX"]), float(r["hydroTorqueY"]), float(r["hydroForceZ"]),
                             float(r["hydroForceX"]), float(r["hydroTorqueX"]), float(r["hydroTorqueZ"])))
    return rows

results = {}
for tag in man["tags"]:
    d = os.path.join(a.archive, tag)
    if not os.path.isdir(d):
        alt = sorted(glob.glob(d + "_attempt*"))
        if not alt: print(f"[{tag}] missing"); continue
        d = alt[-1]
    rows = read_state(d)
    if not rows: print(f"[{tag}] no rows in floater-history"); continue
    rows = [r for r in rows if r[1] == 0 and W0 <= r[0] <= W1]     # one row per main step (the hydro value is fixed within a step)
    uids = sorted({r[2] for r in rows}, key=lambda u: np.mean([r[3] for r in rows if r[2] == u]))
    ang = float(next(c["alpha"] for c in man["cases"] if c["alpha"] is not None and man["tags"][man["angles"].index(c["alpha"])] == tag))
    for bi, uid in enumerate(uids):
        rr = np.array([(r[4], r[5], r[6], r[7], r[8]) for r in rows if r[2] == uid])
        ideal = next(c for c in man["cases"] if c["alpha"] == ang and c["box_index"] == bi)
        results[(ang, DENS[bi])] = dict(My=rr[:, 0].mean(), My_sd=rr[:, 0].std(), Fz=rr[:, 1].mean(), Fz_sd=rr[:, 1].std(), Fx=rr[:, 2].mean(),
                                        Tx=rr[:, 3].mean(), Tz=rr[:, 4].mean(), n=len(rr), ideal_My=-ideal["M_ideal_signed_xz"], ideal_Fz=ideal["F_buoy_ideal"])

if not results: sys.exit("no results")
angles = sorted({k[0] for k in results})
print(f"window steps {W0}..{W1}; torque in LU (about the box centre, about +Y); ideal My = -F*d\n")
print(f"{'heel':>7} | {'rho':>4} | {'My meas':>9} {'+-sd':>7} {'My ideal':>9} {'meas-ideal':>10} | {'Fz meas':>8} {'Fz ideal':>8} {'Fz err %':>8}")
for ang in angles:
    for rho in DENS:
        r = results.get((ang, rho))
        if not r: continue
        print(f"{ang:7.3f} | {rho:4.2f} | {r['My']:9.4f} {r['My_sd']:7.4f} {r['ideal_My']:9.4f} {r['My']-r['ideal_My']:10.4f} | {r['Fz']:8.4f} {r['ideal_Fz']:8.4f} {100*(r['Fz']/r['ideal_Fz']-1):8.2f}")
    print()
print("KEY TEST -- ideal theory says the rho=0.75 and rho=0.25 torque curves are IDENTICAL:")
print(f"{'heel':>7} | {'My(0.75)':>9} {'My(0.25)':>9} {'diff':>8} | {'meas-ideal .75':>14} {'meas-ideal .25':>14}")
for ang in angles:
    a75, a25 = results.get((ang, 0.75)), results.get((ang, 0.25))
    if a75 and a25:
        print(f"{ang:7.3f} | {a75['My']:9.4f} {a25['My']:9.4f} {a75['My']-a25['My']:8.4f} | {a75['My']-a75['ideal_My']:14.4f} {a25['My']-a25['ideal_My']:14.4f}")
# where would each measured curve cross zero (linear interpolation)? the ideal crossing is 26.565 for rho 0.25/0.75, 45 for 0.5
def zero_cross(rho):
    xs = [(ang, results[(ang, rho)]['My']) for ang in angles if (ang, rho) in results]
    for (x0, y0), (x1, y1) in zip(xs, xs[1:]):
        if y0 == 0 or y0 * y1 < 0:
            return x0 - y0 * (x1 - x0) / (y1 - y0)
    return None
crossings = {}
for rho in (0.75, 0.25):
    zc = zero_cross(rho)
    crossings[rho] = None if zc is None else round(float(zc), 2)
print("\nmeasured zero crossing of My (interpolated; ideal 26.565 for 0.75 and 0.25):", crossings)
if a.plot:
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    for rho, c in ((0.75, "C3"), (0.5, "C2"), (0.25, "C0")):
        xs = [ang for ang in angles if (ang, rho) in results]
        ax[0].errorbar(xs, [results[(x, rho)]["My"] for x in xs], [results[(x, rho)]["My_sd"] for x in xs], color=c, marker="o", label=f"MEM rho={rho}")
        ax[0].plot(xs, [results[(x, rho)]["ideal_My"] for x in xs], color=c, ls="--", lw=1)
    ax[0].axhline(0, color="k", lw=.5); ax[0].set_xlabel("heel angle (deg)"); ax[0].set_ylabel("torque about +Y (LU)"); ax[0].legend(); ax[0].set_title("measured (solid) vs ideal (dashed)")
    xs = [ang for ang in angles if (ang, 0.75) in results and (ang, 0.25) in results]
    ax[1].plot(xs, [results[(x, .75)]["My"] - results[(x, .25)]["My"] for x in xs], "ko-"); ax[1].axhline(0, color="k", lw=.5)
    ax[1].set_xlabel("heel angle (deg)"); ax[1].set_ylabel("My(0.75) - My(0.25)  (ideal: 0)"); ax[1].set_title("symmetry test")
    plt.tight_layout(); plt.savefig(a.plot, dpi=130); print("wrote", a.plot)
