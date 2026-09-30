#!/usr/bin/env python3
"""Trace particles through a floater-history folder (ParticleKernels.h FloaterHistoryLogger CSVs).

    python3 analyze_stale_particle.py <floater-history dir> <uid> [<uid> ...] [--from S] [--to S]

Per (mainStep, subCycle) prints the OWNING copy's (isLocal=1) pos / vel / total force / torque / raw MEM hydro force+torque (hydro values are only valid on the owner),
the number of ranks that hold a copy, and every accepted contact/lubrication row (rank_<N>.csv) involving the uid (partner uid, computing rank, detection path, penetration, force on it).
Used for job 4250846 (step 5714, [STALE-PARTICLE-IDX] uids ...727/...726): the stale index is downstream of a runaway, this shows what started it.
"""
import argparse, csv, glob, os, collections
ap = argparse.ArgumentParser(); ap.add_argument("dir"); ap.add_argument("uids", nargs="+", type=int)
ap.add_argument("--from", dest="s0", type=int, default=0); ap.add_argument("--to", dest="s1", type=int, default=10**9)
ap.add_argument("--first-only", action="store_true", help="only print rows from sub-cycle 0..2 and the steps around the first anomaly")
a = ap.parse_args(); U = set(a.uids)
state = collections.defaultdict(list); copies = collections.defaultdict(set); contacts = collections.defaultdict(list)
for f in sorted(glob.glob(os.path.join(a.dir, "rank_*_state.csv"))):
    for r in csv.DictReader(open(f)):
        u = int(r["uid"]); s = int(r["mainStep"])
        if u in U and a.s0 <= s <= a.s1:
            copies[(u, s, int(r["subCycle"]))].add((int(r["observingRank"]), int(r["isGhost"])))
            if r["isLocal"] == "1": state[u].append(r)
for f in sorted(glob.glob(os.path.join(a.dir, "rank_[0-9]*.csv"))):
    if f.endswith("_state.csv"): continue
    for r in csv.DictReader(open(f)):
        s = int(r["mainStep"])
        if not (a.s0 <= s <= a.s1): continue
        for side, other in (("A", "B"), ("B", "A")):
            u = int(r["uid" + side])
            if u in U: contacts[u].append((s, int(r["subCycle"]), r, side, other))
def v(r, *k): return "(" + ",".join(f"{float(r[x]):.4g}" for x in k) + ")"
for u in a.uids:
    print(f"\n================ uid {u}: {len(state[u])} owner-copy rows, {len(contacts[u])} contact rows")
    for r in sorted(state[u], key=lambda r: (int(r["mainStep"]), int(r["subCycle"]))):
        s, c = int(r["mainStep"]), int(r["subCycle"])
        print(f"step {s} sc {c:2d} owner={r['ownerRank']:>3} copies={len(copies[(u,s,c)])} pos={v(r,'posX','posY','posZ')} vel={v(r,'velX','velY','velZ')} angVel={v(r,'angVelX','angVelY','angVelZ')} "
              f"F={v(r,'forceX','forceY','forceZ')} T={v(r,'torqueX','torqueY','torqueZ')} Fhydro={v(r,'hydroForceX','hydroForceY','hydroForceZ')} Thydro={v(r,'hydroTorqueX','hydroTorqueY','hydroTorqueZ')}")
    print("---- contacts (rows are what the computing rank accepted)")
    for s, c, r, side, other in sorted(contacts[u], key=lambda t: (t[0], t[1])):
        fo = v(r, f"forceOn{side}_x", f"forceOn{side}_y", f"forceOn{side}_z")
        print(f"step {s} sc {c:2d} kind={r['kind']} rank={r['computingRank']} partner={r['uid'+other]} (shape {r['shapeType'+other]}, ghost={r['isGhost'+other]}) path={r['detectionPath'][:38]} pen={float(r['penetrationDepth']):.4g} F_on_me={fo}")
