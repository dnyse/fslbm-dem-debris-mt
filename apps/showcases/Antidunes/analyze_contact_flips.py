#!/usr/bin/env python3
"""Scan a floater-history folder's per-rank contact CSVs (rank_<N>.csv) for pairs whose contact overlap or detection path JUMPS between consecutive sub-cycles.

    python3 analyze_contact_flips.py <floater-history dir> [--ratio 4] [--min-pen 0.3] [--top 40]

Overlap = -penetrationDepth (mesa_pd convention: negative = overlapping). A "flip" = the overlap changes by more than --ratio (up or down) with the larger value above --min-pen LU,
or the detection path family changes (closed-form analytical scenario <-> GJK/EPA fallback) while the pair stays in contact.
"""
import argparse, csv, glob, os, collections
ap = argparse.ArgumentParser(); ap.add_argument("dir"); ap.add_argument("--ratio", type=float, default=4.0); ap.add_argument("--min-pen", type=float, default=0.3); ap.add_argument("--top", type=int, default=40)
a = ap.parse_args()
rows = collections.defaultdict(list)          # pair -> [(t, sub, ov, path, rank, kind)]
for f in sorted(glob.glob(os.path.join(a.dir, "rank_[0-9]*.csv"))):
    if f.endswith("_state.csv"): continue
    for r in csv.DictReader(open(f)):
        if r["kind"] != "contact" or r["shapeTypeA"] != "6" or r["shapeTypeB"] != "6": continue     # rod x rod only
        key = (min(int(r["uidA"]), int(r["uidB"])), max(int(r["uidA"]), int(r["uidB"])))
        rows[key].append((int(r["mainStep"]) * 100 + int(r["subCycle"]), int(r["mainStep"]), int(r["subCycle"]), -float(r["penetrationDepth"]), r["detectionPath"], int(r["computingRank"])))
def fam(p): return "GJK" if p.startswith("GJK") else "analytic:" + (p.split("scenario=")[-1].rstrip(")") if "scenario=" in p else p[:20])
events = []; nn = 0; deep = 0; maxov = 0
for key, lst in rows.items():
    lst.sort(); nn += 1
    for (t0, s0, c0, o0, p0, k0), (t1, s1, c1, o1, p1, k1) in zip(lst, lst[1:]):
        maxov = max(maxov, o1)
        if t1 - t0 > 1 and not (s1 == s0 + 1 and c1 == 0 and c0 == 9): continue        # only consecutive sub-cycles
        hi, lo = max(o0, o1), max(min(o0, o1), 1e-9)
        pathflip = fam(p0).startswith("GJK") != fam(p1).startswith("GJK")
        if (hi > a.min_pen and hi / lo > a.ratio) or (pathflip and hi > a.min_pen):
            events.append((hi, key, s0, c0, o0, fam(p0), k0, s1, c1, o1, fam(p1), k1))
print(f"rod-rod pairs with contact rows: {nn}; max overlap seen {maxov:.3f}; flip events (ratio>{a.ratio}, overlap>{a.min_pen}): {len(events)}")
paths = collections.Counter(fam(p) for lst in rows.values() for (_, _, _, _, p, _) in lst)
print("path mix over all rod-rod contact rows:", dict(paths))
for e in sorted(events, reverse=True)[:a.top]:
    hi, key, s0, c0, o0, f0, k0, s1, c1, o1, f1, k1 = e
    print(f"pair {key[0]%100000}/{key[1]%100000}: step {s0}.{c0} ov={o0:.4f} [{f0}, rank {k0}]  ->  step {s1}.{c1} ov={o1:.4f} [{f1}, rank {k1}]")
