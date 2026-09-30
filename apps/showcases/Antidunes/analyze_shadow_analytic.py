#!/usr/bin/env python3
"""Compare the closed-form Kodam+Guo Cylinder x Cylinder result with the GJK/EPA result that was actually used (shadow-analytic CSVs, ModelParameters.shadowAnalyticCylinderLog).

    python3 analyze_shadow_analytic.py <shadow-analytic dir> [--worst 15] [--min-overlap 0.0]

Rows are deduplicated on (mainStep, subCycle, uidA, uidB) (neighbouring ranks compute the same pair). Overlaps are in LU; "err" = closed-form overlap minus GJK overlap (>0: closed form over-estimates).
"""
import argparse, csv, glob, os
ap = argparse.ArgumentParser(); ap.add_argument("dir"); ap.add_argument("--worst", type=int, default=15); ap.add_argument("--tol", type=float, default=0.15, help="|err| above this counts as wrong (5%% of R=3)")
a = ap.parse_args()
import heapq
seen = set(); by = {}; heap = []; nrows = 0; both = gjkonly = cfonly = 0; g = []
def acc(sc):
    return by.setdefault(sc, dict(n=0, wrong=0, sumabs=0.0, under=0.0, over=0.0, flags=0, flagw=0, flipped=0, off=0, sumcf=0.0, sumgjk=0.0, ngjk=0, low=0))
for f in sorted(glob.glob(os.path.join(a.dir, "rank_*.csv"))):
    for r in csv.DictReader(open(f)):
        k = (r["mainStep"], r["subCycle"], r["uidA"], r["uidB"])
        if k in seen: continue                       # neighbouring ranks compute the same pair
        seen.add(k); nrows += 1
        gj, cf = r["gjkContact"] == "1", r["closedFormContact"] == "1"
        both += gj and cf; gjkonly += gj and not cf; cfonly += cf and not gj
        if gj and len(g) < 3000000: g.append(float(r["gjkOverlap"]))
        e = float(r["errCfMinusGjk"]); st = acc(r["scenario"]); st["n"] += 1; st["sumabs"] += abs(e); st["under"] = max(st["under"], -e); st["over"] = max(st["over"], e)
        flag = r["checkFlags"] == "1"; st["flags"] += flag
        if gj:
            go = float(r["gjkOverlap"]); co = float(r["cfOverlap"]) if cf else 0.0
            st["sumcf"] += co; st["sumgjk"] += go
            if go > 0.01: st["ngjk"] += 1; st["low"] += co < 0.5 * go
        if abs(e) > a.tol: st["wrong"] += 1; st["flagw"] += flag
        if cf and gj:
            d = float(r["normalDot"]); st["flipped"] += d < 0; st["off"] += 0 <= d < 0.94
        item = (abs(e), k, r)
        if len(heap) < a.worst: heapq.heappush(heap, (item[0], item[1], item[2]))
        elif item[0] > heap[0][0]: heapq.heapreplace(heap, (item[0], item[1], item[2]))
print(f"{nrows} distinct (step, sub-cycle, pair) rows where GJK or the closed form saw a contact")
print(f"both contact: {both}   GJK contact but closed form finds nothing: {gjkonly}   closed-form contact but GJK none: {cfonly}\n")
print(f"{'closed-form scenario':22s} {'n':>9s} {'|err|>tol':>10s} {'%':>6s} {'mean |err|':>10s} {'max under':>10s} {'max over':>9s} {'check flags':>11s} {'flag among wrong':>16s} {'normal flipped':>14s} {'normal >20deg':>13s}")
for sc, st in sorted(by.items(), key=lambda kv: -kv[1]["n"]):
    print(f"{sc:22s} {st['n']:9d} {st['wrong']:10d} {100*st['wrong']/st['n']:6.2f} {st['sumabs']/st['n']:10.5f} {st['under']:10.3f} {st['over']:9.3f} {st['flags']:11d} {(str(st['flagw'])+'/'+str(st['wrong'])) if st['wrong'] else '-':>16s} {st['flipped']:14d} {st['off']:13d}")
print(f"\nRELATIVE view (contact force ~ overlap): per closed-form scenario  sum(closed-form overlap)/sum(GJK overlap)  and  share of contacts (GJK overlap > 0.01 LU) where the closed form gives < 50 % of GJK's overlap:")
for sc, st in sorted(by.items(), key=lambda kv: -kv[1]["n"]):
    if st["sumgjk"] > 0: print(f"   {sc:22s} overlap ratio {st['sumcf']/st['sumgjk']:7.3f}   <50% of GJK: {st['low']}/{st['ngjk']} ({100*st['low']/max(1,st['ngjk']):.1f}%)")
g.sort()
if g: print(f"\nGJK overlap distribution (rows with a GJK contact): n={len(g)} median={g[len(g)//2]:.5f} p90={g[int(len(g)*.9)]:.5f} p99={g[int(len(g)*.99)]:.4f} max={g[-1]:.4f}")
print(f"\nworst {a.worst} disagreements (largest |closed-form overlap - GJK overlap|):")
for _, k, r in sorted(heap, key=lambda t: -t[0]):
    print(f"  step {r['mainStep']}.{r['subCycle']} pair ..{r['uidA'][-5:]}/..{r['uidB'][-5:]}: GJK overlap {float(r['gjkOverlap']):.4f} | closed form [{r['scenario']}] overlap {float(r['cfOverlap']):.4f} (needs {float(r['neededAlongCfNormal']):.3f} along its own normal; check {'FLAGS' if r['checkFlags']=='1' else 'ok'}) | normal dot {r['normalDot']}")
