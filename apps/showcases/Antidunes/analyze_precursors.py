#!/usr/bin/env python3
"""Summarise the fluid-side precursor logs of a run (class (B) pore instability, THESIS_PLAN.md).

    python3 analyze_precursors.py <run .out file>

Prints: the earliest [POP-PRECURSOR] / [HIGH-VELOCITY-PRECURSOR] / [MEM-GARBAGE-LINK] / STALE-PARTICLE-IDX events (step, phase, cell class, world position, max|f|, rho, fill, pinched-pair count and
obstacle-neighbour count), and how many of the flagged cells were pinched (>= 1 opposite obstacle pair) -- the direct test of the pinched-cell hypothesis.
"""
import re, sys, collections
pop = re.compile(r"\[POP-PRECURSOR t=(\d+) after=(\w+)\] (\w+\s*) cell block-local=\(([\d,]+)\) world=\(([\d.,]+)\) max\|f\|=([\d.e+-]+) at f\[(\d+)\].*rho=([-\d.e+]+) fill=([-\d.e+]+) \| obstacle geometry: pinched opposite pairs=(\d+)/13, obstacle neighbours=(\d+)/26 \| context: former-moving-obstacle neighbours=(\d+), interface neighbours=(\d+), gas neighbours=(\d+), fastest fluid neighbour \|u\|=([\d.e+-]+), obstacle owners \(uids\)=\[ ([\d ]*)\]")
hv = re.compile(r"\[HIGH-VELOCITY-PRECURSOR t=(\d+) after=(\w+)\] (\w+\s*) cell .*world=\(([\d.,]+)\) \|u\|=([\d.]+)")
gl = re.compile(r"\[(\d+)\] \[MEM-GARBAGE-LINK\] particle uid=(\d+) fluid cell block-local \(([\d,]+)\) flags=\[([\w| ]+)\] fill=([-\d.e+]+) density=([-\d.e+]+)")
stale = re.compile(r"\[\s*(\d+)\] \[STALE-PARTICLE-IDX\] rank=(\d+) .*uid=(\d+)")
pops, hvs, gls, stales = [], [], [], []
for line in open(sys.argv[1], errors="replace"):
    if "POP-PRECURSOR" in line:
        m = pop.search(line); pops.append(m.groups()) if m else None
    elif "HIGH-VELOCITY-PRECURSOR" in line:
        m = hv.search(line); hvs.append(m.groups()) if m else None
    elif "MEM-GARBAGE-LINK" in line:
        m = gl.search(line); gls.append(m.groups()) if m else None
    elif "STALE-PARTICLE-IDX" in line:
        m = stale.search(line); stales.append(m.groups()) if m else None
print(f"POP-PRECURSOR lines: {len(pops)}, HIGH-VELOCITY-PRECURSOR: {len(hvs)}, MEM-GARBAGE-LINK: {len(gls)}, STALE-PARTICLE-IDX: {len(stales)}")
for p in sorted(pops, key=lambda p: int(p[0]))[:12]:
    print(f"  POP  t={p[0]} after={p[1]} {p[2].strip()} world=({p[4]}) max|f|={p[5]} f[{p[6]}] rho={p[7]} fill={p[8]}  pinched pairs={p[9]}/13, obstacle neighbours={p[10]}/26, fMO neighbours={p[11]}, interface nbrs={p[12]}, gas nbrs={p[13]}, fastest nbr |u|={p[14]}, owners=[{p[15].strip()}]")
if pops:
    pinched = sum(1 for p in pops if int(p[9]) >= 1); print(f"  flagged cells that are pinched (>=1 opposite obstacle pair): {pinched}/{len(pops)}; obstacle neighbours histogram: {dict(sorted(collections.Counter(int(p[10]) for p in pops).items()))}")
for h in sorted(hvs, key=lambda h: int(h[0]))[:6]: print(f"  HV   t={h[0]} after={h[1]} {h[2].strip()} world=({h[3]}) |u|={h[4]}")
for g in gls[:4]: print(f"  LINK step={g[0]} uid={g[1]} cell={g[2]} flags=[{g[3]}] fill={g[4]} density={g[5]}")
for s in stales[:2]: print(f"  STALE step={s[0]} rank={s[1]} uid={s[2]}")
