"""Bogner-like free box at real gravity: heel-angle history from [BoxForces] (tiltY, pos).  Results in
~/University/Thesis/animations/pinscan_sphere_20260928/bognerlike/<label>/<label>.forces.txt"""
import re, glob, os, sys, numpy as np
ROOT = "/home/dhuber/University/Thesis/animations/pinscan_sphere_20260928/" + (sys.argv[1] if len(sys.argv) > 1 else "bognerlike")
print("case             | n   | tilt at step: 2.5k   10k    25k    50k    75k    end  | last-20k mean  std   min..max  | x range (cells)  y range | target")
for d in sorted(glob.glob(ROOT + "/bl_*")):
    lab = os.path.basename(d); f = f"{d}/{lab}.forces.txt"
    if not os.path.exists(f): continue
    st, tl, xs, ys, zs = [], [], [], [], []
    for l in open(f):
        m = re.search(r"\[BoxForces\].*?step=(\d+).*?pos=<([-\d.e+]+),([-\d.e+]+),([-\d.e+]+)>.*?tiltY=([-\d.e+]+)", l)
        if m: st.append(int(m.group(1))); xs.append(float(m.group(2))); ys.append(float(m.group(3))); zs.append(float(m.group(4))); tl.append(float(m.group(5)))
    st, tl, xs, ys = map(np.array, (st, tl, xs, ys)); rho = int(lab.split("_")[1][1:])/100
    at = lambda s: tl[np.argmin(abs(st - s))]
    end = st.max(); m = st >= end - 20000
    tgt = 45.0 if abs(rho - 0.5) < 1e-9 else 26.565
    print(f"{lab:16s} | {len(st):3d} | {at(2500):6.2f} {at(10000):6.2f} {at(25000):6.2f} {at(50000):6.2f} {at(75000):6.2f} {tl[-1]:6.2f} | {tl[m].mean():6.2f} {tl[m].std():5.2f} {tl[m].min():6.2f}..{tl[m].max():6.2f} | {xs.min():5.1f}..{xs.max():5.1f}    {ys.min():5.1f}..{ys.max():5.1f} | {tgt}")
