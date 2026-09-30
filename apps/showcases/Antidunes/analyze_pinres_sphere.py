"""Pinned-sphere RESOLUTION scan (R1/R2/R4 family, real g_SI): excess buoyancy and pressure-zero offset in cells vs resolution."""
import re, os, numpy as np
ROOT = "/home/dhuber/University/Thesis/animations/pinscan_sphere_20260928"
G = {1: 2.137150014e-05, 2: 2.671437518e-06, 4: 3.339296897e-07}; STEPS = {1: 12000, 2: 24000, 4: 48000}
def cap(zc, r, zfs):
    h = np.clip(zfs - (zc - r), 0, 2*r); return np.pi*h*h*(3*r - h)/3
print("R | r  | g_LU      | n  | mean F_z        std(%)  | z_c      wl(fill)  | ratio    excess   | extra vol/(pi r^2) [cells] | offset z0-wl (exp fit) [cells] | offset/r")
for R in (1, 2, 4):
    base = f"{ROOT}/pinres_R{R}/R{R}_pinres"
    if not os.path.exists(base + ".forces.txt"): continue
    r = 8.0*R; fz, zc = [], None
    for l in open(base + ".forces.txt"):
        m = re.search(r"step=(\d+).*?F_hydro=\(([-\d.e+]+),([-\d.e+]+),([-\d.e+]+)\).*?pos=<[^,]+,[^,]+,([-\d.e+]+)>", l)
        if m and int(m.group(1)) >= STEPS[R]*2//3: fz.append(float(m.group(4))); zc = float(m.group(5))
    fz = np.array(fz); L = open(base + ".column.txt").read().splitlines()
    wl = float([l for l in L if l.startswith("waterline")][0].split()[1])
    p = np.array([[float(x) for x in l.split()[1:3]] for l in L if l.startswith("rho")]); p = p[np.argsort(p[:, 0])]
    m = (p[:, 0] > 5) & (p[:, 0] < wl - 4); s, c = np.polyfit(p[m, 0], np.log(p[m, 1]), 1); off = -c/s - wl
    V = cap(zc, r, wl); ratio = fz.mean()/(G[R]*V); extra = (fz.mean()/G[R] - V)/(np.pi*r*r)
    print(f"{R} | {r:3.0f} | {G[R]:.3e} | {len(fz):2d} | {fz.mean():.6e}  {100*fz.std()/fz.mean():5.2f} | {zc:8.3f} {wl:8.3f} | {ratio:.4f}  {ratio-1:+.4f} | {extra:+.3f}                    | {off:+.3f}                        | {off/r:+.4f}")
