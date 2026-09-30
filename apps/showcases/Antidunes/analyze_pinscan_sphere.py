import re, numpy as np
ROOT = "/home/dhuber/University/Thesis/animations/pinscan_sphere_20260928"
G = 5.06585e-5; ZFS_NOM = 59.8831
def cap(zc, R, zfs):
    h = np.clip(zfs - (zc - R), 0, 2*R); return np.pi*h*h*(3*R - h)/3
import os
for sfx, name in (("", "baseline (simple bounce-back)"), ("_ibb", "interpolated bounce-back (CLI)"), ("_nofw", "useFillWeighting false (R=24 only)"), ("_sig0", "surface tension 0 (real g)"), ("_sig2x", "surface tension doubled (R=24 only)")):
    print(f"== {name}")
    print("R   | n  | mean F_hydro_z  std    | z_c      waterline | V_sub(wl)  ratio    excess  | extra vol/(pi R^2) [cells] | ratio @nominal wl")
    for R in (12, 24, 36):
        base = f"{ROOT}/R{R}{sfx}/R{R}_real_f075{sfx}"
        if not os.path.exists(base + ".forces.txt"): continue
        fz, zc = [], None
        for l in open(base + ".forces.txt"):
            m = re.search(r"step=(\d+).*?F_hydro=\(([-\d.e+]+),([-\d.e+]+),([-\d.e+]+)\).*?pos=<[^,]+,[^,]+,([-\d.e+]+)>", l)
            if m and int(m.group(1)) >= 6000: fz.append(float(m.group(4))); zc = float(m.group(5))
        fz = np.array(fz)
        wl = float([l for l in open(base + ".column.txt") if l.startswith("waterline")][0].split()[1])
        V = cap(zc, R, wl); ratio = fz.mean()/(G*V); extra = fz.mean()/G - V
        print(f"{R:3d} | {len(fz):2d} | {fz.mean():.5f}  {fz.std():.5f} | {zc:.3f}  {wl:.3f} | {V:9.0f}  {ratio:.4f}  {ratio-1:+.4f} | {extra/(np.pi*R*R):+.3f} | {fz.mean()/(G*cap(zc, R, ZFS_NOM)):.4f}")
