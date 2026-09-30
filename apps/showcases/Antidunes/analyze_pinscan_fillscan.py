"""Fill scan (THESIS_PLAN.md 2026-09-28): pinned R=24, 0.75 submerged, real gravity, liquid level shifted so the fill of the
top interface cell (phi = frac of the sum-fill waterline) varies. Prediction of the interface-cell-centre pressure-condition
hypothesis: rho=1 sits (0.5-phi) above the fill waterline, excess extra-volume/(pi R^2) = (r_w^2/R^2) * (0.5-phi)."""
import re, os, numpy as np
ROOT = "/home/dhuber/University/Thesis/animations/pinscan_sphere_20260928"; G = 5.06585e-5; R = 24
def cap(zc, zfs):
    h = np.clip(zfs - (zc - R), 0, 2*R); return np.pi*h*h*(3*R - h)/3
cases = [("baseline (nominal)", "R24", "R24_real_f075")] + [(f"level {lab}", f"R24_{t}", f"R24_real_f075_{t}") for lab, t in (("-0.175", "lvm175"), ("+0.225", "lvp225"), ("+0.625", "lvp625"))]
print("case              | wl(sum fill) phi_top | z0-wl (rho=1 fit) | 0.5-phi | mean F_z   ratio    excess | extra vol/(pi R^2) | predicted (rw^2/R^2)*(z0-wl) | predicted (rw^2/R^2)*(0.5-phi)")
for name, d, base in cases:
    f = f"{ROOT}/{d}/{base}"
    if not os.path.exists(f + ".forces.txt"): continue
    fz, zc = [], None
    for l in open(f + ".forces.txt"):
        m = re.search(r"step=(\d+).*?F_hydro=\(([-\d.e+]+),([-\d.e+]+),([-\d.e+]+)\).*?pos=<[^,]+,[^,]+,([-\d.e+]+)>", l)
        if m and int(m.group(1)) >= 6000: fz.append(float(m.group(4))); zc = float(m.group(5))
    fz = np.array(fz); L = open(f + ".column.txt").read().splitlines()
    wl = float([l for l in L if l.startswith("waterline")][0].split()[1]); phi = wl - np.floor(wl)
    p = np.array([[float(x) for x in l.split()[1:3]] for l in L if l.startswith("rho")]); p = p[np.argsort(p[:, 0])]
    m = (p[:, 0] < wl - 4) & (p[:, 0] > 5); a, b = np.polyfit(p[m, 0], p[m, 1], 1); off = (1 - b)/a - wl
    V = cap(zc, wl); extra = (fz.mean()/G - V)/(np.pi*R*R); rw2 = R**2 - (zc - wl)**2
    print(f"{name:17s} | {wl:7.3f}   {phi:5.3f}    | {off:+.3f}            | {0.5-phi:+.3f}  | {fz.mean():.5f}  {fz.mean()/(G*V):.4f}  {fz.mean()/(G*V)-1:+.4f} | {extra:+.3f}            | {rw2/R**2*off:+.3f}                        | {rw2/R**2*(0.5-phi):+.3f}")
