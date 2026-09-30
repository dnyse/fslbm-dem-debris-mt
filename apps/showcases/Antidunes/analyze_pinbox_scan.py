"""Pinned-box heel scan at real gravity: measured MEM force / torque_y vs ideal hydrostatics at the measured waterline, with and without the
level raised by the measured pressure-zero offset (exponential rho(z) fit).  Results: ~/University/Thesis/animations/pinscan_sphere_20260928/pinbox/"""
import re, glob, os, sys, numpy as np, math
ROOT = "/home/dhuber/University/Thesis/animations/pinscan_sphere_20260928/" + (sys.argv[1] if len(sys.argv) > 1 else "pinbox"); G = 5.06585e-5; S = 30.0; LY = 60.0
def square(zc, a):
    a = math.radians(a); c, sn = math.cos(a), math.sin(a); h = S/2
    return [(x*c + z*sn, -x*sn + z*c + zc) for x, z in ((-h, -h), (h, -h), (h, h), (-h, h))]
def clip_below(poly, zw):
    out = []
    for i in range(len(poly)):
        p, q = poly[i], poly[(i+1) % len(poly)]; pin, qin = p[1] <= zw, q[1] <= zw
        if pin: out.append(p)
        if pin != qin:
            t = (zw - p[1])/(q[1] - p[1]); out.append((p[0] + t*(q[0] - p[0]), zw))
    return out
def area_centroid(poly):
    if len(poly) < 3: return 0.0, 0.0, 0.0
    A = cx = cz = 0.0
    for i in range(len(poly)):
        x0, z0 = poly[i]; x1, z1 = poly[(i+1) % len(poly)]; c = x0*z1 - x1*z0
        A += c; cx += (x0 + x1)*c; cz += (z0 + z1)*c
    A *= 0.5; return A, cx/(6*A), cz/(6*A)
def ideal(zc, a, level, z0):
    A, xb, zb = area_centroid(clip_below(square(zc, a), level))
    rho = 1 + 3*G*(z0 - zb)                              # hydrostatic density at the centroid (compressible, rho=1 at z0)
    F = G*rho*A*LY; return F, -F*xb                      # torque_y = -F*(x_B - x_G)
rows = []
for d in sorted(glob.glob(ROOT + "/box_*")):
    lab = os.path.basename(d); base = f"{d}/{lab}"
    if not os.path.exists(base + ".forces.txt"): continue
    Fz, Ty, pz, tilt = [], [], None, []
    for l in open(base + ".forces.txt"):
        m = re.search(r"\[BoxForces\].*?step=(\d+).*?F_hydro=\(([-\d.e+]+),([-\d.e+]+),([-\d.e+]+)\).*?pos=<[^,]+,[^,]+,([-\d.e+]+)>.*?tiltY=([-\d.e+]+)", l)
        if m and int(m.group(1)) >= 6000: Fz.append(float(m.group(4))); pz = float(m.group(5)); tilt.append(float(m.group(6)))
        m = re.search(r"\[BoxTorque\].*?step=(\d+).*?hydroTorque=\(([-\d.e+]+),([-\d.e+]+),([-\d.e+]+)\)", l)
        if m and int(m.group(1)) >= 6000: Ty.append(float(m.group(3)))
    Cl = open(base + ".column.txt").read().splitlines()
    wl = float([l for l in Cl if l.startswith("waterline")][0].split()[1])
    p = np.array([[float(x) for x in l.split()[1:3]] for l in Cl if l.startswith("rho")]); p = p[np.argsort(p[:, 0])]
    m = (p[:, 0] > 5) & (p[:, 0] < wl - 4); s, c = np.polyfit(p[m, 0], np.log(p[m, 1]), 1); z0 = -c/s; delta = z0 - wl
    a = float(np.mean(tilt)); F, T = np.mean(Fz), np.mean(Ty)
    F0, T0 = ideal(pz, a, wl, z0); F1, T1 = ideal(pz, a, z0, z0)
    rows.append((lab, a, pz, wl, delta, F, np.std(Fz), F/F0, F/F1, T, np.std(Ty), T0, T1))
print("case             | heel   z_c     wl     offset | F_meas    F/F_ideal(wl)  F/F_ideal(wl+off) | tau_meas   std     tau_ideal(wl)  tau_ideal(wl+off) | tau_meas-ideal(wl)  tau_meas-ideal(wl+off)")
for r in rows:
    print(f"{r[0]:16s} | {r[1]:6.2f} {r[2]:7.3f} {r[3]:7.3f} {r[4]:+.3f} | {r[5]:8.4f}  {r[7]:.4f}         {r[8]:.4f}           | {r[9]:+9.4f} {r[10]:6.4f} {r[11]:+9.4f}     {r[12]:+9.4f}          | {r[9]-r[11]:+8.4f}            {r[9]-r[12]:+8.4f}")
