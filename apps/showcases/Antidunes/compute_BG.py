"""Compute (approximate) center of buoyancy B and center of gravity G for the Bogner-like free box, from
the logged pose (position, tiltY) and the measured ambient water level, and compare the resulting analytic
hydrostatic torque to the ACTUALLY measured MEM torque_y. G = box position (homogeneous body, comOffsetZ=0).
B = centroid of the box's square cross-section (side s, in its local X-Z plane) clipped below the ambient
level, rotated by tiltY, expressed in world (x,z). Approximation: assumes the ambient level is uniform along
the box's own Y-extent (no slosh/tilt-induced variation along Y) -- a simplification, not exact.
"""
import re, sys, collections, statistics, math
import numpy as np

def load_forces(path):
    S=[];T=[];X=[];Y=[];Z=[];Fz=[];Fx=[]
    Ts=[];Ty=[]
    for l in open(path):
        m = re.search(r"\[BoxForces\].*?step=(\d+).*?F_hydro=\(([-\d.e+]+),([-\d.e+]+),([-\d.e+]+)\).*?pos=<([-\d.e+]+),([-\d.e+]+),([-\d.e+]+)>.*?tiltY=([-\d.e+]+)", l)
        if m:
            S.append(int(m.group(1))); Fx.append(float(m.group(2))); Fz.append(float(m.group(4)))
            X.append(float(m.group(5))); Y.append(float(m.group(6))); Z.append(float(m.group(7))); T.append(float(m.group(8)))
        m2 = re.search(r"\[BoxTorque\].*?step=(\d+).*?hydroTorque=\(([-\d.e+]+),([-\d.e+]+),([-\d.e+]+)\)", l)
        if m2: Ts.append(int(m2.group(1))); Ty.append(float(m2.group(3)))
    return (np.array(S), np.array(T), np.array(X), np.array(Y), np.array(Z), np.array(Fx), np.array(Fz)), (np.array(Ts), np.array(Ty))

def load_level(path):
    h = collections.defaultdict(list)
    for l in open(path):
        p = l.split()
        if len(p) < 5: continue
        step, gx, gy, ht, n = int(p[0]), int(p[1]), int(p[2]), float(p[3]), int(p[4])
        if n == 0: h[step].append(ht)
    st = sorted(h)
    return st, {s: statistics.median(h[s]) for s in st}

def square(zc, xc, alpha, s):
    a = math.radians(alpha); c, sn = math.cos(a), math.sin(a); hh = s/2
    # local (x,z) about (xc,zc); returns world (x,z) corners
    return [(xc + x*c + z*sn, zc - x*sn + z*c) for x, z in ((-hh,-hh),(hh,-hh),(hh,hh),(-hh,hh))]

def clip_below(poly, zw):
    out = []
    for i in range(len(poly)):
        p, q = poly[i], poly[(i+1) % len(poly)]; pin, qin = p[1] <= zw, q[1] <= zw
        if pin: out.append(p)
        if pin != qin:
            t = (zw - p[1]) / (q[1] - p[1]); out.append((p[0] + t*(q[0]-p[0]), zw))
    return out

def area_centroid(poly):
    if len(poly) < 3: return 0.0, 0.0, 0.0
    A = 0.0; Cx = 0.0; Cz = 0.0
    n = len(poly)
    for i in range(n):
        x0, z0 = poly[i]; x1, z1 = poly[(i+1) % n]
        cross = x0*z1 - x1*z0
        A += cross; Cx += (x0+x1)*cross; Cz += (z0+z1)*cross
    A *= 0.5
    if abs(A) < 1e-12: return 0.0, 0.0, 0.0
    Cx /= (6*A); Cz /= (6*A)
    return A, Cx, Cz

def analyze(label, forces_path, wh_path, s, g_lu, y_len, nprint=25):
    (S,T,X,Y,Z,Fx,Fz), (Ts,Ty) = load_forces(forces_path)
    st, lvl = load_level(wh_path)
    print(f"\n===== {label}  (s={s}, g_LU={g_lu:.4e})")
    print(f"{'step':>7} {'tiltY':>7} {'x_G':>7} {'z_G':>7} {'level':>7} {'x_B':>7} {'z_B':>7} {'lever(xB-xG)':>13} {'V_sub':>8} {'F_ideal':>8} {'tau_ideal':>10} {'F_meas':>8} {'tau_meas':>10}")
    idxs = np.linspace(0, len(S)-1, nprint).astype(int)
    for i in idxs:
        step = S[i]
        j = min(st, key=lambda x: abs(x-step))
        level = lvl[j]
        poly = square(Z[i], X[i], T[i], s)
        sub = clip_below(poly, level)
        A, cx, cz = area_centroid(sub)
        lever = cx - X[i]
        Fideal = g_lu * A * y_len
        tau_ideal = -Fideal * lever
        k = np.argmin(np.abs(Ts - step))
        print(f"{step:7d} {T[i]:7.2f} {X[i]:7.2f} {Z[i]:7.2f} {level:7.2f} {cx:7.2f} {cz:7.2f} {lever:13.3f} {A:8.1f} {Fideal:8.4f} {tau_ideal:10.4f} {Fz[i]:8.4f} {Ty[k]:10.4f}")

if __name__ == "__main__":
    G_LU = 5.06585e-5
    analyze("rho0.75 orig 150k (4273396)",
             "/home/dhuber/University/Thesis/animations/bogner3_deep_fixesoff_20260928/4273396/b3_r075.forces.txt",
             "/home/dhuber/University/Thesis/animations/bogner3_deep_fixesoff_20260928/4273396/b3_r075.waterheight.txt",
             s=18, g_lu=G_LU, y_len=36)
    analyze("s=16 normal seed (4273515)",
             "/home/dhuber/University/Thesis/animations/bogner3_deep_s16_20260929/4273515/b3_r075_s16.forces.txt",
             "/home/dhuber/University/Thesis/animations/bogner3_deep_s16_20260929/4273515/b3_r075_s16.waterheight.txt",
             s=16, g_lu=G_LU, y_len=32, nprint=15)
    analyze("s=16 EQUILIBRIUM seed (4273516)",
             "/home/dhuber/University/Thesis/animations/bogner3_deep_s16_20260929/4273516/b3_r075_s16_eqseed.forces.txt",
             "/home/dhuber/University/Thesis/animations/bogner3_deep_s16_20260929/4273516/b3_r075_s16_eqseed.waterheight.txt",
             s=16, g_lu=G_LU, y_len=32, nprint=20)
    analyze("s=32 FIXED (4273525)",
             "/home/dhuber/University/Thesis/animations/bogner3_deep_s32_20260929/4273525/b3_r075_s32.forces.txt",
             "/home/dhuber/University/Thesis/animations/bogner3_deep_s32_20260929/4273525/b3_r075_s32.waterheight.txt",
             s=32, g_lu=G_LU, y_len=64, nprint=20)
    analyze("s=32 CONTINUOUS 100k (4273570)",
             "/home/dhuber/University/Thesis/animations/bogner3_deep_s32_20260929/4273570/b3_r075_s32_100k.forces.txt",
             "/home/dhuber/University/Thesis/animations/bogner3_deep_s32_20260929/4273570/b3_r075_s32_100k.waterheight.txt",
             s=32, g_lu=G_LU, y_len=64, nprint=12)
