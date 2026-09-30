import numpy as np, glob, re
ROOT="/home/dhuber/University/Thesis/animations/pinscan_sphere_20260928"; G=5.06585e-5
print("case                      | waterline(sum fill) | z where fitted rho=1 | offset (z0-wl) [cells] | slope/(-3g) | rho at wl")
for f in sorted(glob.glob(ROOT+"/*/*.column.txt")):
    L=open(f).read().splitlines(); wl=float([l for l in L if l.startswith("waterline")][0].split()[1])
    p=np.array([[float(x) for x in l.split()[1:3]] for l in L if l.startswith("rho")]); p=p[np.argsort(p[:,0])]
    m=(p[:,0]<wl-4)&(p[:,0]>5)       # liquid interior, away from floor and interface
    a,b=np.polyfit(p[m,0],p[m,1],1); z0=(1-b)/a
    print(f"{f.split('/')[-1][:-11]:26s}| {wl:8.3f}            | {z0:8.3f}             | {z0-wl:+.3f}                 | {a/(-3*G):.3f}       | {a*wl+b:.6f}")
