"""Fully submerged pinned spheres (no waterline): ratio of mean MEM F_hydro_z to the compressible Archimedes
force g * integral(rho(z) dV) over the sphere, rho(z) = measured far-field profile (column.txt)."""
import re, numpy as np
ROOT = "/home/dhuber/University/Thesis/animations/pinscan_sphere_20260928"; G = 5.06585e-5
print("R  | n | mean F_hydro_z   std     | z_c     | g*V (rho=1) | g*int(rho dV) | ratio vs rho=1 | ratio vs compressible (excess) | extra vol/(pi R^2) [cells]")
for R in (12, 24):
    base = f"{ROOT}/R{R}_f100/R{R}_real_f100"
    fz, zc = [], None
    for l in open(base + ".forces.txt"):
        m = re.search(r"step=(\d+).*?F_hydro=\(([-\d.e+]+),([-\d.e+]+),([-\d.e+]+)\).*?pos=<[^,]+,[^,]+,([-\d.e+]+)>", l)
        if m and int(m.group(1)) >= 6000: fz.append(float(m.group(4))); zc = float(m.group(5))
    fz = np.array(fz)
    prof = np.array([[float(x) for x in l.split()[1:3]] for l in open(base + ".column.txt") if l.startswith("rho")]); prof = prof[np.argsort(prof[:, 0])]
    z = np.linspace(zc - R, zc + R, 4001); area = np.pi*np.maximum(R**2 - (z - zc)**2, 0)
    rho = np.interp(z, prof[:, 0], prof[:, 1])
    V = 4/3*np.pi*R**3; comp = G*np.trapezoid(area*rho, z)
    print(f"{R:2d} | {len(fz)} | {fz.mean():.5f}  {fz.std():.5f} | {zc:.3f} | {G*V:.5f} | {comp:.5f} | {fz.mean()/(G*V):.4f} | {fz.mean()/comp:.4f} ({fz.mean()/comp-1:+.2%}) | {(fz.mean()/G - comp/G)/(np.pi*R*R):+.3f}")
