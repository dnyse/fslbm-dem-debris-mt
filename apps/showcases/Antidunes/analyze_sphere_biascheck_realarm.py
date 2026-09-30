import numpy as np, sqlite3, pandas as pd
exec(open("analyze_sphere_biascheck.py").read().split("for tag, d in OLD.items():")[0])
ART_WL = {"d025":58.92, "d050":60.17, "d075":61.26}   # artifact's (real-arm, VTK-measured) waterlines
for tag, d in OLD.items():
    parts = arms(tag); p = parts[0]
    t = p[p.timestep>=p.timestep.max()-10000]
    # fit z(t)=z_inf + A exp(-t/tau) cos(w t+phi) over the last 25k steps by scanning w
    q = p[p.timestep>=25000]; x=q.timestep.values.astype(float); y=q.z.values
    best=None
    for w in np.linspace(2*np.pi/12000, 2*np.pi/2000, 400):
        for tau in (5e3,1e4,2e4,4e4,1e9):
            M = np.c_[np.ones_like(x), np.exp(-(x-x[0])/tau)*np.cos(w*x), np.exp(-(x-x[0])/tau)*np.sin(w*x)]
            c,res,*_ = np.linalg.lstsq(M,y,rcond=None); r=np.sum((M@c-y)**2)
            if best is None or r<best[0]: best=(r,w,tau,c)
    r,w,tau,c = best
    zinf=c[0]; amp=np.hypot(c[1],c[2])
    zfs = ART_WL[tag]
    naive = solve_zc(zfs, lambda z: np.ones_like(z), d['rho'])
    comp  = solve_zc(zfs, lambda z: np.exp(-3*G['real']*(z-zfs)), d['rho'])
    print(f"{tag}: tail-10k range {t.z.min():.3f}..{t.z.max():.3f} (p-p {t.z.max()-t.z.min():.3f}) mean {t.z.mean():.4f} | damped-cos fit: z_inf={zinf:.4f}, period={2*np.pi/w:.0f} steps, tau={tau:.0f}, amp@25k={amp:.3f}")
    print(f"      predicted @ artifact waterline {zfs}: naive {naive:.4f}, compressible {comp:.4f} -> bias of z_inf vs naive {zinf-naive:+.3f}")
