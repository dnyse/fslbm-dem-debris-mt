import numpy as np, pandas as pd, glob
ROOT="/home/dhuber/University/Thesis/animations/"
def fit(z, rho, wl, g):
    m=(z>5)&(z<wl-4)
    a,b=np.polyfit(z[m],rho[m],1); z0_lin=(1-b)/a                       # straight line (what I used before)
    s,c=np.polyfit(z[m],np.log(rho[m]),1); z0_log=-c/s                  # exponential: ln rho = c + s z
    return z0_lin-wl, z0_log-wl, s/(-3*g)
print("== Bogner-gravity R3 free-floating columns (g=7.5e-4)")
for tag in ("d025","d050","d075"):
    fs=sorted(glob.glob(ROOT+f"sphere-r3-biascheck-density-column_20260928/{tag}/density-column/rank_[0-3].csv"))
    c=pd.concat([pd.read_csv(f,usecols=["mainStep","worldZ","density","fill"]) for f in fs]); c=c[c.mainStep>=c.mainStep.max()-10000]
    wl=c.groupby("mainStep").fill.sum().mean(); p=c.groupby("worldZ").density.mean()
    l,e,r=fit(p.index.values,p.values,wl,7.5e-4)
    print(f" {tag}: waterline {wl:.3f} | offset z0-wl: straight-line fit {l:+.3f}  exponential fit {e:+.3f} cell | slope/(-3g) {r:.3f}")
print("== real-gravity pinned columns (g=5.06585e-5)")
for f in sorted(glob.glob(ROOT+"pinscan_sphere_20260928/*/*.column.txt")):
    L=open(f).read().splitlines(); wl=float([l for l in L if l.startswith("waterline")][0].split()[1])
    p=np.array([[float(x) for x in l.split()[1:3]] for l in L if l.startswith("rho")]); p=p[np.argsort(p[:,0])]
    l,e,r=fit(p[:,0],p[:,1],wl,5.06585e-5); print(f" {f.split('/')[-1][:-11]:24s} linear {l:+.3f}  exponential {e:+.3f}")
