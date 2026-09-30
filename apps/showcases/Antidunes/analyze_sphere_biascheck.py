import sqlite3, numpy as np, pandas as pd, glob, sys
ROOT = "/home/dhuber/University/Thesis/animations/sphere-r3-biascheck-density-column_20260928"
R = 24.0; A = 96*135
OLD = {  # artifact tail-10k means (real, bogner) and artifact "target"
 "d025": dict(rho=0.25, old=(67.5228, 67.4640), target=67.26),
 "d050": dict(rho=0.50, old=(60.3748, 60.3977), target=60.17),
 "d075": dict(rho=0.75, old=(53.3056, 53.4827), target=52.92)}
G = {"real": 5.06585e-5, "bogner": 7.5e-4}

def arms(tag):
    con = sqlite3.connect(f"{ROOT}/{tag}/eval/antidunes.sqlite")
    df = pd.read_sql("select timestamp, timestep, maxParticlePos_z z, totalFluidMass m from runs order by timestamp, timestep", con)
    # split arms where timestep resets
    brk = np.where(np.diff(df.timestep.values) < 0)[0]
    idx = [0] + list(brk+1) + [len(df)]
    parts = [df.iloc[idx[i]:idx[i+1]] for i in range(len(idx)-1)]
    return parts

def cap_int(zc, zfs, rho_fun, n=4000):
    """integral over sphere part below zfs of rho(z) dV"""
    z = np.linspace(max(zc-R, 0), min(zc+R, zfs), n)
    area = np.pi*np.maximum(R**2-(z-zc)**2, 0)
    return np.trapezoid(area*rho_fun(z), z)

def solve_zc(zfs, rho_fun, ratio):
    V = 4/3*np.pi*R**3
    lo, hi = zfs-R, zfs+R
    for _ in range(80):
        mid = 0.5*(lo+hi)
        # higher zc -> less submerged -> smaller integral
        if cap_int(mid, zfs, rho_fun) > ratio*V: lo = mid
        else: hi = mid
    return 0.5*(lo+hi)

for tag, d in OLD.items():
    parts = arms(tag)
    print(f"\n===== {tag} (densityRatio {d['rho']}) : {len(parts)} arms in sqlite =====")
    zs = []
    for k, p in enumerate(parts):
        tail = p[p.timestep >= p.timestep.max()-10000]
        zs.append(tail.z.mean())
        print(f" arm{k}: steps {p.timestep.min()}-{p.timestep.max()}  tail-10k z mean={tail.z.mean():.4f}  (old artifact: {d['old'][k]:.4f})  mass tail={tail.m.mean():.0f}")
    # density column (only last arm = bogner)
    dfs = [pd.read_csv(f, usecols=["mainStep","worldZ","flag","density","fill"]) for f in sorted(glob.glob(f"{ROOT}/{tag}/density-column/rank_*.csv")) if len(open(f).read(200).splitlines())>1]
    c = pd.concat(dfs)
    smax = c.mainStep.max(); print(f" density-column steps {c.mainStep.min()}..{smax}, rows {len(c)}")
    tailc = c[c.mainStep >= smax-10000]
    prof = tailc.groupby("worldZ").agg(rho=("density","mean"), fill=("fill","mean"), n=("density","size")).reset_index()
    zfs_t = tailc.groupby("mainStep").fill.sum()          # free-surface height per step (sum of fill over column)
    zfs = zfs_t.mean(); print(f" far-field waterline (sum of fill, tail-10k mean) = {zfs:.4f} (std {zfs_t.std():.4f})")
    liq = prof[(prof.fill>0.999)]
    a, b = np.polyfit(liq.worldZ, liq.rho, 1)
    g = G["bogner"]
    print(f" rho(z) linear fit in liquid: slope={a:.4e} /cell  (theory -3g = {-3*g:.4e}), rho(z=0)={b:.5f}, rho at waterline={a*zfs+b:.5f}")
    rho_fun = lambda z: np.exp(-3*g*(z - zfs))  # isothermal LBM hydrostatics, rho=1 at surface
    naive = lambda z: np.ones_like(z)
    V = 4/3*np.pi*R**3
    zc_comp = solve_zc(zfs, rho_fun, d['rho']); zc_naive = solve_zc(zfs, naive, d['rho'])
    print(f" predicted centre, incompressible Archimedes @ measured waterline : {zc_naive:.4f}")
    print(f" predicted centre, compressible (measured rho(z))                  : {zc_comp:.4f}")
    print(f" MEASURED centre (bogner arm tail-10k)                             : {zs[-1]:.4f}")
    print(f"   bias vs naive = {zs[-1]-zc_naive:+.4f} LU ;  bias vs compressible = {zs[-1]-zc_comp:+.4f} LU ; compressibility explains {(zc_comp-zc_naive):+.4f} LU")
    # mass check: mean-density explanation of the 6% mass difference between arms
    if len(parts)==2:
        mr = parts[1][parts[1].timestep>=parts[1].timestep.max()-10000].m.mean()/parts[0][parts[0].timestep>=parts[0].timestep.max()-10000].m.mean()
        pred = (1+3*G['bogner']*zfs/2)/(1+3*G['real']*zfs/2)
        print(f" mass ratio bogner/real = {mr:.4f} ; predicted from hydrostatic mean density = {pred:.4f}")
    # real-gravity arm: predict with theoretical slope, same waterline (measured position only)
    rho_real = lambda z: np.exp(-3*G["real"]*(z - zfs))
    zc_real = solve_zc(zfs, rho_real, d['rho'])
    print(f" [real-gravity arm, theoretical slope -3g, same waterline] predicted centre {zc_real:.4f}  vs measured {zs[0]:.4f}  (naive {zc_naive:.4f}) -> bias vs compressible {zs[0]-zc_real:+.4f}, vs naive {zs[0]-zc_naive:+.4f}")
