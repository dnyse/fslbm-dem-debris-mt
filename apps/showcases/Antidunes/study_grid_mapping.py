# [Claude, 2026-09-21] Numerical study of the cell-centre particle mapping (MovingParticleMapping: cell = obstacle iff its centre is inside the body): sealed fluid cells between touching bodies and mapped-volume fidelity. Run: python3 study_grid_mapping.py (numpy only, ~2 min). Results: THESIS_PLAN.md 2026-09-21 "GRID MAPPING DEEP DIVE".
import numpy as np, math, sys
rng=np.random.default_rng(1)
N=44
ax=np.arange(N)+0.5
X,Y,Z=np.meshgrid(ax,ax,ax,indexing='ij')
P=np.stack([X,Y,Z],-1)

def rot_random():
    q=rng.normal(size=4); q/=np.linalg.norm(q)
    a,b,c,d=q
    return np.array([[a*a+b*b-c*c-d*d,2*(b*c-a*d),2*(b*d+a*c)],[2*(b*c+a*d),a*a-b*b+c*c-d*d,2*(c*d-a*b)],[2*(b*d-a*c),2*(c*d+a*b),a*a-b*b-c*c+d*d]])

def sphere(c,r): return ((P-c)**2).sum(-1)<=r*r
def cyl(c,R,r,L):   # axis = R[:,2]
    d=P-c; u=d@R  # body-frame coords
    return (u[...,0]**2+u[...,1]**2<=r*r)&(np.abs(u[...,2])<=L/2)
def box(c,R,h):
    u=(P-c)@R
    return (np.abs(u[...,0])<=h[0])&(np.abs(u[...,1])<=h[1])&(np.abs(u[...,2])<=h[2])

def obst_nb(mask):
    m=np.pad(mask.astype(np.int8),1)
    n=np.zeros(mask.shape,np.int8)
    for dx in (-1,0,1):
        for dy in (-1,0,1):
            for dz in (-1,0,1):
                if dx==0 and dy==0 and dz==0: continue
                n+=m[1+dx:1+dx+mask.shape[0],1+dy:1+dy+mask.shape[1],1+dz:1+dz+mask.shape[2]]
    return n
def sealed(mask,thr):
    n=obst_nb(mask); fluid=~mask
    # ignore boundary layer of grid
    s=fluid&(n>=thr); s[:3]=s[-3:]=False; s[:,:3]=s[:,-3:]=False; s[:,:,:3]=s[:,:,-3:]=False
    return int(s.sum())

def rand_offset(): return rng.uniform(0,1,3)
C0=np.array([N/2,N/2,N/2])

def pair_axes(shape, gap, R, r=3.0, L=18.0):
    # two parallel bodies, centres separated perpendicular to the axis by 2r+gap in a random perpendicular direction
    ax_=R[:,2]; perp=R[:,0]*math.cos(rng.uniform(0,2*math.pi))+R[:,1]*math.sin(rng.uniform(0,2*math.pi))
    perp=perp-(perp@ax_)*ax_; perp/=np.linalg.norm(perp)
    return ax_,perp

res={}
def study(name, builder, n=150, gaps=(-0.03,0.0,0.25,0.5,1.0,2.0)):
    out={}
    for g in gaps:
        cnt={20:[],24:[],26:[]}
        for _ in range(n):
            m=builder(g)
            for t in cnt: cnt[t].append(sealed(m,t))
        out[g]={t:(np.mean(v),np.max(v),np.mean(np.array(v)>0)) for t,v in cnt.items()}
    res[name]=out

r=3.0; L=18.0
def two_cyl(g):
    R=rot_random(); ax_,perp=pair_axes('c',g,R)
    c1=C0+rand_offset()-0.5; c2=c1+perp*(2*r+g)
    cm=(c1+c2)/2; c1=c1-(cm-C0); c2=c2-(cm-C0)
    return cyl(c1,R,r,L)|cyl(c2,R,r,L)
def three_cyl(g):
    R=rot_random(); ax_=R[:,2]; e1=R[:,0]; e2=R[:,1]
    d=2*r+g
    off=[np.zeros(3), e1*d, e1*d/2+e2*d*math.sqrt(3)/2]
    cen=sum(off)/3; base=C0+rand_offset()-0.5
    return sum(0 for _ in [])+np.logical_or.reduce([cyl(base+o-cen,R,r,L) for o in off])
def two_sph(g):
    d=rng.normal(size=3); d/=np.linalg.norm(d); c1=C0+rand_offset()-0.5; c2=c1+d*(2*r+g)
    return sphere(c1,r)|sphere(c2,r)
def three_sph(g):
    R=rot_random(); e1=R[:,0]; e2=R[:,1]; d=2*r+g
    off=[np.zeros(3), e1*d, e1*d/2+e2*d*math.sqrt(3)/2]; cen=sum(off)/3; base=C0+rand_offset()-0.5
    return np.logical_or.reduce([sphere(base+o-cen,r) for o in off])
def two_box(g):   # 6x6x18 boxes, side by side, random orientation
    R=rot_random(); e1=R[:,0]; base=C0+rand_offset()-0.5; d=6+g
    return box(base-e1*d/2,R,(3,3,9))|box(base+e1*d/2,R,(3,3,9))
def three_box(g):
    R=rot_random(); e1=R[:,0]; e2=R[:,1]; base=C0+rand_offset()-0.5; d=6+g
    off=[np.zeros(3), e1*d, e1*d/2+e2*d]; cen=sum(off)/3
    return np.logical_or.reduce([box(base+o-cen,R,(3,3,9)) for o in off])

for nm,b in [("2 cylinders",two_cyl),("3 cylinders",three_cyl),("2 spheres",two_sph),("3 spheres",three_sph),("2 boxes 6x6x18",two_box),("3 boxes 6x6x18",three_box)]:
    study(nm,b,n=120)

print("SEALED FLUID CELLS (fluid cell with >=T of its 26 neighbours obstacle), mean / max / fraction of configs with >=1, over 120 random orientation+sub-cell registrations")
for nm,o in res.items():
    print("\n"+nm)
    print(" gap    T>=20 (mean,max,frac)      T>=24 (mean,max,frac)      T>=26")
    for g,d in o.items():
        print(f" {g:5.2f}  "+"   ".join(f"{d[t][0]:5.1f},{d[t][1]:3.0f},{d[t][2]:.2f}" for t in (20,24,26)))

# Study 2: volume fidelity of the cell-centre mapping
print("\nMAPPED VOLUME / ANALYTIC VOLUME (cell-centre test), random orientation + registration, n=400")
def vol_stat(name,fn,Vex,n=400):
    v=np.array([fn().sum()/Vex for _ in range(n)])
    print(f" {name:26s} mean {v.mean():.4f}  std {100*v.std():.2f}%  min {v.min():.3f} max {v.max():.3f}")
vol_stat("sphere r=3",lambda: sphere(C0+rand_offset()-0.5,3.0),4/3*math.pi*27)
vol_stat("sphere r=8",lambda: sphere(C0+rand_offset()-0.5,8.0),4/3*math.pi*512)
vol_stat("cylinder r=3 L=18",lambda: cyl(C0+rand_offset()-0.5,rot_random(),3.0,18.0),math.pi*9*18)
vol_stat("box 6x6x18",lambda: box(C0+rand_offset()-0.5,rot_random(),(3,3,9)),6*6*18)
vol_stat("box 30x60x30 (sweep box, 24^3 grid too small -> skipped)",lambda: np.zeros((1,)),1.0,n=1) if False else None
# aligned (axis-parallel) flat faces: registration sweep of a box, the classic step function
print("\nAXIS-ALIGNED box 6x6x18: mapped volume vs sub-cell shift of one face (x offset 0..1 cell)")
for s in np.linspace(0,1,11):
    m=box(C0+np.array([s,0.0,0.0]),np.eye(3),(3,3,9)); print(f"  shift {s:.1f}: {m.sum()/(6*6*18):.3f}", end=';')
print()
