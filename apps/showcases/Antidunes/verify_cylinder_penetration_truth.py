import json, numpy as np
R, H = 3.0, 9.0
d = json.load(open("pair176.json"))
UA, UB = 4323455642275676176, 4323455642275676180
def axis(r):
    w,x,y,z=[float(r[k]) for k in("quatR","quatI","quatJ","quatK")]; n=np.sqrt(w*w+x*x+y*y+z*z); w,x,y,z=w/n,x/n,y/n,z/n
    return np.array([2*(x*z+w*y),2*(y*z-w*x),1-2*(x*x+y*y)])
def h(n,a): c=abs(n@a); return H*c+R*np.sqrt(max(0,1-c*c))
def s_along(n,a1,a2,dv): return h(n,a1)+h(n,a2)-dv@n          # distance B must move along +n (from A to B) to separate
# dense direction sampling + local refinement for the true minimum translation distance
rng=np.random.default_rng(0)
def true_pd(a1,a2,dv):
    N=200000; v=rng.normal(size=(N,3)); v/=np.linalg.norm(v,axis=1)[:,None]
    c1=np.abs(v@a1); c2=np.abs(v@a2)
    s=H*c1+R*np.sqrt(np.clip(1-c1**2,0,None))+H*c2+R*np.sqrt(np.clip(1-c2**2,0,None))-v@dv
    i=np.argmin(s); n=v[i]; best=s[i]
    for sc in (0.05,0.01,0.002,0.0004):
        for _ in range(300):
            m=n+sc*rng.normal(size=3); m/=np.linalg.norm(m); val=s_along(m,a1,a2,dv)
            if val<best: best,n=val,m
    return best,n
st={}
for r in d["state"]: st[(int(r["mainStep"]),int(r["subCycle"]),int(r["uid"]))]=r
print("step.sc |  reported by contact detection            | ground truth")
for c in sorted(d["contacts"],key=lambda r:(int(r["mainStep"]),int(r["subCycle"]))):
    s,sc=int(c["mainStep"]),int(c["subCycle"])
    if (s,sc,UA) not in st or (s,sc,UB) not in st: continue
    A=st[(s,sc,UA)]; B=st[(s,sc,UB)]
    pa=np.array([float(A[k]) for k in("posX","posY","posZ")]); pb=np.array([float(B[k]) for k in("posX","posY","posZ")]); a1=axis(A); a2=axis(B); dv=pb-pa
    pd,nt=true_pd(a1,a2,dv)
    n=np.array([float(c["normalX"]),float(c["normalY"]),float(c["normalZ"])]); rep=-float(c["penetrationDepth"])
    # required separation along the REPORTED normal (either orientation; keep the smaller)
    need=min(s_along(n,a1,a2,dv),s_along(-n,a1,a2,dv))
    print(f"{s}.{sc} | {c['detectionPath'][22:60]:38s} overlap={rep:8.4f}  needed along its normal={need:7.3f} | true min-translation depth={pd:6.3f}  (centre dist {np.linalg.norm(dv):.2f}, axis angle {np.degrees(np.arccos(min(1,abs(a1@a2)))):.1f} deg)")
