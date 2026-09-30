#!/usr/bin/env python3
"""Census of 'sealed' LIQUID cells in fluid_field VTK frames (needs numpy).
Usage: python3 census_sealed_cells.py <path to vtk-out/fluid_field>   (default: animations/4252146_gjk_sanity_18k/vtk-out/fluid_field)
For every frame: number of liquid cells (mapped_flag 1) with >= 12/16/20/24/25/26 obstacle cells among their 26 neighbours (obstacle = mapped_flag 0 = rods, or >= 4 = walls/pillars;
cells outside the VTK crop count as non-obstacle, so the numbers are lower bounds at the crop edge) and the max |u| among the >=24 ones. [Claude, 2026-09-21] see THESIS_PLAN.md.
"""
import os,re,sys,base64,struct,glob
import numpy as np
ROOT=sys.argv[1] if len(sys.argv)>1 else "animations/4252146_gjk_sanity_18k/vtk-out/fluid_field"
def arr(txt,name=None,dtype=None):
    pat=r'<DataArray type="(\w+)"[^>]*?%sformat="binary">\s*([A-Za-z0-9+/=]+)\s*</DataArray>'%(('Name="%s" [^>]*?'%name) if name else '')
    m=re.search(pat,txt,re.S); typ,b=m.group(1),m.group(2); raw=base64.b64decode(b); n=struct.unpack("<I",raw[:4])[0]
    dt={"Float32":"<f4","Float64":"<f8","UInt32":"<u4"}[typ]; return np.frombuffer(raw[4:4+n],dtype=dt)
def load_frame(step):
    D=f"{ROOT}/simulation_step_{step}"
    blocks=[]
    for fn in os.listdir(D):
        txt=open(D+"/"+fn).read(); i=txt.index("<Points>"); j=txt.index("</Points>")
        pts=arr(txt[i:j]).reshape(-1,3); x0,y0,z0=pts.min(0); x1,y1,z1=pts.max(0); nx,ny,nz=int(round(x1-x0)),int(round(y1-y0)),int(round(z1-z0))
        ci=txt.index("<CellData>"); fl=arr(txt[ci:],"mapped_flag"); ve=arr(txt[ci:],"velocity").reshape(-1,3)
        blocks.append((int(x0),int(y0),int(z0),nx,ny,nz,fl.reshape(nz,ny,nx),np.linalg.norm(ve,axis=1).reshape(nz,ny,nx)))
    X0=min(b[0] for b in blocks); Y0=min(b[1] for b in blocks); Z0=min(b[2] for b in blocks)
    X1=max(b[0]+b[3] for b in blocks); Y1=max(b[1]+b[4] for b in blocks); Z1=max(b[2]+b[5] for b in blocks)
    F=np.full((Z1-Z0,Y1-Y0,X1-X0),255,np.uint8); U=np.zeros(F.shape,np.float32)
    for x0,y0,z0,nx,ny,nz,fl,u in blocks:
        F[z0-Z0:z0-Z0+nz,y0-Y0:y0-Y0+ny,x0-X0:x0-X0+nx]=fl; U[z0-Z0:z0-Z0+nz,y0-Y0:y0-Y0+ny,x0-X0:x0-X0+nx]=u
    return F,U,(X0,Y0,Z0)
def obst_count(F):
    isobs=((F==0)|((F>=4)&(F!=255))).astype(np.uint8)   # rods (0) and walls/pillars (>=4); 255 = outside crop -> not counted
    P=np.pad(isobs,1); cnt=np.zeros(F.shape,np.uint8)
    for dz in (-1,0,1):
        for dy in (-1,0,1):
            for dx in (-1,0,1):
                if dx==dy==dz==0: continue
                cnt+=P[1+dz:1+dz+F.shape[0],1+dy:1+dy+F.shape[1],1+dx:1+dx+F.shape[2]]
    return cnt
steps=sorted(int(d.split("_")[-1]) for d in os.listdir(ROOT) if d.startswith("simulation_step_") and not d.endswith("pvtu"))
if len(sys.argv)>2: steps=[t for t in steps if t in {int(x) for x in sys.argv[2].split(",")}]
print("step | liquid cells | body cells (flag 0) | obstacle-cell count of LIQUID cells with obst nbrs >= 12 / 16 / 20 / 24 / 25 / 26 | max |u| among those with >=24")
out=[]
for s in steps:
    F,U,o=load_frame(s); C=obst_count(F); liq=(F==1)
    ks=[int(((C>=k)&liq).sum()) for k in (12,16,20,24,25,26)]
    sealed=(C>=24)&liq
    mu=float(U[sealed].max()) if sealed.any() else 0.0
    nb=int((F==0).sum())
    print("%5d | %8d | %8d | %s | %.3f | per 1000 body cells >=24: %.2f  >=20: %.1f"%(s,int(liq.sum()),nb," ".join("%5d"%k for k in ks),mu,1000.0*ks[3]/max(nb,1),1000.0*ks[2]/max(nb,1)),flush=True)
    idx=np.argwhere((C>=25)&liq)
    out.append((s,[(int(x+o[0]),int(y+o[1]),int(z+o[2]),int(C[z,y,x]),round(float(U[z,y,x]),3)) for z,y,x in idx[:12]]))
print("\nexamples of cells with >=25 obstacle neighbours (x,y,z cell index, nbrs, |u|):")
for s,ex in out:
    if ex: print(s,ex[:8])
