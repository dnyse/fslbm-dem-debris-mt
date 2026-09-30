import glob,sys,collections,statistics
D=sys.argv[1]
for rho in ("075","050","025"):
    h=collections.defaultdict(float); o=collections.defaultdict(int)
    for f in glob.glob(f"{D}/{rho}/rank_*.csv"):
        for l in open(f).read().splitlines()[1:]:
            p=l.split(",")
            if len(p)<6: continue
            k=(int(p[0]),int(p[2]),int(p[3])); h[k]+=float(p[4]); o[k]+=int(p[5])
    steps=sorted({k[0] for k in h}); print("rho",rho,"steps logged:",len(steps),"last",steps[-1])
    for s in steps[::max(1,len(steps)//10)]+[steps[-1]]:
        v=[h[k] for k in h if k[0]==s and o[k]==0]; allv=[h[k] for k in h if k[0]==s]
        print(f"  step {s:6d} median(no obstacle) {statistics.median(v):7.3f}  min {min(v):7.3f} max {max(v):7.3f}  ncols {len(v)}/{len(allv)}")
