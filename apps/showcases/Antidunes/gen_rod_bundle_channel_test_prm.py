#!/usr/bin/env python3
"""LOCAL test prm: ONE bundle of three mutually touching PARALLEL rods dragged through still water (THESIS_PLAN.md 2026-09-21, 'why rods and not spheres').

In the PID-off GJK crash 4252482 (t=16377) the failing interface cell sat in the thin triangular channel between three parallel rods (..8404, ..8407, ..8408) that stayed locked together (surface gaps -0.03 LU) for >= 7000 steps, moving
at ~0.012 and tilted (axis (-0.74,0.65,-0.15)). Three touching parallel cylinders of radius r enclose a channel of inscribed radius 0.1547 r = 0.46 LU (< half a cell). This prm reproduces ONLY that geometry, at the production unit system
(omega 1.99753, g 5.8258e-5, sigma 1.63945e-3) and production rod size (r=3, L=18), in a small box (48x36x48) periodic in X and Y with a free surface, on 9 MPI ranks (3x3x1 blocks of 16x12x48 cells):

    python3 gen_rod_bundle_channel_test_prm.py [--vx -0.0096 --vy 0.0043 --vz 0 --dz 0 --dx 0 --dy 0 --steps 20000 --axis -0.74,0.65,-0.15 --name rod_bundle_channel_test]
    mpirun -np 9 <build-release>/apps/showcases/Antidunes/Antidunes prms/<name>_gas.prm     (from an empty working directory)

The rods move RIGIDLY with the constant velocity ModelParameters.rodPrescribedVelocityX/Y/Z (zero angular velocity, forces ignored) so that cells are covered / uncovered by the bundle without any DEM dynamics; the water starts at rest.
"""
import argparse, math, pathlib, re
ap = argparse.ArgumentParser()
ap.add_argument("--vx", type=float, default=-0.0096); ap.add_argument("--vy", type=float, default=0.0043); ap.add_argument("--vz", type=float, default=0.0)
ap.add_argument("--dx", type=float, default=0.0, help="sub-cell offset of the bundle centroid (registration)"); ap.add_argument("--dy", type=float, default=0.0); ap.add_argument("--dz", type=float, default=0.0)
ap.add_argument("--steps", type=int, default=20000)
ap.add_argument("--axis", default="-0.74,0.65,-0.15", help="rod axis (default: the axis of the bundle in the crash run at t=17000)")
ap.add_argument("--gap", type=float, default=-0.03, help="surface gap between the rods (crash bundle: -0.01..-0.03 LU)")
ap.add_argument("--name", default="rod_bundle_channel_test")
ap.add_argument("--cpb", default="", help="override cellsPerBlock, e.g. 12,9,48 (16 ranks) -- decomposition-independence test")
ap.add_argument("--guard", type=int, default=0, help="ModelParameters.closeSealedCellsMinNeighbours (0 = off)")
ap.add_argument("--accel", type=float, default=1e-15, help="accelerationXOverride (driving body acceleration); default = still water")
ap.add_argument("--flowtest", action="store_true", help="flow-rate evaluator test: ClogMeasurement.flowRateEnabled, planes at x=15 (through the bundle) and x=40 (rod-free), sampled every 100 steps")
ap.add_argument("--nobundle", action="store_true", help="control: three SPHERES-like isolated rods far apart instead of a bundle (spacing 14 LU)")
ap.add_argument("--pair", action="store_true", help="TWO parallel rods only (matches the real 2-rod parallel-slit geometry of job 4255679's own crash, e.g. --gap -0.069 --axis 0,1,0): lower obstacle-neighbour count than the 3-rod bundle (~9-15/26 vs ~20/26), the sub-mechanism SealedCellGuard.h's default threshold 24 was NOT validated against")
ap.add_argument("--raft", default="", help="R,C: hex-packed raft of R rows x C cols touching rods instead of the 3-rod triangle (log-jam stress test for SealedCellGuard: many long-lived closures + many simultaneous real formerObstacle conversions, unlike the 3-rod case). Bumps the domain/blocks unless --cpb/--dom given.")
ap.add_argument("--dom", default="", help="override domain size 'x,y,z' (default: 48,36,48, or auto-sized for --raft)")
A = ap.parse_args()
base = pathlib.Path(__file__).parent / "prms" / "hopper_clog_cyl_frontblock_pidoff_2node_gjk_sanity_18k_gas.prm"
txt = base.read_text()
if A.raft:
    RAFT_R, RAFT_C = (int(v) for v in A.raft.split(","))
    DOM = tuple(int(v) for v in A.dom.split(",")) if A.dom else (64, 64, 48)
    CPB = tuple(int(v) for v in A.cpb.split(',')) if A.cpb else (16, 16, 48)   # 4x4x1 = 16 ranks
else:
    RAFT_R = RAFT_C = 0
    DOM = tuple(int(v) for v in A.dom.split(",")) if A.dom else (48, 36, 48)
    CPB = tuple(int(v) for v in A.cpb.split(',')) if A.cpb else (16, 12, 48)   # 3x3x1 = 9 ranks
# SetupStages.h:71 needs >= 3 blocks per PERIODIC direction when particles are used; blocks >= the rod interaction radius (9.5)
R_ROD, L_ROD = 3.0, 18.0
ZW = 34.0
# --- geometry: three parallel rods, centres on an equilateral triangle of side 2r+gap in the plane perpendicular to the axis
ax = [float(v) for v in A.axis.split(",")]; n = math.sqrt(sum(a*a for a in ax)); ax = [a/n for a in ax]
up = [0.0, 0.0, 1.0]; d = sum(u*a for u, a in zip(up, ax)); e2 = [u - d*a for u, a in zip(up, ax)]; m = math.sqrt(sum(v*v for v in e2)); e2 = [v/m for v in e2]
e1 = [ax[1]*e2[2]-ax[2]*e2[1], ax[2]*e2[0]-ax[0]*e2[2], ax[0]*e2[1]-ax[1]*e2[0]]
side = 2*R_ROD + A.gap; Rc = side/math.sqrt(3)
C = (DOM[0] / 2.0 + A.dx, DOM[1] / 2.0 + A.dy, 19.6 + A.dz)     # centroid: domain centre (3-rod case: near the original block faces x=16, y=12)
if A.raft:
    # hex-packed raft: RAFT_R rows x RAFT_C cols, row spacing d*sqrt(3)/2 in e2, alternate rows offset by d/2 in e1 (touching neighbours in every direction)
    d = side
    pts = []
    for row in range(RAFT_R):
        for col in range(RAFT_C):
            u = col * d + (d / 2.0 if row % 2 else 0.0)
            v = row * d * math.sqrt(3) / 2.0
            pts.append((u, v))
    cu = sum(p[0] for p in pts) / len(pts); cv = sum(p[1] for p in pts) / len(pts)
    centres = [tuple(C[i] + (p[0] - cu) * e1[i] + (p[1] - cv) * e2[i] for i in range(3)) for p in pts]
elif A.nobundle:
    offs = [(-14.0, 0.0), (0.0, 0.0), (14.0, 0.0)]
    centres = [tuple(C[i] + o[0]*e1[i] + o[1]*e2[i] for i in range(3)) for o in offs]
elif A.pair:
    centres = [tuple(C[i] + s*(side/2.0)*e1[i] for i in range(3)) for s in (-1.0, 1.0)]
else:
    centres = [tuple(C[i] + Rc*(math.cos(math.radians(p))*e1[i] + math.sin(math.radians(p))*e2[i]) for i in range(3)) for p in (90.0, 210.0, 330.0)]
rods = "\n".join(f"   FloatingCylinder {{ position < {c[0]:.4f}, {c[1]:.4f}, {c[2]:.4f} >; radius {R_ROD:g}; length {L_ROD:g}; axis < {ax[0]:.6f}, {ax[1]:.6f}, {ax[2]:.6f} >; densityRatio 1; settlingSteps 0; }}" for c in centres)

def replace_block(t, name, new):
    m = re.search(r"^%s\s*\n\{" % name, t, re.M); assert m, name
    i = m.start(); j = t.index("\n}\n", m.end()) + 3
    return t[:i] + new + t[j:]
def get_block(t, name):
    m = re.search(r"^%s\s*\n\{" % name, t, re.M); i = m.start(); j = t.index("\n}\n", m.end()) + 3; return t[i:j]

t = txt
t = replace_block(t, "BlockForestParameters", f"""BlockForestParameters
{{
   cellsPerBlock                 < {CPB[0]}, {CPB[1]}, {CPB[2]} >;
   periodicity                   < 1, 1, 0 >;
   loadSnapshot                  false;
   storeSnapshot                 false;
   snapshotFrequency             100000;
   snapshotBaseFolder            snapshots;
}}
""")
t = replace_block(t, "DomainParameters", f"""DomainParameters
{{
   domainSize         <{DOM[0]}, {DOM[1]}, {DOM[2]}>;
   wavePeriods        1;
   liquidHeightFactor 2.0;
   floorHeightFactor  {(ZW-16.0031)/8.001943372119138:.6f};   // absoluteLiquidHeight = 16.003 + floor = {ZW} (water surface at z=34, bundle ~14 LU below it)
   initialAmplitude   0;
}}
""")
t = replace_block(t, "PIDParameters", f"""PIDParameters
{{
   targetMeanVelocityMagnitude 0;      // still water
   proportionalGain            2e-4;
   derivativeGain              1e-6;
   integralGain                2e-4;
   maxRamp                     1e-4;
   minActuatingVariable        0;
   maxActuatingVariable        1e-3;
   enabled                     false;
   accelerationXOverride       {A.accel};  // default 1e-15 = effectively no driving force
}}
""")
ph = get_block(t, "PhysicsParameters")
ph_new = f"""PhysicsParameters
{{
   enableWetting     false;
   timesteps         {A.steps};
   Re                3100;
   Fr                1.31;
   We                15.6188;
   useHardcodedUnits true;              // U = 0 (still water): fixed units, then overridden to the PRODUCTION values below
   relaxationRateOverride  1.99753;     // production omega (nu = 2.06502e-4)
   gravityOverride         5.82575e-05; // production gravity
   surfaceTensionOverride  0.00163945;  // production surface tension
   denseNanScanEnabled            true;
   denseNanScanFromStep           0;
   denseNanScanToStep             {A.steps};
   highVelocityPrecursorThreshold 0.5;
   fluidSanityEnabled             true;
   fluidSanityFromStep            0;
   fluidSanityToStep              {A.steps};
   fluidSanityEveryNSteps         1;
   fluidSanityMaxVelocity         0.17;
   fluidSanityDensityMin          0.1;
   fluidSanityDensityMax          1.9;
   fluidSanityMaxPopulation       2.0;
   fluidSanityFillTolerance       1.0;
   fluidSanityMaxExcessMass       1.0;
   fluidSanityMaxNeighbourVelocityJump 0.2;
   fluidSanityMaxNeighbourDensityJump  0.3;
   fluidSanityCheckTopology       true;
   fluidSanityStatsEveryNSteps    100;
   fluidSanityAbortAbove          0;
}}
"""
t = t.replace(ph, ph_new)
# no pillars, three rods
t = replace_block(t, "CylinderPillars", "// (no pillars in this test)\n") if re.search(r"^CylinderPillars", t, re.M) else t
t = replace_block(t, "FloatingCylinders", "FloatingCylinders\n{\n" + rods + "\n}\n")
mp = get_block(t, "ModelParameters")
add = f"""   rodPrescribedVelocityX        {A.vx};    // TEST OPTION: rigid prescribed rod velocity (crash bundle: (-0.0096, 0.0043, 0.0062) at t=16300; z component dropped, the box is small)
   rodPrescribedVelocityY        {A.vy};
   rodPrescribedVelocityZ        {A.vz};
""" + (f"   closeSealedCellsMinNeighbours {A.guard};   // EXPERIMENT: close sub-cell channels (SealedCellGuard.h)\n" if A.guard else "")
mp_new = mp.rstrip("\n")
assert mp_new.endswith("}"); mp_new = mp_new[:-1].rstrip("\n") + "\n" + add + "}\n"
t = t.replace(mp, mp_new)
t = replace_block(t, "BoundaryParameters", "BoundaryParameters\n{\n   Border { direction T;  walldistance -1; NoSlip{} }\n   Border { direction B;  walldistance -1; NoSlip{} }\n}\n")
ev = get_block(t, "EvaluationParameters")
ev_new = re.sub(r"upstreamPos\s*<[^>]*>", "upstreamPos       < 6, 18, 20 >", ev)
ev_new = re.sub(r"downstreamPos\s*<[^>]*>", "downstreamPos     < 42, 18, 20 >", ev_new)
ev_new = re.sub(r"pillarUpstreamPos\s*<[^>]*>", "pillarUpstreamPos    < 6, 6, 20 >", ev_new)
ev_new = re.sub(r"pillarDownstreamPos\s*<[^>]*>", "pillarDownstreamPos  < 42, 6, 20 >", ev_new)
ev_new = re.sub(r"performanceLogFrequency\s+\d+", "performanceLogFrequency 5000", ev_new)
if A.flowtest:
    ev_new = re.sub(r"upstreamPos\s*<[^>]*>", "upstreamPos       < 15, 18, 20 >", ev_new)
    ev_new = re.sub(r"downstreamPos\s*<[^>]*>", "downstreamPos     < 40, 18, 20 >", ev_new)
    ev_new = re.sub(r"(probeHalfWidthXY\s+\d+;)", r"\1\n      flowRateEnabled   true;", ev_new, count=1)
t = t.replace(ev, ev_new)
vt = get_block(t, "VTK")
vt_new = re.sub(r"CellBB_filter_0 \{[^}]*\}", f"CellBB_filter_0 {{ min < 0, 0, 0 >;  max < {DOM[0]}, {DOM[1]}, {DOM[2]} >; }}", vt)
vt_new = re.sub(r"writeFrequency\s+\d+;", "writeFrequency       500;", vt_new, count=1)
t = t.replace(vt, vt_new)
hdr = f"""// {A.name}_gas.prm   (GENERATED by gen_rod_bundle_channel_test_prm.py -- edit the generator, not this file)
// [Claude, 2026-09-21] LOCAL bundle-channel test: {"three ISOLATED rods (control, no channel)" if A.nobundle else ("TWO parallel rods (gap %.3f LU, r=3, L=18) = the real 2-rod slit geometry of job 4255679's own crash" % A.gap if A.pair else "three mutually touching PARALLEL rods (gap %.3f LU, r=3, L=18) = the bundle at the crash site of job 4252482" % A.gap)}, dragged through still water with the constant velocity
// ({A.vx}, {A.vy}, {A.vz}); axis ({ax[0]:.3f}, {ax[1]:.3f}, {ax[2]:.3f}); domain {DOM[0]}x{DOM[1]}x{DOM[2]} periodic in X and Y, water surface at z={ZW:g}, 3x3x1 blocks of {CPB[0]}x{CPB[1]}x{CPB[2]} (9 MPI ranks); centroid start ({C[0]:.2f},{C[1]:.2f},{C[2]:.2f}).
// Production unit system (omega 1.99753, g 5.8258e-5, sigma 1.63945e-3); fluid sweep at |u|>0.17 + dense NaN scan armed from step 0. Everything else as in prms/hopper_clog_cyl_frontblock_pidoff_2node_gjk_sanity_18k_gas.prm.
"""
out = pathlib.Path(__file__).parent / "prms" / f"{A.name}_gas.prm"
out.write_text(hdr + t)
print("wrote", out); print("rod centres:", [tuple(round(x, 3) for x in c) for c in centres])
pd = lambda a, b: math.dist(a, b)
print("pairwise centre distances:", [round(pd(centres[i], centres[j]), 3) for i in range(len(centres)) for j in range(i+1, len(centres))], "(=2r+gap =", 2*R_ROD+A.gap, ")")
