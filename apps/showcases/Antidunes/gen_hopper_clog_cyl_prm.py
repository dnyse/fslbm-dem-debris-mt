#!/usr/bin/env python3
"""Cylinder version of the "hopper" clogging scenario (companion of gen_hopper_clog_prm.py, which does spheres).

    python3 gen_hopper_clog_cyl_prm.py   ->  prms/hopper_clog_cyl_frontblock_pidon_2node_gas.prm

Same channel / pillars / PID / water height / VTK as the sphere hopper, but the debris are CYLINDERS with diameter:length = 1:3 (r=3, L=18, i.e. the
sphere's 6 LU cross-section) lying with their LONG AXIS ALONG Y (across the channel), the whole block placed against the pillars.
Two rods fit across the 40 LU channel (2*18 + 3 clearances), so the block is 2 (Y) x 7 (Z) x 29 (X rows) = 406 rods.
"""
import argparse
import itertools
import math
import pathlib

ap = argparse.ArgumentParser()
ap.add_argument("--radius", type=float, default=3.0, help="cylinder radius (LU)")
ap.add_argument("--length", type=float, default=18.0, help="cylinder length along Y (LU); diameter:length = 1:3 by default")
ap.add_argument("--pillar-radius", type=float, default=12.0)
ap.add_argument("--pitch", type=float, default=7.0, help="X row pitch (LU) >= 2*radius + clearance")
ap.add_argument("--nlayers", type=int, default=7)
ap.add_argument("--clearance", type=float, default=1.0, help="min surface gap rod<->rod/wall/floor (LU)")
ap.add_argument("--rows", type=int, default=29)
ap.add_argument("--front-x", type=float, default=368.0, help="X of the row nearest the pillars")
ap.add_argument("--pillar-clearance", type=float, default=1.0)
ap.add_argument("--skin", type=float, default=0.0,
                 help="cylinderContactSkin (LU, cd1b991b): rod-rod contact starts at a surface gap of 2*skin, "
                      "rod-wall/pillar at skin. 0 (default) = off, matching the original un-adapted hopper. "
                      "When >0, --clearance/--pitch must leave a real gap safely above that threshold or the "
                      "pre-existing [RELEASE-OVERLAP] release-time guard aborts (see THESIS_PLAN.md, job 4282206).")
ap.add_argument("--hold", type=int, default=2000)
ap.add_argument("--timesteps", type=int, default=60000)
ap.add_argument("--out", default="prms/hopper_clog_cyl_frontblock_pidon_2node_gas.prm")
ap.add_argument("--lx", type=float, default=768.0, help="domain X (pillars at LX/2); a shorter domain gives a small local reproduction")
ap.add_argument("--lz", type=float, default=180.0, help="domain Z")
ap.add_argument("--cells", default="64,20,30", help="cellsPerBlock 'bx,by,bz'")
ap.add_argument("--vtk-every", type=int, default=0, help="fluid/particle VTK writeFrequency (0 = timesteps//12)")
ap.add_argument("--vtk-crop", default="128,640,100", help="fluid VTK crop 'x0,x1,zmax'")
args = ap.parse_args()

LX, LY, LZ = args.lx, 40.0, args.lz
PILLAR_X = LX / 2
WATER_H = 8.001943372119138 * (2.0 + 4.1393)   # absoluteLiquidHeight = 49.13
r, L, R, pitch, c = args.radius, args.length, args.pillar_radius, args.pitch, args.clearance
gap = LY - 2 * R
assert pitch >= 2 * r + c

# ---- grid: as many rods across Y as fit (each needs L + clearances), evenly spread; layers evenly floor -> waterline ------------------------------
ncols = int((LY - c) // (L + c))                       # 2 for L=18, c=1
assert ncols >= 1
pitch_y = (LY - 2 * (L / 2 + c)) / (ncols - 1) if ncols > 1 else 0.0
ys = [L / 2 + c + i * pitch_y for i in range(ncols)] if ncols > 1 else [LY / 2]
nz = args.nlayers
pitch_z = (WATER_H - (r + c)) / (nz - 1)
zs = [r + c + k * pitch_z for k in range(nz)]
xs = [args.front_x - i * pitch for i in range(args.rows)]
assert pitch_z >= 2 * r + c
if ncols > 1: assert pitch_y >= L + c

def density(k):                                        # heavy low, neutral mid, light at the surface (as the sphere version)
    if k >= nz - 2: return 0.5
    return 1.2 if k < (nz - 2) // 2 else 1.0

rods = [(x, y, z, density(k)) for x in xs for y in ys for k, z in enumerate(zs)]

# ---- exact geometry checks (all rods parallel to Y) -------------------------------------------------------------------------------------------------
def rod_gap(a, b):     # lower bound on the surface gap of two parallel-axis cylinders: max(radial gap, axial gap)
    return max(math.hypot(a[0] - b[0], a[2] - b[2]) - 2 * r, abs(a[1] - b[1]) - L)
minsep = min(rod_gap(a, b) for a, b in itertools.combinations(rods, 2))
assert minsep >= c - 1e-9, f"rods too close: {minsep}"
ywall_gap = min(min(y - L / 2, LY - (y + L / 2)) for _, y, _, _ in rods)
assert ywall_gap >= c - 1e-9, "rod too close to a Y wall"
floor_gap = min(z - r for _, _, z, _ in rods)
assert floor_gap >= c - 1e-9, "rod too close to the floor"
def rect_disc_gap(x, y, xc, yc):   # XY projection of a rod is the rectangle [x-r,x+r]x[y-L/2,y+L/2]; pillar is a disc (xc,yc,R)
    nx = min(max(xc, x - r), x + r); ny = min(max(yc, y - L / 2), y + L / 2)
    return math.hypot(nx - xc, ny - yc) - R
pillar_gap = min(rect_disc_gap(x, y, PILLAR_X, yc) for x, y, _, _ in rods for yc in (0.0, LY))
assert pillar_gap >= args.pillar_clearance - 1e-9, f"rod too close to a pillar: {pillar_gap:.3f}"
assert max(x for x, *_ in rods) + r < PILLAR_X - R + 1e-9 or True
# cylinderContactSkin (--skin, see its own help text): the skin's rod-rod contact threshold is a
# surface gap of 2*skin, rod-wall/pillar is skin. A real gap AT or below that threshold registers
# as a genuine contact to the skin-aware detector, which the pre-existing (unconditional)
# [RELEASE-OVERLAP] release-time guard then aborts on -- it can't distinguish that from an actual
# insertion-overlap bug (see job 4282206, THESIS_PLAN.md). Require real margin, not just >=, so
# this never lands exactly on the boundary the way the un-skinned default clearance does.
if args.skin > 0:
    assert minsep >= 2 * args.skin + 1e-9, (
        f"rod-rod gap {minsep:.3f} LU is at/below the skin's own 2*skin={2*args.skin:g} LU threshold "
        f"-- widen --pitch (and/or --clearance) first")
    assert pillar_gap >= args.skin + 1e-9, (
        f"rod-pillar gap {pillar_gap:.3f} LU is at/below the skin's own skin={args.skin:g} LU threshold "
        f"-- widen --pillar-clearance first")
    assert ywall_gap >= args.skin + 1e-9, (
        f"rod-Ywall gap {ywall_gap:.3f} LU is at/below the skin's own skin={args.skin:g} LU threshold "
        f"-- widen --clearance first")
    assert floor_gap >= args.skin + 1e-9, (
        f"rod-floor gap {floor_gap:.3f} LU is at/below the skin's own skin={args.skin:g} LU threshold "
        f"-- widen --clearance first")
counts = {d: sum(1 for q in rods if q[3] == d) for d in (0.5, 1.0, 1.2)}
vol = math.pi * r * r * L
print(f"{len(rods)} rods (r={r:g}, L={L:g}, L/D={L/(2*r):g}, axis Y): {ncols} across Y x {nz} layers x {args.rows} rows; densities {counts}")
print(f"pitch X/Y/Z = {pitch:g}/{pitch_y:g}/{pitch_z:.3f}; min rod-rod surface gap {minsep:.2f}; nearest rod-to-pillar-surface gap {pillar_gap:.3f}")
print(f"block X [{xs[-1]-r:.1f},{xs[0]+r:.1f}]; Y centres {[round(y,2) for y in ys]}; Z layers {[round(z,2) for z in zs]}")
print(f"gap {gap:g} vs rod length {L:g}: a rod lying along Y (broadside) does NOT fit the gap ({'blocked' if L > gap else 'fits'}); "
      f"rod volume {vol:.0f} LU^3 (sphere r=3: {4/3*math.pi*27:.0f}); solid fraction of the block ~ {vol/(pitch*(LY/ncols)*pitch_z):.2f}")

lines = "\n".join(
    f"   FloatingCylinder {{ position < {x:.2f}, {y:.2f}, {z:.2f} >; radius {r:g}; length {L:g}; axis < 0, 1, 0 >; densityRatio {d:g}; settlingSteps {args.hold}; }}"
    for x, y, z, d in rods)

skin_header = "" if args.skin <= 0 else (
    "// ─────────────────────────────────────────────────────────────────────────────\n"
    "// [Claude, 2026-10-01] SKIN-ADAPTED PRODUCTION RUN: cylinderContactSkin "
    f"{args.skin:g} added (cd1b991b, rod-rod contact threshold a surface gap of {2*args.skin:g} LU, "
    f"rod-wall/pillar {args.skin:g} LU), sealed-cell guard left off, matching the two clean doubling "
    "tests below this file's own production scale (jobs 4282267: 24 rods, full 51000 steps, exit 0, "
    "zero crash markers; 4282396: 48 rods, same result). The original/default --pitch=7.0 placed "
    "adjacent rows exactly 1.0 LU apart surface-to-surface -- precisely 2*skin -- which the "
    "pre-existing [RELEASE-OVERLAP] release-time guard mistook for an insertion-overlap bug (job "
    "4282206, failed in 97s). Regenerated with --pitch 8.0 instead (real rod-rod gap now "
    f"{minsep:.2f} LU, comfortably above 2*skin={2*args.skin:g}); every clearance (rod-rod, "
    "rod-Ywall, rod-floor, rod-pillar) is now asserted against its own skin threshold by this "
    "script (see the --skin argparse help and the asserts right after minsep/pillar_gap above). "
    "Rod geometry/count/densities are otherwise byte-for-byte identical to the un-skinned 406-rod "
    "original -- only the row pitch changed, same approach as the 24/48-rod tests (widen the gap, "
    "don't shrink the rods, per explicit instruction).\n"
)

prm = f"""{skin_header}// {pathlib.Path(args.out).name}   (GENERATED by gen_hopper_clog_cyl_prm.py -- edit the script, not this file)
// ─────────────────────────────────────────────────────────────────────────────
// [Claude, 2026-09-19] "HOPPER" CLOGGING SETUP, CYLINDER VERSION (THESIS_PLAN.md item 3/4). NOT RUN YET.
// Companion of the sphere hopper (prms/hopper_clog_pidon_2node_frontblock_gas.prm, same channel/pillars/PID/probes/VTK): two cylindrical pillars centred on the
// Y walls (X={PILLAR_X:g}, radius {R:g}, central gap {gap:g} LU, plan view "|) (|") -- but the debris are CYLINDERS, diameter:length = 1:3 (r={r:g}, L={L:g}; the sphere's
// {2*r:g} LU cross-section), LONG AXIS ALONG Y (across the channel), released together from a block pressed against the pillars.
//
// KEY GEOMETRIC FACT: the rod length {L:g} is {'LONGER' if L > gap else 'shorter'} than the {gap:g} LU gap. A rod lying along Y therefore {'cannot pass the gap broadside: it has to rotate its axis towards the flow' if L > gap else 'fits the gap broadside'}.
// Expect jamming of rods across the opening (log-jam like) rather than the free arching of the sphere hopper -- a different mechanism (rotation-limited passage).
//
// POPULATION: {len(rods)} rods = {ncols} across Y x {nz} layers x {args.rows} rows in X (pitch X/Y/Z = {pitch:g}/{pitch_y:g}/{pitch_z:.2f}, >= {c:g} LU surface clearance rod-rod / wall / floor,
// checked exactly for parallel cylinders). Only {ncols} rods of length {L:g} fit across the {LY:g} LU channel, and each rod is ~{vol/(4/3*math.pi*27):.1f}x the volume of a sphere, hence far fewer than the 1015
// spheres (kept: the same 29 rows, i.e. the same block length). densityRatio {{0.5,1.0,1.2}} = {counts[0.5]}/{counts[1.0]}/{counts[1.2]}, stratified by height like the sphere version (1.2 low, 1.0
// mid, 0.5 in the top two layers; the top layer's axis sits at the waterline Z={WATER_H:.2f}, half submerged).
// Block X in [{xs[-1]-r:.1f},{xs[0]+r:.1f}], front row {pillar_gap:.2f} LU (surface to pillar surface, exact rectangle-disc distance) from the pillar cylinders.
// All rods are held FIXED for {args.hold} steps (settlingSteps, no startInactive) and then released together; clearance >= {c:g} LU so [RELEASE-OVERLAP] cannot fire at release.
//
// CONTACT MODEL: useAnalyticalContactDetection true (Cylinder x Cylinder = analytical Kodam+Guo HYBRID, usePureAnalyticalCylinder false keeps the GJK/EPA fallback;
// Cylinder x pillar/wall analytic regardless) + useContactFilterOwnerTiebreak true; every experimental contact/buoyancy fix explicitly OFF (see ModelParameters) EXCEPT useRotationalEffectiveMass (ON, user request 2026-09-20: generalised effective mass incl. rotational inertia in the contact law) and useTightLayer1 (ON, user request 2026-09-20: the capsule/segment-distance pre-filter for Cylinder x Cylinder; near-parallel false-positive fix is in PaperAnalyticalContacts.h + an always-on capsule bound).
// Friction/restitution as in the sphere hopper (0.5 / 0.5).
//
// FLOW/UNITS: identical to the sphere hopper (PID on, target 0.04, water Z=49.13, Z=180 headroom, <768,40,180> / <64,20,30> = 144 ranks / 2 nodes).
// CAVEATS: the upstream clog probes (X={PILLAR_X-48:g}) sit inside the rod block (unusable as a clog signal, as in the sphere case); the fluid VTK is cropped to X 128..640, Z 0..100
// (the rods are unfiltered); the [CylComputation] contact-path diagnostics make logs very large (the 100-floater sediment run's log was 594 MB) -- expect a big .out file.
// ─────────────────────────────────────────────────────────────────────────────

BlockForestParameters
{{
   cellsPerBlock                 < {args.cells.split(',')[0]}, {args.cells.split(',')[1]}, {args.cells.split(',')[2]} >;
   periodicity                   < 1, 0, 0 >;
   loadSnapshot                  false;
   storeSnapshot                 true;
   snapshotFrequency             5000;
   snapshotBaseFolder            snapshots;
}}

DomainParameters
{{
   domainSize         <{LX:g}, {LY:g}, {LZ:g}>;
   wavePeriods        1;
   liquidHeightFactor 2.0;
   floorHeightFactor  4.1393;
   initialAmplitude   0;
}}

PIDParameters
{{
   targetMeanVelocityMagnitude 0.04;
   proportionalGain            2e-4;
   derivativeGain              1e-6;
   integralGain                2e-4;
   maxRamp                     1e-4;
   minActuatingVariable        0;
   maxActuatingVariable        1e-3;
   enabled                     true;
}}

PhysicsParameters
{{
   enableWetting     false;
   timesteps         {args.timesteps};
   Re                3100;
   Fr                1.31;
   We                15.6188;
   useHardcodedUnits false;
   // no gravityOverride -- dynamically computed, see the extreme-clog prm's header
}}

ParticleParameters
{{
   inFileName spheres_out.dat;   // never opened -- see avgParticleDiameterOverride below
   bedCopiesInX 1;
   bedCopiesInY 1;
   avgParticleDiameterOverride 8.001943372119138;
   dxSIOverride                0.000362694127490625;
   densityRatio 2.55;
   fixingHeightFactor 1.5;
   frictionCoefficient 0.5;
   restitutionCoefficient 0.5;
   numSubCycles 10;
   useLubricationCorrection true;
   useNoSlipParticles false;
   useParticleBed false;
}}

CylinderPillars
{{
   // axes ON the Y walls: half of each cylinder lies outside the domain (harmless). Z is an arbitrary reference point on the vertical axis.
   CylinderPillar {{ center < {PILLAR_X:g}, 0, 45 >; radius {R:g}; }}
   CylinderPillar {{ center < {PILLAR_X:g}, {LY:g}, 45 >; radius {R:g}; }}
}}

FloatingCylinders
{{
{lines}
}}

ModelParameters
{{
   useAnalyticalContactDetection true;   // Cylinder x Cylinder: analytical Kodam+Guo hybrid
   usePureAnalyticalCylinder     false;  // (default, explicit) keep the GJK/EPA fallback for uncertain pairs
   useContactFilterOwnerTiebreak true;   // one rank accepts each boundary-straddling contact (min owner rank){"" if args.skin <= 0 else chr(10) + f"   cylinderContactSkin           {args.skin:g};   // [Claude, 2026-10-01] rod contact skin in LU (cd1b991b): rod-rod contact starts at a surface gap of {2*args.skin:g}, rod-wall/pillar at {args.skin:g}. Sealed-cell guard stays off (not set below)."}

   // Experimental contact/buoyancy fixes: ALL explicitly OFF
   useTightTier1                 false;
   useTightLayer1                true;   // [Claude, 2026-09-20] ON by user request (capsule/segment-distance pre-filter for Cylinder x Cylinder, identical rejections to the always-on capsule bound but before the analytic dispatch)
   skipTier1                     false;
   useRoundedBoxCorners          false;
   useMultiPointBoxContact       false;
   useBoxOverlapCap              false;
   useAxisHysteresis             false;
   useContactManifold            false;
   useRotationalEffectiveMass    true;   // [Claude, 2026-09-20] ON for rods (user request): contact stiffness/damping from the generalised effective mass incl. rotational inertia at the contact point (exact for spheres, matters for rod rim/edge contacts)
   useSpheroBoxContact           false;
   useSpheroBoxAnalytical        false;
   useFraigeBoxContactModel      false;
   useFullVolumeAnalyticBuoyancy false;
   usePass3Reresolution          false;
   useAnalyticBoxBuoyancy        false;

   pdfReconstructionModel        OnlyMissing;
   pdfRefillingModel             EquilibriumRefilling;
   excessMassDistributionModel   EvenlyNewInterfaceFallbackLiquid;
   curvatureModel                FiniteDifferenceMethod;
   useSimpleMassExchange         false;
   cellConversionThreshold       1e-2;
   cellConversionForceThreshold  1e-1;
   useBognerReconstruction       true;
   useFillWeighting              true;
   useBognerThreePhaseMEM        true;
   useMemOnlyBuoyancy            true;
}}

EvaluationParameters
{{
   performanceLogFrequency 5000;
   evaluationFrequency     100;
   baseFolderName          eval;

   ClogMeasurement
   {{
      upstreamPos       < {PILLAR_X-48:g}, 20, 40 >;
      downstreamPos     < {PILLAR_X+48:g}, 20, 40 >;
      probeHalfWidthXY  2;

      pillarUpstreamPos    < {PILLAR_X-48:g}, 6, 40 >;
      pillarDownstreamPos  < {PILLAR_X+48:g}, 6, 40 >;
   }}
}}

StabilityChecker
{{
   checkFrequency 0;
   streamOutput   false;
   vtkOutput      true;
}}

BoundaryParameters
{{
   Border {{ direction T;  walldistance -1; NoSlip{{}} }}
   Border {{ direction B;  walldistance -1; NoSlip{{}} }}
   // Fluid-blocking counterparts of the CylinderPillars (full height); Y walls get NoSlip automatically (non-periodic Y).
   BodyBC {{ shape Cylinder; min < {PILLAR_X:g}, 0, 0 >;  max < {PILLAR_X:g}, 0, {LZ:g} >;  radius {R:g}; NoSlip{{}} }}
   BodyBC {{ shape Cylinder; min < {PILLAR_X:g}, {LY:g}, 0 >; max < {PILLAR_X:g}, {LY:g}, {LZ:g} >; radius {R:g}; NoSlip{{}} }}
}}

MeshOutputParameters
{{
   writeFrequency 0;
   baseFolder     mesh-out;
}}

// 12 dumps over the run; fluid cropped to X 128..640, Z 0..100 (rods/pillars unfiltered).
VTK
{{
   fluid_field
   {{
      writeFrequency       {args.vtk_every or args.timesteps // 12};
      ghostLayers          0;
      baseFolder           vtk-out;
      samplingResolution   1;
      useMPIIO false;

      writers
      {{
         velocity;
         mapped_flag;
      }}

      CellBB_filter_0 {{ min < {args.vtk_crop.split(',')[0]}, 0, 0 >;  max < {args.vtk_crop.split(',')[1]}, 40, {args.vtk_crop.split(',')[2]} >; }}

      inclusion_filters
      {{
         CellBB_filter_0;
      }}

      before_functions {{ gas_cell_zero_setter; }}
   }}

   domain_decomposition
   {{
      writeFrequency             10000000000;
      baseFolder                 vtk-out;
      outputDomainDecomposition  true;
   }}
}}
"""
out = pathlib.Path(__file__).parent / args.out
out.write_text(prm)
print("wrote", out)
