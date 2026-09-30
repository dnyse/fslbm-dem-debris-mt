#!/usr/bin/env python3
"""Cuboid (box) version of the "hopper" clogging scenario -- shape-comparison companion of
gen_hopper_clog_cyl_prm.py (cylinders) and gen_hopper_clog_prm.py (spheres).

    python3 gen_hopper_clog_box_prm.py   ->  prms/hopper_clog_box_frontblock_pidon_2node_gas.prm

IDENTICAL channel / pillars / PID / water height / VTK / ModelParameters / population grid
(positions, densities, row/column/layer counts) to the cylinder hopper -- every FloatingCylinder is
replaced 1:1 by a FloatingBox of the SAME bounding footprint (halfExtents <r, length/2, r>, no axis
field needed since an axis-aligned box's long dimension is simply its own Y half-extent), matching
the same shape-comparison convention already used for the smaller 24x/48x doubling-series runs
(boxes_24x_2x4x3_clog_gas.prm's own halfExtents <3,9,3>). This is the FULL production-scale (~406
body) companion to that smaller series: same purpose (does the STALE-PARTICLE-IDX cross-rank
ghost-ownership desync found with cylinders reproduce with box-box SAT contacts instead, at the same
population/layout), but at the scale actually used for the thesis's own hopper-gap clogging scenario
(Table 4.2's "hopper gap" row: domain 768x40x180, 406 floaters), not just the 24/48-body toy series.

Geometry note: since neither shape is ever rotated at setup (both released axis-aligned, long axis
along Y), a box's XY-plan footprint IS the exact same rectangle a cylinder's own bounding box already
occupies -- every geometry check below (Y-wall clearance, floor clearance, pillar clearance via
rect_disc_gap) therefore carries over from the cylinder generator completely unchanged; only the
body-body minimum-separation check needs its own box-specific (AABB) formula in place of the
cylinder's capsule/segment-distance one.
"""
import argparse
import itertools
import math
import pathlib

ap = argparse.ArgumentParser()
ap.add_argument("--radius", type=float, default=3.0, help="box half-extent in X/Z (LU) -- matches the cylinder radius")
ap.add_argument("--length", type=float, default=18.0, help="box full extent along Y (LU); diameter:length = 1:3 by default, as the cylinder")
ap.add_argument("--pillar-radius", type=float, default=12.0)
ap.add_argument("--pitch", type=float, default=7.0, help="X row pitch (LU) >= 2*radius + clearance")
ap.add_argument("--nlayers", type=int, default=7)
ap.add_argument("--clearance", type=float, default=1.0, help="min surface gap box<->box/wall/floor (LU)")
ap.add_argument("--rows", type=int, default=29)
ap.add_argument("--front-x", type=float, default=368.0, help="X of the row nearest the pillars")
ap.add_argument("--pillar-clearance", type=float, default=1.0)
ap.add_argument("--hold", type=int, default=2000)
ap.add_argument("--timesteps", type=int, default=60000)
ap.add_argument("--out", default="prms/hopper_clog_box_frontblock_pidon_2node_gas.prm")
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

# ---- grid: identical to gen_hopper_clog_cyl_prm.py -- as many boxes across Y as fit, evenly spread; layers evenly floor -> waterline ------------------
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

def density(k):                                        # heavy low, neutral mid, light at the surface (as the sphere/cylinder versions)
    if k >= nz - 2: return 0.5
    return 1.2 if k < (nz - 2) // 2 else 1.0

boxes = [(x, y, z, density(k)) for x in xs for y in ys for k, z in enumerate(zs)]

# ---- exact geometry checks (all boxes axis-aligned, half-extents <r, L/2, r>) --------------------------------------------------------------------
def box_gap(a, b):     # AABB surface distance (exact for two axis-aligned, non-rotated boxes)
    dx = max(0.0, abs(a[0] - b[0]) - 2 * r)
    dy = max(0.0, abs(a[1] - b[1]) - L)
    dz = max(0.0, abs(a[2] - b[2]) - 2 * r)
    return math.sqrt(dx * dx + dy * dy + dz * dz)
minsep = min(box_gap(a, b) for a, b in itertools.combinations(boxes, 2))
assert minsep >= c - 1e-9, f"boxes too close: {minsep}"
assert all(y - L / 2 >= c - 1e-9 and LY - (y + L / 2) >= c - 1e-9 for _, y, _, _ in boxes), "box too close to a Y wall"
assert all(z - r >= c - 1e-9 for _, _, z, _ in boxes), "box too close to the floor"
def rect_disc_gap(x, y, xc, yc):   # XY footprint of a box is the SAME rectangle [x-r,x+r]x[y-L/2,y+L/2] a cylinder's bounding box already occupies
    nx = min(max(xc, x - r), x + r); ny = min(max(yc, y - L / 2), y + L / 2)
    return math.hypot(nx - xc, ny - yc) - R
pillar_gap = min(rect_disc_gap(x, y, PILLAR_X, yc) for x, y, _, _ in boxes for yc in (0.0, LY))
assert pillar_gap >= args.pillar_clearance - 1e-9, f"box too close to a pillar: {pillar_gap:.3f}"
counts = {d: sum(1 for q in boxes if q[3] == d) for d in (0.5, 1.0, 1.2)}
vol = (2 * r) * L * (2 * r)
print(f"{len(boxes)} boxes (halfExtents <{r:g},{L/2:g},{r:g}>): {ncols} across Y x {nz} layers x {args.rows} rows; densities {counts}")
print(f"pitch X/Y/Z = {pitch:g}/{pitch_y:g}/{pitch_z:.3f}; min box-box surface gap {minsep:.2f}; nearest box-to-pillar-surface gap {pillar_gap:.3f}")
print(f"block X [{xs[-1]-r:.1f},{xs[0]+r:.1f}]; Y centres {[round(y,2) for y in ys]}; Z layers {[round(z,2) for z in zs]}")
print(f"gap {gap:g} vs box length {L:g}: a box lying along Y (broadside) does NOT fit the gap ({'blocked' if L > gap else 'fits'}); "
      f"box volume {vol:.0f} LU^3 (matches the cylinder's own bounding box exactly); solid fraction of the block ~ {vol/(pitch*(LY/ncols)*pitch_z):.2f}")

lines = "\n".join(
    f"   FloatingBox {{ position < {x:.2f}, {y:.2f}, {z:.2f} >; halfExtents < {r:g}, {L/2:g}, {r:g} >; densityRatio {d:g}; settlingSteps {args.hold}; }}"
    for x, y, z, d in boxes)

prm = f"""// {pathlib.Path(args.out).name}   (GENERATED by gen_hopper_clog_box_prm.py -- edit the script, not this file)
// ─────────────────────────────────────────────────────────────────────────────
// [Claude, 2026-09-30] "HOPPER" CLOGGING SETUP, BOX VERSION -- shape-comparison companion of
// gen_hopper_clog_cyl_prm.py's own cylinder hopper (the crash chased at length in THESIS_PLAN.md,
// section 5 / the doubling-in-X ghost-ownership desync investigation). IDENTICAL channel, pillars,
// PID, population grid (positions/densities/counts) and ModelParameters to the cylinder version;
// every FloatingCylinder is replaced 1:1 by a FloatingBox of the SAME bounding footprint
// (halfExtents <{r:g},{L/2:g},{r:g}>). Two cylindrical pillars centred on the Y walls (X={PILLAR_X:g},
// radius {R:g}, central gap {gap:g} LU) -- but the debris are CUBOIDS, half-extents <{r:g},{L/2:g},{r:g}>
// ({2*r:g}x{L:g}x{2*r:g} full size, the sphere's {2*r:g} LU cross-section), LONG AXIS ALONG Y (across
// the channel), released together from a block pressed against the pillars.
//
// KEY GEOMETRIC FACT: the box's Y-extent {L:g} is {'LONGER' if L > gap else 'shorter'} than the {gap:g} LU gap. A box lying along Y
// therefore {'cannot pass the gap broadside: it has to rotate towards the flow' if L > gap else 'fits the gap broadside'} -- same constraint as the cylinder version.
//
// POPULATION: {len(boxes)} boxes = {ncols} across Y x {nz} layers x {args.rows} rows in X (pitch X/Y/Z = {pitch:g}/{pitch_y:g}/{pitch_z:.2f},
// >= {c:g} LU surface clearance box-box / wall / floor, checked via exact AABB distance). densityRatio {{0.5,1.0,1.2}} =
// {counts[0.5]}/{counts[1.0]}/{counts[1.2]}, stratified by height exactly as the cylinder version (1.2 low, 1.0 mid, 0.5 in the
// top two layers; the top layer's centre sits at the waterline Z={WATER_H:.2f}, half submerged).
// Block X in [{xs[-1]-r:.1f},{xs[0]+r:.1f}], front row {pillar_gap:.2f} LU (surface to pillar surface) from the pillar cylinders.
// All boxes are held FIXED for {args.hold} steps (settlingSteps, no startInactive) and then released together.
//
// CONTACT MODEL: useAnalyticalContactDetection true (Box x Box = SAT, exact for convex polyhedra) +
// useContactFilterOwnerTiebreak true; every experimental contact/buoyancy fix explicitly OFF, matching
// the cylinder hopper's own baseline (useRotationalEffectiveMass true, useTightLayer1 true -- the
// latter is a cylinder-specific pre-filter, a no-op with zero cylinders present, kept only for an
// exact side-by-side diff against the cylinder prm). Friction/restitution as in the cylinder hopper
// (0.5 / 0.5). FLOW/UNITS: identical (PID on, target 0.04, water Z=49.13, Z=180 headroom,
// <768,40,180> / <64,20,30> = 144 ranks / 2 nodes).
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

FloatingBoxes
{{
{lines}
}}

ModelParameters
{{
   useAnalyticalContactDetection true;   // Box x Box: exact SAT
   usePureAnalyticalCylinder     false;  // no-op (no cylinders present); kept for an exact diff against the cylinder prm
   useContactFilterOwnerTiebreak true;   // one rank accepts each boundary-straddling contact (min owner rank)

   // Experimental contact/buoyancy fixes: same baseline as the cylinder hopper
   useTightTier1                 false;
   useTightLayer1                true;   // no-op for Box x Box (cylinder-specific pre-filter); kept for an exact diff
   skipTier1                     false;
   useRoundedBoxCorners          false;
   useMultiPointBoxContact       false;
   useBoxOverlapCap              false;
   useAxisHysteresis             false;
   useContactManifold            false;
   useRotationalEffectiveMass    true;   // same as the cylinder hopper (generalised effective mass incl. rotational inertia)
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

// 12 dumps over the run; fluid cropped to X 128..640, Z 0..100 (boxes/pillars unfiltered).
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
