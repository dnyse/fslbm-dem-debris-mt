#!/usr/bin/env python3
"""[Claude, 2026-09-25] "DROP WAVES" variant of the half-density cylinder hopper (companion of
gen_hopper_clog_cyl_prm.py, which pre-places all rods in a block against the pillars, held FIXED,
then releases them all at once). Per explicit request: instead of pre-placing, release the SAME
210-rod population (--rows 15, 2 across Y x 7 layers x 15 rows -- identical geometry/densities to
the existing halfrows scenario) as ~6 waves of 36 rods each (6 x 2 x 3), dropped from ABOVE the
water surface at the UPSTREAM end of the channel ("top of the stream"), similar to the staggered
startInactive+settlingSteps release ladder already used by the sphere clogging scenario
(clog_measurement_validation_extreme_pidon_gas.prm). Waves fall through the air, land in the flow
far upstream of the pillars, and get carried downstream by the PID-driven current -- clogging (if
any) has to develop through flow-driven transport and pile-up, not from a pre-arranged block.

    python3 gen_hopper_clog_cyl_dropwaves_prm.py   ->  prms/hopper_clog_cyl_dropwaves_gas.prm

Reuses gen_hopper_clog_cyl_prm.py's exact grid layout/validation logic for the 210 rod
X/Y positions and per-layer densities; only the Z placement (above-surface drop height instead of
submerged depth) and the release timing (staggered waves instead of "hold then release together")
differ.
"""
import argparse
import itertools
import math
import pathlib

ap = argparse.ArgumentParser()
ap.add_argument("--radius", type=float, default=3.0)
ap.add_argument("--length", type=float, default=18.0)
ap.add_argument("--pillar-radius", type=float, default=12.0)
ap.add_argument("--pitch", type=float, default=7.0)
ap.add_argument("--nlayers", type=int, default=7)
ap.add_argument("--clearance", type=float, default=1.0)
ap.add_argument("--rows", type=int, default=15)                 # same population as the halfrows scenario (210 rods)
ap.add_argument("--front-x", type=float, default=110.0, help="X of the row nearest the drop zone's downstream edge -- UPSTREAM this time, not pressed against the pillars")
ap.add_argument("--pillar-clearance", type=float, default=1.0)
ap.add_argument("--wave-size", type=int, default=36, help="rods released together per wave (6 x 2 x 3)")
ap.add_argument("--hold", type=int, default=2000, help="steps before the FIRST wave is released (flow spin-up)")
ap.add_argument("--wave-interval", type=int, default=3000, help="steps between successive wave releases")
ap.add_argument("--drop-height-above-surface", type=float, default=6.0, help="lowest drop layer's height above the water surface (LU)")
ap.add_argument("--drop-layer-spacing", type=float, default=7.0, help="vertical spacing between the 7 drop layers (LU)")
ap.add_argument("--timesteps", type=int, default=30000)
ap.add_argument("--out", default="prms/hopper_clog_cyl_dropwaves_gas.prm")
ap.add_argument("--lx", type=float, default=768.0)
ap.add_argument("--lz", type=float, default=180.0)
ap.add_argument("--cells", default="64,20,30")
ap.add_argument("--vtk-every", type=int, default=0)
ap.add_argument("--vtk-crop", default="0,400,140", help="fluid VTK crop 'x0,x1,zmax' -- upstream + pillars, not the old downstream-only block crop")
args = ap.parse_args()

LX, LY, LZ = args.lx, 40.0, args.lz
PILLAR_X = LX / 2
WATER_H = 8.001943372119138 * (2.0 + 4.1393)   # absoluteLiquidHeight = 49.13, same as the halfrows scenario
r, L, R, pitch, c = args.radius, args.length, args.pillar_radius, args.pitch, args.clearance
gap = LY - 2 * R
assert pitch >= 2 * r + c

# ---- X/Y grid: identical logic to gen_hopper_clog_cyl_prm.py (2 across Y, --rows in X) ------------------------------------------------------------
ncols = int((LY - c) // (L + c))
assert ncols >= 1
pitch_y = (LY - 2 * (L / 2 + c)) / (ncols - 1) if ncols > 1 else 0.0
ys = [L / 2 + c + i * pitch_y for i in range(ncols)] if ncols > 1 else [LY / 2]
nz = args.nlayers

# ---- Z: ABOVE the water surface (drop layers), not submerged depths -- the one real geometry change from the block scenario -----------------------
zs = [WATER_H + args.drop_height_above_surface + k * args.drop_layer_spacing for k in range(nz)]
xs = [args.front_x - i * pitch for i in range(args.rows)]
assert zs[-1] < LZ - r - c, f"top drop layer too close to domain ceiling: {zs[-1]}"

def density(k):   # same per-layer stratification as the block scenario, kept only to preserve population composition (60/90/60)
    if k >= nz - 2: return 0.5
    return 1.2 if k < (nz - 2) // 2 else 1.0

rods = [(x, y, z, density(k)) for x in xs for y in ys for k, z in enumerate(zs)]

# ---- exact geometry checks (rods checked at their DROP positions -- horizontal spacing only matters in the air; they fall independently) ----------
def rod_gap(a, b):
    return max(math.hypot(a[0] - b[0], a[2] - b[2]) - 2 * r, abs(a[1] - b[1]) - L)
minsep = min(rod_gap(a, b) for a, b in itertools.combinations(rods, 2))
assert minsep >= c - 1e-9, f"rods too close at drop time: {minsep}"
assert all(y - L / 2 >= c - 1e-9 and LY - (y + L / 2) >= c - 1e-9 for _, y, _, _ in rods), "rod too close to a Y wall"
counts = {d: sum(1 for q in rods if q[3] == d) for d in (0.5, 1.0, 1.2)}
print(f"{len(rods)} rods (r={r:g}, L={L:g}, axis Y): {ncols} across Y x {nz} drop layers x {args.rows} rows; densities {counts}")
print(f"drop zone X [{xs[-1]-r:.1f},{xs[0]+r:.1f}] (upstream; pillars at X={PILLAR_X:g}, clearance {PILLAR_X-r-xs[0]:.1f} LU)")
print(f"drop layers Z {[round(z,2) for z in zs]} (water surface Z={WATER_H:.2f})")

# ---- wave batching: flatten in the natural (x, y, layer) order and chunk into --wave-size groups ---------------------------------------------------
waves = [rods[i:i + args.wave_size] for i in range(0, len(rods), args.wave_size)]
nwaves = len(waves)
last_release = args.hold + (nwaves - 1) * args.wave_interval
print(f"{nwaves} waves of up to {args.wave_size} rods; releases at steps "
      f"{[args.hold + i * args.wave_interval for i in range(nwaves)]}; last release t={last_release}")
if args.timesteps < last_release + 5000:
    print(f"WARNING: --timesteps {args.timesteps} leaves under 5000 steps after the last wave's release (t={last_release}) "
          f"for it to fall/travel/interact -- consider raising --timesteps.")

lines = []
for wave_idx, wave in enumerate(waves):
    release_step = args.hold + wave_idx * args.wave_interval
    lines.append(f"   // wave {wave_idx}: {len(wave)} rods, released (startInactive -> active) at step {release_step}")
    for x, y, z, d in wave:
        lines.append(
            f"   FloatingCylinder {{ position < {x:.2f}, {y:.2f}, {z:.2f} >; radius {r:g}; length {L:g}; axis < 0, 1, 0 >; "
            f"densityRatio {d:g}; settlingSteps {release_step}; startInactive true; }}")
lines = "\n".join(lines)

prm = f"""// {pathlib.Path(args.out).name}   (GENERATED by gen_hopper_clog_cyl_dropwaves_prm.py -- edit the script, not this file)
// ─────────────────────────────────────────────────────────────────────────────
// [Claude, 2026-09-25] "DROP WAVES" cylinder hopper, per explicit request: same 210-rod population
// as the halfrows block scenario (2 across Y x {nz} layers x {args.rows} rows, same densities), but
// released as {nwaves} waves of up to {args.wave_size} rods (6 x 2 x 3) dropped from ABOVE the water
// surface (drop layers Z={[round(z,1) for z in zs]}, surface Z={WATER_H:.2f}) at the UPSTREAM end of
// the channel (drop zone X=[{xs[-1]-r:.1f},{xs[0]+r:.1f}], {PILLAR_X-r-xs[0]:.1f} LU clear of the
// pillars at X={PILLAR_X:g}) -- similar in spirit to the sphere clogging scenario's staggered
// startInactive+settlingSteps release ladder. Each wave falls through the air, lands in the
// PID-driven flow, and is carried downstream toward the pillars -- any clog has to develop from
// flow-driven transport and pile-up, not a pre-arranged block. Wave releases at steps
// {[args.hold + i * args.wave_interval for i in range(nwaves)]} (hold={args.hold}, interval={args.wave_interval});
// {args.timesteps} total timesteps.
//
// THIS RUN IS EXPLICITLY REQUESTED WITH ALL OF THIS SESSION'S 2026-09-25 CRASH-MITIGATION FIXES
// DISABLED (Option A in StreamReconstructAdvectSweep.h, Option B/repairDegenerateObstacleCells in
// SetupStages.h+BuildTimeloop.h, the SimpleBB.h [MEM-GARBAGE-LINK] active guard) -- a genuine
// pre-session baseline for this NEW release mechanic, not validated against any of today's fixes.
// See THESIS_PLAN.md for the disable-patch details and which mirror this was run on.
//
// Same channel/pillars/PID/physics/contact-model as the halfrows block scenario
// (hopper_clog_cyl_frontblock_pidon_2node_halfrows_gas.prm) otherwise -- domain
// <{LX:g},{LY:g},{LZ:g}>, 144 ranks/2 nodes, PID target 0.04, Re/Fr/We unchanged.
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
   inFileName spheres_out.dat;
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
   CylinderPillar {{ center < {PILLAR_X:g}, 0, 45 >; radius {R:g}; }}
   CylinderPillar {{ center < {PILLAR_X:g}, {LY:g}, 45 >; radius {R:g}; }}
}}

FloatingCylinders
{{
{lines}
}}

ModelParameters
{{
   useAnalyticalContactDetection true;
   usePureAnalyticalCylinder     false;
   useContactFilterOwnerTiebreak true;

   useTightTier1                 false;
   useTightLayer1                true;
   skipTier1                     false;
   useRoundedBoxCorners          false;
   useMultiPointBoxContact       false;
   useBoxOverlapCap              false;
   useAxisHysteresis             false;
   useContactManifold            false;
   useRotationalEffectiveMass    true;
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
   BodyBC {{ shape Cylinder; min < {PILLAR_X:g}, 0, 0 >;  max < {PILLAR_X:g}, 0, {LZ:g} >;  radius {R:g}; NoSlip{{}} }}
   BodyBC {{ shape Cylinder; min < {PILLAR_X:g}, {LY:g}, 0 >; max < {PILLAR_X:g}, {LY:g}, {LZ:g} >; radius {R:g}; NoSlip{{}} }}
}}

MeshOutputParameters
{{
   writeFrequency 0;
   baseFolder     mesh-out;
}}

// fluid cropped to X 0..400 (upstream drop zone + pillars), Z 0..140 (drop layers reach ~97).
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
