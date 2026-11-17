#!/usr/bin/env python3
"""Generate the "hopper" clogging scenario prm (THESIS_PLAN.md item 3, second clogging setup).

Two cylindrical pillars whose axes lie ON the channel's Y walls (Y=0 and Y=Ly), i.e. the plan view
    |) (|
leaving one central gap of width  Ly - 2R.  A large population of small spheres, each of which fits
through the gap alone (d < gap), is pre-placed in the water on a regular grid upstream of the
constriction and released together, to see whether they arch across the gap like grains in a hopper.

Everything except the sphere grid/pillars/probes is copied from the extreme-clog family
(clog_measurement_validation_extreme_pidon_2node_gas.prm): same domain <768,40,180>, same
avgParticleDiameterOverride/dxSIOverride water height (absoluteLiquidHeight = 49.13), same PID, same
ModelParameters. Run:  python3 gen_hopper_clog_prm.py   (writes prms/hopper_clog_pidon_2node_gas.prm)
"""
import argparse
import itertools
import math
import pathlib

ap = argparse.ArgumentParser()
ap.add_argument("--radius", type=float, default=3.0, help="sphere radius (LU)")
ap.add_argument("--pillar-radius", type=float, default=12.0, help="pillar radius (LU); gap = Ly - 2R")
ap.add_argument("--pitch", type=float, default=7.0, help="X row pitch (LU); must exceed 2*radius + --clearance")
ap.add_argument("--ncols", type=int, default=5, help="columns across Y (evenly spread between the walls)")
ap.add_argument("--nlayers", type=int, default=7, help="layers in Z (evenly spread floor -> waterline)")
ap.add_argument("--clearance", type=float, default=1.0, help="minimum surface gap sphere<->sphere/wall/floor (LU)")
ap.add_argument("--rows", type=int, default=29, help="number of X-rows of the sphere block")
ap.add_argument("--front-x", type=float, default=290.0, help="X of the row nearest the pillars")
ap.add_argument("--hold", type=int, default=2000, help="settlingSteps: all spheres held FIXED until this step")
ap.add_argument("--timesteps", type=int, default=60000)
ap.add_argument("--out", default="prms/hopper_clog_pidon_2node_gas.prm")
ap.add_argument("--pillar-clearance", type=float, default=1.0, help="minimum sphere-surface to pillar-surface gap (LU)")
args = ap.parse_args()

LX, LY, LZ = 768.0, 40.0, 180.0
PILLAR_X = LX / 2
WATER_H = 8.001943372119138 * (2.0 + 4.1393)  # absoluteLiquidHeight = 49.13 (see extreme-clog prm header)
r, R, pitch = args.radius, args.pillar_radius, args.pitch
gap = LY - 2 * R
assert pitch > 2 * r, "grid pitch must leave a clearance between neighbouring spheres"
assert gap > 2 * r, "gap must be wider than one sphere or nothing can ever pass"

# --- sphere grid ---------------------------------------------------------------------------------
# Y and Z pitches are derived (not free): columns are spread evenly between the two walls and layers
# evenly from just above the floor to the waterline, each keeping >= --clearance to the boundary.
c = args.clearance
ncols, nz = args.ncols, args.nlayers
pitch_y = (LY - 2 * (r + c)) / (ncols - 1)
pitch_z = (WATER_H - (r + c)) / (nz - 1)  # bottom layer at r+c above the floor, top layer AT the waterline
ys = [r + c + i * pitch_y for i in range(ncols)]
zs = [r + c + k * pitch_z for k in range(nz)]
xs = [args.front_x - i * pitch for i in range(args.rows)]
assert pitch >= 2 * r + c and pitch_y >= 2 * r + c and pitch_z >= 2 * r + c, (pitch, pitch_y, pitch_z)

# Density stratified by height so every density starts near where buoyancy would put it (heavy low,
# neutral in the middle, light at the surface) -- avoids a violent buoyancy-driven re-sort at release.
# Bottom third 1.2, middle third 1.0, top two layers (incl. the surface layer) 0.5.
def density(k):
    if k >= nz - 2:
        return 0.5
    return 1.2 if k < (nz - 2) // 2 else 1.0

spheres = [(x, y, z, density(k)) for x in xs for y in ys for k, z in enumerate(zs)]

# --- geometry checks -----------------------------------------------------------------------------
minsep = min(math.dist(a[:3], b[:3]) for a, b in itertools.combinations(spheres, 2))
assert minsep >= 2 * r + c - 1e-9, f"spheres overlap / too tight: min centre distance {minsep}"
assert all(y - r >= c - 1e-9 and LY - y - r >= c - 1e-9 for _, y, _, _ in spheres), "sphere too close to a Y wall"
assert all(z + r <= WATER_H + r for _, _, z, _ in spheres)
clear = (PILLAR_X - R) - (spheres[0][0] + r)
assert min(s[0] for s in spheres) - r > 0 and max(s[0] for s in spheres) + r < PILLAR_X - R, "block not upstream of pillars"
# exact surface-to-surface clearance to each pillar cylinder (axes at (PILLAR_X, 0) and (PILLAR_X, LY), radius R): centre distance in XY minus R minus r
pillar_gap = min(math.hypot(x - PILLAR_X, y - yc) - R - r for x, y, _, _ in spheres for yc in (0.0, LY))
assert pillar_gap >= args.pillar_clearance - 1e-9, f"sphere too close to a pillar: surface gap {pillar_gap:.3f} < {args.pillar_clearance}"
counts = {d: sum(1 for s in spheres if s[3] == d) for d in (0.5, 1.0, 1.2)}
print(f"{len(spheres)} spheres  ({ncols} cols x {nz} layers x {args.rows} rows)  densities {counts}")
print(f"pitch X/Y/Z = {pitch:g}/{pitch_y:g}/{pitch_z:.3f}  (clearance >= {c:g})")
print(f"radius {r} (d={2*r}), pillar radius {R}, gap {gap} -> gap/d = {gap/(2*r):.2f};  min sphere separation {minsep:.3f}")
print(f"nearest sphere-surface to pillar-surface gap {pillar_gap:.2f} LU")
print(f"block X [{xs[-1]-r:.1f},{xs[0]+r:.1f}], free run to pillar surface {clear:.1f} LU;  Z layers {[round(z,2) for z in zs]}")

sphere_lines = "\n".join(
    f"   FloatingSphere {{ position < {x:.2f}, {y:.2f}, {z:.2f} >; radius {r:g}; densityRatio {d:g}; settlingSteps {args.hold}; }}"
    for x, y, z, d in spheres
)

prm = f"""// {pathlib.Path(args.out).name}   (GENERATED by gen_hopper_clog_prm.py -- edit the script, not this file)
// ─────────────────────────────────────────────────────────────────────────────
// [Claude, 2026-09-18] THESIS ITEM 3 (THESIS_PLAN.md): second clogging setup, "hopper" variant.
//
// GEOMETRY  (plan view, flow left -> right, Y up):
//
//        Y=40 ─────)   (─────        two cylindrical pillars, axes at X={PILLAR_X:g}, on the Y walls
//                   ) . (             (Y=0 and Y=40), radius {R:g}  ->  central gap {gap:g} LU wide
//        Y=0  ─────)   (─────        i.e. "|) (|" seen from above; each pillar is a half-cylinder in
//                                    the channel. Pillars span the full height (Z 0..{LZ:g}).
// The previous extreme-clog family used spheres BIGGER than the slits (a plugged aperture); here
// each sphere (radius {r:g}, d={2*r:g}) fits through the {gap:g} LU gap on its own (gap/d = {gap/(2*r):.2f}), so any
// stoppage has to come from ARCHING -- several spheres wedging each other across the opening, as grains
// do in a silo hopper -- not from a single oversized obstacle.
//
// POPULATION: {len(spheres)} spheres ({ncols} across Y x {nz} layers in Z x {args.rows} rows in X, grid pitch X/Y/Z = {pitch:g}/{pitch_y:g}/{pitch_z:.2f} LU,
// >= {c:g} LU surface clearance everywhere), densityRatio {{0.5,1.0,1.2}} = {counts[0.5]}/{counts[1.0]}/{counts[1.2]} (1.5 dropped as requested). ALL are created up front, in the
// water and at the surface, upstream of the pillars (block X in [{xs[-1]-r:.1f},{xs[0]+r:.1f}], {clear:.1f} LU free run to the pillar
// surface) -- no staggered dropping. Density is stratified by height so every sphere starts near its own
// buoyant equilibrium (1.2 low, 1.0 mid-depth, 0.5 in the top two layers -- the top layer sits exactly
// at the waterline Z={WATER_H:.2f}, half submerged) instead of triggering a violent buoyancy re-sort at release.
// settlingSteps {args.hold} (no startInactive): every sphere is held FIXED for the first {args.hold} steps so
// the free surface can form its contact lines around the surface layer (the same reason the surface-start
// floater scenarios hold -- see FloatingSphereParams::settlingSteps), then all are released at once. Clearance
// >= {c:g} LU means no contact on the first sub-cycle, so the [RELEASE-OVERLAP] check cannot fire.
//
// FLOW: PID ON at the extreme-clog family's target 0.04 (fixed discharge). Domain, water height, PID gains,
// ModelParameters, Z=180 headroom and the 144-rank/2-node decomposition are copied unchanged from
// clog_measurement_validation_extreme_pidon_2node_gas.prm so results stay comparable to that family.
// Downstream leg wraps periodically to X=0 and back upstream of the block; the block itself is kept
// well clear of the X=0/{LX:g} seam.
//
// OPEN CHOICES worth revisiting after the first run: gap/d (change --pillar-radius/--radius), hold vs.
// immediate release (--hold 0), PID off with the calibrated force instead, a no-clog reference with a
// wide gap. Time estimate for {args.timesteps} steps: ~2-3 h on 144 ranks (1000 spheres sit on only ~12
// ranks -- DEM load is unbalanced but small; unverified until the first run).
// ─────────────────────────────────────────────────────────────────────────────

BlockForestParameters
{{
   cellsPerBlock                 < 64, 20, 30 >;
   periodicity                   < 1, 0, 0 >;
   loadSnapshot                  false;
   storeSnapshot                 true;
   snapshotFrequency             5000;
   snapshotBaseFolder            snapshots;
}}

DomainParameters
{{
   domainSize         <768, 40, 180>;
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
   // axes ON the Y walls: half of each cylinder lies outside the domain (harmless -- only the part
   // inside the channel can touch a sphere). Z is an arbitrary reference point on the vertical axis.
   CylinderPillar {{ center < {PILLAR_X:g}, 0, 45 >; radius {R:g}; }}
   CylinderPillar {{ center < {PILLAR_X:g}, {LY:g}, 45 >; radius {R:g}; }}
}}

FloatingSpheres
{{
{sphere_lines}
}}

ModelParameters
{{
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
      // Pair 1 on the gap centreline (Y=20), pair 2 off-axis near the Y=6 side (upstream of the pillar's
      // reach, Y<R); both 48 LU either side of the pillar axis -- upstream probes sit in the free run
      // between the sphere block ({xs[0]+r:.0f}) and the pillar surface ({PILLAR_X-R:g}).
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
   // Fluid-blocking counterparts of the CylinderPillars above (CylinderPillars alone only give the DEM
   // contact; the LBM only sees BodyBC) -- both span the FULL domain height. Y walls get their NoSlip
   // automatically (non-periodic Y, see setupBoundaryHandling).
   BodyBC {{ shape Cylinder; min < {PILLAR_X:g}, 0, 0 >;  max < {PILLAR_X:g}, 0, {LZ:g} >;  radius {R:g}; NoSlip{{}} }}
   BodyBC {{ shape Cylinder; min < {PILLAR_X:g}, {LY:g}, 0 >; max < {PILLAR_X:g}, {LY:g}, {LZ:g} >; radius {R:g}; NoSlip{{}} }}
}}

MeshOutputParameters
{{
   writeFrequency 0;
   baseFolder     mesh-out;
}}

// 12 dumps over the run; the fluid field is cropped to the region of interest (X 128..640, Z 0..100)
// to keep each frame small -- sphere/pillar VTK share this writeFrequency and are unfiltered.
VTK
{{
   fluid_field
   {{
      writeFrequency       {args.timesteps // 12};
      ghostLayers          0;
      baseFolder           vtk-out;
      samplingResolution   1;
      useMPIIO false;

      writers
      {{
         velocity;
         mapped_flag;
      }}

      CellBB_filter_0 {{ min < 128, 0, 0 >;  max < 640, 40, 100 >; }}

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
