#!/usr/bin/env python3
"""Generate the thesis "mixed sediment" showcase prm: 100 free floaters (spheres, boxes, cylinders) of three
densities dropped into still water in a CUBE-like periodic domain <120,120,120>.

    python3 gen_mixed_sediment_prm.py     ->  prms/mixed_sediment_showcase_100floaters_gas.prm

Design:
  * 10 x 10 grid, pitch 12, ONE floater per (X,Y) column, so no floater starts above another (no head-on
    free-fall collision between bodies at first water entry).
  * two drop heights by CHECKERBOARD (i+j)%2: the domain is periodic with 10 (even) columns per side, so a
    2-colouring is valid across the wrap; a 3-level (i-j)%3 scheme is not (10 % 3 != 0 -> equal levels meet at
    the seam, where two 14-LU-long boxes on a 12-LU pitch would start overlapping).
  * shape and density: 9 (shape,density) combos, 11 copies each + 1, shuffled with a FIXED seed, so drop
    height is statistically independent of density and shape.
  * shape ASPECT RATIOS as in the equilibrium runs: box 1:2:1 (<15,30,15>), cylinder L/D = 1.25 (r=24, L=60, axis Y); densities 0.5 / 1.0 / 1.5 (float / neutral / sink).
  * contact model: analytical Box-Box SAT + hybrid analytical/GJK Cylinder-Cylinder, owner-rank tiebreak on.
"""
import argparse
import collections
import itertools
import pathlib
import random

ap = argparse.ArgumentParser()
ap.add_argument("--g", type=float, default=9.81, help="gravity in m/s^2 (SI)")
ARGS = ap.parse_args()

BLOCK = (20, 30, 40)                     # cellsPerBlock -> 6 x 4 x 3 = 72 blocks = all 72 cores of one node (domain stays a 120^3 cube)
N, PITCH, OFF = 10, 12.0, 2.0          # columns per side, pitch, offset of column 0 (best X-seam clearance for 20-LU blocks is 2)
L = N * PITCH                            # 120: periodic footprint
LZ = 120.0
Z_LEVELS = (64.0, 76.0)                  # centre heights; waterline is 56.02 (see WATER_H); Z block seams are at 40/80
WATER_H = 8.001943372119138 * (2.862 + 4.1393)   # liquidHeightFactor 2.862 + floorHeightFactor 4.1393
SEED = 1
# LBM gravity from SI g at the fixed scale dx_SI (dxSIOverride) and reference viscosities: g_LU = g_SI * dt^2/dx with
# dt = nu_LU/nu_SI * dx^2  ->  g_LU = g_SI * (nu_LU/nu_SI)^2 * dx^3.  (9.81 reproduces the 2.13715e-05 used elsewhere.)
DX_SI, NU_LU, NU_SI = 0.000362694127490625, 0.000214026, 1.0016e-6
G_SI = ARGS.g
G_LU = G_SI * (NU_LU / NU_SI) ** 2 * DX_SI ** 3
if abs(G_SI - 9.81) < 1e-9:
    assert abs(G_LU / 2.13715005840044e-05 - 1) < 1e-4, G_LU
TIMESTEPS = 100000
VTK_EVERY = 10000
DENS = (0.5, 1.0, 1.5)
SHAPES = ("sphere", "box", "cylinder")
# axis-aligned half-extents (X, Y, Z) at the initial (identity) orientation
CYL_R, CYL_L = 4.5, 11.25   # L/D = 1.25, the equilibrium cylinder (r=24, L=60, axis Y)
HALF = {"sphere": (4.5, 4.5, 4.5), "box": (3.5, 7.0, 3.5), "cylinder": (CYL_R, CYL_L / 2, CYL_R)}   # box 1:2:1 = the equilibrium cuboid <15,30,15>
BOUND_R = {"sphere": 4.5, "box": (3.5**2 + 7**2 + 3.5**2) ** 0.5, "cylinder": (CYL_R**2 + (CYL_L / 2)**2) ** 0.5}

# ---- assignment ------------------------------------------------------------------------------------
combos = [(s, d) for s in SHAPES for d in DENS]
pool = [combos[k % 9] for k in range(N * N)]          # 11 x 9 + 1 -> one combo appears 12 times
random.Random(SEED).shuffle(pool)
cells = []
for i in range(N):
    for j in range(N):
        s, d = pool[i * N + j]
        cells.append((i, j, s, d, Z_LEVELS[(i + j) % 2]))
pos = {c[:2]: (OFF + PITCH * c[0], OFF + PITCH * c[1], c[4]) for c in cells}

# ---- geometry checks (periodic in X and Y) ---------------------------------------------------------
def mind(a, b, L_):
    d = abs(a - b) % L_
    return min(d, L_ - d)

worst_aabb, worst_bound = 1e9, 1e9
for a, b in itertools.combinations(cells, 2):
    (xa, ya, za), (xb, yb, zb) = pos[a[:2]], pos[b[:2]]
    dx, dy, dz = mind(xa, xb, L), mind(ya, yb, L), abs(za - zb)
    ha, hb = HALF[a[2]], HALF[b[2]]
    gaps = (dx - ha[0] - hb[0], dy - ha[1] - hb[1], dz - ha[2] - hb[2])
    worst_aabb = min(worst_aabb, max(gaps))          # separated iff at least one axis gap > 0
    r3 = (dx * dx + dy * dy + dz * dz) ** 0.5
    worst_bound = min(worst_bound, r3 - BOUND_R[a[2]] - BOUND_R[b[2]])
assert worst_aabb >= 1.0, f"initial AABBs overlap / too tight ({worst_aabb})"
for i, j, s, d, z in cells:
    assert z - HALF[s][2] > WATER_H + 1.0, "floater starts in the water"
    assert z + HALF[s][2] < LZ - 10, "too close to the ceiling"
def _clear(v, b):
    return min(v % b, b - v % b)
seam_x = min(_clear(p[0], BLOCK[0]) for p in pos.values())
seam_y = min(_clear(p[1], BLOCK[1]) for p in pos.values())
seam_z = min(_clear(z, BLOCK[2]) for z in Z_LEVELS)
seam_clear = min(seam_x, seam_y, seam_z)
assert 120 % BLOCK[0] == 0 and 120 % BLOCK[1] == 0 and 120 % BLOCK[2] == 0
NBLK = (120 // BLOCK[0], 120 // BLOCK[1], 120 // BLOCK[2])
NRANKS = NBLK[0] * NBLK[1] * NBLK[2]
assert NBLK[0] >= 3 and NBLK[1] >= 3, "need >= 3 blocks per periodic direction"
assert NRANKS == 72
assert len({p[:2] for p in pos.values()}) == 100

cnt_sd = collections.Counter((c[2], c[3]) for c in cells)
cnt_dz = collections.Counter((c[3], c[4]) for c in cells)
print("shapes :", dict(collections.Counter(c[2] for c in cells)))
print("density:", dict(collections.Counter(c[3] for c in cells)))
print("shape x density:", {f"{s[:3]}/{d}": n for (s, d), n in sorted(cnt_sd.items())})
print("density x drop height:", {f"{d}@{int(z)}": n for (d, z), n in sorted(cnt_dz.items())})
print(f"g_SI = {G_SI} m/s^2 -> gravityOverride = {G_LU:.11e} LU (bound check: 120*g = {120*G_LU:.2e} << 1/3)")
print(f"waterline {WATER_H:.2f}; lowest floater bottom {min(c[4]-HALF[c[2]][2] for c in cells):.2f}")
print(f"min AABB axis-gap between any pair (best axis) = {worst_aabb:.2f}; min bounding-sphere gap = {worst_bound:.2f} (<0 is fine, AABB is the test)")
print(f"blocks {BLOCK} -> {NBLK} = {NRANKS} ranks; min clearance of a floater start centre to a block seam: X {seam_x:g}, Y {seam_y:g}, Z {seam_z:g}")

# ---- emit ------------------------------------------------------------------------------------------
def line(c):
    i, j, s, d, z = c
    x, y, _ = pos[(i, j)]
    if s == "sphere":
        return f"   FloatingSphere {{ position < {x:g}, {y:g}, {z:g} >; radius 4.5; densityRatio {d:g}; settlingSteps 0; }}"
    if s == "box":
        return f"   FloatingBox {{ position < {x:g}, {y:g}, {z:g} >; halfExtents < 3.5, 7, 3.5 >; densityRatio {d:g}; settlingSteps 0; }}"
    return f"   FloatingCylinder {{ position < {x:g}, {y:g}, {z:g} >; radius {CYL_R:g}; length {CYL_L:g}; axis < 0, 1, 0 >; densityRatio {d:g}; settlingSteps 0; }}"

blocks = {}
for c in cells:
    blocks.setdefault(c[2], []).append(line(c))
sec = lambda name, s: f"{name}\n{{\n" + "\n".join(blocks[s]) + "\n}\n"

n_sph, n_box, n_cyl = (len(blocks[s]) for s in SHAPES)
prm = f"""// mixed_sediment_showcase_100floaters_gas.prm   (GENERATED by gen_mixed_sediment_prm.py -- edit the script, not this file)
// ─────────────────────────────────────────────────────────────────────────────
// [Claude, 2026-09-18] THESIS SHOWCASE: MIXED SEDIMENT -- 100 free floaters of all three shapes ({n_sph} spheres,
// {n_box} boxes, {n_cyl} cylinders) and three densities (0.5 floats, 1.0 neutral, 1.5 sinks -- Bogner-showcase densities, per request; every shape appears at every density) dropped from rest into
// STILL WATER (targetMeanVelocityMagnitude 0, no forcing) in a CUBE-LIKE domain <120,120,120> (X,Y periodic, floor and
// ceiling walls), g = {G_SI:g} m/s^2. NOT YET RUN -- the first cluster run is the real test.
//
// CONTACT MODEL: useAnalyticalContactDetection true -> Box x Box by the closed-form OBBTree SAT, Cylinder x Cylinder by the
// analytical Kodam+Guo HYBRID (usePureAnalyticalCylinder false: pairs the analytic dispatch flags uncertain fall back to
// GJK/EPA); useContactFilterOwnerTiebreak true -> a contact straddling a rank boundary is accepted by exactly one rank
// (min owner rank), not by two independent geometric tests.
//
// UNITS/GRAVITY: dx_SI = 0.000362694 m (avgParticleDiameterOverride 8.0019 / dxSIOverride -- spheres_out.dat is never opened and
// X=120 does not rescale dx). gravityOverride = g_SI * (nu_LU/nu_SI)^2 * dx^3 = {G_LU:.11e} LU (9.81 m/s^2 gives 2.13715e-05, the
// value used by the other real-gravity showcases). {LZ:g} cells * g = {LZ*G_LU:.1e} << cs^2 = 1/3, so the LBM incompressibility
// bound (Buick & Greated 2000, Bogner Eq. 10) holds with a wide margin. Water depth: liquidHeightFactor 2.862 + floorHeightFactor
// 4.1393 -> waterline Z = {WATER_H:.2f} LU (20.3 mm), ~47% of the cube; the rest is fall height and headroom.
// ModelParameters otherwise the current project set (Bogner three-phase MEM, MEM-only buoyancy, fill weighting).
//
// SHAPES -- aspect ratios taken from the equilibrium (buoyancy) runs: box halfExtents <3.5,7,3.5> = 7x14x7 = 1:2:1, as the equilibrium cuboid
// <15,30,15> (2s x s x s, long axis Y); cylinder radius {CYL_R:g}, length {CYL_L:g}, axis Y -> L/D = {CYL_L/(2*CYL_R):g}, as the equilibrium cylinder
// (r=24, L=60, axis Y -> L/D = 1.25); sphere r=4.5. Volumes: box {7*14*7} , cylinder {3.141592653589793*CYL_R**2*CYL_L:.0f}, sphere {4/3*3.141592653589793*4.5**3:.0f} LU^3. Free to translate and
// rotate, no tilt, no pinning.
//
// LAYOUT -- 10 x 10 grid, pitch 12, ONE floater per (X,Y) column (nothing starts above anything else), column centres at {OFF:g} + 12 i
// (start centres >= {seam_x:g} LU (X) / {seam_y:g} LU (Y) / {seam_z:g} LU (Z) from the {BLOCK} block seams):
//  * two drop heights, Z = {Z_LEVELS[0]:g} / {Z_LEVELS[1]:g}, by CHECKERBOARD (i+j)%2. 10 columns per side is even, so the checkerboard is a
//    valid 2-colouring ACROSS the periodic wrap; (i-j)%3 is not (10 % 3 != 0), and equal levels would meet at the seam, where two
//    14-LU-long boxes on a 12-LU pitch would start overlapping. Lowest floater bottom at Z = {min(c[4]-HALF[c[2]][2] for c in cells):.2f} (waterline {WATER_H:.2f}).
//  * shape and density: the 9 (shape,density) combinations, 11 copies each plus one extra, shuffled with fixed seed {SEED}, so drop
//    height is independent of density and shape. Counts: shapes {n_sph}/{n_box}/{n_cyl}; densities
//    {'/'.join(str(sum(1 for c in cells if c[3]==d)) for d in DENS)} for 0.5/1.0/1.5.
//  * initial AABBs of all 100 floaters (X,Y periodic) are disjoint with >= {worst_aabb:.1f} LU to spare (checked in the generator, all pairs).
//
// EXPECTED (unverified): 0.5 ends half-submerged at the surface, 1.5 sinks to the floor (terminal velocity ~0.016 LU/step, a few
// thousand steps), 1.0 hangs wherever it comes to rest; shape/density interactions where they meet at the surface and on the floor.
//
// timesteps {TIMESTEPS}, VTK every {VTK_EVERY} steps = dumps at 0,{VTK_EVERY},...,{TIMESTEPS-VTK_EVERY} (the loop is `t != timesteps`, so the final step {TIMESTEPS} itself is not written).
// VTK options as in the extreme-clog prms: before_functions {{ gas_cell_zero_setter; }} (without it gas/obstacle cells show stale
// PDF-derived noise), writers velocity + mapped_flag, plus fill_level for the free surface. Particle output (spheres/boxes/cylinders)
// shares this writeFrequency; the fluid is two mid-plane slices, NOT the full volume.
// Domain <120,120,120> / cellsPerBlock <{BLOCK[0]},{BLOCK[1]},{BLOCK[2]}> = {NBLK[0]}x{NBLK[1]}x{NBLK[2]} = {NRANKS} processes = every core of one singlenode node (the domain stays a cube; only the blocks are not). Snapshots every 10000.
// ─────────────────────────────────────────────────────────────────────────────

BlockForestParameters
{{
   cellsPerBlock                 < {BLOCK[0]}, {BLOCK[1]}, {BLOCK[2]} >;
   periodicity                   < 1, 1, 0 >;   // X and Y periodic -- no side walls, only Z floor/ceiling
   loadSnapshot                  false;
   storeSnapshot                 true;
   snapshotFrequency             10000;
   snapshotBaseFolder            snapshots;
}}

DomainParameters
{{
   domainSize         <120, 120, 120>;
   wavePeriods        1;
   liquidHeightFactor 2.862;
   floorHeightFactor  4.1393;
   initialAmplitude   0;
}}

PIDParameters
{{
   targetMeanVelocityMagnitude 0;       // still water -- no flow forcing
   proportionalGain            2e-4;
   derivativeGain              1e-6;
   integralGain                2e-4;
   maxRamp                     1e-4;
   minActuatingVariable        0;
   maxActuatingVariable        1e-3;
}}

PhysicsParameters
{{
   enableWetting     false;
   timesteps         {TIMESTEPS};
   Re                3100;
   Fr                1.31;
   We                15.6188;
   useHardcodedUnits true;              // U = 0: keep the fixed reference viscosity/relaxation rate
   gravityOverride   {G_LU:.14e};   // g_SI = {G_SI:g} m/s^2 at dx_SI = 0.000362694 m (see header)
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
   restitutionCoefficient 0.97;
   numSubCycles 10;
   useLubricationCorrection true;
   useNoSlipParticles false;
   useParticleBed false;
}}

{sec("FloatingSpheres", "sphere")}
{sec("FloatingBoxes", "box")}
{sec("FloatingCylinders", "cylinder")}
ModelParameters
{{
   useAnalyticalContactDetection true;   // Box x Box: closed-form SAT; Cylinder x Cylinder: analytical hybrid
   usePureAnalyticalCylinder     false;  // (default, explicit) keep the GJK/EPA fallback for uncertain cylinder pairs
   useContactFilterOwnerTiebreak true;   // one rank accepts each boundary-straddling contact (min owner rank)

   // Experimental contact/buoyancy fixes: ALL explicitly OFF (they default to false too -- written out so a later
   // default change cannot silently switch one on in this thesis showcase).
   useTightTier1                 false;
   useTightLayer1                false;
   skipTier1                     false;
   useRoundedBoxCorners          false;
   useMultiPointBoxContact       false;
   useBoxOverlapCap              false;
   useAxisHysteresis             false;
   useContactManifold            false;
   useRotationalEffectiveMass    false;
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
   evaluationFrequency     500;
   baseFolderName          eval;
}}

StabilityChecker
{{
   checkFrequency 0;
   streamOutput   false;
   vtkOutput      true;
}}

BoundaryParameters
{{
   // Only Z (top/bottom) walls -- X and Y are periodic, setupBoundaryHandling skips periodic faces.
   Border {{ direction T;  walldistance -1; NoSlip{{}} }}
   Border {{ direction B;  walldistance -1; NoSlip{{}} }}
}}

MeshOutputParameters
{{
   writeFrequency 0;
   baseFolder     mesh-out;
}}

VTK
{{
   fluid_field
   {{
      writeFrequency       {VTK_EVERY};
      ghostLayers          0;
      baseFolder           vtk-out;
      samplingResolution   1;
      useMPIIO false;

      writers
      {{
         velocity;
         fill_level;
         mapped_flag;  // cell type per timestep: 1=liquid 2=interface 3=gas 4=noSlip (VtkWriter.h)
      }}

      // two mid-plane slices (Y=60 and X=60), full height -- not the whole volume
      CellBB_filter_0 {{ min < 0, 60, 0 >;  max < 120, 61, 120 >; }}
      CellBB_filter_1 {{ min < 60, 0, 0 >;  max < 61, 120, 120 >; }}

      inclusion_filters
      {{
         CellBB_filter_0;
         CellBB_filter_1;
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
out = pathlib.Path(__file__).parent / "prms" / "mixed_sediment_showcase_100floaters_gas.prm"
out.write_text(prm)
print("wrote", out)
