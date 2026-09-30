#!/usr/bin/env python3
"""Generate the "Bogner exact" floating-cuboid equilibrium prm at a chosen resolution s (box cross-section in cells).

    python3 gen_bogner_exact_prm.py --s 16    ->  prms/bogner_exact_3box_s16_gas.prm   (his own resolution)
    python3 gen_bogner_exact_prm.py --s 32    ->  prms/bogner_exact_3box_s32_gas.prm   (2x finer)

Reproduces Bogner (2017) Sec. 5.3.3 / Bogner & Ruede (2013): three floating boxes 2s x s x s of density 0.75, 0.5, 0.25 (left to right,
Fig. 5.10) in a 130 x 40 x 40 (for s=16) basin, lower half liquid, tau = 1/1.9, g = 7.5e-4, seed heel 2.86 deg, rotation only about the
longitudinal (Y) axis. Expected stable heel: 26.565 / 45 / 26.565 deg. MEM-only buoyancy (project default), NO analytic buoyancy.
"""
import argparse, pathlib

ap = argparse.ArgumentParser()
ap.add_argument("--s", type=int, default=16, help="box cross-section in cells (Bogner: 16)")
ap.add_argument("--timesteps", type=int, default=0, help="default: 50000 * s/16")
ap.add_argument("--blocks", default="", help="cellsPerBlock override 'bx,by,bz'")
ap.add_argument("--seed", type=float, default=2.86, help="initial heel about Y in degrees (Bogner: 2.86)")
ap.add_argument("--nopin", action="store_true", help="fully free boxes (6 DOF) instead of X,Y-pinned roll-only")
ap.add_argument("--omega", type=float, default=1.9, help="relaxationRateOverride, i.e. lattice omega = 1/tau (Bogner: 1.9); LOWER omega = MORE viscous damping")
ap.add_argument("--evalfreq", type=int, default=500, help="EvaluationParameters.evaluationFrequency (also gates [BoxForces]/[BoxTorque] logging); 1 = every timestep")
ap.add_argument("--only", type=float, default=0, help="keep only the box of this density ratio (e.g. 0.75); 0 = all three")
ap.add_argument("--z0", type=float, default=0, help="override the box centre height (default: upright Archimedes draft); for a large --seed use the ideal heave of that heel (rho 0.75 at 26.565 deg: 16.9223, from the righting scans)")
ap.add_argument("--tag", default="", help="suffix for the output file name (default: derived from --seed/--nopin)")
ARGS = ap.parse_args()

s = ARGS.s
k = s / 16.0
SEED = ARGS.seed
TAG = ARGS.tag or (('_seed' + format(SEED, 'g').replace('.', 'p')) if abs(SEED - 2.86) > 1e-9 else '') + ('_nopin' if ARGS.nopin else '')
LX, LY, LZ = int(round(130 * k)), int(round(40 * k)), int(round(40 * k))
W = LZ / 2.0                                   # lower half liquid
G = 7.5e-4 / k                                 # keep H*g at Bogner's own margin (40 * 7.5e-4 = 0.03) when refining
TIMESTEPS = ARGS.timesteps or int(50000 * k)
VTK_EVERY = max(500, TIMESTEPS // 50)
D_AVG = 8.0 * k                                # avgParticleDiameterOverride chosen so that d * (1.25 + 1.25) = W
assert abs(D_AVG * 2.5 - W) < 1e-9
xs = [LX * (i + 0.5) / 3.0 for i in range(3)]  # three equal thirds along X
dens = [0.75, 0.5, 0.25]                        # left to right, Bogner Fig. 5.10
# cellsPerBlock: as many ranks as fit the node(s), subject to X and Y periodic (>= 3 blocks each) and the ghost-sync reach: a box's half-diagonal must be
# reachable in <= 4 sync hops (setup does 4 rounds), conservatively ceil(halfdiag/block)+1 <= 4 in every direction.
DEFAULT_BLOCKS = {16: (10, 8, 40),    # 13x5x1 = 65 ranks (one node; 72 is impossible: 130 and 40 have no factor 3)
                  32: (20, 16, 40)}   # 13x5x2 = 130 ranks (two nodes; 144 impossible for the same reason)
blocks = tuple(int(v) for v in ARGS.blocks.split(",")) if ARGS.blocks else DEFAULT_BLOCKS[s]
assert LX % blocks[0] == 0 and LY % blocks[1] == 0 and LZ % blocks[2] == 0, (LX, LY, LZ, blocks)
nb = (LX // blocks[0], LY // blocks[1], LZ // blocks[2]); nr = nb[0] * nb[1] * nb[2]
assert nb[0] >= 3 and nb[1] >= 3, "X and Y are periodic: need >= 3 blocks in each"
import math
half_diag = math.sqrt((s / 2.0) ** 2 + s ** 2 + (s / 2.0) ** 2)   # box half-diagonal, halfExtents <s/2, s, s/2>
hops = max(math.ceil(half_diag / b) + 1 for b in blocks)
assert hops <= 4, f"box half-diagonal {half_diag:.1f} needs {hops} ghost-sync hops with blocks {blocks} (setup does 4 rounds)"

def box(x, rho):
    z = ARGS.z0 if ARGS.z0 > 0 else W + s / 2.0 - rho * s   # centre height at the UPRIGHT Archimedes draft d = rho*s (or --z0)
    if ARGS.nopin:
        pins = """      pinPosition          false;         // fully free (6 DOF), as the archived s=30 runs
      pinRotation          false;"""
    else:
        pins = """      pinPosition          true;
      pinAxes              < 1, 1, 0 >;   // hold X,Y; heave (Z) free
      pinRotation          true;
      pinRotationAxes      < 1, 0, 1 >;   // roll about Y free; pitch (X) and yaw (Z) pinned -- Bogner's 'longitudinal axis only'"""
    return f"""   FloatingBox
   {{
      position             < {x:.4f}, {LY/2:g}, {z:g} >;
      halfExtents          < {s/2:g}, {s:g}, {s/2:g} >;
      densityRatio         {rho:g};
      settlingSteps        0;
      initialTiltAngleDeg  {SEED:g};   // seed heel about the free Y axis (Bogner: 2.86)
      comOffsetZ           0;
{pins}
   }}"""

prm = f"""// bogner_exact_3box_s{s}{TAG}_gas.prm   (GENERATED by gen_bogner_exact_prm.py -- edit the script, not this file)
// ─────────────────────────────────────────────────────────────────────────────
// [Claude, 2026-09-19] BOGNER (2017) Sec. 5.3.3 EQUILIBRIUM TEST, GEOMETRY REPRODUCED AS EXACTLY AS THE TEXT ALLOWS, at s = {s} ({k:g}x his resolution).
// Question this answers (THESIS_PLAN.md item 1a): our free-floating cuboid at rho = 0.75 settles ~6 deg above Bogner's 26.565 deg even with his tau and g
// (s=30 box, 96x135x120, 5 deg heel: 32.4 deg), while rho = 0.5 (45.0) and rho = 0.25 (25.3) are fine. Is that offset a property of our setup (geometry,
// heel, free 3D motion) or of the method? Here: his geometry, his tau/g, his heel, his roll-only constraint, ALL THREE densities in one basin.
// (s=16 is his own resolution; s=32 tests his claim that the error shrinks with resolution.)
//{'' if not TAG else chr(10) + '// VARIANT ' + TAG + ': follow-up to the baseline s=16 run (job 4248419), where the 2.86 deg seed decayed to ~0 within 1000 steps and rho 0.5/0.25 stayed upright (target 45 / 26.565): ' + ('seed heel ' + format(SEED, 'g') + ' deg' if abs(SEED-2.86) > 1e-9 else 'seed unchanged') + ('; PINS OFF (fully free 6-DOF boxes)' if ARGS.nopin else '; same pins') + '. Everything else identical to the baseline -- one variable at a time.'}
//
//  * basin {LX} x {LY} x {LZ} cells, lower half liquid (waterline Z = {W:g}); three boxes 2s x s x s = {s}x{2*s}x{s} cells (halfExtents <{s/2:g},{s:g},{s/2:g}>), long axis Y,
//    densities 0.75 / 0.5 / 0.25 left to right, centres on equal thirds of X, Y at mid-width, Z at the upright Archimedes draft (d = rho*s).
//  * tau = 1/1.9 (relaxationRateOverride 1.9), gravityOverride {G:g} ({'his 7.5e-4' if s == 16 else 'his 7.5e-4 / ' + format(k, 'g') + ', keeping H*g = 0.03, his own incompressibility margin'}). The equilibrium angle is static
//    (Archimedes + geometry), so it should not depend on g; g only sets how fast the box gets there.
//  * heel {SEED:g} deg from upright about Y (his value is 2.86; our archived s=30 runs used 5). {'Boxes FULLY FREE (no pins).' if ARGS.nopin else 'Rotation pinned to roll about Y only; X,Y position pinned, heave free.'}
//  * Buoyancy: project default MEM-only (useMemOnlyBuoyancy + useBognerThreePhaseMEM, i.e. his Eq. 18, no separate buoyancy channel). useAnalyticBoxBuoyancy is OFF
//    (explicit) -- no analytic buoyancy anywhere in this run.
//  * Units: avgParticleDiameterOverride {D_AVG:g} / dxSIOverride give absoluteLiquidHeight = d * (1.25 + 1.25) = {W:g}; spheres_out.dat is never opened.
//  * NOT determined by his text, chosen here: side boundaries (X AND Y periodic: the {2*s}-long boxes then leave {LY-2*s} cells to their own image instead of {(LY-2*s)//2} to a wall), his exact timestep count
//    ("until all fluid motions ceased"), and whether his rotation constraint was a pin. So this is not a bit-exact reproduction.
//  * Expected stable heel (his text): 26.565 / 45 / 26.565 deg for rho 0.75 / 0.5 / 0.25. Read the result as tilt folded by the section's 90 deg symmetry
//    (|angle| mod 90, then min(a, 90-a)) from the floaterBox VTK 'rotation' matrix or the [BoxForces] tiltY log lines (one line per box uid, every 500 steps).
//  * Not run yet. Domain {LX}x{LY}x{LZ} / cellsPerBlock <{blocks[0]},{blocks[1]},{blocks[2]}> = {nb[0]}x{nb[1]}x{nb[2]} = {nr} processes; {TIMESTEPS} timesteps, VTK every {VTK_EVERY}.
// ─────────────────────────────────────────────────────────────────────────────

BlockForestParameters
{{
   cellsPerBlock                 < {blocks[0]}, {blocks[1]}, {blocks[2]} >;
   periodicity                   < 1, 1, 0 >;   // X and Y periodic, Z floor/ceiling
   loadSnapshot                  false;
   storeSnapshot                 true;
   snapshotFrequency             10000;
   snapshotBaseFolder            snapshots;
}}

DomainParameters
{{
   domainSize         <{LX}, {LY}, {LZ}>;
   wavePeriods        1;
   liquidHeightFactor 1.25;
   floorHeightFactor  1.25;
   initialAmplitude   0;
}}

PIDParameters
{{
   targetMeanVelocityMagnitude 0;       // still water
   proportionalGain            2e-4;
   derivativeGain              1e-6;
   integralGain                2e-4;
   maxRamp                     1e-4;
   minActuatingVariable        0;
   maxActuatingVariable        1e-3;
}}

PhysicsParameters
{{
   enableWetting          false;
   timesteps              {TIMESTEPS};
   Re                     3100;
   Fr                     1.31;
   We                     15.6188;
   useHardcodedUnits      true;
   relaxationRateOverride {ARGS.omega:g};          // Bogner's tau = 1/1.9 = 0.526; lower omega = more viscous damping ({("Bogner's own value" if abs(ARGS.omega-1.9)<1e-9 else "DAMPING TEST, not Bogner's value")})
   gravityOverride        {G:.6g};
}}

ParticleParameters
{{
   inFileName spheres_out.dat;   // never opened -- see avgParticleDiameterOverride
   bedCopiesInX 1;
   bedCopiesInY 1;
   avgParticleDiameterOverride {D_AVG:g};
   dxSIOverride                0.0004;
   densityRatio 2.55;
   fixingHeightFactor 1.5;
   frictionCoefficient 0.5;
   restitutionCoefficient 0.97;
   numSubCycles 10;
   useLubricationCorrection true;
   useNoSlipParticles false;
   useParticleBed false;
}}

FloatingBoxes
{{
{chr(10).join(box(x, r) for x, r in zip(xs, dens) if not ARGS.only or abs(r - ARGS.only) < 1e-9)}
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
   useAnalyticBoxBuoyancy        false;   // explicitly OFF: no analytic buoyancy in this experiment
}}

EvaluationParameters
{{
   performanceLogFrequency 5000;
   evaluationFrequency     {ARGS.evalfreq};
   dataSliceFrequency      500;   // separate from evaluationFrequency (SimulationParams.h) -- keeps the x-z slice text-file dump throttled even when evaluationFrequency is set to 1 for per-step [BoxForces]/[BoxTorque] logging
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
         mapped_flag;
      }}

      // one mid-width plane (the boxes' cross-section is in X-Z)
      CellBB_filter_0 {{ min < 0, {LY//2}, 0 >;  max < {LX}, {LY//2+1}, {LZ} >; }}

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
out = pathlib.Path(__file__).parent / "prms" / f"bogner_exact_3box_s{s}{TAG}_gas.prm"
out.write_text(prm)
print(f"ghost-sync hops needed {hops} (<=4 ok)")
print(f"s={s}: domain {LX}x{LY}x{LZ}, water {W:g}, g {G:g}, {TIMESTEPS} steps, VTK every {VTK_EVERY}, blocks {blocks} -> {nr} ranks")
print("boxes (x, y, z, rho):", [(round(x, 2), LY / 2, W + s / 2 - r * s, r) for x, r in zip(xs, dens)])
print("wrote", out)
