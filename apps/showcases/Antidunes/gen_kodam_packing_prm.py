#!/usr/bin/env python3
"""Recreate Kodam et al. (2010) Part II, Section 3.3 "Static packing in a cylinder" locally, using
this project's own DropTest.cpp/SetupStages_drop_test.h machinery (fluid-free DEM + floaterCylinder
VTK output, already supports cylinders -- see DropTest.cpp's own header).

Kodam's real test (Table 1 + Section 3.3, checked directly against the paper via Zotero):
  - Cylinder: diameter 8.0mm, length 5.3mm  -> aspect ratio L/D = 0.6625
  - Container: Plexiglas cylinder, internal diameter 50.6mm, height 130mm
  - 250 cylindrical particles, dropped ONE AT A TIME with negligible initial KE and random initial
    orientation, from random positions at the top of the container; the next particle is released
    only once the previous one has come to rest.
  - friction coefficient 0.5, coefficient of restitution 0.85 (Table 1)
  - Reported result (true cylinder DEM vs. experiment): fill height 55.8+-2.2mm vs 53.3+-2.0mm
    (~5% high); bulk solid fraction 0.589+-0.023 vs 0.617+-0.023 (~5% low).

What this recreation keeps EXACT: cylinder aspect ratio (L/D=0.6625), particle count (250), friction
(0.5) and restitution (0.85) from Kodam's own Table 1, and (see below) the SAME real-gravity
calibration and batched-wave release TIMING already proven stable in this project's own
DropTest_cylinders_250_gas.prm.

DELIBERATE approximations, and why:
  - Container cross-section: DropTest's walls are 6 axis-aligned planes (SetupStages_drop_test.h) --
    a true round container wall is not implemented anywhere in this project. Rather than build a new
    curved-wall contact type, the container is a SQUARE footprint with the SAME cross-sectional AREA
    as Kodam's real circular container (bulk solid fraction = volume/[(footprint area)*height] stays a
    fair comparison; only near-wall packing structure right at the boundary differs).
  - Release protocol: NOT strictly one-at-a-time. Two things were tried and found unsafe for a
    one-at-a-time schedule at a LOCAL-run-tractable step count (both found empirically, this session):
      1. Giving every one of 250 particles its own unique, simultaneously-created random (x,y) (all
         250 exist as real bodies from t=0, only "startInactive" distinguishing them) collapsed the
         particle count during the run (250 -> 193 by step 23000, ~1 lost per release event).
      2. Even after fixing that (a small, reused set of well-separated spawn COLUMNS instead), the
         collapse persisted just as badly with a "demo-speed" gravityOverride chosen to make a
         one-at-a-time schedule finish quickly (0.005 LU/step^2, ~250x the value below) -- confirmed
         BY A DIRECT A/B TEST that the established, real-gravity-calibrated DropTest_cylinders_250_gas.prm
         does NOT collapse (numParticles held at 250 across the same step range under normal wave
         release), isolating the cause to gravity magnitude, not particle placement: the fast gravity
         drives impact velocities (~1.3 LU/step over this container's fall height) far outside what
         this project's LinearSpringDashpot contact-time calibration (collisionTimeFactor=4) has ever
         been exercised at elsewhere in this codebase, causing a fast-moving body to tunnel through a
         thin wall within a single DEM sub-step once contact forces spike.
    Given a genuinely safe one-at-a-time schedule at REAL gravity would need ~1.5M timesteps for 250
    particles (impractical for a local interactive run), this recreation instead reuses the EXACT
    proven-stable real-gravity value and the established BATCHED-wave release timing (9 columns fired
    together each wave, matching the working file's own 10-column pattern, just resized for Kodam's
    larger, real cylinder footprint) -- trading strict one-at-a-time for guaranteed stability at a
    tractable local runtime. Small-batch vs. strictly-sequential release is not expected to materially
    change the FINAL static packing outcome this test is actually about.
  - 1 LU = 1 mm (avgParticleDiameterOverride/dxSIOverride copied EXACTLY from the proven file, which
    happen to give this to within 0.02%): every length below reads directly as millimetres, directly
    comparable to the paper's own reported numbers with no conversion.

Run:  python3 gen_kodam_packing_prm.py   ->  prms/DropTest_kodam_packing_gas.prm
"""
import argparse
import math

ap = argparse.ArgumentParser()
ap.add_argument("--radius", type=float, default=4.0, help="cylinder radius (LU=mm); diameter 8mm matches Kodam exactly")
ap.add_argument("--aspect", type=float, default=0.6625, help="length/diameter ratio (Kodam Table 1: 5.3/8.0)")
ap.add_argument("--n", type=int, default=250, help="particle count (Kodam: 250)")
ap.add_argument("--container-diameter-mm", type=float, default=50.6, help="Kodam's real circular container diameter")
ap.add_argument("--container-height-mm", type=float, default=130.0, help="Kodam's real container height (unused for packing, just headroom)")
ap.add_argument("--drop-clearance", type=float, default=30.0, help="extra domain height above the container height for fall clearance (LU)")
ap.add_argument("--wall-margin", type=float, default=2.0, help="min clearance kept between a spawn point and the side walls (LU)")
ap.add_argument("--wave-gap", type=int, default=1500, help="steps between waves -- same as the proven DropTest_cylinders_250_gas.prm")
ap.add_argument("--settle-tail", type=int, default=10000, help="extra steps after the last wave, same convention as the proven file")
ap.add_argument("--friction", type=float, default=0.5, help="Kodam Table 1")
ap.add_argument("--restitution", type=float, default=0.85, help="Kodam Table 1")
ap.add_argument("--seed", type=int, default=42, help="unused (grid layout is deterministic) -- kept for interface symmetry")
ap.add_argument("--out", default="prms/DropTest_kodam_packing_gas.prm")
args = ap.parse_args()

R = args.radius
L = args.aspect * 2.0 * R                      # length from aspect ratio (L/D), D=2R
D = 2.0 * R
container_area_mm2 = math.pi * (args.container_diameter_mm / 2.0) ** 2
side = math.sqrt(container_area_mm2)           # square footprint, area-matched to Kodam's circle
side = round(side)
if side % 3 != 0:
    side += (3 - side % 3)                      # round up to a multiple of 3 for a clean 3x3x1 block decomposition

domain_height = int(round(args.container_height_mm + args.drop_clearance))
spawn_z = domain_height - 10.0                  # 10 LU clearance below the domain's own ceiling plane, matching the
                                                 # proven file's own ~10 LU spawn-to-ceiling gap (170 of 180)

lo = R + L / 2.0 + args.wall_margin
hi = side - lo

# 3x3 = 9 reused spawn columns (resized for Kodam's real, much larger cylinder footprint than the
# original rod-shaped DropTest cylinders) -- same "fixed columns, reused across waves" structure as
# the proven file's own 10 columns, just fewer of them since each column needs far more clearance now.
ncols_grid = 3
step = (hi - lo) / (ncols_grid - 1)
columns = [(lo + gx * step, lo + gy * step) for gy in range(ncols_grid) for gx in range(ncols_grid)]
max_footprint_diag = math.sqrt((2 * R) ** 2 + L ** 2)
min_gap = min(math.hypot(a[0] - b[0], a[1] - b[1]) for i, a in enumerate(columns) for j, b in enumerate(columns) if i != j)
assert min_gap > max_footprint_diag + args.wall_margin, f"columns too close: {min_gap:.2f} vs required >{max_footprint_diag + args.wall_margin:.2f}"

batch = len(columns)
num_waves = math.ceil(args.n / batch)
timesteps = (num_waves - 1) * args.wave_gap + args.settle_tail

print(f"Cylinder: R={R:g} L={L:g} (D={D:g}, L/D={L/D:.4f}, Kodam's own L/D=0.6625)")
print(f"Container: square side={side:g} LU (area {side*side:.1f} vs Kodam's real circle area {container_area_mm2:.1f} mm^2, "
      f"{100*(side*side/container_area_mm2-1):+.1f}% footprint area)")
print(f"{batch} reused spawn columns (min pairwise gap {min_gap:.1f} LU > worst-case footprint diagonal {max_footprint_diag:.1f} LU)")
print(f"Domain: <{side:g},{side:g},{domain_height:g}>, spawn z={spawn_z:g}")
print(f"{args.n} particles in {num_waves} waves of up to {batch}, {args.wave_gap} steps apart -> timesteps={timesteps}")
print(f"Reference: Kodam's true-cylinder fill height 55.8mm = {55.8/D:.3f} diameters -> {55.8/D*D:.1f} LU here; "
      f"solid fraction 0.589 (dimensionless, compare directly)")

lines = []
lines.append(f"// DropTest_kodam_packing_gas.prm  (auto-generated by gen_kodam_packing_prm.py)")
lines.append("// ─────────────────────────────────────────────────────────────────────────────")
lines.append("// Local recreation of Kodam et al. (2010) Part II Sec. 3.3 \"Static packing in a cylinder\":")
lines.append(f"// {args.n} cylinders (R={R:g}, L={L:g} LU=mm, L/D={L/D:.4f}, matching Kodam's real 8.0/5.3mm exactly),")
lines.append(f"// dropped in waves of up to {batch} (one per reused spawn column) into a square, area-matched stand-in")
lines.append("// for Kodam's real 50.6mm-diameter circular container. friction/restitution taken directly from Kodam's")
lines.append("// Table 1 (0.5/0.85). Gravity and wave timing REUSE THE EXACT PROVEN-STABLE VALUES from the existing")
lines.append("// DropTest_cylinders_250_gas.prm -- see this script's own module docstring for why a faster \"demo")
lines.append("// speed\" gravity was tried and found to cause the particle count to collapse (contact-model")
lines.append("// instability at the resulting impact velocity), and why strict one-at-a-time release was dropped")
lines.append("// in favour of the already-validated batched-wave schedule. 1 LU = 1 mm.")
lines.append("// ─────────────────────────────────────────────────────────────────────────────")
lines.append("")
lines.append("BlockForestParameters")
lines.append("{")
lines.append(f"   cellsPerBlock                 < {side//3}, {side//3}, {domain_height} >;  // 3x3x1 = 9 processes")
lines.append("   periodicity                   < 0, 0, 0 >;")
lines.append("   loadSnapshot                  false;")
lines.append("   storeSnapshot                 false;")
lines.append("   snapshotFrequency             10000;")
lines.append("   snapshotBaseFolder            snapshots;")
lines.append("}")
lines.append("")
lines.append("DomainParameters")
lines.append("{")
lines.append(f"   domainSize         <{side}, {side}, {domain_height}>;")
lines.append("   wavePeriods        1;      // unused -- required by parser")
lines.append("   liquidHeightFactor 2.5;    // unused -- required by parser")
lines.append("   floorHeightFactor  4.1393; // unused -- required by parser")
lines.append("   initialAmplitude   0;      // unused -- required by parser")
lines.append("}")
lines.append("")
lines.append("PIDParameters")
lines.append("{")
lines.append("   targetMeanVelocityMagnitude 0;")
lines.append("   proportionalGain            0;")
lines.append("   derivativeGain              0;")
lines.append("   integralGain                0;")
lines.append("   maxRamp                     0;")
lines.append("   minActuatingVariable        0;")
lines.append("   maxActuatingVariable        0;")
lines.append("}")
lines.append("")
lines.append("PhysicsParameters")
lines.append("{")
lines.append("   enableWetting     false;")
lines.append(f"   timesteps         {timesteps};")
lines.append("   Re                3100;")
lines.append("   Fr                1.31;")
lines.append("   We                15.6188;")
lines.append("   useHardcodedUnits true;")
lines.append("   gravityOverride   2.1371500584054327e-05;  // EXACT value from the proven DropTest_cylinders_250_gas.prm")
lines.append("                                              // (\"true g_SI=9.81\"); confirmed empirically stable this session.")
lines.append("}")
lines.append("")
lines.append("ParticleParameters")
lines.append("{")
lines.append("   inFileName                  spheres_out.dat; // never read -- overrides below bypass it")
lines.append("   bedCopiesInX                1;")
lines.append("   bedCopiesInY                1;")
lines.append("   densityRatio                1.0;")
lines.append("   fixingHeightFactor          1.5;")
lines.append(f"   frictionCoefficient         {args.friction:g};  // Kodam Table 1")
lines.append(f"   restitutionCoefficient      {args.restitution:g}; // Kodam Table 1")
lines.append("   numSubCycles                20;")
lines.append("   useLubricationCorrection    false;")
lines.append("   useNoSlipParticles          false;")
lines.append("   useParticleBed               false;")
lines.append("   avgParticleDiameterOverride 8.001943372119138;   // EXACT value from the proven file (gives 1 LU ~ 1mm to 0.02%)")
lines.append("   dxSIOverride                0.000362694127490625; // EXACT value from the proven file")
lines.append("}")
lines.append("")
lines.append("// No Pillars / CylinderPillars / ArchTunnels -- just the 6 auto-created domain walls")
lines.append("// (SetupStages_drop_test.h) standing in for Kodam's cylindrical container (see module docstring).")
lines.append("FloatingCylinders")
lines.append("{")
count = 0
wave = 0
while count < args.n:
    this_batch = min(batch, args.n - count)
    settling = wave * args.wave_gap
    for ci in range(this_batch):
        x, y = columns[ci]
        # deterministic pseudo-random-looking orientation per (column, wave), no RNG needed for reproducibility:
        # a low-discrepancy-ish spread over the unit sphere using the golden angle, distinct per particle index.
        idx = count + ci
        theta = math.acos(1 - 2 * ((idx * 0.6180339887498949) % 1.0))
        phi = 2 * math.pi * ((idx * 0.7548776662466927) % 1.0)
        ax, ay, az = math.sin(theta) * math.cos(phi), math.sin(theta) * math.sin(phi), math.cos(theta)
        inactive = "" if wave == 0 else " startInactive true;"
        lines.append(f"   FloatingCylinder {{ position < {x:.3f}, {y:.3f}, {spawn_z:g} >; radius {R:g}; length {L:g}; "
                     f"axis < {ax:.4f}, {ay:.4f}, {az:.4f} >; densityRatio 1.0; settlingSteps {settling};{inactive} }}")
    count += this_batch
    wave += 1
lines.append("}")
lines.append("")
lines.append("ModelParameters")
lines.append("{")
lines.append("   useAnalyticalContactDetection true;  // Kodam/Guo cylinder-cylinder closed form (PaperAnalyticalContacts.h)")
lines.append("   useMultiPointBoxContact       false;")
lines.append("   useBoxOverlapCap              false;")
lines.append("   useSpheroBoxContact           false;")
lines.append("   spheroBoxRadius               0.05;")
lines.append("   useSpheroBoxAnalytical        false;")
lines.append("   useTightLayer1                true;")
lines.append("   useContactFilterOwnerTiebreak true;")
lines.append("   pdfReconstructionModel        OnlyMissing;")
lines.append("   pdfRefillingModel             EquilibriumRefilling;")
lines.append("   excessMassDistributionModel   EvenlyNewInterfaceFallbackLiquid;")
lines.append("   curvatureModel                FiniteDifferenceMethod;")
lines.append("   useSimpleMassExchange         false;")
lines.append("   cellConversionThreshold       1e-2;")
lines.append("   cellConversionForceThreshold  1e-1;")
lines.append("   useBognerReconstruction       true;")
lines.append("   usePass3Reresolution          true;")
lines.append("   useBognerThreePhaseMEM        true;")
lines.append("   useMemOnlyBuoyancy            true;  // REQUIRED true -- see DropTest.cpp's own header comment")
lines.append("   logFloaterHistory              false;")
lines.append("   logFloaterHistoryStartStep     0;")
lines.append("   logFloaterHistoryEndStep       10000;")
lines.append("   logContactAnomalies            false;")
lines.append("   contactAnomalyStartStep        0;")
lines.append("   contactAnomalyEndStep          10000;")
lines.append("   useHydrodynamicForceTimeAveraging true;")
lines.append("}")
lines.append("")
lines.append("EvaluationParameters")
lines.append("{")
lines.append("   performanceLogFrequency 5000;")
lines.append("   evaluationFrequency     1000;")
lines.append("   baseFolderName          eval;")
lines.append("   logContactDetectionStatsPeriodically true;")
lines.append("}")
lines.append("")
lines.append("VTK")
lines.append("{")
lines.append("   fluid_field")
lines.append("   {")
lines.append("      writeFrequency       500;")
lines.append("      ghostLayers          0;")
lines.append("      baseFolder           vtk-out;")
lines.append("      samplingResolution   1;")
lines.append("   }")
lines.append("")
lines.append("   domain_decomposition")
lines.append("   {")
lines.append("      writeFrequency             10000000000;")
lines.append("      baseFolder                 vtk-out;")
lines.append("      outputDomainDecomposition  true;")
lines.append("   }")
lines.append("}")

with open(args.out, "w") as f:
    f.write("\n".join(lines) + "\n")
print(f"Wrote {args.out}")
