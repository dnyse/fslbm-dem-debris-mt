#!/usr/bin/env python3
"""Fixed-heel RIGHTING-MOMENT scan (Bogner 2017 Sec. 5.3.2 protocol) with the raw MEM torque, for the cube-section box.

    python3 gen_bogner_righting_scan_prm.py --s 16      ->  prms/bogner_righting_scan_s16_a<tag>_gas.prm  (one per heel angle)
                                                            + prms/bogner_righting_scan_s16_manifest.json  (ideal reference values)

Question (THESIS_PLAN.md item 1a): why does rho = 0.75 settle ~6-8 deg above Bogner's 26.565 deg while rho = 0.25 and 0.5 agree? Ideal theory is
SYMMETRIC under rho -> 1-rho: the (1-rho) box is the 180-degree point-reflected image of the rho box (a square maps onto itself), so at the same heel angle (heave set
so that Archimedes holds) the two hydrostatic torques are EQUAL: M_{1-rho}(alpha) = +M_rho(alpha) -- the two stability curves are identical, hence the identical stable
heel (26.565 deg for both 0.25 and 0.75). (An earlier draft of this comment said 'negative'; the generator's own self-check on the exact numbers showed equal.) A FIXED box
at each heel angle isolates the hydrostatic force computation (no dynamics): if the MEM torque of rho=0.75 differs from rho=0.25's, the force computation is biased in one of them.

Per run (one heel angle): the three-box Bogner basin (130x40x40 for s=16, tau=1/1.9, g=7.5e-4) with boxes rho 0.75 / 0.5 / 0.25, ALL held (position and
rotation pinned in every axis) at the same heel angle, each at the heave where its submerged AREA equals rho*s^2 (Archimedes), MEM-only buoyancy, NO analytic
buoyancy. The raw MEM hydrodynamic torque/force is read from the floater-history state CSV over a short window after the fluid has settled.
"""
import argparse, json, math, pathlib

ap = argparse.ArgumentParser()
ap.add_argument("--s", type=int, default=16)
ap.add_argument("--timesteps", type=int, default=0, help="default 20000*s/16")
ap.add_argument("--window", type=int, default=300, help="logging window length in steps (ends at timesteps-700)")
ap.add_argument("--boxes", default="0.75,0.5,0.25", help="densities of the held boxes, left to right (default: all three of Bogner's basin)")
ap.add_argument("--angles", default="0,5,10,15,20,25,26.565,30", help="comma-separated heel angles in degrees")
ap.add_argument("--prefix", default="bogner_righting_scan", help="file-name prefix for the prms and the manifest")
ap.add_argument("--nofill", action="store_true", help="[2026-09-21] write ModelParameters.useFillWeighting false (raw two-phase MEM force, ablation of the waterline phi-weighting); default true")
ap.add_argument("--ibb", action="store_true", help="[2026-09-21] write ModelParameters.useInterpolatedBounceBack true (CLI interpolated bounce-back on the box links, SimpleBB.h); default: simple bounce-back")
ARGS = ap.parse_args()
s = ARGS.s; k = s / 16.0
LX, LY, LZ = int(round(130 * k)), int(round(40 * k)), int(round(40 * k))
W = LZ / 2.0
ZW = W + 0.5                       # actual initial free surface (fill level is initialised at absoluteLiquidHeight + 0.5; measured 56.52 vs 56.02 elsewhere)
G = 7.5e-4 / k
TIMESTEPS = ARGS.timesteps or int(20000 * k)
WIN_START, WIN_END = TIMESTEPS - 1000, TIMESTEPS - 1000 + ARGS.window
D_AVG = 8.0 * k
DENS = [float(v) for v in ARGS.boxes.split(',')]
xs = [LX * (i + 0.5) / len(DENS) for i in range(len(DENS))]
ANGLES = [float(v) for v in ARGS.angles.split(',')]
BLOCKS = {16: (10, 8, 40), 32: (20, 16, 40)}[s]
nb = (LX // BLOCKS[0], LY // BLOCKS[1], LZ // BLOCKS[2]); nr = nb[0] * nb[1] * nb[2]
L_BOX = 2.0 * s                    # box length along Y

def tag(a):                        # a05, a10, ..., a26p565
    return "a" + (f"{int(round(a)):02d}" if abs(a - round(a)) < 1e-9 else format(a, "g").replace(".", "p"))

def square(zc, alpha_deg):        # vertices of the side-s square centred at (0,zc) rotated by alpha, in the (x,z) plane
    a = math.radians(alpha_deg); c, sn = math.cos(a), math.sin(a); h = s / 2.0
    return [(x * c + z * sn, -x * sn + z * c + zc) for x, z in ((-h, -h), (h, -h), (h, h), (-h, h))]

def clip_below(poly, zw):         # Sutherland-Hodgman: part of a convex polygon with z <= zw
    out = []
    for i in range(len(poly)):
        p, q = poly[i], poly[(i + 1) % len(poly)]
        pin, qin = p[1] <= zw, q[1] <= zw
        if pin: out.append(p)
        if pin != qin:
            t = (zw - p[1]) / (q[1] - p[1]); out.append((p[0] + t * (q[0] - p[0]), zw))
    return out

def area_centroid(poly):
    if len(poly) < 3: return 0.0, (0.0, 0.0)
    A = cx = cz = 0.0
    for i in range(len(poly)):
        x0, z0 = poly[i]; x1, z1 = poly[(i + 1) % len(poly)]; cr = x0 * z1 - x1 * z0
        A += cr; cx += (x0 + x1) * cr; cz += (z0 + z1) * cr
    A *= 0.5
    return A, (cx / (6 * A), cz / (6 * A))

def heave_for(rho, alpha):        # centre height so that the submerged AREA = rho*s^2 (bisection; area DEcreases as zc rises)
    lo, hi = ZW - s * 1.5, ZW + s * 1.5
    for _ in range(80):
        mid = 0.5 * (lo + hi); A, _ = area_centroid(clip_below(square(mid, alpha), ZW))
        if A > rho * s * s: lo = mid      # too deep -> raise the box
        else: hi = mid
    return 0.5 * (lo + hi)

manifest = []
for alpha in ANGLES:
    for rho in DENS:
        zc = heave_for(rho, alpha); A, (bx, bz) = area_centroid(clip_below(square(zc, alpha), ZW))
        F = G * A * L_BOX                        # ideal buoyancy force (rho_fluid = 1 LU) = weight at equilibrium heave
        manifest.append(dict(alpha=alpha, rho=rho, z_centre=zc, area_sub=A, area_expected=rho * s * s, x_B_minus_x_G=bx, F_buoy_ideal=F,
                             M_ideal_mag=abs(F * bx), M_ideal_signed_xz=F * bx, box_index=DENS.index(rho)))
# self-check of the symmetry claim on the IDEAL numbers: M_0.25(alpha) == +M_0.75(alpha)
for alpha in ANGLES:
    m = {e["rho"]: e["M_ideal_signed_xz"] for e in manifest if e["alpha"] == alpha}
    if 0.25 in m and 0.75 in m:
        assert abs(m[0.25] - m[0.75]) < 1e-9 * max(1.0, abs(m[0.75])), (alpha, m)
assert all(abs(e["area_sub"] - e["area_expected"]) < 1e-6 for e in manifest)

def box(x, rho, alpha):
    zc = next(e["z_centre"] for e in manifest if e["alpha"] == alpha and e["rho"] == rho)
    return f"""   FloatingBox
   {{
      position             < {x:.4f}, {LY/2:g}, {zc:.4f} >;
      halfExtents          < {s/2:g}, {s:g}, {s/2:g} >;
      densityRatio         {rho:g};
      settlingSteps        0;
      initialTiltAngleDeg  {alpha:g};
      comOffsetZ           0;
      pinPosition          true;
      pinAxes              < 1, 1, 1 >;   // held at the Archimedes heave for this rho and heel
      pinRotation          true;
      pinRotationAxes      < 1, 1, 1 >;   // held at the heel angle
   }}"""

outdir = pathlib.Path(__file__).parent / "prms"
for alpha in ANGLES:
    txt = f"""// {ARGS.prefix}_s{s}_{tag(alpha)}_gas.prm   (GENERATED by gen_bogner_righting_scan_prm.py -- edit the script, not this file)
// ─────────────────────────────────────────────────────────────────────────────
// [Claude, 2026-09-19] FIXED-HEEL RIGHTING-MOMENT SCAN, heel {alpha:g} deg, s = {s} (Bogner 2017 Sec. 5.3.2 protocol, cube section 2s x s x s = {s}x{2*s}x{s}).
// Purpose: measure the raw MEM torque on a HELD box and compare rho = 0.75 with rho = 0.25. Ideal theory: M_{{1-rho}}(alpha) = +M_rho(alpha), i.e. the rho=0.25 and rho=0.75 torque curves are IDENTICAL (checked on the ideal numbers in the
// generator). See gen_bogner_righting_scan_prm.py's docstring. NOT RUN YET.
//  * basin {LX}x{LY}x{LZ}, lower half liquid, initial free surface at Z = {ZW:g} (fill level is initialised half a cell above the nominal {W:g}); X and Y periodic; tau = 1/1.9, g = {G:g}.
//  * {len(DENS)} boxes rho {' / '.join(map(str, DENS))} at X = {', '.join(f'{x:.1f}' for x in xs)}, ALL at heel {alpha:g} deg about Y, each held (position AND rotation pinned in every axis) at the heave where its submerged
//    cross-section area equals rho*s^2: Z = {', '.join(format(next(e['z_centre'] for e in manifest if e['alpha']==alpha and e['rho']==r), '.3f') for r in DENS)} for rho {', '.join(map(str, DENS))} (ideal geometry, waterline {ZW:g}).
//  * MEM-only buoyancy (useMemOnlyBuoyancy + useBognerThreePhaseMEM); useAnalyticBoxBuoyancy explicitly OFF (no analytic buoyancy).
//  * {TIMESTEPS} steps, then the floater-history state CSV (raw MEM hydroForce/hydroTorque of the local copy of each box) is logged for steps {WIN_START}..{WIN_END} only -- the CSV has blown the home
//    quota before, so the window is short and the fluid has {WIN_START} steps to settle (acoustic damping is slow at tau=0.526).
//  * Domain {LX}x{LY}x{LZ} / cellsPerBlock <{BLOCKS[0]},{BLOCKS[1]},{BLOCKS[2]}> = {nb[0]}x{nb[1]}x{nb[2]} = {nr} processes.
// ─────────────────────────────────────────────────────────────────────────────

BlockForestParameters
{{
   cellsPerBlock                 < {BLOCKS[0]}, {BLOCKS[1]}, {BLOCKS[2]} >;
   periodicity                   < 1, 1, 0 >;
   loadSnapshot                  false;
   storeSnapshot                 false;
   snapshotFrequency             100000000;
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
   targetMeanVelocityMagnitude 0;
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
   relaxationRateOverride 1.9;          // Bogner's tau = 1/1.9
   gravityOverride        {G:.6g};
}}

ParticleParameters
{{
   inFileName spheres_out.dat;   // never opened
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
{chr(10).join(box(x, r, alpha) for x, r in zip(xs, DENS))}
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
   useFillWeighting              {'false' if ARGS.nofill else 'true'};   // [2026-09-21] false = ablation (no phi weighting of the MEM force)
   useInterpolatedBounceBack     {'true' if ARGS.ibb else 'false'};   // [2026-09-21] IBB (CLI) on the moving-obstacle links; false = simple bounce-back
   useBognerThreePhaseMEM        true;
   useMemOnlyBuoyancy            true;
   useAnalyticBoxBuoyancy        false;   // explicitly OFF

   logFloaterHistory             true;    // per-rank floater-history/rank_<N>_state.csv (raw MEM hydroForce/hydroTorque), WINDOWED below
   floaterHistoryBaseFolder      floater-history;
   logFloaterHistoryStartStep    {WIN_START};
   logFloaterHistoryEndStep      {WIN_END};
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
      writeFrequency       {TIMESTEPS // 2};
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
    (outdir / f"{ARGS.prefix}_s{s}_{tag(alpha)}_gas.prm").write_text(txt)

(outdir / f"{ARGS.prefix}_s{s}_manifest.json").write_text(json.dumps(dict(s=s, LX=LX, LY=LY, LZ=LZ, z_water=ZW, g=G, box_length=L_BOX, window=[WIN_START, WIN_END],
                                                                                     timesteps=TIMESTEPS, densities=DENS, xs=xs, angles=ANGLES, tags=[tag(a) for a in ANGLES], cases=manifest), indent=1))
print(f"s={s}: {len(ANGLES)} heel angles {ANGLES}, tags {[tag(a) for a in ANGLES]}; {TIMESTEPS} steps, log window {WIN_START}..{WIN_END}; {nr} ranks")
print(f"ideal |M| (LU) by heel [rho {' / '.join(map(str, DENS))}] and heave z_c:")
for alpha in ANGLES:
    row = {e['rho']: e for e in manifest if e['alpha'] == alpha}
    print(f"  {alpha:6.3f} deg: |M| " + " / ".join(f"{row[r]['M_ideal_mag']:7.3f}" for r in DENS) + "   z_c " + " / ".join(f"{row[r]['z_centre']:6.3f}" for r in DENS))
