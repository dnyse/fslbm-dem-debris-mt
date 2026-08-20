# Antidunes — Change Log

---

## Follow-up on DEMValidation: conservation checks + a real Kodam Part II Test 1 replica (2026-08-12)

Two additions to the dry-run below, prompted by the question "does this actually validate against the
papers' data, and does conservation hold?":

**Conservation checks added to every scenario** (`logConservation`/`logMomentumConservation`/
`totalMechanicalEnergy` in `DEMValidationScenarios.h`): kinetic energy must never INCREASE across any
dissipative collision (a model-agnostic hard bound, independent of restitution/friction calibration);
linear momentum must be conserved to numerical precision for the isolated two-body pair tests; total
mechanical energy (KE+PE) must not trend upward between sampling checkpoints in the multi-body gravity
tests. All pass cleanly across every scenario (momentum to ~1e-6, energy-non-increase with wide margin).

**Energy budget (not just non-increase) added** (`logEnergyBudget`, `kineticEnergyOf`): the one-sided
"never increases" bound is necessarily weak for the deliberately-inelastic scenarios (restitution <1 by
design, matching real materials) -- it would pass even if the simulation dissipated a wildly wrong
AMOUNT of energy. Added a companion check comparing the simulated post-impact KE against the EXACT
value the analytical model (or, for pairRestitution at its default frictionless setting, a closed-form
reduced-mass formula) predicts, not just its sign. Wall-impact tests: 1-14% relative error (same order
as the existing velocity-comparison tolerance, since both trace back to the same friction-model gap).

**New `dem_validation_wall_cylinder_kodam_test1.prm`**: a REPLICA of Kodam et al. (2010) Part II's own
Test 1 -- their exact geometry (diameter 8.0mm/length 5.3mm, aspect ratio 0.6625), their exact
restitution (0.85, calibrated by them from Plexiglas sphere drops), their exact protocol (frictionless,
gravity-free, purely-normal impact velocity, zero initial spin, 1 m/s). This is a genuinely stronger
validation than the earlier cylinder test: at tilt=0deg/90deg the simulated rebound ratio is checked
DIRECTLY against the paper's own stated number (0.85), not just against our own re-derived formula --
PASSES at <0.6% error at both points. Intermediate tilt angles show larger *relative* errors at the
tolerance used (0.02, tight since friction is off), but the *absolute* deviations are small (~0.02-0.03
on a ~1 m/s scale) and concentrated near natural zero-crossings of the rebound velocity (a relative-error
metric is inherently unstable near zero regardless of physical correctness) -- not evidence of a real
discrepancy, but flagged honestly rather than tuned away.

---

## NEW: fluid-free DEM validation dry-run, `DEMValidation` (2026-08-12)

**New standalone executable** (`DEMValidation.cpp`, `AnalyticalSolutions.h`, `DEMValidationScenarios.h`,
7 `prms/dem_validation_*.prm` configs, new CMake target) that runs the PRODUCTION DEM code
(`SolidCylinderContactDetection`, `ComputeContacts`, `VelocityVerletPre/PostForceUpdate` from
`ParticleKernels.h`, unmodified) through a minimal fluid-free sub-cycle loop and checks the resulting
DYNAMICS against closed-form analytical solutions and literature values. Existing tests
(`ContactDetectionTest.cpp`, `PaperAnalyticalContactTest.cpp`, `validation/*.py`) only ever checked
static contact geometry (given two placed bodies, is the reported penetration/normal/point correct?) —
nothing previously ran a real time-stepped collision and checked the result.

**Physics**: one general rigid-body single-point-impulse solver (`AnalyticalSolutions.h::rigidBodyImpact`,
the standard Coulomb-friction two-branch impulse-momentum method) applied to every shape rather than
per-shape formulas — specialises to the classic oblique sphere-wall impact, to Kodam et al. (2010) Part
II's cylinder-wall test (contact point taken directly from the not-yet-simulation-wired
`PaperAnalyticalContacts.h` closed forms), and — via an internal self-check
(`selfCheckHousnerRelativeError`) — reproduces Housner (1963)'s rocking-block restitution exactly in the
zero-tangential-velocity limit.

**Results**: 6/7 scenarios PASS —
- `wallImpactSphere`/`wallImpactCylinder`/`wallImpactBox`: single-body wall-bounce vs. the analytical
  solver across a tilt/angle sweep (tolerance 0.25-0.3, calibrated against a real, documented, expected
  gap between mesa_pd's Cundall-Strack tangential-spring friction model and the idealised instantaneous-
  Coulomb-impulse assumption — see `AnalyticalSolutions.h`'s calibration note).
- `pairRestitution`/`energyConservation`: central head-on collisions for all 6 shape-pair combinations
  (incl. Box×Box and Cylinder×Cylinder, the two pairs with no analytic contact overload — GJK/EPA only).
  Restitution matches to <1%; elastic (e=1) energy conservation matches to ~3e-5 — Feng (2023)'s "ECC"
  criterion, satisfied cleanly across every pair. Box-box/cylinder-cylinder additionally get a LIVE
  per-step cross-check of `PaperAnalyticalContacts.h`'s closed forms against the production GJK/EPA
  geometry actually driving the dynamics (`paperAnalyticalMaxOverlapDelta`, measured 0 in this run).
- `packing`: mixed sphere+box+cylinder bed, gravity-settled into a bounded pit; bulk solid fraction 0.415,
  between the sphere (~0.64) and cube (~0.78) RCP literature bounds for the mix.
- `angleOfRepose`: **FAILS** (~1deg, i.e. flat) — root-caused, not just tolerance-missed: mesa_pd's
  `LinearSpringDashpot` implements pure Coulomb sliding friction with **no rolling-resistance torque**, so
  particles roll/slide out to a flat layer instead of mounding up, regardless of shape (sphere vs.
  elongated rod) or containment (open floor vs. bounded pit) — all tried, all failed the same way. A real,
  literature-documented DEM contact-model gap this validation surfaced, not a test bug; see
  `dem_validation_angle_of_repose.prm`'s header and `DEMValidationScenarios.h`'s calibration notes for
  the full diagnosis trail.

**Bug found in shared mesa_pd code** (worked around locally, not patched at the source):
`mesa_pd::data::LinkedCells`'s constructor leaves `infiniteParticles_` at its default-member-initializer
value (0, via `std::atomic<int>`'s C++20 value-init) instead of the `-1` sentinel its own traversal logic
expects — only `.clear()` sets it correctly. `BedGeneration.cpp` already calls `.clear()` every step for
an unrelated reason (reusing one long-lived `LinkedCells` object) and so never hits this; a *freshly
constructed* `LinkedCells` (this app's pattern) needs the same call for a different reason. Manifests as
a genuine self-pairing crash (`AnalyticalContactDetection`'s "colliding with itself" assertion) when
index 0 happens to be the sole finite body in the scene — exactly the single-body wall-impact tests here.
Worth a proper fix in `LinkedCells`'s constructor at some point.

---

## Mixed single-pillar test scaled up: floaters ×2, pillar ×1.5 (2026-07-16)

**`prms/mixed_floaters_single_pillar_gas.prm` resized** — sphere radius 5→10, box `halfExtents`
`<3,3,3>`→`<6,6,6>`, log radius/length 3/10→6/20, pillar radius 8→12 (`CylinderPillars` and the
matching `BodyBC` obstacle radius). Doubling every floater's linear size roughly doubles its reach
too, so the old column/row layout would have collided at t=0 — the whole 3-column × 2-row layout
was **re-solved from scratch** (not scaled by the same 2×) and every pairwise + floater–pillar
clearance re-verified numerically (box reach conservatively bounded by its half-diagonal
`6√3 ≈ 10.39` regardless of the small ±5° tilt):

- New columns: sphere x=0, box x=24, log x=44 (was 20/30/40); rows at y=13/47 (was 15/45).
- Tightest margin in the new layout is the log↔pillar clearance (≈3.19 LU) — the log's own 20 LU
  length brings its nearer end closer to the pillar's centreline than its column offset alone
  would suggest (nearest endpoint sits only `length/2 = 10` LU closer to y=30 than the log's own
  centre). Wall clearance is also tight (≈2.61 LU from the N/S no-slip walls) — do not shrink the
  ±17-ish row half-spacing further without rechecking both numbers.
- Header comment in the prm rewritten with the full clearance derivation so future resizes don't
  have to re-derive it from scratch.

## Floater initial velocity + still-water sphere–sphere collision test (2026-07-16)

**New `FloatingSphereParams::initialVelocity`** (`SimulationParams.h`, default `<0,0,0>`,
backward-compatible) — a linear velocity imparted at creation via `mesa_pd`'s
`setLinearVelocity`, wired into the sphere-creation loop in `SetupStages.h`. General mechanism,
not test-specific: only `FloatingSphereParams` has it so far — `FloatingBoxParams` /
`FloatingCylinderParams` would need the same field + `setLinearVelocity` call added to their own
creation loops (`SetupStages.h` ~line 444 / ~line 502) if a box–box or log–log collision test is
wanted later.

**New `prms/sphere_sphere_collision_gas.prm` + `jobs/run_72_sphere_sphere_collision.sbatch`** —
still water (no pillars, `targetMeanVelocityMagnitude 0`), two spheres (r=6) 40 LU apart on the
domain centreline given opposing `initialVelocity <±0.05,0,0>` so they close and collide head-on
near x=64 (domain centre) instead of relying on flow-driven drift. Gap between surfaces at t=0 is
28 LU (clear); first contact expected ≈280 steps in, ignoring hydrodynamic drag deceleration.

## Buoyancy-plane VTK showcase: `box_stillwater_showcase_gas.prm` (2026-07-16)

**New prm+job pair** exercising `ModelParameters.writeFloaterPlaneVTK` for the first time in any
prm (confirmed via grep beforehand — the flag existed and was wired in `BuildTimeloop.h` but no
prm had ever turned it on) together with `useAnalyticBoxBuoyancy true`, so the exported
`planeOrigin`/`planeNormal` VTK fields are *exactly* the plane the analytic clip (`clipBoxAnalytic`)
cuts against, with no sampling-vs-analytic cross-check mismatch to explain away in the
visualization.

- Same still-water single-cube scenario as `cube_stillwater_gas.prm` (edge 26, released
  `initialTiltAngleDeg 45` — its stable, edge-down orientation), but **X/Y translation pinned**
  (`pinPosition true; pinAxes <1,1,0>`) after observing the un-pinned box visibly wander side to
  side while settling. This is safe here — unlike the earlier-documented "pinning x,y BACKFIRES"
  failure mode (further down this log) — because that failure was specific to a box needing a
  large *righting rotation* from an unstable start; this box already starts at its stable tilt and
  only needs to settle vertically, so pinning XY removes the distracting lateral drift without
  killing any roll-damping the body actually needs.
- `jobs/run_72_box_stillwater_showcase.sbatch` follows the standard template.

## Floater×pillar interaction test suite (2026-07-16)

**7 new isolated single-interaction prm/sbatch pairs**, filling gaps in interaction coverage that
the existing two-pillar jam/clog tests didn't isolate (one case per shape × {glancing contact,
clean gap-pass negative control}, plus rotation/oblique contact where physically meaningful — a
sphere has no rotation case since it's rotationally symmetric):

- `sphere_pillar_glance_gas.prm` / `sphere_pillar_gap_pass_gas.prm`
- `box_pillar_glance_gas.prm` / `box_pillar_gap_pass_gas.prm` / `box_pillar_corner_gas.prm`
  (corner-first contact: redesigned head-on, Δy=0, after realizing `initialTiltAngleDeg` rotates
  about world Y only and doesn't change a box's Y-extent, so a lateral/offset approach would have
  looked identical tilted or not — the tilt-induced increase in X-reach only shows up head-on)
- `cylinder_pillar_gap_pass_gas.prm` / `cylinder_pillar_oblique_gas.prm` (45°-yawed log axis,
  still horizontal so it floats normally — the general skew-axis case for
  `SolidCylinderContactDetection`'s fixed-point iteration, deliberately not a vertical-buoy log,
  which would be a different, unvalidated flotation mode)

**New `prms/mixed_floaters_single_pillar_gas.prm`** — 2 of each shape (6 floaters total) around a
single pillar, all start positions analytically checked pairwise + against the pillar for zero
overlap at t=0 (see the 2026-07-16 resize entry above for the current, post-scale-up numbers).

**8 matching sbatch scripts** added under `jobs/`, all following the established template
(`run_72_sphere_pillar_glance.sbatch`, …, `run_72_mixed_floaters_single_pillar.sbatch`).

## Static-H ablation: `staticReferenceHeightOverride` (2026-07-10)

**New prm knob `ModelParameters.staticReferenceHeightOverride`** (default −1 = disabled). If ≥ 0 it
replaces the flat `absoluteLiquidHeight` (~56.02) as the STATIC MEM reference height H everywhere the
*reconstruction reference* is meant — a pure diagnostic ablation, motivated by the phantom/dynamic-H
A/B results: dynamic-H showed that *raising* H drives the floater down, so this pins H at an extreme
fixed value (domain top, z=180) to isolate the offset effect without the dynamic-H feedback loop.

- `SimulationParams.h`: new `ModelParams::staticReferenceHeightOverride` + convenience
  `SimulationParams::staticReferenceHeight()` (override if ≥ 0, else `computed.absoluteLiquidHeight`).
- Wired into all three consumers of the static H (fluid initialisation is deliberately untouched):
  - `SetupStages.h` — `HydrostaticDensityFtor` for SimpleBB's rhoHydStat (the force-side MEM
    reference; the lever that matters in the baseline configuration) + the `memReferenceHeight` seed.
  - `BuildTimeloop.h` — the Bogner reconstruction manager's hydrostatic-reconstruction H (only active
    with `useHydrostaticReconstruction true`).
  - `ParticleKernels.h` — `surfaceHeight_`, so the `[PhantomDecomp]` deltaH/F_pred diagnostic stays
    consistent with what the reconstruction actually subtracts.

**Physics probed**: shifting H by ΔH adds a constant Δρ = 3·|a_z|·ΔH to the subtracted reference. On a
closed (fully submerged) body that offset integrates to ~zero net force; on the open wetted cap of a
floater straddling the surface it is unbalanced → expected extra DOWNWARD force ≈ |a_z|·ΔH·A_proj. At
ΔH = 124 and A_proj up to π·10² this is O(0.27) — ~13× the test sphere's weight. Prediction: baseline
behaviour while airborne/deep, violent submersion the moment the sphere straddles the waterline —
the mirror image of the levitation baseline.

**Test job**: `prms/sphere_drop_static_h_top_gas.prm` + `jobs/run_72_static_h_top.sbatch` — byte-identical
scenario to `sphere_drop_two_pillars_gas.prm` (the phantom/dynamic-H A/B control, RUN A) except
`staticReferenceHeightOverride 180`. Compare `[FloaterForces]`/`[FloaterBuoyancy]` against RUN A of
`phantom.log` / `dynamic.log`.

---

## Floating CYLINDERS (true logs) + full shape×shape / shape×obstacle contact matrix (2026-07-08)

**New floater shape: `FloatingCylinder` (type 4).** Until now "logs" were Y-elongated boxes; mesa_pd had
no finite-cylinder solid (only the infinite hollow `CylindricalBoundary` used for pillars). A real
cylinder shape was added to mesa_pd core and wired through the whole coupling stack, and the contact
system was rebased so that EVERY floater shape (sphere, box, cylinder) interacts with every other floater
and with every obstacle (halfspace walls/floor, box pillars, cylinder pillars, bed spheres).

**mesa_pd core (new shape, SHAPE_TYPE 6):**
- `src/mesa_pd/data/shape/Cylinder.h` (NEW): finite solid cylinder, **axis = body z-axis**; stores
  `radius_`, `height_` (full length). Correct solid-cylinder inertia (`I_axis = mr²/2`,
  `I_perp = m(3r²+h²)/12`), `getVolume()`, GJK `support()` (cap-rim point), `pack/unpack`.
- Dispatch plumbing: added the Cylinder case to `kernel/SingleCast.h`, the full 7×7 `kernel/DoubleCast.h`
  (regenerated with identical formatting), `mpi/ShapePackUnpack.h` (ghost sync/checkpointing), and
  `common/Contains.h` (`isPointInsideCylinderBF` + `ContainsPointFunctor` overload → used by
  `MovingParticleMapping`, so the cylinder maps onto the LBM grid rotation-aware like the box).

**Contact detection rebased onto GJK/EPA (`ParticleKernels.h`).** `SolidCylinderContactDetection` now
derives from `GeneralContactDetection` instead of `AnalyticContactDetection`:
- **Box×Box is now handled** (GJK/EPA on support functions) — the old defensive skip is removed, so
  box floaters collide with each other and with `FixedBox` pillars. Same for Cylinder×Box,
  Cylinder×Cylinder, Cylinder×Ellipsoid. EPA's contact convention is float-identical to the analytic
  functions (asserted in `tests/mesa_pd/collision_detection/EPA.cpp`), so `LinearSpringDashpot` sees
  the same normal/penetration semantics.
- **Analytic overloads still win** overload resolution (exact-match beats the generic template), so the
  bed-sphere hot path (sphere-sphere/-halfspace/-box) and the battle-tested Box×HalfSpace /
  Box×CylindricalBoundary custom code are unchanged.
- NEW analytic **Sphere×Cylinder** (body-frame radial+axial clamp; interior fallback picks the nearest
  face so the normal never degenerates) — analytic because a log next to thousands of bed spheres
  should not pay GJK per pair per sub-cycle.
- NEW analytic **Cylinder×CylindricalBoundary** (floating log vs solid pillar): pillar axis transformed
  into the log body frame, project-and-clamp fixed-point iteration (radial+axial), same −1.5 LU
  envelop safety clamp + warning as the box version ([CYL-CYL-ENVELOP]).
- Cylinder×HalfSpace (floor/walls) comes from the base's support()-based halfspace overload — exact for
  a finite cylinder. `numParticleTypes` bumped 4→5.

**Config / creation / coupling (app):**
- `SimulationParams.h`: `FloatingCylinderParams` + `FloatingCylinders { FloatingCylinder { … } }` parser.
  Keys: `position`, `radius`, `length`, `axis` (initial world axis direction, default `<0,1,0>` =
  spanwise log; normalised), `densityRatio`, `settlingSteps`, `initialTiltAngleDeg` (about world Y,
  applied after axis alignment), `comOffsetZ`, `pinPosition`, `pinAxes` — all box knobs carry over 1:1.
- `SetupStages.h`: creation (type 4) with rotation aligning body-z to `axis`;
  `interactionRadius = √(r²+(L/2)²)`; included in `blockSyncExtension` / `maxParticleDiameter` (sync
  algorithm choice); `floaterSelector` (initial mapping) + settling release extended to cylinders.
- `BuildTimeloop.h`: `movingObstacleSelector` includes type 4 → the cylinder is a full two-way MEM
  floater (blocks fluid, SimpleBB force, PDF reconstruction), exactly like the box. New
  `floaterCylinder` VTK output: `glyphScale=(2r,L,2r)` + body→world `rotation` tensor (transposed for
  ParaView's Tensor Glyph, same convention as `floaterBox`) + velocity.
- `SimulationObjects.h`: `floaterCylinderUIDs` (indexed-UID + allReduce pattern).

**Buoyancy / righting (`ParticleKernels.h`).** `computeSubmergedVolumesBox()` generalised: boxes and
cylinders share one `clipFloaters_` list (`{uid, isCylinder, halfExtents}`, cylinder = `(r, r, L/2)` on
body z). The plane-fit of the surrounding free surface is unchanged (shape-agnostic shell); the V_sub
clip: boxes keep the analytic 6-tet OBB∩plane option, **cylinders always use the sampling clip** (12³
grid with the out-of-circle points rejected → uniform samples of the cylinder volume; no closed form for
an oblique plane∩cylinder). Centre-of-buoyancy offset → righting torque, `comOffsetZ` pendulum torque,
pinned translation, and the `[BoxBuoyancy]`/`[BoxRightingMoment]` logs (now with `shape=BOX/CYLINDER`)
all work for cylinders; force breakdown logs as `[CylForces]`. `AddExternalForces` treats type 4 as a
floater (density from invMass×volume — `Cylinder::getVolume()` is exact).

**Verified** (no full build possible on this machine): `g++ -fsyntax-only` TUs against the existing
codegen headers compile the entire change set clean, including forced instantiation of the full
`DoubleCast × SolidCylinderContactDetection` dispatch (all 7×7 shape pairs), `ComputeContacts`,
`ComputeLubrication`, `AddExternalForces`, `doDEMSubCycles`, and the setup/timeloop templates.

**New test configs** (`prms/` + `jobs/`, all 72-rank, same domain as the box tests):
- `cylinder_stillwater_gas.prm` / `run_72_cylinder_stillwater.sbatch` — still-water settle of a true log
  (r=7, L=44, axis Y), counterpart to `log_stillwater_gas.prm`. Expect ~50% submersion and a
  roll-NEUTRAL tilt (circular cross-section — contrast with the square log's GM<0 roll).
- `cylinder_two_pillars_gas.prm` / `run_72_cylinder_two_pillars.sbatch` — log (r=3, L=18) swept into the
  two-pillar gap (14 LU) → clogging via the new Cylinder×pillar contact.
- `mixed_floaters_two_pillars_gas.prm` / `run_72_mixed_floaters_two_pillars.sbatch` — sphere + two boxes
  + log swept into the pillars together: exercises every new pair in one run (Cyl×pillar jam, Box×Box
  via GJK/EPA, Box×Cyl / Sphere×Cyl / Sphere×Box pile-up, Box/Sphere×pillar).

---

## Dynamic-H reinstated as an opt-in toggle, re-framed as a CLOGGING test (2026-07-06)

**Context re-read of the two-pillar floater test.** Reviewing `phantom.log`'s per-step `[FloaterForces]`
breakdown (not just the summary tables) surfaced two things:

1. **The settled region's `F_contact+lub` has a persistent, nonzero, DOWNWARD z-component** (~−0.005 to
   −0.011) throughout steps 44000–60000 — not previously examined. Geometry check: pillar clear gap
   (centres y=15/45, radius 8 each) = `30 − 16 = 14` LU; sphere diameter = `20` LU (`radius 10`). **The
   sphere cannot fit through the gap** — it can only wedge against both pillars. Settled position
   `(54.05, 30, ~58)` is at separation ≈0.0007 LU from BOTH pillars simultaneously (exact double contact).
   Per the user, this is intentional: the test is meant to model **clogging** (a body too large for a
   constriction jamming against it), not a free-floating equilibrium. So the mechanical contact/friction
   channel is expected to matter here, not a confound to eliminate.
2. **The `[phantomCorr]` log printed the CLAMPED value, not the raw prediction** — every logged step in
   the settled region showed `F_pred_z == (was)` exactly, meaning the clamp was fully saturated (unclamped
   F_pred ≥ raw MEM) throughout, not marginally exceeded. Fixed: the log now prints
   `F_pred_raw=… F_pred_applied(clamped)=… SATURATED` when they differ, so saturation is visible instead of
   hidden (`ParticleKernels.h`, `AddExternalForces` log line).

Also fixed a latent sign bug found in the same review: the one-sided clamp (`phantomZ =
min(F_pred, hydro_force[2])`) had no floor at 0, so a NEGATIVE F_pred (local trough, `deltaH<0` — e.g. a
wake or antidune trough) would make the subtraction `hydro_force[2] -= phantomZ` ADD to the upward force
instead of leaving it alone. Every test so far kept floaters near the pillar pile-up (`deltaH>0` always),
so this never fired, but it's live for any general run. Fixed to `std::clamp(itp->second, 0,
hydro_force[2])`.

**Dynamic-H reinstated (toggle `useDynamicReferenceHeight`, default false).** Since a clogging body wedged
against the pillars still needs a physically correct hydrodynamic force (mechanical contact resolves the
"does it clog" question; hydrodynamics must still be right so contact isn't confounded by phantom lift),
the dynamic MEM reference height from the 2026-07-01 thread (previously fully reverted — see below) is
brought back as an opt-in toggle, structured cleanly this time:
- `SimulationObjects.h`: `shared_ptr<real_t> memReferenceHeight`.
- `DataInitialization.h`: `HydrostaticDensityFtor` gets a second constructor taking the shared_ptr and
  reading it fresh per call (the constant-height ctor is unchanged, still used for field init).
- `SetupStages.h`: seeds `memReferenceHeight` to `absoluteLiquidHeight`; `setupBoundaryHandling` picks the
  shared_ptr-based functor when the toggle is on, else the constant one (same construction site as before).
- `ParticleKernels.h`: `computeSubmergedVolumes()` accumulates a `wTot`-weighted average of `localSurfaceZ`
  over floaters with a clean fitted ring this step, EMA-blends it into `*memReferenceHeight` (α=0.2) after
  the loop, logs `[DynamicH]` at eval steps. Force-side only (SimpleBB's `rhoHydStat` never touches the
  fluid — confirmed by re-reading `SimpleBB.h` lines 243/269-276 this session), so no PDF/reconstruction
  feedback loop, unlike naive re-derivations of this idea might risk.
- `SimulationParams.h`: `useDynamicReferenceHeight` prm key (`ModelParameters`, default false).

**Test:** `jobs/run_72_dynamich_ab.sbatch` — RUN A baseline (static flat H, `sphere_drop_two_pillars_gas.prm`)
vs RUN B dynamic-H (`sphere_drop_two_pillars_dynamich_gas.prm`), same clogging geometry, both with
`usePhantomForceCorrection false` to isolate the mechanism. Diagnostic, not pass/fail: compare `[DynamicH]`
H(EMA) tracking the pile-up, `F_hydro_z` A vs B at matched submersion (B should show less open-cap
phantom), and `F_contact+lub_z` A vs B (does cancelling the hydrodynamic phantom change how firmly/where
the sphere clogs, given the mechanical wedging is real and expected in both).

**Known limits, carried over from the 2026-07-01 thread (not re-litigated, just re-confirmed applicable):**
dynamic-H is a STATIC/offset correction only — exact for the uncancelled-pile-up phantom, but structurally
cannot remove a genuinely velocity-dependent residual (the suspected Venturi suction), and previously
overshot to 100% submerged once the phantom was cancelled and that residual was exposed. In THIS run the
mechanical contact (clogging) may dominate the vertical force balance regardless — that's the point of
watching `F_contact+lub_z` alongside `F_hydro_z` this time instead of assuming the vertical balance is
purely hydrodynamic.

---

## Option A end-to-end result + one-sided clamp refinement (2026-07-03, phantom.log)

**First end-to-end test (`phantom.log`) — unclamped subtraction OVER-CORRECTED.** Two-pillar floater A/B:
RUN A (OFF) levitated to z≈73 / ~3% submerged (baseline phantom, reproduced). RUN B (ON, unclamped) did NOT
settle at 50% — it **sank to z≈49 / ~100% submerged**. Trace:

| subm | F_pred subtracted | raw MEM_z | corrected MEM_z |
|---|---|---|---|
| 43% | 0.0216 | +0.0257 | +0.0041 (ok) |
| 49% | 0.0216 | +0.0211 | **−0.0005** (F_pred > raw → flips negative) |
| 65% | 0.0169 | +0.0098 | −0.0071 (strongly negative) |

Two combining causes: (1) F_pred marginally EXCEEDS the raw MEM force near 50% (pinned calibration was on the
edge), flipping the corrected force negative; (2) once nudged past ~55% the floater builds downward momentum
and falls into the submerged-regime **−0.02 down-residual** (the secondary artifact the taper leaves
uncorrected) — a trap with no restoring force. The up-phantom had been *preventing* the floater from ever
reaching that trap; removing it unclamped traded levitation for sinking.

**Refinement — ONE-SIDED CLAMP (`AddExternalForces`).** The phantom is strictly an UPWARD artifact, so its
cancellation may only reduce a positive MEM force toward zero, never make it negative:
`phantomZ = min(F_pred, max(0, raw_MEM_z))`, applied only when raw_MEM_z > 0. Structurally prevents the
negative-force tip: at 49% corrected → 0 (not −0.0005), at 52% corrected → 0 with geometric buoyancy
(+0.0009) restoring toward 50%. Removes levitation while the geometric channel supplies the restoring
spring. Re-run `jobs/run_72_phantomfix_ab.sbatch` to confirm RUN B now holds ~50%.

**Bogner free box (`tilt.log`) — settled at −11°, not upright (SEPARATE issue).** The 6:4 box (GM≈+2 ⇒ should
be stable upright at 0°) settled at tiltY≈−11° with a nonzero residual restoring moment M_y≈+0.017 — a free
body should not rest under an unbalanced moment. Unaffected by the phantom work (boxes get no F_pred; this
config has the correction off). Flagged for a separate look at the box righting dynamics.

## Option A phantom-force subtraction + Bogner righting-moment validation (2026-07-03)

**Option A implemented (`usePhantomForceCorrection`, default OFF).** The tapered phantom force F_pred_z (now
computed EVERY step in `computeSubmergedVolumes` and stored in `phantomForceZ_`) is subtracted from the
VERTICAL MEM force on each floating sphere in `AddExternalForces` (`hydro_force[2] -= F_pred_z`), particle
side only — no PDF/fluid feedback, so none of the dynamic-H instability. Keeps horizontal drag and the
vertical fluctuations (only the static double-counted head is removed); self-disables where Δh→0 (no
pile-up). `[FloaterForces]` logs the per-step correction. Toggle in `ModelParameters`.
- End-to-end test: `jobs/run_72_phantomfix_ab.sbatch` — the original two-pillar floater, A) correction OFF
  (baseline, levitates to z≈74 / 1–9% submerged) vs B) ON (`sphere_drop_two_pillars_phantomfix_gas.prm`,
  should sit ~50% submerged at the pile-up surface). PASS = B no longer launches.

**Bogner & Rüde 2013 Fig 5.6 righting-moment validation (setup).** New `[BoxRightingMoment]` log line: for a
box held FIXED at heel α it prints heelDeg, the righting moment M_y = (B−G)×F_buoy (from the analytic clip's
exact submerged centroid), the righting arm GZ = |M|/(V_sub·|g|), and V_sub — one point on the stability
curve M(α). Sweep configs `prms/bogner_righting_6to4_a{00,05,10,15,20,25,30}_gas.prm` (6:4 cuboid: width
b=24 x, height h=16 z, length 2h=32 y; fixed, half-immersed at z=56; analytic buoyancy ON; still water),
job `jobs/run_72_bogner_righting_6to4.sbatch`. Collect M_y vs heelDeg → the stability curve, compare to the
ideal.

**Scope (per user, 2026-07-03): coarse resolution only, goal = the EQUILIBRIUM POSITION** (not the
multi-resolution convergence plot). Primary run: `prms/bogner_righting_6to4_free_gas.prm` +
`jobs/run_72_bogner_free.sbatch` — the 6:4 box FREE-floating (densityRatio 0.5), started 10° off upright; a
wide>tall box is stable upright, so it should right itself: tiltY → ~0°, ~50% submersion, M_y → 0 at
equilibrium. The 7 fixed-angle configs + `run_72_bogner_righting_6to4.sbatch` remain available for the full
stability curve if wanted. The 48 dx / 96 dx resolutions are intentionally NOT built.

---

## FEATURE: analytical (exact) box buoyancy + translation-pin option + tilt test (2026-07-02)

Added an **analytical** submerged-volume / centre-of-buoyancy path for floating boxes as an alternative
to the 12³ point-sampling clip, plus a way to pin a box's translation while leaving rotation free, and a
test that exercises both.

**Analytical clip (`clipBoxAnalytic` in `ParticleKernels.h`).** Closed-form counterpart of the Python PoC
`box_overlap_volume.py` (OBB ∩ half-space), specialised to box-vs-plane so no convex-hull lib is needed:
the oriented box is split into **6 tetrahedra (Kuhn triangulation)** and each tet is clipped against the
fitted free-surface half-space (cases by #vertices-below: 0/1/2/3/4; the 2-in case fans a triangular prism
into 3 tets). Returns V_sub and the world-frame B−G lever, same outputs/frame as the sampling clip, so it
feeds the existing righting-torque path unchanged. Validated in Python against Monte-Carlo (ΔV ≈ 0.02% =
MC noise, centroid to ~1e-4) across rotated-box + tilted-plane and fully-submerged/emerged edge cases.
- Toggle: `ModelParameters { useAnalyticBoxBuoyancy true|false }` (default **false** = legacy sampling).
- The `[BoxBuoyancy]` log always prints a **cross-check**: `method=… [xcheck V_sub: analytic=… sampling=…
  diff=…%]` (sampling still computed at eval steps when analytic is selected), expected to agree <~1%.

**Translation pin (`pinPosition` per box).** `FloatingBoxes/FloatingBox { pinPosition true }` holds only the
box's TRANSLATIONAL DOFs (position reset + linear velocity zeroed each DEM sub-cycle in
`pinBoxTranslations()`), leaving rotation free — unlike the `settlingSteps` FIXED flag, which freezes
everything. Lets us watch the free-rotation equilibrium in isolation, no drift/self-propulsion confound.
Both the old sampling method and free (`pinPosition false`) dynamics remain available — nothing removed.

**Test (`prms/box_tilt_pinned_gas.prm` + `box_tilt_loose_gas.prm`, job `run_72_box_tilt_ab.sbatch`).**
Uniform cube (edge 26, densityRatio 0.5) in still water, started 10° off face-up. A uniform half-submerged
cube is unstable face-up (GM<0) and **stable edge-down (45°)**, so it should roll to **tiltY ≈ 45°**.
RUN A pins translation (clean tilt-to-45° in place); RUN B is fully loose (rolls to ~45° but may drift).
Both use the analytic buoyancy so the cross-check also gets exercised on a live rotating box.

**RESULT (`tilt.log`).**
- **Analytic buoyancy cross-check: PASS.** analytic-vs-sampling V_sub agree to **<0.2% at every step** on a
  live rotating box through all angles (0°→62°) — the 6-tet OBB∩plane clip is correct in production.
- **Loose run (B): PASS.** free to translate+rotate, the cube sank to z≈56, ~46% submerged, and rolled to
  **tiltY ≈ 46–47°** = the stable edge-down orientation. Expected physics confirmed.
- **Pinned run (A): test-setup bug, not code.** All-axes pin at z=62 held the box at only ~22% submerged →
  settled at a different tilt (~62°). The pin feature worked (rotation stayed free); the pin height was
  wrong. **Attempted fix (a):** added a per-axis `pinAxes` mask, pinned `<1,1,0>` (x,y pinned, z free).

### (a) FOLLOW-UP RESULT (`tilt.log`, 2026-07-03): pinning x,y BACKFIRES — reverted recommendation

The z-free depth fix DID work — the pinned box found its own depth (z≈55, ~50% submerged). BUT the tilt then
**spun continuously instead of settling**: tiltY 90°→147°→176°→194°→197°→200°… monotonically rotating
~3°/1000 steps, never stopping. **Cause: pinning x,y removes the rotational damping.** A free floating cube
damps its roll by coupling to sway (horizontal translation) + wave radiation; forbidding horizontal motion
kills that channel, so the buoyancy torque spins the box forever with nothing to dissipate it. Pinning x,y
is therefore NON-PHYSICAL for the righting test.
- **The LOOSE run is the correct tilt test** and it confirms the physics: RUN B settles at **tiltY ≈ 43–45°**
  (edge-down), ~45% submerged, z≈56, with the analytic-vs-sampling V_sub cross-check <0.5% throughout.
- The `pinAxes` mask + translation pin remain in the code (useful for other diagnostics), but the box tilt
  test should use the loose config. box_tilt_pinned_gas.prm kept for reference; not the recommended run.

### (b) FOLLOW-UP RESULT (`pinned.log`, 2026-07-03): submersion taper CONFIRMED live

Re-ran the sweep with the taper active in-solver; it reproduces the offline analysis exactly. Tapered-pred
vs measured: r=10 (43%) 0.0199 vs 0.0198; r=8 (25%) 0.0164 vs 0.0168; r=8 (50%) 0.0138 vs 0.0135; r=8 (75%)
**0.0069 vs 0.0081** (raw was 0.0125 = 54% over → now 15% under); r=5 (71%) 0.0035 vs 0.0064 (small-sphere
under-resolution anomaly, as flagged). Taper works: exact ≤50%, big improvement at high submersion, only the
under-resolved r=5 case remains off. Ready to fold into the Option A subtraction.

---

## DECISIVE: pinned-sphere-at-pillar matrix localizes the phantom (2026-07-02)

`pinned_mat.log` — sphere held **FIXED** at x=53 (the pillar hot-spot), static flat H, hydro-force
submersion scaling DISABLED so the raw MEM force is logged. 2 depths × φ-fill-weighting on/off = 4 runs.
Steady-tail (steps ≥20k) mean **F_hydro_z**:

| run | config | submersion | **mean F_hydro_z** | geometric buoy |
|---|---|---|---|---|
| A | half (z=64), φ-weight **ON**  | ~40%  | **+0.0195 (UP)**   | −0.004 |
| C | half (z=64), φ-weight **OFF** | ~40%  | **+0.0201 (UP)**   | −0.004 |
| B | submerged (z=50), φ **ON**    | 100%  | **−0.0133 (DOWN)** | +0.021 |
| D | submerged (z=50), φ **OFF**   | 100%  | **−0.0135 (DOWN)** | +0.021 |

**Three hard conclusions (this is the clean discriminator — both confounds removed):**

1. **The phantom is REAL and position-locked, NOT a motion/feedback artifact.** Sphere never moves,
   reference never moves, yet the half-submerged sphere feels a steady **+0.02 upward** MEM force for the
   whole run. Kills every "moving body / moving reference / PID chasing itself" story.

2. **It is an INTERFACE effect — the residual FLIPS SIGN with depth** (+0.02 up half-submerged,
   −0.02 down fully submerged). A datum/height mismatch would push the *same* direction at both depths;
   it doesn't. The force is generated where the sphere **crosses the free surface**, not in the bulk —
   which is also why no single dynamic-H reference could fix it (the SLOPE/interface finding, above).

3. **φ fill-weighting is EXONERATED.** A≡C and B≡D to 3–4 sig figs. Toggling the fill-weighting we added
   neither helps nor hurts → not the phantom source; a red herring for this residual.

**Fix direction (now motivated, not a guess): Option A.** The geometric buoyancy channel is correct and
reference-free (Fbuoy behaves properly: −0.004 near half, +0.021 fully submerged). Let geometric buoyancy
own the vertical *static* balance and strip the static vertical component out of the MEM force so MEM
carries only the dynamic (drag/unsteady) part. Verify against the measured residual magnitude (~±0.02).

### Quantitative verification of hypothesis #1 (uncancelled pile-up head) — CONFIRMED

Model: the MEM subtracts a **flat** reference ρ_hyd(z)=1+3·a_z·(z−flatRef), flatRef=56.024. The pile-up
raises the real surface to localSurfaceZ, so an extra head **Δh=localSurfaceZ−flatRef** goes un-subtracted,
leaving a **uniform excess pressure Δp=|a_z|·Δh** on the wetted lower cap. Open surface (gas on top) ⇒
unbalanced ⇒ net up: **F_pred_z = Δp · A_proj**, A_proj = π·h_cap·(2r−h_cap). Closed (submerged) ⇒ uniform
Δp, zero gradient ⇒ F_pred_z→0.

Backed |a_z|=1.018e-5 out of run B (submerged: netF=|a_z|·0.5·vFull=0.02132). Then over run A steady tail
(steps 20k–40k, 21 samples):

| quantity | value |
|---|---|
| mean **F_pred_z** (Δp·A_proj) | **+0.0193** |
| mean **measured F_hydro_z**   | **+0.0195** |
| **ratio** | **1.02** (2% off) |
| submerged (run B) prediction | ≈0 (closed body) — matches the vanished lift |

The mechanism reproduces both the **magnitude** and the **sign flip** with depth. (Step-to-step Pearson
r≈0.45 only, because instantaneous F_hydro also carries drag/LBM fluctuation and localSurfaceZ is a noisy
plane-fit — but the mean is essentially exact.)

**Instrumentation added:** `[PhantomDecomp]` log line in `ParticleKernels.h` (next to `[FloaterBuoyancy]`)
prints Δh, |a_z|, Δp, r_contact, A_proj, openSurface, F_pred_z every eval step, to compare live against the
measured `[FloaterForces]` F_hydro_z. Re-run job: `jobs/run_72_pinned_decomp_ab.sbatch` (pinned, half +
submerged, fill-weight ON — the exonerated OFF variants dropped). **Live re-run confirmed:** mean
F_pred_z=+0.0193 vs measured +0.0195 (RUN A, half); RUN B (submerged) F_pred=0 (closed) vs measured −0.013
(separate secondary residual, NOT the levitation cause).

### Validation sweep — does F_pred = ΔP·A_proj (A_proj∝r²) generalise? (pending)

Before trusting the surgical subtraction (planned Option A: F_hydro_z −= F_pred_z on the particle side only,
downstream of forceMEM, so the fluid/reconstruction is untouched → no dynamic-H feedback, and F_pred→0 on
clean water so it self-disables), validate the formula across size and submersion. Pinned sphere at the
pillar hot-spot (x=53,y=30), density irrelevant for a pinned body. `jobs/run_72_pinned_sweep.sbatch`,
5 runs sharing an r=8/~50% pivot:
- **RADIUS axis** (~50% submerged): r=5, 8, 10 → tests the r² law.
- **SUBMERSION axis** (r=8): ~25%, ~50%, ~75% → tests the A_proj(submersion) shape.
Configs `prms/pinned_sweep_r{05,08,10}_sub50_gas.prm`, `pinned_sweep_r08_sub{25,75}_gas.prm`.
PASS = mean(meas/pred) ≈ 1 over the steady tail for all five.

**RESULT (`predict.log`).** r² law CONFIRMED and formula validated up to half-submerged; a known
over-prediction above ~70% motivated a taper.

| run | subm | A_proj | meas F_hydro_z | raw F_pred_z | ratio |
|---|---|---|---|---|---|
| r=10 (~50%) | 43% | 310 | +0.0198 | +0.0199 | 1.01 |
| r=8  (~25%) | 37% | 193 | +0.0168 | +0.0164 | 1.03 |
| r=8  (~50%) | 51% | 199 | +0.0135 | +0.0145 | 0.97 |
| r=8  (~75%) | 72% | 176 | +0.0081 | +0.0125 | **0.66** (over) |
| r=5  (~50%) | 71% | 68  | +0.0064 | +0.0058 | 1.10 |

- **r² law holds**: across r=5→8→10 pred tracks meas as A_proj spans 68→199→310 (ratio 0.97–1.10).
- **Clean up to ~50% submersion** (25/50%: ratio 0.97–1.03).
- **Over-predicts above ~70%** (r08_s75: 0.66), and the *within-run* trend is a clean monotonic decline
  (ratio 1.09 @54% → 0.83 @62% → 0.65 @72%) → the open-cap A_proj is symmetric about 50% but the
  phantom is not (dies toward the closed-body regime). This motivated the **submersion taper**.

### Submersion taper on F_pred (implemented in ParticleKernels.h [PhantomDecomp])

`F_pred_z = ΔP·A_proj · taper`, `taper(s) = 1 − smoothstep((s − kneeSub)/(1 − kneeSub))`, kneeSub = 0.5,
s = vSub/vFull. Unity for s ≤ 50% (never disturbs the validated regime), smoothsteps to 0 as s → 100%
(closed body ⇒ phantom vanishes, per RUN B). Effect on the sweep (tapered-pred vs measured, the stable
metric — ratios blow up near taper→0 so compare magnitudes): r08_s75 0.0125 (54% over) → **0.0069**
(15% under, vs meas 0.0081); ≤50% cases unchanged. Log now prints `submFrac`, `taper`, `F_pred_raw`, and
the tapered `F_pred_z`.

**Honest caveat:** two runs at ~71% submersion want *different* corrections (r08_s75 needs strong taper,
r05_s50 barely any) → submersion fraction is not the *whole* story; a residual radius/resolution term
remains. But r05 is under-resolved (10 cells across, phantom ~0.006 near the noise floor) where raw AND
tapered both miss — a resolution limit, not a taper failure. kneeSub is a named constant, tunable when more
well-resolved high-submersion runs are available. The taper is PROVISIONAL but strictly improves the
high-submersion regime while leaving the validated low-submersion regime exact.

---

## REVERTED the MEM-reference thread → back to static flat H (2026-07-01)

The pressure-referenced run (`pressref_regression.log`) closed out the whole "adjust the MEM hydrostatic
reference height" line of attack as a **dead end**, and it has been fully removed. Result of the 3-run A/B:

| run | reference | submersion | F_hydro_z | outcome |
|---|---|---|---|---|
| A baseline | flat H=56 | 4.6–12% | +0.018…+0.024 | levitates (phantom up) |
| B pressure-ref | live pressure, H≈61 (53–68) | 96.9–100% | −0.005…−0.035 | **sinks to 100% — same overshoot** |
| C no-pillar regression | dynamic H (surface) | **23.5–100%** | −0.029…+0.009 | **oscillates wildly** |

Two findings killed the approach:
1. **Pressure hypothesis falsified.** The pressure inversion read H≈61 (mean), essentially the SAME as the
   surface geometry (~62–65) — the column near the sphere genuinely IS hydrostatic to ~61–65. The pile-up
   is real static head, not a momentum-supported standing wave. So pressure-ref behaved like surface-ref
   and sank the sphere.
2. **Dynamic H regresses the clean case.** The no-pillar slow case, previously stable at ~48% with static
   H, oscillated 23–100% with dynamic H — the reference tracking the sphere's OWN waterline is a
   destabilizing feedback loop.

**Root diagnosis: it's a SLOPE problem, not an OFFSET problem.** Sweeping H slides F_hydro(submersion) up
or down (+0.02 at H=56, −0.02 at H=65, force-neutral ~60), but the curve's *slope* vs submersion is
destabilizing (F_hydro grows more negative than buoyancy restores). No reference height can fix a slope,
which is why even the physically-correct H≈61 still sank the sphere.

**Reverted (this session's MEM-reference thread, all removed):** `useDynamicReferenceHeight`,
`useDynamicPressureReference`, `hydrostaticReferenceHeight` (params); `SimulationObjects::memReferenceHeight`;
`HydrostaticDensityFtor` back to a const height; ParticleKernels pressure-probe + dynamic-H update +
`pdfFieldID_`; SetupStages/BuildTimeloop back to flat `computed.absoluteLiquidHeight`. Deleted the failed
test configs/jobs (dynamich / localref / pressref). The MEM reference is once again the static flat H the
original working code used.

**Still present (separate, prior mechanism, default OFF, inert):** `useHydrostaticReconstruction` (the
FormerMO reconstruction density reseed, in `src/…/PdfReconstructionManager.h`). It is a different lever
(reconstruction density, not the MEM reference) and untouched by this revert — a candidate for a later
cleanup if we want the tree fully minimal.

**Open question for the next attempt (NOT reference-height):** we still have not isolated whether the
residual vertical MEM force is a numerical artifact or partly real wake downforce, because every run so far
moved both the reference AND the free body at once. The clean discriminator not yet run: a PINNED sphere at
the pillar location, several fixed depths, static reference, developed flow — measure F_hydro_z in
isolation. Leading fix direction (Option A): for floaters, let the reference-free geometric buoyancy own the
vertical static balance and pass only the dynamic part of the MEM force.

The sections below are the now-superseded record of that thread, kept for the reasoning trail.

---

## Phantom-force fix — DYNAMIC MEM reference height H (2026-07-01) [REVERTED — see top]

This is the actual fix for the strong, trapping phantom. It confirms the SimpleBB-flat-H hypothesis
from the reconstruction post-mortem (below) with a constant-H A/B, then replaces the constant with a
live surface tracker.

**Confirmation — constant local-H A/B (`hydrostatic_twopillars_ab.log`).** Two-pillars trap, one build,
identical except the MEM ρ_hyd reference height: RUN A flat global `H=56.02`, RUN B `H=64.8` (the
measured dammed pile-up height near the parked sphere). Only the SimpleBB reference changed — init and
the (already reference-free, local-waterline) geometric buoyancy were untouched.

| | RUN A (flat H=56.02) | RUN B (H=64.8) |
|---|---|---|
| parked z | ~72 (levitating 16 LU above H) | **sinks 73→65→60→55→53** |
| submersion | stuck **3–5%** | **0→32→64→100%** |
| `F_hydro_z` | persistent **+0.020…+0.024** (phantom up) | **flips negative** (+0.014→−0.003…−0.03) |

Correcting *only* the MEM reference height flipped the sphere from levitation to sinking. **That nails
the flat-H MEM reference as the phantom source** — nothing else changed. The geometric buoyancy
(`netF_z = −0.021`, wants to sink) was correct all along; the flat-H MEM phantom (+0.02) was cancelling
it.

**Why a constant is not enough — the overshoot.** RUN B did not settle at ~50% (the `densityRatio 0.5`
float target); it sank to **100% submerged**. The buoyancy log explains it: once submerged, `netF_z`
saturates at +0.0213 (full buoyancy, up), yet the sphere sits centreDepth **+7…+14 below** the surface,
because a fixed `H=64.8` was only right *while the sphere levitated on the hump*. Once it sinks, the true
`localSurfaceZ` drops to ~62, so the fixed 64.8 over-references and pushes **down**. A constant reference
is right at exactly one configuration:
- `H=56`: right far from pillars, too low (phantom **up**) on the hump.
- `H=64.8`: right on the hump, too high (phantom **down**) once sunk.

The only correct reference is the real local surface — which **moves** and which the geometric buoyancy
already computes live (`localSurfaceZ`, plane-fit from the fill field around each sphere). The fix feeds
that same surface into the MEM reference.

**Implementation — live H (toggle `useDynamicReferenceHeight`, default false → legacy flat H).**
- `SimulationObjects.h`: new `shared_ptr<real_t> memReferenceHeight` — the live reference height.
- `DataInitialization.h`: `HydrostaticDensityFtor` now holds a `shared_ptr<real_t>` and reads it fresh
  per call (added a shared-ptr ctor alongside the constant-height one used by init). SimpleBB therefore
  sees the current H every step.
- `SetupStages.h`: seeds `memReferenceHeight` to the configured constant (`absoluteLiquidHeight`, or
  `hydrostaticReferenceHeight` if >0) and passes it to the boundary functor; logs STATIC vs DYNAMIC.
- `ParticleKernels.h`: `computeSubmergedVolumes()` already plane-fits each floater's `localSurfaceZ`;
  now, when the toggle is on, it averages the waterline of floaters that have a CLEAN fitted surface this
  step (ring present, `wTot>1e-6`) and EMA-updates `*memReferenceHeight` (α=0.2) with it. Fully-
  submerged/emerged floaters (no ring) do not contribute — the reference is held rather than snapped to a
  meaningless centre-Z. Runs in the DEM sub-cycles (after this step's FSLBM force, before the next), so it
  applies one LBM step later — negligible, the surface drifts slowly. Logs `[DynamicH]` per eval step.
- `SimulationParams.h`: `useDynamicReferenceHeight` prm key (ModelParameters, default false).
- A/B: `prms/sphere_drop_two_pillars_dynamich_gas.prm` + `jobs/run_72_two_pillars_dynamich_ab.sbatch`
  (RUN A flat H baseline vs RUN B dynamic H; reconstruction reseed left OFF to isolate the reference).

**Scope note:** the reconstruction hydrostatic reseed (`useHydrostaticReconstruction`, below) still uses
its configured constant height; it is a proven ~1e-4 no-op and is left OFF in the dynamic-H config, so no
divergence matters. Only the SimpleBB MEM reference — the proven driver — is made dynamic.

**RESULT (`two_pillars_dynamic.log`).** Dynamic H worked *mechanically* — `[DynamicH]` logged H climbing
58.5→~65 during settling then hovering **62–65** while parked, i.e. it tracked the real waterline live.
And it moved the equilibrium to the right place: `F_hydro_z` now crosses **zero at ~50% submersion** (flat
H put that zero-crossing up at ~5%). Baseline A still levitated (mean F_hydro_z **+0.0226** at 7.5% subm).

| submersion | `F_hydro_z` (RUN B, dynamic H) |
|---|---|
| 40% | +0.005 |
| 41.6% | +0.0016 |
| **54.7%** | **−0.0002 ← zero crossing** |
| 62% | −0.009 |
| 94% | −0.025 |
| 100% | −0.030 |

**But the equilibrium is UNSTABLE — still overshoots to 100%.** Past 50%, `F_hydro_z` keeps growing more
negative and *faster* than the buoyant restoring force (at 62% subm: buoyancy +0.005 up vs F_hydro −0.009
down → net down) → runaway sink, same 100%-submerged endpoint as the constant-64.8 run.

**The real insight — the pile-up is DYNAMIC, not hydrostatic.** This reconciles all three runs: flat H=56
*under*-references (+phantom up, levitates); local surface H≈64 *over*-references (−phantom down, sinks);
truth is **~59–60, in between**. The pillar pile-up is a momentum-supported standing wave (supercritical
Fr=1.31), so part of its 8-LU elevation is DYNAMIC, not static head. Using the fill-field surface geometry
as a *static* reference over-counts the head and injects a spurious depth-growing downforce. **The
reference must track the local hydrostatic PRESSURE, not the surface GEOMETRY.**

---

## Phantom-force fix — PRESSURE-derived dynamic reference height (2026-07-01)

Follow-on to the surface-sourced dynamic H above. Same shared `memReferenceHeight` / same MEM plumbing;
only the SOURCE of the live H changes: instead of the plane-fit surface geometry, invert the MEASURED
ambient pressure.

**Mechanism.** In `computeSubmergedVolumes()`, for each floater sample the fully-LIQUID cells in the
ambient shell just outside the sphere (distC ∈ (r+2, 1.8r], liquid flag, not in another floater's zone)
and read their actual `pdfField->getDensity`. Each inverts the hydrostatic relation for the local datum:
`H = z − (ρ−1)/(3·a_z)`. Averaged over the shell (H is depth-invariant for a hydrostatic column) this
reads the REAL static head — a momentum-supported pile-up carries no extra static pressure, so H lands
below the visible surface, at the ~59–60 the data pointed to, where the MEM force is force-neutral at
~50% submersion. EMA-smoothed (α=0.2) into `memReferenceHeight` exactly as the surface source.

**Code (toggle `useDynamicPressureReference`, only consulted when `useDynamicReferenceHeight` is on):**
- `ParticleKernels.h`: new `pdfFieldID_` member; accum stride 11→13 with `[11] Σ impliedH`, `[12] count`
  over ambient liquid cells; the dynamic-H update picks the PRESSURE target when the toggle is on and any
  ambient liquid was sampled, else falls back to the surface fit, else holds. `[DynamicH]` log tags the
  source (PRESSURE vs surface).
- `SimulationParams.h`: `useDynamicPressureReference` prm key (ModelParameters, default false).

**Tests (`jobs/run_72_pressref_and_regression.sbatch`, three runs on one build):**
- RUN A — two-pillars baseline (flat H) — reference.
- RUN B — `prms/sphere_drop_two_pillars_pressref_gas.prm` — the pressure-referenced fix. Expect
  `F_hydro_z`→0 at ~50% submersion with NO overshoot.
- RUN C — `prms/sphere_drop_slowflow_dynamich_gas.prm` — no-pillar slow control with dynamic H
  (surface source). Regression guard: the surface there is flat, so dynamic H must be a NO-OP —
  must reproduce the static-H control (~48–50% subm, F_hydro_z ~+0.001).

**Status:** implemented, synced; awaiting the run.

---

## Phantom-force fix — depth-consistent (hydrostatic) FormerMO reconstruction (2026-07-01)

Building on the held-static result (the phantom is DYNAMIC, born when the moving boundary uncovers
cells), this traced the FormerMO refill path, confirmed the mechanism with a live diagnostic, ran a
first fix attempt that no-op'd, and landed the working fix.

**Trace.** A cell uncovered by a moving particle (flag `FormerMO`) is reconstructed in
`PdfReconstructionManager::operator()` → `EquilibriumReconstructor::setLocalEquilibrium`
(`Reconstructor.h:139`). Its density comes from `getLocalAverageDensity`: the **unweighted** mean of
the D3Q27 domain (liquid+interface) neighbours' `getDensity`, fallback 1.0. This is **fill-blind** —
in FSLBM an interface cell stores its emptiness in the fill field φ, not in its PDFs, so it reports
ρ≈1 however empty it is. A cell uncovered at/above the waterline thus gets a full-liquid ρ≈1.

**Confirmation run (a)** — `sphere_drop_slowflow_nopillars_gas.prm`, `drop_slowflow.log`. The live
`[RECON-DENSITY-CHECK]` fired 25× during impact, EVERY hit at z≈56.5–57.5 (waterline ≈56) with
`rhoAssigned≈1.0`, `liquidNbrs=0`, all-interface neighbourhood. Since ρ_hyd(z)=1+3·a_z·(z−H) < 1
above the waterline, the MEM reads ρ_assigned − ρ_hyd > 0 → upward lift. Post-release F_hydro_z:
mean **+0.0026**, impact spike **+0.0388** at 3.7% submersion.

**Key surprise:** every flagged cell also had `rhoFillWeighted=0` and **negative `minPhi`**
(−0.004…−0.144) — the waterline interface cells carry NEGATIVE fill (contact-line fill drift), so
Σφ≈0 and any φ-weighted density is degenerate.

**First fix attempt (b₁) — fill-weighted ρ_W = Σφρ/Σφ — CONFIRMED NO-OP.** Run
`fwrecon.log`: the impact spike was **byte-identical** (`+0.0387899` in both runs), mean F_hydro_z
unchanged (+0.0026 → +0.0032, within noise). The Σφ-weighting cannot engage where Σφ≈0 — i.e. on
exactly the cells producing the phantom. Retired.

**Working fix (b₂) — depth-consistent hydrostatic reconstruction.** Re-seed each FormerMO cell to
ρ_hyd(z)=1+3·a_z·(z−H) (the same profile as `HydrostaticDensityFtor`) instead of the fill-blind mean.
f_i^eq ∝ ρ at fixed u and the reconstructor sets pure equilibrium, so scaling every PDF by
ρ_hyd/ρ_assigned lands on f_i^eq(ρ_hyd, u_particle) with velocity preserved exactly. Uses no φ → immune
to the negative fill drift. Near-no-op in the bulk (ρ_assigned ≈ ρ_hyd); at/above the waterline it
removes the ρ≈1 over-density at the MEM's source.

**Code (toggle `useHydrostaticReconstruction`, default false → legacy path unchanged):**
- `PdfReconstructionManager.h`: ctor + `makePdfReconstructionManagerBogner` take
  `useHydrostaticReconstruction`, `hydrostaticAccelerationZ` (a_z), `absoluteLiquidHeight` (H);
  `reseedHydrostatic(x,y,z)` lambda applied right after the reconstructor in both passes (before the
  conserveMomentum read so the booked force uses the corrected density).
- `SimulationParams.h`: `useHydrostaticReconstruction` prm key (ModelParameters, default false).
- `BuildTimeloop.h`: passes the toggle + `(*simObjs->acceleration)[2]` + `computed.absoluteLiquidHeight`.
- A/B: `prms/sphere_drop_slowflow_hydrorecon_gas.prm` + `jobs/run_72_sphere_drop_hydrorecon.sbatch`
  (toggle on; otherwise identical to the (a) baseline).

**A/B RESULTS.**
- Clean slow-flow (no pillars), `hydrorecon.log`: stable, fix active. Both metrics nudged the RIGHT
  way but within noise — mean F_hydro_z +0.00131 → +0.00023, mean submersion 47.5% → 48.9%. Impact
  spike unchanged (+0.0387899). The re-seed only changes ρ by ~1e-4 (the lattice hydrostatic gradient
  is tiny), so this weak-phantom case is a poor amplifier.
- **Two-pillars trap (strong phantom), `hydrorecon_twopillars.log` — NULL.** A (baseline) vs B (fix)
  on one build: meanFz **+0.02063 vs +0.02066**, submersion both trapped at 3–9% (range 0–28%). Near
  step-for-step identical. **The fix does nothing on the case that matters.**

**CONCLUSION: reconstruction over-density is NOT the dominant driver.** It is real (the diagnostic
confirmed it) but the density correction it enables is ~1e-4 — far too small for a +0.02 force. The
strong, trapping phantom is present WITH pillars and weak WITHOUT → it tracks the pillars, not the
uncovered-cell reconstruction. Leading hypothesis: `SimpleBB` subtracts a hydrostatic reference from a
single GLOBAL flat `absoluteLiquidHeight H`; where the pillars dam the flow the local free surface is
ELEVATED above H, so a submerged cell's true pressure exceeds ρ_hyd(z|flat H) → MEM reads ρ_cell −
ρ_hyd > 0 → spurious lift. That is a force-reading defect in SimpleBB, a different locus than
reconstruction — which is why this fix (reconstruction density) left it untouched.

**Status of the ρ_hyd reconstruction fix:** kept as an opt-in toggle (default false). Harmless,
physically correct, stable, directionally right on the clean case — but not the lever for the strong
phantom. Next investigation: SimpleBB's hydrostatic reference vs the local (spatially-varying) surface.

---

## Phantom-force investigation — held-static discriminator (2026-06-30)

Goal: decide whether the phantom lift is a **static** effect (missing gas pressure on an emerged cap
= "point 1", or hydrostatic-staircase error = "point 2") or a **dynamic** effect of the moving
boundary. Two pinned spheres (r=8) in STILL water (U=0, `useHardcodedUnits true`), never released
(`settlingSteps 1e8`), 15000 steps each in one allocation. submFrac OFF, φ-weighting ON. Configs:
`sphere_held_submerged_gas.prm` (centre z=45, fully submerged) and `sphere_held_partial_gas.prm`
(centre z=56, ~54% submerged). Output: `held_static.log`.

| run | submersion | geometric `F_buoy+grav_z` | **MEM `F_hydro_z`** (mean±std, n=30) | verdict |
|---|---|---|---|---|
| A submerged (z=45) | 100% | +0.00753 | **−0.00181 ± 0.00027** (~6.7σ, real) | small DOWNWARD residual |
| B partial (z=56) | ~54% | +0.00056 | **+0.00020 ± 0.00082** (~0.25σ, noise) | no net static force |

FormerMO conversions were **0 on every report in both runs** (pinned sphere → no moving boundary);
`[FILL-WEIGHT]` fired 120×, 0 dry links.

**Findings.**
1. **Point 1 (missing gas pressure on the emerged cap) is REFUTED.** The decision table predicted
   *A≈0 & B>0* for that mechanism. Instead the **partial body shows no net static vertical force**
   (mean = ¼σ, pure oscillation), and the only real residual is on the **fully submerged** body. A
   missing gas down-force would push the partial sphere *up*; it does not.
2. **A small dose of point 2 exists, but it is downward and minor.** The fully-wetted sphere should
   give MEM `F_hydro_z ≈ 0` at hydrostatic rest (the hydrostatic density is subtracted link-by-link
   so buoyancy is not double-counted — buoyancy is the separate geometric +0.00753). The leftover
   **−0.00181 (≈24% of buoyancy, downward)** is the discretization error of that subtraction:
   `rhoHydStat` evaluated at the **cell center** instead of the wall-link intersection on the curved
   surface. It shrinks into the noise for the half-submerged case (fewer wetted links, 3× interface
   noise).
3. **The big UPWARD phantom is DYNAMIC, not static.** With the sphere pinned (zero FormerMO
   conversions) the strong upward lift seen in the drop runs **disappears entirely**. It is therefore
   generated by the **moving boundary**: each cell uncovered as the sphere moves is reconstructed
   (EquilibriumRefilling/Bogner) toward ρ≈1 at a depth where the correct density is ρ_hyd(z) > 1,
   injecting spurious upward momentum per uncovered cell — the **FormerMO reconstruction** mechanism.

**Decision / next levers.**
- Drop the gas-pressure / static-buoyancy line of fixes; it is not the source.
- Primary fix: **depth-consistent reconstruction of uncovered (FormerMO) cells** — refill to ρ_hyd(z)
  instead of ρ≈1, so a freshly-uncovered cell does not seed an upward momentum kick.
- Cheap secondary win for the −0.0018 static residual: evaluate `rhoHydStat` at the **boundary
  link-intersection** rather than the cell center in `SimpleBB.h`.
- φ-weighting and the submFrac retirement remain unchanged.

---

## Phantom-force investigation — two-scenario drop diagnostic (2026-06-29)

Goal: quantify the phantom hydrodynamic lift on a floating sphere and the number of "unknown"
(FormerMO) cells Bogner reconstruction faces. For this study the **submersion scaling in
`ParticleKernels.h` was DISABLED** (`hydro_force *= submFrac` commented out) so the raw MEM force
reaches the particle, and `PdfReconstructionManager` was instrumented with a `[FORMERMO-UNKNOWNS]`
log (total unknown cells + how many of each cell's D3Q27 neighbours are themselves unknown).

**Setup.** r=10 sphere, started ~5 LU above the z≈63 waterline (centre z=78), held FIXED 30000 steps
then released. `useHardcodedUnits true` pins gravity/viscosity/σ so the two scenarios share identical
drop physics. Configs: `sphere_drop_two_pillars_gas.prm`, `sphere_drop_slowflow_nopillars_gas.prm`.

| | Scenario 1 (two pillars, U=0.02) | Scenario 2 (no pillars, U=0.005) |
|---|---|---|
| F_hydro_z at ~3–4% submersion | +0.017 | **+0.039** (step 32000, ~2× weight) |
| F_hydro_z at ~16% submersion | +0.044 | — |
| sphere outcome | **trapped at 2–27% submersion**, never settles | punches through, settles **~43% submerged** |
| steady F_hydro_z | persistently large + | small (+0.006…+0.009), restoring (−) when deep |

**Findings.**
1. **The phantom is intrinsic.** It survives pillar removal *and* a 4× flow reduction — at 3.7%
   submersion the lift is still +0.039 (≈2× weight) where it should be ≈0. So it is a free-surface
   MEM **fill-blindness** artifact, not a pillar/flow effect. Pillars + fast flow merely *trap* the
   floater in the low-submersion shell where the phantom dominates, turning a transient into a
   persistent wrong equilibrium.
2. **Unknowns are well-conditioned, not clustered.** Peak ≈1200 cells/step (scenario 1) / sustained
   ≈1100 cells/step while bobbing (scenario 2), but in BOTH runs ~95%+ have **zero** unknown
   neighbours and the **max is ≤3**. Branch-C ("unclassified FormerMO neighbour → interface
   fallback") is tiny (11–15 cells). The problem is the **force/density per near-boundary cell**, not
   unknown-on-unknown reconstruction.
3. **The φ-weighting fix is already LIVE but insufficient.** `SimpleBB.h` already scales each link's
   force by φ(x_b) (Bogner Eq. 5.9; `[FILL-WEIGHT]` fired 600×, interface cells weighted by φ≈0.52),
   wired in at `AntidunesBoundaryHandling.impl.h:102`. Yet with submFrac OFF the phantom persists, and
   even in clean water the sphere floats ~43% vs the ideal 50% (a ~7% residual upward bias matching
   the steady +0.006…+0.009 F_hydro_z). **Production stability therefore still rests on the disabled
   submFrac geometric scaling**, not on φ-weighting alone.

**Decision.** Keep **φ-weighting (SimpleBB, Eq. 5.9) as the sole MEM correction**; the submFrac
geometric scaling in `ParticleKernels.h` stays disabled (reworded from a temporary diagnostic to a
deliberate retirement — φ-weighting is the physical fix, submFrac was a band-aid). Accepted
behaviour: in clean water the floater settles ~43% vs the ideal 50% (small residual lift); the
pillar/flow trapping case remains imperfect.

**Residual still open (not pursued for now).** The ~7% equilibrium bias likely comes from LIQUID
links (φ=1, unweighted) whose reconstructed density is ρ≈1 regardless of how short the local liquid
column is (cf. `[RECON-DENSITY-CHECK]`), and/or the flat hydrostatic reference on a raised surface.
A z-force decomposition by near-boundary cell type (liquid vs interface) would localize it if revisited.

---

## Physics Analysis

### MEM three-phase correction (Eq. 5.9)

The Momentum Exchange Method force on particles is computed in
`src/lbm_mesapd_coupling/momentum_exchange_method/boundary/SimpleBB.h` (aliased as `MovingObstacle_T` in `AntidunesBoundaryHandling.h:72`).
It implements only the standard two-phase formula (Ladd 1994 / GIMEM Wen 2014):

```
Δj_q = c_q · (f̃ + f')   [+ Galilean correction]
```

When a particle boundary link crosses an **interface or gas cell**, the correct three-phase formula (Eq. 5.9, Bogner & Rüde 2013) is:

```
Δj_q = c_q · φ(x_b) · [f̃(x_b,t) + f(x_b,t+1)]  +  2·c_q · [1−φ(x_b)] · f_q^eq(ρ_g, 0)
```

`SimpleBB` has no access to the fill-level field and applies no such correction.

**Effect in this simulation:** The gas pressure is set to `ρ_gas = 1` (`BubbleModelConstantPressure(real_c(1))` in `Antidunes.cpp:100`), matching the LBM reference liquid density.
With `ρ_gas = ρ_liquid = 1`, the static pressure terms in the error cancel exactly.
The residual error is `O(Ma · (1−φ))` — a velocity-dependent second-order term that is negligible at the low Mach numbers of this simulation.
The missing correction would be significant only with a real air-water density ratio (`ρ_gas/ρ_liquid ≈ 1/1000`).

---

### Buoyancy from submerged volume — spheres and boxes (FIXED)

**Problem (original):** `AddExternalForces::operator()` applied `ρ_f · V_full · g` for all particles,
overestimating buoyancy for partially submerged floaters.

**Fix (`ParticleKernels.h`):** Each step, the submerged volume `V_sub` of every floater is
precomputed and the vertical force becomes `(ρ_P · V_full − V_sub) · g_z`. Non-floaters are
always fully submerged (`V_sub = V_full`) so they are unchanged.

```cpp
// X/Y: full-volume driving force (pressure-gradient body force acts on the whole particle).
// Z:   gravity (rho * vFull) minus correct buoyancy (vSub × g, fluid density = 1).
Vector3< real_t > bodyForce(
   vFull * (*acceleration_)[0],
   vFull * (*acceleration_)[1],
   (rho * vFull - vSub) * (*acceleration_)[2]);
```

- **Sphere** (`computeSubmergedVolumes`): contact-ring plane-fit of the local free surface
  (weighted normal + centroid over interface cells on the sphere surface) → spherical-cap volume.
- **Box** (`computeSubmergedVolumesBox`): plane-fit of the surface in a shell around the box,
  then a 12³ point grid inside the rotated box is clipped against the fitted plane → `V_sub` and
  the centre-of-buoyancy offset (lever arm for the righting torque). Independent of the box-interior
  fill (which becomes obstacle once the box is mapped).

**Cavity-proof surface bound (the sink fix).** A deeply-submerged floater can sit in an entrained
gas cavity; the surrounding "interface" cells are then cavity walls, so the plane-fit under-counts
`V_sub` → buoyancy collapses → the floater sinks to the floor and stays. Fixed with a surface-height
reference `surfRef = max(absoluteLiquidHeight, localSurfaceZ)`, where `localSurfaceZ` is the live
gradient-weighted waterline at the floater:

- floater **top ≤ surfRef** → fully submerged → `V_sub = V_full` (cavity-proof),
- floater **bottom ≥ surfRef** → fully emerged → `V_sub = 0`,
- straddling → the local plane-fit above.

The `max` keeps the reference at the far-field level when a hugging cavity drags the local estimate
down, but still follows a genuinely risen/displaced surface. Validated: the giant sphere and the
giant box both bounce and settle near half-submerged instead of sinking.

The `[FloaterBuoyancy]` / `[BoxBuoyancy]` logs report submersion, the live `localSurfaceZ`, the net
buoyancy+gravity force, and (box) the tilt and B−G lever arm.

---

## Floating box — full two-way LBM coupling

The floating box (type 3) now has the same capabilities as a floating sphere.

- **Mapping (`BuildTimeloop.h`):** the moving-particle mapping selector was generalised from a
  sphere-only filter to `movingObstacleSelector` = *all spheres + floating boxes* (fixed pillars
  excluded). The box is therefore mapped as an LBM obstacle and exchanges momentum via the MEM,
  giving it real hydrodynamic drag (which also damps its rolling — the giant box settles, the small
  one tumbles indefinitely).
- **Buoyancy + righting torque (`ParticleKernels.h`):** `computeSubmergedVolumesBox` feeds both the
  vertical buoyancy force and a righting torque from the centre-of-buoyancy offset
  (`buoyancyCenterOffset_`), so the box rolls/pitches under buoyancy like a real log (GM<0 → rolls).
- **Setup (`SetupStages.h`):** `FloatingBox` creation (`Utility.h` shape with
  `updateMassAndInertia(densityRatio)`), optional `initialTiltAngleDeg` seed, FIXED-on-creation +
  `releaseSettledFloatingBoxes` mirroring the sphere settling mechanism.

### Block-spanning floaters → `SyncGhostOwners`

`blockSyncExtension` and the sync-algorithm choice in `SetupStages.h` now derive from the largest
floater extent (sphere radius / box circumscribed half-diagonal). A floater large enough to span
2+ blocks (e.g. the giant sphere/box) automatically switches from `SyncNextNeighbors` to
`SyncGhostOwners`, which correctly ghost-copies one body across all overlapped blocks.

### Box VTK orientation (`BuildTimeloop.h`, `src/mesa_pd/vtk/WriteOutput.h`)

The `floaterBox` point output writes `edgeLength`, the full 3×3 `rotation` matrix, and `velocity`.
A new `writeOutput` overload for `std::array<real_t,9>` in `WriteOutput.h` lets the rotation matrix
be written as 9 components for a ParaView **Tensor Glyph** (Glyph Type = Box, X/Y/Z Length = edge
lengths, **Extract Eigenvalues OFF**). `getMatrix()` is body→world, but ParaView's Tensor Glyph
applies the tuple with the opposite ordering, so `SelectParticleRotationMatrix` emits the **transpose**
(Rᵀ) — visualization-only; the buoyancy/torque math uses `getMatrix()` directly and is unaffected.

---

## Reconstruction model toggle — Bogner vs all-interface

`makePdfReconstructionManagerBogner` (`BuildTimeloop.h`) drives the FormerMO PDF reconstruction with
`EquilibriumReconstructor`, gated by `ModelParameters.useBognerReconstruction` (`SimulationParams.h`):

- **true** → Bogner & Rüde (2013) per-cell state selection: each uncovered (FormerMO) cell is
  classified by its D3Q27 non-obstacle neighbourhood and reconstructed as gas / liquid / interface.
- **false** → every FormerMO cell is reconstructed uniformly as **interface** (no fill=1
  over-promotion band) — the path that never drifted/NaN'd in the box tests.

`PdfReconstructionManager.h` carries per-branch FormerMO→state counters and an `[OVERPROMOTE-CHECK]`
diagnostic (waterline cells promoted to liquid that should be interface).

---

## Free-surface / MEM NaN guards and diagnostics

- **Stale near-boundary flag (`AntidunesBoundaryHandling.h`):** `BoundaryHandling::refresh()` only
  updates the near-boundary flag on domain (liquid/interface) cells, so when a near-boundary
  interface cell converts to GAS the flag is left stale → the treat-loop would hand that gas cell
  (invalid `-nan` PDF) to `SimpleBB`. The boundary sweep now detects `(gas && nearBoundary)`, warns
  `[STALE-NEARBOUNDARY]`, and clears the flag *before* the treat-loop — removing the cause of the
  SimpleBB dry-link NaN. (0 DRY-LINK / 0 LOW-RHO in the validated runs.)
- **`SimpleBB.h` `[DRY-LINK]`:** finite-check on `pdf_old`; on a non-finite read it logs the cell,
  decoded flags, and raw flag value, then skips the link (the bounce-back victim, not the origin).
- **`ReconstructInterfaceCellABB.h` `[LOW-RHO]`:** before `u /= rho` in the interface
  reconstruction, warns and dumps all 27 PDFs (with directions) when `|rho|` is near zero / non-finite
  — pinpoints a collapsing interface cell before it produces a non-finite velocity.
- **`CurvatureSweep.impl.h`:** the generic "Invalid normal" warning now reports the interface cell's
  global coordinates, its D3Q27 obstacle-neighbour count (→ whether it is a contact-line cell), and
  the direction toward the solid — locating degenerate-normal cells at the particle/free-surface
  junction (no surface-tension force there with wetting disabled).

---

## Experiment configs (`prms/`, `jobs/`)

New `.prm` files and matching 72-rank `.sbatch` scripts (domain `<128,60,180>` /
cellsPerBlock `<32,20,30>`), all with `useBognerReconstruction true`:

| prm | scenario |
|---|---|
| `giant_sphere_gas` / `_free` / `_drop` | giant sphere (r=20) at/over the surface, two pillars / free / dropped |
| `box_free_gas` | **giant** free-floating log `12×36×12` (1:3:1), sweeps the periodic domain |
| `box_pillar_gas` | log wedged between two cylinder pillars |
| `giant_box_drop_gas` | giant 24³ cube dropped from height |
| `giant_box_two_pillars_gas` | giant log placed clearing two pillars (the giant_sphere setup, box-swapped) |
| `ten_floaters_gas` | ten floating spheres |

Each `run_72_*.sbatch` builds + runs via `antidunes.sh`, archiving prior `vtk-out`/`eval`/`snapshots`.

---

## NaN crash fix — PDF reconstruction at the free-surface contact circle

**Symptom:** Simulation crashed with NaN after a few hundred timesteps when a floating sphere
started at or near its equilibrium position at the free surface.

**Root cause:** When the sphere moves and uncovers a contact-circle cell (where the sphere
surface meets the free surface, fill ≈ 0.5), `PdfReconstructionManager` must synthesise fresh
PDFs. `GradsMomentApproximationReconstructor` (Chikatamarla et al., EPL 2006; Dorschner et al.,
JCP 2015) computes `f[i] = f^eq + Pneq` where `Pneq` uses a velocity gradient from the
D3Q27 neighbourhood. In the free-surface context, ~⅓ of those neighbours are gas cells that
pass `isDomain()` but carry near-zero PDFs. This underestimates `ρ_avg` and produces a large
spurious `Pneq` term that pushes individual `f[i]` negative. After one streaming step
`ρ = Σf[i] < 0` and the cumulant collision (ω ≈ 1.998) diverges to NaN.

Note: `GradsMomentApproximationReconstructor` has no fallback for incomplete neighbourhoods —
it was designed for solid-wall contexts where "missing" neighbours are absent, not gas cells
with physically zero PDFs.

### Fix 1 — `BuildTimeloop.h`: switch to `EquilibriumReconstructor`

```cpp
// two-arg overload of makePdfReconstructionManager uses EquilibriumReconstructor internally
lbm_mesapd_coupling::makePdfReconstructionManager< PdfField_T, BoundaryHandling_T >(
   blockForest, pdfFieldID, handlingID, particleFieldID, particleAccessor,
   FormerMO_Flag,
   walberla::free_surface::flagIDs::interfaceFlagID,  // see Fix 2
   strictlyConserveMomentum)
```

`EquilibriumReconstructor` sets `f[i] = f^eq(ρ_avg, u_particle) ≥ 0` for all 27 directions —
guaranteed positive regardless of neighbourhood quality.

### Fix 2 — `BuildTimeloop.h`: mark uncovered cells as `interfaceFlagID`

Contact-circle cells have fill ≈ 0.52 (set by `ObstacleFillLevelSweep` when active, or frozen
from the initial fill otherwise) and are physically interface cells. With `liquidFlagID`,
`StreamReconstructAdvectSweep` performs a full 27-direction streaming step, pulling f=0 from
gas neighbours and corrupted PDFs from obstacle neighbours → `ρ < 0` → NaN.

With `interfaceFlagID`, `StreamReconstructAdvectSweep` *reconstructs* gas-facing PDFs from the
atmosphere density (ρ=1, u=0) rather than streaming zeros. `CellConversionSweep` then converts
the cell to LIQUID or GAS naturally based on fill level within the same timestep.

This flag also applies to submerged and above-surface uncovered cells (always all FormerMO cells),
but those cases are self-correcting within one FSLBM step via `CellConversionSweep`.

---

## Modified waLBerla library files (current diff)

| File | Change |
|---|---|
| `src/lbm_mesapd_coupling/.../reconstruction/PdfReconstructionManager.h` | Bogner per-cell state-selection reconstruction with the `useBognerReconstruction` toggle; FormerMO→state counters; `[OVERPROMOTE-CHECK]` diagnostic |
| `src/lbm_mesapd_coupling/.../boundary/SimpleBB.h` | `[DRY-LINK]` finite-check on `pdf_old` (decodes flags, skips the link) |
| `src/lbm/free_surface/dynamics/functionality/ReconstructInterfaceCellABB.h` | `[LOW-RHO]` guard + 27-PDF dump before `u /= rho` |
| `src/lbm/free_surface/surface_geometry/CurvatureSweep.impl.h` | invalid-normal warning enriched with global coords, obstacle-neighbour count, contact-line detection |
| `src/mesa_pd/vtk/WriteOutput.h` | `writeOutput` overloads for `std::array<real_t,9>` (rotation matrix → ParaView Tensor Glyph) |
| `src/mesa_pd/data/shape/Cylinder.h` | NEW finite solid cylinder shape (SHAPE_TYPE 6, axis = body z): inertia, volume, GJK support(), pack/unpack |
| `src/mesa_pd/kernel/SingleCast.h`, `kernel/DoubleCast.h` | Cylinder added to the shape dispatch (DoubleCast regenerated as the full 7×7 switch) |
| `src/mesa_pd/mpi/ShapePackUnpack.h` | Cylinder case for ghost sync / checkpoint shape (de)serialisation |
| `src/mesa_pd/common/Contains.h` | `isPointInsideCylinderBF` + `ContainsPointFunctor` Cylinder overload (LBM cell mapping) |
| `src/mesa_pd/data/shape/CMakeLists.txt` | register `Cylinder.h` |

These are detailed in the sections above. (Earlier exploratory edits to `CurvatureModel.impl.h`,
`ObstacleFillLevelSweep.impl.h`, `ObstacleNormalSweep.impl.h`, and `ExcessMassDistributionSweep`
were reverted and are no longer part of the diff.)

---

## Meeting 04

### Goal

Fix NaN crashes when a floating sphere touches the free surface, understand the root cause,
test three initial positions (submerged, at surface, in gas), and implement a `settlingSteps`
mechanism for initialisation stabilisation.

---

### Findings from Testing Three Starting Positions

| Variant | Start z | Result | Notes |
|---------|---------|--------|-------|
| `two_cylinder_floater.prm` | 30 (submerged) | ✅ 100 000 steps | Sphere rises through liquid; contact circle grows from 0 — fill-level gradient always well-defined |
| `two_cylinder_floater_surface.prm` | 56 (at surface) | ❌ step 581 | Full-radius contact from t=0; settlingSteps=500 improves from step 12 → 581, but fill-level drift during settling corrupts cells that crash 81 steps after release |
| `two_cylinder_floater_gas.prm` | 90 (in gas) | ❌ step 2668 | Sphere falls ~2259 steps before surface contact; contact circle grows ~409 steps then fails — same NaN mechanism as surface variant |

---

### 1. `SimulationParams.h` — `settlingSteps` in `FloatingSphereParams`

Added `settlingSteps` field (default 0): number of LBM timesteps the sphere is held FIXED
before being released to float freely. Non-zero values let the FSLBM establish a stable
fill-level distribution around the stationary contact circle before the sphere first moves.

---

### 2. `SetupStages.h` — FIXED flag on creation + `releaseSettledFloatingSpheres`

If `fsp.settlingSteps > 0`, the FIXED flag is set on the new particle at creation.
`releaseSettledFloatingSpheres()` is called each timestep and unsets FIXED when
`currentStep == fsp.settlingSteps`.

---

### 3. `Antidunes.cpp` — call `releaseSettledFloatingSpheres`

Called at the top of the main time loop before the FSLBM step.

---

### 4. New `.prm` files — three floating sphere test cases

- **`prms/two_cylinder_floater.prm`**: sphere z=30 (submerged), `settlingSteps 0` — ✅ runs 100 000 steps
- **`prms/two_cylinder_floater_surface.prm`**: sphere z=56 (at surface), `settlingSteps 500` — ❌ crashes at step 581
- **`prms/two_cylinder_floater_gas.prm`**: sphere z=90 (in gas), `settlingSteps 0` — ❌ crashes at step 2668

---

## Meeting 03

### Root Causes Diagnosed

#### 1. `blockSyncExtension` too small

`SyncNextNeighbors` ghost-copies a particle to a neighboring block only if the particle center is within `blockSyncExtension` of the block boundary. The hardcoded `2.5` meant a large floating sphere could physically overlap a block without being synced there → LBM saw stale PDFs → NaN hydro force.

**Fix (`SetupStages.h`):** `blockSyncExtension` is now `max(2.5, largest_floater_radius)`.

#### 2. MPI count mismatch in `computeSubmergedVolumes`

Different vector sizes on different ranks caused `MPI_Allreduce` with mismatched `count` → undefined behaviour → silent corruption → eventual NaN.

**Fix (`ParticleKernels.h`):** Pre-size all vectors to `simParams->floatingSpheres.size()` (same on all ranks). Only the owning rank writes its UID/position; `allReduceInplace(SUM)` broadcasts to all.

#### 3. `phiSum` overcounting

Each block computed `phi_avg = localSum / localCount` then `allReduce SUM` added partial averages — overcounting when ring cells are split across blocks.

**Fix (`ParticleKernels.h`):** Raw numerator and denominator accumulated per Z-slice into flat vectors, both reduced globally, divided after. `phiSum` feeds the spherical-cap formula `V_sub = π·h²·(3r−h)/3`.

#### 4. Sphere top half initialised in gas zone

Sphere center at Z=52, radius=10: top half (Z=52–62) above free surface (≈Z=52.5). `MovingObstacle` cells in gas conflicted with the free-surface PDF scheme → NaN from timestep 0.

**Fix:** Position sphere fully submerged from initialisation.

---

## Meeting 02

### Goal

Fix a NaN quaternion crash when running with a floating sphere obstacle. Replace box pillar with a cylindrical pillar supported in both DEM and LBM.

---

### 1. `ParticleKernels.h` — NaN contact normal guard

Added NaN check in `ComputeContacts::operator()` before applying contact response. Two overlapping particles at exactly the same position produce a zero-length contact normal → NaN after normalization → quaternion abort. Guard logs and skips the contact instead.

---

### 2. `ParticleKernels.h` — `SolidCylinderContactDetection`

Subclass of `AnalyticContactDetection` that overrides `(Sphere, CylindricalBoundary)` contact
detection with exterior (solid cylinder) math. The base class uses hollow/interior math.
`ComputeContacts::acd_` changed to `SolidCylinderContactDetection`; all other shape pairs
fall through to inherited base-class methods unchanged.

**Why subclass the functor:** `DoubleCast` dispatches by `SHAPE_TYPE` and slices shape subclasses,
so overriding the functor is the correct extension point. `using AnalyticContactDetection::operator();`
is required to prevent C++ name hiding of all base-class overloads.

---

### 3. `SimulationParams.h` — `CylinderPillarParams`

Added `CylinderPillarParams` (center + radius) and `parseCylinderPillars()`. `SimulationParams`
now holds both `pillars` (box) and `cylinderPillars` (cylinder) vectors.

---

### 4. `SetupStages.h` — cylinder pillar creation

Radial exclusion check added to `createParticleFct` for cylinder pillars. DEM creation via
`createCylinder()` (using `CylindricalBoundary` shape with `SolidCylinderContactDetection`).

---

## Meeting 01

### Goal

Add a fixed solid box pillar obstacle to the Antidunes simulation.

---

### 1. `Utility.h` — `FixedBox` and `createBox()`

`FixedBox` subclass of `mesa_pd::data::Box` overrides `updateMassAndInertia` to set
`invMass = 0` (same as a wall), so contact stiffness with the pillar is identical to
sphere–wall contacts. `createBox()` sets GLOBAL / INFINITE / FIXED / NON_COMMUNICATING flags.

**Why `INFINITE` flag:** prevents Box↔HalfSpace pairs from reaching `AnalyticContactDetection`
(no overload for that combination). Box↔Sphere is still included because Sphere is not INFINITE.

---

### 2. `SetupStages.h` — `setupParticles()`

`createParticleFct` excludes spheres overlapping the pillar column (nearest-point XY distance).
`createBox()` called after `createPlane()` calls.

---

### 3. `SimulationParams.h` — `parsePillars()`

`PillarParams` with center + halfExtents. `parsePillars()` reads zero-or-more `Pillar` sub-blocks
from the `Pillars` config block. `SimulationParams` holds `std::vector<PillarParams> pillars`.

---

## What Was NOT Changed (and Why)

### Particle mapping (LBM–DEM coupling)
The moving-particle mapping uses `movingObstacleSelector` = spheres + floating boxes (type 3).
**Fixed** box and cylinder pillars are still excluded from MEM mapping — they never conflict with
the NoSlip flags set by the boundary handling. (The floating box, unlike the fixed pillars, *is*
now mapped — see "Floating box — full two-way LBM coupling" above.)

### Fill level inside pillar
`initFlagsFromFillLevels` explicitly skips cells where `handling->isBoundary()` is true.
NoSlip cells keep their flag regardless of fill level.

### Zero velocity in pillar cells
`ZeroSetter` already zeros PDFs in obstacle cells before every VTK write.
