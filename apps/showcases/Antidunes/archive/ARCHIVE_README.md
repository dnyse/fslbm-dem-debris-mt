# Archive

Code and configs built and run during the thesis work in this directory, but that did not make it
into the final thesis results. Kept for reference rather than deleted, mirroring the thesis repo's
own `archive/04_results/toy-bridge/` precedent (see that repo's `CLAUDE.md`). Excluded from the
CMake build (not referenced by any live `CMakeLists.txt` target or `#include`).

## tsunami_pole_forest_driftwood/

`tsunami.cpp` + `TsunamiDataInitialization.h`: a dry-domain/dam-break reproduction of Watanabe et al.
(2023)'s pole-forest/driftwood tsunami-flume experiment, reusing Antidunes' own setup/timeloop/DEM
machinery unchanged except for its own initial-condition file. Built and run four times on the
cluster (`jobs/run_72_toy_pole_forest_driftwood.sbatch`, `run_72_tsunami_pole_forest_driftwood.sbatch`,
`run_288_tsunami_pole_forest_driftwood_proportional.sbatch`,
`run_576_tsunami_pole_forest_driftwood_proportional_elongated.sbatch`, with their `prms/` configs),
but never written up or cited anywhere in the thesis chapters/appendix. Watanabe et al. (2023) is
cited elsewhere in the thesis only for its driftwood DEM contact parameters (restitution/friction),
used by the unrelated, already-in-the-thesis 250-floater drop test (`DropTest.cpp`).

[Claude, 2026-10-29] Archived at the author's request. `CMakeLists.txt`'s `tsunami` executable target
was removed. CORRECTION: initially stated "not `#include`d by anything else" -- missed
`tests/showcases_antidunes/TsunamiDamBreakTest.cpp`, a dedicated unit test for
`TsunamiDataInitialization.h`'s pure geometry functions, which broke when this archival moved that
header. Now also archived here (`tests/TsunamiDamBreakTest.cpp`, include path fixed to
`../TsunamiDataInitialization.h`); its `ANTIDUNES_TsunamiDamBreak` CMake target removed from
`tests/showcases_antidunes/CMakeLists.txt`.

## sealed_cell_guard/

`SealedCellGuard.h` (`ModelParameters.closeSealedCellsMinNeighbours`, `CloseSealedCellsSweep`) --
an experimental guard against sub-cell fluid "channels" sealed between touching rods (introduced
2026-09-21). Superseded by `cylinderContactSkin`/`boxContactSkin` (implemented 2026-09-30), which
*prevents* the sealing by keeping a minimum resolvable gap between bodies' DEM contact detection,
rather than detecting and closing already-sealed cells afterward. Per `THESIS_PLAN.md`
(2026-09-25): "`closeSealedCellsMinNeighbours=24` delays but does NOT eliminate the sealed-cell
mechanism -- and reveals the guard is checking the wrong signature" (total enclosure count, not
specifically "pinched between two separate bodies"). Contact skin became the actual production
approach from 2026-09-30 onward and is the one used in every thesis-cited run since; the guard was
never combined with it in any prm. Confirmed: none of the 9 prms that enabled the guard
(`*_guard21_gas.prm`/`*_guard24_gas.prm`/`*_guard24_rodcap_gas.prm`) are cited in the thesis, and
none were ever git-tracked. Moved here along with its 6 dedicated sbatch jobs, its own unit test
(`tests/Job4255679SealedCellTest.cpp`, which specifically validated whether the guard's threshold
would have caught a real crash -- the one test in the suite that was failing, since the guard
genuinely doesn't reach the threshold for that real geometry), and its prm generator script
(`gen_guard_hopper_prm.py`). Removed from production: `BuildTimeloop.h`'s `#include` and sweep
registration, `SimulationParams.h`'s `closeSealedCellsMinNeighbours` field/parse.

## unused_contact_methods/feng2017_cylinder_cylinder_contact.h

`feng2017::fengCylinderCylinderContact` -- an alternative cylinder-cylinder contact method (Feng,
Han & Owen (2017), *Computer Methods in Applied Mechanics and Engineering* 315, 632-651), implemented
as an experiment to cross-check `paper_analytical::cylinderCylinderContact` (Kodam et al. 2010 /
Guo et al. 2012). Never wired into production (`ParticleKernels.h`) or any `DEMValidationScenarios.h`
test. [Claude, 2026-10-29] CORRECTION: initially reported as having "zero callers anywhere in this
showcase" -- that check was scoped only to `apps/showcases/Antidunes/` and missed
`tests/showcases_antidunes/PaperAnalyticalContactTest.cpp`, which DOES exercise it (6 call sites,
`testFeng*`/`testCylinderCylinderMissPatternAnalysis`) -- already commented out of that file's own
`main()` before this archival ("dropped from active development"), now removed from the file body too
(see `feng2017_validation_tests.cpp.excerpt` below). Moved out of `PaperAnalyticalContacts.h`
2026-10-29; still includes that file for its `CylinderState`/`ContactResult` dependencies, so it stays
compilable by hand if anyone wants to pick it back up.

## unused_contact_methods/feng2017_validation_tests.cpp.excerpt

The `testFeng*` + `testCylinderCylinderMissPatternAnalysis` functions excerpted verbatim from
`tests/showcases_antidunes/PaperAnalyticalContactTest.cpp` (lines 2462-2953 at archive time), which
validated `feng2017::fengCylinderCylinderContact` above. Reference only, not standalone-compilable --
depends on shared test helpers (`FuzzRng`, `GjkResult`, etc.) defined elsewhere in that file.

## unused_contact_methods/spherosimplex_box_contact.h

`boxBoxSpheroAnalytical` + the `eberly_box_distance` namespace (David Eberly's closed-form box-box
distance, ported from GeometricTools) -- the analytical alternative to the GJK/EPA-with-margin
spherosimplex box-contact path (`ModelParameters::useSpheroBoxContact`/`useSpheroBoxAnalytical`).
Unused by any thesis-cited result (confirmed: every thesis-cited `DropTest_250_random_*` prm sets
`useSpheroBoxContact false`; only superseded exploratory prms set it true). Removed entirely from the
live build 2026-10-29, including `ParticleKernels.h`'s dispatch and `SimulationParams.h`'s flags --
see the corresponding commit for the full removal footprint (constructor params, `DEMEngine`,
`CornerScenarioParams`/`Result`, 4 config-table rows). `DropTest_compare_sphero_gjk_gas.prm` and its
sbatch job moved alongside this file.

## unused_contact_methods/fraige_multipoint_box_contact.h

`boxBoxFraigeMultiPointContacts` + its private `clipSegmentToBoxLocal` helper -- a fuller, more
literal replica of Fraige et al. (2008)'s multi-point box-box contact method (adds edge-edge/
edge-face modes and a near-corner normal blend that the live, vertex-only `boxBoxMultiPointContacts`
does not have). Unused by any thesis-cited result. Removed entirely from the live build 2026-10-29,
including `ParticleKernels.h`'s dispatch (collapsed to always call `boxBoxMultiPointContacts`) and
`SimulationParams.h`'s flag. `DropTest_fraige_gas.prm`/`DropTest_fraige_250_gas.prm` and their sbatch
jobs moved alongside this file.

## arch_tunnel_containment/

`runArchTunnelContainment` (DEM-validation scenario `archTunnelContainment`) -- a fluid-free
geometric-containment check for `ComputeTunnelContacts` (`ParticleKernels.h`), the production code
path for the real-scale Pont Ferroviere bridge case study (`chapters/03_implementation/case-study/`,
`chapters/04_results/case-study/` -- intentionally excluded from the thesis per `CLAUDE.md`). Per
this showcase's own `FINAL_THESIS_DEM_VALIDATION.md` (line 120), this test already **FAILED 1/2**
and was explicitly noted there as "out of this plan's cylinder/cuboid-contact-detection focus... not
investigated further." Moved out of `DEMValidationScenarios.h`/`DEMValidation.cpp` 2026-10-29, along
with its 2 dedicated prms (`dem_validation_arch_tunnel.prm`, `dem_validation_arch_tunnel_box.prm` --
neither was referenced by any sbatch job).

**[Claude, 2026-11-16] Fuller cleanup completed at the author's request** -- the whole bridge
case-study code/config set flagged above as a separate decision is now archived too, completing
this entry. Grepped the entire repo root (not just this directory) for every symbol below before
removing anything.

Moved into this same directory:
- `ComputeTunnelContacts.h`: the full feature, self-contained -- `ArchTunnelParams`/`ArchSolidBox`/
  `archDeckFillBoxGeometry`/`archBlockParams`/`ArchHaunchSide`/`archHaunchBoxGeometry`/
  `parseArchTunnels` (all originally `SimulationParams.h`), `createArchTunnel` (originally
  `Utility.h`), and the `ComputeTunnelContacts` class itself (originally `ParticleKernels.h`,
  ~400 lines). Stays compilable by hand: includes `ParticleKernels.h` (for `FloaterHistoryLogger`
  and `closed_form_contacts::detail::farthestThetaOnCircleToLine`) and `Utility.h`.
- `DEMEngine_addArchTunnel.h.excerpt`: `DEMEngine::addArchTunnel` and its two private members
  (`DEMValidationInfra.h`) -- confirmed zero callers anywhere except `arch_tunnel_containment_test.h`
  (this directory), which now needs this excerpt re-added by hand to recompile. Note added to that
  test's own header documenting the dependency.
- `arch_tunnel_box_tests.cpp.excerpt`: `testArchTunnelBoxContainmentHandDerived`/
  `testArchTunnelBoxContainmentNoContact`, excerpted verbatim from
  `tests/showcases_antidunes/PaperAnalyticalContactTest.cpp` (lines 4256-4407 at archive time),
  their 2 call sites removed from that file's `main()`.
- `arch_geometry_tests.cpp.excerpt`: the 7 `testArch*Geometry*`/`testArchDeckFillBox*` functions,
  excerpted verbatim from `tests/showcases_antidunes/ContactDetectionTest.cpp` (lines 1416-1575 at
  archive time), their 7 call sites removed from that file's own `main()`.
- `prms/test_bridge_v1*.prm` (46 files) and `jobs/run_*test_bridge_v1*.sbatch` (48 files), alongside
  the 2 `dem_validation_arch_tunnel*.prm` already here from the first pass.
- `bridge_layout_optimizer.py`, `bridge_voxel_build.py`, `bridge_voxel_render_full.py`,
  `bridge_voxel_view.py`, and their generated outputs (`bridge_layout.html`, `bridge_layout.png`,
  `bridge_voxels_full.npz`, `bridge_voxels_full.png`, `bridge_voxels.png`). Compiled `__pycache__`
  for these scripts was deleted outright, not archived (regeneratable, not source).

Removed from production (verified via repo-wide grep before and after, both targeted symbol
searches and a full rebuild of `Antidunes`/`DEMValidation`/`DropTest`/`BedGeneration`, clean):
- `SimulationParams.h`: the `ArchTunnelParams` struct family and `parseArchTunnels` (moved whole,
  see above), the `archTunnels` field + its constructor init, and the `useArchTunnelAsMesaPDObject`
  field + its constructor init and doc comment.
- `ParticleKernels.h`: the `ComputeTunnelContacts` class (moved, see above) and every wiring site --
  the `computeTunnelContacts_` member, its lazy construction in the `ParticleKernels` constructor
  (gated on `useArchTunnelAsMesaPDObject`), `setHistoryLogger`/`setStepInfo` forwarding, and the
  `forEachParticle` invocation in the main DEM sub-cycle. `FloaterHistoryLogger::logTunnelContact`
  was deliberately left in place -- the archived `ComputeTunnelContacts` still calls back into it,
  same pattern as `feng2017_cylinder_cylinder_contact.h`'s continued use of `ClosedFormContacts.h`.
  Several comments elsewhere in this file compare a *different*, still-live algorithm to
  `ComputeTunnelContacts' own gate/bound/formula` as design precedent -- left untouched, still
  accurate (the comparison doesn't depend on where `ComputeTunnelContacts` physically lives).
- `Utility.h`: the `createArchTunnel` helper (moved, see above). `createContainerCylinder`'s own
  doc comment ("same flags as `createArchTunnel` below") was left as an accurate historical
  reference, not updated.
- `SetupStages.h` and `SetupStages_drop_test.h` (identically duplicated in both): `createArchSolidElements`/
  `createArchBlock` (both overloads) and the `useArchTunnelAsMesaPDObject`-gated construction loop
  in `setupParticles`.
- `DEMValidationInfra.h`: `DEMEngine::addArchTunnel` (excerpted, see above), its two private
  members, and its `if (computeTunnelContacts_) ...` call in `DEMEngine::step()`.
- `DEMValidation.cpp`: corrected the `archTunnelContainment` removal note (first pass) to stop
  claiming `ComputeTunnelContacts` and the `test_bridge_v1_*` family are still live -- they're
  archived now too.
- Two generic prm-template comment strings (`gen_kodam_round_prm.py`, `gen_kodam_packing_prm.py`:
  `"// No Pillars / CylinderPillars / ArchTunnels -- ..."`) were deliberately left alone -- they
  describe the *absence* of these obstruction types in a generated prm, not a dependency on the
  feature, and remain accurate either way.

Left untouched, confirmed NOT part of this feature despite mentioning it by name: two comments in
`tests/showcases_antidunes/ContactDetectionTest.cpp` comparing an unrelated, still-live cylinder
contact gate to "`ComputeTunnelContacts`' own X-extent gate" (still accurate); three comments in
`DEMValidationCornerAndMultiPoint.h` attributing a real, generic corner-contact bug to production
job 4196224, which happened to be a bridge "widerstacks" scenario (the test itself is a generic
box-corner regression check, unrelated to arch tunnels).

## (not archived — see `local_only/`)

[Claude, 2026-11-17] Repo cleanup at the author's request, ahead of the public GitHub mirror
(`github` remote, `dnyse/fslbm-dem-debris-mt`). A long tail of one-off debug/investigation scripts
(19 `analyze_*.py` plus `census_sealed_cells.py`, `visualize_spherobox.py`, and 8 similar one-offs),
23 `gen_*_prm.py`→unrelated raw CSV dumps (`wall_impact_*.csv`, 13 files), old chain-runner shell
scripts (`run_kodam_comparisons_chain.sh`, `run_repro_chain.sh`, `run_shape_comparison.sh`) and their
logs, the stale `changes.md` changelog (superseded by `THESIS_PLAN.md`), `MITIGATION_COMPARISON.md`,
`spherobox_visualization.png`, and the Pont-Ferroviere-only `gen_pont_ferroviere_debris.py` /
`pont_ferroviere_rue_goffontaine.md` were **not** moved here — they were `git rm --cached` and moved
to a sibling top-level `local_only/` folder instead, which is wholly `.gitignore`'d. Same reasoning
as the entries above (superseded/unused, kept for the author's own reference) but the author wants
this specific batch off the public mirror entirely, not just out of the way in a tracked `archive/`.
`THESIS_PLAN.md`, `DEM_VALIDATION_MATH_REFERENCE.md`, and `FINAL_THESIS_DEM_VALIDATION.md` (personal
working documents) got the same treatment, untracked in place and moved into `local_only/` too.

Also in `docs/`: `class-structure.md`, `forces.md`, `forces-diagram.{pdf,png,py}`, and
`material-phases-and-velocities.md` were deleted outright (not archived, not moved to
`local_only/`) at the author's request, after confirming none were referenced by current code or
the thesis; `docs/HOWTO_BUILD_AND_RUN.md` and `docs/clogging-animations/` were kept as still-current.

## unused_contact_methods/fraige_validation_tests.cpp.excerpt

`testBoxBoxFraigeVsGjkVsSatFuzz`/`testFraigeVertexFaceContact`/`testFraigeVertexNearCornerBlend`/
`testFraigeEdgeFaceContact`/`testFraigeEdgeEdgeContact`, excerpted verbatim from
`tests/showcases_antidunes/PaperAnalyticalContactTest.cpp` (lines 4195-4532 at archive time). Unlike
the feng2017 tests, these WERE actively called from that file's own `main()` before this archival;
those 5 calls were removed as part of the same change. Reference only, not standalone-compilable.

## unused_contact_methods/: box-contact mitigations, Tier-1 variants, analytic box buoyancy (2026-11-20)

Eight experimental `ModelParameters` switches that defaulted to off and were never turned on in a
thesis-cited run, removed from the live build together with their companion values:

- `useRoundedBoxCorners`/`boxRoundingRadius`, `useBoxOverlapCap`/`boxOverlapCapFraction`,
  `useAxisHysteresis`/`axisHysteresisMargin`: variants of `closed_form_contacts::boxBoxSAT`
  (rounded edge-edge normal, overlap cap, sticky axis choice) against the SAT axis-flip discontinuity.
- `useMultiPointBoxContact` and `useContactManifold`/`manifoldFadeRate`: one contact point per embedded
  vertex (`BoxVertexContact`, `boxBoxMultiPointContacts`, now in `multipoint_box_contact.h`) and a
  persistent faded vertex set on top of it.
- `useTightTier1`, `skipTier1`: capsule-bound and skip variants of the Tier-1 check in
  `cylinderCylinderResultIsUncertain`. The capsule bound itself stays live as cylinder Layer 1
  (`useTightLayer1`); its conservativeness test is now `testCapsuleBoundConservative10k` in
  `PaperAnalyticalContactTest.cpp` (the original test is in `tight_tier1_accuracy_test.cpp.excerpt`).
- `useAnalyticBoxBuoyancy`: used the analytic OBB/disk clip for the submerged volume of floating boxes and
  cylinders instead of the 12³ sampling. Sampling is now the only method; the analytic clip is still
  computed and logged as a cross-check.

With all switches at their defaults the live code paths are unchanged (verified bit-identical on the unit
tests and the tracked `dem_validation_*.prm` runs). The removed code paths and parameter docs are in
`box_contact_mitigations.cpp.excerpt`. The DEM-validation scenarios that compared these mitigations
(`DEMValidationCornerAndMultiPoint.h`: cornerNormalDiscontinuity, cornerFixComparison,
multiPointDeepOverlapStability/Fuzz, multiPointEdgeCases, contactStiffnessCalibrationSweep,
fallOntoPileStability, multiPointBoxContactFuzz; and `DEMValidationProductionCollisionReplay.h`) moved
here as well; they use the pre-archive `DEMEngine` signature and would need adapting to compile again.
Their prms were never git-tracked.
