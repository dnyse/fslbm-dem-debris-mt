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
was removed; `tsunami.cpp`/`TsunamiDataInitialization.h` are not `#include`d by anything else in this
directory, so nothing else needed updating.

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

## unused_contact_methods/fraige_validation_tests.cpp.excerpt

`testBoxBoxFraigeVsGjkVsSatFuzz`/`testFraigeVertexFaceContact`/`testFraigeVertexNearCornerBlend`/
`testFraigeEdgeFaceContact`/`testFraigeEdgeEdgeContact`, excerpted verbatim from
`tests/showcases_antidunes/PaperAnalyticalContactTest.cpp` (lines 4195-4532 at archive time). Unlike
the feng2017 tests, these WERE actively called from that file's own `main()` before this archival;
those 5 calls were removed as part of the same change. Reference only, not standalone-compilable.
