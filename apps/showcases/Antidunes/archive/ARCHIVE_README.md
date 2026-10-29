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
Guo et al. 2012). Confirmed zero callers anywhere in this showcase outside its own namespace -- never
wired into production (`ParticleKernels.h`) or any `DEMValidationScenarios.h` test. Moved out of
`PaperAnalyticalContacts.h` 2026-10-29; still includes that file for its `CylinderState`/`ContactResult`
dependencies, so it stays compilable by hand if anyone wants to pick it back up.
