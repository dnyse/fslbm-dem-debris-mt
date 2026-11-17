# Antidunes — Reproducing the Thesis Results

This showcase extends waLBerla's Antidunes antidune-bedform code (see `CLAUDE.md` and
`docs/HOWTO_BUILD_AND_RUN.md` for the original baseline scenario) with floating debris
(spheres/cylinders/cuboids), non-spherical contact detection, and a clogging measurement, for
*"Numerical Simulation of Floating Debris Transport and Clogging Using Coupled Free-Surface Lattice
Boltzmann and Discrete Element Methods"*. This file is about reproducing **that thesis's** specific
results. It assumes you've read the thesis; it does not re-explain the physics.

## Reproducibility caveat — read this first

`prms/*.prm` is **gitignored by design**. Hundreds of exploratory variants accumulated there during
development and are not a version-controlled record of what actually ran. A few of the files pointed
to below were themselves *reconstructed* after the fact from a cluster job's own captured log
(`.fulllog.gz`), because the original `.prm` no longer existed on disk — see e.g. the header comment
of `prms/bogner3_r025_cyl_s32_eqfull_vtk_strict_noliqgas_unpinned_gas.prm`. So "byte-identical to the
job that produced the published number" is not guaranteed everywhere.

Where a result below has an external provenance record confirming exactly which job IDs produced it,
this file says so. Where it doesn't, it says **"unconfirmed"** rather than presenting a best guess as
settled fact — check with the author, or look in
`/run/media/dhuber/data2/debris-clogging-archive/` (one `.tar.zst`/folder per run family; some have a
`BACKUP_*.md` naming the exact job IDs and which thesis figure they feed) before trusting it blindly.

Also: `THESIS_PLAN.md` (1 MB, the detailed development log) stops being updated in early October; a
good number of the *final* configs referenced below are dated into November and are not narrated
there. File naming and dates, cross-checked against the archive where possible, are the best
available evidence for what's "final" versus "an earlier abandoned variant."

## Build

From the waLBerla repo root (`/home/dhuber/University/Thesis/mt-dennys`):

```bash
python -m venv .venv && source .venv/bin/activate
pip install "pystencils~=1.3" "lbmpy~=1.3" jinja2

cmake -G Ninja -DWALBERLA_BUILD_WITH_CODEGEN=ON -DWALBERLA_BUILD_SHOWCASES=ON \
  -DWALBERLA_BUILD_WITH_MPI=ON -DPython_ROOT_DIR="$(pwd)/.venv/bin" -B build/ .

cmake --build build/ --target Antidunes BedGeneration DEMValidation DropTest
```

All four targets need `WALBERLA_BUILD_WITH_CODEGEN=ON` except `BedGeneration` (confirmed in
`CMakeLists.txt`). `./antidunes.sh --configure --build` only builds `Antidunes` + `BedGeneration` —
`DEMValidation` and `DropTest` (used below for RQ2) are **not** wired into `antidunes.sh` and need the
plain `cmake --build` above.

On the Fritz cluster (NHR@FAU): `./antidunes.sh --sync` (rsyncs this repo to `fritz.iwia:~/mt-dennys/`),
then `./antidunes.sh --configure --build --jobs 32` there. See `antidunes.sh --help` for every flag.

---

## 1. Floating equilibrium (RQ1)

All `Antidunes` runs (full LBM+DEM+MEM). Config naming: `prms/bogner3_r<RRR>_<shape>_...`, where
`<RRR>` = density ratio ×1000 (`r025`→0.25, `r050`→0.5, `r075`→0.75) and `<shape>` is `sphere`, `cyl`
(cylinder), or **`deep`, which means cuboid** (confirmed by reading the file contents — not an obvious
name).

```bash
mpirun -n <N> build/apps/showcases/Antidunes/Antidunes prms/<file>.prm
```
Check the file's own `BlockForestParameters` (`cellsPerBlock` × block counts) for the exact `<N>` it
needs before running — these equilibrium configs are small (`cellsPerBlock <32,15,15>`-class), not
the 144/288-rank clogging scale below.

### Sphere & cylinder (Fig. 5.1, "within 0.5 LU" headline)

**Confirmed** by `debris-clogging-archive/BACKUP_20261102_bogner3_strict_classification.md`: jobs
4408736–4408741 (the strict-Bogner-classification fix) are what's actually plotted.

| ρ | Cylinder | Sphere |
|---|---|---|
| 0.25 | `prms/bogner3_r025_cyl_s32_eqfull_vtk_strict_noliqgas_unpinned_gas.prm` | *(no local file survives — see below)* |
| 0.50 | `prms/bogner3_r050_cyl_s32_eqfull_vtk_strict_noliqgas_unpinned_gas.prm` | *(same)* |
| 0.75 | `prms/bogner3_r075_cyl_s32_eqfull_vtk_strict_noliqgas_unpinned_gas.prm` | *(same)* |

**Unconfirmed / action needed:** no `..._sphere_..._strict_noliqgas...` file exists in `prms/` for any
density ratio — only the archived job output does
(`bogner3-eqfull-strict-noliqgas_20261102.tar.zst`, folders
`final_anim/440873{7,9,1}_b3_r0{25,50,75}_sphere_s32_eqfull_vtk_strict_noliqgas/`). To rerun the
sphere side, either pull a working `.prm` out of that archive's own job metadata, or reconstruct one
from a job's `.fulllog.gz` the way the cylinder file's header comment describes doing.

### Cuboid heel angle / heave (Fig. 5.3–5.4, "within 0.6°")

The file's own header comment self-identifies as *"the thesis's own primary cuboid equilibrium
run ... used for fig:eq-cuboid-evolution"* — fully unpinned (6 dof), no surface tension:

```bash
mpirun -n <N> build/apps/showcases/Antidunes/Antidunes prms/bogner3_r025_deep_s32_eqfull_vtk_sigma0_water_half_200k_noliqgas_unpinned_gas.prm
# + r050 / r075 siblings, same suffix
```

**Caveat:** `BACKUP_20261102_bogner3_strict_classification.md` records an *earlier* cuboid family
(jobs 4351201–4351203) as feeding "Figure 5.3" as of 2026-11-02. The `unpinned` files above are dated
2026-11-06 and their own header says they supersede that prior (pinned) run per author request — so
they're very likely current, but that supersession has no external BACKUP/PROVENANCE confirmation the
way the sphere/cylinder one does.

Appendix C's ~40× more viscous rerun (Bogner's own liquid, `tau=0.526`): **unconfirmed** — search
`prms/` for the same `deep`/density-ratio family without `sigma0_water` in the name (e.g. a
`..._tau0526_...` variant); not independently traced in this pass.

### MEM-only vs. analytical submerged-rise comparison

**Unconfirmed** — not traced to a specific file in this pass. Candidates by name: the
`gen_buoyancy_bottom_prms.py`-generated `BedGeneration_buoyancy_r{1,2,4}.prm` /
`run_72_bedgen_buoyancy_r{1,2,4}*.sbatch` family, or the `cuboid-buoyancy-sweep-combined.tar.zst` /
`cylinder-buoyancy-sweep-combined.tar.zst` / `sphere-buoyance-sweep-combined.tar.zst` archives. Check
with the author before relying on either guess.

---

## 2. Contact-detection validation (RQ2)

### Wall-impact drop test (cylinder & cuboid, 0°–90° in 5° steps)

Entry point: the standalone `DEMValidation` executable (fluid-free, DEM-only — see its own file
header). The `[DEMValidation]` block's `testCase` key selects the scenario (full list in
`DEMValidation.cpp`'s `main()`); for this result it's `wallImpactCylinder` / `wallImpactBox`.

```bash
build/apps/showcases/Antidunes/DEMValidation prms/dem_validation_wall_cylinder_kodam_sweep_translational.prm
build/apps/showcases/Antidunes/DEMValidation prms/dem_validation_wall_cylinder_kodam_sweep_rotational.prm
build/apps/showcases/Antidunes/DEMValidation prms/dem_validation_wall_box_sweep_translational.prm
build/apps/showcases/Antidunes/DEMValidation prms/dem_validation_wall_box_sweep_rotational.prm
```

`_translational` vs. `_rotational` select the effective-mass formula — these produce the thesis's
34/59 & 31/57 (translational) vs. 59/59 & 57/57 (rotational) pass counts. Single process, runs in
seconds; no SLURM job needed. A `_sweep_fine` variant also exists for both shapes (finer tilt
resolution); its exact role in the thesis is **unconfirmed**.

### Drop-test packing vs. Kodam et al.'s experiment

Entry point: the standalone `DropTest` executable (fluid-free DEM at real block-forest scale — see
its file header). One `.prm` per container, no `testCase` key.

Four containers (the "250 floaters: 4 containers" demo):
```bash
mpirun -n 36 build/apps/showcases/Antidunes/DropTest prms/DropTest_250_random_cylinders_gas.prm
mpirun -n 36 build/apps/showcases/Antidunes/DropTest prms/DropTest_250_random_cuboids_gas.prm
mpirun -n 36 build/apps/showcases/Antidunes/DropTest prms/DropTest_250_random_mixed_gas.prm
```
**Confirmed**: `jobs/run_36_droptest_250_three_shapes.sbatch` runs exactly this, and
`debris-clogging-archive/antidunes-droptest-250-random-three-shapes_20261003.tar.zst` is the archived
output (jobs 4284397/4284398/4287115/4292501), feeding thesis Figure 5.12.

The round Kodam container (fill height / solid fraction vs. the published experiment):
```bash
mpirun -n 9 build/apps/showcases/Antidunes/DropTest prms/DropTest_kodam_round_axial_gas.prm
```
**Unconfirmed which exact file produced the published 53.3/55.8/57.8 mm numbers.** The repo also has
a documented chain script, `run_kodam_comparisons_chain.sh` (repo root), but it currently points at
the *older* `DropTest_kodam_packing_gas.prm`/`_gjk_gas.prm`/`_pureanalytical_gas.prm`, not the
`_round_`/`_round_axial` files. `DropTest_kodam_round_axial_gas.prm` (2026-10-10, the most recent of
the `kodam_round*` family) is the best guess by name and recency — not independently verified against
an archive/provenance record. Check before citing.

### Hand-derived cases & the 100,000-pair fuzz test

These are `ctest` unit tests, not `.prm`-driven runs — confirmed in
`tests/showcases_antidunes/CMakeLists.txt`:

```bash
ctest --test-dir build -R ANTIDUNES_ContactDetection --output-on-failure        # hand-derived fixtures
ctest --test-dir build -R ANTIDUNES_PaperAnalyticalContacts --output-on-failure # the fuzz test
ctest --test-dir build -R ANTIDUNES_ClipAnalyticGeometry --output-on-failure    # geometry clipping
```
`ANTIDUNES_PaperAnalyticalContacts` is confirmed by its own CMake comment to be "compared ... against
SolidCylinderContactDetection's existing GJK/EPA-fallback results" over the same random pairs used for
the thesis's accuracy validation; set `ANTIDUNES_REFERENCE_DUMP=<dir>` to additionally dump the random
pairs. `./antidunes.sh --tests` builds and runs the *entire* `ANTIDUNES_*` ctest battery (which
includes many unrelated development-diagnostic tests, not just these three) in one go.

---

## 3. Clogging (RQ3)

All three layouts are `Antidunes` runs. Every `jobs/*.sbatch` script below follows the same pattern:
`antidunes.sh --build` then `antidunes.sh --run --np <N> --prm prms/<file>.prm`, then archives
`vtk-out/`/`eval/`/`snapshots/` to `$HPCVAULT` (falling back to a local directory if unavailable) —
open any one of them to see the exact commands. **144 ranks / 2 nodes** for centred pillar and hopper,
**288 ranks / 4 nodes** for two pillars (matches the thesis's own performance table).

### Centred pillar (100 spheres, deterministic blockage)
```bash
sbatch jobs/run_144_clog_centred_pillar_startblock_cappid.sbatch
```
Config: `prms/clog_centred_pillar_startblock_cappid_gas.prm`. **Unconfirmed:** no dedicated
clean-reference `.prm`/`.sbatch` exists under a "centred + clean/reference" name; per
`THESIS_PLAN.md`'s 2026-10-02 entry the clean reference reuses earlier `clog_measurement_validation_*`
jobs (4245061/4260083) from before this scenario was renamed to `startblock` — not re-verified here.

### Hopper gap (two wall pillars, 16 LU gap)

| Scenario | Config (`prms/`) | Job (`jobs/`) |
|---|---|---|
| Clean reference | `hopper_clean_reference_150k_gas.prm` | *(no matching sbatch found — unconfirmed)* |
| 406 rods (headline) | `hopper_clog_cyl_startblock_skin0p5_cappid_gas.prm` | `run_144_hopper_clog_cyl_startblock_skin0p5_cappid.sbatch` |
| 406 rods, repeat (scatter check) | `hopper_clog_cyl_startblock_skin0p5_cappid_rep2_gas.prm` | `run_144_hopper_clog_cyl_startblock_skin0p5_cappid_rep2.sbatch` |
| 406 cuboids | `hopper_clog_box_startblock_skin0p5_cappid_gas.prm` | `run_144_hopper_clog_box_startblock_skin0p5_cappid.sbatch` |
| 1015 spheres, D=6 | `hopper_clog_sph_startblock_skin0p5_cappid_gas.prm` | `run_144_hopper_clog_sph_startblock_skin0p5_cappid.sbatch` |
| 600 spheres, D=8 | `hopper_clog_sph4_startblock_cappid_gas.prm` | `run_144_hopper_clog_sph4_startblock_cappid.sbatch` |
| Mixture (140 rods/140 cuboids/315 spheres) | `hopper_clog_mixed_startblock_skin0p5_cappid_gas.prm` | `run_144_hopper_clog_mixed_startblock_skin0p5_cappid.sbatch` |

**Contact-skin sensitivity** (same 406-rod scenario, varying `cylinderContactSkin`): `..._noskin_cappid`,
`..._skin0p05_cappid`, `..._skin0p1_cappid`, `..._skin0p25_cappid`, plus the `skin0p5` headline above —
each with a matching `jobs/run_144_*.sbatch`.

**Other sensitivity runs**, same naming pattern: drive-by-mass (`..._drivebymass_gas.prm`), surface
tension ÷10 / ÷100 (`..._sigmadiv10_gas.prm`, `..._sigmadiv100_gas.prm`).

### Two pillars (812 rods, wide-gap passage)
```bash
sbatch jobs/run_288_twopillar_clog_cyl812_startblock_cappid.sbatch   # 812 rods
sbatch jobs/run_288_twopillar_clean_reference_150k.sbatch            # clean reference
```
Config: `prms/twopillar_clog_cyl812_startblock_cappid_gas.prm`.

### Local smoke test

No reduced-resolution/short-step variant of any clogging scenario is documented anywhere in this
repo. The `jobs/*.sbatch` pattern can be invoked directly without `sbatch` for a quick check, but at
the full process/step count as written. A genuine smoke test means copying one `.prm`, lowering
`timesteps` and the domain scale by hand, and rechecking the process-count constraint described in
`docs/HOWTO_BUILD_AND_RUN.md` ("Domain / Process Count Constraints").

### Post-processing

Figures/tables are produced by scripts in the **thesis repo**, not here:
`overleaf/scripts/clogging/*.py` and `overleaf/scripts/figures/*.py`, reading each run's `eval/`
SQLite output and/or extracted `.dat` files. Cross-reference by run name or job ID over there.

---

## Why this file exists — the state of the rest of the repo

A large amount of exploratory scratch accumulated here over the course of the thesis work: hundreds
of one-off `prms/*.prm` variants (gitignored, not meant to be read as history) and dozens of
`scripts/gen_*.py` prm generators, plus long dated comment trails inside the `.h`/`.cpp` files that
record investigation history rather than current behaviour (one-off `analyze_*.py` debugging scripts
and similarly superseded notes were removed from the tracked tree in a 2026-11-17 cleanup — see
`archive/ARCHIVE_README.md`). This README only points at the configs confirmed (or best-inferred,
explicitly flagged as such) to produce the thesis's published numbers. Most of the rest of
`prms/`/`jobs/`/`scripts/` are superseded dead ends kept only because they're gitignored/untracked
(or, for `scripts/`, simply unused by any current run) rather than cleaned up — not a reliable source
of "what the thesis actually used" on their own.
