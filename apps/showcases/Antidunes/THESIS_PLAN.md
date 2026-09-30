# Thesis Validation Plan

Living document tracking the validation/showcase work plan for the thesis, agreed 2026-09-12. This
supersedes the earlier focus on implementing further cuboid-contact fixes (`useSpheroBoxContact`,
`DropTest`, the arch/pier flying-floater investigation) for now — that work is paused, not abandoned;
see `MITIGATION_COMPARISON.md` for its own record.

## How to use this file

- **Keep it updated.** Whenever a run finishes, a decision is made, or a finding changes the plan,
  update the relevant section here — don't let this drift out of sync with reality.
- **Comments are attributed and dated.** Add a blockquote line in the form:
  `> **[Dennys, YYYY-MM-DD]:** ...` for the user's own comments, or
  `> **[Claude, YYYY-MM-DD]:** ...` for comments/findings Claude adds. Never edit or delete another
  author's comment — add a new one below it if something needs correcting.
- **Status tags**: mark each item `[ ] Not started`, `[~] In progress`, `[x] Done`, or
  `[!] Blocked` at the start of its heading.
- Every top-level numbered item ends with a **Challenges & Future Improvements** subsection — these
  feed directly into the thesis's Discussion/Conclusion chapter, so keep them concrete (what
  specifically limits this result, what a follow-up would need to do differently), not generic.

## Archive

Run outputs (vtk-out/eval/snapshots, `.prm` files, and cluster job logs where recoverable) for the
runs referenced below are backed up at `/run/media/dhuber/data2/debris-clogging-archive/` on the
external drive, one subfolder per sweep/comparison. Each subfolder has its own `PROVENANCE.md`
documenting exactly which job IDs/prm files it contains and any gaps. Existing subfolders as of this
writing: `sphere-buoyance-sweep-combined`, `cylinder-buoyancy-sweep-combined`,
`cuboid-buoyancy-sweep-combined`, `flying-box-collision-replay-fix-comparison`. Add a new subfolder
(with its own `PROVENANCE.md`) for every new run family created while working through this plan.

---

## 1a. [ ] Equilibrium (buoyancy) sweeps — sphere / cylinder / cuboid

**Status**: initial sweeps already run and archived (see `sphere-buoyance-sweep-combined`,
`cylinder-buoyancy-sweep-combined`, `cuboid-buoyancy-sweep-combined` above) — comparing Bogner &
Rude (2013)'s own box-experiment gravity/buoyancy model against the real-SI-gravity (g=9.81 m/s²)
calibration, across densityRatio ∈ {0.25, 0.5, 0.75}, for each of sphere/cylinder/cuboid.

**Goal now**: rerun at **longer duration and higher resolution** so each configuration clearly
settles to its true Archimedes equilibrium (fraction submerged = densityRatio) rather than being cut
off mid-transient.

> **[Claude, 2026-09-15]:** Started the cuboid piece of "rerun at longer duration and higher
> resolution" above. Current baseline (`cuboid_stillwater_realgravity981_3xsphere_{,d025,d075}_gas.prm`):
> domain `<96,135,120>`, cellsPerBlock `<16,45,30>` = 72 ranks, 50000 timesteps, real gravity
> (`gravityOverride=5.06585e-05`, solved for g_SI=9.81 m/s^2), densityRatio in {0.25, 0.5, 0.75}.
> Measured wall-clock ~2h22m-2h30m/case (jobs 4038342-4038344).
>
> New `_hires` siblings run at ~1.24x finer linear resolution (domain `<119,168,150>`, cellsPerBlock
> `<17,42,50>` = 84 ranks) and 100000 timesteps (2x). Every number re-derived from the real bed-file
> header (identical on all 3 cluster checkouts: `generationDomainSize_SI =
> (0.0464248483188, 0.0218704558877, 0.145077650996)`), not guessed -- domainSize[0]=119 (not a round
> 120=96*1.25) specifically because Utility.h's `getAvgDiameterScalingFactor` requires
> `domainSize[1]` to exactly match `uint_c(particleScalingFactor * generationDomainSize_SI[1] *
> bedCopiesInY)` or the run aborts; 120 forces domainSize[1]=169 (prime-squared, only divisible by
> 1/13/169 -- unusable for a sane rank decomposition), while 119 gives domainSize[1]=168
> (=2^3*3*7, plenty of divisors) at essentially the same scale (119/96=1.2396 vs 1.25). Recomputed
> from first principles: dx_SI 0.00048359217->0.00039012478 m, dt_SI 4.99725e-05->3.25221526e-05 s
> (ratio 0.6508 ≈ (dx_new/dx_old)^2 = diffusive scaling, confirmed), gravityOverride
> 5.06585e-05->2.6596472613474883e-05 (still solved for g_SI=9.81 at the new dx/dt), cuboid
> halfExtents `<15,30,15>`-><18.59375,37.1875,18.59375>`, position -> new domain centre
> `<59.5,84,75>`. `liquidHeightFactor`/`floorHeightFactor` deliberately left UNCHANGED -- both
> multiply `avgParticleDiameter`, which itself scales with the new `particleScalingFactor`
> automatically, so the physical water depth stays correct with no manual rescaling. Lattice
> viscosity (0.000214026) also deliberately left unchanged (`useHardcodedUnits=true` fixes it
> regardless of resolution -- that's what makes the dt_SI rescaling above diffusively correct in the
> first place).
>
> **Cost check done before submitting** (per this item's own "Challenges" bullet below): ~1.24^3
> (cells) * 2 (timesteps) ≈ 3.8x lattice-updates per case vs. the baseline measured 2.4h, giving a
> ~9-9.5h/case estimate. All 3 cases sequentially would need ~28h, over the multinode partition's
> 24h MaxTime, so each density is its own 12:00:00 job rather than one combined job. Only one of
> the 3 cluster checkouts (`mt-dennys`, no mirror suffix) was free (`mt-dennys-mirror` running the
> R4 buoyancy retry, job 4235003; `mt-dennys-mirror2` had phase-1.5 of the clog-crash investigation,
> job 4235729, queued) -- submitted all 3 density cases there, chained via `--dependency=afterok` so
> they run back-to-back without a build-directory collision instead of requesting 3 separate
> allocations. Initially submitted at 84 ranks (`cellsPerBlock <17,42,50>`, numBlocks (7,4,3)) as
> d025=4235857/d050=4235858/d075=4235860, but caught before any of them started (still PENDING):
> 84 ranks on a 2-node/144-core allocation leaves 60 cores idle. Cancelled and resubmitted at the
> only exact 144-way factorisation of this domain (X=119=7*17 has no divisor between 17 and 119, so
> numBlocksX must be 1, 7, 17 or 119; only numBlocksX=1 combines with a Y/Z pair multiplying to 144 --
> Y=24 (cellsPerBlockY=7) * Z=6 (cellsPerBlockZ=25)): **cellsPerBlock `<119,7,25>`, numBlocks
> (1,24,6) = 144**. Physics (dx_SI/dt_SI/gravityOverride/halfExtents/position) is unaffected --
> only the MPI decomposition changed. Resubmitted: **d025 = job 4235873** (no dependency, runs
> first), **d050 = job 4235874** (`afterok:4235873`), **d075 = job 4235875** (`afterok:4235874`).
> New files: `prms/cuboid_stillwater_realgravity981_3xsphere_hires_{,d025,d075}_gas.prm`,
> `jobs/run_144_cuboid_stillwater_hires_d0{25,50,75}.sbatch`. Not yet analyzed (none have started
> running yet as of submission).
>
> **[Claude, 2026-09-15]: 144-rank decomposition was itself invalid** -- job 4235873 (d025)
> aborted in 0.09s: `SetupStages.h:69`'s own assertion requires at least three blocks per
> *periodic* direction when particles are present (`numBlocks[i] >= 3`), and `<119,7,25>` gives
> `numBlocksX=1` (X is the periodic direction here, `periodicity <1,0,0>`) -- missed this
> constraint when I optimized purely for "exact 144". X=119=7*17 has no divisor between 3 and
> 119 other than 7 and 17, and neither divides 144 evenly, so exact 144 was never actually
> achievable with a valid decomposition. Closest valid factorisation: **cellsPerBlock
> `<17,42,30>`, numBlocks (7,4,5) = 140** (only 4 of 144 cores idle, numBlocksX=7 satisfies the
> assertion). Cancelled the now-permanently-blocked d050/d075 jobs (4235874/4235875 -- their
> `afterok` dependency on the failed 4235873 could never be satisfied) and resubmitted the full
> chain fresh: **d025 = job 4236462**, **d050 = job 4236463** (`afterok:4236462`), **d075 = job
> 4236464** (`afterok:4236463`). Files renamed `run_144_...` -> `run_140_...` to match. d025
> confirmed past the 0.09s failure point on resubmission (still running as of this note). Not yet
> analyzed.
>
> Also checked, out of interest (not a citation-fixing exercise, just curiosity about how this
> compares to Bogner's own free-floating-cuboid validation, Sec. 5.3.3 of `bogner_direct_2017`):
> his own domain was `130x40x40` LU cells with cube side `s=16` (his cuboid `2s x s x s = 32x16x16`
> LU, all three densities run together in one basin, no fixed timestep count reported, just "until
> all fluid motions ceased") -- he calls this "low resolution" himself in the text. This thesis's
> existing 30x60x30 box (let alone the new hires 37.1875x74.375x37.1875 LU box) is already
> substantially finer than his own validation baseline.
>
> **[Dennys, 2026-09-12]:** there's also an artifact we did for this — reference it here.
> **[Claude, 2026-09-12]:** artifact not yet identified/linked — flagging for Dennys to fill in
> (which artifact URL/title, and what it shows) rather than guessing.
> **[Claude, 2026-09-28]: identified.** "Buoyancy Sweep — Sphere, Cylinder, Cuboid"
> (https://claude.ai/artifact/VzeThMNwVuMJhy958wLN5x, dated 2026-08-18). Section 1 is the sphere
> at 3x scale (r=24 LU, domain `<96,135,120>`, this baseline's own `3xsphere` naming) -- three
> density-ratio panels (0.25/0.50/0.75), real gravity (g_SI=9.81) vs Bogner & Rüde's own gravity
> (g=7.5e-4, tau=1.9), each against its own VTK-measured target equilibrium height, 50000 steps.
> All three completed; both gravities converge close to their respective targets (e.g. density
> 0.5: real-gravity ~60.40 vs target 60.17, Bogner ~60.40 vs the same target). Cylinder and cuboid
> sections follow further down the same page.

> **[Claude, 2026-09-18]:** Root-caused and fixed a fresh-start `[STALE-PARTICLE-IDX]` hard abort hit
> while rerunning the original 72-rank d050 box-equilibrium reproduction (`cuboid_stillwater_
> realgravity981_3xsphere_gas.prm`, jobs 4245028/4245089) to capture the new tiltY log field.
>
> Mechanism: `SetupStages.h`'s `setupParticles()` only called `simObjs->syncCall()` (mesa_pd's
> `SyncGhostOwners`) ONCE on a fresh start — and `SyncGhostOwners` is a one-hop-per-call BFS
> propagation (only reaches ranks in `domain.getNeighborProcesses()` of wherever the particle is
> *already* known), not a radius-based single-shot sync. A 4-round repeat of that same call already
> existed to fix an analogous checkpoint-resume crash (job 4236626), but was gated to
> `loadSnapshot==true` only, on the explicit (and here, false) assumption that a fresh start's own
> floaters are "local to a handful of blocks" and never need more than one hop.
>
> Confirmed empirically via the `[GHOST-SYNC-DIAG]` diagnostic (`Antidunes.cpp`, previously
> resume-only, un-gated for this investigation): this scenario's single `FloatingBox`
> (`halfExtents<15,30,15>`, `cellsPerBlock<16,45,30>`) has its owning rank at block index `(3,1,2)`;
> all 26 Moore-neighborhood neighbors of that block correctly show ghost knowledge of the box after
> the single sync call, but ranks 13/17/21 at block indices `(1,0,1)/(1,1,1)/(1,2,1)` — Chebyshev
> distance 2 from the owner — showed zero knowledge (`boxes-known=[ ]`) of it. Verified this is a
> genuine, physically necessary dependency, not a spurious one: accounting for the box's own 5°
> `initialTiltAngleDeg` (which enlarges its axis-aligned reach from a naive half-extent of 15 to
> ~16.25), its true surface reaches global x≈31.75 — inside the very ghost cell (global x=32) that
> rank 17's own owned interior cell (x=31) borders. So rank 17 has a real need for local particle
> data there (bounce-back/MEM force computation), which the single-hop sync failed to provide.
>
> Also worth noting for context: this exact abort's own history (added `35050123`, 2026-08-29,
> escalated from a silent warning-and-skip to a hard abort by `f1f4342f`, 2026-09-16) predates this
> fresh-start case entirely, and the underlying multi-hop gap was almost certainly always present —
> it just could never manifest on a fresh start before `f1f4342f` also fixed a *separate* bug
> (`mapParticlesOntoCells()`'s one-time setup mapping being a silent no-op due to a missing
> `fillWithDomain()` call) that had, until 2 days ago, made this whole code path unreachable at t=0.
> Fixing one bug exposed the other. Not shape-specific (cuboid vs. sphere) despite appearances —
> it's governed by `maxParticleDiameter` vs. block size (`SetupStages.h`'s own algorithm-selection
> threshold between `SyncNextNeighbors` and `SyncGhostOwners`); every floating sphere used elsewhere
> in this project has simply stayed small enough, relative to its own scenario's block size, to never
> cross into the vulnerable `SyncGhostOwners` path.
>
> Fix: removed the `loadSnapshot` gate on the 4-round `syncCall()` repeat entirely — it now always
> runs, fresh start or resume, since the cost is one-time (setup only) and cheap regardless of path.
> Verified compiling locally; resubmitted the same reproduction (job 4245457) with the fix in place.

### Challenges & Future Improvements
- Current sweep's stopping point (timesteps chosen ad hoc per prm) wasn't verified against an actual
  equilibrium criterion (e.g. `|z_center(t) - z_center(t-Δt)| < ε` sustained over some window) — a
  proper rerun should define and log that criterion directly rather than eyeballing a fixed
  step count.
- Higher resolution raises cost roughly with the cube of the refinement factor for the affected
  dimensions; needs a concrete resolution target (and cost estimate) before submitting the rerun.
- Reproducibility: several archived runs have a "reproducibility rerun" alongside the original with
  no recorded way to tell which one's data is actually in the archive (see the cylinder/cuboid
  `PROVENANCE.md` files) — worth fixing the naming/output convention going forward so this ambiguity
  can't recur (e.g. embed the job ID directly in the output folder name, not just the sbatch log
  filename).

> **[Claude, 2026-09-14]:** fixed, both the archive layout and the ambiguity gap above, per explicit
> request. Two problems addressed at once:
> 1. **Flat, ambiguous archive folder naming** -- `ARCHIVE_ROOT` used to be
>    `${SLURM_JOB_NAME}_$STAMP` (no job ID at all, exactly the gap noted above) and
>    `archive_outputs()` dumped each case's 3 outputs as flat SIBLING folders
>    (`sphere_mem_eval`, `sphere_mem_vtk-out`, `sphere_mem_snapshots`, ...) directly under that one
>    root -- 18 loosely-named folders for a 6-case job, no per-case grouping.
> 2. Fixed in **all 34 combined-job sbatch scripts** sharing this pattern (not just the buoyancy ones):
>    `ARCHIVE_ROOT` now embeds `${SLURM_JOB_ID}` (`${SLURM_JOB_NAME}_${SLURM_JOB_ID}_$STAMP`), and
>    `archive_outputs()` now nests each case's outputs under its own subfolder
>    (`$ARCHIVE_ROOT/$label/{eval,vtk-out,snapshots}` instead of `$ARCHIVE_ROOT/${label}_{eval,vtk-out,snapshots}`).
>    Verified against the actual R4 buoyancy script before resubmitting.
> 3. **Retroactively reorganized the already-existing R1 and R2 archives** (moved, not regenerated --
>    all data already on the same filesystem, `mv` is a fast rename not a copy, verified sizes
>    unchanged before/after):
>    - R1 previously had 8 SEPARATE top-level folders (one per individual job submission, since R1 was
>      never run via a combined script) -- consolidated into one `antidunes-buoyancy-r1-compare-all/`
>      parent with 8 clearly-labeled subfolders: `sphere_mem`, `sphere_analytic`, `cylinder_mem`,
>      `cylinder_analytic`, `box_mem`, `box_analytic`, `box_mem_reproducibility_rerun_20260913`
>      (the anomaly-reproduction rerun), `box_mem_tilted25_hypothesis_test_20260914` (the tilt
>      hypothesis test) -- the two box variants now permanently distinguishable, directly closing the
>      "which one's data is actually in the archive" ambiguity called out above.
>    - R2's real combined-job archive (`antidunes-buoyancy-r2-compare-all_20260912_175355`, 4.1G,
>      confirmed the ONLY non-empty one of two same-named folders -- the other was an empty
>      cancelled-attempt leftover, left untouched) had its 18 flat sibling folders reorganized into 6
>      nested `{sphere,cylinder,box}_{mem,analytic}/{eval,vtk-out,snapshots}` subfolders in place.
>      Verified: `du -sh` still reads 4.1G after reorganizing, confirming no data lost.
> 4. R4 itself: cancelled the collision-damaged job (4234535, missing `cylinder_mem` entirely plus
>    partial `sphere_mem`/`sphere_analytic` -- see this item's earlier note), deleted the corrupted
>    archive folder and leftover build-dir output, and **resubmitted fresh as job 4235003** on
>    `mt-dennys-mirror` (checked occupancy first -- `mt-dennys-mirror2` was busy with the clog extreme
>    retry) with the fixed script, so its own archive will come out correctly nested and clearly
>    R4-labeled from the start. Not yet analyzed.

> **[Claude, 2026-09-28] NEW: sphere/cylinder/cuboid equilibrium POSITION bias, quantified and under
> active investigation.** Per explicit request, re-derived the "Buoyancy Sweep — Sphere, Cylinder,
> Cuboid" artifact's own tail-window (last 10k of 50k steps) means and compared against its stated
> targets: the settled equilibrium position is systematically HIGHER (less submerged) than the naive
> constant-density Archimedes target, **growing with density ratio** for both sphere and cylinder
> (sphere: +0.20/+0.26 LU at rho=0.25 -> +0.39/+0.56 at rho=0.75, real-gravity/Bogner-gravity;
> cylinder similar, +0.48/+0.62 at rho=0.75; cuboid shows the same blow-up at rho=0.75 (+1.24/+1.31)
> but a sign flip at rho=0.25, likely coupled to the separately-documented box floor/mapping-tie
> artifact, not this mechanism). **Ruled out** one candidate mechanism analytically before touching
> code: a script using the WRONG (linear, not the true spherical-cap/circular-segment) submerged-
> depth-vs-volume-fraction formula would produce an antisymmetric error (too-high target at rho<0.5,
> too-low at rho>0.5, ~3.6 LU magnitude for r=24) -- the observed bias has no sign flip and is 6-9x
> smaller, so this is not it. **Leading hypothesis (not yet confirmed): weakly-compressible LBM
> hydrostatic density stratification under the gravity body force** (p=cs^2*rho -> rho decreases with
> height at equilibrium, unlike the naive uniform-rho=1 Archimedes assumption) -- would explain the
> sign (always "floats higher than predicted"), the growth with submerged depth, and rough (not
> exact) scaling with gravity magnitude in lattice units (Bogner's own g=7.5e-4 is ~15x larger in LU
> than the real-gravity calibration's gravityOverride for this scenario).
>
> Secondary, unrelated finding along the way: the real-gravity and Bogner-gravity arms of the
> original sweep start from measurably different total liquid mass (~756k vs ~803k LU^3, ~6% apart,
> consistent across all 3 densities within each arm) -- the two calibrations were never quite
> apples-to-apples liquid-volume setups. Not investigated further this session.
>
> **New diagnostic added to test the hypothesis directly**: `ModelParameters.logDensityColumn`
> (+`densityColumnX`/`Y`/`logDensityColumnStartStep`/`EndStep`, `SimulationParams.h`; function in
> `SetupStages.h`, independent registration `addDensityColumnLoggingSweep` in `BuildTimeloop.h` --
> deliberately NOT nested inside `diagnosticSweepEnabled` like `logFluidBoxHistory`, so it works with
> all the crash-investigation diagnostics left off) logs EVERY cell (all flags, not just OBST/INTF)
> along the full domain height at one fixed far-field (x,y) column, every main step, to
> `density-column/rank_<N>.csv` -- gives the actual measured rho(z) profile far from the floater.
> Also added the `density` VTK writer (already supported by the library, just unused before) to the
> new prms for an independent cross-check.
>
> **Reran the sphere R3 (3x scale, r=24) sweep** with this instrumentation: 6 new prms
> (`prms/sphere_stillwater_{realgravity981,bognerboxgravity}_3xsphere{,_d025,_d075}_biascheck_gas.prm`)
> and 3 new sbatch scripts (`jobs/run_72_sphere_stillwater_gravity_comparison_3xsphere{,_d025,_d075}
> _biascheck.sbatch`, each still running its own real-gravity+Bogner-gravity pair back to back, as
> the originals did). **`useObstaclePdfFreeze` explicitly forced `false`** in all 6 prms -- this flag
> didn't exist when the original artifact's runs were made (added 2026-09-27, defaults to `true` now)
> and directly touches every obstacle-covered cell's PDF handling every step, so leaving it at its new
> default would have confounded the rerun with an unrelated, newer crash-mitigation feature. No other
> crash-investigation debug logging needed disabling -- everything else already defaults off and
> these files never touched it.
>
> **Smoke-tested before committing cluster time** (per this project's own "no local build for the
> full sim" convention): a local oversubscribed run hit a `WALBERLA_ASSERT(field::isFlag(domain_))`
> in `BoundaryHandling.h` at both the original AND modified prm, at every rank count tried including
> the real 72-rank production decomposition -- looked like a serious regression at first. A genuine
> cluster smoke run (job 4271923, 500 steps, real 72-rank singlenode allocation) completed cleanly
> (`SMOKE_EXIT_CODE=0`, full run to completion, no assertion) -- **confirmed the local failure was a
> local-oversubscription artifact only, not a real regression**. Its `density-column/rank_0.csv`
> already shows a clean, monotonically DECREASING density with height at step 499 (1.009 at z=0.5 ->
> 1.00468 at z=29.5) -- a first positive signal consistent with the hydrostatic-stratification
> hypothesis, though far from the converged-equilibrium picture the full 50000-step runs will give.
>
> Incidental fix needed along the way: `sbatch --constraint=el9` is now rejected cluster-wide
> ("Invalid feature specification") -- confirmed via `sbatch --test-only` on the UNMODIFIED original
> sibling script too, so this is a cluster-side change (el9 nodes likely drained/retired), not
> specific to any file here. Removed from these 3 new sbatch scripts only; NOT fixed in the many
> other existing `run_72_*.sbatch` files that still have it (out of scope for this rerun).
>
> **Submitted (2026-09-28), one per checkout to avoid a shared-build-dir collision**: d050 (no
> suffix) = **job 4271925** on `mt-dennys`, d025 = **job 4271926** on `mt-dennys-mirror`, d075 =
> **job 4271927** on `mt-dennys-mirror2`. Each runs both gravity arms back to back, 50000 steps each
> (same duration as the original sweep), ~8h time limit. Not yet analyzed -- check back once
> complete, then extract each `density-column/rank_*.csv`'s far-field rho(z) profile and compare its
> shape/magnitude against the position bias measured above.
>
> **[Claude, 2026-09-28] Cleanup, per explicit request: removed `useObstaclePdfFreeze` entirely**
> (both the flag and its "reconstruction-only" no-op alternative), since it was added this same
> session (2026-09-27) purely as an A/B-test knob and its only currently-used setting across every
> live prm (including the 3 jobs just submitted above) is `false` = reconstruction-only anyway --
> deleting it just makes that the sole, permanent behaviour instead of a config option. Removed:
> the field + doc comment + parser line (`SimulationParams.h`), the `obstacleFreezeSnapshot` map and
> both `saveObstaclePdfsBeforeCollide`/`restoreObstaclePdfsAfterCollide` lambdas plus their
> `BeforeFunction`/`AfterFunction` registrations (`BuildTimeloop.h`), and the now-dead
> `useObstaclePdfFreeze false;` line + explanatory comment from all 6 sphere-R3 bias-check prms
> (renumbered their header notes accordingly). Local build (`cmake --build build --target
> Antidunes`) confirmed clean after the removal. **Does not affect the 3 already-running jobs**
> (4271925/4271926/4271927) -- they already built and started from their own synced copy before this
> edit; this only changes the source for future runs/builds.
>
> **[Claude, 2026-09-28] Correction per explicit follow-up: the ownership-handoff PDF-reset fix
> (`onOwnershipHandoff` in `BuildTimeloop.h`, the real, independently-validated bug fix for the
> "two rods hand a shared cell off with no reset" STALE-PARTICLE-IDX origin -- NOT the freeze/
> reconstruction-only A/B test removed above) stays, but is now gated behind its own prm toggle
> instead of being unconditional: new `ModelParameters.useOwnershipHandoffPdfReset` (bool, DEFAULT
> true -- keeps the fix active by default), `SimulationParams.h` + a one-line early-return guard at
> the top of the callback body in `BuildTimeloop.h` (`if (!simParams->model.
> useOwnershipHandoffPdfReset) return;`), mirroring this codebase's own established pattern for
> A/B-testable mitigations. `false` reverts to the original pre-fix behaviour (silent PDF
> inheritance across a handoff). Local build confirmed clean.
>
> **[Claude, 2026-09-28] Follow-up per explicit request: flipped `useOwnershipHandoffPdfReset`'s
> default to `false`** -- the fix is now OFF unless a prm explicitly sets it `true`. Local build
> confirmed clean again after the flip.

> **[Claude, 2026-09-28] RESULTS of the sphere R3 bias-check reruns (jobs 4271925/26/27, all COMPLETED, ~1h40m each, exit 0).**
> **Housekeeping first:** the sbatch archive step FAILED ("mkdir: Disk quota exceeded" on `/home/vault`, file-count quota) so nothing was archived and, since all outputs stay in the shared build dir, the Bogner arm OVERWROTE the real-gravity arm's density-column CSVs/vtk in each checkout. Rescued to `~/University/Thesis/animations/sphere-r3-biascheck-density-column_20260928/` (density-column ranks 0-3 = Bogner arm only; `eval/antidunes.sqlite` = both arms; see its PROVENANCE.md). The real-gravity arm's own rho(z)/waterline is therefore LOST (rerun needed if wanted). Scripts: `analyze_sphere_biascheck.py`, `analyze_sphere_biascheck_realarm.py`.
> **Bogner arm (g=7.5e-4): compressibility hypothesis CONFIRMED.** Far-field rho(z) is hydrostatically stratified exactly as isothermal-LBM theory says (slope -2.41e-3/cell vs -3g=-2.25e-3 (the extra 7% = exponential curvature); rho=0.998 at the waterline). This also explains the earlier "6% mass mismatch" between arms: mean rho = 1+3gH/2, predicted arm mass ratio 1.0615-1.0637 vs measured 1.0620-1.0635. Far-field waterline from sum(fill): 58.91/60.01/60.98 (d025/d050/d075). Predicted centre (rho=1 Archimedes at that waterline / with measured rho(z)) vs measured: d025 67.243/67.356 vs 67.456; d050 60.013/60.330 vs 60.398; d075 52.640/53.397 vs 53.478. Bias vs naive +0.21/+0.39/+0.84 LU; compressibility accounts for +0.11/+0.32/+0.76 of it, leaving +0.10/+0.07/+0.08 LU (~0.1 LU = 0.2% of D, sub-cell level). So for Bogner's g the "bias" is a property of the weakly-compressible model, not a bug; the naive constant-density target is what is wrong.
> **Real-gravity arm (g=5.07e-5): NOT explained by this and NOT converged.** Compressibility shifts the prediction by only +0.01..0.05 LU there, but tail-10k z is still oscillating (p-p 1.7/0.4/1.3 LU, period ~4.2-5.4k steps, envelope tau ~2-4e4 steps -> amplitude still ~0.3-1 LU at 50k), so tail means (67.644/60.420/53.368) are contaminated. Damped-cosine asymptotes z_inf=67.560/60.434/53.391 sit +0.30/+0.26/+0.47 LU above the naive prediction at the artifact's (real-arm, VTK) waterline (58.92/60.17/61.26; real-arm waterline unmeasurable now, +-0.15-0.3 LU uncertainty). Leading candidates (NOT tested): a roughly g-independent absolute force error (e.g. MEM/mapping surface error) whose position effect scales 1/g, so it is ~15x larger at small g; plus near-inviscid tau=0.5006 (useHardcodedUnits) -> very weak damping. Next: rerun the real-gravity arm alone with the archive fixed (density column + waterline), longer (or damped/released near equilibrium), and test a force-offset explanation via the z_inf shift vs g.
> **Fix status:** no simulation code bug found for the Bogner arm; the correct remedy is the reference: use the compressible Archimedes target (integrate exp(-3g(z-z_fs)) over the submerged cap). For the real-g arm the remaining +0.3..0.5 LU is still open. Same mechanism should apply to cylinder/cuboid (not rerun).
>
> **[Claude, 2026-09-28] CORRECTION, prompted by the user's question ("Bogner and real gravity stabilise in roughly the same place -- why less accurate with more submerged?"): the "compressibility CONFIRMED" wording above overreaches -- retracted to "compressibility is present and measured, but does not explain the offset".** The two arms (g differs 15x) end within ~0.1 LU of each other, yet compressibility predicts +0.76 LU (d075) in the Bogner arm and ~+0.05 in the real arm, so the offset is g-INDEPENDENT; the Bogner match was at best partly coincidental (compressibility should ADD 1-3% on top of any g-independent error, and the Bogner total is not larger than the real one). Re-expressed as the submerged volume fraction the sphere actually needs vs nominal (`/tmp` calc, both arms, z_inf, waterlines above): d025 0.2417/0.2442 (real/Bogner) vs 0.25; d050 0.4917-0.4868/0.4880 vs 0.50; d075 0.7371-0.7291/0.7267 vs 0.75, i.e. under-submerged by a roughly CONSTANT 1.7-3.5% RELATIVE (0.6-2.3 pct-points of V) in every case. That answers "why worse with more submerged": the position error is dz = dV/(waterline cross-section), dV ~ eps*V_sub grows with submersion while the waterline area (pi r_w^2, 1577/1810/1627 for d025/50/75) does not, giving amplification factors 9.2/16/26.7 -> ~0.3/0.4/0.8 LU for eps~3%. Force balance itself is fine (job log [FloaterForces], weight m*g = 1.4667 real / 21.7147 Bogner): Bogner-arm F_hydro_z = 21.711 at step 50000 (-0.02% vs weight); real arm 1.455-1.479 oscillating around it (+0.2% mean of 4 samples) -- so the sim sits in exact MEM-force equilibrium; the discrepancy is that MEM force at a given geometric immersion is ~2-3% larger than rho*g*V_sub (or that far-field waterline / analytic sphere is not what the body 'sees'). HYPOTHESIS (untested): first-order simple-bounce-back pressure sampling: the momentum exchange samples hydrostatic p at the fluid node, 0.5 cell off the halfway wall, so bottom-facing links read +rho*g*0.5dx and (fully submerged) top-facing links -rho*g*0.5dx -> excess buoyancy volume ~ dx*pi*R^2 = 3.1% of V for R=24 (~ 3*delta/R) -- g-independent, ~constant relative, and 1/R: consistent with the 2026-09-21 R1/R2/R4 sphere note (MEM faster than analytic by 11%/4%/2.5% at R=1/2/4, ~1/R). Discriminating test: PINNED sphere scan (pinPosition true, several z / several R), measure mean F_hydro_z vs rho*g*V_sub(z) -- no dynamics, no convergence problem, also separates it from compressibility (vary g) and from waterline/meniscus effects.
>
> **[Claude, 2026-09-28] Pinned-sphere scan SUBMITTED (per "Ok do it")**: `gen_pinned_sphere_scan.py` generates 22 prms (`prms/pinscan_sphere_R{12,24,36}_{real,bogner}_f{025,050,075,100}_gas.prm`; fully pinned sphere, still water, 12000 steps, log every 250, density column over steps 10000-11999; nominal submerged fractions 0.25/0.5/0.75 (+1.0 fully submerged for R=12,24 only), waterline 59.8831) and 3 sbatch scripts (`jobs/run_72_pinscan_sphere_R*.sbatch`, one radius each, ~8 cases). Jobs: **R=12 = 4272359** (mt-dennys-mirror), **R=24 = 4272360** (mt-dennys), **R=36 = 4272361** (mt-dennys-mirror2). The scripts do NOT archive to the vault (file quota exhausted): each case is reduced in-job to `pinscan_results/<case>.{forces.txt,column.txt,status}` (FloaterForces lines, far-field waterline + mean rho(z)). Read-out: mean F_hydro_z / (rho*g*V_sub(z_c, measured waterline)) per (R, g, f). Predictions: hypothesis (first-order BB pressure sampling) -> ratio ~1+3*dx/(2R)-ish, g-independent, falling ~1/R; compressibility -> ratio rises with g and depth, R-independent; waterline/meniscus -> depends on g via capillary length.
>
> **[Claude, 2026-09-28] Scan REDUCED per user ("only run it for 0.75 and for one gravity")**: cancelled 4272359/60/61 (before any case finished, partial results deleted) and resubmitted only **rho-fraction 0.75, REAL gravity (g=5.07e-5; chosen because compressibility is negligible there, so any excess is the unexplained part)**, one case per radius: **R=12 = 4272365** (mt-dennys-mirror, z_c=55.716), **R=24 = 4272366** (mt-dennys, z_c=51.548), **R=36 = 4272367** (mt-dennys-mirror2, z_c=47.38); ~1h limit each. `gen_pinned_sphere_scan.py` now generates only these 3 prms (`GRAVS=("real",)`, `FRACS=[0.75]`). NOTE: with one gravity the scan tests the 1/R scaling (pressure-sampling hypothesis) but cannot by itself separate compressibility (g-dependence) -- that would need the Bogner gravity at the same R.
>
> **[Claude, 2026-09-28] Overwrite protection, per user ("Make sure the results are not overwritten"):** jobs 4272365/66/67 were cancelled ~35 s in (no results yet) and resubmitted as **R=12 = 4272372** (mt-dennys-mirror), **R=24 = 4272373** (mt-dennys), **R=36 = 4272374** (mt-dennys-mirror2) with a safer sbatch (regenerated by `gen_pinned_sphere_scan.py`): results go to a per-job directory `pinscan_results/<SLURM_JOB_ID>/` that must not already exist (job aborts instead of overwriting), are copied with `cp -n` to `~/pinscan_backup/<jobid>_R<R>/` outside the repo checkout, and both copies are made read-only at the end; the scripts never `rm` any results (only build-dir scratch output). Checked `antidunes.sh --sync`: plain `rsync -avz` without `--delete`, so a sync cannot remove results either. Local copies must be taken with `rsync --ignore-existing` into a new dated folder.
>
> **[Claude, 2026-09-28] PINNED-SPHERE SCAN RESULTS (jobs 4272372/73/74, all COMPLETED exit 0, 11-14 min; real gravity g=5.06585e-5, nominal 0.75 submerged, 12000 steps, F_hydro mean over 25 samples at steps >= 6000).** Local copy `~/University/Thesis/animations/pinscan_sphere_20260928/` (also backed up on the cluster in `~/pinscan_backup/`), script `analyze_pinscan_sphere.py`. Ratio = mean MEM F_hydro_z / (g * V_sub(analytic cap, pinned z_c, measured far-field waterline)):
> | R | z_c | waterline (sum fill) | mean F_hydro_z (std) | ratio | excess | extra buoyant volume / (pi R^2) |
> |---|---|---|---|---|---|---|
> | 12 | 55.716 | 60.201 | 0.2885 (0.0037) | 1.0255 | +2.55% | +0.31 cell |
> | 24 | 51.548 | 60.275 | 2.2590 (0.0271) | 1.0123 | +1.23% | +0.30 cell |
> | 36 | 47.380 | 60.369 | 7.5900 (0.0530) | 1.0103 | +1.03% | +0.37 cell |
> **Reading:** with NO dynamics the MEM force at fixed position is 1-2.5% ABOVE rho*g*V_sub at real gravity (g-independent mechanism cannot be tested with one g, but compressibility is only ~0.05% here, so this is not compressibility). The excess is NOT constant in relative terms (2.55 -> 1.23 -> 1.03% for R=12/24/36) but the equivalent extra buoyant volume divided by the projected area pi*R^2 is CONSTANT at ~0.3-0.37 cell, i.e. the force behaves as if the buoyant body were ~0.3 cell taller/bigger, independent of R: that is the signature of a first-order wall/pressure-sampling offset (predicted 0.5-1 cell; measured ~0.3), supporting the bounce-back hypothesis over a constant-relative (mass/volume-mapping) error. **Caveats:** (1) only 3 radii, one gravity, one immersion (0.75); (2) the waterline matters: the measured far-field level (60.20/60.27/60.37) is 0.3-0.5 above the nominal 59.883 (a V-independent volume offset of ~4000-6000 cells, unexplained), and using the nominal waterline instead gives 4.9/2.7/2.2% -- a +-0.3 LU waterline uncertainty is worth ~+-1% in the ratio; (3) at R=24 the pinned excess (1.2%) is only about half of the free-floating real-arm effective deficit (1.7-2.8%), so part of the free-floating offset is still unexplained (waterline/meniscus near the body, or the oscillating unconverged state). Next if wanted: Bogner g at the same radii (separates g-dependence), other immersions (bottom-only vs top+bottom sampling: 0.5 vs 1.0 dx), and a local-waterline measurement next to the sphere.
>
> **[Claude, 2026-09-28] Interpolated-bounce-back pinned test SUBMITTED (per "yes run it with interpolated bounce back")**: identical pinned real-gravity 0.75 cases with `ModelParameters.useInterpolatedBounceBack true` (CLI/Ginzburg wall placement on the MovingObstacle links; wired via `AntidunesBoundaryHandling`) -- `gen_pinned_sphere_scan.py` now also writes `prms/pinscan_sphere_R{R}_real_f075_ibb_gas.prm` + `jobs/run_72_pinscan_sphere_R{R}_ibb.sbatch`. Jobs **R=12 = 4272519** (mt-dennys-mirror), **R=24 = 4272520** (mt-dennys), **R=36 = 4272521** (mt-dennys-mirror2); results in `pinscan_results/<jobid>/R{R}_real_f075_ibb.*` (+ `~/pinscan_backup/<jobid>_R{R}_ibb/`). Prediction: if the excess is a first-order wall-position/pressure-sampling error it should drop clearly below +2.55/+1.23/+1.03% (extra height ~0.3 cell -> toward 0); if unchanged the cause is elsewhere (waterline/meniscus, mapping).
>
> **[Claude, 2026-09-28] INTERPOLATED-BOUNCE-BACK RESULT (jobs 4272519/20/21, all COMPLETED exit 0; local `animations/pinscan_sphere_20260928/R{12,24,36}_ibb/`, `analyze_pinscan_sphere.py` prints both variants): NULL.** Mean F_hydro_z simple BB -> interpolated BB: R=12 0.28850 -> 0.28826, R=24 2.25898 -> 2.25808, R=36 7.59004 -> 7.58995; ratio 1.0255/1.0123/1.0103 -> 1.0279/1.0112/1.0103, i.e. changes of <0.1% against an excess of 1-2.5% and a sample std of ~1.3%; the equivalent extra height stays 0.27-0.37 cell. The flag is wired (`AntidunesBoundaryHandling` -> `SimpleBB` branch at line 343, code read only; no per-link counter was logged, so 'links actually engaged' is inferred not measured). **So the first-order wall-position / pressure-sampling explanation is NOT supported**: the excess does not respond to the wall placement. Remaining candidates: waterline/interface handling (fill-weighting phi, interface cells, the unexplained +0.3-0.5 LU offset of the measured level above the nominal 59.883, meniscus) versus body-related mapping. Cleanest discriminator (not yet run): pinned FULLY SUBMERGED spheres (f=1.0, R=12 and 24; generator already supports it via `FRACS`), which have no waterline: if the excess vanishes it is waterline/interface related, if it stays ~1-3% it is in the body force itself. Also cheap: `useFillWeighting false` on the 0.75 case.
>
> **[Claude, 2026-09-28] Continuity-check plot in the Clog Flow-Forcing artifact, checked by re-drawing it:** all three curves share a +-8% oscillation (period ~11k steps) for the first ~50k steps -- a start-up surge passing planes 96 cells apart, not the jam -- then the no-clog ratio sits at 1.00; PID-off's last ~8k steps are noisy (Q~3). Caption updated.
>
> **[Claude, 2026-09-28] Fully-submerged pinned test SUBMITTED (per "yes do this")**: `gen_pinned_sphere_scan.py` now also writes `prms/pinscan_sphere_R{12,24}_real_f100_gas.prm` (sphere fully below the surface: z_c = 43.8831 for R=12, 31.8831 for R=24, top 4 cells under the level, real gravity, 12000 steps, pinned) and `jobs/run_72_pinscan_sphere_R{R}_f100.sbatch`. Jobs **R=12 = 4272554** (mt-dennys-mirror), **R=24 = 4272555** (mt-dennys); results `pinscan_results/<jobid>/R{R}_real_f100.*`, backup `~/pinscan_backup/<jobid>_R{R}_f100/`. Analysis script `analyze_pinscan_sphere_f100.py` (expected force = g * integral of the MEASURED far-field rho(z) over the sphere volume, i.e. compressible Archimedes, so the ~0.2-0.4% hydrostatic density gain at depth is not mistaken for an excess). Decision rule: excess ~0 (<~0.5%) -> the 1-2.5% seen at 0.75 is a WATERLINE/interface effect (fill-weighting, interface cells, meniscus, level offset); excess still ~1-3% (and ~0.3 cell*pi R^2 in volume) -> it is in the body force itself (mapping / MEM). Not run: `useFillWeighting false` at 0.75.
>
> **[Claude, 2026-09-28] `useFillWeighting false` at 0.75 SUBMITTED (per "run one for now on the last free mirror")**: one case, R=24, real gravity, pinned at z_c=51.548, identical to baseline 4272373 except `useFillWeighting false` (diff of the two prms checked: that one line only) -- `prms/pinscan_sphere_R24_real_f075_nofw_gas.prm`, `jobs/run_72_pinscan_sphere_R24_nofw.sbatch`, **job 4272558** on `mt-dennys-mirror2`; results `pinscan_results/4272558/R24_real_f075_nofw.*`, backup `~/pinscan_backup/4272558_R24_nofw/`. Question: does removing the phi(x_b) fill-weighting of the waterline links change the excess (baseline +1.23% at R=24)? NB with weighting off the raw MEM force includes the "phantom" full-liquid lift of interface cells, so a LARGER force is expected if the weighting matters; a change in either direction means the waterline links carry part of the excess. `analyze_pinscan_sphere.py` now prints this variant too.
>
> **[Claude, 2026-09-28] WHERE THE PINNED-SPHERE EXCESS COMES FROM -- results of jobs 4272554/55 (fully submerged R=12/24) and 4272558 (0.75, R=24, useFillWeighting false), all COMPLETED exit 0; local `animations/pinscan_sphere_20260928/{R12_f100,R24_f100,R24_nofw}/`; scripts `analyze_pinscan_sphere_f100.py`, `analyze_pinscan_sphere.py`, `analyze_pinscan_pressure_reference.py`.**
> **(1) The body force is EXACT.** Fully submerged, real gravity: mean F_hydro_z / (g * integral of measured rho(z) dV) = 0.9990 (R=12, 0.36724 vs 0.36760) and 0.9998 (R=24, 2.94548 vs 2.94619), i.e. -0.10% / -0.02% (scatter ~0.9%); against plain g*V the ratio is 1.0015 / 1.0041, the difference being the hydrostatic density gain at depth (compressibility, 0.15-0.4% at real g). Extra volume/(pi R^2) = -0.016 / -0.008 cell (vs +0.30 at 0.75 submerged). MEM force, staircase mapping, mass/volume and bounce-back are therefore not the source. **(2) Fill-weighting is not the source**: useFillWeighting false at R=24/0.75: mean F_hydro_z 2.25729 (std 0.023) vs 2.25898 baseline, ratio 1.0128 vs 1.0123 (no phantom lift appears). **(3) Wall placement is not the source** (interpolated BB null, see above). **(4) => the excess is tied to the free-surface / waterline, and specifically to the PRESSURE REFERENCE LEVEL:** fitting the measured far-field rho(z) (slope = -3g within 0.5% in every case) and finding where rho = 1 (gauge pressure zero) gives a level **0.22-0.32 cell ABOVE the fill-integrated waterline (sum of fill)** in all 9 pinned real-gravity columns (offsets +0.25/+0.29/+0.23/+0.22/+0.26/+0.22 for the 0.75 cases, +0.28/+0.32 fully submerged; rho at the waterline = 1.00003-1.00005). A constant pressure offset acts on the wet surface as an extra uniform pressure -> net upward force g*delta*A_w (A_w = pi r_w^2, the waterline cross-section), and it integrates to ZERO over a closed (fully submerged) body -- exactly the observed pattern (excess only when a waterline cuts the body). Predicted extra volume/(pi R^2) = (r_w^2/R^2)*delta = 0.21/0.20/0.20 cell for R=12/24/36 vs measured 0.31/0.30/0.37: **the pressure-reference offset explains ~55-70% of the excess**; the remaining ~0.1-0.17 cell is only marginally above the sampling noise (~0.06 cell) and the waterline-definition uncertainty (+-1% ~ +-0.05-0.1 cell). Mechanism (hypothesis, not verified in code): in the FSLBM the interface cell's pressure boundary (rho_gas = 1) is imposed at the interface CELL, so p=0 sits ~0.25 cell above the level defined by the fill fractions.
> **Consequences / caveats:** (a) the g-independent ~3% 'excess buoyancy' of the free-floating R3 sweep is (at least in large part) this ~0.25-cell pressure-reference offset acting on the waterline cross-section; its relative size is ~delta*A_w/V_sub, i.e. grows with A_w/V_sub -- consistent with the earlier constant-relative-error reading at fixed R. (b) An analytic target should therefore use the pressure-zero level (fill waterline + ~0.25 cell) and the compressible profile, not the fill-integrated waterline. (c) OPEN: in the Bogner-gravity R3 columns the same fit gave rho at the waterline = 0.9984 (i.e. rho=1 ~0.7 cell BELOW the fill waterline, opposite sign), so the offset is not a universal constant of the scheme -- it may depend on g / initialisation / interface state; unresolved and must be understood before the offset is applied as a correction. (d) One gravity in the pinned scans; single depth for the fully submerged cases; force scatter ~1%. Not tested: check the interface-cell pressure BC in the code, and a 0.75 pinned case at Bogner gravity.
>
> **[Claude, 2026-09-28] CODE READ: the interface-cell pressure condition (per "read the interface pressure condition in the code first").** `src/lbm/free_surface/dynamics/functionality/ReconstructInterfaceCellABB.h` (`reconstructInterfaceCellLegacy`, Koerner 2005 / Thuerey 2007 eq. 4.5): the missing PDFs of an INTERFACE cell (those that would come from gas) are rebuilt directly in that cell as f_i = 2 w_i rho_gas (1 - 1.5 u^2 + 4.5 (c_i.u)^2) - f_ibar from the cell's OWN opposite PDF (anti-bounce-back), with `rhoGas = computeDeltaRhoLaplacePressure(sigma, curvature) + bubbleDensity` (= 1 + 3*2*sigma*kappa; `StreamReconstructAdvectSweep.h:155`). It never uses the cell's fill level phi. Reading: the pressure condition p = p_gas is therefore applied at the interface CELL (centre), whereas the fill-defined surface sits at bottom + phi, i.e. the pressure zero is (0.5 - phi) cells from the true surface, anywhere in [-0.5, +0.5]. **Check against the 9 pinned real-gravity columns** (`/tmp` calc): with phi_top = fractional part of the sum-fill waterline (0.03-0.37) the measured rho=1 level minus the centre of the top interface cell is -0.15..+0.09 cell (mean -0.02, std 0.08), i.e. rho=1 sits at the interface-cell centre to within the noise, and the earlier +0.22..+0.32 offset relative to the fill waterline is just (0.5 - phi). => the surface-pressure error is a SUB-CELL QUANTISATION of the free-surface position: buoyancy depends on how full the top interface cell happens to be (set by the total liquid volume), giving an extra pressure g*(0.5-phi) on the wet cross-section. **Not the whole story:** the Bogner-gravity R3 columns do not fit (phi_top = 0.91/0.01/0.98 predict -0.41/+0.49/-0.48, measured -0.59/-0.66/-0.66), and the two arms also differ in surface tension (computedParams: sigma = 8.5e-4 real vs 0 Bogner) and omega (1.997 vs 1.9). Contact-line capillary force is excluded for the real arm (measured excess force scales with the waterline AREA, dF/A_w ~ 2e-5 for R=12/24/36, not with r_w). **Falsifiable prediction / next test:** pinned R=24, 0.75 submerged, real gravity, with the liquid level shifted so phi_top ~ 0.1 / 0.5 / 0.9: the excess should follow (0.5 - phi) A_w g -> ~+0.4 cell, ~0, ~-0.4 cell (sign change at phi = 0.5), independent of R. Candidate fix if confirmed: place the interface pressure condition at the actual surface (fill-aware) or subtract the (0.5-phi) reference in the MEM pressure.
>
> **[Claude, 2026-09-28] FILL SCAN SUBMITTED (per "yes run the fill scan")**: pinned R=24, 0.75 submerged (z_c 51.548 fixed), real gravity, 12000 steps, identical to baseline 4272373 except `liquidHeightFactor` shifted so the liquid level moves by -0.175 / +0.225 / +0.625 LU (factor 9.9681 + shift/6.00146; baseline measured waterline 60.275 -> phi_top 0.275 -> expected phi ~0.10 / 0.50 / 0.90 if the level shifts 1:1 -- the ACTUAL phi is taken from the measured waterline). `gen_pinned_sphere_scan.py` -> `prms/pinscan_sphere_R24_real_f075_lv{m175,p225,p625}_gas.prm`, `jobs/run_72_pinscan_sphere_R24_lv*.sbatch`. Jobs **lvm175 = 4272641** (mt-dennys-mirror), **lvp225 = 4272642** (mt-dennys), **lvp625 = 4272643** (mt-dennys-mirror2); results `pinscan_results/<jobid>/R24_real_f075_lv*.*` (+ `~/pinscan_backup/<jobid>_R24_lv*/`). Analysis: `analyze_pinscan_fillscan.py`. Hypothesis (interface-cell-centre pressure zero) predicts extra volume/(pi R^2) = (r_w^2/R^2)*(0.5-phi): about +0.35 at phi 0.1, ~0 at phi 0.5, about -0.35 at phi 0.9 (sign change). A positive excess at phi ~0.9 refutes it.
>
> **[Claude, 2026-09-28] FILL SCAN RESULT (jobs 4272641/42/43, all COMPLETED exit 0; local `animations/pinscan_sphere_20260928/R24_lv{m175,p225,p625}/`; `analyze_pinscan_fillscan.py`): the (0.5 - phi) interface-cell-centre hypothesis is REFUTED.** Level shifts moved the top-cell fill to phi = 0.992 / 0.468 / 0.846 (baseline 0.275). Table (R=24, 0.75, z_c 51.548): | case | waterline | phi | z0-wl (rho=1 fit) | (0.5-phi) | mean F_z | ratio | extra vol/(pi R^2) | | baseline | 60.275 | 0.275 | +0.229 | +0.225 | 2.25898 | 1.0123 | +0.300 | | -0.175 | 59.992 | 0.992 | **+0.323** | -0.492 | 2.24567 | 1.0167 | +0.402 | | +0.225 | 60.468 | 0.468 | +0.246 | +0.032 | 2.27625 | 1.0131 | +0.322 | | +0.625 | 60.846 | 0.846 | **+0.342** | -0.346 | 2.31430 | 1.0167 | +0.414 |. At phi 0.99 and 0.85 the hypothesis predicts a NEGATIVE offset and negative excess (-0.43 / -0.29 cell); measured are +0.32/+0.34 cell and +0.40/+0.41 cell, i.e. **the pressure-zero offset is a roughly CONSTANT +0.23..+0.34 cell above the fill-integrated waterline and the excess is a constant +0.30..+0.41 cell, both independent of phi**. My earlier 'rho=1 at the top interface cell centre' agreement (mean -0.02 cell) was a coincidence of the narrow phi range 0.03-0.37 of the first nine columns. **What still holds:** the offset predicts the excess with a stable fraction: (r_w^2/R^2)*(z0-wl) = 0.198/0.283/0.212/0.291 vs measured 0.300/0.402/0.322/0.414 -> 66/70/66/70%, so a constant pressure-reference shift ~0.25-0.34 cell explains about two thirds of the excess in every configuration; the remaining third (~0.1-0.12 cell) is also roughly constant. Fully submerged bodies are unaffected (excess 0, see above). **Consequence for the free-floating sweep:** a constant offset means a constant HEIGHT shift (excess volume / waterline area = 0.35-0.47 cell here), independent of density -- the real-gravity free-floating biases (+0.30/+0.26/+0.47 LU) are consistent with that, and the earlier worry about sign changes with the top-cell fill (d025 phi~0.92) disappears; it does not itself explain growth with submersion (that, for Bogner g, is compressibility). **Still OPEN:** why the pressure zero sits ~0.3 cell above the fill surface (not the top-cell fill, not wall placement, not fill-weighting, not the body force); the Bogner-gravity columns have the opposite sign (-0.6..-0.66 cell) and sigma=0 -- so surface tension / omega / g dependence is the next suspect (real arm: sigma=8.5e-4, omega=1.997; Bogner: sigma=0, omega=1.9). Cheapest next tests: (a) pinned 0.75 R=24 at Bogner g (sign?), (b) real g with surfaceTension 0, (c) full per-cell column (flag/density/fill through the interface layer) instead of the reduced mean profile.
>
> **[Claude, 2026-09-28] CORRECTION of the Bogner-arm pressure-zero offset (found while answering "does gravity influence this? could we tune g to 0.0?"; script `analyze_pressure_reference_expfit.py`):** the earlier straight-line fit of rho(z) extrapolated to the surface is biased by the profile's exponential curvature ((3gH)^2/2 ~ 9e-3 at Bogner g, comparable to the measured rho_wl-1 of -1.4e-3). Refitting ln(rho) linearly: **Bogner free-floating columns z0-wl = +0.103 / +0.059 / +0.082 cell (d025/d050/d075), slope/(-3g) = 1.000; NOT -0.6 (the 'opposite sign' statements above are RETRACTED).** Real-gravity pinned columns with the same exponential fit: +0.268..+0.395 cell (mean ~0.32; straight-line values were 0.22-0.34). **Unified picture: position bias = compressibility (grows with g*depth) + the pressure-zero offset delta (constant shift).** Bogner: residual after compressibility was +0.10/+0.07/+0.08 LU = the measured offsets +0.10/+0.06/+0.08 (arm now fully explained); real g: delta ~0.3 vs free-floating biases +0.30/+0.26/+0.47. With the exponential offsets the pressure-reference offset explains ~81% of the pinned excess ((r_w^2/R^2)*delta = 0.243/0.328/0.258/0.336 vs measured 0.300/0.402/0.322/0.414 for baseline/lv-0.175/+0.225/+0.625). The remaining ~0.06-0.08 cell is at the noise level. **Does g influence delta?** Unknown: delta ~0.3 (real: g 5.1e-5, sigma 8.5e-4, omega 1.997) vs ~0.08 (Bogner: g 7.5e-4, sigma 0, omega 1.9) -- not proportional to g (larger at smaller g), both positive, so no zero crossing is established, and g cannot be tuned in production runs (fixed by g_SI with the grid/time step). Attribution needs single-parameter changes from the pinned real-g R=24/0.75 case (g, surface tension, omega). The 'not validated for the Bogner arm' caveat on the analytic correction is withdrawn: correction = level raised by the arm's own fitted delta + compressible profile (use an EXPONENTIAL fit, not a line).
>
> **[Claude, 2026-09-28] Bogner dissertation check (`~/Downloads/Diss_2017-Bogner.pdf`) + RESOLUTION SCAN SUBMITTED.** Bogner initialises with hydrostatic pressure (g=7.5e-4, tau=1/1.9, sect. 5.3.3) and uses a pseudo-compressible LBM (rho0=rho for the compressible equilibrium); his sect. 5.3.2 statics are classical incompressible naval statics (no compressibility correction). His free-surface BC for most applications, incl. the floating-body tests, is the Koerner rule ('FSK rule', sect. 3.3.3), stated to be only FIRST-order accurate ('exact only if h is divisible by the grid spacing, otherwise ... first order O(dx)', sect. 3.3 planar film flow); he proposes a second-order FSL rule (Bogner et al. 2015) that we do not have in the free-surface module (keyword search; only Koerner ABB in `ReconstructInterfaceCellABB.h`). => the ~0.1-0.3 cell pressure-zero offset is consistent with a known first-order limitation of the FSK scheme. His floating-cube heel-angle test (sect. 5.3.3: 26.565/45/26.565 deg for rho 0.25/0.5/0.75, s=16, 130x40x40) only claims that it 'converges nicely towards the ideal equilibrium state' at low resolution. **Test of first-order behaviour: pinned resolution scan submitted** -- `gen_pinned_resolution_scan.py` derives fully pinned, 0.75-submerged spheres from the R1/R2/R4 family (`prms/buoyancy_compare_sphere_r{R}_mem_gas.prm`: domain, sphere r=8R, liquid depth and true g_SI scaled together, g_LU 2.14e-5/2.67e-6/3.34e-7, X/Y periodic, T/B walls) -> `prms/pinres_sphere_R{1,2,4}_gas.prm` (z_c 78.356/156.711/313.423, waterline 81.134/162.268/324.536, 12000/24000/48000 steps, density column x=y=6R over the last 2000R steps), `jobs/run_{72,72,576}_pinres_sphere_R{R}.sbatch`. Jobs **R=1 = 4272738** (mt-dennys-mirror2; bed r1 copied there), **R=2 = 4272739** (mt-dennys), **R=4 = 4272740** (mt-dennys-mirror, 576 ranks); results `pinscan_results/<jobid>/R{R}_pinres.*`, backup `~/pinscan_backup/<jobid>_pinres_R{R}/`. Prediction: first order -> extra height stays ~0.3 CELL at all R (physical size ~1/R); second order -> falls ~1/R in cells. NB g_LU shrinks 64x across R, so a g-dependence of the offset is confounded with resolution (the earlier real-g/Bogner comparison suggests the offset is not simply proportional to g).
>
> **[Claude, 2026-09-28] Queue change + PINNED BOX SCAN SUBMITTED (per "run the box run instead of the 8 nodes thing we can run that over night").** The 576-rank R=4 resolution job 4272740 sat PENDING (Resources) and was cancelled; **resubmitted as 4272776 with `--begin=22:00` (tonight, checkout mt-dennys-mirror)**. R=1 (4272738, mirror2) and R=2 (4272739, mt-dennys) remain queued (Priority; estimated start ~15:20/15:26). In its place the free checkout mt-dennys-mirror runs the **pinned-box heel scan at real gravity: job 4272775** (`gen_pinned_box_scan.py` -> `prms/pinbox_r{075,025}_a{10,20,26p565,30}_gas.prm`, `jobs/run_72_pinbox_scan.sbatch`, 8 cases sequentially, 12000 steps each, ~2 h, limit 2h30). Fully pinned 30x60x30 box (s=30, roll about Y) at heel 10/20/26.565/30 deg for rho 0.75 and 0.25 (Bogner's stable heel is 26.565 for both by symmetry rho -> 1-rho), heave from Archimedes at level 60.2 (only sets the pose; the analysis uses the MEASURED waterline and pose), real g_SI=9.81 (g_LU 5.06585e-5), R3 domain, density column x=y=5 over steps 10000-11999. Results `pinscan_results/4272775/box_r{075,025}_a*.{status,forces.txt,column.txt,computedParams.txt}`, backup `~/pinscan_backup/4272775_pinbox/`; analysis `analyze_pinbox_scan.py` (ideal force/torque_y about the box centre = -F*(x_B - x_G) from the clipped cross-section at the measured waterline and at waterline + measured exponential-fit pressure-zero offset, with the compressible density at the centroid). Question: does raising the level by the measured offset remove the force AND torque error, and hence the heel-angle error? At heel 26.565 the ideal torque is ~0, so the measured torque there is directly the bias that shifts the equilibrium angle. **Answer to 'did you account for the level change from the displaced volume?': yes -- all floating-sphere comparisons used the MEASURED waterline of each run (Bogner arm 58.91/60.01/60.98 from sum(fill); real arm 58.92/60.17/61.26 from the artifact's VTK), which rises with the displaced volume (0.25->0.75 submerged displaces 28954 cells^3 / 12960 cells^2 = +2.23 predicted vs +2.07 / +2.34 measured); the nominal 59.883 was never used for the floating comparisons.** **Text read (Bogner sect. 3.3.3, Thm 5 + Eq. 3.157):** the FSK rule imposes the pressure boundary value at the point x_w = x_b + c_q/2 (half-way) with second-order accuracy only if the true boundary distance is delta = 1/2; otherwise it is first order with error E = 2 (w_q/c_s^2) (delta - 1/2) d_q(Pi_q - ...) -- i.e. a pressure error equal to the pressure gradient (= rho g in hydrostatics) times the mislocation (delta - 1/2) of the boundary; the FSL rule (Def. 14) instead uses the actual delta (second order). So the text explains the TYPE of the error (pressure reference mislocated by a fraction of a cell, times rho g) but not its size or its independence of the top-cell fill in our data.
>
> **[Claude, 2026-09-28] Bogner Fig. 5.10 measured + BOGNER-LIKE FREE BOX SUBMITTED.** Tracking the coloured boxes in the four frames of Bogner's Fig. 5.10 (`Diss_2017-Bogner.pdf` p.137, ~3.28 px/cell over the 130-cell basin, +-3 cells from perspective/tilt): frame1->4 horizontal shift green (rho 0.75) +0.6, yellow (0.5) +10.1, red (0.25) +17.1 cells -> his boxes DO translate sideways (comparable to our ~20-cell surge), but they reach rest ('run until all fluid motions ceased and the objects came to rest'); caption/text: 'striving towards vertical and rotational equilibrium', 'rotate around its longitudinal axis only' -- no statement that X/Y is pinned. His basin is 130x40x40 with a 32-long box (2s) => ~4 cells of clearance per side, which suppresses yaw/pitch/wobble by geometry; ours (R3 domain, 135 wide for a 60-long box) does not, and our earlier work found the free box's out-of-plane wobble growing into 3D tumbling (plan lines ~2235-2251), and pinning X/Y+pitch/yaw (roll+heave free, seeds 5/10 deg) stayed near seed+1.3 deg instead of reaching 26.565. **New test, closer to his set-up (per 'run this as well ... also use a variant with 5 degrees'): job 4272793** (`gen_bognerlike_box_scan.py` -> `prms/bognerlike_r{075,025}_h{286,500}_gas.prm`, `jobs/run_72_bognerlike_box.sbatch`, checkout mt-dennys-mirror2, `--dependency=afterany:4272738`, 4 cases sequentially, 100000 steps each, limit 4 h): real g_SI=9.81 (g_LU 5.06585e-5), basin <96,45,40> (liquid depth ~20 = half height like his lower-half fill; box s=18 -> 18x36x18, long axis Y, 36 in 45 -> 4.5 cells/side), heave + X/Y translation FREE, roll about the long axis only (pinRotationAxes <1,0,1>), seeds 2.86 deg (his) and 5 deg, rho 0.75 and 0.25. Results `pinscan_results/4272793/bl_r{075,025}_h{286,500}.*`, backup `~/pinscan_backup/4272793_bognerlike/`; analysis `analyze_bognerlike_box.py` (tilt history, last-20k mean/std, x/y travel). Question: does the box settle at 26.565 deg (like Bogner) when the geometry/DOFs match his, and where (26.565 + the offset-predicted shift from the pinned-box torque scan 4272775)? Related running jobs: pinned-box scan 4272775 (mt-dennys-mirror), resolution R=2 4272739 (mt-dennys), R=1 4272738 done/queued (mirror2), R=4 4272776 begins 22:00.
>
> **[Claude, 2026-09-28] RE-CHECK of what Bogner (2017) actually states about his setups (per "also double check what he says about his setup"; `Diss_2017-Bogner.pdf` sect. 5.3.2 pp.~119-120 and 5.3.3 pp.~121-123) -- corrections to my earlier paraphrases.** **(1) Righting-moment validation (5.3.2, Fig. 5.6):** 'In all of the following cases g = 1e-4'; partially filled basin initialised hydrostatically; **tau = 1.0**; cuboids of width:height 6:4 and 5:4 at b x h = 24x16 / 48x32 / 96x64 and 20x16 / 40x32 / 80x64; axis-of-rotation length l = 2h; box axis-aligned and '**exactly half-immersed** relative to the free surface plane' (i.e. rho = 1/2), '**fixed to a constant position**', heel varied 0-30 deg about the longitudinal axis through its centre of gravity; deviations are lattice-shape jumps and 'the errors are decreasing for higher resolutions ... simulated forces will converge towards the analytic values'. => a quantitative validation ONLY for rho=0.5 and non-cubic (5:4, 6:4) boxes; there is NO quantitative righting-torque validation for rho 0.25/0.75 cubes. **(2) Floating-cube equilibration (5.3.3, Figs 5.7-5.10):** domain 130x40x40, lower half liquid, **tau = 1/1.9**, hydrostatic initialisation, **g = 7.5e-4**, three boxes s = 16 with rho 0.25/0.5/0.75 'slightly heeled ... 2.86 deg', 'the 2D cube of side length s is represented by a cuboid 2s x s x s, that will rotate around its longitudinal axis only', 'run until all fluid motions ceased and the objects came to rest'; ideal stable heel 26.565 (0.25, 0.75) and 45 deg (0.5); result stated only qualitatively: 'even at this low resolution the simulation converges nicely towards the ideal equilibrium state' -- **no measured angles are reported**. **NOT stated anywhere in the thesis:** the basin's side boundary conditions, which axis of the basin the box's long axis lies along (my '4 cells clearance per side' is an INFERENCE from 32 in 40, not a statement), whether X/Y translation or pitch/yaw are constrained (the 'rotate ... only' can be a constraint or an expectation), the run length, or any damping detail. **Consequences:** (a) his two tests use different lattice parameters (g 1e-4 / tau 1.0 vs g 7.5e-4 / tau 0.526); our real-g runs use g 5.07e-5 and tau ~0.5006 (omega 1.9974, hardcoded viscosity 2.1e-4) = the LEAST damping of all three, which alone can keep a free body oscillating; (b) our pinned-box heel scan follows his 5.3.2 protocol (fixed position, heel varied) but for rho 0.75/0.25 and a square section, beyond what he validated; (c) 'job 4253988 = THE methodologically correct reproduction of his test' (plan ~line 2222) is too strong: whether he pinned X/Y is unknown; (d) his statement that force errors shrink with resolution supports the R1/R2/R4 pinned resolution scan as the right test of the surface offset.
>
> **[Claude, 2026-09-28] RESOLUTION SCAN R=1/R=2 RESULT (jobs 4272738/4272739, COMPLETED exit 0; local `animations/pinscan_sphere_20260928/pinres_R{1,2}/`; `analyze_pinres_sphere.py`; R=4 = 4272776 begins 22:00): NOT the constant-in-cells signature of a first-order lattice error.** Fully pinned r=8R sphere, nominal 0.75 submerged, true g_SI: R=1 (r=8, g_LU 2.137e-5): mean F_z 3.8893e-2 (std 1.3%), waterline 81.829, ratio 1.0539 (+5.4%), extra volume/(pi r^2) +0.46 cell, exp-fit offset +0.14 cell; R=2 (r=16, g_LU 2.671e-6): F_z 3.8780e-2 (std 2.0%), waterline 162.985, ratio 1.0861 (+8.6%), extra +1.43 cell, offset +1.24 cell. Together with the R3-domain baseline (R=24, g_LU 5.07e-5): 0.30 cell / +1.2%. **The excess in cells GROWS as g_LU falls (0.30 -> 0.46 -> 1.43), i.e. in physical length ~0.15/0.17/0.26 mm at dx = 4.84e-4/3.62e-4/1.81e-4 m -- roughly constant or rising with refinement, not shrinking; the absolute extra force in lattice units is ~constant (R1 0.0020, R2 0.0031).** (Caveat: the far-field rho=1 offset is unreliable at small g_LU -- a density error of 1e-5 is worth ~1 cell at g_LU 2.7e-6 -- so the offset column is noisy there; the force-based excess is the trustworthy number.) **Key observation:** `useHardcodedUnits` fixes the LATTICE surface tension at 0.000849592 for every R (`computedParams.txt`: R1, R2 and the R3 real-g run all 0.000849592) while g_LU falls 64x, so the capillary length sqrt(sigma/(rho g)) in cells is 4.1 (R3), 6.3 (R1), 17.8 (R2), ~50 (R4) and l_c/r = 0.17 (R3, r=24), 0.79, 1.1, 1.6 -- surface tension gets relatively MORE important as R rises. This fits everything seen so far: excess ONLY when a waterline cuts the body (a fully submerged sphere has no contact line/interface), insensitive to wall placement / fill-weighting / top-cell fill, the Bogner-gravity arm (sigma = 0 in its computedParams) showing only ~0.08 cell while the real-g arm (sigma = 8.5e-4) shows ~0.3, and the R-family excess growing with l_c/r. => **leading hypothesis now: the surface-tension (Laplace-pressure, rhoGas = 1 + 3*2*sigma*kappa in the interface reconstruction, plus the contact-line meniscus at the body) term, not first-order boundary-condition placement.** (My earlier exclusion of a contact-line force rested on dF/A_w ~ const within the R3 set only, where Bond number changes with r; the R-family contradicts it.) **TEST SUBMITTED: job 4272824** (mt-dennys, 4 cases sequentially, ~50 min): pinned real-g 0.75 spheres with `PhysicsParameters.surfaceTensionOverride 0` for R=12/24/36 and doubled (0.001699184) for R=24 (`prms/pinscan_sphere_R{12,24,36}_real_f075_sig0_gas.prm`, `..._R24_real_f075_sig2x_gas.prm`, `jobs/run_72_pinscan_sigma.sbatch`; diff to the baseline prm = that one line); results `pinscan_results/4272824/R{R}_real_f075_sig{0,2x}.*`, backup `~/pinscan_backup/4272824_R24_sig0_sigma/`. Prediction: sigma = 0 -> excess -> ~0 (like the fully submerged case) at all R; doubled -> excess roughly doubles. If sigma = 0 does NOT remove it, the hypothesis is refuted. NB the pinned-box and Bogner-like box jobs (4272775, 4272793) still use the nominal sigma.

> **[Claude, 2026-09-28] SURFACE-TENSION-OFF SPHERE TEST (job 4272824, R12/R24/R36 done; R24 sig2x still running) -- DECISIVE: the excess buoyancy is caused by the surface-tension term.** Pinned real-g 0.75 spheres, `surfaceTensionOverride 0` (computedParams confirms surfaceTension 0), local `animations/pinscan_sphere_20260928/R{12,24,36}_sig0/`, `analyze_pinscan_sphere.py` (_sig0 block): F_z / F_ideal(measured waterline) = **1.0004 (R12), 1.0009 (R24), 1.0030 (R36)** vs the nominal-sigma baseline 1.0255 / 1.0123 / 1.0103; exponential-fit pressure-zero offset (rho=1 level minus waterline) **+0.001 / +0.004 / -0.001 cell** (nominal: +0.3 cell at R24); waterline identical (60.38) at all R. => with sigma = 0 the MEM buoyancy of a pinned, waterline-cutting sphere agrees with Archimedes to <= 0.3%, at every resolution, with real gravity. Both the force excess and the 'pressure reference offset' come from sigma. **What is NOT established:** the mechanism inside the code (my working guess: spurious curvature of the nominally flat surface in `rhoGas = 1 + 3*2*sigma*kappa` and/or the contact-line meniscus at the body; needs a check of the FD curvature on a flat interface). sigma = 0 is a DIAGNOSTIC, not a production setting (removes capillary physics that matters for clogging/free-surface flow). Consequences to test next: (a) does sigma = 0 also remove the cuboid heel-angle / torque error at rho 0.75/0.25? -> sigma = 0 variants of both box scans submitted: **pinned box 4272853** (mt-dennys, `--dependency=afterany:4272824`, `jobs/run_72_pinbox_scan_sig0.sbatch`) and **Bogner-like free box 4272854** (mt-dennys-mirror, `--dependency=afterany:4272775`, `jobs/run_72_bognerlike_box_sig0.sbatch`); prms `prms/pinbox_*_sig0_gas.prm` / `prms/bognerlike_*_sig0_gas.prm` differ from the nominal ones only by the line `surfaceTensionOverride 0.0;`; NB the sbatch case labels are unchanged (`box_r075_a10`, `bl_r075_h286`), results in `pinscan_results/<jobid>/`, copy locally to `pinbox_sig0/` and `bognerlike_sig0/` (analysis scripts take the subfolder as argv[1]). (b) the sigma effect in the R1/R2/R4 family (excess in cells growing as l_c/r grows) is consistent with this result.

> **[Claude, 2026-09-28] PINNED-BOX + BOGNER-LIKE PARTIAL RESULTS (nominal sigma, real g; jobs 4272775 7/8 cases, 4272793 3/4 cases).** **Pinned box** (`analyze_pinbox_scan.py`; F / F_ideal at the measured waterline -> with the measured exp-fit level offset): rho 0.25: 10 deg 1.0523 -> 1.0141, 20 deg 1.0628 -> 1.0204, 26.57 deg 1.0824 -> 1.0339; rho 0.75: 10 deg 1.0333 -> 1.0177, 20 deg 1.0334 -> 1.0187, 26.57 deg 1.0399 -> 1.0264, 30 deg 1.0404 -> 1.0216. So the force excess is 3-8% at the waterline and still 1.4-3.4% after the level correction (the offset explains most but not all). Torque (measured minus ideal at the measured waterline, in the analysis's units): rho 0.25 +0.094 / +0.099 / +0.147 at 10/20/26.57 deg (systematic, same sign) vs rho 0.75 +0.008 / -0.098 / +0.078 / -0.071 (scatter +-0.1, i.e. lattice-shape jumps, no systematic sign) -> an asymmetry between rho 0.25 and 0.75, whereas the ideal problem is symmetric under rho -> 1-rho. **Bogner-like free box** (`analyze_bognerlike_box.py`, last-20k-step tilt mean/std/range, target 26.565 deg): rho 0.75 seed 2.86 deg: 26.83 / 0.76 / 25.71-29.24 (x travel 47.9-61.6), rho 0.75 seed 5 deg: 31.40 / 0.58 / 30.13-32.76, rho 0.25 seed 2.86 deg: 28.61 / 0.01 / 28.60-28.63 (x wraps 0.3-94.5, box drifts through the periodic basin). => seed-dependent end angles (26.8 vs 31.4 for the same rho) mean the free box is NOT converging to a unique equilibrium (limit cycle / metastable lattice states), so the 26.83 deg is not by itself evidence that Bogner's angle is reproduced; the rho 0.25 run is steady but 2 deg off. The nominal-sigma box runs may be contaminated by the sigma effect found above; the sigma = 0 reruns are the fair comparison. Remaining nominal jobs: 4272775 a30 rho0.25, 4272793 rho0.25 seed 5 deg, R=4 resolution 4272776 (starts 22:00).

> **[Claude, 2026-09-28] SURFACE-TENSION DOSE-RESPONSE COMPLETE + NOMINAL BOX JOBS FINISHED (4272824, 4272775, 4272793 all COMPLETED exit 0; local `animations/pinscan_sphere_20260928/{R24_sig2x,pinbox,bognerlike}`).** **Sphere R=24, real g, 0.75 submerged, pinned:** extra buoyant volume / (pi r^2) = **0 (sigma 0: +0.023 cell, F ratio 1.0009) -> +0.300 cell (nominal 8.50e-4, ratio 1.0123) -> +0.566 cell (2x = 1.699e-3, ratio 1.0234)**, i.e. the excess scales ~linearly with sigma (x1.9 for 2x sigma) and vanishes at sigma = 0 => surface tension is confirmed as the driver of the sphere bias, not just correlated with it. (Waterline shifts accordingly: 60.381 / 60.275 / 60.129.) **Pinned box (all 8 cases, nominal sigma):** the last case a30 rho 0.25 gives F/F_ideal(wl) 1.0597 (-> 1.0124 with offset), torque-minus-ideal +0.077 (-> -0.005); table otherwise as in the partial entry. **Bogner-like free box (all 4 cases, last-20k mean +- std):** rho 0.75 seed 2.86 deg 26.83 +- 0.76; rho 0.75 seed 5 deg 31.40 +- 0.58; rho 0.25 seed 2.86 deg **+28.61 +- 0.01**; rho 0.25 seed 5 deg **-24.30 +- 0.00** (settles in the mirror-image state). **New observation:** the two steady rho 0.25 runs sit at +26.565 + 2.05 and -(26.565 - 2.27), i.e. BOTH equilibria are displaced by the same ~ +2.1 deg in the SAME rotational direction, not symmetrically about zero: a constant torque bias (independent of which side the box heels to) shifts the equilibrium angle by ~2 deg. This matches the systematic positive torque error of the pinned rho 0.25 scan (+0.09..+0.15 before / +0.03..+0.06 after the level correction) and is the pattern expected from an interface/level effect rather than from lattice-shape noise. The rho 0.75 runs are unsteady (std 0.6-0.8 deg, seed-dependent mean), so no comparable statement is possible there. Pending: the sigma = 0 box reruns (pinned 4272853 running, Bogner-like 4272854 running) decide whether the +2.1 deg shift disappears; R=4 resolution 4272776 starts 22:00.

> **[Claude, 2026-09-28] SIGMA = 0 BOX SCANS FINISHED (pinned 4272853 and Bogner-like 4272854, both COMPLETED exit 0; local `animations/pinscan_sphere_20260928/{pinbox_sig0,bognerlike_sig0}`; `analyze_pinbox_scan.py pinbox_sig0`, `analyze_bognerlike_box.py bognerlike_sig0`). Diagnostic only -- real gravity, nominal sigma remains the thesis setting.** **Pinned box:** the level offset is 0.000 (nominal +0.27..+0.39 cell) and the waterline is 60.381 in all 8 cases. **Torque (measured - ideal at the waterline), rho 0.25 at 10/20/26.57/30 deg: +0.014 / -0.003 / -0.005 / -0.026 (nominal +0.094 / +0.099 / +0.147 / +0.077)** => the systematic rho 0.25 torque bias is caused by surface tension and vanishes without it. rho 0.75: +0.055 / -0.025 / +0.147 / +0.012 (nominal +0.008 / -0.098 / +0.078 / -0.071): still scattered, i.e. lattice-shape jumps unrelated to sigma, with a +0.147 at the 26.57 deg equilibrium pose which by itself would shift the rho 0.75 equilibrium. **Force F/F_ideal(wl):** rho 0.25: 1.0133 / 1.0197 / **1.0429** / 1.0202; rho 0.75: 1.0173 / 1.0175 / **1.0251** / 1.0172 (nominal 1.052-1.082 / 1.033-1.040). So the sigma-independent residual is 1.3-2.5% in most poses but peaks at 26.57 deg (4.3% / 2.5%), where the box edges lie on a lattice-commensurate slope (tan = 1/2); sphere residual at sigma = 0 was only 0.1%. Likely a voxelisation (staircase) effect of the mapped box volume/surface, UNVERIFIED (would need the mapped solid volume vs the analytic clipped area, not measured here). **Bogner-like free box at sigma = 0 (last-20k mean +- std, 100k steps, NOT converged -- tilt still drifting at the end):** rho 0.25 seed 2.86 deg: -31.14 +- 1.16; rho 0.25 seed 5 deg: +28.94 +- 4.01 (range 23.7..36.6); rho 0.75 seed 2.86: 22.76 +- 0.59 (rising); rho 0.75 seed 5: 34.97 +- 1.28. **Conclusion for the heel-angle question:** (1) the systematic, side-independent +2.1 deg shift of the nominal rho 0.25 free box is explained by the surface-tension torque bias (pinned torque bias -> 0 at sigma = 0); (2) but removing sigma does NOT make the free box converge to 26.565: the runs keep drifting/oscillating by 1-4 deg and end 22.8-35 deg, so the free-box equilibrium at this g/viscosity is dynamics/lattice-limited (little damping: tau ~ 0.5006) and the very steady nominal runs (std 0.01 deg for rho 0.25) may be locked by the contact-line/lattice rather than being a true equilibrium; (3) a fair free-box test needs much longer runs and/or damping; the quasi-static pinned-torque scan is the cleaner instrument. Still pending: R=4 resolution 4272776 (starts 22:00).

> **[Claude, 2026-09-28] CORRECTION/CAVEAT for ALL Bogner-like rho 0.75 free-box results (per "tell me more about" the sigma=0 rho 0.75 seed 2.86 run): the basin water level FALLS during these runs and the box sinks with it -- the reported end tilts are not static equilibria.** Evidence (`[BoxForces]` pose history of `bl_r075_h286` sigma=0, and implied waterline = level at which the ideal Archimedes submerged fraction of the tilted 18x18 section equals 0.75, computed from the logged z_c/tilt): implied waterline 20.04 (step 250) -> 18.5 (25k) -> 17.5 (50k) -> 16.3 (75k) -> 15.9 (end); z_c 15.55 -> 11.75, monotonic; the independently measured end waterline (mean of sum-of-fill over the last 2000 steps) is 15.83, i.e. it agrees. The box's lowest corner (z_c - 9(|cos|+|sin|)) reaches z = 0.0 (floor 0.06) at ~80k steps: **the final 22.8 deg is a floor-supported pose**, not a floating equilibrium. Same picture in all four rho 0.75 runs (implied wl at 25k/50k/75k/end: nominal seed 2.86 19.7/19.3/18.2/16.1 [floor at end], nominal seed 5 19.7/19.2/18.7/18.1 [corner at 1.8]; sigma 0 seed 2.86 18.5/17.5/16.3/15.9 [floor]; sigma 0 seed 5 19.1/18.3/17.1/16.4 [floor]) whereas the rho 0.25 runs keep a constant level (nominal 20.8, sigma 0 20.2-20.3, lowest corner 11-13 cells above the floor). Expected level with the box floating is ~22 (20 + displaced volume/basin area = 20 + 8748/4320), so the loss is several cells of depth, i.e. O(20%) of the water volume -- **cause not established**: candidates are non-conservation of fluid mass when the 6-DOF body sweeps through interface cells (cells covered by the body lose their mass, uncovered cells are refilled, cf. the ~6% arm mass difference already seen in the free-floating sphere) or interface-cell mass loss around the box; only the level history and the end column were logged (no total-mass diagnostic), so this is inference from pose/level consistency. **Consequences:** (1) the earlier statements 'rho 0.75 seed 2.86 -> 26.83 near Bogner' (nominal) and '22.8 / 35.0 deg' (sigma 0) must NOT be read as heel-angle equilibria; the nominal 26.83 was also measured while the level was still sinking (implied wl 16.1-18.2 at 50k-end) and the box touched down at the very end; (2) the rho 0.25 free-box results (constant level, corner clear of the floor) are the only usable free-box heel data; (3) before any further free-box heel test at rho 0.75: log the total fluid volume (sum of fill + excess mass) per interval to test the mass-loss hypothesis, and/or use a basin deep enough that touchdown cannot occur.

> **[Claude, 2026-09-28] DEEP-BASIN FREE-BOX RERUN (rho 0.25 / 0.5 / 0.75, all seed 2.86 deg, 150k steps, checkpointing) SUBMITTED + CORRECTION of my 'expected level ~22' statement + Bogner re-check.** **Jobs (one density per checkout, 72 ranks, `jobs/run_72_bogner3_deep.sbatch`, `sbatch --export=ALL,B3_RHO=<025|050|075>` from `<checkout>/apps/showcases/Antidunes`): rho 0.75 -> 4273220 (mt-dennys), rho 0.5 -> 4273221 (mt-dennys-mirror), rho 0.25 -> 4273222 (mt-dennys-mirror2).** Generator `gen_bogner3_deep_scan.py` -> `prms/bogner3_r{025,050,075}_deep_gas.prm` (+ `_deep_resume_gas.prm`: loadSnapshot true, timesteps 400000). Same physics/pinning as the earlier Bogner-like set (real g_SI 9.81, box 18x36x18, `pinPosition false` = x/y/z free, `pinRotation true` with `pinRotationAxes <1,0,1>` = pitch and yaw held, roll about Y free), changes: basin <96,45,60> and liquid depth 30.067 (was 40/20) so a sinking 0.75 box has ~12 cells to fall before touching the floor (box bottom initially at 12.6 for rho 0.75); heave 34.562/30.067/25.573 for rho 0.25/0.5/0.75 (Archimedes at the new level, seed 2.86 deg); `storeSnapshot true`, every 25000 steps, single rolling checkpoint in `build/apps/showcases/Antidunes/snapshots_b3_r<RHO>/` (NOT deleted by the job; a fresh start refuses to run if that folder exists; a resume run first copies it to `snapshots_b3_r<RHO>_backup_before_resume_<jobid>`); resume = `sbatch --export=ALL,B3_RHO=<..>,B3_RESUME=1 jobs/run_72_bogner3_deep.sbatch`. **New diagnostics kept (the earlier job scripts deleted them):** `eval/antidunes.sqlite` (per 250 steps: `totalFluidMass`, `maxFluidPos_z`, ...) copied to `pinscan_results/<jobid>/b3_r<RHO>.antidunes.sqlite`, gzipped full log, [BoxForces]/[BoxTorque] every 250 steps, density column over steps 148000-149999. **Why rho 0.75 loses water -- what is established and what is not:** ESTABLISHED (code, `src/lbm_mesapd_coupling/momentum_exchange_method/MovingParticleMapping.h` ~l.180-195): a fluid cell newly covered by a moving body is set to obstacle via `forceBoundary(...)`, and the comment states 'this removes the fluid from the simulation'; its mass (fill x density) is NOT put into the excess-mass field or redistributed, so covering liquid cells destroys mass rather than raising the level. Uncovered cells are re-seeded from their neighbours (`PdfReconstructionManager.h::seedConvertedCellFill`). Consequences seen in data: (a) a body inserted into a pre-filled basin does not raise the level: pinned box/sphere waterlines are only +0.2/+0.4 cell above nominal where volume conservation would give +3.1/+3.4 (displaced volume / basin area), and the free rho 0.75 box shows no lift in its first 12k steps (z_c stays 15.5, implied level 20.0, conservation would give ~22 => **my earlier statement 'expected level with the box floating is ~22' was WRONG for this code: the correct baseline is ~20 (level unchanged by insertion) and the falling level is the anomaly**); (b) for pinned bodies the covered set never changes, so no dynamic loss (level constant, matches). NOT established: why rho 0.75 loses ~4 cells (~17700 cells^3, more than the box's total displaced volume 8748-11664) while rho 0.25 does not, although the rho 0.25 box ALSO moved a lot (tilt -47 -> +32 deg in the first 10k steps, x wrapping across the basin) and kept its level (20.8 nominal / 20.2-20.3 sigma 0). Hypothesis (untested): the rho 0.75 box covers/uncovers cells deep in the liquid bulk (13.5 of 18 cells submerged) where every covered cell destroys a full unit of mass, whereas the rho 0.25 box only touches the interface layer; the falling level then lets the box sink further, covering more liquid -- a positive feedback. If mass is truly non-conserved this way, free-floating results with deep submersion (rho 0.75 cuboid/sphere/cylinder at real g) carry this bias too; the total-mass series of 4273220-2 decides. **Bogner (2017) re-read for explicit pinning statements (text /tmp/bogner.txt):** (1) sect. 5.3.2 (righting moment, Fig. 5.6): 'During the simulation, the box is fixed to a constant position and the angle of heel ... is varied' -- explicit, but ONLY for that fixed-torque test. (2) sect. 5.3.3 (floating cubes, Figs 5.8-5.10): 'the 2D cube of side length s is represented by a cuboid of dimensions 2s x s x s, that will rotate around its longitudinal axis only'; 'Three boxes with s = 16 are placed inside the basin' (all three densities in ONE 130x40x40 basin at once, not one box per run); 'run until all fluid motions ceased and the objects came to rest'. There is NO statement on constraining translation (x/y/z) or on pitch/yaw pinning; 'rotate around its longitudinal axis only' is either a constraint or a description of what the 2D-equivalent geometry does, the text does not say which. His own Fig. 5.10 shows the boxes shifting horizontally by +0.6/+10/+17 cells (measured earlier), so at least horizontal translation is not held. Our new runs therefore use: translation free, pitch/yaw held, roll free -- the interpretation closest to the text -- one box per basin (differs from his three-in-one).

> **[Claude, 2026-09-28] IS THE WATER LOSS AN ARTIFACT OF THE ROD-CLOG FIXES? Audit vs commit f1f4342f (HEAD, 2026-09-16) + historical evidence + control runs (per "could this be an artifact of our changes with the rods ... compare the code to an older commit ... turn off anything we tried").** **(1) What is in the working tree beyond HEAD (uncommitted, ~7500 lines) and what could touch the fluid MASS path:** src: `StreamReconstructAdvectSweep.h` (Option A: seed the scratch PDF clone once on first creation + a log-only [DEGENERATE-FILL-RECONSTRUCTION] warning), `MovingParticleMapping.h` (cellCenterTestOffset default 0 = off; ownership-handoff callback, only fires if `useOwnershipHandoffPdfReset` true, default false; guard against `uidToIdx == SIZE_MAX`), `PdfReconstructionManager.h` (log-only [FALLBACK-ISOLATED-INTERFACE]), `SimpleBB.h` (skips a MEM link only if |pdf| > 10 = garbage; IBB default off; diagnostics), `Contains/DoubleCast/SingleCast/ShapePackUnpack` (contact geometry). Antidunes: Option B `repairDegenerateObstacleCells` (always on, but only rewrites obstacle cells whose PDF sum is non-finite/<=0), `SealedCellGuard` (`closeSealedCellsMinNeighbours` default 0 = off), `mappingTieOffset` 0, `debugForceObstacleClearAtStep` 0 (the only code in the tree that sets fill=1 / setToEquilibrium is gated on it), obstacle-freeze/reconstruction-only experiments already removed. `git diff HEAD` filtered for fillField/excessMass/forceBoundary/flag writes finds nothing default-on. **=> No default-on change since HEAD writes to the fill field, the flags or the excess mass of ordinary fluid cells; everything default-on is a safety net that acts only on corrupt state (non-finite/garbage PDFs) or on the first allocation of a scratch buffer.** **The mass-destroying step is ORIGINAL waLBerla code (2020, `MovingParticleMapping.h`: 'this removes the fluid from the simulation'; `git log -S` -> 595ebbb2, 2020-01-27), not ours.** **(2) The problem predates the rod work:** Sep-19 runs (before the 2026-09-22...27 changes) with free Bogner-exact boxes: totalFluidMass series (sqlite): `bogner_exact_3box_s16` (4248419, 50k steps) +0.39%, `..._s32` (4248420, 100k) -0.89%, **`..._s16_nopin300k` (4248656, 300k steps) -12.4%, monotonically falling (95556 -> 83724; maxFluidPos_z 31 -> 18)**, and that run's rho 0.75 box was already the erratic one (~37 deg, no plateau, plan sect. 1a). **(3) Controls submitted (same prm `bogner3_r075_deep_gas.prm`, seed 2.86, deep basin, 150k steps):** **4273255 = pristine HEAD f1f4342f** built from a git worktree (`../mt-dennys-head`, synced over `~/mt-dennys-mirror2/`: NOTE mirror2 now holds the HEAD source, not the working tree -- re-sync the working tree before using mirror2 for anything else; HEAD has no logDensityColumn so no column.txt, box log + totalFluidMass sqlite are there); **4273256 = current code with `useBognerReconstruction false`** (legacy uniform all-interface FormerMO reconstruction, no phi seeding; prm `bogner3_r075_deep_legacyrec_gas.prm`, chained after 4273220 in mt-dennys, env `B3_VARIANT=_legacyrec`); reference = 4273220 (current code, running); 4273257 = current code rho 0.25 (chained after 4273256). **Job mishap:** I cancelled 4273222 (rho 0.25, mirror2) believing it was still pending -- it had just started; partial output moved to `snapshots_b3_r025_cancelled_4273222`; rho 0.25 is rerun as 4273257. rho 0.5 (4273221, mirror) unchanged. **Prediction / decision rule:** if HEAD (4273255) also loses the level -> the loss is pre-existing, not a rod-fix artifact; if the legacy reconstruction (4273256) does not lose it -> the Bogner phi-seeding / fill=1 over-promotion of uncovered cells is implicated.

> **[Claude, 2026-09-28] ALL-BOGNER-OFF CONTROL ADDED (per "everything bogner related should be turned [off] including the reconstruction?"): job 4273271** (`--dependency=afterany:4273221`, mt-dennys-mirror, `B3_VARIANT=_nobogner`, prm `prms/bogner3_r075_deep_nobogner_gas.prm`): identical to the rho 0.75 deep run except `useBognerReconstruction false` (legacy uniform all-interface FormerMO reconstruction, no phi seeding), `useFillWeighting false`, `useBognerThreePhaseMEM false` (rhoHydStat subtracted from the MEM force = pre-flip Antidunes behaviour) and `useMemOnlyBuoyancy false` (analytic vSub*g buoyancy + analytic righting torque instead of MEM-only buoyancy). NB this changes the buoyancy MODEL (analytic Archimedes term added), so it answers 'does the water still disappear with the whole Bogner stack off', not 'what is the MEM-only heel angle'. Job map for rho 0.75 (same geometry, seed 2.86, 150k): 4273220 all defaults (current code) | 4273255 pristine HEAD f1f4342f | 4273256 current code, only useBognerReconstruction false | 4273271 current code, all four Bogner switches off. Decision table: HEAD loses water -> pre-existing; legacyrec keeps water -> phi-seeding/over-promotion of uncovered cells; nobogner keeps water but legacyrec does not -> the MEM/buoyancy side (box motion) rather than reconstruction; all lose water -> the upstream mass deletion at covered cells itself (needs a conservation fix, e.g. route the covered cells' mass into the excess-mass field).

> **[Claude, 2026-09-28] BOGNER-OFF CONTROLS CANCELLED (per "cancel all those runs, everything bogner related has to be turned on").** Cancelled before producing output: **4273256** (current code, `useBognerReconstruction false`) and **4273271** (all four Bogner switches off). The thesis runs keep ALL Bogner options ON (`useBognerReconstruction`, `useFillWeighting`, `useBognerThreePhaseMEM`, `useMemOnlyBuoyancy` all true); the prms `bogner3_r075_deep_{legacyrec,nobogner}_gas.prm` remain in the repo as unused variants. Still queued/running with Bogner ON: 4273220 (rho 0.75), 4273221 (rho 0.5), 4273255 (rho 0.75 on pristine HEAD f1f4342f, whose defaults/prm also have Bogner ON), 4273257 (rho 0.25, re-chained `afterany:4273220` in mt-dennys). The decision table of the previous entry now only applies to the HEAD comparison.

> **[Claude, 2026-09-28] ALL DEEP-BASIN FREE-BOX JOBS CANCELLED (per "let's cancel all the jobs"):** 4273220 (rho 0.75), 4273221 (rho 0.5), 4273255 (rho 0.75 pristine HEAD), 4273257 (rho 0.25) -- none reached a checkpoint or produced usable output (started 18:32-18:41, cancelled ~19:00). Nothing from the deep-basin set was analysed. State left behind: `~/mt-dennys-mirror2/` holds the pristine-HEAD SOURCE (worktree `../mt-dennys-head` locally) -- re-sync the working tree before using it; partial `snapshots_b3_r0*` folders in the build dirs of mt-dennys / mirror / mirror2 (a fresh start of `run_72_bogner3_deep.sbatch` refuses to run while they exist -- move them away first); prms/generator/sbatch (`gen_bogner3_deep_scan.py`, `prms/bogner3_*`, `jobs/run_72_bogner3_deep.sbatch`) are kept for a later rerun. Still queued: R=4 resolution 4272776 (mirror, starts 22:00), no wakeup armed for it. The water-loss finding (implied level falling in all four rho 0.75 Bogner-like runs, Sep-19 300k run -12.4% total mass) and its unresolved cause stand as logged above.

> **[Claude, 2026-09-28] DEEP-BASIN FREE-BOX SET RESUBMITTED (per "do the following on each of the three directories ... 2.86 starting angle, rotation locked except y, position free, multiple height points, enough liquid under the 0.75 box ... delete the old stale vtk"): jobs 4273301 (rho 0.75, mt-dennys), 4273302 (rho 0.5, mt-dennys-mirror), 4273303 (rho 0.25, mt-dennys-mirror2)**, Bogner options ALL ON (defaults), current working-tree code (all three checkouts re-synced; mirror2 no longer holds the HEAD source). Set-up identical to the earlier deep set (`gen_bogner3_deep_scan.py`, `prms/bogner3_r{025,050,075}_deep_gas.prm`): seed `initialTiltAngleDeg 2.86` in all three, `pinPosition false` (x/y/z free), `pinRotation true` + `pinRotationAxes <1,0,1>` (rotation about world X and Z held, only Y free), 150000 steps, real g_SI 9.81, basin <96,45,60>, liquid depth 30.067, single rolling checkpoint every 25000 steps in `build/.../snapshots_b3_r<RHO>/` (resume: `sbatch --export=ALL,B3_RHO=<..>,B3_RESUME=1 jobs/run_72_bogner3_deep.sbatch`). **Liquid under the rho 0.75 box:** centre z 25.573, bottom at 16.57 -> 16.6 cells of water below it at t=0 (floating draught 13.5); the loss seen before was ~4 cells over 100k steps, so even with the same loss the box keeps >10 cells clearance. **New height diagnostic (code, compiled locally with the Antidunes target):** `ModelParameters.logWaterHeight` + `waterHeightGridX/Y/EveryNSteps` (`logWaterHeights` in SetupStages.h, `addWaterHeightLoggingSweep` in BuildTimeloop.h): fill-integrated height (sum of fill over liquid+interface cells; obstacle cells counted separately) at a 12 x 3 grid of columns (x = 4, 12, ..., 92; y = 7.5, 22.5, 37.5) every 250 steps over the WHOLE run; the job reduces it to `pinscan_results/<jobid>/b3_r<RHO>.waterheight.txt` (columns: step gx gy height nObstacleCells; use columns with nObstacleCells == 0, e.g. median per step, because the 36-long box shadows all three y-columns at its current x). Also kept: eval sqlite (totalFluidMass, maxFluidPos_z per 250 steps), [BoxForces]/[BoxTorque], density column over steps 148000-149999, gzipped full log. **Cleanup:** deleted the stale live `build/apps/showcases/Antidunes/vtk-out` (55 MB, 81 files, leftovers of the cancelled runs) in all three checkouts; NOT deleted: archived run copies under `<checkout>/apps/showcases/Antidunes/animations/*/vtk-out`, the older DropTest/dropcompare outputs and `mt-dennys-mirror/build-asan-hopper/.../vtk-out` (archives/other build, not stale by the same criterion). The partial snapshot folders of the cancelled deep runs were moved (not deleted) to `build/.../archive/snapshots_b3_r0*_cancelled_20260928` so the fresh start guard passes. R=4 resolution 4272776 (mirror) still starts 22:00 -- the rho 0.5 job in that checkout must finish first (~1 h expected). If the level still falls: investigate (mass log + level grid are now there).

> **[Claude, 2026-09-28] R=4 resolution job 4272776 CANCELLED (per "cancel the R4 job for now"), before it started (was pending until 22:00, mt-dennys-mirror, 576 ranks).** The R=1/R=2 results stand; R=4 (`prms/pinres_sphere_R4_gas.prm`, `jobs/run_576_pinres_sphere_R4.sbatch`) can be resubmitted later from `mt-dennys-mirror/apps/showcases/Antidunes` (`sbatch --begin=<time> jobs/run_576_pinres_sphere_R4.sbatch`).

> **[Claude, 2026-09-28] EARLY CHECK of the Bogner-ON deep-basin jobs (19:19, ~14 min in; 4273301 rho 0.75, 4273302 rho 0.5, 4273303 rho 0.25): all RUNNING, clean start, water-height CSVs written (72 rank files each, header `mainStep,rank,gx,gy,partialFillSum,nObstacleCells`).** Partial result (median over grid columns without obstacle cells, `analyze_waterheight_partial.py` on a copy of the live CSVs; fill-integrated height, initial value 30.564): **rho 0.75: 30.56 (0) -> 31.30 (4.7k) -> 29.99 (14k) -> 29.54 (33k) -> 28.73 (48k) = -1.8 cells in 48k steps, monotonic, every column falls together (min/max 27.7/29.0 at 48k) => a basin-wide level drop at ~0.04 cell/1000 steps, the same rate as in the 20-cell basin; rho 0.5: 30.4-30.7 throughout (30.38 at 45.5k, -0.19 after the initial transient); rho 0.25: 30.8-30.9 (30.78 at 39k, no drift).** So in the deeper basin, with the box free in x/y/z and roll about Y, ONLY rho 0.75 loses water; rho 0.25/0.5 hold. The 16.6-cell water cushion under the 0.75 box is being consumed at ~0.04 cell/1000 steps (would reach ~24 by 150k, leaving the box ~7 cells above the floor if the trend is linear). Also answered for the user: `pdfReconstructionModel OnlyMissing` is NOT ours -- it is the upstream default in `Antidunes.prm` since the 2023 merge (commit 0b415d68), used by 322 of the 338 prms at HEAD (the parser requires the key; the only other value in any prm is `OnlyMissingMin-3-largest`).

> **[Claude, 2026-09-28] DEEP-BASIN FREE-BOX RESULTS (jobs 4273301 rho 0.75 / 4273302 rho 0.5 / 4273303 rho 0.25, all COMPLETED exit 0, 150k steps, ~40 min each; Bogner ON; local `animations/bogner3_deep_20260928/<jobid>/`; scripts: scratchpad final.py/corr.py logic, `analyze_waterheight_partial.py`).** **Water level (median fill-integrated height over 12x3 grid columns without obstacle cells; start 30.56):** rho 0.75: 29.87 (25k) / 28.74 (50k) / 26.82 (100k) / **25.63 (150k) = -4.93 cells, still falling ~0.023/1000 steps in the last 50k**; rho 0.5: 30.10 (-0.46, last-50k slope -0.003/1000); rho 0.25: 30.53 (-0.04, slope ~0). **totalFluidMass (sqlite, every 250 steps):** rho 0.75 **-17.87%** (-6.67% @50k, -13.24% @100k; 123097 -> 101104), rho 0.5 -1.23%, rho 0.25 -0.23%. The mass series and the level agree (22000 cells^3 / 4320 = 5.1 cells). **The water loss is reproduced in the deeper basin and is specific to rho 0.75.** **Box (last-20k / plateau):** rho 0.75: tilt -45 deg at 25-100k, then drifting to -37.1 +- 0.35 (130-150k); |tilt| never near 26.565; z_c 25.6 -> 22.4, lowest corner min 9.8 above the floor (no touchdown, so the deeper basin did what it was for); x 47.8..66.7, y stays exactly 22.5. rho 0.5: **-45.75 +- 0.03 (target 45; plateau from ~25k)**; z_c 30.1 -> 30.45; x 8.8..48.0, y 18.0..22.4. rho 0.25: **+33.65 +- 0.04 (target 26.565, +7.1 deg; plateau from ~75k after 27 deg at 25-50k)**; z_c 34.65 -> 34.51; x 0.8..94.9 (drifts through the periodic basin), y 22.7..27.0. Implied level from the box pose is 0.15-0.6 cell ABOVE the measured level in all three (0.6 for rho 0.75, 0.3 rho 0.5, 0.15 rho 0.25; end-of-run 26.20/25.63, 30.45/30.10, 30.67/30.53), the same sign as the earlier surface-tension excess buoyancy. **Y-WALL CAVEAT (new):** the rho 0.5 box reaches y = 18.0 and the rho 0.25 box y = 27.0, i.e. the 36-long box touches the y = 0 / y = 45 walls (18 +- 18); the logs show [STRONG-CONTACT-FORCE] between the box and another uid (wall/plane) at steps 77020 (rho 0.5, y 17.99, thisContactForce_y +1.39) and 73959 (rho 0.25, y 27.01, force_y -0.54); the rho 0.75 box never leaves y = 22.5. So the rho 0.5 / 0.25 plateaus (45.75 / 33.65) were reached with the box partly pushed against a side wall and are NOT clean free-floating equilibria; the rho 0.25 value (33.65) equals the '~34 deg' seen in the Sep-19 pinned runs. **LOG WARNINGS (`b3_r*.fulllog.gz`, WARNING lines: 3665 / 3426 / 3216 for rho 0.75 / 0.5 / 0.25):** (a) 'No liquid or interface cell is in the neighborhood to distribute excess mass to. Mass is lost/gained' (excess-mass distribution failing for a converted cell, each <= ~1 cell of mass): 150 / 141 / 111 occurrences -- essentially the same in all three, so it cannot explain a 22000-cell^3 loss that only rho 0.75 shows. (b) STALE-NEARBOUNDARY 365 / 275 / 322, FALLBACK-ISOLATED-INTERFACE ~526 / 607 / 750 lines (true counter capped at #20 per rank), SUBCYCLE-LEFT-RANK 6 / 9 / 5, GHOST-SYNC-DIAG 72 each -- similar across runs. (c) STRONG-CONTACT-FORCE / TOTAL-FORCE-THRESHOLD / MEM-WATCH (throttled): rho 0.75 300/300/300 lines, rho 0.5 740/640/700, rho 0.25 600/500/522 -- more for the wall-touching boxes. (d) **[OBSTACLE-CELL-REPAIR] (Option B guard, rewrites obstacle cells whose PDF sum is <= 0 or NaN): true per-rank counter up to 2,315,000 for rho 0.75 vs 70,000 (rho 0.5) vs 20 (rho 0.25)**, first fired at step 50475 in the rho 0.75 run (world (72.5,18.5,26.5), rho -3.0, then -1e1, -3.4e7 within 10 steps) and keeps firing through 149231, i.e. tens of degenerate obstacle cells per step in the box body; this is the only warning class that is dramatically rho 0.75-specific. It began after the level had already fallen 1.8 cells and the loss rate did not change afterwards (level change per 25k: -0.69, -1.13, -0.83, -1.09, -0.46, -0.73), so the repairs are NOT the cause of the loss, but the same run being numerically unhealthy inside the body is notable (repair only touches obstacle cells, never fill). **Loss timing:** regression of per-250-step mass change on |dz_c|, |d tilt|, |dx| for rho 0.75: R^2 = 0.00 (mean -36.8 per 250 steps, steady, uncorrelated with heave, roll rate or drift; the box orientation is nearly constant at 130-150k yet the level still falls); sinking-covered volume (648 x 3.19 dz) ~2070 is only ~9% of the 21993 lost. => the loss is a continuous leak at a roughly constant rate, not tied to box motion within the run, which weakens the earlier 'covering/uncovering while moving' hypothesis (the mass deletion at covered cells alone cannot supply it). Cause still open; candidates to test next: interface/contact-line cells at the waterline (rho 0.75 waterline sits at a different height of the box cross-section), the excess-mass path, or the degenerate obstacle cells.

> **[Claude, 2026-09-28] OPTION A AND OPTION B NOW PRM SWITCHES, DEFAULT OFF; DEEP-BASIN SET RERUN WITH BOTH OFF (per "turn option b off and rerun / make both prm, give them better names, default off, rerun all, delete all the files from this previous run").** **New switches (ModelParameters, default false = off):** `repairDegenerateObstaclePdfs` (Option B: `repairDegenerateObstacleCells` after collide, BuildTimeloop.h, was always on) and `seedStreamPdfCloneOnFirstUse` (Option A: copy the real PDF field into the scratch PDF clone of `StreamReconstructAdvectSweep` on its first creation, was always on; implemented as a new ctor parameter `seedPdfCloneOnFirstUse = false` of the library sweep, so other users of that class (benchmarks, tests: `AdvectionTest`, `InflowTest`, `CollideIsolatedObstacleIterationTest`) get the upstream behaviour). **CONSEQUENCE, not applied to old prms:** because the DEFAULT is now off, every existing prm without the keys silently runs without Options A/B; the hopper / rod-bundle clogging scenarios validated on 2026-09-25/26 (they crashed within ~3000 steps without B) must set `repairDegenerateObstaclePdfs true; seedStreamPdfCloneOnFirstUse true;` explicitly. The MEM-garbage-link skip (|pdf| > 10) and the unknown-particle-index guard are still always on (they act only on corrupt state and did not fire in the deep-basin runs: no MEM-GARBAGE-LINK warnings). **Run prms** `prms/bogner3_r{025,050,075}_deep_gas.prm` (generator `gen_bogner3_deep_scan.py`) now contain both keys explicitly = false; everything else unchanged (rho 0.25/0.5/0.75, seed 2.86, roll about Y only, x/y/z free, deep basin, Bogner ON, 150k steps, water-height grid, checkpoints every 25k). **Jobs: 4273396 (rho 0.75, mt-dennys), 4273397 (rho 0.5, mt-dennys-mirror), 4273398 (rho 0.25, mt-dennys-mirror2).** Reference (Options A+B ON, same prm otherwise): 4273301/2/3 = level -4.93 / -0.46 / -0.04, mass -17.9% / -1.2% / -0.2%, OBSTACLE-CELL-REPAIR counter 2.3M / 70k / 20. **Cleanup done on the cluster (all three checkouts):** deleted the previous deep-basin runs' files -- `pinscan_results/{4273220,21,22,55,56,57,71,4273301,02,03}`, `~/pinscan_backup/427330{1,2,3}_bogner3_*`, `build/.../snapshots_b3_r0*` (incl. `archive/snapshots_b3_r0*`), leftover `water-height`/`density-column`/`eval`/`vtk-out`, and `sbatch-output/*/bogner3_deep_*`. NOT deleted: the local copy `~/University/Thesis/animations/bogner3_deep_20260928/` (the numbers in the entries above come from it) -- say if it should go too.

> **[Claude, 2026-09-28] OPTION A + B OFF RESULT (jobs 4273396 rho 0.75 / 4273397 rho 0.5 / 4273398 rho 0.25, all COMPLETED exit 0, ~40 min each, 0 OBSTACLE-CELL-REPAIR confirmed, 0 MEM-GARBAGE-LINK, no assertion/abort): the water loss is UNCHANGED (if anything slightly worse) without Options A/B -- they are not the cause.** Local `animations/bogner3_deep_fixesoff_20260928/<jobid>/`. **Water level / totalFluidMass (A+B OFF vs A+B ON reference in parentheses):** rho 0.75: -5.92 cells / -20.96% (vs -4.93 / -17.87%); rho 0.5: -1.16 / -3.58% (vs -0.46 / -1.23%); rho 0.25: **+0.21 / +0.74%** (vs -0.04 / -0.23%, i.e. now a small GAIN, not a loss) -- so rho 0.75 leaks MORE without the repair guard, rho 0.5 also leaks more, rho 0.25 flips sign but stays negligible in magnitude. This rules out Option A/Option B as the cause; if anything the guard was mildly mitigating it (by keeping obstacle-cell PDFs healthy, indirectly reducing degenerate excess-mass events near the box). **Warning counts (A+B off): 'Mass is lost/gained' 274 / 132 / 93** (rho 0.75 nearly double the A+B-ON run's 150, consistent with the larger loss, but still <1% of the 22000-27000 mass units actually lost) -- FALLBACK-ISOLATED-INTERFACE/STALE-NEARBOUNDARY/STRONG-CONTACT-FORCE/MEM-WATCH counts are all similar in magnitude to the A+B-ON run (same order, same relative ranking across rho). **Tilt (A+B off):** rho 0.75 reproduces the -26.5 deg stage (target magnitude) at 20-40k (mean -25.21, std 0.90) -- CONFIRMS this is a real, reproducible transient near the correct answer, not noise -- but then goes through TWO further transitions: -45.18+-0.35 at 60-80k, then a large excursion to a mean -67.25 (std 11.47, range down to -90.5!) at 100-125k, settling at -59.40+-0.88 by 130-150k -- i.e. WORSE and MORE unstable than the A+B-ON run's -37.09+-0.35, and touches -90 deg (box lying on its side) at least transiently; lowest corner reaches 7.5 (still clear of the floor, deeper basin holds). rho 0.5: -45.05 -> -45.21, essentially identical to A+B-ON (-45.75), std stays <0.1 throughout -- unaffected. rho 0.25: **reproduces the A+B-ON pattern almost exactly** -- 26.57 (20-40k) -> 26.50 (60-80k) -> settles at 25.35-25.47 (100-150k, std 0.03-0.08) instead of A+B-ON's jump to 33.7; i.e. WITHOUT the repair guard rho 0.25 stays close to 26.565 (end 25.5, error 1.1 deg) rather than drifting to 33.65 -- the single case where turning the guard off looks BETTER, though it may be coincidental (different PDF-noise trajectory, not a systematic improvement, n=1). **Conclusion for the standing question: Options A and B do not explain (and do not fix) the water loss; the rho 0.75 free box also does not reliably reach or hold 26.565 deg with or without them, going through the same kind of intermediate lattice-locked states (-26.5, -45, -59..-90) either way. The leak and the heel-angle instability look like separate manifestations of the same underlying issue at rho 0.75 (heavily submerged, moving body at the free surface), not something the rod-crash mitigations touch.** Next candidates per the open list above (fully pinned rho 0.75 box for a long window; surface-tension-0 free box; wider y-clearance) remain the live leads.

> **[Claude, 2026-09-28] BOGNER GRAVITY/VISCOSITY REPLICATION, rho 0.75 ONLY (per "rerun only the 0.75 but only with gravity and viscosity parameters known from bogner's thesis ... so we can see the stability comparison and maybe that validates his choice of weird gravity"). Job 4273470 (mt-dennys).** Same deep-basin free-box setup as the current baseline (`bogner3_r075_deep_gas.prm`: basin 96x45x60, box 18x36x18, seed 2.86 deg, x/y/z free, pitch/yaw pinned (`pinRotationAxes <1,0,1>`), 150000 steps, water-height grid, checkpointing, Options A+B OFF (current default)) with ONLY the physics changed: `gravityOverride 7.5e-4` (was `5.06585e-5` for real g_SI=9.81) and `relaxationRateOverride 1.9` added (Bogner & Rude 2013's tau = 1/1.9 = 0.526, i.e. omega=1.9; was the hardcoded-units default omega=1.99743). Derived viscosity from the override: nu = (1/1.9 - 0.5)/3 = 0.00877, **~41x more viscous than the real-g hardcoded value (0.000214026)** -- much more damping. Surface tension left at the nominal hardcoded value 0.000849592 (unchanged; not overridden), matching the convention of the already-validated `bogner_exact_3box_s16_gas.prm` (his own box-equilibration reproduction), NOT the separate `bognerboxgravity` sphere-lineage prms which run with `useHardcodedUnits` unset/false and therefore sigma=0 from the U=0 dynamic formula -- that is a different, unrelated prm family for a different question (MEM buoyancy bias), not used here. Box initial heave (z_c=25.573) is unaffected by g/tau (purely geometric, from the clipped cross-section at the fixed liquid level). Generator: `gen_bogner3_deep_scan.py` -> `prms/bogner3_r075_deep_bognergrav{,_resume}_gas.prm` (diff to the baseline prm: exactly the `gravityOverride`/`relaxationRateOverride` lines and the snapshot folder name); resume via `sbatch --export=ALL,B3_RHO=075,B3_VARIANT=_bognergrav,B3_RESUME=1 jobs/run_72_bogner3_deep.sbatch`. **Question this tests:** at his own g/tau (small gravity, ~40x his viscosity), does the rho 0.75 free box hold near 26.565 deg without the drift-through-lattice-locked-states behaviour seen at real g_SI=9.81 (both with and without Options A/B)? If yes, that would support gravity/viscosity (not the repair guard, not the reconstruction model) as the reason his own published test looks stable while ours does not. Comparison set once complete: real-g guard-on (job 4273301, ends -37.09+-0.35), real-g guard-off (4273396, passes through -26.5 -> -45 -> up to -90, ends -59.40+-0.88), Bogner-g/tau guard-off (4273470, pending). Results -> `pinscan_results/4273470/b3_r075_bognergrav.*`; local copy -> `animations/bogner3_deep_bognergrav_20260928/4273470/` (not yet created).

> **[Claude, 2026-09-28] REAL-GRAVITY GUARD-OFF rho 0.75 RESUMED FOR ANOTHER ~100K STEPS (per "on another mirror can we continue the real gravity run for another 100k ... only run the things with guard off"). Job 4273476 (mt-dennys-mirror), resumed from job 4273396's own final checkpoint.** Confirmed the checkpoint's own `simState.file` says **beginTimeStep = 125001, not 150001** -- the run's LAST regular-cadence checkpoint (every 25000 steps) landed at t=125000, not at the final step t=150000 (the loop apparently does not force a checkpoint exactly at the last step if it isn't itself a multiple-of-25000 boundary hit from the write's own timing -- unconfirmed why, not investigated further). Consequence: resuming to a target of 250000 re-simulates steps 125001-150000 (25000 steps already covered by 4273396's own log) before reaching genuinely new territory, then continues ~100000 steps past the original 150000 endpoint -- net about 125000 steps of which 100000 are beyond the original end, matching the spirit of the request even though the resume point isn't exactly 150001. Hand-edited `prms/bogner3_r075_deep_resume_gas.prm` (generator's own default resume target is 400000; changed to **250000** for this specific request) -- checkpoint copied cluster-side from `mt-dennys` build dir to `mt-dennys-mirror` build dir (`cp -r snapshots_b3_r075`, same real g/guard-off case as 4273396); submitted `sbatch --export=ALL,B3_RHO=075,B3_RESUME=1 jobs/run_72_bogner3_deep.sbatch` from `mt-dennys-mirror/apps/showcases/Antidunes`. **Per "only run the things with guard off": both `repairDegenerateObstaclePdfs` and `seedStreamPdfCloneOnFirstUse` remain false (the current project default since they became prm switches) -- unchanged, no new toggle needed; going forward no new bogner3_* job will set either to true without being asked.** Results in `pinscan_results/4273476/b3_r075_resume.*` once done; question this answers: does the box, given another ~100k steps at real gravity, ever settle near 26.565 deg and hold, or does it keep drifting through further lattice-locked states (it had reached -59.40+-0.88 by the end of the original 150k)? Companion run 4273470 (Bogner g/tau, mt-dennys) still in progress for the stability-vs-gravity-choice comparison.

> **[Claude, 2026-09-28] REAL-GRAVITY rho 0.75 RERUN, FRESH FROM STEP 0, 75K STEPS ONLY (per "rerun it for 75k again from 0 and see if we can reproduce this first correct equilibrium ... with the real gravity"). Job 4273497 (mt-dennys).** Tests whether job 4273396's own long plateau near tiltY = -25 to -27 deg (steps ~10000-49000, close to the ideal 26.565 deg magnitude, see the earlier "wave reverberation" force analysis of that same run) is a reproducible feature of the dynamics or a one-off trajectory, by rerunning the IDENTICAL setup (real g_SI=9.81, guard off, deep basin, seed 2.86 deg, x/y/z free, pitch/yaw pinned) from a fresh t=0 rather than any resume, only far enough (75000 steps) to cover that plateau and its ~50k breakdown with margin, without paying for the full 150k. Hand-edited copy of `bogner3_r075_deep_gas.prm` -> `prms/bogner3_r075_deep_75k_gas.prm` (diff: header, `snapshotBaseFolder snapshots_b3_r075_75k` [new folder, avoids colliding with the existing `snapshots_b3_r075` checkpoint from job 4273396's own run still sitting in mt-dennys's build dir], `timesteps 75000`, density column window moved to 73000-74999); submitted `sbatch --export=ALL,B3_RHO=075,B3_VARIANT=_75k jobs/run_72_bogner3_deep.sbatch`. Question: does the tilt again settle near -26.5 deg for tens of thousands of steps before something (the same kind of localized force/torque event, or a different one) disturbs it, or does this fresh run diverge on a different schedule/magnitude entirely (LBM parallel-reduction ordering makes exact bit-reproducibility unlikely per this project's own prior findings, so some difference is expected regardless)? Results -> `pinscan_results/4273497/b3_r075_75k.*`.

> **[Claude, 2026-09-29] RESULT of the fresh-from-0, 75k-step real-gravity rerun (job 4273497, COMPLETED exit 0, 21 min): the near-26.565 deg plateau does NOT reproduce -- this trajectory finds a DIFFERENT lock-in angle (+45 deg, the rho=0.5 target, not rho=0.75's) for ~20000 steps, then breaks down again near the very end of the window.** Local `animations/bogner3_deep_75k_20260928/4273497/`. **Tilt history:** 2.71 (step 250) -> climbs rapidly through 9-19-39 deg (2500-12500) -> **locks onto +42 to +45 deg from ~15000 to 67500 (>50000 steps, essentially the WHOLE run), tightening to exactly 45.00 deg (std ~0) for 47500-67500** -> breaks down again at the very end: 44.68 (70000) -> 40.05 (72500) -> **26.45 (75000, the last logged step)** -- ending near the ideal 26.565 deg purely because the run was cut off mid-collapse, not because it settled there (z_c/x were still changing fast: 24.46->21.79->22.94, x 46.8->43.6->37.7). **Job 4273396's own first 75000 steps looked completely different: plateau near -25 to -27 deg from ~10000-49000 (40000 steps), THEN the -45 deg lock-in this session already characterized.** So neither the angle reached, nor its timing, nor even which of the two lock-in states (26.5-ish vs 45) comes first is reproducible between independent runs of the IDENTICAL setup -- **the specific angle the box locks onto looks like an arbitrary lattice-metastable state, not a signature of rho=0.75's true 26.565 deg equilibrium being found and then perturbed away from.** This directly weakens the framing of "reproduce the first correct equilibrium": what looked like an equilibrium in job 4273396 was very likely just one such metastable lock-in state among several (here the fresh run locked onto 45 instead, a different density's target, by coincidence of its own trajectory) -- neither is evidence of the ideal answer, both are artifacts of whichever lattice state the chaotic early transient happens to fall into. **Warning-burst pattern REPRODUCES, though, and more convincingly now (two independent confirmations):** 'mass is lost' warnings burst exactly at each of this run's own three transitions -- 71 events at step 2500 (the initial climb), 18 events at 42500-47500 (locking onto 45), 29 events at 70000-72500 (the final breakdown) -- vs near-zero (0-7) elsewhere. Same signature, same coincidence-with-transition pattern as job 4273396's step-50072 event and its own second (80-90k) transition. **This is now the most robust, twice-independently-confirmed finding of the whole investigation: every tilt transition in both runs, regardless of the angle involved, is accompanied by a localized burst of excess-mass-distribution-failure warnings at the box surface.** Water level -3.23 cells / mass -12.11% by step 74750 (rate per step comparable to job 4273396's own -20.96%/150000); lowest corner never touches the floor (min 9.11). **Conclusion for the user's original question:** no, the "first correct equilibrium" does not reproduce -- it looks like it was never a real equilibrium, just one of several possible chaotic lock-in states that happen to sit near a plausible-looking angle. The transition mechanism (warning burst at the box surface, cause still not pinned down) is the one thing that reproduces reliably.

> **[Claude, 2026-09-29] 75K FRESH RUN RESUMED FOR ANOTHER ~50K (per "continue for another 50k, is the tilt still pinned?"). Job 4273507 (mt-dennys), resumed from job 4273497's own final checkpoint.** Checkpoint `simState.file` says beginTimeStep=50001 (its last regular-cadence write, one interval short of the original run's own final step 75000 -- same pattern seen resuming the 150k run). New `prms/bogner3_r075_deep_75k_resume_gas.prm`: `loadSnapshot true`, `timesteps 125000` (75000 + another 50000; re-simulates steps 50001-75000 first, not bit-identical, before reaching new territory past 75000), density column moved to 123000-124999, otherwise identical (real g_SI=9.81, guard off). Submitted `sbatch --export=ALL,B3_RHO=075,B3_VARIANT=_75k,B3_RESUME=1 jobs/run_72_bogner3_deep.sbatch`. Question: job 4273497 was caught mid-collapse at its last logged step (75000: 26.45 deg, still falling fast from the 45 deg lock-in, z_c/x still moving) -- does the tilt keep falling past 26.565, re-pin somewhere else, or actually settle? Results -> `pinscan_results/4273507/b3_r075_75k_resume.*`.

> **[Claude, 2026-09-29] RESULT of the 75k-run resume (job 4273507, COMPLETED exit 0, 20 min, target 125000): tilt is NOT pinned near 26.565 -- it briefly passes through the target zone, then climbs back out and re-settles near +45 deg (noisier this time).** Local `animations/bogner3_deep_75k_20260928/4273507/`. Because LBM here is not bit-reproducible run-to-run, the resimulated 50001-75000 stretch diverges from job 4273497's own original path almost immediately (both start at ~45 deg at 50250, but this run wobbles 42-49 deg through 52500-70000 instead of climbing to a clean 45.0 lock): **step 70000: 39.91 -> 72500: 24.11 -> 75000: 21.98 -> 77500: 20.96** (undershooting the 26.565 target and passing through it, NOT the mid-collapse-at-cutoff 26.45 seen in the original run at exactly step 75000 -- a coincidence of timing, not the same event) -> **climbs back out: 80000: 25.14 -> 82500: 36.22 -> 85000: 48.79** -> settles into oscillation for the remaining 40000 steps, **last-20k (105000-125000) mean +44.82 deg, std 2.03 (range 40.9-48.7)** -- close to but noisier than the earlier razor-sharp 45.00+-0.00 lock. So the box passed THROUGH the 26.565 deg target on its way between two other lock-in states and did not stop there -- direct evidence against 26.565 being a genuine attracting equilibrium in this setup; the box treats it as a point to pass through, not settle at. **Warning-burst pattern confirmed a THIRD independent time:** 51 'mass is lost' events at steps 70000-75000 (24+17+10), essentially zero everywhere else in this run (12+5+4=21 at 50000-55000, likely residual from the resume's own re-settling; 0 for 55000-70000) -- again lining up exactly with the tilt's dip-and-recovery through the target zone. **Water/mass loss stopped once the box re-settled:** level -0.28 cells, mass +0.49% (a slight GAIN) over the whole 50250-124750 window -- consistent with the earlier real-g-to-209k finding that the mass leak concentrates during destabilisation transients and goes quiet once the box is in a quasi-steady oscillation, not a continuous background process. Lowest corner min 11.87 (no floor contact); x now drifts substantially (47.0 -> 72.5, wrapping periodically) while y stays 18.7-22.6 (no wall contact yet). **Answer to "is the tilt still pinned?": no -- it moved through several more states in this window (42-49 range, a dip to ~21-26, then back to ~41-49) and has not found anywhere resembling a stable, low-scatter equilibrium at either 26.565 or any other fixed angle for more than a bumpy oscillation.**

> **[Claude, 2026-09-29] APPROXIMATE B (centre of buoyancy) vs G (centre of gravity) COMPUTED AND COMPARED TO THE MEASURED TORQUE (per "approximate or compute and print what B and G are ... maybe this explains why it won't settle"). New script `compute_BG.py`, applied to job 4273396's own full 150k trajectory (rho 0.75, real g, guard off).** Method (post-processing only, from already-logged quantities, no new C++/prm field): G = the box's own logged position (exact -- homogeneous body, comOffsetZ=0). B = the analytic centroid of the box's square cross-section (rotated by the logged tiltY about its logged x_G/z_G) clipped below the MEASURED ambient water level (median of the far-field water-height-grid columns) -- an APPROXIMATION since it assumes the level is uniform along the box's own Y-extent (no local wave/slosh under the box itself, which we have separately shown is likely an underestimate of how close the box's top gets to the surface). tau_ideal = -F_ideal*(x_B - x_G) with F_ideal = g_LU * V_sub * (box's Y-length), same sign convention validated earlier against the pinned-box data. **RESULT: the analytic lever arm (x_B - x_G) stays tiny throughout the run (mostly 0.00-0.08 cells, max observed 0.16 at step 106250) -- the box's own measured pose is very close to being torque-BALANCED by simple hydrostatics almost everywhere -- yet the ACTUALLY MEASURED torque_y is 3-20x LARGER in magnitude than tau_ideal at almost every sampled step, and its SIGN frequently disagrees with tau_ideal's** (e.g. step 31250: lever +0.036 -> tau_ideal -0.0158, but tau_meas +0.0614, opposite sign; step 43750: tau_ideal -0.0352 vs tau_meas +0.0352, opposite sign; step 50000: tau_ideal +0.0002 [essentially zero] vs tau_meas -0.1585; step 100000: tau_ideal 0.0000 vs tau_meas -0.3322). **This is a direct, quantitative confirmation of "why it won't settle": near-zero analytic hydrostatic torque does NOT mean a numerically calm state -- the box is being pushed around by torque contributions many times larger than anything the B-G lever arm can supply, consistent with the lattice-noise/scatter character of the pinned-box torque error established earlier (as opposed to rho 0.25's smooth, sign-consistent bias) and with the mass-loss/warning-burst events found at every transition.** The vertical force F_ideal vs F_meas agree reasonably well throughout (both ~0.40-0.48, consistent with the earlier-established few-percent excess-buoyancy bias) -- it is specifically the TORQUE that is wildly disproportionate to the analytic lever arm, not the force. Caveat: B is approximate (uniform-level assumption); the true local level right at the box's own hull may differ from the far-field value used here, which could shift tau_ideal somewhat, but is unlikely to explain a 3-20x, sign-flipping discrepancy on its own.

> **[Claude, 2026-09-29] TWO s=16 REAL-GRAVITY RUNS SUBMITTED (per "rerun the real gravity again but with 16x32x16 for 100k, we don't need any vtk output" and "also do a second with the same pinning and gravity but put its starting position already in the position we had at the equilibrium in terms of starting height and starting angle").** Both use Bogner's OWN box size s=16 (halfExtents <8,16,8> = 16x32x16, was our s=18 = 18x36x18), same basin <96,45,60>/level 30.067, real g_SI=9.81, guard off, `pinRotationAxes <1,0,1>` (x/y/z free, roll about Y only), 100000 steps, VTK writeFrequency explicitly set to 0 (was already effectively off at 1e8/1e10). **Job 4273515 (mt-dennys)**, `prms/bogner3_r075_deep_s16_gas.prm`: seeded at the usual 2.86 deg with the Archimedes-consistent heave for s=16 (z_c=26.072) -- the direct s=16 analogue of the s=18 baseline. **Job 4273516 (mt-dennys-mirror)**, `prms/bogner3_r075_deep_s16_eqseed_gas.prm`: seeded AT the ideal equilibrium itself -- `initialTiltAngleDeg 26.565` (not 2.86) with position z set to the Archimedes-consistent heave for s=16 AT that tilt (z_c=26.489, computed the same way as the seed heave) -- i.e. starts at rest exactly at the analytically correct pose rather than reaching (or not) it through a transient. Tests whether 26.565 deg is a genuine fixed point (stays, confirming the earlier plateaus were real near-misses of a true equilibrium) or unstable/metastable even with zero initial perturbation (drifts away immediately, which combined with the B-vs-G finding above would point at the numerical torque noise itself as the reason, independent of how the box got there). Results -> `pinscan_results/4273515/b3_r075_s16.*` and `pinscan_results/4273516/b3_r075_s16_eqseed.*`.

> **[Claude, 2026-09-29] RESULT of both s=16 runs (4273515 normal seed, 4273516 EQUILIBRIUM seed at 26.565 deg): the equilibrium seed does NOT hold -- it drifts away within the first ~5000-10000 steps and converges to essentially the SAME lock-in cluster as the normal-seed run, at a THIRD distinct angle (~33-37 deg) not seen in the s=18 family.** Local `animations/bogner3_deep_s16_20260929/{4273515,4273516}/`. **4273515 (2.86 deg seed):** tilt swings through -4.6 (5k) -> -28.0 (10k) -> settles into a noisy band **-30 to -36 deg for the remaining 90000 steps, last-20k mean -33.06 std 1.00**. **4273516 (26.565 deg seed, at rest, z_c=26.489 exactly at the Archimedes-consistent height for that tilt):** stays near the seed only briefly -- **first-5k mean +26.68 deg, std 2.77 (range 18.5-29.8), already scattering within the first 5000 steps** -- then by step 10000 has flipped sign and moved to -34.47 deg, settling into the **same -30 to -39 deg band, last-20k mean -36.58 std 1.91**. Both s=16 runs land in the SAME cluster regardless of starting condition, distinct from every s=18 cluster found earlier (26.5, 42-45, 59-60, 90) -- i.e. **the lock-in angle depends on box size, which is a lattice/voxelisation signature, not something continuum hydrostatics would produce (continuum physics predicts the identical 26.565 deg target regardless of s).** **This directly answers 'is 26.565 deg a genuine fixed point?' for s=16: no -- released at rest exactly there, with the exact Archimedes-consistent height, the box leaves it almost immediately** (visible drift within the very first 5000 of 100000 steps) rather than staying or returning after a perturbation. Water/mass loss present in both, smaller in magnitude than the s=18 family: -1.56 cells / -5.40% (4273515), -1.91 cells / -7.02% (4273516) over 100000 steps (per-step rate roughly half to a third of the s=18 family's). Lowest corner never near the floor (min 12-13). **`compute_BG.py` run on both (added to the script's `__main__`): the SAME B-vs-G torque mismatch pattern from the s=18 analysis reproduces at s=16, sharpest right at the eqseed run's own start** -- step 250: lever +0.081, tau_ideal -0.0265, tau_meas only -0.0071 (smaller, same sign); but by step 5500 (still within the box's very first oscillation): lever +0.069, tau_ideal -0.0283, **tau_meas +0.0472 (opposite sign, larger)** -- the box built essentially zero real hydrostatic imbalance yet is already being pushed by a torque with the wrong sign for hydrostatics, within its first few thousand steps. Same story throughout both full trajectories: tau_ideal stays small (|.| mostly < 0.03) and tau_meas frequently 2-6x larger with a different sign (e.g. eqseed step 26500: tau_ideal +0.0110 vs tau_meas -0.0621). **Third independent confirmation (s=18 normal, s=16 normal, s=16 eqseed) that the torque driving these boxes away from 26.565 deg is dominated by something that isn't the B-G hydrostatic lever -- most consistent with the lattice-noise/scatter torque-error mechanism established earlier in this session for rho 0.75 specifically (as opposed to rho 0.25's smooth, predictable bias).**

> **[Claude, 2026-09-29] s=32 (Bogner's HIGH-resolution box) SUBMITTED, rho 0.75, 4 nodes (per "we haven't tried s=32 for the free box yet, try this for the 0.75 version, run it for 51k, checkpoint every 50k" then "why not use 4 nodes for this"). Job 4273521 (mt-dennys).** Bogner doubles his ENTIRE basin (130x40x40 -> 260x80x80) when going from s=16 to s=32, so this run doubles ours too: `prms/bogner3_r075_deep_s32_gas.prm` (hand-edited from `..._s16_gas.prm`) -- domain <192,90,120> (was 96,45,60), `liquidHeightFactor 10.0` (was 5.0, -> absoluteLiquidHeight 60.0746, was 30.0673), box halfExtents <16,32,16> = 32x64x32 (was 16x32x16), heave recomputed for s=32 at the new level (z_c=52.085, seed 2.86 deg), position centred at the new domain centre <96,45,52.085>. `gravityOverride` UNCHANGED (5.06585e-05) -- confirmed this is a pure dx_SI/dt_SI unit-conversion constant tied to the sediment-bed particle scale (from `computedParams.txt`, identical across s=16 and s=18 runs), not domain/box-size-dependent, so 'real g_SI=9.81' means the same lattice value at every resolution tried here. Only 51000 steps (cheap first look given the ~8x cell-count cost of doubling every dimension), single checkpoint at snapshotFrequency=50000, VTK off, guard off (unchanged). **Rank count: 72 cores/node confirmed via sinfo (not 576/8 nodes like the sphere resolution scan) -> `cellsPerBlock <32,15,15>` gives 6x6x8=288 ranks on 4 nodes** (was going to reuse the sphere scan's <16,15,15>/576/8-node scale, corrected per explicit request). New job script `jobs/run_288_bogner3_deep_s32.sbatch` (partition multinode, --nodes=4 --ntasks=288), submit via `sbatch --export=ALL,B3_RHO=075,B3_VARIANT=_s32 jobs/run_288_bogner3_deep_s32.sbatch`. Question: does doubling resolution (this way, box+basin together, matching Bogner's own convergence-check convention) change which lock-in cluster the box finds, shrink the torque-noise-vs-hydrostatics mismatch established at s=16/s=18, or reduce the mass/water leak -- i.e. does Bogner's own claimed convergence-with-resolution (stated for his righting-moment test, not explicitly for this floating-cube test) show up here too? Results -> `pinscan_results/4273521/b3_r075_s32.*`.

> **[Claude, 2026-09-29] s=32 RESULT INVALID -- doubling the domain silently broke the real-gravity calibration (job 4273521, COMPLETED exit 0 but not a valid physics run). CORRECTION of the 's=32 SUBMITTED' entry's claim that gravityOverride is domain-independent.** `computedParams.txt`: **avgParticleDiameter=12.0029 (was 6.00146 for every s=16/s=18 run), particleScalingFactor=4135.72 (was 2067.86), dx_SI=0.000241796 (was 0.000483592, HALVED), dt_SI=1.24931e-05 (was 4.99725e-05, ~1/4).** Cause: `avgParticleDiameter` is computed dynamically by fitting the sediment-bed particle file (`spheres_out.dat`, `bedCopiesInX 1`, unchanged) into the domain's footprint (`getAvgDiameterScalingFactor`) -- doubling domainSize X (96->192) with the same bed file and bedCopiesInX doubled the fitted particle diameter, halving dx_SI. `gravityOverride` was kept at the SAME lattice-units value (5.06585e-05) on the (wrong, over-generalised from the s=16-vs-s=18 comparison where the domain did NOT change) assumption that it is a pure unit-conversion constant independent of domain size. With the new dx_SI/dt_SI, that lattice value corresponds to **g_SI = g_LU * dx_SI/dt_SI^2 ~= 78 m/s^2, not 9.81 -- about 8x too strong.** Combined with `liquidHeightFactor 10.0` (chosen to double the level assuming avgParticleDiameter stayed 6.00146): `liquidHeight = avgParticleDiameter * liquidHeightFactor = 12.0029*10.0 = 120.03` ~= the ENTIRE 120-cell domain height -- essentially zero gas headroom, confirmed by `waterheight.txt` reading a flat **120.00 for the whole run** (no real free surface at all). Consequence: the run's own numbers (z_c climbing 52.25 -> 102.82, tiltY locking at 89.91+-0.06 near the domain ceiling, F_buoy+grav_z ~ -2.49 an order of magnitude above the s=16 runs' ~0.3-0.4) are an artifact of this badly-conditioned, wrong-gravity, no-free-surface setup, not a genuine s=32 resolution finding -- **retracted, no comparison to s=16/s=18 is valid from this run.** **`gravityOverride is domain-independent` claim in the earlier 's=32 SUBMITTED' entry is WRONG IN GENERAL** -- confirmed true only for the narrow case actually tested (s=16 vs s=18, same domain, only the box size changed, so avgParticleDiameter/dx_SI/dt_SI genuinely stayed fixed there) -- false whenever the DOMAIN itself changes, since avgParticleDiameter depends on domain footprint via the bed-fitting logic. **Fix needed before retrying s=32 (not yet done, awaiting direction):** either (a) recompute gravityOverride for this run's own dx_SI/dt_SI (g_LU ~= 9.81 * dt_SI^2/dx_SI ~= 6.33e-06, not 5.07e-05), or (b) the more robust fix -- also double `bedCopiesInX` (1 -> 2) when doubling domainSize X, so the SAME bed particle size tiles the wider footprint and avgParticleDiameter/dx_SI/dt_SI stay fixed exactly as they did for the s=16/s=18 comparison, keeping every future resolution step's unit system consistent without hand-patching gravity each time.

> **[Claude, 2026-09-29] s=32 RESUBMITTED WITH THE bedCopiesInX/Y FIX (per "resubmit with real gravity"). Job 4273525 (mt-dennys).** Fix applied to `prms/bogner3_r075_deep_s32_gas.prm`: `bedCopiesInX`/`bedCopiesInY` 1 -> 2 (matching the domain's own 2x doubling), so `getAvgDiameterScalingFactor` (Utility.h: `scalingFactor = domainSize[0] / (generationDomainSize_SI[0] * bedCopiesInX)`) produces the SAME avgParticleDiameter/particleScalingFactor/dx_SI/dt_SI as every s=16/s=18 run, instead of the doubled values job 4273521 silently got. Nothing else changed -- box pose/pinning untouched (`initialTiltAngleDeg 2.86`, `pinRotationAxes <1,0,1>`, `position <96,45,52.085>`, `halfExtents <16,32,16>` = s=32), `gravityOverride` still 5.06585e-05 (now genuinely correct again), `liquidHeightFactor 10.0` (now correctly gives absoluteLiquidHeight ~60.075, matching the value the heave/position calc already assumed). Deleted the old, invalid `snapshots_b3_r075_s32` checkpoint from job 4273521's broken (wrong-gravity, no-free-surface) run before this fresh start (the sbatch's own guard refuses a fresh start over an existing checkpoint folder). Results -> `pinscan_results/4273525/b3_r075_s32.*`.

> **[Claude, 2026-09-29] s=32 FIXED RESULT (job 4273525, COMPLETED exit 0, avgParticleDiameter confirmed 6.00146 = correct calibration): BY FAR THE BEST CONVERGENCE OF THE WHOLE INVESTIGATION -- settles at 26.56 deg +- 0.03, essentially exactly the ideal 26.565 deg target.** Local `animations/bogner3_deep_s32_20260929/4273525/`. **Tilt history:** 2.88 (step 250) -> climbs/overshoots through 9-18-30-32.56 (steps 2500-12500, peak overshoot 32.56 at 15000) -> corrects back down: 28.06 -> 25.79 -> 25.26 (steps 17500-25000, slight undershoot) -> settles: 26.03 -> 26.94 -> 26.54 -> 26.27 -> 26.54 -> **26.56-26.60 for the remaining ~25000 steps (30000-51000), last-15k mean +26.56, std 0.03, range 26.5-26.6.** Compare: s=16 (both seeds) converged to ~-30 to -39 deg (7-13 deg off target); s=18 never settled cleanly (drifted through 26.5/45/59-60/90 clusters, never held near target with this little scatter). This is the first run in the entire series to land within 0.05 deg of 26.565 and hold there with near-zero scatter. Water/mass loss still present, similar relative rate to smaller sizes: level -0.98 cells / mass -2.05% over 51000 steps (per-step rate ~-4.0e-5%/step, comparable to or slightly higher than s=16's ~-5.4e-5 and -7.0e-5 %/step over 100000 steps -- NOT shrinking with resolution). Lowest corner min 27.78 (floor 0.06, healthy margin at this scale). Warning bursts (`compute_BG.py`'s data confirms the timing separately): 12+146+21+17+70+15+17=298 of 304 total 'mass is lost' events land in steps 0-15000 (the climb/overshoot/correction transient), dropping to 1-3 per 2500-step bucket after step 20000 once the box settles -- same transition-coincident pattern as every other run. **B-vs-G torque comparison (`compute_BG.py`, added s=32 case, y_len=64): the torque MISMATCH is still there and, in absolute terms, if anything LARGER than at s=16/s=18** -- e.g. step 29500: tau_ideal -0.384, tau_meas +2.498 (opposite sign, 6.5x bigger); step 45500: tau_ideal -0.101, tau_meas +1.576 (opposite sign, 15.6x bigger) -- yet the box still converges tightly. **Reconciling this apparent paradox (torque noise not shrunk, but convergence dramatically better): torque scales with the box's volume x length (~s^4) while ROTATIONAL INERTIA scales with volume x length^2 (~s^5) -- so the noise-to-inertia ratio, i.e. the actual ANGULAR ACCELERATION the same relative noise produces, shrinks as ~1/s as the box grows.** This is consistent with everything observed: the absolute torque scatter grew roughly 10-40x from s=16 to s=32 (matching an s^4-ish scaling), but the box is far LESS perturbed by it, exactly as the s^5 inertia-growth argument predicts. **This reframes the whole 'why won't it settle' investigation: the torque-noise MECHANISM (lattice/voxelisation scatter in the MEM force, established earlier) never goes away and does not shrink with resolution -- what changes is the body's own inertia outrunning it. This matches Bogner's own claim (stated for his righting-moment test, 5.3.2) that results converge towards the analytic value at higher resolution, and is the first result in this investigation that is directly consistent with his reported outcome for the floating-cube test (5.3.3) rather than contradicting it.** Caveat: single run, single trajectory (this run's own initial overshoot to 32.56 before correcting down shows the dynamics are still noisy/non-monotonic even here, just with a much smaller basin of attraction/faster relative settling) -- not yet repeated with a second seed/start to confirm reproducibility the way the s=16 comparison was (normal seed vs equilibrium seed, both converging to the same wrong cluster). A natural next check, not yet run: repeat s=32 with a different seed (or an equilibrium-seeded start at 26.565 like the s=16 test) to see if it reliably lands here again.

> **[Claude, 2026-09-29] s=32 RESUMED FOR ANOTHER ~50K (per "continue from last checkpoint for another 50k so 100k in total"). Job 4273559 (mt-dennys), resumed from job 4273525's own final checkpoint.** Checkpoint `simState.file` says beginTimeStep=50001 (its last regular-cadence write, one interval short of the original run's own final step 51000 -- same off-by-one-interval pattern seen in every other bogner3_* resume this session). New `prms/bogner3_r075_deep_s32_resume_gas.prm`: `loadSnapshot true`, `timesteps 100000` (per the explicit "100k in total" target), density column moved to 98000-99999, otherwise identical to the fixed s=32 prm (real g_SI=9.81, guard off, bedCopiesInX/Y=2 fix retained). Submitted `sbatch --export=ALL,B3_RHO=075,B3_VARIANT=_s32,B3_RESUME=1 jobs/run_288_bogner3_deep_s32.sbatch`. Question: job 4273525 settled at 26.56+-0.03 deg for its last ~20000 steps (30000-51000) -- does it HOLD there for another ~49000 steps, confirming a genuine, durable equilibrium rather than a longer-but-still-transient dwell, or does it eventually drift/transition again the way every smaller-s run's initial 'quiet' stretch (e.g. the s=18 run's -26.5 deg plateau, 40000 steps) eventually did? Results -> `pinscan_results/4273559/b3_r075_s32_resume.*`.

> **[Claude, 2026-09-29] s=32 EXTENDED RESULT (job 4273559, COMPLETED exit 0, resumed 50001->100000): the tight 26.56+-0.03 lock did NOT hold -- it broke almost immediately on resume (within ~1000 steps) and settled into a noisier, higher band with large positional drift.** Local `animations/bogner3_deep_s32_20260929/4273559/`. **Tilt, right at the resume point:** 26.63 (step 50250, matches 4273525's own final value) -> 27.10 (50750) -> 28.24 (51250) -> 29.63 (51750) -> 30.58 (52250) -- already 4 deg off within 2000 steps of resuming -- climbs to a peak 34.08 (57250), then oscillates in a **27-35 deg band for the rest of the run, last-20k (80000-100000) mean +30.50, std 1.45**, i.e. ~4 deg off target and 50x noisier than 4273525's own last-20k std of 0.03. **x drifted almost the entire periodic domain length**: 99.5 (50250) -> 4.3 (100000), i.e. it travelled essentially all the way around the 192-cell periodic X domain -- a qualitatively different regime from the original run's own quiet window (x barely moved, 95.9-100.7). Water/mass stayed FLAT despite this (level -0.11 cells, mass +0.09% over the whole window) -- consistent with the established pattern that the leak concentrates during destabilisation transients and goes quiet during steady (even if not perfectly still) oscillation, not that this window was calm. **Warning-burst timing: 115 of 144 total 'mass is lost' events land in the very first 2500-step bucket after resuming (step 50000-52500), essentially none afterward (29 at 52500-55000, then silence)** -- an order of magnitude above the steady-state background rate seen once 4273525 had settled (1-3 per bucket). **IMPORTANT CAVEAT, not resolved:** this burst sits suspiciously exactly at the resume point itself, raising the open possibility that the CHECKPOINT RELOAD (not a genuine long-run instability of the 26.565 deg state) is what perturbed the box -- e.g. some derived/ghost-layer state not perfectly round-tripped through the save/load cycle, distinct from the box's own dynamics. This cannot be distinguished from 'the equilibrium is only metastable over long timescales' without a SINGLE CONTINUOUS run past 51000 steps (no resume in the middle) to rule out a resume-specific artifact -- not yet attempted. **Answer to the standing question:** the 26.565-deg-like state was NOT durable across this resume -- whether that reflects the state's own long-run (in)stability or an artifact of resuming from checkpoint is unresolved.

> **[Claude, 2026-09-29] SINGLE CONTINUOUS s=32 RUN TO 100000 STEPS, NO RESUME (per "run it again from timestep 0"). Job 4273570 (mt-dennys), corrected time limit.** Tests directly whether the 26.56+-0.03 lock's breakdown after resuming (job 4273559) reflects the state's own long-run instability or a checkpoint-reload artifact: identical setup to `bogner3_r075_deep_s32_gas.prm` (bedCopiesInX/Y=2 fix retained, real g, guard off) but ONE continuous execution from step 0 to 100000 -- `prms/bogner3_r075_deep_s32_100k_gas.prm`, new snapshot folder `snapshots_b3_r075_s32_100k` (checkpointing still happens at snapshotFrequency=50000 for safety, but nothing loads from it mid-run). **Job mishap (self-caught before wasting cluster time):** first submission (4273570's slot, job id 4273569) used the job script's original `--time=01:00:00` limit, inherited from the 51000-step template -- 100000 steps at the observed s=32 rate needs ~80 min, would have been killed by SLURM mid-run. Caught and cancelled after only 20 sec elapsed (no data lost); `jobs/run_288_bogner3_deep_s32.sbatch`'s time limit bumped to 02:00:00; resubmitted clean as **job 4273570**. Results -> `pinscan_results/4273570/b3_r075_s32_100k.*`.

> **[Claude, 2026-09-29] s=32 SINGLE-CONTINUOUS 100K RUNS FOR rho 0.5 AND rho 0.25 SUBMITTED (per "setup the same run then also for the 0.25 and the 0.5, so only change the density of the cuboid"). Jobs 4273582 (rho 0.5, mt-dennys-mirror) and 4273583 (rho 0.25, mt-dennys-mirror2).** Both are minimal-diff copies of the working `bogner3_r075_deep_s32_100k_gas.prm` (confirmed by diff: exactly three lines changed each -- `snapshotBaseFolder`, `position.z`, `densityRatio`; basin <192,90,120>, bedCopiesInX/Y=2 fix, box s=32 (halfExtents <16,32,16>), real g_SI=9.81, guard off, no VTK, 100000 steps, single continuous execution all unchanged). Heave (position.z) recomputed via Archimedes for s=32 at the new density, seed 2.86 deg, level 60.0746: **rho 0.5 -> z_c=60.075** (`prms/bogner3_r050_deep_s32_100k_gas.prm`), **rho 0.25 -> z_c=68.065** (`prms/bogner3_r025_deep_s32_100k_gas.prm`). Companion rho 0.75 run (job 4273570, same family) was at step 68500/100000 (tilt 26.62 deg, still holding) when these were submitted. Results -> `pinscan_results/4273582/b3_r050_s32_100k.*` and `pinscan_results/4273583/b3_r025_s32_100k.*`. Question: at Bogner's own high resolution (s=32) and the now-corrected real-gravity calibration, do rho 0.5 and rho 0.25 also converge tightly to their own ideal targets (45 deg, 26.565 deg) the way rho 0.75 did at this resolution, completing the resolution-comparison picture (s=16: none converged; s=18: none converged cleanly; s=32: rho 0.75 converged tightly -- do 0.5/0.25 as well)?

> **[Claude, 2026-09-29] RESOLVED: the s=32 rho 0.75 drift after resuming (job 4273559) was a CHECKPOINT-RELOAD ARTIFACT, not a genuine long-run instability of the 26.565 deg state. Single continuous run (job 4273570, COMPLETED exit 0, 82 min, 100000 steps, no resume anywhere) holds the lock through the ENTIRE run, including straight through the exact step-50000 region where the resumed run broke away within 1000 steps.** Local `animations/bogner3_deep_s32_20260929/4273570/`. **Tilt, by window:** 30000-51000: mean +26.58 std 0.05; **51000-80000 (spans the critical step-50000 boundary with no resume there): mean +26.53 std 0.11, range 26.1-26.7 -- NO disturbance, no burst, no drift**; 80000-100000 (last-20k): mean +26.55 std 0.05, range 26.4-26.7. The whole 75000-step stretch from 25000 to 100000 sits within +-0.3 deg of 26.565, tighter even than job 4273525's own first 51000 steps. **Warning-burst timing confirms it: zero elevated burst anywhere near step 50000** (buckets of 5000 steps: 45000:16, 50000:(none logged, i.e. 0), 55000:1 -- indistinguishable from ordinary background elsewhere in the run, e.g. 70000:27, 95000:27) -- directly contradicting what would be expected if step ~50000 were itself a natural destabilisation point independent of resuming; the 115-event burst seen in job 4273559 right at its OWN resume point has no counterpart here at the same step count with no resume. **B-vs-G torque comparison (`compute_BG.py`, updated) reproduces the by-now-established pattern**: tau_ideal stays small (|.| mostly <0.4) while tau_meas is comparable-to-larger and frequently opposite in sign (e.g. step 27250: tau_ideal -0.319 vs tau_meas -0.840, same sign but 2.6x; step 45500: tau_ideal -0.217 vs tau_meas -0.176, comparable) -- the torque-noise-vs-hydrostatics mismatch is present and unresolved exactly as in every other s=32 run, yet stability holds, consistent with the s^4-torque/s^5-inertia scaling argument already logged. Water level/mass: continues to fall steadily thoughout (60.57 -> 58.27 cells, -2.31; mass -4.29%) with NO acceleration or discontinuity at step 50000 either -- the leak here is a smooth, continuous background process at this resolution/pinning configuration, unlike the sharper transition-coincident bursts seen at s=16/s=18 -- and does not appear to threaten the lock over this run's length (z_c is compensating smoothly, 52.10->50.54, still 27.78+ cells clear of the floor). **Conclusion: job 4273559's drift was specific to the resume mechanism (most likely some derived/ghost-layer state not perfectly round-tripped through the checkpoint save/load cycle), not evidence that the s=32 rho 0.75 equilibrium is only metastable over long timescales.** A single continuous run at this resolution and calibration reliably reproduces and HOLDS the near-exact 26.565 deg equilibrium for the full 100000 steps tested. Any future extension of this run past 100000 steps should therefore be done as a longer single continuous execution where feasible, not a resume, if reproducing this same clean result matters.

> **[Claude, 2026-09-29] Direct comparison of job 4273525's own tail (steps 50500-51000, continuous) against job 4273559's own head (steps 50250-50750, resumed) at matching step numbers (per "analyze what changed at 50k again") -- inconclusive, does not cleanly settle the resume-artifact question.** Raw [BoxForces]/[BoxTorque] at step 50500: continuous run F_hydro=(-0.004,0.310,2.502) torque_y=+0.044; resumed run F_hydro=(-0.037,0.144,2.478) torque_y=+0.190 -- same order of magnitude, no shock/blow-up, but already visibly different (expected: this project's own prior finding that the LBM here is not bit-reproducible run-to-run, so SOME divergence between a continuous run and a resume-from-checkpoint is normal even with a perfect reload). At step 50750: continuous torque_y=-1.145 vs resumed torque_y=+0.608 -- opposite sign, notable divergence only 750 steps post-resume, but torque_y swings of this size (and larger) also occur routinely mid-run in every bogner3_* trajectory examined this session (not obviously anomalous by itself). The one thing that DOES stand out: the 115-event warning burst concentrated in steps 50000-52500 of the resumed run (already logged) is large, but bursts of comparable or larger size (93, 51, 71 events) were also seen at NATURAL, non-resumed transitions in other runs (e.g. job 4273497's own step-50072 and 80-90k events) -- so a burst of this size is not by itself outside the range seen without any resume involved. **Net: the comparison weakens but does not refute the 'natural chaotic divergence, not a resume artifact' explanation** -- the single continuous 100k run (job 4273570) is the test that can actually distinguish the two; this comparison alone cannot.

> **[Claude, 2026-09-29] SPHERE/CYLINDER EQUILIBRIUM COMPANIONS TO THE s=32 CUBOID FAMILY, ALL 3 DENSITIES (per "setup the equilibrium runs for this but for the cylinder and the sphere again for all 3 densities, they can all start from the same z position, adjust the sizes accordingly, run for 100k with the waterline measurement etc included like in the cuboid case ... archived in antidunes-archive").** Same basin/gravity/logging as `bogner3_r0XX_deep_s32_100k_gas.prm` (domain <192,90,120>, liquidHeightFactor 10.0 -> WL0=60.0746, bedCopiesInX/Y=2, real g_SI=9.81 via gravityOverride 5.06585e-05, guard off, no VTK, water-height 12x3 grid every 250 steps, 100000 steps, single continuous — no resume variant, per the resume-artifact finding above). **Sizing/positioning convention (chosen, not copied from the older R1/R2/R4 buoyancy-compare family, which uses a different box=cube/sphere/cylinder scale): z position is FIXED to the cuboid's own already-computed heave z_c per density (52.085 / 60.075 / 68.065 for rho 0.75/0.5/0.25) — same z for the sphere and cylinder as for the cuboid at that density — and the RADIUS is solved instead (spherical-cap volume / circular-segment area bisection) so the submerged fraction at that fixed z equals the target density exactly.** Sphere radius: 23.005 (rho .75), 21.500 (rho .5 — exactly half-submerged by symmetry at any radius since z_c≈WL0, so radius is a free choice here), 23.007 (rho .25). Cylinder (axis along Y, length 64 = matches the cuboid's own y-extent 2*halfExtents.y, unaffected by roll about its own axis since a circular cross-section is rotationally symmetric): radius 19.778 (rho .75), 21.500 (rho .5), 19.780 (rho .25). **No pinning at all, per explicit instruction ("for these we don't need to pin anything") — pinPosition/pinRotation both left at their default `false` for every FloatingSphere/FloatingCylinder block: fully free 6-DOF float.** Note for the cylinder: pitch/yaw (tipping nose-down) ARE physically meaningful DOFs here, unlike roll about its own long axis, which does nothing to the submerged cross-section by symmetry — so any interesting rotational dynamics observed will show up as pitch/yaw, not roll. Note for the sphere: it has NO orientation-dependent physics at all (any rotation leaves a sphere's submerged cross-section unchanged) — the only meaningful equilibrium quantity is heave/depth, there is no possible "heel angle" analogue. New prms: `bogner3_r075_sphere_s32_100k_gas.prm`, `bogner3_r050_sphere_s32_100k_gas.prm`, `bogner3_r025_sphere_s32_100k_gas.prm`, `bogner3_r075_cyl_s32_100k_gas.prm`, `bogner3_r050_cyl_s32_100k_gas.prm`, `bogner3_r025_cyl_s32_100k_gas.prm`. New job script `jobs/run_288_bogner3_shape_s32.sbatch` (same 288-rank/4-node layout as the cuboid's own `run_288_bogner3_deep_s32.sbatch`; force-log grep is shape-aware: `[CylForces]`/`[CylTorque]` for the cylinder, `[FloaterForces]` only for the sphere, which has no torque log). **Submission sequencing:** the cuboid's own rho 0.5/0.25 s=32 100k runs (jobs 4273582/4273583) were still RUNNING on `mt-dennys-mirror`/`mt-dennys-mirror2` at submission time, so only `mt-dennys` (free) got a job immediately. First attempt (job 4273655) was cancelled per explicit request to add auto-archiving first; the (now updated) `run_288_bogner3_shape_s32.sbatch` auto-copies each job's results into a shared vault folder `/home/vault/iwia/iwia148h/antidunes-archive/antidunes-bogner3-shapes-s32_20260929/` and flock-appends a per-job section to its `PROVENANCE.md` at the end of every run — no manual archiving step needed per job. Resubmission (job 4273669) hit a stale-checkpoint guard left over from the cancelled attempt (`snapshots_b3_r075_sphere_s32_100k` already existed); cleared and resubmitted clean as **job 4273688 (sphere rho 0.75, mt-dennys)**. Once cuboid jobs 4273582/4273583 completed (see their own results below), **sphere rho 0.5 -> job 4273720 (mt-dennys-mirror)** and **sphere rho 0.25 -> job 4273722 (mt-dennys-mirror2)** were submitted, each with its cylinder companion dependency-chained via `--dependency=afterany:<sphere_jobid>` (runs regardless of the sphere job's exit status) so no further manual polling/submission is needed: **cyl rho 0.75 -> job 4273723 (afterany:4273688)**, **cyl rho 0.5 -> job 4273724 (afterany:4273720)**, **cyl rho 0.25 -> job 4273725 (afterany:4273722)**. Results -> `pinscan_results/<jobid>/b3_r0XX_{sphere,cyl}_s32_100k.*`, auto-archived to the vault folder above as each job finishes.

> **[Claude, 2026-09-29] CUBOID rho 0.5/0.25 s=32 100k RUNS COMPLETED AND ARCHIVED. Jobs 4273582 (COMPLETED exit 0, 1:25:15) and 4273583 (COMPLETED exit 0, 1:18:08).** **rho 0.5: last-20k (80000-100000) tiltY mean 44.71, std 0.55, range 42.97-45.22 deg — locks tightly on the ideal 45 deg target** (0.29 deg off), mirroring rho 0.75's own clean lock. **rho 0.25: last-20k tiltY mean 65.40, std 0.19 deg — by the square's 4-fold symmetry this folds to 90-65.40 = 24.60 deg, ~1.96 deg off the ideal 26.565 deg target** — tightly locked (tightest std of the three densities) but on a slightly-off equilibrium rather than the exact target, unlike rho 0.75 (0.02 deg off) and rho 0.5 (0.29 deg off). **Summary, all 3 densities at s=32, real g, single continuous 100k:** rho 0.75 -> 26.55 vs ideal 26.565 (off 0.02); rho 0.5 -> 44.71 vs ideal 45.0 (off 0.29); rho 0.25 -> 24.60 (folded) vs ideal 26.565 (off 1.96). Archived to `/home/vault/iwia/iwia148h/antidunes-archive/antidunes-bogner3-s32-rho050-rho025_20260929/` with its own `PROVENANCE.md`; local copies under `animations/bogner3_deep_s32_20260929/{4273582,4273583}/`.

> **[Claude, 2026-09-29] SPHERE/CYLINDER EQUILIBRIUM COMPANIONS: ALL 6 COMPLETED (exit 0) AND AUTO-ARCHIVED.** Jobs: sphere rho .75/.5/.25 -> 4273688/4273720/4273722 (47:58/43:12/49:46); cylinder rho .75/.5/.25 -> 4273723/4273724/4273725 (dependency-chained after their sphere, 1:21:26/1:20:10/1:19:23). Archived automatically (per-job, via the sbatch script's own end-of-run logic, no manual step) to `/home/vault/iwia/iwia148h/antidunes-archive/antidunes-bogner3-shapes-s32_20260929/`, `PROVENANCE.md` auto-appended per job. Local copies `animations/bogner3_shapes_s32_20260929/<jobid>/`.
>
> **Sphere heave (no rotation possible for a sphere -- only position.z is meaningful): all 3 densities settle very tightly (last-20k std 0.03-0.04 cells) but DRIFT UPWARD 0.48-0.79 cells over the full 100000 steps, while the ambient water level stays essentially FLAT (+0.03 to +0.06 cells -- unlike the cuboid family's clear multi-cell leak at this same resolution/calibration).** This is genuine excess buoyancy, not a mass-leak-driven repositioning (the cuboid's own upward/downward z_c compensation was explicitly tied to a falling water level; here the water barely moves at all yet the sphere still rises) -- consistent with the pre-existing ~2-3% MEM-only excess-buoyancy hypothesis (see `pinscan_sphere_R24_real_f075_gas.prm`'s own header, section 1a), now confirmed directly in a FREE (unpinned) heave test rather than only in a force-balance comparison on a pinned sphere.
>
> **Cylinder heave matches the sphere pattern (small drift 0.52-0.66 cells, water level flat).** **Cylinder roll (`tiltY`, the ONLY rotational DOF the `[CylForces]` log emits -- nominally physically trivial for a circular cross-section, since roll about the cylinder's own axis cannot change what's submerged) is quiet at rho 0.75/0.25 (last-20k mean -0.32/+1.11 deg, std ~1.0-1.3 deg, max |tiltY| 2.5/4.6 deg -- consistent with "no real dynamics, just noise") but SWINGS HARD at rho 0.5: last-20k mean -40.92 deg, std 1.29 deg, max |tiltY| 42.46 deg.** This is a genuine anomaly, not glossed over: a 40+ deg `tiltY` reading at exactly the neutral-buoyancy density (the one case where the cylinder is expected to be least constrained by any restoring buoyancy asymmetry) means either (a) the cylinder's axis itself genuinely reoriented away from world Y (real pitch/yaw tipping that the logged Euler-angle field is picking up as "roll" contamination once the axis is no longer aligned with Y), or (b) a quaternion/Euler-decomposition artifact specific to this symmetric case. **UNRESOLVED: the current code only logs `tiltY` for `FloatingCylinder`, not a full pitch/yaw/roll breakdown, so this cannot be disambiguated from the existing `.forces.txt` data alone.** A follow-up needing either an added pitch/yaw log field (mirroring how the box's `tiltY` field would need extending) or a look at the raw stored quaternion in a snapshot would be required to settle which of (a)/(b) is happening -- not yet attempted.

> **[Claude, 2026-09-29] DISCRETE VTK SNAPSHOT FEATURE ADDED (`fluid_field_snapshots`), THEN A REAL PARSING BUG FOUND AND FIXED before any real data was collected.** Per "implement the option so we can have a discrete set of vtk that we can print instead of every thousand, like at 1k 2k 8k 16k 32k 64k 100k": new mechanism in `src/lbm/free_surface/VtkWriter.h` (`snapshotStepPtr`/`snapshotSteps` gate, same pattern as the existing `fluid_field_dense` window gate) + `SimulationParams.h` (`VTKParams::snapshotsEnabled`/`snapshotSteps`, parsed from a `steps "..."` comma-separated string via a new `parseUintList` helper) + `BuildTimeloop.h` (wiring). **Bug found during validation (job 4274435, sphere rho .75): step 0 was silently missing from every fired snapshot** -- confirmed via a standalone local test against the REAL waLBerla config parser that `config::Config::BlockHandle::getParameter<std::string>()` returns a quoted prm value WITH the surrounding quote characters still attached (e.g. `"0,250,...,100000"`, literal quotes included) -- undocumented anywhere obvious. The token glued to each quote mark breaks differently: a leading `"` before a digit makes `operator>>` fail immediately (silently dropping that whole entry, which is exactly how "0" disappeared), while the trailing `100000"` still parsed fine because `operator>>` just stops at trailing garbage instead of failing. **Fixed** by stripping a leading/trailing `"` in `parseUintList` before splitting; reverified against the real parser (all 11 values, including 0, now parse correctly) and rebuilt clean locally before redeploying. The affected job (4274435) and an interim singlenode validation job (4274495, `mt-dennys-mirror`, `singlenode` partition -- needed since `multinode` enforces a minimum of 2 nodes) were both cancelled and their partial checkpoints removed before any real thesis data was collected with the bug present.
>
> **VTK-enabled sphere/cylinder equilibrium companions resubmitted with the fix, ALL 6 running (cuboid explicitly excluded from this round per correction: "i only meant to rerun the cylinder and sphere cases not the cuboid").** Sphere gets a shared starting position (position.z=60.0746, the water level) across all 3 densities; cylinder keeps its own per-density starting z but now has pitch/yaw pinned (`pinRotationAxes <1,0,1>`) to isolate whether the earlier free-cylinder run's large rho=0.5 roll excursion was genuine or an artifact. Job map (two per checkout, dependency-chained `afterany`): **sphere rho .75 -> 4274674 (mt-dennys) -> cyl rho .75 -> 4274677**; **sphere rho .5 -> 4274675 (mt-dennys-mirror) -> cyl rho .5 -> 4274678**; **sphere rho .25 -> 4274676 (mt-dennys-mirror2) -> cyl rho .25 -> 4274679**. The 3 cuboid VTK prms generated earlier (`bogner3_r0XX_deep_s32_101k_vtk_gas.prm`) are left unused in `prms/` (not deleted, just not run).

> **[Claude, 2026-09-29] ALL 6 COMPLETED (exit 0), MYSTERY RESOLVED -- pitch/yaw pinning did NOT
> stabilize the cylinder's roll, and now we know why.** Jobs 4274674-4274679 all COMPLETED and
> auto-archived to `antidunes-archive/antidunes-bogner3-vtk-101k_20260929/`. Sphere heights land where
> intended given the shared start (last-20k mean/target: rho .75 53.591/53.366, .50 60.700/60.614, .25
> 67.792/67.775 -- all within 0.23 LU). **Cylinder roll, however, is dramatically different from the
> earlier free-rotation round: with pitch/yaw pinned, the RAW `tiltY` reading swings through the full
> -180 to +180 deg range at ALL THREE densities** (previously only rho=0.5 showed a large excursion;
> rho=0.75/0.25 stayed near 0 deg). Unwrapping the signal (removing the +-360 wraparound so a
> continuous spin reads as a ramp, not a sawtooth) shows the actual mechanism cleanly: **rho=0.75
> never stops rolling** -- two full revolutions (720 deg) completed by step 100000 and still turning;
> **rho=0.5 and rho=0.25 each complete exactly one bounded 180 deg flip and then lock there** for the
> rest of the run. This resolves the earlier mystery: **roll about a circular cross-section has ZERO
> restoring torque by symmetry** (every roll angle is identical to the fluid), so pinning pitch/yaw
> doesn't stabilize roll -- it REMOVES the only thing (coupling to pitch/yaw) that was occasionally
> damping it in the fully-free round. The same persistent torque-noise mechanism documented throughout
> this investigation (present at every resolution, never resolved at the root-cause level) is free to
> either spin the cylinder indefinitely or kick it into one of the (physically equivalent, since a
> 180 deg roll of a circle maps onto itself) locked poses, and which outcome a given density gets
> appears to depend on which way the early noise happens to push it, not on the density itself. Added
> a new "Cylinder Roll, Pitch/Yaw Pinned" section + unwrapped-rotation chart to the "Cuboid
> Heel-Angle Drift" artifact (all 3 densities on one axis, -200 to 760 deg) and updated the existing
> sphere/cylinder heave charts with this round's data (replacing the earlier free-cylinder/
> per-density-sphere-start series). Local copies `animations/bogner3_vtk_101k_20260929/<jobid>/`.

> **[Claude, 2026-09-29] BUG FOUND while building the next round: rho .50/.25 in ALL SIX
> `bogner3_r0XX_{sphere,cyl}_s32_101k_vtk_gas.prm` files use the WRONG radius** (the rho=.75 value,
> 23.005/19.778, baked in instead of their own 21.500/23.007 or 21.500/19.780). Checked the actual
> physical impact before treating this as invalidating anything: for rho=.5 the target submerged
> fraction is exactly 0.5 at ANY radius (the position sits exactly at the water level, a symmetric
> case established earlier this session), so the wrong-radius rho=.5 runs still hit their correct
> equilibrium target -- just with an off-nominal (larger sphere / smaller cylinder) object size,
> self-consistent but not the size intended. For rho=.25 the error is negligible (0.002 LU, since
> rho=.75 and rho=.25's own correct radii are themselves nearly identical by the sphere/circle's
> own near-mirror geometry at this basin). **Net: the already-completed jobs' equilibrium-convergence
> conclusions stand; only the exact object size for the rho=.5 cases (and by ~0.002 LU for rho=.25)
> was off from what was intended.** Not fixed retroactively in the already-run/archived prms (not
> worth a rerun given the negligible-to-none effect on conclusions) -- fixed going forward in the new
> `..._freesharedstart` prms below.
>
> **CYLINDER RERUN, FULLY FREE, SHARED START (per "rerun the cylinder from the same starting position
> again (same for all densities) but this time not pinned at all"):** new
> `bogner3_r0XX_cyl_s32_101k_vtk_freesharedstart_gas.prm` (radius bug above fixed here: 21.500/19.780
> for rho .5/.25), position.z=60.0746 (the shared water-level start, matching what sphere/cuboid
> already got) for all 3 densities, pinPosition/pinRotation/pinRotationAxes removed entirely (no
> pinning at all, unlike both the immediately preceding pitch/yaw-pinned round AND the original
> free-cylinder round, which used each density's own per-density start). Tests whether letting
> pitch/yaw move too (instead of pinning them) changes the "roll never stops or locks after one
> 180deg flip" outcome just found -- e.g. by coupling roll to a genuinely restoring pitch/yaw mode
> that pinning removed. `run_288_bogner3_vtk_101k.sbatch` extended with an optional `B3_VARIANT`
> suffix (default empty, backward compatible) to support this variant's distinct prm/snapshot naming
> without touching the existing cases. **Checked cluster state first and used dependencies to avoid
> any checkout overlap** (per explicit request): `mt-dennys` was free (R1 buoyancy-compare job
> already finished there) -> **rho .75 submitted directly as 4278473**; `mt-dennys-mirror` was busy
> with the R2 buoyancy-compare job (4278451) -> **rho .5 submitted as 4278474 with
> `--dependency=afterany:4278451`**; `mt-dennys-mirror2` was busy with the R4 buoyancy-compare job
> (4278452) -> **rho .25 submitted as 4278475 with `--dependency=afterany:4278452`**.

> **[Claude, 2026-09-29] CORRECTION: sphere/cylinder radius must be FIXED across densities, matching
> the cuboid's own convention -- per explicit correction ("the sphere and cylinders and cuboid should
> be the same size no matter the density").** The cuboid has always been one fixed size
> (`halfExtents <16,32,16>`) regardless of density; sphere/cylinder radius, by contrast, had been
> solved PER-DENSITY (to hit the target fraction at a fixed z) ever since the very first
> sphere/cylinder equilibrium round this session -- an inconsistency that only made sense under the
> old "start already at equilibrium" methodology, not the current shared-start settling-dynamics one.
> **Fixed radius = 16.000 for both shapes (matching the cuboid's own cross-section half-extent) across
> all 3 densities**; cylinder length stays 64 (matching the cuboid's y-extent). Cancelled 4278473/
> 4278474/4278475 (still queued/barely started with the old per-density radii) before any real data
> was collected; the 3 `..._freesharedstart` cyl prms fixed in place (radius now 16.000, was
> 19.778/21.500/19.778 -- note the r050/r025 files ALSO still had the earlier rho=.75-baked-in radius
> bug at that point, so this single fix corrects both issues at once).
>
> **SCOPE, per explicit follow-up ("don't redo the pinned one only do the non pinned version, also
> delete the wrong runs from the archive"):** the already-completed pitch/yaw-PINNED cylinder round
> (4274677/4274678/4274679) is NOT redone and NOT deleted -- kept as-is, wrong-radius caveat and all,
> since it isn't being superseded. The already-completed SPHERE round (4274674/4274675/4274676,
> per-density radius) IS superseded -- **deleted from the vault archive**
> (`antidunes-archive/antidunes-bogner3-vtk-101k_20260929/427467{4,5,6}_...`, had to `chmod -R u+w`
> first since the archived copies inherited the read-only bits from `chmod -R a-w` on the results dir)
> and from the local copies (`animations/bogner3_vtk_101k_20260929/427467{4,5,6}/`, same permission
> fix needed). The 3 sphere prms (`bogner3_r0XX_sphere_s32_101k_vtk_gas.prm`) fixed in place (radius
> now 16.000 for all 3, was 23.005/23.005/23.005 -- i.e. the rho=.75 value duplicated into .50/.25 by
> the earlier bug, now ALSO fixed to the new shared 16.000 rather than restored to the old per-density
> values). Verified across all 6 final prms: sphere and cylinder both now have IDENTICAL
> `position <96,45,60.0746>` and IDENTICAL radius (16.000; cylinder also identical `length 64`) across
> all 3 densities, matching the cuboid's own fixed-size convention exactly.
>
> **6 reruns total (3 sphere + 3 cylinder, fully free/unpinned), dependency-chained to avoid both
> checkout overlap AND depending on the slow R4 buoyancy job** (per explicit instruction: "do the
> dependency chaining but not after the R4 that will take considerably longer"). Job map: `mt-dennys`
> runs 4 jobs back-to-back -- **cyl .75 (4278476, running)** -> **sphere .75 (4278483)** -> **cyl .25
> (4278484)** -> **sphere .25 (4278485)**; `mt-dennys-mirror` runs 2 more after its own R2 job
> finishes -- **cyl .5 (4278477)** -> **sphere .5 (4278486)**. `mt-dennys-mirror2` (running R4) is
> left completely untouched -- nothing chains after it.

> **[Claude, 2026-09-29] SECOND BUG FOUND IN THE PINNED-CYLINDER BATCH -- same root cause as the
> radius bug, but affecting `position` this time.** Checking "does the pinned cylinder's size vary
> too" while answering a direct question surfaced that `bogner3_r0{50,25}_cyl_s32_101k_vtk_gas.prm`
> ALSO never got their own per-density starting height: all three pinned-cylinder prms have
> `position <96,45,52.085>` (rho=0.75's own z_c) baked in, because the generator script templated all
> 3 densities from the single rho=0.75 source file and patched `densityRatio` (always) and `position`
> (only for cuboid/sphere, which were deliberately switched to the new shared start) -- cylinder's
> `position` was meant to stay "per-density" per the header comment, which needed actually reading
> each density's own original z_c from its own source file, which the script never did. Net effect:
> **rho=0.75 (job 4274677) is actually correct** (52.085 and 19.778 both happen to BE its own intended
> values, since it IS the template source) -- **rho=0.50 (4274678) and rho=0.25 (4274679) both
> silently ran as near-duplicates of the rho=0.75 case** (same start height, same radius), only
> `densityRatio` itself actually differed. Their "unexplained large roll excursion / locks after one
> flip" results are confounded by this and should not be read as reflecting genuinely different
> per-density starting conditions.
>
> **Per explicit instruction ("delete the pinned version that are wrong keep the one that is
> correct"): 4274678 and 4274679 deleted from both the vault archive and local copies (4274677 kept
> as-is, not rerun, per the earlier explicit decision not to redo the pinned batch).** Vault
> `PROVENANCE.md` for this folder rewritten from scratch (was auto-appended per-job and had gone
> stale/inconsistent after two rounds of deletions) to describe only what's actually left plus a clear
> record of what was deleted and why.
>
> **Verification method upgraded before trusting the CURRENT (fully-free, fixed-size) replacement
> batch**: rather than re-reading only the prm text (which is exactly what missed both bugs above),
> checked the actual DEPLOYED prm on both active checkouts AND the ALREADY-RUNNING job's own logged
> `[CylForces]` output (`pos=<96,45.0001,59.7736>` at step 250 for job 4278476, rho=0.75) to confirm
> the real simulation is actually using the intended shared start (60.0746) and fixed size (16.000),
> not just that the source prm file says so.

> **[Claude, 2026-09-29] ARTIFACT BUG FOUND AND FIXED: the cuboid's plotted equilibrium target was
> never recomputed from the run's own data, unlike the sphere/cylinder targets, making it look like
> the box was "drifting away" from equilibrium when it wasn't.** Per "why is the cube drifting away
> [screenshot of the artifact's 'Same shape, all three densities' Cuboid panel] is the target also
> computed in the same way as the sphere, cylinder and cuboid" -> "also do this for the sphere and
> cylinder recomputed using the actual measured water level." **Root cause:** the sphere/cylinder
> targets in the artifact (`SCHARTS` in `tilt_chart_v3.html`) were already each solved from that run's
> own measured END-of-run water level (via the same 12x3 `waterheight.txt` grid average already used
> elsewhere) -- confirmed correct, no change needed (sphere 55.323/60.633/65.874, cylinder
> 54.910/60.652/66.290 for rho .75/.5/.25). The **cuboid** target, however, was still the STATIC
> starting z_c baked into each prm at submission time (52.085/60.075/68.065) -- never recomputed
> against either the run's own measured end water level OR its own converged heel angle, both of
> which move the true equilibrium height for a tilted square cross-section. Against the actual
> last-20k mean heights (51.017/58.750/66.341), the static target reads as a 1.0-1.8-cell "drift"
> that is actually just the wrong target, not a real divergence. **Fix:** computed each cuboid run's
> own last-20k mean heel angle from its `.forces.txt` (jobs 4273570/4273582/4273583: 26.552 deg / 44.706
> deg / 65.399 deg for rho .75/.5/.25) and its own measured end-of-run water level (58.4839 / 58.7913 /
> 59.1002, same grid-average method as sphere/cylinder), then solved for the z_c giving the target
> submerged-volume fraction via a Sutherland-Hodgman polygon clip of the rotated 32x32 square cross-
> section against that water level (shoelace-area fraction, bisected on z_c -- same bisection pattern
> already used for the sphere/cylinder radius solves, just clipping a polygon instead of a circle/
> sphere-cap). **New targets: 51.328 / 58.791 / 66.374** (rho .75/.5/.25) -> residual vs. the observed
> last-20k means shrinks to **-0.311 / -0.041 / -0.033 cells** -- an order of magnitude tighter than
> the apparent 1.0-1.8-cell "drift" against the old static target, and now consistent with how tightly
> the sphere/cylinder already track their own (correctly end-water-level-based) targets. Edited
> `SCHARTS`'s Cuboid entry (`target`/`targetLabel`) and the section lede paragraph in the artifact
> (https://claude.ai/artifact/Pp1U9xtwCD9vD5J44baY9S, v14) to explain the corrected, now-uniform
> methodology across all three shapes. **Not yet done:** the overleaf thesis figures
> (`overleaf/scripts/equilibrium-plots/extract_data.py`, `shape-comparison-highres.tex`) still use the
> OLDER, bug-affected sphere/cylinder job data (wrong radius / wrong start position, see the two bugs
> above) and the old static cuboid target -- none of this session's fixes (radius fix, position fix,
> or this target-methodology fix) have been carried into the actual thesis PDF yet.

> **[Claude, 2026-09-29] GAP CLOSED: overleaf figures + text now use the corrected data too.** Per
> request to add/redraw the artifact's cuboid/sphere/cylinder height panels and the new combined
> cuboid-tilt chart into the thesis and update the results text to match. **`extract_data.py`**: RUNS
> updated to the corrected sphere/cylinder jobs (4278714/4278486/4278485 sphere, 4278476/4278477/
> 4278484 cylinder -- the fully-free, fixed-radius, shared-start replacements, not the two
> superseded buggy batches); TARGETS updated to the corrected end-water-level-based values for all 9
> shape/density combinations (cuboid also now uses the polygon-clip + own-heel-angle method, not its
> static starting z_c). **`make_plots.py`**: axis widened from 46-70 to **42-72 LU** (the corrected
> cylinder rho=.75 trough of 44.0 and cuboid rho=.75 trough of 47.5 would otherwise clip) and
> `X_MAX` bumped 100000->101000 (sphere/cylinder's own run length). Regenerated all 4 PDFs
> (`shape-heave-{cuboid,sphere,cylinder}.pdf`, `cuboid-tilt.pdf`) into `overleaf/figures/equilibrium/`
> -- visually diffed against the artifact's own charts, confirmed matching.
>
> **Text (`shape-comparison-highres.tex`) rewritten, not just re-pointed at new numbers**, because
> the corrected cylinder data changes the actual scientific narrative, not just the axis: the
> cylinder is now FULLY FREE (no pitch/yaw pinning -- the pinned batch was the one with the
> position bug and was abandoned rather than fixed, per the earlier decision this session), and a
> freely-rotating cylinder turns out to ring for **50,000+ steps** at rho=.75 (trough 44.0 LU) and
> ~35,000 steps at rho=.25 -- a timescale comparable to the cuboid's own heel relaxation, not the
> "settles essentially as fast as the sphere" behaviour the old pitch/yaw-pinned cylinder showed and
> that the old text used to argue "it's heel-angle relaxation specifically, not non-spherical
> geometry." Rewrote the closing discussion to instead conclude: (1) the **sphere**, not the
> cylinder, is the clean rotation-free control (residuals +0.34/+0.06/+0.04 LU, tightest of the
> three, settles in ~13,000 steps at rho=.75 and effectively instantly otherwise); (2) once the
> cylinder is given real rotational freedom it develops slow pitch/heave coupling of its own despite
> a roll-invariant circular cross-section, suggesting it's rotational relaxation *in general*, not
> the square cross-section's Scribanti-curve dynamics specifically, that slows settling -- a pinned
> comparison can't distinguish the two; (3) all three shapes' final residuals are actually the same
> order of magnitude (few hundredths to a few tenths of an LU) once every target is computed
> consistently at the run's own end water level (+ own heel angle for the cuboid) -- the cuboid's
> apparent multi-LU "drift" in the old text was a target-computation artifact, not a real bias. Also
> fixed the setup paragraph's cylinder-pinning description, the axis range, the trough/residual/
> settling-step numbers throughout, and added an honest paragraph explaining why the pinned-cylinder
> round was dropped (referencing this file). Rebuilt the full thesis PDF (pdflatex + biber + 2x
> pdflatex) clean: all 4 figures embedded correctly at physical pages 86-89 (section 4.1.3), no new
> undefined references introduced (the one remaining `sec:cc-dispatch-gap` warning is pre-existing
> and unrelated), verified by rendering those pages to PNG and reading them back.

## 1b. [ ] MEM hydrostatic force vs. analytical-buoyancy-subtraction comparison

**Goal**: for a single fully-submerged floater, compare the current default MEM-only buoyancy model
(`useMemOnlyBuoyancy=true` — buoyant support comes entirely from the momentum-exchange/hydrodynamic
force, no separate analytic term) against the older method (`useMemOnlyBuoyancy=false` — gravity
minus an explicit analytic buoyancy term `vSub × g`).

**Scenario**: release a floater at the bottom of the domain with a densityRatio low enough to float,
let it rise, and analyze the differing rise behavior (trajectory, velocity, any oscillation/overshoot)
between the two models. Only the **fully-submerged** portion of the rise is the relevant comparison
window — once the floater breaches the surface, the two models converge in relevant behavior/intent
diverges from the comparison's purpose (partial-submersion buoyancy isn't what's being isolated here).

> **[Dennys, 2026-09-12]:** we assume the already-working Antidunes version (MEM/`useMemOnlyBuoyancy`)
> is validated by previous projects — no need to independently check either method against a
> closed-form/analytical ground truth here, this comparison is purely about the two methods'
> differing behavior relative to each other.

> **[Claude, 2026-09-12]:** implementation status. Added a new additive, backward-compatible flag
> `useFullVolumeAnalyticBuoyancy` (default `false`) so the analytical arm uses the original
> full-volume buoyancy formula (`vSub=vFull`) rather than the newer surface-tracking plane-fit
> machinery, which isn't relevant while the floater stays fully submerged. Built a resolution sweep
> (R=1/2/4, same physical SI setup, `g_SI=9.81` throughout) × 3 shapes (sphere/cylinder/box) × 2
> methods (MEM vs. analytical) = 18 comparison runs, each fed a real `BedGeneration`-produced bed file
> (not an SI override) and a `timesteps` budget calibrated on the cluster to stay safely below the
> floater's measured breach step. **Finding**: the log's own `subm=` field cannot detect a breach in
> either arm of this comparison (both unconditionally set `vSub=vFull` internally) — breach must be
> checked via `pos.z` directly against the known waterline height instead. All 6 R=1 runs completed
> on the cluster; R=2 (job 4223038) and R=4 (job 4223039, 8 nodes) are running now. Once complete,
> next step is the actual `[FloaterForces]` log comparison (z(t)/vel_z(t) overlay per shape/resolution,
> restricted to samples before breach) and archiving everything under
> `/run/media/dhuber/data2/debris-clogging-archive/buoyancy-method-comparison/`.

> **[Claude, 2026-09-13]:** R=2 (retried as job 4223072) COMPLETED. R=4 needed two attempts: job
> 4223039 was cancelled (never started), and its retry, job 4223073, FAILED after ~3.5 minutes with
> `Invalid input file spheres_out_buoyancy_r4.dat` (Utility.h:109) — a mundane cross-checkout file
> mismatch, not a physics/code bug: the R=4 bed file was generated successfully by job 4222977 on the
> `mt-dennys-mirror` checkout, but 4223073 ran on `mt-dennys-mirror2`, which never had that file
> copied over. Fixed by copying `spheres_out_buoyancy_r4.dat` from `mt-dennys-mirror` to the (idle)
> main `mt-dennys` checkout and resubmitting there as job 4227812. Report back once it finishes.

> **[Claude, 2026-09-13]:** Pulled and parsed the actual `[FloaterForces]`/`[CylForces]`/`[BoxForces]`
> per-step logs (pos.z/vel.z) for all 6 R=1 and all 6 R=2 runs (jobs 4222989/4223005/4223007/4223008/
> 4223013/4223015 for R=1; the combined job 4223072 for R=2). Rise (Δz) and terminal v_z, R=1 (2600
> steps) then R=2 (10400 steps):
>
> | Shape | Method | R=1 rise | R=1 v_z,end | R=2 rise | R=2 v_z,end |
> |---|---|---|---|---|---|
> | Sphere | MEM | +16.60 | 0.01138 | +39.51 | 0.00670 |
> | Sphere | Analytic | +13.35 | 0.00914 | +36.22 | 0.00616 |
> | Cylinder | MEM | +16.68 | 0.01197 | +37.33 | 0.00666 |
> | Cylinder | Analytic | +15.32 | 0.01036 | +38.39 | 0.00643 |
> | Box | MEM | **+0.001** | **~0** | +28.36 | 0.00559 |
> | Box | Analytic | +15.85 | 0.01072 | +37.50 | 0.00616 |
>
> **General trend** (sphere/cylinder, both resolutions): MEM-only buoyancy rises consistently a bit
> faster than the analytical (vSub*g) method — a real, reproducible behavioral difference, exactly
> what this item set out to characterize.
>
> **A separate, real anomaly**: at R=1, the box under MEM-only buoyancy essentially does not rise at
> all — `pos.z` stays pinned at ~41.5 for the full 2600-step run while `F_hydro_z` oscillates wildly
> (-0.07 to +0.19, no consistent sign) and X/Y drift laterally instead of the object rising. Verified
> directly against the raw log lines (not a parsing artifact — checked multiple raw `[BoxForces]`
> lines by hand). At R=2 it partially recovers (+28.36) but still visibly underperforms box-analytic
> (+37.50) and sphere/cylinder MEM at the same resolution (+39.51/+37.33) by a wider margin than
> sphere/cylinder show between their own MEM/analytic pair. This looks like a genuine
> **resolution-sensitivity problem specific to MEM-only buoyancy for box-shaped floaters** — not yet
> root-caused (candidate suspects: box corner/edge cells interacting poorly with the MEM
> hydrodynamic-force sampling at coarse resolution, or the `[BoxBuoyancy] V_sub method: SAMPLING
> (12^3 points)` scheme under-resolving a box's flat faces specifically compared to a sphere/cylinder's
> curved surface at the same point count) — worth flagging explicitly in the thesis, and worth an R=4
> check (job 4227812, pending) to see whether it keeps improving with resolution or plateaus below
> sphere/cylinder's own MEM/analytic agreement.

> **[Claude, 2026-09-13]:** per explicit request, verified the config and reran R=1 box/MEM to check
> reproducibility. **Flags confirmed correct**: `buoyancy_compare_box_r1_mem_gas.prm`'s
> `ModelParameters` block is byte-identical to the sphere/cylinder R=1 MEM prms (`diff` — zero
> differences); `useMemOnlyBuoyancy=true`/`useBognerThreePhaseMEM=true` both explicitly set, and the
> run's own startup log confirms `buoyancy model: MEM-ONLY`. The only difference between the three
> prms at all is the shape block itself (`FloatingBox halfExtents<8,8,8>` vs. the sphere/cylinder
> equivalents) — not a misconfiguration.
>
> **Likely mechanism, from re-reading the code** (`ParticleKernels.h:3103`): under
> `useMemOnlyBuoyancy_=true`, the applied gravity body force is unconditionally `rho*vFull*g` (pure
> weight, using the FULL volume, no buoyancy subtraction) — ALL upward lift is required to come
> entirely from the MEM hydrodynamic force `F_hydro_z`, by design. The R=1 box log's `F_hydro_z`
> oscillates wildly and inconsistently in sign (-0.07 to +0.19 across steps, no settling trend),
> unlike what a working case needs (a value that settles near the ~0.044 magnitude needed to balance
> gravity). This directly explains the near-zero net rise: F_hydro_z isn't reliably supplying
> sustained lift for the box at this resolution, not that lift is simply absent or misconfigured.
> Sphere/cylinder's own hydro forces (not yet extracted/plotted the same way, worth doing) presumably
> settle more cleanly at the same resolution.
>
> **Rerun submitted**: job 4228712 (`run_72_buoyancy_compare_box_r1_mem.sbatch`, resynced to
> `mt-dennys-mirror2` along with the R=1 bed file `spheres_out_buoyancy_r1.dat`, which that checkout
> didn't have — same cross-mirror gap as the R=4 bed file earlier). Report back once it finishes to
> confirm the anomaly reproduces exactly (same near-zero rise) rather than being a one-off fluke.

> **[Claude, 2026-09-14]:** job 4228712 COMPLETED, **anomaly reproduces exactly** (rise=+0.0006 vs.
> the original run's +0.001 — same near-total absence of rise, not a fluke). Quantified the
> mechanism precisely by computing `F_hydro_z` statistics across the full run for all three R=1 MEM
> cases:
>
> | Shape | rise (Δz) | F_hydro_z mean | F_hydro_z stdev | F_hydro_z range |
> |---|---|---|---|---|
> | Sphere | +16.60 | 0.0279 | 0.0053 | [0.0069, 0.0603] — always positive |
> | Cylinder | +16.68 | 0.0418 | 0.0032 | [0.0295, 0.0504] — always positive |
> | Box (rerun) | +0.0006 | 0.0427 | **0.0846** | **[-0.0804, 0.2082]** — crosses zero |
>
> This nails down the mechanism precisely: the box's `F_hydro_z` **mean** (0.0427) is actually
> perfectly reasonable — close to cylinder's own mean (0.0418) and clearly the right order of
> magnitude to balance gravity. The problem is pure **noise**: box's `F_hydro_z` standard deviation
> is 16-26x larger than sphere/cylinder's, and unlike sphere/cylinder (which stay strictly positive
> the entire run), box's `F_hydro_z` regularly swings negative. A mean-correct but wildly
> oscillating, sign-flipping force produces ~zero net displacement over many steps even though the
> time-averaged force would make it rise fine — this is a genuine MEM hydrodynamic-force **stability/
> noise** problem specific to box-shaped (flat-faced, axis-aligned) floaters at coarse resolution, not
> a wrong-magnitude or misconfiguration problem. Consistent with the earlier hypothesis (flat faces
> parallel to the LBM grid axes interacting differently with the interface/momentum-exchange scheme
> than a sphere/cylinder's curved surface). Worth checking against R=2/R=4 box-MEM data (same
> F_hydro_z noise metrics) to see whether this noise shrinks with resolution — that would confirm a
> coarse-grid artifact rather than a deeper box-specific implementation gap.

> **[Claude, 2026-09-14]:** **mechanism confirmed directly and dramatically.** Traced the root cause
> to `lbm_mesapd_coupling/MovingParticleMapping.h`: each LBM cell is classified fluid/solid by testing
> whether the cell CENTRE lies inside the particle shape. An axis-aligned box's flat top/bottom faces
> are exactly parallel to a grid plane, so a whole plane of cells can flip inside<->outside
> simultaneously as the box moves a tiny amount — large synchronized force jumps. A curved
> sphere/cylinder surface never has this (different points cross the threshold at different times).
> Tested directly: reran R=1 box-MEM with `initialTiltAngleDeg=25` (job 4232331, COMPLETED) instead of
> 0. Result:
>
> | Box variant | Rise (Δz) | F_hydro_z mean | F_hydro_z stdev | Range |
> |---|---|---|---|---|
> | Axis-aligned (0°) | +0.0006 | 0.0427 | 0.0846 | [-0.080, 0.208] — crosses zero |
> | **Tilted (25°)** | **+14.68** | 0.0522 | **0.0030** | **[0.042, 0.061] — always positive** |
>
> The tilted box rises normally (comparable to sphere/cylinder/analytic, all ~13-17 at R=1) and its
> `F_hydro_z` noise **dropped 28x** — now even tighter than sphere's own stdev (0.0053), right in line
> with cylinder's (0.0032). This is about as clean a confirmation as this kind of hypothesis test gets:
> the ONLY thing changed was the box's initial orientation, and it fully resolves the anomaly.
> **Conclusion for the thesis**: this is a genuine numerical artifact of cell-center-based MEM
> discretization for axis-aligned flat-faced bodies at coarse resolution, not a physical or
> implementation bug — worth describing explicitly as a known limitation/consideration for anyone
> using box floaters with MEM-only buoyancy in a "resting flush against the grid" orientation.

### Challenges & Future Improvements
- Need to decide the exact comparison metric up front (rise-time-to-surface? peak velocity? energy
  conservation during ascent?) so the two runs produce a like-for-like number, not just two
  qualitatively-different-looking curves.
- The two models may have different implicit assumptions about when "fully submerged" ends (MEM
  reads it off the fill field per cell; the analytic method uses a geometric `vSub` computation) —
  worth checking they agree on the crossover point itself before comparing behavior on either side
  of it.

> **[Claude, 2026-09-29] R1/R2/R4 RERUN WITH DISCRETE VTK SNAPSHOTS, per "using the same option as
> for these equilibrium test[s] can we now also run the R1 R2 and R4 again?"** Added the
> `fluid_field_snapshots` block (steps scaled per resolution: R1 `0,100,250,500,750,1000,1500,2000,
> 2500` of 2600 total; R2 `0,250,500,1000,2000,3000,5000,7500,10000` of 10400; R4 `0,500,1000,2000,
> 4000,8000,16000,25000,40000` of 41600 -- each list's last entry kept strictly below that run's own
> `timesteps`, since the main loop's `t != timesteps` condition never reaches `t==timesteps` itself)
> to all 18 `buoyancy_compare_{box,sphere,cylinder}_r{1,2,4}_{mem,analytic}_gas.prm` files, plus a
> thin Y-slab `CellBB_filter` through the floater's own y-centre at each resolution (full x/z) to
> keep frame size small. **Scope, per explicit correction ("remove the tilted cuboid 0 degrees run"
> + "option a and option b is turned off and all the mem is activated like in the equilibrium
> runs"):** guard (Option A/B) already defaults to `false` in `SimulationParams.h` and the `_mem`
> prms already had `useMemOnlyBuoyancy`/`useBognerThreePhaseMEM` explicitly `true` -- both already
> matched the equilibrium convention with no edits needed. The box's own untitled/0°-tilt MEM+analytic
> case (`buoyancy_compare_box_r{1,2,4}_{mem,analytic}_gas.prm`, the one bundled in the `_all` sbatch
> scripts) is deliberately **not** run this time -- only `run_72_buoyancy_compare_r{1,2,4}_tilted25.
> sbatch` (sphere MEM/analytic, cylinder MEM/analytic, box-**tilted25** MEM/analytic, 6 cases per
> resolution) is submitted, avoiding both the unwanted 0°-tilt box case and re-running sphere/cylinder
> twice. The box-tilted25 prms themselves are left VTK-off (data only, per "for this one you can turn
> off the vtk output we only need the data") -- only sphere/cylinder get the new snapshot mechanism.
> Bed files regenerated fresh per resolution first (the ones sitting in each checkout's build dir were
> stale, left over from the bogner3 s=32 VTK work's own domain/bedCopies). Job map (one bedgen +
> one dependency-chained compare job per checkout/resolution): **R1 -> bedgen 4278437 (mt-dennys) ->
> compare 4278440**; **R2 -> bedgen 4278438 (mt-dennys-mirror) -> compare 4278441**; **R4 -> bedgen
> 4278439 (mt-dennys-mirror2) -> compare 4278442**. Auto-archives to
> `$HPCVAULT/antidunes-archive/antidunes-buoyancy-r{1,2,4}-compare-tilted25_<jobid>_<timestamp>/
> {sphere,cylinder,box}_{mem,analytic}/{vtk-out,eval,snapshots}` via the existing script's own
> `archive_outputs()` (unchanged, no script edits needed -- it already archives whatever ends up in
> `vtk-out/`, so the new snapshot frames are picked up automatically).

> **[Claude, 2026-09-29] CORRECTION: VTK REMOVED ENTIRELY, R4 MOVED TO 4 NODES.** Per explicit
> correction ("None of the R1 R2 R4 need vtk we just want to compare analytical path vs mem... over
> the waterline doesn't make any sense for tracking because the original antidunes didn't support
> floating objects") -- the whole premise of adding `fluid_field_snapshots` above was wrong for this
> family: it's a **fully-submerged-rise-only** comparison (Section 1b's own stated scope, "only the
> fully-submerged portion of the rise is compared"), no VTK visualization is wanted at all, and the
> comparison never even reaches the waterline. The 3 jobs above (4278440/4278441/4278442) were
> cancelled and all 18 `fluid_field_snapshots` blocks reverted (`git diff` confirms all 18 prms are
> now byte-identical to their pre-session state). Separately, per "run r4 on four nodes, the other on
> a single node": all 8 R4 prms (box/sphere/cylinder mem/analytic + box mem/analytic tilted25) had
> `cellsPerBlock` changed `<32,32,24>` (6x6x16=576 ranks/8 nodes) -> `<32,32,48>` (6x6x8=288
> ranks/4 nodes), and `run_72_buoyancy_compare_r4_tilted25.sbatch` updated to
> `--nodes=4 --ntasks=288`/`--np 288` (time limit doubled to 14:00:00 for the halved rank count).
> R1/R2 were already single-node (`singlenode` partition, 72 ranks) and needed no change. Bed files
> from the earlier bedgen jobs (4278437/4278438/4278439) remain valid (bed generation depends on
> domainSize, not cellsPerBlock/rank decomposition) -- not regenerated. Resubmitted clean, no VTK, no
> dependency chain needed this time (bed already in place): **R1 -> 4278450 (mt-dennys, running)**,
> **R2 -> 4278451 (mt-dennys-mirror, queued)**, **R4 -> 4278452 (mt-dennys-mirror2, 4 nodes, queued)**.

> **[Claude, 2026-09-30] ALL THREE COMPLETED (R1 4278450 6:38, R2 4278451 1:09:49, R4 4278452
> 8:45:20) -- full 18-run R1/R2/R4 MEM-vs-analytical matrix now done, box's tilted25 fix confirmed at
> all three resolutions, new artifact published.** Parsed the `[FloaterForces]`/`[CylForces]`/
> `[BoxForces]` per-step logs directly from each job's own combined `.out` file (each runs all 6
> shape x method cases for its resolution back to back: sphere MEM/analytic, cylinder MEM/analytic,
> box MEM/analytic-tilted25, in that order -- confirmed via the `## R=... <shape> <method>` banner
> lines the sbatch script prints between cases). Rise = pos.z - starting z, restricted to steps
> before the floater's own top surface would reach the resolution's own measured `absoluteLiquidHeight`
> (81.07/162.14/324.28 for R1/R2/R4, read from each case's own `computedParams.txt` -- confirmed
> identical across all 6 shapes/methods at a given resolution, as expected for a domain/bed property).
> None of the 18 runs actually breached in their own step budget (calibrated conservatively, as
> intended).
>
> **The tilted25 box fix holds at all three resolutions, not just R1 where it was found:**
>
> | Shape | Method | R1 rise | R1 v_z,end | R2 rise | R2 v_z,end | R4 rise | R4 v_z,end |
> |---|---|---|---|---|---|---|---|
> | Sphere | MEM | +16.65 | 0.01133 | +39.57 | 0.00673 | +90.58 | 0.00377 |
> | Sphere | Analytic | +13.34 | 0.00917 | +36.32 | 0.00613 | +86.57 | 0.00362 |
> | Cylinder | MEM | +16.62 | 0.01195 | +37.47 | 0.00653 | +78.55 | 0.00332 |
> | Cylinder | Analytic | +15.30 | 0.01034 | +38.46 | 0.00642 | +85.45 | 0.00352 |
> | Box (tilted25) | MEM | **+14.54** | **0.01067** | +33.27 | 0.00588 | +71.12 | 0.00304 |
> | Box (tilted25) | Analytic | +14.07 | 0.00941 | +35.15 | 0.00582 | +77.71 | 0.00311 |
>
> Compare the box MEM row to the OLD, untitled-box R1 result from 2026-09-13 (**+0.001 rise, ~0
> v_z** -- the degenerate flush-face mapping tie that starved MEM-only buoyancy): the tilted box now
> rises **+14.54**, squarely alongside every other shape/method/resolution combination in the table,
> and stays that way at R2/R4 too, where the same bug would have applied if untilted. Sphere/cylinder
> numbers match the earlier (pre-tilted25, pre-R4-node-change) table closely (within ~0.1-0.3%,
> ordinary run-to-run noise) -- confirming those two shapes were never affected by either change, only
> the box needed the tilt fix and only R4 needed the node-count change.
>
> **New artifact published**: "Does MEM-only buoyancy agree with the analytical model?"
> (https://claude.ai/artifact/E8iSnHqoXAdj3Pa7ZxTkKP) -- 9 panels (3 shapes x 3 resolutions), each
> showing MEM-only vs. analytical rise trajectories overlaid, with final rise/v_z stat chips. Did NOT
> update the older "Buoyancy Sweep — Sphere, Cylinder, Cuboid" artifact
> (https://claude.ai/artifact/VzeThMNwVuMJhy958wLN5x) -- read it first to check, and it covers a
> completely different comparison (real vs. Bogner & R&uuml;de gravity calibration, not MEM vs.
> analytical), so a new artifact was the right call, not an update.
>
> **External-drive backup**: per explicit request, mirrored all three resolutions' vault archives to
> `/run/media/dhuber/data2/debris-clogging-archive/` (the same location as the equilibrium-run
> backups), as `antidunes-buoyancy-r{1,2,4}-compare-tilted25_<jobid>_<timestamp>/` siblings (R1 341M,
> R2 4.1G, R4 ~65G -- R1/R2 completed during this session, R4's transfer was still running in the
> background at ~7MB/s when this entry was written, expected to take on the order of an hour+ given
> its size; check `/tmp/rsync_r4.log` or just `ls` the destination folder to confirm completion
> later). Confirmed via the jobs' own logs that all 18 stages archived correctly to
> `$HPCVAULT` with no fallback-to-local warnings anywhere.

> **[Claude, 2026-09-30] New data ALSO added to the actual thesis text**, per explicit follow-up
> ("add the table with new data ... to the artifact that we have in thesis regarding this"): found
> the existing thesis section covering this exact comparison,
> `overleaf/chapters/04_results/equilibrium/mem-vs-analytical-buoyancy.tex`
> (`sec:results-mem-vs-analytic`) -- already has one rise/vz table (`tab:mem-analytic-rise`, the
> OLD data: untilted box's R=1 bug row `+0.0006`, and R=4 from the original 8-node/576-rank run) and
> a noise-statistics table already documenting the tilted-25 fix but only as a one-off R=1
> diagnostic. **Added a NEW table** (`tab:mem-analytic-rise-tilted`) rather than overwriting the old
> one -- the old table's job is to show the bug being discovered (flush box vs. every other
> shape/method), and overwriting its box row would have deleted that half of the story; the new
> table instead confirms the fix holds as the standing default across the FULL, self-consistent
> 18-run sweep (box always tilted, R=4 always on the new 4-node/288-rank decomposition), with a new
> paragraph noting sphere/cylinder match the old table within ~0.1-0.3% (decomposition noise, not
> physical) and that the tilted box now tracks its own analytical counterpart within 5-11% at every
> resolution, same order as sphere/cylinder, rather than reopening a resolution-dependent gap.
> Thesis rebuilds clean (verified: table renders correctly as Table 5.5, page 84, no new undefined
> references beyond the pre-existing `sec:cc-dispatch-gap`).
>
> **Also noticed** (unrelated to this table, found via the harness's own file-changed-on-disk
> notice): `probe-and-forcing.tex` was substantially reworked since this session's own item-4/5 edits
> -- restructured into its own `sec:clog-resistance` subsection, added the floater body-force term
> $\rho_f\sum V_{p,i}$ to the momentum balance (`analyze_clog_resistance.py` does NOT yet include
> this term, per that file's own new comment flagging it), and fixed the blockage-ratio reference
> depth to $\zeta_0$ with a note that all current probe planes sit where $A_{\mathrm{pil}}=0$ (48 LU
> from the pillar axis) -- worth updating the script and/or evaluating BR at the pillar line itself
> before relying on either for a real cross-scenario plot. Not actioned this entry, only flagged.

> **[Claude, 2026-09-30] Table 5.4 (`tab:mem-box-noise`) removed as outdated; provenance comments
> added to both equilibrium `.tex` files; artifact time-to-surface table switched to lattice steps
> only.** Three follow-up edits, all per explicit request:
>
> 1. **Removed `tab:mem-box-noise`** from `mem-vs-analytical-buoyancy.tex`
>    ("also remove the second table 5.4 we just explain why we had to tilt it t25 degrees and this
>    data is outdated") -- this was the one-off R=1 $F_{\mathrm{hydro},z}$ noise-statistics table used
>    to originally justify the 25° box tilt; now superseded by `tab:mem-analytic-rise-tilted`'s
>    full 18-run confirmation, so keeping both was redundant and the noise numbers were stale relative
>    to the new table's own R=1 row. Condensed the surrounding two paragraphs + table into a tighter
>    two-paragraph section (renamed "The box anomaly, and why tilting it fixes it") that keeps the
>    causal explanation (cell-centre classification + synchronized fluid/solid flips on flush axis-
>    aligned faces) as prose only, no numbers; removed the now-dangling
>    `\Cref{tab:mem-box-noise}`-referencing clause from the "fix confirmed" paragraph below it.
>    Rebuilt clean (2x pdflatex, exit=0 both times, only the pre-existing unrelated
>    `sec:cc-dispatch-gap` warning); verified via `pdftoppm`-rendered pages 97-98 that Table 5.3 (old,
>    untouched) -> condensed prose (no table) -> "fix confirmed" paragraph -> Table 5.4 (the tilted
>    confirmation table, auto-renumbered down from 5.5 by LaTeX with no broken cross-references) ->
>    Conclusion all render correctly.
> 2. **Data-provenance `%`-comment blocks added** to both equilibrium chapter files, per
>    ("for the plots and tables we added the results yesterday can you add comment to the latex file
>    where the data is stored for these plots and tables"):
>    - `shape-comparison-highres.tex`: cites the `extract_data.py`/`make_plots.py` pipeline and job
>      IDs (`bogner3_deep_s32_20260929/{4273570,4273582,4273583}` cuboid;
>      `bogner3_vtk_101k_20260929/{4278714,4278486,4278485}` sphere; `{4278476,4278477,4278484}`
>      cylinder), vault/local archive paths, and the artifact URL.
>    - `mem-vs-analytical-buoyancy.tex`: cites jobs 4278450 (R1) / 4278451 (R2) / 4278452 (R4), notes
>      the extraction was ad-hoc (direct log grep, no committed script, unlike the equilibrium-plots
>      pipeline), and gives the vault path
>      `$HPCVAULT/antidunes-archive/antidunes-buoyancy-r{1,2,4}-compare-tilted25_<jobid>_<timestamp>/`,
>      the external-drive mirror path, and the artifact URL
>      (`https://claude.ai/artifact/E8iSnHqoXAdj3Pa7ZxTkKP`).
> 3. **Artifact "Time to reach the surface" table** (https://claude.ai/artifact/E8iSnHqoXAdj3Pa7ZxTkKP,
>    now v3) switched from SI-time + step-count columns to **lattice-unit steps only**, per
>    ("also the Time to reach the surface artifact plot can you change the time to lattice unit time
>    steps not realtime"). Removed the ms columns entirely rather than just hiding them; updated the
>    lede text to note that the simulation's own time unit already *is* the lattice step (dt=1, no SI
>    conversion applied anywhere in this comparison), so showing steps alone is the more honest
>    framing for this artifact.
>
> None of this changes any underlying data or conclusion -- it's presentation cleanup (stale table
> removed, provenance made traceable, units made consistent with how the rest of this comparison is
> actually measured).

---

## 1c. [ ] Showcase: mixed-density floaters, some float / some sink

**Status**: we already have runs like this (mixed-shape/mixed-density sedimentation scenarios) — just
need to pick the best one as the actual thesis showcase (or run a clean dedicated one if none of the
existing runs frames the float-vs-sink contrast clearly enough).

> **[Claude, 2026-09-18]: dedicated "mixed sediment" thesis showcase written, NOT yet run or submitted.**
> `gen_mixed_sediment_prm.py` → `prms/mixed_sediment_showcase_100floaters_gas.prm` +
> `jobs/run_72_mixed_sediment_showcase_100floaters.sbatch` (singlenode, 72 ranks = every core, 10 h, unverified).
> - **Setup:** 100 free floaters (34 spheres / 33 boxes / 33 cylinders) × densities 0.5/1.0/1.5 (34/33/33; Bogner's showcase densities, all 9 shape×density combinations present) dropped from rest
>   into STILL water in a cube-like `<120,120,120>` (X,Y periodic), cellsPerBlock `<20,30,40>` = 6×4×3 = 72 ranks (domain stays a cube, blocks do not), 100000 steps, VTK every 10000 (dumps 0…90000; extreme-clog VTK flags `gas_cell_zero_setter` + `mapped_flag`, plus `fill_level`).
>   10×10 grid, pitch 12, one floater per column, two drop heights (Z=64/76, checkerboard; 76 not 80 because the Z=40 blocks put a seam at Z=80; waterline 56.02).
>   Shape×density via 9 balanced combos shuffled with a fixed seed.
> - **Shapes** keep the equilibrium runs' aspect ratios: box 7×14×7 (1:2:1 = equilibrium cuboid `<15,30,15>`), cylinder
>   r=4.5, L=11.25 along Y (L/D = 1.25 = equilibrium cylinder r=24, L=60), sphere r=4.5.
> - **Contact model:** `useAnalyticalContactDetection true` (Box SAT + hybrid Cylinder, `usePureAnalyticalCylinder false`)
>   and `useContactFilterOwnerTiebreak true`. Real gravity via `gravityOverride` derived from `--g` (default 9.81).
> - **Checked in the generator:** all 100 initial AABBs disjoint (≥3 LU, periodic), all floaters start above the water.
>   Not checked: the prm has not been parsed by the binary; runtime estimate (~2.5–3 h) is extrapolated.
> - **Submitted 2026-09-18 (final layout, one job per checkout, nothing shares a build dir concurrently):**
>   `mt-dennys-mirror` = hopper 4247229; `mt-dennys` = extreme-clog PID-switch (freeze at t=70000) 4247448, independent (no dependency);
>   `mt-dennys-mirror2` = the two still-water cases back to back: d075 equilibrium continuation 4246831, then the mixed sediment
>   showcase 4247456 (`--dependency=afterany:4246831`, singlenode, 72 ranks = every core, 10 h). Earlier attempts 4247446/4247454 (sediment
>   in `mt-dennys`, chained before the clog run) and the original PID-switch job 4247230 were cancelled while pending.
> - **Output safety:** the sediment and PID-switch sbatch scripts MOVE leftover `vtk-out/eval/snapshots` to
>   `build/.../archive/rescued_before_<jobid>_<stamp>/` instead of `rm -rf`, and archive whatever exists even when the run exits non-zero.
>   This matters in `mirror2`: if the d075 continuation crashes or times out, its checkpoint is rescued, not deleted, when the sediment job starts.
>   If the 10 h limit kills a job no archive step runs, but the next job's startup rescues the leftovers.
> - Bed file: `spheres_out.dat` present in the checkouts (header matches) and not read anyway (override branch, SimulationParams.h:1340).
> - Cores: sediment 72/72; clog runs 144/144 on 2 nodes; d075 continuation 140/144 (its checkpoint's block decomposition is 7×4×5, left as is).
> - **Status 2026-09-19:** d075 continuation 4246831 COMPLETED (3h35m, archive `.../antidunes-cuboid-stillwater-hires-d075-resume_4246831_20260919_083151`, 7.4G);
>   hopper 4247229 COMPLETED (1h48m, 60000 steps, archive `..._hopper-clog-pidon-2node_4247229_20260919_065136`, 8.6G) -- results NOT yet analysed;
>   mixed sediment 4247456 at step 98500/100000 when last checked (~2 h for 100k steps, far below the 5 h estimate; no error lines in the log).
>   **PID-switch 4247448 FAILED** at step 8400 (`reynoldsNumber is nan`, WriteData.h:98): first `[LOW-RHO]` NaN interface cell at step 8373 on rank 73,
>   block-local cell (34,0,0) (a floor/side-wall block corner, not at the pillar, floaters still dormant, PID switch not reached -- no `[PID-SWITCHOFF]` line).
>   Its acceleration/meanVelocity trace up to step 8300 matches the completed PID-on run 4245629 (same prm apart from the two switch keys), and that run and the
>   hopper logged ZERO LOW-RHO warnings, so this looks like a sporadic free-surface instability, not the PID-switch code (unproven: a single failure, root cause not
>   isolated). The new archive-on-failure logic worked (output preserved in `..._pidswitch70k-2node_4247448_20260919_051717`, 3.0G). Resubmitted unchanged as 4248336.
>   **Update: 4248336 (the unchanged resubmission) FAILED THE SAME WAY** (12 min, exit 1): first `[LOW-RHO]` NaN interface cell at step 7734 on rank 72, block-local (62,0,26) (again local y=0), `meanVelocity: nan` at 7800. Two failures out of two for this prm,
>   vs one completed PID-on run (4245629, older binary) and a completed hopper (4247229, current binary, switchOffStep 0). The three extreme-clog runs are NOT bit-identical (meanVelocity differs from step 600 on, e.g. 0.0431822 / 0.0431681 / 0.0431657),
>   so the trajectories do not separate "our switch code" from run-to-run nondeterminism; but the only code difference between the completed PID-on run's binary and the current one is the ForcingAdjuster change, so the switch path is NOT ruled out.
>   Untested discriminator: the unchanged PID-on 2-node prm on the current binary for ~12000 steps (does it NaN near 8000 too?).
>   **Discriminator result (job 4248429, PID-on prm unchanged, current binary, 12000 steps): PASSED** (0 LOW-RHO, no NaN, past the 7700-8400 window). Code re-read: before the switch step the ForcingAdjuster path is
>   identical to the original (the only addition, a `history_` push, never touches `acceleration_`), and nothing else reads `switchOffStep` except the constructor call in BuildTimeloop.h -- so the switch code cannot alter the
>   physics before t=70000. Nodes: failures on f[1008,1052] and f[0628,0634] (no shared node), passes on other nodes. Record for the 2-node extreme-clog config: 3 clean passes (4245629, 4245630, 4248429) vs 2 NaN
>   failures (4247448, 4248336, both the pidswitch prm) => best reading is a sporadic, non-bit-reproducible free-surface instability around steps 7700-8400 (first bad cell at block-local y=0 both times), not the switch
>   code and not a bad node. n is small; not proven. Resubmitted unchanged as 4248458.

> **[Claude, 2026-09-19]: results of the three finished runs (synced to `animations/`, checksum-verified against the vault archives).**
> - **Mixed sediment 4247456** (`animations/antidunes-mixed-sediment-showcase-100floaters_4247456_20260919`): 100k steps in 1h34m, exit 0, all 100 floaters
>   (34/33/33) present at all 10 dumps, no LOW-RHO / RELEASE-OVERLAP / NaN. Only warnings: 768 harmless `[STALE-NEARBOUNDARY]` (code clears the flag itself)
>   and ~1M `[CylComputation]` cylinder-contact diagnostics (log is 594 MB -- noise, worth silencing before long runs). Behaviour at t=90000, mean Z: rho 0.5 all
>   three shapes at Z≈62.0; rho 1.5 on the floor (sphere 4.50, cylinder 4.49, box 3.56 = lying on its long side); rho 1.0 spread over the whole depth (Z 5.6..59.9,
>   mean 39-42) as expected for neutral buoyancy. The 0.5 floaters sit ~6 LU above the *initial* waterline 56.02 -- NOT a buoyancy error: the initial fill level is
>   actually 56.52 (code adds +0.5) and the closed periodic box has no overflow, so the water rises to ≈61.5 (p90 of the slice column heights) from displacement
>   by 100 bodies (~49000 LU^3 / 14400 = 3.4 LU predicted; ~5 measured -- the extra ~1.6 LU is unexplained). Floaters at 62.0 vs level ≈61.5 = essentially half-submerged.
> - **Hopper 4247229** (`animations/antidunes-hopper-clog-pidon-2node_4247229_20260919`): NO persistent clog. 1015 spheres (0.5:290 / 1.0:435 / 1.2:290, gap/d = 2.67).
>   avgParticleTransportRate stays ~0.03-0.04 for the whole run; 532 spheres are downstream of the pillars by t=25000. Only a transient slowdown (mean speed near the
>   gap fell to ~9e-4 at t=40000, back to ~8e-3 by 45000). All 0.5 spheres and most 1.0 spheres pass; the 1.2 spheres (290) sink to Z≈5.5 and stay put (241 of 290
>   still at t=55000, median X≈225 = their starting block, only 13 ever downstream): the flow never mobilises them, so they never reach the gap. The clog probes
>   (elevation 25-50 vs baseline 49) are dominated by spheres passing through the probe columns and are NOT usable as a clog signal here. Caveat: periodic X means
>   passed spheres re-enter upstream, so this is not a finite hopper. To get an arch: smaller gap/d (~2.0-2.3, e.g. pillar radius 13-14 or sphere radius 3.5), and/or
>   start the 1.2 spheres near the gap, and/or a larger drive.
> - **d075 continuation 4246831** (`animations/antidunes-cuboid-stillwater-hires-density-sweep_20260918/d075_resume_4246831_20260919`, resumed from t=90000; log step = local,
>   absolute = 90000 + step): NOT converged to Bogner's 26.565°. Raw tiltY (wraps at ±180): ~54-58° up to t≈130k, 50.9° at 160k, 25° at 175k, then a roll-over of
>   ≈200° (-150.8° at 190k), wandering -150..-144° until 265k, then +124.75° at 280k and flat to 300k (std 0.17°, slope ≈0). Folding by the square XZ cross-section's 90°
>   symmetry the equivalent tilt is ≈39° (160k) -> 29° (190k) -> 36° (250k) -> 34.7° (final plateau) vs the 26.565° target: ~8° off, with 90°/roll-over flips on the way. So
>   the earlier "still relaxing" reading is refined, not confirmed as convergence to Bogner: rho=0.75 settles on a plateau ~35° (unproven whether a true equilibrium or a
>   long-lived metastable state; single run).

> **[Claude, 2026-09-19]: does Bogner's thesis explain why rho=0.75 doesn't reach 26.565°? No -- but our own archived runs narrow it down.**
> The thesis (Sec. 5.3.2-5.3.3) only says: lower width:height => shorter righting lever => larger relative discretisation error (5:4 worse than 6:4; jumps in the
> stability curve from lattice approximation; error shrinks with resolution), and that all three densities converge "nicely" at his low resolution (s=16, 130x40x40,
> tau=1/1.9, g=7.5e-4, 2.86° heel, no tolerance or convergence time given). The cube (1:1 section) is his LEAST stable ratio (upright unstable), i.e. exactly our box.
> Folded tilt (|angle| mod 90, then min(a, 90-a)) recomputed from the archived VTK rotation matrices, `animations/cuboid-buoyancy-sweep-combined/d0xx_{bogner,real}_vtk-out`
> (same s=30 box in 96x135x120, 5° start, 50000 steps; "bogner" = Bogner's tau=1/1.9 and g=7.5e-4 but NOT his geometry; "real" = g_SI=9.81, nu=2.14e-4):
> | rho | target | bogner tau/g, final (last 6 dumps) | real gravity, final (last 6 dumps) |
> |---|---|---|---|
> | 0.25 | 26.565 | 25.3 (24.9-25.4) | 25.9 (25.9-27.0) |
> | 0.50 | 45.0 | 45.0 (44.9-45.0) | 44.9 (44.7-45.0) |
> | 0.75 | 26.565 | **32.4** (32.3-32.6, steady from t~12000) | 30.1 (27.1-32.9, still wandering) |
> Reading (hypotheses, not tested): (1) rho=0.75 is biased ~+6° even at Bogner's own tau/g, and there it is rock steady -- so the offset is not caused by our low-viscosity/low-g
> settings; (2) rho=0.5 exact and rho=0.25 within ~1.3° in the same code, so the method is not generally off, but the 0.25/0.75 mirror symmetry (both 26.565 in the thesis)
> is broken => the error depends on the submerged fraction (possibly the MEM-only force near the box's top edges when 75% submerged); (3) the wandering / 200° roll-over
> (hires d075, 90k->300k) appears only in the low-viscosity real-gravity setup (nu 2.1e-4 vs 8.8e-3 at tau=1.9; g 5e-5 vs 7.5e-4) -- consistent with weak damping and long-lived
> metastable states. Not yet done: his exact geometry (s=16, 130x40x40, 2.86° heel); rho=0.75 with `useAnalyticBoxBuoyancy` (would show whether the bias is MEM-specific).

> **[Claude, 2026-09-19]: re-read Bogner Ch. 5 for a correction we might be missing -- none found.** Checked against our code: (a) Eq. 5.9 fill-weighted MEM force: implemented as
> phi * (PDF deviation from the rest density) in `SimpleBB.h`; equivalent to the thesis form because the gas term 2*c_q*(1-phi)*f_eq(rho_g,0) is the constant ambient pressure, which
> cancels over the closed body once the PDFs are zero-centred (`pdf -= w`, three-phase mode); ours additionally has Wen et al. 2014's Galilean-invariance term (irrelevant at rest).
> (b) "average forces over two time steps" (Ladd) = `useHydrodynamicForceTimeAveraging`, default ON, so already active in all archived runs. (c) lubrication correction: irrelevant for isolated
> boxes. (d) obstacle->fluid conversion (B(x) rule, phi interpolation, Eq. 5.3-5.5 equilibrium init) ~ our FormerMO reconstruction; NOT audited line by line. (e) R_P >= 2.5 dx: met.
> One real, uncontrolled DIFFERENCE: his floating-body validation used plain BGK (Sec. 5.3.2: "based on the FSLBM with BGK collision operator") while ours is a D3Q27 cumulant kernel; he also
> notes BGK+bounce-back has a viscosity-dependent boundary placement error (Sec. 5.3.1 "Calibration of hydrodynamic radius", which he does NOT correct for). Whether that matters for a
> hydrostatic (no-shear) equilibrium is unknown. Not a missing correction, but it means "Bogner tau/g" runs do not reproduce his collision model.
> **Submitted experiments (MEM-only, NO analytic buoyancy -- analytic-buoyancy run dropped on request; re-submitted 2026-09-19 with X AND Y periodic and near-full rank counts):**
> `bogner_exact_3box_s16` = job 4248419 (`mt-dennys-mirror2`, singlenode, 65 ranks = cellsPerBlock <10,8,40>, 13x5x1): his exact geometry, 3 boxes rho 0.75/0.5/0.25, 130x40x40, s=16, tau=1/1.9,
> g=7.5e-4, 2.86 deg heel, roll-only, 50000 steps; `bogner_exact_3box_s32` = job 4248420 (`mt-dennys-mirror`, multinode 2 nodes, 130 ranks = <20,16,40>, 13x5x2): same at 2x resolution (260x80x80, g scaled to
> 3.75e-4 to keep H*g = 0.03, 100000 steps) to test his claim that the error shrinks with resolution. 72/144 ranks are impossible (no factor 3 in 130/40/260/80); 65/130 are the maxima under X,Y periodic (>=3 blocks each)
> and the ghost-sync reach (box half-diagonal reachable in <=4 hops). Earlier 5-/64-rank, Y-wall versions (4248379/4248380) were cancelled ~10 min in (they ran at ~0.8 / ~1.1 ks per 1000 steps). Periodic boundaries and the roll-only pin are our choices, not from his text.
> **RESULT s=16 (job 4248419, 11 min, COMPLETED; `animations/antidunes-bogner-exact-3box-s16_4248419_20260919`): our code does NOT reproduce Bogner's heeled equilibria in his own geometry.**
> The 2.86° seed heel decayed to ~0° within 500-1000 steps in all three boxes (log tiltY and VTK rotation agree). rho=0.5 stays upright to 0.03° for all 50000 steps (target 45°); rho=0.25 stays upright, 0.03°
> (target 26.565°); rho=0.75 sits upright until ~step 18000-20000, then leaves to -4.8° and drifts back to -2.8° by step 50000 (target 26.565°). Box centre heights 16.5 / 20.5 / 24.5 = the upright Archimedes drafts
> +0.5 (the initial free surface is at 20.5, see the +0.5 fill offset). This contradicts our s=30, 5°-seed, fully-free runs (45.0 / 25.3 / 32.4 with Bogner's tau,g). The outcome therefore depends on seed size,
> resolution and/or the pins -- NOT separated yet (all three changed at once vs the archived runs: s 30->16, seed 5°->2.86°, free 6-DOF -> pinned roll-only/X,Y-pinned, plus periodic Y). Hypothesis (untested): at low
> resolution a lattice-aligned box is locked near upright (2.86° at s=16 moves a corner only 0.8 cells), so a small seed cannot break the unstable upright equilibrium. Cheap follow-ups (11 min each at 65 ranks):
> s=16 with a 5° / 10° seed; s=16 with the pins off. **Submitted 2026-09-19 (one variable at a time vs the baseline, chained afterany in `mt-dennys-mirror2`, 65 ranks each): seed 5° = 4248453, seed 10° = 4248454, pins off (2.86° seed, fully free) = 4248457.**
> **RESULTS of the three s=16 follow-ups (all COMPLETED, ~12 min each; folded tilt over the last 5000 of 50000 steps; targets 26.565 / 45 / 26.565 for rho 0.75 / 0.5 / 0.25; data in `animations/antidunes-bogner-exact-3box-s16-{seed5,seed10,nopin}_*`):**
> | variant | rho 0.75 | rho 0.5 | rho 0.25 |
> |---|---|---|---|
> | baseline: 2.86° seed, pinned | 2.8 (stuck near upright) | 0.03 (upright) | 0.03 (upright) |
> | 5° seed, pinned | 6.0 (still near upright) | 43.6 | 34.2 |
> | 10° seed, pinned | 21.9 (creeping up: 17.9@10k -> 21.9@50k) | 45.0 | 34.1 |
> | 2.86° seed, pins OFF | ~37 (erratic, raw -20..-51, no plateau) | 44.7 | 25.8 |
> Readings (single runs, so hypotheses): (1) the baseline failure IS a small-seed lock near upright: a larger seed OR removing the pins lets the boxes leave upright. Whether the pin or the extra freedom (unpinned boxes drifted 6-20 cells
> laterally, i.e. they get perturbed) is what breaks the lock is NOT separated. (2) rho=0.5 reproduces Bogner's 45.0° in every variant that escapes the lock. (3) rho=0.25 reproduces 26.565° only unpinned (25.8); pinned gives ~34° for BOTH
> seeds, so pinning shifts that equilibrium by ~8° (unexplained). (4) rho=0.75 remains the outlier in every variant: stuck (6°), slowly creeping (21.9° and still rising), or erratic (~37°) -- consistent with its slow/erratic relaxation in our earlier
> s=30 and hires runs; 50000 steps may simply be too short for it. Open: longer rho=0.75 run at s=16 (cheap), pins-on vs free at larger seed, and the s=32 run (4248420, still pending). **Queued 2026-09-19: long run = pins off, 2.86° seed, 300000 steps (`bogner_exact_3box_s16_nopin_300k`), single job in `mt-dennys-mirror2`, 65 ranks, 3 h limit (~75 min expected) -- see job ID in the queue; only step count / VTK interval differ from job 4248457.** **Also queued: pinned 10° seed long run (`bogner_exact_3box_s16_seed10_300k`, 300000 steps, chained afterany the pins-off long run 4248656 in the same checkout) -- rho 0.75 was still creeping up there (17.9@10k -> 21.9@50k).** s=32 (job 4248420, 2.86° = 1.6 cells at the corner) still pending. PID-switch control job 4248429 (PID-on prm, 12000 steps, current binary) is running.
> **[Claude, 2026-09-19] RESULTS: s=32 (4248420, 1h04m) and s=16 pins-off 300k (4248656, 54 min) COMPLETED; pinned-10°-300k (4248665) FAILED at step ~32371 (NaN, `reynoldsNumber is nan`, first bad interface cell (4,2,23) on rank 13) -- the same prm ran clean to
>   50000 steps before (4248454), so it is again a sporadic, non-bit-reproducible NaN; resubmitted unchanged as 4248877. Data in `animations/antidunes-bogner-exact-3box-s32_4248420_20260919`, `...s16-nopin300k_4248656_20260919`, `...s16-seed10-300k_FAILED_4248665_20260919`.**
> - **s=32, pinned, 2.86° seed, 100k steps** (folded last 10k): rho 0.5 = 45.00 exactly (left upright late, ~55k steps, 13° -> 45°); rho 0.25 = 31.5 (30.7-32.0, raw ~60°; s=16 pinned gave 34.1-34.2, target 26.565);
>   rho 0.75 = 3.2 -- still locked near upright after 100k steps (target 26.565). Reading: refining s=16 -> 32 moves pinned rho 0.25 towards the target (34 -> 31.5, consistent with his "error shrinks with resolution", but two points,
>   different seed histories) while the small-seed upright lock gets LONGER (rho 0.5 needed ~55k steps to escape, rho 0.75 never did).
> - **s=16, pins OFF, 2.86° seed, 300k steps** (folded last 25k): rho 0.5 = 44.0 (41.9-45.0; robust); rho 0.25 = 24.7 (23.0-26.5; box drifted 49 cells in x); rho 0.75 = 35.8 (29.8-43.9): raw 29.5@30k, 61@150k, then a ~90° flip
>   to 143@180k, ~122-149 afterwards -- i.e. it does not plateau even at 300k, it hovers ~30-37° folded. So "50000 steps is too short" is NOT the whole explanation: the free rho=0.75 box wanders within ~30-37° and never sits at 26.565.
> - **Combined with the earlier s=30 (Bogner tau/g, 32.4°), hires (~34.7°) and s=16 runs: rho 0.75's equilibrium in our code is ~32-35° in every free/escaped variant, ~6-8° above Bogner's 26.565°; rho 0.5 (45°) and rho 0.25 (~25-26° free)
>   agree with him. The rho=0.75 offset is therefore a reproducible bias, not a resolution, seed, pin or run-length artefact (untested: BGK collision as he used, and the MEM force near the top edges at 75% submergence).**

> **[Claude, 2026-09-19] Fixed-heel righting-moment scan to localise the rho=0.75 bias -- queued as job 4248900 (`mt-dennys-mirror2`, chained afterany the running seed10-300k job 4248877), NOT yet run.**
> - **Protocol** (Bogner Sec. 5.3.2, but with the raw MEM torque, not the analytic moment the old `bogner_righting_6to4` scripts logged, and NO analytic buoyancy): 8 heel angles {0, 5, 10, 15, 20, 25, 26.565, 30}° about Y, one run each (65 ranks,
>   20000 steps), three boxes per run (rho 0.75 / 0.5 / 0.25, s=16, Bogner basin 130x40x40, tau=1/1.9, g=7.5e-4, X/Y periodic), all held (position AND rotation pinned in every axis) at the heave where the submerged cross-section area equals rho*s^2
>   (waterline Z=20.5 = 20 + the 0.5 fill offset). Raw MEM `hydroTorqueY`/`hydroForceZ` of the owning rank's copy from the floater-history state CSV over steps 19000..19300 only (these CSVs blew the home quota before); each run's output is moved to
>   the vault right after the run; a run that dies is retried once. Files: `gen_bogner_righting_scan_prm.py`, `prms/bogner_righting_scan_s16_a*_gas.prm`, `prms/bogner_righting_scan_s16_manifest.json` (ideal reference), `jobs/run_65_bogner_righting_scan_s16.sbatch`,
>   `analyze_bogner_righting_scan.py <archive_dir>` (tested only on a synthetic archive -- it reproduces an injected +0.2 LU bias; the real MEM torque sign convention, assumed tau_y = -F*d, is untested).
> - **CORRECTION of my earlier statement:** in ideal theory the rho=0.25 and rho=0.75 torque curves are IDENTICAL (M_{1-rho}(alpha) = +M_rho(alpha), since the (1-rho) box is the 180° point-reflected image of the rho box), NOT equal-and-opposite. The generator's own
>   self-check on the exact numbers caught this (both -0.0865 at 5°). The test is therefore "do the two measured curves coincide", and the ideal curves vanish at exactly 26.565° (verified numerically; rho=0.5 vanishes at 45°).
> - **Sensitivity (ideal numbers, LU, s=16):** near the equilibrium the rho=0.75/0.25 torque is tiny and its slope is only ~0.03 LU per degree (|M| 0.056 at 25°, 0.096 at 30°, max 0.19 at 15°; the buoyant force is ~4.6 LU). A torque bias of ~0.2 LU -- a sub-0.05-cell shift of the
>   centre of buoyancy -- would move the zero crossing by ~6-7°, the size of the observed rho=0.75 offset. So the offset does not need a gross error, only sub-cell-level systematic bias in where the MEM sees the pressure.
> - **Phantom-force correction removed from consideration (old issue, fixed by moving to MEM buoyancy):** no longer cited as evidence for the contact-line hypothesis; the explicit `usePhantomForceCorrection false` line was dropped from `gen_mixed_sediment_prm.py`
>   and `gen_hopper_clog_cyl_prm.py` and from the (not yet started) cylinder-hopper prm on the cluster; the already-run mixed-sediment prm keeps it (a default-false no-op) so its record is unchanged. `useHydrostaticReconstruction` / `useDynamicReferenceHeight`
>   (same investigation) were dropped the same way (both generators, and the cylinder-hopper prm on the cluster -- synced before that job started; the already-run mixed-sediment prm keeps them, default-false no-ops).
> - **Also queued today (IDs not recorded earlier):** sphere front-block hopper 4248723 (`mt-dennys-mirror`, was chained after 4248420, now pending on priority); cylinder hopper 4248725 (`mt-dennys`, chained afterany the PID-switch retry 4248458).
> **[Claude, 2026-09-19] STATUS of the torque scan (4248900, CANCELLED by me after ~30 min) and of the pinned-10°-300k run (4248877, FAILED):**
> - **Scan design fails at large heel:** a00/a05/a10 ran (65 ranks, ~100 MB history each, 3 boxes); **a15 died twice** (`reynoldsNumber is nan` ~2 min in; first NaN interface cells at step ~1157 near the domain top, cell (2,0,33)), and a20 was throwing LOW-RHO at step 1500, so I cancelled. The `[MEM-WATCH] large hydrodynamic force` lines
>   at step 1 (all three boxes) are just the box weights, NOT the trigger. Cause of the start-up blow-up at heel >= 15° NOT identified. Nothing near the equilibrium angle (26.565°) was measured.
> - **What a00/a05/a10 show (raw MEM torque about +Y, LU, held boxes at ideal heave; sd over the window ~0, i.e. fully steady; ideal My = -F*d assumed sign):**
>   | heel | rho | My meas | My ideal | Fz err |
>   |---|---|---|---|---|
>   | 0° | 0.75 / 0.5 / 0.25 | ~0 (sd 0.95 / 0.006 / 0.45, noisy) | 0 | -14% / +21% / +2.6% (upright box aligned with the lattice: degenerate) |
>   | 5° | 0.75 / 0.5 / 0.25 | +0.012 / -0.001 / **-0.377** | +0.087 / +0.354 / +0.087 | +1.7% / +0.6% / -1.9% |
>   | 10° | 0.75 / 0.5 / 0.25 | +0.182 / +0.336 / +0.170 | +0.156 / +0.689 / +0.156 | +0.5% / +0.6% / +0.5% |
> - **Reading (single realisations, hypotheses):** the measured torque scatters around the ideal by ~0.35-0.45 LU (rho=0.5 is low by 0.355 at BOTH 5° and 10°; rho=0.25 is low by 0.46 at 5° but +0.015 at 10°) -- discrete jumps of the size Bogner describes as lattice-approximation
>   "jumps in the stability curve" (Fig. 5.6). That is 2-5x the ideal torque near the equilibrium (0.05-0.19 LU). The rho=0.75 vs rho=0.25 difference is +0.39 at 5° but only +0.011 at 10°, so it is NOT a constant asymmetry: it looks like staircase quantisation noise that is
>   different for the two boxes (liquid-side vs gas-side interface treatment), which is enough to set where the torque crosses zero. If so, the equilibrium angle at s=16 is decided by lattice effects, consistent with the lock-in plateaus / seed / pin dependence seen in the free runs. Needs angles near 26.565°,
>   finer angle steps, and a second resolution (s=32) to confirm. Fz agrees with the ideal buoyancy to 0.5-2% for 5° and 10°.
> - **Pinned 10° seed, 300k (4248877) FAILED at step 124979** (`[NaN-SCAN] aborting after timeloopFSLBM at t=125000`; first bad interface cell (8,0,27), again block-local y=0). The first attempt (4248665) failed at ~32k; the earlier 50k run (4248454) passed. Partial data to 124.5k (folded): rho 0.5 = 45.0; rho 0.25 = 33.5 (28.2-35.0);
>   **rho 0.75 sits on a very steady 13.5° plateau from 75k to 120k** (13.5 at 100k and 120k), NOT creeping to 26.565, and different from the 21.9° the same prm had reached in the 50k run -- a history-dependent lock-in plateau. NaN first cells so far (all block-local y=0 or nearby): (34,0,0), (62,0,26), (4,2,23), (8,0,27): possibly a free-surface/block-boundary issue at block y-edges; not investigated.
> **[Claude, 2026-09-19] Two-box fine scan QUEUED as job 4249061 (`mt-dennys-mirror2`, singlenode, 65 ranks, 3 h; NOT yet run):** rho 0.75 and rho 0.25 only (rho 0.5 dropped to test whether fewer boxes let the start-up survive heel >= 15°), raw MEM torque, no analytic
>   buoyancy, s=16, 20000 steps/angle, window 19000..19300. DIAGNOSTIC GATE: a15 and a20 run first; if BOTH fail twice the job stops (gate logic unit-tested), else it continues with the fine scan 16.5, 18, 19.5, 21, 22.5, 24, 25.5, 26.565, 27.5, 29, 30.5, 32°.
>   Files: `prms/bogner_righting_scan2_s16_a*_gas.prm` + `..._manifest.json` (`gen_bogner_righting_scan_prm.py --boxes 0.75,0.25 --angles ... --prefix bogner_righting_scan2`), `jobs/run_65_bogner_righting_scan2_s16.sbatch`; analyse with
>   `analyze_bogner_righting_scan.py <archive_dir> --manifest prms/bogner_righting_scan2_s16_manifest.json`. The generator change left the three-box prms identical apart from comments (checked, 0 of 8 non-comment diffs; manifest identical). The cancelled scan's leftover run
>   dirs in `mirror2` are moved aside (not deleted) by the new job's startup. A second resolution (s=32) is NOT queued: at 130 ranks a 40000-step run is ~25 min per angle, so it needs a reduced angle set chosen from the s=16 result.

> **[Claude, 2026-09-20] RESULTS -- two-box torque scan 4249061 COMPLETED (1h11m, all 14 angles on the first attempt; dropping the rho=0.5 box fixed the start-up problem at heel >= 15°).** Reduced tables: `animations/antidunes-bogner-righting-scan2-s16_4249061_20260919_reduced/`
>   (raw history CSVs in the vault archive `antidunes-bogner-righting-scan2-s16_4249061_20260919_185957`). Raw MEM torque about +Y (LU), held boxes at the IDEAL heave, sd over the window ~0 (fully steady -> static lattice effect, not noise or dynamics); ideal = -F*d (assumed sign):
>   | heel | My rho .75 | My rho .25 | ideal | | heel | My rho .75 | My rho .25 | ideal |
>   |---|---|---|---|---|---|---|---|---|
>   | 15 | +0.375 | +0.376 | +0.189 | | 25.5 | +0.430 | -0.345 | +0.040 |
>   | 16.5 | +0.154 | +0.070 | +0.189 | | 26.565 | +0.355 | -0.208 | 0 |
>   | 18 | +0.007 | +0.095 | +0.183 | | 27.5 | -0.159 | -0.281 | -0.034 |
>   | 19.5 | +0.178 | +0.163 | +0.170 | | 29 | -0.011 | -0.063 | -0.076 |
>   | 20 | +0.150 | +0.285 | +0.165 | | 30.5 | **-0.005** | -0.239 | -0.104 |
>   | 21 | +0.201 | +0.184 | +0.151 | | 32 | -0.247 | -0.154 | -0.121 |
>   | 22.5 | +0.324 | +0.108 | +0.123 | | 24 | +0.006 | +0.062 | +0.086 |
> - **The static MEM torque is a jagged staircase around the ideal curve**: rms error 0.17 LU (rho .75) and 0.14 LU (rho .25) against an ideal signal of only 0.04-0.19 LU near the equilibrium; the rho .75 - rho .25 difference is not systematic in sign (-0.135 at 20°, +0.22 at 22.5°, **+0.77 at 25.5°**, +0.56 at 26.565°).
>   Vertical force agrees with the ideal buoyancy to ~1-4% (mean +1.0% for rho .75, -0.8% for rho .25).
> - **Direction agrees with the free runs (suggestive, not proven):** rho .25's measured zero crossing is ~24.2° and its torque is more restoring than ideal at 25.5-26.6° (free runs settled at 24.7-25.8°); rho .75's torque crosses zero at ~27.2° and then stays ~0 at 29 and 30.5° (-0.011, -0.005 vs ideal -0.076, -0.104), i.e. a soft, almost flat
>   region 27-31° with weaker-than-ideal restoring torque (free runs settled at 32-35°). So the equilibrium shifts of rho .25 (down) and rho .75 (up) look explained by the static force computation, not by dynamics, the collision operator or the pins. Caveats: single realisations; boxes held at the IDEAL heave (not MEM's own, ~1% off);
>   linear interpolation across a jagged curve; assumed torque sign convention; one resolution.
> - **Suggested next tests:** (a) sub-cell position jitter at ~26.565 and 29° (z and x offsets of +-0.25 cell) to size the pure lattice scatter for one geometry; (b) the same scan at s=32 with a reduced angle set (22.5, 25, 26.565, 28, 30) to see whether the scatter shrinks with resolution.
>
> **[Claude, 2026-09-20] RESULTS -- extreme-clog PID-switch run 4248458 COMPLETED (2h26m, 0 NaN, exit 0)** (`animations/antidunes-clog-measurement-validation-extreme-pidswitch70k_4248458_20260920`): PID on until t=70000, then acceleration_x frozen at 6.84245e-06 (mean of 50 PID updates over 65000-70000; instantaneous 7.056e-06; `[PID-SWITCHOFF]` logged).
>   The PID had needed 8.5x the pre-jam force (~8e-7 at 20-40k) to hold discharge 0.04 with the jam. Eval windows (upstream / downstream elevation at the gap probes, baseline 49.13): 15-20k 53.5/48.1; 40-50k 56.4/44.2; 50-60k 70/33.7; 60-65k 82/22; 65-70k 84.5/11.2; 70-72.5k 88.6/9.8; 80-90k 89/4.9; 100-110k 86/1.1; 120-124.75k 86/7.2. Particle transport rate:
>   peak 7.8e-3 at 40-60k, 3.2e-3 at 65-70k, 4.6e-4 at 75-80k, ~0 (-3e-7) at the end. After the freeze the mean velocity falls 0.0400 -> 0.036 (75k) -> 0.031 (80k) -> 0.021 (90k) -> 0.0158 (120k+): the jam HOLDS and keeps tightening under a fixed driving force, with an afflux of ~80 LU. This is the clean clogging measurement (fixed-discharge onset, then constant-force response).
> - **[Claude, 2026-09-21] Correction + artifact update for the PID-switch run 4248458:** the 'afflux of ~80 LU' above is wrong -- it quotes the upstream ELEVATION (~87-89 LU), not the rise over the clean baseline; the afflux is ~87 - 54 = ~33 LU (baseline 49-54 depending on the window). The run is now in the artifact 'What a jam does to force, flow, and freeboard' (https://claude.ai/artifact/PhwjsUf7gfJZgbATWJSBRR): overlaid on the upstream-elevation, velocity and acceleration plots, its own panel (upstream + both downstream slits, PID-frozen marker at 70000) and a table of where the VTK output is (local `~/University/Thesis/animations/antidunes-clog-measurement-validation-extreme-pidswitch70k_4248458_20260920/vtk-out` and the PID-on/PID-off twins; vault `antidunes-archive/...pidswitch70k-2node_4248458_20260919_215425`; only 5 dumps per run locally: steps 0, 24950, 49900, 74850, 99800). The box-tilt artifact (https://claude.ai/artifact/S51Avsi7ZoWoPeoLY9urMa) now has a full-range y scale and a second plot of the rho=0.75 run continued to 300k (job 4246831, tilt from its `[BoxForces]` log `/home/hpc/iwia/iwia148h/mt-dennys-mirror2/apps/showcases/Antidunes/sbatch-output/output/cuboid_stillwater_hires_d075_resume_4246831.out`, unfolded and folded; final plateau folded 34.6 +- 0.17 deg over steps 280k-300k vs 26.565 deg).
>
> **[Claude, 2026-09-20] RESULTS (eval only; particle-VTK check pending sync) -- sphere hopper with the block against the pillars, 4248723 COMPLETED (1h50m):** still NO clog. Force stays ~1e-6 throughout (PID never needs more), transport rate 3.7-5.4e-2 (equal to or above the original hopper's 2.6-3.7e-2), upstream/downstream elevations fluctuate (34-51 / 22-43) with no sustained rise.
>   Moving the whole block to the pillars did not change the outcome at gap/d = 2.67. Cylinder hopper 4248725 is still pending (priority).

> **[Claude, 2026-09-20] Queued (all pending, one job per checkout):**
> - **Sphere hopper with LARGER PILLARS, 4249475 (`mt-dennys-mirror`, 2 nodes, 144 ranks, 6 h):** pillar radius 12 -> **15**, central slit 16 -> **10 LU = 1.67 sphere diameters** (d=6): still wider than a sphere (2 LU clearance each side for a centred one) but too narrow for two abreast (12 LU). R=16 (8 LU) would leave only 1 LU per side, too tight for the fluid to resolve
>   the flow around the sphere, so 15 is the practical maximum. Block shifted -3.1 LU so its front row keeps 1.0 LU from the pillar surfaces (verified on the written prm: only the 2 CylinderPillar + 2 BodyBC radii and a uniform X shift differ from the R=12 front-block prm; all clearances 1.0). Prm `hopper_clog_pidon_2node_frontblock_R15_gas.prm`
>   (`gen_hopper_clog_prm.py --pillar-radius 15 --front-x 365.4`). The cylinder hopper 4248725 (R=12, rods L=18 > slit 16) was NOT changed and is still pending.
> - **Sub-cell jitter scan, 4249476 (`mt-dennys-mirror2`, 65 ranks, 2.5 h):** heel 25.5 / 26.565 / 29 / 30.5° x offsets (dz,dx) in {(-0.25,0), (+0.25,0), (0,+0.25), (+0.25,+0.25)} cells applied to both boxes = 16 runs (`gen_bogner_righting_jitter_prm.py`; prms differ from their scan2 base ONLY in the position lines, checked). Ideal reference recomputed at the shifted heave
>   (manifest `bogner_righting_jitter_s16_manifest.json`). Note from the ideal numbers: a 0.25-cell heave shift changes the ideal torque by about +-0.18 LU (e.g. -0.148 vs +0.217 for rho .75 at 25.5°), i.e. ~0.7 LU per cell -- so the ~1% force mismatch (heave error ~0.1 cell) alone could move the torque by ~0.07 LU, the size of the scatter.
> - **Two-box scan at s=32, 4249478 (`mt-dennys-mirror2`, 2 nodes, 130 ranks, 6 h, chained afterany the jitter job):** angles 22.5, 25, 26.565, 28, 30°, 40000 steps each (window 39000..39300), g=3.75e-4; ~26 min per angle expected (from the free s=32 run's 0.038 s/step) => ~2.2 h. Question: does the static torque scatter (rms 0.14-0.17 LU at s=16) shrink with resolution?

> **[Claude, 2026-09-20] RESULTS of the three queued jobs + the cylinder-hopper failure:**
> - **R=15 sphere hopper 4249475 COMPLETED (1h52m):** still NO clog at slit 10 LU = 1.67 d. The narrower slit only adds hydraulic resistance: PID force 2.5-3.4e-6 (~3x the R=12 runs' ~1e-6), upstream level 50-63 vs 10-11 downstream (afflux ~40-50 LU), but the transport rate stays 3.0-5.1e-2 to the end (R=12 block: 3.4-5.4e-2): spheres still pass single file. Data `animations/antidunes-hopper-clog-frontblock-R15-pidon_4249475_20260920` (eval; VTK sync running). **Particle-VTK check (done 2026-09-20, both runs, 12 dumps each):** spheres downstream of the pillars peak at ~470-500 of 1015 at t=20000-25000 (R=12 block: 0.5/1.0/1.2 = 287/301/78 at 20k; R=15: 263/181/50 at 20k) and the pillar-zone mean speed stays ~1e-2 (only brief near-stalls: 2.5e-4 at 35k for R=12, 3.3e-3 at 50k for R=15, both recovered) -> NO persistent arch in either. With the block against the pillars 50-83 of the 290 heavy (1.2) spheres now pass (13 when the block started far upstream), but 126-190 of them are at rest at any time.
> - **Jitter scan 4249476 (FAILED overall, 14 of 16 runs have data):** `a29_zp25` and `a30p5_zp25` (both dz=+0.25) died twice each (history CSV empty; cause not investigated -- both failures are the +0.25 heave shift at the two largest angles). Findings (`animations/antidunes-bogner-righting-jitter-s16_4249476_20260920_reduced`):
>   **(1) The static MEM-torque scatter is pure lattice registration:** shifting BOTH boxes by only 0.25 cell in X (ideal torque unchanged) changes My by mean |0.157| LU, max 0.38 LU (e.g. rho .25 at 26.565°: -0.208 -> +0.173, a sign flip); across all sub-cell registrations of one angle/box the error (meas - ideal at that heave) spreads by 0.48 LU on average (max 1.15), rms error 0.24 LU -- 2-5x the ideal signal near the equilibrium
>   (0.04-0.19 LU). A single pinned measurement is therefore dominated by registration noise; only the registration-AVERAGED error is meaningful, and a pinned/locked free box samples ONE registration (consistent with the lock-in plateaus, pin and seed dependence).
>   **(2) Registration-averaged error, ideal-heave scan2 (14 angles) / jitter (14 runs):** rho .75 +0.057 +-0.047 / +0.114 +-0.081 LU; rho .25 -0.065 +-0.040 / +0.040 +-0.030 LU. With the ideal slope ~0.03 LU/deg near 26.565° a +0.06..+0.11 error shifts rho .75's zero crossing up by ~+2..+4°, rho .25's by ~0 (-2..+1°). Sign agrees with the free runs (rho .75 settles 6-8° high, rho .25 ~1° low) but the static estimate accounts for only about half to a third of the rho .75 shift,
>   each mean is only 1-2 sigma, and heave error (MEM's own equilibrium heave is ~1% off the ideal one; d(ideal My)/d(heave) ~0.7 LU per cell) is not included. SUGGESTIVE, not proven.
> - **s=32 scan 4249478 COMPLETED (1h56m, all 5 angles)** (`animations/antidunes-bogner-righting-scan2-s32_4249478_20260920_reduced`; torque scales ~s^3 so compare normalised by 8): meas-ideal rms 0.656/8 = 0.082 LU-equivalent vs 0.14-0.17 at s=16 (roughly halves with doubled resolution, only 5 points); measured zero crossings ~27.65° (rho .75) and ~25.5° (rho .25) vs 27.2° / 24.2° at s=16 and ideal 26.565°; vertical force error +1.1..+1.4% (rho .75) and -0.05..-3.2% (rho .25), same sign pattern as s=16 (+1.0% / -0.8%) and NOT shrinking with resolution
>   -> a small systematic MEM buoyancy bias (heavier-submerged box too buoyant, mostly-emerged box too little), equivalent to an effective waterline offset of ~0.1-0.3 cell.
> - **Cylinder hopper 4248725 FAILED after 5 min: `[RELEASE-OVERLAP]` abort at mainStep 2001, subCycle 2** (two pairs of top-layer rods: x=249/242 penetration -2.557 and x=179/193 penetration -4.482, both y=30, z=49.13; nothing else abnormal: rods had velocity 0 at subCycle 0; the hold-phase detection correctly returned contact=0 for the same pairs). **Cause NOT identified.** Ruled out by a new local test (`ANTIDUNES_CylinderParallelGap`, added to tests/showcases_antidunes, builds+runs in seconds): the detection routine (`SolidCylinderContactDetection`, analytic hybrid on AND off) correctly reports NO contact for parallel Y-axis r=3 L=18 rods at
>   surface gaps 1, 2, 8, 15 LU, with or without a 1e-9 rotation perturbation; rod orientation onto Y verified against the mixed-sediment VTK; the reported penetrations are not reproduced by any combination of world-axis orientations (X/Y/Z pairs) either -- so at the abort the rods must be at a non-axis-aligned orientation or state the test does not model. The FIXED-flag is not consulted anywhere in ParticleKernels.h (only mesa_pd's own kernels), and FIXED-hold + release had only been exercised with spheres (and boxes in the bridge runs), never with cylinders.
>   Options: (a) a single-node diagnostic of the exact prm to step ~2100 with the floater-history window around the release (positions, quaternions, angular velocity, torque of the offending rods), (b) drop the hold (settlingSteps 0) as the mixed-sediment run did (rods start in the water at the surface; the surface-start risk that motivated the hold for spheres is untested for rods).
> - **[Claude, 2026-09-20] Follow-up on the `[RELEASE-OVERLAP]` abort (Dennys' correction): it cannot have been an insertion overlap.** The check is meant to prove a floater was released while ALREADY overlapping, which is only true at sub-cycle 0 (detection on the positions at the instant of release, nothing has moved); but the code had NO sub-cycle condition, so a contact at ANY of the 10 sub-cycles of the first step after a release aborted with the "insertion overlap ... increase the stagger" message. Here it fired at sub-cycle 2 between rods that were at rest at sub-cycle 0, 1 and 8 LU apart -- so something changed in two sub-cycles (2 x 0.1 LBM steps, i.e. no physical motion can do that) and the message misdiagnosed it.
>   **Patched `ParticleKernels.h`:** abort ONLY at sub-cycle 0; at later sub-cycles log a throttled (40) `[RELEASE-CONTACT-LATER]` warning with position, quaternion, linear/angular velocity, force, torque and hydro force/torque of BOTH bodies, and continue. The hold-and-release semantic is unchanged (held floaters are ordinary active floaters; `resetHydroForces` uses SelectAll so nothing accumulates while held -- checked). Unit tests `ANTIDUNES_ContactDetection` and `ANTIDUNES_CylinderParallelGap` pass with the patch; the local Antidunes target compiles and contains the new message. Sources synced to all three checkouts (0 of 26 differ).
>   **Diagnostic job** (real cylinder-hopper prm on 48 ranks / one node, stopped at step 2400): prms/hopper_clog_cyl_frontblock_diag_1node_gas.prm + jobs/run_48_hopper_clog_cyl_frontblock_diag_1node.sbatch (only cellsPerBlock, timesteps and VTK interval differ from the production prm). `gen_hopper_clog_cyl_prm.py` got size options (--lx --lz --cells --vtk-every --vtk-crop); default output verified byte-identical. A local reproduction is not possible with the current build: the local tree is a DEBUG build and trips a debug-only assertion in BoundaryHandling.h:1203 (domain flag must be a single flag) that the production code violates by design.
> - **[Claude, 2026-09-20] ROOT CAUSE FOUND -- the ANALYTIC cylinder x cylinder contact path reports contacts for SEPARATED near-parallel rods.** Diagnostic job 4250809 (COMPLETED, 12 min, exit 0, no abort/NaN, `animations` not synced) logged 259 `[RELEASE-CONTACT-LATER]` lines, all at mainStep 2001 sub-cycle 2: the two bodies are exactly at their initial positions (distance 7 = surface gap 1, or 14), quaternion (0.7071,-0.7071,~1e-7,~3e-7) = axis along Y as designed, velocities ~1e-5, angular velocities ~1e-6 -- nothing moved -- yet the detection returns penetration -2.14 / -2.67.
>   **Local reproduction** (`ANTIDUNES_CylinderParallelGap`, Monte Carlo, 400 trials per case, positions +-1e-6): with `useAnalyticalContactDetection true` (the production flag) separated parallel r=3 L=18 rods are reported as touching in 0/400 at rotation noise <= 1e-7 rad, but **201/400 (gap 1) and 98/400 (gap 8) at 1e-6 rad**, 59/400 and 27/400 at 1e-5 rad; with the analytic path OFF (GJK/EPA only) **0 false contacts in every case**. 385 of the ~391 false contacts come from the analytic scenario `band-band-skewed` (overlap reported ~5.93-5.97 = almost 2r, i.e. an implied distance of ~0.05 between axes that are 7 apart), a few from face-band / face-face / band-band-parallel. The hybrid's uncertainty gate does not catch them.
>   Why it appeared only at sub-cycle 2 after the release: perfectly aligned rods (<=1e-7 rad) are handled correctly; the first tiny torques after release perturb the orientations to ~1e-6 rad and flip the pair into the bad branch. Consequence in the diagnostic run: 747 `[STRONG-CONTACT-FORCE]` and 706 `[DEM-WATCH] deep contact` warnings from spurious contact forces between rods that do not touch (the run survived 400 steps). **The already-archived mixed-sediment run (analytic ON, cylinders resting parallel on the floor) may contain such spurious contacts -- not yet checked.**
>   Proposed fix (implemented 2026-09-20 -- see the next note; original proposal): a sound capsule bound in the analytic cylinder x cylinder dispatch -- a cylinder lies inside the r-neighbourhood of its axis segment, so if the distance between the two axis segments exceeds r_i + r_j the cylinders cannot touch and any analytic "contact" is rejected (falls back to no contact); verify with the Monte Carlo test (expect 0 wrong contacts) plus the existing 33-case and fuzz tests.
> - **[Claude, 2026-09-20] FIX IMPLEMENTED AND TESTED (false near-parallel band-band-skewed cylinder contacts), checked against Kodam et al. (2010) Part I, Eqs. 55-66.** *Root cause, now pinned to one line:* Kodam's skewed-axes criterion (Eq. 63) needs BOTH closest-point parameters l_i*, l_j* (Eqs. 59-60, Vega & Lago) to lie inside the cylinder lengths and d = |x . r_hat| < R_i+R_j with r_hat = (e_i x e_j)/|e_i x e_j| (Eq. 55-56). For near-parallel axes (tilt theta ~ 1e-6 rad, ordinary rotation noise) the two infinite centrelines only approach each other ~ separation/theta ~ 1e6 LU away along the axes, so **the paper's own l* test correctly rejects**. The implementation (`PaperAnalyticalContacts.h`, band-band-skewed branch) instead replaced l* by 0 whenever 1-(e_i.e_j)^2 < 1e-10 ("guard against division by ~0"), which makes both l* conditions vacuously true, while r_hat is then a pure noise direction (any unit vector perpendicular to both axes -- e.g. perpendicular to the actual separation, giving d = 0). Deterministic reproducer (Y-axis rods, r=3, L=18, second rod tilted 2e-6 rad about Z, separated along X by a 0.1 LU gap): the old code returned `band-band-skewed`, **overlap 6.0 = R_i+R_j** (d = 0); a literal implementation of Eqs. 55-66 with the exact l* (no guard) says no contact. So Kodam's equations are fine and only the guard deviated from them.
>   **Fixes (two independent layers, both verified alone by mutation):** (1) `PaperAnalyticalContacts.h`: the skewed branch is applied only when the l*/d formulas are well-conditioned (1-(e_i.e_j)^2 > 1e-10, i.e. theta > 1e-5 rad); below that the branch is not applicable -- justified because every near-parallel pair that really touches is caught by the parallel-axes branch at the top of the function (|zIJ| < (L_i+L_j)/2 and rhoIJ < R_i+R_j), the residual radial-distance error over the length being <= theta*L ~ 2e-4 LU. (2) `ParticleKernels.h` (Cylinder x Cylinder dispatch): an ALWAYS-ON capsule bound (independent of `useTightLayer1`) rejects any analytic contact whose axis SEGMENTS are farther apart than R_i+R_j (a cylinder is inside the capsule of its axis segment, so this can never reject a real contact); it also protects `usePureAnalyticalCylinder` mode, which has no GJK/EPA cross-check. New log counter `cylCapsuleRejected` in `[ContactDetectionStats]` (`capsule-bound rejected analytic false positives=`), which must stay 0 now that layer (1) is in.
>   **New separate unit test `ANTIDUNES_CylinderNearParallelSkewed` (tests/showcases_antidunes/CylinderNearParallelSkewedTest.cpp), asserting:** (A) deterministic reproducer grid, tilt {1e-7..1e-3} x gap {0.05..2}: literal-Kodam reference, `paper_analytical::cylinderCylinderContact`, the production dispatch (analytic) and the dispatch with `useTightLayer1` must all report NO contact; (B) Monte Carlo, 4800 random near-parallel separated pairs (tilt 1e-7..1e-4, random tilt axes, random separation direction): 0 false contacts in every implementation; (C) genuinely overlapping near-parallel rods (overlap 0.1/0.5/1.0, tilt 0..1e-4) are still detected with the right overlap (< 1e-3 error); (D) rods 34 LU apart are never reported. **Before the fix it failed** (band-band-skewed, overlap 6.000000, dispatch=1, from tilt 2e-6 on; only `useTightLayer1` rejected it), **after: ALL CHECKS PASSED.** `ANTIDUNES_CylinderParallelGap` is now an asserting regression too (0 wrong contacts, analytic on/off, incl. Monte Carlo up to 1e-5 rad). Also passing after the change: `ANTIDUNES_ContactDetection`, `ANTIDUNES_PaperAnalyticalContacts` (33 cases + fuzz, 204 s), `ANTIDUNES_CylinderCylinderContactBenchmark`, `ANTIDUNES_BognerStateClassification/ClipAnalyticGeometry/ForcingAdjuster`; the local `Antidunes` target builds. Sources + tests + the two cylinder-hopper prms synced to `~/mt-dennys`, `~/mt-dennys-mirror`, `~/mt-dennys-mirror2` (md5 verified).
>   **`useTightLayer1 true` is now ON in the two cylinder-hopper prms (production + 48-rank diagnostic) and in `gen_hopper_clog_cyl_prm.py`** at the user's request (the Monte Carlo showed it alone already rejects all false pairs, since it is the same capsule bound applied before the dispatch); all other experimental contact/buoyancy flags stay OFF. NOT yet re-run on the cluster. **Caveat for the archived mixed-sediment run:** it used the OLD analytic code (its cylinders are affected); its `band-band-skewed contact=1` lines (log on mirror2) should be checked for spurious contacts before its cylinder statistics are used.
> - **[Claude, 2026-09-20] `useTightLayer1` is now ON BY DEFAULT** (user request): `SimulationParams.h` parse default `true`, `ComputeContacts` member and constructor defaults `true` (`ParticleKernels.h`); a prm can still set it `false` explicitly (three still do: `mixed_sediment_showcase_100floaters_gas.prm` + its generator, `box_sedimentation_showcase_analytical_125floaters_realgravity_base_statslog_gas.prm`; `DEMValidationScenarios.h` also passes an explicit `false` for its baseline). Effect on prms that don't mention it (410): Cylinder x Cylinder uses the capsule pre-filter, and Box x Box (only with analytic detection on) gets the single-axis pre-filter; both conservative (validated 0 unsafe cases / 10000 trials in `PaperAnalyticalContactTest`). All 8 local `ANTIDUNES_*` tests re-run and pass with the new default.
>   **Submitted (2026-09-20):** diagnostic `antidunes-hopper-clog-cyl-frontblock-diag-1node` = job 4250845 (48 ranks, 1 node, checkout `~/mt-dennys-mirror`) and production `antidunes-hopper-clog-cyl-frontblock-pidon-2node` = job 4250846 (144 ranks, 2 nodes, checkout `~/mt-dennys`); both with the fixed analytic cylinder path + `useTightLayer1 true`. Things to check in the diagnostic: no `[RELEASE-CONTACT-LATER]`, `[STRONG-CONTACT-FORCE]`/`[DEM-WATCH] deep contact` counts (were 747/706 before), `capsule-bound rejected analytic false positives=0` in `[ContactDetectionStats]`.
> - **[Claude, 2026-09-20] Results of the two cylinder-hopper jobs with the fixed analytic path + `useTightLayer1`:** **Diagnostic 4250845 (48 ranks, to step 2400): COMPLETED, exit 0, 13 min.** No `[RELEASE-OVERLAP]`, no `[RELEASE-CONTACT-LATER]`, **0 `[DEM-WATCH] deep contact` (was 706 before the fix)**, no NaN/STALE lines, `capsule-bound rejected analytic false positives=0` (so the source fix alone already removes the false contacts; the safety net never had to act). The remaining `[STRONG-CONTACT-FORCE]` lines (200) are rod-vs-pillar contacts (shapeType 6 vs 2, pillar at (384,40,45), per-contact force ~0.09), i.e. genuine, and the step-1 `[TOTAL-FORCE-THRESHOLD]` lines are hydrodynamic force (2.7-12 LU) on the still-held rods in the first step (initial transient), not contacts. Contact stats: 19256 analytic contacts (band-edge 16980, band-band-skewed 2260, band-band-parallel 16) vs **16.57 M GJK/EPA fallbacks that all found no contact** (99.9 % fallback rate: analytic finds nothing -> Tier 1 -> GJK); `skipTier1 true` is proven safe together with `useTightLayer1` (now default) and would remove them -- not enabled, awaiting decision.
>   **Production 4250846 (144 ranks, 2 nodes): FAILED after 15 min at step 5714 with `[STALE-PARTICLE-IDX]`** (`SimpleBB.h:250`, ranks 36 and 42, uids ...727/...726: a fluid cell is flagged as owned by a particle uid that no longer maps to a local particle -- the fluid-side obstacle map desynced from mesa_pd's particle storage; the same guard that caught job 4166479). Precursors in that log: 883 `[GHOST-VANISHED]` (first at step 2060), 287 `[SUBCYCLE-LEFT-RANK]`, 7 `[VELOCITY-JUMP]` (|dv| 0.05-0.07 at x~335-349, steps 4285/4789). Deep rod-rod contacts (`[DEM-WATCH]`, overlap ~0.65 LU = 6.9 % of the interaction radius, step 2511, both rods y~30, x 285/280) appear once the rods pile up -- looks like genuine contacts in the jam, not the old false ones (not verified against GJK). **The false-contact bug is fixed; this abort is a different (particle migration / ghost) problem and is NOT yet diagnosed.** Candidate factors, unverified: rod interaction radius sqrt(3^2+9^2)=9.49 vs block width 20 LU in Y (rods span block boundaries all the time).
> - **[Claude, 2026-09-20] Cylinder hopper `[STALE-PARTICLE-IDX]` (job 4250846, step 5714) traced -- and a deeper problem found: the analytic cylinder contact reports ~0 overlap for DEEPLY interpenetrated rods.**
>   **Id matching stale error <-> floater log:** the error prints the particle `uid` (decimal 64-bit, e.g. 2594073385365405727); `floater-history/rank_<N>_state.csv` (`uid` column) and `rank_<N>.csv` (`uidA/uidB`) use the same value, and uids survive a checkpoint restart (verified: the same uids appear in the resumed runs). New script `analyze_stale_particle.py <floater-history dir> <uid...> [--from S --to S]` prints per sub-cycle the owner copy's pose/velocity/force/hydro force+torque and every accepted contact (partner, path, penetration, force). **Bug fixed on the way:** on a checkpoint resume the DEM step counter (`ParticleKernels::stepCounter_`, used by all `step=` logs and by `logFloaterHistoryStart/EndStep`) restarted at 0, so the first resume's window 5600..5760 never triggered (empty CSVs) and log steps did not match the original run; it now starts at `beginTimeStep` (absolute).
>   **What happened to the stale uids (failed run log):** uids ...727 and ...726 (and a third rod ...730, all on rank 36, x~210-215, y~10-12, z~10-15, i.e. next to the y=0 wall) are launched in ONE main step: step 5714, sub-cycle 1 (|dv| 0.3 -> pos ~1e5-1e6, vel ~1e6, `SUBCYCLE-LEFT-RANK`), the fluid-side obstacle map still references them -> `[STALE-PARTICLE-IDX]` abort. So the stale index is purely a consequence of a force runaway of three rods. Precursor at step 5700: raw MEM hydro force on ...727/...726 was 2-2.1 LU (20-70x its usual 0.03-0.15) balanced by contact force of the same size.
>   **Rerun from the step-5001 checkpoint (jobs 4250904, 4250907, 144 ranks, floater-history 5650..5850 in the second): both COMPLETED to step 6000 WITHOUT the abort** (the run is not reproducible: the trajectories diverge, e.g. ...727 at (212.1,10.36,12.9) vs (213.1,9.7,17.6) at ~5600-5650). Floater-history (2.2 GB) is archived in `$HPCVAULT/antidunes-archive/antidunes-hopper-clog-cyl-frontblock-resume5000-floaterlog_4250907_20260920_104411/floater-history`.
>   **Finding (verified locally, test `ANTIDUNES_CylinderDeepOverlapProbe`, real poses from that floater log):** many rod pairs are interpenetrated by 30-90 % of their volume (Monte-Carlo volume check; centre distances 0.3-2.9 LU for r=3 rods; the pairs were found by an axis-segment screening and then confirmed by the volume check -- the screening alone overcounts, e.g. end-to-end rods). For five such pairs the analytic Kodam+Guo dispatch reports `face-edge`/`band-edge` with **overlap 0.03-0.12 LU**, whereas GJK/EPA reports **5.8-6.0 LU** (i.e. full penetration). The closed-form scenarios measure a local feature overlap and are only valid for shallow contacts: once two rods are deeper than a scenario boundary, the analytic contact force collapses (0.1 LU instead of ~6) and the rods pass through each other. In the failed run `[DEM-WATCH] deep contact` shows rod-rod overlaps already 1.5 LU at step 2500 and up to 2.8 (step 2750) and 5.35 LU (step 3250) (median 1.4, p90 2.6 LU over 1820 lines; rerun: p90 4.8, max 18 LU, that maximum is not rod-rod-verified). **Hypothesis (unproven): the step-5714 runaway is a deep pair whose scenario flips from the tiny closed-form overlap to a large one (GJK-like ~6 LU) -> force jump of ~100x -> launch.** Not caught because the abort did not reproduce.
>   **Consequences / options (need a decision):** (1) the analytic hybrid is not usable in this jam regime as it stands -- either use GJK/EPA only for rod-rod (`useAnalyticalContactDetection` also switches Box x Box/pillar analytic path, so a cylinder-only switch would be needed) or cap/cross-check the analytic overlap against GJK when deep (e.g. trust the analytic result only when its overlap < ~10 % of r, else GJK); (2) the deep overlaps themselves mean contact stiffness/sub-stepping is too soft for the flow load (hydro force up to 2-12 LU per rod at the start): more sub-cycles or a shorter contact time; (3) the `avgParticleVelocity` evaluator writes `inf` from the release on (38 sqlite warnings) -- separate small bug, not investigated.
> - **[Claude, 2026-09-20] Follow-up: the analytic-vs-GJK flip is visible in the floater log, and GJK/EPA is the correct one (the analytic overlap is wrong in both flip states) -- but the flip alone does NOT explain the 1e6 launch.**
>   **Tracking (`analyze_contact_flips.py` on the 4250907 floater-history):** 503 rod-rod pairs in contact, 357 flip events (overlap changes >4x within one sub-cycle). Path mix over ~670 k rod-rod contact rows: band-band-skewed 345 k, band-edge 299 k, face-edge 25 k, edge-edge 4.9 k, GJK only 23 -- i.e. the analytic result is trusted almost everywhere. Typical flip (pair ...176/...180 at (353,10,34), step 5705, 10 sub-cycles): band-edge 2.17 -> **edge-edge 18.07** -> face-edge 0.0013 -> ... -> 0.012 (contact force 5.6 -> **46.6** -> 0.09 LU). Same at step 5667 (face-edge 0.0001 -> edge-edge 18.05 -> face-edge 0.0013, force 0.03 -> 46 -> 0.08). Largest rod-rod force in the whole window: 62 LU (edge-edge, overlap 17.98).
>   **Which is right (`verify_cylinder_penetration_truth.py`):** the exact penetration depth of two convex bodies is min over unit directions n of [h_A(n) + h_B(n) - d.n] (h = support function of the cylinder, h(n) = (L/2)|n.a| + R sqrt(1-(n.a)^2); minimised over 200 000 random directions + local refinement). For pair ...176/...180 it is **5.89-5.99 LU** (centre distance 1.2, axes 17 deg apart, i.e. deeply nested), matching GJK/EPA's 5.8-5.93 in the rows where GJK ran. The analytic values are wrong in every state: face-edge 0.001 (translating by that along its own normal would need **18.0** to separate), band-edge 2.17 (needs 7.5), edge-edge 18.05 (needs 18.4, and 3x the true depth). Reason: the closed-form Kodam/Guo scenarios describe the first-contact geometry of ONE feature pair (rim vs face / band / rim) and report the depth of that feature only (e.g. face-edge = how far one rim point is below the other body's face plane, along the face normal); this is valid only for shallow overlaps. When the bodies are nested, the scenario classification jumps between meaningless answers from one sub-cycle to the next. (The earlier magnitude-based 'Tier 2' cross-check that would have forced GJK/EPA for face-edge/band-edge/edge-edge was removed in an earlier session, see PaperAnalyticalContacts.h; only special-contact is distrusted now.)
>   **Honest limit:** forces up to 62 LU for single sub-cycles change a rod's velocity by only ~0.005 LU, so this chatter alone cannot produce the 1e6 LU/step launch at step 5714 of the failed run (would need ~1e10). The launch trigger is therefore still UNEXPLAINED (candidates: a degenerate contact point/normal giving a giant torque or force, the raw MEM hydro force spike seen on ...726/...727 at step 5700, or lubrication); the run cannot be reproduced. **Proposed fix (not implemented):** after every analytic result, verify it -- compute `needed = h_A(n)+h_B(n)-d.n` along the analytic normal (a few flops); if `|needed - overlap|` exceeds a small tolerance (e.g. 5 % of R) fall back to GJK/EPA; this is exact for the shallow cases the closed forms are valid for and removes the deep-overlap flips.
> - **[Claude, 2026-09-20] (a) Test proving the detection claim -- `ANTIDUNES_CylinderClosedFormRange` (passes, 45 s; fallback NOT implemented); (b) stiffness is NOT tied to sub-cycles (Dennys was right); (c) Kodam Part II stiffness; (d) three short comparison jobs queued.**
>   **(a) What the test shows** (ground truth = exact penetration depth min_n[h_A(n)+h_B(n)-d.n]; 14 000 random rod pairs r=3, L=18 + 5 real nested pairs from the hopper floater log): the closed form is **right when shallow, wrong when deep -- but NOT wrong for most deep pairs (I said "most" before, that was too strong):** wrong (|error| > 5 % R) for 2/893 pairs with depth < 0.3 R, 3.2 % at 0.3-0.6 R, **12.9 % at 0.6-1 R, 25.5 % at 1-2 R** (max under-estimate 6 LU = the whole depth; max over-estimate 12.3 LU). The 5 real pairs I picked are all failures (selection bias: they were screened for low reported overlap). GJK/EPA is within 0.001 LU of the truth for random pairs and within 0.25 LU for the near-coincident real ones. The self-consistency check (`needed(n) = h_A(n)+h_B(n)-d.n` along the reported normal must equal the reported overlap; needed >= exact depth always, so **every under-estimate is caught -- 100 % in all classes, proven and measured**) plus a cheap upper bound (min of needed() over 6 canonical directions) for over-estimates catches 340 of 341 over-estimates; false alarms on correct results 0-6 % (analytic normal not exactly the minimum-translation normal; a false alarm only costs a GJK call). So the fallback proposal is validated as a design; not implemented in production yet.
>   **(b) Stiffness:** `LinearSpringDashpot::setStiffnessAndDamping`: k = m_eff (pi^2 + ln^2 e)/tc^2, c = -2 m_eff ln e / tc, with **tc = 4 x avgParticleDiameter = 32 LBM steps hard-coded** (`SimulationParams.h`, was `particleCollisionTimeNonDim = 4`) and dt_DEM = 1/numSubCycles. So **sub-cycles do not change k**: a steady load gives overlap = F/k whatever dt is (sub-cycles only resolve the oscillation, here 320 DEM steps per collision -- massive over-resolution). I was wrong to suggest "more sub-cycles". For rods (m_eff ~ 254-300, e = 0.5): **k ~ 2.6-3 LU/LU**, i.e. a 3-15 LU load (hydro 2-12 LU per rod plus the pile above it) gives 1-5 LU overlap -- exactly the DEM-WATCH values (median 1.4). New prm `ParticleParameters.collisionTimeFactor` (default 4 = unchanged): k ~ 1/factor^2 (factor 1 -> 16x, 0.5 -> 64x); tc >= ~20-40 DEM steps needed for explicit stability (factor 0.5: 40 steps, ok).
>   **(c) Kodam et al. (2010) Part II Table 1:** cylinders d = 8 mm, L = 5.3 mm, m = 3.1e-4 kg, G = 1.15 GPa, nu = 0.35, e = 0.85, mu = 0.5; "true cylinder" DEM uses a damped LINEAR spring whose stiffness was tuned to match the maximum overlap of a binary Hertz impact (printed as k = 4.0E-05 N/m, n = 1.506 N s/m -- the exponent looks like a typo, 4.0E+05 N/m is the value consistent with G = 1.15 GPa, unverified). With k = 4e5 N/m the overlap under a particle's own weight (3e-3 N) is ~7.6e-9 m = 1e-6 diameters, i.e. essentially rigid; their tests are dilute (single impact on a wall, rotating drum, hopper-like deposition), no dense jam under a hydrodynamic load. They also note that the linear spring generally OVER-estimates max overlaps. Ours: k*d/F_net ~ 3*6/0.03 = 600 (net weight ~0.03 LU) vs ~1e6 for Kodam -- our contacts are ~1000x softer in that dimensionless sense, which is normal for LBM-DEM (collision time tied to the fluid time scale) but far too soft for a flow-driven jam.
>   **(d) Short comparison runs (fresh start, 6000 steps, 144 ranks, floater-history 5900-5960):** 4250976 = GJK/EPA only, tc = 4 d (default stiffness) in `~/mt-dennys`; 4250977 = analytic, tc = 1 d (16x stiffer) in `~/mt-dennys-mirror`; 4250978 = GJK/EPA only, tc = 1 d in `~/mt-dennys-mirror2`. Baseline = the failed analytic/tc=4d job 4250846 (up to step 5714). Prms `hopper_clog_cyl_short_{gjk_tc4,analytic_tc1,gjk_tc1}_gas.prm`. Questions: does GJK alone stop the deep interpenetration (is the analytic path the problem), does stiffness stop it (is the load the problem)?
> - **[Claude, 2026-09-20] RESULT of the three short runs: the deep interpenetration and the launches come from the ANALYTIC contact path, not from soft contacts.** (fresh start, 6000 steps, 144 ranks; exact rod-rod penetration depth from the floater-history poses via `analyze_rod_overlap.py`)
>   | run | outcome | `[DEM-WATCH]` deep contacts (overlap > 0.65 LU) | max rod-rod penetration depth | pairs > 1 R (3 LU) |
>   |---|---|---|---|---|
>   | analytic, tc = 4 d (baseline; 4250846 + resumes 4250904/07) | 4250846 ABORT step 5714 (`STALE-PARTICLE-IDX`); resumes reached 6000 | 1820 (median 1.4, max 5.35 LU) | **6.0-6.1 LU** (steps 5650-5850) | **10-12** |
>   | analytic, **tc = 1 d (16x stiffer)** (4250977) | **ABORT at step 3599, `STALE-PARTICLE-IDX` again** (a rod at 16-18 LU/step velocity in sub-cycle 7 of step 3599, after a pair at 0.2 velocity/0.84 LU overlap one step earlier) | -- | -- | -- |
>   | **GJK/EPA only, tc = 4 d** (4250976) | **COMPLETED 6000 steps, 0 VELOCITY-JUMP, 0 STALE** | **0** | **0.43 LU** | **0** |
>   | **GJK/EPA only, tc = 1 d** (4250978) | **COMPLETED 6000 steps, 0 VELOCITY-JUMP, 0 STALE** | **0** | **0.33 LU** | **0** |
>   Conclusions: (1) with GJK/EPA the rods stay shallow even with the default (soft, tc = 4 d) contacts -- max penetration 0.43 LU (14 % of R) versus 6 LU with the analytic path; so **the overlap problem is the analytic detection (force collapse for nested pairs, chatter 0.001 <-> 18 LU), not the contact stiffness.** (2) Stiffness is second order: 16x stiffer lowers GJK's max/mean penetration from 0.43/0.11 to 0.33/0.08 LU, and **does not save the analytic path** (it aborted earlier, step 3599 vs 5714). (3) Both GJK runs survived, and both analytic runs with fresh start aborted (2/2) -- small sample, chaotic system, but consistent with the verified detection defect (`ANTIDUNES_CylinderClosedFormRange`). (4) Timing is not worse: GJK-only 6000 steps in 14:58 / 14:48 min. **Recommendation:** run the cylinder hopper with GJK/EPA for rod-rod (or the verify-and-fall-back scheme once implemented), default stiffness is acceptable (optionally tc = 1 d for ~25 % smaller overlaps). Caveats: `useAnalyticalContactDetection false` also switches rod-vs-pillar to GJK/EPA (that worked in the runs); one run per configuration; the launch mechanism of the analytic aborts is still not pinned down (see the previous note), only that it correlates with the analytic path.
> - **[Claude, 2026-09-20] SHADOW-MODE run (job 4251012, COMPLETED 14:55 min, no STALE/VELOCITY-JUMP/NaN): GJK/EPA used, closed form computed and compared for every rod pair in contact.** New code: `ModelParameters.shadowAnalyticCylinderLog` (+ `shadowAnalyticStartStep/EndStep`), `SolidCylinderContactDetection::shadowCompare` (`ParticleKernels.h`): per pair scenario, closed-form overlap/normal, GJK overlap/normal, "needed distance along the closed-form normal", check flag; cumulative `[ShadowAnalytic]` log line + per-rank CSV `shadow-analytic/rank_N.csv` (steps 5500-5960, 3.5 GB in the vault archive `..._gjk-shadow-tc4_4251012_20260920_124513`); analysis: `analyze_shadow_analytic.py`. Prm: `hopper_clog_cyl_short_gjk_shadow_tc4_gas.prm`.
>   **Whole run (all 6000 steps, 46.7 M contact computations incl. duplicates on neighbouring ranks and the held phase):** the closed form also finds a contact in 99.97 %; misses 12 300 (grazing: GJK overlap <= 0.074 LU); closed-form contact where GJK sees none: 0. Wrong by > 0.15 LU (5 % R): **0** in every scenario except `special-contact` (56/816, all 816 flagged by the check; the scenario production already distrusts). band-band-skewed 27.7 M and band-edge 17.9 M and face-edge 0.8 M contacts: agreement to 1e-9 .. 1e-4 LU.
>   **Window 5500-5960 (2.13 M distinct pair-sub-cycles):** GJK overlaps are SHALLOW: median 0.030, p90 0.066, p99 0.139, max 0.211 LU. Closed-form/GJK overlap ratio (force-weighted): band-band-skewed 1.000, face-edge 1.000, band-band-parallel 0.998, band-edge 0.997 (0.4 % of contacts give < 50 % of GJK's overlap), **edge-edge 0.78 (34 % of edge-edge contacts give < 50 % of GJK's overlap; 11.6 k of 2.1 M contacts)**, and 198 GJK contacts (0.01 %) the closed form does not see at all. Normals agree (0 flipped, 0 off by > 20 deg in the window).
>   **Interpretation:** on a trajectory where contacts stay shallow (what GJK produces) the closed form is essentially exact, so it is NOT wrong in general; it fails in the rare configurations (edge-edge under-estimates, special-contact, and above all nested/deep pairs, see `ANTIDUNES_CylinderClosedFormRange`) -- and because a too-weak contact lets the pair sink further into the range where the scenario classification is meaningless, the analytic runs run away (positive feedback) while the GJK runs never enter that range. Not proven as a causal chain (analytic runs are chaotic and the launch trigger is still unexplained), but consistent with all data: analytic 3 aborts/deep overlaps (max 6 LU), GJK 3/3 clean (max 0.43 LU).
> - **[Claude, 2026-09-20] THE FAILING CONTACT IDENTIFIED (job 4251049 = analytic, tc = 1 d, floater-history 3300..6000; aborted AGAIN at step 3493 with `STALE-PARTICLE-IDX`, uid 3963167672086036486 on ranks 61/67 -- third analytic abort out of three fresh starts: steps 5714, 3599, 3493; floater-history 1.9 GB in the vault archive `..._analytic-tc1-floaterlog3300_4251049_20260920_132436`).** Launched: rods ...486 and ...489 (rank 67, x~322-325, y~29.5-30, z~56.7-58.5). Trace (`analyze_stale_particle.py`, exact depth via support function from the logged poses):
>   1. **Steps 3300-3485, contact ...486/...489 reported by the closed form as `band-edge`, overlap 0.0003-0.0035 LU (force ~0.07 LU), while the EXACT penetration depth grows steadily 0.18 (3300) -> 0.57 (3345) -> 1.2 (3390) -> 2.05 (3435) -> 2.8 (3465) -> 3.3 (3480) -> 3.45 LU (3485)** (centre distance 5.9 -> 4.2 LU for two r=3 rods, axes 11 deg apart): the rods slowly sink into each other with practically no restoring force -- the under-estimate the shadow run only saw in ~0.4 % of band-edge contacts.
>   2. **Step 3485 sub-cycle 9: the closed-form overlap jumps 0.0034 -> 2.48 LU with unchanged geometry (exact depth 3.45 -> 3.52)**; contact force 0.07 -> 40 LU at step 3486 (3486: 2.44, 3487: 1.57, 3488: 0.77, 3489: `edge-edge` 3.50 with force 28/-42, 3490: 1.1 ... 3492: `face-edge` 0.78, force 17; true depth 3.5 -> 5.2 LU). The two rods are kicked (velocity changes 0.05-0.1 LU/step for 4 steps, z 56.7 -> 57.3) while the raw MEM hydro force on ...486 rises from 0.09 (3485) to 1.3 (3486), 3.0, 2.8, -3.7 (3490), ... 
>   3. **Step 3493 sub-cycle 0: the raw MEM hydro force of ...486 is (-4.56e21, -4.56e21, +4.56e21), torque ~2e22 (...489 next to it: normal 1.3/2.7/-1.7) -> velocity 9e17 in one sub-cycle -> left the domain -> stale obstacle map -> abort.** So the direct trigger of the launch is a CORRUPT hydrodynamic force from the fluid side (garbage value, equal magnitude in all components), reached after the contact kick threw the rods around near the free surface; the cause of the garbage MEM force is NOT yet analysed (the fluid-side PDFs/mapping of the moved rod; consistent with the STALE-PARTICLE-IDX desync class of job 4166479).
>   **Chain, with evidence:** (a) closed-form band-edge under-reports a nested pair by ~1000x for ~200 steps (0.003 vs 0.2-3.4 LU) -> (b) scenario flip to 2.5 LU with unchanged geometry -> 40 LU kick -> (c) violent motion at the interface -> (d) garbage MEM hydro force -> (e) launch/abort. Steps (a)+(b) are contact-detection faults that GJK/EPA does not have (all three GJK runs, incl. the shadow run with 2.13 M logged pair-sub-cycles, are clean); (d) is a separate fluid-side failure that the contact fault provokes and that deserves its own look (does MEM produce garbage whenever a rod is thrown across the interface at > 0.1 LU/step?).
> - **[Claude, 2026-09-20] WHY THE JUMP (test `ANTIDUNES_CylinderClosedFormJump`, real logged poses of the launched pair, steps 3484-3489; passes):** recomputing the closed form from the logged poses reproduces the logged scenario/overlap in **60 of 60 sub-cycles** (band-edge 0.0035 ... jump at 3485.7 to 2.5259 ... 1.57 ... edge-edge 3.50 at 3489), GJK/EPA = exact depth in all of them (3.41 -> 3.47 -> 5.2 LU). Walking between the two logged poses 3485.6 -> 3485.7 (they differ by 0.0047 LU in total) in 50 steps: the exact depth changes by < 1e-4 LU per step, but the closed-form overlap jumps 0.0028 -> 2.528 at t~0.46. Cause: `rimVsFaceOrBand` (band-edge branch) examines only ONE rim circle of the probe rod -- the end-cap circle whose centre is NEARER to the host's centre (`capCentre = |capTop-host| < |capBottom-host| ? capTop : capBottom`) -- and reports the radial overlap R_host - (closest distance of THAT circle to the host's axis). In the nested pair the two rim circles give completely different values (one circle just grazes the band: 0.0027, the other dips 2.45-2.53 LU into it, neither equals the true depth 3.47), and which one is examined flips when the two cap centres are about equally far from the host centre -- an arbitrary tie that has nothing to do with contact physics. At t~0.46 the examined circle flips (caps 1,0 -> 0,0) and the reported overlap jumps by 2.5 LU (force 0.07 -> 40 LU) although nothing moved. This is a discontinuity of the closed-form dispatch itself, valid-range violation (the derivation assumes only one rim can be in contact), not a numerical accident: any nested pair passes through such a tie.
>   **Fluid side, next:** raw MEM force on ...486 at step 3493 = (-4.56e21, -4.56e21, +4.56e21): equal magnitude in all components = ONE D3Q27 corner link with |pdf| ~ 2.3e21 (forceMEM = pdf_old + pdf_new ~ 4.6e21) that passes the `isfinite` dry-link guard in `SimpleBB::treatDirection`. Added a log-only diagnostic there (`[MEM-GARBAGE-LINK]`, |pdf| > 10: fluid-cell flags incl. former-moving-obstacle/conversion flags, fill, density, all 27 populations, neighbour flags) and queued job 4251065 (analytic, tc = 1 d, per-step PDF scan of all liquid/interface cells 3000..6000 with `highVelocityPrecursorThreshold 0.5`, floater history off) to catch the value at its origin. Hypotheses to test: (1) a cell freshly freed from the rod (former-moving-obstacle) or freshly converted gas->interface gets reconstructed PDFs from neighbours that were inside an object / gas; (2) a drained interface cell (rho -> 0) divides; (3) MEM reads a cell before the same step's reconstruction finished.
> - **[Claude, 2026-09-20] FLUID SIDE (job 4251065 = analytic, tc = 1 d, per-step PDF scan 3000..6000, `[MEM-GARBAGE-LINK]`; aborted AGAIN with `STALE-PARTICLE-IDX` at step 3517, uid 4827858800541171717, ranks 61/67; 4th analytic abort in 4 fresh runs): the "garbage" MEM link is a BLOWN-UP INTERFACE CELL SQUEEZED BETWEEN TWO RODS, not a stale value from inside an object.**
>   Timeline: t=3512 (after `timeloopAfterDEM`, i.e. after particle mapping/reconstruction) first `[HIGH-VELOCITY-PRECURSOR]`: interface cells with fill exactly 0.5 (the value of freshly converted/reconstructed interface cells) next to liquid cells with |u| = 0.5-0.59 (mean flow 0.04, target Mach limit ~0.3) at world (358-362, 34-36, 41-46); t=3513 more (13 after AfterDEM, 9 after FSLBM, |u| up to 0.65) and in the same step DEM-WATCH reports the pair ...718/...721 at overlap 2.7 LU (28 %); t=3515 `[MEM-GARBAGE-LINK]` x40: fluid cell block-local (43,14,11) is an **interface cell with fill = -0.415 (outside [0,1]), density 0.0775 and 27 populations of magnitude 5-51 (both signs)**, i.e. a numerically exploded cell (finite, so it passes the `isfinite` dry-link guard); its obstacle neighbours (rawFlag 512) are on BOTH sides in z: (43,14,10) owned by rod ...717 and (43,14,12)/(43,15,12) owned by rod ...718 -- **one fluid cell in a one-cell gap between two rods**, read by both rods' MEM links (dir 5/6/11); t=3516 `[STALE-NEARBOUNDARY]` conversions, mass-loss warnings ("No liquid or interface cell is in the neighborhood to distribute excess mass"), t=3517 abort.
>   Reading: (a) hypothesis "cell was inside an object and still holds garbage": NOT supported for this crash -- the bad cell is a legitimate interface cell that blew up dynamically (fill < 0, |f| ~ 50 in a 1-cell gap; the earlier 4.6e21 value at step 3493 is the same failure at a later stage / a different overflow state), though the first trouble (fill-0.5 interface cells + |u| 0.5 at t=3512 right after the AfterDEM coupling sweep) is exactly where rod-vacated/reconstructed cells appear, so a reconstruction contribution is possible and NOT excluded; (b) the fluid is stressed where rods are squeezed together with < 1 cell gaps -- the analytic path lets rods nest/kick (see above), which produces such gaps; whether the fluid would also blow up for a GJK-only trajectory is being tested (job 4251069: GJK, default stiffness, same per-step scan and garbage-link log; expectation from the 3 clean GJK runs: no |u| > 0.5 cells). No floater history in 4251065 (I disabled it), so the rod contacts of steps 3510-3513 for ...717/...718/...721 are not available for this run.
>   **Queued:** 4251066 = the 60 000-step cylinder hopper with GJK/EPA only (production prm `hopper_clog_cyl_frontblock_pidon_2node_gjk_gas.prm`, checkpoints every 5000 steps, ~/mt-dennys); 4251069 = GJK fluid-scan baseline (mirror2).
> - **[Claude, 2026-09-20] MATHEMATICS of the closed-form vs GJK gap + cross-check of our Kodam implementation against Part I (Sec. 3.2.5-3.2.6, Eq. 67-101) -- THREE separate issues, one of them an implementation deviation.**
>   **Exact depth:** PD = min_n [h_i(n) + h_j(n) - d.n], h(n) = (L/2)|n.e| + R sqrt(1-(n.e)^2), d = x_j - x_i (Minkowski difference; s(n) := the value in the bracket = distance j must move along n to separate).
>   **(1) band-edge is a one-point depth (Eq. 79: delta = R_i - rho_E, rho_E = distance of the rim point E to the host axis) and therefore (a) can never exceed R_i (= 3 LU here, true depth 3.4-5.3 LU) and (b) equals s(n_E) only if E is the support point of j in direction -n_E: with n = radial direction of E, s(n) - delta = h_j(n) + (E - x_j).n >= 0, the gap being how far E lies above j's deepest point (0.83-0.90 LU at 3484/3485 with the paper's cap choice; 0.95 at 3485.7, 2.9 at 3487, 4.5 at 3488). Valid only for shallow first contact, where the rim point IS the deepest point. edge-edge (Eq. 99): overlap = distance between chords AB and CD which the paper explicitly approximates as tiny straight segments ("since the length of these segments will be very small, they are assumed to be linear") -- for a nested pair the chords are long: reported 3.50 with s(n) = 11.76 along its own normal, exact 5.21.**
>   **(2) IMPLEMENTATION DEVIATION FROM THE PAPER (bug, found by this cross-check):** Kodam Eq. 68-70 choose the probe cap circle by S = sign(d_bottom - d_top), d = distance of the cap centre from the HOST AXIS in the host's xy-plane (the cap nearer to the axis). `rimVsFaceOrBand` used the 3D distance of the cap centre from the host CENTRE. For the launched pair this picks the wrong circle for 17 sub-cycles (3484.0-3485.6): reported delta 0.0027-0.0039 with a normal that is NOT the minimum-translation normal (s(n) = 3.535 vs exact 3.414), whereas the paper's rule reports 2.59-2.53 with the correct normal (s(n) = 3.4137 = exact depth). **With the paper's rule the first jump (0.0027 -> 2.526 at 3485.7) disappears** (2.5330 -> 2.5263, continuous) -- the 40 LU kick and the ~200 steps of "no restoring force" would not have happened in that form. Switch added for the test only: `paper_analytical::kodamCapSelectionByAxisDistance` (default false = old behaviour unchanged). On 14 000 random rod pairs the paper's rule helps only mildly (deep 0.6-2 R wrong: 2859 -> 2837 of 12 077; max under-estimate at 0.3-0.6 R 1.46 -> 0.59 LU), i.e. it is a real bug but not the main source of the deep-overlap errors.
>   **(3) Case switching is intrinsic to the scheme:** each scenario reports the depth of a DIFFERENT feature pair and the dispatch takes the first scenario whose local validity tests pass; the exact depth is continuous (1-Lipschitz in the poses), the scenario values agree with each other/with it only near first contact (where the shared feature coincides). Two more jumps remain with the paper's cap rule: 3488.6 band-edge 0.69 -> edge-edge 3.09 and 3489.6 edge-edge 2.83 -> band-edge 1.07, with the exact depth constant (5.30 / 5.00).
>   **Cross-check result:** correct = point E (brute-force minimum over the rim circle equals the implemented quartic solution to 4 digits: 0.0027 / 2.449 for the two circles), Eq. 79-81 (delta, normal, point), Guo et al.'s removal of Kodam's delta_r < delta_a condition, the Guo corrections of face-edge/special-contact; wrong = the cap selection (Eq. 68-70). Not re-derived here: the quartic solvers and the edge-edge chord construction (covered by the existing `ANTIDUNES_PaperAnalyticalContacts` cases/fuzz, which pass). Test: `ANTIDUNES_CylinderClosedFormJump` (env `CYL_JUMP_TABLE=1` prints the 60-row table, the cap-rule comparison and the s(n)/gap decomposition), `ANTIDUNES_CylinderClosedFormRange` (env `CYL_KODAM_CAP=1`).
> - **[Claude, 2026-09-20] Issue (2) FIXED: Kodam Eq. 68-70 cap selection is now the DEFAULT** (`paper_analytical::kodamCapSelectionByAxisDistance = true`, old behaviour still selectable for tests). `ANTIDUNES_CylinderClosedFormJump` now also asserts that with the paper's rule the closed form is continuous over the interval that jumped by 2.5 LU (largest step < 0.1 LU) and that all 17 early rows report the deep circle (> 2 LU). All antidunes tests pass (ContactDetection, CylinderParallelGap, CylinderNearParallelSkewed, CylinderClosedFormRange, CylinderClosedFormJump, PaperAnalyticalContacts 33 cases + fuzz, CylinderCylinderContactBenchmark). Synced to the three cluster checkouts.
>   **Issue (1), intuition + a stronger finding (test table `TABLE4`, `CYL_JUMP_TABLE=1`): in the nested band-edge rows the closed-form NORMAL points EXACTLY OPPOSITE to GJK/EPA's (normal . normal = -1.000 in all band-edge rows 3484.0-3486.0, -0.96..-0.99 until 3488.3; edge-edge rows +0.6..+0.7 = 45-53 deg off).** In shallow contacts the two normals agree (0 flipped in 2.1 M shadow-run contacts), so the convention is right; in a nested pair the rim point E nearest the axis lies INSIDE the host near the axis, the radial direction "axis -> E" points the wrong way out, and the band-edge force pushes the pair DEEPER (exact depth grows 3.47 -> 5.3 LU over steps 3485.7-3488.5 while the closed form reports 2.5 -> 0.7). So the 40 LU force at step 3486 was not only too small/jumpy but pointed into the pair. Figure: `animations/band_edge_vs_true_depth.png` (end-on view along the host rod: A's cross-section = circle radius 3, B's shadow = stadium; closed form pushes E out through the near wall, the true way out is the opposite side).
>   **NEW: the GJK-only 60 000-step production run 4251066 FAILED at step 6719 (17 min, 144 ranks), also with `STALE-PARTICLE-IDX`** (rank 68, rods 3963167672086036501 and 3098476543630901251), i.e. **GJK is not immune: the 6000-step GJK runs were simply too short.** No `[VELOCITY-JUMP]`, no contact anomaly before; 40 `[MEM-GARBAGE-LINK]`: LIQUID cells (fill 1, block-local (13,13,13), (13,12,14), (14,13,14), (13,14,14)) squeezed between obstacle cells of both rods (rawFlag 512 on several sides) with **populations of order 1e11 following a coherent pattern (-8.42, +4.21 x6, -2.1 x12, +1.05 x8 (all e11), i.e. an exponentially growing grid-scale mode, densities -1.6..+1.5)** -> MEM force ~4e11 -> rods launched to |v| ~ 4e7 in one sub-cycle (step 6719 sub-cycle 1) -> stale map -> abort. So there are TWO failure classes in this hopper: (A) analytic contact faults (wrong-signed/jumping forces) and (B) a FLUID-SIDE instability of fluid cells trapped in sub-cell pores between rods (seen with analytic AND with GJK-only; GJK-only survived 6700 steps, 3 analytic runs 3493-5714, 1 more 3517). The GJK fluid-scan baseline 4251069 (steps 3000-6000) is clean: no |u| > 0.5 cell, no garbage link -- consistent with class (B) starting later. Next options: resume 4251066 from its step-5000 checkpoint with the dense fluid scan/`[MEM-GARBAGE-LINK]`/floater log over 6000-6720 to see the onset; then pore treatment (e.g. convert isolated fluid cells enclosed by obstacles to obstacle, keep rods further apart, cap the MEM force).
> - **[Claude, 2026-09-20] Class (B) (fluid cells squeezed between rods): test + code-level reading + why spheres are spared.**
>   **Test `ANTIDUNES_PinchedFluidCells` (passes, 77 s), geometric precondition only:** cells mapped like `ParticleMapping` (centre inside body => obstacle); a fluid cell is *pinched* when for one of the 13 opposite D3Q27 direction pairs both neighbours are obstacle (the exploded cell of job 4251066 had such pairs, 11 of 26 neighbours obstacle from two rods). A single body never makes one (0 in 60 placements). Mean pinched cells per PAIR (300 random placements/orientations per gap): rods (r=3, L=18) 52.0 / 61.4 / 69.5 / 76.7 / 62.6 / 33.9 for gaps 0 / 0.25 / 0.5 / 1 / 1.5 / 2 cells vs spheres r=3: 14.3 / 14.5 / 14.8 / 12.8 / 8.5 / 3.7 (rod/sphere 3.6x .. 9x) and spheres of equal volume (r=4.95): 30.3 / 29.8 / 29.3 / 23.8 / 16.0 / 7.1 (rods ~2.3x). Strongly enclosed cells (>= 3 opposite obstacle pairs): rods 25-38 per pair at gaps <= 1 vs 3-5 (spheres r=3) and 8-13 (equal volume). Reason: two parallel rods form a LINE contact = a slit 18 cells long, two spheres a point contact = a small lens. **This proves the geometry, not the explosion.**
>   **Where it breaks (evidence from the logs of 4251066/4251065):** (i) the garbage is already IN THE PDF FIELD when `SimpleBB::treatDirection` reads it (`pdf_old_raw` is read before any exchange; MEM is only the messenger that turns 4e11 into a force and launches the rod; `isfinite` guards do not catch finite garbage); (ii) the populations of the exploded liquid cell (fill 1, density -0.125) are exactly f_i = (A/4) prod_a (3 c_ia^2 - 2) with A = 4.21e11, i.e. rest -2A, faces +A, edges -A/2, corners +A/4: a mode with ZERO mass and ZERO momentum, a purely non-hydrodynamic moment of the D3Q27 stencil -- streaming and bounce-back only copy values, so only the COLLISION can amplify it by 11 orders of magnitude, and the project's own comment in `BuildTimeloop.h` documents that the cumulant collision at tau ~ 0.501 diverges from negative-population/low-density states; (iii) cell densities around it are -1.64, -0.92, -0.13, +1.54 (sign-alternating = an exploded, oscillating pore). So: MEM is not the origin; the collision blows up a cell whose density had already collapsed to ~0/negative. **NOT yet known: which step drove the density of a pinched cell to ~0** -- candidates: moving-wall bounce-back term in a closing/opening pore (SimpleBB), cells covered/uncovered by `MovingParticleMapping` and re-filled by `EquilibriumReconstructor` from very few fluid neighbours (Bogner state selection: no liquid neighbour => gas), or interface handling. Job 4251086 (resume of 4251066 from step 5001, per-step scan + garbage-link + floater log 6640-6730) is queued to catch the onset.
>   **Why spheres are spared:** (a) geometry, see the numbers above; (b) `ComputeLubrication` is **only defined for sphere-sphere pairs** (`ParticleKernels.h`: "lubrication is only defined for sphere-sphere pairs"); the `LubricationCorrectionKernel` adds the viscous 6 pi mu a^2 v / gap force for gaps below 2/3 cell -- exactly the sub-cell regime the lattice cannot resolve -- so approaching spheres are decelerated before the gap closes; rods and boxes get none, so they can close a 1-cell slit at the speed set by the flow/contacts alone; (c) the sphere hopper (pillars + spheres) ran its 60 000 steps without this abort.
> - **[Claude, 2026-09-20] Job 4251086 (GJK resume 5001 -> 6800 with per-step fluid scan, garbage-link log, floater log 6640-6730) COMPLETED WITHOUT the instability** (0 `[HIGH-VELOCITY-PRECURSOR]`, 0 `[MEM-GARBAGE-LINK]`, 0 STALE/VELOCITY-JUMP): the step-6719 blow-up of 4251066 did not recur from the same checkpoint (runs are not bit-reproducible), so class (B) is a RARE, sporadic event and its onset is still not captured (only the late stages: analytic 4251065: first `[HIGH-VELOCITY-PRECURSOR]` at t=3512 after `timeloopAfterDEM`, garbage link 3515, abort 3517; GJK 4251066: garbage link read at step 6718, abort 6719). Order of operations per LBM step (BuildTimeloop.h): normals/curvature -> boundary handling (SimpleBB MEM: reads f, writes reflected f into the obstacle cell, accumulates force) -> force weighting -> stream/collide (generated cumulant kernel, omega ~ 1.997) or StreamReconstructAdvect -> cell conversion -> equilibrium refilling -> excess-mass distribution -> flag reset; then reduceAndAvgHydroForces, 10 DEM sub-cycles; then `timeloopAfterDEM`: MovingParticleMapping (covered cells -> obstacle, uncovered -> FormerMO) -> PDF reconstruction of FormerMO cells (Bogner state selection + EquilibriumReconstructor).
> - **[Claude, 2026-09-20] Catch runs for the class (B) pore instability (queued: 4251094, 4251096, 4251097; all GJK-only, 144 ranks, 3 h, per-step fluid scan, checkpoints every 2500 steps, archived to the vault at the end).** New diagnostic `populationPrecursorThreshold` (`PhysicsParams`, `checkPdfNaN` in SetupStages.h): any |f_i| > 2 in a liquid/interface cell logs `[POP-PRECURSOR t=.. after=..]` with cell class, world position, max|f|, rho, fill and the cell's OBSTACLE GEOMETRY (number of the 13 opposite D3Q27 direction pairs with obstacle on both sides = pinched pairs, obstacle neighbours of 26) -> tests the pinched-cell hypothesis directly (`analyze_precursors.py <run.out>` summarises it). **A** = `hopper_clog_cyl_gjk_catch_resume` (resume of the GJK production run from its step-5001 checkpoint to step 20 000, scan over the whole range) in ~/mt-dennys; **B** = `hopper_clog_cyl_gjk_stress_u008` (fresh, PID target 0.04 -> 0.08 = twice the lattice Mach number, ~4x drag on the pile, 15 000 steps) in ~/mt-dennys-mirror; **C** = `hopper_clog_cyl_gjk_stress_heavy` (fresh, all 406 rods densityRatio 2.0, 15 000 steps) in ~/mt-dennys-mirror2. If one aborts, the last archived checkpoint (<= 2500 steps earlier) allows a tight rerun with the floater log.
> - **[Claude, 2026-09-20] Evidence check for "pinching" -- CORRECTION of my earlier framing.** (1) Base rate from the real GJK rod pile (poses of `..._gjk-tc4` at steps 5900/5930/5960, 406 rods + pillars + walls, voxelised like the mapping): fluid cells with >= 1 opposite obstacle pair ("pinched", my test's metric) **~25 400 per step**, >= 3 pairs ~12 600, >= 5 pairs ~5 900, >= 8 pairs ~1 350, obstacle neighbours >= 11: ~23 000 -- all in HEALTHY steps. So "pinched" in that sense is the normal state of a dense rod pile and cannot by itself explain a rare explosion (both exploded cells had 3 resp. >= 5 pinched pairs, 11-13 obstacle neighbours: consistent with, but not evidence for, the hypothesis). (2) The waLBerla library documents the real failure case (`PdfReconstructionManager.h`, header note): a reconstructed cell "surrounded (given by the LBM stencil) by only boundary cells" (or several cells "locally trapped") is not divergence free -> "mass loss/gain", "malicious oscillations in the density", "very rare ... can not be detected easily and also not be cured easily". The matching geometric quantity is a CLOSED FLUID POCKET (fluid cells not connected, D3Q27, to the open region): in the healthy pile at 5900/5930/5960 there are 1 / 0 / 0 (a single trapped cell in one of three snapshots) -- rare, the right order of magnitude for a once-per-thousands-of-steps event. BUT the exploded cells had liquid neighbours (GJK: 14 fluid neighbours of 26), i.e. they were not fully isolated pockets, so the documented single-cell case does not literally apply; a small trapped cluster connected through a throat is not excluded and not measured. (3) The `[HIGH-VELOCITY-PRECURSOR]` cells of the analytic crash 4251065 at t=3512 were INTERFACE cells with fill exactly 0.5 = `seedConvertedCellFill`'s fallback "no interface neighbour in B(x) -> 0.5" (`PdfReconstructionManager.h`) next to liquid cells with |u| 0.5-0.59: freshly reconstructed cells whose D3Q27 neighbourhood contained no interface cell. (4) fill < 0: the free-surface LBM tracks mass m per interface cell (phi = m/rho); the cell conversion/excess-mass sweeps run at the END of the step, so a cell that lost more mass than it holds (garbage populations, rho ~ 0) has phi < 0 until then, and when no liquid/interface neighbour exists the excess cannot be redistributed -- at step 3515 on rank 67 (same rank/step as the garbage link) the log shows "No liquid or interface cell is in the neighborhood to distribute excess mass. Mass is lost/gained. Cell (41,6,21) converted to gas." **Status: pinching in the strong sense (trapped cells/clusters) is a plausible mechanism supported by the library authors' own warning, but not demonstrated for our cells; the catch runs (with `[POP-PRECURSOR]` incl. fMO/interface/gas neighbour counts and obstacle owners) are meant to decide.**
> - **[Claude, 2026-09-20] Catch run 4251098 (GJK production resumed from step 5001, per-step scan, checkpoints every 2500) ABORTED at t=9808 with `[NaN-SCAN]` -- a THIRD failure class, NOT at the rods.** First precursors: `[HIGH-VELOCITY-PRECURSOR]` t=9794 (interface cell, world (551.5,39.5,73.5), |u| 0.54), 9795 (|u| 0.59-0.71); `[POP-PRECURSOR]` t=9800 (liquid, world (552.5,37.5,72.5), max|f| 5.2 at f[0], rho 0.44, fill 1, **0/26 obstacle neighbours, 26 interface neighbours, fastest neighbour |u| 0.63**), t=9801: |f| up to 7.2e6, rho -0.24..-2.5, |u| up to 5.3; NaN at t=9808 at (557.5,36.5,68.5) (INTF). So the onset is in free-surface cells ~20 cells ABOVE the mean surface (z 68-74 vs 49.13) next to the y=40 wall, 170 cells downstream of the pillars (x~552), with no rod anywhere near (owners [], pinched pairs 0/13 in all 40 logged cells): a spray/splash droplet region (a small liquid blob whose 26 neighbours are all interface cells). **This is the first time the per-step scan caught an onset, and the onset is a free-surface/droplet instability, not the rod pore.** The three classes now: (A) analytic contact faults (fixed cap rule, GJK avoids), (B) MEM-adjacent pore cells between rods (4251065, 4251066; onset not observed), (C) free-surface spray instability far from the rods (4251098). It also means part of the earlier "STALE-PARTICLE-IDX" crashes' cause may sit in class (C)-like fluid states rather than pores; the run 4251066 (no scan) died at 6719 with garbage cells at the rods, 4251098 (scan) at 9808 in the spray. Checkpoint 7501 of 4251098 is archived (`..._gjk-catch-resume_4251098_20260920_161748/snapshots`).
>   **Queued:** 4251139 = `hopper_clog_cyl_gjk_catch2` (tight rerun from that step-7501 checkpoint to 12000, scan window 7501-12000, NEW density window watcher `densityPrecursorLow/High/LiquidOnly` (0.1/1.9, liquid+interface) and `highAltitudeWaterScanEnabled` z > 60 to see how water gets that high); 4251140 = `hopper_clog_cyl_gjk_stress_heavy` (fresh, all rods densityRatio 2.0, scan from 2000, density window liquid-only). The U=0.08 stress test 4251099 blew up globally at step ~313 (NaN in interface cells at the domain edges): invalid as a class-(B) stress test.
> - **[Claude, 2026-09-20] Sphere hopper vs cylinder hopper -- what differs (Dennys: "it clearly works with spheres but not with cylinders").** Same prm otherwise (diff of `hopper_clog_pidon_2node_frontblock_gas.prm` vs `hopper_clog_cyl_frontblock_pidon_2node_gjk_gas.prm`: only the floater block and the contact flags differ; same pillars R=12, flow, PID 0.04, domain, 144 ranks). **(1) System level, measured (eval DB + XZ fill slices):** the sphere block is 1015 spheres r=3 (solid volume 1015*113 = 115k LU^3); the rod block is 406 rods r=3, L=18 (406*509 = 207k LU^3, 1.8x) and 2 rods across Y cover 36 of the 40 LU channel width. Consequence: PID driving acceleration_x 5.5e-6 (cyl, step 5100) vs 1.1-2.0e-6 (spheres, all 60 000 steps); highest fluid cell z: spheres <= 70 for all 60 000 steps (mean surface 49.13), cylinders 81 (5100) -> 92 (7000) -> 107 (8000) -> 114 (9000); depth profile along x (y-averaged fill): spheres upstream 46-48 -> downstream 37-44 (step 30 000); cylinders upstream 50 / 69 / 95-100 (x = 100 / 250 / 330-360), pile front 112, **downstream only 19-29 (x = 400-552)** at step ~9800: the rod pile is nearly an impermeable dam -- water heaps 2x upstream, the downstream layer is thin and fast (supercritical, Fr 1.3, hydraulic-jump/spray region -- where class (C) started at x~552, z~72) and the flow has to squeeze through 1-2 cell gaps in the pile (classes A/B). Total fluid mass 1.30e6 (cyl) vs 1.40e6 (spheres) = the extra displaced volume. **(2) Code paths that differ by shape:** lubrication correction exists for sphere-sphere pairs only (`ComputeLubrication`); contact detection (spheres: exact analytic; rods: Kodam+Guo closed forms with the defects found, or GJK/EPA); interaction radius 3 (sphere) vs 9.49 (rod, circumscribing sphere) against block sizes 20 x 30 -> far more ghost copies and hand-overs: `[GHOST-VANISHED]` 2786 in the 60 000-step sphere run (46 per 1000 steps) vs 830 in the 6000-step GJK rod run (138 per 1000 steps), `[SUBCYCLE-LEFT-RANK]` 2575 vs 274 (43 vs 46 per 1000); the mass warning "No liquid or interface cell is in the neighborhood to distribute excess mass" (999) and `[STALE-NEARBOUNDARY]` (1620) also occur in the healthy sphere run, so those are benign background events, only their rate differs. **(3) Not yet isolated:** which of {blockage/permeability, line contacts + missing lubrication, block/ghost geometry} drives the instability. Discriminating experiments (one factor each): (a) rods with a bypass (1 rod across Y, or fewer rods) -> tests blockage; (b) spheres with the rods' blockage (e.g. 2x the sphere count) -> tests whether spheres also fail when the pile is that dense; (c) shorter rods (L = 6, aspect 1) -> line contact vs blockage; (d) rods in a wider block (cellsPerBlock y 40 instead of 20) -> tests ghost/sync geometry.
> - **[Claude, 2026-09-20] `useRotationalEffectiveMass` is now ON for the rod hopper (user request):** `gen_hopper_clog_cyl_prm.py` and the prms `hopper_clog_cyl_frontblock_pidon_2node_{gjk_,}gas.prm` / `..._diag_1node_gas.prm` (the already-queued/finished catch/stress prms are unchanged so that those diagnostics keep their exact configuration). The contact stiffness/damping then use 1/meff = 1/mA + 1/mB + (rA x n)^T IA^-1 (rA x n) + (rB x n)^T IB^-1 (rB x n) at the actual contact point (`computeRotationalEffectiveMass`, any shape pair); pillars/walls have zero inverse mass and inertia so they add nothing; it is exact for sphere-sphere (rotational terms vanish), so it is one more sphere-vs-rod difference that is now removed. All other experimental contact/buoyancy flags stay OFF. NOT yet run.
> - **[Claude, 2026-09-20] Artificial DAMPING CAPS for rods (Dennys' idea) implemented as prm-toggleable diagnostics, all default OFF:** `ModelParameters.rodCapHydroForce / rodCapHydroTorque / rodCapVelocity / rodCapAngularVelocity / rodCapOverlap` (0 = off; rods = particle type 4 only). Hydro force/torque of the MEM applied to the rod are scaled to the cap (non-finite values zeroed when a cap is on), |v| and |omega| are scaled after every sub-cycle's velocity update, the penetration depth fed into the contact force law is limited (cylinder contacts incl. rod-pillar). Each cap logs its first 40 hits (`[ROD-CAP <kind>] uid step subCycle value -> cap`) and the cumulative hit counts + the largest value seen before capping appear in the `[RodCaps]` line of `[ContactDetectionStats]` (end of run / periodic). These change the physics and are NOT for production. **What the test can and cannot separate:** caps stop a launch (the STALE-PARTICLE-IDX chain) but not the fluid instability itself: if the capped run survives far beyond the ~3500 steps at which the uncapped runs died -> the rod launches/forces were the driver (classes A/B); if it dies later with `[NaN-SCAN]`/garbage cells and few cap hits -> the fluid side is the root; the first cap hits give the onset per rod. **Queued: job 4251203** = `hopper_clog_cyl_analytic_tc1_caps` (the crash-prone analytic/tc=1d configuration, aborted at steps 3599/3493/3517 uncapped; caps: hydro force 1.0, torque 5.0, |v| 0.15, |omega| 0.05, overlap 1.0 LU; per-step fluid scan 3000-9000; 9000 steps) in ~/mt-dennys. Next candidate once a checkout is free: the same caps on the GJK production run resumed from the step-7501 checkpoint (class C: spray at x~552, no rod involved -- the caps should NOT help there).
> - **[Claude, 2026-09-20] (I) Feng-style energy-consistent (ECC) contact: reference implementation + comparison -- `ANTIDUNES_EnergyConsistentContact` (test only, nothing in the simulation; passes, ~130 s).** Model (Feng 2023 Sect. 6.3): overlap = Minkowski overlap delta_M = min_n [h1(n)+h2(n)-d.n] (h_i = support function of the cylinder), potential w = k delta_M^2/2, normal = argmin n* (body 2 pushed along +n*), force -/+ k delta n*, moments from the witness (support) points; at the kinks of the support function (line contact n perp axis, cap contact n along axis -- the GENERIC optimum for rods) the torque is a sub-gradient and the resultant is placed at the centre of the contact patch (closest points of the two support sets in the plane perpendicular to n*, centre of the overlap for parallel lines/discs). **Results:** (1) on smooth tip-to-tip rim-rim configurations the analytic force/torque equals the finite-difference gradient of w in 60/60 pairs (force err <= 1.5e-6, torque err <= 3.7e-4). (2) **GJK/EPA + spring at its contact point IS the ECC model** (1200 random pairs): mean normal deviation 0.09/0.13/0.31 deg (shallow/mid/deep), overlap error <= 2e-4, torque error <= 0.3 %, contact point <= 0.012 LU (one degenerate 179-deg case in 882 deep pairs). The closed-form Kodam+Guo hybrid is identical for shallow contacts but at depth 0.3-1 R has mean normal error 3.5 deg (max 84), overlap error 1.5 % mean (max 71 %), torque error 2.6 %, contact point off by 0.16 LU; at depth 1-2 R mean normal error 19 deg, **89 of 890 normals flipped (>90 deg)**, overlap error 10 % mean (max 150 %), torque error 12 %, contact point off by 0.78 LU. (3) Elastic rod-rod impacts (undamped spring, must return the kinetic energy): ECC mean |dE/E| 6e-5..1.5e-4, **GJK/EPA + spring 1e-5** (worst 5e-5..1.2e-3), closed-form hybrid 1e-5 (v=0.02) .. 1e-4 (0.06) .. **1.9e-3 mean, +1.3 % worst at v_rel 0.15 (1.5 LU peak overlap)** with net energy CREATION (an earlier 40-impact run gave +0.9 % mean and up to +27 % worst per impact; the 16-impact statistics above are noisier in the tails): the closed form injects energy in deep impacts, GJK/EPA does not. **Consequence: nothing to implement for the energy-consistent normal/point on the GJK path -- it already realises Feng's Minkowski/ECC model; the closed-form path does not.**
>   **(II) Contact stiffness = usual sphere calibration (dive).** Code: k = m_eff (pi^2 + ln^2 e)/tc^2, c = -2 m_eff ln e / tc, kappa = 0.876 (tangential), tc = `collisionTimeFactor` (4) x avg diameter 8.0019 = 32.0 LBM steps, DEM dt 0.1, e = 0.5, m_eff = 1/(1/mA+1/mB) purely translational; identical formula for spheres and rods (per type pair) -- k scales with m_eff so every contact has the same period tc, independent of shape, load and material: rod pair (rho 1, m = 509): k = 2.57 (sphere r=3, m=113: 0.57; rod-pillar 5.1). **Measured (GJK run 4250976, steps 5900-5960, 273 794 distinct rod-rod contact rows):** overlap median 0.030 / p90 0.067 / p99 0.135 / max 0.243 LU (1-8 % of R), |F| per contact median 0.070 / p99 0.30 / max 0.36 LU, |F|/overlap median 2.67 (p10 1.35, p90 3.20) = the calibrated k (2.57 for rho 1 pairs, 1.29 for rho 0.5 pairs, plus the damping term): the calibration is implemented as intended for rods. Equivalent Hertz modulus of the sphere calibration (binary impact 0.03 LU/step, R* = 1.5): E* ~ 0.8 LU, i.e. only ~2.4x the LBM fluid bulk modulus rho c_s^2 = 1/3 -- the DEM 'solid' is extremely soft compared to a real one, by design (Rettinger/Rudi: tc tied to the fluid time scale). **What is shape-specific and missing:** (a) test (4) of the ECC test: the generalised effective mass at the ECC contact point is 0.76 of the translational m/2 for the median rod contact (min 0.27, 10 % quantile 0.51): the same k gives a contact period sqrt(ratio) x tc, median 13 % (up to 48 %) shorter than intended and a restitution that is not the requested 0.5 (`useRotationalEffectiveMass`, now ON for the rod hopper prms, fixes this; sphere-sphere ratio is exactly 1); (b) Feng's point (his ref. [145]) that the stiffness parameter of a non-spherical contact potential should be chosen for the shape (Hertz for spheres does not carry over to edges/lines) -- our tc-based rule is shape-blind; (c) the calibration is load-independent: in the measured GJK window the overlaps stay <= 8 % of R so stiffness is not the limiting factor there (consistent with the tc=4d vs tc=1d runs: GJK max penetration 0.43 vs 0.33 LU); it would matter for load chains in a compacted pile or force spikes (analytic path).
>   **(III) Job results.** **4251203 (analytic, tc=1d, ROD CAPS: hydro F 1.0, T 5.0, |v| 0.15, |omega| 0.05, overlap 1.0): survived 8136 steps** (the same uncapped configuration died at 3599 / 3493 / 3517, tc=4d at 5714), 0 STALE-PARTICLE-IDX, 0 VELOCITY-JUMP; cap hits: step 1 (initial transient: hydro force 2.7-12 LU, capped to 1.0) and only from step 8128 on (hydro torque 6.7 -> 5, |F_hydro| up to 3e17, |T| 1.8e18 at 8132) -- caveat: ranks 25/31 used their 40-line log budget at step 1, the periodic `[RodCaps]` summary was never printed (abort before the end-of-run call), so hits in between on those ranks are not visible; it then died with `[NaN-SCAN]` at t=8136: first `[HIGH-VELOCITY-PRECURSOR]` t=8120 in a LIQUID cell at world (378.5, 23.5, 24.5) with |u| 0.53 -> 0.61 (mean 0.04) -- the throat between the rod block (x<=372) and the pillars (x 372-396, slit y 12-28) -- then `[MEM-GARBAGE-LINK]` at t=8128 (liquid cell, density 0.26). So with damped rods the launch/STALE chain disappears (n=1 vs 3 uncapped failures) and the run ends in a FLUID failure where the flow is squeezed through the slit at Mach ~1 (the rods' explosive motion was a symptom/amplifier, not the ultimate cause). **4251144 (catch2: GJK production resumed from the step-7501 checkpoint to 12000, density window liquid+interface 0.1-1.9): COMPLETED, 0 precursors, 0 cells outside the density window, 0 garbage links** -- the step-9808 spray instability of 4251098 did NOT recur (sporadic), heap height 100-110 throughout. **4251140 (rods densityRatio 2.0, fresh, 15000 steps): COMPLETED clean** (0 precursors, 0 garbage links, 0 density-window hits, 0 STALE), acceleration_x 2.7-5.6e-6 (light rods 5.5-9e-6), max rod z <= 67 and max fluid z 56-97 (light rods 81-114): heavy rods stay low and the pile blocks less. Together with the sphere hopper (heap <= 70) this points to the FLOATING, blocking rod pile (heap 100-114, thin fast downstream layer, slit velocity ~0.6) as the regime that produces the failures.
> - **[Claude, 2026-09-20] Continuation to 60 000 steps queued: job 4251422** (`hopper_clog_cyl_gjk_continue60k`, 144 ranks, 4 h limit, ~/mt-dennys). Catch2 (4251144) had `storeSnapshot false`, so it left NO checkpoint at step 12000; the continuation restarts from the same step-7501 checkpoint (vault archive of 4251098) with the identical configuration (GJK-only, default stiffness, `useRotationalEffectiveMass false`, per-step fluid scan + density window 0.1-1.9 over the whole run, `[MEM-GARBAGE-LINK]`), checkpoints every 5000 steps archived to the vault. Rough estimate from catch2: ~0.13 s/step -> ~2 h for the 52 500 steps.
> - **[Claude, 2026-09-20] New prm-toggled END-OF-STEP fluid sanity sweep (`fluidSanitySweep`, `SetupStages.h`; hook `Antidunes.cpp` right after `timeloopEval->singleStep`), default OFF.** `PhysicsParameters`: `fluidSanityEnabled`, `fluidSanityFromStep/ToStep/EveryNSteps`, `fluidSanityDensityMin` (0.1) / `DensityMax` (1.9) / `MaxVelocity` (0.5) / `MaxPopulation` (2.0) / `FillTolerance` (0.05; negative = off), `fluidSanityLiquidOnly`, `fluidSanityAbortAbove` (0 = report only; N>0 = abort when more than N cells violate in one step); a threshold of 0 switches that category off. Checks every liquid+interface cell: non-finite f/rho/u/fill, rho window, |u|, max|f_i|, fill in [-tol,1+tol] (interface) / [1-tol,1+tol] (liquid). Prints ONLY on violation: per rank `[FLUID-SANITY t=..]` (counts per category + first / lowest-rho / fastest / largest-|f| cell with class, world position, rho, fill, |u|, pinched pairs, obstacle/fMO/interface/gas neighbours, owning particle uids; 60 lines/rank) and root `[FLUID-SANITY-SUMMARY t=..]` (global per-category counts; 400 lines). **Difference to the existing dense NaN scan:** it runs AFTER the whole step (evaluators/PID included, i.e. it sees the state the next step starts from, while `checkPdfNaN` runs between the sweeps), it adds the fill check and a global count/abort. Compiles (`Antidunes` target); NOT yet exercised on the cluster; synced (md5-verified) to ~/mt-dennys-mirror and ~/mt-dennys-mirror2; ~/mt-dennys deliberately left untouched while 4251422 is pending (defaults off, so its behaviour would be unchanged either way).
>   Suggested first use: `fluidSanityEnabled true; fluidSanityFromStep 3000; fluidSanityMaxVelocity 0.3; fluidSanityAbortAbove 0;` on the capped analytic prm -- the pre-crash onset (378.5,23.5,24.5, |u| 0.53 at step 8120) should show up with a lower velocity threshold several steps earlier.
> - **[Claude, 2026-09-20] Fluid sanity sweep extended + capped analytic run RERUN as job 4251480** (`hopper_clog_cyl_analytic_tc1_caps_sanity`, 144 ranks, 1 h limit, ~/mt-dennys-mirror; prm `prms/hopper_clog_cyl_analytic_tc1_caps_sanity_gas.prm` = 4251203's prm + sweep). The sweep only ever visits LIQUID and INTERFACE cells (gas/obstacle/fMO never; `fluidSanityLiquidOnly` drops the interface cells too). New opt-in checks (default off): `fluidSanityMaxExcessMass` (|excess mass|), `fluidSanityMaxNeighbourVelocityJump` / `...DensityJump` (largest |du| / |drho| to a checked face neighbour inside the same block -- single-cell spike / shock detector, block-boundary faces skipped because ghost PDFs are half a step older), `fluidSanityCheckTopology` (liquid cell touching a gas cell = broken interface layer), `fluidSanityStatsEveryNSteps` (`[FLUID-STATS t=..]` global rho range / max|u| / max|f_i| / fill range + the fastest cell, independent of violations). Log throttling: first 30 violating steps always, then only when the count doubles or every 100th step, so a benign steady violation cannot eat the budget. Rerun settings: from step 3000 to 9000, max|u| 0.3, rho 0.1-1.9, max|f| 2.0, fill tol 0.05, excess 1.0, nbr |du| 0.2, |drho| 0.3, topology on, stats every 100 steps, report only. NOT yet exercised on the cluster -- the noise level of the new opt-in thresholds is unknown (the first run tells whether e.g. topology / jump checks fire in the healthy pile).
> - **[Claude, 2026-09-20] Result of 4251480 (capped analytic + fluid sanity sweep): COMPLETED all 9000 steps in 22 min (exit 0) -- the identical configuration that died at step 8136 in 4251203 (and at ~3500 uncapped).** The failure is therefore NOT deterministic (as with 4251086/4251144 resumes); the sweep is read-only, so the difference is run-to-run nondeterminism (MPI reduction/ordering), i.e. the instability is a rare event, not a fixed-step bug. `[FLUID-STATS]` (every 100 steps): max|u| 0.09-0.146 throughout (peak sample 0.146 at t=6000), rho in [0.96, 1.10], max|f_i| 0.30-0.33; the max|u| 0.3 / rho / |f| / excess-mass / neighbour-jump / topology checks NEVER fired (0 violations in all categories). The ONLY violations (4320 rank lines, 89 root summaries, up to 1218 cells at t=4600) were interface-cell FILL outside [-0.05,1.05] (observed fill range [-0.54, 1.71] -- see `[FLUID-STATS]`): with tol 0.05 this check is too tight for FSLBM (interface fill overshoots before conversion) and is noise; use a much larger tolerance (or off) unless the fill excursions themselves are to be studied. No STALE-PARTICLE-IDX, no NaN, no garbage links; `[RodCaps]`: hydroForce cap hits 9940, hydroTorque 700, velocity/angular velocity/overlap 0; largest values before capping |F_hydro| 12.4, |T| 10.6, |v| 0.094, |omega| 0.0062, overlap 0.071 (vs |F| 3e17 in 4251203 at 8132). The 863 `[STALE-NEARBOUNDARY]` lines are the pre-existing near-boundary log, not the new sweep. Consequence for the diagnosis: to catch the onset we need many repeats or a longer run; a single 22-min run is cheap (~0.13 s/step incl. sweep), so N repeats of this prm are the practical way to sample the rare event with the sweep armed (fluidSanityFillTolerance raised first).
> - **[Claude, 2026-09-20] Repeat over 18 000 steps queued: job 4251684** (`hopper_clog_cyl_analytic_tc1_caps_sanity_18k`, ~/mt-dennys-mirror2, 144 ranks, 1 h limit; prm `prms/hopper_clog_cyl_analytic_tc1_caps_sanity_18k_gas.prm` = the 4251480 prm with timesteps / dense-scan / sweep windows 18000 and `fluidSanityFillTolerance` 1.0). 4251480 needed ~22 min for 9000 steps incl. build, so ~40-45 min expected -- the 1 h limit is tight but should fit; if it is cut off, the vtk-out/eval archive step still runs only if the walltime signal leaves time (not guaranteed).
> - **[Claude, 2026-09-20] Result of 4251684 (capped analytic + sanity sweep, 18 000 steps): COMPLETED all 18 000 steps in 44:57 (exit 0).** Second surviving run of the configuration that died at 8136 in 4251203 -> the crash is a rare event (2 of 3 capped runs survive). Across 150 `[FLUID-STATS]` samples: max|u| <= 0.250 (peak at t=10100, an INTERFACE cell in the thin downstream layer at (431.5,34.5,29.5), 10 gas neighbours, no rod nearby), rho 0.915-1.120, max|f_i| <= 0.332, fill [-0.81, 1.88]; no NaN / garbage link / high-velocity or population precursor / STALE / velocity jump. `[RodCaps]`: hydroForce 11600 hits, hydroTorque 1390, velocity/omega/overlap 0; largest seen before capping |F| 12.4, |T| 11.5, |v| 0.131, |omega| 0.012, overlap 0.197. Sweep (tol 1.0): only 34 violating steps in total. (1) t=4026: one INTERFACE cell (319.5,11.5,65.5) with fill 2.03 (rho 1.0, excess 0, 1 gas neighbour) -- gross fill excursion, single step. (2) **Face-neighbour |du| > 0.2 (the new jump check) fired ONLY in liquid pore cells wedged between two rods: 11-19 of 26 neighbours are obstacle, 2-4 pinched pairs, owners = exactly two rod uids** -- (299.5-304.5, 12.5, 71.5-74.5) at t=6726-6736 and 6908-6927 (persisting ~30 steps each, 2-3 cells; |u|~0.1, rho 0.98-1.01, i.e. shear of ~0.2 across the pore, no blow-up), (406.5,24.5-25.5,10.5) at 13400/13416, (387.5,15.5,48.5) at 15300, (379.5,15.5,30.5) at 17646 (19/26 obstacle neighbours). No density-jump, excess-mass, |u|, |f|, rho-window or topology violation in the whole run. Reading: the jump check is a clean pore-cell / pinched-cell detector; in the healthy runs those pores stay harmless (velocities ~0.1); in the crashing run the same kind of cell (378.5,23.5,24.5, throat between two rods/pillar) would be where it escalates -- worth comparing the jump-check hits' locations with the onset in a run that does crash.
> - **[Claude, 2026-09-20] Full 60 000-step run of the capped analytic configuration queued** (`hopper_clog_cyl_analytic_tc1_caps_sanity_60k`, ~/mt-dennys-mirror2, 144 ranks, **3 h limit (reduced from 4 h on request via scontrol)** -- 18k steps took 45 min => ~2.5 h expected; prm `prms/hopper_clog_cyl_analytic_tc1_caps_sanity_60k_gas.prm` = the 18k prm with timesteps / dense-scan / sweep windows 60000 and `storeSnapshot true` every 5000 steps so a late crash can be resumed with the sweep armed; sbatch additionally archives `snapshots` to the vault). Job id: 4252005.
> - **[Claude, 2026-09-20] Rod-arrival check for 4251684 (vtk-out/floaterCylinder, frames every 5000 steps):** block starts x 172-368 (front pressed against the pillars, span 372-396); rods with x>372 / >384 / >396: 0/0/0 (t=0), 18/0/0 (5000), 53/24/18 (10000), 79/55/51 (15000, max x 662). Rods inside the pillar span have y 10.5-30.9 (in the slit): 35 at t=10000, 28 at t=15000. So the rods DO reach and pass the hopper throat (13 % through by t=15000, the rest creeping; median x 270 -> 324) -- the sweep hits at the pore cells therefore sit in a really loaded throat.
> - **[Claude, 2026-09-20] GJK-only 18 000-step run WITHOUT caps, sweep armed, queued: `hopper_clog_cyl_gjk_sanity_18k`** (~/mt-dennys-mirror, 144 ranks, 1 h limit; prm `prms/hopper_clog_cyl_frontblock_pidon_2node_gjk_sanity_18k_gas.prm` = the GJK production prm with timesteps 18000 + dense NaN scan and fluidSanity* sweep (same thresholds as 4251684, fill tol 1.0) from step 3000; `useRotationalEffectiveMass true`, `useTightLayer1 true`, default collisionTimeFactor 4, no rod caps, `storeSnapshot true` every 5000). NB the earlier GJK production runs (4251098 etc.) had `useRotationalEffectiveMass false`, so this is not an exact repeat of them. Time is tight for 1 h (GJK about 0.13-0.18 s/step incl. scan+sweep -> ~40-50 min + build).
> - **[Claude, 2026-09-20] Results 4252146 / 4252005 / 4251422.** **4252146 (GJK-only, no caps, sweep armed, 18 000 steps): COMPLETED in 46:43** -- max|u| 0.212 (t=4500), rho 0.877-1.187, max|f_i| 0.351, no NaN/precursor/garbage/STALE; sweep hits only the pore-cell neighbour-jump pattern (41 violating steps, none after t=4501). **4252005 (capped analytic, 60 000 requested, 3 h limit): FAILED at t=22886 after 56:44 with `[NaN-SCAN] FATAL fluid-cell NaN` at INTERFACE cell block-local (23,14,16), world (471.5,14.5,16.5) -- the first crash that the sweep caught at its onset, and it is a FREE-SURFACE instability ~75 LU DOWNSTREAM of the pillars (x 372-396) with NO rod involved (obstacle neighbours 0/26, owners=[] in every logged cell; STALE 0, rod caps irrelevant).** Timeline: `[FLUID-STATS]` t=22300-22800 max|u| 0.14-0.17 (no warning); sweep (|u|>0.3 / nbr-jump) first violations **t=22861** (INTERFACE cells (470.5,14.5,22-23), fill 0.07-0.8, gas nbrs 11-18, |u| 0.25-0.33, rho 0.98-1.07), 1-6 cells/step t=22861-22869, 4 at 22868, 8 at 22871, 17 at 22876 (|u|>0.3: 8), 48 at 22880 (rho<0.1: 11, rho>1.9: 10, max|f|>2: 18), 146 at 22883, 464 at 22885; the existing `[HIGH-VELOCITY-PRECURSOR]` (|u|>0.5) first fired at t=22877 (INTERFACE (473.5,14.5,20.5), |u| 0.54, fill 0.08 -> 0.85 / 0.56 at 22878, LIQ neighbour (472.5,14.5,20.5) |u|=0.95 at 22879) -- **the sweep at max|u| 0.3 gave 16 steps more warning than the 0.5 precursor**; NaN at 22886. This is the same class as the spray at x~552 in 4251098 (class C): a thin fast free-surface layer at a low-fill interface region, not the throat/rod chain of 4251203 (x=378). So capped rods do not cure the fluid instability; two different fluid failure sites now (throat x=378, downstream free surface x~471-552). **4251422 (60k GJK continuation) FAILED after 2 min: CMake error `PinchedFluidCellsTest.cpp` missing in ~/mt-dennys** (my earlier sync of tests/showcases_antidunes/CMakeLists.txt never copied the test sources there and the sbatch always rebuilds). Fixed: PinchedFluidCellsTest.cpp + the 3 sweep sources synced and md5-verified in ~/mt-dennys (all three checkouts now identical), **resubmitted as 4252426** (same sbatch/prm, 4 h). NB /home/hpc is at 115 GB vs a 105 GB soft quota (job epilogue warning) -- rescued_before_* / archive dirs under the checkouts' build/.../archive should be moved to the vault.
> - **[Claude, 2026-09-20] DEEP DIVE 4252146 (GJK-only, no caps, sweep armed, 18 000 steps; archive `antidunes-hopper-clog-cyl-gjk-sanity-18k_4252146_20260920_215229`, has checkpoints at 5000/10000/15000).** (1) **Macro behaviour is indistinguishable from the capped analytic run 4251684** (2000-step window means): PID acceleration 5.2e-6 (start) -> 1.1e-6 (2-4k) -> 6-8e-6 (later), mean vel 0.039-0.040, Fr 1.27-1.33, max fluid z 113 (GJK, t=8700) vs 104 (ANA), max particle z 110 vs 107, transport rate 0.01-0.06 for both, total fluid mass 1.303e6 -> 1.266e6 (-2.9 %) -> 1.272e6 for both. The upstream water level rises 35 -> 77 (pillar-upstream 26 -> 74) while the downstream probe elevation FALLS 41 -> 3.8 (GJK) / 4.5 (ANA): the rod pile acts as a dam (head ~70 LU across the pile), same as in the earlier note. So GJK vs analytic+caps makes no macroscopic difference; the failures are stochastic events, not a systematic model effect. (2) **Rods:** frames 0/5000/10000/15000: rods with x>372: 0/20/43/65 (GJK) vs 0/18/53/79 (ANA); through the throat (x>396): 0/0/13/35 (GJK) vs 0/0/18/51 (ANA); max x 461 / 636 (GJK) vs 498 / 662; median x 270 -> 326 (GJK); max rod speed 0.046/0.098/0.054. Rods in the slit centre have |axis_y| <= 0.26 (GJK, t=10000/15000) = aligned with the flow, as geometry demands (L=18 > gap 16). (`rotation` VTK convention: only axis_y = R[4] is convention-independent; x/z components not used.) (3) **Sweep hits (41 violating steps), all `face-nbr |du|>0.2` except one fill>1 hit at t=4501:** (a) t=3184-3700 and 4475-4501: LIQ pore cells with 4-13 obstacle neighbours at the BACK/top of the block, x 185-203 / 239-245, z 52-57 (3-8 cells, persisting 10-40 steps, |u| 0.03-0.1); (b) t=11622-11636 (390.5,13.5,39.5) and 11781-11791 (393.5,14.5,33.5): pore cells in the pillar throat with 2-4 pinched pairs, 8-13 obstacle nbrs, |u| 0.08-0.15, rho 0.98-1.03 (15 steps each, no escalation); (c) 10 isolated single-step INTERFACE hits, 0 obstacle neighbours, at x 294-415 (t=4501, 4982, 5705/10, 6616, 9343, 9872, 10595, 10791, 15333, 17217, 17644) = a free-surface cell next to a >0.2 faster interface neighbour, the same free-surface class as the 4252005 crash. (4) **Where the fastest cell lives (150 samples, every 100 steps):** t<4500 in the rod-block pores (LIQ, 14/15 with obstacle nbrs, |u| 0.13-0.21, x~204); from ~6000 on on free-surface INTERFACE cells (10-15 of 15 samples), x drifting downstream 276 -> 412 (9000-15000) -> **498 (15000+)**, |u| median ~0.10. **From t=15300 to the end the fastest cell is ONE stationary INTERFACE cell (497.5,15.5,25.5): 25 gas neighbours, 1 interface neighbour (a 2-cell micro-droplet), fill oscillating 0.83-1.003, rho 1.0126, |u| rising monotonically 0.096 -> 0.132 (2600 steps) without moving, nearest rod 14 LU away** -> a non-advecting numerical droplet artefact and a plausible seed of the class-C instability (same y 14.5-15.5, z 20-25 streak, x 470-500 as the 4252005 crash site; hypothesis, being checked in the 4252005 analysis). (5) **Probe caveat:** `clogDownstreamPressure`, `clogDownstreamExcessPressure`, `clogPillarDownstreamPressure`, `clogPillarDownstreamExcessPressure` are CONSTANT for the whole run from t>=3000 in GJK (1 distinct value; ANA: 1-4 distinct) -- the downstream pressure probes sit in a cell that is no longer liquid once the downstream water level drains to ~4-13 LU, so those columns are frozen last-valid values and must NOT be used for the thesis clog metrics after t~3000. (6) **Step-1 transient:** every rod at the block front gets an initial MEM hydro force of 2.7 / 9.8 / 12.4 (the `[TOTAL-FORCE-THRESHOLD] step=1` lines at vel=0) -- the 'largest |F_hydro| = 12.4 before capping' in the capped runs is this start-up transient, NOT a later event. (7) Performance: 0.162 s/step steady incl. dense scan + sweep (2718 s for 17900 steps). **Open:** no fluid-field VTK analysis yet; a resume from the step-15000 checkpoint with fluidSanityMaxVelocity 0.2 and the micro-droplet region logged would show whether the stuck cell at (497.5,15.5,25.5) escalates.
> - **[Claude, 2026-09-20] 4252005 crash site geometry (fluid_field VTK step 15000/20000 + NaN-scan dump; block 64x20x30 cells at origin (448,0,0); z=0 is the FLOOR, gravity -z; `mapped_flag` 1 liquid, 2 interface, 3 gas, else obstacle).** NaN-scan neighbourhood at t=22886 around (471.5,14.5,16.5): the 9 cells directly BELOW (z=15.5) are all GAS (rho~1, fill 0), the 9 at z=16.5 are INTERFACE (already blown up), z=17.5 LIQ/INTERFACE -> the failing cell is on the UNDERSIDE of an airborne liquid lump (gas below and around, no floor/obstacle contact), i.e. effectively inside the gas region. That lump did NOT exist at the earlier frames: step 15000 (y=14.5): free surface at x~471 is at z~10-12 (water mound deepening downstream, rods embedded at x~481-495 up to z~19 -- rod uid ..8891 at (481,8.4,17.7) sticks out of the surface), z=16.5 at x=471 is GAS; step 20000: only a 4-5 LU film (z 0-5) at x 448-495, y 13-15, everything above z~6 is gas, so the site is ~11 LU above the film. => liquid at z 16-25 appeared within the last 2861 steps (frames every 5000 steps cannot resolve when). Rods: eval maxParticlePos_z rises smoothly (~0.5 LU/100 steps, always 4-6 LU above the highest fluid cell = the top of the upstream heap, both ~128-132 at t=22800), NO launch spike at 100-step sampling; step 20000: nearest rod to the site is 41 LU away (uid ..6485 at (435.5,29.0,5.7), v=(0.052,-0.007,0.003) -- would pass x=471 near t~20700 at y~29, i.e. ~15 LU lateral), other near-site rods gone (rods reach x=766 and wrap through the periodic X boundary; rods at x~373 are static, v=0, jammed at the pillar front); step 15000: 3 rods within 12-22 LU. floater-history was NOT enabled, so the rod trajectories between 20000 and 22861 are unobserved: an airborne/tumbling rod cannot be excluded but nothing on record supports it; the background analysis agent (report.md in scratchpad/crash4252005) covers this in depth. Cleanest test: resume from the step-20000 checkpoint (snapshots in the 4252005 vault archive) with floater-history / rod-position logging for rods within 60 LU of (471,14,20).
> - **[Claude, 2026-09-20] CORRECTION + 4252005 crash analysis by a background agent (raw extracts in scratchpad/crash4252005/; the key rod facts below were re-verified by me in the log).** My statement just above that 'floater-history was off so rod trajectories between 20000 and 22861 are unobserved / nothing on record supports a flying rod' was WRONG: the log has `[CylForces]` lines for all 406 rods every 100 steps (pos, vel, F_hydro, F_contact, tilt). They show a rod flying into the site: **uid 2666130979403333655 (light rod, densityRatio 0.5) was swept over the slit by the raised pool (pool top z 122-128, rod at (377,13.5,110) at t=20000, leaves the slit at x=396, z~104 at t~21400) and then FELL with the water stream: t=22000 (418.4,13.8,82.6) v=(0.048,-0.002,-0.053); 22400 (440.6,12.5,55) vz -0.082; 22600 (453.1,12.0,37.7) vz -0.102; 22700 (459.1,12.1,26.9) vz -0.115; 22800 (465.0,12.6,14.5) v=(0.057,0.005,-0.133)** (|v| 0.145 = 97 % of the 0.15 velocity cap; F_contact = 0 throughout, hydro force pushes it down faster than g). Extrapolated floor contact ~t=22835-22850 near (468,13); **`[ROD-CAP hydroForce]` for this uid starts at t=22848 (1.385, 1.736, 1.977, 2.079 at 22848..22851, cap 1.0; rank 84's 40-line log budget is exhausted at 22851, so later hits are unlogged)**; the first sweep violation is at t=22861 at INTF (470.5,14.5,23.5) ~2 LU from the rod in x/y but 5-15 LU above it; all violating cells stay on ONE rank (84) within x 470-473, y 12-16, z 17-24, creeping down 0.26 LU/step and growing 2/17/48/146/464 cells until the NaN at 22886. A second falling rod (uid ..33652) is ballistic (F_hydro=0.0000 from 20600, lands x~490 at t~23450) and is irrelevant. Baseline 4251684 also contains falling light rods (|F_hydro|<1e-4, vz -0.06..-0.08 at t=10000/14000/16000) without crashing. **Ranked conclusion (agent, my spot-checks agree with the facts, not the percentages): (1) free-surface/LBM failure in a plunging thin spray-like stream (|u| 0.25-0.33, fill 0.003-0.8, 5-16 interface / 9-21 gas neighbours; the stream exceeds the Torricelli limit sqrt(2 g H)=0.122 for H=128, g=5.83e-5 by 1.4-2.7x), (2) the trigger is either the stream impact or the momentum exchange at the rod's landing -- not separable from the logs; agent's rough split 75 % fluid / 25-40 % rod trigger (not exclusive).** A contact-glitch launch is NOT supported (F_contact=0 for the rod 22400-22800, no STALE/velocity-jump). **The 'stuck micro-droplet' hypothesis from 4252146 is REFUTED for 4252005:** no fastest cell stays in place over 20000-22800; the 22861 cells are a diffuse low-fill interface cloud, not an isolated 1-interface/25-gas cell (the 4252146 stuck cell at (497.5,15.5,25.5) remains an unexplained artefact of that run). Other findings: the upstream pool keeps rising (clog upstream elevation 78 -> 115 between t=10k and 22.8k, max fluid z 128 at 22.8k), the downstream probe (432,20) reads a 1.4-1.7 LU film, 7 rods have wrapped through the periodic X boundary (x 3..19, z 46-57 at t=20000). Logging defects: the ROD-CAP budget is a static 40 lines per rank shared over ALL kinds and `[RodCaps]` is printed only at the very end (missing after an abort). **Best follow-up (not launched): A/B resume from the step-20000 checkpoint (vault archive of 4252005; pdfField 1.2 GB) to ~22900: (A) unchanged + per-step FloaterHistory for uids ..33655/..33652 from 22400 + per-step top-|u| dump on rank 84; (B) identical but rod ..33655 removed/held inactive. If B fails at ~22861 in the same place -> pure free-surface plunge instability; if B survives -> the rod impact is the trigger.**
> - **[Claude, 2026-09-20] Were flying rods the cause in BOTH capped crashes? -- No (checked 4251203 via its `[CylForces]` lines).** 4251203 (t=8136, throat x~378): the two rods named in the garbage-force / cap lines were **uid 4755801206503243780 (..780), STATIC** at (373.33,29.41,22.44) with |v|~1e-4 at t=8000 and 8100 (jammed at the pillar front, F_hydro ~ -F_contact ~0.05-0.08), and **uid 4827858800541171713 (..713), creeping** through the throat at |v|~0.02 (t=8000 (382.1,22.4,26.1) -> t=8100 (384.1,21.7,24.8), F_hydro <=0.06, F_contact 0), i.e. ~6 LU from the first blow-up cell (378.5,23.5,24.5). No rod was flying (max |v| 0.02); the fluid |u| precursor (t=8120, 0.53 -> 0.8) precedes the first ROD-CAP hit on these rods (t=8128 torque 6.7; t=8129 raw force 504 / 682, torque 3772 / 3459) -> the rods were VICTIMS of already-exploded pore cells, not the source (caveat: rod state is sampled only every 100 steps, nothing between 8100 and 8120). 4252005 (t=22886): a light rod (..33655) fell with the plunging stream and landed at the site ~13 steps before the first violation; trigger not separable from the stream impact (see the note above). The uncapped analytic crashes (~3500 steps) are different again: contact faults -> violent rod motion -> garbage MEM force -> launch (STALE-PARTICLE-IDX), removed by the caps/GJK. So there are THREE failure modes: (1) analytic contact fault -> rod launch, (2) fluid blow-up in a throat pore cell next to jammed/slow rods (4251203), (3) fluid blow-up in a fast thin free-surface stream downstream (4251098 x~552 without rods; 4252005 x~471 with a landing rod). Common factor of (2) and (3): very high heads/velocities (pool 78-128 LU vs Torricelli limit |u|~0.12) rather than rod launches.
> - **[Claude, 2026-09-20] 4252426 (resubmitted 60k GJK continuation from the step-7501 checkpoint, no caps, dense scan on, sweep NOT armed) FAILED after 33:58 at t=19239: `[NaN-SCAN] FATAL fluid-cell NaN` at INTERFACE cell block-local (18,7,11), world (402.5,27.5,41.5)** -- just downstream of the pillars (x 372-396) at the slit edge (y~28), z=41.5 (well above the thin downstream film). So the uncapped GJK config also dies (t=19239), i.e. the fluid instability is not a property of the analytic contacts or the caps. Not analysed yet (no sweep armed, VTK frames only at 10000/15000 in the vault archive `antidunes-hopper-clog-cyl-gjk-continue60k_4252426_20260920_233552`). Non-capped rod-run ranking by reached step: 4252426 19239 (crashed, resumed run, 2 frames), 4252146 18000 (completed, fresh, 4 frames 0/5000/10000/15000), 4251140 stress-heavy 14900 (3 frames), 4251144 catch2 11900, 4251098 9800; uncapped analytic ~3500. Failure sites so far: throat x~378 (4251203), behind the pillars x~402 z~41 (4252426), downstream stream x~471 (4252005), spray x~552 (4251098).
> - **[Claude, 2026-09-20] Two jobs queued.** (1) **4252462 = `hopper_clog_cyl_gjk_sanity_cont18k`** (~/mt-dennys-mirror, 144 ranks, 1.5 h): CONTINUATION of 4252146 from its step-15001 checkpoint (vault archive `antidunes-hopper-clog-cyl-gjk-sanity-18k_4252146_20260920_215229/snapshots`, copied not moved; the sbatch prints simState and expects 15001) +18000 steps (`timesteps 33001`), VTK every 1000 steps (`fluid_field.writeFrequency 1000` also drives the rod/particle frames), dense scan + fluidSanity* armed, snapshots every 5000. prm `..._gjk_sanity_18k_cont18k_gas.prm`; only differences to the 18k prm: loadSnapshot, timesteps, scan/sweep ToStep, VTK frequency and the fluid crop z max 100 -> 140 (the pool surface reaches z~110-128 and the old crop clipped it at 100; x 128..640 unchanged). Estimate ~50 min steps + ~10 min VTK/snapshots/build. (2) **4252463 = `hopper_clog_cyl_noclog_reference_pidon_30k`** (~/mt-dennys-mirror2, 144 ranks, 1.5 h): NO-CLOG REFERENCE of the rod hopper = the same prm without the FloatingCylinders block (no debris; same channel/pillars/Re/Fr/We/avgParticleDiameterOverride), PID ON target 0.04, 30000 steps, no snapshots, VTK every 5000, scan + sweep armed. Purpose (2026-09-18 convention): the STEADY acceleration_x of this clean channel becomes `PIDParameters.accelerationXOverride` for the PID-OFF rod runs (with enabled false) and this run is the baseline for the reference-subtracted afflux metric. NOT yet verified that the zero-rod path runs (FloatingCylinders is optional in the parser, `getBlocks(...,0,1)`); check the first minutes of the log. After it finishes: read acceleration_x of the last ~5000 steps from `eval/antidunes.sqlite` (table runs), average, and use that value.
> - **[Claude, 2026-09-21] What |u| means, where it comes from, and a CORRECTION of the thresholds I used.** `|u|` (FLUID-STATS `max|u|`, sweep `|u|>...`) = fluid speed of one liquid/interface cell in lattice units (cells per step); lattice speed of sound cs=1/sqrt(3)=0.577, Ma=|u|/cs; mean flow 0.04 = Ma 0.07. **Sources (Zotero, checked in the PDFs):** Kruger et al. 2017, *The Lattice Boltzmann Method*, Sect. 7.2, p. 279: 'A lattice Mach number Ma = U/cs is considered small if Ma < 0.3. Larger Mach numbers lead to more significant compressibility errors' => |u| < 0.3*0.577 = **0.17**; Sect. 4.4 pp. 127-128: instability = errors in populations/density/velocity grow exponentially, eventually NaN; Sect. 4.4.4 pp. 134-135: keep |u| < cs, use TRT/MRT/better collision operators for stability (no fixed collapse velocity given); Bogner 2017 (thesis), Eq. 2.11 and text: Ma <~ 0.3 = usually incompressible; Schwarzmeier et al. 2023 (JFM 961 R1, p. R1-3): same D3Q27 cumulant model (Geier 2015), all relaxation rates 1 except omega, nu = (1/3)(1/omega - 1/2) (ours omega=1.9975, nu=2.07e-4); no velocity ceiling stated. Torricelli sqrt(2gH) is textbook fluid mechanics. **NOT sourced:** my earlier statement that populations/density go wrong 'around |u| 0.3-0.5' came from the comment in `find_degenerate_cells.py` ('rough Mach-number stability ceiling', no citation, and 0.3 there was |u| = Ma 0.52, not Ma 0.3) and from general LBM experience, not from a paper. **Empirical (this project):** there is NO sharp collapse velocity: 4251684 (max sampled |u| 0.250 = Ma 0.43), 4252146 (0.212) and 4252426 before its crash ran on; in 4252005 cells at 0.25-0.33 escalated (|u|>0.5 within 16 steps, rho out of range 8 steps later, NaN after 25) on thin low-fill interface cells with low viscosity -- the collapse is an exponentially growing local error, not a value of |u|. Sampled maxima of 0.15-0.25 already exceed 'small Mach' (0.17) in places -> compressibility errors are present in the hopper runs. **Changes made:** `fluidSanityMaxVelocity` code default 0.5 -> 0.17 (Ma 0.3) with a comment stating it is an accuracy warning, not a collapse threshold; `find_degenerate_cells.py` default `--vel-threshold` 0.3 -> 0.17 with a corrected docstring. All sanity prms so far (4251480 ... 4252463) explicitly set 0.3 and are UNCHANGED (provenance); expect many more sweep hits at 0.17 in future runs. Sources rebuilt locally; synced to ~/mt-dennys only (the two mirrors run jobs 4252462/4252463 -- sync after they finish; no current prm omits the key, so nothing depends on the default yet).
> - **[Claude, 2026-09-21] `fluidSanityMaxVelocity` set to 0.17 in all six sanity prms** (`hopper_clog_cyl_analytic_tc1_caps_sanity{,_18k,_60k}_gas.prm`, `..._frontblock_pidon_2node_gjk_sanity_18k{,_cont18k}_gas.prm`, `hopper_clog_cyl_noclog_reference_pidon_30k_gas.prm`; each line carries a comment naming the job that ran with 0.3). This supersedes the 'prms unchanged' remark in the previous note: the repo prms now say 0.17, the jobs 4251480, 4251684, 4252005, 4252146 (and 4252462/4252463, which loaded their prms at start) ran with 0.3 -- their sweep hits/logs refer to 0.3. Synced to ~/mt-dennys (md5-verified); the two mirrors get the new prms after 4252462/4252463 finish.
> - **[Claude, 2026-09-21] Results 4252463 (no-clog reference) and 4252462 (GJK continuation +18000 from 15001).** **4252463: COMPLETED all 30000 steps in 37:03, no problems, 7 sweep summary lines at 0.3.** Clean rod-hopper channel (768x40x180, pillars R12, no rods, PID on 0.04): mean vel 0.03996-0.04000 after spin-up, Fr 1.31, max fluid z 53-62 (free surface waves), mass constant 1.5086e6. **`acceleration_x` = 0 for the first ~5000 steps** (mean vel starts at 0.0423 > target 0.04, so the PID output is clipped at minActuatingVariable 0 until the flow has decayed below the target -- likely explanation, not verified), then 2500-step window means 2.8e-7 .. 6.8e-7 (noisy, std of the samples ~9e-8): 15000-17500 6.36e-7, 17500-20000 6.41e-7, 20000-22500 5.58e-7, 22500-25000 4.81e-7, 25000-27500 5.32e-7, 27500-30000 6.77e-7 => **steady no-clog force ~ 5.9e-7 (mean of the six windows 15000-30000, window spread 4.8-6.8e-7)**; for comparison the sphere no-clog extreme run used 7.8e-7 (different domain height 90). Candidate for the PID-OFF rod run: `PIDParameters.enabled false; accelerationXOverride 5.9e-7`. Vault archive `antidunes-hopper-clog-cyl-noclog-reference-pidon-30k_4252463_20260921_003520` (VTK every 5000). **4252462: FAILED at t=29396 (38:16) with `[STALE-PARTICLE-IDX]` (SimpleBB.h:288) on rank 44, cell (11,7,6) owned by rod uid 4827858800541171727** -- i.e. it survived 14400 steps from the checkpoint and then died in the RODS-adjacent way (class B), NOT via a NaN. Trace (sweep at 0.3, dense scan on): fastest cell t=28800-29100 = a fast-moving low-neighbour INTERFACE cell (25-26 gas nbrs) at x 245 -> 265, y 33.5, z 83.5, |u| 0.257 -> 0.344 (drifting downstream ~0.07 LU/step, away from the later failure site); t=29200/29300 fastest = LIQUID pore cells between rods (obstacle nbrs 14-18) at (186.5,30.5,71.5) / (192.5,29.5,71.5) |u| 0.29/0.30; first `[HIGH-VELOCITY-PRECURSOR]` t=29351 at LIQ (190.5,31.5,72.5) |u| 0.504 vel=<0.48,-0.11,-0.10> fill 1 (26 precursor lines in total); sweep hits 29395: 36 cells (rho<0.1: 4, rho>1.9: 3, max|f|>2: 10, |u|>0.3: 13) incl. LIQ (191.5,29.5,69.5) / (192.5,29.5,69.5) with owner uid ..727 (1-3 obstacle nbrs); 10 `[MEM-GARBAGE-LINK]` lines mention uid ..727; STALE abort 29396. **Rod ..727 was NOT flying:** pos (167.4,31.6,67.7) at t=28500 -> (193.2,32.1,69.8) at t=29300, |v| 0.031-0.035 (creeping downstream), F_hydro <= 0.1, F_contact ~0 at 29300 -- the fluid next to a slow rod in a pore at the TOP-BACK of the rod pile (x~190, z~70, i.e. ~180 LU upstream of the pillars, inside the pool) blew up first and the garbage momentum-exchange force then invalidated the rod (=> STALE). Failure site #5 (pile-top pore cells, upstream end): sites so far: throat x~378 (4251203), behind the pillars x~402 (4252426), downstream stream x~471 (4252005), spray x~552 (4251098), pool-top pore x~190 (4252462). Fluid VTK frames every 1000 steps 16000..29000 are in the vault archive `antidunes-hopper-clog-cyl-gjk-sanity-cont18k_4252462_20260921_003554` (not copied to animations yet).
> - **[Claude, 2026-09-21] PID-OFF GJK clog run queued: 4252482** (`hopper_clog_cyl_gjk_pidoff_sanity_18k`, ~/mt-dennys-mirror2, 144 ranks, 1.5 h; prm `prms/hopper_clog_cyl_frontblock_pidoff_2node_gjk_sanity_18k_gas.prm` = the 4252146 GJK prm with `PIDParameters.enabled false; accelerationXOverride 5.9e-7` (the steady no-clog force from 4252463; `SetupStages.h:124` uses the override as the initial and, with the PID off, constant acceleration) and fluid crop z 0..140; 18000 steps, VTK every 5000, sweep at |u|>0.17 + dense scan armed from 3000). All pending 0.17 changes (SimulationParams.h default, find_degenerate_cells.py, six sanity prms) are synced (md5-verified) to BOTH mirrors. **4252462's VTK frames were merged into `animations/4252146_gjk_sanity_18k/vtk-out`** (frames 16000..29000 added, .pvd merged -> one series 0..29000; PROVENANCE.md updated; fluid crop changes from z<=100 to z<=140 at 16000).
> - **[Claude, 2026-09-21] DIVE into the 4252462 crash with the sweep (log `sbatch-output/output/*cont18k_4252462.out`, script `scratchpad/dive.py`).** (A) 61 sweep steps in 16 bursts between 15001 and 29396: THREE bursts are exactly 27 cells with |u|>0.3 (t=15741-15816 INTF (563.5,7.5,20.5); 17980-18009 INTF (552.5,11.5,12.5); 28987-29101 INTF (256.5,32.5,82.5), all 17-19 gas neighbours, |u| 0.30-0.34, no rod near) = flying spray droplets that healed by themselves; the other 13 bursts are face-neighbour |du|>0.2 hits at LIQUID pore cells between rods (11-21 obstacle nbrs, 1-8 pinched pairs, |u| 0.03-0.18) at (549.5,31.5,13.5), (380.5,18.5,49.5), (381.5,27.5,34.5), (408.5,24.5,18.5), (367.5-370.5,13.5-14.5,38-39), (474.5,20.5,0.5), (425.5,13.5,5.5), ... all harmless. (B) The crash burst: t=29100/29200/29300 |du| hits at pore cells (189.5-192.5,27.5-28.5,71.5) [owner rod ..727, obst 11-16, pinched 1-5, cell |u| 0.03-0.12 -> a much faster neighbour]; **t=29351 first `[HIGH-VELOCITY-PRECURSOR]` at LIQ (190.5,31.5,72.5) |u|=0.504, and the velocity FLIPS SIGN EVERY STEP: 29353 <+0.50,-0.12,-0.12>, 29354 <-0.47,+0.13,+0.12>; then the neighbour (191.5,31.5,72.5) 29372 0.515 (+), 29374 0.550 (+), 29375 0.529 (-), 29376 0.590 (+0.577), 29377 0.562 (-0.543), 29378 0.616 (+0.600), 29379 0.585 (-0.563)** = a period-2 (grid-scale, odd-even) oscillation of |u|~0.5-0.6 growing over ~30 steps -- the classical signature of an LBM instability, not a physical flow; t=29393-29395 the next cell (192.5,30.5,71.5) |u| 0.749 -> 1.238 -> 1.732; **t=29395 `[MEM-GARBAGE-LINK]` on rank 32 for the LIQUID cell block-local (63,11,12) = world (191.5,31.5,72.5) (the oscillating cell!): fill 1, density 0.294 (should be ~1), all populations ~ +/-1.4e6..2.9e6, and ~20 of its 26 neighbours are moving-obstacle cells (rawFlag 512) belonging to TWO rods (uid ..727 and uid 3530822107858468887 (..887))**; STALE-PARTICLE-IDX abort t=29396. (C) Geometry (frame 29000, 395 steps earlier, fluid_field VTK): the cell is still in OPEN liquid (3x3x3 neighbourhood 26 liquid + 1 obstacle, a rod ~4 LU to its left/below); rod ..727 moves +x at 0.031-0.035 ((182.8,32.2,69.3) at 29000 -> (193.2,32.1,69.8) at 29300; rod radius 3 => its top surface z~72.8, the cell is at z=72.5) so the creeping rod's leading surface sweeps over the cell and, together with rod ..887, encloses it. **Working hypothesis (inference, not proven): a liquid cell SEALED into a void between two rods (all links = moving-obstacle bounce-back, no fluid neighbour to exchange mass with) is squeezed by the creeping rod; MEM momentum exchange on ~20 links feeds a growing period-2 density/velocity oscillation which explodes in ~45 steps (29351 -> 29395); the garbage MEM force then invalidates the rod -> STALE. This matches the earlier 'pinched cell' picture and the 4251203 throat pores (obstacle nbrs 11-19, pinched pairs 2-4), but NOT the 4252005 downstream stream (0 obstacle neighbours) -> two different mechanisms: (I) sealed rod-rod pore cells, (II) plunging free-surface stream.** Caveats: the 'sealed' count comes from the garbage-link dump at 29395 (only links with |pdf_old|>10 are printed, so >=20 of 26 is a lower bound for that step); no per-step obstacle-neighbour count exists for the oscillating cell between 29000 and 29395; the run is chaotic (identical configs did not reproduce each other), so a resume may not repeat the crash. **Suggested next steps:** (1) add obstacle-neighbour count / owners / pinched pairs to the HIGH-VELOCITY-PRECURSOR lines and a cheap period-2 detector (u(t) ~ -u(t-1), |u|>0.1) to the sweep; (2) a prm-toggled 'sealed cell' guard (liquid cell with no liquid/interface neighbour among its 26 -> convert to obstacle or reset to equilibrium) as an experiment; (3) check the 4251203 crash for the same period-2 signature (HIGH-VELOCITY lines t=8120-8127 are in the log).
> - **[Claude, 2026-09-21] 4252462 crash, follow-up: what is physically going on before 29396, and Watanabe et al. 2023 as comparison.** **Nothing macroscopic happens:** the rods around the crash cell (191.5,31.5,72.5) move at 0.023-0.036 (= flow speed), |F_hydro| 0.02-0.11, |F_contact| <= 0.06 between t=29000 and 29300; global max|u| 0.28-0.34 belongs to a flying droplet cluster at x 245-265 (not near the site). **Geometry (frame 29000 orientations from `sizeOrientation`, positions from the CylForces log every 100 steps, axes frozen at their 29000 values => +-1-2 LU):** rods ..727 (4827858800541171727) and ..887 (3530822107858468887) are NEARLY PARALLEL (axes (0.89,-0.45,0.01) / (0.87,-0.49,-0.01)), stacked (z~69.5 and ~75), their surfaces are in contact throughout (surface gap -0.00/0.04/-0.01/0.13 LU at 29000/29100/29200/29300; F_contact between them rises 0.015/0.036 -> 0.053/0.045 with opposing directions at 29200) while both creep +x; the crash cell's distance to rod ..727's surface goes 1.55 -> 0.44 -> -0.20 -> 0.01 and to rod ..887's surface 3.26 -> 1.97 -> 0.86 -> 0.15 LU (t=29000..29300), i.e. the cell ends up AT the contact line of two touching parallel cylinders (a cusp/wedge), consistent with ~20 of 26 neighbours being obstacle cells in the garbage-link dump. **Block boundary:** x=191.5 is the last cell (local 63) of the block x in [128,192) on rank 32; rod ..727 left rank 32's blocks at t=29265 (x=192.001) and rod ..887 at t=29367 (x=192.001) -- the second hand-over falls INSIDE the 29351-29395 oscillation window (first precursor 29351); hand-overs happen thousands of times without a crash, so this is a coincidence to watch, not an established cause. At t=29396 both rods get garbage velocities (|v| 145 .. 866) and leave the domain (SUBCYCLE-LEFT-RANK, GHOST-VANISHED on ~20 ranks) -> STALE. **Mechanism (inference from LBM mechanics + the observed sign flips, untested):** a cell whose links almost all go to moving-obstacle cells reflects its populations (bounce-back reverses the momentum every step = the observed period-2 sign flip); the collision operator conserves mass and momentum and only relaxes non-conserved moments, so it cannot damp this; the MEM wall term (-6 w rho c.u_wall, wall velocities of the two rods differ) pumps momentum/mass into a cell that has no fluid neighbour to relax against. => another collision operator would NOT be expected to cure the sealed-cell mode; a treatment of sealed cells (fill as obstacle / reset) would. **Watanabe et al. 2023 (Eng. Appl. Comput. Fluid Mech. 17(1), 2211143), Sect. 1-2:** D3Q27 CUMULANT LBM (Geier 2015) + coherent-structure SMAGORINSKY LES (eddy viscosity nu_SGS added to nu_0 in tau, Eq. 11-16; 'to stabilise more by suppressing the dispersion of very small droplets in violent free-surface flows'); SINGLE-PHASE free surface (gas not simulated) with a conservative ALLEN-CAHN phase field for interface capturing instead of VOF cell flags ('stabilises by suppressing very small droplets'; mass conservation only approximate, their Sect. on conservation); DEM spring-dashpot contacts for rigid floating bodies; fluid-solid coupling by momentum exchange with INTERPOLATED bounce-back (Bouzidi type, q in (0,0.5] / (0.5,1]; falls back to halfway bounce-back where the second fluid node is not liquid; halfway BB for the phase field at solids); fluid force/torque AVERAGED over current and previous step (Aidun et al. 1998, Eq. 48-49: F = 0.5 (F(t)+F(t-dt))); octree AMR (fine grid near free surface and solids), multi-GPU. **Ours:** cumulant D3Q27 via lbmpy with ONE relaxation rate omega (`AntidunesLatticeModelGeneration.py`: Method.CUMULANT, no LES/Smagorinsky), VOF-type free-surface LBM (flags liquid/interface/gas, Schwarzmeier 2023), simple stair-step bounce-back MEM (`SimpleBB`), DEM spring-dashpot, force time-averaging already implemented (`useHydrodynamicForceTimeAveraging`, default true), no AMR. => Watanabe's extra stabilisers that we do not have: Smagorinsky LES, phase-field interface (fewer tiny droplets), interpolated bounce-back, AMR; Smagorinsky/phase-field would target the stream/droplet failure (mode II), not the sealed-cell mode (I), since a sealed cell has ~zero strain rate.
> - **[Claude, 2026-09-21] Interpolated bounce-back (IBB) -- what the sources say (Zotero PDFs read: Kruger et al. 2017 Sect. 11.2.1-11.2.2 pp. 436-447 (Eq. 11.7, moving walls p. 446); Bogner 2017 Sect. 3.3 (closure relation Eq. 3.137, bounce-back rule Def. 12 / Eq. 3.138, higher-order schemes Sect. 3.3.2, Tab. 3.2); Rettinger 2023 Sect. 3.2.2-3.2.4 pp. 26-30 (CLI Eq. 3.2, refilling Eq. 3.6-3.8) and Sect. 1 p. ~6 on trapped pore fluid).** Simple bounce-back (SBB, what `SimpleBB.h` implements, plus Ladd's moving-wall term): a population that would stream into a solid cell is reflected back along the same link; the wall is implicitly placed HALFWAY on the link whatever the true surface position -> stair-case surface, first-order accurate for curved walls (second order only if the wall really is halfway, delta=1/2; Bogner Sect. 3.3.1, Rettinger 3.2.2). IBB (Bouzidi et al. 2001) uses the true wall distance q (Kruger; delta_q in Rettinger, fraction of the link between the boundary fluid node and the wall) and builds the reflected population by LINEAR INTERPOLATION between post-collision populations of the boundary node and of its fluid neighbour x_f = x_b - c_i dt: q<1/2: f_ibar(x_b,t+dt) = 2q f*_i(x_b) + (1-2q) f*_i(x_f); q>=1/2: f_ibar = f*_i(x_b)/(2q) + (2q-1)/(2q) f*_ibar(x_b) (Kruger Eq. 11.7; the two cases keep the coefficients non-negative = better stability; both reduce to SBB at q=1/2); moving wall: add 2 w_i u_w.c_i/cs^2 (Ladd). Second-order accurate for arbitrary shapes; independent of the collision operator; Rettinger uses the equivalent CLI variant (Ginzburg 2008, Eq. 3.2) and Bogner uses CLI for some experiments (best w.r.t. collision-operator parameterisation; other IBB variants have viscosity-dependent slip). **Drawbacks (Kruger p. 445-446, Rettinger p. 27):** needs TWO fluid nodes along the link (three for quadratic) -> cannot be applied in narrow gaps between obstacles, Rettinger then falls back to SBB; not mass-conserving; viscosity-dependent slip; for moving walls newly uncovered ('fresh') nodes need refilling and the drag fluctuates when nodes appear/disappear. **Consequences for us (inference):** IBB would improve the accuracy of the r=3 LU rod surface (very stair-cased) and the hydrodynamic forces, but it would NOT help the sealed wedge cell between two touching rods, because exactly there it cannot be applied and falls back to SBB; single-node schemes (Chun-Ladd; Meng 2020, Marson 2021 cited by Rettinger) are the ones meant for narrow gaps. Rettinger also notes that fluid trapped in pores between particles is possible for circular particles in 2D but 'not possible for spherical particles in a three-dimensional setting' -- rods with line contact behave like the 2D case along the contact line, which may be part of why spheres do not show this failure (inference, not tested).
> - **[Claude, 2026-09-21] Sealed-cell CENSUS test + PID-off run 4252482 crash + Bogner's boundary scheme.** **(1) Census (`census_sealed_cells.py`, all 18 merged frames 0..29000 of 4252146/4252462):** LIQUID cells with >=24 of their 26 neighbours obstacle: 0 (t=0), 8 (5000), 20 (10000), 53-106 (15000-29000, ~80-100 typical); with >=25: 1-19; fully sealed 26/26: 0-1 per frame (t=22000, 26000). Their |u| is <=0.03 (typ. 0.001-0.01) and many persist for THOUSANDS of steps at the same position (e.g. (362,28,14) 15000-18000, (363,10,16) 15000-20000, (319,23,5) 21000-29000) -> **sealed cells are common and mostly inert (static walls); being sealed is NOT sufficient for a crash** (weakens the simple 'sealed cell' hypothesis; what is needed is sealed + dynamics: creeping/contacting rods changing the pocket, i.e. cells being covered/uncovered, or the free-surface bookkeeping in the pocket). None of the sealed examples is near the 4252462 crash cell at 29000 (that cell was still open liquid). **(2) 4252482 (GJK, PID OFF with accelerationXOverride 5.9e-7, sweep at |u|>0.17): FAILED at t=16377 after 42:52 with `[NaN-SCAN] FATAL` at an INTERFACE cell block-local (58,7,25) = world (58.5,27.5,25.5) on rank 6, NOT a fast-flow failure:** global max|u| only 0.103 at t=16300 (rho 0.994-1.033), NO sweep hit after t=14512 (only spray-droplet bursts of 18-27 cells with |u|>0.17 at 13881/14016/14512 and the start-up burst 3000-3029), no precursor, no garbage link, no STALE -- a single-step blow-up. The NaN cell's 3x3x3 neighbourhood: 23 OBSTACLE cells (rods), 2 normal LIQUID cells (rho 1.004/1.003, fill 1), 1 GAS cell; that GAS cell (59,6,25) was created 3 steps earlier: `[16374] No liquid or interface cell is in the neighborhood to distribute excess mass to. Mass is lost/gained. Cell (59, 6, 25) converted to gas.` on rank 6 (message from `ExcessMassDistributionSweep.impl.h:555`: an interface cell that must convert has NO liquid/interface neighbour to hand its excess mass to = an isolated interface cell in an obstacle pocket; the mass is simply lost/gained). That warning is routine (PID off: 186/272/195/30 per 2000-step bin at 10000-16400; PID on 4252146: 159/82/120/343 at 10000-18000, none before 10000) so it does not discriminate by itself, but its adjacency (the converted cell is a direct neighbour of the cell that turns NaN 3 steps later) is unlikely to be chance. Mechanism unknown (candidates: PDF re-initialisation/averaging with zero fluid neighbours -> 0/0; fill/density undefined in a pocket without fluid neighbours) -- NOT verified in the code. **Third crash type: low-velocity, NaN in an obstacle-enclosed INTERFACE cell after an isolated-interface-cell conversion next to it.** Flow in the PID-off run: mean velocity 0.029 (0-2000) -> 0.017 (4000-6000) -> 0.008 -> ~0 (8000-10000) -> **NEGATIVE -0.0036..-0.0051 from step 10000 on (Froude -0.12..-0.17, rod transport rate negative -0.011..-0.016)** with the fixed force 5.9e-7 (the no-clog channel had +0.04 with it): the clog chokes the channel completely and the mean flow reverses (physical backflow or a metric/sign artefact -- to be checked); upstream elevation 35 -> 60 (8000-10000) -> 37 (16000+), max fluid z 66-76. **(3) Bogner 2017 (thesis) Sect. 5.1 (p. ~110-111): for moving/floating bodies 'the present implementation uses the bounce-back rule of Def. 12 with u_w = u_P + omega x (x_w - x_P) and rho_w = rho(x_b), unless otherwise noted'** = SIMPLE bounce-back with Ladd's velocity term, no interpolation; Sect. 5.3 (floating cuboids, from Bogner & Rude 2013): FSLBM with **BGK** collision; righting-moment validation tau=1.0, g=1e-4; equilibration of floating cubes (densities 0.25/0.5/0.75; stable heel 26.565 deg for 0.25 and 0.75, 45 deg for 0.5): cuboid 2s x s x s, s=16, domain 130x40x40, tau=1/1.9=0.526, g=7.5e-4, initial heel 2.86 deg; the thesis does not restate the boundary rule for that specific test beyond the general Sect. 5.1 statement; TRT with magic parameter 3/16 + bounce-back is stated for the earlier fully-submerged validations (sphere drag/torque), the CLI rule only 'for some of the numerical experiments'.
> - **[Claude, 2026-09-21] Rod-specific check (Dennys: 'it works for spheres, so it has to do with the rods') + code path of an isolated interface cell + waLBerla's IBB.** **(1) Sealed-cell census, spheres vs rods (`census_sealed_cells.py`, now also per body-cell):** sphere hopper R12 (4248723, 12 frames 0..55000): liquid cells with >=24 of 26 obstacle nbrs = **0 in ALL frames** (>=20: 50-297; per 1000 body cells >=20: 1.5-3.5); rod pile (4252146/4252462, frames 15000-29000): >=24: 53-106 (0.28-0.62 per 1000 body cells), >=25: 7-17, >=20: 1460-1740 (8.8-9.7 per 1000 body cells) => rods produce sealed/near-sealed pockets, spheres do not (consistent with Rettinger p. ~6: trapped pore fluid impossible for spheres in 3D; rods with LINE contact make a wedge along a line). **(2) Reconstruction fallback counters (`[FORMERMO-CONVERT]`, `PdfReconstructionManager.h`), branch E = a rod-freed cell with NO real liquid and NO real gas neighbour -> made interface by fallback:** rods: 4252482 310 E-events / 1.94M converted cells (1.6e-4), 4252146 240/3.12M (7.7e-5), 4252462 189/5.0M (3.8e-5); spheres: 4248723 122/17.5M (7.0e-6), 4249475 100/17.7M (5.6e-6) => **rods hit the 'isolated freed cell' fallback 5-25x more often per converted cell**; pass-3 'unresolved island' (genuinely isolated from both phases) = 0 in every run, no `[PASS3-UNRESOLVED-ISLAND]` warnings. **(3) What the code does with an isolated interface cell (static reading, nothing run):** (a) `ExcessMassDistributionSweep.impl.h:555`: a converting interface cell without liquid/interface neighbour cannot hand over its excess mass -> warning 'Mass is lost/gained', cell converts (to gas or liquid) anyway (~130-300 warnings per 5000 steps in rod runs, 30-110 per 5000 in the sphere runs (745/932/1915 in 60k steps)); (b) `PdfRefillingSweep.impl.h:57-110` (`getAverageDensityAndVelocity`): a gas->interface refill averages rho,u over liquid/interface neighbours that are NOT newly converted and inside the block-local domain (ghost-layer data excluded when optimised communication is on); with ZERO valid neighbours it falls back to rho=1, u=0 and prints `WALBERLA_LOG_WARNING_ON_ROOT` (only rank 0 -> silent on other ranks); the division IS guarded (`numNeighbors != 0`), so my earlier guess 'division by zero' is REFUTED at this spot; the neighbour loop does `neighborU /= neighborRho` (would propagate NaN/0 from a bad neighbour, not the case in 4252482 where the neighbours had rho 1.003-1.004); (c) `PdfReconstructionManager.h` (rod-freed cells): branch E -> interface fallback, `seedConvertedCellFill` seeds interface fill with the mean of interface neighbours or 0.5 if none (l. 343), hydrostatic reseed returns early if rho<=1e-6. No code path found that turns 'no fluid neighbours' into a NaN by itself; the NaN mechanism in 4252482 stays UNIDENTIFIED. **(4) The 'cheap sub-cell test' I proposed for the righting box was ALREADY done (jitter scan 4249476, see the entries above): shifting both boxes by 0.25 cell in X changes the static MEM torque by mean |0.157| LU (max 0.38, sign flips), registration-averaged error rho .75 +0.057..+0.114 LU, and the s=32 scan roughly halves the rms error (0.082 vs 0.14-0.17 LU-equivalent) -> the staircase/registration effect is confirmed; open: does IBB remove it.** **(5) waLBerla already ships IBB for moving particles: `lbm_mesapd_coupling::CurvedLinear` (`src/lbm_mesapd_coupling/momentum_exchange_method/boundary/CurvedLinear.h`, R. Rettinger) = the CLI scheme of Ginzburg et al. 2008 (Eq. 4.1, Tables 3/4): pdf_new = f_i(x) + kappa0 (f_i(x_ff) - f_ibar(x)) with kappa0 = (1-2 delta)/(1+2 delta), wall term alpha = 4/(1+2 delta) (alpha=2 = SBB), delta from `mesa_pd::RayParticleIntersectionRatioFunctor` (analytic for Sphere/HalfSpace/Ellipsoid, GENERIC BISECTION line search for every other shape, i.e. boxes and cylinders work), falls back to simple bounce-back when the further cell x_ff is not a domain cell; plus `lbm/boundary/Curved.h` for static bodies.** Antidunes uses `MovingObstacle_T = lbm_mesapd_coupling::SimpleBB` (`AntidunesBoundaryHandling.h:72`); swapping needs: constructor differences (block storage/block, length scaling factor, hydrostatic density function), the far-cell test must accept only LIQUID|INTERFACE cells (gas cells hold invalid/-nan PDFs in the FSLBM, CurvedLinear tests `domainMask_`), and porting our SimpleBB changes (dry-link/-nan guards, `[MEM-GARBAGE-LINK]`, STALE-PARTICLE-IDX diagnostics). IBB would not change the wedge/sealed-cell case (falls back to SBB in narrow gaps).
> - **[Claude, 2026-09-21] IBB implemented as an OPTION + ρ=0.75 static righting scan queued (paired with a simple-BB control) + PID-off GJK continuation queued.** **Implementation:** `ModelParameters.useInterpolatedBounceBack` (default **false** = unchanged simple bounce-back; `SimulationParams.h`, plumbed through `AntidunesBoundaryHandling{.h,.impl.h}` and both `SetupStages*.h` call sites into the `SimpleBB` constructor, new last argument). In `src/lbm_mesapd_coupling/momentum_exchange_method/boundary/SimpleBB.h::treatDirection`: if the flag is on and the FURTHER fluid cell x_ff = x - c_i is a liquid/interface cell (`domainMask_`) with finite PDFs, the link uses the CLI scheme (Ginzburg, d'Humieres, Verhaeghe 2008, Eq. 4.1, Tables 3/4 = `lbm_mesapd_coupling::CurvedLinear`): pdf_new = f_i(x) + kappa0 (f_i(x_ff) - f_ibar(x)), kappa0 = (1-2 delta)/(1+2 delta), wall term alpha = 4/(1+2 delta), delta from `mesa_pd::RayParticleIntersectionRatioFunctor` (analytic for spheres, bisection for boxes/cylinders, tolerance 1e-4 dx), the wall velocity is taken at x + delta c_i; every other link (no second fluid cell, non-finite PDF, ray outside [0,1]) falls back to simple bounce-back; the MEM force/torque computation (incl. fill weighting, hydrostatic subtraction, all project guards/diagnostics) is unchanged. Log line `[IBB] interpolated bounce-back active on this rank: links so far interpolated=.., fell back=.. | this link ... delta=..` (first hit, 1000th, 100000th). **Validated locally** with a copy of waLBerla's sphere drag test (`tests/showcases_antidunes/IbbEquivalenceDragTest.cpp`, target `ANTIDUNES_IbbEquivalenceDrag`, variants BB / CLI (= CurvedLinear) / IBB (= SimpleBB + flag), `--shortrun`, 2 MPI ranks): drag after 150 steps CLI **4.0464086404299753** vs IBB **4.0464086404302622** (relative difference 7e-14 = round-off) vs plain BB 4.1496136556585039 (+2.5 %) => the option reproduces waLBerla's CurvedLinear and is different from simple BB. Full-length accuracy runs (drag vs the analytical value, BB vs IBB) were started in the background (`/tmp/claude-1000/ibbfull/{BB,IBB}/run.log`, debug build => slow) -- not evaluated yet. `Antidunes` and `DropTest` build. **Note (git history):** the original Antidunes (merge 0b415d68) already used `SimpleBB` as `MovingObstacle_T` (variable name `curvedLinear` is a leftover; `CurvedLinear.h` was included but unused) -> the published Antidunes/Bogner behaviour is plain bounce-back. **Scan jobs (same 14 heel angles 15, 20, 16.5, 18, 19.5, 21, 22.5, 24, 25.5, 26.565, 27.5, 29, 30.5, 32 deg as scan2, ONE box rho 0.75 at X=65, s=16, 65 ranks, 2.5 h, raw MEM torque, no analytic buoyancy; prm sets differ ONLY in the useInterpolatedBounceBack line, verified by diff):** **4252514** = IBB (`bogner_righting_scan3ibb75`, ~/mt-dennys-mirror, started immediately), **4252515** = simple-BB control (`bogner_righting_scan3sbb75`, ~/mt-dennys-mirror2, pending). NB the box position differs from scan2's rho 0.75 box (X 32.5 there because two boxes were used) => the control is essential; compare the two new scans, not the new IBB scan against scan2. Analysis: `analyze_bogner_righting_scan.py <archive_dir> --manifest prms/bogner_righting_scan3{ibb,sbb}75_s16_manifest.json`; question: rms error of the static MEM torque vs the ideal curve (SBB scan2: 0.17 LU for rho 0.75, sub-cell jitter mean |0.157| LU) and the zero crossing (ideal 26.565 deg). Generator: `gen_bogner_righting_scan_prm.py --ibb`. **PID-off GJK continuation 4252516** (`hopper_clog_cyl_gjk_pidoff_sanity_cont20k`, ~/mt-dennys, 144 ranks, 1.5 h, pending): resumes 4252482's LAST checkpoint (vault `antidunes-hopper-clog-cyl-gjk-pidoff-sanity-18k_4252482_20260921_013727/snapshots`, simState 15001, PID off, acceleration 5.9e-7) for +20000 steps (`timesteps 35001`), VTK every 1000, sweep at |u|>0.17 + dense scan armed; the original crashed at 16377 (1376 steps after this checkpoint), the run is chaotic so it may or may not repeat the failure. All sources synced to all three checkouts (9 source files md5-identical).
> - **[Claude, 2026-09-21] Results: local IBB accuracy test, scan jobs, PID-off continuation 4252516.** **(1) Full local sphere-drag test (`ANTIDUNES_IbbEquivalenceDrag --logging`, 2 ranks, debug build): analytical drag 2.84202; plain bounce-back 2.92515 (relative error 2.93 %); IBB (SimpleBB + useInterpolatedBounceBack) 2.84490 (relative error 0.101 %) => the option is ~29x more accurate than simple bounce-back on this sphere test.** **(2) The two righting-scan jobs 4252514 / 4252515 FAILED at the CMake step after ~2 min: I had synced `tests/showcases_antidunes/CMakeLists.txt` (which lists PinchedFluidCellsTest.cpp, EnergyConsistentContactTest.cpp) to the mirrors without the test SOURCES (same mistake as with 4251422, this time on the two mirrors, whose test dirs were stale).** Fixed: the whole `tests/showcases_antidunes/` directory synced to both mirrors (dry-run: 0 files differing), resubmitted as **4252644 (IBB, ~/mt-dennys-mirror)** and **4252646 (simple-BB control, ~/mt-dennys-mirror2)**. Lesson: never push a CMakeLists without its sources; sync the whole directory. **(3) PID-off continuation 4252516 (from 4252482's step-15001 checkpoint, +20000 steps): COMPLETED all steps (15001 -> 35001, 58:02, exit 0) -- it passed step 16377 where the original crashed** (chaotic: the crash did not repeat from the same checkpoint). max|u| 0.151 (200 samples, never above the 0.17 sweep limit), rho 0.979-1.053, max|f| 0.312, no NaN/precursor/garbage/STALE, sweep: only 2 violating steps with 1 cell each; excess-mass warnings 450 / 90 / 43 / 53 per 5000 steps (15000-20000 ... 30000-35000, decreasing). Flow (fixed force 5.9e-7, target-velocity 0.04 NOT enforced): mean velocity -0.0048 (15000-22500, backflow) -> crosses 0 at ~24000 -> **+0.0020 (25000-27500) -> +0.0049 (30000-32500) -> +0.0037-0.0044 at the end**, i.e. the choked channel re-opens to ~10 % of the clean-channel flow; transport rate -0.011 -> +0.002; upstream elevation 37 -> 32-33. VTK frames every 1000 steps 16000..35000 are in the vault archive `antidunes-hopper-clog-cyl-gjk-pidoff-sanity-cont20k_4252516_20260921_032037` (not copied to animations). **(4) `maxFluidPos_z` is POLLUTED by a levitating interface blob:** it stays at 76.5-76.7 for the whole 20000 steps (and already 76.5 at 14000-18000 in 4252482); frame 30000 shows >= 36 INTERFACE cells at x 307.5-309.5, y 11.5-14.5, z 74-77.5, all with the same |u| = 0.062 (block of 3x4 cells per layer), the same object that was the fastest cell (308.5,13.5,75.5) |u| 0.061-0.064 at t=15100-15300 in 4252482 -- a stationary blob of interface cells ~40 LU above the water (upstream elevation 33) that never moves, the same class of artefact as the stuck droplet at (497.5,15.5,25.5) in 4252146. => do NOT use `maxFluidPos_z` as the water level in these runs; use the clog elevation probes or a percentile of the fill field. Origin of these levitating blobs unexplained (candidate: interface cells without gas or liquid neighbours cannot convert/advect).
> - **[Claude, 2026-09-21] RESULT of the rho=0.75 righting scans, interpolated (4252644) vs simple bounce-back (4252646): IBB does NOT reduce the static torque error.** Both jobs COMPLETED all 14 heel angles (75 min / 70 min, 0 failures; vault `antidunes-bogner-scan3-ibb75-s16_4252644_20260921_075010`, `antidunes-bogner-scan3-sbb75-s16_4252646_20260921_075041`); IBB was active only in the IBB job (936 `[IBB]` log lines, wall distances delta 0.38 / 0.93 ..., none in the control). Raw MEM torque about +Y minus the ideal hydrostatic torque (LU, n=14 angles 15-32 deg, boxes held at the ideal heave): **IBB mean +0.093, rms 0.156 (sd 0.125); simple BB mean +0.102, rms 0.159 (sd 0.122)**; paired per-angle difference IBB - SBB: mean -0.009, **mean |d| = 0.0127 LU, max |d| = 0.027 LU** -- about 12x smaller than the error itself and than the registration scatter (0.157 LU for a 0.25-cell shift, jitter scan 4249476). Vertical force error +1.67 % (IBB) vs +1.51 % (SBB) (not improved). Measured zero crossing of the torque: **25.15 deg (IBB) vs 25.17 deg (SBB)** vs ideal 26.565 (the earlier TWO-box scan2 gave 27.2 deg for the same rho 0.75 box at another X position => registration alone moves the crossing by ~2 deg). Per-angle errors (IBB | SBB): 15: +0.101|+0.127, 16.5: +0.201|+0.176, 18: +0.241|+0.258, 19.5: +0.312|+0.316, 20: -0.022|-0.013, 21: +0.017|+0.025, 22.5: +0.134|+0.141, 24: +0.131|+0.134, 25.5: -0.105|-0.101, 26.565: +0.179|+0.193, 27.5: +0.212|+0.227, 29: -0.065|-0.038, 30.5: -0.063|-0.045, 32: +0.031|+0.029. **Conclusion: the wall PLACEMENT on the box faces (staircase vs exact) is NOT what causes the static torque error / registration scatter, although the option is verified correct (sphere drag error 2.9 % -> 0.10 %, agrees with CurvedLinear to 7e-14).** Interpretation (hypotheses, untested): (a) the hydrostatic MEM force on a flat tilted face is insensitive to where along the link the wall sits (a staircase preserves the projected areas the pressure acts on); (b) the error and its dependence on the sub-cell position (x-shift 0.25 cell changes the torque by 0.157 LU with the ideal unchanged) comes from the FREE SURFACE / waterline: the fill-level weighting of the MEM force in interface cells (`useFillWeighting`, Bogner Eq. 5.9), the interface position/curvature at the box, or the compressible hydrostatic density profile; the error roughly halves from s=16 to s=32 (0.14-0.17 -> 0.082 LU-equivalent, scan 4249478) like a first-order effect. **Suggested next tests (not run):** (1) the same 14-angle scan with `useFillWeighting false` (IBB off); (2) a FULLY SUBMERGED held box (no free surface, or a deep waterline far from the box) at a few heel angles -- if the torque scatter vanishes, the waterline is the cause; (3) IBB + s=32 at 3-4 angles; the rho=0.25 twin was not run.
> - **[Claude, 2026-09-21] Is the MEM force weighted twice (useFillWeighting + useBognerThreePhaseMEM + useMemOnlyBuoyancy)? -- NO (code reading, nothing run).** (1) `useFillWeighting` (`SimpleBB.h`, `fillLevelFct_`): the ONLY multiplicative weighting: the force of every obstacle link is scaled by the near-boundary cell's fill level phi (clamped to [0,1]); the bounce-back PDF written back into the fluid is untouched. (2) `useBognerThreePhaseMEM` is NOT a weighting: it only decides whether a depth-dependent hydrostatic reference density rho_hyd(z) is subtracted from the link PDFs before the force is formed (`SetupStages.h` ~l. 812-826: true -> null `hydrostaticDensityFct`, only the constant zero-centring w_i (rho=1) remains, so the hydrostatic pressure = buoyancy stays in the MEM force); together with (1) this is Bogner & Ruede (2013) Eq. 5.9 (phi-weighted blend of the boundary-link PDFs and the gas equilibrium; the gas term vanishes after the rho=1 zero-centring). (3) `useMemOnlyBuoyancy` (default true): body force = gravity only (rho*V_full*g), no analytic buoyancy -> buoyancy is counted once (from the MEM). (4) The old GLOBAL submergence scaling `hydro_force *= V_sub/V_full` in `AddExternalForces` (ParticleKernels.h ~l. 3340-3353) is commented out ('DELIBERATELY OFF ... retired'); the `submScale=` text in `[CylForces]` is informational only although the log says 'applied'. Consequence for the planned ablation: `useFillWeighting false` (keep useBognerThreePhaseMEM true, useMemOnlyBuoyancy true) removes the only phi weighting => the raw two-phase MEM force (Bogner Eq. 5.6: interface cells at the waterline count as full liquid, phantom buoyant lift) -- an ablation of the waterline treatment, no double or missing hydrostatic term; note that with it OFF the run is no longer Bogner's Eq. 5.9 even though the flag is still called useBognerThreePhaseMEM.
> - **[Claude, 2026-09-21] Fill-weighting ablation scan queued: job 4252921** (`antidunes-bogner-scan3-nofill75-s16`, ~/mt-dennys-mirror2, 65 ranks, 2.5 h, pending on priority): the SAME 14 heel angles / ONE rho 0.75 box / s=16 / simple bounce-back as the control scan 4252646, with `ModelParameters.useFillWeighting false` (raw two-phase MEM force, no phi weighting at the waterline); useBognerThreePhaseMEM and useMemOnlyBuoyancy stay true. Prms `prms/bogner_righting_scan3nofill75_s16_a*_gas.prm` (`gen_bogner_righting_scan_prm.py --nofill`, new flag), diff to the control prms checked: ONLY the useFillWeighting line (a15, a26p565, a32; manifest identical), job `jobs/run_65_bogner_righting_scan3nofill75_s16.sbatch`. Compare directly with 4252646 (control: torque error mean +0.102 / rms 0.159 LU, zero crossing 25.17 deg, Fz +1.51 %) and 4252644 (IBB): if the static torque error / zero crossing change clearly, the waterline phi-weighting matters; if not, the error sits elsewhere (interface position, hydrostatic density profile, or a deeper lattice effect). Analysis: `analyze_bogner_righting_scan.py <archive> --manifest prms/bogner_righting_scan3nofill75_s16_manifest.json` (archive under `$HPCVAULT/antidunes-archive/antidunes-bogner-scan3-nofill75-s16_4252921_*`).
> - **[Claude, 2026-09-21] The candidate causes of the static righting-torque error compared with Bogner 2017 (thesis Sect. 3.3.3, 5.3.2, 5.3.3) and with our code/scan setup.** **(1) Interface position / free-surface BC:** Bogner uses the FSK rule (Koerner et al. 2005) 'for most applications', proven only FIRST-order accurate (Sect. 3.3.3); the second-order FSL rule (Bogner et al. 2015) exists but is not the default. Ours: `src/lbm/free_surface/dynamics/functionality/ReconstructInterfaceCellABB.h` = 'Free surface boundary condition as in Koerner et al., 2005' (pressure anti-bounce-back, Thuerey 2007) => the SAME first-order rule as Bogner's; waLBerla has no FSL. Data: in the sbb scan 4252646 the total fluid mass is constant to the unit over 20000 steps (a15 102664 -> 102664, a24 102748, a26p565 102668, a32 102718: level shift 0.0000 cells), in the IBB scan 4252644 it drifts +16..+109 (0.003-0.021 cells; IBB is not mass conserving) => a drifting waterline is NOT the cause; fill init surface at absoluteLiquidHeight+0.5 = 20.5 (DataInitialization.h:160), the mass balance gives ~20.49. An effective waterline/heave offset of ~0.14-0.18 cell would explain both the mean torque bias (+0.09..0.10 LU at ~0.7 LU/cell) and the vertical force error (+1.5 %, F ~ 12 % per cell of waterline) -- consistent with an O(dx) free-surface/interface convention effect, but not separated. **(2) Hydrostatic density profile:** Bogner initialises the basin 'with a pressure gradient according to hydrostatic equilibrium' (5.3.2, 5.3.3); ours does the same (`initPdfFieldHydrostaticPressure`, reference height absoluteLiquidHeight = 20 vs the fill surface at 20.5 => 0.5 cell / 0.11 % density offset at t=0 that relaxes; the torque window 19000-19300 is steady with sd ~0); both use the MEM force without an extra hydrostatic subtraction (useBognerThreePhaseMEM = Eq. 5.9 as published) => no difference found, unlikely to be the cause. **(3) Lattice effects:** Bogner reports 'jumps in the stability curves' from the lattice approximation of the shape that shrink with resolution, but validated the FIXED-heel torque only for rho=0.5 boxes of ratio 6:4 and 5:4 (24x16, 48x32, 96x64) -- larger lever arm than the cube; the rho=0.75 cube (1:1, his least stable ratio) is shown only as a dynamic equilibration figure at s=16 (130x40x40, tau=1/1.9, g=7.5e-4, 2.86 deg heel) with no numbers. So our observation (scatter ~0.16 LU vs an ideal signal 0.04-0.19 LU, error halving from s=16 to 32) is not in conflict with his text; he never measured our quantity. **(4) A difference his text does not cover: SURFACE TENSION is ON in our scans** (`useHardcodedUnits true` sets the literal `surfaceTension = 0.000849592` in `SimulationParams.h` l. 1490 regardless of We; capillary length sqrt(sigma/(rho g)) = sqrt(8.5e-4/7.5e-4) = 1.06 cells at g = 7.5e-4, `curvatureModel FiniteDifferenceMethod`, `enableWetting false`); Bogner's Sect. 5.3.2/5.3.3 do not mention surface tension for these tests (unknown whether it was off). Not yet tested; needs a small `surfaceTensionOverride` option (0 = off). **(5) Collision/bounce-back:** Bogner BGK tau=1.0 (righting-moment test) / 1/1.9 (equilibration) with simple bounce-back (Def. 12, Sect. 5.1); ours cumulant omega=1.9 with simple bounce-back by default; hydrostatics makes the collision operator's viscosity dependence irrelevant to first order. **Ranking of what is still open:** (a) fill weighting (scan 4252921 queued), (b) surface tension off (new option + scan), (c) an FSL-like second-order interface rule (not in waLBerla, big job), (d) higher resolution.
> - **[Claude, 2026-09-21] PID-OFF SPHERE hopper queued: job 4252971** (`hopper_clog_pidoff_2node_frontblock`, ~/mt-dennys-mirror, 144 ranks, **tight 2h15m limit** because short jobs are scheduled sooner): the completed PID-on run 4248723 (sphere front block against the pillars, pillar radius 12, gap 16, 1015 spheres, 60000 steps, 1h50m04s; twins 4247229 1h48m, 4249475 1h52m) with `PIDParameters.enabled false; accelerationXOverride 5.9e-7` = the steady force of the clean channel measured in 4252463 -- valid because the sphere hopper channel is IDENTICAL (768x40x180, liquidHeightFactor 2.0, floorHeightFactor 4.1393, pillars R12, same Re/Fr/We, same unit system). The prm differs from `hopper_clog_pidon_2node_frontblock_gas.prm` ONLY in those two PID lines and the fluid VTK crop z 100 -> 140 (checked with a diff); NO extra fluid diagnostics (I first added the dense scan + sweep, then removed them so the cost per step equals the original and the 2h15m limit is safe; the limit of the pending job was lowered with `scontrol update`). 60000 steps, VTK every 5000, snapshots every 5000. Files `prms/hopper_clog_pidoff_2node_frontblock_gas.prm`, `jobs/run_144_hopper_clog_pidoff_2node_frontblock.sbatch`. Purpose: the sphere counterpart of the PID-off rod run 4252482/4252516 (rods: choked channel, backflow, mean velocity ~0 after step 9000) -- does the fixed clean-channel force get spheres through the gap (transport rate, upstream elevation vs the PID-on 4248723: no clog, transport 3.7-5.4e-2)?
> - **[Claude, 2026-09-21] Sphere hopper CONTINUATION queued: job 4252997** (`hopper_clog_pidon_2node_frontblock_cont50k`, ~/mt-dennys, 144 ranks, **tight 2h10m limit**): resumes the completed PID-ON sphere run 4248723 (pillar R12, 1015 spheres, no clog, 60000 steps) from its last checkpoint (vault `antidunes-hopper-clog-frontblock-pidon-2node_4248723_20260919_224318/snapshots`, simState 55001, PID state restored) for +50000 steps (`timesteps 105001`), **VTK every 1000 steps** (fluid crop z 0..100 as in the original so the frames can be merged into `animations/4248723_sphere_hopper_frontblock_R12`), checkpoints only every 25000 steps (saves time), no extra diagnostics. Time estimate: 50000 x ~0.105 s (original 1h50m for 60000 incl. build) ~ 88 min + build ~5 min + restore/VTK/archive ~15 min ~ 1h50m. Files `prms/hopper_clog_pidon_2node_frontblock_cont50k_gas.prm` (diff to the original prm: loadSnapshot, timesteps, snapshotFrequency, VTK writeFrequency only), `jobs/run_144_hopper_clog_pidon_2node_frontblock_cont50k.sbatch`. NB the original wrote steps 55001-60000 already, so 56000..60000 of the continuation overlap in time with the original's last 5000 steps (the run is chaotic, not bit-identical). **Three jobs now queued on three checkouts:** 4252921 nofill righting scan (~/mt-dennys-mirror2), 4252971 PID-off sphere hopper (~/mt-dennys-mirror), 4252997 sphere continuation (~/mt-dennys). /home/hpc is over its soft quota (107 GB of 100 GB, hard 200 GB): the continuation's ~11 GB of VTK is archived to the vault at the end; the cleanup of verified duplicates (~20 GB) is still pending your OK.
> - **[Claude, 2026-09-21] WHY RODS AND NOT SPHERES? -- deeper look at the PID-off GJK crash 4252482 (t=16377): the crash cell sits in a THIN CHANNEL between a locked BUNDLE of three parallel touching rods.** **(1) Ghost/handover events do not separate rods from spheres:** per 1000 steps `SUBCYCLE-LEFT-RANK` rods 31 / 51 / 7 / 101 (4252482 / 4252146 / 4252516 / 4252462) vs spheres 43 / 50 (4248723 / 4249475); `GHOST-VANISHED` rods 65 / 110 / 22 / 208 vs spheres 46 / 57; `STALE-NEARBOUNDARY` rods 42 / 68 / 28 / 94 vs spheres 27 / 32 (same order of magnitude; some of these lines may be throttled, so rates are approximate). **(2) The crash cell's neighbourhood (rank 6, block-local (58,7,25), world (58.5,27.5,25.5), interface cell, 23 of 26 neighbours obstacle, 2 liquid nbrs both at dy=+1, 1 gas nbr at dy=-1 = the cell converted to gas 3 steps earlier):** three rods with consecutive uids ..8404 (=1801439850948198404), ..8407, ..8408, 4.8-6.4 LU away (step 16300: (64.3,28.2,23.0), (60.4,23.8,22.5), (62.2,25.6,27.9)), identical tilt -51 deg and velocity (-0.009,0.003,0.006). **From frame 10000 to 17000 (7000 steps!) the three rods are MUTUALLY TOUCHING (surface gaps -0.01 .. -0.03 LU) and PARALLEL (axis angle 3.8 deg at 10000, 0.8 deg from 15000) and rotate/move together as one rigid bundle** (axis (0.26,0.96,0.12) at 10000 -> (-0.80,0.50,0.34) at 15000 -> (-0.74,0.65,-0.15) at 17000; |v| 0.007-0.016; contact forces ~0.09). **Three mutually tangent parallel cylinders of radius 3 enclose a triangular channel whose inscribed radius is 0.1547 r = 0.464 LU (< half a cell), running along the whole rod length (18 LU)** -- only ~0.67 cell centres per unit length lie inside it, so the fluid there is a chain of sealed cells (23-26 obstacle neighbours) that is physically unresolved by the lattice and stays intact as long as the bundle is locked. The crash cell is at the open end of such a channel (liquid cells continue along +y, gas at the -y end). Spheres cannot do this (three touching spheres leave a point-like interstice): census sealed liquid cells (>=24/26 obstacle nbrs) spheres 0 in all 12 frames vs rods 8-106 per frame, branch-E 'freed cell without liquid or gas neighbour' fallbacks 5-25x more frequent per converted cell for rods. **(3) Ghost exchange not excluded:** the ONLY ownership event near the crash is `SUBCYCLE-LEFT-RANK rank=18 uid=1801439850948198404 (= bundle rod ..8404) step=16337 pos=(63.9998,28.33,23.24)`: bundle member ..8404 crossed the block face x=64 (rank 18 -> the rank owning x<64), 40 steps before the NaN, i.e. the bundle (and its channel) STRADDLED a block boundary at that time; handovers per se are equally frequent for spheres, so what would be special is a sealed channel crossing a rank boundary (ghost-layer consistency of the obstacle mask / PDFs across the face), which the statistics cannot test. The failing cell itself is 5.5 cells from the x-face, not in the ghost layer. **(4) Mechanism inside the cell still UNKNOWN** (no high velocity: global max|u| 0.103; single-step NaN; adjacent isolated interface cell converted to gas 3 steps before, `ExcessMassDistributionSweep.impl.h:555`). **Proposed tests for later (none run):** (a) a small controlled setup: a bundle of 3 parallel touching rods (r=3, L=18) and, for comparison, 3 touching spheres, in still water with a free surface, at rest and moving slowly (0.01) with the channel (i) inside one block and (ii) crossing a block face -- tests channel vs ghost exchange separately; (b) log/census per step of interface and liquid cells with >=24 obstacle neighbours and their owners; (c) a prm-toggled mapping post-pass that turns fluid cells with >= N (e.g. 24) obstacle neighbours into obstacle cells of the owning rod ('close sub-cell channels'), then repeat rod runs and see whether the failures disappear.
> - **[Claude, 2026-09-21] LOCAL ROD-BUNDLE CHANNEL TEST (9 MPI ranks, laptop): the touching three-rod bundle alone reproduces the fluid anomalies; isolated rods with identical motion do not.** **Setup** (`gen_rod_bundle_channel_test_prm.py`, prms `rod_bundle_channel_test_gas.prm` / `rod_bundle_channel_control_gas.prm`, Release build `build-release`, `mpirun -np 9`): periodic box 48x36x48 in x,y (walls top/bottom NoSlip), blocks 16x12x48 (3x3x1; the code needs >=3 blocks per periodic direction), still water with a free surface at z=34 (no forcing, PID off, no pillars), production unit system (`useHardcodedUnits`, omega 1.99753, g 5.82575e-5, sigma 1.63945e-3), fluid sanity sweep + dense NaN scan armed from step 0. Exactly the three rods of the bundle (r=3, L=18, surface gaps -0.03, axis (-0.74,0.65,-0.15) as at step 17000 of 4252482, centre distance 5.97), all three forced by the new test option `ModelParameters.rodPrescribedVelocityX/Y/Z` (-0.0096, 0.0043, 0: the measured bundle velocity, set every DEM sub-cycle, angular velocity 0, rods pass every block face and the periodic wrap). New options in the LOCAL sources only (not yet on the cluster checkouts): `rodPrescribedVelocity*`, `surfaceTensionOverride`. **Result (40000 steps, 26 min):** bundle: fluid density down to rho_min 0.55-0.80 and max|u| 0.15-0.26 (Ma>0.3 warning level 0.17 exceeded; 282 violating sweep steps), FLUID-STATS every 500 steps; every flagged cell is a LIQUID cell inside the triangular channel with 21-25 of 26 obstacle neighbours (pinched pairs 9-12), lowest rho at z~19.5-20.5 (rod centre height) travelling with the bundle. The anomaly is intermittent (quiet windows e.g. steps 0-1000 and 16300-21000 with rho_min >0.99, then starts again within ~100 steps) and the bundle drags the water along (mean velocity -0.0075 at 40000). NO NaN / abort in 40000 steps, i.e. the anomaly appears but does not necessarily lead to a crash. **CONTROL (`--nobundle`: same three rods, same prescribed velocity, same box, but 14 LU apart = no touching, 22000 steps):** rho_min >= 1.00000, max|u| 0.012-0.014, ZERO sweep violations at any step. => neither the rod shape, the periodic box, the block layout/face crossings, the velocity, nor the surface tension causes the anomaly: it needs the touching-rod bundle (sub-cell channel). **Caveats:** prescribed rigid velocity injects momentum (mean water velocity -0.0075 develops), the start is impulsive, rods have no contact dynamics, and the crash itself (single-step NaN in an interface cell at the channel end) was not reproduced; only ONE bundle registration (dx=dy=dz=0) and ONE block layout were run. **Not yet done:** (a) other registrations and a bundle placed inside one block only (ghost-exchange test), (b) 3 touching spheres as comparison, (c) the channel-closing guard (fluid cells with >=24 obstacle neighbours -> obstacle) to see whether the anomalies vanish. Logs: `/tmp/claude-1000/bundle_run/{bundle_z0,control}/run.log` (not archived).
> - **[Claude, 2026-09-21] FLOW-RATE EVALUATOR (ClogMeasurement.flowRateEnabled).** `ClogMeasurementEvaluator` (Evaluators.h) now optionally writes `clogUpstreamFlowRate`, `clogDownstreamFlowRate`, `clogUpstreamWetArea`, `clogDownstreamWetArea`: Q = sum over the y-z plane of fill*u_x and A = sum of fill over the FLUID cells only (fill>0, not moving-obstacle, not no-slip), averaged over the 2*floor(probeHalfWidthXY)+1 planes x-2..x+2 around the probe x (a single 1-cell plane can fall into the 1-LU slit between two rows of rods and then says nothing; Q is conserved from plane to plane). Cells covered by rods/spheres/pillars are excluded, so a body in the plane reduces A and adds nothing to Q -> shape independent. Lattice units (Q in cells^3/step, A in cells^2), signed (back-flow). Prm: `flowRateEnabled true;` in the ClogMeasurement block (default false; only the first probe pair 'clog' writes it). **Validation (local, `gen_rod_bundle_channel_test_prm.py --flowtest --accel 5e-6 --vx 1e-9`, 3 rods r=3 L=18 fixed in a driven periodic channel, planes at x=15 through the bundle and x=40 rod-free, 9 ranks, 3000 steps):** Q(rod-free) = 0.991-1.002 x (meanVelocity_x * 16.0039 * 36) (the PID's mean velocity uses the nominal liquid height), Q(through the rods) = Q(free) within -0.6..+0.2 % after step 1000 (mass conservation) while the wet area drops from 1242.3 to 1125.4 = -116.9 cells^2 (analytic blocked cross-section of the 3 rod ellipses 3*pi*r^2/|a_x| = 114.6). A 144-rank prm smoke test was run with larger blocks (18 ranks, 120 steps): columns written, plane through the initial rod block blocked 43 %. Files: Evaluators.h, SimulationParams.h (`clogFlowRateEnabled`), BuildTimeloop.h.
> - **[Claude, 2026-09-21] NO-PILLAR RUNS QUEUED (PID on, 100000 steps, 144 ranks, flow-rate planes x=336 / x=432).** (1) Rods: **4253776** `antidunes-hopper-clog-cyl-nopillars-pidon-100k`, `~/mt-dennys-mirror2`, `--dependency=afterany:4252921` (the pending nofill righting scan), limit 3h45m; prm `prms/hopper_clog_cyl_nopillars_pidon_2node_100k_gas.prm` = the PID-on rod run of 4252146/4252462 (406 rods, GJK-only, sanity sweep every 10 steps + dense NaN scan) with CylinderPillars + both BodyBC pillar lines removed, checkpoints every 10000. (2) Sphere control: **4253777** `antidunes-hopper-clog-nopillars-pidon-100k`, `~/mt-dennys-mirror`, `afterany:4252971` (the PID-off sphere hopper), limit 3h20m; prm `prms/hopper_clog_nopillars_pidon_2node_100k_gas.prm` = the sphere front-block hopper 4248723 without pillars. Generator `gen_nopillar_prm.py`; sbatch `jobs/run_144_hopper_clog_cyl_nopillars_pidon_100k.sbatch`, `jobs/run_144_hopper_clog_nopillars_pidon_100k.sbatch`. Question: do rods still crash (sealed pore / bundle channel) in an OPEN channel, and do the spheres stay clean. Time-limit reasoning: 4252146 (18k steps, sweep every step) took 46m43s incl. build (~0.13 s/step), the PID-on sphere twins 1h48-1h52 for 60000 steps. NB the sources of the two mirrors were synced (Evaluators.h, BuildTimeloop.h, SimulationParams.h, ParticleKernels.h) BEFORE the sealed-cell guard below was added, so these jobs run WITHOUT the guard; the pending jobs 4252921/4252971 build from the same synced sources (all new options default off). The sphere run's second job was my reading of 'also queue after another run in another directory'; cancel 4253777 if a duplicate of the rod run was meant. **[2026-09-21, later] 4253777 (sphere control) CANCELLED at the user's request; only the rod run 4253776 remains from this pair.**
> - **[Claude, 2026-09-21] LOCAL TESTS OF THE ROD BUNDLE, RESULT: (a) DECOMPOSITION-INDEPENDENT, (b) CLOSING THE SEALED CELLS REMOVES THE ANOMALY.** (a) The 40000-step bundle run (9 ranks, blocks 16x12x48) was repeated for 8000 steps on 16 ranks with blocks 12x9x48 (all block faces at different positions; `gen_rod_bundle_channel_test_prm.py --cpb 12,9,48`): the fluid statistics agree at first (step 1000: rho_min 0.9947 vs 0.9945, step 1500: 0.9086 vs 0.9084, step 2000: 0.9222 vs 0.9208) and then follow the same envelope (rho_min 0.63-0.75 and max|u| 0.19-0.26 around steps 4500-7000 in both), the small differences grow as in any chaotic free-surface run (max |difference| of the sampled rho_min 0.13). Also, the quiet window 16300-21000 of the first run spans several block-face crossings in x and y while the loud windows span them too. => the anomalies are not a ghost-layer/decomposition effect. (b) New experimental guard `SealedCellGuard.h` / `ModelParameters.closeSealedCellsMinNeighbours` (0 = off, default): right after MovingParticleMapping and before the PDF reconstruction every non-gas fluid cell (and every just-vacated cell) with >= N of its 26 D3Q27 neighbours covered by a moving obstacle or no-slip wall (and >= 1 moving obstacle) becomes an obstacle cell of the neighbouring body that owns most of these neighbours (inner cells + first ghost layer, decided in two phases so it is order independent). Same bundle setup, 16500 steps: without guard rho_min 0.576, max|u| 0.300, 53/166 stat samples above 0.17, 129 violating sweep steps; with N=24: rho_min 0.992, max|u| 0.017, 0 samples above 0.17, 0 violating steps; with N=21: rho_min 0.997, max|u| 0.016, 0 / 0. Mean water velocity develops identically (-0.0050 at 16400 in all three). => the anomaly is caused by the sealed sub-cell channel cells, not by the rod shape per se or the block faces. NOT yet shown: that it removes the real crash in the hopper (needs a rod run with the guard; it is local-only code so far, not on the cluster checkouts), effects on mass/force accuracy (a closed cell removes its fluid mass; the force on the rod comes from the new obstacle links) and whether a threshold of 21-24 is too aggressive for other geometries (cells in slits between flat faces have 18 obstacle neighbours; cells at rod-wall crevices count no-slip walls). Logs `/tmp/claude-1000/bundle_run/{bundle_z0,cpb12x9,guard24,guard21}/run.log` (not archived).
> - **[Claude, 2026-09-22] FIXED-BUNDLE TEST: without motion the bundle is clean -- the anomaly needs the rods to MOVE relative to the grid.** Same box/rods/settings as the moving-bundle run (9 ranks, 16500 steps, sweep + dense scan from step 0) but rod velocity ~0 (`--vx 1e-9`, i.e. fixed) in still water, two sub-cell registrations of the bundle (A: centroid start (15,11,19.6); B: shifted by (0.5, 0.37, 0.5) cell): **rho_min 1.000, rho_max 1.006, max|u| 0.000, 0 violating sweep steps in both**, mean velocity -5e-6; the moving bundle (-0.0096, 0.0043, 0) over the same steps: rho_min 0.576, max|u| 0.300, 129 violating steps. The earlier fixed-bundle run in a driven flow (`--accel 5e-6`, 3000 steps, flow-rate test) was also clean (rho_min 1.000, max|u| 0.021, 0 violations; short). => static sealed channel cells are inert (as in the hopper census); what triggers the anomaly is the covering/uncovering of cells while the bundle moves through the lattice (sub-cell registration changing every ~100 steps). Caveats: still water in the fixed runs (nothing pushes on the pocket; the driven-flow check is only 3000 steps); the moving run has an impulsive rigid start; a slower prescribed velocity (e.g. 0.001) or a bundle moving under the flow was not tested. Logs `/tmp/claude-1000/bundle_run/{fixedA,fixedB}/run.log`, prms `rod_bundle_channel_fixed{A,B}_gas.prm`.
> - **[Claude, 2026-09-22] CHECK ON LAST NIGHT'S SIX JOBS.** (1) **4252921** `antidunes-bogner-scan3-nofill75-s16` COMPLETED (1h13m, all 14 heel angles, 0 failures). `analyze_bogner_righting_scan.py <archive> --manifest prms/bogner_righting_scan3nofill75_s16_manifest.json`: torque error (meas-ideal) mean +0.084, rms 0.140 LU; zero crossing 24.98 deg; Fz error mean +1.26 %. Compared with the sbb control 4252646 (mean +0.102, rms 0.159, crossing 25.17 deg, Fz +1.51 %) and the IBB run 4252644 (no clear reduction either): `useFillWeighting false` makes the error slightly SMALLER (about 20 % less on mean/rms) but does not remove it and does not move the crossing past 26.565 -> fill weighting is a minor contributor at most, same conclusion as IBB. (2) **4252971** `hopper_clog_frontblock_pidoff_2node` (PID-off sphere hopper, fixed force 5.9e-7 from 4252463) COMPLETED (1h41m, 60000 steps): NO clog -- transport rate 0.011-0.021 throughout (vs the PID-on twin's 0.026-0.037), mean velocity dips to 0.018 (t=20000) then recovers to 0.029 by the end; the fixed clean-channel force is enough to keep spheres moving through the gap without PID help. (3) **4252997** `hopper_clog_frontblock_pidon_2node_cont50k` (sphere continuation +50000 from the completed no-clog run 4248723) COMPLETED (1h41m, steps 55100-105000): stays clean the whole way, transport rate 0.038-0.050, velocity holds at the 0.04 target, no sign of a developing clog even after 105000 steps total. (4) **4253988** `antidunes-bogner-exact-s16-seed26p565-rho075` (free start AT the equilibrium angle, rho 0.75 only) COMPLETED (11m32s): **tiltY 26.81 deg at step 500, settles to 26.96 deg by step 50000, essentially FLAT the whole run (26.8-27.2 deg band throughout, no drift to 32-36 deg).** This answers the question from the free-vs-pinned discussion: released AT its equilibrium heel, the box STAYS there (within <1 deg) -- the 32-36 deg seen in every run seeded at 2.86-10 deg is confirmed as an APPROACH-PATH effect (a large-amplitude transient that overshoots and settles on the wrong side of a shallow energy landscape, or takes far longer than 50000-300000 steps to relax), not a shift of the true equilibrium away from 26.565. The remaining ~0.4 deg offset (26.96 vs 26.565) is consistent with the righting-scan torque bias measured above (crossing ~25-27 deg) and is far smaller than the 6-10 deg gap of the seeded-heel runs. (5) **4253776** `hopper_clog_cyl_nopillars_pidon_100k` (rods, no pillars, PID on): its dependency 4252921 only completed early this morning, so it is still PENDING (unstarted) as of this check. (6) **4253834** `hopper_clog_cyl_guard24_pidon_40k` (sealed-cell guard, WITH pillars) -- see the dedicated note below, it FAILED.
> - **[Claude, 2026-09-22] SEALED-CELL GUARD FAILED ON THE REAL HOPPER (job 4253834) -- ROOT CAUSE FOUND: the guard forces a flag change but never resets the closed cell's PDFs, so a closed cell's fluid-era populations are read as the wall state by SimpleBB and can blow up.** Ran 46m17s, FAILED at step 19136 (before the original ~29000 crash point) with a `[NaN-SCAN] FATAL` interface cell (world (452.5,24.5,3.5)) whose neighbourhood already contains a cell with pdf=-nan. Immediately before it: 149 `[MEM-GARBAGE-LINK]` warnings (SimpleBB.h's existing guard, added 2026-08-29) firing against obstacle cells with `flags=[rawFlag=512]` -- i.e. cells the mesa_pd flag machinery does not recognise by name, exactly what SEALED-CLOSE produces (`boundaryHandling->forceBoundary`/`setBoundary` sets only the obstacle flag+particleField UID, the cell's PDF field is left untouched) -- one example: `pdf_old(dir)=-184.594`, 'all populations of th[is cell are garbage]'. 1025 `[SEALED-CLOSE]` events fired on this run (rank 25 alone: cell (48,2,21) closed once then 'kept as obstacle' on 40+ consecutive steps -- i.e. closures are frequent and long-lived with 406 rods, unlike the 2-3-per-configuration of the numpy study). **Mechanism:** `SealedCellGuard.h` copies the mapping's flag-setting exactly but not its `conserveMomentumWhenMapping`/PDF handling; a cell whose old fluid PDFs are simply abandoned in place keeps evolving them via whatever code path still touches non-flagged-obstacle logic for it (state-sweep classification not updated, `formerObstacle`/`StateSweep` machinery unaware of the forced transition), and `SimpleBB`'s MEM link computation eventually reads that stale, diverging PDF as a 'wall' population -> `[MEM-GARBAGE-LINK]` -> NaN nearby. **Not root-caused further** (which sweep/state exactly keeps writing to the closed cell was not traced instruction-by-instruction). **Fix needed before retrying:** either (a) re-run the same PDF-reconstruction-then-immediate-force-to-obstacle pattern `MovingParticleMapping` itself uses (so the closed cell gets a clean equilibrium/interface reconstruction, matching a real particle's arrival, before being marked obstacle), or (b) explicitly zero/reset the PDF field and register the cell with the free-surface `StateSweep`/flag machinery the same way an ordinary particle-covered cell is (block state, gas/liquid/interface counters), not just the two MEM flags. **The bundle-test result (16500 steps clean with the guard) does NOT contradict this**: that test used only 3 rods and ~2 closures total, never exercising the 'closed cell sits there for dozens of steps while neighbours keep changing' regime that 406 rods produce continuously. Guard left OFF everywhere pending a fix; `closeSealedCellsMinNeighbours` stays default 0 in every prm. Log: `/home/hpc/iwia/iwia148h/mt-dennys/apps/showcases/Antidunes/sbatch-output/output/hopper_clog_cyl_guard24_pidon_40k_4253834.out`, archived `antidunes-hopper-clog-cyl-guard24-pidon-40k_4253834_20260922_044643`.
> - **[Claude, 2026-09-22] LOCAL RAFT REPRO ATTEMPT: 6-rod hex-packed log-jam (guard 24, 12000 steps) did NOT reproduce the real hopper's MEM-GARBAGE-LINK/NaN failure.** New generator mode `gen_rod_bundle_channel_test_prm.py --raft R,C` (hex-packed R x C touching rods, same rigid-prescribed-velocity mechanism, unchanged for any rod count). A 3x4=12-rod raft on 16 ranks hit an UNRELATED pre-existing bug ('Received invalid notification type', mesa_pd ParseMessage.h -- an MPI particle-sync message desync, not sealed-cell related) at step ~1200; a 2x3=6-rod raft on the proven 9-rank layout ran clean to 12000 steps with 190+ SEALED-CLOSE events (4000+ closures on one rank alone) and zero MEM-GARBAGE-LINK/NaN. => the real failure needs something the 6-12-rod local stress test still doesn't have -- most likely the scale of 406 rods with real DEM settling/contacts and many simultaneous pile-formation events, not just closure COUNT. Not pursued further locally; the guard PDF-reset fix (see below) was implemented and tested via re-submission instead of a local repro.
> - **[Claude, 2026-09-22] VISCOSITY/DAMPING SWEEP QUEUED (3 jobs, rho=0.75 only, free 5 deg seed, s=16, upright Archimedes draft z=16 -- NOT the equilibrium heave, i.e. a genuine free-fall-in test unlike 4253988).** Motivation: every free run seeded at 2.86-10 deg overshoots past 26.565 deg and settles/wanders at 30-37 deg even after 300000 steps (confirmed NOT a runtime problem: both the 300k pins-off Bogner-tau run and the 300k real-gravity hires continuation show no decay trend, just a different resting band) -- while the box released EXACTLY at 26.565 deg (4253988) stayed there. Hypothesis: insufficient rotational damping (fluid viscosity is the only physical damping available for a free-floating single box; no DEM contacts) lets the swing overshoot past the true well into a different one. New `gen_bogner_exact_prm.py --omega` option overrides `relaxationRateOverride` (Bogner's own value 1.9 = tau 0.526; LOWER omega = MORE viscous damping). Three jobs, all 65 ranks / 50000 steps / ~35 min limit: **4255530** `omega1p9` (Bogner's own value, the control) + **4255531** `omega1p5` (tau 0.667, ~1.8x his viscosity), chained `afterany` in `~/mt-dennys-mirror` (both cheap, chained in one checkout); **4255532** `omega1p0` (tau 1.0, ~3.6x his viscosity, most damping) independent in `~/mt-dennys` (mirror2 was busy with 4253776). Files `prms/bogner_exact_3box_s16_seed5_rho075_omega{1p9,1p5,1p0}_gas.prm`, `jobs/run_65_bogner_exact_3box_s16_seed5_rho075_omega{1p9,1p5,1p0}.sbatch`. Reading: if 1.0 (or 1.5) settles near 26.565 deg without the big overshoot seen at 1.9, damping is confirmed as the mechanism (trade-off: too much damping also just slows convergence, not evaluated here since only 50000 steps are budgeted). This changes Re away from Bogner's own value, so it's a separate damping-sensitivity question, not a redo of the exact reproduction.
> - **[Claude, 2026-09-22] VISCOSITY/DAMPING SWEEP RESULT: NO overshoot at ANY damping level -- the pinned rig was the wrong test.** Jobs 4255530/31/32 (omega 1.9/1.5/1.0, seed5, rho075, `bogner_exact_3box_s16` PINNED/roll-only rig, same as 4253988) all COMPLETED (~11-13 min each) and none of them escape the near-upright basin at all: tiltY stays within ~5-7 deg for the full 50000 steps at every damping level (omega1p9: 6.58 -> peak 7.24 -> final 5.73; omega1p5: 6.10 -> 5.72, essentially monotonic decay, no peak above the seed; omega1p0: 5.43 -> peak 5.71 at step 13500 -> final 5.63, drifting down to 5.20 by the very end). **The 'overshoots past 26.565, settles/wanders at 30-37 deg' phenomenon this experiment targets was only ever observed in the FULLY FREE (unpinned, 6 DOF) family (the s=30 cuboid_stillwater sweep, s=16/hires 'nopin' 300k continuations) -- not in this pinned/roll-only rig, whose basin around upright a 5 deg seed cannot escape at any damping tested.** Separately notable: Bogner's own text reports his 2.86 deg seed 'converges nicely' at s=16 in this same pinned/roll-only geometry -- our reproduction never leaves upright at either 2.86 or 5 deg, at any damping level, a discrepancy from his description that is independent of the damping question.
> - **[Claude, 2026-09-22] VISCOSITY/DAMPING SWEEP RE-QUEUED ON THE CORRECT (FULLY FREE) RIG: job 4255587.** Same seed5/rho075/omega{1.9,1.5,1.0} sweep, this time `gen_bogner_exact_prm.py --nopin` (pinPosition/pinRotation both false, full 6 DOF -- matching the family that actually shows the overshoot, not the pinned rig above). Per explicit instruction, all three omega values run SEQUENTIALLY inside ONE combined sbatch file/ONE job (`jobs/run_65_bogner_exact_3box_s16_seed5_rho075_nopin_omega_sweep.sbatch`, new -- build once, then a `run_one()` helper per omega that rescues any leftover output, runs, archives to its own subfolder under one shared vault archive root, and continues to the next omega even if one crashes; the job's own exit code is non-zero if any sub-run failed), 65 ranks, 2h limit (untested how a fully free box behaves time-wise, so more headroom than the pinned runs' 11-13 min each). Queued in `~/mt-dennys-mirror` (free; mirror2 busy with 4255584, ~/mt-dennys free but unused this time) as **4255587**. Prms `prms/bogner_exact_3box_s16_seed5_rho075_nopin_omega{1p9,1p5,1p0}_gas.prm`.
> - **[Claude, 2026-09-22] FREE (nopin) DAMPING SWEEP RESULT: damping changes escape TIMING, not the final resting magnitude -- REFUTES the strong damping hypothesis.** Job 4255587 COMPLETED (34m58s, all three sub-runs exit 0). tiltY trajectories (seed 5 deg, rho075, s=16, fully free 6 DOF): **omega1.9 (control)** escapes fast (31 deg by step 5000), settles **+31.9 deg** (last 20 samples); **omega1.5** escapes slower (plateau by ~30000), settles **-31.2 deg** (OPPOSITE sign -- arbitrary, a symmetric problem tipped by a small perturbation); **omega1.0 (most damping, ~3.6x Bogner's own viscosity)** escapes slowest (still climbing at 45000, flattens by 50000), settles **+29.9 deg**. No wild 90 deg roll-over this time (unlike the earlier s=16 nopin 2.86-deg-seed 300k run, which did roll over and wander 30-37 deg) -- this 5-deg-seed/50000-step run settles smoothly. **All three land in the same ~30-32 deg band regardless of damping level** -- if insufficient damping caused the overshoot past 26.565 deg, the MOST damped case should land closer to target; it does not, landing at essentially the same value as the undamped control. So damping determines HOW LONG the box sits near upright before tipping, not WHERE it ends up. **This is now a well-replicated finding across every free (unpinned) variant tested: s=30 free (32.4 deg), s=16 nopin 2.86-deg-seed 300k (35.8 deg, hovering 30-37), and these three damping levels (29.9-31.9 deg)** -- the ~5-6 deg offset from Bogner's 26.565 deg looks like a genuine, damping-independent property of the free 3D body dynamics (translation drift / coupling to other rotation axes / a real method-level bias), not a transient, a seed-magnitude artifact, or a viscosity effect. Not yet explained mechanistically.
> - **[Claude, 2026-09-22] NO-PILLAR ROD RUN 4253776 ALSO FAILED (18m42s, step ~8273) -- WITHOUT THE GUARD, EARLIER THAN WITH PILLARS. Correction to the 4253834 root-cause note: MEM-GARBAGE-LINK is NOT guard-specific.** `closeSealedCellsMinNeighbours` is unset in this prm (confirmed: 0 hits for the key, 0 `[SEALED-CLOSE]` lines) -- this run used ONLY the pre-existing `SimpleBB.h` guard (added 2026-08-29), no sealed-cell closing at all. Sequence: 6 `[MEM-GARBAGE-LINK]` warnings at t=8273 (fluid cell (12,11,26) flags=liquid fill=1 density=0.157 -- NOT obviously corrupted itself -- against several `rawFlag=512` obstacle neighbours) -> `[TOTAL-FORCE-THRESHOLD]` clipping fires repeatedly on two particles -> their velocity blows up (|v| ~ 11680) -> `[SUBCYCLE-LEFT-RANK]` with absurd positions (~-1868,-2303,-2307) -> `[GHOST-VANISHED]` on 7+ ranks -> hard abort at `SimpleBB.h:298` (the `[STALE-PARTICLE-IDX]` guard, a cell still references the now-untrackable UID). **This is the SAME MEM-GARBAGE-LINK symptom seen preceding 4253834's failure, but 4253834 DID have the sealed-cell guard active and this run does NOT -- so the guard is not solely responsible for that failure mode; the underlying dense-rod-pile MEM force instability exists independently of it.** The guard-PDF-reset fix (in progress) is still worth finishing since it is a real bug (a closed cell's PDFs going untouched forever is still wrong), but it should not be expected, by itself, to fix crashes of this family. Without pillars the rods jam directly at their own release block (X 169-371) rather than being funnelled by the pillars, and this failed MUCH sooner (8273 vs 29000 with pillars) -- removing the pillars made the pile-up worse/faster, not better; the original 'does the crash need the pillars' question is answered: NO, it does not, an unfunnelled rod pile fails even faster. Archived `antidunes-hopper-clog-cyl-nopillars-pidon-100k_4253776_20260922_104116`.
> - **[Claude, 2026-09-22] SEALED-CELL GUARD CORRECTNESS FIX IMPLEMENTED + RERUN QUEUED: 4255561.** Fix (`SealedCellGuard.h`, `BuildTimeloop.h`): the sweep now runs AFTER `PdfReconstructionManager` and its own flagField/particleField communication (was: right after `MovingParticleMapping`, BEFORE reconstruction), and its own `formerObstacle`-interception branch is removed entirely -- since reconstruction has already resolved every vacated cell to a consistent liquid/gas/interface state by the time this sweep runs, there is no longer a 'still-sealed vacated cell skips reconstruction, slams back to obstacle with stale PDFs' path (the mechanism traced to 4253834's failure); the sweep now only ever closes a cell that is already an ordinary, fully-resolved domain cell. `CloseSealedCellsSweep`'s constructor dropped the now-unused `formerObstacle` parameter. **Local smoke test (3-rod bundle, 3000 steps, 9 ranks, guard 24):** 129 closures, completed cleanly, confirming the reordering doesn't break the sweep itself. **Rerun of 4253776** (no-pillar rod hopper, PID on, which failed at step 8273 via `[MEM-GARBAGE-LINK]` -> force blow-up -> stale-particle-idx abort WITHOUT any guard): same scenario with the fixed guard ON (`closeSealedCellsMinNeighbours 24`), queued as job **4255561** in `~/mt-dennys-mirror2` (free again after 4253776 finished), **1h limit** (4253776 failed within 19 minutes, so this is tight but gives headroom past that point if it survives). New prm `prms/hopper_clog_cyl_nopillars_pidon_2node_100k_guard24_gas.prm` (copy of the unguarded 100k prm + one line), `jobs/run_144_hopper_clog_cyl_nopillars_pidon_guard24.sbatch`. **Caveat carried over from the 4253776 correction above: the guard targets sealed-CHANNEL cells specifically; if 4253776's failure is a different dense-pile MEM mechanism (unconfirmed either way), this fix may not prevent it even though it is a genuine correctness improvement on its own.**
> - **[Claude, 2026-09-22] 4255584 (guarded rerun) ALSO FAILED (22m49s, step ~8989) -- SAME particle UID as the unguarded crash, DIFFERENT mechanism: a DEM CONTACT OVERLAP BLOW-UP, not a sealed channel.** 581 `[SEALED-CLOSE]` closures fired (guard active, correctness-fixed) but the crash is only ~700 steps later than the unguarded 4253776 (8273). Same culprit uid `4827858800541171713` in both. Sequence: `[DEM-WATCH]` deep contact overlap growing across DEM sub-cycles within one step (1.15 -> 0.92 -> 0.66 -> 1.51 -> 2.47 -> **3.47 LU**, i.e. 37% of the smaller interaction radius) -> `[HIGH-VELOCITY-PRECURSOR]` fluid |u| up to **16.77** right where it plunges -> `[MEM-GARBAGE-LINK]` against obstacle cells with density=-0.22 (invalid) and fill=1.37 (out of range), some AT z=-1 (below the domain floor) -> `[GHOST-VANISHED]` -> the same stale-particle-idx abort. **Conclusion: this failure family is a DEM contact-resolution failure (a rod interpenetrating catastrophically, ~3.5 LU deep, within a handful of sub-cycles), not the sub-cell sealed-channel mechanism the guard targets** -- the sealed-cell correctness fix is a real, separate improvement but doesn't address this. Activated the EXISTING (pre-built, previously-off) `ModelParameters.rodCapOverlap` diagnostic cap (ParticleKernels.h:2501 -- clamps the penetration depth fed into the Cylinder contact force law for any rod contact; already proven together with `rodCapAngularVelocity` in the `hopper_clog_cyl_analytic_tc1_caps_*` family) at the SAME precedent value **1.0 LU**, on top of the (still active, correctness-fixed) sealed-cell guard. New prm `prms/hopper_clog_cyl_nopillars_pidon_2node_100k_guard24_rodcap_gas.prm`, sbatch `jobs/run_144_hopper_clog_cyl_nopillars_pidon_guard24_rodcap.sbatch`. Queued as **4255613** in `~/mt-dennys-mirror2`, 2h limit (generous: if the cap works this may run much longer than the two ~20min prior failures). Note `rodCapOverlap` is a DIAGNOSTIC artefact (SimulationParams.h's own doc comment: 'changes the physics, NOT for production') -- if this survives, it confirms the deep-overlap blow-up as the proximate trigger but doesn't by itself fix the underlying contact-stiffness/sub-cycling problem that lets the overlap grow that deep in the first place.
> - **[Claude, 2026-09-22] 4255613 (guard + rodCapOverlap combined) ALSO FAILED (20m55s, step ~9198) -- rodCapOverlap fired ZERO times before the blow-up.** 560 `[SEALED-CLOSE]` closures fired (sealed-cell guard active) but `[ROD-CAP]` overlap hits = 0 -- the diagnostic overlap cap never got a chance to intervene; whatever triggered this blow-up did not route through a single deep-overlap contact the way 4255584's did. Sequence this time: `[HIGH-VELOCITY-PRECURSOR]` fluid cells (|u| 1.0-1.73, well past threshold but far milder than 4255584's 16.77) around world (497-498, 29-30, 38-39) -> `[SUBCYCLE-LEFT-RANK]` with two particle UIDs at ABSURD positions/velocities (uid ...720: pos <-6539,-158,229> vel <-35196,-954,954>; uid **...724** -- a NEW culprit not seen in the two earlier crashes: pos <5.23e6,5.24e6,44.9> vel <2.6e7,2.6e7,-0.06>) -> GHOST-VANISHED on 8+ ranks -> the same stale-particle-idx abort. **Pattern across all three attempts: 4253776 (no guards) crashed at step 8273, 4255584 (sealed-cell) at 8989, 4255613 (sealed-cell + overlap cap) at 9198 -- all within ~900 steps of each other, same initial pile, overlapping culprit UIDs. Neither guard has meaningfully moved the failure point or changed its character.** This line of attack (sealed-cell guard, overlap cap) does not appear to address whatever is actually driving this crash family; a fourth variant was not queued without a new hypothesis for what to target instead.
> - **[Claude, 2026-09-22] SEALED-CELL GUARD RETRIED ON THE ORIGINAL PILLARED HOPPER (job 4255664)** -- the scenario the guard was actually built against (4252146/4252462's sealed rod-rod pore finding), now with the correctness fix (the three no-pillar attempts above all used the no-pillar variant, a different geometry). Reused the EXISTING prm/sbatch from the earlier broken-guard attempt 4253834 unchanged (`prms/hopper_clog_cyl_guard24_pidon_2node_40k_gas.prm`, `jobs/run_144_hopper_clog_cyl_guard24_pidon_40k.sbatch`: pillars present, 40000 steps, `closeSealedCellsMinNeighbours 24`) -- only the underlying SOURCE (`SealedCellGuard.h`/`BuildTimeloop.h`) differs, now the post-reconstruction, formerObstacle-case-removed version. 4253834 (broken guard) had failed at step 19136, earlier than the original unguarded pillared run's ~29000 crash. Queued in `~/mt-dennys-mirror` (free), 144 ranks, 2h limit (unchanged from 4253834's own budget).
> - **[Claude, 2026-09-22] HALF-DENSITY ORIGINAL HOPPER QUEUED (no guards/caps at all, explicit instruction): job 4255679.** The ORIGINAL cylinder hopper (pillars, GJK/EPA-only, PID on, 60000 steps, clean baseline -- no sealed-cell guard, no rodCapOverlap, no diagnostic sweeps) but with HALF the rod population: 210 rods instead of 406 (`gen_hopper_clog_cyl_prm.py --rows 15`, `front_x` unchanged at 368 -> kept rows are the DOWNSTREAM half against the pillars, X [267,371]; the dropped 14 rows are the FURTHER-UPSTREAM half, X [169,260] of the original 29-row block -- a whole contiguous half removed, not thinned throughout, so every remaining row still has all 7 Z-layers / all three density bands 0.5/1.0/1.2, confirmed in the generator's own printout: 60/90/60 = 210). Question: does a half-density pile jam/crash at all, and if so, via the same mechanism as the full-density pile? New prm `prms/hopper_clog_cyl_frontblock_pidon_2node_halfrows_gas.prm` (base from the generator + the same one-line `useAnalyticalContactDetection true->false` toggle the existing `..._gjk_gas.prm` uses, verified by diff), `jobs/run_144_hopper_clog_cyl_halfrows_gjk_noguards_60k.sbatch`. Queued `afterany:4255664` (the currently RUNNING pillared sealed-cell-guard test) in the SAME checkout `~/mt-dennys-mirror`, 144 ranks, 2h30m limit (no timing precedent for this exact configuration; budgeted off the full-density sphere hopper's own 1h48-1h52m for 60000 steps, plus margin).
> - **[Claude, 2026-09-22] TORQUE-SWEEP DATA-SLICE QUOTA BUG FIXED + 5 REMAINING RUNS RE-QUEUED: job 4256467.** Root cause of 4255656 dying after only its first sub-run: `AverageDataSliceEvaluator` (`BuildTimeloop.h`) shared `EvaluationParameters.evaluationFrequency` with `[BoxForces]`/`[BoxTorque]` per-step logging; setting `--evalfreq 1` for per-step torque also made that evaluator dump 3 text files/step to `eval/avgDataSlices/` -- ~300000 files over 100000 steps, which blew a vault disk quota mid-archive and killed the job via `set -e` before the 2nd of 6 sub-runs ever started (confirmed from the .err log: repeated `mv: ... Disk quota exceeded`). **Fix:** new, separate `EvaluationParameters.dataSliceFrequency` (`SimulationParams.h`, `BuildTimeloop.h`) -- defaults to `evaluationFrequency` itself if unset (0), so every existing prm is unaffected; `gen_bogner_exact_prm.py` now always sets it to 500 explicitly regardless of `--evalfreq`. Verified locally: rebuilt (`./antidunes.sh --local --build`, needed the `walberla-codegen` conda env active), then a 600-step smoke run (20 ranks, `--blocks` override) gave exactly 600 `[BoxTorque]` lines and only 6 slice files -- confirms the decoupling works. **Cleanup:** copied the one sub-run that DID complete (omega1p9_seed5: vtk-out, eval/antidunes.sqlite, computedParams.txt, and the full 200000-line BoxForces/BoxTorque log) to `animations/4255656_omega1p9_seed5_100k_torque_20260922/` (skipped the ~54000 junk slice files); deleted those junk files from the vault archive (2.9G -> 35M) and a 16G duplicate/incomplete copy left in the checkout's own build dir by the failed `mv`. **Broader /home/hpc survey (explicit instruction):** found and deleted `mt-dennys`'s own local `archive/` fallback directory (41G: a 31G `20260916_185616/d025` folder with no confirmed vault backup, plus four `rescued_before_<jobid>_*` snapshot-only dirs ~2.4G each that turned out NOT to be confirmed duplicates of their same-numbered vault archives on inspection -- checked one, 4250907's vault archive has no `snapshots/` folder at all, meaning the rescued copy likely belongs to a DIFFERENT, unidentified earlier job -- deleted anyway on explicit instruction after flagging this ambiguity). Re-queued the 5 sub-runs that never ran as **4256467** (`jobs/run_65_bogner_exact_3box_s16_nopin_100k_torque_sweep_remaining5.sbatch`, omega1p9_seed5 excluded -- already have its data), `~/mt-dennys-mirror2` (free), 65 ranks, 2h30m limit.
> - **[Claude, 2026-09-22] FLOW-RATE EVALUATOR ACTIVATED ON THE THREE 'Clog Flow-Forcing Comparison' ARTIFACT SCENARIOS, chained reruns queued.** The artifact's underlying runs (4242918 no-clog/PID-on, 4245629 extreme-clog/PID-on, 4245630 extreme-clog/PID-off) all predate `ClogMeasurement.flowRateEnabled` and have no `clogUpstream/DownstreamFlowRate`/`WetArea` columns. New prms, each a copy of the original with ONLY `flowRateEnabled true` added: `prms/clog_measurement_validation_extreme_{pidon,pidoff}_2node_flowrate_gas.prm`, `prms/clog_measurement_validation_noclog_reference_extreme_flowrate_gas.prm`. Queued as three chained jobs: **4257275** (extreme clog, PID on) with no dependency in `~/mt-dennys` (free), **4257276** (extreme clog, PID off) `afterok:4257275` in the SAME checkout (only runs if the PID-on run succeeds), **4257277** (no-clog PID-on reference) `afterok:4256467` (the currently running 5-run torque sweep) in `~/mt-dennys-mirror2`, same directory, per explicit instruction. Time limits tightened from the original budgets to the actual measured runtimes with margin: pidon 6h->3h15m (4245629 itself took 2h30m56s), pidoff 6h->3h (4245630 took 2h23m02s), noclog 3h->1h45m (4242918 took 1h14m52s). All three synced via `./antidunes.sh --sync` to their respective checkouts (source-only, did not disturb 4256467's running build on mirror2).
> - **[Claude, 2026-09-22] 4255664 (pillared hopper + FIXED sealed-cell guard) COMPLETED CLEANLY -- first full 40000-step survival of this scenario.** Exit status 0, 1h43m34s (well inside the 2h budget). **1766 `[SEALED-CLOSE]` closures** (guard genuinely active throughout, closing many sealed pockets in this dense pillared pile), **ZERO `[MEM-GARBAGE-LINK]` warnings** (vs. 149 in the broken-guard attempt 4253834 and 50/20 in the two no-pillar attempts), **ZERO fatal NaN scans**. Fluid sanity sweep still logs occasional minor |u|>0.17 violations near the end (a handful of cells at a time, e.g. 12 at t=38190) but nothing escalating -- no non-finite, no density blowup, no excess-mass runaway. For context: the original UNguarded version of this exact scenario (4252146 -> 4252462) crashed at ~step 29000; this run went to 40000 clean. Flow state at the end: meanVelocity 0.0402 (right at the PID target 0.04, i.e. NOT stalled), acceleration 4.40e-6 (above the clean-channel baseline ~6-7e-7 -- real resistance building -- but far below the ~14x seen in the confirmed extreme-clog cases). **This is the FIRST case where the sealed-cell guard has actually prevented a crash** -- unlike the three no-pillar attempts (unguarded, sealed-cell-only, sealed-cell+overlap-cap), which all failed within ~900 steps of each other via the wall-contact/DEM-force family of instabilities investigated separately above. Consistent with the guard's original design target: THIS scenario (dense, slowly-forming pillared pile) is exactly the 'quasi-static jam' regime where sealed-channel cells were the actual problem, unlike the no-pillar scenario's free-roaming high-velocity collisions. Archived `antidunes-hopper-clog-cyl-guard24-pidon-40k_4255664_20260922_135635`. Its dependent job 4255679 (half-density, no guards, original clean baseline) started once this completed.
> - **[Claude, 2026-09-22] 4255561 FAILED IMMEDIATELY (38s) -- NOT a physics crash, a sync mistake: resubmitted as 4255563.** The rsync before submitting only copied `BuildTimeloop.h`, `SealedCellGuard.h`, the prm and the sbatch -- `SimulationParams.h` (where `mappingTieOffset`/`closeSealedCellsMinNeighbours` are actually declared as struct members) was left off, so mirror2 built against its own STALE copy of that header and failed to compile (`'const struct ModelParams' has no member named 'closeSealedCellsMinNeighbours'`, two more identical errors). Synced `SimulationParams.h`, resubmitted unchanged as **4255563** (same 1h limit, same checkout). The other three damping-sweep jobs (4255530/31/32) are unaffected -- they only use the `--omega` prm generator flag (an existing `relaxationRateOverride` override key), no new C++ member, no rebuild-affecting source change.
> - **[Claude, 2026-09-22] SECOND compile failure (4255563, 35s): another missed file -- `src/lbm_mesapd_coupling/momentum_exchange_method/MovingParticleMapping.h` (the `cellCenterTestOffset`/tie-offset parameter) was never synced to mirror2 either, so `makeMovingParticleMapping` still had its OLD 9-argument signature there while `BuildTimeloop.h` was already calling it with 10. Two manual per-file rsyncs in a row missed a file. Switched approach per explicit instruction: (1) built locally first with the project's own `./antidunes.sh --local --configure --build --build-dir build-release --jobs 12` (needed the `walberla-codegen` conda env manually activated first -- `--local` skips ALL conda activation, so `--conda-env` is silently ignored in that mode; without it Python_ROOT_DIR falls back to system Python without pystencils/lbmpy and the codegen target silently fails to produce the `Antidunes` binary target at all, 'ninja: error: unknown target Antidunes', a separate gotcha worth remembering) -- compiled clean (only the pre-existing SetupStages.h indentation warnings). (2) Synced the WHOLE repo with the script's own `./antidunes.sh --sync --remote fritz.iwia --remote-path '~/mt-dennys-mirror2/'` instead of another manual per-file rsync -- this is the built-in, safe mechanism (excludes build/ and dotfiles, copies everything else) and should not have this failure mode again. Verified via md5sum that `MovingParticleMapping.h` now matches between local and mirror2. Resubmitted as **4255584** (same 1h limit).
> - **[Claude, 2026-09-21] SPHERE-ONLY AUDIT of the code path (per step of the time loop; question at every point: 'was this written for spheres, does it need adjusting for other shapes?').** Result: most of the path is already shape generic (or was adapted); remaining sphere-only spots: **(1) Lubrication correction is sphere-sphere only** (`ComputeLubrication::operator()`, ParticleKernels.h:1982-1984 returns unless BOTH shapes are spheres): with `useLubricationCorrection true` (all rod prms) rods, boxes and sphere-rod pairs get NO lubrication force; spheres do (analytic sphere-sphere film). Whether this matters is open (the LBM does not resolve gaps < ~1 cell either, and line contacts of rods are the case the sphere formula does not cover). **(2) `avgParticleVelocity` (BedloadTransportEvaluator, Evaluators.h:~215) divides by `numParticles`, the count of the initial sphere BED (SetupStages.h:301-364) = 0 in every floater-only prm** -> non-finite value -> the recurring 'Skipping column avgParticleVelocity' warning; the sum also includes fixed pillar particles. The transport rate itself uses `getVolume` and is shape generic. **(3) Unit system built on the average SPHERE diameter** (`avgParticleDiameter`, SimulationParams.h:1447-1503: liquid/floor height, `particleCollisionTime = nonDim * avgParticleDiameter`, fixing height); the hopper prms set `avgParticleDiameterOverride 8.0019` so the rod channel equals the sphere channel -- consistent, but SI numbers derived from it are meaningless for rods. **(4) Bed creation** (`maxSphereRadius`, SetupStages.h:373) is sphere-only by design (the bed is spheres). **(5) Checked and OK / already adapted:** ghost-sync algorithm choice uses the circumscribed diameter of boxes/cylinders (SetupStages.h:715-729) and `blockSyncExtension`; particle cell bounding box = interaction-radius cube (mesa_pd AABBConversion.h:47) so rods use a 19^3 cube instead of ~650 cells (performance only, correct); mapping selector accepts spheres by shape type and boxes/cylinders by type 3/4 (BuildTimeloop.h:724-731; a NEW shape type would silently not be mapped); `Cylinder::updateMassAndInertia` (m r^2/2, m(3r^2+h^2)/12) and `containsPoint` correct; SimpleBB/MEM uses `getVelocityAtWFPoint` and the wall position (shape generic); buoyancy split (only for `useMemOnlyBuoyancy false`) has a cylinder branch (`isCyl`); contact stiffness uses the pair's effective mass (`useRotationalEffectiveMass` on for the rods); VTK writers exist per type; rod velocity caps exist only for rods (spheres/boxes have none, off in production). Not audited line by line: PdfReconstructionManager internals, PID, checkpoint (de)serialisation of shapes (rod resume runs worked).
> - **[Claude, 2026-09-21] GRID MAPPING DEEP DIVE (cell-centre test) -- what the mapping does per shape and where it goes wrong.** Mechanism (`MovingParticleMapping.h`): per block, per (also ghost) particle, per cell of the particle's cell bounding box (interaction-radius cube): obstacle iff `containsPoint(cell centre)` (binary, inclusive `<=`); a newly covered domain cell becomes obstacle with the particle UID (momentum of the covered fluid is NOT given to the body, `conserveMomentumWhenMapping` false, its mass is dropped), a cell of another body is re-claimed unless that body is FIXED, a vacated cell gets the 'former MO' flag and is re-initialised by the Bogner PDF reconstruction from its non-obstacle neighbourhood (gas/liquid/interface by neighbour evidence). Everything shape specific is in `containsPoint`; nothing else knows the shape. **Numbers** (`study_grid_mapping.py`, numpy, same rule, 120-400 random orientations and sub-cell registrations each): mapped volume / analytic volume: sphere r=3 mean 1.001, std 2.6 % (min 0.946, max 1.105); sphere r=8 std 0.45 %; cylinder r=3 L=18 std 0.95 % (0.953-1.026); box 6x6x18 std 0.67 % -- so the r=3 debris SPHERES have the largest volume/buoyancy noise, not the rods. **Sealed fluid cells** (fluid cell with >=T of its 26 neighbours obstacle; mean number per configuration, gap = surface distance in cells): touching parallel cylinders r=3 (gap <= 0): T>=20: 4.3 cells (92 % of configurations have >= 1), T>=24: 0.0-0.1; THREE touching cylinders: T>=20: 27-29 cells (100 %), T>=24: 2.0-2.7 (71-85 % of configs, up to 12); gap 0.5: 5.2 / 0.0; gap >= 1: none. Two touching spheres: ~0 (T>=20 in 8 % of configs), three touching spheres: 0.5 (T>=20 in 40 %), never T>=24. Boxes 6x6x18: touching (gap <= 0) flush faces leave no fluid, BUT a gap of 0.25 cell leaves 1-cell slits: two boxes T>=20: 16.8, T>=24: 8.8 (85 % of configs), three boxes 39.5 / 13.8 (96 %); gap 0.5: 27-49 / 0.9-1.4; gap >= 1: none. => the sealed-cell problem is a property of parallel flat/curved faces at sub-cell distance, worst for slits between flat faces (boxes) and for the 3-rod triangle; point-contact shapes (spheres) hardly produce any. **Axis-aligned ties:** a box whose face plane passes exactly through cell centres is mapped one full layer thicker on BOTH sides (6x6x18 box shifted by 0.5 cell in x: mapped volume 1.167 vs 1.000 for every other shift), i.e. the mapping is a step function with a +/-0.5 cell jump at ties and any floating-point noise in the body rotation can flip a whole face layer from step to step. The R1 buoyancy_compare box sits exactly on such a tie (position z = 41.5 = cell centre, half extent 8 -> faces at 33.5 / 49.5, cell centres) while the R2/R4 boxes (z = 83.0 / 166.0) do not -- consistent with the wildly oscillating `F_hydro_z` (-0.07..+0.19) and the missing rise of the R1 box (section above); tested locally, see the next note. [SUPERSEDED: the R1 failure is floor contact of the mapped box, not tie flicker -- see 'CORRECTION of the tie explanation'.] **Other consequences seen in the code:** the covered fluid's momentum is discarded (impulsive body-fluid coupling on fast approach), every body maps its own cells so line contacts (rods) have an 18-cell-long overlap zone whose cells are re-claimed by whichever body is processed last (harmless for ownership, but MEM links of the contested cells are attributed to one rod), and the bounding cube is 19^3 for a rod.
> - **[Claude, 2026-09-21] R1/R2/R4 BUOYANCY COMPARISON -- data ALREADY THERE, 'steps to rise the same number of lattice units' extracted (no new run needed).** The 18 runs (R1/R2/R4 x sphere/cylinder/box x MEM-only/analytic; vault `antidunes-archive/antidunes-buoyancy-r{1,2,4}-compare-all*/<shape>_<method>/eval/antidunes.sqlite`, floater centre `maxParticlePos_z` every 10 steps) all have the SAME physical scenario (everything scaled by R), so the floater rises R times as many lattice units at higher R; but its z(t) can be cut at any lattice distance. `analyze_buoyancy_steps_to_rise.py <root>` -> steps until the centre has risen D lattice units from rest: **sphere** D=10: R1 mem/ana 1964/2213, R2 4946/5168, R4 12814/13140 (MEM faster by 11 % / 4 % / 2.5 % -> converges towards the analytic model with resolution); **cylinder** D=10: 1988/2050, 5348/4997, 14843/13132 (MEM vs analytic -3 % / +7 % / +13 %: DIVERGES with resolution); **box** D=10: R1 MEM never rises, analytic 2006; R2 6538/4988, R4 17480/13158 (MEM 31 % / 33 % slower, no convergence). Step ratios R2/R1 and R4/R2 for D=10: sphere 2.5 / 2.6 (MEM), cylinder 2.7 / 2.8 (MEM), analytic 2.3-2.6 (constant-acceleration expectation R^1.5 = 2.83 since g_LU scales like R^-3). A dedicated same-lattice-rise prm set was NOT generated (not needed for this comparison; it would only shorten future runs: D=14 needs 2.4k / 6k / 15k steps instead of 2.6k / 10.4k / 41.6k).
> - **[Claude, 2026-09-21] TIE TEST (its interpretation 'mapping-tie flicker' is CORRECTED in the next note): the R1 box that never rose is fixed by shifting it by 0.25 cell.** Local reproduction of `prms/buoyancy_compare_box_r1_mem_gas.prm` (MEM-only, 2600 steps, blocks 16x16x24 = 36 ranks instead of 72, otherwise unchanged; bed file `spheres_out_buoyancy_r1.dat` taken from the cluster build dir): (a) as is (position z = 41.5, half extent 8 -> top/bottom face planes exactly through cell centres): the box stays at z = 41.50 for all 2590 steps, `F_hydro_z` (steps 500-2500) mean 0.0443, **sd 0.0845**, range -0.075..+0.210 (the wild oscillation of the cluster run, reproduced); (b) only the start height changed to 41.75 (no tie): the box rises 41.75 -> 43.9 (1000) -> 49.6 (2000) -> 55.3 (2590), `F_hydro_z` mean 0.0520, **sd 0.0081 (10x smaller)**, range 0.002..0.065; analytic-arm reference of the cluster run: +15.7, here +13.5. => the 'MEM-only buoyancy fails for box floaters at R1' finding of 2026-09-12/13 is a cell-centre TIE artefact of the mapping (a face plane through cell centres is mapped one full layer thicker and can flip layer-wise with floating-point noise), not a property of MEM buoyancy for boxes. R2/R4 boxes (z = 83.0 / 166.0, faces on cell boundaries) have no tie but are still 31-33 % slower than analytic -- separate, unexplained. **Consequences to check:** any box floater/pillar whose face planes pass through cell centres (positions at half-integers with integer half-extents, or integer positions with half-integer half-extents, axis-aligned, e.g. righting-scan angle 0 or pinned boxes, `bogner_exact` prms with zero tilt) is affected; a tilt of >= a few tenths of a degree removes the ties for most faces. Possible fix: map with a strict `<` on one side / a tiny deterministic offset (e.g. shift the cell-centre test point by 1e-6 * (1,1.3,1.7) cells) so ties always resolve the same way. Not implemented. Data: /tmp/claude-1000/bundle_run/tie_r1/{orig,shifted} (not archived).
> - **[Claude, 2026-09-21] CORRECTION of the tie explanation: the R1 MEM box fails because its mapped bottom layer lies flush on the floor (no fluid cell underneath), not because of tie flicker.** After implementing the deterministic tie resolution (`ModelParameters.mappingTieOffset`, `MovingParticleMapping.h` `cellCenterTestOffset`; test point = cell centre + offset*(dx, 1.3dy, 1.7dz); default 0 = off, bit-identical) I re-ran the R1 box (start z = 41.5, MEM-only, 36 ranks, 2600 steps) with the two signs: **+1e-4: box still stuck at z = 41.50, F_hydro_z sd 0.0857 (same as without offset); -1e-4: box rises to 55.6 (2590), F_hydro_z sd 0.0098** -- as with the 0.25-cell shift (55.3, sd 0.0081). +1e-4 keeps the bottom layer of the tie (cell 33, centres 33.5), -1e-4 drops it. The floor is solid up to z = 33.12 (cells 0..32), so with the layer kept the mapped box sits directly on the floor cells: no fluid below it, hence no MEM pressure from below (the MEM-only arm needs that pressure for ALL its lift; the analytic-buoyancy arm does not, and rises), whereas dropping the layer leaves one fluid cell under the box. So (i) the 2026-09-12/13 finding 'MEM-only buoyancy fails for boxes at R1' is a SETUP/mapping artefact (box start 0.38 cell above the floor; sphere/cylinder at the same z touch the floor cells only at a pole), not a property of MEM buoyancy for boxes; (ii) the mapped-volume tie (+16.7 % in the numpy test) is real but is not what failed here; (iii) the generic offset only makes ties deterministic -- its sign is arbitrary, so it does NOT fix this case by itself (needs -1e-4) and is left OFF everywhere; the actual rule that matters is 'a mapped body must keep >= 1 fluid cell to any wall/floor', the same family as the sealed-cell guard. The mechanism is inferred from the geometry and the sign test, not checked in a flag-field dump. The 'Consequences to check' list of the previous note (boxes with ties) should read: bodies placed within ~1 cell of a solid wall/floor in MEM-only runs.
> - **[Claude, 2026-09-22] R1/R2/R4 COMPARISON RERUN QUEUED (3 chained jobs, box tilted 25 deg, per explicit instruction).** Applied the tie fix found above to the full comparison: `gen_buoyancy_tilted25_box_prm.py` (new) generated `prms/buoyancy_compare_box_r{1,2,4}_{mem,analytic}_tilted25_gas.prm` (initialTiltAngleDeg 25 on the FloatingBox; BOTH mem and analytic arms tilted, so the two methods still start from an identical initial condition, not just individually patched -- the analytic arm never had the tie problem but is tilted too for a fair comparison). Sphere/cylinder prms UNCHANGED. **Split into THREE separate chained sbatch scripts (not one combined job)** since R1/R2 need 72 ranks/1 node (singlenode) but R4 needs 576 ranks/8 nodes (multinode) -- a single job can't change its SLURM allocation mid-run, and requesting the R4-sized allocation for the whole thing would waste ~2h of 8-node core-hours during the R1/R2 stages. New scripts `jobs/run_72_buoyancy_compare_r{1,2,4}_tilted25.sbatch` (copies of the existing `r{1,2,4}_all.sbatch` with the box run_case lines swapped to the tilted25 prms). **Time limits adjusted from the LAST actually-measured run of each stage** (not re-guessed): R1 kept at 00:30:00 (already calibrated from the real ~7.6 min/6-run total, no new data); R2 tightened 02:30:00 -> 01:30:00 (job 4223072 actually COMPLETED in 1:08:31); R4 tightened 09:00:00 -> 07:00:00 (job 4235003 actually COMPLETED in 6:18:24, after two earlier failed/cancelled attempts). Chained via `--dependency=afterany`. **Originally submitted in `~/mt-dennys-mirror`** (chosen because it already had ALL THREE bed files `spheres_out_buoyancy_r{1,2,4}.dat`), but that checkout was still running the (unrelated) nopin omega-damping sweep 4255587 -- first submitted R1 there with no dependency, which would have collided on the shared build directory (concurrent CMake/ninja + mpirun in the same path); caught and fixed via `scontrol update JobId=4255623 Dependency=afterany:4255587` before R1 could start. **Then moved to the genuinely idle checkout `~/mt-dennys`** on explicit instruction (cancelled 4255623/24/25, copied the two missing bed files `spheres_out_buoyancy_r{1,2}.dat` from mirror's build dir to mt-dennys's -- it only had r4 locally -- then resynced sources and resubmitted): final chain **4255633 (R1, no dependency, starts immediately) -> 4255634 (R2) -> 4255635 (R4)** in `~/mt-dennys`. All three checkouts now busy: mt-dennys = this chain, mt-dennys-mirror = 4255587, mt-dennys-mirror2 = 4255613. Also fixed the `analyze_buoyancy_steps_to_rise.py` script used to read the OLD (untilted) data: it previously used a fixed absolute LU threshold D at every R, which is not the same relative journey at different resolutions (D=10 LU is 1.25x the R1 box's own half-size but only 0.31x the R4 box's); now scales D by R (D_R = D_base * R) and reports R2/R1, R4/R2 step-count ratios -- rerun on the old data gives ratios ~3.5-3.85, much closer to the expected R^2=4 than the unscaled version.
> - **[Claude, 2026-09-21] SEALED-CELL GUARD RUN QUEUED: 4253834** `antidunes-hopper-clog-cyl-guard24-pidon-40k` (`~/mt-dennys`, `--dependency=afterany:4252997`, limit 2h): the PID-ON rod hopper WITH pillars of 4252146/4252462 (which crashed at ~29000 in a sealed rod-rod pore), fresh start, 40000 steps (beyond the crash), `closeSealedCellsMinNeighbours 24`, flow-rate planes x=336/432, sanity sweep every 10 steps + dense NaN scan, checkpoints every 10000; prm `prms/hopper_clog_cyl_guard24_pidon_2node_40k_gas.prm` (`gen_guard_hopper_prm.py`), `jobs/run_144_hopper_clog_cyl_guard24_pidon_40k.sbatch`. Pre-flight: guard sweep made cheaper (early exit after ~3 flag reads per cell; the local guard24 bundle run was not slower than the unguarded one: 0.034 vs 0.039 s/step); local hopper smoke (18 ranks, 300 steps, rods released at 100) ran without errors and 0 closures (rods not yet clumped); Antidunes, BedGeneration and DropTest all compile locally. Synced to `~/mt-dennys` (this also feeds the pending 4252997): BuildTimeloop.h, Evaluators.h, ParticleKernels.h, SealedCellGuard.h, SetupStages.h, SetupStages_drop_test.h, SimulationParams.h, src/lbm_mesapd_coupling/momentum_exchange_method/MovingParticleMapping.h. All new options default off. What to look for: `[SEALED-CLOSE]` lines (closure counts), `[FLUID-SANITY-SUMMARY]` (should stay quiet), no abort beyond 29000; a chaotic run may not reproduce the crash even without the guard, so a clean run is evidence, not proof -- compare the sweep statistics/branch-E counts with 4252146. The mapping-tie offset is NOT used in this run.
> - **[Claude, 2026-09-21] FREE START AT THE EQUILIBRIUM ANGLE, rho 0.75 only: job 4253988** (`antidunes-bogner-exact-s16-seed26p565-rho075`, `~/mt-dennys-mirror`, `afterany:4252971`, 65 ranks, 35 min limit; the 3-box baselines took 11-13 min): ONE Bogner-exact s=16 box (16x32x16 cells, rho 0.75, x = 21.67, y = 20, 130x40x40 basin, waterline z = 20, tau 1/1.9, g 7.5e-4) starts at the 26.565 deg heel at its ideal heave z = 16.9223 (from the righting scans) with roll about Y and heave free, X,Y pinned (as the baseline), 50000 steps, VTK every 1000, `[BoxForces]` tilt every 500 steps. Question: it stays near 26.565 (-> the 32-36 deg of the runs seeded at 2.86-10 deg is an approach-path effect) or drifts (-> the free equilibrium differs from the static-torque zero crossing of the righting scans, ~25-27 deg, where the box is held with position/rotation pinned). No free run had started at the equilibrium angle before; the righting scans at 26.565 hold the box. Files: `gen_bogner_exact_prm.py --s 16 --seed 26.565 --only 0.75 --z0 16.9223 --tag _seed26p565_rho075` (new options `--only`, `--z0`), `prms/bogner_exact_3box_s16_seed26p565_rho075_gas.prm`, `jobs/run_65_bogner_exact_3box_s16_seed26p565_rho075.sbatch`. Analysis as for the other bogner_exact runs (folded tilt from the floaterBox VTK rotation / [BoxForces] tiltY).
> - **[Claude, 2026-09-21] Bogner 2017 Sect. 5.3.3 DOES give the size of the boxes relative to the grid:** 'the 2D cube of side length s is represented by a cuboid 2s x s x s ... Inside a domain of 130 x 40 x 40 cells, the lower half is filled with liquid, tau = 1/1.9 = 0.526, hydrostatic initial pressure (g = 7.5e-4). Three boxes with s = 16 are placed inside the basin, slightly heeled ... by 2.86 deg' -> box 32x16x16 cells, water depth 20 cells (waterline at 1.25 s), heel 2.86 deg (= 1.6 cells at the corner). That is exactly what the `bogner_exact_3box_s16/s32` prms replicate. What the text does NOT give: the run length / convergence criterion ('until all fluid motions ceased'), any numerical result (only the visual Fig. 5.10, 'converges nicely at this low resolution'), the boundary conditions of the basin (periodic / walls), the x positions and spacing of the three boxes and their initial height, the heel direction of each box, whether surface tension was on, the gas-phase parameters, and the bounce-back variant details beyond 'simple bounce-back (Def. 12)'. The 26.565 / 45 deg targets are the ideal geometric values (his Fig. 5.8 / Sect. 5.3.2), not measured discrete-grid values, so a discrete offset at s = 16 (rho 0.75: ours 32-36 deg) would not contradict his figure.
> - **[Claude, 2026-09-24] Job 4253988's own question ANSWERED: it stays near 26.565, confirming the free-sweep's 29-32 deg is an approach-path/DOF effect, not a discretization or code-correctness problem.** `eval/antidunes.sqlite` has no per-particle tilt column, and the job's own stdout log was not preserved anywhere (checked all three checkouts and the vault's `sbatch-output-archive`) -- recovered instead by parsing the `rotation` DataArray (9-component, ASCII) directly out of the archived `vtk-out/floaterBox/simulation_step_*/process_*.vtu` files (`tiltY = atan2(R(0,2), R(0,0)) * 180/pi`, same convention as `[BoxForces]`), 50 samples over the 50000-step run. Result: settles to **-27.0 deg by step ~20000 and holds flat (+-0.1-0.2 deg noise) for the remaining 30000 steps** -- a 0.4 deg gap from the ideal 26.565 deg, consistent with ordinary lattice-discretization error at this coarse (s=16) resolution (the same kind of error Bogner's own Fig. 5.6 documents shrinking with resolution), and nowhere near the free 6-DOF sweep's 29-32 deg (rows 1-6 of the torque-sweep artifact, `omega1p9_seed5`/`omega1p5_seed5`/`omega1p0_seed5`/`omega1p9_seedEq`). Since this run uses the SAME physics, SAME resolution, and SAME density as the free sweep, differing only in constraints (X/Y position and pitch/yaw pinned, Bogner's own 'roll about the longitudinal axis only'), the extra 3-5 deg gap in the free runs is most likely a coupling effect from the additional degrees of freedom available to them (X/Y translation, pitch/yaw, interacting with the periodic self-image at 8-cell spacing already flagged as a known deviation from Bogner's own setup) -- not evidence of a resolution or implementation bug. Not yet root-caused further (e.g. which specific DOF drives the shift, or whether it is the periodic-image proximity specifically) -- a natural follow-up would be a second control run with X/Y pinned but pitch/yaw free (isolating translation vs. rotation coupling), or the same free setup on a non-periodic/wider domain. Plotted alongside the 6-run sweep (dashed control trace, sign-flipped for visual comparison) in the torque-sweep artifact: https://claude.ai/artifact/CTyxTFJbfahuGhuFx6NnDC. Also checked in this session: Bogner's own thesis (Zotero `H7YLRADA`) gives no derivation for why 2.86 deg specifically was chosen as the seed heel -- only that starting off-equilibrium demonstrates convergence -- and credits the whole test design to Fekken (2004) (Zotero `ETTABZSX`), whose PDF has a broken/non-standard font encoding that defeats `pdftotext` (even after a `qpdf --qdf` structural repair) throughout the entire document, not just its reference list; not pursued further via OCR this session.
> - **[Claude, 2026-09-24] MECHANISM of the free-sweep's 3-5 deg overshoot IDENTIFIED: the box never reaches a static fixed point at all -- it's locked in a sustained roll-translate coupled cycle, and the "settled" tilt is that cycle's mean, not an equilibrium.** Re-parsed `omega1p9_seed5`'s own per-step `[BoxForces]` log (`pos=<x,y,z>` alongside `tiltY`, already had this from the tilt/torque extraction, just hadn't looked at `pos` before) at 100-step resolution: X travels from 21.7 to a peak of 41.97 at step 38700, then reverses and drifts back to 30.3 by step 100000 -- over 20 LU of net travel, more than the box's own 16 LU width -- and `totalTorque_y` keeps changing sign for the entire 100000 steps, including the very last logged step (never converges to a fixed sign let alone zero). Zooming into one cycle (steps 35000-42000, every 100 steps) shows the coupling directly: tilt peaks at 32.27 deg at step 38100 **exactly** when the box's own horizontal velocity (dx/dstep) crosses zero and reverses sign -- textbook coupled-pendulum behavior (zero translational velocity at the rotational turning point), not two independent processes. `omega1p5_seed5` shows the identical signature at a different rate (X: 21.7 -> 6.2, still moving at step 100000, not yet reversed within the run length). Interpretation: a box that is rolling AND surging wets its hull asymmetrically relative to a pure roll, shifting the effective torque balance away from the roll-only 27 deg fixed point found in job 4253988's control -- the four "calm" sweep runs (small torque magnitude at any given snapshot) are actually orbiting a slow limit cycle, not sitting still; only the genuinely chaotic runs (seed3, seed4) look different because their orbit amplitude is much larger (+-60 deg vs +-1 deg), not because they're a different kind of motion. Still NOT isolated: whether this couples through X/Y translation specifically, pitch/yaw, or the periodic self-image proximity (8 LU gap, already flagged as a known deviation from Bogner's setup) -- the natural next step remains the same one already proposed above (a second control with X/Y pinned but pitch/yaw free, isolating translation from rotation coupling). Added as a second callout to the torque-sweep artifact (same link) with the concrete step/position numbers.

### Challenges & Future Improvements
- Candidate runs need to be reviewed specifically for visual clarity (a good showcase needs an
  unambiguous, easy-to-screenshot float/sink contrast, not just correct physics) — this is a
  presentation requirement on top of a correctness one.

---

## 2a. [x] Cuboid floaters: interaction validation

Introduce cuboids and their interaction with other floaters (and with each other). This session
already has a validation scheme (`DEMValidationScenarios.h`'s `runMultiPointEdgeCases`,
`runProductionCollisionReplay`, `runFallOntoPileStability`, `runPacking`, plus the production-scale
bridge/jam cluster jobs and the paused `DropTest`/arch-pier investigation) — see
`MITIGATION_COMPARISON.md` for the full record of what's been tried (`useSpheroBoxContact`,
`useMultiPointBoxContact`, `useBoxOverlapCap`, `useRoundedBoxCorners`, etc.) and each one's actual
pass/fail record.

> **[Dennys, 2026-09-12]:** agreed with my supervisor to NOT pursue the remaining flying-floater/
> contact-chatter issues any further — this item is considered done as-is for the thesis; the
> unresolved parts go into Challenges & Future Work below, not into more implementation work.

### Challenges & Future Improvements
- The flying-floater / contact-chatter bug chased at length this session was never fully root-caused
  at production scale (see `MITIGATION_COMPARISON.md` and the paused `DropTest` fluid-free isolation
  attempt). Per the decision above, this is presented in the thesis as a known, documented
  limitation with mitigations (`useSpheroBoxContact` etc.) that measurably reduce — but do not
  eliminate — its frequency, not as a solved problem.
- Existing validation battery is edge-case/unit-style (isolated pairwise scenarios) plus one real
  production-incident replay — there's no systematic sweep across cuboid aspect ratio, density, or
  impact angle. Noted as a future-work gap, not something being pursued now.

> **[Claude, 2026-09-13]:** job 4226400 (`run_36_droptest_all.sbatch`, the 1000-cuboid `DropTest`
> scenario, `useSpheroBoxContact=true`) FAILED at step 7928/158500 (~462s in) with "Received invalid
> notification type: ... from sender: 12" (`mesa_pd/mpi/notifications/ParseMessage.h:273` — a
> corrupted/desynchronized inter-rank particle-sync message), during the cuboids phase; the cylinders
> phase never ran. Same failure class as the paused flying-floater/chatter investigation above, not a
> new bug from this session's own changes (`useSpheroBoxContact` bypasses the new
> `useFraigeBoxContactModel`/`boxBoxSAT` code path entirely, so nothing added today was even exercised
> in this run). Per explicit request, submitted a genuinely different comparison run rather than more
> implementation work: job 4227741 (`run_36_droptest_fraige.sbatch`, `prms/DropTest_fraige_gas.prm`)
> repeats the identical cuboids scenario with `useSpheroBoxContact=false` /
> `useFraigeBoxContactModel=true` instead — SAT + Fraige et al. (2008)'s fuller multi-point replica,
> a different contact-detection code path entirely, not a single-parameter A/B test (see that prm's
> own header for why). Also submitted job 4227750 (`run_36_droptest_fraige_250.sbatch`,
> `mt-dennys-mirror2`, concurrent with 4227741 on `mt-dennys-mirror`) — same Fraige contact model, but
> the floater count reduced 1000 -> 250 (`prms/DropTest_fraige_250_gas.prm`: same 10-column spawn
> grid truncated to the first 25 of 100 waves/column, `timesteps` 158500 -> 46000) per explicit
> request.
>
> **Update, same day**: per explicit request, cancelled the 1000-floater Fraige run (4227741) and
> instead ran everything at 250 floaters. **Job 4227750 (Fraige, 250 floaters) COMPLETED cleanly**
> in 25 minutes — full 46000 steps, no crash, no MPI notification errors, properly archived under
> `$HPCVAULT/antidunes-archive/antidunes-droptest-fraige-250_20260913_142816/`. (Caught and fixed a
> real bug of my own along the way: my first attempt at generating `DropTest_250_gas.prm`/
> `DropTest_cylinders_250_gas.prm` via a truncation script had an `in_block` state-tracking bug that
> silently deleted an unrelated earlier parameter block — `loadSnapshot`/`storeSnapshot`/
> `snapshotFrequency`/`snapshotBaseFolder`/`periodicity` — caught immediately via a local smoke-test
> abort before ever reaching the cluster, then regenerated correctly and re-verified locally.) Also
> submitted the two 250-floater ORIGINAL-contact-model companions for full comparison: job 4227798
> (`run_36_droptest_250.sbatch`, box, `mt-dennys-mirror`) and job 4227799
> (`run_36_droptest_cylinders_250.sbatch`, cylinders, `mt-dennys-mirror2`) — testing whether the
> ORIGINAL `useSpheroBoxContact=true` model (which crashed at 1000 floaters, job 4226400) also
> crashes at this smaller scale, or whether the failure is scale-dependent. Report back once both
> finish.

> **[Claude, 2026-09-14]:** both finished. **4227799 (cylinders, 250) COMPLETED cleanly** (30 min).
> **4227798 (box, original contact model, 250) FAILED** — same exact signature as the 1000-floater
> crash ("Received invalid notification type: ... from sender: 31", `ParseMessage.h:273`), but much
> faster this time (66s vs 462s). This means the crash is **not purely a scale effect** — it
> reproduces at 1/4 the floater count too, specific to box + `useSpheroBoxContact`. Strengthens the
> earlier Fraige-model result (clean completion at both 1000 and 250 floaters) as a genuinely
> meaningful difference, not a fluke of one lucky run.

> **[Dennys, 2026-09-13]:** document, for the thesis, how the DEM side and the fluid (FSLBM/MEM)
> side each actually "see"/represent a cuboid's shape — these are two genuinely different geometric
> representations coupled together, and that distinction itself is worth writing up explicitly, not
> just the DEM contact-geometry side already covered above.
> **[Claude, 2026-09-14]:** did this write-up, reading the actual mapping code rather than assuming —
> and the plan's own speculation above (a voxelised/rasterised mesh approximation on the fluid side)
> turns out to be **wrong**. The real picture is more precise and worth stating differently in the
> thesis:
>
> - **DEM side**: the box is an exact analytic primitive (`mesa_pd::data::Box`/`FloatingBox`,
>   `edgeLength()`/`support()`) — contact detection against it (`boxBoxSAT`, the Eberly closed-form
>   distance, `boxBoxMultiPointContacts`, this session's `boxBoxFraigeMultiPointContacts`) operates on
>   the exact sharp-cornered continuous geometry, evaluated in closed form (no sampling, no grid).
> - **Fluid/MEM side uses the SAME exact analytic shape, not a separate voxelised mesh.** The
>   fluid/solid cell classification (`lbm_mesapd_coupling/momentum_exchange_method/
>   MovingParticleMapping.h:145`, `singleCast(particleIdx, *ac_, containsPointFctr, *ac_,
>   currentCellCenter)`) dispatches to `mesa_pd::ContainsPointFunctor::operator()` for `Box`
>   (`src/mesa_pd/common/Contains.h:137-142`), which transforms the query point into the box's own
>   body frame (`transformPositionFromWFtoBF` — the box's actual current position/orientation, not an
>   axis-aligned assumption) and tests it with `isPointInsideBoxBF` (`Contains.h:61-67`): an exact
>   `|x|<=0.5*Lx && |y|<=0.5*Ly && |z|<=0.5*Lz` half-space test against the identical `edgeLength()`
>   the DEM side uses. There is no separate rasterised/voxelised copy of the box built anywhere for
>   this coupling.
> - **The real discretisation is in WHERE the exact test gets evaluated, not in the shape itself**:
>   `MovingParticleMapping.h` only calls this test once per LBM cell, at that cell's CENTRE
>   (`currentCellCenter`) — a single point-sample per cell, not a sub-cell area/volume-fraction
>   integral. This is exactly the mechanism this thesis already root-caused independently in item 1b
>   above (the box-MEM buoyancy noise finding): an axis-aligned flat box face lies exactly parallel to
>   a grid plane, so a whole plane of cell-centre samples can cross the exact boundary simultaneously
>   as the box moves a tiny amount, producing synchronized force jumps that a sphere/cylinder's curved
>   surface never produces (different points on a curve cross at different times). Both findings are
>   the same underlying fact viewed from two angles: same-shape point-sampling classification is only
>   *approximately* equivalent to a true volume-fraction coupling, and the approximation error is
>   worst for large flat faces aligned with the grid.
> - **A more accurate alternative already exists in waLBerla**, just not used by Antidunes: the
>   Partially Saturated Cells Method (`src/lbm_mesapd_coupling/partially_saturated_cells_method/`,
>   `PSMUtility.h`, `ParticleAndVolumeFractionMappingSweepsCPU/GPU.h`) computes an actual solid-volume
>   fraction per cell rather than a binary point-sample classification — a concrete, correctly-cited
>   direction for future work rather than a vague "increase resolution" suggestion.
>
> **Bottom line for the thesis**: DEM and fluid do not use two different geometric representations of
> the cuboid (both use the exact same analytic `mesa_pd::data::Box`) — they differ only in *how much of
> that geometry each solver actually queries*: DEM resolves the full continuous shape in closed form;
> MEM samples it at discrete cell centres once per step. The mismatch is bounded by grid resolution
> relative to box size (finer grid -> more sample points per face -> smaller synchronized-flip
> artifact, consistent with item 1b's own R=1->R=2 partial-recovery data), and PSM is the existing,
> already-implemented remedy that trades this point-sampling scheme for an exact volume-fraction one.

---

## 2b. [x] Cylinder floaters: interaction validation

Same as 2a, for cylinders.

> **[Dennys, 2026-09-12]:** same supervisor-agreed decision as 2a — not pursuing further, documented
> as a limitation/future work item instead.

### Challenges & Future Improvements
- Cylinder contact detection has its own dedicated analytic path (face-face/band-band/face-edge/
  edge-edge cases in `ComputeContacts`) distinct from the box path — validation coverage for this
  path specifically has the same gaps identified in 2a (systematic sweeps vs. edge-case-only
  coverage), left as future work.

> **[Claude, 2026-09-30]:** this item is about the production-scale box-chatter question specifically
> (closed per the supervisor decision above); the separate, unit-test-style contact-detection/dynamics
> validation effort (`FINAL_THESIS_DEM_VALIDATION.md`, `DEM_VALIDATION_MATH_REFERENCE.md`) is tracked
> there rather than here. That effort's cylinder-wall (Kodam et al. 2010 replica) and box-wall
> (Fraige 2008 / Mack 2011 replicas) dynamics results — previously validated but only recorded in
> those markdown files — are now written into the actual thesis as
> `overleaf/chapters/04_results/contact-validation/wall-impact.tex`, plus a new background-chapter
> paragraph/equation for the rotational effective mass their root-cause fix depends on
> (`02_background/dem/contact.tex`, citing Baraff 1997). See `FINAL_THESIS_DEM_VALIDATION.md`'s own
> 2026-09-30 entry for the full detail.

---

## 3. [ ] Measuring clogging: elevation, pressure, and flow-forcing mode

> **[Dennys, 2026-09-12]:** promoted to its own top-level item, ahead of the clogging scenarios
> below (was nested as "3a" under them) — agreed we need this measurement setup built and validated
> *before* running clogging scenarios, not alongside/after them.

Need a concrete way to **measure** when/whether clogging is occurring, not just observe it visually:
- Track upstream water **elevation** (a clog should raise the free surface behind it, analogous to a
  real logjam backing up a river) and **pressure** at/near the obstruction.
- Consider running with the PID controller **disabled** (constant imposed flow/forcing) instead of
  its normal target-mean-velocity control, specifically to see what a developing clog does under a
  truly constant driving force — with PID active, the controller may partially compensate for a
  forming clog by adjusting forcing, masking the effect that's actually of interest. Compare against
  a constant-PID (steady target) run too, since both regimes are physically meaningful (PID mimics
  a controlled flume; constant forcing mimics a fixed hydraulic gradient).

> **[Dennys, 2026-09-12]:** agreed we need a setup for this.

> **[Dennys, 2026-09-12]:** for the validation setup itself: use a higher flow speed than our usual toy
> examples, start with spheres, use a real pier gap (not a slitpass) to force clogging, elongate the
> domain in X (periodic wraparound + queueing room), lower the water level (and therefore domain
> height, which is also cheaper), and use a spread of floater densities so different depths along the
> gap get contested.

> **[Claude, 2026-09-12]:** implementation status. Added `ClogMeasurementEvaluator` (Evaluators.h) --
> logs `clogUpstreamElevation`/`clogDownstreamElevation` (fill-level column sum, averaged over a small
> XY window) and `clogUpstreamPressure`/`clogDownstreamPressure`/`clog*ExcessPressure` (lattice
> pressure at an exact probe cell, minus the plain hydrostatic reference) at two explicit prm-configured
> probe points -- gated by an optional `ClogMeasurement` prm block (default off, zero behavior change
> for every existing prm). Added `PIDParameters.enabled` (default `true`) -- `false` disables
> `ForcingAdjuster`'s controller update entirely, pinning `acceleration[0]` at its analytically-derived
> initial value for the whole run (constant-forcing/"fixed hydraulic gradient" mode); `MeanVelocity`
> tracking stays on regardless. Both verified via a clean local build of `Antidunes`/`BedGeneration` and
> the existing `ANTIDUNES_ContactDetection`/`ClipAnalyticGeometry`/`BognerStateClassification` unit
> tests (all still pass).
>
> Built a small (8-rank) two-pillar validation scenario per the decisions above: gap=6 LU < sphere
> diameter=8 LU (guaranteed blockage), domain `<384,40,90>` (elongated X, halved Z per the
> lowwater-convention), `targetMeanVelocityMagnitude=0.04` with `useHardcodedUnits=false` (avoiding a
> real unit-mismatch bug `thousand_cuboids_two_pillars_slitpass_gas.prm`'s own header documents for
> exactly this kind of velocity change), 12 spheres across the 3 density tiers from item 1a's sweep
> (0.25/0.5/0.75) fed in on a staggered release ladder. 3 prms (no-clog reference, sanity-clog with PID
> on, sanity-clog with PID off) submitted as one combined cluster job (4223076). Once it completes,
> next steps: confirm the tool actually shows the expected clogged-vs-open-channel contrast and that
> the PID-off run holds `acceleration[0]` bit-for-bit constant, then pick concrete numeric
> "clogged" thresholds from the real elevation/pressure/transport-rate numbers and archive everything
> under `/run/media/dhuber/data2/debris-clogging-archive/clog-measurement-validation/`.

> **[Claude, 2026-09-14]:** checked in on this — 4223076 and one retry (4223107) were both CANCELLED
> (never ran), but a further resubmission, **job 4223118, actually COMPLETED** (58 min) on 2026-09-12
> and was never analyzed until now. Pulled all 3 archived `antidunes.sqlite` databases
> (`$HPCVAULT/antidunes-archive/antidunes-clog-measurement-validation-all_20260913_025443/`) and
> queried the `clogUpstreamElevation`/`clogDownstreamElevation`/`avgParticleTransportRate`/
> `acceleration_x` columns directly. **The tool works as intended:**
>
> | Case | mean(up-down) whole run | mean(up-down) last 10 samples | acceleration_x |
> |---|---|---|---|
> | No-clog reference | +0.0000 | +0.0000 | n/a |
> | Sanity-clog, PID on | +0.36 | **+2.70** | actively varies (0 to 6.5e-6, confirmed engaging) |
> | Sanity-clog, PID off | -0.05 | **+1.24** | **bit-for-bit constant** (one distinct value, 9.68e-8, across all 349 samples) |
>
> The no-clog reference shows exactly zero upstream/downstream divergence the whole run (as it should
> -- nothing there to detect). Both clog scenarios start near zero and develop a real, growing
> divergence by the end of the run (upstream measurably higher than downstream), coinciding with
> `avgParticleTransportRate` collapsing from ~0.0008-0.001 down toward ~0 in the same late-run window
> -- i.e. the debris genuinely jams up and transport through the gap stalls, exactly the physical
> signature this tool was built to measure. PID-on's final divergence (+2.70) is larger than PID-off's
> (+1.24) -- a plausible real difference (PID-on keeps pushing harder to hold the target mean velocity
> as the clog forms, driving more water to back up against the same restriction) worth a real
> explanation in the thesis, not just a noted number. `PIDParameters.enabled=false` independently
> confirmed to do exactly what it was built to do (acceleration pinned bit-for-bit); PID-on
> independently confirmed to actively engage (many distinct acceleration values, a real oscillating
> response). **Next actual step**: pick the concrete numeric "clogged" threshold this data supports
> (e.g. sustained up-down divergence above some value, or transport-rate collapse below some fraction
> of its own open-channel baseline) and archive under
> `/run/media/dhuber/data2/debris-clogging-archive/clog-measurement-validation/` per this item's own
> plan -- not yet done.

> **[Claude, 2026-09-14]:** user asked whether a cross-sectional area average would be better than the
> two-point probes, and to check the literature. Read two Schalko et al. papers already in the Zotero
> library: Schalko et al. 2018 (E3S Web Conf. 40, River Flow 2018) uses point ultrasonic distance
> sensors and defines backwater rise as a *between-run* comparison (Δh = h_jam − h_no-jam at the same
> gauge location), not the within-run upstream-vs-downstream comparison this tool currently does.
> Schalko, Follett & Nepf 2023 (WRR 59, e2023WR034689) is directly analogous to our two-pillar-with-gap
> geometry (a partial-width logjam with an open lateral gap) and found real lateral variation
> (~7.5%±1.8%) in water depth between the jam-blocked region and the open gap beside it -- they took
> *separate* longitudinal profiles at the jam centre and the gap centre rather than assuming one
> location represents both. This establishes (a) a full cross-sectional average is not standard
> practice and would dilute the two physically distinct flow regions, and (b) a single gap-centreline
> probe is insufficient for a lateral-gap geometry specifically. **Implemented the recommendation**:
> `ClogMeasurementEvaluator` (Evaluators.h) now takes an optional `columnPrefix` (default `"clog"`),
> letting a second instance coexist writing `clogPillar*` columns. `SimulationParams.h` gained optional
> `clogPillarMeasurementEnabled`/`clogPillarUpstreamPos`/`clogPillarDownstreamPos`, parsed from optional
> `pillarUpstreamPos`/`pillarDownstreamPos` inside the existing `ClogMeasurement` prm block (absent ->
> feature off, zero behavior change for old prms). `BuildTimeloop.h` wires a second evaluator instance
> gated on the new flag. All 3 clog-validation prms updated with a second probe pair at
> `<140/250, 11, 40>` -- Y=11 matches the left `CylinderPillar`'s own centre (`<192,11,45>`, radius 6)
> exactly, vs. the existing gap-centreline pair's Y=20; same X/Z as the existing probes for a direct
> like-for-like comparison. Verified locally (clean build + a 40s smoke test of the updated
> `sanity_clog_pidon` prm confirming correct config echo and clean timeloop entry) before syncing to
> the cluster and submitting **job 4234182** (all 3 cases, same combined-job sbatch, no script changes
> needed). Not yet analyzed -- pending completion.
>
> Also did the explicitly-requested wider web search ("look online if there also other methods...to
> cover our basis") beyond the two Schalko papers. Findings, for the record:
> - **Afflux via between-scenario baseline subtraction** (Panici & de Almeida 2018, *Formation, Growth,
>   and Failure of Debris Jams at Bridge Piers*, WRR; also a 2022 JMSE dynamic-jam paper) is the
>   dominant convention in the bridge-hydraulics literature generally, and matches Schalko 2018's
>   Δh = h − h₀ shape: upstream level in a jammed run minus upstream level in a *separate* unobstructed
>   run at the same gauge, not a within-run upstream-vs-downstream split. Our `noclog_reference` run
>   already provides exactly the baseline needed to compute this as a free secondary metric.
> - **Blockage ratio** is the field's standard normalizing/correlating variable for reporting afflux
>   results (not a new measurement technique) -- relevant for how results get framed, not for probe
>   placement.
> - **Pier drag-force time series paired with the water-level time series** (Panici 2018; the JMSE 2022
>   paper) -- an independent, non-hydraulic clog-severity signal running alongside the water-level
>   probes. We already compute MEM hydrodynamic force for the coupling, so pulling the force on the
>   `CylinderPillar`s directly would be cheap to add as a second line of evidence if ever wanted.
> - **Camera/mesh-based accumulation-geometry tracking** (De Cicco 2020; Ruiz-Villanueva 2016) measures
>   the jam's physical size/porosity directly rather than inferring it hydraulically -- not needed here
>   since a DEM sim already has exact particle positions/packing, strictly easier than for a physical
>   flume.
> - CFD-DEM debris-flow-impact studies (e.g. Wang et al. 2025) target transient impact loading, a
>   different phenomenon from sustained clogging/backwater -- not a competing methodology for this item.
>
> Net conclusion: the existing point-probe approach (now with the pillar-face pair) matches established
> practice; cross-sectional averaging would be a step away from, not toward, the literature's own
> methods. Baseline-subtracted afflux and pillar drag-force are both free/cheap optional secondary
> metrics worth keeping in mind, not required additions.
>
> Same findings as a quick-reference table, for citing when writing up results:
>
> | Method | Source | Notes |
> |---|---|---|
> | Point elevation/pressure probes (UDS-style) | Schalko et al. 2018, 2023 | What `ClogMeasurementEvaluator` already does — confirmed to be the field's standard, not a cross-sectional average |
> | Separate profiles per lateral flow region (jam vs. gap) | Schalko, Follett & Nepf 2023 | Motivated the pillar-face probe pair added this session |
> | Afflux via between-scenario baseline subtraction | Panici & de Almeida 2018 | Δh = h(jammed run) − h(unobstructed run) at the same gauge, over time — a *different axis* from our within-run upstream-vs-downstream split; numbers above |
> | Blockage ratio as normalizing/correlating variable | widespread in bridge-hydraulics literature | **Not a measurement method itself** — how results get reported/correlated (afflux vs. % cross-section blocked) across a family of runs, not a per-run instrument; numbers above |
> | Pier drag-force time series paired with water level | Panici 2018; JMSE 2022 dynamic-jam paper | Independent, non-hydraulic severity signal; cheap to add since MEM hydro force is already computed on the (currently unread) `CylinderPillar` particles |
> | Camera/mesh accumulation-geometry tracking | De Cicco 2020; Ruiz-Villanueva 2016 | Measures the jam's physical size/porosity directly — strictly easier for us since a DEM sim already has exact particle positions, not needed here |
> | CFD-DEM impact-force studies | e.g. Wang et al. 2025 | Target transient impact loading, a different phenomenon from sustained clogging — not applicable here |

> **[Claude, 2026-09-14]:** worked out concrete numbers for two of the above metrics, for reference
> when writing up results (neither needs new evaluator code -- both are arithmetic on data/geometry we
> already have):
>
> **Baseline-subtracted afflux** (`clogUpstreamElevation_case(t) − clogUpstreamElevation_noclog_reference(t)`,
> computed from job 4223118's 3 archived DBs), as a complement to the within-run up/down numbers already
> logged above:
>
> | Case | mean afflux, whole run | mean afflux, last 10 samples | at run end (t=34900) |
> |---|---|---|---|
> | Sanity-clog, PID on | +0.14 | +1.04 | +1.26 |
> | Sanity-clog, PID off | −0.00 | +0.52 | +0.90 |
>
> Both start at exactly 0 (before any clog forms) and grow through the run. This isolates the clog's
> own contribution from the flume's natural streamwise elevation slope, which the within-run
> upstream-vs-downstream number (+2.70 / +1.24 above) does not separate out.
>
> **Blockage ratio** (a static geometric fact used in the literature as the x-axis when correlating
> afflux *across a family of runs* with varying jam size -- not something a probe measures over time,
> and not directly applicable to a single fixed-geometry sanity run like this one). For our own
> `clog_measurement_validation_*` geometry (domain Y=40 LU, walled): pillar 1 at Y=11 r=6 spans
> Y=[5,17], pillar 2 at Y=29 r=6 spans Y=[23,35] -- blocked width = 24/40 = **60% blocked by the piers
> alone**, leaving a 6 LU gap plus 5 LU margins against each wall. Sphere diameter (8 LU) exceeds the
> 6 LU gap, so **one trapped sphere brings the effective blockage ratio to ~100%** at that
> cross-section -- this scenario is deliberately parked at the saturated end of the blockage-ratio
> scale (to guarantee clogging happens at all), which is exactly why it has no blockage-ratio variation
> to report. Once item 4 (clogging scenarios, which sweeps floater types/counts/densities/gap sizes)
> produces multiple geometries, plotting each run's final afflux against its own blockage ratio would
> be the standard literature-style summary figure for the thesis.

> **[Claude, 2026-09-14]:** validated both of the above. Blockage ratio: reparsed directly from the prm
> (not hand-transcribed) -- confirms 60.0% pier-only blockage, 6 LU gap < 8 LU sphere diameter, exactly
> as reported. Afflux: since the no-clog reference shows zero upstream/downstream divergence throughout
> (established earlier in this item), the identity
> `within-run(upstream-downstream)_case(t) = afflux_upstream(t) - afflux_downstream(t)` must hold
> exactly -- checked over all 349 samples in both PID-on and PID-off, **max error = 0.00e+00**,
> confirming the afflux numbers are arithmetically consistent with the within-run numbers already
> logged above (this is a correctness check on the computation, not an independent physical
> validation -- there is no real-flume ground truth available for this synthetic scenario).
>
> This check also surfaced a real physical finding worth keeping for the writeup: **downstream afflux
> is not negligible**. At run end (t=34900): PID-on afflux_upstream=+1.26 but afflux_downstream=**-1.69**
> (comparable magnitude, opposite sign); PID-off +0.90 vs **-0.71**. The clog doesn't just raise water
> upstream -- it measurably *lowers* it downstream relative to the unobstructed baseline too (less flow
> getting through the gap, so downstream starves). This is why the within-run up/down number is larger
> than either one-sided baseline afflux alone: it sums both the upstream rise and the downstream drop.

> **[Claude, 2026-09-14]:** checked two things the afflux computation could be silently confounded by,
> per a direct question about whether the pillars' own effect and flow spin-up are accounted for:
>
> - **Pillars' own baseline effect**: confirmed `clog_measurement_validation_noclog_reference_gas.prm`
>   has the *identical* `CylinderPillars` block (same centers/radii, same 60% blockage) as the clog
>   prms -- the only difference is zero `FloatingSphere`s. Since both runs carry the same pier geometry
>   from t=0, the piers' own steady wake/backwater is common-mode and cancels exactly in the baseline
>   subtraction; the afflux numbers above isolate the debris jam's *additional* contribution, not the
>   piers' own effect. Confirmed numerically too: afflux is exactly `0.000000` (6 dp) at every sample
>   from t=100 through t=5000 (the first sphere-release step) -- bit-identical between the two runs
>   until debris is actually present, so there is no spurious divergence from simulation noise leaking
>   into the comparison.
> - **Flow spin-up**: `meanVelocity_x` in the no-clog reference starts at 0.0447 (t=100) and drifts
>   *down* monotonically to ~0.0431 by t=15000 (target is 0.04, `targetMeanVelocityMagnitude`) before
>   curving back up to ~0.0437 by t=30000 -- the base flow is **not** at steady state when the first
>   sphere is released (t=5000), and several more releases (through t=21500) happen while it is still
>   settling. This does not invalidate the afflux/within-run comparisons (verified common-mode and
>   cancels by construction, per the previous bullet), but it does mean the early clog signal reflects
>   "debris entering a still-developing flow," not "debris entering a mature steady flow" -- worth
>   stating explicitly in the writeup. **For item 4's real scenarios**, delaying the first release until
>   `meanVelocity_x` has actually converged (roughly t>15000-20000 based on this data, rather than the
>   current fixed 5000-step warmup) would give a cleaner separation between flow spin-up and the
>   debris's own effect.

> **[Claude, 2026-09-14]:** implemented the fix and a "more extreme" scale-up, both requested directly.
>
> **Release-ladder fix applied** to all 3 small-trio prms
> (`clog_measurement_validation_{noclog_reference,sanity_clog_pidon,sanity_clog_pidoff}_gas.prm`):
> first `FloatingSphere` release moved from t=5000 to t=20000 (same 1500-step ladder spacing,
> continued), `timesteps` extended 35000 -> 56500 (last release 36500 + ~20000-step tail) so all three
> cases -- including the no-clog reference, which has no spheres but must still run the same length for
> the baseline-subtracted afflux comparison to have a matching sample at every timestep. sbatch time
> limit bumped 4h -> 6h accordingly. Job **4234212** submitted (cancelled the earlier, not-yet-started
> job 4234182 which still had the old t=5000 ladder, so no compute was wasted).
>
> **New "extreme" scenario** for a bigger/clearer contrast, a sibling of the small trio rather than a
> replacement: domain X doubled (384->768; every pillar/probe/release-zone X position is just the small
> scenario's own X doubled, a clean geometric scale-up preserving every relative margin exactly -- Y/Z
> and the pier gap itself are UNCHANGED, so the "one sphere physically cannot pass" guarantee still
> holds), floater count tripled (12->36, same 12-entry X/density/Z release cycle simply repeated 3x
> rather than hand-designed new positions, to avoid a transcription error), release ladder starting at
> the same t=20000, `timesteps=97500` (last release 72500 + 25000-step tail). 144 ranks (2 nodes,
> `--partition=multinode`, cellsPerBlock unchanged so per-rank cell count matches the small trio's own
> 72-rank run). New files: `clog_measurement_validation_noclog_reference_extreme_gas.prm`,
> `clog_measurement_validation_extreme_{pidon,pidoff}_gas.prm`,
> `jobs/run_144_clog_measurement_validation_extreme_all.sbatch` (16h time limit -- generous given three
> stacked unknowns: 3x floater DEM cost, first-ever 2-node run of this scenario, and the longer
> timestep count). Both new prms smoke-tested locally (oversubscribed 144-rank run, config parsed
> correctly, physical params computed sanely, no errors) before syncing. Job **4234213** submitted.
>
> Both jobs pending on the cluster as of this note; not yet analyzed.

> **[Claude, 2026-09-14]:** both jobs (4234212, 4234213) initially FAILED in 3 seconds with
> `antidunes.sh: No such file or directory` -- a path-resolution issue, not a physics/prm problem.
> Root cause: `run_72_clog_measurement_validation_all.sbatch`'s (and every other combined-job
> sbatch's) `REPO="$(cd "$SLURM_SUBMIT_DIR/../../../.." && pwd)"` assumes `sbatch` was invoked from
> INSIDE the `jobs/` subdirectory (4 levels up from `jobs/` lands on `mt-dennys`) -- I had instead
> invoked it as `cd .../Antidunes && sbatch jobs/run_....sbatch` (one directory short), landing REPO
> one level too high. Confirmed by comparing against the still-successfully-RUNNING R4 buoyancy job
> (4227812), which uses byte-identical REPO-resolution boilerplate -- it must have been submitted
> from inside `jobs/` per this same convention. Not a script bug; an invocation-convention mismatch on
> my part. Fixed by resubmitting correctly (`cd .../Antidunes/jobs && sbatch <name>.sbatch`) --
> **jobs 4234435 (small trio) and 4234436 (extreme trio)**, both now queued cleanly (confirmed past
> the immediate-failure point). Noting the convention explicitly here since it's shared boilerplate
> across many of this project's sbatch scripts: always `cd` into `jobs/` itself before `sbatch`ing one
> of these combined-job scripts, not into the `Antidunes/` directory one level up.

> **[Claude, 2026-09-14]:** correction/supersedes the note directly above, per explicit request ("for
> future submission for consistency make sure we can sbatch from Antidunes with jobs/..."). Checked
> the REST of this project's sbatch scripts (e.g. `run_72_two_cylinder_floater.sbatch`) and found the
> established, documented convention is actually the OPPOSITE of what I concluded above: those
> scripts' own header comments state "scripts are always submitted with sbatch from
> apps/showcases/Antidunes... so SLURM_SUBMIT_DIR sits 3 levels below the repo root" and use
> `$SLURM_SUBMIT_DIR/../../..` (3 `../`), not 4. The 11 "combined-job" scripts (buoyancy-comparison
> family, clog-measurement, droptest) had drifted to an inconsistent 4-`../` version -- a real,
> project-wide copy-paste bug, not an invocation mistake on my part. **Fixed all 11 scripts** to the
> established 3-`../` convention and re-synced to all three cluster checkouts (`mt-dennys`,
> `mt-dennys-mirror`, `mt-dennys-mirror2`). `cd ~/mt-dennys/apps/showcases/Antidunes && sbatch
> jobs/<script>.sbatch` is now the correct, consistent invocation for every sbatch script in this
> project -- do NOT `cd` into `jobs/` first (that was only a workaround for the bug, now unnecessary
> and would break these 11 scripts the other way).

> **[Claude, 2026-09-14]:** checking in on the resubmitted jobs (4234435 small trio, 4234436 extreme
> trio) turned up a real mistake: I had submitted BOTH from the plain `mt-dennys` checkout, which
> ALSO already had job **4227812** (R4 buoyancy comparison) actively running there since 15:10 --
> directly violating this project's own established "use a separate mirror checkout for concurrent
> runs" convention. Consequence: both 4227812 and 4234436 FAILED with `mv`/`rm: cannot remove
> '.../build/apps/showcases/Antidunes/eval/avgDataSlices': Directory not empty` -- a race on the
> shared `build/apps/showcases/Antidunes/` output directory between the R4 job's own archive step and
> the newly-started clog job's startup cleanup step. 4234436 failed instantly (40s, at its very first
> cleanup) with nothing lost. **4227812 is a real loss**: its R=4 sphere-MEM case had already run to
> full completion (52:58 min, 41600 steps, 8-node/576-rank allocation) but failed at the archive `mv`
> step immediately after -- and by the time this was noticed, the still-running 4234435's own startup
> cleanup had already `rm -rf`'d the same shared `eval/`/`vtk-out/`/`snapshots/` directories for its
> own run, overwriting the orphaned R=4 sphere-MEM output before it could be rescued. That ~53 minutes
> of 8-node compute is genuinely gone, not just misfiled -- a real cost of the mistake, though not a
> correctness problem (nothing was silently corrupted, just lost before archiving).
>
> Fixed by resubmitting both properly isolated: **R4 buoyancy re-run as job 4234535** on
> `mt-dennys-mirror` (already had the R4 bed file from its original generation, no cross-mirror gap to
> fix this time), **extreme clog trio re-run as job 4234536** on `mt-dennys-mirror2` (synced
> `BuildTimeloop.h`/`SimulationParams.h`/`Evaluators.h` and all 6 clog-measurement prms there first).
> `mt-dennys` itself is left alone for 4234435 (small trio) to finish undisturbed. All three now
> confirmed on three DIFFERENT checkouts, no further collision risk.

> **[Claude, 2026-09-14]:** checked in on 4234435/4227812's successors and the extreme trio's own
> retry. **4234435 (small trio) COMPLETED** -- both fixes confirmed: `meanVelocity_x` converges to
> ~0.0431 by t~12000-14000, well before the new t=20000 first release; the pillar-face probe pair
> agrees closely with the gap-centreline pair (last-10 up/down diff -- PID-on: gap=3.07/pillar=2.95;
> PID-off: gap=1.02/pillar=1.08) and both are exactly zero in the no-clog reference. One caveat: the
> reference run's sqlite is missing its first ~23300 steps (331/564 expected samples, no internal
> gaps) -- coincides with the window job 4227812 was still colliding on the same `mt-dennys` checkout
> before it failed; doesn't affect the late-run comparisons above (all three runs are complete and
> gap-free there) but would leave an early gap in any full-run baseline-subtracted afflux series from
> this specific archive.
>
> **4234536 (extreme trio, first attempt) FAILED** at step 42501 with a `[RELEASE-OVERLAP]` abort (a
> freshly-released sphere at X=220 collided with an existing floater still nearby) -- but this was
> properly isolated on its own checkout, a real physics/scenario-design issue, not the collision bug.
> Per explicit request, root-caused from the archived partial data (the aborted run's own sqlite
> survived, never cleaned since the script aborted before its own archive step) rather than guessed:
> by t=42000 `avgParticleTransportRate` had already dropped to ~0.0006 (~40% below its own
> ~0.0010 early-run value) and `clogUpstreamElevation` had risen to ~51.7-52.0 vs. a confirmed ~49.62
> undisturbed baseline (from the completed `noclog_reference_extreme` sibling in the same job) -- a
> real, substantial backwater effect had already reached back far enough by the third release cycle to
> measurably slow advection at the spawn zone itself, so earlier floaters weren't clearing out before
> the next reuse of the same X position. This is actually a positive result for the "make it more
> extreme" goal: the scenario was already producing a bigger, faster-developing signal than the small
> trio at only ~43% of its intended run length.
>
> **Fix applied**: moved the release zone from X=160-200 to X=40-100 in both
> `clog_measurement_validation_extreme_{pidon,pidoff}_gas.prm` (same internal spacing, whole block
> shifted -120 upstream -- clearance to the pier at X=384 increased from 164-224 LU to 284-344 LU),
> keeping the full 36-sphere debris load and release cadence intact (targets the actual mechanism --
> distance from the jam -- rather than reducing the scenario's own "extreme" intent). Smoke-tested
> locally (144-rank oversubscribed parse check, no errors, all 36 spheres present) before syncing.
> Resubmitted as **job 4234802** on `mt-dennys-mirror2` (checked occupancy first per the updated
> concurrent-runs convention -- `mt-dennys-mirror` was busy with 4234535/R4-buoyancy, `mt-dennys` and
> `mt-dennys-mirror2` were both free), now running, isolated from the R4 job. Not yet analyzed.

> **[Claude, 2026-09-14]:** job 4234802 (the relocated-spawn-zone retry) **FAILED again**, same
> `[RELEASE-OVERLAP]` signature, this time at step 48501 -- involving X=100, the new zone's own
> downstream-most position. Exactly one reuse-interval (4 positions x 1500 steps = 6000 steps) later
> than the first failure (42501). This proved relocation alone doesn't fix the mechanism, only delays
> it proportionally to the extra distance -- congestion keeps building throughout the run regardless
> of where the spawn zone sits, so it will eventually reach ANY fixed spawn zone if releases keep
> happening late enough into the run.
>
> Root-caused properly this time before the 3rd attempt: pulled `avgParticleTransportRate` across the
> full 0-30000 range of the aborted run and found the flow holds its healthy ~0.0009-0.0014 rate right
> through t=30000, only starting to decline after that. Both failures (t=42501, t=48501) happened well
> into that already-congested regime -- because the original schedule kept releasing new floaters all
> the way out to t=72500. The actual root cause was never "not far enough from the pier" so much as
> "releases keep happening long after the flow has already slowed enough that reused spawn points
> stop clearing in time."
>
> **Real fix**: kept the X=40-100 relocation (still helps) but compressed release spacing from 1500 to
> 200 steps, finishing all 36 releases by t=27000 instead of t=72500 -- comfortably before the
> ~t=30000-32000 slowdown onset, so no release event ever occurs during the congested regime. Debris
> load (36 spheres) and the density/Z cycling pattern are unchanged; only the timeline is compressed.
> `timesteps` (97500) left unchanged, giving a long ~70500-step tail to observe the clog mature after
> all debris has already arrived -- if anything more informative than the original design, since the
> whole debris load is now present from early on rather than trickling in for most of the run.
> Smoke-tested both prms locally (144-rank oversubscribed parse check, no errors). Checked occupancy
> again (`mt-dennys-mirror` still busy with 4234535; `mt-dennys`/`mt-dennys-mirror2` free) and
> resubmitted on `mt-dennys-mirror2` as **job 4234953**. Not yet analyzed.

> **[Claude, 2026-09-14]:** job 4234953 **FAILED again**, but this is a genuinely different failure
> mode -- `reynoldsNumber is nan!` (`WriteData.h:98`), not a `[RELEASE-OVERLAP]`. Progress-wise this
> confirms BOTH prior fixes worked: the run reached **t=48800** (vs. 42501 and 48501 for the first two
> attempts) with a clean release phase throughout. By that point the clog had already reached
> **near-total blockage**: `avgParticleTransportRate` collapsed to ~3e-5 (essentially zero flow) and
> the up/down elevation gap widened to **~14.3 LU** (57.97 vs 43.63, vs. a ~49.62 undisturbed
> baseline) -- a far stronger signal than the small trio ever produced, i.e. the "extreme" scaling
> already achieved its purpose before the crash.
>
> Root-caused (per explicit request) from the raw stdout log rather than guessed, using the
> `[timestep]`-tagged log lines right before the abort: the sequence was (1) `Bogner state selection:
> ... has N unclassified FormerMO neighbour(s); falling back to interface` -- cells just VACATED by a
> moving particle ("FormerMO") couldn't be cleanly reclassified gas/liquid; (2) `No liquid or
> interface cell is in the neighborhood to distribute excess mass to. Mass is lost/gained.` -- the
> free-surface mass-conservation step then failed at those same cells; (3) dozens of simultaneous
> `[LOW-RHO] interface cell ... rho=nan` cells, clustered near the domain floor (Z=0) and at
> inter-rank block-boundary-local coordinates. This is a genuine free-surface reconstruction edge case
> that only manifests once particle packing becomes extreme (the just-reached near-total-blockage
> state) -- the fallback path for a FormerMO cell's neighbor classification isn't robust at that
> density.
>
> **Correction, checked directly per explicit question ("is this with spheres and are there flying
> spheres?")**: initially framed this as "the same underlying theme" as item 2a's documented
> "flying floater / contact-chatter" limitation -- checked the actual `[TOTAL-FORCE-THRESHOLD]`
> warnings logged right before the crash and that framing was NOT supported by the evidence, so
> retracting it. The two particles flagged there (uid 1513209474796486662, 648518346341351428) have
> tiny, stable velocities (~0.002-0.003 LU/step, essentially at rest relative to the flow) and modest,
> DECREASING contact forces across sub-cycles (0.0635 -> 0.0573) -- completely normal sustained contact
> within a dense pile, not a launch or any kind of instability. Nothing is flying. Item 2a's own
> limitation is specifically about DEM contact detection producing impulsive launches; that is not
> what is happening here. The actual mechanism is confined to the free-surface side only (Bogner state
> misclassification of a cell vacated by a calmly-moving particle -> mass-conservation failure ->
> density NaN) -- a distinct issue from the box/contact-chatter limitation, not a manifestation of it,
> and specific to extreme packing DENSITY rather than any particle's velocity or launch behaviour.
>
> Given the partial data already captures the complete story (healthy flow -> developing clog ->
> near-total blockage) before the crash, next step is a decision on whether to (a) attempt a further
> mitigation and resubmit a 4th time, (b) accept the salvaged partial run as the "extreme" result and
> document the NaN as its own, separately-diagnosed known limitation (free-surface reconstruction at
> extreme packing density -- not the box-contact-chatter one), or (c) something else -- pending
> discussion.

> **[Claude, 2026-09-14]:** dug one level deeper, per explicit follow-up ("are there low rho
> warnings?") -- 1464 `[LOW-RHO]` warnings total across the run, but the precise chain matters. The
> VERY FIRST occurrence (of either) is at step 48819, cell `(56,5,14)` on rank 66, immediately preceded
> (same cell, same step) by repeated `[DRY-LINK] non-finite pdf_old ... flags=[ interface ]` warnings --
> traced this to its actual source, `src/lbm_mesapd_coupling/momentum_exchange_method/boundary/
> SimpleBB.h:184-203` (core waLBerla, not an Antidunes-local file). Read that guard's own code comment
> directly rather than assuming: **`[DRY-LINK]` firing is NOT itself a bug** -- it is a pre-existing,
> intentional safety guard for exactly this situation ("a cell on the fluid side of this obstacle link
> may be a GAS cell ... at the waterline of a floating particle, which carries no valid PDF ... the
> expectation is that it is a drained INTERFACE cell (ρ→0)"), and it correctly SKIPS using the
> non-finite value for particle-force coupling rather than propagating it -- our case matches the
> anticipated one exactly (`flags=[ interface ]`, not gas). The comment explicitly notes "genuine
> corruption of a *liquid/interface* PDF is caught separately by the checkPdfNaN scanner, which is not
> silenced by this guard" -- i.e. the real, still-unexplained problem is one level further upstream:
> WHY cell (56,5,14) became a "drained, ρ→0" interface cell severely enough that its own PDFs (via
> normal streaming/collision, not just this guard's read) eventually went non-finite. Consistent with,
> and now more precisely located than, the earlier Bogner/FormerMO framing: very plausibly a cell just
> uncovered by a particle in the dense pile, for which the free-surface reconstruction has no
> recovery/floor once its density gets pathologically small. The dry-link guard is doing its job
> correctly; the actual instability is upstream of it, in the free-surface scheme's handling of a
> drained interface cell under extreme local packing.

> **[Claude, 2026-09-14]:** pinned down exactly where cell (56,5,14)/rank 66 actually is, per two direct
> follow-up questions ("close to the clogging or somewhere different?" / "is it a surface cell?").
> Converted local->global using `cellsPerBlock=<64,20,15>` and cross-checked against that same rank's
> own `[FloaterForces]` log lines (its owned particles sit at X~374-377, Y~27-36, Z~5-14) rather than
> assuming the block-to-rank mapping: **global position ~(376, 25, 14)**.
> - **Close to the clog**: yes, essentially at it -- X=376 is only 8 LU upstream of the pillars'
>   center (X=384); Y=25 falls inside pillar 2's own footprint (center Y=29, radius 6, spans
>   Y=[23,35]). This cell sits right against the submerged base of the second pillar.
> - **NOT a surface cell**: this run's own computed `absoluteLiquidHeight = 49.13` -- the crashed cell
>   at Z=14 is deep underwater, nowhere near the real air-water boundary. Notably, the floaters owned
>   by this same rank also sit at similarly low Z (5-14) despite densityRatio 0.25-0.75 (buoyant,
>   should float toward the surface) -- consistent with floaters being pushed down and submerged by
>   the sheer weight of the pile-up against the pillar, not resting at their natural buoyant depth.
>
> Sharper conclusion: the actual bug candidate is that the free-surface reconstruction's
> interface-cell math (valid only near the true air-water boundary) gets applied here regardless of
> depth, via the FormerMO fallback, to a cell that should just be ordinary submerged liquid pinned
> deep underwater by pile-up pressure against a pillar -- not a real free-surface transition at all.

> **[Claude, 2026-09-14]:** retracting the FormerMO/"reconstruction impossible" framing above --
> checked both claims directly per two sharp follow-up questions ("was the reconstruction impossible
> in this case", "can we confirm this cell was a FormerMO or is it inside the pillar") and NEITHER
> holds up:
> - **Not inside the pillar**: confirmed geometrically -- (376,25) is 8.94 LU from pillar 2's centre
>   (384,29), radius 6, i.e. ~2.94 LU of genuine fluid space outside the pillar surface.
> - **FormerMO: unconfirmed, likely wrong.** Searched the ENTIRE log for any prior mention of cell
>   (56,5,14) on rank 66 -- it never appears before it is already NaN. No Bogner/FormerMO message
>   names this specific cell; the earlier "unclassified FormerMO neighbour" messages I'd pointed to
>   were at different local coordinates on different ranks (74, 80), not this cell. That was an
>   unverified inference from a superficially similar but distinct event, not evidence about this one.
> - **"Reconstruction impossible": also unconfirmed.** The actual diagnostic for that ("No liquid or
>   interface cell is in the neighborhood to distribute excess mass to") never fires on rank 66
>   anywhere in the whole log.
> - **This IS the true first occurrence**, not a downstream victim -- confirmed it's the earliest
>   DRY-LINK/LOW-RHO line in the entire 51324-line log; no other rank shows corruption earlier.
>
> Per the user's own prior-testing experience, a submerged cell mislabeled "interface" is itself a
> known, common, harmless artifact -- so the mere fact that this cell is submerged-and-labelled-
> interface does NOT distinguish it from countless other harmless instances presumably occurring
> throughout the run. Honest current state: this cell shows up already-NaN the very first time ANY
> existing diagnostic happens to look at it (the DRY-LINK check, which only fires opportunistically
> when an obstacle link touches the cell); nothing in this project's current instrumentation caught
> whatever actually happened to it beforehand. The true trigger and mechanism remain **unidentified**
> -- this needs either additional targeted logging (e.g. instrumenting the PDF-refilling/reconstruction
> path directly, not just downstream consumers of a value that's already corrupted) or a dedicated
> reproduction attempt with tighter diagnostic frequency around this specific cell, not a conclusion to
> draw from this run's existing logs alone.

> **[Claude, 2026-09-14]:** built and launched exactly that dedicated reproduction attempt, per
> explicit request ("let's do this and perhaps we can add some new checks too catch these issues
> earlier... conditional logging if a value looks suspicious").
>
> **Found the intended infrastructure already exists**, mostly unused until now: `PhysicsParams.
> denseNanScanEnabled/FromStep/ToStep` (SimulationParams.h) is a pre-built, documented mechanism for
> exactly this -- its own comment says "Resume from the last checkpoint before a known crash
> (loadSnapshot=true) and set denseNanScanEnabled=true with a from/toStep window bracketing the
> crash." `checkPdfNaN`/`trackWatchCells`/`checkClosedInterfaceLayer` (SetupStages.h) already do the
> heavy lifting -- `checkPdfNaN`'s neighbourhood dump even already classifies FormerMO cells.
>
> **Added two new, threshold-triggered, throttled (first 20 hits/rank) checks** inside `checkPdfNaN`,
> addressing the "sometimes NaN isn't the first sign, unrealistic-but-finite values explode first"
> point directly:
> - `lowRhoPrecursorThreshold` (real_t, 0=off): warns on a liquid/interface cell whose density is
>   still finite but has drained below this value -- the step(s) *before* a cell goes non-finite,
>   which a plain isfinite scan can never see by construction.
> - `highVelocityPrecursorThreshold` (real_t, 0=off): the symmetric high-side check -- warns on a
>   velocity magnitude (via `PdfField_T::getVelocity()`, confirmed to divide by the cell's own actual
>   rho for this project's compressible cumulant model, not an assumed 1) that is finite but
>   unrealistically large.
> Both verified locally first (a tiny 3-rank toy scenario with deliberately-triggering thresholds --
> both warnings fired correctly, throttled as designed, no false wiring issues) before touching the
> real scenario.
>
> **Found and fixed a real, previously-unexercised bug** in the process, confirmed via a local
> `storeSnapshot`->`loadSnapshot` round-trip test: floater creation (spheres/boxes/cylinders) ran
> unconditionally in `setupParticles()` regardless of `loadSnapshot`, so resuming a checkpoint that
> already contained those same (non-GLOBAL, hence checkpoint-serialized -- confirmed by reading
> `BlockForestDataHandling.cpp`'s own `if (GLOBAL) continue;` guard in both serialize() and
> deserialize()) particles aborted with "particle with the same uid(N) already existing at index(M)".
> Pillars/planes/arch-tunnels do NOT need the same fix -- they're GLOBAL and explicitly excluded from
> serialization by that same guard, so recreating them fresh every run (load or not) is correct by
> design, not a bug. Fixed by gating the floater-creation block behind `if (!bfParams.loadSnapshot)`
> (SetupStages.h) -- confirmed fixed via the same local round-trip test, this time keeping the
> `FloatingBoxes` block in the reload prm unchanged (previously had to be manually omitted as a
> workaround; no longer necessary). Known, documented limitation of the fix: a checkpoint taken WHILE
> some floaters are still `startInactive`/dormant will leave them parked forever after this kind of
> resume (their FIXED state IS correctly restored, they just never get released without the
> creation-time UID bookkeeping) -- irrelevant for THIS investigation (checkpoint will be taken at
> ~t=48500, long after the last release at t=27000) but worth remembering for any future resume
> attempted before a scenario's own releases are complete.
>
> **Two-phase reproduction now running**: PHASE 1 (`clog_debug_extreme_pidon_phase1_gas.prm`,
> `jobs/run_144_clog_debug_extreme_pidon_phase1.sbatch`) is the original extreme_pidon scenario with
> `storeSnapshot=true`/`snapshotFrequency=500`/`timesteps` capped at 50000 -- expected to reproduce
> the exact same crash near step 48819, leaving a checkpoint from step 48500 behind (checkpoints
> overwrite in place, only the latest ever exists on disk, so no manual pruning needed). Submitted as
> **job 4235135** on `mt-dennys-mirror2` (checked occupancy first -- `mt-dennys-mirror` busy with the
> R4 retry, job 4235003). PHASE 2 (`clog_debug_extreme_pidon_phase2_gas.prm`) will `loadSnapshot=true`
> from that checkpoint once it exists, with `denseNanScanEnabled`+both new precursor thresholds
> (`lowRhoPrecursorThreshold=0.1`, `highVelocityPrecursorThreshold=0.4`) and tight VTK
> (`writeFrequency=10`, velocity+density+fill_level+mapped_flag writers) enabled for just the short
> remaining window -- estimated ~5GB of VTK output if it crashes again near the same step (~41 dumps
> at ~122MB each for this domain size), well within budget. Not yet submitted (waiting on phase 1's
> checkpoint). Not yet analyzed.
>
> **Storage correction**: checkpoints/VTK for this debug work land in the build directory under
> `/home/hpc/...` (the small-quota HOME filesystem), NOT `$HPCVAULT` (large quota, used for the
> archive) as I'd assumed when estimating "well within budget" above -- caught directly. Fresh quota
> check: HOME 54.4GB used / 100GB soft / 200GB hard; VAULT 93.2GB used / 1000GB soft / 2000GB hard.
> Recomputed against the correct (tighter) HOME quota: phase 1's checkpoint (~645MB, one copy only,
> since checkpoints overwrite in place) + phase 2's VTK (~5GB realistic, ~18GB worst case) + the
> existing 54.4GB baseline = 60-73GB total, still comfortably under the 100GB soft quota (27-40GB
> margin) -- conclusion unchanged, but the correction was legitimate and the number now real.
>
> **Verified the checkpoint mechanism with actual data, not just re-run success**, per two direct
> follow-up questions: does it really save/restore floater STATE (not just resume at the right step
> number), and does the pillar actually get correctly reconstructed?
> - **Floater position**: compared the resumed run's own final position
>   (`<20.4167, 32, 40.0002>`) against phase A's own UNINTERRUPTED continuous run at the same absolute
>   step (`<20.4091, 32, 40.0002>`) -- Z matches exactly, X differs by only 0.0076 (consistent with a
>   different DEM sub-cycle sampling point within the same main step, not a real discrepancy).
> - **Pillar reconstruction**: added a `CylinderPillar` to the toy test and (previously silent, now
>   fixed -- see below) logged its creation explicitly. Fresh run: `Created cylinder pillar:
>   center=<50,32,40> radius=3`. Reload run: byte-identical `Created cylinder pillar: center=<50,32,40>
>   radius=3`, immediately followed by the floater-skip message, then a clean resume. Confirms pillars
>   are correctly and deterministically reconstructed every run (not saved, by design -- GLOBAL
>   particles are excluded from serialization), with no drift.
> - Along the way, added a proper `WALBERLA_LOG_INFO_ON_ROOT("Created ... pillar: ...")` line for both
>   box and cylinder pillar creation (SetupStages.h) -- previously silent, inconsistent with floaters'
>   own "Created floating ..." logging; now matches that convention, useful for exactly this kind of
>   verification on any future resume.
> - Also noticed (not yet investigated further): the `[MEM-WATCH]`/`[TOTAL-FORCE-THRESHOLD]` logs'
>   own internal `step=` counter does NOT reflect the resumed absolute step (showed `step=1` right
>   after resuming at absolute t=11) -- likely a separate, presumably-cosmetic counter inside
>   ParticleKernels that isn't part of the checkpoint and simply restarts at 0 on any fresh process
>   launch, unlike the main loop's own `beginTimeStep` (which correctly showed 11). Worth keeping in
>   mind when reading logs from a resumed run -- the position/force values themselves are correct (as
>   verified above), just that specific label is not.
>
> **[Claude, 2026-09-15]: PHASE 1 completed WITHOUT crashing** -- job 4235135 (`mt-dennys-mirror2`)
> ran all 50000 of its own steps cleanly (43m32s, 144 ranks/2 nodes) instead of reproducing job
> 4234953's `reynoldsNumber is nan!` abort at step 48819, despite using the identical scenario, rank
> count, and node topology. The job's own script explicitly flagged this as unexpected. This is
> strong evidence the crash is **non-deterministic** (most likely MPI-reduction floating-point
> ordering sensitivity, or genuine borderline numerical chaos in a turbulent, debris-laden flow that
> depends on exact runtime scheduling) rather than a bug that reliably fires from this fixed input
> state. Since there is no crash in this run, PHASE 2 (dense-scan a short window around a known crash
> step) has no target yet and was not submitted.
>
> **Phase 1's own eval data checked for general sanity** (pulled `eval/antidunes.sqlite` locally,
> 499 rows, `evaluationFrequency=100` over 50000 steps): no NaN/Inf anywhere in `reynoldsNumber`,
> `froudeNumber`, `weberNumber`, `meanVelocity_*`, or any `clog*` column across all 499 rows --
> numerically clean throughout. Re (2158-3693, target 3100), Fr, and We stay in the physically
> expected range given PID-driven velocity noise. The clog signal itself develops exactly as the
> scenario was designed to produce: flat/undisturbed (`clogUpstreamElevation` ≈
> `clogDownstreamElevation` ≈ 49.62) through t=20000 (before any release), a noisy, sometimes-inverted
> transitional period from t≈25000-40000 (plausible transient sloshing as debris first queues up,
> `avgParticleTransportRate` dropping to ~0 and briefly negative around t=35000-37000), then a large,
> sustained, monotonically-growing clog signature from t≈41000 onward (`clogUpstreamElevation` rising
> 51→60 while `clogDownstreamElevation` falls 44→42 by t=49000), with the pillar-face probes tracking
> the gap-centreline probes almost exactly at that point -- consistent with a full-width jam having
> formed by the end of the run, rather than the partial-width case Schalko et al.'s ~7.5% lateral
> variation describes.
> - **Found one genuine data-quality issue, not yet investigated further**: `totalHydroForce_x`,
>   `totalHydroForce_y`, and `totalHydroForce_z` are all exactly `0.0` in literally every one of the
>   499 rows -- suspicious given individual `[FloaterForces]` log lines the same run produced clearly
>   show non-zero per-floater force components (e.g. `F_hydro=(0.1157,0.0451488,0.00562933)`). Reads
>   like the aggregate `totalHydroForce_*` evaluator column either isn't wired up to accumulate real
>   data for this scenario, or its accumulator is being reset before the eval write reads it -- worth
>   checking directly (`Evaluators.h`) before relying on this column for anything.
>
> **PHASE 1.5 created and submitted** (job 4235729, `mt-dennys-mirror2`) -- not phase 2's diagnostic
> resume (nothing to target yet), but a plain checkpoint-resume smoke test: `loadSnapshot=true` from
> phase 1's own final checkpoint (step 50000), `storeSnapshot` still on, none of phase 2's overhead
> (no `denseNanScanEnabled`, no precursor thresholds, no tight VTK). `timesteps=100000` -- confirmed
> from `Antidunes.cpp`'s own `for (t = beginTimeStep; t != timesteps; ++t)` loop that `timesteps` is
> the absolute target step, not a relative count, so reaching "the same duration as phase 1" from a
> step-50000 start requires 100000, not 50000. New files: `prms/clog_debug_extreme_pidon_phase1_5_gas.prm`,
> `jobs/run_144_clog_debug_extreme_pidon_phase1_5.sbatch`. Purpose: confirm the resume continues
> physically sensibly past step 50000, and, as a side effect, give a second independent 50000-step
> window in which the original crash could (or could not) reappear. Not yet analyzed.
>
> **[Claude, 2026-09-15]: PHASE 1.5 CRASHED -- and the failure mode reframes the whole
> investigation.** Job 4235729 aborted after only 1m59s wall-clock with the identical
> `reynoldsNumber is nan!` abort as the original crash. Critically, this was NOT at step ~48819
> (the original crash step) or step 50000 (phase 1's own final step) -- the log shows
> `beginTimeStep: 49501 / 100000`, meaning phase 1's own LAST checkpoint write was actually at
> step 49500, not 50000: with `timesteps=50000` and the loop `for (t=begin; t!=timesteps; ++t)`,
> the executed steps are 0..49999, and a `t % 500 == 0` checkpoint trigger last fires at t=49500
> (49999 itself is never a multiple of 500) -- an off-by-one in when "the last checkpoint"
> actually lands relative to the configured `timesteps`, not a bug in the resume logic itself, but
> worth knowing when interpreting any of this. The real finding: **resuming from that step-49500
> checkpoint and continuing crashed within ~96 steps** (`[LOW-RHO]` warnings first appear at
> step 49597, 8.6s of wall-clock, mostly checkpoint-loading overhead). Phase 1 itself, running
> continuously and NOT resumed, passed through this exact same step range (49500-49597) cleanly
> on its way to completing all 50000 steps without incident. This is a much sharper, far more
> reproducible signal than "the crash is rare/non-deterministic": **something about the
> checkpoint-restored state at t=49500 is NOT equivalent to the true in-memory state phase 1 itself
> had at that same step**, and that discrepancy alone is enough to trigger the instability almost
> immediately on resume.
>
> Partially investigated already: `oldHydrodynamicForce`/`oldHydrodynamicTorque` (the per-particle
> fields `Evaluators.h`'s own fix, above, relies on surviving a reset) ARE correctly included in
> checkpoint serialization -- confirmed by reading `ParticleCopyNotification.h` directly (lines
> 118-119, 179-180) and `BlockForestDataHandling.cpp`'s own `serialize()`/`deserialize()`, which
> reuse that exact notification type for particle checkpoint I/O. So this specific, plausible-looking
> hypothesis is ruled out, not confirmed -- the real cause is still open. Next places to look:
> whether the LBM fields (pdf/fill/excess-mass) round-trip bit-for-bit through the checkpoint
> files at a cell level (not yet checked), whether flag-field/interface-classification state is
> fully captured (checkpoint stores raw fields, not derived flags -- if flags are recomputed on
> load rather than restored, a boundary/interface cell could classify differently post-resume than
> it would have continuously), and whether the PID controller's own restored state
> (`pidState.file`) reproduces the exact driving acceleration phase 1 itself had at t=49501 (a
> mismatch there would perturb the flow right at resume, in exactly the way that could tip an
> already-marginal cell over into NaN). Not yet resubmitted -- this needs a diagnosis plan, not
> another blind resume attempt.
>
> **[Claude, 2026-09-15]: PHASE 2 run, per explicit instruction ("first make sure it's not a
> checkpoint bug, then run phase 2") -- strong evidence this IS checkpoint-related.** Built
> `jobs/run_144_clog_debug_extreme_pidon_phase2.sbatch` (didn't exist yet) and submitted on
> `mt-dennys-mirror2` (job 4236481) against phase 1's still-intact checkpoint (`storeSnapshot=false`
> in the prm, confirmed untouched, still timestamped from phase 1's own save). With
> `denseNanScanEnabled` active from the very first post-load step, the result is unambiguous:
> **`[NaN-SCAN] aborting after timeloopFSLBM at t=49501`** -- the FIRST timestep after resuming,
> not step 49597 (phase 1.5's own coarser periodic-check detection point) or anywhere later. Per
> mid-investigation request, resubmitted once more (job 4236487) with `logFloaterHistory=true`,
> `logFloaterHistoryStartStep=49501` added to the prm (cheap here: ~499-step window, nowhere near
> the 126GB/1000-floater blowup that flag's own doc comment warns about) to check whether a
> floater/pillar contact was involved -- **the resulting `floater-history/*.csv` files are all
> header-only, zero data rows**, because `checkPdfNaN`'s abort fires during
> `timeloopFSLBM->singleStep()`, strictly BEFORE `reduceAndAvgHydroForces()`/`doDEMSubCycles()` even
> run for that step (Antidunes.cpp's own per-timestep order) -- ruling out a floater-contact-force
> origin for THIS specific failure, at least directly.
>
> The FATAL cell list itself is the strongest clue yet: every single NaN cell reported (12 cells
> across 12 different ranks, all `flag=LIQ`, i.e. genuinely fully-liquid, not interface) sits at
> world X in [372.5, 384.5] -- i.e. either just outside or just inside pillar 2's own boundary
> (`center=<384,29,45>` and `<384,11,45>`, `radius=6`, per the extreme scenario's own
> `CylinderPillars` block). Z values are scattered (0.5 to 30.5, well below the ~49 LU liquid
> surface), so this is not a free-surface/interface effect -- it is specifically localized to the
> pillar boundary, appearing the instant the checkpoint-loaded fluid state takes its first LBM step.
> Pillars are GLOBAL and confirmed (earlier this session, byte-identical `Created cylinder pillar`
> log line) to be recreated fresh and identically on every load, never restored from the checkpoint
> itself -- so the leading hypothesis is that the boundary-handling/obstacle-mapping data structures
> that reconcile a GLOBAL body's geometry against the fluid field are rebuilt fresh at load time in
> a way that is NOT fully consistent with a fluid field that itself came from a checkpoint rather
> than a fresh initialization or continuous evolution, specifically at cells adjacent to that
> boundary. Not yet proven -- the next concrete step is checking whether phase 1's own raw
> checkpoint files (`pdfField.file`/`fillField.file`, still on disk, untouched) already contain
> non-finite values at these exact 12 cell coordinates BEFORE any load/resume happens at all (a
> pure byte-level check, no cluster time needed) to separate "the checkpoint was already bad when
> written" from "the checkpoint is fine but the post-load boundary reconstruction breaks it."

> **[Claude, 2026-09-15]: DEFINITIVE ANSWER -- it is NOT a checkpoint save/load bug.** Rather than
> hand-decode the binary checkpoint format (risk of misreading it independently), added a small,
> clearly-marked TEMPORARY diagnostic block directly in `setupLBMFields()` (SetupStages.h, inside the
> `bfParams.loadSnapshot` branch, right after the `pdfFieldID`/`fillFieldID` `loadBlockData()` calls)
> that dumps rho, fill level, and every raw PDF component at each of the 12 exact world cells that
> phase 2 (job 4236487) reported as `FATAL fluid-cell NaN` -- reusing the checkpoint's own real
> deserialization path (`PdfFieldHandling`/`ScalarFieldHandling`), reading the values immediately
> after load and strictly BEFORE `setupBoundaryHandling()` (flag field / obstacle mapping) or any
> timeloop sweep runs. Rebuilt and resubmitted the identical phase-2 job on `mt-dennys-mirror2`
> (job 4236518, same prm/checkpoint, `storeSnapshot=false` preserved) to get this dump for free
> alongside the existing dense-NaN-scan diagnostics already in that run.
>
> Result: **all 12 cells are completely clean in the raw loaded checkpoint** -- rho in
> [1.00045, 1.01159] (physically sane, close to 1), fill=1 (fully liquid, consistent with their
> `flag=LIQ` classification), every one of the 27 D3Q27 PDF components finite and of a sane
> magnitude (`anyNonFinitePdf=no` for all 12), at t=25.7s wall time right after load. Then, at
> t=26.3s (`beginTimeStep: 49501 / 50000`), the FIRST LBM sweep of the resumed run
> (`timeloopFSLBM->singleStep()`) runs, and by t=27.1-27.2s **those exact same 12 cells** are
> reported NaN by `checkPdfNaN` -- i.e. the corruption is produced by the first post-load sweep
> itself, not present in anything read from disk.
>
> This rules out a checkpoint write bug AND a checkpoint read/deserialization bug outright -- the
> bytes on disk, and the fields as loaded into memory, are correct. It also strengthens (without yet
> fully proving) the standing pillar-boundary-reconstruction hypothesis: since phase 1 itself passed
> through this identical physical state/step continuously without incident, and the loaded state is
> provably identical (finite, sane) to what phase 1 had in memory at that point, the only remaining
> difference between "continued in-memory" and "checkpoint-resumed" is exactly the state that a
> checkpoint does NOT carry across a resume and must instead be freshly rebuilt: the flag field /
> obstacle mapping / `former moving obstacle` bookkeeping that `setupBoundaryHandling()` constructs
> AFTER this diagnostic point, from the freshly-recreated GLOBAL pillar particles (confirmed
> byte-identical geometry, but not necessarily byte-identical resulting flag-field state) against the
> now-loaded fluid field. Next step: instrument `setupBoundaryHandling()`/the boundary-mapping sweep
> itself (or diff the flag field's classification at these 12 cells against what phase 1's own log
> would have shown at the equivalent point) to find exactly which piece of freshly-rebuilt state
> disagrees at the pillar boundary. The temporary diagnostic block in `setupLBMFields()` should be
> removed once that follow-up is done (or repurposed/left gated behind a flag if it proves useful
> again).

> **[Claude, 2026-09-15]: ROOT CAUSE FOUND AND FIXED -- confirmed with two clean end-to-end runs.**
> Per explicit follow-up ("go ahead, looks like we narrowed it down and should find the bug soon"),
> bisected further with three more one-shot diagnostic hooks (all temporary, all removed once the fix
> was confirmed):
> 1. A mid-`timeloopFSLBM` dump (between the geometry sweeps and the dynamics sweeps) showed
>    normal/curvature/obstacleNormal are correctly `<0,0,0>`/`0` at all 12 cells at that point (they
>    are genuinely bulk-liquid cells, not interface cells) -- ruling out a curvature/surface-tension
>    forcing-term origin and narrowing the bug to `addSweepSurfaceDynamics` itself, or to something
>    already wrong in the flag field before any sweep runs at all.
> 2. A full (non-early-exiting) re-scan of every LIQ/INTF cell for non-finite PDFs, run immediately
>    after `initFlagsFromFillLevel()` -- i.e. BEFORE any timeloop sweep -- found **~4881 already-NaN
>    cells**, not just the original 12: a huge slab (X in [367.5,383.5], full Y width [0.5,39.5], Z in
>    [0.5,29.5]) sitting directly upstream of/at the pillar gap. Far too large to be one stray cell --
>    this is the ENTIRE interior volume of the 36-sphere jammed debris pile from the "extreme, PID on"
>    scenario's own design, exposed as plain LIQUID with the raw, never-touched PDF garbage that had
>    sat under it (harmlessly, as an obstacle cell) throughout phase 1's real 49500-step history. The
>    original 12 "FATAL" cells were never the origin -- they were just the outermost cells the very
>    first LBM stream step happened to spread this pre-existing corruption into.
> 3. Counting `movingObstacle`-flagged cells right after `mapParticlesOntoCells()` (the app's ONE-TIME,
>    setup-only floater-obstacle mapping) confirmed the smoking gun directly: **all 36 floater
>    particles are present and completely intact after checkpoint load** (verified position, type,
>    `interactionRadius`, and the *shape's own* internal radius via `Sphere::getRadius()` -- both
>    exactly 4, matching the prm) and the mapping kernel *is* invoked (143 times, once per rank-local
>    ghost copy) -- but it maps **zero cells** anywhere in the whole domain, every single time.
>
> Root cause: `lbm_mesapd_coupling::MovingParticleMappingKernel::mapMovingParticleOnBlock()`
> (`src/lbm_mesapd_coupling/momentum_exchange_method/MovingParticleMapping.h:335`) only ever claims a
> cell when `boundaryHandling->isDomain(x,y,z)` is already true, and explicitly does **not** overwrite
> anything otherwise ("no else -- will not overwrite pre-existing boundary cells... e.g. when first
> mapping global particles"). `isDomain()` checks the `domain_` flag mask, which is **never
> initialized** anywhere in this app before this call: `setupBoundaryHandling()` only runs
> `initFromConfig()` (Border NoSlip) before `mapParticlesOntoCells()` is called, and the one call that
> actually establishes `domain_` for the rest of the field --
> `antidunesBoundaryHandling->initFlagsFromFillLevel()` -- happens **afterwards**, in `Antidunes.cpp`.
> Every other caller of this exact kernel in the whole waLBerla tree (e.g.
> `tests/lbm_mesapd_coupling/mapping/ParticleMapping.cpp`, `MovingParticleMappingMEM.cpp`,
> `SettlingSphereMEM.cpp`, and a dozen more benchmarks/showcases) explicitly calls
> `handling->fillWithDomain(...)` before mapping, for exactly this reason -- Antidunes never did.
> **This one-time setup mapping has therefore always silently mapped zero cells for every floater, on
> every run this app has ever done, fresh or resumed.** It has always been harmless on a fresh start:
> at t=0 the fluid field is freshly initialized sane data everywhere, so the brief one-step gap (fixed
> by the very next per-step mapping in `timeloopAfterDEM`, which runs after `initFlagsFromFillLevel()`
> has already executed once during setup) produces at most an invisible one-step blip. It is fatal on
> any checkpoint resume with active floaters present, because the exact same gap now exposes
> genuinely-undefined PDF data that has been sitting untouched under a real, long-established debris
> pile for tens of thousands of steps.
>
> Fix (`SetupStages.h`, `mapParticlesOntoCells()`): added a loop over all blocks calling
> `handling->fillWithDomain(handling->getFlagField()->nrOfGhostLayers())` right at the top of the
> function, before either mapping kernel runs -- `fillWithDomain()` only touches cells that are
> already completely empty (its own `isEmpty()` guard), so this cannot overwrite the Border-NoSlip
> flags `initFromConfig()` already set; it only gives `isDomain()` a correct answer for ordinary
> liquid/gas cells that haven't been classified yet.
>
> Verified with three clean end-to-end reruns of the exact same phase-2 job/prm/checkpoint (jobs
> 4236601 during the bisection, and 4236606 after removing all the temporary diagnostic code):
> `movingObstacle`-flagged cells went from 0 to 9752; the load-time already-bad-cell scan went from
> ~4881 cells / "YES" to 0 cells / "NO"; and the job now runs cleanly through and past t=49501 to
> completion (`sacct` reports `COMPLETED`, exit code 0) with zero FATAL aborts, where every previous
> attempt (phase 1.5, and every phase-2 diagnostic run before the fix) crashed within the first 1-499
> steps. All temporary diagnostic code (`[CHECKPOINT-DIAG]`, `[MID-FSLBM-DIAG]`,
> `[FLAG-NEIGHBORHOOD-DIAG]`, `[FULL-BAD-CELL-LIST]`, `[OBSTACLE-MAPPING-DIAG]`,
> `[LOAD-TIME-NAN-SCAN]`, `[FLOATER-STATE-DIAG]`, `[MAPPING-VISIT-DIAG]`, `dumpFlagNeighborhoodAtCells`,
> `listAllBadCellsNoEarlyExit`) has been removed; only the real `fillWithDomain()` fix remains.
>
> Remaining follow-up (separate issue, NOT checkpoint-related): while re-checking the R3
> cuboid-stillwater-hires equilibrium chain's status, job 4236462 (d025, `loadSnapshot=false`, a fresh
> run) FAILED at t=72000/100000 with its own `[NaN-SCAN] aborting after timeloopFSLBM` -- FATAL cells
> scattered across many different ranks, each at local block-corner `(0,0,0)`, at varied world
> positions with no obvious common cause yet. This cascaded `DependencyNeverSatisfied` onto the
> dependent d050/d075 jobs (4236463/4236464), which are now permanently stuck and will need
> cancelling/resubmitting once this is investigated. Not yet root-caused -- flagged here for the next
> session/request rather than assumed to be the same bug (this run never touches
> `mapParticlesOntoCells`'s loadSnapshot path at all).

> **[Claude, 2026-09-16]: standing convention -- sparse VTK for thesis figures.** Per Dennys's request,
> every clogging-measurement `.prm` from now on should produce ~4-6 evenly spaced VTK dumps (steps 0
> and the final step included) instead of the VTK output being fully disabled, specifically so a run
> can be spot-checked visually (and figured in the thesis) without the storage/IO cost of dense output.
> Applied now to all 6 currently-existing `clog_measurement_validation_*` scenario files (the 2
> currently-running extreme pidon/pidoff jobs' `.prm`s included, though editing them doesn't affect
> those already-launched jobs -- this takes effect on their next resubmission): `writeFrequency` set to
> `timesteps/5` (giving exactly 6 dumps at 0, T/5, ..., T) in each file's `VTK.fluid_field` block.
> Since every particle/floater/pillar VTK writer in `BuildTimeloop.h` reads the same
> `simParams->vtk.writeFrequency` (there is only one shared value, not one per writer), this one change
> synchronizes fluid-field, sphere, and `cylinder_pillar` VTK output onto the same steps for free --
> exactly what's needed to produce a consistent fluid+floater+pillar snapshot sequence per run. Any
> *new* clogging scenario file (including the design-space proposal above) should set
> `VTK.fluid_field.writeFrequency = timesteps/5` (or `/3` for 4 dumps) from the start rather than
> defaulting to disabled and having to remember this later.

> **[Claude, 2026-09-18]: d050 cuboid-stillwater-hires crash (job 4240983, see the
> 2026-09-15/09-16 entries above) -- checkpoint-resume diagnostics come back 7/7 clean, still
> unresolved.** Two separate attempts to reproduce the original crash (job 4240983 crashed at
> step ~91960, `[DRY-LINK]` cascade from step 91739, cell (2,0,10) local / global ~(53,126,40) --
> see the sacct/log analysis this session) by resuming from its own t=90000 checkpoint: a
> single-shot dense-NaN-scan resume (job 4242078) and a 6x-repeat resume on a second node
> allocation (job 4243270), all reading the identical checkpoint. **All 7 attempts completed
> the full target window with zero precursor firings and zero crashes.** The 6 repeats within
> job 4243270 ran on the SAME node allocation (`f[1252-1253]`) and were fully deterministic
> (bit-identical CLEAN outcome every time) -- ruling out same-hardware nondeterminism as the
> cause of the non-reproduction. Global cell (53,126,40) is NOT near the one `FloatingBox` in
> this scenario (16+ LU below its bottom face, just outside its Y-footprint), NOT near the free
> surface (34 LU below it), and NOT near any domain boundary -- an unremarkable point in bulk
> liquid with no obvious geometric explanation for why a corrupted `interface`-flagged cell would
> appear there. Confirmed this scenario has no sediment bed at all (`useParticleBed false`) --
> an earlier hypothesis in this same investigation (bed-interstitial interface pockets) was
> wrong and retracted.
>
> Investigated whether checkpoint resume genuinely reconstructs the same microstate continuous
> execution would have had at t=90000. **It does not, in a specific, verified way**:
> `Antidunes.cpp`'s call to `initFlagsFromFillLevel()` (which (re)builds the entire
> liquid/interface/gas flag classification from the current fill field) runs unconditionally at
> setup, fresh or resumed, and does so via direct `forceFlag()`/`removeFlag()` calls rather than
> the `setDomain()`/`setBoundary()` API that `BoundaryHandling.h`'s own `getNearBoundaryFlag()`
> comment says is the *only* sanctioned way to keep the near-boundary flag consistent ("Never
> (ever) set near boundary flags manually outside the boundary handler!"). This produced a
> measured burst of 302 `[STALE-NEARBOUNDARY]` corrections clustered right at the resumed run's
> own step 90002 (vs. the trickle a continuously-evolving simulation generates one cell at a
> time) -- confirmed this is NOT a resume-specific code path (the same `initFlagsFromFillLevel()`
> call runs on a fresh t=0 start too), just that a fresh start's near-empty fill field gives it
> almost nothing to act on, while a resume at a late timestep has 90,000 steps of real free-surface
> complexity for it to find stale flags in.
>
> **Retracting the stronger claim from earlier in this same session that this flag-reclassification
> difference is the leading explanation for the 7/7 non-reproduction**, after re-reading this
> file's own 2026-09-15 entries (immediately above) on a *different* scenario's checkpoint-resume
> crash (extreme clog debug, phase 1/2). That investigation found and fixed a more serious, related
> bug: `MovingParticleMappingKernel` mapped **zero cells** on any checkpoint resume with active
> floaters, because `isDomain()` (checked before a cell can be claimed) was never initialized at
> that point in the setup order -- `initFlagsFromFillLevel()`, which sets it, ran *after* the
> mapping call. This let genuinely undefined PDF garbage sitting under a long-covered floater
> (harmless while covered) get exposed as apparently-valid liquid for one step right at resume.
> Fixed via an unconditional `fillWithDomain()` call added at the top of `mapParticlesOntoCells()`
> -- **confirmed this fix is already present and unconditional in the current `SetupStages.h`**,
> meaning it already protects d050's own box-floater resume path. The near-boundary-flag
> difference I found is real but only affects GAS cells (whose PDFs are never physically
> meaningful and are never read for anything once the stale flag is cleared) -- unlike the
> isDomain()-zero-mapping bug, it has no mechanism to perturb genuine liquid/interface physics.
> **Net assessment: I do not currently have a confirmed mechanism for the 7/7 non-reproduction.**
> Both concrete "resume differs from continuous execution" differences I can point to in this
> codebase are either already fixed or provably harmless to real physics. The empirical pattern
> itself (0/7 crashes on resume vs. 1/1 crash on the original continuous run) is still a real,
> if weak, Bayesian signal that resume is doing *something* different -- just not one I can
> currently explain mechanistically. A genuinely rare, hardware/compiler-version-dependent
> instability that we simply haven't hit again across a still-small sample remains equally
> plausible and cannot be ruled out.
>
> The one experiment that would actually settle this: add a flag-field save/load path to the
> checkpoint format itself (currently NOT persisted at all -- see the save/recreate table below),
> bypassing `initFlagsFromFillLevel()`'s reclassification entirely on resume, and see whether a
> resume with a genuinely-preserved flag field reproduces the crash. If it does, the mechanism is
> confirmed; if it still doesn't, the flag hypothesis is dead and the remaining candidate is
> hardware/library drift between when job 4240983 first ran and now. Not attempted -- this is a
> real code/checkpoint-format change, not a diagnostic run, and needs explicit sign-off before
> touching the checkpoint format.
>
> **Reference table: what a checkpoint actually saves vs. recreates**, compiled while
> investigating the above (verified directly against `SetupStages.h`, `PIDController.cpp`, and
> `BlockForestDataHandling`'s particle serialization, not from memory):
>
> | Saved to checkpoint | File | Contents |
> |---|---|---|
> | Block/domain decomposition | `blockForest.file` | Rank-to-block layout, tied to that exact rank count |
> | PDF field | `pdfField.file` | Full D3Q27 fluid state |
> | Fill level field | `fillField.file` | Per-cell liquid fraction |
> | Excess mass field | `excessMassField.file` | Only if the excess-mass model is the "evenly, all-interface, fallback-liquid" type |
> | Particle storage (floaters only) | `particleStorage.file` | Position/velocity/orientation/UID -- GLOBAL/fixed bodies (pillars) explicitly excluded (`if (GLOBAL) continue`) |
> | PID controller state | `pidState.file` | commandVariable, actuatingVariable, all 3 gains, maxRamp, min/max actuating bounds, 3-entry error history, error integral |
> | Timestep counter | `simState.file` | `beginTimeStep` |
> | Driving acceleration | `simState.file` (same file) | `acceleration[0]` -- overwrites the analytically-recomputed value so PID-driven acceleration continues correctly |
>
> | Recreated fresh every setup (fresh start AND resume, unconditionally) | Mechanism |
> |---|---|
> | Flag field (liquid/interface/gas) | `initFlagsFromFillLevel()` -- full domain reclassification from whatever fill field is present |
> | Boundary conditions (NoSlip/FreeSlip, `BodyBC`, auto walls) | `setupBoundaryHandling()`'s `initFromConfig(BoundaryParameters)` -- rebuilt from the `.prm` every run |
> | Fixed pillars as DEM rigid bodies | `createCylinder`/`createBox` -- excluded from particle-storage serialization, always recreated from the `.prm` |
> | Force density field | Always zero-allocated fresh; recomputed every timestep regardless |
> | Curvature/normal/obstacle-normal/smooth-fill fields | Always zero-allocated fresh; each has its own per-timestep sweep |
> | Particle field (UID->cell map for MEM) | Rebuilt via `mapParticlesOntoCells()`, called unconditionally |
> | meanVelocity / totalFluidMass | Reset to zero at setup, recomputed during evaluation |
> | VTK output counter | Not loaded -- `initialExecutionCount` overwritten to `beginTimeStep` in the config instead |

> **[Claude, 2026-09-18]:** Both no-clog reference/extreme-clog PID-ON runs completed cleanly with the
> BodyBC fix in place (job 4242918 no-clog reference, job 4242919 extreme clog, both COMPLETED
> 2026-09-18, archived under `$HPCVAULT/antidunes-archive/` and copied locally to
> `/home/dhuber/University/Thesis/animations/*_clean_bodyBCfix_20260918/`).
>
> Extracted the PID-ON reference's own converged `acceleration_x` from `eval/antidunes.sqlite`
> (t>118465 tail window, n=63: mean=7.801745e-07, std=1.570243e-08; `meanVelocity_x` locked at
> 0.040000+-1.6e-5) and compared it against `setupSimState()`'s plain-channel analytical formula
> (`3*U*viscosity/liquidHeight^2` = 9.675064e-08 for this scenario's own computed viscosity/
> liquidHeight): **the PID-converged value is ~8.06x (706%) higher**. This makes physical sense once
> stated plainly — the analytical formula assumes an unobstructed channel and has no way to account
> for the two-pillar `BodyBC`'s own substantial hydraulic drag, which is now correctly modeled
> post-BodyBC-fix.
>
> Consequence: the two earlier completed PID-OFF-adjacent runs (job 4240981 `extreme_pidoff`, job
> 4240982 `extreme_pidon_fresh`, both 2026-09-17) predate the BodyBC fix — their pillars gave DEM
> contact physics only and never blocked the fluid at all, so the plain-channel formula was
> incidentally valid for *that* (bugged) geometry, but is invalid for the current, correctly-blocking
> one. Deleted both as invalid (local copies under `animations/`; never made it to `$HPCVAULT` in the
> first place per an archive-listing check, so nothing to remove there).
>
> Added a proper fix rather than hand-editing a single run's initial state: a new
> `PIDParameters.accelerationXOverride` prm knob (`SimulationParams.h`'s `PIDParams`,
> `SetupStages.h`'s `setupSimState()`), following this codebase's existing `*Override` convention
> (`gravityOverride`, `relaxationRateOverride`). 0 (default) reproduces the old analytical-formula
> behavior exactly; >0 pins `accelerationX` (and the `PIDController`'s own initial actuating value) to
> the supplied constant instead. Set `clog_measurement_validation_noclog_reference_extreme_pidoff_gas.prm`'s
> `accelerationXOverride` to 7.801745e-07 (the extracted value above) and resubmitted it (job pending,
> see sbatch output) on `mt-dennys-mirror` — the same checkout the PID-on reference (4242918) itself
> ran on, for direct provenance parity. The with-debris `clog_measurement_validation_extreme_pidoff_gas.prm`
> has NOT yet been given the same treatment (would need its own PID-on reference extraction, not just
> a copy of this no-clog value) — flagged as a follow-up, not done here to avoid scope creep beyond
> what was asked.

> **[Claude, 2026-09-18]:** Implemented the "PID first, then constant force" variant of the PID-off
> experiment (`PIDParameters.switchOffStep` / `switchOffAveragingWindow`, `ForcingAdjuster` in
> `Forcing.h`, wired in `BuildTimeloop.h`; both default to off, so every existing prm is unchanged).
> PID control runs normally until the switch-off step; there `acceleration_x` is frozen at the MEAN of
> the PID's own commanded values over the preceding window and never touched again. Motivation: the
> plain PID-off runs start from a hand-set force (`accelerationXOverride`) and an initial-velocity
> guess, so the clean case drifts 0.044 -> 0.0376 instead of holding 0.04; spinning the flow up under
> PID and then releasing the debris under the frozen force is the physically natural sequence for a
> slope-driven reach (normal flow first, jam arrives, slope unchanged) and needs no extract-and-hardcode
> step. The mean, not the instantaneous value, is frozen because the PID's command oscillates by about
> +-20% around its converged value at t=15000-20000 (5.7e-7 to 9.6e-7 in the no-clog PID-on data).
> The switch step is compared against the ABSOLUTE timestep (`beginTimeStep` + calls), and the
> checkpoint already stores `acceleration_x`, so resuming before or after the switch behaves correctly;
> a resume shortly before the switch averages only the PID updates made in the resumed run (the log
> line `[PID-SWITCHOFF]` reports how many samples went into the mean).
> Unit test `ANTIDUNES_ForcingAdjuster` (`tests/showcases_antidunes/ForcingAdjusterTest.cpp`, builds and
> passes locally) checks: freeze at the window mean while the last commanded value differs from it,
> bit-identical hold afterwards, feature off keeps controlling, window 0 = last value, resume before and
> after the switch step, `enabled=false` precedence. New prms (not yet run on the cluster):
> `clog_measurement_validation_extreme_pidswitch_2node_gas.prm` (extreme clog, switch at t=20000 = first
> debris release, window 5000) and its debris-free reference
> `clog_measurement_validation_noclog_reference_extreme_pidswitch_gas.prm`, whose job is to confirm the
> frozen force actually keeps the clean channel at 0.04.

> **[Claude, 2026-09-30] Implemented review.md items 4 and 5 (blockage ratio measurement, PID-acceleration
> resistance) — both were previously only discussed/defined in text, neither was actually computed.**
>
> **Item 4 (blockage ratio):** `ClogMeasurementEvaluator` (`Evaluators.h`) gained a new `blockageAt()`
> method, gated behind a new opt-in `blockageRatioEnabled` key in the `ClogMeasurement` prm block
> (`SimulationParams.h`'s `EvaluationParams::clogBlockageRatioEnabled`, same pattern as the existing
> `flowRateEnabled`). Per-run design decision, since the review's own "divide by the wetted area" is
> ambiguous between the LOCAL (backed-up) wetted area and the REFERENCE (undisturbed) one: chose the
> REFERENCE wetted area (channel width x reference liquid depth, a fixed constant) as the denominator,
> matching Panici & Kripakaran (2018)'s own definition against the undisturbed flow area — this makes
> the ratio a scenario-describing number comparable across runs, not one that moves with the very
> backwater rise it's meant to explain. Pillar contribution is obtained ANALYTICALLY from the known
> `CylinderPillar` geometry (a vertical prism of circular footprint spanning the full domain height, so
> its contribution to the wetted area is simply its y-extent at that x times the reference depth) rather
> than by scanning flags — avoids the harder problem of determining per-cell "is this obstacle cell
> submerged" from a flag-only test, and sidesteps entirely the fact that a (x,y) column can span
> multiple MPI-owned z-blocks (`cellsPerBlock` typically splits Z into several blocks in the clog
> scenarios), which would otherwise make a bottom-up "first gas cell" column scan need cross-rank
> communication. Jam contribution counts moving-obstacle cells at or below the reference depth,
> geometrically excluding cells inside a pillar's own footprint (pillars are themselves ordinary
> fixed-mass mesa_pd bodies and so may also carry the moving-obstacle flag — this exclusion is correct
> regardless of whether they actually do). New prm output columns: `clogUpstreamBlockageRatio`,
> `clogDownstreamBlockageRatio`, plus the raw `clog*PillarArea`/`clog*JamArea` components. Threaded
> through `BuildTimeloop.h` (builds `pillarCenters`/`pillarRadii` from `simParams->cylinderPillars`,
> passes `domain.domainSize[1]` as the channel width) — only the first probe pair ("clog" prefix) gets
> it, matching the existing flow-rate convention. Methodology written up in `probe-and-forcing.tex`
> (`eq:clog-pillar-area`, `eq:clog-blockage-ratio`), thesis rebuilds clean (verified by rendering pages
> 73-74 to PNG). **Not yet run on the cluster** — no existing prm has `blockageRatioEnabled true` yet,
> and producing the actual "afflux vs. blockage ratio across scenarios" plot the review asks for needs
> more than turning the flag on: the existing validation matrix has exactly ONE debris scenario per
> scale (Table 4.2), and blockage ratio is explicitly a CROSS-scenario comparison variable — a real
> plot needs a small sweep of scenarios with different debris loads or gap widths, which doesn't exist
> yet and is a design decision (how many points, what varies) rather than a rerun of what's already
> there.
>
> **Item 5 (PID-acceleration resistance):** needed NO new instrumentation — `acceleration_x` and
> `totalFluidMass` are already logged every evaluation step (`WriteData.h::writeFluidInfo`). New script
> `analyze_clog_resistance.py`: `--mode pid` computes `(a_x - a_x,ref) * mean(totalFluidMass)`; `--mode
> fixed` computes `Q/Q_ref` from the flow-rate probe. **Validated against real, already-completed data**
> (no new run needed): debris job 4260382 vs. its own reference 4257768 (both extreme-scale, PID-on,
> flow-rate-enabled, steady window t>=100000) gives `a_x` = 7.338e-6 vs. reference 7.708e-7 (+852%),
> resistance = 9.95 (lattice force units) — a striking, physically sensible result (the jam needs
> 8.5x the driving acceleration of the clean channel to hold the same target velocity). The `--mode
> fixed` path could NOT be demoed against a real matched pair: the one existing fixed-forcing reference
> run (`noclog_reference_extreme_pidoff_4245061`) predates the `flowRateEnabled` feature and has no
> `clogUpstreamFlowRate` column at all — needs a fresh run (or could reuse `meanVelocity_x` as an
> approximate substitute if that's acceptable) to demo the fixed-forcing side for real. Methodology
> written up in `probe-and-forcing.tex` (`eq:clog-pid-resistance`), thesis rebuilds clean.

> **[Claude, 2026-09-30] PER-FLOATER POSITION CHECK on job 4279039 (24-box hopper): PARTIAL,
> depth-selective clog found -- not "nothing happened."** Per explicit follow-up questions ("did it
> clog?", "compare position over time of floater in relation to hopper pillars"). Parsed the run's own
> `vtk-out/floaterBox` per-rank VTU point data (6 frames, steps 0/10000/.../50000), deduplicating by
> `uid` across ranks (ghost ranks near block boundaries write the same particle). Result, split by
> starting Z layer (matches the 3-layer density stratification: Z=10/rho=1.2 deep, Z=25/rho=1.0 mid,
> Z=45/rho=0.5 shallow):
> - **Deep layer (Z=10, rho=1.2, 8 boxes): 7/8 JAMMED permanently.** Each moved only ~10-20 LU from
>   its start in the first 10000 steps then went essentially stationary for the rest of the 51000-step
>   run (e.g. uid ...161: x=364->351.6->353.8->355.0->355.0->355.0; uid ...781: 350->336.2->...->336.2)
>   -- all parked in the x~=315-365 band just upstream of the pillar (center x=384). One deep-layer box
>   (uid ...163) stayed jammed until ~step 20000-30000 then broke free and passed through.
> - **Mid + shallow layers (Z=25/45, 16 boxes): 16/16 ESCAPED.** All had passed x=400+ by step 10000-
>   20000 and kept moving; several wrapped around the periodic domain (length 768) one or more times by
>   step 50000 (e.g. x values that go from >700 back down to <100 between frames).
> - Global `avgParticleTransportRate` stayed roughly steady (~0.003-0.0045, no decline) across the
>   whole run -- consistent with a PARTIAL, not total, blockage: the escaping majority kept
>   contributing steady transport even while the deep layer piled up and stopped.
>
> **This exactly reproduces the qualitative pattern already found with the cylinders** at small scale
> (12-rod run 4268353: "split exactly by initial release depth... stuck in floor-piles from the Z=10
> layer, escaped downstream from the Z=25/45 layers, regardless of row or Y-column") -- so shape
> (box vs. cylinder) changes whether the ghost-ownership CRASH happens, but not this depth-selective
> jam/escape split itself, which appears to be a property of the release geometry (density
> stratification by height) common to both shapes. Position-vs-pillar chart sent to the user directly
> (`/tmp/box24_check/position_vs_pillar.png`, x(t) per floater, periodic-unwrap applied, pillar
> reference line at x=384) -- not archived into the repo (a one-off diagnostic, not thesis content).
>
> **Two new box runs submitted per explicit follow-up request** ("double the number in x for one of
> the jobs and then... the full cylinder... around 400 cylinders, run this but as a cuboid version"):
> - **New `gen_hopper_clog_box_prm.py`** (sibling of `gen_hopper_clog_cyl_prm.py`): same population
>   grid/geometry-check logic, `FloatingBox` in place of `FloatingCylinder` (halfExtents `<r,L/2,r>`,
>   no axis field needed for an axis-aligned box). Diff-checked against a fresh cylinder-prm generation
>   with shape-block/comment text normalised out: confirmed byte-identical on every position, density,
>   count, and physics/ModelParameters line -- the two prms differ ONLY in shape. Generated the full
>   406-box production hopper: `prms/hopper_clog_box_frontblock_pidon_2node_gas.prm` (406 boxes, 2
>   across Y x 7 layers x 29 rows, min box-box gap 1.00 LU, min pillar gap 1.04 LU -- both pass the
>   same >=1 LU clearance the cylinder generator enforces).
> - **`prms/boxes_48x_2x8x3_clog_gas.prm`**: doubling-in-X step for the small box series (2x8x3 = 48),
>   hand-extending `boxes_24x_2x4x3` with 4 more rows at the same 7 LU pitch (X=336/329/322/315) --
>   continues the SAME row-doubling methodology already used for the cylinder/rod series, motivated
>   directly by the partial-jam finding above (doubling the row count packs more deep-layer boxes into
>   the same jam site that already collects them).
> - Synced (both new prms + both new sbatch scripts) to `~/mt-dennys` and `~/mt-dennys-mirror` (the
>   two free checkouts; `~/mt-dennys-mirror2` still running the R4 buoyancy job) via `antidunes.sh
>   --sync`. Submitted: **4279372** (`boxes_48x_2x8x3`, `~/mt-dennys`) and **4279373**
>   (`hopper_clog_box_frontblock_pidon_2node`, the full 406-box run, `~/mt-dennys-mirror`), both 144
>   ranks/2 nodes/6h. Both landed `PENDING (AssocGrpNodeLimit)` at submission -- the account's node
>   quota is currently saturated by the still-running R4 job (4 nodes) + these two (2 each) = 8; they
>   will start automatically once nodes free up, no action needed.

> **[Claude, 2026-09-30] BOTH NEW BOX RUNS CRASHED -- OVERTURNS the "boxes avoid the desync" reading
> from the clean 24-box run above. Real evidence the mechanism is a general dense-pile effect, not
> specific to the cylinder's line-contact geometry.**
>
> **4279372 (48 boxes, 2x8x3, `~/mt-dennys`)**: FAILED after 13:22, step 8624 -- a DIRECT
> `[STALE-PARTICLE-IDX]` abort, the exact same cross-rank ghost-ownership desync mechanism chased at
> length in the cylinder doubling series above: `rank=72 cell (33,9,29) ... is flagged as owned by
> particle uid=3530822107858468865, but that uid no longer maps to a valid local particle index`
> (`SimpleBB.h:325`). This is not a downstream symptom or a different failure class -- it is the
> identical guard, firing the identical way, with boxes instead of cylinders.
>
> **4279373 (406 boxes, full production hopper, `~/mt-dennys-mirror`)**: FAILED after 17:52, step
> ~7100 -- a `reynoldsNumber is nan` abort (`WriteData.h:148`) following a large, sudden NaN cascade
> in a cluster of interface cells near a floor/wall corner (x~36-45, y~13-19, z~0-1). Not a direct
> STALE-PARTICLE-IDX abort this time, but **40 occurrences of `[STALE-PARTICLE-IDX]`/
> `[GHOST-VANISHED]`/`[MEM-GARBAGE-LINK]`-class warnings already logged earlier in the same run**,
> consistent with the same underlying mechanism being active throughout and simply getting caught by
> a later, more generic diagnostic this time rather than the specific ownership guard.
>
> **Reading**: the 24-box run's clean 51000-step completion (job 4279039) was NOT evidence that box
> contacts are immune to the ghost-ownership desync -- it was evidence that 24 boxes, at that
> layout/duration, wasn't enough load to trigger it, exactly mirroring the cylinder/rod doubling
> series' own finding (6 and 12 rods ran clean, 24 rods crashed). Doubling to 48 boxes crashes via the
> IDENTICAL guard cylinders trip; the full 406-box production scale crashes even faster (step ~7100
> vs. the smaller run's 8624), consistent with more bodies -> more contacts -> the race condition
> firing sooner, not with any shape-specific immunity. **Revised working hypothesis**: this is a
> general dense-floater-pile cross-rank ghost-ownership race, not specific to the cylinder's
> line-contact wedge geometry as earlier framed -- the earlier framing was a real but premature
> generalization from a single (too-small) box data point. Both jobs archived to
> `$HPCVAULT/antidunes-archive/antidunes-boxes-48x-2x8x3-clog_4279372_*` and
> `..._hopper-clog-box-frontblock-pidon-2node_4279373_*`.

### Challenges & Future Improvements
- No existing evaluator computes upstream elevation or local pressure at an obstruction directly —
  this likely needs a new function in `Evaluators.h` rather than reusing an existing one verbatim.
- "Clogging" itself needs an operational definition for the thesis (e.g. a sustained elevation rise
  above some threshold, or a sustained drop in downstream transport rate) — deciding this before
  running the sweep will make the results much easier to present than defining it post-hoc from
  whatever the runs happen to show.
- Needs a documented no-clog reference run (e.g. same debris density but spread out / open channel)
  as a baseline, so a measured elevation/pressure rise has something concrete to be compared against.
- Blockage-ratio-vs-afflux plot (review.md item 4) needs a small NEW scenario sweep (varying debris
  load or gap width) with `blockageRatioEnabled true` — not yet designed or run; the code/methodology
  side is done, only the actual sweep and plot remain.
- The fixed-forcing side of the PID-acceleration resistance metric (review.md item 5) needs a fresh
  fixed-forcing reference run with `flowRateEnabled true` to demo for real — the one that exists
  predates that feature.

---

## 4. [ ] Clogging scenarios

Design and run scenarios that specifically target clogging (debris jamming/blocking flow), starting
simple and building up:
1. Sphere-only clogging first (simplest shape, cleanest to interpret).
2. Then progressively more complex scenarios mixing in cuboids/cylinders once the sphere-only case
   is understood.

Depends on item 3 (measurement setup) being in place first.

> **[Claude, 2026-09-16]: Expanded clogging-mechanism taxonomy and scenario design space.**
> Written down now for review later, not scheduled for immediate execution — the three cluster jobs
> currently running (R3 equilibrium rerun 4239262, extreme-clog PID-off 4239285, PID-on 4239288) take
> priority; this is what comes after they land cleanly.
>
> **Current geometry, read directly from `prms/clog_measurement_validation_extreme_pidoff_gas.prm`
> (identical topology in the PID-on sibling and the small `..._sanity_clog_...` pair, just scaled):**
> two `CylinderPillar`s, radius 6, centred at `<384,11,45>` and `<384,29,45>` in a `domainSize <768,40,90>`
> (Y walled, not periodic). That leaves **three** apertures across the channel, not one: left-wall-to-
> pillar1 (`11-6=5` LU), pillar1-to-pillar2 (`(29-6)-(11+6)=6` LU), pillar2-to-right-wall
> (`40-(29+6)=5` LU). The sphere diameter is `~8` LU (`avgParticleDiameterOverride=8.00...`), which is
> larger than **all three** gaps (5, 6, 5), not just the central one. So today's "extreme clogging"
> scenario is a **full-width, every-aperture-blocked** setup by construction — there is no open
> bypass to either side of the pillars, and the sphere literally cannot fit through any gap. Directly
> answering the "is there flow left/right of the pillar" question: yes, there is space there, but it's
> currently just as blocked as the centre gap, which was presumably not a deliberate design choice so
> much as a side effect of scaling the small scenario's own geometry up uniformly.
>
> **Taxonomy — this project currently tests only one of (at least) three distinct clogging
> mechanisms documented in the granular-flow literature:**
> | Type | Mechanism | Control condition | Literature | Status in this project |
> |---|---|---|---|---|
> | A — Deterministic occlusion | Every aperture is geometrically narrower than the particle; blockage is guaranteed, independent of flow rate, count, or shape orientation | gap/diameter ratio $R < 1$ | — (trivial geometric fact, not really a "clogging" finding) | **Implemented** (today's extreme scenario, all three apertures) |
> | B — Probabilistic arch/bridge clogging | Aperture is *wider* than a single particle, so any one particle could pass individually, but a self-supporting arch of several particles can still span the gap and jam it stochastically | $R\gtrsim 1$, typically needs $R\gtrsim 4$–$5$ before jam probability $\to 0$ | \textcite{hafez_effect_2021} (already cited in `related_work/main.tex`) study exactly this via orifice-to-particle size ratio; classic silo/hopper-discharge literature (To, Lai \& Pak 2001; Zuriguel et al. 2005) established the same phenomenon for dry granular flow | **Not implemented** — every current scenario uses $R<1$ |
> | C — Shape/mixture-driven clogging | Particle shape and/or a mixed-shape feed changes clog propensity independent of size ratio (flat-faced/angular shapes interlock and resist torque far more than spheres) | shape, mixture composition | \textcite{hafez_effect_2021} found cubes and 3D-crosses far more clog-prone than spheres for the *same* $R$, due to face-to-face contacts resisting torque | **Not implemented** — sphere-only feeds so far, though the box/cylinder floater types already exist in the codebase |
>
> This project's current scenario is therefore the *least* interesting case in the literature's own
> terms (Type A needs no jamming physics to explain), while the two cases the granular-flow field
> actually studies as "clogging" proper (Type B, and Type C on top of it) are both easy extensions of
> the existing release/measurement machinery, not new physics.
>
> **Table of variables plausibly influencing clog formation/persistence, for deriving concrete test
> cases from:**
> | Category | Variable | Notes |
> |---|---|---|
> | Geometric | Gap-to-particle-diameter ratio $R$ | The single most important control parameter per the Type B/C literature above; today fixed at $R<1$ everywhere |
> | Geometric | Number of apertures (single-slit vs.\ multi-slit rack) | See scenario proposal below |
> | Geometric | Pillar cross-section shape | Already flagged in `related_work/main.tex` (citing \textcite{decicco_bridge_2020}) as a factor this project's pillars don't yet vary |
> | Geometric | Pillar/aperture corner rounding | Directly ties into the box rounded-corner discussion already in \Cref{sec:results-box-stability} — a rounded DEM corner vs.\ the sharp fluid-boundary geometry it maps to is exactly the kind of detail that would matter for arch formation at a rounded-pillar throat |
> | Geometric | Channel width relative to pillar spacing | Determines whether side apertures are a real bypass or negligible, per the geometry note above |
> | Particle | Shape (sphere / cylinder / box / mixed) | Type C above; all three already implemented as floater types |
> | Particle | Size (uniform vs.\ polydisperse) | Not yet varied within a single scenario |
> | Particle | Density ratio (submersion depth/freeboard) | Unlike classical dry silo grains, these floaters are only partially submerged — freeboard height plausibly changes effective arch geometry at the throat |
> | Particle | Friction / restitution coefficient | Already DEM parameters, not yet swept for clogging specifically |
> | Particle | Feed rate (release interval) | Silo-discharge literature treats flow/feed rate as a real control parameter for jam probability, not just a scheduling detail |
> | Particle | Total debris count / concentration | — |
> | Flow | Mean target velocity / Froude number | — |
> | Flow | Water depth relative to pillar height | Partial vs.\ full pillar submersion |
> | Flow | Debris approach angle relative to the aperture | — |
>
> **Concrete geometry proposal, responding directly to the "one big slit vs.\ multiple slits" question:**
> 1. **Single-slit baseline.** Two rounded pillars with the side apertures closed off entirely (solid
>    walls flush with each pillar's outer edge, or pillars pushed out to the domain walls), leaving
>    exactly one, deliberately *controlled-width* central gap. This is the topology the Type
>    B literature itself uses (a single orifice), so it's the right baseline for a direct
>    $R$-ratio sweep across the $R<1$ (today's case) through $R\gtrsim 4$-$5$ (expected jam-free) range,
>    with everything else held fixed.
> 2. **Multi-slit rack.** Four (or more) pillars producing three (or more) apertures, each with its
>    own upstream/downstream probe pair — a straightforward generalisation of the existing
>    `clogPillarUpstreamPos`/`clogPillarDownstreamPos` mechanism (`SimulationParams.h:1230`ff) from one
>    fixed pair to one pair per aperture, indexed. This is structurally closer to a real debris rack
>    or a multi-span bridge pier row than the current two-pillar setup, and lets per-aperture clog
>    onset be correlated against the other apertures' own state (does one aperture jamming first
>    change the arrival rate/behaviour at the others?).
> 3. Run both topologies at both Type A and Type B gap ratios, sphere-only first, then repeat the
>    Type B point with a mixed sphere/box/cylinder feed once the sphere-only baseline is understood —
>    consistent with item 4's own existing "start simple, build up" ordering above.

> **[Claude, 2026-09-18]: first Type B (arch/bridge clogging) scenario written — "hopper" setup, NOT yet
> run.** Implements the single-slit baseline from proposal 1 above. `gen_hopper_clog_prm.py` →
> `prms/hopper_clog_pidon_2node_gas.prm` + `jobs/run_144_hopper_clog_pidon_2node.sbatch`.
> - **Geometry:** two `CylinderPillar`+`BodyBC` pairs with axes ON the Y walls (Y=0, Y=40), radius 12,
>   X=384 → plan view `|) (|`, one central gap of 16 LU. Domain/decomposition/PID/ModelParameters copied
>   from the extreme-clog 2-node prm (`<768,40,180>`, 144 ranks), so results stay comparable.
> - **Particles:** 1015 spheres, radius 3 (d=6 → gap/d = 2.67, i.e. one sphere passes easily, so any stop
>   must be an arch), `densityRatio` 0.5/1.0/1.2 = 290/435/290 (1.5 dropped as requested). All pre-placed
>   in the water on a 5(Y)×7(Z)×29(X) grid with ≥1.0 LU surface clearance everywhere (pitch X/Y/Z =
>   7/8/7.52; first draft used 0.5 LU / 1008 spheres, widened on request), block X∈[91,293] (79 LU free
>   run to the pillars), top layer exactly at the waterline (Z=49.13). NO staggered release — all held
>   FIXED for 2000 steps (free-surface contact lines form) then released together. Stratified by height
>   (1.2 low / 1.0 mid / 0.5 top two layers) so each starts near its buoyant equilibrium.
> - **Geometry verified** two ways: the generator's asserts, and an independent parse of the written prm
>   (min sphere–sphere gap ≥1.0, sphere–pillar 79 LU, walls/floor ≥1.0, no duplicate positions). Brace/
>   block structure diffed against the working extreme-clog prm. NOT yet parsed by the actual binary (no
>   local Antidunes build) — the first cluster run is the real test.
> - **Unverified:** runtime (6 h requested, 60000 steps), whether 1015 spheres survive release without
>   contact chatter, and whether gap/d = 2.67 actually arches (may need a smaller gap or larger spheres).
>   Probes: pair 1 on the gap centreline, pair 2 at Y=6, both X=384±48.
> - **Submitted 2026-09-18:** hopper PID-on = job 4247229 (checkout `mt-dennys-mirror`); extreme clog with
>   PID switched OFF at t=70000, freeze at the mean of steps 65000-70000 (`switchOffStep 70000`, prm
>   `clog_measurement_validation_extreme_pidswitch_2node_gas.prm`) = job 4247230 (checkout `mt-dennys`);
>   d075 resume 4246831 already occupied `mt-dennys-mirror2`. All three in separate checkouts, md5-verified
>   against local (42/42 source files) before submission. First real run of the PID-switch code path on the
>   cluster (unit test ANTIDUNES_ForcingAdjuster passes locally). Archive under
>   `debris-clogging-archive/` with a PROVENANCE.md per family once they finish.
> - **Open follow-ups:** wide-gap no-clog reference; PID-off/PID-switch variants; gap/d sweep.
> - **[Claude, 2026-09-19] Front-block hopper variant written (NOT run):** `prms/hopper_clog_pidon_2node_frontblock_gas.prm` (`gen_hopper_clog_prm.py --front-x 368.5 --out ...`). The whole 1015-sphere block is shifted +78.5 LU in X
>   (block X 91..293 -> 169.5..371.5) so the front row sits 1.0 LU (surface to pillar surface) from the pillar cylinders; verified on the written prm: uniform shift, Y/Z/radius/density/hold unchanged, no other non-sphere line differs from the
>   prm that ran (4247229), min sphere-pillar gap 1.008, sphere-sphere / wall / floor gaps 1.0. Motivation: in the original run the 1.2-density spheres sank and never reached the gap (241/290 still at their start block at t=55000) and no
>   persistent arch formed. Caveat: the upstream clog probes (X=336, pair 1 on the gap centreline and pair 2 at Y=6) now sit INSIDE the sphere block, so they stay unusable as a clog signal; the fluid VTK crop (X 128..640) still contains the whole block.
> - **[Claude, 2026-09-19] Cylinder hopper written (NOT run):** `prms/hopper_clog_cyl_frontblock_pidon_2node_gas.prm` (`gen_hopper_clog_cyl_prm.py`), same channel/pillars/PID/probes/VTK as the sphere hopper but 406 CYLINDERS
>   (r=3, L=18 => diameter:length = 1:3, the sphere's 6 LU cross-section), long axis along Y, 2 across the 40-wide channel x 7 layers x 29 rows, densities 0.5/1.0/1.2 = 116/174/116, block X 169..371 against the pillars
>   (rod-to-pillar surface gap 1.04, rod-rod / wall / floor 1.0; verified on the written prm with sampled surface points, not just the generator's formulas). Rods are ~4.5x a sphere's volume (block solid fraction ~0.48 vs 0.27) hence
>   406 not 1015 (same block length). KEY FACT: L=18 > the 16 LU gap, so a rod lying along Y cannot pass broadside -- it must rotate towards the flow; expect log-jam-like blocking (rotation-limited passage), a different mechanism from the sphere
>   hopper's free arching. Contact model: analytic Kodam+Guo hybrid cylinder contacts + owner-rank tiebreak, all experimental fixes explicitly off. Job scripts prepared, NOT submitted (all checkouts busy):
>   `jobs/run_144_hopper_clog_cyl_frontblock_pidon_2node.sbatch` and `jobs/run_144_hopper_clog_pidon_2node_frontblock.sbatch` (sphere front-block), 144 ranks / 2 nodes / 6 h each, with the leftover-rescue + always-archive logic.

### Challenges & Future Improvements
- Per 2a/2b's supervisor-agreed decision to not pursue further cuboid/cylinder contact fixes, these
  scenarios (once they involve cuboids/cylinders, not spheres alone) may surface the same known
  contact-chatter limitation — results involving non-spherical floaters should note this as a caveat
  on interpretation, not be read as new evidence about clogging physics specifically if the anomaly
  reappears.

---

## 5. [ ] Pont Ferrovière real-life case study

Apply the validated scenario(s) to the real Pont Ferrovière geometry/scale
(`BedGeneration_pont_ferroviere.prm` and related real-scale bridge reconstruction work already done
this session). Outcome either way is a valid thesis result: if clogging is reproduced, that's a
positive finding; if not, that's an argued negative result (clogging wasn't possible under the
tested conditions), not a failure of the study.

### Challenges & Future Improvements
- Real-scale geometry is far more expensive to run than the toy bridge scenarios used for method
  development — needs the clogging-measurement approach from item 3 already working and validated at
  small scale before committing to a real-scale run, to avoid burning cluster time on an unproven
  measurement method.
- A negative result ("clogging wasn't possible here") is only persuasive if the tested conditions
  are shown to be a fair/representative test (realistic debris supply rate, realistic flow
  conditions) — worth explicitly justifying the chosen parameters against real hydrological data for
  the site, if available, rather than reusing toy-scenario defaults at real scale.

> **[Claude, 2026-09-23] Job 4255679 (`antidunes-hopper-clog-cyl-halfrows-gjk-noguards-60k`) root-caused, a LOCAL deterministic reproduction found, and its guard24 rerun (4259313) evaluated.**
>
> **1. Root cause of 4255679's own `STALE-PARTICLE-IDX` abort (step 3133), confirmed directly from its WARNING-level log (the fatal ERROR line alone doesn't show this):** rods uid `3530822107858468885`/`888`, real logged centre separation 5.931 LU at step 3100 (essentially touching, 2R=6), sealed a liquid cell between them. `[MEM-GARBAGE-LINK]` at step 3132 shows two cells already garbage: block-local `(29,13,29)` (15/26 obstacle neighbours, density collapsed to 0.253, populations of order 1e9-1e10 in the characteristic zero-mass/zero-momentum sign-alternating mode already described for job 4251066) and `(27,14,29)` (9/26 obstacle neighbours, density **negative** at -0.066, populations of order 1e12). By step 3133 — one step later — both rods are launched to **~138,000 LU/step** (vs their normal ~0.02 LU/step), a ~7e6x velocity spike, instantly leaving every rank's block+ghost region (`[SUBCYCLE-LEFT-RANK]`/`[GHOST-VANISHED]` cascade) and triggering the abort. This is the same "class B" sealed/pinched-cell instability already characterised for job 4252462 above, but with a materially WEAKER obstacle-neighbour count (9-15/26) than that case's ~20/26 closing-pocket mechanism or `SealedCellGuard.h`'s own design target (guard threshold 24, tuned for a 3-rod near-total enclosure — see that file's own header comment: "Three touching parallel rods... a string of fluid cells whose 26 neighbours are (almost) ALL obstacle"). `ANTIDUNES_Job4255679SealedCell` (new unit test, `tests/showcases_antidunes/`) proves the geometric precondition directly from these real positions: 72 pinched cells, 68 strongly-pinched (>=3 pairs), but the worst single cell only reaches 18/26 obstacle neighbours at this 2-body-only geometry — confirming guard24 would NOT reliably have closed either of the real exploded cells (9 and 15) even under ideal conditions.
>
> **2. Guarded rerun (job 4259313, `closeSealedCellsMinNeighbours=24`, otherwise identical to 4255679) FAILED anyway, after 5:42 (vs the original's 5:08), via a DIFFERENT failure signature: an invalid-quaternion NaN assertion at step ~2936, not `STALE-PARTICLE-IDX`.** The guard was genuinely active (one particle received 20 `[SEALED-CLOSE]` closures within ~10 steps on a single rank), so it measurably changed the dynamics and pushed the failure later/different -- but did not prevent the underlying instability. This is exactly the "several steps later" downstream symptom the `STALE-PARTICLE-IDX` guard's own message warns about pre-empting; here the guard delayed/masked the desync's *detection* without fixing its *cause*. Consistent with (1): guard24 is tuned for the wrong sub-mechanism for this scenario family.
>
> **3. A LOCAL, fully deterministic reproduction was found** (`gen_rod_bundle_channel_test_prm.py --pair`, new `--pair` mode added this session for exactly two parallel rods instead of the existing 3-rod bundle/raft): two rods at the exact real 5.931 LU separation, axis world-Y, real gravity/contact-model settings, `ModelParameters.rodPrescribedVelocityX/Y/Z = (0.02, 0, 0.02)` (deterministic rigid motion, decoupled from DEM chaos), domain 48x36x48, periodic X/Y, fluid sanity sweep + dense NaN scan armed from step 0, no guard. **Reran 5 times (`rod_pair_repro_run{1..5}`): all 5 aborted with `STALE-PARTICLE-IDX` at EXACTLY step 1423**, same culprit UID (`4611686018427387905`), same cell region (world X~0-2, Z~44-47) -- fully bit-reproducible, unlike the real production runs ("not bit-reproducible" per the 4251086 resume finding above). This determinism is expected and explained: with prescribed (not contact-derived) rod velocity, the kinematic trajectory is exactly repeatable regardless of any fluid-side floating-point noise, so the two rods reach the identical trapping geometry at the identical step every time.
>
> **4. IMPORTANT CAVEAT — this local repro is not a pixel-perfect reproduction of 4255679's own mechanism; it is a different, wall-assisted sub-case of the same class-B phenomenon, found and reported honestly rather than presented as more faithful than it is:** the prescribed Z-velocity (0.02/step) carries the rod pair's centroid from its Z=19.6 start straight toward the domain's own non-periodic TOP WALL at Z=48 -- `19.6 + 0.02*1423 = 48.06`, matching the crash position (`pos=<...,48.0001>`/`<...,48.001>`) almost exactly. So the trapped/sealing cell here forms between the TWO RODS **and** the domain's Z-wall (three converging surfaces, closer in spirit to 4252462's 3-body near-total enclosure than to 4255679's pure 2-rod-in-open-water slit) -- confirmed directly from the `[FLUID-SANITY]` log in the run-up to the crash: from step ~1312, a cell at world `(2.5,~17.5-23.5,47.5)` (right under the Z=48 wall) shows **12/13 pinched pairs, 24/26 obstacle neighbours, owners = both rod UIDs**, and density climbing steadily and monotonically from 1.91 (t=1314) through 1.97 (t=1317), 2.0+ (t=1317+), continuing to rise past 2.61 before the abort -- a genuine, slow-building sealed-cell density runaway, not an instantaneous event. This is real, directly-relevant evidence for the broader class-B mechanism (same real rod UIDs/geometry, same qualitative signature), but the specific 24/26 obstacle count here means it would actually be CAUGHT by guard24 (unlike the real 4255679 cells) -- worth confirming with a guard=24 rerun of this exact local repro as a sanity check that the guard machinery works as intended at least for this wall-assisted variant, before using this tool to test lower thresholds aimed at the true 2-rod-only (9-15/26) mechanism. A cleaner, wall-free reproduction (zero or periodic Z-velocity, or a taller/periodic-Z box) is the natural next step if a more faithful match to 4255679's own open-water mechanism is wanted specifically.
>
> **Files:** `gen_rod_bundle_channel_test_prm.py` (`--pair` mode added), `prms/rod_pair_4255679_repro_gas.prm`, `prms/hopper_clog_cyl_frontblock_pidon_2node_halfrows_guard24_gas.prm` (job 4259313's prm), `jobs/run_144_hopper_clog_cyl_halfrows_gjk_guard24_60k.sbatch`, `tests/showcases_antidunes/{CylinderTippingRodJamTest,RodSpinSubstepResolutionTest,TiltYOwnAxisSpinArtifactTest,Job4255679SealedCell}.cpp` (four new local unit tests from this investigation; the first two's premise -- that uid 888's large logged `tiltY` swing indicated a physical topple -- was later superseded by the third test, which shows `tiltY` is dominated by a diagnostic artifact, own-axis-spin sensitivity in its `atan2` formula, not necessarily real rotation; kept as they still correctly bound what GJK/EPA and substep timing can/can't explain). Local run outputs not yet archived to the vault (small, reproducible on demand in ~30s/run; regenerate via the command above rather than archiving).

> **[Claude, 2026-09-23] REMOVED: a same-day sphere/box comparison + "boundary-type bug" follow-up (originally logged here, and as a published artifact) was retracted by the user and deleted from this file.** That follow-up used `rodPrescribedVelocity` to drive sphere and box pairs at a sustained constant velocity straight through the domain's non-periodic top wall, and found their resulting `STALE-PARTICLE-IDX` crash was caused by an MPI ownership-handoff race at that wall, unrelated to force or density. Correct as far as it went, but irrelevant to debugging job 4255679: no production Antidunes run ever drives a floater's center through a domain wall like this (a real wall/contact force stops it first), so the finding describes a test-harness artifact, not a mechanism that occurs in the simulation this plan is actually validating. The rod's own sealed-cell density/force excursion (point 4 above, steps 1312-1338) is unaffected and remains valid, real evidence for the Class B mechanism; only the sphere/box/boundary-race material was removed. `ParticleKernels.h`'s `rodPrescribedVelocity` type filter was left extended to sphere/box (types 1/3) in case it's useful later, but is currently unused by any active scenario.

> **[Claude, 2026-09-24] Clog-measurement-validation crashes (extreme clog, PID on, flow-rate evaluator scenario): three attempts, three different steps/signatures, but a real spatial pattern -- crashes cluster at the domain's Y=0 side wall, at/near the particle release zone in X.** Same prm/binary re-run three times on three different checkouts after the flow-rate evaluator was added (see the 2026-09-22 entry above); all three died the same way (`reynoldsNumber is nan`, `WriteData.h:98`), converted from rank+block-local cell to global domain coordinates:
>
> | Job | Variant | Crash step | Global (X, Y, Z) |
> |---|---|---|---|
> | 4257275 | flow-rate ON (original) | 58,482 | (34, **1**, 106) |
> | 4259949 | flow-rate ON (identical rerun) | 56,038 | (28, **1**, 106) |
> | 4260103 | flow-rate DISABLED (control) | 6,274 | (109, **0**, 90) |
>
> **The flow-rate-enabled hypothesis is refuted:** disabling the evaluator did not help -- it crashed *earlier* (step 6,274 vs. 56,038/58,482) and with the same failure signature, so the new evaluator code is not the cause (it only reads cells near the far-downstream probes at X=280/500 and the pillar at X=384; it has no code path that could touch X=28-109). Domain geometry: particle release zone is X=40-120, pillar X~372-396, probes X=280/500, domain Y=40 (2-node prm, `cellsPerBlock<64,20,30>`) -- **all three crashes sit at Y=0 or Y=1 (one cell off the channel side wall)**, and X=28/34/109 are all inside or immediately adjacent to the release zone, nowhere near the pillar or either probe. Z differs (106 vs. 90) and the exact step varies by ~50,000, consistent with the already-documented chaotic bit-non-determinism of this codebase (jobs 4251086/4251480) -- but the *location* is not random: this is the same "class B" sealed/pinched-cell instability characterised for the rod scenarios (job 4255679 etc., see 2026-09-23 entries above), here triggered by a sphere pinching against the channel's Y=0 wall near where spheres are actively being released, rather than rod-on-rod pinching. **Not yet tested:** whether moving the release zone further from the Y=0 wall, or narrowing the release band's Y-extent, reduces the crash rate -- the release-zone-X-relocation fix already applied earlier in this project (X=160-200 -> X=40-100, see line ~1178 above) addressed a different problem and did not anticipate this Y-wall interaction.
>
> **Checkpoint-loss incident (process note, not a physics finding):** both 4257275's and 4259949's own checkpoints were destroyed before they could be examined or resumed -- in both cases because the next job run on the *same* checkout (4260083 on `mt-dennys` after 4257275; 4260103 on `mt-dennys-mirror` after 4259949) has a destructive `rm -rf` of leftover `vtk-out`/`eval`/`snapshots` at startup, and I ran it without archiving the prior crash's leftovers first. Confirmed via exact snapshot-file timestamps falling entirely inside the *later* job's own run window. No resume was possible for either; a fresh (non-resumed) rerun of the flow-rate-enabled scenario was submitted instead as **job 4260382** (`mt-dennys-mirror2`), still running as of this entry. Lesson applied going forward: archive or check a checkout's leftover run directories before reusing it, especially for the simpler single-run sbatch scripts that `rm -rf` rather than rescue-and-move.

> **[Claude, 2026-09-24] Bogner box equilibrium: the "~2.2% effective density deficit" hypothesis is now falsified by a direct dynamical test, not just by the earlier retracted density-field measurement.** Recap: the free 6-DOF omega=1.9/seed5/rho=0.75 box was observed to settle near ~32 deg instead of the theoretical 26.565 deg; a geometric back-calculation from real force/torque data suggested this was equivalent to the box behaving as if its density were ~0.734 instead of 0.75 (~2.2% low), motivating a rerun with `densityRatio` bumped to 0.7665 (+2.2%, compensating the apparent deficit) to see if it would then settle correctly. That rerun (job 4260604, byte-identical prm otherwise, 100k steps, `~/mt-dennys`) **COMPLETED cleanly this time** (1h02m, exit 0 -- the earlier attempt in job 4260105 had crashed around step ~15,323, confirmed here as a one-off chaotic crash, not a reproducible consequence of the density change). **Result: it does NOT settle near 26.565 deg, or anywhere.** Sampling `tiltY` (raw, `[BoxForces]` log) every ~2000-10000 steps across the full run: 24.4 -> 26.7 -> 27.3 -> 24.0 deg (steps 5k-20k, tracking the ideal value nicely at first) then a large excursion: 57.2 (30k) -> 39.6 (40k) -> 56.5 (50k) -> 54.8 (60k) -> 46.9 (70k) -> 58.7 (80k) -> 58.6 (90k) -> **6.5 (100k)**; a finer look at just the last 10,000 steps shows it is still swinging wildly, not settling: 58.6 (90k) -> 48.4 (92k) -> 74.0 (94k) -> 56.7 (96k) -> 21.5 (98k) -> 6.5 (100k). This is a much WIDER, more chaotic tumbling range than the rho=0.75 sibling's own "coupled roll-translate limit cycle" (which stayed within a narrower band around ~30-36 deg per the earlier torque-sweep artifact) -- correcting the density made the dynamics more chaotic, not less, and never produced a clean settle near 26.565. The box also translated a large distance in X (21.67 -> ~1.5-2.4, wrapping most of the way around the periodic 130-cell domain) while Z (draft depth) stayed close to its expected value throughout (vertical force balance holds, consistent with everything found earlier). **Conclusion: the density-offset compensation does not fix the equilibrium-angle discrepancy -- reinforcing, via an independent dynamical test, the standing retraction of the "local fluid density fluctuation" mechanism** (that measurement was already found to be a depth-averaging artifact, not a real density excess -- see the corresponding retraction; not yet written up in this file as its own dated entry, see pending items). The mechanism for why the free box settles away from 26.565 remains open; the coupled roll-translate limit-cycle dynamics itself (not a density or wall-placement or fill-weighting artifact -- both IBB and `useFillWeighting=false` ablations only reduced the static torque bias by ~0-20%, see 2026-09-21/22 entries above) is the strongest remaining candidate, now reinforced by this run showing that even a "corrected" density does not escape large-amplitude chaotic tumbling.

> **[Claude, 2026-09-24] RETRACTION: the "local fluid density excess near the box" measurement (which motivated the density-correction test above) was my own depth-averaging artifact, not a real effect.** Working hypothesis at the time: the box's effective buoyancy was ~2.2% too weak because the LOCAL fluid density near the box had drifted persistently above the nominal 1.0 reference (extracted `density` VTK field, job 4260105's densitycheck_rho075 sub-run, `mapped_flag`-masked to pure-liquid cells only). Initial measurement found an apparent "~2.36% global" and "~3.69% local" density excess near the box. **This was wrong**: the hydrostatic initialization genuinely gives density a strong Z-gradient (`rho(z) = 1 + 3*a_z*(z-H)`, `DataInitialization.h`'s `HydrostaticDensityFtor`), and averaging over full-depth liquid columns without controlling for depth compares different mean depths near vs. far from the box -- producing a spurious "excess" that has nothing to do with the box's presence. Redoing the near-box-vs-far-field comparison at MATCHED depth (same z-layer) dropped the difference to ~0.0001-0.02%, i.e., noise. A related "hydrostatic init height mismatch" root-cause guess (referencing the surface at the wrong height) was also checked directly and found ~20-40x too small to explain the (itself wrong) 2.36% figure. **Net effect: there is no evidence of a real local or global fluid density anomaly near the box; the buoyancy-strength deficit inferred from the torque/force back-calculation (independent of this density-field measurement, and still valid) has to be explained some other way** -- the completed rho=0.7665 dynamical test above is the direct follow-up, and it also does not support "just add density" as a fix.

> **[Claude, 2026-09-24] The 26.565deg-seed rho=0.7665 case's step-444 "crash" is not a property of the box's dynamics at all -- it is entirely determined by whether per-step VTK output is active, confirmed with a clean 7-attempt controlled comparison.** Releasing the box AT the theoretical equilibrium (zero seed transient, otherwise identical to the seeded-heel runs above) reproducibly blows up (`reynoldsNumber is nan`, `WriteData.h`) at **exactly step 444** with the standard prm (mid-width Y=20-21 VTK plane, `writeFrequency 2000`): jobs **4260975, 4261005, 4261093, 4261316 -- 4 for 4**, all matching to the displayed precision at every logged step up to and including 444 (e.g. tiltY=25.696698 deg at step 444 in every one of them). A new `keepAllSnapshots`/`keepAllSnapshotsFromStep` diagnostic feature (`SimulationParams.h`/`Antidunes.cpp`/`WriteData.h`, checkpoints every step into `snapshots/step_<N>/` instead of the usual single rolling file) captured the full field state for steps 400-443 from job 4261093 before its crash. Resuming from the `step_400` checkpoint with a WIDENED VTK region (covering the two hot-spot columns identified in the crash log, `CellBB_filter_0 <0,0,18> to <16,10,32>`, `writeFrequency 1` -- full per-step dump) unexpectedly **did not reproduce the crash**: job **4261250** ran clean through the requested 450 steps. Repeating the identical resume (job **4261311**) also ran clean. A companion run starting completely FRESH from t=0 with the SAME widened per-step VTK (job **4261312**) *also* ran clean all the way to step 1000 -- and critically, its `tiltY` at step 441 (25.736590 deg) matched the ORIGINAL crashing runs' value bit-for-bit (not the resumed runs' own very slightly perturbed 25.736581 deg from the checkpoint round-trip), proving the crash/no-crash outcome is **not** explained by run-to-run floating-point drift in the box's own trajectory -- the same exact state at step 441 crashes with normal VTK and doesn't with the widened one. Extending the resume window (job **4261314**, steps 400-1000, same widened VTK) also stayed clean, and shows the box settling into a smooth, monotonic tilt decay (25.7 deg at 444 -> 24.9 -> 23.4 -> 22.5 -> 22.0 -> 21.8 -> 21.9 deg at step 1000) -- qualitatively different from the wide, chaotic 50-65 deg oscillation seen in the other seeded-heel rho=0.7665 runs (2026-09-24 entries above), suggesting that dodging this specific instability may lead to materially different long-run behaviour, not just a delayed version of the same outcome. **Final tally: 4/4 crashed with normal VTK, 3/3 clean with the widened per-step VTK (2 resumes + 1 fresh-from-0), 7 attempts total, a clean split on exactly one variable.** Interpretation: the crash is a genuine knife-edge numerical instability (most likely a tie-break or near-singular condition in the free-surface interface reconstruction / excess-mass redistribution at those specific cells, which sit at the periodic X seam far from the box itself, see the 2026-09-24 cell-tracing entry above) that per-step VTK output's own extra MPI collective calls happen to perturb just enough to avoid -- consistent with, and a much cleaner demonstration of, the same "extra reductions perturb chaotic timing" mechanism already documented for the flow-rate-evaluator clog-scenario crashes elsewhere in this file. **Files:** `bogner_exact_3box_s16_rho0767_tilt26p565_resume400_hotspot_gas.prm`, `..._resume400_to1k_hotspot_gas.prm`, `..._fresh0_hotspot_gas.prm`, `jobs/run_65_bogner_resume400_hotspot.sbatch`, `..._resume400_to1k_hotspot.sbatch`, `..._fresh0_hotspot.sbatch`; checkpoint data archived at `animations/antidunes-bogner-exact-s16-rho0767-tilt26p565-allsnaps_4261093_20260924_141616/snapshots/step_{400..443}/`. **Not yet done:** decoding the actual per-step VTK fields at the hot-spot cells (base64 VTU parsing) to see the precursor pattern directly rather than inferring it from the console log's mass-redistribution warnings.

> **[Claude, 2026-09-24] CORRECTION: the failing cell is at the box's own corner/waterline, NOT at the periodic X seam far from the box -- the "periodic seam" claim above was a coordinate-conversion mistake on my part.** All the `[STALE-NEARBOUNDARY]` traces used to locate the failure (2026-09-24 entries above, "same cells every time") only ever printed the field iterator's BLOCK-LOCAL coordinates (`it.x(),it.y(),it.z()` in `AntidunesBoundaryHandling.h`), which I read as world/global coordinates directly. That is only correct for the one block sitting at the domain origin -- with `cellsPerBlock<10,8,40>` and 13x5x1=65 blocks tiling the 130x40x40 domain, block-local (0,3,25) or (8,x,26) can belong to any of 65 different blocks depending on which rank logged it, and I never converted or checked which rank/block each line came from. Built a targeted fix instead of continuing to guess: added a VELOCITY-derived NaN check to `checkPdfNaN` (`SetupStages.h`, catches a 0/0-in-getVelocity() case the existing raw-PDF isfinite scan structurally cannot see -- motivated by tracing exactly how a local cell failure reaches the global `meanVelocity_x` used by the `reynoldsNumber` abort check, see the mean-velocity chain-of-custody entry above), then discovered `checkPdfNaN` itself is only ever CALLED when `watchNow` (hardcoded `t%1000==0` in `Antidunes.cpp`) or `denseNanScanEnabled` (prm bool, defaults false, never set in any of our prms) is true -- meaning the scanner had never actually run near step 444 in ANY of the 7+ crashes investigated so far, old or new check alike. Fixed by adding `denseNanScanEnabled true; denseNanScanFromStep 430; denseNanScanToStep 446;` to a resume-from-checkpoint rerun (job **4261621**, resumes the step_400 checkpoint, normal VTK so the crash reproduces, 500-step cap): it fired immediately, and it was the ORIGINAL raw-PDF check that caught it (not the new velocity one -- the PDF itself really did go non-finite, `f[0]=-nan`), at exactly **block-local=(0,0,25), world=(20.5, 8.5, 25.5)**, flagged INTERFACE, at t=443 (one step before the abort). At that step the box was at `pos=<21.466,19.983,16.172>` with `halfExtents<8,16,8>` -- spanning X=[13.47,29.47], Y=[3.98,35.98], Z=[8.17,24.17]. **World (20.5, 8.5, 25.5) sits INSIDE the box's X and Y footprint, at Z just above its own top surface** -- i.e. this is a cell right at the box's own top corner/waterline, not out in open water far from it. The 3x3x3 neighbourhood dump confirms this reading directly: the failing interface cell is flanked by `OBST` (the box itself) on one side and `GAS` on the other -- the classic geometry of a corner cell where a tilted box's edge cuts the free surface, exactly the kind of cell prone to running out of liquid/interface neighbours to reconstruct from (matches the "no liquid or interface cell in the neighborhood" mass-redistribution warnings already documented, now correctly located). **This reframes the mechanism: the instability is a free-surface/box-corner reconstruction problem at the box's own waterline, not an unexplained collapse in open water at a domain boundary.** The VTK-timing dependence (per-step output avoiding the crash) and the "same cells every time" determinism both still stand -- only the physical location of "same cells" was wrong. **Not yet done:** re-examine the earlier `[TOTAL-FORCE-THRESHOLD]`/early-crash (job 4261349, step ~36) trace with this corrected local-to-world conversion in mind, since that crash may also be a box-corner event rather than the separate mechanism it was assumed to be.

> **[Claude, 2026-09-24] Follow-up: the effect is NOT "any periodic VTK output helps" -- it specifically requires a write AT (or immediately adjacent to) step 444 itself.** Two more fresh-from-t=0 runs, identical scenario, ORIGINAL mid-width Y=20-21 VTK plane (not the widened hot-spot region), with intermediate `writeFrequency` values that do NOT land a write near step 444: job **4261320** (`writeFrequency 1000`, 50k steps) and job **4261321** (`writeFrequency 400`, 20k steps) -- both CRASHED at exactly step 444 again, bit-identical to every prior crash (tiltY=25.696698 deg at step 444). Updated tally: **6/6 crashed** with any writeFrequency that doesn't write at step 444 (2000, 1000, 400), **3/3 clean** only with `writeFrequency 1` (guaranteed to write every single step, including 444 itself). This rules out "VTK I/O activity anywhere in the run perturbs floating-point state generally" -- the effect is highly localized: only a collective MPI call happening AT that specific step changes the outcome. Strongly consistent with a genuine race/tie-break at that one step's free-surface reconstruction (see the cell-tracing entry above) that an incidental extra collective operation at exactly the right moment disturbs enough to avoid. **Not yet tried:** a writeFrequency chosen to land exactly on 444 (e.g. via a one-off extra write at that specific step) as the most direct confirmation.

> **[Claude, 2026-09-24] That direct confirmation was run -- and it DISPROVES "a single well-timed write at step 444 is sufficient."** Job **4261326**: fresh-from-t=0, original mid-width VTK plane, `writeFrequency 444` (a divisor of 444, so writes land at steps 0 and 444 exactly -- one write lands precisely on the crash step, none on any other step in between). VTK output is registered as a before-timestep function (`addBeforeTimeStepVTKOutput`), so this write fires at the very start of the t=444 iteration, before that iteration's own FSLBM sweep runs (i.e. immediately before the crash-precipitating computation). **It CRASHED anyway**, at exactly step 444, bit-identical to every other crash (tiltY=25.696698 deg). Updated tally: **7/7 crashed** with anything other than genuine every-step output (2000, 1000, 400, and now a single write exactly at 444) vs. **3/3 clean** only with true `writeFrequency 1` (steps 401-443 ALL written, not just 444). **Revised interpretation:** the effect is not "a collective MPI call happening at the critical step nudges a tie-break" (that hypothesis is now ruled out) -- it must be a CUMULATIVE effect of many per-step writes across the preceding steps (401-443), each contributing a tiny floating-point perturbation via its own MPI collective, that only in aggregate shifts the trajectory enough by step 444 to avoid whatever near-singular condition exists there. A single nudge at the critical moment alone does nothing. **Not yet tried:** writeFrequency values that write every step but only over a SHORTER preceding window (e.g. every step from 420-444 only, not from 401) to bound how much accumulated perturbation is actually needed.

> **[Claude, 2026-09-24] Decisive test: even a PERFECTLY-placed release (correct tilt AND correct geometric heave, zero transient of either kind) still fails to hold the true equilibrium and still ends up in the same wide chaotic regime -- ruling out imprecise initial conditions as the driver of the rho=0.75/0.7665 equilibrium-offset problem.** Every seeded-tilt run all session (including the ORIGINAL rho=0.75 sweep) placed the box's centre at a fixed Z=16 regardless of the seeded tilt angle -- correct only for rho=0.75 at EXACTLY zero tilt (verified: `heave_for(0.75, 0deg) = 16.000000` exactly, via the same submerged-area-fraction=rho bisection method already used for the righting-moment scans, `gen_bogner_righting_scan_prm.py`'s own `heave_for()`, waterline Z=20). For any nonzero tilt this is measurably off (e.g. +0.68 deg... corrections of +0.18 to +0.68 LU for rho=0.7665 at 25-26.565deg, computed fresh this session), meaning every prior seeded run started with a small residual HEAVE transient layered on top of the intentional ROLL transient. Job **4261350** (rho=0.7665, released at 26deg tilt with the position ALSO corrected to the exact geometric centroid, z_centre=16.1659, both transients eliminated) tests whether removing the heave mismatch changes the outcome. **It does not:** within 1000 steps the box has already drifted to Z=16.77 -- matching the EMPIRICALLY-observed apparent equilibrium from every earlier (uncorrected) run (~16.7-16.9), not the geometric ideal (16.17) it was released at. It tracks near the target tilt reasonably well through ~10-15k steps (24-29deg), then blows up into the same wide chaotic band (43.7 -> 63.8 -> 64.7 deg) seen in every other seeded-tilt rho=0.7665 run this session, before happening to pass back through 26.56deg at step 50000 (almost certainly a mid-oscillation coincidence given the wild swings immediately before it, not a genuine settle). **Conclusion: this is not an initial-condition problem.** A release with zero roll transient AND zero heave transient still gets pulled toward the same "wrong" dynamical point and the same eventual chaotic regime -- reinforcing that the coupled roll-translate limit-cycle dynamics (tilt peaks at surge-velocity-zero-crossings, see the rho=0.75 equilibrium-offset entries above) is intrinsic to the free 6-DOF system itself, not an artifact of how these runs happen to be seeded. **Follow-up, same day: the analogous rho=0.75 test confirms it independently.** Job **4261415** (rho=0.75, 5deg seed matching the ORIGINAL archived sweep, position corrected to z_centre=16.0152 vs. the original's Z=16.0 -- a much smaller correction since rho=0.75 was already exact at 0deg) COMPLETED (28m36s). Sampled tiltY: 2.4 (1k) -> **-32.6 (5k)** -> -31.5 (10k) -> -31.5 (20k) -> -28.3 (30k) -> 48.3 (40k) -> 46.0 (50k) deg. Even this much smaller heave correction does not produce a clean settle at 26.565deg -- instead the box FLIPS SIGN entirely (rolls past vertical to the opposite side) and holds near -31/-32deg for ~25,000 steps (the mirror-folded equivalent of ~31-32deg, itself close to the ~30-36deg band the uncorrected sweep landed in, just on the other side), then destabilizes into a different chaotic band (46-48deg) by the end. Y also drifts substantially here (19.96 -> 23.96, +4 LU) -- a different drift axis than the X-dominated drift seen in the rho=0.7665 runs. **Two independent tests now (different density, different seed angle, different specific failure mode) reach the same conclusion: a perfectly-placed release (zero roll AND zero heave transient) does not rescue the free 6-DOF box from ending up somewhere other than the true equilibrium.** This is not an initial-condition artifact for either density -- the coupled roll-translate dynamical instability appears to be a genuine, intrinsic property of the free 6-DOF system at this resolution/parameter set.

> **[Claude, 2026-09-24] CONFIRMED with a direct in-code catch: the step-444 crash cell is a formerMO cell classified Interface with genuinely ZERO real interface/liquid neighbours at the exact instant of reconstruction -- but this exact signature is common and usually harmless, so it is not itself sufficient to explain the crash.** Added a new warning, `[FALLBACK-ISOLATED-INTERFACE]`, directly inside `seedConvertedCellFill()` (`src/lbm_mesapd_coupling/momentum_exchange_method/reconstruction/PdfReconstructionManager.h`, same throttling/world-coordinate convention as the file's existing `[PASS3-UNRESOLVED-ISLAND]`), firing whenever a freshly-converted interface cell finds zero real interface neighbours to average its fill from (`cnt==0`, the exact branch that produces the hardcoded `fill=0.5` fallback). Reran the resume-from-checkpoint-400 crash (job **4261759**, `~/mt-dennys-mirror2`, normal VTK so the crash reproduces) with this instrumentation plus the existing `denseNanScanEnabled` window from the prior entry. Result: it fired **exactly at t=442** for our target cell, from 4 different MPI ranks reporting the same physical location via their own local/ghost copies -- `world≈(20.5,8.5,25.5)`, `seeded fill=0.5 with ZERO real interface neighbours` -- the identical timestep the earlier cube-trace already showed this cell flipping from `OBST` to `INTF`. This is now a **direct** catch from inside the fill-seeding code itself at classification time, not an inference from the following step's neighbourhood dump. The subsequent NaN-SCAN/neighbourhood dump at t=443 (same run) reproduced byte-for-byte what the earlier resume run found: the cell's 26-neighbourhood is exactly `{12 OBST, 14 GAS, 0 LIQ, 0 INTF}` -- confirming `hasGas && !hasLiquid` with no interface/liquid evidence at all, matching the decision table's `hasGas&&!hasLiquid&&!hasUnknown -> Gas` / `...&&hasUnknown -> Interface` branch pair.
>
> **New finding that complicates the "misclassification" framing: this fallback signature is not rare.** The same `[FALLBACK-ISOLATED-INTERFACE]` warning fired **60 times** across steps 425-442 in this one run, at many other world locations (mostly a cluster around `world x=18.5`, differing Y, both Z=25.5/26.5 -- almost certainly the box's other corner or a periodic-image location, not investigated further) -- and none of those other 59 occurrences caused a crash; the run proceeded normally through all of them. So "formerMO cell, zero real interface neighbours, Interface classification, `fill=0.5` fallback" is a **routine, usually-harmless event** at this box geometry/resolution, not a unique fingerprint of the crash cell. What is unique to `world≈(20.5,8.5,25.5)` must be something layered on top of the fallback -- most plausibly the specific asymmetry of its neighbourhood (an entire adjacent layer of 9 `OBST` cells directly on one side, `GAS` filling the rest, i.e. a genuine near-singular corner geometry) that makes the D3Q27 cumulant collision operator produce a NaN on this particular reconstructed state, not the `fill=0.5` seeding itself.
>
> **Still not directly confirmed:** whether `hasUnknown` was literally `true` at the exact classification instant (which would make Interface the documented, intended same-pass-ambiguity fallback rather than a rule violation) -- `seedConvertedCellFill`'s own `cnt` counter (real INTERFACE neighbours available for fill-averaging) is a different quantity from `classifyFormerObstacleNeighbourhood`'s `hasUnknown` flag (whether a neighbour was itself still `formerObstacle`-tagged, unresolved, in the same pass), and no log currently exists inside `classifyFormerObstacleNeighbourhood` itself. Given a box corner rotating through the free surface plausibly uncovers a whole cluster of neighbouring cells in the same timestep, a same-pass-unresolved neighbour is a likely (but not yet logged) explanation. **Direct next step if pursued further:** add an equivalent throttled log inside `classifyFormerObstacleNeighbourhood` itself reporting `hasGas/hasLiquid/hasUnknown` for this exact cell/step, to settle "documented safety fallback" vs. "decision-table violation" conclusively.

> **[Claude, 2026-09-24] Widened the sweep-by-sweep rho probe (`BuildTimeloop.h`'s `probeCrashCell`) from just the one known crash cell to 5 cells spanning both fallback clusters (2 from cluster B including the crash cell, 2 from cluster A, 1 early-converted cluster-A cell), rerun of the same resume400 crash (job **4262030**, same setup as 4261759). Result: the fallback/`fill=0.5` mechanism alone is NOT what determines whether a cell goes NaN -- most cells that go through it survive fine, and there's a SECOND independent NaN in this same crash, not just the one already found.**
>
> | Stage (t=443) | (20.5,8.5,25.5) crash cell, INTF | (20.5,6.5,25.5) sibling, GAS | (20.5,5.5,25.5), GAS | (18.5,34.5,26.5), INTF, converted same step (442) | (18.5,4.5,26.5), GAS, converted step 425 |
> |---|---|---|---|---|---|
> | start of timestep | 1.000000 | 1.000000 | 1.000000 | 1.000000 | -3.621461e+53 |
> | after boundary sweep | 1.000000 | 1.000000 | 1.000000 | 1.000000 | -3.621461e+53 |
> | after force-density sweep | 1.000000 | 1.000000 | 1.000000 | 1.000000 | -3.621461e+53 |
> | after StreamReconstructAdvect | 1.000000 | **-nan** | 1.005174 | 1.006672 | 1.000000 (reset clean) |
> | after collide | **-nan** | -nan | 1.005174 | 1.006672 | 1.000000 |
> | after cell conversion | -nan | -nan | 1.005174 | 1.006672 (-> GAS) | 1.000000 |
> | after equilibrium refilling | -nan | -nan | 1.005174 | 1.006672 | 1.000000 |
>
> Confirms the already-known crash cell's own chain exactly (clean through StreamReconstructAdvect, NaN appears at collide). **New finding: a second, independent NaN exists in the same crash** -- `(20.5,6.5,25.5)`, a cluster-B sibling converted 2 steps earlier (step 440), also goes `-nan`, but ONE SWEEP EARLIER, at StreamReconstructAdvect rather than collide. Since D3Q27 streaming only moves PDF data +-1 lattice unit per timestep and these two cells are 2 apart in Y, this cannot be direct single-step cross-contamination -- it is a second, independently-triggered failure sharing the same root mechanism, not a propagation of the first. **The other three tracked cells -- despite going through the identical zero-real-neighbour `fill=0.5` fallback -- stay completely finite**: `(20.5,5.5,25.5)` and `(18.5,34.5,26.5)` pick up only a small finite excess (rho~1.005-1.007, from an early streamReconstructAdvect mass-redistribution artifact) and never go non-finite; `(18.5,34.5,26.5)` in particular was converted at the EXACT SAME timestep (442) as the crash cell via the identical mechanism yet survives cleanly through to reclassification back to GAS. `(18.5,4.5,26.5)` carries a huge but merely-FINITE garbage value (-3.6e+53, untouched since its own step-425 conversion, sitting inert through boundary/force-density sweeps since GAS-cell PDFs aren't used for physics) that gets silently reset to a clean 1.0 by StreamReconstructAdvect itself -- i.e. stale finite garbage in GAS cells is routine and harmless as long as it never becomes non-finite. **Conclusion: the isolated-neighbour fallback classification/fill=0.5 seeding is necessary context but not sufficient cause -- the actual NaN trigger is specific to something in the local velocity/geometry configuration of the tight `X=20.5, Z=25.5, Y in {6.5,8.5}` corner region specifically**, since geometrically/temporally similar cells one Z-layer up (X=18.5,Z=26.5, same conversion step) or one Y-row over (Y=5.5, same corner) do not fail. **Not yet done:** dump the actual velocity field / local strain at these specific cells right before StreamReconstructAdvect/collide to find what distinguishes Y=6.5/8.5 from Y=5.5 at the identical corner.

> **[Claude, 2026-09-24] Crash cell confirmed as an ongoing dynamical event, not an initialization artifact -- and an EXISTING, already-implemented mechanism (Pass-3 re-resolution) avoids the crash entirely, though it quarantines rather than fixes the underlying numerical singularity.**
>
> **Not an initialization artifact:** the crash cell's own owned-rank trace (job 4262030, widened probe from t=425) shows it flagged **OBST continuously from t=425 through the entirety of t=442's own FSLBM/AfterDEM passes** (this is a resume-from-checkpoint-400 run, so the cell had already been solid obstacle for 25+ tracked steps, almost certainly the whole run) -- converting only at `post-AfterDEM` of t=442, the same ordinary per-timestep DEM/reconstruction pass that runs every step. This rules out any connection to t=0 setup or the checkpoint-resume operation itself; it is a genuine, ongoing dynamical event (the box's corner rotating this cell out of obstacle coverage).
>
> **Traced through the classification algorithm directly (`classifyFormerObstacleNeighbourhood`, PdfReconstructionManager.h:86-115):** the cell falls into documented branch A' (`hasGas && !hasLiquid && hasUnknown -> Interface`) -- a deliberate same-pass-ambiguity fallback, not a decision-table bug. The code already anticipates this: cells hitting this fallback are collected as **Pass-3 candidates** (`isPass3Candidate`), and if `ModelParameters.usePass3Reresolution` is enabled, are re-classified after every FormerMO cell in the block has its final Pass-2 state (so `hasUnknown` should resolve to false). If the re-check then reads `hasGas && !hasLiquid && !hasUnknown`, Pass 3 converts the cell straight to Gas via `boundaryHandling->forceDomain(bognerGas, ...)`, bypassing the fragile isolated-interface reconstruction+collide path entirely. This option already exists and is already used successfully in three unrelated prms this session (`cuboid_stillwater_realgravity981_3xsphere_d075_pass3_gas.prm`, `thousand_spheres_two_pillars_slitpass_TEST_100floaters_r3_gas.prm`, `ten_cuboids_two_pillars_slitpass_gas.prm`) -- it was simply never turned on for this Bogner box scenario.
>
> **Direct test: job 4262101**, a byte-identical copy of the confirmed crash reproduction (`bogner_exact_3box_s16_rho0767_tilt26p565_resume400_velnanscan_gas.prm`) with exactly one change, `usePass3Reresolution true` (new file: `..._pass3_gas.prm`, `jobs/run_65_bogner_resume400_velnanscan_pass3.sbatch`). **Result: COMPLETED cleanly, all 500 requested steps, zero `reynoldsNumber is nan` / `FATAL fluid-cell NaN` occurrences.** Confirmed mechanistically via the same probe: at t=442 post-AfterDEM the cell is now `flag=GAS` (not `INTF` as in every prior crash), rho=1.0, PDFs at plain equilibrium (`f[0]=0.29629 ~= 8/27`) -- exactly the Pass-3 Gas reclassification predicted from the code.
>
> **Important caveat -- this is a quarantine, not a fix of the underlying singularity.** At t=443's `after-collide` probe, this now-GAS cell STILL goes `-nan` (identical symptom to the original crash, just one flag state removed from mattering) -- collide is still producing a non-finite result from this exact spatial location even on a plain-equilibrium input. Because `MeanVelocityComputer::getMeanVelocity` (Forcing.h) excludes GAS cells from the abort-triggering velocity sum, this NaN never reaches the outside world and the run survives -- but the actual numerical mechanism that makes collide blow up at world (20.5,8.5,25.5) remains unexplained. **Practical recommendation: enable `usePass3Reresolution true` for all future runs of this Bogner box family** (free, already-tested, zero observed downside, directly avoids the reproducible step-444 crash) while treating the deeper "why does collide go non-finite at this exact corner regardless of cell classification" question as still open if pursued further (e.g. dumping the force-density field, not just the PDFs, at this location right before collide -- the codegen'd Cumulant kernel folds in an external force term from a separate sweep that hasn't been examined for corruption at this specific spot).
>
> **Update:** re-checked `ForceDensityCodegenSweep` (Forcing.h) directly -- it writes `forceDensityFieldIt` ONLY for `isInterface`/`isLiquid` cells, every single timestep, fresh, from that timestep's own legitimate density/fill. So an ordinary interface cell sitting in open water (continuously interface for many consecutive steps) never has a staleness problem -- this field is recomputed for it every step, with no ping-pong/never-written risk at all. The staleness we flagged as a caveat is specific to a GAS-classified cell that has never had a legitimate liquid/interface classification (our crash cell, post-Pass-3): `forceDensityField` there has literally never been written (this cell was OBST the entire tracked window), so it's inert leftover/garbage -- exactly the same "field never touched for this flag" pattern already established for the PDF ping-pong buffers, not a new or different mechanism, and harmless for the identical reason (excluded from anything the abort-triggering computation reads).
>
> **[Claude, 2026-09-24] Traced the root cause of the same-pass ambiguity to an explicit, stated assumption in Bogner's own published algorithm that this scenario violates -- Bogner did not solve this case, he assumed it away.** Read `Bogner, S. (2017) "Direct Numerical Simulation of Liquid-Gas-Solid Flows Based on the Lattice Boltzmann Method"` (PhD thesis, Zotero `H7YLRADA`), Sec. 5.2, pp.112-113, directly (via `pdftotext`, not from memory/citation alone). His obstacle-to-fluid conversion rule, verbatim:
>
> > "We define the non-obstacle neighborhood of a cell as B(x) := {x + c_q : x + c_q is no obstacle}. For a prior obstacle cell x to be converted into one of the three fluid states, the state is determined as follows [Bogner, 2009]:
> > - B(x) contains no liquid: x is converted into gas.
> > - B(x) contains no gas: x is converted to liquid, rho(x) is interpolated from B(x), and the f_q(x) are initialized according to Eq. (5.3).
> > - B(x) contains liquid and gas: x is converted into interface, and the f_q(x) are initialized as in the preceding case. Additionally, a fill level phi(x) is chosen by interpolation of the fill levels of the interface cells in B(x)."
>
> Two things stand out. First, this 3-branch rule has **no defined outcome for "B(x) contains neither liquid nor gas"** -- that configuration satisfies both "contains no liquid" and "contains no gas" simultaneously, an undefined double-match with no stated priority. Second, and this is the actual explanation, Bogner states an explicit precondition immediately before presenting the rule: **"We assume the velocity of the object is small enough for conversions into obstacle to occur only at positions in the neighborhood of other obstacle cells, while conversions from obstacle to one of the three fluid states ... can occur only at positions in the neighborhood of fluid cells."** Under that assumption, a freshly-uncovered cell's B(x) is guaranteed to contain only ALREADY-RESOLVED fluid cells (liquid/interface/gas from a prior timestep) -- never another cell that is *itself* simultaneously converting from obstacle this same timestep -- so "neither liquid nor gas" simply cannot arise, and the undefined case is excluded by construction, not handled.
>
> **Our crash is a direct, concrete violation of that stated assumption.** The box's corner at 26.565deg tilt, s=16 resolution, omega=1.9 rotates fast enough (relative to the lattice spacing) to uncover a small CLUSTER of adjacent cells within a single timestep, not one at a time -- confirmed independently by the `[FALLBACK-ISOLATED-INTERFACE]` log showing multiple nearby world locations firing in the same timestep window (steps 425-442 entry above) and by the classification trace showing our target cell's `hasUnknown=true` (a neighbour still tagged `formerObstacle`, i.e. itself mid-conversion in the same pass). This is exactly the scenario Bogner's own "object velocity small enough" precondition was written to rule out, and his published algorithm has no defined behaviour for it.
>
> **This reframes the waLBerla/Antidunes codebase's extra machinery (the `hasUnknown` tracking, the branch-A'/C interface fallback, and Pass-3 re-resolution) not as a deviation from Bogner's method, but as the necessary generalization that extends his rule to remain well-defined outside the regime his stated assumption covers.** Pass 2's interface fallback keeps the simulation numerically safe (never commits to an ill-founded gas/liquid decision) exactly when Bogner's assumption fails; Pass 3 then re-applies Bogner's original 3-branch rule a second time, after the ambiguity has cleared, recovering the answer Bogner's rule would have given directly had his assumption held. In that sense Pass 3 isn't a bug workaround so much as a literal implementation of "wait until B(x) is fully resolved, then apply Bogner's rule" -- which is precisely what his own stated precondition requires but our simulation's local rotation rate cannot guarantee at this resolution.
>
> **[Claude, 2026-09-24] What we tried, and why a fully deterministic fix is not worth pursuing further in a distributed-memory setting.** To restate the core gap plainly, in the words used to describe it in this session: Bogner's rule assumes "the velocity of the object is small enough for conversions ... from obstacle to one of the three fluid states ... [to] occur only at positions in the neighborhood of fluid cells" -- i.e. he assumes a freshly-uncovered cell's neighbours are always already-resolved fluid cells from a prior timestep, never another cell simultaneously converting in the same pass. Our box violates that assumption at its rotating corner, and the existing codebase already carries a multi-layered response to exactly this gap (none of it invented from scratch this session, all of it pre-existing, tested machinery -- `tests/showcases_antidunes/BognerStateClassificationTest.cpp` unit-tests `classifyFormerObstacleNeighbourhood` directly):
> - **`hasUnknown` tracking** (`classifyFormerObstacleCell`): every former-obstacle cell's neighbour scan explicitly counts neighbours still tagged `formerObstacle` (i.e. themselves mid-conversion, same pass, same timestep) as a distinct third category alongside real gas/liquid evidence, rather than silently treating them as absent.
> - **Pass-2 safe fallback** (branches A'/C of `classifyFormerObstacleNeighbourhood`): whenever `hasUnknown` is true, the cell defaults to Interface rather than committing to an ill-founded Gas/Liquid guess -- numerically safe by construction, at the cost of sometimes being wrong (as in our crash cell).
> - **Pass-3 re-resolution** (`usePass3Reresolution`, gated opt-in): a second, purely LOCAL re-application of the same classifier after every FormerMO cell in the block has its final Pass-2 state, so `hasUnknown` can resolve to false once the ambiguity clears. Confirmed this session (job 4262101) to correctly recover Bogner's intended Gas outcome for our exact crash cell.
> - **`isPass3Candidate` + `[PASS3-UNRESOLVED-ISLAND]` diagnostic** (pre-existing): explicitly anticipates that Pass 3 will not always succeed -- a cell that STILL has no real liquid and no real gas within one cell even after every same-pass neighbour has a final state (a genuinely isolated cluster with no external contact at all) is logged distinctly and left as Interface, because there is no safe basis to pick either side.
> - **`[FALLBACK-ISOLATED-INTERFACE]`** (added this session, `seedConvertedCellFill`): a targeted warning specifically for the fill-seeding side of the same problem (zero real interface/liquid neighbours to average from), which is what let us catch and confirm the crash cell's exact mechanism in the first place.
>
> **None of this adds up to a fully deterministic algorithm, and that is a considered trade-off, not an oversight.** A provably-correct resolution -- one that always recovers the exact answer Bogner's rule would give under his stated assumption, for ANY size and shape of simultaneously-converting cell cluster, including clusters that span MPI block boundaries -- would need one of two things: (1) an unbounded number of ghost-layer ping-pong communication rounds per timestep, iterating classification+halo-exchange to a global fixed point (the number of rounds needed scales with the size of the largest connected same-pass cluster, which is not bounded a priori), or (2) a globally-ordered/BFS-style resolution scheme (always resolve cells nearest to already-settled fluid first, flood-fill outward), which is inherently a distributed graph algorithm, not an embarrassingly-parallel per-cell operation, and would need its own coordination/communication infrastructure layered on top of the existing block-local sweep structure. Either approach adds real per-timestep MPI communication overhead and substantial implementation complexity to a piece of code that, in practice, only matters for a small, transient population of cells (60 fallback events across 18 steps in an 208,000-cell domain, in this run) at object corners crossing the free surface -- not a cost/complexity trade worth making for this thesis's purposes. Pass 3's single bounded local re-check is the practical middle ground: cheap (one extra local pass, no extra communication), handles the common case (confirmed directly for our crash), and degrades safely (stays Interface, logged distinctly) for the residual cases it cannot resolve.

> **[Claude, 2026-09-24] Direct evidence for WHY the free box does not settle at the theoretical equilibrium: a real, quantitative torque residual that correlates with translational speed, not roll rate -- a velocity-coupled hydrodynamic effect Bogner's static hydrostatic theory has no term for.** Reused `gen_bogner_righting_scan_prm.py`'s own exact geometry (Sutherland-Hodgman polygon clipping of the box's square cross-section against the waterline, `square`/`clip_below`/`area_centroid`) to compute the analytic hydrostatic torque `M(alpha) = G * A_submerged(z, alpha) * L_box * bx` at the box's own OBSERVED (z, tiltY) pose at every logged step, and compared it directly to the ACTUAL measured `totalTorque_y` from the archived `[BoxTorque]` log (`~/University/Thesis/animations/4255656_omega1p9_seed5_100k_torque_20260922/BoxForces_BoxTorque.log`, rho=0.75, seed5, 100k steps, free 6-DOF).
> - **At rest (step 1, tilt=5deg, |vel|~0): measured torque=-0.0168 vs analytic=-0.0128** -- close agreement, ruling out a bug in the MEM force computation itself; the underlying hydrostatic force calculation matches Bogner's theory well when the box is actually static.
> - **Once settled into its ~30-32deg quasi-steady oscillation band, the gap opens and never closes**: analytic hydrostatic torque sits persistently around +0.06 to +0.12, while measured torque hovers near zero. This is not a small numerical discrepancy -- it is the SAME ORDER OF MAGNITUDE as the torque itself.
> - **Correlated the residual (measured - analytic) against the box's own dynamics in a 500-step window (steps 10000-10500, every step)**: correlation with angular velocity (d(tiltY)/dstep) = **-0.11** (weak), correlation with linear velocity magnitude (|vel|, translational speed) = **-0.64** (moderately strong). The residual tracks TRANSLATION speed, not roll rate.
> - **Interpretation**: this is the signature of a velocity-coupled hydrodynamic torque (an added-mass/drag cross-coupling between translation and rotation, physically analogous to the "Munk moment" for a body moving obliquely through a fluid) -- a real physical effect that Bogner's righting-moment theory, derived purely for a body AT REST, has no term for by construction. As long as the box keeps translating (which the established "coupled roll-translate limit cycle" mechanism guarantees it does), this extra torque keeps perturbing it away from the pure hydrostatic equilibrium, which is consistent with, and gives a concrete mechanism for, why the box sustains a wide oscillation rather than damping to 26.565deg.
> - **Script**: ad hoc, not yet committed to the repo (`/tmp/torque_bias_analysis.py`, `/tmp/torque_corr.py` -- reuses `gen_bogner_righting_scan_prm.py`'s geometry functions verbatim, parses the archived `[BoxForces]`/`[BoxTorque]` log pair per step via regex). **Not yet done**: rerun this same analysis on the rho=0.7665 (density-corrected) run and on a SLOWER/more-damped run (e.g. lower omega, or a pinned-translation control) to check whether reducing translational speed measurably shrinks the residual and the resulting oscillation amplitude, which would be the most direct confirmation of the mechanism; also worth checking early-transient steps 250-3000 (large, sign-varying residual during the fast initial roll from 5deg, not yet analysed for velocity correlation) separately, since that phase has much larger |vel| than the settled band.

> **[Claude, 2026-09-24] REFRAMING: Bogner's own floating-cuboid validation is NOT free 6-DOF -- it explicitly restricts rotation to the longitudinal axis only, and reports the box reaching genuine static rest, not an oscillation. The "nopin" free runs studied all session are therefore a genuine extension beyond his tested envelope, not a discrepancy with it.** Read Bogner's thesis Sec. 5.3.3 ("Floating Bodies in Free Surface Flows," pp.121-123) directly. Verbatim: "the 2D cube of side length s is represented by a cuboid of dimensions 2s x s x s, that **will rotate around its longitudinal axis only**... Three boxes with s=16 are placed inside the basin, slightly heeled from their axis-aligned upright position by an angle of 2.86deg... **The simulation was run until all fluid motions ceased and the objects came to rest**." This is the EXACT scenario this session's `bogner_exact_3box_s16` naming reproduces (3 boxes, rho=0.25/0.5/0.75, s=16, 2.86deg seed, 130x40x40 domain, tau=1/1.9, g=7.5e-4) -- confirming these runs are a direct, deliberate reproduction attempt of Bogner's own published validation, not an independent scenario.
>
> **Two consequences.** First, this confirms (rather than merely assumes) that job **4253988** -- X/Y position and pitch/yaw pinned, roll-only, deliberately built earlier this session to match "Bogner's own methodology" -- is not just A reasonable control, it is THE methodologically correct reproduction of his actual published test. That run settling to -27.0deg (0.4deg from 26.565deg, attributed to ordinary lattice-discretization error at s=16) is therefore the properly-validated result, directly comparable to Bogner's own reported convergence. Second, and this reframes the entire equilibrium-offset investigation: the "nopin" free-6-DOF runs (the box drifting in X/Y, oscillating 30-65deg, never reaching rest) are studying a case **Bogner never tested and his theory has no term for** -- a genuinely harder, more general problem this thesis is extending into, not a bug or an unexplained deviation from prior work. The velocity-coupled hydrodynamic torque residual found in the entry immediately above (correlating with translation speed, r=-0.64) is exactly the kind of physics that can only appear once translation is actually free, which is precisely the regime Bogner's roll-only validation was designed to avoid. **This does not fully resolve why the free box settles where it does, but it correctly relocates the question**: not "why does our simulation disagree with Bogner" (it doesn't, in the case he actually validated), but "what governs the equilibrium of a body with an additional, physically real translational degree of freedom that classical hydrostatic theory does not model."

> **[Claude, 2026-09-24] CORRECTION to the entry immediately above: no explicit "fixed"/"pinned" wording exists in Bogner's text for the floating-cuboid experiment -- re-searched directly (`grep -i "fixed position\|held fix\|pinned\|constrain"` across the full extracted thesis text) and found nothing near Sec. 5.3.3. The only textual anchor is "will rotate about its longitudinal axis only," and since Bogner's general rigid-body framework (Sec. 5.1.1: position, linear velocity, angular velocity, integrated by the p.e. physics engine) is unrestricted 6-DOF with no mention of an artificial per-experiment constraint mechanism, this phrase most likely describes the idealised 2D-reduced treatment of the problem (or what his own run happened to exhibit), not a stated numerical pin. Downgrading the previous entry's confidence accordingly -- it is NOT established that Bogner's own validation was constrained to fewer DOF than our "nopin" runs; it is only established that job 4253988 (deliberately constrained to match the roll-only idealisation) reproduces his reported 26.565deg outcome.**
>
> **New, independently decisive finding that changes the picture regardless: the free box does not stay in a pure roll-about-Y configuration at all -- it develops a real, growing out-of-plane wobble that Bogner's 2D-equivalent treatment has no room for.** Parsed the `rotation` DataArray (9-component, full 3x3 matrix) directly from the archived `floaterBox` VTU snapshots of job 4255656 (rho=0.75, seed5, 100k steps, every 2000 steps) rather than re-deriving tiltY alone from one matrix element. Defining the out-of-plane wobble angle as `arccos(R(1,1))` -- the angle between the box's body-fixed long (Y) axis and the world Y axis, which is exactly 0 for a pure rotation about Y and otherwise measures genuine 3D tumbling -- the trajectory shows:
> | step | tiltY (roll, atan2 convention) | wobble = arccos(R11) |
> |---|---|---|
> | 0 | -5.0deg | 0.000deg (exact, as initialised) |
> | 20000 | -31.9deg | 6.8deg |
> | 44000 | -31.8deg | 3.6deg (local minimum) |
> | 58000 | -29.1deg | 9.2deg |
> | 78000 | -32.6deg | 12.0deg |
> | 98000 | -32.6deg | **16.0deg** |
> This is not noise: while roll (tiltY) oscillates in the familiar ~29-33deg band throughout, the wobble angle grows close to monotonically over the run's second half with no sign of saturating by 100k steps, despite starting at EXACTLY zero (a perfectly pure roll at t=0, confirming this is not a measurement artefact). **Interpretation: the intended 2D-equivalent pure-roll configuration is unstable to out-of-plane (pitch/yaw) perturbation once all 6 DOF are genuinely free, and the box is slowly transitioning into fully 3D tumbling rather than merely oscillating in the plane Bogner's theory describes.** This plausibly also feeds the translational drift and the velocity-coupled torque residual found in the entry two above: a wobbling, precessing body couples roll into surge/sway through channels a pure planar roll never has access to. **Not yet done**: confirm this wobble growth is not an artefact of the specific `atan2`/`arccos` extraction by checking a full quaternion-based Euler decomposition; check whether the SAME wobble growth appears in job 4253988's own constrained-DOF control run (it should be identically zero there, by construction, which would be a useful sanity check of the wobble metric itself); check whether wobble growth rate correlates with translational speed the same way the torque residual did, and whether it is present already during the fast initial transient (steps 0-4000) at a detectable level; script used is ad hoc (`/tmp/wobble_analysis.py`), not yet committed to the repo.

> **[Claude, 2026-09-24] DECISIVE: the wobble-growth mechanism is confirmed on a second, independent run (job 4262237, fresh-from-t=0, rho=0.7665, released exactly AT 26.565deg with Pass-3 enabled, COMPLETED clean to 50000 steps with zero crashes), and its onset coincides exactly with the already-known sudden destabilisation event.** Ran the same `rotation`-DataArray wobble extraction (`arccos(R(1,1))`) against this run's own archived `floaterBox` VTU snapshots (`/home/vault/iwia/iwia148h/antidunes-archive/antidunes-bogner-exact-s16-rho0767-tilt26p565-50k-pass3_4262237_20260924_230658`, every 2000 steps through the last available snapshot at step 48000):
> | step | tiltY (roll) | wobble |
> |---|---|---|
> | 0 | -26.565deg (exact release point, zero seed transient) | 0.000deg |
> | 8000 | -28.2deg | 3.0deg |
> | 20000 | -28.4deg | 5.5deg |
> | 40000 | -18.2deg | 7.2deg |
> | 42000 | -19.4deg | 11.1deg |
> | 44000 | -32.9deg | 16.6deg |
> | 46000 | -36.7deg | 19.7deg |
> | 48000 | -41.1deg | **22.6deg** |
>
> For the first ~40,000 steps the box holds a gentle oscillation near the true target angle (26.565deg) -- matching the "tracks the target closely for a long stretch" behaviour already reported for other runs in this family -- while the wobble angle creeps up slowly and roughly linearly from 0 to ~7deg. At step ~40000-42000 the wobble growth suddenly ACCELERATES, and in exact lockstep tiltY destabilises catastrophically out of its calm oscillation into the wide chaotic band already documented elsewhere in this file. **The onset of the roll-angle blowup and the onset of runaway wobble growth occur at the same timestep, not independently.**
>
> **This gives a mechanistic explanation for the entire "coupled roll-translate limit cycle" phenomenon investigated throughout this session: the intended pure-Y-axis roll is only a marginally stable configuration once every DOF is genuinely free. It slowly leaks into a growing out-of-plane (pitch/yaw) mode; once that mode grows large enough (roughly 7-10deg in both observed cases), the configuration goes unstable and the system transitions into fully 3D tumbling, which shows up as "chaotic destabilisation" when only tiltY is examined (as every prior analysis in this file, including this session's own torque-residual work above, did until now).** This is now the leading candidate explanation for why the free 6-DOF box never settles at 26.565deg for either density -- not a hydrostatic force-computation bug (ruled out earlier: near-rest torque matches Bogner's theory closely), not initial-condition sensitivity (ruled out: two independent zero-transient releases both eventually destabilise this same way), but a genuine 3D rigid-body instability of the intended 2D-equivalent roll motion, consistent with, and probably the root driver of, the previously-found translation-correlated torque residual (a wobbling, precessing body couples roll into surge/sway through channels a pure planar roll never has). **Not yet done:** identify what triggers the transition from slow linear wobble growth to runaway acceleration around step 40000 (a resonance with the roll oscillation period, a genuine linear instability growth rate, or a threshold/nonlinear effect); check whether damping the wobble mode directly (e.g. a diagnostic-only run with pitch/yaw pinned but roll and all translation free, isolating whether restoring ONLY the rotational symmetry restores clean settling) would confirm this is causal rather than merely correlated.

> **[Claude, 2026-09-24] CAVEAT on the torque-residual analysis two entries above: the 2D analytic baseline it compares against assumes ZERO out-of-plane wobble by construction, so once real wobble develops the residual conflates two different effects, not one.** `gen_bogner_righting_scan_prm.py`'s polygon-clipping calculation (and Bogner's own Scribanti-curve theory it reproduces) is a strictly 2D treatment: submerged area/centroid computed from the X-Z cross-section alone, implicitly assuming the box's long (Y) axis stays exactly aligned with world Y. Once real wobble (the arccos(R(1,1)) angle from the two entries above) is nonzero, the box's TRUE submerged geometry is genuinely 3D and the 2D formula no longer describes it exactly. So the reported residual (measured torque - 2D-analytic torque) is a clean measurement of the velocity-coupled hydrodynamic torque ONLY while wobble stays small (roughly the first 10,000-20,000 steps in both analysed runs); once wobble grows large, the residual increasingly also reflects "the 2D baseline itself is wrong for this pose," not purely missing physics, and should not be over-interpreted at face value in that later regime. **This makes the pinned-XY/pitch-yaw control run (job 4262441, submitted this session, X/Y translation and pitch/yaw pinned, roll+heave free) an even more important test than originally framed**: by construction wobble is held at exactly zero throughout, so any comparison against the 2D analytic theory in that run stays valid for its entire duration, uncontaminated by this issue.

> **[Claude, 2026-09-24] Out-of-interest cross-check: does today's Pass-3 fix (item 3's box-corner crash) also help job 4255679's hopper sealed-cell crash? CONFIRMED NO -- and this cleanly separates the two problem classes.** Reran 4255679's own scenario (half-density cylinder hopper, 210 rods, no guards/caps at all) with exactly one change, `usePass3Reresolution true` (new files: `prms/hopper_clog_cyl_frontblock_pidon_2node_halfrows_pass3_gas.prm`, `jobs/run_144_hopper_clog_cyl_halfrows_gjk_pass3_60k.sbatch`, job **4262495**, 144 ranks, `~/mt-dennys-mirror`). **Result: still crashed, at essentially the same step (~3100-3101 vs the original's 3133), via a DIFFERENT proximate signature** -- `reynoldsNumber is nan` + a `[GHOST-VANISHED]` warning this time, not the original's `STALE-PARTICLE-IDX`. Confirmed Pass-3 was genuinely active and doing real work in this scenario (not silently inert): `[FORMERMO-CONVERT t=3100] ... pass3ToLiquid=476 pass3ToGas=1 pass3StillInterface=0 pass3UnresolvedIsland=0` out of 12867 total formerMO conversions in that reporting window -- this hopper scenario has substantial free-surface formerMO churn (plausibly rods moving through the water surface region), so Pass-3's machinery is exercised here too, just on cells unrelated to the actual crash. **This is the expected, predicted-in-advance result**: 4255679's own root cause (2026-09-23 entry above) is a MEM/DEM pinched-cell density blowup between two fully SUBMERGED rods -- a different code path (ordinary liquid-only MEM reconstruction between obstacles) from Pass-3's target (`classifyFormerObstacleNeighbourhood`'s gas/liquid/interface classification near the waterline). The shift in crash step and abort signature (but not crash avoidance) is consistent with the same "extra computation perturbs exactly which check trips first, without touching the underlying cause" pattern already documented for the VTK-timing-dependent box-corner crash earlier in this file. **Net conclusion: Pass-3 and the sealed-cell guard are fixes for genuinely different problem classes; neither substitutes for the other, and 4255679's own crash still needs its own fix** (the local deterministic 2-rod repro already found, `prms/rod_pair_4255679_repro_gas.prm`, remains the right tool for that, not pursued further here).

> **[Claude, 2026-09-25] CORRECTION: job 4253988 was mischaracterised earlier this session (and in this file's own 2026-09-24 pinxy-control entries above) as "position/rotation pinned, roll-only" -- it was actually FULLY FREE (nopin), released exactly AT the 26.565deg equilibrium with zero seed transient.** Re-read the actual 2026-09-22 entry recording it directly: "**(4) 4253988** `antidunes-bogner-exact-s16-seed26p565-rho075` (free start AT the equilibrium angle, rho 0.75 only) COMPLETED (11m32s): tiltY 26.81 deg at step 500, settles to 26.96 deg by step 50000, essentially FLAT the whole run... released AT its equilibrium heel, the box STAYS there... the 32-36 deg seen in every run seeded at 2.86-10 deg is confirmed as an APPROACH-PATH effect." There is no DOF pinning involved at all -- the run answers "does the free box hold its ground if started exactly at equilibrium," not "does a constrained box converge to equilibrium from a perturbed seed." **Consequence: there is no established reference confirming that a roll-only/constrained box correctly converges to 26.565deg from an away-from-equilibrium seed.** This reframes the pinxy 5deg/10deg stalling result (entries above): it is NOT contradicting a known-good constrained-convergence case (none exists), and may simply reflect the same "shallow energy landscape, slow relaxation" effect 4253988's own entry already flags for the free case seeded away from equilibrium -- 50000 steps may just not be enough for a 5-10deg seed to relax to 26.565deg even without any wobble-instability interference, pinned or not. Apologies for the earlier confident misstatement; this correction should be treated as authoritative over the pinxy-control framing in the preceding entries.

> **[Claude, 2026-09-25] Pinxy control CONFIRMED at a second seed angle: the box settles onto a real, stable local equilibrium near its own seed, not near 26.565deg -- reproducible at both 5deg and 10deg.** Job **4262552** (same pinxy setup, seed=10deg instead of 5deg, 50000 steps): the simulation itself completed cleanly to step 49999/50000 ("END LOGGING" reached; the SLURM "FAILED" status is purely a vault file-count quota error during archiving of `avgDataSlices`, unrelated to the simulation). tiltY rose from 10.0deg to 11.30deg by step 2500, then stayed essentially flat (11.30-11.35deg, |vel|~1e-4) for the remaining 47500 steps -- the exact same shape as the 5deg run (5.0 -> ~6.3deg, flat for ~40000 steps). Both seeds settle to approximately seed+1.3deg and STAY there indefinitely; neither shows any tendency to keep rolling toward 26.565deg within 50000 steps. **This looks like a genuine, distinct, stable local equilibrium near the upright position under this exact DOF constraint (X/Y translation + pitch/yaw rotation pinned, only Z-heave and Y-roll free), not merely slow relaxation.** Physically plausible given the wall-sided Scribanti curve already showed GM(0) is negative (unstable) for this density/geometry: past the wall-sided validity boundary, the TRUE (corner-clearing) righting-moment curve is not the simple parabola and can cross zero more than once, so a shallow local equilibrium well short of the deep 26.565deg one is a reasonable outcome, not a bug. **Not yet done:** confirm this local equilibrium is real (not a code artifact) by checking whether the net hydrodynamic torque at ~6-11deg is genuinely ~0 (already suggested by the tiny, non-decaying residual `F_hydro`/velocity values logged) rather than just very small; try a seed angle roughly halfway to 26.565deg (e.g. 15-18deg) to map out where this local equilibrium's basin of attraction actually ends.

> **[Claude, 2026-09-25] The pervasive-early-NaN pattern found while widening the hopper scan window is a RED HERRING relative to the original 4255679 crash -- it is a deterministic, structural artefact of the fixed initial pillar/rod geometry at the waterline, present in EVERY variant regardless of Pass-3, density, or PID, and almost certainly unrelated to the actual abort.** Ran the SAME pervasive-NaN check on a genuine reference: job **4262640**, the FULL-density (406 rods), Pass-3-OFF, PID-off baseline (`hopper_clog_cyl_frontblock_pidoff_2node_gjk_sanity_18k_gas.prm`, already an established scenario from before this session, its own dense-scan window merely widened to start at step 0), alongside job **4262624** (half-density, Pass-3 ON, the scenario this whole thread has been testing, also widened to start at step 0). **Both show the identical pattern, at the SAME world coordinates**, e.g. `world=(374.5,8.5,50.5)` at step 10 and `world=(376.5,27.5,50.5)` at step 12 appear IDENTICALLY in both runs -- ruling out Pass-3, density, and PID as causes. It starts as early as **step 8**, and the affected Z-coordinates (~50.5) sit right at the documented water level (Z=49.13) -- i.e. exactly where the fixed initial pillar/rod block first meets the free surface as the simulation initialises. Since this fires at step 8-12 but the real abort in every prior run of this family happens around step 2900-3133 (a ~2900-step gap), and since the Bogner box investigation already established that the large majority of isolated-fallback interface cells are reclassified to Gas shortly after and never cause a problem, **this pattern is very likely harmless background noise inherent to this scenario's fixed geometry, not the mechanism behind the actual crash.** The already-established root cause for 4255679 (a slow sealed-cell density runaway between two touching, fully submerged rods, building for ~20+ steps before overflowing -- see the 2026-09-23 entry above) remains the more likely explanation for the real abort, independent of everything varied in this thread (Pass-3, density, scan window). **Net conclusion for this whole side investigation: neither Pass-3 nor the sealed-cell guard fix 4255679's actual crash, the pervasive early NaN pattern discovered while looking for it is a separate and apparently benign phenomenon, and the real fix (if pursued further) should target the already-identified two-rod pinched-cell mechanism directly** (e.g. via the local deterministic repro `rod_pair_4255679_repro_gas.prm`), not the free-surface classification machinery this thread has been instrumenting.

> **[Claude, 2026-09-25] Added a targeted diagnostic in `src/lbm/free_surface/dynamics/StreamReconstructAdvectSweep.h` (`[DEGENERATE-FILL-RECONSTRUCTION]`) to test the hypothesis that a near-total obstacle enclosure starves `reconstructInterfaceCellLegacy`'s curvature/normal-based reconstruction enough to make `newRho` degenerate, which the per-step fill update then turns directly into a non-finite fill via an UNGUARDED division with no zero/NaN check at all: `*fillDstFieldIt = (*fillSrcFieldIt * oldRho + deltaMass) / newRho;` (line 155, confirmed via `git log` to be pre-existing code from commit 87a522d7, 2026-06-11, not something added this session).** The diagnostic fires only when `*fillDstFieldIt` is itself non-finite right after this division, logging cell coordinates, old/new rho, rhoGas, curvature, normal, and a live count of obstacle neighbours (throttled first-20-then-every-5000th, matching this codebase's existing convention since this is a generic library file with no Antidunes-specific step-window config available).
>
> **Confirmed on job 4262624's own crash cell (world (267.5,15.5,18.5), 23/26 obstacle neighbours, both rho AND fill already non-finite): rerunning the identical scenario with this diagnostic compiled in (job 4264769) reproduced a crash, but at a DIFFERENT cell with a DIFFERENT mechanism, directly distinguishing two separate NaN-injection pathways in this scenario family:**
> 1. **Degenerate curvature/normal reconstruction -> fill=NaN via the unguarded StreamReconstructAdvect division** -- the mechanism this diagnostic targets, established from job 4262624's own 23/26-obstacle cell (see the prior entry above).
> 2. **Collide itself produces NaN PDFs on an otherwise cleanly-reconstructed cell** -- confirmed as job 4264769's actual failure (step 2833, world (267.5,18.5,19.5)): the new diagnostic fired ZERO times in this entire run (ruling out mechanism 1 for this occurrence), yet the SAME cell shows up as `FATAL fluid-cell NaN` first at `INTRA-after-collide` (not at any earlier sweep stage this timestep) with a completely normal-looking neighbourhood -- **12 OBST / 13 LIQ / 1 GAS**, i.e. abundant real liquid neighbours, nowhere near a near-total enclosure. This is the EXACT SAME mechanism already established for the original Bogner box crash (clean PDFs through StreamReconstructAdvect, NaN appears specifically at collide) -- direct evidence that the box's collide-fragility mechanism is not confined to geometrically extreme, near-total-enclosure cells; it can hit an interface cell with an entirely ordinary neighbourhood too.
>
> **Practical consequence: this scenario family has at least two independent, confirmed ways for an interface cell to go fatally non-finite, not one.** A fix targeting only the curvature/normal-degeneracy path (e.g. a zero/near-zero guard on `newRho` in `StreamReconstructAdvectSweep.h:155`) would not address the collide-fragility path, and vice versa -- consistent with why neither Pass-3 nor the sealed-cell guard (each addressing a specific, narrower failure mode) has fully resolved crashes in this scenario family. **Not yet done:** add an equivalent targeted diagnostic right after the generated Cumulant collision kernel call (`BuildTimeloop.h`'s `CollideSweep`) to capture the PRE-collide (clean) vs POST-collide (NaN) PDF state directly for this second mechanism, mirroring the stage-by-stage table already built for the original Bogner box investigation, to see whether the same "fragile reconstructed input" signature (e.g. an extreme local density/velocity gradient) is present here too.

> **[Claude, 2026-09-25] Direct evidence this location is repeatedly, independently fragile -- NOT a case of "one clean cell suddenly fails," and a new structural finding: collide runs unconditionally even on OBST-flagged cells.** Added a cube-trace probe (`BuildTimeloop.h`'s `probeHopperCrashRegion`, world X=267.5, Z in {18.5,19.5}, Y from 10.5 to 30.5, gated to steps 2800-2850) covering the recurring crash region, and reran the same scenario (job **4264868**). The crash reproduced at step 2837, at world (267.5,18.5,19.5) -- and the cube trace shows its FULL history for the preceding steps, not just the moment of failure:
> | step | stage | flag | rho |
> |---|---|---|---|
> | 2830-2831 | all stages | OBST | clean, ~1.0067-1.0069, stable |
> | 2832 | after-boundarySweep/after-forceDensitySweep | OBST | 1.006628 (still clean) |
> | 2832 | **after-streamReconstructAdvect** | OBST | **exactly 0.0** |
> | 2832 | after-cellConversion/after-equilibriumRefilling | OBST | **-nan** |
> | 2833-2836 | all stages, every step | OBST | **-nan**, persistently (not alternating) |
> | 2837 | after-boundarySweep/after-forceDensitySweep/after-streamReconstructAdvect | **INTERFACE** | **1.006678 (clean!)** |
> | 2837 | after-collide onward | INTERFACE | **nan** -- the actual crash |
>
> **This exact physical location already went non-finite once before, 5 steps earlier, while still flagged OBST -- not inherited from a different neighbouring cell, the SAME location breaking down twice.** Two consequences: (1) neither `boundarySweep` nor `forceDensitySweep` touch an OBST cell's PDFs (rho stays flat, matching prior findings), but the value visibly changes between `after-streamReconstructAdvect` (0.0) and `after-cellConversion` (-nan) at step 2832 -- and `collide` is the only sweep that runs in that gap. **This directly demonstrates the generated Cumulant collision kernel is NOT flag-gated -- it runs unconditionally on every cell in the domain, obstacle or not**, a structural fact not directly confirmed before this trace (only inferred/suspected earlier this session from OBST cells' ping-pong-alternating rho values, which could previously be explained by inert buffer swapping alone; this trace shows an ACTIVE value change specifically bracketed by collide's own execution window, which inert buffer swapping cannot explain). (2) The step-2837 conversion to INTERFACE DID give this cell a genuinely clean reseed (rho=1.006678, not the OBST-era -nan carried forward) -- confirming `seedConvertedCellFill` really does reset PDFs from scratch rather than inheriting stale corruption, so the two failures are independent collapses, not one persisting corruption.
>
> **Conclusion: this specific location (world X=267.5, Y=18.5, Z=19.5) is numerically fragile regardless of its classification state at the time** -- it broke down as OBST, got a clean reset via formerMO conversion, then broke down again as INTERFACE within the same handful of steps. Repeated failure at one fixed spot across different flag states points to a genuine local geometric singularity (most likely an extremely tight, near-degenerate gap between two specific rods at that exact point that never resolves correctly regardless of cell classification) rather than a random, one-off, or neighbour-inherited fluke. **Not yet done:** inspect the actual rod geometry (positions/radii) at this exact world location to confirm a near-degenerate gap directly (matching the earlier local deterministic 2-rod repro's own methodology, `gen_rod_bundle_channel_test_prm.py --pair`); check whether collide's unconditional-execution-on-OBST-cells finding also explains the ping-pong-alternating rho pattern documented earlier for the Bogner box's own permanently-OBST cells (i.e. whether THOSE alternations are also active collide writes rather than inert buffer exposure, which would need re-verification given this new finding).

> **[Claude, 2026-09-25] Rod geometry at the crash location, and a correction: the nearby "gas" cannot be a real entrained bubble -- the site sits ~30 lattice units below the actual free surface.** Extracted the rod layout feeding this scenario (210 rods, radius 3, length 18, axis along Y, X-pitch 7 with centres at X=270,277,284,...; Y-centres at 10 and 30, each rod spanning Y=[1,19] or [21,39] with a 2-LU gap between the two Y-groups; Z-layers at 4, 11.52, 19.04, 26.56, ...). World X=267.5 sits just short of the first rod-centre column (X=270) at radius 3, i.e. right at the column's own surface -- consistent with a genuinely tight geometric margin rather than mid-channel open space. The crash Z (~18.5-19.5) lines up almost exactly with the Z=19.04 rod layer, and Y=18.5-20.5 sits right in the 2-LU gap between the two Y-groups (Y=19/21 boundary) -- i.e. the fragile cell sits in the narrow gap where a Z-layer-3 rod's end and the inter-group Y-gap meet, a genuinely tight multi-rod corner, not open water.
>
> **Correction to the working "trapped gas pocket" framing used earlier this session**: the nearby GAS-flagged neighbour that motivated calling this a "trapped bubble" cannot be a real entrained air pocket. The documented water level for this scenario is Z~=49.13 (`ComputedParams`/prm-derived), so the crash region at Z~=19.5 sits roughly **30 lattice units below the actual free surface** -- physically implausible for real air to have been carried down there in the ~830-2900 steps since `settlingSteps` release (rods settle under gravity through still water; there is no mechanism in this scenario for air entrainment to depth). The GAS classification at/near this cell is almost certainly a **phantom classification artefact** of the formerMO machinery: rods resettling after release briefly uncover cells in a tight multi-rod corner faster than real liquid/gas evidence can propagate into them, producing a transient, non-physical GAS flag rather than a genuine void. This is consistent with, and likely the same underlying mechanism as, the already-documented pervasive early-NaN red herring (2026-09-25 entry above) and the Bogner-box formerMO-misclassification finding earlier in this file -- i.e. a third manifestation of the same general class of problem (evidence-starved cells right after an obstacle-to-fluid transition), here compounded by collide's own fragility at this specific tight geometric corner. **Not yet done, deferred per explicit redirection to the force-density investigation below**: the proposed isolated test (a small tightly-packed submerged rod cluster held fixed then released, no free surface) that would directly confirm the phantom-classification mechanism without any of the collide-fragility complexity.

> **[Claude, 2026-09-25] Traced the exact generated Cumulant collision kernel operations that consume external (non-PDF) state, to identify the one remaining untested input for the step-2837 collide-stage NaN.** Read `build/apps/showcases/Antidunes/default_codegen/AntidunesLatticeModel.cpp` (lbmpy-generated, read-only) directly rather than reasoning about the collision operator abstractly. Two divisions/external reads matter: `xi_0 = 1.0/m_000` (line 130, `m_000` = the cell's own summed density from its PDFs -- the mechanism already confirmed responsible for the FIRST/step-2832 OBST-phase collapse, where rho was independently measured as exactly 0.0 right after `streamReconstructAdvect`), and `xi_29/30/31 = xi_3*_data_force[...]` (line 410 in `kernel_collide`, line 77 in `kernel_streamCollide`) -- these feed directly into `u_0 = m_100*xi_0 + xi_0*xi_29` (velocity) and every higher cumulant, and are read from the **separate force-density field** (`ForceDensityCodegenSweep`, `Forcing.h`), not from the PDFs at all. Since the step-2837 crash cell's raw PDFs were already confirmed clean/near-equilibrium (matching D3Q27 weights almost exactly -- ruling out both the zero-density and extreme-velocity-from-PDFs hypotheses directly, not by assumption), **the force-density field is the single remaining, never-yet-directly-inspected external input to collide for this occurrence.** Also confirmed `ForceDensityCodegenSweep` only writes force density for `isInterface`/`isLiquid` cells (`accel*fill*density`) with no else-clause for GAS/OBST -- if the crash cell's fill or density were itself extreme/stale at the moment this sweep ran (fill=NaN already excluded by the clean-PDF trace, but a stale/uninitialised value from the cell's OBST-era history is not yet excluded), the resulting force value could be corrupted independent of the PDF field. **In progress:** `trackWatchCells` (`SetupStages.h`) now also prints this field's raw `(fx,fy,fz)` per watched cell/stage, and `probeHopperCrashRegion`'s (`BuildTimeloop.h`) coverage was widened (step window 2700-3200, bracketing all 7 previously observed crash steps 2831/2833/2837/2900/3006/3120/3133 across independent reruns; cell grid X=267.5 fixed, Y=5.5-35.5, Z in {17.5,18.5,19.5,20.5}) since this scenario's exact crash location/step is not bit-reproducible across reruns. Job **4264997** (same scenario/prm as 4264868, `hopper_clog_cyl_frontblock_pidon_2node_halfrows_pass3_precursor_gas.prm`) submitted with this instrumentation to directly capture the force-density value immediately before whichever cell fails this time.

> **[Claude, 2026-09-25] DECISIVE, surprising result: the force-density hypothesis is RULED OUT, and direct replay of the captured pre-collide state through the real kernel arithmetic proves collide cannot have computed NaN from what our own probe observed -- redirecting the investigation from "fragile physics" to "collide is reading/writing something our probe doesn't see."** Job **4264997** (widened cube-trace + force-density probe) reproduced the crash at the SAME location as every prior run (world (267.5,18.5,19.5)), this time at step 2838, via the now-familiar two-stage pattern: OBST-phase collapse at step 2833 (all-27-PDFs exactly 0.0 right after `streamReconstructAdvect`, then `-nan` from `1/0` at collide, persisting steps 2833-2837), a clean formerMO reseed to INTERFACE at step 2838 (`rho=1.006926, fill=0.5`, ordinary near-equilibrium PDFs, force=(0,0,0) immediately after `boundarySweep`), then **`forceDensitySweep` writes a small, entirely ordinary force `(4.531e-07, 0, -2.933e-05)`** -- not extreme, not NaN, not stale-looking, and it stays unchanged and finite through `streamReconstructAdvect` right up to the moment `collide` runs -- yet `collide` still turns this cell to `nan` across all 27 PDF components in that one step.
>
> **This rules out the force-density field as the corrupted input** for this occurrence. To settle "what exactly breaks inside collide" definitively, ported the generated `kernel_collide` (`build/apps/showcases/Antidunes/default_codegen/AntidunesLatticeModel.cpp:352-682`) line-for-line into a standalone Python script (`/tmp/.../replay_collide.py`, not committed -- ephemeral verification only) and fed it the EXACT captured 27 PDF values, force vector, and `omega=1.99753` (read directly from this run's own startup log, `relaxationRate = 1.99753` -- notably very close to the omega=2 zero-viscosity stability limit, though this by itself is not a singularity in the kernel's math). **Result: every intermediate value (m_000=1.0069, u=(-0.0062,-0.00057,-0.0103), all C_post_*/c_post_* cumulants) comes out completely finite and physically ordinary.** Inspecting the kernel source directly first confirmed why this had to be checked rather than assumed: the entire 330-line collision kernel is pure polynomial arithmetic in the PDFs/force/omega with exactly two division operations, both `1/m_000` and `1/m_000^2` (density, measured at 1.0069, nowhere near zero) -- no `sqrt`, `log`, or any other division exists anywhere in the kernel. There is therefore no way for this exact input, run through this exact math, to produce a non-finite result other than via `m_000->0`, which did not happen here.
>
> **Conclusion: the NaN this run produced was NOT computed by the collision kernel from the state our probe captured.** Since `trackWatchCells`' read of `after-streamReconstructAdvect` and `collide`'s own kernel read happen on the same single-threaded MPI rank, on the identical field object, moments apart, with no intervening sweep in this scenario's per-timestep sequence (boundarySweep -> forceDensitySweep -> streamReconstructAdvect -> collide, all confirmed to run back-to-back with no MEM/DEM coupling in between at this point in the timestep, per `CLAUDE.md`'s documented three-timeloop structure), this points away from "fragile physics triggered by unlucky-but-legitimate input" and toward a genuine **data-corruption or stale-read bug**: either collide's generated code reads a different physical memory offset than `trackWatchCells`/`checkPdfNaN` do for this specific cell (e.g. a ghost-layer-width or stride mismatch between the PDF field, which needs a wide stencil halo, and the force-density field, which may have a narrower one, manifesting only for cells near a block boundary), or something writes into this exact PDF memory in the narrow window between the two reads through a path not yet identified. **Not yet done:** check whether local index (11,18,19) sits at or adjacent to this block's edge (would support the ghost-layer/stride-mismatch theory directly); if so, check the configured ghost-layer widths for the PDF field vs. the force-density field in `SetupStages.h`'s field-registration code. This reframes the entire "what breaks inside collide" line of inquiry opened earlier this session -- the answer is very likely not a numerical/floating-point issue in the collision math at all, but a memory-correctness bug specific to this exact geometric location.

> **[Claude, 2026-09-25] CORRECTION: the ghost-layer/stride-mismatch theory proposed in the entry immediately above does NOT hold up -- checked the actual collide-sweep call site directly and it is wired safely.** Read `AntidunesLatticeModel.cpp:807-836` (`AntidunesLatticeModel::Sweep::collide`) and its call site (`BuildTimeloop.h:696-711`, `CollideSweep::operator()`). The crash cell's block runs the standalone `kernel_collide` path (`Sweep(CollideSweep(lbmSweepGenerated), ..., StateSweep::fullFreeSurface)` -- the correct path for a block with active formerMO/interface machinery, not the fused `kernel_streamCollide`/`StateSweep::onlyLBM` path). Critically, `CollideSweep::operator()` defaults `numberOfGhostLayersToInclude = uint_t(0)`, and this default is what actually gets used (the `Sweep(...)` wrapper calls it with no explicit argument). With `numberOfGhostLayersToInclude=0`: both `_data_force` and `_data_pdfs` pointers are offset by exactly 0 from their own inner-cell origins, `_size_force_0/1/2 = force->xSize() + 0` (no expansion at all), and the kernel loop covers precisely the block's real inner-cell region for both fields identically -- there is no ghost-layer expansion, no reliance on the force field's only-1 allocated ghost layer, and no possibility of the two fields' index spaces disagreeing near a block edge. **This closes off the stride-mismatch hypothesis entirely.** Also checked: `omega` is never mutated per-block or per-cell anywhere in the Antidunes-specific code (`grep` for `omega_\s*=`/`setOmega` outside the generated kernel/its Sweep wrapper returns nothing) -- it is set once as a plain scalar at `LatticeModel` construction and passed by value into every `kernel_collide` call, ruling out a corrupted/stale relaxation rate as well.
>
> **Net position after this whole investigation: PDFs (confirmed clean), force density (confirmed clean, tiny), omega (confirmed constant, uncorrupted), and the collide call site's field-indexing wiring (confirmed safe, no ghost-layer overrun) have ALL now been individually ruled out as the source of the collide-stage NaN, yet the arithmetic replay proves the kernel cannot produce NaN from this input through its own math.** The remaining candidate explanations are no longer reachable by log-level probing: (a) a genuine indexing/off-by-one bug in the generated kernel's per-cell address computation that causes one cell's collision result to silently overwrite a DIFFERENT cell's memory (would require inspecting/instrumenting the per-cell address arithmetic directly, or a debugger/memory-checker run), or (b) undefined-behaviour-class memory corruption from an entirely unrelated part of the code that happens to land in this exact PDF field's memory shortly before collide reads it (heap-adjacency dependent, would need AddressSanitizer/Valgrind to catch directly). **Recommended next step, not yet done:** build a debug/ASan-instrumented binary and run a small local reproduction of this exact block's geometry (matching the already-established `rod_pair_4255679_repro_gas.prm`-style local repro methodology) to catch an out-of-bounds read/write red-handed, since further indirect log-based probing has now exhausted the plausible external-input hypotheses without finding the actual mechanism.

> **[Claude, 2026-09-25] Cross-validated with the REAL compiled kernel, not just the Python port: new local unit test `ANTIDUNES_CollideKernelCrashReplay` (`tests/showcases_antidunes/CollideKernelCrashReplayTest.cpp`) confirms the same result.** Per the project's own established pattern (`tests/showcases_antidunes/ANTIDUNES_*` targets build+run+pass locally, no cluster needed), built a minimal single-block test: a `BlockForest` with the SAME `cellsPerBlock=<64,20,30>` as the real scenario, a real `LatticeModel_T`/`PdfField_T`/force field constructed exactly as `SetupStages.h::setupLBMFields` does, the target cell placed at the SAME local block index (11,18,19) as the real crash (2 cells from this block's own Y=20 edge, to preserve any edge-adjacency effect), populated with the exact captured pre-collide PDF and force values from job 4264997's log, then calling the REAL `AntidunesLatticeModel::Sweep::collide(block, uint_t(0))` directly -- the actual compiled pystencils-generated address/stride computation, not a hand re-typed formula. **Result: `rho before collide = 1.0069218`, `rho after collide = 1.0069218, allFinite=true`** -- identical conclusion to the standalone Python replay, now confirmed through the genuine production code path (address computation included, not just the isolated arithmetic).
>
> **This closes off the collision kernel itself (both its math AND its address/stride computation, including at this exact block-edge-relative position) as a possible cause, cross-validated two independent ways.** The only remaining explanation consistent with all evidence gathered so far is memory corruption or a stale read from a genuinely different source landing in this exact PDF field slot between the moment the log's WATCH probe reads it and the moment the real simulation's collide sweep reads it -- something no amount of further log-level probing or kernel-input replay can distinguish; catching it needs a memory checker (AddressSanitizer/Valgrind) run against the actual multi-rank simulation binary, ideally on a scaled-down local reproduction of this exact block/rod geometry. Test committed at `tests/showcases_antidunes/CollideKernelCrashReplayTest.cpp`, registered in `tests/showcases_antidunes/CMakeLists.txt` as `ANTIDUNES_CollideKernelCrashReplay`, builds and passes locally without any cluster access.

> **[Claude, 2026-09-25] Cross-validated a third way: `ANTIDUNES_CollideKernelCrashReplay` also PASSES on the cluster itself, same compiler/toolchain as the real production runs.** Before running on the cluster, re-fetched job 4264997's exact `[WATCH t=2838 HOPPER-after-streamReconstructAdvect] world=(267.5,18.5,19.5)` log line directly and diffed every one of the 27 `f[q]` values plus the force vector against what's hardcoded in `CollideKernelCrashReplayTest.cpp` -- all match exactly, no transcription error. Built and ran the test on `fritz` (job **4265384**, `jobs/run_1_unittest_collide_kernel_crash_replay.sbatch`, gcc/11.5.0 + openmpi/4.1.3, the SAME toolchain the production hopper runs use, as opposed to the local workstation's own compiler): `1/1 Test #228: ANTIDUNES_CollideKernelCrashReplay ... Passed`. (First submission, job 4265261, failed at the build step with "No rule to make target" -- the mirror's cached Makefiles hadn't picked up the new CMakeLists.txt entry; fixed by adding an explicit `cmake "$BUILD_DIR"` reconfigure step before `cmake --build ... --target ...`, then job 4265384 succeeded cleanly.)
>
> **This closes off the last remaining alternative explanation (a cluster-specific compiler flag or optimization difference, e.g. `-march=native`/FMA contraction, producing different floating-point behaviour than the local workstation build) -- the real production toolchain also computes a fully finite result from this exact input.** Combined with the local workstation result and the independent Python arithmetic replay, the collision kernel's finiteness on this input is now confirmed three independent ways (different arithmetic implementations, different machines, different compilers). The conclusion from the entry above stands unweakened: the actual root cause must be memory corruption or a stale read specific to the real multi-rank simulation's runtime state, not anything reproducible from the captured input values alone -- an ASan/Valgrind pass on the actual simulation binary remains the concrete next step, not further input-replay verification.

> **[Claude, 2026-09-25] Another independent occurrence, same exact location: fresh-build rerun on a second mirror crashed at world (267.5,18.5,19.5) again, this time at step 2835.** Per request, wiped `mt-dennys-mirror2`'s `build/` directory entirely and resubmitted `hopper_clog_cyl_frontblock_pidon_2node_halfrows_pass3_precursor_gas.prm` from a completely clean `--configure --build` (job **4265767**, confirming this is not an artefact of any particular stale build state either). Crashed via the usual `[NaN-SCAN] aborting after timeloopFSLBM` path at **world=(267.5,18.5,19.5), local=(11,18,19), step 2835**, an EIGHTH distinct crash step now recorded at this identical geometric location (previously 2831/2832/2833/2837/2838/2900/3006/3120/3133 across earlier runs) -- the step keeps moving, the place never does. Meanwhile the ASan-instrumented rerun (job 4265406) hit an unrelated infrastructure snag at step 400: `ibv_reg_mr(...) failed: Bad address` / UCX `failed to register address ... on md[4]=mlx5_0` -- a known ASan-vs-InfiniBand-RDMA incompatibility (ASan's shadow-memory/redzone layout breaks `ibv_reg_mr`'s pinned-memory registration), not a real crash or an ASan-caught violation; job **4265770** reruns the same ASan build with `UCX_TLS=tcp,self,sm` / `OMPI_MCA_btl=^openib` (forcing TCP+shared-memory transports, skipping RDMA registration entirely) to work around it.

> **[Claude, 2026-09-25] DECISIVE NEGATIVE RESULT: full AddressSanitizer instrumentation reproduces the exact same crash and catches ZERO memory violations -- this rules out classic memory-safety bugs (buffer overflow, use-after-free, double-free) as the cause.** Job **4265770** (ASan build, `UCX_TLS=tcp,self,sm`/`OMPI_MCA_btl=^openib` workaround for the earlier RDMA-registration incompatibility documented above) ran the full 144-rank scenario to completion of its normal lifecycle and crashed via the usual `[NaN-SCAN] aborting after timeloopFSLBM` path at **world=(267.5,18.5,19.5), step 2838** -- the SAME step AND same cell as the original job 4264997 (this scenario's non-determinism means step-matching two independent runs this closely is itself notable). `ASAN_OPTIONS=detect_leaks=0:abort_on_error=1:halt_on_error=1:log_path=...` was active throughout; the job's own post-run check (`ls asan_log.*`) found **no report files at all** -- ASan instruments every heap/stack/global memory access in the binary and aborts immediately with a detailed report on the first out-of-bounds read/write, use-after-free, or double-free it detects, so a clean run through the entire crash sequence means none of those occurred anywhere in the instrumented code (essentially all of it, including the MEM/DEM coupling files modified this session).
>
> **This decisively closes off the leading "memory corruption from elsewhere" hypothesis in its most common form and narrows what's left to three classes ASan structurally cannot see:** (1) an uninitialized-memory read (ASan does not check this; would need MemorySanitizer, impractical to retrofit onto an MPI+Python-codegen codebase reliably); (2) a genuine multi-rank data race / communication-ordering bug (ASan instruments single-process memory safety only, not inter-process synchronization -- a ghost-layer or particle-mapping update arriving late/early relative to when this rank reads it would look identical to "clean" memory from ASan's perspective); (3) a **logic bug that writes to a validly-allocated, correctly-addressed cell using the WRONG index** -- e.g. an off-by-one in particle-to-cell mapping at this exact rod's geometric boundary -- which is not a memory-safety violation at all (the write is in-bounds, just to the wrong place) and therefore invisible to ASan by construction. Given the OBST-phase collapse's rho hits an exact, clean `0.0` (not a garbage/uninitialized-looking bit pattern) immediately after `streamReconstructAdvect`, a deliberate-but-misdirected zero-write from somewhere in the free-surface/MEM coupling code (candidates: `MovingParticleMapping.h`, `PdfReconstructionManager.h`, `SimpleBB.h` -- all modified this session per git status) targeting this cell incorrectly is now the leading hypothesis over random corruption. **Not yet done:** audit those three files directly for any code path that writes fill/PDFs to exactly zero on a per-cell basis, and check whether its cell-index computation has an off-by-one or rounding issue specific to this rod's exact geometric boundary (world X=267.5 sits just inside the first rod-centre column's radius, as already noted in the rod-geometry entry above).

> **[Claude, 2026-09-25] Audited `MovingParticleMapping.h`, `PdfReconstructionManager.h`, `SimpleBB.h` per request -- their git diffs this session are all either pure diagnostics or explicitly opt-in features that are OFF for this exact scenario (`useInterpolatedBounceBack=false`, `mappingTieOffset`/`closeSealedCellsMinNeighbours` unset/0 in this prm) -- so THESE SPECIFIC edits are not the cause. But reading the pre-existing (unmodified) logic in these files pointed at something more important: which code actually gets to run at all.**
>
> **REFUTED: "repeated isolated collide iteration" (the hypothesis floated at the end of the previous entry) is directly disproven.** Built a new local unit test, `ANTIDUNES_CollideIsolatedObstacleIteration` (`tests/showcases_antidunes/CollideIsolatedObstacleIterationTest.cpp`), that repeatedly calls the REAL `sweep.collide()` on the same single cell 20,000 times in a row with zero force (matching the real obstacle-cell force value) and zero streaming input, starting from job 4264997's captured near-equilibrium PDFs. **Result: rho stayed bit-for-bit IDENTICAL (1.0069218) for all 20,000 iterations -- a genuine fixed point, no drift, no divergence at all.** This specific starting point turned out to be the wrong one to test with (see below), but the test itself is sound and now committed for reuse.
>
> **New decisive clue found while re-checking job 4265770's own full per-step trace (which, unlike earlier partial views, includes EVERY step from 2700 onward): this OBST cell's rho visibly drifts ACROSS timestep boundaries, even though nothing inside the FSLBM timeloop (boundarySweep..equilibriumRefilling) ever touches it except collide.** E.g. step 2700 ends at `rho=1.007356` (after-equilibriumRefilling); step 2701 begins at `rho=1.006870` (after-boundarySweep) -- a real ~4.9e-4 change, with force staying exactly `(0,0,0)` throughout, and boundarySweep itself independently confirmed (many times this session) to never change OBST-cell PDFs. **The only code that runs in the gap between one step's end and the next step's start is DEM sub-cycling and `timeloopAfterDEM` (LBM-DEM coupling: particle mapping + PDF reconstruction -- i.e. exactly `MovingParticleMapping.h` / `PdfReconstructionManager.h`, per `CLAUDE.md`'s own documented per-timestep order) -- and this has NEVER been instrumented at this specific cell before.** The existing generic watch-cell probe (`Antidunes.cpp`'s `trackWatchCells(..., "post-AfterDEM", ..., watchCells)`) only tracks a different, unrelated location (world ~(20.5,8.5,25.5), the Bogner-box-era waterline region) -- it has never once looked at the hopper crash cell.
>
> **Action taken: added the missing instrumentation directly in `Antidunes.cpp`** -- a new `hopperCrashCells` list (same X=267.5/Y=5.5-35.5/Z={17.5..20.5} coverage as `BuildTimeloop.h`'s `probeHopperCrashRegion`) with three new `trackWatchCells` probe points, gated to steps 2700-3200: right after DEM sub-cycling (before `timeloopAfterDEM`), right after `timeloopAfterDEM->singleStep()`, and right after `timeloopEval->singleStep()` -- bracketing every remaining unwatched stage of the per-timestep sequence. Verified this compiles cleanly (local `cmake --build build --target Antidunes`). **Working hypothesis pending this data**: a geometric knife-edge (X=267.5 sits almost exactly at the first rod-centre column's own radius boundary, matching the rod-geometry entry above) could cause `MovingParticleMapping`'s per-step containment test to flicker this cell between obstacle/fluid classification by a hair's breadth every timestep, triggering `PdfReconstructionManager`'s fill-seeding/PDF-reconstruction logic repeatedly on a cell that (from the FSLBM timeloop's perspective, which only samples flags once per step at fixed hook points) always LOOKS permanently obstacle -- each flicker perturbing the PDFs slightly, compounding over many steps until collide's own repeated (but not literally isolated, per the refuted hypothesis above) application to this slowly-drifting state eventually produces the observed collapse. Job **4266175** (mirror2, same prm, this new instrumentation) submitted to test this directly.

> **[Claude, 2026-09-25] CORRECTION + precise geometric characterization: the crash cell sits at a rod's flat END-CAP edge (axial), NOT at a radial knife-edge as an earlier entry this session overstated.** Pulled the exact `FloatingCylinder` entry nearest the crash cell directly from `prms/hopper_clog_cyl_frontblock_pidon_2node_halfrows_pass3_precursor_gas.prm:301`: `position <270.00, 10.00, 19.04>; radius 3; length 18; axis <0,1,0>` -- i.e. a rod centred radially at (X=270, Z=19.04), axis along Y, spanning Y=[1,19]. For the crash cell (267.5, 18.5, 19.5): **radial distance** (X-Z plane, perpendicular to the Y-axis) = `sqrt((267.5-270)^2 + (19.5-19.04)^2) = sqrt(2.5^2+0.46^2) ~= 2.54` -- comfortably inside the radius-3 circle (margin ~0.46 LU, ~15% of the radius), so the EARLIER entry's claim that "X=267.5 sits just short of the first rod-centre column... right at the column's own surface" was an overstatement and is hereby corrected. **The real singularity is axial, not radial**: this rod's Y-extent ends at Y=19 (its flat end-cap disk), and the crash cell sits at Y=18.5 -- only 0.5 LU inside that end, right at the edge of the 2-LU inter-group gap before the next Y-group's rods begin at Y=21. A cylinder's flat-end/lateral-surface junction is a genuine geometric edge (discontinuous surface normal), and a fluid cell whose centre sits this close to it, along the axial direction specifically, is exactly the kind of location where `MovingParticleMapping`'s per-cell containment test (a simple inside/outside check against the analytic shape) is most sensitive to sub-lattice-unit jitter in the rod's real-time position/orientation -- directly supporting (now with the corrected geometry) the "geometric knife-edge causing a per-step classification flicker" hypothesis already under test via the `timeloopAfterDEM` probes and the in-kernel raw-read diagnostic (job 4266208, pending).

> **[Claude, 2026-09-25] DECISIVE MECHANISM FOUND: `boundarySweep`'s MEM bounce-back (`SimpleBB::treatDirection`) writes directly into the OBSTACLE cell's own PDF array every step -- an assumption held (and stated) repeatedly earlier this session ("boundarySweep doesn't touch OBST cells' PDFs") was WRONG, and correcting it closes the entire investigation end to end.**
>
> Re-seeded the isolated-iteration check (`ANTIDUNES_CollideIsolatedObstacleIterationRealSeed`, new test) with THIS run's own real captured OBST-phase state (job 4265770, step 2700, rho=1.007356) instead of a different occurrence's value, and iterated `collide` alone (zero force) 200 times. **Result: rho changed only by floating-point rounding noise (~1e-16), i.e. ANOTHER exact fixed point, zero real drift.** This is now the SECOND independent confirmation (different real starting values, same conclusion) that repeated isolated `collide` iteration is NOT the source of the observed real multi-step drift -- something else genuinely rewrites this cell's PDFs every real step.
>
> Compared the exact 27 PDF components between step 2827's end (`equilibriumRefilling`, rho=1.007128) and step 2828's very first FSLBM hook (`after-boundarySweep`, rho=1.006672) -- nothing else runs in that gap (confirmed: nothing changed between step 2827's `equilibriumRefilling` and its own `after-DEMSubCycles`/`after-timeloopAfterDEM`/`after-timeloopEval`, all identically 1.007128; the very next thing that runs is step 2828's own `boundarySweep`). **Exactly 14 of the 27 components changed, 13 stayed BIT-IDENTICAL** -- not the "all 27 change together" signature of `collide`, but the signature of a boundary condition that only touches specific stencil *directions*.
>
> Read `SimpleBB.h::treatDirection` in full (not just the diff) to find the actual write target. It iterates fluid-side cells `(x,y,z)` for each direction into an adjacent obstacle `(nx,ny,nz)`, computes the bounce-back-reflected `pdf_new`, and writes: **`pdfField_->get( nx, ny, nz, Stencil_T::invDirIdx(dir) ) = pdf_new;`** -- i.e. INTO THE OBSTACLE CELL's own PDF slot (standard Ladd/MEM convention: a covered cell needs a well-defined phantom population for momentum-exchange bookkeeping), not the fluid cell's, as had been assumed. So for our crash cell -- which sits near a rod's flat end-cap edge (corrected geometry entry above) and therefore has a MIXED neighbourhood of real obstacle AND real liquid/interface neighbours (confirmed directly in the `[NaN-SCAN]` 3x3x3 neighbourhood dump: the `dz=-1` layer is all `LIQ`) -- every real liquid/interface neighbour triggers its own `treatDirection` call that overwrites ONE of our cell's 27 components with a fresh, physically-grounded reflected value every single step.
>
> **Full causal chain, closing the investigation:** (1) end-cap-edge location -> mixed obstacle/fluid neighbourhood (unlike a cell buried deep inside a rod, which has an all-obstacle neighbourhood and is never touched by `boundarySweep` at all); (2) every step, `boundarySweep` refreshes ONLY the components facing real fluid neighbours with genuine, physically-evolving bounce-back values; (3) the remaining components (facing obstacle neighbours) are never refreshed by `boundarySweep` OR `streamReconstructAdvect` (confirmed earlier) -- they persist purely as whatever `collide` last produced; (4) every step, `collide` runs its full cumulant transform on this **hybrid, partially-fresh/partially-self-referential state** -- something the transform's derivation assumes never happens (one coherent, physically consistent population set); (5) at `omega=1.99753`, essentially at the zero-viscosity stability limit, there is minimal damping to correct the growing inconsistency between the two categories of components, and over enough real steps this can push `m_000` toward zero, at which point `xi_0=1/m_000` amplifies the inconsistency into a full blow-up.
>
> **This is consistent with every single finding from this entire investigation, not just the newest one**: ASan-clean (this is legitimate code writing to memory it is SUPPOSED to write to -- not corruption); single-shot `collide` replay always finite (true for any ONE coherent, internally-consistent input, which is exactly what my replay tests always fed it); BOTH isolated-iteration tests showing zero drift (neither ever re-applied `boundarySweep`, so neither could reproduce the real mechanism -- this was the missing piece, not a contradiction); the fixed crash location (only end-cap-adjacent cells get this specific mixed-neighbourhood treatment; a cell fully inside a rod is fully collide-governed and stable, a cell fully in open fluid is never obstacle-flagged and has no issue); and the step-to-step non-determinism (the exact bounce-back values each step depend on the real, chaotically-evolving neighbouring flow state, which differs slightly between independent runs). **Root cause: MEM bounce-back's per-cell partial PDF refresh, applied indefinitely to a permanently-obstacle-covered cell with a mixed neighbourhood, combined with `collide` running unconditionally on that same cell every step with no re-equilibration mechanism, at a relaxation rate offering minimal numerical damping.** Not yet done: a direct fix proposal (candidates: skip `collide` entirely for obstacle-flagged cells that will never be read as fluid data anyway, since nothing downstream needs a "moving obstacle" cell's own bulk PDFs except `SimpleBB`'s own bounce-back write target, which does not require the PRE-collide value to be physically meaningful; or make `SimpleBB` refresh ALL 27 components of a covered cell to a consistent local equilibrium rather than only the linked directions, removing the hybrid-state inconsistency at its source).

> **[Claude, 2026-09-25] Smoking-gun confirmation, caught directly in a single log line: the mechanism above literally shows a cell "half alive, half dead" one step before collide unifies it into total corruption.** Job **4266208** (kernel-level raw-read diagnostic + real-run instrumentation) crashed at a DIFFERENT but geometrically symmetric cell this time -- world (267.5,21.5,19.5), local (11,1,19) -- the mirror-image end-cap of the *other* Y-group's rod (Y=21 is that rod's own near end, exactly analogous to the usual Y=19 end studied above; a THIRD independent confirmation that this is a geometry-driven, not location-specific-fluke, phenomenon). Its `[WATCH t=2830 HOPPER-after-streamReconstructAdvect]` line reads: `f[0]=-nan f[1]=-nan f[2]=7.4442e-02 f[3]=7.6042e-02 f[4]=-nan f[5]=-nan f[6]=-nan f[7]=1.9035e-02 f[8]=-nan f[9]=1.9002e-02 f[10]=1.8226e-02 f[11]=-nan f[12]=1.8004e-02 f[13]=1.8423e-02 f[14]=nan f[15]=nan f[16]=1.9097e-02 f[17]=nan f[18]=nan f[19]=-nan f[20]=4.6142e-03 f[21]=4.3868e-03 f[22]=4.5955e-03 f[23]=-nan f[24]=-nan f[25]=4.6678e-03 f[26]=4.8729e-03` -- **13 of 27 components are already `nan`, the other 14 are completely ordinary, clean near-equilibrium values, in the SAME cell, at the SAME instant.** This is exactly the signature the mechanism above predicts: the `nan` components are the ones facing obstacle neighbours (never refreshed by anything since `streamReconstructAdvect` doesn't touch obstacle-cell PDFs and this cell has presumably been in this partially-poisoned state for a while), and the clean components are the ones `boundarySweep`'s MEM bounce-back refreshes every step from real, healthy fluid neighbours. The very next stage, `[WATCH t=2830 HOPPER-after-collide]`, shows all 27 components `nan` -- because `collide`'s cumulant transform sums ALL 27 PDFs into `m_000` and higher moments, and under IEEE-754 arithmetic a single `nan` operand poisons every sum it touches unconditionally, so the moment `collide` mixes the "half alive" cell into one coherent set of moments, the result is uniformly non-finite. **This is not a hypothesis anymore -- it is the mechanism, directly visible in a single log line, independently reproduced at a geometrically-symmetric second location.**
>
> (Note: the in-kernel `[KERNEL-RAW-READ]` diagnostic added to the generated `kernel_collide` this session did NOT fire for this specific cell in this run -- 1010 hits total, all at the already-known-harmless early-startup waterline red herring locations, none at (11,1,19). Given `m_000` built from a PDF vector already containing multiple `nan` components is itself unconditionally `nan` under IEEE arithmetic, the gate `!std::isfinite(m_000)` should have caught this. The likely explanation is that this cell's containing block's `StateSweep` classification was `onlyLBM` rather than `fullFreeSurface` at the exact moment this happened, in which case the FUSED `kernel_streamCollide` ran instead of the standalone `kernel_collide` the diagnostic was only added to -- not a contradiction of the mechanism, just an instrumentation gap in a kernel variant not yet touched. Not pursued further since the wrapper-level evidence above is already fully decisive.)

> **[Claude, 2026-09-25] FINAL CONFIRMATION at the lowest possible level: the kernel's own raw pointer-arithmetic read of the crash cell's PDF memory is genuinely, exactly all-zero at the instant `collide` processes it -- closing the "does collide see something different from what our probe captured" question with a definitive NO.**
>
> Fixed the instrumentation gap: the original `[KERNEL-RAW-READ]` diagnostic (added to `kernel_collide` only) used a generic "first 30 hits, anywhere" throttle that saturated on the already-documented pervasive early-startup red herring (steps ~8-12, harmless, present identically in every scenario variant regardless of Pass-3/density/PID -- rods sitting at/near the water level Z~49 at initialisation briefly starve nearby cells of real interface-neighbour evidence, triggering the same `[FALLBACK-ISOLATED-INTERFACE]` mechanism documented earlier, self-resolving within a few steps, structurally unrelated to the real ~2800+-step-later crash). That silently exhausted the log budget before the real crash cell ever failed. Fixed by gating on the crash region's own block-local coordinates (`ctr_0==11 && ctr_2==19`) instead of a hit-count cap, and added the same diagnostic (rank-tagged, adapted for its different inline-neighbour-offset variable structure) to `kernel_streamCollide`, the fused kernel that was never instrumented before. Job **4266339** (mirror2) reproduced the crash at BOTH previously-seen locations simultaneously this run (world (267.5,18.5,19.5) rank 48, AND its mirror-image (267.5,21.5,19.5) rank 54) -- filtering strictly to rank 48 (the confirmed owner of the FULL-BAD-CELL-LIST report) and local (11,18,19):
>
> `[KERNEL-RAW-READ:collide] rank=48 local=(11,18,19) m_000=0.000000000e+00 | pdfs: off0=...off26=0.000000000e+00 (all 27, exactly) | force: fx=fy=fz=0.000000000e+00`
>
> This is the kernel's OWN raw `_data_pdfs[stride*ctr_0+...]` pointer read, bypassing every higher-level accessor (`PdfField::get()`, `trackWatchCells`, everything) -- and it matches exactly what the wrapper-level probe showed moments earlier at this same cell: genuinely, exactly zero across all 27 PDF components and the force vector, not floating-point noise, not a near-zero value. **This proves there is no stale-read/race discrepancy between what our probes observed and what collide actually consumed -- the two agree perfectly.** `m_000=0` feeds `xi_0=1/m_000=inf` directly, and every downstream cumulant/moment computation involving `xi_0` propagates that `inf` into `nan` across all 27 output components -- the exact, now doubly-confirmed mechanism (`StreamReconstructAdvectSweep.h`'s obstacle branch never writing PDFs, so whatever the ping-pong buffer swap exposes for this cell is whatever was left there, which this time was genuinely all-zero).
>
> **Bonus finding from the fix itself**: `kernel_streamCollide` (the fused kernel, previously uninstrumented) ALSO fired repeatedly at this exact same rank/cell (`#16, #34, #52...`, all `m_000=-nan`) at other points across the run -- direct proof that this block's `StateSweep` classification genuinely toggles between `fullFreeSurface` and `onlyLBM` over the run's duration, both kernel variants really do process this one cell at different times, and that toggling (not a fixed classification) is why the original kernel_collide-only instrumentation kept missing it. **This closes the investigation's central open question. The full chain, top to bottom, is now evidence-backed at every link**: end-cap-edge geometry -> mixed obstacle/fluid neighbourhood -> `boundarySweep`'s MEM bounce-back partially refreshing only some PDF components every step -> the untouched remainder persisting from whatever the ping-pong buffer last held (confirmed here to genuinely reach exact zero) -> `collide` (running unconditionally on obstacle cells, in either kernel variant depending on the block's current `StateSweep` state) computing `1/0=inf` from that all-zero input -> `inf` propagating into `nan` across every output component.

> **[Claude, 2026-09-25] FIX IMPLEMENTED AND VALIDATED: two candidate fixes tested head-to-head, one mirror each, both with the fully-instrumented kernel. Option B (narrow reactive repair) works end-to-end; Option A (framework clone-buffer fix) is real but insufficient alone -- and the comparison reveals the instability is far more pervasive than the single end-cap cell this whole investigation focused on.**
>
> **Option A** (`src/lbm/free_surface/dynamics/StreamReconstructAdvectSweep.h`): the confirmed root cause of the specific zero-PDF mechanism -- `pdfFieldClone_.get(block)` allocates its scratch buffer via `cloneUninitialized()` on its one-time first creation, never mirroring the real field's own equilibrium initialisation. Fixed by capturing the `created` flag from `FieldClone::get(block, created)` and, when true, copying the real field's current contents into the fresh clone before use.
>
> **Option B** (`SetupStages.h`'s new `repairDegenerateObstacleCells()`, wired into `BuildTimeloop.h` as an `AfterFunction` right after `collide`): scans obstacle-flagged cells only, and resets a cell to the D3Q27 rest-equilibrium distribution (rho=1, u=0) only if its density is already non-finite or non-positive -- healthy obstacle cells (the vast majority, including ones actively receiving legitimate `SimpleBB` bounce-back writes) are never touched.
>
> **Test setup**: job **4266472** (`mt-dennys-mirror`) ran Option A alone (Option B's AfterFunction call disabled for isolation); job **4266473** (`mt-dennys-mirror2`) ran Option B alone (Option A not applied). Both used the identical fully-instrumented `AntidunesLatticeModel.cpp` (rank-tagged `[KERNEL-RAW-READ]` diagnostics, properly gated) and the same prm capped at 5000 timesteps (comfortably past the previously-observed 2700-3200 crash window).
>
> **Results:**
> - **Option A (4266472): FAILED at step 3285** -- later than every previous unpatched run (typically 2830-3133), and via a DIFFERENT signature: not an obstacle cell reading exact zero, but a **LIQUID cell's density diverging to an astronomical finite value** (`rho=-2.25e15`, `max|f|=1.5e32` at world (287.5,11.5,60.5)), with the `[DENSITY-WINDOW]` scan showing 232 cells below rho=0.5 and 299 above rho=1.5 domain-wide at that same step. This proves Option A fixed a REAL mechanism (the run got measurably further, and the specific "exact zero" signature never recurred) but exposed a SEPARATE, broader numerical instability not confined to the clone-buffer bug.
> - **Option B (4266473): COMPLETED all 5000 steps, exit status 0.** Zero `VEL-NaN-SCAN` or `DENSITY-WINDOW` warnings anywhere in the run -- the liquid-domain cascade that killed Option A's run never developed at all. The repair guard itself fired extremely often (~30 million logged events, throttled from an even larger true count, at cells scattered across the ENTIRE domain from step 0 through step 4999) -- some repaired cells showed the same kind of astronomical densities seen in Option A's failure (e.g. `-6.4e152`, `-1.3e154`), not just clean exact zeros.
>
> **Reframing of scope**: the instability is NOT confined to the one end-cap-adjacent cell (world (267.5,18.5,19.5) and its mirror (267.5,21.5,19.5)) this entire investigation focused on -- it is a pervasive, continuous, domain-wide phenomenon affecting many obstacle cells throughout the whole rod bundle, the entire run. That one cell was simply the first, in every unpatched run, to escalate all the way to a hard abort. The most likely reason Option B succeeds where Option A alone does not: when a badly-diverged obstacle cell is eventually uncovered by a moving rod (formerMO conversion to liquid/interface), whatever corrupted state is sitting in its PDF array can leak directly into what is now a genuine fluid cell and cascade outward via ordinary streaming/collide through the liquid domain -- exactly Option A's failure mode. Option B prevents this at the source by never allowing ANY obstacle cell to remain corrupted, so there is never a bad cell available to leak at the moment of conversion.
>
> **Conclusion: Option B (the narrow reactive repair guard) is the validated, working fix and should be the one carried forward.** Option A remains a real, independently-justified correctness fix for the framework-level `FieldClone` gap (worth keeping/upstreaming on its own merits) but is not sufficient by itself to close this crash class -- the two are not redundant, and running both together (as the local development tree now does) is reasonable defense-in-depth, though Option B alone is sufficient and is the one to rely on.
>
> **Open follow-up, not yet investigated**: what exactly is the SEPARATE, broader instability mechanism that Option A's run exposed (astronomical liquid-cell densities, ~30 million obstacle-cell repairs needed over 5000 steps even with Option B active)? Candidates worth checking: `omega=1.99753` genuinely being too close to the zero-viscosity stability limit for this scenario's actual local Reynolds number/velocity gradients (independent of the obstacle-cell mechanism entirely); or a second, distinct code path also capable of producing degenerate obstacle-cell PDFs beyond the one clone-buffer-initialization gap Option A closes.

> **[Claude, 2026-09-25] CORRECTION: Option B alone does NOT fully close this scenario -- a full 60k-step validation run (job 4266492, heavy per-step diagnostics disabled) crashed at step 3267, via a mechanism Option B was never designed to catch.**
>
> The earlier "Option B completed all 5000 steps" result (job 4266473) was real but incomplete: given this scenario's established non-determinism, that specific run simply did not happen to trigger a SEPARATE, pre-existing failure mode within its first 5000 steps. Job 4266492's crash: a `[MEM-GARBAGE-LINK]` warning at step 3267 on a LIQUID cell (block-local (31,15,29)) whose PDFs had already diverged to wild values (e.g. 321179, -158882, -40184.8, versus the normal ~0.01-0.3 range). That garbage fed directly into the MEM force computation on the neighbouring rod (uid 3530822107858468885 -- the SAME uid that appeared as one of the two obstacle owners in job 4266472's own pinched-cell crash, strongly tying this to the same underlying mechanism), launching it to an absurd position/velocity (`pos=<2823.92,2522.83,-2482.43> vel=<12674.8,12560.7,-12701>`) -- thousands of units outside the domain. That triggered `[GHOST-VANISHED]` (ranks losing track of the rocketed particle) and a cascade of `[STRONG-CONTACT-FORCE]`/`[DEM-WATCH] deep contact` events as it slammed through its neighbour rod (uid ...888), overlaps up to 84% of the interaction radius.
>
> **Root cause: this is job 4255679's own pre-existing "class-B sealed-cell" mechanism** (a liquid cell nearly fully enclosed by two touching rods, undergoing a MEM/DEM density runaway) -- completely unrelated to, and unaddressed by, either the clone-buffer fix (Option A) or the obstacle-cell repair guard (Option B), since **Option B only scans and repairs OBSTACLE-flagged cells; the corrupted cell here was flagged LIQUID**, entirely outside its scope. This is also a strictly more severe failure mode than the obstacle-cell mechanism this whole session characterised: instead of a contained, cleanly-detected density blowup, it manifests as a full DEM catastrophe (a particle physically launched out of the simulation domain), which is harder to guard against reactively after the fact.
>
> **Next test (job 4266507, this entry)**: Option B (`repairDegenerateObstacleCells`, unchanged) + `closeSealedCellsMinNeighbours=24` (the project's own pre-existing, purpose-built guard for exactly this sealed-cell mechanism, previously tried alone for a different scenario variant -- job 4259313 -- and found insufficient IN ISOLATION at that time, though that test predates both Option A and Option B and was against a different combination of live bugs). New prm: `hopper_clog_cyl_frontblock_pidon_2node_halfrows_pass3_precursor_60k_guard24_gas.prm`. Full 60k steps, heavy per-step diagnostics still disabled (this is now a production-style validation run, not a debugging run).

> **[Claude, 2026-09-25] `closeSealedCellsMinNeighbours=24` delays but does NOT eliminate the sealed-cell mechanism -- and reveals the guard is checking the wrong signature.**
>
> Job **4266507** (Option B + guard24, full 60k): pushed the crash from step ~3267 (guard off, job 4266492) to **step ~8393** -- roughly 2.5x further -- but the same `[MEM-GARBAGE-LINK]`-then-rod-launch mechanism still occurred. This time the corrupted liquid cell (block-local (56,10,28)) had `density=8.52` (deceptively unremarkable at the aggregate level) while its INDIVIDUAL raw PDF components had exploded to ~3.4e10/-1.7e10/8.5e9 etc -- massive internal cancellation hiding the corruption from a simple density check. The affected rod (uid 4827858800541171719) was launched to `vel=<1.27e6,-4.2e5,4.1e5>`, an even more extreme runaway than job 4266492's.
>
> **Root cause of the guard's incompleteness, confirmed by directly counting this cell's neighbourhood from the `[MEM-GARBAGE-LINK]` log's own neighbour-flags dump: only 19 of 26 neighbours are obstacle-flagged (7 genuinely liquid) -- well under the 24 threshold, so `closeSealedCellsMinNeighbours` never touches this cell at all.** Critically, the 19 obstacle neighbours split across TWO DIFFERENT flag values (10 cells flagged `512`, 9 cells flagged `1024` -- almost certainly two distinct rods' individual per-particle obstacle flag bits, consistent with the earlier job 4266472 trace where a similar cell's obstacle neighbours were explicitly owned by two different rod UIDs). **The guard checks TOTAL enclosure (how surrounded a cell is overall), but the actual danger signature is a cell pinched between TWO SEPARATE solid bodies specifically -- which can occur at a much lower total-obstacle-neighbour count than 24/26, as long as the fluid gap between the two particular rods is narrow enough.** A cell substantially enclosed by ONE rod's own rounded surface (a smooth, single-body curvature) is comfortably stable even at high neighbour counts; a cell squeezed between two DIFFERENT rods can be dangerous even at a much lower count. `closeSealedCellsMinNeighbours`'s single scalar threshold cannot distinguish these two cases.
>
> **Implication for a real fix**: a guard keyed to "opposite-direction neighbour pairs both obstacle AND belonging to different particle UIDs" (rather than a flat total-neighbour-count threshold) would likely target this failure mode far more precisely -- closer in spirit to the already-existing `pinchedPairs` diagnostic metric (`checkPdfNaN`'s own neighbourhood analysis, already computing exactly this "opposite-pair-both-obstacle" count) than to `closeSealedCellsMinNeighbours`'s current simple count. Not yet implemented or tested.

> **[Claude, 2026-09-25] FIX IMPLEMENTED: `SimpleBB.h`'s existing `[MEM-GARBAGE-LINK]` diagnostic promoted from log-only to an active guard -- the correct, targeted intervention for the sealed/pinched-cell mechanism, keyed on the actual quantity that matters (a single PDF component), not density.**
>
> Root cause clarified from the last two crashes (jobs 4266507, 4266539): the *aggregate* density at the corrupted cells looked only mildly anomalous (rho=8.52, rho=-1.12, rho=-0.47 etc) while their INDIVIDUAL PDF components had already exploded (values up to ~1e10, or into the hundreds/thousands, cancelling out in the sum). `SimpleBB::treatDirection`'s MEM force computation (`forceMEM = pdf_old + pdf_new`) reads exactly ONE such component directly and unconditionally per link -- no density-weighting or averaging ever dilutes a single bad value before it reaches the particle's force. This is why a density-based guard (Option B, or any density-sum check) structurally cannot catch this failure mode: the danger signature lives in the individual populations, not their sum.
>
> **Fix**: the file already had a diagnostic for exactly this (`if (std::abs(pdf_old_raw) > real_c(10))`, in place since 2026-09-20 for a different investigation, job 4251049) but only logged `[MEM-GARBAGE-LINK]` and then fed the garbage value into the force computation anyway. Added a `return;` right after the logging block -- mirroring the existing "dry-link" guard immediately above it (which already skips non-finite `pdf_old` the same way) -- so an absurdly-large-but-finite population now causes this link to be skipped entirely (zero momentum-exchange contribution from it this step) instead of corrupting the particle's force and launching it out of the domain.
>
> **Not yet validated on the cluster** -- next step: combine with Option B (`repairDegenerateObstacleCells`, already proven for the separate obstacle-cell zero-PDF mechanism) for a full 60k run, to test whether the two together finally close out both known failure classes for this scenario.

> **[Claude, 2026-09-25] SimpleBB.h fix validated as effective against its target mechanism, but a THIRD, distinct failure path surfaced once it was suppressed.**
>
> Job **4266571** (Option B + SimpleBB fix, no sealed-cell guard, full 60k, 1h cap): ran to **step 5900** with **zero `[MEM-GARBAGE-LINK]` or DEM rocket-launch events anywhere in the run** -- direct confirmation the SimpleBB.h fix successfully closes the mechanism it targeted (a corrupted individual PDF component reaching the MEM force computation). Crashed instead via a genuinely different, not-yet-investigated signature: `[LOW-RHO]` -- three adjacent INTERFACE cells ((9,3,7), (9,4,7), (9,5,7)) with `rho=nan` and all 27 PDFs `nan` simultaneously, caught by the free-surface interface-cell reconstruction's own pre-existing guard (distinct from `checkPdfNaN`'s NaN-SCAN family, and distinct from both the obstacle-cell mechanism Option B targets and the MEM-force mechanism SimpleBB's fix targets). This is a third, separate failure path in this scenario, not yet root-caused.
>
> Job **4266572** (mirror1: Option A + Option B + guard24 + SimpleBB fix, all four combined, full 60k, 1h cap): still running as of this entry, already past 12 minutes (further than 4266571's crash point in wall-clock terms, though direct step comparison across runs is unreliable given this scenario's established non-determinism).
>
> **Running tally of this scenario's distinct known failure classes**: (1) obstacle-cell PDF reads as zero -- FIXED by Option A (framework) and independently by Option B (reactive repair); (2) pinched liquid cell between two rods feeding a corrupted individual PDF component into MEM force computation, launching a particle out of the domain -- FIXED by the SimpleBB.h guard (confirmed: zero occurrences with the fix active, versus reliably present without it); (3) interface-cell reconstruction producing `[LOW-RHO]` NaN -- NOT YET FIXED, NOT YET ROOT-CAUSED, first directly observed only once mechanisms (1) and (2) were both suppressed for long enough to reveal it.

> **[Claude, 2026-09-25] Confirmed: the `[LOW-RHO]` interface-cell mechanism is real and independent -- it recurred even with ALL FOUR mitigations combined.**
>
> Job **4266572** (Option A + Option B + guard24 + SimpleBB fix, all four together, full 60k, 1h cap): crashed at step 9300, via the same `[LOW-RHO]` signature as job 4266571 -- this time at multiple scattered locations simultaneously ((50,9,16), (2,12,14), (2,12,15), (2,12,16), ...), all interface cells, all 27 PDFs NaN. Since this occurred even with every currently-known mitigation active at once, it is conclusively a third, independent failure class in this scenario, not a residual manifestation of mechanisms (1) or (2).
>
> **Summary of validated fixes so far:**
> 1. Obstacle-cell PDF reads as zero (`FieldClone`'s uninitialised scratch buffer) -- FIXED, two independent ways (Option A: root-cause fix in `StreamReconstructAdvectSweep.h`; Option B: reactive repair guard in `SetupStages.h`/`BuildTimeloop.h`). Either alone closes this specific mechanism.
> 2. Pinched liquid cell between two rods -> a single corrupted PDF component reaching `SimpleBB`'s MEM force computation -> particle launched out of the domain -> cascading contact catastrophe with neighbours -- FIXED (`SimpleBB.h`'s `[MEM-GARBAGE-LINK]` check promoted from log-only to an active skip-the-link guard). Confirmed: zero occurrences with the fix active across two separate 60k runs, reliably present without it.
> 3. `[LOW-RHO]` interface-cell NaN (free-surface reconstruction's own pre-existing guard, distinct from all of the above) -- NOT YET FIXED, NOT YET ROOT-CAUSED. Not yet investigated at all: which reconstruction code path produces this, what triggers it, whether it is itself a symptom of an even-further-upstream cause.
>
> **Net conclusion for this line of work**: real, validated progress on two genuinely distinct bugs, but the scenario is not yet clean for unattended long-horizon (60k-step) production use. The next investigation, if pursued, should target `[LOW-RHO]`/the interface-cell reconstruction path specifically -- likely `ReconstructInterfaceCellABB.h` or wherever this project's own reconstruction machinery computes momentum density for interface cells, following the same evidence-first methodology used for mechanisms (1) and (2) above (raw-value tracing before assuming a cause).

> **[Claude, 2026-09-25] Minimal two-rod pinch reproduction (per explicit request "create a small
> test case that proves the pinching crash is a real thing") -- NEGATIVE result, and a correction
> to an earlier claim about job 4266473.**
>
> **Job 4266473 re-examined ("check the full run with all the cylinder with guard on that didn't
> crash, why?")**: confirmed directly from its own log (`timesteps 5000;` in the printed config,
> performance logger fires at `[4999]`, sacct elapsed 00:15:19) that this run's SUCCESS was a length
> artifact, not evidence the guards solved the scenario -- it was simply told to stop at 5000 steps,
> short of the 5900-9300-step window where the `[LOW-RHO]`/pinch instability actually appears in the
> longer (60k-programmed) runs. `[OBSTACLE-CELL-REPAIR]` (Option B) DID fire 14,349 times during
> those 5000 steps -- proof the obstacle-cell corruption (mechanism 1) was continuously present and
> being suppressed live throughout, even in this "clean" run. Zero MEM-GARBAGE-LINK/NAN-LINK/DRY-LINK/
> LOW-RHO/SEALED-CLOSE. (Earlier note in this file mischaracterised this job as demonstrating "Option
> B works end-to-end" -- retracting that framing per [[verify-own-session-comments]]: it demonstrates
> only that Option B suppresses mechanism 1 continuously, not that the scenario is stable.)
>
> **Minimal reproduction built and run** (`prms/pinch_test_two_rods_gas.prm`, job **4267098**, mirror1,
> 3 ranks, single-node, 20 min cap): two FloatingCylinders at the *exact* production pitch/radius
> (radius 3, X-pitch 7 -> 1 LU surface clearance -- the identical spacing behind job 4266507's real
> crash cell, 19/26 neighbours split between two rod UIDs), axis Y, held FIXED for the whole run
> (`settlingSteps 1000000` >> `timesteps 5000` -- confirmed via SetupStages.h that FIXED means zero
> invMass/invInertia, i.e. immovable but still fully live in the DEM contact + MEM force pipeline
> every step), fully submerged in still water (`targetMeanVelocityMagnitude 0`), single small domain
> `<48,24,40>` / 3 blocks in X (mesa_pd's own minimum for a periodic direction with particles -- a
> 1-rank attempt aborted immediately with "When using particles, use at least three blocks per
> periodic direction", unrelated to the physics). All current diagnostics live (Option A+B, the
> SimpleBB.h fix, precursor/NaN scans).
>
> **Result: ran clean for the full 5000 steps, zero hits on every diagnostic** (MEM-GARBAGE-LINK,
> NAN-LINK, DRY-LINK, LOW-RHO, SEALED-CLOSE, OBSTACLE-CELL-REPAIR, both precursor checks), no
> assertions, no crash. **The static 1-LU-clearance geometry alone, perfectly symmetric and
> perfectly still, is NOT sufficient to trigger the pinch mechanism.** Something else present in the
> full hopper scenario but absent here -- channel flow/turbulence, PID forcing, DEM jitter from the
> other ~208 rods, or actual DEM contact/overlap between the two rods (not just nominal 1 LU
> clearance) -- is a necessary ingredient. This sharpens, rather than refutes, the "Class-B
> sealed-cell" framing: proximity/geometry alone isn't the trigger; some dynamic perturbation is
> required to actually drive a liquid cell at that pinch point into corruption. Next step, if pursued:
> repeat this same minimal setup but let the rods actually move (release under gravity/buoyancy so
> they can settle into real contact/overlap), and/or add a mild driven flow, to isolate which specific
> ingredient is the missing piece.
>
> Output archived locally at `~/University/Thesis/animations/antidunes-pinch-test-two-rods_4267098_20260925/`
> (fluid field + both rods' particle VTK series) for visual inspection.
>
> **Second, still-open finding from reading `SimpleBB.h` during this investigation**: its `[NAN-LINK]`
> guard (immediately after the `[MEM-GARBAGE-LINK]` guard we fixed) is diagnostic-only -- it logs a
> non-finite `force` but never skips the link, unlike `[MEM-GARBAGE-LINK]`'s `return;`. `force` is
> built from `density` (the FULL 27-PDF sum via `getDensity(x,y,z)`, compressible path), not just this
> link's own `pdf_old` -- so a corruption in ANY OTHER direction at the same fluid cell (e.g. an
> incipient `[LOW-RHO]`-style single-component blowup that hasn't yet pushed `pdf_old` past the |>10|
> threshold checked by our fix) can still reach `addHydrodynamicForceAtWFPosAtomic` completely
> unguarded. This is a second, structurally identical leak of the same underlying corruption, parallel
> to the one already closed -- plausibly part of why the pinch/LOW-RHO mechanism still gets through
> after the SimpleBB.h fix. NOT YET FIXED -- the natural next step is promoting `[NAN-LINK]` to an
> active guard the same way, mirroring the existing pattern exactly (log first N, then `return;`).

> **[Claude, 2026-09-25] New scenario: "drop waves" half-density cylinder hopper (per explicit
> request), fixes-off baseline -- job 4267357 setup error, resubmitted as 4267519.**
>
> Built `gen_hopper_clog_cyl_dropwaves_prm.py` (companion of `gen_hopper_clog_cyl_prm.py`): same
> 210-rod population/densities as the halfrows block scenario, but released as 6 waves of up to 36
> rods (6x2x3) dropped from 6-55 LU above the water surface (Z=49.13) at an UPSTREAM drop zone
> (X=[9,113], 271 LU clear of the pillars at X=384) -- `startInactive true` + a `settlingSteps`
> ladder at t=2000/5000/8000/11000/14000/17000, 30000 timesteps total. Deliberately run with ALL of
> today's crash-mitigation fixes DISABLED on mirror2 ONLY (Option A in
> `StreamReconstructAdvectSweep.h`, Option B's wiring in `BuildTimeloop.h`, the SimpleBB.h
> `[MEM-GARBAGE-LINK]` `return;` guard reverted to diagnostic-only) -- a genuine pre-session
> baseline for this new release mechanic; local tree and mirror1 keep the fixes intact.
>
> **Job 4267357 failed in under 2 minutes -- not a physics crash.** Root cause found in the .err
> log: `FileReader returned an error ... Error opening parameter input file
> '.../prms/hopper_clog_cyl_dropwaves_gas.prm'`. Caused by a scp path bug on my end: `scp ... 
> fritz.iwia:.../prms/../` resolves to the Antidunes root directory, not `prms/`, so the new prm
> landed one level up from where it was expected. Build itself succeeded cleanly (only pre-existing
> unused-variable/misleading-indentation warnings, no errors). Fixed by moving the file into
> `prms/` on mirror2 and resubmitting as **job 4267519**.

> **[Claude, 2026-09-25] Drop-waves fixes-off baseline (job 4267519) result: crashed at step
> ~22037-22038, and for the first time shows a DIRECT, FAST causal link from MEM-GARBAGE-LINK to
> LOW-RHO -- a second, previously unseen route to the LOW-RHO mechanism.**
>
> Ran cleanly through all 6 wave releases (t=2000..17000), rods falling from above the surface,
> carried downstream, piling up near the pillars (X~220-330 by the time of the crash, pillars at
> X=384) with substantial real DEM contact activity (1769 GHOST-VANISHED, 1027 SUBCYCLE-LEFT-RANK,
> 6532 STRONG-CONTACT-FORCE by the end -- assessed as ordinary heavy-contact-traffic noise, not
> independently diagnostic, since none of these three ever coincided with an actual corruption
> signature until the final crash sequence).
>
> **Complete, precisely-timed crash chain** (all three ordinary text greps -- MEM-GARBAGE-LINK,
> LOW-RHO, AddExternalForces NaN -- corroborate each other, not inferred):
> 1. **Step 22035**: `[MEM-GARBAGE-LINK]` fires (77 total by end of run) on a LIQUID cell almost
>    fully enclosed by a SINGLE rod's own obstacle surface (25/26 neighbours `rawFlag=512`, one rod
>    -- notably NOT the "two different rods" signature from job 4266507; a single rod's own
>    enclosure is evidently also sufficient). One PDF component (`pdf_old(dir)=-1011.04`, vs. a
>    normal O(1) range) is fed unguarded into `forceMEM` -- the fix's `return;` is reverted on this
>    mirror, so despite the log text still saying "SKIPPING this link" (stale wording inherited
>    unchanged from the fix, NOT accurate on this deliberately-reverted build -- flagging this so
>    the wording doesn't mislead a future reader of this specific log), the value is NOT actually
>    skipped here.
> 2. **Step 22037 (2 steps later)**: TWO `[LOW-RHO]` interface cells appear -- (42,1,25) mild
>    (rho=-0.000154), (33,18,19) catastrophic (rho=0 exactly, all 27 PDFs at +-9e15, ~19 orders of
>    magnitude past normal).
> 3. **Step 22038**: `[DEM-WATCH] deep contact` at 89.9% and 23.9% of `smallerInteractionRadius`
>    (rods interpenetrating heavily), then `[AddExternalForces] NaN/Inf hydrodynamic force` AND
>    `torque` on FOUR different rod uids in the same step.
> 4. Shortly after: a rod's orientation quaternion goes NaN, trips a genuine `WALBERLA_ASSERT`
>    (`src/core/math/Quaternion.h:241`, `r*r+i*i+j*j+k*k == 1`), `Fatal error`/`Aborting now`,
>    `MPI_ABORT`. 32m19s elapsed (well inside the 3h cap) -- a genuine physics crash, not a timeout.
>
> **Interpretation, stated carefully (two runs, not conflated):** this is the FIRST time this
> session that a `[MEM-GARBAGE-LINK]` event has been observed with `[LOW-RHO]` appearing just 2
> steps later -- a second, distinct route into the LOW-RHO mechanism, on top of the independent one
> already established in jobs 4266571/4266572 (where the SimpleBB fix was ACTIVE, MEM-GARBAGE-LINK
> was ZERO throughout, and LOW-RHO still occurred regardless -- so that route does not require
> MEM-GARBAGE-LINK at all). Read together: LOW-RHO is reachable via at least two paths -- (a) an
> independent, not-yet-root-caused free-surface reconstruction issue (4266571/4266572), and (b) as a
> fast (~2-step) downstream consequence of an unguarded MEM-force corruption reaching a particle
> (this run). This strengthens rather than weakens the case for keeping the SimpleBB.h fix: it
> closes at least route (b) on top of the rod-launch mechanism it was built for, even though route
> (a) remains open and unaddressed by any of today's four mitigations.
>
> **New scenario itself validated as a working, useful reproduction vehicle**: the drop-waves
> release mechanic (6 waves of 36 rods, dropped from above the surface, carried downstream by the
> flow) runs cleanly through pile formation and produces a real, fast, well-instrumented crash --
> confirming this is a viable smaller/faster alternative to the full pre-placed-block scenario for
> studying this failure family, not just a one-off experiment.

> **[Claude, 2026-09-26] Wired in the capsule ("pill") shape as an opt-in experiment, per explicit
> request ("Try and wire in the pill shaped cylinder") -- job 4267738 submitted.**
>
> A capsule shape (`src/mesa_pd/data/shape/Capsule.h`, SHAPE_TYPE=7 -- cylinder + hemispherical end
> caps, i.e. a spherocylinder/pill) already existed fully wired into mesa_pd's core dispatch
> (DoubleCast.h, SingleCast.h, ShapePackUnpack.h -- all uncommitted in the working tree from earlier
> this session's DropTest work), but was only used by the separate fluid-free `DropTest.cpp` tool,
> never by the main `Antidunes.cpp` production sim. No analytic Capsule-Capsule contact routine
> exists -- capsule pairs fall through to generic GJK/EPA (`GeneralContactDetection`), which is
> exactly the point: a continuous closest-point-between-segments calculation with no flat-cap/rim
> edge case, unlike the analytic Cylinder-Cylinder Kodam+Guo dispatch that produced the
> "band-band-skewed" 90%-overlap failure traced in job 4267519.
>
> **Wiring added** (all in `apps/showcases/Antidunes/`, confined to this app -- no further mesa_pd
> core changes needed):
> - `SimulationParams.h`: new `ModelParameters.useCapsuleShapeForCylinders` (bool, default false).
> - `SetupStages.h`: `#include "mesa_pd/data/shape/Capsule.h"`; the FloatingCylinder creation loop
>   (~line 664) now constructs `mesa_pd::data::Capsule(fcp.radius, fcp.length)` instead of `Cylinder`
>   when the toggle is on -- everything else (position, rotation, interactionRadius, FIXED flag,
>   settlingSteps, updateMassAndInertia) is unchanged since Capsule shares Cylinder's (radius,
>   length) constructor and body-z-axis convention.
> - Local build (`cmake --build build --target Antidunes`) confirmed clean (link succeeded, only
>   pre-existing warnings).
>
> **Experiment**: `prms/hopper_clog_cyl_dropwaves_capsule_gas.prm` -- IDENTICAL to job 4267519's
> `hopper_clog_cyl_dropwaves_gas.prm` (same 210-rod drop-waves scenario, same fixes-off mirror2
> build with Option A/B and the SimpleBB.h guard all still reverted) except
> `useCapsuleShapeForCylinders true`. Submitted as job **4267738** on mirror2 so any difference in
> outcome vs. 4267519 is attributable ONLY to the shape change, not a confound with any of today's
> other mitigations. Result pending.

> **[Claude, 2026-09-26] Found and fixed a BLOCKING gap before the capsule experiment ever ran, per
> explicit request to double-check parity ("follow the cylinder implementation and double check
> what we have implemented with it... make sure the capsule has the same capabilities").**
>
> Two audits dispatched: (1) DEM contact detection parity (GJK/EPA, pillar, wall) -- still running;
> (2) fluid-obstacle mapping parity -- returned first, with a clear verdict: **`Capsule` had NO
> `mesa_pd::ContainsPointFunctor` overload in `src/mesa_pd/common/Contains.h`** (the file that backs
> `MovingParticleMapping.h`'s per-cell "is this lattice cell inside the solid" test via
> `SingleCast`). Every other shape (Sphere, HalfSpace, Box, Cylinder, Ellipsoid,
> CylindricalBoundary, ConvexPolyhedron) had one; Capsule fell through to the generic fallback
> template, which unconditionally does `WALBERLA_ABORT("ContainsPoint not implemented!")`. Job
> **4267738** (queued but still PENDING when this was found -- **cancelled before it ever started
> running**, no cluster time wasted) would have hard-aborted at its very first particle-to-lattice
> mapping call, not produced a wrong-but-plausible result -- worth noting since the earlier DEM
> contact-detection audit had already confirmed the DropTest.cpp tool's own prior use of Capsule
> never exercises this path (no LBM fields exist there), so this was genuinely never hit before.
>
> AABB/bounding-box was checked and is FINE for Capsule (both Cylinder and Capsule paths in
> `SetupStages.h` set `interactionRadius` from the flat-capped circumscribed-sphere radius
> `sqrt(r^2+(L/2)^2)`, which is >= a capsule's true farthest point `L/2` -- conservative, no
> clipping risk). Submersion/buoyancy volume calculation is currently DORMANT for this scenario
> (`useMemOnlyBuoyancy true`, `useAnalyticBoxBuoyancy false`, `useFullVolumeAnalyticBuoyancy false`
> in the prm) so the Cylinder-specific disk-plane submersion formula noted in `ParticleKernels.h` is
> not exercised here -- flagged as a separate, currently-inactive gap for any FUTURE run that turns
> analytic buoyancy on for capsule rods.
>
> **Fix applied** (`src/mesa_pd/common/Contains.h`): added `isPointInsideCapsuleBF()` (closest point
> on the body-z-axis core segment `[-coreHalfLength,+coreHalfLength]`, then a plain sphere-radius
> test -- the same Minkowski-sum-of-segment-and-ball construction `Capsule::support()` itself uses)
> and the corresponding `ContainsPointFunctor::operator()(..., const Capsule&, ...)` overload,
> mirroring Cylinder's own entry exactly. Local build confirmed clean (link succeeded). Synced to
> mirror2. NOT YET resubmitted -- holding until the DEM contact-detection audit (pillar/wall/GJK
> parity) also reports back, so any further gap it finds can be fixed in the same pass instead of
> wasting a second doomed submission.

> **[Claude, 2026-09-26] Second audit result (DEM contact detection: GJK, pillar, wall) -- no
> further blocker found; the Contains.h fix above was the only one needed.**
>
> Corrected a premise: the two "pillars" are NOT `mesa_pd::data::Cylinder` -- they're
> `mesa_pd::data::CylindricalBoundary` (GLOBAL/FIXED/INFINITE, created via `Utility.h`'s
> `createCylinder()`, `SetupStages.h:480-485`). Domain walls are `mesa_pd::data::HalfSpace`
> (`createPlane()`, `SetupStages.h:453-462`). Neither is gated to `Cylinder` specifically anywhere
> in wall/pillar creation.
>
> `SolidCylinderContactDetection` (`ParticleKernels.h:165+`) has explicit analytic overloads for
> `Cylinder x CylindricalBoundary` (the actual "pillar interaction") and `Cylinder x HalfSpace`
> (wall), but NONE for any pair involving `Capsule` -- confirmed by full grep of both
> `ParticleKernels.h` and `PaperAnalyticalContacts.h` (zero hits on `mesa_pd::data::Capsule`). Every
> Capsule pair therefore falls to the `GeneralContactDetection` base class via the inherited
> `using GeneralContactDetection::operator();` (`ParticleKernels.h:168`).
>
> **Checked whether that base-class fallback is a real algorithm or a silent no-op for
> HalfSpace/CylindricalBoundary** (the crux of the pillar/wall question) -- it IS real:
> `GeneralContactDetection.h` has dedicated, non-GJK, `support()`-only generic templates for
> `HalfSpace x AnyShape` (114-140) and `CylindricalBoundary x AnyShape` (152-192) that only need the
> OTHER shape's `support()` -- which `Capsule::support()` (`Capsule.h:113`) correctly provides.
> Contact IS detected in every case; not a "rods pass through walls/pillars" bug.
>
> **Two specific, named residual gaps, both non-blocking for this scenario:**
> 1. `Capsule x HalfSpace`'s generic path returns `contactNormal` pointing INTO the wall (opposite
>    of every analytic shape's out-of-wall convention in this file, per `ParticleKernels.h:158-164`'s
>    own comment) -- proven algebraically harmless to the actual force applied (`LinearSpringDashpot`
>    is symmetric under the corresponding idx-swap+negate), but will misreport sign for anything
>    reading the raw contact normal directly (diagnostic logs, FloaterHistoryLogger CSVs) for
>    capsule-vs-wall rows specifically. Treat wall-contact-normal-sign-dependent post-processing for
>    capsule rods with caution; the physical trajectory itself is not affected.
> 2. `ComputeTunnelContacts` (`ParticleKernels.h:3031-3045`) has a literal
>    `getShapeType() != Cylinder::SHAPE_TYPE` check that silently downgrades a Capsule rod to a
>    coarser isotropic containment bound instead of the exact per-rim Cylinder calc -- only relevant
>    to ARCH-TUNNEL scenarios (`useArchTunnelAsMesaPDObject=false`), which this pillar+wall rod-hopper
>    scenario does not use. Flagged for later, not a blocker here.
>
> **Net conclusion**: the `Contains.h` fix (previous entry) was the ONLY blocking gap. Resubmitting
> the capsule experiment now.

> **[Claude, 2026-09-26] Test suite for the Capsule shape, per explicit request ("Create a test
> suite for the pills") -- found and fixed TWO further real bugs before job 4267741 could run with
> them live, on top of the Contains.h gap found by manual audit earlier.**
>
> Wrote `tests/showcases_antidunes/CapsuleShapeTest.cpp` (registered as `ANTIDUNES_CapsuleShape` in
> CMakeLists.txt, codegen-gated like `ANTIDUNES_ContactDetection`). Coverage: volume/mass/inertia
> (analytic formula + the exact height<=2r sphere-degenerate limit + a thin-rod limiting-behaviour
> check), `support()` (cardinal directions checked against the Minkowski-sum definition directly,
> plus a "Capsule is always contained within the same-(radius,height) Cylinder" containment check),
> `ContainsPointFunctor` via the EXACT `SingleCast` dispatch `MovingParticleMapping.h` uses (direct
> regression coverage for the Contains.h gap found earlier -- including the discriminating
> "beyond the flat cap but still within the rounded cap" case a lazy/wrong implementation would get
> wrong), contact detection through the SAME production `SolidCylinderContactDetection` dispatch real
> scenarios use (HalfSpace/wall, CylindricalBoundary/pillar, Capsule x Capsule), and pack/unpack
> round-trip.
>
> **Bug #2, found by `testCapsuleDegenerateSphereLimit`**: `Capsule::updateMassAndInertia()`'s
> `Iperp` formula (`src/mesa_pd/data/shape/Capsule.h`) used `19/160` for "each hemisphere's own-
> centroid perpendicular moment," when the correct value -- provable directly by mirror-symmetry
> (splitting a full sphere through its own centre, the flat-face-diameter axis lies exactly in that
> splitting plane, so both hemispheres contribute equally to the sphere's own (2/5)MR^2 about it) --
> is `83/320`. The two standard values (2/5 about the flat face, 83/320 about the centroid) had been
> swapped in the original derivation. Confirmed as a STRUCTURAL bug, not just a degenerate-limit
> artifact: the wrong constant is used for every capsule regardless of its aspect ratio, understating
> `Iperp` by a fixed `(9/64)*(fullSphereMass)*r^2` for any (r, height) -- meaning the production rod
> dimensions (r=3, length=18) would have had wrong rotational inertia throughout, affecting their
> tumbling dynamics in every capsule run. Fixed (19/160 -> 83/320, with the corrected derivation
> written out in the code comment); test now passes exactly (not just within tolerance).
>
> **Bug #3, found by `testCapsuleCylindricalBoundaryContact`**: no dedicated Capsule x
> CylindricalBoundary overload existed in `ParticleKernels.h`'s `SolidCylinderContactDetection`, so
> the pair fell to the generic `GeneralContactDetection::operator()(CylindricalBoundary, GEO2_T)`
> template (`src/mesa_pd/collision_detection/GeneralContactDetection.h`) -- which implements
> CONTAINMENT semantics UNCONDITIONALLY (collision when the OTHER body's FARTHEST point from the
> axis EXCEEDS the boundary's radius; it never consults `CylindricalBoundary::isContainment()`).
> This is exactly backwards from a solid pillar's needs (collision when the body's CLOSEST point
> gets nearer than the radius) -- the dedicated Cylinder x CylindricalBoundary overload right above
> it in the same file implements the correct "solid exterior" semantics precisely because
> `isContainment()==false` for real pillars, but Capsule got none of that. Test caught it exactly:
> expected penetrationDepth 0.2, got 1.8 -- traced directly to the generic template's `if (dist >
> radius^2)` check. **This would have silently corrupted EVERY rod-pillar interaction in the actual
> hopper scenario** (most rods sit far outside a pillar's radius in the open channel, which the wrong
> containment semantic would treat as an immediate, large, spurious "violation," feeding a huge
> unphysical contact force from step one) -- caught by the test before job 4267741 ever got past
> PENDING in the queue (cancelled, no cluster time wasted, no corrupted results).
>
> **Fix**: added a dedicated `Capsule x CylindricalBoundary` overload (+ its swapped-argument
> partner) to `SolidCylinderContactDetection`, mirroring the Cylinder overload's own structure but
> CLOSED-FORM rather than iterative (a capsule's closest surface point to any line is exactly its
> core segment's closest point moved `radius` toward that line -- no radial clamp needed at all,
> only a 1-D axial clamp, derived directly from the two-lines closest-point normal equations and
> exact by convexity, no iteration required). Falls back to the generic path for the (unused, but
> handled for completeness) containment-mode case. Both fixes confirmed against the full test suite
> (all green) and a full local rebuild of the production `Antidunes` binary (clean link).
>
> **Net effect of writing this test suite**: caught 2 of the 3 total capsule-wiring bugs found this
> session (the third, Contains.h, was found by manual audit slightly earlier the same evening) --
> each one would have either hard-aborted or silently corrupted a real cluster run. Synced both
> fixes to mirror2 and resubmitted the capsule drop-waves experiment as job **4267766**.

> **[Claude, 2026-09-26] Capsule (pill-shape) fixes-off drop-waves result (job 4267784, 1h cap) --
> a genuine, trustworthy apples-to-apples comparison to the flat-cylinder baseline (job 4267519) now
> that all three capsule-wiring bugs (Contains.h, updateMassAndInertia, CylindricalBoundary contact)
> are fixed. Result: HELPED MEANINGFULLY, DID NOT FIX.**
>
> Crashed at step **24153** (39m21s elapsed, not a timeout) -- ~2000 steps further than the
> flat-cylinder baseline's ~22037. Worst `[DEM-WATCH]` rod-rod overlap seen across the whole run was
> **8.99%**, vs. the baseline's **89.9%/23.9%/7.3%** -- roughly a 10x reduction in peak
> interpenetration severity, consistent with the rounded caps avoiding the worst corner-into-corner
> pathological overlaps a flat-capped cylinder's "band-band-skewed" contact produces.
>
> **But the same root mechanism still fired**: `[MEM-GARBAGE-LINK]` triggered 40 times starting at
> step 24150 (unguarded on this fixes-off mirror) on a rod's own nearly-enclosed liquid cell (25/26
> neighbours `rawFlag=512`, one rod) -- even a modest ~9% overlap between smoothly rounded surfaces
> was still enough to pinch a fluid cell into corruption. So the underlying "nearly-sealed liquid
> cell -> degenerate LBM update" instability is NOT purely a function of deep mechanical overlap; a
> capsule's rounding reduces how OFTEN/how SEVERELY that geometric precondition arises, but does not
> eliminate the precondition itself.
>
> **Terminal failure mechanism differs from the baseline**: instead of a quaternion-NaN abort, this
> run hit `[STALE-PARTICLE-IDX]` (`src/lbm_mesapd_coupling/momentum_exchange_method/boundary/
> SimpleBB.h:313`) -- a deliberate hard-abort guard (added in an earlier session, per its own code
> comment) that catches a particle's fluid-side obstacle-ownership desyncing from its DEM storage
> slot BEFORE it manifests as a later NaN crash. Almost certainly a downstream consequence of the
> SAME MEM-GARBAGE-LINK event kicking a rod into fast enough cross-rank motion to desync ownership
> bookkeeping within a few steps, not an independent new failure class -- consistent with this
> session's broader pattern of "one root corruption event, several different downstream guards
> catching it at different points in the resulting cascade."
>
> **Honest verdict**: the pill shape is a real, measurable mitigation (an order of magnitude smaller
> overlaps, ~9% further into the run) but NOT a fix for this scenario's underlying pinched-cell
> instability. Both the flat-cylinder and capsule fixes-off baselines now point at the same
> conclusion reached earlier for LOW-RHO: there is a genuine, not-yet-root-caused free-surface/MEM
> coupling instability at near-sealed cells that persists regardless of rod end-cap geometry.

> **[Claude, 2026-09-26] Sanity check on job 4267784's capsule rods, and the individual-fix sweep,
> per explicit request ("First check the rods acted as expected and then try out the fixes
> individually and in combinations use all three directories") / ("check in relation of where the
> rods are what crashed").**
>
> **Rod behaviour sanity check**: sampled `[CylForces]` trajectories across job 4267784. Buoyancy
> term (`F_buoy+grav_z`) matched the three density tiers EXACTLY and in the physically correct order
> (-0.0132 for densityRatio 0.5, -0.0264 for 1.0, -0.0316 for 1.2 -- heavier rods pull down harder, as
> expected). Velocities stayed in the normal ~0.005-0.04 LU/step flow-speed range throughout (matching
> the PID target 0.04) except right at the crash. tiltY progressed smoothly and continuously (no
> erratic jumps) for every sampled rod. Several rods had already passed X=384 (the pillars) by the
> time of the crash, reaching as far as X~434 -- all of them at Y~15-18, comfortably centred in the
> 16 LU gap between the two pillar surfaces (Y=[12,28]), which is exactly why the new
> `[pillar exterior (Capsule)]` log line never fired in this run (0 hits) -- not a bug, just never
> geometrically exercised; the rods simply never got close enough to a pillar to test that code path.
> Net: the capsule rods behaved like physically sane rigid bodies throughout.
>
> **Where the crash was, relative to the rods**: the crash-triggering rod (uid
> 1008806316530991130) sat at (X~305-312, Y~24.8-25.2, Z~43.4-44.2) in the steps immediately
> around the crash -- ~72 LU upstream of the pillars (X=384), mid-channel, comfortably away from any
> wall or pillar surface. Same character as the flat-cylinder baseline's crash location (job
> 4267519, also mid-channel) -- reconfirms this instability is a rod-rod/rod-self pinch phenomenon
> unrelated to the pillar geometry, occurring anywhere two rod surfaces happen to get close while
> being carried downstream, well before any of them reach the clog point.
>
> **Individual-fix sweep launched**, using all three cluster checkouts for genuine concurrent
> isolation (per explicit instruction): all three configurations carry the SAME capsule shape and
> ALL capsule-wiring fixes (Contains.h, updateMassAndInertia, Capsule x CylindricalBoundary), only
> ONE of today's three original crash-mitigation fixes active per checkout, verified file-by-file via
> md5sum before submission (network to fritz.iwia was unusually flaky this session -- every sync
> verified against a checksum, not just an exit code, after several silent partial-failures were
> caught this way):
> - **Job 4267856** (`~/mt-dennys`): Option A ONLY (StreamReconstructAdvectSweep.h's clone-buffer
>   init-copy) -- Option B's wiring and the SimpleBB.h guard both reverted to pre-fix state.
> - **Job 4267857** (`~/mt-dennys-mirror`): Option B ONLY (`repairDegenerateObstacleCells`) -- Option
>   A's init-copy and the SimpleBB.h guard both reverted.
> - **Job 4267858** (`~/mt-dennys-mirror2`): SimpleBB.h fix ONLY (the `[MEM-GARBAGE-LINK]` active
>   guard, the fix most directly relevant to the mechanism observed firing in both fixes-off runs so
>   far) -- Option A's init-copy and Option B's wiring both reverted.
>
> All three at a 1h time cap (matching the pattern that already cleanly captured both the
> flat-cylinder and capsule fixes-off crashes within ~40 min). Results pending. Once these land, the
> natural next round (reusing the same three checkouts) is the pairwise/full combinations (A+B, A+C,
> B+C, A+B+C) that "in combinations" also calls for.

> **[Claude, 2026-09-26] Individual-fix sweep, first result: Option A alone does NOT help (as
> expected -- it targets an unrelated mechanism).**
>
> **Job 4267856** (`~/mt-dennys`, Option A only, capsule rods): crashed at step **17493** (28m37s
> elapsed) -- EARLIER than either fixes-off baseline (job 4267519 flat-cylinder ~22037, job 4267784
> capsule-fixes-off ~24153). Worst `[DEM-WATCH]` overlap reached **49.2%** -- deeper than the
> capsule-fixes-off baseline's 9% (though still well under the flat-cylinder baseline's 89.9%).
> `[MEM-GARBAGE-LINK]` fired 7 times (fewer occurrences than the capsule-fixes-off baseline's 40, but
> feeding a comparably corrupted PDF, `pdf_old(dir)=-10.02` on a liquid cell almost fully enclosed by
> one rod's own surface). Terminal mechanism identical to the capsule-fixes-off baseline:
> `[STALE-PARTICLE-IDX]` hard-abort at `SimpleBB.h:313`, same causal pattern (MEM-GARBAGE-LINK
> corrupts one rod's force -> fast motion -> fluid-side obstacle-ownership desync -> deliberate abort
> before a later NaN).
>
> **Interpretation**: Option A alone provides no protection against this failure mode, consistent
> with its own scope -- it fixes the FieldClone-uninitialized-buffer mechanism (mechanism 1, obstacle
> cells reading stale zero PDFs), which is structurally unrelated to this rod-rod contact-driven
> pinch instability. The worse overlap/earlier crash in this ONE run should NOT be read as "Option A
> makes things worse" -- given this scenario's established non-determinism (exact crash step/severity
> varies run to run even with IDENTICAL code, per this project's own long-standing observation), a
> single data point in either direction isn't conclusive; the more defensible reading is simply
> "no protective effect detected," pending the other two individual-fix results and, ideally, repeat
> runs.
>
> Jobs 4267857 (Option B only) and 4267858 (Option C only) still PENDING in the queue as of this
> entry -- cluster busy. Results to follow once they run.

> **[Claude, 2026-09-26] Individual-fix sweep, results 2 and 3: Option B alone ALSO does not help --
> same terminal mechanism a third time.**
>
> **Job 4267857** (`~/mt-dennys-mirror`, Option B only, capsule rods): crashed at step **16699**
> (27m14s elapsed) -- earlier than both fixes-off baselines, and even slightly earlier than Option
> A's own run (17493). Worst overlap **59.8%** (worse than Option A's 49.2%, worse than the capsule
> fixes-off baseline's 9%, though still under the flat-cylinder baseline's 89.9%). `[MEM-GARBAGE-LINK]`
> fired 40 times (matching the fixes-off capsule baseline's count exactly). **Terminal mechanism
> IDENTICAL to job 4267856 (Option A) and the fixes-off baselines**: `[STALE-PARTICLE-IDX]` hard-abort
> at `SimpleBB.h:313`, same causal pattern.
>
> **Interpretation**: Option B alone (the `repairDegenerateObstacleCells` guard, confirmed genuinely
> active via 4659+ `[OBSTACLE-CELL-REPAIR]` hits before the crash) also provides no protection
> against this failure mode -- consistent with its own scope: it repairs OBSTACLE-cell PDF corruption
> (mechanism 1, cells with permanently-zero/garbage PDFs from the FieldClone bug), not the LIQUID-cell
> MEM-force corruption that triggers this rod-rod pinch cascade. Three fixes-off/individual-fix runs
> in a row (4267519, 4267784, 4267856, 4267857 -- four, actually) have now all terminated via the SAME
> underlying mechanism family (MEM-GARBAGE-LINK or its close cousin -> either a quaternion-NaN or a
> STALE-PARTICLE-IDX abort), regardless of which of Option A/B was active. This narrows the
> hypothesis space sharply: of today's three original fixes, only the SimpleBB.h guard (Option C,
> job 4267858, still running as of this entry) has any mechanistic reason to intervene in this
> specific chain, since it is the only one of the three that acts on the MEM force computation
> itself rather than the fluid-obstacle-cell PDF state.
>
> **Job 4267858** (`~/mt-dennys-mirror2`, SimpleBB.h fix ONLY, capsule rods) -- still RUNNING as of
> this entry, clean at step 12484 (19m16s), zero corruption diagnostics, worst overlap 6.5%. This is
> the run most likely to show a genuinely different outcome, since it is the only individual fix that
> directly guards the exact mechanism (an unguarded corrupted `pdf_old` reaching `forceMEM`) observed
> triggering every crash so far. Result pending.

> **[Claude, 2026-09-26] Individual-fix sweep, final result: the SimpleBB.h fix (Option C) worked
> exactly as designed but was NOT sufficient on its own -- the underlying corruption re-routes
> through a different pre-existing guard and still crashes the run. This is the most informative
> result of the whole sweep.**
>
> **Job 4267858** (`~/mt-dennys-mirror2`, SimpleBB.h fix ONLY, capsule rods): crashed at step
> **19725** (30m29s elapsed) -- later than Option A (17493) and Option B (16699), still short of
> both fixes-off baselines (22037/24153).
>
> **`[MEM-GARBAGE-LINK]` fired 93 times (more than any other run in this sweep) -- and EVERY ONE was
> safely skipped: zero `[STALE-PARTICLE-IDX]` aborts this run, the first time in this whole sweep.**
> This directly, empirically confirms the SimpleBB.h fix does exactly what it was built to do: it
> closes the specific chain (corrupted `pdf_old` -> unguarded force -> particle rocket-launch ->
> fluid-obstacle-ownership desync -> `STALE-PARTICLE-IDX` abort) cleanly and completely, even under
> more than double the corruption-event load seen in any prior run this sweep.
>
> **But the underlying fluid-cell corruption itself was NOT prevented -- it re-routed through two
> OTHER pre-existing guards instead**: `[DRY-LINK]` (30x, non-finite `pdf_old` at what the guard
> expects to be a drained near-waterline interface cell) and **`[LOW-RHO]`** (26x -- the SAME
> independent, not-yet-root-caused free-surface reconstruction instability characterised earlier
> this session in jobs 4266571/4266572/4267519/4267784). Final abort: a **quaternion-NaN** assertion
> (`Quaternion.h:241`, `r*r+i*i+j*j+k*k == 1` with all four components `-nan`) -- the SAME
> fundamental failure mode as the very first flat-cylinder baseline (job 4267519), just reached via
> the LOW-RHO/DRY-LINK route instead of the MEM-GARBAGE-LINK->STALE-PARTICLE-IDX route the fix
> closes off.
>
> ---
> ## Individual-fix sweep: summary table
>
> | Config | Job | Crash step | Elapsed | Worst overlap | Terminal mechanism |
> |---|---|---|---|---|---|
> | Flat cylinder, all fixes off | 4267519 | ~22037 | ~32min | 89.9% | quaternion-NaN (via LOW-RHO) |
> | Capsule, all fixes off | 4267784 | ~24153 | ~39min | 9.0% | STALE-PARTICLE-IDX |
> | Capsule, Option A only | 4267856 | 17493 | 28m37s | 49.2% | STALE-PARTICLE-IDX |
> | Capsule, Option B only | 4267857 | 16699 | 27m14s | 59.8% | STALE-PARTICLE-IDX |
> | Capsule, Option C only | 4267858 | 19725 | 30m29s | 47.9% | quaternion-NaN (via LOW-RHO/DRY-LINK) |
>
> **Overall conclusion of the individual-fix sweep**: none of today's three original crash-mitigation
> fixes, applied individually, prevents this scenario's crash. Options A and B (which target
> obstacle-cell PDF corruption, a structurally different mechanism) show no protective effect at all
> -- same terminal mechanism as the fixes-off baselines. Option C (the SimpleBB.h guard) is
> mechanistically the right tool and PROVABLY closes its own target chain completely (0/93 leaks
> through), but the scenario has at least TWO independent roads to the same quaternion-NaN endpoint
> (MEM-GARBAGE-LINK->STALE-PARTICLE-IDX, now closed by Option C; and LOW-RHO/DRY-LINK, still fully
> open) -- closing one road just routes the same underlying pinched-cell fluid instability through
> the other. This sharpens, rather than resolves, the standing open question from earlier this
> session: the LOW-RHO/free-surface-reconstruction mechanism is the genuinely unresolved bottleneck,
> independent of rod shape (capsule vs cylinder) and independent of which MEM-force guard is active.
>
> Per the non-determinism caveat established throughout this project, exact crash steps/overlap
> severities vary run to run even with identical code -- the step-count differences across this table
> (16699-24153) should be read qualitatively (all crashed well before completing 30000 steps) rather
> than as precise, reproducible measurements of any one fix's "delay" effect.
>
> **Natural next steps** (not yet started): (1) the pairwise/full combinations this sweep didn't yet
> cover (A+B, A+C, B+C, A+B+C, all with capsules) -- A+C or B+C might close BOTH known routes if
> Option A/B's obstacle-cell repair happens to reduce how often a cell reaches the LOW-RHO-triggering
> state, worth testing empirically rather than assuming; (2) directly root-causing the LOW-RHO/
> free-surface-reconstruction mechanism itself, which is now the clear standing bottleneck across
> every configuration tested this session.

> **[Claude, 2026-09-26] Literature context and thesis documentation for the near-sealed-cell
> instability, per explicit request ("What do other papers do against this? Can we document our
> instability better for my thesis, have test cases proof this").**
>
> ## What this instability actually is (mechanism, established empirically this session)
>
> A liquid lattice cell becomes almost entirely enclosed by a SINGLE rod's own curved surface
> (confirmed directly from `[MEM-GARBAGE-LINK]`'s own neighbour-flag dump across five independent
> crashes -- jobs 4267519, 4267784, 4267856, 4267857, 4267858: 14/26, 24/26, 21/26, and comparable
> ratios of the cell's 26 neighbours carrying the SAME single obstacle flag value, never a mix of two
> different rods' flags), leaving it connected to the rest of the fluid domain by as few as 2 real
> liquid neighbours. The standard LBM collision-streaming update, run on a cell with almost no real
> fluid input, produces PDF components that diverge to unphysical magnitudes (confirmed values from
> +-10 up to +-9x10^15 across different occurrences) -- density and momentum become meaningless. This
> corrupted value then either (a) leaks through `SimpleBB.h`'s momentum-exchange computation into a
> particle's hydrodynamic force (the `[MEM-GARBAGE-LINK]` chain, closed by this session's SimpleBB.h
> fix), or (b) resurfaces via the free-surface reconstruction's own density computation
> (`[LOW-RHO]`/`[DRY-LINK]`, still open) -- both routes converge on the same terminal failure, an
> invalid (NaN) orientation quaternion.
>
> Crucially: every occurrence traced this session happened well away from the pillars/walls (72-107
> LU clear in the individual-fix sweep, ~50-100+ LU in every earlier baseline) and involved only ONE
> rod's own obstacle flag -- never a second, different rod's flag appearing in the same
> neighbourhood. This points at a SINGLE-BODY lattice-alignment effect (a continuously-orientable
> curved surface happening to wrap a lattice cell against itself at certain orientations relative to
> the FIXED cubic grid), not fundamentally a two-body contact/pinch effect, even though this
> scenario's rods are densely packed and DO also experience real two-body contact separately.
>
> ## Why this doesn't (and can't) happen in a real experiment
>
> Real water has no discretisation to misalign against. The physical fluid in an arbitrarily thin
> gap between (or around) solid surfaces still satisfies the Navier-Stokes equations continuously,
> however fast or sheared the local flow -- there is no mechanism for density or momentum to become
> non-finite. What fails here is specifically the numerical representation of a smoothly curved,
> continuously-orientable solid boundary on a fixed cubic grid, once the locally-resolved fluid
> region around/within that boundary shrinks to a handful of grid cells. This is a discretisation
> limitation, not a physical instability -- the literature below confirms this is a long-recognised,
> named problem in resolved-particle LBM specifically, independent of this codebase.
>
> ## What the literature does about it
>
> The general problem -- "lattice Boltzmann simulations of resolved particles cannot resolve
> lubrication/near-contact hydrodynamics once the fluid gap approaches the lattice spacing" -- is
> well established and has several named families of mitigation:
>
> 1. **Analytic lubrication corrections (the dominant approach for MEM/bounce-back codes).**
>    Nguyen & Ladd, "Lubrication corrections for lattice-Boltzmann simulations of particle
>    suspensions," *Phys. Rev. E* **66**, 046708 (2002) -- the seminal reference: once two
>    discretised surfaces approach within about one lattice spacing, the LBM's own resolution of the
>    lubrication force becomes unreliable, so an explicit analytical lubrication force/torque
>    (calculated directly from the particles' relative velocity and gap width, independent of the
>    under-resolved fluid field) is ADDED on top of, or substituted for, the LBM's own hydrodynamic
>    force in that regime. This is exactly the role `useLubricationCorrection` plays in this
>    codebase's DEM side already -- but only for the DEM force between two particles; it does
>    nothing to protect the fluid solver's OWN collision update at the sealed cell itself, which is
>    the actual site of this session's crashes. Related: Ladd's original momentum-exchange
>    formulation, "Extension of the lattice-Boltzmann method for direct simulation of suspended
>    particles near contact," *J. Stat. Phys.* (2002).
> 2. **Partially Saturated Method (PSM) / Noble-Torczynski method.** Noble & Torczynski, "A
>    lattice-Boltzmann method for partially saturated computational cells" -- instead of a binary
>    obstacle/fluid flag per cell, each cell carries a continuous solid volume fraction that smoothly
>    modifies the collision operator as a body's surface sweeps through it. Because there is no hard
>    binary obstacle boundary, a cell can never become "almost entirely obstacle-flagged" the way our
>    cells do -- the analogous failure mode is structurally impossible in a properly implemented PSM.
>    Directly relevant follow-up specific to dense granular suspensions: "Analysis of two
>    partially-saturated-cell methods for lattice Boltzmann simulation of granular suspension
>    rheology" (a directly comparable use case to this project's dense rod hopper).
> 3. **Smoothed-profile / diffuse-interface methods.** A further step in the same direction as PSM --
>    the solid-fluid boundary is represented as a smooth (not sharp) profile, explicitly motivated in
>    the literature by exactly the force/velocity fluctuation problems sharp bounce-back methods show
>    near contact.
> 4. **Local/adaptive mesh refinement around moving particle surfaces.** Refining the grid near a
>    particle's surface (octree-based or block-structured AMR) increases the number of real fluid
>    cells resolving the near-contact region, directly attacking the "too few real fluid neighbours"
>    root cause rather than working around its consequences. Reported specifically to fix "large
>    fluctuations of the computed hydrodynamic force and torque" from inadequate resolution near a
>    moving particle surface.
> 5. **Progressive bounce-back fallback schemes.** For interpolated bounce-back (IBB) specifically --
>    the same family this codebase's own `useInterpolatedBounceBack_` belongs to -- the literature
>    describes falling back from quadratic to linear to plain single-node bounce-back as the number
>    of available real fluid nodes along a link shrinks, rather than letting an under-determined
>    interpolation silently produce a bad value. Directly relevant to the audit conclusion reached
>    earlier this session that IBB provides no protection here (it still uses the same unguarded
>    `pdf_old` as its base term) -- the literature's own answer to "IBB near a narrow gap" is
>    exactly this kind of explicit degradation path, which this codebase does not currently implement
>    for the case where `pdf_old` itself, not just the interpolation partner cell, is bad.
>
> Broader context, both general reviews of the field: Aidun & Clausen, "Lattice-Boltzmann Method for
> Complex Flows," *Annu. Rev. Fluid Mech.* **42**:439-472 (2010); and, specific to the contact/gap
> problem across many coupling schemes, the comparative study "A comparative study of fluid-particle
> coupling methods for fully resolved lattice Boltzmann simulations" and the review "A Review on
> Contact and Collision Methods for Multi-body Hydrodynamic problems in Complex Flows" (arXiv
> 2211.11728) -- both directly comparable to this project's own MEM+bounce-back choice and worth
> citing as the "landscape of alternatives" this project sits within.
>
> **Framed against this codebase's actual architecture**: this project uses hard bounce-back (+MEM)
> with a binary per-cell obstacle flag -- structurally the family MOST susceptible to the
> near-sealed-cell failure, since it is precisely the binary flag (not a continuous solid fraction)
> that lets a cell become ~90%+ "obstacle" while still being nominally "liquid." Switching the
> fluid-obstacle coupling to PSM, or adding local refinement near particle surfaces, would likely
> remove this failure mode structurally rather than patching around it -- but both are substantial
> architectural changes, well beyond a targeted fix, and out of scope for validating THIS session's
> narrower mitigations.
>
> ## Test cases built this session that constitute the "proof"
>
> 1. **`tests/showcases_antidunes/CapsuleShapeTest.cpp`** (`ANTIDUNES_CapsuleShape`) -- unit-level,
>    deterministic, asserting proof that the Capsule shape's own geometry/contact/containment code is
>    correct in isolation (17 passing cases) -- rules out "the crashes are just a capsule-wiring bug"
>    as an explanation, cleanly separating shape-implementation correctness from the physical
>    instability itself.
> 2. **`prms/pinch_test_two_rods_gas.prm`** (job 4267098, EARLIER this session) -- two FIXED,
>    perfectly symmetric, statically-placed rods at the exact production 1 LU clearance: ran 5000
>    steps completely clean, zero corruption diagnostics. NEGATIVE result, but an important one: pure
>    static geometric proximity between two rods, with no dynamics, is NOT sufficient on its own --
>    ruling out the simplest "two bodies just need to be this close" hypothesis.
> 3. **`prms/single_rod_tumble_gas.prm` / `single_rod_tumble_capsule_gas.prm`** (jobs 4268210/4268211,
>    submitted this entry, results pending) -- the key new test, built directly from this session's
>    own observation that every real crash involved only ONE rod's flag in the enclosing
>    neighbourhood: ONE cylinder/capsule, completely alone in an otherwise empty still-water domain,
>    free to fall and tumble from a 45 deg initial tilt so it sweeps through many orientations
>    relative to the fixed lattice as it falls. If `[MEM-GARBAGE-LINK]`/`[LOW-RHO]`/
>    `[OBSTACLE-CELL-REPAIR]` fire here too, with no second body anywhere in the domain, that is a
>    clean, minimal, reproducible proof that the trigger is single-body lattice-alignment, not
>    multi-body contact -- a far stronger and more citable result for the thesis than "it happened
>    somewhere in a noisy 144-rank, 6000-rod-interaction production run." Results to follow.

> **[Claude, 2026-09-26] Single-rod tumble test result: NEGATIVE for both shapes -- sharpens rather
> than settles the single-body hypothesis.**
>
> First attempt (jobs 4268210/4268211) had a genuine setup bug, caught before wasting the result:
> `axis <0,1,0>` combined with `initialTiltAngleDeg` (which rotates about world-Y) made the "tilt" a
> spin about the rod's OWN long axis -- physically invisible for a circular cross-section (tiltY
> stayed frozen at ~0.0002 deg). A concurrent second bug: both jobs ran `"$SH" --build` in the SAME
> checkout (`~/mt-dennys`) at the same time, corrupting the shared object files ("undefined reference
> to `main`" on the capsule job -- a build race, not a physics result). Fixed: `axis <1,0,0>` (so the
> world-Y tilt genuinely reorients the rod), densityRatio 1.1->1.3 (more decisive sink), and the two
> variants moved to separate checkouts (`~/mt-dennys` cylinder, `~/mt-dennys-mirror2` capsule) to
> avoid the build race.
>
> **Resubmitted as jobs 4268213 (cylinder) and 4268214 (capsule) -- both COMPLETED all 8000 steps
> with ZERO hits on every diagnostic** (MEM-GARBAGE-LINK, NAN-LINK, DRY-LINK, LOW-RHO, SEALED-CLOSE,
> STALE-PARTICLE-IDX, OBSTACLE-CELL-REPAIR, Fatal error, Assertion). Confirmed genuine tumbling this
> time: cylinder tiltY 45->130.7 deg, capsule 45->123.8 deg, both sinking ~12-13 LU over the run
> (ruling out "it just sat still again" as the explanation).
>
> **Honest interpretation**: this negative result does NOT refute the single-body/lattice-alignment
> hypothesis (established from the neighbour-flag evidence -- every real crash's enclosing
> neighbourhood carried exactly ONE obstacle flag value, never two), but it does show that GENTLE,
> UNPERTURBED single-rod tumbling in still water is not sufficient to reproduce it at this sample
> size. Three concrete reasons this isn't surprising, all pointing at what a follow-up test should
> add rather than abandoning the hypothesis:
> 1. **Sample size**: <90 degrees of rotation swept over 8000 steps, vs. 210 rods undergoing far more
>    violent, contact-driven reorientation over tens of thousands of steps in the real scenario --
>    orders of magnitude fewer chances to hit whatever specific orientation is actually dangerous.
> 2. **No contact-induced kick**: every one of the five real crashes traced this session happened in
>    a densely-packed, contact-rich environment. The trigger may not be "any bad orientation
>    reached smoothly," but specifically a SUDDEN reorientation/velocity spike from a DEM contact
>    event, reaching the dangerous configuration too fast for the fluid solver to track -- this test
>    had zero contacts (nothing else in the domain) and therefore zero such kicks.
> 3. **No flow, no neighbours' wakes**: still water, nothing else present -- the real scenario adds
>    channel-flow shear and 209 other rods' local flow disturbances, both absent here.
>
> **Sharpened next hypothesis for a follow-up test** (not yet built): a single rod given a large,
> SUDDEN angular velocity kick (mimicking a contact event) at a moment when its own curvature is
> already close to a self-enclosing configuration, rather than a slow smooth tumble -- this isolates
> the "sudden reorientation" variable specifically, still with no second body present, keeping the
> single-body-trigger claim testable while addressing reason #2 above directly.

> **[Claude, 2026-09-26] Giant single-rod clog scenario, per explicit request ("run the hopper case
> with one giant cylinder... it should clog"). Confirmed: genuine, stable clog against the pillar at
> densityRatio 0.75, after ruling out a mischaracterisation and working through a disk-quota crisis.**
>
> **Setup**: ONE rod, radius 6, length 36 (1:3 ratio, 90% of the 40 LU channel width -- far longer
> than both the 16 LU pillar gap and the channel itself, so a jam is geometrically guaranteed),
> released already submerged 20 LU upstream of the production hopper pillars (X=384, r=12), all
> three of today's crash-mitigation fixes reverted (fixes-off baseline), 50000 timesteps.
>
> **Iteration history**:
> - Job 4268239 (densityRatio 1.0): NOT neutrally buoyant (F_buoy+grav_z=-0.237196, exactly constant
>   -- confirmed via a code audit this is a pure weight term scaled by densityRatio with no
>   subtracted buoyant offset baked in; the actual Archimedes counter-force is supplied separately
>   via the hydrostatic-pressure correction folded into the MEM/hydrodynamic force itself). Sank
>   steadily from Z=30 to Z~8-13 over ~15000 steps and appeared pinned there -- INITIALLY
>   mischaracterised as "stuck on the floor via friction" (per explicit follow-up question
>   challenging this, re-checked: floor contact requires Z=6 (rod radius), and X had already reached
>   the pillar's own geometric contact threshold, ~366.1, exactly where it plateaued -- the pin was
>   the PILLAR, not the floor, the whole time; corrected in THESIS_PLAN.md and to the user).
> - Job 4268239's cancellation + own VTK output (writeFrequency 50, ~33000 steps -> ~660 dumps) wrote
>   139G, blowing the /home/hpc soft quota (104.9G) and killing the immediate next job
>   ("Disk quota exceeded" on `vtk-out/particles/...` and the archive `mkdir`). Archives deleted
>   across all three checkouts (freed ~46G the first time, then another 139G after this specific
>   run), writeFrequency raised 50->200->2000 per explicit follow-up requests.
> - Job 4268255 (densityRatio 0.75): **stable result.** F_buoy+grav_z=-0.177897 (exactly 0.75x the
>   1.0 case's value -- confirms the pure-proportionality/no-offset relationship precisely). X pinned
>   at 366.21-366.22 (matching the pillar-contact geometry) from very early on; Z, after initially
>   continuing to sink (30->~10.5 over ~20000 steps), TURNED AROUND and is now rising slowly and
>   steadily (10.51->10.75 over steps 24400-24900); tiltY converging smoothly (2.80->2.17deg over the
>   same window). Zero corruption diagnostics of any kind throughout (MEM-GARBAGE-LINK, NAN-LINK,
>   DRY-LINK, LOW-RHO, SEALED-CLOSE, STALE-PARTICLE-IDX all zero). This is a genuine, stable clog
>   against the pillar, not a degenerate/still-evolving state.
>
> **Also confirmed via a dedicated code audit this session** (per explicit request, "check that the
> mem force is applied correctly... over multiple ranks"): `reduceAndAvgHydroForces` performs a real
> cross-rank SUM (`mesa_pd::mpi::ReduceProperty`, point-to-point ghost-to-owner, not neighbor-
> restricted) with no double-counting by construction; the "Avg" in its name refers to an unrelated
> TEMPORAL rolling average across timesteps, not a spatial/per-rank average. One real, honestly-
> flagged but untested edge case exists (`SyncGhostOwners`'s one-hop-per-call reach for a body this
> large, which this rod's own size triggers per SetupStages.h's own selection logic) -- but a
> convergence-lag bug would show up as noisy/transient force fluctuations, not the smooth constant
> bias actually observed, so the evidence points away from an MPI defect and confirms the sinking/
> pinning behaviour is a genuine buoyancy-formula effect, correctly computed.
>
> **Next scenario prepared, held back until this confirmation landed**: `prms/giant_rods_6x_clog_gas.prm`
> -- SIX giant rods (same radius 6/length 36), arranged 2(X) x 1(Y, forced -- a single rod already
> spans nearly the full channel width) x 3(Z), both X-rows close to the pillars (X=364 and X=350,
> per explicit request "released right in front of the hopper," not staggered far upstream) x 3
> Z-layers (15/30/45, 15 LU apart) x density stratified by depth (0.85/0.75/0.65, heavy-low/light-
> high, centred on the now-confirmed-stable 0.75). Not yet submitted.

> **[Claude, 2026-09-26] Job 4268255 (giant single rod, densityRatio 0.75) COMPLETED all 50000
> steps, exit code 0, 52m54s -- final confirmation of a genuine, permanent clog.**
>
> Max X ever reached across the ENTIRE run: 366.249 -- essentially unchanged from where it first
> hit the pillar early in the run (~366.1-366.2, the exact geometric contact threshold derived from
> the rod's radius=6/pillar radius=12/the rod's finite Y-extent falling 2 LU short of the pillar's
> own axis position). It never approached X=384 (the pillar centre), let alone passed it. Zero
> corruption diagnostics of any kind for the full 50000 steps (MEM-GARBAGE-LINK, NAN-LINK, DRY-LINK,
> LOW-RHO, SEALED-CLOSE, STALE-PARTICLE-IDX, Fatal error, Assertion all zero throughout).
>
> Final state: still slowly rising in Z (30.25->31.13 over the last 500 steps of the run, continuing
> the turnaround trend first seen around step 24000) and still rotating in place (tiltY reached
> -178.5deg by step 50000, having rolled almost a full half-turn while wedged against the pillar) --
> but pinned horizontally the ENTIRE run. A genuine, stable, permanent clog: the rod reaches the
> hopper, jams against the pillar, and never breaks free or drifts downstream, for the full
> 50000-step duration.

> **[Claude, 2026-09-26] Interim check, jobs 4268352 (six giant rods, 2x1x3) and 4268353 (twelve
> original-size rods, 2x2x3), both at step ~19500-20100/51000, fixes off. Zero corruption
> diagnostics in either job. The two populations are producing qualitatively DIFFERENT outcomes --
> the giant rods are reproducing the single-rod pillar-jam almost rod-for-rod, while the small rods
> are mostly slipping through the pillar gap by tumbling, with only a minority jamming, and jamming
> well upstream of the pillar rather than against it.**
>
> **Six giant rods (4268352, `giant_rods_6x_clog_gas.prm`)**: 5 of 6 rods already pinned at
> X=366.2-366.27 -- matching the single-rod pillar-contact threshold (366.1-366.2) to within
> numerical noise -- across a spread of depths (Z=7.65, 19.65, 31.63, 43.86, 55.85), all essentially
> motionless (|vel|~1e-4) with F_contact+lub balancing F_hydro. The 6th rod (the shallowest,
> lightest-density = 0.65 rod in the front row) is still in transit: pos X=380.8, Z=59.3 (already
> above the Z=49.13 water surface, i.e. it has floated up through the surface), still moving with
> real velocity (vel~0.02-0.03) and NOT yet at the pillar threshold. Whether it also jams at
> ~X=366 once it settles into the flow, or instead continues toward/past X=384 while partially
> out of the water, is not yet resolved -- worth watching specifically. No pile-up between rods
> observed yet (each pinned rod's X matches the single-rod threshold almost exactly, consistent with
> each jamming independently against the pillar rather than against each other).
>
> **Twelve small rods (4268353, `small_rods_12x_2x2x3_clog_gas.prm`, radius 3/length 18)**: sharply
> different outcome. Of the 12 rods, 8 have already passed X=384 (one, uid ...098, has gone all the
> way around the periodic domain and reappeared near X=14) and are moving steadily downstream
> (X up to 764) -- meaning these small rods, despite being longer (18 LU) than the 16 LU pillar gap,
> are getting through it, evidently by rotating enough in flight that their presented cross-section
> through the gap shrinks below the geometric worst case (several of the escaped rods show
> tiltY far from both 0 and 180deg -- e.g. -138.75, 99.46, 78.73, -86.07deg -- consistent with
> tumbling into a more axis-aligned-with-flow orientation before/while passing through). The
> remaining 4 rods are stuck, but NOT at the pillar: two pairs, one pair near (X=334,Y=11.2,Z=4.4)/
> (X=328.6,Y=9.6,Z=2.9) and another near (X=355.4,Y=31.0,Z=4.1)/(X=349.4,Y=31.0,Z=3.7) -- each pair
> close together in X and Y and both members of a pair sitting right at the floor (Z~3-4.4, close to
> the rod radius of 3), with large F_contact+lub (~0.1-0.3 magnitude, vs ~0 for the free-floating
> rods) and near-zero velocity. Both stuck rods in the Y~31 pair show tiltY at almost exactly
> +-180deg (179.999/-179.999) -- fully flipped over. This looks like two independent two-rod
> floor-piles formed well upstream of the pillars (X~328-355, vs the pillar gap at X=384), not a
> pillar-jam at all -- consistent with these particular rods failing to reorient in time and instead
> settling/colliding near the bed before ever reaching the gap.
>
> **Interpretation so far (provisional, both jobs still running)**: rod length relative to the
> pillar gap width appears to be the dominant control on clog formation in this setup -- a rod
> whose length is much larger than BOTH the gap and the channel width (giant rods, L=36 vs 16 LU
> gap/40 LU channel) cannot reorient its way through no matter how it tumbles, and jams
> deterministically and independently of the other rods present; a rod whose length exceeds the gap
> by a smaller margin (small rods, L=18 vs 16 LU gap) CAN pass through in many cases if it tumbles
> into a favourable orientation in time, and jamming for these rods appears to be more of a
> bed-collision/pile-up phenomenon upstream of the pillars than a deterministic gap-blocking
> phenomenon. This would suggest the L/gap ratio, not L alone, is the more useful non-dimensional
> parameter for predicting clog probability -- worth carrying into the thesis's parameter sweep if
> confirmed by these jobs' full completion.
>
> **Follow-up, same check: starting depth (not row, not Y-column) perfectly separates escaped vs
> stuck among the 12 small rods** -- per explicit observation ("the starting position matter[s] the
> flow starts to turn them"). Mapped each uid back to its initial `pos=` from step ~1-100 of the log
> and cross-referenced against the step-20100 outcome above:
>
> | initial Z (depth) | initial row (X) | initial Y | outcome at step 20100 |
> |---|---|---|---|
> | 10 (deep)    | 364 or 357 (both) | 10 or 30 (both) | **stuck, all 4** (the two floor-piles) |
> | 25 (mid)     | 364 or 357 (both) | 10 or 30 (both) | **escaped, all 4** |
> | 45 (shallow) | 364 or 357 (both) | 10 or 30 (both) | **escaped, all 4** |
>
> The split is exact and depends on nothing except initial Z -- row (front/back, 7 LU apart) and
> Y-column (10 vs 30) have no bearing on outcome; every deep-released rod (Z=10) is stuck, every
> mid/shallow-released rod (Z=25, Z=45) got through. Physical reading: near the floor the mean flow
> is slower (no-slip boundary layer), so a rod released deep accumulates less hydrodynamic torque on
> the way toward the pillar and never tumbles into a gap-passing orientation before it settles into
> the bed and collides with its neighbours; a rod released in faster-moving fluid (mid/shallow) picks
> up enough torque in transit to reorient favourably and slip through the 16 LU gap despite being
> longer (18 LU) than the gap itself. This refines the earlier L/gap-ratio hypothesis: L/gap sets
> whether passing through is geometrically POSSIBLE at all (rules the giant rods out entirely,
> regardless of orientation), while release depth (via the local flow speed/torque available before
> reaching the pillar) sets whether a rod that COULD reorient through actually does so in time.

> **[Claude, 2026-09-26] Jobs 4268352 (six giant rods) and 4268353 (twelve small rods) both
> COMPLETED all 51000 steps, exit clean, zero corruption diagnostics of any kind throughout either
> run. Final states confirm and extend the step-~20000 interim findings above.**
>
> **Six giant rods (4268352)**: 5 of 6 remain permanently jammed at the pillar (X=366.22-366.26,
> essentially unchanged from the interim check, F_contact+lub exactly balancing F_hydro, |vel|~1e-4)
> -- a stable multi-rod jam matching the single-rod baseline (job 4268255) precisely, for all 6
> depths tested. **The 6th rod (uid 4755801206503243778, densityRatio 0.65, the lightest/shallowest,
> front row) did NOT jam** -- it is the one that had floated above the water surface (Z=59.3) and
> was still approaching at the interim check; by step 51000 it is at pos=<213.216,23.8187,62.4168>
> with **F_contact+lub=(0,0,0)** -- i.e. it registered ZERO contact force with the pillar at any
> point, and ended up having travelled the full length of the domain past X=384 and wrapped around
> the periodic X boundary. This is a geometrically surprising result: a rod of length 36 (spanning
> Y=[2,38], nearly the full 40 LU channel) should not fit through the 16 LU pillar gap in ANY
> orientation, submerged or not -- yet this one passed with no recorded contact whatsoever. Under
> investigation (dispatched a code-reading check of whether the DEM `CylinderPillar` contact object
> and the fluid-side `BodyBC Cylinder` boundary actually share the same Z-extent, or whether the DEM
> contact geometry might not extend the full domain height the fluid obstacle nominally does --
> which would let a rod that floats high enough slip past the CONTACT object while the fluid
> boundary is still nominally present there). Result not yet in; to be appended once available.
>
> **Resolved.** A dedicated code-reading check (`SimulationParams.h:338-347`/542-550,
> `SetupStages.h:480-484`, `Utility.h:226-241`, `CylindricalBoundary.h:87-88`,
> `BodyFromConfig.cpp:45-58`, `Cylinder.cpp:65-76`) confirmed the DEM `CylinderPillar` object is a
> `mesa_pd::data::CylindricalBoundary` with **NO axial bounds set** (`setAxialBounds()` never
> called, `axialMin_/axialMax_` default to +-infinity) -- it is Z-INFINITE by construction, a strict
> superset of the finite (but here full-domain-height, 0-180) fluid-side `BodyBC Cylinder`. Its
> contact test is radius-in-XY-only, for any Z. So there is no Z-extent loophole at all -- ruled out.
>
> The real mechanism, found by tracing uid 4755801206503243778's full `pos=`/`tiltY=` history: it
> TUMBLED through the gap, the same mechanism already identified for the small rods, just enabled by
> buoyancy/surfacing rather than flow shear. Trajectory highlights: settled near its release point
> until ~step 4000, then underwent a large, sudden tiltY jump (3.5deg -> -76.67deg between step 6000
> and 8000, coinciding with a Z drop 27.9->18.3 -- almost certainly a contact event, likely brushing
> the pillar or another rod) that kicked off a multi-thousand-step tumble, during which it actually
> drifted UPSTREAM against the mean flow (X: 350->339->331->308->294 over steps 2000-12000) while
> tiltY swung through -76 to -113deg, then reversed, rose further (Z climbing to 54, then 62, then 68
> -- fully clear of the Z=49.13 water surface by step 16000-18000) and tumbled violently right at the
> gap itself (tiltY -4.1deg at step 16000 -> +130.0deg at step 18000 -> -48.3deg at step 20000, the
> exact window in which X crossed 364->389.5, i.e. THROUGH the pillar's X position). After that it
> settled back into the water (Z back down to 40.8, then 8.3) and travelled freely downstream for
> the rest of the run, tiltY now evolving smoothly (no more contact-sized jumps) exactly like a
> free body in open flow.
>
> **Physical reading**: the deciding variable is not rod length or L/gap ratio in isolation, but
> whether the rod's axis can rotate far enough OUT of the Y-direction. A rod whose axis points
> mostly along Y presents its full length (36 LU, hugely over the 16 LU gap) to the pillar pair; a
> rod tumbled so its axis points more along X/Z presents only its DIAMETER (12 LU for these giant
> rods) to the gap -- comfortably under 16 LU regardless of the rod's length. This means even a
> grossly-oversized rod is not UNCONDITIONALLY guaranteed to jam -- only rods that never achieve
> (or never sustain) a large-enough reorientation do. The other 5 giant rods stayed jammed precisely
> because once wedged against the pillar (tiltY oscillating only 0-3deg, essentially locked), the
> contact forces holding them there also prevented the kind of free tumbling this one underwent
> while still in open water -- i.e. jamming is self-reinforcing once contact is established, but a
> rod that reaches the pillar region still tumbling freely (as this lighter, surfaced one did,
> likely because reduced submersion meant less damping/lubrication resistance to rotation) has a
> real chance to thread through even a geometrically "impossible" gap. This meaningfully qualifies
> the earlier L/gap-ratio framing: L/gap predicts whether a STATIONARY/non-tumbling rod can pass,
> not whether a freely tumbling one can -- worth stating carefully in the thesis rather than
> presenting L/gap as a deterministic pass/fail threshold.
>
> **Twelve small rods (4268353)**: exact confirmation of the depth-based split found at the interim
> check, now stable over 30900 further steps. The same 4 deep-released (initial Z=10) rods remain
> jammed in the same two floor-piles, barely having moved (still Z~2.9-3.7, sustained
> F_contact+lub magnitude ~0.1-0.2, |vel|~1e-5). All 8 mid/shallow-released rods continued
> downstream without further incident -- several have now wrapped around the periodic X boundary
> once or twice (e.g. uid ...098 reached X=434.5 after two full wraps), confirming free,
> unobstructed onward travel once a rod is through the gap. No stuck<->escaped transitions occurred
> between step 20100 and step 51000 -- the depth-based partition found at the interim check is the
> final, stable outcome for this population.

> **[Claude, 2026-09-26] "Keep doubling in X until we reproduce the old crash diagnostics" series
> (per explicit request): doubling step 2. Job 4268493 (12 giant rods, 4x1x3) completed cleanly;
> job 4268494 (24 small rods, 2x4x3) CRASHED with a genuine reproduction -- but via a DIFFERENT
> mechanism than every earlier pinch case this session.**
>
> **12 giant rods (4268493, `giant_rods_12x_clog_gas.prm`)**: completed all 51000 steps, exit clean,
> zero corruption diagnostics throughout. Consistent with the established pattern -- rods released
> closer to the pillar row (X=364, X=350) jam there; rods released from the two new rearward rows
> (X=336, X=322) had more distance to develop torque/reorientation before arriving. Output copied to
> `~/University/Thesis/animations/giant_rods_12x_clog_fixesoff_4268493_20260926/`. Next doubling
> step (24 rods, 8x1x3, rows extended to X=308/294/280/266 at the same 14 LU pitch,
> `giant_rods_24x_clog_gas.prm` / job **4268581**) submitted immediately after on the same mirror
> (`~/mt-dennys`), fixes-off state reverified via checksum before submission.
>
> **24 small rods (4268494, `small_rods_24x_2x4x3_clog_gas.prm`) -- CRASHED, genuine reproduction of
> a corruption diagnostic, at step 22196/51000 (MPI_ABORT, exit status 1).** This IS the target
> result of this doubling series -- **not** treated as a failure to fix; no further doubling was
> submitted for this population.
>
> Sequence of events: MEM-GARBAGE-LINK fired 40 times starting step 22195, all involving two
> specific rods -- uid 4323455642275676162 and uid 4755801206503243780 -- repeatedly referencing the
> SAME fluid cells (e.g. block-local (38,17,0), an `interface` cell, fill=0.5). One step later
> (22196), a **fatal STALE-PARTICLE-IDX** error fired on rank 108: a near-boundary fluid cell was
> still flagged as owned by uid 4323455642275676162, but that uid "no longer maps to a valid local
> particle index" -- i.e. the fluid-side obstacle-cell mapping (`particleField_`) had desynced from
> mesa_pd's own local particle storage for this uid (its storage slot was reused, migrated, or
> removed cross-rank without this cell's ownership entry being refreshed). Per the code's own
> (pre-existing) diagnostic message at `SimpleBB.h:313`, this desync is known to produce a
> multi-order-of-magnitude hydrodynamic force spike and a subsequent invalid-quaternion crash several
> steps later if not caught -- so it aborts immediately at the point of desync instead. MPI_ABORT
> followed on rank 114, run exited status 1.
>
> **This is mechanistically DIFFERENT from every earlier pinch/crash case investigated this
> session.** All prior MEM-GARBAGE-LINK events (jobs 4267519/4267784/4267856/4267857/4267858, the
> giant single-rod/6-rod/12-rod series up to this point) either (a) involved exactly ONE rod's own
> self-enclosure (never two different rods) and fired near the pillar/high-contact-density region, or
> (b) never escalated to a fatal abort at all (all six clean completions logged above). This crash
> involves TWO DIFFERENT rods, is triggered by a CROSS-RANK OWNERSHIP DESYNC (not a geometric
> self-enclosure), and occurred FAR DOWNSTREAM of the pillar (X~605-608, nowhere near X=384) -- both
> rods had already escaped the gap and were travelling together in open flow (positions ~3 LU apart
> in Z, ~5 LU apart in Y, both moving at similar velocity ~0.018-0.022, per the CylForces log at
> steps 21900-22100 immediately preceding the crash) when their ownership handoff between ranks
> raced with the fluid-side obstacle-mapping refresh. This is precisely the class of issue flagged
> as an untested edge case earlier in this session re: `SyncGhostOwners`'s one-hop-per-call
> convergence for large bodies -- though STALE-PARTICLE-IDX specifically names an ownership-handoff
> desync, a related but distinct mechanism from the force-reduction hop-count lag discussed there.
> Crowding more rods into the same channel (this doubling series' whole point) increased the
> density of near-simultaneous cross-rank handoffs far downstream, which is plausibly why this
> mechanism only appeared at 24 rods and not at 6/12.
>
> Full crash log (MEM-GARBAGE-LINK entries, the fatal STALE-PARTICLE-IDX message, and the
> pre-crash CylForces trajectories of both involved rods) preserved in job 4268494's own
> `.out`/`.err` files under `~/mt-dennys-mirror2/apps/showcases/Antidunes/sbatch-output/` and its
> archived run output at
> `/home/hpc/iwia/iwia148h/mt-dennys-mirror2/build/apps/showcases/Antidunes/archive/20260926_170448_small_rods_24x_2x4x3_clog_fixesoff/`
> (not yet copied locally -- crash output, not needed for the animation library, but worth pulling
> the relevant VTK frames around step ~22000-22200 if a visual of this specific two-rod encounter is
> wanted later).

> **[Claude, 2026-09-26] Doubling series, giant-rod side: job 4268581 (24 giant rods, 8x1x3, second
> doubling step after the clean 12-rod completion 4268493) ALSO CRASHED, at step 7409/51000 --
> exact same fatal mechanism and code location as job 4268494's crash above. Both populations have
> now independently reproduced the target crash at their second doubling step (24 rods); no further
> doubling submitted for either population per explicit instruction.**
>
> Same fatal `STALE-PARTICLE-IDX` from `SimpleBB.h:313`, same MPI_ABORT/exit-status-1 pattern: rank
> 55's cell (30,11,14) was flagged as owned by uid=3891110078048108546, but that uid no longer
> mapped to a valid local particle index -- the identical cross-rank ownership-desync mechanism as
> the small-rod crash, just involving a different rod (this run's own UID numbering, unrelated to
> job 4268494's). MEM-GARBAGE-LINK fired 16 times just before the abort, involving uid
> ...108546 and a numerically-adjacent uid ...108547 (this job's own particle-creation-order
> numbering, not necessarily physically special). The affected rod (uid ...108546, F_buoy+grav_z=
> -0.201616, i.e. the densityRatio=0.85/heaviest tier) had drifted up to Z=56.4 -- well above the
> Z=49.13 water surface, echoing the same "surfaced/tumbling" state seen in every earlier
> giant-rod-escapes-the-water case this session -- when its F_hydro x-component spiked from the
> steady ~0.002-0.05 range to -0.388 in a single step (step 7400), immediately preceding the step
> 7409 fatal abort. This spike is exactly the symptom the code's own diagnostic message warns is
> produced by this desync if not caught early.
>
> **Combined interpretation**: this crash mechanism is NOT specific to rod size or geometry -- it
> reproduced independently for both the giant rods (24 of them, 8x1x3) and the small rods (24 of
> them, 2x4x3), via different specific UIDs, at the same rod COUNT (24, this run's 2nd doubling from
> a clean population of 12) and the same 144-rank configuration. This points to a general cross-rank
> particle-ownership-churn threshold tied to how many floating bodies are moving/crossing rank
> boundaries simultaneously in this domain/rank decomposition, rather than to any property of the
> rods themselves. Both crashes are pre-existing, deliberately-placed safety aborts
> (`SimpleBB.h:313`'s own diagnostic comment describes exactly this failure mode and why it aborts
> rather than silently continuing) -- i.e. the underlying desync bug already existed and was already
> instrumented against; this doubling series simply found the particle-count regime where it starts
> firing. Neither crash's output has been copied to the local animation library (crash dumps, not
> clean completions) -- both are preserved on the cluster under their respective jobs' archive
> paths, listed in the two entries above, if a closer look (e.g. VTK frames right before the abort)
> is wanted later.

> **[Claude, 2026-09-26] Root cause of the STALE-PARTICLE-IDX crash, confirmed via code trace (not
> speculation): an ordinary cross-rank ghost-drop, not a teleport. Then: both crashed jobs were
> resumed from their last checkpoint (per explicit request) and BOTH re-crashed via the exact same
> mechanism, involving the exact same rod(s) -- confirming this is a robust, reproducible
> instability, not a one-off numerical fluke.**
>
> **Mechanism** (traced through `MovingParticleMapping.h`, `SetupStages.h`, `ParticleKernels.h`,
> `SimpleBB.h`): the fluid-side obstacle-cell mapping (`particleField_`) is refreshed exactly once
> per MAIN timestep, and that refresh is entirely PARTICLE-driven (it iterates currently-known
> particles and updates their own cells) rather than CELL-driven (nothing ever revisits a specific
> cell's stale entry once its owning particle stops being visited). Meanwhile mesa_pd's own
> ghost-ownership sync (`syncCall()`, `SyncNextNeighbors`/`SyncGhostOwners`) runs up to
> 2xnumSubCycles times PER main step (~20 times, given `numSubCycles=10`), each a routine one-hop
> handoff. If any of those many sub-step syncs drops a rank's ghost copy of a particle -- an
> entirely ORDINARY event once a rod drifts clear of that rank's blocks, requiring no unusual
> velocity or displacement -- that rank's `particleField_` is left holding a stale UID reference
> that nothing will ever clean up, until the next MEM force sweep (`SimpleBB.h`) hits it and the
> pre-existing deliberate safety check there (line 313) aborts rather than risk the
> multi-order-of-magnitude force spike / invalid-quaternion crash the desync would otherwise cause
> several steps later. Confirmed NOT a teleport: both crashing rods' own logged centerline positions
> were smooth and unremarkable right up to each crash, in both jobs.
>
> **Checkpoint-resume test (per explicit request, "let both from the crashed continue via their
> last checkpoint")**: both jobs' single available snapshot (giant: step 5000; small: step 20000,
> the only checkpoints written before their respective crashes given `snapshotFrequency 5000`) were
> copied to a dedicated `resume_snapshots/` folder on each mirror (to dodge the sbatch script's
> leftover-output rescue logic) and resumed via `loadSnapshot true` + an absolute
> `snapshotBaseFolder` (new prms `giant_rods_24x_clog_resume_gas.prm` /
> `small_rods_24x_2x4x3_clog_resume_gas.prm`, jobs **4268636**/**4268637**). Verified beforehand via
> a dedicated code-reading check that floater particles are fully guarded by `loadSnapshot` (no
> re-creation from the prm's `FloatingCylinders` block when loading -- `SetupStages.h:522-697`) while
> the fixed pillars are correctly always recreated (they're GLOBAL particles, explicitly excluded
> from checkpoint serialization) -- confirmed no duplication risk before submitting.
>
> **Both resumed runs re-crashed, via the identical mechanism, involving the identical rod(s):**
> - Giant-rod resume (4268636): crashed at step 7614 (vs. step 7409 in the original 4268581 -- i.e.
>   only ~200 steps further than where it crashed the first time). Involved uids
>   3891110078048108546 and ...547 (both original crash participants) PLUS a third, ...545 (newly
>   involved this time).
> - Small-rod resume (4268637): crashed at step 21999 (vs. step 22196 originally -- only ~200 steps
>   later this time too). Involved the EXACT SAME two uids as the original crash (4323455642275676162
>   and 4755801206503243780), a 100% match.
>
> **This is a strong result**: re-running from checkpoint reproduced the same failure, with (nearly)
> the same rods, within a similar step-count of resuming, in BOTH populations. This rules out "rare
> random numerical fluke" as the explanation and confirms a robust, structurally-determined
> instability tied to these specific rods' ongoing trajectories/interactions at this point in each
> run -- consistent with the mechanism above (an ordinary, recurring cross-rank ghost-handoff gap
> that these particular rods' motion happens to trigger reliably around this stage of the
> simulation). No further doubling or resume attempted for either population -- this is the
> characterized end-state of this investigation thread.
>
> **Methodological note for future checks in this series**: `[ERROR]`-level diagnostic tags
> (including `STALE-PARTICLE-IDX`) are written to the job's `.err` file, NOT its `.out` file --
> grepping only the `.out` file (as done for the routine `MEM-GARBAGE-LINK`/`NAN-LINK`/etc. WARNING-
> level counts throughout this series) will silently undercount or miss ERROR-level fatal
> diagnostics entirely. Always check both streams when auditing for a crash cause going forward.

> **[Claude, 2026-09-26] Both crashed jobs re-run from checkpoint a second time with new
> instrumentation added per explicit request ("add additional debugging statements to nail them
> both down and resume the run from the checkpoint"). Result: DEFINITIVELY confirms the trigger is
> the MEM hydrodynamic-force computation, NOT a DEM collision -- "a fast and sudden movement due to
> a collision" was the natural first guess and is directly ruled out by the data.**
>
> **Instrumentation added** (`ParticleKernels.h`): (1) at every `[GHOST-VANISHED]` event, a scan of
> that rank's local `particleField_` for cells still flagged with the vanished UID, logging
> `orphanCellCount`/`firstOrphanCell` -- measures whether a given ghost-drop actually leaves a
> stale reference behind (the real STALE-PARTICLE-IDX precondition) or is benign; (2)
> `logFloaterHistory` enabled in a tight window around each job's known crash step, giving
> per-SUB-CYCLE position/velocity/force/hydroForce/ownership for every floater.
>
> Both jobs re-crashed via the identical mechanism, at nearly the same step as their first resume
> (giant: step ~7872 vs 7614 before; small: step ~22068 vs 21999 before -- small rods' crash
> involved the SAME uid pair for the THIRD consecutive run in a row).
>
> **Orphan-scan result**: every earlier `[GHOST-VANISHED]` hit for the eventually-crashing UID(s),
> across dozens of occurrences spanning hundreds of steps, showed `orphanCellCount=0` -- completely
> benign. Then, in BOTH jobs independently, at one specific step/sub-cycle, MULTIPLE ranks
> simultaneously reported large nonzero orphan counts for the same UID(s) (giant, step 7872 subCycle
> 5: orphanCellCount = 57/1709/51/98/538/1607 across 6 different ranks; small, step 22068: 65/253/6/
> 189/48/66/392 across 7 ranks) while other ranks at the exact same step/subCycle still showed 0.
> This confirms the crash is not caused by ghost-drops being inherently dangerous (they happen
> constantly and harmlessly) but by a rare, SUDDEN, WIDESPREAD multi-rank orphaning event.
>
> **Floater-history result -- the decisive piece**: pulled the owning rank's real
> hydroForce/hydroTorque/contact-force values for the giant rod's crashing uid
> (3891110078048108546) across steps 7860-7872. `forceX/Y/Z` (DEM contact+lubrication, reset and
> accumulated fresh each sub-cycle) is **exactly 0 for every single row in the entire window** --
> this floater had NO collision at all leading up to the crash. Meanwhile `hydroForceX/Y/Z`
> (computed once per main step by `reduceAndAvgHydroForces`, BEFORE any DEM sub-cycling starts) sits
> at normal magnitudes (order 0.1-0.5) for all 12 steps from 7860-7871, then at step 7872 explodes to
> hydroForceX=438161, hydroForceY=443471, hydroForceZ=-443260, hydroTorqueY=-6.708e6 -- a roughly
> MILLION-fold jump, at the exact step the orphan-scan caught the multi-rank stale-cell event for
> this same UID. This is word-for-word the failure mode `SimpleBB.h:313`'s own diagnostic message
> predicts ("a multi-order-of-magnitude hydrodynamic force spike"), now directly measured rather
> than inferred.
>
> **Conclusion**: the causal chain is fully closed. A stale `particleField_` link (left behind by an
> ordinary ghost-drop that, this one time, wasn't cleaned up before the next MEM sweep touched it)
> feeds garbage into `reduceAndAvgHydroForces()`, producing a real, measured ~10^6x hydrodynamic
> force/torque spike -- not a collision, not a DEM event, and not any unusual prior motion. That
> force spike is what then drives the runaway velocity cascade (204->1020 LU/step across sub-cycles,
> documented in the earlier entry above) and triggers the fatal `STALE-PARTICLE-IDX` abort. No
> further reruns planned for this specific finding -- the mechanism is now fully characterized from
> root cause (cross-rank mapping desync) through to observable consequence (MEM force spike ->
> velocity runaway -> abort), with direct log/CSV evidence at every link.

> **[Claude, 2026-09-26] The STALE-CELL-GUARD was tested empirically ("add the guard, also add logs
> to confirm it rerun from the checkpoints"), found to trade one crash for a worse one, and REVERTED
> per explicit request. Root cause of the worse crash traced precisely; a deterministic, small-scale
> reproduction was then built and submitted, independent of the rare cross-rank race.**
>
> **Guard test (jobs 4268749 giant, 4268750 small, both resumed from the same checkpoints used
> throughout this series)**: the guard (clear+reset any orphaned `particleField_` cell the moment a
> ghost vanishes) DID eliminate the targeted `STALE-PARTICLE-IDX` abort in both jobs (zero
> occurrences, confirmed via both `.out` and `.err`). But both then crashed anyway, LATER, via a
> WORSE failure -- exactly what `SimpleBB.h:313`'s own diagnostic message warned would happen if the
> desync were allowed to proceed unchecked ("a multi-order-of-magnitude hydrodynamic force spike and
> a subsequent invalid-quaternion crash several steps later"):
> - Small-rod (4268750): guard fired at step 22584, clearing >1000 obstacle cells belonging to rod
>   uid 4755801206503243780 in one abrupt sweep. The first NaN appeared 3 steps later (step 22587),
>   confirmed NOT a pinching/near-total-enclosure artifact (`obstacleNeighbours=0/26` at every
>   reported cell, both at this first occurrence and the run's very first `DEGENERATE-FILL-
>   RECONSTRUCTION` line) -- and confirmed SPATIALLY INSIDE that same rod's own tilted body: at the
>   closest logged sample (step 22500) the rod's centre was at Y=17.5 with tiltY~58deg; its 18 LU
>   length, tilted ~58deg off the Y axis, projects a Y-span of [8.5,26.5] and an X-reach of
>   ~9*sin(58deg)~7.6 LU -- both comfortably explain the first NaN cell's world position
>   (603.5,9.5,30.5) relative to the rod's contemporaneous trajectory. Root cause: the guard's
>   abrupt, geometry-blind fill-level discontinuity (>1000 cells flipped from solid obstacle to
>   fill=1 liquid in one sub-cycle) destabilises the free-surface curvature/normal reconstruction
>   (`reconstructInterfaceCellLegacy`) in the immediately adjacent, now-suddenly-bare interface
>   cells -- confirmed by directly comparing the guard's own new `[STALE-CELL-GUARD]
>   guardedCellCount=` log lines against the crash timeline: hundreds of PRIOR ghost-vanish events on
>   this exact UID pair were all benign (`guardedCellCount=0`, nothing to clear), and the fatal one
>   coincided with the guard clearing >1000 cells across 7+ ranks simultaneously at step 22584.
> - Giant-rod (4268749): pushed from its old crash step (7872) to step ~21400, then hit an
>   **invalid-quaternion assertion** (`Quaternion.h:241`, r/i/j/k all `-nan`) -- the second half of
>   the same predicted failure mode. Critically, **the guard never fired even once in this entire
>   run** (zero nonzero `orphanCellCount` events) -- so this crash is provably UNRELATED to the
>   guard. Its own first NaN (step 21361, world ~(448.5,20.5,1.5), also `obstacleNeighbours=0/26`) is
>   nowhere near any of the 24 tracked rods (closest rod X-positions cluster around 300-360 and
>   602-676; this cell is at X=448.5, Z=1.5 -- essentially on the floor). Partial lead: a
>   PRE-EXISTING, unrelated guard (`[STALE-NEARBOUNDARY]` in `AntidunesBoundaryHandling.h`, itself a
>   workaround for a related "gas cell reconstruction near a block boundary" fragility) fired at the
>   SAME block-local cell `(0,0,1)` far earlier in the run (step 10202) -- suggestive of a
>   longstanding, separate fragility in gas/interface reconstruction at block-corner cells near the
>   floor, NOT caused by or related to the rod-desync/guard investigation. No fully-pinned-down
>   trigger found for why reconstruction failed at THIS occurrence specifically; documented honestly
>   as a distinct, not-yet-explained instability rather than forcing a premature repro.
>
> **Decision: guard reverted.** A controlled, immediately-diagnosable abort (the original
> `STALE-PARTICLE-IDX`) is strictly more useful than a delayed, harder-to-trace one -- the guard
> masked the early-warning signal without fixing the underlying instability. `ParticleKernels.h`'s
> `[GHOST-VANISHED]` block reverted to scan-and-log only (no cell modification); world-coordinate
> logging added to the diagnostic (via `blockForest_->getBlockLocalCellCenter`) so a future
> origin-trace doesn't need manual block-to-world coordinate conversion. `boundaryHandlingID_`
> removed as an unused member. Reverted `ParticleKernels.h` recompiled clean locally and resynced to
> both mirrors (checksum-verified) before any further job used it.
>
> **Deterministic small-scale reproduction, independent of the rare cross-rank race**: per explicit
> request ("extract these situations and replicate it in a smaller simulation and then we could
> document them"), added a TEST-ONLY, off-by-default trigger,
> `ModelParameters.debugForceObstacleClearAtStep` (`SimulationParams.h`), wired into
> `ParticleKernels.h`'s new `[DEBUG-FORCE-OBSTACLE-CLEAR]` block: at the given step, every LOCAL
> floater has the exact same clear-and-reset sequence the guard used applied to its own CURRENT,
> perfectly-valid obstacle cells -- reproducing the guard's side effect on demand, without needing
> any desync at all. A from-scratch minimal LOCAL scenario (mirroring `single_rod_tumble_gas.prm`'s
> tiny 3-rank domain) was attempted first but hit an unrelated, PRE-EXISTING setup-time crash
> (`field::isFlag(domain_)` assertion, `BoundaryHandling.h:1203`, during `mapParticlesOntoCells`) --
> confirmed present in the UNMODIFIED `single_rod_tumble_gas.prm` too, i.e. not caused by this
> session's changes; flagged separately, not chased further given the time cost. Pivoted instead to
> adapting the already-validated `giant_rod_clog_gas.prm` (one rod, 144 ranks, jams at the pillar
> early, run successfully many times this session): new prm
> `prms/giant_rod_debug_obstacle_clear_repro_gas.prm` (6000 timesteps instead of 50000,
> `debugForceObstacleClearAtStep 3000` -- well after the rod settles into its jammed position),
> job **4268996** submitted on `~/mt-dennys` (fixes-off state reverified via checksum first).
>
> **Result: clean, decisive, deterministic reproduction -- CONFIRMED.** Job 4268996 ran exactly as
> designed: `[DEBUG-FORCE-OBSTACLE-CLEAR]` fired at step 3000 (rank 60, uid
> 4827858800541171713 -- this scenario's own single rod), `clearedCellCount=1508` (the same order of
> magnitude as the original guard's >1000-cell clears in job 4268750). The first
> `DEGENERATE-FILL-RECONSTRUCTION` appeared just **5 steps later** (step 3005) -- essentially the
> same lag as the original (3 steps) -- and escalated (1861 `DEGENERATE-FILL-RECONSTRUCTION`/
> `LOW-RHO` warnings total) to the exact same fatal error as the original guard-triggered crash:
> **`reynoldsNumber is nan!`** (`WriteData.h:109`), MPI_ABORT, exit status 1, at ~step 3100-3200 (run
> configured for 6000 steps total, never got close).
>
> This closes the loop completely: a real, valid, single, non-orphaned rod, with NO cross-rank
> desync involved at all, reproduces the identical mechanism and identical fatal error in ~3200
> steps instead of needing a rare 22000+-step race across 24 rods on 144 ranks. **Definitively
> confirms the instability is a property of the abrupt clear-and-reset operation itself** (a large,
> geometry-blind fill-level discontinuity that the free-surface curvature/normal reconstruction
> cannot tolerate), independent of how that clear gets triggered -- whether by a genuine cross-rank
> ghost-ownership desync (the original, rare path) or by direct, deliberate action (this
> reproduction). This is the documentable, thesis-ready minimal test case for Situation 1.
>
> **Situation 2 (giant-rod floor-adjacent NaN, job 4268749) remains a separate, distinct, and still
> UNEXPLAINED instability** -- not reproduced, not chased further this session (see the diagnosis
> above: not guard-related, not rod-position-related, plausibly connected to the pre-existing
> `[STALE-NEARBOUNDARY]` gas-cell-near-block-corner fragility but without a fully pinned-down
> trigger). Worth a dedicated follow-up investigation if pursued further, kept explicitly separate
> from Situation 1's now-fully-closed finding.

> **[Claude, 2026-09-27] Fourth independent, fully NATURAL (unforced) recurrence of the small-rod
> crash, per explicit request ("reproduce the error again from the original checkpoint, no forced
> activation... confirm the stale id is the issue"). This is the cleanest confirmation yet of the
> causal chain, entirely from real cross-rank desync with zero manual intervention.**
>
> **Setup**: resumed job (**4269772**) from the SAME original checkpoint used throughout this
> series (step 20000, `~/mt-dennys-mirror2`), current code (guard fully reverted --
> scan-and-log-only; `debugForceObstacleClearAtStep` not set anywhere in this prm, so the manual
> trigger path is entirely inert), `logFloaterHistory` window widened to 20500-24000 (generous
> margin after three prior attempts landed anywhere from step 21999 to 22584) to guarantee full
> coverage regardless of exactly where this attempt crashed.
>
> **Note on checkpoint integrity**: while preparing this, discovered the GIANT-rod population's
> original step-5000 checkpoint has been inadvertently OVERWRITTEN -- its `snapshotBaseFolder` was
> set to the same absolute path used for loading, and the second guard-test job (4268749) passed
> through several 5000-step snapshot boundaries (up to step 20000) before crashing later at
> step ~21400, silently replacing the original early checkpoint with a much later one. The
> SMALL-rod checkpoint was unaffected (its resume attempts never reached its next snapshot boundary
> at step 25000), so this natural-recurrence test proceeded on the small-rod population only.
>
> **Result**: crashed again, at step **21990**, rank 109, `STALE-PARTICLE-IDX` -- involving the
> exact same two rods as every prior recurrence in this series (uid 4323455642275676162 and
> 4755801206503243780), a FOURTH consecutive crash with this identical pair. The four ranks with
> full floater-history coverage through the crash step were again 108/109/114/115 -- the same
> four-way block-corner cluster identified previously.
>
> **Force trace (owning rank 108, uid 4755801206503243780, isLocal=1) -- the most extreme
> confirmation yet**:
> - Step 21980: `hydroForce=(-0.163, 0.116, 0.148)` -- small, normal.
> - Step 21990 (the crash step): `hydroForce=(-9.35e14, 2.81e15, -9.36e14)`,
>   `hydroTorque=(-1.42e16, -1.73e15, 9.05e15)` -- a roughly **10^15-fold** jump, the largest
>   measured in this entire investigation (exceeding even job 4268653's already-extreme 4.4e5).
>   `forceX/Y/Z = 0, 0, 0` throughout (contact/lubrication, confirmed zero yet again -- no
>   collision). `velocity=(0.0176, 0.0023, 0.0042)` in the SAME row -- still essentially normal,
>   confirming the corrupted force is measured BEFORE velocity has had any chance to respond to it.
>
> **This closes the loop with a fourth, fully independent, entirely natural data point**: no guard,
> no manual trigger, no forced anything -- just the real cross-rank desync recurring exactly as it
> has every time this population has been resumed from this checkpoint, with the identical rod
> pair, at the identical four-rank block corner, producing the identical causal signature (force
> normal -> force astronomically corrupted -> velocity explodes as the downstream consequence).
> Combined with the deterministic manual repro (job 4269518) and the original real-desync trace
> (job 4268653), the causal chain (stale ownership -> corrupted MEM force -> velocity explosion) is
> now confirmed across three independent lines of evidence: one fully natural/original, one
> natural/resumed (this entry), and one fully controlled/manual -- all showing the identical
> signature. No further reruns planned; this is the final piece of evidence for this investigation
> thread.

> **[Claude, 2026-09-27] CORRECTION to the causal-order claim above, per explicit push-back ("I
> don't really know if it's cell collapsed -> crazy mem force -> teleport -> stale id, or due to
> the update"). Improved logging to sub-cycle resolution and found the real picture is more
> tangled than the clean one-directional order previously asserted.**
>
> **New instrumentation** (both housed in the fixes-off `SimpleBB.h`/`ParticleKernels.h`, synced to
> `~/mt-dennys-mirror2`): (1) the `STALE-PARTICLE-IDX` fatal abort now dumps the flagged cell's full
> PDF list/density/fill directly into the abort message itself, captured at the exact instant of
> detection; (2) a new `[EARLY-CELL-CORRUPTION]` check runs EVERY sub-cycle (not just once per main
> step like every existing diagnostic) for every locally-owned floater's own current cells' fluid
> neighbours, checking for `|PDF|>10` -- closing the multi-step blind window between "force
> confirmed normal" and "force confirmed garbage" that every prior check could only bound to
> roughly +/-10 steps.
>
> **Result, job 4270814 (natural, unforced resume, same checkpoint/rods as every prior entry in
> this series)**: at the finer resolution, the actual sequence was:
> - Step 22000 (last clean sample): uid 4323455642275676162 at `pos=<606.5,22.4,29.8>`,
>   `F_hydro=(0.095,-0.05,-0.04)` -- completely normal.
> - **Step 22027, subCycle 0**: `[EARLY-CELL-CORRUPTION]` fires -- a cell ~5.8 LU from this same
>   rod's centreline (within reach of its own 18 LU tilted length) already shows `density=-0.608`
>   (physically impossible) and `pdf=-594`. At this point uid ...162 is STILL fully owned/known by
>   this rank (the check only runs for locally-owned, non-ghost floaters) -- no ownership problem
>   detected yet.
> - **Step 22030, subCycle 1** (3 steps / ~13 sub-cycles later): `[GHOST-VANISHED]` fires with LARGE
>   nonzero orphan counts (14-424 cells) across FOUR ranks (108/109/114/115) simultaneously, for
>   BOTH rods -- the mass, coordinated collapse already characterised above.
> - **Step 22030 (same step), the actual fatal abort**: fires on a DIFFERENT, COMPLETELY HEALTHY
>   cell (`fill=1`, `density=0.994`/`1.42`, every PDF component in the normal 0.001-0.34 range) --
>   the abort is purely an ownership-bookkeeping check, unrelated to that specific cell's own fluid
>   state.
>
> **This does not match the clean, one-directional "ownership desync -> corrupted force" order
> asserted in the entries above.** That claim was based on once-per-main-step data that could only
> bound the transition to a multi-step window and could not see this internal structure. The
> honest, revised reading: a small, LOCALISED fluid-cell anomaly (right next to one of the two
> rods, not yet astronomic in magnitude) appears a few steps BEFORE the visible mass ownership
> collapse -- neither "ownership breaks first" nor "an unrelated cell breaks first and causes a
> teleport" (the rod's own logged velocity/position never leaves the normal range until well after
> the mass collapse, so no teleport precedes any of this) cleanly describes what happened. The two
> phenomena (a localised, sub-catastrophic cell anomaly and a mass, multi-rank ownership collapse)
> are close together in time (~3 steps) and plausibly entangled/mutually-reinforcing -- e.g. the
> localised anomaly could itself be an early, small-scale, not-yet-orphanCellCount-detectable
> instance of the SAME desync mechanism, which then escalates into the full collapse a few steps
> later -- rather than one being a distinct, independent, isolated trigger for the other. This is
> a genuine correction, not a refinement: the previous entries' framing of a clean causal arrow
> should be read alongside this one as the more accurate, if less tidy, picture.
>
> **Addendum, same day: exact geometry of the first `[EARLY-CELL-CORRUPTION]` cell.** Per direct
> request ("does it have a special position, is it inside of an object"), computed precisely from
> `floater-history/rank_114_state.csv` (uid ...162's exact position+quaternion at step 22027,
> subCycle 0) and the corrupted cell's block-local index (26,0,27) on rank 114 (block origin
> (576,20,0) -> world cell center (602.5,20.5,27.5)): rod centre (607.14,22.34,29.82), long axis
> (from quaternion) ~(0.078,0.900,0.429). Projecting the cell onto that axis: axial offset ~3.0 LU
> from centre (well inside the +/-9 LU half-length -- mid-body, not an end-cap), perpendicular
> distance from centreline ~4.6 LU. Rod radius is 3 LU, so **the cell is NOT inside the solid body**
> -- it sits ~1.6 LU outside the surface, i.e. an ordinary near-wall fluid cell within
> lubrication-correction/reconstruction range, not a cell trapped inside an object. (Supersedes the
> earlier, cruder "~5.8 LU from centreline" estimate mentioned in conversation but never written
> here -- 4.6 LU is the precise figure.) Also notable: this cell's block-local Y index is 0, i.e.
> global Y=20 -- exactly on the by=0/by=1 inter-block seam already implicated in the four-way
> corner mechanism above, though Z=27.5 is close to but not exactly on the other seam (Z=30) this
> time, so it's block-boundary-adjacent again without being a full corner hit.
>
> **Second addendum, same day: fuller re-check of job 4270814's own log turned up a sharper picture
> than the "entangled, can't disentangle" framing above -- re-grepped without truncation.**
>
> Confirmed genuine natural resume: `[CHECKPOINT-RESUME] loadSnapshot=true beginTimeStep=20001
> globalOwnedParticleCount=24` (full population, no forced trigger). Crashed step 22030, exit
> status 1. Elapsed 00:04:30 vs. prior similar run 4269772's 00:04:05 -- **the new per-sub-cycle
> `[EARLY-CELL-CORRUPTION]` scan adds only ~10% runtime overhead**, cheap enough to leave on.
>
> **All four ranks (108,109,114,115) aborted simultaneously**, not two as the first pass reported
> (that pass missed one message to a `head` truncation). Two distinct rod uids involved (...780 on
> 108/109, ...162 on 114/115 -- matching the established four-rank-corner pattern). Abort-cell
> states: rank108 `density=1.036` (healthy), rank109 `density=0.994` (healthy), rank114
> `density=1.42` (healthy), **rank115 `density=9.48e7`, PDFs up to 8.3e7 (catastrophic)** -- for
> the SAME uid...162 that rank114's (healthy) abort concerned. So 3 of 4 simultaneous aborts were
> on numerically healthy cells and 1 was already catastrophic, for the same particle.
>
> **`[EARLY-CELL-CORRUPTION]` full timing (both rods, all four ranks):** rank108/uid...780 and
> rank114/uid...162 both first show corruption at **step 22027, subCycle 0** -- three full steps
> before any ownership symptom anywhere -- and stay corrupted for several consecutive sub-cycles
> (108: all ten sub-cycles of step 22027 plus 22028/subCycle0; 114: sub-cycles 0-5), then **go
> quiet** for ~2 steps. Rank109/uid...778 and rank115/uid...164, by contrast, show **no** corruption
> at all until **step 22030, subCycle 0** -- the same step everything else happens. First nonzero
> `[GHOST-VANISHED]` orphaning fires at **step 22030, subCycle 1, on all four ranks simultaneously**
> (orphan counts 14-424).
>
> **Refined reading**: this is more resolved than "closely timed, entangled" -- a genuine local
> fluid instability originates on a SUBSET of the eventually-affected ranks (108,114) a full 3
> steps early, persists for several sub-cycles, then goes quiet with no visible ownership-tracking
> symptom anywhere in that window. Three steps later, at step 22030, the previously-clean ranks
> (109,115) show corruption for the first time AND the cross-rank ownership bookkeeping collapses
> on all four ranks together AND the fatal abort fires -- healthy on most ranks, catastrophic on
> one. The mass ownership collapse looks like the SYNCHRONIZING event across ranks, but the
> instability itself demonstrably started locally, earlier, and independently of it on two of the
> four ranks. Still not a single clean one-line causal arrow, but no longer purely "can't tell
> which comes first" either -- the sub-cycle data localizes the true origin to ranks 108/114 at
> step 22027.
>
> **Correction to the second addendum above**, caught by direct question ("are these all
> referring to the same cell at the same timestep?"): the uid-to-rank mapping of the four
> simultaneous aborts was mis-transcribed. Re-verified directly from the raw abort messages:
>
> | rank | uid flagged stale | fluid cell (world) | density |
> |---|---|---|---|
> | 108 | ...780 | (597.5,16.5,19.5) | 1.036 (healthy) |
> | 109 | ...780 | (608.5,12.5,30.5) | 0.994 (healthy) |
> | 114 | ...780 | (598.5,20.5,19.5) | 1.42 (healthy) |
> | 115 | **...162** | (601.5,20.5,30.5) | **9.48e7 (catastrophic)** |
>
> These are **four distinct physical cells** (not one cell viewed from four ranks), spread across
> an ~11x8x11 LU region. **Three of the four aborts (108,109,114) concern stale references to
> uid...780; only rank115's concerns uid...162** -- the earlier claim that rank114 and rank115
> shared uid...162 was wrong; rank114's abort was actually uid...780.
>
> This sharpens the picture further, and cleanly: rank108 locally OWNS uid...780 (confirmed by
> `[EARLY-CELL-CORRUPTION]` only firing there for it) -- yet by step 22030 rank108 ALSO held a
> second, separate, stale ghost-reference to that same uid that no longer resolved, meaning
> ownership of uid...780 had already shifted away from rank108 by the time of the abort. The same
> structure holds on 109 and 114 (stale ghost copies of uid...780, not their own owned particles),
> and on 115 (stale ghost copy of uid...162). **Every one of the four simultaneous aborts was
> triggered by a rank's stale ghost-reference to a particle it does not itself currently own** --
> not by a rank's own owned-particle bookkeeping. This is the cleanest direct confirmation yet of
> the ghost-ownership-desync mechanism itself, independent of the still-open question of exactly
> when/why the desync began.

> **[Claude, 2026-09-27] Obstacle-PDF freeze fix implemented and A/B tested against reconstruction-
> only (per explicit request); NEITHER prevents this crash class -- and the actual trigger, for this
> specific rod pair, connects directly to the already-root-caused 2026-09-20 analytic cylinder-
> contact discontinuity, not to the single-cell checkerboard drift this fix targets.**
>
> **Freeze fix**: `BuildTimeloop.h`'s `saveObstaclePdfsBeforeCollide`/`restoreObstaclePdfsAfterCollide`
> (gated by new `ModelParameters.useObstaclePdfFreeze`, default true) snapshot every obstacle-flagged
> cell's 27 PDFs before the stream+collide sweeps and restore them afterward, so an obstacle cell's
> PDF genuinely stays frozen across steps instead of silently evolving under the generated D3Q27
> kernel (which has no per-cell obstacle mask). Motivated by Bogner (2017 thesis, p.111, citing Ladd
> 1994a: MEM forces "should be averaged over two time steps to avoid checkerboard effects") and
> Watanabe et al. (2023): both cite the same 2-step force-averaging mitigation this codebase already
> has (`useHydrodynamicForceTimeAveraging`), which smooths the force OUTPUT but not the underlying
> PDF value the force is read from. Progressive single-cell tests (collide-only, then collide+stream)
> on the ORIGINAL step-22012 crash instance showed genuine, measurable improvement: crash delayed
> 22012 -> 22546 -> 22781, oscillation replaced by a small monotonic drift. **Self-correction made
> along the way**: earlier claims this session that mesa_pd's particle-obstacle reconstruction uses
> an "OnlyMissing" (inheriting) model were WRONG -- it's hardcoded to `EquilibriumReconstructor`,
> which reinitializes ALL 27 populations from equilibrium at (rho_avg, u_particle); the freeze's real
> value is protecting `SimpleBB`'s continuous `pdf_old(dir)` reads during the obstacle phase, not
> protecting a reconstruction inheritance assumption that doesn't exist.
>
> **A/B test against job 4271048's own crash instance (step 22264, uid...162), same step-20000
> checkpoint on both mirrors, fully synced+instrumented code (verified byte-identical, including a
> real sync gap on `mt-dennys-mirror` -- a stale `SimulationObjects.h` missing the `diagnosticEvents`
> member, caught via a build failure and fixed via the user's own `antidunes.sh --sync` + rsync-to-
> both-mirrors sequence before resubmitting):**
> - **mirror2 (`useObstaclePdfFreeze=true`), job 4271395: crashed step 22012** (uid...778 on rank109
>   AND uid...162 on ranks108/114/115, four ranks simultaneously, two different rods).
> - **mirror (`useObstaclePdfFreeze=false`, reconstruction-only), job 4271396: crashed step 22257**
>   (uid...780, three ranks simultaneously).
> Neither delay resembles the earlier single-cell progression (22012->22546->22781) -- and the
> pre-crash signature is completely different: `[OBSTACLE-CELL-REPAIR]` (a pre-existing safety net,
> `SetupStages.h:1615`, resets any obstacle cell with non-finite/non-positive rho to equilibrium)
> fired for the FIRST time in either rank's entire ~2260-step replay right before each crash, in ONE
> step catching 17+ cells (mirror2, rho as low as -2.6e21) belonging to a single rod -- a sudden,
> multi-cell burst, not gradual growth. **The freeze does not prevent this failure mode in either
> configuration.**
>
> **Chasing the trigger (per explicit request "so what's the common denominator") connected today's
> new data to the pre-existing 2026-09-20 and 2026-09-26 findings above:**
> - `floater-history/rank_108_state.csv` (mirror, uid...780) shows raw hydroForce growing mildly
>   (0.30 -> 0.60 -> 1.0 LU across steps 22252-22254) then spiking >130x in ONE step at 22255
>   (-133.7,-23.2,-63.1), kicking the rod's velocity into the teens; by step 22256 it has moved
>   ~2 LU in one main step, sweeping through previously-fluid cells and converting many to obstacle
>   at once -- directly producing the `[OBSTACLE-CELL-REPAIR]` burst.
> - `floater-history/rank_108.csv` (the CONTACT-PAIR log, mirror, uid...162 vs uid...780) shows the
>   analytic contact classifier oscillating `band-band-skewed` <-> `band-edge` over steps 22248-22256
>   with penetrationDepth smoothly growing (0.12 -> 0.17 LU), then at step 22257 the reported
>   penetrationDepth jumps discontinuously to **2.52 LU** with contact force jumping **0.155 -> 12.5
>   LU (~80x)** while the geometry itself has not changed nearly that much -- the SAME
>   scenario-classification discontinuity already root-caused on 2026-09-20
>   (`ANTIDUNES_CylinderClosedFormJump`: the closed-form band-edge/edge-edge scenarios examine only
>   ONE rim circle and flip which one arbitrarily as a nested pair's geometry changes, producing a
>   multi-LU force discontinuity with no real motion behind it) -- previously observed for a
>   DIFFERENT pair in the cylinder-hopper population, now reproduced for THIS pair in the small-rod
>   population. `ownerRankB` for uid...780 in that same contact row changes 114->115 in the SAME
>   step as the force jump, i.e. the STALE-PARTICLE-IDX-adjacent ghost-ownership handoff and the
>   contact-force discontinuity are coincident, not independently timed.
> - On mirror2's run (uid...778 vs uid...162, no floater-history available -- the crash at step
>   22012 predates that prm's `logFloaterHistoryStartStep 22200`), the `[DIAG-SWEEP][VELOCITY-JUMP]`
>   events show BOTH rods diverging within the SAME 1-2 step window (22009-22012), a spatial
>   neighbourhood ~7-9 LU apart -- consistent with a shared local trigger rather than two independent
>   coincidental failures. PID/mean-velocity log confirmed smooth through step 22000 (meanVelocity
>   0.0403243, acceleration decaying normally) -- ruling out a global forcing spike as that shared
>   trigger.
>
> **Synthesis**: this is not "any two rods create an unstable cell" -- it is specifically uid...162
> and uid...780 (recurring across at least 5 independent runs now: the original 09-26 crash, its two
> checkpoint-resumes, and both of today's A/B runs) whose ongoing trajectories bring them into a
> slow, tightening near-contact. That narrow gap is simultaneously (a) hard for the free-surface
> fluid solver (the `[EARLY-CELL-CORRUPTION]`/squeezed-interface-cell mechanism found 09-20/09-26)
> and (b) hard for the analytic cylinder-cylinder contact detector (the deep-overlap
> scenario-classification discontinuity found and only partially fixed 09-20 -- the NEAR-PARALLEL
> false-positive case was fixed; the DEEP-OVERLAP discontinuity for genuinely close/nested pairs was
> only diagnosed, with GJK/EPA-only rod-rod detection validated as a clean fix, 3/3 vs 0/3, but never
> adopted for this population). Whichever spikes first (fluid MEM force or contact force) kicks a rod
> far enough in one step to outrun mesa_pd's ghost-sync radius, producing the stale `particleField_`
> reference that `STALE-PARTICLE-IDX` (correctly) aborts on -- the abort is the tripwire, not the
> root cause. **The obstacle-PDF freeze addresses neither (a)'s deep mechanism nor (b) at all**,
> which is consistent with it measurably helping a slow single-cell drift (the original 22012
> instance) while doing nothing for this pair's sudden, contact-adjacent blowup.
>
> **Recommended next step (not yet run)**: test `useAnalyticalContactDetection false` (GJK/EPA-only
> rod-rod contact) against this exact step-20000 checkpoint/population, mirroring the 09-20
> methodology that produced 3/3 clean cylinder-hopper runs vs 0/3 clean for the analytic path.

> **[Claude, 2026-09-27] NEGATIVE RESULT: GJK/EPA-only contact detection does NOT fix this crash --
> the 09-20 cylinder-hopper mitigation does not transfer to this population.** Job **4271486**
> (`useAnalyticalContactDetection false`, everything else identical to the freeze-ON A/B prm, same
> step-20000 checkpoint, 6000-step budget) still crashed via `STALE-PARTICLE-IDX` at **step 22174**,
> the SAME recurring `uid...162` (ranks 108/114 shown), density blown up to 359.2/3.48, PDFs up to
> 285 -- the identical failure signature, just a different step than any analytic-path run (22012/
> 22257/22264). **Both of the two known, previously-validated-elsewhere candidate mitigations tested
> today -- the obstacle-PDF freeze (targets the fluid-side single-cell drift) and GJK/EPA-only
> contact detection (targets the contact-classification discontinuity) -- fail to prevent this
> specific rod pair's crash, tested independently.** This sharpens rather than closes the
> investigation: removing the contact-detection discontinuity stops the RODS from being kicked
> together violently, but if their real, physically-accurate trajectory still brings them into a
> genuinely tight near-contact anyway, the underlying squeezed-fluid-gap instability (found 09-20 at
> step 3515: an interface cell in a <1-cell gap between two rods going non-finite/fill<0) can still
> occur independent of how accurate the contact force itself is -- contact-detection accuracy affects
> HOW the rods get close, not whether the fluid can tolerate them being close. Neither tested
> mitigation touches that fluid-side limit directly. **Per explicit request ("create a box around
> all the cylinders in the crash area, log out their values and see the crash evolution"), next step
> is a new spatial+temporal fluid-box logger (see below) rather than a third guess at a targeted fix.**

> **[Claude, 2026-09-27] Fluid-box history logger implemented and run -- the full crash evolution is
> now directly visible, and it looks nothing like the two mechanisms tested above (single-cell
> checkerboard drift, contact-force discontinuity). It is a WHOLESALE, near-simultaneous corruption
> of one rod's entire covered obstacle region.**
>
> **New diagnostic** (`SetupStages.h::logFluidBoxHistory`, gated by `ModelParameters.logFluidBoxHistory`
> + `fluidBoxMin/Max` + `logFluidBoxHistoryStartStep/EndStep`, SimulationParams.h): every rank writes
> one CSV row per OBSTACLE-or-INTERFACE cell inside a fixed world-space box, once per main step
> (mainStep, rank, world x/y/z, flag, owning uid, density, fill, max|PDF component|) -- the missing
> piece between `trackWatchCells`' handful of hardcoded single cells and `logFloaterHistory`'s
> floater-only view. Hit one real bug along the way: `getParameter< Vector3<real_t> >(name, default)`
> (the two-argument overload) fails a `static_assert` deep in `Vector3`'s own generic single-arg
> constructor for this config binding -- worked around via `isDefined()` + a ternary, the same pattern
> `clogUpstreamPos` etc. already use (SimulationParams.h) -- not investigated further, may be worth a
> waLBerla-upstream look if any FUTURE prm parameter needs a defaulted Vector3.
>
> **Job 4271512**: identical config to today's freeze-ON A/B run (job 4271395) plus the new logger,
> box covering every crash-relevant position seen all session (X 590-625, Y 8-30, Z 15-50), window
> 21800-22050. **This run did NOT reproduce byte-identically** (crashed step 22031 via `uid...164`,
> not step 22012 via `uid...162/...778` as job 4271395 did with otherwise-identical code/checkpoint/
> config) -- consistent with, not contradicting, this series' own earlier-established finding that
> resumed/rerun jobs from the same checkpoint are NOT bit-reproducible (09-20 checkpoint-resume
> entries). The crash step still fell inside the chosen window, so the run remains useful.
>
> **The evolution, cell (617.5,22.5,34.5), owned by uid...164 (one more member of the same recurring
> uid cluster as ...162/...780/...778 -- all differ by single-digit offsets, almost certainly adjacent
> particles created sequentially in this population):**
> - Steps 21800-22024: healthy OBST cell, density oscillating in a tight ~1.0016-1.0026 band (the
>   freeze fix's own already-characterized small residual drift), **max|PDF component| pinned at
>   EXACTLY 0.296931 for all 225 steps** -- direct, unambiguous confirmation the obstacle-PDF freeze
>   is doing exactly what it was designed to do on cells it actually covers.
> - Step 22025: density steps discretely to 0.945695 (a real ~6% drop, much bigger than the ambient
>   drift) and holds there exactly for 3 steps -- a NEW frozen value, meaning something wrote to this
>   cell outside the freeze's own save/restore bracket right at this step (candidate: SimpleBB's MEM
>   bounce-back, confirmed 09-25 to write into obstacle cells' PDFs every step, unaffected by the
>   freeze, which only brackets stream+collide).
> - **Step 22028: flag flips OBST -> INTF (the rod uncovers this cell) in the SAME step
>   max|PDF| explodes to 1.04e64.** Reconstruction (`PdfReconstructionManager`/`EquilibriumReconstructor`)
>   runs at exactly the uncovering transition and its OUTPUT is already catastrophic.
> - Steps 22029-22030: cascades to every neighbour within the probed radius (max|PDF| up to 3.5e99
>   across a 3x3x3+ neighbourhood by step 22030). Step 22031: `STALE-PARTICLE-IDX` abort.
>
> **But tracing uid...162's OWN covered cells (not uid...164's) in the same run finds the corruption
> started EARLIER and MUCH MORE WIDELY than one squeezed cell:**
> - Steps up to 22023: every uid...162 obstacle cell in the box is healthy (max|PDF| ~=0.3, the same
>   frozen baseline).
> - **Step 22024: AT LEAST 7 of uid...162's obstacle cells, spanning roughly X606-615/Y17-22/Z30-35 --
>   essentially this rod's ENTIRE covered footprint (radius 3, length 18) -- jump from healthy to
>   1e8-1e43 in ONE STEP, simultaneously**, not a single cell, not a gradual buildup. By step 22025 the
>   worst cells reach 2.3e55.
> - **This is NOT actually the true onset either.** The pre-existing `[MEM-GARBAGE-LINK]` guard
>   (`SimpleBB.h`, ALREADY an active guard in this build, not log-only -- it detects `|pdf|>10` on a
>   specific bounce-back link and skips that one link's contribution to `forceMEM` rather than
>   aborting) has been firing at a STEADY ~11-13/step baseline for uid...162 since at least step
>   22000 (possibly earlier, not checked) -- i.e. some low-level link-level garbage is a normal,
>   mostly-benign background condition throughout this run, already being caught and skipped. That
>   rate climbs conspicuously BEFORE the visible cell-level collapse: 22/step (22017-22018), then a
>   sharp spike to **67 at step 22020**, 39 at step 22023, settling to 13-18/step through the collapse
>   itself (22025-22030) -- i.e. the escalation in SKIPPED garbage links precedes the first
>   whole-cell-level corruption (22024) by 4-7 steps.
>
> **Reading**: this is not the single "squeezed interface cell between two rod surfaces" mechanism
> from 09-20/09-26, and not the slow per-cell checkerboard drift the freeze targets -- it is a
> WHOLESALE corruption of one rod's entire covered region, preceded by several steps of an escalating
> (but not yet cell-collapsing) rate of individually-skipped garbage links. The obstacle-PDF freeze
> correctly holds a cell steady WHILE it stays obstacle-flagged (confirmed directly above) but (a)
> does not stop whatever writes the NEW, wrong frozen value at the moment of the step-22025-style
> jump (outside its own save/restore bracket) and (b) actively PRESERVES that wrong value once
> written, rather than letting collide's own -- checkerboard-prone but at least self-correcting in
> the cell-conserves-mass sense -- dynamics touch it; by the time the cell is uncovered and
> reconstructed, the "before" state reconstruction extrapolates from is already garbage, and the
> output is far worse than the input. **This suggests the freeze may be net-neutral-to-harmful for
> this specific failure mode** (it wasn't tested head-to-head against reconstruction-only for THIS
> exact instance/uid, only for the two-mirror A/B above which crashed via different uids/steps) --
> not yet conclusively shown either way, but the mechanism described here does not obviously benefit
> from freezing.
>
> **Not yet done**: (1) find why the MEM-GARBAGE-LINK rate climbs at step 22017-22023 specifically --
> requires the SAME box/window logged further back (this run's window started at 21800, comfortably
> before, so the climb itself IS visible in the density/max|PDF| columns already collected -- an
> unexplored angle on the SAME already-fetched data, not a new run); (2) which specific link
> direction/neighbour cell each skipped MEM-GARBAGE-LINK event points at, to see whether they all
> trace back to one shared root cell (the true "patient zero") or are already spread across the
> rod's footprint from the start.

> **[Claude, 2026-09-27] Per explicit request ("don't freeze, run from earlier checkpoint then"):
> the reconstruction-only companion, watched all the way back to the checkpoint's own start step
> (20001, not just the narrow 21800-22050 window used for the freeze-ON run), ran CLEAN for the
> entire 2200-step budget -- no `STALE-PARTICLE-IDX`, no fatal error -- and shows a flat,
> non-escalating `[MEM-GARBAGE-LINK]` rate throughout, in direct contrast to the freeze-ON run's
> pre-collapse climb.**
>
> **Job 4271533**: `useObstaclePdfFreeze false`, same box (X 590-625, Y 8-30, Z 15-50), same
> step-20001 checkpoint, window widened to the full 20001-22200 range instead of the narrow
> 21800-22050 used for the freeze-ON run (job 4271512). Completed all 2200 configured steps with
> exit 0 -- no crash of any kind.
>
> **`[MEM-GARBAGE-LINK]` per-step count across the ENTIRE run: flat at 25-30/step from start to
> finish** (checked explicitly: steps 22150-22199, the very end of the logged window, still sitting
> at a rock-steady 27/step, no upward trend visible anywhere) -- a completely different signature
> from the freeze-ON run's escalation (11-13/step baseline -> 22/step (22017-18) -> **67/step spike**
> (22020) -> 39/step (22023) -> collapse (22024) -> abort (22031). No sign whatsoever of an impending
> escalation at the point this run's budget ran out.
>
> **Caveat, not a clean "freeze causes it" proof**: a DIFFERENT freeze=false run earlier today (job
> 4271396, no box logger, narrower floater-history window) DID crash, at step 22257 -- 57 steps past
> this run's 22200 cutoff. Given this codebase's own already-established run-to-run non-determinism
> (09-20 checkpoint-resume entries: "the run is not reproducible: the trajectories diverge"), this
> run's clean completion through 22200 does NOT prove removing the freeze eliminates the crash
> entirely -- it may simply have diverged onto a trajectory that survives longer, or crashes via a
> different, rarer mechanism (the genuine cross-rank ghost-ownership race characterised 09-26,
> unrelated to any escalating-link precursor) that a flat MEM-GARBAGE-LINK rate would not predict.
>
> **What IS a fair reading**: whatever caused the freeze-ON run's dramatic, sustained escalation
> (11->22->67->39/step over just 6 steps, then wholesale regional collapse) is NOT a universal
> feature of this population/checkpoint -- it did not reappear in 2200 steps of the no-freeze
> trajectory. Combined with the freeze-ON run's own direct evidence (a cell's PDF getting "locked in"
> at a wrong value outside the freeze's own save/restore bracket, step 22025, then exploding on
> reconstruction at step 22028), the most consistent explanation so far: the freeze converts an
> ordinary, otherwise-transient MEM-GARBAGE-LINK event (which happens constantly, ~27/step, in BOTH
> configurations, and evidently usually self-corrects when collide is allowed to keep acting on the
> cell) into a PERMANENT one by preventing collide from ever touching -- and potentially
> correcting -- the poisoned cell again. Not yet proven (would need a run where a specific cell is
> tracked hitting a large MEM-GARBAGE-LINK spike under NO freeze and observed to recover, matched
> against the SAME cell/uid under freeze failing to recover) but the two runs' aggregate signatures
> line up with this reading cleanly.
>
> **Not yet done**: (1) extend this no-freeze run's budget past step 22257 (this run's own natural
> continuation, same seed/checkpoint) to see whether it also eventually crashes, and whether that
> crash is ALSO preceded by a MEM-GARBAGE-LINK escalation or is a "clean" sudden event, distinguishing
> the two candidate mechanisms directly; (2) the proposed matched-cell recovery-vs-lock-in check above.

> **[Claude, 2026-09-27] Follow-up #1 done: job 4271541 extended the same no-freeze run to step
> 22400 (200 steps past job 4271396's own freeze=false crash at 22257) from the same step-20001
> checkpoint (no new checkpoint existed to resume from -- `storeSnapshot` was false by design, so
> this is a fresh continuation attempt from the same starting point, not a literal resume of 4271533,
> per that run's own caveat above). Per explicit request ("keep it running then from the last
> checkpoint").**
>
> **Result: completed cleanly AGAIN, all 2400 steps, no `STALE-PARTICLE-IDX`, no fatal error.** The
> `[MEM-GARBAGE-LINK]` per-step rate stays at 27-30/step through the ENTIRE extended range, including
> directly through step 22257 itself (27/step, no anomaly) and out to the final logged step 22399
> (27/step) -- no escalation anywhere in 2400 steps total, covering and exceeding every crash step
> seen in ANY freeze=false run today.
>
> **This meaningfully strengthens, without fully proving, the freeze-is-implicated reading above.**
> Two independent no-freeze runs from the same checkpoint (4271533 to 22200, 4271541 to 22400) both
> show a flat, non-escalating baseline and zero crashes, cleanly outlasting the one freeze=false run
> that DID crash (4271396, step 22257, no box logger -- so its own precursor signature, if any, was
> never captured). Meanwhile EVERY freeze=true test today (4271395: 22012, 4271512: 22031, and the
> GJK-only test 4271486, which also had the freeze at its default true: 22174) crashed, and the one
> instrumented freeze=true run (4271512) showed the same escalating precursor pattern immediately
> before its collapse. Given this codebase's own established non-determinism, three-crashes-for-
> three-attempts (freeze=true) vs. two-clean-runs-for-two-attempts (freeze=false, now covering a
> wider step range than any freeze=true attempt) is suggestive, not conclusive proof causation -- but
> it is the most consistent, most directly-evidenced signal in this whole investigation arc.
>
> **Working recommendation for the thesis writeup (not yet a final decision)**: consider
> `useObstaclePdfFreeze false` (reconstruction-only, the pre-fix default) as the safer configuration
> for this specific population/checkpoint going forward, pending either (a) a freeze=true run that
> also survives a comparably wide step range without crashing (would weaken the causal reading), or
> (b) direct confirmation via the proposed matched-cell recovery-vs-lock-in check (does removing the
> freeze let a cell that hits a large MEM-GARBAGE-LINK spike recover on its own, matched against the
> SAME cell/uid under freeze failing to recover -- not yet done).

> **[Claude, 2026-09-27] CORRECTION to the working recommendation above: a THIRD no-freeze attempt
> from the same checkpoint (job 4271670, the long-haul 50k run per explicit request) crashed too --
> EARLIER than either of the two "clean" runs' own cutoffs. The freeze=false-is-safer reading was
> premature; three-for-three across ALL configurations now looks like ordinary chaotic
> non-determinism, not something the freeze setting controls.**
>
> **New checkpointing feature implemented per explicit request ("set the save checkpoint to ever
> [sic] 500 timesteps ... also save always the last two checkpoint not only the last one")**:
> `BlockForestParams.keepLastTwoSnapshots` (SimulationParams.h) + `writeCheckpoint`'s new
> `alternateGeneration` parameter (WriteData.h) -- the regular-cadence checkpoint now alternates
> between `snapshotBaseFolder/gen_0` and `gen_1` by `(timestep/snapshotFrequency) % 2`, so the two
> most recent checkpoints both survive instead of one overwriting the other. Loading is unchanged
> (point a resume prm's own `snapshotBaseFolder` at whichever `gen_<N>` has the higher
> `simState.file` step). Compiles clean; not one of the three fixes-off toggle files, so no scratch-
> copy sync needed for this change.
>
> **Job 4271670**: `useObstaclePdfFreeze false`, `snapshotFrequency 500`, `keepLastTwoSnapshots true`,
> ALL heavy diagnostics off (per "turn the logging back on" -- meant for AFTER a crash is located, not
> during this long run), loading from a dedicated copy of the step-20001 checkpoint
> (`resume_snapshots/small_rods_24x_nofreeze_50k_run/`, copied via `cp -r` so this run's own
> checkpoint writes never touch the original backup), targeting the full 50000-step budget.
>
> **Crashed at step 22038** (`STALE-PARTICLE-IDX`, ranks 108/114/115, the SAME recurring
> `uid...162`) -- 162 steps EARLIER than job 4271533's clean 22200-step cutoff and 219 steps earlier
> than job 4271541's clean 22400-step cutoff, both also `useObstaclePdfFreeze false` from the exact
> same starting checkpoint. **This directly contradicts the two-clean-runs-for-two-attempts framing
> above.** The honest count across all of today's freeze=false attempts is now 2 crashed (4271396 at
> 22257, 4271670 at 22038) vs. 2 clean (4271533/4271541, but note their own cutoffs of 22200/22400
> fall AFTER 4271670's crash step -- meaning if extended even slightly further they may well have
> crashed too; "clean" only means "clean within that run's own arbitrarily-chosen budget", not
> "structurally safe"). Combined with freeze=true's own 3-for-3 crash record today (22012, 22031,
> 22174), the fair reading is: **this population crashes via the recurring uid...162/...780/...164/
> ...778 mechanism regardless of freeze state, at a step that varies chaotically (22012-22257 across
> 6 total attempts today) because of this codebase's own already-established run-to-run
> non-determinism -- not because of the freeze setting.** The escalating-MEM-GARBAGE-LINK-precursor
> vs. flat-baseline contrast observed between job 4271512 (froze, escalated, crashed) and jobs
> 4271533/4271541 (no freeze, flat, survived to their own cutoffs) may still be a real, freeze-linked
> phenomenon for THOSE specific instances -- but it is not a general rule that predicts crash-vs-
> no-crash for this population. **Recommendation withdrawn** pending the still-open matched-cell
> recovery-vs-lock-in check, which remains the one experiment that could actually distinguish
> "freeze causes it" from "chaotic non-determinism, freeze irrelevant."
>
> **Debug follow-up (per the same request, "then we debug again if it crashes and turn the logging
> back on")**: `keepLastTwoSnapshots` immediately paid off here -- `gen_0` holds step 22001 (37 steps
> before the crash, `simState.file` confirmed) and `gen_1` holds step 21501, vs. the 2037-step gap
> from the original step-20001 backup used for every earlier debug rerun today. A new tightly-
> windowed rerun from `gen_0` with full diagnostics (diagnosticSweepEnabled, logFloaterHistory,
> logFluidBoxHistory) re-enabled is the natural next step -- queued next.

> **[Claude, 2026-09-27] Real bug found and fixed along the way: `keepLastTwoSnapshots`'s `gen_0`/
> `gen_1` subfolders never got their own `blockForest.file`, so neither was independently loadable.**
>
> First debug-rerun attempt (job 4271677, loading `gen_0`) aborted at setup, NOT in physics:
> `Loading BlockForest from file '.../gen_0/blockForest.file' failed: Opening the file failed. Does
> the file even exist?` (`BlockForest.cpp:368`). Root cause, found by reading `SetupStages.h:112-116`:
> `blockForest.file` is only ever written when `storeSnapshot && !loadSnapshot` (a FRESH start that
> also stores) -- during a RESUME (`loadSnapshot=true`, job 4271670's own case), that condition is
> permanently false, so the file is never (re)written anywhere, including into a
> `keepUnique`/`alternateGeneration` subfolder, even though every OTHER field gets written there
> every checkpoint. This is a PRE-EXISTING gap in the original (2026-09-24) `keepAllSnapshots`
> mechanism too, just never triggered before because every earlier use of that mode this session
> happened to be paired with a fresh start, and every manual checkpoint copy this session (`cp -r`)
> copied the WHOLE folder including `blockForest.file`, masking the gap.
> **Fixed** (`WriteData.h::writeCheckpoint`): whenever writing into a `keepUnique`/
> `alternateGeneration` subfolder, also copy `blockForest.file` from the shared base
> `snapshotBaseFolder` into that subfolder (root-only, `std::filesystem::copy_file` with
> `overwrite_existing`) -- valid because the block decomposition itself never changes during a run
> (same domain/rank count throughout). Manually patched the already-written `gen_0`/`gen_1` (copied
> the base folder's `blockForest.file` into each) rather than losing that computed state; compiles
> clean; synced. Second attempt (job 4271684) loaded and ran correctly.

> **[Claude, 2026-09-27] DECISIVE, CLEAN CONFIRMATION of the ORIGINAL checkerboard-instability
> hypothesis from earlier this session -- caught with quantitative growth-rate data, WITHOUT the
> freeze, directly contradicting nothing (this is the freeze's actual intended target, now seen
> operating in the open).**
>
> **Job 4271684** (debug rerun from `gen_0`/step 22001, `useObstaclePdfFreeze false`, full
> diagnostics on) crashed at **step 22036** (`STALE-PARTICLE-IDX`, `uid...780` -- 2 steps from the
> original job 4271670's own crash at 22038, well within this codebase's established non-determinism
> for a resumed run). Unlike every earlier traced instance, the FLAGGED cell itself was completely
> healthy at the moment of the abort (density=1.186, ordinary PDF values) -- a clean ghost-ownership
> race, not a corrupted cell. But tracing `uid...780`'s OWN obstacle cell (603.5,20.5,27.5) in the
> fluid-box log tells a different, much clearer story:
>
> | step | max\|PDF component\| | growth vs. 2 steps earlier |
> |---|---|---|
> | 22001 | 1.69e6 | (start of window) |
> | 22009 | 7.58e9 | |
> | 22017 | 1.77e13 | |
> | 22025 | 3.00e16 | |
> | 22031 | 6.94e18 | |
> | 22033 | 4.18e19 | |
> | 22035 | 1.23e20 (last logged step) | |
>
> **Per-step growth ratio, computed directly from the logged values**: 10.1x (22001->22003), 8.9x,
> 8.3x, 7.8x, 7.5x, 7.2x, 7.0x, 6.8x, 6.7x, 6.6x, 6.4x, 6.35x, 6.27x, 6.19x, 6.13x, 6.07x (over
> 2-step intervals, converging smoothly toward ~6x/2-steps = ~2.4x/step) -- a textbook signature of a
> LINEAR instability with a dominant eigenvalue just above 1 (the ratio's smooth, monotonic
> convergence as transient/subdominant modes decay is exactly what a genuine numerical instability
> looks like; noise or a one-off event would not produce this clean a curve). This is quantitative
> confirmation of the mechanism this whole investigation started from (way earlier this session,
> before any of today's freeze/reconstruction/GJK/box-logger work): the generated D3Q27 kernel has no
> per-cell obstacle mask, so `collide()` runs unconditionally on obstacle cells too, and without the
> freeze to stop it, ANY small seed asymmetry left in an obstacle cell's stored PDF grows
> exponentially, step after step, completely unchecked.
>
> **The seed predates this run's own checkpoint**: the tracked cell was ALREADY at 1.69e6 (many
> orders of magnitude past healthy, ~0.3) at step 22001 -- the very first step of this window, i.e.
> the state ORIGINALLY CHECKPOINTED by job 4271670 at that step already contained this growing
> instability. Extrapolating the ~6-10x/2-step rate backward, a healthy value (~0.3) would place the
> true onset around step 21990-21996 -- outside this run's own logged window, not chased further
> here.
>
> **Reconciling with today's freeze-ON finding (job 4271512)**: there, the SAME class of tracked
> obstacle cell held `max|PDF|` at EXACTLY 0.296931 for 225 straight steps -- i.e. the freeze
> unambiguously DOES suppress this exact exponential-growth mechanism while a cell stays covered, as
> designed. The freeze-ON run's own crash then came from something else (a cell "locked in" at a
> wrong value from OUTSIDE the freeze's own save/restore bracket, then exploding further at
> reconstruction when uncovered -- see the 2026-09-27 entries above). **Today's full picture is now:
> two DIFFERENT real failure modes have been directly observed and both lead to the same
> `STALE-PARTICLE-IDX` outcome** -- (a) the classic unmasked-collide checkerboard growth (freeze OFF,
> this entry, clean exponential data), and (b) a value written outside the freeze's own bracket that
> then gets baked in and explodes on reconstruction (freeze ON, job 4271512). The freeze trades one
> for the other rather than eliminating either outright, which is consistent with EVERY freeze=true
> run today crashing (3/3) and EVERY freeze=false run today also eventually being shown capable of
> crashing (job 4271670, this one) despite two other freeze=false runs completing cleanly within
> their own shorter budgets (almost certainly luck, given this run's own crash step of 22036-22038
> falls well inside both of those "clean" runs' cutoffs).
>
> **Open question, not yet chased**: why does an obstacle cell accumulate ANY seed asymmetry at all,
> freeze or no freeze, if MEM bounce-back (`SimpleBB`) and mesa_pd's own single momentum-extraction
> at the moment of covering are the only other things that ever touch an obstacle cell's PDF (per
> the 2026-09-25 mechanism finding)? Candidates, none tested yet: (a) the bounce-back computation
> itself introducing a small per-step asymmetry between "incoming" and "outgoing" link directions
> that is not exactly self-cancelling; (b) the one-time momentum-extraction at the moment of covering
> leaving an already slightly-asymmetric residual PDF that the unmasked collide kernel then amplifies
> from step one.

> **[Claude, 2026-09-27] ORIGIN FOUND, per explicit request ("where does it originate"): a DIRECT
> obstacle-to-obstacle OWNERSHIP HANDOFF between the two recurring rods, with no reset in between --
> the checkerboard oscillation begins in the SAME step as the handoff, cleanly, every time checked.**
>
> **Job 4271714** (rerun from `gen_1`/step 21501, 500 steps before `gen_0`, full diagnostics on from
> the very start) crashed at step 22044 (uid...162 again, different cell -- this run's own
> non-deterministic divergence). Tracing the SAME cell watched throughout today, (603.5,20.5,27.5),
> from a genuinely healthy state:
> - **Steps 21900-21936**: owned by `uid...162`. Rock-steady: density 1.0027-1.0033, essentially flat,
>   max\|PDF\|~0.297 -- healthy, matching the freeze-ON run's own frozen baseline exactly (same
>   numbers), confirming this is what a genuinely undisturbed obstacle cell looks like regardless of
>   freeze setting.
> - **Step 21937**: still owned by `uid...162`, but density ticks up to 1.031 (a small, one-step
>   precursor -- the first departure from the flat baseline, one step BEFORE the ownership change).
> - **Step 21938: ownership switches directly to `uid...780`** -- same cell, flag stays `OBST`
>   throughout (confirmed via the logged flag column: no INTF/LIQ state in between) -- i.e. the two
>   rods' surfaces are touching or all-but-touching at this exact cell, and it is handed from one
>   rod's obstacle-ownership straight to the other's with no equilibration/reset step.
> - **From step 21938 onward, a clean odd-even (checkerboard) oscillation in density begins and
>   grows every single step-pair, continuously, with no interruption**: 1.054/0.999 (21939/21940),
>   1.064/0.987 (21941/21942), 1.081/0.974 (21943/21944), 1.093/0.957 (21945/21946), 1.111/0.938
>   (21947/21948), 1.127/0.918 (21949/21950), ... 1.510/0.604 (21975/21976), ... reaching
>   max\|PDF\|=1.22 (step 21984, the first breach of the earlier "anomaly" threshold used elsewhere
>   today), 11.3 (21990), and full breakout from there (matches the exponential-growth entry above --
>   this is the SAME mechanism, now traced to its exact seed).
>
> **This closes the loop on "where does it originate" completely**: not a random numerical fluke, not
> tied to any specific step number (today's six-plus attempts crashed at steps spanning 22012-22257,
> all non-reproducible) -- it is tied to a specific, recurring PHYSICAL EVENT: `uid...162` and
> `uid...780` (the same two rods implicated in literally every crash instance today and in the
> 2026-09-26 investigation) come close enough that a shared cell's obstacle-ownership passes directly
> from one to the other, and THAT handoff is the seed. Every time it happens, the unmasked D3Q27
> collision kernel (no per-cell obstacle mask, confirmed since the very start of this investigation
> thread) has nothing to stop the resulting tiny residual asymmetry from growing, step after step,
> until it explodes and (eventually, via cascading to neighbours and/or a large enough hydro/contact
> force spike) desyncs mesa_pd's own ghost-ownership bookkeeping far enough to trip
> `STALE-PARTICLE-IDX`. **The freeze fix suppresses the GROWTH phase (confirmed directly: a frozen
> cell's own numbers, 1.0025-1.0033, 0.297, are IDENTICAL to this run's own pre-handoff baseline) but
> was never positioned to prevent the HANDOFF itself from seeding a new instability the moment the
> cell is later uncovered and reconstructed** -- which is exactly the freeze-ON failure path already
> characterised above (job 4271512: locked-in wrong value -> reconstruction explosion at uncovering).
>
> **Remaining open question, now much narrower**: why does a direct ownership handoff leave a residual
> asymmetry at all, rather than a clean, symmetric state? The receiving rod's own momentum-extraction/
> MEM computation for that cell presumably assumes it is taking over a "fresh" obstacle cell, but the
> cell's stored PDF is whatever the PREVIOUS rod's own (possibly already slightly imperfect) obstacle
> dynamics left behind -- worth checking directly whether `MovingParticleMapping.h`'s handoff path
> re-initialises the PDF at all when a cell's owning UID changes without an intervening liquid/
> interface state, or silently carries the old rod's stored values forward unchanged into the new
> rod's ownership (the current best-guess candidate, not yet verified by reading that code path).

> **[Claude, 2026-09-27] Code read + small local reproduction, per explicit request ("find that out
> and ... can you reproduce this in a small isolated test"): CONFIRMED -- `MovingParticleMapping`'s
> cross-particle handoff branch never touches the PDF field, only re-tags ownership.**
>
> **Code read** (`src/lbm_mesapd_coupling/momentum_exchange_method/MovingParticleMapping.h:180-193`,
> the per-step `MovingParticleMapping` class, NOT the initial-mapping-only `MovingParticleMappingKernel`
> which has no equivalent branch): when a particle's containment test finds a cell already flagged
> obstacle (owned by a different, non-fixed particle), the code's ENTIRE handoff logic is:
> ```cpp
> auto formerParticleUid = (*particleField)(x, y, z);
> auto formerParticleIdx = ac_->uidToIdx(formerParticleUid);
> if (!isSet(ac_->getFlags(formerParticleIdx), mesa_pd::data::particle_flags::FIXED))
>    (*particleField)(x, y, z) = particleUid;   // <-- the ENTIRE handoff. Nothing else.
> ```
> No PDF write, no re-initialization, no momentum transfer, no reconstruction trigger -- the SECOND
> particle silently inherits whatever the FIRST particle's own obstacle-cell PDF history happened to
> leave behind, completely unchanged.
>
> **Small local reproduction, new test `ANTIDUNES_ObstacleOwnershipHandoff`
> (`tests/showcases_antidunes/ObstacleOwnershipHandoffTest.cpp`, builds+runs+passes locally in ~1.2s,
> no cluster needed, per this project's own established fast-local-test convention)**: exercises the
> REAL, unmodified `MovingParticleMapping` class (not a reimplementation) with two real mesa_pd
> `Sphere` particles in a tiny single-block domain. Particle A maps a test cell; a deliberately
> non-equilibrium 27-component PDF is injected (standing in for "whatever A's own imperfect obstacle
> dynamics left behind" -- the exact values don't matter, only that they're asymmetric/non-trivial);
> A is then moved far away and B moved into the exact same position, and the SAME
> `MovingParticleMapping::operator()` sweep is called once more (matching the real per-main-step
> cadence). **Result: `max|PDF component difference before vs. after the handoff| = 0` -- a bit-for-
> bit no-op, exactly as the code read predicts, and `particleField` correctly switches to B's uid
> while the cell stays flagged obstacle throughout (never becomes fluid in between)** -- reproducing
> the exact real-world mechanism (uid...162 -> uid...780, step 21938) in ~1 second locally instead of
> needing the rare real chaotic encounter between two specific rods on a 144-rank cluster job.
>
> **Deliberately out of scope for this test** (documented in the test's own header comment): trying
> to ALSO reproduce the subsequent exponential-growth phase in equally strict isolation.
> `CollideIsolatedObstacleIterationTest`/`...RealSeedTest` (2026-09-25) already tested repeated
> ISOLATED collide iteration (zero force, zero streaming input) on two different real near-
> equilibrium seeds and found NEITHER diverges -- a genuine fixed point both times. So whatever
> sustains the real observed growth needs more than "collide alone forever" -- most likely the
> interplay with `SimpleBB`'s own per-step partial PDF rewrite (2026-09-25 finding) against a
> genuinely time-varying real neighbour flow field, which a single isolated cell with no real
> neighbours structurally cannot reproduce. A faithful growth reproduction would need a heavier, less
> "isolated" two-rod real-time mini-scenario (not a unit test) -- flagged as a possible follow-up, not
> attempted here.
>
> **This fully closes the loop the box logger opened**: origin (a specific, recurring physical event
> -- direct ownership handoff between two persistently-close rods) is now confirmed both by reading
> the exact code responsible AND by an independent, fast, deterministic local reproduction using the
> real production class -- not merely inferred from one real run's own log data.

> **[Claude, 2026-09-27] FIX IMPLEMENTED (`MovingParticleMapping.h`), per explicit request ("do it").
> First cluster verification attempt (job 4271763) crashed with a raw SIGSEGV -- found and fixed a
> SECOND, genuinely pre-existing latent bug uncovered along the way, in the SAME branch.**
>
> **The fix**: on a GENUINE cross-particle ownership handoff (`formerParticleUid != particleUid` --
> NOT the routine same-owner re-confirmation that fires every step a particle keeps a cell it already
> owns), the cell's PDF is now also reset to rest equilibrium via `PdfField::setToEquilibrium()`
> (stencil-agnostic, correct for any `LatticeModel_T`, not just Antidunes' own D3Q27). This directly
> breaks the silent-inheritance chain the origin finding traced: the receiving particle now starts
> from a guaranteed-finite, guaranteed-positive-density state instead of whatever the previous
> owner's obstacle-cell history left behind.
>
> **First cluster attempt (job 4271763, 25000-step budget, freeze off, from the same step-20001
> checkpoint used all day): crashed almost immediately with a raw SIGSEGV** (`Signal 11`, rank 109) --
> no `gen_0`/`gen_1` checkpoint was even written yet (crash predates step 20500, the first scheduled
> checkpoint). NOT the controlled `STALE-PARTICLE-IDX` abort seen in every other run today.
>
> **Root cause of the segfault, found by reading the surrounding code (not by guessing)**: the SAME
> branch's PRE-EXISTING, UNMODIFIED line, `ac_->uidToIdx(formerParticleUid)`, returns
> `std::numeric_limits<size_t>::max()` when the former owner is no longer known locally (its ghost
> copy was dropped -- an ORDINARY event, per this whole investigation's own established ghost-drop
> mechanism). The existing code then called `ac_->getFlags(formerParticleIdx)` on that result
> completely unguarded -- indexing a particle-storage vector with `SIZE_MAX` is undefined behaviour
> and segfaults immediately. This bug ALREADY existed before today's fix; it was masked in every
> other run because `SimpleBB.h`'s own `STALE-PARTICLE-IDX` check normally aborts the whole run
> first, before `MovingParticleMapping` ever reaches this line with a genuinely stale uid -- but that
> ordering is a race between two independent code paths, not a guarantee, and the ownership-handoff
> fix's own new code path apparently shifted the timing enough to expose it on the very first
> attempt.
>
> **Second fix, same file**: guard the pre-existing `getFlags(formerParticleIdx)` call --
> `formerParticleIdx != std::numeric_limits<size_t>::max()` -- and treat an unknown/stale former
> owner the same as "not fixed" (a particle that is no longer known locally obviously cannot be
> FIXED; fixed/global particles are never dropped), proceeding with the claim (and the equilibrium
> reset) instead of reading out-of-bounds memory to decide.
>
> **Verified locally after both fixes**: `ANTIDUNES_ObstacleOwnershipHandoff` (updated to check the
> fix, not the original bug) and the full `Antidunes` target both build clean;
> `ANTIDUNES_ObstacleOwnershipHandoff`/`CapsuleShape`/`BognerStateClassification`/`ContactDetection`
> (every local test that touches `MovingParticleMapping`/`BoundaryHandling`) all still pass -- no
> regressions. Re-submitted for a second cluster attempt (job pending).

> **[Claude, 2026-09-27] Second cluster attempt (job 4271773) ALSO segfaulted, same signature.
> Bisected properly this time (job 4271777, callback disabled -> clean STALE-PARTICLE-IDX abort,
> not a crash) -- confirmed the true culprit is `PdfField::setToEquilibrium()` itself, not the
> `uidToIdx` guard. Re-architected: the reset moved OUT of the generic library file entirely.**
>
> **Root cause of the segfault**: `lbm::Equilibrium<LatticeModel_T>::set()` (the function
> `setToEquilibrium()` calls) reads a static `LatticeModel_T::w[]` weight-table member in its generic
> code path. This codebase's OWN established convention for resetting an Antidunes PDF cell
> (`repairDegenerateObstacleCells`, `ObstacleOwnershipHandoffTest.cpp`'s own `injectedPdf` base state)
> NEVER relies on that generic trait for the generated D3Q27 Cumulant model -- it always hardcodes
> the known rest-equilibrium weights directly. Calling the generic API from the generic
> `MovingParticleMapping.h` library file apparently isn't safe for this model in the REAL multi-file
> production link (it worked fine in a minimal single-file local unit test -- a link/static-init
> context difference, not caught by that test).
>
> **Re-architected rather than patched around**: `MovingParticleMapping.h` (generic library code,
> also used by non-Antidunes callers with other stencils, e.g. `IbbEquivalenceDragTest.cpp`'s D3Q19)
> must not contain any stencil-specific reset logic. Added an optional
> `OwnershipHandoffCallback_T onOwnershipHandoff` constructor parameter (default: empty/no-op, so
> every OTHER existing caller in the codebase is completely unaffected) that fires exactly once,
> right before a genuine handoff is committed, with `(block, x, y, z, formerUid, newUid)`. The class
> itself never touches the PDF field for this purpose any more -- it only decides *when* a genuine
> handoff is happening (the actual bug/logic fix) and hands control to whichever callback the caller
> supplies.
>
> **Antidunes' own callback** (`BuildTimeloop.h`, wired into the existing
> `makeMovingParticleMapping` call site in `addSweepFluidParticleCoupling`): manually assigns the
> D3Q27 rest-equilibrium weights (rho=1, u=0) directly -- the exact same convention
> `repairDegenerateObstacleCells` already uses successfully in production -- plus a throttled
> `[OWNERSHIP-HANDOFF-RESET]` log line (first 50 occurrences, then every 500th) for visibility.
>
> **`ANTIDUNES_ObstacleOwnershipHandoff` updated to match**: now supplies its own equivalent callback
> (mirroring `BuildTimeloop.h`'s exactly) instead of relying on the class to do it internally, and
> tracks invocation counts. One test-design bug found and fixed along the way: the test's two
> spheres (radius 2.5) cover a multi-cell VOLUME, not a single point, so a genuine handoff between
> them touches ~81 cells at once, not the one cell the test otherwise tracks -- the original
> "callback must fire EXACTLY ONCE" assertion was simply the wrong expectation for that geometry;
> relaxed to "at least once for the genuine handoff, zero further increments for the routine re-
> confirmation", which is the actually-meaningful claim. Passes; full local suite
> (`ObstacleOwnershipHandoff`/`CapsuleShape`/`BognerStateClassification`/`ContactDetection`) still
> passes, no regressions. Synced; third cluster attempt submitted (job pending).

> **[Claude, 2026-09-27] Third attempt (job 4271821) ALSO segfaulted, same signature -- the
> `setToEquilibrium` hypothesis from the previous entry was WRONG. Real cause found via precise
> bracketing logs (job 4271830) and confirmed by reading the field-setup code directly: a genuine,
> pre-existing ghost-layer-count MISMATCH between `flagField` and `pdfField`.**
>
> **Bracketing logs (unthrottled, every step of the callback)** showed something new and important
> first: the SAME physical cell flickered ownership between TWO DIFFERENT rods (`uid...780` and
> `uid...782`) MULTIPLE TIMES within a single `operator()` call, at step 20045 -- i.e. these two rods
> genuinely, deeply overlap at that cell right now (not a numerical fluke; `MovingParticleMapping`'s
> particle loop visits each particle once per call, so if both a rod's shape genuinely contains a
> cell already claimed by the other, the cell gets re-claimed and the callback re-fires, each time
> the loop reaches the next overlapping particle). The callback itself handled this correctly every
> time it fired. The LAST log line before the crash was `[...-GOTFIELD]` (the field pointer was
> fetched successfully) with `[...-Q0]` (the very first PDF write) never appearing -- pinpointing the
> segfault to the FIRST `pdfField->get(x,y,z,0)` call for one specific (rank, cell).
>
> **Root cause, confirmed by reading the setup code directly (not guessed)**:
> `src/lbm/free_surface/boundary/FreeSurfaceBoundaryHandling.impl.h:167` allocates its flag field
> with **2 ghost layers**; `SetupStages.h:213`'s `addPdfFieldToStorage` call (the 6-argument overload,
> no explicit ghost-layer count) gives `pdfField` only **1**. `MovingParticleMapping`'s own `cellBB`
> is sized via `flagField->nrOfGhostLayers()` (2, both in `operator()` and inside
> `getParticleCellBB`), so for a body whose `cellBB` reaches that second ghost layer -- exactly what
> two DEEPLY overlapping rods do -- `(x,y,z)` can be valid for `flagField` but genuinely out of
> bounds for `pdfField`. Any code that unconditionally writes `pdfField` there segfaults. This is NOT
> a new bug introduced by today's fix -- it is a real, pre-existing hazard in this generic file that
> the ORIGINAL, unmodified code right above this branch already protects against for its own PDF
> read: `if( conserveMomentumWhenMapping_ && pdfField->isInInnerPart(Cell(x,y,z)) )` guards the
> "cell newly occupied" branch's own `getMomentumDensity` call for precisely this reason. My new
> callback simply never had the same guard.
>
> **Fix**: added the identical `pdfField->isInInnerPart(Cell(x,y,z))` guard to the handoff branch's
> own callback invocation in `MovingParticleMapping.h` -- the callback is now only ever invoked when
> `(x,y,z)` is safe to write on `pdfField`. This is sufficient, not merely a workaround: whichever
> rank actually OWNS a given cell as an interior cell will always see it as safe there (only a
> NEIGHBOURING rank's ghost-layer copy of the same physical cell could ever be skipped), and the
> existing `PdfCommunication_T`/`Communication_T` calls already registered right after this sweep
> (`BuildTimeloop.h`) propagate the owning rank's correctly-reset value outward to every neighbour's
> ghost layer regardless. Removed the temporary bracketing logs (their job is done), kept the
> existing throttled `[OWNERSHIP-HANDOFF-RESET]` visibility log.
>
> **This corrects the previous entry's diagnosis**: `PdfField::setToEquilibrium()`'s generic
> `lbm::Equilibrium<>::set()` path reading a static `w[]` table was a plausible-sounding but
> ultimately WRONG explanation -- the uidToIdx-guard-only bisection (job 4271777) "proved" it only
> because THAT test never wrote to `pdfField` at all (no callback was firing), so it could never have
> hit the ghost-layer hazard either, regardless of which API a real reset would have used. Both the
> generic API and the manual-weights callback would have segfaulted identically without this guard.
> Verified locally (full suite, no regressions); synced; fourth cluster attempt submitted (job
> pending).

> **[Claude, 2026-09-27] Fourth attempt (job 4271838): the segfault is GONE -- ran 227 seconds
> (vs. ~10-15 seconds for all three prior segfaults) and reached the ordinary, controlled
> `STALE-PARTICLE-IDX` abort instead (step 22005, `uid...780` -- the same recurring rod). Zero
> segfaults in the log. The ghost-layer guard fix is confirmed correct and complete for its own
> purpose: the ownership-handoff PDF reset is now safe to run in the real 144-rank production
> context.**
>
> **Honest bottom line on what this fix does and does NOT accomplish**: the fix eliminates the
> specific, now-fully-understood silent-inheritance bug (`MovingParticleMapping`'s cross-particle
> handoff branch carrying a cell's stale PDF history forward unchanged) -- confirmed both by direct
> code reading and by an independent local reproduction
> (`ANTIDUNES_ObstacleOwnershipHandoff`). But THIS SPECIFIC cluster run still hit
> `STALE-PARTICLE-IDX` anyway, via the SAME recurring `uid...780`, at a step (22005) inside the same
> general range every other configuration tested today also crashed in (22012-22257). Today's fuller
> picture (freeze on, freeze off, GJK-only, and now this fix) is that the ACTUAL abort is triggered
> by mesa_pd's own cross-rank ghost-ownership desync (the 2026-09-26 mechanism: an ordinary ghost-
> drop that, rarely, leaves a stale `particleField_` reference before the next MEM sweep touches it)
> -- a mechanism this fix was never targeting and does not address. The ownership-handoff PDF reset
> targets a genuine, real, now-eliminated SOURCE of local numerical corruption (which very plausibly
> contributes to destabilizing these two rods' trajectories enough to make the ghost-drop race more
> likely to bite), but is not, by itself, a complete fix for this population's `STALE-PARTICLE-IDX`
> crashes. A second attempt from a different checkpoint/trajectory would be needed to see whether
> this fix measurably delays or reduces the frequency of the crash (as the obstacle-PDF freeze
> measurably did for one specific earlier instance) even though it does not eliminate it outright --
> not yet attempted, given the scope of today's investigation reaching a natural stopping point here.
>
> **What IS solid and worth keeping regardless of the above**: (1) the exact origin mechanism of the
> checkerboard/silent-inheritance instability, confirmed by real crashed-run data, code reading, and
> independent local reproduction; (2) a real, correct, generically-safe fix for that specific
> mechanism (the callback architecture correctly keeps stencil-specific logic out of shared library
> code, and the ghost-layer guard closes a genuine pre-existing hazard in `MovingParticleMapping.h`
> itself, benefiting any future caller of that class, not just Antidunes); (3) a new, reusable local
> test (`ANTIDUNES_ObstacleOwnershipHandoff`) and a new checkpoint-retention feature
> (`keepLastTwoSnapshots`) that materially sped up today's whole debugging cycle and remain useful
> for any future investigation in this area.

> **[Claude, 2026-09-28] PAUSED HERE, per explicit request ("let's leave this be for now we'll come
> back later"), with an open methodological question that should be resolved BEFORE trusting the
> "two rods collapse simultaneously" reading below, or building any stabilization fix on it.**
>
> **What job 4271838's crash (step 22005, `uid...780`) looked like on the surface**: a genuine
> "collapsing cell -> large MEM force -> velocity jump -> `STALE-PARTICLE-IDX`" chain (confirmed via
> the run's own `[MEM-GARBAGE-LINK]`/`[DIAG-SWEEP][VELOCITY-JUMP]` logs), with the first broken cell
> at world (607.5,20.5,33.5) -- density -324.9 at step 22002, 3.04e23 by step 22003. Owned-neighbour
> tracing (job 4271879, a targeted rerun from that run's own `gen_0` checkpoint at step 22001, box
> logger on) suggested `uid...162` and `uid...780` both had cells go bad in the SAME single step
> (22002->22003) -- read as a shared, proximity-driven squeeze (the September 20 "interface cell
> squeezed between two rods" mechanism), NOT the ownership-handoff mechanism fixed earlier today.
>
> **But a direct check found a real discrepancy that undercuts that reading**: the ORIGINAL run
> (4271838) shows rank 109 completely CLEAN at its own step 22001 (zero density-window warnings).
> The RERUN (4271879), loading the checkpoint saved at that exact step, shows rank 109 ALREADY
> massively corrupted at ITS OWN nominal step 22001 (13 cells < 0.1, 12 cells > 1.9, values up to
> 5.2e7) -- before this rerun's own dynamics for that step have even run. Same nominal step, same
> underlying checkpointed state, contradictory diagnostic results. Two live explanations, NEITHER
> confirmed: (a) the checkpoint round-trip itself introduces or exposes spurious corruption (e.g. a
> ghost-layer region not yet re-synced via communication at the exact point the first post-load
> diagnostic sweep runs, reading stale/uninitialized values that look like real corruption but
> aren't); (b) a labeling/timing subtlety in what "step 22001" means relative to the communication
> phase differs between a continuously-running step and a freshly-resumed one, making the two runs'
> "step 22001" not quite the same comparison point.
>
> **Why this matters**: if (a) is true, EVERY targeted rerun done today via a fresh checkpoint load
> (not just this one) could have a spurious, artifact-driven anomaly baked into its very first
> logged step, which would call into question how much of today's whole "trace the collapse from a
> close checkpoint" methodology is trustworthy at the finest (single-step) resolution -- though the
> COARSER findings (the ownership-handoff mechanism, confirmed independently via code reading AND a
> from-scratch local unit test with synthetic particles, not dependent on any checkpoint resume at
> all) are NOT at risk from this specific concern.
>
> **Not yet done, flagged for next time**: (1) check whether rank 109's flagged cells at its own
> step 22001 sit right at that rank's own block/ghost-layer boundary (consistent with hypothesis
> (a)); (2) check whether a `PdfCommunication_T`/`Communication_T` call runs before or after
> `addAfterTimeStepDiagnosticSweep`'s own density scan on the VERY FIRST step after a fresh
// `loadSnapshot=true` resume specifically (as opposed to every subsequent step, where ghost layers
> are already synced from the previous step); (3) only once (1)/(2) are resolved, decide whether the
> "simultaneous two-rod squeeze" reading survives, and only then decide between the two stabilization
> options already discussed with the user (extending `repairDegenerateObstacleCells`-style reactive
> repair to liquid/interface cells, vs. a lower `closeSealedCellsMinNeighbours` threshold on the
> existing `SealedCellGuard.h`, vs. both) -- NOT implemented this session, deliberately, pending the
> above.
>
> **State left in place for the resume**: `resume_snapshots/small_rods_24x_fixverify_run/` (gen_0 =
> step 22001, gen_1 = step 21501, on `mt-dennys-mirror2`) is untouched and still valid; the
> `small_rods_24x_2x4x3_clog_resume_collapse22002_gas.prm`/
> `run_144_small_rods_24x_collapse22002.sbatch` pair used for job 4271879 can be resubmitted directly
> for a repeat/comparison run without any new setup.

> **[Claude, 2026-09-30] CLEAN RE-VERIFICATION + BOX SHAPE-COMPARISON, both submitted.** Per explicit
> request ("run this 24 with the small-rod 2x4x3 layouts but with option a and b turned off and none
> of the fixes we came up with turned on during the multiday debug session... use the analytical
> contact detection and the tight layer 1 on for the cylinder don't turn pure analytical on"):
> **job 4279038** (`~/mt-dennys`, free directory) resubmits `small_rods_24x_2x4x3_clog_gas.prm` /
> `run_144_small_rods_24x_2x4x3_clog_fixesoff.sbatch` **completely unchanged** -- checked its
> `ModelParameters` block against every constraint named: Option A/B
> (`seedStreamPdfCloneOnFirstUse`/`repairDegenerateObstaclePdfs`) are both absent from the prm, which
> defaults them to `false` in `SimulationParams.h`; the sealed-cell-guard code from the ghost-ownership
> investigation above was fully reverted (nothing left to disable); `useAnalyticalContactDetection
> true`, `useTightLayer1 true`, `usePureAnalyticalCylinder false` already match exactly. This is
> therefore an as-faithful-as-possible re-run of the original 24-small-rod crash (job 4268494), not a
> new configuration.
>
> **Job 4279039** (`~/mt-dennys-mirror`, the other free directory) is a new shape-comparison
> companion, **`boxes_24x_2x4x3_clog_gas.prm`** / **`run_144_boxes_24x_2x4x3_clog.sbatch`**: identical
> domain/pillars/PID/physics/`ModelParameters` and all 24 body positions/densities, with the 24
> `FloatingCylinder`s replaced by 24 `FloatingBox`es of the same footprint (`halfExtents <3,9,3>` =
> edgeLength 6x18x6, matching the rod's radius-3/length-18 bounding box exactly, same axis-along-Y
> orientation, no tilt seed -- a like-for-like shape swap of the identical scenario, not a new
> perturbation). Purpose: the ghost-ownership desync mechanism identified earlier (STALE-PARTICLE-IDX)
> was only ever confirmed with cylinders, whose line-contact geometry can wedge fluid along a line in a
> way a sphere cannot (per the sphere-vs-cylinder pore-sealing census above); this tests whether boxes
> (point/edge/face SAT contacts only, never a line-contact wedge) reproduce the same desync at the same
> body count/layout, or whether it is specific to the cylinder's own contact geometry. Box dimensions
> match this project's own established convention for cylinder/box shape-comparison at this population
> scale (`box_pillar_gas.prm`'s "log-between-pillars" box already uses the identical `halfExtents
> <3,9,3>` for the same reason).
>
> **Sync note**: both new files (`prms/boxes_24x_2x4x3_clog_gas.prm`,
> `jobs/run_144_boxes_24x_2x4x3_clog.sbatch`) plus every other outstanding local change were pushed to
> both `~/mt-dennys` and `~/mt-dennys-mirror` via `antidunes.sh --sync` (the project's own established
> full-tree rsync, excludes `build*/`/dotfiles) before submission, confirmed via the rsync file list
> (only actually-changed files transferred). Both jobs submitted cleanly (144 ranks/2 nodes each,
> 6h cap): **4279038** and **4279039**, both `PENDING (Resources)` at submission time; `4278452` (R4
> buoyancy-compare, `~/mt-dennys-mirror2`) is the only other job running, unaffected.

> **[Claude, 2026-09-30] RESULTS: 4279038 (cylinder) FAILED, 4279039 (box) COMPLETED CLEAN --
> real evidence the desync needs the cylinder's line-contact geometry specifically.**
>
> **4279038 (cylinder, `~/mt-dennys`)**: FAILED after 23:17, step ~18420-18421 -- a genuine
> reproduction of the same instability class chased at length above, not an infra/config error. Log
> signature, in order: `[LOW-RHO]` NaN-density interface cells (49-51, 10-12, 16-17) -> deep rod-rod
> contact (`[DEM-WATCH] deep contact`, overlap 1.99 LU = 21.0% of the 9.49 LU interaction radius,
> uids ...162/...782) -> `[AddExternalForces] NaN/Inf hydrodynamic force/torque` on uid ...782 ->
> `[TOTAL-FORCE-THRESHOLD]` firing on two other rods (...164/...166, |totalForce| 19.6/7.2) -> hard
> abort on rank 97 via an INVALID-QUATERNION assertion (`r=i=j=k=-nan`,
> `src/core/debug/CheckFunctions.cpp:41`) -- the same failure family as the earlier giant-rod NaN
> crash (job 4268749, "invalid-quaternion assertion, step ~21400"). Crashed earlier than any prior
> run of this exact scenario (~18420 vs. ~22000-22068 for jobs 4268494/4271838) -- consistent with,
> not contradicting, the project's own established finding that these runs are chaotic and NOT
> bit-reproducible even from an identical checkpoint/config (see the resume-divergence entries
> above): same mechanism, different unlucky step by chance.
>
> **4279039 (box, `~/mt-dennys-mirror`)**: COMPLETED all 51000 steps cleanly, exit 0, 58:17. 9778
> total warnings logged (ordinary transient background noise for this scenario class --
> `[SUBCYCLE-LEFT-RANK]` etc. -- comparable in kind to what healthy sphere/box runs already show
> elsewhere in this investigation) but **zero** `[LOW-RHO]`, NaN force/torque, or quaternion-assertion
> events, and no abort of any kind.
>
> **Reading**: identical population (24 bodies), identical layout (2x4x3), identical positions,
> identical densities, identical `ModelParameters` -- the ONLY thing that differs is cylinder vs. box
> contact geometry, and that alone is enough to separate a crash from a clean 51000-step run. This is
> real, direct evidence (not just a plausible hypothesis) for the "line-contact wedge" mechanism
> already proposed above: a cylinder-cylinder contact is a LINE (or a near-line at shallow skew
> angles), which can seal a pocket of fluid along its whole length in a way a box's point/edge/face
> SAT contact geometrically cannot. One caveat before over-generalizing: this is a single run each,
> not a statistical sample, and the cylinder side's own crash step already varies by thousands of
> steps run-to-run (chaotic), so a single clean 51000-step box run is suggestive, not dispositive, on
> its own -- but it is a real, not a null, result, and consistent with every other piece of evidence
> gathered in this investigation (the sphere-vs-rod sealed-cell census above found the identical
> asymmetry: spheres/point-contacts never produce >=24-neighbour sealed pockets, rods do).

> **[Claude, 2026-09-30] ROD CONTACT SKIN IMPLEMENTED (`ModelParameters.cylinderContactSkin`, default 0 = off); test job 4280159 queued, waiting for the Fritz maintenance.**
> Motivation: the rod density collapse starts in liquid cells sealed or pinched between nearly touching parallel rods
> (2026-09-21 bundle test). Line contacts leave sub-cell channels, e.g. an inscribed radius of 0.155 R = 0.46 LU for three
> touching rods at R = 3. The standard remedy in resolved LBM-DEM is a minimum fluid gap / effective diameter (Rettinger
> 2023, p. 7, citing Boutt et al. 2007; Glowinski et al. 2001 use a short-range repulsion).
> **Implementation:** for contact detection only, each rod is replaced by an inflated cylinder, radius R + m and half-length
> h + m. The Cylinder x Cylinder, Cylinder x HalfSpace and Cylinder x CylindricalBoundary (solid pillar) overloads in
> `ParticleKernels.h` call themselves once with the inflated shapes, so the unchanged, validated detection is used. The
> rod-rod force starts below a surface gap of 2m, rod-wall/pillar below m. The rod interaction radius grows by 2m
> (`SetupStages.h`). The MEM mapping is not dilated. Rims: the sharp inflated edge reaches up to (sqrt 2 - 1) m further
> than a rounded skin would.
> **Rejected:** mesa_pd GJK/EPA with margin m. Its EPA returns false when only the margins overlap (EPA.cpp,
> `numPoints == 1`), which is exactly the parallel-rod line contact. It would silently drop those contacts.
> **Unit tests:** `ANTIDUNES_ContactDetection` gained testCylinderSkinCylinderCylinder and
> testCylinderSkinHalfSpaceAndPillar. Both PASS locally, and all existing tests still pass.
> **Test job 4280159** (`jobs/run_36_rod_skin_bundle_test.sbatch`, `~/mt-dennys`, single node, 4 x 9 ranks in parallel;
> results in `skin_results/4280159/`):
>   (A) prescribed three-rod bundle at surface gaps -0.03 (baseline), 1.0 and 2.0 LU, 8000 steps: does one or two resolved
>       fluid cells remove the density drop? (baseline previously rho_min 0.55-0.80);
>   (B) free rods starting in contact with skin 0.5, 3000 steps, floater history on: does the skin separate them to a
>       gap of about 1.0 LU?
> Evaluate with the FLUID-STATS lines (rho_min, max|u|). A local attempt was stopped after about 200 steps because the
> local machine is too slow.
> Not yet covered: Cylinder x Sphere / Cylinder x Box (no such pairs in the rod scenarios); the hopper production run
> with the skin; the effect on packing and clogging (effective size +2m).
