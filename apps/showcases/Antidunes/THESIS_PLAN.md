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

---

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

---

## 1c. [ ] Showcase: mixed-density floaters, some float / some sink

**Status**: we already have runs like this (mixed-shape/mixed-density sedimentation scenarios) — just
need to pick the best one as the actual thesis showcase (or run a clean dedicated one if none of the
existing runs frames the float-vs-sink contrast clearly enough).

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

### Challenges & Future Improvements
- No existing evaluator computes upstream elevation or local pressure at an obstruction directly —
  this likely needs a new function in `Evaluators.h` rather than reusing an existing one verbatim.
- "Clogging" itself needs an operational definition for the thesis (e.g. a sustained elevation rise
  above some threshold, or a sustained drop in downstream transport rate) — deciding this before
  running the sweep will make the results much easier to present than defining it post-hoc from
  whatever the runs happen to show.
- Needs a documented no-clog reference run (e.g. same debris density but spread out / open channel)
  as a baseline, so a measured elevation/pressure rise has something concrete to be compared against.

---

## 4. [ ] Clogging scenarios

Design and run scenarios that specifically target clogging (debris jamming/blocking flow), starting
simple and building up:
1. Sphere-only clogging first (simplest shape, cleanest to interpret).
2. Then progressively more complex scenarios mixing in cuboids/cylinders once the sphere-only case
   is understood.

Depends on item 3 (measurement setup) being in place first.

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
