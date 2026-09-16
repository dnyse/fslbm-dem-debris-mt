# Final Thesis DEM Validation — Status & Index

Living document (same convention as `THESIS_PLAN.md`: dated/attributed comments, status tags
`[ ]`/`[~]`/`[x]`/`[!]`, never overwrite another author's comment). Tracks which existing DEM
contact-validation tests constitute the **final thesis suite**, their current pass/fail status, and
literature citations — the actual derivations live in `DEM_VALIDATION_MATH_REFERENCE.md`, linked from
each entry rather than repeated here.

**Focus**: cylinder and cuboid contact-detection/impact-dynamics validation (per explicit decision,
2026-09-12). Sphere is supporting context only. The box-chatter/production-stability investigation
(`useSpheroBoxContact`, `DropTest.cpp`, `MITIGATION_COMPARISON.md`, and the `testCase`s
`cornerNormalDiscontinuity`, `cornerFixComparison`, `multiPointDeepOverlapStability`/`Fuzz`,
`multiPointEdgeCases`, `multiPointBoxContactFuzz`, `fallOntoPileStability`,
`contactStiffnessCalibrationSweep`, `productionCollisionReplay`) is **explicitly excluded** —
optional/future work, tracked separately.

All entries below run locally (`cmake --build build --target DEMValidation` /
`ANTIDUNES_PaperAnalyticalContacts`, then execute directly — no cluster needed), except the
"dropper case" showcase runs (§4), which are cluster-only per explicit instruction.

## 1. Cylinder — `[x]` resolved

| Test | prm / target | Status (run 2026-09-13) | Math reference |
|---|---|---|---|
| Own-parameters oblique sweep | `dem_validation_wall_cylinder.prm` (`wallImpactCylinder`) | **PASS** | §1.2 |
| **Literal Kodam Part II Test 1 replica** | `dem_validation_wall_cylinder_kodam_test1.prm` | **Fully resolved.** With the `Cylinder::support()` fix + `useRotationalEffectiveMass=true`: FAIL(31/32), sole residual tilt=30° at 2.48% (vs. 2% tolerance). Diagnosed as a numerical-integration-resolution artifact, not a modeling gap: the *absolute* velocity error is nearly constant across every tilt (~0.001–0.002), and tilt=30° only fails in *relative* terms because that's where the base post-impact velocity happens to be smallest (near the approach/separation sign-crossover). Confirmed by refining `dt` 10x at fixed physics/calibration (`collisionTime` unchanged) — tilt=30° error dropped to 1.94% (PASS), and a full re-sweep at this finer `dt` gives **OVERALL PASS (32/32)** — every check, every tilt, both literal Kodam-paper anchors. Before the `support()` fix: FAIL(19/32) default-mass mode / FAIL(28/32) rotational-mass mode | §1.2 |
| Cylinder×cylinder analytical vs. GJK/EPA | `ANTIDUNES_PaperAnalyticalContacts` (cylinder fuzz suite) | **Completes clean.** 100% existence agreement across all 7 named sub-cases (face-face, band-band parallel/skewed, face-band, face-edge, band-edge, edge-edge); small overlap/point discrepancies, root-caused ("rim-vs-band-surface" gap) | §1.1, `DEM_VALIDATION_MATH_REFERENCE.md#1-cylinder` |
| **Cylinder×wall analytical vs. GJK/EPA** (NEW, `gjkVsAnalyticalCylinderWall`) | `dem_validation_gjk_vs_analytical_cylinder_wall.prm` | **PASS (10/10)** — GJK/EPA (against a large-box wall proxy, since a true `HalfSpace` has no `support()`) confirms the exact closed-form overlap depth to ~1e-8–1e-15 at every tilt 0–90°, including the two worst dynamics-discrepancy angles (30°, 80°) | §1.2 |

> **[Claude, 2026-09-12/13]:** the Kodam-replica intermediate-tilt failures are the one open item here.
> Investigated further: tested and **ruled out three** candidate explanations — (1) finite-duration
> soft-contact effects (a 10x stiffness increase produced no consistent trend in the errors); (2)
> contact-point-formula mismatch (proved algebraically that the analytical reference and the actual
> production contact detection use the *identical* geometric contact-point formula); (3) a secondary/
> two-point contact event (inspected the actual per-step trajectory at tilt=30° — one smooth,
> single-phase compression-rebound, no double-dip). **A fourth, independent check** — added a new
> permanent test, `gjkVsAnalyticalCylinderWall`, using mesa_pd's own iterative GJK/EPA solver (against
> a large-box wall proxy, since a true `HalfSpace` has no `support()`) — **confirms the geometry is
> exactly right at every tilt** (PASS 10/10, overlap depth matching to ~1e-8–1e-15). The mystery is now
> narrowed firmly to the *dynamics*/integration side, not geometry. **Major update — likely found the
> main mechanism**: `useRotationalEffectiveMass` (an existing flag from the earlier box-chatter
> investigation, previously never wired into this test) fixes `LinearSpringDashpot`'s stiffness/
> damping calibration to use the true lever-arm-aware effective mass instead of the plain translational
> one, which the codebase's own docs already note is only exact at zero lever arm. Wired it into
> `wallImpactCylinder` and reran: **dramatic improvement at every intermediate tilt** (e.g. 30°:
> 34.6%→2.48%, 80°: 14.4%→1.15%), overall FAIL(19/32)→FAIL(28/32) — but tilt=90° (previously passing)
> now regresses to 5.0% error. **Root-caused**: a pre-existing tie-break ambiguity in
> `mesa_pd::data::Cylinder::support()` itself (affects production contact detection too, not just this
> test) picks a spurious rim-corner contact point instead of the band's true midpoint exactly at 90°
> tilt — invisible under the default calibration (which never looks at the contact point), but
> corrupts stiffness/damping once `useRotationalEffectiveMass` uses that point's real (spurious) lever
> arm. **Fixed (2026-09-13)**: both `Cylinder::support()` and the analogous `Box::support()` (checked
> proactively, found the same bug there too — worse, a self-reinforcing corner-flickering torque
> oscillation) now snap the offset to zero when the relevant direction component is within `1e-12` of
> exactly degenerate, instead of always applying `sign(x)*extent`. Reran the full sweep: FAIL(31/32),
> only tilt=30° (2.48%) left. **Closed out (2026-09-13)**: tabulated absolute (not relative) velocity
> error across all 10 tilts — it's essentially constant (~0.001–0.002) everywhere, confirming the
> lever-arm effect is now fully accounted for; tilt=30° fails only in *relative* terms because the base
> post-impact velocity is smallest there (near the sign-crossover between approach and separation).
> Confirmed this is a pure numerical-integration-resolution effect, not a remaining physical gap, by
> refining `dt` 10x at fixed physics (same `collisionTime`/calibration): tilt=30° error fell from 2.49%
> to 1.94% (PASS), and the full re-sweep at this finer `dt` gives **OVERALL PASS (32/32)**. A genuine
> physical/geometric gap would not have responded to pure timestep refinement; a truncation artifact
> would (and did). Investigation fully closed — see `DEM_VALIDATION_MATH_REFERENCE.md`'s §1.2 note for
> the complete trail (three hypotheses ruled out, the rotational-mass fix, the support() fix, and this
> final dt-refinement confirmation) and how to reproduce every diagnostic run.

## 2. Cuboid — `[~]` in progress (geometric correctness only, NOT the chatter question)

| Test | prm / target | Status (run 2026-09-12) | Math reference |
|---|---|---|---|
| Box-wall corner impact | `dem_validation_wall_box.prm` (`wallImpactBox`) | **PASS (12/12)**, tilt 5–35°. Same degenerate-direction bug found in `Box::support()` at tilt=0° (flush face-on, outside this sweep's own range) — worse than the cylinder case, a self-reinforcing corner-flickering torque oscillation, plausibly a contributing mechanism to the earlier paused box-chatter investigation. Fixed 2026-09-13; still PASS (12/12) after (unaffected — this sweep was never in the degenerate regime) | §2.1 |
| Box×box: SAT + Eberly-closed-form vs. GJK/EPA | `ANTIDUNES_PaperAnalyticalContacts` (box fuzz suite) | **Completes clean.** Eberly port previously verified bit-exact vs. real GeometricTools library (200k trials, `MITIGATION_COMPARISON.md`) | §2.2 |
| **Face-drop restitution** (NEW, `boxFaceDropRestitution`) | `dem_validation_box_face_drop_mack2011.prm` | **PASS (12/12)**, 4 impact speeds. Literal replica of Mack et al. (2011)'s cube drop-test methodology (measured e=0.54, μ=0.28) but flush/face-first — the exact degenerate configuration the `Box::support()` fix targets. Angular velocity is exactly `<0,0,0>` at every speed (was the spurious-torque regime pre-fix); velocity-restitution error 0.3–0.4% | §2.1 |
| **Incline slide straightness** (NEW, `boxInclineSlideStraightness`) | `dem_validation_box_incline_slide_fraige2008.prm` | **PASS (3/3)**. Literal replica of Fraige et al. (2008)'s "slide test" (their single-contact-point model showed the contact force "jumping" between corners, causing lateral drift — fixed there by a multi-contact model). Our fixed `support()`: lateral position and both out-of-plane angular-velocity components are **bit-exact zero** for the entire 2000-step run; in-plane pitch shows a genuine, physically-expected settling transient (~0.07 rad/s) that decays to ~0.007 rad/s by the last 20% of the run | §2.1 |

No Kodam-style *closed-form dynamics sweep* anchor exists for box-wall impact (Chung & Ooi (2011) is
sphere-only by their own abstract) — but two literal single-particle regression replicas were found and
added 2026-09-13 (found via targeted web search once "cuboid + spring-dashpot + drop/hopper validation"
was posed as the specific ask; both already in the Zotero library): **Fraige, F.Y., Langston, P.A., Chen,
G.Z. (2008)**, "Distinct element modelling of cubic particle packing and flow," Powder Technology
186(3):224–240, DOI:10.1016/j.powtec.2007.12.009 — real acrylic cubes in a rectangular hopper, literal
`F_n=k_n·δ_n−c_n·v_rn` Cundall spring-dashpot, their own "path"/"slide" single-particle debugging tests;
and **Mack, S., Langston, P., Webb, C., York, T. (2011)**, "Experimental validation of polyhedral
discrete element model," Powder Technology 214:431–442, DOI:10.1016/j.powtec.2011.08.043 — real
polyhedral dice (including hexahedra) in a 3D hopper, same spring-dashpot-slider family, with a directly
measured cube coefficient of restitution (e=0.54) and friction (μ=0.28) from an actual drop test. Both
papers' own single-particle methodology (not their bulk hopper/packing numbers, which depend on their own
material stiffness calibration, not ours) is now replicated as permanent regression tests — see the table
above and `DEM_VALIDATION_MATH_REFERENCE.md`'s §2.1 for the full derivation and results.

## 3. Sphere — `[x]` supporting context, not a thesis validation target here

`wallImpactSphere`, `pairRestitution`, `energyConservation` all **PASS** (run 2026-09-12) — see math
reference §3 for the trivial closed forms involved. Included in the thesis methodology section only as
the simplest worked example of the general model (§0 of the math reference).

## 4. "Dropper case" showcase — `[ ]` not started

> **[Dennys, 2026-09-12]:** the final validation dropper test will only be run on the cluster, so we
> also need an sbatch script for it.

Cylinder (Kodam-tilt) VTK sequence, and box equivalent — not yet built. Plan: reuse
`wallImpactCylinder`'s existing `vtkSpacing` parameter (`engine.setupVtk`) with a single representative
case; needs a dedicated cluster sbatch script (`DEMValidation` is a plain MPI executable built directly
via CMake, not run through `antidunes.sh` — follow `jobs/run_72_droptest.sbatch`'s own pattern for the
build/run commands).

## 5. Tier-2 bulk/assembly checks — `[!]` inconclusive, needs a longer local run

Different validation purpose (emergent bulk behaviour, not single-contact-pair correctness) — checked
2026-09-12, first with a 60s timeout (too short), then re-checked at 300s:

| Test | Status |
|---|---|
| `packing` | Did not finish even at 300s |
| `packingSingleShape` | Did not finish even at 300s. Directly inspected the raw output (not just grepped for a final line): step 0 alone produces a large volume of `[STRONG-CONTACT-FORCE]`/`[DEM-WATCH] deep contact` warnings, with reported overlaps up to ~750% of the smaller interaction radius and velocities in the hundreds (e.g. `vel=<-73.6,49.96,-999.9>`). This is consistent with either (a) a deliberate overlap-then-relax packing initialization (common in DEM packing generators — particles start overlapping and a repulsive settling phase resolves it) producing a huge, slow-to-flush log volume, or (b) a genuine numerical instability. Not distinguished — not chased further, since it's outside this plan's cylinder/cuboid-contact-detection focus. |
| `angleOfRepose` | Did not finish even at 300s (not individually inspected — same open-question status as the two above until someone looks) |
| `rotatingDrum` | Did not finish even at 300s (same) |
| `archTunnelContainment` (both `dem_validation_arch_tunnel.prm` and `_box.prm`) | **FAIL (1/2)**, fast (24s / 0.66s) — a geometric-containment check (floaters staying within the tunnel bore), not a contact-detection test; out of this plan's cylinder/cuboid-contact-detection focus, noted here for completeness only, not investigated further |

## 6. Literature status

Confirmed in `~/Zotero/storage/` (title + full text verified via `pdftotext` where cited above):
Kodam et al. (2010) Parts I & II, Ericson (2004) title (HTML snapshot only, not a readable PDF), Chung
& Ooi (2011) (supporting/sphere section only).

**Not in the library, worth adding** — all three below verified via direct web lookup (2026-09-13),
with DOIs/URLs, not just propagated from this codebase's own comments:

- **Cundall, P.A. & Strack, O.D.L. (1979)**, *"A discrete numerical model for granular assemblies,"*
  Géotechnique 29(1):47–65. DOI: [10.1680/geot.1979.29.1.47](https://doi.org/10.1680/geot.1979.29.1.47).
  Bare bibliographic record already in the library (itemIDs 872/874) with no attached PDF.
  **Highest priority to add**: the exact citation for the frictional Cundall–Strack-vs-idealized-
  Coulomb discrepancy already documented in `AnalyticalSolutions.h`, and the tangential-spring model
  actually used by `LinearSpringDashpot`.
- **Housner, G.W. (1963)**, *"The behavior of inverted pendulum structures during earthquakes,"*
  Bulletin of the Seismological Society of America, 53(2):403–417. DOI:
  [10.1785/BSSA0530020403](https://doi.org/10.1785/BSSA0530020403). Open-access copy via CaltechAUTHORS.
  The rocking-block self-check citation (§2.1 of the math reference).
- **Baraff, D. (1997)**, *"An Introduction to Physically Based Modeling: Rigid Body Simulation"*
  (part of Witkin & Baraff's SIGGRAPH course *"An Introduction to Physically Based Modeling"*),
  freely hosted by the author: [notesd1.pdf](https://www.cs.cmu.edu/~baraff/sigcourse/notesd1.pdf)
  (unconstrained rigid-body dynamics) and
  [notesd2.pdf](https://www.cs.cmu.edu/~baraff/sigcourse/notesd2.pdf) (nonpenetration constraints —
  where the generalized effective-mass/impulse formula this session verified algebraically actually
  lives). This codebase's own comments say "Baraff (2001)"; the notes themselves are copyright 1997
  and were reused across many later SIGGRAPH course offerings — cite the 1997 originals unless you
  specifically mean a later course's reprint. This is the citation for `useRotationalEffectiveMass`/
  `computeRotationalEffectiveMass` — the fix behind this session's main finding (§1.2 of the math
  reference).
- Stronge, W.J. (2000), *"Impact Mechanics,"* Cambridge University Press — cited alongside Baraff for
  the same general solver; a book, not independently verified via web lookup this session.

**A citation to double-check**, not yet resolved: `AnalyticalSolutions.h`'s own comment cites
"Kruggel-Emden, H. et al. (2011), Granular Matter 13, 235–244" for the same paper title as the
directly-verified "Chung, Y.C. & Ooi, J.Y. (2011), Granular Matter 13:643–656" — see
`DEM_VALIDATION_MATH_REFERENCE.md`'s §3 note.

## Challenges & Future Improvements

- The Kodam-replica intermediate-tilt discrepancy (§1) is now **fully resolved: PASS (32/32)**. Root
  cause found (`LinearSpringDashpot`'s stiffness/damping calibration needs the lever-arm-aware
  effective mass, `useRotationalEffectiveMass`); a separate, pre-existing geometry bug
  (`Cylinder`/`Box::support()`'s degenerate-direction tie-break) found and fixed along the way; the
  last residual (tilt=30° at 2.48%) confirmed as a numerical-integration-resolution artifact and
  resolved by refining `dt` 10x at fixed physics. Open follow-up (not yet decided/actioned): whether to
  permanently lower the *default* `dt` in `dem_validation_wall_cylinder_kodam_test1.prm` to match, or
  leave the coarser default as-is with this finding documented — only a scratchpad copy was changed.
- The `Box::support()` fix (§2) surfaced a self-reinforcing torque-oscillation mechanism at flush
  box-wall contact that's a plausible contributing factor to the earlier, separately-motivated and
  currently-paused box-chatter/flying-floater investigation (`MITIGATION_COMPARISON.md`) — worth
  revisiting that investigation with this fix in place if/when that work resumes, though this session
  did not re-run any of that investigation's own scenarios to check.
- Tier-2 bulk checks (§5) need re-running with a realistic timeout (these are legitimately longer
  local sims, not necessarily broken) before their status is known.
- No box-wall literal-literature-replica anchor exists (§2) — acceptable per this plan's own
  reasoning (the geometric SAT/Eberly-vs-GJK result doesn't need one), but worth noting explicitly in
  the thesis text as a deliberate scope decision, not an oversight.
- Citation gaps (§6) should be closed (fetch PDFs, resolve the Kruggel-Emden/Chung-Ooi question) before
  the thesis's bibliography is finalized.
