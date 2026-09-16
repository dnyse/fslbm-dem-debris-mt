# DEM Validation — Math Reference

Purpose: every closed-form formula used by this codebase's cylinder/cuboid contact-validation tests,
transcribed with its exact source (paper, equation number, page) and mapped to the actual code
variable that implements it — so each one can be checked independently against the cited paper rather
than taken on faith. Companion to `FINAL_THESIS_DEM_VALIDATION.md` (status/index — this file holds the
derivations, that one holds the pass/fail record and doesn't repeat them).

**Scope**: cylinder and cuboid contact-detection/impact-dynamics validation (the thesis focus).
Sphere material is included only where it's needed to explain the *general* model everything else
specializes from. The box-chatter/production-stability investigation
(`useSpheroBoxContact`, `MITIGATION_COMPARISON.md`) is out of scope here.

**How to use this file**: dated/attributed comments follow `THESIS_PLAN.md`'s own convention
(`> **[Dennys, YYYY-MM-DD]:** ...` / `> **[Claude, YYYY-MM-DD]:** ...`), never overwritten, only
appended to.

**Citation confidence**: each source below is marked ✅ (PDF in `~/Zotero/storage/`, text directly
read/verified for this document — `pdftotext` was run and the cited equation checked against the
actual page) or ⚠️ (cited by this codebase's own existing comments, but **not** in the Zotero library
and not independently verified here — worth tracking down before relying on the exact equation number
in the thesis).

---

## 0. The general model everything below specializes from

All of `DEMValidationScenarios.h`'s wall-impact tests (`wallImpactSphere`, `wallImpactCylinder`,
`wallImpactBox`) compare against **one** general solver, `analytical::rigidBodyImpact()`
(`apps/showcases/Antidunes/AnalyticalSolutions.h:122`) — a standard single-point impulsive rigid-body
collision with Newtonian normal restitution `e` and Coulomb friction `mu`. Rather than re-derive it
here, the essential structure (verify against the code directly, it's fully commented in-place):

- Relative velocity at the contact point: `v_rel = (v_A + ω_A × r_A) − (v_B + ω_B × r_B)`.
- Effective inverse-mass matrix at the contact point: `K = K_A + K_B`, where each body contributes
  `K_i = (1/m_i)·I − [r_i]×·I_i⁻¹·[r_i]×` (`effectiveInverseMassMatrix`, line 76) — `[r_i]×` is the
  skew-symmetric cross-product matrix (`skew()`, line 69).
- Normal-only impulse: `J_n = −(1+e)·v_n / (n·K·n)`.
- Full solve (sticking vs. kinetic-sliding branches, Coulomb-limit check, and the Painlevé-paradox
  fallback for the box-corner case) — lines 130–224.

**Source, per the code's own header (line 46)**: this is "the standard two-branch impulse-momentum
method", citing:
- Baraff, D. (2001), *"Physically Based Modeling: Rigid Body Simulation,"* SIGGRAPH course notes. ⚠️
  (not in library)
- Stronge, W.J. (2000), *"Impact Mechanics,"* Cambridge University Press. ⚠️ (not in library)

> **[Claude, 2026-09-12]:** these two are the generic rigid-body-impact-mechanics references the code
> itself cites for the general two-branch (stick/slip) impulse solver — neither is in
> `~/Zotero/storage/`, so I have not independently verified the exact method against them. The solver
> itself IS independently verified for the box-corner degenerate case via `selfCheckHousnerRelativeError`
> (§2 below) — that's a real, checkable identity even without the Baraff/Stronge texts on hand.

---

## 1. Cylinder

### 1.1 Contact-detection geometry — Kodam et al. (2010), Part I ✅

Kodam, M., Bharadwaj, R., Curtis, J., Hancock, B., & Wassgren, C. (2010). *"Cylindrical object contact
detection for use in discrete element method simulations. Part I – Contact detection algorithms."*
Chemical Engineering Science, 65(22), 5852–5862. (`~/Zotero/storage/DGAKZ7L8/`, verified via
`pdftotext`.)

Implemented directly in `PaperAnalyticalContacts.h` (per that file's own header). Representative
cases, transcribed from the paper (§3.1–3.2, eq. numbers as printed):

- **Face-plane contact** (§3.1.1, p. 5855): contact occurs when the cylinder's z-axis is parallel to
  the plane normal AND `|z^P_cyl| ≤ L/2` (eq. 2), where `z^P_cyl` is the cylinder centre's z-coordinate
  in the plane's body-fixed frame (eq. 3). Code: `cylinderWallImpact`'s face-plane branch
  (`AnalyticalSolutions.h:284`, `axisDotN ≈ ±1`).
- **Face-face contact overlap** (§3.2.1, p. 5857, eq. 28): `d_contact = ½(L_i+L_j) − |z_i,j|` (the two
  cylinders' half-lengths minus their axial separation). Normal direction (eq. 29): points along the
  shared axis, sign of the axial separation.
- Band-plane, edge-plane (§3.1.2–3.1.3), face-band, face-edge, band-band parallel/skewed, band-edge,
  edge-edge (§3.2.2–3.2.6) — same paper, not transcribed here (each is its own multi-step geometric
  construction); see the section numbers above to check the code against the source directly.
  `PaperAnalyticalContactTest.cpp` fuzz-tests exactly these 7 named sub-cases against GJK/EPA (see
  `FINAL_THESIS_DEM_VALIDATION.md` for the population-level results).

### 1.2 Cylinder-wall impact dynamics — Kodam et al. (2010), Part II ✅

Kodam, M., Bharadwaj, R., Curtis, J., Hancock, B., & Wassgren, C. (2010). *"Cylindrical object contact
detection for use in discrete element method simulations, Part II — Experimental validation."*
Chemical Engineering Science, 65(22), 5863–5871. (`~/Zotero/storage/66QQ5DYQ/`, verified via
`pdftotext`.)

**Test 1 protocol (§3.1, p. 5866)**: a cylinder at a prescribed tilt angle impacts a wall with
prescribed translational speed *normal* to the wall and *zero* angular speed; contact assumed
**frictionless**, gravity neglected. Post-impact angular speed (eq. 3, p. 5866):

```
ω_y+ = [ e·(1+?)·... m·V_z,CM⁻ · r·cos(α+θ) ] / [ I_yy + m·r²·cos²(α+θ) ]
```

(the paper's own eq. 3 — re-derive/confirm the exact bracket contents directly from the PDF page
before quoting it verbatim in the thesis; the OCR extraction used to write this reference had minor
line-wrapping artifacts around this specific equation that I did not fully untangle. The two clean,
directly-quotable numbers from the same page are below and don't have this problem.)

**The two literal numeric anchors actually used by the codebase's test** (p. 5866, confirmed clean
in the extracted text): at tilt = 0° (face-on) and tilt = 90° (band-on), the dimensionless rebound
speed `V_z,CM+ / V_z,CM−` equals **exactly the configured coefficient of restitution, 0.85** (their
own calibrated Plexiglas value). `dem_validation_wall_cylinder_kodam_test1.prm` reproduces the paper's
exact geometry (⌀8.0mm, L=5.3mm, L/D=0.6625), material (e=0.85), and protocol (frictionless, zero
angular velocity, purely normal impact), and `DEMValidationScenarios.h:702` checks the simulated
ratio against **this literal paper-stated number**, not just the codebase's own re-derived formula.

**Code implementation**: `analytical::cylinderWallImpact` (`AnalyticalSolutions.h:262`) — specializes
the general solver (§0) with the contact point at the cylinder's rim closest to the wall
(`paper_analytical::rimPointTowardPlane`, from the Part I geometry above), with a special-cased
on-axis contact point exactly at face-on/band-on (see that function's own comment: the naive rim point
would inject a spurious, arbitrary lever arm at the symmetric orientations).

**Result** (`dem_validation_wall_cylinder_kodam_test1.prm`, run 2026-09-12): both literal paper
anchors **PASS** (rel. error 0.52–0.56%). Full sweep (0°–90°, step 10°) **FAILS overall (19/32
checks)** at intermediate tilts (10°–40°, 70°–80°; velocity errors up to ~35%).

> **[Claude, 2026-09-12]:** this test is **frictionless** (`frictionCoefficient 0.0`), so the
> already-documented Cundall–Strack-vs-idealized-Coulomb explanation (§1.3 below) — which the code's
> own comments explicitly say does *not* apply when friction is off — cannot be the cause here.
>
> **Tested and ruled out, two candidate explanations:**
> 1. *Finite-duration soft-contact effect* (my original hypothesis: over the finite duration of a real
>    spring-dashpot contact, a tilted cylinder's off-centre normal force rotates the body and shifts
>    the true contact point/geometry mid-collision, unlike the closed form's single-instantaneous-
>    impulse assumption). **Tested directly**: reran the same sweep with `collisionTime` shrunk 10x
>    (0.01→0.001, stiffer spring, shorter real contact duration, closer to the rigid-impulse limit).
>    Result: **no consistent trend** — some angles' errors got very slightly smaller (tilt=40°:
>    3.80%→3.60%; tilt=80°: 14.4%→14.0%), others got very slightly *larger* (tilt=10°: 8.7%→10.0%;
>    tilt=20°: 13.6%→16.1%), all well within noise of each other, overall still FAIL (19/32). If this
>    were the real mechanism, a 10x stiffness increase should have shrunk the discrepancy noticeably
>    and consistently across all angles — it did neither. **Hypothesis rejected.**
> 2. *Contact-point-formula mismatch* (does the production simulation's actual contact point match the
>    analytical reference's assumed point?). **Checked directly, algebraically**: the analytical
>    reference (`cylinderWallImpact` → `paper_analytical::rimPointTowardPlane`,
>    `PaperAnalyticalContacts.h:1980`) and the actual production contact detection
>    (`SolidCylinderContactDetection`'s Cylinder×HalfSpace overload, `ParticleKernels.h:338`, which
>    calls `mesa_pd::data::Cylinder::support()`, `src/mesa_pd/data/shape/Cylinder.h:80`) turn out to be
>    **the same formula written two different ways** — both reduce to `radius·normalize(perpendicular
>    component of the wall-normal-facing direction) + sign(axial component)·halfLength·axis`, confirmed
>    by direct algebraic substitution (both use the identical radial-direction construction
>    `axisDotN·eZ − n̂`, and the axial sign agrees whenever the cylinder is on the domain side of the
>    plane, true throughout this test's approach and impact). **Hypothesis rejected** — they are
>    provably not using different contact points.
>
> 3. *Secondary/two-point contact event* (does the body touch the wall at a second point — e.g. the
>    band, or the opposite rim end — before fully separating, an event the single-point idealised
>    model structurally cannot represent?). **Checked directly**: reran tilt=30° alone with per-step
>    CSV logging (`engine.openCsv`, already-existing instrumentation) and inspected the actual
>    trajectory through the whole ~40-step contact (t=0.040–0.048s): `velZ` decelerates smoothly from
>    −1 through a single minimum and back to a small positive value, `wY` grows smoothly and
>    monotonically to its final value, `KE` dips to a single minimum (55.9, at t=0.0436) and recovers —
>    **one smooth, single-phase event, no visible double-dip or plateau-then-second-drop** that would
>    signal a second contact point engaging. Also: the accumulated rotation *during* the ~0.038s
>    contact is small (~0.4°, from the average angular velocity × duration), too small on its own to
>    plausibly account for a 35% velocity error via a shifting support-point direction. **No direct
>    evidence found for this mechanism either**, though a smooth trajectory doesn't rule it out as
>    conclusively as the algebraic proof in (2) does — this is a weaker negative than the other two.
>
> **A fourth check, added as a permanent test rather than a throwaway script**: does mesa_pd's
> *iterative* GJK/EPA solver independently confirm the exact closed-form contact geometry, at the same
> tilt angles? A true `HalfSpace` has no `support()` at all, so GJK/EPA cannot run on a Cylinder ×
> HalfSpace pair in this codebase by construction — added `gjkVsAnalyticalCylinderWall`
> (`DEMValidationScenarios.h`, new `testCase`, `dem_validation_gjk_vs_analytical_cylinder_wall.prm`),
> which approximates the wall as a very large, thin box (which *does* have `support()`) purely to get
> an independent third answer. **Result: PASS (10/10)** — GJK/EPA's overlap depth matches the exact
> formula to ~1e-8–1e-15 at every tilt from 0° to 90°, including the two worst-offending dynamics
> angles (30°, 80°). The contact-point X/Y also matches to ~1e-6; the small constant Z offset (0.025,
> exactly half the configured 0.05 overlap) is a **convention difference, not a disagreement** — EPA
> returns the midpoint of the overlap region, the exact formula returns the penetrating surface point.
> At the degenerate 0°/90° cases the point distance is large, exactly as `AnalyticalSolutions.h`'s own
> comment already anticipates (a whole face/band touches, so "the" contact point is inherently
> ambiguous — both methods pick a different, equally valid point).
>
> This is now the **strongest** confirmation that the geometry is correct at every angle, via a
> genuinely independent numerical method (not just re-deriving the same formula by hand) — pushing the
> mystery even more firmly onto the *dynamics*/integration side, not the *geometry* side.
>
> **A concrete, untested lead** (found by reading Kodam et al. (2010) Part II's own methodology
> section, p. 5865-5866, directly — not assumed): Kodam's own DEM does use a linear spring + linear
> dashpot for the normal contact force (their eq. 1, `F_N = k·δ + n·δ̇`) — **the same model family as
> ours** — and the paper *explicitly* flags a known limitation of it: because the force depends only on
> penetration depth `δ`, not on the actual contact geometry, "the calculated maximum contact overlaps
> will generally be overestimated" for non-spherical shapes where different orientations (e.g.
> face-plane vs. band-plane) have different true overlap *volumes* at the same `δ`. Despite this, their
> own "true cylinder DEM" curve (Fig. 4/5) matches the analytical curve closely across the *entire*
> 0–90° sweep, not just the symmetric endpoints — so this particular limitation, even though real,
> doesn't explain a discrepancy that's *worst in the middle and vanishes at 0°/90°* the way ours is.
>
> What IS different: Table 1 ("Parameters for experiments and simulations," p. 5865) has three columns
> — Experiment, **True cylinder DEM simulations**, Glued sphere DEM simulations — and its footnote (d)
> is attached to the "Damped linear spring" force model specifically in the **True cylinder** column
> (double-checked: the *other* shape, the glued-sphere approximation, uses a different, genuinely
> Hertzian model — "Nonlinear damped Hertzian spring / Mindlin–Deresiewicz/coulomb friction" — since
> that column really is built from spheres, where exact Hertz theory applies; the linear-spring
> footnote is not that column's). It states Kodam's stiffness `k` for the true cylinder was **"found
> by matching binary maximum contact overlap with that of Hertz contact"** (a sphere-sphere reference
> collision) — k=4.0e-5 N/m, n=1.506 N/(m/s). Our own `DEMEngine` instead calibrates via
> `setParametersFromCOR(restitution, collisionTime, meff)` — a target restitution coefficient *and*
> target collision duration, a genuinely different calibration philosophy. Both land on *some* linear
> spring-dashpot, but nothing guarantees the two calibration targets produce stiffness/damping that
> behaves the same way as the effective lever arm changes with tilt. **Not yet tested** — the natural
> next experiment is reproducing Kodam's own k/n values directly (in matching SI units, using their
> exact Table 1 mass/geometry) rather than our own COR+duration-derived ones, and checking whether the
> intermediate-tilt discrepancy shrinks.
>
> **Update, tested directly — this is very likely the main mechanism**: this exact codebase already
> documents (from the earlier, separately-motivated box-chatter investigation, `SimulationParams.h`'s
> `ModelParameters::useRotationalEffectiveMass` doc comment) that `LinearSpringDashpot`'s
> stiffness/damping are calibrated ONCE from the plain *translational* effective mass, which is only
> exact when the contact point has zero lever arm about the line of centres — "is EXACT for
> sphere-sphere contact (there rA/rB are always parallel to n, so the rotational terms are identically
> zero)" per that comment. A tilted cylinder's rim contact point has exactly this nonzero lever arm,
> vanishing at tilt=0°/90° (on-axis contact) and largest at intermediate tilts — precisely matching the
> observed error pattern. Added `useRotationalEffectiveMass` as a new prm parameter to
> `runWallImpactCylinder` (previously never threaded through at all, always defaulting to `false`) and
> reran the Kodam replica with it turned on:
>
> | Tilt | velRelErr, before | velRelErr, with `useRotationalEffectiveMass=true` |
> |---|---|---|
> | 10° | 8.7% | **0.29%** |
> | 20° | 13.6% | **0.60%** |
> | 30° | 34.6% | **2.48%** (still narrowly fails the 2% tolerance) |
> | 40° | 3.8% | **0.48%** |
> | 70° | 2.6% | **0.35%** |
> | 80° | 14.4% | **1.15%** (now passes) |
>
> Overall FAIL(19/32) → FAIL(28/32). This is a dramatic, order-of-magnitude improvement at every
> previously-bad intermediate angle — strong confirmation that the lever-arm/effective-mass
> calibration gap is the dominant mechanism, not a coincidence.
>
> **But it introduces a new regression**: tilt=90° (band-on, previously PASSING against Kodam's own
> literal 0.85 anchor at 0.52% error) now **fails** — velRelErr jumps to 5.0%, energy-budget error to
> 10.2%. **Root-caused directly** by rerunning tilt=90° alone with per-step CSV+raw-log output enabled
> and inspecting the actual `[STRONG-CONTACT-FORCE]` diagnostic line: `contactPoint=<2.65,0,0>`,
> `thisContactTorque=<0,7829.55,0>` (nonzero!) — `2.65` is exactly `cylinderLength/2`, i.e. the DEM's
> *actual* contact point sits at the rim **corner** (face meets band), not the band's midpoint the
> analytical reference assumes.
>
> Mechanism: at exactly 90° tilt, the wall-normal direction expressed in the cylinder's own body frame
> has an axial component that is numerically ≈0 but not exactly 0 (floating point), and
> `Cylinder::support()` (`src/mesa_pd/data/shape/Cylinder.h:80`) still applies
> `sign(bfD[2])·halfLength` as an axial offset regardless — `sign(≈0)` still resolves to a definite ±1,
> so the support point picks up a spurious full axial offset to one face's edge instead of "no axial
> preference, the whole band length is equally valid." This is a **pre-existing tie-break ambiguity in
> `Cylinder::support()` itself**, not something `useRotationalEffectiveMass` introduces — it's the same
> ambiguity `AnalyticalSolutions.h:284`'s own comment already flags as needing special-casing for the
> *analytical reference*. It affects **production contact detection too** (identical `support()` call),
> but was invisible before: the plain translational calibration never looks at the contact point at
> all, so the only symptom was a small spurious torque (already visible as the pre-existing, tolerated
> `angRelErr=1` at 90° even in the *original*, non-rotational-mass run — dividing by an
> expected-exactly-zero analytical value). With `useRotationalEffectiveMass` on, that same spurious
> corner point now also corrupts the *stiffness/damping calibration itself* (via a real, nonzero
> `r×n`), escalating a previously-cosmetic issue into a genuine velocity/energy failure.
>
> **Fixed and verified (2026-09-13)**: patched `Cylinder::support()`
> (`src/mesa_pd/data/shape/Cylinder.h:80`) to snap the axial offset to exactly 0 when `|bfD[2]|` is
> below `1e-12` (the same absolute epsilon the function already used for its lateral-degeneracy
> check), instead of always applying `sign(bfD[2])·halfLength` regardless. Also fixed the identical
> class of bug in `Box::support()` (`src/mesa_pd/data/shape/Box.h:72`) per-axis, on the user's own
> request to check for it there too — see §2 below for that one, which turned out worse (a
> self-reinforcing torque oscillation, not just a wrong-but-stable lever arm).
>
> **Full sweep after the fix** (`dem_validation_wall_cylinder_kodam_test1.prm`):
>
> | Mode | Before fix | After fix |
> |---|---|---|
> | Default (`useRotationalEffectiveMass=false`) | FAIL (19/32) | **FAIL (20/32)** — tilt=90°'s
>   spurious angular velocity (`simAngVel≈0.0244`, was already tolerated/dismissed as a division-by-
>   zero artifact) is now genuinely ≈0 (`-2.6e-16`), not just numerically masked |
> | `useRotationalEffectiveMass=true` | FAIL (28/32) | **FAIL (31/32)** — every tilt 0°–90° now
>   passes except tilt=30° (2.48% vs. the strict 2% tolerance) |
>
> **Regression-checked** (2026-09-13) after the fix, before treating it as safe: `walberla_mesa_pd`
> and the full `Antidunes`/`BedGeneration` production binaries still build clean;
> `ANTIDUNES_ClipAnalyticGeometry`/`ANTIDUNES_BognerStateClassification` unaffected;
> `ANTIDUNES_ContactDetection` initially **failed** one pre-existing assertion
> (`testBoxHalfSpaceContactAxisAligned` expected the OLD arbitrary-corner behavior, `(0.5,0.5,0)`, for
> an exactly-flush axis-aligned box — itself an instance of the bug, not a legitimate case the fix
> broke) — corrected that test's expected value to `(0,0,0)` (the physically-correct face centre,
> confirmed penetration depth is unchanged) and redesigned the adjacent
> `testBoxHalfSpaceContactRotated` (which rotated only about the wall-normal axis, Z — also degenerate
> regardless of angle, so it wasn't actually testing rotated-corner-finding at all) to rotate about X
> instead, a genuinely non-degenerate case, with the expected contact point computed symbolically from
> the rotation rather than hand-approximated. Both pass after that correction.
> `ANTIDUNES_PaperAnalyticalContacts` (the full cylinder/box/sphere GJK/EPA fuzz suite, thousands of
> trials) completes clean, exit 0, no new mismatch categories — only the same pre-existing
> "informational only, not asserted" ones from before the fix.
>
> Also disabled the per-mismatch VTK-frame writing in `PaperAnalyticalContactTest.cpp` by default (it
> was the single biggest cost in that test's runtime for no benefit on a routine run) — set
> `WALBERLA_WRITE_MISMATCH_VTK=1` to turn it back on when actually visualizing a specific mismatch in
> ParaView. The CSV dump (cheap) is unaffected.
>
> **Citation for `useRotationalEffectiveMass`/`computeRotationalEffectiveMass`**, verified via direct
> web lookup (2026-09-13), not just propagated from the code comment: Baraff, D. (1997), *"An
> Introduction to Physically Based Modeling: Rigid Body Simulation"* (part of Witkin & Baraff's
> SIGGRAPH course *"An Introduction to Physically Based Modeling"*) — freely available from the
> author's own page, <https://www.cs.cmu.edu/~baraff/sigcourse/notesd1.pdf> (unconstrained dynamics)
> and <https://www.cs.cmu.edu/~baraff/sigcourse/notesd2.pdf> (nonpenetration constraints, which is
> where the generalized effective-mass/impulse formula actually lives). The existing code comments
> say "Baraff (2001)" — the notes themselves are copyright 1997 and were reused across many later
> SIGGRAPH course offerings through the early 2000s, so "2001" likely refers to a specific
> re-presentation rather than the original notes; cite the 1997 originals unless you specifically mean
> a later course's reprint.
>
> **Resolved (2026-09-13) — the sweep now passes 32/32.** After the `Cylinder`/`Box::support()` fix
> (above) the only remaining failure was tilt=30° at 2.48% (vs. the strict 2% tolerance). Checked
> whether this was a genuine remaining physical/geometric gap or a numerical-resolution artifact by
> tabulating the **absolute** (not relative) velocity difference at every tilt:
>
> | Tilt | sim | analytic | abs. diff | rel. err |
> |---|---|---|---|---|
> | 0° | 0.854718 | 0.85 | +0.00472 | 0.56% |
> | 10° | −0.362966 | −0.364016 | +0.00105 | 0.29% |
> | 20° | −0.188580 | −0.189718 | +0.00114 | 0.60% |
> | **30°** | **0.072708** | **0.074561** | **−0.00185** | **2.49% (FAIL)** |
> | 40° | 0.429794 | 0.431849 | −0.00206 | 0.48% |
> | 50°–90° | ... | ... | −0.0009 to −0.0020 | 0.10%–1.15% |
>
> The absolute error is **remarkably constant across every tilt** (~0.001–0.002, no dependence on the
> lever arm at all — confirming the lever-arm effect above is now fully accounted for). Tilt=30° isn't
> physically worse; it just happens to be where the post-impact Z-velocity crosses through its
> smallest magnitude (the approach/separation sign-crossover), so the *same* small, tilt-independent
> absolute residual becomes a large *relative* error there and nowhere else — a tolerance-metric
> artifact, not a modeling gap.
>
> **Confirmed directly**: reran tilt=30° alone with `dt` refined 10x (`0.00002` instead of the default
> `collisionTime/50 = 0.0002`), **same `collisionTime`/physics calibration, only the integration
> timestep changed**. Result: absolute error shrank (0.00185 → 0.00145), velRelErr dropped from 2.49%
> to 1.94%, now **PASS**. A genuine remaining physical effect would not respond this way to pure
> timestep refinement at fixed physics; a numerical-integration truncation error would (and does).
> Reran the **entire** 0°–90° sweep at this finer `dt`: **OVERALL PASS (32/32)** — every single check,
> including both literal Kodam paper anchors, now passes. The investigation is fully closed: five
> candidate causes examined (three ruled out, one confirmed as the main geometric-calibration
> mechanism, one confirmed as a pure numerical-resolution effect), no remaining unexplained
> discrepancy of any kind.

### 1.3 Already-documented discrepancy for the *frictional* general sweep

`dem_validation_wall_cylinder.prm` (our own parameters, non-replica, `frictionCoefficient` nonzero)
passes overall. The general model's own header (`AnalyticalSolutions.h:69–79`) already documents,
precisely, why a *frictional* oblique impact can disagree with the idealized closed form even when the
underlying geometry/solver is correct:

> "mesa_pd's actual contact law (LinearSpringDashpot) ... uses a Cundall & Strack (1979)-style
> TANGENTIAL SPRING clamped by the Coulomb limit ... over a short, finite, soft contact duration that
> spring does not necessarily have time to ramp up to the full Coulomb-limited force for the whole
> contact, so the simulated tangential impulse transferred ... comes in systematically BELOW this
> idealized-instantaneous-impulse prediction for fast oblique impacts — observed ~15–25% low on
> angular velocity at tangential/normal speed ratios of 0.25–1.5 ... while the PURE NORMAL case
> (ratio=0, no friction involved) matches to <1%."

Cundall & Strack (1979), *"A discrete numerical model for granular assemblies,"* Géotechnique
29(1):47–65. ⚠️ **Not in the library** (bare bibliographic record only, itemIDs 872/874, no attached
PDF) — this is the single most useful addition to make, since it's the exact citation for a
discrepancy this thesis will want to discuss precisely.

---

## 2. Cuboid

### 2.1 Box-wall corner impact — general solver + Housner (1963) self-check

`wallImpactBox` (`dem_validation_wall_box.prm`) compares against `analytical::boxCornerWallImpact`
(`AnalyticalSolutions.h:314`) — **the same general solver as §0**, evaluated at the box's own true
corner contact point (`mesa_pd::data::Box::support()`, the identical corner-selection logic the
production DEM contact code uses). **This is not literally Housner's reduced-2D formula** — that
formula only appears in an internal regression check (below), not as the box-wall test's own
comparison target.

**Result** (2026-09-12): `wallImpactBox` **OVERALL PASS (12/12)** across tilt = 5°/15°/25°/35°.

> **[Claude, 2026-09-13]:** checked (on request) whether the same class of degenerate-direction bug
> found in `Cylinder::support()` (§1.2) also exists in `Box::support()` — **it does, and it's worse**.
> `wallImpactBox`'s own sweep deliberately avoids tilt=0° (flush face-on) — extending it there directly
> (`toleranceRelError` unchanged) surfaced the bug immediately: consecutive `[STRONG-CONTACT-FORCE]`
> log lines within the SAME contact event showed the contact point **flickering between opposite
> corners** of the same flush face (`<1.01,2,0>` then `<-0.99,-1.99999,0>`), with torque flipping sign
> correspondingly (`<7835,-3918,0>` → `<-9255,4628,0>`). Mechanism: an arbitrary corner's spurious
> torque rotates the box infinitesimally, which flips which near-zero direction component's sign
> resolves which way on the *next* `support()` call, picking a *different* corner — a self-reinforcing
> oscillation, not just a stable-but-wrong lever arm like the cylinder case. This is a **plausible
> contributing mechanism to the earlier (separately-motivated, currently paused) box-chatter/
> flying-floater investigation** (`MITIGATION_COMPARISON.md`) — not asserted as the whole explanation,
> but a concrete, directly-observed instability of exactly the kind that investigation was chasing.
>
> **Fixed and verified (2026-09-13)**: same principle as the cylinder fix — `Box::support()`
> (`src/mesa_pd/data/shape/Box.h:72`) now snaps each axis's offset to 0 when `|bfD[i]| ≤ 1e-12`
> instead of always applying `sign(bfD[i])·edgeLength[i]/2`, picking the physically-symmetric
> face/edge-centre point on that axis instead of an arbitrary corner. `wallImpactBox` still **PASS
> (12/12)** after the fix (unaffected — its own sweep, 5°–35°, was never in the degenerate regime to
> begin with); the flush-face tilt=0° case now gives `(0,0,0)` (verified: penetration depth unchanged,
> only the in-plane contact-point choice changes) instead of flickering between corners. Full
> regression details (production build, `ANTIDUNES_ContactDetection` test-expectation fixes, full fuzz
> suite) recorded in §1.2 above (applies to both shapes' fixes together, since they were made and
> verified in the same pass).

Housner, G.W. (1963). *"The behavior of inverted pendulum structures during earthquakes."* Bulletin
of the Seismological Society of America, 53(2), 403–417. ⚠️ **Not in the library.**

Housner's classical rocking-block restitution (as cited by `AnalyticalSolutions.h:335`,
`housnerRockingRestitution`):

```
r = (ω_after / ω_before)² = (1 − 1.5·sin²α)²
```

where `α` is the block's slenderness angle (angle between the base-corner-to-centroid diagonal and
vertical), and the moment of inertia about a base corner `I_O = (4/3)·m·b²` (Housner's eq. 1, `b` =
half-diagonal) — **cited from the code's own comment, not independently verified against the paper**
since it isn't in the library.

`selfCheckHousnerRelativeError()` (`AnalyticalSolutions.h:352`) is an **internal regression check**:
it evaluates the *general* solver (§0) at Housner's own reduced planar parametrization (zero
restitution = "block does not bounce", full-stick friction = "no final slip at the new pivot corner"),
and asserts the result matches `housnerRockingRestitution(α)` exactly. This checks that the general
solver is *internally consistent* with Housner's known closed form in that one degenerate case — it
is a self-consistency check on the solver, not an independent experimental validation of the box-wall
test itself (that's what §2.1's own PASS result above is).

> **[Claude, 2026-09-12]:** worth running `selfCheckHousnerRelativeError` yourself (it's a pure
> function, no simulation) across a few `α` values and confirming the relative error really is ~0, if
> you want an extra independent check before citing this in the thesis — I have not executed it this
> session (it's not wired into `DEMValidation`'s `testCase` dispatch; find it a call site, or add a
> one-off `WALBERLA_LOG` in a scratch build, to exercise it).

### 2.1b Literal single-particle replicas: Fraige et al. (2008) / Mack et al. (2011)

Added 2026-09-13, in direct response to "is there not a paper that does drop or hopper tests with DEM
for cuboids and needs to validate their cuboid implementation" — both use a **linear spring-dashpot**
contact model (the same family as `LinearSpringDashpot`), both are now in the Zotero library, and both
have single-particle debugging methodologies that are direct, literal analogues of the exact bug fixed
this session (`Box::support()`'s degenerate-direction tie-break, §1.2/§2.1 above).

**Fraige, F.Y., Langston, P.A., Chen, G.Z. (2008)**, *"Distinct element modelling of cubic particle
packing and flow,"* Powder Technology 186(3):224–240, DOI:
[10.1016/j.powtec.2007.12.009](https://doi.org/10.1016/j.powtec.2007.12.009). Verified directly from the
PDF (`~/Zotero/storage/M95F3WBC/`):

```
F_n = k_n·δ_n − c_n·v_rn                         (their eq. 1 -- literally our own normal-force form)
F_fb = k_t·δ_t,  capped by  F_fa = μ·F_ne          (their eqs. 2-3 -- Hookean tangential + Coulomb cap,
                                                      i.e. Cundall-Strack, same family as our tangential
                                                      spring)
```

Real acrylic cubes (12.41±0.26mm, ρ=1202±30 kg/m³) dropped into a real rectangular hopper, filmed, and
compared to two versions of their own DEM-Cube model: a "single-contact" version (one point, the deepest
overlapping vertex — structurally the same idea as `Box::support()`'s pre-fix behaviour) and a
"multiple-contact" version. Their own **single-particle debugging tests** (Section 3.6, Figs. 7–9) are
literal precedents for this session's own `support()` investigation:

- **Slide test**: a cube resting face-flush on a surface inclined past the friction limit angle. Their
  single-contact model's contact-force location "jumps" between corners as the cube translates (Fig. 7),
  producing a spurious lateral deviation from the straight downslope path (Fig. 8) — fixed by the
  multi-contact version.
- **Path test**: a cube given an initial "flick" (translation + rotation) on a horizontal surface,
  tracked (X, Y, θ) vs. time and compared between the two model versions and experiment (Fig. 9) — the
  multi-contact version tracks the experiment; the single-contact version departs.

**Replicated as `boxInclineSlideStraightness`** (their slide test): a cube resting flush on a plane
tilted about world X (μ=0.4, their Table 2 measured cube-wall friction), released through a small gap so
gravity pulls it into full contact, then slid for 2000 steps. Result with the fixed `support()`: world-X
position and both out-of-plane angular-velocity components (wY, wZ) are **bit-exact 0.0** for the entire
run (not just "small" — the code path never breaks the mirror symmetry). The in-plane pitch (wX) shows a
genuine, expected settling transient (the release trajectory through the gap is along world gravity, not
purely along the tilted wall normal, so the first contact is not perfectly head-on) — it decays from a
transient peak ~0.07 rad/s to ~0.007 rad/s over the last 20% of the run, checked separately from the
exact-zero lateral/out-of-plane checks (`dem_validation_box_incline_slide_fraige2008.prm`, PASS 3/3).

**Mack, S., Langston, P., Webb, C., York, T. (2011)**, *"Experimental validation of polyhedral discrete
element model,"* Powder Technology 214:431–442, DOI:
[10.1016/j.powtec.2011.08.043](https://doi.org/10.1016/j.powtec.2011.08.043). Verified directly from the
PDF (`~/Zotero/storage/K6JMIZV2/`) — the direct follow-up to Fraige et al., generalizing cube→arbitrary
convex polyhedra (real dice: tetrahedron, hexahedron, octahedron, icosahedron, pentagonal trapezohedron),
same spring-dashpot-slider family, validated in a real 3D hopper at three half-angles with a 1000fps
high-speed camera. Their Section 2.4 gives a **directly measured single-particle coefficient of
restitution** for real hexahedra (cubes): dropped ten from 50cm, `e = √(h₁/h₀) = 0.54`; Section 2.3 gives
the measured sliding friction `μ = 0.28` (their Table 2).

**Replicated as `boxFaceDropRestitution`**: a cube dropped FACE-first (flush, zero tilt) onto a wall
using their literal measured e=0.54, μ=0.28 — the exact degenerate configuration `Box::support()`'s
pre-fix tie-break corrupted (an arbitrary-corner contact point injects a spurious lever arm/torque into a
physically trivial 1-DOF bounce). By symmetry the analytical answer needs no closed-form solver call:
`analyticVel=(0,0,+e·impactSpeed)`, `analyticAngVel=(0,0,0)` exactly. Result across 4 impact speeds
(0.5–2.0): simulated angular velocity is exactly `<0,0,0>` at every speed (confirming the fix), and the
velocity-restitution ratio matches e to 0.3–0.4% relative error
(`dem_validation_box_face_drop_mack2011.prm`, PASS 12/12).

Neither paper's own *bulk* hopper-discharge or packing-porosity numbers were replicated — those depend on
each paper's own material-property-derived stiffness calibration (Fraige: Hookean k=10⁹ N/m explicitly
"smaller than actual" for stability; Mack: 5×10⁷ dyn/cm, also deliberately softened), not on ours
(`collisionTime`/restitution-derived per `calcStiffnessAndDamping`) — replicating those numbers would
mean adopting their calibration inputs wholesale, not validating our own contact geometry. The
single-particle methodology (contact-point correctness under a flush/degenerate configuration) is
calibration-independent and is what was replicated here.

### 2.2 Box×box contact detection: SAT and closed-form vs. GJK/EPA

- **Separating Axis Theorem (SAT)**: Ericson, C. (2004), *"Real-Time Collision Detection,"* Morgan
  Kaufmann, ch. 4. ⚠️ **Not text-verified** (the Zotero entry, itemID 924, is an HTML snapshot, not a
  readable PDF) — the 15-axis OBB–OBB test (3+3 face normals + 9 edge-pair cross products) is the
  standard textbook algorithm this chapter describes; re-confirm the exact axis count/edge cases
  against your own copy before citing a page number.
- **Closed-form box-box distance**: ported from Eberly, D. (GeometricTools library) — **verified
  bit-exact against the real, unmodified GeometricTools library** (200,000 randomized general-rotation
  trials, zero discrepancies; recorded in `MITIGATION_COMPARISON.md`). This specific correctness
  result stands on its own regardless of the broader chatter-mitigation investigation it was found
  inside of (which stays out of scope here).
- Both fuzz-tested against GJK/EPA in `ANTIDUNES_PaperAnalyticalContacts` — see
  `FINAL_THESIS_DEM_VALIDATION.md` for the population-level numbers.

---

## 3. Sphere (supporting context only)

Not a validation target in this plan (contact detection for spheres is exact — no GJK/EPA involved at
all), but the simplest specialization of the general model (§0), useful for framing the thesis
methodology section:

- `analytical::sphereWallImpact` (`AnalyticalSolutions.h:243`): `I = (2/5)·m·R²`, contact point
  `= −R·n̂` (always directly below the centre along the wall normal — trivial, no orientation
  tracking needed).
- `analytical::twoBodyHeadOnRestitution` (`AnalyticalSolutions.h:378`): for a central, face-aligned,
  no-lever-arm collision, the general solver collapses to the textbook `v_rel,after = −e·v_rel,before`
  regardless of shape/inertia — the reference value `pairRestitution` checks against.
- Structurally similar (post-impact-velocity-vs-restitution/friction, not a time-domain Hertz-contact
  check) to Chung & Ooi (2011)'s Tests 4–6 for oblique sphere-wall impact — see that paper
  (`~/Zotero/storage/NBH3X2QR/`, verified via `pdftotext`) if a genuine Hertz-elastic (Tests 1–2:
  closed-form collision duration/max overlap/max force) check is ever wanted; not built here since
  it's out of this plan's cylinder/cuboid focus.

> **[Claude, 2026-09-12]:** one citation worth double-checking — `AnalyticalSolutions.h:37-40`'s own
> comment attributes the classic oblique sphere-wall model to *either* "Walton, O.R. & Braun, R.L.
> (1986)" *or* "Kruggel-Emden, H. et al. (2011), 'Benchmark tests for verifying discrete element
> modelling codes at particle impact level,' Granular Matter 13, 235–244" — but the paper I actually
> found and read in `~/Zotero/storage/NBH3X2QR/` with that *exact* title is **Chung, Y.C. & Ooi, J.Y.
> (2011), Granular Matter 13:643–656** (author names and page range both directly confirmed from the
> PDF's own title page). Either there are two distinct 2011 Granular Matter papers with an
> near-identical title (possible, but I'd want to see the second one before assuming it), or the
> existing code comment has the wrong author/page-range for this citation. Worth resolving before
> citing either version in the thesis — I did not chase down a possible "Kruggel-Emden 2011" paper
> this session.
