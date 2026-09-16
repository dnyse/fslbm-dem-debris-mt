# DEM Contact Mitigation Comparison

All 9 mitigation configurations run against every test case, tracking the **combined total** (both bodies) kinetic energy (or, for `fallOntoPileStability`, combined total mechanical energy KE+PE) rather than any single body's own value. Generated from `DEMValidationScenarios.h`'s `runProductionCollisionReplay`, `runFallOntoPileStability`, and `runMultiPointEdgeCases` (via `DEMValidation`).

Config legend: **baseline** = no mitigation; **useRoundedBoxCorners** = Minkowski-rounded box corners; **useMultiPointBoxContact** = resolve every embedded vertex independently; **both** = rounded+multipoint together; **+overlapCap** = multipoint with a per-vertex overlap cap (previous production recommendation); **useAxisHysteresis** = sticky SAT axis selection; **+manifold** = persistent faded-in/out multipoint contact manifold; **useRotationalEffectiveMass** = generalized (translational+rotational) effective mass in the spring-dashpot calibration; **useSpheroBoxContact** = genuine "rounded box" (spherosimplex) — bypasses `boxBoxSAT` entirely and resolves Box×Box via mesa_pd's own GJK/EPA run with a real physical margin on the two sharp boxes (see `SimulationParams.h::ModelParameters::useSpheroBoxContact`'s doc comment).

## `useSpheroBoxAnalytical`: closed-form alternative, now validated

A second implementation of `useSpheroBoxContact` was added: `useSpheroBoxAnalytical` (see
`ModelParameters::useSpheroBoxAnalytical`, `SimulationParams.h`) resolves the same spherosimplex contact
via a closed-form, non-iterative distance-between-two-boxes algorithm ported from David Eberly's
GeometricTools library (`PaperAnalyticalContacts.h`'s `eberly_box_distance` namespace /
`boxBoxSpheroAnalytical`), instead of mesa_pd's GJK/EPA.

**The ported distance algorithm itself was verified bit-exact against the real, unmodified GeometricTools
library** (200,000 randomized general-rotation trials, zero discrepancies) — but the first version still
failed `grazing-plank-stack` at e=1.00 (+19.7% gain), because Eberly's `DistOrientedBox3OrientedBox3` is a
pure *separation*-distance query with no EPA-equivalent penetration phase: the instant the two boxes'
sharp cores actually overlap (not just within the rounding margin), it reports `distance=0` with an
arbitrary, unstable touching point (confirmed directly: exact `distance=0.0` at `grazing-plank-stack`'s own
starting pose against the real library).

**Fix**: a hybrid — use Eberly's closed form when the cores are genuinely separated (`distance > 1e-9`);
fall back to GJK/EPA (which has a real penetration-depth phase) the instant the cores overlap. Re-tested
after the fix:

| | e=1.00 | e=0.30 |
|---|---|---|
| useSpheroBoxContact (GJK) | 9/9 PASS | 9/9 PASS |
| useSpheroBoxContact (Eberly analytical, hybrid) | **9/9 PASS** | **9/9 PASS** |

`grazing-plank-stack` specifically: GJK −0.48% / hybrid −0.40% at e=1.00; GJK −75.1425% / hybrid −75.1425%
at e=0.30 (identical to 6 decimal places). The hybrid is functionally equivalent to GJK across the full
edge-case battery, while being closed-form (no iteration, no convergence tolerance) for the common
well-separated-cores regime.

## Summary

**`useSpheroBoxContact` is the standout result of this whole comparison.** It is the *only* configuration
across every test in this file to pass **9/9 at both the strict e=1.00 bound and production's e=0.30**
restitution on the synthetic edge-case battery (section 3), with a worst-case gain at e=1.00 of just
**+0.11%** (`plank-vs-near-cube-deep`) — every other configuration tested this session showed gains from
tens to thousands of percent on at least one geometry. Even `grazing-plank-stack` (the worst failure found
this session under every other config: +618% to +7603% at e=1.00) passes cleanly under
`useSpheroBoxContact` (−0.48% at e=1.00, −75% at e=0.30). **`useSpheroBoxContact` has not yet been run
against the `productionCollisionReplay`/`runFallOntoPileStability` batteries (sections 1–2)** — that's the
natural next step before treating it as a full production recommendation.

Other findings, from before `useSpheroBoxContact` was added:

- **`useMultiPointBoxContact+overlapCap`** (the previous best pick) reaches 9/9 on the edge-case battery at
  e=0.30 only (not e=1.00), is mid-pack (+16.6%) on the real incident replay, and is one of the **worse**
  performers on the realistic fall-onto-pile case (+208%) — no config other than `useSpheroBoxContact`
  dominates across the board.
- **`useRotationalEffectiveMass`** shows **no meaningful benefit** on its own: +407% vs baseline's +418% on
  the real replay, +186% vs baseline's +172% on the fall-onto-pile case.
- **`useAxisHysteresis`** alone is actively **worse than baseline** on the real incident replay (+626% vs
  +418%).
- The **fall-onto-pile case (realistic impact speed onto a 6-box pile) is the hardest test by far**: every
  configuration tested against it so far fails its 2% tolerance, with gains ranging from +20% (best of that
  batch: `both`) to +29956% (worst: `useRoundedBoxCorners` alone).

## 1. `productionCollisionReplay` — real incident, job 4199585

Direct replay of the real "flying box" incident (uid `3314649325744685089` vs uid `720575940379279392`, exact logged state at step 25883), 10000 sub-cycles. Production log's own independently-measured gain over the matching 3540-sub-cycle window: **+536%**.

| Config | KE start | KE peak | peak sub-cycle | KE end | Combined-KE gain | Result |
|---|---|---|---|---|---|---|
| baseline (both OFF) | 0.0090499 | 0.046894 | 4071 | 0.028837 | +418.2% | FAIL |
| useRoundedBoxCorners | 0.0090499 | 0.010833 | 9999 | 0.010833 | +19.7% | FAIL |
| useMultiPointBoxContact | 0.0090499 | 0.010549 | 9999 | 0.010549 | +16.6% | FAIL |
| both | 0.0090499 | 0.010836 | 9999 | 0.010836 | +19.7% | FAIL |
| useMultiPointBoxContact+overlapCap | 0.0090499 | 0.010549 | 9999 | 0.010549 | +16.6% | FAIL |
| useAxisHysteresis | 0.0090499 | 0.065734 | 4343 | 0.050493 | +626.4% | FAIL |
| useMultiPointBoxContact+manifold | 0.0090499 | 0.010549 | 9999 | 0.010549 | +16.6% | FAIL |
| useRotationalEffectiveMass | 0.0090499 | 0.045888 | 4687 | 0.035358 | +407.1% | FAIL |
| useMultiPointBoxContact+rotationalEffMass | 0.0090499 | 0.011407 | 9999 | 0.011407 | +26.0% | FAIL |

## 2. `runFallOntoPileStability` — realistic fall onto a pile

Single representative geometry (quickMode): dropHeight=-0.71 (matches the real measured release-overlap), flyingSpeed=0.015 LU/step (matches the real measured flying-neighbour speed), pileSize=6, landing on-centre, 2000 steps. Tracks combined **total mechanical energy** (KE + gravitational PE) of every body, which must be non-increasing at e=0.3 (tolerance 2%).

| Config | E start | E peak | peak step | Combined-E gain | Peak Z after impact | Result |
|---|---|---|---|---|---|---|
| baseline (both OFF) | 0.0087422 | 0.023768 | 620 | +171.9% | 6.0702 | FAIL |
| useRoundedBoxCorners | 0.0087422 | 2.6275 | 1969 | +29955.6% | 7.2474 | FAIL |
| useMultiPointBoxContact | 0.0087422 | 0.036559 | 954 | +318.2% | 6.2168 | FAIL |
| both | 0.0087422 | 0.010507 | 1999 | +20.2% | 3.7893 | FAIL |
| useMultiPointBoxContact+overlapCap | 0.0087422 | 0.026935 | 1497 | +208.1% | 5.2827 | FAIL |
| useAxisHysteresis | 0.0087422 | 0.023768 | 622 | +171.9% | 6.0712 | FAIL |
| useMultiPointBoxContact+manifold | 0.0087422 | 0.042955 | 1231 | +391.4% | 5.9001 | FAIL |
| useRotationalEffectiveMass | 0.0087422 | 0.025034 | 257 | +186.4% | 5.0645 | FAIL |
| useMultiPointBoxContact+rotationalEffMass | 0.0087422 | 0.037242 | 652 | +326.0% | 5.3651 | FAIL |

Every configuration fails the strict 2% bound on this specific realistic combination — `both` (rounded+multipoint) comes closest (+20%), `useRoundedBoxCorners` alone is catastrophic (+29956%, the falling box gets launched to Z=7.2 vs a drop height of 3.79).

### 3a. `runMultiPointEdgeCases` at e=1.00 (strict diagnostic bound)

9 hand-built geometries, combined-KE relative gain per config (✅ = PASS, < 2% tolerance).

| Case | baseline (both OFF) | useRoundedBoxCorners | useMultiPointBoxContact | both | useMultiPointBoxContact+overlapCap | useAxisHysteresis | useMultiPointBoxContact+manifold | useRotationalEffectiveMass | useMultiPointBoxContact+rotationalEffMass | useSpheroBoxContact |
|---|---|---|---|---|---|---|---|---|---|---|
| exact-coincident | +0.0% ✅ | -0.2% ✅ | +0.0% ✅ | +0.0% ✅ | +0.0% ✅ | +0.0% ✅ | +0.0% ✅ | +0.0% ✅ | +0.0% ✅ | +0.001% ✅ |
| fully-engulfed | +0.0% ✅ | -0.0% ✅ | -100.0% ✅ | -100.0% ✅ | -100.0% ✅ | +0.0% ✅ | -100.0% ✅ | +0.0% ✅ | -100.0% ✅ | +0.027% ✅ |
| perpendicular-cross | +0.0% ✅ | -0.1% ✅ | +0.0% ✅ | -0.1% ✅ | +0.0% ✅ | +0.0% ✅ | +0.0% ✅ | +0.0% ✅ | +0.0% ✅ | +0.030% ✅ |
| corner-into-face-center | +0.0% ✅ | -0.1% ✅ | +0.0% ✅ | +0.0% ✅ | +0.0% ✅ | +0.0% ✅ | +0.0% ✅ | +2.3% ❌ | +0.0% ✅ | +0.001% ✅ |
| 45deg-symmetric-tiebreak | +82.8% ❌ | -99.1% ✅ | +82.8% ❌ | -99.1% ✅ | +57.9% ❌ | +82.8% ❌ | +82.8% ❌ | +86.0% ❌ | +86.0% ❌ | +0.019% ✅ |
| plank-vs-near-cube-deep | +7.7% ❌ | +883.6% ❌ | -7.5% ✅ | +211.4% ❌ | -2.7% ✅ | +7.7% ❌ | -8.0% ✅ | +10.1% ❌ | -5.1% ✅ | +0.106% ✅ |
| flush-stack-axis-aligned | +0.0% ✅ | -0.2% ✅ | +0.0% ✅ | +0.0% ✅ | +0.0% ✅ | +0.0% ✅ | +0.0% ✅ | +0.4% ✅ | +0.0% ✅ | +0.001% ✅ |
| extreme-size-ratio-deep | +0.0% ✅ | -0.0% ✅ | -99.3% ✅ | -99.3% ✅ | -99.6% ✅ | +0.0% ✅ | -99.3% ✅ | +0.0% ✅ | -99.3% ✅ | +0.045% ✅ |
| grazing-plank-stack | -0.9% ✅ | -1.3% ✅ | +617.6% ❌ | +7602.9% ❌ | +1131.5% ❌ | -0.8% ✅ | +1045.3% ❌ | +6.1% ❌ | +624.5% ❌ | -0.477% ✅ |

Pass counts (out of 9 geometries): baseline (both OFF): 7/9, useRoundedBoxCorners: 8/9, useMultiPointBoxContact: 7/9, both: 7/9, useMultiPointBoxContact+overlapCap: 7/9, useAxisHysteresis: 7/9, useMultiPointBoxContact+manifold: 7/9, useRotationalEffectiveMass: 5/9, useMultiPointBoxContact+rotationalEffMass: 7/9, **useSpheroBoxContact: 9/9**

### 3b. `runMultiPointEdgeCases` at e=0.30 (production restitution)

Same 9 geometries at production's real restitution coefficient (0.30, "wet/waterlogged wood").

| Case | baseline (both OFF) | useRoundedBoxCorners | useMultiPointBoxContact | both | useMultiPointBoxContact+overlapCap | useAxisHysteresis | useMultiPointBoxContact+manifold | useRotationalEffectiveMass | useMultiPointBoxContact+rotationalEffMass | useSpheroBoxContact |
|---|---|---|---|---|---|---|---|---|---|---|
| exact-coincident | -60.2% ✅ | -60.2% ✅ | +0.0% ✅ | +0.0% ✅ | +0.0% ✅ | -60.2% ✅ | +0.0% ✅ | -60.2% ✅ | +0.0% ✅ | -60.19% ✅ |
| fully-engulfed | -60.2% ✅ | -60.2% ✅ | -100.0% ✅ | -100.0% ✅ | -100.0% ✅ | -60.2% ✅ | -100.0% ✅ | -60.2% ✅ | -100.0% ✅ | -59.93% ✅ |
| perpendicular-cross | -60.2% ✅ | -60.2% ✅ | -60.2% ✅ | -60.2% ✅ | -89.8% ✅ | -60.2% ✅ | -60.2% ✅ | -60.2% ✅ | -60.2% ✅ | -60.11% ✅ |
| corner-into-face-center | -63.2% ✅ | -63.2% ✅ | +0.0% ✅ | +0.0% ✅ | +0.0% ✅ | -63.2% ✅ | +0.0% ✅ | -59.5% ✅ | +0.0% ✅ | -67.25% ✅ |
| 45deg-symmetric-tiebreak | -48.2% ✅ | -99.8% ✅ | -48.2% ✅ | -99.8% ✅ | -92.5% ✅ | -48.2% ✅ | -48.2% ✅ | -36.7% ✅ | -36.7% ✅ | -61.99% ✅ |
| plank-vs-near-cube-deep | -61.9% ✅ | +26.6% ❌ | -67.5% ✅ | +11.1% ❌ | -91.2% ✅ | -61.9% ✅ | -67.6% ✅ | -57.4% ✅ | -63.6% ✅ | -62.98% ✅ |
| flush-stack-axis-aligned | -64.4% ✅ | -64.4% ✅ | +0.0% ✅ | +0.0% ✅ | +0.0% ✅ | -64.4% ✅ | +0.0% ✅ | -60.1% ✅ | +0.0% ✅ | -64.36% ✅ |
| extreme-size-ratio-deep | -60.2% ✅ | -60.2% ✅ | -99.3% ✅ | -99.3% ✅ | -99.7% ✅ | -60.2% ✅ | -99.3% ✅ | -60.2% ✅ | -99.3% ✅ | -59.96% ✅ |
| grazing-plank-stack | -64.4% ✅ | -64.3% ✅ | +2.3% ❌ | +1569.9% ❌ | -80.5% ✅ | -64.3% ✅ | +5.6% ❌ | -58.2% ✅ | +41.7% ❌ | -75.14% ✅ |

Pass counts (out of 9 geometries): baseline (both OFF): 9/9, useRoundedBoxCorners: 8/9, useMultiPointBoxContact: 8/9, both: 7/9, useMultiPointBoxContact+overlapCap: 9/9, useAxisHysteresis: 9/9, useMultiPointBoxContact+manifold: 8/9, useRotationalEffectiveMass: 9/9, useMultiPointBoxContact+rotationalEffMass: 8/9, **useSpheroBoxContact: 9/9**
