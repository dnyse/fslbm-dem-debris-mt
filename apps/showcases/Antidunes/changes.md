# Antidunes — Change Log

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

### Buoyancy uses submerged volume (FIXED)

**Problem (original):** `AddExternalForces::operator()` applied `ρ_f · V_full · g` for all particles,
overestimating buoyancy for partially submerged floating spheres (type 1).

**Fix (active in `ParticleKernels.h`):** For floating spheres, the precomputed submerged volume
`V_sub` from `computeSubmergedVolumes()` (spherical-cap integration over the fill-level field)
is used instead of the full sphere volume. Non-floaters are always fully submerged so their
behaviour is unchanged (`V_sub = V_full`).

```cpp
// X/Y: full-volume driving force (pressure-gradient body force acts on the whole particle).
// Z:   gravity (rho * vFull) minus correct buoyancy (vSub × g, fluid density = 1).
Vector3< real_t > bodyForce(
   vFull * (*acceleration_)[0],
   vFull * (*acceleration_)[1],
   (rho * vFull - vSub) * (*acceleration_)[2]);
```

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

## Active source changes (waLBerla library files)

### `src/lbm/free_surface/surface_geometry/CurvatureModel.impl.h`

`ObstacleNormalSweep` and `ObstacleFillLevelSweep` are now **unconditional** in
`FiniteDifferenceMethod::addSweeps` (previously gated on `enableWetting_`). Also
`NormalSweep` always uses `includeObstacleNeighbors=true`.

This affects other showcases that use the library class. Antidunes uses `BuildTimeloop.h`
directly and manages these sweeps independently (currently commented out — see note below).

### `src/lbm/free_surface/surface_geometry/ObstacleFillLevelSweep.impl.h`

Removed the `WALBERLA_CHECK_GREATER(normal.length(), 0, ...)` assertion that aborted when an
obstacle cell's neighbourhood is exactly symmetric (e.g., sphere equator with liquid below,
gas above, producing a zero net normal). Replaced with an equal-weight fallback:

```cpp
const real_t weight = obstNormalSqrLen > 0
   ? std::abs(obstNormal * dirVector_normalized)  // normal-weighted (original)
   : real_c(1);                                   // equal-weight fallback when normal is zero
```

### `src/lbm/free_surface/surface_geometry/ObstacleNormalSweep.impl.h`

Added `WALBERLA_LOG_DEVEL_ON_ROOT` when an obstacle cell's averaged direction vector is zero
(triggers the equal-weight fallback in `ObstacleFillLevelSweep`). Visible with `WALBERLA_LOGLEVEL=DEVEL`.

### `src/lbm/free_surface/dynamics/ExcessMassDistributionSweep.h/.impl.h`

Added `newGas` flag to the isolated-cell warning in `distributeMassInterfaceAndLiquid` so the
log message reports whether the orphaned cell converted to gas or liquid.

### `src/lbm_mesapd_coupling/momentum_exchange_method/boundary/SimpleBB.h`

Added finite-check on `pdf_old` before the bounce-back computation:

```cpp
if (!std::isfinite(pdf_old))
   WALBERLA_LOG_WARNING("[SimpleBB] NaN/Inf pdf_old at fluid cell (" << x << "," << y << "," << z
      << ") -> obstacle (" << nx << "," << ny << "," << nz << ") dir=" << dir << " pdf=" << pdf_old);
```

`SimpleBB` is a victim, not the NaN origin — this diagnostic identifies which fluid cell was
already corrupt when the bounce-back tried to read it.

### `src/lbm_mesapd_coupling/momentum_exchange_method/reconstruction/PdfReconstructionManager.h`

Added finite-check after `(*reconstructor_)(...)` to detect if reconstruction itself produced NaN
(distinguishes fill-level collapse from reconstruction failure as two different failure modes).

---

## Note on `ObstacleNormalSweep` / `ObstacleFillLevelSweep` in Antidunes

These sweeps are currently **commented out** in `BuildTimeloop.h`, and `NormalSweep` has
`includeObstacleNeighbors=false` to match. The two are a coupled pair:
`ObstacleFillLevelSweep` writes fictitious fill values into obstacle cells of `smoothFillField`;
`NormalSweep` reads them when `includeObstacleNeighbors=true`. Without `ObstacleFillLevelSweep`
running, obstacle cells have fill=0 (indistinguishable from gas), so enabling
`includeObstacleNeighbors` would corrupt normals near pillars and sphere.

Re-enable both together if wetting / contact-angle effects at obstacle surfaces are needed.

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
Both `mapParticlesOntoCells` and `MovingParticleMapping` use a `sphereSelector` filtering by
`Sphere::SHAPE_TYPE`. Box and cylinder pillars are excluded from MEM mapping — they never
conflict with the NoSlip flags set by the boundary handling.

### Fill level inside pillar
`initFlagsFromFillLevels` explicitly skips cells where `handling->isBoundary()` is true.
NoSlip cells keep their flag regardless of fill level.

### Zero velocity in pillar cells
`ZeroSetter` already zeros PDFs in obstacle cells before every VTK write.
