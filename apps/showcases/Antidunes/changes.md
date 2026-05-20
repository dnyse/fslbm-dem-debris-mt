# Antidunes — Change Log

---

## Meeting 02

### Goals
- Fix a NaN quaternion crash when running with a floating sphere obstacle
- Replace box pillar with a cylindrical pillar supported in both DEM and LBM

---

### 1. `ParticleKernels.h` — NaN contact normal guard

**Added NaN check in `ComputeContacts::operator()`** before applying the contact response:

```cpp
const auto& n = acd_.getContactNormal();
if (std::isnan(n[0]) || std::isnan(n[1]) || std::isnan(n[2]))
{
   WALBERLA_LOG_WARNING_ON_ROOT("[ComputeContacts] NaN contact normal detected! ...");
   return; // skip this contact to avoid NaN propagation
}
```

**Why:** Two overlapping particles at exactly the same position produce a zero-length contact normal that collapses to NaN after normalization inside `detectSphereSphereCollision`. Passing NaN forces into the velocity Verlet integrator propagates to the quaternion update (`Quaternion.h:241`), which then aborts. The guard logs the offending pair and skips the contact instead.

---

### 2. `ParticleKernels.h` — `SolidCylinderContactDetection`

**Added a subclass of `AnalyticContactDetection`** that overrides `(Sphere, CylindricalBoundary)` contact detection with exterior (solid cylinder) math:

```cpp
class SolidCylinderContactDetection : public mesa_pd::collision_detection::AnalyticContactDetection
{
public:
   using AnalyticContactDetection::operator(); // prevent name hiding of all base overloads

   template< typename Accessor >
   bool operator()(const size_t idx1, const size_t idx2,
                   const mesa_pd::data::Sphere& s,
                   const mesa_pd::data::CylindricalBoundary& cb,
                   Accessor& ac)
   {
      // Radial distance from cylinder axis; contact fires when sphere is too close from outside
      const mesa_pd::Vec3 delta  = posSphere - posCyl;
      const mesa_pd::Vec3 radial = delta - (delta * axis) * axis;
      const real_t dist          = radial.length();
      const real_t separation    = dist - cb.getRadius() - s.getRadius();
      if (separation < getContactThreshold()) { ... }
   }
   // symmetric overload delegates to above
};
```

**Why `AnalyticContactDetection` is subclassed, not `CylindricalBoundary`:** `DoubleCast` dispatches by `SHAPE_TYPE` integer and `static_cast`s to the base type — a subclass with the same `SHAPE_TYPE` gets sliced and always invokes the base-class (hollow/interior) detection. Subclassing the *functor* instead works because `DoubleCast` deduces the concrete functor type, so C++ overload resolution picks the subclass method.

**Why `using AnalyticContactDetection::operator();` is required:** Defining any `operator()` in the subclass hides *all* base-class `operator()` overloads by C++ name hiding rules. Without the `using` declaration, `DoubleCast`'s calls for all other shape pairs (Sphere-Sphere, Sphere-Box, Ellipsoid-*, etc.) would fail to compile.

**`ComputeContacts::acd_`** changed from `AnalyticContactDetection` to `SolidCylinderContactDetection`. All existing contact pairs (sphere-sphere, sphere-box, sphere-halfspace) are unaffected — they fall through to the inherited base-class methods.

---

### 3. `Utility.h` — `createCylinder()`

**Added `createCylinder()` helper** alongside `createBox()`, using `CylindricalBoundary` as the shape:

```cpp
inline auto createCylinder(mesa_pd::data::ParticleStorage& ps, const mesa_pd::Vec3& center,
                           real_t radius, const mesa_pd::Vec3& axis = mesa_pd::Vec3(0, 0, 1))
```

Sets the same GLOBAL / INFINITE / FIXED / NON_COMMUNICATING flags as `createBox`. Type is 2. Interaction radius is `infinity()` (same as `createBox`).

**Why `CylindricalBoundary` for a solid cylinder:** It stores exactly the data needed (radius + axis) and its `SHAPE_TYPE` is already handled by `DoubleCast`. The hollow/interior contact math it normally triggers is replaced by `SolidCylinderContactDetection`. End caps are handled by the existing floor/ceiling `HalfSpace` planes — no finite-length cylinder shape is needed.

---

### 4. `SimulationParams.h` — `CylinderPillarParams`

**Added `CylinderPillarParams`** (center + radius) and **`parseCylinderPillars()`** alongside the existing `PillarParams` / `parsePillars()`. `SimulationParams` now holds:

```cpp
std::vector< PillarParams >         pillars;         // box pillars (existing)
std::vector< CylinderPillarParams > cylinderPillars; // cylinder pillars (new)
```

Config block key: `CylinderPillars` → `CylinderPillar { center <x,y,z>; radius r; }`. Zero entries is valid.

---

### 5. `SetupStages.h` — cylinder pillar creation

**`createParticleFct` lambda** extended with a radial exclusion check for cylinder pillars:

```cpp
for (const auto& pillar : simParams->cylinderPillars)
{
   const real_t dx = pos[0] - pillar.center[0];
   const real_t dy = pos[1] - pillar.center[1];
   if (std::sqrt(dx * dx + dy * dy) < pillar.radius + maxSphereRadius)
      return false;
}
```

**DEM creation loop** added after the box pillar loop:

```cpp
for (const auto& pillar : simParams->cylinderPillars)
   createCylinder(*simObjs->particleStorage, pillar.center, pillar.radius);
```

---

### 6. `Antidunes.prm` — cylinder pillar config

**DEM side:** Added `CylinderPillars` block (active) and kept `Pillars` block (commented out):

```
CylinderPillars
{
   CylinderPillar
   {
      center < 64, 30, 60 >;
      radius 10;
   }
}
```

**LBM side:** Added `BodyBC` initializer in `BoundaryParameters` (the `Body` key maps to `OverlapFieldFromBody` for fill-level init; `BodyBC` maps to `BoundaryFromBody` for LBM boundary flags):

```
BodyBC { shape Cylinder; min < 64, 30, 0 >; max < 64, 30, 120 >; radius 10; NoSlip{} }
```

`min`/`max` are the cylinder axis endpoints (bottom and top). End caps are the existing floor/ceiling `Border` NoSlip lines.

---

### 7. Floating sphere — stability fix

The floating sphere (`FloatingSphereParams`, `AddExternalForces` with type-1 check, optional `FloatingSphere` config block) was already implemented. A NaN quaternion crash was diagnosed and fixed by correcting the position and radius in `Antidunes.prm`:

- **Root cause:** Sphere placed at `z=60` with `radius=10` — the center sat exactly at the free surface (`absoluteLiquidHeight ≈ 60 LU`). MEM computes momentum exchange only over fluid-covered cells; with the sphere half in gas the torque is large, asymmetric, and accumulates until angular velocity diverges.
- **Fix:** Reduce radius to 4 LU and position fully below the free surface (`z=48`, top of sphere at `z=52`). The sphere is currently commented out in the `.prm` and can be re-enabled.

### Goal
Fix a SIGSEGV crash (fault address `0x18`) that occurred on the first call to `timeloopFSLBM->singleStep()` at `beginTimeStep = 0`.

### Root cause (diagnosed)
Progressive debug logging narrowed the crash to the box pillar VTK output added in Meeting 01. The execution order is:
1. fluid field VTK ✓ (logged "done")
2. sphere particle VTK ✓ (logged "done")
3. box pillar VTK ✗ → SIGSEGV at address `0x18` before "done" was logged

Address `0x18` = offset 24 from a null pointer = `inertiaBF_` in `BaseShape`. The crash occurred inside `vtk::writeFiles(boxVTKOutput)` — the exact internal call site could not be pinpointed without a debug build with symbols on the cluster (GCC 8.5.0, no `-g` flag in the release build).

### Fix

**`BuildTimeloop.h`** — removed the entire box pillar VTK block (`SelectBoxEdgeLength`, `boxPointData`, `boxVTKOutput`, and the before-step registration). The pillar appears as a zero-velocity void in every fluid field VTK slice, and its geometry is fully specified in `Antidunes.prm`. Also removed the `#include "mesa_pd/data/shape/Box.h"` that was only needed for this block.

**`SetupStages.h`** — removed the debug `firstCall` shared_ptr and associated per-rank log from the `createParticleFct` lambda.

**`Antidunes.cpp`** — removed the post-VTK debug marker `addFuncBeforeTimeStep` that was added to identify the crash phase.

---

### Goal
Make pillar geometry configurable at runtime (no recompile) and support an arbitrary number of pillars.

---

### Changed Files

#### 1. `SimulationParams.h`

**Added `parsePillars()` helper:**
```cpp
inline std::vector< PillarParams > parsePillars(const shared_ptr< config::Config >& cfg)
{
   std::vector< PillarParams > result;
   config::Config::Blocks blocks;
   cfg->getOneBlock("Pillars").getBlocks("Pillar", blocks, 0, SIZE_MAX);
   for (const auto& b : blocks)
      result.emplace_back(b);
   return result;
}
```

**Changed `PillarParams pillar` → `std::vector<PillarParams> pillars`** in `SimulationParams`, initialized via `parsePillars(cfg)`.

---

#### 2. `SetupStages.h` — `setupParticles()`

**`createParticleFct` lambda** now loops over all pillars:
```cpp
for (const auto& pillar : simParams->pillars)
{
   const real_t dx = std::max(real_t(0), std::abs(pos[0] - pillar.center[0]) - pillar.halfExtents[0]);
   const real_t dy = std::max(real_t(0), std::abs(pos[1] - pillar.center[1]) - pillar.halfExtents[1]);
   if (std::sqrt(dx * dx + dy * dy) < maxSphereRadius)
      return false;
}
```

**`createBox()` call** replaced by a loop:
```cpp
for (const auto& pillar : simParams->pillars)
{
   auto p = createBox(*simObjs->particleStorage, pillar.center, real_t(2) * pillar.halfExtents);
   p->setType(2);
}
```

---

#### 3. `Antidunes.prm`

Replaced the single `PillarParameters` block with a `Pillars` block containing one or more `Pillar` sub-blocks:

```
Pillars
{
   Pillar
   {
      center      < 200, 30, 60 >;
      halfExtents < 16, 10, 60 >;
   }
}
```

Zero `Pillar` blocks is valid — the simulation runs without any pillar. For a box pillar, a matching `CellInterval` is also needed in `BoundaryParameters`:
- `CellInterval min = center - halfExtents`
- `CellInterval max = center + halfExtents - 1`

---

## Meeting 01

This section describes all changes made to add a fixed solid box pillar to the Antidunes simulation.

---

## Overview

A fixed rectangular pillar (`Box` shape in MESA-PD) was added to the middle of the simulation domain, geometrically matching the `CellInterval` already defined in `BoundaryParameters` of `Antidunes.prm`. The pillar acts as a solid obstacle for both the DEM particle bed and the LBM fluid.

**Pillar geometry (LBM lattice units):**
- CellInterval (from `Antidunes.prm`): `min <184, 20, 0>`, `max <215, 39, 119>`
- Box center: `(200, 30, 60)`
- Box half-extents: `(16, 10, 60)`
- Box edge lengths: `(32, 20, 120)`

---

## Changed Files

### 1. `Utility.h`

**Added `FixedBox` subclass** (before `createBox()`):
```cpp
class FixedBox : public mesa_pd::data::Box {
public:
   explicit FixedBox(const mesa_pd::Vec3& edgeLength) : Box(edgeLength) {}
   void updateMassAndInertia(const real_t /*density*/) override {
      mass_         = std::numeric_limits<real_t>::infinity();
      invMass_      = real_t(0);
      inertiaBF_    = mesa_pd::Mat3(std::numeric_limits<real_t>::infinity());
      invInertiaBF_ = mesa_pd::Mat3(real_t(0));
   }
};
```

**Why:** The contact force stiffness is computed as `k_n = m_eff * (π² + ln²(e)) / t_c²` where
`m_eff = 1 / (invMass_sphere + invMass_obstacle)`. For a wall (`HalfSpace`), `invMass = 0` so
`m_eff = mass_sphere`. `FixedBox` overrides `updateMassAndInertia` to set `invMass = 0`, making
the contact response with the pillar identical to wall contacts.

**Added `createBox()` function:**
```cpp
inline auto createBox(mesa_pd::data::ParticleStorage& ps,
                      const mesa_pd::Vec3& center,
                      const mesa_pd::Vec3& edgeLength)
```

Sets GLOBAL / INFINITE / FIXED / NON_COMMUNICATING flags. Type defaults to 0 and is overridden to 2 in `setupParticles`. Interaction radius is `infinity()`.

**Why `INFINITE` flag:** `ExcludeInfiniteInfinite` in the contact loop excludes pairs where both particles are `INFINITE`. Without it, Box↔HalfSpace pairs would reach `AnalyticContactDetection`, which has no overload for that combination and aborts. With `INFINITE` set:
- Box↔HalfSpace → both INFINITE → excluded ✓
- Box↔Sphere → Sphere not INFINITE → still included ✓

---

### 2. `SetupStages.h` — `setupParticles()`

**Modified `createParticleFct`** to exclude spheres overlapping the pillar column using nearest-point distance in XY (Z is omitted because the pillar spans the full domain height).

**Added `createBox()` call** after `createPlane()` calls.

---

### 3. `Antidunes.prm`

Added `CellInterval` in `BoundaryParameters` for the LBM NoSlip boundary matching the box geometry.

---

## What Was NOT Changed (and Why)

### Particle mapping (LBM–DEM coupling)
Both `mapParticlesOntoCells` and `MovingParticleMapping` use a `sphereSelector` filtering by `Sphere::SHAPE_TYPE`. Box and cylinder pillars are excluded from MEM mapping — they never conflict with the NoSlip flags set by the boundary handling.

### Fill level inside pillar
`initFlagsFromFillLevels` explicitly skips cells where `handling->isBoundary()` is true. NoSlip cells keep their flag regardless of fill level.

### Zero velocity in pillar cells
`ZeroSetter` already zeros PDFs in obstacle cells before every VTK write.

---

## Verified Contact Detection Correctness

The sphere–box contact normal always points **outward** via `detectSphereBoxCollision` in `AnalyticCollisionFunctions.h`. `DoubleCast` routes `Sphere ↔ Box` pairs to `AnalyticContactDetection`, which calls `detectSphereBoxCollision`. `LinearSpringDashpot` applies the contact force using the same parameter set (type 2 vs 0, same friction/restitution/kappa) as sphere–wall contacts.