# Antidunes — Box Pillar Obstacle Changes

This document describes all changes made to add a fixed solid box pillar to the Antidunes simulation and the reference code each change is based on.

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

**Added includes:**
```cpp
#include "mesa_pd/data/shape/Box.h"
#include "mesa_pd/data/shape/HalfSpace.h"
```
`HalfSpace.h` was previously only transitively included and is needed by `createPlane()`.

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
`m_eff = mass_sphere`. A plain `Box` with finite volume would give a slightly different `m_eff`.
`FixedBox` overrides `updateMassAndInertia` to set `invMass = 0` (same as `HalfSpace`), making
the contact response with the pillar identical to the contact response with the floor/ceiling walls.

**Reference:** `src/mesa_pd/data/shape/HalfSpace.h` — `HalfSpace::updateMassAndInertia` (sets
`mass_ = infinity`, `invMass_ = 0`). The `protected` members `mass_`, `invMass_`, `inertiaBF_`,
`invInertiaBF_` are defined in `src/mesa_pd/data/shape/BaseShape.h`.

---

**Added `createBox()` function:**
```cpp
inline auto createBox(mesa_pd::data::ParticleStorage& ps,
                      const mesa_pd::Vec3& center,
                      const mesa_pd::Vec3& edgeLength)
{
   auto p0 = ps.create(true);
   p0->setPosition(center);
   p0->setBaseShape(std::make_shared< FixedBox >(edgeLength));
   p0->getBaseShapeRef()->updateMassAndInertia(real_t(1));
   p0->setOwner(walberla::mpi::MPIManager::instance()->rank());
   p0->setType(0);
   p0->setInteractionRadius(real_t(0.5) * edgeLength.length());
   mesa_pd::data::particle_flags::set(p0->getFlagsRef(), mesa_pd::data::particle_flags::GLOBAL);
   mesa_pd::data::particle_flags::set(p0->getFlagsRef(), mesa_pd::data::particle_flags::FIXED);
   mesa_pd::data::particle_flags::set(p0->getFlagsRef(), mesa_pd::data::particle_flags::INFINITE);
   mesa_pd::data::particle_flags::set(p0->getFlagsRef(), mesa_pd::data::particle_flags::NON_COMMUNICATING);
   return p0;
}
```

**Key differences from `createPlane()`:**
- Uses `FixedBox` instead of `HalfSpace`
- `interactionRadius` = half-diagonal of the box (bounding sphere radius), not infinity
- Has `INFINITE` flag — required to prevent Box↔HalfSpace contact detection (see below)

**Reference:** `createPlane()` in `Utility.h` — same flag pattern (`GLOBAL`, `INFINITE`, `FIXED`,
`NON_COMMUNICATING`). The `INFINITE` flag is required: `ExcludeInfiniteInfinite`
in the contact loop (`ParticleKernels.h:237`) excludes pairs where **both** particles are `INFINITE`.
Without it, Box↔HalfSpace pairs would be passed to `AnalyticContactDetection`, which has no
overload for that combination and aborts at the generic fallback (line 166). With `INFINITE` set:
- Box↔HalfSpace → both INFINITE → excluded ✓
- Box↔Sphere → Sphere not INFINITE → still included ✓

---

### 2. `SetupStages.h` — `setupParticles()`

**Added box geometry constants** (before the `if (bfParams.loadSnapshot)` block):
```cpp
const Vector3< real_t > boxCenter(real_t(200), real_t(30), real_t(60));
const Vector3< real_t > boxHalfExtents(real_t(16), real_t(10), real_t(60));
```

Defined once and used in both the particle exclusion filter and the `createBox()` call, so the
geometry is specified in a single place.

---

**Modified `createParticleFct`** to exclude spheres overlapping the pillar column:
```cpp
const real_t maxSphereRadius = real_t(0.5) * simParams->computed.avgParticleDiameter * real_t(1.1);
auto createParticleFct = [&simParams, &boxCenter, &boxHalfExtents, maxSphereRadius](Vector3<real_t> pos) {
   if (pos[2] >= initializationProfile(...)) return false;
   // pillar spans full z height; use nearest-point distance in x-y
   const real_t dx = std::max(real_t(0), std::abs(pos[0] - boxCenter[0]) - boxHalfExtents[0]);
   const real_t dy = std::max(real_t(0), std::abs(pos[1] - boxCenter[1]) - boxHalfExtents[1]);
   return std::sqrt(dx * dx + dy * dy) >= maxSphereRadius;
};
```

**Why:** `initSpheresFromFile` calls `createParticleFct(pos)` with the sphere center before
radius is known. Without this filter, spheres from the bed file that fall inside (or touching)
the pillar volume would be created, causing their centers to be inside the box. The
`detectSphereBoxCollision` function has a special `!outside` fallback (line 96–121 in
`AnalyticCollisionFunctions.h`) for sphere centers inside the box — avoiding this entirely
is cleaner.

The exclusion check uses the same nearest-point distance formula as `detectSphereBoxCollision`
itself: clamp sphere center to box bounds, measure distance. The z-axis check is omitted because
the pillar spans the full domain height (z half-extent = 60 = domain half-height), making that
term always zero.

**Reference:** `src/mesa_pd/collision_detection/AnalyticCollisionFunctions.h` —
`detectSphereBoxCollision` lines 56–142. The `margin` uses the `1.1 × avgParticleDiameter`
factor already used at `SetupStages.h:336` for lubrication correction safety.

---

**Added `createBox()` call** (after `createPlane()` calls):
```cpp
createBox(*simObjs->particleStorage, boxCenter, real_t(2) * boxHalfExtents);
```

**Reference:** `createPlane()` calls directly above at `SetupStages.h:331–332`.

---

### 3. `BuildTimeloop.h` — `addBeforeTimeStepVTKOutput()`

**Added include:**
```cpp
#include "mesa_pd/data/shape/Box.h"
```

**Added custom VTK output selector and writer for the pillar:**
```cpp
struct SelectBoxEdgeLength {
   using return_type = mesa_pd::Vec3;
   mesa_pd::Vec3 operator()(const mesa_pd::data::ParticleStorage::iterator& pIt) const {
      return static_cast< mesa_pd::data::Box* >(pIt->getBaseShape().get())->getEdgeLength();
   }
};
auto boxPointData = std::make_shared< mesa_pd::vtk::ParticleVtkOutput >(simObjs->particleStorage);
boxPointData->addOutput< SelectBoxEdgeLength >("edgeLength");
boxPointData->setParticleSelector([](const mesa_pd::data::ParticleStorage::iterator& pIt) {
   return pIt->getBaseShape()->getShapeType() == mesa_pd::data::Box::SHAPE_TYPE;
});
auto boxVTKOutput = vtk::createVTKOutput_PointData(boxPointData, "pillar",
   simParams->vtk.writeFrequency, simParams->vtk.baseFolder, ...);
timeloop->addFuncBeforeTimeStep(vtk::writeFiles(boxVTKOutput), "VTK: pillar data");
```

Output file: `vtk-out/pillar_*.vtu` — one point at the box center with attribute `edgeLength = (32, 20, 120)`.

**ParaView visualization:** Load `pillar_*.vtu` → Filters → Glyph → glyph type **Box** →
Scale Array: `edgeLength`, Scale Mode: `vector`. This renders the pillar as a correctly-sized
solid box at the correct position.

**Reference:** Existing sphere particle VTK output at `BuildTimeloop.h:109–122`. The
`SelectBoxEdgeLength` selector follows the `OutputSelector` pattern in
`src/mesa_pd/vtk/OutputSelector.h` and the `SelectParticleInteractionRadius` pattern in
`src/mesa_pd/data/ParticleStorage.h:1480`.

---

## What Was NOT Changed (and Why)

### LBM boundary condition
The `CellInterval` in `Antidunes.prm` already defines the pillar as a `NoSlip` boundary for the
fluid. `BoundaryFromCellInterval` (registered in `AntidunesBoundaryHandling.impl.h:295`) handles
this at startup via `initFromConfig`.

### Zero velocity in pillar cells
`ZeroSetter` in `src/lbm/free_surface/VtkWriter.h:186` already zeros PDFs in cells where
`flagInfo_.isGas() || flagInfo_.isObstacle()` before every VTK write. NoSlip cells are obstacle
cells → velocity is always shown as zero in VTK output. For the simulation physics, the NoSlip
bounce-back operates on adjacent fluid cell PDFs and never reads the obstacle cell's own PDFs.

### Particle mapping (LBM–DEM coupling)
Both `mapParticlesOntoCells` (`SetupStages.h:390`) and `MovingParticleMapping`
(`BuildTimeloop.h:359`) use a `sphereSelector` that filters by `Sphere::SHAPE_TYPE`. The
`FixedBox` (`Box::SHAPE_TYPE = 3`) is excluded from mapping, so it never conflicts with the
`CellInterval` NoSlip flags already set by the boundary handling.

### Fill level inside pillar
`initFlagsFromFillLevels` in `src/lbm/free_surface/InitFunctions.h:86` explicitly skips cells
where `handling->isBoundary()` is true. NoSlip cells are boundary cells → they keep their NoSlip
flag regardless of fill level. `initFillLevelsInBoundaries` additionally sets fill_level = 1 in
all obstacle cells to prevent the bubble model from treating them as gas.

---

## Verified Contact Detection Correctness

The sphere–box contact normal always points **outward** (from box surface toward sphere center)
for spheres approaching from outside, via `detectSphereBoxCollision` in
`src/mesa_pd/collision_detection/AnalyticCollisionFunctions.h`:

```cpp
const Vec3 q( R * p );  // nearest point on box surface (world frame, relative to box center)
const Vec3 n( d - q );  // sphere_center - nearest_surface_point  →  always points outward
```

The `DoubleCast` kernel (`src/mesa_pd/kernel/DoubleCast.h:88–92`) routes `Sphere ↔ Box` pairs
to `AnalyticContactDetection`, which calls `detectSphereBoxCollision`. `LinearSpringDashpot`
then applies the contact force using the computed normal — identical parameter set (type 0 vs 0,
same friction/restitution/kappa) as sphere–wall contacts.