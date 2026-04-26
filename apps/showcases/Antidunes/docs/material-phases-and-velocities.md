# Material Phases and Velocity Overview

## The Three Phases

The simulation tracks three distinct material phases simultaneously. Each occupies cells in the lattice, and a **flag field** tells the LBM which type every cell is at any moment.

```
z
^
|  [gas]          [gas]          [gas]         <- above free surface, passive
|  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~          <- free surface (interface cells, ε ∈ [0,1])
|  [liquid]  [liquid]   [liquid]  [liquid]     <- active LBM fluid
|  [liquid]  [●particle●] [liquid] [liquid]    <- particles displace fluid cells
|  [●●●●●●●●●●●●●●●●●●●●●●●●●●●●●●●●●●●]     <- fixed particle bed (floor)
+-----------------------------------------> x (flow direction)
```

---

## Gas

Gas cells sit above the free surface and are entirely passive. The bubble model runs with constant pressure (ρ = 1) and assigns no velocity to gas. The FSLBM does not stream PDFs through gas cells; it only maintains a density boundary condition at the interface so the free surface can move freely.

**Code:** The bubble model is instantiated in `Antidunes.cpp:100` as `BubbleModelConstantPressure(1.0)`. Cell type flags (liquid / interface / gas) are registered and managed via the `FlagInfo` struct in the waLBerla free-surface module; `AntidunesBoundaryHandling::initFlagsFromFillLevel()` (`AntidunesBoundaryHandling.h:103`) sets the initial gas/liquid/interface flags from the fill level field after all solid obstacles are placed.

---

## Liquid

Liquid is the active LBM phase. Velocity and density live **implicitly in the PDFs** of each cell — recovered via the zeroth and first moments of the distribution.

### How velocity is driven

A body-force acceleration vector `acceleration[0]` acts in the x-direction on every liquid and interface cell. The force density applied per cell is:

```
liquid cell:    f = acceleration * ρ
interface cell: f = acceleration * ρ * ε      (scaled by fill level ε)
```

**Code:** `ForceDensityCodegenSweep::operator()` in `Forcing.h:52–74`. The `globalAcceleration_` shared pointer is the same object adjusted by the PID controller.

### PID flow-rate control

The acceleration is not fixed — a PID controller adjusts it every `evaluationFrequency` steps so the domain-averaged liquid velocity tracks `targetMeanVelocityMagnitude`.

- Mean velocity computation: `MeanVelocityComputer::getMeanVelocity()` (`Forcing.h:121–159`) — sums `pdfField->getVelocity(cell)` over all liquid + interface cells and MPI-reduces.
- PID update: `ForcingAdjuster::operator()` (`Forcing.h:192–209`) — calls `pidController_->update(meanVelocity->length())` on rank 0, broadcasts the new `acceleration[0]` to all ranks.
- PID implementation: `PIDController.h` / `PIDController.cpp`.

### Pressure

Initialized hydrostatically: `ρ = 1 + 3 * g * (z − H)`, giving higher density near the bed and ρ = 1 at the free surface.

**Code:** `initPdfFieldHydrostaticPressure()` in `DataInitialization.h:206–233`, using `HydrostaticDensityFtor` defined at `Utility.h:96–108`.

### Initial velocity profile

A Poiseuille profile:  
`u_x(z) = (acceleration_x / 2ν) * h * (2H − h)`  
where `h` is height above the local (sinusoidal) bed surface and `H` is the mean fluid height.

**Code:** `initPdfFieldPoiseuille()` in `DataInitialization.h:167–196`.

---

## Free Surface (Interface Cells)

Interface cells sit between liquid and gas. Each holds a **fill level** ε ∈ [0, 1] representing the fraction of the cell occupied by liquid. The FSLBM moves mass between adjacent cells each step, converting cells from interface → liquid (ε → 1) or interface → gas (ε → 0) when thresholds are crossed (`cellConversionThreshold` in `Antidunes.prm`). The sinusoidal bedform shape visible in the output is encoded in the spatial pattern of these fill levels.

**Initial profile:** A Monte-Carlo estimation samples 101×101 sub-points per cell and counts how many lie below the sinusoidal surface `cos(2πx/λ + π) * A + offset`.

**Code:** `initFillFieldSinus()` in `DataInitialization.h:117–164`. The profile function itself is `initializationProfile()` at `Utility.h:111–114`.

---

## Particles

Particles are Lagrangian spheres tracked by the DEM (mesa_pd). Each has an explicit position, linear velocity, and angular velocity — there is no per-cell particle velocity field.

### Velocity update

Velocity Verlet integration, run `numSubCycles` times (default 10) per LBM step with a fractional timestep `dt = 1 / numSubCycles`:

```
pre-force:   x(t+dt) = x(t) + v(t)*dt + 0.5*a(t)*dt²   (move positions)
force phase: compute all forces (see below)
post-force:  v(t+dt) = v(t) + 0.5*(a(t) + a(t+dt))*dt  (update velocities)
```

**Code:** `ParticleKernels::doDEMSubCycles()` in `ParticleKernels.h:206–267`. Pre-force uses `VelocityVerletPreForceUpdate`, post-force uses `VelocityVerletPostForceUpdate` (both mesa_pd kernels).

### Forces on each particle

| Force | Direction | Description |
|---|---|---|
| Hydrodynamic force | 3D | MEM: transferred from LBM boundary (see coupling section) |
| Driving body force | x, y | `volume * acceleration[x,y]` — same acceleration as the fluid |
| Effective weight | z | `volume * (densityRatio − 1) * acceleration[z]` — gravity minus buoyancy |
| Contact force | 3D | Linear spring-dashpot (particle–particle and particle–wall) |
| Lubrication correction | 3D | Squeeze-film force for near-contact pairs in viscous fluid |

**Code:** `AddExternalForces::operator()` in `ParticleKernels.h:134–146`. Note the comment on line 135 explains why buoyancy is subtracted only from z.  
Contact forces: `ComputeContacts::operator()` at `ParticleKernels.h:97–120`, using `LinearSpringDashpot`.  
Lubrication: `ComputeLubrication::operator()` at `ParticleKernels.h:59–70`, using `LubricationCorrectionKernel` with a gap-width-dependent formula.

**Why `(densityRatio − 1)` for z only:**  
The Momentum Exchange Method already transfers the hydrostatic pressure contribution to the particle in z. Applying the full gravitational weight would double-count buoyancy in z, so only the net submerged weight is added in that direction.

### Fixed vs. mobile particles

Particles with centers below `particleFixingHeight` (a multiple of the average particle diameter, set via `fixingHeightFactor` in `Antidunes.prm`) are given **infinite mass** at setup. The velocity Verlet kernel skips infinite-mass particles automatically, making them a rigid floor. All particles above that threshold are mobile and form the evolving bedform.

**Code:** The fixed-particle flag is set during `setupParticles()` in `SetupStages.h`. The threshold is `simParams->computed.particleFixingHeight`, derived in `ComputedParams` at `SimulationParams.h:222`.

---

## Phase Coupling: the MEM Handshake

This is the heart of the three-way interaction. It happens in a fixed order each timestep, explicit in the main loop (`Antidunes.cpp:137–169`):

```
① FSLBM step  (timeloopFSLBM->singleStep)
   LBM streams PDFs through the domain.
   At every fluid cell adjacent to a particle surface, the
   MovingObstacle boundary (SimpleBB) computes the momentum
   exchanged between fluid and particle — this accumulates
   the hydrodynamic force/torque on each particle.

② reduceAndAvgHydroForces  (particleKernels.reduceAndAvgHydroForces)
   MPI all-reduce sums hydrodynamic forces across ranks.
   Temporal averaging (current + previous LBM step) smooths
   the force before passing to DEM.

③ DEM sub-cycles  (particleKernels.doDEMSubCycles, ×10)
   Particles integrate positions and velocities using the
   averaged hydrodynamic force + contact + weight.
   The fluid field is NOT updated yet during these cycles.

④ resetHydroForces  (particleKernels.resetHydroForces)
   Clear hydrodynamic accumulators so the next LBM step starts fresh.

⑤ MovingParticleMapping + PDF reconstruction  (timeloopAfterDEM->singleStep)
   Now the fluid grid is updated to match new particle positions:
   - Cells newly covered by a particle  →  flagged solid, PDFs discarded
   - Cells newly vacated by a particle  →  flagged liquid/interface,
     PDFs reconstructed from surrounding fluid
```

**Velocity consistency at the particle surface:**  
The `SimpleBB` bounce-back boundary uses the **particle's own DEM velocity** to set the no-slip condition on the fluid at the particle surface. So the fluid sees the correct wall velocity of a moving particle during each LBM stream step.

**Code:**  
- `SimpleBB` (Moving Obstacle boundary): `lbm_mesapd_coupling/momentum_exchange_method/boundary/SimpleBB.h`  
- `MovingObstacle_T` registered in `AntidunesBoundaryHandling.h:74`  
- `MovingParticleMapping` and PDF reconstruction sweeps added in `BuildTimeloop.h` via `addSweepFluidParticleCoupling()`  
- Hydrodynamic force reduction and averaging: `ParticleKernels.h:186–203`  
- Force reset: `ParticleKernels.h:177–183`

---

---

## Boundary Conditions

The domain is **periodic in x and y** (streamwise and spanwise) and **bounded in z** (vertical). This is set by `periodicity < 1, 1, 0 >` in `BlockForestParameters` of `Antidunes.prm:5`.

### Domain walls (z direction)

Both the top and bottom walls use **NoSlip** — the standard full bounce-back boundary condition, which enforces zero fluid velocity at the wall surface.

```
Antidunes.prm:
  Border { direction T;  walldistance -1; NoSlip{} }   // top wall (z = max)
  Border { direction B;  walldistance -1; NoSlip{} }   // bottom wall (z = min)
```

`walldistance -1` means the boundary flag is placed in the ghost layer (one cell outside the physical domain), which is the standard waLBerla convention for solid walls.

**Code:** `setupBoundaryHandling()` in `SetupStages.h:353–366` calls `antidunesBoundaryHandling->initFromConfig()`, which reads the `BoundaryParameters` block and delegates to `BoundaryFromDomainBorder`. The `NoSlip_T` condition itself is instantiated in `AntidunesBoundaryHandling.impl.h:87`.

> **Note:** x and y NoSlip borders are present in the config but commented out (`Antidunes.prm:82–88`). They can be re-enabled to turn the periodic channel into a closed duct.

### Particle surfaces

How particle surfaces interact with the fluid depends on the `useNoSlipParticles` flag (`Antidunes.prm:51`, default `false`):

| `useNoSlipParticles` | Boundary type | Hydrodynamic forces | Typical use |
|---|---|---|---|
| `false` (default) | **MovingObstacle** (SimpleBB) | Computed via MEM | Mobile particles + fixed bed |
| `true` | Static **NoSlip** | Not computed | Fully fixed beds only |

With `useNoSlipParticles = false` (the normal case), every cell overlapping a particle is flagged `MovingObstacle`. The `SimpleBB` boundary imposes the **particle's own DEM velocity** as the no-slip condition at the surface during each LBM stream step — this is both a wall BC and the MEM force measurement in one.

For the fixed bed particles (below `particleFixingHeight`), this is also `MovingObstacle`, but their DEM velocity is always zero, so the effect is equivalent to a static no-slip wall.

**Code:**
- `mapParticlesOntoCells()` in `SetupStages.h:368–393` — chooses between `ParticleMappingKernel` (NoSlip path) and `MovingParticleMappingKernel` (MovingObstacle path).
- `MovingObstacle_T` defined as `lbm_mesapd_coupling::SimpleBB<...>` in `AntidunesBoundaryHandling.h:74`.
- Safety check: if `useNoSlipParticles = true` but not all particles are fixed, the simulation aborts (`SetupStages.h:301–303`).

### DEM plane walls

In addition to the LBM no-slip walls, two **infinite-mass DEM planes** are created at the bottom and top of the domain for particle–wall contact detection in the DEM:

```cpp
createPlane(*simObjs->particleStorage, bf->getDomain().minCorner(), Vector3<real_t>(0, 0,  1)); // floor
createPlane(*simObjs->particleStorage, bf->getDomain().maxCorner(), Vector3<real_t>(0, 0, -1)); // ceiling
```

These are mesa_pd `Plane` particles with infinite mass (never moved by force). They participate in the same `LinearSpringDashpot` contact detection as particle–particle contacts, preventing particles from leaving the domain in z.

**Code:** `SetupStages.h:325–326`.

### Free surface (gas–liquid interface)

The gas–liquid interface is not a hard wall but a **pressure boundary** managed by the FSLBM. At interface cells adjacent to gas, `SimplePressureWithFreeSurface` applies a constant reference density ρ = 1 (atmospheric pressure), derived from `BubbleModelConstantPressure(1.0)`.

**Code:** `Antidunes.cpp:100` (bubble model creation); `AntidunesBoundaryHandling.impl.h:88–97` (pressure BC instantiation); the `Pressure_T` alias is `walberla::free_surface::SimplePressureWithFreeSurface` (`AntidunesBoundaryHandling.h:69`).

### All registered boundary types

`AntidunesBoundaryHandling` registers all of the following types at startup (even if not all are used in the default config):

| Type | Flag ID | Used by default | Purpose |
|---|---|---|---|
| `NoSlip` | `NoSlip` | Yes (walls) | Zero-velocity wall |
| `MovingObstacle` (SimpleBB) | `MovingObstacle` | Yes (particles) | Particle surface + MEM |
| `SimplePressureWithFreeSurface` | `Pressure` | Yes (free surface) | Atmospheric pressure at interface |
| `FreeSlip` | `FreeSlip` | No (commented out) | Symmetry / slip wall |
| `UBB` | `UBB` | No | Prescribed velocity wall |
| `UBB_Inflow` | `UBB_Inflow` | No | Inflow with interface cell generation |
| `Pressure` (outflow) | `PressureOutflow` | No | Pressure outlet |
| `Outlet` | `Outlet` | No | Convective outlet |

**Code:** Full boundary handling type definition at `AntidunesBoundaryHandling.h:77–80`; all flag UIDs defined as static constants in `AntidunesBoundaryHandling.impl.h:197–275`.

---

## Summary: velocity ownership per phase

| Phase | Velocity representation | Updated by |
|---|---|---|
| Gas | None | — |
| Liquid / Interface | PDFs in `pdfField` (implicit) | LBM collision + stream kernel |
| Free surface position | Fill level `ε` in `fillField` | FSLBM mass transfer sweeps |
| Particles (mobile) | Explicit `velocity` per particle in `ParticleStorage` | DEM velocity Verlet |
| Particles (fixed bed) | Zero (infinite mass) | Never updated |
