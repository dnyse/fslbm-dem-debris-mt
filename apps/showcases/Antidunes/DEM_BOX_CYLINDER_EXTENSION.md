# DEM extension for box and cylinder floaters

Summary of what had to change in the DEM to support floating boxes and cylinders, on top of a
solver originally built around rotationally-symmetric spheres.

## What DEM does, in three pieces

1. **Equations of motion** — Newton–Euler, integrated with `mesa_pd::kernel::VelocityVerlet`,
   sub-cycled at a much finer `Δt_p ≪ Δt` than the fluid solver.
2. **Contact detection + force law** — for any touching pair, obtain a contact point / normal /
   penetration depth, then apply the same spring-dashpot + Coulomb friction law regardless of how
   those three numbers were produced.
3. **Lattice mapping** — which fluid cells belong to a given body, for the momentum-exchange
   coupling.

## Why spheres were the easy case

A sphere only ever needs a **scalar** moment of inertia (`I = 2/5 m R²`), and its footprint,
lattice mapping, and contact geometry against anything are all independent of how it's oriented.
Rotation is tracked (friction imparts spin) but nothing downstream cares what that rotation *is*.

## What changes for box and cylinder

The integrator and force law stayed generic and untouched — `mesa_pd`'s rigid-body machinery
(quaternion-based orientation `Rot3`, the `ContainsPointFunctor` point-in-shape dispatch) already
supported arbitrary rigid bodies before this work. What had to be built or adapted was everything
*downstream* that used to get away with a sphere's rotational symmetry:

- **Anisotropic inertia tensor.** Box: `I = m/12 · diag(b²+c², a²+c², a²+b²)`. Cylinder:
  `I_axis = ½ m r²`, `I_⊥ = m(3r²+h²)/12`. Since these principal values differ, orientation now
  affects how torque converts to angular acceleration.

- **Orientation-aware lattice mapping.** Cell membership now requires transforming the cell into
  the body frame (via `Rᵀ`) before testing against the box's half-extents or the cylinder's
  radius/half-length, instead of a single distance check.

- **New contact-detection overloads** (`ParticleKernels.h`, `SolidCylinderContactDetection`):
  - Sphere–cylinder: radial + axial clamp in the cylinder's body frame (closed form).
  - Box–wall, box–pier, cylinder–pier: rotation-aware — transform the wall/pier geometry into the
    floater's body frame, find the closest point there, transform back. The pier cases need a
    small fixed-point iteration (2–3 steps) since a pier is an infinite line, not a point or plane.
  - Box–box, box–cylinder, cylinder–cylinder: no closed form exists, so these fall back to the
    generic GJK/EPA path, which works off any shape's support function.

- **Numerical safety.** A sphere touching a thin pillar can never "swallow" it. A large box or log
  rotating near a pillar can end up with the pillar's axis passing through its cross-section —
  geometrically valid but numerically degenerate (near-zero contact normal → unbounded spring
  force → NaN in the coupled fluid). The box/cylinder–pier overloads clamp the reported
  penetration depth to a bounded maximum overlap to guard against this; sphere contacts never
  needed it.

- **Buoyancy no longer falls out for free.** For a small sphere, Archimedes' force emerges
  implicitly from the hydrostatic pressure field the MEM already sees. For a large, slowly-tilting
  box, the submerged volume and centre of buoyancy have to be computed explicitly against the
  locally fitted free-surface plane — exactly (tetrahedral clip) for boxes, by sampling for
  cylinders (no closed form for an obliquely cut cylinder) — because orientation now directly
  determines how much volume is submerged and where the buoyant force acts, producing a righting
  or capsizing torque a sphere could never generate.

- **New configuration knobs that are meaningless for a sphere**: `initialTiltAngleDeg` (release
  away from equilibrium orientation to watch it right itself), `comOffsetZ` (tune metacentric
  height, i.e. ballast), `pinPosition`/`pinAxes` (isolate translational vs. rotational behaviour).

## In short

The integrator and force law stayed generic and untouched; what this thesis actually added was
everything that has to know a rigid body's orientation — contact geometry against flat walls and
cylindrical piers, the submerged-volume/buoyancy calculation, and the numerical guards that keep
that orientation-dependence from blowing up near a thin obstacle.
