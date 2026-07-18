# Independent cross-checks for ContactDetectionTest.cpp

Three independent validations of `SolidCylinderContactDetection` (`ParticleKernels.h`), none of
which share code with the C++ under test or with each other:

| Script | What it is | Coverage |
|---|---|---|
| `validate_pybullet.py` | Bullet Physics, a separately-developed rigid-body engine | All 12 shape pairs, 27 checks |
| `validate_lammps_granular.py` | LAMMPS's GRANULAR package, a real DEM/granular contact-mechanics implementation | Sphere x Sphere, Sphere x HalfSpace only (LAMMPS's granular package is sphere-only) |
| `validate_pure_python.py` | From-scratch closest-point/separating-axis geometry, numpy only | All 12 shape pairs, 26 checks |

All three currently pass in full (57/57 checks).

## Why not PyGran / LIGGGHTS

PyGran has no physics of its own — its only working backend (`engine.liggghts`) is a ctypes wrapper
around a separately-built LIGGGHTS shared library; PyGran does not bundle or build it, and its other
backend (`engine.simple`) is an unimplemented stub. Building LIGGGHTS from source (available in the
AUR as `liggghts-git`) is a real, multi-minute, system-wide install with its own risk of failure —
disproportionate to what we're actually checking here. LAMMPS's own official `lammps` package (pip,
precompiled) already provides the same underlying granular contact physics LIGGGHTS forked from
(`pair_style granular`, `fix wall/gran`), so `validate_lammps_granular.py` calls it directly instead.

## Setup

**pybullet + LAMMPS** (one venv covers both):
```bash
python3 -m venv /path/to/venv
/path/to/venv/bin/pip install pybullet numpy lammps
```
The pip `lammps` wheel is built against the MPICH ABI (`libmpi.so.12`), which may not match the
system's OpenMPI (`libmpi.so.40`) — point `LD_LIBRARY_PATH` at any MPICH-ABI-compatible `libmpi.so.12`
already on the machine (e.g. an Intel oneAPI or bundled ParaView MPI install) when running the LAMMPS
script:
```bash
LD_LIBRARY_PATH="/opt/intel/oneapi/<version>/lib:$LD_LIBRARY_PATH" /path/to/venv/bin/python3 validate_lammps_granular.py
```
`validate_pybullet.py` and `validate_pure_python.py` need no such workaround.

**pure Python**: just needs `numpy` — no physics engine, no LD_LIBRARY_PATH workaround.

## Two real findings from writing these

1. **`Cylinder x HalfSpace`'s `contactNormal` sign was wrong before this validation pass.** It used
   to fall through to `GeneralContactDetection`'s generic HalfSpace-vs-anything template, which
   returns `-wallNormal` (into the wall) — the opposite sign from every other HalfSpace pair
   (`Sphere`/`Box` both return `+wallNormal`, confirmed correct by `pybullet` and, for `Sphere`, by
   `LAMMPS`'s own force direction). Fixed in `ParticleKernels.h` by giving `Cylinder x HalfSpace` its
   own overload using the same formula as `Box x HalfSpace`. All three scripts here confirm the fixed
   sign.
2. **Several `contactPoint` mismatches against `pybullet` are genuine mathematical degeneracies, not
   bugs**: a flat cylinder cap against a flat plane, or a flat box face against an infinite pillar's
   uniform-radius surface, has an entire *region* of equally-deep points, not a unique one. Our C++
   code picks one valid point from that region; `pybullet`'s independent GJK/EPA implementation picks
   a different, equally-valid one. `validate_pybullet.py` skips the point check for exactly these
   cases (documented inline), matching the precedent `ContactDetectionTest.cpp` itself already sets
   for `Box x Box` / `Cylinder x Box` (only the well-defined axis component is asserted there, for the
   same reason).
