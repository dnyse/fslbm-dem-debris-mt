# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is the **Antidunes** showcase — a parallel 3D numerical simulation of antidune bedforms that couples three physics solvers:
- **Free-Surface Lattice Boltzmann Method (FSLBM)** for two-phase fluid flow
- **Discrete Element Method (DEM)** for granular particle dynamics
- **Momentum Exchange Method (MEM)** for fluid-particle coupling

It lives inside the **waLBerla** monorepo at `apps/showcases/Antidunes/`. There is no standalone build; it compiles as part of the waLBerla CMake project (root: `/home/dhuber/University/Thesis/mt-dennys`).

Published results: https://doi.org/10.1017/jfm.2023.262

## Build

**Prerequisites:** C++20 compiler, CMake >= 3.26, MPI library, Python >= 3.10 with `pystencils~=1.3`, `lbmpy~=1.3`, `jinja2` (needed for LBM kernel code generation at build time).

```bash
# From the waLBerla repo root
python -m venv .venv && source .venv/bin/activate
pip install "pystencils~=1.3" "lbmpy~=1.3" jinja2

cmake -G Ninja \
  -DWALBERLA_BUILD_WITH_CODEGEN=ON \
  -DWALBERLA_BUILD_SHOWCASES=ON \
  -DWALBERLA_BUILD_WITH_MPI=ON \
  -DPython_ROOT_DIR="$(pwd)/.venv/bin" \
  -B build/ .

cmake --build build/ --target Antidunes
cmake --build build/ --target BedGeneration
```

`WALBERLA_BUILD_WITH_CODEGEN=ON` is required — the D3Q27 Cumulant LBM collision kernel is generated at build time by `AntidunesLatticeModelGeneration.py` via lbmpy.

## Running

The number of MPI processes must exactly equal `(domainX/cellsPerBlockX) * (domainY/cellsPerBlockY) * (domainZ/cellsPerBlockZ)`. The default config requires 768 processes (64 × 3 × 4). `setupBlockForest` enforces this with an abort.

```bash
# Step 1: generate settled particle bed (prerequisite for main sim)
mpirun -np <N> ./BedGeneration BedGeneration.prm

# Step 2: run the main simulation
mpirun -np <N> ./Antidunes Antidunes.prm
```

Four experimental configurations (E1–E4) are documented as comments in `Antidunes.prm` with their Re, Fr, We values.

## Architecture

### Key Files

| File | Role |
|---|---|
| `Antidunes.cpp` | Entry point; orchestrates setup stages and the main time loop |
| `BedGeneration.cpp` | Standalone preprocessing tool; settles particles via DEM, writes `spheres_out.dat` |
| `SimulationParams.h` | All configuration parameters parsed from `.prm` file; `ComputedParams` derives all physical LBM units from dimensionless numbers |
| `SimulationObjects.h` | Holds all mutable runtime objects (block forest, field IDs, particle storage, PID controller, etc.) |
| `SetupStages.h` | Initialization functions (`setupBlockForest`, `setupLBMFields`, `setupParticles`, `setupBoundaryHandling`, etc.) |
| `BuildTimeloop.h` | Adds sweeps and evaluators to the three separate time loops |
| `ParticleKernels.h` | DEM sub-cycling: velocity Verlet, lubrication correction, contact forces, MPI sync |
| `AntidunesBoundaryHandling.h` | Extends waLBerla's free-surface boundary handling with NoSlip, FreeSlip, UBB, Pressure, MovingObstacle |
| `PIDController.h/.cpp` | Discrete PID controller maintaining target mean flow velocity via body force adjustment |
| `Forcing.h` | Body force density sweep, mean velocity computation, PID→acceleration field connection |
| `Evaluators.h` | XZ-slice averages, bedload transport, total hydrodynamic force evaluators |
| `WriteData.h` | Checkpoint writing (tmp-then-rename for atomicity) and fluid info to SQLite |
| `DataInitialization.h` | Field initialization: sinusoidal fill level, Poiseuille velocity profile, hydrostatic pressure |
| `AntidunesLatticeModelGeneration.py` | lbmpy script generating the D3Q27 Cumulant LBM kernel; run automatically by CMake |

### Per-Timestep Execution Order

```
timeloopFSLBM->singleStep()       # Free-surface LBM sweeps
particleKernels.reduceAndAvgHydroForces()
particleKernels.doDEMSubCycles(N=10)
particleKernels.resetHydroForces()
timeloopAfterDEM->singleStep()    # LBM-DEM coupling (particle mapping + PDF reconstruction)
timeloopEval->singleStep()        # Diagnostics, PID update
```

### Key Conventions

- **`SimulationParams` is immutable after construction**; all mutable state lives in `SimulationObjects`. Both are passed as `shared_ptr` everywhere.
- **Header-only implementation**: almost all logic is in `.h` files. Only `PIDController.cpp` and the two entry-point `.cpp` files are non-header.
- **All unit conversions are centralized in `ComputedParams`** — viscosity, relaxation rate, gravity, surface tension, SI scaling are all derived once from the dimensionless numbers Re, Fr, We.
- **Three separate `SweepTimeloop`s** (`timeloopFSLBM`, `timeloopAfterDEM`, `timeloopEval`) run sequentially via `.singleStep()` each timestep, making the coupling order explicit.
- **MPI-only parallelism** — no OpenMP by default. All file I/O uses `WALBERLA_ROOT_SECTION()` (rank 0 only).
- All code is under namespace `walberla::antidunes`.

### Outputs

| Output | Path | Contents |
|---|---|---|
| Checkpoints | `snapshots/` | Block forest, PDF/fill/excess-mass fields, particle storage, PID state |
| Evaluation | `eval/` | `antidunes.sqlite`, computed params, XZ-slice averages per timestep |
| Visualization | `vtk-out/` | Fluid fields and particle point data for ParaView |
| Bed file | `.` | `spheres_out.dat` (particle positions + radii in SI units) |

Post-processing notebooks: `slice_evaluation.ipynb` and `pyvista.ipynb`.
