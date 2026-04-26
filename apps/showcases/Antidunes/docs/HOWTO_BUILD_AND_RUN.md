# Antidunes Showcase — Build & Run Guide

## Prerequisites

Load the required modules on the supercomputer before doing anything:

```bash
module load gcc/14.2.0
module load openmpi/4.1.3-gcc14.2.0
module load cmake
module load python
```

---

## Python / lbmpy Setup (required for code generation)

The `Antidunes` target uses pystencils/lbmpy to generate a lattice model at build time.
The system Python 3.9 (anaconda) is pinned to lbmpy 1.3.3, which is too old.
Create a conda environment with a newer version:

```bash
conda create -n walberla-codegen python=3.12 -y
conda activate walberla-codegen
pip install "lbmpy>=1.4"
pip install jinja2
```

---

## Configure (from the repo root)

Run all CMake commands from `/path/to/walberla/` (the repository root).

```bash
cmake -B build/ \
  -DCMAKE_C_COMPILER=$(which gcc) \
  -DCMAKE_CXX_COMPILER=$(which g++) \
  -DPython_ROOT_DIR=$CONDA_PREFIX \
  -DWALBERLA_BUILD_SHOWCASES=ON \
  -DWALBERLA_BUILD_WITH_CODEGEN=ON \
  -DCMAKE_BUILD_TYPE=Release \
  .
```

**Why explicit compilers?** Without them CMake picks up `/usr/bin/cc` (GCC 8.5.0),
which does not support C++20. The loaded `gcc/14.2.0` must be used.

To force a clean reconfiguration without deleting compiled objects:

```bash
cmake --fresh -B build/ \
  -DCMAKE_C_COMPILER=$(which gcc) \
  -DCMAKE_CXX_COMPILER=$(which g++) \
  -DPython_ROOT_DIR=$CONDA_PREFIX \
  -DWALBERLA_BUILD_SHOWCASES=ON \
  -DWALBERLA_BUILD_WITH_CODEGEN=ON \
  -DCMAKE_BUILD_TYPE=Release \
  .
```

### Key CMake flags

| Flag | Value | Reason |
|------|-------|--------|
| `WALBERLA_BUILD_SHOWCASES` | `ON` | Enables showcase apps |
| `WALBERLA_BUILD_WITH_CODEGEN` | `ON` | Required for `Antidunes` target |
| `Python_ROOT_DIR` | `$CONDA_PREFIX` | Points to the conda env with lbmpy>=1.4 |
| `CMAKE_C_COMPILER` / `CMAKE_CXX_COMPILER` | `$(which gcc/g++)` | Forces GCC 14, not system GCC 8 |

---

## Build

```bash
cmake --build build/ --target Antidunes    -- -j$(nproc)
cmake --build build/ --target BedGeneration -- -j$(nproc)
```

`BedGeneration` does **not** require `WALBERLA_BUILD_WITH_CODEGEN` and can be built independently.

---

## Running

All commands below must be run from the build output directory:

```bash
cd build/apps/showcases/Antidunes
```

### Step 1 — Generate the particle bed

`Antidunes` requires a `spheres_out.dat` file. Generate it first:

```bash
mpirun -n 9 ./BedGeneration BedGeneration.prm
```

`BedGeneration.prm` hardcodes `blocks <3, 3, 1>` → exactly 9 MPI processes required.
The output file `spheres_out.dat` is written to the current directory.

### Step 2 — Run the simulation

```bash
mpirun -n 72 ./Antidunes Antidunes.prm
```

---

## Domain / Process Count Constraints

The number of MPI processes must equal the number of blocks exactly:

```
numProcesses == (domainSize[0]/cellsPerBlock[0])
              * (domainSize[1]/cellsPerBlock[1])
              * (domainSize[2]/cellsPerBlock[2])
```

With `cellsPerBlock <50, 20, 40>` the original domain `<3200, 60, 160>` gives
`64 × 3 × 4 = 768` blocks (the production setup).

### Scaling down to 72 processes

The particle bed also imposes a compatibility constraint driven by a scaling factor:

```
scalingFactor  = domainSize[0] / (bedDomainSize_SI_X * bedCopiesInX)
               = domainSize[0] / (0.8 * bedCopiesInX)

# Must hold exactly (integer check):
scalingFactor * bedDomainSize_SI_Y * bedCopiesInY == domainSize[1]
              = scalingFactor * 0.015 * bedCopiesInY
```

For **72 processes** the valid configuration is:

| Parameter | Value |
|-----------|-------|
| `domainSize` | `<400, 60, 120>` |
| `cellsPerBlock` | `<50, 20, 40>` (unchanged) |
| `bedCopiesInX` | `1` |
| `bedCopiesInY` | `8` |
| blocks | `8 × 3 × 3 = 72` ✓ |
| scalingFactor | `400 / 0.8 = 500` |
| Y check | `500 × 0.015 × 8 = 60` ✓ |

The bed was generated with `domainSize_SI_Y = 0.015 m`. At scalingFactor=500,
one bed copy spans only 7.5 lattice cells — not an integer — so 8 copies are needed
to exactly fill the 60-cell Y dimension.

### General formula for other process counts

Given a target process count `N` and `cellsPerBlock <50, 20, 40>`:

1. Find integer factorization `bx * by * bz = N`
2. Set `domainSize = <50*bx, 20*by, 40*bz>`
3. Find `bedCopiesInY` such that `(50*bx / 0.8) * 0.015 * bedCopiesInY` is an integer equal to `20*by`
4. Adjust `bedCopiesInX` analogously if needed

---

## Output

| Directory | Contents |
|-----------|----------|
| `vtk-out/` | VTK files (fluid field: fill level, velocity, density, mapped flag) every 1000 steps |
| `eval/` | `bedload.txt`, `fluidInfo.txt` evaluation data every 1000 steps |
| `snapshots/` | Checkpoint snapshots every 10000 steps |

---

## SLURM Job Template (production: 768 processes)

```bash
#!/bin/bash
#SBATCH --job-name=antidunes
#SBATCH --nodes=11
#SBATCH --ntasks=768
#SBATCH --time=48:00:00
#SBATCH --output=antidunes_%j.out

module load gcc/14.2.0
module load openmpi/4.1.3-gcc14.2.0
conda activate walberla-codegen

cd build/apps/showcases/Antidunes
mpirun ./Antidunes Antidunes.prm
```