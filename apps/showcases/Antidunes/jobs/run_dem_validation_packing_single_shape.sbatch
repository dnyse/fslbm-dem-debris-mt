#!/bin/bash -l
#SBATCH --job-name=antidunes-dem-validation-packing-single-shape
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --output=sbatch-output/output/antidunes_dem_validation_packing_single_shape_%j.out
#SBATCH --error=sbatch-output/error/antidunes_dem_validation_packing_single_shape_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# Fluid-free DEM validation dry-run, Kodam et al. (2010) Fig. 10 / Section 3.3 recreation at REAL,
# paper-matched scale (mm particles, real 9.81 gravity, 250 particles) -- see the prm file's header for
# the exact value-by-value correspondence. This scale needs a short timestep for stable stiff-particle
# contacts, so this is a genuinely long run (requested the full 24h singlenode partition limit; user
# explicitly chose full-fidelity-over-runtime for this comparison). ntasks=1: no MPI domain decomposition
# needed (single block, DEM-only, no LBM fields).
set -euo pipefail

if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    REPO="$(cd "$SLURM_SUBMIT_DIR/../../.." && pwd)"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    REPO="$(cd "$SCRIPT_DIR/../../../.." && pwd)"
fi
BIN_DIR="$REPO/build/apps/showcases/Antidunes"

module load gcc/11.5.0
module load openmpi/4.1.3-gcc11.5.0
module load cmake/3.31.8
module load python/3.13-miniforge
CONDA_BASE="$(conda info --base 2>/dev/null)"
# shellcheck source=/dev/null
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate walberla

cmake --build "$REPO/build" --target DEMValidation --parallel "$(nproc)"

cd "$BIN_DIR"
rm -rf vtk-out/packing_single_shape

# `|| true`: a scenario FAIL (nonzero exit) is a real, useful result, not a script error -- must not skip
# the archive step below. If SLURM kills this job on the 24h wall-time limit instead, the archive step
# never runs either, but the (still-being-written, vtkSpacing=2000) vtk-out/packing_single_shape/ tree is
# left in place under build/ and can be rescued manually afterward -- not lost, same recovery already done
# once for job 4048639.
./DEMValidation "$REPO/apps/showcases/Antidunes/prms/dem_validation_packing_single_shape.prm" || true

STAMP="$(date +%Y%m%d_%H%M%S)"
if [[ -n "${HPCVAULT:-}" ]]; then
    ARCHIVE_DIR="$HPCVAULT/antidunes-archive/${SLURM_JOB_NAME}_$STAMP"
else
    ARCHIVE_DIR="$BIN_DIR/archive/$STAMP"
fi
mkdir -p "$ARCHIVE_DIR"
[[ -d vtk-out/packing_single_shape ]] && mv vtk-out/packing_single_shape "$ARCHIVE_DIR/packing_single_shape"
