#!/bin/bash -l
#SBATCH --job-name=antidunes-dem-validation-rotating-drum
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --output=sbatch-output/output/antidunes_dem_validation_rotating_drum_%j.out
#SBATCH --error=sbatch-output/error/antidunes_dem_validation_rotating_drum_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# Fluid-free DEM validation dry-run, Kodam et al. (2010) Fig. 11 / Section 3.4 recreation at REAL,
# paper-matched scale (mm particles, real 9.81 gravity, 250 particles, 25rpm) -- see the prm file's header
# for the exact value-by-value correspondence. Genuinely long run (24h singlenode partition limit
# requested; user explicitly chose full-fidelity-over-runtime). ntasks=1: no MPI domain decomposition
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
rm -rf vtk-out/rotating_drum

# `|| true`: a scenario FAIL (nonzero exit) is a real, useful result, not a script error -- must not skip
# the archive step below. If SLURM kills this job on the 24h wall-time limit instead, the archive step
# never runs either, but the (still-being-written, vtkSpacing=2000) vtk-out/rotating_drum/ tree is left in
# place under build/ and can be rescued manually afterward -- not lost.
./DEMValidation "$REPO/apps/showcases/Antidunes/prms/dem_validation_rotating_drum.prm" || true

STAMP="$(date +%Y%m%d_%H%M%S)"
if [[ -n "${HPCVAULT:-}" ]]; then
    ARCHIVE_DIR="$HPCVAULT/antidunes-archive/${SLURM_JOB_NAME}_$STAMP"
else
    ARCHIVE_DIR="$BIN_DIR/archive/$STAMP"
fi
mkdir -p "$ARCHIVE_DIR"
[[ -d vtk-out/rotating_drum ]] && mv vtk-out/rotating_drum "$ARCHIVE_DIR/rotating_drum"
