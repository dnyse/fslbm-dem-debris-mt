#!/bin/bash -l
#SBATCH --job-name=antidunes-dem-validation-rotating-drum-diag
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=00:15:00
#SBATCH --output=sbatch-output/output/antidunes_dem_validation_rotating_drum_diag_%j.out
#SBATCH --error=sbatch-output/error/antidunes_dem_validation_rotating_drum_diag_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# Short diagnostic run (settleSteps=20000, tumbleRevolutions=0 -- see dem_validation_rotating_drum_diag.prm)
# to localize the settle-phase energy-monotonicity FAIL job 4049208 reported for the full-length
# rotatingDrum run (worst checkpoint-to-checkpoint relative energy increase = 4.32). Uses the new
# per-checkpoint diagnostic logging added to runRotatingDrum's settle-phase loop. 15-minute budget is
# generous -- the full 400000-step settle phase itself only took ~3032s on this partition, so 20000 steps
# should finish in well under a minute of actual compute.
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
rm -rf vtk-out/rotating_drum_diag

./DEMValidation "$REPO/apps/showcases/Antidunes/prms/dem_validation_rotating_drum_diag.prm" || true

STAMP="$(date +%Y%m%d_%H%M%S)"
if [[ -n "${HPCVAULT:-}" ]]; then
    ARCHIVE_DIR="$HPCVAULT/antidunes-archive/${SLURM_JOB_NAME}_$STAMP"
else
    ARCHIVE_DIR="$BIN_DIR/archive/$STAMP"
fi
mkdir -p "$ARCHIVE_DIR"
[[ -d vtk-out/rotating_drum_diag ]] && mv vtk-out/rotating_drum_diag "$ARCHIVE_DIR/rotating_drum_diag"
