#!/bin/bash -l
#SBATCH --job-name=antidunes-dem-validation-packing-single-shape-sphere-smoke
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=00:15:00
#SBATCH --output=sbatch-output/output/antidunes_dem_validation_packing_single_shape_sphere_smoke_%j.out
#SBATCH --error=sbatch-output/error/antidunes_dem_validation_packing_single_shape_sphere_smoke_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# Short smoke test (10 particles, dropTimeoutSteps 5000) to check the sphere shape path works before
# committing to the full 24h, 250-particle Kodam-exact run -- see
# dem_validation_packing_single_shape_sphere_smoke.prm.
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
rm -rf vtk-out/packing_single_shape_sphere_smoke

./DEMValidation "$REPO/apps/showcases/Antidunes/prms/dem_validation_packing_single_shape_sphere_smoke.prm" || true

STAMP="$(date +%Y%m%d_%H%M%S)"
if [[ -n "${HPCVAULT:-}" ]]; then
    ARCHIVE_DIR="$HPCVAULT/antidunes-archive/${SLURM_JOB_NAME}_$STAMP"
else
    ARCHIVE_DIR="$BIN_DIR/archive/$STAMP"
fi
mkdir -p "$ARCHIVE_DIR"
[[ -d vtk-out/packing_single_shape_sphere_smoke ]] && mv vtk-out/packing_single_shape_sphere_smoke "$ARCHIVE_DIR/packing_single_shape_sphere_smoke"
