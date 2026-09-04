#!/bin/bash -l
#SBATCH --job-name=antidunes-dem-validation-packing-single-shape-sphere
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --output=sbatch-output/output/antidunes_dem_validation_packing_single_shape_sphere_%j.out
#SBATCH --error=sbatch-output/error/antidunes_dem_validation_packing_single_shape_sphere_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# Sphere-shape variant of run_dem_validation_packing_single_shape.sbatch -- see
# dem_validation_packing_single_shape_sphere.prm for why (between-shapes comparison + a simpler-contact
# control against the cylinder run's create=100/250, bulkSolidFraction=1.34 anomaly).
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
rm -rf vtk-out/packing_single_shape_sphere

./DEMValidation "$REPO/apps/showcases/Antidunes/prms/dem_validation_packing_single_shape_sphere.prm" || true

STAMP="$(date +%Y%m%d_%H%M%S)"
if [[ -n "${HPCVAULT:-}" ]]; then
    ARCHIVE_DIR="$HPCVAULT/antidunes-archive/${SLURM_JOB_NAME}_$STAMP"
else
    ARCHIVE_DIR="$BIN_DIR/archive/$STAMP"
fi
mkdir -p "$ARCHIVE_DIR"
[[ -d vtk-out/packing_single_shape_sphere ]] && mv vtk-out/packing_single_shape_sphere "$ARCHIVE_DIR/packing_single_shape_sphere"
