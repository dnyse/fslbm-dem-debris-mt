#!/bin/bash -l
#SBATCH --job-name=antidunes-dem-validation-packing-angleofrepose
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=00:30:00
#SBATCH --output=sbatch-output/output/antidunes_dem_validation_%j.out
#SBATCH --error=sbatch-output/error/antidunes_dem_validation_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# Fluid-free DEM validation dry-run (DEMValidation.cpp / DEMValidationScenarios.h) -- no MPI domain
# decomposition needed (single block, DEM-only, no LBM fields), so ntasks=1. Runs the existing
# "packing" (Kodam et al. 2010 Fig. 10 analogue: mixed sphere/box/cylinder bed settling into a bounded
# container) and "angleOfRepose" (Fig. 11 analogue: static pour-pile repose angle) scenarios back to
# back, each with VTK output enabled (vtkSpacing>0 in the prm), then archives vtk-out/{packing,
# angle_of_repose} out of the shared build/ tree so a later run on this or another mirror can't clobber
# them -- same tmp-then-archive pattern the other jobs use, just DEM-scoped subfolders instead of the
# whole vtk-out/eval/snapshots set (this app writes no eval/ or snapshots/ at all).
set -euo pipefail

if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    REPO="$(cd "$SLURM_SUBMIT_DIR/../../.." && pwd)"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    REPO="$(cd "$SCRIPT_DIR/../../../.." && pwd)"
fi
BIN_DIR="$REPO/build/apps/showcases/Antidunes"

# antidunes.sh's own module-loading block (it doesn't build the DEMValidation target itself, so this
# job calls cmake directly and needs the same environment set up first).
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
rm -rf vtk-out/packing vtk-out/angle_of_repose

# `|| true` on each run: DEMValidation exits nonzero on a scenario FAIL (e.g. angleOfRepose's known,
# already-documented no-rolling-friction limitation -- see DEMValidationScenarios.h), which is a real
# and useful result, not a script error. Under `set -e` an unguarded nonzero exit here would abort the
# script BEFORE the archive step below ever runs, silently losing BOTH scenarios' VTK output even when
# the other one passed cleanly (hit exactly this the first time this job ran, job 4048639).
./DEMValidation "$REPO/apps/showcases/Antidunes/prms/dem_validation_packing.prm" || true
./DEMValidation "$REPO/apps/showcases/Antidunes/prms/dem_validation_angle_of_repose.prm" || true

STAMP="$(date +%Y%m%d_%H%M%S)"
if [[ -n "${HPCVAULT:-}" ]]; then
    ARCHIVE_DIR="$HPCVAULT/antidunes-archive/${SLURM_JOB_NAME}_$STAMP"
else
    ARCHIVE_DIR="$BIN_DIR/archive/$STAMP"
fi
mkdir -p "$ARCHIVE_DIR"
for dir in vtk-out/packing vtk-out/angle_of_repose; do
    [[ -d "$dir" ]] && mv "$dir" "$ARCHIVE_DIR/$(basename "$dir")"
done
