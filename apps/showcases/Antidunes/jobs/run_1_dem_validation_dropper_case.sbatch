#!/bin/bash -l
#SBATCH --job-name=antidunes-dem-validation-dropper-case
#SBATCH --partition=singlenode
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=00:15:00
#SBATCH --output=sbatch-output/output/dem_validation_dropper_case_%j.out
#SBATCH --error=sbatch-output/error/dem_validation_dropper_case_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# Final thesis DEM validation "dropper case" showcase (see FINAL_THESIS_DEM_VALIDATION.md item 4) --
# one clean cylinder-wall and one box-wall oblique impact, VTK-enabled, for a thesis figure. Cluster-
# only per explicit instruction, even though this is a near-instantaneous run (both cases together ran
# in well under a minute locally without VTK) -- 1 rank is genuinely enough: DEMValidation's
# wallImpactCylinder/wallImpactBox test cases are plain mesa_pd single/two-body impacts (DEMEngine,
# DEMValidationScenarios.h), with no block-forest domain decomposition at all -- unlike DropTest.cpp
# (jobs/run_72_droptest.sbatch), which specifically needs 72 ranks to replicate production's own
# block-forest decomposition. Time limit 15min is a large safety buffer for what's actually a
# sub-minute job.
#
# Built directly via cmake/mpirun rather than antidunes.sh, since that script is hardcoded to the
# Antidunes/BedGeneration/tsunami target names (same reasoning as run_72_droptest.sbatch's own header).
set -euo pipefail

module load gcc/11.5.0
module load openmpi/4.1.3-gcc11.5.0
module load cmake/3.31.8
module load python/3.13-miniforge

if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    REPO="$(cd "$SLURM_SUBMIT_DIR/../../.." && pwd)"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    REPO="$(cd "$SCRIPT_DIR/../../../.." && pwd)"
fi
BUILD_DIR="$REPO/build"
BIN_DIR="$BUILD_DIR/apps/showcases/Antidunes"
PRM_DIR="$REPO/apps/showcases/Antidunes/prms"

cmake --build "$BUILD_DIR" --target DEMValidation -- -j"$(nproc)"

cd "$BIN_DIR"

STAMP="$(date +%Y%m%d_%H%M%S)"
if [[ -n "${HPCVAULT:-}" ]] && mkdir -p "$HPCVAULT/antidunes-archive/${SLURM_JOB_NAME}_${SLURM_JOB_ID}_$STAMP" 2>/dev/null; then
    ARCHIVE_ROOT="$HPCVAULT/antidunes-archive/${SLURM_JOB_NAME}_${SLURM_JOB_ID}_$STAMP"
else
    ARCHIVE_ROOT="$BIN_DIR/archive/$STAMP"
    mkdir -p "$ARCHIVE_ROOT"
fi

run_case() {  # $1 = banner, $2 = prm file, $3 = vtk output dir name (outputBaseName + the code's own
              # "_tilt<value>" suffix, e.g. "dropper_case_cylinder_tilt45.000000" -- confirmed by an
              # actual local run, not assumed)
    echo "############################################################"
    echo "## $1"
    echo "############################################################"
    mpirun -n 1 ./DEMValidation "$PRM_DIR/$2"
    if [[ -d "vtk-out/$3" ]]; then
        mkdir -p "$ARCHIVE_ROOT"
        mv "vtk-out/$3" "$ARCHIVE_ROOT/$3"
    fi
}

run_case "Dropper case: cylinder" dem_validation_dropper_case_cylinder.prm "dropper_case_cylinder_tilt45.000000"
run_case "Dropper case: box"      dem_validation_dropper_case_box.prm      "dropper_case_box_tilt25.000000"

echo "## DONE. Archived under $ARCHIVE_ROOT/{dropper_case_cylinder,dropper_case_box}"
