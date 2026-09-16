#!/bin/bash -l
#SBATCH --job-name=antidunes-buoyancy-box-r1-mem-tilted25
#SBATCH --partition=singlenode
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=00:45:00
#SBATCH --output=sbatch-output/output/buoyancy_compare_box_r1_mem_tilted25_%j.out
#SBATCH --error=sbatch-output/error/buoyancy_compare_box_r1_mem_tilted25_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# Diagnostic: same as run_72_buoyancy_compare_box_r1_mem.sbatch but initialTiltAngleDeg=25 (about
# world Y) on the box -- see prms/buoyancy_compare_box_r1_mem_tilted25_gas.prm's own header for the
# full hypothesis being tested (does tilting the box off-axis reduce the F_hydro_z noise the
# axis-aligned box shows at R=1, THESIS_PLAN.md item 1b).
# --constraint=el9 deliberately omitted (see run_72_bedgen_buoyancy_r1.sbatch's own header comment).
set -euo pipefail

if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    REPO="$(cd "$SLURM_SUBMIT_DIR/../../.." && pwd)"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    REPO="$(cd "$SCRIPT_DIR/../../../.." && pwd)"
fi
SH="$REPO/apps/showcases/Antidunes/antidunes.sh"
PRM_DIR="$REPO/apps/showcases/Antidunes/prms"
BIN_DIR="$REPO/build/apps/showcases/Antidunes"

for dir in vtk-out eval snapshots; do
    if [[ -d "$BIN_DIR/$dir" ]]; then
        rm -rf "${BIN_DIR:?}/$dir"
    fi
done

"$SH" --build --jobs "$(nproc)" --run --np 72 --prm "$PRM_DIR/buoyancy_compare_box_r1_mem_tilted25_gas.prm"

STAMP="$(date +%Y%m%d_%H%M%S)"
if [[ -n "${HPCVAULT:-}" ]] && mkdir -p "$HPCVAULT/antidunes-archive/${SLURM_JOB_NAME}_${SLURM_JOB_ID}_$STAMP" 2>/dev/null; then
    ARCHIVE_DIR="$HPCVAULT/antidunes-archive/${SLURM_JOB_NAME}_$STAMP"
else
    ARCHIVE_DIR="$BIN_DIR/archive/$STAMP"
    mkdir -p "$ARCHIVE_DIR"
fi
for dir in vtk-out eval snapshots; do
    if [[ -d "$BIN_DIR/$dir" ]]; then
        mv "$BIN_DIR/$dir" "$ARCHIVE_DIR/$dir"
    fi
done
