#!/bin/bash -l
#SBATCH --job-name=antidunes-buoyancy-box-r2-mem
#SBATCH --partition=singlenode
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=01:30:00
#SBATCH --output=sbatch-output/output/buoyancy_compare_box_r2_mem_%j.out
#SBATCH --error=sbatch-output/error/buoyancy_compare_box_r2_mem_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# Thesis item 1b (THESIS_PLAN.md) -- MEM vs. analytical buoyancy comparison, box, R=2, mem
# method only. The original R=2 combined-job log (job 4223072) that this data came from was
# rotated off the cluster before its own F_hydro_z noise statistics (mean/stdev/range) were ever
# extracted -- only rise/v_z had been pulled at the time. This standalone rerun exists purely to
# regenerate that one log so the box-MEM noise trend (R=1 -> R=2 -> R=4) can be completed for the
# thesis results section (chapters/04_results/equilibrium/mem-vs-analytical-buoyancy.tex).
# Same prm, same np=72/single-node as the original R=2 combined job (run_72_buoyancy_compare_r2_all.sbatch).
# --constraint=el9 deliberately omitted (see run_72_bedgen_buoyancy_r2.sbatch's own header comment).
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
        echo "WARNING: removing leftover $BIN_DIR/$dir from an unarchived previous run" >&2
        rm -rf "${BIN_DIR:?}/$dir"
    fi
done

"$SH" --build --jobs "$(nproc)" --run --np 72 --prm "$PRM_DIR/buoyancy_compare_box_r2_mem_gas.prm"

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

echo "## DONE. Archived under $ARCHIVE_DIR"
