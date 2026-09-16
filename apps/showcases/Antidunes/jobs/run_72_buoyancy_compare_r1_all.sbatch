#!/bin/bash -l
#SBATCH --job-name=antidunes-buoyancy-r1-compare-all
#SBATCH --partition=singlenode
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=00:30:00
#SBATCH --output=sbatch-output/output/buoyancy_compare_r1_all_%j.out
#SBATCH --error=sbatch-output/error/buoyancy_compare_r1_all_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# Thesis item 1b (THESIS_PLAN.md) -- MEM vs. analytical buoyancy comparison, R=1, all 3 shapes x 2
# methods run back-to-back in ONE job (matches item 1a's own run_case() sbatch precedent) so the
# node allocation only needs to be requested once, not 6 separate times. R=1: 72 ranks / 1 node,
# cellsPerBlock <16,16,12> (true baseline). Time limit 00:30:00 = ~4x the measured ~7.6min total
# (6 x ~76s each, per THESIS_PLAN.md's own timing data) as a safety buffer.
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
        echo "WARNING: removing leftover $BIN_DIR/$dir from an unarchived previous run" >&2
        rm -rf "${BIN_DIR:?}/$dir"
    fi
done

STAMP="$(date +%Y%m%d_%H%M%S)"
if [[ -n "${HPCVAULT:-}" ]] && mkdir -p "$HPCVAULT/antidunes-archive/${SLURM_JOB_NAME}_${SLURM_JOB_ID}_$STAMP" 2>/dev/null; then
    ARCHIVE_ROOT="$HPCVAULT/antidunes-archive/${SLURM_JOB_NAME}_${SLURM_JOB_ID}_$STAMP"
else
    ARCHIVE_ROOT="$BIN_DIR/archive/$STAMP"
    mkdir -p "$ARCHIVE_ROOT"
fi

"$SH" --build --jobs "$(nproc)"

archive_outputs() {  # $1 = label
    for d in vtk-out eval snapshots; do
        [[ -d "$BIN_DIR/$d" ]] && mkdir -p "$ARCHIVE_ROOT/$1" && mv "$BIN_DIR/$d" "$ARCHIVE_ROOT/$1/$d"
    done
}

run_case() {  # $1 = banner, $2 = prm file, $3 = label
    echo "############################################################"
    echo "## $1"
    echo "############################################################"
    "$SH" --run --np 72 --prm "$PRM_DIR/$2"
    archive_outputs "$3"
}

run_case "R=1 sphere MEM"        buoyancy_compare_sphere_r1_mem_gas.prm        sphere_mem
run_case "R=1 sphere analytic"   buoyancy_compare_sphere_r1_analytic_gas.prm   sphere_analytic
run_case "R=1 cylinder MEM"      buoyancy_compare_cylinder_r1_mem_gas.prm      cylinder_mem
run_case "R=1 cylinder analytic" buoyancy_compare_cylinder_r1_analytic_gas.prm cylinder_analytic
run_case "R=1 box MEM"           buoyancy_compare_box_r1_mem_gas.prm           box_mem
run_case "R=1 box analytic"      buoyancy_compare_box_r1_analytic_gas.prm      box_analytic

echo "## DONE. Archived under $ARCHIVE_ROOT/{sphere,cylinder,box}_{mem,analytic}_{vtk-out,eval,snapshots}"
