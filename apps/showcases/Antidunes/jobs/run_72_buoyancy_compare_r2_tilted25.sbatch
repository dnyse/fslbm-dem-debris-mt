#!/bin/bash -l
#SBATCH --job-name=antidunes-buoyancy-r2-compare-tilted25
#SBATCH --partition=singlenode
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=01:30:00
#SBATCH --output=sbatch-output/output/buoyancy_compare_r2_tilted25_%j.out
#SBATCH --error=sbatch-output/error/buoyancy_compare_r2_tilted25_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# [Claude, 2026-09-22] TILTED25 VARIANT (part of a chained R1->R2->R4 rerun, THESIS_PLAN.md): the box runs (both mem and analytic
# arms) use the tilted25 prms -- initialTiltAngleDeg 25 on the FloatingBox -- instead of upright. Root cause: the untilted box sits
# on a mapping tie (axis-aligned face planes through cell centres), which for R1/MEM left it flush on the floor with no fluid cell
# underneath, starving MEM-only buoyancy (traced and fixed 2026-09-22, 'CORRECTION of the tie explanation'; proven once for R1/box/
# MEM only, job 4232331, never applied to R2/R4 or the analytic arm before now). Sphere/cylinder prms are UNCHANGED (no tie problem
# was ever observed for them). See gen_buoyancy_tilted25_box_prm.py.
# Chained: R1 (no dependency) -> R2 (afterany R1) -> R4 (afterany R2), same checkout (mt-dennys-mirror, which already has all
# three bed files spheres_out_buoyancy_r{1,2,4}.dat -- avoids a repeat of the earlier cross-mirror bed-file gap).
# Time limit ADJUSTED from the last real run: job 4223072 (the R2 compare-all that actually COMPLETED) took 1:08:31 total for the 6 untilted runs; the
# previous script's own 02:30:00 budget was padded well beyond that. Tightened to 01:30:00 (~22 min margin over the last measured total) so this gets scheduled sooner.

# Thesis item 1b (THESIS_PLAN.md) -- MEM vs. analytical buoyancy comparison, R=2, all 3 shapes x 2
# methods run back-to-back in ONE job (matches item 1a's own run_case() sbatch precedent) so the
# node allocation only needs to be requested once, not 6 separate times. R=2: 72 ranks / 1 node,
# cellsPerBlock <32,32,24> (8x baseline per-rank cells, single-node compromise). Time limit 02:30:00 = ~2.3x the extrapolated ~66min total
# (6 x ~11min each, per THESIS_PLAN.md's own timing data) as a safety buffer.
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

run_case "R=2 sphere MEM"        buoyancy_compare_sphere_r2_mem_gas.prm        sphere_mem
run_case "R=2 sphere analytic"   buoyancy_compare_sphere_r2_analytic_gas.prm   sphere_analytic
run_case "R=2 cylinder MEM"      buoyancy_compare_cylinder_r2_mem_gas.prm      cylinder_mem
run_case "R=2 cylinder analytic" buoyancy_compare_cylinder_r2_analytic_gas.prm cylinder_analytic
run_case "R=2 box MEM tilted25"           buoyancy_compare_box_r2_mem_tilted25_gas.prm      box_mem
run_case "R=2 box analytic tilted25"      buoyancy_compare_box_r2_analytic_tilted25_gas.prm box_analytic

echo "## DONE. Archived under $ARCHIVE_ROOT/{sphere,cylinder,box}_{mem,analytic}_{vtk-out,eval,snapshots}"
