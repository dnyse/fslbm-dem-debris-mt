#!/bin/bash -l
#SBATCH --job-name=antidunes-buoyancy-r4-compare-tilted25
#SBATCH --partition=multinode
#SBATCH --nodes=4
#SBATCH --ntasks=288
#SBATCH --time=14:00:00
#SBATCH --output=sbatch-output/output/buoyancy_compare_r4_tilted25_%j.out
#SBATCH --error=sbatch-output/error/buoyancy_compare_r4_tilted25_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# [Claude, 2026-09-22] TILTED25 VARIANT (part of a chained R1->R2->R4 rerun, THESIS_PLAN.md): the box runs (both mem and analytic
# arms) use the tilted25 prms -- initialTiltAngleDeg 25 on the FloatingBox -- instead of upright. Root cause: the untilted box sits
# on a mapping tie (axis-aligned face planes through cell centres), which for R1/MEM left it flush on the floor with no fluid cell
# underneath, starving MEM-only buoyancy (traced and fixed 2026-09-22, 'CORRECTION of the tie explanation'; proven once for R1/box/
# MEM only, job 4232331, never applied to R2/R4 or the analytic arm before now). Sphere/cylinder prms are UNCHANGED (no tie problem
# was ever observed for them). See gen_buoyancy_tilted25_box_prm.py.
# [Claude, 2026-09-29] Switched from 576 ranks/8 nodes to 288 ranks/4 nodes (cellsPerBlock <32,32,24>
# -> <32,32,48> in all R4 prms, 6x6x16 -> 6x6x8) per explicit request ("run r4 on four nodes, the
# other on a single node"); time limit doubled to 14:00:00 accordingly (half the ranks, roughly
# double the per-case wall time, plus the same margin as the original 576-rank estimate).
# Chained: R1 (no dependency) -> R2 (afterany R1) -> R4 (afterany R2), same checkout (mt-dennys-mirror, which already has all
# three bed files spheres_out_buoyancy_r{1,2,4}.dat -- avoids a repeat of the earlier cross-mirror bed-file gap).
# Time limit ADJUSTED from the last real run: job 4235003 (the R4 compare-all that actually COMPLETED, after two earlier failed/cancelled attempts --
# 4234535 4:30:33 CANCELLED, 4227812 0:52:58 FAILED, 4223073 0:03:32 FAILED) took 6:18:24 total for the 6 untilted runs. Tightened from 09:00:00 to 07:00:00
# (~42 min margin over the last measured total) so this gets scheduled sooner.

# Thesis item 1b (THESIS_PLAN.md) -- MEM vs. analytical buoyancy comparison, R=4, all 3 shapes x 2
# methods run back-to-back in ONE job (matches item 1a's own run_case() sbatch precedent) so the
# node allocation only needs to be requested once, not 6 separate times. R=4: 576 ranks / 8 nodes,
# cellsPerBlock <32,32,24> (8x baseline per-rank cells, partial multi-node scale-up, not the
# full R^3/64-node scaling). Time limit 09:00:00 = ~2.2x the extrapolated ~4.1h total
# (6 x ~41min each, per THESIS_PLAN.md's own timing data) as a safety buffer.
# --constraint=el9 deliberately omitted (see run_72_bedgen_buoyancy_r4.sbatch's own header comment).
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
    "$SH" --run --np 288 --prm "$PRM_DIR/$2"
    archive_outputs "$3"
}

run_case "R=4 sphere MEM"        buoyancy_compare_sphere_r4_mem_gas.prm        sphere_mem
run_case "R=4 sphere analytic"   buoyancy_compare_sphere_r4_analytic_gas.prm   sphere_analytic
run_case "R=4 cylinder MEM"      buoyancy_compare_cylinder_r4_mem_gas.prm      cylinder_mem
run_case "R=4 cylinder analytic" buoyancy_compare_cylinder_r4_analytic_gas.prm cylinder_analytic
run_case "R=4 box MEM tilted25"           buoyancy_compare_box_r4_mem_tilted25_gas.prm      box_mem
run_case "R=4 box analytic tilted25"      buoyancy_compare_box_r4_analytic_tilted25_gas.prm box_analytic

echo "## DONE. Archived under $ARCHIVE_ROOT/{sphere,cylinder,box}_{mem,analytic}_{vtk-out,eval,snapshots}"
