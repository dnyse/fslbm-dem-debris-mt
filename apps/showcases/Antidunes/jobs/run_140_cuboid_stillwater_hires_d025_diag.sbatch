#!/bin/bash -l
#SBATCH --job-name=antidunes-cuboid-stillwater-hires-d025-diag
#SBATCH --partition=multinode
#SBATCH --nodes=2
#SBATCH --ntasks=140
#SBATCH --time=01:00:00
#SBATCH --output=sbatch-output/output/cuboid_stillwater_hires_d025_diag_%j.out
#SBATCH --error=sbatch-output/error/cuboid_stillwater_hires_d025_diag_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# Thesis item 1a (THESIS_PLAN.md) -- diagnostic resume of job 4236462's own checkpoint (t=70000,
# still on disk in build/apps/showcases/Antidunes/snapshots/ since the job crashed before its own
# archive-move step ran) to observe the onset of the `[NaN-SCAN] aborting after timeloopFSLBM at
# t=72000` crash with denseNanScanEnabled bracketing 70000-72500. See
# cuboid_stillwater_realgravity981_3xsphere_hires_d025_diag_gas.prm's own header. 2500 steps at
# 4236462's own measured ~1.03s/step (72000 steps in 1h14m19s) is ~43min -- 1h requested for margin.
#
# 140 ranks (matches 4236462's own decomposition exactly -- the checkpoint's blockForest.file is
# tied to that exact rank/block layout). Deliberately does NOT touch snapshots/ (unlike the original
# run's own sbatch, which removes leftover snapshots/ at the start) -- this job LOADS that checkpoint
# and (storeSnapshot=false in the diag prm) must not overwrite or delete it.
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

for dir in vtk-out eval; do
    if [[ -d "$BIN_DIR/$dir" ]]; then
        echo "WARNING: removing leftover $BIN_DIR/$dir from an unarchived previous run" >&2
        rm -rf "${BIN_DIR:?}/$dir"
    fi
done
# NOTE: deliberately NOT touching snapshots/ -- this job loads job 4236462's own t=70000 checkpoint.

"$SH" --build --jobs "$(nproc)"

echo "############################################################"
echo "## cuboid stillwater hires d025 -- DIAGNOSTIC resume from t=70000, dense NaN scan 70000-72500"
echo "############################################################"
"$SH" --run --np 140 --prm "$PRM_DIR/cuboid_stillwater_realgravity981_3xsphere_hires_d025_diag_gas.prm"

echo "## DIAG finished -- check output for [NaN-SCAN]/[LOW-RHO-PRECURSOR]/[HIGH-VELOCITY-PRECURSOR] onset."
