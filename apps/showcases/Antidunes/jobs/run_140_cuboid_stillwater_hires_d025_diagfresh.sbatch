#!/bin/bash -l
#SBATCH --job-name=antidunes-cuboid-stillwater-hires-d025-diagfresh
#SBATCH --partition=multinode
#SBATCH --nodes=2
#SBATCH --ntasks=140
#SBATCH --time=01:30:00
#SBATCH --output=sbatch-output/output/cuboid_stillwater_hires_d025_diagfresh_%j.out
#SBATCH --error=sbatch-output/error/cuboid_stillwater_hires_d025_diagfresh_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# Thesis item 1a (THESIS_PLAN.md) -- fresh-run diagnostic for the R3 equilibrium NaN crash (job
# 4237630 FAILED at t=6000/100000, ~84 FATAL cells all at local (0,0,0) -- same first-per-block
# reporting artifact as the earlier checkpoint-resume investigation). No checkpoint exists yet
# (crash is before the first scheduled snapshot at t=10000), so this reruns from scratch with
# denseNanScanEnabled bracketing 5500-6500 tightly (meanVelocity was smooth through the last logged
# point, t=5500, so the onset is somewhere in that 500-step window) plus low-rho/high-velocity
# precursor checks and listAllBadCellsNoEarlyExit (temporary diagnostic in SetupStages.h) to get the
# TRUE extent of the corruption instead of one cell per block.
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

"$SH" --build --jobs "$(nproc)"

echo "############################################################"
echo "## cuboid stillwater hires d025 -- FRESH diagnostic, dense NaN scan 5500-6500"
echo "############################################################"
"$SH" --run --np 140 --prm "$PRM_DIR/cuboid_stillwater_realgravity981_3xsphere_hires_d025_diagfresh_gas.prm"

echo "## DIAGFRESH finished -- check output for [FULL-BAD-CELL-LIST]/[LOW-RHO-PRECURSOR]/[HIGH-VELOCITY-PRECURSOR]."
