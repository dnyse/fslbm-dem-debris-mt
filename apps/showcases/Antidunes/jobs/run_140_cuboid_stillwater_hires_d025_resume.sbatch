#!/bin/bash -l
#SBATCH --job-name=antidunes-cuboid-stillwater-hires-d025-resume
#SBATCH --partition=multinode
#SBATCH --nodes=2
#SBATCH --exclude=f0220,f0258
#SBATCH --ntasks=140
#SBATCH --time=08:00:00
#SBATCH --output=sbatch-output/output/cuboid_stillwater_hires_d025_resume_%j.out
#SBATCH --error=sbatch-output/error/cuboid_stillwater_hires_d025_resume_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# Thesis item 1a (THESIS_PLAN.md) -- RESUME #2 of the R3 higher-resolution/longer-duration cuboid
# equilibrium rerun (densityRatio=0.25). First resume attempt (job 4238851) reached step 30000/100000
# at ~469ms/step (matching the pre-fix box-MEM rate exactly) before being CANCELLED deliberately --
# not a crash -- specifically to restart with the just-applied logging fix (removed the
# addIgnoreRegex-based suppression AND commented out the underlying noisy per-cell framework log
# statements in ObstacleNormalSweep.impl.h/ObstacleFillLevelSweep.impl.h/PdfReconstructionManager.h/
# CurvatureSweep.impl.h -- see git history/session notes). That fix measured an ~8.3x per-step
# speedup on the R2 box-MEM comparison (440ms/step -> 53ms/step, job 4238846 vs 4239167). Resuming
# again from the step-30000 checkpoint (freshly written just before cancellation) via
# loadSnapshot=true (already set in the prm).
#
# Unlike the fresh-run sbatch, this script does NOT delete an existing snapshots/ directory on
# startup -- that directory IS the checkpoint being resumed from. Only vtk-out/eval leftovers (not
# needed to resume) are cleaned.
#
# Time budget: remaining 70000 steps at the OLD 469ms/step rate would need ~9.1h; even a conservative
# 3x (not the full ~8.3x measured on box-MEM) speedup brings that under 3.1h. Requested 8h for
# margin without requesting an excessive allocation.
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
        echo "WARNING: removing leftover $BIN_DIR/$dir from the timed-out previous attempt" >&2
        rm -rf "${BIN_DIR:?}/$dir"
    fi
done

if [[ ! -f "$BIN_DIR/snapshots/simState.file" ]]; then
    echo "FATAL: no snapshots/simState.file found in $BIN_DIR -- refusing to run with loadSnapshot=true" >&2
    exit 1
fi

STAMP="$(date +%Y%m%d_%H%M%S)"
if [[ -n "${HPCVAULT:-}" ]] && mkdir -p "$HPCVAULT/antidunes-archive/${SLURM_JOB_NAME}_${SLURM_JOB_ID}_$STAMP" 2>/dev/null; then
    ARCHIVE_ROOT="$HPCVAULT/antidunes-archive/${SLURM_JOB_NAME}_${SLURM_JOB_ID}_$STAMP"
else
    ARCHIVE_ROOT="$BIN_DIR/archive/$STAMP"
    mkdir -p "$ARCHIVE_ROOT"
fi

"$SH" --build --jobs "$(nproc)"

echo "############################################################"
echo "## cuboid stillwater hires, densityRatio=0.25 -- RESUME from checkpoint"
echo "############################################################"
"$SH" --run --np 140 --prm "$PRM_DIR/cuboid_stillwater_realgravity981_3xsphere_hires_d025_gas.prm"

for d in vtk-out eval snapshots; do
    [[ -d "$BIN_DIR/$d" ]] && mkdir -p "$ARCHIVE_ROOT/d025" && mv "$BIN_DIR/$d" "$ARCHIVE_ROOT/d025/$d"
done

echo "## DONE. Archived under $ARCHIVE_ROOT/d025/{eval,vtk-out,snapshots}"
