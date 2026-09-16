#!/bin/bash -l
#SBATCH --job-name=antidunes-clog-debug-extreme-pidon-phase1
#SBATCH --partition=multinode
#SBATCH --nodes=2
#SBATCH --ntasks=144
#SBATCH --time=12:00:00
#SBATCH --output=sbatch-output/output/clog_debug_extreme_pidon_phase1_%j.out
#SBATCH --error=sbatch-output/error/clog_debug_extreme_pidon_phase1_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# Thesis item 3 (THESIS_PLAN.md) -- PHASE 1 of a two-phase debug rerun of job 4234953's
# `reynoldsNumber is nan!` crash (traced to a [DRY-LINK]/[LOW-RHO] cell at global ~(376,25,14),
# right against pillar 2, first appearing at step 48819). Single-case job, NOT the 3-case
# run_case() combined pattern -- this only needs to reproduce ONE case (extreme, PID on) with
# checkpointing enabled, to produce a snapshot near the crash for PHASE 2
# (run_144_clog_debug_extreme_pidon_phase2.sbatch) to resume from with full diagnostics + tight VTK.
#
# clog_debug_extreme_pidon_phase1_gas.prm is otherwise identical to
# clog_measurement_validation_extreme_pidon_gas.prm except storeSnapshot=true
# (snapshotFrequency=500) and timesteps capped at 50000 (a safety margin past the expected
# crash, not the original 97500). Expected outcome: this job crashes with the same
# `reynoldsNumber is nan!` abort as job 4234953, near step 48819 -- leaving behind a checkpoint
# from step 48500 (the last multiple of 500 before that) in snapshots/, ready for phase 2.
#
# 144 ranks / 2 nodes -- same resource shape as the real extreme scenario (cellsPerBlock
# <64,20,15>, numBlocks (768/64, 40/20, 90/15) = 144).
#
# Time limit 12:00:00: the original crashed run (job 4234953) took ~1h13m of wall-clock to reach
# step 48819 including the "no-clog reference" case that ran first in that combined job; this is a
# single-case job (extreme_pidon only), so expect well under that -- 12h is a generous ceiling, not
# an estimate.
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
# NOTE: deliberately NOT clearing snapshots/ here -- this job's whole point is to WRITE one there
# for phase 2 to consume afterward. A stale snapshot from an unrelated earlier run would still be
# overwritten cleanly by this job's own first checkpoint (storeSnapshot always replaces, never
# appends), so no explicit pre-clear is needed.

"$SH" --build --jobs "$(nproc)"

echo "############################################################"
echo "## clog debug: extreme, PID on -- PHASE 1 (checkpointed reproduction)"
echo "############################################################"
"$SH" --run --np 144 --prm "$PRM_DIR/clog_debug_extreme_pidon_phase1_gas.prm"

echo "## PHASE 1 finished without crashing (unexpected -- see this script's own header). Check eval/ for results."
