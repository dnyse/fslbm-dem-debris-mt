#!/bin/bash -l
#SBATCH --job-name=antidunes-clog-debug-extreme-pidon-phase1-5
#SBATCH --partition=multinode
#SBATCH --nodes=2
#SBATCH --ntasks=144
#SBATCH --time=02:30:00
#SBATCH --output=sbatch-output/output/clog_debug_extreme_pidon_phase1_5_%j.out
#SBATCH --error=sbatch-output/error/clog_debug_extreme_pidon_phase1_5_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# Thesis item 3 (THESIS_PLAN.md) -- PHASE 1.5: a plain checkpoint-resume smoke test, run after
# PHASE 1 (job 4235135) completed all 50000 of its own steps WITHOUT reproducing job 4234953's
# `reynoldsNumber is nan!` crash (unexpected -- see clog_debug_extreme_pidon_phase1_gas.prm's own
# header). With no crash to dense-scan, PHASE 2 has nothing to target yet. This job instead just
# resumes PHASE 1's own final checkpoint (step 50000, in snapshots/) and continues for the SAME
# DURATION phase 1 ran (another 50000 steps -> timesteps=100000 in the prm, since `timesteps` is
# the absolute target step), with NONE of phase 2's diagnostic overhead (no denseNanScanEnabled, no
# precursor thresholds, no tight VTK) -- just confirming loadSnapshot resumes correctly and the run
# continues physically sensibly, while also giving a second independent window (steps 50000-100000)
# in which the original crash could, or could not, reappear.
#
# 144 ranks / 2 nodes -- same resource shape as phase 1.
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
# NOTE: deliberately NOT touching snapshots/ here -- this job's whole point is to LOAD the
# checkpoint phase 1 left there, and (with storeSnapshot also on) to keep overwriting it as this
# run progresses past step 50000.

"$SH" --build --jobs "$(nproc)"

echo "############################################################"
echo "## clog debug: extreme, PID on -- PHASE 1.5 (checkpoint-resume smoke test)"
echo "############################################################"
"$SH" --run --np 144 --prm "$PRM_DIR/clog_debug_extreme_pidon_phase1_5_gas.prm"

echo "## PHASE 1.5 finished -- check eval/ for whether the run continued sensibly past step 50000."
