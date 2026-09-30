#!/bin/bash -l
#SBATCH --job-name=antidunes-bogner-exact-s16-rho0767-tilt26p565-allsnaps
#SBATCH --partition=singlenode
#SBATCH --nodes=1
#SBATCH --ntasks=65
#SBATCH --time=0:30:00
#SBATCH --output=sbatch-output/output/bogner_exact_3box_s16_rho0767_tilt26p565_allsnaps_%j.out
#SBATCH --error=sbatch-output/error/bogner_exact_3box_s16_rho0767_tilt26p565_allsnaps_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# [Claude, 2026-09-24] Third attempt at the 26.565deg-seed rho=0.7665 case (jobs 4260975, 4261005 both
# FAILED identically: reynoldsNumber nan at exactly step 444, a fully deterministic crash -- see
# THESIS_PLAN.md and this run's own prm header for the full trace). This run adds `keepAllSnapshots
# true` / `keepAllSnapshotsFromStep 400` (new BlockForestParameters fields) so every step from 400
# onward is checkpointed into its own snapshots/step_<N>/ subfolder instead of being overwritten,
# letting us load and inspect the pdf/fill/excess-mass fields cell-by-cell at every step leading up to
# the crash. timesteps cut to 1000 (crash is deterministic at ~444; no need to run further). Expect a
# crash around step 444 again -- that IS the expected/successful outcome for this diagnostic run, not
# a failure to react to. ~45-50 checkpoints of the full 130x40x40 domain expected; small enough to keep
# all of them (no rescue-before-delete needed here -- this run's checkpoints ARE the point).
# Queue: independent, confirmed clean/free before submission.
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

"$SH" --build --jobs "$(nproc)"

echo "############################################################"
echo "## Bogner exact box equilibrium, s=16, rho=0.7665, tilt=26.565deg, all-snapshots diagnostic"
echo "############################################################"
set +e
"$SH" --run --np 65 --prm "$PRM_DIR/bogner_exact_3box_s16_seed5_rho0767_nopin_omega1p9_tilt26p565_allsnaps_gas.prm"
RC=$?
set -e
echo "## run exited with status $RC (a crash around step 444 is the EXPECTED outcome here -- the point is the per-step checkpoints leading up to it, not a clean finish)"

STAMP="$(date +%Y%m%d_%H%M%S)"
VAULT_CANDIDATE="${HPCVAULT:-}/antidunes-archive/${SLURM_JOB_NAME}_${SLURM_JOB_ID}_${STAMP}"
if [[ -n "${HPCVAULT:-}" ]] && mkdir -p "$VAULT_CANDIDATE" 2>/dev/null \
   && touch "$VAULT_CANDIDATE/.write_test" 2>/dev/null; then
    rm -f "$VAULT_CANDIDATE/.write_test"
    ARCHIVE_ROOT="$VAULT_CANDIDATE"
else
    ARCHIVE_ROOT="$BIN_DIR/archive/${STAMP}_rho0767_tilt26p565_allsnaps"
    mkdir -p "$ARCHIVE_ROOT"
    echo "WARNING: \$HPCVAULT unavailable at archive time -- falling back to LOCAL archive at $ARCHIVE_ROOT (move this off manually to avoid quota pressure)" >&2
fi
for d in vtk-out eval snapshots; do
    [[ -d "$BIN_DIR/$d" ]] && mv "$BIN_DIR/$d" "$ARCHIVE_ROOT/$d"
done

echo "## DONE. Archived under $ARCHIVE_ROOT"
exit 0
