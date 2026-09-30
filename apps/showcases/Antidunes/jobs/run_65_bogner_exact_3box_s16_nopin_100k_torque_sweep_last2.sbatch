#!/bin/bash -l
#SBATCH --job-name=antidunes-bogner-exact-s16-nopin-100k-torque-last2
#SBATCH --partition=singlenode
#SBATCH --nodes=1
#SBATCH --ntasks=65
#SBATCH --time=2:30:00
#SBATCH --output=sbatch-output/output/bogner_exact_3box_s16_nopin_100k_torque_last2_%j.out
#SBATCH --error=sbatch-output/error/bogner_exact_3box_s16_nopin_100k_torque_last2_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# [Claude, 2026-09-23] RERUN of the LAST 2 sub-runs from job 4256467 (the "remaining5" job), which itself
# TIMED OUT (2h30m limit) after completing only 3 of its 5: omega1p5_seed5, omega1p0_seed5, omega1p9_seed3
# all finished cleanly (exit 0, archived under
# /home/vault/iwia/iwia148h/antidunes-archive/antidunes-bogner-exact-s16-nopin-100k-torque-remaining5_4256467_.../
# {omega1p5_seed5,omega1p0_seed5,omega1p9_seed3}/), but omega1p9_seed4 was killed mid-run (no valid output)
# and omega1p9_seedEq never started.
#
# ROOT CAUSE of the timeout: the 2h30m budget was sized from an earlier estimate of ~23 min/run (scaling
# 4255587's 50000-step runs linearly to 100000 steps), but each 100k-step run's OWN internal wall-clock
# (per its END LOGGING line) actually took ~2629s = ~43.8 min -- roughly 2x the estimate -- plus real
# per-run archive/rescue overhead (~6-7 min each, from mv-to-vault of vtk-out/eval/snapshots). 3 full runs
# alone (~131 min) + build + 3x archiving already consumed essentially the whole 150-minute budget.
#
# This job re-sizes the budget from that OBSERVED per-run cost instead: 2 runs x ~44 min + 2 x ~10 min
# archive overhead (rounded up for margin) + ~2 min build = ~110 min -> 2h30m limit kept for headroom.
#
# 2 runs, both rho=0.75, fully free (nopin), 100000 steps, evaluationFrequency 1 / dataSliceFrequency 500:
#   omega 1.9/seed4, omega 1.9/seedEq (free start AT the equilibrium heel)
# Queue: independent, ~/mt-dennys-mirror2 (free -- 4256467 already ended there, one way or another).
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

"$SH" --build --jobs "$(nproc)"

STAMP="$(date +%Y%m%d_%H%M%S)"
VAULT_CANDIDATE="${HPCVAULT:-}/antidunes-archive/${SLURM_JOB_NAME}_${SLURM_JOB_ID}_${STAMP}"
if [[ -n "${HPCVAULT:-}" ]] && mkdir -p "$VAULT_CANDIDATE" 2>/dev/null \
   && touch "$VAULT_CANDIDATE/.write_test" 2>/dev/null; then
    rm -f "$VAULT_CANDIDATE/.write_test"
    ARCHIVE_ROOT="$VAULT_CANDIDATE"
else
    ARCHIVE_ROOT="$BIN_DIR/archive/${STAMP}_bogner_exact_s16_nopin_100k_torque_last2"
    mkdir -p "$ARCHIVE_ROOT"
    echo "WARNING: \$HPCVAULT unavailable at archive time -- falling back to LOCAL archive at $ARCHIVE_ROOT (move this off manually to avoid quota pressure)" >&2
fi

OVERALL_RC=0

run_one() {
    local tag="$1" prm="$2"

    local RESCUE="$BIN_DIR/archive/rescued_before_${tag}_${SLURM_JOB_ID:-manual}_$(date +%Y%m%d_%H%M%S)"
    for dir in vtk-out eval snapshots; do
        if [[ -d "$BIN_DIR/$dir" ]]; then
            mkdir -p "$RESCUE"
            echo "WARNING: leftover $BIN_DIR/$dir before $tag -- moving to $RESCUE (NOT deleting)" >&2
            mv "$BIN_DIR/$dir" "$RESCUE/$dir"
        fi
    done

    echo "############################################################"
    echo "## Bogner exact box equilibrium, s=16 nopin 100k torque (last2), $tag"
    echo "############################################################"
    set +e
    "$SH" --run --np 65 --prm "$prm"
    local rc=$?
    set -e
    echo "## $tag exited with status $rc -- archiving whatever output exists"
    [[ $rc -ne 0 ]] && OVERALL_RC=1

    for d in vtk-out eval snapshots; do
        [[ -d "$BIN_DIR/$d" ]] && mkdir -p "$ARCHIVE_ROOT/$tag" && mv "$BIN_DIR/$d" "$ARCHIVE_ROOT/$tag/$d"
    done
    echo "## DONE $tag. Archived under $ARCHIVE_ROOT/$tag/{eval,vtk-out,snapshots}"
}

run_one omega1p9_seed4 "$PRM_DIR/bogner_exact_3box_s16_seed4_rho075_nopin_omega1p9_100k_torque_gas.prm"
run_one omega1p9_seedEq "$PRM_DIR/bogner_exact_3box_s16_seedEq_rho075_nopin_omega1p9_100k_torque_gas.prm"

echo "## ALL TWO DONE. Archive root: $ARCHIVE_ROOT"
exit "$OVERALL_RC"
