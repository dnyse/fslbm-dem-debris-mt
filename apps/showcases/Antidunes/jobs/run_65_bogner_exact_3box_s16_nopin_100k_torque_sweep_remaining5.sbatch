#!/bin/bash -l
#SBATCH --job-name=antidunes-bogner-exact-s16-nopin-100k-torque-remaining5
#SBATCH --partition=singlenode
#SBATCH --nodes=1
#SBATCH --ntasks=65
#SBATCH --time=2:30:00
#SBATCH --output=sbatch-output/output/bogner_exact_3box_s16_nopin_100k_torque_remaining5_%j.out
#SBATCH --error=sbatch-output/error/bogner_exact_3box_s16_nopin_100k_torque_remaining5_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# [Claude, 2026-09-22] RERUN of the 5 sub-runs from job 4255656 that never ran (that job's 6-run combined sweep completed only its FIRST
# sub-run, omega1p9_seed5, then died archiving it: AverageDataSliceEvaluator was sharing evaluationFrequency with [BoxForces]/[BoxTorque]
# logging, so --evalfreq 1 also made it dump 3 text files/step -- ~300000 files over 100000 steps blew a vault disk quota mid-archive).
# FIX (SimulationParams.h/BuildTimeloop.h): new, separate EvaluationParameters.dataSliceFrequency (defaults to evaluationFrequency itself
# if unset, so nothing else changes); gen_bogner_exact_prm.py now sets it to 500 explicitly regardless of --evalfreq. Verified locally: a
# 600-step smoke test gave exactly 600 [BoxTorque] lines (genuine per-step) and only 6 slice files (properly throttled), rc=0.
# omega1p9_seed5 itself is NOT rerun here -- it already completed successfully in 4255656 and was copied to
# animations/4255656_omega1p9_seed5_100k_torque_20260922/ (vtk-out + eval/antidunes.sqlite + the full per-step BoxForces/BoxTorque log).
# 5 runs, all rho=0.75, fully free (nopin), 100000 steps, evaluationFrequency 1 / dataSliceFrequency 500:
#   omega 1.5/seed5, omega 1.0/seed5, omega 1.9/seed3, omega 1.9/seed4, omega 1.9/seedEq (free start AT the equilibrium heel)
# Time budget: 4255587's own 3 runs at 50000 steps took 34m58s total (~11.6 min/run); 100000 steps ~2x that -> ~23 min/run x 5 = ~115 min
# + build/archive overhead (now minimal, no more huge file dumps) -- 2h30m limit for margin.
# Queue: independent, ~/mt-dennys-mirror2 (free).
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
    ARCHIVE_ROOT="$BIN_DIR/archive/${STAMP}_bogner_exact_s16_nopin_100k_torque_remaining5"
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
    echo "## Bogner exact box equilibrium, s=16 nopin 100k torque (remaining5), $tag"
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

run_one omega1p5_seed5 "$PRM_DIR/bogner_exact_3box_s16_seed5_rho075_nopin_omega1p5_100k_torque_gas.prm"
run_one omega1p0_seed5 "$PRM_DIR/bogner_exact_3box_s16_seed5_rho075_nopin_omega1p0_100k_torque_gas.prm"
run_one omega1p9_seed3 "$PRM_DIR/bogner_exact_3box_s16_seed3_rho075_nopin_omega1p9_100k_torque_gas.prm"
run_one omega1p9_seed4 "$PRM_DIR/bogner_exact_3box_s16_seed4_rho075_nopin_omega1p9_100k_torque_gas.prm"
run_one omega1p9_seedEq "$PRM_DIR/bogner_exact_3box_s16_seedEq_rho075_nopin_omega1p9_100k_torque_gas.prm"

echo "## ALL FIVE DONE. Archive root: $ARCHIVE_ROOT"
exit "$OVERALL_RC"
