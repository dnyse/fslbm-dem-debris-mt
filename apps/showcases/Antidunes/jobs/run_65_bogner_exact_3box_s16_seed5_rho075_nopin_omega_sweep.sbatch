#!/bin/bash -l
#SBATCH --job-name=antidunes-bogner-exact-s16-seed5-rho075-nopin-omega-sweep
#SBATCH --partition=singlenode
#SBATCH --nodes=1
#SBATCH --ntasks=65
#SBATCH --time=2:00:00
#SBATCH --output=sbatch-output/output/bogner_exact_3box_s16_seed5_rho075_nopin_omega_sweep_%j.out
#SBATCH --error=sbatch-output/error/bogner_exact_3box_s16_seed5_rho075_nopin_omega_sweep_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# [Claude, 2026-09-22] Thesis item 1a (THESIS_PLAN.md): VISCOSITY/DAMPING TEST, rerun of the seed5/rho075 omega{1.9,1.5,1.0} sweep, this time with the box FULLY FREE (--nopin: 6 DOF,
# no position/rotation pins) instead of the position-pinned/roll-only Bogner-exact rig the first attempt (4255530/31/32) used. That first attempt showed NO overshoot at all at any
# damping level -- all three stayed locked near the 5 deg seed for the full 50000 steps, because the PINNED rig has a basin around upright that a 5 deg seed apparently can't escape,
# regardless of damping. The actual "overshoots past 26.565, settles/wanders at 30-37 deg" behaviour this experiment is meant to explain was only ever seen in the FULLY FREE family
# (the s=30 cuboid_stillwater sweep, and the s=16/hires 'nopin' 300k continuations) -- this rerun matches that setup instead.
#
# All THREE omega values (1.9 Bogner's own value / 1.5 / 1.0, most damping) run SEQUENTIALLY in this ONE job/ONE sbatch file (explicit instruction), each to its own archive
# subfolder, build done once up front. A crashed/blown-up run does not stop the sweep -- each run's exit status is captured, its output is archived regardless, and the NEXT omega
# still runs; the job's own exit code is non-zero if ANY of the three failed. 65 ranks, 50000 steps each; the pinned siblings took 11-13 min each, so budget similarly per run here
# (untested: a fully free box could behave very differently, e.g. drift/rotate into a slower- or faster-converging state) -- 2h limit for headroom across all three plus build.
# See prms/bogner_exact_3box_s16_seed5_rho075_nopin_omega{1p9,1p5,1p0}_gas.prm and gen_bogner_exact_prm.py's --nopin/--omega options.
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

# Prefer $HPCVAULT, decided and verified HERE (right before use), one archive root shared by all three sub-runs of this job.
STAMP="$(date +%Y%m%d_%H%M%S)"
VAULT_CANDIDATE="${HPCVAULT:-}/antidunes-archive/${SLURM_JOB_NAME}_${SLURM_JOB_ID}_${STAMP}"
if [[ -n "${HPCVAULT:-}" ]] && mkdir -p "$VAULT_CANDIDATE" 2>/dev/null \
   && touch "$VAULT_CANDIDATE/.write_test" 2>/dev/null; then
    rm -f "$VAULT_CANDIDATE/.write_test"
    ARCHIVE_ROOT="$VAULT_CANDIDATE"
else
    ARCHIVE_ROOT="$BIN_DIR/archive/${STAMP}_bogner_exact_s16_seed5_nopin_omega_sweep"
    mkdir -p "$ARCHIVE_ROOT"
    echo "WARNING: \$HPCVAULT unavailable at archive time -- falling back to LOCAL archive at $ARCHIVE_ROOT (move this off manually to avoid quota pressure)" >&2
fi

OVERALL_RC=0

run_one() {
    local tag="$1" prm="$2"

    # Never DELETE a previous sub-run's output -- move it aside instead (same filesystem -> instant rename).
    local RESCUE="$BIN_DIR/archive/rescued_before_${tag}_${SLURM_JOB_ID:-manual}_$(date +%Y%m%d_%H%M%S)"
    for dir in vtk-out eval snapshots; do
        if [[ -d "$BIN_DIR/$dir" ]]; then
            mkdir -p "$RESCUE"
            echo "WARNING: leftover $BIN_DIR/$dir before $tag -- moving to $RESCUE (NOT deleting)" >&2
            mv "$BIN_DIR/$dir" "$RESCUE/$dir"
        fi
    done

    echo "############################################################"
    echo "## Bogner exact 3-box equilibrium test, s=16 seed5 nopin, $tag"
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

run_one omega1p9 "$PRM_DIR/bogner_exact_3box_s16_seed5_rho075_nopin_omega1p9_gas.prm"
run_one omega1p5 "$PRM_DIR/bogner_exact_3box_s16_seed5_rho075_nopin_omega1p5_gas.prm"
run_one omega1p0 "$PRM_DIR/bogner_exact_3box_s16_seed5_rho075_nopin_omega1p0_gas.prm"

echo "## ALL THREE DONE. Archive root: $ARCHIVE_ROOT"
exit "$OVERALL_RC"
