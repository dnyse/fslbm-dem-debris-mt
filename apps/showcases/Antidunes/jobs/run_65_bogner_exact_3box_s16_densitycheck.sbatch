#!/bin/bash -l
#SBATCH --job-name=antidunes-bogner-exact-s16-densitycheck
#SBATCH --partition=singlenode
#SBATCH --nodes=1
#SBATCH --ntasks=65
#SBATCH --time=2:00:00
#SBATCH --output=sbatch-output/output/bogner_exact_3box_s16_densitycheck_%j.out
#SBATCH --error=sbatch-output/error/bogner_exact_3box_s16_densitycheck_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# [Claude, 2026-09-24] Two follow-ups to the torque-sweep equilibrium-shift investigation (THESIS_PLAN.md
# "MECHANISM of the free-sweep's 3-5 deg overshoot" entry): the fully-free omega1p9_seed5 run (rho=0.75, 5
# deg seed) converges to ~32deg instead of Bogner's 26.565, while an X/Y+pitch/yaw-pinned control at the
# SAME rho/seed correctly holds at 27.0deg. Geometric back-calculation showed the discrepancy is
# equivalent to the box behaving as if its effective density were ~0.734 instead of the configured 0.75
# (~2.2% low), with the vertical force balance itself confirmed essentially exact (buoyancy cancels
# gravity to 5 sig figs) -- i.e. B is displaced sideways, not under/over-strength.
#
# Run 1 (densitycheck_rho075): BYTE-IDENTICAL replica of omega1p9_seed5 (same rho=0.75, same 5deg seed,
# same everything) with "density;" added to the fluid_field VTK writers list (not present before) --
# tests whether the LOCAL fluid density near the box has drifted persistently above the nominal 1.0
# reference (the mechanism that would explain an apparent buoyancy-strength deficit without the box's
# own configured mass changing).
#
# Run 2 (densitycheck_rho0767): same everything, but densityRatio bumped from 0.75 to 0.7665 (+2.2%,
# compensating the effective deficit) -- if the deficit is roughly density-independent, this should
# settle near the CORRECT 26.565deg instead of ~32deg. Also has density; VTK output for a second,
# independent look at the local density field.
#
# Both prms already carry dataSliceFrequency=500 (the omega1p9_seed5 prm they were copied from predates
# that fix and was the one that caused job 4255656's original vault file-count quota failure -- added
# explicitly here to not repeat it).
#
# Timing: prior 100k-step s16 single-box-family runs took ~44 min actual wall time (last2 job's own
# observed cost) + ~5-7 min archive overhead each; two runs -> ~110 min; 2h limit for margin.
# Queue: independent, ~/mt-dennys-mirror2 (confirmed free before submission).
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
    ARCHIVE_ROOT="$BIN_DIR/archive/${STAMP}_bogner_exact_s16_densitycheck"
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
    echo "## Bogner exact box equilibrium, s=16 density-check, $tag"
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

run_one densitycheck_rho075 "$PRM_DIR/bogner_exact_3box_s16_seed5_rho075_nopin_omega1p9_100k_densitycheck_gas.prm"
run_one densitycheck_rho0767 "$PRM_DIR/bogner_exact_3box_s16_seed5_rho0767_nopin_omega1p9_100k_densitycheck_gas.prm"

echo "## ALL TWO DONE. Archive root: $ARCHIVE_ROOT"
exit "$OVERALL_RC"
