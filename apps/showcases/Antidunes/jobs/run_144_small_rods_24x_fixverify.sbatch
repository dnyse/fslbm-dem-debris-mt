#!/bin/bash -l
#SBATCH --job-name=antidunes-small-rods-fixverify
#SBATCH --partition=multinode
#SBATCH --nodes=2
#SBATCH --ntasks=144
#SBATCH --time=01:00:00
#SBATCH --output=sbatch-output/output/small_rods_24x_fixverify_%j.out
#SBATCH --error=sbatch-output/error/small_rods_24x_fixverify_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# [Claude, 2026-09-27] Per explicit request ("do it"): verifies the MovingParticleMapping ownership-
# handoff fix against the same step-20001 checkpoint used throughout today's investigation, freeze
# OFF, 25000-step budget (comfortably past every crash step seen today: 22012-22257). Checkpoints
# every 500 steps with keepLastTwoSnapshots, into a fresh dedicated folder that never touches the
# original backup.
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

RESCUE="$BIN_DIR/archive/rescued_before_${SLURM_JOB_ID:-manual}_$(date +%Y%m%d_%H%M%S)"
for dir in vtk-out eval snapshots; do
    if [[ -d "$BIN_DIR/$dir" ]]; then
        mkdir -p "$RESCUE"
        echo "WARNING: leftover $BIN_DIR/$dir from a previous run -- moving to $RESCUE (NOT deleting)" >&2
        mv "$BIN_DIR/$dir" "$RESCUE/$dir"
    fi
done

"$SH" --build --jobs "$(nproc)"

echo "############################################################"
echo "## FIX VERIFICATION: MovingParticleMapping ownership-handoff fix, freeze OFF, from step 20001"
echo "## toward step 25000 -- comfortably past every crash step observed today."
echo "############################################################"
set +e
"$SH" --run --np 144 --prm "$PRM_DIR/small_rods_24x_2x4x3_clog_resume_fixverify_gas.prm"
RUN_RC=$?
set -e
echo "## run exited with status $RUN_RC -- archiving whatever output exists"

STAMP="$(date +%Y%m%d_%H%M%S)"
VAULT_CANDIDATE="${HPCVAULT:-}/antidunes-archive/${SLURM_JOB_NAME}_${SLURM_JOB_ID}_${STAMP}"
if [[ -n "${HPCVAULT:-}" ]] && mkdir -p "$VAULT_CANDIDATE" 2>/dev/null \
   && touch "$VAULT_CANDIDATE/.write_test" 2>/dev/null; then
    rm -f "$VAULT_CANDIDATE/.write_test"
    ARCHIVE_ROOT="$VAULT_CANDIDATE"
else
    ARCHIVE_ROOT="$BIN_DIR/archive/${STAMP}_small_rods_24x_fixverify"
    mkdir -p "$ARCHIVE_ROOT"
    echo "WARNING: \$HPCVAULT unavailable at archive time -- falling back to LOCAL archive at $ARCHIVE_ROOT (move this off manually to avoid quota pressure)" >&2
fi
for d in vtk-out eval; do
    [[ -d "$BIN_DIR/$d" ]] && mv "$BIN_DIR/$d" "$ARCHIVE_ROOT/$d"
done
# NOTE: snapshots are NOT moved -- this run's checkpoints (resume_snapshots/small_rods_24x_fixverify_run/
# gen_0|gen_1) are left in place for a future resume/debug, not archived away.

echo "## DONE. Archived under $ARCHIVE_ROOT (checkpoints remain in resume_snapshots/small_rods_24x_fixverify_run/)"
exit "$RUN_RC"
