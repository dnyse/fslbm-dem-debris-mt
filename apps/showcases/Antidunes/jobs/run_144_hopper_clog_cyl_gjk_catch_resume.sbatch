#!/bin/bash -l
#SBATCH --job-name=antidunes-hopper-clog-cyl-gjk-catch-resume
#SBATCH --partition=multinode
#SBATCH --nodes=2
#SBATCH --ntasks=144
#SBATCH --time=03:00:00
#SBATCH --output=sbatch-output/output/hopper_clog_cyl_gjk_catch_resume_%j.out
#SBATCH --error=sbatch-output/error/hopper_clog_cyl_gjk_catch_resume_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# [Claude, 2026-09-20] RESUME + FLOATER-HISTORY rerun of job 4250846 (GJK-only cylinder hopper 4251066, aborted at step 6719): restores that job's step-5001 checkpoint
# from the vault, re-runs to step 6000 with logFloaterHistory on for steps 5650..5850 (see prms/hopper_clog_cyl_gjk_catch_resume_gas.prm). Same binary/code as the failed run.
# Restored checkpoint is COPIED (the vault archive stays untouched). Output (vtk-out, eval, floater-history) is archived even if the run aborts again (expected).
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

# [Claude, 2026-09-18] Never DELETE a previous run's output. These jobs are chained in ONE checkout with --dependency, so a
# predecessor that crashed or hit its time limit leaves vtk-out/eval/snapshots here -- the old `rm -rf` would have destroyed it.
# Move them aside instead (same filesystem -> instant rename); archive/rescued_* is then moved off by hand.
RESCUE="$BIN_DIR/archive/rescued_before_${SLURM_JOB_ID:-manual}_$(date +%Y%m%d_%H%M%S)"
for dir in vtk-out eval snapshots floater-history; do
    if [[ -d "$BIN_DIR/$dir" ]]; then
        mkdir -p "$RESCUE"
        echo "WARNING: leftover $BIN_DIR/$dir from a previous run -- moving to $RESCUE (NOT deleting)" >&2
        mv "$BIN_DIR/$dir" "$RESCUE/$dir"
    fi
done

CKPT="${HPCVAULT:-/nonexistent}/antidunes-archive/antidunes-hopper-clog-cyl-frontblock-pidon-2node-gjk_4251066_20260920_143250/snapshots"
if [[ ! -f "$CKPT/simState.file" ]]; then echo "FATAL: checkpoint $CKPT/simState.file not found" >&2; exit 1; fi
mkdir -p "$BIN_DIR/snapshots"
cp -a "$CKPT/." "$BIN_DIR/snapshots/"
echo "## restored checkpoint: simState = $(cat "$BIN_DIR/snapshots/simState.file")  (expect 5001)"
du -sh "$BIN_DIR/snapshots"

"$SH" --build --jobs "$(nproc)"

echo "############################################################"
echo "## hopper cylinders front-block RESUME from step 5001 with floater-history (steps 5650..5850)"
echo "############################################################"
# a crashed run must still reach the archive step below (set -e would skip it), so capture the exit status
set +e
"$SH" --run --np 144 --prm "$PRM_DIR/hopper_clog_cyl_gjk_catch_resume_gas.prm"
RUN_RC=$?
set -e
echo "## run exited with status $RUN_RC -- archiving whatever output exists"

# Prefer $HPCVAULT, decided and verified HERE (right before use) -- see the 288-rank sibling job's
# own 2026-09-16 comment for the full rationale (autofs staleness, local-only is a last resort).
STAMP="$(date +%Y%m%d_%H%M%S)"
VAULT_CANDIDATE="${HPCVAULT:-}/antidunes-archive/${SLURM_JOB_NAME}_${SLURM_JOB_ID}_${STAMP}"
if [[ -n "${HPCVAULT:-}" ]] && mkdir -p "$VAULT_CANDIDATE" 2>/dev/null \
   && touch "$VAULT_CANDIDATE/.write_test" 2>/dev/null; then
    rm -f "$VAULT_CANDIDATE/.write_test"
    ARCHIVE_ROOT="$VAULT_CANDIDATE"
else
    ARCHIVE_ROOT="$BIN_DIR/archive/${STAMP}_hopper_clog_cyl_gjk_resume5000_fluidscan"
    mkdir -p "$ARCHIVE_ROOT"
    echo "WARNING: \$HPCVAULT unavailable at archive time -- falling back to LOCAL archive at $ARCHIVE_ROOT (move this off manually to avoid quota pressure)" >&2
fi
for d in vtk-out eval snapshots floater-history; do
    [[ -d "$BIN_DIR/$d" ]] && mv "$BIN_DIR/$d" "$ARCHIVE_ROOT/$d"
done

echo "## DONE. Archived under $ARCHIVE_ROOT"
exit "$RUN_RC"
