#!/bin/bash -l
#SBATCH --job-name=antidunes-hopper-clog-cyl-halfrows-gjk-noguards-60k
#SBATCH --partition=multinode
#SBATCH --nodes=2
#SBATCH --ntasks=144
#SBATCH --time=02:30:00
#SBATCH --output=sbatch-output/output/hopper_clog_cyl_halfrows_gjk_noguards_60k_%j.out
#SBATCH --error=sbatch-output/error/hopper_clog_cyl_halfrows_gjk_noguards_60k_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# [Claude, 2026-09-22] NEW TEST (explicit instruction): the ORIGINAL cylinder hopper -- WITH pillars, GJK/EPA-only contacts, PID on, 60000 steps -- with NO guards/caps at all
# (no sealed-cell guard, no rodCapOverlap, no fluid sanity sweep, no dense NaN scan: the clean baseline, per prms/hopper_clog_cyl_frontblock_pidon_2node_halfrows_gas.prm's own
# absence of any of those keys) but HALF the rod population: 210 rods instead of 406 (`gen_hopper_clog_cyl_prm.py --rows 15`, front_x UNCHANGED at 368 so the kept rows are the
# DOWNSTREAM half against the pillars, X [267,371]; the dropped 14 rows are the FURTHER-UPSTREAM half, X [169,260] in the original 29-row block -- a whole contiguous half removed,
# not thinned throughout, so all 7 Z-layers / all three density bands (0.5/1.0/1.2) are still present in every remaining row). Question: does a half-density pile even jam/crash at
# all, and if it does, does it fail via the same mechanism (same rough step, same culprit-uid family) as the full-density pile despite half the rods? 144 ranks, 2h30m limit
# (no solid timing precedent for this exact half-density configuration; budgeted similar to the full-density sphere hopper's own 1h48-1h52m for 60000 steps, plus margin).
# Queue: afterany:4255664 (the currently RUNNING pillared-hopper sealed-cell-guard test), SAME checkout ~/mt-dennys-mirror.
# See prms/hopper_clog_cyl_frontblock_pidon_2node_halfrows_gas.prm and gen_hopper_clog_cyl_prm.py's own --rows option.
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
for dir in vtk-out eval snapshots; do
    if [[ -d "$BIN_DIR/$dir" ]]; then
        mkdir -p "$RESCUE"
        echo "WARNING: leftover $BIN_DIR/$dir from a previous run -- moving to $RESCUE (NOT deleting)" >&2
        mv "$BIN_DIR/$dir" "$RESCUE/$dir"
    fi
done

"$SH" --build --jobs "$(nproc)"

echo "############################################################"
echo "## clog measurement: rod hopper WITH pillars, PID on, sealed-cell guard 24, 40000 steps, 2-node/144-rank"
echo "############################################################"
# a crashed run must still reach the archive step below (set -e would skip it), so capture the exit status
set +e
"$SH" --run --np 144 --prm "$PRM_DIR/hopper_clog_cyl_frontblock_pidon_2node_halfrows_gas.prm"
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
    ARCHIVE_ROOT="$BIN_DIR/archive/${STAMP}_hopper_clog_cyl_halfrows_gjk_noguards_60k"
    mkdir -p "$ARCHIVE_ROOT"
    echo "WARNING: \$HPCVAULT unavailable at archive time -- falling back to LOCAL archive at $ARCHIVE_ROOT (move this off manually to avoid quota pressure)" >&2
fi
for d in vtk-out eval snapshots; do
    [[ -d "$BIN_DIR/$d" ]] && mv "$BIN_DIR/$d" "$ARCHIVE_ROOT/$d"
done

echo "## DONE. Archived under $ARCHIVE_ROOT"
exit "$RUN_RC"
