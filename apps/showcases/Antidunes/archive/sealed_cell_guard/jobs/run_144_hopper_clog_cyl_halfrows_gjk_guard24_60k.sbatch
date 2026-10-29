#!/bin/bash -l
#SBATCH --job-name=antidunes-hopper-clog-cyl-halfrows-gjk-guard24-60k
#SBATCH --partition=multinode
#SBATCH --nodes=2
#SBATCH --ntasks=144
#SBATCH --time=02:30:00
#SBATCH --output=sbatch-output/output/hopper_clog_cyl_halfrows_gjk_guard24_60k_%j.out
#SBATCH --error=sbatch-output/error/hopper_clog_cyl_halfrows_gjk_guard24_60k_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# [Claude, 2026-09-23] GUARDED RERUN of job 4255679 (antidunes-hopper-clog-cyl-halfrows-gjk-noguards-60k), which
# aborted at step 3133 via STALE-PARTICLE-IDX. Root cause confirmed from that run's own WARNING-level log: a liquid
# fluid cell sealed between rods uid ...885/...888 (real centre separation 5.93 LU, essentially touching) developed
# a garbage population mode ([MEM-GARBAGE-LINK], step 3132) and by step 3133 both rods were launched to ~138,000
# LU/step (vs ~0.02 LU/step normal), triggering the abort -- the same "class B" sealed/pinched-cell instability
# already root-caused (job 4252462) and already fixed (job 4255664, pillared full-density hopper: guard24 took that
# scenario from a ~step-29000 crash to a full clean 40000-step survival) via SealedCellGuard.h. This run uses the
# IDENTICAL 4255679 config (210-rod halfrows population, GJK-only, no rodcap/fluid-sanity-sweep) with
# closeSealedCellsMinNeighbours=24 added -- see prms/hopper_clog_cyl_frontblock_pidon_2node_halfrows_guard24_gas.prm
# for the full writeup. Submitted standalone (no --dependency): this checkout (~/mt-dennys-mirror) was confirmed
# free of running/pending jobs before submission (the other two checkouts, ~/mt-dennys and ~/mt-dennys-mirror2, each
# have their own pending clog-measurement/bogner jobs queued -- this run does not touch either). Checkpointing is ON
# (storeSnapshot true, snapshotFrequency 1000, finer than the usual 5000 so a checkpoint sits close to the original
# crash region for cheap inspection/resume either way this goes).
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

# [Claude, 2026-09-18] Never DELETE a previous run's output -- move it aside instead (same filesystem -> instant
# rename); archive/rescued_* is then moved off by hand. (Confirmed clean before this submission; kept for safety.)
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
echo "## clog measurement: rod hopper WITH pillars, PID on, sealed-cell guard 24, GJK-only, 60000 steps, 2-node/144-rank"
echo "############################################################"
# a crashed run must still reach the archive step below (set -e would skip it), so capture the exit status
set +e
"$SH" --run --np 144 --prm "$PRM_DIR/hopper_clog_cyl_frontblock_pidon_2node_halfrows_guard24_gas.prm"
RUN_RC=$?
set -e
echo "## run exited with status $RUN_RC -- archiving whatever output exists"

# Prefer $HPCVAULT, decided and verified HERE (right before use).
STAMP="$(date +%Y%m%d_%H%M%S)"
VAULT_CANDIDATE="${HPCVAULT:-}/antidunes-archive/${SLURM_JOB_NAME}_${SLURM_JOB_ID}_${STAMP}"
if [[ -n "${HPCVAULT:-}" ]] && mkdir -p "$VAULT_CANDIDATE" 2>/dev/null \
   && touch "$VAULT_CANDIDATE/.write_test" 2>/dev/null; then
    rm -f "$VAULT_CANDIDATE/.write_test"
    ARCHIVE_ROOT="$VAULT_CANDIDATE"
else
    ARCHIVE_ROOT="$BIN_DIR/archive/${STAMP}_hopper_clog_cyl_halfrows_gjk_guard24_60k"
    mkdir -p "$ARCHIVE_ROOT"
    echo "WARNING: \$HPCVAULT unavailable at archive time -- falling back to LOCAL archive at $ARCHIVE_ROOT (move this off manually to avoid quota pressure)" >&2
fi
for d in vtk-out eval snapshots; do
    [[ -d "$BIN_DIR/$d" ]] && mv "$BIN_DIR/$d" "$ARCHIVE_ROOT/$d"
done

echo "## DONE. Archived under $ARCHIVE_ROOT"
exit "$RUN_RC"
