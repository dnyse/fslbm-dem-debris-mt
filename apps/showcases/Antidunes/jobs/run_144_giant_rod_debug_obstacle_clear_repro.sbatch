#!/bin/bash -l
#SBATCH --job-name=antidunes-giant-rod-debug-obstacle-clear-repro
#SBATCH --partition=multinode
#SBATCH --nodes=2
#SBATCH --ntasks=144
#SBATCH --time=06:00:00
#SBATCH --output=sbatch-output/output/giant_rod_debug_obstacle_clear_repro_%j.out
#SBATCH --error=sbatch-output/error/giant_rod_debug_obstacle_clear_repro_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# [Claude, 2026-09-26] DETERMINISTIC REPRO of the small-rod job's (4268750) NaN mechanism, per
# explicit request ("extract these situations and replicate it in a smaller simulation"). Adapted
# from the already-validated giant_rod_clog_gas.prm (one rod, jams at the pillar early, 144 ranks) --
# see prms/giant_rod_debug_obstacle_clear_repro_gas.prm's own header for the full derivation. At step
# 3000 (well after the rod has settled into its jammed position), ParticleKernels.h's
# [DEBUG-FORCE-OBSTACLE-CLEAR] hook deliberately clears this rod's own current obstacle cells --
# reproducing the reverted STALE-PARTICLE-IDX guard's side effect without needing the rare cross-rank
# desync that originally triggered it. Fixes-off, matching every other scenario in this series.
# 6000 timesteps (short -- the point is the step-3000 trigger, not a long run), 6h cap (generous).
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
echo "## giant single-rod debug-obstacle-clear repro, fixes OFF (Option A/B, SimpleBB.h guard all reverted)"
echo "############################################################"
set +e
"$SH" --run --np 144 --prm "$PRM_DIR/giant_rod_debug_obstacle_clear_repro_gas.prm"
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
    ARCHIVE_ROOT="$BIN_DIR/archive/${STAMP}_giant_rod_debug_obstacle_clear_repro"
    mkdir -p "$ARCHIVE_ROOT"
    echo "WARNING: \$HPCVAULT unavailable at archive time -- falling back to LOCAL archive at $ARCHIVE_ROOT (move this off manually to avoid quota pressure)" >&2
fi
for d in vtk-out eval snapshots; do
    [[ -d "$BIN_DIR/$d" ]] && mv "$BIN_DIR/$d" "$ARCHIVE_ROOT/$d"
done

echo "## DONE. Archived under $ARCHIVE_ROOT"
exit "$RUN_RC"
