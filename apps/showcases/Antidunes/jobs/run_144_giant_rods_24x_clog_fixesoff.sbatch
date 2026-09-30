#!/bin/bash -l
#SBATCH --job-name=antidunes-giant-rods-24x-clog-fixesoff
#SBATCH --partition=multinode
#SBATCH --nodes=2
#SBATCH --ntasks=144
#SBATCH --time=06:00:00
#SBATCH --output=sbatch-output/output/giant_rods_24x_clog_fixesoff_%j.out
#SBATCH --error=sbatch-output/error/giant_rods_24x_clog_fixesoff_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# [Claude, 2026-09-26] Twenty-four giant rods (8x1x3, second doubling step in this series, after
# the 12-rod scenario 4268493), per explicit request ("keep doubling them in x until into our old
# issues again"). See prms/giant_rods_24x_clog_gas.prm's own header for the full derivation: same
# rod geometry (radius 6, length 36, 1:3 ratio) and 3 Z-layer densities (0.85/0.75/0.65) as every
# scenario in this series, but 8 rows in X (364/350/336/322/308/294/280/266, 14 LU pitch) instead of
# 4. Run on THIS mirror checkout with all three of 2026-09-25's crash-mitigation fixes (Option A/B,
# SimpleBB.h guard) reverted -- a fixes-off baseline, so if crowding this many rods into the channel
# reproduces one of this session's earlier-diagnosed corruption mechanisms, it's a genuine,
# undefended reproduction. 51000 timesteps, 6h cap.
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
echo "## twenty-four giant rods (8x1x3), fixes OFF (Option A/B, SimpleBB.h guard all reverted)"
echo "############################################################"
set +e
"$SH" --run --np 144 --prm "$PRM_DIR/giant_rods_24x_clog_gas.prm"
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
    ARCHIVE_ROOT="$BIN_DIR/archive/${STAMP}_giant_rods_24x_clog_fixesoff"
    mkdir -p "$ARCHIVE_ROOT"
    echo "WARNING: \$HPCVAULT unavailable at archive time -- falling back to LOCAL archive at $ARCHIVE_ROOT (move this off manually to avoid quota pressure)" >&2
fi
for d in vtk-out eval snapshots; do
    [[ -d "$BIN_DIR/$d" ]] && mv "$BIN_DIR/$d" "$ARCHIVE_ROOT/$d"
done

echo "## DONE. Archived under $ARCHIVE_ROOT"
exit "$RUN_RC"
