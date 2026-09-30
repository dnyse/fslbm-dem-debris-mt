#!/bin/bash -l
#SBATCH --job-name=antidunes-small-rods-48x-2x8x3-clog-fixesoff
#SBATCH --partition=multinode
#SBATCH --nodes=2
#SBATCH --ntasks=144
#SBATCH --time=06:00:00
#SBATCH --output=sbatch-output/output/small_rods_48x_2x8x3_clog_fixesoff_%j.out
#SBATCH --error=sbatch-output/error/small_rods_48x_2x8x3_clog_fixesoff_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# [Claude, 2026-09-26] Forty-eight small rods (2x8x3, second doubling step in this series, after
# the 24-rod scenario 4268494), per explicit request ("keep doubling them in x until into our old
# issues again"). See prms/small_rods_48x_2x8x3_clog_gas.prm's own header for the full derivation:
# same rod geometry (radius 3, length 18, 1:3 ratio), same 2 Y-columns and 3 Z-layer densities
# (1.2/1.0/0.5) as every scenario in this series, but 8 rows in X
# (364/357/350/343/336/329/322/315, 7 LU pitch) instead of 4. Run on THIS mirror checkout with all
# three of 2026-09-25's crash-mitigation fixes (Option A/B, SimpleBB.h guard) reverted -- a
# fixes-off baseline, so if crowding this many rods into the channel reproduces one of this
# session's earlier-diagnosed corruption mechanisms, it's a genuine, undefended reproduction.
# 51000 timesteps, 6h cap.
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
echo "## forty-eight small rods (2x8x3), fixes OFF (Option A/B, SimpleBB.h guard all reverted)"
echo "############################################################"
set +e
"$SH" --run --np 144 --prm "$PRM_DIR/small_rods_48x_2x8x3_clog_gas.prm"
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
    ARCHIVE_ROOT="$BIN_DIR/archive/${STAMP}_small_rods_48x_2x8x3_clog_fixesoff"
    mkdir -p "$ARCHIVE_ROOT"
    echo "WARNING: \$HPCVAULT unavailable at archive time -- falling back to LOCAL archive at $ARCHIVE_ROOT (move this off manually to avoid quota pressure)" >&2
fi
for d in vtk-out eval snapshots; do
    [[ -d "$BIN_DIR/$d" ]] && mv "$BIN_DIR/$d" "$ARCHIVE_ROOT/$d"
done

echo "## DONE. Archived under $ARCHIVE_ROOT"
exit "$RUN_RC"
