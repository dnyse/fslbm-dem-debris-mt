#!/bin/bash -l
#SBATCH --job-name=antidunes-boxes-24x-2x4x3-clog
#SBATCH --partition=multinode
#SBATCH --nodes=2
#SBATCH --ntasks=144
#SBATCH --time=06:00:00
#SBATCH --output=sbatch-output/output/boxes_24x_2x4x3_clog_%j.out
#SBATCH --error=sbatch-output/error/boxes_24x_2x4x3_clog_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# [Claude, 2026-09-30] Shape-comparison companion to run_144_small_rods_24x_2x4x3_clog_fixesoff.sbatch,
# per explicit request ("run this one of the free directory with cylinders and on the other free one
# with cubes"): identical scenario (domain, pillars, PID, physics, ModelParameters, all 24
# positions/densities) with the 24 FloatingCylinders replaced by 24 FloatingBoxes of the same
# footprint (halfExtents <3,9,3>, matching the rod's radius-3/length-18 bounding box) -- see
# prms/boxes_24x_2x4x3_clog_gas.prm's own header. Tests whether the STALE-PARTICLE-IDX cross-rank
# ghost-ownership desync found in the rod doubling series (THESIS_PLAN.md) is specific to the
# cylinder's line-contact geometry or a general dense-body-pile effect. Option A/B and the (fully
# reverted) sealed-cell guard are off, matching the rod run exactly. 51000 timesteps, 6h cap.
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
echo "## twenty-four boxes (2x4x3), shape-comparison companion to the small-rod run"
echo "############################################################"
set +e
"$SH" --run --np 144 --prm "$PRM_DIR/boxes_24x_2x4x3_clog_gas.prm"
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
    ARCHIVE_ROOT="$BIN_DIR/archive/${STAMP}_boxes_24x_2x4x3_clog"
    mkdir -p "$ARCHIVE_ROOT"
    echo "WARNING: \$HPCVAULT unavailable at archive time -- falling back to LOCAL archive at $ARCHIVE_ROOT (move this off manually to avoid quota pressure)" >&2
fi
for d in vtk-out eval snapshots; do
    [[ -d "$BIN_DIR/$d" ]] && mv "$BIN_DIR/$d" "$ARCHIVE_ROOT/$d"
done

echo "## DONE. Archived under $ARCHIVE_ROOT"
exit "$RUN_RC"
