#!/bin/bash -l
#SBATCH --job-name=antidunes-bridge-lv3-DEBUG-widerstacks-analytical-softdrop-clog-boxes-ownertiebreak-layer1only-multipoint
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=06:00:00
#SBATCH --output=sbatch-output/output/antidunes_test_bridge_v1_raisedwater_bigjam_mixed_level3_DEBUG_widerstacks_analytical_softdrop_clog_boxes_ownertiebreak_layer1only-multipoint_%j.out
#SBATCH --error=sbatch-output/error/antidunes_test_bridge_v1_raisedwater_bigjam_mixed_level3_DEBUG_widerstacks_analytical_softdrop_clog_boxes_ownertiebreak_layer1only-multipoint_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# Rerun with useMultiPointBoxContact=true to test the SAT-normal-flip fix against this exact
# RELEASE-OVERLAP-crashing scenario (see the multipoint prm's own header for full context).
# BOX counterpart of the cylinder layer1only-multipoint rerun: useAnalyticalContactDetection=true (boxBoxSAT,
# already the base file's own setting), useTightLayer1=true (NEW Box x Box single-axis pre-filter),
# useContactFilterOwnerTiebreak=true (matching the cylinder sibling). NOTE: this exact scenario is
# already known to have crashed twice with RELEASE-OVERLAP (jobs 4175745/4179445) -- a pre-existing,
# unrelated floater-release-timing issue; may recur here regardless of today's work. Run from a
# SEPARATE cluster checkout so it can't race any other run's build or output directories.
set -euo pipefail

if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    REPO="$(cd "$SLURM_SUBMIT_DIR/../../.." && pwd)"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    REPO="$(cd "$SCRIPT_DIR/../../../.." && pwd)"
fi
BIN_DIR="$REPO/build/apps/showcases/Antidunes"

for dir in vtk-out eval snapshots floater-history; do
    if [[ -d "$BIN_DIR/$dir" ]]; then
        echo "WARNING: removing leftover $BIN_DIR/$dir from an unarchived previous run" >&2
        rm -rf "${BIN_DIR:?}/$dir"
    fi
done

"$REPO/apps/showcases/Antidunes/antidunes.sh" \
    --build --jobs "$(nproc)" \
    --run --np 72 \
    --prm "$REPO/apps/showcases/Antidunes/prms/test_bridge_v1_raisedwater_bigjam_mixed_level3_DEBUG_widerstacks_analytical_softdrop_clog_boxes_ownertiebreak_layer1only_multipoint_gas.prm"

STAMP="$(date +%Y%m%d_%H%M%S)"
if [[ -n "${HPCVAULT:-}" ]] && mkdir -p "$HPCVAULT/antidunes-archive/${SLURM_JOB_NAME}_${SLURM_JOB_ID}_$STAMP" 2>/dev/null; then
    ARCHIVE_DIR="$HPCVAULT/antidunes-archive/${SLURM_JOB_NAME}_$STAMP"
else
    ARCHIVE_DIR="$BIN_DIR/archive/$STAMP"
    mkdir -p "$ARCHIVE_DIR"
fi
for dir in vtk-out eval snapshots floater-history; do
    if [[ -d "$BIN_DIR/$dir" ]]; then
        mv "$BIN_DIR/$dir" "$ARCHIVE_DIR/$dir"
    fi
done
