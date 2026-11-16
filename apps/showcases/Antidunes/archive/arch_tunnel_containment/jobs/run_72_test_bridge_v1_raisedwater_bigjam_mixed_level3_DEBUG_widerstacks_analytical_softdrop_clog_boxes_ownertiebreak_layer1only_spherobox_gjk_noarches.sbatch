#!/bin/bash -l
#SBATCH --job-name=antidunes-bridge-lv3-DEBUG-widerstacks-analytical-softdrop-clog-boxes-ownertiebreak-layer1only-spherobox_gjk_noarches
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=06:00:00
#SBATCH --output=sbatch-output/output/antidunes_test_bridge_v1_raisedwater_bigjam_mixed_level3_DEBUG_widerstacks_analytical_softdrop_clog_boxes_ownertiebreak_layer1only-spherobox_gjk_noarches_%j.out
#SBATCH --error=sbatch-output/error/antidunes_test_bridge_v1_raisedwater_bigjam_mixed_level3_DEBUG_widerstacks_analytical_softdrop_clog_boxes_ownertiebreak_layer1only-spherobox_gjk_noarches_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# Same config as the spherobox_gjk rerun (useSpheroBoxContact=true, useSpheroBoxAnalytical=false),
# EXCEPT the arches are REMOVED ENTIRELY (both the DEM ArchTunnels block and the 3 arch BodyBC
# blocks -- not just a flag flip; see the noarches prm's own header for why
# useArchTunnelAsMesaPDObject=false alone is insufficient, that just swaps to a different DEM
# primitive, ComputeTunnelContacts, which still enforces containment). Only the piers remain as
# solid geometry. See the noarches prm's own header for full context: job 4212104 (spherobox_gjk,
# arches present) still FAILED at STALE-PARTICLE-IDX step 115258, and log analysis found every
# [VELOCITY-JUMP]/the pre-crash [MEM-WATCH] event clustering at z~60 (the channel ceiling) in the
# bridge/arch X-span (47-72) -- pointing at the arch tunnels themselves, not Box x Box contact, as
# the actual energy-injection mechanism. Run from a SEPARATE cluster checkout so it can't race any
# other run's build/output.
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
    --prm "$REPO/apps/showcases/Antidunes/prms/test_bridge_v1_raisedwater_bigjam_mixed_level3_DEBUG_widerstacks_analytical_softdrop_clog_boxes_ownertiebreak_layer1only_spherobox_gjk_noarches_gas.prm"

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
