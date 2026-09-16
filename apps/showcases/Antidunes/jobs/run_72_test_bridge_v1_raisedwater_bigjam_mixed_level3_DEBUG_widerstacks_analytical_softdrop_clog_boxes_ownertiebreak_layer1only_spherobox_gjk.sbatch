#!/bin/bash -l
#SBATCH --job-name=antidunes-bridge-lv3-DEBUG-widerstacks-analytical-softdrop-clog-boxes-ownertiebreak-layer1only-spherobox_gjk
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=06:00:00
#SBATCH --output=sbatch-output/output/antidunes_test_bridge_v1_raisedwater_bigjam_mixed_level3_DEBUG_widerstacks_analytical_softdrop_clog_boxes_ownertiebreak_layer1only-spherobox_gjk_%j.out
#SBATCH --error=sbatch-output/error/antidunes_test_bridge_v1_raisedwater_bigjam_mixed_level3_DEBUG_widerstacks_analytical_softdrop_clog_boxes_ownertiebreak_layer1only-spherobox_gjk_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# Rerun with useSpheroBoxContact=true, useSpheroBoxAnalytical=false (GJK/EPA-with-margin path, NOT the
# Eberly closed-form alternative -- that gets its own separate rerun). useMultiPointBoxContact and
# useBoxOverlapCap both TURNED OFF (superseded). See the spherobox_gjk prm's own header for full context:
# jobs 4201656 (useMultiPointBoxContact alone) and 4202194 (+ useBoxOverlapCap) both FAILED at
# RELEASE-OVERLAP step 122001 / STALE-PARTICLE-IDX step 106903 -- both downstream consequences of the
# same already-diagnosed flying-floater instability. useSpheroBoxContact is this session's
# best-validated Box x Box mitigation (9/9 PASS at both e=1.00 and e=0.30 on the local DEMValidation
# edge-case battery, best/tied-best on the real job-4199585 incident replay) but this is its FIRST run
# at actual production scale. logFloaterHistory is ON for the ENTIRE run (window [0,250000]) since a
# recurrence's step, if any, is not yet known the way prior incidents' were. Run from a SEPARATE cluster
# checkout so it can't race any other run's build/output.
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
    --prm "$REPO/apps/showcases/Antidunes/prms/test_bridge_v1_raisedwater_bigjam_mixed_level3_DEBUG_widerstacks_analytical_softdrop_clog_boxes_ownertiebreak_layer1only_spherobox_gjk_gas.prm"

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
