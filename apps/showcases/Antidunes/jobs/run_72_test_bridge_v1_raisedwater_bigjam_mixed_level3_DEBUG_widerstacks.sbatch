#!/bin/bash -l
#SBATCH --job-name=antidunes-bridge-lv3-DEBUG-widerstacks
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=02:30:00
#SBATCH --output=sbatch-output/output/antidunes_test_bridge_v1_raisedwater_bigjam_mixed_level3_DEBUG_widerstacks_%j.out
#SBATCH --error=sbatch-output/error/antidunes_test_bridge_v1_raisedwater_bigjam_mixed_level3_DEBUG_widerstacks_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# PARALLEL run (separate mirror checkout) alongside the Z-height debug run on mt-dennys-mirror2.
# Same confirmed root cause (a floater released while overlapping an ACTIVE neighbour -- now
# detected exactly, via assertNoOverlapAtRelease/ComputeContacts's freshlyActivatedFloaterUIDs hook,
# SetupStages.h/ParticleKernels.h), but testing a DIFFERENT lever: job 4166006 showed a real overlap
# between two DIFFERENT, nearby aim points (row X=15 vs row X=22, only 7 LU apart), not just
# same-point sharing. THIS run widens the row-to-row separation (X=22 -> 40 for every row-2 aim
# point, all shapes) instead of raising spawn height. See the prm's own header for the full
# derivation.
set -euo pipefail

if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    REPO="$(cd "$SLURM_SUBMIT_DIR/../../.." && pwd)"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    REPO="$(cd "$SCRIPT_DIR/../../../.." && pwd)"
fi
BIN_DIR="$REPO/build/apps/showcases/Antidunes"

"$REPO/apps/showcases/Antidunes/antidunes.sh" \
    --build --jobs "$(nproc)" \
    --run --np 72 \
    --prm "$REPO/apps/showcases/Antidunes/prms/test_bridge_v1_raisedwater_bigjam_mixed_level3_DEBUG_widerstacks_gas.prm"

STAMP="$(date +%Y%m%d_%H%M%S)"
if [[ -n "$HPCVAULT" ]]; then
    ARCHIVE_DIR="$HPCVAULT/antidunes-archive/${SLURM_JOB_NAME}_$STAMP"
else
    ARCHIVE_DIR="$BIN_DIR/archive/$STAMP"
fi
mkdir -p "$ARCHIVE_DIR" || ARCHIVE_DIR="$BIN_DIR/archive/$STAMP"
[[ -d "$ARCHIVE_DIR" ]] || mkdir -p "$ARCHIVE_DIR"
for dir in vtk-out eval snapshots floater-history; do
    if [[ -d "$BIN_DIR/$dir" ]]; then
        mv "$BIN_DIR/$dir" "$ARCHIVE_DIR/$dir"
    fi
done
