#!/bin/bash -l
#SBATCH --job-name=antidunes-bridge-lv3-DEBUG-widerstacks-nogjk-lowwater-clog
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=05:00:00
#SBATCH --output=sbatch-output/output/antidunes_test_bridge_v1_raisedwater_bigjam_mixed_level3_DEBUG_widerstacks_nogjk_lowwater_clog_%j.out
#SBATCH --error=sbatch-output/error/antidunes_test_bridge_v1_raisedwater_bigjam_mixed_level3_DEBUG_widerstacks_nogjk_lowwater_clog_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# GUARANTEED-CLOGGING run: boxes-only (no spheres/cylinders), 270 floating boxes across 18 spawn
# points (both upstream rows x all 9 Y aim-positions, covering every arch's own width, not just its
# centre), 15 waves per point at the proven-safe 9000-step gap, liquidHeightFactor=2.0 (lowwater,
# job 4167764's setting), 150000 steps so the jam has time to develop and hold after the last
# release (step 131000). See the prm's own header for the full derivation from the lowwater file.
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
    --prm "$REPO/apps/showcases/Antidunes/prms/test_bridge_v1_raisedwater_bigjam_mixed_level3_DEBUG_widerstacks_nogjk_lowwater_clog_gas.prm"

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
