#!/bin/bash -l
#SBATCH --job-name=antidunes-bridge-lv3-DEBUG-crashrepro-floaterlog
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=01:00:00
#SBATCH --output=sbatch-output/output/antidunes_test_bridge_v1_raisedwater_bigjam_mixed_level3_DEBUG_crashrepro_floaterlog_%j.out
#SBATCH --error=sbatch-output/error/antidunes_test_bridge_v1_raisedwater_bigjam_mixed_level3_DEBUG_crashrepro_floaterlog_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# Reproduces job 4166006's exact confirmed crash (Z=74 spawn height, ORIGINAL/untripled release
# timing) with ModelParameters.logFloaterHistory turned on (unfiltered -- safe here since the crash
# happens very early, ~step 8001, so data volume stays small). Goal: fine-grained per-sub-cycle
# force/torque/velocity breakdown for uid=576460752303423513 (the box that fell at only 65% of pure
# free-fall speed with an unexplained Y-drift -- lubrication ruled out, box-only pair) in the steps
# leading up to its collision with uid=1441151880758558745.
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
    --prm "$REPO/apps/showcases/Antidunes/prms/test_bridge_v1_raisedwater_bigjam_mixed_level3_DEBUG_crashrepro_floaterlog_gas.prm"

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
