#!/bin/bash -l
#SBATCH --job-name=antidunes-bridge-lv3-DEBUG-release34000
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=02:30:00
#SBATCH --output=sbatch-output/output/antidunes_test_bridge_v1_raisedwater_bigjam_mixed_level3_DEBUG_release34000_%j.out
#SBATCH --error=sbatch-output/error/antidunes_test_bridge_v1_raisedwater_bigjam_mixed_level3_DEBUG_release34000_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# FOURTH debug attempt for the "flying floater" chain traced from job 4164460. assertNoOverlapAtRelease()
# (SetupStages.h) hard-aborts, with both floaters' uid/shapeType/pos/interactionRadius/overlap in the
# message, if a floater is released while overlapping an already-active neighbour. Job 4165898 (tripled
# release timing) aborted, but on what looks like a bounding-sphere false positive for an elongated box
# pair -- see the prm's own header for the full geometric check. THIS run reverts to the ORIGINAL
# (untripled) release timing (3000/9000-step gaps) already shown insufficient, at timesteps=100000, to
# get a clean baseline: expect the assertion to fire early, and each firing needs the same manual
# axis-extent sanity check before being treated as a confirmed real overlap.
# See prms/test_bridge_v1_raisedwater_bigjam_mixed_level3_DEBUG_release34000_gas.prm's own header.
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
    --prm "$REPO/apps/showcases/Antidunes/prms/test_bridge_v1_raisedwater_bigjam_mixed_level3_DEBUG_release34000_gas.prm"

STAMP="$(date +%Y%m%d_%H%M%S)"
if [[ -n "$HPCVAULT" ]]; then
    ARCHIVE_DIR="$HPCVAULT/antidunes-archive/${SLURM_JOB_NAME}_$STAMP"
else
    ARCHIVE_DIR="$BIN_DIR/archive/$STAMP"
fi
mkdir -p "$ARCHIVE_DIR" || ARCHIVE_DIR="$BIN_DIR/archive/$STAMP"
[[ -d "$ARCHIVE_DIR" ]] || mkdir -p "$ARCHIVE_DIR"
# `if` form (not `[[ -d ]] && mv`) so a legitimately-absent floater-history/ (logFloaterHistory is
# off) doesn't trip `set -e` into a false "FAILED" -- see run_72_test_bridge_v1_lowwater_bigjam_mixed.sbatch's own comment for the full story.
for dir in vtk-out eval snapshots floater-history; do
    if [[ -d "$BIN_DIR/$dir" ]]; then
        mv "$BIN_DIR/$dir" "$ARCHIVE_DIR/$dir"
    fi
done
