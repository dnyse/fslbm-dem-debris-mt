#!/bin/bash -l
#SBATCH --job-name=antidunes-toy-trunk-first-staggered-rotated90
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=03:00:00
#SBATCH --output=sbatch-output/output/antidunes_test_bridge_v1_toy_trunk_first_staggered_rotated90_%j.out
#SBATCH --error=sbatch-output/error/antidunes_test_bridge_v1_toy_trunk_first_staggered_rotated90_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# 90-degree-rotated (broadside-to-flow) variant of the trunk-first-staggered toy scenario (job 4176138,
# COMPLETED cleanly at 150000 steps). See the prm's own header for the full rotation derivation and the
# geometric tension it introduces (a 22 LU broadside trunk cannot fit 3 non-overlapping same-time-release
# groups in this 60-wide domain with wall clearance -- mitigated via repositioning + staggered trunk
# release times, not fully eliminated). EXPECT this may hard-abort with RELEASE-OVERLAP; if so, that is
# anticipated, not a new bug -- see the prm's own header.
#
# ISOLATION: submit this job from THIS repo copy specifically (not the widerstacks-spheres job's copy --
# see that job's own comment) so the two jobs' build directories, running concurrently, cannot clobber
# each other's binaries/output.
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
    --prm "$REPO/apps/showcases/Antidunes/prms/test_bridge_v1_toy_trunk_first_staggered_rotated90_gas.prm"

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
