#!/bin/bash -l
#SBATCH --job-name=antidunes-box-sedimentation-showcase-analytical-1000floaters-realgravity
#SBATCH --partition=multinode
#SBATCH --constraint=el9
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=72
#SBATCH --time=10:00:00
#SBATCH --output=sbatch-output/output/antidunes_box_sedimentation_showcase_analytical_1000floaters_realgravity_%j.out
#SBATCH --error=sbatch-output/error/antidunes_box_sedimentation_showcase_analytical_1000floaters_realgravity_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# Genuine 10x10x10 grid (1000 floaters) sibling of
# run_144_box_sedimentation_showcase_analytical_125floaters_realgravity.sbatch: same box size
# (halfExtents <3,9,3>), same real Y-neighbour-contact design, same true g_SI=9.81 gravity. Domain
# scaled proportionally in X/Y (160,160 -> 320,320), Z unchanged (180); cellsPerBlock scaled
# proportionally too (40,40,20 -> 80,80,20), so process count stays IDENTICAL at 144 (4x4x9) -- no
# new node/core budget versus the 125-floater run. Time limit raised (4h -> 10h) since both the
# fluid domain (~4x cells, X and Y doubled) and the particle/DEM workload (8x floaters) scale up
# substantially. See prms/box_sedimentation_showcase_analytical_1000floaters_realgravity_gas.prm's
# header for the full derivation.
if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    REPO="$(cd "$SLURM_SUBMIT_DIR/../../.." && pwd)"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    REPO="$(cd "$SCRIPT_DIR/../../../.." && pwd)"
fi
BIN_DIR="$REPO/build/apps/showcases/Antidunes"

"$REPO/apps/showcases/Antidunes/antidunes.sh" \
    --build --jobs "$(nproc)" \
    --run --np 144 \
    --prm "$REPO/apps/showcases/Antidunes/prms/box_sedimentation_showcase_analytical_1000floaters_realgravity_gas.prm"

STAMP="$(date +%Y%m%d_%H%M%S)"
if [[ -n "$HPCVAULT" ]]; then
    ARCHIVE_DIR="$HPCVAULT/antidunes-archive/${SLURM_JOB_NAME}_$STAMP"
else
    ARCHIVE_DIR="$BIN_DIR/archive/$STAMP"
fi
mkdir -p "$ARCHIVE_DIR"
# `if` (not `[[ ]] &&`) so a missing directory's false test doesn't leak into this being the
# script's own final exit status -- see the 125-floater job's own comment (job 4104656) for the
# full mechanism this avoids.
for dir in vtk-out eval snapshots floater-history; do
    if [[ -d "$BIN_DIR/$dir" ]]; then mv "$BIN_DIR/$dir" "$ARCHIVE_DIR/$dir"; fi
done
