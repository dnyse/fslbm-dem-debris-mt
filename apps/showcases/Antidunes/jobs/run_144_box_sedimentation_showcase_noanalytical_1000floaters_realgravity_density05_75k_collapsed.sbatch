#!/bin/bash -l
#SBATCH --job-name=antidunes-box-sedimentation-showcase-noanalytical-1000floaters-realgravity-density05-75k-collapsed
#SBATCH --partition=multinode
#SBATCH --constraint=el9
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=72
#SBATCH --time=08:00:00
#SBATCH --output=sbatch-output/output/antidunes_box_sedimentation_showcase_noanalytical_1000floaters_realgravity_density05_75k_collapsed_%j.out
#SBATCH --error=sbatch-output/error/antidunes_box_sedimentation_showcase_noanalytical_1000floaters_realgravity_density05_75k_collapsed_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# useAnalyticalContactDetection=false sibling of
# run_144_box_sedimentation_showcase_analytical_1000floaters_realgravity_density05_75k_collapsed_v3.sbatch
# (job 4125031, analytical ON, crashed at step 18188 with the same isolated-cell [LOW-RHO]
# mechanism, plus an unexplained long tail of floaters still well above the waterline at step
# 18000). Identical setup otherwise (spawn Z=118, 1000-step release interval, densityRatio 0.5,
# timesteps 75000, VTK writeFrequency 2000 velocity-only, logFloaterHistory off) -- pure GJK/EPA
# fallback for every Box x Box pair instead of the closed-form analytical path, to see whether that
# changes the crash timing or the high-floater tail. See the prm's own header.
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
    --prm "$REPO/apps/showcases/Antidunes/prms/box_sedimentation_showcase_noanalytical_1000floaters_realgravity_density05_75k_collapsed_gas.prm"

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
