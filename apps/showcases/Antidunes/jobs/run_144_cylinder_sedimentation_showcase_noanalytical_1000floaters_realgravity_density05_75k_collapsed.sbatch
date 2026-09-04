#!/bin/bash -l
#SBATCH --job-name=antidunes-cylinder-sedimentation-showcase-noanalytical-1000floaters-realgravity-density05-75k-collapsed
#SBATCH --partition=multinode
#SBATCH --constraint=el9
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=72
#SBATCH --time=08:00:00
#SBATCH --output=sbatch-output/output/antidunes_cylinder_sedimentation_showcase_noanalytical_1000floaters_realgravity_density05_75k_collapsed_%j.out
#SBATCH --error=sbatch-output/error/antidunes_cylinder_sedimentation_showcase_noanalytical_1000floaters_realgravity_density05_75k_collapsed_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# NOT YET SUBMITTED as of writing -- prepared to run only after confirming
# run_144_box_sedimentation_showcase_analytical_1000floaters_realgravity_density05_75k_collapsed_v3.sbatch
# (job 4125031) completes/behaves as expected.
#
# Cylinder sibling of that box scenario: same collapsed-plane, staggered-drop setup (spawn Z=118,
# 1000-step release interval, all released by step 9000, densityRatio 0.5, timesteps 75000, VTK
# writeFrequency 2000 velocity-only, logFloaterHistory off) but with FloatingCylinder (radius 3,
# length 18) instead of FloatingBox, and useAnalyticalContactDetection false (pure GJK/EPA for
# every Cylinder x Cylinder pair, matching the existing *_noanalytical_* convention -- typically
# slower per-step than the closed-form analytical path, hence the longer --time budget here).
# See the prm's own header for the full rationale.
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
    --prm "$REPO/apps/showcases/Antidunes/prms/cylinder_sedimentation_showcase_noanalytical_1000floaters_realgravity_density05_75k_collapsed_gas.prm"

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
