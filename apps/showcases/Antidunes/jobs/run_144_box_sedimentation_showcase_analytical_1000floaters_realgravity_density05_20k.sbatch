#!/bin/bash -l
#SBATCH --job-name=antidunes-box-sedimentation-showcase-analytical-1000floaters-realgravity-density05-20k
#SBATCH --partition=multinode
#SBATCH --constraint=el9
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=72
#SBATCH --time=03:00:00
#SBATCH --output=sbatch-output/output/antidunes_box_sedimentation_showcase_analytical_1000floaters_realgravity_density05_20k_%j.out
#SBATCH --error=sbatch-output/error/antidunes_box_sedimentation_showcase_analytical_1000floaters_realgravity_density05_20k_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# Uniform-density (all 0.5, all floating) comparison sibling of
# run_144_box_sedimentation_showcase_analytical_1000floaters_realgravity.sbatch, which mixes 5
# density classes and reliably crashes at steps 3315-4204 (jobs 4106364, 4112811, 4114357, 4115699,
# 4116319). This run has no logFloaterHistory and no dense-VTK window -- classic VTK output only,
# to see how far a uniform, non-sinking floater population gets within 20000 steps. See the prm's
# own header for the full rationale.
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
    --prm "$REPO/apps/showcases/Antidunes/prms/box_sedimentation_showcase_analytical_1000floaters_realgravity_density05_20k_gas.prm"

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
