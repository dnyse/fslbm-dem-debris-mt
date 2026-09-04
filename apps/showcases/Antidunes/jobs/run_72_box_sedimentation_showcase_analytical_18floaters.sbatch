#!/bin/bash -l
#SBATCH --job-name=antidunes-box-sedimentation-showcase-analytical-18floaters
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=03:00:00
#SBATCH --output=sbatch-output/output/antidunes_box_sedimentation_showcase_analytical_18floaters_%j.out
#SBATCH --error=sbatch-output/error/antidunes_box_sedimentation_showcase_analytical_18floaters_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# Scaled-up variant of run_72_box_sedimentation_showcase_analytical.sbatch: same v6 layout
# methodology (single Y position per column, no Y-neighbour contacts, gravity-driven Z-stacking
# only) but with box size matched to thousand_cuboids_two_pillars_slitpass_gas.prm's
# halfExtents <3, 9, 3> instead of v6's own <3.5, 7, 3.5>, and the floater count doubled
# (9 -> 18) by doubling the Z-stack depth per column (3 rows -> 6 rows), not by adding
# Y-neighbour contacts. See prms/box_sedimentation_showcase_analytical_18floaters_gas.prm's
# header for the full rationale.
# domainSize <128, 60, 180> / cellsPerBlock <32, 20, 30> = 4x3x6 = 72 processes
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
    --prm "$REPO/apps/showcases/Antidunes/prms/box_sedimentation_showcase_analytical_18floaters_gas.prm"

STAMP="$(date +%Y%m%d_%H%M%S)"
if [[ -n "$HPCVAULT" ]]; then
    ARCHIVE_DIR="$HPCVAULT/antidunes-archive/${SLURM_JOB_NAME}_$STAMP"
else
    ARCHIVE_DIR="$BIN_DIR/archive/$STAMP"
fi
mkdir -p "$ARCHIVE_DIR"
for dir in vtk-out eval snapshots; do
    [[ -d "$BIN_DIR/$dir" ]] && mv "$BIN_DIR/$dir" "$ARCHIVE_DIR/$dir"
done
