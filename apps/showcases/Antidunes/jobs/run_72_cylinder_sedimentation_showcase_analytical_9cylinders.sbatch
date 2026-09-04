#!/bin/bash -l
#SBATCH --job-name=antidunes-cylinder-sedimentation-showcase-analytical-9cylinders
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=03:00:00
#SBATCH --output=sbatch-output/output/antidunes_cylinder_sedimentation_showcase_analytical_9cylinders_%j.out
#SBATCH --error=sbatch-output/error/antidunes_cylinder_sedimentation_showcase_analytical_9cylinders_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# Cylinder conversion of run_72_box_sedimentation_showcase_analytical.sbatch: same v6 9-shape
# Z-stacking-only layout (single Y position per column, no Y-neighbour contacts) that never
# crashed for boxes, with the 9 FloatingBoxes replaced by matching-footprint FloatingCylinders
# (radius 3.5, length 14, axis Y). useAnalyticalContactDetection true routes every
# Cylinder x Cylinder contact through PaperAnalyticalContacts.h's cylinderCylinderContact hybrid
# path (analytical + GJK/EPA fallback gate) instead of pure GJK/EPA. NOT the same scenario as
# run_72_cylinder_sedimentation_showcase_analytical.sbatch (a much denser, unrelated 100-log
# grid). See prms/cylinder_sedimentation_showcase_analytical_9cylinders_gas.prm's header for the
# full rationale.
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
    --prm "$REPO/apps/showcases/Antidunes/prms/cylinder_sedimentation_showcase_analytical_9cylinders_gas.prm"

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
