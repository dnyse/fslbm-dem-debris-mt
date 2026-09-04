#!/bin/bash -l
#SBATCH --job-name=antidunes-cylinder-sedimentation-showcase-analytical
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=03:00:00
#SBATCH --output=sbatch-output/output/antidunes_cylinder_sedimentation_showcase_analytical_%j.out
#SBATCH --error=sbatch-output/error/antidunes_cylinder_sedimentation_showcase_analytical_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# Analytical-contact-detection smoke test: identical to run_72_cylinder_sedimentation_showcase.sbatch
# (same dense 100-log 5x4x5 grid, 12 LU spacing, "end-on and side-on collisions ... expected")
# except the prm sets useAnalyticalContactDetection true, routing every Cylinder x Cylinder
# contact through PaperAnalyticalContacts.h's HYBRID dispatcher (Kodam et al. 2010 + Guo et al.
# 2012 analytical result, gated by cylinderCylinderResultIsUncertain() which falls back to
# GJK/EPA whenever the analytical magnitude can't be trusted) instead of pure GJK/EPA. This is
# the densest, most collision-heavy scenario in the project -- a real stress test of the fallback
# gate under load. See prms/cylinder_sedimentation_showcase_analytical_gas.prm's header.
# domainSize <128, 60, 180> / cellsPerBlock <32, 20, 30> = 4x3x6 = 72 blocks
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
    --prm "$REPO/apps/showcases/Antidunes/prms/cylinder_sedimentation_showcase_analytical_gas.prm"

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
