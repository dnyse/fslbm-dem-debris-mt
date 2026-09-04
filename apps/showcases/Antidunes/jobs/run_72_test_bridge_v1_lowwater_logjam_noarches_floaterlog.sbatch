#!/bin/bash -l
#SBATCH --job-name=antidunes-test-bridge-v1-lowwater-logjam-noarches-floaterlog
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=06:00:00
#SBATCH --output=sbatch-output/output/antidunes_test_bridge_v1_lowwater_logjam_noarches_floaterlog_%j.out
#SBATCH --error=sbatch-output/error/antidunes_test_bridge_v1_lowwater_logjam_noarches_floaterlog_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# Geometry-control mirror of job 4134552 (test_bridge_v1_lowwater_logjam_floaterlog, FloatingBoxes):
# same 9 log-shaped FloatingBoxes at the same positions/settlingSteps, same piers, but the arches
# (DEM ArchTunnels + LBM arch BodyBC, all 3 of them) are removed entirely -- only the two piers
# remain as solid geometry, deck/arch volume is fully open. Mirror of the sphere-shape control
# (run_72_test_bridge_v1_lowwater_logjam_spheres_floaterlog.sbatch): that one held geometry fixed
# and swapped floater shape; this one holds floater shape fixed and removes the geometry that
# produced the snag, to separate "box shape artifact" from "arch tunnel geometry artifact". See
# prms/test_bridge_v1_lowwater_logjam_noarches_floaterlog_gas.prm's own header for full rationale.
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
    --prm "$REPO/apps/showcases/Antidunes/prms/test_bridge_v1_lowwater_logjam_noarches_floaterlog_gas.prm"

STAMP="$(date +%Y%m%d_%H%M%S)"
if [[ -n "$HPCVAULT" ]]; then
    ARCHIVE_DIR="$HPCVAULT/antidunes-archive/${SLURM_JOB_NAME}_$STAMP"
else
    ARCHIVE_DIR="$BIN_DIR/archive/$STAMP"
fi
mkdir -p "$ARCHIVE_DIR"
for dir in vtk-out eval snapshots floater-history; do
    [[ -d "$BIN_DIR/$dir" ]] && mv "$BIN_DIR/$dir" "$ARCHIVE_DIR/$dir"
done
