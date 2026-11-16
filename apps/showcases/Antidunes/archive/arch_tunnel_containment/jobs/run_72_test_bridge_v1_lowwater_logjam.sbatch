#!/bin/bash -l
#SBATCH --job-name=antidunes-test-bridge-v1-lowwater-logjam
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=06:00:00
#SBATCH --output=sbatch-output/output/antidunes_test_bridge_v1_lowwater_logjam_%j.out
#SBATCH --error=sbatch-output/error/antidunes_test_bridge_v1_lowwater_logjam_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# Log-jam test on the toy bridge (recessed-arch layout, restitutionCoefficient=0.5, low-water free
# surface, targetMeanVelocityMagnitude=0.01): 30000 steps, no floaters at all for the first 5000
# (flow field develops under the new PID target first), then 9 log-shaped FloatingBoxes
# (halfExtents<6,1,1.5>, densityRatio=0.6) release in a staggered X-Y grid (2 upstream rows,
# settlingSteps 5000-15000) and drift downstream to (hopefully) jam at the bridge. See
# prms/test_bridge_v1_lowwater_logjam_gas.prm's own header for the full layout/timing rationale.
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
    --prm "$REPO/apps/showcases/Antidunes/prms/test_bridge_v1_lowwater_logjam_gas.prm"

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
