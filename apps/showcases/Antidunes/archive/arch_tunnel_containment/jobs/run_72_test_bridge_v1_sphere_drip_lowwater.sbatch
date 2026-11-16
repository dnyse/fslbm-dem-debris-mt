#!/bin/bash -l
#SBATCH --job-name=antidunes-test-bridge-v1-sphere-drip-lowwater
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=06:00:00
#SBATCH --output=sbatch-output/output/antidunes_test_bridge_v1_sphere_drip_lowwater_%j.out
#SBATCH --error=sbatch-output/error/antidunes_test_bridge_v1_sphere_drip_lowwater_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# Long-duration bridge floater-contact test: 130000 steps, low-water free surface
# (liquidHeightFactor 2.0), 100 radius-3 densityRatio=0.5 spheres released one every 1000 steps
# starting at step 25000 (after the flow field develops), cycling through a 10-position Y-row via
# startInactive. useBognerThreePhaseMEM/useMemOnlyBuoyancy both explicitly true (the buoyancy-flag
# fix -- see prms/test_bridge_v1_sphere_drip_lowwater_gas.prm's header for why this matters: without
# it, every floater in the earlier mixed_75k run sank to the domain floor regardless of density).
# Time budget generous (6h vs the mixed_75k run's 56min for 75k steps/8 floaters): up to 100
# concurrent floaters plus a longer run raises DEM contact cost meaningfully.
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
    --prm "$REPO/apps/showcases/Antidunes/prms/test_bridge_v1_sphere_drip_lowwater_gas.prm"

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
