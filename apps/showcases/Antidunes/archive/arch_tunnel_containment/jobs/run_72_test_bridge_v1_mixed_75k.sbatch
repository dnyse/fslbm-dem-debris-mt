#!/bin/bash -l
#SBATCH --job-name=antidunes-test-bridge-v1-mixed-75k
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=04:00:00
#SBATCH --output=sbatch-output/output/antidunes_test_bridge_v1_mixed_75k_%j.out
#SBATCH --error=sbatch-output/error/antidunes_test_bridge_v1_mixed_75k_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# Extended bridge floater-contact test: 75000 steps (vs. the 30k baseline in
# run_72_test_bridge_v1.sbatch), full mix of all three floater shapes
# (sphere/box/cylinder) sized to match mixed_shape_sedimentation_showcase_gas.prm
# (radius/halfExtents 4.5/<3.5,7,3.5>/4.5, up from the original radius-2 spheres),
# all densityRatio 0.5. Same piers+arches geometry, same radial-gate fix in
# ComputeTunnelContacts (verified clean on job 4065061 at the smaller scale) --
# this run checks the fix holds at floater sizes close to the arch radii
# themselves (sphere r=4.5 vs left/right arch r=5) and over a much longer window.
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
    --prm "$REPO/apps/showcases/Antidunes/prms/test_bridge_v1_mixed_75k_gas.prm"

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
