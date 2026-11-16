#!/bin/bash -l
#SBATCH --job-name=antidunes-test-bridge-v1-overtop
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=04:00:00
#SBATCH --output=sbatch-output/output/antidunes_test_bridge_v1_overtop_%j.out
#SBATCH --error=sbatch-output/error/antidunes_test_bridge_v1_overtop_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# Rerun of the overtop bridge variant (job 4064455's overtop half, which already confirmed clean
# sustained overtopping) with two geometry/output fixes applied to test_bridge_v1_overtop_gas.prm:
#   - Piers now capped at Z=76 (the deck top) instead of running the full domain height -- they stop
#     flush with the top of the bridge instead of poking up into the open air above it.
#   - Added mapped_flag to the fluid_field VTK writers, so the bridge geometry (piers/deck/arches,
#     all NoSlip=4) can be isolated directly in ParaView independent of the water.
# See the prm's own header for the full water-level derivation (liquidHeightFactor=5.9).
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
    --prm "$REPO/apps/showcases/Antidunes/prms/test_bridge_v1_overtop_gas.prm"

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
