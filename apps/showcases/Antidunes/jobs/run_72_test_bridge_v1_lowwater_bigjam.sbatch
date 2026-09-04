#!/bin/bash -l
#SBATCH --job-name=antidunes-test-bridge-v1-lowwater-bigjam
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=03:00:00
#SBATCH --output=sbatch-output/output/antidunes_test_bridge_v1_lowwater_bigjam_%j.out
#SBATCH --error=sbatch-output/error/antidunes_test_bridge_v1_lowwater_bigjam_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# Big-jam variant of the confirmed-working low-water log-jam scenario (job 4134882): 54 logs
# (6 stacked per upstream aim point, overlapping spawn positions, staggered release) instead
# of 9, 100000 timesteps instead of 30000, floater history CSV logging off, and domain height
# cut from 180 to 90 LU (still 72 blocks/ranks -- cellsPerBlock.z halved to match). See
# prms/test_bridge_v1_lowwater_bigjam_gas.prm's own header for the full derivation of every
# change.
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
    --prm "$REPO/apps/showcases/Antidunes/prms/test_bridge_v1_lowwater_bigjam_gas.prm"

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
