#!/bin/bash -l
#SBATCH --job-name=antidunes-thousand-spheres-collapsed-144
#SBATCH --partition=multinode
#SBATCH --constraint=el9
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=72
#SBATCH --time=08:00:00
#SBATCH --output=sbatch-output/output/antidunes_thousand_spheres_collapsed_144_%j.out
#SBATCH --error=sbatch-output/error/antidunes_thousand_spheres_collapsed_144_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# 2-node / 144-core (multinode partition, MinNodes=2 on this cluster) rerun of the 1000-sphere
# radius-2.0 scenario: half the target mean velocity, 25000-step flow-development warmup, all 10
# Z-layers collapsed onto the lowest layer's Z (fully overlapping via startInactive), released one
# layer (~100 spheres) every 2000 steps, 200000 steps total. domainSize UNCHANGED (160,75,180);
# only cellsPerBlock changed (Z: 30->15) to get 4x3x12=144 blocks for the extra node. See
# prms/thousand_spheres_two_pillars_slitpass_collapsed_gas.prm's header for the full derivation.
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
    --run --np 144 \
    --prm "$REPO/apps/showcases/Antidunes/prms/thousand_spheres_two_pillars_slitpass_collapsed_gas.prm"

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
