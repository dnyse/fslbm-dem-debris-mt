#!/bin/bash -l
#SBATCH --job-name=antidunes-sphere-impact-sedimentation-block-atsurface
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=03:00:00
#SBATCH --output=sbatch-output/output/antidunes_sphere_impact_sedimentation_block_atsurface_%j.out
#SBATCH --error=sbatch-output/error/antidunes_sphere_impact_sedimentation_block_atsurface_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# Same 1000-sphere Bogner Fig 5.11 replication block (diagonal density weave)
# as sphere_impact_sedimentation_block_gas.prm, but released essentially AT
# the water surface (0.5 LU gap, vs. the base file's 30 LU free-fall drop) --
# isolates whether the block "punching through" the surface in the dropped
# version is driven by impact momentum built up during the fall, or would
# happen anyway. 20000 timesteps (diagnostic, not the full run). See prm
# file header for full derivation; unchanged from the base file otherwise.
#
# set -e so a fatal in-simulation crash actually fails this SLURM job instead
# of being masked by the archiving step's own clean exit.
set -euo pipefail

# domainSize <150, 141, 240> / cellsPerBlock <50, 47, 30> = 3x3x8 = 72 blocks
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
    --prm "$REPO/apps/showcases/Antidunes/prms/sphere_impact_sedimentation_block_atsurface_gas.prm"

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
