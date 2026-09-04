#!/bin/bash -l
#SBATCH --job-name=antidunes-thousand-cuboids-slitpass-75k
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=04:00:00
#SBATCH --output=sbatch-output/output/antidunes_thousand_cuboids_slitpass_75k_%j.out
#SBATCH --error=sbatch-output/error/antidunes_thousand_cuboids_slitpass_75k_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# Box-shaped sibling of run_72_thousand_logs_two_pillars_slitpass.sbatch: same domain/pillars/timing,
# 1000 floating boxes (halfExtents <3,9,3> -- same 1:3 aspect ratio and overall size as the cylinder
# scenario's radius-3/length-18 logs) instead of cylinders. Box x Box uses boxBoxSAT, which is
# provably symmetric under argument-order swap (unlike Cylinder x Cylinder's face-edge/band-edge
# branch), so this is a comparison run to check the mesa_pd::mpi::ContactFilter double-accept race
# root-caused for cylinders does NOT reproduce here. See
# prms/thousand_cuboids_two_pillars_slitpass_gas.prm's header for the full derivation.
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
    --prm "$REPO/apps/showcases/Antidunes/prms/thousand_cuboids_two_pillars_slitpass_gas.prm"

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
