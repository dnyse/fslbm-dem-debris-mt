#!/bin/bash -l
#SBATCH --job-name=antidunes-thousand-logs-slitpass-75k
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=04:00:00
#SBATCH --output=sbatch-output/output/antidunes_thousand_logs_slitpass_75k_%j.out
#SBATCH --error=sbatch-output/error/antidunes_thousand_logs_slitpass_75k_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# Precursor test (single node, 75000 steps) before a big multi-node log run, mirroring the
# thousand_spheres_two_pillars_slitpass_collapsed_gas.prm approach but with 1000 floating
# cylinders ("logs": radius 2, length 12, axis <0,1,0>) instead of spheres. Same halved mean
# velocity (0.01), 25000-step warmup, 10 batches of 100 staggered every 2000 steps. See
# prms/thousand_logs_two_pillars_slitpass_gas.prm's header for the full derivation.
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
    --prm "$REPO/apps/showcases/Antidunes/prms/thousand_logs_two_pillars_slitpass_gas.prm"

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
