#!/bin/bash -l
#SBATCH --job-name=antidunes-ten-cuboids-slitpass
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=02:00:00
#SBATCH --output=sbatch-output/output/antidunes_ten_cuboids_slitpass_%j.out
#SBATCH --error=sbatch-output/error/antidunes_ten_cuboids_slitpass_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# Minimal, non-ASan, non-instrumented reproduction of the fluid rho->NaN / invalid-quaternion crash:
# only 10 floating boxes (vs the full scenario's 100), otherwise identical (same domain/pillars/
# timing/units). See prms/ten_cuboids_two_pillars_slitpass_gas.prm's header for the full derivation
# and why this reduction was created (an ASan RelWithDebInfo run of the 100-box version confirmed
# the SimpleBB.h memory-safety fix holds, but the scenario still crashed via the same
# invalid-quaternion assertion afterward -- this run is for fast, cheap iteration on THAT remaining
# bug without ASan's overhead or 100 floaters' worth of log volume).
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
    --prm "$REPO/apps/showcases/Antidunes/prms/ten_cuboids_two_pillars_slitpass_gas.prm"

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
