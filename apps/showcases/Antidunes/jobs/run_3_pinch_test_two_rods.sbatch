#!/bin/bash -l
#SBATCH --job-name=antidunes-pinch-test-two-rods
#SBATCH --partition=singlenode
#SBATCH --nodes=1
#SBATCH --ntasks=3
#SBATCH --time=00:20:00
#SBATCH --output=sbatch-output/output/pinch_test_two_rods_%j.out
#SBATCH --error=sbatch-output/error/pinch_test_two_rods_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# [Claude, 2026-09-25] MINIMAL REPRODUCTION, per explicit request ("create a small test case
# that proves the pinching crash is a real thing"). See prms/pinch_test_two_rods_gas.prm's own
# header for full derivation. Two static (FIXED for the whole run) cylinders at the exact
# production pitch/radius (radius 3, X-pitch 7 -> 1 LU clearance), fully submerged, still water,
# 3 ranks (mesa_pd's own minimum block count along a periodic direction when using particles),
# 5000 timesteps, all current diagnostics (Option A+B, SimpleBB fix, precursor/NaN scans) live.
set -euo pipefail

module load gcc/11.5.0
module load openmpi/4.1.3-gcc11.5.0
module load cmake/3.31.8
module load python/3.13-miniforge

if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    REPO="$(cd "$SLURM_SUBMIT_DIR/../../.." && pwd)"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    REPO="$(cd "$SCRIPT_DIR/../../../.." && pwd)"
fi
BUILD_DIR="$REPO/build"
BIN_DIR="$BUILD_DIR/apps/showcases/Antidunes"

for dir in vtk-out eval snapshots; do
    if [[ -d "$BIN_DIR/$dir" ]]; then
        echo "WARNING: removing leftover $BIN_DIR/$dir from an unarchived previous run" >&2
        rm -rf "${BIN_DIR:?}/$dir"
    fi
done

cmake --build "$BUILD_DIR" --target Antidunes -- -j"$(nproc)"

cd "$BIN_DIR"
mkdir -p eval snapshots
mpirun -n 3 ./Antidunes "$REPO/apps/showcases/Antidunes/prms/pinch_test_two_rods_gas.prm"

STAMP="$(date +%Y%m%d_%H%M%S)"
if [[ -n "${HPCVAULT:-}" ]] && mkdir -p "$HPCVAULT/antidunes-archive/${SLURM_JOB_NAME}_$STAMP" 2>/dev/null; then
    ARCHIVE_DIR="$HPCVAULT/antidunes-archive/${SLURM_JOB_NAME}_$STAMP"
else
    ARCHIVE_DIR="$BIN_DIR/archive/$STAMP"
    mkdir -p "$ARCHIVE_DIR"
fi
for dir in vtk-out eval snapshots; do
    if [[ -d "$BIN_DIR/$dir" ]]; then
        mv "$BIN_DIR/$dir" "$ARCHIVE_DIR/$dir"
    fi
done
