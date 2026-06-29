#!/bin/bash -l
#SBATCH --job-name=antidunes-giant-log-righting
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=02:00:00
#SBATCH --output=sbatch-output/output/antidunes_giant_log_righting_%j.out
#SBATCH --error=sbatch-output/error/antidunes_giant_log_righting_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# build + run: CLEAN righting test for the giant log (12x36x12, densityRatio 0.5). Starts
# HALF-SUBMERGED at z~56 with a 15-deg tilt (no drop/splashdown), settled then released —
# isolates the lowered-COM righting torque (comOffsetZ=-3 -> GM=+2). Watch [BoxBuoyancy]
# tiltY decay 15 -> 0. Set comOffsetZ 0 in the prm for the rolling control case.
# domainSize <128, 60, 180> / cellsPerBlock <32, 20, 30> = 4x3x6 = 72 blocks
REPO="$HOME/mt-dennys"
BIN_DIR="$REPO/build/apps/showcases/Antidunes"

# archive previous output so nothing is overwritten
STAMP="$(date +%Y%m%d_%H%M%S)"
ARCHIVE_DIR="$BIN_DIR/archive/$STAMP"
mkdir -p "$ARCHIVE_DIR"
for dir in vtk-out eval snapshots; do
    [[ -d "$BIN_DIR/$dir" ]] && mv "$BIN_DIR/$dir" "$ARCHIVE_DIR/$dir"
done

"$REPO/apps/showcases/Antidunes/antidunes.sh" \
    --build --jobs "$(nproc)" \
    --run --np 72 \
    --prm "$REPO/apps/showcases/Antidunes/prms/giant_log_righting_gas.prm"
