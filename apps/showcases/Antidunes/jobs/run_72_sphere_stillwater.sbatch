#!/bin/bash -l
#SBATCH --job-name=antidunes-sphere-stillwater
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=02:00:00
#SBATCH --output=sbatch-output/output/antidunes_sphere_stillwater_%j.out
#SBATCH --error=sbatch-output/error/antidunes_sphere_stillwater_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# build + run: STILL-WATER drop test — one floating sphere (r=8, densityRatio 0.5) dropped ~6 LU
# into quiescent water (targetMeanVelocityMagnitude 0, useHardcodedUnits true → valid gravity/
# viscosity/surface tension despite U=0). No flow, no pillars; should settle to 50% submersion.
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
    --prm "$REPO/apps/showcases/Antidunes/prms/sphere_stillwater_gas.prm"
