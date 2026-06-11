#!/bin/bash -l
#SBATCH --job-name=antidunes-box-drop
#SBATCH --partition=singlenode
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=02:00:00
#SBATCH --output=sbatch-output/output/antidunes_box_drop_%j.out
#SBATCH --error=sbatch-output/error/antidunes_box_drop_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# build + run: box equilibrium / free-settling test.
# One freely-rotating cube (edge 24 LU, ratio 0.5) released from rest above the surface with a
# 15 deg tilt — observe rest depth in Z (~50% submersion) and the orientation it settles into
# (centre-of-buoyancy righting/capsizing torque, Bogner & Ruede 2013 Fig. 5.4).
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
    --prm "$REPO/apps/showcases/Antidunes/prms/giant_box_drop_gas.prm"
