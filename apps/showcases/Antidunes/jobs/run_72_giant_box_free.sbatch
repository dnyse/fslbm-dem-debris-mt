#!/bin/bash -l
#SBATCH --job-name=antidunes-giant-box-free
#SBATCH --partition=singlenode
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=02:00:00
#SBATCH --output=sbatch-output/output/antidunes_giant_box_free_%j.out
#SBATCH --error=sbatch-output/error/antidunes_giant_box_free_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# build + run: GIANT free-floating log (12×36×12 LU, halfExtents 6×18×6, 1:3:1 ratio) — box analog
# of the giant sphere (circumscribed half-diagonal ≈ 19.9 ≈ sphere radius 20). Starts x=10, sweeps the
# periodic domain; densityRatio 0.5 → 50% submerged equilibrium; rolls (GM<0), 5° tilt seed.
# Body spans 2+ blocks → auto-selects SyncGhostOwners, same as the giant sphere.
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
    --prm "$REPO/apps/showcases/Antidunes/prms/box_free_gas.prm"
