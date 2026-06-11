#!/bin/bash -l
#SBATCH --job-name=antidunes-giant-box-two-pillars
#SBATCH --partition=singlenode
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=02:00:00
#SBATCH --output=sbatch-output/output/antidunes_giant_box_two_pillars_%j.out
#SBATCH --error=sbatch-output/error/antidunes_giant_box_two_pillars_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# build + run: the giant_sphere scenario with the sphere swapped for the GIANT box
# (12×36×12 LU, halfExtents 6×18×6, 1:3:1 ratio). Placed at (30,30,62) clearing the two
# cylinder pillars (centres y=15,45 r=8), bottom at surface, settle-and-release.
# densityRatio 0.5 → 50% submerged; GM<0 → rolls; body spans 2+ blocks → SyncGhostOwners.
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
    --prm "$REPO/apps/showcases/Antidunes/prms/giant_box_two_pillars_gas.prm"
