#!/bin/bash -l
#SBATCH --job-name=antidunes-9stacked
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=02:00:00
#SBATCH --output=sbatch-output/output/antidunes_nine_stacked_%j.out
#SBATCH --error=sbatch-output/error/antidunes_nine_stacked_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# build + run: 9 floating spheres (r=8, densityRatio 0.5) stacked 3x3 in the gas above the
# surface at x~9, all held FIXED and released TOGETHER at t=30k — the column falls as a stack
# (lower rows cushion the upper ones, avoiding a hard far-fall surface crossing) and drifts
# into the two pillars (x=64) as a clogging mass.
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
    --prm "$REPO/apps/showcases/Antidunes/prms/nine_stacked_release_pillar_gas.prm"
