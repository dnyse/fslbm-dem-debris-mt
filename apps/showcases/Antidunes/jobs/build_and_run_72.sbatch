#!/bin/bash -l
#SBATCH --job-name=antidunes-build-run-72
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=10:00:00
#SBATCH --output=sbatch-output/output/antidunes_build_run_72_%j.out
#SBATCH --error=sbatch-output/error/antidunes_build_run_72_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# configure + build + bed generation + run in one job
# domainSize <128, 60, 120> / cellsPerBlock <32, 20, 20> = 4x3x6 = 72 blocks
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
    --configure \
    --build --jobs "$(nproc)" \
    --bed --bed-np 9 \
    --run --np 72
