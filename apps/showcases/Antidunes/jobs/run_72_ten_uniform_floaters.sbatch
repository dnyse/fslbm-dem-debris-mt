#!/bin/bash -l
#SBATCH --job-name=antidunes-10unif-floaters
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=02:00:00
#SBATCH --output=sbatch-output/output/antidunes_ten_uniform_floaters_%j.out
#SBATCH --error=sbatch-output/error/antidunes_ten_uniform_floaters_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# build + run: 10 uniform floating spheres (r=8, densityRatio 0.5 — same size as the
# single sphere in two_cylinder_floater_gas.prm) released upstream (low X) of the two
# cylinder pillars at x=64; swept downstream by the +X flow and jammed against / between
# the piers (multi-body clogging, uniform size).
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
    --prm "$REPO/apps/showcases/Antidunes/prms/ten_floaters_uniform_pillar_gas.prm"
