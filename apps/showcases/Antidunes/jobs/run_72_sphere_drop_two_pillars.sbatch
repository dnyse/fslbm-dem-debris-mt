#!/bin/bash -l
#SBATCH --job-name=antidunes-sphere-drop-2pillars
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=03:00:00
#SBATCH --output=sbatch-output/output/antidunes_sphere_drop_2pillars_%j.out
#SBATCH --error=sbatch-output/error/antidunes_sphere_drop_2pillars_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# Phantom-force / unknown-PDF diagnostic run.
# Normal sphere (r=10) starts ~5 LU above the surface, held FIXED for 30000 steps, then
# released. Two cylinder pillars (64,15)/(64,45), r=8. Hydro-force submersion scaling is
# DISABLED in ParticleKernels.h (raw MEM force). Watch the log for [FORMERMO-CONVERT] and
# [FORMERMO-UNKNOWNS] (count of unknown cells + how clustered their unknown neighbours are).
# domainSize <128, 60, 180> / cellsPerBlock <32, 20, 30> = 4x3x6 = 72 blocks
REPO="$HOME/mt-dennys"
BIN_DIR="$REPO/build/apps/showcases/Antidunes"

"$REPO/apps/showcases/Antidunes/antidunes.sh" \
    --build --jobs "$(nproc)" \
    --run --np 72 \
    --prm "$REPO/apps/showcases/Antidunes/prms/sphere_drop_two_pillars_gas.prm"

# archive this run's own output so the next job doesn't overwrite it
STAMP="$(date +%Y%m%d_%H%M%S)"
ARCHIVE_DIR="$BIN_DIR/archive/$STAMP"
mkdir -p "$ARCHIVE_DIR"
for dir in vtk-out eval snapshots; do
    [[ -d "$BIN_DIR/$dir" ]] && mv "$BIN_DIR/$dir" "$ARCHIVE_DIR/$dir"
done
