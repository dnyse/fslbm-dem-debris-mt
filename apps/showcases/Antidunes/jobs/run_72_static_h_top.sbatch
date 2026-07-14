#!/bin/bash -l
#SBATCH --job-name=antidunes-static-h-top
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=02:00:00
#SBATCH --output=sbatch-output/output/antidunes_static_h_top_%j.out
#SBATCH --error=sbatch-output/error/antidunes_static_h_top_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# build + run: STATIC-H ABLATION — identical to the phantom/dynamic-H A/B baseline
# (sphere_drop_two_pillars) but with the static MEM reference height H raised from the flat
# still-water level (~56.02) to the DOMAIN TOP (z=180) via staticReferenceHeightOverride.
# Isolates the pure reference-offset effect: a closed (fully submerged) body should feel ~nothing,
# an open wetted cap at the waterline takes an unbalanced ~13x-weight DOWNWARD push. Compare
# [FloaterForces]/[FloaterBuoyancy] against RUN A of phantom.log / dynamic.log (H=56.02 control).
# domainSize <128, 60, 180> / cellsPerBlock <32, 20, 30> = 4x3x6 = 72 blocks
REPO="$HOME/mt-dennys"
BIN_DIR="$REPO/build/apps/showcases/Antidunes"

"$REPO/apps/showcases/Antidunes/antidunes.sh" \
    --build --jobs "$(nproc)" \
    --run --np 72 \
    --prm "$REPO/apps/showcases/Antidunes/prms/sphere_drop_static_h_top_gas.prm"

# archive this run's own output so the next job doesn't overwrite it
STAMP="$(date +%Y%m%d_%H%M%S)"
# if $HPCVAULT is set, archive to persistent vault storage; otherwise the local build dir
if [[ -n "$HPCVAULT" ]]; then
    ARCHIVE_DIR="$HPCVAULT/antidunes-archive/${SLURM_JOB_NAME}_$STAMP"
else
    ARCHIVE_DIR="$BIN_DIR/archive/$STAMP"
fi
mkdir -p "$ARCHIVE_DIR"
for dir in vtk-out eval snapshots; do
    [[ -d "$BIN_DIR/$dir" ]] && mv "$BIN_DIR/$dir" "$ARCHIVE_DIR/$dir"
done
