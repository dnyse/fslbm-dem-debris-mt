#!/bin/bash -l
#SBATCH --job-name=antidunes-sphere-drop-nophi
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=03:00:00
#SBATCH --output=sbatch-output/output/antidunes_sphere_drop_nophi_%j.out
#SBATCH --error=sbatch-output/error/antidunes_sphere_drop_nophi_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# Phantom-force ABLATION reference: useFillWeighting false -> SimpleBB applies the RAW two-phase
# MEM force (no Eq. 5.9 fill-weighting). submFrac scaling is also off, so this run has NO phantom
# correction at all. Identical otherwise to the slow-flow / no-pillar scenario 2. Compare F_hydro_z
# vs subm and the settled submersion against the phi-weighted run to see what phi-weighting buys.
# This changes the boundary-handling source, so the build will recompile.
# domainSize <128, 60, 180> / cellsPerBlock <32, 20, 30> = 4x3x6 = 72 blocks
REPO="$HOME/mt-dennys"
BIN_DIR="$REPO/build/apps/showcases/Antidunes"

"$REPO/apps/showcases/Antidunes/antidunes.sh" \
    --build --jobs "$(nproc)" \
    --run --np 72 \
    --prm "$REPO/apps/showcases/Antidunes/prms/sphere_drop_slowflow_nophi_gas.prm"

# archive this run's own output so the next job doesn't overwrite it
STAMP="$(date +%Y%m%d_%H%M%S)"
ARCHIVE_DIR="$BIN_DIR/archive/$STAMP"
mkdir -p "$ARCHIVE_DIR"
for dir in vtk-out eval snapshots; do
    [[ -d "$BIN_DIR/$dir" ]] && mv "$BIN_DIR/$dir" "$ARCHIVE_DIR/$dir"
done
