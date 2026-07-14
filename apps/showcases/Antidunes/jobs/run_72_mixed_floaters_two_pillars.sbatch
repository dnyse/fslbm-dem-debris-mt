#!/bin/bash -l
#SBATCH --job-name=antidunes-mixed-floaters-two-pillars
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=03:00:00
#SBATCH --output=sbatch-output/output/antidunes_mixed_floaters_two_pillars_%j.out
#SBATCH --error=sbatch-output/error/antidunes_mixed_floaters_two_pillars_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# build + run: ALL-SHAPE interaction test — one sphere + two boxes + one cylinder log swept
# into two cylinder pillars. Exercises every newly enabled contact pair in one run:
# Cylinder x pillar (log jams the gap), Box x Box (the box pair meets at the centreline, via
# GJK/EPA), Box x Cylinder / Sphere x Cylinder / Sphere x Box (drift pile-up behind the log),
# and Box/Sphere x pillar. Watch F_contact in [BoxForces]/[CylForces]/[FloaterForces].
# domainSize <128, 60, 180> / cellsPerBlock <32, 20, 30> = 4x3x6 = 72 blocks
REPO="$HOME/mt-dennys"
BIN_DIR="$REPO/build/apps/showcases/Antidunes"

"$REPO/apps/showcases/Antidunes/antidunes.sh" \
    --build --jobs "$(nproc)" \
    --run --np 72 \
    --prm "$REPO/apps/showcases/Antidunes/prms/mixed_floaters_two_pillars_gas.prm"

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
