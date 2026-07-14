#!/bin/bash -l
#SBATCH --job-name=antidunes-sphere-held-static
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=02:00:00
#SBATCH --output=sbatch-output/output/antidunes_sphere_held_static_%j.out
#SBATCH --error=sbatch-output/error/antidunes_sphere_held_static_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# STATIC PHANTOM-FORCE TEST (point 1 vs point 2). Two held-fixed spheres in STILL water,
# read the steady MEM [FloaterForces] F_hydro_z at hydrostatic rest:
#   A) fully submerged  (sphere_held_submerged_gas.prm,  centre z=45)
#   B) ~50% submerged   (sphere_held_partial_gas.prm,    centre z=56)
# A≈0 & B>0 -> missing gas pressure on the emerged cap (point 1).
# both >0   -> hydrostatic-staircase error (point 2).   both ≈0 -> no static phantom.
# Both runs in one allocation; grep the .out for [FloaterForces] under each banner below.
# domainSize <128, 60, 180> / cellsPerBlock <32, 20, 30> = 4x3x6 = 72 blocks
REPO="$HOME/mt-dennys"
BIN_DIR="$REPO/build/apps/showcases/Antidunes"
SH="$REPO/apps/showcases/Antidunes/antidunes.sh"
PRM_DIR="$REPO/apps/showcases/Antidunes/prms"
STAMP="$(date +%Y%m%d_%H%M%S)"

# build once (recompiles the useFillWeighting boundary-handling change if not already built)
"$SH" --build --jobs "$(nproc)"

archive_outputs() {  # $1 = label
    for d in vtk-out eval snapshots; do
        [[ -d "$BIN_DIR/$d" ]] && mkdir -p "$BIN_DIR/archive/$STAMP" && mv "$BIN_DIR/$d" "$BIN_DIR/archive/$STAMP/${1}_$d"
    done
}

echo "############################################################"
echo "## RUN A: FULLY SUBMERGED  (sphere_held_submerged_gas.prm) ##"
echo "############################################################"
"$SH" --run --np 72 --prm "$PRM_DIR/sphere_held_submerged_gas.prm"
archive_outputs submerged

echo "############################################################"
echo "## RUN B: ~50%% SUBMERGED   (sphere_held_partial_gas.prm)   ##"
echo "############################################################"
"$SH" --run --np 72 --prm "$PRM_DIR/sphere_held_partial_gas.prm"
archive_outputs partial

echo "## DONE. grep -E 'RUN A|RUN B|FloaterForces' on this .out to compare F_hydro_z."
