#!/bin/bash -l
#SBATCH --job-name=antidunes-sphere-drop-hydrorecon
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=03:00:00
#SBATCH --output=sbatch-output/output/antidunes_sphere_drop_hydrorecon_%j.out
#SBATCH --error=sbatch-output/error/antidunes_sphere_drop_hydrorecon_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# FIX TEST (b, revised): identical to run_72_sphere_drop_slowflow.sbatch (the (a) baseline) but with
# useHydrostaticReconstruction true -> each FormerMO cell is re-seeded to the depth-consistent
# density rho_hyd(z) = 1 + 3*a_z*(z - H) instead of the fill-blind unweighted mean. Supersedes the
# fill-weighted (sum(phi*rho)/sum(phi)) attempt, which no-op'd because the waterline interface cells
# have negative fill drift (sum(phi) ~ 0). A/B the settled submersion and [FloaterForces] F_hydro_z
# at low submersion against the baseline run.
# domainSize <128, 60, 180> / cellsPerBlock <32, 20, 30> = 4x3x6 = 72 blocks
REPO="$HOME/mt-dennys"
BIN_DIR="$REPO/build/apps/showcases/Antidunes"

"$REPO/apps/showcases/Antidunes/antidunes.sh" \
    --build --jobs "$(nproc)" \
    --run --np 72 \
    --prm "$REPO/apps/showcases/Antidunes/prms/sphere_drop_slowflow_hydrorecon_gas.prm"

# archive this run's own output so the next job doesn't overwrite it
STAMP="$(date +%Y%m%d_%H%M%S)"
ARCHIVE_DIR="$BIN_DIR/archive/$STAMP"
mkdir -p "$ARCHIVE_DIR"
for dir in vtk-out eval snapshots; do
    [[ -d "$BIN_DIR/$dir" ]] && mv "$BIN_DIR/$dir" "$ARCHIVE_DIR/$dir"
done
