#!/bin/bash -l
#SBATCH --job-name=antidunes-cylinder-stillwater-analytic
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=02:00:00
#SBATCH --output=sbatch-output/output/antidunes_cylinder_stillwater_analytic_%j.out
#SBATCH --error=sbatch-output/error/antidunes_cylinder_stillwater_analytic_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# build + run: direct A/B partner of run_72_cylinder_stillwater.sbatch — identical partial-start
# still-water settle scenario (cylinder radius 7, length 44, axis Y, densityRatio 0.5), single
# variable changed: useAnalyticBoxBuoyancy true, which selects clipCylinderAnalytic() (exact
# disk-plane clip) instead of the 12^3 Cartesian sampling grid that showed a spurious roll-torque
# drift (tiltY away from its 0 deg seed) on this axisymmetric body. Watch [BoxBuoyancy]
# shape=CYLINDER method=ANALYTIC, tiltY (should stay near 0), submersion (should converge ~50%).
# domainSize <128, 60, 180> / cellsPerBlock <32, 20, 30> = 4x3x6 = 72 blocks
REPO="$HOME/mt-dennys"
BIN_DIR="$REPO/build/apps/showcases/Antidunes"

"$REPO/apps/showcases/Antidunes/antidunes.sh" \
    --build --jobs "$(nproc)" \
    --run --np 72 \
    --prm "$REPO/apps/showcases/Antidunes/prms/cylinder_stillwater_analytic_gas.prm"

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
