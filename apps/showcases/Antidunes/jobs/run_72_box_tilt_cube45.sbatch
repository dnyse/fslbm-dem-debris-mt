#!/bin/bash -l
#SBATCH --job-name=antidunes-box-tilt-cube45
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=03:00:00
#SBATCH --output=sbatch-output/output/antidunes_box_tilt_cube45_%j.out
#SBATCH --error=sbatch-output/error/antidunes_box_tilt_cube45_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# CUBE 45° edge-down test — reproduce the classic righting test with the NEW ANALYTIC buoyancy (exact
# OBB∩plane clip), LONG run so the left↔right rocking damps. Free cube (edge 26, densityRatio 0.5) in still
# water, started 10° off face-up (unstable) → should roll to the STABLE edge-down orientation, tiltY → ±45°.
# 250k steps (was 75k, still rocking). Watch [BoxBuoyancy] tiltY converge to ~±45° and submersion ~50%,
# and the [BoxRightingMoment] M_y → 0 at rest. Analytic-vs-sampling V_sub cross-check should stay <~1%.
# grep -E 'BoxBuoyancy|BoxRightingMoment' on the .out.
# domainSize <128, 60, 180> / cellsPerBlock <32, 20, 30> = 4x3x6 = 72 blocks
REPO="$HOME/mt-dennys"
BIN_DIR="$REPO/build/apps/showcases/Antidunes"
SH="$REPO/apps/showcases/Antidunes/antidunes.sh"
PRM_DIR="$REPO/apps/showcases/Antidunes/prms"
STAMP="$(date +%Y%m%d_%H%M%S)"

"$SH" --build --jobs "$(nproc)"

echo "############################################################"
echo "## CUBE 45° edge-down — analytic buoyancy, long run"
echo "############################################################"
"$SH" --run --np 72 --prm "$PRM_DIR/box_tilt_cube45_analytic_gas.prm"
for d in vtk-out eval snapshots; do
    [[ -d "$BIN_DIR/$d" ]] && mkdir -p "$BIN_DIR/archive/$STAMP" && mv "$BIN_DIR/$d" "$BIN_DIR/archive/$STAMP/cube45_$d"
done

echo "## DONE. grep -E 'BoxBuoyancy|BoxRightingMoment' on this .out; tiltY→±45°, submersion~50%, M_y→0."
