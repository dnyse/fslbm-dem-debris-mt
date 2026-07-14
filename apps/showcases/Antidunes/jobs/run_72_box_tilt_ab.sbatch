#!/bin/bash -l
#SBATCH --job-name=antidunes-box-tilt-ab
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=10:00:00
#SBATCH --output=sbatch-output/output/antidunes_box_tilt_ab_%j.out
#SBATCH --error=sbatch-output/error/antidunes_box_tilt_ab_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# BOX ROTATION-EQUILIBRIUM A/B with the new ANALYTICAL buoyancy (OBB∩plane 6-tet clip) + cross-check.
# Uniform cube (edge 26, densityRatio 0.5) in still water, started 10° off face-up. It should roll to
# the STABLE edge-down orientation, tiltY → ~45°.
#   A) PINNED  (pinPosition true)  -> translation held, rotation free: clean tilt-to-45° in place.
#   B) LOOSE   (pinPosition false) -> both free: rolls to ~45° but may drift while tumbling.
# Both use useAnalyticBoxBuoyancy=true; [BoxBuoyancy] prints tiltY and the analytic-vs-sampling V_sub
# cross-check (xcheck V_sub: analytic=… sampling=… diff=…%), expected to agree to <~1%.
# grep -E 'RUN A|RUN B|BoxBuoyancy|xcheck' on the .out.
# domainSize <128, 60, 180> / cellsPerBlock <32, 20, 30> = 4x3x6 = 72 blocks
REPO="$HOME/mt-dennys"
BIN_DIR="$REPO/build/apps/showcases/Antidunes"
SH="$REPO/apps/showcases/Antidunes/antidunes.sh"
PRM_DIR="$REPO/apps/showcases/Antidunes/prms"
STAMP="$(date +%Y%m%d_%H%M%S)"

"$SH" --build --jobs "$(nproc)"

archive_outputs() {  # $1 = label
    for d in vtk-out eval snapshots; do
        [[ -d "$BIN_DIR/$d" ]] && mkdir -p "$BIN_DIR/archive/$STAMP" && mv "$BIN_DIR/$d" "$BIN_DIR/archive/$STAMP/${1}_$d"
    done
}

run_case() {  # $1 = banner, $2 = prm, $3 = label
    echo "############################################################"
    echo "## $1"
    echo "############################################################"
    "$SH" --run --np 72 --prm "$PRM_DIR/$2"
    archive_outputs "$3"
}

run_case "RUN A: PINNED translation, free rotation (box_tilt_pinned)" box_tilt_pinned_gas.prm pinned
run_case "RUN B: LOOSE (both free)                (box_tilt_loose)"   box_tilt_loose_gas.prm  loose

echo "## DONE. grep -E 'RUN A|RUN B|BoxBuoyancy|xcheck' on this .out to read tiltY→45° and the V_sub cross-check."
