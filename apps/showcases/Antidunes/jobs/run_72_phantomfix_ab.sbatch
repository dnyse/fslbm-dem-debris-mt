#!/bin/bash -l
#SBATCH --job-name=antidunes-phantomfix-ab
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=10:00:00
#SBATCH --output=sbatch-output/output/antidunes_phantomfix_ab_%j.out
#SBATCH --error=sbatch-output/error/antidunes_phantomfix_ab_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# OPTION A end-to-end test — the original two-pillar FLOATER (r=10 sphere, densityRatio 0.5, released at
# step 30000 into developed U=0.02 flow with the pile-up), run twice on ONE build:
#   A) BASELINE  usePhantomForceCorrection OFF (sphere_drop_two_pillars_gas.prm)          -> levitates (phantom)
#   B) FIXED     usePhantomForceCorrection ON  (sphere_drop_two_pillars_phantomfix_gas.prm)-> should sit ~50%
# The correction subtracts the tapered F_pred_z (validated live to ratio ≈ 1) from the VERTICAL MEM force on
# the particle side only (no PDF feedback); self-disables where there is no pile-up. PASS = B settles near
# 50% submersion at the pile-up surface instead of launching up to z≈74 / ~1–9% submersion like A.
# The [FloaterForces] log in B prints the per-step correction: [phantomCorr F_pred_z=… subtracted …].
# grep -E 'RUN A|RUN B|FloaterForces|FloaterBuoyancy|phantom-force correction' on the .out.
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

run_case "RUN A: BASELINE (phantom correction OFF)" sphere_drop_two_pillars_gas.prm          baseline
run_case "RUN B: FIXED    (phantom correction ON)"  sphere_drop_two_pillars_phantomfix_gas.prm phantomfix

echo "## DONE. grep -E 'RUN A|RUN B|FloaterForces|FloaterBuoyancy' on this .out. PASS = B stays ~50% submerged (no levitation)."
