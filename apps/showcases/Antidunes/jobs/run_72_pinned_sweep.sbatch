#!/bin/bash -l
#SBATCH --job-name=antidunes-pinned-sweep
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=12:00:00
#SBATCH --output=sbatch-output/output/antidunes_pinned_sweep_%j.out
#SBATCH --error=sbatch-output/error/antidunes_pinned_sweep_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# PHANTOM-FORMULA VALIDATION SWEEP — does F_pred = ΔP·A_proj (A_proj∝r²) hold across sphere size and
# submersion? Pinned sphere at the pillar hot-spot (x=53,y=30) in the two-pillar pile-up (U=0.02),
# static flat H; each run logs live [PhantomDecomp] F_pred_z next to measured [FloaterForces] F_hydro_z.
# 5 runs on ONE build, sharing the r=8 / ~50% pivot:
#   RADIUS axis (~50% submerged):  r=5, r=8(pivot), r=10   -> tests the r² law
#   SUBMERSION axis (r=8):         ~25%, ~50%(pivot), ~75% -> tests the A_proj(submersion) shape
# Density is irrelevant for a PINNED sphere (MEM force depends on geometry+fluid only), so it is fixed.
# PASS = ratio meas/pred ≈ 1 across ALL five (mean over the steady tail, steps ≥ 20000).
# Compare: grep -E 'RUN |PhantomDecomp|FloaterForces' on the .out.
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

run_case "RUN r05_sub50: RADIUS axis, r=5, ~50%"    pinned_sweep_r05_sub50_gas.prm r05_sub50
run_case "RUN r08_sub50: PIVOT,       r=8, ~50%"    pinned_sweep_r08_sub50_gas.prm r08_sub50
run_case "RUN r10_sub50: RADIUS axis, r=10, ~50%"   pinned_sweep_r10_sub50_gas.prm r10_sub50
run_case "RUN r08_sub25: SUBMERSION axis, r=8, ~25%" pinned_sweep_r08_sub25_gas.prm r08_sub25
run_case "RUN r08_sub75: SUBMERSION axis, r=8, ~75%" pinned_sweep_r08_sub75_gas.prm r08_sub75

echo "## DONE. grep -E 'RUN |PhantomDecomp|FloaterForces' on this .out. PASS = meas/pred ≈ 1 for all 5."
