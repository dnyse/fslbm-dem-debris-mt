#!/bin/bash -l
#SBATCH --job-name=antidunes-pinned-decomp-ab
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=06:00:00
#SBATCH --output=sbatch-output/output/antidunes_pinned_decomp_ab_%j.out
#SBATCH --error=sbatch-output/error/antidunes_pinned_decomp_ab_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# PINNED-SPHERE DECOMPOSITION A/B — re-run of the two heights that matter (fill-weight ON only, since
# the prior matrix proved fill-weighting ON≡OFF), now with the new [PhantomDecomp] log line that prints
# the terms of hypothesis #1 every eval step, right next to the measured [FloaterForces] F_hydro_z:
#   A) HALF-submerged (z=64, open surface) -> expect F_pred_z ≈ +0.02 ≈ measured (phantom present)
#   B) SUBMERGED      (z=50, closed body)  -> expect F_pred_z ≈ 0     (uniform head, no net; measured ~-0.02 secondary)
# [PhantomDecomp] prints: deltaH (localSurfaceZ−flatRef), |a_z|, deltaP (=|a_z|·deltaH), r_contact,
#   A_proj (=π·h_cap·(2r−h_cap)), openSurface, F_pred_z (=deltaP·A_proj).
# Read side-by-side: grep -E 'RUN A|RUN B|PhantomDecomp|FloaterForces' on the .out.
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

run_case "RUN A: HALF z=64 (open surface), fill-weight ON"      pinned_pillar_half_fw_gas.prm      half_fw
run_case "RUN B: SUBMERGED z=50 (closed body), fill-weight ON"  pinned_pillar_submerged_fw_gas.prm submerged_fw

echo "## DONE. grep -E 'RUN A|RUN B|PhantomDecomp|FloaterForces' on this .out to compare F_pred_z vs measured F_hydro_z."
