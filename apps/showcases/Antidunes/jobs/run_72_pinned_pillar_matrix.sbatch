#!/bin/bash -l
#SBATCH --job-name=antidunes-pinned-pillar-matrix
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=10:00:00
#SBATCH --output=sbatch-output/output/antidunes_pinned_pillar_matrix_%j.out
#SBATCH --error=sbatch-output/error/antidunes_pinned_pillar_matrix_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# PINNED-SPHERE-AT-PILLAR matrix — the clean discriminator. Sphere held FIXED at x=53,y=30 (phantom
# hot-spot) in developed two-pillar flow (U=0.02), STATIC flat H. Removes both confounds (moving body,
# moving reference). 2 heights x fill-weighting on/off = 4 runs on ONE build:
#   A) HALF-submerged  (z=64) + φ fill-weighting ON   (pinned_pillar_half_fw_gas.prm)
#   B) SUBMERGED       (z=50) + φ fill-weighting ON   (pinned_pillar_submerged_fw_gas.prm)
#   C) HALF-submerged  (z=64) + φ fill-weighting OFF  (pinned_pillar_half_nofw_gas.prm)
#   D) SUBMERGED       (z=50) + φ fill-weighting OFF  (pinned_pillar_submerged_nofw_gas.prm)
# Read the STEADY tail of [FloaterForces] F_hydro_z and [FloaterBuoyancy] localSurfaceZ/submersion.
# Questions: (1) is F_hydro_z ~0 at a fixed config, or a steady residual (artifact/real wake force)?
#            (2) does the φ fill-weighting (ON vs OFF) move the residual — i.e. is it helping or hurting?
# grep -E 'RUN A|RUN B|RUN C|RUN D|FloaterForces|FloaterBuoyancy' on the .out to compare.
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

run_case "RUN A: HALF z=64, fill-weight ON"  pinned_pillar_half_fw_gas.prm        half_fw
run_case "RUN B: SUBMERGED z=50, fill-weight ON"  pinned_pillar_submerged_fw_gas.prm   submerged_fw
run_case "RUN C: HALF z=64, fill-weight OFF" pinned_pillar_half_nofw_gas.prm      half_nofw
run_case "RUN D: SUBMERGED z=50, fill-weight OFF" pinned_pillar_submerged_nofw_gas.prm submerged_nofw

echo "## DONE. grep -E 'RUN A|RUN B|RUN C|RUN D|FloaterForces|FloaterBuoyancy' on this .out to compare."
