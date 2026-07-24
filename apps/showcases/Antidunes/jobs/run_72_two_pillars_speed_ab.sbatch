#!/bin/bash -l
#SBATCH --job-name=antidunes-two-pillars-speed-ab
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=08:00:00
#SBATCH --output=sbatch-output/output/antidunes_two_pillars_speed_ab_%j.out
#SBATCH --error=sbatch-output/error/antidunes_two_pillars_speed_ab_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# SPEED sweep on the two-pillars trap, both runs on ONE build, STATIC flat H (reverted baseline):
#   A) targetMeanVelocity 0.02  (sphere_drop_two_pillars_gas.prm)      -> fast, current phantom
#   B) targetMeanVelocity 0.005 (sphere_drop_two_pillars_slow_gas.prm) -> 4x slower -> smaller pile-up
# Single-variable change (only the PID target speed). Fills the missing 2x2 cell (pillars x speed).
# Prediction (double-count/pile-up theory): B's phantom shrinks toward the clean ~50% behaviour as the
# pile-up collapses; F_hydro_z should track the measured localSurfaceZ pile-up height. If B keeps the
# full +0.02 phantom, it is the pillars' presence (BC/mapping), not the pile-up.
# grep -E 'RUN A|RUN B|FloaterForces|FloaterBuoyancy|localSurfaceZ' on the .out to compare.
# domainSize <128, 60, 180> / cellsPerBlock <32, 20, 30> = 4x3x6 = 72 blocks
# Resolve the repo root programmatically (mirrors antidunes.sh's own SCRIPT_DIR/REPO_ROOT
# derivation) -- no hardcoded path, so this works across mirror clones (mt-dennys,
# mt-dennys-mirror, mt-dennys-mirror2, ...) regardless of which one submitted the job. Jobs
# are always submitted with `sbatch` from apps/showcases/Antidunes (same dir as antidunes.sh),
# so SLURM_SUBMIT_DIR (guaranteed set by SLURM to that submission directory) sits 3 levels
# below the repo root; fall back to deriving it from this script's own location (four levels
# up from apps/showcases/Antidunes/jobs/) for direct/non-sbatch execution.
if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    REPO="$(cd "$SLURM_SUBMIT_DIR/../../.." && pwd)"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    REPO="$(cd "$SCRIPT_DIR/../../../.." && pwd)"
fi
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

echo "############################################################"
echo "## RUN A: FAST  (U=0.02, two_pillars_gas.prm)             ##"
echo "############################################################"
"$SH" --run --np 72 --prm "$PRM_DIR/sphere_drop_two_pillars_gas.prm"
archive_outputs fast

echo "############################################################"
echo "## RUN B: SLOW  (U=0.005, two_pillars_slow_gas.prm)       ##"
echo "############################################################"
"$SH" --run --np 72 --prm "$PRM_DIR/sphere_drop_two_pillars_slow_gas.prm"
archive_outputs slow

echo "## DONE. grep -E 'RUN A|RUN B|FloaterForces|FloaterBuoyancy' on this .out to compare."
