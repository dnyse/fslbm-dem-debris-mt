#!/bin/bash -l
#SBATCH --job-name=antidunes-two-pillars-ab
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=06:00:00
#SBATCH --output=sbatch-output/output/antidunes_two_pillars_ab_%j.out
#SBATCH --error=sbatch-output/error/antidunes_two_pillars_ab_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# STRONG-PHANTOM A/B on the two-pillars trap, both runs on ONE build:
#   A) useHydrostaticReconstruction FALSE (sphere_drop_two_pillars_gas.prm)          -> baseline
#   B) useHydrostaticReconstruction TRUE  (sphere_drop_two_pillars_hydrorecon_gas.prm) -> fix
# Baseline: F_hydro_z +0.017..+0.044, sphere trapped at 2-27% submersion, never settles.
# If the depth-consistent reconstruction is the right lever, B should break free -> ~50%.
# grep -E 'RUN A|RUN B|FloaterForces|RECON-DENSITY-CHECK' on the .out to compare.
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

# build once (picks up the useHydrostaticReconstruction code if not already built)
"$SH" --build --jobs "$(nproc)"

archive_outputs() {  # $1 = label
    for d in vtk-out eval snapshots; do
        [[ -d "$BIN_DIR/$d" ]] && mkdir -p "$BIN_DIR/archive/$STAMP" && mv "$BIN_DIR/$d" "$BIN_DIR/archive/$STAMP/${1}_$d"
    done
}

echo "############################################################"
echo "## RUN A: BASELINE  (two_pillars_gas.prm, fix OFF)        ##"
echo "############################################################"
"$SH" --run --np 72 --prm "$PRM_DIR/sphere_drop_two_pillars_gas.prm"
archive_outputs baseline

echo "############################################################"
echo "## RUN B: FIX       (two_pillars_hydrorecon_gas.prm, ON)  ##"
echo "############################################################"
"$SH" --run --np 72 --prm "$PRM_DIR/sphere_drop_two_pillars_hydrorecon_gas.prm"
archive_outputs hydrorecon

echo "## DONE. grep -E 'RUN A|RUN B|FloaterForces' on this .out to compare submersion & F_hydro_z."
