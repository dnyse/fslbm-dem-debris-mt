#!/bin/bash -l
#SBATCH --job-name=antidunes-sphere-waterlevel-sweep
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=01:00:00
#SBATCH --output=sbatch-output/output/antidunes_sphere_waterlevel_sweep_%j.out
#SBATCH --error=sbatch-output/error/antidunes_sphere_waterlevel_sweep_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# WATER-LEVEL SWEEP -- three short (5000-step) runs on ONE build, back to back, testing
# whether a free-floating densityRatio=0.5 sphere's equilibrium submersion depends on the
# absolute water depth once gravity is a FIXED constant (Bogner's gravityOverride=7.5e-4,
# relaxationRateOverride=1.9) instead of this project's own Fr-derived gravity (which is
# inversely tied to liquidHeight by construction: gravity=(U/Fr)^2/liquidHeight).
#
# domainSize <96,45,40> and gravity/viscosity are IDENTICAL across all 3 runs -- only
# liquidHeightFactor (hence liquidHeight, hence where the waterline sits) differs:
#   low   liquidHeight ~14.0
#   mid   liquidHeight ~19.941 (matches sphere_density_showcase_bogner_domain_matched_gas.prm)
#   high  liquidHeight ~26.0
#
# Expected result (basic Archimedes' principle): the sphere settles to ~50% submerged in
# all 3 runs regardless of water depth -- i.e. liquidHeight has NO measurable influence on
# where it floats, once gravity itself no longer depends on liquidHeight. Check the eval
# output / [FloaterBuoyancy] log for each run's converged submersion fraction and compare.
#
# domainSize <96, 45, 40> / cellsPerBlock <16, 15, 10> = 6x3x4 = 72 processes
#
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

run_case() {  # $1 = banner, $2 = prm, $3 = label
    echo "############################################################"
    echo "## $1"
    echo "############################################################"
    "$SH" --run --np 72 --prm "$PRM_DIR/$2"
    archive_outputs "$3"
}

run_case "RUN low: liquidHeight ~14.0" \
    sphere_waterlevel_sweep_low_bogner_gas.prm low
run_case "RUN mid: liquidHeight ~19.941" \
    sphere_waterlevel_sweep_mid_bogner_gas.prm mid
run_case "RUN high: liquidHeight ~26.0" \
    sphere_waterlevel_sweep_high_bogner_gas.prm high

echo "## DONE. Compare the converged submersion fraction ([FloaterBuoyancy] log / eval sqlite)"
echo "## across low/mid/high -- expect ~50% in all three if liquidHeight has no influence."
