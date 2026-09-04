#!/bin/bash -l
#SBATCH --job-name=antidunes-sphere-stillwater-gravity-comparison-d025
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=02:00:00
#SBATCH --output=sbatch-output/output/antidunes_sphere_stillwater_gravity_comparison_d025_%j.out
#SBATCH --error=sbatch-output/error/antidunes_sphere_stillwater_gravity_comparison_d025_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# densityRatio-0.25 sibling of run_72_sphere_stillwater_gravity_comparison.sbatch (job 3947614,
# densityRatio 0.5): same three-way still-water gravity comparison for a single free-floating
# sphere (radius 8, unpinned in every axis), Bogner three-phase MEM buoyancy throughout, run back
# to back in a single job on one build:
#   1. sphere_stillwater_oldhardcoded_d025_gas.prm      -- this project's own g=7.02618e-06
#   2. sphere_stillwater_bognerboxgravity_d025_gas.prm  -- Bogner & Rude (2013)'s own g=7.5e-4
#   3. sphere_stillwater_fixed_dynamicgravity_d025_gas.prm -- job 3947080's dynamically computed
#      g=1.16888e-05, taken here as a fixed constant
# Archimedes predicts 25% submerged at equilibrium here (vs. 50% for the densityRatio-0.5 set) --
# this run tests whether the gravity-magnitude-dependent bias observed at 0.5 (weaker gravity ->
# bigger deviation from the target fraction) behaves the same way off the halfway point, and
# whether the spherical-cap curve's shallower slope away from 50% submersion changes how visible
# a given position bias is in the submerged-fraction metric. See the three prm files for full
# rationale. 30000 timesteps each.
# domainSize <96, 45, 40> / cellsPerBlock <16, 15, 10> = 6x3x4 = 72 processes

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

# if $HPCVAULT is set, archive to persistent vault storage (same convention as the
# single-run job scripts: ${SLURM_JOB_NAME}_$STAMP); otherwise fall back to the local build
# dir, which is scratch/work storage and not guaranteed to persist.
if [[ -n "$HPCVAULT" ]]; then
    ARCHIVE_ROOT="$HPCVAULT/antidunes-archive/${SLURM_JOB_NAME}_$STAMP"
else
    ARCHIVE_ROOT="$BIN_DIR/archive/$STAMP"
fi

"$SH" --build --jobs "$(nproc)"

archive_outputs() {  # $1 = label
    for d in vtk-out eval snapshots; do
        [[ -d "$BIN_DIR/$d" ]] && mkdir -p "$ARCHIVE_ROOT" && mv "$BIN_DIR/$d" "$ARCHIVE_ROOT/${1}_$d"
    done
}

run_case() {  # $1 = banner, $2 = prm, $3 = label
    echo "############################################################"
    echo "## $1"
    echo "############################################################"
    "$SH" --run --np 72 --prm "$PRM_DIR/$2"
    archive_outputs "$3"
}

run_case "RUN old hardcoded gravity (g=7.02618e-06), densityRatio 0.25" \
    sphere_stillwater_oldhardcoded_d025_gas.prm old_hardcoded
run_case "RUN Bogner box gravity (g=7.5e-4, tau=1.9), densityRatio 0.25" \
    sphere_stillwater_bognerboxgravity_d025_gas.prm bogner_box_gravity
run_case "RUN fixed dynamic-value gravity (g=1.16888e-05), densityRatio 0.25" \
    sphere_stillwater_fixed_dynamicgravity_d025_gas.prm fixed_dynamic_gravity

echo "## DONE. Archived under $ARCHIVE_ROOT/{old_hardcoded,bogner_box_gravity,fixed_dynamic_gravity}_{vtk-out,eval,snapshots}"
