#!/bin/bash -l
#SBATCH --job-name=antidunes-cylinder-gravity-comparison
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=02:00:00
#SBATCH --output=sbatch-output/output/antidunes_cylinder_gravity_comparison_%j.out
#SBATCH --error=sbatch-output/error/antidunes_cylinder_gravity_comparison_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# Cylinder analogue of the sphere dynamic-gravity vs. Bogner-box-gravity comparison
# (jobs 3947080 / 3947124), run back to back in a single job on one build:
#   1. cylinder_dynamic_gravity_bogner_gas.prm     -- this project's own Fr-derived gravity
#   2. cylinder_free_flow_bognerboxgravity_gas.prm -- Bogner & Rude (2013)'s floating-cube gravity
# Same domain, same real mean flow, same Bogner three-phase MEM buoyancy, same fully free
# (unpinned) floater and densityRatio 0.5 -- but the floater is now a FloatingCylinder with its
# long axis along X (parallel to the flow, axis <1,0,0>; radius 8, length 32) instead of a sphere.
# timesteps 50000 (up from the sphere runs' 30000), VTK fluid_field writeFrequency 1000 (up from
# 2000). See prms/cylinder_dynamic_gravity_bogner_gas.prm and
# prms/cylinder_free_flow_bognerboxgravity_gas.prm for full rationale.
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

run_case "RUN dynamic gravity (Fr-derived, g~1.169e-5)" \
    cylinder_dynamic_gravity_bogner_gas.prm dynamic_gravity
run_case "RUN Bogner box gravity (g=7.5e-4, tau=1.9)" \
    cylinder_free_flow_bognerboxgravity_gas.prm bogner_box_gravity

echo "## DONE. Archived under $ARCHIVE_ROOT/{dynamic_gravity,bogner_box_gravity}_{vtk-out,eval,snapshots}"
