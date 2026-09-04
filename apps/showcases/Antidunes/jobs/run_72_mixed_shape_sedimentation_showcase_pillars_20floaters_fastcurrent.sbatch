#!/bin/bash -l
#SBATCH --job-name=antidunes-mixed-shape-sedimentation-showcase-pillars-20floaters-fastcurrent
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=04:00:00
#SBATCH --output=sbatch-output/output/antidunes_mixed_shape_sedimentation_showcase_pillars_20floaters_fastcurrent_%j.out
#SBATCH --error=sbatch-output/error/antidunes_mixed_shape_sedimentation_showcase_pillars_20floaters_fastcurrent_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# build + run: MIXED-SHAPE, WORKING-DENSITY-ARRANGEMENT, PILLARS variant (see
# prms/mixed_shape_sedimentation_showcase_pillars_20floaters_fastcurrent_gas.prm for the
# full rationale) -- NaN-crash ISOLATION VARIANT ("faster current, original 20-floater
# grid"): run_72_mixed_shape_sedimentation_showcase_pillars_27floaters.sbatch filled the
# 20-floater grid out to 27, enlarged spheres 1.25x, AND doubled
# targetMeanVelocityMagnitude all at once and blew up with a global fluid NaN at t=39000
# (mixed_27f.log). The other isolation variants (_27floaters_fastcurrent/_halfstep/_baseline)
# all still carry the 20->27 floater-count change; this one doesn't -- it goes back to the
# ORIGINAL 20-floater grid (7 sphere + 7 box + 6 cylinder) and original sphere radius 4.5,
# changing ONLY targetMeanVelocityMagnitude (doubled to 0.01), the purest possible test of
# the faster current alone. Same two static CylinderPillars (DEM contact body, Utility.h
# createCylinder) + matching BodyBC cylinders (BoundaryParameters, the actual fluid-blocking
# geometry) + walled Y (periodicity <1,0,0>) + settlingSteps 10000 on every floater,
# unchanged from the base 27-floater version.
# domainSize <128, 60, 180> / cellsPerBlock <32, 20, 30> = 4x3x6 = 72 processes
# (unchanged from every other showcase in this project)
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

"$REPO/apps/showcases/Antidunes/antidunes.sh" \
    --build --jobs "$(nproc)" \
    --run --np 72 \
    --prm "$REPO/apps/showcases/Antidunes/prms/mixed_shape_sedimentation_showcase_pillars_20floaters_fastcurrent_gas.prm"

# archive this run's own output so the next job doesn't overwrite it
STAMP="$(date +%Y%m%d_%H%M%S)"
# if $HPCVAULT is set, archive to persistent vault storage; otherwise the local build dir
if [[ -n "$HPCVAULT" ]]; then
    ARCHIVE_DIR="$HPCVAULT/antidunes-archive/${SLURM_JOB_NAME}_$STAMP"
else
    ARCHIVE_DIR="$BIN_DIR/archive/$STAMP"
fi
mkdir -p "$ARCHIVE_DIR"
for dir in vtk-out eval snapshots; do
    [[ -d "$BIN_DIR/$dir" ]] && mv "$BIN_DIR/$dir" "$ARCHIVE_DIR/$dir"
done
