#!/bin/bash -l
#SBATCH --job-name=antidunes-mixed-shape-sedimentation-showcase-pillars-27floaters-fastcurrent-density05
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=04:00:00
#SBATCH --output=sbatch-output/output/antidunes_mixed_shape_sedimentation_showcase_pillars_27floaters_fastcurrent_density05_%j.out
#SBATCH --error=sbatch-output/error/antidunes_mixed_shape_sedimentation_showcase_pillars_27floaters_fastcurrent_density05_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# build + run: MIXED-SHAPE, WORKING-DENSITY-ARRANGEMENT, PILLARS variant (see
# prms/mixed_shape_sedimentation_showcase_pillars_27floaters_fastcurrent_density05_gas.prm
# for the full rationale) -- NaN-crash ISOLATION VARIANT ("fast current, uniform low
# density"): comparing jobs 3885722/3885723/3885677/3882075 (all 27-floater, all crashed)
# against the ORIGINAL WORKING 20-floater run (job 3863469) shows every floater drifts
# downstream (+X) under the nonzero current, and the Y=10/Y=50 anchors line up with the two
# CylinderPillars (Y=15/45, radius 8) -- in the working run that floater grazes past the
# pillar, in every crashing 27-floater run it instead gets PINNED in sustained contact
# shortly before the crash. This variant tests whether floater density affects
# pinned-vs-grazing: all 27 floaters are densityRatio 0.5 instead of the normal
# [0.5,1.0,2.0] per-anchor mix, keeping the 27-floater grid, sphere radius 4.5, and
# targetMeanVelocityMagnitude doubled to 0.01. Same two static CylinderPillars (DEM contact
# body, Utility.h createCylinder) + matching BodyBC cylinders (BoundaryParameters, the
# actual fluid-blocking geometry) + walled Y (periodicity <1,0,0>) + settlingSteps 10000 on
# every floater, unchanged from the base 27-floater version.
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
    --prm "$REPO/apps/showcases/Antidunes/prms/mixed_shape_sedimentation_showcase_pillars_27floaters_fastcurrent_density05_gas.prm"

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
