#!/bin/bash -l
#SBATCH --job-name=antidunes-mixed-shape-sedimentation-showcase-pillars-27floaters-fastcurrent-bogner-computed
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=04:00:00
#SBATCH --output=sbatch-output/output/antidunes_mixed_shape_sedimentation_showcase_pillars_27floaters_fastcurrent_bogner_computed_%j.out
#SBATCH --error=sbatch-output/error/antidunes_mixed_shape_sedimentation_showcase_pillars_27floaters_fastcurrent_bogner_computed_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# build + run: BOGNER MODEL (computed physics) variant of
# run_72_mixed_shape_sedimentation_showcase_pillars_27floaters_fastcurrent.sbatch --
# same mixed-shape 27-floater grid, pillars, walled Y, fast current
# (targetMeanVelocityMagnitude 0.01), domain as the base file, and
# useBognerThreePhaseMEM + useMemOnlyBuoyancy both true (no rhoHydStat subtraction,
# no separate analytic vSub*g buoyancy channel), but WITHOUT the
# gravityOverride/relaxationRateOverride used by
# run_72_mixed_shape_sedimentation_showcase_pillars_27floaters_fastcurrent_bogner.sbatch.
# That full-override variant was found (from a real run screenshot at 99k steps)
# to produce NO sediment buildup in front of either pillar at all: overriding
# gravity/viscosity to Bogner's still-water validation values collapses
# Re 3100->~52 and Fr 1.31->~0.15, which (a) destroys the near-critical Fr
# regime antidune bedforms need to couple the free surface to the bed, and
# (b) raises submerged particle weight ~107x, pushing the critical shear stress
# needed to entrain a grain far above what the flow produces. This variant keeps
# Re/Fr similarity intact (gravity/relaxationRate computed normally) so bedload
# transport/mound formation in front of the pillars remains physically possible,
# while still using the Bogner three-phase MEM force formulation. See
# prms/mixed_shape_sedimentation_showcase_pillars_27floaters_fastcurrent_bogner_computed_gas.prm
# for the full rationale.
# domainSize <128, 60, 180> / cellsPerBlock <32, 20, 30> = 4x3x6 = 72 processes
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
    --prm "$REPO/apps/showcases/Antidunes/prms/mixed_shape_sedimentation_showcase_pillars_27floaters_fastcurrent_bogner_computed_gas.prm"

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
