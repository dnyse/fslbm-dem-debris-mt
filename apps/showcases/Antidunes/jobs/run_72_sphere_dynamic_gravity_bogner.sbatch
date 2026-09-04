#!/bin/bash -l
#SBATCH --job-name=antidunes-sphere-dynamic-gravity-bogner
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=01:00:00
#SBATCH --output=sbatch-output/output/antidunes_sphere_dynamic_gravity_bogner_%j.out
#SBATCH --error=sbatch-output/error/antidunes_sphere_dynamic_gravity_bogner_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# build + run: single free-floating densityRatio=0.5 sphere, Bogner three-phase
# MEM buoyancy (useBognerThreePhaseMEM + useMemOnlyBuoyancy), but UNLIKE
# run_72_sphere_density_showcase_bogner*.sbatch this run does NOT override
# gravity/relaxationRate -- gravity is left to the project's own dynamic
# Fr-derived formula (gravity = (U/Fr)^2 / liquidHeight), and
# targetMeanVelocityMagnitude is set to 0.02 (E1) so the PID controller is
# actually driving a real mean flow instead of still water. Tests whether the
# sphere still converges to ~50% submerged (basic Archimedes' principle) under
# dynamically-computed gravity plus a real background flow, not just in the
# still-water / gravityOverride cases already checked elsewhere. See
# prms/sphere_dynamic_gravity_bogner_gas.prm for full rationale. 30000
# timesteps.
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

"$REPO/apps/showcases/Antidunes/antidunes.sh" \
    --build --jobs "$(nproc)" \
    --run --np 72 \
    --prm "$REPO/apps/showcases/Antidunes/prms/sphere_dynamic_gravity_bogner_gas.prm"

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
