#!/bin/bash -l
#SBATCH --job-name=antidunes-sphere-free-flow-bognerboxgravity
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=01:00:00
#SBATCH --output=sbatch-output/output/antidunes_sphere_free_flow_bognerboxgravity_%j.out
#SBATCH --error=sbatch-output/error/antidunes_sphere_free_flow_bognerboxgravity_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# build + run: direct sibling of run_72_sphere_dynamic_gravity_bogner.sbatch
# (job 3947080) -- identical free-floating densityRatio=0.5 sphere, Bogner
# three-phase MEM buoyancy, and real mean flow (targetMeanVelocityMagnitude
# 0.02) -- except gravity/relaxationRate are now OVERRIDDEN to Bogner & Rude
# (2013)'s own floating-cube validation values (g=7.5e-4, tau=1/1.9=0.526)
# instead of the dynamically computed Fr-derived gravity (~1.169e-5). Domain
# height (40) already matches Bogner's own box-experiment domain, so his g is
# used unscaled at the same Eq. 10 stability margin he himself validated.
# Tests whether the much stronger gravity changes the submersion-oscillation
# behaviour seen in the dynamic-gravity sibling run once the sphere is fully
# free (job 3947080: oscillated in a tight band around the correct 50%
# submerged once unpinned). See prms/sphere_free_flow_bognerboxgravity_gas.prm
# for full rationale. 30000 timesteps.
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
    --prm "$REPO/apps/showcases/Antidunes/prms/sphere_free_flow_bognerboxgravity_gas.prm"

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
