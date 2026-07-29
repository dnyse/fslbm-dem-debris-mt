#!/bin/bash -l
#SBATCH --job-name=antidunes-sphere-density-showcase-bogner-domain-matched
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=00:45:00
#SBATCH --output=sbatch-output/output/antidunes_sphere_density_showcase_bogner_domain_matched_%j.out
#SBATCH --error=sbatch-output/error/antidunes_sphere_density_showcase_bogner_domain_matched_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# build + run: same Bogner three-phase-MEM/buoyancy model as
# run_72_sphere_density_showcase_bogner.sbatch (job 3905990, which converged
# stably but ~8-12% under-submerged at all three densities), except gravity
# and viscosity are now overridden to match Bogner & Rude (2013)'s own
# Fekken-cube validation exactly (g=7.5e-4, tau=1/1.9=0.526, thesis section
# 5.3.3), instead of Antidunes' own Fr-derived antidune-flow scaling
# (g=7.02618e-06, ~107x smaller). Domain reshaped to <96,45,40> -- the closest
# footprint to Bogner's own <130,40,40> that is compatible with the existing
# spheres_out.dat bed file's fixed X:Y generation ratio (see the prm header
# for why the ratio can't be freely chosen). Sphere radius reduced 12->8 to
# match Bogner's own cube side length s=16. See
# prms/sphere_density_showcase_bogner_domain_matched_gas.prm for full
# rationale. 30000 timesteps (equilibrium was already established within the
# first ~800 steps in the prior run).
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
    --prm "$REPO/apps/showcases/Antidunes/prms/sphere_density_showcase_bogner_domain_matched_gas.prm"

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
