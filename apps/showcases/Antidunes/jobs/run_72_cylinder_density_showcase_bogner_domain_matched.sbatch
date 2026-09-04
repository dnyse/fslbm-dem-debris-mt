#!/bin/bash -l
#SBATCH --job-name=antidunes-cylinder-density-showcase-bogner-domain-matched
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=00:45:00
#SBATCH --output=sbatch-output/output/antidunes_cylinder_density_showcase_bogner_domain_matched_%j.out
#SBATCH --error=sbatch-output/error/antidunes_cylinder_density_showcase_bogner_domain_matched_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# build + run: same treatment as run_72_sphere_density_showcase_bogner_domain_matched.sbatch
# (job 3908647, converged to within 0.6 percentage points of target density ratios),
# applied to the cylinder/log buoyancy showcase instead of spheres. Same three-phase
# MEM/buoyancy model (useBognerThreePhaseMEM + useMemOnlyBuoyancy), same gravity/viscosity
# override matching Bogner & Rude (2013)'s own validated regime (g=7.5e-4,
# tau=1/1.9=0.526), and the SAME domain <96,45,40> reused verbatim from the sphere run
# (same X:Y generation-ratio lock applies, so it's the same closest-compatible footprint).
# Cylinder radius (8) already matched what the sphere run had to shrink down to, so no
# geometry change was needed there; only the log length (24, along Y) was checked against
# the smaller Y domain (60->45) and still clears the N/S walls by 10.5 LU. See
# prms/cylinder_density_showcase_bogner_domain_matched_gas.prm for full rationale.
# 30000 timesteps, matching the sphere run's convergence timescale at this gravity.
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
    --prm "$REPO/apps/showcases/Antidunes/prms/cylinder_density_showcase_bogner_domain_matched_gas.prm"

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
