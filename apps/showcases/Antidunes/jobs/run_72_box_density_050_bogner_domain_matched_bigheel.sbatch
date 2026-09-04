#!/bin/bash -l
#SBATCH --job-name=antidunes-box-density-050-bogner-domain-matched-bigheel
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=01:00:00
#SBATCH --output=sbatch-output/output/antidunes_box_density_050_bogner_domain_matched_bigheel_%j.out
#SBATCH --error=sbatch-output/error/antidunes_box_density_050_bogner_domain_matched_bigheel_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# build + run: same Fekken/Bogner densityRatio-0.5 box equilibration as
# run_72_box_density_050_bogner_domain_matched.sbatch, but seeded with a larger
# initial heel (10 deg instead of 2.86 deg) and only 50000 timesteps -- the original
# 2.86 deg seed decayed straight back to ~0 deg upright and stayed there for the full
# 300k steps, never reaching Bogner's predicted 45 deg heeled equilibrium. densityRatio
# 0.5 is the symmetric/marginal case (0 deg and 45 deg both potential equilibria,
# separated by an unstable point expected around ~22.5 deg); 10 deg is still short of
# that expected barrier, so this checks whether upright's basin of attraction is
# smaller than the naive symmetric estimate. Only need to see which basin the
# trajectory falls into, not full convergence, hence the shorter walltime/timestep count.
# See prms/box_density_050_bogner_domain_matched_bigheel_gas.prm for full rationale.
# domainSize <128, 60, 54> / cellsPerBlock <32, 20, 9> = 4x3x6 = 72 processes
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
    --prm "$REPO/apps/showcases/Antidunes/prms/box_density_050_bogner_domain_matched_bigheel_gas.prm"

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
