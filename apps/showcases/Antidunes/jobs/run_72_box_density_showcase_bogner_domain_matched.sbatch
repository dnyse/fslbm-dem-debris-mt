#!/bin/bash -l
#SBATCH --job-name=antidunes-box-density-showcase-bogner-domain-matched
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=01:00:00
#SBATCH --output=sbatch-output/output/antidunes_box_density_showcase_bogner_domain_matched_%j.out
#SBATCH --error=sbatch-output/error/antidunes_box_density_showcase_bogner_domain_matched_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# build + run: side-by-side combined counterpart to the three separate
# run_72_box_density_0{25,50,75}_bogner_domain_matched.sbatch jobs -- all three
# densities (0.25/0.5/0.75) in ONE domain/run, analogous to
# run_72_sphere_density_showcase_bogner_domain_matched.sbatch and
# run_72_cylinder_density_showcase_bogner_domain_matched.sbatch. Same gravity/
# viscosity override matching Bogner & Rude (2013) Sec. 3.3.1 (g=7.5e-4,
# tau=1/1.9=0.526), same domain <128,60,54> reused verbatim from the
# per-density Bogner-domain-matched box runs (X/Y kept at original scale so
# avgParticleDiameter and box halfExtents 8x16x8 stay unchanged; Z shrunk to
# 54). X positions (32/64/96) copied from the pre-existing non-Bogner combined
# box_density_showcase_gas.prm -- same 32-LU spacing, same clearance margins.
# NOTE: 2.86 deg seed heel (matching Bogner's own paper) -- per this project's
# own per-density testing, densityRatio 0.5 and 0.75 are expected to decay back
# to upright at this small seed (only 0.25 sustains a heel); see
# prms/box_density_showcase_bogner_domain_matched_gas.prm header for the
# follow-up (a 10 deg seed was needed for 0.5 to reach its predicted 45 deg
# branch in isolated testing).
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
    --prm "$REPO/apps/showcases/Antidunes/prms/box_density_showcase_bogner_domain_matched_gas.prm"

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
