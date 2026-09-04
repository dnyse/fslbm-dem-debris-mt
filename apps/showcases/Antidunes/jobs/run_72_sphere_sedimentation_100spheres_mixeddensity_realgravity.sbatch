#!/bin/bash -l
#SBATCH --job-name=antidunes-sphere-sed-100spheres-mixeddensity-realgravity
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=02:00:00
#SBATCH --output=sbatch-output/output/antidunes_sphere_sed_100spheres_mixeddensity_realgravity_%j.out
#SBATCH --error=sbatch-output/error/antidunes_sphere_sed_100spheres_mixeddensity_realgravity_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# Mixed-density diagnostic: full 100 spheres, same domain/gravity as
# sphere_sedimentation_showcase_gas.prm, but densities reassigned via
# densities[(i+j+k)%3] (v2 3D-mixed scheme) instead of v3's column-uniform
# densities[i%3] -- breaks up the spatially-coherent same-density impact
# front to test whether that coherence (not the gravity value itself) drove
# the 100-sphere real-gravity crash (job 3996320, ~step 9524/50000). Real
# 9.81 m/s^2 gravity (gravityOverride=2.13715005840044e-05). timesteps
# reduced 50000 -> 30000. See prm file header for full derivation.
#
# set -e so a fatal in-simulation crash (nonzero exit from antidunes.sh)
# actually fails this SLURM job instead of being masked by the archiving
# step's own clean exit -- sacct reporting COMPLETED/ExitCode=0:0 on a run
# that crashed mid-simulation is a real, previously-hit failure mode here.
set -euo pipefail

# domainSize <128, 60, 180> / cellsPerBlock <32, 20, 30> = 4x3x6 = 72 blocks
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
    --prm "$REPO/apps/showcases/Antidunes/prms/sphere_sedimentation_100spheres_mixeddensity_realgravity_gas.prm"

STAMP="$(date +%Y%m%d_%H%M%S)"
if [[ -n "$HPCVAULT" ]]; then
    ARCHIVE_DIR="$HPCVAULT/antidunes-archive/${SLURM_JOB_NAME}_$STAMP"
else
    ARCHIVE_DIR="$BIN_DIR/archive/$STAMP"
fi
mkdir -p "$ARCHIVE_DIR"
for dir in vtk-out eval snapshots; do
    [[ -d "$BIN_DIR/$dir" ]] && mv "$BIN_DIR/$dir" "$ARCHIVE_DIR/$dir"
done
