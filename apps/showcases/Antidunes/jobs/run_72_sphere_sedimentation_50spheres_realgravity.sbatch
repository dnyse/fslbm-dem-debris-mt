#!/bin/bash -l
#SBATCH --job-name=antidunes-sphere-sed-50spheres-realgravity
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=02:00:00
#SBATCH --output=sbatch-output/output/antidunes_sphere_sed_50spheres_realgravity_%j.out
#SBATCH --error=sbatch-output/error/antidunes_sphere_sed_50spheres_realgravity_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# Half-count diagnostic: 50 spheres (dropped Y=24/48 rows, kept Y=12/36),
# same domain/gravity/density scheme as sphere_sedimentation_showcase_gas.prm,
# real 9.81 m/s^2 gravity (gravityOverride=2.13715005840044e-05). Isolates
# impact density (fewer spheres hitting the surface at once) as a variable
# against the 100-sphere crash (job 3996320, ~step 9524/50000). timesteps
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
    --prm "$REPO/apps/showcases/Antidunes/prms/sphere_sedimentation_50spheres_realgravity_gas.prm"

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
