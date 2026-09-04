#!/bin/bash -l
#SBATCH --job-name=antidunes-sphere-sed-100spheres-mixeddensity-realgravity-pass3off
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=04:00:00
#SBATCH --output=sbatch-output/output/antidunes_sphere_sed_100spheres_mixeddensity_realgravity_pass3off_%j.out
#SBATCH --error=sbatch-output/error/antidunes_sphere_sed_100spheres_mixeddensity_realgravity_pass3off_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# Pass-3 A/B comparison, OFF side. Same mixed-density 100-sphere scenario as
# run_72_sphere_sedimentation_100spheres_mixeddensity_realgravity.sbatch, but
# full 50000 timesteps and usePass3Reresolution=false (prm:
# sphere_sedimentation_100spheres_mixeddensity_realgravity_pass3off_gas.prm).
# Run on a separate mirror from the _pass3on_ sibling job so both binaries
# build/run concurrently without sharing build/vtk-out directories.
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
    --prm "$REPO/apps/showcases/Antidunes/prms/sphere_sedimentation_100spheres_mixeddensity_realgravity_pass3off_gas.prm"

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
