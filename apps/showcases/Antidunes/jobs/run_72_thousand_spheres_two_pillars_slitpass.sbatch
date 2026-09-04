#!/bin/bash -l
#SBATCH --job-name=antidunes-thousand-spheres-two-pillars-slitpass
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=12:00:00
#SBATCH --output=sbatch-output/output/antidunes_thousand_spheres_two_pillars_slitpass_%j.out
#SBATCH --error=sbatch-output/error/antidunes_thousand_spheres_two_pillars_slitpass_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# build + run (v2): 1000 small spheres (radius 1.0, all densityRatio 0.5, dynamically-computed
# gravity -- no gravityOverride), ALL starting outside/above the free surface, in 10 Z-layers
# of 100 spheres each. Each layer's settlingSteps is offset by 1000 from the last (300, 1300,
# ..., 9300) so the total release staggers into 10 events of 100 bodies instead of one event of
# 1000 -- v1 (all 1000 released simultaneously at settlingSteps=300, straddling the surface)
# crashed with a fatal NaN-quaternion assertion at step 321 (job 4056593). Domain widened in
# X/Y (128->160, 60->90; block COUNT unchanged, still 72 processes/one node) with pillars
# pushed further right (X=138, centres Y=23/67, r=8, gaps 15/28/15 LU) and spheres upstream
# (X=[10,100], Y=[5,85]).
#
# set -euo pipefail: v1's script lacked this, so the mid-run crash still let this script reach
# the archive step and exit 0 -- sacct reported COMPLETED/0:0 for a run that actually crashed at
# step 321/100000. Fixed here so a future crash reports as FAILED, matching every other job
# script in this project that guards against exactly this failure mode.
# domainSize <160, 90, 180> / cellsPerBlock <40, 30, 30> = 4x3x6 = 72 blocks
set -euo pipefail

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
    --prm "$REPO/apps/showcases/Antidunes/prms/thousand_spheres_two_pillars_slitpass_gas.prm"

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
