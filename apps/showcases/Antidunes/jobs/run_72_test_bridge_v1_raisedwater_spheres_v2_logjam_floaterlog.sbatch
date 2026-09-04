#!/bin/bash -l
#SBATCH --job-name=antidunes-test-bridge-v1-raisedwater-spheres-v2-logjam-floaterlog
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=03:00:00
#SBATCH --output=sbatch-output/output/antidunes_test_bridge_v1_raisedwater_spheres_v2_logjam_floaterlog_%j.out
#SBATCH --error=sbatch-output/error/antidunes_test_bridge_v1_raisedwater_spheres_v2_logjam_floaterlog_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# V2 of the sphere raised-water run (job 4135949): water level lowered further (liquidHeightFactor
# 3.0 -> 2.9, ~57.13 -> ~56.33) since visual review showed the previous level still submerged too
# much of the arch openings, and sphere count doubled 9 -> 18 (2 per aim point, Z-staggered spawn
# to avoid the spawn-overlap collision bug found in the bigjam box run). See
# prms/test_bridge_v1_raisedwater_spheres_v2_logjam_floaterlog_gas.prm's own header for the full
# derivation.
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
    --prm "$REPO/apps/showcases/Antidunes/prms/test_bridge_v1_raisedwater_spheres_v2_logjam_floaterlog_gas.prm"

STAMP="$(date +%Y%m%d_%H%M%S)"
if [[ -n "$HPCVAULT" ]]; then
    ARCHIVE_DIR="$HPCVAULT/antidunes-archive/${SLURM_JOB_NAME}_$STAMP"
else
    ARCHIVE_DIR="$BIN_DIR/archive/$STAMP"
fi
mkdir -p "$ARCHIVE_DIR"
for dir in vtk-out eval snapshots floater-history; do
    [[ -d "$BIN_DIR/$dir" ]] && mv "$BIN_DIR/$dir" "$ARCHIVE_DIR/$dir"
done
