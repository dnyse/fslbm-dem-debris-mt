#!/bin/bash -l
#SBATCH --job-name=antidunes-test-bridge-v1-raisedwater-spheres-v3-logjam-floaterlog
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=03:00:00
#SBATCH --output=sbatch-output/output/antidunes_test_bridge_v1_raisedwater_spheres_v3_logjam_floaterlog_%j.out
#SBATCH --error=sbatch-output/error/antidunes_test_bridge_v1_raisedwater_spheres_v3_logjam_floaterlog_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# V3 of the sphere raised-water run (job 4136666): water level dropped BELOW the springing line
# (liquidHeightFactor 2.9 -> 2.7, ~56.33 -> ~54.73, ~1.27 LU below the springing line at 56)
# since visual review showed the previous level still looked too submerged. Relies on the
# current damming up against the piers/deck to push the local water level near the bridge back
# above the springing line dynamically. See
# prms/test_bridge_v1_raisedwater_spheres_v3_logjam_floaterlog_gas.prm's own header for the full
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
    --prm "$REPO/apps/showcases/Antidunes/prms/test_bridge_v1_raisedwater_spheres_v3_logjam_floaterlog_gas.prm"

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
