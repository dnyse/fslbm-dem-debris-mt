#!/bin/bash -l
#SBATCH --job-name=antidunes-thousand-spheres-two-pillars-TEST-100floaters-res2x
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=02:00:00
#SBATCH --output=sbatch-output/output/antidunes_thousand_spheres_two_pillars_TEST_100floaters_res2x_%j.out
#SBATCH --error=sbatch-output/error/antidunes_thousand_spheres_two_pillars_TEST_100floaters_res2x_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# Isolation test 4 of 4: same as TEST_100floaters (100 spheres, radius=1.0, Pass-3 off) but
# domain resolution doubled (X 160->320, Y 75->150), sphere radius kept literally at 1.0. See
# thousand_spheres_two_pillars_slitpass_TEST_100floaters_res2x_gas.prm's header for the full
# derivation and why this differs from the r2/r3 "bigger sphere" tests. 2h budget (4x the cell
# count of the 1h tests).
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
    --prm "$REPO/apps/showcases/Antidunes/prms/thousand_spheres_two_pillars_slitpass_TEST_100floaters_res2x_gas.prm"

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
