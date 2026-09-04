#!/bin/bash -l
#SBATCH --job-name=antidunes-thousand-spheres-two-pillars-TEST-0floaters
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=01:00:00
#SBATCH --output=sbatch-output/output/antidunes_thousand_spheres_two_pillars_TEST_0floaters_%j.out
#SBATCH --error=sbatch-output/error/antidunes_thousand_spheres_two_pillars_TEST_0floaters_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# Isolation test 0 of 4: same widened domain + pillars as the real scenario, ZERO floating
# spheres, 10000 timesteps. See thousand_spheres_two_pillars_slitpass_TEST_0floaters_gas.prm's
# header for why -- checking whether job 4057248's t=2000 NaN (mid-water, far from any sphere)
# comes from the domain/pillar/gravity setup alone.
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
    --prm "$REPO/apps/showcases/Antidunes/prms/thousand_spheres_two_pillars_slitpass_TEST_0floaters_gas.prm"

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
