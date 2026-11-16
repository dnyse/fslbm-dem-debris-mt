#!/bin/bash -l
#SBATCH --job-name=antidunes-test-bridge-v1-overtop-and-lowwater
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=08:00:00
#SBATCH --output=sbatch-output/output/antidunes_test_bridge_v1_overtop_and_lowwater_%j.out
#SBATCH --error=sbatch-output/error/antidunes_test_bridge_v1_overtop_and_lowwater_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# Two water-level variants of the baseline bridge test (test_bridge_v1_gas.prm, job 4060404),
# run back-to-back in one job: same 72-process canal domain, same piers/arches/Re/Fr/We/PID/
# timesteps=30000, only liquidHeightFactor differs.
#
#   1. test_bridge_v1_overtop_gas.prm    (liquidHeightFactor 5.9,  absoluteLiquidHeight~80.3):
#      baseline already piles up water ~9 LU right at the piers (confirmed from job 4060404's
#      own velocity-field data) but never reached the deck top (76). This variant starts the
#      undisturbed water level itself 4.3 LU ABOVE the deck top, so overtopping is expected from
#      the start rather than relying on pile-up alone.
#   2. test_bridge_v1_lowwater_gas.prm   (liquidHeightFactor 2.0,  away-baseline absoluteLiquidHeight~49.1):
#      the baseline run's arches turned out to be close to fully flooded at steady state (near-
#      bridge level 71.5 exceeds even the tallest crown, 66), so this variant LOWERS the whole
#      water column instead, aiming (via the baseline's own ~9 LU pile-up estimate, not re-derived
#      for this shallower regime) to land the free surface genuinely INSIDE the arch tunnels
#      instead of filling them.
#
# See each prm's own header for the full derivation. Builds once (same source for both), then runs
# each prm in turn, archiving each run's own vtk-out/eval/snapshots separately before the next one
# starts (so the second run doesn't overwrite the first's output).
set -euo pipefail

if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    REPO="$(cd "$SLURM_SUBMIT_DIR/../../.." && pwd)"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    REPO="$(cd "$SCRIPT_DIR/../../../.." && pwd)"
fi
BIN_DIR="$REPO/build/apps/showcases/Antidunes"

echo "==> Building Antidunes once (shared by both runs below)"
"$REPO/apps/showcases/Antidunes/antidunes.sh" \
    --build --jobs "$(nproc)"

archive_run() {
    local run_name="$1"
    local stamp
    stamp="$(date +%Y%m%d_%H%M%S)"
    local archive_dir
    if [[ -n "${HPCVAULT:-}" ]]; then
        archive_dir="$HPCVAULT/antidunes-archive/${SLURM_JOB_NAME}-${run_name}_$stamp"
    else
        archive_dir="$BIN_DIR/archive/${run_name}_$stamp"
    fi
    mkdir -p "$archive_dir"
    for dir in vtk-out eval snapshots; do
        [[ -d "$BIN_DIR/$dir" ]] && mv "$BIN_DIR/$dir" "$archive_dir/$dir"
    done
    echo "==> Archived $run_name output to $archive_dir"
}

echo "==> Run 1/2: overtop variant"
"$REPO/apps/showcases/Antidunes/antidunes.sh" \
    --run --np 72 \
    --prm "$REPO/apps/showcases/Antidunes/prms/test_bridge_v1_overtop_gas.prm"
archive_run "overtop"

echo "==> Run 2/2: lowwater variant"
"$REPO/apps/showcases/Antidunes/antidunes.sh" \
    --run --np 72 \
    --prm "$REPO/apps/showcases/Antidunes/prms/test_bridge_v1_lowwater_gas.prm"
archive_run "lowwater"

echo "==> Both runs complete."
