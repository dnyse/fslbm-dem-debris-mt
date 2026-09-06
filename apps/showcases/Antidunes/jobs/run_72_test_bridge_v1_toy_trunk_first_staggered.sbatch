#!/bin/bash -l
#SBATCH --job-name=antidunes-toy-trunk-first-staggered
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=03:00:00
#SBATCH --output=sbatch-output/output/antidunes_test_bridge_v1_toy_trunk_first_staggered_%j.out
#SBATCH --error=sbatch-output/error/antidunes_test_bridge_v1_toy_trunk_first_staggered_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# Toy-scale "main trunk released first, then staggered secondary debris" release scheme, before
# applying it to the real-scale pont_ferroviere scenarios. 33 total floaters (3 aim points x [1
# trunk + 10 secondary pieces]), useContactFilterOwnerTiebreak enabled (box x cylinder pairs are
# unconditionally GJK). liquidHeightFactor=2.0 (validated safe tier). The 65000-step validation
# (job 4176072) COMPLETED cleanly; this reruns at 150000 steps for a longer observation window
# after the last release (t=45000). --time bumped 01:30:00 -> 03:00:00 (65000 steps took 37:50;
# 150000 steps scales to ~1h27m, doubled here for margin). See the prm's own header for the full
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
    --prm "$REPO/apps/showcases/Antidunes/prms/test_bridge_v1_toy_trunk_first_staggered_gas.prm"

STAMP="$(date +%Y%m%d_%H%M%S)"
if [[ -n "$HPCVAULT" ]]; then
    ARCHIVE_DIR="$HPCVAULT/antidunes-archive/${SLURM_JOB_NAME}_$STAMP"
else
    ARCHIVE_DIR="$BIN_DIR/archive/$STAMP"
fi
mkdir -p "$ARCHIVE_DIR" || ARCHIVE_DIR="$BIN_DIR/archive/$STAMP"
[[ -d "$ARCHIVE_DIR" ]] || mkdir -p "$ARCHIVE_DIR"
for dir in vtk-out eval snapshots floater-history; do
    if [[ -d "$BIN_DIR/$dir" ]]; then
        mv "$BIN_DIR/$dir" "$ARCHIVE_DIR/$dir"
    fi
done
