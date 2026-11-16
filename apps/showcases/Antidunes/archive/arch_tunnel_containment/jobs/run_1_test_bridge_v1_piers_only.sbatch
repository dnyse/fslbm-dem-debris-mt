#!/bin/bash -l
#SBATCH --job-name=antidunes-test-bridge-v1-piers-only
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=01:00:00
#SBATCH --output=sbatch-output/output/antidunes_test_bridge_v1_piers_only_%j.out
#SBATCH --error=sbatch-output/error/antidunes_test_bridge_v1_piers_only_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# First test-bridge shakedown: two round piers, real fluid flow, no deck/arch/floaters yet.
# Single block (1 process) -- small, fast-iterating scale, not the real bridge dimensions. See
# test_bridge_v1_piers_only_gas.prm's header for the BodyBC architectural constraint (no Box/
# Plane shape, no subtraction) that's why there's no arch/deck yet.
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
    --run --np 1 \
    --prm "$REPO/apps/showcases/Antidunes/prms/test_bridge_v1_piers_only_gas.prm"

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
