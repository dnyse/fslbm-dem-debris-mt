#!/bin/bash -l
#SBATCH --job-name=antidunes-single-rod-tumble-cylinder
#SBATCH --partition=singlenode
#SBATCH --nodes=1
#SBATCH --ntasks=3
#SBATCH --time=00:30:00
#SBATCH --output=sbatch-output/output/single_rod_tumble_%j.out
#SBATCH --error=sbatch-output/error/single_rod_tumble_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# [Claude, 2026-09-26] Minimal single-body reproduction, per explicit request ("document our
# instability better for my thesis, have test cases proof this"). ONE cylinder, free to fall and
# tumble from a 45 deg initial tilt, in an otherwise empty still-water domain -- no second rod, no
# flow, no pillars. See prms/single_rod_tumble_gas.prm's own header for full derivation. Tests
# whether the near-total-enclosure precondition behind this session's [MEM-GARBAGE-LINK]/[LOW-RHO]
# crashes arises from a SINGLE rod's own curvature/lattice-alignment alone.
set -euo pipefail

if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    REPO="$(cd "$SLURM_SUBMIT_DIR/../../.." && pwd)"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    REPO="$(cd "$SCRIPT_DIR/../../../.." && pwd)"
fi
SH="$REPO/apps/showcases/Antidunes/antidunes.sh"
PRM_DIR="$REPO/apps/showcases/Antidunes/prms"
BIN_DIR="$REPO/build/apps/showcases/Antidunes"

for dir in vtk-out eval snapshots; do
    if [[ -d "$BIN_DIR/$dir" ]]; then
        echo "WARNING: removing leftover $BIN_DIR/$dir from an unarchived previous run" >&2
        rm -rf "${BIN_DIR:?}/$dir"
    fi
done

"$SH" --build --jobs "$(nproc)"

echo "############################################################"
echo "## single rod tumble test: single_rod_tumble_gas"
echo "############################################################"
set +e
"$SH" --run --np 3 --prm "$PRM_DIR/single_rod_tumble_gas.prm"
RUN_RC=$?
set -e
echo "## run exited with status $RUN_RC"

STAMP="$(date +%Y%m%d_%H%M%S)"
ARCHIVE_ROOT="$BIN_DIR/archive/${STAMP}_single_rod_tumble_gas"
mkdir -p "$ARCHIVE_ROOT"
for d in vtk-out eval snapshots; do
    [[ -d "$BIN_DIR/$d" ]] && mv "$BIN_DIR/$d" "$ARCHIVE_ROOT/$d"
done

echo "## DONE. Archived under $ARCHIVE_ROOT"
exit "$RUN_RC"
