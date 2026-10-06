#!/bin/bash -l
#SBATCH --job-name=antidunes-small-rods-24x-2x4x3-clog-skin0p5
#SBATCH --partition=multinode
#SBATCH --nodes=2
#SBATCH --ntasks=144
#SBATCH --time=06:00:00
#SBATCH --output=sbatch-output/output/small_rods_24x_2x4x3_clog_skin0p5_%j.out
#SBATCH --error=sbatch-output/error/small_rods_24x_2x4x3_clog_skin0p5_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# [Claude, 2026-10-01] Rerun of the 24-rod hopper that FAILED as job 4279038 (invalid-quaternion
# NaN cascade, step 18420), with exactly one change: cylinderContactSkin 0.5 (cd1b991b's rod
# contact skin) added, sealed-cell guard left off. Isolates whether the skin alone fixes this
# hopper's cross-rank ghost-ownership desync, before combining it with any other fix. 51000
# timesteps, 6h cap (same as 4279038, which itself used well under this on a clean run).
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

RESCUE="$BIN_DIR/archive/rescued_before_${SLURM_JOB_ID:-manual}_$(date +%Y%m%d_%H%M%S)"
for dir in vtk-out eval snapshots; do
    if [[ -d "$BIN_DIR/$dir" ]]; then
        mkdir -p "$RESCUE"
        echo "WARNING: leftover $BIN_DIR/$dir from a previous run -- moving to $RESCUE (NOT deleting)" >&2
        mv "$BIN_DIR/$dir" "$RESCUE/$dir"
    fi
done

"$SH" --build --jobs "$(nproc)"

echo "############################################################"
echo "## twenty-four small rods (2x4x3) hopper, cylinderContactSkin 0.5, sealed-cell guard OFF"
echo "## (rerun of FAILED job 4279038 with the skin fix alone added)"
echo "############################################################"
set +e
"$SH" --run --np 144 --prm "$PRM_DIR/small_rods_24x_2x4x3_clog_skin0p5_gas.prm"
RUN_RC=$?
set -e
echo "## run exited with status $RUN_RC -- archiving whatever output exists"

STAMP="$(date +%Y%m%d_%H%M%S)"
VAULT_CANDIDATE="${HPCVAULT:-}/antidunes-archive/${SLURM_JOB_NAME}_${SLURM_JOB_ID}_${STAMP}"
if [[ -n "${HPCVAULT:-}" ]] && mkdir -p "$VAULT_CANDIDATE" 2>/dev/null \
   && touch "$VAULT_CANDIDATE/.write_test" 2>/dev/null; then
    rm -f "$VAULT_CANDIDATE/.write_test"
    ARCHIVE_ROOT="$VAULT_CANDIDATE"
else
    ARCHIVE_ROOT="$BIN_DIR/archive/${STAMP}_small_rods_24x_2x4x3_clog_skin0p5"
    mkdir -p "$ARCHIVE_ROOT"
    echo "WARNING: \$HPCVAULT unavailable at archive time -- falling back to LOCAL archive at $ARCHIVE_ROOT (move this off manually to avoid quota pressure)" >&2
fi
for d in vtk-out eval snapshots; do
    [[ -d "$BIN_DIR/$d" ]] && mv "$BIN_DIR/$d" "$ARCHIVE_ROOT/$d"
done

echo "## DONE. Archived under $ARCHIVE_ROOT"
exit "$RUN_RC"
