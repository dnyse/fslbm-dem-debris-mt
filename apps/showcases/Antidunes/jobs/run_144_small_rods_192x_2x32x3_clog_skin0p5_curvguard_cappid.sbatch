#!/bin/bash -l
#SBATCH --job-name=antidunes-small-rods-192x-2x32x3-clog-skin0p5-curvguard-cappid
#SBATCH --partition=multinode
#SBATCH --nodes=2
#SBATCH --ntasks=144
#SBATCH --time=04:00:00
#SBATCH --output=sbatch-output/output/small_rods_192x_2x32x3_clog_skin0p5_curvguard_cappid_%j.out
#SBATCH --error=sbatch-output/error/small_rods_192x_2x32x3_clog_skin0p5_curvguard_cappid_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# Hopper gap, 192 rods, capped PID (prms/small_rods_192x_2x32x3_clog_skin0p5_curvguard_cappid_gas.prm): 2 across
# x 32 rows x 3 layers of rods R 3, L 18 (axis y) upstream of the 16 LU gap, cylinderContactSkin 0.5 with an
# 8 LU row pitch, guardCurvatureInputs on, PID cap 2.65e-6. 100000 steps, 144 ranks / 2 nodes; 4 h requested
# (the 51000-step uncapped version took ~1.2 h including the build).
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
echo "## 192 rods, CAPPED PID 2.65e-6"
echo "## (doubling-in-X step after the clean 48-rod skin test, job 4282396)"
echo "############################################################"
set +e
"$SH" --run --np 144 --prm "$PRM_DIR/small_rods_192x_2x32x3_clog_skin0p5_curvguard_cappid_gas.prm"
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
    ARCHIVE_ROOT="$BIN_DIR/archive/${STAMP}_small_rods_192x_2x32x3_clog_skin0p5_curvguard_cappid"
    mkdir -p "$ARCHIVE_ROOT"
    echo "WARNING: \$HPCVAULT unavailable at archive time -- falling back to LOCAL archive at $ARCHIVE_ROOT (move this off manually to avoid quota pressure)" >&2
fi
for d in vtk-out eval snapshots; do
    [[ -d "$BIN_DIR/$d" ]] && mv "$BIN_DIR/$d" "$ARCHIVE_ROOT/$d"
done

echo "## DONE. Archived under $ARCHIVE_ROOT"
exit "$RUN_RC"
