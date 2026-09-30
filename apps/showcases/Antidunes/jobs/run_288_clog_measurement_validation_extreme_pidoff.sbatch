#!/bin/bash -l
#SBATCH --job-name=antidunes-clog-measurement-validation-extreme-pidoff
#SBATCH --partition=multinode
#SBATCH --nodes=4
#SBATCH --ntasks=288
#SBATCH --time=04:00:00
#SBATCH --output=sbatch-output/output/clog_measurement_validation_extreme_pidoff_%j.out
#SBATCH --error=sbatch-output/error/clog_measurement_validation_extreme_pidoff_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# Thesis item 3/4 (THESIS_PLAN.md) -- the fourth cell of the 2x2 clog-measurement matrix (no-clog x
# extreme-clog, PID-on x PID-off). accelerationXOverride=7.801745e-07 (the no-clog PID-on run's own
# converged force, job 4242918) pins a fixed driving force instead of letting PID compensate for the
# clog -- see this prm's own PIDParameters comment for why that specific value and not anything
# derived from this scenario's own PID-on sibling.
#
# [Claude, 2026-09-18]: RENAMED from run_144_..._pidoff.sbatch (144 ranks / 2 nodes) to run_288_...
# (288 ranks / 4 nodes) for the same reason as the PID-on sibling's own rename -- domainSize's Z
# doubled (90 -> 180) to fix a hard ceiling-clipping issue found in the PID-on run's data (see
# clog_measurement_validation_extreme_pidoff_gas.prm's own header), which at unchanged
# cellsPerBlock<64,20,15> doubles numBlocks in Z (6 -> 12) and total ranks (144 -> 288). This job was
# also cancelled and relaunched here specifically because its first attempt (job 4245570) had
# already started under the OLD Z=90 domain when the ceiling issue was found -- cheap to redo this
# early rather than let it run for hours on a value already known to be truncated.
#
# Built with today's logging fix (noisy per-cell ObstacleNormalSweep/ObstacleFillLevelSweep/Bogner
# state selection/CurvatureSweep log statements removed at the source) -- measured ~8x per-step
# speedup on the R2 box-MEM comparison. Per-rank work is unchanged from the old 144-rank job (same
# cellsPerBlock), so its own ~1.4h estimate should still hold; 4h kept for the same margin as before.
# A 2-node/144-rank fallback (coarser Z resolution, same physical domain) exists as
# run_144_clog_measurement_validation_extreme_pidoff_2node.sbatch if 4 nodes aren't available.
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
echo "## clog measurement validation: extreme clog, PID off (4-node/288-rank, Z=180)"
echo "############################################################"
"$SH" --run --np 288 --prm "$PRM_DIR/clog_measurement_validation_extreme_pidoff_gas.prm"

# [Claude, 2026-09-16]: prefer $HPCVAULT, decided and verified HERE (right before use) -- see the
# PID-on sibling job's own comment for the full rationale (autofs staleness, local-only is a last
# resort).
STAMP="$(date +%Y%m%d_%H%M%S)"
VAULT_CANDIDATE="${HPCVAULT:-}/antidunes-archive/${SLURM_JOB_NAME}_${SLURM_JOB_ID}_${STAMP}"
if [[ -n "${HPCVAULT:-}" ]] && mkdir -p "$VAULT_CANDIDATE" 2>/dev/null \
   && touch "$VAULT_CANDIDATE/.write_test" 2>/dev/null; then
    rm -f "$VAULT_CANDIDATE/.write_test"
    ARCHIVE_ROOT="$VAULT_CANDIDATE"
else
    ARCHIVE_ROOT="$BIN_DIR/archive/${STAMP}_extreme_pidoff"
    mkdir -p "$ARCHIVE_ROOT"
    echo "WARNING: \$HPCVAULT unavailable at archive time -- falling back to LOCAL archive at $ARCHIVE_ROOT (move this off manually to avoid quota pressure)" >&2
fi
for d in vtk-out eval snapshots; do
    [[ -d "$BIN_DIR/$d" ]] && mv "$BIN_DIR/$d" "$ARCHIVE_ROOT/$d"
done

echo "## DONE. Archived under $ARCHIVE_ROOT"
