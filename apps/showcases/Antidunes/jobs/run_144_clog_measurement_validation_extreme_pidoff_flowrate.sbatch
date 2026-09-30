#!/bin/bash -l
#SBATCH --job-name=antidunes-clog-measurement-validation-extreme-pidoff-flowrate
#SBATCH --partition=multinode
#SBATCH --nodes=2
#SBATCH --ntasks=144
#SBATCH --time=03:00:00
#SBATCH --output=sbatch-output/output/clog_measurement_validation_extreme_pidoff_flowrate_%j.out
#SBATCH --error=sbatch-output/error/clog_measurement_validation_extreme_pidoff_flowrate_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# [Claude, 2026-09-18] 2-NODE FALLBACK of run_288_clog_measurement_validation_extreme_pidoff.sbatch
# -- same scenario (extreme clog, PID off, accelerationXOverride=7.801745e-07, Z=180 domain with the
# doubled headroom fix), but using clog_measurement_validation_extreme_pidoff_2node_flowrate_gas.prm
# (cellsPerBlock<64,20,30> instead of <64,20,15> -- half the Z resolution, same physical domain) so
# numBlocks (12,2,6)=144 fits 2 nodes instead of the canonical file's 288-rank/4-node requirement.
# Use this one when only 2 nodes are available on the multinode partition; use the 288-rank job
# otherwise for the full Z resolution.
#
# Time: doubling cellsPerBlock's Z extent doubles the per-rank cell count (same rank count as the OLD
# 144-rank/Z=90 job, which ran in ~53ms/step per that job's own timing) -- roughly 2x the per-step
# cost per rank. The old (pre-headroom-fix) 144-rank PID-on sibling completed in 1h39m (job 4242919);
# 6h requested here for margin against the ~2x per-rank slowdown from the coarser Z.
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
echo "## clog measurement validation: extreme clog, PID off (2-node/144-rank fallback, Z=180)"
echo "############################################################"
"$SH" --run --np 144 --prm "$PRM_DIR/clog_measurement_validation_extreme_pidoff_2node_flowrate_gas.prm"

# Prefer $HPCVAULT, decided and verified HERE (right before use) -- see the 288-rank sibling job's
# own comment for the full rationale (autofs staleness, local-only is a last resort).
STAMP="$(date +%Y%m%d_%H%M%S)"
VAULT_CANDIDATE="${HPCVAULT:-}/antidunes-archive/${SLURM_JOB_NAME}_${SLURM_JOB_ID}_${STAMP}"
if [[ -n "${HPCVAULT:-}" ]] && mkdir -p "$VAULT_CANDIDATE" 2>/dev/null \
   && touch "$VAULT_CANDIDATE/.write_test" 2>/dev/null; then
    rm -f "$VAULT_CANDIDATE/.write_test"
    ARCHIVE_ROOT="$VAULT_CANDIDATE"
else
    ARCHIVE_ROOT="$BIN_DIR/archive/${STAMP}_extreme_pidoff_2node"
    mkdir -p "$ARCHIVE_ROOT"
    echo "WARNING: \$HPCVAULT unavailable at archive time -- falling back to LOCAL archive at $ARCHIVE_ROOT (move this off manually to avoid quota pressure)" >&2
fi
for d in vtk-out eval snapshots; do
    [[ -d "$BIN_DIR/$d" ]] && mv "$BIN_DIR/$d" "$ARCHIVE_ROOT/$d"
done

echo "## DONE. Archived under $ARCHIVE_ROOT"
