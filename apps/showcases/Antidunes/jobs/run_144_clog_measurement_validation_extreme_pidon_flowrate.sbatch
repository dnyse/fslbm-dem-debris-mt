#!/bin/bash -l
#SBATCH --job-name=antidunes-clog-measurement-validation-extreme-pidon-flowrate
#SBATCH --partition=multinode
#SBATCH --nodes=2
#SBATCH --ntasks=144
#SBATCH --time=03:15:00
#SBATCH --output=sbatch-output/output/clog_measurement_validation_extreme_pidon_flowrate_%j.out
#SBATCH --error=sbatch-output/error/clog_measurement_validation_extreme_pidon_flowrate_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# [Claude, 2026-09-22] RERUN of the extreme-clog PID-ON scenario (job 4245629, source of the 'Clog Flow-Forcing Comparison' artifact) with the ClogMeasurement
# flow-rate evaluator turned on (flowRateEnabled true, Evaluators.h flowRateAt -- didn't exist when 4245629 ran). Everything else identical: 144 ranks, extreme clog,
# PID on, Z=180 domain. Time limit 3h15m (4245629 itself COMPLETED in 2h30m56s -- tightened from the original 6h budget to match). Queue: independent, ~/mt-dennys (free).
# See prms/clog_measurement_validation_extreme_pidon_2node_flowrate_gas.prm.
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
echo "## clog measurement validation: extreme clog, PID on (fresh, 2-node/144-rank fallback, Z=180)"
echo "############################################################"
"$SH" --run --np 144 --prm "$PRM_DIR/clog_measurement_validation_extreme_pidon_2node_flowrate_gas.prm"

# Prefer $HPCVAULT, decided and verified HERE (right before use) -- see the 288-rank sibling job's
# own 2026-09-16 comment for the full rationale (autofs staleness, local-only is a last resort).
STAMP="$(date +%Y%m%d_%H%M%S)"
VAULT_CANDIDATE="${HPCVAULT:-}/antidunes-archive/${SLURM_JOB_NAME}_${SLURM_JOB_ID}_${STAMP}"
if [[ -n "${HPCVAULT:-}" ]] && mkdir -p "$VAULT_CANDIDATE" 2>/dev/null \
   && touch "$VAULT_CANDIDATE/.write_test" 2>/dev/null; then
    rm -f "$VAULT_CANDIDATE/.write_test"
    ARCHIVE_ROOT="$VAULT_CANDIDATE"
else
    ARCHIVE_ROOT="$BIN_DIR/archive/${STAMP}_extreme_pidon_2node"
    mkdir -p "$ARCHIVE_ROOT"
    echo "WARNING: \$HPCVAULT unavailable at archive time -- falling back to LOCAL archive at $ARCHIVE_ROOT (move this off manually to avoid quota pressure)" >&2
fi
for d in vtk-out eval snapshots; do
    [[ -d "$BIN_DIR/$d" ]] && mv "$BIN_DIR/$d" "$ARCHIVE_ROOT/$d"
done

echo "## DONE. Archived under $ARCHIVE_ROOT"
