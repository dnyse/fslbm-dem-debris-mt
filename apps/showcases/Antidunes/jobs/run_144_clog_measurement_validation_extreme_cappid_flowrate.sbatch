#!/bin/bash -l
#SBATCH --job-name=antidunes-clog-measurement-validation-extreme-cappid-flowrate
#SBATCH --partition=multinode
#SBATCH --nodes=2
#SBATCH --ntasks=144
#SBATCH --time=03:30:00
#SBATCH --output=sbatch-output/output/clog_measurement_validation_extreme_cappid_flowrate_%j.out
#SBATCH --error=sbatch-output/error/clog_measurement_validation_extreme_cappid_flowrate_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# Centred pillar, extreme scale, capped PID, staggered release, with the ClogMeasurement flow-rate evaluator
# (prms/clog_measurement_validation_extreme_cappid_2node_flowrate_gas.prm): 768 x 40 x 180 channel, pillar R 12
# leaving two 8 LU slits, 100 spheres R 5 in 10 release batches, PID cap 2.65e-6, 124750 steps.
# 144 ranks / 2 nodes; 3.5 h requested (~2.5 h measured for this scenario).
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
"$SH" --run --np 144 --prm "$PRM_DIR/clog_measurement_validation_extreme_cappid_2node_flowrate_gas.prm"

# Archive to $HPCVAULT if it is writable right now (checked here, not at job start: autofs mounts can go stale);
# the local archive under build/ is only a last-resort fallback.
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
