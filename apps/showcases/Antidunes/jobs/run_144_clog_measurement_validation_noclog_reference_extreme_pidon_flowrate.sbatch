#!/bin/bash -l
#SBATCH --job-name=antidunes-clog-measurement-validation-noclog-reference-extreme-pidon-flowrate
#SBATCH --partition=multinode
#SBATCH --nodes=2
#SBATCH --ntasks=144
#SBATCH --time=01:45:00
#SBATCH --output=sbatch-output/output/clog_measurement_validation_noclog_reference_extreme_pidon_flowrate_%j.out
#SBATCH --error=sbatch-output/error/clog_measurement_validation_noclog_reference_extreme_pidon_flowrate_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# No-clog reference for the extreme centred-pillar runs, PID on, with the ClogMeasurement flow-rate evaluator
# (prms/clog_measurement_validation_noclog_reference_extreme_flowrate_gas.prm): the same pillar layout without
# floaters. The free-surface spin-up transient (up to +-11 LU between the upstream and downstream probes) is the
# same as in the jam runs, so subtracting this run isolates the debris effect (baseline-referenced afflux).
# 124750 steps, as the jam runs, so the late-time windows line up; domain height Z = 90 (the jam runs use 180).
# 144 ranks / 2 nodes, as the rest of the extreme family; 1 h 45 min requested (~1.25 h measured).
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
echo "## clog measurement validation: no-clog reference, extreme, PID on"
echo "############################################################"
"$SH" --run --np 144 --prm "$PRM_DIR/clog_measurement_validation_noclog_reference_extreme_flowrate_gas.prm"

# Archive to $HPCVAULT if it is writable right now (checked here, not at job start: autofs mounts can go stale);
# the local archive under build/ is only a last-resort fallback.
STAMP="$(date +%Y%m%d_%H%M%S)"
VAULT_CANDIDATE="${HPCVAULT:-}/antidunes-archive/${SLURM_JOB_NAME}_${SLURM_JOB_ID}_${STAMP}"
if [[ -n "${HPCVAULT:-}" ]] && mkdir -p "$VAULT_CANDIDATE" 2>/dev/null \
   && touch "$VAULT_CANDIDATE/.write_test" 2>/dev/null; then
    rm -f "$VAULT_CANDIDATE/.write_test"
    ARCHIVE_ROOT="$VAULT_CANDIDATE"
else
    ARCHIVE_ROOT="$BIN_DIR/archive/${STAMP}_noclog_reference_extreme_pidon"
    mkdir -p "$ARCHIVE_ROOT"
    echo "WARNING: \$HPCVAULT unavailable at archive time -- falling back to LOCAL archive at $ARCHIVE_ROOT (move this off manually to avoid quota pressure)" >&2
fi
for d in vtk-out eval snapshots; do
    [[ -d "$BIN_DIR/$d" ]] && mv "$BIN_DIR/$d" "$ARCHIVE_ROOT/$d"
done

echo "## DONE. Archived under $ARCHIVE_ROOT"
