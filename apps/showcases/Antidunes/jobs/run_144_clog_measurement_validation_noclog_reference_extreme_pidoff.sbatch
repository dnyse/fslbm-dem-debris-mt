#!/bin/bash -l
#SBATCH --job-name=antidunes-clog-measurement-validation-noclog-reference-extreme-pidoff
#SBATCH --partition=multinode
#SBATCH --nodes=2
#SBATCH --ntasks=144
#SBATCH --time=03:00:00
#SBATCH --output=sbatch-output/output/clog_measurement_validation_noclog_reference_extreme_pidoff_%j.out
#SBATCH --error=sbatch-output/error/clog_measurement_validation_noclog_reference_extreme_pidoff_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# Thesis item 3/4 (THESIS_PLAN.md) -- standalone run of the extreme no-clog reference (PID-off), needed
# to properly baseline the completed extreme_pidoff_gas.prm run (job 4240981, COMPLETED, 115000 steps).
# See run_144_clog_measurement_validation_noclog_reference_extreme_pidon.sbatch's own header for the
# full motivation (elevation-difference metric confounded by a flow-spinup transient of up to +-11 LU
# pre-release, unrelated to clogging) -- this is its PID-off twin, using
# clog_measurement_validation_noclog_reference_extreme_pidoff_gas.prm (PIDParameters.enabled=false,
# otherwise identical to the PID-on reference prm).
#
# 144 ranks / 2 nodes, matching every other job in this "extreme" clog-measurement family. No debris at
# all, so no [RELEASE-OVERLAP] risk; 3h requested for margin (see PID-on sibling's own timing estimate).
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
echo "## clog measurement validation: no-clog reference, extreme, PID off"
echo "############################################################"
"$SH" --run --np 144 --prm "$PRM_DIR/clog_measurement_validation_noclog_reference_extreme_pidoff_gas.prm"

# Prefer $HPCVAULT, decided and verified HERE (right before use), not at job start -- see
# run_144_clog_measurement_validation_extreme_pidoff.sbatch's own 2026-09-16 comment for the full
# rationale (autofs staleness, local-only is a last resort not a default).
STAMP="$(date +%Y%m%d_%H%M%S)"
VAULT_CANDIDATE="${HPCVAULT:-}/antidunes-archive/${SLURM_JOB_NAME}_${SLURM_JOB_ID}_${STAMP}"
if [[ -n "${HPCVAULT:-}" ]] && mkdir -p "$VAULT_CANDIDATE" 2>/dev/null \
   && touch "$VAULT_CANDIDATE/.write_test" 2>/dev/null; then
    rm -f "$VAULT_CANDIDATE/.write_test"
    ARCHIVE_ROOT="$VAULT_CANDIDATE"
else
    ARCHIVE_ROOT="$BIN_DIR/archive/${STAMP}_noclog_reference_extreme_pidoff"
    mkdir -p "$ARCHIVE_ROOT"
    echo "WARNING: \$HPCVAULT unavailable at archive time -- falling back to LOCAL archive at $ARCHIVE_ROOT (move this off manually to avoid quota pressure)" >&2
fi
for d in vtk-out eval snapshots; do
    [[ -d "$BIN_DIR/$d" ]] && mv "$BIN_DIR/$d" "$ARCHIVE_ROOT/$d"
done

echo "## DONE. Archived under $ARCHIVE_ROOT"
