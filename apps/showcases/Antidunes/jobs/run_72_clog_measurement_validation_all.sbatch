#!/bin/bash -l
#SBATCH --job-name=antidunes-clog-measurement-validation-all
#SBATCH --partition=singlenode
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=06:00:00
#SBATCH --output=sbatch-output/output/clog_measurement_validation_all_%j.out
#SBATCH --error=sbatch-output/error/clog_measurement_validation_all_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# Sanity check of the ClogMeasurement evaluator (Evaluators.h) and of the PIDParameters.enabled constant-forcing
# flag (Forcing.h) on a cheap two-pillar scenario (384 x 40 x 90, pillars R 6 leaving a 6 LU gap, 12 spheres R 4
# released from step 20000 on, once the flow has spun up). All three prms run back-to-back in one job so the
# allocation is requested only once:
#   1. no-clog reference (no floating spheres at all)
#   2. sanity clog, PID enabled (controlled-flume mode)
#   3. sanity clog, PID disabled (constant forcing, fixed hydraulic gradient)
# 72 ranks (cellsPerBlock <64,20,15>, 6 x 2 x 6 blocks): --nodes=1 on singlenode reserves a whole 72-core node
# anyway, so fewer ranks would only leave cores idle and give each rank more cells.
# 6 h requested: a comparable run (21600 cells/rank, 72 ranks) took 56.4 ms per step, so 56500 steps at ~19200
# cells/rank take ~53 min per run, ~2.7 h for all three, doubled for the pillar BodyBC and floater DEM overhead.
set -euo pipefail

# Submit with `sbatch jobs/...` from apps/showcases/Antidunes: SLURM_SUBMIT_DIR is then three levels below the
# repo root.
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

STAMP="$(date +%Y%m%d_%H%M%S)"
if [[ -n "${HPCVAULT:-}" ]] && mkdir -p "$HPCVAULT/antidunes-archive/${SLURM_JOB_NAME}_${SLURM_JOB_ID}_$STAMP" 2>/dev/null; then
    ARCHIVE_ROOT="$HPCVAULT/antidunes-archive/${SLURM_JOB_NAME}_${SLURM_JOB_ID}_$STAMP"
else
    ARCHIVE_ROOT="$BIN_DIR/archive/$STAMP"
    mkdir -p "$ARCHIVE_ROOT"
fi

"$SH" --build --jobs "$(nproc)"

archive_outputs() {  # $1 = label
    for d in vtk-out eval snapshots; do
        [[ -d "$BIN_DIR/$d" ]] && mkdir -p "$ARCHIVE_ROOT/$1" && mv "$BIN_DIR/$d" "$ARCHIVE_ROOT/$1/$d"
    done
}

run_case() {  # $1 = banner, $2 = prm file, $3 = label
    echo "############################################################"
    echo "## $1"
    echo "############################################################"
    "$SH" --run --np 72 --prm "$PRM_DIR/$2"
    archive_outputs "$3"
}

run_case "no-clog reference"        clog_measurement_validation_noclog_reference_gas.prm    noclog_reference
run_case "sanity-clog, PID on"      clog_measurement_validation_sanity_clog_pidon_gas.prm   sanity_clog_pidon
run_case "sanity-clog, PID off"     clog_measurement_validation_sanity_clog_pidoff_gas.prm  sanity_clog_pidoff

echo "## DONE. Archived under $ARCHIVE_ROOT/{noclog_reference,sanity_clog_pidon,sanity_clog_pidoff}_{vtk-out,eval,snapshots}"
