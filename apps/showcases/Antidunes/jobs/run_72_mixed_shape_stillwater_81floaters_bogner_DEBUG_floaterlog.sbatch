#!/bin/bash -l
#SBATCH --job-name=antidunes-mixed-stillwater-81f-DEBUG-floaterlog
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=00:30:00
#SBATCH --output=sbatch-output/output/antidunes_mixed_shape_stillwater_81floaters_bogner_DEBUG_floaterlog_%j.out
#SBATCH --error=sbatch-output/error/antidunes_mixed_shape_stillwater_81floaters_bogner_DEBUG_floaterlog_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# Reproducing job 4101985's [ContactHistory DUPLICATE] with logFloaterHistory on this time, to
# check whether it's fatal with current code (that historical job reported COMPLETED despite
# logging the duplicate -- unexplained, possibly an exit-code propagation quirk) and, if it fires,
# get the fine-grained state/contact CSVs for uid=2017612633061982216 (a sphere at <16,50,64>,
# nowhere near any periodic or block boundary -- direct evidence this ContactFilter race isn't
# periodic-wrap-specific).
set -euo pipefail

if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    REPO="$(cd "$SLURM_SUBMIT_DIR/../../.." && pwd)"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    REPO="$(cd "$SCRIPT_DIR/../../../.." && pwd)"
fi
BIN_DIR="$REPO/build/apps/showcases/Antidunes"

"$REPO/apps/showcases/Antidunes/antidunes.sh" \
    --build --jobs "$(nproc)" \
    --run --np 72 \
    --prm "$REPO/apps/showcases/Antidunes/prms/mixed_shape_stillwater_81floaters_bogner_DEBUG_floaterlog_gas.prm"

STAMP="$(date +%Y%m%d_%H%M%S)"
if [[ -n "$HPCVAULT" ]]; then
    ARCHIVE_DIR="$HPCVAULT/antidunes-archive/${SLURM_JOB_NAME}_$STAMP"
else
    ARCHIVE_DIR="$BIN_DIR/archive/$STAMP"
fi
mkdir -p "$ARCHIVE_DIR" || ARCHIVE_DIR="$BIN_DIR/archive/$STAMP"
[[ -d "$ARCHIVE_DIR" ]] || mkdir -p "$ARCHIVE_DIR"
for dir in vtk-out eval snapshots floater-history; do
    if [[ -d "$BIN_DIR/$dir" ]]; then
        mv "$BIN_DIR/$dir" "$ARCHIVE_DIR/$dir"
    fi
done
