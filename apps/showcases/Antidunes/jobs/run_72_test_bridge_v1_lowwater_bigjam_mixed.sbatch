#!/bin/bash -l
#SBATCH --job-name=antidunes-test-bridge-v1-lowwater-bigjam-mixed
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=03:00:00
#SBATCH --output=sbatch-output/output/antidunes_test_bridge_v1_lowwater_bigjam_mixed_%j.out
#SBATCH --error=sbatch-output/error/antidunes_test_bridge_v1_lowwater_bigjam_mixed_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# Mixed-shape variant of the confirmed-working bigjam scenario (job 4141536, box-only): same
# domain/physics/piers/arches/ModelParameters, but the 162 floaters are split three ways (54 box +
# 54 sphere + 54 cylinder, equal volume, wave-index-mod-3 cycling per aim point). See
# prms/test_bridge_v1_lowwater_bigjam_mixed_gas.prm's own header for the full derivation.
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
    --prm "$REPO/apps/showcases/Antidunes/prms/test_bridge_v1_lowwater_bigjam_mixed_gas.prm"

STAMP="$(date +%Y%m%d_%H%M%S)"
if [[ -n "$HPCVAULT" ]]; then
    ARCHIVE_DIR="$HPCVAULT/antidunes-archive/${SLURM_JOB_NAME}_$STAMP"
else
    ARCHIVE_DIR="$BIN_DIR/archive/$STAMP"
fi
mkdir -p "$ARCHIVE_DIR"
# NOTE: unlike the original bigjam sbatch, this uses an if-statement (not `[[ -d ]] && mv`) for the
# archival loop -- under `set -e`, a bare `&&` list whose test half is false (e.g. floater-history/
# legitimately not existing, since logFloaterHistory is off) exits the whole script with a nonzero
# code even though the simulation and archival both actually succeeded. This is the exact false-
# "FAILED" alarm pattern hit repeatedly earlier this session (e.g. job 4141536, which completed
# cleanly at step 100000 but sacct still reported FAILED for this reason). An `if` statement's exit
# status doesn't trip `set -e` regardless of the test's truth value, so this reports success
# correctly whether or not floater-history/ exists.
for dir in vtk-out eval snapshots floater-history; do
    if [[ -d "$BIN_DIR/$dir" ]]; then
        mv "$BIN_DIR/$dir" "$ARCHIVE_DIR/$dir"
    fi
done
