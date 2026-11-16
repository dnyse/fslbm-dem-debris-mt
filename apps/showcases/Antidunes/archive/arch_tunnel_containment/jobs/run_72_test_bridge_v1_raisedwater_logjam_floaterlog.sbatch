#!/bin/bash -l
#SBATCH --job-name=antidunes-test-bridge-v1-raisedwater-logjam-floaterlog
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=06:00:00
#SBATCH --output=sbatch-output/output/antidunes_test_bridge_v1_raisedwater_logjam_floaterlog_%j.out
#SBATCH --error=sbatch-output/error/antidunes_test_bridge_v1_raisedwater_logjam_floaterlog_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# Raised-water variant of the log-jam scenario: the lowwater run (job 4134882, with both the
# effR-vs-exact-corner fix and the springing-line semicircle fix already applied to
# ComputeTunnelContacts) confirmed 0 tunnel-contact rows across the whole 30000-step run and no
# floater ever reaching above Z=55 -- meaning that scenario's water is too low for any log to
# ever reach the real (upper-semicircle) arch structure at all, so it couldn't exercise the
# fixed containment logic against a genuine interaction. This run raises liquidHeightFactor
# 2.0 -> 3.2 (absoluteLiquidHeight ~49.1 -> ~58.73, ~2.7 LU above the springing line Z=56) and
# the floater spawn height 55 -> 60 to match, so the logs should genuinely reach into the arch
# openings this time. See prms/test_bridge_v1_raisedwater_logjam_floaterlog_gas.prm's own header
# for the full derivation.
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
    --prm "$REPO/apps/showcases/Antidunes/prms/test_bridge_v1_raisedwater_logjam_floaterlog_gas.prm"

STAMP="$(date +%Y%m%d_%H%M%S)"
if [[ -n "$HPCVAULT" ]]; then
    ARCHIVE_DIR="$HPCVAULT/antidunes-archive/${SLURM_JOB_NAME}_$STAMP"
else
    ARCHIVE_DIR="$BIN_DIR/archive/$STAMP"
fi
mkdir -p "$ARCHIVE_DIR"
for dir in vtk-out eval snapshots floater-history; do
    [[ -d "$BIN_DIR/$dir" ]] && mv "$BIN_DIR/$dir" "$ARCHIVE_DIR/$dir"
done
