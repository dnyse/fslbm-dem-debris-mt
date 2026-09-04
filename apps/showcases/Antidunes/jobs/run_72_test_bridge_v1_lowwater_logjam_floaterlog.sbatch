#!/bin/bash -l
#SBATCH --job-name=antidunes-test-bridge-v1-lowwater-logjam-floaterlog
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=06:00:00
#SBATCH --output=sbatch-output/output/antidunes_test_bridge_v1_lowwater_logjam_floaterlog_%j.out
#SBATCH --error=sbatch-output/error/antidunes_test_bridge_v1_lowwater_logjam_floaterlog_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# Rerun of job 4134495 (log-jam test on the toy bridge) with ModelParameters.logFloaterHistory
# turned on -- that first run completed clean on the fluid side (30000/30000 steps, zero
# LOW-RHO/MIGRATION-NAN/LARGE-OVERLAP) but produced one genuine anomaly: uid
# 1008806316530991107 flew from Z~61 to Z~178.9 near the domain ceiling between steps
# 28000-29000 while approaching the central arch, with 100 [DEM-WATCH]/[STRONG-CONTACT-FORCE]
# hits (both saturated) right before the launch. Only the coarse 1000-step [BoxForces] samples
# caught this last time -- per-subcycle logFloaterHistory should show the actual triggering
# contact directly. Only 9 floaters total, so unfiltered full-run logging is cheap (no
# logFloaterHistoryReleaseStepFilter needed). Everything else identical to job 4134495 --
# see prms/test_bridge_v1_lowwater_logjam_floaterlog_gas.prm's own header.
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
    --prm "$REPO/apps/showcases/Antidunes/prms/test_bridge_v1_lowwater_logjam_floaterlog_gas.prm"

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
