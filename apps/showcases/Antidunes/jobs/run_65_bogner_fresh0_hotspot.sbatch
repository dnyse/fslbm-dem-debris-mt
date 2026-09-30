#!/bin/bash -l
#SBATCH --job-name=antidunes-bogner-fresh0-hotspot
#SBATCH --partition=singlenode
#SBATCH --nodes=1
#SBATCH --ntasks=65
#SBATCH --time=0:30:00
#SBATCH --output=sbatch-output/output/bogner_fresh0_hotspot_%j.out
#SBATCH --error=sbatch-output/error/bogner_fresh0_hotspot_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# [Claude, 2026-09-24] Companion to run_65_bogner_resume400_hotspot.sbatch: same 26.565deg-seed
# rho=0.7665 scenario, same widened hot-spot VTK region (writeFrequency=1), but run completely FRESH
# from t=0 (loadSnapshot false) instead of resuming from the step_400 checkpoint. Run alongside the
# resume version on a different checkout so both a from-scratch and a resumed attempt exist side by
# side. timesteps capped at 1000 (safe margin past the step-444 crash point seen in 3 of 4 prior
# resume attempts) to keep the per-step VTK I/O bounded.
# Queue: independent, confirmed clean/free before submission.
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
echo "## Bogner box, s=16, rho=0.7665, tilt=26.565deg, FRESH from t=0, hot-spot VTK region"
echo "############################################################"
set +e
"$SH" --run --np 65 --prm "$PRM_DIR/bogner_exact_3box_s16_rho0767_tilt26p565_fresh0_hotspot_gas.prm"
RC=$?
set -e
echo "## run exited with status $RC (a crash around step 444 would MATCH the resume-run pattern; a clean finish would NOT)"

STAMP="$(date +%Y%m%d_%H%M%S)"
VAULT_CANDIDATE="${HPCVAULT:-}/antidunes-archive/${SLURM_JOB_NAME}_${SLURM_JOB_ID}_${STAMP}"
if [[ -n "${HPCVAULT:-}" ]] && mkdir -p "$VAULT_CANDIDATE" 2>/dev/null \
   && touch "$VAULT_CANDIDATE/.write_test" 2>/dev/null; then
    rm -f "$VAULT_CANDIDATE/.write_test"
    ARCHIVE_ROOT="$VAULT_CANDIDATE"
else
    ARCHIVE_ROOT="$BIN_DIR/archive/${STAMP}_fresh0_hotspot"
    mkdir -p "$ARCHIVE_ROOT"
    echo "WARNING: \$HPCVAULT unavailable at archive time -- falling back to LOCAL archive at $ARCHIVE_ROOT (move this off manually to avoid quota pressure)" >&2
fi
for d in vtk-out eval snapshots; do
    [[ -d "$BIN_DIR/$d" ]] && mv "$BIN_DIR/$d" "$ARCHIVE_ROOT/$d"
done

echo "## DONE. Archived under $ARCHIVE_ROOT"
exit 0
