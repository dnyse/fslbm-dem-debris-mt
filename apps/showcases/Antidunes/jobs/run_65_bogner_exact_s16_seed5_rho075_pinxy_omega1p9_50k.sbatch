#!/bin/bash -l
#SBATCH --job-name=antidunes-bogner-exact-s16-seed5-rho075-pinxy-omega1p9-50k
#SBATCH --partition=singlenode
#SBATCH --nodes=1
#SBATCH --ntasks=65
#SBATCH --time=1:00:00
#SBATCH --output=sbatch-output/output/bogner_exact_3box_s16_seed5_rho075_pinxy_omega1p9_50k_%j.out
#SBATCH --error=sbatch-output/error/bogner_exact_3box_s16_seed5_rho075_pinxy_omega1p9_50k_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# [Claude, 2026-09-24] Direct test of the wobble-instability hypothesis (THESIS_PLAN.md): same box
# scenario as job 4255656 (rho=0.75, seed5, omega=1.9), but with X/Y translation pinned (pinAxes
# <1,1,0>, Z free) and pitch/yaw rotation pinned (pinRotationAxes <1,0,1>, roll/Y free). If this
# restores clean settling near 26.565deg like the roll-only control (job 4253988) already did, that
# confirms the extra free DOFs (translation + pitch/yaw) -- not the hydrostatic force computation --
# are what drive the destabilisation documented throughout this session.
# 50000 steps (half of 4255656's 100k; the wobble-growth transition in the comparable Pass-3 run
# showed up by step ~40000, so 50k should be enough to see whether this constrained case avoids it).
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
echo "## Bogner exact box, s=16, rho=0.75, seed5, omega=1.9, X/Y+pitch/yaw pinned, roll+heave free, 50k"
echo "############################################################"
"$SH" --run --np 65 --prm "$PRM_DIR/bogner_exact_3box_s16_seed5_rho075_pinxy_omega1p9_50k_gas.prm"

STAMP="$(date +%Y%m%d_%H%M%S)"
VAULT_CANDIDATE="${HPCVAULT:-}/antidunes-archive/${SLURM_JOB_NAME}_${SLURM_JOB_ID}_${STAMP}"
if [[ -n "${HPCVAULT:-}" ]] && mkdir -p "$VAULT_CANDIDATE" 2>/dev/null \
   && touch "$VAULT_CANDIDATE/.write_test" 2>/dev/null; then
    rm -f "$VAULT_CANDIDATE/.write_test"
    ARCHIVE_ROOT="$VAULT_CANDIDATE"
else
    ARCHIVE_ROOT="$BIN_DIR/archive/${STAMP}_seed5_rho075_pinxy_omega1p9_50k"
    mkdir -p "$ARCHIVE_ROOT"
    echo "WARNING: \$HPCVAULT unavailable at archive time -- falling back to LOCAL archive at $ARCHIVE_ROOT (move this off manually to avoid quota pressure)" >&2
fi
for d in vtk-out eval snapshots; do
    [[ -d "$BIN_DIR/$d" ]] && mv "$BIN_DIR/$d" "$ARCHIVE_ROOT/$d"
done

echo "## DONE. Archived under $ARCHIVE_ROOT"
