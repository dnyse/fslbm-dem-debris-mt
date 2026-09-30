#!/bin/bash -l
#SBATCH --job-name=antidunes-bogner-exact-s16-rho0767-tilt26p565-50k-pass3
#SBATCH --partition=singlenode
#SBATCH --nodes=1
#SBATCH --ntasks=65
#SBATCH --time=1:00:00
#SBATCH --output=sbatch-output/output/bogner_exact_3box_s16_rho0767_tilt26p565_50k_pass3_%j.out
#SBATCH --error=sbatch-output/error/bogner_exact_3box_s16_rho0767_tilt26p565_50k_pass3_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# [Claude, 2026-09-24] Fresh-from-t=0 rerun of the confirmed 4-for-4 step-444 crash scenario
# (bogner_exact_3box_s16_seed5_rho0767_nopin_omega1p9_50k_tilt26p565_gas.prm, jobs
# 4260975/4261005/4261093/4261316), with exactly one change: ModelParameters.usePass3Reresolution true
# (see the pass3 prm's own header and THESIS_PLAN.md). The 500-step resume-from-checkpoint-400
# reproduction (job 4262101) confirmed Pass-3 avoids the known step-444 crash by correctly reclassifying
# the formerMO cell to Gas instead of Interface. This run checks the FULL fresh 50000-step trajectory --
# does it run clean throughout, crash somewhere else, and where does the box actually settle.
# Timing: same domain/box/steps as the original 50k run (~31min observed) + archive overhead; 1h limit.
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
echo "## Bogner exact box equilibrium, s=16, rho=0.7665, tilt=26.565deg, 50k, Pass-3 enabled"
echo "############################################################"
"$SH" --run --np 65 --prm "$PRM_DIR/bogner_exact_3box_s16_seed5_rho0767_nopin_omega1p9_50k_tilt26p565_pass3_gas.prm"

STAMP="$(date +%Y%m%d_%H%M%S)"
VAULT_CANDIDATE="${HPCVAULT:-}/antidunes-archive/${SLURM_JOB_NAME}_${SLURM_JOB_ID}_${STAMP}"
if [[ -n "${HPCVAULT:-}" ]] && mkdir -p "$VAULT_CANDIDATE" 2>/dev/null \
   && touch "$VAULT_CANDIDATE/.write_test" 2>/dev/null; then
    rm -f "$VAULT_CANDIDATE/.write_test"
    ARCHIVE_ROOT="$VAULT_CANDIDATE"
else
    ARCHIVE_ROOT="$BIN_DIR/archive/${STAMP}_rho0767_tilt26p565_50k_pass3"
    mkdir -p "$ARCHIVE_ROOT"
    echo "WARNING: \$HPCVAULT unavailable at archive time -- falling back to LOCAL archive at $ARCHIVE_ROOT (move this off manually to avoid quota pressure)" >&2
fi
for d in vtk-out eval snapshots; do
    [[ -d "$BIN_DIR/$d" ]] && mv "$BIN_DIR/$d" "$ARCHIVE_ROOT/$d"
done

echo "## DONE. Archived under $ARCHIVE_ROOT"
