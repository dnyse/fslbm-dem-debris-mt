#!/bin/bash -l
#SBATCH --job-name=antidunes-bogner-exact-s16-densitycheck-rho0767-only
#SBATCH --partition=singlenode
#SBATCH --nodes=1
#SBATCH --ntasks=65
#SBATCH --time=1:15:00
#SBATCH --output=sbatch-output/output/bogner_exact_3box_s16_densitycheck_rho0767_only_%j.out
#SBATCH --error=sbatch-output/error/bogner_exact_3box_s16_densitycheck_rho0767_only_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# [Claude, 2026-09-24] REPEAT of just the rho=0.7665 (+2.2% density correction) sub-run from job 4260105,
# which crashed at step ~15323 (reynoldsNumber nan, a FormerMO reconstruction cascade -- 10 unknown cells,
# 90% with unknown neighbours, in 2-3 steps -- NOT the density-instrumented rho=0.75 sibling, which
# completed cleanly using the SAME binary in the SAME job). Verified before this rerun: the reconstruction
# algorithm itself (PdfReconstructionManager.h) has zero uncommitted changes, and both prior sub-runs used
# one binary built once for both -- so the crash is a genuine consequence of the density change (or the
# chaotic dynamics it put the system on), not a code difference. This rerun checks whether that crash was
# a one-off (chaotic, as with every other crash pattern seen this session) or reproduces.
#
# Same prm as before: bogner_exact_3box_s16_seed5_rho0767_nopin_omega1p9_100k_densitycheck_gas.prm
# (densityRatio 0.7665, density; already in the VTK writers list, dataSliceFrequency=500 already set).
# 100k steps, ~44 min expected (per job 4260105's own observed per-run cost) + ~5-10 min archive; 1h15m
# limit for margin.
# Queue: independent, ~/mt-dennys (confirmed clean/free before submission).
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
echo "## Bogner exact box equilibrium, s=16 density-check, densitycheck_rho0767_only"
echo "############################################################"
"$SH" --run --np 65 --prm "$PRM_DIR/bogner_exact_3box_s16_seed5_rho0767_nopin_omega1p9_100k_densitycheck_gas.prm"

STAMP="$(date +%Y%m%d_%H%M%S)"
VAULT_CANDIDATE="${HPCVAULT:-}/antidunes-archive/${SLURM_JOB_NAME}_${SLURM_JOB_ID}_${STAMP}"
if [[ -n "${HPCVAULT:-}" ]] && mkdir -p "$VAULT_CANDIDATE" 2>/dev/null \
   && touch "$VAULT_CANDIDATE/.write_test" 2>/dev/null; then
    rm -f "$VAULT_CANDIDATE/.write_test"
    ARCHIVE_ROOT="$VAULT_CANDIDATE"
else
    ARCHIVE_ROOT="$BIN_DIR/archive/${STAMP}_densitycheck_rho0767_only"
    mkdir -p "$ARCHIVE_ROOT"
    echo "WARNING: \$HPCVAULT unavailable at archive time -- falling back to LOCAL archive at $ARCHIVE_ROOT (move this off manually to avoid quota pressure)" >&2
fi
for d in vtk-out eval snapshots; do
    [[ -d "$BIN_DIR/$d" ]] && mv "$BIN_DIR/$d" "$ARCHIVE_ROOT/$d"
done

echo "## DONE. Archived under $ARCHIVE_ROOT"
