#!/bin/bash -l
#SBATCH --job-name=antidunes-bogner-exact-s16-rho0767-tilt26-zempirical-50k
#SBATCH --partition=singlenode
#SBATCH --nodes=1
#SBATCH --ntasks=65
#SBATCH --time=1:00:00
#SBATCH --output=sbatch-output/output/bogner_exact_3box_s16_rho0767_tilt26_zempirical_50k_%j.out
#SBATCH --error=sbatch-output/error/bogner_exact_3box_s16_rho0767_tilt26_zempirical_50k_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# [Claude, 2026-09-24] Starting-angle sensitivity test at the "corrected" density 0.7665, omega=1.9
# (see THESIS_PLAN.md): job 4260604 (seeded at 5deg) tracked 26.565deg closely for ~15,000 steps then
# grew into a wide chaotic orbit. This is the 26.565deg-seed sibling (released AT the theoretical equilibrium, zero seed transient), run
# alongside the 2.9deg and 10deg siblings on the other two checkouts, all three simultaneously.
# 50000 steps only (half of 4260604's length); ~31 min expected (half of 4260604's 62 min) + archive
# overhead; 1h limit for margin.
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
echo "## Bogner exact box equilibrium, s=16, rho=0.7665, tilt=2.9deg, 50k"
echo "############################################################"
"$SH" --run --np 65 --prm "$PRM_DIR/bogner_exact_3box_s16_seed5_rho0767_nopin_omega1p9_50k_tilt26_zempirical_gas.prm"

STAMP="$(date +%Y%m%d_%H%M%S)"
VAULT_CANDIDATE="${HPCVAULT:-}/antidunes-archive/${SLURM_JOB_NAME}_${SLURM_JOB_ID}_${STAMP}"
if [[ -n "${HPCVAULT:-}" ]] && mkdir -p "$VAULT_CANDIDATE" 2>/dev/null \
   && touch "$VAULT_CANDIDATE/.write_test" 2>/dev/null; then
    rm -f "$VAULT_CANDIDATE/.write_test"
    ARCHIVE_ROOT="$VAULT_CANDIDATE"
else
    ARCHIVE_ROOT="$BIN_DIR/archive/${STAMP}_rho0767_tilt26_zempirical_50k"
    mkdir -p "$ARCHIVE_ROOT"
    echo "WARNING: \$HPCVAULT unavailable at archive time -- falling back to LOCAL archive at $ARCHIVE_ROOT (move this off manually to avoid quota pressure)" >&2
fi
for d in vtk-out eval snapshots; do
    [[ -d "$BIN_DIR/$d" ]] && mv "$BIN_DIR/$d" "$ARCHIVE_ROOT/$d"
done

echo "## DONE. Archived under $ARCHIVE_ROOT"
