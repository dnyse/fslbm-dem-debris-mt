#!/bin/bash -l
#SBATCH --job-name=antidunes-bogner-exact-3box-s16-seed10-300k
#SBATCH --partition=singlenode
#SBATCH --nodes=1
#SBATCH --ntasks=65
#SBATCH --time=3:00:00
#SBATCH --output=sbatch-output/output/bogner_exact_3box_s16_seed10_300k_%j.out
#SBATCH --error=sbatch-output/error/bogner_exact_3box_s16_seed10_300k_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# [Claude, 2026-09-19] Thesis item 1a (THESIS_PLAN.md): LONG version (300000 steps) of the s=16 pinned, 10 deg-seed Bogner-exact run (job 4248454). There rho 0.5 reached 45.0
# and rho 0.25 ~34.1, but rho 0.75 was still creeping up (17.9@10k -> 21.9@50k, target 26.565): does it converge with more time? Same prm as 4248454 except timesteps
# 50000 -> 300000 (VTK every 6000). ~12-13 min per 50000 steps on 65 ranks => ~75 min; 3 h requested. Chained after the pins-off long run (same checkout,
# --dependency=afterany); the rescue-instead-of-delete startup logic below moves a predecessor's leftover output aside.
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

# [Claude, 2026-09-18] Never DELETE a previous run's output. These jobs are chained in ONE checkout with --dependency, so a
# predecessor that crashed or hit its time limit leaves vtk-out/eval/snapshots here -- the old `rm -rf` would have destroyed it.
# Move them aside instead (same filesystem -> instant rename); archive/rescued_* is then moved off by hand.
RESCUE="$BIN_DIR/archive/rescued_before_${SLURM_JOB_ID:-manual}_$(date +%Y%m%d_%H%M%S)"
for dir in vtk-out eval snapshots; do
    if [[ -d "$BIN_DIR/$dir" ]]; then
        mkdir -p "$RESCUE"
        echo "WARNING: leftover $BIN_DIR/$dir from a previous run -- moving to $RESCUE (NOT deleting)" >&2
        mv "$BIN_DIR/$dir" "$RESCUE/$dir"
    fi
done

"$SH" --build --jobs "$(nproc)"

echo "############################################################"
echo "## Bogner exact 3-box equilibrium test, s=16 seed10 300k steps, 65 ranks"
echo "############################################################"
# a crashed run must still reach the archive step below (set -e would skip it), so capture the exit status
set +e
"$SH" --run --np 65 --prm "$PRM_DIR/bogner_exact_3box_s16_seed10_300k_gas.prm"
RUN_RC=$?
set -e
echo "## run exited with status $RUN_RC -- archiving whatever output exists"

# Prefer $HPCVAULT, decided and verified HERE (right before use), not at job start -- see
# run_140_cuboid_stillwater_hires_d050.sbatch's own comment for the full rationale.
STAMP="$(date +%Y%m%d_%H%M%S)"
VAULT_CANDIDATE="${HPCVAULT:-}/antidunes-archive/${SLURM_JOB_NAME}_${SLURM_JOB_ID}_${STAMP}"
if [[ -n "${HPCVAULT:-}" ]] && mkdir -p "$VAULT_CANDIDATE" 2>/dev/null \
   && touch "$VAULT_CANDIDATE/.write_test" 2>/dev/null; then
    rm -f "$VAULT_CANDIDATE/.write_test"
    ARCHIVE_ROOT="$VAULT_CANDIDATE"
else
    ARCHIVE_ROOT="$BIN_DIR/archive/${STAMP}_bogner_exact_s16_seed10_300k"
    mkdir -p "$ARCHIVE_ROOT"
    echo "WARNING: \$HPCVAULT unavailable at archive time -- falling back to LOCAL archive at $ARCHIVE_ROOT (move this off manually to avoid quota pressure)" >&2
fi

for d in vtk-out eval snapshots; do
    [[ -d "$BIN_DIR/$d" ]] && mkdir -p "$ARCHIVE_ROOT/bogner_exact" && mv "$BIN_DIR/$d" "$ARCHIVE_ROOT/bogner_exact/$d"
done

echo "## DONE. Archived under $ARCHIVE_ROOT/bogner_exact/{eval,vtk-out,snapshots}"
exit "$RUN_RC"
