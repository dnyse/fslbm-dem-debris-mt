#!/bin/bash -l
#SBATCH --job-name=antidunes-bogner-exact-s16-seed5-rho075-omega1p9
#SBATCH --partition=singlenode
#SBATCH --nodes=1
#SBATCH --ntasks=65
#SBATCH --time=0:35:00
#SBATCH --output=sbatch-output/output/bogner_exact_3box_s16_seed5_rho075_omega1p9_%j.out
#SBATCH --error=sbatch-output/error/bogner_exact_3box_s16_seed5_rho075_omega1p9_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# [Claude, 2026-09-22] Thesis item 1a (THESIS_PLAN.md): VISCOSITY/DAMPING TEST, rho=0.75 only, seeded at 5 deg (Bogner's own upright Archimedes draft z=16, NOT the equilibrium heave --
# a genuine free-fall-in test, unlike 4253988 which started AT the equilibrium). Motivation: every free run seeded at 2.86-10 deg overshoots past the 26.565 deg target and settles/wanders at 30-37 deg
# even after 300000 steps (not a runtime problem -- it genuinely lands in a different band). Hypothesis: insufficient rotational damping (fluid viscosity is the only physical damping available; the box
# has no DEM contacts) lets the swing overshoot past the true well. This run: relaxationRateOverride 1.9 (Bogner's own value, omega unchanged) -- the CONTROL of the 3-run sweep (1.5 and 1.0 in sibling
# jobs). 65 ranks, 50000 steps, ~11-13 min expected (same cost as the s=16 baselines). Queue: independent (no dependency), ~/mt-dennys-mirror2.
# See prms/bogner_exact_3box_s16_seed5_rho075_omega1p9_gas.prm and gen_bogner_exact_prm.py's new --omega option.
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
echo "## Bogner exact 3-box equilibrium test, s=16 seed5, 65 ranks"
echo "############################################################"
# a crashed run must still reach the archive step below (set -e would skip it), so capture the exit status
set +e
"$SH" --run --np 65 --prm "$PRM_DIR/bogner_exact_3box_s16_seed5_rho075_omega1p9_gas.prm"
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
    ARCHIVE_ROOT="$BIN_DIR/archive/${STAMP}_bogner_exact_s16_seed5"
    mkdir -p "$ARCHIVE_ROOT"
    echo "WARNING: \$HPCVAULT unavailable at archive time -- falling back to LOCAL archive at $ARCHIVE_ROOT (move this off manually to avoid quota pressure)" >&2
fi

for d in vtk-out eval snapshots; do
    [[ -d "$BIN_DIR/$d" ]] && mkdir -p "$ARCHIVE_ROOT/bogner_exact" && mv "$BIN_DIR/$d" "$ARCHIVE_ROOT/bogner_exact/$d"
done

echo "## DONE. Archived under $ARCHIVE_ROOT/bogner_exact/{eval,vtk-out,snapshots}"
exit "$RUN_RC"
