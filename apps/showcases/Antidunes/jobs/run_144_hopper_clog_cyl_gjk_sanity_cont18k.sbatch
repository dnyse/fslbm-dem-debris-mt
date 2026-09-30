#!/bin/bash -l
#SBATCH --job-name=antidunes-hopper-clog-cyl-gjk-sanity-cont18k
#SBATCH --partition=multinode
#SBATCH --nodes=2
#SBATCH --ntasks=144
#SBATCH --time=01:30:00
#SBATCH --output=sbatch-output/output/hopper_clog_cyl_gjk_sanity_cont18k_%j.out
#SBATCH --error=sbatch-output/error/hopper_clog_cyl_gjk_sanity_cont18k_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# [Claude, 2026-09-20] CONTINUATION of 4252146 from its step-15001 checkpoint, +18000 steps, VTK every 1000 (1.5 h limit); see prms/hopper_clog_cyl_frontblock_pidon_2node_gjk_sanity_18k_cont18k_gas.prm's header. Output (vtk-out, eval, floater-history) archived to the vault.
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
for dir in vtk-out eval snapshots floater-history; do
    if [[ -d "$BIN_DIR/$dir" ]]; then
        mkdir -p "$RESCUE"
        echo "WARNING: leftover $BIN_DIR/$dir from a previous run -- moving to $RESCUE (NOT deleting)" >&2
        mv "$BIN_DIR/$dir" "$RESCUE/$dir"
    fi
done

CKPT="${HPCVAULT:-/nonexistent}/antidunes-archive/antidunes-hopper-clog-cyl-gjk-sanity-18k_4252146_20260920_215229/snapshots"
if [[ ! -f "$CKPT/simState.file" ]]; then echo "FATAL: checkpoint $CKPT/simState.file not found" >&2; exit 1; fi
mkdir -p "$BIN_DIR/snapshots"
cp -a "$CKPT/." "$BIN_DIR/snapshots/"
echo "## restored checkpoint: simState = $(head -1 "$BIN_DIR/snapshots/simState.file")  (expect 15001)"
du -sh "$BIN_DIR/snapshots"

"$SH" --build --jobs "$(nproc)"

echo "############################################################"
echo "## hopper_clog_cyl_gjk_sanity_cont18k"
echo "############################################################"
# a crashed run must still reach the archive step below (set -e would skip it), so capture the exit status
set +e
"$SH" --run --np 144 --prm "$PRM_DIR/hopper_clog_cyl_frontblock_pidon_2node_gjk_sanity_18k_cont18k_gas.prm"
RUN_RC=$?
set -e
echo "## run exited with status $RUN_RC -- archiving whatever output exists"

# Prefer $HPCVAULT, decided and verified HERE (right before use) -- see the 288-rank sibling job's
# own 2026-09-16 comment for the full rationale (autofs staleness, local-only is a last resort).
STAMP="$(date +%Y%m%d_%H%M%S)"
VAULT_CANDIDATE="${HPCVAULT:-}/antidunes-archive/${SLURM_JOB_NAME}_${SLURM_JOB_ID}_${STAMP}"
if [[ -n "${HPCVAULT:-}" ]] && mkdir -p "$VAULT_CANDIDATE" 2>/dev/null \
   && touch "$VAULT_CANDIDATE/.write_test" 2>/dev/null; then
    rm -f "$VAULT_CANDIDATE/.write_test"
    ARCHIVE_ROOT="$VAULT_CANDIDATE"
else
    ARCHIVE_ROOT="$BIN_DIR/archive/${STAMP}_hopper_clog_cyl_gjk_sanity_cont18k"
    mkdir -p "$ARCHIVE_ROOT"
    echo "WARNING: \$HPCVAULT unavailable at archive time -- falling back to LOCAL archive at $ARCHIVE_ROOT (move this off manually to avoid quota pressure)" >&2
fi
for d in vtk-out eval snapshots floater-history; do
    [[ -d "$BIN_DIR/$d" ]] && mv "$BIN_DIR/$d" "$ARCHIVE_ROOT/$d"
done

echo "## DONE. Archived under $ARCHIVE_ROOT"
exit "$RUN_RC"
