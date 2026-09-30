#!/bin/bash -l
#SBATCH --job-name=antidunes-hopper-clog-cyl-halfrows-gjk-pass3-precursor
#SBATCH --partition=multinode
#SBATCH --nodes=2
#SBATCH --ntasks=144
#SBATCH --time=00:30:00
#SBATCH --output=sbatch-output/output/hopper_clog_cyl_halfrows_gjk_pass3_precursor_%j.out
#SBATCH --error=sbatch-output/error/hopper_clog_cyl_halfrows_gjk_pass3_precursor_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# [Claude, 2026-09-24] Out-of-interest variant of job 4255679 (the half-density cylinder hopper, no
# guards/caps at all: no sealed-cell guard, no rodCapOverlap, no fluid sanity sweep, no dense NaN scan)
# with exactly ONE change: ModelParameters.usePass3Reresolution true. 4255679's own crash was already
# root-caused (THESIS_PLAN.md, 2026-09-23): a MEM/DEM pinched-cell density blowup between two fully
# submerged rods (STALE-PARTICLE-IDX at step 3133, preceded by a negative-density [MEM-GARBAGE-LINK]
# cell with only 9-15/26 obstacle neighbours -- below SealedCellGuard.h's own guard24 threshold, and
# guard24's own rerun (job 4259313) failed anyway via a different signature). Pass-3 re-resolution only
# affects classifyFormerObstacleNeighbourhood's free-surface (gas/liquid/interface) classification for
# cells becoming fluid near the waterline -- a different code path from ordinary submerged-liquid MEM
# reconstruction between two rods well below the surface -- so it is NOT expected to fix this specific
# crash mechanism. Run anyway, out of interest, to confirm directly rather than assume.
# Same half-density population as 4255679 (210 rods, downstream half, `gen_hopper_clog_cyl_prm.py
# --rows 15`), 144 ranks, 2h30m limit (same budget as the original).
# See prms/hopper_clog_cyl_frontblock_pidon_2node_halfrows_pass3_precursor_gas.prm.
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
echo "## clog measurement: half-density rod hopper WITH pillars, PID on, NO guards, Pass-3 ON, denseNanScan armed 3000-3200, capped 5000 steps, 2-node/144-rank"
echo "############################################################"
# a crashed run must still reach the archive step below (set -e would skip it), so capture the exit status
set +e
"$SH" --run --np 144 --prm "$PRM_DIR/hopper_clog_cyl_frontblock_pidon_2node_halfrows_pass3_precursor_gas.prm"
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
    ARCHIVE_ROOT="$BIN_DIR/archive/${STAMP}_hopper_clog_cyl_halfrows_gjk_pass3_precursor"
    mkdir -p "$ARCHIVE_ROOT"
    echo "WARNING: \$HPCVAULT unavailable at archive time -- falling back to LOCAL archive at $ARCHIVE_ROOT (move this off manually to avoid quota pressure)" >&2
fi
for d in vtk-out eval snapshots; do
    [[ -d "$BIN_DIR/$d" ]] && mv "$BIN_DIR/$d" "$ARCHIVE_ROOT/$d"
done

echo "## DONE. Archived under $ARCHIVE_ROOT"
exit "$RUN_RC"
