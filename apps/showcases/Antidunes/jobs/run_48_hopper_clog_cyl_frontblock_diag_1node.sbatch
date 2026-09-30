#!/bin/bash -l
#SBATCH --job-name=antidunes-hopper-clog-cyl-frontblock-diag-1node
#SBATCH --partition=singlenode
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --time=01:00:00
#SBATCH --output=sbatch-output/output/hopper_clog_cyl_frontblock_diag_1node_%j.out
#SBATCH --error=sbatch-output/error/hopper_clog_cyl_frontblock_diag_1node_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# [Claude, 2026-09-20] DIAGNOSTIC for the cylinder-hopper abort (job 4248725: [RELEASE-OVERLAP] at mainStep 2001, subCycle 2, between rods that were 1 and 8 LU apart and at rest at sub-cycle 0).
# The real prm on 48 ranks / one node (cellsPerBlock <128,20,45>), stopped at step 2400 (400 steps past the release at 2000). ParticleKernels.h now aborts only at sub-cycle 0 (a real insertion
# overlap); a contact at a LATER sub-cycle of the first step after a release is logged as [RELEASE-CONTACT-LATER] with position, quaternion, velocities, force, torque and hydro force/torque of BOTH bodies.
# Look for those lines: they should show which rod moved/rotated (or which copy has stale state) between sub-cycle 0 and 2. Expected ~12 min; 1 h requested. Little output (VTK only at step 0).
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
echo "## clog measurement validation: hopper cylinders front-block, PID on, 2-node/144-rank"
echo "############################################################"
# a crashed run must still reach the archive step below (set -e would skip it), so capture the exit status
set +e
"$SH" --run --np 48 --prm "$PRM_DIR/hopper_clog_cyl_frontblock_diag_1node_gas.prm"
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
    ARCHIVE_ROOT="$BIN_DIR/archive/${STAMP}_hopper_clog_cyl_frontblock_diag_1node"
    mkdir -p "$ARCHIVE_ROOT"
    echo "WARNING: \$HPCVAULT unavailable at archive time -- falling back to LOCAL archive at $ARCHIVE_ROOT (move this off manually to avoid quota pressure)" >&2
fi
for d in vtk-out eval snapshots; do
    [[ -d "$BIN_DIR/$d" ]] && mv "$BIN_DIR/$d" "$ARCHIVE_ROOT/$d"
done

echo "## DONE. Archived under $ARCHIVE_ROOT"
exit "$RUN_RC"
