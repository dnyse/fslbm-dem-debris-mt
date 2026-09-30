#!/bin/bash -l
#SBATCH --job-name=antidunes-hopper-clog-cyl-dropwaves-capsule-optionC-only
#SBATCH --partition=multinode
#SBATCH --nodes=2
#SBATCH --ntasks=144
#SBATCH --time=01:00:00
#SBATCH --output=sbatch-output/output/hopper_clog_cyl_dropwaves_capsule_optionC_only_%j.out
#SBATCH --error=sbatch-output/error/hopper_clog_cyl_dropwaves_capsule_optionC_only_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# [Claude, 2026-09-25] NEW scenario per explicit request: same 210-rod half-density cylinder
# population as the halfrows hopper, but released as ~6 waves of 36 rods (6x2x3) dropped from ABOVE
# the water surface at the UPSTREAM end of the channel ("top of the stream"), instead of pre-placed
# in a block against the pillars -- see gen_hopper_clog_cyl_dropwaves_prm.py and
# prms/hopper_clog_cyl_dropwaves_capsule_gas.prm for the full derivation. THIS RUN IS DELIBERATELY WITH ALL
# OF THIS SESSION'S 2026-09-25 CRASH-MITIGATION FIXES DISABLED (Option A in
# StreamReconstructAdvectSweep.h, Option B/repairDegenerateObstacleCells in BuildTimeloop.h, the
# SimpleBB.h [MEM-GARBAGE-LINK] active guard) -- a genuine pre-session baseline for this new release
# mechanic, run on THIS mirror checkout ONLY (the local tree and other mirrors keep the fixes intact).
# See THESIS_PLAN.md for exactly which lines were reverted and why.
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
echo "## drop-waves half-density rod hopper, ALL 2026-09-25 crash-mitigation fixes DISABLED (baseline)"
echo "############################################################"
set +e
"$SH" --run --np 144 --prm "$PRM_DIR/hopper_clog_cyl_dropwaves_capsule_gas.prm"
RUN_RC=$?
set -e
echo "## run exited with status $RUN_RC -- archiving whatever output exists"

STAMP="$(date +%Y%m%d_%H%M%S)"
VAULT_CANDIDATE="${HPCVAULT:-}/antidunes-archive/${SLURM_JOB_NAME}_${SLURM_JOB_ID}_${STAMP}"
if [[ -n "${HPCVAULT:-}" ]] && mkdir -p "$VAULT_CANDIDATE" 2>/dev/null \
   && touch "$VAULT_CANDIDATE/.write_test" 2>/dev/null; then
    rm -f "$VAULT_CANDIDATE/.write_test"
    ARCHIVE_ROOT="$VAULT_CANDIDATE"
else
    ARCHIVE_ROOT="$BIN_DIR/archive/${STAMP}_hopper_clog_cyl_dropwaves_capsule_optionC_only"
    mkdir -p "$ARCHIVE_ROOT"
    echo "WARNING: \$HPCVAULT unavailable at archive time -- falling back to LOCAL archive at $ARCHIVE_ROOT (move this off manually to avoid quota pressure)" >&2
fi
for d in vtk-out eval snapshots; do
    [[ -d "$BIN_DIR/$d" ]] && mv "$BIN_DIR/$d" "$ARCHIVE_ROOT/$d"
done

echo "## DONE. Archived under $ARCHIVE_ROOT"
exit "$RUN_RC"
