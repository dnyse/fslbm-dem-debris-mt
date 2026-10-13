#!/bin/bash -l
#SBATCH --job-name=antidunes-skin-pocket-test
#SBATCH --partition=singlenode
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --time=01:30:00
#SBATCH --output=sbatch-output/output/skin_pocket_test_%j.out
#SBATCH --error=sbatch-output/error/skin_pocket_test_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# [Claude, 2026-10-07] Cheap test of the sub-cell pocket mechanism behind the skin-0.05 crash 4311067 (author: "create a simple
# and cheap test run proofing this. small domain only floaters"): prms/skin_pocket_test_m{0p05,0p1,0p25,0p5}_gas.prm
# (gen_hopper_stress_prms.py --skin-pocket-test-only; 84 rods of the 4293495 start block, 128 LU periodic channel, 8000 steps).
# The four skins run IN PARALLEL on one node, 16 ranks each, in this one checkout; their outputs cannot collide because every
# prm writes to its own folders (vtk-out_m<tag>, eval_m<tag>, snapshots_m<tag>) and its own log file. --bind-to none so that
# the four mpiruns do not all pin to cores 0-15.
set -uo pipefail

if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    REPO="$(cd "$SLURM_SUBMIT_DIR/../../.." && pwd)"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    REPO="$(cd "$SCRIPT_DIR/../../../.." && pwd)"
fi
SH="$REPO/apps/showcases/Antidunes/antidunes.sh"
PRM_DIR="$REPO/apps/showcases/Antidunes/prms"
BIN_DIR="$REPO/build/apps/showcases/Antidunes"
LOG_DIR="$REPO/apps/showcases/Antidunes/sbatch-output/output"
TAGS="0p05 0p1 0p25 0p5"

# never delete a previous run's output -- move aside
RESCUE="$BIN_DIR/archive/rescued_before_${SLURM_JOB_ID:-manual}_$(date +%Y%m%d_%H%M%S)"
for t in $TAGS; do for d in vtk-out_m$t eval_m$t snapshots_m$t; do
    if [[ -d "$BIN_DIR/$d" ]]; then mkdir -p "$RESCUE"; mv "$BIN_DIR/$d" "$RESCUE/$d"; fi
done; done

"$SH" --build --jobs "$(nproc)"
[[ -f "$BIN_DIR/Antidunes" ]] || { echo "binary missing"; exit 1; }

# [Claude, 2026-10-07] the runs call mpirun directly (four in parallel), so load the same modules antidunes.sh loads for itself
# (first attempt 4311135 failed with "mpirun: command not found")
module load gcc/11.5.0
module load openmpi/4.1.3-gcc11.5.0
cd "$BIN_DIR"
declare -A PIDS
for t in $TAGS; do
    mpirun -n 16 --bind-to none ./Antidunes "$PRM_DIR/skin_pocket_test_m${t}_gas.prm" > "$LOG_DIR/skin_pocket_test_m${t}_${SLURM_JOB_ID:-manual}.out" 2>&1 &
    PIDS[$t]=$!
done
RC=0
for t in $TAGS; do
    wait "${PIDS[$t]}"; rc=$?
    echo "## skin $t exited with status $rc"
    [[ $rc -ne 0 ]] && RC=1
done

STAMP="$(date +%Y%m%d_%H%M%S)"
VAULT_CANDIDATE="${HPCVAULT:-}/antidunes-archive/${SLURM_JOB_NAME:-skin-pocket-test}_${SLURM_JOB_ID:-manual}_${STAMP}"
if [[ -n "${HPCVAULT:-}" ]] && mkdir -p "$VAULT_CANDIDATE" 2>/dev/null; then ARCHIVE_ROOT="$VAULT_CANDIDATE"
else ARCHIVE_ROOT="$BIN_DIR/archive/${STAMP}_skin_pocket_test"; mkdir -p "$ARCHIVE_ROOT"; fi
for t in $TAGS; do for d in vtk-out_m$t eval_m$t snapshots_m$t; do
    [[ -d "$BIN_DIR/$d" ]] && mv "$BIN_DIR/$d" "$ARCHIVE_ROOT/$d"
done; cp "$LOG_DIR/skin_pocket_test_m${t}_${SLURM_JOB_ID:-manual}.out" "$ARCHIVE_ROOT/" 2>/dev/null; done
echo "## DONE. Archived under $ARCHIVE_ROOT"
exit "$RC"
