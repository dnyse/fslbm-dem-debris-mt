#!/bin/bash -l
#SBATCH --job-name=antidunes-hopper-clog-cyl-restart25001-densevtk
#SBATCH --partition=multinode
#SBATCH --nodes=2
#SBATCH --ntasks=144
#SBATCH --time=01:00:00
#SBATCH --output=sbatch-output/output/hopper_clog_cyl_restart25001_densevtk_%j.out
#SBATCH --error=sbatch-output/error/hopper_clog_cyl_restart25001_densevtk_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# [Claude, 2026-10-02] DENSE-VTK RESTART of the crashed job 4284723 (406-rod hopper,
# cylinderContactSkin 0.5 + guardCurvatureInputs true), from its own checkpoint at step 25001
# (copied -- not moved -- from the vault archive's snapshots/ into
# build/apps/showcases/Antidunes/snapshots_restart_4284723/ ahead of this job). 4284723 crashed at
# step 27344 ([STALE-PARTICLE-IDX], a downstream symptom of 80 [MEM-GARBAGE-LINK] warnings starting
# at step 27338 and the first NaN normal at 27339). This restart adds diagnoseNaNNormals + the full
# precursor-scan suite windowed around 27200-27400, and a DENSE VTK snapshot schedule (every single
# step from 27334 through 27344) cropped to the hopper + full downstream channel, to see exactly
# where cylinder uid 2666130979403333644 (riding the drained downstream sheet since step ~23000)
# and the MEM-GARBAGE-LINK cell interact. 2400 steps (25001->27400) at ~6 steps/s plus VTK I/O, 1h
# cap. Deliberately does NOT rebuild -- same checkout/binary that produced 4284723's own crash, per
# explicit instruction ("do NOT resync or change code first: the restart should use the code that
# crashed").
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

# [Claude, 2026-09-18] Never DELETE a previous run's output -- move aside instead (same filesystem
# -> instant rename); archive/rescued_* is then moved off by hand. Only touches vtk-out/eval/
# snapshots (NOT snapshots_restart_4284723, the checkpoint copy this job reads from).
RESCUE="$BIN_DIR/archive/rescued_before_${SLURM_JOB_ID:-manual}_$(date +%Y%m%d_%H%M%S)"
for dir in vtk-out eval snapshots; do
    if [[ -d "$BIN_DIR/$dir" ]]; then
        mkdir -p "$RESCUE"
        echo "WARNING: leftover $BIN_DIR/$dir from a previous run -- moving to $RESCUE (NOT deleting)" >&2
        mv "$BIN_DIR/$dir" "$RESCUE/$dir"
    fi
done

if [[ ! -d "$BIN_DIR/snapshots_restart_4284723" ]]; then
    echo "FATAL: $BIN_DIR/snapshots_restart_4284723 is missing -- this job reads the checkpoint from there (step 25001, copied from 4284723's own vault archive). Aborting before even trying to run." >&2
    exit 1
fi

# [Claude, 2026-10-02] Deliberately NOT building: this restart must use the exact binary that
# produced 4284723's own crash, not a freshly (re)compiled one.
echo "## SKIPPING build -- using the existing ./Antidunes binary, unchanged since job 4284723"

echo "############################################################"
echo "## DENSE-VTK RESTART of crashed job 4284723, from checkpoint step 25001, to step 27400"
echo "############################################################"
set +e
"$SH" --run --np 144 --prm "$PRM_DIR/hopper_clog_cyl_frontblock_pidon_2node_skin0p5_curvguard_restart25001_densevtk_gas.prm"
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
    ARCHIVE_ROOT="$BIN_DIR/archive/${STAMP}_hopper_clog_cyl_restart25001_densevtk"
    mkdir -p "$ARCHIVE_ROOT"
    echo "WARNING: \$HPCVAULT unavailable at archive time -- falling back to LOCAL archive at $ARCHIVE_ROOT (move this off manually to avoid quota pressure)" >&2
fi
for d in vtk-out eval snapshots; do
    [[ -d "$BIN_DIR/$d" ]] && mv "$BIN_DIR/$d" "$ARCHIVE_ROOT/$d"
done

# [Claude, 2026-10-02] Note the checkpoint's own provenance in the archive instead of duplicating
# the 2.4G of snapshot files a second time -- it's already sitting in 4284723's own vault archive.
cat > "$ARCHIVE_ROOT/PROVENANCE_restart.txt" <<EOF
This run (job ${SLURM_JOB_ID}, ${SLURM_JOB_NAME}) restarted from job 4284723's own checkpoint at
step 25001, copied (not moved) from:
  /home/vault/iwia/iwia148h/antidunes-archive/antidunes-hopper-clog-cyl-frontblock-pidon-2node-skin0p5-curvguard_4284723_20261002_154014/snapshots/
into build/apps/showcases/Antidunes/snapshots_restart_4284723/ ahead of this job.
That checkpoint is NOT duplicated here -- it already lives in 4284723's own vault archive above.
Same checkout/binary as 4284723 (no resync, no rebuild before this run).
EOF

echo "## DONE. Archived under $ARCHIVE_ROOT"
exit "$RUN_RC"
