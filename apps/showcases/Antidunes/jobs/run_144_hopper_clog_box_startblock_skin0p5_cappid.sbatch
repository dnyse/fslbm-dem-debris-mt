#!/bin/bash -l
#SBATCH --job-name=antidunes-hopper-clog-box-startblock-skin0p5-cappid
#SBATCH --partition=multinode
#SBATCH --nodes=2
#SBATCH --ntasks=144
#SBATCH --time=06:00:00
#SBATCH --output=sbatch-output/output/hopper_clog_box_startblock_skin0p5_cappid_%j.out
#SBATCH --error=sbatch-output/error/hopper_clog_box_startblock_skin0p5_cappid_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# [Claude, 2026-10-01] Skin-adapted full-scale cylinder hopper: same 406 rods (r=3, L=18, D:L=1:3,
# long axis along Y) pressed against the pillars, PID on, analytic cylinder contacts + owner-rank
# tiebreak as the original `hopper_clog_cyl_frontblock_pidon_2node_gas.prm`/job, with
# cylinderContactSkin 0.5 added (cd1b991b) and the row pitch widened 7->8 LU (via
# `gen_hopper_clog_cyl_prm.py --pitch 8.0 --skin 0.5`) so the skin's own 1.0 LU rod-rod contact
# threshold doesn't collide with the [RELEASE-OVERLAP] release-time guard the way it did at the
# 24-rod test scale (job 4282206, failed in 97s). That guard-vs-skin interaction was then confirmed
# fixed at 24 and 48 rods (jobs 4282267, 4282396: both full-length, exit 0, zero crash markers) --
# this is the same fix applied directly to the real ~400-rod production scenario, per explicit
# request ("adapt the original run 400 cylinder run to also work with this skin workaround").
# See hopper_clog_box_startblock_skin0p5_cappid_gas.prm's own generated header for the full
# geometry derivation. The un-skinned original sphere hopper took 1h48m for 60000 steps on 144
# ranks; rods add contact cost and very verbose [CylComputation] diagnostics (expect a
# multi-hundred-MB log), so 6h requested (same as the original cylinder job).
# [Claude, 2026-10-03] capped-PID copy of run_144_hopper_clog_cyl_frontblock_pidon_2node_skin0p5_curvguard.sbatch (406 rods, 100000 steps;
# 4284723 ran ~6.6 steps/s, so ~4.3 h plus VTK output)
# [Claude, 2026-10-04] DRIFTING CUBOIDS (gen_hopper_stress_prms.py --box-startblock-only): cuboid twin of the start-block
# rod run 4293495 (block x 8..232), 406 boxes 6x18x6, boxContactSkin 0.5, 150000 steps, release 500, VTK 0..149000,
# wallOuterGhostLayer, probes z 20, blockage ratio. 6 h requested (4293494 took 4h39).
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
# -> instant rename); archive/rescued_* is then moved off by hand.
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
echo "## clog measurement: hopper cylinders front-block, PID on, 2-node/144-rank, SKIN 0.5 adapted"
echo "############################################################"
set +e
"$SH" --run --np 144 --prm "$PRM_DIR/hopper_clog_box_startblock_skin0p5_cappid_gas.prm"
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
    ARCHIVE_ROOT="$BIN_DIR/archive/${STAMP}_hopper_clog_box_startblock_skin0p5_cappid"
    mkdir -p "$ARCHIVE_ROOT"
    echo "WARNING: \$HPCVAULT unavailable at archive time -- falling back to LOCAL archive at $ARCHIVE_ROOT (move this off manually to avoid quota pressure)" >&2
fi
for d in vtk-out eval snapshots; do
    [[ -d "$BIN_DIR/$d" ]] && mv "$BIN_DIR/$d" "$ARCHIVE_ROOT/$d"
done

echo "## DONE. Archived under $ARCHIVE_ROOT"
exit "$RUN_RC"
