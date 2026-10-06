#!/bin/bash -l
#SBATCH --job-name=antidunes-hopper-clog-box-frontblock-pidon-2node-skin0p5
#SBATCH --partition=multinode
#SBATCH --nodes=2
#SBATCH --ntasks=144
#SBATCH --time=06:00:00
#SBATCH --output=sbatch-output/output/hopper_clog_box_frontblock_pidon_2node_skin0p5_%j.out
#SBATCH --error=sbatch-output/error/hopper_clog_box_frontblock_pidon_2node_skin0p5_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# [Claude, 2026-10-01] Skin-adapted full-scale BOX hopper, run alongside the cylinder version of the
# same rerun (jobs/run_144_hopper_clog_cyl_frontblock_pidon_2node_skin0p5.sbatch), per explicit
# request ("on the other free mirror run the same for 400 cuboids, but also enable the skin").
# boxContactSkin 0.5 (new: box analog of cd1b991b's cylinderContactSkin, added this session --
# ParticleKernels.h/SimulationParams.h/SetupStages.h) + row pitch widened 7->8 LU via
# `gen_hopper_clog_box_prm.py --pitch 8.0 --skin 0.5` so the boxes have enough space between them to
# not trigger [RELEASE-OVERLAP] at release (same mechanism/fix as the cylinder version's job 4282206
# -> 4282267). Sealed-cell guard left off. 50000 timesteps (down from the original 60000), VTK on
# the discrete snapshot schedule 500/1000/2000/4000/8000/16000/32000/49000 (fluid_field_snapshots +
# the new addSnapshotAwareVTKWrite wiring for the per-box/pillar point data too -- BuildTimeloop.h).
# See hopper_clog_box_frontblock_pidon_2node_skin0p5_gas.prm's own generated header for the full
# geometry derivation. 6h requested (same as the original box job, which itself used well under this
# -- resize down with `scontrol update` once this run's own actual rate is known).
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
echo "## clog measurement: hopper boxes front-block, PID on, 2-node/144-rank, SKIN 0.5 adapted"
echo "############################################################"
set +e
"$SH" --run --np 144 --prm "$PRM_DIR/hopper_clog_box_frontblock_pidon_2node_skin0p5_gas.prm"
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
    ARCHIVE_ROOT="$BIN_DIR/archive/${STAMP}_hopper_clog_box_frontblock_pidon_2node_skin0p5"
    mkdir -p "$ARCHIVE_ROOT"
    echo "WARNING: \$HPCVAULT unavailable at archive time -- falling back to LOCAL archive at $ARCHIVE_ROOT (move this off manually to avoid quota pressure)" >&2
fi
for d in vtk-out eval snapshots; do
    [[ -d "$BIN_DIR/$d" ]] && mv "$BIN_DIR/$d" "$ARCHIVE_ROOT/$d"
done

echo "## DONE. Archived under $ARCHIVE_ROOT"
exit "$RUN_RC"
