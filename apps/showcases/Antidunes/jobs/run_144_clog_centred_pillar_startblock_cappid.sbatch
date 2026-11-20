#!/bin/bash -l
#SBATCH --job-name=antidunes-clog-centred-pillar-startblock-cappid
#SBATCH --partition=multinode
#SBATCH --nodes=2
#SBATCH --ntasks=144
#SBATCH --time=06:00:00
#SBATCH --output=sbatch-output/output/clog_centred_pillar_startblock_cappid_%j.out
#SBATCH --error=sbatch-output/error/clog_centred_pillar_startblock_cappid_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# Centred pillar, extreme scale, capped PID, start block (prms/clog_centred_pillar_startblock_cappid_gas.prm,
# from scripts/gen_hopper_stress_prms.py --centred-startblock-only): as run_144_clog_centred_pillar_cappid_current,
# but the 100 spheres R 5 start as one block at the upstream end and are all released at step 500, like the
# hopper start blocks. 150000 steps, 144 ranks / 2 nodes; 6 h requested (the staggered-release base scenario
# took ~2.5 h for 124750 steps).
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

for dir in vtk-out eval snapshots; do
    if [[ -d "$BIN_DIR/$dir" ]]; then
        echo "WARNING: removing leftover $BIN_DIR/$dir from an unarchived previous run" >&2
        rm -rf "${BIN_DIR:?}/$dir"
    fi
done

"$SH" --build --jobs "$(nproc)"

echo "############################################################"
echo "## clog measurement validation: extreme clog, PID on (fresh, 2-node/144-rank fallback, Z=180)"
echo "############################################################"
"$SH" --run --np 144 --prm "$PRM_DIR/clog_centred_pillar_startblock_cappid_gas.prm"

# Archive to $HPCVAULT if it is writable right now (checked here, not at job start: autofs mounts can go stale);
# the local archive under build/ is only a last-resort fallback.
STAMP="$(date +%Y%m%d_%H%M%S)"
VAULT_CANDIDATE="${HPCVAULT:-}/antidunes-archive/${SLURM_JOB_NAME}_${SLURM_JOB_ID}_${STAMP}"
if [[ -n "${HPCVAULT:-}" ]] && mkdir -p "$VAULT_CANDIDATE" 2>/dev/null \
   && touch "$VAULT_CANDIDATE/.write_test" 2>/dev/null; then
    rm -f "$VAULT_CANDIDATE/.write_test"
    ARCHIVE_ROOT="$VAULT_CANDIDATE"
else
    ARCHIVE_ROOT="$BIN_DIR/archive/${STAMP}_clog_centred_pillar_startblock_cappid"
    mkdir -p "$ARCHIVE_ROOT"
    echo "WARNING: \$HPCVAULT unavailable at archive time -- falling back to LOCAL archive at $ARCHIVE_ROOT (move this off manually to avoid quota pressure)" >&2
fi
for d in vtk-out eval snapshots; do
    [[ -d "$BIN_DIR/$d" ]] && mv "$BIN_DIR/$d" "$ARCHIVE_ROOT/$d"
done

echo "## DONE. Archived under $ARCHIVE_ROOT"
