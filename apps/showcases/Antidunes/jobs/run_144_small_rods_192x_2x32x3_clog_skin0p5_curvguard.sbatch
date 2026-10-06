#!/bin/bash -l
#SBATCH --job-name=antidunes-small-rods-192x-2x32x3-clog-skin0p5-curvguard
#SBATCH --partition=multinode
#SBATCH --nodes=2
#SBATCH --ntasks=144
#SBATCH --time=02:00:00
#SBATCH --output=sbatch-output/output/small_rods_192x_2x32x3_clog_skin0p5_curvguard_%j.out
#SBATCH --error=sbatch-output/error/small_rods_192x_2x32x3_clog_skin0p5_curvguard_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# [Claude, 2026-10-01] Doubling-in-X step after the clean 48-rod skin test (job 4282396, full
# 51000 steps, exit 0, zero crash markers -- same result as the 24-rod test 4282267). Same
# cylinderContactSkin 0.5, sealed-cell guard still off, same widened 8 LU row pitch -- now
# 32 rows instead of 8 (2x32x3 = 192 rods), per explicit request ("we know the 48
# cylinder works, then on two new freed directory run the 96 and the 192 version with the skin
# activated"). 2h cap (4282396 ran 51000 steps well under 1h on 144 ranks; this grows the
# floater count but not the domain/fluid work, so should stay well inside 2h).
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
echo "## one hundred ninety-two small rods (2x32x3) hopper, cylinderContactSkin 0.5, sealed-cell guard OFF"
echo "## (doubling-in-X step after the clean 48-rod skin test, job 4282396)"
echo "############################################################"
set +e
"$SH" --run --np 144 --prm "$PRM_DIR/small_rods_192x_2x32x3_clog_skin0p5_curvguard_gas.prm"
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
    ARCHIVE_ROOT="$BIN_DIR/archive/${STAMP}_small_rods_192x_2x32x3_clog_skin0p5_curvguard"
    mkdir -p "$ARCHIVE_ROOT"
    echo "WARNING: \$HPCVAULT unavailable at archive time -- falling back to LOCAL archive at $ARCHIVE_ROOT (move this off manually to avoid quota pressure)" >&2
fi
for d in vtk-out eval snapshots; do
    [[ -d "$BIN_DIR/$d" ]] && mv "$BIN_DIR/$d" "$ARCHIVE_ROOT/$d"
done

echo "## DONE. Archived under $ARCHIVE_ROOT"
exit "$RUN_RC"
