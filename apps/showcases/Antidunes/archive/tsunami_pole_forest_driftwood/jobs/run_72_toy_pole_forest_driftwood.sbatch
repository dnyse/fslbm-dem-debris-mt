#!/bin/bash -l
#SBATCH --job-name=antidunes-toy-pole-forest-driftwood
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=02:30:00
#SBATCH --output=sbatch-output/output/antidunes_toy_pole_forest_driftwood_%j.out
#SBATCH --error=sbatch-output/error/antidunes_toy_pole_forest_driftwood_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# NEW toy scenario: a staggered 8-pole "tsunami protection forest" hit by drifting cuboid debris,
# based on reference imagery (a laboratory tsunami/driftwood-vs-coastal-forest study). Water/flow
# setup reused verbatim from the validated toy-bridge family (steady PID-driven flow, NOT a tsunami
# wave/dry-bed runup -- see the prm's own header for why). Cuboids released only after an initial
# settling window so the flow has time to establish around the forest before debris arrives. 24
# cuboids, 4 aim points x 6 waves, 9000-step gap (the exact gap/box-size/velocity combo already
# validated this session in the toy trunk-first-staggered scenario). 80000 steps. See the prm's own
# header for the full derivation.
set -euo pipefail

if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    REPO="$(cd "$SLURM_SUBMIT_DIR/../../.." && pwd)"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    REPO="$(cd "$SCRIPT_DIR/../../../.." && pwd)"
fi
BIN_DIR="$REPO/build/apps/showcases/Antidunes"

"$REPO/apps/showcases/Antidunes/antidunes.sh" \
    --build --jobs "$(nproc)" \
    --run --np 72 \
    --prm "$REPO/apps/showcases/Antidunes/prms/toy_pole_forest_driftwood_gas.prm"

STAMP="$(date +%Y%m%d_%H%M%S)"
if [[ -n "$HPCVAULT" ]]; then
    ARCHIVE_DIR="$HPCVAULT/antidunes-archive/${SLURM_JOB_NAME}_$STAMP"
else
    ARCHIVE_DIR="$BIN_DIR/archive/$STAMP"
fi
mkdir -p "$ARCHIVE_DIR" || ARCHIVE_DIR="$BIN_DIR/archive/$STAMP"
[[ -d "$ARCHIVE_DIR" ]] || mkdir -p "$ARCHIVE_DIR"
for dir in vtk-out eval snapshots floater-history; do
    if [[ -d "$BIN_DIR/$dir" ]]; then
        mv "$BIN_DIR/$dir" "$ARCHIVE_DIR/$dir"
    fi
done
