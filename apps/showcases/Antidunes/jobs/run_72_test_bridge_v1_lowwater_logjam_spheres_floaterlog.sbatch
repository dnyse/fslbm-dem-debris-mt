#!/bin/bash -l
#SBATCH --job-name=antidunes-test-bridge-v1-lowwater-logjam-spheres-floaterlog
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=06:00:00
#SBATCH --output=sbatch-output/output/antidunes_test_bridge_v1_lowwater_logjam_spheres_floaterlog_%j.out
#SBATCH --error=sbatch-output/error/antidunes_test_bridge_v1_lowwater_logjam_spheres_floaterlog_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# Shape-control mirror of job 4134552 (test_bridge_v1_lowwater_logjam_floaterlog, FloatingBoxes):
# same geometry/timing/restitution/PID, but the 9 log-shaped FloatingBoxes are replaced 1:1 with
# radius-3 FloatingSpheres. That box run showed every log's rotation peaking at X~55, right where
# its 6-LU nose reaches the recessed arch collar's sharp X=61 edge -- logs consistently "stand up"
# there instead of floating through submerged. Spheres have no corner/edge to catch on that collar,
# so this isolates whether the standing-up behaviour is a box-shape artifact of the collar geometry
# or a more general property of the arch-tunnel contact model. See
# prms/test_bridge_v1_lowwater_logjam_spheres_floaterlog_gas.prm's own header for full rationale.
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
    --prm "$REPO/apps/showcases/Antidunes/prms/test_bridge_v1_lowwater_logjam_spheres_floaterlog_gas.prm"

STAMP="$(date +%Y%m%d_%H%M%S)"
if [[ -n "$HPCVAULT" ]]; then
    ARCHIVE_DIR="$HPCVAULT/antidunes-archive/${SLURM_JOB_NAME}_$STAMP"
else
    ARCHIVE_DIR="$BIN_DIR/archive/$STAMP"
fi
mkdir -p "$ARCHIVE_DIR"
for dir in vtk-out eval snapshots floater-history; do
    [[ -d "$BIN_DIR/$dir" ]] && mv "$BIN_DIR/$dir" "$ARCHIVE_DIR/$dir"
done
