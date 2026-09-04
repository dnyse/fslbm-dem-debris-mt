#!/bin/bash -l
#SBATCH --job-name=antidunes-test-bridge-v1-lowwater-bigjam-cylinders-pureanalytical
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=04:00:00
#SBATCH --output=sbatch-output/output/antidunes_test_bridge_v1_lowwater_bigjam_cylinders_pureanalytical_%j.out
#SBATCH --error=sbatch-output/error/antidunes_test_bridge_v1_lowwater_bigjam_cylinders_pureanalytical_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# CYLINDER-FLOATER, PURE-ANALYTICAL variant of the 108-floater bigjam scenario: same 9-aim-point/
# 12-wave layout as job 4137147, but FloatingBox -> FloatingCylinder (radius 1.38, length 12,
# matched to the box's 2x3 cross-sectional area), 50000 timesteps (first-test smoke run, not a
# full jam-duration run), useAnalyticalContactDetection=true + usePureAnalyticalCylinder=true (no
# GJK/EPA fallback at all -- exercises the real Kodam eq. 82-101 edge-edge implementation and the
# rest of the analytical dispatch chain under real coupled-physics load). Run from a separate repo
# checkout (mt-dennys-mirror) so this doesn't collide with the main checkout's concurrent
# 162-floater cuboid run (job 4139454) or the no-analytical sibling variant on mt-dennys-mirror2.
# See prms/test_bridge_v1_lowwater_bigjam_cylinders_pureanalytical_gas.prm's own header for the
# full derivation.
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
    --prm "$REPO/apps/showcases/Antidunes/prms/test_bridge_v1_lowwater_bigjam_cylinders_pureanalytical_gas.prm"

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
