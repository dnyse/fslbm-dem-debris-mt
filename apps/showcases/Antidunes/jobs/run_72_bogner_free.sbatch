#!/bin/bash -l
#SBATCH --job-name=antidunes-bogner-free
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=06:00:00
#SBATCH --output=sbatch-output/output/antidunes_bogner_free_%j.out
#SBATCH --error=sbatch-output/error/antidunes_bogner_free_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# BOGNER 6:4 cuboid — EQUILIBRIUM POSITION (coarse resolution only). Free-floating 6:4 box (densityRatio
# 0.5), started 10° off upright, in still water. A wide>tall box is stable UPRIGHT, so it should right
# itself: tiltY → ~0°, ~50% submersion, with the righting moment M_y → 0 at equilibrium. Analytic buoyancy
# ON (exact submerged centroid drives the righting torque).
# grep -E 'BoxBuoyancy|BoxRightingMoment' on the .out; watch tiltY and M_y settle.
# domainSize <128, 60, 180> / cellsPerBlock <32, 20, 30> = 4x3x6 = 72 blocks
# Resolve the repo root programmatically (mirrors antidunes.sh's own SCRIPT_DIR/REPO_ROOT
# derivation) -- no hardcoded path, so this works across mirror clones (mt-dennys,
# mt-dennys-mirror, mt-dennys-mirror2, ...) regardless of which one submitted the job. Jobs
# are always submitted with `sbatch` from apps/showcases/Antidunes (same dir as antidunes.sh),
# so SLURM_SUBMIT_DIR (guaranteed set by SLURM to that submission directory) sits 3 levels
# below the repo root; fall back to deriving it from this script's own location (four levels
# up from apps/showcases/Antidunes/jobs/) for direct/non-sbatch execution.
if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    REPO="$(cd "$SLURM_SUBMIT_DIR/../../.." && pwd)"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    REPO="$(cd "$SCRIPT_DIR/../../../.." && pwd)"
fi
BIN_DIR="$REPO/build/apps/showcases/Antidunes"
SH="$REPO/apps/showcases/Antidunes/antidunes.sh"
PRM_DIR="$REPO/apps/showcases/Antidunes/prms"
STAMP="$(date +%Y%m%d_%H%M%S)"

"$SH" --build --jobs "$(nproc)"

echo "############################################################"
echo "## Bogner 6:4 FREE — equilibrium position"
echo "############################################################"
"$SH" --run --np 72 --prm "$PRM_DIR/bogner_righting_6to4_free_gas.prm"
for d in vtk-out eval snapshots; do
    [[ -d "$BIN_DIR/$d" ]] && mkdir -p "$BIN_DIR/archive/$STAMP" && mv "$BIN_DIR/$d" "$BIN_DIR/archive/$STAMP/free_$d"
done

echo "## DONE. grep -E 'BoxBuoyancy|BoxRightingMoment' on this .out; tiltY→~0°, M_y→0 at equilibrium."
