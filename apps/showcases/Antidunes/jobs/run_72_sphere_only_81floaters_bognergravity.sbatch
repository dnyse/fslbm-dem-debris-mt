#!/bin/bash -l
#SBATCH --job-name=antidunes-sphere-only-81floaters-bognergravity
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=04:00:00
#SBATCH --output=sbatch-output/output/antidunes_sphere_only_81floaters_bognergravity_%j.out
#SBATCH --error=sbatch-output/error/antidunes_sphere_only_81floaters_bognergravity_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# 81 spheres (same all-sphere conversion of the proven-clean 81-floater grid
# used by sphere_only_81floaters_realgravity981_gas.prm), but with Bogner &
# Rude's own gravity+viscosity pair (gravityOverride=1.66667e-4,
# relaxationRateOverride=1.9) -- identical physics to the concurrently-running
# mixed_shape_stillwater_81floaters_bogner_gas.prm (job 3996932). Isolates
# shape (spheres vs. mixed) as the variable, holding gravity/viscosity fixed
# at the combination already shown to run stably.
#
# domainSize <128, 60, 180> / cellsPerBlock <32, 20, 30> = 4x3x6 = 72 blocks
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
    --prm "$REPO/apps/showcases/Antidunes/prms/sphere_only_81floaters_bognergravity_gas.prm"

STAMP="$(date +%Y%m%d_%H%M%S)"
if [[ -n "$HPCVAULT" ]]; then
    ARCHIVE_DIR="$HPCVAULT/antidunes-archive/${SLURM_JOB_NAME}_$STAMP"
else
    ARCHIVE_DIR="$BIN_DIR/archive/$STAMP"
fi
mkdir -p "$ARCHIVE_DIR"
for dir in vtk-out eval snapshots; do
    [[ -d "$BIN_DIR/$dir" ]] && mv "$BIN_DIR/$dir" "$ARCHIVE_DIR/$dir"
done
