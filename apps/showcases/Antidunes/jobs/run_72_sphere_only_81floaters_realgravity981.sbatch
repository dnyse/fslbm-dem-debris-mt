#!/bin/bash -l
#SBATCH --job-name=antidunes-sphere-only-81floaters-realgravity981
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=04:00:00
#SBATCH --output=sbatch-output/output/antidunes_sphere_only_81floaters_realgravity981_%j.out
#SBATCH --error=sbatch-output/error/antidunes_sphere_only_81floaters_realgravity981_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# 81 spheres (all-sphere conversion of the proven-clean 81-floater "clean-
# factorial" grid, same 9 X-anchors x 3 Y-replicates x 3 Z-levels, same
# per-X-anchor density arrangement) dropped into still water, using TRUE
# Earth gravity g_SI=9.81 m/s^2 (gravityOverride=2.13715005840044e-05, same
# derivation as sphere_sedimentation_showcase_gas.prm's -- same domain). No
# relaxationRateOverride (standard useHardcodedUnits default, unlike the
# Bogner-gravity source file's tau=1/1.9). useBognerThreePhaseMEM true.
#
# Purpose: isolate whether the earlier 100-sphere crash (job 3996320, step
# ~9524/50000, NaN hydrodynamic force at a particle straddling the free
# surface under this same 9.81 gravity) was specific to that dense 100-sphere
# packing, or a more general 9.81-gravity + useBognerThreePhaseMEM
# interaction -- using the layout already proven collision-free (job
# 3856318) as the geometry instead. Complementary to the concurrently-running
# Bogner-gravity 81-floater job (3996932): same layout family, different
# gravity magnitude, to separate "layout" from "gravity" as variables.
#
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

"$REPO/apps/showcases/Antidunes/antidunes.sh" \
    --build --jobs "$(nproc)" \
    --run --np 72 \
    --prm "$REPO/apps/showcases/Antidunes/prms/sphere_only_81floaters_realgravity981_gas.prm"

# archive this run's own output so the next job doesn't overwrite it
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
