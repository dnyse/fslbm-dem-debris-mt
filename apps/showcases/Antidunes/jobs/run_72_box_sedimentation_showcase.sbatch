#!/bin/bash -l
#SBATCH --job-name=antidunes-box-sedimentation-showcase
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=03:00:00
#SBATCH --output=sbatch-output/output/antidunes_box_sedimentation_showcase_%j.out
#SBATCH --error=sbatch-output/error/antidunes_box_sedimentation_showcase_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# build + run: CUBOID variant of the sedimentation showcase, v6. v1 (100-cube
# 5x4x5 grid) crashed -- a WALBERLA_CHECK "entry already present" in mesa_pd's
# ContactHistoryNotification.h (the contact-ownership race, traced to ~90 of
# 100 cubes packed near rank boundaries simultaneously). v2 cut this to 9
# cubes; v3 briefly tried thin sticks; v4 introduced Y-neighbour contacts for
# the first time (18-box 3(X)x2(Y)x3(Z) grid, halfExtents <3.5,7,3.5>, only
# 2 LU Y-clearance) and crashed differently -- a NaN rotation quaternion
# (core/math/Quaternion.h) at 174 sec. v5 widened the Y-clearance to 10 LU on
# the theory that the tight clearance itself was the cause -- but v5 crashed
# AGAIN with the identical signature at almost the same time (179 sec), now on
# two ADJACENT ranks simultaneously, which points instead at the same
# still-unpatched ContactFilter race from v1 (a double-applied force large
# enough to blow up angular velocity to NaN within one DEM sub-cycle, before
# the once-per-step duplicate-history check ever runs).
#
# This version (v6, see prms/box_sedimentation_showcase_gas.prm for the full
# rationale) is a controlled test: it returns to EXACTLY v2's layout (single Y
# position per column, 9 boxes, no Y-neighbour contacts at all -- only
# gravity-driven Z-stacking, which never crashed in this file's history) while
# keeping v4/v5's box shape. If this runs cleanly, it isolates the Y-neighbour
# contact type as necessary for the crash; if it crashes too, that points
# squarely at the underlying mesa_pd race regardless of further showcase-level
# tuning.
# domainSize <128, 60, 180> / cellsPerBlock <32, 20, 30> = 4x3x6 = 72 processes
# (unchanged from every other showcase in this project)
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
    --prm "$REPO/apps/showcases/Antidunes/prms/box_sedimentation_showcase_gas.prm"

# archive this run's own output so the next job doesn't overwrite it
STAMP="$(date +%Y%m%d_%H%M%S)"
# if $HPCVAULT is set, archive to persistent vault storage; otherwise the local build dir
if [[ -n "$HPCVAULT" ]]; then
    ARCHIVE_DIR="$HPCVAULT/antidunes-archive/${SLURM_JOB_NAME}_$STAMP"
else
    ARCHIVE_DIR="$BIN_DIR/archive/$STAMP"
fi
mkdir -p "$ARCHIVE_DIR"
for dir in vtk-out eval snapshots; do
    [[ -d "$BIN_DIR/$dir" ]] && mv "$BIN_DIR/$dir" "$ARCHIVE_DIR/$dir"
done
