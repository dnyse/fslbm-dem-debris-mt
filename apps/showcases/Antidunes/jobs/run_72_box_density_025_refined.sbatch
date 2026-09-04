#!/bin/bash -l
#SBATCH --job-name=antidunes-box-density-025-refined
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=08:00:00
#SBATCH --output=sbatch-output/output/antidunes_box_density_025_refined_%j.out
#SBATCH --error=sbatch-output/error/antidunes_box_density_025_refined_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# build + run: Fekken [2004] floating-cube equilibration, densityRatio 0.25, isolated
# single-box run (X,Y periodic, no side walls) -- REFINED-RESOLUTION / SHORTER-DOMAIN rerun
# of run_72_box_density_025.sbatch (see prms/box_density_025_refined_gas.prm for the full
# rationale): 1.5x finer resolution in X/Y (domainSize/cellsPerBlock scaled, box geometry
# scaled to match) plus a much shorter Z domain (liquidHeightFactor/floorHeightFactor cut
# down, domainSize.z 180->96) since this no-bed test never needed the original 124 LU of
# gas headspace. Net ~1.2x MORE total cells than the original (the X/Y refinement costs more
# than the Z cut saves) -- walltime bumped from 6h to 8h accordingly.
# domainSize <192, 90, 96> / cellsPerBlock <48, 30, 16> = 4x3x6 = 72 processes
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
    --prm "$REPO/apps/showcases/Antidunes/prms/box_density_025_refined_gas.prm"

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
