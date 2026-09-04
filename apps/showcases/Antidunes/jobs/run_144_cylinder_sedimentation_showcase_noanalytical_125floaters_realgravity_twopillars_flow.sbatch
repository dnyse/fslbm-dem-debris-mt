#!/bin/bash -l
#SBATCH --job-name=antidunes-cylinder-sedimentation-showcase-noanalytical-125floaters-realgravity-twopillars-flow
#SBATCH --partition=multinode
#SBATCH --constraint=el9
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=72
#SBATCH --time=05:00:00
#SBATCH --output=sbatch-output/output/antidunes_cylinder_sedimentation_showcase_noanalytical_125floaters_realgravity_twopillars_flow_%j.out
#SBATCH --error=sbatch-output/error/antidunes_cylinder_sedimentation_showcase_noanalytical_125floaters_realgravity_twopillars_flow_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# Cylinder sibling of run_144_box_sedimentation_showcase_analytical_125floaters_realgravity_twopillars_flow.sbatch
# (job 4112794, completed cleanly, full 120000 steps, no crash, floaters stayed under Z=63 after
# release). Every FloatingBox swapped for a FloatingCylinder (radius 3, length 18) and
# useAnalyticalContactDetection set false (pure GJK/EPA for every Cylinder x Cylinder pair) --
# everything else (domain, pillars, density mix, release schedule, timesteps) identical. Tests
# whether the collapsed-plane experiments' oversized MEM "phantom lift" force (SimpleBB.h's force
# computation has zero shape-type branching -- confirmed by code inspection) is shape-agnostic, by
# checking whether this known-good, non-collapsed scenario stays equally clean with cylinders. See
# the prm's own header.
if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    REPO="$(cd "$SLURM_SUBMIT_DIR/../../.." && pwd)"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    REPO="$(cd "$SCRIPT_DIR/../../../.." && pwd)"
fi
BIN_DIR="$REPO/build/apps/showcases/Antidunes"

"$REPO/apps/showcases/Antidunes/antidunes.sh" \
    --build --jobs "$(nproc)" \
    --run --np 144 \
    --prm "$REPO/apps/showcases/Antidunes/prms/cylinder_sedimentation_showcase_noanalytical_125floaters_realgravity_twopillars_flow_gas.prm"

STAMP="$(date +%Y%m%d_%H%M%S)"
if [[ -n "$HPCVAULT" ]]; then
    ARCHIVE_DIR="$HPCVAULT/antidunes-archive/${SLURM_JOB_NAME}_$STAMP"
else
    ARCHIVE_DIR="$BIN_DIR/archive/$STAMP"
fi
mkdir -p "$ARCHIVE_DIR"
# `if` (not `[[ ]] &&`) so a missing directory's false test doesn't leak into this being the
# script's own final exit status -- see the 125-floater box job's own comment (job 4104656) for the
# full mechanism this avoids.
for dir in vtk-out eval snapshots floater-history; do
    if [[ -d "$BIN_DIR/$dir" ]]; then mv "$BIN_DIR/$dir" "$ARCHIVE_DIR/$dir"; fi
done
