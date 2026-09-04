#!/bin/bash -l
#SBATCH --job-name=antidunes-bridge-lv3-DEBUG-widerstacks-nogjk
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=02:30:00
#SBATCH --output=sbatch-output/output/antidunes_test_bridge_v1_raisedwater_bigjam_mixed_level3_DEBUG_widerstacks_nogjk_%j.out
#SBATCH --error=sbatch-output/error/antidunes_test_bridge_v1_raisedwater_bigjam_mixed_level3_DEBUG_widerstacks_nogjk_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# ANALYTIC-ONLY variant of the widerstacks scenario (job 4166652 crashed on a mesa_pd::mpi::ContactFilter
# double-accept race, Cylinder x Box, which is an unconditional-GJK/EPA shape pairing). This run removes
# every FloatingCylinder (54 entries) while keeping FloatingBoxes/FloatingSpheres, piers, and arches
# identical, so every DEM contact in the scenario goes through analytic (closed-form) detection only --
# Sphere x Sphere, Sphere x Box, Box x Box (useAnalyticalContactDetection=true), and Sphere/Box vs the
# pillars/arches (CylindricalBoundary, this showcase's own analytic override). See the prm's own header
# for the verification against mesa_pd/collision_detection/AnalyticContactDetection.h +
# GeneralContactDetection.h before this file was built. Goal: test whether removing GJK/EPA's iterative
# contact-point computation avoids (or reduces the frequency of) the ContactFilter ownership-arbitration
# race, while the real fix (the documented-but-unimplemented smallest-witnessing-rank tie-break in
# ContactFilter.h) is still pending.
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
    --prm "$REPO/apps/showcases/Antidunes/prms/test_bridge_v1_raisedwater_bigjam_mixed_level3_DEBUG_widerstacks_nogjk_gas.prm"

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
