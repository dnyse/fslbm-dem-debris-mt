#!/bin/bash -l
#SBATCH --job-name=antidunes-box-sedimentation-showcase-analytical-125floaters-realgravity
#SBATCH --partition=multinode
#SBATCH --constraint=el9
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=72
#SBATCH --time=04:00:00
#SBATCH --output=sbatch-output/output/antidunes_box_sedimentation_showcase_analytical_125floaters_realgravity_%j.out
#SBATCH --error=sbatch-output/error/antidunes_box_sedimentation_showcase_analytical_125floaters_realgravity_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# Genuine 5x5x5 grid (125 floaters) sibling of run_72_box_sedimentation_showcase_analytical_18floaters.sbatch:
# same box size (halfExtents <3,9,3>), domain grown <128,60,180> -> <160,160,180> to give real Y-neighbour
# spacing (deliberately reintroducing Y-neighbour contacts, unlike the 9-/18-floater files' single-Y-position
# design). cellsPerBlock <40,40,20> = 4x4x9 = 144 processes -- matches the project's existing 144-core
# precedent (run_144_thousand_spheres_two_pillars_slitpass_collapsed.sbatch): multinode partition, 2 nodes,
# 72 tasks/node. gravityOverride set to true g_SI=9.81 m/s^2 (checked against the LBM incompressibility
# bound and confirmed stable at this domain's scale -- see prm header). See
# prms/box_sedimentation_showcase_analytical_125floaters_realgravity_gas.prm's header for the full derivation.
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
    --prm "$REPO/apps/showcases/Antidunes/prms/box_sedimentation_showcase_analytical_125floaters_realgravity_gas.prm"

STAMP="$(date +%Y%m%d_%H%M%S)"
if [[ -n "$HPCVAULT" ]]; then
    ARCHIVE_DIR="$HPCVAULT/antidunes-archive/${SLURM_JOB_NAME}_$STAMP"
else
    ARCHIVE_DIR="$BIN_DIR/archive/$STAMP"
fi
mkdir -p "$ARCHIVE_DIR"
# `if` (not `[[ ]] &&`) so a missing directory's false test doesn't leak into this being the
# script's own final exit status -- bash's `for` returns its last command's exit code, and a bare
# `[[ -d ]] && mv` as the LAST statement of the LAST iteration (floater-history, only created when
# logFloaterHistory is on) would make the whole job report FAILED/ExitCode 1 to SLURM even when the
# simulation itself completed cleanly (confirmed happening exactly this way for job 4104656).
for dir in vtk-out eval snapshots floater-history; do
    if [[ -d "$BIN_DIR/$dir" ]]; then mv "$BIN_DIR/$dir" "$ARCHIVE_DIR/$dir"; fi
done
