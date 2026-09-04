#!/bin/bash -l
#SBATCH --job-name=antidunes-cylinder-sedimentation-showcase-pureanalytical-10k-floaterhistory
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=01:00:00
#SBATCH --output=sbatch-output/output/antidunes_cylinder_sedimentation_showcase_pureanalytical_10k_floaterhistory_%j.out
#SBATCH --error=sbatch-output/error/antidunes_cylinder_sedimentation_showcase_pureanalytical_10k_floaterhistory_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# "Pure Kodam+Guo" sibling of run_72_cylinder_sedimentation_showcase_analytical_10k_floaterhistory.sbatch:
# identical dense 100-log pile/settings, but usePureAnalyticalCylinder true -- removes the
# Cylinder-Cylinder GJK/EPA safety net entirely, trusting the closed-form Kodam et al. (2010) +
# Guo et al. (2012) dispatch as the sole answer, exactly as those papers' own algorithm stands on
# its own (neither paper uses GJK/EPA or any iterative solver). See
# prms/cylinder_sedimentation_showcase_pureanalytical_10k_floaterhistory_gas.prm's header.
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
    --prm "$REPO/apps/showcases/Antidunes/prms/cylinder_sedimentation_showcase_pureanalytical_10k_floaterhistory_gas.prm"

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
