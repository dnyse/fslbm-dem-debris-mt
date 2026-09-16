#!/bin/bash -l
#SBATCH --job-name=antidunes-tsunami-pole-forest-driftwood-proportional
#SBATCH --partition=multinode
#SBATCH --constraint=el9
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=72
#SBATCH --time=08:00:00
#SBATCH --output=sbatch-output/output/tsunami_pole_forest_driftwood_proportional_%j.out
#SBATCH --error=sbatch-output/error/tsunami_pole_forest_driftwood_proportional_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# 4-node / 288-core (multinode partition) run of prms/tsunami_pole_forest_driftwood_proportional_gas.prm
# -- a proportionally-scaled version of the tsunami dam-break scenario: every dimension in that prm is
# Watanabe et al. 2023's own Section 4.2 number * 1.2 (the factor that keeps pole diameter at the same
# 6 LU/6 cm resolution already used in the non-proportional scenario). domainSize <1872,72,72>
# (9,704,448 LU cells, ~14x the non-proportional scenario's 691,200) with cellsPerBlock <26,36,36> gives
# exactly (1872/26)*(72/36)*(72/36) = 72*2*2 = 288 blocks, matching this job's rank count. See the prm's
# own header for the full proportionality derivation and what is/isn't rigorously Watanabe-scaled.
#
# --time=08:00:00 is a rough first guess: the non-proportional scenario (691,200 cells, 100000 steps) ran
# in ~35 min on 72 ranks; this domain is ~14x more cells at the same step count (120000, only ~1.2x more
# steps), so per-step cost should dominate -- with ~4x the ranks (288 vs 72) roughly offsetting the ~14x
# cell increase, a few hours seems like the right order of magnitude, but this has NOT actually been run
# yet, so treat this walltime as generous headroom, not a measured estimate.
#
# Do NOT run this without first checking: (a) available disk quota on the target filesystem -- the prm's
# own VTK/eval comments estimate ~40GB of vtk-out alone for this run, learned the hard way from run 4181010
# crashing mid-checkpoint on a nearly-full /home/hpc quota; (b) that $BIN_DIR/archive and any stale prior
# run output has been cleared first (see the pre-run cleanup block below, which handles the common case
# automatically, but a full quota check is still worth doing manually before a multi-hour, multi-node
# allocation).
set -euo pipefail

if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    REPO="$(cd "$SLURM_SUBMIT_DIR/../../.." && pwd)"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    REPO="$(cd "$SCRIPT_DIR/../../../.." && pwd)"
fi
BIN_DIR="$REPO/build/apps/showcases/Antidunes"

for dir in vtk-out eval snapshots floater-history; do
    if [[ -d "$BIN_DIR/$dir" ]]; then
        echo "WARNING: removing leftover $BIN_DIR/$dir from an unarchived previous run" >&2
        rm -rf "${BIN_DIR:?}/$dir"
    fi
done

"$REPO/apps/showcases/Antidunes/antidunes.sh" \
    --configure \
    --build-tsunami --jobs "$(nproc)" \
    --run-tsunami --np 288 \
    --tsunami-prm "$REPO/apps/showcases/Antidunes/prms/tsunami_pole_forest_driftwood_proportional_gas.prm"

# Archive destination: $HPCVAULT when set and not over its own quota, else $BIN_DIR/archive on /home/hpc.
# See run_72_tsunami_pole_forest_driftwood.sbatch's own comment here for the disk-quota incident this
# guards against -- keep archive/ on /home/hpc trimmed between runs; nothing in this script does that
# automatically.
STAMP="$(date +%Y%m%d_%H%M%S)"
if [[ -n "${HPCVAULT:-}" ]] && mkdir -p "$HPCVAULT/antidunes-archive/${SLURM_JOB_NAME}_${SLURM_JOB_ID}_$STAMP" 2>/dev/null; then
    ARCHIVE_DIR="$HPCVAULT/antidunes-archive/${SLURM_JOB_NAME}_$STAMP"
else
    ARCHIVE_DIR="$BIN_DIR/archive/$STAMP"
    mkdir -p "$ARCHIVE_DIR"
fi
for dir in vtk-out eval snapshots floater-history; do
    if [[ -d "$BIN_DIR/$dir" ]]; then
        mv "$BIN_DIR/$dir" "$ARCHIVE_DIR/$dir"
    fi
done
