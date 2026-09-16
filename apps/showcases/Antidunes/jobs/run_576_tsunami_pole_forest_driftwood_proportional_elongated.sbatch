#!/bin/bash -l
#SBATCH --job-name=antidunes-tsunami-pole-forest-driftwood-proportional-elongated
#SBATCH --partition=multinode
#SBATCH --constraint=el9
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=72
#SBATCH --time=10:00:00
#SBATCH --output=sbatch-output/output/tsunami_pole_forest_driftwood_proportional_elongated_%j.out
#SBATCH --error=sbatch-output/error/tsunami_pole_forest_driftwood_proportional_elongated_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# 8-node / 576-core (multinode partition) run of
# prms/tsunami_pole_forest_driftwood_proportional_elongated_gas.prm -- built on top of the (still unrun)
# proportional scenario with three changes: domainX doubled (1872->3744, giving 2420 LU of run-out past the
# forest instead of 548, so debris that clears the forest can keep travelling instead of immediately
# bouncing off the far wall), half the driftwood logs rotated 90 degrees (broadside-to-flow, for a more
# chaotic pile-up against the poles), and logFloaterHistory enabled (diagnostic -- see the prm's own header
# point 3). cellsPerBlock <26,36,36> on domain <3744,72,72> gives exactly
# (3744/26)*(72/36)*(72/36) = 144*2*2 = 576 blocks, matching this job's rank count.
#
# THIS RUN IS EXPECTED TO REPRODUCE A KNOWN BUG, NOT DEMONSTRATE FLOATING DEBRIS: analysis of the completed
# non-proportional 100k run (job 4181190) found the driftwood boxes do not actually float -- they sit flush
# on the floor under many LU of standing water, with the MEM hydrodynamic force (F_hydro_z, meant to supply
# ALL buoyant support under useMemOnlyBuoyancy=true) coming out NEGATIVE instead of positive for every box
# still touching the floor. The one box that got knocked clear of the floor (onto/against a pole) showed the
# correct sign. Root cause not yet fixed -- suspected to be that these boxes start resting on the dry floor
# from step 0 and get engulfed in place as water arrives, never receiving the initial impact/detachment
# event that every other validated floater scenario in this project gets (all of them drop boxes from the
# air into already-established water). logFloaterHistory is turned on specifically to get a much finer-
# grained diagnostic trail of this than the once-per-500-step BoxForces log used so far. Do not be surprised
# if the driftwood again just sits on the floor and gets shoved along it rather than floating.
#
# --time=10:00:00 is a rough first guess, scaled from the (also never-run) 4-node/9.7M-cell proportional
# scenario's own untested 8h guess: this domain is 2x the cells (19.4M) with 2x the ranks (576 vs 288), so
# per-rank workload is unchanged, but cross-node communication overhead across 8 nodes instead of 4 adds
# some uncertain overhead -- treat this as generous headroom, not a measured estimate.
#
# Do NOT run this without first checking available disk quota -- see run_72_tsunami_pole_forest_driftwood.sbatch's own comment for the /home/hpc quota-exhaustion incident this guards against. The pre-run
# cleanup below handles the common case (leftover output from an unarchived previous run) automatically.
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
    --run-tsunami --np 576 \
    --tsunami-prm "$REPO/apps/showcases/Antidunes/prms/tsunami_pole_forest_driftwood_proportional_elongated_gas.prm"

# Archive destination: $HPCVAULT when set and not over its own quota, else $BIN_DIR/archive on /home/hpc.
# See run_72_tsunami_pole_forest_driftwood.sbatch's own comment for the disk-quota incident this guards
# against -- keep archive/ on /home/hpc trimmed between runs; nothing in this script does that
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
