#!/bin/bash -l
#SBATCH --job-name=antidunes-tsunami-pole-forest-driftwood
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=03:00:00
#SBATCH --output=sbatch-output/output/tsunami_pole_forest_driftwood_%j.out
#SBATCH --error=sbatch-output/error/tsunami_pole_forest_driftwood_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# First run of the NEW `tsunami` executable (tsunami.cpp) -- a SEPARATE binary from Antidunes, added
# to recreate Watanabe et al. 2023's dry-domain dam-break tsunami test case (a reservoir behind a gate
# collapsing under gravity across a dry floor into a pole forest / driftwood field), closing the gap
# the wet-channel toy scenario (prms/toy_pole_forest_driftwood_gas.prm) explicitly flagged in its own
# header. See TsunamiDataInitialization.h / tsunami.cpp for the mechanism; prms/
# tsunami_pole_forest_driftwood_gas.prm for this scenario's own geometry/parameter derivation.
#
# antidunes.sh gained --build-tsunami/--run-tsunami flags (additive, alongside its existing
# --build/--run for the Antidunes target) specifically so this job can go through the same shared
# build/run helper every other scenario uses, instead of duplicating the module-load/cmake/mpirun
# sequence by hand in this script.
#
# --configure is included because CMakeLists.txt changed (new `tsunami` target) since this mirror's
# build directory was last configured -- cmake -B picks up the new target from the changed
# CMakeLists.txt; it is idempotent/harmless to run again for every other target too.
#
# SECOND submission: run 4180567 (first submission, 15000 steps) was thoroughly checked against both
# the log and VTK output (wave_front_check.py / find_degenerate_cells.py) and came back clean --
# reservoir-collapse, wave transit, and pole/floater interaction all matched expectations. That check
# also found (a) two FloatingBox entries in the prm genuinely overlapped at t=0 -- fixed in the prm now,
# and (b) the true wave-front speed is ~2x slower than the original shallow-water estimate -- timesteps
# raised from 15000 to 100000 (with writeFrequency scaled the same ~6.7x, so VTK volume stays ~9.5GB
# instead of ballooning to ~56GB) to actually watch the tank settle across multiple reflections. See the
# prm's own header for the full derivation.
#
# Also: run 4180567's OWN archived eval/avgDataSlices was found contaminated with ~267 leftover files
# (steps up to 103000) from some earlier, unrelated, never-cleaned-up run that predated it in this same
# $BIN_DIR -- cleaned out by hand after the fact. The pre-run cleanup below removes any such leftover
# vtk-out/eval/snapshots/floater-history from $BIN_DIR before this job starts, so a run that manages to
# skip its own archiving step (e.g. a crash) can never again contaminate the NEXT run's output the same
# way -- a real completed run always ends by archiving these away (see the tail of this script), so
# anything still present at start-of-job is leftover cruft, never live data to preserve.
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
    --run-tsunami --np 72 \
    --tsunami-prm "$REPO/apps/showcases/Antidunes/prms/tsunami_pole_forest_driftwood_gas.prm"

# Archive destination: $HPCVAULT when set and not over its own quota, else $BIN_DIR/archive on
# /home/hpc. run 4181010 crashed mid-checkpoint-write because /home/hpc had been driven to 207.7G of a
# 209.7G hard byte quota by earlier runs' archives piling up here after vault rejected them (vault was
# hard-blocking ALL new files that day due to a file-count quota breach: `Filec 201K > FileQ 200K` in the
# job epilogue's quota report, unrelated to its own byte usage, which was nowhere near its limit) -- fixed
# by clearing out /home/hpc directly (not by adding another storage tier here). Keep archive/ on
# /home/hpc trimmed between runs; nothing in this script does that automatically.
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
