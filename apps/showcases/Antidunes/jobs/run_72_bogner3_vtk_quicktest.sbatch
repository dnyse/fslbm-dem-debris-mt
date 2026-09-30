#!/bin/bash -l
#SBATCH --job-name=antidunes-vtk-quicktest
#SBATCH --partition=singlenode
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=00:15:00
#SBATCH --output=sbatch-output/output/vtk_quicktest_%j.out
#SBATCH --error=sbatch-output/error/vtk_quicktest_%j.err

# [Claude, 2026-09-29] THROWAWAY validation job, NOT thesis data -- per "on a different mirror just
# try on a single node 10 min run, easier to get resources to test this": confirms the new discrete
# fluid_field_snapshots VTK feature fires at the right steps with real data, cheaply. Single node (72
# ranks, cellsPerBlock <32,30,30>), 3000 steps, steps list 0,250,...,3000. No vault archiving (this
# run isn't meant to be kept).
#   sbatch jobs/run_72_bogner3_vtk_quicktest.sbatch
set -uo pipefail
if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    REPO="$(cd "$SLURM_SUBMIT_DIR/../../.." && pwd)"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    REPO="$(cd "$SCRIPT_DIR/../../../.." && pwd)"
fi
BIN_DIR="$REPO/build/apps/showcases/Antidunes"
SH="$REPO/apps/showcases/Antidunes/antidunes.sh"
PRM="$REPO/apps/showcases/Antidunes/prms/_quicktest_vtk_snapshots_sphere_gas.prm"
OUT="$REPO/apps/showcases/Antidunes/pinscan_results/${SLURM_JOB_ID:-manual_$(date +%Y%m%d_%H%M%S)}"
mkdir -p "$(dirname "$OUT")"
mkdir "$OUT" || { echo "results dir $OUT already exists -- refusing to overwrite"; exit 1; }
rm -rf "$BIN_DIR/density-column" "$BIN_DIR/eval" "$BIN_DIR/vtk-out" "$BIN_DIR/water-height" "$BIN_DIR/snapshots_vtk_quicktest"

"$SH" --build --jobs "$(nproc)" --bed --bed-np 72

TMPLOG="$OUT/quicktest.fulllog.tmp"
"$SH" --run --np 72 --prm "$PRM" > "$TMPLOG" 2>&1
echo "exit=$?" > "$OUT/quicktest.status"
grep -E "\[VTK\] fluid_field_snapshots|Assertion|FATAL" "$TMPLOG" > "$OUT/quicktest.vtk_and_errors.txt"
find "$BIN_DIR/vtk-out" -type f 2>/dev/null | sort > "$OUT/quicktest.vtk_file_listing.txt"
gzip -c "$TMPLOG" > "$OUT/quicktest.fulllog.gz" && rm -f "$TMPLOG"
echo "## DONE. Results in $OUT (vtk-out listing: $OUT/quicktest.vtk_file_listing.txt)"
