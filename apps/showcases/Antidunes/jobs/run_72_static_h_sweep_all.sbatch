#!/bin/bash -l
#SBATCH --job-name=antidunes-static-h-sweep-all
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=04:00:00
#SBATCH --output=sbatch-output/output/antidunes_static_h_sweep_all_%j.out
#SBATCH --error=sbatch-output/error/antidunes_static_h_sweep_all_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# STATIC-H SENSITIVITY SWEEP -- combined run of all 8 variants (h1 / h44 / h48 / default /
# h64 / h68 / h180 / dynamic, see prms/static_h_sweep_*_gas.prm) in a single job, one
# floating sphere (radius 12, densityRatio 0.5) starting ALREADY at 50% submersion (centre
# z=56, right at the undisturbed waterline, no drop/fall phase) in quiescent still water, no
# flow, no pillars, for 25000 timesteps each. Builds the Antidunes binary ONCE (identical
# across all 8 variants -- only the prm's ModelParameters differ), then runs each variant
# in turn, archiving eval/vtk-out/snapshots into its own subdirectory so nothing gets
# overwritten between runs. At the end, parses every variant's own saved stdout log for the
# per-eval-step [FloaterForces] line (uid=... step=N subm=X% ... pos=<x,y,z> ...) into a
# clean per-variant CSV (step,submersion_pct,centreZ) for the time-vs-submersion% and
# time-vs-position graphs.
#
# domainSize <128, 60, 180> / cellsPerBlock <32, 20, 30> = 4x3x6 = 72 processes
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
ANTIDUNES_DIR="$REPO/apps/showcases/Antidunes"
BIN_DIR="$REPO/build/apps/showcases/Antidunes"

VARIANTS=(h1 h44 h48 default h64 h68 h180 dynamic)

echo "==> Building Antidunes (once, shared by all ${#VARIANTS[@]} variants)"
"$ANTIDUNES_DIR/antidunes.sh" --build --jobs "$(nproc)"

STAMP="$(date +%Y%m%d_%H%M%S)"
if [[ -n "$HPCVAULT" ]]; then
    ARCHIVE_ROOT="$HPCVAULT/antidunes-archive/${SLURM_JOB_NAME}_$STAMP"
else
    ARCHIVE_ROOT="$BIN_DIR/archive/$STAMP"
fi
mkdir -p "$ARCHIVE_ROOT"

for v in "${VARIANTS[@]}"; do
    echo
    echo "=================================================================="
    echo "==> Running variant: $v"
    echo "=================================================================="
    VARIANT_DIR="$ARCHIVE_ROOT/$v"
    mkdir -p "$VARIANT_DIR"

    "$ANTIDUNES_DIR/antidunes.sh" \
        --run --np 72 \
        --prm "$ANTIDUNES_DIR/prms/static_h_sweep_${v}_gas.prm" \
        2>&1 | tee "$VARIANT_DIR/run.log"

    for dir in vtk-out eval snapshots; do
        [[ -d "$BIN_DIR/$dir" ]] && mv "$BIN_DIR/$dir" "$VARIANT_DIR/$dir"
    done
done

echo
echo "=================================================================="
echo "==> All variants done. Parsing [FloaterForces] logs into CSVs"
echo "=================================================================="
for v in "${VARIANTS[@]}"; do
    LOG="$ARCHIVE_ROOT/$v/run.log"
    CSV="$ARCHIVE_ROOT/$v/floater_timeseries.csv"
    echo "step,submersion_pct,centreZ" > "$CSV"
    grep '\[FloaterForces\]' "$LOG" \
        | grep -oP 'step=\K\d+|subm=\K[0-9.]+|pos=<\K[0-9.eE+-]+,[0-9.eE+-]+,[0-9.eE+-]+' \
        | paste -d',' - - - \
        | awk -F',' '{print $1","$2","$5}' >> "$CSV"
    echo "   wrote $CSV ($(($(wc -l < "$CSV") - 1)) data points)"
done

echo
echo "==> Archive root: $ARCHIVE_ROOT"
ls -la "$ARCHIVE_ROOT"
