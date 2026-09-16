#!/bin/bash -l
#SBATCH --job-name=antidunes-box-stability-curve-all
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=08:00:00
#SBATCH --output=sbatch-output/output/antidunes_box_stability_curve_all_%j.out
#SBATCH --error=sbatch-output/error/antidunes_box_stability_curve_all_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# STABILITY CURVE (Bogner & Rude 2013, Fig. 5.4/5.5 methodology) -- combined run of all 10
# heel-angle variants (heel00/05/10/15/20/25/30/35/40/45, see
# prms/box_stability_curve_heel*_gas.prm) in a single job. Each holds the Fekken densityRatio
# 0.5 cuboid (b:h=4:4 square cross-section, refined 1.5x resolution, shortened domain -- see
# box_density_050_refined_gas.prm) at a FIXED heel angle (pinRotationAxes <1,1,1>, ALL
# rotation pinned) while Z translation is free to find the buoyant depth at that angle, for
# 25000 timesteps. Builds the Antidunes binary ONCE (identical across all 10 variants -- only
# initialTiltAngleDeg differs), then runs each variant in turn, archiving eval/vtk-out/
# snapshots into its own subdirectory so nothing gets overwritten. At the end, parses every
# variant's own saved stdout log for the per-eval-step [BoxRightingMoment] line (heelDeg=...
# M_y=... GZ=... submersion=...) into a clean per-variant CSV
# (log_line,heelDeg,M_y,GZ,submersion_pct) -- log_line is the log file's own line number for
# that entry (NOT the simulation timestep -- [BoxRightingMoment] has no step= field -- but it
# is monotonically increasing and evenly spaced by evaluationFrequency, a valid ordinal proxy
# for time). The settled tail of each CSV gives one (heelDeg, M_y) point on the
# Torque-vs-Heeling-Angle stability curve; the CSVs' own time series also confirm heelDeg
# stays pinned and M_y/GZ converge rather than drift.
#
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
ANTIDUNES_DIR="$REPO/apps/showcases/Antidunes"
BIN_DIR="$REPO/build/apps/showcases/Antidunes"

ANGLES=(00 05 10 15 20 25 30 35 40 45)

echo "==> Building Antidunes (once, shared by all ${#ANGLES[@]} heel-angle variants)"
"$ANTIDUNES_DIR/antidunes.sh" --build --jobs "$(nproc)"

STAMP="$(date +%Y%m%d_%H%M%S)"
if [[ -n "$HPCVAULT" ]]; then
    ARCHIVE_ROOT="$HPCVAULT/antidunes-archive/${SLURM_JOB_NAME}_${SLURM_JOB_ID}_$STAMP"
else
    ARCHIVE_ROOT="$BIN_DIR/archive/$STAMP"
fi
mkdir -p "$ARCHIVE_ROOT"

for a in "${ANGLES[@]}"; do
    echo
    echo "=================================================================="
    echo "==> Running heel angle: $a deg"
    echo "=================================================================="
    VARIANT_DIR="$ARCHIVE_ROOT/heel${a}"
    mkdir -p "$VARIANT_DIR"

    "$ANTIDUNES_DIR/antidunes.sh" \
        --run --np 72 \
        --prm "$ANTIDUNES_DIR/prms/box_stability_curve_heel${a}_gas.prm" \
        2>&1 | tee "$VARIANT_DIR/run.log"

    for dir in vtk-out eval snapshots; do
        [[ -d "$BIN_DIR/$dir" ]] && mv "$BIN_DIR/$dir" "$VARIANT_DIR/$dir"
    done
done

echo
echo "=================================================================="
echo "==> All heel angles done. Parsing [BoxRightingMoment] logs into CSVs"
echo "=================================================================="
for a in "${ANGLES[@]}"; do
    LOG="$ARCHIVE_ROOT/heel${a}/run.log"
    CSV="$ARCHIVE_ROOT/heel${a}/righting_moment_timeseries.csv"
    echo "log_line,heelDeg,M_y,GZ,submersion_pct" > "$CSV"
    grep -n '\[BoxRightingMoment\]' "$LOG" \
        | awk -F: '{print $1}' \
        | paste -d',' - <(grep '\[BoxRightingMoment\]' "$LOG" \
            | grep -oP 'heelDeg=\K[0-9.eE+-]+|M_y=\K[0-9.eE+-]+|GZ=\K[0-9.eE+-]+|submersion=\K[0-9.eE+-]+' \
            | paste -d',' - - - -) >> "$CSV"
    echo "   wrote $CSV ($(($(wc -l < "$CSV") - 1)) data points)"
done

echo
echo "==> Archive root: $ARCHIVE_ROOT"
ls -la "$ARCHIVE_ROOT"
