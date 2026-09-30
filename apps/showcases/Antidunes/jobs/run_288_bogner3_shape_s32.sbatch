#!/bin/bash -l
#SBATCH --job-name=antidunes-bogner3-shape-s32
#SBATCH --partition=multinode
#SBATCH --nodes=4
#SBATCH --ntasks=288
#SBATCH --time=02:00:00
#SBATCH --output=sbatch-output/output/bogner3_shape_s32_%j.out
#SBATCH --error=sbatch-output/error/bogner3_shape_s32_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# [Claude, 2026-09-29] Companion to run_288_bogner3_deep_s32.sbatch, for the free-floating SPHERE/CYLINDER
# equilibrium runs at the same s=32 basin/resolution/gravity calibration as the bogner3 cuboid family
# (per "setup the equilibrium runs for this but for the cylinder and the sphere again for all 3
# densities, they can all start from the same z position, adjust the sizes accordingly, run for 100k
# with the waterline measurement etc included like in the cuboid case, archived in antidunes-archive").
# Single continuous 100000-step run only (no resume variant -- the cuboid family's own resume-artifact
# finding argues against resuming these at all). Fully free 6-DOF float, no pinning (prms have no
# pinPosition/pinRotation fields at all -- both default false).
#   sphere rho .75: sbatch --export=ALL,B3_RHO=075,B3_SHAPE=sphere jobs/run_288_bogner3_shape_s32.sbatch
#   cyl    rho .5 : sbatch --export=ALL,B3_RHO=050,B3_SHAPE=cyl    jobs/run_288_bogner3_shape_s32.sbatch
# Results in pinscan_results/<jobid>/ (never overwritten, backup in ~/pinscan_backup), snapshots kept
# in build/.../snapshots_b3_r<RHO>_<SHAPE>_s32_100k (NOT deleted), eval sqlite kept.
set -uo pipefail
RHO="${B3_RHO:?set B3_RHO=075|050|025 via --export=ALL,...}"
SHAPE="${B3_SHAPE:?set B3_SHAPE=sphere|cyl via --export=ALL,...}"
if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    REPO="$(cd "$SLURM_SUBMIT_DIR/../../.." && pwd)"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    REPO="$(cd "$SCRIPT_DIR/../../../.." && pwd)"
fi
BIN_DIR="$REPO/build/apps/showcases/Antidunes"
SH="$REPO/apps/showcases/Antidunes/antidunes.sh"
PRM_DIR="$REPO/apps/showcases/Antidunes/prms"
OUT="$REPO/apps/showcases/Antidunes/pinscan_results/${SLURM_JOB_ID:-manual_$(date +%Y%m%d_%H%M%S)}"
SNAP="$BIN_DIR/snapshots_b3_r${RHO}_${SHAPE}_s32_100k"
mkdir -p "$(dirname "$OUT")"
mkdir "$OUT" || { echo "results dir $OUT already exists -- refusing to overwrite"; exit 1; }
BACKUP="$HOME/pinscan_backup/$(basename "$OUT")_bogner3_r${RHO}_${SHAPE}_s32_100k"
PRM="bogner3_r${RHO}_${SHAPE}_s32_100k_gas.prm"; LAB="b3_r${RHO}_${SHAPE}_s32_100k"
[[ -f "$PRM_DIR/$PRM" ]] || { echo "FATAL: $PRM_DIR/$PRM does not exist"; exit 1; }
[[ -e "$SNAP" ]] && { echo "FATAL: $SNAP already exists -- refusing to start a fresh run over an existing checkpoint (move it away first)"; exit 1; }
rm -rf "$BIN_DIR/density-column" "$BIN_DIR/eval" "$BIN_DIR/vtk-out" "$BIN_DIR/water-height"

"$SH" --build --jobs "$(nproc)" --bed --bed-np 288

TMPLOG="$OUT/$LAB.fulllog.tmp"
"$SH" --run --np 288 --prm "$PRM_DIR/$PRM" > "$TMPLOG" 2>&1
echo "exit=$?" > "$OUT/$LAB.status"
if [[ "$SHAPE" == "cyl" ]]; then
    grep -E "\[CylForces\]|\[CylTorque\]" "$TMPLOG" | cut -c1-460 > "$OUT/$LAB.forces.txt"
else
    grep -E "\[FloaterForces\]" "$TMPLOG" | cut -c1-460 > "$OUT/$LAB.forces.txt"   # sphere: no torque log, no orientation
fi
grep -E "Assertion|FATAL" "$TMPLOG" | head -5 >> "$OUT/$LAB.status"
grep -c "WARNING" "$TMPLOG" | sed 's/^/warnings=/' >> "$OUT/$LAB.status"
cat "$BIN_DIR"/density-column/rank_*.csv 2>/dev/null | awk -F, '$1 ~ /^[0-9]+$/ { s[$1]+=$8; r[$5]+=$7; c[$5]++ } END { for (k in s) { n++; t+=s[k] } if (n) printf "waterline_mean %.5f nsteps %d\n", t/n, n; for (z in r) printf "rho %s %.6f\n", z, r[z]/c[z] }' | sort -k1,1 -k2,2n > "$OUT/$LAB.column.txt"
cat "$BIN_DIR"/water-height/rank_*.csv 2>/dev/null | awk -F, '$1 ~ /^[0-9]+$/ { k=$1" "$3" "$4; h[k]+=$5; o[k]+=$6 } END { for (k in h) print k, h[k], o[k] }' | sort -k1,1n -k2,2n -k3,3n > "$OUT/$LAB.waterheight.txt"   # step gx gy height nObstacleCells
cp "$BIN_DIR/eval/computedParams.txt" "$OUT/$LAB.computedParams.txt" 2>/dev/null
cp "$BIN_DIR/eval/antidunes.sqlite" "$OUT/$LAB.antidunes.sqlite" 2>/dev/null
gzip -c "$TMPLOG" > "$OUT/$LAB.fulllog.gz" && rm -f "$TMPLOG"
rm -rf "$BIN_DIR/density-column" "$BIN_DIR/eval" "$BIN_DIR/vtk-out" "$BIN_DIR/water-height"

mkdir -p "$BACKUP" && cp -n "$OUT"/* "$BACKUP"/ && chmod -R a-w "$OUT" "$BACKUP"

# [Claude, 2026-09-29] AUTO-ARCHIVE to the HPC vault (per "it should automatically archive into
# antidunes[-archive] ... so we can dependency chain them" -- no manual archiving step after each
# job). One shared family folder + PROVENANCE.md for all 6 sphere/cylinder equilibrium runs;
# flock-guarded append since jobs from different checkouts may finish close together.
VAULT_DIR="/home/vault/iwia/iwia148h/antidunes-archive/antidunes-bogner3-shapes-s32_20260929"
if mkdir -p "$VAULT_DIR" 2>/dev/null; then
    JOBARCH="$VAULT_DIR/${SLURM_JOB_ID:-manual}_${LAB}"
    mkdir -p "$JOBARCH" && cp "$OUT"/* "$JOBARCH"/ 2>/dev/null
    EXITCODE=$(head -1 "$OUT/$LAB.status" 2>/dev/null | sed 's/exit=//')
    WARNCOUNT=$(grep '^warnings=' "$OUT/$LAB.status" 2>/dev/null | sed 's/warnings=//')
    GEOM=$(grep -E "position|radius|length|densityRatio" "$PRM_DIR/$PRM" | grep -v "//" | sed 's/^[[:space:]]*/  /')
    (
        flock -x 200
        {
            echo ""
            echo "## ${SLURM_JOB_ID:-manual}_${LAB}"
            echo ""
            echo "- shape=$SHAPE rho=0.${RHO:1} exit=${EXITCODE:-unknown} warnings=${WARNCOUNT:-unknown}"
            echo "- prm: prms/$PRM"
            echo "$GEOM"
            echo "- archived: $JOBARCH"
        } >> "$VAULT_DIR/PROVENANCE.md"
    ) 200>"$VAULT_DIR/.provenance.lock"
    echo "## Also archived to $JOBARCH"
else
    echo "WARNING: could not create $VAULT_DIR -- NOT auto-archived to vault, only in $OUT / $BACKUP" >> "$OUT/$LAB.status"
fi

echo "## DONE. Results in $OUT (backup: $BACKUP); checkpoints kept in $SNAP"
