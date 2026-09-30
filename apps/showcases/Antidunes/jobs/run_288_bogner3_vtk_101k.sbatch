#!/bin/bash -l
#SBATCH --job-name=antidunes-bogner3-vtk-101k
#SBATCH --partition=multinode
#SBATCH --nodes=4
#SBATCH --ntasks=288
#SBATCH --time=02:10:00
#SBATCH --output=sbatch-output/output/bogner3_vtk_101k_%j.out
#SBATCH --error=sbatch-output/error/bogner3_vtk_101k_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# [Claude, 2026-09-29] VTK-enabled thesis-visualization companions to the bogner3 s=32 equilibrium
# family. Cuboid/sphere share ONE starting position across all 3 densities (position.z=60.0746, the
# water level) per "run the spheres and cuboids ... have them start ... from the same starting
# position"; cylinder keeps its own per-density starting z but now has pitch/yaw pinned
# (pinRotationAxes <1,0,1>) per "in the cylinder case have the tilt x and z pinned". VTK uses the new
# discrete fluid_field_snapshots mechanism (SimulationParams.h/VtkWriter.h, added this session) at
# steps 0,250,500,1000,2000,4000,8000,16000,32000,64000,100000 -- dense early (settling transient),
# one frame exactly at the intended end. 101000 total steps (not 100000) for safety margin. Single
# continuous run (no resume), real g_SI=9.81, guard off.
#   cuboid rho .75 : sbatch --export=ALL,B3_RHO=075,B3_SHAPE=deep   jobs/run_288_bogner3_vtk_101k.sbatch
#   sphere rho .5  : sbatch --export=ALL,B3_RHO=050,B3_SHAPE=sphere jobs/run_288_bogner3_vtk_101k.sbatch
#   cyl    rho .25 : sbatch --export=ALL,B3_RHO=025,B3_SHAPE=cyl    jobs/run_288_bogner3_vtk_101k.sbatch
# Results in pinscan_results/<jobid>/ (never overwritten, backup in ~/pinscan_backup), vtk-out COPIED
# into the results dir too (small: 11 frames x 1 XZ-slab each), snapshots kept in build/.../
# snapshots_b3_r<RHO>_<SHAPE>_s32_101k_vtk (NOT deleted), eval sqlite kept. Auto-archives to the vault
# on completion (same pattern as run_288_bogner3_shape_s32.sbatch).
set -uo pipefail
RHO="${B3_RHO:?set B3_RHO=075|050|025 via --export=ALL,...}"
SHAPE="${B3_SHAPE:?set B3_SHAPE=deep|sphere|cyl via --export=ALL,...}"
VARIANT="${B3_VARIANT:-}"   # [Claude, 2026-09-29] optional suffix, e.g. "_freesharedstart" -- default empty = unchanged behaviour
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
SNAP="$BIN_DIR/snapshots_b3_r${RHO}_${SHAPE}_s32_101k_vtk${VARIANT}"
mkdir -p "$(dirname "$OUT")"
mkdir "$OUT" || { echo "results dir $OUT already exists -- refusing to overwrite"; exit 1; }
BACKUP="$HOME/pinscan_backup/$(basename "$OUT")_bogner3_r${RHO}_${SHAPE}_s32_101k_vtk${VARIANT}"
PRM="bogner3_r${RHO}_${SHAPE}_s32_101k_vtk${VARIANT}_gas.prm"; LAB="b3_r${RHO}_${SHAPE}_s32_101k_vtk${VARIANT}"
[[ -f "$PRM_DIR/$PRM" ]] || { echo "FATAL: $PRM_DIR/$PRM does not exist"; exit 1; }
[[ -e "$SNAP" ]] && { echo "FATAL: $SNAP already exists -- refusing to start a fresh run over an existing checkpoint (move it away first)"; exit 1; }
rm -rf "$BIN_DIR/density-column" "$BIN_DIR/eval" "$BIN_DIR/vtk-out" "$BIN_DIR/water-height"

"$SH" --build --jobs "$(nproc)" --bed --bed-np 288

TMPLOG="$OUT/$LAB.fulllog.tmp"
"$SH" --run --np 288 --prm "$PRM_DIR/$PRM" > "$TMPLOG" 2>&1
echo "exit=$?" > "$OUT/$LAB.status"
if [[ "$SHAPE" == "cyl" ]]; then
    grep -E "\[CylForces\]|\[CylTorque\]" "$TMPLOG" | cut -c1-460 > "$OUT/$LAB.forces.txt"
elif [[ "$SHAPE" == "sphere" ]]; then
    grep -E "\[FloaterForces\]" "$TMPLOG" | cut -c1-460 > "$OUT/$LAB.forces.txt"
else
    grep -E "\[BoxForces\]|\[BoxTorque\]" "$TMPLOG" | cut -c1-460 > "$OUT/$LAB.forces.txt"
fi
grep -E "\[VTK\] fluid_field_snapshots" "$TMPLOG" > "$OUT/$LAB.vtk_schedule.txt"
grep -E "Assertion|FATAL" "$TMPLOG" | head -5 >> "$OUT/$LAB.status"
grep -c "WARNING" "$TMPLOG" | sed 's/^/warnings=/' >> "$OUT/$LAB.status"
cat "$BIN_DIR"/density-column/rank_*.csv 2>/dev/null | awk -F, '$1 ~ /^[0-9]+$/ { s[$1]+=$8; r[$5]+=$7; c[$5]++ } END { for (k in s) { n++; t+=s[k] } if (n) printf "waterline_mean %.5f nsteps %d\n", t/n, n; for (z in r) printf "rho %s %.6f\n", z, r[z]/c[z] }' | sort -k1,1 -k2,2n > "$OUT/$LAB.column.txt"
cat "$BIN_DIR"/water-height/rank_*.csv 2>/dev/null | awk -F, '$1 ~ /^[0-9]+$/ { k=$1" "$3" "$4; h[k]+=$5; o[k]+=$6 } END { for (k in h) print k, h[k], o[k] }' | sort -k1,1n -k2,2n -k3,3n > "$OUT/$LAB.waterheight.txt"
cp "$BIN_DIR/eval/computedParams.txt" "$OUT/$LAB.computedParams.txt" 2>/dev/null
cp "$BIN_DIR/eval/antidunes.sqlite" "$OUT/$LAB.antidunes.sqlite" 2>/dev/null
mkdir -p "$OUT/vtk-out" && cp -r "$BIN_DIR"/vtk-out/* "$OUT/vtk-out/" 2>/dev/null
gzip -c "$TMPLOG" > "$OUT/$LAB.fulllog.gz" && rm -f "$TMPLOG"
rm -rf "$BIN_DIR/density-column" "$BIN_DIR/eval" "$BIN_DIR/vtk-out" "$BIN_DIR/water-height"

mkdir -p "$BACKUP" && cp -rn "$OUT"/* "$BACKUP"/ && chmod -R a-w "$OUT" "$BACKUP"

# Auto-archive to the vault (same pattern as run_288_bogner3_shape_s32.sbatch)
VAULT_DIR="/home/vault/iwia/iwia148h/antidunes-archive/antidunes-bogner3-vtk-101k_20260929"
if mkdir -p "$VAULT_DIR" 2>/dev/null; then
    JOBARCH="$VAULT_DIR/${SLURM_JOB_ID:-manual}_${LAB}"
    mkdir -p "$JOBARCH" && cp -r "$OUT"/* "$JOBARCH"/ 2>/dev/null
    EXITCODE=$(head -1 "$OUT/$LAB.status" 2>/dev/null | sed 's/exit=//')
    WARNCOUNT=$(grep '^warnings=' "$OUT/$LAB.status" 2>/dev/null | sed 's/warnings=//')
    GEOM=$(grep -E "position|radius|length|halfExtents|densityRatio|pinRotation" "$PRM_DIR/$PRM" | grep -v "//" | sed 's/^[[:space:]]*/  /')
    (
        flock -x 200
        {
            echo ""
            echo "## ${SLURM_JOB_ID:-manual}_${LAB}"
            echo ""
            echo "- shape=$SHAPE rho=0.${RHO:1} exit=${EXITCODE:-unknown} warnings=${WARNCOUNT:-unknown}"
            echo "- prm: prms/$PRM"
            echo "$GEOM"
            echo "- archived: $JOBARCH (includes vtk-out/)"
        } >> "$VAULT_DIR/PROVENANCE.md"
    ) 200>"$VAULT_DIR/.provenance.lock"
    echo "## Also archived to $JOBARCH"
else
    echo "WARNING: could not create $VAULT_DIR -- NOT auto-archived to vault, only in $OUT / $BACKUP" >> "$OUT/$LAB.status"
fi

echo "## DONE. Results in $OUT (backup: $BACKUP); checkpoints kept in $SNAP"
