#!/bin/bash -l
#SBATCH --job-name=antidunes-bogner3-deep-s32
#SBATCH --partition=multinode
#SBATCH --nodes=4
#SBATCH --ntasks=288
#SBATCH --time=02:00:00
#SBATCH --output=sbatch-output/output/bogner3_deep_s32_%j.out
#SBATCH --error=sbatch-output/error/bogner3_deep_s32_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# [Claude, 2026-09-29] Bogner-like free box at Bogner's HIGH-resolution box size s=32 (doubled basin 192x90x120, gen_bogner3_deep_scan.py
# family, hand-edited prm bogner3_r075_deep_s32_gas.prm). 288 ranks (6x6x8, cellsPerBlock <32,15,15>), 4 nodes (72 cores/node) -- not the sphere
# resolution scan's own 576/8-node scale, per explicit request to use 4 nodes here. Cheap first look: 51000 steps only, single checkpoint at snapshotFrequency=50000.
#   fresh : cd <checkout>/apps/showcases/Antidunes && sbatch --export=ALL,B3_RHO=075,B3_VARIANT=_s32 jobs/run_288_bogner3_deep_s32.sbatch
#   resume: add B3_RESUME=1 (needs bogner3_r075_deep_s32_resume_gas.prm, not yet created)
# Results in pinscan_results/<jobid>/ (never overwritten, backup in ~/pinscan_backup), snapshots kept in build/.../snapshots_b3_r<RHO>_s32/ (NOT deleted), eval sqlite kept.
set -uo pipefail
RHO="${B3_RHO:?set B3_RHO=075 via --export=ALL,B3_RHO=075,B3_VARIANT=_s32}"; RESUME="${B3_RESUME:-0}"; VAR="${B3_VARIANT:-_s32}"
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
SNAP="$BIN_DIR/snapshots_b3_r${RHO}${VAR}"
mkdir -p "$(dirname "$OUT")"
mkdir "$OUT" || { echo "results dir $OUT already exists -- refusing to overwrite"; exit 1; }
BACKUP="$HOME/pinscan_backup/$(basename "$OUT")_bogner3_r${RHO}${VAR}"
if [[ "$RESUME" == "1" ]]; then
    PRM="bogner3_r${RHO}_deep${VAR}_resume_gas.prm"; LAB="b3_r${RHO}${VAR}_resume"
    [[ -f "$SNAP/simState.file" ]] || { echo "FATAL: no $SNAP/simState.file to resume from"; exit 1; }
    cp -r "$SNAP" "${SNAP}_backup_before_resume_${SLURM_JOB_ID:-manual}" || { echo "FATAL: checkpoint backup failed"; exit 1; }
else
    PRM="bogner3_r${RHO}_deep${VAR}_gas.prm"; LAB="b3_r${RHO}${VAR}"
    [[ -e "$SNAP" ]] && { echo "FATAL: $SNAP already exists -- refusing to start a fresh run over an existing checkpoint (move it away or resume)"; exit 1; }
fi
rm -rf "$BIN_DIR/density-column" "$BIN_DIR/eval" "$BIN_DIR/vtk-out" "$BIN_DIR/water-height"

"$SH" --build --jobs "$(nproc)" --bed --bed-np 288

TMPLOG="$OUT/$LAB.fulllog.tmp"
"$SH" --run --np 288 --prm "$PRM_DIR/$PRM" > "$TMPLOG" 2>&1
echo "exit=$?" > "$OUT/$LAB.status"
grep -E "\[BoxForces\]|\[BoxTorque\]" "$TMPLOG" | cut -c1-460 > "$OUT/$LAB.forces.txt"
grep -E "Assertion|FATAL" "$TMPLOG" | head -5 >> "$OUT/$LAB.status"
grep -c "WARNING" "$TMPLOG" | sed 's/^/warnings=/' >> "$OUT/$LAB.status"
cat "$BIN_DIR"/density-column/rank_*.csv 2>/dev/null | awk -F, '$1 ~ /^[0-9]+$/ { s[$1]+=$8; r[$5]+=$7; c[$5]++ } END { for (k in s) { n++; t+=s[k] } if (n) printf "waterline_mean %.5f nsteps %d\n", t/n, n; for (z in r) printf "rho %s %.6f\n", z, r[z]/c[z] }' | sort -k1,1 -k2,2n > "$OUT/$LAB.column.txt"
cat "$BIN_DIR"/water-height/rank_*.csv 2>/dev/null | awk -F, '$1 ~ /^[0-9]+$/ { k=$1" "$3" "$4; h[k]+=$5; o[k]+=$6 } END { for (k in h) print k, h[k], o[k] }' | sort -k1,1n -k2,2n -k3,3n > "$OUT/$LAB.waterheight.txt"   # step gx gy height nObstacleCells
cp "$BIN_DIR/eval/computedParams.txt" "$OUT/$LAB.computedParams.txt" 2>/dev/null
cp "$BIN_DIR/eval/antidunes.sqlite" "$OUT/$LAB.antidunes.sqlite" 2>/dev/null
gzip -c "$TMPLOG" > "$OUT/$LAB.fulllog.gz" && rm -f "$TMPLOG"
rm -rf "$BIN_DIR/density-column" "$BIN_DIR/eval" "$BIN_DIR/vtk-out" "$BIN_DIR/water-height"

mkdir -p "$BACKUP" && cp -n "$OUT"/* "$BACKUP"/ && chmod -R a-w "$OUT" "$BACKUP"
echo "## DONE. Results in $OUT (backup: $BACKUP); checkpoints kept in $SNAP"
