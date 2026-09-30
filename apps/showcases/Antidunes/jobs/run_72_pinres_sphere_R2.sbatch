#!/bin/bash -l
#SBATCH --job-name=antidunes-pinres-sphere-R2
#SBATCH --partition=singlenode
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=01:30:00
#SBATCH --output=sbatch-output/output/pinres_sphere_R2_%j.out
#SBATCH --error=sbatch-output/error/pinres_sphere_R2_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# [Claude, 2026-09-28] Pinned-sphere RESOLUTION scan, R=2 (see gen_pinned_resolution_scan.py). Submit as
# `cd <checkout>/apps/showcases/Antidunes && sbatch jobs/run_72_pinres_sphere_R2.sbatch`. Needs spheres_out_buoyancy_r2.dat in the build dir.
# Results (per-job dir, never overwritten, read-only backup in ~/pinscan_backup, no vault archive -- its file quota is exhausted).
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
mkdir -p "$(dirname "$OUT")"
mkdir "$OUT" || { echo "results dir $OUT already exists -- refusing to overwrite"; exit 1; }
BACKUP="$HOME/pinscan_backup/$(basename "$OUT")_pinres_R2"
[[ -f "$BIN_DIR/spheres_out_buoyancy_r2.dat" ]] || { echo "missing bed file spheres_out_buoyancy_r2.dat"; exit 1; }

"$SH" --build --jobs "$(nproc)"

L=R2_pinres
rm -rf "$BIN_DIR/density-column" "$BIN_DIR/eval" "$BIN_DIR/vtk-out" "$BIN_DIR/snapshots"
TMPLOG="$OUT/$L.fulllog.tmp"
"$SH" --run --np 72 --prm "$PRM_DIR/pinres_sphere_R2_gas.prm" > "$TMPLOG" 2>&1
echo "exit=$?" > "$OUT/$L.status"
grep "\[FloaterForces\]" "$TMPLOG" | cut -c1-420 > "$OUT/$L.forces.txt"
grep -E "Assertion|FATAL" "$TMPLOG" | head -5 >> "$OUT/$L.status"
cat "$BIN_DIR"/density-column/rank_*.csv 2>/dev/null | awk -F, '$1 ~ /^[0-9]+$/ { s[$1]+=$8; r[$5]+=$7; c[$5]++ } END { for (k in s) { n++; t+=s[k] } if (n) printf "waterline_mean %.5f nsteps %d\n", t/n, n; for (z in r) printf "rho %s %.6f\n", z, r[z]/c[z] }' | sort -k1,1 -k2,2n > "$OUT/$L.column.txt"
cp "$BIN_DIR/eval/computedParams.txt" "$OUT/$L.computedParams.txt" 2>/dev/null
rm -f "$TMPLOG"
rm -rf "$BIN_DIR/density-column" "$BIN_DIR/eval" "$BIN_DIR/vtk-out" "$BIN_DIR/snapshots"
mkdir -p "$BACKUP" && cp -n "$OUT"/* "$BACKUP"/ && chmod -R a-w "$OUT" "$BACKUP"
echo "## DONE. Results in $OUT (backup: $BACKUP)"
