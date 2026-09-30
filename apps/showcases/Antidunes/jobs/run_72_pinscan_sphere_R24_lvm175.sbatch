#!/bin/bash -l
#SBATCH --job-name=antidunes-pinscan-sphere-R24_lvm175
#SBATCH --partition=singlenode
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=01:00:00
#SBATCH --output=sbatch-output/output/pinscan_sphere_R24_lvm175_%j.out
#SBATCH --error=sbatch-output/error/pinscan_sphere_R24_lvm175_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# [Claude, 2026-09-28] Pinned-sphere MEM buoyancy scan, radius 24 (see gen_pinned_sphere_scan.py). Submit as
# `cd <checkout>/apps/showcases/Antidunes && sbatch jobs/run_72_pinscan_sphere_R24_lvm175.sbatch`.
# NO vault archive (its file quota is exhausted, which lost data on 2026-09-28): each case is reduced in-job to
# a few KB under pinscan_results/<jobid>/ (unique per job, never overwritten, backed up to ~/pinscan_backup, read-only) (FloaterForces lines, waterline, mean rho(z)); big outputs are deleted.
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
BACKUP="$HOME/pinscan_backup/$(basename "$OUT")_R24_lvm175"

"$SH" --build --jobs "$(nproc)" --bed --bed-np 72

run_case() {  # $1 banner, $2 prm, $3 label
    echo "## $1"
    rm -rf "$BIN_DIR/density-column" "$BIN_DIR/eval" "$BIN_DIR/vtk-out" "$BIN_DIR/snapshots"
    TMPLOG="$OUT/$3.fulllog.tmp"
    "$SH" --run --np 72 --prm "$PRM_DIR/$2" > "$TMPLOG" 2>&1
    echo "exit=$?" > "$OUT/$3.status"
    grep "\[FloaterForces\]" "$TMPLOG" | cut -c1-420 > "$OUT/$3.forces.txt"
    grep -E "Assertion|FATAL|NaN" "$TMPLOG" | head -5 >> "$OUT/$3.status"
    cat "$BIN_DIR"/density-column/rank_*.csv 2>/dev/null | awk -F, '$1 ~ /^[0-9]+$/ { s[$1]+=$8; r[$5]+=$7; c[$5]++ } END { for (k in s) { n++; t+=s[k] } if (n) printf "waterline_mean %.5f nsteps %d\n", t/n, n; for (z in r) printf "rho %s %.6f\n", z, r[z]/c[z] }' | sort -k1,1 -k2,2n > "$OUT/$3.column.txt"
    cp "$BIN_DIR/eval/computedParams.txt" "$OUT/$3.computedParams.txt" 2>/dev/null
    rm -f "$TMPLOG"
    rm -rf "$BIN_DIR/density-column" "$BIN_DIR/eval" "$BIN_DIR/vtk-out" "$BIN_DIR/snapshots"
}

run_case "R24 real f0.75_lvm175" pinscan_sphere_R24_real_f075_lvm175_gas.prm R24_real_f075_lvm175

mkdir -p "$BACKUP" && cp -n "$OUT"/* "$BACKUP"/ && chmod -R a-w "$OUT" "$BACKUP"   # 2nd copy outside the checkout; read-only
echo "## DONE. Results in $OUT (backup: $BACKUP)"
