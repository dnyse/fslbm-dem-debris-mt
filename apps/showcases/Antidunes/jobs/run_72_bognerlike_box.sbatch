#!/bin/bash -l
#SBATCH --job-name=antidunes-bognerlike-box
#SBATCH --partition=singlenode
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=04:00:00
#SBATCH --output=sbatch-output/output/bognerlike_box_%j.out
#SBATCH --error=sbatch-output/error/bognerlike_box_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# [Claude, 2026-09-28] Bogner-like free box at real gravity (gen_bognerlike_box_scan.py): 4 cases sequentially. Submit as
# `cd <checkout>/apps/showcases/Antidunes && sbatch [--dependency=afterany:<job in the same checkout>] jobs/run_72_bognerlike_box.sbatch`.
# Per-job results dir (never overwritten), read-only backup in ~/pinscan_backup, NO vault archive (file quota exhausted).
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
BACKUP="$HOME/pinscan_backup/$(basename "$OUT")_bognerlike"

"$SH" --build --jobs "$(nproc)" --bed --bed-np 72

run_case() {  # $1 banner, $2 prm, $3 label
    echo "## $1"
    rm -rf "$BIN_DIR/density-column" "$BIN_DIR/eval" "$BIN_DIR/vtk-out" "$BIN_DIR/snapshots"
    TMPLOG="$OUT/$3.fulllog.tmp"
    "$SH" --run --np 72 --prm "$PRM_DIR/$2" > "$TMPLOG" 2>&1
    echo "exit=$?" > "$OUT/$3.status"
    grep -E "\[BoxForces\]|\[BoxTorque\]" "$TMPLOG" | cut -c1-460 > "$OUT/$3.forces.txt"
    grep -E "Assertion|FATAL" "$TMPLOG" | head -5 >> "$OUT/$3.status"
    cat "$BIN_DIR"/density-column/rank_*.csv 2>/dev/null | awk -F, '$1 ~ /^[0-9]+$/ { s[$1]+=$8; r[$5]+=$7; c[$5]++ } END { for (k in s) { n++; t+=s[k] } if (n) printf "waterline_mean %.5f nsteps %d\n", t/n, n; for (z in r) printf "rho %s %.6f\n", z, r[z]/c[z] }' | sort -k1,1 -k2,2n > "$OUT/$3.column.txt"
    cp "$BIN_DIR/eval/computedParams.txt" "$OUT/$3.computedParams.txt" 2>/dev/null
    rm -f "$TMPLOG"
    rm -rf "$BIN_DIR/density-column" "$BIN_DIR/eval" "$BIN_DIR/vtk-out" "$BIN_DIR/snapshots"
}

run_case "rho0.75 seed 2.86" bognerlike_r075_h286_gas.prm bl_r075_h286
run_case "rho0.75 seed 5.0" bognerlike_r075_h500_gas.prm bl_r075_h500
run_case "rho0.25 seed 2.86" bognerlike_r025_h286_gas.prm bl_r025_h286
run_case "rho0.25 seed 5.0" bognerlike_r025_h500_gas.prm bl_r025_h500

mkdir -p "$BACKUP" && cp -n "$OUT"/* "$BACKUP"/ && chmod -R a-w "$OUT" "$BACKUP"
echo "## DONE. Results in $OUT (backup: $BACKUP)"
