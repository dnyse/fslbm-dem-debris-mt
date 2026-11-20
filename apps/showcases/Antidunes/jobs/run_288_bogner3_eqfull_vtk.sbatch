#!/bin/bash -l
#SBATCH --job-name=b3-eqfull
#SBATCH --partition=multinode
#SBATCH --nodes=4
#SBATCH --ntasks=288
#SBATCH --time=03:00:00
#SBATCH --output=sbatch-output/output/bogner3_eqfull_vtk_%j.out
#SBATCH --error=sbatch-output/error/bogner3_eqfull_vtk_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# Floating-body equilibrium runs after Bogner (2017) Sec. 5.3.3, full length with whole-domain VTK: a free
# cuboid ("deep"), sphere or cylinder at densityRatio 0.75 / 0.50 / 0.25, s = 32, in a still-water basin
# periodic in x and y (prms/bogner3_r{075,050,025}_{deep,sphere,cyl}_s32_eqfull_vtk<VARIANT>_gas.prm, from
# prms/gen_bogner3_eqfull_vtk_prms.py). One case per job, each with its own job name, e.g.
#   sbatch --job-name=b3-eqfull-deep-r075 --export=ALL,B3_RHO=075,B3_SHAPE=deep jobs/run_288_bogner3_eqfull_vtk.sbatch
# B3_SHAPE is deep (default), sphere or cyl; B3_VARIANT is an optional prm suffix. The checkpoint folder is read
# from the prm (snapshotBaseFolder), and the job refuses to start over an existing one.
# Results go to pinscan_results/<jobid>/ (never overwritten) with a read-only copy in ~/pinscan_backup; the job
# exits with the simulation's exit code, so Slurm reports a crash.
# 100000 / 101000 steps take ~85 min on 288 ranks (82 min measured for 100000 steps), plus VTK and build: 3 h.
set -uo pipefail
RHO="${B3_RHO:?set B3_RHO=075|050|025 via --export=ALL,...}"
SHAPE="${B3_SHAPE:-deep}"   # deep (cuboid, default) | sphere | cyl (with B3_VARIANT=_freesharedstart)
VARIANT="${B3_VARIANT:-}"   # optional suffix, e.g. "_freesharedstart" -- default empty = unchanged behaviour
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
SNAP=""   # read from the prm below
mkdir -p "$(dirname "$OUT")"
mkdir "$OUT" || { echo "results dir $OUT already exists -- refusing to overwrite"; exit 1; }
BACKUP="$HOME/pinscan_backup/$(basename "$OUT")_bogner3_r${RHO}_${SHAPE}_s32_eqfull_vtk${VARIANT}"
PRM="bogner3_r${RHO}_${SHAPE}_s32_eqfull_vtk${VARIANT}_gas.prm"; LAB="b3_r${RHO}_${SHAPE}_s32_eqfull_vtk${VARIANT}"   # VARIANT e.g. _ypin
[[ -f "$PRM_DIR/$PRM" ]] || { echo "FATAL: $PRM_DIR/$PRM does not exist"; exit 1; }
SNAP="$BIN_DIR/$(sed -n 's#//.*##; s/^[[:space:]]*snapshotBaseFolder[[:space:]]\+\([^;[:space:]]*\);.*/\1/p' "$PRM_DIR/$PRM")"
[[ "$SNAP" != "$BIN_DIR/" ]] || { echo "FATAL: no snapshotBaseFolder in $PRM"; exit 1; }
[[ -e "$SNAP" ]] && { echo "FATAL: $SNAP already exists -- refusing to start a fresh run over an existing checkpoint (move it away first)"; exit 1; }
# Move leftovers aside instead of deleting them (never delete a previous run's output).
RESCUE="$BIN_DIR/archive/rescued_before_${SLURM_JOB_ID:-manual}_$(date +%Y%m%d_%H%M%S)"
for d in density-column eval vtk-out water-height; do
    if [[ -d "$BIN_DIR/$d" ]]; then mkdir -p "$RESCUE"; mv "$BIN_DIR/$d" "$RESCUE/$d"; fi
done

CONF=""; [[ -f "$REPO/build/CMakeCache.txt" ]] || CONF="--configure"   # fresh checkout: configure first
"$SH" $CONF --build --jobs "$(nproc)" --bed --bed-np 288

TMPLOG="$OUT/$LAB.fulllog.tmp"
"$SH" --run --np 288 --prm "$PRM_DIR/$PRM" > "$TMPLOG" 2>&1
SIMEXIT=$?   # also returned by the job (see the end), so Slurm reports a crash
echo "exit=$SIMEXIT" > "$OUT/$LAB.status"
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

# Whole-domain VTK is large (GBs, thousands of block files per run): the home backup leaves vtk-out out, the
# vault gets it as ONE tar file (vault file quota), and once the tar is verified vtk-out is removed from $OUT.
mkdir -p "$BACKUP" && rsync -a --exclude=vtk-out "$OUT"/ "$BACKUP"/ && chmod -R a-w "$BACKUP"

# Auto-archive to the vault, with one PROVENANCE.md entry per job
VAULT_DIR="/home/vault/iwia/iwia148h/antidunes-archive/antidunes-bogner3-eqfull-vtk_20261008"   # separate folder: the earlier one is packed as a tarball
if mkdir -p "$VAULT_DIR" 2>/dev/null; then
    JOBARCH="$VAULT_DIR/${SLURM_JOB_ID:-manual}_${LAB}"
    mkdir -p "$JOBARCH" && rsync -a --exclude=vtk-out "$OUT"/ "$JOBARCH"/
        # vtk-out as one tar; vtk-out is removed from $OUT only after the tar lists completely
    if tar -cf "$JOBARCH/vtk-out.tar" -C "$OUT" vtk-out && tar -tf "$JOBARCH/vtk-out.tar" > /dev/null \
       && [[ "$(tar -tf "$JOBARCH/vtk-out.tar" | grep -vc '/$')" -eq "$(find "$OUT/vtk-out" -type f | wc -l)" ]]; then
        rm -rf "$OUT/vtk-out"
        echo "vtk-out archived as $JOBARCH/vtk-out.tar (removed from $OUT)" >> "$OUT/$LAB.status"
    else
        echo "WARNING: vtk-out tar to the vault FAILED -- vtk-out kept in $OUT" >> "$OUT/$LAB.status"
    fi
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

chmod -R a-w "$OUT"
echo "## DONE. Results in $OUT (backup: $BACKUP); checkpoints kept in $SNAP"
exit "$SIMEXIT"
