#!/bin/bash -l
#SBATCH --job-name=antidunes-cylcrash
#SBATCH --partition=multinode
#SBATCH --nodes=4
#SBATCH --ntasks=288
#SBATCH --time=00:40:00
#SBATCH --output=sbatch-output/output/cylcrash_%j.out
#SBATCH --error=sbatch-output/error/cylcrash_%j.err

# [Claude, 2026-10-07] crash diagnosis of the cylinder equilibrium rerun 4313587 (NaN at step 8733). Runs one prm of
# prms/gen_cylcrash_prms.py:  sbatch --export=ALL,CC_PRM=cylcrash_base_gas.prm jobs/run_288_cylcrash.sbatch
# Results (status, forces, full log) in pinscan_results/<jobid>/; the checkpoint snapshots_cylcrash stays in the build dir.
# Never deletes earlier output: leftovers of a previous run are moved aside.
set -uo pipefail
PRM="${CC_PRM:?set CC_PRM}"; LAB="${PRM%_gas.prm}"
REPO="$(cd "$SLURM_SUBMIT_DIR/../../.." && pwd)"
BIN_DIR="$REPO/build/apps/showcases/Antidunes"; SH="$REPO/apps/showcases/Antidunes/antidunes.sh"; PRM_DIR="$REPO/apps/showcases/Antidunes/prms"
OUT="$REPO/apps/showcases/Antidunes/pinscan_results/${SLURM_JOB_ID}"; mkdir -p "$(dirname "$OUT")"; mkdir "$OUT" || exit 1
RESCUE="$BIN_DIR/archive/rescued_before_${SLURM_JOB_ID}_$(date +%Y%m%d_%H%M%S)"
for d in density-column eval vtk-out water-height; do
    if [[ -d "$BIN_DIR/$d" ]]; then mkdir -p "$RESCUE"; mv "$BIN_DIR/$d" "$RESCUE/$d"; fi
done
if [[ "$PRM" == cylcrash_base_gas.prm && -e "$BIN_DIR/snapshots_cylcrash" ]]; then mkdir -p "$RESCUE"; mv "$BIN_DIR/snapshots_cylcrash" "$RESCUE/"; fi
"$SH" --build --jobs "$(nproc)" --bed --bed-np 288
"$SH" --run --np 288 --prm "$PRM_DIR/$PRM" > "$OUT/$LAB.fulllog" 2>&1
echo "exit=$?" > "$OUT/$LAB.status"
grep -E "\[CylForces\]|\[CylTorque\]" "$OUT/$LAB.fulllog" | cut -c1-460 > "$OUT/$LAB.forces.txt"
gzip "$OUT/$LAB.fulllog"
echo "## DONE $OUT"
