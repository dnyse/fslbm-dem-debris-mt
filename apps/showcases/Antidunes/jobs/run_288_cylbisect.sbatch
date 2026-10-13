#!/bin/bash -l
#SBATCH --job-name=cylbisect
#SBATCH --partition=multinode
#SBATCH --nodes=4
#SBATCH --ntasks=288
#SBATCH --time=01:30:00
#SBATCH --output=sbatch-output/output/cylbisect_%j.out
#SBATCH --error=sbatch-output/error/cylbisect_%j.err

# [Claude, 2026-10-08] crash bisect of the r075 cylinder (NaN at step 8750 with code a4b3b50b). Runs prms/cylbisect_r075_gas.prm
# (thesis prm of 4278476, 9000 steps) in the checkout it is submitted from: ~/mt-dennys-bisectA (code 5683c85f) or
# ~/mt-dennys-bisectB (a4b3b50b + legacy seeding). Fresh build dir, hence 1:30 h. Exits with the simulation's exit code.
#   cd ~/mt-dennys-bisectA/apps/showcases/Antidunes && sbatch --job-name=cylbisect-A jobs/run_288_cylbisect.sbatch
set -uo pipefail
PRM=cylbisect_r075_gas.prm; LAB="${PRM%_gas.prm}"
REPO="$(cd "$SLURM_SUBMIT_DIR/../../.." && pwd)"
BIN_DIR="$REPO/build/apps/showcases/Antidunes"; SH="$REPO/apps/showcases/Antidunes/antidunes.sh"; PRM_DIR="$REPO/apps/showcases/Antidunes/prms"
OUT="$REPO/apps/showcases/Antidunes/pinscan_results/${SLURM_JOB_ID}"; mkdir -p "$(dirname "$OUT")"; mkdir "$OUT" || exit 1
"$SH" --configure --build --jobs "$(nproc)" --bed --bed-np 288   # [Claude, 2026-10-08] fresh checkout: configure first (Release, codegen, as the other checkouts)
"$SH" --run --np 288 --prm "$PRM_DIR/$PRM" > "$OUT/$LAB.fulllog" 2>&1
SIMEXIT=$?
echo "exit=$SIMEXIT" > "$OUT/$LAB.status"
grep -E "\[CylForces\]|\[CylTorque\]" "$OUT/$LAB.fulllog" | cut -c1-460 > "$OUT/$LAB.forces.txt"
grep -c "FALLBACK-ISOLATED-INTERFACE" "$OUT/$LAB.fulllog" | sed 's/^/fallback_isolated_interface=/' >> "$OUT/$LAB.status"
gzip "$OUT/$LAB.fulllog"
echo "## DONE $OUT"
exit "$SIMEXIT"
