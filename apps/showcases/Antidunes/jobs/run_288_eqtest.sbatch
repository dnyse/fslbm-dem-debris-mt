#!/bin/bash -l
#SBATCH --job-name=eqtest
#SBATCH --partition=multinode
#SBATCH --nodes=4
#SBATCH --ntasks=288
#SBATCH --time=01:15:00
#SBATCH --output=sbatch-output/output/eqtest_%j.out
#SBATCH --error=sbatch-output/error/eqtest_%j.err

# [Claude, 2026-10-08] short diagnostic runs of the equilibrium set-up (prms/eqtest_*_gas.prm), e.g.
#   sbatch --job-name=eqtest-deep025-walls-oldseed --export=ALL,EQ_PRM=eqtest_deep025_walls_oldseed_gas.prm jobs/run_288_eqtest.sbatch
# Results (status, forces of all floater types, fallback count, eval sqlite, full log) in pinscan_results/<jobid>/; no VTK kept.
# Configures a fresh checkout first. Exits with the simulation's exit code.
set -uo pipefail
PRM="${EQ_PRM:?set EQ_PRM}"; LAB="${PRM%_gas.prm}"
REPO="$(cd "$SLURM_SUBMIT_DIR/../../.." && pwd)"
BIN_DIR="$REPO/build/apps/showcases/Antidunes"; SH="$REPO/apps/showcases/Antidunes/antidunes.sh"; PRM_DIR="$REPO/apps/showcases/Antidunes/prms"
OUT="$REPO/apps/showcases/Antidunes/pinscan_results/${SLURM_JOB_ID}"; mkdir -p "$(dirname "$OUT")"; mkdir "$OUT" || exit 1
CONF=""; [[ -f "$REPO/build/CMakeCache.txt" ]] || CONF="--configure"
"$SH" $CONF --build --jobs "$(nproc)" --bed --bed-np 288
"$SH" --run --np 288 --prm "$PRM_DIR/$PRM" > "$OUT/$LAB.fulllog" 2>&1
SIMEXIT=$?
echo "exit=$SIMEXIT" > "$OUT/$LAB.status"
grep -E "\[CylForces\]|\[CylTorque\]|\[BoxForces\]|\[BoxTorque\]|\[FloaterForces\]" "$OUT/$LAB.fulllog" | cut -c1-460 > "$OUT/$LAB.forces.txt"
grep -c "FALLBACK-ISOLATED-INTERFACE" "$OUT/$LAB.fulllog" | sed 's/^/fallback_isolated_interface=/' >> "$OUT/$LAB.status"
cp "$BIN_DIR/eval/antidunes.sqlite" "$OUT/$LAB.antidunes.sqlite" 2>/dev/null
gzip "$OUT/$LAB.fulllog"
rm -rf "$BIN_DIR/density-column" "$BIN_DIR/eval" "$BIN_DIR/vtk-out" "$BIN_DIR/water-height"
echo "## DONE $OUT"
exit "$SIMEXIT"
