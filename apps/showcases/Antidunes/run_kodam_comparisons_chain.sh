#!/bin/bash
set -uo pipefail
BASE=/home/dhuber/University/Thesis/mt-dennys/apps/showcases/Antidunes
BIN=/home/dhuber/University/Thesis/mt-dennys/build-release/apps/showcases/Antidunes/DropTest
DATE=$(date +%Y%m%d)

run_one() {
  local name="$1" np="$2" prm="$3" outdir="$4" extra="$5"
  echo "=== [$(date '+%F %T')] starting $name (np=$np) ==="
  cd "$BASE/$outdir" || exit 1
  mpirun -np "$np" $extra "$BIN" "$BASE/prms/$prm" > run.log 2>&1
  echo "=== [$(date '+%F %T')] finished $name, exit=$? ==="
  cd "$BASE"
}

run_one "SAT-cuboid (250)"        36 DropTest_250_sat_gas.prm                        "DropTest_250_sat_output_${DATE}"               "--oversubscribe"
run_one "Kodam packing (closed-form)" 9 DropTest_kodam_packing_gas.prm               "DropTest_kodam_packing_output_${DATE}"         ""
run_one "Kodam packing (GJK-only)"    9 DropTest_kodam_packing_gjk_gas.prm           "DropTest_kodam_packing_gjk_output_${DATE}"     ""
run_one "Kodam packing (pure-analytical)" 9 DropTest_kodam_packing_pureanalytical_gas.prm "DropTest_kodam_packing_pureanalytical_output_${DATE}" ""
echo "=== ALL DONE ==="
