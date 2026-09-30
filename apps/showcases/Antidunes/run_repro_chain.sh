#!/bin/bash
set -uo pipefail
BASE=/home/dhuber/University/Thesis/mt-dennys/apps/showcases/Antidunes
BIN=/home/dhuber/University/Thesis/mt-dennys/build-release/apps/showcases/Antidunes/Antidunes
for i in 1 2 3 4 5; do
  echo "=== [$(date '+%F %T')] starting run $i ==="
  cd "$BASE/rod_pair_repro_run${i}"
  rm -rf vtk-out eval snapshots
  mpirun -np 9 "$BIN" "$BASE/prms/rod_pair_4255679_repro_gas.prm" > run.log 2>&1
  RC=$?
  echo "=== [$(date '+%F %T')] run $i finished, exit=$RC ==="
  cd "$BASE"
done
echo "=== ALL REPRO RUNS DONE ==="
