#!/bin/bash
# IR-1's MEQ arm re-run on c57ae78, then IR-2, one configuration per invocation.
# Stages, logs and chain.status go to runs/freegs4e-benchmark/race-runs/.
bench=$(cd "$(dirname "$0")/.." && pwd)
stages=${MEQ_RUNS_DIR:-$bench/../../../runs}/freegs4e-benchmark/race-runs
mkdir -p "$stages" && cd "$stages" || exit 1
"$bench/venv/bin/python" "$bench/race_inverse.py" run ir1-f-meq --config machine-f-diiid-shaped --codes meq --starts cold,warm --repeats 3 > ir1-f-meq.log 2>&1
echo "ir1-f-meq exit $?" >> chain.status
for c in c-mast d-tcv e-diamagnetic g-mastu; do
  "$bench/venv/bin/python" "$bench/race_inverse.py" run ir2-${c%%-*} --config machine-$c-shaped --starts cold,warm --repeats 3 > ir2-${c%%-*}.log 2>&1
  echo "ir2-${c%%-*} exit $?" >> chain.status
done
echo CHAIN DONE >> chain.status
