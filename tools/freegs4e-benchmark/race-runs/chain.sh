#!/bin/bash
# IR-1's MEQ arm re-run on c57ae78, then IR-2, one configuration per invocation.
cd /home/ian/projects/meq/tools/freegs4e-benchmark
./venv/bin/python race_inverse.py run race-runs/ir1-f-meq --config machine-f-diiid-shaped --codes meq --starts cold,warm --repeats 3 > race-runs/ir1-f-meq.log 2>&1
echo "ir1-f-meq exit $?" >> race-runs/chain.status
for c in c-mast d-tcv e-diamagnetic g-mastu; do
  ./venv/bin/python race_inverse.py run race-runs/ir2-${c%%-*} --config machine-$c-shaped --starts cold,warm --repeats 3 > race-runs/ir2-${c%%-*}.log 2>&1
  echo "ir2-${c%%-*} exit $?" >> race-runs/chain.status
done
echo CHAIN DONE >> race-runs/chain.status
