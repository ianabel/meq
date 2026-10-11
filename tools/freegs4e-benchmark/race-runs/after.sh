#!/bin/bash
# C shaped's MEQ arm, after halfdisc's coil-ownership fix; waits for the chain.
cd /home/ian/projects/meq/tools/freegs4e-benchmark
for i in $(seq 1 6000); do grep -q 'CHAIN DONE' race-runs/chain.status && break; sleep 60; done
./venv/bin/python race_inverse.py run race-runs/ir2-c-meq --config machine-c-mast-shaped --codes meq --starts cold,warm --repeats 3 > race-runs/ir2-c-meq.log 2>&1
echo "ir2-c-meq exit $?" >> race-runs/chain.status
echo "AFTER DONE" >> race-runs/chain.status
