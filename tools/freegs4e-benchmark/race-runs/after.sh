#!/bin/bash
# C shaped's MEQ arm, after halfdisc's coil-ownership fix; waits for the chain.
bench=$(cd "$(dirname "$0")/.." && pwd)
stages=${MEQ_RUNS_DIR:-$bench/../../../runs}/freegs4e-benchmark/race-runs
mkdir -p "$stages" && cd "$stages" || exit 1
for i in $(seq 1 6000); do grep -q 'CHAIN DONE' chain.status 2>/dev/null && break; sleep 60; done
"$bench/venv/bin/python" "$bench/race_inverse.py" run ir2-c-meq --config machine-c-mast-shaped --codes meq --starts cold,warm --repeats 3 > ir2-c-meq.log 2>&1
echo "ir2-c-meq exit $?" >> chain.status
echo "AFTER DONE" >> chain.status
