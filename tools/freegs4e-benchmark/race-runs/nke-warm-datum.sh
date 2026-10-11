#!/bin/bash
# Waits for D shaped's orchestrator (26758), then re-runs the warm freegsnke datum.
cd /home/ian/projects/meq/tools/freegs4e-benchmark
while kill -0 26758 2>/dev/null; do sleep 30; done
echo "ir2-d exit (orchestrator gone)" >> race-runs/chain.status
./venv/bin/python race-runs/nke-warm-datum.py > race-runs/ir2-d-nke-warm.log 2>&1
echo "ir2-d-nke-warm exit $?" >> race-runs/chain.status
