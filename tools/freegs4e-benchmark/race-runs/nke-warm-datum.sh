#!/bin/bash
# Waits for D shaped's orchestrator (26758), then re-runs the warm freegsnke datum.
bench=$(cd "$(dirname "$0")/.." && pwd)
stages=${MEQ_RUNS_DIR:-$bench/../../../runs}/freegs4e-benchmark/race-runs
mkdir -p "$stages" && cd "$stages" || exit 1
while kill -0 26758 2>/dev/null; do sleep 30; done
echo "ir2-d exit (orchestrator gone)" >> chain.status
"$bench/venv/bin/python" "$bench/race-runs/nke-warm-datum.py" > ir2-d-nke-warm.log 2>&1
echo "ir2-d-nke-warm exit $?" >> chain.status
