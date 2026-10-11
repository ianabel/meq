"""D shaped's warm freegsnke datum, re-run on race_arm_nke.py's DATUM_GS_TOL
after the 1e-10 run was stopped; one record into ir2-d-nke-warm/."""
import argparse, json, os, sys
HERE = "/home/ian/projects/meq/tools/freegs4e-benchmark"
sys.path.insert(0, HERE)
import race_inverse as ri
out = os.path.join(HERE, "race-runs/ir2-d-nke-warm")
os.makedirs(out, exist_ok=True)
plan = json.load(open(os.path.join(HERE, "race-runs/ir2-d/plan.json")))
args = argparse.Namespace(load_max=plan["load_max"], quiet_wait=plan["quiet_wait"],
                          row_tol=plan["row_tol"])
state = os.path.join(HERE, "race-runs/ir2-d/machine-d-tcv-shaped/cold/freegsnke/r02/rep0/state")
with open(os.path.join(out, "records.jsonl"), "a") as log:
    ri.run_arm("freegsnke", plan["problem"], plan["config"], 2, "warm",
               os.path.join(out, "machine-d-tcv-shaped/warm/freegsnke/r02/rep0"),
               state, args, log)
