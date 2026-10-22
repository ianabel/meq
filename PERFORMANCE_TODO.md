# Performance runs owed

**MEQ on C shaped (`machine-c-mast-shaped`) converges on two of its six
rungs**, k2 h1.00 and k4 h0.50, cold and warm. With no accepted finest rung, it
has no self-error on C. Every cold run on C, and on F, fails the first bordered
Newton and is solved by the Picard-then-Newton retry; the two defects below are
where that retry fails too. Measured at `0e7177a`, the race's `ir3-c-meq`
stage; the doc is the Claude Docs artifact *Inverse race: F, C and D shaped*
(`claude.ai/artifact/WAq3D1udMQakGMYyvHdL1g`).

Each failed run's inputs are in
`runs/freegs4e-benchmark/race-runs/ir3-c-meq/machine-c-mast-shaped/cold/meq/r0N/rep0/`
(`meq.toml`, `meq.msh`, `guess.gf`, `meq.log`); `rN` is the rung, 0 to 5 in the
ladder's order.

## 1. The middle rungs stall, and the retry trips the X-point relocation cap

k2 h0.50 (`r01`), k3 h1.00 (`r02`) and k3 h0.50 (`r03`), all three repeats and
the seeded start: the bordered Newton exits "has stopped moving: 8 steps in a
row took less than 1/64 of the Newton step". The Picard-then-Newton retry then
finds no damping with a finite residual, because every trial moves the X-point
past the 0.1225 m relocation cap: 0.123573 m (`r01`), 0.122586 m (`r02`),
0.124438 m (`r03`). The rungs either side, `r00` and `r04`, also hit the cap
(3 and 5 line-search trials rejected) but converge.

## 2. The datum crawls to the wrong equilibrium

k4 h0.35 (`r05`): the retry converges after the bordered Newton CRAWLED 40
steps in a row below 1/64 of the Newton step, with 559 line-search trials
rejected by the cap. The answer is a different equilibrium: psi_ax 0.2769 Wb/rad
against freegs4e's finest 0.0901, P2 and P3 of opposite sign to freegs4e's,
and the plasma containing the symmetry axis. The driver refuses it (exit 1,
"THE PLASMA CONTAINS THE SYMMETRY AXIS"), which is right; reaching it is the
defect. Its seeded start exits 2.

## Before any timing that follows a fix

* **A quiet machine.** The race enforces this itself (`--load-max`,
  `--quiet-wait`); `runs/inverse-race/queue.sh` adds a five-minute
  `/proc/stat` gate before each configuration.
* **Commit what is being measured, then build and install from it**, so
  `meq --version` and every record's `provenance.tree` name a commit rather
  than `-dirty`.
* **The race arms at the default `OMP_WAIT_POLICY`**, with
  `OMP_NUM_THREADS=8 MKL_NUM_THREADS=8`, the physical core count
  (`tools/freegs4e-benchmark/cores.py`).
* **Re-run C's MEQ arm whole**, cold and warm, three repeats, every rung, and
  merge it with `ir2-c`'s other codes' records as `queue.sh` does. Check every
  record's `newton_longest_crawl`, and MEQ against freegs4e's finest: k4 h0.50
  is 3.7e-6 from it at `0e7177a`.

## Not in this round, deliberately

* **Copying `ir3-{f,c,d}-final`** from `runs/freegs4e-benchmark/race-runs/`
  into the tracked `tools/freegs4e-benchmark/race-runs/`, beside the stages
  the doc first quoted. Not decided.
* **Why F's cold start got about 10× cheaper** between `c57ae78` and `0e7177a`
  (k3 h0.50: 307 s to 31 s, same 2 iterations and 3 sweeps). Not measured.
* **E and G shaped**, which have never been run for any code. A separate,
  all-codes stage.
* **The exterior settings.** Moving F's `[boundary.exterior] Radius` from 3.4 m
  to 3.6 m and its modes from 10 to 20 cuts F's error 4–9× at k3–k5, but it
  changes the problem MEQ's file poses, so it is an accuracy decision with its
  own write-up, not a change to slip in under a timing re-run.
