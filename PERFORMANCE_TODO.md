# Performance runs owed

**F shaped still fails its first bordered Newton on every cold rung** and is
solved by the Picard-then-Newton retry, to the same currents. The cause is
measured ([M-315](MEASUREMENTS.md#m-315)): the support is frozen at the guess's
X-point flux, 0.0914, while the `psi_bnd` unknown starts at zero and the first
step carries it to −0.188, so `F` jumps at the frozen edge and the Jacobian
cannot see it. Starting `psi_bnd` at the edge cures F and is a branch control;
so are three ways of making the source continuous. None is landed, and the
candidate left is a continuous source as a rung reached only on failure.

**C's k2 h1.00 rung exits 4**, a best fit at a target residual of 4.7e-09, with
the fix and without it: the inverse best fit's stationarity stop (`TODO`, the
single-precision round's item 3), not a branch question.

## The timing runs owed

None. The inverse race's MEQ arms on F, C and D shaped are
[M-343](MEASUREMENTS.md#m-343) and [M-344](MEASUREMENTS.md#m-344), and IH-0's
five repeats are [M-314](MEASUREMENTS.md#m-314).

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

## Not in this round, deliberately

* **Copying `ir4-{f,c,d}-final`** from `runs/freegs4e-benchmark/race-runs/`
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
