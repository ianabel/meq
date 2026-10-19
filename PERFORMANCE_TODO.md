# Performance runs owed

Two pieces of measurement, in this order: **profile D and F shaped first**,
then **re-run every MEQ arm the inverse-race doc reports** (F, C and D shaped).
The profiles come first because they decide whether anything changes in the
solver or in MFEM before the race numbers are taken, and the race numbers
should be taken once, on the code they will be quoted for.

The doc is the Claude Docs artifact *Inverse race: F, C and D shaped*
(`claude.ai/artifact/WAq3D1udMQakGMYyvHdL1g`). Its MEQ rows predate the stall
test, and C has no MEQ row at all.

## Before either

* **A quiet machine.** Nothing else running — the last profiles shared the box
  with an `xcql3d` run and a MaNTA build, and the k = 3 cold profile read twice
  its uncontended wall. `uptime` under 1 before starting; the race enforces
  this itself (`--load-max`, `--quiet-wait`), the profiles do not.
* **Commit what is being measured, then build and install from it**, so
  `meq --version` and every record's `provenance.tree` name a commit rather
  than `-dirty`. The working tree carries the stall test and crawl report,
  the threaded cold-start transfer and `--trace-solver-stats`, none of which
  any recorded race number has seen.
* **`OMP_WAIT_POLICY=passive` for every profile**, or the `cores` column is
  meaningless: a spinning team makes a serial leg read high
  (`docs/linear_solvers.rst`). Leave the race arms at their default, as every
  run already in the doc was.
* **`OMP_NUM_THREADS=8 MKL_NUM_THREADS=8`**, the physical core count
  (`tools/freegs4e-benchmark/cores.py`), as the race uses.

## 1. Profiled D and F shaped

**What it answers.**

| question | read from | what it decides |
| --- | --- | --- |
| Does PARDISO refine? Per step the ~20-column solve costs ~3× the factorisation | `--trace-solver-stats`: perturbed pivots per factorisation, refinement steps per solve | item 1 of `../../mfem-hdg-dev/doc/HDG-TRACE-SOLVE-AND-NPC-REDUCE-FROM-MEQ.md`: closed if it never refines |
| How threaded are the legs, honestly? | `--profile` under `passive`: `cores` per leg | whether reduce+recover (15–21%) and the backsolve (15–19%) are thread-bound; item 2 of the same request |
| Did threading the cold-start transfer pay? | `--profile`: `outside solve() → makeSolver`, cold runs | was 0.86 s (k = 2) and 1.45 s (k = 3) at 1.00 core on D; warm reads 0.008 s |
| What does the stall test cost F's cold start now? | `--profile` legs, and the log's "stopped moving" line | single runs read 14.8 s / 20.8 s / 43.5 s at k2 h1 / k3 h1 / k3 h0.5 against 80 / 128 / 307 s before |
| Does a converged answer crawl? | the converged report's crawl line; `newton_longest_crawl` in the `.nc` | any warning here is a defect report, not a pass (`TODO`, *Bordered solves that CRAWL*) |

**The runs.** Both configurations, both starts, at the rungs the doc's
1e-3 and 1e-4 rows are decided on: k2 h1.00 and k3 h1.00 for both, and k3
h0.50 for F, which is F's 1e-4 rung. That is 2 × 2 × 2 + 2 = 10 configurations.

Each configuration gets a directory from the race arm — a cold run first, then
a warm one restarted from it — and is then re-run by hand under the
instruments, so the profile is of exactly the file the race times:

```sh
cd source/tools/freegs4e-benchmark
P=../../../runs/freegs4e-benchmark/ref-interp/problem.json
OUT=../../../runs/perf-profile          # under runs/, never in source/
export OMP_NUM_THREADS=8 MKL_NUM_THREADS=8
for cfg in machine-d-tcv-shaped machine-f-diiid-shaped; do
  for r in 0 2 $( [ $cfg = machine-f-diiid-shaped ] && echo 3 ); do
    c=$OUT/$cfg/r$r/cold; w=$OUT/$cfg/r$r/warm; mkdir -p $c $w
    venv/bin/python race_arm_meq.py --problem $P --config $cfg --rung $r --start cold --out $c
    venv/bin/python race_arm_meq.py --problem $P --config $cfg --rung $r --start warm --state $c --out $w
  done
done
```

Then in each of those directories, three times each and interleaved across
the ten so drift spreads over all of them, with `$MEQ` the installed binary the
arm itself runs (`venv/bin/python -c 'import benchpaths; print(benchpaths.meq())'`):

```sh
$MEQ --mesh-ready meq.toml                                        # plain wall
OMP_WAIT_POLICY=passive $MEQ --mesh-ready --profile meq.toml      # the legs
$MEQ --mesh-ready --trace-solver-stats meq.toml                   # once is enough
```

A warm directory's `meq.toml` names its mesh and guess by absolute path into
the cold one, so re-running in place is right; copy a directory only with its
paths rewritten.

**Not measurable from MEQ, and asked upstream instead:** reduce and recover
separately. `GetNPCTraversalTime()` times them together, and the request's
item 2 asks for the split and the colour count.

**Record** the leg tables as an `M-nn` entry in `MEASUREMENTS.md`, and put the
two numbers the upstream request asks for — refinement steps, and the passive
`cores` of the backsolve and of reduce+recover — into that request as a reply
section rather than a new file.

## 2. Re-run the MEQ arms the doc reports

**All of it**: F, C and D shaped, cold and warm, three repeats, every rung of
the ladder including the k4 h0.35 datum. One configuration per invocation, as
the race protocol requires:

```sh
cd source/tools/freegs4e-benchmark
for c in f-diiid c-mast d-tcv; do
  venv/bin/python race_inverse.py run ir3-${c%%-*}-meq \
      --config machine-$c-shaped --codes meq --starts cold,warm --repeats 3 \
      > ../../../runs/freegs4e-benchmark/race-runs/ir3-${c%%-*}-meq.log 2>&1
done
```

Sequentially, not chained behind anything that can die: the last chain stopped
at D when its orchestrator went away, and the C re-run waiting on its
`CHAIN DONE` never ran.

**Then merge with the other codes' records**, which are not being re-run, and
summarise. The doc's F summary was made this way — `ir1-f-final` is `ir1-f`'s
non-MEQ records plus `ir1-f-meq`'s:

```sh
R=../../../runs/freegs4e-benchmark/race-runs
for pair in f:ir1-f-final c:ir2-c d:ir2-d-final; do
  k=${pair%%:*} old=race-runs/${pair##*:}; new=$R/ir3-$k-final; mkdir -p $new
  cp $old/plan.json $new/
  python3 -c "import json,sys; [print(l,end='') for l in open('$old/records.jsonl') if json.loads(l)['code']!='meq']" > $new/records.jsonl
  cat $R/ir3-$k-meq/records.jsonl >> $new/records.jsonl
  venv/bin/python race_inverse.py summarise $new
done
```

`ir2-c` holds C's failed MEQ runs as `meq` records with exit 1; the filter
drops them with the rest of the old MEQ rows. Copy the three `ir3-*-final`
stages into the tracked `race-runs/` beside the old ones, which stay as the
record of what the doc first said.

**What to check before quoting anything:**

* every MEQ run accepted, and on C that means exit 4 ("best fit") with the
  rows inside `--row-tol` — C's targets are rank-deficient by symmetry;
* `newton_longest_crawl` in each record: report any converged answer at 8 or
  more in the doc rather than averaging it into a median;
* MEQ against freegs4e's finest still near 1.5e-5 on F and 3.5e-6 on D. The
  stall test changes iterations, not answers, and the single F runs agreed to
  every printed digit; a moved floor means something else changed;
* F's warm rungs unchanged within noise (3–17 s): the stall test does not fire
  on a warm start.

**Then update the doc**: the time-to-accuracy chart, the who-finished table
(C's MEQ cell), the defects table (the stall test as a fix that changes timed
numbers), and the caveats (F's cold start no longer lands in a 200-iteration
stall).

## Not in this round, deliberately

* **E and G shaped**, which have never been run for any code. A separate,
  all-codes stage.
* **The exterior settings.** Moving F's `[boundary.exterior] Radius` from 3.4 m
  to 3.6 m and its modes from 10 to 20 cuts F's error 4–9× at k3–k5, but it
  changes the problem MEQ's file poses, so it is an accuracy decision with its
  own write-up, not a change to slip in under a timing re-run.
