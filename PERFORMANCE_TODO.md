# Performance runs owed

**Re-run every MEQ arm the inverse-race doc reports** (F, C and D shaped). The
profiles that were to come first are done, M-300: PARDISO refines every solve
to its cap and the backsolve and reduce+recover legs run on three to four of
eight cores, both now with `../../mfem-hdg-dev` as a request. Nothing in the
solver changed on them, so the race numbers can be taken on the code as it is.

The doc is the Claude Docs artifact *Inverse race: F, C and D shaped*
(`claude.ai/artifact/WAq3D1udMQakGMYyvHdL1g`). Its MEQ rows predate the stall
test, and C has no MEQ row at all.

## Before it

* **A quiet machine.** Nothing else running; the race enforces this itself
  (`--load-max`, `--quiet-wait`).
* **Commit what is being measured, then build and install from it**, so
  `meq --version` and every record's `provenance.tree` name a commit rather
  than `-dirty`.
* **The race arms at the default `OMP_WAIT_POLICY`**, as every run already in
  the doc was; passive is for a profile's `cores` column.
* **`OMP_NUM_THREADS=8 MKL_NUM_THREADS=8`**, the physical core count
  (`tools/freegs4e-benchmark/cores.py`), as the race uses.

## Re-run the MEQ arms the doc reports

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
