# The inverse race: the arm contract

`INVERSE-RACE-PLAN.md` is the design; this is the interface every arm
implements so `race_inverse.py` can drive four codes the same way. An arm is a
command. The orchestrator runs it once per (configuration, rung, start,
repeat), in a fresh directory, with the environment it chooses (threads, MKL),
and reads one JSON record back. **An arm never decides whether it won**: it
reports what it computed and the orchestrator applies every acceptance test.

## The problem file

`inverse_problem.py` writes it from a reference set
(`runs/freegs4e-benchmark/ref-interp/` by default: `fgsref.py --shaped
--interpolating-profiles`, so every code solves the analytic profile). It is a
list of configurations in the format M-254's `problem.json`
introduced, plus:

| field | what |
|---|---|
| `vacuum_seed` | `{ supply: A }` for the FREE supplies: the minimum-norm currents meeting the targets against the design blob and the fixed supplies' own field — `freegs4e`'s first control step. **The cold start's free currents** |
| `design` | the reference's `design` block: `R0`, `Ip`, the design X-points and isoflux pairs — the blob a cold start represents, each code in its own way |
| `shifted_xpoint_targets`, `shifted_isoflux_targets` | **the warm start's targets**: every target point moved 0.01 m along the unit vector from the reference's magnetic axis to it |
| `meq_example` | the shipped `examples/machine-*-shaped.toml` the MEQ arm derives from |
| `symmetric`, `rank`, `rows` | the targets are up–down symmetric; the vacuum response's rank and the row count. **Where `rank < rows` the row set is consistent only by that symmetry** (C, E, G), so an arm that can make its discretisation symmetric should: an asymmetric one cannot meet all the rows and reports a best fit (M-255, E) |

Everything else — `conductors`, `supplies` with `free` and `reference_current`,
`xpoint_targets`, `isoflux_targets`, `Ip`, `profile_shape`, the two amplitudes,
`reference` — is as before.

## The command

```
<arm> --problem P.json --config NAME --rung R --start START [--state DIR] --out DIR
<arm> --problem P.json --config NAME --list-rungs
```

* `--rung R` indexes the arm's OWN ladder, 0 the coarsest. `--list-rungs`
  prints a JSON list, one object per rung, `{"rung": R, "label": "...",
  "finest": bool}` — the last rung is the accuracy datum, run but not timed
  as a contestant.
* `--start`:
  * `cold` — fixed supplies at `reference_current`, free supplies at
    `vacuum_seed`, the plasma as the code's ordinary representation of the
    design blob carrying `Ip`, the DESIGN targets. Nothing solved for.
  * `seeded` — the same, free supplies at `reference_current`. Untimed: the
    orchestrator uses it only to obtain a converged state when `cold` does
    not finish.
  * `warm` — restart from `--state DIR` (this arm's own converged design-target
    run at the same rung): its field, its plasma, its currents; the SHIFTED
    targets. This is the timed shape step.
* `--out DIR` exists and is empty. The arm writes `DIR/record.json`, and
  anything it needs for a later `warm` start under `DIR` (its `state`).

## The record

```json
{
  "code": "meq", "config": "machine-f-diiid-shaped", "rung": 2,
  "rung_label": "k2 h1.0", "unknowns": 12345, "start": "cold",
  "finished": true, "status": "converged",
  "wall_s": 12.3, "solve_s": 8.1,
  "currents": {"F4A": 1.7e5, "...": 0.0},
  "psi_axis": 0.3758, "psi_bndry": 0.0708,
  "axis": [1.768, -0.0002], "xpoints": [[1.20, -1.00], [1.20, 1.00]],
  "rows": [1e-12, 3e-12, 2e-13, 4e-12, 1e-14, 2e-14],
  "iterations": 4, "state": "DIR/state",
  "provenance": {"code_version": "...", "tree": "sha", "notes": "..."}
}
```

* `wall_s` is **end to end**, from the arm's start to its record — meshing,
  boundary matrices, the solve, its own outputs. `solve_s` is the solve alone
  where the code can say, else `null`.
* `currents` carries every FREE supply in amperes, in the problem's supply
  convention: a supply's current, member conductor `c` carrying
  `multiplier × I`.
* `rows` is the residual of each target row in the problem's order — X-point
  rows as `B_R, B_Z` per point in T, then one `psi1 − psi2` per isoflux pair in
  Wb/rad — **evaluated from the arm's converged FIELD, not taken from the
  solver's stopping test**. A solver that reports success with its rows at
  0.1 T (freegsnke's step ratchet, M-255) is caught here and nowhere else.
  For `warm` the rows are against the SHIFTED targets.
* `finished` false with a `status` saying why — a DNF is a row, not a gap.
* `unknowns` is the code's own measure of size (MEQ's dofs, NICE's P1 nodes,
  the grid points), for the Pareto plot's other axis.

## The rules an arm keeps

* No other code's answer, ever, beyond what the problem file holds.
* Its own ladder, from coarse to the finest it can afford in memory;
  `INVERSE-RACE-PLAN.md` §4 lists what each code reached in M-255.
* Threads from the environment it is given; it does not set its own.
* It never edits another code's source tree; patches live in the harness and
  are recorded in `provenance.notes`.
