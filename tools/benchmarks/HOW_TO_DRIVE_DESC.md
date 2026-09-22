# How to drive DESC

DESC is an inverse spectral equilibrium code in JAX — it solves for the flux
surfaces rather than for `psi(R, z)`, and its profiles are labelled by
**toroidal** flux. It is MEQ's second independent comparison.

**Harness:** `tools/desc-benchmark/`  ·  **Interpreter:**
`tools/desc-benchmark/venv/bin/python` (**not** the freegs4e venv — that one has
no `desc` and fails with `ModuleNotFoundError: No module named 'desc'`).

## What works, and what to spend your time on

**THE FIXED-BOUNDARY LADDER WORKS AND MATCHES.** This is the arm to use.

```sh
tools/desc-benchmark/venv/bin/python tools/desc-benchmark/race_desc.py fixed-h-circular
```

At matched accuracy MEQ `k = 2 r = 0` reads 9.7266e-05 against DESC `M = 12` at
9.3509e-05 — **MEQ 4.80× faster cold, DESC 11.9× faster warm**. Both clocks
matter: the cold one is what solving an equilibrium once costs and the warm one
is what a parameter scan pays. The mutual floor is 1.1815e-04 and the harness
refuses targets below it. M-157, M-163 section 3.

**THE FREE-BOUNDARY ARM MOSTLY DOES NOT CONVERGE. DO NOT BUILD ON IT.**

```sh
tools/desc-benchmark/venv/bin/python tools/desc-benchmark/descfreeb.py \
    -M 12 --boundary circle
```

Of 11 runs taken with a convergence check in place, **4 converged and 7 did
not**, and the ones that converge land **0.7% to 4.6%** from the reference —
nowhere near the 1e-04 the fixed ladder reaches. Every small number the first
version of M-162 published came from a run that had **not moved from its seed**.
M-162.

## Things that will cost you an afternoon

**READ `result["success"]`.** `scipy`'s least-squares result carries the final
iterate whether or not the optimiser got anywhere, so a failed run hands back
its own starting point and looks like an answer. Seed from the reference's LCFS
and a failed run reports a tiny error *because it did not move*. `descfreeb.py`
now prints the exit message and refuses to call the result an answer; this was
the single most expensive mistake in the campaign.

**`A bad approximation caused failure to predict improvement.`** is the message
you will see most. Starting AT a stationary point — a reference seed — is the
usual cause: the trust region cannot predict improvement from a point where the
gradient already vanishes.

**DESC's DEFAULT `ftol` IS 1e-2**, which is very loose. `--desc-defaults` passes
none of the three tolerances so DESC uses its own; that is the control for a
trust-region collapse, not the setting to publish from.

**THE TOROIDAL FIELD IS NOT OPTIONAL AND ITS ABSENCE IS NOT AN ERROR.**
`BoundaryError`'s second residual is `B_out² − B_in² − 2 μ₀ p = 0`, and `B_in`
carries `B_phi = g/R` while poloidal field coils produce none. Leave out the
`ToroidalMagneticField` and the optimiser trades a good residual against an
impossible one and walks the boundary **0.81 m on a plasma of minor radius
0.34**. `g_at_gamma` is the scalar that supplies it.

**AN LCFS WITH A CORNER CANNOT BE REPRESENTED.** `FourierRZToroidalSurface` is a
truncated Fourier series and cannot turn an X-point corner, so **only a LIMITED
case can be posed free-boundary at all**. Every other free-boundary machine in
`examples/` is diverted. The same fact is why the `examples/fixed-*.toml` ladder
is posed at `psi_n = 0.95`.

**GIVE IT THREADS WITH `taskset`, NOT WITH ENVIRONMENT VARIABLES.** DESC's CPU
work goes through XLA's own Eigen thread pool, which reads neither
`OMP_NUM_THREADS` nor `MKL_NUM_THREADS`: a solve left alone was measured at
**1257% CPU** — twelve and a half cores — with `OMP_NUM_THREADS=8` set.
`race_desc.py` pins with `taskset` to one logical cpu per physical core, which
is a cap the process cannot talk its way out of, and gives MEQ the same set.

**THE PROFILES ARE LABELLED BY TOROIDAL FLUX, WHICH IS NOT KNOWN UNTIL THE
EQUILIBRIUM IS.** So `race_desc.py` poses each resolution self-consistently
**untimed** and clocks only the one forward cold-start solve that follows.

## Three harness defects worth not re-introducing

All three were the same shape — **the harness recorded something other than what
ran** — and all three are fixed:

| | recorded | actual |
|---|---|---|
| no convergence check | the final iterate, unlabelled | `success` was never read |
| `boundary_start` | the `--boundary` flag | `--blend` overrides it, so `blend = 0.000` is the CIRCLE and was stored as `"reference"` |
| `boundary_moved` on a circle seed | distance from `+a sin θ` | the seed surface is `−a sin θ`, so it measured against the circle's own **reflection** — every circle run read about `2a = 0.68 m` |

The seed **surfaces** are fine and were checked rather than assumed:
`--boundary circle` builds an exact two-mode circle, `--blend 0.0` fits the same
curve to `M` modes, and the two agree to **9.2e-16 m**. It was only the
diagnostic tuple returned beside the surface that carried the reflection.

**AND THE `.npz` NOW RECORDS `ftol`, `xtol`, `gtol` AND `maxiter`.** Without
them an archived file from the tolerance study cannot be told from one from the
seed study, and neither can be compared with a run taken later — which is what
made the first free-boundary sweep unrescuable by re-reading it, and forced it
to be re-run rather than re-scored.
