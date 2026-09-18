# The coils out of the mesh: `psi = psi_coils_analytic + psi_solved`

Ian's idea, 2026-09-18, in his words: *"possibly we could have inside-coil
approximation that uses the fact that the near field of a rectangular coil is
solvable? approximate psi_coil as psi_analytic_coil + polynomial"*, against the
requirement that *"we should be faster at 1e-3 and 1e-4 not just at very high
accuracy, high order should allow us to not use very many elements"*.

`CLAUDE.md` is the index and `CLAUDE_FB.md` owns the free boundary; this defers
to both. It is a DESIGN, not a measurement — nothing here has been built.

---

## 0. What would falsify it, and the cheapest experiment

Put first, because this plan has one load-bearing assumption and it is testable
in an afternoon without writing any of it.

**The assumption**: that MEQ's error on a machine case is limited by the coils —
either by their REGULARITY (the source jumps across a rectangle's edge, so `psi`
has limited smoothness there and a degree-`k` polynomial cannot pay) or by their
MODEL (MEQ's uniform rectangles against freegs4e's discrete filaments).

**M-111 already says it is one of those two and does not say which.** F DIII-D's
relative error is **flat at 5.785e-03 across a 16x range in dofs**, off-coil flat
at 6.3e-04 to 7.6e-04, and M-111's own verdict is *"that is the conductor model
... it is not a discretisation error and no rung buys it down"*.

**The experiment that separates them**: compute the reference's `psi` twice on
the same grid — once from freegs4e's filaments, once from the same currents
spread uniformly over MEQ's rectangles — and difference them, with no MEQ solve
anywhere. Both are Green's-function sums and neither involves a mesh.

* If that difference is about 5.785e-03, the error is the MODEL and this plan's
  value is that it lets MEQ adopt freegs4e's model exactly. Regularity is a
  bystander.
* If it is far smaller, the error is REGULARITY and the subtraction is the
  repair.
* If it is far larger, something else is wrong and this plan is premature.

**Do this before building anything.** It is one script, it needs no MEQ, and it
decides which half of the argument below is the real one.

---

## 1. The split

Write the total flux as

```
psi = psi_c + psi_p
```

where `psi_c` is the field of the CONDUCTORS alone, computed analytically, and
`psi_p` is what MEQ solves for. Since `Delta*` is linear and the conductor
currents are prescribed inputs of a forward solve,

```
Delta* psi_p = -mu0 r J_plasma( psi_c + psi_p )
```

with the coil term gone from the right-hand side entirely. **The coils need not
be in the mesh at all.**

`psi_c` is a Green's-function integral over each conductor. It needs no closed
form — a closed form for a rectangle is not required and probably not worth
chasing — because the integrand is smooth except for a log at the source point
and a tensor-product Gauss rule reaches machine precision on it cheaply. MEQ
already owns that quadrature in `meq::Coil`.

**AND IT IS COMPUTED ONCE PER MESH, NOT PER NEWTON STEP.** The conductor
currents do not move during a forward solve; only the profile scale does. So
`psi_c` at every quadrature point and every boundary node is a precompute, and
its cost is amortised over every iteration of every support sweep.

## 2. What it buys, in the order the evidence supports

**(a) The conductor model stops being a discretisation question.** This is the
strongest claim and the one M-111 most directly supports. With `psi_c` computed
by quadrature, MEQ can use whatever conductor model the reference uses —
filaments at freegs4e's own positions, or exact rectangles — and the 5.785e-03
that no rung buys down becomes a modelling CHOICE rather than an error floor.
A benchmark that currently cannot resolve below 1.4e-04 because the two codes
model conductors differently becomes one where they model them identically.

**(b) Elements.** On `examples/mastu-nke.msh`: **1024 of 9361 triangles, 10.9%,
carry a coil attribute** — and that undersells it, because `CoilSize = 0.03`
against `Size = 0.30` forces a 10x graded refinement AROUND each coil as well.
The honest estimate is a fifth to a third of the mesh, not an order of
magnitude. **Do not sell this plan on the element count.**

**(c) Regularity, IF §0's experiment says so.** A jump in the source across a
rectangle's edge limits `psi`'s smoothness there, and a rectangle CORNER limits
it further; neither is something degree `k` recovers. Remove the coils from the
solved field and `psi_p`'s only irregularity is the plasma edge, which is
`CLAUDE_FB.md`'s `k <= j` result and is a separate known cap.

## 3. What it costs, and the two that are not obvious

* **The Dirichlet datum and the exterior coupling both move.** On `Gamma` the
  condition becomes `psi_p = psi_datum - psi_c`, and FB-5's DtN coupling then
  sees the plasma alone. That is arguably SIMPLER — the exterior problem's
  source is one connected region instead of a plasma plus twenty-three
  conductors — but it is a real change to `meq::ExteriorDtN`'s setup and to
  `setBoundaryData()`, and the two must move together or the datum is
  double-counted.
* **Every place that reads `psi` must read `psi_c + psi_p`.** The critical-point
  finder, the plasma-support fill, the limiter and X-point borders, the output
  writers, the flux-surface tracer. **This is the real cost of the plan**: the
  solved field stops being the physical field, and every consumer that forgets
  is a silent wrong answer rather than a failure. `CLAUDE.md`'s own record has
  three separate defects of exactly that shape — `psi*` against `psi_h`, the
  band continuation reading the foot, `scaledF` against `f`.
* **FB-7 already does a piece of this** — a conductor outside `Gamma` that the
  mesh does not carry — so the machinery for "a coil MEQ does not mesh" exists
  and this generalises it inward rather than inventing it.

## 4. What it does NOT fix

The plasma edge. `ConstrainPaxisIp` at `alpha_n = 1.2` makes the edge a
fractional power and caps the order at `k <= 1.2` whatever the coils do, so
**MAST-U cannot demonstrate this plan's benefit** and a smooth-edge case is
needed to test it. Nor does it touch anything in M-115, M-117 or M-119 — the
constraint discontinuities and the cold failures are about the border, not the
conductors.

## 5. Staging

| | |
|---|---|
| **CS-0** | §0's two-Green's-function difference. Decides model against regularity, needs no MEQ, and can kill or redirect the plan |
| **CS-1** | `meq::ConductorField`: `psi_c` and `grad psi_c` by quadrature, against an analytic single-loop check |
| **CS-2** | the split on a FIXED-boundary case with coils, where nothing else moves |
| **CS-3** | the Dirichlet datum and the DtN coupling |
| **CS-4** | every consumer of `psi`, with a test per consumer that the total is read |
| **CS-5** | re-take M-111 with the conductor models matched |
