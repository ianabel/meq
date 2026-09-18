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

**MEASURED, 2026-09-18, and the answer is split.**
`tools/freegs4e-benchmark/conductor_model.py` builds both fields from the same
18 DIII-D currents on the reference's own grid — freegs4e's point filaments
against MEQ's uniform rectangles, 24² Gauss per rectangle — with no solver and
no mesh anywhere:

| region | model difference | M-111's irreducible error | |
|---|---|---|---|
| everywhere | **6.081e-03** | **5.785e-03** | **1.05×** |
| off the conductors | 2.315e-04 | 6.3e-04 to 7.6e-04 | 0.33× |
| inside them | 4.143e-02 | — | upper bound, log kernel |

**THE GLOBAL IRREDUCIBLE ERROR IS THE CONDUCTOR MODEL.** A figure that refused
to refine across a 16× range in dofs is reproduced to **5 per cent** by two
Green's-function sums. M-111 inferred that from the error refusing to refine;
this measures it.

**OFF the conductors it explains only a third**, so about two thirds of that
floor is something else — regularity, the plasma edge, or the reference's own
accuracy. *The first version of the verdict in that script read only the
off-conductor row and concluded the opposite; both rows have to be read.*

---

## 0b. THE OBJECTIVE IS ACCELERATION, NOT A CHANGE OF ANSWER

Ian, on reading the above: *"the objective of the coil subtraction is to
accelerate, not to give up solving for the exact psi including the plasma
effect"*.

That rules out the reading this plan first leaned on — "adopt freegs4e's
filament model and the floor goes away". **`psi_c` must be MEQ's own conductor
model evaluated exactly**, not somebody else's; the split is a change of
REPRESENTATION and must reproduce the same equilibrium. The benchmark agreement
that follows is a side effect of both codes then resolving their own conductors
exactly, not of MEQ adopting freegs4e's.

## 0c. AND THE VACUUM IS LINEAR, WHICH IS THE BIGGER HALF

Ian again: *"as we're outside the plasma boundary, that's linear, so we have a
way to cheat"*.

Outside the plasma `Delta* psi = 0`. The vacuum region carries no source, so
nothing there needs a non-linear solve and nothing there needs the resolution
the plasma does — it needs only to carry a harmonic extension. MEQ already
exploits this OUTSIDE `Gamma`, which is FB-5's exterior DtN coupling.

**The two ideas compound, and that is the real acceleration argument.**
`Gamma` currently has to sit outside every conductor, because the conductors are
meshed: `examples/mastu-nke.toml` needs a disc of radius **3.4** for a machine
whose plasma reaches 1.8. Take the conductors out of the mesh and `Gamma` can
move inward to just outside the plasma's reach — a half-disc of radius 2.1
instead of 3.4 is **2.6× less area** before a single element is saved on the
coils themselves or on the 10× grading around them.

So the ordering of the wins is the opposite of what §2 first said: the element
count is the point, and it comes from moving `Gamma` in rather than from
deleting 1024 coil triangles.

**What still has to be meshed**: the plasma AND the vacuum between its moving
edge and `Gamma`, because the free boundary moves during the solve and a DtN map
on a moving surface would have to be rebuilt every iteration. `Gamma` is fixed
and always contains the plasma, which is what makes its DtN a precompute.

---

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

**(a) THE DOMAIN SHRINKS, which is the acceleration.** See §0c: `Gamma` is at
radius 3.4 only because the conductors are meshed, and can go to about 2.1 once
they are not — 2.6× less area. On top of that, **1024 of 9361 triangles, 10.9%,
carry a coil attribute** on `examples/mastu-nke.msh`, and `CoilSize = 0.03`
against `Size = 0.30` forces a 10× graded refinement around each of them. The
two together are the case for the plan.

**(b) The conductors are resolved EXACTLY rather than to the mesh's order**, at
no mesh cost — which is where M-111's global 6.081e-03 goes. Not by adopting
anybody else's model (§0b) but by integrating MEQ's own rectangle to machine
precision instead of representing its field in a degree-k space.

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
