# Cut elements under the interpolatory reaction — measured, and what it does NOT settle

From `hdgdev` (`gf-interp-hdg-dev`, commit `cf573a2522`), answering CCSZ §6.8
item 5. It bears on `INTERPOLATORY-HDG-PLAN.md` §0.2, §5.2, §5.3 and §7.2.

**Read the scope paragraph first.** What I measured is a **jump** in `F` — your
`j = 0` — and your production case is `j ≥ 1`, where `F` is continuous and only
its derivative kinks. **So none of this settles §0.2.** It bounds the worst
case, it corrects one thing we had written down about your plan, and it turns
up a hazard that is specific to the interpolatory route and that §5.3 does not
name. Take the numbers as an upper bound on what a cutting discontinuity can
cost, not as an estimate of what yours will.

---

## 0. The four findings, shortest form

1. **A discontinuity lying ON a mesh line costs exactly nothing**, at every
   degree tried. Your **per-element flood-fill mask is of this kind** and is
   therefore free.
2. **A discontinuity that CUTS costs `O(h)` in `u*`, independent of `k`.** The
   degree buys nothing back.
3. **Quadrature is not a safe harbour on a cut element.** At `j = 0` the two
   routes are a coin toss — the quadrature arm is smaller at 15 of 29 cut
   interface positions and larger at 14 — because each is accurate exactly
   where the interface misses **its own** point set, and the two point sets do
   not coincide.
4. **The interpolatory route can be BLIND to an interface**, returning
   bit-identically the answer for an interface snapped to the cell boundary,
   while reporting full `k+2` superconvergence. This is a correctness hazard
   rather than an accuracy one and it is the item we would most want you to
   read.

---

## 1. Scope, stated before the numbers

We added two gates to `convdiff`'s problem 10 (CCSZ Example 4.1,
`F(u) = u³ − u` on triangles at `τ = 1`):

| flag | gate | your language |
|---|---|---|
| `-rcm 1 -rc a` | `F = 0` where `x₀ ≥ a` — a **plane** | a jump, `j = 0`, on a fixed surface |
| `-rcm 2 -rc c` | `F = 0` where `u ≥ c` — a **level set of the solution** | a jump, `j = 0`, on a moving surface |

The manufactured source carries the same gate, so the continuous problem stays
consistent and the two arms of `-rx` solve the same thing.

**Both gates are discontinuous. Neither is your `j ≥ 1`.** Your §5.2 is about
`F = O(d^j)` with `d` the distance to the edge, which for `j ≥ 1` is continuous
across it — a kink, not a jump. An interpolant of a continuous function is not
blind to a feature it does not have a node on; an interpolant of a jump is. So
finding 4 below is a `j = 0` phenomenon and finding 2's `O(h)` is the `j = 0`
corner of your `h^{j+1/2}` estimate, not a contradiction of it.

**Why the plane gate is the useful one:** on `data/inline-tri.mesh` the plane
`a = 0.5` is a mesh line at **every** refinement, so it gives a genuine uncut
control — same discontinuity, same magnitude, resolved exactly — and any other
value cuts a column of cells. That control is what makes the rest of this mean
anything, and it is the thing your `PlasmaEdgeConvergence` fixture cannot
easily produce, the plasma edge being a curve.

---

## 2. An aligned discontinuity costs nothing — and this is about your mask

Rates of `u / q / u*` over four refinements, `-rx 1`:

| `k` | ungated `F` | `F` gated ON a mesh line |
|---|---|---|
| 1 | 2.02 2.01 **3.01** | 2.02 2.01 **3.01** |
| 2 | 3.02 3.01 **4.00** | 3.02 3.01 **4.00** |
| 3 | 4.02 4.01 **5.00** | 4.02 4.01 **5.00** |

Six digits apart at every mesh, and **bit-identical** at the finest for `k = 2`.
The jump is not the problem; cutting is.

**The consequence for you is in §7.2 rather than §5.2.** Your
`sourceValue()` consults two gates and they are different animals:

* `elementCarriesPlasma( tr.ElementNo )`, the XP-1 flood-fill mask, zeroes `F`
  on **whole elements**. Its discontinuity lies on element boundaries by
  construction — it is exactly the aligned case above, and it costs nothing.
* `insidePlasma( ψ )`, the pointwise gate, is the one that cuts.

So the gate you **cannot currently express** through
`NodalReactionFunction::Eval()` is the free one, and the gate you can express is
the expensive one. That is an argument for §7.2's ask that your plan does not
make, and we think it is the strongest one available: the element mask is not a
cost you would be paying to keep, it is a cost you would be paying to lose.

**§7.2 is noted and NOT actioned here.** `Eval( const Vector &x, const Vector
&u, Vector &F )` still receives only the physical coordinate. The integrator
does have the element in hand at both entry points, and it is a signature
change rather than a mechanism, as you say. It has not been raised as a
separate ask on our side yet and it is not in this commit.

---

## 3. A cutting discontinuity costs O(h), and the degree buys nothing

**First, a methodological warning that cost us a sweep.** A FIXED interface
samples a different phase of the cell at every refinement, so a rate read off
consecutive meshes is not a rate. At `-rc 1/3` — a mesh line at no level — the
`u*` errors run 9.7585e-04, 1.8514e-04, 2.3842e-04, 2.8307e-05, i.e. the
sequence **2.40, −0.36, 3.07**. Every number is correct and the sequence means
nothing. This applies to `PlasmaEdgeConvergence` too: as `h` halves, the
plasma edge's position **within** the cells it cuts is uncontrolled, so a
four-point rate there has an uncontrolled constant by construction. We would
not trust a single such sequence to decide §0.2.

Holding the interface at a **fixed fraction of a cell** keeps the geometry
similar as `h` halves. Then `r = 3 → 4`:

| fraction of a cell | k=1 interp | k=1 quad | k=2 interp | k=2 quad |
|---|---|---|---|---|
| 0.000 (mesh line) | **3.00** | **3.00** | **4.00** | **4.00** |
| 0.125 | 0.96 | 0.96 | 1.00 | 1.00 |
| 0.250 | 1.20 | 1.17 | 1.00 | 1.50 |
| 0.375 | 1.04 | 1.04 | 1.00 | 1.00 |
| 0.500 (cell centre) | **2.96** | **2.95** | 2.01 | 2.12 |
| 0.625 | 0.97 | 0.97 | 1.01 | 1.01 |
| 0.750 | 0.88 | 0.89 | 1.01 | 1.61 |
| 0.875 | 1.05 | 1.06 | 1.01 | 1.01 |
| 1.000 (mesh line) | **3.00** | **3.00** | **4.00** | **4.00** |

**1.00 at every cut position and at both degrees.** At the fixed mesh `h = 1/32`
that is 1.87e-07 → 1.84e-04 at `k = 2`, a factor of 990; at `k = 1`,
1.33e-05 → 1.95e-04, a factor of 15. At `k = 1` the flux and potential are
barely touched; at `k = 2` they are dragged down with `u*`, their own error
having fallen below the pollution.

**Against §5.2.** Your arithmetic gives `h^{j+1/2}`, so `h^{0.5}` at `j = 0`.
We measure `h^{1.0}`. You already note the bound is optimistic by an order on
the quadrature route where you can check it; at `j = 0` on ours it is
**pessimistic by half an order**. Both observations say the same thing about
the bound — it is an estimate whose sign of error is not fixed — and neither
replaces §0.2.

**The cell centre is a genuine outlier and we do not explain it.** At
fraction 0.5 the rate is 2.96 at `k = 1` (as good as aligned) and 2.01 at
`k = 2`. It reproduces at two meshes and in both arms. Whatever it is, it means
a coarse position sweep can mislead badly, which is why the next section is
sampled at 32 positions and not 8.

---

## 4. Quadrature is not a safe harbour — and this re-aims §0.2's fallback

**First, a correction to something we wrote about your plan.** Our CCSZ §6.8
item 5 said "meq's plan is to interpolate on uncut elements and keep quadrature
where the edge cuts, which would need a per-element opt-out". **That is not
what §4.1 says.** Your switch is `enum class SourceTerm { Quadrature,
Interpolatory }` at the configuration level — fixed boundary interpolatory,
free boundary quadrature — not a per-element choice. We withdraw the
paraphrase. The measurement below still bears on §0.2, because §0.2's fallback
is "refuse the interpolatory route on the free-boundary path and keep
quadrature", and the question is whether the quadrature route is a place worth
retreating **to** on the elements that cut.

Sweeping the interface across one cell in 32 steps, `k = 2`, `h = 1/32`, the
error is a **staircase**: consecutive samples repeat to every printed digit and
then jump. The plateau boundaries are each arm's own evaluation points,
measured on the mesh rather than inferred — as a fraction of a cell width:

| | positions |
|---|---|
| enriched nodes, `k = 2` (degree 3) | 0.065 0.130 0.309 0.333 0.374 0.627 0.667 0.691 0.870 0.935 |
| the control's rule, order 6 | 0.053 0.063 0.126 0.249 0.310 0.364 0.499 0.501 0.637 0.690 0.751 0.874 0.937 0.947 |

A representative slice, `u*` relative error:

| fraction | interpolatory | quadrature | quad/interp |
|---|---|---|---|
| 0.0000 | 1.8681e-07 | 1.8681e-07 | 1.000 |
| 0.0625 | **1.8681e-07** | 6.7961e-05 | **363.8** |
| 0.1250 | 5.4097e-05 | 7.2167e-05 | 1.334 |
| 0.2500 | 7.7245e-05 | **9.8536e-07** | **0.013** |
| 0.3750 | 1.8355e-04 | 1.3557e-04 | 0.739 |
| 0.5000 | 1.3462e-06 | 8.3158e-07 | 0.618 |
| 0.7500 | 7.8002e-05 | **1.1568e-06** | **0.015** |
| 0.9375 | **1.8681e-07** | 6.8984e-05 | **369.3** |

**Of the 29 genuinely cut positions the quadrature arm is smaller at 15 and
larger at 14**, and the ratio spans 0.013 to 369. At the worst positions in the
cell it is 26% smaller; at the next-worst, 25% larger.

So at `j = 0`, retreating to quadrature on the cut elements is a coin toss
decided by where the interface happens to fall relative to a Gauss rule's
points — not a recovery of the order. **We are not telling you §0.2 will go
this way**: you measure the quadrature route at `min(k+1, j+1.5)` for `j ≥ 1`
and it demonstrably works there. What we are saying is that "quadrature is the
safe fallback" is a property of your `j ≥ 1` measurements and not of quadrature,
and it does not extend to a genuine jump.

**What would actually fix a cut element is a rule that knows where the cut is**
— a subdivision along the interface, or a moment-fitted rule. That is an
unfitted-FEM machine, it serves the **load** as much as the reaction, and
nobody has asked for it. We are not building it.

---

## 5. The blindness hazard — the part we would most want you to read

At interface positions inside the **outermost enriched node** — 0.065 of a cell
width at `k = 2` — the interpolatory arm returns `1.8681e-07`, which is the
**uncut** value. Not close to it: we diffed two entire runs, one with the
interface on the mesh line and one 1/16 of a cell inside it, with timing lines
filtered, and **every printed digit of every norm is identical**. The load's
quadrature is blind there too, not just the reaction's.

So for a band of interface positions about 6.5% of a cell wide at each end, the
method silently solves the problem with the interface **snapped to the cell
boundary**, and reports full `k+2` superconvergence while doing it. The
quadrature arm is not blind at the same position (6.7961e-05, 364× worse),
because its rule carries points at 0.053 and 0.063.

**Why this matters more for you than for us.** On a fixed-boundary problem this
is an accuracy question and a small one. On a **free-boundary** problem the
edge position is part of the answer, and a residual that cannot see the edge
move cannot be differentiated with respect to that motion. Your §5.3 looks at
the adjacent hazard — a node landing *on* the edge, and the support count moving
by one element — and concludes correctly that for `j ≥ 1` the residual stays
continuous. The hazard here is the opposite one: not a node **on** the edge but
**no node across it**, where the residual is not merely continuous but
**constant** in the edge's position. For `j ≥ 1` the pointwise gate is
continuous and we would expect this to weaken a great deal or vanish — but
`elementCarriesPlasma` is a genuine jump, and if any future mask is evaluated
pointwise rather than per element it would reintroduce exactly this.

**A concrete check, cheap in your fixture**: hold the mesh and `k` fixed, move
`Ψ_edge` by a fraction of a cell, and see whether `ψ*` moves at all. If it is
constant over a band, that band is blind. That is one sweep and it does not
need a convergence study.

---

## 6. What this does and does not license

**Does:**

* An aligned discontinuity — your element mask — is free. Measured at `k = 1, 2, 3`.
* A cutting jump costs `O(h)` in `u*`, whatever the degree.
* A four-point rate with an uncontrolled interface phase is not a rate. Use a
  fixed relative position, or many more refinements, or report the phase.
* At a jump, neither route dominates and the winner is an accident of geometry.

**Does not:**

* It does **not** answer §0.2. Your `j ≥ 1` is a kink; ours is a jump.
* It does **not** separate the source's quadrature from the reaction's. Both
  integrate a discontinuous integrand on the same elements and our two arms
  share a bit-identical load. A quadrature-order sweep on each would separate
  them; we did not run one.
* It does **not** cover your `1/R` weight, your `extraOrder`, or Remark 2.2's
  `ℓ` (§7.4) — all of which could change the constants and none of which we varied.

## 7. How to reproduce, and the cost of the instrument

```
convdiff -p 10 -dg -hb -nls 3 -tau0 1.0 -pp -rtol 1e-11 -prec 1 \
         -m ../../data/inline-tri.mesh -r <1..4> -o <k> \
         -rx {1|2} -rcm {0|1|2} -rc <threshold>
```

`-rcm 0` is CCSZ's own `F` and is the default, so nothing moves unless asked;
the five serial and three parallel `p10_*` references pass unchanged. The
manufactured source carries the gate, and `-rcm` on a problem with no reaction
is refused rather than ignored.

One timing note, since it is the arm you would retreat to: with a **fixed**
interface the quadrature control costs about 1.2× the interpolatory route at
`k = 2`, `h = 1/32`. With the interface a **level set of the solution**, so
that it moves with the iterate, it becomes far worse — at `k = 1` on the finest
mesh the interpolatory arm takes 8.8 s and the quadrature arm had not finished
after **16 minutes**. We did not chase why, and the quadrature integrator was
never meant for production; but a free-boundary retreat to it is not obviously
cheap.

Newton survives the gate in all our arms: the gate's derivative is a delta and
is dropped, so Newton sees a one-sided Jacobian, and every gated arm converged
in 3 to 6 outer steps against 4 to 5 ungated. We note this only to say the
delta did not have to be modelled here — your `j = 0` Newton failures are about
the free-boundary system and are not reproduced or contradicted by this.
