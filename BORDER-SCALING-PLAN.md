# The bordered merit's ruler: one scale for six kinds of constraint

Written 2026-09-15, the same day as `BORDERED-GLOBALISATION-PLAN.md` and after
it. **That file stands as the record for the globalisation route** — its §2.2
(the NPC Jacobian is solve-only, so there is no `Jᵀ`, no merit gradient, no
trust region and no JFNK), its §3.5 (`buildForms()` replaces the `DarcyForm`
outright, so nothing may switch a globalisation inside
`solveWithNormalisation()`) and its §1.3 (the axis border degenerates to an
empty row when no O-point can be located) are all load bearing here and are
cited rather than restated. Its **BG-0 gate survives into this document as
BS-0**, sharpened by three more hypotheses.

`CLAUDE.md` is the maintainer's index, `CLAUDE_FB.md` owns the borders,
`CLAUDE_HDGGS.md` owns Newton, and `MEASUREMENTS.md` M-103 is the standing
measurement. This file defers to all four.

**Line numbers are as read on 2026-09-15 at 22:30 and `GradShafranov.{hpp,cpp}`
are under active uncommitted edit** — `BorderStep` and `borderSteps()` were
added to the working tree during the afternoon this was written and moved the
`solveWithNormalisation()` body by about thirty lines. Every citation names the
function or the lambda; `grep` for the name.

**The change in one sentence.** `augmentedNorm()` weights six structurally
different constraints with **one** scale taken from the first of them, plus one
hand-derived mesh-dependent factor for the X-point; replace it with a per-border
scale derived from each border's own column and corner, and decide the X-point's
scale on a physical argument rather than on `r·h`.

---

## STATUS: BS-0 IS RUN AND ITS GATE FIRES BOTH WAYS

→ **[M-117](MEASUREMENTS.md#m-117)**. The instrument gaps §0.2 names are closed:
`reportBorderSteps()` is reachable **on failure** and not only past a converged
return, `directionFinite` is recorded and flagged per step, and the line
search's **trial ladder** prints wherever it halved more than twice — which is
the ratio that separates *a direction a reweighting would have taken* from *a
direction reweighting cannot help*.

**H2 IS EXCLUDED OUTRIGHT.** `directionFinite` is true on every step of all six
exit-2 configurations, so §0.1's row that would have made this whole route moot
does not hold and the elimination is not the defect.

**H1 AND H3 SPLIT THE SIX, WHICH IS WHY THE GATE POINTS TWO WAYS AT ONCE.**
Three cases are **X-point excursions** — TCV's null travels 1.6 m in `z` from
its seed and E's flips sign in one step — which is H3 and makes **BS-5 the
plan** for them, exactly as §7's gate says. The other three are step failures
with a finite direction, which is H1.

**AND §1.3's CENTRAL CLAIM IS MEASURED AND IS STRONGER THAN THE DOCUMENT PUT
IT.** M-117 reads the merit's composition off the same traces: on C MAST the
border dominates, `g·axis` **0.75** of a merit of 0.843, while on D TCV the
field does, `g·axis` −1.31e-04 against a `‖R‖` of **2.953**. **The border's
share of the merit is three orders of magnitude apart between two machines in
one benchmark**, so there is no single `gamma` that can be right for both — and
that unifies M-113's inversion, M-114's local optimum and `BorderMeritWeight`'s
one-sided saturation, all three of which were taken on **one** machine at one
end of that range.

**WHAT HAS SINCE MOVED THE CASES IS NOT A RULER.** `setBorderRegularisation()`
removes the singular-Jacobian throw (**[M-123](MEASUREMENTS.md#m-123)**) and an
initial guess spreading `I_p` over an **ellipse** converges MAST-U from its own
file (**[M-124](MEASUREMENTS.md#m-124)**), **necessary together and neither
sufficient**. So BS-1 to BS-4 are unstarted and are no longer the first thing to
reach for; BS-5, which this plan makes conditional, is the stage M-117 selects.

## 0. BS-0 — THE GATE, AND THE INSTRUMENT IS ALREADY BUILT

Put first because it decides whether anything below is the right repair, and
because the run that settles it now costs one printer and six invocations.

### 0.1 Three hypotheses fit every observation to date

| | hypothesis | what it predicts | what fixes it |
|---|---|---|---|
| **H1** | the **merit is badly scaled**, so Armijo rejects steps that genuinely improve the field | the accepted `damping` is small at every iteration; the field term `‖R‖` *falls* across a step while the total merit does not, or a border term dominates the total | this document |
| **H2** | the **direction is non-finite** — `y` and `z` are unchecked where the bordered solution `solved(i)` is checked (`GradShafranov.cpp:6963`) | `directionFinite` false at the failing iteration, and the twelve halvings fail for a reason unrelated to step length | the bordered elimination, not the line search |
| **H3** | `refreshXPoint()` **relocates the X-point discontinuously**, so the merit is not a continuous function of the step at all | `xR`, `xZ` jump between iterations or between trials; a line search on a discontinuous merit cannot succeed and its failure looks like H1 | a constraint on the relocation — `BORDERED-GLOBALISATION-PLAN.md` §6.2 |

They are not exclusive and H3 can *cause* H1. The point of BS-0 is that they
want three different repairs and nothing published separates them.

### 0.2 The instrument exists and NOTHING reads it

`GradShafranovSolver::BorderStep` (`GradShafranov.hpp:2609-2643`) already carries
everything all three hypotheses need, and `borderStepHistory` is populated at
three sites in `solveWithNormalisation()`:

* the per-constraint decomposition, weights applied, pushed at the top of each
  iteration (`GradShafranov.cpp:6376-6390`) — `fieldNorm`, `axis`, `boundary`,
  `current`, `exterior`, `xPoint`, `xR`, `xZ`;
* `directionNorm` and `directionFinite`, taken **before any step length is
  tried** (`GradShafranov.cpp:7047-7057`) — H2's discriminator;
* `damping`, `trials` and `armijo`, filled in when the line search finishes
  (`GradShafranov.cpp:7153-7155` on acceptance, `:7166-7168` on the fallback to
  `bestDamping`) — H1's.

**A printer landed in `apps/meq.cpp` while this was being written**, at
`:3744-3774`, gated on `--profile` and emitting exactly the §0.3 table. So BS-0
is not "build an instrument" and is now not even "print it".

**IT IS "REACH IT ON THE RUNS THAT MATTER", AND TODAY THE SIX CANNOT.** The
printer sits at `:3744`, after the solve block; the bordered failure path
returns `SolveFailed` from the catch at **`apps/meq.cpp:2762`**, roughly a
thousand lines earlier. **So on every one of M-103's six exit-2 configurations
the table is never printed** — it is reachable only on a run that succeeded,
which is the one case where nobody needs it. Moving the print above that
`return`, or printing it from the catch, is the whole of BS-0's code.

Three further gaps, all small and all worth closing in the same edit:

* **`z` is recorded only through `y`.** The doxygen says *"`‖y‖` and `‖z‖`"*
  (`GradShafranov.hpp:2634-2636`) and the code stores `y.Norml2()` alone. `z`,
  `zB`, `zL` and the exterior `z`s are separate backsolves and any one of them
  can come back non-finite while `y` does not.
* **`bestNorm` is not recorded.** On the fallback branch the step taken is the
  least-bad one and the residual is allowed to *rise*; `armijo == false` says
  that happened but not by how much.
* **The table is behind `--profile`**, which is right for a successful run and
  wrong for a failing one: a diagnosis a reader has to know to ask for is a
  diagnosis nobody gets. Print it unconditionally on failure.

### 0.3 The ONE run, and exactly what is printed

`build/meq examples/machine-b-ffprime.toml`, shipped configuration, nothing
tuned. One table per solve, one row per iteration, written to stdout by the
driver beside `reportResiduals()` (`apps/meq.cpp:555`):

```
 it        ||R||      g*Gax     g*Gbnd      g*GIp     g*Gext      g*GX    ||y||  fin  damp  tri  arm        xR        xZ
  0    1.357e-02  ...
```

with `||R||` the unweighted field residual and the five constraint columns
already carrying their weights, so that

```
sqrt( fieldNorm^2 + axis^2 + boundary^2 + current^2 + exterior^2 + xPoint^2 )
```

reproduces `newtonResiduals()[ it ]` exactly. **That identity is the check that
the decomposition is the merit and not a second calculation of it**, and it is
the BS-0 regression (§7).

**What each hypothesis looks like in that output:**

| | `fin` | `damp` | the constraint columns | `xR, xZ` |
|---|---|---|---|---|
| **H1** | true throughout | `0.125`, `0.0625`, `0.002` — small and falling | **one column dominates** the total while `‖R‖` falls | move smoothly |
| **H2** | **false** at the failing iteration | 0, `arm` false, `tri` 12 | meaningless at that row | irrelevant |
| **H3** | true | small | the `g*GX` column spikes | **jump** between rows |

If `fin` is true, `damp` is small at every row, and no column jumps, H1 is the
answer and this document is the plan. If `fin` goes false, stop and open the
elimination. If `xR, xZ` jump, the merit is discontinuous and no reweighting of
a discontinuous function helps.

### 0.4 What would falsify the whole route before any of it is built

**If `damp` is 1.0 at most iterations on machine B, evidence 1's reading is
wrong** and the crawl is something else. §1.1 says why that reading is an
inference from a linear model rather than an observation, and this is the
measurement that turns it into one.

---

## 1. What the evidence says, and one reading of it I believe is wrong

### 1.1 The damping ladder, inferred — and it is an inference

`scratchpad/lam/machine-b-ffprime-L1000.meqlog` is the shipped machine B at full
current, warm-started, and its residual history is

```
 0  7.270102e-02    1  6.394657e-02    2  5.813085e-02    3  5.088189e-02
 4  3.717791e-02    5  3.482153e-02    6  3.261605e-02    7  3.255232e-02
 8  4.417186e-02    9  4.151532e-02   10  4.023869e-02   11  3.992551e-02
12  3.976984e-02   13  3.975043e-02      then exit 2
```

Ratios: **0.8796, 0.9091, 0.8753, 0.7307, 0.9366, 0.9367, 0.9980**, then a 36%
rise, then a stall at 0.547 relative.

For Newton on an `ℓ²` merit the direction is a descent direction and
`merit(α) ≈ merit(0)(1 − α)` — `CLAUDE_HDGGS.md`, *Why it fails, measured*,
states exactly that and upstream's swept `α = 1.2e-4` is the same arithmetic. So
those ratios read as `α ≈ 0.120, 0.091, 0.125, 0.269, 0.063, 0.063, 0.002`, i.e.
three to nine halvings, and **not one full step in seven.** That is the
coordinator's reading and I agree with its conclusion.

**But three of the seven are not within 1% of a ladder value** — `0.8796`
against `1 − 1/8 = 0.875`, `0.9091` against `0.9375` or `0.875`, `0.7307`
against `0.75` — so the linear model is being asked to resolve a factor of two
in `α` at 3% accuracy. It is evidence for *"the accepted damping is small"* and
not proof of *"the accepted damping is `1/8`"*. **`damping` is recorded in
`BorderStep` and is not printed anywhere**; print it. One run replaces the
inference with the number.

Two further readings of that log, both of which matter later:

* **The rise at iteration 8 is the `bestDamping` fallback, not divergence.**
  `GradShafranov.cpp:7160-7168` takes the least-bad trial when Armijo accepts
  none, and `armijo == false` is exactly that. A 36% rise with a reported order
  of −156 is the signature and nothing else in the loop can produce it.
* **The exit-2 throw comes after a five-iteration stall at 0.547**, not out of a
  healthy descent. Whatever kills iteration 13 has been building since 8.

### 1.2 Removing the X-point rows: what it actually produced

This is the reading I think is wrong, and the correction strengthens rather than
weakens the conclusion.

The experiment is clean — `diff examples/machine-b-ffprime.toml
scratchpad/pin/machine-b-ffprime-pinned.toml` is **one section name**,
`[boundary.xpoint]` → `[boundary.limiter]` at identical `R = 1.1, Z = −0.6`,
which `Config.cpp:1462-1463` makes mutually exclusive precisely because both pin
`ψ_bnd` and only one makes it an unknown. Six pinned runs, all converging in 2
or 3 Newton iterations. **What they converge TO is the finding:**

| | pinned `ψ_ax` | reference `ψ_ax` | O-point located? | support |
|---|---|---|---|---|
| A testtokamak | 1.406280e-01 | 8.271751e-02 | **NO** | settled |
| B ffprime | 1.869864e-01 | 9.723774e-02 | **NO** | did not settle |
| D tcv | 2.743790e-01 | 2.253393e-02 | **NO** | settled |
| E diamagnetic | 1.044338e-01 | 8.376562e-02 | yes, at ( 1.5473, **1.3740** ) | did not settle |
| G mastu shaped | **−5.062356e+00** | 9.255544e-02 | **NO** | settled |
| C mast | 5.895752e-02 | 5.795820e-02 | yes, and **refused** | did not settle |

**Four of the six converge to a wall-hugging annulus with no interior extremum
of `ψ_h` anywhere on the mesh** — the driver's own warning fires on A, B, D and
G-shaped. That is not "the wrong equilibrium" in a general sense. It is
**exactly** the degenerate branch `BORDERED-GLOBALISATION-PLAN.md` §1.3
identifies: `locateAxisPoint()` fails, `peakAt()` falls through to the nodal
maximum with `constraintLocated` false (`GradShafranov.cpp:5872-5905`), the
border row is built **empty** (`:6483-6486`), and the axis constraint degenerates
to `ψ_ax ← max ψ_h` — a Picard update on a system with one fewer real
constraint. Machine B pinned reports `constraint psi_ax − max psi_h = 1.260e-11`
and converges beautifully **to a solution of that degenerate system**.

And E's "located axis" is at `z = +1.374 m` where the reference's is at
`z = +0.000000`; on a machine with `R₀ ≈ 1.26, a ≈ 0.68` that is a critical
point somewhere near the top of the domain, not the magnetic axis. So five of
six are wrong and the sixth, C, was already right **unpinned** — M-103's C
filament row reads `ψ_ax = 5.896474e-02`, `ψ_bnd = −2.356703e-02` against the
pinned run's `5.895752e-02` and `−2.356983e-02`, agreeing to five figures. C is
a null control, not a success.

**So the correction is this.** Evidence 2 does not show that the X-point rows
are where the trouble is *concentrated*. It shows that **the X-point rows are
the only thing standing between these configurations and the annulus branch**,
and that removing them lets the iteration fall into it cleanly and quickly.
`ψ_bnd` at a *prescribed* limiter is a datum; `ψ_bnd` at a *located* X-point is
an unknown that couples the topology to the solve, and it is the coupling that
is hard and also the coupling that is doing the physics.

**Two consequences, and they point the same way as the coordinator's ranking for
a different and stronger reason.**

1. **Pin-then-release is correctly withdrawn**, and the reason is sharper than
   "it converges to the wrong equilibrium": the pinned problem's solution has no
   O-point, so it is not on the same solution branch as the target at all. A
   homotopy needs its endpoints on one branch.
2. **The X-point rows' scale is the highest-value question in this document**,
   because they carry the topology and because §3.2 shows the natural scaling
   rule is *undefined* for exactly those two rows.

It also **unifies two of M-103's four outcome classes**.
`BORDERED-GLOBALISATION-PLAN.md` §1.2 sorts exit-2 and exit-0-no-O-point into
different bins, with only the first
a globalisation candidate. Evidence 2 shows a one-word change moving B, D and E
from the first bin into the second. **They are one failure seen from two sides:
the iteration is being pulled towards the annulus branch, and the X-point rows
are what makes that pull expensive rather than what causes it.** That is a real
correction to my earlier taxonomy and it is recorded here rather than by editing
that file.

### 1.3 One `gamma` for six kinds — confirmed, and it is narrower than it looks

`augmentedNorm()` (`GradShafranov.cpp:6184-6195`) is

```
sqrt( ||R||^2
      + gamma^2 ( cAx^2 + cBnd^2 + cCurrent^2 )
      + gamma^2 sum_m cModes_m^2
      + gamma^2 xScale^2 ( cXr^2 + cXz^2 ) )
```

and `gamma = initialColumn.Norml2()` (`:6160`), where `initialColumn` is
`∂R/∂ψ_ax` obtained by a **central difference** at the first iterate
(`:6151-6158`). The comment at `:6167-6183` concedes the exterior case in its
own words —
*"there is no second natural scale to give it"* — and derives `xScale = r·h`
from `∇̄ψ = r q` for the X-point.

**And `gamma` touches nothing but the merit.** `grep -n '\bgamma\b'` over
`GradShafranov.cpp`, excluding `gammaHMarker`, gives eleven hits: the two
comments, the definition, the five uses inside `augmentedNorm`, and the five
`BorderStep` records. It never enters `rowDot`, `cornerEntry`, `constraintAt`,
`columnZ` or the dense solve. The same is true of `xScale`
(`:5317, 5476-5482, 6180, 6194, 6385`).

**That is provable and it is the single most useful fact in this document.** The
bordered elimination solves `M δp = f` with `M_ij = D_ij − b_i·z_j` and
`f_i = b_i·y − G_i` (`GradShafranov.cpp:6626-6644`). Scaling constraint row `i`
consistently — `G_i → s G_i`, `b_i → s b_i`, `D_ij → s D_ij` — scales row `i` of
`[ M | f ]` by `s` and leaves `δp` **exactly unchanged**. So:

* **the Newton direction is independent of any consistent row scaling**, and in
  MEQ it is independent of `gamma` outright, because `gamma` is not applied to
  the rows at all;
* **`gamma` enters only two decisions: which trial the Armijo test accepts, and
  when `norm <= target` fires**;
* so a rescaling is a pure line-search-and-stopping change with **zero** effect
  on the direction — which bounds the risk precisely and is what makes §5's
  question answerable rather than rhetorical.

The code already knows the first half: `constraintAt()`'s comment at `:6847-6849`
says an X-point scale applied *there* *"would be a row scaling of the dense
system, which changes the step rather than only the yardstick"*. It is warning
against scaling `G_i` **alone**, which is inconsistent and does change `δp`. A
scale applied to the whole row does not, and `gamma` is applied to neither.

### 1.4 AND CONTINUATION IN `μ₀ I_p` IS DEAD — my own BG-7, falsified

`scratchpad/lam.log` is a plasma-current ramp over
`λ ∈ { 1.00, 0.90, 0.75, 0.60, 0.50, 0.40, 0.30, 0.20, 0.10, 0.05 }`:

* **machine B: all ten values fail identically**, *"no finite damped step"*, from
  `I_p = 2.5e+05 A` down to `1.25e+04 A`. There is no `λ` at which the shipped
  machine B is solvable.
* **machine E: the first six fail the same way**; `λ = 0.30` and `λ = 0.20`
  converge in 13 and 12 iterations to `ψ_ax = 7.4e-04` and `8.3e-04` with
  `ψ_bnd = −6.87e-02` and `−6.68e-02` — **refused by the topology guard**, axis
  flux ABSENT, X drift 1.093 m. So the ramp reaches a degenerate branch rather
  than the equilibrium.

`BORDERED-GLOBALISATION-PLAN.md` §5 ranks `μ₀ I_p` continuation fourth and
argues it is sanctioned by `CLAUDE_FB.md` (*"Ramping `I_p` is the exception and
is fine"*). **That argument is about legitimacy and it is still right; the route
is nonetheless measured dead on the two machines it was aimed at.** Its own
falsifier — *"if `θ` has to fall below about 1e-2 to take the first step, there
is a limit point"* — is met in the strongest possible form: at `λ = 0.05` there
is still no first step. **BG-7 is withdrawn.** Recorded rather than deleted,
because that is the standing rule and because the *shape* of the finding is
informative: the failure is not amplitude-dependent, so it is not an eigenvalue
balance and it is not stiffness.

### 1.5 Two side findings from the same sweep, neither of them scaling

Both are outside this plan's scope and both are cheap to check, so they are
recorded here rather than lost.

**`machine-c-mast` may be a false refusal.** MEQ solves it — axis at
`( 1.0156, −0.0000 )` against the reference's `( 1.006896, 0.000000 )`, `ψ_ax`
`5.896474e-02` against `5.795820e-02` (**1.7%**), `ψ_bnd` `−2.356703e-02`
against `−2.359584e-02` (**0.12%**), X-point `9.473e-03 m` from its seed — and
then refuses it, because `−ψ_bnd/span = +0.2856 > 0` puts `r = 0` on the plasma
side of the level set. **The reference's own `ψ_bndry` is negative too**, so the
guard is refusing the equilibrium `freegs4e` reports. `checkAxisSource()`
evaluates `nonlinearSource->f( r, z, 0.0 )` directly
(`GradShafranov.cpp:2777-2779`) and **does not consult `plasmaComponentMask`**,
so it asks a *pointwise* question where the solve under
`PlasmaConnectivity::Component` places current by *connectivity*. Whether that
makes it a false refusal on C is one call to `plasmaComponentHolds()` on the
element `worstR, worstZ` names. **Not established, and not a scaling question.**

**The X-point runs to a spurious null on A, and the number is specific.** A's
shipped run converges in 2 iterations with `ψ_bnd = 1.100691e-03` at a located
X-point of `( 0.605600, +0.377959 )`, **1.096 m from the seed** and on the wrong
side of the midplane, against a reference `ψ_bndry` of `3.240412e-02`. E's ramp
rows report an X drift of `1.093 m`. So the X-point does not wander — it
**relocates**, and it takes `ψ_bnd` with it, which collapses the span and
produces the annulus. That is H3, and BS-0 is what tells H3 from H1.

---

## 2. The scaling problem, stated

### 2.1 Six kinds, six units

| border | unknown `u_i` | constraint `G_i` | corner `D_ii` | column `c_i = ∂R/∂u_i` |
|---|---|---|---|---|
| `ψ_ax` | a flux | `s − ψ_h( x* )`, a flux | `1` exactly (`:6398`) | assembled, `assembleNormalisationColumn` (`:4951`) |
| `ψ_bnd` | a flux | `s_B − ψ_h( x_lim )`, a flux | `1` (`:6833`) | assembled, same routine (`:6894-6899`) |
| current | a **dimensionless** profile scale `λ` | `∫F/r − μ₀I_p`, **an ampere-metre** | `( ∫F/r )/λ` (`:6831-6832`) | assembled, `assembleCurrentColumn` (`:4599`) |
| exterior, `N` of them | a Gegenbauer coefficient | `T_m`, a transmission integral | `blockEntry( m )`, diagonal (`:6836`) | assembled exactly, once per mesh (`:6054`), and **constant** |
| X-point, 2 of them | a **length** | `q_i( x_X )`, a flux gradient over `r` | `∇q`, the only non-diagonal block (`:6771-6777`) | **exactly zero** (`xZeroColumn`, `:6838-6842`) |

`augmentedNorm` gives the first four the *same* `gamma`, which is `‖∂R/∂ψ_ax‖`,
and the fifth `gamma·r·h`.

**The current row is the clearest unit error.** `gamma` has units of
`[R]/[flux]`; `G_L` has units of `μ₀I_p`. `gamma·G_L` is therefore
`[R]·[μ₀I_p]/[flux]`, which is not `[R]`. On DIII-D `μ₀I_p` is `1.26` in SI and
`ψ` is `O(1e-1)`, so the mismatch is a factor of ten there and would be a factor
of `1e6` in a code that spoke amperes — which is one of the reasons
`setPlasmaCurrent()` takes `μ₀I_p` rather than `I_p` (`CLAUDE_FB.md`, *THE
ARGUMENT IS `μ₀ I_p`, NOT `I_p`*). The convention hides the error; it does not
remove it.

### 2.2 What the right scale is, derived

`gamma_i G_i` must have `R`'s units, and the quantity that converts them is the
**field residual's response per unit of the constraint's own residual**:

```
gamma_i = || c_i || / | D_ii |
```

Read it as: `G_i/D_ii` is the perturbation of `u_i` that would zero constraint
`i` on its own, and `‖c_i‖` converts that perturbation into `R`'s units. Check
it against the five kinds:

* **`ψ_ax`**: `D = 1`, so `gamma_ax = ‖c_ax‖` — **exactly today's `gamma`**. The
  present code is this rule, correctly derived, for one border and then copied.
* **`ψ_bnd`**: `D = 1`, so `gamma_bnd = ‖c_bnd‖`. Different from `gamma_ax`
  whenever the two columns differ, which they do — `∂F/∂ψ_ax` and `∂F/∂ψ_bnd`
  are the two expressions in `CLAUDE_HDGGS.md`, *THE COLUMNS ARE AVAILABLE IN
  CLOSED FORM*, and they differ by `(Ψ − 1)` against `Ψ`.
* **current**: `gamma_L = ‖c_L‖·λ/|∫F/r|`, which cancels the ampere-metre
  exactly. **The unit error goes away by construction.**
* **exterior**: `gamma_m = ‖c_m‖/|blockEntry(m)|`, per mode. The DtN block is
  diagonal and its entries vary strongly with mode number, so the `N` modes
  currently share a scale that is right for none of them.
* **X-point**: `‖c_X‖ = 0`, so `gamma_X = 0`. **The rule is not merely
  inconvenient there — it is degenerate.** §3.2.

### 2.3 The degeneracy is a true statement, not a gap in the rule

`c_X = ∂R/∂x_X = 0` **exactly**, and the code says why: `( r_X, z_X )` reach the
field residual *only* through `ψ_bnd`, whose column is `zB`
(`GradShafranov.cpp:6838-6842`) — so `z_X = J⁻¹0 = 0` and XP-3 costs no
backsolve at all. Moving the X-point does not move the field residual.

So **no weighting of the X-point constraint derived from the field residual can
be anything but zero**, and the scale for those two rows has to come from
outside the algebra. `xScale = r·h` is exactly such an outside choice, honestly
made. The question is not whether to replace it with a derived quantity — there
is none — but whether a *better* outside choice exists.

---

## 3. What `gamma_i` should be, and what it costs

### 3.1 Cost: four of the five kinds are free

Every column the rule needs is **already computed inside the loop**:

| | where it is already computed | extra work for `gamma_i` |
|---|---|---|
| `c_ax` | `sourceColumn()` into `column` (`:5933-5960`) | one `Norml2()` |
| `c_bnd` | `columnB` (`:6894-6905`) | one `Norml2()` |
| `c_L` | `assembleCurrentColumn( unknown, columnL )` (`:6889`) | one `Norml2()` |
| `c_m` | `exteriorColumns`, assembled **once per mesh** and constant (`:6054`, `:6080`) | `N` `Norml2()`s, once |
| `D_ii` | `cornerEntry( i, i )` (`:6746-6839`) | nothing |

So per solve the rule costs `3 + N` vector norms **captured at iteration 0**, and
**no extra residual evaluation whatever**. Against a step that already spends
`N + 4` backsolves and one `ComputeH()`, that is unmeasurable.

**And it can be made cheaper than today.** `gamma` is currently taken from
`initialColumn`, computed by a central difference before the loop — **two full
`fieldResidual()` calls** (`GradShafranov.cpp:6151-6158`) — while the step
itself uses the *assembled* column from `sourceColumn()`. Reading `gamma_ax` off
iteration 0's own column instead removes those two evaluations and makes the
yardstick consistent with the Jacobian. **Do not do it in the same stage** — see
§3.4, it moves `HighBetaConvergence` — but it is available and it is a saving,
not a cost.

### 3.2 The X-point: three candidates, and `r·h` is not obviously the worst

Today: `xScale = |r_X|·h( element )`, frozen at the first refresh
(`GradShafranov.cpp:5476-5482`). It converts `q` into the change of `ψ` across
the X-point's own element, which is defensible and is **mesh-dependent**: halve
`h` and the X-point's contribution to the merit halves, so a mesh refinement
silently reweights the line search.

| | scale | mesh-free? | degenerate when | already assembled? |
|---|---|---|---|---|
| **X1** | `gamma·r·h·‖q‖` — today | **no** | never | yes |
| **X2** | from the row, `gamma/‖b_X‖` | no — `b_X` is a vector of shape functions and scales with the element | never | `xFluxDofsR/Z`, `xFluxShape` (`:6718-6725`) |
| **X3** | **the implied displacement**: `gamma_ax·‖(∇q)⁻¹ q( x_X )‖ / L`, with `L` a length from the configuration | **yes** | `∇q` singular — a degenerate null | `xFluxJacobian`, exact (`:5458-5462`) |

**X3 is the one derived from the physics rather than from the discretisation.**
`(∇q)⁻¹ q( x_X )` is the Newton displacement of the X-point itself — *how far
this point is from a true null of `q_h`*, in metres, computed from a 2×2 matrix
MEQ already assembles exactly and which `cornerEntry()` already uses as the
corner block. Dividing by a configuration length `L` — `[boundary.exterior]
Radius`, which every machine case has, or the axis-to-limiter distance — makes
the term dimensionless, and multiplying by `gamma_ax` puts it in the merit's
units. The reading is then *"the X-point is this fraction of a machine radius
from a null"*, which is what a physicist means and is what M-103's **1.096 m**
excursion needs the merit to be able to see.

**Three things about X3 that are not settled and that the document must not
pretend are:**

1. **`(∇q)⁻¹` needs a floor.** A degenerate `∇q` is exactly the case where the
   X-point is merging with another critical point, i.e. the interesting failure,
   and a pseudo-inverse threshold is a new heuristic replacing an old one. The
   cheap variant `gamma_ax·‖q‖/(‖∇q‖ L)` has no inverse and no singularity of
   its own, at the cost of being wrong when `∇q` is strongly anisotropic — which
   near a separatrix it is.
2. **Which `L`.** `[boundary.exterior] Radius` is the obvious one and is always
   present on these cases, but it is the *domain's* size rather than the
   plasma's. The distance from the located axis to the X-point seed is more
   physical and is unavailable on an iterate with no located axis — which is
   four of the six pinned runs.
3. **X1 may be right for a reason X3 is not.** `r·h` says *"a misplacement
   smaller than one element is not a misplacement"*, which is a discretisation
   statement and is arguably what a line search on a discrete problem should
   compare against. X3 says *"a misplacement is a misplacement"*. These are
   different claims about what the merit is for, and **only the measurement
   decides**: BS-2 runs both on the same configurations.

**The frozen-at-first-refresh property must survive whichever is chosen.**
`xScaleFrozen` (`:5476-5482`) exists for exactly the reason `gamma` is frozen,
and its comment says so.

### 3.3 The merit that results

```
gamma_i = || c_i || / | D_ii |          i in { psi_ax, psi_bnd, current, modes }
gamma_X = an external choice             X1 or X3

merit = sqrt( ||R||^2 + sum_i ( gamma_i G_i )^2 + ( gamma_X . X-point term )^2 )
```

All `gamma_i` **frozen at iteration 0**, exactly as one `gamma` is today.

### 3.4 The frozen-`gamma` argument survives, and here is exactly what moves

`GradShafranov.cpp:6167-6175` freezes `gamma` so that *"a `gamma` recomputed each
step would put the Jacobian's own variation into the convergence order and
manufacture one"*, and `:4077` documents the printed residual as `‖( R, γG )‖`.
**A vector of frozen scales has that property identically** — it is still one
fixed yardstick, with more than one entry. The argument is preserved rather than
replaced, and the doxygen at `:4077` needs one sentence rather than a rewrite.

**What moves, precisely:**

| | published in | borders | moves? |
|---|---|---|---|
| **M-30** `ψ_ax` vs the dimensional estimate | `HighBetaConvergence` | `ψ_ax` **alone** | **no — bit identical**, because `gamma_ax = ‖c_ax‖/1` is today's `gamma` unchanged |
| **M-31** coupled vs decoupled | `HighBetaConvergence` | `ψ_ax` alone | **no**, and the control keeps its comparability: `:6145-6148` gives `Decoupled` the same `gamma` deliberately, and the same rule gives it the same `gamma_ax` |
| **M-24** first-step blow-up | unbordered NPC | none | **no — not this path** |
| **M-32, M-33** `ψ_bnd` tables | `FreeBoundaryCoupling` | two or more | **yes**, in iteration count and in the last digits of the answer |
| **M-43** the 2×2 factorial | `FreeBoundaryCoupling` | two or more | **yes**, same |
| **XP-3, FB-7** | `XPointBorder`, `FreeBoundaryCoupling` | four to six | **yes** |

**The rule is: one border is bit-unchanged, more than one is not.** That is a
consequence of `gamma_ax` being the definition today, and it is the strongest
available regression — `HighBetaConvergence` must not move by one digit, and if
it does, the implementation is wrong rather than the measurement being stale.

**Which is why §3.1's saving must NOT be taken in the same stage.** Reading
`gamma_ax` off the assembled column instead of the differenced `initialColumn`
changes it in the eleventh digit and breaks that regression for a reason that
has nothing to do with the scaling. Take it later, on its own, with its own
before-and-after.

**And the histories stop being comparable across configurations**, which is a
real loss and should be said in the driver's output rather than discovered:
today two runs of different machines both print `‖(R, γG)‖` with one `γ`; after
this they print norms in different composite scales. The `BorderStep`
decomposition is the compensation — it prints `‖R‖` unweighted, which **is**
comparable across configurations and is what a reader should compare.

---

## 4. The affine-covariant alternative, and why it is not the recommendation

There is a formulation in which the scaling question does not arise: measure
the step rather than the residual. Deuflhard's affine-**covariant** Newton uses
`‖δx‖` — or, for a monotonicity test, the simplified-Newton correction
`‖ J(x_k)⁻¹ R( x_k + α δx ) ‖` — which is invariant under **any** row scaling of
the equations, so `gamma_i` disappears entirely.

**It fits MEQ's structural constraints better than most alternatives.**
`BORDERED-GLOBALISATION-PLAN.md` §2.2 rules out every globalisation needing
`Jᵀ`, because `DarcyNPCOperator::Jacobian::Mult` aborts
(`darcyhybridization.hpp:3756-3761`). This one needs no `Jᵀ`: the step is
already computed, and the simplified correction is **one more backsolve against
a factorisation MEQ already holds**, which is exactly what the queue-and-flush
machinery of M-98 blocks through `ArrayMult`.

**Three reasons it is not the recommendation today:**

1. **Cost.** One extra backsolve *per trial*, and the log shows up to nine
   trials per iteration. That is up to nine extra trace backsolves per step
   against today's `N + 4`.
2. **It changes the merit on every case, including the four that converge.**
   §5's hazard applies at full strength, where §3's per-border rule leaves
   single-border cases provably untouched.
3. **It does not answer the X-point question, it dissolves it** — which is
   attractive, and which also means the X-point's misplacement would be weighted
   by whatever the linear algebra implies rather than by a stated physical
   choice. Given §1.2, a *stated* choice about the X-point is worth having.

**Recorded as the fallback if BS-2 cannot settle X1 against X3.** It is the
principled answer and it is a larger change.

---

## 5. Would rescaling change which equilibrium is reported?

**It can, and the mechanism is exact rather than speculative.** §1.3 shows
`gamma` enters only the Armijo acceptance and the stopping test. Both decide
*which* step is taken, and M-26 measures three solve routes reaching discrete
solutions **9.4% apart** on an under-resolved mesh. So a reweighting is
**not** a diagnostic change: it is capable of moving the answer on a case that
converges today.

`CLAUDE.md`'s standing rule therefore binds, and the document turns it into an
assertion rather than a hope. **The four configurations that converge and agree
in M-103 are the regression:**

| | `ψ_ax` today |
|---|---|
| F DIII-D, filament | 3.759851e-01 (`machine-f-diiid.nc`, `axis_normalised_flux = 1.00000000002119`) |
| F DIII-D, shaped | M-103's 9.102e-04 relative `L2` against the shaped reference |
| G MAST-U, filament | 8.804164e-02 (`machine-g-mastu.nc`, `axis_normalised_flux = 0.999999999279991`) |
| D TCV, shaped | 2.250051e-02 (`machine-d-tcv-shaped.nc`, `axis_normalised_flux = 1.00000000003387`) |

**The acceptance is that all four reach the same equilibrium, not the same
history.** The iteration count may fall — that is the hoped-for outcome — and
`ψ_ax`, `ψ_bnd`, the located axis position and the relative `L2` against the
reference must not move beyond the solve's own tolerance. A moved answer on any
of the four is a **correctness regression and stops the stage**, whatever it does
for the six that fail.

---

## 6. Recommendation, ranked

**1. BS-0, and it is now cheap.** The instrument is built and unread. Print it,
run machine B, and separate H1 from H2 from H3. Nothing below is worth building
until the merit is known to be the binding constraint, and §1.5's X-point
relocation is a live H3 candidate with a specific number attached.

**2. BS-1, `gamma_i = ‖c_i‖/|D_ii|` for the four kinds where it is defined.**
Free — `3 + N` vector norms at iteration 0, no extra residual evaluation — and
it is the rule the code already applies correctly to `ψ_ax` and then copies. It
fixes a genuine unit error on the current row and gives the `N` exterior modes
their own diagonal. **Single-border cases are provably bit-unchanged**, which
makes it the rare change with a free, exact regression.

**3. BS-2, the X-point's scale, measured rather than argued.** It is the highest
value item and the least determined: §2.3 shows no derived scale exists, §1.2
shows these two rows carry the topology, and §3.2 gives three candidates with
three different claims about what a merit is for. Run X1 against X3 on the same
configurations and let the measurement choose.

**4. BS-3, the no-answer-moves regression**, run *before* BS-4 rather than after
— §5 is a stop condition, not a report.

**5. The affine-covariant merit (§4) is the fallback**, not the plan.

**6. Withdrawn, with the measurement recorded**: `μ₀ I_p` continuation
(§1.4, `BORDERED-GLOBALISATION-PLAN.md` §5 / BG-7), dead on both machines at
every `λ` down to 0.05.

**7. Still standing from the other file, unchanged**: BG-1's `PicardSweeps`
control is in flight and is orthogonal to all of this — it changes where the
iteration starts, where BS-1..BS-3 change which steps it accepts. If the
in-flight sweep closes B, E or D, **BS-0 still has to run**, because a basin
change that hides a mis-scaled merit leaves it to be rediscovered on the next
machine.

---

## 7. The staged plan

Stages are BS-0 to BS-6, not colliding with BG-*, FB-*, XP-*, IN-* or PE-*. Each
ends in a **measured** acceptance and names its file.

### BS-0 — Print the border decomposition, and classify the failure — **DONE, M-117**

**Do.** Make the printer at `apps/meq.cpp:3744-3774` reachable from the
bordered failure path — it is unreachable past the `return SolveFailed` at
`:2762`, which is exactly the six runs it exists for. Close the three gaps in
§0.2: print it unconditionally on failure rather than behind `--profile`, record
`max(‖z‖, ‖zB‖, ‖zL‖, ‖z_m‖)` and its finiteness beside `directionNorm`, and
record `bestNorm` on the fallback branch. Run the six M-103 exit-2
configurations and, as controls, machine A and machine F.

**Acceptance, and it is two things.**

1. **An identity, asserted**: for every iteration of every bordered solve,
   `sqrt( fieldNorm² + axis² + boundary² + current² + exterior² + xPoint² )`
   equals `newtonResiduals()[ i ]` to **1e-14 relative**. That is what makes the
   decomposition the merit rather than a second opinion about it.
2. **A classification**, as a new `MEASUREMENTS.md` anchor (M-104 if free): each
   of the eight runs against H1/H2/H3, carrying the accepted `damping` per
   iteration, which constraint column dominates, and the X-point's position per
   iteration.

**File.** `tests/convergence/SolverContract.cpp` for the identity — it is a
contract on the solver's own reporting and needs one small bordered fixture.

**Gate.** H2 dominant → stop, open the elimination. H3 dominant → the X-point
relocation is the defect and BS-5 becomes the plan. H1 → continue.

### BS-1 — Per-border `gamma_i` for the four kinds where it is defined

**Do.** Replace the scalar `gamma` with a `std::vector<double>` sized to the
border, filled at iteration 0 from `‖c_i‖/|D_ii|` using the columns the loop
already computes, frozen thereafter. `gamma_X` keeps `xScale` untouched in this
stage. **Do not** change how `gamma_ax` is obtained (§3.4).

**Acceptance, three assertions:**

1. **`HighBetaConvergence` is bit identical** — every digit of M-30's table and
   M-31's two histories, because a single-border problem has `gamma_ax` equal to
   today's `gamma` by construction. **This is the stage's own correctness
   proof**: if it moves, the implementation is wrong.
2. **The weights actually differ.** On a configuration with all six kinds live —
   `examples/machine-f-diiid.toml`, or `tests/convergence/DivertedMachine.hpp` —
   report `gamma_i/gamma_ax` per border. **The claim being tested is that the
   spread is large**; if every ratio is within a factor of two of 1, the defect
   is smaller than evidence 3 suggests and the stage buys little.
3. **The current row's unit error is gone**, checked by a scale invariance: run
   the same configuration with the profile tables scaled by `10` and the target
   `μ₀I_p` unchanged, so `λ` moves by `1/10`. `gamma_L G_L` must be **unchanged**
   under that, where today it moves by `10`. That is an exact, dimensionless
   check needing no reference.

**File.** `tests/convergence/HighBetaConvergence.cpp` for (1);
`tests/convergence/XPointBorder.cpp` for (2) and (3), which already solves the
diverted machine five times and already has every border live.

### BS-2 — The X-point's scale, X1 against X3

**Do.** Add `gamma_X` as a choice — `XPointScale::ElementFlux` (today's `r·h`)
and `XPointScale::Displacement` (§3.2's X3, with the no-inverse variant as a
third if the 2×2 solve needs a floor). A library setter, **not a TOML key**:
it changes which equilibrium is reported and `CLAUDE.md` refuses that as a file
option. Default stays `ElementFlux` until (2) below says otherwise.

**Acceptance, two measurements:**

1. **Mesh independence, which is the thing X1 provably lacks.** Solve one
   converging configuration at `h`, `h/2` and `h/4` and report
   `gamma_X·(X-point term)` at the first iterate. **X1's must fall by about two
   per refinement and X3's must not move beyond the solve's own accuracy.** If
   X3 moves with `h` the implementation is wrong; if X1 does not, the mesh
   dependence is not where §3.2 says it is and the argument for X3 weakens.
2. **The failing six**, both scales, reporting for each: converged or not, the
   accepted damping per iteration from BS-0's printer, the located axis, and the
   X-point's distance from its seed. **The number that decides it is how many of
   the six converge with a located O-point under each scale.**

**File.** `tests/convergence/XPointBorder.cpp` for (1) — it is XP-3's
acceptance and the natural home for a claim about XP-3's rows.

### BS-3 — No answer moves. RUN BEFORE BS-4

**Do.** Nothing. Re-run.

**Acceptance.** The four converging M-103 configurations of §5 reach the same
equilibrium under the chosen scaling: `ψ_ax`, `ψ_bnd` and the located axis
position within the solve's own tolerance of today's, and the relative `L2`
against each reference unmoved in its first two digits. Separately, every
driver-against-library pin holds — `theDriverRunsTheAdaptiveLoop` at 4.4e-14,
`theDriverSolvesOnACurvedBoundary` at 1.6e-16, `DriverAcceptance`'s 1.189e-16 —
since none of those is bordered and none may move at all.

**File.** `tests/convergence/DriverAcceptance.cpp` and
`tests/convergence/FreeBoundaryCoupling.cpp`. **A moved answer stops the
stage**, whatever BS-2 did for the six.

**And it re-publishes rather than re-baselines.** M-32, M-33 and M-43 are
multi-border and their iteration counts will move (§3.4). Re-take them under the
same anchors' successors; **do not edit the existing anchors**, per the standing
rule.

### BS-4 — Re-take M-103

**Do.** All seven machines, both conductor models, shipped configurations,
`OMP = MKL = 16`, both arms — the control M-96 insists on.

**Acceptance.** A new anchor with M-103's columns plus the accepted-damping
summary from BS-0's printer, refusing `absent` as well as `bad` per M-103's own
harness note. **The number that decides the route is how many of the six exit-2
configurations close with a located O-point.**

### BS-5 — The X-point relocation, conditional on BS-0 finding H3 — **SELECTED: M-117 finds H3 on three of six**

**Do.** Only if BS-0 says the X-point jumps. `refreshXPoint()` already rejects a
trial that carries the point out of the mesh; the change is to reject one that
moves it further than a stated fraction of `L`. **It is a branch control and not
a step rule** — `BORDERED-GLOBALISATION-PLAN.md` §6.2 — so it is opt-in and
reactive and it carries BS-3's regression with it.

**Acceptance.** On machine A, the located X-point stays within a stated distance
of its seed **and** an O-point is located, against today's 1.096 m and no
O-point at all.

**File.** `tests/convergence/XPointBorder.cpp`.

### BS-6 — The two side findings, scheduled separately

**Do.** (a) Check whether `checkAxisSource()` should consult
`plasmaComponentMask` before refusing: one call to `plasmaComponentHolds()` on
the element it names, on `machine-c-mast`. (b) Nothing else.

**Acceptance.** For (a), a stated answer: either C's refusal is correct and the
guard stays, or the guard is asking a pointwise question about a solve posed by
connectivity and the refusal is qualified. **Either way MEQ's C solve agrees
with `freegs4e` to 1.7% in `ψ_ax` and 0.12% in `ψ_bnd`, and that belongs in a
measurement whichever way the guard goes.**

**File.** `tests/unit/PlasmaComponentTests.cpp` if the guard moves; otherwise a
`MEASUREMENTS.md` note and no code.

---

## 8. Falsification

1. **BS-0 finds `damp = 1.0` at most iterations on machine B.** §1.1's reading
   is an inference from a linear model and this is the measurement that can
   refute it. The merit is then not what is rejecting the steps and this
   document is void.
2. **BS-0 finds `directionFinite == false`.** H2, and the defect is in the
   bordered elimination. Void, and `BORDERED-GLOBALISATION-PLAN.md` §0.1's gate
   fires instead.
3. **BS-0 finds the X-point jumping.** H3: the merit is discontinuous, and
   reweighting a discontinuous function is not a repair. BS-5 becomes the plan
   and BS-1 becomes hygiene.
4. **BS-1's assertion 2 finds every `gamma_i/gamma_ax` within a factor of two of
   1.** Then the six kinds happen to share a scale on these problems and evidence
   3 is a correctness argument without a numerical consequence. Land BS-1 anyway
   — the unit error is real — but stop expecting it to converge anything.
5. **BS-3 finds a moved answer on any of the four.** A correctness regression
   that would trade four right answers for some number of new ones. Stop.
6. **BS-4 closes fewer than two of the six.** The merit was mis-scaled and that
   was not what was blocking these machines. Record it, keep BS-1 for the unit
   error, and go back to `BORDERED-GLOBALISATION-PLAN.md`'s ladder.
7. **A quiet one, and the most likely.** §1.2 says four of six pinned runs
   converge to a branch with no O-point. **If the rescaled merit makes the six
   converge to the same annulus**, the scaling was the barrier and the annulus
   was always the attractor — a *worse* outcome than exit 2, since M-103 already
   records that *"a run that used to fail now reports a non-equilibrium as a
   success, which is strictly the worse of the two."* BS-4's acceptance is
   therefore **converged with a located O-point**, never converged alone.

---

## 9. What this document does NOT resolve

No number below exists anywhere in this tree and none is invented here.

1. **Which of H1, H2, H3 it is.** BS-0. Everything is conditional on it, and
   §1.5's 1.096 m relocation makes H3 a live contender rather than a formality.
2. **What `L` should be in X3**, and whether `(∇q)⁻¹` needs a floor or the
   no-inverse variant is enough. §3.2's three open items.
3. **Whether X1 is actually wrong.** Its mesh dependence is real; whether a
   merit *should* be mesh dependent is a genuine question and §3.2 states both
   sides. BS-2's first measurement is the only thing that settles it.
4. **Whether the four converging machines survive.** §5 makes it an assertion,
   not a prediction, and there is no reason to expect either outcome.
5. **What happens to the comparability of printed histories.** §3.4 says the
   composite scale stops being comparable across configurations and that
   `‖R‖` unweighted is the replacement, but nothing has been decided about what
   the driver should print by default or what `newtonResiduals()` should return.
6. **Whether the affine-covariant merit is affordable.** §4 estimates up to nine
   extra trace backsolves per step from the log's trial counts; it has not been
   timed, and M-98's blocked `ArrayMult` may make it much cheaper than that.
7. **Whether the annulus branch is the attractor for all six.** §8 item 7. The
   pinned runs say it is for four; nothing says what the rescaled merit would
   find.
8. **Whether `machine-c-mast` is a false refusal.** §1.5, BS-6.
9. **Whether the `‖r₀‖` denominators are comparable at all** between the pinned
   and unpinned runs quoted in §1.2 — the pinned B starts at `1.356916e-02` and
   the unpinned at `7.270102e-02`, which is two borders' difference in the merit
   as well as a different problem, and nothing separates those two causes.
10. **What the in-flight `PicardSweeps` control says.** `scratchpad/psweep.log`
    was still running when this was written and carries only the `P0` baseline.
    It is BG-1's acceptance, not BS-*'s, and it is orthogonal — but if it closes
    B, E or D it changes the priority of everything here.
