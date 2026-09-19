# The inverse solve, and equilibrium reconstruction

Written 2026-09-16. `CLAUDE.md` is the maintainer's index, `CLAUDE_FB.md` owns
the borders and the `freegs4e` benchmark, `CLAUDE_HDGGS.md` owns Newton, and
`MEASUREMENTS.md` holds the published tables; this file defers to all four and
restates none of them. `BORDERED-GLOBALISATION-PLAN.md` and
`BORDER-SCALING-PLAN.md`, both written 2026-09-15, attack the same cold-start
failure from the solver side and §8 says how this relates to them.

**THIS FILE HAS TWO HALVES AND THEY ARE NOT THE SAME KIND OF DOCUMENT.**

| | |
|---|---|
| **Part One**, §§1–9 | a **design** for an inverse solve — given shape targets, find the coil currents. Staged, with measured acceptance criteria naming their test files, in the manner of `FREE-BOUNDARY-PLAN.md` |
| **Part Two**, §§10–15 | **exploratory research** on equilibrium reconstruction from magnetic measurements. It ends in open questions rather than stages, deliberately, and §15 says which |

**Line numbers in `src/meq/GradShafranov.{hpp,cpp}` are as read on 2026-09-16
and both files are under uncommitted edit at that moment** — the same caveat
the two globalisation plans carry, for the same reason. Every citation also
names the function or the class, which is the durable half; `grep` for the name
rather than trusting the number.

**NO NUMBER IN THIS FILE IS NEW.** Every measurement is cited to an `M-nn`
anchor, a `file:line`, or a paper with its section and equation number. Where
the design needs a number that does not exist, the experiment that would produce
it is named and the gap is left open. Nothing here has been run.

---

## 0. WHAT WOULD FALSIFY PART ONE, AND THE CHEAPEST EXPERIMENT

Put first because it decides whether the rest is worth building, and because in
this case the experiment costs no MEQ code at all.

### 0.1 The load-bearing claim

**An inverse solve converges from a cold start where the forward solve does
not, because its geometric targets pin the topology while the currents move.**
That is the claim the whole of §2.2 rests on, and the whole reason this is
worth building before it is worth having.

It is established for `freegs4e` and **not** for MEQ. `forward.py:166-188`
records it as the reason `solve_from_cold()` exists, and the evidence is that
all seven references in `tools/freegs4e-benchmark/` were *made* by the inverse
solve from cold, while the same code run forward from cold on machine A reaches
`psi_axis` 1.0441e-01 against 8.2718e-02 — **2.6e-01 relative, with three
O-points and two X-points** (`forward.py:169-176`). MEQ's own forward failures
are M-103.

**What does not follow is that MEQ's inverse solve would behave the same way.**
`freegs4e` is finite-difference Picard with the currents re-solved on every pass
(`picard.py:61`, `picard.py:136`); MEQ would be one Newton with the currents as
border unknowns. The mechanism by which the targets help — the topology cannot
run away while `ψ` at the target points is pinned — is plausibly common to both
and is not the same mechanism.

### 0.2 The experiment, and it needs no MEQ code

`solve_from_cold()` already exists and already runs both stages. What has not
been done, and what would settle the gate, is **to run it on all seven machines
and record which converge and what the released forward stage costs**.

`forward.py:178-187` asserts the inverse stage converges on all seven "because
that is how the references were made", and records for machine A that releasing
the constraint re-converges forward in **one iteration to 7.9e-08 of the
reference**. That is one machine. The experiment is `python3 forward.py` over
the seven, tabulating per machine: inverse-stage iterations, forward-stage
iterations after release, and the final agreement.

**If the released forward stage takes many iterations, or lands on a different
branch, on more than one or two of the seven, the basin-finder claim is weaker
than §2.2 states** and the motivation falls back to feature parity alone — which
is a real motivation and a much smaller one.

### 0.3 The second gate, which is MEQ's own and is cheaper still

**Is MEQ's cold-start failure actually a topology failure?** M-103's six exit-2
rows report *"no damping of the bordered Newton step gave a finite residual —
`psi_ax`"*, and `BORDERED-GLOBALISATION-PLAN.md` §0.1 establishes that **nothing
in the code distinguishes a genuine excursion from a non-finite direction**, and
that its BG-0 experiment has not been run.

**If BG-0 comes back "a non-finite direction" on most of the six, an inverse
capability does not help them either** — the defect would be in the bordered
elimination, and pinning the shape would not repair a `NaN` in `y`. So **BG-0 is
a dependency of this document's motivation**, not merely a neighbour of it. §8
says so again.

---

# PART ONE — AN INVERSE SOLVE FOR MEQ

## 1. What "inverse solve" means, in two codes that mean different things

The two reference implementations agree on the objective and disagree on almost
everything else, including on whether an X-point is a target at all. Reading
only one of them would produce a design that looks complete and is not.

### 1.1 `freegs4e`: a linear least-squares on the coil response, every Picard pass

`../freegs4e/freegs4e/control.py`, class `constrain` (`:15-123`). Read the code
rather than the docstring: the docstring says "Adjust coil currents using
constraints" and the mechanism is narrower and more interesting than that.

Each target contributes rows to one matrix `A` and one right-hand side `b`:

| target | rows | `b` | `A` |
|---|---|---|---|
| **X-point** `(R, Z)` | **two** (`:61-77`) | `−B_r(R,Z)`, `−B_z(R,Z)` at the current equilibrium | `tokamak.controlBr(R,Z)`, `tokamak.controlBz(R,Z)` — the per-coil field response |
| **isoflux** `(R₁,Z₁,R₂,Z₂)` | **one** (`:80-91`) | `ψ(R₂,Z₂) − ψ(R₁,Z₁)` | the **difference** `controlPsi(R₁,Z₁) − controlPsi(R₂,Z₂)` |
| **psival** `(R,Z,ψ*)` | **one** (`:94-98`) | `ψ* − ψ(R,Z)` | `controlPsi(R,Z)` |

and the step is solved by Tikhonov regularisation (`:109-121`), verbatim from
the source:

```
# minimise || Ax - b ||^2 + ||gamma x ||^2
#   x = (A^T A + gamma^2 I)^{-1}A^T b
current_change = dot( inv( dot( transpose(A), A ) + self.gamma**2*eye(ncontrols) ),
                      dot( transpose(A), b ) )
tokamak.controlAdjust( current_change )
```

**Four things the code says that the docstring does not.**

* **`A` is the VACUUM coil response and carries no plasma at all.**
  `controlPsi`/`controlBr`/`controlBz` (`machine.py:274-341`, `:910-976`) sum
  Green's functions over the coils of a circuit at one ampere. So the linear
  system is exact for the coil field and treats the plasma's own response to the
  current change as zero. That is what makes it cheap and it is also why it must
  be re-run every Picard pass.
* **`γ` penalises the STEP, not the current.** The unknown is
  `current_change`, applied by `controlAdjust()`. So `γ` is a damping on the
  adjustment — a trust region — and **not** the zeroth-order Tikhonov on current
  magnitude that the comment's notation suggests. The default is `1e-12`
  (`:43`), which is a minimum-norm least squares with the regularisation doing
  essentially nothing except conditioning the normal equations.
* **The X-point condition is `∇ψ = 0` and nothing else.** It says *a critical
  point is here*; it does not say *a saddle*. An O-point satisfies it equally.
* **It runs inside the Picard loop, before the first pass and after every one**
  (`picard.py:59-61`, `:135-136`), so the currents and the field converge
  together rather than in nested loops.

`tools/freegs4e-benchmark/fgsref.py:163-184` is MEQ's own `SyncConstrain`
subclass, which exists only to work around a `current_vec` staleness bug and
changes nothing about the algorithm.

### 1.2 CEDRES++: a full-space SQP whose only target is isoflux

`refs/CEDRES.pdf` — Heumann *et al.*, *J. Plasma Phys.* **81** (2015) 905810301,
open access. The paper was read in full for this document and the reading
**confirms `refs/Refs.md`'s claim with three corrections of detail**; §3.4 below
is the one that matters to MEQ.

**Table 1 (p. 5) is the 2×2 of modes** — static/evolution × direct/inverse. But
the paper numbers **eight** problems, not four, because each mode splits again on
whether `I_P` is prescribed (Problems 1–8, §§2.1–2.4), and then a second,
fully discrete set as **Inverse Problems 17–20** (§3.4). `refs/Refs.md` says
"four" and means the modes; the problems are eight.

**The cost functional is (2.12), §2.2 p. 9**, transcribed exactly:

```
K( psi ) := ( 1/2 ) sum_{i=1}^{N_desi} ( psi( r_i, z_i ) - psi( r_desi, z_desi ) )^2
```

with `(r_i, z_i)` points on a desired boundary curve `Γ_desi` and
`(r_desi, z_desi)` **one further point on the same curve whose flux is not
prescribed**. So this is a **gauge-free isoflux** — it evaluates to zero exactly
when `Γ_desi` is a `ψ`-isoline, at whatever level. MEQ has met the gauge-free
idea before, in `SurfaceFit`; `CLAUDE_INVERSION.md` owns it.

**The regulariser is (2.13)**:

```
R( I_{1,1}, ..., I_{L,N_L} ) := sum_{i=1}^{L} sum_{j=1}^{N_i} ( w_{i,j}/2 ) I_{i,j}^2
```

— a **zeroth-order Tikhonov penalty on the current magnitudes themselves**, with
a zero reference, one weight per coil, and `H` diagonal. **The word "Tikhonov"
never appears in the paper**; it calls `w_{i,j} > 0` "regularization weights",
and gives **no rule for choosing them and no example value anywhere**. Its
entire well-posedness argument is one sentence, §2.2 p. 8: *"A regularization
term ensures well posedness of the inverse problem."*

**The optimisation is FULL-SPACE and the GS equation is a genuine equality
constraint**, (3.30) §3.4 p. 21: `min_{u,y} ½yᵀKy + ½uᵀHu s.t. B(y) = F(u)`, with
`y` the finite-element coefficients **and** `λ`, and `u` the `N` individual coil
currents. The paper **explicitly rejects** the reduced-space alternative:

> "This is the main reason for which gradient type methods for a corresponding
> unconstrained optimization problem for the reduced cost function are too
> expensive: one evaluation of the gradient requires the very expensive solution
> of the nonlinear problem in the third line of (3.32)."

The SQP iteration is (3.33), a 3×3 block KKT system in `(y, u, p)` with `D_y B`
and **`D_y Bᵀ`** both appearing, solved directly; second derivatives of `B` and
`F` are dropped, so — the paper's own words, p. 22 — *"The quadratic convergence
of Newton's method deteriorates to super-linear convergence"*, on a system
*"roughly twice as large as the algebraic system of a Newton iteration of the
direct problem"*.

**Targets: isoflux ONLY.** There is no `∇ψ = 0` condition, no saddle constraint,
no strike point — the word "strike" does not occur. An X-point is obtained by
asking for a `Γ_desi` that *"has at least one corner"* (§2.2 p. 8) and letting
the isoflux fit produce it. The paper acknowledges the gap: *"Clearly, it is also
possible to define other cost functions forcing the plasma to have other
characteristics. CEDRES++ can be easily extended in this direction."*

**Inequality constraints on currents and forces are PLANNED, not implemented**
(§2.2 p. 9: *"We are planning to include for example upper and lower bounds on
the currents and the forces in the coils"*), and limits are checked **post hoc**
by inspecting the answer (§4.2.3).

**AND THE PAPER PUBLISHES NO MEASUREMENT OF ANY INVERSE SOLVE.** §4.1's Table 2
— 5 Newton iterations to 3.935226e-12 on 577,415 unknowns, 368 s — is the
**direct** problem, Problem 2. The only inverse example is §4.2.3, a WEST
CRONOS post-processing whose entire quantitative content is *"It can be seen
that the CRONOS boundary (black dots) is well matched"* against Figure 8. The
one performance claim for the SQP is qualitative and unbacked: *"the overall
computing time in practical examples has about the same magnitude as the
computing time for solving the constraint for given control parameters"* — an
inverse solve costs about one direct solve, with no number behind it. **Anyone
sizing MEQ's version against CEDRES++ has nothing to size against.**

### 1.3 Where they agree, and where they differ

| | `freegs4e` | CEDRES++ |
|---|---|---|
| **objective** | agree — hit geometric shape targets by moving coil currents | agree |
| **isoflux** | one row per pair, `ψ(P₁) = ψ(P₂)`, hard | **the only** target, soft, in a sum of squares, gauge-free |
| **X-point** | **explicit**, `B_r = B_z = 0` at a prescribed point, two rows | **absent as a condition** — implied by a cornered `Γ_desi` |
| **prescribed `ψ` at a point** | `psivals`, one row | absent — the gauge is deliberately free |
| **formulation** | a **linear** least squares in the currents, on the frozen vacuum response, re-solved every Picard pass | a **nonlinear** full-space SQP with the discrete GS equation as an equality constraint |
| **the plasma's response to a current change** | **neglected** within a pass; recovered by the outer Picard | **included**, through `D_y B` |
| **regularisation** | on the **step**, `γ` default `1e-12` | on the **current magnitude**, weights `w_{i,j}`, no value published |
| **convergence** | Picard's, i.e. linear | superlinear, second derivatives dropped |
| **coil bounds** | none | planned, not built |
| **published inverse measurements** | none in the code | **none in the paper** |

**The disagreement that matters to MEQ is the X-point row**, and MEQ is in a
position neither code is in: XP-3 already makes `(R_x, Z_x)` unknowns of the
Newton with exact rows (§3.2), so MEQ can write `freegs4e`'s condition, and can
also write the condition **neither** code has — *the X-point, wherever the
solution puts it, is to be at this target position* — as an ordinary residual on
an unknown it already carries. §5 is the table.

### 1.4 Lackner §2.3: the problem is badly posed, and the falsifier

`refs/LacknerFreeBoundary.pdf`, §2.3, as indexed in `refs/Refs.md:492`: the
inverse problem is **badly posed in the sense of Hadamard**, with Fourier
truncation and Zakharov's Tikhonov regularisation as the practical remedies,
and a cautionary result in his figures 3–4 that **too many Fourier components
leave the plasma surface essentially unchanged while making the field near the
conductors wild**.

**That last result is this design's falsifier and it is the reason §4.3 puts the
regularisation where it does.** Restated in MEQ's terms: the map from coil
currents to plasma shape has a large near-null space, so a shape that is hit to
high accuracy says nothing about whether the currents are sensible. A design
that reports only the shape residual will report success on a current set no
machine could carry.

**The concrete consequence for the acceptance criteria**: every stage in §9 that
asserts a shape residual **also** asserts something about the currents — their
norm, or their agreement with an independently-known answer. A shape-only
acceptance is exactly the test Lackner's figures say will pass on a bad answer.

**I have not read the Lackner PDF myself for this document**; §1.4 is taken from
`refs/Refs.md:492`, which is a reading of the paper made by this project, and
the PDF is on disk at `refs/LacknerFreeBoundary.pdf`. **The subagent tasked with
§2.3 reported on it; see §15 for what remains unverified.**

---

## 2. Three motivations, and the second is measured

### 2.1 Feature parity, which is the weakest of the three

Every production free-boundary code has this and MEQ does not. `[[coils]]`
`Current` is an input (`src/meq/Config.cpp:933-985`), so MEQ can answer *what
equilibrium do these currents make* and cannot answer *what currents make this
equilibrium*, which is the question a scenario designer actually asks. That is a
real gap and it is not on its own an argument for doing it now.

### 2.2 THE INVERSE SOLVE IS A BASIN FINDER, AND THIS IS THE PRIMARY DRIVER

**MEQ fails to converge from a cold start on four of seven machine cases.**
M-103: B, D and E on both conductor arms and C under shaped conductors exit 2
with *"no damping of the bordered Newton step gave a finite residual — `psi_ax`"*;
A and G-shaped exit 0 while writing a `.nc` carrying **no
`axis_normalised_flux` at all**, which M-103 is emphatic is not convergence.
Only DIII-D closes on both arms.

**AND THE REFERENCE CODE HAS THE SAME PROBLEM IN THE SAME MODE**, which is the
finding that reframes this from a MEQ weakness into a property of the problem.
`tools/freegs4e-benchmark/forward.py` was written for exactly this and its
docstring records it (`:56-70`): with the current control on, so it is not the
`Ip_logic` trap, **a cold forward `freegs4e` solve on machine A converges in 87
iterations to `psi_axis` 1.0441e-01 against 8.2718e-02 — 2.6e-01 relative — with
three O-points and two X-points**, i.e. a multi-lobed plasma, profile rescaled by
`L = 0.765`. The file's own words: *"MEQ's cold-start failures in M-103 are not
a MEQ weakness: the reference code, in the same mode, misses the same branch."*

**And the way through is the inverse solve.** `solve_from_cold()`
(`forward.py:166-235`) is a two-stage chain:

1. **inverse** — hand the coils back their `control` flag, run Picard with
   `SyncConstrain( xpoints, isoflux )`. `forward.py:178-182`: *"It converges
   from cold on all seven machines — that is how the references were made —
   because its geometric targets pin the topology while the currents move."*
2. **release** — freeze the currents where stage 1 left them, drop the
   constraint, re-solve forward.

On machine A the released stage **re-converges in one iteration to 7.9e-08 of
the reference**. So *cold → inverse → release* is a complete cold-start recipe,
and **an inverse capability in MEQ would be a globalisation for MEQ's own
forward problem as much as it is a feature.**

**The mechanism, stated so it can be argued with.** A forward solve from cold
has nothing holding the topology: the plasma may grow lobes, the X-point may
collapse onto the axis — M-103 measures exactly that, *"both no-axis rows put
the located X-point more than a metre from its seed"* — and the current
constraint only fixes the *total*, which a multi-lobed plasma satisfies as
happily as a single one. An inverse solve adds equations that a wrong topology
cannot satisfy: an isoflux pair straddling the midplane, and a saddle at a named
place. The currents absorb the mismatch, and the currents are exactly the
variables a cold start has no information about.

**`forward.py:184-187` also records a continuation knob worth carrying**:
walking `γ` **up** before releasing is continuation in how hard the coil circuit
is allowed to respond, since large `γ` damps the adjustment towards doing
nothing — which is what forward *is*. So `γ` is a homotopy parameter between the
inverse and forward problems, and that is a second use for it beyond §4.4's.

**AND THERE IS A MEASURED WARNING AGAINST THE OBVIOUS ALTERNATIVE.**
`ROADMAP.md:86-90` records, from §7.18's re-measurement, that **continuation in
the coil current "works and is not needed", reaching the failing currents on a
degenerate branch**. So a plain current-homotopy is not the same thing as an
inverse solve and has already been tried; what settled that episode was
**choosing consistent inputs** — Shafranov's `B_v` for the vertical field a given
`I_p` needs. An inverse solve is the automatic version of that reasoning, which
is a better argument for it than "continuation".

### 2.3 It closes the benchmark's standing mode mismatch from the other side

`CLAUDE_FB.md:2304`, *STANDING RULE: COMPARE IN THE SAME MODE MEQ SOLVES IN*:
every `freegs4e` reference in the tree was made by an **inverse** solve, the
currents are an output of the reference and an input to MEQ, *"so the two codes
are not being asked the same question"*. The measurement that forces it, machine
A at `k = 2`:

| | `psi_axis` |
|---|---|
| reference, filaments | 8.2717514448e-02 |
| reference, `ShapedCoil` | 8.2716849781e-02 |
| **the two references against each other** | **8e-06** |
| MEQ forward, either current set | **2.4e-02 from both** |

The two references differ in conductor model and agree to 8e-06 **because the
inverse solve re-tunes to the same targets**. `forward.py` is the fix that rule
prescribes — run `freegs4e` forward too.

**An inverse capability in MEQ is the other fix, and it is a sharper test.** An
inverse-against-inverse comparison is insensitive to the conductor model by the
same 8e-06 argument, so what it measures is MEQ's **discretisation** rather than
the forward sensitivity of the equilibrium to a current set. The existing
fourteen inverse references become directly comparable, which they are not
today. This is a benefit of the feature and **not** a reason to build it before
`forward.py`'s route is exhausted.

---

## 3. What MEQ already has, and why the answer is "most of it"

### 3.1 The border is already exactly the shape an inverse solve needs

`GradShafranovSolver::solveWithNormalisation()` (`GradShafranov.cpp:5068`)
carries five kinds of border unknown today, against **one** factorisation of the
field Jacobian per Newton step:

| border | row | column | corner |
|---|---|---|---|
| `psi_ax` | potential shape functions at the located axis, **exact** | `assembleNormalisationColumn()` (`:4969`), analytic | exactly `1` under NPC |
| `psi_bnd` | a covector on the potential at the contact | differenced or analytic | — |
| `(R_x, Z_x)` — XP-3's two | `q_r`, `q_z` shape functions on the flux dofs, **exact** | **exactly zero** (`:6886`, `xZeroColumn`) | exact 2×2 |
| `mu0 I_p` | `assembleCurrentRow()` (`:4738`), analytic | `assembleCurrentColumn()` (`:4617`) | `assembleCurrentNormalisationCorner()` (`:4679`) |
| `a_n`, the Gegenbauer coefficients | the transmission integrals | `assembleExteriorColumns()` (`:4917`) — **constant in the iterate**, assembled once per mesh | **diagonal** |

and the elimination, from the comment at `:6671-6691`:

```
[ J  C ] [ dx ]   [ -R ]
[ B  D ] [ dp ] = [ -G ],   M dp = f,  M_ij = D_ij - b_i . z_j,  f_i = b_i . y - G_i
```

with `y = J^-1 R`, `z_j = J^-1 c_j`. **One factorisation and one backsolve per
border.** The dense corner is solved at `:7004-7005` by
`mfem::DenseMatrixInverse`, on a matrix of order `nBorderTotal` — for DIII-D
today, 1 + 1 + 2 + 1 + 10 modes = **15**.

**An inverse solve adds one more kind of border unknown and one more kind of
border row. Nothing else about that structure changes.**

### 3.2 The envelope theorem, which is the same trick an isoflux row wants

`GradShafranov.hpp:487-550`, `AxisConstraint::LocatedAxis`. `psi_ax = psi_h(x*)`
with `x*` a zero of `q_h`, so `x*` moves with the solution and the chain rule
gives

```
dG/dlambda = -[ dpsi_h/dlambda |_x*  +  grad( psi_h )( x* ) . dx*/dlambda ]
```

and **`grad_bar( psi ) = r q`, so `grad( psi_h )( x* ) = 0` at a zero of `q_h`
identically**. The position term vanishes, no sensitivity of the root find is
needed, and the row is the potential shape functions of `x*`'s element evaluated
at `x*` — `(k+1)(k+2)/2` entries, exact, one element.

**Two consequences for the inverse design, and the second is not obvious.**

* **An isoflux row at a PRESCRIBED point needs no envelope argument at all**,
  because the point does not move. `ψ_h(P₁) − ψ_h(P₂) = 0` is a difference of
  two potential-shape-function covectors, exact, two elements, and it is
  structurally the **same object** as the `psi_ax` row minus the root find.
  That is the cheapest possible target row and it is the one both reference
  codes use most.
* **An isoflux row at a MOVING point — "wherever the separatrix is, put it
  through here" — gets the envelope theorem back**, and so does a target on the
  X-point's own position. §5 works this through.

`GradShafranov.hpp:733-793` records that `setLimiterSurface()` and
`setBoundaryFluxPoint()` each invoke this argument once more, so the pattern is
established rather than novel.

### 3.3 THE COIL TERM IS EXACTLY LINEAR IN THE CURRENT, AND ITS COLUMN IS CONSTANT

**This is the finding that makes the whole design cheap, and it is read out of
the source rather than assumed.**

`CoilAugmentedNormalisedSource::f()` (`Coils.cpp:1219-1222`) is

```cpp
return plasmaSource->f( r, z, psi ) + coilSet->f( r, z );
```

and `CoilSet::f()` (`Coils.cpp:886-897`) is

```cpp
double density = 0.0;
for ( Coil const &c : coilList )
    if ( c.contains( r, z ) )
        density += c.currentDensity();
return permeability*r*density;
```

with `currentDensity() = I/area` and `area = 4·halfWidth·halfHeight` fixed
geometry. Therefore, **exactly**:

```
dF/dI_c  =  mu0 * r / area_c   on coil c's rectangle,   0 elsewhere.
```

**Three properties follow and each removes a cost.**

* **It does not depend on `psi`.** `CoilAugmentedNormalisedSource::dFdPsi()`
  (`Coils.cpp:1224-1227`) forwards to the plasma source alone — the coils
  contribute nothing to it. So the coil column is **independent of the
  iterate**, exactly as the exterior Gegenbauer columns are
  (`GradShafranov.cpp:4917`, and `BORDERED-GLOBALISATION-PLAN.md` §2.3 lists it
  among the four border pieces that need no `∂F/∂ψ`). **Assembled once per mesh,
  not once per Newton step.**
* **It does not depend on the plasma support.** The coil term is added
  unconditionally, outside every `elementInPlasma()` gate; `assembleCurrentRow()`
  says so in its own comment (`:4755-4759`): *"Every quantity assembled in this
  loop is about the PLASMA term … Only `meq::SourceIntegrator` needs the
  `fOutsidePlasma()` branch, because only `f()` carries the coils."* So
  `ConfineToPlasma` and the support sweeps leave the coil column alone. §6.
* **It is a `DomainLFIntegrator` with a constant coefficient on a rectangle.**
  The weak form carries the `1/r` that turns `F` into the right-hand side, so
  the coefficient the integrator sees is `mu0/area_c`, a **constant**. This is
  the simplest column in the whole border.

**The exterior conductors are a second and different case, and it is also
linear.** A conductor outside `Γ` is an `ExteriorCoilSet` member with **no
`f()` at all** (`Coils.hpp:836-935`) — *"An exterior conductor contributes
nothing to the interior equation"* — and reaches the residual only through
`conductorNormalFlux()` (`GradShafranov.cpp:1949`), which sums the set's own
closed-form fields. That is linear in the currents too and independent of the
iterate, so a **free exterior current has a constant column through the
transmission rows** instead of through the source. Both cases are constant; they
are assembled by different code. §6.3.

### 3.4 AND MEQ HAS NO `Jᵀ`, WHICH DECIDES THE ARCHITECTURE BEFORE ANY TRADE-OFF

`BORDERED-GLOBALISATION-PLAN.md` §2.2, verified against the installed header:
`mfem::DarcyNPCOperator::Jacobian::Mult` **aborts**
(`../mfem/install/include/mfem/fem/darcy/darcyhybridization.hpp:3756-3761`) —
after `ComputeH()` the local arrays hold **factored** blocks, so the Jacobian can
be *solved with* and cannot be *applied*. **There is no `Jᵀ` either.**

**CEDRES++'s architecture is therefore not available to MEQ, and this is not a
preference.** Its SQP iteration (3.33) is a 3×3 block KKT system containing
`D_y Bᵀ(yⁱ)` explicitly in the first block row. MEQ cannot form that block, cannot
apply it, and cannot hand it to a direct solver. A full-space KKT solve is
**structurally unbuildable** on the NPC path.

**What IS available is exactly what the bordered elimination uses**: `J⁻¹`
applied to a right-hand side and to each border column. That is `NPCReduce()`,
the trace solve, and `NPCRecover()`, which `DarcyNPCSolver::ArrayMult` already
blocks across several right-hand sides in one pass — worth **1.34×** at 14
columns, → **[M-98](MEASUREMENTS.md#m-98)**.

**So the recommendation in §4 is not chosen over CEDRES++'s; it is the only one
of the two MEQ can build.** That is a stronger statement than a trade-off and it
should be made first, because a reader who knows the CEDRES++ paper will
otherwise ask why its formulation was not adopted.

---

## 4. The recommended architecture

**In one sentence: the free coil currents become additional border unknowns of
the same bordered Newton, the shape targets become additional border rows, and
the regularisation lives in the dense corner — which is a matrix of order tens.**

### 4.1 The unknowns and the equations

Write `I` for the `N_free` free coil currents (or circuit currents, §4.6), and
keep every border MEQ has. The augmented system Newton closes is

```
R( x, p, I )        = 0      the hybridized field residual, as today
G_j( x, p )         = 0      the existing borders: psi_ax, psi_bnd, (R_x,Z_x), mu0 I_p, a_n
T_m( x, p )         = 0      the NEW target rows, m = 1 .. N_targ
```

with `x` the NPC unknown `(q, psi, psihat)` and `p` the existing border scalars.
The coil currents appear in `R` and **in nothing else** — no target row contains
`I` directly, because a target is a statement about the field. So:

| | |
|---|---|
| **coil column** `∂R/∂I_c` | the constant load vector of §3.3. One backsolve per free current per step; **no re-assembly** |
| **coil corner** `∂G_j/∂I_c`, `∂T_m/∂I_c` | **exactly zero**, for every existing border and every target row |
| **target row** `∂T_m/∂x` | a covector on the potential or the flux — §5 |
| **target column** | **exactly zero**, the mirror of XP-3's: a target row introduces no unknown of its own |

**The structure is the same economy XP-3 found**, `GradShafranov.hpp:1993-1999`:
*"The bordered elimination pays one backsolve per border COLUMN … So both
columns are exactly zero … An X-point costs two rows of a 4×4 and nothing
else."* Here the targets are free in the same way and the currents are not: each
free current is a genuine column and costs a backsolve, but a **constant** one.

**The cost per Newton step is therefore `N_free` extra backsolves and nothing
else** — no extra factorisation, no extra assembly after the first, and the
backsolves block through `ArrayMult` with the ones already queued.

### 4.2 Why not an outer loop, and the number that settles it

**MEQ has already built an outer loop around this solver and measured it.**
`CLAUDE_INVERSION.md`, *Driving by `q(ψ)`*: the `q(ψ)`-driven solve is a KINSOL
outer Newton whose unknowns are the coefficients of a fitted `g²`, wrapping the
bordered Newton. `DriverAcceptance::theDriverSolvesForTheToroidalField` reads
**42.5 s, 12 outer iterations, 43 equilibria** to recover a three-coefficient
closed form.

**43 equilibria for 3 outer unknowns.** A differenced outer Jacobian costs one
equilibrium per unknown per outer step, so at DIII-D's **18 coils** — and
`examples/machine-*.toml` carries 4, 4, 11, 24, 4, 18 and 26 — the same
architecture costs of order **19 equilibrium solves per outer step**, against
`N_free` **backsolves** per Newton step for the bordered route. The two are not
close, and the gap widens with the coil count.

**And that is only the cost argument.** The house objection is stronger and is
about the answer rather than the clock: `CLAUDE.md` records that
`PicardThenNewton` and `Globalisation` are not exposed as keys because three
solve routes reach discrete solutions **9.4% apart** on an under-resolved mesh,
and `CLAUDE_INVERSION.md` records that *"a converged outer solve is not evidence
the root is the answer"* — both failing configurations of the `q(ψ)` loop
converged, to 3.5e-07, on a fixed point that was not the answer. **An outer loop
adds a second convergence criterion, and a second thing that can converge to the
wrong place quietly.**

**The `q(ψ)` loop is an outer loop for a reason that does not apply here**, and
saying so is the honest version of this argument: `g²` enters the source through
`setGGPrime()`, a whole profile object, and there is no assembled route to
`∂R/∂(spline coefficient)`. **There is an assembled route to `∂R/∂I_c` — it is a
constant vector, §3.3.** That is precisely the difference between a quantity that
belongs in a border and one that does not.

**What an outer loop would buy, and it is not nothing.** It would reuse
`freegs4e`'s algorithm unchanged, need no change to `solveWithNormalisation()`,
and let the regularisation be somebody else's well-tested least squares. §9's
**IV-0** is exactly that, as a throwaway instrument rather than as the product.

### 4.3 THE REGULARISATION GOES IN THE DENSE CORNER, AND THE CORNER IS TINY

**With `N_targ` target rows and `N_free` current columns, the dense corner stops
being square**, and that is the whole of the new numerical content.

Today the corner `M` is `nBorderTotal × nBorderTotal` and is inverted outright
(`GradShafranov.cpp:7004-7005`):

```cpp
mfem::DenseMatrixInverse inverse( dense );
inverse.Mult( right, solved );
```

With the inverse borders added, `M` has `nBorderTotal + N_targ` rows and
`nBorderTotal + N_free` columns. On DIII-D with its 18 coils, 2 X-point targets
and 2 isoflux pairs that is **21 × 33** — a matrix small enough to solve by any
means at all, inside a Newton step that costs a sparse factorisation.

**So the regularisation is a change to one dense solve of order tens, and to
nothing else.** Replace the inverse with a regularised least squares

```
minimise  || M dp - f ||^2  +  || Gamma ( dp )_I ||^2
```

where `(dp)_I` is the current block of the step and `Gamma` is diagonal with one
weight per free current. Normal equations give `(MᵀM + ΓᵀΓ) dp = Mᵀf`, which is
**algebraically `freegs4e`'s `x = (AᵀA + γ²I)⁻¹Aᵀb`** (`control.py:113-118`)
applied to the Newton step in the corner rather than to a decoupled vacuum
response. Nothing in `J`, the factorisation, the backsolves or the NPC machinery
changes.

**AND ONE IMPLEMENTATION COVERS BOTH REFERENCE FORMULATIONS, WHICH IS WORTH
SAYING PLAINLY.**

* If the targets are **attainable** — `N_targ ≤ N_free` and the shape is
  reachable — the least squares has a consistent solution, the target residuals
  go to zero, and this is a genuine Newton on an underdetermined system. **The
  quadratic convergence survives**, because nothing has been dropped; `Γ`
  selects *which* root, not *whether* there is one. This is `freegs4e`'s
  formulation.
* If the targets are **not** attainable, the least squares returns the best-fit
  step, which is the Gauss–Newton step for `min ½‖T‖² + ½‖ΓI‖²`. This is
  CEDRES++'s formulation, and it degrades to superlinear for the reason CEDRES++
  gives at p. 22 — second-order terms dropped.

**The discriminator is measurable and should be reported**: whether `‖T‖` goes
to round-off or stalls. A run that stalls has been handed infeasible targets and
must say so rather than reporting a converged Newton. §4.5.

**Take a library for the dense solve.** `CLAUDE.md`'s standing preference —
*"the chance that you found the perfect Householder implementation and are able
to maintain it indefinitely into the future is low"* — and `src/meq/SurfaceFit.cpp`
already hand-rolls a Householder QR and a one-sided Jacobi SVD that the same
rule says should become Eigen. A rank-revealing decomposition here is the right
instrument twice over: it gives the regularised solve **and** it gives the
numerical rank, which is the only honest report of how many of the currents the
targets actually determine. **`mfem::DenseMatrixSVD` exists and should be
measured against Eigen before either is chosen**; that is IV-2's experiment.

### 4.4 WHAT `γ` MEANS PHYSICALLY, AND HOW TO SET IT

**It is an inverse mutual inductance, and that is not a metaphor.**

The two terms being traded are a flux residual, in webers, and a current, in
amperes. For them to be commensurable `Γ` must carry units of Wb/A — an
inductance — and the relevant one is the **mutual inductance between a coil and
a target point**, `ψ_target = M I`. So

```
Gamma_c  ~  1 / M_c,      M_c = the flux at the target points per ampere in coil c
```

and `M_c` is **already computed**: it is the constant column of §3.3, or
equivalently `CoilSet::psiOf( c, r, z )` (`Coils.hpp:741`) evaluated at the
targets. **The natural weight is available for free and needs no tuning**, which
is a better position than either reference code is in — CEDRES++ publishes no
weight at all and `freegs4e` defaults to `1e-12`, a number with no units
attached to it that works because it is effectively zero.

**What the choice means operationally:**

| `Γ → 0` | the minimum-norm Newton step. Hits the targets; says nothing about whether the currents are sensible. **This is Lackner's failure mode**, §1.4 |
| `Γ ~ 1/M` | one ampere of coil current is worth one coil-point-flux of shape error. The scale at which the two terms are the same size |
| `Γ → ∞` | the current step is damped to zero and the solve degenerates to the **forward** problem. This is `forward.py:184-187`'s continuation knob, §2.2 |

**AND THIS IS THE SAME QUESTION `BORDER-SCALING-PLAN.md` IS ABOUT.** That
document's one-sentence statement is that `augmentedNorm()` *"weights six
structurally different constraints with one scale taken from the first of
them"*, and proposes a per-border scale derived from each border's own column
and corner. An inverse border makes that **eight** kinds of constraint and adds
one whose scale has a closed-form physical answer. So BS-1's machinery is a
dependency and the inverse border is a **new datum for it**, not a new problem:
the per-border scale BS-1 wants for a coil current is `1/M_c`, derivable rather
than tuned. §8.

### 4.5 Underdetermined, infeasible, and which is the normal case

**Underdetermined is the normal case and the numbers say so.** Counted from
`examples/machine-*.toml` and `fgsref.py`'s `CASES` (`:694-773`):

| machine | coils in the TOML | targets in `fgsref.py` | rows | slack |
|---|---|---|---|---|
| A testtokamak | 4 | 2 X-points, 1 isoflux | 5 | **−1** |
| C MAST | 11 | 2 X-points, 3 isoflux | 7 | +4 |
| D TCV | 24 | 2 X-points, 2 isoflux | 6 | **+18** |
| F DIII-D | 18 | 2 X-points, 2 isoflux | 6 | **+12** |
| G MAST-U | 26 | — | — | — |

**Machine A is OVERdetermined and that is not a mistake in the table.** With four
coils and five rows there is in general no exact solution, and `freegs4e` gets an
answer anyway because its formulation is a least squares throughout. **A design
that treats the targets as hard equations fails on machine A**, which is the
first machine in the benchmark. That alone settles the choice in §4.3: the corner
solve must be a least squares, not a square solve, and the two regimes must be
handled by one code path.

**AND THE COIL COUNTS IN THAT TABLE ARE THE FLATTENED ONES, WHICH OVERSTATES THE
FREEDOM.** §4.6.

**What to do when the targets are infeasible.** Report it, do not hide it. The
least-squares corner will return a step regardless, and the run will converge to
a best match that is not the requested shape. Per `CLAUDE.md`'s refusal
discipline — *"what the driver refuses rather than approximates"* — the honest
behaviour is:

* **converged, targets met**: `‖T‖` at round-off. Report the shape residual and
  the currents.
* **converged, targets NOT met**: `‖T‖` stalls above a threshold. This is a
  legitimate answer to *what is the best this machine can do*, and it must be
  **labelled as such in the output and in the exit status**, not reported as a
  plain success. It is the inverse-solve analogue of M-103's *"EXIT 0 IS NOT
  CONVERGED"*, and the failure mode is identical: a reader excuses a shape
  residual they would refuse if it were named.
* **rank deficient**: the decomposition's numerical rank is below `N_free`.
  Report the rank and the null-space dimension. This is Lackner's warning
  made quantitative, and it is information the user needs rather than an error.

**What is NOT resolved here**: whether inequality constraints on the currents —
coil limits, force limits — belong in this design at all. CEDRES++ plans them
and does not have them (§1.2); `freegs4e` has nothing. Adding them turns the
corner solve from a least squares into a bounded least squares, which is a
different and much larger dependency. **Recommended out of scope for IV-0..IV-6
and named in §15.**

### 4.6 CIRCUITS, AND A FLAT COIL LIST IS THE WRONG INPUT FOR THIS

**This is measured and it is a hard requirement rather than a refinement.**
`CLAUDE_FB.md:79-81`, on converting the seven machines:

> "Three things the conversion has to get right, each of which silently converges
> if got wrong. The CIRCUIT TOPOLOGY, because the references come from an
> INVERSE solve and **flattening a circuit into independent coils hands it more
> freedom than the machine has**. The `control` flag per conductor, which says
> whether that current may move at all."

`freegs4e`'s `Circuit` is a list of `(label, coil, multiplier)` and its
`controlPsi` sums `multiplier*coil.controlPsi` over the group
(`machine.py:274-295`) — **one current per circuit, with fixed per-coil
multipliers**, which is how a real machine is wired. `freegs4e`'s `TestTokamak`
(`machine.py:1965-1985`) has four coils in two symmetric pairs; MEQ's
`examples/machine-a-testtokamak.toml` has **four independent `[[coils]]`
blocks**. So MEQ's flat list would give an inverse solve four freedoms where the
machine has two, and it would use them — producing an up–down asymmetric current
set for an up–down symmetric machine, converging perfectly.

**So `[[circuits]]` is not optional for the inverse solve**, and `[[coils]]` as
it stands is a correct input for the *forward* problem and an incorrect one for
this. §7 has the schema. **It is also the single largest piece of new
configuration in this design**, and it is why IV-1 is a config stage rather than
a solver stage.

---

## 5. Target types, and what each contributes as a Newton row

Ordered by how much MEQ already has. **Every row below is a covector on `x`
with a zero column**, so each costs corner entries and no backsolve.

| target | equation | row | exact? | MEQ has |
|---|---|---|---|---|
| **isoflux pair**, prescribed points | `psi_h(P1) - psi_h(P2) = 0` | difference of two potential shape-function covectors, two elements | **exact** — polynomial evaluation | the `psi_ax` row minus the root find (`GradShafranov.hpp:487-550`) |
| **flux value at a point** | `psi_h(P) - psi* = 0` | one potential covector | **exact** | the `psi_bnd` row at a prescribed contact (`setBoundaryFluxPoint`) |
| **X-point AT a named place** | `q_r(P) = 0`, `q_z(P) = 0` | two flux shape-function covectors on one element | **exact** | **XP-3's own rows** (`GradShafranov.cpp:6745-6771`), with the point held fixed instead of unknown |
| **X-point MOVED TO a target** | `R_x - R_target = 0`, `Z_x - Z_target = 0` | `-e` on the X-point border unknowns | **exact, one entry each** | XP-3 makes `(R_x,Z_x)` unknowns; this pins them |
| **strike point** | `psi_h(S) - psi_bnd = 0`, `S` on the wall | one potential covector minus the `psi_bnd` border unknown | **exact** | nothing names a wall; `[boundary.limiter] SurfaceAttribute` is the nearest |
| **prescribed separatrix** | `psi_h(P_i) - psi_bnd = 0` for `i = 1..N` | `N` rows of the strike-point kind | **exact** | as above |

**Two observations that are not obvious and are the reason to use MEQ's border
rather than to port `freegs4e`'s.**

**FIRST: THE TWO WAYS TO ASK FOR AN X-POINT ARE DIFFERENT QUESTIONS, AND ONLY
MEQ CAN ASK THE SECOND.** Row three is `freegs4e`'s: *there shall be a critical
point here*. Row four is *the X-point this solve is already tracking shall be
here* — and it is available only because XP-3 made `(R_x, Z_x)` unknowns of the
Newton (`GradShafranov.hpp:1971-2039`). CEDRES++ has neither and gets its
X-point from a cornered isoflux curve (§1.2).

They differ in a way that matters at a cold start. Row three is satisfied by
**any** critical point, an O-point included — which is precisely the degeneracy
M-103 reports, *"the X-point constraint collapsing onto a spurious point"*. Row
four inherits XP-3's discipline instead: `setXPointBoundary()` *"follows ONE
saddle, from the value given here, exactly as the axis constraint follows one
O-point"*, and its element is **frozen within a Jacobian and re-decided at each
accepted step**. **Row four is the one to build first**, on the basis that it
carries a topology guarantee row three does not.

**SECOND: A MOVING ISOFLUX POINT GETS THE ENVELOPE THEOREM BACK.** *"Put the
separatrix through this point"* is row five, `psi_h(S) = psi_bnd`, which is
exact. But *"put the plasma edge 5 cm from this wall tile"* — a gap constraint,
which is what a real shape controller asks for — makes the point itself a
function of the solution, and then the row needs `d(point)/dx`. **It is NOT free
by the envelope theorem**, because `∇ψ` does not vanish on the separatrix the
way it does at a critical point. That is a genuine extra cost and a gap
constraint is therefore **deliberately out of scope** for IV-0..IV-6. Named in
§15.

---

## 6. Interaction with `ConfineToPlasma`, the support sweeps, and the exterior DtN

### 6.1 `ConfineToPlasma`: the coil column is untouched, the target rows are not

Established in §3.3: the coil term is added outside every `elementInPlasma()`
gate, so **the coil columns are unaffected by the support**. The target rows are
covectors on `psi_h` and `q_h`, which are affected by it exactly as every other
border row already is — no new interaction.

**The live hazard is an ordering one that has already bitten twice, and the
inverse solve makes it reachable a third time.** `CLAUDE.md` records that
`XPointOuter` was red for sessions because a fixture called `setPlasmaSupport()`
on the plasma source **before** constructing the coil wrapper, so the wrapper's
flag stayed false and **XP-1's flood fill never ran** — 333 elements of private
flux region carrying a current nobody asked for. `CoilAugmentedNormalisedSource::
setPlasmaSupport()` (`Coils.cpp:1277-1286`) sets the wrapped source *and* the
base for exactly that reason.

**An inverse solve constructs and reconfigures the coil wrapper more than a
forward one does**, since the currents move. Whatever IV-3 builds must set the
support **through the wrapper**, and the acceptance must assert
`plasmaComponentWanted()` rather than assuming it — which is the lesson M-82
records and the one a 2×2 cross was needed to see.

### 6.2 The plasma support sweeps: the inverse border sits INSIDE the loop

`[solver] PlasmaSupportSweeps` is the outer loop XP-3 deliberately left standing
(`CLAUDE_FB.md:1566-1568`): freeze the support, solve, refreeze at the answer,
solve again, stop when a sweep changes nothing. **The inverse borders belong
inside a sweep, with the currents carried across sweeps as `psiAxisGuess` and the
X-point seed already are.**

**The hazard is that the support now depends on the currents**, so a sweep that
converges the currents against a frozen support may be frozen against a support
the new currents would not produce. That is the same hazard the sweep already
manages for `ψ`, and there is a measured warning that it can fail to settle:
M-103 records machine A reporting *"the plasma support took 4 sweeps and DID NOT
SETTLE"*. **An inverse solve gives the support more ways to move**, so IV-5's
acceptance must report the sweep count and assert it settles, not merely that
the targets were met.

**And the standing trap applies unchanged**: `CLAUDE.md`'s *A RESIDUAL
EVALUATION MAY NOT SIT BETWEEN `NPCGradient()` AND THE SOLVES IT ENTITLES* —
`fieldResidual()` calls `refreshPlasmaComponent()`, so a coil column solved
across a residual evaluation is solved against a Jacobian belonging to a
different plasma support. **The coil columns are constant and can be assembled
before the Jacobian exists**, which sidesteps this entirely; that is a second
reason to assemble them once per mesh rather than per step.

### 6.3 The exterior DtN: two kinds of free current, and one of them is refused

`meq::ExteriorDtN` (`src/meq/ExteriorDtN.hpp:117`) is the semicircle's
Dirichlet-to-Neumann map in the Gegenbauer basis, and it rests on the exterior
being **`Δ*`-harmonic** — source-free — outside `Γ`.

| conductor | where its current enters | free current's column |
|---|---|---|
| **inside `Γ`** — every coil in all seven machines | `CoilSet::f()`, the source | the constant load vector of §3.3 |
| **outside `Γ`** | `conductorNormalFlux()` (`GradShafranov.cpp:1949`), the transmission rows | a constant column through the transmission rows |

Both are linear and constant, and both are buildable. **The mixed case is
already known to be unreachable on these machines**: M-95 measures that putting
`Γ` between the plasma and the coils is impossible on all seven, *"the margin is
negative at every centre on every machine, −0.26 m at best"*, and
`ROADMAP.md:172-176` records that at a usable clearance the conductors cascade
and nothing can be externalised.

**So IV-3 should build the interior case only**, and `ExteriorCoilSet::clearance()`
(`Coils.hpp:924`) is the check that refuses the other — *"A conductor inside
`Γ` is not a small error: `psi_coil` is then not `Δ*`-harmonic where the
expansion assumes it is, and the run converges to a machine nobody described."*
That refusal exists and needs no new code.

**The Gegenbauer coefficients are themselves border unknowns and stay so.** A
current change moves the plasma, which moves the transmission integrals, which
moves `a_n` — all within one Newton, against one factorisation, exactly as
`setExteriorCoupling()` already arranges (`GradShafranov.hpp:2089-2127`). **No
interaction to design; the border was built to compose.**

---

## 7. The TOML schema

House convention: **TOML keys are `UpperCamelCase`**, deliberately unlike the C++
rule (`CLAUDE.md`, *Layout*). `rejectUnknownKeys` runs first, so any key that is
to get a diagnostic of its own must be listed as accepted
(`src/meq/Config.cpp:428`).

### 7.1 A coil whose current may move

Follow `[boundary.xpoint]`'s established pattern exactly: the value given is
**an initial value for an unknown, not a prescription**
(`GradShafranov.hpp:1974-1977`).

```toml
[[coils]]
Name      = "F1A"
CentreR   = 0.860800
CentreZ   = 0.168300
HalfWidth = 0.050000
HalfHeight= 0.050000
Current   = -1.3773984082e+05    # the SEED when Free = true
Free      = true                  # NEW
```

`Free` joins the accepted list at `Config.cpp:933-934`. **`Current` or
`CurrentDensity` stays required even when `Free` is true** — a Newton unknown
needs a starting value, and the existing "exactly one of the two" refusal
(`:960-973`) is unchanged. A file with no `Free` coil and an `[inverse]` block is
refused at parse, naming the contradiction.

### 7.2 Circuits, which §4.6 makes mandatory rather than optional

```toml
[[circuits]]
Name = "P1"
Free = true
Current = 1.0e5                       # the seed for the circuit's own current
Coils = [ { Name = "P1U", Multiplier =  1.0 },
          { Name = "P1L", Multiplier =  1.0 } ]
```

A circuit's members are named `[[coils]]` blocks, which then **must not**
themselves carry `Free`; naming a coil in a circuit and marking it free is
refused rather than resolved, on the same principle as `Current` and
`CurrentDensity`. A coil in no circuit is its own one-member circuit, which keeps
every existing file valid and unchanged.

**The `Multiplier` is signed and is what makes an anti-series pair expressible**
— the thing `CLAUDE_FB.md:79-81` says the flattening loses.

### 7.3 The targets

```toml
[inverse]
Gamma          = 0.0          # 0 selects the derived 1/M weight of 4.4
MaxTargetResidual = 1.0e-08   # above this, "converged but targets not met"

[[inverse.isoflux]]
R1 = 1.2   Z1 = -1.0   R2 = 1.2   Z2 = 1.0

[[inverse.xpoint]]
R = 1.2    Z = -1.0
Track = true                   # row 4 of section 5: pin the TRACKED X-point.
                               # false gives row 3, a critical point here.

[[inverse.flux]]
R = 1.05   Z = 0.0   Psi = 0.0   # the gauge, when one is wanted
```

**`Gamma = 0` meaning "derive it" rather than "no regularisation" is a
deliberate choice and should be argued with.** The alternative — `Gamma` absent
meaning derived, `Gamma = 0` meaning literally zero — is more honest and more
typing. **Not resolved**; §15.

**`[[inverse.xpoint]] Track` and `[boundary.xpoint]` interact and the interaction
must be refused or specified.** `[boundary.xpoint]` already makes `(R_x, Z_x)`
unknowns and pins `psi_bnd` to them. `Track = true` additionally pins those
unknowns to a target, which over-determines them unless a coil current is free to
absorb it — which is the whole point. **A file with `[[inverse.xpoint]] Track =
true` and no `[boundary.xpoint]` is a parse error**, since there is no tracked
X-point to pin.

### 7.4 The reporting side

An inverse run must write out what it chose. The `.nc` gains the solved currents
per circuit, the target residual per target, and the numerical rank of the corner
(§4.5). **Per §1.4 this is not optional**: a file that reports only the shape is
the report Lackner's figures 3–4 say will look perfect on a bad answer.

---

## 8. Relationship to the two globalisation plans

Both were written 2026-09-15 and both attack M-103 from the solver side. **The
verdict is: complementary, with one of them a hard dependency, and this document
is not a superset of either.**

**THE DEPENDENCY HAS SINCE BEEN DISCHARGED, AND IT DISCHARGES IN THIS PLAN'S
FAVOUR.** BG-0 and BS-0 are the same gate and it has been run →
**[M-117](MEASUREMENTS.md#m-117)**: `directionFinite` is true on **every step of
all six** exit-2 configurations, so hypothesis H2 — the non-finite direction
that would have made both plans and this one moot — **does not hold**. The six
are one fatal step (two), chronic (one) and **X-point excursions** (three),
all of which are failures of the path or the problem rather than of the
elimination.

**And the solver-side campaign that followed came back negative**, which
sharpens the case for changing the problem rather than the step:
`Globalisation::BorderedPicardThenNewton` closes **0 of 6** and makes four of
them worse (**[M-119](MEASUREMENTS.md#m-119)**), and `PicardSweeps` fixes none
and turns three into silent wrong answers (**[M-116](MEASUREMENTS.md#m-116)**).
What *has* moved cases is not a step rule at all —
`setBorderRegularisation()` plus an `I_p` spread over an **ellipse**, M-123 and
M-124, necessary together and neither sufficient — and an inverse solve is the
same species of intervention: it changes what is being solved.

**What none of that touches is §0.1.** The load-bearing claim here is still
established for `freegs4e` and **assumed for MEQ**, and the experiment that
would settle it still costs no MEQ code and has still not been run.

| | | |
|---|---|---|
| **`BORDERED-GLOBALISATION-PLAN.md`** | **complementary, and its BG-0 gates this document's motivation** | It changes **how far the step goes**; an inverse solve changes **what problem is being solved**. Its §2.2 is load bearing here — no `Jᵀ` — and §3.4 above uses it to rule out CEDRES++'s architecture outright. And its **BG-0 experiment decides whether M-103's six exit-2 rows are path failures at all**: if they are non-finite directions rather than genuine excursions, neither a globalisation nor an inverse solve helps them, and both plans need re-aiming. §0.3 |
| **`BORDER-SCALING-PLAN.md`** | **complementary, and a dependency in the other direction** | It replaces one merit scale over six constraint kinds with a per-border scale. An inverse border makes it **eight** kinds — and, unusually, brings a scale with a closed-form physical answer, `1/M_c` (§4.4). So BS-1's machinery should land first and the inverse border is a new and easy datum for it, not a new problem. Its **BS-0 gate is BG-0 restated** and the same caveat applies |

**Is an inverse capability a superset of either? No, and the claim should not be
made.** Three reasons:

* **It cannot rescue a broken elimination.** If BG-0/BS-0's hypothesis H2 holds —
  a non-finite direction out of the bordered elimination — the inverse solve
  inherits the same elimination and the same defect, on a **larger** border.
* **It makes the merit-scaling problem worse before it makes anything better.**
  Eight constraint kinds in `augmentedNorm()` against six.
* **Its own convergence is not established.** §0.1: the basin-finder claim is
  measured for `freegs4e` and assumed for MEQ.

**What it does offer that neither plan does** is a route that changes the
*problem* rather than the *step*, which is the family
`BORDERED-GLOBALISATION-PLAN.md` §2.2 says is forced on MEQ anyway: *"Three of
the four families below therefore change the operator or the problem rather than
the step rule, and that is forced rather than chosen."* **An inverse solve is a
fifth member of that family and arguably the most physical one** — it adds
information about the machine rather than damping the arithmetic.

---

## 9. Stages

Each ends in a **measured** acceptance criterion naming its test file, per
`CLAUDE.md`'s testing stance. **A test asserts the behaviour that is wanted and
fails until it is there.**

### IV-0 — the gate, in Python, costing no MEQ code

**Run `solve_from_cold()` on all seven machines** and tabulate inverse-stage
iterations, released-forward iterations, and final agreement with the reference.
§0.1. Run `BORDERED-GLOBALISATION-PLAN.md`'s BG-0 beside it.

**Acceptance**: a table in `MEASUREMENTS.md` under a **new** anchor — the highest
today is **M-103** and anchors are never renumbered. **If fewer than five of
seven release cleanly, the basin-finder motivation is not established and IV-3
onward should not start on that argument.** Harness: `tools/freegs4e-benchmark/
forward.py`, already written.

### IV-1 — circuits, and `Free`, in the configuration layer alone

`[[circuits]]`, `Free` on coils and circuits, and the refusals of §7. **No
solver change.** `Config` is MFEM-free, so this is CI-gated
(`CLAUDE.md`, *Coverage*).

**Acceptance**: `tests/unit/ConfigTests.cpp` — every refusal in §7.1–7.3 fires
with a message naming the key; a circuit's multipliers reach `meq::CoilSet` such
that the seven machines' flattened and circuited descriptions produce **the same
`CoilSet::f()` to machine precision** when the circuit currents are set to the
flattened values. That last assertion is what stops the circuit layer changing
the forward answer, and it is the `theDriverReproducesTheLibrary` pattern.

### IV-2 — the regularised corner solve, in isolation

The rectangular least squares of §4.3, with rank reporting, **as a free function
with no solver attached**. Measure `mfem::DenseMatrixSVD` against Eigen at the
sizes §4.3 gives (21 × 33 for DIII-D) and take the library `CLAUDE.md`'s standing
preference points at.

**Acceptance**: a new `tests/unit/` case — against a closed-form rank-deficient
system, the returned step is the minimum-norm solution to `1e-14`; the reported
rank is exact; and `Γ → ∞` drives the current block of the step to zero
monotonically. **Mutation-test it**, per `CLAUDE.md`: a dropped transpose and a
`γ` on the wrong block must both be caught.

### IV-3 — coil currents as border unknowns, one machine, isoflux only

The constant column of §3.3, the isoflux row of §5, interior conductors only
(§6.3). DIII-D, which is the one machine that converges forward today (M-103), so
that a failure is attributable.

**Acceptance**: `tests/convergence/InverseSolve.cpp` —
`theInverseSolveRecoversAKnownCurrentSet`. Take `examples/machine-f-diiid.toml`,
solve it forward, read the isoflux targets **off its own converged separatrix**,
perturb the currents, and assert the inverse solve returns to the original
currents. **Two assertions, and the second is §1.4's**: the shape residual at
round-off, **and the recovered currents within a stated tolerance of the
originals**. A shape-only assertion is the test Lackner's figures say will pass
on a wrong answer.

**This is a self-consistency test and its weakness should be stated in the
file**: the targets come from the answer, so it cannot detect a target the
machine cannot reach. IV-6 is the test that can.

### IV-4 — the X-point rows, both kinds

Row three and row four of §5. Row four first, per §5's argument.

**Acceptance**: `InverseSolve.cpp` — on the same machine, moving the X-point
target by a stated distance moves the solved X-point to it within a stated
tolerance, and the currents move by an amount that is **reported rather than
asserted** (there is no independent answer to assert against). Plus a control:
the coil-free arm must fail, so that the currents are shown to be doing the work
— the `theDriverSolvesOnACurvedBoundary` pattern, which pins against the library
**and** against a do-nothing arm.

### IV-5 — through the driver, with the support sweeps

`apps/meq.cpp`, `ConfineToPlasma`, `PlasmaSupportSweeps`, the exterior coupling.
§6.

**Acceptance**: `tests/convergence/DriverAcceptance.cpp` —
`theDriverSolvesAnInverseProblem` pins the driver against the same solve driven
through the library, in the manner of `theDriverRunsTheAdaptiveLoop`'s
**1.666e-16 over 2694 dofs**. And separately: the support sweep **settles**, with
its count reported — §6.2's hazard, and M-103's *"took 4 sweeps and DID NOT
SETTLE"* is the failure it guards.

### IV-6 — the seven machines, and the claim in §2.2 tested on MEQ

Run all seven from cold in inverse mode, then release and re-solve forward with
the currents frozen, which is `solve_from_cold()`'s chain in MEQ.

**Acceptance**: a table under a new `M-nn` anchor, one row per machine per
conductor arm, reporting inverse convergence, released-forward convergence, and
agreement with the `freegs4e` reference. **The headline number is how many of
the fourteen M-103 configurations this converts from red to green**, and it is
the only acceptance in this list that tests the motivation rather than the
mechanism.

**It may legitimately come out at zero.** §0.3: if M-103's failures are not
topology failures, this changes nothing about them, and that would be a finding
worth recording rather than a stage that failed.

---

# PART TWO — EQUILIBRIUM RECONSTRUCTION

**THIS HALF IS RESEARCH AND ENDS IN OPEN QUESTIONS RATHER THAN STAGES**, which
is a statement about how well the problem is understood here and not about its
importance. §15 is the question list and §16 says what was not resolved.

**The conclusion, stated first so the rest can be read against it: a
reconstruction is data-limited, not discretisation-limited, and MEQ's two
strongest properties — high order, and a solved flux `q` — buy less than they
buy anywhere else in this project.** §14 is the argument and it is quantitative.

---

## 10. The four references, and what is on disk

**THREE OF THE FOUR ARE ON DISK AND THE FOURTH IS NOT.** Checked 2026-09-16 by
listing `refs/`, not by remembering:

| file | what it actually is | in `refs/Refs.md`? |
|---|---|---|
| `refs/MHD Equilibrium Reconstruction in the DIII-D Tokamak.pdf` | **Lao, St. John, Peng, Ferron, Strait, Taylor, Meyer, Zhang & You**, *Fusion Science and Technology* **48**:2 (2005) 968–977, `10.13182/FST48-968`. 11 pp. **The primary reference and the one that answers the structural questions** | **NO** |
| `refs/EFIT2006.pdf` | **Appel, Huysmans, Lao, McCarthy, Muir, Solano, Storrs, Taylor, Zwingmann** + EFDA ITM Task Force + JET-EFDA Contributors, *A Unified Approach to Equilibrium Reconstruction*, EFDA–JET–CP(06)03-07, 33rd EPS, Rome 2006. **4 pages of content** | **NO** |
| `refs/EFIT-revisited-for-MAST.pdf` | **Appel & Lupelli**, *Equilibrium reconstruction in an iron core tokamak using a deterministic magnetisation model*, *Comput. Phys. Commun.* **223** (2018) 1–17, `10.1016/j.cpc.2017.09.016`, JET Contributors. **THE FILENAME IS WRONG ABOUT THE MACHINE** | **NO** |
| `refs/Fitzgerald-AnisotropicEFIT.pdf` | **ABSENT.** `refs/Refs.md:297` indexes it as an *"unread seed"* — Fitzgerald, Appel & Hole, *EFIT tokamak equilibria with toroidal flow and anisotropic pressure using the two-temperature guiding-centre plasma*, *Nucl. Fusion* **53** (2013) 113040, **`10.1088/0029-5515/53/11/113040`** | indexed, PDF absent |

**`EFIT2006.pdf` IS NOT THE PAPER ITS TITLE SUGGESTS AND THE DOCUMENT SHOULD NOT
BE PLANNED AROUND IT.** It is an EPS conference preprint about **software
architecture and data plumbing**: no numbered equations, no statement of the
Grad–Shafranov equation, no measurement list, no `p′`/`ff′` basis, no benchmark
numbers, nothing on ill-posedness. Its "unified" means *one executable for many
machines* — EFIT-f95 and EFIT2006 sharing an ITM hierarchical data structure
reached through IDAM, MDS+ and XML — **not one physics model reconciled across
codes.** Its §4 is a six-item list of what users suffer, and every item is an
interface complaint: *"There is no standard equilibrium reconstruction code …
There is no consensus on methods for accessing experimental databases … There is
no standard definition of units."* **It does not compare reconstruction codes at
all.**

That is worth recording rather than glossing, because the natural expectation —
a paper by EFIT's own author with "unified" in the title is the best statement of
what reconstruction requires — is wrong, and acting on it would put a
load-bearing citation on a paper that cannot carry it. **The structural answers
come from the DIII-D paper and from the iron-core paper's §2**, and §§11–12 say
which at each point.

**`refs/Refs.md` HAS NO ENTRY FOR ANY OF THE THREE PDFs.** §17 lists what to add.
Per `CLAUDE.md`, this file does not edit `Refs.md`.

**And `EFIT-revisited-for-MAST.pdf` is misnamed.** The machine is **JET**
throughout — discharge 82 937 for validation, 91 576 for the application, JET's
iron core, JET's `Nomatil RD` material, compared against COMSOL and CREATE. A
search of the paper finds **three** occurrences of "MAST", none of them a
comparison: one sentence noting that Fig. 1b and 1c plot MAST and JET flux
surfaces, the Fig. 1 caption, and the word "Master" in a thesis reference. There
is no MAST-U mention and no "revisiting" anywhere. **The nearest defence for the
filename is in the OTHER paper** — `EFIT2006.pdf` §2 says EFIT-f95 is *"based on
the EFIT source code originally developed for MAST"* — so whoever named the file
may have attached that sentence to the wrong PDF. **Cite it by its real title and
machine.**

---

## 11. The measurements, and what each contributes mathematically

**Every measurement contributes one row of one least-squares system, and the
distinction that matters is what it constrains — the boundary, or the interior.**

### 11.1 DIII-D's set, and the two counts differ

The **installed** set (§II.C, p. 971): *"67 magnetic probes located in two
toroidal planes, 41 flux loops, 6 saddle loops, and 19 Rogowski loops to measure
plasma and external coil currents"*, around *"18 external shaping coils"*
(§II.A).

The set **actually used** in the paper's own reconstructions (§IV.B, pp. 973–974)
is smaller: *"39 flux loops, 55 magnetic probes, 1 full Rogowski loop, 34
channels of MSE measurements, and kinetic profile data from Thomson scattering,
ECE, and CER"*. **The paper does not explain the reduction**; the mechanism is
the uncertainty vector, §I: *"The use of the uncertainty vector allows all
available data to be flexibly included or entirely left out in the
reconstruction"* — setting `σ → ∞` drops a channel.

| measurement | what it reads | what it constrains |
|---|---|---|
| **flux loop** | total toroidal flux at a poloidal position | `ψ` at a point outside the plasma — a **linear** functional of all currents |
| **magnetic probe** | local `B` component | `∇ψ` at a point, same linearity |
| **Rogowski** | enclosed current | `I_p` and coil currents — an integral constraint |
| **saddle loop** | *"the difference in magnetic flux between two toroidal loops"* (Appel & Lupelli §6) | a flux **difference**, which is an isoflux-shaped row |
| **MSE** | *"the local pitch angles of the magnetic field lines inside the plasma using a spectroscopic technique"* (§I) | `B_z/B_φ` **inside** the plasma — the `q` profile |
| **Li beam** (Zeeman) | edge pitch angle, *"not sensitive to the radial electric field"* | the **edge** current density |
| **Thomson / ECE / CER** | `n_e`, `T_e`, `T_i`, rotation | `P(ψ)` directly, plus a **topological** constraint from `T_e` being a flux function |

### 11.2 JET's set is different, and the difference is worth keeping

Appel & Lupelli §6: *"The magnetic data set in the standard JET intershot has **35
magnetic probes, 4 flux loops, 27 saddle loops, and 10 measured power supply
currents**."*

**Against DIII-D's 55 probes and 39 flux loops, JET uses 4 flux loops and 27
saddle loops** — the opposite balance. **THIS IS THE SPREAD, AND IT SHOULD NOT BE
SMOOTHED OVER**: there is no canonical magnetic measurement set, the numbers
differ by an order of magnitude between machines in the same quantity, and a
synthetic-diagnostic layer built against one machine's balance is not a general
tool. It is also a hint about the iron: a machine whose flux is routed through
ferromagnetic material measures flux *differences* rather than absolute flux.

### 11.3 The coil currents are FIT PARAMETERS, not inputs — and that is a design lesson

**This is the single most transferable sentence in the DIII-D paper** for anyone
coming from a forward solver, §II.B:

> "Note that in EFIT, rather than treating the external coil currents `I_ej` as
> given, `α⃗` consists of all the unknown plasma current profile parameters `α_n`
> and `γ_n` as well as the external coil currents `I_ej`. This approach allows
> all the magnetic and external coil current data to be consistently treated with
> their respective measurement uncertainties."

The currents **are** measured — that is what the 19 Rogowskis are for — and they
enter `χ²` as measurements with their own `σ`, rather than as constraints.
**So Part One's inverse solve and a reconstruction are the same solve with
different weights**: a free current with `σ → ∞` is Part One's `Free = true`, and
a current with a finite `σ` is a soft version of the same thing. That is an
architectural connection worth noticing and §13.1 returns to it.

And the vacuum-shot precondition, §II.C, which MEQ has no equivalent of:

> "A necessary condition for equilibrium reconstruction is that in the absence of
> the plasma, the predicted magnetic responses from the external coils must agree
> with the experimental measurements within the measurement uncertainties. This is
> routinely done in DIII-D using vacuum shots."

---

## 12. The structure of the problem: exactly what is linear

**The claim in the brief is confirmed explicitly and in the papers' own words,
and the mechanism is more specific than "the coefficients are linear".**

### 12.1 The equation and the basis

DIII-D Eq. (1):

```
Delta* psi = -mu0 R J_phi ,   J_phi = R P'( psi ) + mu0 F F'( psi ) / ( 4 pi^2 R )
```

**Transcribe that second term carefully**: this paper's `F` is `2 pi R B_phi/mu0`,
so `mu0` is in the **numerator**. It is *not* the `ff'/(mu0 R)` form, and copying
that form while attributing it to this paper is a unit error waiting to happen.
Appel & Lupelli's (2) uses the other convention, `(1/(mu0 R)) f f'`, with
`f = R B_phi`. **Two of the four references on disk write this term differently
and both are right in their own convention.**

The basis, DIII-D Eq. (4):

```
P'( psi ) = sum_n alpha_n y_n( x ) ,    F F'( psi ) = sum_n gamma_n y_n( x )
x = ( psi - psi_0 )/( psi_1 - psi_0 )
```

with, §II.B, *"A polynomial, a spline representation, and a local representation
are available in EFIT"* — `y_n(x) = x^n` (Eq. 15) and a **variable-tension
spline** (Eq. 16, after Cline, *Commun. ACM* **17** 218). Appel & Lupelli (3)–(4)
agree: *"a typical choice is to use `c_i = d_i = psibar_p^{i-1}`"*, with tension
splines available and declared out of scope.

### 12.2 THE LINEARISATION IS A FROZEN FLUX MAP, AND THAT IS THE WHOLE TRICK

DIII-D §II.B, first sentence:

> "Since the unknown parameters `α_n` and `γ_n` appear linearly, Eqs. (4) and (5)
> can be combined to explicitly relate the unknown parameter vector `α⃗` directly
> to the measurement vector `M⃗` through the response matrix `R̄`."

and the mechanism is Eq. (5), whose **superscripts are the entire content**:

```
C_i^{m+1}( r_i ) = sum_j G_Ci( r_i, r_ej ) I_ej
                 + int_{V^m} dR' dZ' G_Ci( r_i, r' ) J_phi[ R', psi^m( r' ), alpha^{m+1} ]
```

**The plasma volume `V^m` and the flux map `psi^m` are at cycle `m`; the
parameters `alpha^{m+1}` are the unknowns.** So `x = (psi^m - psi_0^m)/(psi_1^m -
psi_0^m)` is a **fixed function of position** within a cycle, the basis functions
`y_n(x)` are fixed fields, and the measurement response is a fixed linear
functional of the coefficients. **Fixed `psi` ⇒ fixed normalised flux ⇒ fixed
boundary ⇒ a linear least squares.**

The objective, Eq. (6): `chi^2 = sum_i ( ( M_i - C_i )/sigma_i )^2`, giving the
linear system Eq. (7): `R̄ alpha = M`, with `R̄` *"an appropriate rearrangement of
the Green induction functions and their integration over the plasma volume
weighted by the uncertainties `sigma_i`"*. Appel & Lupelli (7)–(8) are the same
object written out per current source, ending at `chi^2 = || A I - b ||`.

### 12.3 It is INTERLEAVED, not nested, and that matters

**EFIT does not converge a GS solve and then fit.** The abstract:

> "A response function formalism and a Picard linearization scheme are used to
> **efficiently combine the equilibrium and the fitting iterations** and search for
> the optimum solution vector."

and §I: *"EFIT retains the computational efficiency of the filament code approach
by **interleaving** the equilibrium and the fitting iterations with a Picard
linearization scheme."* The equilibrium Picard Eq. (3) and the fitting
linearisation Eq. (5) **share the index `m`** — one linear fit per equilibrium
cycle. Appel & Lupelli §2 describe the same loop with a third stage for the iron:
*"a loop which iterates towards a converged equilibrium force balance solution by
successively invoking the magnetisation model, the linearised Grad–Shafranov
solver and a least-squares algorithm."*

**Neither paper gives an iteration count.** DIII-D gives only a criterion —
`eps = Max|( psi^{m+1} - psi^m )/( psi_0^{m+1} - psi_1^{m+1} )| ~ 1e-8` at
129×129, `~1e-4` at 65×65. **That is a gap in both papers, not in the reading.**

### 12.4 What is NOT linear, and the split is architectural

DIII-D §IV, on a localised ECCD current channel `J_local = gamma_local
cos^2( kappa x~ )/R`:

> "`γ_local` appears linearly and can be determined similarly as other current
> profile parameters … `x_local` and `Δ_local` appear nonlinearly and are
> determined in a **separate external optimization loop**."

and the same for spline knots: *"An option is available in EFIT to optimize the
knot locations to best fit the pressure and MSE data using an **external
nonlinear optimization driver**."*

**Amplitudes go in the inner linear solve; LOCATIONS go in an outer nonlinear
driver.** That is EFIT's own architecture, stated twice, and it is the shape any
MEQ reconstruction would inherit.

### 12.5 The solve itself, and the absence worth naming

DIII-D §IV.A: *"At each iteration, the linear minimization problem is solved
using the **singular value decomposition** method by decomposing `R̄` into a
diagonal matrix."* Appel & Lupelli §2: *"In general `A` will not be square and
the **minimum norm** of Eq. (8) is computed using singular value decomposition."*

**Neither paper states a truncation threshold, a condition number, or which
singular directions are dropped.** That is the operative regularisation decision
in the whole method and it is in neither document. Whoever builds this must go to
Lao *et al.*, *Nucl. Fusion* **25** (1985) 1611 — **not on disk** — or to EFIT
itself.

**What the DIII-D paper does give is a weak-constraint mechanism with a measured
insensitivity**, §IV.A, the JT representation:

```
P'( psi )  = alpha_0 + alpha_1 x ,              0.1 alpha_0 + 0.1 alpha_1 = 0
FF'( psi ) = gamma_0 + gamma_1 x + gamma_2 x^2 , 0.1 gamma_0 + 0.1 gamma_1 + 0.1 gamma_2 = 0
```

> "The coefficients 0.1 appearing in these constraint equations describe the weak
> weighting of these constraints relative to the other magnetic measurements when
> `R̄` is decomposed. … **Increasing the coefficients from 0.1 to 0.2 only changes
> the results slightly.**"

**A soft edge condition with a weight, whose value the answer is insensitive to.**
That is the same device as Part One's `Γ` and the same reassurance CEDRES++ never
offers about its own `w_{i,j}`.

---

## 13. What MEQ has, and what MEQ lacks entirely

### 13.1 What helps, and honestly how much

| | |
|---|---|
| **`q` is a SOLVED variable** | critical points are roots of `q_h` rather than of a differentiated `ψ_h` (`CLAUDE_INVERSION.md`), so `psi_0` and `psi_1` — which define `x` and hence the whole frozen map of §12.2 — come from a field carrying the potential's own order. EFIT re-finds them by searching a 65×65 or 129×129 grid every cycle. **Genuinely better, and §14 argues it is not the binding constraint** |
| **the flux-surface machinery** | `FluxSurfaceFamily`, the tracer, the averages, the Zernike disc basis and the gauge-free fit are all built and all are things a reconstruction needs downstream — `beta_p`, `l_i`, `q(psi)`, the separatrix. `CLAUDE_INVERSION.md` |
| **a bordered Newton with exact derivatives** | `psi_ax`, `psi_bnd`, `(R_x, Z_x)`, `mu0 I_p` and the exterior coefficients are already unknowns with analytic rows. §3.1. **A profile coefficient is the same kind of object as `mu0 I_p`** — a global scalar entering through the source — so the border generalises to it |
| **coil currents as unknowns** | Part One. §11.3 shows EFIT treats them the same way, so Part One is a **prerequisite** of a reconstruction rather than a neighbour of it |
| **toroidal flow** | `CLAUDE_FLOW.md`. DIII-D §V reconstructs with rotation, Eqs. (24)–(33), adding a **third** stream function `P_omega'(psi) = sum_n omega_n y_n(x)` and keeping the same linear structure. MEQ's `meq::RotatingSource` is the forward half of exactly this. §13.3 |
| **high order** | and **§14 concludes it buys little here**, which is the one place this document contradicts the usual argument for MEQ |

### 13.2 What MEQ lacks entirely — three layers, none of them started

**(a) A SYNTHETIC DIAGNOSTIC FORWARD MODEL.** Given an equilibrium, compute what
each probe and loop *would* read. **MEQ has the pieces for exactly one kind of
sensor and nothing for the rest.** `meq::GridSampler` samples `psi` and `B` on a
grid (`src/meq/Sampler.{hpp,cpp}`), so a flux loop or a probe **inside the mesh**
is a sample. What is missing:

* sensors **outside `Gamma`**, which is where most of them are. `ExteriorDtN`
  can evaluate the exterior field from the converged coefficients
  (`ExteriorDtN.hpp`, *"the field at an exterior point"*), so this is reachable —
  but nothing calls it for this purpose and nothing has measured its accuracy at
  a sensor position.
* a **sensor description** — position, orientation, type. DIII-D's probes have
  an orientation angle; Appel & Lupelli's Table 1 lists 18 JET probes as
  `( R, Z, angle )`. MEQ has no type for this.
* **Green's functions from a source to a sensor**, which is what makes EFIT's
  response matrix precomputable: *"The Green induction functions can be
  precomputed and stored in lookup tables."* MEQ has `CoilSet::psi()` and
  `flux()`, which are coil-to-point Green's functions already
  (`Coils.hpp:737-767`) — **so the coil half exists** and the plasma half, an
  integral of `G` over the plasma volume against a basis function, does not.

**(b) A MEASUREMENT AND UNCERTAINTY DATA MODEL.** No `sigma`, no channel
enable/disable, no file format, no notion of a time slice. EFIT's uncertainty
vector is not a detail — it is the mechanism by which 67 probes become 55.

**(c) A FIT AND REGULARISATION LAYER.** Part One's §4.3 corner solve is a
regularised least squares of order tens; this is one of order hundreds against
thousands, with an SVD truncation policy neither paper states. **`SurfaceFit`'s
truncation threshold and Levenberg–Marquardt damping are the nearest thing MEQ
has**, and `CLAUDE.md` already records that they *are* the gauge that makes the
gauge-free fit well posed — the same role, one layer up.

### 13.3 The Fitzgerald overlap, which cannot be assessed

`refs/Refs.md:297` records Fitzgerald, Appel & Hole (2013) as *"the only
implementation in the list and the only entry carrying both generalisations at
once"* — toroidal flow **and** anisotropic pressure, in EFIT. **MEQ has toroidal
flow and not anisotropy**, so the overlap is genuine and partial.

**The PDF is not on disk and I have not read it, so nothing further is claimed
about it here.** `doi:10.1088/0029-5515/53/11/113040` — and `Refs.md` records
that `TODO` transcribed the article number as 113408, which 404s. What **can** be
said from the DIII-D paper is that EFIT's rotating reconstruction (§V) adds a
third stream function and keeps the inner fit linear, and that MEQ's
`RotatingSource` computes the forward direction of the same physics.

---

## 13.4 THE IRON CORE: A CAPABILITY GAP, NAMED PRECISELY

**MEQ has no ferromagnetic material model whatever**, and the interesting part is
*which* assumption breaks, because the answer is not the obvious one.

### The obvious answer is wrong: the Grad–Shafranov operator is NOT modified

Appel & Lupelli state the equation once, for the whole domain, **iron included**,
and `mu` does not appear in it. Their (1) and (2):

```
Delta* psi_p = -2 mu0 R J_phi
J_phi = R dp/dpsi_p + ( 1/( mu0 R ) ) f df/dpsi_p + J_ext( R, Z )
J_ext = J_pf + J_induced + J_iron
```

**FLAG: the factor of 2 in (1) is as printed and is not MEQ's convention.** With
`psi_p = R A_phi` the conventional statement carries no 2. Either it is a
misprint or their `Delta*` or `psi_p` carries a factor the paper does not state.
**Do not port it.**

So the iron is **equivalent currents on the right-hand side**, on the same
footing as the poloidal-field coils, in an otherwise unchanged vacuum
formulation. Confirmed structurally by their solver: the field is solved *"using
the fast 2-D cyclic reduction finite difference method of Buneman"* — a
Cartesian-grid Poisson solver, which **cannot represent a jump in `mu`**. There
is no interface condition reaching the discretisation at all.

The currents come from the magnetisation `M = ( mu_r - 1 ) H`, `J_m = curl M`,
split three ways: `J_b` on the iron–air interface (13), `J_f = ( mu_r - 1 ) curl H`
in the bulk (15), and `J_mu = grad mu_r x H` (16). *"The magnetisation model
described in the current work includes the `J_mu` and `J_b` components but
excludes `J_f`. The assumption `J_f = 0` is reasonable for tokamaks that have a
laminated iron-core as is the case for JET."* And the axisymmetric reduction
leaves one scalar: *"the only finite components are `J_phi^b` and `J_phi^mu`."*

### So which of the three is it? — **(ii), extra source terms, with a caveat**

Extra source terms in an unchanged formulation: a **line** current density on the
iron–air interface and a **volume** current density on triangles inside the iron.
**Not** a modification of a Dirichlet-to-Neumann operator — the boundary integral
lives on the *physical* material interface, not on an artificial boundary, and its
output is fed back as a source. MEQ's exterior treatment would be structurally
unchanged.

**THE CAVEAT IS THE ONE THAT MATTERS TO MEQ, AND IT KILLS THE ECONOMY PART ONE
RESTS ON.** The iron's response is **nonlinear and field-dependent**:

* `mu_r` is a **tabulated B–H curve**, not a formula — *"based on the available
  data of the dominant material used for the JET iron core, Nomatil RD"*. There
  is no closed-form permeability model anywhere in the paper.
* The boundary-integral matrix and right-hand side both carry `mu_r` pointwise,
  so both **must be rebuilt whenever the field changes**.
* Saturation is observed, not hypothetical: *"the large corner current resulting
  in a local saturation with `mu_r ≈ 1`"*, and stripes of current appearing
  *"above a threshold magnitude of solenoid current |I_p1| ≳ 10 kA"*.
* The paper's own framing, §1: *"the presence of the iron core represents a
  complication by **adding an additional non-linearity** to the calculation of
  equilibrium force-balance."*

**Part One's whole cheapness argument is that a coil column is CONSTANT in the
iterate (§3.3). An iron column would not be.** It would have to be re-assembled
every Newton step, through a dense boundary-integral solve of its own. That is a
different and much more expensive object, and it is the honest reason iron is a
capability gap rather than an afternoon's work.

### What specifically breaks in `src/meq`, named

| assumption | where it lives | what iron does to it |
|---|---|---|
| **the exterior is `Delta*`-harmonic** | `src/meq/ExteriorDtN.hpp:117` and its header comment — the Gegenbauer separation assumes `Delta* psi = 0` outside `Gamma`, and the mode functions `rho^(1-n)` are the decaying solutions of the **source-free** equation | Iron **inside `Gamma`** breaks it outright. `ExteriorCoilSet::clearance()` (`Coils.hpp:924`) already names this failure for a conductor: *"psi_coil is then not `Delta*`-harmonic where the expansion assumes it is … the run converges to a machine nobody described."* **Iron inside `Gamma` is the same failure and nothing would detect it** |
| **an exterior conductor's current is KNOWN** | `meq::ExteriorCoilSet` — it sums `psi` and `flux` over conductors whose currents are data (`Coils.hpp:836-935`) | Magnetisation currents are **not known**; they are determined by the field. Outside `Gamma` they are `Delta*`-harmonic away from the iron, so they *could* be carried as **border unknowns** beside the Gegenbauer coefficients — but with an iterate-dependent column, above |
| **a conductor is a rectangle of UNIFORM current density** | `meq::Coil` (`Coils.hpp:303-355`) — *"a real winding has turns and this does not"*, and `CoilSet::f()` sums `mu0 r I/area` over rectangles | The iron needs a **line** density on a curve and a **non-uniform volume** density on a triangulation. Appel & Lupelli use **1442 piecewise-constant boundary elements and 4570 internal elements** at high resolution. `meq::Coil` cannot express either |
| **a conductor does not reach the axis** | `meq::Coil`'s constructor refuses `centreR - halfWidth <= 0` (`Coils.hpp:316`, mirrored in `Config.cpp:956`), because the operator's `1/r` is not integrable through `r = 0` | **A transformer core's central limb is ON the axis.** So MEQ cannot express JET's iron geometry even as a static conductor, before any magnetisation question arises. This is the sharpest single statement of the gap |

**AND THERE IS A MEASURED RESULT IN THAT PAPER THAT BEARS DIRECTLY ON PART ONE
AND ON §14, AND IT IS THE BEST THING IN EITHER OF THE NEW PDFs.** Appel &
Lupelli compare their **deterministic** model against the **fitting** model of
O'Brien *et al.* (1992), in which the boundary magnetisation currents are free
parameters fitted through quasi-measurements:

| | fitting model | deterministic model |
|---|---|---|
| free parameters | **110** (10 power supply + 2 plasma + 98 iron) | **12** (10 + 2) |
| `⟨chi²⟩` — the magnetic fit | **0.25** | **≈ 1** |
| probe residual | better | up to 15 mT / 5% |
| outboard separatrix vs HRTS | **~40 mm** | **15 mm** |
| outboard strike point vs IR | **~30 mm** | **0–10 mm** |

**THE MODEL WITH NINE TIMES FEWER PARAMETERS FITS THE MAGNETICS WORSE AND GETS
THE PHYSICAL TOPOLOGY RIGHT.** The paper's own conclusion, §7: *"the large number
of free parameters in the fitting magnetisation model have the possibility of
absorbing inconsistencies in the data making systematic errors difficult to
identify."*

**That is Lackner's figures 3–4 (§1.4), measured on a real machine in 2018**, and
it is the same statement: extra degrees of freedom that the data cannot see get
set by noise, the misfit improves, and the answer gets worse. **It is evidence
for Part One's §4.5 rank reporting and against ever accepting a fit residual as
the acceptance criterion.**

---

## 14. THE HONEST OBSTACLE: THIS IS DATA-LIMITED, NOT DISCRETISATION-LIMITED

**Asked whether high-order HDG buys anything for a problem whose data is sparse
and noisy, the answer this reading supports is: very little, and the numbers are
in the papers.**

### 14.1 How many parameters the magnetics actually support

DIII-D §IV.A, on the production magnetics-only parametrisation:

> "**There are five parameters. Three of them describe `β_p`, `ℓ_i`, and `I_p`,
> and the remaining two are determined by the two boundary conditions on
> `P′(ψ₁)` and `FF′(ψ₁)`.**"

Appel & Lupelli §6, independently:

> "the plasma parameterisations used in the equilibrium calculations describe the
> plasma in terms of **just two degrees of freedom**"

**Five coefficients, of which three carry information — or two, on JET.** MEQ
resolves `psi` at `k+1` over thousands of dofs. **The mismatch is three orders of
magnitude**, and no amount of discretisation accuracy adds a sixth parameter the
magnetics can constrain.

### 14.2 And the failure of magnetics-only is a TOPOLOGY error, not a magnitude one

DIII-D §IV.B, the measured demonstration on shot 92043:

> "As expected, reconstruction using external magnetic data does not provide
> accurate information on the `P` and `q` profiles. **`q₀` and the minimum `q`
> value `q_min`, from the magnetically reconstructed case, are both 0.76, whereas
> the values from the full reconstruction case are 1.23 and 1.01.**"

From Table I, the same four reconstructions:

| | magnetic | + MSE | + kinetic | + rotation |
|---|---|---|---|---|
| `I_p` (MA) | 1.482 | 1.481 | 1.486 | 1.492 |
| `kappa` | 1.84 | 1.86 | 1.87 | 1.87 |
| `beta_p` | 1.63 | 1.65 | 1.62 | 1.56 |
| **`q_0`** | **0.76** | **1.13** | **1.23** | **1.32** |
| **`q_min`** | **0.76** | **1.00** | **1.01** | **1.00** |

**The shape and the global parameters agree to a per cent or two across all four;
`q_0` moves by 74%.** And `q_0 = q_min` in the magnetics-only column means it
returned a **monotonic** `q` where the truth is flat or reversed — the wrong
topology, not a wrong number. §I states the limit outright: external magnetics
*"alone can only yield global current and pressure profile parameters"* and *"to
fully reconstruct the pressure and current profiles … internal current and
kinetic profile measurements must be used."*

**MEQ cannot improve that column by discretising better.** It is an information
limit, and the remedy is MSE and kinetic data — 34 channels of it.

### 14.3 What high order might still buy, stated fairly

Three things, none of them the headline, and **none measured**:

* **Fewer unknowns for the same forward accuracy.** EFIT runs 65×65 between shots
  at ~1.8 CPU-s per reconstruction, 129×129 for stability analysis. If `k = 3`
  on a coarse mesh matches 129×129, the *inner* solve is cheaper — and the inner
  solve is already not the bottleneck at 250 reconstructions in 25 seconds over
  24 CPUs.
* **`psi_0` and `psi_1` from a solved `q`.** These define `x` and therefore the
  entire frozen map of §12.2, and EFIT finds them by searching a grid. An `O(h)`
  error in the axis position perturbs every basis function `y_n(x)`. **Whether
  that matters at the data's noise level is exactly the unanswered question**, and
  it is the one experiment in this half worth naming: reconstruct the same
  synthetic dataset with the axis located to `O(h^{k+1})` and to `O(h)`, and see
  whether the fitted coefficients move by more than their own uncertainty.
* **The separatrix and strike points.** Appel & Lupelli measure model differences
  at **15 mm, 30 mm, 6 mm** — and say the HRTS and IR uncertainties are
  comparable. **A discretisation that resolves the separatrix to microns is
  resolving far below the measurement.** This cuts against high order rather
  than for it.

**The honest summary: MEQ's high order is worth having for the forward problem
inside a reconstruction and is not what makes a reconstruction good.** Saying so
is more useful than the alternative, because it says where the effort would
actually go — the synthetic diagnostics and the fit layer of §13.2, neither of
which is a discretisation question.

---

## 15. Is the `q(psi)`-driven solve a partial precedent? — **Yes, and a closer one than it looks**

`CLAUDE_INVERSION.md`, *Driving by `q(ψ)`*, and `src/meq/SafetyFactorSolve.{hpp,cpp}`.

**What matches, point for point:**

| reconstruction | the `q(psi)` solve |
|---|---|
| infer profile coefficients from data | infer `g²` coefficients from a target `q(psi)` |
| the data are functionals of the whole solution | `V'` and `⟨R⁻²⟩` are functionals of the whole solution |
| the profile is **fitted**, few coefficients | *"`g²` is carried as a handful of least-squares coefficients rather than a table, so a differenced Jacobian costs a few map evaluations a step"* |
| an outer loop around a forward solve | KINSOL outer Newton around the bordered Newton |
| one solver serves every evaluation | *"Only `gg′` changes between steps, so the mesh, the spaces, the forms and the trace solver's symbolic factorisation all survive"* |

**It is a reconstruction with one measurement type, no noise and no uncertainty
model** — and every structural lesson it learned transfers:

* **The fit is what makes the outer Jacobian affordable**, and it was put there
  for conditioning: *"the innermost surface's inverted `g` moves 13% for a 1.3%
  change in the profile."* That is an ill-conditioning statement of exactly
  §14's kind, met and handled.
* **A converged outer solve is not evidence the root is the answer.** Both
  failing configurations converged to 3.5e-07 on a fixed point that was not the
  answer. **In a reconstruction, with noise, this is the normal case rather than
  a bug** — and it is why §13.2(b)'s uncertainty model is not optional.
* **The knot span must cover the whole of `psi`, not the data's range.** Laying
  knots only over `[0.05, 0.95]` left *"the outer residual AT THE ANSWER 1.3e-01
  rather than zero"*, and both methods then converged to the same wrong answer to
  six digits. **EFIT's edge boundary conditions (Eqs. 19–23) are the same
  problem**, and DIII-D §IV.A says getting them wrong is not subtle: *"Forcing
  `P′(ψ₁)`, `FF′(ψ₁) = 0` in H-mode configurations often leads to an inaccurate
  reconstruction of the plasma boundary."*

**Where it is NOT a precedent**: it is an outer loop, which §4.2 argues against
for the inverse solve — and for a reconstruction an outer loop may well be
right, because EFIT's own architecture puts *locations* (knots, channel centres)
in an outer nonlinear driver and only *amplitudes* in the inner linear fit
(§12.4). **The `q(psi)` solve is the outer driver; what MEQ does not have is
EFIT's cheap inner linear fit.**

---

## 16. Open questions, and what I did NOT resolve

**Part Two ends here deliberately. These are the questions, not a stage list.**

1. **Does the frozen-flux-map linearisation survive HDG?** EFIT's inner fit is
   linear because `psi^m` freezes the basis fields `y_n(x)`. In MEQ the profile
   reaches the residual through `meq::SourceIntegrator` at quadrature points, so
   the same freezing is available in principle — **and nobody has written down
   what `dR/d alpha_n` is at fixed `psi`, or whether it is as cheap as
   `assembleCurrentColumn()`.** This is the first thing to settle and it decides
   whether a reconstruction is a border problem or an outer loop.
2. **Is the exterior field at a sensor position accurate enough?** `ExteriorDtN`
   can evaluate it; nothing has measured it against anything. Most sensors are
   outside `Gamma`.
3. **What is EFIT's SVD truncation policy?** Not in any of the three PDFs on
   disk. It is the operative regularisation decision in the method.
4. **How many Newton/Picard iterations does a reconstruction take?** Neither
   paper says. It is the number that would size the work.
5. **Does the axis-location accuracy matter at the data's noise level?** §14.3.
   The one experiment in this half that would settle whether MEQ's distinctive
   property is worth anything here.
6. **Anisotropic pressure** — unassessable; the PDF is absent (§13.3).
7. **Inequality constraints** — coil and force limits, deferred out of Part One
   in §4.5 and unresolved for both halves.
8. **Gap constraints** — §5's last paragraph. A moving isoflux point does **not**
   get the envelope theorem, and what its row costs is not worked out.
9. **`Gamma = 0` meaning "derive it"** — §7.3, a schema question left open.
10. **Whether the `.nc` is the right place for the inverse solve's report**, or
    whether a fit wants a format of its own.

### What I did not verify

* **I did not open any PDF directly.** All four readings are through
  `pdftotext` extractions, three of them by delegated readers whose reports flag
  their own mangled passages. The DIII-D PDF's maths font decodes to Latin
  letters and required a substitution key; **the Lackner PDF is a scanned 1976
  page whose displayed equations (17) and (18) are reconstructions, not
  transcriptions**, and are labelled as such in §1.4's source. **The prose of
  Lackner §2.3 is intact and the quotes in §1.4 are verbatim; the equations are
  not.**
* **`refs/Refs.md:492`'s summary of Lackner §2.3 agrees with the paper** on every
  point checked — Hadamard, Fourier truncation, Zakharov/Tikhonov, and the
  figures 3–4 result. **Where `Refs.md` is imprecise is CEDRES++**: "four problem
  statements" is four *modes* and eight *problems*, and "Tikhonov" is the right
  mathematics under a name the paper never uses. **Trust the paper.**
* **I have not run anything.** No number in this file was produced here.
* **`GradShafranov.{hpp,cpp}` were under uncommitted edit while this was
  written.** Line numbers may have moved; the names have not.

---

## 17. `refs/Refs.md` needs three entries, and this file does not add them

Per `CLAUDE.md`, `refs/Refs.md` is the index and is written by hand. **None of
the three reconstruction PDFs on disk has an entry**, and a PDF in `refs/` that
the index does not know about is the condition this project's own rules exist to
prevent. Recommended, with the dois verified as far as the documents themselves
allow:

| file | entry |
|---|---|
| `MHD Equilibrium Reconstruction in the DIII-D Tokamak.pdf` | Lao, St. John, Peng, Ferron, Strait, Taylor, Meyer, Zhang & You, *Fusion Science and Technology* **48**:2 (2005) 968–977, **`10.13182/FST48-968`**. The EFIT reference. Note its `F = 2πRB_φ/μ₀` convention, which puts `μ₀` in the numerator of the `FF′` term |
| `EFIT2006.pdf` | Appel, Huysmans, Lao, McCarthy, Muir, Solano, Storrs, Taylor, Zwingmann + EFDA ITM Task Force + JET-EFDA Contributors, *A Unified Approach to Equilibrium Reconstruction*, **EFDA–JET–CP(06)03-07**, 33rd EPS Rome 2006. **Architecture and data plumbing, not physics** — the entry should say so, since the title promises otherwise. **Its reference [1] miscites Lao 1985 as Nucl. Fusion vol 22; it is vol 25** |
| `EFIT-revisited-for-MAST.pdf` | Appel & Lupelli, *Equilibrium reconstruction in an iron core tokamak using a deterministic magnetisation model*, *Comput. Phys. Commun.* **223** (2018) 1–17, **`10.1016/j.cpc.2017.09.016`**. **THE FILENAME IS WRONG — the machine is JET, not MAST.** Consider renaming. Carries the 12-vs-110 parameter result of §13.4 and the `Δ*ψ = −2μ₀RJ_φ` factor-of-2 flag |

and `refs/Refs.md:297`'s Fitzgerald row should record that **the PDF is still not
on disk**, since it currently reads as an unread seed without saying whether it
was ever fetched.
