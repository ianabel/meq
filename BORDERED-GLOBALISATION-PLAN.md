# Globalising the bordered Newton

Written 2026-09-15. `CLAUDE.md` is the maintainer's index, `CLAUDE_HDGGS.md`
owns Newton and the globalisation ladder, and `CLAUDE_FB.md` owns the borders;
this file is the design for one change that sits across all three and defers to
each of them on anything they already say. `FREE-BOUNDARY-PLAN.md` §7 and §10
are the stage record this one continues from, and `MEASUREMENTS.md` M-103 is the
measurement that makes it worth writing.

**Line numbers in `src/meq/GradShafranov.cpp` and `apps/meq.cpp` are as read on
2026-09-15, and `GradShafranov.{hpp,cpp}` are under uncommitted edit at that
moment** — the working tree deletes `NonlinearOrdering` entirely, so NPC is the
only ordering and every "under the condensation" branch below is dead code
waiting to be removed. Every citation also names the function or the lambda,
which is the durable half; `grep` for the name rather than trusting the number.

**The change in one sentence.** Make a globalisation reachable when `ψ_ax` is an
unknown of the bordered Newton, so that the four of seven diverted machines that
do not converge from a cold start have a second rung to fall to — and do it
without letting a rung change which equilibrium a converging run reports.

---

## 0. WHAT WOULD FALSIFY THIS, AND THE CHEAPEST EXPERIMENT

Put first because it decides whether the rest is worth reading, and because in
this case the first experiment may kill the whole route.

### 0.1 The gate: it is NOT established that these are path failures

Six of M-103's fourteen configurations exit 2 with

```
no damping of the bordered Newton step gave a finite residual -- psi_ax
through zero, most often, which is the branch leaving the physical one
```

(`GradShafranov.cpp:7115`). That message is reached only when `bestDamping` is
still `0.0` after **twelve halvings** (`GradShafranov.cpp:7013-7018`), i.e. when
*every* trial down to `α = 2⁻¹¹` produced a non-finite `trialNorm`
(`GradShafranov.cpp:7086-7106`).

**That is a strange thing for a step-length failure to be.** The trial state is
`savedState + α·(direction)` and `savedState`'s own augmented norm was finite at
the top of the same iteration (`GradShafranov.cpp:6358-6372`), so as `α → 0` the
trial tends to a state whose norm is known finite. Three things can produce a
non-finite trial at `α = 2⁻¹¹` and they want three different repairs:

| cause | what it is | does a globalisation help? |
|---|---|---|
| **a non-finite direction** | `y` or `z` from the trace solve carries `NaN`/`Inf`, so `unknown.Add( -damping, y )` is `NaN` at every `α`. Only the *border* solution is checked for finiteness (`GradShafranov.cpp:6959-6961`); `y`, `z`, `zB`, `zL` and the exterior columns are not | **No.** The Jacobian or its elimination is the defect |
| **a throwing constraint** | `refreshXPoint()` carrying the point out of the mesh, `peakAt()`/`locateAxisPoint()`, or `setNormalisation()` refusing a `ψ_ax` through zero. All are caught and turned into `trialNorm = ∞` (`GradShafranov.cpp:7097-7106`) | **Only if** the throw is caused by the step. If the *saved* state already throws, damping cannot rescue it |
| **a genuine excursion** | the direction is finite, the constraints evaluate, and the augmented norm is `Inf` because the residual overflowed | **Yes.** This is the case a globalisation is for |

**Nothing in the code distinguishes them and nothing has measured which it is.**
The message asserts the third ("psi_ax through zero, most often") and the code
cannot know that.

**The experiment.** In the `!accepted` branch (`GradShafranov.cpp:7112-7116`),
record per trial which of the three fired — `std::isfinite` over `y`, `z` and
every queued column before the loop; the exception's `what()` when the catch
fires; and `trialNorm` when neither. Re-run the six exit-2 configurations
(`machine-b-ffprime`, `machine-b-ffprime-shaped`,
`machine-e-diamagnetic`, `machine-e-diamagnetic-shaped`, `machine-d-tcv`,
`machine-c-mast-shaped`) and tabulate. It is one diagnostic block, six runs and
no new algorithm.

**If the answer is "a non-finite direction" on most of the six, this whole plan
is the wrong plan** and what is wanted is an investigation of the bordered
elimination — most likely `M_ij = D_ij − b_i·z_j` going singular as the X-point
row and the `ψ_bnd` row collapse onto each other, which the dense inverse
(`GradShafranov.cpp:6955-6958`) does not detect until after it has produced
`Inf`.

### 0.2 The control that already exists and has NEVER been run on these cases

`[solver] PicardSweeps` is family (b) below, built and committed in `c12501d`
under the title *"a basin-finder that does not yet work"*. **No example, no
test and no benchmark script in the tree sets it** — `grep -rl PicardSweeps`
finds it only in `Config.{hpp,cpp}` and `apps/meq.cpp` — and M-103 was taken on
the shipped configurations with "nothing tuned per machine", so **every row of
M-103 was measured at `PicardSweeps = 0`.**

The commit message records it measured on **case A alone**, where it fails in
two specific ways (§4.3). Case A is one of the three rows that converges anyway
and is a *wrong-root* failure rather than an exit-2 one, so the existing
measurement is against the failure mode the pre-stage is least likely to fix.

**The experiment.** Run the six exit-2 configurations at
`PicardSweeps ∈ { 2, 4, 8 }`, `PicardBlend = 0.5`, everything else shipped. It
costs fourteen `meq-run` invocations and no code at all. **If two or more of the
six close, family (b) is already the answer and family (a) is redundant work.**

### 0.3 The argument that survives both

**The refusal is stated as a structural impossibility and it is not one.**
Both throws say the Picard paths *cannot* be bordered:

```
psi_ax as an unknown is implemented for Globalisation::None only -- the
KINSOL paths drive a residual of their own and the Picard ones do not build
a Jacobian at all                                   GradShafranov.cpp:5053

... and the Picard ones build no Jacobian to border GradShafranov.cpp:7292
```

§2.3 shows that **three of the four bordered Jacobian pieces need no Newton
Jacobian at all** — the row is exact by the envelope theorem, the corner is
exactly `1`, and the column is assembled from the source's own
`normalisationDerivatives()`. What a Picard path lacks is not the border, it is
the *operator MEQ borders*. That is a much weaker obstacle and it is removable.

So even if §0.1 and §0.2 both come back negative, the comment is wrong and
should be corrected; and if either comes back positive, §2.3 is the reason the
work is small.

---

## 1. The problem

### 1.1 The two throws

`GradShafranovSolver::solve()` refuses every globalisation but `None` once a
`NormalisedSource` is set (`GradShafranov.cpp:7291-7292`), and
`solveWithNormalisation()` refuses it again on entry
(`GradShafranov.cpp:5052-5053`). The second is unreachable through the first and
is there as a contract on a function that can also be entered for a pure
exterior-coupled solve (`GradShafranov.cpp:7347-7350`).

The driver knows it: `apps/meq.cpp:2751-2763` catches the failure, prints
*"no globalisation MEQ has can carry that border — so there is no fallback to
try"*, and returns `SolveFailed` **without attempting a second rung**, where the
unnormalised path falls to `Globalisation::PicardThenNewton`
(`apps/meq.cpp:2765-2790`).

### 1.2 M-103, sorted by failure mode rather than by machine

M-103's table is by machine and conductor model. Sorted by what actually
happened, its fourteen configurations are:

| n | outcome | configurations |
|---|---|---|
| 4 | **converged and agrees** | F filament, F shaped, D shaped, G filament |
| 6 | **exit 2**, no finite damped step | B ×2, E ×2, D filament, C shaped |
| 3 | **exit 0, no O-point anywhere** | A ×2, G shaped |
| 1 | **exit 1**, the plasma contains the symmetry axis | C filament |

**Only the six exit-2 rows are candidates for a globalisation.** The other four
failures are converged solves of the system as posed — a globalisation changes
the path, and a converged solve has already finished travelling.

The scale of what is at stake, from M-103: DIII-D closes on both arms and
reproduces M-96 to every digit (4.593e-03 against the filament reference,
9.102e-04 against the shaped one); the three no-axis rows carry 3.060e-01,
4.146e-01 and **1.440e+01** relative `L2`; and G under shaped conductors spends
**187 s of solve**, exits 0, and reports that 1.44e+01 — *"A run that used to
fail now reports a non-equilibrium as a success, which is strictly the worse of
the two."*

### 1.3 Three of the fourteen are not path failures, and the mechanism is known

**The axis border silently decouples when no axis can be found, and then the
solve converges on a system with one fewer real constraint.**

`peakAt()` (`GradShafranov.cpp:5872-5905`) tries `locateAxisPoint()` under
`AxisConstraint::LocatedAxis` and, failing, falls through to the nodal maximum
over the potential block with `constraintLocated` left **false**. The border row
is then built by the `else` branch at `GradShafranov.cpp:6428-6494`, which — for
`argDof < 0` — sets an **empty row**:

```
if ( argDof < 0 )
{
	borderDofs.SetSize( 0 );
	border.SetSize( 0 );
}
```

`rowDot()` returns exactly `0.0` for an empty row, so `b = 0`, the corner is
`1`, and the dense solve degenerates to `step[ 0 ] = −constraint`, i.e.
`ψ_ax ← max ψ_h`. **That is a Picard update, not a Newton one**, and the
comment beside it says so.

On a normal solve this fires at iteration 0 only, because the cold reference
evaluates the constraint on a state whose flux and potential blocks are zeroed
(`GradShafranov.cpp:6203-6280`) — `located 0, rowSize 0` at iteration 0 and
`located 1, rowSize 6` after it, measured on `examples/diverted-tokamak.toml`.
**On a run with no interior extremum of `ψ_h` anywhere it fires at every
iteration.** The equation being converged is then "`ψ_ax` equals the largest
nodal value of `ψ_h`", which a wall-hugging annulus and an `r = 0` boundary
layer both satisfy — both of them documented failure shapes in `CLAUDE_FB.md`,
*`ψ_ax` was the weak point of the whole free-boundary path*.

**No globalisation can touch this.** The iteration converges; it converges to a
solution of a degenerate system. The repair is a refusal or a different
constraint, and it belongs to `FREE-BOUNDARY-PLAN.md` §11 rather than here.
§6.1 says what the shape of it would be.

### 1.4 What a globalisation is allowed to change

`CLAUDE.md`'s standing rule and `MEASUREMENTS.md` M-26 together set the bar.
M-26 has three solve routes reaching discrete solutions **9.4% apart in
`max ψ_h`** on an under-resolved mesh, all three converged to `rel_tol = 1e-12`.
So:

* a globalisation is **not** a performance knob and must not become a default;
* it must be **reactive** — tried only after an observed failure, never chosen
  from a property of `F` or of `‖r₀‖`;
* and the run must **say which rung produced the answer**, because a reader
  cannot otherwise tell M-26's hazard from a converged equilibrium.

`Globalisation` is deliberately not a TOML key for exactly this reason, and `PicardSweeps` is opt-in for exactly this reason
(`Config.hpp:1000-1001`: *"a key that silently changed which equilibrium those
runs report would be exactly the kind of thing CLAUDE.md refuses to expose"*).
Nothing below may weaken that.

**And there is no predictive trigger to be had.** Two candidates are measured
and dead: `max|∂F/∂ψ|/λ₁` (`CLAUDE_HDGGS.md`, *Why MEQ's Newton struggles*; two
points, no threshold, and it needs the range of `ψ` which is not known before
solving) and `‖r₁‖/‖r₀‖`, which M-24 measures **anti**-correlated. **A third
observation, from this tree's own shipped outputs, says the iteration count
carries no signal either**: `machine-a-testtokamak.nc` and
`machine-a-testtokamak-shaped.nc` both record `newton_iterations = 2` and no
`axis_normalised_flux` at all, while `machine-f-diiid.nc` — the one that is
right — also records `newton_iterations = 2` with
`axis_normalised_flux = 1.00000000002119`. The wrong root and the right one
take the same number of steps. **That reading is an observation and not a
measurement** — those three files carry `meq_version` 969b626, 969b626 and
d4a203b, so they are not one run of one tree, and the count is the last support
sweep's rather than the solve's. It is worth re-taking inside BG-0, where the
instrumentation is being added anyway. Do not propose a fourth detector without
saying how it would be falsified on all fourteen of M-103's rows.

---

## 2. What the bordered path already has, and what it structurally cannot have

### 2.1 There is already a globalisation in there, and it is an Armijo backtrack

`solveWithNormalisation()` carries a backtracking line search on the augmented
norm (`GradShafranov.cpp:6998-7130`): twelve halvings from `α = 1`, accepted on
`trialNorm < ( 1 − 1e-4·α )·norm`, with the best non-improving trial taken as a
fallback. It damps **every** unknown together — the fields, the trace, `ψ_ax`,
`ψ_bnd`, the plasma-current scale, the exterior coefficients and the X-point's
two coordinates — which is the half of a line search a trace-only operator
cannot express.

So the question is not "should the bordered path have a globalisation"; it has
one. The question is what to do when a *step-length* rule is not enough, which
is the same question `CLAUDE_HDGGS.md`, *Why it fails, measured — and why a line
search does not fix it* already answered for the unbordered path: **what works
is fixing where the iterate is, not how far it steps.**

### 2.2 The NPC Jacobian is solve-only, so most textbook globalisations are out

`mfem::DarcyNPCOperator::Jacobian::Mult` **aborts**
(`../mfem/install/include/mfem/fem/darcy/darcyhybridization.hpp:3756-3761`):
after `ComputeH()` the local arrays hold factored blocks, so the Jacobian can be
solved with and cannot be applied. There is no `Jᵀ` either.

That removes, structurally:

* a **merit-function gradient** `∇‖R‖² = 2JᵀR`, hence any Armijo test with a
  real directional derivative rather than the `(1 − 1e-4α)` surrogate that is
  there;
* a **dogleg** or any trust region needing a Cauchy point;
* **Newton–Krylov / JFNK**, which needs `J·v`;
* and anything that wants to *compare* `J` between iterates.

What is available is: the step direction, the residual and the constraints at
trial points, the dense `(N+4)` border system, and the freedom to change what
operator is assembled. Three of the four families below therefore change the
operator or the problem rather than the step rule, and that is forced rather
than chosen.

### 2.3 The borders do not need a Newton Jacobian, and the throw's reason is wrong

This is the load-bearing observation of the whole document. Taking the four
pieces of the bordered system in turn, at the state `x` with border unknowns
`p`:

| piece | how it is obtained today | does it need `∂F/∂ψ`? |
|---|---|---|
| **row `b`**, `−∂ψ_ax/∂x` | the potential shape functions of `x*`'s element, exact and undifferenced, because `∇ψ_h( x* ) = 0` at a zero of `q_h` — the envelope theorem (`GradShafranov.cpp:6391-6427`) | **No.** It needs only that `ψ` is an unknown of whatever linear system is being solved |
| **corner `d`**, `∂G/∂ψ_ax` | exactly `1` under NPC, `G = s − ψ_j` and `s` does not appear in `ψ` (`GradShafranov.cpp:6379-6389`) | **No** |
| **column `c`**, `∂R/∂ψ_ax` | `assembleNormalisationColumn()` (`GradShafranov.cpp:4951`) runs `SourceIntegrator`'s own quadrature loop over `∂F/∂s` from `NormalisedSource::normalisationDerivatives()`, into the potential block and nowhere else | **No.** It is `∂F/∂ψ_ax` at fixed `ψ`, a different derivative |
| **the field block `J`** | `npc->GetGradient( unknown )` (`GradShafranov.cpp:6613`) | **Yes** |

The same holds for the other borders: the plasma-current row, column and corner
are all analytic (`assembleCurrentRow`, `assembleCurrentColumn`,
`assembleCurrentNormalisationCorner`); the exterior rows are transmission
integrals and their columns are *constant* in the iterate; and XP-3's two rows
are functions of `q_h` alone with a zero column
(`GradShafranov.cpp:6828-6842`).

**So exactly one of the four pieces needs the Newton Jacobian, and it is the one
a Picard method replaces on purpose.** The throw's sentence — *"the Picard ones
build no Jacobian to border"* — is true of the field block and false of the
border, and the design below turns on that distinction.

### 2.4 A border-Picard already exists and is already measured failing

`Normalisation::Decoupled` (`setNormalisationCoupling()`) zeroes the column
(`sourceColumn`, `GradShafranov.cpp:5934-5940`) and the row (`coupled` false,
`GradShafranov.cpp:6425`), leaving `step[0] = −constraint` — a Picard update on
`ψ_ax` wrapped round a full Newton on the field. **M-31 measures it**: coupled
reaches `4.447e-15` in 4 iterations, decoupled sits at `8.239e-02` after 15 and
does not converge.

That is the *opposite* split from family (a): Newton on the field, Picard on the
border. It is worth stating because it is the obvious first idea and it is
already dead, and because it is the control that says the coupling is doing real
work. **Do not re-derive it.**

---

## 3. Family (a): a genuinely bordered Picard

### 3.1 The mechanism

Keep the bordered elimination, every border, the Armijo loop and the augmented
norm exactly as they are. Change only what the **field block** is:

```
Newton    J = A_lin − (∂F/∂ψ)/r mass term,   residual R( x ) with F( r, z, ψ_h )
Picard    J = A_lin,                          residual R( x ) with F( r, z, ψ^k )
```

where `ψ^k` is the last accepted iterate's potential and `A_lin` is the
hybridized HDG operator with no reaction term. Each step is then **Picard in the
field and Newton in the borders, solved simultaneously** rather than in
alternation — which is what distinguishes it from `Normalisation::Decoupled`
(§2.4) and from the driver's outer sweep (§4).

The fixed point it drives is the papers' own — `CLAUDE_HDGGS.md`, *Picard,
keeping the local problems linear* — and it inherits that method's property:
undamped Picard stalls and needs relaxation, Anderson at depth 1 converges and
depth 2 fails. See §3.6 on whether Anderson is reachable here.

### 3.2 It goes in `SourceIntegrator`, NOT in `buildForms()`

The obvious implementation is the existing one: `usesNonlinearForms()` returns
false, `buildForms()` puts the potential block on the **linear** form
(`GradShafranov.cpp:3447-3487`), and `prepare()` assembles `FrozenSource` onto
the right-hand side (`GradShafranov.cpp:3809-3828`). **That route is wrong here,
for four independent reasons:**

1. `solveWithNormalisation()` builds an `mfem::DarcyNPCOperator` from
   `darcy->GetHybridization()` (`GradShafranov.cpp:5133-5136`). With the
   potential block on the linear form, `EnableHybridization()` takes the `M_p`
   branch (`darcyform.cpp:367`) rather than the `Mnl_p` one, and whether
   `DarcyNPCOperator::Mult` is even defined on that configuration is **not
   established**. It would have to be read out of `darcyhybridization.cpp`
   before anything is built on it.
2. The exterior coupling is refused outright when `usesNonlinearForms()` is
   false (`GradShafranov.cpp:5078-5084`), and every machine case in M-103 has
   one.
3. `FrozenSource::Eval` (`GradShafranov.cpp:788-794`) calls `source.f()`
   directly. It does **not** consult the element mask, so XP-1's connectivity
   fill is inert on that path and `fOutsidePlasma()` — the thing that keeps the
   coils alive outside the plasma component — is never reached. That is the
   same defect shape as §4.4 below.
4. `buildForms()` **replaces the `DarcyForm` outright**
   (`GradShafranov.cpp:3394`), which is fatal mid-solve. See §3.5.

**The right mechanism is one flag and one pointer on `meq::SourceIntegrator`**,
which MEQ owns:

```
void setFrozenPotential( mfem::GridFunction const *frozen );
```

With it set, `AssembleElementVector` evaluates `F` at `frozen`'s value at the
quadrature point instead of at `shape*elfun`, and `AssembleElementGrad` returns
zero for the reaction term. Nothing else moves: the source stays on the
non-linear potential mass form, `usesNonlinearForms()` stays true, the element
mask and `fOutsidePlasma()` keep working through `sourceValue()`
(`GradShafranov.cpp:340-351`), `DarcyNPCOperator` is used unchanged, the
exterior guard is untouched, and **no form is rebuilt**.

### 3.3 What each border piece costs under it

Unchanged, item for item:

* `b` — the shape functions at `x*`. `locateAxisPoint()` roots `q_h`, which is
  still a solved field under the frozen source. Exact.
* `d` — still exactly `1`.
* `c` — `assembleNormalisationColumn()` evaluates `∂F/∂ψ_ax` **at the state it
  is handed**. Under the frozen mechanism it must be handed the *frozen*
  potential, or the column will be the derivative of a source the residual is
  not using. That is a one-line consistency requirement and it is exactly the
  trap the code already guards for the plasma support
  (`GradShafranov.cpp:6604-6611`: *"a Jacobian assembled on a different support
  from its residual is the defect this guards"*).
* the current, exterior and X-point pieces — unchanged.

### 3.4 The Jacobian is constant, and whether that can be cashed is OPEN

With `∂F/∂ψ` dropped, the field block depends on the mesh, `τ`, `r` and the
extension and on **nothing that moves between iterations** — `F` is on the
residual only. Mathematically that is one factorisation for the whole Picard
phase, against one per Newton step.

**Whether MFEM can be told that is not established.** `GetGradient()` runs
`ComputeH()` unconditionally, and `ComputeH()` is the item M-80 and
`DarcyHybridization::GetComputeHTime()` exist to measure. There is a
condensation cache (`CanCacheCondensation()`, refused under
`LocalFactorMode::Batched`) and `SetReuseSymbolic()` already keeps the symbolic
factorisation across Newton steps. **Do not assume the prize.** If it cannot be
cashed, a bordered Picard iteration costs about what a bordered Newton
iteration costs, and the route buys robustness only — which is what M-36
already says Picard is for: *"a robustness route, not a faster one"*.

### 3.5 `usesNonlinearForms()` and `built`, and why this is worse here

`setGlobalisation()` resets **`built`** as well as `prepared`
(`GradShafranov.cpp:1223-1233`), which is the fix for the latent defect
`CLAUDE_HDGGS.md`, *Picard, then Newton* records. On the unbordered path that is
sufficient, because `solveByPicardThenNewton()` switches globalisation
**between** two calls to `solve()` (`GradShafranov.cpp:3907`, `3923`).

**On the bordered path a mid-solve switch is not merely wrong, it is undefined
behaviour.** `buildForms()` does

```
darcy = std::make_unique<mfem::DarcyForm>( fluxFes.get(), potentialFes.get() );
```

(`GradShafranov.cpp:3394`), destroying the `DarcyForm` that the live `npc`,
`reduced`, `traceX`, `traceB`, `darcySolution`, `darcyRhs` and the
`TimedSolver`'s bound operator all refer to. `reprepare()`'s own comment
(`GradShafranov.cpp:5145-5155`) records what happens when only the
*hybridization handle* goes stale: *"the process dies inside this function with
nothing in the trace to say so."* Replacing the whole form is strictly worse.

**Two design rules follow and both are absolute:**

* **A globalisation may never be switched inside `solveWithNormalisation()`.**
  The Picard-to-Newton handoff of §3.7 flips a local flag that `SourceIntegrator`
  reads; it does not call `setGlobalisation()` and it does not touch `built`.
  This is a second, independent argument for the mechanism of §3.2.
* **`usesNonlinearForms()` has to learn about the normalised source.** Today it
  reads

  ```
  return globalisationChoice != Globalisation::AndersonPicard
      && globalisationChoice != Globalisation::PicardOnly;
  ```

  (`GradShafranov.cpp:3379-3385`). If a bordered `PicardOnly` is to keep the
  source on the non-linear form, that predicate must become
  `normalisedSource != nullptr || ( the two comparisons )`. **That is not a
  cosmetic edit**: `usesNonlinearForms()` is read by `buildForms()`
  (`:3447`), by `solve()`'s Picard dispatch (`:7317`), by `fieldsAreState`
  (`:7326`) and by the exterior-coupling guard (`:5078`), and it decides which
  branch `EnableHybridization()` takes — which decides whether MEQ's face
  constraint goes on `c_bfi_p` or `c_nlfi_p`, i.e. whether the batched face
  route is reachable at all (`CLAUDE.md`, *MEQ is in the 15 that the batched
  face kernel takes*). Any change to it needs the bit-exactness cases in
  `SolverContract.cpp` re-run.

### 3.6 Anderson, KINSOL, and a rule this tree keeps

M-36 and `CLAUDE_HDGGS.md` both say plain Picard needs damping and Anderson at
depth 1 does not. MEQ gets Anderson from `mfem::KINSolver( KIN_FP )`
(`solveByPicard()`, `GradShafranov.cpp:1309-1370`).

**KINSOL cannot be handed this fixed point.** `KINSolver::Mult` drives
`oper( x ) = 0` or `u = G( u )` over one vector; the bordered iteration's state
is `(x, p)` with `p`'s equations being the constraints, and MEQ's residual is
not zero (hence `ShiftedResidual`). Wiring it would re-create exactly the
"drives a residual of its own" problem the throw complains about.

So a bordered Picard is a loop inside `solveWithNormalisation()`, and Anderson —
if wanted — is a small dense least squares over the last `m` increments of
`(x, p)`. **That collides with `CLAUDE.md`'s standing preference for a
maintained library over a hand-rolled algorithm**, and the collision should be
stated rather than resolved by taste: KINSOL is the maintained implementation
and it structurally cannot take this problem; Eigen is already a sanctioned
dependency for dense least squares and would carry the QR. **Recommendation:
build the damped fixed point first with no acceleration (BG-3), measure it, and
only add Anderson if the iteration count is the thing that kills it** — M-36
reads 122 to 290 Picard iterations on the unbordered stiff cases, and at M-90's
1.1 s per non-linear step on DIII-D that is two to five minutes for a reactive
rung, which is affordable.

### 3.7 The handoff

`Globalisation::PicardThenNewton` on the bordered path is then: run the frozen
iteration to a **loose tolerance**, clear the flag, and continue the same loop
with the real Jacobian. It is one function, one state vector and no rebuild.

Two rules carry over from `CLAUDE_HDGGS.md`, *Picard, then Newton*, and both
were learned the expensive way:

* **Stage 1 failing to converge is not an error.** §4.5 converges at both orders
  from a Picard state that never met its own tolerance.
* **Do not replace the tolerance with an iteration budget.** The handoff is not
  monotone in Picard effort — budgets of 400 and 3 converge where 40 and 10
  fail. Drive it to a tolerance.

### 3.8 What family (a) fails to fix

* The three **no-axis** rows (§1.3). The constraint has degenerated; the path is
  irrelevant.
* C under filaments, which converges in 3 Newton steps with a properly located
  axis and is then refused by the topology guard for a plasma containing
  `r = 0`. A globalisation that found a *different* branch might help, but
  nothing in a Picard says it would find the right one, and M-26 says a
  different route is entitled to a different solution.
* Anything that is a bad **input**. M-103 runs the shipped configurations;
  `CLAUDE_FB.md`'s §7.15 conclusion — *choose consistent inputs*, Shafranov's
  `B_v` for the prescribed `I_p` — is a different and possibly cheaper lever.

---

## 4. Family (b): freeze the borders and hand off

### 4.1 It exists, and it is `[solver] PicardSweeps`

`runPicardPrestage()` (`apps/meq.cpp:2412-2617`) is family (b), built and
shipped. Per sweep it:

* freezes `(ψ_ax, ψ_bnd)` on the source via `setNormalisation()`, and the plasma
  edge via `freezePlasmaEdge()`;
* builds an **unbordered** solver — `makeSolver( mesh, first, /*bordered=*/false )`
  (`apps/meq.cpp:1714-1732`), which calls `setSource( Source const & )` rather
  than the `NormalisedSource` overload, and skips the limiter, the X-point and
  the plasma-current borders (`apps/meq.cpp:1744-1754`, `:1784-1789`);
* sets `Globalisation::PicardThenNewton` on it — **legal there and nowhere else
  in the run**, because with the normalisation frozen there is no `ψ_ax`
  unknown (`apps/meq.cpp:2476-2477`);
* solves at a deliberately loose `1e-4` relative tolerance, treating a stall as
  usable (`apps/meq.cpp:2495-2540`);
* under-relaxes the **field** by `PicardBlend` and re-reads `ψ_ax` at the
  located O-point, `ψ_bnd` at the bounding point, and rescales to the target
  `I_p`.

MEQ therefore already runs **three** nested iterations on a machine case: this
pre-stage, then `PlasmaSupportSweeps` (4 on every failing machine file), then
the bordered Newton with its Armijo. Family (b) is not a new pattern; it is the
outer-freeze pattern MEQ already applies to the plasma support, applied to the
normalisation as well.

### 4.2 It changes the PROBLEM, and here is exactly which one

This must be said plainly, because it is the difference between a globalisation
and a two-solver pipeline. What each sweep solves is

```
−∇̄·( (1/r) ∇̄ψ ) = F( r, z, ψ ; ψ_ax*, ψ_bnd*, λ* ) / r
```

with `ψ_ax*`, `ψ_bnd*` and the profile scale `λ*` **constants**. That is an
ordinary semi-linear Grad–Shafranov problem — `meq::NormalisedSource` reduced to
a `meq::Source` — and it has a different solution set from the bordered system
for every `(ψ_ax*, ψ_bnd*, λ*)` that is not the answer. `CLAUDE_HDGGS.md`,
*What did not work, in the order it was tried*, item 1 is the measurement:
*"A fixed `ψ_ax` is not a simplification of a normalised profile; unless the
value happens to be right it is a **different problem**"*, and
`theFixedNormalisationIsADifferentProblem` is the standing control.

So the only question that matters is §4.3, and the answer is not reassuring.

### 4.3 Is the handoff state in the bordered Newton's basin?

**Measured on case A only, and the answer there is no.** From `c12501d`'s commit
message, which is the only record:

* reading `ψ_bnd` at the **fixed seed** runs away — `I_p` goes from `3.2e+05` to
  `4.4e+06` in one sweep — because the seed sits inside the plasma and `ψ` there
  rises as the current concentrates;
* **tracking the saddle instead**, and under-relaxing the field and the scale as
  `freegs4e` does, **inverts the span at the first sweep**;
* each sweep's inner solve **stalls rather than diverging** — `8.985e-01` to
  `1.26e-01` in two steps, then a limit cycle at about `1e-03`.

The third bullet is the informative one: **the inner solve is not the obstacle.**
A 0.1% solve of the frozen problem is available and the outer map is what does
not contract. `tools/freegs4e-benchmark/mkpicardguess.py` records the same
runaway from a standalone reimplementation: the core goes from 310 cells to 1
over 150 sweeps while `ψ_ax` and `ψ_bnd` chase each other.

**And the fixed point being iterated is not the equation being solved.** The
outer map is `(ψ_ax, ψ_bnd, λ) ↦ read them off the solution of the frozen
problem`, which is a *Gauss–Seidel* split of a system the bordered Newton solves
*simultaneously*. §2.4's M-31 is the same split measured from the other side and
it does not converge either. That is not proof that family (b) cannot work —
`freegs4e` runs this map successfully on all seven — but it is a reason to
expect the blend, not the inner solver, to be where the work is.

### 4.4 A defect in it, READ FROM THE CODE AND NOT MEASURED

**The connectivity fill never runs in the pre-stage, so the pre-stage and the
bordered solve are posed on two different supports.**

`plasmaComponentWanted()` (`GradShafranov.cpp:4133-4141`) requires
`normalisedSource != nullptr`. The pre-stage solver received its source through
`setSource( Source const & )` (`apps/meq.cpp:1731`, with `source = normalised`
at `apps/meq.cpp:1021`), so that member is **null** and the predicate is false.
Therefore:

* `refreshPlasmaComponent()` returns immediately at its first line
  (`GradShafranov.cpp:4197`), so the driver's explicit
  `step->refreshPlasmaComponent( *picardState )` at `apps/meq.cpp:2521` is
  **inert**;
* `setPlasmaSupportFrozen( true )` beside it freezes nothing;
* the pointwise `Ψ > 0` test still applies, because it lives inside
  `NormalisedMHDSource::f()` (`Source.cpp:100-101`);
* so the pre-stage solves the **level-set** support problem and hands its answer
  to a solve posed on the **connected component** — which is precisely the
  distinction XP-1 exists to make, and `setPlasmaConnectivity`'s default is
  `Component` (`GradShafranov.hpp:3177`).

This is the same shape as M-82's defect — *"`setPlasmaSupport()` only reaches the
fill if it is called on the wrapper"* — and it hides the same way: every
experiment on the pre-stage has had the fill off in **both** arms.

**It is a code reading and not a measurement.** The experiment is one `printf` of
`plasmaComponentWanted()` and `componentCount()` inside the pre-stage on
`machine-a-testtokamak`, and it costs one run.

### 4.5 What family (b) fails to fix

Everything family (a) fails to fix (§3.8), plus one more: **it cannot use the
borders to help itself.** The prescribed plasma current is what
`CLAUDE_FB.md`, *The plasma current is a border unknown*, measures as turning an
ill-conditioned non-linear eigenvalue balance into an ordinary unknown — and the
pre-stage discards it (`apps/meq.cpp:1784-1789` gates `setPlasmaCurrent()` on
`bordered`), replacing it with an explicit rescale of `λ` between sweeps. That
is the amplitude-fixed problem §7.13 measured failing *from everywhere*, re-posed
once per sweep.

---

## 5. Family (c): continuation in `μ₀ I_p`

Not in the brief, cheap, and **explicitly sanctioned** by `CLAUDE_FB.md`:

> *Ramping `I_p` is the exception and is fine: the plasma is present throughout,
> and it is prescribed input rather than something recovered from a black-box
> `F`.*

The general objection to continuation — `CLAUDE_HDGGS.md`, *THIS MUST NOT GO
INTO THE DRIVER* — is that `meq::Source` exposes `f()` and `dFdPsi()` and has no
amplitude parameter, so `F_λ = λF` is the only black-box homotopy and it
degenerates to the harmonic problem at `λ = 0`. **`μ₀ I_p` is not that.** It is
`setPlasmaCurrent()`'s argument, a number the file supplies, and at any value on
the path the plasma exists and every border means what it means.

The rung would be: on observed failure, solve at `μ₀ I_p` scaled by
`θ ∈ (0, 1]`, walk `θ` to 1 with the adaptive step control that already worked
on §4.4's `c₃` ramp (halve on failure, grow by 1.3 on success — 9 solves with 2
retreats), warm-starting each solve from the last. Every machine that fails in
M-103 prescribes a current, so the lever exists on all six.

**What it does not address**: the no-axis rows again, and the possibility that
the *target* current is itself inconsistent with the conductor currents, which
is §7.15's finding and is an input problem.

**Its own falsifier**: if `θ` has to fall below about `1e-2` to take the first
step, there is a limit point on the branch and the equilibrium being asked for
does not continuously connect to a low-current one. §4.4's ramp records that the
step "never had to fall below 1e-3", which is the shape of a path with no limit
point.

---

## 6. What none of them fixes

### 6.1 The degenerate axis row

§1.3's mechanism. Three of M-103's fourteen rows converge on a system whose
`ψ_ax` equation is `ψ_ax ← max ψ_h` with an empty Jacobian row. The candidate
repairs, none of which is a globalisation:

* **refuse.** If `constraintLocated` is false at every iteration of a converged
  solve, the run did not solve the problem it was given. `checkAxis()` already
  warns (`absent`) and deliberately does not refuse, because a wall-hugging
  annulus is a real thing to look at — but *converging with a decoupled border*
  is a different statement from *the answer has no O-point*, and only the first
  is a solver fact.
* **report it.** The `.nc` carries `axis_normalised_flux` when the axis is
  located and nothing when it is not, which is why M-103's harness note says
  *"the verdict is the datum and the norm is not"*. A second attribute saying
  the axis border was decoupled for `n` of `m` iterations costs one integer.
* **constrain the topology.** Nothing in the bordered system says the plasma is
  a core (`CLAUDE_FB.md`: *"Every constraint is satisfied by it ... and none
  says the plasma is a core"*). That is a real design question and it is
  `FREE-BOUNDARY-PLAN.md` §11's, not this file's.

### 6.2 The X-point excursion

M-103: both no-axis rows put the located X-point **more than a metre** from its
seed — A at 1.096 m, G-shaped landing at `( 0.053, −0.008 )`, essentially the
symmetry axis. The X-point's two unknowns move only themselves inside the
damping loop (`GradShafranov.cpp:6977-6980`, against the zero column at
`:6836`), so a cap on their per-step excursion is a three-line change *inside*
the existing Armijo — reject a trial that moves `x_X` further than some multiple
of the local element size, exactly as `refreshXPoint()` already rejects one that
leaves the mesh.

**But it is a branch control, not a step-length rule**, and CLAUDE.md's standing
rule applies with full force: it would change which equilibrium is reported on
runs that converge today. It must be opt-in or reactive, and it needs its own
measurement. Listed here so it is not confused with a globalisation.

### 6.3 The topology guard's refusal

C under filaments converges in 3 Newton steps to a well-located axis and is
refused because `−ψ_bnd/( ψ_ax − ψ_bnd ) = +2.8555e-01` puts the plasma over
`r = 0`. That is the guard working. Whether a different branch exists on that
machine is a question for a continuation or for the inputs, and M-103 does not
answer it.

---

## 7. Recommendation, ranked

**1. Diagnose before building (BG-0).** It is not established that the exit-2
mode is a path failure at all, and if it is a non-finite direction no
globalisation in this document helps. Six runs and a diagnostic block.

**2. Repair and measure family (b), which is already built (BG-1).** Fix §4.4's
connectivity gap, then run `PicardSweeps ∈ { 2, 4, 8 }` on the six exit-2
configurations. This is the **control** any new rung has to beat, and it costs
one bug fix and fourteen runs. A "does not yet work" measured on one machine,
with the fill off, is not a verdict.

**3. Build family (a) — the bordered Picard (BG-2 to BG-4).** It is the route
that is *right* rather than the route that is cheap: it keeps every border
exact, keeps every guard live, changes no form so it cannot trip §3.5's
dangling-`DarcyForm` hazard, and it fixes the thing family (b) structurally
cannot — the normalisation and the field moving **together** rather than in
alternation. §2.4 and §4.3 are two independent measurements of what happens when
they alternate.

**4. Family (c), `μ₀ I_p` continuation, as a second rung (BG-7).** Cheap,
sanctioned, orthogonal to (a) and (b), and aimed squarely at the
prescribed-current family that fails. Deferred behind (a) only because a rung
that is tried after a rung is worth less than the first rung.

**5. Rejected outright, with reasons**: KINSOL on the bordered residual (§3.6 —
it cannot take the border); a dogleg, trust region or JFNK (§2.2 — the NPC
Jacobian is solve-only); `Normalisation::Decoupled` as a rung (§2.4 — measured
failing, M-31); any predictive trigger (§1.4 — two signals measured useless,
one of them anti-correlated, and a third that looks the same way).

---

## 8. The staged plan

Each stage ends at a **measured** acceptance criterion and names the file the
assertion lives in. Stages are BG-0 to BG-7 so they do not collide with FB-*,
XP-*, IN-* or PE-*.

### BG-0 — Classify the exit-2 failure

**Do.** Instrument the `!accepted` branch (`GradShafranov.cpp:7112-7116`) to
record, per trial: whether `y`, `z` and every queued column were finite before
the loop; the `what()` of any caught exception; and `trialNorm` otherwise.
Extend the thrown message to name which of §0.1's three causes fired. Run the
six exit-2 configurations both arms.

**Acceptance.** A new `MEASUREMENTS.md` anchor (M-104 if free) classifying all
fourteen M-103 configurations into { non-finite direction, throwing constraint,
non-improving step, converged-degenerate }, with the failing iteration number and
the cause for each. **Not "it runs": the table is the deliverable.**

**Regression.** `tests/convergence/SolverContract.cpp` — a case asserting the
three causes produce three distinguishable messages, provoked by (i) a source
that refuses the normalisation, (ii) an X-point seed outside the mesh, and (iii)
a deliberately singular border. The existing
`"the bordered Jacobian is singular in ( psi_ax, psi_bnd, a )"` throw
(`GradShafranov.cpp:6964`) is the model.

**Gate.** If the modal cause is a non-finite direction, stop here and open a
different plan.

### BG-1 — Repair `PicardSweeps` and measure it as the control

**Do.** Give the pre-stage solver a real component fill. The cleanest route is
to let the unbordered solver know its source is a `NormalisedSource` without
making `ψ_ax` an unknown — a `setSource( NormalisedSource &, double )` overload
that takes the normalisation as **fixed**, or an explicit
`setNormalisedSourceFixed()`. Then `plasmaComponentWanted()` is true, the
driver's existing `refreshPlasmaComponent()` / `setPlasmaSupportFrozen()` calls
stop being inert, and the pre-stage and the bordered solve are posed on the same
support.

**Acceptance.** Two numbers, in the same anchor as BG-0 or the next one:
`componentCount()` in the pre-stage before and after (it is 1 or more after, and
the fill is reported as having run at all); and the six exit-2 configurations at
`PicardSweeps ∈ { 2, 4, 8 }`, reporting how many close and what `ψ_ax` each
reaches against its reference.

**Regression.** `tests/convergence/PlasmaConnectivity.cpp` — the fill runs on a
solver whose normalisation is fixed, asserted by `componentCount()` and by the
element count the fill reaches, against the pointwise control.

**Gate.** If four or more of the six close, stop: family (b) is the answer and
BG-2 onward is optional work.

### BG-2 — The frozen source term, unbordered

**Do.** `meq::SourceIntegrator::setFrozenPotential( mfem::GridFunction const * )`
per §3.2: `F` evaluated at the frozen field, the reaction term returning zero,
the element mask and `fOutsidePlasma()` untouched. Reachable on the unbordered
path first, where there is a Newton to compare against.

**Acceptance.** On `tests/analytic/ManufacturedNonlinear.hpp` at `k = 2`, the
frozen iteration driven to `1e-12` and plain Newton reach the same discrete
solution to **better than 1e-10 in `max ψ_h`** — the same assertion shape as
`andersonPicardReachesTheSameSolutionAsNewton`. Separately,
`GetNumLocalNLIterations()` still reads 0, and `SolovievConvergence`'s table is
**bit-identical** (a Solov'ev source has `∂F/∂ψ ≡ 0`, so freezing must change
nothing at all — that is the sharpest available check that the frozen path is
the same discretisation).

**File.** `tests/convergence/SolverContract.cpp`, beside the existing ordering
and source contracts.

### BG-3 — The bordered Picard

**Do.** Lift the two throws for `Globalisation::PicardOnly` alone. In
`solveWithNormalisation()`, hold a local `bool frozen` that drives
`SourceIntegrator::setFrozenPotential()` and is updated to the accepted iterate's
potential after each step; hand `assembleNormalisationColumn()` the same frozen
state (§3.3). Add relaxation on the field increment, default 1.0 per
`setPicardDamping()`'s existing default and for the same measured reason.
**Do not** call `setGlobalisation()` inside the loop (§3.5). **Do not** wire
KINSOL (§3.6). `usesNonlinearForms()` gains the `normalisedSource` clause and
`SolverContract.cpp`'s bit-exactness cases are re-run because of it.

**Acceptance.** On `tests/analytic/HighBetaPoloidal.hpp` at `ν = 2, A = 1`,
`k = 2`, `n = 8` — the case M-30 publishes — the bordered Picard reaches the
same `ψ_ax` as the bordered Newton to **at least eight significant figures**
(M-30's published value is `3.059006e-01`), and `HighBetaConvergence`'s own
table is unmoved. Report the iteration count; it is expected to be one to two
orders larger than Newton's 4.

**File.** `tests/convergence/HighBetaConvergence.cpp`, as
`theBorderedPicardReachesTheSameSolutionAsTheBorderedNewton`.

### BG-4 — The handoff, `PicardThenNewton` on the bordered path

**Do.** Lift the throws for `PicardThenNewton`; run BG-3's iteration to a loose
relative tolerance (start at `1e-4`, which is what the driver's pre-stage
already uses and for the same stated reason), clear the frozen flag, continue.
Stage 1 failing to meet its tolerance is not an error. Carry
`picardIterationCount` as the unbordered path does.

**Acceptance.** Two assertions, and the second is the one with teeth:

1. On a configuration where `Globalisation::None` **throws** — the first of
   BG-0's exit-2 cases reduced to a test-sized mesh — the rung converges, and
   the converged state satisfies every border: `|ψ_ax − ψ_h( x* )|`,
   `|ψ_bnd − ψ_h( x_limiter )|` and `|∫F/r − μ₀I_p|/μ₀I_p` all below `1e-10`,
   with a **located** axis (`constraintLocated` true at the last iteration).
2. On a configuration where `None` **converges** — the diverted fixture in
   `tests/convergence/DivertedMachine.hpp` — the rung reaches the **same**
   equilibrium: `ψ_ax`, `ψ_bnd` and the X-point position agreeing with the
   `None` answer to `1e-9`. M-26's hazard is that it will not, and this is the
   assertion that finds out.

**File.** A new `tests/convergence/BorderedGlobalisation.cpp`, registered as a
ctest. **Not** `FreeBoundaryCoupling.cpp`, which is already the longest case in
the suite (M-14's discussion has it at 281–433 s across three readings of one
tree).

### BG-5 — The driver rung, and saying which rung answered

**Do.** Replace the immediate `return SolveFailed` at `apps/meq.cpp:2752-2763`
with one retry on the bordered `PicardThenNewton`, rebuilding the solver first
(a caught `ErrorException` leaves one unusable — `CLAUDE_HDGGS.md`, *On
SUNDIALS*), and keep the existing message as the text printed when the retry
also fails. Record the rung in the `.nc` as a global attribute beside
`newton_iterations`, e.g. `globalisation = "none"` or `"picard-then-newton"`,
and in the `Cycle` record the driver already keeps (`globalised`).

**Acceptance.** `DriverAcceptance.cpp`: a configuration that fails cold and
succeeds on the rung, with the `.nc` attribute naming the rung; **and every
shipped example unchanged** — the driver-against-library pins
(1.189e-16 over 15,360 dofs, 1.6e-16 on the curved path, 4.4e-14 on the adaptive
loop) all hold to the digit, because nothing on a converging path may move.

### BG-6 — Re-take M-103

**Do.** Re-run all seven machines, both conductor models, shipped configurations
plus the new rung, `OMP = MKL = 16`, both arms — which is the only thing that
can tell a conversion defect from a pre-existing one, per M-96's own note.

**Acceptance.** A new anchor with M-103's columns plus a rung column and the
`compare.py` verdict, refusing `absent` as well as `bad` per M-103's harness
note. **The number that decides it is how many of the six exit-2 configurations
close with a located axis.**

### BG-7 — Continuation in `μ₀ I_p`, conditional on BG-6

**Do.** Only if BG-6 leaves exit-2 rows standing. A `θ`-ramp on
`setPlasmaCurrent()`'s argument with §4.4's adaptive step control, warm-started
from the previous `θ`, as a third rung after `PicardThenNewton`.

**Acceptance.** The number of solves, the number of retreats, and the smallest
`θ`-step taken, on each configuration it is run on — the same three numbers
§4.4's `c₃` ramp reports. If the step has to fall below `1e-2`, the branch has a
limit point and the ramp is the wrong tool; say so rather than tuning it.

---

## 9. Keeping it reactive, and how a reader tells which branch was reached

Three commitments, all of them enforceable by a test:

* **No new default.** `Globalisation::None` stays the default in the class and
  the rung is entered only from a `catch`. BG-5's acceptance includes every
  shipped example being bit-unchanged, which is what enforces it.
* **No new TOML key that selects a globalisation.** `Globalisation` is not
  exposed today and must not become exposed — M-26 is the reason and it is
  unchanged by anything here. `PicardSweeps` stays what it is: an opt-in
  initialiser, off by default, documented as changing which equilibrium is
  reported.
* **The answer says how it was reached.** The `.nc` gains one attribute; the
  driver already prints the residual history on failure
  (`reportResiduals`) and already tracks `globalised` per cycle. Without this a
  reader cannot tell M-26's 9.4% from a converged equilibrium, and M-103's own
  harness note — *"the verdict is the datum and the norm is not"* — is the
  precedent for making the verdict machine-readable.

---

## 10. Falsification

Any one of these kills the route, and each is cheaper than the stage it gates.

1. **BG-0 finds the modal exit-2 cause is a non-finite direction.** Then the
   defect is in the bordered elimination or the trace solve and no step rule or
   operator swap reaches it. The plan is void and the investigation is a
   different one.
2. **BG-1's repaired `PicardSweeps` closes four or more of the six.** Then
   family (b) works, family (a) is redundant, and the correct outcome is a
   documentation change plus the connectivity fix.
3. **BG-3's bordered Picard does not contract on any failing configuration
   within 300 iterations**, or limit-cycles the way the driver's pre-stage does
   at `1e-3`. Then the bordered fixed point is no better conditioned than the
   alternating one and family (a)'s central claim — that simultaneous beats
   alternating — is false.
4. **BG-4's second assertion fails**: the rung converges on a case `None`
   already solves and reports a *different* equilibrium. Then the rung is a
   branch changer on working problems and may never be entered automatically,
   which removes most of its value even where it converges.
5. **BG-6 closes fewer than two of the six.** A rung that rescues one
   configuration out of fourteen is not worth a permanent second code path;
   record the measurement and delete the rung.
6. **Cost.** At M-90's 1.1 s per non-linear step on DIII-D, a rung needing
   M-36-sized Picard counts (122–290) costs two to five minutes per failed
   machine. That is affordable reactively. **Ten times that is not**, and if the
   bordered fixed point needs thousands of iterations the route dies on the wall
   clock regardless of whether it converges.

---

## 11. What this document does NOT resolve

Eleven things, each with the measurement that would settle it. **No number
below appears anywhere in this tree and none is invented here.**

1. **Whether the exit-2 failures are path failures.** BG-0. Everything else is
   conditional on it.
2. **Whether `DarcyNPCOperator::Mult` is defined when the potential block is on
   the linear form.** Read `darcyhybridization.cpp`'s `NPCResidual` and its
   preconditions. It decides whether §3.2's rejected route was even available.
3. **Whether the constant Picard Jacobian can be cashed as one factorisation.**
   Read `CanCacheCondensation()` and `SetReuseSymbolic()`, then measure
   `GetComputeHTime()` over a frozen run against a Newton run of the same length
   with `tests/performance/NewtonStepProfile`. It is the difference between a
   rung that costs the same per step and one that costs a fraction.
4. **Whether Anderson is needed on the bordered fixed point, and where it would
   come from.** KINSOL cannot take it (§3.6). Measure the undamped and damped
   iteration counts in BG-3 first; only if they are the binding constraint does
   the Eigen-versus-hand-rolled question arise.
5. **What `γ` should be for a Picard phase.** It is frozen at `‖c‖` from the
   first iterate so the printed history compares like with like
   (`GradShafranov.cpp:6150-6182`). A phase whose first iterate is a different
   state gets a different `γ`, and the handoff would then print two histories in
   two scales. Decide it before BG-4 prints anything an order is read off.
6. **Whether to freeze at the last accepted iterate or the last trial.** The
   Armijo loop evaluates up to twelve trials per step; freezing at the accepted
   one is the defensible choice and costs one extra copy. Not measured.
7. **Whether §4.4's connectivity gap is actually the reason `PicardSweeps` fails
   on case A.** It is a code reading. One `printf` of
   `plasmaComponentWanted()` and `componentCount()` settles it.
8. **Whether the pre-stage belongs inside or outside the support loop.** Today
   it runs once, before it (`apps/meq.cpp:2724-2733`), while
   `PlasmaSupportSweeps = 4` on every failing machine — so there are two nested
   outer freezes and nobody has measured the ordering. A 2×2 cross over
   `{ pre-stage first, pre-stage per support sweep } × { fill on, fill off }` is
   the experiment, and M-82's lesson is that a one-key experiment separates two
   hypotheses only if everything else is where you think it is.
9. **Whether the degenerate axis row should refuse.** §6.1. It is a
   branch-selection decision, it is `FREE-BOUNDARY-PLAN.md` §11's, and it is the
   only thing that addresses three of M-103's fourteen rows.
10. **Whether the X-point wants an excursion cap.** §6.2. Three lines, and it
    changes which equilibrium converging runs report, so it needs its own
    measurement on `DivertedMachine.hpp` before it is entertained.
11. **Whether the shipped machine inputs are themselves consistent.** M-103 runs
    them untuned and `CLAUDE_FB.md` §7.15's conclusion is that consistent inputs
    are the cheap fix — Shafranov's `B_v` for the prescribed `I_p`. `PsiAxis`
    seeding is known not to help on case A and has not been tried on the six.
    If the inputs are the problem, a globalisation is being asked to find an
    equilibrium that the conductors and the current target do not jointly
    describe.
