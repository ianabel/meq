# The HDG Grad-Shafranov solver: the equation, the discretisation, the solve

**`CLAUDE.md`'s core-solver half** — stages 0 to 6 of the port, and everything
about how MEQ solves the equation rather than about what it is pointed at.
`CLAUDE.md` is still the maintainer's index and carries the build, the commands,
the traps and the layout. The two campaign files are `CLAUDE_FB.md` (free
boundary) and `CLAUDE_INVERSION.md` (flux-surface inversion); both inherit
everything below.

**`refs/Refs.md` indexes the papers.** The two marked ✔ there are not background
reading — they *are* the method, and `src/meq` is an implementation of them.

**Measurement tables live in `MEASUREMENTS.md`** under stable `M-nn` anchors.
**Do not renumber the anchors.**

**Toroidal flow is `CLAUDE_FLOW.md`**, being a campaign rather than part of the
core method — though a rotating equilibrium is the same equation with a different
`F`, so everything below is what it is driven through.

## The equation being solved

Fixed boundary: the plasma boundary `Γ` is known and taken to be the level set
`ψ = 0`, so this is an interior Dirichlet problem. From
`refs/HDG-GradShafranov.pdf` eqs (1)–(4):

```
-∇̄·( (1/R) ∇̄ψ ) = F(R,z,ψ) / R     in Ω ⊂ R²
ψ = 0                               on Γ = ∂Ω

F(R,z,ψ) := μ₀ R² dp/dψ + g dg/dψ
```

`∇̄ := (∂_r, ∂_z)` acts formally like a vector of partial derivatives independent
of the coordinate system — it is *not* the cylindrical gradient, and the
distinction is the whole content of the `1/R` and `R` weights below. `p(ψ)` is
the plasma pressure and `g(ψ)/R` the toroidal field function; both are user
input, and it is their `ψ`-dependence that makes the problem semi-linear.

Recast as a first-order system by introducing the **flux** `q = (1/R)∇̄ψ`:

```
q − (1/R)∇̄ψ = 0,      −∇̄·q = F/R,      ψ = 0 on Γ
```

Introducing `q` is not a numerical convenience. The physically interesting output
is the magnetic field, which is built from `∇ψ`, so discretising `q` directly is
the reason to prefer HDG here — it gives the derivative at the same order as the
potential rather than one order down.

### The two papers disagree about the sign of the Solov'ev source

Checked, resolved, and worth not rediscovering.

`refs/HDG-GradShafranov.pdf` eq. (10) gives `F = −((1−A)r² + A)`.
`refs/HDG-GradShafranov-Adaptive.pdf` eq. (21) gives `F = +((1−A)r² + A)`.

**The first is right**, and the second contradicts its own eq. (1). Applying
`Δ*` to the particular solution *both* papers publish settles it analytically:

```
Δ*( r⁴/8 )            =  r²
Δ*( (A/2) R² ln R )   =  A
Δ*( −A r⁴/8 )         = −A r²
                        ─────────────────
Δ*( ψ_P )             =  (1−A) r² + A
```

and since both papers define `−Δ*ψ = F`, `F = −((1−A)r² + A)`. The twelve
homogeneous terms contribute nothing, being `Δ*`-harmonic. This is confirmed
numerically: `Soloviev.hpp`'s `deltaStarFD()` recomputes `Δ*ψ` by central
differences, and `SolovievConvergence.cpp`'s `solovievSourceMatchesTheOperator`
asserts it against `−f()` over the benchmark rectangle, so the whole twelve-term
transcription is checked rather than trusted, and checked against the same `f()`
the solver is actually fed.

Take this as the standing warning about the benchmark: **the published
coefficients are not self-checking**, and a sign error here would show up as a
solver that converges beautifully to the wrong equilibrium.

### A third erratum: eq. (20)'s `η₅` does not vanish on the exact solution

`ψ̂_h ∈ P_k(e)`, but `ψ*`'s trace has degree `k+1`, which no element of `M_h` can
represent. Orthogonality splits the term:

```
‖ψ̂ − ψ*‖²_e = ‖ψ̂ − P_M ψ*‖²_e + ‖(I − P_M) ψ*‖²_e
```

The second piece survives even when you substitute the *exact* solution and the
best possible trace, at `O(h^{k+3/2})`, which `h_e^{-1}` and `O(h^{-2})` edges
turn into an `O(h^k)` floor. Measured, the printed term **is** that floor —
agreeing with it to between 0.05% and 1.8% at every `k` and every mesh — so it
converges at `k`, not `k+1`, and drags the total to 1.44/2.32/3.55 instead of
1.99/2.99/3.97.

Taking the difference *inside* `M_h` restores `k+1` and reproduces the rates
GS-2's own Table 1 reports. `TraceComparison::Projected` is the default;
`Literal` is kept so the suite keeps measuring the difference.

**A separate `η₅` problem on the extension path**, and this one was nearly fatal
to the adaptive loop: on `Γ_h` the term compares `ψ*` against a trace pinned to
*zero* rather than the `φ_h` actually imposed, so the difference is
`O(dist(Γ_h, Γ)) = O(h)`. Unmitigated, `η = 4.09e-1` where `η₁ = 2.12e-3`,
converging at ~0.5 — the loop would have run, produced plausible pictures, and
refined the wrong elements.

**THAT WAS AN OMISSION AND IT IS A REPAIR.**
**Excluding those faces is the wrong repair** — it restores `k+1` by deleting
the term exactly where the geometry error lives. `setTransferredBoundary()` takes
the datum as a second argument instead — `GradShafranovSolver::transferredDatum()`, built on
`mfem::TransferredDatumCoefficient` — so `η₅` compares `ψ*` against the `φ_h`
actually imposed and the faces stay in. Measured on the extension benchmark at
`k = 2` over `h = 0.213 / 0.106 / 0.053`, all three treatments on one solve:

→ **[M-15](MEASUREMENTS.md#m-15)** — coarsest · finest · rate

So the term is now **two percent of `η₁`** rather than four orders larger than
it, and it converges at `k+1` rather than at a half.
`theTransferredDatumRestoresEtaFive` in `ExtensionConvergence.cpp` is that table,
and it keeps the pinned column as its control: if that ever converges, the
comparison is empty.

**η BARELY MOVES, AND THAT IS THE POINT RATHER THAN A DISAPPOINTMENT.** The
excluded and datum columns agree to four figures, because a correctly evaluated
`η₅` on `Γ_h` is small. What changes is that the boundary elements now *have* an
`η₅` contribution instead of none, so the marking sees them. The adaptive loop is
unchanged — 97 → 254 → 342 → 449, `η` 4.7352e-04 → 6.8668e-05 — which is what
says the repair did not perturb what already worked.

**THE SIGN IS THE THING TO GET RIGHT, AND IT IS TESTED SEPARATELY.**
`transferredDatum()` must hand `mfem::PathLiftCoefficient` the **raw flux
block**, which holds `−q`, and not `flux()`, whose sign is undone. Feed it the
wrong one and `φ_h` comes back as `−ψ` rather than `ψ`, because the `g(a(x))`
that would otherwise survive is zero on `Γ`. So the failure is the answer with
its sign reversed, not a small bias.
`theTransferredDatumReproducesTheImposedCondition` checks `φ_h` against the
exact `ψ` it transfers — 1.65e-3 → 4.0e-5 relative, converging at 3.25 then 3.63
— and a sign error reads about 2.

## The discretisation

HDG, in the LDG-H form of `refs/HDG-GradShafranov.pdf` eq. (8) — restated in
`refs/HDG-GradShafranov-Adaptive.pdf` eq. (13) with the block structure made
explicit, which is the more useful version when writing assembly.

```
(R q_h, v)_Th + (ψ_h, ∇̄·v)_Th − ⟨ψ̂_h, v·n⟩_∂Th   = 0
(q_h, ∇̄w)_Th − ⟨q̂_h·n, w⟩_∂Th                    = (F/R, w)_Th
⟨q̂_h·n, μ⟩_∂Th\Γh                                 = 0
ψ̂_h = φ_h on Γh

q̂_h·n := q_h·n + τ(ψ_h − ψ̂_h)        on ∂Th
```

Spaces, all of the **same degree** `k` — hybridization removes the inf-sup
compatibility condition that a classical mixed method would impose:

```
V_h = [P_k(K)]²    flux            L2_FECollection, vdim 2
W_h = P_k(K)       potential       L2_FECollection
M_h = P_k(e)       hybrid trace    DG_Interface_FECollection
```

**The volume spaces are on the closed Gauss-Lobatto basis, and that is a
convention rather than a requirement** — `convdiff.cpp`'s "as it is customary for
HDG to match trace DOFs", quoted and inherited until it was measured. A nodal
basis does not change the *space*: Lobatto and Legendre both span `P_k(K)`, so
the discretisation cannot see the choice, and nothing in the hybridization needs
volume dofs to sit on faces — every coupling is a face **integral**, computed by
quadrature against both bases. Measured, with both volume spaces switched to
`GaussLegendre` and nothing else touched: rates 1.996/3.001/3.998 in `ψ`, the
curved benchmark's `L2` at 2.742813e-05 down both paths agreeing to 6.6e-15, and
the `k = 3` Newton history reproducing the one recorded below to seven digits.
The finest Solov'ev error reads 5.510484e-10 against 5.510486e-10 — the last
place.

So it is kept for alignment with the miniapp MEQ was ported from, and for nothing
else. **What it costs** is that a dof is a point value *on* the element boundary,
where an L2 field is discontinuous, so reading `W_h` by nodal interpolation at
another mesh's dof points is ambiguous — measured, 9% to 28% wrong.
`meq::FieldTransfer` projects rather than interpolates, which is basis-agnostic
and is what a non-nested transfer needs anyway.

### The assembled flux is −q, and the papers' τ carries the wrong sign

Two conventions settled by measurement rather than argument, both of which cost
time if rediscovered. `tests/convergence/SolovievConvergence.cpp` records what
every alternative produces.

**`DarcyForm` holds `−q`, not `q`.** It is built for `u = −k ∇p`, the opposite
sign to `q = (1/R)∇̄ψ`, and the two integrators that make the hybridization
consistent — `NormalTraceJumpIntegrator` and the trace rows of
`HDGDiffusionIntegrator` — have that sign baked in and take no scaling argument,
so there is no way to flip it in the assembly. `GradShafranovSolver::flux()`
undoes it once, into a separate GridFunction; the block vector stays in
`DarcyForm`'s convention, because a stage-4 Newton residual assembled by
`DarcyForm` will expect it there. Do not change one without the other.

The same convention makes the potential right-hand side `−(F/R, w)`: with the
default `bsymmetrize = true` the second block row is assembled as
`−B q − M_p ψ = b_p`. `bsymmetrize = false` is not an escape — it does not work
here at all, giving an error flat at 3e-1 for every combination of signs.

**And `τ` carries the opposite sign to `refs/HDG-GradShafranov.pdf` eq (8e)**,
which prints `q̂·n := q·n + τ(ψ − ψ̂)`. With that paper's own `q = (1/R)∇̄ψ` and
`−∇̄·q = F/R`, testing (8a)–(8d) against `(v,w,μ) = (q_h, ψ_h, ψ̂_h)` gives

```
(R q, q) − τ‖ψ − ψ̂‖²_∂Th = 0
```

which is indefinite, and the local solves are not guaranteed invertible. The
stable sign is `−τ`, which is exactly what assembling in `DarcyForm`'s convention
with a positive `τ` produces. This is a second sign slip in the same pair of
papers as the Solov'ev source one above. It is a well-posedness argument, not an
observed failure: **both** signs converge at `k+1` on this benchmark, with
`τ > 0` in `DarcyForm`'s convention giving slightly the smaller error at `k = 1`.

**`τ = 1`.** Both papers set it there and note that optimal order needs only
`τ = O(1)`. The pre-port code used `τ = 5.0` with no recorded reason; do not
reinstate that without a measurement to justify it.

**And MFEM will not give you a constant `τ` unless you ask.**
`HDGDiffusionIntegrator`'s built-in stabilisation is

```
τ± = ( β ± (α/2)(u·n)/|u·n| ) { h⁻¹ Q }
```

— the LDG choice, scaled by the inverse local mesh size *and* by the diffusion
coefficient, and not what the papers use. Measured on the Solov'ev benchmark it
costs a full order in the flux: `q` converges at `k`, not `k+1`, while `ψ` still
converges at `k+1`, so **a study of `ψ` alone would have passed it**.

The designed way out is the `HDGStabilization` hook in
`fem/darcy/bilininteg_hdg.hpp`: subclass it, return the constant from `Eval()`,
install it with `SetStabilization()`. `meq::ConstantStabilization` is that
subclass. Read `StabValue()` there — with a hook installed it divides the
quadrature weight out before calling and multiplies it back after, so returning
a bare `τ` yields exactly `⟨τψ, w⟩`. The MFEM tree has no ready-made constant
implementation; the only subclasses are the test fixtures in
`tests/unit/fem/test_bilininteg_hdg.cpp` and `test_darcy_degenerate.cpp`, which
are the idiom to copy.

**Keeping `IsConstant()` true also keeps MEQ out of an `EvalGrad` trap.** That
header warns that omitting `EvalGrad` for a *non-constant* stabilisation gives
"no wrong answer, only slow Newton convergence — a failure that survives a
passing regression suite". A constant `τ` never calls it. One more reason to
keep `τ` constant unless something measured says otherwise, which is also the
sibling project's standing advice: **do not derive `τ` from the local
coefficient.**

### Post-processing is back, and it was free

Stage 3 dropped `ψ*` on the grounds that `q` already gives the physical quantity
at `k+1`, and noted that stage 6's estimator would need it again. It does — eq.
(20) uses `ψ*` in **four of its five terms** (`η₁`, `η₂`, `η₄`, `η₅`).

**`DarcyForm::Reconstruct()` supplies it, measured at `k+2`**: rates 3.03, 4.03,
5.00 for `k = 1, 2, 3`, and it survives the extension path and Newton. No
hand-written local solve was needed and the old `GSSolver::Postprocess()` stays
deleted. `thePostProcessedPotentialSurvivesNewton` is that measurement: MFEM
lifts the non-linear potential integrators as a Jacobian frozen at the computed
potential, which on Example 5 delivers 3.05, 4.05, 5.03 and is 47×, 113×, 125×
smaller than `ψ_h` on the finest mesh.

The measurement that justifies *requiring* `ψ*` rather than quoting it: building
`η₂` on raw `ψ_h` loses **exactly one order at every `k`** — 2.002 vs 0.998,
3.000 vs 2.002, 3.992 vs 2.981 — and is 124× to 407× larger on the finest mesh.
That is the defect the pre-modernisation estimator had.

**IT WAS BROKEN, SILENTLY AND PER ELEMENT, AND THE SHAPE OF THAT IS THE PART
WORTH KEEPING.** The local post-processing problem is a pure **Neumann** problem
by construction — the total flux driving it is in H(div), so the element's flux
balance is already satisfied and the potential is determined only up to a
constant. The element average is what closes it, and that is *unconditional*.
MFEM had been treating it as one of three choices: a flag set on the mere
presence of a non-linear integrator skipped the mean-value branch, and a
singular matrix was then factored anyway with `Factor( m, TOL = 0.0 )` testing
`|pivot| <= 0.0` and the return value discarded. Where `∂F/∂ψ` vanishes the
integrator's Jacobian *is* the zero matrix, so the term **is** singular.

**The flag is per element, so the corruption was too**, and that is what made it
invisible. With `∂F/∂ψ` vanishing on an eighth of the domain — a tabulated
profile with a flat segment, which is ordinary — elements were **20× wrong**
while the whole-domain norm read 1.87; at half the domain, 61× against 7.57.
**Any global check misses it.** And a floor of **1e-12** on `∂F/∂ψ` — twelve
orders below everything else in the problem, incapable of moving a solution —
took the ratio from 7.565317 to 0.998525. That is a singular matrix and nothing
else.

The fix is `gf-hdg-dev`'s *"The postprocessing closes on the element average,
always"*. Measured from MEQ's side on a case that reads 20.3, 64.1 and 61.6 on
dead elements without it: **1.0069 with the fix, against 1.0048 where `∂F/∂ψ`
does not vanish at all.** The driver's adaptive loop is back on the published estimator — `η`
reaching 6.87e-5 on **449** elements where the degraded one needed 1069 to reach
4.77e-4.

**MEQ carried no runtime check for this, deliberately.** A solver should not
stand permanently on guard against its dependency. The state of the defect lived
in the suite instead, as a test asserting the *correct* behaviour and failing.
`thePostProcessedPotentialIsCorrectWhereTheJacobianVanishes` still runs, now as
a regression rather than a tripwire.

**And the defect did not reach the solve — measured, not argued.** The forward
solve's local problem takes its potential block from the HDG stabilisation on
`Mnl_p`'s **interior faces**, a fixed bilinear form whatever the source does,
while `Reconstruct()` builds a different local problem whose regularisation is
what failed. `theReconstructionDefectDoesNotReachTheSolve` compares both paths
in `ψ_h` **and** `q_h` over `k = 1…4` and three meshes: worst relative
difference **1.6e-13** and **2.6e-13**, against discretisation errors from 2.2e-3
down to 7.3e-11. It also checks that post-processing leaves both bit-identical,
since the driver writes them afterwards.

**Two things to keep.** The fix is on `gf-hdg-dev` and on no other branch, so a
`meq-integration` rebuilt without it silently loses this — see the merge recipe.
And the failure mode is the transferable part: **a per-element defect that a
global norm cannot see, in a quantity four of the estimator's five terms are
built on.**

## The NPC port, and the parity gap that came back the other way

**MFEM deleted `NLOrdering::LineariseThenCondense` **, as *"a
condensation in disguise"*: it was an `Operator` on the trace alone that kept the
linearisation as hidden state, and NPC's fields are Newton state, which a
trace-only operator has nowhere to put. `SetNonlinearOrdering()` went with it.
That was a breaking change for MEQ, whose default ordering it was.

**MEQ's default is now `NonlinearOrdering::NPC`** — `mfem::DarcyNPCOperator`
with `mfem::DarcyNPCSolver`, Newton on the full `(q, ψ, ψ̂)` system with the
Jacobian solved by hybridized elimination, Nguyen, Peraire & Cockburn eqs
(14)–(18). `CondenseThenLinearise` is kept as the **backup**, and the reference
is `../mfem-hdg-dev/doc/HDG-ORDERING-API.md`.

**The port cost MEQ almost nothing structurally, and the reason is an accident
worth knowing.** `solution` has always been a three-block
`{ flux, potential, trace }` `BlockVector` on a four-entry `blockOffsets`, with
`darcyFlux`, `potentialGf` and `traceGf` all `MakeRef`'d into it. **That is the
NPC unknown, block for block.** So the fields are already in place when the
iteration returns, and `RecoverFEMSolution()` leaves the Newton path rather than
needing rework. `meq::Relinearised` — the wrapper that existed only to pair the
residual with the gradient under the deleted mode — was deleted with it.

### What NPC bought, measured

**`GetNumLocalNLIterations()` is exactly zero**, which is the acceptance signal
that this is NPC and not a condensation wearing the name.
`SolverContract::theOrderingsAgreeAndOnlyOneIteratesLocally` asserts it on
Example 5 and prints the comparison:

→ **[M-16](MEASUREMENTS.md#m-16)** — `k` · NPC local NL its · condense local NL its · L2 vs exact · relative in `ψ`

Same L2 to seven figures, same 4 Newton iterations, ψ agreeing to round-off.
That is what says the backup is a backup and not a different problem.

**AND AN NPC ITERATION IS THREE TO FOUR TIMES CHEAPER THAN A CONDENSATION
ONE**, which is what the iteration-count tables everywhere else in this file do
*not* show and is where the suite's 872 s → 499 s actually came from. Measured
Measured on two cases where both orderings converge and agree:

→ **[M-17](MEASUREMENTS.md#m-17)** — NPC · CondenseThenLinearise

**4.3× and 3.7× on wall clock at the same or fewer iterations.** So reading the
parity table as "NPC needs 17 where the condensation needs 11" understates NPC
by most of a factor of four; the honest currency here is seconds, and the local
non-linear iteration counts in the third column are what NPC deletes.

**The bordered Newton got simpler and stronger, and this is the real prize.**
With `ψ` an unknown rather than a function of the trace, two of the three
bordered quantities stop being finite differences:

→ **[M-18](MEASUREMENTS.md#m-18)** — condensation · **NPC**

and the whole frozen-seed apparatus goes away — `formSystem()` is no longer
re-called once per accepted step, because there is no element-local non-linear
solve whose initial guess could go stale. `HighBetaConvergence` reproduces
**every `ψ_ax` in the table under *The measurement* to every digit printed**, and
the constraint is now satisfied to machine zero:

**These are `AxisConstraint::NodalMaximum` numbers**, which is what the tree ran
— the point of the table is what the NPC port did to the border,
and re-measuring it under the located-axis constraint would change the variable
being studied. `ψ_ax − max ψ_h` is that constraint's own residual.

→ **[M-19](MEASUREMENTS.md#m-19)** — `ν` · `A` · `ψ_ax` · `ψ_ax − max ψ_h`, was · now · Newton, was · now

`HighBetaConvergence` went **17 s → 3.3 s**. The printed residual is a different
quantity — `‖(R, γG)‖` over the full system rather than the trace — so the
histories are not comparable across the port and `γ` is a different number; the
converged `ψ_ax` is what is comparable, and it is unchanged.

### What NPC cost, and it is the parity gap in the other direction

**THIS IS THE FINDING TO READ BEFORE ASSUMING THE PORT WAS FREE.** The
expectation going in was that NPC would remove the parity gap the deleted mode
had. It does not. It replaces it with a mirror image, and a wider one. Measured
Raw undamped Newton, cap 60:

**WIDENED FROM 7 CASES TO 14, AND ONE SENTENCE OF IT WAS
WRONG.** The seven-case version said NPC was *never* better on any case
measured. It is better on two, and they are the same two corner: well resolved,
higher order.

→ **[M-20](MEASUREMENTS.md#m-20)** — case · **NPC** · CondenseThenLinearise

**The pattern is resolution, not stiffness.** NPC loses where the mesh is too
coarse for the source and wins where it is not — `k = 2, n = 24` and
`k = 3, n = 16` go to NPC. And since an NPC step carries **no element-local
non-linear solve at all**, an equal iteration count is already an NPC win on
wall clock.

**The §4.3 rows are new and they matter, because they are the control.** The
barrier defeats *both* orderings at `k = 1` and `k = 2` on `n = 16`, so this is
not a table about one ordering being sound and the other not. (These rows are
run **without** `setInitialGuess()`, unlike the §4.3 rows under *Picard, then
Newton*, which seed the ramp — that is why `k = 2` reads a failure here and
`ok, 10` there. The guess is doing the work in that table, and it is worth not
confusing the two.)

**AND WHEN NPC FAILS IT FAILS CHEAPLY, WHICH IS NOT A SMALL THING.** Same
barrier case, `k = 1, n = 16`, both reaching the cap of 60:

→ **[M-21](MEASUREMENTS.md#m-21)** — wall clock · element-local non-linear iterations

A factor of 33 in the cost of finding out that neither works. That is exactly
the uniformity §6 promises, and it holds on an iteration going nowhere as much
as on one that is working — which makes NPC the better thing to *try first*
even on the cases it loses.

That NPC is not a free win is not a surprise to upstream, whose own §6
says *"Reach for [the reduced trace operator] unless you have a reason not to"*
and *"**NPC is not automatically faster.** Its advantage is uniformity of the
local work … not fewer floating-point operations."*

### MFEM's own test suite reproduces the parity gap, and two of its tests are red

**`MFEM_ENABLE_TESTING` IS NOW ON, and turning it on found this in the first
run.** MEQ had been building with `-DMFEM_ENABLE_TESTING=OFF` because of a CMake
bug in serial CUDA builds; that bug is fixed (see `CLAUDE.md`, *Which MFEM*), so MFEM's own
131 registered tests are available. Run over the Darcy and HDG tags on
`meq-integration` at `974871c456`:

```
MKL_NUM_THREADS=1 \
  ./tests/unit/unit_tests "[DarcyForm],[NPC],[HDG],[DarcyHybridization],[NonlinearDarcy]"
```

**72 of 74 test cases pass; 98,434 of 98,436 assertions.** Both failures are in
`tests/unit/fem/test_darcy_npc.cpp`, and both are NPC on a stiff pedestal:

| test | what fails |
|---|---|
| *"A stiff source converges by condensation and by NPC alike"* | the condensation section **passes**, the NPC section **fails** |
| *"NPC solves stiff problems LineariseThenCondense cannot"* | fails outright |

Both at `n = 32`, order 1, `σ² = 0.003`, stalling at **1.1590e-03** after 41
iterations with `local_nl_iters = 0` — so it really is NPC, and the residual is
sitting at 1e-3 exactly as MEQ's §4.5 does at its own `n = 32`. **And the call is
`RunNPC( P, 40, true, GM::Assembled )`, where the `true` is the backtracking line
search** — upstream's recommended globalisation, on, and still failing.

**This is MEQ's parity gap reproduced inside MFEM's own suite**, on their fixture
rather than MEQ's, in two tests whose *names* encode the claim it contradicts. It
is independent evidence for everything in the two sections above, and it was
free: it took one flag and one run.

**Do not read it as a regression report.** MEQ has never run these tests before —
testing was off — so there is no earlier green to compare
against, and `974871c456` predates upstream's own response commit. What it
establishes is that the finding is not about MEQ's discretisation, since it
reproduces on a fixture MEQ did not write.

### Why it fails, measured — and why a line search does not fix it

**IT IS A PROPERTY OF THE TWO ALGORITHMS, NOT A DEFECT IN EITHER, AND THE
MECHANISM IS MEASURED RATHER THAN ARGUED.**

Both orderings take a bad first step — §4.2 at `k = 1, n = 16` goes 1.14e-01 →
6.25e+00 under NPC and 7.78e+00 → 5.57e+01 under the condensation. The
condensation recovers and reaches its quadratic endgame in 14; NPC wanders
non-monotonically between 1e-2 and 1e+1 for sixty steps and never enters a
basin. Not a stall, not divergence to NaN — undamped Newton outside its basin.

**Where the Newton remainder goes, per block.** At the cold iterate
`x = [0, 0, 2.771]` — flux and potential zero, trace carrying the Dirichlet
datum — one full step gives

→ **[M-22](MEASUREMENTS.md#m-22)** — flux · potential · trace

**The flux and trace rows of the NPC residual are LINEAR in `(q, ψ, ψ̂)`** — the
flux equation and the flux constraint carry no `F` — so a Newton step annihilates
them *exactly*, to 1e-16 and 1e-13. **All the non-linearity is in the potential
row**, and so is the entire Newton remainder: 8.1e-02 → 6.25e+00. And the
correction is `O(10²)` against a solution where `ψ ~ 0.3`, because at `ψ ≡ 0` the
linearised operator `−∇̄·((1/R)∇̄·) − (∂F/∂ψ)/R` is at its most indefinite — the
pedestal sits at `max|∂F/∂ψ|/λ₁ ≈ 7`.

**What the condensation does instead is nonlinear elimination, and that is a
preconditioner.** Its element-local solves take `(q, ψ)` to the *exact* solution
of each element's non-linear problem at the current trace, so the fields never
take an unphysical linear excursion and the outer residual is a tamer function of
a much smaller unknown. That is a real and well-understood robustness mechanism,
and it is precisely what NPC gives up in exchange for uniform local work.

**A LINE SEARCH ON THE FULL RESIDUAL DOES NOT FIX IT, AND THE TABLE ABOVE SAYS
WHY.** Upstream recommends backtracking on the full residual as NPC's
globalisation — `NewtonSolver::ComputeScalingFactor`, a dozen lines,
`miniapps/hdg/navierstokes.cpp` carries one. **It was implemented in MEQ and
measured, and it made every case worse**, including the five that converge
undamped: all of them went to sixty iterations, creeping by about 1% a step.

**THE MECHANISM OF THE CREEP IS A DEFECT IN THE REFERENCE IMPLEMENTATION, AND
UPSTREAM FOUND IT AFTER MEQ REPORTED THE SYMPTOM.** MEQ's first reading — that
`α` collapses because the `ℓ²` merit charges for restoring `(1 − α)` of the flux
and trace residuals — describes the pressure correctly and gets the mechanism
wrong. `NSBacktrackingNewton` accepts on `Norm(rt) < n0`, **a monotone test with
no sufficient-decrease constant**. For Newton on an `ℓ²` merit the direction is
always a descent direction, so `merit(α) ≈ merit(0)(1 − α)` and *any* small
enough `α` passes: upstream swept it and `α = 1.2e-4` still "succeeds",
improving the merit by 1e-4 relative. **So where `α = 1` is rejected the search
does not fail — it creeps**, which is exactly the 1%-a-step crawl MEQ measured.
An Armijo test would reject those steps and fail honestly. Fixing it is open
upstream, and **it would not make MEQ's problem converge**; it would stop the
globalisation disguising its own failure. **KINSOL's `KIN_LINESEARCH` does use
an Armijo condition and fails on the same cases**, which is why the two rows
agree.

**AND THE BLOCK STRUCTURE IS NOT WHAT SEPARATES MEQ'S PROBLEM FROM UPSTREAM'S —
SEVERITY IS.** This file claimed the failure follows from two rows being linear
and one non-linear. Upstream measured the same shape on their own Darcy pedestal
— a full step takes the flux row to 4e-17 and the trace row to 3e-14 with the
whole remainder in the potential row — and there the line search *works*: 3 of 4
configurations fail undamped and converge with backtracking. **What differs is
that a full step IMPROVES their potential residual, 4.7e-01 → 1.7e-01, and
MULTIPLIES MEQ's BY 77.** The block structure is common to both; the size of the
first correction is not. So the recommendation in their §6 is supported there
and wrong here, and MEQ's original explanation was over-general.

MEQ **asked whether §6's baseline was the deleted mode, and it was not** — it
was plain undamped `NewtonSolver` on NPC, so the disagreement is real rather
than an artefact of what was being compared against. On the strength of MEQ's
measurements upstream has withdrawn two §6 claims: *"NPC is not automatically
faster … not fewer floating-point operations"*, and *"reach for [the reduced
trace operator] unless you have a reason not to"* as a general default. §6 now
says the choice turns on what the element-local non-linear solve costs on your
problem.

So the line search was **not kept**, and the comment at the `Globalisation::None`
switch in `GradShafranov.cpp` records why. What a step-length rule would need
here is to damp the non-linear block without un-doing the exact annihilation of
the linear ones — different, and unattempted.

**What DOES work is fixing WHERE the iterate is rather than how far it steps**,
which is the reactive ladder MEQ already has, and this also explains why: Picard
hands NPC a physically sensible `(q, ψ, ψ̂)` instead of `(0, 0, g_D)`, and the
huge first correction never arises. It is the same observation as the warm-start
one below — **under NPC the fields are state, so what they start at matters**,
and `(0, 0)` is the worst place to start.

**Both failures are cured by the reactive ladder, which is what the driver
runs**, so nothing in production is affected:

→ **[M-23](MEASUREMENTS.md#m-23)** — NPC, plain · NPC, `PicardThenNewton`

and both are cured by refinement as well, which is the standing reading of these
sources: *a difficulty measured at one resolution is not a property of the
problem*.

### There is no cheap discriminator, and the obvious one is ANTI-correlated

**ASKED AND ANSWERED BY MEASUREMENT.** The reactive ladder pays for
a failure by running Newton to its cap before handing off. The obvious
improvement is to notice earlier — and the obvious signal, the first Newton
step making the residual worse, **is not a signal at all**. NPC, plain Newton,
cap 60, cold start:

→ **[M-24](MEASUREMENTS.md#m-24)** — case · out · its · `‖r₀‖` · `‖r₁‖/‖r₀‖`

**The largest first-step blow-up in the table converges, and the two cases
where the residual came DOWN at step one both fail.** 1169× at
`k = 1, n = 24` finishes in 26; the barrier's 0.86 and 0.87 run to the cap.
Any threshold that catches the pedestal at `n = 16` throws away `n = 24`, and
no threshold whatever catches the barrier. Under the condensation it is the
same story from the other side — the barrier drops to 4.3e-02 and 1.1e-02 at
step one and still fails.

**And `‖r₀‖` carries no information either, because under NPC it is very nearly
source-INDEPENDENT.** Look down that column: the pedestal and the barrier at
`k = 1, n = 16` both read 1.136e-01, and the pedestal and the layer at
`k = 2, n = 16` both read 1.124e-01 — to four figures, on different sources. At
the cold iterate `(0, 0, g_D)` the residual is dominated by the Dirichlet datum
entering the flux row through `⟨ψ̂, v·n⟩` and the potential row through the
stabilisation, and `(F/R, w)` is a small correction to both. `‖r₀‖` is a
measurement of the boundary data and the mesh, not of the difficulty.

**So the ladder stays reactive and stays triggered by observed failure**, which
is what *Why MEQ's Newton struggles* already concluded from `max|∂F/∂ψ|/λ₁`.
Two independent candidate detectors have now been measured and both are
useless; the difference is that this one is *anti*-correlated, which is worse
than uninformative. **Do not add a predictive trigger without a measurement
that separates these fourteen rows**, and note that the failure is cheap to
observe anyway — see the 1.6 s against 53.5 s above.

### Should `PicardThenNewton` simply be the default? No, and cost is the weakest of the three reasons

Plain Newton against `PicardThenNewton`, both under
NPC, Picard capped at the same 60 iterations, wall clock on an otherwise idle
machine:

→ **[M-25](MEASUREMENTS.md#m-25)** — case · plain · `PicardThenNewton` · slowdown · rel. `max ψ_h`

**Three reasons, in increasing order of how decisive they are.**

**1. It costs 1.4× to 3.1× on everything that already works**, which is real
but is the reason you would trade away for robustness if the other two did not
exist. Note also that this is with Picard *capped at 60*; the cap is a tuning
parameter and *Picard, then Newton* has already measured the handoff to be
**non-monotone** in Picard effort, so 60 working on these fourteen rows is not
a property to build a default on.

**2. It does not cover everything.** The barrier fails under plain Newton and
fails under `PicardThenNewton`, and the `picard 0` in those rows is Picard
itself throwing on the first iteration. A default that is 2× slower and still
needs a fallback is not a default.

**3. IT SILENTLY CHANGES WHICH DISCRETE SOLUTION YOU GET, AND THIS IS THE ONE
THAT SETTLES IT.** Look at the last column on the rows that converge both ways:
6.1e-06, 1.6e-05, and **1.0e-02** on §4.5 at `k = 2, n = 16`. These are
*converged* solves of the same discrete system to `rel_tol = 1e-12`, differing
by up to one percent in `max ψ_h`. Run three ways rather than two, the same
case gives three answers:

→ **[M-26](MEASUREMENTS.md#m-26)** — route · Newton steps · `max ψ_h`

**A spread of 9.4%.** This is the multiple-solution finding recorded under
*The suite is green* — coarse discretisations of these sources
carry more than one solution — meeting the solver from a third direction. It
means changing the default globalisation is **not** a performance decision:
it changes the equilibrium MEQ reports on an under-resolved mesh, which is
precisely the mesh an adaptive run starts from.

So `Globalisation::None` stays the default and the ladder stays reactive. What
the ladder is for is a solve that would otherwise not happen at all, and paying
for it only then is what keeps the answer on the branch plain Newton finds.

**AND THE OTHER ORDERING IS NOT A CHEAPER RUNG THAN PICARD — MEASURED, BECAUSE
IT LOOKED LIKE IT WOULD BE.** `CondenseThenLinearise` converges both cases NPC
fails on, in 13–14 iterations against Picard's ~200, so switching ordering on
failure looks like the obvious cheap rung. It is not: an iteration is not the
currency, and a condensation iteration carries thousands of element-local
non-linear solves.

→ **[M-27](MEASUREMENTS.md#m-27)** — case · NPC plain · `CondenseThenLinearise` · NPC `PicardThenNewton`

**`PicardThenNewton` is between 2.0× and 2.6× cheaper than the ordering switch
on every case either of them cures**, and on the cases neither cures the
ordering switch costs 25× to 30× more to find that out. The ladder MEQ already
has is the right one and nothing needs adding to it.

The two cures also do **not** agree at these resolutions — 2.5e-13 on the
pedestal, but 1.0e-03 on §4.5 at `n = 30` and 9.4e-02 at `k = 2, n = 16` — which
is reason 3 again, and a second argument against reaching for whichever cure is
convenient.

**A corollary for anyone comparing the two orderings**: do it on a *resolved*
mesh. A disagreement at `k = 2, n = 16` is this, not a defect.
`SolverContract::theOrderingsAgreeAndOnlyOneIteratesLocally` runs on Example 5
at `n = 8` where both converge in 4 steps, and asserts 1e-9 — that bound is
only reachable because the discretisation there has one solution.

### KINSOL: the two orderings are exact opposites

Measured over §4.2 at `k = 1, n = 32` and `n = 40`, `k = 2, n = 24`, `k = 3,
n = 32`:

→ **[M-28](MEASUREMENTS.md#m-28)** — NPC · CondenseThenLinearise

**The wiring is ruled out, and sharply.** Under NPC `KIN_NONE` goes through the
same `ShiftedResidual`, the same `DarcyNPCOperator` and the same
`DarcyNPCSolver` and reproduces plain Newton's iteration count **exactly** — 9
against 9, 8 against 8 — so the residual, the Jacobian and the sign convention
are all correct and the fault is in the globalisation alone.
`kinsolAgreesWithNewtonWhereBothConverge` now picks the strategy that converges
for the ordering under test, which is a *sharper* adapter check than before: a
line search may take a different path to the same answer, an undamped KINSOL
may not.

**`KIN_LINESEARCH` and MEQ's own backtracking fail on the same cases**, both
driving an `ℓ²` merit over the whole residual — KINSOL's is `½‖fscale·F‖²` with
`fscale = 1`. **This is one finding, not two.** They fail *differently*, though,
and the difference is instructive: KINSOL applies an Armijo sufficient-decrease
condition and so gives up honestly, while the hand-written one has only a
monotone test and creeps instead. See the mechanism under *Why it fails*.

### The flux is seeded, and it buys a diagnostic and not one iteration

Projecting the guess onto the potential and the trace while leaving the **flux
block at zero** is the obvious thing and is wrong under NPC, where `q` is an
unknown: the guessed state is then inconsistent in exactly the row that couples
them, and `‖r₀‖` goes **up** with a good guess — at `k = 3`, cold 1.771e-01
against warm 2.638e-01. A guess for `ψ` says nothing about `q` without
differentiating it, and a bare `mfem::Coefficient` cannot be differentiated.

`setInitialGuess( mfem::GridFunction const & )` now seeds it, and the projection
is **weighted**:

```
( R q_h, v ) = ( ∇̄ψ_g, v )     on each element,   block := −q_h
```

which is the flux row of (8a) itself rather than an approximation of it — so the
seeded state SATISFIES the row instead of merely being near it. Writing
`q = (1/R)∇̄ψ_g` and interpolating at the nodes is the obvious alternative and is
wrong twice: the closed Gauss-Lobatto basis puts nodes ON element boundaries, so
on FB-A's domain some sit at `R = 0` exactly where `1/R` divides one numerical
zero by another; and nodal interpolation is not what the residual asks for.
**The weight `R` removes the singularity rather than guarding it.** `V_h` is
discontinuous so the solve is element-local — one small dense factorisation each.

**Measured, `‖r₀‖` went from 0.67× the cold value to 83×:**

→ **[M-29](MEASUREMENTS.md#m-29)** — cold · warm

**AND IT DOES NOT CHANGE THE ITERATE, WHICH IS THE FINDING RATHER THAN A
DISAPPOINTMENT.** The whole residual is **affine in `q`** — the flux row, the
potential row and the trace row are all linear in it, and `F` depends on `ψ`
alone — so an undamped Newton step lands on the *same* point whatever `q` it
started from. Measured on the free-boundary limited tokamak of
`FREE-BOUNDARY-PLAN.md` §7.16, seeded and unseeded: **bit-identical residual
histories from iteration 1 through 17**, seven digits, with only `‖r₀‖` moving.
The warm-start iteration counts are unmoved too — 2 against 4 cold, and
4 / 2 / 2 against 4 / 4 / 4 across adaptive cycles.

So the standing suggestion was right about the mechanism and wrong about the
prize: what it buys is an **honest initial residual**, which matters because
`‖r₀‖` is what a reader uses to judge whether a guess was any good, and it
matters on the free-boundary path where a guess is part of the problem
statement. It buys no iterations, and it cannot: **the same argument says
seeding the exterior coefficients would buy none either**, the residual being
affine in `a` as well. What decides the branch is `ψ`, `ψ̂`, `ψ_ax`, `ψ_bnd` and
the scale, and nothing else.

**Where `‖r₀‖` still goes up is where the state genuinely cannot satisfy the
row**: on the coupled free-boundary path the exterior coefficients start at zero,
so the datum on `Γ_h` is absent, and no `q` makes the flux row vanish there. The
limited case reads 5.697e-01 unseeded against 7.424e-01 seeded for exactly that
reason. `aWarmStartCutsTheWorkAndNotTheAnswer` asserts strictly fewer iterations
and an unmoved answer, which is the property that survives all of this.

### One more test moved from the stopping rule to the property

`MillerConvergence::diagnosticExactSolutionOnThePolygon` asserted
`newtonIterations() <= 1` on the Solov'ev source, whose `∂F/∂ψ` is zero, so the
system is affine and step one must be exact. **It was measuring the stopping
rule.** MFEM stops at `max(rel_tol·‖r₀‖, abs_tol)`, so whether an exact step is
also the *last* step depends on where `‖r₁‖/‖r₀‖` — round-off over the residual's
scale — falls relative to `rel_tol`. NPC's residual is the full system rather
than the trace, the ratio moved from just under 1e-12 to just over, and a solve
that had always had a step to spare started taking it:

```
1.466e-02  ->  1.567e-14  ->  2.708e-17
```

Twelve orders in step one, against a target of 1.466e-14 that 1.567e-14 misses by
7%. It now asserts the **drop**, `‖r₁‖/‖r₀‖ < 1e-8`, which the sweep meets
between 4.8e-13 and 3.1e-11 — the floor *degrades* with refinement, because
`‖r₁‖` sits at round-off while `‖r₀‖` shrinks with the mesh — against the
`O(1)` an inexact step on an affine system would give. Rates were unmoved
throughout: 1.995 / 2.999 / 3.999.

## Newton, and the obligation it creates

MEQ uses **Newton**. Both papers use Anderson-accelerated Picard. This is a
deliberate departure, and it has a cost that must be respected everywhere a
source term is written.

`refs/HDG-GradShafranov.pdf` §4 states the design MEQ is reversing: the authors
keep `F` as opaque problem data so the solver "relies only on the discretization
of the toroidal operator `Δ*`", and pay for it by iterating on *every* source,
even one linear in `ψ` that could have been folded into the bilinear form.

Newton makes the opposite trade. It puts `∂F/∂ψ` into the operator — a mass term
`−(∂F/∂ψ)/R` on the potential block, carried through hybridization by
`DarcyForm::GetGradient`. So:

**Every `Source` must supply `dFdPsi`, and every `Profile` must supply `Prime`.**
A profile that cannot differentiate itself cannot be used. This is why
`src/meq/Profiles.hpp` keeps the Hermite cubic — it gives `f` and `f'` from the
same data — and why `SourceTests.cpp` checks `dFdPsi` against a finite difference
of `F`. That test is not box-ticking; it is what stands between a typo in a
derivative and a Newton iteration that quietly fails to converge quadratically
while still converging.

**And a source whose profiles are functions of NORMALISED flux must be a
`meq::NormalisedSource`, not a `meq::Source`** — because `ψ_ax` is then a
functional of the solution and the Jacobian acquires non-local terms through it,
which `dFdPsi` cannot carry and the finite-difference check on `dFdPsi` cannot
see. That is a whole section of its own, below.

There is independent support for the choice, from the free-boundary literature
rather than the HDG one. CEDRES++ (`refs/CEDRES.pdf`) reports that fixed-point
iterations "usually suffer from very slow convergence or even fail to converge,
which made researchers move towards Newton-type methods", and names **vertically
unstable plasmas** as a case where Picard does not converge at all. Lackner's own
review says the same from the other end: his plain Picard, eq. (3), "will
converge to the physically trivial solution ψ ≡ 0 if admitted by the formulation
of the problem".

Three further things from CEDRES++, all of which bite before free boundary does.

**Differentiate the discrete residual, not the continuous one.** A continuous-
level Newton derivative for the plasma-current term exists (Blum 1989) and
CEDRES++ deliberately does not use it: *"there is no theoretical evidence that
this formula holds also for plasma equilibria with boundaries that contain
X-points. In particular the second term on the right-hand side seems to blow up
if ψ reaches a critical point."* They differentiate the Galerkin form instead.
MEQ gets this right for free, since `DarcyForm::GetGradient` differentiates the
assembled operator — but know that the shortcut fails exactly where the physics
is interesting.

**Normalised flux will make the Jacobian non-local.** Their profiles, like
`meq::Profile`, are functions of `ψ_N = (ψ − ψ_ax)/(ψ_bnd − ψ_ax)` on `[0,1]`.
`ψ_ax` and `ψ_bnd` are global functionals of the solution, so `∂F/∂ψ` acquires
terms through them, and those "lead to non-local entries in the stiffness
matrix". Fixed boundary with `ψ = 0` on a known `Γ` does not need the
normalisation and so does not have the problem — but the day the profiles are
driven by normalised flux, `MHDSource::dFdPsi` is incomplete, **and the existing
finite-difference test will not catch it**, because `f()` and `dFdPsi()` would be
missing the same terms.

**That day arrived, and it is done: `ψ_ax` is now an unknown of the non-linear
system.** `HighBetaConvergence` was the standing red test that said it was not,
and it is green. This is the section to read before touching normalised
profiles or starting free boundary, because three cheaper answers were measured
and killed on the way to the one that works.

`refs/GourdainContour.pdf` §V eq (39) specifies profiles the way equilibrium
codes actually do — `p = Σ aᵢ Ψⁱ`, `F² = F²_ax − (F²_ax − F²_bnd) Σ bᵢ Ψⁱ`, with
`Ψ` the normalised flux — and
`refs/CowleyHighBeta.pdf` / `refs/HsuCowleyHighBeta.pdf` construct high-`β_p`
equilibria of that kind asymptotically. `tests/analytic/HighBetaPoloidal.hpp` is
that source.

### What did not work, in the order it was tried

1. **With `ψ_ax` fixed, the profile is inert.** The solution reaches
   `Ψ = 0.0013`, so a peaked `p ~ Ψ^ν` contributes `Ψ^(ν−1) ~ 1e-9` of itself.
   Pressure amplitudes of 1 and 512 gave `ψ ∈ [2.570e-07, 1.259e-03]` — identical
   to every digit. A fixed `ψ_ax` is not a simplification of a normalised
   profile; unless the value happens to be right it is a **different problem**.
   `theFixedNormalisationIsADifferentProblem` still asserts that, as the control.
2. **So "is it a stiff case?" was unanswerable, and the answer looked like no
   for the wrong reason.** Every configuration converged in 1–2 Newton steps at
   reaction ratios `max|∂F/∂ψ|/λ₁` from 0.09 to **2523** — the current hole is
   unsolvable at 26. Nothing was stiff because no non-linearity was switched on.
   It also showed the ratio was being sampled over `[0, ψ_ax]`, a range the
   solve never visits: **that diagnostic needs the solution's actual range.** It
   now gets it, and the answer is **1.9 at `ν = 2` and 13–14 at `ν = 4`** — real
   but not extreme, between the pedestal's 7 and the hole's 26.
3. **Closing the loop outside the solver does not rescue it, and the reason is a
   pole.** Iterating `ψ_ax ← max ψ` with relaxation converges — to
   `ψ_ax ~ 1e-12`, a degenerate fixed point where `ψ` and `ψ_ax` shrink together,
   `Ψ` stays `O(1)`, and the pressure gradient `νA/ψ_ax` runs to 1e12 while the
   solution it drives sits at 1e-12. Mapped afterwards, the cause is plain: the
   outer map `ψ_ax ↦ max ψ` has a **pole beside its own fixed point**. At `ν = 2`,
   `A = 1` the fixed point is 0.3059 and the pole is at 0.2996, six parts in a
   thousand away. Relaxing harder is not a fix for a pole.

**And there is a second solution that is not the equilibrium, which is separate
from the normalisation and cost its own afternoon.** At a *fixed* `ψ_ax` this
equation has a small positive solution and a large one — measured at `ν = 4`,
`A = 1`, `ψ_ax = 0.42` they are 3.0e-3 and 6.2e-1. Only the large one can satisfy
`max ψ = ψ_ax`, so the constraint takes the small branch out of the solution set
— **but not out of the iteration's reach**, and Newton from the Dirichlet datum
walks straight onto it. So this path needs `setInitialGuess()` with a bump of
about the right height, and that is part of the problem statement rather than an
optimisation. It is *not* the trivial-branch trap: unlike every GS-2 §4.2–4.5
source this one does not vanish at `ψ = 0`, the smallest `|F(R,z,0)|` over the box
being 0.2.

### What works: `ψ_ax` inside the residual, as a bordered Newton

`meq::NormalisedSource` is the interface — a `Source` that also has
`setNormalisation()` — and
`GradShafranovSolver::setSource( NormalisedSource &, double )` takes `ψ_ax` as an
unknown. The system Newton closes is

```
R( λ, s ) = 0                     the hybridized trace residual, source normalised by s
G( λ, s ) = s − max ψ_h( λ, s )   = 0
```

with the Jacobian bordered,

```
[  A    c  ]     A = ∂R/∂λ            the existing hybridized Jacobian
[  bᵀ   d  ]     c = ∂R/∂s            dense
                 b = −∂(max ψ_h)/∂λ   sparse
                 d = 1 − ∂(max ψ_h)/∂s
```

**`c` and `b` are the non-local terms**, in exactly the sense CEDRES++ means when
it warns that a normalised profile "leads to non-local entries in the stiffness
matrix": `ψ_ax` is a functional of the whole solution, so a perturbation of the
trace near the magnetic axis moves the source *everywhere*.

**Why it cannot be a rank-one update inside the element blocks**, which is what a
CG code would do. In an `H¹` discretisation `ψ_ax` is one entry of the global
unknown and the Jacobian simply acquires a rank-one term. Hybridization
eliminates flux and potential **element by element**, and a term coupling every
element to the one element holding the axis is precisely what that elimination
cannot represent. The border is where it goes instead — and it costs **one
factorisation and two backsolves**, not a second matrix: solve `A y = R` and
`A z = c`, then `δs = (b·y − G)/(d − b·z)` and `δλ = −y − z δs`.

**`c` is dense and `b` is not, and the asymmetry is structural.** `s` enters
every element's source, so `∂R/∂s` has an entry on every trace dof — one central
difference in a *scalar*, two residual evaluations. `max ψ_h` is one nodal value
of one element, and under hybridization that element's recovered potential
depends only on the trace dofs of its own faces, so `b` has at most `3(k+1)`
entries and the rest are **exactly zero**. Measured, not assumed:
`theAxisSensitivityIsLocalToItsElement` perturbs trace dofs off that element and
`max ψ_h` moves by **0.0000e+00**, against 7.0e-5 for the element's own dofs.

**THE COLUMNS ARE AVAILABLE IN CLOSED FORM, AND THE PARAGRAPH
BELOW — "there is no choice about it" — IS NOW HALF WRONG.** It is right about
`b`, the row, which is a derivative of the *condensed* residual and would need
the elimination's sensitivity. It is wrong about `c`, the column: `s` reaches the
residual **only through the source**, so `∂R/∂s` is the assembly of `∂F/∂s`, and
that is one line of algebra. With `F = g(Ψ)/σ`, `σ = ψ_ax − ψ_bnd` and
`Ψ = (ψ − ψ_bnd)/σ`:

```
∂F/∂ψ_ax  = −( g′(Ψ)·Ψ + g(Ψ) )/σ²
∂F/∂ψ_bnd =  ( g′(Ψ)·(Ψ − 1) + g(Ψ) )/σ²
```

`g(Ψ) = μ₀r²p′(Ψ) + gg′(Ψ)` is what `f()` already builds, and `g′` is one
`Profile::prime()` of each stored profile — the second derivative level that
exists because the rotating source needed it.
`meq::NormalisedSource::normalisationDerivatives()` is the interface, non-pure
so a source that has not implemented it keeps working, and
`the_normalisation_derivatives_are_analytic` pins it against a
Richardson-extrapolated difference at **1.5e-12 to 5e-12** over four values of
`ψ`, both derivatives.

**WHY IT MATTERS IS THE MOVING SUPPORT AND NOT THE SPEED.** A differenced column
floors the iteration — measured on the half-disc, a coupled solve descends to
about **3e-09 and then sits there** for as many iterations as it is given, which
is the difference's accuracy and not the discretisation's. And with
`ConfineToPlasma` on it is worse than a floor: perturbing `s` **moves the edge**,
so the two evaluations have different supports and the difference straddles a
kink. The closed form returns **exactly zero** outside the plasma, asserted.

**IT IS WIRED IN AND THE SIGN WAS CHECKED RATHER THAN ARGUED.**
`GradShafranovSolver::assembleNormalisationColumn()` runs `SourceIntegrator`'s
own quadrature loop over `∂F/∂s`, into the potential block and nowhere else, with
the same `−w F/R` sign — the flux and trace rows carry no `F`.
`setBorderColumn()` keeps the differenced route as a control and
`BorderColumn::Analytic` is the default, falling back silently for a source that
does not supply the derivatives. **NPC only**: under the condensation the
residual is the reduced trace one and this assembly is not it.

**WHAT IT BUYS, MEASURED BOTH WAYS ON ONE PROBLEM.**
`theAnalyticColumnAgreesWithTheDifferencedOne` reads `ψ_ax` and `ψ_bnd` agreeing
to **ten digits** between the two routes — which is what says the sign and the
block are right — and a final residual of **5.70e-16 against 1.69e-13**, a factor
of **296**. The floor is gone. `HighBetaConvergence` is **bit identical**, every
digit of the published table.

**AND IT IS NOT A UNIVERSAL SPEED-UP.** On the `( N + 1 )` coupled system with a
non-linear source the two routes are indistinguishable — 4 Newton steps each,
the same residual to every digit. The assembled column earns its place where the
*difference* is poor: a stiff border, and structurally wherever a moving support
makes the two evaluations straddle the plasma edge.

**AND THE WRAPPER HAD TO FORWARD IT, WHICH IS THE SAME TRAP `setPlasmaSupport`
ALREADY DOCUMENTS.** `CoilAugmentedNormalisedSource` overrides
`normalisationDerivatives()` to delegate — the coil term is amperes and
contributes exactly zero to either derivative — because the base's default
returns `false`, so without the override **every run carrying a `[[coils]]` block
would have fallen back to a differenced column** and nothing would have said so.

**AND THE AXIS POSITION NEEDS NO DERIVATIVE, WHICH IS A THEOREM AND NOT A
SHORTCUT — AND IT IS WHAT MADE OPTION 3 AFFORDABLE.** `ψ_ax` is a *stationary* value, so if the axis position moves with
the solution the chain rule gives `dψ/dλ = ∂ψ/∂λ + ∇ψ·∂(R*,z*)/∂λ` and
**`∇ψ = 0` at an interior extremum**. The second term vanishes identically, so
the largest *nodal* value loses nothing the Jacobian would have used — the
envelope theorem, and a better reason for that definition than the one recorded
above. **It does not extend to an X-point**, where the constraint is `q = 0`
rather than a stationary value and the corner block is `∇q`; see
`FREE-BOUNDARY-PLAN.md` §10.4.

**AND THE SAME THEOREM IS WHY OPTION 3 COSTS NOTHING IN THE JACOBIAN.** With
`ψ_ax` constrained at the LOCATED axis, `G = ψ_ax − ψ_h( x* )` and `x*` moves
with the solution — but `∇̄ψ = R q` and `q_h( x* ) = 0` by definition, so
`∇ψ_h( x* ) = 0` and the position term vanishes there too. No sensitivity of the
root find is needed, and under NPC the row is just the potential shape functions
of `x*`'s element evaluated at `x*`: `( k+1 )( k+2 )/2` entries, exact, one
element. `−e_j` is the special case where `x*` lands on a node. See *Option 3*.

**Both borders are differenced rather than assembled, and there is no choice
about it.** They are derivatives of the *condensed* residual. Assembling them
would need the sensitivity of the element-local eliminations — for `c` the
derivative of each local solve with respect to a parameter of its own source, for
`b` the derivative of the recovered potential with respect to the trace — and
`DarcyHybridization` exposes neither. Differencing the assembled residual is the
same principle CEDRES++ states for the local term: differentiate the **discrete**
residual, never the continuous equation.

**`ψ_ax` IS THE FLUX AT THE LOCATED MAGNETIC AXIS, AND THE NODAL MAXIMUM IS THE
CONTROL.** `setAxisConstraint()` chooses; `AxisConstraint::LocatedAxis` is the
default and constrains `ψ_ax` at a zero of `q_h`.

What follows describes `AxisConstraint::NodalMaximum`, which is what the control
does and is why it is worth having: the largest nodal value differs from the
maximum of the polynomial by `O(h^{k+1})`, both converge to `max ψ`, and the
nodal one is differentiable in a form the border can use directly. See *Option
3* below for what separates the two definitions in practice.

**A backtracking line search is not optional here.** The full step converges for
the mild profiles and wanders for `ν = 4` at amplitude 10, the augmented residual
reading 2.7e-1, 6.6e-2, 2.0e0, 7.7e-1, 2.5e2, 6.7e4, 7.7e6, 3.8e8 before `ψ_ax`
crosses zero and the source refuses the normalisation. That is not the Jacobian
being wrong — the same Jacobian finishes the milder cases at observed order 2 —
it is the equilibrium being a **mountain-pass solution of a superlinear problem**,
where the linearised operator is indefinite and an undamped step leaves the
basin.

**The printed residual is `‖(R, γG)‖` with `γ` frozen at `‖c‖` from the first
iterate.** `G` is a flux and `R` is a trace residual, so the two cannot simply be
concatenated; `‖c‖` is the factor that converts a perturbation of `ψ_ax` into the
units `R` is measured in. Freezing it keeps the history a comparison of like with
like — a `γ` recomputed each step would put the Jacobian's own variation into the
convergence history and manufacture orders out of it. The **control** path gets
the same `γ`, computed the same way, or the two histories would not be
comparable.

### The measurement, and the control that makes it mean something

`k = 2`, `n = 8`, converged `ψ_ax` against the dimensional estimate
`√(νA/λ₁)`:

→ **[M-30](MEASUREMENTS.md#m-30)** — `ν` · `A` · `ψ_ax` · ratio to estimate · `max\|∂F/∂ψ\|/λ₁` · Newton

**AND THE `|dF|/λ₁` COLUMN'S TWO `ν = 4` ENTRIES WERE STALE, CORRECTED
 TO 13.35 AND 14.17 FROM 13.14 AND 14.13.** Caught while checking that
an unrelated change had not perturbed this table: the run at `HEAD` prints 13.35
and 14.17, and so does the changed tree, so the doc had drifted from the code
rather than either moving. **Every `ψ_ax` is unaffected and was bit-identical
throughout**, which is the point — the stale entries are in a sampled diagnostic
and not in an answer, and that is exactly why nobody noticed.

**RE-MEASURED UNDER `AxisConstraint::LocatedAxis`, AND THE `ν = 4`
ROWS MOVED BY 9%.** The `ν = 2` rows are unchanged to six figures; the `ν = 4`
ones read 3.068708e-01 against 2.834510e-01, and the iteration counts went
8/11/10 → 10/9/8. **The old column `ψ_ax − max ψ_h` is deleted rather than
updated**: it was the constraint's own residual under a definition that no longer
holds, and the quantity that replaces it is `normalisationResidual()`, which the
test now asserts directly.

**AND WHY `ν = 4` MOVED IS THE CLEAREST DEMONSTRATION OF THE DEFECT IN THE
TREE.** `sweep()` on the converged `ν = 4, n = 8` field finds **exactly one**
maximum, at `( 1.1112, −0.0083 )` with `|q| = 2.8e-14`, carrying
`ψ = 3.068708e-01`. The largest NODAL value is **3.348506e-01 at
`( 1.1500, 0.0000 )` — 9% ABOVE the field's only magnetic axis** — because
`ψ_h` is discontinuous and that node sits on an element boundary where one
element's polynomial overshoots. **The old constraint was pinning `ψ_ax` to an
overshoot**, on a benchmark that has been in this suite since stage 4. The gap
converges away: 7.63e-02 → 1.94e-03 → 7.63e-04 at `n = 8, 16, 32`.

**AND THE LOCATED VALUE IS NOT MORE ACCURATE ON A COARSE MESH, WHICH IS WORTH
SAYING BECAUSE IT WOULD BE THE COMFORTABLE CLAIM.** The `n = 32` answer is about
2.841e-01, and at `n = 8` the located reading (3.069e-01) is *further* from it
than the nodal one (2.835e-01) — the overshoot happened to compensate. **The
argument for the definition is meaning, not accuracy.**

At `ν = 2` the source is linear in `ψ` and the problem is a linear eigenvalue
problem, which is why `ψ_ax` is converged in the mesh to six figures between
`n = 8` and `n = 16`. At `ν = 4` it is genuinely non-linear and under-resolved at
`n = 8`, which is what the paragraph above measures.

**`Normalisation::Decoupled` is the control, and it is what a test of this can
actually assert on.** It is the same solver, mesh, guess, line search and
stopping rule with exactly three quantities zeroed — `c`, `b` and `d − 1` — so
the step in `ψ_ax` reduces to `ψ_ax ← max ψ_h` and the trace step is a Newton
step that does not know `ψ_ax` is about to move. That is `ψ_ax` outside the
residual, done as favourably as possible. Measured at `ν = 2`, `A = 1`:

→ **[M-31](MEASUREMENTS.md#m-31)** — residual · iterations

**It does not converge slowly. It does not move.** That is the test the ROADMAP
asked for when it said the work "needs a test that can *see* the missing terms",
and it is needed for the reason recorded under *A wrong Jacobian is invisible to
a convergence table*: `SourceTests`' finite-difference check on `dFdPsi`
structurally cannot see this one, because `f()` and `dFdPsi()` are both evaluated
at whatever normalisation is set and agree with each other however wrong it is.

### The trap that cost the most, and it is in MFEM's local solves

**`DarcyHybridization` captures the element-local Newton's initial guess at
`FormLinearSystem()` time and keeps it for the life of the reduced system.**
`EliminateVDofsInRHS` copies the flux and potential blocks into `darcy_u` and
`darcy_p`, and every local solve in every subsequent residual evaluation,
gradient and recovery starts there — however far the trace has since travelled.
`ComputeSolution()` does **not** use the output vector as a guess, so seeding
that is inert.

Left alone this is a performance problem on an ordinary Newton path and a
**correctness** problem on this one. Measured at `ν = 4`, amplitude 10, `n = 16`,
`k = 2`: 40,000 to 60,000 element-local iterations per outer step, most of them
hitting the cap of 100 — and a local solve that ran out of iterations returns
whatever it had reached, which is not a function of anything. Differencing
`ψ_ax` by 9e-6 then moved `max ψ_h` from 0.8961 to **2.04**, and on the next step
to **3.84**. The corner of the border read 1.6e5 where it should read about 1,
the step in `ψ_ax` collapsed to 1e-8 against a constraint residual of 3e-3, and
the iteration stalled — looking exactly like a singular border and being nothing
of the kind.

`solveWithNormalisation()` therefore re-forms the system from the recovered state
once per accepted step (`formSystem()`, factored out of `prepare()` for this).
Every local solve is then within an iteration or two of its answer, so it
converges, so it is continuous in `ψ_ax`, so the difference is a derivative. The
same run then finishes in five steps.

**Anything else that differentiates a hybridized residual by differencing it will
hit this.** It is the reason a solver-level finite difference is not simply "one
more residual evaluation".

### The production normalised source, and what a working Newton looks like

`meq::NormalisedMHDSource` is the production source built on two `meq::Profile`s
in normalised flux. **It is wired**: `makeNormalisedSource` serves it as well as
the rotating source, so `[source] Type = "mhd"` with `Normalised = true` reaches
the bordered Newton. `makeSource` **throws** on a normalised configuration rather than
returning one, since a `NormalisedSource` is-a `Source` and the plain path would
converge with `ψ_ax` frozen at the guess.

**What working Newton looks like.** CEDRES++ Table 2, on 577k unknowns: relative
residual `2.7e0 → 9.2e-2 → 1.8e-3 → 5.3e-6 → 3.9e-12` in five iterations. That is
the shape stage 4 should produce. A run that grinds down linearly means the
Jacobian disagrees with the residual.

MEQ's own, measured at `k = 3`, `h = 0.1`, 832 trace dofs:

```
   it            ||r||    ||r||/||r_0||    order
    0     1.121994e+01     1.000000e+00        -
    1     4.085930e-02     3.641669e-03        -
    2     5.744224e-04     5.119658e-05    0.759
    3     1.411243e-07     1.257800e-08    1.949
    4     4.142741e-14     3.692303e-15    1.810
```

### The Solov'ev coefficients were wrong twice, and are now checked

`Soloviev.hpp`'s `nstx()` does **not** use the coefficients printed in
`refs/HDG-GradShafranov-Adaptive.pdf` eq. (22c); `nstxAsPublished()` keeps those.
It also does not use this file's *first* correction of them. Both errors are
recorded in the fixture's header, and the shape of the second is the more
instructive.

**Error one: the published set satisfies none of Cerfon & Freidberg's twelve
constraints.** `refs/CerfonFreidberg.pdf` §IX eq. (28) gives them for an up-down
asymmetric single null. At the four points where `ψ` must vanish the printed
coefficients give `−7.5e-3`, `+3.9e-4`, `+2.9e-2`, `−9.8e-3` — the third being
11% of the axis flux.

**Error two: `α` is not `δ`.** The re-solve substituted `sin α = δ` into eq.
(11)'s `N₁ = −(1+α)²/(εκ²)` and `N₂ = (1−α)²/(εκ²)`. But those mean `α` itself,
which for `δ = 0.35` is `arcsin(0.35) = 0.3576`. `N₃ = −κ/(ε cos²α)` is immune,
since `cos²α = 1 − δ²` either way — so the error hit two of three conditions and
nothing else at all.

Settled by differentiating C&F's model surface eq. (9) directly at `τ = 0, π,
π/2`: `N₁ = −0.5907049043`, `N₂ = +0.1322804125`, `N₃ = −2.9220542041`, each
matching the `α` form exactly and the `δ` form by 1.1% and 2.4% at the two
equatorial points.

**The moral, which is the transferable part: checking a solve against the
formula it used cannot detect a misread formula.** The first correction verified
the constraints, and passed, because it verified them against its own wrong `N`.
Only an independent quantity catches that — here the curvature of the surface the
coefficients are supposed to reproduce. The corrected set matches it to 2e-11,
3e-10 and 2e-12.

**And they are no longer asserted nowhere.**
`tests/convergence/SolovievGeometryConvergence.cpp` evaluates all of C&F's
conditions on every set in the fixture, and `nstxUsesTheCorrectAlpha` checks both
readings of `α` side by side so a silent revert fails rather than passing both
ways. Any statement elsewhere in the tree that the coefficients are unchecked is
out of date.

The reason this needed a dedicated test at all: every `ψ_i` is `Δ*`-harmonic, so
**any** coefficients leave `F`, `Δ*ψ` and every convergence rate exact. Nothing
else in the suite can see a wrong one. The absolute-error ceilings move, which is
why they are recorded beside each `checkOrder()` call with the value they were
set from.

`ExtensionConvergence` still takes `Γ` to be the interior surface `ψ = −0.03`
rather than `ψ = 0`: with correct coefficients `ψ = 0` *is* the separatrix, which
passes through an X-point — a **corner** of `Γ`, where both transfer-path
families give out and the Cockburn–Solano analysis does not reach.

**A tooling warning that has now cost time three times, and the third is the
worst.** `pdftotext` silently drops this paper's minus signs, and an `ε` in
another. **Read the rendered page.**

**AND ON `refs/CEDRES.pdf` IT DOES SOMETHING DIFFERENT AND MORE DANGEROUS: IT
DROPS RADICALS AND DISPLACES EXPONENTS.** Checked against the page rendered at
900 dpi, on the kernel of their boundary form (3.5):

| quantity | the page | `pdftotext -layout` |
|---|---|---|
| `M`'s denominator | `2π(r₁r₂)^{3/2}` | `2π(r1 r2 ) 2`, with the `3` on the **next line** |
| `k` | `k = √( 4r_j r_k / ((r_j+r_k)² + (z_j−z_k)²) )` | **the `√` is absent entirely** |
| `δ_±` | `√( r₁² + (ρ_Γ ± z₁)² )` | `r12 + (ρΓ ± z1 )2` — **no `√`** |
| (3.1)'s weights | `∫ψ²r`, `∫(∇ψ)² R^{−1}` | exponents displaced to the line above |

The minus signs came through correctly this time, which is the point: **the
failure mode is the tool's and it is not the same failure twice.** Two of those
four are silently fatal. Losing the radical on `k` turns the **modulus** into
the **parameter** — `E(k)` against `E(k²)`, which is the classic elliptic-integral
error and converges to a wrong answer rather than failing; and `(r₁r₂)^{3/2}`
read as `(r₁r₂)²` or `(r₁r₂)^{1/2}` changes the weight the whole of
`FREE-BOUNDARY-PLAN.md` §3.2 turns on. **Render the page for any equation
carrying a radical, a fractional exponent or a weight**, which in this subject
is most of them.

### A wrong Jacobian is invisible to a convergence table

The single most useful thing measured in stage 4, and the justification for
asserting on Newton's *order* rather than merely on the rates.

Perturbing `∂F/∂ψ` by **+5%** and re-running the whole Example 5 study leaves
**every error and every convergence rate unchanged to six significant figures**.
The discretisation is untouched, because Newton converges to the same discrete
solution whatever Jacobian carried it there. What changes is only the path:
observed order drops to exactly `1.000`, and it takes 10 iterations instead of 4.
At +50% Newton diverges to NaN; drop the Jacobian mass term entirely and McCarthy
takes 87 iterations instead of 1, at a constant contraction of 0.7711 per step,
with the errors again unchanged.

So a rate table cannot see a Jacobian error at all. Three things can, and all
three are in the suite: the finite-difference check on the assembled Jacobian
(4e-11 relative, at the `O(step²)` floor), the assertion on observed Newton
order, and the McCarthy rung, whose affine source must finish in exactly one
step.

### Why MEQ's Newton struggles where other Newton solvers do not

**READ THIS FIRST: THREE OF THE FOUR "STIFF" SOURCES WERE MERELY
UNDER-RESOLVED.** The rest of this section is kept because its falsified
hypotheses are worth not repeating, but its premise is largely wrong. Raw
Newton, undamped, cold start, iterations to converge over `n = 16, 24, 32, 48`:

→ **[M-34](MEASUREMENTS.md#m-34)** — `k = 1` · `k = 2` · `k = 3`

**§4.2, §4.3 and §4.5 are ordinary problems on a resolved mesh** — 7 to 17
iterations — and both refinement paths cure them independently: `h`-refinement
takes §4.2 at `k = 1` from 42 to 9, and `p`-refinement takes `n = 16` from 42 to
10 without touching the mesh. `k = 3, n = 48` on the pedestal takes **7**. The
"MEQ is doing the easier problem and finding it harder" red flag below was
raised from a benchmark run at a single under-resolved point, `k = 1, h = 0.05`.

**§4.4 is the exception and is genuinely unsolved.** It fails at *every* order
and every mesh tried, including `k = 3, n = 48`, spending 1.8M element-local
iterations to do it. Refinement does nothing, and neither does Picard.

**And the reason is measured, not guessed: the Jacobian's reaction term has swept
past ~26 eigenvalues of the operator it is added to.** `∂F/∂ψ` for eq. (26)
ranges over `[−579.5, +565.7]`, against a first Dirichlet eigenvalue of
`λ₁ = π²(1/w² + 1/h²) = 22.3` on the benchmark box. The linearised operator
`−∇̄·((1/R)∇̄·) − (∂F/∂ψ)/R` is therefore strongly indefinite, and the continuous
problem is **multi-valued**. The driver is the `c₃(1 − e₂)cos(c₄ψ)` term at
`c₃ = −18`, `c₄ = 10π` — the very feature that empties the core of current —
whose derivative carries `c₃c₄ ≈ 566`.

That is why refinement is powerless: there is no discretisation error to remove.
Ratio against `λ₁`, by amplitude:

→ **[M-35](MEASUREMENTS.md#m-35)** — `c₃` · 0 · −2 · −4 · −9 · **−18**

Note the pedestal itself already sits at 7 and converges, so exceeding `λ₁` is
not by itself fatal; 26 is. **This supersedes an earlier reading of §4.4 as a
trivial-branch problem.** `F(R, 0) = 0` is true and the trivial branch is real —
see `CLAUDE.md`, *Traps* — but it is not what defeats the iteration here, since the runs above
carry non-homogeneous ramp data that keeps `ψ` away from zero and fail anyway.

Do not file §4.4 under stiffness.

**§4.4 does have a solution, and continuation reaches it.** Adaptive
continuation in the added term's amplitude — halve the step on failure, grow it
by 1.3 on success — walked `c₃` from 0 to −18 at `k = 2, n = 32` in **9 solves
with 2 retreats**, the last taking 8 Newton iterations. Uniform steps do not do
it (ten of them stall at `c₃ = −10.8`), so the step control is the point, and
there is **no limit point on the branch** — the step never had to fall below
1e-3, it simply had to shrink over the last sixth.

**THIS MUST NOT GO INTO THE DRIVER, AND THE REASON IS NOT TASTE.** What was
ramped is `c₃`, a constructor argument of the `CurrentHole` *test fixture*. The
`Source` interface exposes `f( R, z, ψ )` and `dFdPsi( R, z, ψ )` and nothing
else — **there is no amplitude parameter in it and no way to recover one from a
black-box `F`**. The continuation is unavailable to a driver on principle, not
merely inadvisable. The black-box analogue, `F_λ = λF`, *is* expressible but is
**untested and a different homotopy**: at `λ = 0` it degenerates to the harmonic
problem, where the path above starts from the converged *pedestal*.

**And it cannot be predicted that continuation is needed.** `max|∂F/∂ψ|/λ₁` is
computable black-box, `dFdPsi` being mandatory — but the pedestal converges at 7
and the hole fails at 26 (two points, no threshold), and the ratio needs the
range of `ψ`, which is not known before solving. A detector calibrated on that
would be fitting noise. **The driver gets a reactive ladder, never a predictive
one**; see also *There is no cheap discriminator*, where a second candidate was
measured and turned out to be anti-correlated. If continuation is ever wanted
generally the prerequisite is a `Source` that can parameterise itself, which is
an interface change to argue on its own merits.

So the remedy for a production run is **resolution** — the adaptive loop rather
than a globalisation. What globalisation buys is the **coarse start** an
adaptive loop necessarily begins from, which is exactly where raw Newton fails.

What survives below: hybridization really does put a non-linear solve inside
every element elimination under `CondenseThenLinearise`, and that is why an
under-resolved §4.2 grinds where a CG code would not. What does **not** survive
is the conclusion that this makes the problems unreachable — it makes an
under-resolved discretisation expensive, and NPC removes the local solves
without removing the difficulty.

**Two hypotheses below are falsified and one is now known to be incomplete.**

Three Grad–Shafranov codes solve this equation by Newton and report it robust:

| | discretisation | nonlinear systems per Newton step |
|---|---|---|
| Serino, Tang, Tang, Kolev & Lipnikov (`refs/MFEM-GS-Newton.pdf`) | CG, `ψ ∈ H¹`, MFEM | **1, global** |
| CEDRES++ (`refs/CEDRES.pdf`) | CG finite elements | **1, global** |
| FreeGSNKE | 4th-order finite differences, Newton–Krylov | **1, global** |
| **MEQ** | **hybridized HDG** | **1 global + one per element** |

**In a CG or finite-difference discretisation `ψ` is one global unknown vector
and `F(ψ)` enters only the global residual**, so Newton linearises once — Serino
et al. eq (3.1) and eq (3.7) are the whole nonlinear structure. **Hybridization
changes that.** Static condensation expresses `(q, ψ)` element by element in
terms of `ψ̂`, and when `F` depends on `ψ` that elimination is *itself a
nonlinear solve, once per element per residual evaluation* — none of them
globalised, and any one failing poisons the whole residual. That is what
`el: N not convered in 100 iters` is, and it is why globalising the outer
iteration does nothing; see *On SUNDIALS*.

**Two comfortable explanations are wrong, and the papers say so directly.**
Newton is not the problem: Serino et al. built their Newton solver precisely
because "conventional Picard-based solvers fail to converge" on the Taylor
state, and report the residual reaching 1e-6 "in a small handful of iterations".
Free boundary is not the problem either — they note the fixed problem is
"significantly easier", so **MEQ is doing the easier problem and finding it
harder**. It also explains the GS papers' choice: keeping `F` as opaque problem
data leaves **every local solve linear**, which in a hybridized method is
coherence rather than fastidiousness.

**AND THE STRUCTURAL STORY WAS ITSELF FALSIFIED.** This section predicted that
applying Newton to the full `(q, ψ, ψ̂)` system — `refs/HDG-NPC-2.pdf` §2.6,
which is what `NonlinearOrdering::NPC` is and what MEQ now defaults to — would
fix the stiff sources. It does not. NPC has **no element-local non-linear
iteration at all**, confirmed at `GetNumLocalNLIterations() == 0`, and it loses
on precisely the under-resolved cases this section is about. What it buys is
uniform local work and 3.7×–4.3× of wall clock, not robustness. **The
element-local non-linear solves were never the cause.** See *The NPC port*.

**What follows, in order of how targeted it is.**

1. **The local solver was never chosen.** `SetLocalNLSolver` offers `Newton`,
   `LBFGS` and `LBB`, and MEQ hardcoded `Newton` — undamped, on exactly the
   problems that were failing. Measured on §4.2 at `k = 1, h = 0.05`: Newton 42
   outer iterations, LBFGS 36, **LBB 25**. `setLocalSolver()` exposes it. **It
   is inert under NPC**, which has no local non-linear solve to configure.
2. **The control that pointed at the local solves, and why it did not prove what
   it looked like it proved.** The variable was isolated by **Picard on MEQ's
   own linear path** — identical mesh, spaces, hybridization and `τ`, but `F`
   evaluated at the previous iterate and handed to `setSource( Coefficient & )`,
   which makes the potential block linear and every local elimination a linear
   solve. On §4.2 at `k = 1, h = 0.05`, the case that fails: Newton **fails**,
   undamped Picard **stalls at 3.5e-1**, and Picard at `ω = 0.5` reaches
   **2.8e-8 and falling**. So the mesh is fine, the discretisation is fine and
   the problem *is* solvable there.

   **But that control differed in GLOBALISATION as well as in local linearity**,
   and the later globalisation cross says which half was doing the work — the
   damping, not the local linearity. It is why the `ω = 1` row matters: undamped
   Picard stalls, which is the weakness CEDRES++ and Serino et al. both report
   and the reason the GS papers use **Anderson**-accelerated Picard rather than
   plain. **A diagnosis, not a recommendation**: relaxed Picard took 200
   iterations to reach 2.8e-8 where Newton takes 42 on the mesh Newton manages.

**AN ORDERING DOES NOT FIX A STIFF SOURCE, AND THAT PREDICTION HAS NOW BEEN
FALSIFIED TWICE.** The argument was that condensing before linearising is the
*cause* of MEQ's stiff-source trouble and that reversing the order would fix it.
Reversing it was confirmed genuinely in effect by
`GetNumLocalNLIterations() == 0` and was **strictly worse on every stiff
case**. NPC — which is the canonical version of the same
idea and is MEQ's default today — reproduces that verdict; see *The NPC port*.
**Do not expect an ordering to fix a stiff source.** It is the standing warning
against reasoning from structure to robustness in this file.

**A measurement technique worth reusing: cross-linearisation residuals.**
Evaluate the reduced residual twice at one trace, relinearising *at a different
trace* in between. Under condense-first the difference is exactly 0.000e+00 at
every local stiffness — it must be, since the local problems are solved to
convergence and nothing about how they were reached survives. Any nonzero
answer says the operator carries hidden state and no Newton can converge on it.
That is what condemned the deleted mode, at 149% of the residual's own value.

**AND IT FOUND A LIVE DEFECT IN THE ORDERING MEQ STILL SHIPS AS THE BACKUP.**
Under condense-first the residual is a perfect function of the trace, but the
assembled gradient disagrees with a central difference of it **by a factor of
100** once the source's local width reaches `σ² ≤ 0.01`. That is the
element-local solves hitting their 100-iteration cap: a residual computed from
an unconverged local solve is not smooth, so the difference is meaningless and
the Jacobian may or may not be. **It is the measured explanation for the
pedestal's wandering iterations at `k = 1, n = 16` under
`CondenseThenLinearise`**, and it is the same defect that corrupts a differenced
border — see the seed note under `CLAUDE.md`, *Traps*. Refinement cures it, because a
resolved local problem converges inside its cap.

**Globalising the outer iteration makes the local solves worse, not better**,
which is the other half of *On SUNDIALS*. On §4.2 at `k = 1, n = 16` under
condense-first: plain Newton converges in 31 with 133,168 element-local
iterations; `KIN_LINESEARCH` fails at 45 having spent **1,381,527**; `KIN_NONE`
fails at 19. A line search chooses how far to move the trace and cannot make
the local problems at that trace well posed.

3. **Picard, keeping the local problems linear — implemented, as a bridge.**
   `Globalisation::AndersonPicard` is `KINSolver(KIN_FP)` over the fixed point
   `ψ^{k+1} = G(ψ^k)`, where `G` freezes `F` at the previous iterate, puts it on
   the right-hand side and does one linear solve. The papers' own method, and it
   works where Newton grinds. On §4.2 at `k = 1, h = 0.05`:

   | method | outcome |
   |---|---|
   | Newton | 42 iterations |
   | Picard, undamped | stalls |
   | Picard, `ω = 0.5` | 248 iterations |
   | **Anderson, depth 1, undamped** | **162 iterations** |
   | Anderson, depth 2 and above | fails |

   **Two surprises, both defaults now set from measurement rather than from the
   papers.** Plain Picard *needs* damping and Anderson does not — so
   `setPicardDamping` defaults to 1.0, which is right for one path and wrong for
   the other. And **HDG-GS-1's `m = 2` fails here** where `m = 1` converges, so
   `setAndersonDepth` defaults to 1; whether that is this fixed point's
   conditioning or KINSOL's implementation is not established, and raising it
   expecting the papers' behaviour will not work.

   It is a **robustness route, not a faster one**: 162 iterations against
   Newton's 42, each one a full linear solve.
   `andersonPicardReachesTheSameSolutionAsNewton` pins that both reach the same
   discrete solution — `ψ` agreeing to seven figures — which is what makes it an
   alternative rather than a different problem.
4. **Continuation in the source amplitude**, which also addresses the trivial
   branch.

### Picard, then Newton — the route for a coarse mesh

**This is what reaches the hard cases without refining them, and it is
`Globalisation::PicardThenNewton`.**
Anderson-accelerated Picard walks the iterate into Newton's basin; plain Newton
takes it from there and supplies the quadratic endgame Picard structurally cannot.
Measured, with `setInitialGuess()` seeding the ramp on every row:

→ **[M-36](MEASUREMENTS.md#m-36)** — case · Newton alone · Picard · **Newton from Picard**

Three of the four unreachable cases become reachable. §4.3 at `k = 1, h = 0.05`
is the one to quote: nothing else in the solver touches it **at that
resolution** — plain Newton fails at 60, a line search fails at 24 having spent
5.4M element-local iterations, and linearise-first aborts — and the handoff
finishes in four Newton steps,
`8.3e-01 → 1.5e-03 → 1.0e-04 → 7.4e-07 → 3.7e-11`, agreeing with Picard's own
answer to **4.7e-10**. That agreement is what makes it a handoff rather than a
change of problem, and `picardThenNewtonRecoversQuadraticOrder` asserts it.

**Be precise about what this is worth, because refinement reaches the same three
cases.** Raw Newton solves §4.2, §4.3 and §4.5 perfectly well once resolved — see
the table under *Why MEQ's Newton struggles* — so the handoff is **not** the only
route to them, and it is not the route to prefer when refining is available. What
it is for is the **coarse mesh**: an adaptive run must solve on its initial mesh
before it has an estimator to refine with, and that first solve is exactly the
under-resolved regime where raw Newton fails. That is a real job, and it is a
narrower one than "the route for stiff sources".

**It does not rescue §4.4.** The current hole fails under Picard, under the
handoff, and at every order and mesh up to `k = 3, n = 48`. Its problem is the
trivial branch, not the iteration.

**Picard's job is not to solve the problem.** §4.5 converges at both orders from
a Picard state that never met its own tolerance, so stage 1 stopping short is an
expected outcome, not an error, and `solveByPicardThenNewton()` swallows that
throw deliberately. This is a globalisation, not a two-solver pipeline.

**Do not replace the tolerance with an iteration budget.** The handoff is *not
monotone* in Picard effort — on §4.5 at `k = 1`, budgets of 400 and 3 converge
while 40 and 10 fail, and on §4.3 at `k = 1` a budget of 3 diverges to `1e4`. A
budget tuned on one mesh will betray you on the next. Picard's own tolerance is
the trigger that worked wherever it was reached.

It is not cheap: stage 1 spent 122 to 290 iterations, each a full linear solve.
Reach for it when `Globalisation::None` fails, not before.

**A caution on reading the printed order.** Taking the best observed order over
any triple of a short or non-monotone history manufactures values of 3.2, 3.6 and
9.36 out of runs that are not converging at all. Only a monotone tail supports an
order claim; the 2.02, 2.01 and 1.96 above are those.

**And it exposed a latent defect, now fixed.** `setGlobalisation()` reset
`prepared` but not `built`, while `buildForms()` branches on
`usesNonlinearForms()` — which reads `globalisationChoice` — to decide whether the
potential block goes on the linear or the nonlinear form. Switching a live solver
between a Picard path and a Newton one therefore reused the other path's blocks.
`setNonlinearOrdering()` and `setLocalSolver()` beside it always did reset
`built`. Latent until now because every caller built a fresh solver per
globalisation; `PicardThenNewton` switches twice inside one `solve()`.

### On SUNDIALS

`mfem::KINSolver` **derives from `mfem::NewtonSolver`** (`linalg/sundials.hpp`)
and is reachable through `setGlobalisation()`; `SetJFNK` and `EnableAndersonAcc`
come with it. SUNDIALS 7.5.0 is built in at `../sundials/cuda-install`.

**IT IS NOT A DROP-IN, AND THE DIFFERENCE IS SILENT.** `NewtonSolver::Mult(b,x)`
forms `R = oper(x) − b`, while **`KINSolver::Mult` declares its first argument
without a name** and solves `oper(x) = 0`. MEQ's trace right-hand side is not
zero, so handing it straight to KINSOL converges — to the solution of a
different problem. `ShiftedResidual` in `GradShafranov.cpp` is the adapter, and
it reproduces `NewtonSolver`'s residual exactly, which is what makes any
comparison between the two paths mean anything.
`kinsolAgreesWithNewtonWhereBothConverge` would catch it being dropped or its
sign flipped. For the NPC-era behaviour see *KINSOL: the two orderings are exact
opposites*.

**Globalising the outer iteration does not globalise the inner ones**, which is
the structural point and is why a line search rescues nothing under
`CondenseThenLinearise`. The failure there is
`el: N not convered in 100 iters` — MFEM's **element-local** non-linear solve,
one per element per residual evaluation. A line search chooses how far to move
the *trace*; it cannot make the local problems at that trace well posed, and
KINSOL never sees them. What would: damping the **local** solves, which is
MFEM's to offer; continuation, so each solve starts from the previous answer; or
**Picard on the outer loop**, which leaves every local problem linear. That last
is what both papers do. **NPC removes the whole question** by having no
element-local non-linear solve at all.

**`MFEM_USE_LAPACK` APPEARS TO FIX THE PEDESTAL AND DOES NOT — IT TIPS A
MARGINAL ITERATION BY CHANGING THE ROUNDING.** Against a LAPACK build §4.2 at
`k = 1, h = 0.05` converges in 42; against `install-nolapack`, identical but for
that flag, it fails at 60. **But at `h = 0.0333` both builds converge in 23, the
same number.** Only the marginal mesh moves. There is no mechanism for it to be
anything else: MEQ sets `LPrecType::LU`, both implementations partial-pivot on
the same rule, and `Factor( m, TOL = 0.0 )` is `dgetrf`'s singularity condition
exactly. What is left is blocked BLAS-3 against unblocked scalar loops —
identical arithmetic, different summation order, `O(1e-16)`. **And the BLAS here
is threaded MKL**, whose reduction order depends on the thread count, so whether
that local Newton converges is machine- and environment-dependent. Do not treat
42 iterations as reproducible. An intermediate claim in this file that LAPACK
had fixed a third of the globalisation problem was wrong.

**`MFEM_USE_EXCEPTIONS` is enabled** and does what the driver's exit codes need:
`MFEM_ERROR_THROW` is the default error action and `mfem::ErrorException`
derives from `std::exception`, so **MEQ needs nothing of its own** — a plain
`catch ( std::exception const & )` catches it. §4.4 at `k = 1, n = 16` reports a
failure with a usable iteration count rather than taking the process down with
SIGABRT, so a driver can return exit code 2 rather than dying.

**A caveat worth keeping in view**: the throw unwinds out of the middle of MFEM,
and the objects are left as the throw found them. MEQ's paths construct a fresh
solver per solve and so do not care, but **do not assume a
`GradShafranovSolver` is reusable after a caught `ErrorException`.**

### The suite is green, and the last red one was a mesh chosen for a dead solver

The NPC port dropped the whole suite from
872 s to 499 s — the element-local non-linear solves that cost `PedestalConvergence`
its time are gone, 213 s against 572 s, and `HighBetaConvergence` 3.3 s against
17 s. It left one red assertion, and that is now fixed rather than tolerated.

**What was red**: `internalLayerSelfConverges`, GS-2 §4.5 at `k = 1`, on the
**coarsest** mesh of its self-convergence sweep, `n = 24`.

**WHAT IT ACTUALLY WAS: A MESH SEQUENCE CALIBRATED AGAINST AN ORDERING MFEM HAS
DELETED.** `meshesFor()`'s table recorded `h = 0.0333: 26 it` for that point,
measured under the old default. The discretisation reason was there all along and
nobody had connected it: §4.5's ridge is **about 0.025 wide in space**, and
`n = 24` gives `h = 0.0333` — **the feature was thinner than a cell**. The
self-difference at that point was measuring the approach to the asymptotic
regime, not the rate, which is the same objection the file already records
against starting `k = 2, 3` coarser than 16.

Moving the `k = 1` sequence to `{ 48, 96, 192 }` fixes both at once, and the
rates improve rather than being relaxed to fit:

→ **[M-37](MEASUREMENTS.md#m-37)** — §4.5, `k = 1` · `ψ` · `q`

Design order for `k = 1` is 2, and the corner control caps `q` at about 2.2, so
the finer sequence sits at design order in `ψ` and at the cap in `q` — which is
what the study is for. Plain NPC converges at all three points, 12/13/12.
`PedestalConvergence` went 216 s → 262 s, and the whole suite only 509 s →
**522 s**: the 60 wasted iterations at `n = 24` paid for most of `n = 192`.

**And the failure had two distinct characters, which is worth keeping**:

→ **[M-38](MEASUREMENTS.md#m-38)** — `n` · `h`/0.025 · plain NPC

The `n = 32` row is not divergence: the residual descends 7.99e-02 → 2.69e-03 in
ten steps and then *orbits*, period 4 to a mean relative mismatch of 0.041
against 0.32–0.44 at every other period, eleven orders short of its target.
Newton with a stable periodic orbit, which is a different animal from the
wandering at `n = 24` and from anything else recorded in this file.

**NONE OF IT IS THE DISCRETISATION.** Every route that converges reaches the same
discrete solution — `CondenseThenLinearise` and NPC-with-`PicardThenNewton` agree
to between 2.9e-11 and 2.3e-09 in L2 at every mesh from 24 to 192, and their
self-differences are **bit-identical**, so the rate does not depend on the route.
Both failing points are reachable: `PicardThenNewton` takes 11 Newton steps at
`n = 24` and **4** at `n = 32`.

**AND THE BASIN FOLLOWS THE FEATURE, NOT THE DOF COUNT.** Marking on `|∂F/∂ψ|`
at the datum ramp — an a-priori marker using nothing the solve is not given — and
refining that band once from `n = 24`:

→ **[M-39](MEASUREMENTS.md#m-39)** — trace dofs · plain NPC

Inside the basin on **54%** of uniform `n = 48`'s dofs, and with 21% more than
the uniform mesh that fails. The band must be generous, though: marking at half
the peak takes only 102 of 1152 elements and still fails, because the datum is
only a proxy for where the ridge sits. This does not go into the
self-convergence study — a rate needs a uniform `h` — but it is the right move
for anyone who needs a coarse solve on a localised-feature source.

**THE TRANSFERABLE LESSON, AND IT IS THE ONE THIS FILE KEEPS RELEARNING**: a
"coarsest usable mesh" is a property of the **solver** as much as of the
benchmark. That table was calibrated against a solver that no longer exists, and
the red assertion was the calibration going stale rather than anything about
MEQ. Re-measure it after any change to the ordering or the globalisation.

**`pedestalConvergenceIsAResolutionThreshold` ASSERTS THE CURES, NOT THE KNIFE
EDGE**, and nothing is entitled to assert the edge in either direction: 42
iterations against a LAPACK build and a failure at 60 against
`install-nolapack` is one flag's difference in threaded-MKL reduction order, so
an assertion on it is an assertion about this machine today. It asserts the two independent refinement cures, each with a
factor of four in it:

→ **[M-40](MEASUREMENTS.md#m-40)** — `k = 1, n = 16`

That two independent paths both cure it is what identifies the cause as
under-resolution. §4.4 is the control: it fails at every order and mesh tried,
because its trouble is that the problem is multi-valued. The threshold itself
is still asserted, as `iterations(n = 16) ≥ 2 × iterations(n = 32)`, which holds
whichever side the rounding falls on.

**THESE COARSE DISCRETISATIONS CARRY MORE THAN ONE SOLUTION, and it is measured
from three directions now.** `andersonPicardReachesTheSameSolutionAsNewton`
cannot gate one mesh at 1e-6; over a sweep it reads 9.1e-05, 1.2e-13, 3.3e-06
and 4.9e-13 at `n = 16, 24, 32, 48` — round-off on some meshes and 1e-5 on
others, **with no trend in the mesh**. It is not a stopping tolerance:
tightening rtol from 1e-8 to 1e-12 leaves the `n = 16` figures *bit identical*.
Both iterations are fully converged and their fixed points differ. The test now
sweeps three meshes and asserts the two things actually entailed — the **best**
agreement is at round-off, which is what says the Picard path solves MEQ's
problem rather than a neighbouring one, and the **worst** is bounded well below
a different problem. A per-mesh gate at 1e-6 is not reinstatable. See
*Should `PicardThenNewton` simply be the default?* for the same phenomenon at
9.4% across three solve routes.

**`NewtonConvergence`'s finite-difference Jacobian check reads 4e-11**, the
`O(step²)` floor, and MEQ keeps its hoisted `GetGradient()` regardless of
whether the library requires it: holding the linearisation fixed across a
difference is the right thing to write, and upstream's own finding is that
re-taking the gradient *after* a difference silently buys extra corrections and
can hide exactly this class of defect.

## Toroidal flow → `CLAUDE_FLOW.md`

`meq::RotatingSource` solves RoPP (136) closed by its (96) and (97), reachable
from a TOML file as `[source] Type = "rotating"`. It is a campaign of its own and
lives in `CLAUDE_FLOW.md`: two species in closed form, `n` species by a
safeguarded root find, normalised flux through the bordered Newton above, the
third derivative level it cost `meq::Profile`, the three errors found in Li & Zhu
— all of which converge beautifully — and Maschke & Perrin as the one exact
solution MEQ has with `T′ ≠ 0` and `ω′ ≠ 0`.

## The linear solves, and what they should be

A hybridized HDG scheme needs exactly two linear solvers: one for the global
face-coupled trace system, one for the small dense per-cell systems. MEQ has a
**third**, and it is not an oversight — it is the price of Newton.

| | what MEQ uses | |
|---|---|---|
| global trace | **selectable**: `UMFPackSolver` (METIS ordering, the default), `PardisoSolver` (`REAL_STRUCTURE_SYMMETRIC`), `CuDSSSolver` (`NONSYMMETRIC` + `FULL`) | three unsymmetric sparse LUs, agreeing to 5e-14; `setTraceSolver()` |
| *fallback, no direct solver at all* | GMRES | **unpreconditioned** on the Newton path, `GSSmoother` on the linear path |
| per-cell dense | MFEM `LUFactors`, partial-pivot LU | as a 2×2 block: LU on the flux block `A`, local Schur `S`, LU on `S` |
| **per-cell nonlinear** | element-local `NewtonSolver`, 100 iters, rtol 1e-12, `LPrecType::LU` | one iteration *per element per residual evaluation* — **and MEQ's default ordering no longer reaches it**, see below |

And a **fourth** when `ψ_ax` is an unknown, which is not a fourth solver: the
bordered system is solved by block elimination against the same factorisation of
the trace Jacobian, so it costs one extra backsolve per Newton step and nothing
else. Assembling the border into an `(n+1)` matrix would put a dense row and a
dense column into the factorisation for no gain, which is why it is not done.

**The third one exists because MEQ uses Newton rather than Picard, and
`NonlinearOrdering::NPC` removes it.** Picard evaluates `F` at the previous
iterate and leaves every local problem linear — one dense factorisation and done.
Condense-first Newton puts `∂F/∂ψ` inside the local problem and makes it
nonlinear. NPC, which is **MEQ's default**, differentiates the full
`(q, ψ, ψ̂)` system first and hybridizes the linear system that results, so every
element-local operation is a linear solve and `GetNumLocalNLIterations()` stays
at **zero** — asserted, not assumed, by
`SolverContract::theOrderingsAgreeAndOnlyOneIteratesLocally`, which reads 0
against the condensation's 3644/3560/3412 at `k = 1, 2, 3`.

So the row above is what MEQ *has*, and it is one `setNonlinearOrdering()` call
away rather than unreachable — `CondenseThenLinearise` is the **backup**, and on
stiff under-resolved meshes it is the one that works. See *The NPC port*.

**And there is a fifth thing that is not a solver at all but decides how fast
the third and fourth run**: `setAssemblyMode()`, which threads the element loop
that builds all of them. Measured under *Threading, measured*.

**The trace solver is a run-time choice**: `setTraceSolver()` picks among
whichever of the three the build has, and
`traceSolverAvailable()` says which those are without throwing. All four
construction sites — Picard, Newton, the bordered Newton and the linear path —
now go through one factory, so the *only* thing that differs between them is
whether the symbolic analysis is retained. It is, everywhere except the linear
path, which factorises once and destroys the object.

**Choosing one is never a numerical decision.** `theTraceSolversAgree` drives
each of them through the whole solver and pins the recovered `ψ` against the
first — 5.0e-14 for PARDISO against UMFPack — and cuDSS is checked the same way
by the `TraceSolverCuDSS` ctest, which is separate only because cuDSS needs an
`mfem::Device` and that is global state a Boost test case should not be
configuring.

### The trace matrix is symmetric on the fitted path and is not on the extension path

**And that is not a defect — it is what the extension technique is.** This is the
single most important thing in this section, and it was nearly missed by
measuring only the easy configuration.

`GradShafranov.cpp` justifies GMRES over CG with "the hybridized trace system is
small but not symmetric positive definite in this sign convention". Measured on
the **fitted** benchmark that is true as written and misleading:

→ **[M-54](MEASUREMENTS.md#m-54)** — `k` · `n` · nnz/row · rel `\|A_ij − A_ji\|` · Rayleigh quotients, free subspace

**Symmetric to round-off and negative definite** — so on a fitted mesh `−A` is
SPD and every symmetric method applies. The only positive diagonals are exactly
the essential trace dofs (64/96/128 of them, all exactly `1.0`), which is
`DIAG_ONE` putting a unit row in.

The **Newton Jacobian is the same** on that path: free–free block symmetric to
5e-16, every Rayleigh quotient in `[−1.91, −1.30]`, for `∂F/∂ψ` of either sign.
Its whole-matrix asymmetry reads 0.52, but that is entirely `GetGradient`
zeroing the essential rows without eliminating the matching columns — a
boundary-condition artefact, not a property of the operator. **Measure the
free–free block, not the assembled matrix**, or this looks like a wildly
unsymmetric problem.

**On the extension path none of that holds.** Same measurement, same code, a
curved `Γ`:

→ **[M-55](MEASUREMENTS.md#m-55)** — fitted · extension

The cause is `HDGExtensionIntegrator`, and it is structural rather than a bug.
Its element matrix is

```
elmat( dof*di + i, dof*dj + j ) += w * nor(di) * shape(i) * L(j, dj)
```

— an outer product of the normal trace of basis `i` against the **path lifting**
of basis `j`, two unrelated vectors. Measured on `Γ_h` faces its relative
asymmetry is exactly **1.0**, and it is **16.8× larger** than the `(R q, v)`
domain term it sits beside, which is itself symmetric to the last bit
(`max|A_ij − A_ji| = 0`). `DarcyForm::AssembleFluxMassBdrFaces()` deposits it
straight into the hybridization's per-element flux block, so it reaches both the
local factorisation and, through the Schur complement, the global trace matrix.

**So the transfer technique costs self-adjointness.** The continuous
Grad–Shafranov operator is self-adjoint; `Δ*` with a transferred Dirichlet datum
is not, because the datum on a face depends on the flux along a path leaving the
element. That is a property of Cockburn–Solano transfer, not of this
implementation, and it means **MEQ's headline configuration is genuinely
non-symmetric**. UMFPACK's unsymmetric LU is the right solver for it.

What survives on both paths is that the **symmetric part is negative definite** —
every Rayleigh quotient measured is negative on either path. So `−A` has positive
definite symmetric part, which is what gives GMRES a convergence bound and what
makes a preconditioned Krylov method reasonable at all.

### A quarter of every Newton step was thrown away — fixed, and switched on

`NewtonSolver::Mult` calls `prec->SetOperator( *grad )` every iteration while
the sparsity pattern does not change between Newton steps, so a
`UMFPackSolver::SetOperator` holding `void *Symbolic` as a **local variable**
recomputes and discards the symbolic analysis each time — and MEQ asks for METIS
ordering, which makes it dearer than the default.

**`Symbolic` is now a member and MEQ calls `SetReuseSymbolic()`.** Verified by
count rather than by clock — `theSymbolicAnalysisIsReusedAcrossNewtonSteps`
prints and asserts it:

```
  Newton took 4 iterations: 1 symbolic analyses, 4 numeric factorisations
```

One analysis, one factorisation per step. **A count, not a timing, on purpose**:
a timing here would be a measurement about this machine, and this is a
measurement about the code. It is also the *only* thing that could notice the
reuse lapsing — the pattern check is exact and re-analyses whenever it fails, so
a lapse costs speed and nothing else. No wrong answer, no failed convergence,
nothing a rate table or an error norm could ever see.

**Where it is switched on, and where deliberately not:**

| site | reuse | why |
|---|---|---|
| Newton path | **yes** | `NewtonSolver::Mult` re-`SetOperator`s the same object every iteration |
| Picard path | **yes**, and the solver is **hoisted to a member** | it was built inside `picardStep()`, so it had nothing to reuse across calls; Picard runs 122 to 290 full factorisations, a larger absolute win than Newton's |
| linear path | no | one factorisation, object destroyed immediately: retaining the analysis buys nothing and costs a copy of the pattern |

The Picard hoist is safe across `prepare()` rebuilding `reduced` every iteration,
because `SetReuseSymbolic()` documents comparing the **pattern**, entry by entry,
rather than the object — it "accepts a matrix reassembled in place and a matrix
rebuilt into a fresh object with the same structure alike", and re-analyses
whenever the check fails. Reuse is a request, never an assumption.

The measurement that motivated it:

→ **[M-56](MEASUREMENTS.md#m-56)** — `n` · symbolic · numeric · backsolve · symbolic share

### Threading, measured — and the two axes are not the same axis

**`tests/performance/` is the harness and is deliberately NOT a ctest**:
everything it reports is a timing, and this project's standing rule is that a
threaded timing on this machine is a measurement about the machine.
`TraceSolverScaling` returns non-zero only for the two *correctness* properties
that make its timings mean anything — threaded assembly reproducing serial
assembly bit for bit, and PARDISO reproducing UMFPACK. `scan.sh` drives it one
process per point, because MKL fixes its threading at first use and an
in-process sweep would measure the first setting five times.

**Two knobs, and sweeping them together is actively misleading:**

| | drives |
|---|---|
| `OMP_NUM_THREADS` | `DarcyHybridization`'s threaded element assembly |
| `MKL_NUM_THREADS` | UMFPACK's BLAS, PARDISO's internals, **and MFEM's element-local dense LU** |

That third entry is the surprise, and it is why the axes must be separated: swept
together, the `k = 3` assembly blow-up swamps every other column. Full account
under `CLAUDE.md`, *Traps*; the short of it is `MKL_NUM_THREADS=1`, non-negotiable, and the
culprit is `ComputeH()`'s dense LU rather than UMFPACK.

**THREADED ASSEMBLY IS NOW EVERY ELEMENT-LOCAL LOOP, NOT JUST `ComputeH()`,
AND THAT INVERTS THE DEFAULT.** `SetAssemblyMode( Threaded )`
threads more than the assembly: `MultNL()` too — the residual and the Jacobian
assembly, and so `NPCResidual()` and `NPCGradient()`,
which is **every NPC step** — plus `ReduceRHS()`, `ComputeSolution()` and the
two RHS eliminations. The two kinds of loop are threaded differently and the
difference is the scatter *target*: `ComputeH()` writes an unfinalized
`SparseMatrix`, which carries one `current_row` and one RowNode allocator for
the whole matrix, so its scatter stays serial and in element order; `MultNL()`
writes a `Vector`, where disjoint indices really are independent, so it is
walked in element **colour** order. Both are bit for bit — a colouring only
changes the order in which a face's two elements accumulate, and `a + b = b + a`.
MFEM notes this would *not* survive an `H1_Trace` (EDG) trace space, where a dof
takes more than two contributions; MEQ's is `DG_Interface_FECollection`, so it
holds.

**EXCEPT FOR TWO, AND ON THE BORDERED PATH THEY ARE THE LARGEST LEG OF A
THREADED STEP.** `NPCReduce()` and `NPCRecover()` — the element loops either
side of the trace solve, which is how `J^-1` is applied — carry no
`omp parallel`, single-vector or blocked. That is deliberate upstream and the
doxygen on `NPCGradient()` gives the reason: they are *"under 6% of the step,
flat in mesh size and order, and Amdahl caps any gain there"*, on the pedestal
problem at four resolutions.

**THAT MEASUREMENT IS FIXED BOUNDARY AND ONE RIGHT-HAND SIDE, AND MEQ'S
BORDERED STEP IS NEITHER.** `J^-1` is applied to the field residual **and to
every border column** — `psi_ax`, `psi_bnd`, the prescribed current, the
X-point's two coordinates and the `N` exterior modes, all against the one
factorisation — through `DarcyNPCSolver::ArrayMult()`. So the traversal is
`O( elements x columns )` where the integrator-bound loops are
`O( elements )` and already threaded, and the share inverts with the border and
again with every exterior mode. Measured on DIII-D at 10 modes, **0.212 s at one
thread and 0.212 s at eight**, which is 30.6% of the threaded step and the
largest leg in it.

`StepProfile` has a leg for it, and that leg is where `otherSeconds()` went:
`TimedSolver` decorates the TRACE solver, so it times the middle of
`NPCReduce -> solve -> NPCRecover` and neither end, and the remainder was 45% of
a threaded step and named nothing. It is timed as the wall clock of
`ArrayMult()` **minus** the trace solve inside it, so it is disjoint from the
backsolve leg rather than containing it, and `other` falls to 14%.

→ **[M-126](MEASUREMENTS.md#m-126)** — the leg split at 1 and 8 threads · the
two instruments that found it before the leg existed · the barrier spin · what
threading it would be worth

**AND A THIRD OF THE `%CPU` ON THAT RUN IS BARRIER SPIN RATHER THAN WORK**,
which matters here because it is what a utilisation figure is read off.
`OMP_WAIT_POLICY=passive` takes the CPU from about 308% to about 188% for
roughly 1.5% of the wall clock; MEQ's own OpenMP regions and MKL's contribute
independently and about equally. **Set it before reading anything into a
`%CPU`**, and do not read a low one as idle cores — here it is the honest number
and the high one is the artefact.

**AND THE ONE INTEGRATOR MEQ EVALUATES PER ELEMENT UNDER THREADING IS ITS
OWN, WHICH IS WHY THIS IS SAFE.** MFEM's stock integrators are not uniformly
thread safe — `VectorMassIntegrator`, `MassIntegrator` and `DiffusionIntegrator`
hold their scratch as plain members, and 67 of the 99 classes carrying scratch
are unguarded. MEQ installs the first of those. It is never reached from a
threaded loop, and the reason is structural rather than lucky:

* `CopyLinearGradBlocks()` **always declines for MEQ**, on
  `lop_type == LocalOpType::PotNL`.
* On that branch `ConstructGrad()` skips the flux block outright —
  `else if ( !ad_done && !A_empty && lop_type != PotNL )` — because a
  specialised local operator keeps its block factored in place.
* What it does evaluate per element is `m_nlfi_p`, and that is
  **`meq::SourceIntegrator`**, MEQ's own, made reentrant when it first met this
  loop. See *Traps* in `CLAUDE.md`.

The bit-exactness assertions between assembly modes are the independent check
and have been green throughout; a shared-scratch integrator on eight threads
does not survive them.

**AND UPSTREAM IS ABOUT TO REFUSE THIS CONFIGURATION ON THE WRONG PREDICATE.**
A proposed abort keys on *did `CopyLinearGradBlocks()` decline* — which MEQ does
on every run, for the second of two unrelated reasons. Filed back as
`../mfem-hdg-dev/doc/HDG-THREADED-REFUSAL-FROM-MEQ.md`; **expect an MFEM update
to abort every threaded MEQ run until it lands.**

**AND WITHOUT IT `--profile`'s `cores` COLUMN INVERTS.** A serial leg reads
HIGH, not low, because the other seven threads are spinning through it: on one
run `constraint location` reads 1.01 cores under `passive` and **5.76** under
the default, and `I_p` reads 1.03 against **8.95**.

**MEQ'S DEFAULT IS THEREFORE `Threaded` NOW, AND WAS `Serial`.** The measurement
that settled it before was taken against the option as it then was, and the
inversion is the whole justification:

→ **[M-57](MEASUREMENTS.md#m-57)** — Serial · **Threaded**

**The bordered Newton is the case that argues hardest FOR the flag**, and the
reason is structural rather than lucky: it is dominated by residual evaluations
rather than by assembly, and residuals are what `MultNL()` threading covers. It
does not cost anything at one thread either, which is what makes this safe as a
default rather than merely faster on average.

`defaultAssemblyMode()` is **build-conditional**, and it has to be:
`buildForms()` passes the mode straight to MFEM, whose abort-rather-than-fall-back
is only guarded by `setAssemblyMode()`, and the constructor bypasses the setter.
A build without `MFEM_USE_OPENMP` or `MFEM_THREAD_SAFE` therefore defaults to
`Serial`. `Config`'s default is a plain `Threaded` because that header is
MFEM-free; `apps/meq.cpp` downgrades it when the build cannot honour it, and
**refuses** only when the file said `"threaded"` explicitly.

**THE BIT-EXACTNESS HOLDS AT `MKL_NUM_THREADS=1` AND NOT ABOVE IT**, which is
new and is not a defect. `threadedAssemblyReproducesSerialAssemblyExactly` and
its nonlinear sibling both require `0.000e+00`, and both get it under ctest,
which sets that variable. At `MKL_NUM_THREADS=8` the two modes differ by
**1.3e-15 in `ψ` and 1.3e-13 in the flux** — because the serial path hands MKL
eight threads and the threaded path hands it one, so a blocked BLAS-3 sums in a
different order from an unblocked loop. Arithmetic reassociation inside MKL,
not a race in MFEM or in MEQ. It is also independent confirmation of the nesting
result below: the only way the two modes can round differently is if MKL is
doing something different in each.

**AND THE NONLINEAR CASE HAD TO BE ADDED, BECAUSE THE EXISTING ONE STRUCTURALLY
COULD NOT SEE THE LOOP THAT MATTERS.**
`threadedAssemblyReproducesSerialAssemblyExactly` builds its source as an
`mfem::FunctionCoefficient`, which takes MEQ's **linear** path:
`usesNonlinearForms()` is false, `meq::SourceIntegrator` is never installed and
`MultNL()` is never called. So it exercises `ComputeH()` and nothing else — which
was the whole option when it was written.
`threadedAssemblyReproducesSerialAssemblyOnANonlinearSource` is the missing case
and it asserts three things: `ψ`, the flux, **and the Newton iteration count**.
The count is the one with teeth — a racing *residual* moves the answer and the
first two catch it, while a racing *gradient* reaches the same discrete solution
by a longer path, which is exactly the failure *A wrong Jacobian is invisible to
a convergence table* says no error norm and no rate can see.

**IT FOUND A REAL RACE, IN THE ONE INTEGRATOR MEQ INSTALLS ON THAT LOOP.**
MFEM states the obligation plainly and says it cannot check it: *"Any integrator
the caller installs … sits on this loop and must be thread-safe too. An
integrator holding per-point scratch as a plain member will race, silently."*
`meq::SourceIntegrator` held exactly that — `mfem::Vector shape`, resized and
refilled per quadrature point, shared across every element of a colour — and
`SourceIntegrator` **is** the whole semi-linear term. It is now guarded on
`MFEM_THREAD_SAFE`, the way MFEM guards its own (`HDGDiffusionIntegrator`'s
dozen members are guarded the same way, and MEQ's build has
`MFEM_THREAD_SAFE = YES`, so they do not exist here).
`meq::PoloidalFieldCoefficient` carried the same defect, latent because nothing
constructs it, and is guarded rather than deleted.

**The test was checked to DISCRIMINATE, not merely to pass.** With the guard
deliberately removed and the race reinstated, it fails 3 runs out of 3 — and it
fails *loudly*, `NewtonSolver` aborting on a NaN residual rather than returning
a quietly wrong answer, because a concurrent `Vector::SetSize` is a reallocation
and so is memory corruption rather than merely a stale read. **Do not read that
as a guarantee the failure would always be loud**: the same race on a different
mesh or thread count can just as easily return finite nonsense.

**One obligation is recorded and NOT closed**: `meq::RotatingSource` throws from
`f()` and `dFdPsi()` when a species temperature goes non-positive or the
quasineutrality root find fails, and an exception escaping an OpenMP structured
block is undefined behaviour. It is left recorded because the loop is MFEM's —
MFEM has the same exposure through its own `MFEM_VERIFY` under
`MFEM_USE_EXCEPTIONS` — and inventing an error path the library does not support
would replace a crash with a silent NaN. Prefer `Serial` for a rotating source
until the iterate is known good.

**PARDISO beats UMFPACK on the trace solve, and beats it even sequentially** —
1.50x on analyse-plus-factor and 1.41x on the backsolve at 37,248 trace dofs
with `MKL_NUM_THREADS=1` on both sides, agreeing to 1.0e-14 or better at every
point. It scales to about 8 threads (1.87x setup, 1.96x solve; 16 buys nothing
more).

**AND THOSE THREADS ARE NOW SPENDABLE, WHICH THIS FILE SAID THEY WERE NOT.
THREADED ASSEMBLY IS WHAT UNLOCKS THEM, AND IT NEEDED NO NEW CODE.** The
paragraph that stood here said `MEQ cannot have that`, because `MKL_NUM_THREADS`
is process-wide and the setting that makes PARDISO fast is the setting that
makes `ComputeH()` forty times slower at `k = 3`; and that item 0 of *What to
do* — `mkl_set_num_threads_local()` around the trace solve, or
`LocalFactorMode::Batched` — was the way out. **Neither is needed.**

**MKL SUPPRESSES ITS OWN THREADING INSIDE AN ACTIVE OpenMP REGION.** Once the
element loop is itself such a region, the element-local dense work is *nested*,
MKL runs it sequentially, and `MKL_NUM_THREADS` costs it nothing — while the
trace solve, which runs on the master thread **outside** any parallel region,
still takes all of them. Measured on a whole nonlinear solve, `k = 3, n = 16`,
`OMP_NUM_THREADS=8`:

→ **[M-58](MEASUREMENTS.md#m-58)** — `MKL=1` · `MKL=8`

**1285x between the two modes at `MKL=8`.** The recipe is therefore
`AssemblyMode::Threaded` + `TraceSolver::Pardiso` + `OMP_NUM_THREADS` and
`MKL_NUM_THREADS` set to the **same** value above one. What PARDISO's threads
themselves add on top is modest — 0.1313 s to 0.1163 s at `k = 2, n = 24`, about
**1.13x** — so the threading win is overwhelmingly the assembly and PARDISO is
the part that makes it *safe* to ask for MKL threads at all.

**TWO CONFIGURATIONS TO REFUSE OR WARN ABOUT, AND BOTH WERE MET WHILE
MEASURING.**

* **UMFPACK can never take MKL threads, whatever the assembly mode.** Its BLAS
  calls are in the trace solve, on the master thread, outside any parallel
  region — so they are *not* nested and get the full count. Measured in situ: a
  `NpcThreadScaling` sweep at `MKL=8` that included UMFPack rows sat at **266%
  CPU making no progress**, which is precisely the barrier-spinning signature
  `CLAUDE.md`'s *Traps* records for `SolovievConvergence` at 1.24 s against
  177.83 s. If
  `MKL_NUM_THREADS > 1`, use PARDISO.
* **Threaded assembly at `OMP_NUM_THREADS=1` with MKL threads on is
  catastrophic**, and it is a configuration a user can reach by accident. A team
  of one thread does not get the nested-region suppression. Measured on the
  isolated element-local kernels: **12.4 s against 0.069 s serial** at `k = 3`.
  `apps/meq.cpp` warns about exactly this combination at startup, and only when
  it is actually present.

**AND THE `ComputeH()` STORY THIS FILE TOLD WAS WRONG IN ONE DETAIL WORTH
CORRECTING.** Item 0 names "`ComputeH()`'s element-local dense **LU**". Measured
on the kernels at MEQ's own block sizes, `dgetrf` **does not degrade at all** —
0.0137 s at `MKL=1` against 0.0138 s at `MKL=8` at `k = 3`, and 0.0256 against
0.0262 at `k = 4`. MKL does not thread a `dgetrf` that small. What degrades is
the **back-substitution** (`dgetrs` with `nrhs` = the element's trace dofs) and
the **Schur-complement `dgemm`**: at `k = 3`, 0.0070 to 0.0331 and 0.0020 to
0.0294. The `k = 2` / `k = 3` threshold `CLAUDE.md`'s *Traps* records is
confirmed — `k = 2`
does not move at any thread count — but the **40x magnitude is not reproduced**
on the clean single-MKL link line, which post-dates that measurement; the
kernels give 3.3x. The 383x above is in situ and on the whole nonlinear solve,
where `MultNL()`'s element-local work is hit on every residual and every
Jacobian rather than once per assembly.

**cuDSS is correct here and not measurable here.** It agrees with UMFPACK to
**3.5e-14 or better** from 9,408 to 148,224 trace dofs, which is the question
that needed answering locally. The timings are irreproducible **by a factor of
thirty at fixed size** — three runs of the same binary on the same 9,408-dof
problem gave 1.05 s, 1.13 s and 2.81 s — because WSL2 shares the GPU with the
Windows host, whose load MEQ does not control, and because a consumer card runs
FP64 at 1/32 of FP32 where a datacentre part runs it at about 1/2. It won nothing
on this hardware. **Do not re-time it here**; it needs a different machine, not
another afternoon.

**THE MEMORY FIGURE THIS PARAGRAPH USED TO CARRY IS DELETED RATHER THAN
UPDATED.** It read *"`nvidia-smi` reports 7.5 of 8 GB used with no compute
processes"*, and today the same command reports **1.1 of 8 GB** — the host's
load, not a property of the card or of MEQ. The FP64 ratio is the part of the
argument that does not move, and it is the part that decides.

**AND ONE TRAP THAT WOULD HAVE POISONED ALL OF IT.** cuDSS queues work on a
stream and returns. Timed without a device synchronise, the warm setup of a
**148,224-unknown factorisation reads 2.0e-04 s** — four orders of magnitude
out, and *plausible enough to publish*, because "the reordering was reused, so
of course it is fast" is a story that fits. `TraceSolverScaling` calls
`cudaDeviceSynchronize()` inside the timing loop for **every** solver, CPU ones
included, so the sync can never be the thing that was forgotten. **Any future
device timing in this project must do the same.**

**What a fresh `GradShafranovSolver` does:** `setTraceSolver()` is `UMFPack`
(the only backend present in every build, and what every rate in the suite was
measured with); `setAssemblyMode()` is **`Threaded` wherever the build can
honour it** and `Serial` otherwise — build-conditional, see
above for why the default moved; `MKL_NUM_THREADS` is 1, set on every ctest.

**The fastest reachable set is PARDISO plus threaded assembly at 8**, and the
1.24x figure is what a **linear** solve gives, where the flag threads assembly
and nothing else. On a **nonlinear** solve — which is what MEQ does — the same
set is worth **2.8x to 3.0x**, because the flag now
threads the residual and the Jacobian too. PARDISO's own contribution within
that is about **1.13x**; the rest is the assembly mode.

Threaded assembly **is** the default now. PARDISO is not, because
`MFEM_USE_MKL_PARDISO` is off in most builds and oneMKL's terms are not
everybody's to accept; `setTraceSolver()`, `traceSolverAvailable()` and
`[solver] TraceSolver` are how a caller who has it takes it. **The reason to
take it is not its 1.13x** — it is that PARDISO is the only trace solver that
tolerates `MKL_NUM_THREADS > 1` at all, and UMFPACK collapses there whatever the
assembly mode.

### What to do, in order of value

**The element-local dense work does not need getting off threaded MKL, which is
which is the obvious place to start.** `AssemblyMode::Threaded` settles it: MKL
suppresses its own threading inside an active OpenMP region, so once the element
loop *is* such a region the local work is nested and free — **1285× between the
two modes at `MKL=8`, `k = 3`**. Two details worth having, because both are easy
to get wrong: it is not the dense **LU** that degrades — `dgetrf` does not move
at these block sizes at all — but the back-substitutions and the
Schur-complement `dgemm`; and **`LocalFactorMode::Batched` is not the answer and
should not be taken**, since it batches `InvertA()`/`InvertD()`, which run once
from `Finalize()`, while the per-linearisation factorisation is inside
`ComputeElementH()`, and it aborts on an exact zero pivot where the serial loop
carries on.

**Symbolic reuse and a run-time trace solver are both in.**
`theSymbolicAnalysisIsReusedAcrossNewtonSteps` asserts the reuse by count —
about 23% off each step for no numerical change, on both paths — and
`setTraceSolver()` picks among three backends that agree to 5e-14 through the
whole solver. **The default stays UMFPack**, because oneMKL's licence is not
everybody's to accept; a caller who has PARDISO gets it for one line, and its
*sequential* advantage is 1.50× on the factorisation with no MKL threads at all.

**And the ranking below is the one the ASYMMETRY finding gives.** Ranking
Cholesky highly is what measurements taken only on the fitted path suggest, and
it is wrong wherever `Γ` is curved.

1. **Precondition the Newton path's GMRES fallback.** It has none, while the
   linear path's fallback has a `GSSmoother`. That inconsistency is a plain bug
   and is worth fixing whatever else happens; without SuiteSparse the Newton
   path is currently running unpreconditioned GMRES on every step.
2. **Symmetric methods — fitted path only.** CG or MINRES on `−A`, and a CHOLMOD
   Cholesky, are all correct on a fitted mesh and all invalid on the extension
   path. CHOLMOD is already linked (`-lcholmod` via SuiteSparse) but MFEM wraps
   only UMFPACK and KLU. Worth having only if fitted-domain runs become a
   workload in their own right; the driver's target is curved boundaries, where
   none of it applies.
3. **Cholesky on the per-cell flux block — do not.** `(R q, v)` is SPD, so this
   looks like a free 2×. It is wrong on the extension path for exactly the
   reason above: `AssembleFluxMassBdrFaces()` puts a term with relative
   asymmetry 1.0, sixteen times larger than the mass term, into the very block
   `LU_A` factors. A correct version would have to branch on whether an
   extension is configured, for at best a small gain on 9×9 to 30×30 blocks.
   Not worth a conditional in that inner loop. **`LUFactors` is the right
   choice and should stay.**
4. **Do not reach for AMG.** At 2D serial sizes direct wins: 50k trace dofs
   factorise in 0.4 s. AMG earns its place in 3D or in parallel, and serial MFEM
   has none anyway — `HypreBoomerAMG` needs MPI.
5. **Do not re-time cuDSS here.** Correctness is established, the timings are
   irreproducible on this machine by a factor of thirty, and the reason is the
   machine. It needs a datacentre part, not another afternoon.

### PARDISO and the MKL link line: what is still true

*Threading, measured* carries the numbers; three things belong here.

**oneAPI MKL 2026.1 fixed a real defect, and it was packaging.** Debian's
`intel-mkl` 2020.4.304 returned error `-3` from PARDISO at 12,544 and 28,032
trace dofs — the sizes MEQ actually runs. Against oneAPI at
`/opt/intel/oneapi/mkl/latest` it runs at every size and agrees with UMFPACK to
round-off. **Both paths must keep shipping**: oneMKL's licence is not
everybody's to accept and `MFEM_USE_MKL_PARDISO` is off in most builds, so
`tests/convergence/TraceSolverComparison.cpp` compiles and passes either way,
skipping the PARDISO columns when they are absent. It asserts only the
**agreement**; the timings are printed, per the standing rule.

**`SetReuseSymbolic()` exists on `PardisoSolver` too**, with the same
pattern-comparison contract as `UMFPackSolver`'s, so *A quarter of every Newton
step* carries over without re-arguing.

**MEQ BUILDS ITS OWN SuiteSparse SO THAT EXACTLY ONE MKL IS LOADED, AND THE LINK
LINE IS ONLY HALF OF WHY.** Take Debian's and the line carries oneAPI's
`mkl_gnu_thread` from the explicit `BLAS_LIBRARIES` beside Debian's
`mkl_intel_thread` plus `iomp5` behind SuiteSparse — two MKL *versions*, two MKL
threading layers and two OpenMP runtimes, working only because oneAPI comes
first in the link order.

**THE HALF THAT NO CMAKE VARIABLE COULD FIX**: Debian's `libumfpack.so` and
`libcholmod.so` carry a hard `NEEDED` on `libblas.so.3`, which on this machine
is a Debian alternatives symlink to **`libmkl_rt.so`**. So Debian's MKL 2020
loaded at runtime whatever the link line said. **Editing the link line alone
would have looked like a fix and left the real one in place.**

**The fix is MEQ's own SuiteSparse**, at `../suitesparse`, exactly as
`../sundials/cuda-install` is MEQ's own SUNDIALS. It is built from **v7.12.2 —
the same version Debian ships** — so the only variable that changes is the BLAS,
not the numerics:

```sh
MKL="-L/opt/intel/oneapi/mkl/latest/lib;-Wl,-rpath,/opt/intel/oneapi/mkl/latest/lib;\
-lmkl_intel_lp64;-lmkl_gnu_thread;-lmkl_core;-lgomp;-lpthread;-lm;-ldl"
cmake -S src -B build -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX=/home/ian/projects/suitesparse/install \
  -DSUITESPARSE_ENABLE_PROJECTS="suitesparse_config;amd;btf;camd;ccolamd;colamd;cholmod;klu;umfpack" \
  -DBLAS_LIBRARIES="$MKL" -DLAPACK_LIBRARIES="$MKL" -DBLA_VENDOR=Intel10_64lp \
  -DBLAS_INCLUDE_DIRS=/opt/intel/oneapi/mkl/latest/include \
  -DCMAKE_C_FLAGS=-I/opt/intel/oneapi/mkl/latest/include \
  -DBUILD_STATIC_LIBS=OFF -DSUITESPARSE_USE_CUDA=OFF -DSUITESPARSE_USE_FORTRAN=OFF
```

**Three snags, none of them obvious, all cost a configure round trip.**

* SuiteSparse takes `BLAS_LIBRARIES` as-is if you set it — but its threading
  probe then does `string(REGEX MATCH "^Intel" ... ${BLA_VENDOR})` on a
  `BLA_VENDOR` that is **empty**, because supplying the libraries skips
  `find_package(BLAS)`. Set `BLA_VENDOR` anyway; it is informational.
* That probe's `try_run` passes `LINK_LIBRARIES` but **never
  `BLAS_INCLUDE_DIRS`**, so it fails on a missing `mkl.h`. The include path has
  to go in `CMAKE_C_FLAGS` as well.
* Only the nine projects MFEM needs are built. `SUITESPARSE_ENABLE_PROJECTS`
  keeps GraphBLAS — which dwarfs everything else — out of the build entirely.

**AND ON THE MFEM SIDE THE STALE CACHE HAD TO GO.** `mfem_find_package` resolves
SuiteSparse's `BLAS` requirement through **CMake's own `FindBLAS`**, which had
cached Debian's paths in `BLAS_mkl_*_LIBRARY`. Deleting those
(`cmake -U "SuiteSparse_*" -U "BLAS_mkl_*" -DSuiteSparse_DIR=...`) is what lets
MFEM's explicit oneAPI `BLAS_LIBRARIES` win.

**Why `FindBLAS` chose the WRONG threading layer, which is worth knowing because
it will do it again**: it selects `mkl_gnu_thread` + `gomp` only when a
**Fortran** compiler is loaded and is GNU; otherwise it takes
`mkl_intel_thread` + `iomp5`. SuiteSparse and MFEM are C/C++ builds, so the
default is Intel threading — mismatched with everything else here. That is
`FindBLAS.cmake:505-513`, and it is why the libraries are handed over
explicitly rather than discovered.

**VERIFIED THREE WAYS, because "one MKL" is a claim about the loader and not
about a text file.** The installed `config.mk` carries oneAPI and nothing else;
`ldd` on a test binary resolves `libmkl_intel_lp64`, `libmkl_gnu_thread` and
`libmkl_core` to `/opt/intel/oneapi/...` with **no `libmkl_rt`, no `libblas.so.3`,
no `libiomp5`, and `libgomp` as the only OpenMP runtime**; and `LD_DEBUG=libs`
shows the only MKL objects *initialised* are oneAPI's. Debian's `libmetis` is
still on the line and is **not even loaded** — CHOLMOD now bundles its own,
prefixed `SuiteSparse_*`/`cholmod_*`, so there is no symbol collision and
UMFPACK keeps its METIS ordering. Suite green twice — **23/23 as it then was**, 490 s and 540 s on
the same code, which is the run-to-run spread on this machine and is worth
knowing before anyone reads a 10% change in a suite time as a result.

**`MKL_THREADING_LAYER=GNU` is inert and has been dropped everywhere.** It only
ever configured the `libmkl_rt` **dispatcher**, and there is no dispatcher any
more: `SolovievConvergence`, `NewtonConvergence` and `FieldConvergence` produce
**bit-identical output** with it set and with it unset. It was briefly kept as a
free guard, and that argument was rejected for a reason that generalises —
**most machines will not have MKL at all**, so it named a library the majority
of readers do not have, about a failure they cannot suffer. A guard that is
usually inapplicable is a false instruction to whoever reads it next, and it
spends the credibility of the ones that are real. The driver's matching runtime
warning went for the same reason. **Choosing a threading layer is CMake's job**,
which can look at what the build actually links.

**The hazard went off once, and in the opposite direction to the one predicted.**
The prediction was wrong answers from mixed threading layers. What actually
happens is that a stray `libmkl_sequential` is **load bearing**: it makes MKL
resolve sequential for everything, so timings taken against it are fast and
nobody notices that threaded MKL is ruinous here. Removing it — the *correct*
thing to do — exposes a 140x regression.
**The suite was never deliberately sequential; it was accidentally so.** Same
species of finding as the threaded-BLAS one: a property of the link line
masquerading as a property of the code.

**And the honest caveat on all of it**: on a hard case the dominant cost is not
the global solve. Globalisation is a bigger lever than anything in this
section, and under `CondenseThenLinearise` the element-local *non-linear*
iteration dominates outright.

**The transferable lesson**, which is the same one the Solov'ev coefficients
taught: a property measured on the easy configuration is not a property of the
code. Symmetry held to 2e-16 on a fitted rectangle and failed at 5.4e-1 on the
geometry MEQ is actually for.
