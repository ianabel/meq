# CCSZ interpolatory HDG for MEQ's source term

Written 2026-09-15. `CLAUDE.md` is the maintainer's index and
`CLAUDE_HDGGS.md` owns the discretisation; this file is the design for one
change to it and defers to both on anything they already say.
`refs/Refs.md`, *Interpolatory HDG*, is the bibliography and both papers are
on disk.

**Line numbers in `src/meq/GradShafranov.cpp` are as read on 2026-09-15 and
that file is under active edit** — it moved by about twenty lines in the
`buildForms()` region during the afternoon this was written. Every citation
below also names the function, which is the durable half; `grep` for the name
rather than trusting the number.

**The change in one sentence.** Stop integrating `F( r, z, ψ )` at quadrature
points inside `meq::SourceIntegrator` and instead **interpolate it into the
postprocessing space `Z_h = P^{k+1}`**, so that the source's contribution to
the residual and to the Jacobian become fixed matrices applied to a vector of
nodal `F` and `dF/dψ` values.

**THE SINGLE BIGGEST FACT, AND IT CHANGES THE SHAPE OF THE WHOLE JOB: MFEM HAS
ALREADY BUILT THIS, AND MEQ LINKS IT TODAY.**
`fem/darcy/reaction_hdg.hpp` and `fem/darcy/postprocess_hdg.hpp` in
`../mfem/install/include/mfem/` carry `HDGPostprocessBlocks`,
`NodalReactionFunction`, `HDGReactionIntegratorBase`,
`HDGInterpolatoryReactionIntegrator` and `HDGQuadratureReactionIntegrator`, and
`nm` finds every one of them defined in `libmfem.a`. `DarcyHybridization` carries
`Bg_data`, the Jacobian's `(1,0)` block, written specifically for this method and
named after it in its own doxygen. So this is **not** an MFEM feature request
with an MEQ consumer attached; it is an MEQ wiring job against a built,
documented and unit-tested library facility, plus two small things upstream does
not yet offer (§7).

---

## 0. WHAT WOULD FALSIFY THIS, AND THE CHEAPEST EXPERIMENT

Put first because it decides whether the rest is worth reading.

### 0.1 The falsifier for the PERFORMANCE case, and it costs about two hours

**`meq::SourceIntegrator` already has the knob that bounds the prize.**
`GradShafranovSolver::setSourceQuadratureOrder( int )`
(`src/meq/GradShafranov.cpp:1469`) sets `sourceQuadratureExtra`, default **4**,
and the rule is `2k + tr.OrderW() + extraOrder`
(`src/meq/GradShafranov.cpp:353-364`). On a straight-sided triangle
`IsoparametricTransformation::OrderW()` is `(1-1)·2 = 0`, so the shipped rule at
`k = 2` is order 8, which `mfem::IntRules` serves with **16 points**; at
`extraOrder = 0` it is order 4, served with **6 points**.

**So MEQ can delete 10 of its 16 source quadrature points today, with one
existing setter and no new discretisation.** Run `tests/performance/NewtonStepProfile`
— which already splits a step into residual / gradient / factor / backsolve /
plasma-fill / other, and already reports `ComputeH()` as a sub-split of the
gradient — over `extraOrder ∈ { 4, 2, 0 }` on all three of its cases
(`example5`, `pedestal`, `highbeta`), medians of five, `MKL_NUM_THREADS=1`, on
an idle machine.

**Read it as a derivative.** Let `T( N_q )` be the step time. Going 16 → 6
points removes 62.5% of the source quadrature work. The interpolatory method
removes what is left of it *and then charges* 10 nodal `F` evaluations plus two
GEMMs (§3.3), so

```
prize( interpolatory )  <  ( 16 / 10 ) × [ T(16) − T(6) ]
```

is a hard upper bound in the same units, measured on MEQ's own cases with MEQ's
own profiler. **If `T(16) − T(6)` is inside the 1.4%–14% run-to-run scatter
M-80 records, the performance case is dead and the port must stand on §0.3
alone.**

It needs one flag added to `tests/performance/NewtonStepProfile.cpp` (a
`--source-extra N` passed through to `setSourceQuadratureOrder()`), a `-j6`
rebuild and three runs. **No number for it is printed anywhere in this
document, because it has not been run.**

**A second reading comes free from the same runs**: the rates in
`SolovievConvergence` and `NewtonConvergence` at `extraOrder = 0` say how much
of the shipped over-integration is buying accuracy rather than time.

### 0.2 The falsifier for the ACCURACY case, and MEQ already owns the instrument

`tests/convergence/PlasmaEdgeConvergence.cpp` with `tests/analytic/PlasmaEdge.hpp`
produces M-10 — `ψ*`'s rate over `j ∈ { 0, 1, 2, 3 }` and `k ∈ { 1, 2, 3 }`,
with `k+2` reached exactly on the entries where `k ≤ j`. Re-run that table with
the interpolatory term in place of the quadrature one. **If the bold entries
stop reaching `k+2`, the interpolatory route is refused on the free-boundary
path and kept on the fixed-boundary one** — which costs nothing, because the two
integrators coexist by construction (§4.1). §5.2 gives the arithmetic that says
this is the likely outcome and why it is an estimate rather than a result.

### 0.3 The argument that survives both falsifiers

**The interpolatory form is the only shape in which MEQ's own integrator can
ever reach a device.** `CLAUDE.md`, *Which integrators need device offload*,
ranks `meq::SourceIntegrator` **first** and says it "is the only one on this
list MEQ can write itself". A quadrature loop containing a virtual
`meq::Source::f()` call per point has no kernel in it. `A9 F( γ )` and
`A9 diag( F'( γ ) ) [ B11 B12 ]` are a batched GEMM and a pointwise map, and
`DarcyHybridization::LocalFactorMode::Batched` already carries the `(1,0)` block
through the batched routes — upstream's
*"The batched local routes carry the (1,0) gradient block"* is a test for
exactly that. This argument is structural and is not decided by §0.1.

**And it deletes a load-bearing constraint.** `GradShafranov.cpp:3396-3418`
documents at length why *the whole potential block* has to sit on the non-linear
form — "the one structural decision in this file that is not obvious from the
weak form". With the source on the **block** non-linear form instead, the HDG
face stabilisation goes back on the **linear** `M_p`, which is precisely what
MFEM's own end-to-end test does
(`tests/unit/fem/test_darcy_reaction.cpp`, `RunAffine()`). That removes the
silent-drop hazard the comment warns about rather than documenting it better.

---

## 1. What CCSZ interpolatory HDG is, stated for MEQ's equation

### 1.1 The papers, and which one

`refs/SuperconvergentHDG-I.pdf` — Chen, Cockburn, Singler & Zhang,
*Superconvergent Interpolatory HDG Methods for Reaction Diffusion Equations I:
An HDG_k Method*, J. Sci. Comput. **81** (2019) 2188–2212,
doi `10.1007/s10915-019-01081-3`. **This is the one to implement**, and
`refs/Refs.md` already says so. Its spaces are MEQ's spaces: `HDG_k`, all of
`V_h`, `W_h`, `M_h` at degree `k`, `τ` constant and `O(1)`, simplices.

Its predecessor (Cockburn, Singler & Zhang, JSC **79** (2019) 1777–1800)
interpolates `F( u_h )` into `W_h = P^k` and **loses the superconvergence of the
postprocessed solution**. Remark 2.1 is the whole difference and it is one line:
interpolate `F( u*_h )` into `Z_h = P^{k+1}` instead. For MEQ that difference is
not academic — `ψ*` at `k+2` is what four of the estimator's five terms, the
adaptive loop and `OutputConvergence` all rest on.

`refs/SuperconvergentHDG-II.pdf` is context only: HHO-flavoured spaces and
`τ = 1/h`, both of which MEQ would have to change.

### 1.2 The equation, in MEQ's convention

From `CLAUDE_HDGGS.md`:

```
q − (1/r) ∇̄ψ = 0,      −∇̄·q = F( r, z, ψ ) / r,      ψ = 0 on Γ
F := μ₀ r² dp/dψ + g dg/dψ
```

`DarcyForm` holds **`−q`**, and the potential right-hand side is assembled as
`−( F/r, w )`. Both conventions are load bearing and both are settled in
`CLAUDE_HDGGS.md`, *The assembled flux is −q*.

### 1.3 The three spaces MEQ has, and the fourth it needs

Built at `src/meq/GradShafranov.cpp:592-600`:

| | |
|---|---|
| `V_h` flux | `L2_FECollection( k, 2, GaussLobatto )`, `vdim = 2`, `Ordering::byNODES` |
| `W_h` potential | `L2_FECollection( k, 2, GaussLobatto )`, `vdim = 1` |
| `M_h` trace | `DG_Interface_FECollection( k, 2 )`, `vdim = 1` |

**New**: `Z_h`, the enriched/postprocessing space,
`L2_FECollection( k+1, 2, GaussLobatto )`, `vdim = 1`. MEQ does **not** own one
today — `DarcyForm::Reconstruct()` builds its own internally on first use
("the primary collections cloned at `k+1`", `GradShafranov.hpp:3165-3171`), and
that one is private to it.

### 1.4 The interpolation operator, and what is interpolated

`I_h` is **elementwise Lagrange interpolation at the finite element nodes of
`Z_h`**. It is not a projection and it is not a quadrature refinement: for
`g` continuous on each element, `I_h g ∈ Z_h` is the unique degree-`k+1`
polynomial agreeing with `g` at those `ns = (k+2)(k+3)/2` nodes.

`u*_h ∈ Z_h` is the **classic** HDG postprocessing — Nguyen, Peraire &
Cockburn eq (25) — solved on each element `K`:

```
( ∇̄ψ*, ∇̄z )_K = − ( iK · (flux block), ∇̄z )_K     for all z ∈ P^{k+1}(K)
( ψ*, 1 )_K    =   ( ψ_h, 1 )_K
```

**For MEQ `iK = r`**, and the sign works out with no flip: the stored block is
`−q`, and `∇̄ψ = r q = −r · (stored)`, so `( ∇̄ψ*, ∇̄z ) = −( r · stored, ∇̄z )`,
which is the formula with `iK = r`. `HDGPostprocessBlocks::SetDiffusionInverse(
Coefficient & )` is where that goes. **The sign must be checked and not
argued** — this is the same trap as `transferredDatum()` in `CLAUDE_HDGGS.md`,
where feeding the wrong flux gives the answer with its sign reversed rather than
a small bias.

Because this problem is a **pure Neumann** one closed by the element average,
and because `F` does not appear in it at all, the map from element coefficients
to the nodal values of `ψ*` is **linear and constant in the state**:

```
γ_K = B11_K α_K + B12_K β_K
```

with `α_K` the element's flux coefficients, `β_K` its potential coefficients,
`B11_K` of shape `ns × na`, `B12_K` of shape `ns × nd` and **rank one** (`ψ*`
sees the potential only through its element average). Both are geometry and
`iK` alone, so they are assembled once per mesh.

**This is NOT what MEQ calls `ψ*` today.** `GradShafranovSolver::postProcess()`
(`src/meq/GradShafranov.cpp:7632-7692`) uses `DarcyForm::Reconstruct()`, the
richer mixed reconstruction, which **lifts the non-linear potential integrators
into its local problem and takes their gradient at the computed potential** — so
its output is not a linear function of the unknowns, and feeding it back into
the source would make every element's reconstruction an implicit local fixed
point. `postprocess_hdg.hpp` says so in as many words. §6.4 says what MEQ does
about having two postprocessed potentials.

### 1.5 The residual and the Jacobian, for MEQ's `F`

The paper's `F` is `F( u )`. MEQ's is `F( r, z, ψ )` with an explicit spatial
dependence and a `1/r` weight. **The explicit spatial dependence costs nothing**:
the `Z_h` nodes have fixed reference positions, so their physical images `x_n`
are geometry and are precomputed alongside `A9` —
`HDGReactionIntegratorBase` stores them as `x_data` and
`NodalReactionFunction::Eval( const Vector &x, const Vector &u, Vector &F )`
takes the physical coordinate. So MEQ's `F( r_n, z_n, γ_n )` drops straight in.

With `A9_K( i, j ) = ( χ_j, φ_i )_K`, shape `nd × ns`:

```
residual, potential row :   R_K = A9_K · [ F( x_n, γ_n ) ]_{n=1..ns}
Jacobian, (1,1) block   :   A9_K diag( F'( x_n, γ_n ) ) B12_K
Jacobian, (1,0) block   :   A9_K diag( F'( x_n, γ_n ) ) B11_K
```

exactly the expressions after eq (9) of paper I. **`A9`, `B11` and `B12` are
assembled once per mesh; nothing about the source is reassembled between Newton
iterates.**

**THE `(1,0)` BLOCK IS NEW TO MEQ AND IT IS THE STRUCTURAL CHANGE.** Today the
source contributes a potential-row × potential-column mass term and nothing
else. Under CCSZ the potential equation depends on the **flux**, because `u*` is
reconstructed from it. `DarcyHybridization::Bg_data` exists for precisely this
and its doxygen names the method; `GetGradBMatrix( el, gradient, B )` is the one
funnel that returns the Jacobian's `(1,0)` block when `Bg_data` is live and the
linear divergence form otherwise.

Upstream's *"A linear interpolatory reaction converges in one Newton step"* is
the falsifier for it, and is worth copying: with `F` affine the whole discrete
problem is linear, so a correct Jacobian converges in **exactly one** step, and
one with the `(1,0)` block missing was measured at five steps at a constant
factor of 4.2e-3. **The residual is bit-identical either way**, so the
discriminating quantity is the iteration count and not the answer.

### 1.6 The `1/r`, and why it must not be interpolated

MEQ assembles `−( F/r, w )`. The shipped `HDGReactionIntegratorBase` builds
`A9( i, j ) = ( χ_j, φ_i )_K` with **no coefficient**, so as it stands MEQ would
have to hand it `F/r` as the nodal function. That is wrong twice:

* **It has a pole on the symmetry axis.** `Z_h`'s Lagrange node set is closed,
  so an element touching `r = 0` — which the free-boundary half-disc mesh of
  `tools/mesh/halfdisc.py` has — carries a node **at** `r = 0` exactly. A Gauss
  rule on a triangle has strictly interior points and never meets it. MEQ
  already refuses a source that does not vanish on the axis, so the value there
  is `0/0`, which is a `NaN` and not a limit.
* **It destroys the one exactly-representable case.** `F` for Solov'ev is
  `−( (1−A) r² + A )`, a quadratic, so `I_h F = F` exactly for every `k ≥ 1`.
  `F/r = −( (1−A) r + A/r )` is not a polynomial at any degree.

**So the right form is to weight `A9`**: `A9( i, j ) = ( χ_j / r, φ_i )_K`,
still constant, still assembled once, with `F` itself interpolated. That is
§7.1, a small upstream request.

---

## 2. Verdict

**Decided by the experiment in §0.1, with a strong independent reason to do it
anyway.**

* The **accuracy** question is answered for the smooth fixed-boundary ladder and
  the answer is *no order is lost* — proved in paper I and measured by MFEM on
  triangles at `1.00, 3.00, 4.01, 5.00` for `k = 0..3` with an `O(1)` `τ`, which
  is `k+2` for `k ≥ 1`. It is **open** at a plasma edge that cuts an element,
  where §5.2's arithmetic estimates a loss of one order and §0.2 is the test.
* The **performance** question is open and §0.1 bounds it. §3.3's arithmetic
  says the prize is much smaller than the abstract's "all quadratures performed
  once" suggests, because MEQ's `Z_h` node count is within a factor of two of
  its quadrature point count.
* The **cost** is low and concentrated: the library half is built, and MEQ's
  half is one adapter class, one branch in `buildForms()`, and — the real work —
  **five** existing quadrature loops that must acquire interpolatory twins
  (§6.2).
* The **structural** case (§0.3) stands whatever §0.1 says.

**Not worth doing** only if §0.1 shows no measurable time in the source
quadrature **and** the device path is abandoned.

---

## 3. What it buys

### 3.1 The dimensions, exactly

Straight-sided triangle, `neq = 1`, `dim = 2`. `nd` and `na` agree with M-80's
level-3 table at `k = 2` (`na = 12`, `nd = 6`), which is the cross-check that
this arithmetic is about MEQ.

| `k` | `nd = (k+1)(k+2)/2` | `na = 2 nd` | `ns = (k+2)(k+3)/2` | rule order `2k+4` | `N_q` | `N_q / ns` |
|---|---|---|---|---|---|---|
| 1 | 3 | 6 | 6 | 6 | 12 | 2.00 |
| 2 | 6 | 12 | 10 | 8 | 16 | 1.60 |
| 3 | 10 | 20 | 15 | 10 | 25 | 1.67 |
| 4 | 15 | 30 | 21 | 12 | 33 | 1.57 |

*Arithmetic, not measurement.* `N_q` is read out of MFEM's own triangle table in
`fem/intrules.cpp` (`case 6: IntegrationRule(12)`, `case 8: (16)`,
`case 10: (25)`, `case 12: (33)`); `OrderW() = ( order − 1 ) · dim = 0` for a
straight triangle from `fem/eltrans.cpp:493`.

**On a CURVED mesh the geometry is order 2, so `OrderW() = 2`, the rule order is
`2k + 6` and `N_q` at `k = 2` is 25** — `N_q / ns = 2.50`. The interpolatory
route is relatively better there, and the curved and free-boundary paths are
where MEQ actually runs.

### 3.2 THE UNCOMFORTABLE ARITHMETIC: at `extraOrder = 0` INTERPOLATION EVALUATES `F` MORE OFTEN THAN QUADRATURE DOES

Same table with the minimal `2k` rule:

| `k` | rule order `2k` | `N_q` | `N_q / ns` |
|---|---|---|---|
| 1 | 2 | 3 | **0.50** |
| 2 | 4 | 6 | **0.60** |
| 3 | 6 | 12 | 0.80 |
| 4 | 8 | 16 | 0.76 |

**So the "fewer evaluations of `F`" argument holds only because MEQ deliberately
over-integrates by four orders.** `SourceIntegrator`'s own comment gives the
reason — `F` is not a polynomial in `ψ` at all, Example 5 carries `exp(−ψ)` —
and `setSourceQuadratureOrder()` exists because `PlasmaEdgeConvergence` measured
half an order of `ψ*` riding on it. This is recorded because it is the single
most likely way to be wrong about the prize, and it is the reason §0.1's
experiment is phrased as a derivative rather than as a yes/no.

### 3.3 Work per element per evaluation, at `k = 2`

`nd = 6`, `na = 12`, `ns = 10`, `N_q = 16`. *Arithmetic, not measurement.*

| | quadrature (today) | interpolatory |
|---|---|---|
| **residual** | 16 × `Source::f()`; 16 × `CalcShape` (6 values); 16 × `Tr.Transform`; 96 MAC of `AXPY` | 10 × `Source::f()`; `γ = B11 α + B12 β` ≈ 136 MAC (`B12` is rank one); `A9 F` = 60 MAC. **No `CalcShape`, no `Transform`** |
| **Jacobian** | 16 × `Source::dFdPsi()`; 16 × `CalcShape`; 16 × `Tr.Transform`; 16 × `AddMult_a_VVt` = 576 MAC | 10 × `Source::dFdPsi()`; `γ` ≈ 136 MAC; `A9 diag(F') B11` = 720 MAC; `A9 diag(F') B12` ≈ 96 MAC |
| **extra** | — | an `nd × na` = 72-entry `Bg_data` block per element per linearisation, **replacing** rather than adding to the shared linear `B` at `ComputeElementH()` and `MultInv()` |

**Read honestly, this says three things.**

1. The **residual** is unambiguously cheaper: fewer `F` calls *and* no geometry
   work at all.
2. The **Jacobian** does **more** dense flops — roughly 950 against 576 — and
   fewer `dF/dψ` calls, 10 against 16. Whether that is a win depends entirely on
   the cost of one `meq::Source::dFdPsi()` relative to a multiply-add. For
   `meq::SplineProfile`, with an interval search and a Hermite cubic behind two
   virtual calls, it is tens of MACs and the trade is good; for the analytic
   `Soloviev` source it is a handful and the trade is bad.
3. **`ComputeH()` is untouched**, and M-80 measures it at **68–73% of the
   gradient leg on every row, in both arms, across three sources and three mesh
   sizes**. The interpolatory method removes assembly, not the element-local
   dense factorisation and Schur complement.

### 3.4 The Amdahl ceiling, from M-80's published shares

*Arithmetic on measured shares, not a new measurement.* Threaded, MKL unpinned:

| case | residual | gradient | `ComputeH` as % of gradient | gradient − `ComputeH` | ceiling if the source term became free |
|---|---|---|---|---|---|
| example5 | 16.0% | 33.3% | 71.8% | 9.4% | **≤ 25.4%** |
| pedestal | 14.5% | 33.4% | 73.4% | 8.9% | **≤ 23.4%** |
| highbeta (bordered) | 23.8% | 14.4% | 70.2% | 4.3% | **≤ 28.1%** |

**And these are ceilings with a lot of slack in them**, because the residual leg
is `NPCResidual` in full — the flux mass, the divergence form, the HDG face
terms, the trace row and the local dense algebra — of which the source is one
part. No attempt is made here to guess that fraction; §0.1 measures it.

**The bordered row is the one MEQ runs**, and it is the most favourable: M-80
records the border spending **four residual evaluations per step against one
gradient and one factorisation**, so a residual-side saving is worth four times
a gradient-side one there.

### 3.5 What is NOT bought, said plainly

* **No reduction in trace-solve time**, factorisation or backsolve.
* **No reduction in `ComputeH()`**, which is the largest single item in the
  gradient leg.
* **No reduction in the once-per-mesh assembly.** `prepare()` gains
  `HDGPostprocessBlocks::Assemble()` and `HDGReactionIntegratorBase::Assemble()`;
  both are new work paid once per mesh, and on an adaptive run that is once per
  cycle.
* **No change to the answer where both converge** — which is a *claim to be
  tested*, not a property. Unlike `AssemblyMode` and `TraceSolver`, this **is a
  discretisation change** and may legitimately move the answer at `O(h^{k+2})`.
  IH-8 draws the consequence for the TOML schema.

---

## 4. The design

### 4.1 Both routes coexist, and that is not a hedge

MFEM ships `HDGInterpolatoryReactionIntegrator` and
`HDGQuadratureReactionIntegrator` on one base **deliberately**, so that the two
differ in the treatment of `F` alone and a disagreement is a statement about the
interpolation. MEQ should do the same at its own level: keep
`meq::SourceIntegrator` exactly as it is and add the interpolatory route behind

```
enum class SourceTerm { Quadrature, Interpolatory };
```

with `Quadrature` the default until §0.1 and §0.2 have both reported. Every
stage below is then a controlled comparison against a running control, which is
the only way any of its acceptance criteria mean anything — and it is what lets
the free-boundary path keep the quadrature route if §0.2 goes badly while fixed
boundary takes the interpolatory one.

### 4.2 Where the integrators go, and the constraint that disappears

Today (`src/meq/GradShafranov.cpp:3395-3434`), on the Newton path everything
potential-shaped is on `Mnl_p`:

| | today |
|---|---|
| source | `Mnl_p->AddDomainIntegrator( SourceIntegrator )` |
| HDG face stabilisation, interior | `Mnl_p->AddInteriorFaceIntegrator( HDGDiffusionIntegrator )` |
| HDG face stabilisation, boundary | `Mnl_p->AddBdrFaceIntegrator( HDGDiffusionIntegrator, fittedMarker )` |

and `GradShafranov.cpp:3396-3418` explains that the two forms cannot be mixed:
face integrators left on `M_p` beside a non-linear source **abort**, and domain
integrators left there are **silently dropped**.

Under CCSZ, following MFEM's own `RunAffine()`:

| | interpolatory |
|---|---|
| source | `darcy->GetBlockNonlinearForm()->AddDomainIntegrator( reaction )` |
| HDG face stabilisation, interior | `darcy->GetPotentialMassForm()->AddInteriorFaceIntegrator( ... )` — the **linear** `M_p` |
| HDG face stabilisation, boundary | `darcy->GetPotentialMassForm()->AddBdrFaceIntegrator( ..., fittedMarker )` |

`EnableHybridization()` tests `if (M_p)` first, so the face integrators become
`c_bfi_p`, the linear constraint assembled once per solve — which is the route
`CLAUDE.md` records MEQ already takes through `FaceIntegratorsAreLinear()`, so
this is at worst neutral and at best removes the predicate from the picture
entirely. **Verify it rather than assume it**: `batchedPotFaceAssemblyTaken()`
already exists as an accessor and should read the same in both arms.

### 4.3 The flux mass stays factored once, and the promise is how

`HDGReactionIntegratorBase::GetBlockRowMask()` returns `1 << 1` — the potential
row alone. `DarcyHybridization::FluxMassIsPrefactored()` is true exactly when
nothing makes the flux row non-linear, **and a block non-linear integrator whose
mask promises the potential row alone is explicitly admitted**. So MEQ keeps the
`LocalOpType::PotNL` regime, where the flux mass is factored once at
`Finalize()`.

*"The row READS the flux — `gamma = B11 u_e + B12 p_e` — and that is a different
statement: the promise is about what is written."* Upstream's own test
`GetBlockRowMask() keeps the flux mass factored once` runs the identical problem
with the promise withdrawn and requires the two to agree to 1e-12, so the mask
is a performance declaration and not a discretisation one. MEQ should assert
`FluxMassIsPrefactored()` in the interpolatory arm for the same reason MEQ
asserts `GetNumLocalNLIterations() == 0` for NPC: it is the signal that the
intended route was taken.

### 4.4 `HDGPostprocessBlocks` and the freshness contract

One instance, owned by the solver, constructed on `( fluxFes, potentialFes,
enrichedFes )`, `SetDiffusionInverse( radius )`, `Assemble()`.

**It records `Mesh::sequence` and aborts on a mismatch**, and it exposes
`CoefficientsMoved()` for a moving diffusivity. MEQ's `iK = r` never moves, so
only the mesh matters — and the adaptive loop refines in place, which bumps the
sequence. `buildForms()` is called from `prepare()` and `prepare()` runs after
every refinement, so assembling the blocks there is correct. **The abort is a
feature**: it turns "stale blocks after a refinement" from a wrong answer into a
named failure.

---

## 5. Accuracy

### 5.1 On smooth problems there is no order to lose, and this is proved AND measured

Paper I, Corollary 3.15 and Theorem 3.19:

```
‖q − q_h‖ ≤ C h^{k+1},   ‖u − u_h‖ ≤ C h^{k+1},   ‖u − u*_h‖ ≤ C h^{k+1+min{k,1}}
```

i.e. `k+2` for `k ≥ 1` and **no superconvergence at `k = 0`**, which the linear
case shares. Table 1 reads 2.95, 3.02, 3.02, 3.01 for `u*` at `k = 1`.

**Why the consistency term does not cost an order.** The interpolation enters
the error bound through `‖F(u) − I_h F(u)‖_{T_h}` in eq (27), linearly. For a
smooth `F(u)`, interpolation into `P^{k+1}` gives `C h^{k+2} |F(u)|_{k+2}` —
**exactly the target order for `u*`, and one order better than the target for
`u_h` and `q_h`.** That is the whole content of Remark 2.1: interpolating into
`W_h = P^k` instead gives `h^{k+1}`, which pollutes `u*` at `h^{k+1}` and kills
the superconvergence. The choice of `Z_h` is not a convenience, it is the
theorem.

MFEM has already measured it on the configuration MEQ would use — triangles,
`O(1)` `τ`, `−Δu + u³ − u = g` — at **1.00, 3.00, 4.01, 5.00** for `k = 0..3`,
against **−0.25, 1.96, 3.08** with the default `1/h` stabilisation. MEQ runs
`meq::ConstantStabilization` at `τ = 1`, so it is on the right side of that.

**Two coincidences to avoid.** First, on a tensor-product element with MFEM's
default Gauss-Legendre L2 basis the interpolatory term **is** the quadrature
term identically — interpolation degenerates to collocation, measured at a
relative gap of 4.3e-16 on a box against 1.2e-03 on a triangle. MEQ's meshes are
triangles and its basis is Gauss-Lobatto, so it is safe from this, but any
comparison written on a quadrilateral fixture would pass for any implementation
whatsoever. Second, on a Gauss-Legendre box the gap converges at `h^{k+4}` by
orthogonality; a Gauss-Lobatto basis recovers `h^{k+2}` exactly. **Build `Z_h`
Gauss-Lobatto**, matching MEQ's other volume spaces, and expect the full
`h^{k+2}` consistency error rather than a lucky cancellation.

### 5.2 THE PLASMA EDGE IS WHERE THE ORDER IS AT RISK, AND THE ESTIMATE SAYS ONE ORDER

*This subsection is arithmetic in the style of M-09, not a measurement.*

Let `j` be the order to which the profiles vanish at the plasma edge, `p' ~ Ψ^j`
— a modelling choice, `j = 1` for FreeGS's default. `meq::SplineProfile` clamps
outside its knots and `NormalisedSource::insidePlasma()` is a hard pointwise
gate, so across the edge `F` behaves like `c d_+^j` with `d` the signed
distance.

On an element the edge crosses, `d = O(h)`, so `F = O(h^j)` there and its
`P^{k+1}` Lagrange interpolation error is `O(h^j)` pointwise — the *constant*
falls with `k`, the order does not, which is the same mechanism M-09 records for
the best approximation of `|d|^m`. Then

```
‖F(u) − I_h F(u)‖_{L²(K)}      = O( h^{j} · h )      = O( h^{j+1} )
‖F(u) − I_h F(u)‖_{L²(band)}   = O( h^{j+1} / h^{1/2} ) = O( h^{j+1/2} )
```

over the `O(1/h)` cut elements. Against M-09's measured `min( k+1, j+1.5 )` for
`ψ_h` under the plain Gauss rule, **that is one order lower.** At `j = 1,
k = 2` it is `h^{1.5}` against a measured 2.58.

**Two honest qualifications, and they point in opposite directions.**

* The bound is crude. Theorem 3.14 carries the consistency term in `L²` with no
  duality argument, and the interpolation error of a kinked function oscillates
  in sign, so the observed rate can be better than the bound. §0.2 is the
  experiment and `PlasmaEdgeConvergence` already has fixtures at
  `j = 0, 1, 2, 3`.
* The same arithmetic applied to the quadrature route would also predict
  `h^{j+1/2}`, and MEQ measures `h^{j+1.5}`, so the method is known to be
  optimistic by an order **on the route it can be checked against**. That is an
  argument for running the test and not for trusting either number.

**The consequence if it goes badly is bounded**, because §4.1 keeps both
integrators: fixed boundary and the smooth fixed-`Ψ_edge` cases take the
interpolatory route, the confined free-boundary path keeps the quadrature one,
and `SourceTerm` is the switch. That is a worse outcome than a clean win and a
much better one than a dead port.

### 5.3 A node landing on the plasma edge is a new, small hazard

Under quadrature, `insidePlasma()` is evaluated at `N_q = 16` interior Gauss
points; under interpolation, at `ns = 10` nodes, **some of which lie on element
boundaries** because `Z_h` is a closed Lagrange space. As Newton moves the edge
across a node, one entry of `F(γ)` switches between the profile value and zero.

For `j ≥ 1` this is continuous — `F` itself is continuous across the edge, only
its derivative kinks — so the residual stays continuous in the unknowns and
there is a Jacobian to iterate with. That is the same condition
`Source.hpp:304-311` already imposes and already measures: violate it and Newton
fails at every degree and mesh. At `j = 0` `CLAUDE_FB.md` records that Newton
cannot solve the problem at all, so this changes nothing there.

**What does change is that a boundary node is shared between two elements**, and
the two elements may disagree about which side of the edge it is on if their
`γ` differ there — which they do, `Z_h` being discontinuous. That is correct
behaviour (each element interpolates its own `F`) and is called out only so it
is not mistaken for a bug when the plasma-support count moves by one element.

### 5.4 The Lipschitz condition, and whether MEQ's sources satisfy it

Paper I proves superconvergence under a **global** Lipschitz condition
(Assumption 3.4) and then under a **local** one (Assumption 3.16) for `k ≥ 1`
on a quasi-uniform mesh with `h` small enough. Locally Lipschitz is satisfied by
every smooth `F`, so the analytic ladder — `Soloviev`, `McCarthy`,
`ManufacturedNonlinear`, `SimilarityExponential` — is covered.
**`NormalisedMHDSource` with `ConfineToPlasma` is locally Lipschitz in `ψ` for
`j ≥ 1` and is not for `j = 0`**, which is one more statement of the same
threshold.

The **quasi-uniform** hypothesis is worth noting against the adaptive loop,
which grades the mesh deliberately. That is a gap in the theory and not a reason
to expect a failure; it is why AD-2 in §6.5 measures the adaptive rates rather
than assuming them.

---

## 6. What changes in MEQ, file by file

Nothing below is written yet. Function names are the existing ones.

### 6.1 New, and small

| | |
|---|---|
| `src/meq/GradShafranov.cpp:592-600` | build `enrichedColl` = `L2_FECollection( order+1, dim, GaussLobatto )` and `enrichedFes` beside the three existing spaces |
| new member | `std::unique_ptr<mfem::HDGPostprocessBlocks> postprocessBlocks`, `SetDiffusionInverse( radius )`, `Assemble()` from `buildForms()` |
| **new class** `meq::NodalSource : public mfem::NodalReactionFunction` | the adapter. `Eval()` returns MEQ's source at one node; `EvalJacobian()` returns `dF/dψ`. Holds `Source const *`, the cached `NormalisedSource const *`, and the plasma-component mask — see §7.2 for why the last one does not fit the interface as shipped |
| `buildForms()`, `GradShafranov.cpp:3395-3440` | a third branch on `SourceTerm::Interpolatory`, laying the integrators out as §4.2 says. **`react->Assemble()` before `AddDomainIntegrator()`**, because `BlockNonlinearForm` takes ownership |
| `GradShafranov.hpp` / `.cpp` | `setSourceTerm( SourceTerm )`, `sourceTerm()`, and an accessor for `FluxMassIsPrefactored()` so a test can assert the route |

**The sign.** MFEM's `AssembleElementVector` writes `+A9 F(γ)` into the
potential row; MEQ's `SourceIntegrator` writes `−w F/r`. So `meq::NodalSource`
must return the negated, `r`-weighted quantity, and **which negation is correct
must be established by running both**, exactly as upstream's own fixture comment
says ("The sign is the form's, established by running both"). IH-2's `McCarthy`
acceptance is what settles it: `F` linear in `ψ` means a sign error is a
different equation and a wrong answer, not merely slow convergence.

### 6.2 THE REAL WORK: FIVE EXISTING QUADRATURE LOOPS THAT DIFFERENTIATE THE ASSEMBLED RESIDUAL

`src/meq/GradShafranov.cpp:4480-4489` states the rule and it is the thing most
likely to be missed:

> ALL THREE USE `meq::SourceIntegrator`'S OWN QUADRATURE RULE, and they have
> to: these are derivatives of the ASSEMBLED residual, not of the continuous
> one, so a different rule would differentiate a different function.

There are in fact **five**, all sharing the free function `sourceRule()` at
`GradShafranov.cpp:4494-4500`:

| | |
|---|---|
| `assemblePlasmaCurrent()` | `:4533` |
| `assembleCurrentColumn()` | `:4585` |
| `assembleCurrentNormalisationCorner()` | `:4646` |
| `assembleCurrentRow()` | `:4704` |
| `assembleNormalisationColumn()` | `:4919`, rule at `:4975` |

**Under the interpolatory term the assembled residual is a different function,
so every one of these must acquire an interpolatory twin or the borders are
differentiating something the solve is not solving.** The twins are simpler than
the originals, not harder: the normalisation column becomes
`A9 · [ ∂F/∂s ( x_n, γ_n ) ]`, using the closed-form
`NormalisedSource::normalisationDerivatives()` that already exists, and the
plasma current becomes `b2ᵀ`-weighted nodal sums.

**This is the hazard the project already has a name for.** A wrong border column
does not produce a wrong answer, it produces a Newton that stalls — and
`CLAUDE_HDGGS.md`, *A wrong Jacobian is invisible to a convergence table*, is
the standing warning. Every twin needs the treatment
`the_normalisation_derivatives_are_analytic` already gives the analytic column:
pinned against a Richardson-extrapolated difference of the **interpolatory**
residual.

### 6.3 `postProcess()` changes without being edited, and that must be measured

`GradShafranovSolver::postProcess()` (`:7632-7692`) calls
`DarcyForm::Reconstruct()`. Today that lifts `Mnl_p`'s non-linear potential
integrators into its local problem and takes their gradient at the computed
potential, which is why `thePostProcessedPotentialSurvivesNewton` reads
3.05 / 4.05 / 5.03. Under the interpolatory branch **`Mnl_p` has no integrators
at all**: the source is on the block non-linear form, which
`ReconstructFluxAndPot()`'s whitelist does not consult, and the face terms are
back on `M_p`.

**Two consequences, in opposite directions.**

* `ψ*` becomes a **different function**. Its rate must be re-measured; there is
  no reason to expect it to fall, since CCSZ's own `u*` has no reaction term
  either and reaches `k+2`, but "no reason to expect" is not an acceptance.
* The singular-matrix defect `CLAUDE_HDGGS.md`, *Post-processing is back*,
  records **cannot arise**, because it needed a non-linear potential integrator
  whose Jacobian could be the zero matrix and there is now no such integrator.
  `thePostProcessedPotentialIsCorrectWhereTheJacobianVanishes` should be run in
  the interpolatory arm and is expected to be trivially green there.

### 6.4 Two postprocessed potentials, and they are not the same function

After the port MEQ holds:

* `ψ*` from `DarcyForm::Reconstruct()` — reported in the output, used by
  `meq::Estimator` in four of eq (20)'s five terms, and used by the adaptive
  loop;
* `u*` from `HDGPostprocessBlocks` — internal to the source term, never
  written out, existing only as the `γ` vector inside an element.

`postprocess_hdg.hpp` and `refs/Refs.md` both say the two answer different
questions. **Do not unify them**, and in particular do not feed `Reconstruct()`'s
`ψ*` back into the source: it is not linear in the unknowns, so `B11` and `B12`
would not be constant and each element's reconstruction would become an implicit
local fixed point. The plan is to keep both and to measure, once, how far apart
they are on a converged Solov'ev solution — a number worth knowing and cheap
(IH-1).

### 6.5 What does not change

`src/meq/Source.{hpp,cpp}` — no change to `meq::Source`,
`meq::NormalisedSource` or any concrete source. `f()` and `dFdPsi()` are already
pointwise in `( r, z, ψ )` and that is exactly what a `NodalReactionFunction`
wants. This is the payoff of `CLAUDE.md`'s rule that `Profiles` and `Source`
keep MFEM out.

`src/meq/Config.{hpp,cpp}` — **no new key**, see IH-8.

`src/meq/Estimator.hpp`, `Field`, `Sampler`, `Output`, `WarmStart`,
`CriticalPoints`, `FluxSurfaces` — untouched.

---

## 7. What is needed from MFEM

Per `CLAUDE.md`, upstream work is requested by **writing a document into
`../mfem-hdg-dev/doc/` and doing nothing else in that tree** — no branches, no
commits, no builds. Both items below are small and both are in the interpolatory
machinery upstream has just built, so per the standing rule *do not file findings
against unfinished work* they should be raised only when the port reaches them
and can say what it measured.

### 7.1 A weight coefficient on `A9`

**Ask**: `HDGReactionIntegratorBase::SetWeight( Coefficient & )`, so that
`A9( i, j ) = ( c χ_j, φ_i )_K` with `c` defaulting to 1.

**Why**: §1.6. Without it MEQ must interpolate `F/r`, which has a pole at `r = 0`
on the half-disc mesh that the free-boundary path uses, and which turns the one
exactly-representable source in the tree — Solov'ev, a quadratic — into an
inexactly-represented one at every degree.

**Size**: the coefficient is evaluated inside the existing quadrature loop in
`HDGReactionIntegratorBase::Assemble()`, once per mesh. It cannot move any
existing answer, `c = 1` being the default.

**Workaround if refused**: MEQ can special-case `r < ε` in
`meq::NodalSource::Eval()` and return zero there, which is correct wherever the
source vanishes on the axis — which MEQ already refuses to run without. That
loses the Solov'ev exactness and is a worse place to put the knowledge.

### 7.2 The element index in `NodalReactionFunction::Eval()`

**Ask**: pass the element number, either as an argument or through a
`SetElement( int )` the integrator calls before its node loop.

**Why, and this one BLOCKS the diverted free-boundary path.**
`meq::SourceIntegrator::sourceValue()` consults **two** gates: the pointwise
`NormalisedSource::insidePlasma( ψ )`, and the **per-element** XP-1 flood-fill
mask `elementCarriesPlasma( tr.ElementNo )`, which is what keeps a private-flux
region from carrying a current channel nobody asked for —
`CLAUDE.md` and M-82 record what happened the one time that fill did not run.
`NodalReactionFunction::Eval( const Vector &x, const Vector &u, Vector &F )`
receives only the physical coordinate, so **the element mask cannot be
expressed through the interface as shipped**, and neither can
`fOutsidePlasma( r, z )`, which is selected by the same mask.

The integrator has `Tr.ElementNo` in hand at both entry points; this is a
signature change and not a mechanism.

**Workaround if refused**: none that is acceptable. Locating a node in the mesh
to recover its element is absurd for a quantity the caller already has. Until
this lands, the interpolatory route is **fixed-boundary and limited-plasma
only**, which is a sensible place to stop anyway and is where §8's stages stop.

### 7.3 What is NOT needed

* **No new integrator.** `HDGInterpolatoryReactionIntegrator` is it.
* **No new assembly route.** `Bg_data`, `GetGradBMatrix()`, `GetBlockRowMask()`
  and the batched twins are all built.
* **No NPC change** — but see §9.1, which is a testing gap rather than a code
  gap.

---

## 8. The stages

House format: each ends at a **measured convergence rate or a measured timing**,
not at "it runs". Every stage runs against the quadrature control of §4.1.

### IH-0 — the falsifier (§0.1), before anything is written

**Do**: add `--source-extra N` to `tests/performance/NewtonStepProfile.cpp`,
rebuild `-j6`, run `extraOrder ∈ { 4, 2, 0 }` × three cases, medians of five,
`MKL_NUM_THREADS=1`, idle machine, load gate as M-80 uses.

**Acceptance**: `T(16) − T(6)` reported per case with its scatter, and
`( 16/10 ) × [ T(16) − T(6) ]` published as the port's upper bound. **A result
inside the scatter is a decision to stop** unless §0.3 is taken up on its own.
Also report the `SolovievConvergence` and `NewtonConvergence` rates at
`extraOrder = 0`, since the same runs produce them.

### IH-1 — `u*` alone, no solver change

**Do**: construct `HDGPostprocessBlocks` on a **converged** MEQ Solov'ev
solution, with `iK = r`, and compare its `u*` against the exact `ψ`. Also report
`‖u* − ψ*‖` against `Reconstruct()`'s output (§6.4).

**Acceptance**: `u*` converges at **`k+2`** for `k = 1, 2, 3` on the dyadic
Solov'ev ladder — 3, 4, 5 to within the suite's usual slack. **A rate of `k+1`
means the sign or `iK` is wrong**, which is the failure mode
`CLAUDE_HDGGS.md` records for `transferredDatum()`, and is the whole reason this
stage exists separately: it isolates the postprocessing from the source before
either can hide the other.

### IH-2 — the interpolatory source, fixed boundary, the analytic ladder

**Do**: `meq::NodalSource`, the `buildForms()` branch, `SourceTerm`. Run the
ladder in `tests/analytic/`, both arms.

**Acceptance**, in the order the ladder was built for:

| fixture | `∂F/∂ψ` | what it must show |
|---|---|---|
| `Soloviev` | `0` | **one** Newton step, and with §7.1 in place `I_h F = F` exactly at `k ≥ 1`, so `ψ_h` and `q_h` must agree with the quadrature arm to **round-off** rather than merely to the same rate. Without §7.1 the acceptance weakens to the rate |
| `McCarthy` | constant `T` | **one** Newton step and the closed form to the published tolerance — this is upstream's own affine falsifier and is what settles the sign and the `(1,0)` block |
| `ManufacturedNonlinear` | varies | rates `k+1` in `ψ_h` and `q_h`, `k+2` in `ψ*`, `k = 1, 2, 3`, and quadratic Newton |
| `SimilarityExponential` | `nF` | the same against an *exact* nonlinear solution |

Plus: `FluxMassIsPrefactored()` true, `GetNumLocalNLIterations()` zero,
`batchedPotFaceAssemblyTaken()` unchanged between arms.

### IH-3 — the timing, on the same cases as IH-0

**Do**: `NewtonStepProfile` again, both arms, same protocol.

**Acceptance**: the measured step-time change per case, published beside IH-0's
bound, with the residual and gradient legs split out. **IH-3 must come in under
IH-0's bound**; if it does not, one of the two measurements is wrong and the
discrepancy is the finding.

### IH-4 — the bordered Newton (§6.2)

**Do**: the interpolatory twins of the five loops. Pin each against a
Richardson-extrapolated difference of the interpolatory residual.

**Acceptance**: `HighBetaConvergence` green in the interpolatory arm; `ψ_ax`
agreeing with the quadrature arm to the discretisation difference and not
merely to the solver tolerance; the final residual reaching the analytic
column's floor (M's `5.70e-16` order) and not the differenced column's
(`1.69e-13`). **A residual that floors three orders high means a twin is
differentiating the old assembly.**

### IH-5 — the curved boundary and the adaptive loop

**Do**: `miller-curved` and `miller-adaptive` in both arms.

**Acceptance**: the extension benchmark's rates unmoved within
`ExtensionConvergence`'s own 0.30-per-pair slack; `η` monotone over the adaptive
cycles; the element counts within one refinement of the control. Note that
`HDGPostprocessBlocks` must have been re-`Assemble()`d after every refinement —
its `Mesh::sequence` check turns a failure to do so into an abort, so this stage
also tests §4.4.

### IH-6 — the plasma edge (§0.2), and it may be a refusal

**Do**: `PlasmaEdgeConvergence`'s M-10 table in the interpolatory arm.

**Acceptance**: `k+2` retained on every `k ≤ j` entry. **If not, the finding is
recorded and the free-boundary path keeps `SourceTerm::Quadrature`** — which is
not a failure of the stage, it is the stage doing its job. Either way the table
is published under a new `M-nn` anchor beside M-10.

### IH-7 — free boundary, and it is gated on §7.2

**Do**: nothing until the element index reaches
`NodalReactionFunction::Eval()`. Then `FreeBoundaryCoupling` and the
`freegs4e` benchmark in both arms.

**Acceptance**: the limited tokamak reproducing `freegs4e` to the same 1.3e-04,
and the diverted case's flood fill covering the same element set. Until §7.2
lands, `setSourceTerm( Interpolatory )` must **throw** when the source is
confined to the plasma, rather than silently running the fill's gate off.

---

### IH-8 — the schema question: no TOML key until every stage has reported

`CLAUDE.md` draws the line at whether a key can change the answer:
`AssemblyMode` and `TraceSolver` are exposed because they cannot;
`Globalisation` and `NonlinearOrdering` are withheld because they report
different discrete solutions. **`SourceTerm` is a discretisation change and
therefore changes the answer by construction**, at `O(h^{k+2})` on a smooth
problem and possibly by an order at a plasma edge. It stays a library setter
until IH-2 through IH-6 have all reported, and it is exposed then only if they
agree — with the same "and it was measured" gloss `docs/` gives the other
choices.

---

## 9. Interactions with the three things that make MEQ unusual

### 9.1 NPC, and the element-local Newton's frozen guess

**The frozen-guess trap does not apply on the NPC path, and MEQ's default is
NPC.** `CLAUDE_HDGGS.md`, *The trap that cost the most*, is about
`DarcyHybridization` capturing the element-local Newton's initial guess at
`FormLinearSystem()` time; under NPC there is no element-local non-linear solve
at all, which is what `GetNumLocalNLIterations() == 0` asserts. The whole
`formSystem()`-per-accepted-step apparatus went away with the port to NPC and
nothing here brings it back.

**It DOES apply to `CondenseThenLinearise`**, the backup, where the block
non-linear integrator makes every element's elimination non-linear — and
upstream's `RunAffine()` even had to raise the local solve tolerances
(`SetLocalNLSolver( Newton, 1000, 1e-13, 1e-16 )`) to stop the local rtol of
1e-6 capping the outer Newton at ~1e-5. **MEQ must do the same in the
interpolatory arm if it runs the condensation at all**, and `SolverContract`'s
`theOrderingsAgreeAndOnlyOneIteratesLocally` is where that is checked.

**AND THE GAP WORTH NAMING: UPSTREAM HAS NO NPC TEST FOR THE INTERPOLATORY
TERM.** `tests/unit/fem/test_darcy_reaction.cpp` drives everything through
`FormLinearSystem()` and a `NewtonSolver` on the reduced operator — the
condensation. Reading the code, `DarcyNPCOperator::Mult` →
`DarcyHybridization::NPCResidual` → `MultNL( MultNlMode::AtFields )` →
`LocalResidual` → `LocalNLOperator::AddMultBlock`, which begins
`if (!dh.m_nlfi && !dh.c_nlfi) { return; }` and then runs the block integrator —
so the path exists. **But "the path exists" is not "the path is tested", and
MEQ would be its first user.** That is exactly the position
`CLAUDE.md` describes for `SetEssentialBC` with a non-linear reduced operator,
and the standing advice applies: if IH-2 produces a converged-but-wrong answer,
look here before looking at MEQ's assembly.

### 9.2 The bordered Newton on `ψ_ax`

Three separate points.

* **The five loops (§6.2) are the work**, and getting one wrong is invisible to
  a convergence table.
* **The border is where the prize lands.** M-80: the bordered row spends four
  residual evaluations per step against one gradient and one factorisation, and
  the residual is the leg the interpolatory form unambiguously improves (§3.3).
  `highbeta` is therefore the case to measure, not `example5`.
* **`b`, the row, is untouched.** `∂(max ψ_h)/∂λ` is the potential shape
  functions of the axis element evaluated at `x*`, and `∇ψ_h( x* ) = 0` by the
  envelope argument. That argument is about `ψ_h`, not about `ψ*`, so the
  interpolatory term does not reach it. **`c`, the column, is entirely rebuilt**
  as `A9 · [ ∂F/∂s ]`, and it gets simpler: the closed-form derivatives already
  exist and return exactly zero outside the plasma.
* **The moving-support argument for an analytic column gets STRONGER.**
  `CLAUDE_HDGGS.md` records that a differenced column straddles a kink when the
  edge moves and floors the iteration at ~3e-09. Under interpolation the
  kink is sampled at 10 nodes instead of 16 points, so a difference is if
  anything worse. `BorderColumn::Analytic` must stay the default and the
  differenced route stays as the control.

### 9.3 Free boundary, the moving support, and the kinked profile

* **§7.2 is a hard gate**: the per-element flood-fill mask cannot be expressed
  through `NodalReactionFunction` as shipped, so the diverted path is blocked
  until upstream passes the element index. `setSourceTerm( Interpolatory )` must
  refuse a confined source until then rather than silently dropping the fill.
* **§5.2 is the accuracy question** and IH-6 is its test. The likely outcome is
  a per-route decision rather than a single answer.
* **§5.3** is the node-on-the-edge hazard: real, small, continuous for `j ≥ 1`,
  and already fatal at `j = 0` for reasons that predate this.
* **`freezePlasmaEdge()` composes unchanged.** XP-2 freezes the edge threshold
  so it does not move mid-iteration; the interpolatory term reads
  `supportAxis()`/`supportBoundary()` through the same
  `NormalisedSource::insidePlasma()` and inherits the freeze.
* **`fOutsidePlasma( r, z )`** — the coil term — is `ψ`-independent, so under
  the interpolatory form it is a constant vector `A9 · [ fOut( x_n ) ]` per
  element and could move to the right-hand side entirely. **Do not do that in
  the same stage as the port**; it is a separate simplification and mixing them
  makes a disagreement unattributable.

---

## 10. Hypotheses considered here and rejected, so they are not re-derived

* **"Interpolate `F( ψ_h )` into `W_h`."** That is the 2019 predecessor, and it
  loses the superconvergence of the postprocessed solution — which for MEQ means
  `ψ*` at `k+1`, the estimator a full order worse, and the adaptive loop
  refining on a degraded indicator. `Z_h = P^{k+1}` is not a refinement of the
  idea, it is the theorem.
* **"Feed `DarcyForm::Reconstruct()`'s `ψ*` to the source, since MEQ already has
  it."** It lifts the non-linear potential integrators and takes their gradient
  at the computed potential, so it is not linear in the unknowns, `B11`/`B12`
  would not be constant, and each element's reconstruction would become an
  implicit local fixed point. §6.4.
* **"Interpolate `F/r`, since that is what the weak form carries."** A pole at
  `r = 0` on the half-disc, and it destroys the Solov'ev exactness. §1.6, §7.1.
* **"The saving is that quadrature of `F` disappears."** At MEQ's degrees
  `ns` is within a factor of two of `N_q`, and at `extraOrder = 0` it is
  **larger**. §3.2. The paper's framing is a time-stepping one — matrices
  assembled once before the integration — and MEQ already reassembles nothing
  but the source.
* **"It removes `ComputeH()`'s cost."** It does not; `ComputeH()` is 68–73% of
  the gradient leg on every row of M-80 and is element-local dense algebra, not
  assembly. §3.3, §3.4.
* **"Paper I measures the speed-up."** It does not. Grepped: the paper reports
  Table 1's convergence history and pattern-formation figures, and **no timing
  of any kind**. Remark 2.3's claim of computational advantage is asserted. The
  only measurement of the gap this project has access to is MFEM's own, in
  `reaction_hdg.hpp`'s doxygen.
* **"MEQ's `AssemblyMode` or `TraceSolver` precedent means `SourceTerm` can be a
  TOML key."** Those cannot change the answer; this one changes it by
  construction. IH-8.

---

## 11. A build-system trap that will outlive this plan

**`gf-interp-hdg-dev` IS A FIFTH BRANCH AND `CLAUDE.md`'S MERGE RECIPE DOES NOT
LIST IT.** The recipe under *`meq-integration`: the branch MEQ builds from*
names four: `gf-hdg-subdomains-dev`, `direct-solver-symbolic-reuse`,
`gf-hdg-linearise-first`, `gf-hdg-dev`. Checked in `../mfem/mfem-src` on
2026-09-15, `meq-integration` contains **six** — those four plus
`gf-interp-hdg-dev`, which introduced `fem/darcy/reaction_hdg.*` in
`db9b64cef9` and `c31d059ff0`, and `gf-hdg-p-adaptivity`. The installed
`reaction_hdg.hpp` is byte-identical to `meq-integration`'s.

**So re-creating `meq-integration` by following the recipe as written would
delete the entire library half of this port**, with the symptom being a
compile failure on a missing header rather than anything naming a branch. If
this plan is taken up, the recipe and its containment loop must grow those two
entries. Recorded here rather than acted on, because this file changes no source
and `CLAUDE.md` is not this document's to edit.

---

## 12. Summary of what is measured, what is cited, and what is arithmetic

| claim | status |
|---|---|
| MFEM implements CCSZ-I; the symbols are in `libmfem.a`; `Bg_data` carries the `(1,0)` block | **verified on disk**, headers + `nm` + `git`, 2026-09-15 |
| `gf-interp-hdg-dev` is contained in `meq-integration` and absent from the recipe | **verified**, `git merge-base --is-ancestor`, 2026-09-15 |
| Paper I reports no timing of any kind | **verified**, full-text search of `refs/SuperconvergentHDG-I.pdf` |
| `u*` rates `k+2` for `k ≥ 1`, no superconvergence at `k = 0`; `‖F(u) − I_h F(u)‖` enters linearly at `h^{k+2}` | **cited**, paper I Cor. 3.15, Thm 3.19, Table 1, eq (27) |
| 1.00 / 3.00 / 4.01 / 5.00 on triangles at `O(1)` `τ`; the Gauss-Legendre box collocation coincidence | **cited**, MFEM `postprocess_hdg.hpp` and `reaction_hdg.hpp` doxygen |
| `N_q` = 12 / 16 / 25 / 33 at `k = 1..4`; `ns` = 6 / 10 / 15 / 21; `OrderW() = 0` on a straight triangle | **arithmetic**, from MFEM's `fem/intrules.cpp` and `fem/eltrans.cpp:493` and MEQ's rule at `GradShafranov.cpp:353-364` |
| the flop counts in §3.3 and the Amdahl ceilings in §3.4 | **arithmetic**, the second on M-80's published shares |
| the plasma-edge consistency estimate `O(h^{j+1/2})` | **arithmetic**, in M-09's style, and known to be an order pessimistic on the one route it can be checked against |
| anything in §0.1, IH-0 or IH-3 | **not run.** No number is printed for any of it |
