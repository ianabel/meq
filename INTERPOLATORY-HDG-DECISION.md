# Where CCSZ-I goes: MFEM machinery or MEQ code

**Written 2026-09-09.** A decision about *placement*, not about whether to adopt
the method — `ROADMAP.md` §12.1 is the plan and the acceptance criteria, and
`refs/Refs.md`'s *Interpolatory HDG* is the bibliographic account. Neither says
where the code lives. This does.

**Line numbers into `fem/darcy/` are against `meq-integration` at
`2613ee1b3b`**, the merge MEQ builds from, read out of `../mfem/mfem-src`.
`../mfem-hdg-dev`'s working copy of `darcyhybridization.cpp` is over a thousand
lines away from it, so the function names are the durable reference.

Nothing is measured here. Every number quoted is read out of a paper or a source
file; the rate claims MEQ would accept the method on **do not exist yet**, and
the section *What MEQ can do today* is what would produce them.

## The decision

**Split, with two asks on `../mfem-hdg-dev` and everything physical on MEQ's
side.** The load-bearing reason is one line of somebody else's private code:
`darcyhybridization.cpp:3208` sets `grad_arr(1,0) = NULL` and the comment above
it says why — *"Block (1,0) stays NULL: the divergence form is linear, so B is
already exact in `Bf_data`"*. CCSZ-I is a counterexample to that sentence inside
the same family of methods, because `u*` depends on `q_h`, so the potential
equation acquires a Jacobian block against the flux; and every place that block
would have to reach — `Bf_data`, `GetBnlMatrix()`, `ComputeElementH()`,
`MultInv()`, `LocalNLOperator::B` — is private to `DarcyHybridization`. **MEQ
cannot supply it from outside at any price short of reimplementing the
hybridized elimination.** So the two asks are (1) `HDGPotentialPostprocessor`
exposing its per-element `B11`/`B12` as matrices rather than only applying them,
and (2) a solution-dependent (1,0) gradient block, mirroring the (0,1) one that
already exists for the mirror-image reason. Everything else — `I_h`, the choice
of nodal basis for `Z_h`, the `1/r` weight, the explicit `(r, z)` dependence,
`ConfineToPlasma`'s cut elements, the plasma-component mask, the acceptance
ladder — is MEQ's, and none of it belongs in a general library. The request is
`../mfem-hdg-dev/doc/CCSZ-INTERPOLATORY-HDG-FROM-MEQ.md`.

## What §2 actually says, and four refinements to the brief

§2 of `refs/SuperconvergentHDG-I.pdf` is as summarised. `V_h`, `W_h`, `M_h` at
degree `k`, `Z_h = P^{k+1}`, `τ` nonnegative, elementwise constant and `O(1)`
(eq 3); the nonlinear term is `(I_h F(u*_h), v_h)` with `I_h` the elementwise
Lagrange interpolation onto `Z_h`'s finite element nodes; the postprocessing is
eqs (4a)–(4b), rewritten with a `P^0` multiplier as (7a)–(7b); `γ = B11 α +
B12 β`; the term is `A9 𝓕(γ)`; the Jacobian blocks are `A9 diag(𝓕'(γ)) B11`
and `A9 diag(𝓕'(γ)) B12`, and they sit in the **potential row** of `F'(x)`,
against the flux and potential columns respectively. Four things to add.

* **`B12` is rank one per element.** The right-hand side of (7) against `β` is
  `[0; b2^ℓ β^ℓ]`, and `b2^ℓ β^ℓ` is the scalar `∫_K u_h`: `u*` sees `u_h`
  **only through its element average**. So `B12^ℓ = c^ℓ (b2^ℓ)^T` for one column
  `c^ℓ`, and the potential–potential Jacobian block of the interpolatory
  reaction term has rank one where MEQ's present one is full. Remark 2.2 says
  `ℓ = 0, …, k−1` all work, which raises the rank to `dim P^ℓ(K)` and is the
  one dial in the method MEQ has a free choice of.
* **Table 1 is the new method, not the old one.** It exhibits the *restored*
  `k+2` — 2.95, 3.02, 3.02, 3.01 for `u*` at `k = 1` — and at `k = 0` it reads
  0.97, no better than `u_h`. That is the theory: Theorem 3.14 gives
  `‖u − u*_h‖ ≤ C h^{k+1+min(k,1)}`, and §3's own sentence is *"superconvergence
  is only obtained for `k ≥ 1`"*. The loss under the earlier method is asserted
  in Remark 2.1 and in §5; no table of it appears.
* **THE PAPER IS PARABOLIC AND HAS NO STEADY THEOREM.** Every estimate is
  `L^∞(0,T; L^2)` for eq (1)'s `∂_t u − Δu + F(u) = f`, and `Δt^{-1} M`
  contributes a positive-definite term to the (1,1) block that MEQ's steady
  problem does not have. There is no elliptic result in paper I, so
  transferring the rates to Grad–Shafranov is an inference, and MEQ's acceptance
  is a measurement rather than a check against a theorem. The paper's remaining
  assumptions bind too: the dual-regularity assumption before eq (10) is
  *"satisfied if Ω is convex"*, which MEQ's domains and its graded `Γ_h` are
  not, and Theorem 3.19's local-Lipschitz version assumes the mesh is
  quasi-uniform, which an adaptive MEQ mesh is not.
* **MFEM closes the local problem by row replacement, not by the paper's
  border.** `postprocess_hdg.cpp:204-205` zeroes row 0 of the `P^{k+1}` Neumann
  stiffness and writes the mass row `(χ_j, 1)_K` into it, then puts
  `(u_h, 1)_K` in that entry of the right-hand side (line 215). That is
  equivalent to (4) in exact arithmetic — the rows of the stiffness sum to zero
  for a partition-of-unity nodal basis and `Σ_i (q, ∇χ_i) = 0` makes the
  right-hand side consistent, so dropping any single row loses nothing — but it
  is a different `nd_s`-square system from the paper's `(nd_s+1)`-square
  bordered one, and the `B11`/`B12` that come out of it are that system's.
  MaNTA's 1-D implementation (`../MaNTA/Postprocessing.cpp`) takes the paper's
  bordered form with an Eigen `FullPivLU`. Either is fine; they must not be
  mixed up when the blocks are checked against each other.

Everything else in the brief is confirmed, including the two facts that decide
the placement: `DarcyForm::Reconstruct()` lifts the nonlinear potential
integrators as a Jacobian frozen at the computed potential
(`darcyform.cpp:1370-1404`), so its `ψ*` is not a linear function of the
unknowns; and `HDGPotentialPostprocessor`'s local matrix is
`AddMult_a_AAt(w, dshape_s, A)` and nothing else (`postprocess_hdg.cpp:164`),
with the physics entering only through the right-hand side, which is what makes
`B11` and `B12` constant.

## The seam

### On MFEM's side, and why it cannot be anywhere else

**1. `B11` and `B12` as matrices.** `HDGPotentialPostprocessor::Compute()`
already forms them and throws them away: the element loop builds `A`, replaces
the mean row, factors it, and back-solves per equation. Their extraction is
`Ai` applied to the two right-hand-side operators instead of to one vector —
`B11^ℓ = −Ã⁻¹ Ã2^ℓ` where `Ã2` is the gradient–flux coupling with row `i_c`
zeroed, and `B12^ℓ = (Ã⁻¹ e_{i_c}) (mass_p)^T`. What makes this the library's
business rather than MEQ's is not the algebra, which is thirty lines, but
everything the class already owns around it: reading the flux layout out of its
space rather than assuming it, the `ik`/`iK` diffusion-inverse coefficient — the
hook MEQ's `r` goes through — the `neq` blocking, and the enriched-space
construction. A MEQ copy would have to own all four, and the standing rule in
`CLAUDE.md` is a maintained library over a hand-rolled algorithm even when the
hand-rolled one measures no slower.

**2. A solution-dependent (1,0) gradient block.** `A9 diag(𝓕'(γ)) B11` is
`d_dofs × a_dofs` — **exactly `Bf_data`'s shape and orientation** — and it must
be added to `B` wherever `B` acts as the (1,0) block of the Jacobian. There are
three such sites and all are private:

| site | role |
|---|---|
| `ComputeElementH()`, `darcyhybridization.cpp:1446` | `S = D + B·AiBt`, and `BAiCt = B·AiCt` in the trace loop |
| `MultInv()`, `:3140` | the eliminated Jacobian applied, reached with `with_bnl = true` from `NPCReduce()`, `NPCRecover()` and `MultNL(GradMult)` |
| `LocalNLOperator::B`, `:4763`, used at `:5203` | `grad.SetBlock(1, 0, &B)`, hard-wired to the linear block |

The elimination itself needs **no** new algebra: it is a 2×2 block Schur
complement that already keeps the two off-diagonal blocks distinct — `AiBt`
carries the (0,1) block and `B` the (1,0) one — and the (0,1) side is already
solution-dependent through `Bnl_data`, whose header comment (`:322-342`) states
the mirror case in as many words: *"For a flux law `q = D(p) u` the flux
equation depends on the potential, so the (0,1) block of the local Jacobian is
not simply the transpose of the linear divergence form."* There is even a
precedent for the exact defect the omission causes: the comment at `:2181`
records that passing the linear `∓B^T` alone *"made the matrix-free gradient
disagree with the assembled one"*. The same funnel means one change to
`MultInv()` covers `GradientMode::MatrixFree` as well as `Assembled`.

**3. A slot that does not force `FullNL` — a question, not a demand.**
MEQ installs its source on `m_nlfi_p`, which selects
`LocalOpType::PotNL` (`:2463`), and `PotNL` is the branch where `InvertA()`
factors the flux block once and `ComputeElementH()` then skips re-factoring it
(`:1440`). Any `m_nlfi` at all falls through to `FullNL` (`:2485`), which gives
that up: one dense LU of `A` per element per Newton step. **Unmeasured**, and it
may be small. What MEQ would rather have is a slot that contributes the potential
residual and the (1,1) and (1,0) gradient blocks while leaving the flux block
linear. Whether that is expressible without a fourth `LocalOpType` is upstream's
call.

### On MEQ's side

* **`I_h` needs nothing new and MEQ should build it.**
  `L2_TriangleElement` is a `NodalFiniteElement`, so `I_h g` is `g` evaluated at
  `fe->GetNodes()` and nothing more. **But the basis matters and MEQ's present
  choice is wrong for this.** MEQ builds its spaces
  `L2_FECollection(k, dim, BasisType::GaussLobatto)` (`GradShafranov.cpp:533`),
  and `L2_TriangleElement`'s barycentric node construction with Lobatto points
  puts dofs on vertices and edges — MEQ's own comment says so. On a domain
  reaching the axis, which `examples/limited-tokamak.toml` and
  `tools/mesh/halfdisc.py` both are, an element vertex sits at `r = 0`, and the
  integrand MEQ interpolates is `F/r`. Quadrature never meets this because Gauss
  points are interior; **nodal interpolation does.** The fix is to build `Z_h`
  with an open basis — `BasisType::GaussLegendre` — whose triangle nodes are
  strictly interior. `u*` is the same `P^{k+1}` function either way, so this
  changes `I_h` and not the postprocessing, and `HDGPotentialPostprocessor`
  already accepts a caller-supplied space.
* **The axisymmetric weight and the `(r, z)` dependence are notational, as the
  roadmap says, but the choice of *where* to put `1/r` is not free.**
  Interpolating the whole integrand, `I_h(F(r,z,ψ*)/r)`, keeps `A9 = [(χ_j,
  φ_i)]` a pure mass coupling. Keeping `1/r` under the quadrature instead makes
  `A9 = [(χ_j/r, φ_i)]`, which diverges logarithmically on an element touching
  `r = 0` unless `χ_j` vanishes there. Take the first.
* **`ConfineToPlasma` is MEQ's, and the roadmap overstates it by one word.**
  `F` is discontinuous across `{Ψ = 0}` **at `j = 0` only**; at `j ≥ 1` it is
  `|d|^j`, a kink, and paper I's local-Lipschitz assumption is satisfied.
  `PLASMA-EDGE-PLAN.md` §1 measures MEQ's cut-element caps at
  `min(k+2, j+2.5)` for `ψ*`, so cut elements are already the accuracy
  bottleneck and an interpolant that is `O(h^{j+1})` there would make them worse
  rather than differently bad. The hybrid — interpolate on uncut elements, keep
  the quadrature where the edge cuts — is the answer, and it needs the same
  per-element mask `SourceIntegrator::setPlasmaComponent()` already carries.
  Nothing about it belongs in MFEM.

### The API MEQ would need, in the order it would be used

```cpp
// HDGPotentialPostprocessor, per element, in the space Compute() would use.
void GetLocalBlocks( int el, DenseMatrix &B11, DenseMatrix &B12 ) const;
void Assemble();                        // factor every element's A, once
void Compute( GridFunction &p_s ) const;   // as now, but an apply

// DarcyHybridization: the mirror of the existing (0,1) facility.
//   grad_arr( 1, 0 ) read rather than forced NULL, in ConstructGrad() and
//   LocalNLOperator::AddGradBlock(); stored per element in Bf_data's own
//   d_dofs x a_dofs orientation; added to B in the three sites above,
//   under the existing with_bnl flag.
```

Nothing in either is Grad–Shafranov-shaped.

## The case against, honestly

* **MaNTA already implements CCSZ-I without any of this**, in 1-D on Eigen:
  `../MaNTA/Postprocessing.cpp` stores `B11_`, `B12_`, `A9_` and `b1_` per cell
  from the paper's bordered system, `SystemSolver.cpp:1209-1264` accumulates
  both Jacobian blocks including the `B11` one, and `Tests/UnitTests/
  SolveJacTests.cpp:289-293` names the flux coupling as *"the only genuinely new
  coupling the scheme"* has. So the method demonstrably does not require library
  support — it requires being able to *reach* the block, which in MaNTA is
  MaNTA's own matrix and in MEQ is somebody else's private array. The argument
  for MFEM here is entirely about `DarcyHybridization`'s encapsulation, not
  about the difficulty of the algebra.
* **Swapping `ψ*` is not free, and it moves the estimator by construction.**
  `Estimator.hpp` computes all five η terms at `ψ*`, and η₂ is
  `‖q_h − (1/r)∇̄ψ*‖`. Under the classic postprocessing `∇u*` *is* the `L2`
  projection of `−iK q_h` onto `∇P^{k+1}(K)`, so η₂ becomes a measure of the
  non-gradient part of `q_h` — a defensible indicator, and a different function
  from what `Reconstruct()`'s richer `ψ*` gives. η₁, η₄ and η₅ move too.
  Every adaptive number in `CLAUDE.md` — `miller-adaptive`'s 4.73e-4 → 6.87e-5,
  the extension path's η₅ rate of 2.78 — would need re-measuring, and a change
  in any of them is a finding rather than a new baseline.
* **The second facility is not independently needed, and claiming otherwise
  would be wrong.** MEQ does not need `B11`/`B12` *as matrices* for anything it
  does today: `postProcess()` is called once per solve and the driver caches it.
  It becomes necessary the moment `ψ*` is wanted inside a residual evaluation —
  which is CCSZ-I, and which is also `ROADMAP.md` §12.2's plasma edge located on
  `ψ*`. Until one of those is built, the matrices are an optimisation of a call
  that happens once.
* **`FullNL` may cost more than the interpolation saves.** The saving is
  quadrature of `F` that the roadmap declines to apportion out of the 46–53% the
  device-offload plan gives all integrators; the cost is a per-element flux LU
  per step, plus one more element-local postprocessing solve per residual
  evaluation than MEQ does now. Neither side is measured. It is entirely
  possible that the honest answer after measurement is §12.3's Shamanskii
  iteration, which needs nothing from anybody.
* **A negative decision would have been available on a narrower reading.** If
  MEQ were willing to run the interpolatory source with the (1,0) block simply
  omitted, no library change would be needed at all, ever — see the next
  section. What that gives up is quadratic Newton, and MEQ's whole Newton
  chapter is about how little margin there is on stiff sources.

## What MEQ can do today, with no library change

The residual side needs nothing. `m_nlfi` is a
`BlockNonlinearFormIntegrator` that already receives `(u_l, p_l)` and writes to
both `bu` and `bp` (`AddMultBlock`, `:4846`), so a MEQ integrator can compute
`u*` locally, interpolate `F/r` onto `Z_h`'s nodes and write `A9 𝓕(γ)` into the
potential residual with the code as it stands. Under NPC the flux is Newton
state, so both arguments are the current iterate.

The Jacobian is then exact **except for the flux coupling**: the block
integrator can fill `grad_arr(1,1)` with `A9 diag(𝓕'(γ)) B12` — a
`d_dofs × d_dofs` block, which the existing plumbing accepts — and leave (1,0)
empty. That is a well-defined and useful intermediate:

1. **The discretisation's rates do not depend on the Jacobian.** A converged
   answer is the same whatever route reached it, so the `k+1` / `k+1` / `k+2`
   ladder on `ManufacturedNonlinear` and `SimilarityExponential` at `k = 1…3`
   can be established with (1,0) missing. That is the acceptance criterion, and
   it is reachable today.
2. **The cost of the missing block is then a measurement rather than a
   prediction** — iteration counts with and without, on the cases MEQ already
   has. That is the evidence a request should carry, and
   `HDG-BEM-COUPLING-FROM-MEQ.md`'s own rule is that *"a version that has been
   used is a better request"*.
3. `B11`/`B12` can be formed in MEQ meanwhile by duplicating
   `HDGPotentialPostprocessor`'s loop, which is short. If the blocks arrive from
   the library, MEQ's copy goes away; if the rates fail, MEQ deletes eighty
   lines and nobody upstream spent an afternoon.

Independently of all of it, and cheaper than any of it: **`ψ*` can move to
`HDGPotentialPostprocessor` now.** It is public, it takes MEQ's `r` through
`SetDiffusionInverse`, and the only `Reconstruct()` output anything reads is
`postProcessedPotential()` — `SolverContract.cpp:499` records that `trace()`,
`totalFlux()` and `postProcessedFlux()` *"are computed on every solve and
post-process and were checked by nothing"*. The classic postprocessing is also
the more robust of the two: `darcyform.cpp:1755` calls it *"immune to all of it,
structurally rather than by luck"* about the closure question, and it carries
none of `Reconstruct()`'s refusals — no `neq == 1` cap on an H(div) flux, no
requirement that the face constraint's gradient be trace-independent, no
dependence on the `gf-hdg-dev` reconstruction fix that
`thePostProcessedPotentialIsCorrectWhereTheJacobianVanishes` guards. The price
is the estimator re-measurement in *The case against*.

## Open questions

* **What fraction of an NPC step is quadrature of `F`?** Everything about
  whether this is worth building rests on it, and it is one `perf` run on
  `SourceIntegrator::AssembleElementVector` and `AssembleElementGrad` against
  the whole step. Do this before anything else.
* **Does it retire `setSourceQuadratureOrder()`?** The roadmap's second prize.
  Under interpolation the only quadrature of a nonpolynomial integrand left is
  the postprocessing's own, and `A9` is degree `2k+1`, inside the tabulated
  triangle rules — so the sweep into Grundmann–Möller territory
  (`QUADRATURE-HIGH-ORDER-TRIANGLES-FROM-MEQ.md`, weights to −1.9e+07) would
  have nothing to be swept for. Unverified: the cut-element hybrid keeps a
  quadrature, and it is exactly the case where the sweep mattered.
* **Is `HDGPotentialPostprocessor` reentrant enough to run inside
  `MultNL`'s threaded element loop?** `postprocess_hdg.cpp` takes
  `mesh->GetElementTransformation(z)`, the shared-scratch overload that
  `CLAUDE.md`'s *Traps* records as a silent wrong answer under threading. It is
  correct where it is used today, which is a serial post-processing call; a
  residual evaluation is not that. This is a requirement on the new path, not a
  defect in the existing one.
* **Which `ℓ` in Remark 2.2?** `ℓ = 0` is the paper's choice and gives a
  rank-one `B12`; `ℓ = k−1` matches `u*` to `u_h`'s lower moments and may
  behave better where `F` varies strongly inside an element, which for MEQ is
  the plasma edge. Unmeasured either way, and the analysis covers all of them.
* **Does `k = 0` matter to MEQ?** Paper I gives no superconvergence there and
  paper II's method (B) does, at the cost of different spaces and `τ = 1/h`.
  MEQ's convergence suite runs `k = 1…4`, so probably not — but the adaptive
  loop's coarse first cycle is where a `k = 0` question would come from.
* **Does the curved extension path survive it?** `Ω_h` is a union of background
  elements and the postprocessing is a plain element loop over them, so nothing
  obviously breaks; but `ψ*` through the transfer already reads `k+1` to `k+1.5`
  rather than `k+2` (`PLASMA-EDGE-PLAN.md` §1), and a method whose whole point
  is `ψ*`'s extra order meets that shortfall head on. Unmeasured.
