The papers behind meq's numerics. The PDFs are gitignored (publisher material,
several MB each); these tables are the tracked part, so fetch each one from its
doi and save it under the file name given.

**Every doi here has been checked against Crossref**, and the two arXiv ids
against the arXiv API. This is a rule, not a formality: a doi recalled from
memory is as likely to point at an unrelated paper as at the right one.

**Entries marked ✔ have been read against meq's formulation**; their annotations
say what the paper contributes and what it costs. Unmarked entries are seeds —
the annotation comes from the abstract or from how the read papers cite them, so
treat it as a reason to fetch the paper, not as a summary of what it says about
meq, and **replace it once you have read it**.

## The method meq implements

meq solves the **fixed-boundary** Grad–Shafranov equation: the plasma boundary
`Γ` is known, taken without loss of generality to be the level set `ψ = 0`, and
the problem is an interior Dirichlet one.

```
-∇̄·( (1/r) ∇̄ψ ) = F(r,z,ψ)/r   in Ω,        ψ = 0 on Γ
F(r,z,ψ) := μ₀ r² dp/dψ + g dg/dψ
```

recast as a first order system by introducing the **flux** `q = (1/r)∇̄ψ`, which
is not a numerical convenience — `q` is what the magnetic field is made of, so
discretising it directly is the point.

**These two papers are the method, not background to it.** Between them they fix
every choice in `src/meq`: the spaces, the numerical flux, the stabilisation, the
post-processing, the error estimator, and the treatment of the curved boundary.
Read HDG-GS-1 first and completely; HDG-GS-2 is the same method plus adaptivity
and is the one to work from thereafter.

**Both are open access on arXiv**, which is worth knowing because the versions of
record are paywalled at Elsevier, and both are in *Computer Physics
Communications* despite `pii` links that look like a different Elsevier journal.

The single most consequential number in both: **`τ = 1`**. The old meq code used
`τ = 5.0` with no recorded justification.

| Reference | URL (doi or arxiv) | Short Description | File Name |
| --- | --- | --- | --- |
| ✔ Computer Physics Communications 235 (2019) 120–132 | https://doi.org/10.1016/j.cpc.2018.09.013 | Sánchez-Vizuet & Solano, **the method meq implements**. The LDG-H formulation of the mixed system, eq. (8): spaces `V_h = [P_k]²`, `W_h = P_k`, `M_h = P_k(e)` — **all the same degree, with no inf-sup compatibility condition to respect**, which is the freedom hybridization buys. Numerical flux `q̂·n = q·n + τ(ψ_h − ψ̂_h)` with **`τ = 1`**; they note optimal order needs only `τ = O(1)`. Measured `k+1` for both `ψ` and `∇ψ` (Tables 1–5) on ITER, NSTX, ASDEX and FRC geometries, with round-off dominating beyond `k ≈ 5` or `h/8` — so a convergence study that stops improving there is behaving correctly, not failing. Two things to notice that are easy to miss: **they describe the superconvergent post-processing but did not implement it**, because the quantity of interest is `B ∝ ∇ψ`, which the method already gives at full order; and their §4 states the design decision meq reverses — keeping `F` as opaque problem data so the solver "relies only on the discretization of the toroidal operator `Δ*`", at the price of treating every source iteratively. Also carries the Solov'ev right-hand side `F = −((1−A)r² + A)` and the twelve homogeneous solutions. Paywalled; arXiv:1712.04148 | HDG-GradShafranov.pdf |
| ✔ Computer Physics Communications 255 (2020) 107239 | https://doi.org/10.1016/j.cpc.2020.107239 | Sánchez-Vizuet, Solano & Cerfon, the same method **plus adaptivity**, and the better statement of the discrete system — eqs (15)–(18) give the local/global block split explicitly, including which linear solves are shared between the two sides of the condensed system and therefore computed once. Introduces the **companion mesh** `T_c^h`, a minimal cover of `Ω̄` that exists only to update the computational domain during refinement; no computation is ever done on it. §2.7 gives a *nonlinear* local post-processing but then records that **dropping the `F` terms and using the plain linear post-processing is also effective, even in the semi-linear case** — which is what meq does, and this sentence is the licence for it. §3.1 is the residual estimator, and §3.3 the refinement strategy that keeps `dist(Γ_h, Γ) = O(h_loc)`; without that update the transfer degrades as the mesh refines, which is the counter-intuitive failure this paper exists to fix. **Read the warning at the end of §3.3**: refinement *appears* to crowd the boundary, but that is the proximity condition at work, not the error indicator, and reading it as the latter will send you chasing a bug that is not there. §4.1 is meq's headline benchmark — an exact NSTX Solov'ev solution with `A = −0.52` and twelve coefficients given to 15 digits. **Two errata found by checking it**: eq. (21) gives the Solov'ev source as `F = (1−A)r² + A`, which contradicts the paper's own eq. (1) — applying `Δ*` to the particular solution of eq. (22a) gives `(1−A)r² + A`, so `F = −Δ*ψ` carries a minus sign and HDG-GS-1 eq. (10) is the correct statement; and eq. (22c) prints `c₇` and `c₁₀` as the same value, `−0.000044132956899`, which may be a transcription slip. `CLAUDE.md` records the first in full. Paywalled; arXiv:1903.01724 | HDG-GradShafranov-Adaptive.pdf |
| Journal of Computational Physics 93 (1991) 1–107 | https://doi.org/10.1016/0021-9991(91)90074-u | Takeda & Tokuda, the standing review of MHD equilibrium computation for tokamaks. A hundred pages of what everyone else does and why; the reference for placing meq's approach among the alternatives rather than for anything meq implements. Paywalled | EquilibriumReview.pdf |
| SIAM Journal on Numerical Analysis 47 (2009) 1319–1365 | https://doi.org/10.1137/070706616 | Cockburn, Gopalakrishnan & Lazarov, the unified hybridization framework — where the trace equation and `τ` come from, and the reference for why the equal-order spaces above are admissible. Paywalled; free copy on PDXScholar | HDG-UnifiedHybridization.pdf |

## Curved boundaries by extension from polygonal subdomains

`Γ` is a smooth curve, and the plasma boundary of anything interesting is not a
polygon. The usual answers are to fit the mesh to `Γ` with isoparametric
elements, or to immerse `Γ` in a background mesh and accept low order. **Both
GS papers take a third route**, and it is the reason the mesh generator meq used
to carry has been deleted rather than ported: the computational domain `Ω^h` is
simply the union of background elements lying entirely inside `Ω`, and the
Dirichlet data is carried from `Γ` to `Γ_h` by a **line integral of the extended
flux along a transfer path**.

The consequence for meq is structural. A *uniform, shape-regular triangulation of
a box* is the whole meshing requirement — `Mesh::MakeCartesian2D` plus
refinement. There is no geometry-conforming mesh anywhere in the pipeline.

**This is already implemented in the MFEM branch meq builds on.**
`fem/darcy/extension_hdg.{hpp,cpp}` in `../mfem-hdg-dev` provides `TransferPath`,
`ClosestPointPath`, `LevelSetPath`, `VertexConePath`, `ElementExtension`,
`PathTraceCoefficient` and `HDGExtensionIntegrator`, with
`miniapps/hdg/extension.cpp` as a worked driver. So this section is here to
explain what that machinery is doing and why, not to be implemented from scratch.

| Reference | URL (doi or arxiv) | Short Description | File Name |
| --- | --- | --- | --- |
| SIAM J. Sci. Comput. 34 (2012) A497–A519 | https://doi.org/10.1137/100805200 | Cockburn & Solano, **the transfer-path technique itself** — Dirichlet data on the true boundary carried to the mesh boundary by a 1-D line integral of the extension, reducing a curved-boundary problem to quadrature along a segment. Cited as [28] in HDG-GS-1 and the origin of its eqs (5)–(7). The result that matters: optimal order is retained **even when `dist(Γ_h, Γ)` is only `O(h)`**, where earlier work needed `O(h^(k+1))`. Paywalled | HDG-CurvedExtensions.pdf |
| Mathematics of Computation 83 (2013) 665–699 | https://doi.org/10.1090/s0025-5718-2013-02747-0 | Cockburn, Qiu & Solano, the a priori analysis for the extension technique — HDG-GS-1's [34], and the source of its claim that the rates it measures "agree with the theory". Read this one before trusting an observed rate on a curved geometry, since it says which rates are actually predicted. Paywalled | HDG-ExtensionAnalysis.pdf |
| Journal of Scientific Computing 59 (2013) 512–543 | https://doi.org/10.1007/s10915-013-9776-y | Cockburn & Solano, the same for convection–diffusion. Listed because it is the nearest analysis to a Grad–Shafranov operator with a non-constant coefficient, `1/r` here. Paywalled | HDG-CurvedExtensionsConvDiff.pdf |

## Adaptivity, and the estimator meq's old code got wrong

The estimator is eq. (20) of HDG-GS-2, in five terms `η₁…η₅`: the strong residual
of the divergence equation, the residual of the constitutive equation
`q − (1/r)∇ψ`, the jumps of `q·n` and of `ψ` across faces, and the mismatch
between the hybrid unknown and the solution on element boundaries.

**Two of those terms use the post-processed `ψ*_h`, not `ψ_h`, and that is not
cosmetic.** The paper is explicit that `η₂` would converge at reduced order if it
used the raw `ψ_h`, because it differentiates the approximation; substituting
`ψ*_h` is what preserves order `k+1`. The estimator in meq's pre-modernisation
`CockburnEstimator.hpp` used raw `ψ_h` in both places — it was a degraded version
of the published estimator, and reproducing it would waste the post-processing
the solver already computes.

`η₅` likewise pairs the hybrid unknown against `ψ*_h`, where the old code paired
it against `ψ_h`.

| Reference | URL (doi or arxiv) | Short Description | File Name |
| --- | --- | --- | --- |
| Numerische Mathematik 148 (2021) 919–958 | https://doi.org/10.1007/s00211-021-01221-8 | Sánchez, Sánchez-Vizuet & Solano, a priori **and a posteriori** analysis of the unfitted HDG method for **semi-linear** elliptic problems. This is HDG-GS-2's reference [19] — the analysis it defers its estimator to as "the subject of a separate communication" — since published. **The load-bearing reference for anything meq does with the estimator**, because it is the only one that treats the semi-linearity and the unfitted domain together, which is exactly meq's situation. Open access | HDG-UnfittedSemilinear.pdf |
| Journal of Scientific Computing 90 (2022) | https://doi.org/10.1007/s10915-022-01767-1 | Sánchez, Sánchez-Vizuet & Solano, the sequel, for a broader class of nonlinear elliptic problems. Read after the Numer. Math. paper; relevant if the profile nonlinearity turns out to sit outside the earlier paper's hypotheses. Open access | HDG-UnfittedNonlinear.pdf |
| SIAM Journal on Numerical Analysis 51 (2013) 676–693 | https://doi.org/10.1137/120866269 | Cockburn & Zhang, the residual estimator HDG-GS-2's eq. (20) is built from, for **linear** elliptic problems on fitted polygonal meshes. The place to look for what each `η` term is estimating before the Grad–Shafranov complications are added. Paywalled | HDG-APosteriori.pdf |
| Mathematics of Computation 85 (2015) 1113–1141 | https://doi.org/10.1090/mcom/3014 | Cockburn, Nochetto & Zhang, the contraction property of adaptive HDG — the convergence proof for the solve→estimate→mark→refine loop, **assuming Dörfler marking**. Worth knowing because HDG-GS-2 actually used *maximum* marking with `γ = 0.3` in its experiments, for which the analysis is acknowledged as still open, so meq should be able to run both. Open access | AHDG-Contraction.pdf |
| Stenberg, ESAIM M2AN 25 (1991) 151–167 | https://doi.org/10.1051/m2an/1991250101511 | The local post-processing idea both GS papers use: solve a small problem per element in `P_(k+1)` and gain an order. The `ψ*_h` that `η₁`, `η₂` and `η₅` above depend on. Open access | Stenberg-Postprocessing.pdf |

## The nonlinearity — where meq leaves the papers

**meq uses Newton. Both papers use Anderson-accelerated Picard.** This is a
deliberate divergence and the entries below are here to make the trade explicit
rather than to be followed.

The papers' choice follows from their stated goal, quoted in §1 above: keep `F`
as problem data the user supplies, so the solver only ever discretises `Δ*`. The
price is that even a source *linear* in `ψ` — which could be folded into the
bilinear form as a mass matrix — is iterated. HDG-GS-1 §4 says so directly, of
the ASDEX example.

**Newton reverses that trade.** It puts `∂F/∂ψ` into the operator, which buys
quadratic convergence and costs the requirement that the source expose its
derivative. `src/meq/Source.hpp` therefore has `dFdPsi` alongside `F`, and
profile splines must supply `f'` as well as `f`. The obligation is real: a source
that cannot differentiate itself cannot be used.

Anderson is not thereby ruled out — MFEM's `KINSolver` offers it via
`EnableAndersonAcc`, and `KINSolver` derives from `mfem::NewtonSolver`, so the two
are interchangeable at the call site. Keep these papers to hand for the day a
pedestal makes Newton's globalisation an issue.

**And there is support for the choice from the free-boundary side**: CEDRES++
(§ *Free boundary* below) reports that fixed-point iterations converge very
slowly or fail outright, vertically unstable plasmas being the case that forced
the field to Newton-type methods. That is a stronger statement than "Newton is
faster", and it is the one to cite if the choice is ever questioned.

| Reference | URL (doi or arxiv) | Short Description | File Name |
| --- | --- | --- | --- |
| SIAM Journal on Numerical Analysis 53 (2015) 805–819 | https://doi.org/10.1137/130919398 | Toth & Kelley, convergence analysis for Anderson acceleration. HDG-GS-1's [45], and the source of its **depth `m = 2`**: they report no gain beyond `m ≥ 3`, and HDG-GS-1's own experiments agreed. The reference for what Anderson does and does not guarantee, which matters if it is ever used as a fallback when Newton stalls. Paywalled | AndersonAcceleration.pdf |
| SIAM Journal on Numerical Analysis 49 (2011) 1715–1735 | https://doi.org/10.1137/10078356x | Walker & Ni, Anderson acceleration for fixed-point iterations, and its relation to GMRES on the linearised problem. HDG-GS-2's [27]. Paywalled | AndersonWalkerNi.pdf |

## Analytic benchmarks

Every convergence claim in `tests/convergence` is measured against one of these.
The ordering matters: a wrong sign convention converges, at the right rate, to
the wrong function, so **an order-of-accuracy study alone cannot catch it** — only
a comparison against a closed form can. That is why the Solov'ev case with
published coefficients comes before any self-convergence study.

| Reference | URL (doi or arxiv) | Short Description | File Name |
| --- | --- | --- | --- |
| Physics of Plasmas 17 (2010) 032502 | https://doi.org/10.1063/1.3328818 | Cerfon & Freidberg, *"One size fits all"* analytic solutions to the Grad–Shafranov equation. **The parametrization behind both GS papers' test cases** — their `(ε, δ, κ)` for aspect ratio, triangularity and elongation, and the procedure that fixes the twelve coefficients `c₁…c₁₂` from geometric constraints. Cited as [47] in HDG-GS-1 and [41] in HDG-GS-2. Needed to generate a Solov'ev benchmark for any geometry other than the one HDG-GS-2 §4.1 tabulates. Paywalled | Soloviev-CerfonFreidberg.pdf |
| Physics of Plasmas 14 (2007) 112508 | https://doi.org/10.1063/1.2803759 | Guazzotto & Freidberg, a family of analytic equilibria with **dissimilar source functions** — `p = (S/μ₀)ψ` and `g² = Tψ² + 2Uψ + g₀²`, giving `F = Tψ + Sr² + U`, a source **linear in `ψ`** whose solution is the Bessel-function combination of HDG-GS-1 eq. (14). HDG-GS-1's [48] and its ASDEX Example 4. The useful benchmark for a source that depends on `ψ` at all, and the one where the papers' iterate-everything design is most obviously paying a price Newton does not. Paywalled | AnalyticEquilibria-Bessel.pdf |

## Free boundary: the direction after the fixed-boundary target

meq's first target is fixed boundary, and the free-boundary code it used to have
is unported in `attic/free-boundary/`. These two are the seeds for bringing it
back properly rather than restoring what was there.

The difficulty is that free boundary is an **exterior** problem: the vacuum field
extends to infinity, so the domain is unbounded and the plasma boundary is itself
an unknown, defined by a limiter contact or a separatrix. The standard answer,
and the one the old attic code implemented by hand, is Lackner's: reduce the
unbounded domain to a boundary condition through the Green's function, and couple
that to a finite element solve inside. Doing it as a *coupled BIM/FEM* rather
than as a hand-rolled Green's-function quadrature is the improvement to aim for.

**These are annotated from their abstracts and introductions only** — enough to
say why they are here, not enough to be read as summaries. None carries a ✔.

### What meq takes from Lackner, and what it does not

Lackner 1976 is the origin of the whole approach and is cited as such by
CEDRES++ and by the old attic code alike. **It is important to be precise about
which half of it meq is adopting, because the two halves have opposite verdicts.**

**Adopt: the reduction of the unbounded domain.** Representing the exterior
solution through the analytic Green's function, so the unbounded vacuum region
becomes a boundary condition on a finite computational domain, is sound, is
universal, and is what CEDRES++ means when it cites Lackner. Nothing here
disputes it.

**Reject: the outer fixed-point iteration it is embedded in.** Lackner's §2.1 is
a catalogue of iteration schemes for the nonlinearity, and read as a catalogue of
their failure modes it is damning:

* Plain Picard, his eq. (3) `Δ*ψ^(n+1) = −f(r, ψ^n)`, **"will converge to the
  physically trivial solution ψ ≡ 0 if admitted by the formulation of the
  problem."** Not slow — converged, and to the wrong answer.
* The Marder–Weitzner three-level scheme, eq. (5), fixes that by damping with a
  parameter `α`, and he is explicit about the cost: *"The price to be paid for
  this consists in a slow convergence rate, as α weights the correction to ψ^n at
  every three-step cycle."*
* The schemes that actually get used in practice keep chosen "conserved
  quantities" fixed cycle to cycle and pin the magnetic axis by adjusting
  vertical and radial field coefficients each iteration — an explicit analogue of
  feedback position control, needed because the physical configuration is
  vertically unstable.

That last point is the same one CEDRES++ makes when it says fixed-point
iterations "fail to converge — a very important example is vertically unstable
plasmas". Lackner is not contradicting it; he is describing the machinery a
fixed-point scheme needs in order to survive, and meq's answer is to not be a
fixed-point scheme.

**So the consequence for the boundary treatment is structural, not a preference.**
In Lackner's arrangement the Green's function enters as an *explicit boundary-
condition update evaluated inside the outer loop*, so its cost is paid every
cycle. That is exactly what `attic/free-boundary/` implements —
`GreensFunctionBoundaryCoefficient::Eval` calls `BoundaryPsi`, which loops over
every boundary face with singular quadrature, once per quadrature point, and the
whole `O(N²)` sweep is redone at each Picard iteration.

A coupled BIM/FEM formulation instead puts the boundary integral operators **into
the system**, where Newton differentiates through them and the interior and
exterior are solved together. Gatica & Hsiao is what makes that cheap: take the
artificial coupling boundary to be a circle, invert the boundary integral
operators exactly, and one weakly singular term survives. So the route is
Lackner's reduction, expressed as a monolithic coupled system, solved by Newton —
which is what CEDRES++ does and what meq should do.

**This reasoning is from §§1–2.2 of Lackner and the introduction of CEDRES++.**
The rest of both papers is unread, and the argument should be re-checked against
Lackner's §2.2 in particular, which is where the linearised solve in the
unbounded domain is actually described.

**One thing already worth extracting, because it bears on a decision meq has
already taken.** meq uses Newton where both Grad–Shafranov papers use
Anderson-accelerated Picard, and §4 above frames that as a trade. CEDRES++'s
introduction makes the stronger claim from the free-boundary side: *"Simple
fixed-point iterations usually suffer from very slow convergence or even fail to
converge, which made researchers move towards Newton-type methods... they can
converge in cases where fixed-point iterations do not converge — a very important
example is vertically unstable plasmas."* So Newton is not merely a preference
here; it is what the free-boundary literature moved to, and the case where Picard
fails outright is a physically important one. Confirm this against the body of
the paper before leaning on it.

| Reference | URL (doi or arxiv) | Short Description | File Name |
| --- | --- | --- | --- |
| Journal of Plasma Physics 81 (2015) 905810301 | https://doi.org/10.1017/s0022377814001251 | Heumann, Blum, Boulbe, Faugeras, Selig, Ané, Brémond, Grandgirard, Hertout & Nardon, **CEDRES++** — a survey of the computational methods in a production quasi-static *free-boundary* equilibrium code, on ITER and WEST. The nearest thing in the literature to what meq would become. Piecewise-linear FEM for the flux map coupled to a **boundary element method** for the unbounded exterior, and a **Newton method for the discretised nonlinear problem** covering all three nonlinearities at once — the current profile, the free plasma boundary, and the ferromagnetic permeability where there is an iron transformer. Also covers the *inverse* problems (find the coil currents giving a desired plasma shape), for which they say Newton is the main building block. Read for the coupling and the Newton formulation rather than for the discretisation, which is low order where meq's is not. Open access on Cambridge Core | CEDRES.pdf |
| Computer Physics Communications 12 (1976) 33–44 | https://doi.org/10.1016/0010-4655(76)90008-4 | Lackner, "Computation of ideal MHD equilibria" — **the origin of the Green's-function reduction of the unbounded domain**, and the method the code in `attic/free-boundary/` implements. Read the section above for which half of it meq adopts. Its §2.1 is the reason for the other half: plain Picard "will converge to the physically trivial solution ψ ≡ 0 if admitted by the formulation of the problem", the damped three-level alternative pays for stability with an explicitly slow convergence rate, and the schemes used in practice need conserved quantities and per-cycle magnetic-axis pinning to hold the column in place. Also worth knowing for the inverse problem: he classes formulation III — find the applied currents producing a prescribed plasma surface — as **badly posed in the sense of Hadamard**, which is the caveat to carry into any inverse work later. Paywalled | LacknerFreeBoundary.pdf |
| Journal of Mathematical Analysis and Applications 189 (1995) 442–461 | https://doi.org/10.1006/jmaa.1995.1029 | Gatica & Hsiao, the **uncoupling** of boundary integral and finite element methods for *nonlinear* boundary value problems. The trick: choose the artificial coupling boundary to be a **circle** (or a sphere in 3-D), which lets the boundary integral operators be inverted *exactly*, so the weak formulation retains only one boundary term — the weakly singular single-layer operator. They report the coding and computational work more than halved against standard coupling, and the quadrature made much easier because what survives is only weakly singular. Their model problem is exactly the shape of the free-boundary vacuum region: a nonlinear second-order elliptic equation inside, becoming Laplace in the unbounded exterior. Note the authorship: Gatica is at Universidad de Concepción, the same department as Solano of the two Grad–Shafranov papers. Paywalled | DecouplingBIM-FEM.pdf |

## The finite element library

| Reference | URL (doi or arxiv) | Short Description | File Name |
| --- | --- | --- | --- |
| Computers & Mathematics with Applications 81 (2021) 42–74 | https://doi.org/10.1016/j.camwa.2020.06.009 | Anderson et al., "MFEM: A modular finite element methods library" — the library meq is built on, and the citation to use. Note that meq does **not** build against `mfem/master`: it needs the HDG work in `../mfem-hdg-dev`, whose `DarcyForm` replaced the older `HDGBilinearForm` API meq used to depend on. `CLAUDE.md` says which branch and why. Open access | MFEM.pdf |
