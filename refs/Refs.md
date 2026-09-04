The papers behind MEQ's numerics. The PDFs are gitignored (publisher material,
several MB each); these tables are the tracked part, so fetch each one from its
doi and save it under the file name given.

**Every doi here has been checked against Crossref**, and the two arXiv ids
against the arXiv API. This is a rule, not a formality: a doi recalled from
memory is as likely to point at an unrelated paper as at the right one.

**CHECK BY TITLE, NEVER BY VOLUME AND PAGE.** The Maschke & Perrin row below
carried a wrong doi for months and it had *passed* an earlier check, because the
volume and the pages agreed — they belong to the right issue, and the doi
pointed at a different article inside it. Volume and page agreeing is what a
wrong doi in a correctly transcribed citation looks like.

**AND THIS FILE'S OWN CLAIM ABOUT THE ARXIV IDS WAS FALSE UNTIL 2026-09-02.**
The sentence above asserted the check as a standing rule while `TODO` recorded,
correctly, that it had never been done. **Both ids have now been checked against
the arXiv API and both hold**: arXiv:1605.01538 is Kaltsas & Throumoulopoulos
with the title and the published doi this file gives it, and arXiv:2606.11821 is
real, is VEQ, and has no version of record. So the claim was true when it was
finally tested, which is luck rather than diligence — **a claim that a check was
performed is not a check**, and that is the same lesson as the volume-and-page
one directly above.

**Entries marked ✔ have been read against MEQ's formulation**; their annotations
say what the paper contributes and what it costs. Unmarked entries are seeds —
the annotation comes from the abstract or from how the read papers cite them, so
treat it as a reason to fetch the paper, not as a summary of what it says about
MEQ, and **replace it once you have read it**.

## The method MEQ implements

MEQ solves the **fixed-boundary** Grad–Shafranov equation: the plasma boundary
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

The single most consequential number in both: **`τ = 1`**. The old MEQ code used
`τ = 5.0` with no recorded justification.

| Reference | URL (doi or arxiv) | Short Description | File Name |
| --- | --- | --- | --- |
| ✔ Computer Physics Communications 235 (2019) 120–132 | https://doi.org/10.1016/j.cpc.2018.09.013 | Sánchez-Vizuet & Solano, **the method MEQ implements**. The LDG-H formulation of the mixed system, eq. (8): spaces `V_h = [P_k]²`, `W_h = P_k`, `M_h = P_k(e)` — **all the same degree, with no inf-sup compatibility condition to respect**, which is the freedom hybridization buys. Numerical flux `q̂·n = q·n + τ(ψ_h − ψ̂_h)` with **`τ = 1`**; they note optimal order needs only `τ = O(1)`. Measured `k+1` for both `ψ` and `∇ψ` (Tables 1–5) on ITER, NSTX, ASDEX and FRC geometries, with round-off dominating beyond `k ≈ 5` or `h/8` — so a convergence study that stops improving there is behaving correctly, not failing. Two things to notice that are easy to miss: **they describe the superconvergent post-processing but did not implement it**, because the quantity of interest is `B ∝ ∇ψ`, which the method already gives at full order; and their §4 states the design decision MEQ reverses — keeping `F` as opaque problem data so the solver "relies only on the discretization of the toroidal operator `Δ*`", at the price of treating every source iteratively. Also carries the Solov'ev right-hand side `F = −((1−A)r² + A)` and the twelve homogeneous solutions. Paywalled; arXiv:1712.04148 | HDG-GradShafranov.pdf |
| ✔ Computer Physics Communications 255 (2020) 107239 | https://doi.org/10.1016/j.cpc.2020.107239 | Sánchez-Vizuet, Solano & Cerfon, the same method **plus adaptivity**, and the better statement of the discrete system — eqs (15)–(18) give the local/global block split explicitly, including which linear solves are shared between the two sides of the condensed system and therefore computed once. Introduces the **companion mesh** `T_c^h`, a minimal cover of `Ω̄` that exists only to update the computational domain during refinement; no computation is ever done on it. §2.7 gives a *nonlinear* local post-processing but then records that **dropping the `F` terms and using the plain linear post-processing is also effective, even in the semi-linear case** — which is what MEQ does, and this sentence is the licence for it. §3.1 is the residual estimator, and §3.3 the refinement strategy that keeps `dist(Γ_h, Γ) = O(h_loc)`; without that update the transfer degrades as the mesh refines, which is the counter-intuitive failure this paper exists to fix. **Read the warning at the end of §3.3**: refinement *appears* to crowd the boundary, but that is the proximity condition at work, not the error indicator, and reading it as the latter will send you chasing a bug that is not there. §4.1 is MEQ's headline benchmark — an exact NSTX Solov'ev solution with `A = −0.52` and twelve coefficients given to 15 digits. **Two errata found by checking it**: eq. (21) gives the Solov'ev source as `F = (1−A)r² + A`, which contradicts the paper's own eq. (1) — applying `Δ*` to the particular solution of eq. (22a) gives `(1−A)r² + A`, so `F = −Δ*ψ` carries a minus sign and HDG-GS-1 eq. (10) is the correct statement; and eq. (22c) prints `c₇` and `c₁₀` as the same value, `−0.000044132956899`, which is **almost certainly a typesetting duplicate**: solving Cerfon–Freidberg's ten NSTX single-null conditions with the other eleven fixed and `c₁₀` free gives `−2.87e-3`, 65× larger, and drops the residual of those conditions from 4.6e-2 to 8.9e-5. Corroborated independently — a Solov'ev normalisation should put the separatrix at `ψ = 0`, and this set does not: the X-point sits at `(0.6958, −1.8069)` with `ψ = −8.7e-3`, so the zero level set is not a closed curve. **Note that nothing in MEQ's suite can catch this**, since `ψ₁…ψ₁₂` are `Δ*`-harmonic and any `c_i` leaves both `Δ*ψ = −F` and every convergence rate exact; `CLAUDE.md` records the consequences. A third report, that `c₁`'s digits also disagree, is **false** — it came from `pdftotext`, which silently drops this paper's minus signs. Read the rendered page. `CLAUDE.md` records the first in full. Paywalled; arXiv:1903.01724 | HDG-GradShafranov-Adaptive.pdf |
| Journal of Computational Physics 93 (1991) 1–107 | https://doi.org/10.1016/0021-9991(91)90074-u | Takeda & Tokuda, the standing review of MHD equilibrium computation for tokamaks. A hundred pages of what everyone else does and why; the reference for placing MEQ's approach among the alternatives rather than for anything MEQ implements. Paywalled | EquilibriumReview.pdf |
| ✔ Journal of Computational Physics 228 (2009) 3232–3254 | https://doi.org/10.1016/j.jcp.2009.01.030 | Nguyen, Peraire & Cockburn, implicit high-order HDG for **linear** convection–diffusion. The baseline statement of the method MEQ's MFEM branch implements — `fem/darcy/`'s formulation is this one, and `../mfem-hdg-dev/doc/HDG-ROADMAP.md` names it as the target: fully discontinuous spaces for both flux and scalar, coupled only through a trace space and a stabilisation `τ`, with no Raviart–Thomas space anywhere. Read it for what `τ` is *for* before reading anything about what value to give it. Paywalled | HDG-NPC-1.pdf |
| ✔ Journal of Computational Physics 228 (2009) 8841–8855 | https://doi.org/10.1016/j.jcp.2009.08.030 | Nguyen, Peraire & Cockburn, the same for **nonlinear** convection–diffusion, and **the paper that says how a nonlinear HDG problem is meant to be solved**. §2.6 is the section to read: Newton–Raphson is applied to the **full** `(q, u, û)` system first, giving eq (14), a *linear* system in the increments `(δq, δu, δû)`; the hybridization is then applied to **that**, eqs (16)–(18), giving `K δΛ = F` for the trace increment alone. The local elimination is therefore a matrix inverse of block-diagonal `A`, `B`, `D` — "the inverse can be computed on each element independently … it results from applying the LDG method to solve the **linearized** PDE". **Every local operation is a linear solve; there is no element-local nonlinear iteration anywhere in the method.** MEQ does it the other way round — hybridize first, so the elimination is itself nonlinear and runs a Newton per element per residual evaluation — and `CLAUDE.md` records what that costs. **Linearise then condense, not condense then linearise.** Also the reference for a solution-dependent `τ(u_h, û_h)`, its positivity condition (6)–(7), and the `∂₁τ`, `∂₂τ` terms that make the Jacobian consistent — which is the same warning `fem/darcy/bilininteg_hdg.hpp` gives about `EvalGrad`. Paywalled | HDG-NPC-2.pdf |
| SIAM Journal on Numerical Analysis 47 (2009) 1319–1365 | https://doi.org/10.1137/070706616 | Cockburn, Gopalakrishnan & Lazarov, the unified hybridization framework — where the trace equation and `τ` come from, and the reference for why the equal-order spaces above are admissible. Paywalled; free copy on PDXScholar | HDG-UnifiedHybridization.pdf |

## Curved boundaries by extension from polygonal subdomains

`Γ` is a smooth curve, and the plasma boundary of anything interesting is not a
polygon. The usual answers are to fit the mesh to `Γ` with isoparametric
elements, or to immerse `Γ` in a background mesh and accept low order. **Both
GS papers take a third route**, and it is the reason the mesh generator MEQ used
to carry has been deleted rather than ported: the computational domain `Ω^h` is
simply the union of background elements lying entirely inside `Ω`, and the
Dirichlet data is carried from `Γ` to `Γ_h` by a **line integral of the extended
flux along a transfer path**.

The consequence for MEQ is structural. A *uniform, shape-regular triangulation of
a box* is the whole meshing requirement — `Mesh::MakeCartesian2D` plus
refinement. There is no geometry-conforming mesh anywhere in the pipeline.

**This is already implemented in the MFEM branch MEQ builds on.**
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

## Adaptivity, and the estimator MEQ's old code got wrong

The estimator is eq. (20) of HDG-GS-2, in five terms `η₁…η₅`: the strong residual
of the divergence equation, the residual of the constitutive equation
`q − (1/r)∇ψ`, the jumps of `q·n` and of `ψ` across faces, and the mismatch
between the hybrid unknown and the solution on element boundaries.

**Two of those terms use the post-processed `ψ*_h`, not `ψ_h`, and that is not
cosmetic.** The paper is explicit that `η₂` would converge at reduced order if it
used the raw `ψ_h`, because it differentiates the approximation; substituting
`ψ*_h` is what preserves order `k+1`. The estimator in MEQ's pre-modernisation
`CockburnEstimator.hpp` used raw `ψ_h` in both places — it was a degraded version
of the published estimator, and reproducing it would waste the post-processing
the solver already computes.

`η₅` likewise pairs the hybrid unknown against `ψ*_h`, where the old code paired
it against `ψ_h`.

| Reference | URL (doi or arxiv) | Short Description | File Name |
| --- | --- | --- | --- |
| Numerische Mathematik 148 (2021) 919–958 | https://doi.org/10.1007/s00211-021-01221-8 | Sánchez, Sánchez-Vizuet & Solano, a priori **and a posteriori** analysis of the unfitted HDG method for **semi-linear** elliptic problems. This is HDG-GS-2's reference [19] — the analysis it defers its estimator to as "the subject of a separate communication" — since published. **The load-bearing reference for anything MEQ does with the estimator**, because it is the only one that treats the semi-linearity and the unfitted domain together, which is exactly MEQ's situation. Open access | HDG-UnfittedSemilinear.pdf |
| Journal of Scientific Computing 90 (2022) | https://doi.org/10.1007/s10915-022-01767-1 | Sánchez, Sánchez-Vizuet & Solano, the sequel, for a broader class of nonlinear elliptic problems. Read after the Numer. Math. paper; relevant if the profile nonlinearity turns out to sit outside the earlier paper's hypotheses. Open access | HDG-UnfittedNonlinear.pdf |
| SIAM Journal on Numerical Analysis 51 (2013) 676–693 | https://doi.org/10.1137/120866269 | Cockburn & Zhang, the residual estimator HDG-GS-2's eq. (20) is built from, for **linear** elliptic problems on fitted polygonal meshes. The place to look for what each `η` term is estimating before the Grad–Shafranov complications are added. Paywalled | HDG-APosteriori.pdf |
| Mathematics of Computation 85 (2015) 1113–1141 | https://doi.org/10.1090/mcom/3014 | Cockburn, Nochetto & Zhang, the contraction property of adaptive HDG — the convergence proof for the solve→estimate→mark→refine loop, **assuming Dörfler marking**. Worth knowing because HDG-GS-2 actually used *maximum* marking with `γ = 0.3` in its experiments, for which the analysis is acknowledged as still open, so MEQ should be able to run both. Open access | AHDG-Contraction.pdf |
| Stenberg, ESAIM M2AN 25 (1991) 151–167 | https://doi.org/10.1051/m2an/1991250101511 | The local post-processing idea both GS papers use: solve a small problem per element in `P_(k+1)` and gain an order. The `ψ*_h` that `η₁`, `η₂` and `η₅` above depend on. Open access | Stenberg-Postprocessing.pdf |

## The nonlinearity — where MEQ leaves the papers

**MEQ uses Newton. Both papers use Anderson-accelerated Picard.** This is a
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

**But the comparison that actually matters is not Newton against Picard — it is
what the discretisation does to the nonlinear structure.** Serino et al. below
build a Newton solver *in MFEM* and report it robust where Picard fails, on a
harder (free-boundary) problem than MEQ's. The difference is that their `ψ` is a
single global `H¹` unknown, so `F(ψ)` enters one global residual and Newton
linearises once. MEQ's hybridization eliminates flux and potential element by
element, which turns a nonlinear `F` into a nonlinear solve **per element per
residual evaluation**. That is where MEQ's stiff cases fail, and it is why the
HDG-GS papers' Picard is a coherent choice rather than a conservative one:
Picard evaluates `F` at the previous iterate and leaves every local solve
linear. `CLAUDE.md` carries the measurements.

| Reference | URL (doi or arxiv) | Short Description | File Name |
| --- | --- | --- | --- |
| ✔ arXiv:2407.03499 (2025) | https://arxiv.org/abs/2407.03499 | Serino, Q. Tang, X.-Z. Tang, Kolev & Lipnikov, **an adaptive Newton-based free-boundary Grad–Shafranov solver, built on MFEM** — Kolev is an MFEM core developer. The reference that reframes MEQ's nonlinear difficulty. §3: `ψ ∈ H¹(Ω)`, standard shape functions, eq (3.1) one global weak form, eq (3.7) one linearised system per Newton step, with the Jacobian of the domain-dependent plasma term obtained by **shape calculus** rather than by differentiating a discrete residual — the opposite choice from MEQ's and CEDRES++'s, and worth reading against CEDRES++ §5, which rejects the analytic route because it "seems to blow up if ψ reaches a critical point". Two statements bear directly on MEQ: Newton succeeds where "conventional Picard-based solvers fail to converge" on a Taylor state, so Picard is not the safe default; and their own earlier fixed-boundary solver is described as "significantly easier" than this one, so MEQ cannot explain its difficulty by the boundary condition. Also the reference for preconditioning a GS Newton system — block factorisation with AMG on the elliptic sub-blocks | MFEM-GS-Newton.pdf |

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
| Physics of Plasmas 6 (1999) 3554–3560 | https://doi.org/10.1063/1.873630 | McCarthy, analytic solutions with **dissimilar source functions** — the pressure and the toroidal field function not proportional to one another. `p = (S/μ₀)ψ` and `g² = Tψ² + 2Uψ + g₀²`, giving `F = Tψ + Sr² + U`: a source **linear in `ψ`**, solved by the eighteen-term Bessel and Neumann combination of HDG-GS-1 eq. (14). This is HDG-GS-1's [48] and its ASDEX Upgrade Example 4. **Implemented as `tests/analytic/McCarthy.hpp`**, whose eighteen-term transcription is checked numerically: `−Δ*ψ` reproduces `F` to 3.6e-6 and the analytic gradients match finite differences to 5e-9. It is the middle rung of MEQ's benchmark ladder — `dF/dψ = T`, a nonzero *constant*, so a Newton Jacobian either carries the mass term or it does not and no algebra hides the difference, unlike Solov'ev where `dF/dψ` vanishes and the nonlinear manufactured case where it is messy. **Bibliographic trap**: an earlier version of this file credited eq. (14) to Guazzotto & Freidberg, *A family of analytic equilibrium solutions for the Grad–Shafranov equation*, Phys. Plasmas 14 (2007) 112508 — a real paper, but not this one and not HDG-GS-1's [48]. The giveaway was in MEQ's own code, where the fixture had been named `McCarthyEquilibrium` all along. Paywalled; **PDF not in refs/** | McCarthy-DissimilarSources.pdf |

## Toroidal rotation, and the generalised Grad–Shafranov equation

`FLOW-PLAN.md` is the design. The first of these **is the equation MEQ will
solve**; the other two are prior art on solving something close to it, and the
difference between "close" and "the same" is the whole hazard of this section.

**Rotating-equilibrium papers use at least three different closures** — each
species isothermal on a flux surface, adiabatic with a ratio of specific heats,
or constant density on a surface. All three produce
`Δ*ψ = −μ₀r² ∂p/∂ψ|_r − gg′` with `p` varying on the surface, and they are easy
to mistake for one another. Only the first is (136). Check any borrowed closed
form with `deltaStarFD()` before trusting it, exactly as the Solov'ev
coefficients are checked.

| Reference | URL (doi or arxiv) | Short Description | File Name |
| --- | --- | --- | --- |
| Rep. Prog. Phys. 76 (2013) 116201 | https://doi.org/10.1088/0034-4885/76/11/116201 | Abel, Plunk, Wang, Barnes, Cowley, Dorland & Schekochihin, *Multiscale gyrokinetics for rotating tokamak plasmas*. **The equation MEQ implements for flow**: eq (136), printed page 19, is the generalized Grad–Shafranov equation for a rigidly rotating plasma, closed by (96) for the poloidal density variation and (97) for the electrostatic potential `φ₀` that holds quasineutrality against it. Each species is a local Maxwellian with `T_s(ψ)`, so the closure is **isothermal on a flux surface** and the density is *not* a flux function. Three things to read beside it: (59), where `⟨φ₀⟩_ψ = 0` is stated to be an arbitrary convention rather than a closure — which is what lets MEQ use a local gauge instead; (128), the force balance, against which (136) reduces to `Δ*ψ = −4πR² ∂p/∂ψ\|_R − II′`; and §11.3 with (243), the low-Mach limit, which **is the static equation MEQ already solves** and is therefore a free regression test. Gaussian units, and `ψ` is poloidal flux per radian — the same convention as MEQ's `g g′`. Open access | RotatingGK.pdf |
| Comput. Phys. Commun. 260 (2021) 107264 | https://doi.org/10.1016/j.cpc.2020.107264 | Li & Zhu, *Solving the Grad–Shafranov equation using spectral elements for tokamak equilibrium with toroidal rotation* — the extension of NIMEQ, NIMROD's equilibrium solver. **The closest prior art to what MEQ is doing**: a high-order FEM Grad–Shafranov solver with toroidal rotation, and the source of the exact benchmark. Their (8), `P = P₀(ψ)exp[m_iΩ²R₀²/2T (R²/R₀² − 1)]` with `T = T_i + T_e`, **is** RoPP's isothermal closure with the gauge referenced to the magnetic axis rather than to `⟨R²⟩_ψ` — independent support for the local gauge. §3.1 gives a **new closed-form rotating Solov'ev solution**, (14)–(15), whose `M₀ → 0` limit (16) is the static Solov'ev particular solution; that is MEQ's acceptance test for FL-4. **Two cautions.** Its source is constant in `ψ`, so `∂F/∂ψ = 0` and it says nothing about the Jacobian. And **its eq (9) carries two sign errors** — the `dΩ/dψ` and `dT/dψ` corrections to `∂P/∂ψ` are both reversed relative to differentiating their own (8), found independently twice and confirmed to 40 digits. Harmless for both of their benchmarks, because each has `dC/dψ = 0`; live for anyone implementing (9) with profiled `T(ψ)` or `Ω(ψ)`. Their (6) also omits the `μ₀` their (9) carries, and their §5 parameters contradict their own Fig. 2 caption by a factor of 1.3e5 — use the text's values Also reports that the optimal range of their Picard relaxation parameter *narrows* with rotation, which is a hint the problem stiffens | SpectralElementGSRotation.pdf |
| Plasma Physics 22 (1980) 579–594 | https://doi.org/10.1088/0032-1028/22/6/007 | Maschke & Perrin, *Exact solutions of the stationary MHD equations for a rotating toroidal plasma*. **A second exact benchmark for (136), and the citation everyone gives for it is wrong** — Li & Zhu's [48] and `TODO` both name a Phys. Lett. A 102 (1984) 106 paper; the 1980 one is what is actually reproduced. **The 1984 paper is real, and its doi is `10.1016/0375-9601(84)90790-4`** — checked 2026-09-02 and recorded here so that the confusable one is pinned rather than looked up again by whoever meets the trap next. The two differ by a single word of title, *An analytic solution* (1984, Phys. Lett. A 102, 106–108) against *Exact solutions* (1980, Plasma Physics 22, 579–594), and the authors do not differ at all. That is the whole mechanism of the confusion. **And this row's own doi was wrong until 2026-09-01**: it read `.../22/6/009`, which Crossref resolves to *Propagation of guided electron plasma waves on a plasma cylinder* — a different article in the same issue. The right one ends `/007`. So the rule in this file's preamble caught a doi that had been checked against the *volume and page*, which agree, rather than against the *title*, which did not. It carries **two** solutions in two closures and only one is ours. **§4** takes the *temperature* as a surface quantity (`B·∇T = 0`) and is (136)'s isothermal closure; **§3** is a genuine polytrope `p = A(S)ρ^γ` and is **not**. The trap is that `γ` appears in §4 too — but only inside `γΩ²`, where its whole job is to make `Ω` the adiabatic Mach number, so it cancels out of the solution and **every `γ` is usable**, not just `γ = 1`. Verified by direct substitution into MEQ's equation set rather than by re-derivation: 8e-26 relative by 50-digit finite differences, exactly zero by sympy, with `g g′ ≠ 0` exercised. Restore `μ₀` by `p_SI = p_M&P/μ₀`, since they work in `j = ∇×B`. Two caveats: their (4.7) forces `C` constant, so this fixture **cannot see** the `C′(ψ)` term — the very term Li & Zhu got the sign of wrong — and `∂F/∂ψ = 0`, so it sits beside `Soloviev.hpp` rather than `McCarthy.hpp` | MaschkePerrin.pdf |
| Phys. Plasmas 11 (2004) 604–614 | https://doi.org/10.1063/1.1637918 | Guazzotto, Betti, Manickam & Kaye, *Numerical study of tokamak equilibria with arbitrary flow* — the FLOW code. **A different and harder problem than MEQ's**: poloidal *as well as* toroidal flow, where the equation changes type across the poloidal sonic surface and the free functions multiply. Kept as context and as a possible cross-check in a pure-toroidal isothermal limit, **not as a test of (136)**. Its value to MEQ is the survey of what an equilibrium-with-flow code has to get right, and the reminder that "equilibrium with flow" names a family rather than an equation. Paywalled | GuazzottoFLOW.pdf |

## Anisotropic pressure

`p_∥ ≠ p_⊥`, so the scalar `p(ψ)` in the source is replaced by a pair. **This is
the other generalisation, and it is independent of the rotation one above.**
`TODO` names the question that decides how large a job it is, and none of these
rows answers it yet: whether the anisotropy enters as a **modified source under
an unchanged `Δ*`**, in which case only `src/meq/Source.{hpp,cpp}` changes and
the discretisation is untouched — or whether it **changes the operator**, in
which case the HDG assembly changes too and this is a much larger piece of work.

**Every row here is an unread seed, and the provenance is the reason to be
careful.** They are transcribed from CEDRES++'s bibliography — the set a
production free-boundary code thought worth citing as future work, not a
literature search — which is exactly the provenance this file's opening rule
exists for. Checked against Crossref on 2026-09-02, that transcription turned
out to carry **one dead article number** and **one structural surprise**.

**THE LIST IS NOT SEVEN ANISOTROPIC-PRESSURE PAPERS.** It is a *mixed* list of
flow papers and anisotropy papers, because CEDRES++ names both as planned
extensions and the transcription copied the list whole. Sorted by what the
authoritative titles actually say:

| | |
|---|---|
| **anisotropy** | Zwingmann 2001, Cooper 2009, Pustovitov 2010 |
| **anisotropy *and* flow** | Fitzgerald 2013 |
| **flow only** | Goedbloed & Lifschitz 1997, Maschke & Perrin 1984, Guazzotto 2004 |
| **neither directly** | Grad 1967 — foundational, and the origin of the guiding-centre plasma Fitzgerald's title invokes |

Three of those seven are therefore about **rotation**, which is
`FLOW-PLAN.md`'s subject and not this one. `TODO` had already caught two of them
— it flags Maschke & Perrin as a wrong citation and Guazzotto as a different and
harder problem — so the pattern was half seen; the check is what shows it covers
Goedbloed & Lifschitz as well. Maschke & Perrin and Guazzotto stay in *Toroidal
rotation* above, where they belong, and are not duplicated here.

**So Fitzgerald et al is the one to read first**, on two independent counts: it
is the only *implementation* in the list, and it is the only entry that treats
anisotropy and toroidal rotation together — which, if it holds up, makes this
section and `FLOW-PLAN.md` one piece of work rather than two, exactly as `TODO`
suspected from the title. Note what the bibliographic check does and does not
establish: that the paper exists and says so, not that the physics composes.

| Reference | URL (doi or arxiv) | Short Description | File Name |
| --- | --- | --- | --- |
| Nuclear Fusion 53 (2013) 113040 | https://doi.org/10.1088/0029-5515/53/11/113040 | Fitzgerald, Appel & Hole, *EFIT tokamak equilibria with toroidal flow and anisotropic pressure using the two-temperature guiding-centre plasma*. **The place to start**, being the only implementation in the list and the only entry carrying both generalisations at once; the reference most likely to answer whether `Δ*` survives the anisotropy. **`TODO` transcribed the article number as 113408, and that doi is a 404** — the correct one is **113040**, a transposition of the last two digits, which is precisely the class of error a page-and-volume check waves through. Free to read at IOP, no open licence. Unread seed | Fitzgerald-AnisotropicEFIT.pdf |
| Plasma Physics and Controlled Fusion 43 (2001) 1441–1456 | https://doi.org/10.1088/0741-3335/43/11/302 | Zwingmann, Eriksson & Stubberfield, *Equilibrium analysis of tokamak discharges with anisotropic pressure*. The nearest of these to MEQ's own setting — axisymmetric and tokamak, equilibrium analysis rather than a 3D code — so the one to read beside Fitzgerald for what anisotropy does to an axisymmetric equilibrium. Free to read at IOP, no open licence. Unread seed | Zwingmann-AnisotropicEquilibrium.pdf |
| Plasma Physics and Controlled Fusion 52 (2010) 065001 | https://doi.org/10.1088/0741-3335/52/6/065001 | Pustovitov, *Anisotropic pressure effects on plasma equilibrium in toroidal systems*. Theory rather than an implementation, and therefore the one to read for what the anisotropy does to the equilibrium relations themselves before reading a code that has already committed to a closure. Free to read at IOP, no open licence. Unread seed | Pustovitov-AnisotropicPressure.pdf |
| Computer Physics Communications 180 (2009) 1524–1533 | https://doi.org/10.1016/j.cpc.2009.04.006 | Cooper, Hirshman, Merkel, Graves, Kisslinger, Wobig, Narushima, Okamura & Watanabe, *Three-dimensional anisotropic pressure free boundary equilibria*. **Read the first word of the title before planning around this one**: it is 3D and, from the author list, the VMEC/ANIMEC stellarator line — so it is evidence about the anisotropic closure and about what a free-boundary anisotropic solve has to get right, and **not** a template for an axisymmetric HDG discretisation. Hybrid open access at Elsevier. Unread seed | Cooper-Anisotropic3D.pdf |
| The Physics of Fluids 10 (1967) 137–154 | https://doi.org/10.1063/1.1761965 | Grad, *Toroidal Containment of a Plasma*. The foundational entry, and the one whose relevance is least direct: its abstract is about the nonexistence of flux surfaces, uncontained particle drifts and the nonexistence of self-consistent equilibria at small gyroradius, and does not mention anisotropy. It is here because Grad is the origin of the guiding-centre plasma that Fitzgerald's title invokes, so it is background to the **closure** rather than to the numerics. The journal registers under its pre-1994 name, *The Physics of Fluids*, which is why a search for "Physics of Fluids" can miss it. Paywalled. Unread seed | Grad-ToroidalContainment.pdf |
| Physics of Plasmas 4 (1997) 3544–3564 | https://doi.org/10.1063/1.872251 | Goedbloed & Lifschitz, *Stationary symmetric magnetohydrodynamic flows*. **A flow paper, not an anisotropy one**, and the title alone will not tell you: its abstract is about translation-symmetric self-similar flow, four velocity regimes separated by limiting-line and Alfvén singularities, and the three MHD shocks that connect them — anisotropic pressure does not appear in it. Kept in this section because that is where the transcription put it, flagged because it is the row most likely to be read as something it is not. Its nearest relative in `refs/` is `GuazzottoFLOW.pdf` and the same caution governs both: *flow* names a family of equations rather than one, and the difference between "close" and "the same" is this project's standing hazard. `TODO` abbreviates the title as "Stationary symmetric MHD flows"; the record spells out *magnetohydrodynamic*. Paywalled. Unread seed | GoedbloedLifschitz-StationaryFlows.pdf |

**Guazzotto, Betti, Manickam & Kaye 2004** is the seventh entry of the
transcribed list and is **already in this file**, under *Toroidal rotation*, as
`GuazzottoFLOW.pdf` — fetched and annotated. It is not repeated here.

## Free boundary: the direction after the fixed-boundary target

MEQ's first target is fixed boundary, and the free-boundary code it used to have
is unported in `attic/free-boundary/`. These two are the seeds for bringing it
back properly rather than restoring what was there.

The difficulty is that free boundary is an **exterior** problem: the vacuum field
extends to infinity, so the domain is unbounded and the plasma boundary is itself
an unknown, defined by a limiter contact or a separatrix. The standard answer,
and the one the old attic code implemented by hand, is Lackner's: reduce the
unbounded domain to a boundary condition through the Green's function, and couple
that to a finite element solve inside. Doing it as a *coupled BIM/FEM* rather
than as a hand-rolled Green's-function quadrature is the improvement to aim for.

### What MEQ takes from Lackner, and what it does not

Lackner 1976 is the origin of the whole approach and is cited as such by
CEDRES++ and by the old attic code alike. **It is important to be precise about
which half of it MEQ is adopting, because the two halves have opposite verdicts.**

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
fixed-point scheme needs in order to survive, and MEQ's answer is to not be a
fixed-point scheme.

**A correction to an earlier version of this section.** It said the attic code
implemented the wasteful branch of Lackner and that its `O(N²)` sweep was
intrinsic. That was wrong on both counts, and §2.2 is what settles it. Lackner
gives *three* ways to handle the unbounded domain, and ranks them:

1. Direct integration of the Green's function over the source region, his
   eq. (11) — `O(N²M²)`, and he calls it "the most straightforward, but also
   most wasteful".
2. Expansion of the source in Jacobi polynomials — cheaper, and it continues
   naturally outside the region so that torques on the coils fall out, but "for
   large aspect ratios, the expansion into Jacobi polynomials converges badly".
3. **Von Hagenow's method**, his eqs. (14)–(16), which he calls the most
   efficient: solve once with a fictitious conducting shell (`Δ*ψ̂ = −g`, `ψ̂ = 0`
   on `∂R`), recover the true boundary values from a **boundary** integral of the
   Green's function against `∂ψ̂/∂n`, then solve a second time with those values.
   Total cost of order **two fast solver steps**.

`attic/free-boundary/` implements the third, not the first — which the old
config file was telling us all along, with `BoundaryCondition = "VonHagenow"`,
and which `meq.cpp` performed as literally Lackner's three steps: solve with zero
BC, build `GreensFunctionBoundaryCoefficient` from the resulting flux, solve
again. So the *algorithm* was the good one.

What the attic code gets wrong is narrower and more ordinary: Lackner's cost
estimate assumes the `4(N+M)²` Green's function evaluations are done **once per
grid**, tabulated and reused. `BoundaryPsi` instead calls `GreensFunction(r, r*)`
inside the quadrature loop every time it is invoked, and it is invoked per
quadrature point per boundary face, per outer iteration. That is an
implementation defect, not a defect of the method.

**The real argument is still structural, but it is a different one, and Lackner
supplies it himself.** Of the capacitance-matrix variant he warns: *"If many
calculations have to be carried out for a given ∂R (e.g. if it corresponds to the
copper shell of an actually existing device), this algorithm will therefore be
quite efficient, but it will probably not be competitive with iteration methods
if the geometry of R is changed after each calculation."* The whole efficiency
argument rests on amortising a precomputation over a **fixed** boundary.

**MEQ's stage 6 is adaptive mesh refinement.** The geometry changes every cycle,
which is precisely the case Lackner excludes. So the amortisation that makes
von Hagenow cheap is the thing MEQ's adaptivity destroys — and that, rather than
any inefficiency in the method as published, is why MEQ should not build on it.

A coupled BIM/FEM formulation instead puts the boundary integral operators
**into the system**, where Newton differentiates through them and interior and
exterior are solved together, with no precomputation to invalidate. Gatica &
Hsiao is what makes that affordable.

**And this is confirmed rather than inferred.** CEDRES++ p. 13 says it outright:
*"The bilinear form c(·,·) follows basically from the so-called uncoupling
procedure (Gatica and Hsiao 1995) for the usual coupling of boundary integral and
finite element methods… The Green's function that is used in the derivation of
the boundary integral method for our problem was used earlier in finite
difference methods for the Grad–Shafranov–Schlüter equations (Lackner 1976)."*
Their `Γ` is a **semi-circle** of radius `ρ_Γ` enclosing the iron, coils and
passive structures — Gatica & Hsiao's circular coupling boundary, exactly. The
resulting form `c` is added to the operator in their variational formulation
(3.6), which is what "inside the system" means concretely.

### Four things CEDRES++ teaches about Newton, which apply to MEQ now

Not later, when free boundary arrives — now, while `Source` and the stage-4
Newton are being written.

**1. Differentiate the discrete residual, not the continuous one.** This is
CEDRES++'s stated distinguishing feature and the reasoning is a warning. The
*continuous* Newton derivative of the plasma-current term exists in the
literature (their eq. 3.25, from Blum 1989) and rests on shape calculus — and
they distrust it: *"there is no theoretical evidence that this formula holds also
for plasma equilibria with boundaries that contain X-points. In particular the
second term on the right-hand side seems to blow up if ψ reaches a critical
point."* Their §3.3 concludes such approaches "are not very trustworthy" and they
differentiate the Galerkin formulation instead. **MEQ gets this right by
construction** — `DarcyForm::GetGradient` differentiates the assembled operator —
but the reason is worth knowing, because it says the shortcut is wrong precisely
where the physics is interesting.

**2. Normalised flux makes the Jacobian non-local.** Their profiles are functions
of `ψ_N = (ψ − ψ_ax)/(ψ_bnd − ψ_ax)` on a fixed domain `[0,1]` — the same
convention as `meq::Profile`. But `ψ_ax` and `ψ_bnd` are *global functionals of
the solution*, so `∂F/∂ψ` picks up terms through them, and CEDRES++ p. 20 is
explicit: those terms *"lead to non-local entries in the stiffness matrix"*,
connecting the axis and boundary coefficients to every coefficient near the
plasma domain. **MEQ does not have this yet** — fixed boundary with `ψ = 0` on a
known `Γ` needs no normalisation — but the moment the profiles are driven by
normalised flux, `MHDSource::dFdPsi` as it stands is incomplete, and a
finite-difference check against `f()` alone will not reveal it, because both
would be missing the same terms.

**3. Higher order costs you quadrature derivatives.** Their §5 lists what stops
them going beyond `P1`, and one item lands squarely on MEQ: *"we need to
implement sufficiently accurate quadrature rules for polygonal domains with
non-straight boundaries. On top of this, we need to implement for the Newton
method the derivatives of such quadrature rules."* MEQ's stage 5 introduces
exactly such quadrature — the transfer-path integrals over the extension region —
and stage 4's Newton will have to differentiate through it. Expect that to be the
hard part.

**4. What "Newton is working" looks like.** Their Table 2, for a 577,415-unknown
ITER case: relative residual `2.67e0 → 9.16e-2 → 1.78e-3 → 5.25e-6 → 3.94e-12`
over five iterations. Textbook quadratic. That is the shape to expect from MEQ's
stage 4, and a run that grinds down linearly instead means a Jacobian that
disagrees with the residual — the failure this project's finite-difference test
on `dFdPsi` exists to catch.

**A footnote on their §5, which reads like a description of MEQ.** The obstacles
CEDRES++ lists to going higher order are: needing `hp` refinement because the
solution is non-regular near material interfaces; needing quadrature over
polygonal domains with non-straight boundaries; and, as a "promising
alternative", switching to *"separate meshes and separate polynomial degrees for
the representation of the flux in the plasma domain and its exterior"* so that
the magnetic axis and X-point are not confined to mesh vertices. That last is
close to what the extension-from-subdomains construction gives, and the first two
are what stages 5 and 6 are for. It is worth knowing that a production code's
open problems are this project's architecture, in both directions: it is
encouraging about the choice, and it is a warning about how much of the work is
in the quadrature.

**One thing already worth extracting, because it bears on a decision MEQ has
already taken.** MEQ uses Newton where both Grad–Shafranov papers use
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
| ✔ Journal of Plasma Physics 81 (2015) 905810301 | https://doi.org/10.1017/s0022377814001251 | Heumann, Blum, Boulbe, Faugeras, Selig, Ané, Brémond, Grandgirard, Hertout & Nardon, **CEDRES++** — a production quasi-static *free-boundary* code, and the nearest thing in the literature to what MEQ would become. P1 Lagrangian FEM on the interior coupled to a boundary integral form `c(·,·)` on a **semi-circular** artificial boundary, built by Gatica & Hsiao's uncoupling procedure and using Lackner's Green's function — see the two sections above, which are largely drawn from this paper. **Its distinguishing feature is a Newton method on the *discretised* equations**, chosen because the continuous-level derivative used by its predecessors SCED and Proteus has no theoretical backing at X-points and appears to blow up at critical points of `ψ`. Measured: perfect quadratic convergence, five iterations to 4e-12 on 577k unknowns (Table 2); linear convergence in the number of unknowns, as P1 implies. Also carries the four static/evolution × direct/inverse problem statements, the inverse problems as SQP with Tikhonov regularisation reusing *the same* derivatives as the Newton solve, and the flux-surface-average and geometric-coefficient post-processing in ITM conventions. Two cautions worth carrying: profiles in normalised flux make the Jacobian non-local, and their §5 names quadrature over polygonal domains with curved boundaries — plus **the derivatives of that quadrature** — as the obstacle to higher order. No analytic free-boundary solution exists, so they validate by convergence to a fine-mesh reference. Open access on Cambridge Core | CEDRES.pdf |
| ✔ Computer Physics Communications 12 (1976) 33–44 | https://doi.org/10.1016/0010-4655(76)90008-4 | Lackner, "Computation of ideal MHD equilibria" — a review, and **the origin of the Green's-function reduction of the unbounded domain**. §2.1 is iteration schemes for the nonlinearity and is the reason MEQ uses Newton: plain Picard "will converge to the physically trivial solution ψ ≡ 0 if admitted by the formulation of the problem"; the damped Marder–Weitzner three-level alternative buys stability with an explicitly slow convergence rate; and what was used in practice needed conserved quantities plus per-cycle magnetic-axis pinning, which is feedback position control in disguise. §2.2 is the unbounded-domain treatment and ranks three methods — direct Green's-function integration ("most straightforward, but also most wasteful"), a Jacobi-polynomial expansion that "converges badly" at large aspect ratio, and **von Hagenow's**, eqs (14)–(16), costing only two fast-solver steps. `attic/free-boundary/` implements the third. **The load-bearing sentence for MEQ is his caveat** that the precompute-and-reuse structure "will probably not be competitive with iteration methods if the geometry of R is changed after each calculation" — which is what adaptive refinement does every cycle. §2.3 is the inverse problem, classed as **badly posed in the sense of Hadamard**, with Fourier truncation and Zakharov's Tikhonov regularisation as the practical remedies, and a nice cautionary figure (his figs 3–4) showing that too many Fourier components leave the plasma surface unchanged while making the field near the conductors wild. Paywalled | LacknerFreeBoundary.pdf |
| ✔ SIAM Journal on Scientific Computing 34 (2012) A28–A47 | https://doi.org/10.1137/110823237 | Cockburn, Sayas & Solano, **Coupling at a Distance HDG and BEM** — the method MEQ's free boundary is built on, and the reason free boundary is approachable at all. An exterior Dirichlet problem solved by HDG on a **polyhedral** subdomain `D_h` coupled to a spectral method on a **smooth** artificial boundary `Γ` that the mesh does not fit, the two joined by Dirichlet-to-Neumann operators in the unmeshed region between. **Its reference [5] is Cockburn & Solano, which is MEQ's stage 5**: the family of paths, the extension `E_h(q_h)` and the lifting `L_h(g) = g∘a + ∫_σ C E_h(q_h)·m` are `mfem::TransferPath`, `ElementExtension` and `TransferredDatumCoefficient` term for term — so the free-boundary coupling is stage 5 with the datum `g` unknown instead of zero. Optimal orders measured with `dist(Γ_h, Γ) = O(h)`, which is assumption P.1 and which `meq::AdaptiveDomain` already maintains. §4.2 is the circular interface, where the boundary integrals are explicit — Gatica & Hsiao's uncoupling. **Two things MEQ does differently and deliberately**: their §4.3 couples by a Richardson fixed point between interior and exterior, which is outside the residual and is exactly what MEQ rejects everywhere else; and for the *axisymmetric* operator on a **semicircular** `Γ` the exterior map is diagonal in a Gegenbauer basis, so the BEM collapses to a symbol and there are no integral operators at all — derived and measured in `FREE-BOUNDARY-PLAN.md` §3, not in this paper. Paywalled | CouplingAtADistance.pdf |
| Journal of Mathematical Analysis and Applications 189 (1995) 442–461 | https://doi.org/10.1006/jmaa.1995.1029 | Gatica & Hsiao, the **uncoupling** of boundary integral and finite element methods for *nonlinear* boundary value problems. The trick: choose the artificial coupling boundary to be a **circle** (or a sphere in 3-D), which lets the boundary integral operators be inverted *exactly*, so the weak formulation retains only one boundary term — the weakly singular single-layer operator. They report the coding and computational work more than halved against standard coupling, and the quadrature made much easier because what survives is only weakly singular. Their model problem is exactly the shape of the free-boundary vacuum region: a nonlinear second-order elliptic equation inside, becoming Laplace in the unbounded exterior. Note the authorship: Gatica is at Universidad de Concepción, the same department as Solano of the two Grad–Shafranov papers. Paywalled | DecouplingBIM-FEM.pdf |

## Solution inversion: psi(R,z) to R(psi,l), z(psi,l)

**This is `ROADMAP.md` item 10's missing machinery**, and the same machinery
`MANTA-COUPLING.md` needs for both of its "genuinely hard" items --
flux-surface-averaged geometry and its psi-derivative. The prior art in this
tree is `v0-legacy:FluxSurfaces.cpp`, a predictor-corrector tracer that
**has never compiled** (`main()` has an uninitialised declaration and an
undeclared `psi_values`); the `Trace` function itself is sound and is textbook
Allgower-Georg.

**EVERY ENTRY BELOW IS A SEED.** The surveys they came from quoted them at
length, but none has been read against MEQ's own formulation to the standard
this file's preamble means by the tick, so none carries one. Replace the
annotation when you read it.

**Four findings the surveys established, recorded here because they decide the
design and are spread across several of the papers below.**

* **Subdivide-then-march -- what ParaView, VisIt, VTK and MFEM's exporters all
  do -- is second order in the SUB-CELL size with NO `k`-dependence at all.**
  Remacle et al.'s error reduces by `2^{2r}` for `r` levels of subdivision and
  the polynomial degree never enters. The naive route throws away `k-1` orders,
  which is the same structural gap `CLAUDE.md` records for `B` in the band,
  arriving from the other direction.
* **The "chords give second order" objection is NOT binding.** The
  representation error is governed by the point spacing `ds`, which is a free
  parameter decoupled from `h` and `k`. Chords need `ds <~ h^{(k+1)/2}`; cubic
  Hermite using the tangent from **`q`** needs only `ds <~ h^{(k+1)/4}` --
  about one point per element at `k = 3` -- and costs nothing, since the tracer
  already evaluates `q`.
* **For a CLOSED surface, a periodic parametrisation plus the equispaced
  trapezoidal rule converges GEOMETRICALLY** (Trefethen & Weideman). **The trap
  is the metric**: build a spectral rule and then difference neighbouring node
  positions for `|dx/dtheta|` and the whole scheme silently reverts to second
  order. `q` gives the arc-length element pointwise at full order. Not seen
  stated anywhere; this file is the record.
* **An exact rational contour exists only for `k = 2`.** A level set of a
  degree-`d` bivariate polynomial is an algebraic curve of degree `d`, and a
  rational parametrisation needs genus 0 -- a generic smooth plane cubic has
  genus 1. So Wiley et al.'s exactness is a `k = 2` technique in principle, not
  for want of effort, and MEQ runs to `k = 4`. Note also that MEQ's solve mesh
  is **straight-sided**, so their rational *quartic* (which comes from allowing
  curved element geometry) collapses to a conic here.

**THE DESC PAIR, READ 2026-09-03, AND THEY ANSWER THE `ψ`-ELEMENT QUESTION
MORE DIRECTLY THAN ANYTHING ELSE IN THIS LIST.**

| | |
|---|---|
| `DESC-Dudt.pdf` | Dudt & Kolemen, *DESC: A Stellarator Equilibrium Solver*, 2020. §III.B is the axis treatment; §III.A the Fourier-Zernike basis; (5) and (23) are `λ` and the gauge fixing |
| `DESC-Panici.pdf` | Panici et al., *The DESC stellarator code suite Part I*. Figure 5 is the spectral-width comparison against VMEC |

**Three things they establish, all of which MEQ needed.**

* **The axis regularity condition is Lewis & Bellan**, *J. Math. Phys.* **31**,
  2592-2596 (1990): a Fourier coefficient `a_m( rho )` of an analytic function on
  a disc must equal `rho^m` times an even series. `src/meq/Zernike.hpp` derived
  the same constraint independently and did not know the citation. Fourier-Zernike
  satisfies it inherently, which is the whole argument for the basis.
* **DESC reports MEQ's own axis finding**: *"more modes are required in the core,
  resulting in increased error near the magnetic axis for a given resolution"*,
  and attributes it to being *"restricted to operating in straight field-line
  coordinates, which may appear to be a disadvantage compared to the optimal
  poloidal angle of the VMEC formulation"* -- **Hirshman & Breslau**,
  spectral condensation.

  **BOTH ARE NOW IN `refs/` AND BOTH RECORDS ARE VERIFIED FROM THE PDFs,
  2026-09-04**, which is what this file's preamble means by checking a citation
  rather than inheriting one:

  | | |
  |---|---|
  | `LewisBellan.pdf` | *Physical constraints on the coefficients of Fourier expansions in cylindrical coordinates*, H. R. Lewis & P. M. Bellan, *J. Math. Phys.* **31**, 2592-2596 (1990), `10.1063/1.529009` |
  | `HirshmanBreslau.pdf` | *Explicit spectrally optimized Fourier series for nested magnetic surfaces*, S. P. Hirshman & J. Breslau, *Phys. Plasmas* **5**, 2664-2675 (1998), `10.1063/1.872954` |
  | `HirshmanVMEC.pdf` | *Steepest-descent moment method for three-dimensional magnetohydrodynamic equilibria*, S. P. Hirshman & J. C. Whitson, *Phys. Fluids* **26**, 3553-3568 (1983), `10.1063/1.864116` -- VMEC itself, and the origin of the optimal poloidal angle DESC compares itself against |

  Both titles taken from DESC's reference list survived the check unchanged, and
  so did the volumes and the pages. **What could NOT be verified is Panici**:
  `refs/DESC-Panici.pdf` is `arXiv:2203.17173v3` and carries no journal DOI at
  all, only a Zenodo one for the software. Its journal, volume and article
  number had been supplied from memory, so the bibliography now cites the
  preprint, which is what the document itself can be checked against.
* **And Panici measures the mechanism that fixes it**: *"DESC, while not
  explicitly enforcing any poloidal angle constraints, ends up finding an optimal
  representation through the course of the optimization procedure."* A solve
  slides the angle to what its basis can represent; a fit cannot. That is the
  structural difference between MEQ's IN-3 and DESC, and it is the argument for
  IN-4's gauge-free fit.

**AND THE GAP THAT IS MEQ'S TO FILL.** No paper in any of the three surveys
treats iso-contouring a field whose **gradient is an independently solved
unknown at the same order**. The nearest is the gradient-augmented level set
line (Nave, Rosales & Seibold), where the gradient is *advected* rather than
co-solved. Separately, nobody extracts from a `P_k` finite-element psi at
`O(h^{k+1})` and demonstrates that `q(psi)` and the flux-surface averages
converge at that rate. MEQ is the only code in either survey with `q` in hand.

| Reference | URL (doi or arxiv) | Short Description | File Name |
| --- | --- | --- | --- |
| Computer Physics Communications 190 (2015) 72–88 | https://doi.org/10.1016/j.cpc.2015.01.015 | Lee & Cerfon, **ECOM** — **the central reference for this item, and the only paper that treats extraction as a numerical-analysis problem in its own right.** Chebyshev-local Newton along radial rays (bracket the Chebyshev subinterval, build a local expansion, Newton, typically `m <= 5`), derivatives taken from the *same* Chebyshev representation so `dpsi/dr` on the contour is spectral too, then resample to equispaced `theta` and apply the trapezoidal rule — spectral **because** the integrand is periodic. Contours are placed on a **global Chebyshev grid in psi**, which is what protects the magnetic *shear* and not merely `q`. Their §IV Fig. 5 compares `q` and `s-hat` error against CHEASE **at equal run time**: CHEASE wins on `q` at small grids (the conformal map's crowding), shear is comparable, and ECOM wins outright as the grid refines. **Limitation that matters: the conformal map cannot handle an X-point.** | ECOM.pdf |
| Computer Physics Communications 97 (1996) 219–260 | https://doi.org/10.1016/0010-4655(96)00046-x | Lütjens, Bondeson & Sauter, **CHEASE** — the existing precedent for **order-preserving extraction from a finite-element psi**, which is exactly MEQ's situation. `ISOFIND` traces surfaces by intersecting them with the equilibrium mesh and **bisecting the element cubic** ("performance tests have shown that numerical evaluation of the roots by a bissection method is much cheaper ... than using the analytical Cardan formulas"), and the flux-surface integrals are done by a method that **"preserves the accuracy of the bicubic finite element solution"**. Also the reference for the star-shaped assumption being built into a *discretisation* rather than into an extractor: the polar mesh is re-centred on the magnetic axis after a coarse solve so that the innermost surfaces are guaranteed to surround the grid centre. §2.1 restricts to a **single magnetic axis**. | CHEASE.pdf |
| IEEE TVCG 15 (2009) 1227–1234 | https://doi.org/10.1109/tvcg.2009.194 | Etiene, Scheidegger, Nonato, Kirby & Silva — **where the order-of-accuracy numbers come from**, and the methodological match to this project. They apply the **method of manufactured solutions** to six isosurfacing codes and measure observed order against derived formal order: **vertex position `O(h^2)` (measured 1.94), normals `O(h)` (0.93), curvature does not converge at all (measured −3.35)**. Two things to carry: the geometric error carries a constant proportional to `1/\|grad f\|`, so it is **worst exactly at the magnetic axis and the X-point**; and their order study **found real bugs** (a hard-coded iteration cap in a root finder) that no visual inspection or unit test had caught — the same argument `CLAUDE.md` makes for rate tables over unit suites, reached independently. | Verifiable_Visualization_for_Isosurface_Extraction.pdf |
| Computer Methods in Applied Mechanics and Engineering 313 (2017) 759–784 | https://doi.org/10.1016/j.cma.2016.10.019 | Fries, Omerović, Schöllhammer & Steidl, **higher-order meshing of implicit geometries, part I** — the **rigorous `O(h^{k+1})` reconstruction of a curved zero level set**, from the unfitted-FEM literature rather than the visualization one, and **simplex-friendly** where Saye's hyperrectangle method is not. Reports optimal convergence on the automatically reconstructed surface meshes, and records that **node placement optimised for integration beats equally-spaced nodes** significantly. Part II (`10.1016/j.cma.2017.07.037`) is approximation and integration *on* the extracted manifold and is the half that bears on flux-surface averages. **The local copy is the arXiv preprint (1706.00578v1), not the version of record.** | HigherOrderMeshing.pdf |
| The Visual Computer 24 (2008) 187–200 | https://doi.org/10.1007/s00371-007-0184-x | Reuter, Mikkelsen, Sherbrooke, Maekawa & Patrikalakis — polynomial system solving in the **barycentric Bernstein basis**, so the subdivision primitive is a **triangle** rather than a bounding box. **MEQ's mesh is native to it**, with no box-splitting of triangular elements, plus degree elevation to tighten the convex hull. This is the route to finding the magnetic axis and any X-point **exhaustively**: the hull test discards only provably-empty cells, so nothing can be missed. Displaces the box-based literature (Sherbrooke & Patrikalakis' Projected Polyhedron, `10.1016/0167-8396(93)90019-y`; Mourrain & Pavone, `10.1016/j.jsc.2008.04.016`) for this application. **Pair it with a degree test, and note that topological degree gives the SUM of indices, never a count — degree zero does not imply no root, so it is a certification test and never an exclusion test.** | SolvingNonlinearPolynomialsBarycentric.pdf |
| Joint EUROGRAPHICS–IEEE TCVG Symposium on Visualization (VisSym 2003) 167–176 | https://doi.org/10.2312/VisSym/VisSym03/167-176 | Wiley, Childs, Gregorski, Hamann & Joy, *Contouring Curved Quadratic Elements*. **The exact-contour result, and its ceiling.** "A contour in a bivariate quadratic function defined over a triangle in parameter space is a conic section and can be represented by a rational-quadratic function, while in physical space it is a rational quartic." The quartic arises only from allowing **curved element geometry**; MEQ's solve mesh is straight-sided (`curveBoundaryOnto` is output-only and runs last), so the physical-space contour of a `k = 2` field is itself an exact conic. **It does not extend past `k = 2`** — see the genus argument above. Predecessor: Worsey & Farin, `10.1016/0167-8396(90)90041-o`. **Not in Crossref** — Eurographics registers its own handles; the doi resolves to `diglib.eg.org` and encodes the page range. | Contouring_Curved_Quadratic_Elements.pdf |
| Dagstuhl Follow-Ups 2 (2011) 276–291 | https://doi.org/10.4230/DFU.Vol2.SciViz.2011.276 | Pagot, Vollrath, Sadlo, Weiskopf, Ertl & Comba, *Interactive Isocontouring of High-Order Surfaces* — **the best short survey of the direct-extraction field**; its §2 is the map. Also the honest caveat on interval arithmetic for FEM data ("polynomials with hundreds of coefficients and thousands of polynomials ... for each frame") and on point-based methods ("every isovalue change typically takes several minutes"). Its Fig. 1 is a **DG** simulation. **Not in Crossref** — Dagstuhl registers with DataCite; the doi resolves to `drops.dagstuhl.de`. | IsocontouringHOSurfaces.pdf |
| IEEE TVCG 12 (2006) 114–125 | https://doi.org/10.1109/tvcg.2006.12 | Nelson & Kirby — ray-tracing spectral/hp elements for isosurfaces, the origin of the **ElVis** line. Restrict the field to a ray **within one element**, giving a univariate polynomial, and root-find; the element partitioning is what respects the DG jump, since roots are sought segment by segment and never interpolated across a face. Later members: `10.1109/tvcg.2011.206` (cut surfaces on the GPU), `10.1109/tvcg.2012.218` (ElVis itself), `10.1109/tvcg.2013.96` (volume rendering). Note the 2011 paper's contouring **assumes `C^0` across faces**, which is the assumption MEQ's field does not satisfy. | ray-tracing-polymorphic-multidomain-spectral-hp-elements.pdf |
| Computer Graphics Forum 21 (2002) 19–32 | https://doi.org/10.1111/1467-8659.00563 | Theisel, *Exact Isosurfaces for Marching Cubes* — trilinear contours as trimmed rational cubic Bézier patches, i.e. the exact level set of the interpolant marching cubes only approximates. Relevant as the tensor-cell analogue of Wiley et al., and as the demonstration that "exact" is available for low-degree interpolants and not beyond. *Crossref gives pp. 19–32; the PDF's own header says 19–31.* | IsosurfacesForMarchingCubes.pdf |
| SIAM Journal on Scientific Computing 35 (2013) A212–A230 | https://doi.org/10.1137/120874059 | Mirzaee, King, Ryan & Kirby — **SIAC filtering on unstructured triangular meshes**, which is the principled answer to "the contour does not close up across a face". Convolving a B-spline kernel against a DG field raises inter-element regularity to `C^{k-1}` and **increases** the order of accuracy rather than smearing it, so you contour the filtered field. **The unstructured-triangle case is the one that matters for MEQ**; the uniform-mesh superconvergence result is the classical one. | SmoothnessIncreasingFilters.pdf |
| Journal of Scientific Computing 38 (2009) 164–184 | https://doi.org/10.1007/s10915-008-9230-8 | Walfisch, Ryan, Kirby & Haimes — the **one-sided** SIAC variant, which is what makes the filter usable near a domain boundary where the symmetric kernel's support runs off the edge. Directly relevant on the extension path, where `Gamma_h` is where the interesting surfaces are. Companion: Steffen, Curtis, Kirby & Ryan, `10.1109/tvcg.2008.9`, the streamline-integration original. | OneSidedSmoothnessFilters.pdf |
| Journal of Computational Physics 375 (2018) 1179–1204 | https://doi.org/10.1016/j.jcp.2018.09.017 | Maunoury, Besse, Mouysset, Pernet & Haas — builds a **representation mesh adapted to the `hp` solution** rather than tessellating uniformly, which is the middle road between subdivide-then-march and true high-order extraction. Successor: Feuillet, Maunoury & Loseille, `10.1016/j.jcp.2020.109860`, "pixel-exact" rendering of high-order meshes *and* solutions. | PostProcessingForVis.pdf |

**FOUR OF THESE ARRIVED 2026-09-02 AND ARE NOW ON DISK**; the file name column
says which. Allgower & Georg is on disk only in part and outside `refs/` — see
its row.

| Reference | URL (doi or arxiv) | Short Description | File Name |
| --- | --- | --- | --- |
| SIAM Classics in Applied Mathematics 45 (2003) | https://doi.org/10.1137/1.9780898719154 | Allgower & Georg, *Introduction to Numerical Continuation Methods* — **the textbook for what `v0-legacy:FluxSurfaces.cpp` already does.** Ch. 3 Newton as corrector (17–27); ch. 5 convergence of Euler-Newton methods (37–43); **ch. 6 steplength adaptation (44–60), seventeen pages, the longest topic in the book**; ch. 8 detection of bifurcation points (75–90), which is what an X-point is in this language; and **ch. 15, "Approximating Implicitly Defined Manifolds" (233–251)**, the under-cited one and the closest to covering a whole family of surfaces rather than one path. The load-bearing fact: **the predictor's order buys step length, not accuracy** — the corrector puts every accepted point on the level set to solver tolerance regardless, which is why Euler-plus-Newton is the norm. Original: Springer SCM 13 (1990), `10.1007/978-3-642-61257-2`. **The free PDF that used to sit on the author's page is now a SOFT 404 — HTTP 200 with an HTML body — so a status-code check passes and hands you a 274-byte "PDF".** **THE WHOLE BOOK IS ON DISK at `~/AllgowerBook/1.9780898719154.ch{1..16}.pdf`**, outside `refs/` because it is per-chapter rather than one document. **Ch. 15 is the one to read before choosing a representation**, and not for the reason its title suggests: with `H(R, z, c) = psi(R,z) − c` the entire family of flux surfaces is ONE implicitly defined 2-manifold, and since `dH/dc = −1` everywhere that manifold has **no critical points at all** — it is smooth over the axis and over an X-point, where the individual level sets are not curves. See `INVERSION-PLAN.md` §3.1. | `~/AllgowerBook/` (ch. 1–16) |
| SIAM Review 56 (2014) 385–458 | https://doi.org/10.1137/130932132 | Trefethen & Weideman, *The Exponentially Convergent Trapezoidal Rule* — why a closed flux surface should be integrated in a **periodic** parameter, and the justification for the whole Fourier-in-angle approach. Geometric rather than algebraic convergence for a smooth periodic integrand. | TrefethenTrapezoidal.pdf |
| Computational Geometry 24 (2003) 75–94 | https://doi.org/10.1016/s0925-7721(02)00093-7 | Carr, Snoeyink & Axen, *Computing contour trees in all dimensions* — **the canonical contour-tree construction**: join tree and split tree by one union-find sweep each, then merged. `O(n log n + N alpha(N))`, the `n log n` being a global sort of vertex heights. Stated for simplicial complexes in any dimension, so unstructured triangles are native. **For MEQ only the JOIN tree is needed** — "nested closed surfaces bounded by the separatrix" is purely a superlevel-set question, which halves the work. The output-sensitive version, Chiang, Lenz, Lu & Rote, `10.1016/j.comgeo.2004.05.002`, sorts only the *critical* vertices at `O(N + t log t)`, which for `t ~ 3` removes the sort entirely. | ComputingContourTrees.pdf |
| IEEE TVCG 26 (2020) 162–172 | https://doi.org/10.1109/tvcg.2019.2934338 | Jallepalli, Levine & Kirby — **the DG obstruction, measured.** Evaluates the three standard ways of getting a high-order FEM field into a topological analysis (regular-grid sampling, element subdivision, L-SIAC) and finds all lossy, all distorting of feature shape, and — counterintuitively — **sampling and subdivision can perform WORSE as sampling resolution increases**. Read before building anything topological on `psi_h`. **MEQ has an asset this literature does not assume: HDG carries a single-valued trace `psihat` on every face, which is a ready-made gluing rule — topology consistent by construction rather than by filtering.** That last sentence is this project's inference, not the paper's. | The_Effect_of_Data_Transformations_on_Scalar_Field_Topological_Analysis_of_High-Order_FEM_Solutions.pdf |
| Computer Graphics Forum 36 (2017) 23–33 | https://doi.org/10.1111/cgf.13165 | Nucha, Bonneau, Hahmann & Natarajan, *Computing Contour Trees for 2D Piecewise Polynomial Functions* — **the only paper on this exact problem**, per-cell monotone path tracing then stitched. A Crossref search for follow-up work by that group returns none. | ComputingContourTrees2D.pdf |
| Discrete & Computational Geometry 37 (2007) 103–120 | https://doi.org/10.1007/s00454-006-1276-5 | Cohen-Steiner, Edelsbrunner & Harer, *Stability of Persistence Diagrams* — the reason persistence is the right tool for spurious critical points **and** the reason its threshold need not be tuned: the diagram is 1-Lipschitz in the sup norm, so **every spurious feature has persistence at most `2\|\|psi_h - psi\|\|_inf`**, a quantity MEQ already measures to convergence-rate precision. A threshold derived from the discretisation error rather than guessed. Foundational paper: Edelsbrunner, Letscher & Zomorodian, `10.1007/s00454-002-2885-2`. | — |

**A scoping result that belongs with these, established 2026-09-02 and not from any
single paper.** For MEQ's **fixed-boundary** problem the whole open-versus-closed
and separatrix-crossing apparatus is **not needed**. Poincare-Hopf on a disc gives
`#max + #min - #saddle = chi = 1`, and the *maximum principle* — not topology —
supplies `#min = 0`, since single-signed `F` makes `psi` a supersolution of a
uniformly elliptic operator so its minimum sits on `Gamma`. Hence one axis implies
**zero interior saddles**. **Free boundary is where this changes, and it changes by
force**: the vacuum region makes the domain an **annulus**, `chi = 0`, so
`#max + #min - #saddle = 0` and with one axis a saddle **must** exist. X-points are
forced by the topology of the domain rather than being incidental. Two carve-outs:
GS-2 §4.4, whose reversed core current destroys the supersolution property and is
the one fixed-boundary case where an interior saddle is possible; and rotating
runs where the object of interest is `p` rather than `psi`, since `p` is not a flux
function and has its own critical points in its own places.

**And a trap worth one line: the Morse-Smale complex gives the WRONG curves.** Its
"separatrices" are integral curves of `grad psi`, **orthogonal** to the level sets;
the plasma separatrix is the level **set** through the X-point. Reaching for an MS
complex expecting one yields the perpendicular family.

## Other analytic benchmarks, and what each one actually tests

MEQ's benchmark ladder is a ladder in how the source depends on `ψ`, because
that is what decides what a Newton Jacobian has to get right: `Soloviev.hpp`
(constant, `∂F/∂ψ = 0`), `McCarthy.hpp` (linear, `∂F/∂ψ = T`),
`ManufacturedNonlinear.hpp` and `SimilarityExponential.hpp` (nonlinear).
`CLAUDE.md` has the table.

**The distinction worth keeping in mind is manufactured versus exact.** A
manufactured solution picks a convenient `ψ` and *builds* `F` to fit it, so the
resulting source is not of a form any physical profile produces — HDG-GS-1's
Example 5 is one of these. An exact solution picks the free functions `p(ψ)` and
`g(ψ)` and *derives* `ψ`. Both are valid verification, but only the second tests
the solver against an equation somebody might actually pose.

**The best single find of a literature search for nonlinear cases was that four
of them were already on disk** — `HDG-GradShafranov-Adaptive.pdf` §§4.2–4.5 give
a pressure pedestal, an internal transport barrier, a current hole and an
internal layer, every parameter specified, all differentiating in closed form.
They have no exact solutions, so the paper judges them by its error estimator;
but the pedestal at `σ² = 0.005` is much the stiffest thing available and is what
the Newton path has not yet been pushed on. `TODO` has the equations transcribed.

| Reference | URL (doi or arxiv) | Short Description | File Name |
| --- | --- | --- | --- |
| ✔ Physics Letters A 380 (2016) 3373–3377 | https://doi.org/10.1016/j.physleta.2016.08.011 | Kaltsas & Throumoulopoulos, exact solutions by **similarity reduction**, and the source of `tests/analytic/SimilarityExponential.hpp`. Their eq. (1) is `Δ*u + f(u) + g(u)r² = 0`, so MEQ's `F = f(ψ)(1 + εr²)`; the reduction needs `g = εf` with `ε = a²/b²`, and `u = ψ` when there is no poloidal flow. §3.2 gives two nonlinear cases: `f = f₀ + f₁u²` solved by a **Weierstrass elliptic function**, and `f = f₀e^{nu}` solved by **sech²** (eq. 22) — the latter being, as far as this search found, **the only exact solution in elementary functions of a genuinely nonlinear Grad–Shafranov equation**. Verified here: the sech² form satisfies their eq. (11) identically by hand, and `−Δ*ψ = F` to 1.7e-7 numerically. **Read the rendered page, not an extraction**: `pdftotext` renders `g(u) = εf(u)` as `g(u) = f(u)`, dropping the `ε`, the same failure mode that ate every minus sign in the Solov'ev coefficients. One caveat the paper states itself — its nonlinear solutions "do not form closed surfaces", so there is no `ψ = 0` plasma boundary; irrelevant for a fixed-boundary benchmark on a rectangle with exact Dirichlet data, which is how MEQ uses it. Open access, arXiv:1605.01538 | GS-SimilarityReduction.pdf |
| ✔ Physics of Plasmas 17 (2010) 032502 | https://doi.org/10.1063/1.3328818 | Cerfon & Freidberg, *"One size fits all"*. §II gives the seven-term up-down **symmetric** expansion; **§IX gives the twelve-term asymmetric one** and, in its eq. (28), the twelve constraints that determine `c₁ … c₁₂` for a single null. That section is what showed the coefficients printed in HDG-GS-2 eq. (22c) to be wrong throughout — see `CLAUDE.md`. Also carries the `N₁, N₂, N₃` curvature coefficients, the X-point placement `x_sep = 1 − 1.1δε`, `y_sep = −1.1κε`, the figures of merit, and worked ITER, ST, spheromak and FRC cases. Paywalled | CerfonFreidberg.pdf |
| Physics of Plasmas 11 (2004) 3510–3518 | https://doi.org/10.1063/1.1756167 | Atanasiu, Günter, Lackner & Miron. Two families of exact solutions, with two free parameters each in the pressure and poloidal-current profiles. **Not a new rung**: the families are chosen so that the Grad–Shafranov equation becomes *linear* inhomogeneous, which is McCarthy's case. Worth having for a different reason — it constructs X-points by **superposition** and is explicitly aimed at benchmarking near them, which is exactly where the extension technique of stage 5 gives out, `Γ` acquiring a corner. Paywalled | Atanasiu-AnalyticalSolutions.pdf |
| Physics of Plasmas 14 (2007) 112508 | https://doi.org/10.1063/1.2803759 | Guazzotto & Freidberg, a family of analytic equilibrium solutions. Listed for completeness and to keep it distinct from McCarthy: an earlier version of this file wrongly credited HDG-GS-1's eq. (14) to this paper. Paywalled | GuazzottoFreidberg.pdf |
| Journal of Computational Physics 243 (2013) 28–45 | https://doi.org/10.1016/j.jcp.2013.02.045 | Pataki, Cerfon, Freidberg, Greengard & O'Neil, a fast high-order Grad–Shafranov solver. For cross-comparison of error levels rather than as a test case. Paywalled | GS-FastHighOrder.pdf |
| Journal of Computational Physics 316 (2016) 63–93 | https://doi.org/10.1016/j.jcp.2016.04.002 | Palha, Koren & Felici, a **mimetic spectral element** solver. The most directly comparable numbers in the literature: it reports convergence for Solov'ev, FRC and spheromak analytic solutions, so its tables can be read against MEQ's. Paywalled | GS-MimeticSpectral.pdf |
| arXiv:2606.11821 (2026) | https://arxiv.org/abs/2606.11821 | Zhang, Xie, Li, Meng, Wang & Wang, *VEQ: a fast parametric Grad–Shafranov solver for fixed-boundary tokamak equilibria with flexible source profiles*. Listed for one idea worth stealing regardless of the solver: it reports **both** a projected residual norm and a *pointwise strong-form* residual, and uses the pair to tell genuine convergence from a projection that has merely stopped moving. Its own abstract is the argument for the pair — enriching the representation improves *interior* force balance while the global RMS and maximum stay "dominated by near-boundary contributions", which is a diagnosis a projected norm alone cannot make. Its unknowns are MXH flux-surface harmonics and shifted-Chebyshev profile coefficients, so the parametrisation is nothing like MEQ's and the diagnostic is the transferable part. **Checked against the arXiv API 2026-09-02**: the id resolves, the title and the six authors are as given, submitted 2026-06-10 — and there is **no published doi**, the entry carrying neither `journal_ref` nor `doi` where arXiv:1605.01538 carries both, with a Crossref title search finding no version of record either. Preprint only, and **not a future-dated id**: 2606 is June 2026, which is in the past. Unread seed | VEQ.pdf |

## The finite element library

| Reference | URL (doi or arxiv) | Short Description | File Name |
| --- | --- | --- | --- |
| Computers & Mathematics with Applications 81 (2021) 42–74 | https://doi.org/10.1016/j.camwa.2020.06.009 | Anderson et al., "MFEM: A modular finite element methods library" — the library MEQ is built on, and the citation to use. Note that MEQ does **not** build against `mfem/master`: it needs the HDG work in `../mfem-hdg-dev`, whose `DarcyForm` replaced the older `HDGBilinearForm` API MEQ used to depend on. `CLAUDE.md` says which branch and why. Open access | MFEM.pdf |
