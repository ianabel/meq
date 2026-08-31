# Normalised flux under linearise-then-condense

A design, not an implementation. Nothing here has been built. Written 2026-08-30.

The task was to design a mechanism carrying `ψ_ax` as an unknown under
`NLOrdering::LineariseThenCondense`, extending to `ψ_bnd` and the `N` Gegenbauer
coefficients of `FREE-BOUNDARY-PLAN.md` §3, **without reintroducing
condense-first by the back door**: no element-local nonlinear solves, no
`GetNumLocalNLIterations()` above zero, no ordering switched quietly mid-solve.

**The design turned out to be much smaller than expected, because the constraint
it was written around is false.** §1 is that measurement, and everything after it
follows from it.

## 1. The premise is false, and the header doc is what misled us

The brief said, and `CLAUDE.md` and `HDG-BEM-COUPLING-FROM-MEQ.md` now both say,
that under linearise-first an auxiliary unknown "enters only through `r_lin`,
which `GetGradient()` caches and `Mult()` does not refresh, so differencing the
reduced residual with respect to it returns **exactly zero**."

**Measured, 2026-08-30. It returns nothing of the kind.** One `GetGradient()` to
set a linearisation point, then repeated `Mult()` at the *same* trace with the
source's normalisation moved:

| case | ordering | `‖∂R/∂ψ_ax‖`, central FD, at `h/ψ_ax` = 1e−4, 1e−5, 1e−6, 1e−7 |
|---|---|---|
| `ν = 2, A = 1, n = 8` | **linearise-first** | 1.104842, 1.104842, 1.104842, 1.104842 |
| | condense-first | 1.106144, 1.106144, 1.106144, 1.106144 |
| `ν = 4, A = 10, n = 16` | **linearise-first** | 1.554135, 1.554134, 1.554134, 1.554134 |
| | condense-first | 1.560985, 1.560985, 1.560985, 1.560985 |
| `ν = 4, A = 100, n = 16` | **linearise-first** | 1.570001, 1.570001, 1.570001, 1.570001 |
| | condense-first | 1.577094, 1.577094, 1.577094, 1.577094 |

Stable to six or seven digits across four decades of step, which is the signature
of a derivative and not of noise, and agreeing with condense-first to 0.4% —
which is the difference between two discretisations of the same quantity, not an
artefact.

**Why the documentation says otherwise.** `darcyhybridization.hpp:647` still
documents the substitution as

```
(q, u)(L)  =  (q, u)_lin + M^-1 ( -r_lin - [C; E] (L - L_lin) )
```

and that is no longer what the code does. `Relinearise()` says so in terms:

> *The local residual at the linearisation point is deliberately not retained. It
> used to be, and `MultInvLin()` used it to predict; that made the retained
> residual enter the substitution twice and cost the gradient its exactness. It
> now enters once, where it belongs, as the correction `MultInvLin()` evaluates
> at the fields it is given.*

So `r_lin` is **not cached at all**. What `MultInvLin()` does is a linear
prediction followed by `corrections` local Newton corrections — one when
evaluating, two when relinearising — and each correction calls `LocalResidual()`,
which calls `LocalNLOperator::Mult()`, which calls
`m_nlfi_p->AssembleElementVector()`, which is meq's `SourceIntegrator`, which
calls `source->f( r, z, ψ )` **freshly**. Any change to the source's internal
state between two `Mult()` calls therefore reaches the residual.

**The header doc is stale rather than wrong-in-spirit** — it describes the
scheme before the fix that made `GetGradient()` exact — and it is worth telling
the MFEM tree, because it is what cost this project a wrong architectural
conclusion.

## 2. What is actually frozen, and what the border is

| frozen between `GetGradient()` calls | re-evaluated at every `Mult()` |
|---|---|
| `lin_trace`, `lin_u`, `lin_p` | the local residual `r`, at the current fields **and the current source state** |
| the factored local Jacobian `M` and its Schur complement | the trace row, at the corrected fields |

With one correction, the fields a trace implies are

```
(u, p)(L, s) = (u, p)_lin + M⁻¹( −[C; E](L − L_lin) ) + M⁻¹( −r( ·, L; s ) )
```

and the reduced residual is the trace row evaluated at them. Differentiating in
`s` with `M`, `C`, `E` and the linearisation point all fixed:

```
∂R/∂s  =  −[ C' E' ] M⁻¹ ∂r/∂s .
```

**That is the whole border, and every factor in it is something this ordering
already forms.** `M⁻¹` is the stored factorisation; `[C' E']` is the trace row;
`∂r/∂s` is element-local and is derived in §3. Compare condense-first, where the
same column is the derivative of a *converged nonlinear local solve* — it exists,
but only through the implicit function theorem and only if that solve converged.

**So linearise-first is not merely compatible with an auxiliary unknown; it is
the better ordering for one.** Under condense-first the border's smoothness is
contingent on every element-local Newton converging, and `CLAUDE.md` records what
happens when they do not: differencing `ψ_ax` by 9e−6 moved `max ψ_h` from 0.896
to 3.84 and the border read 1.6e5 where it should read 1. Under linearise-first
there is no local iteration to fail, so the residual is a smooth function of the
auxiliary unknown **by construction**.

*Not established:* the experiment above was run at the initial-guess trace, where
condense-first is also well behaved and emitted no local non-convergence
warnings. It demonstrates that the border exists and is smooth under
linearise-first; it does **not** demonstrate the robustness advantage deep in a
Newton iteration. The argument for that is structural, not measured.

## 3. `∂r/∂s` needs nothing new from `meq::Source`

A normalised source has the form `F = H( r, z, Ψ )/( ψ_ax − ψ_bnd )` with
`Ψ = ( ψ − ψ_bnd )/( ψ_ax − ψ_bnd )`. Differentiating,

```
∂F/∂ψ_ax   =  − [ F + ( ψ − ψ_bnd ) ∂F/∂ψ ] / ( ψ_ax − ψ_bnd )
∂F/∂ψ_bnd  =  + [ F + ( ψ − ψ_ax  ) ∂F/∂ψ ] / ( ψ_ax − ψ_bnd )
```

**Both are built from `f()` and `dFdPsi()` alone**, which every `meq::Source`
already supplies and which `SourceTests` already checks against a finite
difference. No interface change, no new virtual, nothing for a source author to
get wrong.

Verified symbolically 2026-08-30 against a generic `H`, together with the
identity that makes them checkable at runtime for free:

```
∂F/∂ψ + ∂F/∂ψ_ax + ∂F/∂ψ_bnd  =  0
```

which holds because shifting `ψ`, `ψ_ax` and `ψ_bnd` together leaves `Ψ`
unchanged. **That is the cheapest unit test in this document** and it needs no
solver at all.

The element-local vector is then the potential-block entry

```
( ∂r/∂s )_p  =  − ( (∂F/∂s)( r, z, ψ_h ) / r , w )_K ,
```

on the same quadrature rule `SourceIntegrator` already uses for the residual and
the Jacobian, so the three cannot drift apart through a rule change. The flux
block of `∂r/∂s` is zero: the source touches only the potential equation.

## 4. What MFEM has to expose — and it is much less than was asked for

`HDG-BEM-COUPLING-FROM-MEQ.md` §3 asks for **auxiliary globally-coupled unknowns
carried through the condensation**: `SetNumAuxiliaryUnknowns`,
`AssembleAuxFluxMatrix`, a bordered `GetGradient()`. That is an architectural
change to `DarcyHybridization`.

**§2's algebra says it is not necessary.** All meq needs is the two directions of
the local elimination, exposed as operators:

```cpp
/** @brief The trace-space image of a per-element local perturbation:
        y = -[C' E'] M^-1 g,   g given per element as (g_u, g_p).

    This is the operation MultInvLin()'s correction already performs on the
    local residual, stopped one step earlier and handed back. Exposing it is
    what lets a caller ASSEMBLE an auxiliary column of the bordered Jacobian
    instead of differencing the reduced residual for it.

    Valid only where a linearisation exists, i.e. after GetGradient(). */
void ApplyLocalToTrace(const BlockVector &g_local, Vector &y_trace) const;

/** @brief The fields a trace increment implies, per element:
        (du, dp) = -M^-1 [C; E] dL.

    The prediction half of MultInvLin(), without the correction. It is what an
    auxiliary ROW needs: the sensitivity of a functional of the fields (a
    maximum, a flux through a surface) to the trace. */
void ApplyTraceToLocal(int el, const Vector &dL_faces,
                       Vector &du, Vector &dp) const;
```

Both are element-local, both reuse `MultInv()` and the stored factors, and
**neither requires the hybridization to know that auxiliary unknowns exist.** meq
assembles the border itself and does the bordered solve outside — for which it
already has working code.

**How this refines the existing request.** §3 of `HDG-BEM-COUPLING-FROM-MEQ.md`
should be replaced, not amended: the ask changes from *"carry `M` extra unknowns
through the elimination"* to *"expose the elimination in both directions"*. That
is smaller, more general, has no bearing on the reduced system's size or type,
and serves any global constraint rather than only this one. §2 of that document —
the two rectangular integrators — is unaffected and still stands.

**And nothing is blocked in the meantime.** §1 measured that differencing works
under linearise-first, stably. So the assembled border is an *accuracy and cost*
improvement over a route that already functions, not a prerequisite. That is the
opposite of what the request currently claims.

## 5. Free boundary: the same mechanism, `N + 2` times

| unknown | how it reaches `r` | `∂r/∂·` |
|---|---|---|
| `ψ_ax` | the source's normalisation | §3, exact |
| `ψ_bnd` | the same | §3, exact |
| `a_n`, `n = 1…N` | the transferred datum's data half, `⟨ g∘a, v·n ⟩_e` on `Γ_h` faces | the rectangular block of `HDG-BEM-COUPLING-FROM-MEQ.md` §2.1 — **already an assembled matrix**, and it lands in the *flux* block of `∂r/∂a_n` |

So all `N + 2` columns are `−[C' E'] M⁻¹` applied to an element-local vector, and
one `ApplyLocalToTrace` call each. The transmission rows `T_m` need
`E_h(q_h)·ν` on `Γ` and its sensitivity to the trace, which is `ApplyTraceToLocal`
composed with §2.2's boundary quadrature.

**Where it stops is cost, not correctness.** `N + 2` applications of
`ApplyLocalToTrace` per Newton step is cheap — one local backsolve per element
each, no global solve. What is not free is the bordered *solve*. Two routes, and
the second is better than the fallback `FREE-BOUNDARY-PLAN.md` §4 assumes:

* **Block elimination**: `N + 2` backsolves of `A` against one factorisation.
* **Assemble the bordered sparse matrix** of size `n_trace + N + 2` and solve
  once. Possible precisely *because* the columns are explicit vectors rather than
  the output of a differencing procedure, and because they are sparse — supported
  on `Γ_h`-adjacent trace dofs for `a`, and dense but cheap for `ψ_ax`, `ψ_bnd`.

The second is what to build. It is also what makes the `N ≈ 20–40` of
`FREE-BOUNDARY-PLAN.md` §3.3 affordable without argument.

## 6. The semismooth part, and a hazard specific to this ordering

`ψ_ax = max ψ_h` over the potential dofs, so the argmax element can change
between iterations. That is an ordinary semismooth Newton and the existing
condense-first implementation already handles it: recompute the argmax at each
Jacobian build, use it for the border row, and evaluate the constraint with the
true maximum.

**What is new under linearise-first, and what I have not seen stated anywhere:
the fields at the linearisation point and the fields at an evaluation of the same
trace are one local correction apart.** `MultInvLin()` takes two corrections when
relinearising and one when evaluating, so `lin_u` carries a fields state that a
subsequent `Mult()` at `L = L_lin` will correct once more. Two consequences:

* **Never read `ψ_h` from the linearisation state.** The constraint `G = ψ_ax −
  max ψ_h` must be evaluated from the fields a residual evaluation produced, or
  it is measuring a different function from the one the border differentiates.
* **The argmax may differ between the two.** The gap is `O( M⁻¹ r )` — small
  once the local residual is small, and not small early in the iteration, which
  is exactly when the argmax is most likely to be moving anyway.

*Not established.* Whether that gap ever flips the argmax in practice, and
whether it costs the quadratic tail when it does, is unmeasured. It is the part
of this design I would expect to bite first, and it is cheap to instrument:
record the argmax element at the linearisation and at the evaluation and count
disagreements.

## 7. What this does *not* fix

**It does not make linearise-first solve the pedestal.** Re-measured
2026-08-30, that ordering still fails at 60 on §4.2 at every resolution tried,
including ones where condense-first takes five iterations. This document removes
the *border* as an obstacle to adopting linearise-first; it says nothing about
the obstacle that actually blocks it.

So the honest statement of what has changed: `setSource( NormalisedSource &,
double )`'s refusal of `LineariseThenCondense` is **not justified by the reason
its message gives**, and can be lifted today with no MFEM change at all. Whether
meq should *run* on that ordering remains a separate question, decided by the
pedestal diagnosis and not by this.

## 8. The cheapest experiment that would falsify this design

In increasing cost:

1. **The sum identity of §3**, `∂F/∂ψ + ∂F/∂ψ_ax + ∂F/∂ψ_bnd = 0`, as a unit
   test on a `NormalisedSource`. Needs no solver, no mesh, no MFEM. If it fails,
   §3's algebra is wrong and §5 collapses with it. Minutes.
2. **Lift the ordering guard and run the existing bordered Newton under
   linearise-first**, on `HighBetaConvergence`'s six cases, with the border still
   obtained by differencing. **This is the decisive one.** If §2's substitution
   structure is what this document claims, the observed Newton order is 2 and the
   converged `ψ_ax` matches the condense-first answer to discretisation error. If
   the structure differs, the order drops to 1 — which is the standing result
   that a wrong Jacobian converges to the right answer at the wrong rate, and is
   precisely the signal this test is for. An hour, because the solver already
   exists and only the guard is in the way.
3. **Once MFEM exposes `ApplyLocalToTrace`**, assemble the column and compare it
   against the central-difference column measured in §1. They must agree to the
   difference's own accuracy. That is the check that the assembled border and the
   differenced one are the same object, and it is the one that would justify
   retiring the differencing.

Experiment 2 is the one to run, and it can be run today. Its most likely failure
mode is not the border at all but §6's argmax gap, which is why §6 says to
instrument that count while running it.
