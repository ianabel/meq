# What MaNTA needs from MEQ

A note on the interface a magnetic-field solver has to present in order to be
coupled to [MaNTA](https://github.com/ianabel/MaNTA) as a `FieldModel`. Written
from MaNTA's side, at the point where the coupling machinery is complete and no
field model is registered in that tree: MaNTA supplies the block structure, the
coupled Jacobian and adjoint solves, the DOF bookkeeping and the serialisation,
and something else has to supply the equilibrium. MEQ is the obvious candidate.

Nothing here is a decision MEQ has to accept as given. It is a statement of what
the existing socket is shaped like, and where the shape is negotiable it says so.

The authoritative version of the contract is `FieldModel.hpp`,
`FieldModelSpec.hpp` and `docs/field_coupling.rst` in the MaNTA tree. Where this
note and those disagree, they are right and this is stale.

---

## Which coupling is this? — and the one already in `TODO`

MEQ's `TODO` already carries an entry (raised 2026-08-24) for reading `p'(psi)`
and `FF'(psi)` from a NetCDF file that MaNTA may have written. That is a
**different and looser coupling** from the one described here, and both are
worth having:

| | File-based (`TODO` entry) | `FieldModel` (this note) |
|---|---|---|
| Direction | MaNTA → MEQ only | two-way |
| Convergence | outer Picard loop, or none at all | one Newton on the coupled system |
| Coupling strength | MEQ consumes a frozen profile | equilibrium and transport solved together |
| MEQ runs | standalone, as it does now | in-process, as a plugin MaNTA dlopens |
| Cost | one equilibrium per outer iteration | one residual per transport residual |
| Effort | small; a reader and a config key | substantial; everything below |

The file-based route is the right first step and is useful on its own — it gives
MEQ realistic profiles to solve against without any of the machinery below. It is
also the natural way to *validate* the tight coupling later: the same problem
solved both ways should agree, and the loose one is much easier to debug.

Nothing here supersedes that entry. Treat this as the description of where the
road goes after it.

### One correction to that `TODO` entry while it is open

Its question 3 asks "which group is which?", on the premise that MaNTA's output
groups are generic `Var0`, `Var1`, … with no metadata saying which is pressure.
**That premise is out of date.** MaNTA physics cases name their own variables
now, and the netCDF group takes that name:

```cpp
nc_output.AddGroup(problem->getVariableName(i), problem->getVariableDescription(i));
nc_output.AddVariable(problem->getVariableName(i), "u", "Value",
                      problem->getVariableUnits(i), y.u(i));
```

So a modern case writes groups called `n`, `p_i`, `T_e` — whatever it declared —
each carrying a `description` attribute, with `units` on every variable inside.
`Var0`/`Var1` survive only for the three cases whose *width* comes from their
config and which therefore have no names to give (`MatrixDiffusion`,
`MatrixDiffusionTest`, `LinearDiffSourceTest`), via `numberedFields`. The
reference file that entry was established from happens to be one of those.

That does not remove the need for a config key naming the mapping — a group
called `p_i` still has to be *declared* to be the pressure MEQ wants — but it
does mean the information is there to be matched against rather than absent, and
the units attribute is a free consistency check.

Questions 1 and 2 of that entry stand unchanged, and question 1 (get the MaNTA
case to output `p'` directly, so `u = p'` and `q = p''` line up with
`meq::Knot`) is worth settling early: **it is a decision on the MaNTA side that
costs almost nothing there and saves an order of accuracy here.**

---

## 1. The headline: MaNTA does not want a converged equilibrium

This is the point most likely to be got backwards, so it comes first.

The instinct when coupling an equilibrium solver to a transport code is to give
the transport code a function `equilibrium(pressure_profile) -> geometry`, which
internally runs a Newton solve to convergence. That is the Picard-style outer
loop, and it is *not* what MaNTA is asking for.

MaNTA asks for MEQ's **constraint rows**, not their solution. The Grad–Shafranov
residual joins the DAE that IDA/KINSOL is already solving, MEQ's unknowns join
the solution vector, and **MEQ's Newton is subsumed into MaNTA's**. The
equilibrium is converged only when the whole coupled system is converged, once,
rather than being re-converged inside every transport residual evaluation.

Concretely: `FieldResidual` is called on *every* residual evaluation, at states
that are nowhere near an equilibrium, and it must return the residual at that
state without complaint. Its Jacobian is applied once per Newton iteration
rather than iterated to convergence per call.

That is what makes the coupling affordable, and it is also what makes the demands
on MEQ different from the demands on a standalone solver. A standalone solver
needs to converge. A coupled one needs to be **differentiable and evaluable off
the solution manifold**.

---

## 2. The coupled system

A field model supplies `nFieldDOF` unknowns $\psi$ with one residual row each,

$$R_m(\psi, \dot\psi, u, q, \sigma, \phi, t) = 0, \qquad m = 1 \ldots \texttt{nFieldDOF}$$

together with a map from $\psi$ to `nGeometry` **geometry slots** $g_s(\psi, x, t)$
that the transport physics reads. The coupling is two-way and both directions go
through a declared interface: the field rows may read the transport solution, the
transport hooks may read the geometry, and neither reaches into the other's
unknowns directly.

The Jacobian is

$$\begin{pmatrix} A & A_1 \\ A_2 & B \end{pmatrix}
  \begin{pmatrix} \delta y \\ \delta\psi \end{pmatrix}
  = \begin{pmatrix} r_y \\ r_\psi \end{pmatrix}$$

| Block | Is | Supplied by |
|---|---|---|
| $A$ | the HDG transport operator | MaNTA (static condensation, as always) |
| $A_1 = \partial(\text{transport})/\partial\psi$ | chain rule through geometry | physics case's `d*_dGeometry` × MEQ's `dGeometry_dpsi` |
| $A_2 = \partial R/\partial y$ | how the equilibrium sees transport | **MEQ** |
| $B = \partial R/\partial\psi + \alpha\,\partial R/\partial\dot\psi$ | MEQ's own block | **MEQ** |

$A_1$ factorises across the interface — MaNTA's physics case reports
$\partial(\text{row})/\partial g_s$ and MEQ reports $\partial g_s/\partial\psi_m$ —
so neither side needs to know the other's internals.

---

## 3. The interface, hook by hook

Nine pure virtuals and six with defaults. Deriving from `FieldModel`:

### Required

| Hook | Signature sketch | What MEQ must do |
|---|---|---|
| `FieldResidual` | `(out, psi, dpsidt, states, points, weights, t)` | The GS residual rows at this $\psi$ and this transport state. `out` is length `nFieldDOF`, arrives zeroed. |
| `Geometry` | `(out, psi, x, t)` | The metric at **one** point $x$. Length `nGeometry`, arrives zeroed. See §5 — this is called a lot. |
| `dGeometry_dpsi` | `(out, psi, x, t)` | Shape `(nGeometry, nFieldDOF)`, arrives zeroed. The hard one; see §6. |
| `FieldResidualPrime` | `(dR, dRdot, dRdpsi, dRddpsidt, psi, dpsidt, states, points, weights, t)` | Every row's derivatives at once. All four arrive zeroed. |
| `InitialFieldValue` | `(out)` | The starting guess for $\psi$. Not called on a restart. |

### Overridable, and the real seam

| Hook | Default | Why MEQ overrides it |
|---|---|---|
| `updateFieldJacobian(dRdpsi, dRddpsidt, alpha)` | forms $B$ densely, partial-pivot LU | MEQ's block is large and structured; do not form it |
| `applyB` / `applyBTranspose` | dense matvec | ditto |
| `solveB` / `solveBTranspose` | dense LU solve | **this is MEQ's existing hybridized solve** |
| `resetForRun()` | no-op | mandatory if MEQ caches anything; see §8 |

**MaNTA never needs $B$ itself — only the ability to apply and invert it, in
both directions.** That is the whole design of the seam, and it is why a
large-DOF solver can plug in at all.

### The transport state MEQ receives

`states` is a `GlobalState`: the transport solution ($u$, $q$, $\sigma$, aux)
sampled on the physics nodes, `points` are the abscissae, and `weights` is one
quadrature weight per node, so

$$\int f \, \mathrm{d}x = \texttt{weights.dot(f\_at\_nodes)}.$$

**Use those weights rather than a rule of your own.** This is a rule MaNTA
learnt the hard way on its scalar hooks: a physics case there used a global
adaptive Kronrod rule over a piecewise polynomial, which is not a smooth function
of the coefficients, and disagreed with its own Jacobian by 8%. The weights are
what makes $\partial(\int u)/\partial u_j = w_j$ exactly, which is what
`FieldResidualPrime` has to report.

---

## 4. What is $\psi$? — choosing the DOF vector

This is the first real design decision and it is MEQ's to make.

`nFieldDOF` unknowns are appended to MaNTA's solution vector, after the global
scalars:

```
[ sigma | q | u | aux ]   per cell, for each cell
[ lambda ]                all face traces
[ mu ]                    all global scalars
[ psi ]                   MEQ's unknowns
```

They are carried in the IDA/KINSOL vector, written to the netCDF output and the
restart file, and each declares itself differential or algebraic (which reaches
`IDASetId`).

**The natural choice looks like MEQ's condensed trace unknowns.** MEQ already
hybridizes and condenses element interiors away, leaving a global system in the
trace alone; that global system is exactly what $B$ has to apply and invert, and
`solveB` would then be MEQ's existing hybridized solve rather than anything new.
The interior recovery MEQ already does element-by-element is then what `Geometry`
needs in order to evaluate at an arbitrary $x$.

Three consequences to weigh before committing:

* **Every field DOF is in IDA's local error test** unless `SuppressAlgebraicError`
  is set, and MaNTA's measurements say turning that on is a trade rather than an
  improvement (a restart round trip degrades from 1.9e-6 to 8.6e-4). With a large
  trace space this could plausibly dominate step-size control. Untested — there is
  no fixture in MaNTA anywhere near this regime.
* **The restart file carries all of them.** Fine, but sizeable.
* **The exact coupled solve costs `nField + 1` transport solves.** See §9.

**Grad–Shafranov is algebraic**: there is no $\dot\psi$ in the equation, so every
DOF declares `differential = false`. Getting this wrong is a specific, documented
failure — declaring a DOF differential when its row carries no time derivative
gives `IDA_LINESEARCH_FAIL` (-13), because `IDA_YA_YDP_INIT` holds every
differential *value* fixed and a row that can reach no unknown it may move is
irreducible. MaNTA refuses it by name at `initialize()` rather than letting IDA
report it, but only because it checks the residual, not the declaration.

If the evolving-equilibrium case ever wants a $\dot\psi$ term, the slot
(`dRddpsidt`) is already there and already weighted by $\alpha$.

---

## 5. Geometry: 2-D equilibrium to 1-D metric

MaNTA is 1-D in a flux label; MEQ is 2-D in $(r, z)$. The geometry slots are
where the dimensional reduction happens, and MEQ owns it.

`FieldModelSpec::label` names the spatial coordinate the geometry is expressed
against — MaNTA does not interpret it, but records it in the output so a run says
what its $x$ meant. So MEQ declares its own flux label and supplies the metric on
it.

**The concrete slot list is negotiated with the transport physics case, not
fixed by MaNTA.** For orientation, the interface MaNTA's mirror physics used
before this feature existed was, on a volume coordinate $V$:
$\psi(V)$, $B(V,z)$, $R(\psi,s)$, $\mathrm{d}R/\mathrm{d}V$, $V'(V)$, the mirror
ratio, $L_\parallel(V)$, $R_\mathrm{max}$, $R_\mathrm{min}$. For a tokamak the
expected set is the usual 1-D transport metric on $\rho$:

$$V'(\rho),\quad \langle|\nabla\rho|\rangle,\quad \langle|\nabla\rho|^2\rangle,\quad
\langle R^{-2}\rangle,\quad \langle B^2 \rangle,\quad q(\rho),\quad F = R B_\phi,\quad S(\rho)$$

The *shape* of the requirement is what matters: a set of scalar functions of the
flux label, each differentiable with respect to every $\psi$ DOF.

### The call pattern is the performance problem

`Geometry` is **pointwise**: MaNTA calls it once per physics node, per residual
evaluation, handing the whole $\psi$ vector each time
(`SystemSolver::evaluateGeometry` is a plain loop over `points`).

For a flux-surface average that means locating the surface through $x$ and
integrating on it — per node, per residual. That is not affordable done naively,
so **MEQ will need to cache per $\psi$**: compute the surface geometry once when
$\psi$ changes and serve individual points from it. The cache invalidation
contract is §8.

If the pointwise signature turns out to be the wrong shape — if MEQ would much
rather be asked for all nodes at once — that is worth raising on the MaNTA side.
Every other physics hook in MaNTA exists in both pointwise and batched forms; the
geometry hook is the exception, and not obviously deliberately so.

---

## 6. `dGeometry_dpsi`, and the moving-surface term

This is the requirement most likely to be underestimated.

The geometry slots are flux-surface averages. A flux-surface average is an
integral over a surface **whose location is itself determined by $\psi$**. So

$$\frac{\partial}{\partial \psi_m} \langle f \rangle_\rho
 = \underbrace{\left\langle \frac{\partial f}{\partial \psi_m} \right\rangle}_{\text{integrand moves}}
 + \underbrace{(\text{terms from the surface moving})}_{\text{easy to forget}}$$

Both terms are needed. Dropping the second gives a Jacobian that is wrong but
*plausible* — and MaNTA's architecture has a specific, unpleasant property here
worth knowing about:

> **An error in the forward Jacobian does not produce a wrong answer, only slow
> Newton convergence.** MaNTA never assembles its Jacobian, so a missing block
> shows up as a solver that takes more iterations, not as a bad number. Several
> defects in MaNTA survived a passing regression suite for months for exactly
> this reason.

The corollary is that **`dGeometry_dpsi` must be finite-difference tested
directly**, against `Geometry`, as a unit test. Convergence of the coupled solve
is not evidence that it is right. MaNTA's own equivalent tests are
`SolveJacTests.cpp` (finite-difference the residual, require $J\,\delta y = g$)
and MEQ already has the same discipline — its existing unit layer differences
$\partial F/\partial\psi$ against a finite difference, which is precisely the
right instinct applied to the right quantity.

The **adjoint** raises the stakes: there, a missing block *does* produce a
silently wrong gradient beside a perfectly good objective. MaNTA has been bitten
by exactly this (a `dSigma/dPhi` block absent from the adjoint assembly, costing
nothing visible until a test was written for it).

---

## 7. $A_2 = \partial R/\partial y$: how the equilibrium sees transport

The physical content of the coupling in this direction is that the GS source term

$$F(r, z, \psi) = \mu_0 r^2 \frac{\mathrm{d}p}{\mathrm{d}\psi} + g\frac{\mathrm{d}g}{\mathrm{d}\psi}$$

depends on the pressure, and the pressure is what MaNTA is computing. So $p(\psi)$
stops being a prescribed profile function and becomes a functional of the
transport solution — $p = \sum_s n_s T_s$ over whatever species the physics case
carries, evaluated on the flux label.

`FieldResidualPrime` must therefore report $\partial R_m / \partial(\text{transport DOF})$
for every field row and every node. In MaNTA's data layout that is written
through `dR[m].Variable()`, the whole `(nVars, nPoints)` matrix — **not** through
`dR[m][j].u(0)`, which compiles, modifies a temporary, and silently discards the
write. (`GlobalState::operator[]` returns a `State` *by value*. This trap is
documented in MaNTA's test fixtures because it was hit.)

Two notes:

* **`dRdot` cannot be filled and leaving it zero is correct.** `FieldResidual`
  receives `states` and no `states_dot`, so a field row has no way to depend on
  transport time derivatives in the first place. The slot exists because the
  assembly already weights it by $\alpha$. What must *not* happen is putting
  $\partial R/\partial\dot\psi$ there — that belongs in `dRddpsidt`, and written
  into `dRdot` it lands in the coupling row at entirely the wrong DOFs with
  nothing to say so.
* **How $p$ is built from $(n, T)$ is a physics-case question**, not MEQ's. The
  cleanest split is probably that MaNTA's physics case owns the profile mapping
  and MEQ receives $p(\rho)$ and $\mathrm{d}p/\mathrm{d}\psi$ — but that is a
  conversation, and it determines which side owns which half of the chain rule
  in $A_2$.

---

## 8. Two contracts that are easy to miss

### Fail by throwing, never by fudging

A model that **cannot evaluate at this state** — no x-point, separatrix left the
domain, Newton inside a nested solve diverged — **must throw** from
`FieldResidual`. MaNTA catches it and returns 1 to IDA, which treats that as a
*recoverable* error and retries with a smaller step.

Returning a plausible-looking nonsense value instead converts a recoverable step
into a wrong answer. This matters more than usual here because §1 guarantees MEQ
*will* be called at states far from equilibrium, routinely, as a normal part of
the Newton iteration.

### `resetForRun` is not optional if MEQ caches anything

MaNTA's `initialize()` skips matrix setup when the solver has already been
initialised, so anything computed once per *object* rather than once per *run* is
stale on the second run. Given §5 says MEQ will certainly cache, this applies.

`resetForRun()` is called from the unconditional part of `initialize()` for
exactly this reason. MaNTA pins it with a test asserting that a reused coupled
solver matches a fresh one **bit for bit** — zero tolerance, deliberately,
because the last defect of this kind left the second run completing, plausible,
and wrong in the eleventh digit.

---

## 9. Cost, and why the iterative path exists

`FieldSolve` picks between two routes to the same answer:

* **`exact`** — form the Schur complement onto $\psi$ by applying the transport
  inverse to every column of $A_1$. That is `nField + 1` transport solves per
  Jacobian solve. It is a verification tool and the oracle the other path is
  checked against.
* **`iterative`** (default) — block Gauss–Seidel between the blocks with
  Irons–Tuck acceleration, one transport solve per sweep.

Break-even is `#sweeps < nField + 1`. **No fixture in MaNTA is on the winning
side of it** — iterative measures 1.5× more expensive at `nField = 1` and
2.2–6.3× at `nField = 5`, for the same answer.

That is not a criticism of the iterative path; it is a statement that MaNTA has
never yet run the regime the iterative path was built for. **MEQ is that regime.**
With $N_\text{magnetics} \gg N_\text{HDG}$ the exact route's `nField + 1`
transport solves and $O(\texttt{nField}^3)$ dense factorisation are hopeless, and
the iterative path becomes the only viable one. MEQ coupling is the first real
test of whether the bet was right.

`iterative` is a cost choice and never an accuracy one: a sweep that exhausts its
cap escalates to the exact solve rather than returning an under-converged answer,
in both the forward and adjoint directions.

---

## 10. Adjoints

If gradients are wanted through the coupled system, MEQ needs **both directions**
of its block operations — `applyBTranspose` and `solveBTranspose` beside the
forward pair. A model supplying only one direction cannot be silently
accommodated, which is why they are declared separately rather than inferred.

The coupled adjoint is the transpose block for block; the elimination runs the
other way and the Schur complement onto $\psi$ becomes

$$\left(B^T - A_1^T A^{-T} A_2^T\right) z_\psi = G_\psi - A_1^T A^{-T} G_y$$

For a hybridized solver the transpose solve is usually cheap to provide — the
condensed system is what it is — but it is worth confirming early rather than
discovering late that only the forward direction was ever built.

---

## 11. What MaNTA cannot do yet

Honest list, so none of these arrive as surprises:

* **`nScalars > 0` alongside a field model is refused** at `setFieldModel`. The
  reason is a branch disagreement rather than a missing feature (one
  `dSources_dScalars` path builds a `State` with no geometry rows), so it is
  fixable, but it is refused today. If MEQ coupling needs global scalars — a
  total current constraint, a plasma current control loop — this needs doing
  first.
* **A field model cannot be written in Python.** `FieldModel` has no pybind11
  class, and `FieldModel` is a `ProblemSelection` key, so a coupled run is a
  config-file run. Not a problem for MEQ, which is C++.
* **An objective whose integrand reads geometry directly loses its
  $\mathrm{d}G/\mathrm{d}\psi$ term** — `AdjointProblem` reports derivatives with
  respect to $u$, $q$, $\sigma$ and $\phi$, and geometry is not among them. Zero
  by construction today rather than by assumption.
* **A field model cannot depend on an adjoint parameter**, so
  $\partial R/\partial p$ is zero. Relevant the moment anyone wants to optimise a
  coil current.
* **`DegreeAdaptation` with a field model is refused** — the adaptive driver
  builds a solver per degree and has no way to carry the model.
* **`Superconvergent = true` with spatial adjoint parameters** throws, as it does
  without a field model.

---

## 12. Integration mechanics

The good news: **MaNTA's build does not need to know about MFEM.**

`FieldModel.hpp` and `FieldModelSpec.hpp` are both in MaNTA's `INSTALL_HEADERS`,
and `runManta` dlopens everything in the config's `PhysicsPlugins` array *before*
instantiating the field model — deliberately, so that a model can come from a
plugin. So the integration path is:

1. `make install PREFIX=...` in MaNTA.
2. MEQ builds an additional shared object containing a `FieldModel` subclass,
   compiled against the installed MaNTA headers, linked against MEQ and MFEM.
3. A MaNTA config names it:

```toml
[configuration]
PhysicsPlugins = ["/path/to/libmeq-fieldmodel.so"]
FieldModel = "MEQ"
```

Registration is a static-initialisation side effect
(`REGISTER_FIELD_MODEL_HEADER` / `_IMPL`), the same pattern as a physics case.

**Two traps, neither of which is a link error:**

* **The plugin must be compiled with the flags `pkg-config --cflags manta`
  reports.** Eigen types cross the interface (`VectorRef`, `MatrixRef`, `Vector`
  are all Eigen), Eigen aligns to the widest vector unit the compiler knows about,
  and it inlines expression templates into both sides of the boundary. A plugin
  built without the core's `-march=` faults inside an aligned AVX-512 load the
  first time the solver touches its state. `manta.pc` records the *concrete*
  architecture so a mismatch is a compile error rather than a run-time crash.
  `-DEIGEN_USE_BLAS` travels for the same reason. **This needs reconciling with
  MEQ's CMake flags and with how MFEM was built** — worth checking early, it is
  the kind of thing that costs an afternoon.
* **The plugin must not link `-lmanta`.** MaNTA links its core objects directly,
  so a plugin pulling in `libmanta.so` gets a second copy of the registry and
  registers into a map the solver never reads — silently. Compile against the
  headers alone and let the loader bind MaNTA's symbols to the host process,
  which is linked `-rdynamic` for exactly that.

---

## 13. Checklist

What MEQ must be able to do, reduced to a list:

- [ ] Evaluate its residual at an **arbitrary, non-converged** $\psi$ and an
      arbitrary transport state, without solving anything to convergence
- [ ] Expose a DOF vector $\psi$ (probably the condensed trace) that MaNTA can
      own, carry in the DAE vector and serialise
- [ ] Report $\partial R/\partial\psi$ **as an operator**: apply, solve, and both
      transposes — never assembled
- [ ] Report $\partial R/\partial(\text{transport DOFs})$ on the physics nodes,
      using MaNTA's quadrature weights
- [ ] Evaluate flux-surface-averaged geometry at an arbitrary point, cheaply
      enough to serve one node per call (so: cached, keyed on $\psi$)
- [ ] Differentiate that geometry with respect to every $\psi$ DOF, **including
      the moving-surface term**, and finite-difference test it
- [ ] Throw rather than return nonsense when the state is not one it can evaluate
- [ ] Discard all caches on `resetForRun()`
- [ ] Declare its DOFs algebraic (no $\dot\psi$ in Grad–Shafranov)

The two genuinely hard ones are **`dGeometry_dpsi` with the surface-motion term**
and **making `Geometry` fast enough for the pointwise call pattern**. Everything
else is bookkeeping.

---

## 14. Open questions for the MaNTA side

Things this note could not settle, listed so they are not lost:

1. **Should `Geometry` gain a batched form?** Every other physics hook in MaNTA
   has one. The pointwise-only geometry hook looks like an oversight rather than
   a decision, and §5 says it is the one that will hurt.
2. **Does the error test need to exclude the field block?** `SuppressAlgebraicError`
   is all-or-nothing across $\sigma$, $q$, $\lambda$, $\phi$ *and* the field DOFs.
   A large trace space may need to be dropped from the error test without
   dropping the rest — nothing in MaNTA has ever run at a scale where this is
   testable.
3. **Who owns the $p(\psi)$ mapping** — the MaNTA physics case or the MEQ field
   model? This determines which side owns which half of the $A_2$ chain rule
   (§7).
4. **Does `nScalars > 0` need unblocking first?** A current constraint or a
   control loop would need it (§11).
