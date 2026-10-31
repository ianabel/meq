# Anisotropic pressure, in MEQ, against geq

**DRAFT.** Nothing here is built. No stage is scheduled, and `ROADMAP.md` has
not yet placed this plan in its order.

**The target**: MEQ solves the anisotropic generalisation of RoPP (136) on the
free-boundary rotating mirror. The answer is checked against geq's
`PrescribedAnisotropy`, which freegs4e's `ProfilesCentrifugalMirror` folds into
its `Jtor`. `MIRROR-PLAN.md` is the machine, the gauge and the comparison
protocol, and this plan reuses all three. `CLAUDE_FLOW.md` is the rotating
source. `MEASUREMENTS.md` will hold the tables.

**The reference is pinned**: `refs/AlphaGK.pdf`, eq (6.18). It is an
unpublished draft, and `refs/Refs.md` carries what it says, its two typos
(confirmed by the first author), and the check that exposes them. geq and freegs4e cite a later version's numbering,
Appendix B, (B18).

## The equation, in MEQ's convention

```
−∇̄·( (1 − Δ)/R ∇̄ψ ) = F/R

F = μ₀R² Σ_s ∫d³w m_s w_∥² ∂F_s/∂ψ|_{ε_s, μ_s}
  + μ₀R⁴ ω ω′ Σ_s m_s n_s
  + Ī Ī′/( 1 − Δ )

Δ = μ₀ Σ_s ( p_∥s − p_⊥s )/B²
```

The left-hand side is (6.18)'s, because `R ∇̄·( (1−Δ)/R ∇̄ψ ) =
(1 − Δ)Δ*ψ − ∇̄Δ·∇̄ψ`. `F` is (6.18)'s right-hand side, negated and in SI, with
the rotation term's `R⁴` that the draft drops. **At `Δ = 0`, with every species
Maxwellian, `F` is `docs/rotation.rst`'s `F` term for term**: (6.19) is the
bracket, and `Ī` is `g`.

**Every species may carry a `Δ_s`.** The draft sums over its energetic species
because of its ordering. An equilibrium has no such restriction: a species is
whatever `F_s` it has, a Maxwellian is the case `Δ_s = 0` with (6.19) as its
integral, and nothing in the equation tells a thermal species from a fast one.
**So MEQ has no "fast species" concept.** A population with its own
distribution is its own species. geq's `fast_fraction` is bookkeeping for
carrying fast ions inside the bulk density, and it is translated once, at the
reference (AN-4), rather than built into MEQ.

**The source is at fixed `( ε, μ )`, not at fixed `B`** (page 14), because `ε`
carries the centrifugal energy. MEQ names it for what it is. freegs4e's
`dp_par/dpsi|_{B,Psi}` is that misnaming.

**`B` enters through `Δ`.** With `g = 0`, `B = |∇̄ψ|/R`. With a toroidal field,
`B² = |∇̄ψ|²/R² + Ī²/( R²(1 − Δ)² )`, so a `Δ` that depends on `B` is implicit
in itself. **The mirror has `g = 0` and avoids this entirely**, which is one
reason the mirror comes first and the tokamak is not on this plan's path.

## What changes in MEQ

### 1. The flux mass, and with it the meaning of `q`

`GradShafranovSolver::buildForms()` puts `( R q, v )` in `DarcyForm`'s flux
mass, the inverse of the diffusion coefficient `1/R`, with `q = (1/R)∇̄ψ`. The
anisotropic operator makes that `( R/(1 − Δ) q, v )`, with
`q = ( (1 − Δ)/R )∇̄ψ`.

* **`Δ` depends on `ψ`**, through `Ψ` and therefore through both labels, and
  later on `B = |q|/(1 − Δ)`. **So the flux mass becomes non-linear and couples
  to the potential.** MEQ uses only the potential-mass non-linear slot today
  (`GetPotentialMassNonlinearForm()`). The installed MFEM has the block slot,
  `DarcyHybridization::SetBlockNonlinearIntegrator()`, and a block form forces
  `LocalOpType::FullNL`. **`FullNL` has no factored Schur complement**
  (`darcyhybridization.hpp`, beside `LocalOpType::FluxNL`'s), so every element
  solve gets dearer. That cost is measured at AN-2, not assumed.
* **The linear `( R q, v )` must leave the bilinear form whole** once the block
  law is on. *A nonlinear potential mass and a linear one do not mix* is the
  trap on the potential side (`CLAUDE.md`, *Traps*). Whether the flux side has
  the same trap is unmeasured, so it is assumed until shown otherwise.
* **Every consumer of `q` assumes `q = (1/R)∇̄ψ`.**
  * Unaffected, because `1 − Δ > 0`: `CriticalPointFinder`'s zeros of `q`.
  * Unaffected, because `Δ = 0` there: the transfer on `Γ_h` and the DtN.
  * Affected, or not yet known:
    * **B** in the `.nc`, and the band's continuation of it;
    * the estimator's flux terms;
    * `Reconstruct()`, which reads the flux mass;
    * the warm start's `q` seed;
    * the axis search;
    * `FluxSurfaces`.
  * **AN-1 starts with that audit**, and the list above is a starting point,
    not its result.

**The route not taken: geq's.** geq keeps `Δ*` and moves the anisotropy into
the current, `∇Δ` lagged one Picard step. Its divergence form is the same idea
with no division. MEQ's whole method is a Newton on an exact Jacobian, and
freegs4e measures the lag failing at `Δ = +0.2` (a limit cycle). So MEQ puts `Δ`
where the operator has it. **This is declined on argument, not on
measurement.** If AN-2's `FullNL` cost proves prohibitive, the lagged source is
the fallback to measure against it.

### 2. Coil subtraction stops being free

The shipped mirror subtracts geq's 20 filaments. Under the vacuum operator
`ψ_c` is annihilated. Under the anisotropic one it is not, wherever `Δ ≠ 0`.

Write `q = q_c + q_r`, with `q_c = (1/R)∇̄ψ_c` the vacuum flux, divergence-free
away from the filaments. The divergence equation of the remainder is unchanged.
Its flux equation gains a known load:

```
( R/(1 − Δ) ) q_r − ∇̄ψ_r = −( Δ/(1 − Δ) ) ∇̄ψ_c
```

The load is supported where `Δ` is, which is inside the plasma and away from
every filament. It is non-linear in `ψ` through `Δ`. **`Δ` and `B` read the
total flux**: `B = |q_r + q_c|/(1 − Δ)`.

This is derived here and not yet checked. **AN-0 checks it by substitution**,
and `DarcyForm`'s `−q` sign convention (`CLAUDE_HDGGS.md`) decides the sign in
the assembly.

### 3. The source

`NormalisedRotatingSource` gains, per species:

* **`Δ_s`**, prescribed;
* **the species' source integral**, prescribed, in place of its Maxwellian
  bracket. A Maxwellian species keeps the bracket, and `Δ_s = 0`.

The rotation term keeps every species' density. That matches the draft's
(6.17), which carries `m_α n_α R² ω ω′` for the energetic species too. It also
matches freegs4e, which leaves the rotation term unscaled on purpose.

Prescribed means a function of `( Ψ, R, z )`, and of `B` from AN-6. A kinetic
closure, geq's `KineticClosure` tables from mirrorfp, is a later consumer of
the same interface and not on this plan.

### 4. The stabilisation

A box support defaults to `τ = 3/max(R, 0.1 m)`, measured for the coefficient
`1/R` (M-293). The coefficient becomes `(1 − Δ)/R`, and `τ(1 − Δ)` is the
obvious candidate. **The `τ` precedents (M-286, M-292, M-293) say to measure it
and not to assume it**, and AN-3 does.

### 5. Refusals

* **`1 − Δ ≤ 0` anywhere in the support is the firehose limit**, where the
  operator stops being elliptic. MEQ refuses it before a solve and at every
  residual. It does not clip, which is where geq differs:
  `PrescribedAnisotropy` clips at `delta_max = 0.8`, with a warning.
* **A `Δ` that does not vanish at a support edge is a coefficient jump.** HDG
  represents a jump exactly when the jump lies on element faces, which
  `[mesh.generate] AlignPlasma` arranges for the box's `z` edges. The labels
  are flux surfaces and are not faces, so a jump there costs order the way a
  stepped source does (M-09). **Whether to refuse it, or to warn and report it,
  is decided at AN-3** by what it measurably costs.

## The reference: what geq solves

**The pins already carry it.** freegs4e `9a1ee0c` is the anisotropic-`Jtor`
commit itself, and geq `73b9a72` has `PrescribedAnisotropy`. So the benchmark
stays on the pins every reference so far was regenerated at
(`tools/geq-benchmark/README.md`). geq's HEAD is 89 commits on, and nothing
here needs it.

**Six facts about geq set the reference.**

1. **`Δ` is lagged, in one of two rearrangements.**
   * The historical one divides by `1 − Δ` and differences `∇Δ` with
     `np.gradient`.
   * `aniso_divergence_form=True` adds a lagged `−∇·((Δ/R)∇ψ)/μ₀`.
   * They agree only at convergence. **Their difference is a control** on
     geq's own discretisation of `∇Δ`.
2. **`Δ > 0` fails in the historical form**: `+0.2` is a Picard limit cycle,
   measured in freegs4e's docstring on the GF1 mach-6 case. The divergence
   form at `Δ > 0` is unmeasured. **The comparison is posed at `Δ < 0`**, the
   pickup-ion sign, where freegs4e measured `−0.1` to `−0.6` converging.
3. **`∇Δ` is taken on the unmasked field.** A `Δ` that does not decay at the
   mask edge draws a warning, and the jump is not represented.
4. **The box end differs between the codes.** geq's support ends at its domain
   edge, `|z| = 1.5`. MEQ's box ends at `z_edge = 1.5`, inside `Γ = 3.0`. A
   `Δ` that is non-zero there is a coefficient jump in MEQ, and invisible to
   geq's Green's-function boundary. **So `Δ = δ₀ g( ψ_n ) h( z )`**, with `g`
   vanishing at both labels to an order and `h` vanishing before `|z| = 1.5`.
   That is the label precondition MR-4 established, extended to the box end.
5. **geq's input file cannot pose it.** At HEAD, `[physics.anisotropy]` offers
   only a constant or a Gaussian `Δ( ψ_n )`; a Gaussian does not vanish at the
   labels, and the section is not at the pin anyway. **The cases are per-case
   scripts with a callable `Δ`**, in the shape of
   `cases/geq_<shape>.py`, calling only geq's and freegs4e's interfaces.
6. **`fast_fraction` scales the named species' pressure bracket by `1 − f`
   and leaves its rotation term and quasineutrality alone.**
   `dp_par_dpsin` is then the fast population's full source term. The two are
   given together or not at all. **The first comparison sets both to zero**,
   which isolates the operator and needs no translation. The translation into
   MEQ's species is AN-6's.

**The flux variable flips.** geq's `ψ_n` is MEQ's `1 − Ψ`. `Δ`'s shape is
mirrored, and `dp_par_dpsin`'s sign is negated, exactly as MR-4's profiles were.

## Stages

### AN-0 — the equation, pinned and checked — DONE

**DONE.** `tools/anisotropic_equation.py` runs all four checks and exits
non-zero on any failure. It **derives** (6.19) from a Maxwellian at fixed
`( ε, μ )` rather than transcribing it, and each typo check runs in both
directions: the corrected form holds and the printed one does not.
`tests/analytic/ManufacturedAnisotropic.hpp` carries two `Δ` models on Example
5's flux and box. The position model is AN-1's. The flux model,
`Δ = d ψ²`, is AN-2's, built now because the same script checks it at no cost.
`AnisotropicConvergence` checks the fixture four ways:

* against a central difference of its own operator;
* against the script's symbolic `F`, to 1e-13;
* against Example 5 at `Δ = 0`;
* that `1 − Δ > 0.5` on the box.

Flipping the sign of the `∇Δ·∇ψ` term fails both the operator check and the
symbolic check.


* **The reduction**: `Δ = 0`, every species Maxwellian, reduces the equation
  above to `docs/rotation.rst`'s `F` term for term, by sympy. This is the check
  that catches both of the draft's typos, and it is asserted with them
  corrected.
* **The identity** `R ∇̄·( (1−Δ)/R ∇̄ψ ) = (1 − Δ)Δ*ψ − ∇̄Δ·∇̄ψ`, by sympy.
* **The subtraction load** of *What changes* §2, by substitution.
* **A manufactured fixture**: an exact `ψ`, a prescribed `Δ( R, z )`, and the
  `F` the anisotropic operator makes of them. It is checked by 50-digit finite
  differences, as `MaschkePerrin.hpp` was, and lives beside `Soloviev.hpp`
  under `tests/analytic/`.

**Acceptance**: all four exact. No solver code changes in this stage.

### AN-1 — the operator, with `Δ` a function of position — DONE

**DONE** → **[M-328](MEASUREMENTS.md#m-328)**.

**The interface.** `meq::Anisotropy`, MFEM-free, gives `Δ( R, z )`.
`GradShafranovSolver::setAnisotropy()` weights the flux mass by
`R/(1 − Δ)`, as `fluxMassRadius` times a `1/(1 − Δ)` coefficient. **So `q` is
`( (1 − Δ)/R )∇̄ψ`**, and `flux()` says so. The flux mass stays a bilinear form,
so the linear path and the Newton path both take it.

**Refusals.** `buildForms()` refuses `1 − Δ ≤ 0` at every vertex and every
point of an order-`2k + 2` rule, serially, before any assembly; an exception
cannot leave the threaded loops. It also refuses an anisotropy beside any of:

* `setExtension()`;
* `setExteriorCoupling()`;
* `setConductorField()`;
* a Device.

**The audit of `q`'s consumers**, as built:

| consumer | under an anisotropy |
|---|---|
| the warm start's `q` seed | weighted like the flux mass, so the seed is the `q` the solve carries |
| `Reconstruct()`, `ψ*` | consistent with no change, since it reads the forms: `ψ*` converges at `k+2` |
| `CriticalPointFinder` | unchanged: the zeros of `q` are `∇ψ`'s, and at a zero the Hessian is scaled by `1/(1 − Δ) > 0`, so the classification stands |
| `ResidualEstimator` | **refused**: its terms are the isotropic operator's |
| `ContourTracer` and the flux surfaces | **refused**: they read `|q|` as `B_p` |
| the extension on `Γ_h`, the DtN, `ConductorField` | refused with their setters, until AN-3 |
| `Field.hpp`'s B and the driver's output | **unreachable**: the driver cannot pose an anisotropy. AN-6 divides `1 − Δ` out |

**The constant-`Δ` statement, corrected.** The draft of this stage claimed
`ψ_Δ = ψ₀/(1 − Δ)`. **That is wrong twice.** It ignores the boundary data, and
in HDG the stabilisation does not scale with the operator. Substituting
`q = (1 − Δ) q′` divides the conservation and constraint rows by `1 − Δ`. So a
constant `Δ` **is** the isotropic problem with `F/(1 − Δ)` and `τ/(1 − Δ)`,
with `q_Δ = (1 − Δ) q₀`, exactly, and that is what is asserted. **It is also
the first evidence for *What changes* §4**: to keep the stabilisation in
proportion to the diffusion, `τ` scales with `1 − Δ`.

### AN-2 — `Δ` a function of the flux

**Split in two**, so the non-linear flux law is measured before the labels
are added on top of it.

#### AN-2a — `Δ( R, z, ψ )`, physical `ψ`, plain Newton — DONE

**DONE** → **[M-329](MEASUREMENTS.md#m-329)**.

**The interface.** `meq::Anisotropy` is now general: `delta( R, z, ψ )`,
`dDeltaDPsi()` and `dependsOnPsi()`. `meq::PositionAnisotropy` is AN-1's case,
and keeps the linear flux mass bit for bit.

**The flux law.** A `Δ` that reads `ψ` puts **the whole** flux mass on
`mfem::MixedConductionNLFIntegrator( FunctionDiffusionFlux )`, with weight
`R/(1 − Δ)` and derivative `R Δ′/(1 − Δ)²`; there is no linear `M_u`.
`DarcyForm::Reconstruct()` lifts `M_u` whenever it exists and would drop a law
beside it, so `M_u` is never asked for on this branch. Under Picard the
derivative is zero, as the source's reaction term is.

**Refusals.**

* A `ψ`-independent coefficient source is refused; pose it as a `meq::Source`
  with `dFdPsi = 0`.
* A normalised source is refused until AN-2b.
* **An explicit `AssemblyMode::Threaded` is refused, and a default one becomes
  Serial.** `Reconstruct()`'s threaded loop shares one
  `mfem::FrozenDualFluxCoefficient`, whose `Eval()` writes its own members, so
  MEQ's `SetIntegratorsThreadSafe()` promise would be false. Filed upstream as
  `../../mfem-hdg-dev/doc/FROZEN-FLUX-LAW-THREADING-FROM-MEQ.md`.
* **The firehose limit is checked on the solve's last `ψ`.** A solve that ends
  past it is refused. A solve that fails with its last iterate past it names
  the limit rather than the iteration.

**No finite-difference Jacobian check.** The condensed gradient aborts under
`FullNL`, and NPC's is solve-only, so the instrument is the observed Newton
order. Zeroing the derivative is caught: 10 steps in place of 4, best order
1.22.

#### AN-2b — `Δ( Ψ )`, through the labels — DONE

**DONE** → **[M-330](MEASUREMENTS.md#m-330)**.

**The interface.** `Anisotropy::argument()` names what `delta()` reads:
`Position`, `Flux` or `NormalisedFlux`. Under `NormalisedFlux` the solver
hands it `Ψ = (ψ − ψ_bnd)/(ψ_ax − ψ_bnd)`, from the normalised source's labels
as the bordered Newton last set them. It applies `1/(ψ_ax − ψ_bnd)` to
`dΔ/dΨ` itself, in `anisotropyAt()`, the one place that mapping is made.

**The label columns are differenced beside a `Δ` in `Ψ`.** The source's
analytic `normalisationDerivatives()` are `dF/ds`, but the flux row moves with
the labels too. The difference perturbs the labels `anisotropyAt()` reads, so
it carries that term. Measured as a mutation: on the source's analytic column,
Newton is linear, order 1.000, 16 steps in place of 4.

**The axis label at a point, and what it costs `ψ*`.**

* Under the located axis, `ψ_ax`'s irregular error, MEQ's open defect, goes
  straight into the operator.
* With the label at a prescribed point, as a mirror reads its labels, Newton
  is quadratic.
* **Read from `ψ_h`, `ψ*` is not at `k+2`.** The label is `ψ_h` at a point,
  so `Δ` carries `ψ_h`'s pointwise error, `k+1` with a constant that moves
  with where the point sits in its element. That bounds `ψ*`, and it holds for
  any operator coefficient read through the labels.
* **Read from `ψ*`, it is**, → **[M-331](MEASUREMENTS.md#m-331)**.
  `setLabelReading( LabelReading::PostProcessed )` reads each prescribed
  label from `ψ*` by defect correction. Each pass's label row is
  `ψ_label − ψ_h(x) − c = 0`, with `c` the last pass's `(ψ* − ψ_h)(x)`, so the
  bordered Newton within a pass is exactly the one it always was. At the
  fixed point the label is `ψ*(x)`. `ψ*` is back at `k+2` for `k = 1…3`, in 2
  to 6 passes. The corrections are exactly zero by default, so every other
  path is unchanged.

**Refused.**

* `Flux` beside a normalised source, which is not measured.
* `NormalisedFlux` without one, since there are no labels to read.
* **`setPlasmaCurrent()` beside any anisotropy.** The constraint integrates
  `F/R`, and under `Δ`, `Δ*ψ = (∇ψ·∇Δ − F)/(1 − Δ)`, so `F/R` is not the
  toroidal current.

### AN-3 — the mirror, through MEQ

**Split in two, because coil subtraction under `Δ` needs upstream.** The
remainder's flux row carries `−(Δ/(1 − Δ))∇ψ_c` (*What changes* §2, verified
at AN-0). That is an **affine** flux law, `f u + g`, and `g` depends on `ψ`.
MFEM's flux law carries it in the solve. `DarcyForm::Reconstruct()` freezes
the law at zero flux, though, and its local problem has a zero flux
right-hand side. So `g` would reach the solve and not `ψ*`, which would be a
different function and not merely a lower order. Filed as
`../../mfem-hdg-dev/doc/AFFINE-FLUX-LAW-RECONSTRUCTION-FROM-MEQ.md`.

#### AN-3a — the anisotropic mirror with meshed packs

On `RotatingMirror`'s machine, through `MirrorMachine.hpp`, with the packs
meshed, so that there is no load:

* **Lift the extension and DtN refusals** where `Δ` vanishes on `Γ_h`, checked
  on the solve's last `ψ`, since both read `q` there as `∇ψ/R`.
* **Measure `τ`**: the box default against `τ(1 − Δ)`, which AN-1's
  constant-`Δ` identity says keeps the stabilisation in proportion.
* **Answer the jump question** of *Refusals*, at the box's `z` edges, which
  are faces.
* The labels from `ψ*` (AN-2b) on both labels.

**Acceptance**:

* **`Δ = 0` reproduces the meshed-pack mirror to round-off.**
* At `δ₀ = −0.1` and `−0.3`, Newton converges from the vacuum guess, and MEQ
  converges against itself at its design rate (MR-4c's form).

#### AN-3b — the subtracted mirror

**Waits on the request above.** It needs an affine flux law, `g` included in
the reconstruction, and the refusal beside `setConductorField()` lifted.

**Acceptance**: the subtracted machine against AN-3a's meshed one, as
COIL-SUBTRACTION-PLAN's instrument does.

### AN-4 — geq's reference

`tools/geq-benchmark/export_anisotropic.py`, beside `export_rotating.py`. It
follows MR-4a's protocol:

* direct `freeBoundary`;
* `φ₀` tightened;
* `N_D` held;
* `rtol = 2e-9`;
* the labels exported by `psiRZ` at the points.

To that it adds `Δ = δ₀ g( ψ_n ) h( z )` at `δ₀ = −0.1` and `−0.3`, with
`dp_par_dpsin = fast_fraction = 0`, in **both** rearrangements. The profile and
`Δ` tables are tracked; geq's answer is regenerated, not committed.

**Acceptance**:

* geq's own self-convergence at 129², 257² and 513², re-taken.
* The two rearrangements' difference, falling as geq refines.
* `δ₀ = 0` reproduces M-261's tables.

### AN-5 — MEQ against geq

MR-4d's matrix:

* `ψ` at geq's grid nodes, nodes outside `D_h` dropped, geq's four corners
  excluded;
* both labels;
* against 129², 257² and 513².

**Acceptance**: the difference falls at geq's own observed rate. At 513² it
sits within a small multiple of geq's own 257-against-513 difference.

### AN-6 — beyond geq, and through the front end

* **`Δ( B )` through `|q|`**, the mirror's case.
* **`Δ > 0` towards the firehose limit.** MEQ against itself only, since geq
  has no converged reference there.
* **geq's paired `dp_par_dpsin` and `fast_fraction`**, translated into a
  species split and compared.
* **The driver's keys and the `.nc`**: `p_∥`, `p_⊥`, `Δ`, and the firehose and
  mirror-mode margins (geq's `anisotropy_stability`).
* **A race through each code's own front end**: geq's case script and MEQ's
  `meq.toml`, by `race_geq.py`'s protocol.

## Open, before any of it is built

* **What `q`'s consumers assume.** AN-1 begins with the audit.
* **`FullNL`'s cost, and whether `FullNL` runs under a Device.** Neither is
  measured.
* **Whether geq's divergence form converges at `Δ > 0`.** If it does, AN-6's
  firehose leg gets a reference after all.
* **What a prescribed `Δ( ψ, z )` is physically.** Parallel force balance ties
  `∂p_∥/∂s` to `( p_∥ − p_⊥ ) ∂ ln B/∂s`. A `Δ` prescribed independently of
  `B`, with an `h( z )` cut, is a statement of a coefficient and not of a
  distribution. That is enough for a code benchmark, and it is called that.
  A self-consistent `Δ` is AN-6's `Δ( B )` and, beyond it, a kinetic closure.
