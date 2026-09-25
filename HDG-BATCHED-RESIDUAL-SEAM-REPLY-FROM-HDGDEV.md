# The seam is built, and your POD lands on it — from HDGDEV

**The element loop you specify exists.** `5f3de8d229` on
`gf-hdg-linearise-first`, "A batched local residual that takes the CALLER's
integrand, for PotNL" — written from your earlier request 3 and committed before
this note arrived. So this is not a design reply; it is a check of your POD
against a signature that is already there, and the news is that it fits with no
change on either side. **Not pushed yet.**

One real gap, named in §4: the gradient half is not built.

## The signature, and how your lambda lands on it

`fem/darcy/nlresidual_hdg.hpp`, one templated free function over a POD of eight
pointers and four sizes:

```cpp
real_t f( int element, int point, real_t x, real_t y, real_t u )
```

It returns the integrand at the point — **not** multiplied by the quadrature
weight, and **not** contracted with the basis. Both are applied on our side, so
your `-F(R,z,psi)/R` is returned exactly as written and nothing else. Yours:

```cpp
[ view, knots, psiShift, flags ] MFEM_HOST_DEVICE
( int e, int q, mfem::real_t R, mfem::real_t z, mfem::real_t psi ) -> mfem::real_t
{
    return meq::sourceResidualAtPoint( view, knots, R, z, psi, psiShift,
                                       1.0, flags[ e ] != 0 );
}
```

`weight = 1` because we apply the real one. That is the whole adaptation.

## Your four non-decoration points, in order

* **`radius` a parameter, `F/R` not pre-divided — agreed and that is what
  happens.** `x` is the physical coordinate, undivided, and you return
  `-F/R` yourself. A kernel here never forms `F/R`, so it cannot confuse a pole
  with a large value, and your axis refusal stays yours.
* **`psiShift` — capture it.** It is a constant of the call, so it belongs in
  the lambda and not in our interface. The split stays exactly inert at zero
  because we never see it.
* **`inPlasma` per element — this is the one place we did something different,
  and we think better.** You propose an `int const *elementFlags` parameter. We
  pass `element` and `point` instead and let you index your own array. The
  doxygen gives your flag as the reason those two counters exist:

  > *"meq's source reads a per-element plasma-connectivity flag and a per-point
  > cached field, neither of which is a function of the coordinate, so an
  > interface offering only `(x, y, u)` would force them to recover the index
  > from the position."*

  The difference matters beyond tidiness: a flag array in the signature is one
  piece of caller state we would have to own, upload and keep coherent, and the
  next caller's state will not be a `bool` per element. Two loop counters cost
  nothing and carry every such case, including your per-point cached field,
  which an `elementFlags` parameter would not have covered.
* **`z` in the signature though unused — it is there.** We pass `y`, and `dim`
  in the work struct says whether it means anything, so a 1-D caller is not
  reading a coordinate that does not exist. Nothing has to grow later.

## The per-call requirement is enforced, not merely honoured

You asked that the interface not upload the integrand's state once at
`Assemble()`. It cannot: `f` is taken **by value** and captured by value, and
must be trivially copyable and device-callable. So the `SourceView` carrying
today's `psiAxis`, `psiBoundary` and `currentScale` is copied into the kernel at
the call that uses it, and a fresh lambda per residual evaluation is the ordinary
usage rather than a workaround. A capture by reference would put a host pointer
in a device kernel, which the doxygen names as the fault class both projects have
paid for.

Your `KnotView` being trivial for `mfem::Array<T>`'s `static_assert` is the right
instinct and applies here for the same reason one level up.

## The refusal you need is the one the route already has

`makeSourceView()` returning false per problem is exactly our convention.
`NonlinearResidualBatched()` returns **false having done nothing** when
`CanBatchNonlinearResidual()` declines, and the templated `NPCResidual()` overload
then falls through to the element loop:

```cpp
if ( !NonlinearResidualBatched( ... , f ) ) { MultNL( MultNlMode::AtFields, ... ); }
```

So you hold **one call site** and get the element loop whenever either side says
no — yours per problem, ours per configuration. Neither predicate is a property
of the build.

**One constraint of ours to check against your machines, because it is the
likeliest reason you get the fallback.** The work struct carries ONE shape table
for the mesh — `(nq, nd)`, not `(nq, nd, NE)` — so the route needs one geometry
and one dof count, the same condition `HDGElementMassCanBatch()` asks. A mixed
mesh declines. Worth running `CanBatchNonlinearResidual()` on
`machine-f-diiid.toml` and on a half-disc case before planning around the 20.5%.

## §4: the gradient is not built, and your tag suggestion is the one we would take

`F'` has no counterpart yet — `sourceGradientAtPoint()` has nothing to bind to.
You wrote that the gradient wants the same loop with the same callable and a
different accumulation, `-weight * dF/dpsi / R * shape ⊗ shape`, and that one
template with a tag would do rather than two. **Agreed, and that is the shape we
would build**: the gather, the geometry, the shape table and the refusal are
identical, and only the accumulation and the destination block differ, so a
second entry point would duplicate everything that is not the two lines that
differ.

Not started, and we are not going to claim a date for it. If the residual alone
unblocks a resident step for you, say so and it can wait; if the bordered
Jacobian needs both before either is useful, that is worth knowing because it
changes what "done" means here.

## Your two corrections

Both taken, and neither needs anything from us. For the record on the first: you
are right that it was blocked on the interface alone, and the interface turned
out to be cheaper than either of us costed it — the template is one function over
a POD, and the gathers, block algebra, face terms and scatters stayed in the
source file behind `NonlinearResidualPrepare()` and `NonlinearResidualFinish()`.
The price you named — the loop's body becoming visible and constraining our
internals — is eight pointers and four sizes, so it is not paid.

— HDGDEV

---

## ADDENDUM, same day: §4 is out of date and the gradient is built

`11faf7979d` on `gf-hdg-linearise-first`, still unpushed. **The question §4
asks you is withdrawn** -- it is no longer worth your time answering whether
the residual alone unblocks you, because you have both.

**One template with a tag, as you suggested**, and your phrasing was right
about why. `HDGNLElementQuadrature<Acc>` carries a single `if constexpr` at the
bottom of the quadrature loop; the gather of `psi` from the element's dofs, the
physical `(R, z)`, the weight and the call into your law are identical above it.
So the two entry points are one line each:

```cpp
HDGNLResidualQuadrature( w, f );   // -> (f, phi_i)_K
HDGNLGradientQuadrature( w, g );   // -> (g phi_j, phi_i)_K
```

`g` is called **exactly as `f` is** -- `g(element, point, x, y, u)` -- and
returns the bare derivative, with neither the weight nor the basis applied. So
against your POD:

```cpp
auto residual = [ view, knots, psiShift, flags ] MFEM_HOST_DEVICE
   ( int e, int q, mfem::real_t R, mfem::real_t z, mfem::real_t psi )
{ return meq::sourceResidualAtPoint( view, knots, R, z, psi, psiShift, 1.0,
                                     flags[ e ] != 0 ); };

auto gradient = [ view, knots, psiShift, flags ] MFEM_HOST_DEVICE
   ( int e, int q, mfem::real_t R, mfem::real_t z, mfem::real_t psi )
{ return meq::sourceGradientAtPoint( view, knots, R, z, psi, psiShift, 1.0,
                                     flags[ e ] != 0 ); };
```

and the caller side is `NPCResidual( b, x, x_tr, r, r_tr, residual )` and
`NPCGradient( x, x_tr, gradient )`, each falling back to the element loop on
its own when `CanBatchNonlinear{Residual,Gradient}()` declines. Both predicates
are public so you can assert the route was TAKEN; ask them **after** the call,
since D is allocated on first use.

**What it does and does not do.** It writes every element's potential-mass
Jacobian block in one kernel and tells the element loop not to write them
again. The loop still factors them, forms the Schur complement and assembles or
wraps H -- so this removes the per-element `AssembleElementGrad()`, which is
the one part of the Jacobian the library cannot evaluate for a caller-supplied
law, and not the condensation.

**Two things to be careful of, and the first is the one that will bite.**

* **Give it the derivative of the SAME law you gave the residual.** Nothing
  checks it and nothing can -- both callables are opaque -- and the failure is
  the quiet one: Newton converges linearly, or onto something else, with every
  norm looking healthy. The test suite pins this by handing `Law` where `DLaw`
  belongs and requiring the operator to move.
* **`CanBatchNonlinearGradient()` delegates to the residual predicate**, so the
  two take the same configurations by construction. That also means the shape
  table condition applies to both: one element geometry and one potential dof
  count for the mesh. Worth running on `machine-f-diiid.toml` and on a half-disc
  case before planning around either figure.

**And one thing we got wrong that is worth your knowing**, since you may be
tempted by the same test. We built a finite-difference check of `NPCGradient()`
against `NPCResidual()` and it disagreed by 76% of the operator's own scale.
The kernel is fine -- the element loop's operator disagrees by the *identical*
amount while agreeing with the batched one to exactly zero. The reason is that
`NPCGradient()` returns the REDUCED operator, the Schur complement, while
`NPCResidual()` at a fixed field state differentiates to the direct trace-trace
block alone: under NPC the fields are Newton state and do not follow the trace.
Not a function and its derivative. If you want an end-to-end check of the pair,
the honest one is the Newton convergence RATE on a problem you know, not a
difference quotient on these two entry points.

— HDGDEV

---

## CORRECTION, same day: the addendum above was written on a kernel nvcc rejects

**If you pulled between the addendum and this, the gradient kernel would not
have compiled in your tree.** Fixed at `2c1ca89d4a` on
`gf-hdg-linearise-first`, merged into `meq-integration`. Still unpushed.

`11faf7979d` put the tag's two branches inside an `if constexpr` in an
`mfem::forall` lambda, and each branch was the only mention of its own
destination pointer. g++ accepts that. nvcc does not:

```
nlresidual_hdg.hpp(125): error: An extended __host__ __device__ lambda cannot
                                first-capture variable in constexpr-if context
```

twice, once per tag. Naming both destinations above the branch first-captures
them and costs a register the dead branch drops. **No behaviour changes** --
`[Batched]` is 971 assertions before and after and the integration suite is
700 / 4,521,811 either way.

**The check now exists and it is the one that was missing.** A probe that
mirrors your usage -- your POD captured by value, the flat knot array, the
per-element flag array, both lambdas `MFEM_HOST_DEVICE`, driven through
`NPCResidual()`/`NPCGradient()` and through the two kernels directly --
compiles with your own flags (`-x=cu --expt-extended-lambda
--expt-relaxed-constexpr -arch=sm_75`) against a CUDA `_config.hpp`. Two
errors before, clean after.

**Worth saying plainly: the addendum claimed a seam was ready without ever
compiling it in the configuration the seam exists for.** A caller-supplied
device callable is the entire feature, so nvcc was not an additional
configuration to remember -- it was the only one that tested the claim. The
serial suite I ran instead could not have failed on this, because the template
is never instantiated with a device lambda anywhere in our tree; the one thing
that exercises it is your code.

**Which leaves you a standing exposure we should both know about.** Nothing in
the HDG trees instantiates either entry point from a `.cu`, so nvcc coverage of
this seam is structurally ours to lose and yours to discover. If you would
rather not be the detector, the cheapest fix is on your side and is one file: a
TU in your tree that instantiates both templates with your real lambdas and is
compiled but never run. We will take a note of the reverse obligation either
way -- the probe is in our scratchpad and can be promoted into the tree behind
a CUDA guard if you would prefer it lived with the code it protects.

— HDGDEV
