# Free-boundary machinery (not built)

This directory holds the free-boundary half of the original meq: the von
Hagenow / Lackner scheme, which computes the boundary data for psi from a
Green's function integral over the plasma and coil currents, together with the
coil-current model and the special quadratures that integral needs.

None of it is compiled. It is not in `CMakeLists.txt` and no target references
it. It is kept because it is the only copy of a working implementation of the
scheme, and because parts of it are worth more than the effort of rewriting
them.

## Why it is not built

Two separate reasons, one of scope and one of API.

**Scope.** meq's first target is the *fixed*-boundary Grad-Shafranov problem:
psi = 0 on a prescribed Gamma. That problem needs no boundary integral at all,
so nothing here is on the path to it. Free boundary is a later question, and
when it is asked it will be asked of a solver built on `DarcyForm`, not on the
one this code was written against.

**API.** The code targets the MFEM branch meq was originally developed on, and
that branch is gone:

* `HDGBilinearForm` no longer exists. The Waterloo HDG interface meq was built
  on was replaced in the 4.9.1 development branch by `DarcyForm`, which
  assembles the same hybridised system through a different set of objects. Every
  file here that assembles anything is written against the old one.
* `GridFunction::GetValueFacet`, which the old error estimator used to evaluate
  the trace unknown on a face, was a patch carried by that branch. It is not in
  MFEM 4.9.1, and there is no drop-in replacement -- the equivalent has to be
  reached through `FaceElementTransformations` instead.

A port is therefore a rewrite of the assembly, not a matter of fixing includes.

## What is here

| File | What it is |
| --- | --- |
| `FreeBoundary.{hpp,cpp}` | The Green's function, the boundary integral `BoundaryPsi` that evaluates psi on Gamma from the normal flux there, `GreensFunctionPsi` for the coil contribution, and two `mfem::Coefficient` subclasses that present the result to a linear form. |
| `BoundaryTraceIntegrators.{hpp,cpp}` | `HDGBoundaryTraceIntegrator` and `HDGBoundaryNormalTraceIntegrator`: the linear form integrators that put the computed boundary data into the HDG system. |
| `utility.hpp` | `Coil`, `Jtor` (the sum of the coil currents at a point), and the beginnings of a domain/boundary-condition description. |
| `GaussLog.hpp`, `GaussLog10.hpp` | 5- and 10-point extended Gauss-log quadrature rules -- nodes and weights for integrands with a logarithmic endpoint singularity, which is what the Green's function has when the evaluation point lies on the face being integrated over. |
| `VacuumGFSoln.cpp`, `VacuumMeshSoln.cpp` | Standalone drivers that compute the exact vacuum field by integrating the Green's function over the coil rectangles with Gauss-Kronrod quadrature, for comparison against a solve. |
| `vacuum-testing/` | Configurations, meshes and reference output for those two drivers. |

`BoundaryTraceIntegrators.{hpp,cpp}` did not exist under those names in the
original tree. They were extracted from `HDGGSIntegrator.{hpp,cpp}` -- which
also contained the fixed-boundary integrators, and which was deleted in the
move to `DarcyForm` -- so that the free-boundary parts survived the deletion.
The rest of that file, and everything else removed in the reorganisation, is
recoverable from the `v0-legacy` tag:

    git show v0-legacy:HDGGSIntegrator.cpp

## Worth keeping

`GreensFunction` in `FreeBoundary.cpp` is the reason this directory exists
rather than a note saying the code was deleted. The Green's function for the
Grad-Shafranov operator is written in terms of the complete elliptic integrals
K(k) and E(k) with

    k^2 = 4 R R* / ( (R + R*)^2 + (Z - Z*)^2 )

and k -> 1 as the source point approaches the field point, where K has a
logarithmic singularity. Evaluating the standard expression there loses all
precision: `1 - k^2` is a difference of nearby quantities, and by the time it is
small enough to matter it has no significant digits left. The function has a
separate branch for that regime which expands K in the small parameters
`dR = (R* - R)/R` and `dZ = (Z* - Z)/R` -- keeping the linear terms, which
doubles can still resolve, and dropping the quadratic ones, which they cannot --
and takes the logarithm in a scaled form chosen to avoid underflow when one of
dR, dZ is much smaller than the other. That analysis is the expensive part of
this code, and it is independent of which MFEM API assembles the system.

The singularity handling in `BoundaryPsi` is the same kind of content: it
detects whether the evaluation point lies on the face currently being
integrated, and if so splits the face at that point and applies the Gauss-log
rule to each half with the singularity moved to the endpoint, rather than
letting an ordinary Gauss rule integrate a logarithm.

## Known defect, to fix if this is ever ported

`FreeBoundary.cpp`, in `GreensFunctionBoundaryCoefficient::Eval`:

    throw new std::logic_error( "This only works with Segements!" );

This throws a *pointer* to an exception, not an exception. `catch
(std::exception&)` will never catch it, so it propagates past every handler that
looks correct and terminates the program; if something does catch it as
`std::logic_error*`, the object leaks unless that handler deletes it. Drop the
`new`. The two other throws in the same file are already correct, which is what
makes this one easy to miss.

## If this is revived, do not revive the algorithm

Decided 2026-08-24; the reasoning is in `refs/Refs.md`, section *Free boundary*.

This code implements Lackner's 1976 scheme, and there are two separable things in
it. The **reduction of the unbounded domain** through the analytic Green's
function is sound and is what any free-boundary solver does, CEDRES++ included.
The **outer fixed-point iteration it is embedded in** is what meq is moving away
from, and it is the reason the code here has the shape it has.

`GreensFunctionBoundaryCoefficient::Eval` calls `BoundaryPsi`, which sweeps every
boundary face with singular quadrature — once per quadrature point, so `O(N²)` —
and the whole sweep is repeated at every Picard iteration, because in that
arrangement the Green's function is an explicit boundary-condition *update*
computed inside the loop.

meq uses Newton. The corresponding structure is a **coupled BIM/FEM system**, in
which the boundary integral operators sit inside the operator Newton
differentiates, rather than in an outer loop around it. Gatica & Hsiao's
uncoupling makes that cheap by taking the artificial coupling boundary to be a
circle, so the boundary integral operators invert exactly and a single weakly
singular term survives.

So the parts of this directory worth carrying forward are the quadratures and the
Green's function evaluation — the numerics listed under *Worth keeping* — and not
the control flow around them.
