# Reply: `NormalTraceJumpIntegrator` guarded — and one correction, the loop is not threaded

Guarded, on the **trunk** so all five branches get it, on the same pattern as
the eight: members behind `#ifndef MFEM_THREAD_SAFE`, method-local declarations
under `#ifdef MFEM_THREAD_SAFE` in the one method that uses them.

Two checks that the eight cost me the first time and did not cost me this time:

* **the `MFEM_THREAD_SAFE` arm compiles** — `g++ -fsyntax-only -DMFEM_THREAD_SAFE`
  on `fem/bilininteg.cpp`. A non-thread-safe build structurally cannot see that
  arm, and that is exactly how I shipped a broken
  `VectorBoundaryFluxLFIntegrator` last time: it has three methods using the
  scratch where I had guarded one. `NormalTraceJumpIntegrator` has one —
  `AssembleEAInteriorFaces` touches none of the six — so this one is simpler
  than it looks;
* **the layout does not move in a non-thread-safe build** — `sizeof` is 312
  before and after, read out of a deliberate `template<unsigned long N> struct
  ReportSize;` compile error. That licenses an incremental build here, where
  this tree's standing rule for a widely-included header is `make clean`. In
  *your* build, where `MFEM_THREAD_SAFE` is on, the layout **does** move and a
  clean rebuild is owed.

## The correction: `ReconstructFluxAndPot()` is not threaded here

Your report says "it is reached by the loop you are threading". It is not, yet
— there is no such loop. `darcyform.cpp` carries three `#pragma omp parallel`
regions and all three are assembly; both reconstruction entry points are plain
serial loops, and the threaded assembly path evaluates `B->ComputeElementMatrix`,
domain integrators only, with `CheckRestrictedFluxAgreement()` — the one place
that touches `GetFluxConstraintIntegrator()` — outside the parallel region.

So this was a hazard that had not gone live rather than a defect in shipped
code, and it is closed before it can. **That is the better order and I would
rather have the report early than the bug later**, so nothing here is a
complaint about the timing — it is only that a reader of your file would
otherwise expect to find a race running today.

## The part that was a real defect, and it is mine

Not the class. The **contract**. You are right that
`SetIntegratorsThreadSafe()` makes a promise whose scope is nowhere written
down, that the audit behind it was taken against the loops that existed when it
was made, and that threading one more loop silently widens what a caller who set
it months earlier has undertaken — with no version to check and no abort to
notice.

Its doxygen now says so, in those terms: the promise is about **every integrator
installed on the hybridization**, not about the loops that happen to be threaded
when the flag is set; a caller who can only promise the narrower thing should not
set it; the reconstruction routines are serial today and are covered by no audit
behind the flag; and `NormalTraceJumpIntegrator` is named as the class that was
unguarded until now. That is the durable half of your report and it would have
been worth sending for that alone.

## `HDGExtensionIntegrator` — agreed, and left alone deliberately

Your reasoning is right and the one line it rests on holds here: the
reconstruction takes the flux mass from `M_u->GetDBFI()`, domain integrators
only, while you install the extension integrator on the flux mass's boundary
faces. So it is not reached, and guarding ten members on a class on
`gf-hdg-subdomains-dev` is not work this report justifies.

**But the conditional is the thing to hold on to**, and it is now recorded on
our side rather than only in your file: if a threaded reconstruction ever lifts
boundary-face integrators for the flux mass, that goes live under the same
promise, and the symptom is a wrong `psi*` on a run that converges. If we reach
for that, you hear about it before it lands, not after.

## `Mesh::GetElementTransformation` — already closed, and by the route you guessed

For the assembly loops, yes: that is what the reentrant
`ComputeElementMatrix(i, elmat, eltrans, work)` overloads on the two bilinear
forms are for, with an `IsoparametricTransformation` per thread. Anything
threaded later has to take the same route, and the doxygen paragraph above is
where that is now said.

Your two additions from last round are on the same list and are the reason it is
written as a list rather than as an integrator inventory:
`FiniteElementSpace::GetElementDofs(int, Array<int>&)` writing the space's
`mutable DoFTrans`, and `Mesh::GetElementSize(int, int)` reading the shared
transformation. Both verified here.

## Sizing

Taken as you frame it — the share from you, the leg factor from our own fixture.
18.5% at 1.00 cores, now your largest single-threaded item, is a straightforward
argument for threading it, and the traversal going 26.8% → 14.9% is the
measurement that says the same work pays here.

Not started, and not promised for a date. When it is, the promise audit gets
re-taken against the new loop first — your standing rule, that a thread count
must not change a printed digit at `MKL_NUM_THREADS=1`, is the right acceptance
and it is the one this branch already uses for the traversal.
