# `ReconstructFluxAndPot()` is threaded — and one thing you have to call

Follow-up to `HDG-RECONSTRUCT-THREADING-REPLY-FROM-HDGDEV.md`, whose last
section says "the reconstruction routines are serial today and are covered by
no audit behind the flag". Both halves of that have changed.

## What you have to do

**Call `SetIntegratorsThreadSafe()` as well as `SetAssemblyMode(Threaded)`.**
Without it this routine stays serial, whatever the mode. That is deliberate
and it is the thing your own report asked for: the loop evaluates two LINEAR
face integrators — `c_bfi`'s `AssembleFaceMatrix()` and `c_bfi_p`'s
`AssembleHDGFaceMatrix()` — which `MultNL()`'s audit behind that flag never
covered, so a caller who promised for that loop has not promised for this one.
Making the flag license it silently is exactly the asymmetry you filed against.

The whole reach of the promise is now written on the setter, loop by loop, so
you can check what you are undertaking rather than infer it.

## What it is worth

Measured here with `MKL_NUM_THREADS=1` and `OMP_WAIT_POLICY=passive`, speedup
of the whole routine against the serial mode:

| | nt=1 | nt=2 | nt=4 | nt=8 |
|---|---|---|---|---|
| 2-D 64² quads, k=2, 4096 el | 0.99 | 1.74 | 2.63 | **4.65x** |
| 2-D 96² quads, k=1, 9216 el | 0.92 | 1.79 | 2.29 | **3.47x** |
| 3-D 12³ hexes, k=2, 1728 el | 0.98 | 1.61 | 3.03 | **4.54x** |

Against your 0.620 s of a 3.302 s run at 1.00 cores, 4.5x on that leg is about
**1.17x on the run** — a little more than the 1.10x the earlier note predicted,
because the prediction assumed it would thread no better than the residual leg.

**The one-thread column is a tax, not noise**, 0.92x to 0.99x, and it is the
same tax you measured independently on the NPC traversal (0.95x on the leg,
0.98x on the run). Entering the region and building the per-thread workspace is
not free. The k=1 row scales worst because the serial replay below is a larger
share of a cheaper solve.

## Why it is not simply a parallel loop, which is the part worth reading

**The enriched trace is whichever element visited the face last.** Each element
solves a local problem whose trace is free on *every* face, so the two elements
either side of an interior face produce two different values for that face's
trace dofs — and `tr.SetSubVector()` assigns. Nothing in the tree said so.
Measured by reversing the element order on one fixed solution:

| | forwards | backwards |
|---|---|---|
| `‖u*‖` | 7.9321522926801267 | 7.9321522926801267 |
| `‖p*‖` | 5.7077841841727244 | 5.7077841841727244 |
| `‖tr*‖` | 5.4581137616313953 | **5.3635341954733962** |

So **no colouring reproduces the serial answer**: a colouring changes which
element writes last, which moves the trace by 2% of its own norm. The routine
therefore SOLVES in parallel and REPLAYS the writes serially in element order,
a bounded chunk of local solutions at a time. Bit for bit against serial at
every thread count — `max|Δ| == 0.0e+00` on the total flux, the enriched flux,
the enriched potential and the enriched trace.

If you would rather have the extra ~15% of the leg than bit-for-bit agreement,
say so and the replay can become a colouring behind an opt-in; it is a
one-line change to the arrangement and a large one to the contract.

## Two things you will want to know

**Your `HDGExtensionIntegrator` reasoning still holds and is now moot anyway.**
The loop still takes the flux mass from `M_u->GetDBFI()`, domain integrators
only, so the boundary-face extension integrator is not reached. Its eleven
scratch members are behind `#ifndef MFEM_THREAD_SAFE` on
`gf-hdg-subdomains-dev` regardless. Note for your build: that guard **does**
move the layout where `MFEM_THREAD_SAFE` is on, so a clean rebuild is owed.

**One hazard your inventory found and we mis-scoped.** You flagged
`FiniteElementSpace::GetElementDofs(int, Array<int>&)` writing the space's
shared `mutable DofTransformation`, and our note filed it as "benign by
accident of the space" because for L2 the only write is a null pointer. That
accident does not hold here: this routine also serves an **RT** flux, and
`DofTransformation::SetFaceOrientations()` is `Fo_ = Fo`, an `Array<int>` copy
assignment — so two threads gathering dofs from one space can be allocating and
freeing that member's buffer at once. Heap corruption, not a wrong answer. The
loop takes the three-argument overload with per-thread storage.

It was found by reading the class **after** the RT test had already passed at
eight threads, twice. Worth repeating wherever you are threading your own
loops: a green threaded run is evidence about one interleaving.

`ReconstructTotalFlux()`, the other half of `Reconstruct()`, is **still
serial**. It is the smaller half and nobody has taken the split; say if it
matters to you.
