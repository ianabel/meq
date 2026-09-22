# How to drive NICE

NICE is a free-boundary equilibrium code, C++ and P1 finite elements on an
unstructured triangulation, with a **direct Newton** and an `I_p` Lagrange
multiplier. It is the fastest code in MEQ's comparison set.

**Tree:** `/home/ian/projects/nice`  ·  **Harness:** `meq-race/`

## The one command

```sh
/home/ian/projects/nice/meq-race/race.sh limited 1.0     # 1933 P1 dofs
/home/ian/projects/nice/meq-race/race.sh limited 0.25    # 6674
/home/ian/projects/nice/meq-race/race.sh limited 0.0625  # 26166
/home/ian/projects/nice/meq-race/race.sh diiid   1.0
```

Each run is from scratch: the script regenerates the deck, wipes `input/`,
`output/` and `restart/`, and never warm-starts. The second argument scales
**every** `*TriangleArea` target, so it is the resolution knob and areas go as
`h²` — `0.25` is `h/2`.

## Things that will cost you an afternoon

**`MKL_THREADING_LAYER=GNU` IS REQUIRED.** It is set inside `race.sh`. Without
it UMFPACK returns silent garbage and the Newton produces `-nan`. Debian's
`libblas.so.3` is an alternatives symlink to `libmkl_rt.so` on this machine, so
NICE gets MKL whether or not anybody asked for it.

**THE OpenMP BUILD IS PREFERRED AND BUYS NOTHING, AND BOTH HALVES MATTER.**
The stock `build/CMakeCache.txt` reads `CHOOSE_OPENMP:BOOL=OFF` while `src/`
carries **270 `#pragma omp`** directives, which looks like the benchmark denying
NICE its parallelism. It measures as worth nothing: `build-omp/` is within
**1%** of the serial build at every resolution while burning 2.6 to 3.3 cores,
and the *non*-OpenMP build is actively **slower** with eight MKL threads (1.51 s
against 1.37 s at 6674 dofs). M-164 section 1.

`race.sh` therefore prefers `build-omp/nice_dir` when it exists and falls back to
`build/nice_dir`, which answers identically; `NICE_BIN` overrides both. NICE is
represented with its parallelism *available* rather than compiled out, and the
conclusion does not depend on which is picked.

**BUT AN OpenMP BUILD MUST BE GIVEN A THREAD COUNT, OR IT TAKES EVERY LOGICAL
CPU.** Unset, OpenMP grabs all 16 here on 8 physical cores, SMT siblings share
an FPU, and the solve goes from **108 ms to 387 ms** — so switching to the
OpenMP build without pinning makes NICE look 3.6× worse than it is. `race.sh`
defaults to the physical core count, which is also what MEQ is given;
`NICE_THREADS` overrides it. At 6674 dofs it reads 412 ms at one thread against
414 ms at eight, which is the finding above restated from the harness.

To build the OpenMP one, both flags are needed — the second is not optional and the
failure is `fatal error: Sparse: No such file or directory`:

```sh
cmake -B build-omp -DCMAKE_BUILD_TYPE=Release -DCHOOSE_OPENMP=ON \
      -DEIGEN3_INCLUDE_DIR=/usr/include/eigen3 \
      -DCMAKE_CXX_COMPILER=g++ -DCMAKE_Fortran_COMPILER=gfortran
cmake --build build-omp -j6
```

`-DEIGEN_DONT_PARALLELIZE` is unconditional in NICE's own `CMakeLists.txt`.
Leave it: it is the author's choice, and `CHOOSE_OPENMP` switches on NICE's own
pragma regions, not Eigen's.

**COUNT P1 DOFS, NOT ELEMENTS**, when matching resolution against MEQ:
`nodes − nodes on r = 0 + 1`. `psi = 0` is imposed on the axis and the `+1` is
the `I_p` Lagrange multiplier. `race.sh` prints it.

**SWEEP THE WHOLE `*TriangleArea` FAMILY TOGETHER.** The semicircular ABB
boundary's segment length is `sqrt( 2 · airTriangleArea )`, so the accuracy of
the boundary-integral (DtN) condition is set by the **air** knob while the
plasma is set by the **vacuum** one. Moving only one is a half-refinement.

## What NICE is NOT doing that MEQ is

**NICE FINDS ITS OWN LIMITER CONTACT.** It takes the limiter as a polygon
(`limiter.txt`) and puts `psi_bnd` at the contact it locates — on the limited
case, `( 0.7359, 0.2296 )`, which is on the true circle of radius 0.35. MEQ's
shipped file is handed a point instead, and that point is `freegs4e`'s own
staircase grid-node artefact. **The two are not solving the same boundary
condition** and the difference is most of the 1–3% between them. See M-164
section 3 before comparing any `psi_bnd`.

`λ`, which NICE solves for as `I_p / I_p( λ = 1 )`, comes out within **0.3% of
1** on both cases. That is the free check that the profile conversion — the flux
span, the `r0` gauge and the `1/μ0` — is right: an error there shows as `λ` off
by that factor while the equilibrium still converges.
