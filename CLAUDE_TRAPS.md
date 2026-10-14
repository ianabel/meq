# Traps: the build, the code, the shell, and a red test

**The long form of `CLAUDE.md`'s *Traps* checklist** — each entry with the measurement and the story of how it hid, which is the half worth keeping. Traps specific to a subject live with it: the solver's in `CLAUDE_HDGGS.md`, MFEM's API and build in `CLAUDE_MFEM.md`, the device in `CLAUDE_DEVICE.md`, the configuration parser in `CLAUDE_DRIVER.md`.

**The standing rules are `CLAUDE.md`'s.** Measurement tables live in `MEASUREMENTS.md` under stable `M-nn` anchors, which are not renumbered. A hypothesis that was measured and falsified stays; the history of the file does not — present tense, no "used to", no dated "now fixed". What was decided against is `CLAUDE_DECLINED.md`.

## Reading a red test: the fixture, and the coin toss

**WHERE A CASE LIKE `theTwoBorderSolveReportsATrueMagneticAxis` GOES RED, SUSPECT
THE FIXTURE BEFORE THE SOLVER.** `examples/limited-tokamak.toml` has **every** ingredient the case
has — a domain reaching `R = 0`, `[boundary.limiter]`, `[boundary.exterior]` and
`ConfineToPlasma` — and converges in **11 Newton steps** to a `Ψ` of 1.0000,
reproducing freegs4e to 1.3e-04. What it had that the fixture did not is **coils
and a prescribed current**. The fixture now has both, plus the confinement, and
it is **green at five limiter radii in a row**, 1.08 to 1.18, with the coil-free
control converging to the wrong topology at three of them; see `CLAUDE_FB.md`,
*§7.12b's fixture was the defect*, and `FREE-BOUNDARY-PLAN.md` §11.3 for the diagnosis and §11.7 for
the repair.

Neither ingredient alone repairs the fixture, and confinement with a prescribed
current does not either; it took the full 2 x 2 to say so
([M-134](MEASUREMENTS.md#m-134)).
→ **[CLAUDE_DECLINED.md](CLAUDE_DECLINED.md)**, *Confinement, clamping and a prescribed current: the M-134 two-by-two*

**IT MUST NOT ASSERT `checkAxis().agrees`, WHICH IS A TAUTOLOGY UNDER THE
LOCATED-AXIS CONSTRAINT**: `Ψ` there reads 1 by construction — measured, that
very configuration reports `Ψ = 1.0000 AGREES` on a solve whose axis source is a
pole. A test that cannot fail is worse than no test. It asserts four things, and
`checkAxisSource().bounded` is the one with teeth.

**AND IT HAS NOW HAPPENED TWICE, ON A CASE THAT HAD BEEN RED FOR SESSIONS.**
`XPointOuter` was failing at its bootstrap solve with a long, careful list of
ruled-out causes and two surviving suspects, both of them about the solver. It
was neither: the fixture called `setPlasmaSupport()` on the plasma source
**before** constructing the coil wrapper, so the wrapper's own flag stayed false,
`plasmaComponentWanted()` read false, and **XP-1's flood fill never ran on the
only diverted case in the tree** — 333 elements of private flux region, 44% of
the candidates, carrying a current channel nobody asked for. `apps/meq.cpp` has
always done it the other way round and is correct. See `CLAUDE_FB.md`,
*`setPlasmaSupport()` only reaches the fill if it is called on the wrapper*, and
→ **[M-82](MEASUREMENTS.md#m-82)**.

**THE PART WORTH GENERALISING IS HOW IT HID.** Every experiment on that case
changed one key and read the outcome, which is the right instinct and is exactly
what could not find this: the fill was off in **both** arms of every one of them,
so a real measurement — "turning `ConfineToPlasma` off makes it converge" —
supported a conclusion that was wrong. It took a **2 × 2 cross** to see that the
fill and the frozen plasma edge are each necessary and neither sufficient. *A
one-key experiment separates two hypotheses only if everything else is where you
think it is.*

Per *Testing stance*, a red suite is the intended signal while a defect stands.
**Do not "fix" it by relaxing an assertion** — the repair in both cases above was
to give a fixture something it was missing, and the assertions never moved.

**AND A THIRD CASE WENT RED WHERE HEAD WAS GREEN, AND IT WAS NEITHER THE FIXTURE
NOR THE SOLVER — IT WAS A COIN TOSS, AND `HEAD` LOSES IT TOO.**
`XPointBorder`'s `theXPointMeritWeightChangesTheWorkAndNotTheAnswer` failed at
`XPointMeritWeight = 4` with *"weight 4 did not solve"* while a `HEAD` worktree
passed it in 81 Newton iterations against a cap of 200. **A worktree comparison
says *this is not how the tree was*; it does not say *this is a regression***,
and the step that separates them is a **numerically neutral perturbation applied
to BOTH binaries**. `MKL_NUM_THREADS` is that perturbation — it reassociates a
blocked BLAS-3 and changes nothing else. Bootstrap iterations at weight 4,
against the file's cap of 200:

| | MKL=1 | MKL=2 | MKL=4 | MKL=8 | MKL=16 |
|---|---|---|---|---|---|
| `HEAD` | 81 | 65 | 115 | **over the cap** | **over the cap** |
| working tree | **over the cap** | 84 | 75 | 103 | 103 |

**`HEAD` FAILS THE SAME CASE IN TWO OF FIVE CONFIGURATIONS**, and every arm that
converges reports `psi_ax` to eleven digits and the X-point to six decimals —
the invariant the case exists for, which never moved. `MKL_NUM_THREADS=1` — what
`tests/CMakeLists.txt` pins, and the only value ctest ever runs — is one draw,
not a property of a tree. → **[M-186](MEASUREMENTS.md#m-186)**.

**The repair is the same shape as `PedestalConvergence`'s and the assertions did
not move.** The case took `BOOST_TEST_REQUIRE( got.sweeps > 0 )` per weight,
which is a **precondition** on a chaotic plateau — the residual sits near
`1.8e-04` for sixty iterations, RISES to `4.3e-04`, then falls to `1e-12` in
five — and in the more brittle direction, since a `REQUIRE` abandons the weights
after it as well. It now sweeps four weights either side of one, prints and
skips a weight that does not solve inside the cap, and requires three of the four
to compare. **What was relaxed is a precondition, not a tolerance**, and that
distinction is the whole of what makes this a repair rather than the thing this
section forbids. **That the majority is a requirement and not a second coin toss
had to be measured before it could be asserted**: the other three weights read
`13, 13, 14` iterations at every `MKL_NUM_THREADS` tried, none of them failing,
so weight 4 is the only one near the plateau and at most one dropout is
reachable.

**AND THE ROUND-OFF HAD NO SEMANTIC CAUSE TO FIND, WHICH IS WHY BISECTING WOULD
HAVE BURNED THE SESSION.** Diffing the two binaries' full transcripts, 192 lines
of 377 differ and **every one is a residual at `1e-12` or below** differing in
its fifth figure — `1.5e-17` absolute against a problem scale of `7e-2`. No
iteration count, no `psi_ax`, no support count differs anywhere else in the
binary. The tolerances that had been rewritten as
`std::max( <the double value>, <something following real_t's epsilon> )` each
**are** their historical value in double, checked one at a time. What moves the
last bits is the compilation: a literal moved behind a function call or a
`std::max` changes which FMA contractions `g++` emits at `-O3`. *Identity of a
constant is not identity of the emitted code, and on a plateau that is enough.*

## Provenance: what CMake reads once

**AND ANYTHING A CMakeLists READS WITH `execute_process()` IS READ ONCE, WHEN
THE BUILD DIRECTORY IS CONFIGURED — WHICH IS HOW A BINARY STOPS KNOWING WHICH
TREE IT CAME FROM** → **[M-224](MEASUREMENTS.md#m-224)**. The commit was
captured that way and `configure_file`d into `meq_version.hpp`, so a long-lived
build directory stamped **every** binary it produced, and every `.nc` those
binaries wrote, with whatever `HEAD` was the day somebody first configured
it. Measured at five commits stale, and the string is *plausible*
— a real commit of this repository, merely not this one.

**IT IS A PROVENANCE DEFECT AND THREE THINGS IN THESE NOTES REST ON IT**: the
`git worktree` recipe for separating *this was already failing* from *I broke
it*; the device programme, which is a comparison between two builds throughout;
and `docs/running.rst`'s promise that the answer to *which equilibrium did this
plot come from* is a file that can be re-run. `cmake/GitCommit.cmake` runs as a
**build step** and writes through `copy_if_different`, so the relink happens
when the commit moves rather than once per build, and a `-dirty` suffix comes
from `git status --porcelain --untracked-files=no` — tracked files only, because
this repository's root holds hundreds of untracked run products and a flag that
is always set means nothing.

**AND THE FILE SAID NOTHING ABOUT THE BUILD AT ALL, WHICH IS THE LARGER HALF.**
`meq --version` now reads `MEQ 0.1.0 (6a078e2) (MFEM 4.10.1, double, cuda)` and
the `.nc` carries `meq_precision` and `meq_device_backends` beside
`meq_version`, because a double build and a single build are otherwise two
binaries with one name answering to floors three orders apart
([M-204](MEASUREMENTS.md#m-204)). `meq::realPrecision()` and
`meq::deviceBackends()` are in `src/meq/BuildInfo.cpp` **and not in its header**:
they are answers about MFEM's own preprocessor, and a precision macro tested in
a translation unit that has not included MFEM is not false but UNDEFINED — which
is `MachineFixedBoundary.cpp`'s trap, in `CLAUDE_DEVICE.md`, *Single precision*.

## Code hygiene no compiler checks

**AND A BRACELESS `for` ADOPTS THE NEXT BRACED BLOCK AS ITS BODY, WITH NOTHING
IN `-Wall -Wextra` SAYING SO — NOT EVEN `-Wmisleading-indentation`.** Left in
`GradShafranov.cpp` by a session that was experimenting with the bordered rows:

```cpp
        flush();

        for ( int mode = 0; mode < nModes && mode < 3; ++mode )


        // rowDot() AND THE INVERSE TOGETHER, ...
        {
                LegTimer const denseTimer( ... );
                ...                       // the WHOLE bordered dense solve
        }
```

The loop's body is that block, so the bordered step's row contractions and its
dense solve ran `min( nModes, 3 )` times — and **zero times whenever
`nModes == 0`**, which is every bordered case without exterior DtN coupling.
`rotating-normalised` then stalls with `||r||` CREEPING UP at observed order
1.000, and the q-driven outer loop and the meshed limiter go with it.

**Tested rather than assumed: `g++ -Wall -Wextra`, `g++ -Wmisleading-indentation`
and `clang++ -Wall -Wextra` all compile the shape silently.** A braced compound
statement as a loop body is ordinary C++ and the indentation heuristic does not
fire on it, so **there is no compiler flag that would have caught this**.

**HOW IT SURVIVED IS THE HALF WORTH KEEPING.** The session that left it verified
on three example solves — `machine-f-diiid`, `limited-tokamak`,
`soloviev-nstx` — and five unit binaries, and every one of those is immune:
the first two carry exterior modes so the loop runs at least once, and the third
takes a single Newton step. **Three cases chosen for being quick to compare were
chosen against exactly the configurations the defect reaches.** Same shape as the
six-of-fifty-eight sample under the rename entry above.

**AND IT IS GREPPABLE, WHICH IS THE CHEAP HALF.** The shape is a braceless
`for`/`if`/`while` head whose next *code* line is a bare `{` with blank or
comment lines between them — Allman style puts the brace on the very next line,
so any gap at all is the tell, whatever the indentation. Thirty lines of python
over `src/meq`, `apps` and `tests` reads **zero** today. Run it after any
session that edited a control-flow head; it costs a second and no compiler
offers it.

**The instrument that settled it is the one these notes already prescribe**:
`git worktree add <dir> HEAD --detach`, configure with `-DMFEM_DIR=` and
`-DTOML11_DIR=`, build `meq` alone, and run the failing example against both
binaries. HEAD reached `psi_ax = 1.039325e-01` in 7 Newton steps and the working
tree did not converge in 40 — which is the difference between *this was already
red* and *I broke it*, and no amount of reading the diff establishes it as
cheaply.

**AND NEVER COMMIT A FILE A SUBAGENT OWNS ON THE STRENGTH OF ITS COMPILING.**
an agent was searching for a fixture in
`FreeBoundaryCoupling.cpp`, the file was staged after checking it built, and the
commit captured 99 lines of the agent's temporary `zzExperiment` scaffolding.
**Scaffolding compiles.** The check that was wanted is a `git diff` against what
the file is supposed to contain, or a `grep` for the marker names an agent uses
for its own workings — not a build. Staging early does lock the content against
a later write, which was the right instinct; it locks in whatever is there at
that moment, which was the wrong assumption. **Ask the agent to confirm the file
is clean, and read the diff.**

**AND NEVER `cmake --build` WHILE `ctest` IS RUNNING.** 
for a comment-only header change, which relinked ten test binaries *including
the one ctest was starting*. Nothing errored — replacing a running executable's
file leaves the running process on its old inode — so the run continued and
would have reported a pass over a mix of old and new binaries. A suite result is
only an acceptance measurement if every binary in it is the one being committed;
that run was killed and re-run rather than believed.

**A MECHANICAL IDENTIFIER RENAME IS NOT VERIFIED BY A BUILD AND A LINT, AND
`-Wshadow` IS THE INSTRUMENT THAT SAYS SO.** Renaming the cylindrical radius
from `r` to `radius` across the tree touched about 1,300 sites; the tree built
clean, the `naming` ctest passed, `sphinx-build -W` passed, and **six sites had
changed an answer**. C++ is what makes it silent: an inner declaration may
shadow an outer one of the same type, and a declaration may initialise itself.
So `double const r = radius*sin( theta )` inside a loop over `radius` became
`double const radius = radius*sin( theta )` — legal, and `-nan` at runtime — and
a lambda parameter renamed onto a captured circle radius quietly redirected
`phi = radius*radius - dr*dr - dz*dz` from the circle to the field point, which
is a different level set that still converges.

Configure a second tree with `-DCMAKE_CXX_FLAGS="-Wshadow -Winit-self"` and grep
its build log for the renamed identifier. **Build it CLEAN**: an incremental
rebuild only re-warns about the files just edited, so the second look at a log
is systematically narrower than the first and says so nowhere.

**THE WORST ONE WAS IN THE LOG FROM THE FIRST RUN AND WAS FILTERED OUT, WHICH IS
THE PART WORTH KEEPING.** The warnings split into `shadows a previous local` and
`shadows a member`. The member ones were checked in `GradShafranovSolver`, whose
`radius` is an `mfem::FunctionCoefficient` — no conversion to or from `double`
exists, so a stolen binding cannot compile and the shadowing is provably inert —
and that finding was then applied to the CATEGORY. `meq::analytic::PlasmaEdge`
has a member `radius` too and it is a **`double`**: the disc's radius `a`. Its
fifteen sites were in the lines already filtered away, so
`levelSet()` returned `R² − dr² − dz²` instead of `a² − dr² − dz²` and **the
analytic fixture stopped being the equilibrium it claims to be** — `f()` off
from `−Δ*ψ` by 3.2e+01 against a 1.0e-08 floor. *A reason that makes one
instance safe is not a reason about the category, and the grep that encodes it
is applied to the category.*

**AND THE FIXTURE'S OWN TEST COULD NOT REPORT IT, BECAUSE THAT BINARY WAS
ALREADY THE DELIBERATE RED.** `PlasmaEdgeConvergence` is red by design while the
staircase `Γ_{p,h}` stands, so "PlasmaEdgeConvergence fails" carried no
information. What settled it was **building `HEAD` in a `git worktree` and
diffing the failure inventories**: 4 errors in one case against 318 across seven.
`git worktree add <dir> HEAD --detach`, configured out of it into
`../builds/<name>` with `-DMFEM_DIR=` and `-DTOML11_DIR=` pointed at the main
checkout and `-DMEQ_RUNS_DIR=` at `../runs` — all three defaults are derived
from the worktree's own parent — is the whole recipe, and it is the only check that
distinguishes *this was already failing* from *I broke it*. **A red test cannot
also serve as a regression on something else**, and a tree-wide rename touches
every test including the red ones.

**WHAT NO COMPILER CAN SEE IS `r` THAT NEVER MEANT A RADIUS**, and the tree had
nine: loop counters for a **repeat**, a **refinement level** and a table **row**,
a range-for over `Reading`s, two `r = rate( ... )`, and the `||r||` residual norm
— which the driver had been printing inconsistently as `||R||` in the bordered
step table while the four test harnesses printed `||r||`. Ask `git` for the
originals instead: `git grep -nE "\b(double|int|auto)\s+(const\s+)?r\s*="
HEAD` and read the initialisers.

**THE PYTHON HALF FAILS DIFFERENTLY AND `py_compile` IS FREE.** `!r` in an
f-string is the **repr conversion**, not an identifier, so `{float(v)!r}` became
`{float(v)!R}` and **two tools stopped parsing** — `SyntaxError: f-string:
invalid conversion character`. Nothing in `ctest` builds the tools, so nothing
would have said so. The same pass also rewrote the `k3r0` case-stem scheme that
`race.py` and `MEASUREMENTS.md`'s published tables share, and renamed two result
dicts to `R`. `for f in $(git ls-files '*.py'); do python3 -m py_compile "$f";
done` catches the first class outright.

**AND A FENCED BLOCK HOLDING CAPTURED OUTPUT IS DATA, NOT PROSE.** Three
transcripts were rewritten into records of something no program ever printed —
a Newton table's `||r||` header, a mesh check's `minimum r over the mesh`, an
X-point locator's `r = 0.699700`. Extract every fenced block from `HEAD` and
from the working tree and diff them pairwise; the maths blocks SHOULD move under
a convention change and the transcripts must not, so the two are told apart by
asking which program prints the line.

**A VERIFICATION THAT COMPARES PRINTED VALUES MUST NAME WHICH BINARIES IT RAN.**
This rename was reported verified on "six value-printing tests byte-identical to
baseline", which was true and did not cover any of the cases that had moved. Six
of fifty-eight is a sample, and a sample chosen for being easy to diff is chosen
against exactly the tests whose output is geometry.

## Paths: read beside the file, written beside the run

**A RELATIVE PATH A CONFIGURATION READS IS RELATIVE TO THE TOML FILE; ONE THE
RUN WRITES IS RELATIVE TO THE WORKING DIRECTORY.** Read: the profile tables,
`[initialguess] File` and `MeshFile`, and `[mesh] File` when nothing generates
it. Written: `[output] Directory`, `[conductors] CacheFile`, and the `[mesh]
File` a `[mesh.generate]` block makes. `docs/configuration.rst`, *Relative
paths*, is the list. Two consequences bite:

* **A copy of a configuration moves away from its data.** A test or a tool that
  edits a shipped TOML and writes the variant into its own directory has
  changed what every relative input resolves against; `meq::test::shippedConfig()`
  (`tests/SourceTree.hpp`) re-anchors them, and a hand-rolled copy must do the
  same. `Configuration::fromString()` has no file at all, so its inputs resolve
  against the working directory.
* **Run from a directory under `../runs/`, never from `source/`**: the output,
  a generated mesh and its `.meq-mesh` stamp land in the working directory, and
  untracked output in `git status` is a bug in whatever wrote it.

## Shell and tooling

**Four tooling traps that cost one session six idle polling loops**, and the
reason they are worth writing down is that a polling loop which can never
terminate is *invisible* when the thing it polls for is also reported some other
way.

* **`cd` persists between Bash tool calls**, so a `cd x && ...` that assumes
  `source/` silently runs somewhere else. Use absolute paths.
* **Shell VARIABLES do not persist**, only the working directory does — the same
  trap read the other way round, and the more dangerous half. `L=...` in one call
  is gone by the next, so `until grep -q '^EXIT=' "$LOG"` with an unset `LOG`
  becomes `grep -q '^EXIT=' ""`, and **`grep -q` with no file argument reads
  stdin**. In a detached background process stdin never reaches EOF, so the loop
  blocks silently for the life of the session.
* **A completion marker appended to a log is not necessarily at the start of a
  line.** `echo "EXIT=$?" >> log` after a Boost test binary lands on the end of
  its final ANSI reset sequence, because Boost's coloured output does not
  terminate with a newline — so `grep -q '^EXIT='` never matches. ctest's output
  *does* end with a newline, so the same waiter works on `ctest` and hangs on a
  bare test binary, which gets diagnosed as "the run is slow". Drop the `^`.
* **AND THE WAITER MUST GREP FOR THE MARKER THE WRITER ACTUALLY WROTE**, which
  is the same trap one level up and bites hardest when two parties agree on a
  log file but not on its marker. Measured here: writers emitting
  `echo "CTEST EXIT=$?"` against waiters greping `^EXIT=` — the line begins
  `CTEST `, so the pattern never matches and **fifteen shells polled a finished
  run indefinitely**, at a `sleep` each, for as long as the session lasted. The
  run had succeeded; nothing was wrong but the pattern. **A shared log file is a
  shared protocol**: write the marker and the pattern in one place, or have the
  waiter match the loosest thing that cannot appear early — `EXIT=` unanchored.
* **`grep -c` exits 1 on zero matches**, so a background check reports failure
  spuriously.
* **`pgrep -f` AND `pkill -f` MATCH THE INVOKING SHELL'S OWN COMMAND LINE.**
  This is the one that has now fired three times in a day, twice on the same
  afternoon, and it fails in two different directions. A waiter built on
  `pgrep -f "cmake --install build"` matches *itself*, so it reports the install
  running for as long as the session lasts and the thing it waits for is never
  seen — the install had in fact finished minutes earlier. And `pkill -f ctest`
  in a shell whose own command line contains `ctest` **kills that shell**, which
  surfaces as an unexplained exit code 144 (128 + SIGTERM) rather than as
  anything mentioning the pattern. Use `pgrep -x`/`pkill -x`, which match the
  process *name*, or check for the artefact — a file's timestamp, an exit
  marker — rather than for a process at all.
* **AND `-x` SILENTLY MATCHES NOTHING WHEN THE BINARY'S NAME IS LONGER THAN 15
  CHARACTERS**, which is what the advice above sends you into. `-x` compares
  against the kernel's `comm` field, which is capped at 15 bytes, so
  `pgrep -x PlasmaEdgeConvergence` — 23 characters — never matches and never
  errors: it just answers *not running*. A waiter built on it falls through
  immediately, and what follows reads a **stale binary** and concludes something
  false about it. Measured here: a `--run_test` filter reported "no test cases
  matching" for a case that was in the source, because the waiter had returned
  before the link finished.

  **So both forms fail and in opposite directions** — `-f` never stops, `-x`
  never starts — and `MEQ`'s test binaries are nearly all over 15 characters, so
  the shorter form is the one that looks right here and is not. **Wait on the
  ARTEFACT**: `until grep -q '<marker>' out.txt; do sleep 5; done`, or let the
  harness's own completion notification do it. Five waiters spun for up to four
  hours in one session on `until ! pgrep -f 'build/tests/PlasmaEdge...'`, every
  one of them matching its own command line.

* **AND THE WRITER MAY NEVER WRITE AT ALL, WHICH NO PATTERN CAN FIX.** The
  bullet above is about a waiter and a writer disagreeing on the marker; this
  is the case where the writer is *gone* — killed mid-flight, timed out, or
  restarted under a different filename because the first attempt was being
  redone. The file then sits at its last byte forever and an unbounded waiter
  sleeps on it for the life of the session. Measured here: **three waiters on
  one such file, alive at 19.5 hours**, after the run producing it was killed
  and restarted writing somewhere else. They cost nothing but a `sleep` each
  and were invisible because every answer they were waiting for had already
  arrived by another route.

  **So give a waiter a BOUND**: `for i in $( seq 1 120 ); do <check> && break;
  sleep 30; done` surfaces a dead producer as *gave up* rather than as
  silence. **And when a long run is killed and restarted, kill its waiters
  too** — they are not addressed to the run, they are addressed to a filename.

**If a waiter is used, check it actually fired.**
