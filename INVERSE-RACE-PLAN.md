# The inverse race: four codes, one posed problem, time to accuracy

**What this is.** A time-to-accuracy race of four free-boundary inverse solvers
on the configurations where all four converge to the same answer: MEQ, NICE,
`freegs4e` and freegsnke. The problem is fixed by M-254 and M-255, and those two
entries are the precondition: they establish that the four codes are posed
identically, and how far apart their answers are when each is converged. A race
is only a comparison of costs once it is known that the codes agree on what they
are computing.

`MEASUREMENTS.md` M-155 (MEQ forward against `freegs4e`, 6.2× to 19.8× slower)
and M-157 (MEQ against DESC on the fixed-boundary rung, 2.0× to 9.4× faster) are
the two races already run. This is the first race on the inverse problem and
the first with more than two codes.

---

## 1. The contestants, and what each is

| code | discretisation | inverse | driven by |
|---|---|---|---|
| **MEQ** | HDG, degree `k`, curved half-disc, exterior DtN | the supplies as border unknowns of one Newton; constrained least squares in the corner | `build/meq-run` on a TOML with `[inverse]`; `tools/freegs4e-benchmark/inverse_compare.py`'s arm |
| **NICE** | P1 FEM on its own triangulation, boundary integral at the semicircle | reduced SQP, Tikhonov always present (set to 1e12 A, i.e. off) | `nice/meq-race/inverse/nice_inv_rows`, which is NICE's SQP with the brief's rows |
| **`freegs4e`** | 4th-order finite differences on a rectangle; Green's-function boundary | a vacuum-response least squares every Picard pass | `fgsref.py`'s `SyncConstrain` |
| **freegsnke** | `freegs4e`'s grid, operator and Green's functions | Newton-Krylov with a frozen-plasma then full-Jacobian current update | `tools/freegs4e-benchmark/inverse_nke.py` |

**freegsnke and `freegs4e` are one discretisation and two solvers** (M-255:
they agree to 1e-06). For accuracy they are one vote; for cost they are two
entries, and the difference between them is a solver comparison on a fixed
discretisation, which is a result in its own right.

## 2. The cases: the five shaped machines

**C, D, E, F and G shaped**, against references built by `fgsref.py --shaped
--interpolating-profiles` with `shaped.py`'s turns and quadrature fixes:

* every conductor a uniform rectangle, which all four codes represent exactly
  — NICE natively, `freegs4e` and freegsnke through `rectangle_rule`, MEQ by
  meshing it;
* the analytic profile, which all four solve, one scale closing `I_p`;
* the free supplies those `problem.json` names — as many as the target rows
  determine — and every other supply held at the reference's value;
* targets `fgsref.py`'s own: each design X-point a critical point, each
  isoflux pair levelled.

**Out, each for a stated reason:**
* **A and B** are overdetermined, and the four codes' best fits are four
  different stationary points (M-255: freegsnke's ends worse than it starts).
  A race needs one answer.
* **The plain arms** carry filaments, which neither NICE nor MEQ's inverse can
  represent; their gaps are the conductor model (M-255), not cost.
* **E** races on `[mesh.generate] Symmetric = true` for MEQ, without which
  MEQ's generated mesh breaks the up-down symmetry its five rows rely on.

## 3. The accuracy measure

**Each code against its own finest run**, as `race.py` and
`time_to_accuracy.py` already do, because the cross-code floor is not any one
code's error:

* **the free currents**, max |ΔI| over the free supplies over the largest free
  current — the quantity the inverse solve exists to produce;
* `psi_ax` relative, the located X-point and axis, for the record.

**The cross-code floor is reported separately and is the ceiling on any claim
about accuracy between codes**: 3e-04 to 4e-03 on these five (M-255). A
"faster at 1e-04" claim compares each code with itself.

Accuracy targets: **1e-02, 1e-03, 1e-04, 1e-05** in the currents.

## 4. The ladders

| code | knob | rungs |
|---|---|---|
| MEQ | `PolynomialDegree` × `PlasmaSize`/`CoilSize` | `k` = 1–4 on meshes from coarse to the shipped, as `time_to_accuracy.py` builds them — the coarse-mesh-high-order corner is where MEQ's claim lives |
| NICE | `--area-scale` over the whole `*TriangleArea` family | 1, 1/4, 1/16, 1/64 — 1/64 is 7–15 GB, and D shaped has been killed for memory three times |
| `freegs4e` | `--nx` | 65, 129, 257, 513 — the cached boundary matrix is 4.3 GB at 513 |
| freegsnke | grid | 65, 129, 257 — 513 was killed at 16 GB on F |

**The finest rung of each is the accuracy datum and is NOT on the clock as a
contestant.** A ladder that cannot afford one rung beyond its last timed one
cannot report its last timed one's error.

## 5. The two starts, and they answer different questions

**COLD — the machine and the targets, and nothing solved for.** A common cold
state every code can be handed: the design plasma footprint carrying `I_p`
(`mkcoldguess.py`'s blob), with the free supplies at the **vacuum stage** —
`freegs4e`'s first control step, the minimum-norm currents meeting the targets
against that blob (`inverse.py`'s `vacuum_currents`). Every code gets those
currents as its starting supplies and its own ordinary representation of the
blob. **A configuration a code does not solve from cold is a DNF, with the
reason**, not a gap: M-253 and M-255 already say MEQ will not finish C shaped
cold because its support is frozen per sweep, and that is part of the result.

**WARM — a shape step, the scan workload.** Each code starts from ITS OWN
converged answer for the design targets; the targets then move — the X-points
1 cm outward, each isoflux point 1 cm — and the re-solve is timed. This is what
a shape scan or a controller does, it needs no shared guess, and it is the start
where a Newton method's advantage, if it has one, shows.

## 6. The protocol, inherited, every clause of which was learned once

From `race.py`, `race_nke.py` and M-102, M-111, M-155, M-157:

* **Interleaved per rung**, never arm by arm, so a machine that drifts drifts
  through all four.
* **Quiet first, checked before every rung**: the load average sampled just
  before each rung, against a threshold that allows one busy arm. No other
  agent or build on the machine for the duration.
* **One case per invocation**, the machine settling between them — 9.49 s
  inside a three-case run against 5.41 s settled, M-111.
* **Median of three repeats**; best-of answers a different question.
* **Threads matched** at the physical core count (`cores.py`, 8 here):
  `OMP_NUM_THREADS = MKL_NUM_THREADS = 8` for MEQ and NICE, the BLAS's threads
  for the Python arms.
* **End to end on the clock**, and the solve alone reported beside it where the
  code can say: MEQ's includes gmsh and its output writers, NICE's its mesher,
  `freegs4e`'s its boundary matrix.
* **Provenance on every row**: every tree's SHA, `libmfem.a`'s mtime, NICE's
  build, the freegs4e and freegsnke versions (freegsnke 3.1.0 pins
  `freegs4e >= 0.14`; the clone is 0.12).
* **A run is counted only if its rows are met** — to 1e-08 T and Wb/rad,
  checked by the harness, not by the code's own exit status. freegsnke reports
  SUCCESS with rows at 0.1 T (M-255), and its stopping rule is a step
  ratchet, not a residual: its arm sets its tolerances so the rows decide.

## 7. The infrastructure — BUILT

| piece | file | what it does |
|---|---|---|
| the contract | `tools/freegs4e-benchmark/RACE-ARMS.md` | the command every arm takes and the record it writes; rows evaluated from the converged field, never the solver's stopping test |
| the problem | `tools/freegs4e-benchmark/inverse_problem.py` → `ref-interp/problem.json` | the determined free subset, the vacuum-stage cold seed, the 1 cm shifted warm targets, the symmetry flag, from `fgsref.py --shaped --interpolating-profiles` references in `ref-interp/` |
| MEQ | `race_arm_meq.py` | the shipped example with `[inverse]`, rung = degree × every mesh size; warm restarts on the cold run's own mesh and field, so nothing is re-meshed on the clock; `Symmetric` meshing where the rank is below the row count |
| NICE | `nice/meq-race/inverse/race_arm.py` | area-scale rungs; warm reuses the converged nodal ψ and currents through NICE's restart route but re-meshes and re-assembles, NICE having no saved-operator restart |
| `freegs4e` | `race_arm_fgs.py` | `fgsref.py`'s solve without its blend and reseed retries; grid rungs 65² to 513² |
| freegsnke | `race_arm_nke.py` | grid rungs 65² to 257²; rows, not the step ratchet, decide `finished` |
| the orchestrator | `race_inverse.py` | `run` and `summarise`, the protocol of §6 |

**Two things the building found.**

* **MEQ's ladder has no `k = 1` rung.** At `k = 1` the X-point border's Newton
  rests on a `∇q` that is constant per element, and on DIII-D shaped the first
  bordered Newton throws the X-point out of the mesh even restarted from its
  own converged answer with the same targets — 41 iterations on the retry rung
  against 2 at `k = 2`. A `k = 1` rung would time that defect.
* **MEQ's cold start is its expensive one, and the race will say so.** On F
  shaped at `k = 2` the vacuum-seeded cold start runs 153 Newton iterations on
  the first rung without converging and the retry finishes it: 352 s, against
  11.2 s for the warm shape step (2 Newton iterations) and `freegs4e`'s 13.1 s
  cold at 129². M-253's frozen support, measured as a cost.

**What NICE cannot honour**: a symmetric mesh (so C, E and G are best fits to
the discretisation's level there), a restart without re-meshing, and threads
beyond MKL's (it is built without OpenMP). Its walls, which it requires, are
the v3 reference outline dilated clear of the coils for C, D, E and G and a
hand-drawn polygon for F (`race_walls.json`): geometry entering no row, and the
one input not drawn from the problem file.

## 8. Stages

| | | acceptance | |
|---|---|---|---|
| **IR-0** | the four arms on F shaped at one rung each, untimed | every arm's record has rows met and currents within M-255's floor of the others | **done** |
| **IR-1** | F shaped, every rung, both starts, timed | a Pareto front per code; the finest rungs self-converged | **done**, `race-runs/ir1-f-final`; MEQ re-run as `ir3-f`, under `../runs/freegs4e-benchmark/race-runs/` and untracked |
| **IR-2** | C, D, E, G shaped, one invocation each | the same, and the DNFs with their reasons | **C and D done**, `race-runs/ir2-c`, `ir2-d-final`, MEQ re-run as `ir3-c`, `ir3-d` beside `ir3-f`; on C MEQ converges at two rungs of six (`PERFORMANCE_TODO.md`) and NICE reaches only a best fit. **E and G not run** |
| **IR-3** | the publication | a new `M-nn` per start: per case and target the cheapest run of each code, the cross-code floor beside it | **open**: the results are the Claude Docs artifact *Inverse race: F, C and D shaped* (`claude.ai/artifact/WAq3D1udMQakGMYyvHdL1g`) and no `M-nn` is written |

## 9. What the result can and cannot say

**It can say** what each code costs to reach a given self-consistent accuracy on
an identically posed inverse problem, and which codes can start from cold.

**It cannot say** which code is more ACCURATE below the cross-code floor,
3e-04 to 4e-03, and the publication must not let a reader infer it. What sets
that floor is open: about 3e-04 between MEQ and the finite-difference family on
F shaped with both converged, NICE between them within its own resolution.

**Expect honestly**: M-155 puts MEQ 6× to 20× slower than `freegs4e` on the
forward problem at the accuracies `freegs4e` reaches, and nothing about an
inverse border changes the per-unknown cost that set that. M-157 says MEQ's
margin is at high accuracy and with few elements. If MEQ wins anything here it
is the 1e-04 and 1e-05 targets, and the warm start, where a quadratic Newton
pays and a Picard iteration does not. The race is worth the machine time
because nothing else measures that, not because MEQ is expected to win it.
