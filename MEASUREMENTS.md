# MEQ's published tables

**Every number MEQ claims, in one place.** These are MEQ's measurement tables,
kept here rather than inline because a maintainer does not need them on every
reading and does need to know where they are. **The argument each one supports
stays beside a pointer to the row**, in whichever of the four maintainer files
owns the subject: `CLAUDE.md` (the index, the build, the traps),
`CLAUDE_HDGGS.md` (the equation, the discretisation, the solve, the linear
solves), `CLAUDE_FB.md` (free boundary and the `freegs4e` benchmark),
`CLAUDE_INVERSION.md` (the flux-surface work) or `CLAUDE_FLOW.md` (toroidal
flow). **Do not renumber the anchors** — all five point at them by `M-nn`.

**THE RULES THEY WERE TAKEN UNDER HAVE NOT CHANGED, AND THEY ARE WHAT MAKES A
TABLE HERE WORTH ANYTHING.**

* **Every stage ends at a measured number, not at "it runs".** That is what
  caught the Solov'ev coefficients, the `τ` sign, the `ψ*` regression, the
  hybridization ordering, the under-resolution and the inert normalisation.
* **A property measured on the easy configuration is not a property of the
  code.** Symmetry held to 2e-16 on a fitted rectangle and failed at 5.4e-1 on
  the geometry MEQ exists for.
* **A difficulty measured at one resolution is not a property of the problem.**
* **A measurement taken on this machine may not be a measurement about the
  code** — threaded MKL decides a marginal Newton, a consumer GPU decides
  nothing about a production one, and a suite time is only a measurement on an
  idle machine.
* **A column that does not move when something that should move it changes is
  measuring the instrument.** This tree has mistaken an instrument for a result
  at least five times; the tell is always the same.
* **A control is not optional.** Most tables here carry one, and where a table
  has an "off" or "removed" column that column is the whole point: without it
  the rest is compatible with the code doing nothing.

**RE-MEASURE BEFORE QUOTING.** A number here was true when it was taken. Several
have already moved once — under the NPC port, under the `ψ_bnd` repair, and under
`AxisConstraint::LocatedAxis` — and this file records what a stale number costs
rather than pretending it cannot happen.

Anchors are stable: `CLAUDE.md` points at `M-nn`, so do not renumber.

---


## Status: the solver works, and MEQ is a program

### M-01

| | per cycle | total | final L2 |
|---|---|---|---|
| cold | 4, 4, 4 | 12 | 3.652757e-06 |
| **warm** | **4, 2, 2** | **8** | 3.652757e-06 |


## Free boundary: MEQ SOLVES A MACHINE CASE, AND IT AGREES WITH freegs4e

### M-02

| | freegs4e, 513² | MEQ | apart |
|---|---|---|---|
| `ψ_ax` | 9.308752e-02 | **9.307342e-02** | **1.5e-04** |
| `ψ_bnd` | 2.622462e-02 | **2.622580e-02** | **4.5e-05** |
| profile scale | 0.939672, predicted from `L` | **0.9400254** | 3.8e-04 |
| `I_p` | 3.0e+05 A | 3.000003e+05 A | it is the constraint |


## One argument where two were meant: psi_bnd was zeroed inside the Jacobian window

### M-03

| | |
|---|---|
| the Newton direction against its own linearised system | **109% wrong**; 1.167e-07 once repaired |
| the current column `∂R/∂λ` | **46.5%** wrong; 1.2e-10 repaired |
| the corner entry `D(λ, ψ_bnd)` | **92.8%** wrong — a factor of **13.9**; 1.4e-10 repaired |
| `∫F/r` fed to the border | 4.4217e-01 against a true 3.3768e-01, **and the OPPOSITE SIGN** against the target: `+6.55e-02` where the truth is `−3.93e-02` |


## FB-5: the exterior coupling is an unknown of the same Newton

### M-04

| `n` | `\|a − exact\|` | `\|a − superposition\|` | Newton |
|---|---|---|---|
| 12 | 1.32e-03 | **2.34e-15** | **1** |
| 24 | 1.32e-04 | **1.67e-15** | **1** |

### M-05

| cycle | `Γ_h` faces, off | `\|a − exact\|`, off | `Γ_h` faces, **on** | `\|a − exact\|`, **on** | `η₆` |
|---|---|---|---|---|---|
| 0 | 34 | 1.3194e-03 | 34 | 1.3194e-03 | 8.5849e-04 |
| 1 | 34 | 1.3194e-03 | **46** | 1.1259e-03 | 5.3416e-04 |
| 2 | 34 | 1.3194e-03 | **57** | 2.8643e-04 | 2.0841e-04 |
| 3 | 34 | 1.3194e-03 | **64** | **1.5297e-04** | 1.2697e-04 |

### M-06

| `Γ_h` faces | 71 | 142 | 282 | 570 |
|---|---|---|---|---|
| `η₅` on `Γ_h` | 2.864e-01 | 4.262e-01 | 6.151e-01 | 8.787e-01 |
| ratio, against `√2 = 1.4142` | — | 1.488 | 1.443 | **1.429** |


## FB-7: a conductor outside `Γ`, and the type that carries it has no `f()`

### M-07

| | |
|---|---|
| the datum given, `ψ` against `ψ_coil` | **1.980 / 3.409 / 4.069** at `k = 1,2,3`, against a datum-removed control **148,165×** larger |
| the coupled path | **one Newton step** at `n = 12/24/48` — the residual is affine in `( x, a )` and an exact Jacobian must finish it in one — with `\|a\|` at **4.076** and `ψ` at **4.486** |
| a rectangle against a filament, same current and centre | **4.242e-04** on `Γ` and **4.076e-04** in the domain, about **3×** §7.16's published 1.3e-04 agreement in `ψ_ax` |
| the same, through the **coupled** path | **0.9998** of the datum-given answer |
| **the two conductor routes against each other** | **1.951 / 3.727 / 3.879** at `k = 1,2,3`, against a control **6306×** larger |


## FB-6 is beaten, and what was in the way was the limiter, pinned to a dof

### M-08

| | | | | order | limit | vs the 129² reference |
|---|---|---|---|---|---|---|
| **NearestDof** | 9.466087e-02 | 9.464282e-02 | 9.464039e-02 | 2.89 | 9.46400e-02 | 2.0e-03 |
| **ExactPoint** | 9.484400e-02 | 9.482652e-02 | 9.482423e-02 | 2.93 | **9.482388e-02** | **7.9e-05** |


## The plasma edge caps the order, and it is the PROFILE that sets the cap

### M-09

| | best approximation, ANY method | MEQ, plain Gauss rule |
|---|---|---|
| `ψ_h` | `min( k+1, j+2.5 )` | `min( k+1, **j+1.5** )` |
| `q_h` | `min( k+1, j+1.5 )` | `min( k+1, j+1.5 )` — **already at its own bound** |
| `ψ*` | `min( k+2, j+2.5 )` | `min( k+2, j+2 )` at the shipped rule, **`j+2.5`** once it is raised |

### M-10

| `j` | `k = 1` | `k = 2` | `k = 3` | `k = 4` |
|---|---|---|---|---|
| 0 | 1.89 | 1.95 | 1.70 | — |
| 1 | **3.00** | 2.87 | 2.68 | — |
| 2 | **3.00** | **3.99** | 3.88 | — |
| 3 | **3.00** | **4.00** | **4.99** | ~5, lost |


## At `j = 0` the question does not arise: Newton cannot solve it at all

### M-11

| | 401 samples | 801 | 1601 |
|---|---|---|---|
| `j = 0`, plain rule | 2.871e-04 | 2.909e-04 | **2.928e-04** |
| `j = 0`, composite | 1.798e-05 | 1.225e-05 | **8.794e-06** |
| `j = 1`, plain | — | 2.18e-07 | — |
| `j = 2`, plain | — | 8.26e-09 | — |

### M-12

| q8 | q12 | q20 | q40 | q80 |
|---|---|---|---|---|
| 2.61e-04 | 1.01e-04 | 1.34e-05 | 6.68e-08 | **5.40e-10** |

### M-13

| `k` | `ψ`, axis / control | `q`, axis / control | conditioning ratio grows by |
|---|---|---|---|
| 1 | **2.000 / 2.000** | 1.79 / 2.00 | 7.169 |
| 2 | **3.000 / 3.000** | 2.51 / 3.00 | 7.536 |
| 3 | **4.000 / 4.000** | 3.999 / 3.00 | 7.631 |


## Commands

### M-14

| | wall | CPU |
|---|---|---|
| `-j1`, `OMP=16` — the old default | **710.5 s** | — |
| `-j16`, `OMP=1` | 265.2 s | 350% |
| **`-j4`, `OMP=4`** | **223.2 s** | 514% |


## A third erratum: eq. (20)'s `η₅` does not vanish on the exact solution

### M-15

| | coarsest | finest | rate |
|---|---|---|---|
| `η₁`, for scale | 2.1161e-03 | 4.7256e-05 | — |
| **`η₅` on the datum** | **9.5889e-05** | **2.0281e-06** | **2.78** |
| `η₅` on the pinned zero | 4.0688e-01 | 2.3226e-01 | 0.40 |
| `η` with the datum | 2.1746e-03 | 4.8375e-05 | 2.75 |
| `η` with the faces excluded | 2.1745e-03 | 4.8375e-05 | 2.75 |


## What NPC bought, measured

### M-16

| `k` | NPC local NL its | condense local NL its | L2 vs exact | relative in `ψ` |
|---|---|---|---|---|
| 1 | **0** | 3644 | 7.614902e-03 both | 2.9e-15 |
| 2 | **0** | 3560 | 2.494695e-04 both | 1.6e-11 |
| 3 | **0** | 3412 | 6.000146e-06 both | 9.0e-12 |

### M-17

| | NPC | CondenseThenLinearise | |
|---|---|---|---|
| §4.2 pedestal `k = 2, n = 24` | **9 its, 0.92 s**, 0 local NL | 10 its, **3.93 s**, 75,490 local NL | agree to 4.1e-11 |
| similarity `k = 2, n = 16` | **4 its, 0.19 s**, 0 local NL | 4 its, **0.70 s**, 10,898 local NL | agree to 3.6e-11 |

### M-18

| | condensation | **NPC** |
|---|---|---|
| `b = −∂(max ψ_h)/∂x` | `3(k+1)` central differences over one element's trace dofs | **exactly `−e_j`**, one entry, not differenced |
| `d = ∂G/∂s` | `1 −` a central difference | **exactly `1`** |
| `c = ∂R/∂s` | one central difference in a scalar | the same |

### M-19

| `ν` | `A` | `ψ_ax` | `ψ_ax − max ψ_h`, was | now | Newton, was | now |
|---|---|---|---|---|---|---|
| 2 | 1 | 3.058984e-01 | 1.7e-16 | 0.0 | 4 | 4 |
| 2 | 10 | 9.607537e-01 | −4.4e-16 | 0.0 | 4 | 4 |
| 2 | 100 | 3.036075e+00 | 4.4e-16 | 0.0 | 4 | 4 |
| 4 | 1 | 2.834510e-01 | **3.1e-12** | −5.6e-17 | 8 | **6** |
| 4 | 10 | 8.643745e-01 | 8.9e-16 | 0.0 | 11 | **7** |
| 4 | 100 | 2.720222e+00 | **1.0e-10** | 0.0 | 10 | **7** |


## What NPC cost, and it is the parity gap in the other direction

### M-20

| case | **NPC** | CondenseThenLinearise |
|---|---|---|
| §4.2 pedestal `k = 1, n = 16` | **fails at 60** | ok, 14 |
| §4.2 pedestal `k = 1, n = 24` | ok, 26 | ok, 23 |
| §4.2 pedestal `k = 1, n = 32` | ok, 10 | ok, 9 |
| §4.2 pedestal `k = 2, n = 16` | ok, 17 | ok, 11 |
| §4.2 pedestal `k = 2, n = 24` | ok, **9** | ok, 10 |
| §4.2 pedestal `k = 3, n = 16` | ok, **8** | ok, 12 |
| §4.5 layer `k = 1, n = 30` | **fails at 60** | ok, 13 |
| §4.5 layer `k = 1, n = 60` | ok, 10 | ok, 10 |
| §4.5 layer `k = 2, n = 16` | ok, **51** | ok, 10 |
| §4.5 layer `k = 2, n = 24` | ok, 13 | ok, 12 |
| §4.3 barrier `k = 1, n = 16` | fails at 60 | **fails at 60** |
| §4.3 barrier `k = 2, n = 16` | fails at 60 | **fails at 60** |
| similarity solution `k = 2, n = 16` | ok, 4 | ok, 4 |
| similarity solution `k = 3, n = 16` | ok, 4 | ok, 4 |

### M-21

| | wall clock | element-local non-linear iterations |
|---|---|---|
| **NPC** | **1.6 s** | **0** |
| CondenseThenLinearise | **53.5 s** | **2,027,130** |


## Why it fails, measured — and why a line search does not fix it

### M-22

| | flux | potential | trace |
|---|---|---|---|
| residual at `x` | 7.93e-02 | 8.13e-02 | **0.00e+00** |
| residual after a full step | **1.14e-16** | **6.25e+00** | 1.82e-13 |
| the correction `c` | **3.50e+02** | **1.51e+02** | 1.48e+01 |

### M-23

| | NPC, plain | NPC, `PicardThenNewton` |
|---|---|---|
| §4.2 pedestal `k = 1, n = 16` | fails at 60 | **ok, 5** |
| §4.5 layer `k = 1, n = 30` | fails at 60 | **ok, 5** |


## There is no cheap discriminator, and the obvious one is ANTI-correlated

### M-24

| case | out | its | `‖r₀‖` | `‖r₁‖/‖r₀‖` |
|---|---|---|---|---|
| §4.2 pedestal `k = 1, n = 16` | **FAIL** | 60 | 1.136e-01 | 5.50e+01 |
| §4.2 pedestal `k = 1, n = 24` | ok | 26 | 9.242e-02 | **1.17e+03** |
| §4.2 pedestal `k = 1, n = 32` | ok | 10 | 7.988e-02 | 2.28e+02 |
| §4.2 pedestal `k = 2, n = 16` | ok | 17 | 1.124e-01 | 2.91e+02 |
| §4.2 pedestal `k = 2, n = 24` | ok | 9 | 9.172e-02 | 5.04e+01 |
| §4.2 pedestal `k = 3, n = 16` | ok | 8 | 9.542e-02 | 5.72e+01 |
| §4.5 layer `k = 1, n = 30` | **FAIL** | 60 | 8.253e-02 | 3.79e+02 |
| §4.5 layer `k = 1, n = 60` | ok | 10 | 5.817e-02 | 3.55e+01 |
| §4.5 layer `k = 2, n = 16` | ok | 51 | 1.124e-01 | 2.91e+02 |
| §4.5 layer `k = 2, n = 24` | ok | 13 | 9.172e-02 | 5.04e+01 |
| §4.3 barrier `k = 1, n = 16` | **FAIL** | 60 | 1.136e-01 | **8.56e-01** |
| §4.3 barrier `k = 2, n = 16` | **FAIL** | 60 | 1.124e-01 | **8.67e-01** |
| similarity `k = 2, n = 16` | ok | 4 | 1.152e-01 | 5.18e-02 |
| similarity `k = 3, n = 16` | ok | 4 | 9.784e-02 | 5.23e-02 |


## Should `PicardThenNewton` simply be the default? No, and cost is the weakest of the three reasons

### M-25

| case | plain | | `PicardThenNewton` | | slowdown | rel. `max ψ_h` |
|---|---|---|---|---|---|---|
| | out, its | s | out, picard+newton | s | | |
| §4.2 pedestal `k = 1, n = 16` | FAIL 60 | 1.40 | **ok** 60+5 | 0.67 | **0.48** | — |
| §4.2 pedestal `k = 1, n = 24` | ok 26 | 1.37 | ok 60+11 | 1.85 | 1.35 | **6.1e-06** |
| §4.2 pedestal `k = 1, n = 32` | ok 10 | 0.98 | ok 60+3 | 2.60 | 2.65 | 3.7e-16 |
| §4.2 pedestal `k = 2, n = 16` | ok 17 | 0.61 | ok 60+2 | 0.88 | 1.43 | 0.0 |
| §4.2 pedestal `k = 2, n = 24` | ok 9 | 0.78 | ok 60+1 | 1.98 | 2.52 | 5.6e-16 |
| §4.2 pedestal `k = 3, n = 16` | ok 8 | 0.52 | ok 60+2 | 1.61 | 3.10 | 5.6e-16 |
| §4.5 layer `k = 1, n = 30` | FAIL 60 | 5.02 | **ok** 60+5 | 2.61 | **0.52** | — |
| §4.5 layer `k = 1, n = 60` | ok 10 | 3.72 | ok 60+3 | 11.41 | 3.07 | 8.9e-16 |
| §4.5 layer `k = 2, n = 16` | ok 51 | 1.81 | ok 60+4 | 0.93 | 0.52 | **1.0e-02** |
| §4.5 layer `k = 2, n = 24` | ok 13 | 1.09 | ok 60+3 | 2.15 | 1.97 | **1.6e-05** |
| §4.3 barrier `k = 1, n = 16` | FAIL 60 | 1.52 | **FAIL** 0+60 | 2.22 | 1.46 | — |
| §4.3 barrier `k = 2, n = 16` | FAIL 60 | 2.28 | **FAIL** 0+60 | 3.14 | 1.38 | — |
| similarity `k = 2, n = 16` | ok 4 | 0.16 | ok 13+1 | 0.30 | 1.87 | 1.8e-16 |
| similarity `k = 3, n = 16` | ok 4 | 0.28 | ok 13+1 | 0.54 | 1.93 | 1.8e-16 |

### M-26

| route | Newton steps | `max ψ_h` |
|---|---|---|
| NPC, plain | 51 | 3.1831e-01 |
| CondenseThenLinearise, plain | 10 | 3.4779e-01 |
| NPC, `PicardThenNewton` | 4 | 3.1514e-01 |

### M-27

| case | NPC plain | `CondenseThenLinearise` | NPC `PicardThenNewton` |
|---|---|---|---|
| §4.2 pedestal `k = 1, n = 16` | FAIL 60, 1.48 s | ok 14, 1.91 s | **ok 5, 0.74 s** |
| §4.5 layer `k = 1, n = 30` | FAIL 60, 5.31 s | ok 13, 5.75 s | **ok 5, 2.84 s** |
| §4.5 layer `k = 2, n = 16` | ok 51, 1.87 s | ok 10, 1.99 s | **ok 4, 1.04 s** |
| §4.3 barrier `k = 1, n = 16` | FAIL 60, 1.60 s | FAIL 60, **53.5 s** | FAIL 60, 2.27 s |
| §4.3 barrier `k = 2, n = 16` | FAIL 60, 2.36 s | FAIL 60, **96.6 s** | FAIL 60, 3.77 s |


## KINSOL: the two orderings are exact opposites

### M-28

| | NPC | CondenseThenLinearise |
|---|---|---|
| plain Newton | ok, 8–10 | ok, 7–10 |
| `KIN_NONE` | **ok, 8–10** | **fails at 60, every case** |
| `KIN_LINESEARCH` | **fails at 60, every case** | ok, 24–31 |


## ~~A warm start no longer shows up in `‖r₀‖`, and that is structural~~ — THE FLUX IS SEEDED NOW, AND IT BUYS A DIAGNOSTIC AND NOT ONE ITERATION

### M-29

| | cold | warm |
|---|---|---|
| before | 1.771e-01 | 2.638e-01 — *worse than cold* |
| **after** | 1.770849e-01 | **2.130472e-03** |


## The measurement, and the control that makes it mean something

### M-30

| `ν` | `A` | `ψ_ax` | ratio to estimate | `max\|∂F/∂ψ\|/λ₁` | Newton |
|---|---|---|---|---|---|
| 2 | 1 | 3.059006e-01 | 1.0209 | 1.88 | 4 |
| 2 | 10 | 9.607543e-01 | 1.0139 | 1.91 | 4 |
| 2 | 100 | 3.036075e+00 | 1.0132 | 1.91 | 4 |
| 4 | 1 | **3.068708e-01** | 0.7242 | 13.35 | **10** |
| 4 | 10 | **9.462726e-01** | 0.7061 | 14.17 | **9** |
| 4 | 100 | **2.981136e+00** | 0.7035 | 14.27 | **8** |

### M-31

| | residual | iterations |
|---|---|---|
| coupled | 8.310e-02 → **4.447e-15** | 4 |
| decoupled | 8.310e-02 → 8.239e-02 | 15, not converged |


## `ψ_bnd`: settable, an unknown, and reachable from a file

### M-32

| | |
|---|---|
| table over `[0, 1]` | **does not converge** — creeps at 0.99 a step |
| table over `[−0.6, 1.4]` | **7 Newton steps**, `ψ_ax = 9.76e-02` |
| the same plus `ConfineToPlasma = true` | 92 steps, `ψ_ax = **4.08e-01**` |

### M-33

| limiter `R` | Newton | final residual | `ψ_ax` | `ψ_bnd` |
|---|---|---|---|---|
| 1.05 | 5 | 3.2276e-15 | 1.088199362e-01 | 8.557804281e-03 |
| 1.15 | 9 | 4.2489e-15 | 1.117391871e-01 | 1.405147764e-02 |
| 1.20 | 6 | 1.0735e-16 | 1.091632931e-01 | 9.210810991e-03 |
| 1.30 | 5 | 2.9677e-15 | 1.088257139e-01 | 8.568770649e-03 |


## Why MEQ's Newton struggles where other Newton solvers do not

### M-34

| | `k = 1` | `k = 2` | `k = 3` |
|---|---|---|---|
| §4.2 pedestal | **42**, **22**, 9, 9 | 10, 9, 8, 6 | 11, 9, 7, **7** |
| §4.3 barrier | fail, fail, **7**, 8 | 10, 7, 7, 8 | —, 8, 8, 7 |
| §4.5 layer | fail, **17**, 11, 11 | —, 12, 14, 11 | —, 10, 21, 12 |
| §4.4 current hole | **abort, abort, fail** | abort, fail, fail | —, **fail, fail, fail** |

### M-35

| `c₃` | 0 | −2 | −4 | −9 | **−18** |
|---|---|---|---|---|---|
| `max\|∂F/∂ψ\| / λ₁` | 7 | 7 | 8 | 13 | **26** |


## Picard, then Newton — the route for a coarse mesh

### M-36

| case | Newton alone | Picard | **Newton from Picard** |
|---|---|---|---|
| §4.2 pedestal `k = 1, n = 16` | ok, 31 | ok, 197 | **ok, 4** |
| §4.2 pedestal `k = 1, n = 24` | ok, 7 | ok, 290 | **ok, 4** — order 2.02 |
| §4.3 barrier `k = 1, n = 16` | **fails at 60** | ok, 122 | **ok, 4** — order 2.01 |
| §4.3 barrier `k = 2, n = 16` | ok, 10 | ok, 120 | **ok, 3** — order 1.96 |
| §4.5 layer `k = 1, n = 16` | **fails at 60** | 400, not converged | **ok, 28** |
| §4.5 layer `k = 2, n = 16` | **fails at 60** | 400, not converged | **ok, 5** |
| §4.4 hole `k = 1, n = 16` | aborts | not converged | **aborts** |


## The suite is green, and the last red one was a mesh chosen for a dead solver

### M-37

| §4.5, `k = 1` | `ψ` | `q` |
|---|---|---|
| `{ 24, 48, 96 }` | 1.786 | 2.480 |
| **`{ 48, 96, 192 }`** | **2.160** | **2.184** |

### M-38

| `n` | `h`/0.025 | plain NPC |
|---|---|---|
| 24 | 1.33 | fails at 60, **wandering** to `ψ ∈ [−1.51, +1.44]` |
| 32 | 1.00 | fails at 60, in a **PERIOD-4 LIMIT CYCLE** at ~4e-3 |
| 48 | 0.67 | ok, 12 |

### M-39

| | trace dofs | plain NPC |
|---|---|---|
| uniform `n = 32` | 6,272 | **fails at 60** |
| **band-refined from `n = 24`** | **7,588** | **ok, 22** |
| uniform `n = 48` | 14,016 | ok, 12 |

### M-40

| | | |
|---|---|---|
| `k = 1, n = 16` | 42 iterations | the knife edge, recorded and *not* asserted on |
| `k = 1, n = 32` | 9 iterations | `h`-refinement cures it |
| `k = 2, n = 16` | 11 iterations | `p`-refinement cures it, mesh untouched |


## What is measured

### M-41

| | |
|---|---|
| `φ₀` against a **brentq root of (97)**, independent Python | 1.9e-14 |
| `Σ_s Z_s n_s = 0` at every radius, over a Mach sweep | 1e-13 |
| `p` against Li & Zhu (8) | 1e-14 |
| `dFdPsi` against a central difference, five Mach scales, two steps | the `O(h²)` floor |
| `ω = 0` against `MHDSource`, pointwise | 1e-13 |
| **`ω = 0` through the solver** | reproduces `SolovievConvergence`'s errors **to every printed digit** |
| rotating Solov'ev, `k = 1,2,3` | 1.995 / 2.995 / 3.995 in `ψ`, 1.980 / 2.985 / 3.984 in `q`, Newton = 1 |
| source against fixture, `M² = 0, 1, 4` | asserted at 1e-12 — two independent implementations of the same physics |
| assembled Jacobian vs a difference of the assembled residual | 3.1e-11 |
| Newton order on the manufactured nonlinear case | 1.980, and `k+1` at 2.007 / 2.999 / 4.002 |
| root find vs closed form at two species | `φ₀` 1e-12, `∂φ₀/∂ψ` 1e-11, `F` 1e-11, `∂F/∂ψ` 1e-9 |
| three species (D, C⁶⁺, e), `Σ_s Z_s n_s = 0` | 1e-12 at every radius |
| bordered Newton, the constraint residual | **0.000e+00** on three meshes; tail order 2.000. Measured as `ψ_ax − max ψ_h` before 2026-09-07 and as `normalisationResidual()` since, the definition having moved |


## IN-A: the axis as a root of `q`, and a degree is never a count

### M-42

| | `k = 1` | `k = 2` | `k = 3` |
|---|---|---|---|
| position of the axis, rate over the whole sequence | **2.340** | **3.484** | **4.447** |


## `CriticalPointFinder`'s axis and `GradShafranovSolver::psiAxis()` are THE SAME POINT since 2026-09-07

### M-43

| limiter | `gg′` | `ψ_bnd` | `Ψ_axis` | `F( 0, z )` | verdict |
|---|---|---|---|---|---|
| no | 0.05 | 0 | 0 | 0.0e+00 | AGREES |
| **yes** | **0.05** | 2.74e-02 | **−2.80e-01** | **1.43e-01** | **REFUSES** |
| no | 0 | 0 | 0 | 0.0e+00 | AGREES |
| yes | 0 | 3.67e-02 | −3.96e-01 | 0.0e+00 | AGREES |


## §7.12b's fixture was the defect, and the suite is 41/41

### M-44

| removed | what happens |
|---|---|
| `ConfineToPlasma` | `\|F\|` on `r = 0` at **1.0e-02 of scale**; with it, **exactly 0.0** |
| the prescribed current | **does not converge at any of the four radii** — §7.14's non-linear eigenvalue problem, `Λ = A/span²` on a region that is itself unknown, so scaling `A` changes nothing |
| the vertical field | converges at limiter 1.05 to §7.14's **annulus**, axis at `r = 1.38` on a domain reaching 1.50; **does not converge at all** at 1.15 or 1.20 |

### M-45

| limiter | coil `μ₀I` | coils | Newton | `ψ_ax` | `ψ_bnd` | `\|F\|` on `r = 0` | `ψ_ax` attained at | `Ψ` at the O-point |
|---|---|---|---|---|---|---|---|---|
| 1.08 | −7.8072e-02 | yes | 11 | 1.212249e-02 | 1.366962e-03 | **0.0000e+00** | ( 0.779, 0.000 ) | **1.0000** |
| 1.08 | 0 | **NO** | 96 | 6.222586e-02 | 3.854703e-02 | 0.0000e+00 | **( 1.381, 0.000 )** | 1.0000 |
| 1.10 | −7.6159e-02 | yes | 11 | 1.203052e-02 | 1.156123e-03 | **0.0000e+00** | ( 0.779, 0.000 ) | **1.0000** |
| 1.10 | 0 | **NO** | 11 | 6.361819e-02 | 4.050822e-02 | 0.0000e+00 | **( 1.381, 0.000 )** | 1.0000 |
| 1.12 | −7.4351e-02 | yes | 24 | 1.189041e-02 | 9.539203e-04 | **0.0000e+00** | ( 0.815, 0.000 ) | **1.0000** |
| 1.12 | 0 | **NO** | 15 | 6.523134e-02 | 4.264655e-02 | 0.0000e+00 | **( 1.381, 0.035 )** | 1.0000 |
| 1.15 | −7.1816e-02 | yes | 8 | 1.177335e-02 | 7.029152e-04 | **0.0000e+00** | ( 0.815, 0.000 ) | **1.0000** |
| 1.18 | −6.9463e-02 | yes | 8 | 1.167206e-02 | 4.794681e-04 | **0.0000e+00** | ( 0.850, 0.000 ) | **1.0000** |

### M-46

| | `Ψ` at the located axis | |
|---|---|---|
| projected paraboloid | **1.0000**, axis to 1e-8 | agrees |
| one dof at 29× the peak | **3.4483e-02** = 1/29 | **refuses** |
| one dof at 3.2× | **3.1250e-01** | **refuses** |
| monotone `ψ`, the annulus branch | no O-point of that sense | **refuses** |
| `HighBetaConvergence`'s real solve | **1.0005** | agrees |


## IN-1: a spectral rule fed a second-order Jacobian is a second-order scheme

### M-47

| | |
|---|---|
| **from `q`, pointwise** | **7.03** |
| the trap: the same trapezoid, `ρ′` by differencing positions | 1.97 |
| chord sum | 1.99 |


## The band, and the extension that was chosen by measuring both

### M-48

| `k` | flux Taylor | transfer lift | lift closer at `n = 64` |
|---|---|---|---|
| 1 | 1.797 | **2.298** | **40×** |
| 2 | 2.138 | **3.995** | **1,610×** |
| 3 | 2.138 | **4.848** | **84,695×** |


## The largest lever is one integer, and it is free

### M-49

| `setWalkDepth()` | seconds | fallbacks |
|---|---|---|
| 2 | 0.0895 | 367 |
| **4, the default** | 0.0531 | 183 |
| 8 | 0.0482 | 30 |
| 12 | **0.0262** | **0** |
| 16 | 0.0255 | 0 |


## Threading: available, and blocked by one non-reentrant function

### M-50

| threads | over surfaces | over rays |
|---|---|---|
| 2 | 1.89× | 1.89× |
| 4 | 3.16× | 3.36× |
| 16 | 3.31× | **7.55×** |


## Three things in the kernel survey worth not re-doing

### M-51

| | measured | |
|---|---|---|
| Zernike by Kintner recurrence against per-mode evaluation | **6.3×**, agreeing to 4.4e-15 | **not worth it** — 11.5% of a 0.6% stage |
| a fit at 20,000 points by Vandermonde + GEMM against per-point | **34.7×**, agreeing to 4.4e-16 m | **worth it for IN-6**, which is the many-point case, and not for anything today |
| batched field evaluation by element | **0.90–1.14×** | **do not** |


## §8's separatrix cut: nothing gives out, and that is the finding

### M-52

| `Ψ_N` | 0.50 | 0.80 | 0.90 | **0.95** | 0.98 | 0.99 | 0.995 |
|---|---|---|---|---|---|---|---|
| band share, `h = 0.0425` | 0% | 1% | 15% | **41%** | 80% | 94% | 98% |
| band share, `h = 0.0212` | 0% | 0% | 2% | **28%** | 46% | 78% | 91% |


## The coarea formula is a third cross-check and it was free

### M-53

| | best over four steps |
|---|---|
| plain central difference | 5.161e-06 |
| **Richardson** | **6.835e-09** — a factor of **755** |
| Richardson with `ρ'` **differenced** | 4.695e-05, **6870×** worse |


## The trace matrix is symmetric on the fitted path and is not on the extension path

### M-54

| `k` | `n` | nnz/row | rel `\|A_ij − A_ji\|` | Rayleigh quotients, free subspace |
|---|---|---|---|---|
| 1 | 416 | 9.4 | 2.2e-16 | 200/200 negative, largest −1.44 |
| 2 | 624 | 14.1 | 3.6e-16 | 200/200 negative, largest −1.30 |
| 3 | 832 | 18.8 | 6.9e-16 | 200/200 negative, largest −1.28 |

### M-55

| | fitted | extension |
|---|---|---|
| free–free rel `\|A_ij − A_ji\|` | 2e-16 | **5.4e-1** |
| Rayleigh quotients | 200/200 negative | 200/200 negative, largest −0.898 |


## A quarter of every Newton step was thrown away — fixed, and switched on

### M-56

| `n` | symbolic | numeric | backsolve | symbolic share |
|---|---|---|---|---|
| 12,544 | 20.8 ms | 72.1 ms | 10.9 ms | **22%** |
| 49,664 | 104.2 ms | 325.6 ms | 54.6 ms | **24%** |


## Threading, measured — and the two axes are not the same axis

### M-57

| | Serial | **Threaded** | |
|---|---|---|---|
| `HighBetaConvergence` | 3.22, 3.19 s | **1.34, 1.33 s** | **2.4x faster**, where the gate once made it **1.8x slower** |
| `PedestalConvergence` | 233.3 s | **138.3 s** | 1.69x, the heaviest test in the suite |
| example5 `k = 2, n = 24`, whole solve, 8 threads | 0.3690 s | **0.1313 s** | 2.81x |
| example5 `k = 3, n = 16`, whole solve, 8 threads | 0.2797 s | **0.0927 s** | 3.02x |
| the same at **one** thread | 0.3686 s | 0.3664 s | **1.01x — a wash** |

### M-58

| | `MKL=1` | `MKL=8` | |
|---|---|---|---|
| **Serial** assembly | 0.2797 s | **107.19 s** | 383x, and this is the trap this file records |
| **Threaded** assembly | 0.0927 s | **0.0834 s** | immune |


## Traps

### M-59

| assembly + reduction | `MKL=1` | `MKL=2` |
|---|---|---|
| `k = 2, n = 64` | 0.469 s | 0.440 s — untouched |
| `k = 3, n = 16` | 0.056 s | **2.306 s** |
| `k = 3, n = 32` | 0.216 s | **8.973 s** |


## MEQ against freegs4e, and root selection is the whole difficulty

### M-60

| case | machine | rel `L2` | rel `L∞` |
|---|---|---|---|
| A | test tokamak, classic | **1.5e-04** | 5.4e-04 |
| G | MAST-U | 6.4e-04 | 1.6e-03 |
| B | test tokamak, `ff'`-dominated | 7.1e-04 | 1.5e-03 |
| E | test tokamak, diamagnetic | 2.7e-03 | 5.8e-03 |
| F | DIII-D | 2.9e-03 | 5.2e-03 |
| D | TCV | 4.4e-03 | 7.9e-03 |
| C | MAST | 6.8e-03 | 1.3e-02 |

### M-61

| start | `max ψ` | vs reference |
|---|---|---|
| ramp, amplitude ≤ 0.1 | 1.9313e-03 | −95.7% |
| **seeded from the reference** | **4.525154e-02** | **−0.067%** |
| ramp, amplitude ≥ 0.2 | 5.225179e-02 | +15.4% |
| freegs4e | 4.528205e-02 | — |

### M-62

| changed | MXH fit | MEQ rel `L2` | limited by |
|---|---|---|---|
| 129², 10 harmonics, `k=2 r=2` | 2.060e-04 | 1.519e-04 | the contour |
| 257² | 1.063e-04 | 6.751e-05 | the contour |
| 513² | **1.032e-04** | 5.726e-05 | **the FITTER** — 3% for 4× the grid |
| 16 harmonics | **2.206e-05** | 5.825e-05 | **MEQ** — `L2` did not move |
| MEQ `k=3` | | **8.802e-06** | ? |
| `k=3 r=3` | | 8.744e-06 | not MEQ |
| output grid 257², 513² | | 8.797e-06, 8.800e-06 | not the sampling |

### M-63

| reference | MXH shape fit | MEQ rel `L2` | rel `L∞` |
|---|---|---|---|
| 129² | 2.060e-04 m | **1.519e-04** | 5.403e-04 |
| 257² | 1.063e-04 m | **6.751e-05** | 3.140e-04 |
| | 1.94× | **2.25×** | 1.72× |

### M-64

| `k` | refine | elements | dofs | dofs/ref pt | wall | rel L2 |
|---|---|---|---|---|---|---|
| 1 | 2 | 1735 | 15,615 | 0.94 | 5.25 s | 3.15e-04 |
| 1 | 3 | 7185 | 64,665 | 3.89 | 13.58 s | 1.46e-04 |
| 2 | 1 | 397 | 7,146 | 0.43 | 2.69 s | 1.01e-03 |
| 2 | 2 | 1735 | 31,230 | 1.88 | 6.01 s | 1.52e-04 |
| **3** | **1** | **397** | **11,910** | **0.72** | **3.00 s** | **1.72e-04** |
| 3 | 2 | 1735 | 52,050 | 3.13 | 6.45 s | 1.44e-04 |

### M-65

| grid | `ψ_ax` | `ψ_bnd` |
|---|---|---|
| 129² | 9.483141e-02 | 2.781829e-02 |
| 257² | 9.337971e-02 | 2.649111e-02 |
| 513² | 9.308752e-02 | 2.622462e-02 |

### M-66

| `n` | points | `ψ_ax` | Picard | wall | peak RSS |
|---|---|---|---|---|---|
| 129 | 16,641 | 0.09483140885 | 41 | 4.85 s | 340 MB |
| 257 | 66,049 | 0.09337971414 | 31 | 12.84 s | 1000 MB |
| **513** | **263,169** | **0.09308752051** | **27** | **72.75 s** | **5663 MB** |

### M-67

| | freegs4e, 513² | MEQ, `k = 3`, 26375 el | apart |
|---|---|---|---|
| `ψ_ax` | 9.308752e-02 | **9.307342e-02** | **1.5e-04** |
| `ψ_bnd` | 2.622462e-02 | **2.622580e-02** | **4.5e-05** |


### M-68

`meq::RotatingSource` against `tests/analytic/VaryingCentrifugal.hpp`, over
`r ∈ [0.6, 1.4] × z ∈ {−0.4, 0, 0.4} × ψ ∈ [−0.2, 1.0]`, relative with a floor
of 1. Both closures, because at two species the general path must give the
closed form's answer:

| | `Closure::ClosedForm` | `Closure::RootFind` |
|---|---|---|
| `F` | **7.134e-16** | 1.070e-15 |
| `dF/dψ` | **5.260e-16** | 8.801e-16 |
| `p` | 4.44e-16 | 5.55e-16 |
| `e φ₀` | 2.22e-16 | 5.00e-15 |
| `d( e φ₀ )/dψ` | 1.11e-16 | 6.05e-15 |
| `n_s` | 4.44e-16 | 2.94e-15 |

### M-69

The same fixture's two routes to its own pressure, which is what says its
algebra follows from (96) and (97) rather than agreeing with MEQ's transcription
of the same algebra. Route A bisects (97) with `C` appearing nowhere; route B is
the closed form through `C`, `C′` and `C″`:

| | |
|---|---|
| (97)'s residual at the bisected root | 6.94e-17 |
| `p`, closed form against bisection | 3.33e-16 |
| `dp/dψ`, closed form against a difference of the bisection | 4.04e-12 |
| `d²p/dψ²`, the same | 4.53e-10 |

### M-70

The same fixture mutated, worst deviation from `meq::RotatingSource`. `C` drifts
**2.50×** over the sweep and `C″ = 2.00` exactly, so `C′` is never small:

| mutation | `F` | `dF/dψ` | `p` |
|---|---|---|---|
| conforming | 7.13e-16 | 5.26e-16 | — |
| `C′, C″ → 0` | **7.90e-01** | **8.72e-01** | 4.44e-16 |
| `C′, C″ → −C′, −C″`, which is Li & Zhu's (9) | **1.58e+00** | **1.02e+00** | 4.44e-16 |

`p` does not move under either, which is what says the mutation is confined to
the derivatives rather than being caught upstream of them.

### M-71

`mfem::Mesh::FindPoints` calls, counted on a breakpoint rather than by
`Contour::fallbackLocations`, which cannot see a seed — before and after
`ContourTracer::traceFromAxis()` takes `CriticalPoint::element` as its hint:

| | before | after |
|---|---|---|
| `theTracerClosesAndTheElementWalkDoesNotFallBack`, 6 surfaces | **6** | **0** |
| `SurfaceAverageConvergence`, whole binary | **196** | **1** |
| `FluxSurfaceConvergence`, whole binary | 650 | 527 |
| `FluxGridConvergence`, whole binary | 737 | 64 |

The accounting closes: `SurfaceAverageConvergence` makes 195 `traceFromAxis`
calls and 196 − 195 = 1, that one a genuine counted walk failure in
`fitByAngle()`. What is left elsewhere is the band tests' own deliberate
`Mesh::FindPoints` calls, the test-only three-argument `sampleAt()`, and counted
walk failures — no seeds.

### M-72

XP-0: the X-point of `Soloviev::nstx()` as a zero of `q_h`, located by `sweep()`
with no seed, on `[0.35, 1.15] × [−2.10, −1.30]`. Exact saddle
`( 0.699700, −1.716000 )`. `ratio` is the position error over the pointwise
error of `q` at that same point; `dq/dx` is symmetric with eigenvalues
`+0.899099` and `−0.578731`, so it is trapped in `[ 1.112, 1.728 ]`.

| k | h | `\|x − x_X\|` | rate | L2(q) | `\|q_h − q\|(x_X)` | rate | ratio |
|---|---|---|---|---|---|---|---|
| 1 | 0.20000 | 3.499506e-03 | – | 3.761508e-03 | 2.619010e-03 | – | 1.34 |
| 1 | 0.10000 | 1.226291e-03 | 1.513 | 9.902667e-04 | 9.173171e-04 | 1.514 | 1.34 |
| 1 | 0.05000 | 2.955408e-04 | 2.053 | 2.539655e-04 | 2.178164e-04 | 2.074 | 1.36 |
| 1 | 0.02500 | 1.052146e-04 | 1.490 | 6.429865e-05 | 7.698971e-05 | 1.500 | 1.37 |
| 1 | 0.01250 | 2.290133e-05 | 2.200 | 1.617600e-05 | 1.498590e-05 | 2.361 | 1.53 |
| **1** | **whole** | | **1.814** | | | **1.862** | |
| 2 | 0.20000 | 1.916478e-04 | – | 2.171810e-04 | 1.732927e-04 | – | 1.11 |
| 2 | 0.10000 | 2.096125e-05 | 3.193 | 2.984338e-05 | 1.612502e-05 | 3.426 | 1.30 |
| 2 | 0.05000 | 1.595953e-06 | 3.715 | 3.869881e-06 | 1.022963e-06 | 3.978 | 1.56 |
| 2 | 0.02500 | 5.718299e-07 | 1.481 | 4.904573e-07 | 3.322901e-07 | 1.622 | 1.72 |
| 2 | 0.01250 | 9.738147e-09 | 5.876 | 6.164823e-08 | 8.285009e-09 | 5.326 | 1.18 |
| **2** | **whole** | | **3.566** | | | **3.588** | |
| 3 | 0.20000 | 1.714588e-05 | – | 2.046713e-05 | 9.938342e-06 | – | 1.73 |
| 3 | 0.10000 | 3.585823e-07 | 5.579 | 1.531698e-06 | 3.144368e-07 | 4.982 | 1.14 |
| 3 | 0.05000 | 5.112740e-08 | 2.810 | 1.024806e-07 | 3.316105e-08 | 3.245 | 1.54 |
| 3 | 0.02500 | 1.479874e-09 | 5.111 | 6.562879e-09 | 1.330497e-09 | 4.639 | 1.11 |
| 3 | 0.01250 | 1.408137e-10 | 3.394 | 4.140062e-10 | 9.577455e-11 | 3.796 | 1.47 |
| **3** | **whole** | | **4.223** | | | **4.166** | |

Wanted 1.75 / 2.75 / 3.75, that is `k+1` less the axis study's 0.25 of slack.
The whole-sequence rate over the axis study's shorter `{ 4, 8, 16, 32 }` is
1.685 / 2.796 / 4.500, and at `n = 128` k = 1 reads 1.916 from `n = 4` and
2.017 from `n = 8` while k = 2 drops to 2.884.

### M-73

XP-0's reference, checked before it is used: the prescribed X-point of each
Solov'ev fixture against the saddle of its own `ψ`.

| fixture | prescribed | `\|∇ψ\|` there | saddle of `ψ` | `\|∇ψ\|` | det | moved by |
|---|---|---|---|---|---|---|
| `nstx` | ( 0.699700, −1.716000 ) | 6.40e-17 | ( 0.699700, −1.716000 ) | 9.81e-17 | −0.2547 | 1.11e-16 |
| `iterExample2` | ( 0.883840, −0.704000 ) | 4.48e-16 | ( 0.883840, −0.704000 ) | 1.78e-16 | −0.9098 | 3.33e-16 |
| `nstxExample3` | ( 0.712570, −1.458600 ) | 1.78e-16 | ( 0.712570, −1.458600 ) | 4.16e-17 | −0.5773 | 2.22e-16 |
| **`nstxAsPublished`** | ( 0.699700, −1.716000 ) | **2.97e-02** | **( 0.695811, −1.806937 )** | 6.99e-17 | −0.1712 | **9.10e-02** |

And the audit over `[0.35, 1.55] × [−2.10, 0.30]`, `nstx`, n = 16, which holds
both critical points:

| k | degree | defect | χ | loops | found | Σ index | worst turn | transverse |
|---|---|---|---|---|---|---|---|---|
| 1 | 0 | 0.0e+00 | 1 | 1 | 2 | 0 | 0.088 | no (0.00) |
| 2 | 0 | 2.8e-16 | 1 | 1 | 2 | 0 | 0.070 | no (0.00) |
| 3 | 0 | 4.2e-16 | 1 | 1 | 2 | 0 | 0.070 | no (0.00) |

### M-74

PE-0: `PlasmaEdge` on `Ω_{p,h}`, the elements lying entirely inside the disc,
`n = 8, 16, 32, 64` (26 → 2076 elements), sequence rates. Each row is run
twice — with `λ` the exact trace on `Γ_{p,h}`, and with `λ` on `Γ_p` and
transferred, which is what PE-3 must do. The cut column is M-10, the same
equilibrium on a mesh the plasma edge cuts.

| | `λ` on `Γ_{p,h}` | | | `λ` on `Γ_p`, transferred | | | CUT (M-10) | | |
|---|---|---|---|---|---|---|---|---|---|
| `j` | `k=1` | `k=2` | `k=3` | `k=1` | `k=2` | `k=3` | `k=1` | `k=2` | `k=3` |
| **`ψ*`** | | | | | | | | | |
| 0 | 2.876 | 3.942 | **4.932** | 2.230 | 3.161 | 4.000 | 1.89 | 1.95 | 1.70 |
| 1 | 2.943 | 3.761 | 4.851 | 2.834 | 2.615 | 3.875 | 3.00 | 2.87 | 2.68 |
| 2 | 2.945 | 3.942 | 4.644 | 3.057 | 3.238 | 3.946 | 3.00 | 3.99 | 3.88 |
| **`ψ_h`** | | | | | | | | | |
| 0 | 1.875 | 2.925 | 3.922 | **1.997** | **3.043** | **4.009** | 1.50 | 1.50 | 1.50 |
| 1 | 1.945 | 2.753 | 3.832 | 1.979 | 2.789 | 3.913 | 2.00 | 2.50 | 2.50 |
| 2 | 1.946 | 2.945 | 3.585 | 1.975 | 3.003 | 3.889 | 2.00 | 3.00 | 3.50 |

`ψ*`'s spread over `j`: **0.069 / 0.181 / 0.288** on `Γ_{p,h}` against the cut
mesh's **1.11 / 2.04 / 2.18**. The control — the same solve with `λ` deleted, so
`ψ = 0` on `Γ_p` — is flat at **−0.064** and 1.0e+06 times larger at the finest
mesh (2.9028e-01 against 2.7834e-07).

### M-75

PE-0's premise with the geometry taken out: the same equilibrium on a FITTED
rectangle `r ∈ [0.85, 1.15] × z ∈ [−0.15, 0.15]`, strictly inside the disc, so
no element is cut and no boundary is a staircase.

| `j` | `k = 2`: `ψ` / `q` / `ψ*` | `k = 3`: `ψ` / `q` / `ψ*` |
|---|---|---|
| 0 | 2.998 / 2.998 / **4.003** | 3.999 / 3.996 / **4.998** |
| 1 | 2.987 / 2.983 / **3.989** | 3.992 / 3.992 / **4.993** |
| 2 | 3.000 / 2.991 / **3.995** | 3.981 / 3.982 / **4.944** |

### M-76

PE-0: the transferred route is INTERMITTENT on a smooth circle, so the
mesh-dependent fragility `PLASMA-EDGE-PLAN.md` §9.4 attributes to the X-point's
corner is met with no corner present. `j = 0`, `k = 3`,
`VertexConePath::NumWidened() = 0` and `dist( Γ_{p,h}, Γ_p )/h ∈ [1.02, 1.33]`
at every mesh.

| `n` | 16 | 24 | 32 | 48 | 64 | 96 | 128 |
|---|---|---|---|---|---|---|---|
| `L2(ψ)` | 3.388e-07 | **5.886e-07** | 2.066e-08 | 4.075e-09 | 1.025e-09 | 3.083e-10 | 9.539e-11 |
| `L2(q)` | 4.317e-06 | **1.757e-05** | 3.170e-07 | 1.078e-07 | 2.339e-08 | **1.785e-08** | 7.010e-09 |

Raising the path quadrature to order 24 moves these in the sixth figure —
2.339475e-08 against 2.339478e-08 at `n = 64` — so it is not the quadrature.

### M-77

**What the CUDA build costs MEQ on the CPU path.** Two MFEM installs identical
in every option, TPL and path — SUNDIALS, SuiteSparse, GSLIB, OpenMP,
`MFEM_THREAD_SAFE`, LAPACK, PARDISO, `MFEM_PRECISION=double` — differing only in
`MFEM_USE_CUDA` and the `MFEM_USE_CUDSS` that has to follow it. Nothing in the
run uses cuDSS: the driver refuses that key and UMFPack is the default. The case
is `examples/limited-tokamak.toml`'s geometry at 1601 elements, `k = 3`,
`MKL_NUM_THREADS=1`.

**The two builds produce byte-identical output** — 19 Newton iterations,
`ψ_ax = −1.232304e-01`, the same axis-inside-a-conductor warning — so the
comparison is of cost and not of answers.

| | CUDA = YES | CUDA = NO | |
|---|---|---|---|
| wall, `OMP=4`, five paired runs | 21.39 · 21.46 · 21.56 · 21.75 · 22.32 | 19.84 · 20.06 · 20.39 · 20.44 · 20.55 | **1.07×** |
| user, `OMP=1` | 26.75 · 28.63 | 25.33 · 25.63 | ~1.09× |
| `operator new` calls | **226,710,640** | **171,344,397** | **−55.4 M, −24%** |

The allocation count is exact and load-invariant; the timings are not, and the
run-to-run scatter is 3.6% and 4.3% against a 7% effect — which is why five
pairs are quoted rather than one, and no-CUDA is faster in **five out of five**.

**A first pair taken while another build was finishing read 31.05 s against
26.09 s, 1.19×, and that number is wrong.** It is the same trap this tree
records for every suite time: a timing under contention is a measurement about
the machine.

**The mechanism is only PARTLY the CUDA lambda wrapper.** `operator new`
interposed at its call site shows the wrapper accounting for 24% of the
allocation traffic — `mfem::forall` on the host path constructs an
`__nv_hdl_wrapper_t`, whose `manager::do_call` is separately visible in `perf` at
1.9–2.9% — but **171 M allocations survive with CUDA off**, so the bulk of the
traffic is MFEM's ordinary per-element temporaries and not the wrapper.

### M-78

**Why PARDISO is the default trace solver, and what its MKL threads actually
buy on a whole run.** Two problems on `examples/diverted-tokamak.toml`, both
with `AssemblyMode = "threaded"` and `OMP_NUM_THREADS=8`, best of two, on a
machine gated at `load < 0.7` for four consecutive samples and watched for
foreign processes throughout every run — `intruders=[none]` on all sixteen.
`dv` is the fixture as it ships, 2642 elements at `k = 2`; `rf` is
`RefinementLevels = 1` at `k = 3`, 10701 elements, so the trace system is the
larger share of the run.

Wall seconds, whole driver run including all three output formats:

| | `dv`, MKL=1 | `dv`, MKL=8 | `rf`, MKL=1 | `rf`, MKL=8 |
|---|---|---|---|---|
| UMFPack | 15.33 | 16.22 | 51.07 | 51.68 |
| **PARDISO** | 13.92 | **12.81** | 47.37 | **44.43** |

**Two separate effects, and only one of them is free.** Changing the solver at
`MKL=1` is worth **1.10× and 1.08×** and costs nothing — user time falls too,
35.25 s → 34.12 s and 72.10 s → 70.04 s. Adding MKL threads to PARDISO is worth
a further **1.09× and 1.07×** and is **not** free: user time goes 35.62 s →
62.42 s on `dv`, so a 9% wall gain costs **75% more CPU**. On a machine running
one job that is a good trade and on a shared one it is not.

**UMFPack moves the other way, which is the whole reason the pair of defaults is
chosen together.** MKL threads make it *slower* — 15.33 → 16.22 s, −6% — while
raising its user time 35.25 s → 45.70 s. Best against worst across the table is
**1.27×** on `dv` and **1.16×** on `rf`.

Per-stage, from `NpcThreadScaling --orders 2 --sizes 32 --repeats 2`, threaded
assembly's `solve` column alone — which is where the difference lives, `prep`
being flat at 0.98–1.04 in every configuration:

| case | UMFPack MKL=1 | UMFPack MKL=8 | PARDISO MKL=1 | PARDISO MKL=8 |
|---|---|---|---|---|
| example5 | 0.2070 | 0.2493 (**−20%**) | 0.1634 | **0.1330** (+23%) |
| pedestal | 0.2665 | 0.3430 (**−29%**) | 0.2090 | **0.1697** (+23%) |

So on the solve in isolation the two ends of the table are **1.87×** and
**2.02×** apart, against 1.27× and 1.16× on the whole run — the difference being
everything else a driver run does. **The 1.9× that `setTraceSolver()`'s own
documentation quotes for PARDISO's MKL threads is the isolated
trace-solve figure and does not survive to the wall clock**; 7–9% does.

**Correctness at `MKL=1` is exact and at `MKL=8` it is not, and that is
expected rather than a finding.** At `MKL=1`: threaded against serial assembly
0.000e+00 in both ψ and the flux, identical Newton counts, trace solvers
agreeing to 1.110e-15. At `MKL=8`: 1.443e-15 in ψ and 1.839e-13 in the flux,
because the serial element loop gets every MKL thread while the threaded one
gets none — MKL suppresses itself inside an active OpenMP region — so a blocked
BLAS-3 reassociates against an unblocked loop. `NpcThreadScaling` says so in
the message and does not count it as a failure above `MKL=1`.

**AND THE CHANGE COSTS ONE TEST, WHICH QUALIFIES A CLAIM THIS PROJECT MAKES
ABOUT WHY `TraceSolver` IS SAFE TO EXPOSE.** `PedestalConvergence`'s
`andersonPicardReachesTheSameSolutionAsNewton` goes **red** on PARDISO and green
on UMFPack, isolated by rebuilding with that one line changed and nothing else.
The failure is not a disagreement — it is a **precondition**:

```
n = 16: Newton did not converge, so there is nothing to compare Picard against
```

GS-2 section 4.2 at `k = 1`, `n = 16`, unglobalised Newton, 500 iterations
allowed, tolerance 1e-8 — the coarsest and stiffest case in the suite. Under
UMFPack Newton reaches it; under PARDISO it does not.

**Nothing about this contradicts the 1e-15 agreement, and that is the point.**
The trace solvers agree to 1.110e-15 *per solve*, and they do here too. What
differs is which side of a convergence boundary several hundred compounded
1e-15 differences land on, on a problem where Newton is marginal to begin with —
the same binary prints `Newton alone : FAILED in 60 iterations` for a
neighbouring case under *both* solvers, with `PicardThenNewton` converging in
one stage-2 step, so unglobalised non-convergence here is a property of the
problem rather than of the package.

So the standing claim — *"`AssemblyMode` and `TraceSolver` cannot change the
answer, so a file may choose them freely"* — needs one word: they cannot change
the answer **where both converge**. On a marginal Newton they can change
**whether** it converges. That is weaker than the `Globalisation` hazard it was
contrasted against, which reaches two different discrete solutions 9.4% apart
and reports one of them, but it is not nothing.

### M-79

**A whole NPC solve does not survive an `mfem::Device`, and cuDSS is a
bystander.** `NpcThreadScaling --device cuda`, `k = 2`, `n = 32`, on an RTX 2070
SUPER with CUDA 13.3. The Device configures correctly — `Device configuration:
cuda,cpu`, `Memory configuration: host-std,cuda` — and then:

| | cuDSS | PARDISO | UMFPack |
|---|---|---|---|
| example5 Newton steps, needs 4 | **0/0** | **0/0** | **0/0** |
| pedestal Newton steps | 4/4 | 6/7 | 4/0 |
| flux, threaded vs serial | **4.006e-02** | **4.006e-02** | **4.006e-02** |

**The same 4.006e-02 to four figures on all three, and 0/0 on all three**, which
is what says there is one common fault rather than three — and it says the fault
is not cuDSS, since UMFPack and PARDISO are host solvers that never touch the
device. At `OMP_NUM_THREADS=8` it does not merely mis-solve but aborts, from
several threads at once, in `MemoryManager::CheckHostMemoryType_` reached
through `Vector::AddElementVector` inside `DarcyHybridization::MultNL`'s own
OpenMP region: *"host pointer is not registered"*. `MFEM_USE_EXCEPTIONS` turns
that into a throw, which cannot leave a parallel region, so it lands as
`terminate called recursively`.

**Not localised to `BlockVector`, which was the first hypothesis and is wrong.**
A standalone program doing `BlockVector::Update()` over four offsets, an
`AddElementVector` into a block, and one into a `MakeRef` view of a block gives
identical results under `"cpu"` and `"cuda"`, exit 0 both ways. The trigger is
inside the HDG path.

**The first backtrace taken was of a breakpoint on
`CheckHostMemoryType_` and named MEQ's own constructor.** That function is
called on every `Register_`, most of them successful, so the breakpoint fired on
a benign call and the site it reported was innocent. `catch throw` names the
real one. Same species as every other instrument-not-answer finding in this
tree.

Nothing here is filed upstream: `doc/HDG-DEVICE-OFFLOAD.md` says the integrator
kernels are not built, and the standing rule is not to report findings against
work that has not landed. It also means `apps/meq.cpp`'s refusal of
`TraceSolver = "cudss"` is right for a second and stronger reason than the trade
argument written beside it — **the solve does not survive a Device at all.**

**THE DISCRIMINATING EXPERIMENT, RUN AT UPSTREAM'S REQUEST: A LINEAR SOURCE
FAILS TOO, SO A NON-LINEAR INTEGRAND IS NOT THE CAUSE.** Upstream built their own
reproduction of M-79 and it does **not** reproduce — not under `Device("debug")`
and not under real CUDA on this same RTX 2070 SUPER — and they asked for
`--device cuda` at `OMP_NUM_THREADS=1` with MEQ's source replaced by a linear
one, to separate their leading hypothesis (MEQ's genuinely non-linear
`SourceIntegrator` running a Newton loop) from the other three differences.

`NpcThreadScaling --soloviev`, `k = 2`, `n = 32`, UMFPack, `OMP_NUM_THREADS=1`,
`MKL_NUM_THREADS=1`. **`SolovievEquilibrium` is the right substitution and a
truly `ψ`-independent source would not be**: `dFdPsi` is identically zero so the
problem is linear and Newton takes one step, but it still arrives through
`setSource( Source const & )`, so `nonlinearSource` is non-null,
`usesNonlinearForms()` is true, and the **form routing is unchanged** — the same
`c_bfi_p` branch, `meq::SourceIntegrator` still on the non-linear form's domain,
`AssembleElementVector` still running per element per residual. A `ψ`-independent
source would take the other branch of `buildForms()`, construct `M_p`, and change
which arm of `EnableHybridization()` fires, so a pass would implicate the routing
as readily as the integrand.

| | Newton its | L2 vs exact | max &#124;ψ_h&#124; |
|---|---|---|---|
| `--device cpu` | 2/2 | **1.916e-07** | 2.662508e-01 |
| `--device cuda` | **0/0** | **1.823e-01** | **0.000000e+00** |

**The potential comes back IDENTICALLY ZERO**, and that is a stronger statement
than "wrong". Newton stops at iteration zero because the residual it evaluated
was zero, and `ψ_h` is left at the initial iterate — whose potential block is
zero under NPC, the Dirichlet datum riding in the trace block. So 1.823e-01 is
just ‖ψ_exact‖ and **nothing was computed at all**. The 10⁶ gap is not a
numerical difference.

**Which kills hypothesis 1 and hypothesis 2 together.** The source is linear, so
a Newton loop over a non-linear integrand is not required to produce this. And
`OMP_NUM_THREADS=1` with serial assembly reproduces it, so threaded assembly is
not required either — upstream's reproduction runs the serial routes and passes,
and so does this one on the same hardware, differing only in what MEQ puts on the
forms.

**AND UPSTREAM'S CONTROL CANNOT SEE THIS, BY CONSTRUCTION.** Their reproduction's
domain integrator on the non-linear form "carries a zero coefficient, so it is
arithmetically inert while the routing is yours" — but a contribution whose
correct value is exactly zero is indistinguishable from a contribution that was
never read back. That is the same blind spot they correctly identified in their
own two-arm device cases, surviving into the fix for it. It is also the signature
they already found once on hardware and called suspicious: a trace residual
"non-zero on the host and **exactly zero** on the device… the signature of a
value never read back".

**AND THE CAUSE IS AN ALIAS LEAVING ITS BASE'S VALIDITY FLAGS STALE, FOUND WITH
`mfem::Device( "debug" )` AND FIXED AT FOUR SITES.** That device has device
memory semantics with host arithmetic and `mprotect`s the host page, so a raw
host read of a device-valid buffer is a named fault with a backtrace rather than
a wrong number. `NpcThreadScaling --soloviev --device debug` faults on the first
Newton residual, and each fix moves the fault **forward** to the next site —
which is what distinguishes a fix from a coincidence:

| # | site | whose | the fault it produces |
|---|---|---|---|
| 1 | `DarcyNPCOperator::Mult` | MFEM | `NewtonSolver::Mult`'s `r -= b`, iteration 0 |
| 2 | `DarcyNPCSolver::Mult` | MFEM | `NewtonSolver::Mult`'s `add( x, -c_scale, c, x )` |
| 3 | `GradShafranovSolver::solve()` | **MEQ** | `GridFunction::operator=` |
| 4 | `NpcThreadScaling`'s `copyOf()` | **MEQ** | a raw `g( i )` |

`GetAliasDevicePtr()` ends in `AliasProtect()` on the **base's** host range and
never touches the base's flags — `Memory::SyncAlias`'s own comment says so and
says what is owed. So the base still reads `VALID_HOST`, the next `ReadWrite_`
computes `copy = !( flags & VALID_DEVICE )` = true, and `GetDevicePtr()` memcpys
h2d **out of the page it just protected**. Under CUDA there is no `mprotect`, so
it silently copies the stale host copy over the device buffer the solve had
correctly written.

**Two symptoms, two different links, which is why it looked like one
inexplicable fault.** Site 1 destroys the first residual and produces the
**0/0**; site 3 is MEQ holding three long-lived `MakeRef` aliases —
`darcyFlux`, `potentialGf`, `traceGf` — onto one `BlockVector solution`, and
produces the **identically zero potential**, the solve by then being correct and
nobody reading it back. The evidence for the mechanism is the faulting call
itself: `operator-=` takes `v` first and `*this` second, and

```
    GetDevicePtr h_ptr=0x7fffbd1a5000 bytes=370176 copy=1     <- b, fine
    GetDevicePtr h_ptr=0x7fffbd14a000 bytes=370176 copy=1     <- r, faults
FAULT at 0x7fffbd14a000
```

the fault address is the residual's base `h_ptr` **to the byte**, block 0 being
at offset 0.

With all four fixed, `k = 2`, `n = 32`, UMFPack, `OMP_NUM_THREADS=1`:

| | Newton its | L2 vs exact | max &#124;ψ_h&#124; |
|---|---|---|---|
| `--device cpu` | 2/2 | 1.916e-07 | 2.662508e-01 |
| `--device cuda` | 2/2 | 1.916e-07 | 2.662508e-01 |
| `--device debug` | 2/2 | 1.916e-07 | 2.662508e-01 |

Every digit of the CPU row on both devices, and the `debug` row is the one that
carries the confidence: the whole solve runs with the host pages protected and
never faults. **The MFEM half is two hunks and four lines**, reported as
`HDG-DEVICE-ALIAS-CHAIN-FROM-MEQ.md`; the `xb.Neg()` in site 2 re-protects the
base, so the base sync has to come **after** it and `SyncFromBlocks()` alone
leaves the fault where it was.

**`DarcyHybridization::ReducedGradient()` is reached 12 times in this run and
does not fault**, on a library predating upstream's own fix to it — so upstream's
prediction that MEQ's `GradientMode::Assembled` plus a Dirichlet trace list
would meet it is right about reachability and not borne out here.

**AND NONE OF IT CHANGES THE TRADE.** The integrators are still 46–53% of a step
with no kernels, so `apps/meq.cpp`'s refusal of `TraceSolver = "cudss"` stands —
on the trade argument alone now, the stronger reason beside it having been this
entry.

**Finiteness was checked before any norm and passed**, on both arms and both
assembly modes — `mfem::Vector::CheckFinite()` on the potential and the flux,
taken ahead of every reduction, because `Norml2()` guards with `fabs(v) > 0` and
so reports zero for an all-NaN vector. So the numbers in this entry are real
numbers, including the 4.006e-02 above.

### M-80

**Where a MEQ Newton step's time goes** — the four-leg split
`HDG-NEWTON-STEP-PROFILE-FROM-HDGDEV.md` asks for at its level 1.
`tests/performance/NewtonStepProfile`, **median of five**, `OMP_NUM_THREADS=8`,
`MKL_NUM_THREADS=1`, PARDISO, machine gated at `load < 0.7` for four
consecutive samples and watched throughout: `intruders=[none]`, loadavg 0.10 at
the gate and 0.23 at the end.

The legs are split by call site and do not overlap; `other` is the **remainder**,
total minus the five timed legs. `prepare()` is outside `total` entirely and is
reported beside it, being paid once per mesh rather than once per step.

| case | assembly | its | total s | residual | gradient | factor | backsolve | other |
|---|---|---|---|---|---|---|---|---|
| example5 k=2 n=24 | serial | 4 | 0.2030 | 29.1% | **45.8%** | 15.1% | 2.8% | 7.3% |
| example5 k=2 n=24 | threaded | 4 | 0.0899 | 13.8% | 25.3% | **37.4%** | 6.6% | 16.8% |
| pedestal k=2 n=32 | serial | 5 | 0.4764 | 26.0% | **43.1%** | 14.0% | 3.2% | 12.4% |
| pedestal k=2 n=32 | threaded | 5 | 0.2033 | 12.4% | 24.9% | **34.4%** | 7.3% | 21.2% |
| highbeta k=2 n=16 | serial | 4 | 0.1596 | **48.2%** | 24.7% | 6.5% | 2.7% | 17.5% |
| highbeta k=2 n=16 | threaded | 4 | 0.0725 | 23.1% | 13.4% | 15.3% | 6.5% | **40.9%** |

**THE ANSWER TO THE QUESTION AS ASKED: THE TRACE SOLVE IS NOT OVER HALF IN ANY
ROW.** `factor + backsolve` reads **17.9%, 17.2%, 9.2%** serial and **44.0%,
41.7%, 21.8%** threaded, against upstream's own 54–59% on 128×128 quads and
gffp's ~1%. MEQ is between them and nearer gffp. Element-local work —
`residual + gradient`, which is where both candidate pieces of work live — is
**74.9%, 69.1%, 72.9%** serial and **39.1%, 37.3%, 36.5%** threaded.

**AND THREADING HAS ALREADY TAKEN MOST OF WHAT (a) WOULD TAKE, WHICH IS THE
FINDING THAT BEARS ON THE DECISION RATHER THAN ON THE PROFILE.** The gradient
leg — `GetGradient`, which MFEM documents as "assemble and factor the Jacobian
at x", so `ComputeH()` is inside it — falls from 45.8% to 25.3% on example5 and
43.1% to 24.9% on the pedestal when the element loop becomes an OpenMP region.
(a) removes *work* and threading removes *time*; they compose, but (a)'s ceiling
against a threaded run is `0.7 × gradient` ≈ **9–18%**, where against a serial
one it is 17–32%.

**THE BORDERED ROW IS A DIFFERENT SHAPE AND IT IS THE ONE MEQ ACTUALLY RUNS.**
With `ψ_ax` an unknown the residual leg dominates — 48.2% serial — because the
border spends **four residual evaluations per step against one gradient and one
factorisation**: measured call counts over four steps are residual ×16, gradient
×4, factor ×4, backsolve ×8. So on MEQ's headline configuration the residual
integrators matter more than `ComputeH()` does, which points at (b)'s half of
the split rather than (a)'s. It forces **no re-assembly**: `setNormalisation()`
changes what the source integrator evaluates and the source sits on the
non-linear form, so the next residual picks it up without `darcy->Assemble()`.

**`other` at 40.9% on that row is honest and is not a leg.** It is the bordered
Newton's own cost — the located-axis search, the peak scan, the augmented norm,
the border's dense algebra — plus MFEM's vector arithmetic inside the loop. It
grows as a *share* under threading because the legs shrink and it does not: in
seconds it is 0.0280 serial against 0.0296 threaded, i.e. flat.

**AND THE PINNED `MKL_NUM_THREADS=1` MADE THE TRACE SOLVE'S SHARE AN UPPER
BOUND, NOT A PRODUCTION NUMBER.** Upstream asked for the threaded rows again with
MKL unpinned, on the reasoning that pinning is right for the *serial* rows — it
stops the gradient leg being a measurement of MKL's threading threshold, the
factor of forty at `k = 3` — while also stopping PARDISO threading, which is the
only leg that would use it. Same medians of five, same process, loadavg 0.75 to
1.37:

| threaded, k=2 | trace solve | gradient | residual | other | total s |
|---|---|---|---|---|---|
| example5, pinned | **42.8%** | 26.5% | 13.2% | 17.2% | 0.0919 |
| example5, **unpinned** | **30.1%** | 33.3% | 16.0% | 21.5% | 0.0826 |
| pedestal, pinned | **41.0%** | 25.9% | 12.1% | 20.8% | 0.2249 |
| pedestal, **unpinned** | **23.7%** | 33.4% | 14.5% | 28.0% | 0.1841 |
| highbeta, pinned | 22.7% | 13.0% | 21.3% | 41.8% | 0.0828 |
| highbeta, **unpinned** | 14.1% | 14.4% | 23.8% | 47.2% | 0.0760 |

In seconds the trace solve falls **1.58× / 2.11× / 1.74×**, the whole solve only
8–18%, and residual and gradient are within a few percent in absolute time across
the two arms — which is the threaded-assembly prediction, MKL suppressing its own
threading inside an active OpenMP region. **The trace solve is no longer the
largest leg on any row**; the gradient is, on both fixed-boundary cases. The
Amdahl ceiling on taking the element-local legs to zero rises with it, from
1.66× / 1.61× / 1.52× to **1.97× / 1.92× / 1.62×**.

**LEVEL 2 — `ComputeH()` INSIDE THE GRADIENT LEG — IS THE STABLEST NUMBER IN
EITHER TABLE.** Upstream's `GetComputeHTime()` / `ResetComputeHTime()` /
`GetComputeHCalls()`, statics over a function-local accumulator, wired into
`StepProfile` as a **sub-split of the gradient leg rather than a sixth leg**:
`ComputeH()` runs inside `GetGradient()`, so it is already counted there and
adding it to the others double-counts.

| threaded | ComputeH s | calls | % of gradient | % of whole solve |
|---|---|---|---|---|
| example5, pinned | 0.0177 | 4 | **72.7%** | 19.2% |
| example5, **unpinned** | 0.0197 | 4 | **71.8%** | 23.9% |
| pedestal, pinned | 0.0426 | 5 | **73.0%** | 18.9% |
| pedestal, **unpinned** | 0.0452 | 5 | **73.4%** | 24.6% |
| highbeta, pinned | 0.0074 | 4 | **68.5%** | 8.9% |
| highbeta, **unpinned** | 0.0077 | 4 | **70.2%** | 10.1% |

**68–73% on every row, in both arms, across three sources and three mesh
sizes** — far more stable than anything else here, and it is the discriminator:
the gradient leg is `ComputeH()` plus a quarter of other things. At upstream's
70–75% of `ComputeElementH` by flop count, (a) is bounded at `0.725 × 23.9%` ≈
**17%** of a whole solve on example5, **18%** on the pedestal and **7%** on
highbeta — the top of the 9–18% estimated from the serial rows.

**AND (a)'s PRECONDITION SURVIVES THE EXTENSION PATH, WHICH THIS HARNESS SAID IT
DID NOT.** `NewtonStepProfile`'s level-3 note claimed the curved and
free-boundary cases are the other side of a line and that the shares do not speak
for them. Wrong, and withdrawn: the error is **asymmetry read as state
dependence**, the same one CLAUDE.md already records. `Bnl_data` is written at one
site, `ConstructGrad()`, and only from `grad_Aup`; `HDGExtensionIntegrator` is a
`BilinearFormIntegrator` landing in `A`, assembled once per mesh. So the 17–18%
carries to the problems MEQ exists for, and the note that said otherwise would
have argued against (a) on every one of them.

**AND THE BOUND MAY BE LOW.** The 17–18% is MEQ's 72.5% measured share of
`ComputeElementH`; upstream's own flop count at these dimensions — `na=12`,
`nd=6`, `nc=3`, `nf=3`, `T=9`, 3888 of 4770 — puts it at **82%**. Not worth
arguing to a third digit, but it says which end of the range a delivered saving
would come from, so read 17–18% as a floor rather than a target.

Level 3, the dimensions, so the cacheable fraction can be computed elsewhere:

| case | k | elements | na | nd | nc | nf | trace dofs | ess. trace | numfact/its |
|---|---|---|---|---|---|---|---|---|---|
| example5 | 2 | 1152 | 12 | 6 | 3 | 3 | 5328 | 288 | 4/4 |
| pedestal | 2 | 2048 | 12 | 6 | 3 | 3 | 9408 | 384 | 5/5 |
| highbeta | 2 | 512 | 12 | 6 | 3 | 3 | 2400 | 192 | 4/4 |

`na` is flux dofs per element counting both components, `nd` potential dofs per
element, `nc` trace dofs per face, `nf` faces per element. **The extension path
is OFF on all three**, which is (a)'s stated precondition — fitted meshes, datum
on the mesh boundary, so `HDGExtensionIntegrator` deposits nothing into the flux
block and `Bnl` is not live. MEQ's curved and free-boundary cases are the other
side of that line and this table does not speak for them.

**A note on the statistic.** These are MEDIANS where `NpcThreadScaling` reports
a best-of, and the difference is deliberate: a minimum is the right centre for a
wall time and the wrong one for a share. The spread over five runs is 1.4% to
14% of the total (0.1945–0.2227 on the widest row, example5 serial), which is
the scatter upstream warned makes a single-run share meaningless.

### M-81

**PE-0's shortfall is not the measuring window, and refining discriminates
against BOTH columns.** `tests/convergence/PlasmaEdgeConvergence`'s
`theUncutPlasmaKeepsItsOrderAtEveryVanishingOrder` reports rates endpoint to
endpoint over `{ 8, 16, 32, 64 }`. The obvious reading of its reds is that the
coarse end is pre-asymptotic — at `h = 0.1` the plasma subdomain is **26
elements**, its `dist/h` reads 0.851 against 1.02–1.33 everywhere else, and the
transferred column's per-pair `psi*` rate climbs monotonically, 1.441, 2.415,
2.836 at `j = 0, k = 1`. The whole study was therefore re-run at `n = 128` as
well and every window read out of the five levels:

| window | failing assertions |
|---|---|
| `{ 8, 16, 32, 64 }` — **the one in the tree** | **17** |
| `{ 8, 16, 32, 64, 128 }` | 20 |
| `{ 16, 32, 64, 128 }` | 21 |
| `k <= 2` `{ 16, 32, 64, 128 }`, `k = 3` `{ 16, 32, 64 }` | 22 |
| `k <= 2` `{ 16, 32, 64, 128 }`, `k = 3` `{ 8, 16, 32, 64 }` | 18 |

**Every window is worse, and the two columns run out of room at opposite ends.**

The **fitted** column — the exact trace on `Gamma_{p,h}`, so the staircase with
every other variable removed — floors at `n = 128` exactly as that file's own
comment predicts: `L2( psi* )` reads **3.63e-14** at `j = 0, k = 3` and
**1.89e-14** at `j = 2, k = 3`, and the last-pair rate collapses to 4.764 and
3.147 against `k+2 = 5`. Where it is NOT at the floor its per-pair rates reach
the design order — **2.952** at `j = 0, k = 1` and **3.909** at `j = 1, k = 2`
against `k+2 = 3` and `4`. So the staircase costs a slowly vanishing
perturbation and not an order.

The **transferred** column stalls for a different reason and three orders of
magnitude above the floor. At `j = 0, k = 3`:

| | `n = 64` | `n = 128` | last-pair rate |
|---|---|---|---|
| `L2( q )`, transferred | 2.339e-08 | 7.010e-09 | **1.739** |
| `L2( psi* )`, transferred | 5.477e-10 | 7.842e-11 | **2.804** |
| `L2( psi_h )`, transferred | 1.025e-09 | 9.539e-11 | 3.426 |
| `L2( psi* )`, fitted | 9.866e-13 | 3.631e-14 | 4.764 (the floor) |

`psi*` and `psi_h` land within 20% of each other, so **the post-processing has
bought nothing there** — the transfer's own error is what both are measuring.

**So the two halves of PE-0's red are different findings and one instrument
cannot separate them.** What fixes the fitted column's REPORTED rate is a
`Gamma_{p,h}` that is not a staircase — `meq::AdaptiveDomain`, the companion mesh
of GS-2 §3.3.

**AND THE SENTENCE THAT STOOD HERE — "what caps the transferred column is the
transfer, and no geometry on the plasma side touches it" — IS WRONG.**
→ **[M-83](MEASUREMENTS.md#m-83)**: that column loses no order at all. It carries
a CONSTANT, amplified by the lifting extrapolating a degree-`k` polynomial
outside its element, and the constant is set by `dist/h` — which is precisely
what the companion mesh controls. **The same lever moves both columns, for
unrelated reasons**, and the reading that they needed different cures came from
reading a sequence rate as an order.

### M-82

**XP-2 needs TWO things and neither is sufficient alone: the connectivity fill
actually running, and the plasma edge held fixed within each Newton.**
`tests/convergence/XPointOuter` on `examples/diverted-tokamak.toml` —
freegs4e's `A_testtokamak_classic`, `k = 2` on 2870 elements, the stored
Green's-function guess, the bootstrap pin at `( 1.143144, −0.553965 )`.

| fill live | edge frozen | outcome |
|---|---|---|
| ✗ | ✗ | **FAILED** at the 200 cap, `‖r‖` stalled at 4.79e-04 |
| ✗ | ✓ | converged in 17 — to the **wrong branch**, `ψ_ax` = −7.99e-02 |
| ✓ | ✗ | **FAILED** at the 200 cap |
| ✓ | ✓ | **CONVERGED in 7**, and XP-2 is met |

**EVERY SINGLE-KEY EXPERIMENT LANDS IN A FAILING ROW**, which is why this
survived three sessions of changing one thing at a time.

**THE FILL WAS NEVER RUNNING, AND THAT WAS A FIXTURE DEFECT.**
`setPlasmaSupport()` was called on the plasma source *before* the coil wrapper
was constructed, so the wrapper's own flag stayed false — and the wrapper is
what the solver holds, so `plasmaComponentWanted()` read false and XP-1's flood
fill never ran. `F` was still confined, the inner source's pointwise test being
live, so nothing failed loudly. `apps/meq.cpp` has always done it the other way
round and is correct. With it fixed, the pointwise test carries **752 elements
over 3 components** and the fill keeps the **419** holding the axis: the other
**333, 44% of the candidates, are the private flux region** — a current channel
nobody asked for, on the only diverted case in the tree. §10.3 predicted exactly
that and XP-1 was built for it.

**THE FROZEN EDGE IS §10.5'S OWN PRESCRIPTION APPLIED TO THE SUPPORT** rather
than only to the bounding point: `meq::NormalisedSource::freezePlasmaEdge` holds
the values `insidePlasma()` tests against for the whole of one solve, and
`GradShafranovSolver::setPlasmaSupportFrozen` holds the fill's mask the same
way. `ψ_bnd` is pinned at `ψ_h` at a **saddle**, so it is a non-local functional
of the iterate; letting the support edge follow it sweeps `F`'s on/off region
across the plasma edge every residual, for which the Jacobian carries no surface
term. Freezing the edge alone converges (7, 4, 4, 4, 4); freezing the mask too
costs nothing and buys the tail (7, 4, 3, 2, 2).

**XP-2 IS MET.** The outer fixed point contracts quadratically and finds both
nulls of the double-null machine:

| sweep | its | X-point | step | `ψ_ax` | `ψ_bnd` |
|---|---|---|---|---|---|
| boot | 7 | ( 1.093, −0.602 ) | — | 8.270827e-02 | 3.238153e-02 |
| 1 | 4 | ( 1.093, −0.604 ) | 1.909e-03 | 8.266009e-02 | 3.237931e-02 |
| 2 | 3 | ( 1.093, −0.604 ) | 3.732e-06 | 8.266004e-02 | 3.237932e-02 |
| 3 | 2 | ( 1.093, −0.604 ) | 9.554e-10 | 8.266004e-02 | 3.237932e-02 |
| 4 | 2 | ( 1.093, −0.604 ) | 2.255e-13 | 8.266004e-02 | 3.237932e-02 |

Against freegs4e: the X-point sits **4.378e-04 m** away, `ψ_bnd` agrees to
**0.08%** and `ψ_ax` to **0.07%**; the upper null is found at ( 1.109, 0.796 )
carrying 2.8931e-02 against the reference's 2.891019e-02. `ψ_h` at the
reference axis reads 8.265989e-02 against this solve's own `ψ_ax` of
8.266004e-02, so the border closed on the field it constrains.

**AND A MEASUREMENT THAT WAS TRUE SUPPORTED A CONCLUSION THAT WAS WRONG.**
Turning `ConfineToPlasma` off does make the solve converge — 21 steps, to
`ψ_ax` = −8.09e-02, a *different* equilibrium with `F` on in the vacuum.
Reading that as "the moving support is the hazard" was an inference from a real
one-key experiment to the wrong cause, because a third variable — the fill —
was silently off in **both** of its rows. A one-key experiment separates two
hypotheses only if everything else is where you think it is.

### M-83

**What caps `psi*` on the transferred route is the LIFTING'S EXTRAPOLATION, not
the quadrature — and the quadrature hypothesis is falsified cleanly.** The
transferred datum is `phi_h(x) = g(a(x)) + int_sigma C E_h(u_h).m ds`, entering
the flux equation as `< phi_h, v.n >_e` on each face of `Gamma_{p,h}`. That is two
quadratures: the inner one along the path integrates a polynomial of the flux
degree on a straight-sided element and is exact, and `mfem::HDGExtensionIntegrator`'s
own header says the **outer** one is not, because the foot map `a(x)` is not
polynomial in `x`. MFEM defaults it to `2k+2`, and the data half through
`VectorBoundaryFluxLFIntegrator` to `2k`. Neither had ever been varied — the
experiment on record raised the *path* rule, which is the one that is already
exact.

`GradShafranovSolver::setExtensionQuadratureOrder()` is the knob. Swept at fixed
`h` on `PlasmaEdge( 0 )`, `n = 32 -> 64`, rule orders `{ default, 12, 24, 48, 80 }`:

| | `L2(q)` at `n=64` | rate | `L2(psi*)` at `n=64` | rate |
|---|---|---|---|---|
| `k=3`, transferred, default | 2.339479e-08 | 3.760 | 5.476803e-10 | 4.802 |
| `k=3`, transferred, rule 12 | **8.673375e-09** | **4.637** | 4.335030e-10 | 4.790 |
| `k=3`, transferred, rule 24/48/80 | 8.673373e-09 | 4.637 | 4.335022e-10 | 4.790 |
| `k=3`, **fitted**, every rule | 1.451937e-09 | 3.937 | 9.866429e-13 | 4.949 |

**Three things at once.** The rule **converges by order 12** — 12, 24, 48 and 80
agree to six or seven digits. The fitted column is **bit-identical at every
rule**, which is the null control, since it installs no extension. And the
default really is inexact: `L2(q)` improves **2.7x** at `k=3` and its pair rate
goes 3.760 -> 4.637, so MEQ's flux on the extension path is costing accuracy for
nothing.

**BUT IT DOES NOT CLOSE THE `psi*` GAP**, which is what the sweep was for:
4.335e-10 converged against the fitted column's 9.866e-13, still **440x**. So a
variational crime in the boundary quadrature is NOT what costs `psi*` its order.

**THE GAP IS THE EXTRAPOLATION, AND THE EVIDENCE IS THAT IT TRACKS BOTH `d/h` AND
`k`.** `E_h` is the element's own degree-`k` polynomial evaluated OUTSIDE its
element, along a path of length `d ~ 1.3 h`. Extrapolating a degree-`k`
polynomial a fixed multiple of its own support amplifies its error by a
Chebyshev/Markov factor that explodes in `k`. Ratio of transferred to fitted
`L2(psi*)` at `j = 0`, against the measured `dist(Gamma_p, Gamma_{p,h})/h`:

| `d/h` | `k=1` | `k=2` | `k=3` |
|---|---|---|---|
| 0.851 | 7.8 | 18.4 | 80.0 |
| 1.023 | 21.6 | 62.0 | 289.6 |
| 1.264 | 28.5 | 88.5 | 501.3 |
| 1.256 | 30.0 | 93.2 | 555.1 |

Monotone in `d/h` at fixed `k`, and explosive in `k` at fixed `d/h`. For
comparison `T_k( 1 + 2 d/h )` at `d/h = 1.26` reads **3.5 / 23.8 / 163.9** — the
same shape and the same order of magnitude, which is as much as a bound of that
kind should be expected to give.

**WHY THE RATE LOOKS RIGHT LOCALLY AND WRONG ACROSS THE SEQUENCE.** The
amplification is a CONSTANT, not an order: the transferred `psi*` pair rate at
`n = 32 -> 64` reads 2.84 / 3.87 / 4.79 against `k+2` of 3 / 4 / 5. What the
constant does is start the error one to two orders higher, so a rate read
endpoint to endpoint over a finite range of `h` reads about `k+1` — and at `k=3,
n=128` the transferred error stalls at 7.8e-11 where the fitted column has
reached the double-precision floor.

**THE RATIO DOES NOT PLATEAU CLEANLY ENOUGH TO GATE ON, AND THAT WAS PROPOSED
AND MEASURED OUT.** Since the amplification is `( 1 + 2 d/h )^k` and P.1 bounds
`d/h`, the ratio ought to stop growing, and over the last refinement most rows
oblige -- 1.05 at `j = 0, k = 1` and `2`, 1.08 to 1.17 across `j = 1`, 1.30 at
`j = 2, k = 3`. **`j = 2, k = 2` grows by 2.40**, its transferred pair rates
reading 3.500, 3.520, 2.693, 2.037 while the fitted column is clean at 3.97. So
the ratio inherits this route's own intermittency -- which PE-0's header already
records on a geometry with no corner anywhere -- and any ceiling loose enough to
pass that row is too loose to catch a transfer leaving its regime. It is printed
per mesh and not asserted.

**AND IT PREDICTS THE FIX, WHICH IS A LEVER THIS TREE ALREADY HAS.** If the gap
is set by `d/h`, then shrinking `d/h` shrinks it — and `d/h` is exactly what
GS-2 section 3.3's companion mesh, `meq::AdaptiveDomain`, exists to control.
**[M-81](MEASUREMENTS.md#m-81)** names that same lever for the FITTED column's
staircase, for an unrelated reason. Untested here, and it is the experiment to
run next.

### M-84

**`Gamma_{p,h}`'s re-entrant corners are 225° AND 270°, not "a staircase of
270°", and the band refinement that removes the transfer penalty does NOT make
PE-0's rates pass.** Two measurements, both from
`tests/convergence/PlasmaEdgeConvergence`.

**THE CORNERS.** The background is `Element::TRIANGLE` — right triangles of
45/45/90, six meeting at an interior vertex — so a boundary vertex's interior
angle is a sum drawn from `{ 45°, 90° }` and the reachable re-entrant angles are
225°, 270° and 315°. Summing the owning triangles' angles at every vertex of
`Gamma_{p,h}`:

| interior angle | `n = 16` | `n = 32` | `n = 64` |
|---|---|---|---|
| 90° | 6 | 12 | 22 |
| 135° | 4 | 12 | 16 |
| 180° | 14 | 20 | 54 |
| **225°** | 0 | **8** | **12** |
| **270°** | **4** | **10** | **20** |
| 315° | 0 | 0 | 0 |

**315° never occurs**, and the reason is the selection rule: an element is
dropped if ANY vertex is outside, so a drop propagates to every triangle sharing
that outside vertex and a lone 45° triangle can never be removed while its five
neighbours stay. So the worst angle is **270° at every resolution** and the
singular exponent is `pi/omega = 2/3` — which is what the old "270°" phrasing
assumed, correctly, for the wrong reason. **The re-entrant count grows like
`1/h`** — 4, 18, 32 over `n = 16, 32, 64` — so the corners are a fixed fraction
of the interface however fine the mesh, which is why the fitted column's
shortfall is a slowly vanishing perturbation rather than something refinement
removes.

**THE BAND REFINEMENT.** `meq::AdaptiveDomain`'s constructor is bit-for-bit what
`makePlasmaSubdomain()` already did — same `MarkLevelSetSubdomain`, same
`extraRefine`, same SubMesh — so only `refine()` can move anything. Marking the
elements of `T_h` that own a face of `Gamma_{p,h}`, so step 3's second half
pushes refinement into the band, at `j = 0`, `k = 3`, `n = 32`:

| band | elem | `dist/h` | transferred `psi*` | fitted `psi*` | ratio |
|---|---|---|---|---|---|
| 0 | 490 | 0.894 | 1.527669e-08 | 3.047638e-11 | **501** |
| 1 | 904 | 0.444 | 5.483084e-10 | 2.531174e-11 | 21.7 |
| 2 | 1792 | 0.301 | 1.176300e-09 | 2.530864e-11 | 46.5 |
| 3 | 3524 | 0.147 | **3.801957e-11** | 2.530864e-11 | **1.50** |

**IT REMOVES THE TRANSFER PENALTY, WHICH CONFIRMS M-83'S MECHANISM.** `dist/h`
falls 0.894 -> 0.147 because `refine()` grows `T_h` TOWARDS `Gamma_p` as band
children fall inside, and the amplification collapses with it: the
transferred/fitted ratio goes 28.5 -> 1.03 at `k = 1`, 88.5 -> 1.00 at `k = 2`
and 501 -> 1.50 at `k = 3`. The transferred column stops being distinguishable
from the fitted one, which is exactly what PE-0 set out to show.

**AND IT DOES NOT MAKE THE RATES PASS — IT MAKES THEM WORSE.** Run as a
sequence, three band refinements at every background mesh, `n = 8 -> 64`:

| | `psi` | `q` | `psi*` | wanted |
|---|---|---|---|---|
| `j=0, k=1` fitted | 1.355 | 1.323 | 2.337 | 2 / 2 / 3 |
| `j=0, k=2` fitted | 2.507 | 2.506 | 3.541 | 3 / 3 / 4 |
| `j=0, k=3` fitted | 3.475 | 3.460 | 4.493 | 4 / 4 / 5 |

Every rate falls about 0.5 below `k+1`, **including the fitted column's, which
read 1.875 / 2.925 / 3.922 without the band refinement**. The cause is that a
fixed band count is NOT a self-similar family: element counts go 688 -> 8388
over an eightfold background refinement, 12.2x where a quasi-uniform 2-D family
would give 64x, because conforming triangle refinement PROPAGATES and floods a
large fraction of a coarse domain while staying near the boundary on a fine one.
The coarse member is therefore disproportionately resolved and flattens the
measured rate. Neither `h` nor `1/sqrt(NE)` is an honest abscissa for it.

**AND "THE FAMILY IS NOT SELF-SIMILAR" IS THE WRONG DIAGNOSIS OF THAT**, which
is worth recording because it is the reading the element counts invite. The
family IS fixed-rule -- uniform background at `1/n`, three band refinements,
`h_band/h_interior = 1/8` at every level -- so a rate against the interior `h` is
well defined, and the count growing 12.2x rather than 64x only says the band is a
1-D feature. Read the ERRORS and not the sequence rate and the banded pair rates
are **rising monotonically** toward target:

| `j=0, k=1` fitted | `L2(psi)` | pair | `L2(psi*)` | pair |
|---|---|---|---|---|
| `n=8`, 688 elem | 3.869755e-04 | — | 7.665502e-06 | — |
| `n=16`, 1604 | 2.015991e-04 | 0.941 | 1.976824e-06 | 1.955 |
| `n=32`, 3524 | 7.625434e-05 | 1.403 | 3.834074e-07 | 2.366 |
| `n=64`, 8388 | 2.310155e-05 | **1.723** | 5.943014e-08 | **2.690** |

against `k+1 = 2` and `k+2 = 3`; the same shape at `k = 3`, 3.062 / 3.572 /
**3.792** and 4.078 / 4.594 / **4.806**. That is a **pre-asymptotic head**, not a
broken abscissa: three band refinements at `n = 8` give **688 elements against a
nominal `h` of 0.141**, so the coarse end is far better resolved than its `h`
implies and the error range is compressed.

**AND THE REASON NOT TO REACH FOR IT ANYWAY IS COST, WHICH IS THE SHARP
FINDING.** Band refinement is not a general accuracy improvement — it is
specifically an antidote to the transfer's extrapolation amplification, and it is
a large net loss wherever that amplification is absent. At `k = 3`, `psi*`:

| | elements | `psi*` |
|---|---|---|
| fitted, band 0, `n=32` | 490 | 3.0476e-11 |
| fitted, band 3, `n=8` | 688 | 1.0329e-08 |
| fitted, band 0, `n=64` | 2076 | 9.8664e-13 |
| fitted, band 3, `n=64` | 8388 | 9.0468e-13 |

At comparable cost the **uniform mesh is 339x better**, and at `n = 64` band 3
buys 1.09x the accuracy for 4.0x the elements. It is spending elements resolving
a boundary where the solution is smooth. On the TRANSFERRED route the same
spend pays: band 3 at `n = 32` (3524 elements) beats plain `n = 64` (2076) by
**14.4x**.

**So the two routes want different meshes**, and PE-0 as staged measures them on
one. That, rather than self-similarity, is the design question.

### M-85

**On the route the coupled method must use, band refinement is a clear win, and
one refinement is the efficient point.** `Interface::Fitted` was renamed
`Interface::ArtificialExactTraceOnGammaH` while taking this, because the old name
invited exactly the error the previous round made: reading it as an
implementation one could choose instead of the transfer, and then comparing costs
against a column nobody can run. From PE-3 on `lambda` IS the interface unknown
and exists nowhere but on `Gamma_p` — there is no function of position to
evaluate on `Gamma_{p,h}`. It is a diagnostic, not a route.

Transferred route, `j = 0`, values at `n = 64` over the sequence `{ 8, 16, 32, 64 }`:

| `k` | band | elem | `dist/h` | `L2(psi)` | `L2(psi*)` | `r(psi)` | `r(psi*)` |
|---|---|---|---|---|---|---|---|
| 1 | 0 | 2076 | 1.256 | 2.798137e-05 | 1.972225e-06 | **1.997** | 2.230 |
| 1 | **1** | 2976 | 0.471 | 2.329046e-05 | **2.769459e-07** | 1.515 | **2.694** |
| 1 | 2 | 4752 | 0.281 | 2.310986e-05 | 7.564871e-08 | 1.366 | 2.577 |
| 1 | 3 | 8388 | 0.148 | 2.310168e-05 | 6.088313e-08 | 1.356 | 2.341 |
| 2 | 0 | 2076 | 1.256 | 2.783392e-07 | 4.529251e-08 | **3.043** | 3.161 |
| 2 | **1** | 2976 | 0.471 | 2.411700e-07 | **3.241408e-09** | 2.535 | **3.694** |
| 3 | 0 | 2076 | 1.256 | 1.025308e-09 | 5.476803e-10 | **4.009** | 4.000 |
| 3 | **1** | 2976 | 0.471 | 7.589940e-10 | **8.240617e-11** | 3.493 | 3.876 |

**ONE REFINEMENT BUYS MOST OF IT**: `psi*` improves **7x / 14x / 6.6x** at
`k = 1, 2, 3` for **1.43x the elements**, and its RATE improves with it, 2.230 ->
2.694 and 3.161 -> 3.694. Band 3 costs 4x the elements for a further 4.5x at
`k = 1` and returns worse rates. `psi` itself saturates at band 1 and is
interior dominated thereafter — 2.329e-05, 2.311e-05, 2.310e-05 — so the whole
benefit is to the post-processing, which is what M-83 says the transfer was
polluting.

**AND `psi`'s MEASURED RATE FALLS, WHICH IS THE UNCOMFORTABLE HALF.** 1.997 ->
1.515, 3.043 -> 2.535, 4.009 -> 3.493. Band refinement removes the geometric
error — `Omega_{p,h}` against `Omega_p` — and that benefit is **3.2x at `n = 8`
and 1.2x at `n = 64`**, so it decays FASTER than `psi`'s own `O( h^(k+1) )`.
Removing it lowers the coarse end and leaves the sequence measuring the true
behaviour from a lower baseline. **So band 0's `psi` rate of 1.997 is partly
flattered**: it clears `k+1` with help from a higher-order coarse-end term, and
the cleaner geometry exposes that rather than causing it.

**NO BAND COUNT MAKES PE-0 GREEN.** Band 0 passes `psi` and fails `psi*`; band 1
does the reverse. The errors are uniformly better at band 1 on both variables,
which is what matters for the METHOD; what does not survive is the assumption
that one uniform-sequence rate gate can score a geometry whose error has two
terms converging at different orders.

**WHAT PE-0 NOW GATES, AND WHAT IT STOPPED GATING.** Acting on the above:
`psi` and `q` are reported and no longer rate-gated on `Omega_{p,h}` — they are
interior dominated, and `theSamePlasmaOnAFittedDomainIsCleanAtEveryVanishingOrder`
already scores that claim on a FITTED rectangle at a tighter slack, `k+1` less
0.10 against this file's 0.15. The monotonicity check per pair stays, since that
is the part a coarse-end contamination cannot flatter. `psi*` keeps its rate
gates, on two rows answering two questions: the artificial column at `k+2`, which
is the premise, and the TRANSFERRED column at **one band refinement** at `k+1`,
which is the geometry the coupled method would run on.

PE-0 goes from **10 failures to 4**, and the four are two causes rather than a
mixed bag:

| failing | rate | wanted | cause |
|---|---|---|---|
| artificial, `j=1, k=2` | 3.761 | 3.85 | the staircase corners, M-84 |
| artificial, `j=2, k=3` | 4.643 | 4.85 | the same |
| banded transferred, `j=2, k=2` | 2.741 | 2.85 | `j = 2`, the stiffest edge |
| banded transferred, `j=2, k=3` | 3.372 | 3.85 | the same |

**Both remaining transferred failures are at `j = 2` and neither is at `j = 0`
or `j = 1`**, which is worth knowing before reading them as the transfer: `j` is
the order to which the profiles vanish at the edge, so this is the corner of the
grid where the solution is least regular across `Gamma_p`.

### M-86

**XP-3: the X-point as two unknowns of the same Newton reaches XP-2's answer to
1.9e-14 m, and costs no backsolve.** `tests/convergence/XPointBorder` on
`examples/diverted-tokamak.toml` — the same fixture as M-82, `k = 2`, the stored
Green's-function guess, the bootstrap pin at `( 1.143144, −0.553965 )`, which is
7.07e-02 m from the answer and about two thirds of an element.

The border adds `q_r( x ) = q_z( x ) = 0` to `ψ_bnd − ψ_h( x ) = 0` and makes
`( r_X, z_X )` unknowns. The support is still an outer state — §10.5's
combinatorial one, which is not differentiable — so the sweeps below move the
support and nothing else:

| sweep | its | X-point | moved | `ψ_ax` | `ψ_bnd` | tail order | support |
|---|---|---|---|---|---|---|---|
| 0 | 14 | ( 1.093, −0.603 ) | 7.043e-02 | 8.265144e-02 | 3.237253e-02 | **1.664** | 419/752 |
| 1 | 3 | ( 1.093, −0.604 ) | 7.063e-05 | 8.266004e-02 | 3.237932e-02 | 1.088 | 440/797 |
| 2 | 2 | ( 1.093, −0.604 ) | 3.041e-10 | 8.266004e-02 | 3.237932e-02 | 0.606 | 440/797 |
| 3 | 2 | ( 1.093, −0.604 ) | 1.739e-14 | 8.266004e-02 | 3.237932e-02 | 0.606 | 440/797 |

**THE BORDER CLOSES ON A SADDLE OF THE FIELD IT SOLVED**, which is the whole
stage: an independent `CriticalPointFinder` root find on the converged field puts
the saddle **6.355e-15 m** from where the border left the point, and `ψ_bnd`
agrees with `ψ_h` there to **5.83e-16**. The point sits 0.067 into its element in
reference coordinates, so it is not on a mesh line. Against freegs4e it is
**4.378e-04 m** away — XP-2's own figure to every digit.

**AGAINST XP-2 IN ONE PROCESS, ON ONE FIXTURE:**

| | X-point | `ψ_ax` | `ψ_bnd` |
|---|---|---|---|
| border, XP-3 | ( 1.093103, −0.603529 ) | 8.266004e-02 | 3.237932e-02 |
| outer loop, XP-2 | ( 1.093103, −0.603529 ) | 8.266004e-02 | 3.237932e-02 |
| apart | **1.934e-14 m** | 1.278e-13 of the span | 3.616e-14 of the span |

**IT IS NOT CHEAPER IN NEWTON STEPS AND THE WIN IS STRUCTURAL.** XP-3 spends
14 + 3 + 2 + 2 = **21**; XP-2 spends 7 + 4 + 3 + 2 + 2 = **18** on the same
fixture in the same process. The bootstrap is where it goes: XP-2's is pinned at
a fixed point and is an easy smooth problem, and XP-3's carries the X-point
7.07e-02 m with the line search damping hard for its first eight steps. What XP-3
buys is that **the X-point is no longer an outer state** — one discrete state
where XP-2 has two — and that the answer is a fixed point of the solve rather
than of a loop around it.

**THE OBSERVED ORDER IS THE FIXTURE'S AND NOT THE BORDER'S, AND THE CONTROL IS
WHAT SAYS SO.** XP-3's bootstrap tail reads **1.664** over 14 iterations and
XP-2's — the same bordered system less exactly these two rows — reads **1.667**
over 7. Neither is quadratic. Every warm sweep on this fixture reaches the solve's
floor in two steps, so an order read off one of those measures where the floor is:
both routes read **0.606** on their last sweep, from histories that are three
points long.

**AND THE FIRST STEP OF EVERY BORDERED SOLVE IS TAKEN WITH THE `ψ_ax` BORDER
DECOUPLED, WHICH IS A DEFECT WHOSE REPAIR CHANGES WHICH EQUILIBRIUM IS
REPORTED.** `peakAt()` writes `constraintLocated` and it was last called on the
COLD state when the convergence target was computed — where the flux and
potential blocks are zero, so no axis is found. Measured with a trace at the
border build: `located 0, rowSize 0` at iteration 0 of every solve and
`located 1, rowSize 6` at every iteration after it. Re-establishing the state
before the loop is three lines and does what it should:

| | unrepaired | repaired |
|---|---|---|
| XP-3's bordered solve | 14 steps, orders 6.75 · 0.28 · 4.39 · 1.66 | **10 steps**, orders 1.88 · 1.37 · 1.33 · 1.51, **same answer every digit** |
| XP-2's limiter-pinned solve | 7 steps, `ψ_ax` 8.266004e-02, X-point 4.4e-04 m from freegs4e | **3 steps**, `ψ_ax` **8.052272e-02**, X-point **1.26e-02 m** away and 1.1e-02 from freegs4e |

So the decoupled first step is **conservative**, and on this machine it is what
keeps the iteration inside the physical branch's basin while a full Newton step
leaves it. That makes it a branch-selection question rather than a Jacobian one,
CLAUDE.md's standing rule is that nothing may silently change which equilibrium
is reported, and it is left alone — recorded here so that the next person to
find the stale flag does not repair it without measuring the second row.

### M-87

**XP-4: the diverted machine from a TOML file, against freegs4e.** The driver on
`examples/diverted-tokamak-xpoint.toml` — the same machine as M-82 and M-86, `k = 2`
on 2642 elements, the stored Green's-function guess — and the reference is
freegs4e's own `A_testtokamak_classic` at 129², an up-down asymmetric double null
whose coil currents its control system solved for.

**Nothing about the null is handed over.** MEQ finds it as three rows of its
Newton, from a seed 7.07e-02 m away that the file states deliberately wrong;
freegs4e finds it by a critical-point search on a finite-difference field. The
agreement in its POSITION is therefore a result rather than a precondition, which
is the structural difference from the limited case — `theDriverSolvesALimitedTokamak`
hands BOTH codes the same limiter point, because freegs4e's own contact is a grid
artefact that moves to the other side of the machine between 129² and 513².

| | freegs4e, 129² | MEQ, `k = 2`, 2642 el | apart |
|---|---|---|---|
| X-point | ( 1.093144118, −0.603965084 ) | **( 1.093103369, −0.603529197 )** | **4.4e-04 m** |
| `ψ_ax` | 8.271751445e-02 | **8.266003630e-02** | **6.9e-04** |
| `ψ_bnd` | 3.240412551e-02 | **3.237931762e-02** | **7.7e-04** |
| `I_p` | 2.0e+05 A | 2.0e+05 A | it is the constraint, 8.9e-13 |
| profile scale | 1, by construction | 9.985753e-01 | 1.4e-03 |

**THE SEED IS THE CONTROL.** It sits 7.071e-02 m from the reference and the
border moves the point 7.043e-02 m, ending **161×** closer than it started. Every
other row above is satisfiable by a border that never ran, since the seed is
itself within 7.1e-02 m of the answer.

**THE SUPPORT LOOP SETTLES IN THREE SWEEPS** — 14 + 3 + 2 Newton steps — and the
driver reproduces `XPointBorder`'s in-process answer to **every printed digit**,
including the per-sweep iteration counts. So the TOML route and the API route are
the same computation.

**THE FIELD, AND THE CONDUCTOR MODEL IS THE WHOLE OF THE DISAGREEMENT.** MEQ's
`.nc` — `ψ*` at degree 3 on 129², the band dropped — against the reference
bilinearly interpolated onto it. 5610 of MEQ's 12612 interior nodes are
comparable, the rest lying outside freegs4e's own box; `scale` is
`max |ψ_ref| = 1.051e-01`.

| region | nodes | rel L2 | rel L∞ | worst at |
|---|---|---|---|---|
| everything comparable | 5610 | 7.116e-03 | **2.973e-01** | ( 1.751, −0.615 ) |
| inside P2L's box | 10 | 1.385e-01 | 2.973e-01 | ( 1.751, −0.615 ) |
| inside P2U's box | 10 | 9.474e-02 | 1.924e-01 | ( 1.751, +0.584 ) |
| **minus those two and a 5 cm collar** | 5478 | **5.268e-04** | **1.818e-03** | ( 1.639, −0.034 ) |
| the plasma, `ψ_ref > ψ_bnd` | 1342 | 7.176e-04 | 1.818e-03 | ( 1.639, −0.034 ) |

**Every bit of the 2.97e-01 is inside two conductors**, and it is the conductor
MODEL rather than either solver: `P2L` and `P2U` are freegs4e FILAMENTS — a point
source with a logarithmic singularity — and MEQ's are 0.1 × 0.1 m rectangles
meshed into the domain carrying a uniform current density. No refinement of
either code closes that. `P1L` and `P1U` are `ShapedCoil`s and agree by
construction; they contribute no nodes here at all, sitting at `z = ±1.10`
outside the reference's `[ −1, 1 ]` box. Drop the two filaments and their collar
and the two codes agree over 5478 nodes at **5.3e-04**, with the worst point on
the outboard midplane near the plasma edge rather than anywhere near a coil.

**AND ~7e-04 IS THE FLOOR THIS COMPARISON CAN MEAN.** `fgsref.py` fits a
`UnivariateSpline` to its own analytic profile shape before solving and the fit
MOVES it — 2.514e-05 of the amplitude in `p'` and **1.707e-02** in `ff'` — while
the tables the TOML carries are the analytic shape, for the reasons
`examples/diverted-tokamak.toml`'s header records. So the two codes are solving
sources that differ at the per-cent level in `ff'`, and an agreement much tighter
than what is measured would be evidence of a shared mistake rather than of two
right answers.

The comparison is
`tools/freegs4e-benchmark/compare.py --free-boundary --exclude-box 1.75,-0.60,0.10,0.10 --exclude-box 1.75,0.60,0.10,0.10`;
`--free-boundary` is what says there is no gauge shift and no `-meta.json`, both
codes solving `ψ → 0` at infinity, and `--exclude-box` reports the conductors as
a row of its own rather than dropping them silently.

### M-88

**THE REFERENCE'S OWN ACCURACY, AND IT IS THE FLOOR THIS BENCHMARK HAS ALWAYS
SAT ON.** `tools/README.md` and `docs/validation.rst` both record that MEQ
*saturates* the `freegs4e` comparison — its error stops falling at about
1.4e-04 because that is the reference's accuracy, not MEQ's — and "refine the
reference, not MEQ" has been the standing next step since. This is that
refinement, run: every one of the seven diverted cases at 129², 257² and 513²,
the fine grids seeded from the coarse one (`fgsref.py --seed-from=auto`).

**The grids nest point for point**, which is the reason for the `2ⁿ+1`
convention, so `ψ` is differenced directly on the 129² points all three share.
No interpolation enters, and the coarse grid *is* a subset of the fine one.

| case | \|129−257\| | \|129−513\| | \|257−513\| | ratio |
|---|---|---|---|---|
| A TestTokamak | 3.054e-05 | 1.578e-04 | 1.547e-04 | **0.20** |
| B ff′-dominated | 1.661e-05 | 2.302e-05 | 1.333e-05 | 1.25 |
| C MAST | 5.059e-04 | 7.299e-04 | 2.367e-04 | 2.14 |
| D TCV | 1.076e-03 | 1.172e-03 | 1.691e-04 | 6.36 |
| E diamagnetic | 4.716e-05 | 7.141e-05 | 2.434e-05 | 1.94 |
| F DIII-D | 1.345e-04 | 1.348e-04 | 7.082e-05 | 1.90 |
| G MAST-U | 5.668e-02 | 6.811e-02 | 4.811e-02 | 1.18 |

relative `L2` in `ψ`; `ratio` is `|129−257| / |257−513|`, which for a scheme
converging at order `p` reads `2^p`.

**IT IS ABOUT FIRST ORDER, ON A FOURTH-ORDER OPERATOR.** C, E and F read 1.9 to
2.1 and B reads 1.25; only D reaches 6.4. `GSsparse4thOrder` is what discretises
`Δ*`, so the operator is not what caps this — the **plasma edge** is.
`psi_bndry` comes from a critical-point search on the grid and `core_mask` is a
per-cell boolean, so which cells carry current is resolved to `O( h )` however
accurate the stencil is. That is the same `k ≤ j` edge story FB-4 tells about
MEQ's own moving support, met in a finite-difference code.

**A's 513² IS FURTHER FROM ITS 257² THAN ITS 129² IS** — ratio 0.20, the only
one below 1 — so on that case the finest run is the outlier rather than the
truth, and its `ψ_ax` moves 8.271794e-02 → 8.272642e-02 where 129² → 257² moved
it by 4e-10. **G does not converge at all**: 4.8e-02 between its two finest
grids, with `ψ_bndry` wandering 5.569e-04 → 9.630e-04 → 7.416e-04. MAST-U's
boundary flux is within a factor of a hundred of zero, so its normalised flux is
ill-conditioned in a way none of the others are.

**WHAT IT COSTS TO GET THERE, seeded**, which understates a cold run:

| | 129² cold | 257² seeded | 513² seeded |
|---|---|---|---|
| A | 4.9 s | 13.8 s | 98.1 s |
| C | 3.3 s | 12.2 s | 73.5 s |
| D | 4.5 s | 14.1 s | 89.2 s |
| E | 4.6 s | 18.1 s | **838.8 s** |
| F | 2.7 s | 12.2 s | 69.4 s |
| G | 7.6 s | 18.6 s | 260.9 s |

so a grid doubling costs 3× and then a further 5–7×, for a factor of about two
in accuracy. E's 838.8 s is its Picard hitting the iteration cap at 401 steps
twice rather than converging.

**THE TRANSFERABLE PART.** "Refine the reference" was the right instinct and
does not work: at 20× the cost the reference is twice as good, and on two of
seven cases it is not better at all. **A comparison against `freegs4e` means
something at about 1e-04 relative and nothing below it**, at any grid — so a
MEQ run that agrees to 1e-04 has reached the floor, and the way past it is a
different reference rather than a finer one.

### M-89

**THE RACE, BOTH CODES COLD, ON DIII-D.** `F_diiid_conventional` is the one of
the seven diverted machines MEQ solves from a genuinely cold start — the
conductors at their given currents, a current blob shaped by the design
profiles, `[boundary.xpoint]` seeded at the null the control system was *asked*
for, and nothing from a converged equilibrium. So it is the one case on which a
wall clock compares two solvers rather than two starting positions.

Held equal: both on all 16 threads; MEQ on its own defaults, the PARDISO trace
solver and threaded assembly, which M-78 measures as the fast pair. Both timed
**end to end** — MEQ's includes generating its mesh with gmsh, the solve and
four output formats; `freegs4e`'s includes its boundary matrix, the Picard loop
and its own diagnostics. Neither is a solve time.

| freegs4e, cold | wall | Picard | `ψ_ax` | `ψ_bnd` |
|---|---|---|---|---|
| 129² | 6.4 s | 33 | 3.758545305e-01 | 7.081839835e-02 |
| 257² | 16.0 s | 32 | 3.758520294e-01 | 7.081519830e-02 |

MEQ, cold, measured against the 513² reference of M-88. `rel L2` is over every
comparable node and `no coils` excludes a 5 cm collar round each conductor,
which is M-87's distinction between two solvers and two conductor models:

| rung | elements | dofs | Newton | solve | wall | rel L2 | no coils | `ψ_ax` | `ψ_bnd` |
|---|---|---|---|---|---|---|---|---|---|
| k=1 | 4848 | — | — | — | — | **FAILED** | | | |
| k=1, h/2 | 19507 | — | — | — | — | **FAILED** | | | |
| **k=2** | **4848** | **87264** | **2** | **12.1 s** | **14.2 s** | 5.846e-03 | **7.456e-04** | 3.49e-04 | 5.36e-05 |
| k=2, h/2 | 19507 | 351126 | 2 | 62.4 s | 67.6 s | 5.984e-03 | 6.710e-04 | 3.13e-04 | 2.84e-04 |
| k=3 | 4848 | 145440 | 2 | 22.4 s | 25.5 s | 5.824e-03 | 6.104e-04 | 2.90e-04 | 3.47e-04 |
| k=3, h/2 | 19507 | 585210 | 2 | 113.5 s | 122.7 s | 5.986e-03 | 6.818e-04 | 3.23e-04 | 2.68e-04 |
| k=2, 3 adaptive cycles | 8234 | 148212 | 2 | 45.6 s | 47.4 s | 5.977e-03 | 6.681e-04 | 3.14e-04 | 2.73e-04 |
| k=3, 3 adaptive cycles | 5752 | 172560 | 2 | 69.7 s | 71.5 s | 5.979e-03 | 6.825e-04 | 3.25e-04 | 2.70e-04 |

**EVERY CONVERGED RUNG GIVES THE SAME ANSWER.** 6.7× the degrees of freedom and
8.7× the wall clock, between the cheapest rung and the dearest, move the
agreement by nothing: 6.1e-04 to 7.5e-04 outside the conductors and 2.9e-04 to
3.5e-04 in `ψ_ax`, with no trend in `k` or in `h`. **MEQ is not the error
here** — M-88 measures `freegs4e`'s own accuracy at about 1e-04 relative and
about first order, so this is the floor between two codes and two conductor
models, met at the very first rung.

**SO THE ANSWER TO "WHAT DOES EQUIVALENT ACCURACY COST" IS 14.2 s AGAINST
6.4 s**, at `k = 2` on 4848 elements against 129², with MEQ's 12.1 s solve being
two Newton steps after three plasma-support sweeps. Read it as a factor of two
and not as a ratio: MEQ is C++ against Python and writes four output files where
`freegs4e` writes one, and 2.1 s of MEQ's wall is sampling `ψ` onto a 513² grid
purely so this table could be made.

**ADAPTIVITY BUYS NOTHING HERE, AND THAT IS NOT A CRITICISM OF IT.** Three
cycles at `k = 2` grow the mesh from 4848 to 8234 elements and cost 3.8× the
wall for an agreement that moves from 7.456e-04 to 6.681e-04 — inside the
scatter of the column. A residual estimator refines where the DISCRETISATION
error is, and on this problem there is no discretisation error left to chase
above the floor. Adaptivity is for a problem whose answer is still moving.

**`k = 1` DOES NOT CONVERGE AT ALL**, at either mesh — "no damping of the
bordered Newton step gave a finite residual" on the coarse one and a plain
non-convergence on the refined one. The cold start's margin is thinner than the
converged answer suggests, and degree is part of what buys it.

### M-90

**WHERE MEQ'S TIME GOES ON THE DIII-D MACHINE CASE, AND EVERY STEP TAKEN IN A
DUMBER WAY THAN `freegs4e`.** M-89 puts MEQ at 14.2 s against 6.4 s for a
Python finite-difference Picard on the same equilibrium — a factor of two the
wrong way for a compiled high-order code, which is a bug list rather than a
property of the method.

**THE UNIT THAT MATTERS IS THE NON-LINEAR STEP.** The run takes three plasma
support sweeps of 7, 3 and 2 Newton steps — **twelve steps in ~13 s, 1.1 s
each** — against `freegs4e`'s 33 Picard steps in 6.4 s, **0.19 s each**. So MEQ
is about **six times slower per non-linear step** on a comparable system: 87264
dofs hybridized to a ~22k trace against a 16641-point grid.

`perf record -F 199 --call-graph dwarf`, `OMP_NUM_THREADS=1`, one process. The
profile has no hot spot — the top entry is 8.2% — which is itself the finding:
the cost is per-element dense work and allocation, run too many times.

**1. `reprepare()` RE-SEEDED THE INITIAL GUESS, 49 TIMES. FIXED, AND IT IS 29%
OF THE SOLVE.** The bordered loop re-prepares whenever the exterior
coefficients move, and every call site assigns the iterate from a saved state
on the **next line** — so `projectOntoTrace()` and `seedFluxFromGuess()` were
computed and discarded every time. Counted on this case at `Modes = 10`: 20, 15
and 14 preparations in the three solves.

| | solve | wall |
|---|---|---|
| before | 10.44 s | 11.78 s |
| after | **7.41 s** | **8.40 s** |

`psi_ax` is 3.759851e-01 either way, the normalisation constraint moving in its
last digit only (−6.465e-12 against −6.466e-12). `prepare( bool )` is the fix.

**2. FOURTEEN LINEAR SOLVES PER NEWTON STEP WHERE `freegs4e` DOES ONE.** One for
the Newton direction and thirteen for the borders — ten Gegenbauer modes,
`psi_ax`, `psi_bnd` and the current. Each is a `DarcyNPCSolver::Mult`: an
element loop to reduce, a trace backsolve, an element loop to recover. Twelve
steps is **168 element-loop pairs**. They are solved one at a time and the
thirteen border right-hand sides are all known at once, so one pass over
thirteen vectors and one PARDISO call at `nrhs = 13` replaces thirteen passes.
This is the largest structural item on the list.

**3. THE EXTERIOR RESPONSE COLUMNS ARE FINITE-DIFFERENCED, `Modes + 2`
PREPARATIONS AND RESIDUALS PER SOLVE.** Measured exactly — 6, 8 and 12
preparations per solve at `Modes` 4, 6 and 10. **The coupling is LINEAR in `a`**:
it reaches the residual as a load term, so `dr/da_m` is independent of the
iterate, of the Newton step and of the plasma support. It is already hoisted out
of the Newton loop; it is still rebuilt for each of the three support sweeps,
and could be built once per mesh.

**4. THE CUDA-ENABLED MFEM COSTS 7% OF THE CPU WALL.** M-77, already measured
against an otherwise identical `../mfem/install-nocuda`, byte-identical output:
`mfem::forall`'s host path builds nvcc's host-lambda wrapper, and 55.4 million
of a 30-second run's 226.7 million allocations go away with CUDA off. A
production build has a measured reason to be a separate install.

**5. `L2_TriangleElement::CalcShape` IS 6% OF THE PROFILE.** Basis functions
re-evaluated at every quadrature point of every element on every residual and
every Jacobian. MFEM's `DofToQuad` cache exists for exactly this and the
hybridized path cannot reach it: `DarcyForm` assembles through
`BilinearForm::ComputeElementMatrix()`, the dense per-element host route.

**6. ALLOCATION IS ABOUT 11%** — `operator new` 4.6%, `malloc` 3.5%, `free`
2.7% — in per-element temporaries. Item 4 removes a quarter of it for free.

**7. `AssemblyMode::Batched` IS UNREACHABLE FROM MEQ.** Upstream already has a
batched local factorisation and a batched flux-mass domain assembly; MEQ's enum
carries `Serial` and `Threaded` and never calls `SetLocalFactorMode`. The
element-local dense work this would reach is `mkl_lapack__dgetrs_` 13.4%,
`mfem::Mult( DenseMatrix, ... )` 10.5% and `MultNL` 12.7%. CLAUDE.md's offload
list already ranks it second and calls it "a smaller job than a new kernel and
entirely in this tree".

**8. THREE SUPPORT SWEEPS ARE THREE FULL SOLVES**, 7 + 3 + 2 steps, and the
third exists only to observe that the support did not move. A settling test that
did not need a full re-solve would be worth about 15% of the run.

**9. `Modes = 10` IS MORE THAN THIS CASE NEEDS.** The retained spectrum's tail
reads 1.38e-02 at ten modes and 6.18e-02 at six, against the 1e-01 the run
itself advises; `psi_ax` moves 4.5e-04 between them, which is at M-88's floor.
Six modes costs four fewer preparations and residuals per solve and four fewer
backsolves per step. Four modes is too coarse — 1.8e-03 in `psi_ax`, and the run
says so.

**10. THE `.nc` SAMPLING IS 0.97 s OF 8.40 s** at 129², and was 2.1 s at the
513² M-89 used. It is output rather than solve, and it is reported separately
for that reason.

**THE ARITHMETIC OF THE REMAINING GAP.** After item 1 the solve is 7.41 s for
twelve steps, 0.62 s each, against `freegs4e`'s 0.19 s. Items 2 and 3 are the
per-step multiplier — fourteen solves and `Modes + 2` residuals where the
reference does one of each — and items 4 to 7 are the constant factor on every
element loop. None of them is about the discretisation.

### M-91

**THE BATCHED DEFAULTS, AND A DEVICE FAULT THE DEVICE DEFAULT CANNOT REACH.**
MFEM carries three batched axes and MEQ now reaches all three. What each
defaults to is a different answer, because what each is worth is different:

| | default | why |
|---|---|---|
| `AssemblyMode` | `Batched` **on a device**, `Threaded` otherwise | it is a device mode; its `D` accumulation goes through `AtomicAdd`, which costs on a host where the per-face loop's plain `+=` does not |
| `LocalFactorMode` | `Serial` | upstream's in-situ figure is **10–12% faster at order 2 and 24% slower at order 6**, so the sign depends on the problem |
| `TraceAssemblyMode` | `Serial` | the one mode here that is **not bit exact** |

The device question is `Device::Allows( Backend::DEVICE_MASK )`, asked at
construction — `mfem::Device` is process-wide state configured before the first
`Vector`, so a solver built after it sees it and one built before it could not
have used it. In the driver it applies only where the file named no mode: an
asked-for mode is honoured, an inherited one is chosen.

**AND `AssemblyMode::Batched` NEEDS NEITHER `MFEM_USE_OPENMP` NOR
`MFEM_THREAD_SAFE`**, which a first version of this got wrong by grouping it
with `Threaded`. Read out of `DarcyHybridization::SetAssemblyMode` rather than
assumed: it aborts for `Threaded` and for nothing else. The conservative
grouping was not free — a device build without OpenMP would have been refused
its own default.

**THE FREE-BOUNDARY PATH DOES NOT SURVIVE A DEVICE, AND THAT IS NOT THIS
CHANGE'S DOING.** `--device debug` on the DIII-D machine case faults in the
first support sweep:

    An illegal memory access was made!
    MFEM abort: Error while accessing address 0x...
     ... in function: void mfem::internal::MmuError( int, siginfo_t *, void * )

**The control says it is pre-existing**: with `AssemblyMode` forced to `serial`
and again to `threaded` — so no batched path is taken at all — the same case
faults identically. So the device default is correct and currently unreachable
end to end on a free-boundary run.

It is M-79's class met at a new site. That entry fixed four unsynced host reads
for the **unbordered** NPC solve; the bordered path has host arithmetic of its
own — `rowDot()`'s dot products, the dense bordered matrix, the backtracking —
reading solve outputs through `operator()`, which neither syncs nor
invalidates. `solveWithNormalisation()` already syncs at several producers and
says so; the fault says the list is not complete. **`--device debug` is the
instrument and it names the site with a backtrace**, which is why it costs one
command to take further.

**TAKEN FURTHER, AND IT WAS TWO SITES.** Each fix moved the fault FORWARD,
which is the only thing that distinguishes a fix from a coincidence:

| | fault at | fixed by |
|---|---|---|
| 1 | `refreshXPoint`, XP-3's border | `state.HostRead()`, and `HostReadWrite()` on the two fields it writes one element of |
| 2 | `refreshPlasmaComponent`, called from the driver's support loop | `state.HostRead()` |
| 3 | — the solve converges | |

`mfem::Vector::operator()` is a RAW accessor — it neither syncs nor
invalidates — and both sites read an iterate a device-side trace solve had just
written. **Read and not read-write for the fill**, which is a pure function of
the state; **read-WRITE for `refreshXPoint`'s two fields**, because only that
element's dofs are written there and `HostWrite()`, whose contents are
undefined, would silently lose every other element's.

**THE FIX IS IN THE CONSUMER AND NOT AT THE PRODUCER**, which is the argument
`rowDot()` already carries in the same function: both are funnels that every
path into their border goes through — the Jacobian, each line-search trial, the
fallback step — so a sync there cannot be outgrown by a new caller, where a list
of producers can.

**WITH BOTH IN, THE DIII-D CASE SOLVES UNDER `--device debug` AND THE ANSWER IS
THE CPU'S**: `psi_ax = 3.759851e-01`, `psi_bnd = 7.081394e-02`, the X-point at
( 1.200929, −0.999491 ), two Newton steps, the support settled in three sweeps.

**WHAT IS LEFT IS UPSTREAM'S AND IS NOT ON THE SOLVE PATH.** After the solve,
`postProcess()` aborts:

    Verification failed: ( it != maps->memories.end() ) is false:
     --> host pointer is not registered: h_ptr = 0x7fffaff0c000
     ... in mfem::MemoryManager::CheckHostMemoryType_

    #16 mfem::DarcyHybridization::ReconstructTotalFlux [ clone .cold ]
    #17 mfem::DarcyForm::ReconstructTotalFlux
    #18 meq::GradShafranovSolver::postProcess

`0x7fff...` is the STACK, so the unregistered pointer is a `Vector` over a stack
buffer, and the `.cold` frame says a second throw — from
`MmuHostMemorySpace::Dealloc`, during unwinding — is what turns it into a
`terminate`. MEQ hands that call only heap-backed grid functions and a block
view over its own solution, and the smaller
`examples/free-boundary-halfdisc.toml` runs the same post-processing clean. So
it is recorded here rather than filed, the device offload being explicitly under
construction.


### M-92

**THE OUTPUT STAGE, AND WHERE ITS TIME ACTUALLY GOES.** The question was whether
the NetCDF write loops efficiently. It does not cost anything worth looking at —
one `putVar` per variable, a bulk call, **5 ms** — and the phase around it costs
**0.95 s of a 7.98 s run, 12%**. Timed per writer on
`examples/machine-f-diiid.toml`, `k = 2`, 4848 elements, `OMP_NUM_THREADS=4`,
at two grid sizes because only some of it is grid-shaped:

| | 129² before | 129² after | 513² before | 513² after |
|---|---|---|---|---|
| `postProcess()` | 0.637 | 0.632 | 0.618 | 0.622 |
| `.mesh` + two `.gf` + `_psistar.gf` | 0.049 | 0.053 | 0.048 | 0.054 |
| **`GridSampler` constructor** | **0.118** | **0.0017** | **0.567** | **0.0070** |
| **the three sampling passes** | **0.0135** | **0.0085** | **0.211** | **0.134** |
| the NetCDF write | 0.005 | 0.016 | 0.018 | 0.024 |
| the `.vtu` | 0.127 | 0.130 | 0.126 | 0.131 |
| **output total** | **0.949** | **0.842** | **1.587** | **0.972** |

**THE LOCATOR IS 69x AND 81x, AND IT IS ONE LINE OF ARITHMETIC REPLACING A
NEWTON SOLVE.** `ElementTransformation::TransformBack()` constructs an
`InverseElementTransformation`, searches a point set for a starting guess and
iterates — **1.45 us per call**, measured identically at both grid sizes
(78k calls for 118 ms, 393k for 567 ms). A straight-sided triangle's map is
affine, so its inverse is a 2x2 solve. `Mesh::GetNodes() == nullptr` is exactly
the condition, and a curved mesh keeps the Newton route.

**AND MOST OF THOSE CALLS WERE REJECTIONS**, which is why the win is larger than
the arithmetic suggests. The element's index range is padded by one cell in each
direction against round-off, so at a grid spacing near the mesh's own only about
one candidate in five is inside: at 129² the constructor ran ~78k inversions to
locate 12520 nodes.

**THE SAMPLING PASSES ARE 1.6x AND THE REMAINDER IS `CalcShape`.**
`GridFunction::GetValue()` and `GetVectorValue()` take an element and a point,
so each call re-fetches what only the element decides — an `Array<int>` of dofs,
a gathered `Vector`, a `Vector` for the shape. **Three heap allocations and a
gather per node, for a dot product of ten numbers**, measured at 0.35 us. The
located nodes are now grouped by element (a counting sort into CSR), so that
half is hoisted out of the node loop and what is left is `CalcShape` — which for
`L2_TriangleElement` is a Vandermonde solve and is genuinely per point.

**THE GROUPING IS ALSO WHAT WOULD MAKE THE PASSES PARALLEL.** Grouped by
element, each element writes a disjoint set of node indices, so the outer loop is
a `forall` with no reduction and no contention; a node-major loop over a
scattered element map is neither. Not done — at 129² the whole grid path is now
27 ms.

**THE ANSWER DOES NOT MOVE, AND THE MASK DOES NOT EITHER.** Against the same
binary with the previous sampler, on the same solve:

| | 129² | 513² |
|---|---|---|
| nodes located | 12520/16641, identical | 199651/263169, identical |
| `inside`, `extrapolated` | bit for bit | bit for bit |
| `psi` | 1.5e-14 relative | 5.5e-15 |
| `B_R`, `B_Z` | 1.2e-14, 1.8e-14 | 2.2e-14, 2.4e-14 |

That is round-off from a different order of operations, in the inverse map and
in the dot product, and nothing else.

**THE ACCEPTANCE IS THE NEWTON ROUTE ITSELF**, run over every element in index
order — which is what the constructor would have done — with ownership read out
through a `P_0` field whose value is the element index, so it is reported exactly
rather than approximated. `theAffineInverseLocatesWhatTheNewtonInverseDoes`:
**1681 nodes located identically, 783 of them claimed by more than one element**,
worst sampled difference 4.66e-15. The grid is 41 nodes across 8 cells so that
every fifth line falls on a mesh line: a node in the middle of an element cannot
distinguish the two inverses and a node on a face can, and the case asserts that
it met some.

**AND THE BIGGEST ITEM IN THE OUTPUT PHASE IS NOT GRID-SHAPED AT ALL.**
`postProcess()` is **0.62 s, 67% of the output at the default grid and 8% of the
whole run**, and it does not move with the grid because it has nothing to do with
it. It is one call into `DarcyForm::Reconstruct()`, and the profile says it is
real per-element work rather than one-off setup: `VectorMassIntegrator`,
`VectorDivergenceIntegrator`, `HDGDiffusionIntegrator` and
`NormalTraceJumpIntegrator` all re-assembled at the enriched order, plus a local
solve, for every element. `ReconstructTotalFlux` is 29% of it — and `totalFlux()`
is read by `SolverContract` and by nothing else, but it is an INPUT to
`ReconstructFluxAndPot()`, so it is not droppable. This is the cost of reporting
`psi*` rather than `psi_h`, which is a decision taken on its merits elsewhere;
it is recorded here because a reader timing the output stage will meet it first.
