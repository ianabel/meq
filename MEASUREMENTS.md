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

**THE SECONDS IN THIS TABLE ARE SUPERSEDED BY
[M-97](MEASUREMENTS.md#m-97)**, which re-takes the race on an idle machine and
after M-90's `reprepare()` fix: 7.93 s against 3.99 s, a ratio of 1.99 rather
than 2.22. The accuracy columns are unaffected — the answers are identical to
every digit — and the rows below stand as a measurement of the code and the
conditions of the day.

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
`mfem::Mult( DenseMatrix, ... )` 10.5% and `MultNL` 12.7%.

**THOSE THREE PERCENTAGES ARE `perf report --children`, i.e. INCLUSIVE, AND THE
PARAGRAPH ABOVE THEM QUOTES A *SELF* FIGURE — TWO VIEWS IN ONE ANCHOR.** *"The
profile has no hot spot, the top entry is 8.2%"* is `--no-children`, which is
why item 7 can carry a 13.4% that appears to exceed the stated maximum: nothing
is wrong with either number and they are not comparable. Re-taken on the current
library at `OMP = MKL = 8` to settle it, since upstream asked which view this
was: the top **self** entry is **7.66%** (`mkl_blas_def_dgemm_pst`), while
`mkl_lapack__dgetrs_` reads **self 0.89%, children 8.38%**. **An anchor that
mixes two views of one profile is a measurement about the instrument**, and it
is recorded here rather than silently corrected because the numbers themselves
stand. CLAUDE.md's offload
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

### M-93

**THE TWO LOCATORS, AND WHAT A NODAL FIELD COSTS THE FAST ONE.** MEQ finds
points in a mesh by two unrelated routes — `meq::GridSampler`, which inverts the
element loop and then inverts the element map, and `mfem::FindPointsGSLIB`,
which `meq::FieldTransfer` goes through and whose simplex path splits triangles
into quads internally. They had never been asked the same question.
`theTwoLocatorsAgreeAtTheSamePoints` in `tests/convergence/SamplerConvergence.cpp`
is the comparison: `ψ_h` at `k = 2` on 2048 triangles, sampled at one set of
points by both, with `SetL2AvgType( NONE )` on the gslib side.

| grid | nodes | interior | differing | border | differing | different element | worst border |
|---|---|---|---|---|---|---|---|
| generic, 65² | 4225 | 3648 | **0** | 577 | 60 | 60 | 1.167e-06 |
| face-aligned, 33² | 1089 | 0 | — | 1089 | **0** | 0 | 2.776e-17 |

**WHERE THE POINT IS UNAMBIGUOUS THE TWO AGREE EXACTLY** — 3648 of 3648 interior
nodes, worst difference **7.216e-16**. Two implementations sharing no locating
machinery landing on the same element every time is the check no comparison
against an exact solution can make, because a smooth field is nearly right in
the neighbouring element too.

**EVERY DISAGREEMENT IS A BORDER NODE, AND BORDER NODES ARE COMMON**: 577 of
4225 on a grid with no alignment intended, **one in seven**, because a uniform
grid over a uniform mesh is commensurate with it far more often than a generic
point would be. The worst difference there is **1.167e-06** against **2.032e-06**
of error in `ψ_h` itself, which is the face jump of an L2 field and not an error
in either locator. On the fully aligned grid the two make the *same* tie-break
on all 1089.

**THE DIFFERENT-ELEMENT COUNT IS NOT A STABLE NUMBER AND MUST NOT BE GATED.**
Building the same mesh with an extent of `0.8` rather than `1.4 − 0.6` — one
ulp, `0.7999999999999999` against `0.8` — moves it from **202 to 331** on an
otherwise identical grid. It is decided by which side of a face a node falls on
at the last bit.

**AND THE DEFAULT L2 AVERAGING MEASURES SOMETHING ELSE ENTIRELY.**
`FindPointsGSLIB::Interpolate()` on an L2 field re-interpolates every border
point through an H1 projection unless told not to. Left at the default, **all
577 border nodes differ**, worst 1.93e-06, because the value returned belongs to
neither element's polynomial. **This settles `TODO`'s first sampling question in
the negative**: `SetL2AvgType` is a real rule and it is the wrong one for a
discontinuous field, which is why `meq::FieldTransfer` already sets `NONE` and
records the default costing it 28%, 10% and 12% in L2 at `k = 1, 2, 3`. Picking
a side, which is what `GridSampler` does, is the better rule.

**A MESH ACQUIRES NODES WITHOUT BECOMING CURVED, AND TESTING THE POINTER FOR
NULL COSTS TWO ORDERS OF MAGNITUDE.** `GridSampler`'s affine fast path — M-92's
69× and 81× — keyed on `Mesh::GetNodes() == nullptr`, and
`meq::FieldTransfer`'s constructor calls `Mesh::EnsureNodes()` because
`FindPointsGSLIB` refuses a mesh without one. Constructing a `FieldTransfer`
therefore dropped every later `GridSampler` on that mesh onto the Newton route,
for the life of the process, with the answer unchanged:

| grid | nodes null | after a `FieldTransfer` | | repaired |
|---|---|---|---|---|
| 129² | 1.14 ms | 181.89 ms | **159.4×** | 0.73 ms |
| 513² | 13.54 ms | 1129.31 ms | **83.4×** | 13.08 ms |

The repair is to test the geometry's **degree** rather than the pointer: an
order-1 nodal field on a triangle is the same affine map the vertices are, and
the values are compared against the vertex array rather than assumed equal to
it, since MFEM lets a caller deform a mesh through its nodes. **It was latent in
the driver rather than live** — both `FieldTransfer`s in `apps/meq.cpp` are
built on the stored or previous-cycle mesh, and `previousMesh` is a *copy* — so
nothing enforced it and nothing would have reported it.

### M-94

**WHAT LOCATING THE PLASMA EDGE EXACTLY IS WORTH — `ROADMAP.md` ITEM 12.2's
CEILING.** Item 12.2 proposes locating `{ Ψ > 0 }` with `ψ*` rather than `ψ_h`.
The two fixtures in `tests/analytic/PlasmaEdge.hpp` bracket it without any of it
being built: `MovingPlasmaEdge` reads its support off `ψ_h`, `PlasmaEdge` cuts
at the **exact** edge, and they share one exact solution. `ψ*` is a better
approximation to `ψ` than `ψ_h` is and the fixed cut *is* `ψ`, so the oracle
column is an upper bound on the item.
`locatingTheEdgeExactlyBuysNoOrder` is the case.

| `j` | `k` | ratio moving/fixed in `ψ_h`, `n` = 8 → 32 | in `ψ*` |
|---|---|---|---|
| 1 | 1 | 0.950 · 0.952 · 0.953 | 0.541 · 0.582 · 0.600 |
| 1 | 2 | 1.197 · 1.164 · 1.118 | 0.824 · 0.945 · 0.924 |
| 1 | 3 | 1.035 · 1.002 · 1.002 | 1.015 · 1.092 · 1.065 |
| 2 | 1 | 0.949 · 0.952 · 0.953 | 0.575 · 0.587 · 0.595 |
| 2 | 2 | 1.261 · 1.279 · 1.283 | 0.507 · 0.499 · 0.493 |
| 2 | 3 | **2.983 · 2.434 · 2.025** | 1.178 · 1.139 · 1.033 |

**NO ORDER IS AVAILABLE.** If the level set limited the rate the oracle's
advantage would grow as `h` fell. It does not: over a fourfold refinement the
ratio grows by at most **1.018** in `ψ_h` and **1.121** in `ψ*`, and in the one
row where the oracle is clearly ahead it *falls*, 2.98 → 2.03.

**AND IN THE REGIME MEQ OPERATES IN IT IS WORTH NOTHING.** At `k ≤ j`, where
`ψ*` keeps `k+2`, the exact edge is worth between **0.95× and 1.28×** — the
moving edge is the better of the two in half those rows. The one corner with a
real factor is `j = 2, k = 3`, about **2×**, which is `k > j`, the regime where
the cut already caps the order and which `PLASMA-EDGE-PLAN.md` exists to address
by a different route. `ψ*` would capture only the part of that gap lying between
`ψ_h` and `ψ`.

### M-95

**CAN `Γ` GO BETWEEN THE PLASMA AND THE COILS? NO, ON ALL SEVEN DIVERTED
MACHINES, AND THE OBSTRUCTION IS GEOMETRIC RATHER THAN NUMERICAL.** The
proposal is FB-7's domain reduction taken to its limit: put the artificial
boundary inside every conductor, mesh only the plasma, and carry the whole coil
set as exterior filaments. `meq::ExteriorDtN`'s `Γ` is a semicircle centred at
`( 0, zCentre )` on the symmetry axis, so the question is whether any
`( zCentre, R )` encloses the last closed flux surface and excludes every
conductor. Measured over `zCentre ∈ [ −60, 60 ]`, with each conductor's
**nearest** point and each LCFS point's distance:

| | machine | best `z_c` | margin, m | at `|z_c| = 50` | blocking conductor |
|---|---|---|---|---|---|
| A | testtokamak classic | 0.000 | **−0.368** | −1.753 | P1L, `r` = 1.000, `z` = −1.100 |
| B | testtokamak peaked ff′ | 0.000 | **−0.270** | −1.754 | P1L, 1.000, −1.100 |
| C | MAST spherical | 0.000 | **−1.320** | −2.912 | P1, 0.150, 0.000 |
| D | TCV | −0.002 | **−0.746** | −1.761 | OH.OH1, 0.430, 0.000 |
| E | testtokamak diamagnetic | 0.000 | **−0.255** | −1.749 | P1L, 1.000, −1.100 |
| F | DIII-D | −0.002 | **−1.500** | −2.651 | F1A, 0.861, +0.168 |
| G | MAST-U | 0.000 | **−1.354** | −3.151 | Pc, 0.067, 0.000 |

The margin is `min over conductors of distance − max over the LCFS of distance`
and has to be positive. It is negative everywhere, by 0.26 m at best, and
**sliding the centre away makes it worse, not better** — the last column is the
same quantity at `|z_c| = 50` m. A semicircle centred on the axis that reaches
the plasma's outboard midplane necessarily swallows anything at smaller radius
near the same height, and every one of these machines has such a conductor:
a central solenoid on C, D, F and G, and on A, B and E a divertor coil tucked
under the plasma at `r` = 1.0 where the plasma reaches `r` = 1.78.

**AND THE MIXED CASE — SOME COILS IN, SOME OUT — FAILS ON A DIFFERENT
CONSTRAINT, NOT ON GEOMETRY.** `Γ` may legally sit in any gap between
conductors, keeping the ones inside meshed and carrying the ones outside as
`meq::ExteriorCoilSet`. Taking the smallest such gap above the plasma gives
1.5× to 4.3× in mesh area and puts up to 15 conductors outside — and it is
unusable, because the exterior expansion converges as `( ρ_src/ρ_Γ )^n` and
those radii sit at a **clearance of 1.01 to 1.02**. Requiring a usable
clearance closes it:

| required clearance | coils externalised, A…G | mesh area saved |
|---|---|---|
| 1.10 | 0, 0, 0, 0, 0, 0, **2** | 1.30× (2.08× on G) |
| **1.25** | 0, 0, 0, 0, 0, 0, 0 | **1.00–1.03×** |
| 1.40 | 0, 0, 0, 0, 0, 0, 0 | 0.80× — worse than today |

**The coils cascade.** Clearing the outermost *kept* conductor by 25% puts `Γ`
beyond the next conductor out, which must then also be kept, and so on — the
gaps between conductors in distance-from-origin are smaller than the clearance
the DtN needs. At clearance 1.25 the best admissible `Γ` is the one already
configured: 2.39 against 2.40 on A, 3.38 against 3.40 on F, 3.49 against 3.50
on G. **The existing radii are already at the optimum**, which was not known
before and is the useful half of this measurement.

### M-96

**THE REFERENCES' CONDUCTORS REBUILT AS `ShapedCoil`, AND IT IS WORTH A FACTOR
OF FIVE TO SEVEN ON THE DIII-D COMPARISON.** M-87 attributes the whole of that
comparison's `2.97e-01` relative `L∞` to the conductor model: `freegs4e`'s
`Coil` is a filament with a log singularity in `ψ` and MEQ's is a rectangle of
uniform current density. `tools/freegs4e-benchmark/shaped.py` rebuilds any of
the seven machines with every filament and solenoid replaced by
`freegs4e.shaped_coil.ShapedCoil` on **`conductors.py`'s own rectangle** — the
same one MEQ meshes — preserving circuit topology, `control` flags and the
`current × turns` convention. `fgsref.py --shaped` is the switch.

All seven converge. The rectangles are `shrink_to_fit`'s, so TCV narrows two
conductors and MAST-U six, which is what keeps them from intersecting.

**MEQ against the two references, same solver, same rungs** — `k = 2`, 4848
elements, the reference's own 129² grid, 3672 comparable nodes of which 80 lie
inside a conductor:

| reference | rel `L2` | rel `L∞` | outside `L2` | outside `L∞` | inside `L2` | inside `L∞` |
|---|---|---|---|---|---|---|
| filaments | 4.593e-03 | 1.004e-01 | 9.697e-04 | 4.567e-03 | 3.043e-02 | 1.004e-01 |
| **ShapedCoil** | **9.102e-04** | **1.403e-02** | **5.132e-04** | **2.047e-03** | **5.119e-03** | **1.403e-02** |
| | **5.0×** | **7.2×** | 1.9× | 2.2× | **5.9×** | **7.2×** |

**It improves the nodes OUTSIDE the conductors too, by about a factor of two,
which M-87 did not predict.** The finite-size correction to a filament's field
goes as `( w/d )²` and does not stop at the conductor's edge; `fgsref.py`'s own
note on the limited case estimated it at 1.6e-02 on the near field. So the
conductor model was not only the `L∞` — it was a fifth of the `L2` as well.

**THE INVERSE SOLVE MOVES THE COIL CURRENTS, so the MEQ side must be
regenerated.** These references are produced by `freegs4e`'s inverse solve
against X-point and isoflux constraints, and giving the conductors extent
changes the currents it chooses — up to 38% relative on DIII-D's F3A, though
that is 1.8 kA on a coil set whose others carry 170 kA. A MEQ configuration
built from the filament reference and compared against the shaped one would be
measuring that, so `make_diverted_case.py` is re-run from the shaped `.npz`.

**THE END-TO-END LOOP CLOSES ON DIII-D AND ON NO OTHER MACHINE, AND THAT IS
MEQ'S COLD START RATHER THAN THE CONVERSION.** Every one of the seven was run
both ways, same driver settings, same rungs — which is the only thing that can
tell a conversion defect from a pre-existing one:

| | machine | filament | shaped |
|---|---|---|---|
| A | testtokamak classic | converged, **axis guard fires** | converged, **axis guard fires** |
| B | testtokamak peaked ff′ | FAILED at 3.975e-02 | FAILED at 3.977e-02 |
| C | MAST spherical | FAILED at 5.857e-03 | FAILED at 3.266e-01 |
| D | TCV | FAILED at 9.897e-01 | **converged** |
| E | testtokamak diamagnetic | FAILED at 2.347e-02 | FAILED at 2.339e-02 |
| F | DIII-D | converged | converged |
| G | MAST-U | converged | FAILED at 2.945e-08 |

**Five of the seven behave identically on the two arms.** B and E fail at
residuals agreeing to three figures either way — 3.975 against 3.977e-02, 2.347
against 2.339e-02 — which is about as clean a control as this lineup offers, and
A converges to a `ψ_ax` the axis guard rejects on both. So those failures are
M-89's finding restated: DIII-D is the one of the seven MEQ solves from a
genuinely cold start, and the conversion does not change that.

**The two that differ, differ in OPPOSITE directions**, which is what says this
is branch selection and not accuracy: TCV gains a solve and MAST-U loses one.
Those are two of the three machines carrying a `Solenoid`, and they are exactly
where the conversion stops being a perturbation — the inverse solve moves to a
different root, with **4 sign flips and Δψ_ax of 0.39 of the span on MAST, 6 and
0.10 on MAST-U**, against no sign flips and 1.3e-05 to 7.1e-05 of span on the
four machines without one. Turning a 100-to-324-filament solenoid stack into one
uniform rectangle is a far larger modelling change than giving a filament a 5 cm
box, and `conductors.py` says so in its own header. **MAST-U's "failure" is also
not a divergence** — it is cycling at 2.9e-08 against a target the cold `‖r₀‖`
put just below it.

**So the reference side is finished for all seven and the comparison is
finished for one.** What blocks the other six is MEQ's cold-start basin, which
is a separate open item, and the honest reading of the DIII-D factor is that it
is one machine's number until a seeded start lets the others be measured.

**RE-TAKEN** — → **[M-103](MEASUREMENTS.md#m-103)**. DIII-D's two numbers
reproduce to every digit and the cold-start basin still blocks the rest, but
three rows changed character: A and G-shaped now exit 0 with **no located
O-point at all**, and C is refused on plasma topology rather than failing.

**AND THE TWO REFERENCES MUST NOT SHARE A FILENAME.** They describe different
machines. `--shaped` writes `<case>_shaped.npz` and stamps `conductor_model`
into the file and the machine label — added after a first run overwrote the
filament baseline that every number in M-87, M-88 and M-89 was taken against,
with a file that looked exactly like it.

### M-97

**THE DIII-D RACE RE-TAKEN ON A QUIET MACHINE, BEST OF FIVE.** M-89's race was
run while the machine was not idle and, more importantly, before M-90's
`reprepare()` fix landed, so its absolute seconds are a measurement of different
code under different conditions. Re-taken here with nothing else running —
nothing above 1.6% CPU, load decayed to 1.4 after the test suite was stopped —
`OMP_NUM_THREADS = MKL_NUM_THREADS = 16` on both arms, both cold, both timed end
to end, and MEQ's mesh stamp deleted before every run so gmsh is inside the
number as M-89 defined it.

| | best | median | worst | spread |
|---|---|---|---|---|
| **freegs4e**, 129², cold | **3.985** | 4.055 | 4.795 | 1.20× |
| **MEQ**, `k = 2`, 4848 elements | **7.927** | 8.322 | 8.743 | 1.10× |
| — of which solve | 6.375 | 6.752 | 6.777 | 1.06× |
| — of which gmsh | 0.716 | 0.730 | 1.105 | 1.54× |

**Same work and the same answers, checked rather than assumed**: freegs4e
reports `psi_axis` 0.3758545305 and `psi_bndry` 0.07081839835, every digit
M-89's, and MEQ converges in 2 Newton steps on 4848 elements at degree 2 to
`psi_ax` 3.759851e-01, which is M-90's post-fix value to every digit. So the
accuracy columns of M-89 carry over unchanged and only the clock moved.

**EQUIVALENT ACCURACY NOW COSTS 7.93 s AGAINST 3.99 s, A FACTOR OF 1.99**,
where M-89 recorded 2.22×. The ratio is the robust part; the absolute seconds
were not.

**RE-TAKEN AGAIN AFTER THE THREE SOLVE-LEG COMMITS** —
→ **[M-102](MEASUREMENTS.md#m-102)**, which reads **1.34×** on the same
protocol, with freegs4e's arm moving 1.6% and MEQ's solve leg 1.71×.

**AND `freegs4e` IS THE CONTROL THAT SEPARATES THE TWO CAUSES**, which is the
transferable half. Its code has not changed since M-89, so its whole 6.4 s →
3.99 s — a factor of **1.61** — is the machine and nothing else. MEQ's 14.2 s →
7.93 s is 1.79×, and M-90 already measured its own fix at 11.78 s → 8.40 s wall
on that day's machine, so the decomposition is *the reprepare fix, then the
machine*, in that order of size. **A timing arm whose code did not change is
worth keeping in a race for exactly this reason**: without it there is no way to
tell a code win from a quiet afternoon, and this file's standing rule — that a
suite time is only a measurement on an idle machine — is otherwise unenforceable
after the fact.

**The spreads are also the argument for best-of-five rather than one run.**
freegs4e's first run is the outlier at 4.795 s against 3.985 s for its best,
which is cold file cache; MEQ's spread is 1.10× and its solve leg 1.06×. A
single run of either could have been quoted 20% high.

### M-98

**THE BORDERED STEP'S BACKSOLVES, BLOCKED INTO ONE `ArrayMult`: 1.34× ON THE
SOLVE LEG.** Every backsolve of a bordered Newton step is against one
factorisation — the residual, `ψ_ax`'s column, `ψ_bnd`'s, the prescribed
current's and the `N` exterior Gegenbauer columns — which is the case
`mfem::DarcyNPCSolver::ArrayMult()` exists for and, by its own documentation,
was built for at MEQ's request. On `examples/machine-f-diiid.toml`, `Modes = 10`
at `k = 2` on 4848 elements, that is **14 columns against one Jacobian**.

Controlled by rebuilding with the queue's dispatch flipped to the single-vector
route and nothing else changed, and **run as three interleaved blocks** because
two other agents' jobs were on the machine:

| block | best | median | worst |
|---|---|---|---|
| **A1** blocked | **8.760** | 8.940 | 9.241 |
| **B** single-vector control | 11.531 | 11.863 | 12.133 |
| **A2** blocked | **8.580** | 9.388 | 9.665 |

**1.34× best-of-five** (11.531 / 8.580), 1.26–1.33× on the medians. The blocked
and control distributions do not overlap at all — the worst blocked run, 9.665,
is faster than the best control run, 11.531 — so the effect is larger than the
drift, which is the whole reason for interleaving rather than running A then B.

**The absolute seconds are contended and are not comparable with M-97's 6.375 s
solve**, which was taken on an idle machine. `rate3d` and an `extension` miniapp
were running throughout. Read the ratio.

**IT BEATS THE 1.24× MFEM MEASURED FOR THE SAME ROUTINE, AND THE DIFFERENCE IS
THE TRACE SOLVE.** Upstream's figure is on 4802 triangles at order 2 with 13
columns — very nearly this problem — with the trace leg unblocked *by
construction*, since `UMFPackSolver` does not override `ArrayMult`. MEQ's
default trace solver is PARDISO, which does, so MEQ gets a third leg upstream's
number excludes: one walk of the factors instead of thirteen.

**AND THAT LEG IS ONLY REACHED BECAUSE `TimedSolver` OVERRIDES `ArrayMult`
TOO.** `mfem::Operator::ArrayMult()`'s base implementation loops `Mult()`, so a
decorator that inherits it silently un-blocks whatever it wraps — same answers,
same call counts, no diagnostic. MEQ wraps its trace solver in a timing
decorator on every bordered path, so without the override the blocked trace
solve would have been bought and then thrown away at the wrapper. **A
pass-through decorator is not pass-through for a method it does not name.**

**The answers do not move**: all fifteen runs report `psi_ax` 3.759851e-01 and
`psi_bnd` 7.081394e-02, to every digit the driver prints, and converge in 2
Newton steps. That is agreement to seven figures and not a bitwise claim —
upstream's documentation is explicit that with LAPACK the blocked route's dense
products become GEMM where the single-vector route had GEMV, so **a bitwise
assertion here would pass on a build without LAPACK and fail on this one**.

**WHAT IS NOT BLOCKED, AND WHY IT IS DELIBERATE.** `ψ_bnd`'s column has a
differenced fallback that runs `fieldResidual()`, which reaches `npc->Mult()`
and, on a moving support, `refreshPlasmaComponent()`. `NPCResidual()` is
non-`const` on the hybridization where `NPCReduce()` and `NPCRecover()` are
`const`, and the support decides the source, so neither the factored local
blocks nor the problem they were factored for survive it by contract. That
branch therefore flushes the queue *before* differencing and takes its own
single-column solve afterwards — the exact order the path had before there was
a queue. The analytic column is the default and never reaches it, so the
ordinary path still blocks all 14.

### M-99

**THE THREE BATCHED AXES CROSSED, AND ONLY ONE OF THEM PAYS.**
`[solver] AssemblyMode` × `LocalFactorMode` × `TraceAssemblyMode`, eight cells,
on `examples/machine-f-diiid.toml` — 4848 elements at `k = 2`, `Modes = 10`,
three plasma-support sweeps — at `OMP = MKL = 8` against
`../mfem/install-nocuda`, `RelWithDebInfo`. Each round runs all eight cells back
to back and is **normalised to its own baseline**, so a machine that drifts
between rounds shows as spread within a column rather than as a difference
between rows.

| AssemblyMode | LocalFactor | TraceAssembly | r1* | r2 | r3 | r4* |
|---|---|---|---|---|---|---|
| threaded | serial | serial | 1.000 | 1.000 | 1.000 | 1.000 |
| threaded | serial | **batched** | **0.906** | **0.929** | **0.910** | **0.928** |
| threaded | batched | serial | 0.991 | 0.978 | 1.010 | 1.023 |
| threaded | batched | batched | 0.990 | 0.965 | 0.954 | 0.982 |
| **batched** | serial | serial | 1.246 | 1.171 | 1.180 | — |
| **batched** | serial | batched | 1.183 | 1.151 | 1.192 | — |
| **batched** | batched | serial | 1.268 | 1.174 | 1.264 | — |
| **batched** | batched | batched | 1.735 | 1.159 | 1.222 | — |

\* r1 and r4 began on a verified-idle machine; r2 and r3 ran against another
agent's build, and r4's own last four cells were contaminated by a suite
starting mid-round, which is why the `batched` half of that column is withheld
rather than printed. **The ratios agree across contaminated and clean rounds
alike**, which is the argument for normalising within a round.

* **`TraceAssemblyMode = "batched"` is worth 7 to 9 per cent**, in all four
  rounds, and it was OFF by default. It is now on.
* **`AssemblyMode = "batched"` is a 15 to 27 per cent LOSS**, in all three
  rounds that measured it. It is not a faster `threaded`; it is the mode whose
  atomics pay for themselves on a device and charge rent on a host, which is
  what `apps/meq.cpp` already says and what nothing had yet measured here.
* **`LocalFactorMode = "batched"` is neutral** — 0.978 to 1.023 — and it makes
  the trace win *worse* rather than better when both are on, 0.973 against
  0.918 for the trace axis alone. The two are not additive and the pair is not
  the configuration to ship.

**Every cell converges in 2 Newton iterations to `psi_ax` 3.759851e-01**, to
every digit the driver prints, across all 32 runs. The modes move the clock and
not the answer — with the documented exception that the batched trace assembly
orders a row's columns differently, so the last bits of the trace solve differ
and a bitwise comparison needs `TraceAssemblyMode = "serial"` named explicitly.

**AND TURNING IT ON BY DEFAULT COST A CONVERGENCE — IN THE CONTROL, NOT IN THE
EQUILIBRIUM, WHICH IS NOT WHERE IT WAS LOOKED FOR.** With the member default
flipped too, `theDriverReachesTheExteriorCoupling` throws *"the non-linear
iteration did not converge"*. The obvious reading is that MFEM's batched trace
matrix is delicate on a bordered free-boundary solve. It is not. Swept over
( n = 24, 32, 40 ) × ( degree 2, 3 ) × ( serial, batched ), replicating that
test's library arm exactly:

| n | deg | coupled, serial | coupled, batched | control, serial | control, batched |
|---|---|---|---|---|---|
| 24 | 2 | 8 · 1.006474394e-01 | 8 · 1.006474394e-01 | 39 · 7.907820021e-02 | **did not converge** |
| 24 | 3 | 9 · 1.002834451e-01 | 9 · 1.002834451e-01 | 69 · 1.269081676e-01 | 69 · 1.269081676e-01 |
| 32 | 2 | 9 · 1.003226437e-01 | 9 · 1.003226437e-01 | 12 · 8.724452563e-02 | 12 · 8.724452563e-02 |
| 32 | 3 | 8 · 1.001930645e-01 | 8 · 1.001930645e-01 | 9 · 1.269125268e-01 | 9 · 1.269125268e-01 |
| 40 | 2 | 9 · 1.002709869e-01 | 9 · 1.002709869e-01 | 16 · −9.773981741e-02 | 16 · −9.773981743e-02 |
| 40 | 3 | 8 · 1.001858541e-01 | 8 · 1.001858541e-01 | 8 · 1.269131498e-01 | 8 · 1.269131498e-01 |

Newton steps · `psi_ax`. **The coupled arm is robust in all twelve cells** — 8
or 9 steps, and serial against batched agreeing to every printed digit. The
fragility is entirely the **control**, the same problem with the coupling
removed, and it is not surprising once named: that is a zero datum on an
artificial boundary in the vacuum, posed only to be different from the real
problem, so it is not an equilibrium and nothing makes it well posed. At degree
2 its `psi_ax` reads 7.91e-02, 8.72e-02, **−9.77e-02** across the three meshes —
it changes **sign**, so it is picking a different branch at each resolution
rather than converging — and the n = 40 pair differs between the trace modes in
its last digit, which is the ulp sensitivity visible directly. At degree 3 it is
stable to four figures on all three meshes.

So the fixture moved to **( n = 32, degree 3 )**, where both arms converge
comfortably under both modes, and **both defaults are now `batched`**. The case
reads `7.257e-17` relative over 24130 dofs against the driver, with the control
still `1.069` relative away.

**THE GUESS AMPLITUDE WAS THE WRONG AXIS AND THE EXAMPLE SAYS SO IN WRITING.**
The first attempt at de-fragilising this swept the bump amplitude and required a
majority to converge, by analogy with
`andersonPicardReachesTheSameSolutionAsNewton`'s mesh sweep. That analogy fails
here: `examples/free-boundary-halfdisc.toml` records that this problem has two
non-trivial roots and that the guess *"does not merely start the iteration, it
CHOOSES WHICH EQUILIBRIUM IS REPORTED"*. Measured — at amplitude 0.06 the
library converges in 15 steps to `psi_ax` 9.432e-02 against the driver's
1.006e-01, **105% apart in L2**, which is the documented multiplicity met from
the inside and not a solver defect. Sweeping a branch-selecting parameter in one
arm of a two-arm pin compares different problems. The amplitude is part of the
configuration and stays at the file's 0.1.

**AND THE GUARD AGAINST A SILENT FALLBACK IS ITSELF PARTLY MISLABELLED, WHICH IS
THE FINDING WORTH KEEPING.** Every one of these paths falls back silently in
MFEM, so `GradShafranovSolver` exposes four predicates and the driver prints
`batched paths taken: …`. But `batchedLocalFactorTaken()` returns
`h->CanBatchLocalFactor()` — **can**, not **did** — and reads `yes` in every
cell, including `LocalFactorMode = "serial"`. Its sibling
`batchedLocalSolveTaken()` does track the mode, `NO` to `yes`, so the four are
not consistent with one another and the header's promise that they answer "a
much narrower question than whether it was asked for" holds for some and not
others. Reading that row without opening the implementation would have supported
a wrong conclusion about which path ran. The driver also prints the line **only
when a batched mode is requested**, so the all-serial baseline — the one cell
where what MFEM does by default is most worth knowing — reports nothing.

### M-100

**WHERE A THREADED BORDERED STEP ACTUALLY GOES, AND IT IS NOT WHERE THE FIVE-LEG
SPLIT SAID.** `otherSeconds()` was **70% of a step at eight threads** against 43%
at one — a remainder that GROWS with the thread count is measuring the wrong
thing, because every named leg threads and whatever is left does not. The legs
below are the decomposition that followed, on
`examples/machine-f-diiid.toml`, `../mfem/install-nocuda`, `RelWithDebInfo`,
last solve of three support sweeps.

**The five original legs, 1 thread against 8:**

| leg | `N = 1` | `N = 8` | speedup |
|---|---|---|---|
| residual | 0.324 | 0.051 | **6.4×** |
| gradient | 0.114 | 0.046 | 2.5× |
| trace factorisation | 0.092 | 0.053 | 1.7× |
| trace backsolve | 0.135 | 0.072 | 1.9× |
| **other** | **0.512** | **0.522** | **1.00×** |

**Dead flat.** Inverting Amdahl on the 1.48× the whole solve gets from eight
physical cores puts the serial fraction near 63%.

**AND THE REASON IT WAS UNNAMED IS THAT THE WORK IS IN `const` MEMBERS CALLED
FROM THE LINE SEARCH.** `assemblePlasmaCurrent()`, the row assemblers and
`locateLimiterContact()` are const, and the backtracking re-evaluates
`refreshXPoint`, the plasma current and all ten `transmissionConstraint` sweeps
of Γ per trial damping. Timing at the call sites in the Newton loop missed
exactly those calls; `profile` is now `mutable` and the timers live in the
functions. `refreshXPoint` had no timer at all.

**The constraint leg, decomposed — one slice is the whole of it:**

| slice | s | share | calls |
|---|---|---|---|
| constraint location | 0.230 | 24.1% | 52 |
| — **axis** | **0.204** | **21.3%** | 4 |
| — — **cold full sweep** | **0.199** | **20.9%** | **2** |
| — X-point | 0.004 | 0.4% | 3 |
| — limiter | 0.000 | 0.0% | 0 |
| — `I_p` | 0.015 | 1.7% | 5 |
| — transmission | 0.000 | 0.0% | 0 |

**Two `CriticalPointFinder::sweep()` calls are a fifth of a threaded step.** A
sweep roots every one of 4848 elements; the seeded `tryFindAxisFrom()` is
bounded by a ring count. **100 ms against 2.5 ms, a factor of 40.**

**A HYPOTHESIS MEASURED AND KILLED ON THE WAY**: that
`transmissionConstraint()` and `rowDot()` were the target, because the exterior
rows are stored dense over the whole unknown (~110k) while being supported on Γ
alone. The sparsity is real; the cost is not. Transmission reads **0.000 s**.

**THE FIX IS TO STOP TAKING THE SWEEP, NOT TO MAKE IT FASTER.**
`previousAxisR`, `previousAxisZ` and `havePreviousAxis` were **locals** of
`solveWithNormalisation()`, so the warm route was reached from the second
evaluation of a solve onwards and discarded at the closing brace — and
`PlasmaSupportSweeps` calls `solve()` again on the same solver for an axis that
moves about a millimetre between sweeps. They now seed from `axisLocatedValue`
and its two coordinates, which are already members and already mean *where the
last solve put the axis*.

| | before | after |
|---|---|---|
| cold sweeps per solve | 2 | **1** |
| axis leg | 0.204 s, 21.3% | **0.116 s, 13.5%** |
| constraint location | 0.230 s, 24.1% | **0.146 s, 16.9%** |

`psi_ax` 3.759851e-01 either way, 2 Newton iterations, same support sweeps. The
constraint residual moves in its last printed digit, −6.465e-12 → −6.466e-12:
the ring search and the full sweep reach the same critical point by different
root-find paths.

**THE REMAINING SWEEP CANNOT BE FIXED BY A SEED, AND IT DID NOT NEED ONE.** It
is `peakAt( coldState, … )`, the cold reference for the convergence target,
where the flux and potential blocks are zero and no axis exists — so it roots
4848 elements in order to fail.

**Warm starting it is impossible twice over**, which is worth recording because
a seed is the obvious thing to reach for after the fix above. The seeded
`tryFindAxisFrom()` needs an extremum to follow and this field has none, so it
fails and falls through to the same cold sweep; and on a fresh solve
`havePreviousAxis` is false anyway, since nothing has yet succeeded on a real
iterate. A seed cannot help a search whose field is identically zero.

**What it needed was an early-out, and the justification is local and
provable.** Under NPC `traceBase` is `blockOffsets[ 2 ]`, so the loop three
lines above the call has zeroed the flux block *and* the potential block
outright. `q_h` and `psi_h` are therefore identically zero on that state by
construction, the axis search cannot succeed, and `peakAt()`'s fallback — a
nodal maximum over an all-zero potential block — is exactly `0`. The call is
replaced by its own answer, with `constraintLocated = false`,
`plasmaSeedElement = -1` and the source's normalisation set, which is precisely
the state `peakAt()` leaves behind when the locate fails.

| | with the sweep | early-out |
|---|---|---|
| cold sweeps per solve | 1 | **0** |
| axis leg | 0.113–0.125 s | **0.002 s** |
| constraint location | 0.120 s, 17.7% | **0.023 s, 3.9%** |
| `psi_ax` | 3.759851e-01 | 3.759851e-01 |
| constraint residual | −6.466e-12 | −6.466e-12 |

Four interleaved rounds, `build-nocuda` at `OMP = MKL = 8`. **The axis leg
saving is ~0.116 s and reproducible to the millisecond; the solve-level delta
is consistent in sign across all four rounds but reads 0.04 to 0.48 s**, so on a
machine at load 6.8 its magnitude is not resolvable and only the leg figure is
claimed. Against a ~4 s solve that is about 3%. Note the legs sum to ~0.74 s
against that solve: they instrument the Newton loop, not the assembly ahead of
it.

**AND THIS IS NOT THE BRANCH-SELECTION CHANGE THE SAME CALL SITE WARNS ABOUT**,
which is the distinction that makes it safe to take. That warning is about
adding `peak = peakAt( unknown, s, … )` *after* this block so the flag describes
the ITERATE — which makes `constraintLocated` TRUE at iteration 0, couples the
axis row on the first step, and moves the diverted machine to a different
equilibrium. This change keeps the flag FALSE and therefore keeps the decoupled
first step exactly; it declines to spend a full-mesh sweep discovering that
falseness. The suite is 50/51 with only `PlasmaEdgeConvergence`'s deliberate red,
so every bordered case in the tree — `HighBetaConvergence`,
`FreeBoundaryCoupling`, `XPointBorder`, `LimiterCurve`, `DriverAcceptance` —
confirms the locate was failing there too.

**WHAT THIS SAYS ABOUT CCSZ.** `meq::SourceIntegrator` lives in the residual
leg, which is **8.0% of a threaded step on this case**. The 23–28% Amdahl
ceiling in `INTERPOLATORY-HDG-PLAN.md` §3.4 comes off M-80's high-beta fixture —
fixed boundary, no border, no exterior modes, no line-search constraint
re-evaluation — and does not carry to a diverted free-boundary run. Making the
source integrator free would buy under 8% here.

### M-101

**MEQ TAKES THE CONDENSATION CACHE, AND `LocalFactorMode = "batched"` IS WHAT
TURNS IT OFF.**

`DarcyHybridization` holds a `CondensationCache` — `A⁻¹Bᵀ`, `B A⁻¹Bᵀ`, `A⁻¹Cᵀ`,
`B A⁻¹Cᵀ − E`, `C A⁻¹Bᵀ + G` and `C A⁻¹Cᵀ`, kept across a Newton loop — and
upstream sizes it against MEQ's own M-80: 1.85× on `ComputeH()` at order 2, about
**11% of a whole MEQ solve at k = 2**. It applies only under
`LocalOpType::PotNL`, whose member `lop_type` defaults to `FullNL` and is
inferred at `Finalize()`. MEQ asked neither predicate, so whether it was getting
any of this was unknown rather than assumed.

`GradShafranovSolver::fluxMassIsPrefactored()` and `condensationCacheTaken()`
now expose the two, and `meq --profile` prints them. On the DIII-D machine case,
`build-nocuda` at `OMP = MKL = 8`:

| configuration | flux mass prefactored (PotNL) | condensation cache |
|---|---|---|
| as shipped | **yes** | **yes** |
| `LocalFactorMode = "batched"` | yes | **NO** |

**Both read `yes` by default, so MEQ was already getting it** — the question was
worth asking and the answer is that nothing is being left on the table.

**AND THE SECOND ROW IS THE FINDING.** `CanCacheCondensation()` has seven ways
to return false, so it is not a tautology, and one of them is
`lfac_mode == LocalFactorMode::Batched` — *"the batched factorisation owns AiBt
and the Schur complement itself and hands them in; caching underneath it would
be two owners of one buffer."* So that mode is a **trade** and not a knob that
does nothing: it buys a batched local factorisation and pays the condensation
cache. That is the mechanism behind M-99's flat `LocalFactorMode` axis, which
was recorded there as neutral-to-negative with no explanation — two effects of
opposite sign, not one absent effect. It also means the two modes must never be
read as independent when the cross above is re-run.

The other six refusals are a block non-linear (0,1) gradient (`Bnl_data`), a
non-linear flux mass, a non-linear face constraint on either the potential row
or the full system, and an empty mesh. MEQ trips none of them, which is the
positive half of the first row.

### M-102

**THE DIII-D RACE, RE-TAKEN AFTER THE THREE SOLVE-LEG COMMITS: 1.34× WHERE
M-97 READ 1.99×.** Same protocol as M-97 — both arms cold, timed end to end,
`OMP_NUM_THREADS = MKL_NUM_THREADS = 16` on both, MEQ's mesh stamp deleted
before every run so gmsh is inside the number, best of five — with one
addition: the rounds are **interleaved**, freegs4e then MEQ five times, so a
drift in the machine lands on both arms rather than on whichever ran second.
The machine was idle, load 0.02 at the start and nothing else above 8% CPU.

| | best | median | worst | spread | M-97's best |
|---|---|---|---|---|---|
| **freegs4e**, 129², cold | **3.922** | 4.121 | 4.251 | 1.08× | 3.985 |
| **MEQ**, `k = 2`, 4848 elements | **5.245** | 5.344 | 5.473 | 1.04× | 7.927 |
| — of which solve | 3.733 | 3.827 | 3.940 | 1.06× | 6.375 |
| — of which gmsh | 0.684 | 0.700 | 0.715 | 1.04× | 0.716 |
| — of which output | 0.788 | 0.803 | 0.815 | 1.03× | — |
| — of which setup | 0.019 | 0.019 | 0.020 | 1.05× | — |

**The answers did not move, which is what makes this a timing measurement.**
freegs4e reports `psi_axis` 0.3758545305 and `psi_bndry` 0.07081839835 in 33
Picard steps, every digit of M-89's and M-97's; MEQ converges in **2 Newton
steps on 4848 elements at degree 2** to `psi_ax` 3.759851e-01, which is M-90's
post-fix value to every digit.

**AND `freegs4e` IS AGAIN THE CONTROL THAT SEPARATES CODE FROM MACHINE.** Its
code has not changed since M-89, and its arm moved 3.985 → 3.922, **1.6%** —
so the machine is the same machine. MEQ's moved 7.927 → 5.245, a factor of
**1.51**, and its solve leg 6.375 → 3.733, a factor of **1.71**. All of that is
code.

**The 1.71 decomposes about as the three commits predict.** M-98 measures its
blocked `ArrayMult` at **1.34×** on this exact leg of this exact case; the
remaining 1.71/1.34 = **1.28×** is the axis carry, the cold sweep answered in
closed form and the batched trace assembly — and the axis work was sized at
about 3% of a solve on its own, so the batched trace assembly is carrying most
of that residue. The two legs that did **not** change are the check on the
decomposition: gmsh reads 0.684 against 0.716 and neither commit touched it.

**The spreads are tighter than M-97's and the reason is the file cache.**
M-97 records freegs4e's first run as an outlier at 4.795 s against 3.985 s for
its best; here the smoke round that warmed the cache was run and discarded
before the five, so freegs4e spreads 1.08× and MEQ 1.04×. A single run of
either would now be quoted within 8%.

### M-103

**THE SEVEN DIVERTED MACHINES, BOTH CONDUCTOR MODELS, RE-TAKEN: DIII-D IS
STILL THE ONLY ONE THAT CLOSES, AND THREE ROWS CHANGED CHARACTER.** M-96's
end-to-end loop re-run on the shipped configurations, nothing tuned per
machine, `OMP = MKL = 16`, both arms — which is the only thing that can tell a
conversion defect from a pre-existing one.

| | machine | M-96 filament | filament | M-96 shaped | shaped |
|---|---|---|---|---|---|
| A | testtokamak classic | converged, axis guard fires | **exit 0, NO AXIS**, 3.060e-01 | converged, axis guard fires | **exit 0, NO AXIS**, 4.146e-01 |
| B | testtokamak peaked ff′ | FAILED at 3.975e-02 | **exit 2**, no finite damped step | FAILED at 3.977e-02 | **exit 2**, same |
| C | MAST spherical | FAILED at 5.857e-03 | **exit 1, REFUSED: plasma contains the symmetry axis** | FAILED at 3.266e-01 | **exit 2**, no finite damped step |
| D | TCV | FAILED at 9.897e-01 | **exit 2**, no finite damped step | converged | **converged**, 3.908e-02 |
| E | testtokamak diamagnetic | FAILED at 2.347e-02 | **exit 2**, same | FAILED at 2.339e-02 | **exit 2**, same |
| F | DIII-D | converged | **converged**, 4.593e-03 | converged | **converged**, 9.102e-04 |
| G | MAST-U | converged | **converged**, 1.394e-02 | FAILED at 2.945e-08 | **exit 0, NO AXIS**, 1.440e+01 |

**M-96's headline survives and so does its control.** DIII-D closes on both
arms and reproduces M-96's two numbers **to every digit** — 4.593e-03 against
the filament reference and 9.102e-04 against the shaped one, the 5.0× M-96
attributes to the conductor model. B and E still behave identically on the two
arms, and TCV still gains a solve under the shaped conductors. What blocks the
other machines is still MEQ's cold-start basin.

**"EXIT 0" IS NOT "CONVERGED", AND THREE ROWS ARE THE DIFFERENCE.** A, and
G under shaped conductors, exit 0, print a Newton convergence and write a
`.nc` that carries **no `axis_normalised_flux` at all** — the driver's own
warning is *"no interior extremum of `psi_h` was found anywhere on the mesh,
so `psi_ax` is not the flux at a magnetic axis"*. `compare.py`'s axis gate has
three verdicts and this is `absent`, which **warns and lets the comparison
run**, by the deliberate choice recorded in `axis_status()`: a wall-hugging
annulus is a real thing to look at, and an old `.nc` must not become
unreadable. So these rows are refused by a reader and not by the harness.

| `.nc` | `axis_normalised_flux` | located at |
|---|---|---|
| `machine-f-diiid` | 1.0000 | ( 1.769, −0.000 ) |
| `machine-f-diiid-shaped` | 1.0000 | ( 1.769, −0.000 ) |
| `machine-g-mastu` | 1.0000 | ( 1.004, 0.000 ) |
| `machine-d-tcv-shaped` | 1.0000 | ( 0.908, −0.002 ) |
| `machine-a-testtokamak` | **absent** | — |
| `machine-a-testtokamak-shaped` | **absent** | — |
| `machine-g-mastu-shaped` | **absent** | — |

**AND THE SHARED SYMPTOM IS THE X-POINT, NOT THE AXIS.** Both no-axis rows put
the located X-point more than a metre from its seed — A at 1.096 m, G-shaped at
1.257 m, the latter landing at ( 0.053, −0.008 ), which is essentially the
symmetry axis. A also reports *"the plasma support took 4 sweeps and DID NOT
SETTLE"*. The absent O-point is the consequence; the X-point constraint
collapsing onto a spurious point is the cause, and it is the thing to
instrument.

**G UNDER SHAPED CONDUCTORS IS THE ONE THAT GOT WORSE.** M-96 records it
*failing* at 2.945e-08 and reads that as cycling against a target the cold
`‖r₀‖` put just below it — an honest failure. It now spends **187 s of solve**,
exits 0, and reports a field **1.44e+01** relative `L2` from the reference. A
run that used to fail now reports a non-equilibrium as a success, which is
strictly the worse of the two. Whether the solve changed or the criterion now
accepts what it used to cycle on is a bisect against M-96's tree and is not
answered here.

**C's REFUSAL IS THE ONE THAT GOT BETTER, AND IT IS A DIFFERENT FINDING FROM
M-96's.** Under filaments MAST now **converges in 3 Newton steps** with a
properly located axis — ( 1.0156, −0.0000 ), normalised flux 1.0000 — and is
then refused, exit 1, by the topology guard: *"THE PLASMA CONTAINS THE SYMMETRY
AXIS. This is not a large error, it is the wrong topology."* `psi_bnd` came out
**−2.356703e-02** against `psi_ax` +5.896474e-02, so −`psi_bnd`/( `psi_ax` −
`psi_bnd` ) = +2.8555e-01 is positive and the plasma reaches `r = 0`. M-96 had
this as a plain non-convergence at 5.857e-03; the solver now gets far enough
for the guard to name what is wrong with the answer.

**THE FAILURES ARE ALL ONE FAILURE AND IT IS THE BORDERED NEWTON.** B, D and E
on both arms, and C under shaped conductors, exit **2** with
*"no damping of the bordered Newton step gave a finite residual — `psi_ax`"*.
That is a single mode, not five, and it is `[source] Normalised = true` making
`psi_ax` an unknown that no globalisation MEQ has can carry from these cold
starts. M-96 quoted a residual for each; this mode reaches none, since no
damped step is finite.

**A HARNESS NOTE, because it cost a wrong table before it was caught.** Reading
`compare.py`'s `absent` verdict as acceptable turns the three no-axis rows into
plain "converged" rows carrying 3.06e-01, 4.15e-01 and 1.44e+01 — and the
relative errors are the only thing that looks wrong, which on a benchmark whose
conductor model is known to cost 2.97e-01 is exactly the number a reader
excuses. **The verdict is the datum and the norm is not.** Anything sweeping
these cases must refuse `absent` as well as `bad`.

### M-104

**THE FORWARD FREE-BOUNDARY SOLVE IS VERTICALLY UNSTABLE, AND ONE SCALAR OF
FEEDBACK COSTING LESS THAN A MILLIAMP MAKES IT SOLVABLE.** Every reference in
`tools/freegs4e-benchmark/` is an INVERSE solve, whose currents are an output;
MEQ solves forward, so the standing rule wants forward against forward. Running
freegs4e forward — its own converged currents FROZEN, its own profiles, `I_p`
prescribed, **no constraint**, seeded at its own converged equilibrium — is the
same problem MEQ is given, and on the elongated machines it does not converge.

**THE RESIDUAL HISTORY CANNOT SEE WHY, WHICH IS THE POINT OF WATCHING THE
AXIS.** `rel` measures how far `psi` moved between passes and says nothing
about where the plasma went, so a solve failing to resolve and one walking out
of the machine look identical in it. DIII-D, filament conductors, 129²:

| it | rel | `Ip_logic` L | nO | nX | `psi_axis` | `Raxis` | **`Zaxis`** |
|---|---|---|---|---|---|---|---|
| 0 | 5.51e-08 | 1.000000 | 5 | 2 | 3.758545e-01 | 1.768397 | **−0.000155** |
| 20 | 6.85e-09 | 1.000000 | 5 | 2 | 3.758545e-01 | 1.768397 | **−0.000155** |
| 40 | 6.89e-08 | 1.000001 | 5 | 2 | 3.758547e-01 | 1.768397 | **−0.000156** |
| 60 | 6.81e-07 | 1.000007 | 5 | 2 | 3.758563e-01 | 1.768399 | **−0.000164** |
| 80 | 7.10e-06 | 1.000075 | 5 | 2 | 3.758736e-01 | 1.768411 | **−0.000248** |
| 100 | 7.84e-05 | 1.000811 | 5 | 2 | 3.760626e-01 | 1.768545 | **−0.001148** |
| 120 | 8.50e-04 | 1.008999 | 5 | 2 | 3.781240e-01 | 1.770031 | **−0.011280** |
| 140 | 6.56e-03 | 1.083508 | 5 | 3 | 3.952126e-01 | 1.779643 | **−0.104115** |
| 160 | 2.54e-02 | 1.446969 | 6 | 5 | 4.519223e-01 | 1.749739 | **−0.555220** |
| 180 | 4.65e-01 | 0.307380 | 5 | **0** | −1.164425e-01 | 0.859546 | **+0.171989** |

**A FACTOR OF TEN EVERY TWENTY PASSES IN BOTH COLUMNS AT ONCE, WHICH IS ONE
UNSTABLE EIGENVALUE AND NOT A CONVERGENCE PROBLEM.** `Zaxis` and `rel` grow in
lockstep at **1.122 per iteration** while `Raxis` moves by 1.4e-05 m over the
first eighty — a rigid VERTICAL displacement. The seeded state is a fixed point
to eight digits and holds for forty passes before the mode is large enough to
see. At about iteration 170 it falls off, and what it lands on has **no
X-points at all**. Machine A, run identically, reports `Zaxis` span **0.000000**
over four passes and converges: it is neither elongated nor up–down symmetric.

**SO THE FAILURE IS THE EQUILIBRIUM'S AND NOT THE ITERATION'S**, and the
inverse solve suppresses it because its geometric targets pin the shape. The
minimum a real machine has that this problem does not is ONE antisymmetric
current combination driven by the axis's own height — no flux targeted
anywhere, nothing re-optimised, the equilibrium still the forward one for
whatever currents result. Gain in amps per metre of axis error, 220 passes
allowed, DIII-D:

| gain | status | its | `Zaxis` span | `psi_axis` | vs the inverse reference | max ∣d*I*∣ |
|---|---|---|---|---|---|---|
| 0 | maxits | 221 | 1.69e+00 | −1.164772e-01 | 1.310e+00 | 0 |
| **+1e5** | **rtol** | **12** | **5.30e-09** | **3.758545e-01** | **9.980e-08** | **2.9434e-04 A** |
| −1e5 | maxits | 221 | 1.64e+00 | −1.819602e-01 | 1.484e+00 | 4.0446e+04 |
| **+1e6** | **rtol** | **12** | **1.02e-09** | **3.758545e-01** | **9.753e-08** | **8.0980e-04 A** |
| −1e6 | maxits | 221 | 1.80e+00 | 8.967487e+02 | 2.385e+03 | 2.1287e+08 |
| +1e7 | rtol | 31 | 1.22e+00 | −1.082351e+03 | 2.881e+03 | 2.6687e+08 |
| −1e7 | rtol | 65 | 1.32e+00 | 2.683694e+03 | 7.139e+03 | 6.4519e+08 |

**EIGHT HUNDRED MICROAMPS AGAINST COIL CURRENTS OF 1.75e+05 A**, and the
forward answer then reproduces the inverse reference to **1e-07**. Read the
`max |dI|` column before the agreement: it is what says the pinned equilibrium
IS the unpinned one rather than a differently-controlled machine. The sign is
load bearing and the wrong one drives the instability rather than damping it;
the working window is about a decade wide and 1e+07 overshoots into a different
failure.

**THE WHOLE SET, FORWARD, PINNED AT EACH REFERENCE'S OWN AXIS HEIGHT**, at a
fixed 1e+06 A/m — which is what exposes that a fixed gain does not travel:

| case | status | its | L | `psi_axis` | vs inverse | max ∣d*I*∣ |
|---|---|---|---|---|---|---|
| A testtokamak classic | rtol 4 | | 1.000000 | 8.2717509271e-02 | **6.260e-08** | 1.18e-05 |
| A shaped | rtol 4 | | 1.000000 | 8.2716844494e-02 | **6.392e-08** | 1.66e-05 |
| B peaked ff′ | rtol 34 | | 1.000000 | 9.7237762126e-02 | **1.948e-07** | 5.44e-03 |
| B shaped | rtol 38 | | 1.000000 | 9.7234142830e-02 | **2.015e-07** | 6.46e-03 |
| C MAST | rtol 168 | | 0.421707 | −1.7087234940e+02 | 2.949e+03 | 7.92e+07 |
| D TCV | maxits | 401 | 165.6 | 2.5053895020e+00 | 1.102e+02 | 9.20e+05 |
| D shaped | rtol 52 | | 0.146006 | −9.5142389060e+00 | 4.387e+02 | 1.64e+07 |
| E diamagnetic | rtol 17 | | 1.000000 | 8.3765613708e-02 | **9.409e-08** | 3.12e-05 |
| E shaped | rtol 4 | | 1.000000 | 8.3767813752e-02 | **8.369e-08** | 3.19e-05 |
| F DIII-D | rtol 12 | | 1.000000 | 3.7585449384e-01 | **9.753e-08** | 8.10e-04 |
| F shaped | rtol 14 | | 1.000000 | 3.7586544435e-01 | **1.018e-07** | 5.36e-04 |
| G MAST-U | rtol 116 | | 0.318944 | −2.5506552801e+02 | 3.041e+03 | 7.27e+07 |

**EIGHT OF TWELVE CLOSE TO 1e-07 AT A CURRENT INCREMENT OF MICROAMPS TO
MILLIAMPS**, against a table where every unpinned row but A's and E's was
2e-01 to 1e+02 away. The four that do not are over-driven, not unpinnable:
`max |dI|` reads 1e+07 to 1e+08 A, so 1e+06 A/m is simply the wrong gain for a
conductor table that resolves solenoid windings into individual turns. The
harness therefore CALIBRATES by default — one probe increment and the axis's
response to it, then a Newton step on the axis height — and a fixed gain is the
override rather than the mode. **The instrument here is the axis position and
nothing else could have found this**: every one of these rows reports a
plausible residual history.

### M-105

**MACHINE A'S 2.6e-02 IS THE PRESCRIBED-CURRENT ROW AND NOTHING ELSE, AND THE
ROW IS NEARLY SINGULAR BECAUSE `I_p(λ)` HAS A MAXIMUM WHERE THE TARGET IS.**

**READ [M-107](#m-107) FIRST: THE 2.6e-02 IS CONDITIONAL ON A STARTING VALUE.**
Every run below has `[source] PsiAxis` = 8.271751e-02, the reference's own
`psi_axis`, which was set while investigating something else. That lies inside
the basin of the WORSE of two roots. From the shipped `PsiAxis` = 6.023264e-02
the same machine and the byte-identical seed give `psi_ax` 8.248908e-02 at
λ = 0.995386 — **2.76e-03**. What follows is a correct and complete account of
the branch reached from the reference's own value, and the fold it measures is
real; it is not an account of what the shipped configuration does.
Nine causes were eliminated one at a time — MEQ's h and p refinement,
freegs4e's grid in both modes, the conductor model, the exterior mode count,
the exterior RADIUS, the profile conversion, the plasma support, the X-point
row, and the source's own arithmetic. Dropping `[source] PlasmaCurrent`, which
fixes the profile amplitude at one and removes that row from the bordered
Newton, moves the disagreement by a factor of **seventeen**:

| | `psi_ax` | `psi_bnd` | `Raxis` | X-point | `I_p` | vs the reference |
|---|---|---|---|---|---|---|
| freegs4e, 129², both modes | 8.271751e-02 | 3.240413e-02 | 1.351273 | ( 1.093144, −0.603965 ) | 2.0000000e+05 | — |
| MEQ, λ = 1, k=2, 2039 el | 8.2514479e-02 | 3.2414841e-02 | 1.348055 | ( 1.092655, −0.604123 ) | 1.997272e+05 | **2.46e-03** |
| MEQ, λ = 1, k=3, 8267 el | 8.2594336e-02 | 3.2391442e-02 | 1.349135 | ( 1.092787, −0.603935 ) | 1.998794e+05 | **1.49e-03** |
| MEQ, `I_p` = 2.0e5, k=2 | 8.068175e-02 | 3.137674e-02 | 1.326374 | ( 1.082839, −0.597793 ) | 2.000007e+05 | 2.46e-02 |
| MEQ, `I_p` = 2.0e5, k=3 | 8.056770e-02 | 3.133561e-02 | 1.325200 | ( 1.082380, −0.597556 ) | 2.000000e+05 | 2.60e-02 |

**AND IT REFINES WHERE THE CONSTRAINED ONE DOES NOT.** At fixed amplitude the
X-point agrees to **3.6e-04 m**, `psi_bnd` to 3.9e-04, and every column improves
from `k=2` to `k=3`; with the current prescribed the gap is 2.5e-02 at BOTH
resolutions and at every exterior radius.

**THE CURVE THE ROW SOLVES ON**, measured by scaling the profile tables — which
is exactly what the row's unknown does — and integrating each converged field's
own source over its own support. The reconstruction reproduces MEQ's reported
`I_p` to **2.3e-06** where MEQ prescribes one, which is what makes this a
measurement of MEQ rather than of a grid:

| λ | `psi_ax` | `psi_bnd` | core area | `Raxis` | **`I_p`** |
|---|---|---|---|---|---|
| 0.90 | 7.3606623e-02 | 2.9000581e-02 | 0.728590 | 1.254147 | 1.922488e+05 |
| 0.94 | 7.7468288e-02 | 3.0334437e-02 | 0.756819 | 1.290669 | 1.968365e+05 |
| **0.967085** | 8.0788118e-02 | 3.1398313e-02 | 0.784376 | 1.327892 | **2.001277e+05** |
| 0.97 | — | — | — | — | **exit 2, no convergence** |
| 1.00 | 8.2514479e-02 | 3.2414841e-02 | 0.764213 | 1.348055 | 1.997272e+05 |
| 1.03 | 8.2581941e-02 | 3.3067558e-02 | 0.718508 | 1.341968 | 1.981757e+05 |
| 1.06 | 8.2889176e-02 | 3.3622904e-02 | 0.689606 | 1.338809 | 1.972577e+05 |
| 1.10 | 8.3444675e-02 | 3.4379539e-02 | 0.655327 | 1.335586 | 1.964680e+05 |

**`I_p(λ)` PEAKS AT 2.0013e+05 AND THE PRESCRIBED CURRENT IS 2.0000000e+05 —
0.065% BELOW THE MAXIMUM.** So the row `I_p(λ) = I_target` has **two** roots,
one either side of the peak, and the amplitude is free to move three per cent
between them.

**THIS SAID "A NEAR TANGENCY, `∂I_p/∂λ ≈ 0` BETWEEN THEM" AND THAT IS
FALSIFIED** — see [M-108](#m-108), where the derivative is computed from the
bordered solve's own Schur complement rather than inferred from the sweep's
shape. The two roots carry **+0.664** and **−0.296** of `I_p` per unit relative
scale: opposite signs, neither near zero, each individually well conditioned.
The defect is NON-UNIQUENESS and not singularity. Everything below about the
curve and its maximum stands; what was wrong was reading a close pair of roots
as a merged one. That is the physical
equilibrium current limit and not a numerical artefact — raising the amplitude
raises the current density and shrinks the plasma, and past the peak the second
wins, core area falling 0.784 → 0.655 while `I_p` falls with it. **The one
amplitude that fails to converge outright is 0.97, at the top of the curve.**

**AND THE REFERENCE'S OWN OPERATING POINT IS NOT A ROOT AT ALL.** freegs4e
reports `Ip_logic` L = 1.000000, i.e. λ = 1 with `I_p` = 2.0000000e+05 exactly;
MEQ at λ = 1 carries **1.997272e+05**, 1.36e-03 short — so 2.0e5 is
unreachable there and the constrained Newton leaves for the tangency. That
deficit is MEQ's own discretisation of the current integral and it HALVES:

| | `I_p` at λ = 1 | short of the target by |
|---|---|---|
| k = 2, 2039 elements | 1.997272e+05 | 1.364e-03 |
| k = 3, 8267 elements | 1.998794e+05 | 6.03e-04 |
| freegs4e | 2.0000000e+05 | — |

**TWO PREDICTIONS FOLLOW AND BOTH HOLD.** Ask for a current safely off the
tangency and MEQ returns to the reference's branch; ask for one nearer it and
it does not:

| prescribed `I_p` | profile scale | `psi_ax` | vs the reference |
|---|---|---|---|
| 2.0000e+05 | 0.967085 | 8.068175e-02 | 2.46e-02 |
| 1.9950e+05 | 0.962807 | 8.010147e-02 | 3.16e-02 |
| **1.9900e+05** | **1.013095** | **8.251901e-02** | **2.42e-03** |
| **1.9800e+05** | **1.034809** | **8.262044e-02** | **1.17e-03** |

**A HALF PER CENT LESS CURRENT IS A FACTOR OF TEN CLOSER, AND ONE PER CENT LESS
IS A FACTOR OF TWENTY-ONE.** The scale flips from 0.963 to 1.013 between those
two rows, which is the solver crossing the peak — and at 1.99e+05 and below
there is only one root on the falling branch, `∂I_p/∂λ` = −5.2e+04 rather than
zero, and the row is ordinary.

**THE PROFILE SCALE IS THEREFORE A DIAGNOSTIC MEQ ALREADY PRINTS AND NOBODY
READS.** λ far from one on a case whose profiles came from a converged
reference means the current row is working hard, and here it means it is at a
tangency. The control is DIII-D, seeded identically: **λ = 0.9993353**, and its
constrained answer is **3.48e-04** from its reference — against its own
fixed-amplitude arm at 1.88e-03. M-103's "DIII-D is the only machine that
closes" is that number.

### M-106

**THE WHOLE DIVERTED SET WITH THE SEED FIXED, BOTH ARMS OF THE CURRENT ROW —
AND THE PROFILE SCALE IS A CLEAN SCREEN ON THE ANSWER.** M-103 runs the shipped
configurations, which start COLD, so every failure in it is a failure to find
the basin and the table cannot separate that from anything else. Here each
machine is seeded from its OWN reference reconstructed by Green's functions over
its own half-disc at `n = 128` — `mkexactguess.py`, the same construction
machine A's `a128.gf` uses — and then run TWICE from that one seed: with
`[source] PlasmaCurrent`, and without it, which fixes the profile amplitude at
one. Everything else is the shipped configuration.

| machine | reference `psi_axis` | prescribed | λ | rel | fixed, λ = 1 | rel |
|---|---|---|---|---|---|---|
| A testtokamak classic | 8.2717514448e-02 | 8.248908e-02 | 0.995386 | **2.762e-03** | no damped step | |
| A shaped | 8.2716849781e-02 | 8.246192e-02 | 0.994454 | **3.082e-03** | no damped step | |
| B peaked ff′ | 9.7237743182e-02 | 1.640031e-01 | **1.230204** | 6.866e-01 **refused** | no damped step | |
| B shaped | 9.7234123240e-02 | 1.640122e-01 | **1.230356** | 6.868e-01 **refused** | no damped step | |
| C MAST | 5.7958202646e-02 | no damped step | | | no damped step | |
| C shaped | 9.0113316844e-02 | no damped step | | | 9.018843e-02 | **8.335e-04** |
| D TCV | 2.2533925689e-02 | no damped step | | | no damped step | |
| D shaped | 2.1736011974e-02 | 2.250051e-02 | 0.993500 | 3.517e-02 | 2.280990e-02 | 4.941e-02 |
| E diamagnetic | 8.3765621589e-02 | 1.025203e-01 | **0.414259** | 2.239e-01 | no damped step | |
| E shaped | 8.3767820763e-02 | 1.025271e-01 | **0.414124** | 2.239e-01 | no damped step | |
| F DIII-D | 3.7585453049e-01 | 3.759851e-01 | 0.999335 | **3.474e-04** | 3.765616e-01 | 1.881e-03 |
| F shaped | 3.7586548261e-01 | 3.758642e-01 | 0.998925 | **3.412e-06** | 3.767981e-01 | 2.481e-03 |
| G MAST-U | 8.3895051962e-02 | 8.395672e-02 | 0.999756 | **7.351e-04** | 8.399722e-02 | 1.218e-03 |
| G shaped | 9.2555437956e-02 | no damped step | | | no damped step | |

**SIX OF FOURTEEN CLOSE AT 3.5e-02 OR BETTER ON THE PRESCRIBED ARM, AGAINST TWO
IN M-103** — and every one of M-103's *seven* `exit 2` rows that now converges
does so because of the seed, not because anything in the solver changed. **The
cold start is the larger half of M-103.**

**AND λ SEPARATES THE ANSWERS PERFECTLY.** Every row with ∣λ − 1∣ < 0.007 lands
between **3.4e-06 and 3.5e-02**; every row outside it — both B's at 1.2302 and
both E's at 0.4142 — is **2.2e-01 to 6.9e-01** out. Nothing in between. On a
case whose profiles came from a converged reference the amplitude has no reason
to move, so a scale far from one says the current row went somewhere else, and
it says so BEFORE any reference is consulted. The axis guard independently
refuses both B rows, so the two disagreeing screens agree.

**DROPPING THE CURRENT ROW IS NOT A GENERAL CURE AND MACHINE A'S FOLD IS A'S.**
It is the whole of the answer on **C shaped**, where the prescribed arm has no
finite damped step and the fixed arm reaches **8.335e-04**; it is the opposite on
A, F and G, where the fixed arm is worse or fails outright — F shaped goes from
3.4e-06 to 2.5e-03 and both A rows stop converging. M-105's tangency is a
property of machine A at a particular starting point, not a mode of the set.

**THE SEED RESOLUTION IS A VARIABLE AND IT IS NOT SMALL.** `mkexactguess.py`
defaults to `n = 32`, measured on `examples/limited-tokamak.toml`; it does not
carry here. At `n = 32` machine A fails on BOTH arms, and at `n = 128` — the
identical construction — both converge. A sweep built on the default would have
measured the seed.

### M-107

**`[source] PsiAxis` IS DOCUMENTED AS A STARTING VALUE AND IT SELECTS WHICH
EQUILIBRIUM IS REPORTED — AND THE WRONG ROOT'S BASIN IS A WINDOW AROUND THE
RIGHT ANSWER.** Machine A, one seed field (byte-identical), one configuration,
one key changed:

| `PsiAxis` | `psi_ax` | λ | vs the reference 8.271751e-02 |
|---|---|---|---|
| 3.0e-02 | exit 2 | | |
| 4.0e-02 | exit 2 | | |
| 4.5e-02 | 8.248908e-02 | 0.995386 | **2.762e-03** |
| 5.0e-02 | 8.248908e-02 | 0.995386 | **2.762e-03** |
| **6.023264e-02** — what `make_diverted_case.py` writes | 8.248908e-02 | 0.995386 | **2.762e-03** |
| 7.0e-02 | 8.248908e-02 | 0.995386 | **2.762e-03** |
| 7.2e-02 | 8.248908e-02 | 0.995386 | **2.762e-03** |
| 7.4e-02 | 8.068175e-02 | 0.967085 | 2.463e-02 |
| 7.6e-02 | 8.068175e-02 | 0.967085 | 2.463e-02 |
| 8.0e-02 | 8.068175e-02 | 0.967085 | 2.463e-02 |
| **8.271751e-02** — **the reference's own answer** | 8.068175e-02 | 0.967085 | **2.463e-02** |
| 8.4e-02 | 8.068175e-02 | 0.967085 | 2.463e-02 |
| 8.6e-02 | 8.248908e-02 | 0.995386 | **2.762e-03** |
| 9.0e-02 | 8.248908e-02 | 0.995386 | **2.762e-03** |
| 1.0e-01 | 8.248908e-02 | 0.995386 | **2.762e-03** |

**TWO ROOTS, EACH REPRODUCED TO SEVEN FIGURES FROM EVERY START IN ITS BASIN**,
so these are genuine separated solutions of the constrained problem and not a
convergence artefact. The good one carries λ = 0.995386 and sits 2.8e-03 from
the reference; the bad one carries λ = 0.967085 and sits 2.5e-02 away. **The bad
root's basin is the interval ( 7.2e-02, 8.6e-02 ), and `psi_axis` = 8.271751e-02
— the value a reader would call the right answer — is inside it.** Starting 25
per cent LOW or 20 per cent HIGH gives the better equilibrium; starting at the
truth gives the worse one by a factor of nine.

**THIS IS THE FOLD OF M-105 SEEN FROM THE OTHER SIDE.** `I_p( λ )` peaks 0.065
per cent below the prescribed current on the branch reached from the reference's
own `psi_axis`, which is why that branch cannot satisfy the row; the branch
reached from further away can, at λ = 0.995386. Both facts are one fact.

**AND IT IS A TRAP FOR EXACTLY THE PERSON WHO KNOWS THE ANSWER.** M-105's
headline — machine A disagreeing by 2.46e-02 — was taken on a configuration
whose `PsiAxis` had been set to the reference's own value while investigating
something else. The shipped configurations write the COLD GUESS's peak and are
outside the bad basin, so **M-96, M-103 and M-87 are not affected**; what was
affected was the diagnosis. *A key documented as "only a starting value" is
still an input when the problem has more than one solution.*

### M-108

**`d I_p / d lambda` ALONG THE MANIFOLD, WHICH THE BORDERED SOLVE ALREADY HAS —
AND IT FALSIFIES M-105'S "TANGENCY" READING.** The dense border matrix the
solve forms and factors every iteration IS the Schur complement `S` of the
border, so the sensitivity of the current residual to the profile scale with
every other border residual held at zero is `1/( S^-1 )_II`: one solve against
a factorisation that already exists, on a matrix of order at most `4 + modes`.
`GradShafranovSolver::plasmaCurrentSensitivity()` reports it and the driver
prints it beside the profile scale.

**IT AGREES WITH AN INDEPENDENT FINITE DIFFERENCE**, which is the only check
that means anything for a derivative. Machine A, both roots, against the
amplitude sweep of M-105 read as a secant:

| root | scale | Schur complement | M-105's secant | |
|---|---|---|---|---|
| from `PsiAxis` = 8.271751e-02 | 0.967085 | **+1.3729e+05 A** | +1.2151e+05 A ( λ 0.94 → 0.967 ) | rising branch |
| from `PsiAxis` = 6.023264e-02 | 0.995386 | **−5.9473e+04 A** | −5.172e+04 A ( λ 1.00 → 1.03 ) | falling branch |

**AND THE TWO ROOTS ARE NOT A DOUBLE ROOT.** M-105 read the amplitude sweep as
a near tangency — the target 0.065 per cent below the maximum, `∂I_p/∂λ ≈ 0`,
two roots nearly on top of each other. The derivative says otherwise: the two
solutions carry **+0.664** and **−0.296** of `I_p` per unit relative scale,
opposite signs and neither near zero. `I_p( λ )` turns over BETWEEN them and each
root is individually well conditioned. **The defect is non-uniqueness, not
singularity** — the prescribed current is attainable twice, and `[source]
PsiAxis` picks which. M-107's basin map is the same statement.

**SO THE USEFUL WARNING IS THE SIGN AND NOT THE MAGNITUDE.** A negative
elasticity means more profile would carry LESS current — the plasma shrinking
faster than its current density rises — so `I_p( λ )` has already turned over
below this scale and a second equilibrium carries the same current at a lower
one. The driver says exactly that, and keeps a separate near-zero warning for a
genuine tangency, which this case turns out not to be. On machine A the warning
fires on the root at 0.995386 and correctly asserts the existence of the one at
0.967085.

**AND IT DOES NOT FALSE-ALARM ON THE TWO MACHINES THAT AGREE BEST.** The same
runs from M-106, seeded identically:

| machine | scale | `d I_p/d scale` | elasticity | warning | vs its reference |
|---|---|---|---|---|---|
| F DIII-D | 0.999335 | +1.3067e+06 A | **+1.306** | none | 3.47e-04 |
| G MAST-U | 0.999756 | +4.6585e+05 A | **+0.776** | none | 7.35e-04 |
| A, from `PsiAxis` 8.2718e-02 | 0.967085 | +1.3729e+05 A | +0.664 | none | 2.46e-02 |
| A, from `PsiAxis` 6.0233e-02 | 0.995386 | −5.9473e+04 A | **−0.296** | **fires** | 2.76e-03 |

**An elasticity of about one is what an ordinary current row looks like** — the
current very nearly proportional to the amplitude — and the two machines that
reproduce their references to 1e-04 read 1.31 and 0.78. **The warning is a
non-uniqueness warning and not an accuracy one**, which is why it fires on the
BETTER of machine A's two roots: that root is real, is the one nearer the
reference, and does have a second equilibrium at a lower scale carrying the same
current. Both facts are true at once and the diagnostic says the one it can
know.

### M-109

**THE AXIS GUARD WAS REFUSING A CORRECT ANSWER ON A READING IT PRINTED THE
REFUTATION OF, AND MACHINE C IS THE CASE.** M-103 has MAST under filament
conductors converging in 3 Newton steps and then refused, exit 1, *"THE PLASMA
CONTAINS THE SYMMETRY AXIS"*. The driver tested
`checkAxisSource().axisInsidePlasma` — the pure level set `-psi_bnd/span > 0`,
which `insidePlasma()` asks POINTWISE — while printing *"plasma component kept
it: NO"* on the line above it.

**THE LEVEL SET AND THE SUPPORT GENUINELY DISAGREE THERE, AND THE REFERENCE
SETTLES WHICH IS RIGHT.** `C_mast_spherical.npz` converges with
`psi_bndry` = **−2.3596e-02**, so `Psi( r = 0 ) = +0.2863` and the level set does
contain the axis. Its near-axis lobe is nevertheless a SEPARATE component:
walking the reference's own midplane inboard from the magnetic axis, `psi` falls
to **−3.0876e-02 at R = 0.1594**, below `psi_bnd`, before rising to −1.3990e-02
at `R = 0.100`. A flood fill run on `{ Psi_N > 0 }` says otherwise only because
it leaks through the X-point — the failure this tree already records for its own
fill. Erode the threshold and it closes:

| ε | components | the axis component reaches |
|---|---|---|
| 0 | 1 | R = 0.1000 |
| 1e-04 | 2 | **R = 0.3078** |
| 1e-02 | 2 | R = 0.3227 |

against the reference's own `core_mask`, which reaches **R = 0.308**.
`G_mastu_simple` behaves identically: 0.1000 at `ε = 0`, **0.2781** at 1e-04,
`core_mask` 0.278.

**AND freegs4e IS NOT USING A WALL THERE AT ALL.** `C_mast_spherical.npz`
carries **zero** wall points, so `mask_inside_limiter` is `None` and nothing
masks `Jtor` geometrically. Connectivity is the whole of what keeps that
reference off the axis, which is what says MEQ's fill was right and the guard
was not.

**BOTH GUARDS MOVED ONTO THE SUPPORT, AND THE SECOND ONE MATTERS AS MUCH AS THE
FIRST.** `bounded` reads `worstOnAxis`, which is `f()` evaluated pointwise and
knows nothing about elements, so correcting only the topology guard would have
moved the same false refusal one branch down. `AxisSourceCheck` now carries
`supportReachesAxis`, `worstOnAxisInSupport` and `boundedInSupport` beside the
level-set readings, which are kept and REPORTED — a run where the two disagree
is one where the fill is the only thing between the load and a `1/r` pole.

| | before | after |
|---|---|---|
| machine C, filament, shipped configuration | **exit 1, refused** | **exit 0** |
| `psi_ax` | 5.896474e-02 | 5.896474e-02, unchanged |
| `psi_bnd` | −2.356703e-02 | −2.356703e-02, unchanged |
| axis | ( 1.0156, −0.0000 ) | ( 1.0156, −0.0000 ) |
| against its reference 5.795820e-02 | — | **1.737e-02** |

**IT IS THE SAME SOLVE TO EVERY DIGIT** — only the verdict changed — and it now
emits a warning naming the disagreement: `| F |` on `r = 0` is assembled as zero
where the pointwise test reads 9.208398e-02.

**A THREADING NOTE, BECAUSE IT COST A WRONG CONCLUSION FOR TEN MINUTES.** That
case is on the knife edge CLAUDE.md documents: at `MKL_NUM_THREADS=1` it fails
in the bordered Newton with no finite damped step, and at `OMP = MKL = 16` it
converges in 3. M-103 was taken at 16. A guard fix verified at 1 would have read
as no change at all.

### M-110

**`[source] ExcludeAttributes` — A REGION THAT CAN NEVER BE PLASMA, BY MESH
ATTRIBUTE, AND IT IS INERT WHERE IT SHOULD BE.** `psi_bnd` confines the plasma
and `PlasmaConnectivity` separates the lobes a level set leaves joined; what
neither can do is know that a lobe is behind a wall, which matters where several
O-points sit across a saddle. An ATTRIBUTE rather than a polygon because the
support is re-decided on every residual evaluation and a vessel does not move.

**THE WHOLE PATH, END TO END**: `[mesh.generate] Vessel` → `meq --mesh-command`
→ `halfdisc.py --vessel` → element attribute **30** outside the polygon →
`[source] ExcludeAttributes = [ 30 ]` → the fill never traverses, seeds at or
counts those elements, and meq::PlasmaComponent::holds() is false there ahead of
its own constant-true shortcut, so the source is off with or without
`ConfineToPlasma`.

| test | result |
|---|---|
| standalone mesh, 2 coils + a hexagonal vessel | attributes **1: 71, 10: 4, 11: 4, 30: 538** |
| machine A, coil interiors excluded ( 10–13 ) | **104 of 2039**, converges, `psi_ax` unchanged to 7 digits |
| machine A, real vessel, `[mesh.generate] Vessel` | **1373 of 2105**, converges, `psi_ax` 8.255892e-02 |
| the same mesh and start, exclusion REMOVED | `psi_ax` 8.255892e-02, λ 0.9938978 — **identical to every printed digit** |
| an attribute the mesh does not carry | **refused**, exit 1, listing the attributes present |

**INERT IS THE REQUIRED PROPERTY AND IT IS WHAT WAS MEASURED**: on a case where
the fill was already right, excluding 1373 of 2105 elements changes nothing at
all. It buys work, not answers, until the configuration is one connectivity
cannot settle.

**AND FRAGMENTING A VESSEL PERTURBS THE MESH, WHICH ON MACHINE A COSTS
CONVERGENCE FROM ONE PARTICULAR START.** With `PsiAxis` at the shipped
6.023264e-02 the vessel-fragmented mesh fails in the bordered Newton — **and so
does the same mesh with the exclusion REMOVED**, which is what says this is the
mesh and not the feature. Regenerating the same mesh WITHOUT the vessel
converges. Moving `PsiAxis` to 5.0e-02, which M-107 shows is the same basin,
converges on the vessel mesh at **1.917e-03** against the committed mesh's
2.762e-03. Machine A sits near a fold with two roots; a three per cent
perturbation of the mesh is exactly what such a case is fragile to.

### M-111

**THE RACE ON THREE MACHINES, AND DIII-D'S 1.34× IS THE BEST CASE RATHER THAN
THE TYPICAL ONE.** M-97 and M-102 race one machine. The cold-forward audit
behind [M-112](#m-112) says which others can be raced at all — MEQ must close
from a cold start and freegs4e must reach the same equilibrium — and that is
**C MAST** and **G MAST-U** beside **F DIII-D**. Interleaved, best of five, one
case per invocation, each gated on the one-minute load falling below 0.5:

| case | freegs4e 129², cold | MEQ, cold | ratio | MEQ's rung |
|---|---|---|---|---|
| **F DIII-D** | **3.903** | **5.408** | **1.39×** | `k = 2`, 4848 el |
| **G MAST-U** | **8.612** | **26.298** | **3.05×** | `k = 3`, 6707 el |
| **C MAST** | **4.666** | **18.044** | **3.87×** | `k = 2`, 4198 el |

Spreads 1.03× on freegs4e's arm and 1.04× on MEQ's, and **`psi_ax` is identical
across all five rounds of every case** — 5.896474e-02, 3.759851e-01,
8.403188e-02 — which is what makes these timings rather than solves.

**M-102 REPRODUCES, WHICH IS THE CHECK THAT LICENSES THE OTHER TWO ROWS.**
freegs4e reads 3.903 against M-102's 3.922, **1.005×**; MEQ reads 5.408 against
5.245, 1.03×; the ratio 1.39× against 1.34×. Both inside M-102's own spreads, on
a day's worth of intervening commits.

**AND THE ACCURACY THOSE SECONDS BUY**, from the same rungs, each measured
against its own code's finest run:

| case | rel `L2` | ... off the coils | `psi_ax` |
|---|---|---|---|
| F DIII-D | 5.785e-03 | **7.599e-04** | 3.54e-04 |
| C MAST | 5.342e-03 | **3.268e-03** | 1.76e-02 |
| G MAST-U | 1.919e-02 | **1.542e-02** | 1.29e-04 |

**`k = 1` FAILS ON ALL THREE CASES, BOTH RUNGS**, which is uniform and is new.
**C converges on ONE rung of eight** — `k2r0` — and every other one fails,
so its row is a single point and not a ladder. **F's error does not refine**:
5.785e-03, 5.802e-03, 5.763e-03, 5.803e-03 across a 16× range in dofs, with the
off-coil column flat at 6.3e-04 to 7.6e-04. That is the conductor model, MEQ's
rectangles against freegs4e's filaments, and M-87 sizes it; it is not a
discretisation error and no rung buys it down. **Ignore G's `psi_bnd`**: the
reference's is 5.6e-04, so a relative error on it is meaningless.

**A DEFECT IN `race.py` ITSELF, AND IT IS WHY THIS IS NOT ITS OUTPUT.** The
script runs its cases back to back and its rungs back to back, so on a
multi-case invocation every case is timed on a machine the previous ones have
been holding at sixteen threads. Measured on DIII-D's `k2r0`, one binary, one
configuration, one answer:

| | wall | solve |
|---|---|---|
| inside a three-case `race.py C F G` | 9.49 | 8.01 |
| **settled, interleaved, best of five** | **5.41** | **≈3.8** |

**A factor of two of pure machine state**, with `psi_ax` 3.759851e-01 either
way. M-97 and M-102 never met this because they are single-case and interleaved.
Anything multi-case must settle between cases and interleave its control, or it
is measuring the box.

**AND THE CONTROL HAS TO BE INTERLEAVED WITH WHAT IT CONTROLS, WHICH IS THE
PART THAT WENT WRONG.** freegs4e's arm was read from inside the race and MEQ's
from after it, and the pair was compared: freegs4e moved 1.05% and MEQ 2.2×, so
the machine looked innocent and MEQ looked like a regression. It was neither.
M-102 states the interleaving as protocol and the reason is exactly this.

### M-112

**NEITHER CODE'S COLD FORWARD RUN WORKS ANYWHERE, AND freegs4e'S FAILS ALL
FOURTEEN.** A COLD FORWARD run is the machine's currents frozen, the profiles,
`I_p` prescribed, no constraint, and nothing that was solved for. It is exactly
what MEQ's shipped configurations are asked -- `make_diverted_case.py` builds
their `[initialguess]` with `mkcoldguess.py`, from the conductors, the target
current and the shape the operator asked for.

`forward.py --cold-forward`, every reference:

| | rel | | rel |
|---|---|---|---|
| A classic | 2.623e-01 | A shaped | 2.622e-01 |
| B peaked ff′ | 2.065e-01 | B shaped | 2.063e-01 |
| C MAST | 3.709e+00 | C shaped | 2.067e+00 |
| D TCV | 1.093e+02 | D shaped | 1.300e+00 |
| E diamagnetic | 4.442e-02 | E shaped | 1.480e+00 |
| **F DIII-D** | **1.310e+00** | F shaped | 1.272e+00 |
| G MAST-U | 6.797e-01 | G shaped | 5.580e+01 |

**Nine of the fourteen run to the 401-iteration cap**; the five reporting `rtol`
converge to a different equilibrium, with `Ip_logic` L reaching 6.9, 27.8 and
0.093. **DIII-D — the one machine that closes every other way — is 1.310e+00 out
with `psi_axis` NEGATIVE.**

**MEQ's cold table, taken on the same tree**, `psi_ax` against each reference:

| | rel | | rel |
|---|---|---|---|
| A classic | 5.687e-01 | A shaped | 5.690e-01 |
| B ff′ | exit 2 | B shaped | exit 2 |
| **C MAST** | **1.737e-02** | C shaped | exit 2 |
| D TCV | exit 2 | **D shaped** | **3.517e-02** |
| E diamagnetic | exit 2 | E shaped | exit 2 |
| **F DIII-D** | **3.474e-04** | **F shaped** | **3.412e-06** |
| **G MAST-U** | **4.943e-02** | G shaped | 4.171e+00 |

**Five of fourteen, against freegs4e's zero, so the intersection is EMPTY.**
There is no case in this benchmark where both codes solve the same problem from
a cold start. C is a pass only since M-109 moved the axis guard onto the
support; M-103 refused it.

**SO THE RACE HAS NEVER BEEN FORWARD AGAINST FORWARD.** `race.py` drives
freegs4e through `fgsref.py`, which builds with `SyncConstrain` -- the INVERSE
solve, whose geometric targets are what carry it from cold. That is why its arm
converges on all seven there while its forward arm fails everything, and it is
the mode mismatch `CLAUDE_FB.md`'s standing rule names, sitting inside M-97,
M-102 and M-111. The wall clocks compare different amounts of work and the
accuracy columns compare MEQ-forward against freegs4e-inverse.

**AND `make_diverted_case.py` WROTE A COMMENT NAMING THE WRONG GUESS BUILDER**,
telling every generated file that `mkexactguess.py` had built its guess while
the code called `mkcoldguess.py`. The two make this a different benchmark:
`mkexactguess` sums Green's functions over the reference's own converged `Jtor`
and puts the answer in the starting position. Believing the comment is what made
M-103's failures look like something other than cold-start failures.

### M-113

**`[solver] XPointMeritWeight`: IT BUYS ITERATIONS, IT NEVER MOVES THE ANSWER,
AND ITS SENSITIVITY INVERTS BETWEEN CASES.** M-111 finds MAST spending 26 of its
43 bootstrap iterations on a plateau where the line search halves eight or nine
times per step and accepts 1/128 to 1/256 of the direction. The decomposition
says why: the axis row carries about **0.75** of `augmentedNorm` and the X-point
rows about **0.07**, so a step that would fix the null is never worth taking and
it crawls a millimetre at a time -- 9 cm away from a seed that was 1 cm from the
answer, and back. The key multiplies the length `r h` that puts those two rows
into the merit, and into nothing else: the border still solves `q_r = q_z = 0`.

**MAST, `k = 2`, 4198 elements, every sweep counted:**

| weight | 0.1 | 0.3 | **1.0** | 3.0 | 6.0 | 10.0 | 15.0 | **20.0** | 30.0 | 100.0 |
|---|---|---|---|---|---|---|---|---|---|---|
| iterations | 57 | 57 | **56** | 53 | 44 | 39 | 35 | **35** | **fails** | **fails** |
| wall / s | 20.44 | 20.30 | 19.61 | 18.56 | 14.10 | 13.09 | 12.01 | **10.91** | exit 2 | exit 2 |

**`psi_ax` is 5.896474e-02 and the X-point ( 0.708468, −1.095754 ) at EVERY
weight that converges**, so this is work and not an answer. MAST goes from 3.87×
freegs4e to **2.34×**. Below one it does nothing, which fits the diagnosis: the
rows already counted for almost nothing, so starving them further changes
nothing.

**AND IT IS INERT WHERE THERE IS NO PLATEAU**, which is the other half of the
claim:

| | weight 1 | weight 10 |
|---|---|---|
| F DIII-D, `k2r0` | 12 its, 4.925 s | 12 its, 4.949 s |
| G MAST-U, `k3r0` | 19 its, 26.886 s | 18 its, 26.619 s |

So G's 3.05× really is dofs — 201,210 at `k = 3` against DIII-D's 87,264 — and
not a line search fighting itself.

**THE SENSITIVITY INVERTS, AND THAT IS WHY THE DEFAULT STAYS AT THE NATURAL
SCALE.** On `examples/diverted-tokamak.toml`, XP-3's own acceptance fixture:

| weight | 0.25 | **1.0** | 4.0 |
|---|---|---|---|
| bootstrap iterations | 13 | **14** | **82** |

`psi_ax` 8.2660036302e-02 and X ( 1.093103, −0.603529 ) at 0.25 and at 4.0,
identically -- the invariant holds -- but **raising the weight is 5.9× WORSE
there** where on MAST it is 1.6× better. So the natural scale `r h` is right on
one case and badly wrong on another, there is no value to promote to a default,
and **a weight must not be carried between machines.** Above about 30 MAST fails
outright with *"no damping gave a finite residual -- psi_ax through zero"*: the
other border rows starved in the merit, which is the predicted failure and is
loud rather than silent.

**`theXPointMeritWeightChangesTheWorkAndNotTheAnswer` is the regression**, and it
asserts the invariant at 0.25 AND at 4.0 -- both sides of one, since an
invariant tested only in the direction that pays is half tested.

### M-114

**THE MERIT'S THREE WEIGHTS, MEASURED — AND "THE AXIS ROW IS OVER-WEIGHTED" IS
FALSIFIED.** M-113 leaves the X-point weight helping MAST and nothing explaining
why. Two further weights were built to find out, and both change the MERIT and
nothing else: `gamma` and `xScale` appear only in `augmentedNorm()` and in the
printed record, never in the step, so **the Newton direction is bit-identical at
every weight and only the accept/reject verdict moves.** That is the fact the
rest of this rests on.

**FIRST, THE LINE SEARCH'S TRIAL LADDER, which is what says whether the
direction or the merit is at fault.** MAST's plateau, iteration 15, every trial
recorded:

```
 current merit   8.4313e-01
 full step       4.5122e+00 -> 3.638 -> 2.162 -> 1.244 -> 0.9328 -> 0.8609 -> 0.8455 -> 0.8427
```

The full step is **3.8 to 5.4× worse than the iterate** right across the plateau,
and the halving sequence descends monotonically, creeping under the current
merit only at the seventh halving. **The Newton step buys almost no first-order
decrease in this merit.** The discriminator proposed for this — *full step ≈
current means the direction is good* — is WRONG, and measured: the ratio says
"bad direction" while reweighting the same direction demonstrably halves the
work. What matters is the SHAPE of the halving sequence, not the full step's
ratio.

**SECOND, THE THREE WEIGHTS ON MAST**, `k = 2`, 4198 elements, all four support
sweeps counted, `psi_ax` **5.896474e-02** in every run that converged:

| | 0.01 | 0.03 | 0.1 | 0.3 | **1.0** | 3.0 | 10 | 20 | 30 |
|---|---|---|---|---|---|---|---|---|---|
| `XPointMeritWeight` | | | 57 | 57 | **56** | 53 | 39 | **35** | fails |
| `BorderMeritWeight` | | **38** | 52 | 55 | **56** | 56 | 56 | | 56 |
| `AxisMeritWeight` | **fails** | **210** | 106 | 79 | **56** | **fails** | | | |

**`BorderMeritWeight` SATURATES ABOVE ONE AND THAT IS A CORRECTNESS SIGNAL.**
3, 10 and 30 give bit-identical sweep counts, as they must: once the border
terms dominate, `gamma` is an overall factor and Armijo's test is scale
invariant. Raising the whole border is structurally a no-op; only lowering it —
which is raising the FIELD — does anything.

**AND THE AXIS ROW'S DEFAULT IS A LOCAL OPTIMUM.** Both directions degrade and
both ends fail: 0.03 costs 210 iterations against 56, and 0.01 and 3.0 do not
converge at all. `AxisMeritWeight = 0.03` WITH `XPointMeritWeight = 20` also
fails, so the two do not stack — they interfere.

**WHAT THAT KILLS.** The hypothesis was that `gamma cAx` carrying 0.75 of the
merit against the X-point rows' 0.07 made the axis row over-weighted, since the
two interventions that help both cut its share. It is wrong. Lowering the WHOLE
border by 33× helps and lowering the AXIS ALONE by 33× is four times worse, and
what separates them is that the first preserves the ratios among border rows and
the second destroys them. So the balance among the border rows is about right as
it stands, and the coincidence that both helpful knobs reduce the axis share was
a coincidence.

**WHAT SURVIVES.** `XPointMeritWeight` pays on MAST and is measured; nothing
explains WHY, and M-113's inversion on `examples/diverted-tokamak.toml` — 14
iterations to 82 at weight 4 — still stands unexplained beside it. `gamma` being
frozen at the first iterate cannot be the cause of any of this: it multiplies
every border row identically and CANCELS in any ratio between two of them, so it
can only move border against field. That was asserted here before it was checked
and it is wrong.

### M-115

**WHY MAST'S LINE SEARCH PLATEAUS: SIX HYPOTHESES MEASURED AND KILLED, AND WHAT
THE INSTRUMENTS SAY INSTEAD.** M-111 records MAST spending 26 of 43 bootstrap
iterations at damping 1/128 to 1/256 while the residual falls 0.7% a step, and
M-113 records a 20× weight on the X-point rows fixing it there and making
`examples/diverted-tokamak.toml` **5.9× worse**. This is the account of what
that is and is not.

**THE ONE STRUCTURAL FACT EVERYTHING RESTS ON**: `gamma` and `xScale` appear
only in `augmentedNorm()` and the printed record, never in the step. **The
Newton direction is bit-identical at every weight and only the accept/reject
verdict moves.**

| # | hypothesis | killed by |
|---|---|---|
| 1 | `gamma` frozen at the first iterate mis-balances the rows | it multiplies every border row identically and CANCELS in any ratio between two of them; it can only move border against field |
| 2 | the axis row is over-weighted ( 0.75 of the merit against the X-point's 0.07 ) | M-114: its default is a LOCAL OPTIMUM — 0.03 costs 210 iterations against 56, and 0.01 and 3.0 fail outright |
| 3 | the merit is mis-equilibrated against the Schur complement | `‖S_axis‖/‖S_xpoint‖` = 1.00/79.3 predicts a weight of **0.26** where the measured optimum is **20** — wrong by 77× and pointing the wrong way |
| 4 | the residual is discontinuous because the X-point is re-located per trial | the located point halves EXACTLY with the damping — 0.0972, 0.0487, 0.0243, 0.0121, 0.0061, 0.0030, 0.0016 — and extrapolates onto the iterate's own |
| 5 | the Jacobian is inexact, so the direction is not `-J^-1 R` | the FIELD block's deviation from `( 1 - a )‖R_0‖` falls by **3.3 to 4.0** per halving, which is the `O( a^2 )` an exact direction owes |
| 6 | so drop the constraints: judge Armijo on `‖R‖` alone | **MAST FAILS**, exit 2. The augmented residual wanders 6.19 → 108.8 → 155 → 6.01 → 47.3 → 8.42 → 24.1 and never closes |

**HYPOTHESIS 5 WAS KILLED BY A TEST PROPOSED AGAINST IT, AND THE ARITHMETIC
BEHIND IT WAS WRONG TWICE.** `augmentedNorm` returns `‖W F‖`, not `½‖W F‖²`, so
the exact-Newton prediction is `phi( t ) = ( 1 - t ) phi( 0 )` with slope
`-phi( 0 )` — not `-phi( 0 )²`, which is what was first compared against. And
the augmented merit cannot answer this question at all: it carries constraints
evaluated after `peakAt()`, `refreshLimiterContact()` and `refreshXPoint()` have
re-located the axis, the contact and the saddle inside every trial, and nothing
linearises those. **Only the field block is both smooth and linearised**, and it
is clean.

**WHAT NUMBER 6 ESTABLISHES IS WORTH MORE THAN THE FIX IT REFUTES**: the border
constraints are in the merit for a reason. Remove them and the line search takes
steps that improve the field while wrecking the constraints, and the solve never
closes. So the repair cannot be to exclude the un-linearised part; it has to be
to linearise it.

**AND THE SEVENTH IS SUPPORTED BUT NOT PROVEN. `psi_h` IS AN L2 FIELD AND
`psi_bnd - psi_h( x_X )` IS A POINT EVALUATION OF IT.** A broken space jumps
across faces, so that constraint can be discontinuous in the X-point's position
while the position itself is smooth. Measured on MAST at iteration 15, the
constraint and the element the point was found in:

```
 +3.736e-03/e1518  +2.609e-03/e1707  +2.138e-03/e3905  +1.994e-03/e2413
 +1.767e-03/e1132  +1.857e-03/e1132  +1.895e-03/e1132  +1.913e-03/e1132
```

Successive differences ×1e-6: **1127, 471, 144, 227**, 90, 38, 18 — ratios 2.39,
3.27, **0.63**, 2.52, 2.37, 2.11. **Within one element it halves cleanly; across
the face from e3905 into e2413 the ratio INVERTS.** The anomaly is about 5% of
the constraint, which is `gamma cBnd` ≈ 0.10 of a merit of 0.843, against an
Armijo decrease requirement of ~1e-07 at that damping — large enough to matter
by four orders. **This is a correlation with face crossings and not a
demonstration of cause**, and it is recorded as such: six hypotheses before it
also fitted the data they were built from.

**WHY MAST AND NOT DIII-D.** MAST's X-point traverses about 9 cm within a single
plateau iteration, crossing many elements from a seed 1 cm from the answer;
DIII-D's converges 1.06e-03 m from its seed and never leaves its element. The
X-point weight helps MAST because it changes how much a face crossing costs in
the merit, which is why it transfers to no other machine and inverts on the
fixture.

**THE INSTRUMENTS THIS LEAVES**, all of which paid for themselves: the trial
ladder ( every trial's merit, the field block alone, and the deviation from
`( 1 - a )‖R_0‖` ), the located X-point and the `psi_bnd` constraint per trial
with the element index, and the Schur complement's row norms and `cond_1`. They
are how six hypotheses were killed in a day rather than argued about.

### M-116

**`[solver] PicardSweeps` ON THE SIX COLD FAILURES: IT FIXES NONE OF THEM AND
TURNS THREE INTO SILENT WRONG ANSWERS.** `BORDERED-GLOBALISATION-PLAN.md` §0.2
records that **nothing in the tree has ever set this key** — no example, no test,
no benchmark — so every row of M-103 was taken at zero, and it was the cheapest
untried thing available: fourteen `meq-run` invocations and no code. It has now
been run. The pre-stage freezes `psi_ax`, `psi_bnd` and the plasma edge and
solves the UNBORDERED problem, where `PicardThenNewton` is legal because there
is no `psi_ax` unknown.

| case | sweeps 2 / 4 / 8 | `psi_ax` | vs its reference |
|---|---|---|---|
| B peaked ff′ | exit 2, exit 2, exit 2 | — | no change |
| B shaped | exit 2, exit 2, exit 2 | — | no change |
| **C MAST shaped** | **exit 0** | 3.519343e-01 | **2.905e+00** |
| D TCV | exit 2, exit 2, exit 2 | — | no change |
| **E diamagnetic** | **exit 0** | 1.030822e-01 | **2.306e-01** |
| **E shaped** | **exit 0** | 1.030887e-01 | **2.306e-01** |

**ZERO OF SIX FIXED, AND THREE OF SIX MADE WORSE IN THE ONE WAY THIS PROJECT
CARES ABOUT MOST.** Those three exit **0** with **no warning of any kind** —
the axis guard is silent, the source guard is silent, nothing in the output
says the answer is a different machine — and C shaped is out by a factor of
**3.9**. A run that fails is actionable; a green run reporting a different
equilibrium is not, and *Testing stance* exists for exactly this.

**AND THE SWEEP COUNT DOES NOT MATTER**: 2, 4 and 8 give bit-identical answers
on all three, so the pre-stage reaches its own fixed point and the bordered
solve inherits it whatever it is.

**E's PROBLEM IS ROOT SELECTION AND NOT GLOBALISATION, WHICH THIS SETTLES.**
The pre-stage reaches 1.030822e-01 and M-106's exact-seed run reached
1.025203e-01 — 0.5% apart, both 23% from the reference. Two routes that share
nothing find the same wrong branch, so no amount of step-length management is
the repair there.

### M-117

**WHAT THE SIX COLD FAILURES ACTUALLY ARE, AND NOT ONE IS A NON-FINITE
DIRECTION.** `BORDERED-GLOBALISATION-PLAN.md` §0.1 asks for this classification
and names the condition that would make the whole globalisation campaign moot:
if the bordered step's DIRECTION is non-finite, a line search can only choose
how far along it to go and no globalisation is the repair. **It does not hold.**
The border trace prints on failure under `--profile` and records
`directionFinite` per step; it is true on every step of all six.

| case | category | evidence |
|---|---|---|
| B peaked ff′ | **one fatal step** | 12 halvings exhausted on **1 of 27** steps |
| B shaped | one fatal step | 1 of 26 |
| C MAST shaped | **chronic** | 55 of 288 |
| D TCV | **X-point EXCURSION** | seed ( 0.7500, −0.5500 ) → ( 0.9296, −1.9875 ) → ( 1.0288, −2.1393 ) |
| E diamagnetic | X-point excursion | ( 1.1000, −0.6000 ) → ( 1.4909, **+0.9299** ) → ( 1.8960, +1.3010 ) |
| E shaped | X-point excursion | as E |

**TCV's NULL TRAVELS 1.6 m IN Z FROM ITS SEED AND E's FLIPS SIGN IN ONE STEP.**
That is not a step-length failure at all; it is the iteration leaving the saddle
XP-3 was told to follow, which `setXPointBoundary`'s own error message names.

**AND THE MERIT'S COMPOSITION EXPLAINS BOTH ENDS AT ONCE, WHICH IS THE FINDING
THAT UNIFIES M-111 THROUGH M-115.** Read off the same traces:

| | merit dominated by | what the X-point does |
|---|---|---|
| **C MAST** | the BORDER — `g·axis` **0.75** of a merit of 0.843 | crawls a millimetre an iteration; the plateau |
| **D TCV** | the FIELD — `g·axis` −1.31e-04 and `g·bnd` −1.74e-03 against ‖R‖ **2.953** | runs away unrestrained |

**The border constraints are three orders of magnitude apart in their share of
the merit between two machines in the same benchmark.** So there is no single
weight that can be right for both, and M-113's inversion, M-114's local optimum
and `BorderMeritWeight`'s one-sided saturation are all the same fact seen from
different machines: **MAST is the case where the border already dominates, so
raising it CANNOT do anything, and it is the only case the border weight was
measured on.** That deletion was taken on one machine at one end of the range.

### M-118

**Three defects in the bordered Newton's derivative supply, found by an audit
of what is analytic and what is differenced, and what each costs.**

The audit's own headline first, because it corrects the premise the audit was
commissioned under: **MEQ's bordered Jacobian is already almost entirely
analytic.** Of roughly twenty derivative quantities entering `J`, three are
differenced — `∂R/∂ψ_ax` at iteration 0 (always), the two normalisation columns
when the source does not supply `normalisationDerivatives()`, and the exterior
columns under `BorderColumn::Differenced`, which is the deliberate control.
Every constraint row, every corner entry and the current column are closed form.

**And the L2 jump is not an obstruction to an analytic derivative.** Six entries
in the whole system need `d/dx` of a point evaluation — the two `ψ_bnd` corner
entries against the X-point coordinates, and XP-3's four `∇q` ones — and all six
are exact arithmetic on the element's own polynomial today. The jump obstructs
the **constraint being differentiable at all**, which is a different problem and
one that analytic differentiation cannot touch. M-115's plateau hypothesis is
about the second, not the first.

| | defect | reachable | what it cost |
|---|---|---|---|
| 1 | `NormalisedRotatingSource` never applied `currentScale()` | `[source] Type = "rotating"` + `Normalised` + `PlasmaCurrent`, no guard anywhere | `∫F/r` independent of λ, so the current constraint cannot be satisfied by the unknown that exists to satisfy it — while `assembleCurrentColumn()` and `cornerEntry( I, I )` report a sensitivity of `scaledF/λ` that does not exist |
| 2 | `NormalisedRotatingSource::normalisationDerivatives()` unwritten | every bordered rotating solve | the two normalisation columns differenced, and a difference perturbs the normalisation, which moves the edge, so it straddles a kink |
| 3 | `assembleCurrentNormalisationCorner()` returned a silent **zero** on a source that refused | any source without (2) | see below |

**(3) IS THE ONE WORTH A NUMBER, AND THE NUMBER IS NOT WHAT WAS PREDICTED.** Its
sibling `assembleNormalisationColumn()` falls back to a central difference in
exactly the same situation: one degraded and the other deleted. The repair is a
difference of `assemblePlasmaCurrent()` on the shared step, and the capability
probe is hoisted out of the quadrature loop — asking inside it would return a
partial sum over the elements already walked if a source ever refused late.

Measured on `HighBetaConvergence`'s rectangle with `setPlasmaCurrent`, a
refusing source against the same source unwrapped:

| corner | n | Newton | final residual | `psi_ax` | scale |
|---|---|---|---|---|---|
| assembled | 8 | **4** | 2.0478e-15 | 5.242280976e-03 | 7.713043e-05 |
| differenced | 8 | **5** | 1.1617e-16 | 5.242280976e-03 | 7.713043e-05 |
| assembled | 16 | **4** | 2.1119e-15 | 5.242827733e-03 | 7.708217e-05 |
| differenced | 16 | **5** | 1.7912e-17 | 5.242827733e-03 | 7.708217e-05 |

**AND THE SAME TABLE WITH THE DEFECT PUT BACK, WHICH IS THE HALF THAT MATTERS:**
the zeroed corner reads **42 iterations at both mesh sizes**, against 4 — and
**still converges**, to 6.28e-13 and 2.81e-13. So on this problem a deleted
corner costs a factor of ten in work and nothing in the answer.

**That mutation changed what the test asserts.** The case was written asserting
that a deleted corner fails to converge within its cap, on the arithmetic that
a geometric 0.8 a step needs about 130 iterations to cross twelve decades —
which is `CLAUDE_FB.md`'s measurement on a harder problem and is not this one.
**The case passed against the defect it was written to catch**, and only running
the mutation found that out. It now asserts on the iteration count, which sits
an order of magnitude clear of both behaviours: 5 against a bound of 10, where
the defect gives 42.

**A fourth finding is recorded and NOT repaired**, because repairing it needs an
instrument that does not exist. The axis row drops its position term by the
envelope theorem, justified in `GradShafranov.cpp` by *"`grad_bar( psi ) = r q`,
so `grad( psi_h )( x* ) = 0` at a zero of `q_h` **IDENTICALLY**"*. That relation
is continuous: the discrete flux equation makes `r q_h − ∇̄ψ_h` the local lifting
of the trace jump, so `∇ψ_h( x* )` is `O( h^k )` rather than zero — which is the
same fact as `q_h` converging a full order better than `∇ψ_h`, and is why the
mixed method exists. `cornerEntry()`'s XP-3 arm says exactly this about the same
identity (*"whose identity is only weak"*) 380 lines away, and the two cannot
both be right. The missing piece is
`∇ψ_h( x* )ᵀ ( ∇q_h )⁻¹ ( flux shape at x* )`, and **every factor is already
assembled** — `xFluxJacobian` is `∇q_h`, `xFluxShape` is `∂q_h/∂u`. It is a
fourth candidate for XP-3's observed order of 1.664 against XP-2's 1.667,
alongside the three `CLAUDE_FB.md` already lists. It is left alone because a
Jacobian term moves which branch a marginal case selects — M-26 measures 9.4% in
`max ψ_h` — and the border has no Jacobian-against-difference case, where the
field block has had one in `NewtonConvergence.cpp` all along.

**The transferable part**: `[source] Type = "rotating"` has now produced this
shape three times — `ConfineToPlasma`, `PlasmaCurrent` and
`normalisationDerivatives()`. Every accepted key on Config's rotating branch is
worth checking against what `meq::NormalisedRotatingSource` actually reads.

### M-119

**`Globalisation::BorderedPicardThenNewton` on the six cold failures: 0 of 6
close, and the prediction that it would close three is FALSIFIED.**

The shipped configurations, one `build/meq-run` each, `MKL=1 OMP=8`, machine
idle (load average 0.03 before the first). **Both arms come from ONE run**: the
driver's ladder tries `Globalisation::None` first and falls to the rung on the
observed failure, so the comparison is the same binary, the same mesh and the
same guess, with no possibility of the attribution error M-102 records.

| case | Newton arm | | Picard arm | | |
|---|---|---|---|---|---|
| | its | min ‖r‖/‖r₀‖ | its | min ‖r‖/‖r₀‖ | |
| `machine-b-ffprime` | 13 | **4.4776e-01** | 21 | 5.2024e-01 | worse |
| `machine-b-ffprime-shaped` | 12 | **4.4773e-01** | 23 | 5.2084e-01 | worse |
| `machine-c-mast-shaped` | 163 | 6.6615e-01, ending **3.0183e+06** | 35 | **3.1903e-01** | **much better, still fails** |
| `machine-d-tcv` | 6 | **3.3509e-01** | 8 | 3.3742e-01 | level |
| `machine-e-diamagnetic` | 9 | **2.5903e-01** | 10 | 6.1105e-01, ending 6.7623e+00 | worse |
| `machine-e-diamagnetic-shaped` | 9 | **2.5806e-01** | 9 | 6.1116e-01, ending 6.7468e+00 | worse |

**All six still exit 2.** Five fail with the same message in both arms — *no
damping of the bordered Newton step gave a finite residual*.

**WHAT WAS PREDICTED AND WHAT HAPPENED ARE DIFFERENT SETS, WHICH IS THE POINT OF
HAVING WRITTEN THE PREDICTION DOWN FIRST.** M-117's classification said B and
B-shaped (one fatal step) and C-shaped (chronic) were the shapes a field-block
preconditioner could address, and D/E/E-shaped were X-point excursions it could
not. Three of six was recorded as the *maximum*. The result is **one** of six
showing improvement, and it is not one of the two the argument leaned on: B and
B-shaped are made WORSE, and both arms there plateau at about 0.5 of the initial
residual and then wander without descending, which is neither a fatal step nor
something a linearisation choice reaches.

**THE ONE INFORMATIVE ROW IS `c-mast-shaped`, AND ITS VALUE IS THE CHANGED
FAILURE MODE RATHER THAN THE BETTER RESIDUAL.** Under Newton it spends 163
iterations and diverges to 3.0e+06. Under Picard it reaches 3.19e-01 in 35 and
then throws

```
the bordered Jacobian is singular in ( psi_ax, psi_bnd, a )
```

— a guard that has been in the tree since FB-3 (`64cca6e`) and that the Newton
arm never reaches. **With the field block replaced by `A_lin`, which is
unconditionally invertible, what is left singular is the BORDER'S Schur
complement.** That is `BORDERED-GLOBALISATION-PLAN.md` §0.1's first row — *a
non-finite direction; the Jacobian or its elimination is the defect* — reached
by measurement rather than by argument, and it says the obstruction on this case
is the border and not the field. §6.1's degenerate axis row and M-115 are where
it points.

**AND THE HANDOFF NEVER HAPPENS THE WAY IT WAS DESIGNED TO.** The tolerance is
1e-3 relative and no case gets below 2.5e-01, so on all six the phase change is
the budget or the line-search rescue, never the tolerance. A handoff tolerance
is the right control for a case that is merely slow; none of these is.

**The transferable part.** The bordered Picard is now a measured instrument
rather than a candidate repair: replacing the field block with something that
cannot be singular is a way of ASKING whether a failure is the field's, and on
five of six the answer is no. Four of six getting worse is itself a finding —
a Picard direction is longer and less well aimed, so on the X-point excursion
cases it walks the null further, which is E's residual growing 6.7× where
Newton's merely stalls.

### M-120

**The border's Jacobian-against-difference case, and what it found: the axis
row's envelope argument is false, and the omission is 0.03 to 0.16 per cent of
the row it sits in.**

`NewtonConvergence.cpp` has checked the assembled FIELD Jacobian against a
central difference of the assembled residual since stage 4. **The border never
had the equivalent**, and M-118 records two defects one would have caught.
`tests/convergence/BorderJacobian.cpp` is it.

**THE IDENTITY THE AXIS ROW RESTS ON IS FALSE, AND `q_h` IS NOT THE REASON.**
`C_Ax = s - psi_h( x* )` with `x*` a root of `q_h`; the position term is dropped
by the envelope theorem, which is exact only if `grad psi_h( x* ) = 0`. At the
located axis, on `HighBetaConvergence`'s rectangle at `k = 2`:

| n | h | \|q_h(x*)\| | \|grad psi_h(x*)\| |
|---|---|---|---|
| 8 | 0.1000 | 7.9213e-17 | **4.9651e-03** |
| 12 | 0.0667 | 1.5293e-16 | **1.1467e-02** |
| 16 | 0.0500 | 1.2493e-14 | **8.3029e-04** |
| 24 | 0.0333 | 1.6318e-13 | **2.5528e-03** |

`q_h` is at round-off there — it *is* the root, so that column is the control
saying the finder converged. `grad psi_h` is eleven to thirteen orders larger.

**AND `psi*` IS THE RIGHT FIELD TO ASK, WHICH BUYS AN ORDER AND NOT AN
IDENTITY.** `psi*` is the field whose gradient *is* the solved flux, so if the
envelope argument holds anywhere it is there:

| n | h | \|grad psi_h(x*)\| | \|grad psi*(x*)\| | ratio |
|---|---|---|---|---|
| 8 | 0.1000 | 4.9651e-03 | 1.2341e-03 | 0.249 |
| 12 | 0.0667 | 1.1467e-02 | 9.6484e-04 | 0.084 |
| 16 | 0.0500 | 8.3029e-04 | 1.5539e-04 | 0.187 |
| 24 | 0.0333 | 2.5528e-03 | 1.1366e-04 | 0.045 |

The ratio falls with `h`, which is the extra order showing — and `grad psi*` is
still 1e-04, nowhere near round-off. The reason is structural rather than a
matter of resolution: the local post-processing solves
`( grad psi*, grad v )_K = ( r q_h, grad v )_K`, so `grad psi*` is the L2
projection of `r q_h` onto the GRADIENTS of `P^(k+2)( K )` — the nearest
gradient field to `r q_h` and not `r q_h` itself. `r q_h` is not a discrete
gradient, so a residual survives however fine the mesh.

**So there are three repairs and they are not the same one.** (a) evaluate
`psi*` at the root of `q_h`: one order smaller, still wants the correction.
(b) add the term to the row as it stands: exact for what `psi_ax` currently
MEANS, one 2x2 solve from pieces already assembled. (c) define `x*` as a
critical point of `psi*` and evaluate `psi*` there: then `grad psi*( x* ) = 0`
is the DEFINING equation and the envelope theorem is exact by construction, with
no correction term at all. (c) is the self-consistent design and it is not free
— it redefines `psi_ax`, which `recoverPeak()` documents as chosen for
differentiability, so it carries M-26's branch-selection risk; and the row
becomes `-( d psi*/du )( x* )`, the local reconstruction operator applied to the
shape functions, which is linear and per-element but is machinery MEQ does not
expose. `postProcess()`'s 0.62 s is NOT the objection it looks like: only the
axis element's reconstruction is wanted, not the mesh's.

**None of the three removes the face jump below.** `psi*` is element-wise, so
the jump shrinks by an order — `O( h^(k+2) )` against `O( h^(k+1) )` — and does
not go away.
**No rate is quoted and none should be**: where in its element the axis lands
changes with the mesh, and the column is not monotone. The finding is that it is
not zero, which an identity does not survive.

**WHAT THE OMISSION IS WORTH, COMPUTED RATHER THAN DIFFERENCED.** MEQ's flux
space is a scalar collection at `vdim = 2`, so `dq_d/dc_(j,d') = phi_j
delta_(dd')` and, from `q_h( x* ) = 0`, the dropped entry on flux dof `(j, d')`
is `w_d' phi_j( x* )` with `w = ( grad q_h )^-T grad psi_h`. The row MEQ does
assemble has norm `| phi( x* ) |` over the potential dofs, so the ratio is
`| w |` up to a shape-norm ratio the case reports rather than assumes:

| n | \|grad q_h\| | \|neglected\| | of the kept row |
|---|---|---|---|
| 8 | 7.6937e+00 | 1.0364e-03 | **1.107e-03** |
| 12 | 7.6584e+00 | 1.4468e-03 | **1.560e-03** |
| 16 | 7.6806e+00 | 2.4331e-04 | **2.730e-04** |
| 24 | 7.6701e+00 | 3.0546e-04 | **3.509e-04** |

Every factor is already assembled for XP-3's own rows — `xFluxJacobian` is
`grad q_h`, `xFluxShape` is `dq_h/du` — so the correction is one 2x2 solve and a
scatter. **The case is RED until it is written**, per the testing stance, and
reads `GradShafranovSolver::axisRowCarriesEnvelopeTerm` so that writing the term
turns it green rather than needing the case edited.

**AND THE FIRST VERSION OF THIS MEASUREMENT WAS WRONG IN A WAY THAT BECAME THE
THIRD FINDING.** It differenced the term: perturb the flux, relocate `x*`,
difference `psi_h( x* )`. The readings were -5.760e-04, -1.158e-03, -2.321e-03
as the step halved — **doubling**, so the NUMERATOR was constant. A constant
numerator is a discontinuity, not a derivative. Printing the element index
settled it:

| step | numerator | \|x* moved\| | element at -a | element at +a |
|---|---|---|---|---|
| 9.4434e-03 | -1.087822e-05 | 1.8698e-04 | **275** | **242** |
| 4.7217e-03 | -1.093199e-05 | 1.0292e-04 | **275** | **242** |
| 2.3609e-03 | -1.095862e-05 | 6.4947e-05 | **275** | **242** |

The two arms land in **different elements at every step**, and `psi_h` is L2, so
the difference reports the face jump however finely it is probed. **That is
M-115's surviving hypothesis — a point evaluation of a broken field jumping
across a face — observed on the AXIS row**, which M-117 notes carries 0.75 of
MAST's merit and has no instrument at all where the X-point has four. It is now
an asserted property: the numerator does not fall with the step.

**The transferable part**: differencing a constraint whose evaluation point is
*located* rather than prescribed measures the location's discontinuity, not its
derivative, and the two are told apart by the element index rather than by the
number. The closed form has no such trouble because it never leaves one element.

### M-121

**The axis row's envelope term, written — and assembling it is measured
HARMFUL, so it is off by default.**

M-120 established that the axis row drops a term worth 0.03 to 0.16 per cent of
the row it sits in, and that every factor was already assembled for XP-3's own
rows. It is now written. `GradShafranovSolver::AxisRow::WithEnvelope` turns it
on and `PositionDropped` is the default.

**The sign trap is `DarcyForm`'s `-q` and it bites twice**, which is worth
setting out because a sign error here does not give a wrong answer — it costs
convergence and nothing else. The unknown's flux block holds `-q`, so with
`q_d = -sum_j u_(j,d) phi_j`:

```
dq_d/du_(j,d')    = -phi_j delta_(dd')
d( x* )/du_(j,d') = -( grad q )^-1 dq/du = +( grad q )^-1 e_d' phi_j
b_(j,d')          = -w_d' phi_j( x* ),   w := ( grad q )^-T grad psi_h
```

**IT IS CORRECT, BY THE ONLY TEST THAT CAN SEE A JACOBIAN** — the equilibrium
does not move. `psi_ax` 8.266004e-02, `psi_bnd` 3.237932e-02 and the X-point at
(1.093, -0.604) to every printed digit, on both libraries and under both
settings. A wrong sign would raise the iteration count; the counts below do not
behave like that.

**BUT WHAT IT IS WORTH IS LIBRARY-DEPENDENT, AND THE FIRST READING OF THIS TABLE
CREDITED THE TERM WITH SOMETHING THE LIBRARY DID.** Taken first against
`libmfem.a` of 2026-09-15 and re-taken against 2026-09-17, both arms each time:

| XP-3 | 09-15 dropped | 09-15 envelope | 09-17 dropped | 09-17 envelope |
|---|---|---|---|---|
| bootstrap iterations | 15 | 15 | **14** | **15** |
| sweep 1 iterations | **4** | **3** | **3** | **3** |
| sweep 3, final residual | 1.702e-13 | 5.551e-16 | 1.511e-14 | 2.532e-15 |
| bootstrap, final residual | 1.273e-14 | 1.736e-16 | — | — |

**The 4-to-3 at sweep 1 was real on 09-15 and does not reproduce on 09-17,
because the newer library reaches 3 on its own.** On today's library the term
COSTS a bootstrap iteration, 14 against 15, and buys a final residual six times
deeper. That is a much more marginal trade than the first reading of this table
suggested, and it is recorded rather than quietly restated: the comparison was
internally valid on the library it was taken on, and it did not survive the
library moving under it. *A measurement of a Jacobian change is a measurement
against one library, and this file already says that about suite times.*

**AND IT COSTS A PHYSICAL CASE, WHICH IS WHY IT IS NOT THE DEFAULT.**

| | dropped | with envelope |
|---|---|---|
| `FreeBoundaryCoupling`, 09-15 | **5 of 5** limiter radii | **4 of 5** — R = 1.18 no longer converges |
| `FreeBoundaryCoupling`, 09-17 | **passes** | **fails** — the coil-free control runs 2 rows of 3 |
| `HighBetaConvergence` assembled corner, final residual (09-15) | 2.048e-15 / 2.112e-15 | 3.017e-13 / 3.166e-13 |
| XP-3 bootstrap, residual by step 8 (09-15) | 6.0e-05 | 1.8e-03 |

**The harm is the part that DOES reproduce across the library change**, with a
different symptom — 09-15 lost a limiter radius outright, 09-17 loses a row of
the coil-free control — and it is what decides the default.

The endgame improves and the approach gets worse. That shape, and the lost
radius, are what a more exact derivative of a function that is **not
differentiable** would produce — and M-120's third finding is exactly that
`C_Ax` is discontinuous in the flux, jumping by the L2 face jump whenever the
axis crosses a face. Sharpening the smooth part's derivative aims the step more
precisely at a residual that does not follow it there, and the line search
rejects.

**So the two halves of M-120 are in tension and the second wins.** The term is
real, it is computable, it is correctly signed, and the constraint it
differentiates is the wrong kind of function to be differentiating harder. **The
route worth taking is the one that makes the constraint continuous first** —
M-115's local reconstruction, or the `psi*` variants set out at the end of
M-120, of which (c) is the self-consistent one: define `x*` as a critical point
of `psi*` and evaluate `psi*` there, where the envelope theorem is exact by
construction and there is no correction term to add at all.

**The transferable part**: *this was the third time in this campaign that a
repair which is right by derivation was measured to be worse in practice*, after
`PicardSweeps` (M-116) and the bordered Picard rung (M-119). All three share a
shape — the derivation assumes a smoothness the discrete problem does not have.

### M-122

**The MAST-U forward reference over a resolution ladder, on a quiet machine —
and the accuracy ceiling it sets for anything compared against it.**

freegsnke's `example02` with the notebook removed and a grid sweep round it:
MAST-U, **forward**, `constrain=None`, twelve prescribed active currents,
`ConstrainPaxisIp` at `paxis = 8e3`, `I_p = 6e5`, `alpha_m = 1.8`,
`alpha_n = 1.2`, `gs_operator_order = 4`, tolerance 1e-9. Median of three
repeats per level; load 0.08 before, held under 0.5 for a sustained 300 s first.

| grid | median s | peak MB | Picard | NK | rel | `psi_axis` | `I_p` |
|---|---|---|---|---|---|---|---|
| 33×65 | 0.05 | 357 | 15 | 11 | 6.08e-10 | 9.18724696e-02 | 6.00000e+05 |
| 65×129 | 0.25 | 468 | 15 | 11 | 3.06e-10 | 9.18720993e-02 | 6.00000e+05 |
| 129×257 | 1.05 | 1291 | 15 | 10 | 1.17e-10 | 9.18649766e-02 | 6.00000e+05 |
| 257×513 | 6.88 | 6598 | 15 | 12 | 2.36e-10 | 9.18641969e-02 | 6.00000e+05 |

**THE ITERATION COUNT IS MESH-INDEPENDENT — EXACTLY 15 PICARD AT EVERY LEVEL**,
and 10 to 12 Newton–Krylov. That is the property that makes this a fair
comparison point: the reference's work per level is flat, so a ratio against
MEQ's is a statement about the codes rather than about how each degrades.

**AND THE REFERENCE'S OWN ACCURACY CEILING IS ABOUT 1e-05, WHICH IS THE NUMBER
TO READ THIS TABLE FOR.** `psi_axis` moves **7.8e-05** relative between 65×129
and 129×257 and **8.5e-06** between 129×257 and 257×513 — non-monotone, so no
clean order can be read off it. **Agreement better than about 1e-05 on this
case is therefore meaningless**, exactly as `tools/README.md` records for the
freegs4e benchmark saturating at 1.4e-04. Do not report a MEQ figure below that
as accuracy; it is the reference's noise.

**257×513 IS THE TOP RUNG ON THIS MACHINE AND NOT A CHOICE.** Peak resident
memory is 6.6 GB there, against 1.3 GB one level down — very nearly ×5 per
level, so 513×1025 would want about 26 GB and this is WSL2.

**A DEFECT IN THE HARNESS, FOUND BY RUNNING IT.** `race_nke.py`'s
`--require-quiet` first compared the load average before AND after and would
have **refused every successful race**: the ladder above left the load at 1.36
having started at 0.08, entirely its own doing, so an after-check cannot tell a
busy machine from a busy measurement. The gate now samples just **before each
level**, where the only thing running is the harness, with one busy arm's worth
of slack. *A gate that always fires is worse than no gate, and only running it
shows which it is.*

### M-123

**`BorderRegularisation` on MAST-U: it removes the failure it was built for, and
the case still does not converge cold.**

`examples/mastu-nke.toml`, freegsnke's `example02` machine as MEQ solves it,
9361 elements at `k = 2`, `OMP = MKL = 8`, **no `[initialguess]`**.

**WITHOUT IT the solve dies at iteration 0** with *"the bordered Jacobian is
singular in ( psi_ax, psi_bnd, a )"* — on both the Newton arm and the bordered
Picard one, 0 iterations each.

**WITH `BorderRegularisation = 1e-3` (and 1e-2, and 1e-1) the throw is gone**,
and what is behind it is this:

| arm | behaviour |
|---|---|
| `Globalisation::None` | stalls **completely** at `7.539822e-01`, 5.44× the initial residual — **identical to seven digits at iterations 40, 80, 120, 160 and 200**. The iterate stops moving at all |
| bordered Picard → Newton | reaches **1.51e-02** relative at best, then creeps UPWARD at a steady 1.001 per iteration to the 400 cap |

**THE OPTION DID WHAT IT WAS DESIGNED TO DO.** It is the first thing this
campaign has built that did — M-116, M-119 and M-121 are three repairs that were
right by derivation and measured worse — and the difference is that it was aimed
at a mechanism M-119 MEASURED (`c-mast-shaped` throwing on a singular border once
the field block was made non-singular) rather than at one derived from first
principles.

**AND TWO THINGS IN THAT TABLE ARE NEW.**

**The bordered Picard rung is the arm that makes progress here**, reaching
1.51e-02 where plain Newton does not move at all. M-119 has it making four of
six cold failures WORSE, so this is the first case where it is the better arm —
one data point, and it is the opposite sign from that measurement.

**The creep is an instability signature, not a stall.** A steady 1.001 per
iteration is ×1.22 over 200, which is M-104's vertical mode's shape — that
measured DIII-D's forward Picard growing ×10 per 20 passes in `Zaxis` — at about
a hundredth of the rate. A residual that sits still is a dead iteration; one
that grows geometrically is a mode being amplified.

**WHAT IS NOT YET TESTED IS THE THING MOST LIKELY TO FIX IT.** The config has no
`[initialguess]`, and `examples/machine-g-mastu.toml`'s own comment is *"NOT A
NICETY. A cold start on a free-boundary machine case wanders and does not
converge; the guess is part of the problem statement."* M-112's cold-forward
audit measures MEQ failing cold on 9 of 14 machine configurations. **So this
result is consistent with "no guess" and is not evidence about the solver**;
`mkexactguess.py` is the bounded work that would settle it.

### M-124

**MAST-U converges from its own configuration file, and what it took was
spreading `I_p` over an ELLIPSE.**

`examples/mastu-nke.toml`, freegsnke's `example02` machine, 9052 elements at
`k = 2`, `OMP = MKL = 8`, `[initialguess] Type = "conductors"`. Nothing is
supplied from outside the file — no `.gf`, no `mkexactguess.py`.

**IT TAKES `BorderRegularisation` AS WELL, AND THAT IS NOT A DETAIL.** Every row
below carries `BorderRegularisation = 1.0e-3`. Without it this machine fails
whatever the guess: cold it dies at iteration 0 with a singular bordered
Jacobian, which is M-123; WITH the ellipse it dies differently, at *"no damping
of the bordered Newton step gave a finite residual — psi_ax through zero"*. So
M-123's repair and this one are necessary together and neither is sufficient,
and the committed example was measured failing with only one of them.

| guess | Newton its | `psi_ax` | `psi_bnd` | X-point | from seed | wall |
|---|---|---|---|---|---|---|
| conductors alone | 17+8+7+9 = 41, support NEVER settled | **3.257887e-01** | 1.326955e-02 | ( 0.521, 0.475 ) | **1.574 m** | 318.4 s |
| **+ ellipse ( 0.50, 0.90 )** | 7+3+2 = 12, settled in 3 | **9.162567e-02** | 2.864064e-02 | ( 0.598637, −1.097187 ) | **1.578e-04 m** | **292.8 s** |
| + ellipse, DEFAULT semi-axes ( 0.425 round ) | 5+3+2+2 = 12, settled in 4 | **9.162567e-02** | 2.864064e-02 | ( 0.598637, −1.097187 ) | 1.578e-04 m | 314.5 s |
| + FILAMENT at the same point | — | **fails** | | | | |
| freegsnke reference | 25–27 Picard | 9.187209931e-02 | 2.867695588e-02 | ( 0.598480, −1.097170 ) | — | 0.25 s at 65×129 |

**Against the reference: `psi_ax` to 2.68e-03 relative, `psi_bnd` to 1.27e-03,
the X-point to 1.6e-04 m.** The residual history of the last sweep is
3.010466e-02 → 4.401979e-11 → 1.953785e-15, and the axis lands at
( 0.9489, −0.0000 ).

**THE CONDUCTORS' FIELD ALONE FINDS A BRANCH THAT IS NOT ONE.** Its axis is at
( 2.6509, −0.9134 ) — outside the machine — with `psi_ax` 3.5× the reference.
That is M-26's hazard reached by a guess with the right scale and no plasma in
it, and it is the clearest instance in this tree of a converged solve being the
wrong answer rather than a failed one.

**AND THE COLUMN'S SHAPE DOES NOT MATTER, WHICH IS THE RESULT THAT MAKES THIS
USABLE.** A column tuned to MAST-U's own plasma — semi-axes ( 0.50, 0.90 )
against a separatrix reaching `r = 0.26` to `1.4` and `|z| = 1.1` — and a round
one at the DEFAULT `0.5*CentreR = 0.425` reach the same `psi_ax`, `psi_bnd` and
X-point **to every printed digit**, in the same 12 Newton iterations. So the
file costs the user ONE number, the guessed major radius, and MEQ supplies the
rest.

**THE FILAMENT IS WHAT THIS REPLACES AND IT FAILS OUTRIGHT.** `meq::filamentPsi`
diverges logarithmically at the filament, so `I_p` carried on one at the guessed
axis is an unbounded spike exactly where `locateAxisPoint()` has to look:

| distance from the guessed axis | filament | ellipse ( 0.45, 0.70 ) |
|---|---|---|
| 0.300 m | 1.846806e-01 | 1.551829e-01 |
| 0.030 m | 4.111789e-01 | 1.367698e-01 |
| 0.003 m | 6.666827e-01 | 1.329906e-01 |
| 0.0003 m | 9.281104e-01 | 1.325976e-01 |

against a reference `psi_axis` of 9.187e-02. **A factor of 5.0 across four
decades of approach for the filament against 1.17 for the ellipse**, and the
filament is still climbing.

**AND AN ELLIPSE RATHER THAN A RECTANGLE, WHICH `meq::Coil` ALREADY IS AND
WHICH WOULD ALSO BE BOUNDED.** A uniform current density over a rectangle
carries a logarithm in its second derivatives at each of four corners. Those are
artefacts of the shape and nothing in the equilibrium puts them there; an
ellipse is also what a plasma column actually is.

**`meq::ellipsePsi()`, AND ITS QUADRATURE IS TWO RULES RATHER THAN ONE.**
Measured in `tests/unit/CoilsTests.cpp`:

| | |
|---|---|
| `Delta* psi` inside the column against `-mu0 r I/( pi a b )` | **1.5e-05 to 4.7e-05** relative — the second difference's own truncation |
| `Delta* psi` outside it | **8.4e-07 to 1.8e-05** of the interior value |
| a shrinking column against its filament | rate **1.98, 1.99, 2.00** — the quadrupole moment, `O( a^2 )` |
| the step across the boundary, extrapolated from each side | **8.7e-06 to 9.2e-04** relative |
| cost | **66 us** per field point at the default order |

**ONE RULE WAS TRIED FIRST AND FAILS OUTSIDE, SEVERELY RATHER THAN GRADUALLY.**
Sweeping chords from the field point is what kills the kernel's `log( 1/rho )`
inside the ellipse, and it is the right rule there. Seen from OUTSIDE, the
ellipse subtends a CONE, so a uniform rule in `theta` spends its nodes mostly on
rays that miss — a distant or small ellipse gets a handful, sometimes none. With
that rule alone the far field stopped converging to its filament at 1e-04 and
then went BACKWARDS (errors 1.520e-03, 1.059e-04, 3.577e-03, 1.505e-02 over four
halvings), and `Delta*` outside the column read **order the interior value with
an erratic sign**. The second rule sweeps the ellipse's own disc, where the
integrand is analytic.

**AND HALVING THE ANGLES IS A REAL 2x THAT THE SEAM REFUSES.** 66 us to 34 us,
with the far field and `Delta*` unmoved — both are spectral out there — and the
step across the boundary going from 9.2e-04 to **3.5e-03**.

### M-125

**`UpDownSymmetry` is exact on a mirror-symmetric mesh, and MAST-U's mesh is not
one.**

The projection averages every dof with the dof at its own reflection, so it
needs a mesh whose dofs are mirror-paired. `tests/convergence/UpDownSymmetry.cpp`
is the acceptance and it has two cases.

**ON A SYMMETRIC MESH IT IS THE IDENTITY, TO ROUND-OFF, IN EVERY BLOCK.** A
bordered solve of an entirely even problem on the standard box in
QUADRILATERALS, run with the projection and without it:

| `k` | `n` | dofs | `psi_ax` | `d psi_ax` | `d psi` | `d flux` |
|---|---|---|---|---|---|---|
| 1 | 8 | 256 | 9.816290e-01 | 6.4e-15 | 7.9e-15 | 7.4e-15 |
| 1 | 12 | 576 | 9.617180e-01 | 6.9e-16 | 1.8e-15 | 1.6e-15 |
| 2 | 8 | 576 | 1.045507e+00 | 2.8e-15 | 3.4e-15 | 3.6e-15 |
| 2 | 12 | 1296 | 1.042921e+00 | 2.1e-16 | 1.9e-15 | 3.2e-15 |

The FLUX column is the one with teeth: `q_r` is even in `z` and `q_z` is odd, so
the map carries a sign per component, and with it backwards the projection lands
on the ANTIsymmetric subspace where the only even field is zero.

**ON MAST-U IT REFUSES, AND IT IS RIGHT TO.** `examples/mastu-nke.msh` has
**3624 of 4735 vertices with no mirror partner**. The machine is not the
problem: its 23 `[[coils]]` are mirror-paired to the last digit, one pair
excepted whose currents are **equal and opposite at 1.9107e-04 A** against a
total of 2.310105e+06 A — 8e-11 of it, which is the 2.4e-10 the applied field
was measured symmetric to. **The GEOMETRY is symmetric and gmsh's triangulation
of it is not.** `tools/mesh/halfdisc.py --symmetric` meshes one half and
reflects it, which is M-127; the refusal below is what a mesh made without it
still gets.

**AND `MakeCartesian2D`'s TRIANGLES ARE NOT THE SYMMETRIC MESH EITHER**, which
is the trap to know before reaching for one: it splits every cell along ONE
diagonal, so a box symmetric in `z` has a triangulation that is not. The refusal
fires on it, naming the element at ( 0.633333, −0.500000 ). The QUADRILATERAL
variant has no diagonal to choose.

**TWO DEFECTS IN THE PROJECTION WERE FOUND BY BUILDING THAT CASE AND NEITHER
COULD HAVE BEEN FOUND BY RUNNING MAST-U**, which throws before reaching them:

* **`ProjectCoefficient( x -> x(0) )` is not the route to a dof's coordinates.**
  It carries `MFEM_VERIFY( VectorDim() == 1 )`, so the FLUX space — the one
  space with vdim 2 and the only one whose map needs a sign — throws out of it;
  and the trace space's dofs are on FACES, which an element-wise projection does
  not reach at all.
* **A coordinate does not identify a dof.** Both volume spaces are L2 on the
  closed Gauss-Lobatto basis, so a dof sits ON the element boundary and every
  element meeting a vertex has its own dof there. A table keyed on position
  holds one of them. The element is now matched first, by centroid, and the
  dofs paired inside the matched pair.

**THE FIRST VERSION OF THE SYMMETRIC CASE READ 0.000e+00 IN EVERY BLOCK AND WAS
VACUOUS.** It used a plain linear `solve()`, and the projection is wired into
`solveWithNormalisation()` and nowhere else — `psi_ax` being an unknown is what
makes the branch a question at all. **An exact zero from a routine that
reassociates sums is not a pass; it is a routine that did not run.**

### M-126

**WHERE THE OTHER FIVE CORES GO ON A BORDERED FREE-BOUNDARY SOLVE.** The
question was why a MAST-U run at `OMP = MKL = 8` sits at 250–290% of a possible
800%. Taken on `examples/machine-f-diiid.toml` as a five-second proxy for the
same path — free boundary, X-point, 10 exterior modes, 4848 elements, `k = 2`,
PARDISO, `AssemblyMode::Threaded` — on an 8-core Ryzen 7 3800X.

**FIRST, A THIRD OF THAT PERCENTAGE IS NOT WORK.** `OMP_WAIT_POLICY=passive`
deletes two fifths of the CPU time for about 1.5% of the wall clock — four
interleaved pairs on a quiet machine, medians 4.965 s default against 5.040 s
passive, 14.6 s of user time against 8.6 s, `passive` slower in all four pairs
and the spread within each arm about 2%. **So it is very nearly free and not
exactly free**, and a first reading that said "free" was six runs on a
contended machine where the difference is below the noise:

| | wall | user | %CPU |
|---|---|---|---|
| `OMP=8 MKL=8` | 5.45, 5.44 s | 16.51, 16.27 s | 316%, 312% |
| `OMP=8 MKL=8` **passive** | 5.39, 5.46 s | 9.61, 9.57 s | **197%, 193%** |
| `OMP=8 MKL=1` | 6.12, 6.48 s | 12.08, 12.63 s | 208%, 207% |
| `OMP=8 MKL=1` **passive** | 6.38, 6.33 s | 8.84, 9.11 s | 153%, 157% |
| `OMP=1 MKL=8` | 7.37, 7.12 s | 12.00, 11.94 s | 173%, 178% |
| `OMP=1 MKL=8` **passive** | 7.07, 7.13 s | 8.23, 8.36 s | 126%, 126% |

**6.8 s of a 16.4 s CPU budget is barrier spin**, and the two pools contribute
independently and about equally — 3.3 s from MEQ's own OpenMP regions, 3.7 s
from MKL's, and 3.3 + 3.7 = 7.0 against the 6.8 measured together.
`GOMP_SPINCOUNT=1000` reads the same as `passive`. **So the utilisation figure
that prompted the question is inflated by the runtime waiting, and the honest
number is about 190%.**

**AND IT IS WHAT MAKES A PER-LEG `cores` READING MEAN ANYTHING**, which is the
reason to set it that has nothing to do with the trade above. Under the default
policy the same run reports `constraint location` at **5.76 cores** and `I_p` at
**8.95** — legs that read 1.01 and 1.03 under `passive`, and that no thread but
one enters. The CPU is the other seven spinning through them. A serial leg
therefore reads HIGH rather than low, which is the opposite of the mistake the
column invites.

**SECOND, THE WALL CLOCK: EIGHT CORES BUY 1.43x AND ONE LEG OWNS THE GAP.**

**THE WHOLE-RUN BUDGET, WHICH CLOSES AGAINST THE WALL CLOCK**, quiet machine,
`passive`, three support sweeps. `cores` is `cpu/wall` over the leg, so it is
the mean number of cores busy in it; the legs sum to the solve phase and the
three phases sum to the wall:

```
MEQ: wall 4.753 s = setup 0.020 + solve 3.927 + output 0.806
MEQ: where the run went -- wall 4.753 s, cpu 9.138 s, 1.92 cores on average
                              wall s   share   cores      calls
  setup                       0.020    0.4%    1.00          1
  solve                       3.927   82.6%    2.12          3
  output                      0.806   17.0%    1.00          1
    of which postProcess      0.609   12.8%    1.00          1
  residual                    0.305    6.4%    6.75         28
  gradient                    0.254    5.3%    4.78         12
    of which ComputeH         0.152    3.2%       -         12
  trace factorisation         0.230    4.8%    3.26         12
  trace backsolve             0.545   11.5%    3.06        168
  NPC reduce+recover          1.320   27.8%    1.00         12
  constraint location         0.185    3.9%    1.01        253
  border assembly             0.227    4.8%    1.00         66
  border dense solve          0.184    3.9%    1.00         12
  re-assembly                 0.198    4.2%    1.10         16
  other (remainder)           0.143    3.0%    1.00          0
  outside solve()             0.335    7.0%    1.00          3
    of which makeSolver       0.061    1.3%    1.00          1
    of which axis checks      0.096    2.0%    1.00          1
```

**Every leg that threads reads 3 to 6.8 cores and every leg that does not reads
1.00**, which is the column doing its job. At `OMP = MKL = 1` the same run is
6.775 s with **every leg at exactly 1.00** — the control that says the column is
not measuring noise — so eight cores buy **1.43x**.

**THE PER-SOLVE SPLIT AT 1 AND 8 THREADS**, which is how the flat leg was found
before the run-level budget existed. Last solve of three support sweeps, 2
Newton iterations, 14 border columns:

| leg | `OMP = MKL = 1` | `OMP = MKL = 8` | speedup |
|---|---|---|---|
| residual | 0.384 s | 0.086 s | **4.5x** |
| gradient | 0.123 s | 0.045 s | 2.7x |
| — of which `ComputeH` | 0.040 s | 0.028 s | 1.4x |
| trace factorisation | 0.105 s | 0.058 s | 1.8x |
| trace backsolve | 0.145 s | 0.081 s | 1.8x |
| **`NPCReduce` + `NPCRecover`** | **0.212 s** | **0.212 s** | **1.00x** |
| constraint location | 0.023 s | 0.023 s | 1.00x |
| border assembly | 0.046 s | 0.047 s | 0.98x |
| border dense solve | 0.034 s | 0.039 s | 0.87x |
| re-assembly | 0.002 s | 0.003 s | — |
| other (remainder) | 0.094 s | 0.097 s | 0.97x |
| **total** | **1.169 s** | **0.692 s** | 1.69x |
| whole run | 7.631 s | 5.146 s | 1.48x |

**Identical to the millisecond across a factor of eight in threads, and 30.6%
of the threaded step** — the largest leg there, ahead of the residual's 12.5%.
Serial legs together are 61% of a threaded step.

**THE LEG IS NEW AND `other` IS WHAT IT CAME OUT OF.** `TimedSolver` decorates
the TRACE solver, so it times the middle of
`NPCReduce -> trace solve -> NPCRecover` and neither end; the two element loops
landed in the remainder, which is why that remainder was **45% of a threaded
step and named nothing**. Timing `DarcyNPCSolver::ArrayMult()` and subtracting
the trace solve inside it leaves the traversal alone, and `other` falls from
0.357 s to **0.097 s**. Same `psi_ax` to every digit, 3.759851e-01, constraint
−6.466e-12: the change is a clock read.

**TWO INSTRUMENTS FOUND IT BEFORE THE LEG EXISTED, WHICH IS WHY THE LEG IS
BELIEVED.**

* **`perf` with DWARF call graphs** — `.eh_frame` unwinds a Release build with
  no `-g`. Of the main thread's 4.71 s, `NPCRecover` is 0.63 s and `NPCReduce`
  0.59 s: **26% together, with every sample of both on the main thread.**
  Beside them, also entirely serial: the trace solve 0.57 s,
  `DarcyForm::Reconstruct()` for `psi*` 0.60 s, assembly 0.30 s, the critical
  point searches 0.22 s. **`Reconstruct()` is the second item and is the same
  shape** — `ReconstructTotalFlux()` and all three `ReconstructFluxAndPot()`
  overloads are plain element loops with no `pragma omp` in their bodies,
  checked by extracting each whole function rather than a line range. It is
  harder than the traversal rather than easier: the scatter is disjoint, but
  the loop re-assembles integrators at the enriched order, so it needs
  `MFEM_THREAD_SAFE` and the reentrant `ElementTransformation` route. Worth
  about 1.10x against the traversal's 1.24x.
* **Column scaling** — sweeping `[boundary.exterior] Modes` at 10, 6, 3, which
  is 14, 10 and 7 border columns, the then-unnamed remainder read 0.330, 0.283
  and 0.232 s. A straight line in the column count, slope 0.014 s per column
  and intercept 0.134 s, so **59% of it scaled with the columns** — which is
  what said the traversal rather than the vector arithmetic around it.

**THE MECHANISM IS A MEASUREMENT OF SOMEBODY ELSE'S THAT DOES NOT CARRY.**
Neither routine has an `omp parallel` — single-vector or blocked — where every
other element-local loop in `DarcyHybridization` does. That is deliberate, and
upstream's doxygen on `NPCGradient` says why: the two loops are *"under 6% of
the step, flat in mesh size and order, and Amdahl caps any gain there"*, on the
pedestal problem at four resolutions. **Those four cases are fixed boundary and
ONE right-hand side.** A bordered step applies `J^-1` to `N + 4` columns against
one factorisation, so the traversal is `O( elements x columns )` where the
integrator-bound loops are `O( elements )` and already threaded. The share
inverts as soon as there is a border, and again with every exterior mode.

**WHAT IT WOULD BE WORTH, AND WHO CAN DO IT.** At the residual leg's own
efficiency on the same mesh in the same process — 4.5x — the traversal falls to
about 0.047 s: the step from 0.692 s to 0.527 s, **1.31x**, and the whole DIII-D
run about **1.24x**. `NPCRecover` is embarrassingly parallel as written, its
only writes being the calling element's own L2 dofs; `NPCReduce` accumulates at
face dofs shared by two elements and needs atomics or a face colouring. Both
are MFEM's, and the request is
`../mfem-hdg-dev/doc/HDG-NPC-TRAVERSAL-FROM-MEQ.md`.

**WHAT IS NOT THE ANSWER, MEASURED.** `OMP_NUM_THREADS=4` is **slower** than 8
— 5.8 s against 5.5 s — so the gap is not oversubscription of the 8 physical
cores by a 16-thread hyperthreaded topology. And `AssemblyMode`,
`LocalFactorMode` and `TraceAssemblyMode` are already at the settings M-99 and
M-101 chose; nothing in that group moves this leg, which reaches none of them.

### M-127

**`halfdisc.py --symmetric` MESHES `z >= 0` AND REFLECTS IT, AND THE PAIRING IS
EXACT RATHER THAN CLOSE.** M-125 established that `[solver] UpDownSymmetry` is
the identity on a symmetric mesh and refuses MAST-U's; this is the mesh it
stops refusing. The flag cuts the disc with a second half-plane, clips the
conductors and the limiter to the upper half analytically, meshes, and then
reflects.

**THE REFLECTION ACTS ON THE WRITTEN FILE, AND THAT IS WHY THE PAIRING IS
EXACT.** gmsh writes the half exactly as an ordinary run does and a pass then
rewrites `$Nodes` and `$Elements` in place. **The image coordinate is the
decimal TOKEN with its sign flipped** — drop a leading `-` or add one — which is
exact at whatever precision gmsh chose, where negating the double and
reformatting is exact only if the formatter is. `$MeshFormat` and
`$PhysicalNames` pass through untouched, so nothing about the format is
re-implemented.

| | |
|---|---|
| MAST-U's 23 conductors paired about `z = 0` | worst discrepancy **2.220e-16 m**, 1 self-mirror (the solenoid) |
| the symmetric mesh | 4753 nodes, 9396 triangles, 108 boundary segments, 43 seam nodes |
| unpaired nodes · involution failures · elements with no image · non-positive areas | **0 · 0 · 0 · 0** |
| worst node pairing | **0.000e+00 m**, at a tolerance of 4.808e-09 m |
| the committed `mastu-nke.msh`, same checker | 4735 nodes, **3588 unpaired**, 8620 of 9468 elements with no image |
| edge manifold | 14040 edges used twice, 108 once, the 108 exactly the tagged boundary |
| order 2, mid-edge offset from the chord | **3.212e-03 in each half, identically** |
| the default path, 4 option families | **byte-identical**; `diverted-tokamak-generated` 2854 elements, 1467 nodes |

3588 unpaired by VERTEX against the 3624 `buildMirrorMaps()` reports by DOF is a
different question with the same answer.

**THE SEAM IS SHARED AND NOT DUPLICATED, AND THE NEAR MISS IS INVISIBLE.** A
node at `z == 0.0` **exactly** keeps its tag and is used by both halves; a node
at `0 < |z| <= tol` is REFUSED, because duplicating it gives a mesh that looks
perfect and is two disconnected halves. Exact zero is available because the cut
is a straight geometric edge, the same property the axis relies on. Checked
topologically rather than by coordinates: every edge is used by exactly two
triangles or exactly one, and the ones used once are precisely the tagged
boundary segments.

**AND THE `z = 0` LINE GETS NO PHYSICAL GROUP.** It bounds the meshed half and
is INTERIOR to the finished file, so the boundary classification grows a third
case. Without it the seam joins Gamma — which would impose Gamma's transferred
exterior datum straight across the midplane, and is this module's own recorded
defect one geometry further on.

**THREE THINGS THAT COULD HAVE BEEN TABULATED AND ARE DERIVED INSTEAD.** The
orientation permutation comes from gmsh's own reference-element table rather
than from the two cases somebody thought to write down, so it is right at any
geometric order — checked at order 2 by the mid-edge offset above. The conductor
pairing is greedy AND CONSUMING, so two conductors at the same place pair with
each other rather than both claiming the first match, which is what makes the
table an INVOLUTION — the property `buildMirrorMaps()` itself asserts and the
one a partner table can get wrong while leaving nothing unmatched. And a
reflected conductor takes its PARTNER's attribute, which is not the attribute of
the half-plane copy.

**WHAT IT REFUSES**, each an `argparse` error at exit 2: an unpaired conductor,
a `--limiter` off the midplane (that flag names ONE circle and a mirror pair is
two), and a `--vessel` outline that is not its own image — tested on the CURVE
and not the vertex set, since a symmetric vertex set joined up asymmetrically
clips to an upper half whose reflection is a different outline.

**`--plasma` IS NOT REFUSED WHEN ASYMMETRIC, AND THE LINE IS DRAWN
DELIBERATELY.** It is a size field and not geometry: the machine is unchanged,
the meshed half is refined exactly as asked and the other half gets the
reflection. MAST-U's box IS asymmetric — `z` in [−1.5099, 1.4987] — so this is
live, and the tool prints a note rather than leaving it to be discovered.

**MEQ TAKES THE RESULT.** On the full symmetric mesh with
`[solver] UpDownSymmetry = true`, at `k = 1` capped at two Newton iterations,
`buildMirrorMaps()` built the potential, trace and flux maps, `projectUpDown()`
ran on the guess, and the bordered Newton iterated 1.053e-01 → 3.814e-02 →
2.825e-03 before the cap stopped it. The control on the committed asymmetric
mesh, same configuration, throws *"the potential space's element at
( 0.196750, −1.560556 ) with no mirror partner about z = 0"*.

**TWO DEFECTS FOUND THAT HAVE NOTHING TO DO WITH SYMMETRY.**

* **A `--vessel` reaching past Gamma silently meshes a bigger domain.**
  `occ.fragment` keeps the parts of a tool lying OUTSIDE the shape it is
  fragmented into, so a vessel polygon poking past `rho` becomes extra surfaces.
  Measured at `--rho 1.5` with a polygon out to `r = 2.5`: the written mesh
  reaches **r = 2.5**, with **32 nodes carrying Gamma's attribute** on edges
  inside it — the transferred exterior datum imposed in the middle of the
  domain. `--check` catches it, but only with `--check` and only after the mesh
  is built. It is now a parse-time refusal, which `--coil` and `--limiter`
  already had.
* **Eleven conductors collide with the limiter's attribute and twenty-one with
  the vessel's.** The conductors are `10 + i` while `LIMITER_ATTRIBUTE = 20` and
  `OUTSIDE_VESSEL_ATTRIBUTE = 30`, so the eleventh coil IS the limiter's
  attribute and the twenty-first is the outside of the vessel — and
  `[source] ExcludeAttributes = [ 30 ]` would then exclude a conductor. No
  example in the tree reaches it: MAST-U has 23 conductors and neither a
  limiter nor a vessel, and the limiter cases have two conductors. **Refused
  rather than renumbered**, since moving the constants would change every mesh
  that already has one and every config that names it, to buy a case nobody has
  posed.

**AND A CORRECTION TO WHERE THE REFUSAL HAPPENS.** `buildMirrorMaps()` is NOT
in `prepare()`. It is called from inside the bordered Newton driver, after the
mesh, the spaces, the assembly and the guess — so a run that is going to be
refused takes **about 55 s at 2322 elements and `k = 1`**, and printed nothing
at 150 s on the full mesh at `k = 2`. That is why `[solver] UpDownSymmetry` on a
file that also carries `[mesh.generate]` without `Symmetric = true` is now a
PARSE error: where the file says how the mesh is made, MEQ knows at parse time
that the answer will be no.

**WHAT IS NOT MEASURED, PLAINLY.** `k = 2` on the full symmetric mesh was not
run to convergence. The entity pairing is `k`-independent — elements match by
centroid, faces by centroid, and the seam faces are their own mirrors; only the
dof count INSIDE a matched entity changes, and those are symmetric iff the
entity is — so `k = 2` is a sound inference and not a measurement. And no
equilibrium was compared between the symmetric and committed MAST-U meshes:
they are different meshes of the same machine, so nothing bit-comparable exists
and only the reference would arbitrate.

### M-128

**MAST-U SPENT 93.6% OF ITS RUN EVALUATING THE INITIAL GUESS, AND THE LEG
PROFILE IS WHAT FOUND IT.** M-126 asked where a free-boundary run's cores go and
answered it on DIII-D. Asked of MAST-U, the same instrument gave a completely
different shape — and the difference was not the solver.

**Before**, `examples/mastu-nke.toml`, `OMP = MKL = 8`, `OMP_WAIT_POLICY=passive`,
quiet machine (loadavg 0.24 at the gate):

```
MEQ: wall 275.437 s = setup 0.038 + solve 273.924 + output 1.474
MEQ: where the run went -- wall 275.437 s, cpu 284.650 s, 1.03 cores on average
  re-assembly               133.518   48.5%    1.00         16
  outside solve()           133.447   48.4%    1.00          4
  NPC reduce+recover          2.403    0.9%    1.00         12
  trace backsolve             0.891    0.3%    3.32        168
  residual                    0.629    0.2%    7.20         32
```

**97% in two rows at 1.00 cores, and the leg M-126 filed upstream is 0.9% of
it.** `prepare()` cost **8.3 s a call here against 12 ms on DIII-D** — 670x for
1.87x the elements, so not a mesh-size effect.

**THE FLAT PROFILE NAMED IT IN ONE READING**, 56k samples over the whole
converging run:

| symbol | share |
|---|---|
| `meq::geometricKernel` | **50.23%** |
| Boost `ellint_rd_imp` (Carlson `R_D`) | **38.10%** |
| `meq::coilPsi` | **5.29%** |
| `meq::ellipsePsi` | 0.18% |

and the call graph put **all of it under `prepare()`** — 57% through
`GridFunction::ProjectCoefficient` onto the potential space and 43% through
`projectOntoTrace()`. What is projected is the `Type = "conductors"` guess: **23
`coilPsi` calls plus one `ellipsePsi` at every nodal point of BOTH spaces**,
about 95,000 points.

**THE ORDER IS THE WHOLE COST, AND IT IS A REFERENCE ORDER USED FOR A GUESS.**
`coilPsi` runs a panelled tensor Gauss rule; `defaultCoilQuadratureOrder = 32`
is chosen to make the interior usable as a reference. Measured over 1176 points
of MAST-U's disc against order 32, with `max |psi_coil| = 1.515779e-01 Wb/rad`,
all 23 conductors per point:

| order | max abs err | relative | us / point | vs 32 |
|---|---|---|---|---|
| 2 | 4.317102e-03 | 2.848e-02 | 6.6 | 211x |
| 3 | 1.544457e-03 | 1.019e-02 | 13.4 | 104x |
| 4 | 6.521601e-04 | 4.302e-03 | 23.2 | 60x |
| **6** | **9.055783e-05** | **5.974e-04** | **50.0** | **28x** |
| 8 | 1.168790e-05 | 7.711e-05 | 88.1 | 16x |
| 12 | 9.744078e-07 | 6.428e-06 | 197.0 | 7x |
| 16 | 4.775856e-08 | 3.151e-07 | 342.7 | 4x |
| 32 | — | — | **1392.4** | 1x |

**1392 us a point x 95,000 points x two projections is 265 s**, against the 267 s
the two rows measured — which is what says the mechanism is understood rather
than merely correlated.

**TWO CHANGES, AND NEITHER TOUCHES THE EQUILIBRIUM.**

* **`meq::guessCoilQuadratureOrder = 6`.** At 6 the guess is wrong by 9.1e-05
  Wb/rad against a `psi_ax` of 9.2e-02 — a thousandth of the quantity being
  guessed, and far below the 1.6e-04 m the X-point seed is already allowed.
  **The guess is not the problem statement**; nothing about the answer depends
  on it, only which basin Newton starts in.
* **`prepare()` caches the seed projection.** The driver prepares once to read
  `potential()` for its support sweep and `solve()` prepares again a moment
  later; both seed, with the same guess, and nothing between them changes it.
  The cached vectors are the projection's own output copied back, so this is
  bit-exact. It tightens one contract: a caller mutating a guess **in place**
  must call `setInitialGuess()` again. Both overloads drop the cache, and every
  caller in the tree already does exactly that.

| | wall | re-assembly | outside solve() |
|---|---|---|---|
| before | **275.437 s** | 133.518 s | 133.447 s |
| order 6 alone | **31.394 s** | 10.439 s | 10.534 s |
| order 6 + seed cache | **20.939 s** | **0.388 s** | 11.054 s |

**13.2x, and the answer does not move**: `psi_ax` 9.162567e-02, `psi_bnd`
2.864064e-02, the X-point at ( 0.598637, −1.097187 ) 1.578e-04 m from its seed,
the axis at ( 0.9489, −0.0000 ), 2 Newton iterations in the last sweep — every
printed digit as before. The constraint residual moves in its last, −1.990e-14
to −1.986e-14, which is the guess having changed. DIII-D is unmoved on all of
`psi_ax`, `psi_bnd`, the X-point and the constraint, warm-started and never
touching this path.

**WHAT IS LEFT IS ONE PROJECTION**, 11.054 s of 20.939 s, in `outside solve()` —
the driver's own `prepare()`, which is the first and cannot be cached from
anything. Lowering the order further, threading the projection, or evaluating
the conductors once per distinct point are the remaining levers and none has
been taken.

**THE TRANSFERABLE PART IS ABOUT THE ESTIMATE AND NOT THE NUMBER.** The comment
beside this guess read *"meq::ellipsePsi() is 66 us a point ... about 6 s of
MAST-U's 291 s"*. That was measured and it was right: `ellipsePsi` is 0.18% of
the run. What it did not say is that the same loop evaluates **23 conductors at
the same points**, and those are the other 93.4%. **A measurement of one term in
a loop is not a measurement of the loop.**

### M-129

**WHAT AN `mfem::Device` COSTS MEQ'S NON-LINEAR SOLVE, AND THE MECHANISM THE
DRIVER CLAIMED IS NOT THE ONE THE CROSS SUPPORTS.** `apps/meq.cpp`'s `--device`
comment asserted that the Newton step count rises *"because the device path's
element-local evaluation is inexact where `dF/dpsi` is non-zero"*, at a
*"1.7x to 2.0x whole-solve loss"*, with no anchor. Both halves are now
measured. RTX 2070 SUPER, CUDA 13.3, `MFEM_USE_DOUBLE`; **this card runs FP64 at
1/32 to 1/64 of FP32 where a datacentre part runs about 1/2, so the seconds are
this machine's and the ITERATION COUNTS are the transferable half.**

`OMP_NUM_THREADS = MKL_NUM_THREADS = 1`, `--device cpu` against `--device cuda`
in one binary:

| case | path | `dF/dpsi` | cpu | cuda |
|---|---|---|---|---|
| `soloviev-nstx`, k=3, 1536 el | plain | **0** | 1 it | 1 it |
| `mhd-rectangle`, k=2, 768 el | plain | non-zero | 5 it, 0.54 s | **9 it, 1.20 s** |
| `limited-tokamak`, k=3, 1601 el | **bordered** | non-zero | 12 it, 4.94 s | **12 it, 5.43 s** |
| `machine-f-diiid`, k=2, 4848 el | **bordered** | non-zero | 2 it, 7.60 s | **fails, falls back, 13 it, 8.82 s** |

**THE STATED MECHANISM IS FALSIFIED BY THE THIRD ROW.** `limited-tokamak` and
`mhd-rectangle` are both `Type = "mhd"` with non-zero `dF/dpsi`, and one pays
1.8x the steps while the other pays none. *"Where `dF/dpsi` is non-zero"* is
therefore not the discriminator. The Solov'ev row says nothing either way: one
step is one step.

**AND THE 1.7x TO 2.0x IS NOT A RANGE, IT IS A SPREAD ACROSS CASES** — 2.2x,
1.10x and 1.16x on the three that converge, with the third reaching a different
answer.

**THE DEVICE CAN CHANGE THE ANSWER AND CAN MAKE THE BORDERED NEWTON FAIL.**
`machine-f-diiid` under `--device cuda` does not converge on the bordered
Newton at all; the driver's reactive ladder catches it and
`bordered-picard-then-newton` finishes in 13 steps. **`psi_ax` reads
3.759873e-01 against the host's 3.759851e-01** — 2.2e-6 relative, in the sixth
digit — with the constraint residual 9.37e-14 against −6.466e-12. Twice
reproduced, to the same digits both times. **AND IT IS A DEFECT RATHER THAN A
PROPERTY: M-130 walks it to two missing syncs, after which this row reads 2 it
and the host's every digit.** The row stands because the failure it records is
what a device does to an unsynced read, which is the reason to run one. `limited-tokamak` by contrast agrees
exactly: 9.484400e-02 either way. So *"reproduces the CPU answer to every
printed digit"* is true of the cases it was measured on and **not** a property
of the device path.

**AT `OMP_NUM_THREADS=8` THE DEVICE STILL ABORTS, BUT NOT EVERYWHERE**, which
is the other half nobody had separated:

| | `OMP=8`, `--device cuda` |
|---|---|
| `mhd-rectangle`, plain | **aborts** — `alias pointer is not registered` / `host pointer is not registered` in `MemoryManager::CheckHostMemoryType_`, from several threads at once, landing as `terminate called recursively` |
| `limited-tokamak`, bordered | **converges**, 12 iterations, `psi_ax` 9.484400e-02, the host's answer |

`--device debug` tracks `--device cuda` on every row above.

**WHAT M-79 REPAIRED IS REPAIRED, AND IT IS A DIFFERENT SYMPTOM FROM THESE.**
The `0/0 Newton iterations on a case that needs four` and the 4.006e-02 flux
disagreement are gone: the same `OMP=1` configuration now runs a real Newton to
a residual floor of 2.14e-15 against the host's 7.03e-17. What remains is a
degraded rate, not a solve that does not run.

**THE TRANSFERABLE PART.** A claim with no anchor survived in the tree because
it was plausible and nobody could cheaply contradict it — and the commit that
introduced it recorded the OPPOSITE on its own case, *"no iteration penalty"*
on `limited-tokamak`, which this table reproduces. Both observations were
right; the sentence that generalised one of them into a mechanism was not.
**Two cases that disagree are a cross to run, not a discrepancy to pick a side
of.**

### M-130

**`machine-f-diiid` SOLVES ON THE DEVICE, AND M-129'S WORST ROW WAS TWO MISSING
SYNCS RATHER THAN A PROPERTY OF THE DEVICE PATH.** M-129 measured the bordered
Newton failing under `--device cuda`, the reactive ladder finishing in 13 steps
against 2, and `psi_ax` differing in its sixth digit. Every digit of that is
reproducible and none of it is the device: it is M-79's alias chain with two
links nobody had walked. Same machine, same binary, `OMP = MKL = 1`:

| case | host | cuda, M-129 | cuda, with both syncs |
|---|---|---|---|
| `soloviev-nstx`, k=3 | 1 it | 1 it | 1 it |
| `mhd-rectangle`, k=2 | 5 it | 9 it | **9 it — unchanged, and still open** |
| `limited-tokamak`, k=3 | 12 it | 12 it | 12 it |
| `machine-f-diiid`, k=2 | 2 it | fails, ladder, 13 it | **2 it** |

`machine-f-diiid` under `--device cuda` now reads `psi_ax = 3.759851e-01`,
`psi_bnd = 7.081394e-02`, X-point `( 1.200929, -0.999491 )` and constraint
residual `-6.466e-12` — **the host's answer to every printed digit** — and runs
to the same numbers under `--device debug`.

**THE CHASE, BECAUSE THE METHOD IS THE TRANSFERABLE PART.**
`mfem::Device( "debug" )` named the first fault on the first run, and each fix
moved it FORWARD, which is the only thing separating a fix from a coincidence:

| | where | what |
|---|---|---|
| 1 | `DarcyHybridization::EliminateVDofsInRHS`, `darcy_u = x.GetBlock(0)` | host memcpy out of an mprotected page |
| 2 | `solveWithNormalisation()`, `essentialTrace[ i ]` | raw `Array<int>::operator[]` on device state |
| 3 | `DarcyForm::ReconstructTotalFlux` in `postProcess()` | **MFEM's, still open**, and `debug`-only |

**Link 1 is the same gap M-79 records, running the other way.** `prepare()`
seeds the iterate by writing THROUGH `potentialGf` and `traceGf`, which are
`MakeRef` aliases of `solution`. A device write through an alias marks the
ALIAS device-valid, calls `AliasProtect` on the base's host range, and — by
`Memory::SyncAlias`'s own account — leaves the base's validity flags exactly as
they were. `formSystem()` then builds `darcySolution` and `traceX` from those
flags, so both views inherit the lie. Caught by breaking on `mprotect` with the
faulting page as the condition, which names the protector directly:

```
AliasProtect <- GetAliasDevicePtr <- Vector::operator= <- GridFunction::operator=
              <- meq::GradShafranovSolver::prepare(bool)
```

`solve()` already synced base to aliases at the end; nothing synced aliases to
base at the start. Three `SyncAliasMemory()` calls in `formSystem()`, inert
without a Device because `Memory::SyncAlias` returns on its first line when the
base is not registered.

**Link 2 is MEQ's own and needed no MFEM at all.** `FormLinearSystem()` has
just eliminated on `GetEssentialTrueDofs()`, so that `Array<int>` is device
state, and `Array<int>::operator[]` neither syncs nor invalidates. One
`HostRead()`.

**WHAT THE TWO SYMPTOMS SHARE.** Under `debug` link 1 is a named fault with a
backtrace. Under `cuda` nothing is protected, the stale host copy is used, the
border columns are differenced against a state that is not the iterate, and the
run reports *"the bordered Jacobian is singular in ( psi_ax, psi_bnd, a )"*.
**One fault, two faces**, and the CUDA face names a border and not a memory.

**IT IS NOT A REGRESSION FROM THE GUESS-SEED CACHE, WHICH IS WHERE THE
INVESTIGATION SPENT ITS TIME.** The `mprotect` backtrace lands in `prepare()`'s
cached branch — `potentialGf = guessSeedPotential`, 229,376 bytes, the potential
block — so the cache looked like the trigger and the immediately preceding
commit looked like the cause. Building the parent commit and running it under
`--device debug` killed that: it faults too, and harder, terminating where the
current tree recovers through the ladder. **The cached assignment is a device
write where the projection it replaces was a host write, so it CHANGES WHICH
SITE PROTECTS FIRST without being the defect.** A backtrace names the
instruction, not the bug.

**THE TRANSFERABLE PART, AND IT IS AGAINST M-129 RATHER THAN AGAINST ANYBODY
ELSE.** M-129 measured correctly and concluded *"a device is an instrument for
finding unsynced reads, not a second way to get the answer"* — while the thing
it had in hand WAS an unsynced read, of exactly the kind that sentence names.
The measurement was treated as a property to be recorded rather than a symptom
to be chased, and the difference between those two is one `gdb` session.
**A result that is stable, reproducible and in the sixth digit is still a bug
until something explains it.**

### M-131

**THE `q`-DRIVEN RUN NOW REPORTS THE `q` IT REACHED, AND IT IS THE TARGET TO
4.5e-07.** MEQ computed the safety factor, carried a `safety_factor` column for
it in `Output.cpp`, and wrote it on no run at all: the column is gated on
`FluxSurfaceFamily::safetyFactorAvailable`, which is
`static_cast<bool>( FluxFamilyOptions::toroidalField )`, and `apps/meq.cpp` set
that at neither extraction site. **Found by closing `INVERSION-PLAN.md` rather
than by using the code.**

**THE REFUSAL WAS RIGHT FOR ONE CASE AND WAS NEVER REVISITED FOR THE OTHER.** A
`meq::Source` carries `g g'` and not `g`, so a PRESCRIBED-field run cannot
recover `g = sqrt( g_edge² + 2∫ g g' dΨ )` without a constant of integration no
`[source]` key supplies — and a column of zeroes would be indistinguishable from
a machine with no toroidal field, so absent is the correct signal and stays.
**`[source] SafetyFactorFile` is the exception**: the outer Newton solves for
the coefficients of `g²` itself, and `apps/meq.cpp` already held them in scope
at the writer and already printed them to the `.nc` as an attribute. So the one
run whose entire purpose is to reach a target `q` was the one run that could not
show the `q` it reached.

`examples/q-driven.toml`, 512 elements at `k = 2`, degree-2 `g²`, 24 surfaces
over `Ψ_N ∈ [0.05, 0.95]`, one run, 44 s:

| `Ψ_N` | `Ψ` source | `q` reached | `q` asked for | relative |
|---|---|---|---|---|
| 0.0500 | 0.9500 | 1.072306 | 1.072307 | 4.769e-08 |
| 0.1255 | 0.8745 | 1.063448 | 1.063448 | 7.777e-08 |
| 0.2351 | 0.7649 | 1.058562 | 1.058562 | 3.130e-08 |
| 0.3788 | 0.6212 | 1.065494 | 1.065494 | 3.516e-08 |
| 0.5567 | 0.4433 | 1.097155 | 1.097156 | 7.568e-08 |
| 0.7686 | 0.2314 | 1.188384 | 1.188384 | 1.714e-07 |
| 0.9500 | 0.0500 | 1.432489 | 1.432488 | 1.586e-07 |

**Worst 4.547e-07 over all 24 surfaces**, against the same run's coefficient
recovery of 5.037e-06 — the `q` reads an order better than the `g` it is built
from because the target was measured on this mesh, so the extraction error is
common to both sides.

**THE HAZARD IS A REFLECTION AND IT WOULD HAVE BEEN INVISIBLE.**
`extractFluxSurfaces` calls the callback with the **raw** `ψ`, and the
coefficients are in the **source's** `Ψ` — `fitToroidalFieldSquared()` fits
against `1 − normalisedFlux` for exactly this reason, because `meq::
normalisedFlux` is zero on the axis and the source's `Ψ` is one there. A wiring
that reached for the family's own label would have evaluated the polynomial
reversed and reported a `q` that is smooth, plausible and the profile backwards.
`examples/q-driven-q.dat`'s own header already records the same trap on the way
IN: *"Handing a table written the other way round does NOT fail: it converges,
at full order, to an equilibrium with its shear reversed"*. **The same
convention pair, the same failure mode, met twice at opposite ends of one
loop.**

**TWO GUARDS, AND NEITHER IS DECORATION.** `g²` is a *fitted* polynomial, so a
degree the family does not determine can dip below zero inside the cut while
every coefficient looks ordinary; `sqrt` of that is a NaN in a file whose only
per-node mask is about the band. Positivity is therefore checked over the
family's own range before the callback is installed, and a failure leaves the
variable **absent** — the same signal a prescribed run gives — and says so on
stdout. And the acceptance asserts the column is **not flat**, because a lambda
that ignored its argument would write one repeated value and pass a tolerance
against a slowly varying target.

`DriverAcceptance::theQDrivenRunReportsTheSafetyFactorItReached`. It shares one
run with `theDriverSolvesForTheToroidalField`: the two are halves of one round
trip, the run is 44 s, and running it twice would make the second case's numbers
those of a different equilibrium.

### M-132

**THE INITIAL `psi_bnd` SELECTS AMONG DISCRETE EQUILIBRIA ON THE DIVERTED
MACHINE, AND REFINEMENT DOES NOT MERGE THEM.** `FREE-BOUNDARY-PLAN.md` §10.5
measured two answers from two starting values and nominated *"whether the gap
closes under refinement"* as the next thing to measure — FB-R. It is now
measured, on the **library** fixture rather than through the driver, so that the
located axis, the support sweeps and the guess transfer are not variables.

`tests/convergence/DivertedMachine.hpp` at three uniform refinements of
`examples/diverted-tokamak.msh`, `k = 2`, the X-point pinned at freegs4e's
`( 1.093144, -0.603965 )`, everything else identical. The seed is
`GradShafranovSolver::setBoundaryFluxInitial()`, which exists because `psi_bnd`
had no setter — it starts at `psiBoundary()`, zero on a fresh solver.

| elements | seed | `psi_ax` | `psi_bnd` | its | final `‖r‖` |
|---|---|---|---|---|---|
| 2642 | 0 | 8.48629920e-02 | 3.22975097e-02 | 26 | 5.53e-12 |
| 2642 | 3.240e-02 | 8.26600366e-02 | 3.23793176e-02 | 7 | 3.95e-12 |
| 10701 | 0 | 7.79347348e-02 | **4.63459237e-03** | 11 | 1.63e-14 |
| 10701 | 3.240e-02 | 8.25753259e-02 | 3.23486098e-02 | 7 | 2.15e-12 |
| 43026 | 0 | 7.56865206e-02 | **4.33770950e-03** | 12 | 2.88e-14 |
| 43026 | 3.240e-02 | **1.29008101e-01** | 2.71353946e-02 | **97** | 2.21e-14 |

**EVERY ROW IS A CONVERGED SOLVE** — the worst final residual in the table is
5.5e-12 and four of six are below 3e-12. None of this is a failure to converge.

**THE ARMS DO NOT MERGE.** At 43,026 elements they stand at `psi_bnd`
2.714e-02 against 4.338e-03, a factor of **6.3**, having been 0.25% apart on
the coarsest mesh. Whatever else is true, *the gap does not close*, which is the
question FB-R asked.

**AND NEITHER ARM IS A CLEAN SEQUENCE, WHICH IS WHY THIS IS NOT THE WHOLE
ANSWER.** The unseeded arm is the coarse mesh that is the outlier — 3.230e-02
once, then 4.635e-03 and 4.338e-03, settling from the middle level. The seeded
arm is the opposite: stable to **0.09%** across the first two levels and then
**breaking** on the third, at 97 iterations and a `psi_ax` 56% above the
reference. Two refinements are not enough to call either a limit.

**THE AXIS COLUMN WAS DROPPED AND THE REASON IS WORTH MORE THAN THE COLUMN.**
The probe reported a magnetic axis from
`CriticalPointFinder::findAxis( AxisSense::Maximum )`, and on levels 0 and 1 it
returned `( 1.006240, -1.099327 )` for every row. That is **inside coil P1L** —
`CentreR = 1.00`, `CentreZ = -1.10`, half-extent 0.05, so `r ∈ [ 0.95, 1.05 ]`
and `z ∈ [ -1.15, -1.05 ]`. `findAxis()` with no argument refuses outright here,
reporting *"3 maxima, 2 minima and 2 saddles"*, and forcing a sense does not
find the plasma — it picks a conductor's O-point.

**This is §10.7's defect 2 arriving from the library side.** `apps/meq.cpp` has
a three-tier conductor exclusion on the fill seed for exactly this reason; the
library's `findAxis()` has none, because it is a general critical-point finder
and a coil O-point is a real critical point of `psi`. So **any axis position
taken from `findAxis()` on a machine with meshed conductors is about a
conductor unless something excludes them**, and a displacement quoted in metres
from such a call — §10.5's own `6.9e-02 m` among them — needs checking against
the coil list before it is believed.

**WHAT FB-R NOW NEEDS**, stated so the next attempt does not repeat this one: a
conductor-excluded axis location, a topology diagnostic that says which
equilibrium each row is, and a fourth level. What it has established is that the
phenomenon is **not** an artefact of the driver's guess handling — it reproduces
in a bare library fixture — and that `setBoundaryFluxInitial()` must stay a
library setter and never become a TOML key, since a file that silently moved
`psi_bnd`'s start by 3e-02 would move the reported `psi_bnd` by a factor of six.

### M-133

**AN X-POINT ON THE PLASMA EDGE COSTS THE APPROXIMATION NOTHING, AND THE
SUPPORT-CORNER COUNTER-ARGUMENT DOES NOT BITE.** `FREE-BOUNDARY-PLAN.md` §10.1
predicts this and says plainly that it cuts against expectation: at a null `Ψ`
vanishes **quadratically**, so a source `F ~ Ψ^j` vanishes to order `2j` there
and the crossing is *smoother* than the branches rather than rougher. §10.1
nominated XP-1 to stop it being a prediction; XP-1 went elsewhere — to the fill
leaking through the band of elements straddling the separatrix — so FB-S has
been an argument and not a number for the whole campaign.

**THE DESIGN NEEDS NO SOLVER AND NO EXACT DIVERTED EQUILIBRIUM**, which is why
it could have been done at any point in the last month.
`theCutCapsTheOrderBeforeAnyMethodIsChosen` already measures the plasma edge's
cap as a property of **best approximation** — element-local L2 projection onto
`P_k` over a dyadic mesh sequence, quadrature at `2k + 14` so the rule is never
what is being read. FB-S is that case with **the cut as the only thing that
differs**:

| arm | `φ` | zero set |
|---|---|---|
| smooth | `a² − dr² − dz²` | a circle |
| crossed | `dz² − dr²` | two lines through a null |

with `ψ = (φ/a²)^m`, `m = j + 2`, over `n = 8, 16, 32, 64` and `k = 1, 2, 3`.
The rate is taken across the whole sequence rather than per pair, for
`ExtensionConvergence`'s reason: which elements the cut passes through is not a
smooth function of `h`.

```
  j = 0 ( m = 2, cap m + 1/2 = 2.5 )
    k       smooth     rate      crossed     rate   difference
    1    2.709e-04    1.918    6.081e-03    1.988     +0.070
    2    2.333e-05    2.484    8.069e-05    2.771     +0.287
    3    7.583e-06    2.532    2.424e-05    2.513     -0.020

  j = 1 ( m = 3, cap m + 1/2 = 3.5 )
    1    2.631e-04    1.877    6.463e-02    1.963     +0.086
    2    9.032e-06    2.943    8.194e-04    2.977     +0.033
    3    7.166e-07    3.388    7.113e-06    3.908     +0.520

  j = 2 ( m = 4, cap m + 1/2 = 4.5 )
    1    2.813e-04    1.825    5.947e-01    1.924     +0.099
    2    9.629e-06    2.995    1.082e-02    2.940     -0.055
    3    3.254e-07    3.795    1.227e-04    3.960     +0.165
```

**THE WORST DROP OVER NINE `( j, k )` PAIRS IS 0.055 OF AN ORDER**, against the
0.25 the case allows for a moving cut, and the crossed arm is **faster** in
seven of the nine. So the prediction holds, and so does the reason for caring
about it: the competing argument — that the support's **corner** costs an order
the way `meq-corner-costs-the-extension-its-order` measures it costing the
extension — would have capped a diverted plasma below FB-4's `|d|^{j+2}`
wherever the null sits inside an element. It does not.

**ONLY THE RATES ARE COMPARABLE AND THE TEST SAYS SO IN ITS ASSERTION.** The
smooth arm's `φ` is bounded by `a²` and the crossed arm's grows to the corner of
the box, so the two carry constants up to **2000×** apart — 2.813e-04 against
5.947e-01 at `j = 2, k = 1`. A reader comparing the two error columns is reading
the normalisation and not the cut, which is why `worstDrop` is a difference of
rates.

**WHAT IT DOES NOT SETTLE.** This is best approximation, so it bounds what any
method on these spaces can do and says nothing about whether MEQ's solve
attains it on a diverted machine — the fill, the moving support and the X-point
border are all outside it. It is the half of §10.1 that was resting on an
argument; the solved half is XP-3's and is already green.

The case is `anXPointOnThePlasmaEdgeCostsTheApproximationNothing` in
`tests/convergence/PlasmaEdgeConvergence.cpp`.

### M-134

**A PRESCRIBED CURRENT MOVES CONFINEMENT OFF ZERO AND DOES NOT CURE IT, SO THE
FIXTURE IS NOT REPAIRABLE THE WAY §7.12b's WAS.** `FREE-BOUNDARY-PLAN.md` §11.2
retracts `theTwoBordersConvergeTogether`'s `ψ_ax` and `ψ_bnd` as physics because
the case closes only with a source that carries current in the **vacuum**, and
records two repairs each tried alone and each failing 4 of 4. FB-T is the
missing fourth cell — confinement **and** a prescribed current — and the
attribution turns on it: if pinning the amplitude repairs confinement, the
failure was §7.14's amplitude-fixed-with-a-moving-support difficulty, whose
known cure is exactly the current border.

**THE WHOLE 2 × 2 IS RE-RUN RATHER THAN ONE CELL BOLTED ONTO REMEMBERED
NUMBERS**, which is M-82's lesson in this same file: a one-key experiment
separates two hypotheses only if everything else is where you think it is, and
there the fill was off in *both* arms of every experiment run on it. `k = 2`,
`n = 24` on the half-disc, four limiter radii, `PowerProfile` at `j = 1`.

```
  arm                          R=1.05  R=1.15  R=1.20  R=1.30   converged
  unconfined, no current        ok5     ok8     ok6     ok5      4/4   control
  CONFINED, no current          --      --      --      --       0/4
  CONFINED, mu0 Ip = 0.02       --      --      --      --       0/4
  CONFINED, mu0 Ip = 0.05       --      --      --      --       0/4
  CONFINED, mu0 Ip = 0.10       --      --      --      --       0/4
  CONFINED, mu0 Ip = 0.20       ok20    --      --      --       1/4
  CONFINED, mu0 Ip = 0.50       ok25    --      ok14    --       2/4
  CONFINED, mu0 Ip = 1.00       ok40    --      --      --       1/4
  CONFINED, mu0 Ip = 2.00       --      ok19    ok32    --       2/4
  CONFINED, mu0 Ip = 5.00       --      --      --      --       0/4
```

**THE ANSWER IS PARTIAL SUPPORT AND AN EXPLICIT REFUSAL OF THE REST.** A
prescribed current does move confinement off zero — 0 of 4 becomes 2 of 4 — so
§7.14's diagnosis is doing real work and the amplitude is part of the
difficulty. It is **not the whole difficulty**, and three things in that table
say so:

* **The clamped arm is cured outright by the same lever.** §11.2's own row reads
  *clamped profiles with a prescribed current converge in 22 steps*, 4 of 4.
  Confinement gets 2 of 4 at its best over a **250× range** of target current.
* **The maximum is interior, so this is not a truncated sweep.** 0.50 and 2.00
  both give 2 of 4 and 5.00 gives 0 of 4. The first run of this case peaked at
  its own largest target and would have been published as "2 of 4" when the
  sweep had simply stopped too early; the case now prints a warning when the
  best cell is the largest target tried, and on this run it does not fire.
* **WHICH radii close depends on the current, and not monotonically.** 0.50
  closes `R = 1.05` and `1.20`; 2.00 closes `1.15` and `1.20`; 1.00, between
  them, closes only `1.05`. Over the whole sweep three of the four radii close
  at *some* current and `R = 1.30` never does, but **no single current closes
  more than two** — 6 of 32 confined cells in all. A lever that cured the
  underlying ill-posedness would not behave like this; a lever that shifts which
  discrete branch a given radius lands on would.

**SO THE FIXTURE IS NOT REPAIRABLE BY GIVING IT THE INGREDIENT IT LACKS**, and
that is the opposite of how §7.12b's fixture and M-82's fill came out. §11.2's
retraction stands unchanged, and what FB-T adds is that it is not going to be
lifted by the current border.

**THE ONLY ASSERTION IS THE CONTROL.** `baseline == 4` is what
`theTwoBordersConvergeTogether` already asserts, and it is in this case so that a
zero elsewhere in the cross is a statement about the key that moved rather than
about the fixture having drifted underneath it. Asserting on the confined arms
would be asserting a known defect, which the testing stance forbids.

The case is `confinementWithAPrescribedCurrentIsTheFourthCell` in
`tests/convergence/FreeBoundaryCoupling.cpp`.

### M-135

**THE THREADED NPC TRAVERSAL IS WORTH 1.25× ON THE WHOLE DIII-D RUN, AND IT IS
THE ONLY THING THAT MOVED.** [M-126](MEASUREMENTS.md#m-126) measured
`NPCReduce()`/`NPCRecover()` at 0.212 s whether given one thread or eight and
called it MEQ's largest single serial leg at 30.6% of a threaded step; upstream
threaded both overloads, and `../mfem/install` is rebuilt from `meq-integration`
at `3424f33bd5` to take it.

**THE PAIR IS INTERLEAVED AND THE MACHINE WAS QUIET**, which both matter here
and neither is decoration. A binary built before the upgrade was kept, and
`libmfem.a` is static, so the two arms differ in the library and in **nothing
else** — same MEQ commit, same mesh, same configuration. Runs alternate
`before, after, before, after`: a batched pair is what caught upstream's own
`OMP_WAIT_POLICY` figure, and this machine had three agents on it that day.
Three pairs per thread count, medians below, `machine-f-diiid.toml`,
`OMP_WAIT_POLICY=passive`.

**`OMP = MKL = 8`:**

| leg | before | cores | after | cores | |
|---|---|---|---|---|---|
| **NPC reduce+recover** | **1.315** | **1.00** | **0.503** | **3.32** | **2.61×** |
| residual | 0.301 | 6.76 | 0.301 | 6.74 | flat |
| gradient | 0.244 | 4.71 | 0.250 | 4.67 | flat |
| of which ComputeH | 0.148 | — | 0.152 | — | flat |
| trace factorisation | 0.224 | 3.29 | 0.225 | 3.31 | flat |
| trace backsolve | 0.490 | 3.18 | 0.494 | 3.20 | flat |
| border assembly | 0.113 | 2.53 | 0.110 | 2.53 | flat |
| of which driver prepare | 0.182 | 1.01 | 0.146 | **1.30** | the threaded assembly |
| of which postProcess | 0.616 | 1.00 | 0.620 | 1.00 | flat, and see below |
| **solve** | **3.316** | 2.42 | **2.471** | 3.39 | **1.34×** |
| **wall** | **4.148** | | **3.310** | | **1.25×** |

**TWO THINGS MOVED, NOT ONE.** The traversal is the headline; the other is
`driver prepare`, 0.182 → 0.146 s with its `cores` going 1.01 → 1.30, which is
upstream's threading of `DarcyForm`'s hybridized assembly loops arriving in MEQ.
It is small here and it is the reason MEQ has to make the thread-safety promise
at all — those loops evaluate the **linear** domain integrators, so the promise
is wider than the one MEQ asked for.

**AND THE ONE-THREAD TAX IS REAL, WHICH IS THE HALF A SPEEDUP TABLE HIDES.**
At `OMP = MKL = 1` the same pair reads:

| | before | after | |
|---|---|---|---|
| NPC reduce+recover | 1.375 | 1.446 | **0.95×** |
| wall | 6.983 | 7.121 | 0.98× |

**About 5% slower on the leg and 2% on the run**, from the OpenMP region and
the colouring that a single thread pays for and cannot use. Upstream's own
table shows the same sign at the same size, 0.199 → 0.218 s at one thread. So
the threaded traversal is **not free below two threads**, and a serial
deployment of MEQ is very slightly worse off for this change.

**`psi_ax = 3.759851e-01` IN ALL TWELVE RUNS**, with the X-point at
`( 1.200929, −0.999491 )` and the constraint at `−6.466e-12`. Upstream warned to
expect the last bits to move, since the blocked traversal reassociates; on this
problem it did not cost even those.

**UPSTREAM PREDICTED 1.24× FROM MEQ'S OWN ARITHMETIC AND THE ANSWER IS 1.25×.**
That is the part worth keeping: the prediction was made by taking M-126's table,
applying a measured 4.35× to one leg and leaving every other leg alone. The leg
factor here is 2.61× rather than 4.35× — MEQ's traversal is twelve calls of
fourteen columns and upstream's fixture is not — and the *run* figure still
lands where the arithmetic said, because the share was right even though the
factor was not.

**WHAT IS NOW THE LARGEST SERIAL ITEM.** `postProcess` — `DarcyForm::Reconstruct()`,
four integrators re-assembled at the enriched order per element — is **0.620 s
of a 3.302 s run, 18.8%, at 1.00 cores**, and it did not move. With the
traversal threaded it is the biggest single-threaded thing left in a MEQ run
after the direct trace solve.

### M-136

**`LocalFactorMode::Batched` NO LONGER COSTS THE CONDENSATION CACHE, SO
[M-101](MEASUREMENTS.md#m-101)'s TRADE IS GONE.** M-101 read that key as *"a
trade and not a knob that does nothing — it buys a batched local factorisation
and pays the condensation cache"*, because `CanCacheCondensation()` refused
outright under it. Upstream has since given `FactorElementsBatched()` the
cache's own buffer as its destination and done the same for the batched
face-pair kernel, and told MEQ that `condensationCacheTaken()` should now report
`yes` on that row. **It does.**

`machine-f-diiid.toml`, the new library, quiet machine, `OMP = MKL = 8`,
`OMP_WAIT_POLICY=passive`, three interleaved pairs, medians:

| `LocalFactorMode` | wall | solve | gradient | ComputeH | cache taken |
|---|---|---|---|---|---|
| `serial` — the default | 3.454 | 2.579 | 0.264 | 0.159 | **yes** |
| `batched` | 3.418 | 2.578 | 0.256 | **0.152** | **yes** |

**IT IS NOW A SMALL CONSISTENT WIN WHERE IT WAS A WASH, AND IT IS NOT WORTH A
DEFAULT.** `ComputeH` is lower in **all three** pairs — 0.164/0.156,
0.159/0.152, 0.157/0.152 — for about 4%, and the gradient leg follows it at 3%.
The wall moves 1.0%, which is **inside the serial arm's own spread** of 3.431 to
3.511 s. So the honest statement is that the batched local factorisation now
buys 4% of a leg that is 6% of the solve, and nothing a user would see.

**MEQ'S DEFAULT STAYS `Serial`**, on the ground that a 1% wall change is not a
reason to move a default that has never been asserted bit-for-bit against the
other route — which is the distinction `AssemblyMode` and `TraceSolver` are held
to under *the two performance keys* in `CLAUDE.md`. What has changed is the
**reason**: M-101 said the key was a trade, and it is now simply a small gain
that MEQ does not take.

### M-137

**THE A+B+C PAYOFF: `1.16×` ON THE WHOLE DIII-D RUN AT EIGHT THREADS, AGAINST A
PREDICTION OF `1.14×`.** `THREADING-PLAN.md`'s items A, B and C were built,
merged and asserted bit-exact, and the one thing the plan never had was the
clock. This is it.

**BOTH ARMS ARE BUILT AGAINST THE SHIPPING `libmfem.a`, WHICH IS THE DESIGN
DECISION THAT MAKES THIS A MEASUREMENT ABOUT MEQ.** `before` is `60d8505` — the
commit before `0108f49` landed A, B, C and D — plus **one line** cherry-picked
from `c3af7e5`, `SetIntegratorsThreadSafe()`. Without that line the baseline
exits 2 against today's library and there is no arm at all; it is a promise the
library reads before deciding whether a threaded loop is permitted, and it
changes no arithmetic. `after` is the working tree. `libmfem.a` is static and
both were linked after the install, so **MEQ's own source is the only
variable** — the mirror image of [M-135](MEASUREMENTS.md#m-135), where the
library was the only variable.

`machine-f-diiid.toml`, `OMP_WAIT_POLICY=passive`, quiet machine verified by the
harness itself rather than on a claim, interleaved `before, after, before,
after`:

| | pair ratios | median | best-of | arm medians |
|---|---|---|---|---|
| **`OMP = MKL = 8`**, six pairs | 1.133 1.136 1.143 1.170 1.196 1.200 | **1.156×** | 1.134× | 1.174× |
| **`OMP = MKL = 1`**, three pairs | 1.020 1.032 1.033 | **1.032×** | 1.023× | 1.028× |

The `before` arm's own spread is 1.065× and the `after` arm's is 1.032×, which
is why the pair ratio rather than a ratio of medians is the statistic quoted:
both arms drift together across a session and an interleaved pair cancels it.
**`1.16×` against the plan's `1.14×`, and the `1.106×` floor it quoted for the
measured parts alone is cleared.**

**WHERE IT CAME FROM, LEG BY LEG, WHICH IS THE PART AN ESTIMATE CANNOT GIVE.**
One pair, both arms, same run:

| leg | before | cores | after | cores | | item |
|---|---|---|---|---|---|---|
| **border dense solve** | 0.205 | 1.00 | **0.005** | 1.00 | **41×** | A |
| **border assembly** | 0.247 | 1.00 | **0.111** | 2.60 | **2.2×** | A **and** B |
| **constraint location** | 0.190 | 1.01 | **0.057** | 4.04 | **3.3×** | B, C |
| — of which cold full sweep | 0.076 | 1.00 | 0.014 | **7.79** | 5.4× | C |
| — of which `I_p` | 0.069 | 1.03 | 0.017 | 5.30 | 4.1× | B |
| — of which axis | 0.083 | 1.00 | 0.022 | 5.52 | 3.8× | C |
| other (remainder) | 0.154 | 1.00 | 0.116 | 1.06 | 1.3× | — |
| outside solve() | 0.317 | 1.14 | 0.239 | 1.67 | 1.3× | D instruments it |
| residual | 0.332 | 6.76 | 0.318 | 6.75 | flat | — |
| gradient | 0.278 | 4.66 | 0.270 | 4.75 | flat | — |
| trace factorisation | 0.256 | 3.38 | 0.238 | 3.31 | flat | — |
| trace backsolve | 0.554 | 3.17 | 0.523 | 3.22 | flat | — |
| NPC reduce+recover | 0.548 | 3.39 | 0.510 | 3.41 | flat | — |
| postProcess | 0.628 | 1.00 | 0.632 | 1.00 | flat | — |
| **wall** | **4.114** | 2.49 | **3.427** | 2.84 | **1.20×** | |

**WHAT CAN AND CANNOT BE ATTRIBUTED TO ONE ITEM, AND THE LEG BOUNDARIES DECIDE
IT.** `border dense solve` is **item A alone**: 0.205 → 0.005 s, **0.200 s and
4.9% of the before run**, from 5.0% of the budget to 0.1%. That is the sparsity
argument landing exactly — 984 nonzeros in a vector of about 109,000, and the
dot products no longer walk the zeros. A dot product that skips 108,000 zeros
does not get cheaper, it stops existing, and "~0.2%" was the wrong shape of
estimate for it. `constraint location` is **B and C**: 0.133 s, 3.2%, against
the 6.2% the two were sized at together.

**`border assembly` CANNOT BE SPLIT BY THIS RUN AND MUST NOT BE READ AS A's.**
That `LegTimer` is installed at six sites — `exteriorTransmissionRows()`, which
is item A's, and `assembleCurrentColumn()`, `assembleCurrentNormalisationCorner()`
and `assembleCurrentRow()`, which are item B's threaded plasma loops — so its
0.247 → 0.111 s carries both, and its **2.60 cores is B**, not a sparsity item
acquiring parallelism. Separating them needs a third binary with one item in it.

**THE LEGS NO ITEM TOUCHED ARE THE CONTROL.** residual, gradient, trace
factorisation, trace backsolve and NPC reduce+recover are flat to 2–6%, one
pair's noise, and `postProcess` is flat to **0.6%** — which is what says the
difference is where the items are and not spread over the run.

**THE `cores` COLUMN IS THE INDEPENDENT WITNESS.** Every leg an item threaded
reads 1.00 in the `before` arm and 4.04 to 7.79 in the `after` one, and every
leg no item touched reads the same in both. `cold full sweep` at **7.79 cores of
8** is item C's critical-point sweep doing what it was written to do.

**THE ANSWER SURVIVES IT, AND THE EXACTNESS CLAIM NEEDS ONE CORRECTION.**
`psi_ax = 3.759851e-01` in all eighteen runs, with `psi_bnd`, the X-point
`( 1.200929, −0.999491 )` and the constraint `−6.466e-12` identical throughout.
The `psi.gf` md5 is constant within each arm at each thread count, so the runs
are deterministic. **It is not constant ACROSS the arms**: at `OMP = MKL = 1`,
`ψ` moves by `2.6e-15` absolute against a `max |ψ|` of `3.76e-01`, i.e. **6.9e-15
relative**, and exactly one printed line moves with it — the iteration-2 Newton
residual, `2.992906e-12` against `2.992877e-12`, itself at the round-off floor.
Item B's per-element partial sums are a different grouping of the same
quadrature contributions at any thread count, which is what that is.

**BUT `OMP_NUM_THREADS` ITSELF STILL CHANGES NO PRINTED DIGIT, AND THE
SEPARATION MATTERS.** Moving the thread count from 1 to 8 with **`MKL_NUM_THREADS`
held at 1** gives byte-identical printed output; moving both together moves that
same residual line. So the cross-cutting rule the plan states is intact and the
mover is MKL's blocked BLAS-3 reassociating, which `CLAUDE_HDGGS.md` already
records at 1.3e-15 and which the registered tests pin at `MKL_NUM_THREADS=1`.
**A one-axis claim needs a one-axis experiment**: the first reading moved both
and would have blamed OpenMP for MKL's arithmetic.

**AND THE LIBRARY CHANGE IS INVISIBLE HERE, WHICH STRENGTHENS M-135.** The
preserved pre-upgrade binary and today's give the same printed output to every
digit at both thread counts, including that residual line — so M-135's
"`psi_ax` did not move" is the weaker statement, and the true one is that the
whole printed run did not move.

### M-138

**THE WHOLE-RUN BUDGET ON THE SHIPPING LIBRARY, AND IT IS NOT THE ONE
`THREADING-PLAN.md` WAS WRITTEN AGAINST.** Same run as
[M-137](MEASUREMENTS.md#m-137)'s `after` arm — `machine-f-diiid.toml`,
`OMP = MKL = 8`, `OMP_WAIT_POLICY=passive`, quiet machine. This supersedes
[M-126](MEASUREMENTS.md#m-126)'s budget as the sizing table for anything on this
path; M-126 stays as the measurement that found the flat NPC leg.

```
MEQ: where the run went -- wall 3.427 s, cpu 9.745 s, 2.84 cores on average
                              wall s   share   cores      calls
  setup                       0.020    0.6%    1.00          1
  solve                       2.577   75.2%    3.45          3
  output                      0.830   24.2%    1.00          1
    of which postProcess      0.632   18.5%    1.00          1
  residual                    0.318    9.3%    6.75         28
  gradient                    0.270    7.9%    4.75         12
    of which ComputeH         0.160    4.7%       -         12
  trace factorisation         0.238    6.9%    3.31         12
  trace backsolve             0.523   15.2%    3.22        168
  NPC reduce+recover          0.510   14.9%    3.41         12
  constraint location         0.057    1.7%    4.04        253
    of which axis             0.022    0.6%    5.52         19
      cold full sweep         0.014    0.4%    7.79          1
    of which X-point          0.017    0.5%    1.00         16
    of which I_p              0.017    0.5%    5.30         28
  border assembly             0.111    3.2%    2.60         66
    of which transmission     0.033    1.0%    1.00          3
  border dense solve          0.005    0.1%    1.00         12
  re-assembly                 0.191    5.6%    1.08         16
  other (remainder)           0.116    3.4%    1.06          0
  outside solve()             0.239    7.0%    1.67          3
    of which makeSolver       0.055    1.6%    1.00          1
    of which axis checks      0.020    0.6%    6.91          1
    of which driver prepare   0.159    4.7%    1.26          1
```

**THE 30/70 SPLIT IS INVERTED: ABOUT 59% OF THE RUN NOW THREADS AND 41% DOES
NOT.** The seven legs reading 2.60 cores or better — residual, gradient, trace
factorisation, trace backsolve, NPC reduce+recover, constraint location, border
assembly — are **2.027 s of 3.427, 59.1%**. Three things moved it there and only
one of them is MEQ's: upstream threaded the NPC traversal
([M-135](MEASUREMENTS.md#m-135)) and the hybridized assembly loops, and MEQ's
own A, B and C took `constraint location` and `border assembly`
([M-137](MEASUREMENTS.md#m-137)).

**THE LARGEST SINGLE-THREADED ITEM IS `Reconstruct()` AND IT IS MFEM'S.**
`postProcess` is **0.632 s, 18.5%, at 1.00 cores** — four integrators
re-assembled at the enriched order per element, the price of reporting `ψ*`
rather than `ψ_h`. It did not move under any of this and it is now larger than
the NPC traversal it used to sit behind.

**THE LARGEST SINGLE-THREADED ITEM MEQ OWNS IS `re-assembly`, AT 5.6%.** That is
`prepare()` under the line search, and it is the subject of
`BORDERED-GLOBALISATION-PLAN.md` §12. Everything else MEQ owns and does not
thread is under 3.5%: `other (remainder)` 3.4%, the non-`postProcess` part of
`output` 5.8%, `makeSolver` 1.6%.

**WHAT AN AMDAHL SUM SAYS IS LEFT.** Threading `Reconstruct()` perfectly would
be `1/(1 − 0.185 × (1 − 1/6.5)) = 1.19×`; removing `re-assembly` outright is
`1.06×`. Those are the two largest single levers on this path, and **both are
larger than anything remaining in `THREADING-PLAN.md`**, which is the honest
reason that plan closes rather than continues.

### M-139

**CS-0b: MATCHING THE CONDUCTOR MODEL TAKES DIII-D's BENCHMARK FLOOR FROM
5.785e-03 TO 7.7e-04, AND LANDS IT ON THE QUADRATURE — WHICH IS WHERE
`COIL-SUBTRACTION-PLAN.md` §0a SAID IT WOULD LAND.** M-111 reads MEQ against
freegs4e's **filament** DIII-D at `5.785e-03` and records that it **does not
refine**: 5.785e-03, 5.802e-03, 5.763e-03, 5.803e-03 across a 16× range in dofs.
§0a measured the two conductor models directly and predicted that matching them
would leave `6.580e-04`, the quadrature difference alone. This is that
experiment run end to end.

`race.py --shaped`, `F_diiid_conventional_shaped`, MEQ cold from
`mkcoldguess.py`, measured against freegs4e's own 257² shaped reference:

| rung | dofs | newton | rel `L2` | ... off the coils | `psi_ax` | X-point |
|---|---|---|---|---|---|---|
| k1r0 | 43,632 | **87** | **7.107e-01** | 7.318e-01 | 1.45e-01 | **( 1.0942, 0.1593 )** |
| k1r1 | 175,563 | 14 | 9.059e-04 | 6.199e-04 | 8.23e-05 | ( 1.2000, **+0.9997** ) |
| **k2r0** | 87,264 | 2 | **8.484e-04** | 5.273e-04 | **2.18e-06** | ( 1.2004, −0.9996 ) |
| **k2r1** | 351,126 | 2 | **7.692e-04** | 3.740e-04 | 4.09e-05 | ( 1.2004, −0.9996 ) |
| **k3r0** | 145,440 | 2 | **7.680e-04** | 3.718e-04 | 6.55e-05 | ( 1.2004, −0.9995 ) |
| **k3r1** | 585,210 | 2 | **7.750e-04** | 3.867e-04 | 3.12e-05 | ( 1.2004, −0.9996 ) |
| k2r0a3 | 147,726 | 2 | 7.754e-04 | 3.741e-04 | 4.00e-05 | ( 1.2004, −0.9995 ) |
| k3r0a3 | 172,590 | 2 | 7.830e-04 | 3.909e-04 | 3.01e-05 | ( 1.2004, −0.9996 ) |

**7.5× BETTER GLOBALLY AND 2.0× OFF THE COILS**, against §0a's predicted 9.2×
and 3.4× on the model difference itself. The prediction was of the *difference
between two conductor models*; this is what removing it from a *solve* is worth,
and the two agreeing to within a factor of 1.2 is the part that says §0a
measured the right thing.

**AND IT STILL DOES NOT REFINE, WHICH IS THE RESULT RATHER THAN A
DISAPPOINTMENT.** 8.484e-04 → 7.692e-04 → 7.680e-04 → 7.750e-04 over a **6.7×
range in dofs** and from `k = 2` to `k = 3`: nine per cent, and not monotone.
**The floor moved down by 7.5× and it is still not MEQ's discretisation.**
§0a named what would be left — MEQ integrates its rectangles at 24² tensor
Gauss and `freegs4e.shaped_coil.ShapedCoil` caps at **6 points per triangle**,
so the two evaluate the same model to different precision — and predicted
`6.580e-04`. **Measured, 7.7e-04.** That is agreement to 17% on a number
predicted before the solve existed.

**THE REFERENCE IS NOT THE LIMITER, WHICH HAD TO BE CHECKED BEFORE ANY OF THE
ABOVE MEANT ANYTHING.** The 129² grid is a point-for-point subset of the 257²
one, so the two references difference without interpolating either:
**rel `L2` 9.864e-05**, rel `L∞` 1.567e-04, `psi_axis` 1.181e-06. An order of
magnitude below the 7.7e-04 floor. So the floor is a real disagreement between
the codes and not the reference's own resolution — and a 513² truth, which this
run did not have, would not move it.

**`k = 1` IS WORSE THAN FAILING NOW, AND THAT IS A CHANGE WORTH NAMING.** M-111
records `k = 1` **failing** on all three cases. Here both `k = 1` rungs *return
an answer*: `k1r0` takes **87** Newton steps to a completely different
equilibrium — X-point at `( 1.0942, 0.1593 )` against the reference's active
`( 1.1999, −1.0000 )`, `psi_ax` out by 14.5% — and `k1r1` converges to the
**upper** null, `( 1.2000, +0.9997 )`, which on this nearly up-down symmetric
machine is very nearly degenerate with the active one and still reads
9.059e-04. **A rung that converges to the wrong branch is worse than one that
fails**, and `psi_ax` alone would excuse `k1r1` at 8.23e-05. The X-point column
is what separates them; M-111's own note about `k3r0` on MAST-U is the same
trap.

**EVERY OTHER RUNG TAKES TWO NEWTON STEPS AND FINDS THE SAME NULL** to
`( 1.2004, −0.9996 )` against the reference's `( 1.1999, −1.0000 )`, which is
4e-04 m in `R` and 4e-04 m in `Z`, and `k2r0` reports `psi_ax` to **2.18e-06**.

**THE TIMING COLUMNS OF THIS RUN ARE NOT MEASUREMENTS AND ARE NOT QUOTED.**
A peer's build started partway through the ladder. The accuracy columns are
rates and norms and do not care; the wall clock does, and re-racing DIII-D
shaped against M-111's filament row wants a quiet machine.

### M-140

**MEQ's "EXACT RESTART" WAS ONE BIT OUT ON A QUARTER OF ITS COEFFICIENTS, AND A
TEST NAMED `theMfemFilesRoundTripExactly` PASSED THROUGHOUT.** Found while
asking whether a `GridFunction` could be written to NetCDF for
`COIL-SUBTRACTION-PLAN.md` CS-6, which turns on whether the existing format is
exact.

**16 IS THE NUMBER OF DIGITS A DOUBLE IS ACCURATE TO; 17 IS THE NUMBER NEEDED TO
RECOVER ONE.** `std::numeric_limits<double>::max_digits10` is **17**.
`src/meq/Output.cpp` set `out.precision( 16 )` at both `.gf` sites, under a
comment reasoning that 16 is exact where MFEM's default 8 is not — right about
the 8 and wrong about the 16.

200,000 doubles drawn uniformly from `[ −1, 1 ]`, formatted and parsed back:

| precision | not bit-recovered |
|---|---|
| **16** | **50,204 of 200,000 — 25.1%** |
| 17 | **0** |

**AND IT IS NOT A PROPERTY OF RANDOM DOUBLES ONLY — IT IS MEASURED ON MEQ'S OWN
FIELD.** `theMfemFilesRoundTripExactly` solves a Solov'ev case at `k = 2`,
writes the pair, reads them back with no knowledge of what wrote them, and takes
the worst coefficient difference:

| `writeMfem()` precision | worst coefficient difference |
|---|---|
| 16 | **5.551e-17** — one ulp at that scale |
| **17** | **0.000e+00** |

**THE TEST'S NAME CLAIMED THE PROPERTY AND ITS ASSERTION DID NOT CHECK IT.** It
read `BOOST_TEST( worst < 1.0e-13 )`, so one ulp passed by four orders of
magnitude. A tolerance where the property is an identity — and the name made it
look like coverage, which is worse than no test. It now asserts `worst == 0.0`,
and **that is falsified rather than assumed**: reverting `Output.cpp` to 16 and
rebuilding makes it fail with the 5.551e-17 above, and restoring 17 makes it
pass at `0.000e+00`.

**WHY NOTHING ELSE CAUGHT IT.** `WarmStartConvergence` is the other test of the
restart and asserts on **iteration counts** — a warm start must change the work
and not the answer — so a one-ulp perturbation sits far below the Newton
tolerance and the restart still finishes in one iteration. A legitimate test
that is structurally blind to this: it checks the consequence, not the property.

**AND IT STRENGTHENS CS-6 RATHER THAN ONLY FIXING A BUG.** §9 proposes NetCDF
for the exact field because a `.gf` cannot carry metadata; NetCDF stores
**binary** doubles, so the round-trip question does not arise there at all — and
the file is smaller than the ASCII it replaces. The precision fix makes the
current format correct; the format change makes the question disappear.

### M-141

**THE TWO CONDUCTOR ROUTES AGREE, AND MESHING TO THE COIL IS WORTH 43× AT THE
COARSEST LEVEL.** `COIL-SUBTRACTION-PLAN.md` §0a-pre is a standing requirement —
MEQ must always be able to solve finite-sized coils accurately the old way, and
the split is an option rather than a replacement. This is the measurement that
gives it teeth, and it is only possible because CS-1b put **rectangles** in
`meq::ConductorField`: a rectangle is the conductor that can go either way, so
it is the only one that can be the cross-check. A filament can only be
subtracted.

**THE SUBTRACTED ARM IS THE TRUTH, WHICH IS THE RIGHT WAY ROUND.** With the
coil's own field as the Dirichlet datum and no plasma, the remainder is
identically zero and the total is `psi_c` **exactly**
→ **[§7.4](COIL-SUBTRACTION-PLAN.md)**. The meshed arm solves
`Δ* psi = −mu0 r j_phi` with the same datum and a **top-hat** source, so it
approximates that same field to the mesh's order. The difference between them
therefore *is* the meshed route's discretisation error.

`k = 2`, box `[ 0.6, 1.4 ] × [ −0.4, 0.4 ]`, one rectangle at `( 1.00, 0.00 )`
carrying `1.0e5 A`, compared at quadrature points outside the conductor:

| elements | dofs | **coil edges ON vertices** (±0.10) | rate | **coil edges INSIDE elements** (±0.06) | rate |
|---|---|---|---|---|---|
| 128 | 768 | **1.585e-05** | — | 6.756e-04 | — |
| 512 | 3,072 | **1.678e-06** | 3.24 | 1.276e-03 | **−0.92** |
| 2,048 | 12,288 | **2.612e-07** | 2.68 | 1.917e-04 | 2.73 |

**43× AT THE COARSEST LEVEL, AND THE UNFITTED COLUMN IS NOT EVEN MONOTONE.**
The source is a top hat, so `psi` is not `C²` across the coil's edge; when that
edge cuts element interiors, **which** elements it cuts is not a smooth function
of `h` and the error wanders. That is exactly the behaviour `CLAUDE.md`'s
*Unfitted convergence needs a two-tier rate assertion* records for the extension
path, met here for the same reason — **a fact about cutting a discontinuity
rather than about either conductor route**.

**SO THIS IS A MEASUREMENT OF WHY `tools/mesh/halfdisc.py` MESHES TO THE COILS**,
which that tool does by construction and which `CLAUDE.md` records as the reason
`[[coils]]` emits exact doubles — *"so the mesh aligns to the rectangle
`meq::Coil`'s quadrature uses to the ulp"*. The alignment was known to matter;
this puts 43× on it.

**AND THE RATES SAY THE TWO ROUTES CONVERGE TO ONE FIELD.** 3.24 and 2.68 at
`k = 2`, where `k+1 = 3` is what a smooth solution gives — outside the conductor
`Δ* psi = 0`, so smooth is what it is there. The acceptance asserts only that
the difference **falls**, not that it reaches any particular order, because the
order argument is about the top hat and the requirement is about the two routes
solving the same problem.
