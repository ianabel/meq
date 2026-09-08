# MEQ's published tables

**Every number MEQ claims, in one place.** These are the measurement tables that
used to sit inline in `CLAUDE.md`. They were moved out because a maintainer does
not need them on every reading and does need to know where they are — the
argument each one supports is still in `CLAUDE.md`, beside a pointer to the row.

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

