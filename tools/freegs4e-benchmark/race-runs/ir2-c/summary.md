## machine-c-mast-shaped

### cold start

**freegs4e**

| rung | size | median wall s | repeats | accepted | dI vs own finest | psi_ax vs own finest |
|---|---|---|---|---|---|---|
| 65^2 | 4225 | 7.35 | 3 | 3/3 | 1.14e-05 | 3.69e-05 |
| 129^2 | 16641 | 11.40 | 3 | 3/3 | 4.40e-06 | 8.05e-07 |
| 257^2 | 66049 | 33.33 | 3 | 3/3 | 1.49e-07 | 2.10e-07 |
| 513^2 (datum) | 263169 | 174.98 | 1 | 1/1 | -- | -- |

**freegsnke**

| rung | size | median wall s | repeats | accepted | dI vs own finest | psi_ax vs own finest |
|---|---|---|---|---|---|---|
| nx65 | 4225 | 11.79 | 3 | 3/3 | 6.73e-06 | 3.74e-05 |
| nx129 | 16641 | 58.93 | 3 | 3/3 | 4.52e-06 | 5.39e-07 |
| nx257 (datum) | 66049 | 237.97 | 1 | 1/1 | -- | -- |

**meq** -- the finest rung has no accepted run: no self-error can be reported

| rung | size | median wall s | repeats | accepted | dI vs own finest | psi_ax vs own finest |
|---|---|---|---|---|---|---|
| k2 h1.00 | None | -- | 3 | 0/3 | -- | -- |
| k2 h0.50 | None | -- | 3 | 0/3 | -- | -- |
| k3 h1.00 | None | -- | 3 | 0/3 | -- | -- |
| k3 h0.50 | None | -- | 3 | 0/3 | -- | -- |
| k4 h0.50 | None | -- | 3 | 0/3 | -- | -- |
| k4 h0.35 (datum) | None | -- | 1 | 0/1 | -- | -- |

**nice** -- the finest rung has no accepted run: no self-error can be reported

| rung | size | median wall s | repeats | accepted | dI vs own finest | psi_ax vs own finest |
|---|---|---|---|---|---|---|
| area x1 | 3271 | -- | 3 | 0/3 | -- | -- |
| area x1/4 | 12432 | -- | 3 | 0/3 | -- | -- |
| area x1/16 | 48651 | -- | 3 | 0/3 | -- | -- |
| area x1/64 (datum) | 192694 | -- | 1 | 0/1 | -- | -- |

**Cheapest run reaching each self-consistent accuracy (currents), median wall s:**

| target | freegs4e | freegsnke | meq | nice |
|---|---|---|---|---|
| 1e-02 | 7.35 (65^2) | 11.79 (nx65) | -- | -- |
| 1e-03 | 7.35 (65^2) | 11.79 (nx65) | -- | -- |
| 1e-04 | 7.35 (65^2) | 11.79 (nx65) | -- | -- |
| 1e-05 | 11.40 (129^2) | 11.79 (nx65) | -- | -- |

**Cross-code floor** (finest accepted runs, dI over the free supplies):

* freegs4e - freegsnke: 2.02e-07

No accuracy target below 2.0e-07 distinguishes the codes' ANSWERS; the table above compares each code with itself.

### warm start

**freegs4e**

| rung | size | median wall s | repeats | accepted | dI vs own finest | psi_ax vs own finest |
|---|---|---|---|---|---|---|
| 65^2 | 4225 | 6.29 | 3 | 3/3 | 5.45e-05 | 4.79e-05 |
| 129^2 | 16641 | 10.32 | 3 | 3/3 | 7.58e-07 | 4.86e-06 |
| 257^2 | 66049 | 30.65 | 3 | 3/3 | 1.85e-07 | 6.22e-07 |
| 513^2 (datum) | 263169 | 201.61 | 1 | 1/1 | -- | -- |

**freegsnke**

| rung | size | median wall s | repeats | accepted | dI vs own finest | psi_ax vs own finest |
|---|---|---|---|---|---|---|
| nx65 | 4225 | 7.59 | 3 | 3/3 | 5.48e-05 | 4.73e-05 |
| nx129 | 16641 | 21.60 | 3 | 3/3 | 1.16e-06 | 5.59e-06 |
| nx257 (datum) | 66049 | 141.61 | 1 | 1/1 | -- | -- |

**Cheapest run reaching each self-consistent accuracy (currents), median wall s:**

| target | freegs4e | freegsnke |
|---|---|---|
| 1e-02 | 6.29 (65^2) | 7.59 (nx65) |
| 1e-03 | 6.29 (65^2) | 7.59 (nx65) |
| 1e-04 | 6.29 (65^2) | 7.59 (nx65) |
| 1e-05 | 10.32 (129^2) | 21.60 (nx129) |

**Cross-code floor** (finest accepted runs, dI over the free supplies):

* freegs4e - freegsnke: 1.41e-07

No accuracy target below 1.4e-07 distinguishes the codes' ANSWERS; the table above compares each code with itself.
