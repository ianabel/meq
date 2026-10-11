## machine-f-diiid-shaped

### cold start

**freegs4e**

| rung | size | median wall s | repeats | accepted | dI vs own finest | psi_ax vs own finest |
|---|---|---|---|---|---|---|
| 65^2 | 4225 | 6.55 | 3 | 3/3 | 2.88e-05 | 7.75e-06 |
| 129^2 | 16641 | 9.19 | 3 | 3/3 | 2.00e-06 | 1.53e-06 |
| 257^2 | 66049 | 22.69 | 3 | 3/3 | 3.57e-07 | 9.51e-09 |
| 513^2 (datum) | 263169 | 132.86 | 1 | 1/1 | -- | -- |

**freegsnke**

| rung | size | median wall s | repeats | accepted | dI vs own finest | psi_ax vs own finest |
|---|---|---|---|---|---|---|
| nx65 | 4225 | 13.30 | 3 | 3/3 | 1.39e-05 | 6.62e-06 |
| nx129 | 16641 | 66.44 | 3 | 3/3 | 3.54e-06 | 1.16e-06 |
| nx257 (datum) | 66049 | 24052.66 | 1 | 1/1 | -- | -- |

**meq**

| rung | size | median wall s | repeats | accepted | dI vs own finest | psi_ax vs own finest |
|---|---|---|---|---|---|---|
| k2 h1.00 | 29088 | 80.26 | 3 | 3/3 | 3.07e-04 | 6.63e-05 |
| k2 h0.50 | 67398 | 69.19 | 3 | 3/3 | 2.95e-04 | 1.51e-05 |
| k3 h1.00 | 48480 | 128.33 | 3 | 3/3 | 1.33e-04 | 4.33e-06 |
| k3 h0.50 | 112330 | 306.93 | 3 | 3/3 | 4.76e-05 | 2.72e-06 |
| k4 h0.50 | 168495 | 440.60 | 3 | 3/3 | 2.11e-05 | 2.06e-06 |
| k4 h0.35 (datum) | 289050 | 795.41 | 1 | 1/1 | -- | -- |

**nice**

| rung | size | median wall s | repeats | accepted | dI vs own finest | psi_ax vs own finest |
|---|---|---|---|---|---|---|
| area x1 | 4279 | 1.56 | 3 | 3/3 | 1.44e-02 | 5.74e-04 |
| area x1/4 | 16663 | 5.51 | 3 | 3/3 | 2.57e-03 | 1.44e-04 |
| area x1/16 | 65533 | 25.81 | 3 | 3/3 | 6.42e-04 | 1.60e-05 |
| area x1/64 (datum) | 259737 | 102.05 | 1 | 1/1 | -- | -- |

**Cheapest run reaching each self-consistent accuracy (currents), median wall s:**

| target | freegs4e | freegsnke | meq | nice |
|---|---|---|---|---|
| 1e-02 | 6.55 (65^2) | 13.30 (nx65) | 69.19 (k2 h0.50) | 5.51 (area x1/4) |
| 1e-03 | 6.55 (65^2) | 13.30 (nx65) | 69.19 (k2 h0.50) | 25.81 (area x1/16) |
| 1e-04 | 6.55 (65^2) | 13.30 (nx65) | 306.93 (k3 h0.50) | -- |
| 1e-05 | 9.19 (129^2) | 66.44 (nx129) | -- | -- |

**Cross-code floor** (finest accepted runs, dI over the free supplies):

* freegs4e - freegsnke: 2.66e-07
* freegs4e - meq: 1.46e-05
* freegs4e - nice: 4.26e-04
* freegsnke - meq: 1.43e-05
* freegsnke - nice: 4.26e-04
* meq - nice: 4.19e-04

No accuracy target below 4.3e-04 distinguishes the codes' ANSWERS; the table above compares each code with itself.

### warm start

**freegs4e**

| rung | size | median wall s | repeats | accepted | dI vs own finest | psi_ax vs own finest |
|---|---|---|---|---|---|---|
| 65^2 | 4225 | 4.51 | 3 | 3/3 | 2.53e-05 | 1.46e-05 |
| 129^2 | 16641 | 6.63 | 3 | 3/3 | 1.31e-06 | 2.79e-07 |
| 257^2 | 66049 | 18.14 | 3 | 3/3 | 2.64e-07 | 1.31e-07 |
| 513^2 (datum) | 263169 | 104.72 | 1 | 1/1 | -- | -- |

**freegsnke**

| rung | size | median wall s | repeats | accepted | dI vs own finest | psi_ax vs own finest |
|---|---|---|---|---|---|---|
| nx65 | 4225 | 5.66 | 3 | 3/3 | 2.01e-05 | 9.57e-06 |
| nx129 | 16641 | 24.76 | 3 | 3/3 | 1.31e-06 | 3.28e-07 |
| nx257 (datum) | 66049 | 17212.78 | 1 | 1/1 | -- | -- |

**meq**

| rung | size | median wall s | repeats | accepted | dI vs own finest | psi_ax vs own finest |
|---|---|---|---|---|---|---|
| k2 h1.00 | 29088 | 3.06 | 3 | 3/3 | 6.21e-04 | 6.12e-05 |
| k2 h0.50 | 67398 | 6.79 | 3 | 3/3 | 1.33e-04 | 4.82e-06 |
| k3 h1.00 | 48480 | 5.08 | 3 | 3/3 | 1.45e-04 | 4.80e-06 |
| k3 h0.50 | 112330 | 11.59 | 3 | 3/3 | 5.11e-05 | 3.56e-06 |
| k4 h0.50 | 168495 | 16.85 | 3 | 3/3 | 2.19e-05 | 2.26e-06 |
| k4 h0.35 (datum) | 289050 | 30.95 | 1 | 1/1 | -- | -- |

**nice**

| rung | size | median wall s | repeats | accepted | dI vs own finest | psi_ax vs own finest |
|---|---|---|---|---|---|---|
| area x1 | 4279 | 1.34 | 3 | 3/3 | 1.59e-02 | 7.16e-04 |
| area x1/4 | 16663 | 4.13 | 3 | 3/3 | 2.28e-03 | 1.42e-04 |
| area x1/16 | 65533 | 18.52 | 3 | 3/3 | 6.47e-04 | 1.84e-05 |
| area x1/64 (datum) | 259737 | 71.13 | 1 | 1/1 | -- | -- |

**Cheapest run reaching each self-consistent accuracy (currents), median wall s:**

| target | freegs4e | freegsnke | meq | nice |
|---|---|---|---|---|
| 1e-02 | 4.51 (65^2) | 5.66 (nx65) | 3.06 (k2 h1.00) | 4.13 (area x1/4) |
| 1e-03 | 4.51 (65^2) | 5.66 (nx65) | 3.06 (k2 h1.00) | 18.52 (area x1/16) |
| 1e-04 | 4.51 (65^2) | 5.66 (nx65) | 11.59 (k3 h0.50) | -- |
| 1e-05 | 6.63 (129^2) | 24.76 (nx129) | -- | -- |

**Cross-code floor** (finest accepted runs, dI over the free supplies):

* freegs4e - freegsnke: 1.98e-07
* freegs4e - meq: 1.44e-05
* freegs4e - nice: 4.45e-04
* freegsnke - meq: 1.44e-05
* freegsnke - nice: 4.45e-04
* meq - nice: 4.38e-04

No accuracy target below 4.5e-04 distinguishes the codes' ANSWERS; the table above compares each code with itself.
