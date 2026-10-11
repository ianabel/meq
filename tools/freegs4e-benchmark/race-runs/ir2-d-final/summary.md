## machine-d-tcv-shaped

### cold start

**freegs4e**

| rung | size | median wall s | repeats | accepted | dI vs own finest | psi_ax vs own finest |
|---|---|---|---|---|---|---|
| 65^2 | 4225 | 20.48 | 3 | 3/3 | 4.57e-05 | 6.54e-06 |
| 129^2 | 16641 | 21.70 | 3 | 3/3 | 3.97e-06 | 6.79e-07 |
| 257^2 | 66049 | 43.47 | 3 | 3/3 | 3.91e-07 | 1.68e-07 |
| 513^2 (datum) | 263169 | 438.35 | 1 | 1/1 | -- | -- |

**freegsnke**

| rung | size | median wall s | repeats | accepted | dI vs own finest | psi_ax vs own finest |
|---|---|---|---|---|---|---|
| nx65 | 4225 | 13.16 | 3 | 3/3 | 6.06e-05 | 1.01e-05 |
| nx129 | 16641 | 243.97 | 3 | 3/3 | 3.72e-06 | 3.77e-07 |
| nx257 (datum) | 66049 | 70459.36 | 1 | 1/1 | -- | -- |

**meq**

| rung | size | median wall s | repeats | accepted | dI vs own finest | psi_ax vs own finest |
|---|---|---|---|---|---|---|
| k2 h1.00 | 30558 | 6.78 | 3 | 3/3 | 6.57e-04 | 1.15e-04 |
| k2 h0.50 | 83598 | 16.52 | 3 | 3/3 | 1.95e-04 | 3.75e-05 |
| k3 h1.00 | 50930 | 10.88 | 3 | 3/3 | 6.12e-05 | 2.59e-06 |
| k3 h0.50 | 139330 | 29.31 | 3 | 3/3 | 2.93e-06 | 4.66e-07 |
| k4 h0.50 | 208995 | 51.97 | 3 | 3/3 | 4.19e-07 | 2.18e-07 |
| k4 h0.35 (datum) | 391365 | 74.10 | 1 | 1/1 | -- | -- |

**nice**

| rung | size | median wall s | repeats | accepted | dI vs own finest | psi_ax vs own finest |
|---|---|---|---|---|---|---|
| area x1 | 6825 | 1.93 | 3 | 3/3 | 5.96e-03 | 5.97e-04 |
| area x1/4 | 25473 | 6.90 | 3 | 3/3 | 1.65e-03 | 3.04e-05 |
| area x1/16 (datum) | 99994 | 37.51 | 1 | 1/1 | -- | -- |

**Cheapest run reaching each self-consistent accuracy (currents), median wall s:**

| target | freegs4e | freegsnke | meq | nice |
|---|---|---|---|---|
| 1e-02 | 20.48 (65^2) | 13.16 (nx65) | 6.78 (k2 h1.00) | 1.93 (area x1) |
| 1e-03 | 20.48 (65^2) | 13.16 (nx65) | 6.78 (k2 h1.00) | -- |
| 1e-04 | 20.48 (65^2) | 13.16 (nx65) | 10.88 (k3 h1.00) | -- |
| 1e-05 | 21.70 (129^2) | 243.97 (nx129) | 29.31 (k3 h0.50) | -- |

**Cross-code floor** (finest accepted runs, dI over the free supplies):

* freegs4e - freegsnke: 2.58e-07
* freegs4e - meq: 3.53e-06
* freegs4e - nice: 6.73e-04
* freegsnke - meq: 3.32e-06
* freegsnke - nice: 6.73e-04
* meq - nice: 6.74e-04

No accuracy target below 6.7e-04 distinguishes the codes' ANSWERS; the table above compares each code with itself.

### warm start

**freegs4e**

| rung | size | median wall s | repeats | accepted | dI vs own finest | psi_ax vs own finest |
|---|---|---|---|---|---|---|
| 65^2 | 4225 | 14.34 | 3 | 3/3 | 4.72e-05 | 2.08e-05 |
| 129^2 | 16641 | 18.22 | 3 | 3/3 | 3.12e-06 | 9.24e-07 |
| 257^2 | 66049 | 39.85 | 3 | 3/3 | 1.43e-06 | 1.12e-07 |
| 513^2 (datum) | 263169 | 252.44 | 1 | 1/1 | -- | -- |

**freegsnke**

| rung | size | median wall s | repeats | accepted | dI vs own finest | psi_ax vs own finest |
|---|---|---|---|---|---|---|
| nx65 | 4225 | 7.06 | 3 | 3/3 | 5.90e-05 | 1.54e-05 |
| nx129 | 16641 | 76.34 | 3 | 3/3 | 1.53e-06 | 9.91e-08 |
| nx257 (datum) | 66049 | 127.94 | 1 | 1/1 | -- | -- |

**meq**

| rung | size | median wall s | repeats | accepted | dI vs own finest | psi_ax vs own finest |
|---|---|---|---|---|---|---|
| k2 h1.00 | 30558 | 3.46 | 3 | 3/3 | 2.28e-03 | 3.36e-04 |
| k2 h0.50 | 83598 | 8.85 | 3 | 3/3 | 4.19e-05 | 9.28e-06 |
| k3 h1.00 | 50930 | 6.18 | 3 | 3/3 | 4.25e-05 | 3.58e-06 |
| k3 h0.50 | 139330 | 16.67 | 3 | 3/3 | 2.58e-06 | 2.05e-06 |
| k4 h0.50 | 208995 | 23.76 | 3 | 3/3 | 3.20e-07 | 2.42e-08 |
| k4 h0.35 (datum) | 391365 | 48.72 | 1 | 1/1 | -- | -- |

**nice**

| rung | size | median wall s | repeats | accepted | dI vs own finest | psi_ax vs own finest |
|---|---|---|---|---|---|---|
| area x1 | 6825 | 1.82 | 3 | 3/3 | 5.35e-03 | 4.80e-04 |
| area x1/4 | 25473 | 6.35 | 3 | 3/3 | 2.07e-03 | 1.40e-04 |
| area x1/16 (datum) | 99994 | 33.53 | 1 | 1/1 | -- | -- |

**Cheapest run reaching each self-consistent accuracy (currents), median wall s:**

| target | freegs4e | freegsnke | meq | nice |
|---|---|---|---|---|
| 1e-02 | 14.34 (65^2) | 7.06 (nx65) | 3.46 (k2 h1.00) | 1.82 (area x1) |
| 1e-03 | 14.34 (65^2) | 7.06 (nx65) | 6.18 (k3 h1.00) | -- |
| 1e-04 | 14.34 (65^2) | 7.06 (nx65) | 6.18 (k3 h1.00) | -- |
| 1e-05 | 18.22 (129^2) | 76.34 (nx129) | 16.67 (k3 h0.50) | -- |

**Cross-code floor** (finest accepted runs, dI over the free supplies):

* freegs4e - freegsnke: 1.36e-06
* freegs4e - meq: 3.47e-06
* freegs4e - nice: 7.52e-04
* freegsnke - meq: 2.11e-06
* freegsnke - nice: 7.53e-04
* meq - nice: 7.53e-04

No accuracy target below 7.5e-04 distinguishes the codes' ANSWERS; the table above compares each code with itself.
