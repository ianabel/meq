# The q-driven Jacobian by shape derivative

A design, not yet built. It replaces the differenced extraction in
`meq::SafetyFactorEquilibrium::differentiate()` with the derivative of each
surface's `q/g` along each knot's tangent, computed once on the base family.
Two problems motivate it:

* **The floor.** Differencing a single-precision extraction leaves the
  Jacobian's rows noise-dominated. No step does better than about 4e-5 of
  `max |J|` on the fixed box, or 1e-3 on the limited and diverted boxes
  (MEASUREMENTS.md M-QS-1).
* **The cost.** The `2m` extractions a step are most of the loop's time
  (TODO, *The q-driven Jacobian's cost is the extraction*).

## What is differentiated

The extraction delivers, per surface `i` at label `Psi_N,i`,

    ( q/g )_i = V'_i < R^-2 >_i / 4 pi^2 = ( 1/4 pi^2 ) sum_k w_k / R_k^2,
    w_k = 2 pi R_k | dx/dtheta |_k dtheta / | grad psi |_k

on the contour `psi~ = c_i`, with `c_i = psi_b + Psi_N,i ( psi_ax - psi_b )`.
`psi~` is the kept continuous interpolant. `| grad psi |` is read from the
flux `R | q~ |`, as the tracer's weights read it, and it is not
`| grad psi~ |`. The derivative has to be the derivative of THIS map: the
level read on `psi~`, the metric read on `q~`. A clean shape derivative of
`int dl/( R | grad psi~ | )` would disagree with re-solved differences at the
discretisation's `O( h^( k+1 ) )`.

Along a tangent `t` -- the fields `psi~ + t dpsi~` and `q~ + t dq~`, the
labels moving by `dc_i = dpsi_b + Psi_N,i ( dpsi_ax - dpsi_b )` -- write the
surface integral as `I( c ) = closed-int_{psi~ = c} f dl`, with
`f = 1/( R^2 | q~ | )` (so that `( q/g ) = I/2 pi`). Then

    dI/dt = closed-int [ df + ( d_n f + kappa f ) V_n ] dl,

where:

* `df = -f ( q~ . dq~ )/| q~ |^2` is the integrand's own change;
* `V_n = ( dc - dpsi~ )/| grad psi~ |` is the surface's normal velocity,
  with `n = grad psi~/| grad psi~ |`;
* `kappa = div n` is its curvature, from the Hessian of `psi~`;
* `d_n f` is `f`'s normal derivative, from the gradient of `q~`, plus the
  explicit `1/R^2`.

Everything is an element polynomial evaluated at the base family's nodes:

| quantity | source |
|---|---|
| `psi~` value, gradient and Hessian | the continuous interpolant |
| `q~` value and gradient | the continuous interpolant |
| `dpsi~` value and `dq~` value | the tangent's interpolant |

Conductors, where the solve subtracts them, add `psi_c`'s gradient and
Hessian to `psi~`'s, and `q_c` and its gradient to `q~`. They are fixed along
every tangent, so they add nothing to `dpsi~` or `dq~`. `ConductorField` gives
`psi_c`, its gradient and `q_c` analytically; the Hessian and the gradient of
`q_c` are differenced in double, as the X-point border's corner already is.

## What the quadrature costs it

The family's quadrature is the trapezoid rule in the ray angle about the
located axis, and its nodes move with the axis. The continuous shape
derivative drops two things:

* the quadrature's own dependence on the origin, which M-391's ray-origin
  term measured at 1e-5 to 4e-5 of the Jacobian;
* the trapezoid error's derivative.

Both are `O( h^( k+1 ) )`-small differences between the discrete map and its
continuous idealisation. Re-solved differences in double resolve them
(`SafetyFactorSolver`'s 1e-6 bound held to 1.2e-7 with the origin moved). So
the acceptance is not "equal to the differenced Jacobian". It is:

1. In double, the shape-derivative Jacobian is re-solved differences to
   within the same 1e-6 of its size where the reference is differentiable
   (`QDrivenFreeEdge`'s precondition). If the dropped terms exceed that,
   carry the origin term explicitly: it is a sum over nodes of
   `d( node )/d( origin ) . grad( integrand )`, computable from the same
   samples.
2. In single, the rows' noise falls from 4e-5 to the evaluation's own
   (1.5e-6 on the fixed box, 1e-5 on the free-edge boxes, M-QS-1).
3. The loops of `QDrivenFreeEdge`, the limited-box diagonal and DeviceParity's
   q-driven card arm close in single.

## Where the band breaks it

On a curved boundary the outer surfaces are traced partly in the band, where
the fields are `BandExtension::TransferLift` continuations (`ContourTracer`)
and have no element polynomial to take a Hessian of. For band nodes the
derivative falls back to central differences of the tracer's own
`sampleField()` gradient over a fixed physical step. That carries the
single-precision noise on those nodes only. The share of band nodes is
printed per family.

## Cost

Per Jacobian:

* **One base family**, already traced for the residual.
* **`m` interpolant applies**, for the tangents (`continuousInterpolant()`).
* **The node evaluations**:
  * **The base family's samples:** `m surfaces x A angles` (`A` = 128 by
    default), each a value, a gradient and a Hessian of `psi~` and of each
    component of `q~`, taken once.
  * **Each direction's samples:** `dpsi~` and `dq~` at the same nodes,
    `m x m x A` evaluations of one element polynomial each.
  * **The work:** at `m = 19`, `A = 128` that is 46,000 element evaluations a
    Jacobian, against today's `2m = 38` full extractions of `m` traced
    surfaces each (each trace a corrector walk plus `A` ray solves).
  * **The expected gain:** the extraction's 39 s of the example's 60 s (M-399)
    should fall to under a second.
* **Device route:** a kernel over ( node, direction ) pairs gathering each
  node's element dofs, the same shape as `meq::ContinuousInterpolant::apply()`'s
  gather. The node set and its element and reference coordinates are fixed
  per Jacobian and uploaded once. Parity arms on `debug` and `cuda` against the
  host loop, bit-equal on `debug`.

## Stages

1. **Host, interior nodes only.** Held to `SafetyFactorSolver`'s re-solved
   differences in double on the box. The ray-origin term is added if needed.
2. **Free edge.** `dc` carries `dpsi_b`, and conductors enter. Held on
   `QDrivenFreeEdge` (judged entries) and on the limited machine.
3. **The band fallback.** `QDrivenCurvedBoundary`.
4. **Single precision.** The `QDrivenFreeEdge` loops, the limited-box
   diagonal, the card arm's repeatability over three runs.
5. **The device kernel**, with its parity arms.

The differenced Jacobian stays as `Options::differenced` and as the test that
checks this one.
