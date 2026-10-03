# Reply: `FirstCrossing()` stops when the midpoint stops moving

Your fix is in, as you wrote it. It is `gf-hdg-subdomains-dev` `238f719f50`,
and in `meq-integration` as `5c2324e921`.

```cpp
const real_t tm = 0.5 * (ta + tb);
if (tm <= ta || tm >= tb) { break; }
```

## What we checked

* **Your diagnosis holds.** The bracket target `tol * length` is absolute,
  and the `1e-13` default is below what single precision can represent.
* **Double is unchanged.** There the bracket can shrink to about
  `2e-16 * length`, so the `1e-13` test always fires first. The full serial
  suite moves by only the new case (578 / 4,312,100). The extension
  miniapp's airfoil, which runs `VertexConePath` through this bisection,
  reproduces every printed digit.
* **The pin runs in double.** Our single-precision unit tests skip, so the
  new case puts a double into your situation instead: a tolerance no double
  bracket can reach (`1e-30`), a cap of 10,000, and a count of level-set
  calls. **It failed with 10,002 calls before the fix and passes with under
  100 after.** A second section checks that the default tolerance still
  stops at the `1e-13` bracket, with the end point within `1e-13`.

## One correction to the scope

The `tol_ = 1e-13, max_iter_ = 100` pair is on three classes, not four:
`LevelSetPath`, `VertexConePath` and `TiltedPath`. `ClosestPointPath` is a
closed form with no bisection. The line you cited as `:476` is
`TiltedPath`'s constructor. All three go through the one `FirstCrossing()`,
so the fix covers every caller.

## On your tolerance sweep

We agree with your conclusion and thank you for measuring it. The tolerance
decides where a path ends, so any value loose enough to save time moves an
answer. Your `DriverAcceptance` 620 against 612 is a clean demonstration of
that. Staying on the default is right, and with this fix the default is now
also cheap in single.

## What `5c2324e921` also carries, if you have not taken it yet

* `HDGExtensionIntegrator::SetProjectedLift()` (opt-in): lifts
  `Π_K(C u_h)` instead of integrating `C` along the path. It is the fix
  where `C` diverges like 1/x towards a Dirichlet boundary.
* `CheckPathAngles()`, and the opt-in `TiltedPath` and path-checking
  `MarkLevelSetSubdomain()`, for paths that run along their own face.
* The datum's two halves now share one face rule by default
  (`ExtensionDatumLFIntegrator` for new code).
* `VertexConePath`'s brief now states that the cone is off unless
  `use_cone_` is set.
