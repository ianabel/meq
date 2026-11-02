Conductor models
================

A machine's poloidal field coils are prescribed currents, and MEQ
can carry them three ways. ``[conductors] Model`` chooses between them — see
:ref:`config-conductors` for the keys — and this
chapter is what each one does to the equation, what it is exact for, and what it
costs.

.. list-table::
   :header-rows: 1
   :widths: 18 42 40

   * - ``Model``
     - The conductor is
     - and the mesh
   * - ``"meshed"`` (the default)
     - a rectangle of uniform current density, a **domain source** in the
       interior equation
     - must resolve each rectangle, and is graded around it
   * - ``"subtracted"``
     - the same rectangle, its flux :math:`\psi_c` **integrated exactly** and
       taken out of the equation
     - need not contain it at all
   * - ``"filament"``
     - a point filament at each rectangle's centre — or a stack of them,
       :ref:`below <conductors-filament-stack>` — carrying the same total
       current, subtracted the same way
     - need not contain it at all

The first two describe the **same** machine and must agree to the
discretisation. The third describes a **different** one: a filament set and a
rectangle set carrying the same currents are two machines, and the ``.nc`` file
records which produced it in its ``conductor_model`` attribute so that two
outputs can be told apart before they are differenced.

.. _conductors-meshed:

The meshed conductor is a first-class route, and stays one
----------------------------------------------------------

**MEQ can always solve finite-sized coils accurately with the coils in the
mesh, and the split is an option beside that rather than a replacement for
it.** ``"meshed"`` is the default for a rectangle, it is exercised and measured
like everything else, and nothing about the split is allowed to make it pay:
every place the split adds :math:`\psi_c` does so through a pointer that is null
unless a subtracting model was asked for, so a run that does not ask is
bit-for-bit the run it would have been without the split existing — and the
test suite asserts that rather than assuming it.

Three reasons make that a requirement and not a courtesy:

* **It is the subtraction's instrument.** The same machine solved with a
  rectangle meshed and with it subtracted must reach the same equilibrium to
  the discretisation. That comparison is the only check that the split changed
  the *representation* and not the *answer*, and it exists only while both
  routes work. A rectangle is the one conductor that can go either way, which
  is what makes the comparison possible; a filament can only be subtracted, so
  it can never be the cross-check.
* **It needs nothing evaluated per quadrature point.** The meshed source is a
  constant current density on the coil's elements; the split evaluates a
  Green's-function integral at every point a consumer asks about.
* **It is the only route that puts the conductor's own current into the
  interior equation**, which is where a force on the conductor would be read.

.. _conductors-split:

The split
---------

Write the poloidal flux as the conductors' field plus a remainder,

.. math::

   \psi = \psi_c + \psi_p ,

where :math:`\psi_c` is the field of the conductors **alone**, in free space,
and :math:`\psi_p` is what MEQ solves for. The operator is linear and a forward
solve's conductor currents are prescribed, so with the coil's source written as
in :doc:`configuration` (``[[coils]]``),

.. math::

   \dstar\psi_p = \dstar\psi - \dstar\psi_c
       = -\mu_0 R \left( J_{\mathrm{plasma}} + J_{\mathrm{coil}} \right)
         - \left( -\mu_0 R\, J_{\mathrm{coil}} \right)
       = -\mu_0 R\, J_{\mathrm{plasma}}(\psi_c + \psi_p) .

**The conductor's source cancels exactly.** It is in neither the mesh nor the
right-hand side, and the remainder sees a bounded source supported on the
plasma alone. For a filament that is the whole trick: :math:`\psi_c` carries the
logarithm a line current puts into :math:`\psi`, so :math:`\psi_p` is smooth
where :math:`\psi` is not.

Two things in that last line are easy to get wrong and MEQ does not:

* **The plasma current is a function of the physical flux.** The source is
  evaluated at :math:`\psi_c + \psi_p`, never at the solved :math:`\psi_p`. A
  source evaluated at the remainder converges at the full rate to a different
  equilibrium.
* **So is its Jacobian.** :math:`\partial F/\partial\psi_p` at :math:`\psi_p`
  is :math:`\partial F/\partial\psi` *evaluated at* :math:`\psi_c + \psi_p`. The
  factor :math:`\mathrm{d}(\psi_c + \psi_p)/\mathrm{d}\psi_p` is exactly one,
  and that is true and irrelevant: the argument still moves. A Jacobian taken at
  the wrong argument converges either way, to a different problem or merely
  more slowly, which is why neither shows in a convergence table.

:math:`\psi_c` needs no closed form. For a filament it is the textbook
current-loop flux, which MEQ evaluates through Carlson's symmetric elliptic
integrals because they stay accurate arbitrarily close to the loop, where the
complete-integral form loses every digit; for a rectangle it is that kernel
integrated over the cross-section by a tensor Gauss rule, finite everywhere
including inside the coil. ``[conductors] QuadratureOrder`` sets the rule.
Measured clear of the conductors, the shipped order agrees with a far higher
one to round-off, and raising it moves a solved split by orders of magnitude
less than the discretisation error — so the order is not a knob that changes
the answer, which you can check on your own machine the same way.

**And it is computed once per mesh, not per Newton step.** The conductor
currents do not move during a forward solve; only the profile scale does. So
:math:`\psi_c` at every potential degree of freedom and every source
quadrature point is a precompute, amortised over every iteration — and,
through ``CacheFile``, over every *run* on the same machine
(:ref:`below <conductors-cache>`).

.. _conductors-exact:

What the split is exact for, and what it cannot be
--------------------------------------------------

The split resolves the conductors to machine precision rather than to the
mesh's order, at no mesh cost — but only at a point where :math:`\psi_c` is
**evaluated**. A finite-element function cannot hold it: projecting
:math:`\psi_c` onto the discrete space is :math:`\psi_c` to the mesh's order,
which is exactly the error the split exists to remove, and for a filament no
polynomial space contains a logarithm at all. So MEQ never forms
:math:`\psi_c + \psi_p` as one grid function, and the outputs divide along that
line:

.. list-table::
   :header-rows: 1
   :widths: 40 16 44

   * - Where
     - Exact?
     - Why
   * - the source and its Jacobian, the critical points, the limiter and
       X-point borders, the plasma-current integral, the flux-surface tracer
       and the averages over its surfaces, the gridded ``.nc``
     - **yes**
     - each evaluates :math:`\psi_c` in closed form at the point it is reading
   * - ``<stem>_psi.gf``, the VTK, the ``.nc`` coefficient variables
     - **no, by design**
     - each is the field's own coefficients, so it holds :math:`\psi_p`
   * - ``<stem>_psi_total.gf``
     - at the nodes
     - :math:`\psi_h + I_h(\psi_c)`, interpolated between them, for looking at

The second row is why the ``.nc`` carries its conductor table rather than a
flag: a remainder is meaningful only together with the conductors it is a
remainder *from*, and :math:`\psi_c` is recoverable from neither the mesh, nor
the space, nor the coefficients. :doc:`output` has the variables.

The cleanest statement of the exactness is a test the suite runs. Put a
filament inside the domain, give the solve no plasma, and hand it the
filament's own field as the boundary datum. Then
:math:`\dstar\psi_p = 0` with :math:`\psi_p = 0` on :math:`\Gamma`, so the
remainder is identically zero and the reported total is the filament's field
**exactly** — at every point, including arbitrarily close to the conductor,
which no meshed solve of any degree can do. The suite asserts that at zero
tolerance, beside a control that skips the subtraction and gets a smooth,
convergent and entirely plausible wrong answer: the harmonic extension of the
datum.

.. _conductors-total:

Every reader of the flux reads the total
----------------------------------------

Under the split the solved field is no longer the physical field, and every
consumer that forgets is a silent wrong answer rather than a failure — a level
set of the remainder is a smooth closed curve that is not a flux surface, and a
null of the remainder's gradient need not lie anywhere near the physical
X-point. So the library makes the choice explicit at every read. A field view
is constructed by name — ``meq::PotentialView::physical( solved, conductors )``
for :math:`\psi_c + \psi_p`, or ``meq::PotentialView::remainder( solved )`` for
:math:`\psi_p` alone, with ``meq::FluxVectorView`` the same for :math:`q` — and
there is no default and no implicit conversion from a grid function. Without
conductors the two are bit-identical to reading the grid function directly.

**The physical flux** is what the source and its Jacobian, the magnetic-axis and
X-point search, the support fill, the limiter contact and its border, the
X-point border, the plasma-current integral and its derivatives, the
flux-surface tracer and the gridded output read.

**The remainder** is read deliberately in three places, each for a reason:

* **The residual error estimator** measures the discretisation of the problem
  actually solved, and that problem is the remainder's. Its one exception is the
  exterior transmission indicator, which compares a normal flux across
  :math:`\Gamma` and so adds :math:`q_c` back — without it the indicator reports
  a perfectly matched boundary as mismatched by the conductors' own flux, and
  the adaptive loop refines :math:`\Gamma` chasing an error that is not there.
* **The jump of** :math:`\psi_h` **across a face** — :math:`\psi_c` is the same
  on both sides and cancels.
* **A border row's derivative.** The limiter contact enters Newton twice: as a
  value in the residual, which must be the physical :math:`\psi_c + \psi_p`,
  and as a functional applied to search directions, which must not carry the
  constant :math:`\psi_c` — a directional derivative of a constant is zero. MEQ
  keeps the two as separate functions, since a residual and a Jacobian that
  describe different functions make Newton drive every other border while that
  one stands still.

**On a half-disc that reaches the axis**, :math:`q_c = (1/R)\,\gradbar\psi_c`
is :math:`0/0` on :math:`R = 0`, so the gridded output and the tracer take
:math:`B_c` through the closed-form axis limit. The critical-point search
evaluates an element's polynomial slightly outside the element on purpose, and
an iterate that leaves the half-plane :math:`R \ge 0` is abandoned rather than
continued on the remainder alone: a point off the half-plane is not a point of
the machine. :math:`\psi_c` itself is even in :math:`R`, so where a boundary
datum lands a hair across the axis it is evaluated by reflection, exactly; its
gradient is not even, and is refused there instead.

**Where a candidate is screened for being inside a conductor** — the axis
search and the support fill's seed both refuse a point inside copper, because a
rectangle's field has its own O-point in its middle — the screen reads the
machine's conductors, not the source's, since a subtracting model removes them
from the source. Filaments are not screened: a filament has no interior, and
excluding its neighbourhood would need a clearance rule MEQ does not make.

.. _conductors-exterior:

The boundary datum and the exterior coupling
--------------------------------------------

On a fixed boundary the datum is for the **physical** flux, as
``[boundary]`` says, and what MEQ imposes on :math:`\psi_p` is
:math:`\psi|_\Gamma - \psi_c|_\Gamma`. A subtracted conductor may then sit
anywhere.

Under ``[boundary.exterior]`` two terms move, and they must move together or
the conductor is counted twice:

* **The transferred datum** becomes :math:`g - \psi_c` at each transfer-path
  target, rather than :math:`g`.
* **The transmission condition gains the subtracted conductors' moment on**
  :math:`\Gamma`. The condition is that the *physical* interior normal flux
  match the exterior one. The interior side MEQ assembles is
  :math:`\int_\Gamma (q_h\cdot\nu)\, C_m` over the solved flux, which is
  :math:`q_p`, while the physical flux is :math:`q_p + q_c` — so
  :math:`\int_\Gamma (q_c\cdot\nu)\, C_m` is **added**.

An exterior conductor, one carried outside :math:`\Gamma` by the exterior
expansion itself, enters the same integral with the **opposite** sign: its
field is part of the exterior solution, so it appears on both sides and is
removed from the interior one. The sign is decided by which side of
:math:`\Gamma` the conductor is on, and that makes the geometry a precondition:

* a **subtracted** conductor must lie strictly **inside** :math:`\Gamma`, since
  the exterior Gegenbauer series stands for sources within :math:`\Gamma` and a
  conductor beyond it puts a singularity in the region that series describes;
* an **exterior** conductor must lie strictly **outside** it.

Both are measured to the conductor's **farthest** point, so one straddling
:math:`\Gamma` belongs to neither and is refused by both, whichever is set
first.

In a tokamak this is not a choice between the two. The coils hug the vessel and
the plasma fills it, so at any :math:`\Gamma` large enough to contain the plasma
almost every conductor is inside it, and the exterior expansion can carry only
the few that clear it. Taking the conductors out of the mesh shrinks the
**meshed area inside** :math:`\Gamma` — the conductor elements and the grading
around them — and not :math:`\Gamma` itself, which the plasma's own extent
bounds.

.. _conductors-filaments:

Filaments, and a filament on a mesh point
-----------------------------------------

A point filament has no finite-element representation as a current density, so
**subtraction is the only way MEQ carries one**, and ``"filament"`` is that
route. It is also the cheapest field there is: one elliptic-integral evaluation
per point per filament, against a whole cross-section quadrature per point per
rectangle.

What a filament costs is accuracy **at** the conductor. Far away, one point at
a rectangle's centre is exact to leading order — the collapse conserves the
total current — and near it the line current's logarithm differs from the
rectangle's finite field by as much as the field itself. Measured on a whole
solve, the difference falls by orders of magnitude within a few conductor
widths, so ``"filament"`` is a percent-level answer near the coils and a far
better one away from them. That is a modelling choice to make deliberately, not
a defect to hide, and :doc:`configuration` has the figures for one machine.

**A mesh vertex or node exactly on a filament is refused, and a point near one
is not.** The flux is evaluated, not approximated, and the elliptic-integral
form stays accurate as a point approaches the loop — it holds a constant offset
from :math:`-\tfrac12 \ln\varepsilon` at distances where the textbook form
returns NaN — so a large value at one evaluation point is not a large error. What the
refusal catches is *coincidence*: the mesher putting a node at the same
coordinates, which is reachable by writing the obvious file, since a meshed
coil's centre is where a filament model puts its filament. It is not a
clearance rule, and nothing measured asks for one. MEQ refuses at setup, naming
the filament and the mesh point, rather than moving the node or taking a finite
part, either of which would report a run that looks as if it worked.

.. _conductors-filament-stack:

Dividing a rectangle is the midpoint rule
-----------------------------------------

``[conductors] FilamentSize`` divides each ``[[coils]]`` block into a stack of
filaments, one at the centre of each cell of a grid over the rectangle, each
carrying its share of the current; ``[[coils]] FilamentsR`` and
``FilamentsZ`` set one block's stack directly. :ref:`dividing-a-rectangle` has
the keys.

**The stack is the midpoint rule for the integral the rectangle's own field
is**, and that is what makes it a model rather than a heuristic. It therefore
converges to the rectangle at second order in the cell size, which the suite
asserts as a measured *rate* — near misses such as filaments at cell corners, a
current density where a share belongs, or a transposed index all give small
errors and none of them gives two. The size is a **length** rather than a count
because a cell's error at a field point goes as the square of its size over the
point's standoff, so one length buys a predictable accuracy on conductors of
every aspect ratio.

The convergence holds **clear of the winding only**. Inside it a stack has one
logarithmic singularity per filament where the rectangle has none, and no cell
size removes them.

**The default is one filament per block, and that is deliberate.** Dividing
converges towards the rectangle — and therefore *away* from a code that models
its coils as point filaments at the same centres. Against such a reference the
undivided model is the like-for-like one and dividing makes the agreement worse;
against the winding the reverse holds. Neither setting is the safe one, so a run
prints which it took and the ``.nc`` conductor table lists every filament
individually.

.. _conductors-cache:

The conductor flux is cached at MEQ's own points
------------------------------------------------

The precompute is per mesh, and a run driven by a transport code pays it per
*solve*: the profiles move between solves and the machine does not.
``[conductors] CacheFile`` keeps it on disk (see :ref:`caching-psi-c`), and
the equilibrium ``.nc`` carries the same data as its ``psi_coil`` group
whenever the split is in use.

What is stored is :math:`\psi_c` at every potential degree of freedom and every
source quadrature point, and :math:`q_c` at every node of the flux space — the
second because the critical points are roots of the *physical* flux and the
search tabulates :math:`q_c` at every element node, which is the larger cost of
the two.

**It holds values at the points MEQ uses them, not a grid to interpolate
from.** The familiar alternative samples a coil set's field on a uniform grid
and interpolates it back — mesh-free, reusable at any resolution, and
approximate. MEQ takes the other trade because the split's identities are exact:
a :math:`\psi_c` that came back to within an interpolation error would turn
each of them into a tolerance. So a reload *is* the recompute, every printed
number is identical between a cold run and a warm one, and the cache is
**bound to one mesh, one polynomial degree and one quadrature rule**.

A stale cache is not a slow run but a plausible field belonging to another
machine, so three checks guard it, at the three moments it can go stale:

* **its signature** — the mesh, the conductors, the degree and the element and
  degree-of-freedom counts — before anything is installed, naming the first that
  differs, conductors first;
* **its quadrature points against the rule this run would use**, when the forms
  are built, since raising the quadrature order moves the points while the mesh,
  degree and conductors stay put;
* **its signature again at every solve**, since an adaptive cycle refines the
  mesh without rebuilding the solver.

**A cache never fails a run.** Absent, stale, corrupt or another machine's, it
costs a rebuild, changes no number, and is reported, because a run that has
silently stopped using its cache is a slowdown with no other symptom.

The one thing that is not cached is the **output grid**: the gridded ``.nc``
evaluates :math:`\psi_c` and :math:`B_c` at grid nodes, and the grid is not the
mesh. ``[conductors] InterpolatedFlux`` is the related option for the
critical-point search, and it is off by default because it can change which
equilibrium a cold solve selects — :ref:`interpolated-flux` says when it is
safe.

.. _conductors-costs:

What each model costs and buys
------------------------------

.. list-table::
   :header-rows: 1
   :widths: 16 42 42

   * - Model
     - Buys
     - Costs
   * - ``"meshed"``
     - No per-point field evaluation. The conductor current is in the interior
       equation. The reference the other two are checked against.
     - Conductor elements, and a mesh graded to each conductor's edges. The
       conductor's field resolved to the mesh's order; a source discontinuous
       at the coil edge, which costs the convergence *rate* wherever an element
       straddles it.
   * - ``"subtracted"``
     - The same machine with a smaller mesh: no conductor elements and no
       grading. The conductors' field exact rather than of the mesh's order,
       and the remainder smooth across where the coil edges were.
     - A cross-section quadrature per field point — a large share of a cold
       run on a rectangle machine, which ``CacheFile`` removes from every run
       after the first. The ``.gf`` and the VTK hold the remainder.
   * - ``"filament"``
     - Everything ``"subtracted"`` buys, at one elliptic integral per point per
       filament, so caching it barely matters. The only way to carry a
       point-filament machine, and so the like-for-like comparison with a code
       that uses one.
     - A different machine: a percent-level difference from the rectangles at
       the conductors, falling quickly away from them. A single filament per
       block is a poor stand-in for a conductor long beside the plasma, such as
       a solenoid, which is what ``FilamentSize`` addresses.

The mesh saving is what the split is for, and what it is worth depends on how
much of a machine's mesh its conductors occupy. ``meq-run`` re-meshes when
``Model`` changes, since under a subtracting model the generator is not told
about the conductors at all, so the comparison on your own machine is one key
and two runs.

**What the split does not change** is the plasma edge. A profile whose current
density vanishes at the edge as a fractional power caps the convergence order
there whatever the conductors do, so on such a machine the split buys the mesh
and the exactness of the conductors' field but not a higher order in the
plasma.
