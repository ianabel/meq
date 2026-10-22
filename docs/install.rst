Building MEQ
============

MEQ is a C++17 library and a single executable. It is built with CMake and
depends on a development branch of MFEM; everything else is either vendored, or
optional, or a package your distribution already has.

.. _install-dependencies:

What you need
-------------

.. list-table::
   :header-rows: 1
   :widths: 22 14 64

   * - Dependency
     - Required?
     - Notes
   * - A C++17 compiler
     - yes
     - Developed against recent GCC. No compiler extensions are used.
   * - CMake ≥ 3.20
     - yes
     - Boost is found in ``CONFIG`` mode, because CMake 4 removed
       ``FindBoost``.
   * - **MFEM**, the HDG branch
     - see below
     - The single hard dependency of the solver. ``mfem/master`` **will not
       work**; see :ref:`install-mfem`.
   * - toml11
     - yes
     - Vendored as the git submodule ``extern/toml11``. A clone without
       ``--recursive`` fails at configure time, and the error message says so.
   * - netcdf-cxx4
     - optional
     - Found with ``pkg-config``. Without it the ``.nc`` writer is dropped from
       the library and the driver cannot produce its interchange output. Debian
       calls it ``libnetcdf-c++4-dev``.
   * - Boost
     - **yes**
     - **Header-only, and needed for the library rather than only for the
       tests.** ``Zernike.cpp`` takes its Jacobi polynomials from Boost.Math —
       the Zernike radial polynomial is a Jacobi polynomial under a coordinate
       change, so this replaces a hand-rolled recurrence and supplies an exact
       derivative with it. Nothing is linked.
   * - Eigen
     - **yes**
     - **Header-only, and needed for the library rather than only for the
       tests.** ``SurfaceFit.cpp`` takes its least-squares solve from Eigen's
       ``JacobiSVD``, which is a column-pivoted QR followed by a one-sided
       Jacobi sweep — the same algorithm that file used to carry by hand. The
       swap was made for maintenance rather than for speed; the two are within
       a factor of two of each other on MEQ's own workload. Nothing is linked
       and no MEQ header mentions it. Debian calls it ``libeigen3-dev``.
   * - Boost.Test
     - for the tests
     - ``unit_test_framework``, the shared-library build, which is what supplies
       ``main()`` for each test binary.
   * - clang-tidy
     - optional
     - Gates the ``naming`` test only. Absent, that test is skipped with a
       message rather than failing.
   * - Python: ``numpy``, ``matplotlib``, ``netCDF4``
     - for ``tools/``
     - Only :ref:`plot_equilibrium.py <tools-plot>` needs them; CMake does not
       check for them.

SuiteSparse, SUNDIALS, GSLIB, LAPACK, MKL and CUDA are all reached *through*
MFEM rather than found separately — MEQ reads which of them are present out of
MFEM's own configuration and adapts. :ref:`install-mfem` says which of them MEQ
actually uses and for what.

.. _install-mfem:

MFEM, and why not ``master``
----------------------------

MEQ is built on ``DarcyForm`` and the HDG integrators in MFEM's ``fem/darcy/``
directory, and for curved boundaries on the transfer-path machinery in
``fem/darcy/extension_hdg.*``. Neither is in a released MFEM. The branch is
where the HDG work for this family of solvers lives, and MEQ is an early user of
it.

.. important::

   **What MEQ needs is a local merge of seventeen MFEM development branches, and it
   is published on no remote.** This page assumes you already have it. If you do
   not, ``INSTALL.md`` in the source tree is the account of which branches they
   are, what each one contributes, what breaks silently if one is missing, and
   the unresolved question of where that merge ought to live.

   The most easily lost of the four is the one carrying the post-processing
   fix, because it currently arrives as an *ancestor* of two of the others
   rather than as a merge of its own. Without it :math:`\psi^\star` is wrong
   per element wherever :math:`\partial F/\partial\psi` vanishes — see
   :ref:`postprocessing-singular`.

The MFEM features MEQ reads out of the configuration, and what each buys:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - MFEM option
     - What MEQ does with it
   * - ``MFEM_PRECISION``
     - ``double`` or ``single``, and it decides ``mfem::real_t`` throughout the
       library. **MEQ has no key for this and cannot have one**: it is a
       property of the MFEM an install links, inherited through ``MFEM_DIR``,
       so choosing it means choosing which MFEM to point at. Both build and
       both pass the test suite. See :ref:`install-precision`.
   * - ``MFEM_USE_SUITESPARSE``
     - UMFPACK, one of the direct solvers for the global trace system. **MEQ
       needs this or** ``MFEM_USE_MKL_PARDISO``, and does not configure without
       one of them; there is no iterative fallback. See :doc:`linear_solvers`.
       SuiteSparse ships no single-precision build, so a
       ``MFEM_PRECISION = single`` install has no UMFPACK and MEQ's default
       trace solver there is PARDISO.
   * - ``MFEM_USE_SUNDIALS``
     - KINSOL, which backs the Anderson–Picard and Picard-then-Newton
       globalisations. See :doc:`nonlinear`.
   * - ``MFEM_USE_GSLIB``
     - ``FindPointsGSLIB``. Optional: MEQ interpolates a warm start from one
       mesh onto another with its own point locator, which reproduces gslib to
       every printed digit and has no floating-point abstraction to lose in a
       single-precision build.
   * - ``MFEM_USE_EXCEPTIONS``
     - Makes an MFEM error a C++ exception rather than an abort. **This is what
       makes the driver's exit code 2 reachable**: without it a failed solve
       takes the process down with ``SIGABRT`` instead of reporting itself.
   * - ``MFEM_USE_MKL_PARDISO``
     - A trace solver, selectable at run time, and the default wherever it is
       present.
   * - ``MFEM_USE_CUDA``, ``MFEM_USE_CUDSS``
     - The device path. It buys three things: the ``cuDSS`` trace solver,
       MFEM's batched element and face kernels — which are reached through
       ``[solver] AssemblyMode``, ``LocalFactorMode`` and ``TraceAssemblyMode``
       and are the device route — and MEQ's own device translation units, which
       CMake hands to ``nvcc`` at the architecture it reads out of MFEM's
       ``config.mk``. Without CUDA those compile to host stubs and everything
       still builds and runs.
   * - ``MFEM_USE_OPENMP`` + ``MFEM_THREAD_SAFE``
     - Both together enable threaded element assembly, which MEQ exposes as an
       opt-in. Threaded assembly *aborts* rather than falling back if only one
       of the two is set.
   * - ``MFEM_USE_MPI``
     - MEQ is serial and does not use it. It is detected only so that a
       parallel MFEM fails intelligibly rather than at link time.

.. _install-precision:

Choosing the precision, which means choosing an MFEM
----------------------------------------------------

``MFEM_PRECISION`` decides ``mfem::real_t`` for the whole library, so it is
settled when MFEM is configured and is not something a MEQ run can ask for.
There is no ``Precision`` key in the configuration file, and there will not be.
Building MEQ at one precision or the other is therefore a matter of which MFEM
``MFEM_DIR`` points at:

.. code-block:: sh

   cmake -S source -B builds/double -DMFEM_DIR=/path/to/mfem-double/install
   cmake -S source -B builds/single -DMFEM_DIR=/path/to/mfem-single/install

Each installs to its own prefix, ``install/double`` and ``install/single``.

Both work and both pass the test suite. What differs:

* **Double is the reference.** Every published measurement is taken there, and
  it is the only build in which MEQ's own acceptance criteria — measured
  convergence rates, :math:`k+1` in :math:`\psi` and :math:`k+2` in
  :math:`\psi^\star` — can be read at all: at a round-off floor near
  :math:`10^{-7}` there is no mesh sequence left above it to measure a rate on.
  A single build is validated by *agreement with the double build*, not by its
  own rate tables.
* **Double keeps two libraries single cannot have.** SuiteSparse and GSLIB have
  no single-precision build. Neither is required — MEQ's default trace solver
  is PARDISO where it is present, and the warm start uses MEQ's own point
  locator — but a single install is not a double install minus one flag, and a
  comparison between the two has to account for that.
* **Single is for a device.** The case for it is arithmetic throughput on a
  part whose double rate is a small fraction of its single one; on a
  workstation card that is most of them, and on a datacentre part it is about
  half, which is why a **double** device build is the right one there. The
  device path itself is identical in both: it keys on whether an
  ``mfem::Device`` is configured, never on ``real_t``.
* **Accuracy is qualified per problem, not per build.** Measured across the
  shipped machines, a single-precision equilibrium's ``psi_ax`` agrees with the
  double build's over a range of three orders of magnitude — so the question
  *is single precision good enough* has an answer per case and no answer per
  install.

Whichever is built, ``meq --version`` and every output file's
``meq_precision`` attribute say which one produced a given answer.

.. _install-layout:

The project directory
---------------------

The git checkout is one of four directories in a *project directory*, and
nothing is built, tested or run inside it:

.. code-block:: text

   <project>/source/            the git checkout
   <project>/builds/<name>/     one CMake build tree per configuration
   <project>/install/<name>/    what cmake --install builds/<name> puts there
   <project>/runs/              everything a run of meq or of tools/ writes

The last three are untracked and can be deleted whenever their contents are no
longer wanted. **CMake refuses a build tree inside the source tree**, rather
than leaving one to be gitignored, so untracked output in ``git status`` is a
bug in whatever wrote it.

Building
--------

From the project directory:

.. code-block:: sh

   git -C source submodule update --init --recursive     # extern/toml11
   cmake -S source -B builds/double
   cmake --build builds/double -j4
   cmake --install builds/double

The default install prefix is ``install/<name>`` for ``builds/<name>``, so the
last line puts ``meq``, ``meq-convert`` and ``meq-run`` in ``install/double/bin``, and the mesh
generator ``meq-run`` calls under ``install/double/share/meq/tools``. The
installed ``meq-run`` finds both relative to itself, so the build tree can be
deleted once it is installed.

``MFEM_DIR`` is how the MFEM install is found. It is a cache variable, it also
reads the environment, and it defaults to ``mfem/install`` beside the project
directory — ``../../mfem/install`` from ``source/`` — so on a machine laid out
that way the ``-D`` is unnecessary:

.. code-block:: sh

   cmake -S source -B builds/double -DMFEM_DIR=/path/to/mfem/install

The finder accepts either layout an MFEM build leaves behind — an installed tree
with ``include/``, ``lib/`` and ``share/mfem/config.mk``, or an in-source make
build with ``config/config.mk`` — and needs no flag to tell them apart. It reads
the enabled features out of ``config.mk`` rather than probing for them, so the
configure summary reports exactly what MFEM was built with.

.. warning::

   **Never run a bare** ``make -j`` **or** ``cmake --build -j`` **with no
   number.** With no argument the job count is unbounded and goes to the host's
   core count, which on a container or a WSL2 virtual machine is not the same
   thing as the memory available to it. Give a number: 4 to 8 for MEQ, and 4 for
   the MFEM tree, whose translation units are large.

Configure-time options
----------------------

.. list-table::
   :header-rows: 1
   :widths: 26 12 62

   * - Option
     - Default
     - Meaning
   * - ``MEQ_BUILD_APP``
     - ``ON``
     - Build the ``meq`` executable. A fatal error if it is on and the library
       has no sources — which is what happens without MFEM.
   * - ``MEQ_BUILD_TESTS``
     - ``ON``
     - Build the test binaries and register them with CTest.
   * - ``MEQ_ENABLE_COVERAGE``
     - ``OFF``
     - Instrument with ``--coverage -O0 -g``. See :ref:`testing-coverage` for
       why this is an option rather than a build type.
   * - ``MFEM_DIR``
     - ``../../mfem/install`` from ``source/``
     - Where to find MFEM. Cache variable, environment variable, or ``-D``.
   * - ``TOML11_DIR``
     - ``extern/toml11``
     - Where to find toml11. Cache variable or environment variable.
   * - ``CMAKE_BUILD_TYPE``
     - ``Release``
     - Forced to ``Release`` for single-configuration generators when unset.

.. note::

   **A MEQ debug build buys you the debugger and nothing else** — and it is
   worth saying so, because the opposite was written down here and in the
   project's own notes for some time.

   MFEM's internal consistency checks are ``MFEM_ASSERT``, which is gated on
   ``MFEM_DEBUG`` — **not** on ``NDEBUG``, and nothing derives one from the
   other. A normal MFEM install sets ``MFEM_DEBUG = NO``, so every one of those
   assertions is already dead in every MEQ build, release or debug, and no flag
   passed to MEQ can revive one. MEQ's own sources contain no ``assert()``
   either, so ``-DNDEBUG`` disables nothing of MEQ's.

   Recovering those checks means building a **second MFEM** configured
   ``MFEM_DEBUG=YES`` and pointing ``MFEM_DIR`` at it. That is a real option
   when chasing a suspected contract violation between MEQ and the library, and
   it is a much larger undertaking than a build type.

Building without MFEM
---------------------

``find_package(MFEM)`` is deliberately **not** ``REQUIRED``. Without it, CMake
builds the MFEM-free half of the library — the configuration parser, the
profiles, the source terms, the boundary shapes — skips every file that includes
``mfem.hpp``, and skips every test that needs a solver, printing a status line
for each thing it left out.

This is not a convenience. It is what makes continuous integration possible at
all: the MFEM branch MEQ needs is not published anywhere a hosted runner can
fetch it, so CI builds and tests only that half.

.. important::

   **Read a green CI badge as evidence about the configuration and profile
   layers and about nothing else.** The HDG assembly, Newton, the extension
   technique, the estimator and every convergence claim are exercised only by a
   local build against a real MFEM. :doc:`testing` says what that covers.

Checking the build
------------------

.. code-block:: sh

   cd builds/double && ctest --output-on-failure

The full suite takes several minutes, runs in the build tree, and writes nothing
into ``source/``. Nothing needs to be set in the environment: the test
registration puts ``MKL_NUM_THREADS=1`` on every test itself, which is the one
variable MEQ needs and a no-op on a machine without MKL. Running a test binary
by hand needs the same one thing, from a directory it may write into:

.. code-block:: sh

   cd builds/double/tests && MKL_NUM_THREADS=1 ./SolovievConvergence

.. warning::

   That variable is not cosmetic. MEQ's hot inner loop factorises a great many
   *small* dense matrices — one or two per element, per residual evaluation —
   and on a threaded MKL each of those calls pays for a fork and a barrier that
   dwarfs the arithmetic. The effect was measured on the development machine and
   it was large, growing rapidly with polynomial degree. It is a property of
   where the threshold sits in your BLAS, so **measure it on yours** before
   concluding it does not apply; the safe setting is 1.

   This is in tension with the threaded trace solvers, which want MKL threads,
   and the tension is real rather than an oversight — ``MKL_NUM_THREADS`` is
   process-wide. :doc:`linear_solvers` discusses it.

Building the documentation
--------------------------

These pages are Sphinx, and their dependencies are pinned in
``docs/requirements.txt``:

.. code-block:: sh

   cd source
   python3 -m venv .venv-docs
   .venv-docs/bin/pip install -r docs/requirements.txt
   .venv-docs/bin/sphinx-build -W -b html docs ../builds/docs/html

or, equivalently, ``cmake --build builds/double --target docs``, which creates
the same virtual environment if it is missing and writes to
``builds/double/docs/html``. The ``-W`` is deliberate and matches
what Read the Docs runs, so a broken cross-reference fails the build here rather
than publishing a page with a dead link.
