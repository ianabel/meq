Building meq
============

meq is a C++17 library and a single executable. It is built with CMake and
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
MFEM rather than found separately — meq reads which of them are present out of
MFEM's own configuration and adapts. :ref:`install-mfem` says which of them meq
actually uses and for what.

.. _install-mfem:

MFEM, and why not ``master``
----------------------------

meq is built on ``DarcyForm`` and the HDG integrators in MFEM's ``fem/darcy/``
directory, and for curved boundaries on the transfer-path machinery in
``fem/darcy/extension_hdg.*``. Neither is in a released MFEM. The branch is
where the HDG work for this family of solvers lives, and meq is an early user of
it.

The MFEM features meq reads out of the configuration, and what each buys:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - MFEM option
     - What meq does with it
   * - ``MFEM_USE_SUITESPARSE``
     - UMFPACK, the default direct solver for the global trace system. Without
       any direct solver meq falls back to GMRES; see :doc:`linear_solvers`.
   * - ``MFEM_USE_SUNDIALS``
     - KINSOL, which backs the Anderson–Picard and Picard-then-Newton
       globalisations. See :doc:`nonlinear`.
   * - ``MFEM_USE_GSLIB``
     - ``FindPointsGSLIB``, which is how a warm start is interpolated from one
       mesh onto another.
   * - ``MFEM_USE_EXCEPTIONS``
     - Makes an MFEM error a C++ exception rather than an abort. **This is what
       makes the driver's exit code 2 reachable**: without it a failed solve
       takes the process down with ``SIGABRT`` instead of reporting itself.
   * - ``MFEM_USE_MKL_PARDISO``
     - An alternative trace solver, selectable at run time.
   * - ``MFEM_USE_CUDA``, ``MFEM_USE_CUDSS``
     - A third trace solver, ``cuDSS``. A prerequisite rather than a speed-up:
       ``fem/darcy/`` has no partial-assembly kernels and meq issues no device
       kernels of its own.
   * - ``MFEM_USE_OPENMP`` + ``MFEM_THREAD_SAFE``
     - Both together enable threaded element assembly, which meq exposes as an
       opt-in. Threaded assembly *aborts* rather than falling back if only one
       of the two is set.
   * - ``MFEM_USE_MPI``
     - meq is serial and does not use it. It is detected only so that a
       parallel MFEM fails intelligibly rather than at link time.

Building
--------

.. code-block:: sh

   git submodule update --init --recursive     # extern/toml11
   cmake -B build
   cmake --build build -j4

``MFEM_DIR`` is how the MFEM install is found. It is a cache variable, it also
reads the environment, and it defaults to ``../mfem/install`` relative to the
source tree, so on a machine laid out that way the ``-D`` is unnecessary:

.. code-block:: sh

   cmake -B build -DMFEM_DIR=/path/to/mfem/install

The finder accepts either layout an MFEM build leaves behind — an installed tree
with ``include/``, ``lib/`` and ``share/mfem/config.mk``, or an in-source make
build with ``config/config.mk`` — and needs no flag to tell them apart. It reads
the enabled features out of ``config.mk`` rather than probing for them, so the
configure summary reports exactly what MFEM was built with.

.. warning::

   **Never run a bare** ``make -j`` **or** ``cmake --build -j`` **with no
   number.** With no argument the job count is unbounded and goes to the host's
   core count, which on a container or a WSL2 virtual machine is not the same
   thing as the memory available to it. Give a number: 4 to 8 for meq, and 4 for
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
     - ``../mfem/install``
     - Where to find MFEM. Cache variable, environment variable, or ``-D``.
   * - ``TOML11_DIR``
     - ``extern/toml11``
     - Where to find toml11. Cache variable or environment variable.
   * - ``CMAKE_BUILD_TYPE``
     - ``Release``
     - Forced to ``Release`` for single-configuration generators when unset.

.. note::

   **A debug build is worth doing occasionally**, and not only for the debugger.
   Several of the contracts between meq and MFEM's block structures are checked
   by ``MFEM_ASSERT``, which is compiled out under ``NDEBUG``. At least one real
   defect in meq — a block vector passed with the wrong number of blocks — was
   invisible in a release build and caught immediately by an assert.

Building without MFEM
---------------------

``find_package(MFEM)`` is deliberately **not** ``REQUIRED``. Without it, CMake
builds the MFEM-free half of the library — the configuration parser, the
profiles, the source terms, the boundary shapes — skips every file that includes
``mfem.hpp``, and skips every test that needs a solver, printing a status line
for each thing it left out.

This is not a convenience. It is what makes continuous integration possible at
all: the MFEM branch meq needs is not published anywhere a hosted runner can
fetch it, so CI builds and tests only that half.

.. important::

   **Read a green CI badge as evidence about the configuration and profile
   layers and about nothing else.** The HDG assembly, Newton, the extension
   technique, the estimator and every convergence claim are exercised only by a
   local build against a real MFEM. :doc:`testing` says what that covers.

Checking the build
------------------

.. code-block:: sh

   cd build && ctest --output-on-failure

The full suite takes several minutes. Nothing needs to be set in the
environment: the test registration puts ``MKL_NUM_THREADS=1`` on every test
itself, which is the one variable meq needs and a no-op on a machine without
MKL. Running a test binary by hand needs the same one thing:

.. code-block:: sh

   MKL_NUM_THREADS=1 ./tests/SolovievConvergence

.. warning::

   That variable is not cosmetic. meq's hot inner loop factorises a great many
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

   python3 -m venv /tmp/docsvenv
   /tmp/docsvenv/bin/pip install -r docs/requirements.txt
   /tmp/docsvenv/bin/sphinx-build -W -b html docs docs/_build/html

or, equivalently, ``make -C docs html``. The ``-W`` is deliberate and matches
what Read the Docs runs, so a broken cross-reference fails the build here rather
than publishing a page with a dead link.
