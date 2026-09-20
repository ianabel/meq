# Unit tests for meq_core.
#
# One executable per test file, each its own Boost.Test module, each registered
# with ctest. Keeping them separate means a test that crashes the process --
# which, in a library that calls into MFEM, is a realistic failure mode -- takes
# down only its own suite and is named precisely by ctest.

# CONFIG mode: FindBoost was removed in CMake 4, so this relies on the
# BoostConfig.cmake that Boost itself installs (1.90 here). The imported target
# defines BOOST_UNIT_TEST_FRAMEWORK_DYN_LINK when it resolves to the shared
# library, which is what makes Boost.Test emit its own main(); without it the
# link fails with an undefined reference to main.
find_package(Boost CONFIG REQUIRED COMPONENTS unit_test_framework)

# Each entry is the base name of a file in unit/, and by convention tests the
# unit of the same name: ConfigTests.cpp tests src/meq/Config.cpp.
set(MEQ_UNIT_TESTS
	BoundaryShapeTests
	CoilsTests
	ConductorFieldTests
	ConfigTests
	ExteriorDtNTests
	FluxFamilyTests
	PlasmaComponentTests
	ProfilesTests
	RotatingSourceTests
	SafetyFactorTests
	SourceTests
	SourceFactoryTests
	ZernikeTests)

foreach(_test IN LISTS MEQ_UNIT_TESTS)
	set(_test_source "${CMAKE_CURRENT_SOURCE_DIR}/unit/${_test}.cpp")

	# The test files are being written alongside the units they cover, so
	# both halves have to be there before there is anything to build: a test
	# whose unit is missing from meq_core would configure happily and then
	# fail at link with undefined references.
	string(REGEX REPLACE "Tests$" "" _unit "${_test}")
	list(FIND MEQ_CORE_SOURCES_PRESENT "src/meq/${_unit}.cpp" _unit_index)

	if(NOT EXISTS "${_test_source}")
		message(STATUS "meq: no ${_test}.cpp yet, test not registered")
	elseif(NOT TARGET meq_core OR _unit_index EQUAL -1)
		message(STATUS "meq: ${_test}.cpp exists but src/meq/${_unit}.cpp does not, test not registered")
	else()
		add_executable(${_test} "unit/${_test}.cpp")
		target_link_libraries(${_test} PRIVATE
			meq_core
			Boost::unit_test_framework)

		# So that a unit test can reach the closed-form solutions in
		# analytic/. SourceFactoryTests uses one to check that the library's
		# copy of the Example 5 source has not drifted from the fixture's.
		target_include_directories(${_test} PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")

		# WORKING_DIRECTORY is load-bearing, not tidiness: ConfigTests opens
		# examples/*.toml by a path relative to the source tree, so without
		# this it passes from build/ and fails from anywhere else. Pinning it
		# makes an out-of-tree build behave the same as an in-tree one.
		add_test(NAME ${_test} COMMAND ${_test}
			WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}")

		# MKL_THREADING_LAYER=GNU USED TO BE SET HERE AND IS DELIBERATELY
		# NOT ANY MORE. It guarded a real trap: libblas.so.3 on the
		# development machine was an alternatives symlink to libmkl_rt.so,
		# the MKL dispatcher, whose default threading layer silently
		# returned wrong BLAS-3 results to UMFPACK. meq no longer loads
		# that dispatcher -- it builds against its own SuiteSparse, which
		# links oneAPI's threading layer directly -- so the variable became
		# inert, measured as bit-identical output across three
		# value-asserting suites with it set and unset.
		#
		# It is dropped rather than kept as a cheap guard because MOST
		# MACHINES WILL NOT HAVE MKL AT ALL, and a variable naming a
		# library that is not there is a false instruction to whoever reads
		# it next. Choosing a threading layer belongs in CMake, which can
		# look at what the build actually links, and that is where it goes
		# when meq runs on more than one machine.
		#
		# MKL_NUM_THREADS=1 STAYS, and for a completely different reason:
		# not correctness but a factor of 140. It is a no-op wherever MKL
		# is absent, so it costs a non-MKL machine nothing. UMFPACK
		# calls BLAS on the small dense frontal matrices of a sparse LU,
		# thousands of times, and each threaded call costs a fork, a join
		# and a barrier that dwarf the arithmetic. Measured on this
		# machine, one process, nothing else running:
		#
		#     SolovievConvergence   MKL=1    1.24 s
		#     SolovievConvergence   MKL=16 177.83 s   at 273% CPU
		#
		# 273% of sixteen cores is the shape of the answer -- the threads
		# are spinning on barriers, not working. Tests that never reach
		# UMFPACK are untouched (SamplerConvergence: 0.64 s against
		# 0.63 s), which is what identifies the cause.
		#
		# This did not use to bite, and the reason is worth knowing: the
		# link line used to carry libmkl_sequential behind SuiteSparse,
		# so MKL resolved to the sequential layer whatever it was asked
		# for. Pointing MFEM's PARDISO at oneAPI's threaded layer removed
		# it and revealed the real behaviour. The suite was never
		# deliberately sequential; it was accidentally so.
		#
		# The culprit there is MFEM's element-local dense LU rather than
		# UMFPACK, which degrades by only about 40% across the whole
		# thread range. See CLAUDE.md under Traps.
		#
		# It is set per test rather than expected of whoever runs ctest,
		# because a benchmark or a suite that takes two hours instead of
		# one minute is not a failure anybody reads as a threading
		# problem.
		set_tests_properties(${_test} PROPERTIES
			ENVIRONMENT "MKL_NUM_THREADS=1")
	endif()
endforeach()

# ---------------------------------------------------------------------------
# convergence -- a measured rate, not a tolerance
# ---------------------------------------------------------------------------

# The acceptance criterion for each stage of the port. These assert an order of
# accuracy against a closed-form solution, which is what catches the failure the
# unit tests above cannot see: a sign or weight convention that converges, at the
# right rate, to the wrong function. See CLAUDE.md, "Testing stance".
#
# They are registered the same way as the unit tests -- one executable, one
# Boost.Test module, one ctest -- and differ only in living under convergence/
# and in needing tests/ on the include path, since they pull the closed-form
# solutions out of analytic/.
set(MEQ_CONVERGENCE_TESTS
	ConductorSubtraction
	FieldConvergence
	OutputConvergence
	SamplerConvergence
	SolovievConvergence
	NewtonConvergence
	SolverContract
	TraceSolverComparison
	WarmStartConvergence
	HighBetaConvergence
	ExtensionConvergence
	EstimatorConvergence
	AdaptiveRefinement
	PedestalConvergence
	MillerConvergence
	SolovievGeometryConvergence
	RotatingNormalisedConvergence
	RotatingSolovievConvergence
	RotatingSourceConvergence
	RotatingNewtonConvergence
	# Maschke & Perrin's section 4, the second exact rotating equilibrium. It
	# solves the SAME equation as the rotating Solov'ev case above, so nothing
	# about the discretisation is new; what is new is that its (4.7) constrains
	# only omega^2/T, leaving both an arbitrary function of psi. It is therefore
	# the only closed form MEQ has that can drive meq::RotatingSource with
	# psi-DEPENDENT temperature and rotation, where the alternative check is a
	# central difference of MEQ's own f().
	MaschkePerrinConvergence
	# The one thing the flow campaign was missing: meq::RotatingSource at
	# C'( psi ) != 0 against an INDEPENDENT closed form. C is the exponent both
	# species share and C' is non-zero only where omega^2/T varies from surface
	# to surface -- which is the term Li & Zhu print with the wrong sign, and
	# which no exact rotating solution can reach, C constant being precisely what
	# makes the equation solvable in closed form. So the fixture is a manufactured
	# one, tests/analytic/VaryingCentrifugal.hpp, deriving (96), (97) and (136)
	# from the paper and calling no MEQ code. It carries a bisection route of its
	# own beside the algebra, so that two transcriptions agreeing is not the whole
	# of the evidence.
	VaryingCentrifugalConvergence
	CriticalPointConvergence
	FluxSurfaceConvergence
	FluxGridConvergence
	SurfaceAverageConvergence
	SurfaceFitConvergence
	# FB-A, the first stage of FREE-BOUNDARY-PLAN.md and the only one that needs
	# no free boundary: a vacuum solve on a mesh whose inner edge IS the axis,
	# where the flux mass ( r q, v ) degenerates and the operator's 1/r is not
	# integrable. It is an ordinary convergence test because it has an ordinary
	# exact solution -- see tests/analytic/VacuumHarmonic.hpp for why one exists
	# at all, which the plan did not know.
	AxisConvergence
	# XP-1, the connectivity of the confined plasma support. It is a convergence
	# test rather than a unit one because the flood fill's INPUT is a mesh and a
	# solved field; the fill itself is MFEM-free and is unit tested as
	# PlasmaComponentTests.
	PlasmaConnectivity
	# XP-2: psi_bnd from the located X-point, as an outer fixed point.
	XPointOuter
	# XP-3: the same statement made INSIDE the Newton -- the X-point is two more
	# unknowns and q_r = q_z = 0 are two more rows. It shares XPointOuter's
	# fixture through convergence/DivertedMachine.hpp, deliberately: its headline
	# assertion is agreement with that case, so the two must solve one machine.
	XPointBorder
	# THE BORDER'S JACOBIAN AGAINST A DIFFERENCE OF ITS CONSTRAINTS, which is the
	# check NewtonConvergence.cpp has always made on the FIELD block and which
	# the border has never had. M-118 records two defects one would have caught.
	# What it settles is the envelope assumption the axis row rests on, in the
	# one form that is decidable from outside the assembly: the row is supported
	# on potential dofs, so it asserts that the constraint cannot respond to the
	# FLUX -- and a flux perturbation moves the axis while leaving psi_h alone.
	BorderJacobian
	# THREADING-PLAN.md items A and B, which are one case each and share a
	# fixture. Item A compresses the eleven border rows that walked the whole
	# 109,200-entry vector onto the indices they are non-zero on; item B threads
	# the five plasma element loops that build them. Both claim EXACTNESS rather
	# than agreement -- dropping a term that is exactly zero cannot move a sum,
	# and an L2 scatter is disjoint by construction -- so both assert at
	# 0.000e+00 rather than at a tolerance. It calls the five assemblers
	# directly, which is what they were made public for: a race that does not
	# fire on this mesh at this thread count leaves a whole-solve pin looking
	# perfect.
	BorderAssembly
	# [solver] UpDownSymmetry: the dof-for-dof projection onto fields even in z,
	# which exists so a DOUBLE NULL need not be told which of its two saddles to
	# follow. Two cases, and the second is the one MAST-U meets: the projection
	# is the identity on a problem that already has the symmetry -- the flux's
	# per-component sign being what that can fail on -- and it REFUSES a mesh
	# whose dofs are not mirror-paired, which gmsh's triangulation of a
	# perfectly symmetric machine is not.
	UpDownSymmetry
	# FB-1's coupling matrix P: the trace projection of the exterior modes onto
	# Gamma_h. It needs a curved boundary and so needs MFEM, which is why it is
	# a convergence test rather than a unit one -- ExteriorDtN's own tests, which
	# CI can run, stop at the operator.
	FreeBoundaryCoupling
	# FB-4: what the plasma edge costs the order. The source stops dead on a
	# curve no mesh line follows and that moves while Newton runs, and the
	# answer is that the PROFILE decides the rate rather than the quadrature --
	# psi* keeps k+2 exactly when p' vanishes at the edge to order j >= k. It
	# also carries the measurement that says why no cut quadrature was built.
	PlasmaEdgeConvergence
	# psi_ax is the largest NODAL value of psi_h, and nothing in that definition
	# says the largest nodal value is a magnetic axis. This is the guard --
	# CriticalPointFinder::checkAxis(), which reads the normalised flux at a zero
	# of q_h and must find 1 there -- and it is a convergence test rather than a
	# unit one because half of it drives the bordered Newton.
	AxisAgreement
	# FB-6's remaining capability gap: the limiter as a CURVE rather than a
	# prescribed point, restricted to the POLYGON the mesh is fitted to. A
	# convergence test rather than a unit one because the contact is found on a
	# mesh and the constraint is closed by the bordered Newton -- and because the
	# assertion with teeth is Newton's own observed ORDER, which is the only
	# thing that can see whether the border row really is exact.
	LimiterCurve
	# ROADMAP.md item 10: driving the equilibrium by q( psi ) rather than by
	# g( psi ). A convergence test because the loop is solve -> extract ->
	# invert -> re-solve, and what it asserts is that the fixed point is the
	# equilibrium whose target q was built from a KNOWN g.
	# INVERSION-PLAN.md stage IN-5: surfaces that do not close. A convergence
	# test because what it asserts is that a Chebyshev series in normalised arc
	# length converges on an open curve while a periodic basis cannot.
	OpenSurfaces
	SafetyFactorSolver)

# DriverAcceptance is registered separately below: it runs the meq BINARY rather
# than linking the library, so it needs the app to have been built and needs the
# path to it, and it must not be registered when MEQ_BUILD_APP is OFF.

foreach(_test IN LISTS MEQ_CONVERGENCE_TESTS)
	set(_test_source "${CMAKE_CURRENT_SOURCE_DIR}/convergence/${_test}.cpp")

	# Every convergence test drives the solver, so there is nothing to build
	# until GradShafranov is in meq_core.
	list(FIND MEQ_CORE_SOURCES_PRESENT "src/meq/GradShafranov.cpp" _solver_index)

	if(NOT EXISTS "${_test_source}")
		message(STATUS "meq: no ${_test}.cpp yet, test not registered")
	elseif(NOT TARGET meq_core OR _solver_index EQUAL -1)
		message(STATUS "meq: ${_test}.cpp exists but src/meq/GradShafranov.cpp does not, test not registered")
	else()
		add_executable(${_test} "convergence/${_test}.cpp")

		# For "analytic/Soloviev.hpp" and
		# "analytic/ManufacturedNonlinear.hpp". meq_core already carries
		# src/ and src/meq/, so "meq/GradShafranov.hpp" needs nothing here.
		target_include_directories(${_test} PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")

		target_link_libraries(${_test} PRIVATE
			meq_core
			Boost::unit_test_framework)

		# WORKING_DIRECTORY is load-bearing, not tidiness: ConfigTests opens
		# examples/*.toml by a path relative to the source tree, so without
		# this it passes from build/ and fails from anywhere else. Pinning it
		# makes an out-of-tree build behave the same as an in-tree one.
		add_test(NAME ${_test} COMMAND ${_test}
			WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}")

		# Not optional, and least optional here of anywhere in the tree: the
		# hybridized trace system is solved by UMFPACK, whose BLAS-3 calls MKL
		# silently corrupts under its default threading layer. The symptom is a
		# convergence table full of plausible numbers that are wrong, which is
		# exactly what this test exists to detect and exactly what it would then
		# fail to detect. See the comment on the unit tests above.
		set_tests_properties(${_test} PROPERTIES
			ENVIRONMENT "MKL_NUM_THREADS=1")
	endif()
endforeach()

# ---------------------------------------------------------------------------
# The performance harness. BUILT BUT NOT REGISTERED AS A TEST, deliberately.
#
# Everything it reports is a timing, and this project's standing rule is that a
# threaded timing on this machine is a measurement about the machine rather than
# about the code. A ctest that fails when the machine is busy teaches people to
# ignore the suite. It returns non-zero only for the two CORRECTNESS properties
# that make its timings meaningful -- threaded assembly reproducing serial
# assembly bit for bit, and PARDISO reproducing UMFPACK -- so it is still usable
# in CI as a guard, just not as a benchmark gate.
#
# Run it through tests/performance/scan.sh.
# ---------------------------------------------------------------------------

list(FIND MEQ_CORE_SOURCES_PRESENT "src/meq/GradShafranov.cpp" _solver_index)
if(TARGET meq_core AND NOT _solver_index EQUAL -1
   AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/performance/TraceSolverScaling.cpp")
	add_executable(TraceSolverScaling "performance/TraceSolverScaling.cpp")
	target_include_directories(TraceSolverScaling PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
	target_link_libraries(TraceSolverScaling PRIVATE meq_core)

	# ONE ctest out of the performance harness, and it is not a benchmark gate.
	#
	# cuDSS agreement is the one thing about the device path that CAN be asserted
	# on this machine: its timings vary by a factor of thirty run to run, but it
	# either reaches UMFPACK's answer or it does not. SolverContract cannot check
	# it, because cuDSS reads its data through MFEM's device-aware accessors and
	# so needs an mfem::Device, which is global state that has to be configured
	# before any other MFEM object exists -- doing that inside a Boost test case
	# would change how every other test in that binary allocates.
	#
	# So the check runs here instead, on the smallest sizes that exercise the
	# path, with --repeats 1 because this is a correctness run and the timings it
	# prints are noise. The binary exits 0 with a "skipped" line when the build
	# has no CUDA, so this passes on a machine with no GPU rather than failing
	# for not having one.
	add_test(NAME cuDSSTraceSolver
		COMMAND TraceSolverScaling --device cuda --orders 2 --sizes 8,16
		        --no-assembly-timing --repeats 1
		WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}")
	set_tests_properties(cuDSSTraceSolver PROPERTIES
		ENVIRONMENT "MKL_NUM_THREADS=1")
endif()

# InversionScaling -- the same rules, for the solution-inversion chain.
#
# NOT a ctest, and unlike TraceSolverScaling it contributes no ctest either:
# every property it asserts on is already asserted by FluxSurfaceConvergence,
# CriticalPointConvergence, SurfaceFitConvergence and SurfaceAverageConvergence
# on the serial path, and the ones that are NOT -- threaded extraction
# reproducing serial extraction bit for bit, and the candidate kernels
# reproducing the code they would replace -- are properties of a THREAD COUNT
# and of a build, which is exactly the kind of thing this project keeps out of
# the suite. It exits non-zero for them all the same, so it is usable as a
# guard by whoever runs it.
#
# Eigen is OPTIONAL and header only. It is used by section 3 alone, to PRICE
# the least-squares solve that src/meq/SurfaceFit.cpp hand-rolls -- nothing in
# src/meq includes it and nothing in meq depends on it, so a build without it
# prints "Eigen is not available" and measures everything else. If it ever
# earns its place in the library, INSTALL.md and the install chapter of docs/
# both enumerate dependencies and would have to say so.
if(TARGET meq_core AND NOT _solver_index EQUAL -1
   AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/performance/InversionScaling.cpp")
	add_executable(InversionScaling "performance/InversionScaling.cpp")
	target_include_directories(InversionScaling PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
	target_link_libraries(InversionScaling PRIVATE meq_core)

	find_package(Eigen3 3.3 QUIET NO_MODULE)
	if(TARGET Eigen3::Eigen)
		target_link_libraries(InversionScaling PRIVATE Eigen3::Eigen)
		target_compile_definitions(InversionScaling PRIVATE MEQ_HAVE_EIGEN)
		message(STATUS "meq: Eigen ${Eigen3_VERSION} found, InversionScaling will price the least-squares solve against it")
	else()
		message(STATUS "meq: Eigen not found, InversionScaling skips its least-squares comparison")
	endif()

	# OpenMP for section 4, which is the only place in meq that opens a parallel
	# region of its own. MFEM's own -fopenmp comes through MFEM_COMPILE_OPTIONS
	# already -- see cmake/FindMFEM.cmake -- but asking for it here is what makes
	# the dependency explicit rather than inherited.
	find_package(OpenMP QUIET COMPONENTS CXX)
	if(TARGET OpenMP::OpenMP_CXX)
		target_link_libraries(InversionScaling PRIVATE OpenMP::OpenMP_CXX)
	endif()
endif()

# NpcThreadScaling -- what AssemblyMode::Threaded is worth on a NONLINEAR solve,
# and whether PARDISO's MKL threads are reachable at the same time.
#
# NOT a ctest, and it contributes none, for the standing reason: every number in
# it is a timing. It exits non-zero only for the correctness properties that make
# those timings mean anything, so it is usable as a guard by whoever runs it.
#
# It exists BESIDE TraceSolverScaling rather than inside it because that harness
# times a LINEAR solve, and MFEM has since threaded MultNL() -- the residual and
# the Jacobian assembly, and so every NPC step. A linear problem never calls
# MultNL at all, so the option's value has moved to a loop the older harness
# structurally cannot reach. The bit-identity check here is likewise the
# nonlinear counterpart of SolverContract's, and it is the one that sees
# meq::SourceIntegrator.
if(TARGET meq_core AND NOT _solver_index EQUAL -1
   AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/performance/NpcThreadScaling.cpp")
	add_executable(NpcThreadScaling "performance/NpcThreadScaling.cpp")
	target_include_directories(NpcThreadScaling PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
	target_link_libraries(NpcThreadScaling PRIVATE meq_core)

	# omp_get_max_threads() only, to report the thread count it is running at --
	# a threaded timing that does not say how many threads it used is not a
	# measurement. MFEM's own -fopenmp arrives through MFEM_COMPILE_OPTIONS
	# already; asking here makes the dependency explicit rather than inherited.
	find_package(OpenMP QUIET COMPONENTS CXX)
	if(TARGET OpenMP::OpenMP_CXX)
		target_link_libraries(NpcThreadScaling PRIVATE OpenMP::OpenMP_CXX)
	endif()
endif()

# NewtonStepProfile -- WHERE a Newton step's time goes, split into legs that do
# not overlap. The measurement HDG-NEWTON-STEP-PROFILE-FROM-HDGDEV.md asks for.
#
# NOT a ctest, same standing reason: every number is a timing. It exits non-zero
# only when a case fails to converge, since a leg share taken from a diverged
# solve is a share of the wrong thing.
#
# It is separate from NpcThreadScaling rather than a flag on it because the two
# answer different questions and want different statistics. That harness asks
# whether threading and PARDISO are worth taking and reports a BEST-of, which is
# right for a minimum wall time; this one asks what fraction of a step each leg
# is and reports a MEDIAN of five, which is right for a share -- and upstream's
# request says so explicitly, their own 1.36-1.59x on one loop having been
# unresolvable against a 5% end-to-end scatter.
if(TARGET meq_core AND NOT _solver_index EQUAL -1
   AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/performance/NewtonStepProfile.cpp")
	add_executable(NewtonStepProfile "performance/NewtonStepProfile.cpp")
	target_include_directories(NewtonStepProfile PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
	target_link_libraries(NewtonStepProfile PRIVATE meq_core)

	find_package(OpenMP QUIET COMPONENTS CXX)
	if(TARGET OpenMP::OpenMP_CXX)
		target_link_libraries(NewtonStepProfile PRIVATE OpenMP::OpenMP_CXX)
	endif()
endif()

# ---------------------------------------------------------------------------
# The driver, end to end. Everything above calls the library directly; this is
# the only test that goes through a config file and out to files, which is the
# path a user takes.
# ---------------------------------------------------------------------------

if(TARGET meq AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/convergence/DriverAcceptance.cpp")
	add_executable(DriverAcceptance "convergence/DriverAcceptance.cpp")
	target_include_directories(DriverAcceptance PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
	target_link_libraries(DriverAcceptance PRIVATE
		meq_core
		Boost::unit_test_framework)

	# Passed rather than guessed. A test that assumed "../meq" would keep
	# passing while running nothing at all under a different generator.
	target_compile_definitions(DriverAcceptance PRIVATE
		MEQ_DRIVER_PATH="$<TARGET_FILE:meq>"
		# And the wrapper beside it. Not a target -- file(GENERATE) writes a
		# script, so there is no $<TARGET_FILE:> for it -- but the path is
		# fixed by the same file(GENERATE) call in the top-level CMakeLists.
		MEQ_RUN_PATH="${CMAKE_BINARY_DIR}/meq-run")

	# The source tree, so examples/*.toml resolve -- and so the driver's output
	# lands somewhere predictable for the read-back.
	add_test(NAME DriverAcceptance COMMAND DriverAcceptance
		WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}")
	set_tests_properties(DriverAcceptance PROPERTIES
		ENVIRONMENT "MKL_NUM_THREADS=1")
else()
	message(STATUS "meq: DriverAcceptance not registered (needs MEQ_BUILD_APP=ON)")
endif()

# MachineFixedBoundary runs the SHIPPED examples/fixed-*.toml through the
# shipped binary, so it needs the same two things DriverAcceptance does: the
# driver's path, and the source tree as its working directory because those
# files name their profile tables as `examples/...`.
#
# IT IS A DRIVER TEST ON PURPOSE. What rots about an example is the FILE -- a
# key renamed, a default moved, a table's path -- and a library-level fixture
# would re-implement the curved path and then notice none of that.
if(TARGET meq AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/convergence/MachineFixedBoundary.cpp")
	add_executable(MachineFixedBoundary "convergence/MachineFixedBoundary.cpp")
	target_include_directories(MachineFixedBoundary PRIVATE
		"${CMAKE_CURRENT_SOURCE_DIR}")
	target_link_libraries(MachineFixedBoundary PRIVATE
		meq_core
		Boost::unit_test_framework)
	target_compile_definitions(MachineFixedBoundary PRIVATE
		MEQ_DRIVER_PATH="$<TARGET_FILE:meq>")
	add_test(NAME MachineFixedBoundary COMMAND MachineFixedBoundary
		WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}")
	set_tests_properties(MachineFixedBoundary PROPERTIES
		ENVIRONMENT "MKL_NUM_THREADS=1")
else()
	message(STATUS "meq: MachineFixedBoundary not registered (needs MEQ_BUILD_APP=ON)")
endif()

# ---------------------------------------------------------------------------
# naming -- the identifier convention, enforced
# ---------------------------------------------------------------------------

# The house convention (UpperCamelCase types and enum values, lowerCamelCase
# functions, methods, variables, parameters and members) is written down in
# ../.clang-tidy as readability-identifier-naming options. This test is what
# makes it more than a document.
#
# clang-tidy is deliberately NOT wired into CMAKE_CXX_CLANG_TIDY: that would run
# it on every compile, and it is an opt-in check rather than something to pay
# for on every build.
find_program(CLANG_TIDY_EXE NAMES clang-tidy)
mark_as_advanced(CLANG_TIDY_EXE)

if(NOT CLANG_TIDY_EXE)
	message(STATUS "meq: clang-tidy not found, naming test not registered")
elseif(NOT MEQ_CORE_SOURCES_PRESENT)
	message(STATUS "meq: no ported sources yet, naming test not registered")
else()
	# clang-tidy reads the flags for each file out of the compile database
	# that CMAKE_EXPORT_COMPILE_COMMANDS writes into the build tree, which is
	# also what limits the check to the ported code: the compile database
	# contains exactly what some target compiles. apps/meq.cpp is in no target
	# because it still speaks the old MFEM API, so it is not in the database and
	# is not checked. That is deliberate and temporary -- it violates the
	# convention wholesale, and it is renamed when it is ported, at which point
	# it joins meq_core and is checked like everything else.
	# src/meq/Estimator.{hpp,cpp} made that journey in stage 6 and is checked.
	# attic/ used to be excluded here permanently; it was removed 2026-09-07.
	#
	# HeaderFilterRegex in ../.clang-tidy is what pulls src/meq/*.hpp into the
	# check: a header is not a translation unit, so GradShafranov.hpp is seen
	# only because GradShafranov.cpp includes it.
	set(_naming_sources "")
	foreach(_source IN LISTS MEQ_CORE_SOURCES_PRESENT)
		list(APPEND _naming_sources "${CMAKE_SOURCE_DIR}/${_source}")
	endforeach()

	# --warnings-as-errors is what makes this a test rather than a report:
	# clang-tidy exits 0 while printing complaints otherwise. Only the naming
	# check is promoted, so a compiler warning from some header does not fail
	# the run -- though a compile *error* still does, which is what should
	# happen, since a file that does not parse has not been checked.
	add_test(NAME naming
		COMMAND "${CLANG_TIDY_EXE}"
			-p "${CMAKE_BINARY_DIR}"
			"--config-file=${CMAKE_SOURCE_DIR}/.clang-tidy"
			--quiet
			--warnings-as-errors=readability-identifier-naming
			${_naming_sources})
endif()
