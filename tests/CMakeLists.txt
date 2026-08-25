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
	ConfigTests
	ProfilesTests
	SourceTests)

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

		# WORKING_DIRECTORY is load-bearing, not tidiness: ConfigTests opens
		# examples/*.toml by a path relative to the source tree, so without
		# this it passes from build/ and fails from anywhere else. Pinning it
		# makes an out-of-tree build behave the same as an in-tree one.
		add_test(NAME ${_test} COMMAND ${_test}
			WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}")

		# MKL_THREADING_LAYER=GNU is not optional on this machine.
		# /usr/lib/x86_64-linux-gnu/libblas.so.3 is an alternatives
		# symlink pointing at libmkl_rt.so, and MKL's default threading
		# layer silently returns wrong BLAS-3 results to UMFPACK -- which
		# MFEM's direct solver uses. Wrong answers, no error, no warning.
		# The same trap is recorded in the MFEM tree's CLAUDE.md, where
		# every run of its own test suite needs the same variable.
		set_tests_properties(${_test} PROPERTIES
			ENVIRONMENT "MKL_THREADING_LAYER=GNU")
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
	SolovievConvergence)

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

		# For "analytic/Soloviev.hpp". meq_core already carries src/ and
		# src/meq/, so "meq/GradShafranov.hpp" needs nothing here.
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
			ENVIRONMENT "MKL_THREADING_LAYER=GNU")
	endif()
endforeach()

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
	# contains exactly what some target compiles. src/meq/Estimator.hpp and
	# apps/meq.cpp are in no target because they still speak the old MFEM API,
	# so they are not in the database and are not checked. That is deliberate
	# and temporary -- they violate the convention wholesale, and they are
	# renamed when they are ported, at which point they join meq_core and are
	# checked like everything else. attic/ is excluded permanently; see its
	# README.
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
