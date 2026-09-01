# meq -- the Maryland Equilibrium Solver.
#
# A hybridizable discontinuous Galerkin solver for the fixed-boundary
# Grad-Shafranov equation, built on MFEM's DarcyForm.
#
# This is the stage-2 build of the port onto MFEM 4.9.1. Stages 1 and 2 cover
# configuration, profiles, the source term and the linear HDG Grad-Shafranov
# operator, and build them into meq_core with unit tests and a convergence test.
# The residual estimator (src/meq/Estimator.hpp) and the driver (apps/meq.cpp)
# still speak the old Waterloo HDG API, which no longer exists in 4.9.1, so they
# are not compiled by anything here yet; see MEQ_BUILD_APP below and
# attic/free-boundary/README.md for what happened to the rest.

cmake_minimum_required(VERSION 3.20)

project(meq
	VERSION 0.1.0
	DESCRIPTION "Hybridizable DG solver for the Grad-Shafranov equation"
	LANGUAGES CXX)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# For clangd and everything else that wants to know how a file is compiled.
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

option(MEQ_BUILD_TESTS "Build the unit tests" ON)

# ON since stage 7e: apps/meq.cpp is the driver now, not the legacy vacuum-field
# program built on MFEM 4.5.1's HDGBilinearForm. It needs MFEM, so the guard
# below still applies -- a build without it gets the MFEM-free half and no
# executable, which is what CI does.
option(MEQ_BUILD_APP "Build the meq executable" ON)

# gcov instrumentation, off by default. When on it is attached to meq_core
# PUBLIC, so the flags reach the test executables as well -- gcov needs both
# the compile and the link side on every object it reports, and a test binary
# that merely links an instrumented library produces no .gcda of its own.
#
# Deliberately an option rather than a build type. A coverage build wants -O0,
# and forcing that through CMAKE_BUILD_TYPE would silently drop NDEBUG for
# anyone who had set it, changing which MFEM_ASSERTs are live. Stage 4 found a
# BlockVector size mismatch that only an assert catches, so that difference is
# real and should be chosen rather than inherited.
option(MEQ_ENABLE_COVERAGE "Build with gcov instrumentation" OFF)

if(MEQ_ENABLE_COVERAGE)
	if(NOT CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
		message(FATAL_ERROR "MEQ_ENABLE_COVERAGE needs gcc or clang")
	endif()
	add_library(meq_coverage INTERFACE)
	target_compile_options(meq_coverage INTERFACE --coverage -O0 -g)
	target_link_options(meq_coverage INTERFACE --coverage)
	message(STATUS "meq: coverage instrumentation ON (-O0, --coverage)")
endif()

# A single-configuration generator with no build type gets no optimisation at
# all, which for a finite element code is the difference between a test suite
# that runs and one that is waited on.
get_property(_meq_multi_config GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)
if(NOT _meq_multi_config AND NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
	set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
		Debug Release RelWithDebInfo MinSizeRel)
endif()

# ---------------------------------------------------------------------------
# Dependencies
# ---------------------------------------------------------------------------

# MODULE mode explicitly: cmake/FindMFEM.cmake reads the flags out of the make
# build's config/config.mk, which is the only description an in-source MFEM
# build leaves behind. See that file for why there is no config package to use.
# NOT REQUIRED, deliberately.
#
# meq's solver needs MFEM, but half the library does not: Config, Profiles,
# Source and SourceFactory are free of it by design, which is what lets the
# configuration and profile layers be unit tested without the numerics.
#
# That split is what makes continuous integration possible at all here. The
# MFEM meq builds against is a local branch (gf-hdg-subdomains-dev in
# ../mfem-hdg-dev) that is not published on any remote, so a hosted CI runner
# cannot obtain it. Without this, CMake could not configure there and nothing
# could be tested. With it, CI builds and tests the MFEM-free half and runs
# the naming check over it.
#
# A developer build finds MFEM and gets everything; nothing about the full
# build changes. If MFEM is missing the configure output says exactly what is
# being skipped rather than quietly building less than you asked for.
find_package(MFEM MODULE)

if(NOT MFEM_FOUND)
	message(STATUS "meq: MFEM not found -- building the MFEM-free subset only.")
	message(STATUS "meq:   excluded: the solver, the estimator, every convergence test")
	message(STATUS "meq:   set MFEM_DIR to a built MFEM to get the whole library")
endif()

# NetCDF, for the driver's output grid. Optional in the same spirit as MFEM: a
# tree without it still builds and tests everything else, and the writer is the
# only thing that goes missing. Found through pkg-config because that is how
# netcdf-cxx4 ships -- it has no CMake config -- and the C++4 API is the one
# MaNTA writes with, which matters because the two codes are meant to exchange
# these files.
find_package(PkgConfig QUIET)
if(PkgConfig_FOUND)
	pkg_check_modules(NETCDFCXX QUIET IMPORTED_TARGET netcdf-cxx4)
endif()

if(NETCDFCXX_FOUND)
	message(STATUS "meq: netcdf-cxx4 ${NETCDFCXX_VERSION}")
else()
	message(STATUS "meq: netcdf-cxx4 not found -- the NetCDF writer is excluded.")
	message(STATUS "meq:   install libnetcdf-c++4-dev to get it")
endif()

# toml11, header only, used by Config. The checkout carries two layouts:
#
#   include/toml.hpp        + include/toml11/*.hpp   -- the normal one
#   single_include/toml.hpp                          -- one amalgamated header
#
# include/ is the one to put on the search path. It resolves both #include
# "toml.hpp" (which is what src/meq/Config.hpp uses) and the individual
# #include "toml11/value.hpp" style, whereas single_include/ resolves only the
# first. Putting both on the path would give two different toml.hpp for the
# preprocessor to choose between, so only one is added.
# toml11 is a submodule at extern/toml11, pinned to a commit, as MaNTA does it.
# MFEM deliberately is not: its history is enormous and it needs an out-of-tree
# build anyway, so it stays a hand-built library found through MFEM_DIR.
if(NOT DEFINED TOML11_DIR AND DEFINED ENV{TOML11_DIR})
	set(_toml11_dir_default "$ENV{TOML11_DIR}")
else()
	set(_toml11_dir_default "${CMAKE_CURRENT_SOURCE_DIR}/extern/toml11")
endif()

set(TOML11_DIR "${_toml11_dir_default}" CACHE PATH
	"Root of the toml11 checkout (header only)")
unset(_toml11_dir_default)

find_path(TOML11_INCLUDE_DIR
	NAMES toml.hpp
	HINTS "${TOML11_DIR}"
	PATH_SUFFIXES include single_include
	NO_DEFAULT_PATH
	DOC "Directory containing toml11's toml.hpp")

if(NOT TOML11_INCLUDE_DIR)
	if(TOML11_DIR STREQUAL "${CMAKE_CURRENT_SOURCE_DIR}/extern/toml11")
		# Much the most likely cause: a fresh clone without --recursive.
		message(FATAL_ERROR
			"toml11 not found at ${TOML11_DIR}.\n"
			"It is a git submodule and has not been checked out. Run:\n"
			"    git submodule update --init --recursive\n"
			"Or point -DTOML11_DIR=... at a toml11 checkout elsewhere.")
	else()
		message(FATAL_ERROR
			"toml11 not found under TOML11_DIR=${TOML11_DIR}.\n"
			"Point -DTOML11_DIR=... at a toml11 checkout containing "
			"include/toml.hpp or single_include/toml.hpp.")
	endif()
endif()

mark_as_advanced(TOML11_INCLUDE_DIR)

# An IMPORTED INTERFACE target so the include directory is treated as a system
# one: toml11 is a deep template library and its warnings are not ours to fix.
add_library(toml11::toml11 INTERFACE IMPORTED)
set_target_properties(toml11::toml11 PROPERTIES
	INTERFACE_INCLUDE_DIRECTORIES "${TOML11_INCLUDE_DIR}")

# ---------------------------------------------------------------------------
# meq_core
# ---------------------------------------------------------------------------

# Listed explicitly rather than globbed: a glob is not re-run when a file
# appears, so a new source silently fails to build until someone re-runs cmake
# by hand, and a stale file left in the tree silently starts building again.
#
# src/meq/GradShafranov.cpp joined this list in stage 2, ported from
# HDGBilinearForm to DarcyForm. src/meq/Estimator.{hpp,cpp} joined it in stage 6:
# the residual estimator of refs/HDG-GradShafranov-Adaptive.pdf eq (20), rewritten
# rather than ported. The old header called GridFunction::GetValueFacet, which was
# a patch carried by the old MFEM branch and is not in 4.9.1, and it also used the
# raw psi_h where the published estimator uses the post-processed psi*.
#
# src/meq/Solution.hpp is gone rather than ported: it wrapped hand-rolled block
# offsets and a second set of finite element spaces around what DarcyForm now
# owns itself. What was left of its job lives on GradShafranovSolver.
set(MEQ_CORE_SOURCES
	src/meq/BoundaryShape.cpp
	src/meq/Config.cpp
	src/meq/Estimator.cpp
	src/meq/Field.cpp
	src/meq/GradShafranov.cpp
	src/meq/Output.cpp
	src/meq/Profiles.cpp
	src/meq/RotatingSource.cpp
	src/meq/Sampler.cpp
	src/meq/Source.cpp
	src/meq/SourceFactory.cpp
	src/meq/WarmStart.cpp)

# Sources that will not compile without MFEM. Filtered out when it is absent;
# see the find_package note above for why that case exists at all.
set(MEQ_CORE_SOURCES_NEEDING_MFEM
	src/meq/Estimator.cpp
	src/meq/Field.cpp
	src/meq/GradShafranov.cpp
	src/meq/Output.cpp
	src/meq/Sampler.cpp
	src/meq/WarmStart.cpp)

# And the one that additionally needs netcdf-cxx4. Kept separate because the two
# dependencies are independent: a tree with MFEM and no NetCDF builds everything
# but the grid writer.
set(MEQ_CORE_SOURCES_NEEDING_NETCDF
	src/meq/Output.cpp)

# The port is being done file by file, so a source in the list above may not
# have been written yet. Configure with what is there and say what is missing,
# rather than failing and blocking everything else in the tree.
set(MEQ_CORE_SOURCES_PRESENT "")
set(_meq_core_missing "")
foreach(_source IN LISTS MEQ_CORE_SOURCES)
	list(FIND MEQ_CORE_SOURCES_NEEDING_MFEM "${_source}" _needs_mfem)
	list(FIND MEQ_CORE_SOURCES_NEEDING_NETCDF "${_source}" _needs_netcdf)
	if(NOT MFEM_FOUND AND NOT _needs_mfem EQUAL -1)
		# Not "missing" -- deliberately excluded. Kept out of _meq_core_missing
		# so the two reasons a source is absent stay tellable apart.
		list(APPEND _meq_core_no_mfem "${_source}")
	elseif(NOT NETCDFCXX_FOUND AND NOT _needs_netcdf EQUAL -1)
		list(APPEND _meq_core_no_netcdf "${_source}")
	elseif(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${_source}")
		list(APPEND MEQ_CORE_SOURCES_PRESENT "${_source}")
	else()
		list(APPEND _meq_core_missing "${_source}")
	endif()
endforeach()

if(_meq_core_missing)
	string(REPLACE ";" ", " _meq_core_missing_text "${_meq_core_missing}")
	message(STATUS "meq: not yet present, omitted from meq_core: ${_meq_core_missing_text}")
endif()

if(_meq_core_no_netcdf)
	string(REPLACE ";" ", " _meq_core_no_netcdf_text "${_meq_core_no_netcdf}")
	message(STATUS "meq: no netcdf-cxx4, omitted from meq_core: ${_meq_core_no_netcdf_text}")
endif()

if(MEQ_CORE_SOURCES_PRESENT)
	add_library(meq_core STATIC ${MEQ_CORE_SOURCES_PRESENT})

	# Both src/ and src/meq/ are on the include path, so a header can be
	# reached either as "Config.hpp" (which is how the existing sources and
	# apps/meq.cpp include each other) or as "meq/Config.hpp".
	target_include_directories(meq_core PUBLIC
		"${CMAKE_CURRENT_SOURCE_DIR}/src"
		"${CMAKE_CURRENT_SOURCE_DIR}/src/meq")

	if(MEQ_ENABLE_COVERAGE)
		# PUBLIC so the flags reach the test executables as well. gcov needs
		# both the compile and the link side on every object it reports, and a
		# test binary that links an instrumented library but is not itself
		# instrumented produces no .gcda for its own translation unit.
		target_link_libraries(meq_core PUBLIC meq_coverage)
	endif()

	# PUBLIC, not PRIVATE: meq's own headers include mfem.hpp and toml.hpp,
	# so anything that includes them needs both on its compile line.
	target_link_libraries(meq_core PUBLIC toml11::toml11)

	# Guarded, because an unfound package leaves no target and CMake passes the
	# name through to the linker verbatim -- the failure is
	# "cannot find -lMFEM::MFEM" from ld, at link time, long after configure
	# said everything was fine.
	if(MFEM_FOUND)
		target_link_libraries(meq_core PUBLIC MFEM::MFEM)
	endif()

	# PRIVATE, unlike MFEM: netcdf appears only inside Output.cpp and in no
	# header, so nothing that includes meq's headers needs it on its line.
	if(NETCDFCXX_FOUND)
		target_link_libraries(meq_core PRIVATE PkgConfig::NETCDFCXX)
		target_compile_definitions(meq_core PUBLIC MEQ_USE_NETCDF)
	endif()

	# Guarded by compiler id rather than assumed: these are the GCC/Clang
	# spellings, and nothing here is worth breaking a configure over on a
	# compiler that spells them differently.
	target_compile_options(meq_core PRIVATE
		$<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-Wall;-Wextra>)
else()
	message(WARNING
		"meq: none of the meq_core sources exist yet, so the library target "
		"was not created. Nothing depending on it will be built.")
endif()

# ---------------------------------------------------------------------------
# Targets
# ---------------------------------------------------------------------------

if(MEQ_BUILD_APP)
	if(NOT TARGET meq_core)
		message(FATAL_ERROR "MEQ_BUILD_APP is ON but meq_core has no sources to build from.")
	endif()

	# The commit, so a run can be traced back to the tree that produced it. Not
	# fatal when git is absent or this is an unpacked tarball -- the version
	# falls back to the project version alone.
	set(MEQ_GIT_COMMIT "unknown")
	find_package(Git QUIET)
	if(GIT_FOUND AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/.git")
		execute_process(
			COMMAND "${GIT_EXECUTABLE}" rev-parse --short HEAD
			WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
			OUTPUT_VARIABLE MEQ_GIT_COMMIT
			OUTPUT_STRIP_TRAILING_WHITESPACE
			ERROR_QUIET)
	endif()

	add_executable(meq apps/meq.cpp)
	target_link_libraries(meq PRIVATE meq_core)
	configure_file(
		"${CMAKE_CURRENT_SOURCE_DIR}/cmake/Version.hpp.in"
		"${CMAKE_CURRENT_BINARY_DIR}/meq_version.hpp"
		@ONLY)
	target_include_directories(meq PRIVATE "${CMAKE_CURRENT_BINARY_DIR}")
endif()

if(MEQ_BUILD_TESTS)
	enable_testing()
	add_subdirectory(tests)
endif()

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------

message(STATUS "meq ${PROJECT_VERSION}: ${CMAKE_BUILD_TYPE} build with ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "  MFEM      ${MFEM_VERSION} at ${MFEM_DIR}")
if(MFEM_FEATURES)
	string(REPLACE ";" " " _meq_mfem_features "${MFEM_FEATURES}")
	string(REPLACE "MFEM_USE_" "" _meq_mfem_features "${_meq_mfem_features}")
	message(STATUS "            options: ${_meq_mfem_features}")
endif()
message(STATUS "  toml11    ${TOML11_INCLUDE_DIR}")
message(STATUS "  app       ${MEQ_BUILD_APP} (MEQ_BUILD_APP)")
message(STATUS "  tests     ${MEQ_BUILD_TESTS} (MEQ_BUILD_TESTS)")
