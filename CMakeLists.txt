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

# Still off after stage 2. apps/meq.cpp constructs an HDGBilinearForm and drives
# the old estimator, neither of which exists in MFEM 4.9.1, so switching this on
# today produces a wall of compiler errors rather than a program. The solver core
# is ported; the driver is a separate step.
option(MEQ_BUILD_APP "Build the meq executable (needs the stage-2 solver port)" OFF)

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
find_package(MFEM MODULE REQUIRED)

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
# HDGBilinearForm to DarcyForm. src/meq/Estimator.hpp has not -- it calls
# GridFunction::GetValueFacet, which was a patch carried by the old MFEM branch
# and is not in 4.9.1 -- and is stage 6 work in any case.
#
# src/meq/Solution.hpp is gone rather than ported: it wrapped hand-rolled block
# offsets and a second set of finite element spaces around what DarcyForm now
# owns itself. What was left of its job lives on GradShafranovSolver.
set(MEQ_CORE_SOURCES
	src/meq/Config.cpp
	src/meq/GradShafranov.cpp
	src/meq/Profiles.cpp
	src/meq/Source.cpp)

# The port is being done file by file, so a source in the list above may not
# have been written yet. Configure with what is there and say what is missing,
# rather than failing and blocking everything else in the tree.
set(MEQ_CORE_SOURCES_PRESENT "")
set(_meq_core_missing "")
foreach(_source IN LISTS MEQ_CORE_SOURCES)
	if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${_source}")
		list(APPEND MEQ_CORE_SOURCES_PRESENT "${_source}")
	else()
		list(APPEND _meq_core_missing "${_source}")
	endif()
endforeach()

if(_meq_core_missing)
	string(REPLACE ";" ", " _meq_core_missing_text "${_meq_core_missing}")
	message(STATUS "meq: not yet present, omitted from meq_core: ${_meq_core_missing_text}")
endif()

if(MEQ_CORE_SOURCES_PRESENT)
	add_library(meq_core STATIC ${MEQ_CORE_SOURCES_PRESENT})

	# Both src/ and src/meq/ are on the include path, so a header can be
	# reached either as "Config.hpp" (which is how the existing sources and
	# apps/meq.cpp include each other) or as "meq/Config.hpp".
	target_include_directories(meq_core PUBLIC
		"${CMAKE_CURRENT_SOURCE_DIR}/src"
		"${CMAKE_CURRENT_SOURCE_DIR}/src/meq")

	# PUBLIC, not PRIVATE: meq's own headers include mfem.hpp and toml.hpp,
	# so anything that includes them needs both on its compile line.
	target_link_libraries(meq_core PUBLIC
		MFEM::MFEM
		toml11::toml11)

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

	add_executable(meq apps/meq.cpp)
	target_link_libraries(meq PRIVATE meq_core)
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
message(STATUS "  app       ${MEQ_BUILD_APP} (MEQ_BUILD_APP; still unported)")
message(STATUS "  tests     ${MEQ_BUILD_TESTS} (MEQ_BUILD_TESTS)")
