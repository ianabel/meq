# FindMFEM.cmake -- locate MFEM and describe it as the imported target MFEM::MFEM.
#
# meq builds against an *installed* MFEM: ../mfem/install, produced by a CMake
# build in ../mfem/build from sources in ../mfem/mfem-src, which sits on the
# gf-hdg-subdomains-dev branch. That library is built with everything meq wants
# -- SUNDIALS for KINSol, GSLIB for FindPoints, SuiteSparse for UMFPACK, LAPACK
# -- and is meant to stay put while ../mfem-hdg-dev continues to be rebuilt for
# development. Pointing meq at a tree somebody is actively editing is how you
# get a suite that fails for reasons that have nothing to do with meq.
#
# It still works against an in-source make build, and that is deliberate: a
# hand-built ../mfem-hdg-dev has libmfem.a and mfem.hpp together at the top of
# the tree, and this module handles both layouts. Set MFEM_DIR to whichever is
# wanted.
#
# What both layouts leave behind is config.mk -- in config/ for a source tree
# and in share/mfem/ for an install -- which records every option the library
# was compiled with and every third-party library it must be linked against.
# This module parses that file rather than hardcoding the flags, because the
# same MFEM gets rebuilt with different options and each change moves the link
# line. Anything absent from config.mk is simply left empty; nothing here
# requires a particular variable to be there. (A CMake install also ships
# MFEMConfig.cmake, so find_package(MFEM CONFIG) would work -- but only for the
# install, and supporting one route for both layouts is worth more than using
# the more idiomatic one for half of them.)
#
# Input variables:
#
#   MFEM_DIR       Root of the MFEM tree. Cache variable, so -DMFEM_DIR=... on
#                  the command line wins; otherwise the environment variable
#                  MFEM_DIR is used if set; otherwise ../mfem/install next to
#                  this project.
#
# Result variables:
#
#   MFEM_FOUND             true if both the header and the library were found
#   MFEM_VERSION           e.g. 4.9.1
#   MFEM_INCLUDE_DIRS      MFEM's own include dir plus the third-party ones
#   MFEM_LIBRARIES         libmfem plus the third-party link line
#   MFEM_USE_<FEATURE>     YES/NO for each build option in config.mk, so that a
#                          caller can react to MFEM_USE_MPI, MFEM_USE_SUNDIALS
#                          and friends without parsing anything itself
#
# Imported target:
#
#   MFEM::MFEM     carries the include directories, the link line, any compile
#                  options MFEM's third-party flags demand, and the C++
#                  standard MFEM itself was compiled with.

include(FindPackageHandleStandardArgs)

# Command line beats environment beats the sibling-checkout default. Note that
# after the first configure MFEM_DIR is in the cache, so the environment is not
# consulted again -- change it with -DMFEM_DIR=... or by editing the cache.
if(NOT DEFINED MFEM_DIR AND DEFINED ENV{MFEM_DIR})
	set(_mfem_dir_default "$ENV{MFEM_DIR}")
else()
	set(_mfem_dir_default "${CMAKE_SOURCE_DIR}/../mfem/install")
endif()

set(MFEM_DIR "${_mfem_dir_default}" CACHE PATH
	"Root of the MFEM tree to build against (source tree or install prefix)")
unset(_mfem_dir_default)

# Both layouts: an in-source make build has mfem.hpp and libmfem.a at the root,
# an installed one has include/ and lib/.
find_path(MFEM_INCLUDE_DIR
	NAMES mfem.hpp
	HINTS "${MFEM_DIR}"
	PATH_SUFFIXES include
	DOC "Directory containing mfem.hpp")

find_library(MFEM_LIBRARY
	NAMES mfem
	HINTS "${MFEM_DIR}"
	PATH_SUFFIXES lib lib64
	DOC "The MFEM library")

# config.mk lives in config/ in a source tree and in share/mfem/ in an install.
# Deliberately NO_DEFAULT_PATH: a config.mk picked up from somewhere else on the
# system would describe a different library than the one we just found.
find_file(MFEM_CONFIG_MK
	NAMES config.mk
	HINTS "${MFEM_DIR}"
	PATH_SUFFIXES config share/mfem
	NO_DEFAULT_PATH
	DOC "MFEM's generated config.mk, which records its build options")

# ---------------------------------------------------------------------------
# Parse config.mk
# ---------------------------------------------------------------------------

set(_mfem_config_text "")
if(MFEM_CONFIG_MK)
	file(READ "${MFEM_CONFIG_MK}" _mfem_config_text)
endif()

# Read one `NAME = value` assignment out of config.mk. A macro rather than a
# function so the result lands in the caller's scope. The name must be followed
# by whitespace and then '=', which is what stops MFEM_USE_METIS from matching
# the MFEM_USE_METIS_5 line and MFEM_CXX from matching MFEM_CXXFLAGS. Absent
# variables yield an empty string.
macro(_mfem_config_value _out _name)
	set(${_out} "")
	if(NOT _mfem_config_text STREQUAL "")
		string(REGEX MATCH "(^|\n)[ \t]*${_name}[ \t]*=[^\n]*" _mfem_line "${_mfem_config_text}")
		if(NOT _mfem_line STREQUAL "")
			string(REGEX REPLACE "^[ \t\n]*${_name}[ \t]*=[ \t]*" "" ${_out} "${_mfem_line}")
			string(STRIP "${${_out}}" ${_out})
		endif()
		unset(_mfem_line)
	endif()
endmacro()

_mfem_config_value(MFEM_VERSION_STRING MFEM_VERSION_STRING)
_mfem_config_value(_mfem_cxxflags       MFEM_CXXFLAGS)
_mfem_config_value(_mfem_cppflags       MFEM_CPPFLAGS)
_mfem_config_value(_mfem_tplflags       MFEM_TPLFLAGS)
_mfem_config_value(_mfem_ext_libs       MFEM_EXT_LIBS)
_mfem_config_value(MFEM_GIT_STRING      MFEM_GIT_STRING)

# Every MFEM_USE_* switch, verbatim, so callers can test for a feature without
# reading config.mk themselves. Read out of the file rather than from a fixed
# list, so options added by a future MFEM release appear without editing this.
string(REGEX MATCHALL "(^|\n)[ \t]*MFEM_USE_[A-Z0-9_]+[ \t]*=[^\n]*" _mfem_use_lines "${_mfem_config_text}")
set(MFEM_FEATURES "")
foreach(_line IN LISTS _mfem_use_lines)
	string(REGEX REPLACE "^[ \t\n]*(MFEM_USE_[A-Z0-9_]+)[ \t]*=.*$" "\\1" _name "${_line}")
	string(REGEX REPLACE "^[^=]*=[ \t]*" "" _value "${_line}")
	string(STRIP "${_value}" _value)
	set(${_name} "${_value}")
	if(_value STREQUAL "YES")
		list(APPEND MFEM_FEATURES "${_name}")
	endif()
endforeach()
unset(_mfem_use_lines)

# The version. config.mk is authoritative; config.hpp is the fallback for an
# install that ships no config.mk.
if(MFEM_VERSION_STRING)
	set(MFEM_VERSION "${MFEM_VERSION_STRING}")
elseif(MFEM_INCLUDE_DIR AND EXISTS "${MFEM_INCLUDE_DIR}/config/_config.hpp")
	file(STRINGS "${MFEM_INCLUDE_DIR}/config/_config.hpp" _mfem_version_line
		REGEX "^#define[ \t]+MFEM_VERSION_STRING[ \t]+")
	string(REGEX REPLACE ".*\"([^\"]+)\".*" "\\1" MFEM_VERSION "${_mfem_version_line}")
	unset(_mfem_version_line)
else()
	set(MFEM_VERSION "")
endif()

# ---------------------------------------------------------------------------
# Turn the parsed make flags into CMake usage requirements
# ---------------------------------------------------------------------------

set(MFEM_INCLUDE_DIRS "")
set(MFEM_COMPILE_DEFINITIONS "")
set(MFEM_COMPILE_OPTIONS "")

if(MFEM_INCLUDE_DIR)
	list(APPEND MFEM_INCLUDE_DIRS "${MFEM_INCLUDE_DIR}")
endif()

# MFEM_TPLFLAGS and MFEM_CPPFLAGS carry the third-party include directories --
# on this machine, -I/usr/include/suitesparse, which is exactly what a program
# using DarcyForm needs in order to see the UMFPACK headers MFEM's own headers
# include. Sort them into include dirs, defines and everything else.
separate_arguments(_mfem_tpl_list UNIX_COMMAND "${_mfem_tplflags} ${_mfem_cppflags}")
foreach(_flag IN LISTS _mfem_tpl_list)
	if(_flag MATCHES "^-I(.+)$")
		list(APPEND MFEM_INCLUDE_DIRS "${CMAKE_MATCH_1}")
	elseif(_flag MATCHES "^-D(.+)$")
		list(APPEND MFEM_COMPILE_DEFINITIONS "${CMAKE_MATCH_1}")
	elseif(NOT _flag STREQUAL "")
		list(APPEND MFEM_COMPILE_OPTIONS "${_flag}")
	endif()
endforeach()
unset(_mfem_tpl_list)
list(REMOVE_DUPLICATES MFEM_INCLUDE_DIRS)

# MFEM_CXXFLAGS is mostly optimisation settings, which are the consumer's
# business and not MFEM's -- CMAKE_BUILD_TYPE decides those. Two things in it do
# have to be propagated, because they change the ABI of MFEM's headers rather
# than just the code generation: the C++ standard, and OpenMP.
set(MFEM_CXX_STANDARD "")
if(_mfem_cxxflags MATCHES "-std=[a-z+]+([0-9][0-9])")
	set(MFEM_CXX_STANDARD "${CMAKE_MATCH_1}")
endif()

if(_mfem_cxxflags MATCHES "(^| )-fopenmp( |$)")
	list(APPEND MFEM_COMPILE_OPTIONS "-fopenmp")
endif()

# The link line. Items beginning with '-' (the -L and -l of MFEM_EXT_LIBS) are
# passed to the linker verbatim and in order, which is what a static libmfem.a
# needs: it comes first, its dependencies after.
set(MFEM_LIBRARIES "")
if(MFEM_LIBRARY)
	list(APPEND MFEM_LIBRARIES "${MFEM_LIBRARY}")
endif()
separate_arguments(_mfem_ext_list UNIX_COMMAND "${_mfem_ext_libs}")
foreach(_flag IN LISTS _mfem_ext_list)
	if(NOT _flag STREQUAL "")
		list(APPEND MFEM_LIBRARIES "${_flag}")
	endif()
endforeach()
unset(_mfem_ext_list)

# A parallel MFEM needs MPI on the consumer's compile line too, since mfem.hpp
# includes mpi.h. meq itself is serial for now; this is here so that pointing
# MFEM_DIR at a parallel build fails at find_package(MPI) with something
# intelligible rather than at the first #include.
if(MFEM_USE_MPI STREQUAL "YES")
	find_package(MPI REQUIRED COMPONENTS CXX)
	list(APPEND MFEM_LIBRARIES MPI::MPI_CXX)
endif()

find_package_handle_standard_args(MFEM
	REQUIRED_VARS MFEM_LIBRARY MFEM_INCLUDE_DIR
	VERSION_VAR MFEM_VERSION)

if(MFEM_FOUND AND NOT TARGET MFEM::MFEM)
	# UNKNOWN rather than STATIC: the same tree can be built either way
	# (config.mk's MFEM_SHARED says which), and UNKNOWN handles both.
	add_library(MFEM::MFEM UNKNOWN IMPORTED)
	set_target_properties(MFEM::MFEM PROPERTIES
		IMPORTED_LOCATION "${MFEM_LIBRARY}"
		INTERFACE_INCLUDE_DIRECTORIES "${MFEM_INCLUDE_DIRS}")

	# The library itself is IMPORTED_LOCATION; INTERFACE_LINK_LIBRARIES is
	# everything it in turn depends on.
	list(REMOVE_AT MFEM_LIBRARIES 0)
	if(MFEM_LIBRARIES)
		set_target_properties(MFEM::MFEM PROPERTIES
			INTERFACE_LINK_LIBRARIES "${MFEM_LIBRARIES}")
	endif()
	list(INSERT MFEM_LIBRARIES 0 "${MFEM_LIBRARY}")

	if(MFEM_COMPILE_DEFINITIONS)
		set_target_properties(MFEM::MFEM PROPERTIES
			INTERFACE_COMPILE_DEFINITIONS "${MFEM_COMPILE_DEFINITIONS}")
	endif()

	if(MFEM_COMPILE_OPTIONS)
		set_target_properties(MFEM::MFEM PROPERTIES
			INTERFACE_COMPILE_OPTIONS "${MFEM_COMPILE_OPTIONS}")
	endif()

	# Whatever MFEM was compiled with is the floor for anything that includes
	# its headers, so it travels with the target.
	if(MFEM_CXX_STANDARD)
		set_target_properties(MFEM::MFEM PROPERTIES
			INTERFACE_COMPILE_FEATURES "cxx_std_${MFEM_CXX_STANDARD}")
	endif()
endif()

mark_as_advanced(MFEM_INCLUDE_DIR MFEM_LIBRARY MFEM_CONFIG_MK)
