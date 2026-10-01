# Regenerate meq_version.hpp from the commit the tree is on RIGHT NOW.
#
# Run in script mode (`cmake -P`) as a build step, not at configure time.
#
# WHY NOT execute_process() IN CMakeLists.txt, WHICH IS THE OBVIOUS PLACE.
# It runs once, when the build directory is configured, so the commit it reads
# is frozen there -- and a long-lived build directory then stamps every binary
# it produces, and every output file that binary writes, with whatever HEAD
# happened to be the day somebody first ran `cmake -B build`. Measured: a `meq`
# linked from 6a078e2 reported `0.1.0 (6d4b68a)` and the netCDF it wrote carried
# that same wrong string as its `meq_version` attribute -- five commits stale,
# with nothing anywhere saying so.
#
# THAT IS A PROVENANCE DEFECT AND NOT A COSMETIC ONE. `CLAUDE.md`'s own method
# for telling *this was already failing* from *I broke it* is to build HEAD in a
# `git worktree` and compare two binaries; the device programme is a comparison
# between two builds throughout; and `docs/running.rst` tells a user that the
# answer to "which equilibrium did this plot come from" is a file that can be
# re-run. All three rest on a binary knowing which tree it came from.
#
# copy_if_different IS WHAT KEEPS THIS CHEAP. The header is rewritten only when
# the commit actually moves, so `meq` relinks when the answer changes and not
# once per build.
#
# Inputs, all -D on the command line: MEQ_SOURCE_DIR, MEQ_BINARY_DIR,
# MEQ_PROJECT_VERSION.

# The template spells the version `@PROJECT_VERSION@`, which is a project()
# variable and does not exist in script mode. Restated from the -D so that
# cmake/Version.hpp.in stays a single template both callers can use.
set(PROJECT_VERSION "${MEQ_PROJECT_VERSION}")

set(MEQ_GIT_COMMIT "unknown")

find_package(Git QUIET)
if(GIT_FOUND AND EXISTS "${MEQ_SOURCE_DIR}/.git")
	execute_process(
		COMMAND "${GIT_EXECUTABLE}" rev-parse --short HEAD
		WORKING_DIRECTORY "${MEQ_SOURCE_DIR}"
		OUTPUT_VARIABLE MEQ_GIT_COMMIT
		OUTPUT_STRIP_TRAILING_WHITESPACE
		ERROR_QUIET)
	if(MEQ_GIT_COMMIT STREQUAL "")
		set(MEQ_GIT_COMMIT "unknown")
	endif()

	# A BINARY BUILT FROM AN EDITED TREE IS NOT THE COMMIT IT NAMES, and saying
	# so is the difference between a provenance string and a guess. Every
	# measurement this project publishes is taken from a working tree at some
	# point, and `CLAUDE.md` records a session committing a subagent's
	# scaffolding because a build was taken as evidence of content.
	#
	# Tracked files only: untracked output lands in the repository root by the
	# hundred -- `.mesh`, `.gf`, `.nc` from every example anyone has run -- so
	# including it would mark every binary dirty and the flag would mean
	# nothing.
	execute_process(
		COMMAND "${GIT_EXECUTABLE}" status --porcelain --untracked-files=no
		WORKING_DIRECTORY "${MEQ_SOURCE_DIR}"
		OUTPUT_VARIABLE _meq_git_dirty
		OUTPUT_STRIP_TRAILING_WHITESPACE
		ERROR_QUIET)
	if(NOT _meq_git_dirty STREQUAL "")
		set(MEQ_GIT_COMMIT "${MEQ_GIT_COMMIT}-dirty")
	endif()
endif()

configure_file(
	"${MEQ_SOURCE_DIR}/cmake/Version.hpp.in"
	"${MEQ_BINARY_DIR}/meq_version.hpp.candidate"
	@ONLY)

execute_process(
	COMMAND "${CMAKE_COMMAND}" -E copy_if_different
		"${MEQ_BINARY_DIR}/meq_version.hpp.candidate"
		"${MEQ_BINARY_DIR}/meq_version.hpp")
