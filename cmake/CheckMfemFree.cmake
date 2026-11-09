# The MfemFreeIncludes ctest: no file of the MFEM-free half reaches an MFEM
# header, directly or through any header it includes.
#
# CI builds that half and nothing else (.github/workflows/ci.yml), because the
# MFEM branch meq needs is published nowhere. A developer's tree always has
# MFEM, so an MFEM-free file that picks up mfem.hpp through some meq header
# still compiles here, and the first anyone hears of it is a CI run that no
# longer compiles. This check makes that a red test in the normal suite.
#
# Each listed file is PREPROCESSED with exactly the command the build compiles
# it with, read out of compile_commands.json, plus -H, which prints every
# header opened and how deep. MFEM's include directories are on that command
# line in a tree with MFEM, so an MFEM header IS reached, and -H names it and
# the chain of meq headers that led there. In a tree without MFEM the
# preprocessor fails on the missing header instead, which is the same verdict.
#
# Run as cmake -P with
#   COMPILE_COMMANDS  the build tree's compile_commands.json
#   LIST              a file of absolute source paths, one per line

foreach(_var COMPILE_COMMANDS LIST)
	if(NOT DEFINED ${_var})
		message(FATAL_ERROR "CheckMfemFree.cmake needs -D${_var}=")
	endif()
endforeach()

file(READ "${COMPILE_COMMANDS}" _database)
file(STRINGS "${LIST}" _sources)

# An MFEM header is mfem.hpp, or a header in one of the subdirectories of an
# installed MFEM's include/mfem/ -- general/, linalg/, config/ and the rest.
# By name rather than by MFEM's include directories: FindMFEM's
# MFEM_INCLUDE_DIRS carries MFEM's third-party directories too, and one of
# those may be /usr/include, where every standard header lives.
function(_is_mfem_header _path _out)
	if(_path MATCHES "/mfem\\.hpp$" OR _path MATCHES "/mfem/[a-z_]+/[^/]+\\.h(pp)?$")
		set(${_out} TRUE PARENT_SCOPE)
	else()
		set(${_out} FALSE PARENT_SCOPE)
	endif()
endfunction()

# One pass over the database, keeping the entries for the listed files: each
# string(JSON) on the whole document parses all of it, so a lookup per file
# would cost files x entries parses.
string(JSON _entries LENGTH "${_database}")
math(EXPR _last "${_entries} - 1")
foreach(_i RANGE ${_last})
	string(JSON _entry GET "${_database}" ${_i})
	string(JSON _file GET "${_entry}" file)
	list(FIND _sources "${_file}" _at)
	if(NOT _at EQUAL -1)
		string(JSON _command_${_at} GET "${_entry}" command)
		string(JSON _directory_${_at} GET "${_entry}" directory)
	endif()
endforeach()

set(_failures "")
set(_checked 0)
set(_at 0)
foreach(_source IN LISTS _sources)
	set(_command "${_command_${_at}}")
	set(_directory "${_directory_${_at}}")
	math(EXPR _at "${_at} + 1")
	if(NOT _command)
		list(APPEND _failures "${_source}: not in ${COMPILE_COMMANDS}, so not checked")
		continue()
	endif()

	# The compile command, less what writes an output, then preprocess only.
	separate_arguments(_args UNIX_COMMAND "${_command}")
	set(_preprocess "")
	set(_skip_next FALSE)
	foreach(_arg IN LISTS _args)
		if(_skip_next)
			set(_skip_next FALSE)
		elseif(_arg MATCHES "^-(o|MT|MF|MQ)$")
			set(_skip_next TRUE)
		elseif(_arg MATCHES "^-(c|MD|MMD)$" OR _arg MATCHES "^-(o|MT|MF|MQ).")
		else()
			list(APPEND _preprocess "${_arg}")
		endif()
	endforeach()
	list(APPEND _preprocess -E -H -o /dev/null)

	execute_process(COMMAND ${_preprocess}
		WORKING_DIRECTORY "${_directory}"
		RESULT_VARIABLE _rc
		ERROR_VARIABLE _tree)
	math(EXPR _checked "${_checked} + 1")

	if(NOT _rc EQUAL 0)
		# The compiler's diagnostics, without -H's tree around them.
		string(REPLACE "\n" ";" _lines "${_tree}")
		list(FILTER _lines EXCLUDE REGEX "^\\.+ |^/[^ :]*$|^Multiple include guards")
		string(REPLACE ";" "\n" _errors "${_lines}")
		list(APPEND _failures "${_source}: does not preprocess:\n${_errors}")
		continue()
	endif()

	# The common case, no MFEM header anywhere, without the walk below.
	if(NOT _tree MATCHES "/mfem(\\.hpp|/[a-z_]+/[^/\n]+\\.h(pp)?)(\n|$)")
		continue()
	endif()

	# -H's lines are "<dots> <path>", one dot per level of inclusion. Keep the
	# open chain, so the first MFEM header is reported with how it was reached.
	string(REPLACE "\n" ";" _lines "${_tree}")
	set(_chain "")
	foreach(_line IN LISTS _lines)
		if(NOT _line MATCHES "^(\\.+) +(.+)$")
			continue()
		endif()
		string(LENGTH "${CMAKE_MATCH_1}" _depth)
		set(_header "${CMAKE_MATCH_2}")
		get_filename_component(_header "${_header}" ABSOLUTE BASE_DIR "${_directory}")
		list(LENGTH _chain _open)
		while(_open GREATER_EQUAL _depth)
			list(POP_BACK _chain)
			list(LENGTH _chain _open)
		endwhile()
		list(APPEND _chain "${_header}")
		_is_mfem_header("${_header}" _hit)
		if(_hit)
			string(REPLACE ";" "\n    includes " _route "${_chain}")
			list(APPEND _failures "${_source}\n    includes ${_route}")
			break()
		endif()
	endforeach()
endforeach()

if(_failures)
	list(LENGTH _failures _count)
	string(REPLACE ";" "\n\n" _report "${_failures}")
	message(FATAL_ERROR
		"${_count} file(s) of the MFEM-free half reach an MFEM header, and CI, "
		"which has none, will not compile them. Move the MFEM-dependent part "
		"behind a header of its own, as meq/StageChecksumMfem.hpp is.\n\n"
		"${_report}")
endif()

message(STATUS "MFEM-free: ${_checked} files checked, none reaches an MFEM header")
