# The SharedTestCore ctest: one copy of MFEM, and of meq, per test process.
#
# Every test executable links meq_core_tests, a shared library carrying meq's
# objects and all of libmfem (CMakeLists.txt). A test that ALSO links
# libmfem.a or libmeq_core.a ahead of it -- MFEM::MFEM or meq_core -- takes
# its own copies of whatever it calls first, and so has two of MFEM's globals,
# or meq's, in the process: a silent ODR violation the linker does not report.
# What gives it away is object code: a non-inline function is a strong global
# text symbol (nm's `T`), whereas what a test instantiates from the headers is
# inline or a template, and weak (`W`). So:
#
#   * the shared library defines MFEM's functions -- a few thousand, or it is
#     not the library this check is about;
#   * no test executable defines, strong, an mfem:: or meq:: function the
#     library also defines. (A test may define meq:: functions of its own --
#     ConfigTests' operator<< for the Config enums -- and those are not in it.)
#
# Run as cmake -P with NM, CORE (the library) and LIST (a file of executables,
# one per line).

foreach(_var NM CORE LIST)
	if(NOT DEFINED ${_var})
		message(FATAL_ERROR "CheckTestLinkage.cmake needs -D${_var}=")
	endif()
endforeach()

# The mangled names of the strong text symbols in namespaces mfem and meq.
function(_strong_text _file _out)
	execute_process(COMMAND "${NM}" --defined-only "${_file}"
		OUTPUT_VARIABLE _symbols RESULT_VARIABLE _rc)
	if(NOT _rc EQUAL 0)
		message(FATAL_ERROR "${NM} failed on ${_file}")
	endif()
	string(REGEX MATCHALL "[0-9a-f]+ T _ZNK?(4mfem|3meq)[^\n]*" _hits "${_symbols}")
	list(TRANSFORM _hits REPLACE "^[0-9a-f]+ T " "")
	set(${_out} "${_hits}" PARENT_SCOPE)
endfunction()

_strong_text("${CORE}" _core)
set(_core_mfem "${_core}")
list(FILTER _core_mfem INCLUDE REGEX "^_ZNK?4mfem")
list(LENGTH _core_mfem _core_count)
if(_core_count LESS 1000)
	message(FATAL_ERROR "${CORE} defines ${_core_count} of MFEM's functions; "
		"it should carry all of libmfem")
endif()
message(STATUS "${CORE}: ${_core_count} MFEM functions")

file(STRINGS "${LIST}" _executables)
set(_bad 0)
list(LENGTH _executables _checked)
foreach(_exe IN LISTS _executables)
	_strong_text("${_exe}" _hits)
	set(_twice "")
	foreach(_symbol IN LISTS _hits)
		list(FIND _core "${_symbol}" _index)
		if(NOT _index EQUAL -1)
			list(APPEND _twice "${_symbol}")
		endif()
	endforeach()
	list(LENGTH _twice _count)
	if(_count GREATER 0)
		list(GET _twice 0 _first)
		message(SEND_ERROR "${_exe} defines ${_count} functions ${CORE} also "
			"defines (${_first}, ...): it links a second copy of MFEM or meq")
		math(EXPR _bad "${_bad} + 1")
	endif()
endforeach()
if(_bad GREATER 0)
	message(FATAL_ERROR "${_bad} of ${_checked} test executables carry their own MFEM or meq")
endif()
message(STATUS "${_checked} test executables, none carrying MFEM's or meq's object code")
