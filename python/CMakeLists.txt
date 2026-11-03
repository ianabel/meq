# The `meq` Python package: python/meq/, and the _meq extension over meq::Run.
#
#   ${MEQ_PYTHON_ROOT}/meq/     what `import meq` finds
#     *.py                      symlinks to python/meq/ in the source
#     _meq<abi>.so              built here
#
# MEQ_PYTHON_ROOT is a sys.path entry: PYTHONPATH pointing at it is how the
# PythonDriverParity test, and anything run by hand against this build, imports
# the package. Symlinks rather than copies, so an edit to a .py file is seen
# without a rebuild; a .py file added to python/meq/ is picked up at the next
# configure. This is MaNTA's layout (../../manta/src/python/CMakeLists.txt).

# ---- the interpreter, once ------------------------------------------------
#
# THE SYSTEM PYTHON, AND NOT A VENV PREFERRED OVER IT: MEQ and MaNTA are built
# on one machine for one interpreter, so that a transport workflow imports both
# (PYTHON-PLAN.md). The headers, the ABI suffix of the module and the
# interpreter the test runs are all derived from this one find, so the three
# cannot disagree; -DPython3_EXECUTABLE= names another.
find_package(Python3 3.11 REQUIRED COMPONENTS Interpreter Development.Module)
message(STATUS "meq: the Python module is built for ${Python3_EXECUTABLE} (${Python3_VERSION})")
# For tests/, a sibling directory, where this find's variables are not in scope.
set(MEQ_PYTHON_EXECUTABLE "${Python3_EXECUTABLE}" CACHE INTERNAL
	"The interpreter the meq module is built for")

# pybind11 from the submodule, pinned to MaNTA's commit so that the two
# modules share one pybind11 internals ABI in one interpreter. It takes the
# Python3 found above rather than searching again.
set(PYBIND11_FINDPYTHON ON)
add_subdirectory("${PROJECT_SOURCE_DIR}/extern/pybind11"
	"${PROJECT_BINARY_DIR}/extern/pybind11" EXCLUDE_FROM_ALL)

# ---- meq_core inside a shared object ----------------------------------------
#
# POSITION-INDEPENDENT, BECAUSE THE MODULE IS A SHARED OBJECT AND meq_core IS
# STATIC. The compiler's default here is -fPIE, which is position independent
# for an executable but reaches thread_local storage by the local-exec model,
# and meq_core has thread_locals (CriticalPoints.cpp's pointOf(), among
# others) that a dlopen'd object cannot reach that way. -fno-semantic-
# interposition keeps -fPIC from turning every non-inline function into an
# interposable one that GCC may not inline: nothing interposes meq_core.
#
# MFEM's static archive needs nothing of the kind: it carries no 32-bit
# absolute relocations, and --exclude-libs below makes its symbols local, so
# its PC-relative ones resolve inside the module.
set_target_properties(meq_core PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_compile_options(meq_core PRIVATE
	$<$<COMPILE_LANGUAGE:CXX>:-fno-semantic-interposition>)

# ---- the extension ----------------------------------------------------------
#
# NO_EXTRAS for MaNTA's reasons: no second LTO, and no strip, because a
# backtrace from inside MKL or MFEM needs its symbols.
#
# --exclude-libs,ALL MAKES EVERY SYMBOL FROM A STATIC ARCHIVE LOCAL TO THE
# MODULE: MFEM's and meq_core's. That is what lets MFEM's PIE objects link into
# a shared object at all, and it keeps them out of the process's global symbol
# table, where a second module carrying its own MFEM or toml11 -- MaNTA's --
# would otherwise bind to whichever was loaded first.
set(MEQ_PYTHON_ROOT "${PROJECT_BINARY_DIR}/python" CACHE INTERNAL
	"The sys.path entry holding this build's meq package")
set(_meq_package_dir "${MEQ_PYTHON_ROOT}/meq")

pybind11_add_module(_meq NO_EXTRAS _meq.cpp)
target_link_libraries(_meq PRIVATE meq_core toml11::toml11)
target_link_options(_meq PRIVATE "LINKER:--exclude-libs,ALL")
# ---- MKL and SuiteSparse of its own -----------------------------------------
#
# STATIC AND HIDDEN, SO THAT NOTHING ELSE IN THE PROCESS CAN BIND TO THEM OR
# THEM TO IT. Debian's numpy loads the system MKL into the process's global
# scope; a module that links MKL shared then binds its BLAS, LAPACK and PARDISO
# to that MKL if numpy came first, and MKL's lazily loaded kernels mix with it
# even if numpy came after -- measured, the answers move about 1e-13
# (PYTHON-PLAN.md, PY-1). Linked from static archives under --exclude-libs,
# the module carries its own copy, which loads no kernels and exports nothing.
# SuiteSparse with it, because UMFPACK and CHOLMOD call BLAS and, left shared,
# would pull the shared MKL back in.
#
# HOW: MFEM's exported link line names these as -L<dir> -l<name>, and GNU ld
# takes the first -L directory that has the name. A directory of this build's
# own, searched first, answers every one of those names with a linker script
# that GROUPs all the static archives -- MKL's three depend on one another
# and on nothing in order, which a group is what resolves. The meq binary's
# link is untouched: it keeps the shared libraries, and the two are the same
# versions, which PythonDriverParity holds to the bit.
#
# OpenMP cannot be made private the same way: libgomp.a holds local-exec TLS a
# shared object cannot, and RTLD_DEEPBIND, the other route, makes libgomp read
# a stale environment and ignore OMP_NUM_THREADS. So libgomp stays shared and
# is bound when meq is imported, and meq/__init__.py refuses a process that
# already holds another OpenMP: import meq before numpy.
set(MEQ_PYTHON_MKL_STATIC_DIR "/opt/intel/oneapi/mkl/latest/lib" CACHE PATH
	"Where libmkl_intel_lp64.a, libmkl_gnu_thread.a and libmkl_core.a are, for the module")
set(MEQ_PYTHON_SUITESPARSE_STATIC_DIR "${MEQ_PROJECT_DIR}/../suitesparse/install-static/lib"
	CACHE PATH "A static, position-independent SuiteSparse, for the module")
set(_meq_static_names mkl_intel_lp64 mkl_gnu_thread mkl_core)
set(_meq_suitesparse_names umfpack klu amd btf cholmod colamd camd ccolamd suitesparseconfig)
set(_meq_static_archives "")
foreach(_name IN LISTS _meq_static_names)
	list(APPEND _meq_static_archives "${MEQ_PYTHON_MKL_STATIC_DIR}/lib${_name}.a")
endforeach()
foreach(_name IN LISTS _meq_suitesparse_names)
	list(APPEND _meq_static_archives "${MEQ_PYTHON_SUITESPARSE_STATIC_DIR}/lib${_name}.a")
endforeach()
set(_meq_static_missing "")
foreach(_archive IN LISTS _meq_static_archives)
	if(NOT EXISTS "${_archive}")
		list(APPEND _meq_static_missing "${_archive}")
	endif()
endforeach()
if(_meq_static_missing)
	message(STATUS "meq: the Python module links MKL and SuiteSparse SHARED, since "
		"these are missing: ${_meq_static_missing}. It then shares them with "
		"anything else in the process (python/CMakeLists.txt says what that costs)")
	set(MEQ_PYTHON_STATIC_LIBRARIES OFF CACHE INTERNAL "")
else()
	set(_meq_static_view "${PROJECT_BINARY_DIR}/python-static-libraries")
	file(MAKE_DIRECTORY "${_meq_static_view}")
	list(JOIN _meq_static_archives " " _meq_group)
	foreach(_name IN LISTS _meq_static_names _meq_suitesparse_names)
		file(WRITE "${_meq_static_view}/lib${_name}.a"
			"/* python/CMakeLists.txt: every static archive _meq links, as one group */\n"
			"GROUP ( ${_meq_group} )\n")
	endforeach()
	message(STATUS "meq: the Python module links MKL and SuiteSparse static and hidden")
	set(MEQ_PYTHON_STATIC_LIBRARIES ON CACHE INTERNAL "")
	# First among the -L directories, so every MKL and SuiteSparse name finds
	# its GROUP script here before the shared library elsewhere.
	target_link_options(_meq BEFORE PRIVATE "-L${_meq_static_view}")
endif()

# meq_version.hpp: the module stamps the commit it was linked at, as meq does.
target_include_directories(_meq PRIVATE "${PROJECT_BINARY_DIR}")
add_dependencies(_meq meq_version_header)
set_target_properties(_meq PROPERTIES
	LIBRARY_OUTPUT_DIRECTORY "${_meq_package_dir}"
	INSTALL_RPATH_USE_LINK_PATH ON)

# WHICH LIBRARIES THE MODULE RESOLVES TO, recorded after every link, so that
# `import meq` can refuse a process that already holds another file of the same
# name -- another MKL, another SUNDIALS (meq/__init__.py says why).
add_custom_command(TARGET _meq POST_BUILD
	COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/record_libraries.py"
		"$<TARGET_FILE:_meq>" "${_meq_package_dir}/_libraries.py"
	COMMENT "Recording the libraries _meq resolves to"
	VERBATIM)

file(GLOB _meq_python_sources CONFIGURE_DEPENDS
	"${CMAKE_CURRENT_SOURCE_DIR}/meq/*.py")
file(MAKE_DIRECTORY "${_meq_package_dir}")
foreach(_source IN LISTS _meq_python_sources)
	get_filename_component(_name "${_source}" NAME)
	file(CREATE_LINK "${_source}" "${_meq_package_dir}/${_name}" SYMBOLIC)
endforeach()

# ---- the install ------------------------------------------------------------
#
# lib/python3.X/site-packages under the prefix, CPython's own layout for one;
# install/<name> is not on any interpreter's path, so a user puts this on
# PYTHONPATH, as for the build tree.
set(MEQ_INSTALL_PYTHONDIR
	"lib/python${Python3_VERSION_MAJOR}.${Python3_VERSION_MINOR}/site-packages"
	CACHE STRING "Where the meq package is installed, relative to the prefix")
install(FILES ${_meq_python_sources} "${_meq_package_dir}/_libraries.py"
	DESTINATION "${MEQ_INSTALL_PYTHONDIR}/meq")
install(TARGETS _meq LIBRARY DESTINATION "${MEQ_INSTALL_PYTHONDIR}/meq")
