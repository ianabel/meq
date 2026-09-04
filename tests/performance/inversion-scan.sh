#!/bin/sh
#
# The solution-inversion scan, INVERSION-PLAN.md stage IN-P.
#
# ONE PROCESS PER THREAD COUNT, which is not fastidiousness and is the same
# reason scan.sh gives beside it: MKL fixes its threading at first use and the
# OpenMP runtime forks its team once, so a sweep inside one process would measure
# the first setting several times and report it as a flat scaling curve.
#
# MKL_NUM_THREADS=1 THROUGHOUT AND IT IS NOT AN AXIS HERE. Everything measured
# in this file is OpenMP over independent work -- surfaces, rays, psi DOFs --
# and none of it goes near a BLAS call worth threading. CLAUDE.md records
# ComputeH()'s element-local dense LU degrading by a factor of forty at k = 3
# under threaded MKL, and the variable is process-wide, so raising it here would
# make the SOLVE that feeds every one of these measurements slower and change
# nothing about the extraction. If a run appears to gain from it, that is what
# is happening.
#
# WHY THE THREAD SECTION BUILDS ONE SOLVE PER THREAD, which is the first thing
# a reader of the output should know: mfem::Mesh::FindPoints is not reentrant --
# it loops over every element through the SHARED GetElementTransformation( int )
# overload -- and every entry point of the tracer reaches it. So the sections
# below measure the parallelism that is AVAILABLE, with per-thread state, rather
# than pretending a shared const ContourTracer is safe. It is not: the first
# version of that section aborted.
#
# Usage:  tests/performance/inversion-scan.sh [build-dir] [repeats]

BUILD="${1:-build}"
REPEATS="${2:-2}"
BIN="$BUILD/tests/InversionScaling"

[ -x "$BIN" ] || BIN="$BUILD/InversionScaling"
if [ ! -x "$BIN" ]; then
	echo "inversion-scan.sh: cannot find InversionScaling under $BUILD" >&2
	exit 1
fi

echo "inversion-scan.sh: $BIN, best of $REPEATS, $(nproc) cores available"

echo
echo "############################################################"
echo "# The breakdown, the scaling and the kernels -- all serial."
echo "# These are the numbers a design decision is traded against."
echo "############################################################"
OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 "$BIN" --sections 1,2,3 --repeats "$REPEATS"
STATUS=$?

echo
echo "############################################################"
echo "# The algorithmic levers: continuation against parallelism,"
echo "# the pointwise call pattern, and dGeometry_dpsi."
echo "############################################################"
OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 "$BIN" --sections 5 --repeats "$REPEATS"
[ $? -eq 0 ] || STATUS=1

echo
echo "############################################################"
echo "# Thread scaling, over surfaces and over rays."
echo "# The internal sweep is over num_threads( t ) inside ONE"
echo "# process, which is legitimate HERE and nowhere else in this"
echo "# directory: no MKL and no lazily built global table is"
echo "# reached inside the parallel regions -- the serial reference"
echo "# pass runs first and warms every one of them."
echo "############################################################"
MKL_NUM_THREADS=1 OMP_PROC_BIND=close OMP_PLACES=cores \
	"$BIN" --sections 4
[ $? -eq 0 ] || STATUS=1

echo
if [ $STATUS -eq 0 ]; then
	echo "inversion-scan.sh: every correctness property held."
else
	echo "inversion-scan.sh: a correctness property FAILED. The timings above" \
	     "are not a result until that is fixed." >&2
fi
exit $STATUS
