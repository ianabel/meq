#!/bin/sh
#
# The thread scan. One process per point, which is not fastidiousness: MKL fixes
# its threading at first use, so a sweep inside one process would measure the
# first setting several times and report it as a flat scaling curve.
#
# TWO AXES, SWEPT SEPARATELY, because they drive different things and sweeping
# them together produces a table that means nothing:
#
#   OMP_NUM_THREADS  drives DarcyHybridization's threaded element assembly.
#   MKL_NUM_THREADS  drives UMFPACK's BLAS and PARDISO's internal parallelism.
#
# Holding the other at 1 is what isolates each. Sweeping both at once was tried
# and is actively misleading -- UMFPACK's BLAS collapses under threading (see
# below) by a factor of 140, which swamps every other column in the row and
# makes the assembly speedup unreadable.
#
# MKL_THREADING_LAYER=GNU is set throughout for the reason every ctest sets it:
# /usr/lib/x86_64-linux-gnu/libblas.so.3 resolves to libmkl_rt, the dispatcher,
# which silently corrupts UMFPACK's BLAS-3 without it. Silently -- you get
# numbers, and they are wrong. A benchmark is exactly where that goes unnoticed.
#
# WHAT THE MKL AXIS IS ACTUALLY MEASURING, so nobody reads the collapse as a
# mistake: UMFPACK calls BLAS on the small dense frontal matrices of a sparse
# LU, thousands of times per factorisation, and a threaded call costs a fork, a
# join and a barrier that dwarf the arithmetic. PARDISO is a genuinely parallel
# sparse solver and should go the other way. The scan exists to show both.
#
# Usage:  tests/performance/scan.sh [build-dir] [repeats]

BUILD="${1:-build}"
REPEATS="${2:-3}"
BIN="$BUILD/tests/TraceSolverScaling"

[ -x "$BIN" ] || BIN="$BUILD/TraceSolverScaling"
if [ ! -x "$BIN" ]; then
	echo "scan.sh: cannot find TraceSolverScaling under $BUILD" >&2
	exit 1
fi

echo "scan.sh: $BIN, best of $REPEATS, $(nproc) cores available"

echo
echo "############################################################"
echo "# AXIS 1: assembly threads (OMP), MKL pinned to 1"
echo "############################################################"
for T in 1 2 4 8 16; do
	echo
	echo "======== OMP_NUM_THREADS=$T  MKL_NUM_THREADS=1 ========"
	OMP_NUM_THREADS="$T" MKL_NUM_THREADS=1 \
	MKL_THREADING_LAYER=GNU OMP_PROC_BIND=close OMP_PLACES=cores \
	"$BIN" --repeats "$REPEATS"
done

echo
echo "############################################################"
echo "# AXIS 2: direct-solver threads (MKL), OMP pinned to 1"
echo "############################################################"
for T in 1 2 4 8 16; do
	echo
	echo "======== MKL_NUM_THREADS=$T  OMP_NUM_THREADS=1 ========"
	OMP_NUM_THREADS=1 MKL_NUM_THREADS="$T" \
	MKL_THREADING_LAYER=GNU \
	"$BIN" --repeats "$REPEATS"
done
