#!/bin/sh
#
# The NPC thread scan, and it sweeps the two axes TOGETHER -- which is the
# opposite of what scan.sh does, deliberately, and the reason is the whole point
# of this script.
#
# scan.sh holds one axis at 1 while sweeping the other, because sweeping both was
# measured to be actively misleading: threaded MKL collapsed the element-local
# dense work by a large factor, which swamped every other column in the row and
# made the assembly speedup unreadable. That rule is correct FOR SERIAL ASSEMBLY
# and it is what the two axes mean there:
#
#   OMP_NUM_THREADS  drives DarcyHybridization's threaded element loops.
#   MKL_NUM_THREADS  drives UMFPACK's BLAS, PARDISO's internals, AND -- the
#                    surprise -- the element-local dense factorisations.
#
# THAT THIRD ENTRY IS WHAT CHANGES UNDER THREADED ASSEMBLY, and it changes
# because of how MKL behaves rather than because of anything MEQ or MFEM does.
# MKL suppresses its own threading inside an ACTIVE OpenMP parallel region. So
# once the element loop is itself an OpenMP region, the element-local dgetrs and
# dgemm inside it are nested, MKL runs them sequentially, and MKL_NUM_THREADS
# costs them nothing. Measured on the isolated kernels at k = 3: the same work in
# a serial loop went 0.0109 s to 0.0648 s from MKL=1 to MKL=8, while in an OpenMP
# loop it read 0.0019 s and 0.0020 s -- identical.
#
# The trace solve, meanwhile, runs on the master thread OUTSIDE any parallel
# region, so it gets all of them. That is what makes PARDISO's MKL threads
# spendable at all, and it is why the interesting configuration is OMP = MKL = N
# rather than either axis alone. CLAUDE.md's *What to do* item 0 -- "get
# ComputeH()'s element-local dense LU off threaded MKL" -- is about the serial
# mode, and threaded assembly is a route around it that needs no new code.
#
# ONE COMBINATION IS PATHOLOGICAL AND THE SCAN DELIBERATELY INCLUDES IT.
# OMP_NUM_THREADS=1 with MKL_NUM_THREADS>1 under threaded assembly: a team of
# ONE thread does not get the nested-region suppression, and the isolated
# kernels read 12.4 s against 0.069 s serial at k = 3. The binary warns about it;
# the scan runs it so the warning has a number beside it. Expect that row to be
# slow and do not read it as a defect.
#
# One process per point, as scan.sh does, because MKL fixes its threading at
# first use and an in-process sweep would measure the first setting repeatedly.
#
# Usage:  tests/performance/npc-scan.sh [build-dir] [repeats]

BUILD="${1:-build}"
REPEATS="${2:-3}"
BIN="$BUILD/tests/NpcThreadScaling"

[ -x "$BIN" ] || BIN="$BUILD/NpcThreadScaling"
if [ ! -x "$BIN" ]; then
	echo "npc-scan.sh: cannot find NpcThreadScaling under $BUILD" >&2
	exit 1
fi

echo "npc-scan.sh: $BIN, best of $REPEATS, $(nproc) cores available"

echo
echo "############################################################"
echo "# BASELINE: everything serial. What MEQ does today."
echo "############################################################"
OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 "$BIN" --repeats "$REPEATS"

echo
echo "############################################################"
echo "# TOGETHER: OMP = MKL = N. The configuration this script is for."
echo "#"
echo "# Threaded assembly gets the element loops; the trace solve gets the"
echo "# MKL threads; the element-local dense work is nested and pays nothing"
echo "# for them. If the nesting argument holds in situ, these rows improve"
echo "# monotonically and the k = 3 rows do NOT collapse."
echo "############################################################"
for T in 2 4 8 16; do
	echo
	echo "======== OMP_NUM_THREADS=$T  MKL_NUM_THREADS=$T ========"
	OMP_NUM_THREADS="$T" MKL_NUM_THREADS="$T" \
	OMP_PROC_BIND=close OMP_PLACES=cores \
	"$BIN" --repeats "$REPEATS"
done

echo
echo "############################################################"
echo "# THE CONTROL: OMP = N, MKL = 1."
echo "#"
echo "# Threaded assembly with the trace solve left sequential. The"
echo "# difference between this and the block above is what the MKL threads"
echo "# bought -- which is to say, what PARDISO bought -- and it is the only"
echo "# honest way to attribute it."
echo "############################################################"
for T in 2 4 8 16; do
	echo
	echo "======== OMP_NUM_THREADS=$T  MKL_NUM_THREADS=1 ========"
	OMP_NUM_THREADS="$T" MKL_NUM_THREADS=1 \
	OMP_PROC_BIND=close OMP_PLACES=cores \
	"$BIN" --repeats "$REPEATS"
done

echo
echo "############################################################"
echo "# THE PATHOLOGICAL ROW, run so the warning has a number."
echo "# OMP_NUM_THREADS=1 with MKL threads, under threaded assembly."
echo "############################################################"
OMP_NUM_THREADS=1 MKL_NUM_THREADS=8 "$BIN" --repeats 1 --orders 3 --sizes 32
