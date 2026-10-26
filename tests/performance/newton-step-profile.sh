#!/bin/sh
#
# The Newton-step leg profile, run the way the measurement has to be run.
#
# THREE THINGS THIS SCRIPT IS FOR, all of them named by
# HDG-NEWTON-STEP-PROFILE-FROM-HDGDEV.md as ways the measurement can be spoiled:
#
#   MKL_NUM_THREADS IS PINNED TO 1. Two reasons, and they are different reasons.
#   MEQ's is that the SERIAL assembly rows would otherwise collapse -- the
#   element-local dense LU sits outside any OpenMP region there and pays full MKL
#   threading per call, measured at forty times slower at k = 3, which would make
#   the gradient leg's share a measurement of MKL's threading threshold rather
#   than of ComputeH(). Upstream's is that omp_set_num_threads() retunes MKL, so
#   a leg share taken while both axes are moving measures their interaction.
#
#   MEDIANS OF FIVE. The binary's default; stated here so a reader of the script
#   does not have to know it. A share from one run is not a share.
#
#   ONE PROCESS PER THREAD SETTING, because MKL fixes its threading at first use
#   and an in-process sweep would measure the first setting repeatedly.
#
# It does NOT gate on an idle machine, deliberately: that belongs to whoever runs
# it, and a script that sleeps for four load samples is a script people work
# around. Read /proc/loadavg first, or run it under a gate of your own. A leg
# SHARE is more robust to contention than a wall time is -- every leg is slowed
# together -- but the serial/threaded comparison is not.
#
# Usage:  tests/performance/newton-step-profile.sh [build-dir] [repeats] [source-extra]
#
# source-extra, when given, is handed to the binary as --source-extra: the
# order added to 2k for meq::SourceIntegrator's rule, INTERPOLATORY-HDG-PLAN.md's
# IH-0 knob. Absent, the binary leaves the solver's default alone.
#
# Runs from any working directory. build-dir defaults to the project's
# builds/double, found from this script's own location (source/tests/performance
# under the project); a relative build-dir is taken from the working directory.
# Nothing is written: the binary reads no file and prints its tables.

HERE=$(cd "$(dirname "$0")" && pwd -P)
PROJECT=$(cd "$HERE/../../.." && pwd -P)
BUILD="${1:-$PROJECT/builds/double}"
REPEATS="${2:-5}"
SOURCE_EXTRA=""
[ -n "$3" ] && SOURCE_EXTRA="--source-extra $3"
BIN="$BUILD/tests/NewtonStepProfile"

[ -x "$BIN" ] || BIN="$BUILD/NewtonStepProfile"
if [ ! -x "$BIN" ]; then
	echo "newton-step-profile.sh: cannot find NewtonStepProfile under $BUILD" >&2
	exit 1
fi

echo "newton-step-profile.sh: $BIN, median of $REPEATS, $(nproc) cores"
echo "loadavg at start: $(cat /proc/loadavg)"

# OMP_NUM_THREADS drives the threaded rows; the serial rows ignore it. Both are
# run in ONE process because AssemblyMode is a per-solver setting rather than an
# environment one -- unlike MKL, which is why that one gets a pinned value and
# not a sweep.
echo
echo "############################################################"
echo "# OMP_NUM_THREADS=8, MKL_NUM_THREADS=1 -- MEQ's own configuration"
echo "############################################################"
OMP_NUM_THREADS=8 MKL_NUM_THREADS=1 \
OMP_PROC_BIND=close OMP_PLACES=cores \
"$BIN" --repeats "$REPEATS" $SOURCE_EXTRA

# THE SECOND ARM ANSWERS A DIFFERENT QUESTION AND IS THREADED-ONLY ON PURPOSE.
# Pinning MKL to 1 is right for the serial rows -- it stops the gradient leg
# being a measurement of MKL's threading threshold, which is the factor of forty
# at k = 3 recorded in CLAUDE.md -- but it ALSO stops UMFPACK's dense frontal
# kernels threading, so the trace solve's share in the pinned arm is an UPPER
# bound on its production share rather than a production number. Under threaded
# assembly the element-local work is nested in an active OpenMP region where MKL
# suppresses its own threading anyway, so this arm moves the trace solve and
# very little else. It is the arm that says whether Amdahl caps the element-local
# programme. The serial rows are NOT re-run unpinned: there the element loop is
# outside every region and the collapse would swamp the thing being measured.
echo
echo "############################################################"
echo "# OMP_NUM_THREADS=8, MKL_NUM_THREADS UNPINNED -- threaded rows only"
echo "############################################################"
env -u MKL_NUM_THREADS \
OMP_NUM_THREADS=8 OMP_PROC_BIND=close OMP_PLACES=cores \
"$BIN" --repeats "$REPEATS" --threaded $SOURCE_EXTRA

echo
echo "loadavg at end: $(cat /proc/loadavg)"
