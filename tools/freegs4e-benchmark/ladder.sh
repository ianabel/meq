#!/bin/bash
# Refine the freegs4e REFERENCE, which is what tools/README.md says the next
# work on this benchmark is: MEQ saturates the comparison at about 1.4e-04
# because that is the reference's own accuracy, not MEQ's.
#
# Runs one case over a nested ladder of grids and records wall clock, peak
# resident memory and the Picard count at each, so that "how does freegs4e
# behave under pressure" is answered with numbers rather than with the two
# endpoints.  The grids are 2^n + 1 so each is a point-for-point subset of the
# next, which is what a reference refinement study needs.
#
#   FGSREF_OUT=<dir> ladder.sh <case-prefix> [n ...]
set -u
CASE="${1:-A_}"
shift || true
LEVELS=${*:-"129 257 513 1025 2049"}
HERE="$( cd "$( dirname "$0" )" && pwd )"
OUT="${FGSREF_OUT:?set FGSREF_OUT to a scratch directory}"
VENV="${FGS_VENV:?set FGS_VENV to the python that can import freegs4e}"
FGS="${FREEGS4E:-/home/ian/projects/freegs4e}"

mkdir -p "$OUT"
printf "%6s %12s %10s %10s %8s\n" n wall_s peak_MB picard status | tee "$OUT/ladder.txt"
for n in $LEVELS; do
   log="$OUT/run-n$n.log"
   /usr/bin/time -v -o "$OUT/time-n$n.txt" \
      env FGSREF_OUT="$OUT" PYTHONPATH="$FGS" "$VENV" "$HERE/fgsref.py" \
      --nx=$n "$CASE" > "$log" 2>&1
   rc=$?
   wall=$(awk -F': ' '/Elapsed \(wall clock\)/{print $2}' "$OUT/time-n$n.txt")
   peak=$(awk -F': ' '/Maximum resident set size/{printf "%.0f", $2/1024}' "$OUT/time-n$n.txt")
   pic=$(grep -oE "picard[^0-9]*[0-9]+" "$log" | tail -1 | grep -oE "[0-9]+$")
   st=$(grep -E "^(OK|CHECK|FAILED)" "$log" | tail -1 | awk '{print $1}')
   printf "%6s %12s %10s %10s %8s\n" "$n" "${wall:-?}" "${peak:-?}" "${pic:-?}" "${st:-rc$rc}" \
      | tee -a "$OUT/ladder.txt"
done
