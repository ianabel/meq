#!/usr/bin/env bash
# The correctness gate before any timing: each case solved on the host and on
# the GPU must print the same Newton counts, support sweeps and psi_ax to every
# digit -- DeviceParity's criterion (MEASUREMENTS.md M-228). A timing of an arm
# that reached a different answer is not a timing of the same problem.
#
#   ./check.sh            mf and ff1 (a few minutes on a login node's A100)
#   CASES="mf ff1 ff2" ./check.sh
#
# Writes check.txt and check/<case>-<device>.log; exits non-zero on a mismatch.

set -uo pipefail
KIT=$(cd "$(dirname "$0")" && pwd)
# shellcheck disable=SC1091
. "$KIT/env.sh" || { echo "no env.sh: run build.sh first"; exit 1; }
# shellcheck disable=SC1091
. "$KIT/cases.sh"
CASES=${CASES:-"mf ff1"}
THREADS=${THREADS:-8}
OUT="$KIT/check"
mkdir -p "$OUT"
exec > >(tee "$KIT/check.txt") 2>&1

echo "check, $(date -Iseconds), $(hostname), $("$MEQ" --version)"
status=0
for c in $CASES; do
	for device in cpu cuda; do
		dir="$OUT/$c-$device"
		prepare_case "$c" "$dir"
		( cd "$dir" && OMP_NUM_THREADS=$THREADS MKL_NUM_THREADS=$THREADS \
			"$MEQ" --device "$device" $(case_flags "$c") "$(case_toml "$c")" \
			> "$OUT/$c-$device.log" 2>&1 )
		echo "  $c on $device: exit $?, $(grep -m1 '^MEQ: wall' "$OUT/$c-$device.log")"
	done
	# The answer's lines: the support sweeps, the Newton count, and each
	# printed label. Not the constraints or the residual table, which sit at
	# round-off and differ between PARDISO and cuDSS by design (M-229).
	answer() {
		grep -E '^ +[0-9]+ +[0-9]+ +[0-9.e+-]+ +[0-9.e+-]+ +\(|converged in|psi_ax =|psi_bnd =|the axis, as a zero' "$1" \
			| sed -E 's/ *\(via [^)]*\)//; s/, constraint .*//; s/^ +//'
	}
	if diff <(answer "$OUT/$c-cpu.log") <(answer "$OUT/$c-cuda.log") > "$OUT/$c.diff"; then
		echo "  $c: PASS -- $(grep -m1 'psi_ax =' "$OUT/$c-cuda.log" | sed 's/^ *//')"
	else
		echo "  $c: FAIL -- host and GPU disagree:"
		sed 's/^/      /' "$OUT/$c.diff"
		status=1
	fi
done
[ $status -eq 0 ] && echo "check: PASS" || echo "check: FAIL -- do not time this build; send check.txt back"
exit $status
