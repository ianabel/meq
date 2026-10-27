#!/usr/bin/env bash
# The timing matrix: every case on every arm, in interleaved rounds, each run
# with --profile. Run it on ONE EXCLUSIVE GPU NODE (job.slurm), after check.sh
# has passed on the same build.
#
#   CASES        default "mf ff1 ff2 ff3"
#   ARMS         default "cpu8 cpu64 cuda": cpuN is the host at N threads
#                (OMP = MKL = N), cuda is one GPU
#   CUDA_THREADS host threads beside the GPU (default 8; MEQ's host loops are
#                serial under a Device either way, M-188)
#   ROUNDS       default 3; FF3_ROUNDS default 2 (ff3 is minutes a run)
#   NSYS=1       afterwards, one cuda run per case under nsys, for the bus
#
# Writes results/<case>/<arm>/r<k>/{run.log,exit,time.txt} and results/env.txt,
# then results-<host>-<date>.tar.gz to bring back.

set -uo pipefail
KIT=$(cd "$(dirname "$0")" && pwd)
# shellcheck disable=SC1091
. "$KIT/env.sh" || { echo "no env.sh: run build.sh first"; exit 1; }
# shellcheck disable=SC1091
. "$KIT/cases.sh"
CASES=${CASES:-"mf ff1 ff2 ff3"}
ARMS=${ARMS:-"cpu8 cpu64 cuda"}
CUDA_THREADS=${CUDA_THREADS:-8}
ROUNDS=${ROUNDS:-3}
FF3_ROUNDS=${FF3_ROUNDS:-2}
RESULTS="$KIT/results"
mkdir -p "$RESULTS"
# GNU time for peak memory and CPU seconds, where the node has it.
TIMER=()
[ -x /usr/bin/time ] && TIMER=(/usr/bin/time -v -o time.txt)

# What the numbers were taken on.
{
	echo "date        $(date -Iseconds)"
	echo "host        $(hostname)"
	echo "meq         $("$MEQ" --version)"
	echo "job         ${SLURM_JOB_ID:-none} nodes=${SLURM_JOB_NODELIST:-} cpus=${SLURM_CPUS_ON_NODE:-}"
	echo "arms        $ARMS (cuda with $CUDA_THREADS host threads); rounds $ROUNDS, ff3 $FF3_ROUNDS"
	echo
	lscpu 2>/dev/null | grep -E '^(Model name|Socket|Core|Thread|NUMA node)'
	echo
	nvidia-smi --query-gpu=index,name,memory.total,driver_version,clocks.max.sm,power.limit --format=csv 2>/dev/null
	echo
	(module list 2>&1 || true)
	echo
	cat "$KIT/PROVENANCE.txt" 2>/dev/null
} > "$RESULTS/env.txt"

run_one() # case, arm, round
{
	local c=$1 arm=$2 r=$3 dir device threads
	dir="$RESULTS/$c/$arm/r$r"
	case "$arm" in
		cpu*) device=cpu; threads=${arm#cpu} ;;
		cuda) device=cuda; threads=$CUDA_THREADS ;;
		*) echo "unknown arm $arm"; return 1 ;;
	esac
	prepare_case "$c" "$dir"
	local start=$SECONDS
	( cd "$dir" && OMP_NUM_THREADS=$threads MKL_NUM_THREADS=$threads OMP_PROC_BIND=close OMP_PLACES=cores \
		"${TIMER[@]}" "$MEQ" --device "$device" --profile $(case_flags "$c") \
		"$(case_toml "$c")" > run.log 2>&1 )
	local code=$?
	echo "$code" > "$dir/exit"
	# The products are large and nobody reads them here; the logs are the data.
	rm -f "$dir"/*.nc "$dir"/*.gf "$dir"/*.msh
	printf '  %-4s %-6s r%d  exit %d  %4ds  %s\n' "$c" "$arm" "$r" "$code" $((SECONDS - start)) \
		"$(grep -m1 '^MEQ: wall' "$dir/run.log" | sed 's/^MEQ: //')"
}

echo "timing: cases [$CASES], arms [$ARMS], into $RESULTS"
for r in $(seq 1 "$ROUNDS"); do
	echo "round $r"
	for c in $CASES; do
		if [ "$c" = ff3 ] && [ "$r" -gt "$FF3_ROUNDS" ]; then continue; fi
		for arm in $ARMS; do
			run_one "$c" "$arm" "$r"
		done
	done
done

if [ "${NSYS:-0}" = 1 ] && command -v nsys >/dev/null; then
	echo "nsys: one cuda run per case"
	for c in $CASES; do
		dir="$RESULTS/nsys/$c"
		prepare_case "$c" "$dir"
		( cd "$dir" && OMP_NUM_THREADS=$CUDA_THREADS MKL_NUM_THREADS=$CUDA_THREADS \
			nsys profile --trace=cuda --force-overwrite true -o "$c" \
			"$MEQ" --device cuda $(case_flags "$c") "$(case_toml "$c")" > run.log 2>&1 \
			&& nsys export --type sqlite --force-overwrite true -o "$c.sqlite" "$c.nsys-rep" > /dev/null 2>&1 )
		rm -f "$dir"/*.nc "$dir"/*.gf "$dir"/*.msh
		echo "  $c: $(ls "$dir" | tr '\n' ' ')"
	done
fi

bundle="$KIT/results-$(hostname -s)-$(date +%Y%m%d-%H%M).tar.gz"
tar -czf "$bundle" -C "$KIT" results $(cd "$KIT" && ls preflight.txt check.txt check 2>/dev/null)
echo "bring back: $bundle"
command -v python3 >/dev/null && python3 "$KIT/summarize.py" "$RESULTS"
