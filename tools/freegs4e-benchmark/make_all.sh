#!/bin/sh
# Every diverted freegs4e reference -> a MEQ free-boundary configuration.
#
# ONE STEM PER REFERENCE, and the letter is the reference's own, so a row of a
# comparison table and a file in examples/ name the same machine.  A is
# TestTokamak classic, which examples/diverted-tokamak*.toml also carries by
# hand -- kept, because the two agreeing is the check that this script is
# writing the same machine somebody once wrote out themselves.
#
#     sh tools/freegs4e-benchmark/make_all.sh [examples-dir]
#
# Re-runnable: it overwrites the TOML, the two profile tables and the guess for
# each case, all of which are derived from the .npz beside this script.
set -e
here=$(cd "$(dirname "$0")" && pwd)
out=${1:-$here/../../examples}

for pair in \
	"A_testtokamak_classic machine-a-testtokamak" \
	"B_testtokamak_peaked_ffdom machine-b-ffprime" \
	"C_mast_spherical machine-c-mast" \
	"D_tcv_conventional machine-d-tcv" \
	"E_testtokamak_diamagnetic machine-e-diamagnetic" \
	"F_diiid_conventional machine-f-diiid" \
	"G_mastu_simple machine-g-mastu"
do
	set -- $pair
	ref=$1
	stem=$2
	echo "============================================================"
	echo "$ref -> $stem"
	"$here/venv/bin/python" "$here/make_diverted_case.py" \
		"$here/$ref.npz" "$here/$ref.json" "$out" "$stem"
done
