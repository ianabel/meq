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
#     sh tools/freegs4e-benchmark/make_all.sh [examples-dir] --shaped
#
# Re-runnable: it overwrites the TOML, the two profile tables and the guess for
# each case, all of which are derived from the .npz beside this script.
#
# --shaped READS THE OTHER SET OF REFERENCES and writes the other set of stems:
# <ref>_shaped.npz -> machine-x-...-shaped.*, the machines whose conductors are
# ShapedCoils on MEQ's own rectangles rather than filaments.  Those references
# come from `fgsref.py --shaped` and have to exist first; see README.md.  The
# two sets are tracked side by side in examples/ because they describe
# DIFFERENT MACHINES -- freegs4e's inverse solve picks different coil currents
# once the conductors have extent -- so neither is derivable from the other.
set -e
here=$(cd "$(dirname "$0")" && pwd)
out=${1:-$here/../../examples}
suffix=""
case "${2:-}" in
	--shaped) suffix="_shaped" ;;
	"") ;;
	*) echo "unknown option: $2" >&2; exit 2 ;;
esac

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
	echo "$ref$suffix -> $stem${suffix:+-shaped}"
	"$here/venv/bin/python" "$here/make_diverted_case.py" \
		"$here/$ref$suffix.npz" "$here/$ref$suffix.json" "$out" \
		"$stem${suffix:+-shaped}"
done
