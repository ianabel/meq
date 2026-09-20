#!/bin/sh
# Half a dozen MACHINE-RELEVANT FIXED-BOUNDARY cases -> examples/.
#
#     sh tools/freegs4e-benchmark/make_fixed.sh [examples-dir] [reference-dir]
#
# THE RECIPE IS IAN'S: solve a machine FREE boundary in freegs4e, fit MXH to its
# psi_n = 0.95 surface, and hand MEQ that curve with psi = 0 on it.  psi is
# defined only up to a constant, so the result is the SAME equilibrium in a
# different gauge -- a fixed-boundary problem with real machine geometry, real
# tabulated profiles, and an independent code's answer to check against.
# make_case.py is the generator and its docstring is the argument.
#
# THE REFERENCES ARE NOT IN THE REPOSITORY AT THIS RESOLUTION AND DO NOT NEED TO
# BE.  The shipped TOML carries the MXH coefficients and the two profile tables,
# which is everything a run needs; the .npz is needed only to COMPARE against
# freegs4e afterwards.  Rebuild them with
#
#     FGSREF_OUT=/some/where PYTHONPATH=../freegs4e \
#         venv/bin/python fgsref.py --shaped --nx=257 --seed-from=auto
#
# which is about six minutes for all nine machines because --seed-from spends
# the grid nesting.  257^2 rather than the committed 129^2 because the MXH fit
# is extracted from that grid and IS the floor of the comparison: at 129^2 with
# 10 harmonics it reads 2.7e-04 to 7.6e-04 m, and at 257^2 with 20 it reads
# 4.4e-05 to 1.8e-04.  Above 20 harmonics nothing improves, which is what says
# the remaining residual is the grid rather than the truncation.
#
# WHY THESE SIX.  A range of SHAPES and a range of SIZES, because both were
# missing: every machine case in this tree is free boundary, and five of the
# eight shipped meshes are under 4096 elements.
#
#   h  limited circular   kappa 1.01, delta 0.05   the nearly-circular control,
#                                                  and LIMITED, so no X-point
#                                                  exists anywhere in it
#   a  TestTokamak        kappa 1.19, delta 0.33   the plain diverted case, and
#                                                  the same machine that
#                                                  examples/machine-a-*.toml
#                                                  solves FREE boundary
#   e  TestTokamak diam.  kappa 1.31, delta 0.35   exactly up-down symmetric,
#                                                  and gg' of the other sign
#   d  TCV                kappa 1.75, delta 0.23   the most elongated
#   g  MAST-U             aspect 1.56              spherical, the lowest aspect
#   f  DIII-D             kappa 1.47, a 0.59 m     the large conventional one,
#                                                  near double null
set -e
here=$(cd "$(dirname "$0")" && pwd)
out=${1:-$here/../../examples}
refs=${2:-$here/ref-n257}

if [ ! -d "$refs" ]; then
	echo "no reference directory $refs -- see the header for how to build it" >&2
	exit 2
fi

# stem                     reference                          elements
for row in \
	"fixed-h-circular        H_limited_circular_shaped          1000" \
	"fixed-a-testtokamak     A_testtokamak_classic_shaped       2500" \
	"fixed-e-diamagnetic     E_testtokamak_diamagnetic_shaped   6000" \
	"fixed-d-tcv             D_tcv_conventional_shaped         15000" \
	"fixed-g-mastu           G_mastu_simple_shaped             36000" \
	"fixed-f-diiid           F_diiid_conventional_shaped       90000"
do
	set -- $row
	"$here/venv/bin/python" "$here/make_case.py" "$refs/$2.npz" \
		-o "$out" --stem "$1" --elements "$3"
done
