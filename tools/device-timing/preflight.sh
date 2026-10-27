#!/usr/bin/env bash
# What this machine provides for build.sh, and what build.sh will supply itself.
# Read-only: it builds and downloads nothing. Writes preflight.txt beside itself.
#
# On Perlmutter, load the environment build.sh will use first:
#   module load PrgEnv-gnu cudatoolkit/13.2 cray-hdf5 cray-netcdf cmake python
# and any MKL module `module spider mkl` offers.

set -u
KIT=$(cd "$(dirname "$0")" && pwd)
OUT="$KIT/preflight.txt"
exec > >(tee "$OUT") 2>&1

have() { command -v "$1" >/dev/null 2>&1; }
row() { printf '  %-14s %-8s %s\n' "$1" "$2" "$3"; }

echo "MEQ device-timing preflight, $(date -Iseconds), $(hostname)"
echo
echo "== the machine"
have lscpu && lscpu | grep -E '^(Model name|Socket|Core|Thread|NUMA node\(s\))' | sed 's/^/  /'
have nvidia-smi && nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv,noheader | sed 's/^/  GPU: /'
echo
echo "== modules loaded"
(module list 2>&1 || true) | sed 's/^/  /'
echo
echo "== the build environment"

if have cmake; then row cmake ok "$(cmake --version | head -1)"; else row cmake MISSING "module load cmake"; fi
CXX=${CXX:-g++}
if have "$CXX"; then row compiler ok "$($CXX --version | head -1)"; else row compiler MISSING "$CXX not found: module load PrgEnv-gnu"; fi
if have nvcc; then
	row nvcc ok "$(nvcc --version | grep release | sed 's/.*release //')"
else
	row nvcc MISSING "module load cudatoolkit"
fi
if have python3; then row python3 ok "$(python3 --version 2>&1)"; else row python3 MISSING "module load python"; fi

if [ -n "${CUDSS_DIR:-}" ] && [ -d "$CUDSS_DIR/lib/cmake/cudss" ]; then
	row cuDSS ok "CUDSS_DIR=$CUDSS_DIR"
else
	row cuDSS fetch "build.sh downloads 0.8.0.10 for this CUDA major from developer.download.nvidia.com"
	# The NVIDIA HPC SDK ships cuDSS too; noted, not used, so the version is pinned.
	ls -d "${NVIDIA_PATH:-/nonexistent}"/math_libs/*/lib64/libcudss.so* 2>/dev/null | head -2 | sed 's/^/                    (also in the HPC SDK: /; s/$/)/'
	if have curl; then
		curl -fsI -m 10 https://developer.download.nvidia.com/compute/cudss/redist/redistrib_0.8.0.json >/dev/null \
			&& row "" "" "  ... and that host is reachable" \
			|| row "" WARN "  ... and that host is NOT reachable from here: set CUDSS_DIR"
	fi
fi

if [ -n "${MKLROOT:-}" ] && [ -f "$MKLROOT/include/mkl_pardiso.h" ]; then
	row oneMKL ok "MKLROOT=$MKLROOT"
else
	row oneMKL fetch "no MKLROOT: build.sh runs 'pip install mkl-devel' into the prefix"
	(module spider mkl 2>&1 || true) | grep -iE 'mkl|intel' | head -5 | sed 's/^/                    module: /'
fi

if [ -n "${NETCDF_DIR:-}" ] && [ -f "$NETCDF_DIR/include/netcdf.h" ]; then
	row netCDF-C ok "NETCDF_DIR=$NETCDF_DIR"
	if [ -f "$NETCDF_DIR/include/ncFile.h" ] || [ -f "$NETCDF_DIR/include/netcdf" ]; then
		row netcdf-cxx4 ok "in NETCDF_DIR"
	else
		row netcdf-cxx4 build "from the bundle, against NETCDF_DIR"
	fi
elif have nc-config; then
	row netCDF-C ok "nc-config --prefix = $(nc-config --prefix)"
	row netcdf-cxx4 build "from the bundle"
else
	row netCDF-C MISSING "module load cray-hdf5 cray-netcdf (MEQ's .nc output needs it)"
fi

SUNDIALS_DIR=${SUNDIALS_DIR:-${SUNDIALS_ROOT:-}}
# The install PREFIX, holding include/ and lib/ or lib64/. The patterns are
# anchored: a SUNDIALS without CUDA still carries the name, as
# "/* #undef SUNDIALS_CUDA_ENABLED */", and every one carries all three
# precisions in a comment.
cfg="${SUNDIALS_DIR:-}/include/sundials/sundials_config.h"
if [ -n "${SUNDIALS_DIR:-}" ] && [ -f "$cfg" ]; then
	v=$(sed -n 's/#define SUNDIALS_VERSION "\(.*\)"/\1/p' "$cfg")
	libdir=$(ls -d "$SUNDIALS_DIR"/lib64 "$SUNDIALS_DIR"/lib 2>/dev/null | while read -r d; do
		ls "$d"/libsundials_kinsol.* >/dev/null 2>&1 && echo "$d" && break; done)
	if ! grep -q '^#define SUNDIALS_CUDA_ENABLED' "$cfg"; then
		row SUNDIALS BAD "$v at $SUNDIALS_DIR was built without ENABLE_CUDA=ON"
	elif ! grep -q '^#define SUNDIALS_DOUBLE_PRECISION 1' "$cfg"; then
		row SUNDIALS BAD "$v at $SUNDIALS_DIR is not double precision"
	elif [ -z "$libdir" ]; then
		row SUNDIALS BAD "$v at $SUNDIALS_DIR has no libsundials_kinsol in lib/ or lib64/"
	else
		row SUNDIALS ok "$v, CUDA, double, at $SUNDIALS_DIR (libraries in ${libdir##*/}/)"
	fi
elif [ -n "${SUNDIALS_DIR:-}" ]; then
	row SUNDIALS BAD "SUNDIALS_DIR=$SUNDIALS_DIR has no include/sundials/sundials_config.h: it is the install prefix"
else
	row SUNDIALS build "7.5.0 with CUDA, from the bundle (or export SUNDIALS_DIR, the install prefix)"
fi
if [ -n "${BOOST_ROOT:-}" ]; then row Boost ok "BOOST_ROOT=$BOOST_ROOT"; else row Boost build "1.86 Boost.Math headers, from the bundle"; fi
if [ -n "${EIGEN3_DIR:-}" ]; then row Eigen ok "EIGEN3_DIR=$EIGEN3_DIR"; else row Eigen build "3.4.0, from the bundle"; fi

echo
echo "== what else is on offer (for whoever reads this)"
for m in cudss mkl intel sundials suitesparse boost eigen netcdf; do
	(module spider "$m" 2>&1 || true) | grep -E '^\s+[A-Za-z].*: ' | head -3 | sed "s/^/  [$m] /"
done
echo
echo "written to $OUT"
