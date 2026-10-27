#!/usr/bin/env bash
# Build SUNDIALS, MFEM (meq-integration) and MEQ for a device timing run, from
# the bundle, into ./install. Run from anywhere; everything is relative to the
# kit. Each step logs to build-logs/<step>.log and is skipped when its install
# already exists, so a failed build resumes where it stopped.
#
#   PREFIX     install root                      (default: <kit>/install)
#   CUDA_ARCH  compute capability, no dot        (default: 80, the A100)
#   JOBS       parallel compile jobs             (default: 16)
#   CXX, CC    host compilers                    (default: g++, gcc)
#   CUDSS_DIR, MKLROOT, NETCDF_DIR, BOOST_ROOT, EIGEN3_DIR, SUNDIALS_DIR (or SUNDIALS_ROOT)
#              use the machine's instead of fetching or building one.
#              A SUNDIALS_DIR must be built with ENABLE_CUDA=ON (for this
#              CUDA_ARCH), SUNDIALS_PRECISION=DOUBLE, SUNDIALS_INDEX_SIZE=64 and
#              KINSOL; MFEM refuses a CUDA build against one without CUDA.
#
# On Perlmutter:  module load PrgEnv-gnu cudatoolkit/13.2 cray-hdf5 cray-netcdf cmake python

set -euo pipefail
KIT=$(cd "$(dirname "$0")" && pwd)
SRC="$KIT/src"
PREFIX=${PREFIX:-$KIT/install}
BUILD=${BUILD:-$KIT/build}
LOGS="$KIT/build-logs"
CUDA_ARCH=${CUDA_ARCH:-80}
JOBS=${JOBS:-16}
export CXX=${CXX:-g++} CC=${CC:-gcc}
mkdir -p "$PREFIX" "$BUILD" "$LOGS"

say() { echo "[build] $*"; }
die() { echo "[build] FAILED: $*" >&2; exit 1; }
step() # name, then the command; the log is the record
{
	local name=$1; shift
	say "$name ..."
	if ! "$@" > "$LOGS/$name.log" 2>&1; then
		tail -30 "$LOGS/$name.log" >&2
		die "$name -- see $LOGS/$name.log"
	fi
}
for tool in cmake nvcc "$CXX" python3 curl tar; do
	command -v "$tool" >/dev/null || die "$tool not found (see preflight.sh)"
done

CUDA_VERSION=$(nvcc --version | sed -n 's/.*release \([0-9][0-9]*\)\.\([0-9]*\).*/\1.\2/p')
CUDA_MAJOR=${CUDA_VERSION%%.*}
say "nvcc $CUDA_VERSION, sm_$CUDA_ARCH, $CXX $($CXX -dumpfullversion), prefix $PREFIX"

# --------------------------------------------------------------------------
# cuDSS: the machine's, or NVIDIA's redistributable for this CUDA major.
# --------------------------------------------------------------------------
if [ -z "${CUDSS_DIR:-}" ]; then
	CUDSS_DIR="$PREFIX/cudss"
	if [ ! -d "$CUDSS_DIR/lib/cmake/cudss" ]; then
		case "$CUDA_MAJOR" in
			12) sum=29716a05c64b28531309afcb7a9dd5044081a43197e72bdd0fbf29e9b4ee4707 ;;
			13) sum=ba18f5fd80dcbbe905d158caac5b3061d848442bb5abd477b5f296b4257a4937 ;;
			*) die "no cuDSS 0.8.0 build for CUDA $CUDA_MAJOR; set CUDSS_DIR" ;;
		esac
		archive="libcudss-linux-x86_64-0.8.0.10_cuda${CUDA_MAJOR}-archive.tar.xz"
		url="https://developer.download.nvidia.com/compute/cudss/redist/libcudss/linux-x86_64/$archive"
		say "fetching cuDSS 0.8.0.10 (CUDA $CUDA_MAJOR)"
		curl -fsSL -o "$BUILD/$archive" "$url" || die "cannot fetch $url; set CUDSS_DIR"
		echo "$sum  $BUILD/$archive" | sha256sum -c - >/dev/null || die "cuDSS archive checksum mismatch"
		mkdir -p "$CUDSS_DIR"
		tar -xJf "$BUILD/$archive" -C "$CUDSS_DIR" --strip-components=1
	fi
fi
[ -d "$CUDSS_DIR/lib/cmake/cudss" ] || die "no lib/cmake/cudss under CUDSS_DIR=$CUDSS_DIR"
CUDSS_LIB=$(dirname "$(find "$CUDSS_DIR" -name 'libcudss.so*' | head -1)")

# --------------------------------------------------------------------------
# oneMKL: PARDISO is the host trace solver the comparison is made against,
# on MKL's GNU OpenMP layer (CLAUDE_MFEM.md, "Rebuilding MFEM").
# --------------------------------------------------------------------------
if [ -z "${MKLROOT:-}" ]; then
	MKLROOT="$PREFIX/mkl"
	if [ ! -f "$MKLROOT/include/mkl_pardiso.h" ]; then
		# mkl-devel pulls mkl and mkl-include: lib/libmkl_*.so.N and include/.
		# Pinned to the MKL MEQ's own measurements are taken on (M-303).
		step mkl python3 -m pip install --no-warn-script-location --prefix "$MKLROOT" mkl-devel==2026.2.0
	fi
fi
[ -f "$MKLROOT/include/mkl_pardiso.h" ] || die "no mkl_pardiso.h under MKLROOT=$MKLROOT"
if [ -d "$MKLROOT/lib/intel64" ]; then MKL_LIBDIR=lib/intel64; else MKL_LIBDIR=lib; fi
MKL_LIB="$MKLROOT/$MKL_LIBDIR"
# A pip install ships libmkl_*.so.N only, without the unversioned link names.
for lib in mkl_core mkl_intel_lp64 mkl_gnu_thread; do
	if [ ! -e "$MKL_LIB/lib$lib.so" ]; then
		versioned=$(ls "$MKL_LIB"/lib$lib.so.* 2>/dev/null | sort -V | tail -1)
		[ -n "$versioned" ] || die "no lib$lib.so in $MKL_LIB"
		ln -s "$(basename "$versioned")" "$MKL_LIB/lib$lib.so"
	fi
done
BLAS="-L$MKL_LIB;-Wl,-rpath,$MKL_LIB;-lmkl_intel_lp64;-lmkl_gnu_thread;-lmkl_core;-lgomp;-lpthread;-lm;-ldl"

# --------------------------------------------------------------------------
# Header libraries.
# --------------------------------------------------------------------------
if [ -z "${EIGEN3_DIR:-}" ]; then
	if [ ! -d "$PREFIX/eigen" ]; then
		tar -xzf "$SRC/eigen-3.4.0.tar.gz" -C "$BUILD"
		step eigen sh -c "cmake -S '$BUILD/eigen-3.4.0' -B '$BUILD/eigen-build' \
			-DCMAKE_INSTALL_PREFIX='$PREFIX/eigen' -DBUILD_TESTING=OFF \
			&& cmake --install '$BUILD/eigen-build'"
	fi
	EIGEN3_DIR=$(dirname "$(find "$PREFIX/eigen" -name Eigen3Config.cmake | head -1)")
fi
if [ -z "${BOOST_ROOT:-}" ]; then
	if [ ! -d "$PREFIX/boost" ]; then
		tar -xJf "$SRC/boost-1.86.0-cmake.tar.xz" -C "$BUILD"
		step boost sh -c "cmake -S '$BUILD/boost-1.86.0' -B '$BUILD/boost-build' \
			-DBOOST_INCLUDE_LIBRARIES=math -DCMAKE_INSTALL_PREFIX='$PREFIX/boost' \
			-DBUILD_TESTING=OFF && cmake --build '$BUILD/boost-build' -j $JOBS \
			&& cmake --install '$BUILD/boost-build'"
	fi
	BOOST_ROOT="$PREFIX/boost"
fi

# --------------------------------------------------------------------------
# netCDF: MEQ's .nc output is not optional, and it is C++4.
# --------------------------------------------------------------------------
if [ -z "${NETCDF_DIR:-}" ] && command -v nc-config >/dev/null; then
	NETCDF_DIR=$(nc-config --prefix)
fi
[ -n "${NETCDF_DIR:-}" ] && [ -f "$NETCDF_DIR/include/netcdf.h" ] \
	|| die "no netCDF-C: module load cray-hdf5 cray-netcdf, or set NETCDF_DIR"
export PKG_CONFIG_PATH="$PREFIX/pkgconfig:$PREFIX/netcdf-cxx4/lib/pkgconfig:$PREFIX/netcdf-cxx4/lib64/pkgconfig:${PKG_CONFIG_PATH:-}"
mkdir -p "$PREFIX/pkgconfig"
if pkg-config --exists netcdf-cxx4 2>/dev/null; then
	say "netcdf-cxx4: $(pkg-config --modversion netcdf-cxx4) via pkg-config"
elif [ -f "$NETCDF_DIR/include/ncFile.h" ] && ls "$NETCDF_DIR"/lib*/libnetcdf_c++4.* >/dev/null 2>&1; then
	# The C++ bindings are in the netCDF install but nothing describes them.
	ncxx_lib=$(dirname "$(ls "$NETCDF_DIR"/lib*/libnetcdf_c++4.* | head -1)")
	cat > "$PREFIX/pkgconfig/netcdf-cxx4.pc" <<EOF
Name: netcdf-cxx4
Description: netCDF C++4 bindings found in $NETCDF_DIR
Version: 4.3.1
Cflags: -I$NETCDF_DIR/include
Libs: -L$ncxx_lib -Wl,-rpath,$ncxx_lib -lnetcdf_c++4 -lnetcdf
EOF
	say "netcdf-cxx4: the bindings in NETCDF_DIR, described by $PREFIX/pkgconfig/netcdf-cxx4.pc"
else
	if [ ! -d "$PREFIX/netcdf-cxx4" ]; then
		tar -xzf "$SRC/netcdf-cxx4-4.3.1.tar.gz" -C "$BUILD"
		nclib=$(dirname "$(ls "$NETCDF_DIR"/lib*/libnetcdf.so* | head -1)")
		step netcdf-cxx4 sh -c "cd '$BUILD/netcdf-cxx4-4.3.1' && ./configure \
			--prefix='$PREFIX/netcdf-cxx4' CPPFLAGS='-I$NETCDF_DIR/include' \
			LDFLAGS='-L$nclib -Wl,-rpath,$nclib' && make -j $JOBS && make install"
	fi
	pkg-config --exists netcdf-cxx4 || die "netcdf-cxx4 built but pkg-config cannot see it"
fi

# --------------------------------------------------------------------------
# SUNDIALS with CUDA: MFEM refuses a CUDA build against one without it.
# --------------------------------------------------------------------------
# gx's Makefile.perlmutter names the same thing SUNDIALS_ROOT.
SUNDIALS_DIR=${SUNDIALS_DIR:-${SUNDIALS_ROOT:-}}
if [ -n "${SUNDIALS_DIR:-}" ]; then
	# Anchored, as in preflight.sh: unanchored, both match a SUNDIALS without
	# CUDA in single precision.
	[ -f "$SUNDIALS_DIR/include/sundials/sundials_config.h" ] \
		|| die "SUNDIALS_DIR=$SUNDIALS_DIR has no include/sundials/sundials_config.h: it is the install prefix"
	grep -q '^#define SUNDIALS_CUDA_ENABLED' "$SUNDIALS_DIR/include/sundials/sundials_config.h" \
		|| die "SUNDIALS_DIR=$SUNDIALS_DIR is not a SUNDIALS built with ENABLE_CUDA=ON"
	grep -q '^#define SUNDIALS_DOUBLE_PRECISION 1' "$SUNDIALS_DIR/include/sundials/sundials_config.h" \
		|| die "SUNDIALS_DIR=$SUNDIALS_DIR is not double precision"
	say "SUNDIALS: $(sed -n 's/#define SUNDIALS_VERSION "\(.*\)"/\1/p' "$SUNDIALS_DIR/include/sundials/sundials_config.h") at $SUNDIALS_DIR"
else
	SUNDIALS_DIR="$PREFIX/sundials"
fi
if [ ! -d "$SUNDIALS_DIR" ]; then
	step sundials sh -c "cmake -S '$SRC/sundials' -B '$BUILD/sundials' \
		-DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX='$PREFIX/sundials' \
		-DENABLE_CUDA=ON -DCMAKE_CUDA_ARCHITECTURES=$CUDA_ARCH \
		-DCMAKE_CUDA_HOST_COMPILER=$CXX -DENABLE_OPENMP=ON \
		-DENABLE_LAPACK=ON -DLAPACK_LIBRARIES='$BLAS' -DBLAS_LIBRARIES='$BLAS' \
		-DSUNDIALS_INDEX_SIZE=64 -DSUNDIALS_PRECISION=DOUBLE \
		-DEXAMPLES_ENABLE_C=OFF -DEXAMPLES_ENABLE_CXX=OFF -DEXAMPLES_INSTALL=OFF \
		&& cmake --build '$BUILD/sundials' -j $JOBS && cmake --install '$BUILD/sundials'"
fi

# --------------------------------------------------------------------------
# MFEM: meq-integration, as CLAUDE_MFEM.md configures it, at sm_$CUDA_ARCH.
# --------------------------------------------------------------------------
# MFEM's FindSUNDIALS looks under SUNDIALS_DIR/lib only, so a SUNDIALS that
# installs to lib64/ (SUSE's default, Perlmutter's) is found through
# CMAKE_LIBRARY_PATH instead.
if [ ! -f "$PREFIX/mfem/include/mfem/mfem.hpp" ] && [ ! -f "$PREFIX/mfem/include/mfem.hpp" ]; then
	step mfem sh -c "cmake -S '$SRC/mfem' -B '$BUILD/mfem' \
		-DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX='$PREFIX/mfem' \
		-DCMAKE_CXX_COMPILER=$CXX -DMFEM_USE_CUDA=YES -DCUDA_ARCH=sm_$CUDA_ARCH \
		-DCMAKE_CUDA_ARCHITECTURES=$CUDA_ARCH -DCMAKE_CUDA_HOST_COMPILER=$CXX \
		-DMFEM_USE_CUDSS=YES -DCUDSS_DIR='$CUDSS_DIR' \
		-DMFEM_USE_SUNDIALS=YES -DSUNDIALS_DIR='$SUNDIALS_DIR' \
		-DCMAKE_LIBRARY_PATH='$SUNDIALS_DIR/lib64;$SUNDIALS_DIR/lib' \
		-DMFEM_USE_MKL_PARDISO=YES -DMKL_PARDISO_DIR='$MKLROOT' -DMKL_LIBRARY_DIR=$MKL_LIBDIR \
		-DMKL_PARDISO_MKL_SEQUENTIAL_LIBRARY='$MKL_LIB/libmkl_gnu_thread.so' \
		-DMFEM_USE_LAPACK=YES -DBLAS_LIBRARIES='$BLAS' -DLAPACK_LIBRARIES='$BLAS' \
		-DMFEM_USE_OPENMP=YES -DMFEM_THREAD_SAFE=YES -DMFEM_USE_EXCEPTIONS=YES \
		-DMFEM_USE_MEMALLOC=ON -DMFEM_USE_ZLIB=YES -DMFEM_USE_MPI=NO \
		-DMFEM_USE_SUITESPARSE=NO -DMFEM_USE_GSLIB=NO -DMFEM_PRECISION=double \
		-DMFEM_ENABLE_TESTING=OFF -DMFEM_ENABLE_EXAMPLES=OFF -DMFEM_ENABLE_MINIAPPS=OFF \
		&& cmake --build '$BUILD/mfem' -j $JOBS && cmake --install '$BUILD/mfem'"
fi

# --------------------------------------------------------------------------
# MEQ. Tests are off: the correctness gate is check.sh, through the driver.
# --------------------------------------------------------------------------
step meq sh -c "cmake -S '$SRC/meq' -B '$BUILD/meq' -DCMAKE_BUILD_TYPE=Release \
	-DCMAKE_INSTALL_PREFIX='$PREFIX/meq' -DCMAKE_CXX_COMPILER=$CXX \
	-DCMAKE_CUDA_ARCHITECTURES=$CUDA_ARCH -DCMAKE_CUDA_HOST_COMPILER=$CXX \
	-DMFEM_DIR='$PREFIX/mfem' -DTOML11_DIR='$SRC/meq/extern/toml11' \
	-DMEQ_BUILD_TESTS=OFF -DMEQ_RUNS_DIR='$KIT/results' \
	-DBoost_ROOT='$BOOST_ROOT' -DEigen3_DIR='$EIGEN3_DIR' \
	&& cmake --build '$BUILD/meq' -j $JOBS && cmake --install '$BUILD/meq'"

# The run scripts source this; it is what the libraries above need at run time.
cat > "$KIT/env.sh" <<EOF
# Written by build.sh, $(date -Iseconds). Sourced by check.sh and run.sh.
export MEQ="$PREFIX/meq/bin/meq"
export LD_LIBRARY_PATH="$CUDSS_LIB:$MKL_LIB:$SUNDIALS_DIR/lib:$SUNDIALS_DIR/lib64:$PREFIX/netcdf-cxx4/lib:$PREFIX/netcdf-cxx4/lib64:\${LD_LIBRARY_PATH:-}"
EOF
# shellcheck disable=SC1091
. "$KIT/env.sh"
grep -q 'cuda' <("$MEQ" --version) || die "meq --version does not report a CUDA build"
say "done: $("$MEQ" --version)"
say "next: ./check.sh"
