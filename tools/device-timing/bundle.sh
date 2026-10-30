#!/usr/bin/env bash
# Pack the device timing kit into one tarball for a machine with no Claude Code
# and possibly no internet: MEQ's source as it stands, MFEM meq-integration,
# SUNDIALS, the header libraries, netcdf-cxx4, the case inputs and these
# scripts. Run on the development machine:
#
#   tools/device-timing/bundle.sh [OUTDIR]       default: <project>/runs/device-timing
#
# Inputs, overridable:
#   MFEM_SRC, MFEM_REF          ../../mfem/mfem-src, meq-integration
#   SUNDIALS_SRC, SUNDIALS_REF  ../../sundials/src, v7.5.0 (what MFEM's CUDA build links)
#   DOWNLOADS                   <project>/runs/device-timing/downloads, holding
#                               eigen/eigen-3.4.0.tar.gz,
#                               netcdf-cxx4/netcdf-cxx4-4.3.1.tar.gz and
#                               boost/boost-1.86.0-cmake.tar.xz
#   MESH                        a machine-f-diiid.msh made by meq-run --mesh-only

set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
SOURCE=$(cd "$HERE/../.." && pwd)
PROJECT=$(cd "$SOURCE/.." && pwd)
OUT=${1:-$PROJECT/runs/device-timing}
MFEM_SRC=${MFEM_SRC:-$PROJECT/../mfem/mfem-src}
MFEM_REF=${MFEM_REF:-meq-integration}
SUNDIALS_SRC=${SUNDIALS_SRC:-$PROJECT/../sundials/src}
SUNDIALS_REF=${SUNDIALS_REF:-v7.5.0}
DOWNLOADS=${DOWNLOADS:-$PROJECT/runs/device-timing/downloads}
MESH=${MESH:-$PROJECT/runs/device-timing/mesh/machine-f-diiid.msh}

die() { echo "bundle: $*" >&2; exit 1; }
for f in eigen/eigen-3.4.0.tar.gz netcdf-cxx4/netcdf-cxx4-4.3.1.tar.gz boost/boost-1.86.0-cmake.tar.xz; do
	[ -f "$DOWNLOADS/$f" ] || die "missing $DOWNLOADS/$f (see this script's header)"
done
if [ ! -f "$MESH" ]; then
	# The mesh is a build product of [mesh.generate]; make it once, here,
	# where gmsh is, so the GPU machine needs none.
	mkdir -p "$(dirname "$MESH")"
	( cd "$(dirname "$MESH")" && "$PROJECT/install/double/bin/meq-run" --mesh-only \
		"$SOURCE/examples/machine-f-diiid.toml" ) || die "cannot make the machine-f-diiid mesh"
fi

STAGE=$(mktemp -d)
trap 'rm -rf "$STAGE"' EXIT
KIT="$STAGE/meq-device-timing"
mkdir -p "$KIT/src/meq" "$KIT/cases"

# The scripts, at the kit's top level.
cp "$HERE"/{README.md,preflight.sh,build.sh,check.sh,cases.sh,run.sh,job.slurm,summarize.py,update.sh} "$KIT/"

# MEQ as it stands: tracked files from the working tree, toml11 with them, and
# new files not yet added (not ignored ones), which the diff cannot carry.
( cd "$SOURCE" && git ls-files -z --recurse-submodules | tar --null -T - -cf - ) | tar -xf - -C "$KIT/src/meq"
untracked=$(git -C "$SOURCE" ls-files --others --exclude-standard)
[ -z "$untracked" ] || ( cd "$SOURCE" && git ls-files -z --others --exclude-standard | tar --null -T - -cf - ) | tar -xf - -C "$KIT/src/meq"
git -C "$SOURCE" diff HEAD > "$KIT/src/meq.diff"
# What meq --version and every .nc will say, in cmake/GitCommit.cmake's spelling.
commit=$(git -C "$SOURCE" rev-parse --short HEAD)
[ -n "$(git -C "$SOURCE" status --porcelain --untracked-files=normal)" ] && commit="$commit-dirty"
echo "$commit" > "$KIT/src/meq/GIT_COMMIT"

# MFEM and SUNDIALS from their repositories, at the named refs.
git -C "$MFEM_SRC" archive --prefix=mfem/ "$MFEM_REF" | tar -xf - -C "$KIT/src"
git -C "$SUNDIALS_SRC" archive --prefix=sundials/ "$SUNDIALS_REF" | tar -xf - -C "$KIT/src"
cp "$DOWNLOADS/eigen/eigen-3.4.0.tar.gz" "$DOWNLOADS/netcdf-cxx4/netcdf-cxx4-4.3.1.tar.gz" \
	"$DOWNLOADS/boost/boost-1.86.0-cmake.tar.xz" "$KIT/src/"

# The cases: every file the two configurations read, and the mesh.
for f in machine-f-diiid.toml machine-f-diiid-pprime.dat machine-f-diiid-ggprime.dat \
         machine-f-diiid-guess.gf machine-f-diiid-guess.mesh \
         fixed-f-diiid.toml fixed-f-diiid-pprime.dat fixed-f-diiid-ggprime.dat; do
	cp "$SOURCE/examples/$f" "$KIT/cases/"
done
cp "$MESH" "$MESH.meq-mesh" "$KIT/cases/"

{
	echo "MEQ device timing kit, bundled $(date -Iseconds) on $(hostname)"
	echo "meq       $(git -C "$SOURCE" rev-parse HEAD) $( [ -s "$KIT/src/meq.diff" ] && echo "+ src/meq.diff ($(wc -l < "$KIT/src/meq.diff") lines, uncommitted)")"
	[ -z "$untracked" ] || echo "          + new, not in the diff: $(echo $untracked)"
	echo "mfem      $MFEM_REF at $(git -C "$MFEM_SRC" rev-parse "$MFEM_REF") (github.com/ianabel/mfem-meq-integration)"
	echo "sundials  $SUNDIALS_REF at $(git -C "$SUNDIALS_SRC" rev-parse "$SUNDIALS_REF^{commit}")"
	echo "eigen 3.4.0, netcdf-cxx4 4.3.1, boost 1.86.0 (cmake), sha256:"
	( cd "$KIT/src" && sha256sum eigen-3.4.0.tar.gz netcdf-cxx4-4.3.1.tar.gz boost-1.86.0-cmake.tar.xz | sed 's/^/  /' )
} > "$KIT/PROVENANCE.txt"

chmod +x "$KIT"/*.sh "$KIT/summarize.py"
mkdir -p "$OUT"
name="meq-device-timing-$(date +%Y%m%d-%H%M).tar.gz"
tar -czf "$OUT/$name" -C "$STAGE" meq-device-timing
echo "bundle: $OUT/$name ($(du -h "$OUT/$name" | cut -f1))"
cat "$KIT/PROVENANCE.txt"
