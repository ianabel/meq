#!/usr/bin/env bash
# Bring an unpacked kit up to date from git instead of a new bundle, then
# build. Run on a login node (it needs github.com), in place of build.sh:
#
#   ./update.sh                     MEQ main and MFEM main, the newest of each
#   MEQ_GIT_REF=607e43b ./update.sh any branch, tag or full commit of either
#
#   MEQ_GIT_URL, MEQ_GIT_REF    https://github.com/ianabel/meq.git, main
#   MFEM_GIT_URL, MFEM_GIT_REF  https://github.com/ianabel/mfem-meq-integration.git, main
#
# The first run replaces the bundle's src/meq and src/mfem with clones (the
# bundled trees are kept as src/*.bundled). After that, each run fetches and
# checks out the refs. MFEM is rebuilt only when its source moved: build.sh
# skips any library already installed, so this removes MFEM's install when the
# commit it was built from is not the one checked out. MEQ always rebuilds,
# incrementally. SUNDIALS, netcdf-cxx4, Eigen and Boost are left as they are.
#
# The kit's scripts and case inputs are refreshed from MEQ's
# tools/device-timing and examples/, except any you have edited (job.slurm,
# with your account in it): those are kept, and named.
#
# The body is one function, read whole before it runs, so that refreshing this
# script mid-run cannot change what is running.

main()
{
	set -euo pipefail
	KIT=$(cd "$(dirname "$0")" && pwd)
	SRC="$KIT/src"
	PREFIX=${MEQ_PREFIX:-$KIT/install}
	MEQ_GIT_URL=${MEQ_GIT_URL:-https://github.com/ianabel/meq.git}
	MEQ_GIT_REF=${MEQ_GIT_REF:-main}
	MFEM_GIT_URL=${MFEM_GIT_URL:-https://github.com/ianabel/mfem-meq-integration.git}
	MFEM_GIT_REF=${MFEM_GIT_REF:-main}
	SCRIPTS=(README.md preflight.sh build.sh check.sh cases.sh run.sh job.slurm summarize.py update.sh)
	CASES=(machine-f-diiid.toml machine-f-diiid-pprime.dat machine-f-diiid-ggprime.dat
	       machine-f-diiid-guess.gf machine-f-diiid-guess.mesh
	       fixed-f-diiid.toml fixed-f-diiid-pprime.dat fixed-f-diiid-ggprime.dat)

	say() { echo "[update] $*"; }
	die() { echo "[update] FAILED: $*" >&2; exit 1; }

	# What the kit's scripts were shipped as, before anything moves: a script
	# that differs from it was edited here, and is kept.
	OLD=$(mktemp -d)
	trap 'rm -rf "$OLD"' EXIT
	for f in "${SCRIPTS[@]}"; do
		[ -f "$SRC/meq/tools/device-timing/$f" ] && cp "$SRC/meq/tools/device-timing/$f" "$OLD/$f"
	done

	checkout "$SRC/meq" "$MEQ_GIT_URL" "$MEQ_GIT_REF"
	# toml11 only; its own submodules are its docs and tests.
	git -C "$SRC/meq" submodule -q update --init
	checkout "$SRC/mfem" "$MFEM_GIT_URL" "$MFEM_GIT_REF"

	local mfem_commit built=""
	mfem_commit=$(git -C "$SRC/mfem" rev-parse HEAD)
	# The first clone, of the same MFEM the bundle carried (the mirror adds
	# PROVENANCE.md and leaves out .github/): the install stands.
	if [ -d "$PREFIX/mfem" ] && [ ! -f "$PREFIX/mfem/KIT_SOURCE_COMMIT" ] && [ -d "$SRC/mfem.bundled" ] \
		&& diff -rq --exclude=.git --exclude=.github --exclude=PROVENANCE.md "$SRC/mfem.bundled" "$SRC/mfem" > /dev/null; then
		echo "$mfem_commit" > "$PREFIX/mfem/KIT_SOURCE_COMMIT"
	fi
	[ -f "$PREFIX/mfem/KIT_SOURCE_COMMIT" ] && built=$(cat "$PREFIX/mfem/KIT_SOURCE_COMMIT")
	if [ -d "$PREFIX/mfem" ] && [ "$built" != "$mfem_commit" ]; then
		say "MFEM's install is from ${built:-the bundle}, not ${mfem_commit:0:10}: rebuilding MFEM"
		rm -rf "$PREFIX/mfem"
	fi

	local kept=()
	for f in "${SCRIPTS[@]}"; do
		local new="$SRC/meq/tools/device-timing/$f"
		[ -f "$new" ] || continue
		if [ -f "$KIT/$f" ] && [ -f "$OLD/$f" ] && ! cmp -s "$KIT/$f" "$OLD/$f"; then
			cmp -s "$OLD/$f" "$new" || kept+=("$f")
			continue
		fi
		cp "$new" "$KIT/$f"
	done
	chmod +x "$KIT"/*.sh "$KIT/summarize.py"
	for f in "${CASES[@]}"; do
		cp "$SRC/meq/examples/$f" "$KIT/cases/$f"
	done
	# cases.sh writes the refinement levels' configurations from
	# fixed-f-diiid.toml when they are missing; they go so that it does again.
	rm -f "$KIT"/cases/fixed-f-diiid-L*.toml

	say "meq  $MEQ_GIT_REF at $(git -C "$SRC/meq" rev-parse --short HEAD)"
	say "mfem $MFEM_GIT_REF at ${mfem_commit:0:10} ($(git -C "$SRC/mfem" log -1 --format=%s))"
	"$KIT/build.sh"
	echo "$mfem_commit" > "$PREFIX/mfem/KIT_SOURCE_COMMIT"

	for f in "${kept[@]}"; do
		say "KEPT your edited $f; the new one is src/meq/tools/device-timing/$f"
	done
	# The machine-f mesh in cases/ is gmsh's, made on the development machine,
	# and mf runs it with --mesh-ready: a claim MEQ takes on trust. Its
	# sidecar holds the generator command it was made with, which is what
	# --mesh-command prints for a configuration whose mesh is current.
	if ! diff -q <("$PREFIX/meq/bin/meq" --mesh-command "$KIT/cases/machine-f-diiid.toml") \
		"$KIT/cases/machine-f-diiid.msh.meq-mesh" > /dev/null; then
		die "machine-f-diiid.toml's mesh changed, so cases/machine-f-diiid.msh is stale: mf needs a new bundle's mesh (or gmsh here)"
	fi
	say "done; next: ./check.sh"
}

checkout() # directory, url, ref: a clone at ref, made from the bundle's tree if it has none
{
	local dir=$1 url=$2 ref=$3
	if [ ! -d "$dir/.git" ]; then
		say "cloning $url into $dir"
		rm -rf "$dir.new"
		git clone -q "$url" "$dir.new" || die "cannot clone $url (a private repository needs your GitHub credentials here)"
		if [ -d "$dir" ]; then
			rm -rf "$dir.bundled"
			mv "$dir" "$dir.bundled"
		fi
		mv "$dir.new" "$dir"
	fi
	[ -z "$(git -C "$dir" status --porcelain --untracked-files=no)" ] \
		|| die "$dir has local changes; commit or discard them first"
	git -C "$dir" fetch -q origin "$ref" || die "no $ref at $url"
	git -C "$dir" checkout -q --detach FETCH_HEAD
}

main "$@"
exit
