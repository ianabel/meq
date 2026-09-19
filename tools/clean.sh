#!/usr/bin/env bash
#
# Remove the temporary files a MEQ checkout accumulates.
#
# WHAT MAKES THIS SAFE IS NOT THE PATTERN LIST, IT IS `git check-ignore`.
# Every candidate is confirmed IGNORED before it is removed, so a file git
# tracks -- or one that is merely untracked, which is how work in progress and
# incoming correspondence look -- is skipped whatever it is named. That is why
# `examples/limited-tokamak-guess.{mesh,gf}` survive a rule that deletes
# `*.mesh` and `*.gf`: .gitignore negates them, check-ignore says so, and this
# script never has to know they are special.
#
# A directory is removed only when git tracks nothing under it AND nothing
# under it is untracked-but-not-ignored. One stray `??` anywhere inside and the
# whole directory is left alone.
#
# THREE THINGS ARE NEVER REMOVED even though they are ignored, because they are
# expensive or impossible to recreate rather than merely large:
#
#   refs/*.pdf                     fetched by doi one at a time; refs/Refs.md
#                                  is the index and the PDFs are not in git
#   .venv-docs/                    the pinned sphinx 7.4.7 / sphinx-material
#                                  0.0.36 the docs build needs; the system
#                                  sphinx is 8.2.3 and above the pin
#   */venv/, */ref-n*/             the freegs4e benchmark's interpreter and its
#                                  reference solutions
#
# Usage:
#   tools/clean.sh [options]
#
#   (no options)     run artefacts and editor scratch -- the safe default
#   --builds         also every build-*/ directory EXCEPT ./build
#   --all-builds     also ./build itself
#   --docs           also docs/_build, docs/doctrees and the LaTeX aux files
#   --all            everything above except ./build
#   -n, --dry-run    say what would go, remove nothing
#   -q, --quiet      totals only
#   -h, --help       this
#
set -uo pipefail

dryRun=0
quiet=0
doBuilds=0
doLiveBuild=0
doDocs=0

while [ $# -gt 0 ]; do
	case "$1" in
		--builds)      doBuilds=1 ;;
		--all-builds)  doBuilds=1; doLiveBuild=1 ;;
		--docs)        doDocs=1 ;;
		--all)         doBuilds=1; doDocs=1 ;;
		-n|--dry-run)  dryRun=1 ;;
		-q|--quiet)    quiet=1 ;;
		-h|--help)     sed -n '3,40p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
		*)             printf 'clean.sh: unknown option %s (try --help)\n' "$1" >&2; exit 2 ;;
	esac
	shift
done

root=$( git rev-parse --show-toplevel 2>/dev/null )
if [ -z "$root" ]; then
	printf 'clean.sh: not inside a git checkout\n' >&2
	exit 1
fi
# Refuse anywhere that is not MEQ. A cleaner that runs in the wrong tree is the
# one kind of bug this script must not have, and two files settle it.
if [ ! -d "$root/src/meq" ] || [ ! -f "$root/CLAUDE.md" ]; then
	printf 'clean.sh: %s does not look like a MEQ checkout; refusing\n' "$root" >&2
	exit 1
fi
cd "$root" || exit 1

reclaimed=0
removedCount=0
keptCount=0

# Expensive to recreate, so kept whatever the patterns say.
protected()
{
	case "$1" in
		refs/*.pdf|.venv-docs|.venv-docs/*) return 0 ;;
		*/venv|*/venv/*|*/.venv|*/.venv/*) return 0 ;;
		tools/freegs4e-benchmark/ref-*)    return 0 ;;
	esac
	return 1
}

# A file goes only if git says it is ignored.
fileIsRemovable()
{
	[ -f "$1" ] || [ -L "$1" ] || return 1
	protected "$1" && return 1
	git check-ignore -q -- "$1"
}

# A directory goes only if git tracks nothing inside it and nothing inside it
# is untracked-but-not-ignored.
dirIsRemovable()
{
	[ -d "$1" ] || return 1
	protected "$1" && return 1
	[ -n "$( git ls-files -- "$1" | head -1 )" ] && return 1
	[ -n "$( git status --porcelain --ignored=matching -uall -- "$1" 2>/dev/null \
	         | grep -m1 '^??' )" ] && return 1
	return 0
}

remove()
{
	local path=$1 why=$2 size
	size=$( du -sb "$path" 2>/dev/null | cut -f1 )
	[ -z "$size" ] && size=0
	reclaimed=$(( reclaimed + size ))
	removedCount=$(( removedCount + 1 ))
	[ "$quiet" -eq 1 ] || printf '  %-10s %s  (%s)\n' \
		"$( numfmt --to=iec --suffix=B "$size" 2>/dev/null || echo "$size" )" \
		"$path" "$why"
	[ "$dryRun" -eq 1 ] || rm -rf -- "$path"
}

seen=" "

consider()
{
	local path=$1 why=$2
	# A path can be reached twice -- <stem>_cycles/ matches both the ParaView
	# sweep and the explicit list below it -- and in a dry run, where nothing
	# is actually removed, that would count its bytes twice and report a
	# total nobody can reconcile against `du`.
	case "$seen" in *" $path "*) return ;; esac
	seen="$seen$path "
	if [ -d "$path" ] && [ ! -L "$path" ]; then
		if dirIsRemovable "$path"; then remove "$path" "$why"
		else keptCount=$(( keptCount + 1 )); fi
	elif fileIsRemovable "$path"; then
		remove "$path" "$why"
	else
		keptCount=$(( keptCount + 1 ))
	fi
}

heading()
{
	[ "$quiet" -eq 1 ] || printf '\n%s\n' "$1"
}

[ "$dryRun" -eq 1 ] && printf 'DRY RUN -- nothing will be removed\n'

# ---------------------------------------------------------------------------
# 1. Solve outputs.
#
# apps/meq writes the same equilibrium three times and optionally a fourth
# reduction, so one run leaves a .mesh, three .gf, a .nc, a ParaView directory
# and sometimes a _surfaces.nc. They land in whatever directory the run was
# started from, which is usually the repo root.
# ---------------------------------------------------------------------------
heading 'Solve outputs'
while IFS= read -r -d '' f; do
	consider "$f" 'solve output'
done < <( find . -maxdepth 2 \
	-path ./.git -prune -o \
	-path './build*' -prune -o \
	-path ./extern -prune -o \
	-path ./.venv-docs -prune -o \
	-path './tools/freegs4e-benchmark/venv' -prune -o \
	-path './tools/freegs4e-benchmark/ref-*' -prune -o \
	-type f \( -name '*.mesh' -o -name '*.gf' -o -name '*.nc' \
	        -o -name '*.vtk' -o -name '*.vtu' -o -name '*.pvtu' -o -name '*.pvd' \
	        -o -name '*.meq-mesh' -o -name '*.msh' \
	        -o -name 'driver-acceptance-*.toml' \) -print0 )

# The ParaView collections are directories -- <stem>/ holding <stem>.pvd and
# Cycle??????/, and <stem>_cycles/ for an adaptive run. Identified by what is
# inside them rather than by name, so a directory that merely shares a stem
# with an example is not at risk.
heading 'ParaView collections'
for d in */ ; do
	d=${d%/}
	case "$d" in
		.git|build|build-*|extern|apps|src|tests|tools|docs|examples|refs|cmake|.venv-docs|.claude) continue ;;
	esac
	if compgen -G "$d/*.pvd" > /dev/null || compgen -G "$d/Cycle??????" > /dev/null; then
		consider "$d" 'ParaView collection'
	fi
done
for d in *_cycles visit_dump* ParaView ; do
	[ -e "$d" ] && consider "$d" 'ParaView collection'
done

# ---------------------------------------------------------------------------
# 2. Editor, tooling and WSL scratch.
# ---------------------------------------------------------------------------
heading 'Editor and tooling scratch'
for p in .cache .clangd compile_commands.json tags tmp gmon.out core ; do
	[ -e "$p" ] && consider "$p" 'scratch'
done
while IFS= read -r -d '' f; do
	consider "$f" 'scratch'
done < <( find . -path ./.git -prune -o -path './build*' -prune -o \
	-path ./extern -prune -o -path ./.venv-docs -prune -o \
	-type f \( -name '*~' -o -name '*.swp' -o -name '*.swo' \
	        -o -name '*:Zone.Identifier' -o -name 'core.[0-9]*' \) -print0 )

# ---------------------------------------------------------------------------
# 3. Documentation build products -- only with --docs, because a docs build is
#    minutes rather than seconds and is often the thing being iterated on.
# ---------------------------------------------------------------------------
if [ "$doDocs" -eq 1 ]; then
	heading 'Documentation build products'
	for p in docs/_build docs/doctrees ; do
		[ -e "$p" ] && consider "$p" 'sphinx output'
	done
	while IFS= read -r -d '' f; do
		consider "$f" 'latex aux'
	done < <( find docs/manual -maxdepth 1 -type f \
		\( -name '*.aux' -o -name '*.bbl' -o -name '*.blg' -o -name '*.fdb_latexmk' \
		-o -name '*.fls' -o -name '*.idx' -o -name '*.ind' -o -name '*.log' \
		-o -name '*.nav' -o -name '*.out' -o -name '*.snm' -o -name '*.toc' \
		-o -name '*.synctex.gz' \) -print0 2>/dev/null )
fi

# ---------------------------------------------------------------------------
# 4. Build directories -- only when asked, and ./build only when asked twice.
#
# ./build is the one every command in CLAUDE.md names, so losing it costs a
# full rebuild on the next thing anybody does. The build-*/ siblings are
# experiments and are what actually fills the disk.
# ---------------------------------------------------------------------------
if [ "$doBuilds" -eq 1 ]; then
	heading 'Build directories'
	for d in build-* ; do
		[ -d "$d" ] || continue
		consider "$d" 'build directory'
	done
	if [ "$doLiveBuild" -eq 1 ] && [ -d build ]; then
		consider build 'build directory (the live one)'
	elif [ -d build ]; then
		[ "$quiet" -eq 1 ] || printf '  %-10s %s  (kept -- pass --all-builds)\n' \
			"$( du -sh build 2>/dev/null | cut -f1 )" build
	fi
fi

# ---------------------------------------------------------------------------
printf '\n%s %d path(s), %s%s\n' \
	"$( [ "$dryRun" -eq 1 ] && echo 'Would remove' || echo 'Removed' )" \
	"$removedCount" \
	"$( numfmt --to=iec --suffix=B "$reclaimed" 2>/dev/null || echo "$reclaimed bytes" )" \
	"$( [ "$keptCount" -gt 0 ] && printf '; %d candidate(s) kept because git tracks them or they are not ignored' "$keptCount" )"

if [ "$doBuilds" -eq 0 ]; then
	printf 'Build directories were not touched. tools/clean.sh --builds removes the build-*/ experiments.\n'
fi
