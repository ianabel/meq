#!/usr/bin/env bash
# Refresh github.com/ianabel/mfem-meq-integration: meq-integration's source as
# a flattened snapshot, so a GPU machine can `git pull` MFEM (update.sh)
# instead of unpacking a new bundle. Run on the development machine:
#
#   tools/device-timing/mirror.sh
#
# Each refresh is ONE new commit on the mirror's main, holding the whole tree
# at MFEM_REF; none of MFEM's history goes out. A commit on top, never a
# rewrite, so a clone fast-forwards. meq-integration is a local merge made only
# to be built against, so this is the only place it is published: it never goes
# to MFEM.
#
# Inputs, overridable:
#   MFEM_SRC, MFEM_REF  ../../mfem/mfem-src, meq-integration
#   MFEM_BASE           the MFEM master meq-integration is merged on (45799ebea8)
#   UPSTREAM            the remote holding the HDG branches (hdgdev)
#   MIRROR_URL          github:ianabel/mfem-meq-integration.git
#   MIRROR_CACHE        <project>/runs/device-timing/mirror, its working clone
#   MIRROR_PUSH         yes; "no" commits in the clone and stops, to look first

set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
PROJECT=$(cd "$HERE/../../.." && pwd)
MFEM_SRC=${MFEM_SRC:-$PROJECT/../mfem/mfem-src}
MFEM_REF=${MFEM_REF:-meq-integration}
MFEM_BASE=${MFEM_BASE:-45799ebea8}
UPSTREAM=${UPSTREAM:-hdgdev}
MIRROR_URL=${MIRROR_URL:-github:ianabel/mfem-meq-integration.git}
CACHE=${MIRROR_CACHE:-$PROJECT/runs/device-timing/mirror}

die() { echo "mirror: $*" >&2; exit 1; }
ref=$(git -C "$MFEM_SRC" rev-parse --verify "$MFEM_REF^{commit}") || die "no $MFEM_REF in $MFEM_SRC"
short=${ref:0:10}
date=$(git -C "$MFEM_SRC" log -1 --format=%cs "$ref")
base=$(git -C "$MFEM_SRC" rev-parse --short=10 "$MFEM_BASE")
git -C "$MFEM_SRC" merge-base --is-ancestor "$MFEM_BASE" "$ref" || die "$MFEM_BASE is not under $MFEM_REF: set MFEM_BASE"

if [ -d "$CACHE/.git" ]; then
	git -C "$CACHE" fetch -q origin
	git -C "$CACHE" checkout -q main
	git -C "$CACHE" reset -q --hard origin/main
else
	git clone -q "$MIRROR_URL" "$CACHE"
fi

# The tree, replaced whole: a file meq-integration dropped must go too.
# .github/ stays out: it is MFEM's CI, which GitHub would run on every push.
git -C "$CACHE" rm -rq --ignore-unmatch .
git -C "$MFEM_SRC" archive "$ref" | tar -xf - -C "$CACHE" --exclude=.github

# Which of the HDG line's branches it holds, at which tip. Every one listed is
# wholly in; one only partly in is named as such, because a fix on its tip
# would then be missing.
branches=$(
	for b in $(git -C "$MFEM_SRC" for-each-ref --format='%(refname:short)' "refs/remotes/$UPSTREAM" \
	           | grep -v -E "/backup/|/wip/|/HEAD$|^$UPSTREAM$"); do
		tip=$(git -C "$MFEM_SRC" rev-parse "$b")
		git -C "$MFEM_SRC" merge-base --is-ancestor "$tip" "$MFEM_BASE" && continue
		if git -C "$MFEM_SRC" merge-base --is-ancestor "$tip" "$ref"; then
			printf '    %-38s %s\n' "${b#"$UPSTREAM"/}" "${tip:0:10}"
		elif ! git -C "$MFEM_SRC" merge-base --is-ancestor "$(git -C "$MFEM_SRC" merge-base "$tip" "$ref")" "$MFEM_BASE"; then
			printf '    %-38s %s, PARTLY: %s commits on its tip are not in\n' "${b#"$UPSTREAM"/}" "${tip:0:10}" \
				"$(git -C "$MFEM_SRC" rev-list --count "$ref..$tip")"
		fi
	done
)

cat > "$CACHE/PROVENANCE.md" <<EOF
# mfem-meq-integration — a flattened snapshot, not a fork

This tree is **MFEM 4.10.1 as MEQ builds against it**: the source of branch
\`meq-integration\` at \`$ref\`
($date), exported with \`git archive\`. Each commit here is one such snapshot,
whole, with no MFEM history; a clone fast-forwards from one to the next.
\`.github/\` is left out, so that GitHub does not run MFEM's CI here.

\`meq-integration\` is a local merge, made only to be built against, of MFEM
\`master\` at \`$base\` and the following branches of the HDG development line
(\`mfem-hdg-dev\`), each at the tip it is in at:

\`\`\`
$branches
\`\`\`

**Do not send this upstream.** It is a custom blend; each branch goes upstream on
its own. When they are all merged into MFEM, this repository is retired.

Built for MEQ with CMake, at minimum:

    -DMFEM_USE_CUDA=YES -DCUDA_ARCH=sm_80 (A100; sm_75 on the RTX 2070 SUPER)
    -DMFEM_USE_CUDSS=YES -DMFEM_USE_SUNDIALS=YES (SUNDIALS built with CUDA)
    -DMFEM_USE_MKL_PARDISO=YES and/or -DMFEM_USE_SUITESPARSE=YES
    -DMFEM_USE_OPENMP=YES -DMFEM_THREAD_SAFE=YES -DMFEM_USE_LAPACK=YES
    -DMFEM_USE_EXCEPTIONS=YES -DMFEM_USE_MPI=NO

MEQ's \`CLAUDE_MFEM.md\` is the full record, and MEQ's
\`tools/device-timing/build.sh\` the recipe used off the development machine.
EOF

git -C "$CACHE" add -A
if git -C "$CACHE" diff --cached --quiet; then
	echo "mirror: already at $MFEM_REF $short"
	exit 0
fi
git -C "$CACHE" commit -q -m "meq-integration at $short, flattened

MFEM 4.10.1 master at $base with the HDG development branches MEQ
builds against merged in; see PROVENANCE.md. The whole tree, with no
MFEM history.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
if [ "${MIRROR_PUSH:-yes}" = no ]; then
	echo "mirror: committed in $CACHE, not pushed (MIRROR_PUSH=no)"
	exit 0
fi
git -C "$CACHE" push -q origin main
echo "mirror: $MIRROR_URL main is now meq-integration $short ($(git -C "$CACHE" rev-parse --short HEAD))"
