# The timing cases, shared by check.sh and run.sh. Sourced, not run.
#
#   mf    examples/machine-f-diiid.toml: the bordered free-boundary machine,
#         4,848 elements; its gmsh mesh is shipped in cases/, so no gmsh
#   ffL   examples/fixed-f-diiid.toml at RefinementLevels = L (1, 2, 3):
#         6,245, 25,487 and 102,840 elements; it meshes itself
#
# MEASUREMENTS.md M-235 and M-242 are the same matrix on an RTX 2070 SUPER.

CASEDIR="$KIT/cases"

case_toml() # case -> the configuration to run
{
	case "$1" in
		mf) echo "$CASEDIR/machine-f-diiid.toml" ;;
		ff[1-9]) echo "$CASEDIR/fixed-f-diiid-L${1#ff}.toml" ;;
		*) echo "unknown case $1" >&2; return 1 ;;
	esac
}

case_flags() # case -> extra driver flags
{
	case "$1" in
		mf) echo "--mesh-ready" ;;
		*) echo "" ;;
	esac
}

prepare_case() # case, run directory: make the directory a run can start in
{
	local c=$1 dir=$2
	rm -rf "$dir"
	mkdir -p "$dir"
	case "$c" in
		# A generated mesh is read from the WORKING directory (CLAUDE.md,
		# "Where things go"), so it goes where the run is.
		mf) cp "$CASEDIR/machine-f-diiid.msh" "$CASEDIR/machine-f-diiid.msh.meq-mesh" "$dir/" ;;
		ff[1-9])
			local toml
			toml=$(case_toml "$c")
			[ -f "$toml" ] || sed -E "s/^RefinementLevels = .*/RefinementLevels = ${c#ff}/" \
				"$CASEDIR/fixed-f-diiid.toml" > "$toml"
			;;
	esac
}
