"""Where MEQ's tools find the solver, its sources and their own output.

THE PROJECT DIRECTORY HOLDS FOUR THINGS, AND ONLY ONE OF THEM IS TRACKED:

    <project>/source/            this repository -- nothing is written here
    <project>/builds/<name>/     one CMake build tree per configuration
    <project>/install/<name>/    `cmake --install builds/<name>`: bin/meq, bin/meq-run
    <project>/runs/              everything a run writes, a tool's under runs/<tool>/

builds/ and runs/ can be deleted whenever their contents are no longer wanted,
so a tool writes into runs/ and runs the INSTALLED binary, never one in a build
tree.  The sibling checkouts this project benchmarks against -- ../freegs4e,
../geq, ../mfem -- sit beside the project directory, not beside source/.

Each location reads an environment variable first:

    MEQ_RUNS_DIR       runs/
    MEQ_INSTALL        an install prefix; default install/$MEQ_BUILD
    MEQ_BUILD          which install/<name>; default "double"
    MEQ_SIBLINGS       the directory holding ../freegs4e and friends

Import it from anywhere under tools/ with

    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
    import meqpaths

adjusting the ".." to the depth of the script.
"""

import os

SOURCE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PROJECT = os.path.dirname(SOURCE)
SIBLINGS = os.environ.get("MEQ_SIBLINGS", os.path.dirname(PROJECT))
RUNS = os.environ.get("MEQ_RUNS_DIR", os.path.join(PROJECT, "runs"))


def source(*parts):
    """A path in the source tree: source("examples", "mhd-rectangle.toml")."""
    return os.path.join(SOURCE, *parts)


def sibling(name, *parts):
    """A checkout beside the project directory: sibling("geq")."""
    return os.path.join(SIBLINGS, name, *parts)


def runs(*parts, create=True):
    """A directory under runs/, made if it is not there: runs("geq-benchmark")."""
    path = os.path.join(RUNS, *parts)
    if create:
        os.makedirs(path, exist_ok=True)
    return path


def install(build=None):
    """The install prefix of `build`, default $MEQ_BUILD or "double"."""
    if build is None and os.environ.get("MEQ_INSTALL"):
        return os.environ["MEQ_INSTALL"]
    name = build or os.environ.get("MEQ_BUILD", "double")
    return os.path.join(PROJECT, "install", name)


def meq(build=None):
    """The installed meq binary."""
    return os.path.join(install(build), "bin", "meq")


def meq_run(build=None):
    """The installed meq-run wrapper."""
    return os.path.join(install(build), "bin", "meq-run")
