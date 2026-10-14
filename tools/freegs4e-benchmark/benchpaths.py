"""Where this benchmark reads and writes: tools/meqpaths.py, plus this
directory's own layout and the one rewrite every variant TOML needs.

NOTHING IS WRITTEN BESIDE THESE SCRIPTS.  Every product goes under

    OUT = runs/freegs4e-benchmark/
        <case>.npz, <case>.json     fgsref.py's 129^2 references, and its
        summary.json, summary_shaped.json
        ref-n<N>/                   fgsref.py --nx=N, unless FGSREF_OUT says
        ref-n257-fixed-boundary/    the 257^2 set examples/fixed-*.toml came from
        ref-interp/                 fgsref.py --interpolating-profiles, and the
                                    inverse race's problem.json
        race-runs/<stage>/          race_inverse.py run <stage>
        <script>/                   each other script's default output

which is the old layout of this directory with the scripts taken out.  The
venv stays here: it is an environment, not output.

THE TRACKED PROBLEM STATEMENTS ARE THE EXCEPTION, and are source:
make_all.sh, make_fixed.sh and the generators they drive write
examples/*.toml, their .dat tables and their guesses into source/examples.
"""

import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import meqpaths                                                 # noqa: E402

SOURCE = meqpaths.SOURCE
EXAMPLES = meqpaths.source("examples")
OUT = meqpaths.runs("freegs4e-benchmark")
VENV = os.path.join(HERE, "venv", "bin", "python")
FREEGS4E = meqpaths.sibling("freegs4e")
FREEGSNKE = meqpaths.sibling("freegsnke")
NICE = meqpaths.sibling("nice")
LIBMFEM = meqpaths.sibling("mfem", "install", "lib", "libmfem.a")


def out(*parts):
    """A directory under OUT, made if it is not there."""
    return meqpaths.runs("freegs4e-benchmark", *parts)


def meq():
    """The installed meq."""
    return meqpaths.meq()


def meq_run():
    """The installed meq-run."""
    return meqpaths.meq_run()


# THE KEYS meq::Configuration RESOLVES AGAINST THE FILE'S OWN DIRECTORY.  A
# variant written anywhere but beside its original reads them from the wrong
# place, so absolutise() pins them first.  [mesh] File is an input only when
# there is no [mesh.generate]: a generated mesh is the run's product and stays
# relative to the working directory, as do [output] Directory and
# [conductors] CacheFile.
_ANYWHERE = ("PPrimeFile", "GGPrimeFile", "SafetyFactorFile", "OmegaFile",
             "TemperatureFile", "DensityFile")
_TABLE = re.compile(r"^\s*\[\[?\s*([^\]\s]+)\s*\]\]?")
_STRING = re.compile(r'^(\s*)([A-Za-z]+)(\s*=\s*)"([^"]*)"(.*)$')


def absolutise(text, directory):
    """`text`, a TOML read from `directory`, with every relative INPUT path
    made absolute, so the file can be written and run from anywhere."""
    directory = os.path.abspath(directory)
    generated = re.search(r"(?m)^\s*\[mesh\.generate\]", text) is not None
    table, lines = "", []
    for line in text.split("\n"):
        heading = _TABLE.match(line)
        if heading:
            table = heading.group(1)
        m = _STRING.match(line)
        if m and m.group(4) and not os.path.isabs(m.group(4)):
            key = m.group(2)
            if (key in _ANYWHERE
                    or (table == "initialguess" and key in ("File", "MeshFile"))
                    or (table == "mesh" and key == "File" and not generated)):
                line = '%s%s%s"%s"%s' % (
                    m.group(1), key, m.group(3),
                    os.path.normpath(os.path.join(directory, m.group(4))),
                    m.group(5))
        lines.append(line)
    return "\n".join(lines)


def read_example(path):
    """A TOML's text with its input paths absolute; see absolutise()."""
    with open(path) as handle:
        return absolutise(handle.read(), os.path.dirname(os.path.abspath(path)))
