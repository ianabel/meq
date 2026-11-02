"""record_libraries.py <module> <output>: the shared libraries a built _meq
resolves to, as meq/_libraries.py.

Run after every link of the module. meq/__init__.py compares them with what
the process already has mapped before it loads the module: a library of the
same name from another file means the module would run on a library it was
not built against (PYTHON-PLAN.md, PY-6), and the import is refused.
"""

import os
import subprocess
import sys

module, output = sys.argv[1:3]
listing = subprocess.run(["ldd", module], capture_output=True, text=True, check=True).stdout

linked = {}
for line in listing.splitlines():
    if "=>" not in line:
        continue
    path = line.split("=>", 1)[1].strip().split(" (", 1)[0].strip()
    if path.startswith("/"):
        linked[os.path.basename(path)] = os.path.realpath(path)

with open(output, "w") as out:
    out.write('"""The shared libraries _meq was linked against, written by\n'
              'python/record_libraries.py at build time. Not hand-edited."""\n\n')
    out.write("LINKED = {\n")
    for name in sorted(linked):
        out.write(f"    {name!r}: {linked[name]!r},\n")
    out.write("}\n")
