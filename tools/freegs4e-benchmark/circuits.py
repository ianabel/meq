"""Each machine's conductors twice: flattened, and wired as [[circuits]].

INVERSE-AND-RECONSTRUCTION-PLAN.md IV-1's acceptance, written out as fixtures.
For every reference this emits two TOML fragments over the SAME rectangles:

    <stem>-flattened.toml   [[coils]] with every conductor's total current,
                            which is what examples/machine-*.toml carry
    <stem>-circuited.toml   [[coils]] with no current on a circuit member, and
                            [[circuits]] carrying one current per supply and a
                            signed Multiplier per member

and tests/unit/ConfigTests.cpp asserts the two parse to the same coil set.
That is the test that the circuit layer cannot change a FORWARD answer: a
circuit is a notation there, and only an inverse solve -- which frees the
circuit's one current rather than each conductor's -- can tell them apart.

WHERE THE MULTIPLIER COMES FROM.  conductors.flatten() gives a circuit member
the total `turns * multiplier * I_circuit`, so MEQ's Multiplier is
`turns * multiplier`: the share of the supply's current through that
rectangle.  A solenoid is one rectangle carrying `Ns * I`, so a solenoid inside
a circuit has Multiplier `Ns * multiplier`, and a solenoid on its own supply is
its own circuit with Multiplier `Ns`.

WHAT IS A CIRCUIT HERE.  Exactly what the machine definition calls one: a
freegs4e `Circuit`.  A top-level Coil, ShapedCoil or Solenoid is its own supply
and is written as a plain coil with its own Current, because a coil in no
circuit is its own one-member circuit -- which is what keeps every existing
file valid.  So TestTokamak and DIII-D, which have no Circuit, produce two
identical fragments, and the comparison is non-trivial on MAST, TCV and MAST-U.

THE CIRCUIT CURRENT IS RECOVERED FROM THE REFERENCE, from its first member, as
`total / Multiplier`.  Every number is written with seventeen significant
digits, so the round trip through a double is exact and the two fragments
differ only in the order of one multiplication.

USAGE

    python3 circuits.py OUTDIR            # every reference's filament arm
"""

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
os.environ.setdefault("FGSREF_OUT", "/tmp")

from fgsref import CASES                       # noqa: E402

MACHINES = [("A_testtokamak_classic", "machine-a-testtokamak"),
            ("B_testtokamak_peaked_ffdom", "machine-b-ffprime"),
            ("C_mast_spherical", "machine-c-mast"),
            ("D_tcv_conventional", "machine-d-tcv"),
            ("E_testtokamak_diamagnetic", "machine-e-diamagnetic"),
            ("F_diiid_conventional", "machine-f-diiid"),
            ("G_mastu_simple", "machine-g-mastu")]


def wiring(machine):
    """( label, circuit name or None, MEQ multiplier ) per rectangle, in
    conductors.flatten()'s depth-first order."""
    out = []

    def share(item):
        cls = type(item).__name__
        if cls == "Solenoid":
            return float(item.Ns)
        return float(getattr(item, "turns", 1.0))

    for label, item in machine.coils:
        if type(item).__name__ == "Circuit":
            for name, coil, multiplier in item.coils:
                if type(coil).__name__ == "Circuit":
                    raise SystemExit("%s: a circuit inside a circuit is not "
                                     "something [[circuits]] expresses" % label)
                out.append(("%s.%s" % (label, name), label,
                            float(multiplier)*share(coil)))
        else:
            out.append((label, None, share(item)))
    return out


def coil_block(w, label, R, Z, hw, hh, current=None):
    w("[[coils]]")
    w('Name = "%s"' % label)
    w("CentreR = %.17g" % R)
    w("CentreZ = %.17g" % Z)
    w("HalfWidth = %.17g" % hw)
    w("HalfHeight = %.17g" % hh)
    if current is not None:
        w("Current = %.17g" % current)
    w("")


def main():
    if len(sys.argv) != 2:
        raise SystemExit(__doc__)
    outdir = sys.argv[1]
    os.makedirs(outdir, exist_ok=True)
    cases = {c["name"]: c for c in CASES}

    for ref, stem in MACHINES:
        d = np.load(os.path.join(HERE, ref + ".npz"), allow_pickle=True)
        labels = [str(x) for x in d["coil_labels"]]
        R, Z = d["coil_R"], d["coil_Z"]
        hw, hh = d["coil_half_width"], d["coil_half_height"]
        total = np.asarray(d["coil_currents"], float)

        wires = wiring(cases[ref]["make"]())
        if [x[0] for x in wires] != labels:
            raise SystemExit("%s: the machine's own order %s is not the "
                             "reference's %s" % (ref, [x[0] for x in wires],
                                                  labels))

        flat, circ = [], []
        circuits = {}
        for i, (label, circuit, mult) in enumerate(wires):
            coil_block(flat.append, label, R[i], Z[i], hw[i], hh[i], total[i])
            if circuit is None:
                coil_block(circ.append, label, R[i], Z[i], hw[i], hh[i],
                           total[i])
                continue
            coil_block(circ.append, label, R[i], Z[i], hw[i], hh[i])
            if circuit not in circuits:
                circuits[circuit] = dict(current=total[i]/mult, members=[])
            circuits[circuit]["members"].append((label, mult))

        for name, c in circuits.items():
            circ.append("[[circuits]]")
            circ.append('Name = "%s"' % name)
            circ.append("Current = %.17g" % c["current"])
            members = ", ".join('{ Name = "%s", Multiplier = %.17g }' % m
                                for m in c["members"])
            circ.append("Coils = [ %s ]" % members)
            circ.append("")

        head = ("# %s, from tools/freegs4e-benchmark/%s.npz by circuits.py.\n"
                "# Generated: edit the script, not this.\n\n" % (stem, ref))
        for kind, lines in (("flattened", flat), ("circuited", circ)):
            path = os.path.join(outdir, "%s-%s.toml" % (stem, kind))
            with open(path, "w") as f:
                f.write(head + "\n".join(lines))
        print("%-24s %2d conductors, %2d circuits" % (stem, len(wires),
                                                     len(circuits)))


if __name__ == "__main__":
    main()
