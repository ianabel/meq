"""IV-6: MEQ's inverse solve from cold on the fourteen benchmark configurations,
then RELEASED -- re-solved forward with the currents frozen where the inverse
left them. That is forward.py's solve_from_cold() chain, run on MEQ.

INVERSE-AND-RECONSTRUCTION-PLAN.md section 9. IV-0 (MEASUREMENTS.md M-244) ran
the same chain in freegs4e and found the basin on thirteen of fourteen
configurations where a cold forward solve found it on none; this is section
2.2's claim tested on MEQ.

WHAT EACH RUN IS GIVEN, AND IT IS WHAT IV-0 GAVE freegs4e:

    the machine, WIRED            examples/machine-*.toml's own rectangles, with
                                  [[circuits]] as freegs4e's Circuit wires them
                                  -- circuits.py's wiring() -- and EVERY supply
                                  Free, as forward.py sets control = True on
                                  every top-level item
    the currents: ZERO, or the    freegs4e's constructors default every current
    VACUUM STAGE's                to zero, and IV-0 started there; that is
                                  `--seed zero`. The default, `--seed vacuum`,
                                  is what freegs4e's control does on its first
                                  step from there: the minimum-norm currents
                                  that best meet the targets against the
                                  guess's own plasma, by linear least squares
                                  on the circuits' Green's functions -- B_R and
                                  B_Z at each design X-point, psi1 - psi2 per
                                  isoflux pair, freegs4e's rows. Nothing solved
                                  for is read either way
    the targets                   fgsref.py's own case table: each design
                                  X-point as [[inverse.xpoint]] Track = false --
                                  a critical point there, which is freegs4e's
                                  X-point constraint -- and each isoflux pair
    the guess                     the blob plus the coils at the seed currents,
                                  mkcoldguess.py --currents; at zero seeds
                                  that is the blob alone, --remainder

WHY THE VACUUM STAGE EXISTS, MEASURED: from zero currents the guess has no
separatrix anywhere, so MEQ's first support sweep freezes EVERY candidate
element as plasma (2124 of 2124 on DIII-D), and the inverse converges -- every
target met to 1e-17 -- on a topology whose axis is a vacuum O-point beyond
the outboard coils, ( 2.97, -1.38 ). freegs4e re-decides its core every Picard
step and never freezes one, which is why IV-0's zero start does not meet
this. MEASUREMENTS.md M-253.

Everything else -- mesh, profiles, plasma current, the X-point seed, the
exterior, the support sweeps -- is the shipped file's, unchanged.

THE RELEASE copies the shipped file again, writes the inverse's circuit
currents in as data with nothing Free, and restarts from the inverse's own
answer on the same mesh. A release that converges and keeps its axis is the
basin found; the reference columns say whether it is the REFERENCE's basin.

USAGE

    python3 inverse.py OUTDIR [--seed zero|vacuum] [--gamma G] [A B ...]

`--gamma` is [inverse] Gamma, the regularisation of the step in every free
current in flux per ampere; zero, the default, is none.

Writes OUTDIR/<config>-{inverse,release}.{toml,log} and the outputs beside
them, each run made with OUTDIR as its working directory, and prints one row
per configuration -- so OUTDIR belongs under runs/freegs4e-benchmark/.
"""

import os
import re
import subprocess
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
os.environ.setdefault("FGSREF_OUT", "/tmp")

import benchpaths                              # noqa: E402
from fgsref import CASES                       # noqa: E402
from circuits import MACHINES, wiring          # noqa: E402
from make_diverted_case import disc_radius     # noqa: E402
from mkcoldguess import plasma_filaments, with_design  # noqa: E402
from mkexactguess import greens                # noqa: E402
import conductors                              # noqa: E402
import json                                    # noqa: E402

NUM = r"([-+]?\d\.\d+e[-+]\d+)"
TIMEOUT = 2400


def coil_blocks(text):
    """( start, end, label ) of every [[coils]] block, in file order."""
    out = []
    heads = list(re.finditer(r"(?m)^\[\[coils\]\][ \t]*\n", text))
    for head in heads:
        s = head.start()
        nxt = re.search(r"(?m)^\[", text[head.end():])
        e = head.end() + nxt.start() if nxt else len(text)
        label = re.search(r'(?m)^Name = "([^"]+)"', text[s:e]).group(1)
        out.append((s, e, label))
    return out


def rewire(text, wires, currents, free):
    """The shipped file's [[coils]] with every current replaced: a circuit
    member carries none, a coil on its own supply carries `currents[label]`,
    and [[circuits]] blocks are appended with theirs. `free` marks every
    supply Free when True, none when False, and those it names when it is a
    set of supply names."""
    if free is True or free is False:
        free = {label for label, _, _ in wires} | {c for _, c, _ in wires if c} \
            if free else set()
    owner = {label: (circuit, mult) for label, circuit, mult in wires}
    pieces, last = [], 0
    for s, e, label in coil_blocks(text):
        block = text[s:e]
        block = re.sub(r"(?m)^(Current|CurrentDensity|Free) = .*\n", "", block)
        circuit, _ = owner[label]
        if circuit is None:
            extra = "Current = %.17g\n" % currents[label]
            if label in free:
                extra += "Free = true\n"
            body = block.rstrip("\n") + "\n" + extra + "\n"
        else:
            body = block.rstrip("\n") + "\n\n"
        pieces.append(text[last:s])
        pieces.append(body)
        last = e
    pieces.append(text[last:])
    out = "".join(pieces)

    circuits = {}
    for label, circuit, mult in wires:
        if circuit is not None:
            circuits.setdefault(circuit, []).append((label, mult))
    lines = []
    for name, members in circuits.items():
        lines.append("[[circuits]]")
        lines.append('Name = "%s"' % name)
        lines.append("Current = %.17g" % currents[name])
        if name in free:
            lines.append("Free = true")
        lines.append("Coils = [ %s ]" % ", ".join(
            '{ Name = "%s", Multiplier = %.17g }' % m for m in members))
        lines.append("")
    return out.rstrip("\n") + "\n\n" + "\n".join(lines) + "\n"


def set_key(text, table, key, value):
    """Replace `key = ...` inside `[table]`, which must hold it."""
    m = re.search(r"(?m)^\[%s\]\s*$" % re.escape(table), text)
    if not m:
        raise SystemExit("no [%s] table" % table)
    nxt = re.search(r"(?m)^\[", text[m.end():])
    end = m.end() + nxt.start() if nxt else len(text)
    body, n = re.subn(r"(?m)^%s = .*$" % re.escape(key), "%s = %s" % (key, value),
                      text[m.end():end])
    if n != 1:
        raise SystemExit("[%s] %s found %d times" % (table, key, n))
    return text[:m.end()] + body + text[end:]


def drop_key(text, table, key):
    """Remove `key = ...` from `[table]` if it is there. A `.nc` restart is
    refused beside MeshFile and Content -- the file carries its mesh and says
    what it holds -- so a shipped .gf restart's keys have to go."""
    m = re.search(r"(?m)^\[%s\]\s*$" % re.escape(table), text)
    if not m:
        return text
    nxt = re.search(r"(?m)^\[", text[m.end():])
    end = m.end() + nxt.start() if nxt else len(text)
    body = re.sub(r"(?m)^%s = .*\n?" % re.escape(key), "", text[m.end():end])
    return text[:m.end()] + body + text[end:]


def target_response(npz_path, case, wires):
    """The vacuum response of freegs4e's target rows to one ampere in each
    supply -- B_R and B_Z at each design X-point, psi1 - psi2 per isoflux
    pair -- and the same rows of the design blob's own field. Returns
    ( A, supplies, blob rows ), A having one column per supply in `supplies`'
    order, each supply a ( name, [ ( label, multiplier ) ] )."""
    d = with_design(np.load(npz_path), npz_path)
    flat = {c["label"]: c for c in conductors.from_npz(d)}
    # THE SAME FILAMENTS mkcoldguess.py sums, so the guess written from these
    # currents is the field the currents were chosen against.
    filaments = {}
    for label, c in flat.items():
        nz = max(1, int(round(c["half_height"]/max(c["half_width"], 1e-9))))
        edges = np.linspace(c["Z"] - c["half_height"], c["Z"] + c["half_height"],
                            nz + 1)
        filaments[label] = [(c["R"], float(z), 1.0/nz,
                             max(c["half_width"], c["half_height"]/nz))
                            for z in 0.5*(edges[:-1] + edges[1:])]
    src_R, src_Z, src_I, soft, _ = plasma_filaments(d, npz_path)

    def blob(R, Z):
        return float(np.asarray(greens(src_R, src_Z, R, Z, soft=soft)) @ src_I)

    supplies = []
    for label, circuit, mult in wires:
        if circuit is None:
            supplies.append((label, [(label, 1.0)]))
        elif not supplies or supplies[-1][0] != circuit:
            supplies.append((circuit, [(label, mult)]))
        else:
            supplies[-1][1].append((label, mult))

    def unit(members):
        def f(R, Z):
            return sum(m*w*float(greens(rc, zc, R, Z, soft=s))
                       for label, m in members
                       for rc, zc, w, s in filaments[label])
        return f

    def rows(f):
        h = 1.0e-5
        out = []
        for R, Z in case["xpoints"]:
            out.append((f(R + h, Z) - f(R - h, Z))/(2*h)/R)
            out.append((f(R, Z + h) - f(R, Z - h))/(2*h)/R)
        for R1, Z1, R2, Z2 in case["isoflux"]:
            out.append(f(R1, Z1) - f(R2, Z2))
        return out

    A = np.array([rows(unit(m)) for _, m in supplies]).T
    return A, supplies, np.array(rows(blob))


def vacuum_currents(npz_path, case, wires):
    """freegs4e's first control step from zero currents: the minimum-norm
    supply currents that best meet the case's targets against the design
    blob's field. Returns ( per supply, per conductor total )."""
    A, supplies, blob = target_response(npz_path, case, wires)
    b = -blob
    x = np.linalg.lstsq(A, b, rcond=None)[0]
    per_supply = {name: float(v) for (name, _), v in zip(supplies, x)}
    per_conductor = {}
    for (name, members), v in zip(supplies, x):
        for label, m in members:
            per_conductor[label] = float(m*v)
    return per_supply, per_conductor


def run(toml, log):
    """meq-run on `toml`, made in the directory the file is written to: its
    Prefix, and any mesh [mesh.generate] makes, are relative to that."""
    env = dict(os.environ, MKL_NUM_THREADS="1", OMP_WAIT_POLICY="passive")
    t0 = time.time()
    with open(log, "w") as fh:
        try:
            code = subprocess.run(
                [benchpaths.meq_run(), toml],
                cwd=os.path.dirname(os.path.abspath(toml)), env=env,
                stdout=fh, stderr=subprocess.STDOUT, timeout=TIMEOUT).returncode
        except subprocess.TimeoutExpired:
            code = "T/O"
    return code, time.time() - t0, open(log).read()


def reference_rows(text, npz):
    """psi_ax relative error, axis and X-point distances against the npz, and
    whether the run reported a located axis -- the cold table's columns."""
    rpa = float(npz["psi_axis"])
    m = re.search(r"psi_ax = " + NUM, text)
    pa = float(m.group(1)) if m else float("nan")
    m = re.search(r"the axis, as a zero of q_h: psi = " + NUM
                  + r" at \(\s*([-\d.]+),\s*([-\d.]+)\s*\)", text)
    dax = float(np.hypot(float(m.group(2)) - float(npz["Raxis"]),
                         float(m.group(3)) - float(npz["Zaxis"]))) \
        if m else float("nan")
    x = re.search(r"psi_bnd = " + NUM + r".*?\(\s*([-\d.]+),\s*([-\d.]+)\s*\)",
                  text)
    if x:
        refx = np.asarray(npz["xpoints"], float)
        dX = float(np.min(np.hypot(refx[:, 0] - float(x.group(2)),
                                   refx[:, 1] - float(x.group(3))))) \
            if refx.size else float("nan")
        xR, xZ = float(x.group(2)), float(x.group(3))
    else:
        dX, xR, xZ = float("nan"), float("nan"), float("nan")
    return dict(pa=pa, rel=abs(pa - rpa)/abs(rpa), dax=dax, dX=dX,
                located=m is not None, xR=xR, xZ=xZ)


def main():
    if len(sys.argv) < 2:
        raise SystemExit(__doc__)
    out = os.path.abspath(sys.argv[1])
    seed = "vacuum"
    args = sys.argv[2:]
    if "--seed" in args:
        at = args.index("--seed")
        seed = args[at + 1]
        del args[at:at + 2]
    if seed not in ("zero", "vacuum"):
        raise SystemExit("--seed is zero or vacuum")
    gamma = 0.0
    if "--gamma" in args:
        at = args.index("--gamma")
        gamma = float(args[at + 1])
        del args[at:at + 2]
    want = args
    # The driver writes Directory/Prefix, and the runs are made in `out`.
    os.makedirs(out, exist_ok=True)
    cases = {c["name"]: c for c in CASES}

    print("%-28s | %4s %3s %9s %4s %9s %6s | %4s %4s %10s %8s %8s %6s" % (
        "config", "inv", "its", "target", "null", "cond", "wall",
        "rel", "its", "psi_ax rel", "axis m", "X m", "wall"), flush=True)
    for ref, stem in MACHINES:
        letter = ref[0]
        if want and letter not in want:
            continue
        case = cases[ref]
        wires = wiring(case["make"]())
        for suffix in ("", "-shaped"):
            name = "machine-%s%s" % (stem[len("machine-"):], suffix)
            npz_name = ref + suffix.replace("-", "_")
            npz_path = os.path.join(benchpaths.OUT, npz_name + ".npz")
            npz = np.load(npz_path, allow_pickle=True)
            shipped = benchpaths.read_example(
                os.path.join(benchpaths.EXAMPLES, name + ".toml"))

            # THE SEEDS AND THE GUESS.
            zero = {label: 0.0 for label, circuit, _ in wires if circuit is None}
            zero.update({circuit: 0.0 for _, circuit, _ in wires if circuit})
            gmesh = os.path.join(out, name + "-guess.mesh")
            ggf = os.path.join(out, name + "-guess.gf")
            command = [sys.executable, os.path.join(HERE, "mkcoldguess.py"),
                       npz_path, gmesh, ggf, "%g" % disc_radius(npz), "32"]
            if seed == "zero":
                seeds = zero
                command.append("--remainder")
            else:
                seeds, totals = vacuum_currents(npz_path, case, wires)
                jpath = os.path.join(out, name + "-vacuum.json")
                json.dump(totals, open(jpath, "w"), indent=1)
                command += ["--currents", jpath]
            made = subprocess.run(command, capture_output=True, text=True)
            if made.returncode != 0:
                raise SystemExit(made.stderr)
            peak = [l for l in made.stdout.splitlines() if l.startswith("PsiAxis ")]
            psi_axis = float(peak[0].split()[1])

            # THE INVERSE.
            inv = rewire(shipped, wires, seeds, free=True)
            inv = set_key(inv, "source", "PsiAxis", "%.10e" % psi_axis)
            inv = set_key(inv, "initialguess", "File", '"%s"' % ggf)
            inv = set_key(inv, "initialguess", "MeshFile", '"%s"' % gmesh)
            inv = set_key(inv, "output", "Prefix",
                          '"%s"' % (name + "-inverse"))
            targets = ["", "[inverse]", "Gamma = %r" % gamma, ""]
            for R, Z in case["xpoints"]:
                targets += ["[[inverse.xpoint]]", "R = %r" % R, "Z = %r" % Z,
                            "Track = false", ""]
            for R1, Z1, R2, Z2 in case["isoflux"]:
                targets += ["[[inverse.isoflux]]", "R1 = %r" % R1, "Z1 = %r" % Z1,
                            "R2 = %r" % R2, "Z2 = %r" % Z2, ""]
            inv = ("# IV-6: %s posed INVERSE from %s currents by inverse.py.\n"
                   "# Generated: edit the script, not this.\n\n" % (name, seed)
                   + inv.rstrip("\n") + "\n" + "\n".join(targets))
            itoml = os.path.join(out, name + "-inverse.toml")
            open(itoml, "w").write(inv)
            icode, iwall, itext = run(itoml, os.path.join(out, name + "-inverse.log"))

            m = re.search(r"largest target residual " + NUM, itext)
            target = float(m.group(1)) if m else float("nan")
            m = re.search(r"corner rank (\d+) nullity (\d+)(?:, condition " + NUM + r")?",
                          itext)
            null = m.group(2) if m else "-"
            cond = float(m.group(3)) if m and m.group(3) else float("nan")
            m = re.search(r"converged in (\d+) Newton", itext)
            iits = m.group(1) if m else "-"
            solved = {}
            for label, value in re.findall(
                    r"(?m)^\s{7}(\S+)\s+([-+]\d\.\d+e[-+]\d+) A  from", itext):
                solved[label] = float(value)
            irow = reference_rows(itext, npz)

            row = "%-28s | %4s %3s %9.2e %4s %9.2e %6.1f |" % (
                name, icode, iits, target, null, cond, iwall)
            if icode not in (0, 4) or len(solved) != len(zero):
                print(row + "  inverse did not finish: no release", flush=True)
                continue

            # THE RELEASE: the inverse's currents as data, from its own answer.
            rel = rewire(shipped, wires, solved, free=False)
            rel = set_key(rel, "source", "PsiAxis", "%.10e" % irow["pa"])
            # The inverse's own answer, from its .nc: MEQ's restart, carrying
            # its mesh and saying what it holds.
            rel = set_key(rel, "initialguess", "File",
                          '"%s"' % os.path.join(out, name + "-inverse.nc"))
            rel = drop_key(rel, "initialguess", "MeshFile")
            rel = drop_key(rel, "initialguess", "Content")
            if np.isfinite(irow["xR"]):
                rel = set_key(rel, "boundary.xpoint", "R", "%.10f" % irow["xR"])
                rel = set_key(rel, "boundary.xpoint", "Z", "%.10f" % irow["xZ"])
            rel = set_key(rel, "output", "Prefix",
                          '"%s"' % (name + "-release"))
            rel = ("# IV-6: %s RELEASED -- the inverse's currents frozen.\n"
                   "# Generated by inverse.py: edit the script, not this.\n\n"
                   % name + rel)
            rtoml = os.path.join(out, name + "-release.toml")
            open(rtoml, "w").write(rel)
            rcode, rwall, rtext = run(rtoml, os.path.join(out, name + "-release.log"))
            m = re.search(r"converged in (\d+) Newton", rtext)
            rits = m.group(1) if m else "-"
            rrow = reference_rows(rtext, npz)
            tag = "" if rrow["located"] else "  NO AXIS"
            print(row + " %4s %4s %10.3e %8.2e %8.2e %6.1f%s" % (
                rcode, rits, rrow["rel"], rrow["dax"], rrow["dX"], rwall, tag),
                flush=True)
    print("INVERSE DONE", flush=True)


if __name__ == "__main__":
    main()
