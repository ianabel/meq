"""MEQ's inverse solve against freegs4e's, IN THE SAME MODE: both codes given
the same targets, the same fixed coils and the same free supplies, and asked
for the currents.

Every reference fgsref.py writes is an INVERSE solve -- freegs4e re-tunes its
coil currents every Picard pass to put a critical point on each design X-point
and to level each isoflux pair -- and the benchmark's other columns compare
MEQ's FORWARD solve at those currents against it. That measures how sensitive
the equilibrium is to a current set, not how close the two codes are
(CLAUDE_FB.md, *Standing rule: compare in the same mode MEQ solves in*). MEQ
now has an inverse solve of its own (INVERSE-AND-RECONSTRUCTION-PLAN.md Part
One), so the comparison can be made in the mode the reference was made in.

WHICH SUPPLIES ARE FREE, AND WHY NOT ALL OF THEM. freegs4e frees every supply
from zero and its targets do not determine them: DIII-D has 18 supplies and
six target rows, MAST-U 14 and a rank of five. What freegs4e reports in the
null space is whatever its path left there -- on C and G almost the whole
current vector lies in the null space of the vacuum target response -- and
MEQ's inverse makes a different choice there, so all-free currents differ by a
selection, not by an error. So each configuration frees a WELL-POSED subset:
the supplies a column-pivoted QR of the vacuum target response picks first, as
many as its rank, with every other supply held at freegs4e's own value. The
reference is then an exact solution of the restricted problem too, and a
locally unique one, so both codes are answering one question with one answer.

THREE ARMS, each through the shipped driver on the shipped file with only its
currents and [inverse] table changed:

    forward     every supply at freegs4e's current, nothing free: the
                benchmark's existing comparison, re-run on this binary
    inverse     the subset free, SEEDED at freegs4e's currents, so the solve
                starts where the forward one does and the corner moves the
                currents only as far as MEQ's own equilibrium needs them moved
                to meet the targets. The headline arm
    perturbed   the same, with the free supplies seeded --perturb away (10% by
                default): does it come back to the same currents

The targets are fgsref.py's case table: each design X-point as
[[inverse.xpoint]] Track = false, freegs4e's B_R = B_Z = 0 condition, and
each isoflux pair.

THE CURRENT COLUMN is max |I_MEQ - I_ref| over the free supplies, relative to
the largest free |I_ref|. On a shaped configuration the ref-n257 reference --
the same freegs4e solve at twice the resolution -- gives freegs4e's own
resolution sensitivity over the same supplies as a yardstick, `fgs 257`.

USAGE

    venv/bin/python inverse_compare.py OUTDIR [--perturb F] [--arms a,b,c]
                                       [--set table.Key=value ...] [A B ...]

`--set` overrides one key of the shipped file in every arm -- the resolution
study is `--set discretisation.PolynomialDegree=3`,
`--set mesh.generate.PlasmaSize=0.0634` and `--set boundary.exterior.Modes=20`.
A [mesh.generate] override sends the mesh to OUTDIR, so the shipped .msh is
never regenerated in place.

`--guess exact` starts every arm from the REFERENCE's own equilibrium,
rebuilt by Green's functions over MEQ's half-disc (mkexactguess.py at n = 128,
the resolution M-106 needed), with [source] PsiAxis at the reference's
psi_axis. It separates "MEQ cannot solve this problem" from "MEQ's cold start
lands in another basin": the answer is in the starting position.

`--reference DIR` reads the freegs4e references -- currents, psi, axis,
X-points -- from DIR rather than from beside this script.
`fgsref.py --interpolating-profiles` writes a set solved with the analytic
profile MEQ's tables carry, where the committed set solved a smoothing
spline of it that moves ff' by about 1.3e-02; THAT is the same-profile
comparison, and the committed set is the same-mode one with a known per-cent
profile difference left in.

with PYTHONPATH reaching freegs4e, as every script here needs.
"""

import json
import os
import re
import sys

import numpy as np
import scipy.linalg

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, HERE)
os.environ.setdefault("FGSREF_OUT", "/tmp")

from fgsref import CASES                                       # noqa: E402
from circuits import MACHINES, wiring                          # noqa: E402
from inverse import (NUM, rewire, run, reference_rows,         # noqa: E402
                     set_key, target_response)
from make_diverted_case import disc_radius                     # noqa: E402
import subprocess                                              # noqa: E402


def supply_currents(json_path, supplies):
    """Each supply's current in MEQ's sense from a reference's conductors:
    the conductor's own total for a coil alone, the circuit current (total
    over multiplier) for a circuit."""
    cond = {c["label"]: c["current"]
            for c in json.load(open(json_path))["conductors"]}
    return {name: cond[members[0][0]]/members[0][1]
            for name, members in supplies}


def free_subset(A):
    """The supplies a column-pivoted QR of the row-normalised vacuum response
    takes first, as many as its numerical rank."""
    Ar = A/np.linalg.norm(A, axis=1, keepdims=True)
    s = np.linalg.svd(Ar, compute_uv=False)
    rank = int(np.sum(s > 1.0e-8*s[0]))
    _, _, piv = scipy.linalg.qr(Ar, pivoting=True)
    return sorted(piv[:rank].tolist()), rank, s


def put_key(text, table, key, value):
    """set_key, or the key appended to the table when the file has none."""
    m = re.search(r"(?m)^\[%s\]\s*$" % re.escape(table), text)
    if m and not re.search(r"(?m)^%s = " % re.escape(key),
                           text[m.end():m.end() + (re.search(r"(?m)^\[",
                                text[m.end():]) or re.search(r"\Z",
                                text[m.end():])).start()]):
        return text[:m.end()] + "\n%s = %s" % (key, value) + text[m.end():]
    return set_key(text, table, key, value)


def psi_bnd_rel(text, npz):
    m = re.search(r"psi_bnd = " + NUM, text)
    ref = float(npz["psi_bndry"])
    return abs(float(m.group(1)) - ref)/abs(ref) if m else float("nan")


def solved_currents(text):
    return {label: float(value) for label, value in re.findall(
        r"(?m)^\s{7}(\S+)\s+([-+]\d\.\d+e[-+]\d+) A  from", text)}


def main():
    if len(sys.argv) < 2:
        raise SystemExit(__doc__)
    out = os.path.abspath(sys.argv[1])
    args = sys.argv[2:]
    perturb = 0.10
    if "--perturb" in args:
        at = args.index("--perturb")
        perturb = float(args[at + 1])
        del args[at:at + 2]
    arms = ["forward", "inverse", "perturbed"]
    if "--arms" in args:
        at = args.index("--arms")
        arms = args[at + 1].split(",")
        del args[at:at + 2]
    refdir = HERE
    if "--reference" in args:
        at = args.index("--reference")
        refdir = os.path.abspath(args[at + 1])
        del args[at:at + 2]
    guess = "cold"
    if "--guess" in args:
        at = args.index("--guess")
        guess = args[at + 1]
        del args[at:at + 2]
    if guess not in ("cold", "exact"):
        raise SystemExit("--guess is cold or exact")
    overrides = []
    while "--set" in args:
        at = args.index("--set")
        key, value = args[at + 1].split("=", 1)
        table, key = key.rsplit(".", 1)
        overrides.append((table, key, value))
        del args[at:at + 2]
    want = args
    rel_out = os.path.relpath(out, ROOT)
    os.makedirs(out, exist_ok=True)
    cases = {c["name"]: c for c in CASES}
    results = []

    print("%-26s %-9s | %4s %3s %9s | %9s %8s | %9s %9s %8s %8s" % (
        "config", "arm", "exit", "its", "target", "dI/I", "fgs 257",
        "psi_ax", "psi_bnd", "axis m", "X m"), flush=True)
    for ref, stem in MACHINES:
        letter = ref[0]
        if want and letter not in want:
            continue
        case = cases[ref]
        wires = wiring(case["make"]())
        for suffix in ("", "-shaped"):
            name = "machine-%s%s" % (stem[len("machine-"):], suffix)
            npz_name = ref + suffix.replace("-", "_")
            npz_path = os.path.join(refdir, npz_name + ".npz")
            if not os.path.exists(npz_path):
                continue        # a --reference set need not carry every case
            npz = np.load(npz_path, allow_pickle=True)
            shipped = open(os.path.join(ROOT, "examples", name + ".toml")).read()

            A, supplies, _ = target_response(npz_path, case, wires)
            current = supply_currents(os.path.join(refdir, npz_name + ".json"),
                                      supplies)
            pick, rank, sv = free_subset(A)
            free = [supplies[i][0] for i in pick]
            scale = max(abs(current[f]) for f in free)
            fine = os.path.join(refdir, "ref-n257", npz_name + ".json")
            fgs257 = float("nan")
            if os.path.exists(fine):
                c257 = supply_currents(fine, supplies)
                fgs257 = max(abs(c257[f] - current[f]) for f in free)/scale

            if guess == "exact":
                gmesh = os.path.join(out, name + "-exact.mesh")
                ggf = os.path.join(out, name + "-exact.gf")
                made = subprocess.run(
                    [sys.executable, os.path.join(HERE, "mkexactguess.py"),
                     npz_path, gmesh, ggf, "%g" % disc_radius(npz), "128"],
                    capture_output=True, text=True)
                if made.returncode != 0:
                    raise SystemExit(made.stderr)

            for arm in arms:
                seeds = dict(current)
                if arm == "perturbed":
                    for k, f in enumerate(free):
                        seeds[f] = current[f] + (perturb if k % 2 == 0
                                                 else -perturb)*scale
                text = rewire(shipped, wires, seeds,
                              free=set() if arm == "forward" else set(free))
                if guess == "exact":
                    text = set_key(text, "initialguess", "File", '"%s"' % ggf)
                    text = set_key(text, "initialguess", "MeshFile",
                                   '"%s"' % gmesh)
                    text = set_key(text, "source", "PsiAxis",
                                   "%.10e" % float(npz["psi_axis"]))
                for table, key, value in overrides:
                    text = put_key(text, table, key, value)
                if any(t == "mesh.generate" for t, _, _ in overrides):
                    text = set_key(text, "mesh", "File", '"%s"' % os.path.join(
                        rel_out, "%s-%s.msh" % (name, arm)))
                text = set_key(text, "output", "Prefix", '"%s"' % os.path.join(
                    rel_out, "%s-%s" % (name, arm)))
                head = ("# %s, %s arm, by inverse_compare.py: edit the script,"
                        " not this.\n# Free: %s\n\n" % (
                            name, arm, ", ".join(free) if arm != "forward"
                            else "nothing"))
                if arm != "forward":
                    targets = ["", "[inverse]", ""]
                    for R, Z in case["xpoints"]:
                        targets += ["[[inverse.xpoint]]", "R = %r" % R,
                                    "Z = %r" % Z, "Track = false", ""]
                    for R1, Z1, R2, Z2 in case["isoflux"]:
                        targets += ["[[inverse.isoflux]]", "R1 = %r" % R1,
                                    "Z1 = %r" % Z1, "R2 = %r" % R2,
                                    "Z2 = %r" % Z2, ""]
                    text = text.rstrip("\n") + "\n" + "\n".join(targets)
                toml = os.path.join(out, "%s-%s.toml" % (name, arm))
                open(toml, "w").write(head + text)
                code, wall, log = run(toml, os.path.join(
                    out, "%s-%s.log" % (name, arm)))

                m = re.search(r"largest target residual " + NUM, log)
                target = float(m.group(1)) if m else float("nan")
                m = re.search(r"converged in (\d+) Newton", log)
                its = m.group(1) if m else "-"
                row = reference_rows(log, npz)
                solved = solved_currents(log)
                dI = float("nan")
                if arm != "forward" and all(f in solved for f in free):
                    dI = max(abs(solved[f] - current[f]) for f in free)/scale
                result = dict(config=name, arm=arm, exit=code, its=its,
                              target=target, dI=dI, fgs257=fgs257,
                              psi_ax=row["rel"], psi_bnd=psi_bnd_rel(log, npz),
                              axis=row["dax"], X=row["dX"], wall=wall,
                              located=row["located"], free=free, rank=rank,
                              solved={f: solved.get(f) for f in free},
                              reference={f: current[f] for f in free})
                results.append(result)
                print("%-26s %-9s | %4s %3s %9.2e | %9.2e %8.1e | %9.2e %9.2e "
                      "%8.2e %8.2e%s" % (
                          name, arm, code, its, target, dI, fgs257,
                          row["rel"], result["psi_bnd"], row["dax"], row["dX"],
                          "" if row["located"] else "  NO AXIS"), flush=True)
            json.dump(results, open(os.path.join(out, "results.json"), "w"),
                      indent=1)
    print("COMPARE DONE", flush=True)


if __name__ == "__main__":
    main()
