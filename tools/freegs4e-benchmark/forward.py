"""FORWARD-MODE freegs4e inputs: the same problem MEQ is given.

WHY THIS EXISTS, AND IT IS A STANDING RULE RATHER THAN A CONVENIENCE.

`fgsref.py` builds every reference through `SyncConstrain( xpoints, isoflux )`
-- freegs4e's INVERSE solve, which re-solves the coil currents on every Picard
pass to hit its geometric targets.  The currents are an OUTPUT of the reference
and an INPUT to MEQ, which solves FORWARD.  So the two codes are not answering
the same question, and the residual between them is not a measurement of either
solver.  `CLAUDE_FB.md`, *Standing rule: compare in the same mode MEQ solves
in*, has the full argument and the measurement that forced it.

The tell: machine A's filament and ShapedCoil references differ in conductor
model and agree with EACH OTHER on `psi_axis` to 8e-06, because the inverse
solve re-tunes to the same targets either way, while MEQ sits 2.4e-02 from
both.

WHAT A FORWARD INPUT IS.  The reference's own converged coil currents, frozen;
its own profile arrays; its own grid and wall; and no constraint object.  The
plasma current is still PRESCRIBED, because MEQ prescribes it -- `[source]
PlasmaCurrent` is a row of the bordered Newton -- so `Ip_logic` stays on.

EVERYTHING COMES OUT OF THE `.npz`, which is what makes these inputs
self-contained and what makes them describe the same machine MEQ meshes:
`coil_R/Z/half_width/half_height/currents/kind/labels` are the conductor table
`conductors.py` hands to `make_diverted_case.py`, `R`/`Z`/`order` are the grid,
`wall_R/Z` the wall, and `psi_n`/`pprime`/`ffprime` the source with stage 3's
`Ip_logic` scale already folded in.

THREE TRAPS, ALL THREE MET WHILE WRITING THIS, ALL THREE SILENT.

  1. `Machine.getPsitokamak()` -- what `Equilibrium.psi()` uses -- reads
     `Machine.current_vec`, NOT `coil.current`.  Setting the currents without
     `getCurrentsVec()` leaves the grid psi with no coil field at all;
     `find_critical()` then reports no X-points and Picard converges happily to
     a small plasma against the outboard wall.  Measured here before the fix:
     `psi_axis` 5.94e-01 against 8.27e-02, `Ip` 1.31e+06 against 2.0e+05.
     `fgsref.py`'s WORKAROUND 1 documents it and predicts exactly that symptom.
  2. `Ip_logic=False` IS NOT THE FORWARD PROBLEM.  The saved arrays are the
     source with no hidden rescaling, so turning the logic off looks right and
     is not: with no current control the forward problem runs away, a bigger
     core carrying more current which grows the core.  Measured: `Ip` 5.69e+05
     against 2.0e+05 and `psi_axis` 3.03e-01.  MEQ controls `I_p`, so this must
     too, and a healthy run then reports `L` within rounding of 1.
  3. A COLD FORWARD START DOES NOT FIND THE DIVERTED BRANCH, AND CURRENT
     CONTROL DOES NOT RESCUE IT.  With `Ip_logic=True` -- so this is not trap 2
     -- a cold forward solve on machine A converges in 87 iterations to
     `psi_axis` 1.0441e-01 against the reference's 8.2718e-02, **2.6e-01
     relative**, with THREE O-points and two X-points and the profile rescaled
     by `L = 0.765` to hold the prescribed current.  That is a multi-lobed
     plasma, not a near miss.  So the seed is part of the problem statement
     exactly as `mkexactguess.py` records it is for MEQ -- and it is worth
     saying plainly that **MEQ's cold-start failures in M-103 are not a MEQ
     weakness**: the reference code, in the same mode, misses the same branch.
     `solve_from_cold()` is the way through and uses the inverse solve as a
     basin finder.

USAGE

    python3 forward.py                       # every case, agreement table
    python3 forward.py A F                   # named cases
    python3 forward.py --write out/          # a standalone input per case
    python3 forward.py --nx 129 257 A        # a refinement ladder

`--write` emits one runnable script per case that reads only its `.npz` and
this directory's two documented workarounds.
"""

import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
os.environ.setdefault("FGSREF_OUT", "/tmp")

# GeneralProfile and picard_loop come from fgsref because both are WORKAROUNDS
# for this fork rather than choices of ours: freegs4e split Jtor into
# Jtor_part1/part2 and left no Profile.Jtor() for its own solve path
# (WORKAROUND 2), and freegs4e.picard.solve has no status return.
from fgsref import CASES, GeneralProfile, picard_loop          # noqa: E402

from freegs4e import equilibrium, machine                       # noqa: E402


def case_for(npz_name):
    """The CASES entry a reference belongs to, by name, shaped or not."""
    stem = npz_name[:-len("_shaped")] if npz_name.endswith("_shaped") \
        else npz_name
    for c in CASES:
        if c["name"] == stem:
            return c
    raise KeyError("no CASES entry for %r" % npz_name)


def build_machine(npz):
    """The reference's conductors, at the reference's currents, as freegs4e
    objects -- and the current vector synced, which is trap 1."""
    labels = [str(x) for x in npz["coil_labels"]]
    coils = []
    for i, label in enumerate(labels):
        kind = str(npz["coil_kind"][i])
        R = float(npz["coil_R"][i])
        Z = float(npz["coil_Z"][i])
        current = float(npz["coil_currents"][i])
        if kind == "shaped":
            from freegs4e.shaped_coil import ShapedCoil
            hw = float(npz["coil_half_width"][i])
            hh = float(npz["coil_half_height"][i])
            shape = [(R - hw, Z - hh), (R + hw, Z - hh),
                     (R + hw, Z + hh), (R - hw, Z + hh)]
            coil = ShapedCoil(shape, current=current)
        else:
            coil = machine.Coil(R, Z, current=current)
        # FORWARD: nothing here is ever re-solved.
        coil.control = False
        coils.append((label, coil))

    wall = machine.Wall(np.asarray(npz["wall_R"], float),
                        np.asarray(npz["wall_Z"], float))
    tok = machine.Machine(coils, wall)
    tok.getCurrentsVec()                    # TRAP 1: without this, no coil field
    return tok


def build_equilibrium(npz, tok, nx=None, seeded=True):
    """The reference's own grid, seeded with its own plasma flux (trap 3).

    `seeded=False` leaves Equilibrium.__init__'s own Gaussian bump, which is
    what solve_from_cold() wants: there the basin is found by the inverse
    stage rather than supplied."""
    R = np.asarray(npz["R"], float)
    Z = np.asarray(npz["Z"], float)
    n = int(nx) if nx else R.size
    eq = equilibrium.Equilibrium(tokamak=tok, order=int(npz["order"]),
                                 Rmin=float(R[0]), Rmax=float(R[-1]),
                                 Zmin=float(Z[0]), Zmax=float(Z[-1]),
                                 nx=n, ny=n)
    # WORKAROUND 3: __init__ leaves psi_func/mask/psi_axis unset.
    eq.mask_outside_limiter = None
    eq.mask_inside_limiter = None

    if not seeded:
        eq._updatePlasmaPsi(eq.plasma_psi)
        return eq

    seed = np.asarray(npz["plasma_psi"], float)
    if seed.shape != (n, n):
        from scipy.interpolate import RectBivariateSpline
        seed = RectBivariateSpline(R, Z, seed)(eq.R[:, 0], eq.Z[0, :])
    eq._updatePlasmaPsi(seed)
    return eq


def build_profile(npz, case):
    """The reference's own source, with the plasma current still PRESCRIBED."""
    return GeneralProfile(Ip=float(npz["Ip"]), fvac=float(npz["fvac"]),
                          psi_n=np.asarray(npz["psi_n"], float),
                          pprime_data=np.asarray(npz["pprime"], float),
                          ffprime_data=np.asarray(npz["ffprime"], float),
                          Raxis=1.0, Ip_logic=True)   # TRAP 2


def solve_from_cold(npz, case, nx=None, rtol=1e-9, maxits=500, gammas=None):
    """COLD to forward, using the INVERSE solve as the basin finder.

    THE PROBLEM THIS SOLVES.  A forward solve from cold does not find this
    branch.  Measured on machine A with the current control ON, so it is not
    trap 2: it converges in 87 iterations to `psi_axis` 1.0441e-01 against the
    reference's 8.2718e-02 -- 2.6e-01 relative -- with THREE O-points and two
    X-points, i.e. a multi-lobed plasma, and the profile rescaled by L = 0.765
    to hold the prescribed current.  That is the same cold-start basin problem
    MEQ has, in the reference code, and it is worth stating plainly: MEQ's
    failures in M-103 are not a MEQ weakness.

    THE WAY THROUGH IS THAT THE INVERSE SOLVE IS A BASIN FINDER.  It converges
    from cold on all seven machines -- that is how the references were made --
    because its geometric targets pin the topology while the currents move.  So
    stage 1 runs it, and stage 2 RELEASES the constraint and re-solves forward
    from that state with the currents frozen at whatever stage 1 chose.

    `gammas` optionally walks the constraint's Tikhonov regularisation up
    before releasing it, which is continuation in how hard the coil circuit is
    allowed to respond: large gamma damps the current adjustment towards doing
    nothing, which is what forward IS.  None releases in one step.
    """
    from fgsref import SyncConstrain

    tok = build_machine(npz)
    eq = build_equilibrium(npz, tok, nx, seeded=False)
    prof = build_profile(npz, case)

    stages = []
    started = time.perf_counter()

    # THE INVERSE STAGE NEEDS COILS IT MAY MOVE.  build_machine() freezes every
    # coil because that is what FORWARD means, so stage 1 has to hand them back
    # -- otherwise the constraint solve has nothing to adjust and this
    # degenerates silently into the cold forward solve it exists to avoid.
    for _, coil in tok.coils:
        coil.control = True
    tok.getCurrentsVec()

    # ---- stage 1: inverse, the basin finder ---------------------------------
    for gamma in (gammas or [1e-12]):
        ctrl = SyncConstrain(xpoints=case["xpoints"], isoflux=case["isoflux"],
                             gamma=gamma)
        st = picard_loop(eq, prof, ctrl, rtol=rtol, atol=1e-14, maxits=maxits)
        stages.append(("inverse gamma=%.1e" % gamma, st["status"],
                       st["iterations"]))

    # ---- stage 2: release, and the currents stay where stage 1 left them ----
    for _, coil in tok.coils:
        coil.control = False
    tok.getCurrentsVec()
    st = picard_loop(eq, prof, None, rtol=rtol, atol=1e-14, maxits=maxits)
    stages.append(("forward", st["status"], st["iterations"]))
    wall = time.perf_counter() - started

    out = dict(nx=eq.R.shape[0], status=st["status"], its=st["iterations"],
               wall=wall, rel=float(st["rel_history"][-1]),
               L=float(getattr(prof, "L", float("nan"))),
               Ip=float(np.sum(prof.jtor)*eq.dR*eq.dZ), stages=stages)
    for key in ("psi_axis", "psi_bndry"):
        v = getattr(eq, key, None)
        out[key] = float(v) if v is not None else float("nan")
    try:
        out["Raxis"] = float(eq.Rmagnetic())
        out["Zaxis"] = float(eq.Zmagnetic())
    except Exception:
        out["Raxis"] = out["Zaxis"] = float("nan")
    out["currents"] = {l: float(c.current) for l, c in tok.coils}
    return out


def trace_forward(npz, case, nx=None, rtol=1e-9, maxits=400, seeded=True,
                  pin_z=None, pin_coils=None, pin_gain=None):
    """The forward solve WITH THE AXIS WATCHED, iteration by iteration.

    WHY A TRACE AND NOT A RESIDUAL HISTORY.  `rel_history` measures how far psi
    moved between passes; it cannot distinguish a solve that is failing to
    resolve from one whose plasma is WALKING VERTICALLY out of the machine.
    Those are different diagnoses with different repairs, and every forward
    failure in the table -- C, D, F, G at `maxits` -- looks the same in a
    residual column.  The axis position separates them in one run.

    `pin_z` optionally applies the minimum vertical feedback that a real
    machine has and this forward problem does not: an up-down ANTISYMMETRIC
    current increment on `pin_coils`, proportional to the axis's departure from
    `pin_z`.  It is deliberately not a constraint solve -- nothing is
    re-optimised, no target flux is imposed, and the equilibrium is still the
    forward one for whatever currents result.  See the note beside it in
    `main()` for what a converged pinned run does and does not establish.
    """
    from freegs4e import critical

    tok = build_machine(npz)
    eq = build_equilibrium(npz, tok, nx, seeded=seeded)
    prof = build_profile(npz, case)

    pairs = _antisymmetric_pairs(tok, pin_coils) if pin_z is not None else []
    base = [float(c.current) for pr in pairs for _, c in pr]
    if pin_z is not None:
        print("    feedback on %d antisymmetric pairs: %s"
              % (len(pairs), ", ".join("%s/%s" % (u[0], l[0])
                                       for u, l in pairs)), flush=True)

    hist = []

    # THE PIN'S CALIBRATION STATE.  `s` is metres of axis per amp of
    # antisymmetric increment: known in advance only if a gain was given, and
    # otherwise measured by one probe.  `cap` keeps a bad first estimate from
    # emptying the machine into one step.
    state = dict(s=(1.0/pin_gain if pin_gain else None), probe=None,
                 applied=0.0, err=0.0,
                 cap=(0.1*max((abs(c.current) for _, c in tok.coils),
                              default=1.0) or 1.0))

    def watch(it, eqn, profiles, arel):
        psi = eqn.psi()
        try:
            opt, xpt = critical.find_critical(eqn.R, eqn.Z, psi,
                                              profiles.mask_inside_limiter,
                                              profiles.Ip)
        except Exception:
            opt, xpt = [], []
        row = dict(it=it, rel=float(arel),
                   L=float(getattr(profiles, "L", float("nan"))),
                   n_o=len(opt), n_x=len(xpt),
                   psi_axis=float(getattr(eqn, "psi_axis", float("nan"))
                                  or float("nan")),
                   psi_bndry=float(getattr(eqn, "psi_bndry", float("nan"))
                                   or float("nan")))
        if len(opt):
            row["Raxis"] = float(opt[0][0])
            row["Zaxis"] = float(opt[0][1])
        else:
            row["Raxis"] = row["Zaxis"] = float("nan")
        row["Ip"] = float(np.sum(profiles.jtor)*eqn.dR*eqn.dZ)
        hist.append(row)

        if pin_z is not None and pairs and np.isfinite(row["Zaxis"]):
            # THE FEEDBACK, APPLIED BETWEEN PASSES.  An antisymmetric increment
            # makes a vertical field and little else, so this moves the axis
            # without retuning the shape -- which is what makes it feedback
            # rather than an inverse solve: no flux is targeted anywhere.
            #
            # A FIXED GAIN DOES NOT TRAVEL BETWEEN MACHINES, MEASURED.  1e+06
            # A/m holds DIII-D at a current increment of 8e-04 A and destroys
            # MAST, whose conductor table resolves solenoid windings and whose
            # currents are two orders larger.  So the default CALIBRATES: one
            # probe increment, the axis's response to it, and a Newton step on
            # the axis height thereafter.  The probe is scaled by the machine's
            # own largest current, which is the only length-free scale here.
            err = row["Zaxis"] - pin_z
            if state["s"] is None:
                if state["probe"] is None:
                    state["probe"] = 1.0e-6*max(abs(c.current)
                                                for _, c in tok.coils) or 1.0
                    step = state["probe"]
                else:
                    ds = err - state["err"]
                    # a probe that moved the axis by less than the drift it was
                    # measured against tells us nothing; probe harder
                    if abs(ds) < 1.0e-3*abs(err) or ds == 0.0:
                        state["probe"] *= 10.0
                        step = state["probe"]
                    else:
                        state["s"] = ds/state["applied"]
                        step = -err/state["s"]
            else:
                step = -err/state["s"]
            step = float(np.clip(step, -state["cap"], state["cap"]))
            for (_, upper), (_, lower) in pairs:
                upper.current += step
                lower.current -= step
            state["err"], state["applied"] = err, step
            tok.getCurrentsVec()
        row["pin"] = float(np.max(np.abs(
            [c.current - i0 for (_, c), i0 in zip(
                [p for pr in pairs for p in pr], base)]))) if pairs else 0.0

    started = time.perf_counter()
    status = picard_loop(eq, prof, None, rtol=rtol, atol=1e-14, maxits=maxits,
                         watch=watch)
    wall = time.perf_counter() - started

    out = dict(nx=eq.R.shape[0], status=status["status"],
               its=status["iterations"], wall=wall,
               rel=float(status["rel_history"][-1]),
               L=float(getattr(prof, "L", float("nan"))),
               Ip=float(np.sum(prof.jtor)*eq.dR*eq.dZ), history=hist)
    for key in ("psi_axis", "psi_bndry"):
        v = getattr(eq, key, None)
        out[key] = float(v) if v is not None else float("nan")
    out["Raxis"] = hist[-1]["Raxis"] if hist else float("nan")
    out["Zaxis"] = hist[-1]["Zaxis"] if hist else float("nan")
    out["currents"] = {l: float(c.current) for l, c in tok.coils}
    return out


def _antisymmetric_pairs(tok, want=None, tol_r=0.05, tol_z=0.05, z_min=0.05):
    """Coils that sit as APPROXIMATE mirror images in Z, as ( upper, lower ).

    Approximate is the operative word and a real machine is why: DIII-D's F1A
    and F1B sit at Z = +0.1683 and -0.1737, so an exact match finds no pairs at
    all and the feedback silently does nothing.  The default tolerances are
    5 cm, which pairs every DIII-D F-coil with its opposite number and pairs
    nothing across banks.

    Every pair carries the same increment, so `pin_gain` is amps of
    antisymmetric current per metre of axis error, shared equally.  That is a
    CHOICE and a crude one -- a machine's vertical control uses a designed
    combination -- and `pin_coils` is how to name a better one.
    """
    coils = [(l, c) for l, c in tok.coils]
    pairs = []
    used = set()
    for i, (li, ci) in enumerate(coils):
        if li in used or (want and li not in want):
            continue
        Ri, Zi = float(getattr(ci, "R", np.nan)), float(getattr(ci, "Z", np.nan))
        if not np.isfinite(Ri) or abs(Zi) < z_min:
            continue
        for lj, cj in coils[i + 1:]:
            if lj in used or (want and lj not in want):
                continue
            Rj = float(getattr(cj, "R", np.nan))
            Zj = float(getattr(cj, "Z", np.nan))
            if abs(Rj - Ri) <= tol_r and abs(Zj + Zi) <= tol_z:
                upper, lower = ((li, ci), (lj, cj)) if Zi > 0 \
                    else ((lj, cj), (li, ci))
                pairs.append((upper, lower))
                used.add(li)
                used.add(lj)
                break
    return pairs


def solve_forward(npz, case, nx=None, rtol=1e-9, maxits=400):
    tok = build_machine(npz)
    eq = build_equilibrium(npz, tok, nx)
    prof = build_profile(npz, case)

    started = time.perf_counter()
    status = picard_loop(eq, prof, None, rtol=rtol, atol=1e-14, maxits=maxits)
    wall = time.perf_counter() - started

    out = dict(nx=eq.R.shape[0], status=status["status"],
               its=status["iterations"], wall=wall,
               rel=float(status["rel_history"][-1]),
               L=float(getattr(prof, "L", float("nan"))),
               Ip=float(np.sum(prof.jtor)*eq.dR*eq.dZ))
    for key, attr in (("psi_axis", "psi_axis"), ("psi_bndry", "psi_bndry")):
        v = getattr(eq, attr, None)
        out[key] = float(v) if v is not None else float("nan")
    try:
        out["Raxis"] = float(eq.Rmagnetic())
        out["Zaxis"] = float(eq.Zmagnetic())
    except Exception:
        out["Raxis"] = out["Zaxis"] = float("nan")
    return out


TEMPLATE = '''"""FORWARD-mode freegs4e input for {name}.

Written by tools/freegs4e-benchmark/forward.py.  Everything below comes out of
{name}.npz, so this describes the same machine MEQ meshes and solves the same
problem MEQ is given: the reference's own coil currents FROZEN, its own profile
arrays, its own grid and wall, the plasma current PRESCRIBED, and NO constraint.

The reference itself is an INVERSE solve, so its currents are an output; this
is the forward re-solve of them, and it must reproduce the reference.  See
CLAUDE_FB.md, *Standing rule: compare in the same mode MEQ solves in*.

    PYTHONPATH=/home/ian/projects/freegs4e python3 {stem}_forward.py [nx]
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, {bench!r})
import forward

NPZ = os.path.join({bench!r}, "{name}.npz")


def main():
    nx = int(sys.argv[1]) if len(sys.argv) > 1 else None
    npz = np.load(NPZ, allow_pickle=True)
    case = forward.case_for("{name}")
    r = forward.solve_forward(npz, case, nx)

    print("  {name}, FORWARD")
    print("    grid        {{:d}}^2".format(r["nx"]))
    print("    status      {{status}} in {{its}} iterations, rel {{rel:.2e}}"
          .format(**r))
    print("    Ip_logic L  {{L:.9f}}   (1 means the saved arrays are the source)"
          .format(**r))
    print("    Ip          {{Ip:.6e}}".format(**r))
    print("    psi_axis    {{psi_axis:.10e}}".format(**r))
    print("    psi_bndry   {{psi_bndry:.10e}}".format(**r))
    print("    axis        ( {{Raxis:.6f}}, {{Zaxis:.6f}} )".format(**r))
    print()
    print("    reference (INVERSE): psi_axis {{:.10e}}  psi_bndry {{:.10e}}"
          .format(float(npz["psi_axis"]), float(npz["psi_bndry"])))
    rel = abs(r["psi_axis"] - float(npz["psi_axis"]))/abs(float(npz["psi_axis"]))
    print("    forward vs inverse in psi_axis: {{:.3e}} relative".format(rel))


if __name__ == "__main__":
    main()
'''


def write_input(name, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    path = os.path.join(out_dir, "%s_forward.py" % name)
    with open(path, "w") as fh:
        fh.write(TEMPLATE.format(name=name, stem=name, bench=HERE))
    return path


def references():
    """Every reference .npz beside this script, filaments and shaped."""
    found = []
    for c in CASES:
        for suffix in ("", "_shaped"):
            name = c["name"] + suffix
            if os.path.exists(os.path.join(HERE, name + ".npz")):
                found.append((name, c))
    return found


def main():
    argv = list(sys.argv[1:])
    out_dir = None
    if "--write" in argv:
        i = argv.index("--write")
        out_dir = argv[i + 1]
        del argv[i:i + 2]
    grids = []
    if "--nx" in argv:
        i = argv.index("--nx")
        j = i + 1
        while j < len(argv) and argv[j].isdigit():
            grids.append(int(argv[j]))
            j += 1
        del argv[i:j]
    cold = "--cold" in argv
    if cold:
        argv.remove("--cold")
    trace = "--trace" in argv
    if trace:
        argv.remove("--trace")
    pin_z = None
    if "--pin-z" in argv:
        i = argv.index("--pin-z")
        pin_z = float(argv[i + 1])
        del argv[i:i + 2]
    pin_gain = None                 # None means CALIBRATE, which is the default
    pin_auto = "--pin" in argv
    if pin_auto:
        argv.remove("--pin")
    pin_coils = None
    if "--pin-coils" in argv:
        i = argv.index("--pin-coils")
        pin_coils = set(argv[i + 1].split(","))
        del argv[i:i + 2]
    if "--pin-gain" in argv:
        i = argv.index("--pin-gain")
        pin_gain = float(argv[i + 1])
        del argv[i:i + 2]
    every = 1
    if "--every" in argv:
        i = argv.index("--every")
        every = int(argv[i + 1])
        del argv[i:i + 2]
    want = [a.upper() for a in argv]

    rows = [(n, c) for n, c in references()
            if not want or any(n.upper().startswith(w) for w in want)]

    if out_dir:
        for name, _ in rows:
            print("  wrote %s" % write_input(name, out_dir), flush=True)
        return

    if trace:
        for name, case in rows:
            npz = np.load(os.path.join(HERE, name + ".npz"), allow_pickle=True)
            for nx in (grids or [None]):
                print("\n  %s   FORWARD, axis watched%s"
                      % (name,
                         "" if pin_z is None
                         else "   PINNED to Zaxis = %.4f, gain %s"
                              % (pin_z, "calibrated" if pin_gain is None
                                 else "%g" % pin_gain)), flush=True)
                r = trace_forward(npz, case, nx, pin_z=pin_z,
                                  pin_coils=pin_coils, pin_gain=pin_gain)
                ref = "    reference (INVERSE): psi_axis %.10e" \
                      % float(npz["psi_axis"])
                if "Raxis" in npz.files:
                    ref += "  axis ( %.6f, %.6f )" % (float(npz["Raxis"]),
                                                      float(npz["Zaxis"]))
                print(ref, flush=True)
                extra = "  %11s" % "max |dI|" if pin_z is not None else ""
                print("    %5s %10s %10s %4s %4s %12s %12s %10s %11s%s"
                      % ("it", "rel", "L", "nO", "nX", "psi_axis",
                         "psi_bndry", "Raxis", "Zaxis", extra), flush=True)
                h = r["history"]
                for row in h:
                    if row["it"] % every and row is not h[-1]:
                        continue
                    extra = ("  %11.4e" % row["pin"]) if pin_z is not None \
                        else ""
                    print("    %5d %10.2e %10.6f %4d %4d %12.6e %12.6e "
                          "%10.6f %11.6f%s"
                          % (row["it"], row["rel"], row["L"], row["n_o"],
                             row["n_x"], row["psi_axis"], row["psi_bndry"],
                             row["Raxis"], row["Zaxis"], extra), flush=True)
                zs = np.array([x["Zaxis"] for x in h])
                good = np.isfinite(zs)
                print("    %s in %d iterations;  Zaxis %+.6f -> %+.6f, "
                      "span %.6f" % (r["status"], r["its"],
                                     zs[good][0] if good.any() else np.nan,
                                     zs[good][-1] if good.any() else np.nan,
                                     (zs[good].max() - zs[good].min())
                                     if good.any() else np.nan), flush=True)
        return

    if cold:
        print("\n  COLD -> INVERSE -> FORWARD, no stored answer used",
              flush=True)
        print("  stage 1 finds the basin with the geometric targets; stage 2\n"
              "  RELEASES them and re-solves forward with the currents frozen\n"
              "  at whatever stage 1 chose.  A forward stage that converges in\n"
              "  one or two iterations is the finding: the inverse solve's\n"
              "  answer IS a fixed point of the forward problem.\n", flush=True)
        print("  %-38s %6s %8s %7s %10s %16s %11s"
              % ("case", "nx", "inverse", "forward", "L", "psi_axis",
                 "vs inverse"), flush=True)
        for name, case in rows:
            npz = np.load(os.path.join(HERE, name + ".npz"), allow_pickle=True)
            try:
                r = solve_from_cold(npz, case, grids[0] if grids else None)
            except Exception as exc:
                print("  %-38s  FAILED %r" % (name, exc), flush=True)
                continue
            ref = float(npz["psi_axis"])
            rel = abs(r["psi_axis"] - ref)/abs(ref)
            inv = "/".join(str(x[2]) for x in r["stages"][:-1])
            print("  %-38s %6d %8s %7d %10.6f %16.10e %11.3e"
                  % (name, r["nx"], inv, r["stages"][-1][2], r["L"],
                     r["psi_axis"], rel), flush=True)
        return

    print("\n  FORWARD AGAINST THE INVERSE REFERENCE IT CAME FROM", flush=True)
    print("  currents frozen, Ip prescribed, no constraint, seeded from the "
          "reference", flush=True)
    if pin_auto:
        print("  VERTICALLY PINNED: one antisymmetric current combination "
              "driven by the axis's own\n  height -- %s -- and NO flux targeted "
              "anywhere.  `max |dI|` is what that\n  cost, and it is the number "
              "to read before reading the agreement."
              % ("gain calibrated from the machine's own response" if
                 pin_gain is None else "%g A per metre of error" % pin_gain),
              flush=True)
    print(flush=True)
    print("  %-38s %6s %8s %6s %10s %16s %11s %11s"
          % ("case", "nx", "status", "its", "L", "psi_axis", "vs inverse",
             "max |dI|"), flush=True)
    for name, case in rows:
        npz = np.load(os.path.join(HERE, name + ".npz"), allow_pickle=True)
        for nx in (grids or [None]):
            try:
                if pin_auto or pin_z is not None:
                    z0 = float(npz["Zaxis"]) if pin_auto else pin_z
                    r = trace_forward(npz, case, nx, pin_z=z0,
                                      pin_coils=pin_coils, pin_gain=pin_gain)
                    dI = r["history"][-1].get("pin", 0.0)
                else:
                    r = solve_forward(npz, case, nx)
                    dI = 0.0
            except Exception as exc:
                print("  %-38s %6s  FAILED %r"
                      % (name, nx or "ref", exc), flush=True)
                continue
            ref = float(npz["psi_axis"])
            rel = abs(r["psi_axis"] - ref)/abs(ref)
            print("  %-38s %6d %8s %6d %10.6f %16.10e %11.3e %11.3e"
                  % (name, r["nx"], r["status"], r["its"], r["L"],
                     r["psi_axis"], rel, dI), flush=True)


if __name__ == "__main__":
    main()
