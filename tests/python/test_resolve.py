"""PythonResolve: an Equilibrium solved again as its profiles change (PY-3).

PYTHON-PLAN.md, PY-3: a re-solve is a new meq::Run restarted exactly from the
last state, on the same mesh, and replaces it only if it solves. Asserted:
an unchanged re-solve is the answer in no Newton step; a changed profile
re-solves warm to the cold solve of the same profiles, in fewer steps; a
re-solve that fails leaves the Equilibrium as it was; profiles set before the
first solve are the constructor's; and a configuration with [adaptivity]
solves once. In this process, which imports meq first (meq/__init__.py).
"""

import os

import meq

import netCDF4
import numpy
import pytest

from harness import EXAMPLES, compare, fresh

CASE = "rotating-normalised"
TOML = str(EXAMPLES / f"{CASE}.toml")

# The attributes that describe how a run solved rather than what it solved
# to: a re-solve's are its own.
HOW = {"final_residual", "newton_iterations", "warm_start"}


def omega(value):
    """rotating-normalised's constant Omega, as a table: its slot, in Psi,
    across the range its own density table covers."""
    return meq.Table([-0.5, 1.5], [value, value], [0.0, 0.0], abscissa="Psi")


class inside:
    """Run in `where`, which is where an Equilibrium writes."""

    def __init__(self, where):
        self.where = where

    def __enter__(self):
        self.previous = os.getcwd()
        os.chdir(self.where)
        return self.where

    def __exit__(self, *exc):
        os.chdir(self.previous)


def written(where, attribute):
    with netCDF4.Dataset(where / f"{CASE}.nc") as nc:
        return nc.getncattr(attribute)


def psi(where):
    with netCDF4.Dataset(where / f"{CASE}.nc") as nc:
        nc.set_auto_mask(False)
        return numpy.asarray(nc["psi"][...])[nc["inside"][...] == 1]


def test_an_unchanged_resolve_is_the_answer():
    """The same profiles solved again start at the answer: no more than one
    Newton step, and every field of the file the same, bit for bit."""
    first, again = fresh("unchanged/first"), fresh("unchanged/again")
    eq = meq.Equilibrium(TOML)
    with inside(first):
        eq.solve()
        eq.write()
    with inside(again):
        eq.solve()
        eq.write()
    assert written(again, "newton_iterations") <= 1
    compare(first, again, CASE, skip=HOW)


def test_a_changed_profile_resolves_warm_to_the_cold_solve():
    """Omega 10% up: warm from the last state, the answer is the cold solve's
    of the same profiles to the solver's tolerance, in fewer Newton steps."""
    warm, cold = fresh("changed/warm"), fresh("changed/cold")
    eq = meq.Equilibrium(TOML)
    with inside(fresh("changed/start")):
        eq.solve()
    eq.set_profile("Omega", omega(4.4e5))
    with inside(warm):
        eq.solve()
        eq.write()
    with inside(cold):
        cold_eq = meq.Equilibrium(TOML, profiles={"Omega": omega(4.4e5)})
        cold_eq.solve()
        cold_eq.write()
    gap = numpy.max(numpy.abs(psi(warm) - psi(cold))) / numpy.max(numpy.abs(psi(cold)))
    steps = written(warm, "newton_iterations"), written(cold, "newton_iterations")
    print(f"\n  warm against cold: {gap:.2e} of max |psi|, Newton {steps[0]} against {steps[1]}")
    assert gap < 1e-9
    assert steps[0] < steps[1], "the warm re-solve took no fewer steps than the cold one"


def test_a_failed_resolve_leaves_the_equilibrium_as_it_was():
    """Omega at 3e6 does not solve. The SolveError leaves the Equilibrium
    holding the last state: what it writes is the file a solve of the first
    profiles writes, every attribute and variable, and it solves again
    afterwards from there."""
    reference, after = fresh("failed/reference"), fresh("failed/after")
    with inside(reference):
        once = meq.Equilibrium(TOML)
        once.solve()
        once.write()

    eq = meq.Equilibrium(TOML)
    with inside(fresh("failed/start")):
        eq.solve()
    eq.set_profile("Omega", omega(3.0e6))
    with inside(fresh("failed/attempt")) as attempt:
        with pytest.raises(meq.SolveError):
            eq.solve()
        assert not list(attempt.glob("*.nc")), "the failed re-solve wrote output"
    with inside(after):
        eq.write()
    compare(reference, after, CASE)

    eq.set_profile("Omega", omega(4.4e5))
    with inside(fresh("failed/recovered")):
        eq.solve()
        eq.write()


def test_profiles_set_before_the_first_solve_are_the_constructors():
    """set_profile() before any solve is profiles= at construction, to the
    bit: the run is prepared again, cold."""
    given, set_after = fresh("before/given"), fresh("before/set")
    with inside(given):
        a = meq.Equilibrium(TOML, profiles={"Omega": omega(4.4e5)})
        a.solve()
        a.write()
    with inside(set_after):
        b = meq.Equilibrium(TOML)
        b.set_profile("Omega", omega(4.4e5))
        b.solve()
        b.write()
    compare(given, set_after, CASE)


def test_an_adaptive_configuration_solves_once():
    """miller-adaptive's mesh is the adaptive loop's, not its file's, so a
    re-solve could not restart on it: refused, and the Equilibrium still
    holds and writes its one solve."""
    toml = str(EXAMPLES / "miller-adaptive.toml")
    with inside(fresh("adaptive")) as where:
        eq = meq.Equilibrium(toml)
        eq.solve()
        with pytest.raises(meq.ConfigError):
            eq.solve()
        assert eq.write() is True
        assert list(where.glob("*.nc"))
