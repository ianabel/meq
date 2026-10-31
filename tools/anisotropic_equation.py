#!/usr/bin/env python3
"""AN-0: THE ANISOTROPIC GRAD-SHAFRANOV EQUATION, CHECKED BEFORE ANY OF IT IS BUILT.

ANISOTROPY-PLAN.md, AN-0. refs/AlphaGK.pdf eq (6.18) is the equation, and
refs/Refs.md records its two typos. Each check below exits non-zero if it fails.
None of them reads MEQ's code, because they are about the equation and the
fixture's closed form, not about the solver.

  1. THE MAXWELLIAN BRACKET (6.19), derived rather than transcribed. The
     velocity integral of m w_par^2 dF/dpsi at fixed ( eps, mu ), taken over
     a Maxwellian in the rotating frame, equals the corrected bracket. The
     bracket as printed, without m_s on the centrifugal energy, does not.
  2. THE REDUCTION TO (136). At Delta = 0 and I-bar = g, (6.18)'s right-hand
     side, corrected, negated and in SI, is docs/rotation.rst's F term for
     term. As printed, with R^2 rather than R^4 on the rotation term, it is
     not.
  3. THE OPERATOR. R div( ( 1 - Delta )/R grad psi ) equals
     ( 1 - Delta ) Delta* psi - grad psi . grad Delta for any psi and Delta.
  4. THE SUBTRACTION LOAD. With q = ( 1 - Delta )/R grad( psi_r + psi_c ) and
     q_c = grad( psi_c )/R, the remainder q_r = q - q_c satisfies
     ( R/( 1 - Delta ) ) q_r - grad psi_r = -( Delta/( 1 - Delta ) ) grad psi_c,
     and div q_r = div q wherever Delta* psi_c = 0.
  5. THE MANUFACTURED FIXTURE, tests/analytic/ManufacturedAnisotropic.hpp.
     Its F is the closed form the header states, for both Delta models: the
     symbolic operator applied to the exact flux, and independently a 50-digit
     central difference at sample points. It prints the reference values that
     tests/convergence/AnisotropicConvergence.cpp pins f() to.

Every equation is in MEQ's convention:
    -div( ( 1 - Delta )/R grad psi ) = F/R.

USAGE, from source/:    python3 tools/anisotropic_equation.py
"""

import sys

import mpmath
import sympy as sp

FAILURES = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + (f"  ({detail})" if detail else ""))
    if not ok:
        FAILURES.append(name)


def is_zero(expression):
    return sp.simplify(sp.expand(expression)) == 0


# ---------------------------------------------------------------------------
# 1. (6.19), from a Maxwellian at fixed ( eps, mu )
# ---------------------------------------------------------------------------
def check_maxwellian_bracket():
    print("1. the Maxwellian bracket (6.19)")
    m, T, N, e, Z, R, omega = sp.symbols("m T N e Z R omega", positive=True)
    phi0 = sp.Symbol("phi0", real=True)
    dlnN, dlnT = sp.symbols("dlnN dlnT", real=True)
    w, theta = sp.symbols("w theta", positive=True)

    # The energy in the rotating frame, with the centrifugal potential: RoPP
    # (59)'s eps_s. At fixed eps the only psi-dependence of a Maxwellian is
    # through N( psi ) and T( psi ); phi0 and omega enter through eps alone,
    # which is held, so their derivatives never appear.
    potential = Z*e*phi0 - m*omega**2*R**2/2
    eps = m*w**2/2 + potential
    F = N*(m/(2*sp.pi*T))**sp.Rational(3, 2)*sp.exp(-eps/T)

    # d F/d psi at fixed eps: F [ N'/N - 3/2 T'/T + eps T'/T^2 ]
    dF = F*(dlnN - sp.Rational(3, 2)*dlnT + eps*dlnT/T)

    # d^3 w = 2 pi w^2 sin( theta ) dw dtheta, w_par = w cos( theta )
    integrand = m*(w*sp.cos(theta))**2*dF*2*sp.pi*w**2*sp.sin(theta)
    moment = sp.integrate(sp.integrate(integrand, (theta, 0, sp.pi)), (w, 0, sp.oo))

    n = N*sp.exp(-potential/T)
    corrected = n*(T*dlnN + (Z*e*phi0 - m*omega**2*R**2/2 + T)*dlnT)
    printed = n*(T*dlnN + (Z*e*phi0 - omega**2*R**2/2 + T)*dlnT)

    check("the moment equals the bracket with m_s on 1/2 omega^2 R^2", is_zero(moment - corrected))
    check("the bracket as printed, without m_s, is not the moment", not is_zero(moment - printed))


# ---------------------------------------------------------------------------
# 2. (6.18) at Delta = 0 is docs/rotation.rst's F
# ---------------------------------------------------------------------------
def check_reduction():
    print("2. the reduction of (6.18) to (136)")
    R, mu0, omega, domega, g, dg = sp.symbols("R mu0 omega domega g dg", real=True)
    species = 2
    m = sp.symbols(f"m0:{species}", positive=True)
    Z = sp.symbols(f"Z0:{species}", real=True)
    T = sp.symbols(f"T0:{species}", positive=True)
    N = sp.symbols(f"N0:{species}", positive=True)
    dlnN = sp.symbols(f"dlnN0:{species}", real=True)
    dlnT = sp.symbols(f"dlnT0:{species}", real=True)
    e, phi0 = sp.symbols("e phi0", real=True)

    def density(s):
        return N[s]*sp.exp(m[s]*omega**2*R**2/(2*T[s]) - Z[s]*e*phi0/T[s])

    def bracket(s):
        return density(s)*(T[s]*dlnN[s]
                           + (Z[s]*e*phi0 - m[s]*omega**2*R**2/2 + T[s])*dlnT[s])

    mass_density = sum(m[s]*density(s) for s in range(species))

    # docs/rotation.rst, "The equation, in MEQ's convention"
    rotation_rst = (mu0*R**2*sum(bracket(s) for s in range(species))
                    + mu0*R**4*omega*domega*sum(m[s]*density(s) for s in range(species))
                    + g*dg)

    # (6.18)'s right-hand side in SI (4 pi -> mu0, c drops), at Delta = 0 and
    # I-bar = g. MEQ's F is its NEGATIVE: (1 - Delta) Delta* psi - ... = -F.
    def rhs_618(rotation_power):
        return (-mu0*R**2*sum(bracket(s) for s in range(species))
                - mu0*mass_density*R**rotation_power*omega*domega
                - g*dg/(1 - 0))

    # (6.17)'s rotation current is c m n R^2 omega omega', and (6.15) makes
    # Delta* psi = -( 4 pi R^2/c ) j.grad phi, so the term carries R^2 * R^2.
    check("corrected (6.18), R^4 on the rotation term, is (136)",
          is_zero(-rhs_618(4) - rotation_rst))
    check("(6.18) as printed, R^2 on the rotation term, is not",
          not is_zero(-rhs_618(2) - rotation_rst))


# ---------------------------------------------------------------------------
# 3. the operator identity
# ---------------------------------------------------------------------------
R, z = sp.symbols("R z", positive=True)


def delta_star(f):
    return R*sp.diff(sp.diff(f, R)/R, R) + sp.diff(f, z, 2)


def anisotropic_operator(f, delta):
    """R div_bar( ( 1 - Delta )/R grad_bar f ), the left-hand side of (6.18)."""
    return R*(sp.diff((1 - delta)/R*sp.diff(f, R), R)
              + sp.diff((1 - delta)/R*sp.diff(f, z), z))


def check_operator():
    print("3. the operator")
    psi = sp.Function("psi")(R, z)
    delta = sp.Function("Delta")(R, z)
    paper = (1 - delta)*delta_star(psi) - (sp.diff(psi, R)*sp.diff(delta, R)
                                           + sp.diff(psi, z)*sp.diff(delta, z))
    check("R div( (1 - Delta)/R grad psi ) = (1 - Delta) Delta* psi - grad psi . grad Delta",
          is_zero(anisotropic_operator(psi, delta) - paper))


# ---------------------------------------------------------------------------
# 4. the subtraction load
# ---------------------------------------------------------------------------
def check_subtraction():
    print("4. the subtraction load")
    psi_r = sp.Function("psi_r")(R, z)
    psi_c = sp.Function("psi_c")(R, z)
    delta = sp.Function("Delta")(R, z)

    def grad(f):
        return sp.Matrix([sp.diff(f, R), sp.diff(f, z)])

    def div(v):
        return sp.diff(v[0], R) + sp.diff(v[1], z)

    q = (1 - delta)/R*grad(psi_r + psi_c)
    q_c = grad(psi_c)/R
    q_r = q - q_c

    load = -(delta/(1 - delta))*grad(psi_c)
    residual = (R/(1 - delta))*q_r - grad(psi_r) - load
    check("( R/(1 - Delta) ) q_r - grad psi_r = -( Delta/(1 - Delta) ) grad psi_c",
          all(is_zero(c) for c in residual))

    # div q_c = div( grad psi_c/R ) = Delta*( psi_c )/R, which vanishes for a
    # vacuum field away from its filaments.
    check("div q_r = div q - Delta*( psi_c )/R",
          is_zero(div(q_r) - (div(q) - delta_star(psi_c)/R)))


# ---------------------------------------------------------------------------
# 5. the manufactured fixture
# ---------------------------------------------------------------------------
# Example 5's flux, as tests/analytic/ManufacturedNonlinear.hpp and
# ManufacturedAnisotropic.hpp state it.
R0 = sp.Rational(-1, 2)
KR = sp.Rational(115, 100)*sp.pi
KZ = sp.Rational(115, 100)
PSI = sp.sin(KR*(R + R0))*sp.cos(KZ*z)

# The two Delta models, with the header's parameters.
POSITION_AMPLITUDE = sp.Rational(-3, 10)
POSITION_A = sp.Rational(5, 2)
POSITION_B = sp.Integer(2)
FLUX_AMPLITUDE = sp.Rational(-1, 2)

POINTS = [(sp.Rational(7, 10), sp.Rational(-1, 2)),
          (sp.Integer(1), sp.Integer(0)),
          (sp.Rational(5, 4), sp.Rational(3, 10)),
          (sp.Rational(13, 10), sp.Rational(-1, 5))]


def position_delta():
    return POSITION_AMPLITUDE*sp.cos(POSITION_A*(R - 1))*sp.cos(POSITION_B*z)


def flux_delta(psi_value):
    return FLUX_AMPLITUDE*psi_value**2


def check_fixture():
    print("5. the manufactured fixture")
    s = sp.Symbol("s", real=True)

    f0 = -delta_star(PSI)  # Example 5's F at its own solution
    header_f0 = (KR**2 + KZ**2)*PSI + (KR/R)*sp.cos(KR*(R + R0))*sp.cos(KZ*z)
    check("Example 5's F at Delta = 0 is the header's F_0", is_zero(f0 - header_f0))

    grad_psi_sq = sp.diff(PSI, R)**2 + sp.diff(PSI, z)**2
    delta_p = position_delta()
    models = {
        "position": (delta_p,
                     (1 - delta_p)*header_f0
                     + sp.diff(delta_p, R)*sp.diff(PSI, R) + sp.diff(delta_p, z)*sp.diff(PSI, z)),
        "flux": (flux_delta(PSI),
                 (1 - flux_delta(PSI))*header_f0
                 + sp.diff(flux_delta(s), s).subs(s, PSI)*grad_psi_sq),
    }

    mpmath.mp.dps = 50
    h = mpmath.mpf("1e-12")

    for name, (delta, header_f) in models.items():
        exact_f = -anisotropic_operator(PSI, delta)
        check(f"{name}: -R div( (1 - Delta)/R grad psi ) is the header's closed form",
              is_zero(exact_f - header_f))

        # Independently of sympy's derivatives: a 50-digit central difference
        # of the operator, evaluated through mpmath at each point.
        psi_n = sp.lambdify((R, z), PSI, "mpmath")
        delta_n = sp.lambdify((R, z), delta, "mpmath")
        f_n = sp.lambdify((R, z), header_f, "mpmath")

        def flux_r(r, zz):
            return (1 - delta_n(r, zz))/r*(psi_n(r + h, zz) - psi_n(r - h, zz))/(2*h)

        def flux_z(r, zz):
            return (1 - delta_n(r, zz))/r*(psi_n(r, zz + h) - psi_n(r, zz - h))/(2*h)

        worst = mpmath.mpf(0)
        for (r, zz) in POINTS:
            r, zz = mpmath.mpf(sp.N(r, 60)), mpmath.mpf(sp.N(zz, 60))
            divergence = ((flux_r(r + h, zz) - flux_r(r - h, zz))/(2*h)
                          + (flux_z(r, zz + h) - flux_z(r, zz - h))/(2*h))
            worst = max(worst, abs(-r*divergence - f_n(r, zz)))
        check(f"{name}: 50-digit central difference agrees", worst < mpmath.mpf("1e-20"),
              f"worst {mpmath.nstr(worst, 3)}")

        print(f"     {name}: reference F for AnisotropicConvergence.cpp")
        for (r, zz) in POINTS:
            value = sp.N(header_f.subs({R: r, z: zz}), 20)
            print(f"       {{ {float(r)!r}, {float(zz)!r}, {sp.N(value, 17)} }},")


def main():
    check_maxwellian_bracket()
    check_reduction()
    check_operator()
    check_subtraction()
    check_fixture()
    if FAILURES:
        print(f"\n{len(FAILURES)} FAILED: " + "; ".join(FAILURES))
        return 1
    print("\nall checks hold")
    return 0


if __name__ == "__main__":
    sys.exit(main())
