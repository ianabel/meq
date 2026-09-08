#ifndef MEQ_TESTS_MASCHKEPERRIN_HPP
#define MEQ_TESTS_MASCHKEPERRIN_HPP

/*
 * The Maschke & Perrin rotating equilibrium: the SECOND exact solution of the
 * generalised Grad-Shafranov equation for a plasma in sonic toroidal rotation
 * with an isothermal closure, and the only one in this directory whose
 * temperature and rotation frequency are FUNCTIONS OF psi.
 *
 * Source:
 *   refs/MaschkePerrin.pdf  section 4 -- E. K. Maschke and H. Perrin, "Exact
 *                  solutions of the stationary MHD equations for a rotating
 *                  toroidal plasma", Plasma Physics 22 (1980) 579-594.
 *                  doi 10.1088/0032-1028/22/6/007
 *   refs/Refs.md   which pins the doi, and the 1984 Phys. Lett. A paper that
 *                  is confusable with it
 *   docs/rotation.rst  on the three closures that look alike on the page
 *
 * WHY IT IS HERE AT ALL, GIVEN THAT RotatingSoloviev.hpp EXISTS. The two
 * fixtures solve the SAME partial differential equation: Li & Zhu's (12) and
 * Maschke & Perrin's (4.10) are
 *
 *     Delta*( psi ) = -p1 r^2 exp[ M2 ( r^2/R0^2 - 1 ) ] - F0            Li & Zhu
 *     L F           = -( P/R0^4 ) R^2 exp[ m R^2/2R0^2 ] - M/R0^2        (4.10)
 *
 * which differ by renaming and by a constant absorbed into the amplitude, and
 * the solutions differ no more than that -- (4.17)'s harmonic terms 1, R^2 and
 * X^2 R^2 - R^4/4 are three of Li & Zhu's four, and their particular solutions
 * are the same function. SO NOTHING ABOUT THE DISCRETISATION IS NEW HERE, and a
 * convergence study driven from f() below measures what
 * RotatingSolovievConvergence.cpp already measures.
 *
 * WHAT IS NEW IS THE SOURCE, AND IT IS THE WHOLE POINT OF THE FIXTURE. Li &
 * Zhu's rotating Solov'ev case holds T_0 and Omega_0 CONSTANT, so every check of
 * meq::RotatingSource against a closed form in this tree is a check at
 * T' = omega' = 0. Maschke & Perrin's (4.7) constrains only the RATIO,
 *
 *     omega^2/( Rbar T ) = constant,
 *
 * leaving T( F ) an arbitrary surface function and omega( F ) with it. So this
 * is the one exact solution MEQ has that can be handed to meq::RotatingSource
 * with FIVE profiles varying in psi -- two temperatures, two densities and
 * omega -- whose variations must cancel exactly to leave a source independent of
 * psi. MaschkePerrinConvergence.cpp is where that cancellation is measured.
 *
 * The alternative check is a finite difference of MEQ's own f(), which is what
 * RotatingSourceTests.cpp does, and CLAUDE.md records why that is not the same
 * thing: a term missing from both f() and dFdPsi() is invisible to it. A closed
 * form is not.
 *
 * WHAT IT STILL CANNOT SEE, AND THE RECORD SHOULD NOT BE READ AS SAYING
 * OTHERWISE. (4.7) makes the shared exponent coefficient
 * C = omega^2( Z_1 m_2 - Z_2 m_1 )/( Z_1 T_2 - Z_2 T_1 ) a CONSTANT -- that is
 * precisely what collapses (4.6) to (4.8) and makes the equation solvable in
 * closed form at all. So C'( psi ) is zero here exactly as it is for Li & Zhu,
 * and the C' term -- the one Li & Zhu print with the wrong sign -- is touched by
 * nothing in this tree but RotatingSourceTests' dFdPsi sweep. This fixture does
 * not close that gap and no exact rotating solution can, since a varying C is
 * what makes the equation unsolvable in closed form.
 *
 * dF/dpsi IS IDENTICALLY ZERO, so this sits on Soloviev.hpp's rung of the ladder
 * tests/analytic keeps and not McCarthy.hpp's:
 *
 *   Soloviev.hpp               F constant in psi and in r beyond the r^2
 *   RotatingSoloviev.hpp       F constant in psi, exponential in r^2
 *   MaschkePerrin.hpp          the same, from psi-DEPENDENT profiles
 *   McCarthy.hpp               F linear in psi
 *   ManufacturedNonlinear.hpp  F nonlinear in psi
 *
 * A Newton solve on it must therefore finish in one step, and the fixture says
 * nothing whatever about a Jacobian. RotatingNewtonConvergence.cpp is what does.
 */

#include <cmath>
#include <stdexcept>

namespace meq
{
namespace analytic
{

/*
 * FOUR TRAPS IN READING THIS PAPER, ALL OF WHICH CONVERGE BEAUTIFULLY TO THE
 * WRONG EQUILIBRIUM.
 *
 * 1. THE CITATION. Everyone -- Li & Zhu's [48] included -- cites Maschke &
 *    Perrin, Phys. Lett. A 102 (1984) 106, "An analytic solution ...". That
 *    paper is real and is not this one. This is Plasma Physics 22 (1980) 579,
 *    "Exact solutionS ...", same authors, one word of title apart. refs/Refs.md
 *    pins both dois so the confusable one need not be looked up again.
 *
 * 2. SECTION 3 IS NOT THIS CLOSURE. The paper carries TWO solutions. Section 3
 *    takes the ENTROPY as a surface quantity, p = A( S )rho^gamma, which is a
 *    genuine polytrope and is not (136)'s closure; its g_S carries
 *    ( 1 + Omega^2 R^2/2R0^2 )^{eta+2}, a power law. Section 4 takes the
 *    TEMPERATURE as a surface quantity, B.grad( T ) = 0, and IS the isothermal
 *    closure meq::RotatingSource implements. The two sections' equations are
 *    laid out identically and differ in one exponent, which is the whole
 *    mechanism of the trap.
 *
 * 3. gamma CANCELS, SO EVERY gamma IS USABLE. gamma appears in section 4 only
 *    inside the group gamma Omega^2. Its job there is to convert between the
 *    ADIABATIC sound speed, in which (4.11) defines Omega as the Mach number at
 *    R0, and the ISOTHERMAL one, which is the only one the equation knows:
 *    (4.7) reads omega^2/( Rbar T ) = gamma Omega^2/R0^2, so gamma Omega^2 is
 *    R0^2 omega^2/( Rbar T ), the isothermal Mach number squared at R0. That
 *    group is the single parameter, it is what machSquared() returns below, and
 *    a reader who fixes gamma = 1 to make the fixture "isothermal" has changed
 *    nothing.
 *
 * 4. THE PAPER IS IN j = curl( B ) UNITS. Its (2.13) is
 *    ( L F + JJ' )grad( F ) = -R^2 grad( p ) + rho R^3 omega^2 grad( R ) with no
 *    mu0, against MEQ's F = mu0 r^2 dp/dpsi|_r + g g'. So p_SI = p_M&P/mu0, and
 *    a fixture that carried the paper's p as an SI pressure would be out by
 *    1.26e-6. Everything below is in the paper's units with mu0 = 1, which is
 *    also what the rest of this directory uses.
 *
 * AND pdftotext IS NOT SAFE ON IT. CLAUDE.md records that tool silently dropping
 * minus signs on one paper in refs/ and dropping RADICALS and displacing
 * EXPONENTS on another. On this one it renders "( r1 r2 ) 2" for
 * ( r1 r2 )^{3/2} and loses the square root off an elliptic modulus. Every
 * equation transcribed below was read off the page rendered at 200 dpi.
 */

/*
 * SIGN AND WEIGHT CONVENTION.
 *
 * The paper's operator is defined immediately below its (2.13):
 *
 *     L F := d2F/dX2 + d2F/dR2 - ( 1/R ) dF/dR,
 *
 * with X the vertical coordinate and R the major radius, which is MEQ's Delta*
 * term for term. Its (4.10) is
 *
 *     L F + M/R0^2 + ( P/R0^4 ) R^2 exp[ gamma Omega^2 R^2/2R0^2 ] = 0
 *
 * and MEQ writes the equation as -Delta*( psi ) = F, so
 *
 *     F( r, z, psi ) = ( P/R0^4 ) r^2 exp[ m r^2/2R0^2 ] + M/R0^2            (*)
 *
 * with m := gamma Omega^2. POSITIVE, and F is the full right hand side numerator
 * with no 1/r applied, as everywhere else in this directory. Read against MEQ's
 * own F = mu0 r^2 dp/dpsi|_r + g g' at mu0 = 1, (*) says
 *
 *     dp/dpsi|_r = ( P/R0^4 ) exp[ m r^2/2R0^2 ],      g g' = M/R0^2,
 *
 * and the first of those is the paper's (4.8) with its (4.9) substituted.
 *
 * F IS POSITIVE HERE, SO THE MAGNETIC AXIS IS AN INTERIOR MAXIMUM of psi. That
 * is the opposite of Soloviev.hpp, whose F is single-signed negative and whose
 * axis is an interior MINIMUM -- see CLAUDE.md on why meq::CriticalPointFinder
 * seeds from both nodal extremes. It is the same sign as the high-beta source.
 *
 * deltaStarFD() below recomputes Delta*( psi ) by central differences and
 * MaschkePerrinConvergence.cpp asserts it against -f() over the benchmark box,
 * so the transcription of (3.16), (4.16) and (4.18) is checked rather than
 * trusted.
 */

/**
 * Maschke & Perrin's section 4 equilibrium, (4.17) with (3.16)'s poloidal
 * current term restored:
 *
 *   psi = C P R^2/R0^2 - M X^2/( 2 R0^2 )
 *         + ( eps_a - 1 ) P R^2/( 4 R0^4 ) ( X^2 - R^2/4 )
 *         + F_0
 *         + ( P/m^2 ) [ 1 + m R^2/2R0^2 - exp( m R^2/2R0^2 ) ]
 *
 * The first three terms are Delta*-harmonic but for the M one, which supplies
 * the constant M/R0^2 of (4.10); the last is (4.16)'s g_T, the particular
 * solution of the exponential part. (4.17) itself is printed for M = 0, which
 * the paper adopts in section 3.3 as "the toroidal magnetic field inside the
 * plasma is equal to the vacuum magnetic field"; the M term is (3.16)'s and is
 * carried here because it is what puts g g' under the source check.
 *
 * C IS NOT FREE. Section 4.4 fixes it by requiring an extremum of psi at
 * ( X = 0, R = R_a ), which is what makes R_a the magnetic axis, and (4.18) is
 * the result. This class computes it rather than storing it, so a caller cannot
 * hand it a C inconsistent with its own r_a.
 *
 * Lengths are arbitrary but consistent, and unlike RotatingSoloviev.hpp there is
 * no ln r anywhere, so r = 0 is not excluded by the expansion -- only by the
 * 1/r of the operator itself.
 *
 * NUMERICAL STABILITY AS m -> 0, WHICH IS THE SAME 0/0 RotatingSoloviev.hpp
 * MEETS AND IS WORTH SPELLING OUT AGAIN BECAUSE THE ALGEBRA IS DIFFERENT.
 * (4.16) has a prefactor 1/( gamma^2 Omega^4 ) = 1/m^2 on a brace that vanishes
 * like m^2, so as printed it loses every digit it has for small m. Writing
 *
 *     w := m R^2/( 2 R0^2 ),     G2( w ) := ( e^w - w - 1 )/w^2,
 *
 * the whole of g_T is exactly
 *
 *     ( P/m^2 )( 1 + w - e^w ) = -( P/m^2 ) w^2 G2( w )
 *                              = -( P R^4/( 4 R0^4 ) ) G2( w ),
 *
 * in which m has cancelled ALGEBRAICALLY and nothing is divided by it. G2 is
 * entire with G2( 0 ) = 1/2, so m = 0 needs no branch and falls out as
 * -P R^4/( 8 R0^4 ), which is the limit the paper prints under its own (4.16).
 * The radial derivative goes the same way through G1( w ) := ( e^w - 1 )/w, and
 * so does (4.18)'s second term, which is ( r_a^2/4 ) G1( m r_a^2/2 ).
 *
 * What still needs care is small w, which happens for any m whenever R is small.
 * G1 and G2 switch to their Taylor series below |w| = seriesThreshold().
 */
class MaschkePerrinEquilibrium
{
	public:
		/// @param majorRadiusIn  R0, the radius the rotation is referenced to.
		///                       Strictly positive.
		/// @param machSquaredIn  m = gamma Omega^2 = R0^2 omega^2/( Rbar T ),
		///                       the ISOTHERMAL Mach number squared at R0. See
		///                       trap 3 above on why this is not called Omega^2.
		///                       Zero is legal and gives the static equilibrium.
		/// @param pressureIn     P, the pressure amplitude of (4.9). At mu0 = 1
		///                       the source is ( P/R0^4 ) r^2 exp( ... ).
		/// @param currentIn      M of (4.9), so that g g' = M/R0^2. Zero is
		///                       (4.17) as the paper prints it.
		/// @param ellipticityIn  eps_a of (3.16), which fixes the cross-section
		///                       through (4.19). eps_a = 0 gives a circular
		///                       cross-section at m = 0; it may be negative.
		/// @param axisRadiusIn   r_a = R_a/R0, where the magnetic axis is put.
		///                       Strictly positive.
		/// @param fluxOffsetIn   F_0, the additive constant of (3.16). It is what
		///                       chooses WHICH surface is psi = 0 and nothing
		///                       else: F, Delta*( psi ) and every convergence
		///                       rate are blind to it.
		MaschkePerrinEquilibrium( double majorRadiusIn, double machSquaredIn,
		                          double pressureIn, double currentIn,
		                          double ellipticityIn, double axisRadiusIn,
		                          double fluxOffsetIn )
			: majorRadius( majorRadiusIn ), machSquared( machSquaredIn ),
			  pressure( pressureIn ), current( currentIn ),
			  ellipticity( ellipticityIn ), axisRadius( axisRadiusIn ),
			  fluxOffset( fluxOffsetIn )
		{
			if ( !( majorRadius > 0.0 ) )
			{
				throw std::invalid_argument( "the major radius must be positive" );
			}
			if ( !( axisRadius > 0.0 ) )
			{
				throw std::invalid_argument( "the axis radius must be positive" );
			}
			if ( !( machSquared >= 0.0 ) )
			{
				throw std::invalid_argument( "the isothermal Mach number squared must not be negative" );
			}
		}

		/*
		 * THE THREE CONFIGURATIONS BELOW, AND HOW THEIR CONSTANTS WERE FIXED.
		 *
		 * The paper prints no numbers: its figures are drawn for families of
		 * eps_a and Omega and its section 4 leaves P, M and F_0 to the reader.
		 * So the geometry here is chosen, as RotatingSoloviev.hpp's is, to put a
		 * closed psi = 0 contour well inside MEQ's standard benchmark box
		 * [0.6,1.4] x [-0.6,0.6]:
		 *
		 *     R0 = 1, r_a = 1        the magnetic axis at the geometric centre,
		 *                            which (4.18) then enforces exactly
		 *     eps_a = 0              a circular cross-section in the static
		 *                            limit, so that (4.19)'s elongation is
		 *                            entirely the rotation's doing
		 *     P = 2                  the pressure amplitude
		 *     F_0 from psi( 1.2, 0 ) = 0
		 *
		 * Only ONE geometric condition is available, unlike RotatingSoloviev's
		 * four, and that is a property of the solution family rather than a
		 * choice: eps_a, r_a and m fix the SHAPE of the level sets outright, and
		 * F_0 chooses only which of them is psi = 0. There is no r^2 ln r term
		 * to trade against, so a second condition would over-determine it.
		 *
		 * MEASURED, evaluated in double precision over the benchmark box:
		 *
		 *                              stationary   rotating   withPoloidalCurrent
		 *   C, from (4.18)               0.125000   0.199361      0.199361
		 *   psi at the axis              0.024200   0.061787      0.061787
		 *   min psi on the box          -0.443800  -0.624511     -0.678511
		 *   max psi on its boundary     -0.027000  -0.039187     -0.039187
		 *   ( b/a )^2 at the axis, (4.19) 1.000000   2.297443      1.767263
		 *
		 * The fourth row is the property that matters and the one to re-measure
		 * after any edit: psi is strictly negative on the whole boundary of the
		 * box and positive at the axis, so the psi = 0 level set is a closed
		 * curve strictly inside it and the Dirichlet data is non-homogeneous.
		 *
		 * THE LAST ROW IS ONLY (4.19) FOR THE FIRST TWO COLUMNS. (4.19) is
		 * derived under M = 0, and the M term changes d2psi/dX2 at the axis
		 * without changing d2psi/dR2 -- so with M = 0.3 the true axis elongation
		 * is 1.767 against (4.19)'s 2.297, a 23% disagreement that is the
		 * formula's scope and not an error in either.
		 */

		/// m = 0: no rotation, so g_T collapses to -P R^4/( 8 R0^4 ) and this is
		/// the static polynomial equilibrium the paper notes under (4.16).
		///
		/// It is the control for the whole fixture: anything measured on
		/// rotating() that is not also visible here belongs to the rotation
		/// rather than to the geometry or to the solver. It is also the
		/// configuration in which meq::RotatingSource is handed a null omega, so
		/// it exercises the closure's own omega = 0 branch against a closed form.
		static MaschkePerrinEquilibrium stationary()
		{
			return MaschkePerrinEquilibrium( 1.0, 0.0, 2.0, 0.0, 0.0, 1.0,
			                                 -0.1008 );
		}

		/// m = 1: the representative rotating case, and the one the convergence
		/// study is run on. (4.17) exactly, M being zero.
		///
		/// The isothermal Mach number reaches 1 at R0 and 1.4 at the outboard
		/// edge of the box, so the pressure gradient varies by exp( 0.98 - 0.18 )
		/// = 2.2 across it -- a real exponential rather than a perturbation --
		/// and (4.19) puts the axis elongation at sqrt( 2.297443 ) = 1.5157,
		/// which at eps_a = 0 is produced by the rotation alone.
		static MaschkePerrinEquilibrium rotating()
		{
			return MaschkePerrinEquilibrium( 1.0, 1.0, 2.0, 0.0, 0.0, 1.0,
			                                 -0.16449220852040905 );
		}

		/// m = 1 with M = 0.3, so that g g' is non-zero and (3.16)'s
		/// -M X^2/( 2 R0^2 ) term is under the Delta* check.
		///
		/// It is the only configuration here in which the poloidal current
		/// enters, and it is the one (4.19) does NOT describe -- see the table
		/// above. F_0 is unchanged from rotating(), because the M term vanishes
		/// on the midplane where F_0 is pinned.
		static MaschkePerrinEquilibrium withPoloidalCurrent()
		{
			return MaschkePerrinEquilibrium( 1.0, 1.0, 2.0, 0.3, 0.0, 1.0,
			                                 -0.16449220852040905 );
		}

		/// The poloidal flux function, evaluated in the form that is stable as
		/// m -> 0. See the class comment.
		double psi( double r, double z ) const
		{
			double const r0Sq = majorRadius*majorRadius;
			double const r0Fourth = r0Sq*r0Sq;
			double const r2 = r*r;
			double const z2 = z*z;

			return coefficientC()*pressure*r2/r0Sq
			       - 0.5*current*z2/r0Sq
			       + ( ellipticity - 1.0 )*pressure*r2*( z2 - 0.25*r2 )/( 4.0*r0Fourth )
			       + fluxOffset
			       + particular( r );
		}

		/// grad_bar( psi ) = ( d_r psi, d_z psi ), differentiated by hand rather
		/// than differenced. Not the HDG flux: that is this divided by r, see
		/// flux().
		void gradPsi( double r, double z, double &dPsiDr, double &dPsiDz ) const
		{
			double const r0Sq = majorRadius*majorRadius;
			double const r0Fourth = r0Sq*r0Sq;
			double const r2 = r*r;
			double const z2 = z*z;

			// d/dr of ( eps_a - 1 ) P r^2( z^2 - r^2/4 )/( 4 R0^4 ) is
			// ( eps_a - 1 ) P ( 2 r z^2 - r^3 )/( 4 R0^4 ).
			dPsiDr = 2.0*coefficientC()*pressure*r/r0Sq
			       + ( ellipticity - 1.0 )*pressure*( 2.0*r*z2 - r2*r )/( 4.0*r0Fourth )
			       + particularPrime( r );

			dPsiDz = -current*z/r0Sq
			       + ( ellipticity - 1.0 )*pressure*r2*z/( 2.0*r0Fourth );
		}

		/// The HDG flux q = grad_bar( psi )/r.
		void flux( double r, double z, double &qR, double &qZ ) const
		{
			gradPsi( r, z, qR, qZ );
			qR /= r;
			qZ /= r;
		}

		/// The Grad-Shafranov source of (4.10), read through MEQ's
		/// -Delta*( psi ) = F:
		///
		///     F = ( P/R0^4 ) r^2 exp[ m r^2/( 2 R0^2 ) ] + M/R0^2.
		///
		/// Returns F, not F/r. Independent of psi, which is what makes this
		/// equilibrium linear -- and exponential in r^2, which is the whole
		/// structural consequence of sonic rotation.
		double f( double r, double /*z*/, double /*psiValue*/ ) const
		{
			double const r0Sq = majorRadius*majorRadius;
			double const r0Fourth = r0Sq*r0Sq;

			return pressure*r*r*std::exp( machSquared*r*r/( 2.0*r0Sq ) )/r0Fourth
			       + current/r0Sq;
		}

		/// dF/dpsi. Identically zero, and NOT because the profiles are constant:
		/// T( psi ), omega( psi ) and both densities all vary, and it is (4.7)
		/// plus the linearity of p_T in F that makes their variations cancel.
		/// See the class comment, and theRotatingSourceReproducesTheClosedForm
		/// in MaschkePerrinConvergence.cpp, where the cancellation is measured
		/// through meq::RotatingSource rather than asserted here.
		double dFdPsi( double, double, double ) const
		{
			return 0.0;
		}

		/// Delta*( psi ), by central differences of psi(). Used to verify the
		/// sign convention documented above, and to catch a mistyped term in
		/// (3.16) or (4.16), without depending on the hand-derived gradients.
		///
		/// The implementation RotatingSoloviev.hpp and McCarthy.hpp carry,
		/// deliberately: this is the one check every fixture in this directory
		/// shares, and it is worth being able to diff them and see nothing.
		double deltaStarFD( double r, double z, double h = 1.0e-4 ) const
		{
			// r d_r( ( 1/r ) d_r psi ) as a second difference of the inner
			// quantity, plus d_zz psi.
			auto innerR = [ & ]( double rr )
			{
				return ( psi( rr + h, z ) - psi( rr - h, z ) ) / ( 2.0 * h ) / rr;
			};

			double const dRInner = ( innerR( r + h ) - innerR( r - h ) ) / ( 2.0 * h );
			double const dZZ = ( psi( r, z + h ) - 2.0 * psi( r, z ) + psi( r, z - h ) )
			                 / ( h * h );

			return r * dRInner + dZZ;
		}

		/// C of (4.18), computed rather than stored:
		///
		///     C = ( ( eps_a - 1 )/8 ) r_a^2 + ( 1/2m )( exp( m r_a^2/2 ) - 1 )
		///       = ( ( eps_a - 1 )/8 ) r_a^2 + ( r_a^2/4 ) G1( m r_a^2/2 ),
		///
		/// the second form being the one evaluated, in which m has cancelled and
		/// m = 0 needs no branch. It is what makes r_a the magnetic axis, and
		/// theAxisConditionIsThePapersOwn asserts the condition it comes from.
		double coefficientC() const
		{
			double const raSq = axisRadius*axisRadius;
			return ( ellipticity - 1.0 )*raSq/8.0 + 0.25*raSq*g1( 0.5*machSquared*raSq );
		}

		/// ( b/a )^2 at the magnetic axis from (4.19),
		///
		///     ( b/a )_T^2 = 2 exp( m r_a^2/2 )/( 1 - eps_a ) - 1.
		///
		/// VALID FOR M = 0 ONLY, and this function does not check that, because
		/// the disagreement at M != 0 is itself worth measuring -- see the table
		/// in the class comment. It is the fixture's one independent geometric
		/// check: a published closed form for a quantity that can be recovered
		/// from psi's own Hessian, which is the same service Cerfon &
		/// Freidberg's conditions do for Soloviev.hpp.
		double axisEllipticitySquared() const
		{
			return 2.0*std::exp( 0.5*machSquared*axisRadius*axisRadius )
			       /( 1.0 - ellipticity ) - 1.0;
		}

		/// R0.
		double getMajorRadius() const
		{
			return majorRadius;
		}

		/// m = gamma Omega^2, the isothermal Mach number squared at R0. See trap
		/// 3 in the file comment on why it is not called Omega^2.
		double getMachSquared() const
		{
			return machSquared;
		}

		/// P of (4.9). At mu0 = 1, dp/dpsi|_r is ( P/R0^4 ) exp( m r^2/2R0^2 ).
		double getPressure() const
		{
			return pressure;
		}

		/// M of (4.9), so that g g' = M/R0^2.
		double getCurrent() const
		{
			return current;
		}

		/// eps_a of (3.16).
		double getEllipticity() const
		{
			return ellipticity;
		}

		/// r_a = R_a/R0, the magnetic axis.
		double getAxisRadius() const
		{
			return axisRadius;
		}

		/// F_0 of (3.16), the additive constant that chooses which surface is
		/// psi = 0.
		double getFluxOffset() const
		{
			return fluxOffset;
		}

		/// |w| below which G1 and G2 are summed rather than evaluated, with
		/// w = m r^2/( 2 R0^2 ). Public because it is the one number a test of
		/// the crossover has to know: a sweep that never straddles it measures
		/// one branch twice.
		///
		/// 0.5 costs G2 about one significant digit to cancellation on the
		/// closed-form side and the series needs a dozen terms to reach round-off
		/// on the other, so there is a wide plateau either side and the exact
		/// value is not critical. Both functions switch here even though only G2
		/// needs to, so that the value and the radial derivative change branch at
		/// the same place and ONE continuity test covers both.
		static double seriesThreshold()
		{
			return 0.5;
		}

	private:
		double majorRadius;
		double machSquared;
		double pressure;
		double current;
		double ellipticity;
		double axisRadius;
		double fluxOffset;

		/// ( e^w - 1 )/w, with its removable singularity filled: G1( 0 ) = 1.
		static double g1( double w )
		{
			if ( std::abs( w ) < seriesThreshold() )
			{
				// sum_{k >= 0} w^k/( k + 1 )!, by the ratio w/( k + 1 ). At
				// |w| <= 0.5 the twenty-fifth term is below 1e-40, so the loop
				// bound is a fixed count rather than a tolerance -- a
				// deterministic number of flops, and no branch that could
				// behave differently on two builds.
				double term = 1.0;
				double sum = 1.0;
				for ( int k = 1; k < 25; ++k )
				{
					term *= w/( k + 1 );
					sum += term;
				}
				return sum;
			}
			return std::expm1( w )/w;
		}

		/// ( e^w - w - 1 )/w^2, with its removable singularity filled:
		/// G2( 0 ) = 1/2.
		static double g2( double w )
		{
			if ( std::abs( w ) < seriesThreshold() )
			{
				// sum_{k >= 0} w^k/( k + 2 )!, by the ratio w/( k + 2 ).
				double term = 0.5;
				double sum = 0.5;
				for ( int k = 1; k < 25; ++k )
				{
					term *= w/( k + 2 );
					sum += term;
				}
				return sum;
			}
			return ( std::exp( w ) - w - 1.0 )/( w*w );
		}

		/// (4.16)'s g_T, written as -( P r^4/( 4 R0^4 ) ) G2( w ). m has
		/// cancelled: this is finite and accurate at m = 0, where G2( 0 ) = 1/2
		/// recovers the paper's own -P( r/R0 )^4/8.
		double particular( double r ) const
		{
			double const r0Sq = majorRadius*majorRadius;
			double const w = 0.5*machSquared*r*r/r0Sq;
			return -0.25*pressure*r*r*r*r*g2( w )/( r0Sq*r0Sq );
		}

		/// d/dr of particular(). Differentiating (4.16) gives
		/// -( P r/( m R0^2 ) )( e^w - 1 ), and the same cancellation applies:
		/// e^w - 1 = w G1( w ) = ( m r^2/2R0^2 ) G1( w ), leaving
		/// -( P r^3/( 2 R0^4 ) ) G1( w ).
		double particularPrime( double r ) const
		{
			double const r0Sq = majorRadius*majorRadius;
			double const w = 0.5*machSquared*r*r/r0Sq;
			return -0.5*pressure*r*r*r*g1( w )/( r0Sq*r0Sq );
		}
};

}
}

#endif
