#ifndef MEQ_TESTS_CURRENTLOOP_HPP
#define MEQ_TESTS_CURRENTLOOP_HPP

/*
 * The exact poloidal flux of a single circular current loop coaxial with the
 * axis, for FB-0 and FB-1.
 *
 * FREE-BOUNDARY-PLAN.md section 3.3 names this field as the independent
 * check on the exterior Dirichlet-to-Neumann map -- "a circular current loop,
 * psi from complete elliptic integrals, which is NOT the formula any of the
 * above came from" -- and section 8's table makes it FB-1's acceptance as
 * well: a vacuum
 * solve driven by coils has the sum of its loops' fields as its exact answer.
 * That is a closed form rather than a self-convergence study, which is what
 * CLAUDE.md's testing stance asks for and what a sign or scaling error can
 * actually be caught by.
 *
 * It is a VACUUM field: Delta* psi = 0 everywhere except on the loop itself, so
 * f() and dFdPsi() are identically zero and a Newton solve on it is affine --
 * one step, like Soloviev.hpp and VacuumHarmonic.hpp. Unlike VacuumHarmonic.hpp
 * it is not a polynomial: it is what a real coil produces, it has the right
 * far-field decay, and it is singular at the conductor, all three of which are
 * properties FB-1 has to survive.
 *
 *
 * THE FORMULA, AND WHAT WAS CHECKED RATHER THAN TRUSTED
 * ----------------------------------------------------
 *
 * For a loop of radius a at height z0 carrying current I, with
 *
 *     d^2 = ( a + r )^2 + ( z - z0 )^2,      k^2 = 4 a r / d^2,
 *
 * the flux is
 *
 *     psi = ( mu0 I / pi )   sqrt( a r ) / k  [ ( 1 - k^2/2 ) K( k ) - E( k ) ]
 *         = ( mu0 I / 2 pi ) d               [ ( 1 - k^2/2 ) K( k ) - E( k ) ].
 *
 * The two lines are the SAME expression: sqrt( a r ) / k = d / 2 identically,
 * straight from the definition of k. Measured over 4732 points, they agree to a
 * worst relative difference of 2.6e-16. The second is what psi() evaluates, and
 * the reason is the axis -- see below. This is the only change made to the
 * transcription; CLAUDE.md records two separate occasions in this project where
 * a transcribed formula converged beautifully to the wrong function, so it
 * was established by substitution rather than by re-reading, and the
 * substitution is reported under "What was measured".
 *
 * K and E are the COMPLETE ELLIPTIC INTEGRALS OF THE MODULUS k, not of the
 * parameter m = k^2. C++17's std::comp_ellint_1 and std::comp_ellint_2 take the
 * modulus, so the argument below is k and not k*k -- which is the classic error
 * with these functions, and is why it is said here twice. Feeding m to them
 * gives a field that is smooth and plausible and is not Delta*-harmonic, and
 * the Delta* check catches it enormously: MEASURED, the slip reads
 * |Delta* psi| = 1.233e+01 against its own |psi( 1, 1 )| = 8.851e-02, a ratio
 * of 1.39e+02 where the correct field reads 4.18e-05. A factor of three
 * million, not a marginal call -- which is the argument for having the check
 * at all rather than proof-reading the argument list.
 *
 *
 * THE CONSTANT, WHICH Delta* CANNOT FIX
 * -------------------------------------
 *
 * Delta* is linear, so Delta* psi = 0 says nothing whatever about the overall
 * multiplicative constant: any multiple of this field passes that test. The
 * constant is therefore stated rather than checked by it, and then checked
 * separately against a quantity that does fix it.
 *
 * The convention is psi = r A_phi -- poloidal flux PER RADIAN, which is the
 * total flux through a circle of radius r divided by 2 pi. It is MEQ's own
 * convention: with q = ( 1/r ) grad_bar psi, CLAUDE.md writes the magnetic
 * field as B = ( -q_z, +q_r ), i.e.
 *
 *     B_r = -( 1/r ) d_z psi,      B_z = +( 1/r ) d_r psi,
 *
 * and in that convention a loop's field at its own centre must be the textbook
 * B_z = mu0 I / ( 2 a ). Measured through flux(), it is, to a worst 2.8e-12
 * relative over a = 0.5, 1.0 and 2.0. THAT is what pins the constant -- and the
 * sign with it -- and mu0 appears explicitly below as 4 pi x 10^-7 rather than
 * being folded into a coefficient.
 *
 * So a, z0 and r are in metres, I is in amperes, and psi is in weber per radian
 * (T m^2). The far field is the loop's dipole,
 * psi -> mu0 I a^2 r^2 / ( 4 R^3 ), which it reaches to 8.25e-7 relative at
 * R = 1000 a, at 8.25e-5 at R = 100 a and 8.23e-3 at R = 10 a. psi > 0
 * everywhere for I > 0.
 *
 * A caller who wants numbers of order one rather than of order 1e-7 should use
 * unitFlux(), which sets the current so that the prefactor mu0 I / pi is
 * exactly 1 -- a loop current of 2.5 MA, which pi/mu0 works out to exactly.
 * That is a choice of units, made in one visible place, and it is the loop
 * every measurement below was taken on.
 *
 *
 * WHERE IT IS VALID, AND WHAT HAPPENS AT THE LOOP
 * -----------------------------------------------
 *
 * k -> 1 as ( r, z ) -> ( a, z0 ), where K( k ) diverges logarithmically.
 * BOTH psi AND ITS GRADIENT DIVERGE THERE -- this is a filamentary current,
 * so it is the two-dimensional line-current logarithm, exactly as
 * A_phi ~ -( mu0 I / 2 pi ) ln( distance ). Measured on the equator with
 * mu0 I / pi = 1, psi at a
 * distance eps outboard of the loop reads 3.4956, 5.7962 and 8.0828 at
 * eps = 1e-3, 1e-5 and 1e-7, against -( 1/2 ) ln eps of 3.4539, 5.7565 and
 * 8.0590: a logarithm, not a finite limit. The gradient goes as 1 / eps,
 * reading -4.9825e+02 and -4.9997e+04 at eps = 1e-3 and 1e-5 against the
 * -1/( 2 eps ) it is approaching.
 *
 * Numerically it gives out well before the mathematics does. 1 - k^2 is formed
 * by cancellation, and on the equator 1 - k ~ eps^2 / 8, so k rounds to exactly
 * 1.0 at eps of about 1e-8 -- where std::comp_ellint_1( 1.0 ) is NaN (not
 * infinity, which is worth knowing: it propagates rather than saturating,
 * while comp_ellint_2( 1.0 ) is a perfectly ordinary 1). The gradient loses its
 * accuracy first, because dA/dk divides by that same cancelled 1 - k^2: at
 * eps = 1e-7 psi is still good to three figures and dPsiDr reads -4.6912e+06
 * against an expected -5.0e+06, six percent out.
 *
 * PRACTICAL DOMAIN: keep evaluation points at least 1e-5 of a loop radius from
 * every loop if the gradient is wanted, 1e-7 if only psi is. For FB-0's
 * semicircle at rho_Gamma = 2.5 or 4, and for FB-1's coils outside the
 * computational domain, that is automatic.
 *
 * There is deliberately NO asymptotic branch handling here:
 * attic/free-boundary/FreeBoundary.cpp has one, it is out of scope for a test
 * fixture, and a fixture that quietly switches formula near the loop is a
 * fixture that can disagree with itself.
 *
 * THE AXIS r = 0 IS THE ONE CASE THAT IS EXACT RATHER THAN MERELY BOUNDED.
 * There k = 0, K( 0 ) = E( 0 ) = pi/2, so the bracket is the difference of
 * two identical doubles and psi( 0, z ) is 0.0 BIT EXACTLY at every z tried
 * -- which is the boundary condition the free-boundary problem imposes on
 * the axis, so it is worth having exactly. That is the reason psi() evaluates
 * the d/2 form: the printed sqrt( a r ) / k form is 0/0 there and returns NaN.
 * Near the axis
 * psi ~ ( mu0 I a^2 / 4 d^3 ) r^2, so it vanishes quadratically -- measured,
 * psi( 1e-3, 0 )/r^2 = 0.785398 against mu0 I a^2 / 4 d^3 = 0.785398 on the
 * unitFlux() loop.
 *
 * THE DERIVATIVES DO NOT ALL FAIL THERE, AND WHICH ONES DO IS WORTH STATING.
 * dPsiDr carries a 1/( 2 r ) and returns NaN at exactly r = 0; dPsiDz does
 * not and returns 0.0, which is also the correct limit, since psi ~ r^2 g( z )
 * makes d_z psi ~ r^2 g'( z ). flux() is NaN in BOTH components, because it
 * divides
 * the pair by r. So gradPsi( 0, z ) is ( NaN, 0 ) and flux( 0, z ) is
 * ( NaN, NaN ): only psi() reaches the axis, and a caller must not read the
 * one component that happens to be finite as evidence that the others are. That
 * is the opposite arrangement to VacuumHarmonic.hpp, whose whole purpose is a
 * flux that is bounded at r = 0, and the two together are the pair FB-A wants.
 *
 *
 * MULTIPLE LOOPS
 * --------------
 *
 * CurrentLoop is ONE loop and holds no vector; a caller with several sums them,
 * because Delta* is linear and a sum of vacuum fields is a vacuum field.
 * CurrentLoopSet below is that sum, written once so that FB-1 -- which drives a
 * vacuum solve with a coil SET -- does not write it again.
 *
 *
 * WHAT WAS MEASURED
 * -----------------
 *
 * On the unitFlux() loop, over r in [0.10, 3.00) and z in [-2.00, 2.00] at a
 * spacing of 0.05, excluding a disc of radius 0.2 about the loop -- 4732 points
 * -- against a representative |psi( 1, 1 )| = 1.9659e-01:
 *
 *  1. Delta* psi = 0, by deltaStarFD():
 *
 *         h = 1e-3   worst 8.504e-04   4.33e-03 of |psi( 1, 1 )|
 *         h = 1e-4   worst 8.218e-06   4.18e-05          "
 *         h = 1e-5   worst 2.967e-04   1.51e-03          "
 *
 *     The worst point is on the excluded disc's own edge -- ( 1.20, 0.00 ) at
 *     h = 1e-3, ( 1.00, 0.20 ) at h = 1e-4 -- which is where the fourth
 *     derivative the stencil truncates is largest. Excluding a disc of radius
 *     0.5 instead, h = 1e-4 gives worst 6.299e-07, 3.20e-06 of the same
 *     reference. h = 1e-5 being WORSE than h = 1e-4 is the difference floor of
 *     the nested stencil, not the field: this is a check on the transcription
 *     at the level a difference can see, and no better.
 *
 *  2. The analytic dPsiDr and dPsiDz against a central difference of psi(), on
 *     the same points, worst relative to |dPsiDr| + |dPsiDz|:
 *
 *         h = 1e-3   9.900e-06   6.859e-06
 *         h = 1e-4   9.898e-08   6.860e-08      clean O( h^2 )
 *         h = 1e-5   1.202e-09   9.308e-10
 *         h = 1e-6   1.470e-08   1.154e-08      round-off, going back up
 *
 *     THAT FLOOR IS THE INSTRUMENT AND NOT THE DERIVATIVE, which is the point
 *     src/meq/Zernike.hpp makes and this file makes again: a central difference
 *     carries its own O( h^2 ) truncation, so the column above measures the
 *     difference until round-off takes over at h = 1e-5. Richardson
 *     extrapolation ( 4 D( h/2 ) - D( h ) ) / 3 from h = 1e-3 sees past it and
 *     reads 3.999e-11 and 2.363e-10 -- one to two orders below anything the
 *     plain difference reaches, and as close to exact as this measurement gets.
 *
 *  3. psi( 0, z ) == 0.0 exactly at z = -2, -0.5, 0, 0.5 and 2. Not to a
 *     tolerance -- the comparison is against the literal 0.0.
 *
 *  4. CurrentLoopSet holding one loop reproduces that loop bit for bit --
 *     0.000e+00 in psi, in both gradient components and in deltaStarFD -- and
 *     a two-loop set ( a = 0.8, z0 = -0.4, I = 3.0 MA and a = 1.6, z0 = 0.7,
 *     I = -1.1 MA ) reproduces the sum of its members' separate fields to
 *     0.000e+00 likewise. That set is itself a vacuum field: Delta* psi worst
 *     7.252e-07 against |psi( 1.2, 0 )| = 2.197e-01, at the same difference
 *     floor as one loop.
 *
 *  5. B_z at the loop centre against the textbook mu0 I / ( 2 a ), taken as the
 *     r -> 0 limit of +q_r and Richardson-extrapolated: 2.8e-12, 1.3e-13 and
 *     2.2e-13 relative at a = 0.5, 1.0 and 2.0. This is the check that fixes
 *     the constant, which item 1 structurally cannot, and it fixes the SIGN
 *     convention with it.
 */

#include <cmath>
#include <cstddef>
#include <vector>

namespace meq
{
namespace analytic
{

/// The poloidal flux of one circular current loop of radius a at height z0,
/// coaxial with the z axis.
///
/// psi() is valid everywhere including r = 0. The derivatives -- dPsiDr(),
/// dPsiDz(), gradPsi() and flux() -- carry a 1/r and require r > 0. Nothing
/// here may be evaluated ON the loop, where the field is genuinely singular;
/// see the file comment for how close is close enough.
class CurrentLoop
{
	public:
		/// pi, to the precision the rest of tests/analytic/ uses.
		static constexpr double pi = 3.14159265358979323846;

		/// The permeability of free space, in H/m. Written out rather than
		/// folded into a prefactor, because the overall constant is exactly
		/// what Delta* psi = 0 cannot check -- see the file comment.
		static constexpr double mu0 = 4.0e-7*pi;

		/// @param radiusIn   the loop radius a, in metres. Strictly positive.
		/// @param heightIn   the loop height z0, in metres.
		/// @param currentIn  the loop current I, in amperes. Signed: reversing
		///                   it reverses psi, since the field is linear in it.
		CurrentLoop( double radiusIn, double heightIn, double currentIn )
			: radiusValue( radiusIn ), heightValue( heightIn ),
			  currentValue( currentIn )
		{
		}

		/// The plain physical case: a one-metre loop on the midplane carrying
		/// one ampere. psi is then of order 1e-7 Wb/rad, which is correct and
		/// inconvenient -- see unitFlux() if the scale matters.
		static CurrentLoop unitCurrent()
		{
			return CurrentLoop( 1.0, 0.0, 1.0 );
		}

		/// The same geometry with the current chosen so that the prefactor
		/// mu0 I / pi is exactly 1, giving psi of order one.
		///
		/// A CHOICE OF UNITS, MADE IN ONE VISIBLE PLACE rather than by scaling
		/// the formula. Every measurement quoted in the file comment was taken
		/// on this loop, so psi( 1, 1 ) = 1.9659e-01 there.
		static CurrentLoop unitFlux()
		{
			return CurrentLoop( 1.0, 0.0, pi/mu0 );
		}

		/// The distance d = sqrt( ( a + r )^2 + ( z - z0 )^2 ) from the field
		/// point to the loop's MIRROR ring at -a, which is what sets the scale
		/// of the elliptic argument. Exposed because psi is d/2 times a bracket
		/// and a reader checking the algebra wants both halves.
		double distance( double r, double z ) const
		{
			double const dz = z - heightValue;
			double const sum = radiusValue + r;
			return std::sqrt( sum*sum + dz*dz );
		}

		/// The elliptic MODULUS k = 2 sqrt( a r ) / d, not the parameter k^2.
		///
		/// Exposed so that a caller can test how close it is standing to the
		/// loop: k -> 1 there, and k rounding to exactly 1.0 is what turns psi
		/// into a NaN. 1 - k of about 1e-16 is the edge.
		double modulus( double r, double z ) const
		{
			return 2.0*std::sqrt( radiusValue*r )/distance( r, z );
		}

		/// The poloidal flux psi = r A_phi, in weber per radian.
		///
		/// Evaluated as ( mu0 I / 2 pi ) d [ ( 1 - k^2/2 ) K - E ], which is
		/// identically the sqrt( a r )/k form and unlike it is exact at r = 0.
		/// std::comp_ellint_1 and _2 take the MODULUS, which is why the
		/// argument is k and the bracket carries k*k.
		double psi( double r, double z ) const
		{
			double const k = modulus( r, z );
			double const kK = std::comp_ellint_1( k );
			double const kE = std::comp_ellint_2( k );

			return ( mu0*currentValue/( 2.0*pi ) )*distance( r, z )
			       *( ( 1.0 - 0.5*k*k )*kK - kE );
		}

		/// d psi / d r. ANALYTIC, not a difference.
		///
		/// With A( k ) := ( 1 - k^2/2 ) K - E and psi = ( mu0 I / 4 pi ) 2 d A,
		/// the chain rule needs dA/dk, and the standard derivatives
		/// dK/dk = ( E - ( 1 - k^2 ) K ) / ( k ( 1 - k^2 ) ) and
		/// dE/dk = ( E - K ) / k collapse it to
		///
		///     dA/dk = ( k / 2 ) [ E / ( 1 - k^2 ) - K ],
		///
		/// which is worth writing down because the three-term form it comes
		/// from does not obviously simplify. The geometry supplies
		/// d d / d r = ( a + r ) / d and
		/// d k / d r = k [ 1/( 2 r ) - ( a + r ) / d^2 ], the latter from
		/// differentiating ln k = ln 2 + ( ln a + ln r )/2 - ln d.
		///
		/// The 1/( 2 r ) is why this is NaN at r = 0. It is a real 1/r and not
		/// an artefact: psi ~ r^2 there, so d psi / d r ~ r and the limit
		/// exists, but the expression as written does not reach it.
		double dPsiDr( double r, double z ) const
		{
			double const dz = z - heightValue;
			double const sum = radiusValue + r;
			double const d2 = sum*sum + dz*dz;
			double const d = std::sqrt( d2 );
			double const k = 2.0*std::sqrt( radiusValue*r )/d;
			double const kK = std::comp_ellint_1( k );
			double const kE = std::comp_ellint_2( k );

			double const bracket = ( 1.0 - 0.5*k*k )*kK - kE;
			double const dAdk = 0.5*k*( kE/( 1.0 - k*k ) - kK );

			double const dDdr = sum/d;
			double const dkdr = k*( 1.0/( 2.0*r ) - sum/d2 );

			return ( mu0*currentValue/( 2.0*pi ) )*( dDdr*bracket + d*dAdk*dkdr );
		}

		/// d psi / d z. Analytic, by the same route, with
		/// d d / d z = ( z - z0 ) / d and d k / d z = -k ( z - z0 ) / d^2.
		///
		/// It has no 1/r of its own, but it shares dAdk's 1/( 1 - k^2 ) and so
		/// gives out at the loop with everything else. On the midplane of a
		/// loop at z0 = 0 it is exactly zero by symmetry.
		double dPsiDz( double r, double z ) const
		{
			double const dz = z - heightValue;
			double const sum = radiusValue + r;
			double const d2 = sum*sum + dz*dz;
			double const d = std::sqrt( d2 );
			double const k = 2.0*std::sqrt( radiusValue*r )/d;
			double const kK = std::comp_ellint_1( k );
			double const kE = std::comp_ellint_2( k );

			double const bracket = ( 1.0 - 0.5*k*k )*kK - kE;
			double const dAdk = 0.5*k*( kE/( 1.0 - k*k ) - kK );

			double const dDdz = dz/d;
			double const dkdz = -k*dz/d2;

			return ( mu0*currentValue/( 2.0*pi ) )*( dDdz*bracket + d*dAdk*dkdz );
		}

		/// grad_bar( psi ) = ( d_r psi, d_z psi ). Not the HDG flux: that is
		/// this divided by r, see flux().
		void gradPsi( double r, double z, double &dR, double &dZ ) const
		{
			dR = dPsiDr( r, z );
			dZ = dPsiDz( r, z );
		}

		/// The HDG flux q = grad_bar( psi ) / r, in the convergence harness's
		/// own signature.
		///
		/// Written as gradPsi() then divided by r, which is what every fixture
		/// here except VacuumHarmonic.hpp does. There is no cancellation to
		/// exploit at r = 0: q_r ~ ( mu0 I a^2 / 2 d^3 ) is a finite limit but
		/// the expression reaching it is not, so this is NaN on the axis. A
		/// caller needing the axis wants VacuumHarmonic.hpp.
		void flux( double r, double z, double &qR, double &qZ ) const
		{
			gradPsi( r, z, qR, qZ );
			qR /= r;
			qZ /= r;
		}

		/// The source. Identically zero: this is a vacuum field, singular only
		/// on the loop, which is not in the computational domain.
		double f( double /*r*/, double /*z*/, double /*psi*/ ) const
		{
			return 0.0;
		}

		/// And so is its derivative, which makes a Newton solve affine -- one
		/// step, exactly, like Soloviev.hpp and VacuumHarmonic.hpp.
		double dFdPsi( double /*r*/, double /*z*/, double /*psi*/ ) const
		{
			return 0.0;
		}

		/// Delta*( psi ), by central differences of psi(), in exactly the
		/// arrangement Soloviev.hpp, ManufacturedNonlinear.hpp and
		/// VacuumHarmonic.hpp use.
		///
		/// It must come out zero, and that is the check that established the
		/// transcription rather than trusting it -- see the file comment for
		/// the numbers and for why h = 1e-5 is worse than h = 1e-4. Keep the
		/// evaluation off the loop by several times h: the stencil reaches
		/// r +/- 2h, and it is the fourth derivative near the singularity that
		/// sets the worst residual.
		double deltaStarFD( double r, double z, double h = 1.0e-4 ) const
		{
			auto innerR = [ & ]( double rr )
			{
				return ( psi( rr + h, z ) - psi( rr - h, z ) )/( 2.0*h )/rr;
			};

			double const dRInner = ( innerR( r + h ) - innerR( r - h ) )/( 2.0*h );
			double const dZZ = ( psi( r, z + h ) - 2.0*psi( r, z ) + psi( r, z - h ) )
			                   /( h*h );

			return r*dRInner + dZZ;
		}

		/// The loop radius a.
		double radius() const { return radiusValue; }

		/// The loop height z0.
		double height() const { return heightValue; }

		/// The loop current I, in amperes.
		double current() const { return currentValue; }

	private:
		double radiusValue;
		double heightValue;
		double currentValue;
};

/// A set of coaxial loops, whose flux is the sum of theirs.
///
/// FB-1 drives a vacuum solve with a coil SET, and the sum of vacuum fields is
/// a vacuum field, so this is arithmetic and nothing more -- it holds no
/// geometry of its own and adds no approximation. It exists so that the sum is
/// written once.
///
/// Every caveat of CurrentLoop applies to every member: psi() reaches the axis
/// and the derivatives do not, and no point may sit on any loop.
class CurrentLoopSet
{
	public:
		/// An empty set. psi() is then zero, which is a vacuum field too.
		CurrentLoopSet() = default;

		/// Add a loop. Order is irrelevant to the answer up to the
		/// associativity of the sum.
		void addLoop( CurrentLoop const &loop )
		{
			loopValues.push_back( loop );
		}

		/// Add a loop by its parameters, for a caller building a set inline.
		void addLoop( double radiusIn, double heightIn, double currentIn )
		{
			loopValues.emplace_back( radiusIn, heightIn, currentIn );
		}

		/// The total flux.
		double psi( double r, double z ) const
		{
			double total = 0.0;
			for ( CurrentLoop const &loop : loopValues )
			{
				total += loop.psi( r, z );
			}
			return total;
		}

		/// d psi / d r of the total field.
		double dPsiDr( double r, double z ) const
		{
			double total = 0.0;
			for ( CurrentLoop const &loop : loopValues )
			{
				total += loop.dPsiDr( r, z );
			}
			return total;
		}

		/// d psi / d z of the total field.
		double dPsiDz( double r, double z ) const
		{
			double total = 0.0;
			for ( CurrentLoop const &loop : loopValues )
			{
				total += loop.dPsiDz( r, z );
			}
			return total;
		}

		/// grad_bar( psi ) of the total field.
		void gradPsi( double r, double z, double &dR, double &dZ ) const
		{
			dR = dPsiDr( r, z );
			dZ = dPsiDz( r, z );
		}

		/// The HDG flux q = grad_bar( psi ) / r of the total field.
		void flux( double r, double z, double &qR, double &qZ ) const
		{
			gradPsi( r, z, qR, qZ );
			qR /= r;
			qZ /= r;
		}

		/// The source of the total field: still identically zero.
		double f( double /*r*/, double /*z*/, double /*psi*/ ) const
		{
			return 0.0;
		}

		/// And its derivative.
		double dFdPsi( double /*r*/, double /*z*/, double /*psi*/ ) const
		{
			return 0.0;
		}

		/// Delta*( psi ) of the total field, summed over the members.
		///
		/// Delta* is linear and so is the difference stencil, so this is the
		/// same quantity as differencing psi() above; it is written as a sum
		/// because that is what makes it obviously zero when each member is.
		double deltaStarFD( double r, double z, double h = 1.0e-4 ) const
		{
			double total = 0.0;
			for ( CurrentLoop const &loop : loopValues )
			{
				total += loop.deltaStarFD( r, z, h );
			}
			return total;
		}

		/// How many loops there are.
		std::size_t size() const { return loopValues.size(); }

		/// The i-th loop.
		CurrentLoop const &loop( std::size_t i ) const { return loopValues[ i ]; }

	private:
		std::vector<CurrentLoop> loopValues;
};

}
}

#endif // MEQ_TESTS_CURRENTLOOP_HPP
