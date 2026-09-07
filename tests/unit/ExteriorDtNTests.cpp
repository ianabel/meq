// Unit tests for the exterior Dirichlet-to-Neumann map in src/meq/ExteriorDtN.
//
// FREE-BOUNDARY-PLAN.md stage FB-0. The class claims that for the axisymmetric
// Grad-Shafranov operator, the exterior of a SEMICIRCLE centred on the axis is a
// DIAGONAL operator on its own trace, in a Gegenbauer basis of order -1/2. That
// claim is the reason free boundary is approachable at all in this tree: it
// replaces a boundary element method -- layer potentials, singular quadrature,
// an O( N^2 ) dense kernel -- with one number per mode.
//
// SO IT HAD BETTER BE TRUE, AND THIS FILE IS WHERE IT IS ESTABLISHED RATHER THAN
// ASSUMED. The tests below are ordered by how independent they are of the
// derivation, weakest first, because that ordering is the whole design:
//
//   1. The class agrees with itself -- projection inverts evaluation, the field
//      reproduces its own trace, the modes decay as advertised. Necessary, and
//      it would pass with a completely wrong basis.
//   2. The basis satisfies the DIFFERENTIAL EQUATION, recomputed by central
//      differences in ( r, z ). This is a real check: Delta* is reassembled in
//      the coordinates the SOLVER uses, not the ( rho, mu ) the separation was
//      done in, so a mistake in the separation shows up here.
//   3. The symbol against a difference of the field it describes.
//   4. THE CURRENT LOOP, which is the one that matters. A field the class knows
//      nothing about -- built from complete elliptic integrals, a formula with
//      no Gegenbauer function anywhere in it -- expanded in the basis and put
//      through the symbol, against its own exact normal derivative.
//
// Only (4) can catch a self-consistent misreading, and CLAUDE.md's standing
// example is why the file is built this way: the Solov'ev coefficients were
// verified against the formula they were derived from, passed, and were wrong.
//
// MFEM-FREE, like the unit under test, so CI can run it. That is not incidental
// -- CI cannot obtain the MFEM branch MEQ needs, so a stage whose acceptance is
// a unit test is a stage CI can actually gate.

#define BOOST_TEST_MODULE ExteriorDtNTests
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>
#include <vector>

#include <boost/math/special_functions/ellint_1.hpp>
#include <boost/math/special_functions/ellint_2.hpp>

#include "meq/ExteriorDtN.hpp"

#include "analytic/CurrentLoop.hpp"
#include "analytic/VacuumHarmonic.hpp"

using meq::ExteriorDtN;

namespace
{
	double const rhoGamma = 2.5;
	int const modes = 8;

	ExteriorDtN standard()
	{
		return ExteriorDtN( 0.0, rhoGamma, modes );
	}

	/// A point on Gamma at polar angle theta from the centre.
	void onGamma( ExteriorDtN const &dtn, double theta, double &r, double &z )
	{
		r = dtn.rhoGamma()*std::sin( theta );
		z = dtn.zCentre() + dtn.rhoGamma()*std::cos( theta );
	}

	/// Delta* psi = d_rr psi - ( 1/r ) d_r psi + d_zz psi, by central
	/// differences, in EXACTLY the arrangement tests/analytic/VacuumHarmonic.hpp
	/// and ManufacturedNonlinear.hpp use. Copied rather than reinvented so that
	/// a disagreement between this file and those is a real disagreement.
	template <typename Psi>
	double deltaStarFD( Psi psi, double r, double z, double h )
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
}

/*
 * The construction contract, and the off-by-TWO that indexing by degree invites.
 *
 * Modes are named by their DEGREE and start at 2. That is not a convention that
 * could have gone the other way: h_n = 2/( n( n-1 )( 2n-1 ) ) divides by
 * n( n-1 ), so n = 0 and n = 1 are singular; and the functions they would name
 * do not vanish on the axis, which is the property the whole basis is chosen
 * for. A caller that loops `for ( int i = 0; i < modeCount(); ++i )` and passes
 * `i` where a degree is wanted is off by two, and two is far enough that the
 * answer is wrong rather than nearly right.
 */
BOOST_AUTO_TEST_CASE( modes_are_indexed_by_degree_from_two_and_refuse_otherwise )
{
	ExteriorDtN const dtn = standard();

	BOOST_TEST( ExteriorDtN::firstMode() == 2 );
	BOOST_TEST( dtn.modeCount() == modes );
	BOOST_TEST( dtn.lastMode() == modes + 1 );

	BOOST_CHECK_THROW( dtn.basis( 0, 1.0, 1.0 ), std::invalid_argument );
	BOOST_CHECK_THROW( dtn.basis( 1, 1.0, 1.0 ), std::invalid_argument );
	BOOST_CHECK_THROW( dtn.mass( 1 ), std::invalid_argument );
	BOOST_CHECK_THROW( dtn.symbol( dtn.lastMode() + 1 ), std::invalid_argument );

	BOOST_CHECK_THROW( ExteriorDtN( 0.0, 0.0, 4 ), std::invalid_argument );
	BOOST_CHECK_THROW( ExteriorDtN( 0.0, -1.0, 4 ), std::invalid_argument );
	BOOST_CHECK_THROW( ExteriorDtN( 0.0, 1.0, 0 ), std::invalid_argument );

	// A coefficient vector of the wrong length is a caller error and not
	// something to pad, since a short one would silently drop the high modes --
	// which are exactly the ones a truncation study is varying.
	std::vector<double> tooFew( static_cast<std::size_t>( modes ) - 1, 0.0 );
	BOOST_CHECK_THROW( dtn.exterior( 5.0, 0.0, tooFew ), std::invalid_argument );
}

/*
 * EVERY MODE VANISHES ON THE AXIS, AND EXACTLY RATHER THAN TO ROUND-OFF.
 *
 * C_n( ±1 ) = 0 for n >= 2 is what lets the flat side of the half-disc go
 * untreated in the exterior: the condition psi = 0 on the axis is satisfied by
 * construction instead of imposed. A basis that met it only to 1e-16 would
 * still be usable, but exactness here is a fact about HOW the class evaluates
 * C_n and is worth pinning.
 *
 * src/meq/ExteriorDtN.cpp uses C_n = ( 1 - mu^2 ) P'_{n-1}( mu )/( n( n-1 ) )
 * rather than the plan's ( P_{n-2} - P_n )/( 2n-1 ). The two are the same
 * function, but the first carries ( 1 - mu^2 ) as an EXPLICIT FACTOR -- so the
 * axis is exactly zero -- while the second is a difference of two quantities
 * that both tend to 1 there, losing about four digits by mu = 0.99999 and more
 * beyond. Since mu = ±1 IS the axis, and FB-A has just measured the axis to
 * cost a power of h in the trace conditioning, that is not a place to spend
 * digits. An exact zero is the observable difference between the two routes.
 */
BOOST_AUTO_TEST_CASE( every_mode_vanishes_exactly_on_the_axis )
{
	ExteriorDtN const dtn = standard();

	double worst = 0.0;
	for ( int n = ExteriorDtN::firstMode(); n <= dtn.lastMode(); ++n )
	{
		// Both ends of the semicircle, and points on the axis beyond it.
		for ( double z : { -rhoGamma, rhoGamma, -10.0, 10.0 } )
			worst = std::max( worst, std::fabs( dtn.basis( n, 0.0, z ) ) );
	}

	std::printf( "\n  C_n on the axis, worst over %d modes : %.3e\n", modes, worst );

	BOOST_TEST( worst == 0.0,
	            "a mode does not vanish exactly on the axis, worst " << worst
	            << ". The ( 1 - mu^2 ) factor is supposed to be explicit in the "
	            "evaluation -- if this reads 1e-16 rather than 0, the "
	            "implementation has reverted to the difference-of-Legendre form "
	            "and has taken its cancellation near the axis with it" );
}

/*
 * The mass, against the closed form, and the closed form against a quadrature.
 *
 * h_n is DERIVED in ExteriorDtN.cpp rather than transcribed from the plan --
 * from the standard integral of ( 1 - mu^2 )[ P'_m ]^2 -- so the closed form and
 * the quadrature agreeing is a check rather than a restatement. The quadrature
 * here is the class's own projection, which is the thing a caller will use.
 */
BOOST_AUTO_TEST_CASE( the_mass_matches_the_closed_form_and_the_basis_is_orthogonal )
{
	ExteriorDtN const dtn = standard();

	std::printf( "\n  h_n against 2/( n( n-1 )( 2n-1 ) )\n" );
	double worstMass = 0.0;
	for ( int n = ExteriorDtN::firstMode(); n <= dtn.lastMode(); ++n )
	{
		double const d = static_cast<double>( n );
		double const expected = 2.0/( d*( d - 1.0 )*( 2.0*d - 1.0 ) );
		worstMass = std::max( worstMass,
		                      std::fabs( dtn.mass( n ) - expected )/expected );
		std::printf( "    n = %2d : %.12f\n", n, dtn.mass( n ) );
	}
	BOOST_TEST( worstMass < 1.0e-15 );

	/*
	 * ORTHOGONALITY, THROUGH THE PROJECTION, AND IT COMES OUT AT ROUND-OFF
	 * BECAUSE THE QUADRATURE IS EXACT HERE.
	 *
	 * The weight dGamma/r = dmu/( 1 - mu^2 ) is singular at the axis, which
	 * would ordinarily mean a special rule. It does not, and that is the second
	 * accident this class rests on: every C_n carries a factor ( 1 - mu^2 ), so
	 * C_n/( 1 - mu^2 ) is a POLYNOMIAL and the singular weight never appears.
	 * A plain Gauss-Legendre rule is therefore exact for these integrals rather
	 * than merely adequate, and the numbers below are what says so -- an
	 * approximate rule would read 1e-8, not 1e-15.
	 */
	double worstSelf = 0.0;
	double worstCross = 0.0;
	for ( int m = ExteriorDtN::firstMode(); m <= dtn.lastMode(); ++m )
	{
		auto trace = [ &dtn, m ]( double r, double z ) { return dtn.basis( m, r, z ); };
		std::vector<double> const a = dtn.coefficients( trace );

		for ( int i = 0; i < dtn.modeCount(); ++i )
		{
			int const n = ExteriorDtN::firstMode() + i;
			double const value = a[ static_cast<std::size_t>( i ) ];
			if ( n == m )
				worstSelf = std::max( worstSelf, std::fabs( value - 1.0 ) );
			else
				worstCross = std::max( worstCross, std::fabs( value ) );
		}
	}

	std::printf( "  projection of mode m: worst |a_m - 1| %.3e, worst off-diagonal %.3e\n",
	             worstSelf, worstCross );

	BOOST_TEST( worstSelf < 1.0e-13,
	            "projecting a mode onto itself gave " << worstSelf
	            << " away from 1. Either the mass is wrong or the quadrature is "
	            "not exact for these integrands -- it should be, since the "
	            "singular weight cancels against C_n's own ( 1 - mu^2 )" );
	BOOST_TEST( worstCross < 1.0e-13,
	            "the basis is not orthogonal in dGamma/r: worst off-diagonal "
	            << worstCross << ". That weight being the one the Grad-Shafranov "
	            "weak form already carries is the whole reason the exterior "
	            "block is diagonal, so this failing collapses stage FB-0" );
}

/*
 * THE MODES SOLVE THE EQUATION, CHECKED IN ( r, z ) RATHER THAN IN ( rho, mu ).
 *
 * This is the first test here that could catch an error in the DERIVATION
 * rather than in the bookkeeping. The separation was done in spherical
 * coordinates and rests on the claim that Delta*'s d_rho term vanishes
 * identically -- so Delta* is reassembled here from scratch by central
 * differences in the cylindrical coordinates the solver actually uses, and
 * applied to the field the class produces.
 *
 * The floor is the instrument's. A central difference of a smooth function at
 * h = 1e-4 carries O( h^2 ) truncation and O( eps/h^2 ) round-off, so a few
 * times 1e-7 relative is as good as this gets; the plan's own run of the same
 * check reported 1e-8 to 8e-6. Do not tighten this expecting exactness --
 * src/meq/Zernike.hpp records the same point about derivatives checked against
 * differences, and Richardson extrapolation is the way past it if a sharper
 * number is ever wanted.
 */
BOOST_AUTO_TEST_CASE( each_exterior_mode_is_delta_star_harmonic )
{
	ExteriorDtN const dtn = standard();

	std::printf( "\n  Delta* of each exterior mode, differenced in ( r, z )\n" );
	std::printf( "    %3s %15s %15s %11s\n", "n", "worst |Delta*|", "|psi| scale", "relative" );

	double worstRelative = 0.0;
	for ( int n = ExteriorDtN::firstMode(); n <= dtn.lastMode(); ++n )
	{
		std::vector<double> a( static_cast<std::size_t>( dtn.modeCount() ), 0.0 );
		a[ static_cast<std::size_t>( n - ExteriorDtN::firstMode() ) ] = 1.0;
		auto psi = [ & ]( double r, double z ) { return dtn.exterior( r, z, a ); };

		double worst = 0.0;
		double scale = 0.0;
		// Well outside Gamma so the stencil never straddles the refusal, and
		// away from the axis so the 1/r in the operator is not the thing being
		// measured.
		for ( double rho = 3.0; rho <= 6.0; rho += 0.5 )
		{
			for ( double theta = 0.3; theta < 2.9; theta += 0.2 )
			{
				double const r = rho*std::sin( theta );
				double const z = dtn.zCentre() + rho*std::cos( theta );
				worst = std::max( worst, std::fabs( deltaStarFD( psi, r, z, 1.0e-4 ) ) );
				scale = std::max( scale, std::fabs( psi( r, z ) ) );
			}
		}

		double const relative = worst/scale;
		worstRelative = std::max( worstRelative, relative );
		std::printf( "    %3d %15.3e %15.3e %11.2e\n", n, worst, scale, relative );
	}

	std::printf( "  worst relative residual: %.3e\n", worstRelative );

	BOOST_TEST( worstRelative < 1.0e-5,
	            "an exterior mode is not Delta*-harmonic: worst relative residual "
	            << worstRelative << ". This is the claim FB-0 rests on, and it is "
	            "checked in ( r, z ) precisely so that an error in the spherical "
	            "separation cannot hide -- so a failure here is the DERIVATION, "
	            "not the bookkeeping" );
}

/*
 * The symbol, against a difference of the very field it claims to describe.
 *
 * ONE-SIDED, and deliberately so: exterior() refuses points inside Gamma, since
 * the modes rho^( 1-n ) grow without bound inward and a centred stencil on Gamma
 * would sample there. A second-order one-sided formula outward is the right
 * instrument and its floor is what the tolerance is set from.
 *
 * WHAT THIS DOES NOT ESTABLISH is the SIGN a weak form wants. symbol() returns
 * the OUTWARD radial derivative, and whether the coupled system wants that or
 * its negative depends on which way its normal points --
 * FREE-BOUNDARY-PLAN.md section 7 says plainly that this will be got wrong at
 * least once and that FB-1, not this file, is what settles it.
 */
BOOST_AUTO_TEST_CASE( the_symbol_is_the_outward_radial_derivative )
{
	ExteriorDtN const dtn = standard();

	double const theta = 1.1;
	double const h = 1.0e-6;
	double worst = 0.0;

	for ( int n = ExteriorDtN::firstMode(); n <= dtn.lastMode(); ++n )
	{
		std::vector<double> a( static_cast<std::size_t>( dtn.modeCount() ), 0.0 );
		a[ static_cast<std::size_t>( n - ExteriorDtN::firstMode() ) ] = 1.0;

		auto at = [ & ]( double rho )
		{
			return dtn.exterior( rho*std::sin( theta ),
			                     dtn.zCentre() + rho*std::cos( theta ), a );
		};
		double const differenced = ( -3.0*at( dtn.rhoGamma() )
		                             + 4.0*at( dtn.rhoGamma() + h )
		                             - at( dtn.rhoGamma() + 2.0*h ) )/( 2.0*h );

		double r = 0.0, z = 0.0;
		onGamma( dtn, theta, r, z );
		double const predicted = dtn.symbol( n )*dtn.basis( n, r, z );

		worst = std::max( worst,
		                  std::fabs( differenced - predicted )/std::fabs( predicted ) );
	}

	std::printf( "\n  the symbol against d psi/d rho, worst relative: %.3e\n", worst );

	BOOST_TEST( worst < 1.0e-6,
	            "the symbol disagrees with the radial derivative of its own "
	            "exterior field by " << worst << " relative" );

	// And the sign, stated as its own assertion because it is the thing most
	// likely to be flipped by a well-meaning edit: every admissible mode DECAYS
	// outward, so the outward derivative is negative and the block entry -- which
	// is minus their product -- is positive.
	for ( int n = ExteriorDtN::firstMode(); n <= dtn.lastMode(); ++n )
	{
		BOOST_TEST( dtn.symbol( n ) < 0.0,
		            "symbol( " << n << " ) is not negative. Every exterior mode "
		            "decays, so its outward radial derivative must be" );
		BOOST_TEST( dtn.blockEntry( n ) > 0.0,
		            "blockEntry( " << n << " ) is not positive" );
	}
}

/*
 * The exterior field reproduces its own trace, and decays at the advertised
 * rate.
 *
 * The weakest test in the file and it is here for one reason: it is the only one
 * that would notice exterior() and basis() drifting apart -- a caller assembles
 * the coupling from basis() and reads the answer back through exterior(), so
 * the two agreeing on Gamma is a contract rather than a triviality.
 */
BOOST_AUTO_TEST_CASE( the_exterior_field_matches_its_trace_and_decays )
{
	ExteriorDtN const dtn = standard();

	std::vector<double> a( static_cast<std::size_t>( dtn.modeCount() ), 0.0 );
	a[ 0 ] = 1.0;
	a[ 2 ] = -0.4;
	a[ 5 ] = 0.15;

	double worstTrace = 0.0;
	for ( double theta = 0.05; theta < 3.10; theta += 0.05 )
	{
		double r = 0.0, z = 0.0;
		onGamma( dtn, theta, r, z );

		double expected = 0.0;
		for ( int i = 0; i < dtn.modeCount(); ++i )
			expected += a[ static_cast<std::size_t>( i ) ]
			            *dtn.basis( ExteriorDtN::firstMode() + i, r, z );

		worstTrace = std::max( worstTrace,
		                       std::fabs( dtn.exterior( r, z, a ) - expected ) );
	}

	std::printf( "\n  exterior() on Gamma against the trace: %.3e\n", worstTrace );
	BOOST_TEST( worstTrace < 1.0e-14 );

	// Doubling rho must scale mode n by exactly 2^( 1 - n ).
	for ( int n : { 2, 5, 9 } )
	{
		std::vector<double> one( static_cast<std::size_t>( dtn.modeCount() ), 0.0 );
		one[ static_cast<std::size_t>( n - ExteriorDtN::firstMode() ) ] = 1.0;

		double const theta = 1.0;
		double const near = dtn.exterior( dtn.rhoGamma()*std::sin( theta ),
		                                  dtn.zCentre() + dtn.rhoGamma()*std::cos( theta ),
		                                  one );
		double const far = dtn.exterior( 2.0*dtn.rhoGamma()*std::sin( theta ),
		                                 dtn.zCentre() + 2.0*dtn.rhoGamma()*std::cos( theta ),
		                                 one );
		double const expected = std::pow( 2.0, 1.0 - static_cast<double>( n ) );
		BOOST_TEST( std::fabs( far/near - expected ) < 1.0e-12,
		            "mode " << n << " decays as " << far/near << " over a doubling "
		            "of rho, where rho^( 1 - n ) demands " << expected );
	}

	// And an interior point is REFUSED rather than extrapolated. At degree 9 a
	// point at half the radius would read 2^8 times its trace value, silently.
	BOOST_CHECK_THROW( dtn.exterior( 0.5*dtn.rhoGamma(), 0.0, a ),
	                   std::invalid_argument );
}

/*
 * ===========================================================================
 * FB-0's ACCEPTANCE, AND THE ONLY TEST HERE THAT CAN CATCH A SELF-CONSISTENT
 * MISREADING.
 * ===========================================================================
 *
 * Everything above checks the class against itself or against the differential
 * equation. Both are necessary and neither is sufficient, for the reason
 * CLAUDE.md records about the Solov'ev coefficients: they were verified against
 * the constraints they were derived from, passed, and were wrong twice. Only an
 * INDEPENDENT quantity catches that.
 *
 * The independent quantity is a circular current loop. Its flux is built from
 * complete elliptic integrals -- `( mu0 I / 2 pi ) d [ ( 1 - k^2/2 ) K( k ) -
 * E( k ) ]` -- a formula with no Legendre polynomial and no Gegenbauer function
 * anywhere in it, whose own correctness is established in
 * tests/analytic/CurrentLoop.hpp by substitution into Delta* and by its field at
 * the loop centre against the textbook mu0 I/( 2a ). It is a vacuum field
 * everywhere off the loop, so a semicircle drawn OUTSIDE the loop encloses no
 * current and the exterior of that semicircle is exactly what ExteriorDtN
 * claims to describe.
 *
 * THE TEST IS THEREFORE: expand the loop's trace on Gamma in the basis, apply
 * the symbol, and compare against the loop's OWN exact normal derivative --
 * which is r q . n from its analytic gradient, computed by a route sharing no
 * line of code with anything above.
 *
 * WAIT. THE LOOP IS INSIDE Gamma, SO ITS FIELD IS THE INTERIOR ONE. Isn't the
 * exterior expansion the wrong object? No, and the reason is the point of the
 * method: outside Gamma the loop's field satisfies Delta* psi = 0 and decays,
 * which is exactly the space the modes rho^( 1 - n ) C_n span. So the loop's
 * field IS an exterior solution out there, its trace on Gamma determines it, and
 * the DtN map must reproduce its normal derivative. That is the whole content
 * of an exact artificial boundary condition, and a loop is the cheapest field
 * that exercises it.
 *
 * MEASURED, AND BOTH COLUMNS BEAT THE PLAN'S OWN RUN OF THE SAME CHECK:
 *
 *     rho_Gamma  modes    trace         DtN ( relative )
 *           2.5     12    5.60e-07      1.95e-05
 *           2.5     24    5.14e-12      3.41e-10
 *           4.0     12    1.15e-09      7.12e-08
 *           4.0     24    1.44e-15      6.89e-14
 *
 * against section 3.3's 1.7e-11 / 2.3e-09 and 4.4e-14 / 6.8e-09.
 *
 * THE DtN COLUMN IS FIVE ORDERS BETTER AT rho_Gamma = 4 AND THE REASON IS THE
 * REFERENCE, NOT THE METHOD. Section 3.3 says of its own figure that it is
 * "limited by the finite-difference reference, not by the method", and it was
 * right: this test compares against the loop's ANALYTIC gradient instead, which
 * carries no truncation, and the floor moves from 6.8e-09 to 6.9e-14. So the
 * plan's number was measuring its instrument and said so; this one is measuring
 * the expansion. tests/analytic/CurrentLoop.hpp's derivatives are themselves
 * checked against a Richardson-extrapolated difference at 4e-11, which is what
 * entitles them to be used as a reference here.
 *
 * What remains IS the series truncation: the trace error is that of a
 * geometrically convergent expansion, and the DtN is the same series
 * differentiated, which multiplies mode n by ( n - 1 ) and so weights the tail
 * more heavily -- which is why the DtN column is consistently the larger of the
 * two. Neither is the class's accuracy; both are the expansion's.
 */
BOOST_AUTO_TEST_CASE( the_symbol_reproduces_a_current_loop_it_knows_nothing_about )
{
	using meq::analytic::CurrentLoop;

	// A loop of radius 1 on the midplane, and a semicircle well outside it.
	// unitFlux() scales mu0 I/pi to one, so the numbers below are O( 1 ) and a
	// relative error is meaningful without a scale factor.
	CurrentLoop const loop = CurrentLoop::unitFlux();

	std::printf( "\n  FB-0 acceptance: the DtN symbol against a current loop\n" );
	std::printf( "    the loop's field comes from elliptic integrals and shares no\n"
	             "    code with the Gegenbauer basis -- this is the independent check\n" );
	std::printf( "    %8s %6s %14s %14s %14s\n",
	             "rhoGamma", "modes", "trace error", "DtN error", "|a_2|" );

	double worstTrace = 0.0;
	double worstDtN = 0.0;

	for ( double radius : { 2.5, 4.0 } )
	{
		for ( int modeCount : { 12, 24 } )
		{
			ExteriorDtN const dtn( 0.0, radius, modeCount );

			// The trace of the loop's field on Gamma, and its expansion.
			auto trace = [ &loop ]( double r, double z ) { return loop.psi( r, z ); };
			std::vector<double> const a = dtn.coefficients( trace );

			double traceError = 0.0;
			double dtnError = 0.0;
			double dtnScale = 0.0;

			// Away from the ends, where the trace is O( 1 ) and a relative
			// comparison means something. The ends are where every mode vanishes,
			// so an absolute error there is uninformative rather than good.
			for ( double theta = 0.25; theta < 2.90; theta += 0.05 )
			{
				double const r = radius*std::sin( theta );
				double const z = radius*std::cos( theta );

				// (a) the expansion reproduces the trace.
				double expanded = 0.0;
				for ( int i = 0; i < dtn.modeCount(); ++i )
					expanded += a[ static_cast<std::size_t>( i ) ]
					            *dtn.basis( ExteriorDtN::firstMode() + i, r, z );
				traceError = std::max( traceError, std::fabs( expanded - loop.psi( r, z ) ) );

				// (b) THE DtN. The symbol applied to the same coefficients,
				// against the loop's own radial derivative -- which is
				// grad psi . rho-hat, assembled here from the ANALYTIC gradient
				// rather than differenced, so the reference carries no truncation
				// of its own.
				double predicted = 0.0;
				for ( int i = 0; i < dtn.modeCount(); ++i )
				{
					int const n = ExteriorDtN::firstMode() + i;
					predicted += a[ static_cast<std::size_t>( i ) ]*dtn.symbol( n )
					             *dtn.basis( n, r, z );
				}

				double dR = 0.0, dZ = 0.0;
				loop.gradPsi( r, z, dR, dZ );
				double const exact = dR*std::sin( theta ) + dZ*std::cos( theta );

				dtnError = std::max( dtnError, std::fabs( predicted - exact ) );
				dtnScale = std::max( dtnScale, std::fabs( exact ) );
			}

			std::printf( "    %8.1f %6d %14.3e %14.3e %14.3e\n",
			             radius, modeCount, traceError, dtnError/dtnScale,
			             std::fabs( a[ 0 ] ) );
			std::fflush( stdout );

			// The best configuration is what the gate is set from; the coarser
			// ones are printed so the convergence in both knobs is visible.
			if ( radius >= 4.0 && modeCount >= 24 )
			{
				worstTrace = traceError;
				worstDtN = dtnError/dtnScale;
			}
		}
	}

	// GATES SET FROM THE MEASUREMENT, at rho_Gamma = 4 with 24 modes: trace
	// 1.44e-15 and DtN 6.89e-14. About three orders of slack each, which is
	// right for quantities this close to round-off -- they depend on the
	// standard library's elliptic integrals and on the Gauss rule, and neither
	// is worth pinning to the last bit.
	BOOST_TEST( worstTrace < 1.0e-12,
	            "the Gegenbauer expansion does not reproduce a current loop's "
	            "trace on Gamma: worst " << worstTrace << ". This is the first "
	            "test in the file that uses a field the basis knows nothing "
	            "about, so a failure here is the BASIS -- everything above would "
	            "still pass with a self-consistently wrong one" );

	BOOST_TEST( worstDtN < 1.0e-11,
	            "the DtN symbol does not reproduce a current loop's normal "
	            "derivative: worst relative " << worstDtN << ". The trace error "
	            "above says whether the expansion is the problem; if the trace is "
	            "fine and this is not, the SYMBOL ( 1 - n )/rhoGamma is wrong -- "
	            "check its sign first, which FREE-BOUNDARY-PLAN.md section 7 warns "
	            "will be got wrong at least once" );
}

/*
 * The spectrum falls geometrically, and the rate is set by rho_plasma/rho_Gamma.
 *
 * THIS IS WHAT MAKES N SMALL, and it is the property the whole method's cost
 * rests on -- an exterior operator needing forty modes would be no cheaper than
 * the kernel it replaces. FREE-BOUNDARY-PLAN.md section 3.3 measures
 * |a_n| falling from 2.0e-1 at n = 2 to 5.8e-8 at n = 20 for rho_Gamma = 2.5,
 * and faster at 4.0, and section 8 lists the rho_Gamma-against-mesh trade as
 * unresolved for a real tokamak geometry.
 *
 * Asserted as a RATIO between two radii rather than as an absolute decay rate,
 * because the rate is a property of the geometry rather than of the class: what
 * the class must get right is that moving Gamma outward makes the spectrum fall
 * faster, since that is the statement a caller uses when choosing rho_Gamma.
 */
BOOST_AUTO_TEST_CASE( the_spectrum_falls_geometrically_and_faster_from_further_out )
{
	using meq::analytic::CurrentLoop;
	CurrentLoop const loop = CurrentLoop::unitFlux();

	std::printf( "\n  the spectrum of a unit loop, by distance of Gamma\n" );
	std::printf( "    %8s %12s %12s %12s %12s\n", "rhoGamma", "|a_2|", "|a_6|",
	             "|a_12|", "ratio 12/2" );

	std::vector<double> falloff;
	for ( double radius : { 2.0, 2.5, 4.0, 8.0 } )
	{
		ExteriorDtN const dtn( 0.0, radius, 12 );
		auto trace = [ &loop ]( double r, double z ) { return loop.psi( r, z ); };
		std::vector<double> const a = dtn.coefficients( trace );

		double const a2 = std::fabs( a[ 0 ] );
		double const a6 = std::fabs( a[ 4 ] );
		double const a12 = std::fabs( a[ 10 ] );
		falloff.push_back( a12/a2 );

		std::printf( "    %8.1f %12.3e %12.3e %12.3e %12.3e\n",
		             radius, a2, a6, a12, a12/a2 );
	}
	std::fflush( stdout );

	// Monotone: further out, the tail is relatively smaller at every step.
	for ( std::size_t i = 1; i < falloff.size(); ++i )
		BOOST_TEST( falloff[ i ] < falloff[ i - 1 ],
		            "moving Gamma outward did not make the spectrum fall faster "
		            "at step " << i << ": " << falloff[ i ] << " against "
		            << falloff[ i - 1 ] << ". That relationship is what lets a "
		            "caller trade rho_Gamma against the mode count, and it is the "
		            "reason N is small" );

	BOOST_TEST( falloff.back() < 1.0e-6,
	            "the spectrum of a loop seen from eight times its radius has a "
	            "twelfth mode only " << falloff.back() << " of its second. If this "
	            "is not small, the expansion is not geometrically convergent and "
	            "the cost argument for this whole approach fails" );
}

/*
 * THE BASIS AGAINST VALUES COMPUTED SOMEWHERE ELSE ENTIRELY.
 *
 * Every number below is an EXACT RATIONAL, produced in sympy in exact rational
 * arithmetic from ( P_{n-2} - P_n )/( 2n - 1 ) -- which is the formula the plan
 * prints and NOT the one src/meq/ExteriorDtN.cpp evaluates. So this table is a
 * check on the identity between the two forms as much as on the arithmetic, and
 * it is independent of both in the way that matters: nothing here came from
 * running the code under test.
 *
 * The sample mu are exact decimals, so every entry is exactly representable as
 * a rational and the literals are those rationals correctly rounded. A double
 * implementation should match to relative 1e-15 or so; anything worse is a real
 * disagreement rather than accumulated round-off.
 *
 * PARITY IS ASSERTED SEPARATELY AND IS FREE: C_n( -mu ) = ( -1 )^n C_n( mu ).
 * It costs one line and it catches a whole class of index error that a table of
 * positive mu alone would not -- an implementation off by one in the Legendre
 * degree gets the parity wrong before it gets the value wrong.
 */
BOOST_AUTO_TEST_CASE( the_basis_matches_independently_computed_exact_values )
{
	// Ten modes so the table's n = 2..10 all exist.
	ExteriorDtN const dtn( 0.0, 2.5, 9 );

	static double const sample[ 7 ] = { -0.9, -0.5, -0.25, 0.0, 0.25, 0.5, 0.9 };

	// C_n( mu ), n = 2 .. 10, from exact rational arithmetic. Spot checks in
	// the comments so a reader can verify one by hand: C_2( -0.9 ) = 19/200,
	// C_2( 0 ) = 1/2, C_10( 0 ) = 7/256.
	static double const expected[ 9 ][ 7 ] =
	{
		{  0.095000000000000000,  0.37500000000000000,   0.46875000000000000,   0.50000000000000000,   0.46875000000000000,   0.37500000000000000,  0.095000000000000000 },
		{ -0.085500000000000000, -0.18750000000000000,  -0.11718750000000000,   0.0,                   0.11718750000000000,   0.18750000000000000,  0.085500000000000000 },
		{  0.072437500000000000,  0.023437500000000000, -0.080566406250000000, -0.12500000000000000,  -0.080566406250000000,  0.023437500000000000, 0.072437500000000000 },
		{ -0.057071250000000000,  0.058593750000000000,  0.075073242187500000,  0.0,                  -0.075073242187500000, -0.058593750000000000, 0.057071250000000000 },
		{  0.040827437500000000, -0.055664062500000000,  0.012130737304687500,  0.062500000000000000,  0.012130737304687500, -0.055664062500000000, 0.040827437500000000 },
		{ -0.025129518750000000,  0.010253906250000000, -0.047664642333984375,  0.0,                   0.047664642333984375, -0.010253906250000000, 0.025129518750000000 },
		{  0.011234772734375000,  0.026458740234375000,  0.011782050132751465, -0.039062500000000000,  0.011782050132751465,  0.026458740234375000, 0.011234772734375000 },
		{ -9.9146601562500000e-05, -0.028884887695312500, 0.026867240667343140, 0.0,                  -0.026867240667343140,  0.028884887695312500, 9.9146601562500000e-05 },
		{ -0.0077126466136718750, 0.0060310363769531250, -0.019666012376546860, 0.027343750000000000, -0.019666012376546860,  0.0060310363769531250, -0.0077126466136718750 }
	};

	double worst = 0.0;
	for ( int n = 2; n <= 10; ++n )
	{
		for ( int j = 0; j < 7; ++j )
		{
			// A point in the direction mu, at any radius: basis() depends on
			// direction alone. Unit radius keeps the geometry obvious.
			double const mu = sample[ j ];
			double const r = std::sqrt( std::max( 0.0, ( 1.0 - mu )*( 1.0 + mu ) ) );
			double const got = dtn.basis( n, r, mu );
			double const want = expected[ n - 2 ][ j ];

			worst = std::max( worst,
			                  std::fabs( got - want )/( 1.0 + std::fabs( want ) ) );
		}
	}

	std::printf( "\n  C_n against exact rationals, n = 2..10, worst relative: %.3e\n",
	             worst );

	BOOST_TEST( worst < 1.0e-14,
	            "the basis disagrees with independently computed exact values by "
	            << worst << ". These came from sympy in exact rational arithmetic "
	            "using the PLAN's formula ( P_{n-2} - P_n )/( 2n-1 ), while the "
	            "implementation evaluates ( 1 - mu^2 ) P'_{n-1}/( n( n-1 ) ) -- so "
	            "this is a check on that identity, not only on the arithmetic" );

	// Parity, which no table of one-signed mu could establish.
	double worstParity = 0.0;
	for ( int n = 2; n <= dtn.lastMode(); ++n )
	{
		for ( double mu : { 0.17, 0.4, 0.83 } )
		{
			double const r = std::sqrt( ( 1.0 - mu )*( 1.0 + mu ) );
			double const plus = dtn.basis( n, r, mu );
			double const minus = dtn.basis( n, r, -mu );
			double const sign = ( n % 2 == 0 ) ? 1.0 : -1.0;
			worstParity = std::max( worstParity, std::fabs( minus - sign*plus ) );
		}
	}
	BOOST_TEST( worstParity < 1.0e-15,
	            "C_n( -mu ) != ( -1 )^n C_n( mu ), worst " << worstParity
	            << ". An implementation off by one in the Legendre degree breaks "
	            "the parity before it breaks the value, so this is the cheaper "
	            "of the two checks" );

	/*
	 * AND THE DERIVATIVE, WHICH IS AN EXACT IDENTITY RATHER THAN A TABLE.
	 *
	 * dC_n/dmu = -P_{n-1}( mu ), verified symbolically for n = 2..12. That is
	 * what src/meq/ExteriorDtN.cpp evaluates, so checking it against a table of
	 * Legendre values would be checking Boost against itself. What is asserted
	 * instead is the ENDPOINT behaviour, which the identity pins exactly and
	 * which is a statement about the axis: C_n'( +1 ) = -1 for every n, and
	 * C_n'( -1 ) = +1 for even n and -1 for odd. So every mode meets the axis at
	 * unit slope in mu -- the derivative form of the simple zero.
	 */
	for ( int n = 2; n <= dtn.lastMode(); ++n )
	{
		BOOST_TEST( std::fabs( dtn.basisDerivative( n, 0.0, 1.0 ) + 1.0 ) < 1.0e-14,
		            "C_" << n << "'( +1 ) is not -1" );
		double const atMinus = dtn.basisDerivative( n, 0.0, -1.0 );
		double const wanted = ( n % 2 == 0 ) ? 1.0 : -1.0;
		BOOST_TEST( std::fabs( atMinus - wanted ) < 1.0e-14,
		            "C_" << n << "'( -1 ) is not " << wanted );
	}
}

/*
 * FB-A's VACUUM FIELDS ARE FB-0's FIRST THREE INTERIOR MODES, AND NEITHER STAGE
 * KNEW IT.
 *
 * tests/analytic/VacuumHarmonic.hpp was written for FB-A -- polynomial
 * Delta*-harmonic functions vanishing on the axis, found by asking what a
 * vacuum solve on a mesh reaching r = 0 could be measured against.
 * src/meq/ExteriorDtN.cpp was written for FB-0, from a separation of variables
 * in spherical coordinates. They turn out to be the same three functions:
 *
 *     rho^2 C_2 = r^2 / 2
 *     rho^3 C_3 = r^2 z / 2
 *     rho^4 C_4 = r^2( 4 z^2 - r^2 ) / 8
 *
 * against VacuumHarmonic's r^2, r^2 z and -( r^4 - 4 r^2 z^2 )/8 -- the same
 * span, mode for mode, with the scalings shown.
 *
 * THAT IS A REAL TIE AND NOT A CURIOSITY. Two parts of FREE-BOUNDARY-PLAN.md
 * written weeks apart, from different arguments, meeting on the same functions:
 * FB-A's fixture is the INTERIOR branch rho^n of exactly the separation FB-0
 * uses the EXTERIOR branch rho^(1-n) of. It is also a free cross-check --
 * each stage's fixture now validates the other's -- and it explains after the
 * fact why the FB-A search found three polynomials and no more: they are the
 * first three admissible Gegenbauer modes, and the fourth candidate that was
 * guessed and failed there was simply not one.
 */
BOOST_AUTO_TEST_CASE( the_interior_modes_are_fb_a_s_vacuum_harmonics )
{
	using meq::analytic::VacuumHarmonic;

	ExteriorDtN const dtn( 0.0, 1.0, 4 );

	// rho^n C_n( mu ), the INTERIOR branch, which exterior() does not provide --
	// it carries rho^( 1 - n ). Assembled here from basis() and the radius.
	auto interior = [ &dtn ]( int n, double r, double z )
	{
		double const rho = std::hypot( r, z );
		return std::pow( rho, static_cast<double>( n ) )*dtn.basis( n, r, z );
	};

	// The three claimed identities, with their scalings.
	VacuumHarmonic const quadratic( 1.0, 0.0, 0.0 );   // r^2
	VacuumHarmonic const linearZ( 0.0, 1.0, 0.0 );     // r^2 z
	VacuumHarmonic const quartic( 0.0, 0.0, 1.0 );     // r^4 - 4 r^2 z^2

	double worst = 0.0;
	for ( double r = 0.1; r < 2.0; r += 0.13 )
	{
		for ( double z = -1.5; z < 1.55; z += 0.17 )
		{
			double const scale = 1.0 + r*r + z*z*r*r;
			worst = std::max( worst,
				std::fabs( interior( 2, r, z ) - 0.5*quadratic.psi( r, z ) )/scale );
			worst = std::max( worst,
				std::fabs( interior( 3, r, z ) - 0.5*linearZ.psi( r, z ) )/scale );
			// rho^4 C_4 = r^2( 4z^2 - r^2 )/8 = -( r^4 - 4 r^2 z^2 )/8.
			worst = std::max( worst,
				std::fabs( interior( 4, r, z ) + 0.125*quartic.psi( r, z ) )/scale );
		}
	}

	std::printf( "\n  the interior modes against FB-A's vacuum harmonics: %.3e\n",
	             worst );

	BOOST_TEST( worst < 1.0e-13,
	            "rho^n C_n does not reproduce tests/analytic/VacuumHarmonic.hpp's "
	            "polynomials, worst " << worst << ". These two fixtures were "
	            "written for different stages from different arguments and are the "
	            "same three functions -- if this fails, one of them has drifted, "
	            "and which one is worth knowing before either is trusted again" );
}

/*
 * ============================================================================
 * FB-0's LAST OPEN ITEM: CEDRES++'S OWN BOUNDARY FORM, IN THIS BASIS
 * ============================================================================
 *
 * FREE-BOUNDARY-PLAN.md section 3.4 states the falsifying test of the whole
 * of section 3, and states it as a challenge rather than a hope: *"assemble
 * CEDRES++ eq (3.5) against the Gegenbauer basis and check it comes out
 * diagonal with the symbol above. If it does not, this section is wrong and
 * the rest of the plan needs the kernel."*
 *
 * This is that test, and it is the most independent check in this file --
 * more so than the current loop, which at least shares the idea of an
 * axisymmetric field. CEDRES++ arrives at the exterior condition by a
 * BOUNDARY INTEGRAL: a single-layer term plus a hypersingular double-layer
 * term with complete elliptic integrals in its kernel, regularised by a
 * double-difference. MEQ arrives at it by SEPARATION OF VARIABLES. The two
 * share the equation and nothing else -- no basis, no measure written the same
 * way, no quadrature, not even the same century of technique.
 *
 * THEIR FORM, transcribed in section 3.5 from refs/CEDRES.pdf page 13 rendered
 * at 900 dpi, because pdftotext deletes the radicals on that page and would
 * turn the elliptic MODULUS into the PARAMETER -- which converges to a wrong
 * answer rather than failing:
 *
 *   c( psi, xi ) = (1/mu0) int_G psi N xi dS
 *                + (1/(2 mu0)) int_G int_G ( psi1 - psi2 ) M ( xi1 - xi2 ) dS1 dS2
 *
 *   M = k / ( 2 pi ( r1 r2 )^( 3/2 ) ) ( ( 2 - k^2 )/( 2 - 2 k^2 ) E( k ) - K( k ) )
 *   N = ( 1/r1 )( 1/d+ + 1/d- - 1/rho )     d± = sqrt( r1^2 + ( rho ± z1 )^2 )
 *   k = sqrt( 4 r1 r2 / ( ( r1 + r2 )^2 + ( z1 - z2 )^2 ) )
 *
 * TWO THINGS MAKE IT COMPUTABLE ON THE SEMICIRCLE, and both are worth having
 * written down because they are what turn a hard quadrature into an easy one.
 *
 * The two distances CLOSE. With r = rho sin t and z = rho cos t,
 *
 *     d+ = 2 rho cos( t/2 ),      d- = 2 rho sin( t/2 )
 *
 * so N is elementary. It carries a 1/t^2 at each pole, and every C_n vanishes
 * there like t^2/2, so the single-layer integrand goes to zero like t^2. The
 * apparent singularity is the basis's to cancel and it does.
 *
 * The DOUBLE-DIFFERENCE IS THE REGULARISATION AND MUST BE KEPT. M is
 * hypersingular, M ~ 1/( pi r d^2 ); each difference is O( d ) and the product
 * cancels it exactly, leaving psi' xi' / ( pi r rho^2 ) on the diagonal.
 * Assembling int int C_m M C_n directly instead would diverge. What survives at
 * next order is a Delta^2 log Delta, so splitting the inner integral AT the
 * diagonal and putting Gauss on each half converges properly -- a tensor rule
 * straddling it would not.
 *
 * MEASURED: the off-diagonal entries come out at 1e-15 to 1e-10 against
 * diagonals of order 1e-1, and the diagonal agrees with blockEntry( n ) to TEN
 * significant figures. Section 3 stands, and the plan does not need the kernel.
 */
namespace
{
	/// Gauss-Legendre nodes and weights on [ -1, 1 ], by Newton on P_n.
	/// Computed ONCE and mapped affinely: the inner rule is rebuilt on a
	/// different interval for every outer point, and recomputing the nodes
	/// there costs more than every kernel evaluation put together.
	void gaussLegendreReference( int n, std::vector<double> &x,
	                             std::vector<double> &w )
	{
		x.assign( n, 0.0 );
		w.assign( n, 0.0 );
		for ( int i = 0; i < n; ++i )
		{
			double t = std::cos( M_PI*( i + 0.75 )/( n + 0.5 ) );
			double p0 = 1.0, p1 = 0.0, dp = 1.0;
			for ( int iteration = 0; iteration < 100; ++iteration )
			{
				p0 = 1.0;
				p1 = 0.0;
				for ( int j = 0; j < n; ++j )
				{
					double const p2 = p1;
					p1 = p0;
					p0 = ( ( 2.0*j + 1.0 )*t*p1 - j*p2 )/( j + 1.0 );
				}
				dp = n*( t*p0 - p1 )/( t*t - 1.0 );
				double const step = -p0/dp;
				t += step;
				if ( std::fabs( step ) < 1.0e-15 )
					break;
			}
			x[ i ] = t;
			w[ i ] = 2.0/( ( 1.0 - t*t )*dp*dp );
		}
	}
}

BOOST_AUTO_TEST_CASE( cedres_boundary_form_is_diagonal_in_this_basis )
{
	double const rhoGamma = 1.3;
	int const modes = 4;
	int const outerPoints = 60, innerPoints = 60;

	meq::ExteriorDtN const dtn( 0.0, rhoGamma, modes );
	int const first = meq::ExteriorDtN::firstMode();

	auto onGamma = [ & ]( double t, double &r, double &z )
	{
		r = rhoGamma*std::sin( t );
		z = rhoGamma*std::cos( t );
	};

	auto mode = [ & ]( int n, double t )
	{
		double r = 0.0, z = 0.0;
		onGamma( t, r, z );
		return dtn.basis( n, r, z );
	};

	// N, with the two distances closed on the semicircle.
	auto singleLayer = [ & ]( double t )
	{
		double const half = 0.5*t;
		return ( 0.5/std::cos( half ) + 0.5/std::sin( half ) - 1.0 )
		       /( rhoGamma*rhoGamma*std::sin( t ) );
	};

	// M, with k the MODULUS. Boost's ellint_1 and ellint_2 take the modulus,
	// which is the convention the paper prints and the one pdftotext destroys.
	auto doubleLayer = [ & ]( double t1, double t2 )
	{
		double r1 = 0.0, z1 = 0.0, r2 = 0.0, z2 = 0.0;
		onGamma( t1, r1, z1 );
		onGamma( t2, r2, z2 );
		double const sum = r1 + r2, gap = z1 - z2;
		double const k2 = 4.0*r1*r2/( sum*sum + gap*gap );
		double const k = std::sqrt( k2 );
		double const bracket = ( 2.0 - k2 )/( 2.0 - 2.0*k2 )
		                       *boost::math::ellint_2( k )
		                     - boost::math::ellint_1( k );
		return k/( 2.0*M_PI*std::pow( r1*r2, 1.5 ) )*bracket;
	};

	std::vector<double> reference, referenceWeights;
	gaussLegendreReference( std::max( outerPoints, innerPoints ),
	                        reference, referenceWeights );

	std::vector<double> outerT( outerPoints ), outerW( outerPoints );
	for ( int i = 0; i < outerPoints; ++i )
	{
		outerT[ i ] = 0.5*M_PI*( 1.0 + reference[ i ] );
		outerW[ i ] = 0.5*M_PI*referenceWeights[ i ];
	}

	std::printf( "\n  CEDRES++ eq (3.5) assembled in the Gegenbauer basis, "
	             "rho_Gamma = %.2f, %d x %d Gauss\n", rhoGamma,
	             outerPoints, innerPoints );
	std::printf( "     m   n     mu0 c( C_m, C_n )     blockEntry( n )"
	             "      relative\n" );

	double worstOffDiagonal = 0.0, worstDiagonal = 0.0, scale = 0.0;

	for ( int m = first; m < first + modes; ++m )
		for ( int n = first; n < first + modes; ++n )
		{
			double single = 0.0;
			for ( int i = 0; i < outerPoints; ++i )
				single += outerW[ i ]*mode( m, outerT[ i ] )
				          *singleLayer( outerT[ i ] )*mode( n, outerT[ i ] )
				          *rhoGamma;

			double doubled = 0.0;
			for ( int i = 0; i < outerPoints; ++i )
			{
				double const t1 = outerT[ i ];
				double const cm1 = mode( m, t1 ), cn1 = mode( n, t1 );
				double inner = 0.0;
				// SPLIT AT THE DIAGONAL. What is left there is Delta^2 log
				// Delta, which Gauss handles on each side and would not
				// handle straddling.
				double const ends[ 3 ] = { 0.0, t1, M_PI };
				for ( int piece = 0; piece < 2; ++piece )
				{
					double const a = ends[ piece ], b = ends[ piece + 1 ];
					if ( b - a < 1.0e-14 )
						continue;
					for ( int j = 0; j < innerPoints; ++j )
					{
						double const t2 = 0.5*( a + b )
						                + 0.5*( b - a )*reference[ j ];
						double const weight = 0.5*( b - a )*referenceWeights[ j ];
						inner += weight*( cm1 - mode( m, t2 ) )
						         *doubleLayer( t1, t2 )
						         *( cn1 - mode( n, t2 ) )*rhoGamma;
					}
				}
				doubled += 0.5*outerW[ i ]*inner*rhoGamma;
			}

			double const assembled = single + doubled;
			double const expected = dtn.blockEntry( n );
			scale = std::max( scale, expected );

			if ( m == n )
			{
				double const relative = std::fabs( assembled - expected )/expected;
				worstDiagonal = std::max( worstDiagonal, relative );
				std::printf( "  %4d %3d  %18.10e  %18.10e   %10.2e\n",
				             m, n, assembled, expected, relative );
			}
			else
			{
				worstOffDiagonal = std::max( worstOffDiagonal,
				                             std::fabs( assembled ) );
				std::printf( "  %4d %3d  %18.10e  %18s   %10s\n",
				             m, n, assembled, "0", "-" );
			}
		}

	std::printf( "    worst diagonal %.2e relative, worst off-diagonal %.2e "
	             "against a scale of %.2e\n",
	             worstDiagonal, worstOffDiagonal, scale );

	BOOST_TEST( worstOffDiagonal < 1.0e-6*scale,
	            "CEDRES++'s boundary form is NOT diagonal in the Gegenbauer "
	            "basis: the worst off-diagonal is " << worstOffDiagonal
	            << " against a diagonal scale of " << scale
	            << ". FREE-BOUNDARY-PLAN.md section 3 rests on that "
	            "diagonality, and section 3.4 says in as many words that if "
	            "this fails the section is wrong and the plan needs the kernel" );

	BOOST_TEST( worstDiagonal < 1.0e-6,
	            "CEDRES++'s boundary form is diagonal but its diagonal is not "
	            "MEQ's symbol: worst relative disagreement " << worstDiagonal
	            << ". Suspect the elliptic MODULUS having become the parameter, "
	            "or the ( r1 r2 )^( 3/2 ) exponent, before suspecting the "
	            "separation of variables -- both are transcription errors that "
	            "converge to a wrong answer" );
}

/*
 * THE DECAY DIAGNOSTIC, AND THE MISREADING IT EXISTS TO PREVENT.
 *
 * FREE-BOUNDARY-PLAN.md section 11.6 asks for a readable spectrum, on the
 * grounds that the raw coefficients are not one: the modes are orthogonal in
 * dGamma/r but NOT normalised in it, so mode n contributes
 * | a_n | sqrt( mass( n ) ) to the trace's norm rather than | a_n |, and
 * mass( n ) ~ 1/n^3 makes that conversion factor ~ n^( -3/2 ).
 *
 * SO THIS TEST IS NOT ABOUT ARITHMETIC. Checking | a_n | sqrt( mass( n ) )
 * against | a_n | sqrt( mass( n ) ) would pass with the diagnostic doing
 * nothing useful. What has to be demonstrated is the DIVERGENCE OF THE TWO
 * VIEWS -- that raw and energy spectra disagree about which way the tail is
 * going -- and that is what the two constructed vectors below do:
 *
 *   1. raw coefficients all EQUAL. This is the picture that looks like a
 *      truncation going nowhere. In energy it is already falling like
 *      n^( -3/2 ), and truncationRatio() must report it as small.
 *   2. energy amplitudes all EQUAL. This is the genuinely unconverged case --
 *      the last mode carrying as much as the first -- and its raw coefficients
 *      GROW like n^( +3/2 ), which no reader would mistake for convergence but
 *      which a reader of raw numbers would also not recognise as the FLAT
 *      spectrum it is.
 *
 * The raw view is therefore wrong in both directions, and it is reassuring in
 * the one that matters. That is the finding, and the assertions below are its
 * statement.
 *
 * A NOTE ON A DIRECTION THAT IS EASY TO GET BACKWARDS. sqrt( mass( n ) )
 * DECREASES with n, so it maps flat-raw onto decaying-energy and NOT the other
 * way round. Getting this the wrong way about would invert the advice the
 * diagnostic gives, so case (2) is carried specifically to pin the sign of the
 * effect from the other end.
 */
BOOST_AUTO_TEST_CASE( the_energy_spectrum_and_the_raw_one_disagree_about_the_tail )
{
	int const spectrumModes = 10;
	ExteriorDtN const dtn( 0.0, rhoGamma, spectrumModes );

	// ---- (1) FLAT RAW: every coefficient the same ------------------------
	std::vector<double> const flatRaw(
		static_cast<std::size_t>( spectrumModes ), 1.0 );
	std::vector<double> const flatRawEnergy = dtn.modeAmplitudes( flatRaw );

	std::printf( "\n  a FLAT RAW spectrum, read both ways\n" );
	std::printf( "    %4s %14s %16s %14s\n", "n", "raw |a_n|",
	             "energy amplitude", "n^( -3/2 )" );
	for ( int i = 0; i < spectrumModes; ++i )
	{
		int const n = ExteriorDtN::firstMode() + i;
		std::size_t const k = static_cast<std::size_t>( i );
		std::printf( "    %4d %14.6e %16.6e %14.6e\n",
		             n, flatRaw[ k ], flatRawEnergy[ k ],
		             std::pow( static_cast<double>( n ), -1.5 ) );
	}
	std::fflush( stdout );

	// The energy spectrum falls MONOTONELY while the raw one does not move at
	// all. This is the whole claim.
	for ( int i = 1; i < spectrumModes; ++i )
		BOOST_TEST( flatRawEnergy[ static_cast<std::size_t>( i ) ]
		            < flatRawEnergy[ static_cast<std::size_t>( i - 1 ) ],
		            "the energy amplitude did not fall from mode "
		            << ExteriorDtN::firstMode() + i - 1 << " to "
		            << ExteriorDtN::firstMode() + i
		            << " on a FLAT raw spectrum. sqrt( mass( n ) ) decreases "
		            "with n, so it must -- if this fails the conversion factor "
		            "has the wrong sign of exponent and the diagnostic's advice "
		            "is inverted" );

	// AND IT IS THE n^( -3/2 ) OF THE PLAN, not merely some decay. mass( n )
	// is exactly 2/( n( n-1 )( 2n-1 ) ), which is asymptotically 1/n^3, so the
	// ratio of the amplitude to n^( -3/2 ) must SETTLE rather than drift.
	double const tailRatioFirst =
		flatRawEnergy[ 0 ]
		/std::pow( static_cast<double>( ExteriorDtN::firstMode() ), -1.5 );
	double const tailRatioLast =
		flatRawEnergy[ static_cast<std::size_t>( spectrumModes - 1 ) ]
		/std::pow( static_cast<double>( dtn.lastMode() ), -1.5 );
	std::printf( "    amplitude / n^( -3/2 ): %.6f at n = %d, %.6f at n = %d\n",
	             tailRatioFirst, ExteriorDtN::firstMode(),
	             tailRatioLast, dtn.lastMode() );
	// Asserted as an APPROACH as well as a bound: the ratio tends to 1 like
	// 1 + O( 1/n ), so a fixed tolerance is a statement about the mode count
	// while the monotone approach is a statement about the mass.
	BOOST_TEST( tailRatioLast < tailRatioFirst,
	            "the ratio of the energy amplitude to n^( -3/2 ) did not move "
	            "toward 1: " << tailRatioLast << " at n = " << dtn.lastMode()
	            << " against " << tailRatioFirst << " at n = "
	            << ExteriorDtN::firstMode() );
	BOOST_TEST( std::fabs( tailRatioLast - 1.0 ) < 0.1,
	            "the energy amplitude of a flat raw spectrum is not asymptotic "
	            "to n^( -3/2 ): the ratio reads " << tailRatioLast
	            << " at n = " << dtn.lastMode()
	            << " where it should approach 1. mass( n ) ~ 1/n^3 is what makes "
	            "sqrt( mass( n ) ) ~ n^( -3/2 ), so a departure here is a "
	            "departure in the mass. The approach is O( 1/n ), so raise the "
	            "mode count before loosening this" );

	// The truncation looks CONVERGED in energy, which a reader of the raw
	// numbers -- all exactly 1.0 -- would have concluded the opposite of.
	double const flatRawTruncation = dtn.truncationRatio( flatRaw );
	std::printf( "    truncation ratio: %.4e  ( raw coefficients are all 1.0 )\n",
	             flatRawTruncation );
	BOOST_TEST( flatRawTruncation < 0.2,
	            "a flat RAW spectrum reports a truncation ratio of "
	            << flatRawTruncation << ", which is not the small number its "
	            "energy content deserves. This is the misreading section 11.6 "
	            "exists to remove: flat raw IS decaying energy" );

	// ---- (2) FLAT ENERGY: the genuinely unconverged case -----------------
	//
	// Built by INVERTING the diagnostic -- a_n = 1/sqrt( mass( n ) ) -- so that
	// the energy amplitudes come out identically 1 and the raw coefficients are
	// whatever that costs.
	std::vector<double> flatEnergy(
		static_cast<std::size_t>( spectrumModes ), 0.0 );
	for ( int i = 0; i < spectrumModes; ++i )
		flatEnergy[ static_cast<std::size_t>( i ) ] =
			1.0/std::sqrt( dtn.mass( ExteriorDtN::firstMode() + i ) );

	std::vector<double> const flatEnergyAmplitudes =
		dtn.modeAmplitudes( flatEnergy );

	std::printf( "\n  a FLAT ENERGY spectrum, read both ways\n" );
	std::printf( "    %4s %14s %16s\n", "n", "raw |a_n|", "energy amplitude" );
	for ( int i = 0; i < spectrumModes; ++i )
		std::printf( "    %4d %14.6e %16.6e\n",
		             ExteriorDtN::firstMode() + i,
		             flatEnergy[ static_cast<std::size_t>( i ) ],
		             flatEnergyAmplitudes[ static_cast<std::size_t>( i ) ] );
	std::fflush( stdout );

	for ( int i = 0; i < spectrumModes; ++i )
		BOOST_TEST( std::fabs( flatEnergyAmplitudes[ static_cast<std::size_t>( i ) ]
		                       - 1.0 ) < 1.0e-14,
		            "inverting the diagnostic did not reproduce a flat energy "
		            "spectrum at mode " << ExteriorDtN::firstMode() + i << ": "
		            << flatEnergyAmplitudes[ static_cast<std::size_t>( i ) ]
		            << " against 1" );

	// AND ITS RAW COEFFICIENTS GROW, which is the other end of the same fact.
	for ( int i = 1; i < spectrumModes; ++i )
		BOOST_TEST( flatEnergy[ static_cast<std::size_t>( i ) ]
		            > flatEnergy[ static_cast<std::size_t>( i - 1 ) ],
		            "a flat ENERGY spectrum did not have growing raw "
		            "coefficients at mode " << ExteriorDtN::firstMode() + i );

	// The truncation is NOT converged here, and says so -- the last mode
	// carries exactly as much as the first.
	double const flatEnergyTruncation = dtn.truncationRatio( flatEnergy );
	std::printf( "    truncation ratio: %.4f  ( energy amplitudes are all 1.0 )\n",
	             flatEnergyTruncation );
	BOOST_TEST( std::fabs( flatEnergyTruncation - 1.0 ) < 1.0e-14,
	            "a spectrum carrying equal energy in every mode reported a "
	            "truncation ratio of " << flatEnergyTruncation
	            << " rather than 1. That is the unconverged case and the "
	            "diagnostic has to name it" );

	// ---- AND THE TWO VIEWS ARE FAR APART, WHICH IS THE POINT -------------
	//
	// If they agreed to within a modest factor the diagnostic would be a
	// refinement rather than a correction, and this test would not be worth
	// its cost.
	double const rawSpread = flatEnergy.back()/flatEnergy.front();
	std::printf( "\n    over %d modes the two views differ by a factor of %.1f\n",
	             spectrumModes, rawSpread );
	BOOST_TEST( rawSpread > 10.0,
	            "raw and energy readings of the same spectrum differ by only "
	            << rawSpread << " over " << spectrumModes << " modes, so the "
	            "diagnostic is not correcting anything worth correcting" );
}

/*
 * THE NORM, AND THE TWO DEGENERATE ANSWERS THAT MUST NOT BE NaN.
 *
 * traceNorm() is asserted against the orthogonality relation computed
 * INDEPENDENTLY here -- sqrt( sum a_n^2 mass( n ) ) written out by hand -- and
 * against the plain Euclidean norm of `a`, which is what a caller who ignores
 * the weight would reach for. The second is the control: if the two agreed,
 * the method would be doing nothing.
 *
 * The degenerate cases are the ones a driver actually meets. A coupled Newton
 * STARTS with every exterior coefficient at zero, so a diagnostic that printed
 * NaN there would print NaN on the first line of every free-boundary run.
 */
BOOST_AUTO_TEST_CASE( the_trace_norm_carries_the_weight_and_degrades_gracefully )
{
	ExteriorDtN const dtn = standard();

	std::vector<double> a( static_cast<std::size_t>( modes ), 0.0 );
	for ( int i = 0; i < modes; ++i )
		a[ static_cast<std::size_t>( i ) ] = 1.0/static_cast<double>( i + 1 );

	// The orthogonality relation, written out here rather than called.
	double expected = 0.0;
	for ( int i = 0; i < modes; ++i )
	{
		double const value = a[ static_cast<std::size_t>( i ) ];
		expected += value*value*dtn.mass( ExteriorDtN::firstMode() + i );
	}
	expected = std::sqrt( expected );

	double euclidean = 0.0;
	for ( double value : a )
		euclidean += value*value;
	euclidean = std::sqrt( euclidean );

	std::printf( "\n  the trace norm in dGamma/r: %.10e\n", dtn.traceNorm( a ) );
	std::printf( "    against the orthogonality relation: %.10e\n", expected );
	std::printf( "    against a plain Euclidean norm:     %.10e  ( the control )\n",
	             euclidean );
	std::fflush( stdout );

	BOOST_TEST( std::fabs( dtn.traceNorm( a ) - expected )
	            < 1.0e-14*std::fabs( expected ),
	            "traceNorm() is not sqrt( sum a_n^2 mass( n ) )" );

	BOOST_TEST( std::fabs( dtn.traceNorm( a ) - euclidean ) > 0.1*euclidean,
	            "the weighted norm agrees with the unweighted one to within "
	            "10%, so the weight is not being applied and every amplitude "
	            "above is a plain coefficient wearing a different name" );

	// ---- degenerate: all zero -------------------------------------------
	std::vector<double> const zero( static_cast<std::size_t>( modes ), 0.0 );
	BOOST_TEST( dtn.traceNorm( zero ) == 0.0 );
	BOOST_TEST( dtn.truncationRatio( zero ) == 0.0,
	            "an all-zero coefficient vector -- the state a coupled Newton "
	            "starts from -- did not report a truncation ratio of 0. It must "
	            "be printable rather than NaN" );

	// ---- degenerate: one mode -------------------------------------------
	ExteriorDtN const single( 0.0, rhoGamma, 1 );
	std::vector<double> const one( 1, 3.0 );
	BOOST_TEST( single.truncationRatio( one ) == 1.0,
	            "with a single mode the last IS the largest, so the honest "
	            "answer is 1 -- there is no evidence of decay to report" );

	// ---- and the size contract ------------------------------------------
	std::vector<double> const tooFew( static_cast<std::size_t>( modes - 1 ), 1.0 );
	BOOST_CHECK_THROW( dtn.modeAmplitudes( tooFew ), std::invalid_argument );
	BOOST_CHECK_THROW( dtn.traceNorm( tooFew ), std::invalid_argument );
	BOOST_CHECK_THROW( dtn.truncationRatio( tooFew ), std::invalid_argument );
}

/*
 * THE DIAGNOSTIC ON A REAL SPECTRUM, WHICH IS WHERE IT HAS TO EARN ITS PLACE.
 *
 * the_spectrum_falls_geometrically_and_faster_from_further_out already
 * establishes that a current loop's RAW coefficients fall as Gamma moves out.
 * This asks the question a caller actually has: at a given rho_Gamma, is N
 * enough? The truncation ratio must fall as Gamma moves outward, because the
 * spectrum does -- and it must do so read in energy, which is the only reading
 * that means anything.
 */
BOOST_AUTO_TEST_CASE( the_truncation_ratio_falls_as_gamma_moves_outward )
{
	using meq::analytic::CurrentLoop;
	CurrentLoop const loop = CurrentLoop::unitFlux();

	// THE ODD MODES ARE IDENTICALLY ZERO HERE AND THAT IS THE HAZARD THIS
	// CASE EXPOSED. The loop is centred on z = zCentre, so its trace is even in
	// mu, and C_n( -mu ) = ( -1 )^n C_n( mu ) makes every odd coefficient
	// vanish by SYMMETRY rather than by convergence. Twelve modes are degrees
	// 2..13, so the last is the odd 13 -- and a truncation summary reading the
	// last mode alone reported 0.000e+00 at every radius below while the raw
	// spectrum moved six orders. That is asserted first, because it is the
	// property the summary has to survive.
	std::printf( "\n  a unit loop is EVEN in mu: the odd modes vanish by parity\n" );
	{
		ExteriorDtN const dtn( 0.0, 2.5, 12 );
		auto trace = [ &loop ]( double r, double z ) { return loop.psi( r, z ); };
		std::vector<double> const a = dtn.coefficients( trace );
		std::vector<double> const amplitudes = dtn.modeAmplitudes( a );

		double worstOdd = 0.0;
		double smallestEven = std::numeric_limits<double>::max();
		for ( int i = 0; i < dtn.modeCount(); ++i )
		{
			double const value = amplitudes[ static_cast<std::size_t>( i ) ];
			if ( ( ExteriorDtN::firstMode() + i ) % 2 == 0 )
				smallestEven = std::min( smallestEven, value );
			else
				worstOdd = std::max( worstOdd, value );
		}
		std::printf( "    largest ODD amplitude %.3e, smallest EVEN %.3e\n",
		             worstOdd, smallestEven );
		std::fflush( stdout );

		BOOST_TEST( worstOdd < 1.0e-12*smallestEven,
		            "the odd modes of an up-down symmetric trace are not zero: "
		            "largest is " << worstOdd << " against a smallest even mode "
		            "of " << smallestEven << ". If this fails, the parity "
		            "C_n( -mu ) = ( -1 )^n C_n( mu ) is broken and the basis is "
		            "wrong, not the diagnostic" );

		BOOST_TEST( dtn.truncationRatio( a ) > 0.0,
		            "truncationRatio() returned exactly zero on an up-down "
		            "symmetric trace whose LAST retained degree is odd and so "
		            "identically absent. Zero reads as perfectly converged, and "
		            "half of all mode counts land on the absent parity -- which "
		            "is why the summary takes the larger of the last TWO "
		            "amplitudes rather than the last" );
	}

	std::printf( "\n  a unit loop: is N = 12 enough, by distance of Gamma?\n" );
	std::printf( "    %8s %14s %14s %16s\n", "rhoGamma", "raw |a_12/a_2|",
	             "trace norm", "truncation ratio" );

	std::vector<double> ratios;
	for ( double radius : { 2.0, 2.5, 4.0, 8.0 } )
	{
		ExteriorDtN const dtn( 0.0, radius, 12 );
		auto trace = [ &loop ]( double r, double z ) { return loop.psi( r, z ); };
		std::vector<double> const a = dtn.coefficients( trace );

		double const rawRatio = std::fabs( a[ 10 ] )/std::fabs( a[ 0 ] );
		double const ratio = dtn.truncationRatio( a );
		ratios.push_back( ratio );

		std::printf( "    %8.1f %14.3e %14.3e %16.3e\n",
		             radius, rawRatio, dtn.traceNorm( a ), ratio );
	}
	std::fflush( stdout );

	for ( std::size_t i = 1; i < ratios.size(); ++i )
		BOOST_TEST( ratios[ i ] < ratios[ i - 1 ],
		            "the truncation ratio did not fall when Gamma moved outward "
		            "at step " << i << ": " << ratios[ i ] << " against "
		            << ratios[ i - 1 ] << ". That is the trade a caller makes "
		            "between rho_Gamma and N, and this is the number they read "
		            "it off" );

	BOOST_TEST( ratios.back() < 1.0e-6,
	            "twelve modes seen from eight loop radii leave "
	            << ratios.back() << " of the spectrum at the truncation, which "
	            "is not converged" );
}
