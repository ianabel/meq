// Unit tests for meq::Source and its implementations.
//
// The Grad-Shafranov right hand side
//
//     F( r, z, psi ) = mu0 r^2 p'( psi ) + ( g g' )( psi )
//
// is where a sign or a factor goes wrong silently: with a consistent-but-wrong F
// the solver still converges, to the wrong equilibrium, and nothing complains.
// Two families of test guard against that here.
//
//   * The Solov'ev source is checked against the closed form of HDG-GS-1 eq (10)
//     term by term, and against an MHDSource built from the constant profiles
//     that eq (10) comes from. If MHDSource's convention -- which of p' and g g'
//     carries the mu0, and which the r^2 -- were wrong, those two would disagree.
//   * dFdPsi is checked against a central difference of f(). meq closes the
//     nonlinearity with Newton rather than the papers' Picard iteration, so the
//     derivative is load-bearing, and it is exactly the quantity that can be
//     wrong without changing the converged answer. That test is deliberately the
//     longest thing in this file.

#define BOOST_TEST_MODULE SourceTests
#ifndef BOOST_TEST_DYN_LINK
#define BOOST_TEST_DYN_LINK
#endif

#include <boost/test/unit_test.hpp>

#include "meq/Source.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <memory>
#include <stdexcept>
#include <type_traits>
#include <vector>

namespace
{

	// A profile given by a pair of closed-form functions, so that a test can hand
	// a Source a derivative that is exactly right by construction. Also the
	// smallest possible check that Profile really is usable as an interface.
	class AnalyticProfile : public meq::Profile
	{
		public:
			// The second derivative is optional because meq::MHDSource never asks
			// for one -- it stores the pre-multiplied products, so F costs one
			// evaluation and dF/dpsi one prime(). A caller that does ask and did
			// not supply it gets a throw rather than a plausible zero.
			AnalyticProfile( std::function<double( double )> value, std::function<double( double )> derivative,
				std::function<double( double )> secondDerivative = {} )
				: valueFunction( std::move( value ) ), derivativeFunction( std::move( derivative ) ),
				  secondDerivativeFunction( std::move( secondDerivative ) )
			{
			}

			double operator()( double psi ) const override
			{
				return valueFunction( psi );
			}

			double prime( double psi ) const override
			{
				return derivativeFunction( psi );
			}

			double doublePrime( double psi ) const override
			{
				if ( !secondDerivativeFunction )
					throw std::logic_error( "AnalyticProfile: this profile was built without a second derivative" );

				return secondDerivativeFunction( psi );
			}

		private:
			std::function<double( double )> valueFunction;
			std::function<double( double )> derivativeFunction;
			std::function<double( double )> secondDerivativeFunction;
	};

	// p'( psi ) and p''( psi ): non-polynomial, non-monotonic, and with a
	// derivative that changes sign inside [ 0, 1 ], so that a dropped or
	// mis-scaled term cannot hide.
	double pressureDerivative( double psi )
	{
		return 1.0 + 0.8*std::sin( 4.0*psi ) - 0.5*psi*psi;
	}

	double pressureSecondDerivative( double psi )
	{
		return 3.2*std::cos( 4.0*psi ) - psi;
	}

	// ( g g' )( psi ) and its derivative, deliberately a different shape again.
	double ggPrimeValue( double psi )
	{
		return 0.4*std::exp( -2.0*psi ) + 0.25*std::tanh( 3.0*psi - 1.0 );
	}

	double ggPrimeDerivative( double psi )
	{
		double const coshTerm = std::cosh( 3.0*psi - 1.0 );
		return -0.8*std::exp( -2.0*psi ) + 0.75/( coshTerm*coshTerm );
	}

	std::shared_ptr<meq::Profile const> analyticPressureProfile( double scale = 1.0 )
	{
		return std::make_shared<AnalyticProfile>(
			[ scale ]( double psi ) { return scale*pressureDerivative( psi ); },
			[ scale ]( double psi ) { return scale*pressureSecondDerivative( psi ); } );
	}

	std::shared_ptr<meq::Profile const> analyticGGPrimeProfile()
	{
		return std::make_shared<AnalyticProfile>( ggPrimeValue, ggPrimeDerivative );
	}

	std::shared_ptr<meq::Profile const> constantProfile( double value )
	{
		return std::make_shared<meq::ConstantProfile>( value );
	}

	// Absolute comparison with a floor of one, so a quantity passing through zero
	// does not make a relative tolerance meaningless.
	void checkClose( double actual, double expected, double tolerance, char const * what, double where )
	{
		double const scale = std::max( 1.0, std::fabs( expected ) );
		BOOST_CHECK_MESSAGE( std::fabs( actual - expected ) <= tolerance*scale,
			what << " at " << where << ": got " << actual << ", expected " << expected
			<< " (difference " << actual - expected << ", allowed " << tolerance*scale << ")" );
	}

	// dF/dpsi by central differences of f() at fixed ( r, z ).
	double differenceDFdPsi( meq::Source const & source, double r, double z, double psi, double h )
	{
		return ( source.f( r, z, psi + h ) - source.f( r, z, psi - h ) )/( 2.0*h );
	}

	// The check that protects the Newton solve: the analytic derivative against a
	// difference of the residual's own source term, over a range of radii and
	// fluxes.
	void checkDerivativeAgainstDifferences( meq::Source const & source, std::vector<double> const & radii,
		double psiLower, double psiUpper, int samples, double h, double tolerance, char const * what )
	{
		for ( double r : radii )
		{
			for ( int i = 0; i <= samples; ++i )
			{
				double const psi = psiLower + ( psiUpper - psiLower )*i/samples;
				double const analytic = source.dFdPsi( r, 0.0, psi );
				double const difference = differenceDFdPsi( source, r, 0.0, psi, h );

				BOOST_CHECK_MESSAGE( std::fabs( analytic - difference ) <= tolerance*std::max( 1.0, std::fabs( difference ) ),
					what << ": dFdPsi disagrees with a central difference at r = " << r << ", psi = " << psi
					<< ": analytic " << analytic << ", difference " << difference
					<< " (error " << analytic - difference << ")" );
			}
		}
	}

	std::vector<double> const testRadii{ 0.25, 0.6, 1.0, 2.5, 6.2 };

}

BOOST_AUTO_TEST_SUITE( soloviev_source_tests )

BOOST_AUTO_TEST_CASE( f_matches_the_closed_form )
{
	// HDG-GS-1 eq (10): with mu0 p' = -C, g g' = -A and the flux normalised so
	// that A + C = 1, F = -( ( 1 - A ) r^2 + A ).
	for ( double a : { -0.52, -0.115, 0.0, 0.25, 1.0, 2.0 } )
	{
		meq::SolovievSource const source( a );

		BOOST_CHECK_EQUAL( source.a(), a );
		checkClose( source.a() + source.c(), 1.0, 1e-15, "A + C", a );

		for ( double r : testRadii )
		{
			for ( double z : { -1.7, 0.0, 0.9 } )
			{
				for ( double psi : { -0.5, 0.0, 0.37, 1.0, 2.0 } )
				{
					double const expected = -( ( 1.0 - a )*r*r + a );
					checkClose( source.f( r, z, psi ), expected, 1e-14, "Solov'ev F", r );
				}
			}
		}
	}
}

BOOST_AUTO_TEST_CASE( f_does_not_depend_on_z_or_psi )
{
	meq::SolovievSource const source( -0.115 );
	double const reference = source.f( 1.3, 0.0, 0.0 );

	for ( double z : { -3.0, -0.4, 0.0, 2.2 } )
		for ( double psi : { -2.0, 0.0, 0.5, 11.0 } )
			BOOST_CHECK_EQUAL( source.f( 1.3, z, psi ), reference );
}

BOOST_AUTO_TEST_CASE( dfdpsi_is_identically_zero )
{
	// Not "small": exactly zero. The Solov'ev problem is linear in psi, and a
	// Newton solve on it has to converge in a single step.
	for ( double a : { -0.52, 0.0, 0.5, 1.0 } )
	{
		meq::SolovievSource const source( a );

		for ( double r : testRadii )
			for ( double z : { -1.0, 0.0, 1.0 } )
				for ( double psi : { -1.0, 0.0, 0.5, 1.0, 4.0 } )
					BOOST_CHECK_EQUAL( source.dFdPsi( r, z, psi ), 0.0 );
	}
}

BOOST_AUTO_TEST_CASE( agrees_with_the_equivalent_mhd_source )
{
	// The Solov'ev profiles are mu0 p' = -C and g g' = -A with C = 1 - A. Building
	// an MHDSource from exactly those constants and getting the same F is a check
	// on MHDSource's convention: which profile carries the mu0 r^2, and which
	// stands alone.
	for ( double a : { -0.52, 0.0, 0.3, 1.0 } )
	{
		meq::SolovievSource const soloviev( a );
		double const c = 1.0 - a;

		// In normalised units, mu0 = 1, so p' = -C directly.
		meq::MHDSource const normalised( constantProfile( -c ), constantProfile( -a ), 1.0 );

		// And in SI, where p' = -C/mu0 and the g g' profile is untouched.
		meq::MHDSource const si( constantProfile( -c/meq::vacuumPermeability ), constantProfile( -a ) );

		for ( double r : testRadii )
		{
			checkClose( normalised.f( r, 0.4, 0.7 ), soloviev.f( r, 0.4, 0.7 ), 1e-14, "normalised MHDSource against Solov'ev", r );
			checkClose( si.f( r, 0.4, 0.7 ), soloviev.f( r, 0.4, 0.7 ), 1e-13, "SI MHDSource against Solov'ev", r );
			BOOST_CHECK_EQUAL( normalised.dFdPsi( r, 0.4, 0.7 ), 0.0 );
			BOOST_CHECK_EQUAL( si.dFdPsi( r, 0.4, 0.7 ), 0.0 );
		}
	}
}

BOOST_AUTO_TEST_CASE( rejects_a_non_finite_parameter )
{
	BOOST_CHECK_THROW( meq::SolovievSource( std::nan( "" ) ), std::invalid_argument );
	BOOST_CHECK_THROW( meq::SolovievSource( std::numeric_limits<double>::infinity() ), std::invalid_argument );
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( mhd_source_tests )

BOOST_AUTO_TEST_CASE( f_is_mu0_r_squared_pprime_plus_ggprime )
{
	// Isolate the two terms. The pressure term, and only the pressure term,
	// carries mu0 r^2; the g g' term is added as it stands. Getting this backwards
	// is the single most likely way to produce a converged, wrong equilibrium.
	{
		meq::MHDSource const pressureOnly( constantProfile( 2.0 ), constantProfile( 0.0 ), 3.0 );
		for ( double r : testRadii )
			checkClose( pressureOnly.f( r, 0.0, 0.5 ), 3.0*r*r*2.0, 1e-14, "pressure term", r );
	}

	{
		meq::MHDSource const currentOnly( constantProfile( 0.0 ), constantProfile( -1.25 ), 3.0 );
		for ( double r : testRadii )
		{
			// No mu0 and no r^2 on this one, at any radius.
			BOOST_CHECK_EQUAL( currentOnly.f( r, 0.0, 0.5 ), -1.25 );
		}
	}

	// Both together, with profiles that actually vary.
	meq::MHDSource const source( analyticPressureProfile(), analyticGGPrimeProfile(), 1.0 );
	for ( double r : testRadii )
	{
		for ( double psi : { 0.0, 0.15, 0.5, 0.9, 1.0 } )
		{
			double const expected = r*r*pressureDerivative( psi ) + ggPrimeValue( psi );
			checkClose( source.f( r, -0.3, psi ), expected, 1e-14, "F against its definition", psi );
		}
	}
}

BOOST_AUTO_TEST_CASE( mu0_defaults_to_the_si_value_and_can_be_overridden )
{
	meq::MHDSource const standard( constantProfile( 1.0 ), constantProfile( 0.0 ) );
	meq::MHDSource const normalised( constantProfile( 1.0 ), constantProfile( 0.0 ), 1.0 );

	// Written out rather than taken from the header, so that a change to the
	// constant has to be made deliberately in two places.
	double const siPermeability = 4.0e-7*3.14159265358979323846;

	checkClose( standard.mu0(), siPermeability, 1e-15, "default mu0", 0.0 );
	BOOST_CHECK_EQUAL( normalised.mu0(), 1.0 );
	checkClose( standard.f( 2.0, 0.0, 0.5 ), siPermeability*4.0, 1e-14, "F with the default mu0", 2.0 );
	BOOST_CHECK_EQUAL( normalised.f( 2.0, 0.0, 0.5 ), 4.0 );
}

BOOST_AUTO_TEST_CASE( dfdpsi_matches_the_profile_derivatives_term_by_term )
{
	// With one profile constant, dFdPsi must be exactly the other profile's
	// derivative, scaled the same way its value is scaled in F.
	{
		meq::MHDSource const pressureOnly( analyticPressureProfile(), constantProfile( 7.0 ), 2.5 );
		for ( double r : testRadii )
			for ( double psi : { 0.0, 0.2, 0.55, 1.0 } )
				checkClose( pressureOnly.dFdPsi( r, 0.0, psi ), 2.5*r*r*pressureSecondDerivative( psi ), 1e-14, "pressure contribution to dFdPsi", psi );
	}

	{
		meq::MHDSource const currentOnly( constantProfile( 7.0 ), analyticGGPrimeProfile(), 2.5 );
		for ( double r : testRadii )
			for ( double psi : { 0.0, 0.2, 0.55, 1.0 } )
				checkClose( currentOnly.dFdPsi( r, 0.0, psi ), ggPrimeDerivative( psi ), 1e-14, "g g' contribution to dFdPsi", psi );
	}

	// And the pressure contribution scales as r^2 while the other does not.
	meq::MHDSource const source( analyticPressureProfile(), analyticGGPrimeProfile(), 1.0 );
	double const psi = 0.42;
	double const atOne = source.dFdPsi( 1.0, 0.0, psi );
	double const atTwo = source.dFdPsi( 2.0, 0.0, psi );
	double const ggContribution = ggPrimeDerivative( psi );
	checkClose( atTwo - ggContribution, 4.0*( atOne - ggContribution ), 1e-13, "r^2 scaling of the pressure term", 2.0 );
}

BOOST_AUTO_TEST_CASE( dfdpsi_agrees_with_a_central_difference )
{
	// The test that protects the Newton solve. Both terms of F vary with psi here,
	// and mu0 = 1 so that neither term is numerically negligible against the
	// other: with the SI mu0 and an O(1) pressure profile the pressure term is
	// 1e-6 of the total, and an error in it would hide inside the tolerance.
	meq::MHDSource const source( analyticPressureProfile(), analyticGGPrimeProfile(), 1.0 );

	checkDerivativeAgainstDifferences( source, testRadii, 0.0, 1.0, 200, 1e-5, 1e-6, "analytic profiles, mu0 = 1" );

	// Again with a smaller step, to show the agreement is limited by the
	// difference formula rather than by a discrepancy that would survive h -> 0.
	checkDerivativeAgainstDifferences( source, testRadii, 0.0, 1.0, 50, 1e-4, 1e-6, "analytic profiles, h = 1e-4" );

	// And beyond the nominal [ 0, 1 ], where a Newton iterate can stray.
	checkDerivativeAgainstDifferences( source, testRadii, -0.5, 1.5, 100, 1e-5, 1e-6, "analytic profiles, psi outside [0,1]" );
}

BOOST_AUTO_TEST_CASE( dfdpsi_agrees_with_a_central_difference_in_si_units )
{
	// The same check with the real mu0, and a pressure profile of a realistic
	// size ( p' ~ 1/mu0 ) so that both terms of F are of comparable magnitude and
	// the pressure term is genuinely being tested.
	meq::MHDSource const source( analyticPressureProfile( 1.0/meq::vacuumPermeability ), analyticGGPrimeProfile() );

	checkDerivativeAgainstDifferences( source, testRadii, 0.0, 1.0, 200, 1e-5, 1e-6, "SI units" );
}

BOOST_AUTO_TEST_CASE( dfdpsi_agrees_with_a_central_difference_for_spline_profiles )
{
	// The configuration-file path: both profiles tabulated and interpolated. Here
	// dFdPsi is a derivative of the interpolants, so the difference is taken away
	// from the knots, where the interpolant is only C^1.
	unsigned int const intervals = 32;
	auto const pPrimeSpline = std::make_shared<meq::SplineProfile>( pressureDerivative, pressureSecondDerivative, intervals );
	auto const ggPrimeSpline = std::make_shared<meq::SplineProfile>( ggPrimeValue, ggPrimeDerivative, intervals );

	meq::MHDSource const source( pPrimeSpline, ggPrimeSpline, 1.0 );

	for ( double r : testRadii )
	{
		for ( unsigned int i = 0; i < intervals; ++i )
		{
			double const psi = ( i + 0.5 )/intervals;
			double const analytic = source.dFdPsi( r, 0.0, psi );
			double const difference = differenceDFdPsi( source, r, 0.0, psi, 1e-6 );
			checkClose( analytic, difference, 1e-6, "spline dFdPsi against a central difference", psi );
		}
	}

	// The spline source should also be close to the analytic one it was sampled
	// from -- a much weaker statement than the derivative check above, but it
	// catches a profile wired into the wrong slot.
	meq::MHDSource const exact( analyticPressureProfile(), analyticGGPrimeProfile(), 1.0 );
	for ( double r : testRadii )
	{
		for ( int i = 0; i <= 40; ++i )
		{
			double const psi = i/40.0;
			checkClose( source.f( r, 0.0, psi ), exact.f( r, 0.0, psi ), 1e-4, "spline F against analytic F", psi );
		}
	}
}

BOOST_AUTO_TEST_CASE( an_out_of_range_flux_is_clamped_rather_than_thrown )
{
	// A Newton iterate overshoots; the profiles clamp, so F stays finite and
	// nothing escapes into the quadrature loop.
	auto const pPrimeSpline = std::make_shared<meq::SplineProfile>( pressureDerivative, pressureSecondDerivative, 8 );
	auto const ggPrimeSpline = std::make_shared<meq::SplineProfile>( ggPrimeValue, ggPrimeDerivative, 8 );
	meq::MHDSource const source( pPrimeSpline, ggPrimeSpline, 1.0 );

	for ( double psi : { -12.0, -1.0, 2.0, 40.0 } )
	{
		BOOST_CHECK_NO_THROW( source.f( 1.5, 0.0, psi ) );
		BOOST_CHECK( std::isfinite( source.f( 1.5, 0.0, psi ) ) );
		// Beyond the table the profiles are constant, so the source contributes
		// nothing to the Jacobian there.
		BOOST_CHECK_EQUAL( source.dFdPsi( 1.5, 0.0, psi ), 0.0 );
	}

	double const clamped = source.f( 1.5, 0.0, 5.0 );
	checkClose( clamped, source.f( 1.5, 0.0, 1.0 ), 1e-14, "clamped F above the table", 5.0 );
}

BOOST_AUTO_TEST_CASE( rejects_null_profiles_and_a_non_finite_mu0 )
{
	std::shared_ptr<meq::Profile const> const good = constantProfile( 1.0 );
	std::shared_ptr<meq::Profile const> const none;

	BOOST_CHECK_THROW( meq::MHDSource( none, good ), std::invalid_argument );
	BOOST_CHECK_THROW( meq::MHDSource( good, none ), std::invalid_argument );
	BOOST_CHECK_THROW( meq::MHDSource( none, none ), std::invalid_argument );
	BOOST_CHECK_THROW( meq::MHDSource( good, good, std::nan( "" ) ), std::invalid_argument );
	BOOST_CHECK_NO_THROW( meq::MHDSource( good, good ) );
}

BOOST_AUTO_TEST_CASE( shares_ownership_of_its_profiles )
{
	auto pressure = std::make_shared<meq::ConstantProfile>( 3.0 );
	auto current = std::make_shared<meq::ConstantProfile>( -1.0 );

	meq::MHDSource const source( pressure, current, 1.0 );
	BOOST_CHECK_EQUAL( pressure.use_count(), 2 );

	pressure.reset();
	current.reset();

	// The source keeps them alive.
	BOOST_CHECK_EQUAL( source.f( 2.0, 0.0, 0.5 ), 4.0*3.0 - 1.0 );
	BOOST_CHECK_EQUAL( source.pPrime()( 0.5 ), 3.0 );
	BOOST_CHECK_EQUAL( source.ggPrime()( 0.5 ), -1.0 );
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( normalised_mhd_source_tests )

/*
 * THE NORMALISATION IS ONE FACTOR IN f() AND TWO IN dFdPsi(), and getting the
 * second one wrong is the classic error: it leaves the converged answer exactly
 * where it was and costs only the order of the Newton convergence to it. See
 * CLAUDE.md, "A wrong Jacobian is invisible to a convergence table".
 */
BOOST_AUTO_TEST_CASE( f_carries_one_factor_of_the_normalisation )
{
	double const psiAxis = 0.37;
	meq::NormalisedMHDSource source( analyticPressureProfile(), analyticGGPrimeProfile(),
	                                 psiAxis, 1.0 );

	for ( double r : testRadii )
	{
		for ( int i = 0; i <= 10; ++i )
		{
			double const psi = psiAxis*i/10.0;
			double const psiN = psi/psiAxis;
			double const expected = ( r*r*pressureDerivative( psiN ) + ggPrimeValue( psiN ) )
			                        /psiAxis;
			checkClose( source.f( r, 0.0, psi ), expected, 1e-14, "normalised F", psi );
		}
	}
}

BOOST_AUTO_TEST_CASE( dfdpsi_carries_two )
{
	double const psiAxis = 0.37;
	meq::NormalisedMHDSource source( analyticPressureProfile(), analyticGGPrimeProfile(),
	                                 psiAxis, 1.0 );

	for ( double r : testRadii )
	{
		for ( int i = 0; i <= 10; ++i )
		{
			double const psi = psiAxis*i/10.0;
			double const psiN = psi/psiAxis;
			double const expected = ( r*r*pressureSecondDerivative( psiN )
			                          + ggPrimeDerivative( psiN ) )/( psiAxis*psiAxis );
			checkClose( source.dFdPsi( r, 0.0, psi ), expected, 1e-14,
			            "normalised dF/dpsi", psi );
		}
	}
}

/*
 * AND THE FINITE DIFFERENCE, AT SEVERAL NORMALISATIONS -- because this check is
 * the one that CANNOT see the thing the bordered Newton exists for. f() and
 * dFdPsi() are both evaluated at whatever normalisation is currently set, so
 * they agree with each other however wrong that normalisation is. What a
 * difference of f() in psi cannot detect is a missing dF/dpsi_ax, and nothing in
 * a unit test can: only a difference of the ASSEMBLED residual can, which is
 * what tests/convergence/HighBetaConvergence.cpp measures.
 */
BOOST_AUTO_TEST_CASE( dfdpsi_agrees_with_a_central_difference_at_every_normalisation )
{
	meq::NormalisedMHDSource source( analyticPressureProfile(), analyticGGPrimeProfile(),
	                                 1.0, 1.0 );

	for ( double psiAxis : { 0.05, 0.5, 1.0, 4.0, 40.0 } )
	{
		source.setNormalisation( psiAxis );
		BOOST_CHECK_EQUAL( source.normalisation(), psiAxis );
		checkDerivativeAgainstDifferences( source, testRadii, 0.05*psiAxis, 0.95*psiAxis,
		                                   20, 1e-7*psiAxis, 1e-6,
		                                   "NormalisedMHDSource" );
	}
}

/*
 * A NORMALISATION OF ZERO IS A THROW AND NOT AN INFINITY. The degenerate fixed
 * point the outer iteration used to find is psi_ax running to zero with psi
 * shrinking beside it, so an iterate that reaches it must say so rather than
 * quietly returning 1/0 -- the bordered Newton's line search catches this throw
 * and rejects the step that caused it.
 */
BOOST_AUTO_TEST_CASE( rejects_a_normalisation_of_zero )
{
	BOOST_CHECK_THROW( meq::NormalisedMHDSource( analyticPressureProfile(),
	                                             analyticGGPrimeProfile(), 0.0, 1.0 ),
	                   std::invalid_argument );

	meq::NormalisedMHDSource source( analyticPressureProfile(), analyticGGPrimeProfile(),
	                                 1.0, 1.0 );
	BOOST_CHECK_THROW( source.setNormalisation( 0.0 ), std::invalid_argument );
	BOOST_CHECK_THROW( source.setNormalisation( std::numeric_limits<double>::quiet_NaN() ),
	                   std::invalid_argument );

	// And it is unchanged by a rejected attempt, so a caught throw leaves a
	// usable source behind.
	BOOST_CHECK_EQUAL( source.normalisation(), 1.0 );
}

BOOST_AUTO_TEST_CASE( rejects_null_profiles )
{
	BOOST_CHECK_THROW( meq::NormalisedMHDSource( nullptr, analyticGGPrimeProfile(), 1.0, 1.0 ),
	                   std::invalid_argument );
	BOOST_CHECK_THROW( meq::NormalisedMHDSource( analyticPressureProfile(), nullptr, 1.0, 1.0 ),
	                   std::invalid_argument );
}

/*
 * THE SAME SOURCE AS MHDSource WHEN THE NORMALISATION IS ONE, which is the
 * cross-check that no sign or factor drifted between the two classes. At
 * psi_ax = 1 the chain rule contributes nothing and the two must agree exactly.
 */
BOOST_AUTO_TEST_CASE( agrees_with_the_unnormalised_source_at_unit_normalisation )
{
	meq::MHDSource plain( analyticPressureProfile(), analyticGGPrimeProfile(), 1.0 );
	meq::NormalisedMHDSource normalised( analyticPressureProfile(), analyticGGPrimeProfile(),
	                                     1.0, 1.0 );

	for ( double r : testRadii )
	{
		for ( int i = 0; i <= 10; ++i )
		{
			double const psi = i/10.0;
			BOOST_CHECK_EQUAL( normalised.f( r, 0.0, psi ), plain.f( r, 0.0, psi ) );
			BOOST_CHECK_EQUAL( normalised.dFdPsi( r, 0.0, psi ), plain.dFdPsi( r, 0.0, psi ) );
		}
	}
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( source_interface_tests )

BOOST_AUTO_TEST_CASE( implementations_are_usable_polymorphically )
{
	// The class this replaces inherited *privately* from its abstract base, so it
	// could never be held or called as one -- which is why the plasma model was
	// never wired up. These conversions are the regression test: they do not
	// compile under private inheritance.
	static_assert( std::is_base_of_v<meq::Source, meq::MHDSource>, "MHDSource must be a Source" );
	static_assert( std::is_base_of_v<meq::Source, meq::SolovievSource>, "SolovievSource must be a Source" );
	static_assert( std::is_base_of_v<meq::Source, meq::NormalisedSource>, "NormalisedSource must be a Source" );
	static_assert( std::is_base_of_v<meq::NormalisedSource, meq::NormalisedMHDSource>, "NormalisedMHDSource must be a NormalisedSource" );
	static_assert( std::is_convertible_v<meq::MHDSource *, meq::Source *>, "MHDSource must convert to its base" );
	static_assert( std::is_convertible_v<meq::SolovievSource *, meq::Source *>, "SolovievSource must convert to its base" );
	static_assert( std::has_virtual_destructor_v<meq::Source>, "a Source is deleted through a base pointer" );

	std::vector<std::unique_ptr<meq::Source>> sources;
	sources.push_back( std::make_unique<meq::SolovievSource>( -0.115 ) );
	sources.push_back( std::make_unique<meq::MHDSource>( analyticPressureProfile(), analyticGGPrimeProfile(), 1.0 ) );
	sources.push_back( std::make_unique<meq::NormalisedMHDSource>( analyticPressureProfile(), analyticGGPrimeProfile(), 0.8, 1.0 ) );

	for ( auto const & source : sources )
	{
		BOOST_CHECK( std::isfinite( source->f( 1.4, 0.2, 0.6 ) ) );
		BOOST_CHECK( std::isfinite( source->dFdPsi( 1.4, 0.2, 0.6 ) ) );

		// Every Source, whatever it is, must have a derivative consistent with its
		// own value: this is the loop the Newton solve will run inside.
		checkDerivativeAgainstDifferences( *source, { 1.4 }, 0.05, 0.95, 20, 1e-5, 1e-6, "through Source&" );
	}
}

BOOST_AUTO_TEST_SUITE_END()
