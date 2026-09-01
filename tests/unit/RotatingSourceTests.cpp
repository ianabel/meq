#define BOOST_TEST_MODULE MeqRotatingSourceTests

#include <boost/test/unit_test.hpp>

#include <cmath>
#include <cstddef>
#include <functional>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

#include "meq/RotatingSource.hpp"
#include "meq/Source.hpp"

/*
 * FL-0 to FL-3 of FLOW-PLAN.md, which is to say everything about the rotating
 * source that can be asserted without a solver.
 *
 * The order below is the order of the plan, and it is also the order of blame:
 * if the charge-neutrality helpers are wrong then the closed form for phi_0 is
 * derived from a false premise, if phi_0 is wrong then the pressure is, and if
 * the pressure is wrong then so is F. Each group therefore assumes the one above
 * it has passed.
 */

namespace
{

	// Deuterium and electrons, in SI. Real masses rather than round numbers,
	// because the ratio m_e/m_i is what decides whether keeping the electron mass
	// in the closed form was worth the term it costs.
	double const deuteronMass = 3.3435837768e-27;
	double const electronMass = 9.1093837015e-31;
	double const keV = 1.602176634e-16;
	double const referenceRadius = 1.0;

	// omega0 chosen so that the exponent reaches about 1 at the outboard edge of
	// the benchmark box: m_i omega^2 ( r^2 - rRef^2 )/2( T_i + T_e ) with
	// r = 1.4, rRef = 1 and T_i + T_e = 3.7 keV.
	double const sonicOmega = 6.0e5;

	// A profile given by three closed forms, so a test can hand the source
	// derivatives that are exactly right by construction rather than differenced.
	class AnalyticProfile : public meq::Profile
	{
		public:
			AnalyticProfile( std::function<double( double )> value,
				std::function<double( double )> derivative,
				std::function<double( double )> secondDerivative )
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
				return secondDerivativeFunction( psi );
			}

		private:
			std::function<double( double )> valueFunction;
			std::function<double( double )> derivativeFunction;
			std::function<double( double )> secondDerivativeFunction;
	};

	std::shared_ptr<meq::Profile const> analytic( std::function<double( double )> f,
		std::function<double( double )> fPrime, std::function<double( double )> fDoublePrime )
	{
		return std::make_shared<AnalyticProfile const>( std::move( f ), std::move( fPrime ), std::move( fDoublePrime ) );
	}

	// The profiles. Non-polynomial, non-monotonic and all different shapes, so
	// that a dropped or mis-scaled term cannot hide behind another. Every one is
	// positive on [ 0, 1 ]: a species with a non-positive temperature has no
	// Maxwellian, and a negative density is not a plasma.

	double ionTemperature( double psi ) { return keV*( 2.0 - psi + 0.3*std::sin( 3.0*psi ) ); }
	double ionTemperaturePrime( double psi ) { return keV*( -1.0 + 0.9*std::cos( 3.0*psi ) ); }
	double ionTemperatureDoublePrime( double psi ) { return keV*( -2.7*std::sin( 3.0*psi ) ); }

	double electronTemperature( double psi ) { return keV*( 1.5 - 0.8*psi + 0.2*std::cos( 2.0*psi ) ); }
	double electronTemperaturePrime( double psi ) { return keV*( -0.8 - 0.4*std::sin( 2.0*psi ) ); }
	double electronTemperatureDoublePrime( double psi ) { return keV*( -0.8*std::cos( 2.0*psi ) ); }

	double ionDensity( double psi ) { return 1.0e19*( 1.2 - 0.5*psi*psi + 0.1*std::exp( -psi ) ); }
	double ionDensityPrime( double psi ) { return 1.0e19*( -psi - 0.1*std::exp( -psi ) ); }
	double ionDensityDoublePrime( double psi ) { return 1.0e19*( -1.0 + 0.1*std::exp( -psi ) ); }

	double rotation( double psi ) { return sonicOmega*( 1.0 - 0.4*psi + 0.15*std::sin( 2.0*psi ) ); }
	double rotationPrime( double psi ) { return sonicOmega*( -0.4 + 0.3*std::cos( 2.0*psi ) ); }
	double rotationDoublePrime( double psi ) { return sonicOmega*( -0.6*std::sin( 2.0*psi ) ); }

	double ggPrimeValue( double psi ) { return 0.4*std::exp( -2.0*psi ) + 0.25*std::tanh( 3.0*psi - 1.0 ); }
	double ggPrimeDerivative( double psi )
	{
		double const coshTerm = std::cosh( 3.0*psi - 1.0 );
		return -0.8*std::exp( -2.0*psi ) + 0.75/( coshTerm*coshTerm );
	}
	double ggPrimeSecondDerivative( double psi )
	{
		double const coshTerm = std::cosh( 3.0*psi - 1.0 );
		return 1.6*std::exp( -2.0*psi ) - 4.5*std::tanh( 3.0*psi - 1.0 )/( coshTerm*coshTerm );
	}

	std::shared_ptr<meq::Profile const> ionTemperatureProfile()
	{
		return analytic( ionTemperature, ionTemperaturePrime, ionTemperatureDoublePrime );
	}

	std::shared_ptr<meq::Profile const> electronTemperatureProfile()
	{
		return analytic( electronTemperature, electronTemperaturePrime, electronTemperatureDoublePrime );
	}

	std::shared_ptr<meq::Profile const> ionDensityProfile( double scale = 1.0 )
	{
		return analytic( [ scale ]( double psi ) { return scale*ionDensity( psi ); },
			[ scale ]( double psi ) { return scale*ionDensityPrime( psi ); },
			[ scale ]( double psi ) { return scale*ionDensityDoublePrime( psi ); } );
	}

	std::shared_ptr<meq::Profile const> rotationProfile( double scale = 1.0 )
	{
		return analytic( [ scale ]( double psi ) { return scale*rotation( psi ); },
			[ scale ]( double psi ) { return scale*rotationPrime( psi ); },
			[ scale ]( double psi ) { return scale*rotationDoublePrime( psi ); } );
	}

	std::shared_ptr<meq::Profile const> ggPrimeProfile()
	{
		return analytic( ggPrimeValue, ggPrimeDerivative, ggPrimeSecondDerivative );
	}

	/// Deuterium and electrons with Z_i = 1, so charge neutrality on the
	/// reference curve makes the electron density equal to the ion one. Built
	/// through neutralisingDensity() rather than by hand, which is the point of
	/// that helper.
	std::vector<meq::Species> hydrogenicSpecies( double densityScale = 1.0 )
	{
		std::vector<meq::Species> species( 2 );

		species[ 0 ].mass = deuteronMass;
		species[ 0 ].charge = 1.0;
		species[ 0 ].temperature = ionTemperatureProfile();
		species[ 0 ].density = ionDensityProfile( densityScale );

		species[ 1 ].mass = electronMass;
		species[ 1 ].charge = -1.0;
		species[ 1 ].temperature = electronTemperatureProfile();
		species[ 1 ].density = meq::neutralisingDensity( species, 1 );

		return species;
	}

	meq::RotatingSource makeSource( double omegaScale = 1.0 )
	{
		return meq::RotatingSource( hydrogenicSpecies(),
			omegaScale == 0.0 ? nullptr : rotationProfile( omegaScale ),
			ggPrimeProfile(), referenceRadius );
	}

	std::vector<double> const testRadii{ 0.6, 0.85, 1.0, 1.2, 1.4 };
	std::vector<double> const testFluxes{ 0.0, 0.15, 0.4, 0.62, 0.85, 1.0 };

	void checkClose( double actual, double expected, double tolerance, char const * what, double where )
	{
		double const scale = std::max( 1.0, std::fabs( expected ) );
		BOOST_CHECK_MESSAGE( std::fabs( actual - expected ) <= tolerance*scale,
			what << " at " << where << ": got " << actual << ", expected " << expected
			<< " (difference " << actual - expected << ", allowed " << tolerance*scale << ")" );
	}

	// A relative comparison, for quantities whose natural size is nowhere near
	// one -- a pressure in Pa, a density in m^-3.
	void checkRelative( double actual, double expected, double tolerance, char const * what, double r, double psi )
	{
		double const scale = std::max( std::fabs( expected ), std::fabs( actual ) );
		double const allowed = tolerance*std::max( scale, 1.0e-300 );
		BOOST_CHECK_MESSAGE( std::fabs( actual - expected ) <= allowed,
			what << " at r = " << r << ", psi = " << psi << ": got " << actual
			<< ", expected " << expected << " (relative error "
			<< std::fabs( actual - expected )/scale << ")" );
	}

}

BOOST_AUTO_TEST_SUITE( rotating_source_tests )

/*
 * FL-0: THE SPECIES CONTAINER AND THE CHARGE-NEUTRALITY SOLVE.
 *
 * Fixing the gauge phi_0( rRef ) = 0 removes exactly one function's worth of
 * freedom from the densities, so for n species there are n - 1 independent
 * density flux functions. neutralisingDensity() is how that last one is
 * obtained rather than guessed.
 */

BOOST_AUTO_TEST_CASE( theNeutralityResidualIsTheChargeWeightedSum )
{
	std::vector<meq::Species> const species = hydrogenicSpecies();

	for ( double psi : testFluxes )
	{
		double const byHand = 1.0*( *species[ 0 ].density )( psi ) - 1.0*( *species[ 1 ].density )( psi );
		checkClose( meq::chargeNeutralityResidual( species, psi ), byHand, 1.0e-14, "the neutrality residual", psi );
	}
}

BOOST_AUTO_TEST_CASE( aNeutralisingDensityClosesTheConstraint )
{
	std::vector<meq::Species> const species = hydrogenicSpecies();

	for ( double psi : testFluxes )
	{
		double const scale = std::fabs( ( *species[ 0 ].density )( psi ) );
		BOOST_CHECK_MESSAGE( std::fabs( meq::chargeNeutralityResidual( species, psi ) ) <= 1.0e-14*scale,
			"charge neutrality is violated at psi = " << psi << ": sum of Z_s n_s0 is "
			<< meq::chargeNeutralityResidual( species, psi ) << " against a species density of " << scale );
	}
}

BOOST_AUTO_TEST_CASE( aNeutralisingDensityIsExactAtEveryDerivativeLevel )
{
	// A combination differenced rather than combined would pass the value check
	// above and fail here, and the Jacobian is built on doublePrime().
	std::vector<meq::Species> const species = hydrogenicSpecies();
	meq::Profile const & electrons = *species[ 1 ].density;

	for ( double psi : testFluxes )
	{
		checkRelative( electrons.prime( psi ), ionDensityPrime( psi ), 1.0e-14, "the neutralising density's prime", 0.0, psi );
		checkRelative( electrons.doublePrime( psi ), ionDensityDoublePrime( psi ), 1.0e-14,
			"the neutralising density's doublePrime", 0.0, psi );
	}
}

BOOST_AUTO_TEST_CASE( aSpeciesSetThatViolatesNeutralityIsRefused )
{
	// The closed form for phi_0 is DERIVED from Z_1 n_10 = -Z_2 n_20. Without it
	// the answer is not approximate, it is wrong, so this has to throw rather
	// than warn.
	std::vector<meq::Species> species = hydrogenicSpecies();
	species[ 1 ].density = ionDensityProfile( 1.05 );

	BOOST_CHECK_THROW( meq::RotatingSource( species, rotationProfile(), ggPrimeProfile(), referenceRadius ),
		std::invalid_argument );
}

BOOST_AUTO_TEST_CASE( aZeroChargeSpeciesCannotBeNeutralised )
{
	std::vector<meq::Species> species = hydrogenicSpecies();
	species[ 1 ].charge = 0.0;

	BOOST_CHECK_THROW( meq::neutralisingDensity( species, 1 ), std::invalid_argument );
	BOOST_CHECK_THROW( meq::neutralisingDensity( species, 5 ), std::out_of_range );
}

BOOST_AUTO_TEST_CASE( theConstructorRefusesWhatItCannotSolve )
{
	std::vector<meq::Species> const good = hydrogenicSpecies();

	// Three species: the general phi_0 root find is FL-6 and is not written, so
	// refusing is the honest behaviour rather than silently using two of them.
	std::vector<meq::Species> three = good;
	three.push_back( good[ 0 ] );
	BOOST_CHECK_THROW( meq::RotatingSource( three, rotationProfile(), ggPrimeProfile(), referenceRadius ),
		std::invalid_argument );

	// Same-sign charges: quasineutrality then has no root at all.
	std::vector<meq::Species> sameSign = good;
	sameSign[ 1 ].charge = 1.0;
	BOOST_CHECK_THROW( meq::RotatingSource( sameSign, rotationProfile(), ggPrimeProfile(), referenceRadius ),
		std::invalid_argument );

	std::vector<meq::Species> badMass = good;
	badMass[ 0 ].mass = -1.0;
	BOOST_CHECK_THROW( meq::RotatingSource( badMass, rotationProfile(), ggPrimeProfile(), referenceRadius ),
		std::invalid_argument );

	std::vector<meq::Species> noTemperature = good;
	noTemperature[ 0 ].temperature = nullptr;
	BOOST_CHECK_THROW( meq::RotatingSource( noTemperature, rotationProfile(), ggPrimeProfile(), referenceRadius ),
		std::invalid_argument );

	BOOST_CHECK_THROW( meq::RotatingSource( good, rotationProfile(), nullptr, referenceRadius ),
		std::invalid_argument );
	BOOST_CHECK_THROW( meq::RotatingSource( good, rotationProfile(), ggPrimeProfile(), 0.0 ),
		std::invalid_argument );
	BOOST_CHECK_NO_THROW( meq::RotatingSource( good, rotationProfile(), ggPrimeProfile(), referenceRadius ) );
}

/*
 * FL-1: phi_0 AND THE PRESSURE.
 *
 * With two species (97) is linear in phi_0 after taking logarithms, so there is
 * no root find and no inner tolerance. These are the properties that says the
 * closed form is the solution of (97) rather than merely a plausible formula.
 */

BOOST_AUTO_TEST_CASE( thePotentialVanishesOnTheReferenceCurve )
{
	// Not "small": exactly zero. This is the gauge, and it is imposed by
	// construction rather than converged to.
	meq::RotatingSource const source = makeSource();

	for ( double psi : testFluxes )
		BOOST_CHECK_EQUAL( source.potential( referenceRadius, psi ), 0.0 );
}

BOOST_AUTO_TEST_CASE( quasineutralityHoldsAtEveryRadius )
{
	// THE CENTRAL PROPERTY OF THE CLOSURE. Charge neutrality is imposed on
	// r = rRef by construction; (97) is the statement that it survives to every
	// other radius, and it does because both species carry the same exponent.
	for ( double omegaScale : { 0.0, 0.5, 1.0, 2.0 } )
	{
		meq::RotatingSource const source = makeSource( omegaScale );

		for ( double r : testRadii )
		{
			for ( double psi : testFluxes )
			{
				double const ions = source.density( 0, r, psi );
				double const electrons = source.density( 1, r, psi );
				double const residual = ions - electrons;

				BOOST_CHECK_MESSAGE( std::fabs( residual ) <= 1.0e-13*ions,
					"quasineutrality fails at r = " << r << ", psi = " << psi
					<< ", omega scale " << omegaScale << ": n_i = " << ions
					<< ", n_e = " << electrons << ", relative residual " << residual/ions );
			}
		}
	}
}

BOOST_AUTO_TEST_CASE( theDensitiesReduceToTheirReferenceValuesOnTheReferenceCurve )
{
	meq::RotatingSource const source = makeSource();
	std::vector<meq::Species> const species = hydrogenicSpecies();

	for ( double psi : testFluxes )
	{
		checkRelative( source.density( 0, referenceRadius, psi ), ( *species[ 0 ].density )( psi ),
			1.0e-15, "the ion density on the reference curve", referenceRadius, psi );
		checkRelative( source.density( 1, referenceRadius, psi ), ( *species[ 1 ].density )( psi ),
			1.0e-15, "the electron density on the reference curve", referenceRadius, psi );
	}
}

BOOST_AUTO_TEST_CASE( thePressureMatchesTheIsothermalClosedForm )
{
	// Li & Zhu, Comput. Phys. Commun. 260 (2021) 107264, eq (8):
	// P = P0( psi ) exp[ m_i omega^2 ( r^2 - rRef^2 )/2T ] with T = T_i + T_e.
	// Written out here with the electron mass kept, which is the exact two
	// species answer rather than the m_e -> 0 one their paper quotes.
	meq::RotatingSource const source = makeSource();

	for ( double r : testRadii )
	{
		for ( double psi : testFluxes )
		{
			double const tSum = ionTemperature( psi ) + electronTemperature( psi );
			double const p0 = ionDensity( psi )*tSum;
			double const w = rotation( psi );
			double const exponent = ( deuteronMass + electronMass )*w*w
				*( r*r - referenceRadius*referenceRadius )/( 2.0*tSum );

			checkRelative( source.pressure( r, psi ), p0*std::exp( exponent ), 1.0e-14,
				"the pressure against the isothermal closed form", r, psi );
			checkRelative( source.densityExponent( r, psi ), exponent, 1.0e-14,
				"the shared density exponent", r, psi );
		}
	}
}

BOOST_AUTO_TEST_CASE( thePotentialIsTheOneThatBalancesTheCentrifugalDrift )
{
	// e phi_0 = omega^2 ( r^2 - rRef^2 )( m_1 T_2 - m_2 T_1 )/2( Z_1 T_2 - Z_2 T_1 ),
	// and substituting it back into (96) for each species separately must give
	// the same density this source reports. That is the check that the potential
	// and the densities are solutions of the SAME equation.
	meq::RotatingSource const source = makeSource();
	std::vector<meq::Species> const species = hydrogenicSpecies();

	for ( double r : testRadii )
	{
		for ( double psi : testFluxes )
		{
			double const y = source.potential( r, psi );
			double const delta = r*r - referenceRadius*referenceRadius;
			double const w = rotation( psi );

			for ( std::size_t s = 0; s < species.size(); ++s )
			{
				double const t = ( *species[ s ].temperature )( psi );
				double const exponent = species[ s ].mass*w*w*delta/( 2.0*t ) - species[ s ].charge*y/t;

				checkRelative( source.density( s, r, psi ),
					( *species[ s ].density )( psi )*std::exp( exponent ), 1.0e-13,
					"eq (96) rebuilt from the reported potential", r, psi );
			}
		}
	}
}

/*
 * FL-2: THE omega -> 0 COLLAPSE.
 *
 * The stage to protect. It exercises the species container, the chain rule,
 * phi_0, the Gaussian-to-SI conversion and the sign of Delta* against an answer
 * meq already has, and it is the reason a sign or a factor of 4 pi cannot hide.
 */

BOOST_AUTO_TEST_CASE( withoutRotationItReproducesTheMhdSource )
{
	// At omega = 0 every exponential is one, the densities become flux functions
	// again and p collapses to P0 = sum_s n_s0 T_s. So this source must equal
	// meq::MHDSource built on p' = P0', to round-off.
	auto const p0Prime = analytic(
		[]( double psi )
		{
			return ionDensityPrime( psi )*( ionTemperature( psi ) + electronTemperature( psi ) )
				+ ionDensity( psi )*( ionTemperaturePrime( psi ) + electronTemperaturePrime( psi ) );
		},
		[]( double psi )
		{
			return ionDensityDoublePrime( psi )*( ionTemperature( psi ) + electronTemperature( psi ) )
				+ 2.0*ionDensityPrime( psi )*( ionTemperaturePrime( psi ) + electronTemperaturePrime( psi ) )
				+ ionDensity( psi )*( ionTemperatureDoublePrime( psi ) + electronTemperatureDoublePrime( psi ) );
		},
		[]( double ) { return 0.0; } );

	meq::MHDSource const mhd( p0Prime, ggPrimeProfile() );

	// Both routes to "no rotation": a null omega profile, and one that is
	// identically zero. They must agree with each other as well as with MHDSource.
	meq::RotatingSource const noProfile( hydrogenicSpecies(), nullptr, ggPrimeProfile(), referenceRadius );
	meq::RotatingSource const zeroProfile( hydrogenicSpecies(), rotationProfile( 0.0 ), ggPrimeProfile(), referenceRadius );

	for ( double r : testRadii )
	{
		for ( double psi : testFluxes )
		{
			checkRelative( noProfile.f( r, 0.0, psi ), mhd.f( r, 0.0, psi ), 1.0e-13,
				"F with no rotation against MHDSource", r, psi );
			checkRelative( noProfile.dFdPsi( r, 0.0, psi ), mhd.dFdPsi( r, 0.0, psi ), 1.0e-13,
				"dF/dpsi with no rotation against MHDSource", r, psi );
			checkRelative( zeroProfile.f( r, 0.0, psi ), mhd.f( r, 0.0, psi ), 1.0e-13,
				"F with omega identically zero against MHDSource", r, psi );
			checkRelative( zeroProfile.dFdPsi( r, 0.0, psi ), mhd.dFdPsi( r, 0.0, psi ), 1.0e-13,
				"dF/dpsi with omega identically zero against MHDSource", r, psi );
		}
	}
}

BOOST_AUTO_TEST_CASE( theSourceIsTheRadialPressureGradientPlusGgPrime )
{
	// F = mu0 r^2 dp/dpsi|_r + g g', which is the collapsed form of RoPP (136).
	// Checking f() against its own two pieces is what localises a failure to one
	// of them rather than to "the source".
	meq::RotatingSource const source = makeSource();

	for ( double r : testRadii )
	{
		for ( double psi : testFluxes )
		{
			double const expected = meq::vacuumPermeability*r*r*source.dPressureDPsi( r, psi )
				+ ggPrimeValue( psi );

			checkRelative( source.f( r, 0.0, psi ), expected, 1.0e-14, "F against its two pieces", r, psi );
		}
	}
}

/*
 * FL-3: THE JACOBIAN.
 *
 * The load-bearing test of this whole stage. CLAUDE.md's standing finding is
 * that a wrong dF/dpsi leaves every error and every convergence rate unchanged
 * and moves only Newton's observed order -- so a rate table cannot see this and
 * a central difference is what can. Swept over Mach number, because the terms
 * being checked are exactly the ones that vanish at omega = 0.
 */

BOOST_AUTO_TEST_CASE( dPressureDPsiAgreesWithACentralDifference )
{
	for ( double omegaScale : { 0.0, 0.5, 1.0, 1.5, 2.0 } )
	{
		meq::RotatingSource const source = makeSource( omegaScale );
		double const h = 1.0e-6;

		for ( double r : testRadii )
		{
			for ( double psi : testFluxes )
			{
				double const difference = ( source.pressure( r, psi + h ) - source.pressure( r, psi - h ) )/( 2.0*h );

				checkRelative( source.dPressureDPsi( r, psi ), difference, 1.0e-7,
					"dp/dpsi against a central difference", r, psi );
			}
		}
	}
}

BOOST_AUTO_TEST_CASE( dFdPsiAgreesWithACentralDifference )
{
	// The step and tolerance are the O( h^2 ) floor of the difference formula,
	// as in SourceTests: agreement at two steps is what says the residual is the
	// difference formula's error rather than a discrepancy that survives h -> 0.
	for ( double omegaScale : { 0.0, 0.5, 1.0, 1.5, 2.0 } )
	{
		meq::RotatingSource const source = makeSource( omegaScale );

		for ( double h : { 1.0e-5, 1.0e-4 } )
		{
			for ( double r : testRadii )
			{
				for ( double psi : testFluxes )
				{
					double const difference = ( source.f( r, 0.0, psi + h ) - source.f( r, 0.0, psi - h ) )/( 2.0*h );

					BOOST_CHECK_MESSAGE(
						std::fabs( source.dFdPsi( r, 0.0, psi ) - difference )
							<= 1.0e-6*std::max( 1.0, std::fabs( difference ) ),
						"omega scale " << omegaScale << ", h = " << h << ": dFdPsi disagrees with a central "
						"difference at r = " << r << ", psi = " << psi << ": analytic "
						<< source.dFdPsi( r, 0.0, psi ) << ", difference " << difference
						<< " (error " << source.dFdPsi( r, 0.0, psi ) - difference << ")" );
				}
			}
		}
	}
}

BOOST_AUTO_TEST_CASE( theMachNumberReachesOrderOneOnTheBenchmarkBox )
{
	// A guard on the fixture rather than on the code: if the profiles above were
	// to drift to a Mach number of 0.01 then every test in this file would still
	// pass and none of them would be testing rotation. The exponent is M^2/2 in
	// the usual sense, so this asks for a real centrifugal effect at the
	// outboard edge.
	meq::RotatingSource const source = makeSource();
	double const exponent = source.densityExponent( 1.4, 0.0 );

	BOOST_CHECK_MESSAGE( exponent > 0.5,
		"the density exponent at the outboard edge is only " << exponent
		<< ", so these profiles are barely rotating and this file is not testing what it claims" );
	BOOST_CHECK_MESSAGE( source.density( 0, 1.4, 0.0 ) > 1.5*source.density( 0, 0.6, 0.0 ),
		"the centrifugal density asymmetry is negligible, so these profiles test nothing" );
}

BOOST_AUTO_TEST_SUITE_END()
