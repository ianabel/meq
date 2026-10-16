#define BOOST_TEST_MODULE MeqRotatingSourceTests

#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <cstdio>
#include <limits>
#include <string>
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
 * Everything about the rotating source that can be asserted without a solver.
 *
 * The order below is the order the pieces were built in, and it is also the
 * order of blame:
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
	// the benchmark box: m_i omega^2 ( R^2 - R_ref^2 )/2( T_i + T_e ) with
	// R = 1.4, R_ref = 1 and T_i + T_e = 3.7 keV.
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

	double carbonTemperature( double psi ) { return keV*( 1.9 - 0.7*psi + 0.25*std::sin( 2.5*psi ) ); }
	double carbonTemperaturePrime( double psi ) { return keV*( -0.7 + 0.625*std::cos( 2.5*psi ) ); }
	double carbonTemperatureDoublePrime( double psi ) { return keV*( -1.5625*std::sin( 2.5*psi ) ); }

	double carbonDensity( double psi ) { return 1.0e17*( 4.0 - 1.5*psi*psi + 0.5*std::cos( psi ) ); }
	double carbonDensityPrime( double psi ) { return 1.0e17*( -3.0*psi - 0.5*std::sin( psi ) ); }
	double carbonDensityDoublePrime( double psi ) { return 1.0e17*( -3.0 - 0.5*std::cos( psi ) ); }

	/// Deuterium, fully stripped carbon and electrons, with the electron density
	/// closed by neutrality. Three species, so (97) is transcendental and the
	/// general root find is the only route.
	std::vector<meq::Species> impuritySpecies()
	{
		std::vector<meq::Species> species( 3 );

		species[ 0 ].mass = deuteronMass;
		species[ 0 ].charge = 1.0;
		species[ 0 ].temperature = ionTemperatureProfile();
		species[ 0 ].density = ionDensityProfile();

		species[ 1 ].mass = 6.0*deuteronMass;
		species[ 1 ].charge = 6.0;
		species[ 1 ].temperature = analytic( carbonTemperature, carbonTemperaturePrime, carbonTemperatureDoublePrime );
		species[ 1 ].density = analytic( carbonDensity, carbonDensityPrime, carbonDensityDoublePrime );

		species[ 2 ].mass = electronMass;
		species[ 2 ].charge = -1.0;
		species[ 2 ].temperature = electronTemperatureProfile();
		species[ 2 ].density = meq::neutralisingDensity( species, 2 );

		return species;
	}

	meq::RotatingSource makeImpuritySource( double omegaScale = 1.0 )
	{
		return meq::RotatingSource( impuritySpecies(),
			omegaScale == 0.0 ? nullptr : rotationProfile( omegaScale ),
			ggPrimeProfile(), referenceRadius );
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
	//
	// THE FLOOR IS NOT OPTIONAL FOR A QUANTITY THAT VANISHES SOMEWHERE. phi_0 is
	// identically zero on R = R_ref, so at that radius both sides are round-off
	// and a purely relative comparison asks whether one 1e-33 equals another --
	// which is always false and never interesting. Callers pass the quantity's
	// own scale on this problem, so the test reads "agree to `tolerance`
	// relatively, or to `tolerance` of the scale, whichever is looser".
	void checkRelative( double actual, double expected, double tolerance, char const * what, double radius, double psi,
		double floor = 0.0 )
	{
		double const scale = std::max( std::max( std::fabs( expected ), std::fabs( actual ) ), std::fabs( floor ) );
		double const allowed = tolerance*std::max( scale, 1.0e-300 );
		BOOST_CHECK_MESSAGE( std::fabs( actual - expected ) <= allowed,
			what << " at R = " << radius << ", psi = " << psi << ": got " << actual
			<< ", expected " << expected << " (relative error "
			<< std::fabs( actual - expected )/scale << ")" );
	}

}

BOOST_AUTO_TEST_SUITE( rotating_source_tests )

/*
 * FL-0: THE SPECIES CONTAINER AND THE CHARGE-NEUTRALITY SOLVE.
 *
 * Fixing the gauge phi_0( R_ref ) = 0 removes exactly one function's worth of
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
	// R = R_ref by construction; (97) is the statement that it survives to every
	// other radius, and it does because both species carry the same exponent.
	for ( double omegaScale : { 0.0, 0.5, 1.0, 2.0 } )
	{
		meq::RotatingSource const source = makeSource( omegaScale );

		for ( double radius : testRadii )
		{
			for ( double psi : testFluxes )
			{
				double const ions = source.density( 0, radius, psi );
				double const electrons = source.density( 1, radius, psi );
				double const residual = ions - electrons;

				BOOST_CHECK_MESSAGE( std::fabs( residual ) <= 1.0e-13*ions,
					"quasineutrality fails at R = " << radius << ", psi = " << psi
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
	// P = P0( psi ) exp[ m_i omega^2 ( R^2 - R_ref^2 )/2T ] with T = T_i + T_e.
	// Written out here with the electron mass kept, which is the exact two
	// species answer rather than the m_e -> 0 one their paper quotes.
	meq::RotatingSource const source = makeSource();

	for ( double radius : testRadii )
	{
		for ( double psi : testFluxes )
		{
			double const tSum = ionTemperature( psi ) + electronTemperature( psi );
			double const p0 = ionDensity( psi )*tSum;
			double const w = rotation( psi );
			double const exponent = ( deuteronMass + electronMass )*w*w
				*( radius*radius - referenceRadius*referenceRadius )/( 2.0*tSum );

			checkRelative( source.pressure( radius, psi ), p0*std::exp( exponent ), 1.0e-14,
				"the pressure against the isothermal closed form", radius, psi );
			checkRelative( source.densityExponent( 0, radius, psi ), exponent, 1.0e-14,
				"the ion density exponent", radius, psi );

			// AND THE TWO SPECIES SHARE IT. That is the property, not an
			// implementation detail: it is what makes Sum_s Z_s n_s vanish at
			// every R once it vanishes at R_ref.
			checkRelative( source.densityExponent( 1, radius, psi ), source.densityExponent( 0, radius, psi ),
				1.0e-15, "the electron exponent against the ion one", radius, psi );
		}
	}
}

BOOST_AUTO_TEST_CASE( thePotentialIsTheOneThatBalancesTheCentrifugalDrift )
{
	// e phi_0 = omega^2 ( R^2 - R_ref^2 )( m_1 T_2 - m_2 T_1 )/2( Z_1 T_2 - Z_2 T_1 ),
	// and substituting it back into (96) for each species separately must give
	// the same density this source reports. That is the check that the potential
	// and the densities are solutions of the SAME equation.
	meq::RotatingSource const source = makeSource();
	std::vector<meq::Species> const species = hydrogenicSpecies();

	for ( double radius : testRadii )
	{
		for ( double psi : testFluxes )
		{
			double const y = source.potential( radius, psi );
			double const delta = radius*radius - referenceRadius*referenceRadius;
			double const w = rotation( psi );

			for ( std::size_t s = 0; s < species.size(); ++s )
			{
				double const t = ( *species[ s ].temperature )( psi );
				double const exponent = species[ s ].mass*w*w*delta/( 2.0*t ) - species[ s ].charge*y/t;

				checkRelative( source.density( s, radius, psi ),
					( *species[ s ].density )( psi )*std::exp( exponent ), 1.0e-13,
					"eq (96) rebuilt from the reported potential", radius, psi );
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

	for ( double radius : testRadii )
	{
		for ( double psi : testFluxes )
		{
			checkRelative( noProfile.f( radius, 0.0, psi ), mhd.f( radius, 0.0, psi ), 1.0e-13,
				"F with no rotation against MHDSource", radius, psi );
			checkRelative( noProfile.dFdPsi( radius, 0.0, psi ), mhd.dFdPsi( radius, 0.0, psi ), 1.0e-13,
				"dF/dpsi with no rotation against MHDSource", radius, psi );
			checkRelative( zeroProfile.f( radius, 0.0, psi ), mhd.f( radius, 0.0, psi ), 1.0e-13,
				"F with omega identically zero against MHDSource", radius, psi );
			checkRelative( zeroProfile.dFdPsi( radius, 0.0, psi ), mhd.dFdPsi( radius, 0.0, psi ), 1.0e-13,
				"dF/dpsi with omega identically zero against MHDSource", radius, psi );
		}
	}
}

BOOST_AUTO_TEST_CASE( theSourceIsTheRadialPressureGradientPlusGgPrime )
{
	// F = mu0 R^2 dp/dpsi|_r + g g', which is the collapsed form of RoPP (136).
	// Checking f() against its own two pieces is what localises a failure to one
	// of them rather than to "the source".
	meq::RotatingSource const source = makeSource();

	for ( double radius : testRadii )
	{
		for ( double psi : testFluxes )
		{
			double const expected = meq::vacuumPermeability*radius*radius*source.dPressureDPsi( radius, psi )
				+ ggPrimeValue( psi );

			checkRelative( source.f( radius, 0.0, psi ), expected, 1.0e-14, "F against its two pieces", radius, psi );
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

		for ( double radius : testRadii )
		{
			for ( double psi : testFluxes )
			{
				double const difference = ( source.pressure( radius, psi + h ) - source.pressure( radius, psi - h ) )/( 2.0*h );

				checkRelative( source.dPressureDPsi( radius, psi ), difference, 1.0e-7,
					"dp/dpsi against a central difference", radius, psi );
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
			for ( double radius : testRadii )
			{
				for ( double psi : testFluxes )
				{
					double const difference = ( source.f( radius, 0.0, psi + h ) - source.f( radius, 0.0, psi - h ) )/( 2.0*h );

					BOOST_CHECK_MESSAGE(
						std::fabs( source.dFdPsi( radius, 0.0, psi ) - difference )
							<= 1.0e-6*std::max( 1.0, std::fabs( difference ) ),
						"omega scale " << omegaScale << ", h = " << h << ": dFdPsi disagrees with a central "
						"difference at R = " << radius << ", psi = " << psi << ": analytic "
						<< source.dFdPsi( radius, 0.0, psi ) << ", difference " << difference
						<< " (error " << source.dFdPsi( radius, 0.0, psi ) - difference << ")" );
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
	double const exponent = source.densityExponent( 0, 1.4, 0.0 );

	BOOST_CHECK_MESSAGE( exponent > 0.5,
		"the density exponent at the outboard edge is only " << exponent
		<< ", so these profiles are barely rotating and this file is not testing what it claims" );
	BOOST_CHECK_MESSAGE( source.density( 0, 1.4, 0.0 ) > 1.5*source.density( 0, 0.6, 0.0 ),
		"the centrifugal density asymmetry is negligible, so these profiles test nothing" );
}

/*
 * FL-6: n SPECIES.
 *
 * Above two, (97) is transcendental and phi_0 comes from a safeguarded scalar
 * Newton, with its two psi-derivatives by implicit differentiation. The first
 * test below is the one that matters: at two species the general path and the
 * closed form are the same problem solved twice, so they must agree, and if they
 * do not then one of them is wrong and nothing else here can say which.
 */

BOOST_AUTO_TEST_CASE( theRootFindReproducesTheClosedFormAtTwoSpecies )
{
	for ( double omegaScale : { 0.0, 0.5, 1.0, 2.0 } )
	{
		meq::RotatingSource const closed( hydrogenicSpecies(),
			omegaScale == 0.0 ? nullptr : rotationProfile( omegaScale ),
			ggPrimeProfile(), referenceRadius, meq::vacuumPermeability,
			meq::RotatingSource::Closure::ClosedForm );
		meq::RotatingSource const found( hydrogenicSpecies(),
			omegaScale == 0.0 ? nullptr : rotationProfile( omegaScale ),
			ggPrimeProfile(), referenceRadius, meq::vacuumPermeability,
			meq::RotatingSource::Closure::RootFind );

		// phi_0 vanishes identically at R_ref, and identically everywhere when
		// omega is zero, so the floor has to be the problem's own energy scale
		// rather than anything read off this configuration. e phi_0 is in Joules,
		// so a keV is the scale, and 1e-12 of a keV is not a potential.

		for ( double radius : testRadii )
		{
			for ( double psi : testFluxes )
			{
				checkRelative( found.potential( radius, psi ), closed.potential( radius, psi ), 1.0e-12,
					"phi_0 by root find against the closed form", radius, psi, keV );
				checkRelative( found.dPotentialDPsi( radius, psi ), closed.dPotentialDPsi( radius, psi ), 1.0e-11,
					"dphi_0/dpsi by implicit differentiation against the closed form", radius, psi, keV );
				checkRelative( found.f( radius, 0.0, psi ), closed.f( radius, 0.0, psi ), 1.0e-11,
					"F by root find against the closed form", radius, psi );
				checkRelative( found.dFdPsi( radius, 0.0, psi ), closed.dFdPsi( radius, 0.0, psi ), 1.0e-9,
					"dF/dpsi by root find against the closed form", radius, psi );
			}
		}
	}
}

BOOST_AUTO_TEST_CASE( quasineutralityHoldsForThreeSpecies )
{
	// D, C6+ and electrons. The carbon is what makes this a real test: it is six
	// times as charged and six times as heavy as the deuterium, so it feels the
	// centrifugal force hardest and is held back by the field hardest.
	meq::RotatingSource const source = makeImpuritySource();

	for ( double radius : testRadii )
	{
		for ( double psi : testFluxes )
		{
			double charge = 0.0;
			double magnitude = 0.0;
			for ( std::size_t sp = 0; sp < 3; ++sp )
			{
				double const n = source.density( sp, radius, psi );
				charge += source.species()[ sp ].charge*n;
				magnitude += std::fabs( source.species()[ sp ].charge*n );
			}

			BOOST_CHECK_MESSAGE( std::fabs( charge ) <= 1.0e-12*magnitude,
				"three-species quasineutrality fails at R = " << radius << ", psi = " << psi
				<< ": sum of Z_s n_s is " << charge << " against " << magnitude );
		}
	}
}

BOOST_AUTO_TEST_CASE( dPotentialDPsiAgreesWithACentralDifference )
{
	// The inner analogue of the dFdPsi check, and the thing that catches an
	// implicit-differentiation slip. Differencing the root find from outside is
	// what this must NOT be doing, so a disagreement here means the closed-form
	// derivative and the numerical one have parted company.
	meq::RotatingSource const source = makeImpuritySource();
	double const h = 1.0e-7;

	for ( double radius : testRadii )
	{
		for ( double psi : testFluxes )
		{
			double const difference = ( source.potential( radius, psi + h ) - source.potential( radius, psi - h ) )/( 2.0*h );

			checkRelative( source.dPotentialDPsi( radius, psi ), difference, 1.0e-6,
				"dphi_0/dpsi against a central difference", radius, psi, keV );
		}
	}
}

BOOST_AUTO_TEST_CASE( dFdPsiAgreesWithACentralDifferenceForThreeSpecies )
{
	// The load-bearing test of the general path: dF/dpsi runs the chain rule
	// twice through the profiles AND through phi_0's two implicit derivatives.
	meq::RotatingSource const source = makeImpuritySource();

	for ( double h : { 1.0e-5, 1.0e-4 } )
	{
		for ( double radius : testRadii )
		{
			for ( double psi : testFluxes )
			{
				double const difference = ( source.f( radius, 0.0, psi + h ) - source.f( radius, 0.0, psi - h ) )/( 2.0*h );

				BOOST_CHECK_MESSAGE(
					std::fabs( source.dFdPsi( radius, 0.0, psi ) - difference )
						<= 1.0e-6*std::max( 1.0, std::fabs( difference ) ),
					"three species, h = " << h << ": dFdPsi disagrees with a central difference at R = "
					<< radius << ", psi = " << psi << ": analytic " << source.dFdPsi( radius, 0.0, psi )
					<< ", difference " << difference );
			}
		}
	}
}

BOOST_AUTO_TEST_CASE( theHeavyImpurityIsCentrifugallyEnriched )
{
	// The physics the whole exercise is for, and a control on the sign of phi_0:
	// carbon is twelve times as heavy as deuterium, so it is flung outboard
	// harder, and its density ratio to the bulk must RISE with radius. If phi_0
	// carried the wrong sign this would fall instead, and every algebraic test
	// above would still pass.
	meq::RotatingSource const source = makeImpuritySource();

	double const inboard = source.density( 1, 0.6, 0.0 )/source.density( 0, 0.6, 0.0 );
	double const outboard = source.density( 1, 1.4, 0.0 )/source.density( 0, 1.4, 0.0 );

	BOOST_CHECK_MESSAGE( outboard > 1.2*inboard,
		"the carbon-to-deuterium ratio is " << outboard << " outboard against " << inboard
		<< " inboard, so the impurity is not being centrifugally enriched" );
}

BOOST_AUTO_TEST_CASE( theClosedFormIsRefusedWhereItDoesNotApply )
{
	std::vector<meq::Species> const three = impuritySpecies();

	BOOST_CHECK_THROW( meq::RotatingSource( three, rotationProfile(), ggPrimeProfile(), referenceRadius,
		meq::vacuumPermeability, meq::RotatingSource::Closure::ClosedForm ), std::invalid_argument );

	// And Automatic picks the general path for three, the cheap one for two.
	meq::RotatingSource const auto3( three, rotationProfile(), ggPrimeProfile(), referenceRadius );
	meq::RotatingSource const auto2 = makeSource();
	BOOST_CHECK( auto3.closure() == meq::RotatingSource::Closure::RootFind );
	BOOST_CHECK( auto2.closure() == meq::RotatingSource::Closure::ClosedForm );
}

/*
 * FL-7: NORMALISED FLUX.
 *
 * The normalised source is a wrapper, so what needs testing is the SCALING and
 * nothing else -- and the scaling is asymmetric, one factor of 1/psi_ax in f()
 * and two in dFdPsi(). Getting that wrong gives a source that is perfectly
 * self-consistent and a Jacobian that is off by psi_ax, which is exactly the
 * defect a rate table cannot see.
 */

BOOST_AUTO_TEST_CASE( theNormalisedSourceCarriesOneFactorInFAndTwoInTheJacobian )
{
	for ( double psiAxis : { 0.35, 1.0, -0.8, 2.4 } )
	{
		meq::NormalisedRotatingSource const normalised( hydrogenicSpecies(), rotationProfile(),
			ggPrimeProfile(), referenceRadius, psiAxis );
		meq::RotatingSource const plain = makeSource();

		for ( double radius : testRadii )
		{
			for ( double psiN : testFluxes )
			{
				double const psi = psiN*psiAxis;

				checkRelative( normalised.f( radius, 0.0, psi ), plain.f( radius, 0.0, psiN )/psiAxis,
					1.0e-14, "the normalised F against the plain one rescaled", radius, psi );
				checkRelative( normalised.dFdPsi( radius, 0.0, psi ),
					plain.dFdPsi( radius, 0.0, psiN )/( psiAxis*psiAxis ),
					1.0e-14, "the normalised dF/dpsi against the plain one rescaled", radius, psi );
			}
		}
	}
}

BOOST_AUTO_TEST_CASE( theNormalisedJacobianAgreesWithACentralDifferenceInPhysicalFlux )
{
	// THE CHECK THAT DOES NOT ASSUME THE ANSWER. The test above compares the
	// wrapper against the formula it was written from, which cannot see a wrong
	// number of factors if the same wrong number is in both. This differences
	// f() in PHYSICAL psi, so it is the derivative the solver actually needs.
	for ( double psiAxis : { 0.35, 1.0, -0.8 } )
	{
		meq::NormalisedRotatingSource const source( hydrogenicSpecies(), rotationProfile(),
			ggPrimeProfile(), referenceRadius, psiAxis );

		double const h = 1.0e-5*std::fabs( psiAxis );

		for ( double radius : testRadii )
		{
			for ( double psiN : testFluxes )
			{
				double const psi = psiN*psiAxis;
				double const difference = ( source.f( radius, 0.0, psi + h ) - source.f( radius, 0.0, psi - h ) )/( 2.0*h );

				BOOST_CHECK_MESSAGE(
					std::fabs( source.dFdPsi( radius, 0.0, psi ) - difference )
						<= 1.0e-6*std::max( 1.0, std::fabs( difference ) ),
					"psi_ax = " << psiAxis << ": dFdPsi disagrees with a central difference in physical psi at R = "
					<< radius << ", psi = " << psi << ": analytic " << source.dFdPsi( radius, 0.0, psi )
					<< ", difference " << difference << " -- a wrong number of factors of psi_ax presents "
					"exactly like this" );
			}
		}
	}
}

BOOST_AUTO_TEST_CASE( theNormalisationCanBeMovedAndIsRefusedAtZero )
{
	meq::NormalisedRotatingSource source( hydrogenicSpecies(), rotationProfile(),
		ggPrimeProfile(), referenceRadius, 0.5 );

	BOOST_CHECK_EQUAL( source.normalisation(), 0.5 );

	// The solver moves it once per residual evaluation, so it has to answer for
	// the new value immediately.
	double const before = source.f( 1.2, 0.0, 0.25 );
	source.setNormalisation( 0.25 );
	BOOST_CHECK_EQUAL( source.normalisation(), 0.25 );
	BOOST_CHECK( source.f( 1.2, 0.0, 0.25 ) != before );

	BOOST_CHECK_THROW( source.setNormalisation( 0.0 ), std::invalid_argument );
	BOOST_CHECK_THROW( source.setNormalisation( std::nan( "" ) ), std::invalid_argument );

	// A refused value must not have been half-applied.
	BOOST_CHECK_EQUAL( source.normalisation(), 0.25 );
}

/**
 * `ConfineToPlasma` IS SILENTLY INERT ON A ROTATING SOURCE, AND THIS CASE
 * ASSERTS THE BEHAVIOUR THAT IS WANTED RATHER THAN THE ONE THAT IS THERE.
 *
 * meq::NormalisedSource::setPlasmaSupport() switches `F` off wherever `Psi <= 0`
 * and both normalised sources honour it at all three gates -- f(), dFdPsi() and
 * normalisationDerivatives(). This case is what keeps the rotating one there.
 *
 * IT IS REACHABLE FROM A FILE, which is what made the omission a defect rather
 * than an unimplemented corner. Config.cpp's Rotating branch lists
 * "ConfineToPlasma" among its accepted keys and calls readPlasmaSupport(), and
 * apps/meq.cpp calls setPlasmaSupport( true ) on whatever normalised source it
 * built -- so `[source] Type = "rotating"` with `Normalised = true` and
 * `ConfineToPlasma = true` parses, reports nothing, and would solve the
 * unconfined problem.
 *
 * WORSE THAN A MISSING ZERO, BECAUSE THE PROFILES CLAMP. A meq::SplineProfile
 * holds its endpoint value outside the knot range, so at `Psi < 0` the source
 * does not tail off: it sits at its plasma-edge value out to the wall, which is
 * a current density filling the vacuum region. And the element-level half of the
 * support DOES run -- GradShafranovSolver::plasmaComponentWanted() reads
 * plasmaSupport(), which is true -- so XP-1's flood fill masks whole elements
 * while the pointwise test inside them does nothing. Half a confinement is a
 * harder thing to recognise in an answer than none.
 *
 * THIS IS THE SHAPE CLAUDE.md ALREADY RECORDS UNDER "A RESERVED KEY IS ONLY
 * RESERVED ON THE PATHS THAT CALL THE REFUSAL": one source type reads a key
 * directly where another goes through a shared helper, and nothing tested the
 * pairing. **AND IT HAS NOW HAPPENED THREE TIMES ON THIS ONE PAIR OF SOURCE
 * TYPES** -- `ConfineToPlasma` here, `PlasmaCurrent` in
 * theCurrentScaleMultipliesTheRotatingPlasmaTerm, and
 * normalisationDerivatives() in theRotatingNormalisationDerivativesAreAnalytic.
 * Every accepted key on the rotating branch of Config.cpp is worth checking
 * against what meq::NormalisedRotatingSource actually reads; a fourth is more
 * likely than not.
 */
BOOST_AUTO_TEST_CASE( aConfinedRotatingSourceVanishesOutsideThePlasma )
{
	// psi_ax = 1, psi_bnd = 0, so Psi = psi and psi < 0 is outside the plasma.
	meq::NormalisedRotatingSource source( hydrogenicSpecies(), rotationProfile(),
		ggPrimeProfile(), referenceRadius, 1.0 );
	source.setPlasmaSupport( true );

	BOOST_TEST_REQUIRE( source.plasmaSupport() );
	// The base class's own test agrees the point is outside; it is only f() and
	// dFdPsi() that do not ask it.
	BOOST_TEST_REQUIRE( !source.insidePlasma( testRadii[ 0 ], 0.0, -0.5 ) );

	for ( double radius : testRadii )
		for ( double psi : { -0.5, -0.1, -1.0e-6 } )
		{
			BOOST_TEST( source.f( radius, 0.0, psi ) == 0.0,
			            "F = " << source.f( radius, 0.0, psi ) << " at R = " << radius
			            << ", Psi = " << psi << ", outside the plasma of a source "
			            "with setPlasmaSupport() on. ConfineToPlasma is reachable "
			            "for Type = \"rotating\" from a TOML file and is inert: "
			            "NormalisedRotatingSource::f() never calls insidePlasma()" );
			BOOST_TEST( source.dFdPsi( radius, 0.0, psi ) == 0.0,
			            "dF/dpsi = " << source.dFdPsi( radius, 0.0, psi ) << " at R = "
			            << radius << ", Psi = " << psi << ", outside the plasma; see above" );
		}

	// And the confinement must not have switched anything off INSIDE, which is
	// the half a bare "return 0" would break.
	for ( double radius : testRadii )
	{
		meq::NormalisedRotatingSource unconfined( hydrogenicSpecies(), rotationProfile(),
			ggPrimeProfile(), referenceRadius, 1.0 );
		BOOST_TEST( source.f( radius, 0.0, 0.5 ) == unconfined.f( radius, 0.0, 0.5 ) );
	}
}

BOOST_AUTO_TEST_SUITE_END()

/*
 * setCurrentScale() MULTIPLIES THE ROTATING PLASMA TERM, AND IT IS REACHABLE
 * FROM A FILE THAT DID NOT GET IT.
 *
 * meq::NormalisedSource::setCurrentScale()'s contract is "SCALE THE PLASMA
 * TERM, so that the TOTAL PLASMA CURRENT can be prescribed instead of the
 * profile amplitude", and the bordered Newton makes that scale an unknown.
 * meq::NormalisedMHDSource applies it in f() and in dFdPsi(); this class did
 * not apply it anywhere.
 *
 * WHAT THAT COSTS IS NOT A MISSING FEATURE BUT A JACOBIAN DESCRIBING A
 * DIFFERENT PROBLEM FROM THE RESIDUAL. With `F` independent of lambda:
 *
 *   - assemblePlasmaCurrent() returns an `int F/R` that does not respond to
 *     lambda at all, so the constraint cannot be satisfied by the unknown that
 *     exists to satisfy it;
 *   - assembleCurrentColumn() computes dR/dlambda as scaledF/lambda, which is
 *     not zero -- so the column is a fabrication;
 *   - cornerEntry( I, I ) returns ( int F/R )/lambda, likewise.
 *
 * Newton is then driven towards a current it can never deliver by a
 * sensitivity that does not exist.
 *
 * AND IT IS REACHABLE: Config.cpp's Rotating branch lists "PlasmaCurrent"
 * among its accepted keys, readPlasmaCurrent() refuses it only when
 * Normalised is false, and apps/meq.cpp calls setPlasmaCurrent() with no
 * source-type guard. THE SAME SHAPE CLAUDE.md RECORDS UNDER "A RESERVED KEY IS
 * ONLY RESERVED ON THE PATHS THAT CALL THE REFUSAL", and the second time this
 * exact pair of source types has produced it.
 *
 * THE ASSERTION IS TO A FEW ULP AND NOT TO EXACT EQUALITY, and the reason is
 * arithmetic rather than a concession. f() evaluates `( lambda H )/span` where
 * this case computes `lambda ( H/span )` -- the division sits between the two
 * multiplications, so they genuinely reassociate and differ in the last bit.
 * The grouping matches meq::NormalisedMHDSource::f()'s, which is worth more
 * than making one test's arithmetic exact. The property being asserted is
 * `f( lambda ) = lambda f( 1 )`, and four ulp is the right tolerance for it.
 */
BOOST_AUTO_TEST_CASE( theCurrentScaleMultipliesTheRotatingPlasmaTerm )
{
	meq::NormalisedRotatingSource source( hydrogenicSpecies(), rotationProfile(),
		ggPrimeProfile(), referenceRadius, 1.0 );
	source.setNormalisation( 0.9, -0.2 );

	BOOST_TEST_REQUIRE( source.currentScale() == 1.0,
	                    "the scale must default to one, or every solve that "
	                    "does not prescribe a current changes" );

	for ( double const lambda : { 0.5, 1.0, 2.5 } )
		for ( double const radius : testRadii )
			for ( double const psi : { 0.1, 0.4, 0.7 } )
			{
				source.setCurrentScale( 1.0 );
				double const bareF = source.f( radius, 0.0, psi );
				double const bareD = source.dFdPsi( radius, 0.0, psi );

				source.setCurrentScale( lambda );
				double const scaledFValue = source.f( radius, 0.0, psi );
				double const scaledDValue = source.dFdPsi( radius, 0.0, psi );

				auto closeEnough = []( double got, double want )
				{
					return std::abs( got - want )
					       <= 4.0*std::numeric_limits<double>::epsilon()
					          *std::max( std::abs( want ), 1.0e-300 );
				};

				BOOST_TEST( closeEnough( scaledFValue, lambda*bareF ),
				            "f() at lambda = " << lambda << ", R = " << radius
				            << ", psi = " << psi << " reads " << scaledFValue
				            << " against " << lambda*bareF );
				BOOST_TEST( closeEnough( scaledDValue, lambda*bareD ),
				            "dFdPsi() at lambda = " << lambda << ", R = " << radius
				            << ", psi = " << psi << " reads " << scaledDValue
				            << " against " << lambda*bareD );
			}

	// AND scaledF() IS THE SAME THING HERE, because this is the plasma source
	// itself and not a wrapper around one. A coil-augmented source overrides
	// these to forward to what it holds; an unwrapped one must not, or the
	// current border would count the scale twice.
	//
	// THIS PAIR IS EXACT, and legitimately so: the base implementation forwards
	// the same call with the same arguments, so there is no arithmetic between
	// them to reassociate.
	source.setCurrentScale( 2.5 );
	for ( double const radius : testRadii )
	{
		BOOST_TEST( source.scaledF( radius, 0.0, 0.4 ) == source.f( radius, 0.0, 0.4 ) );
		BOOST_TEST( source.scaledDFdPsi( radius, 0.0, 0.4 )
		            == source.dFdPsi( radius, 0.0, 0.4 ) );
	}
}

/*
 * THE ROTATING SOURCE'S NORMALISATION DERIVATIVES ARE ANALYTIC, AND THE SPECIES
 * CLOSURE NEVER HAS TO BE DIFFERENTIATED AGAIN TO GET THEM.
 *
 * The mirror of SourceTests.cpp's the_normalisation_derivatives_are_analytic,
 * and it exists for a stronger reason than symmetry: without an override the
 * base returns false, and TWO callers degrade rather than one.
 * solveWithNormalisation()'s psi_ax column falls back to a central difference
 * -- which perturbs the normalisation, which moves the plasma edge, so it
 * straddles a kink rather than measuring a derivative -- and
 * assembleCurrentNormalisationCorner() used to fall back to a silent ZERO,
 * which CLAUDE_FB.md measures turning the current-bordered Newton from
 * quadratic into a clean geometric contraction of about 0.8 a step.
 *
 * WHY IT IS SIX LINES AND NOT A DERIVATION. This class is a pure wrapper:
 * F = S H( Psi )/span with H := inner.f, and inner.dFdPsi IS dH/dPsi. So the
 * normalisation enters only through the argument and the overall 1/span, and
 * meq::NormalisedMHDSource's algebra applies verbatim. The header said this
 * was not written because "the rotating closure's version" would have to be;
 * it does not.
 *
 * Checked against a difference of f() ITSELF, one level below any assembly, so
 * a failure here is the formula and can be nothing else -- and Richardson
 * extrapolated, because a plain central difference would floor the comparison
 * at its own O( h^2 ) truncation rather than at the derivative.
 */
BOOST_AUTO_TEST_CASE( theRotatingNormalisationDerivativesAreAnalytic )
{
	double const z = 0.0;

	meq::NormalisedRotatingSource source( hydrogenicSpecies(), rotationProfile(),
		ggPrimeProfile(), referenceRadius, 1.0 );

	// A SCALE OTHER THAN ONE, deliberately: the scale multiplies F and so
	// multiplies both derivatives, and at lambda = 1 a formula that dropped it
	// would pass.
	source.setCurrentScale( 1.7 );

	auto difference = [ & ]( double radius, double psi, bool axis )
	{
		auto at = [ & ]( double step )
		{
			source.setNormalisation( axis ? 0.9 + step : 0.9,
			                         axis ? -0.2 : -0.2 + step );
			double const plus = source.f( radius, z, psi );
			source.setNormalisation( axis ? 0.9 - step : 0.9,
			                         axis ? -0.2 : -0.2 - step );
			double const minus = source.f( radius, z, psi );
			return ( plus - minus )/( 2.0*step );
		};
		double const h = 1.0e-4;
		return ( 4.0*at( h/2.0 ) - at( h ) )/3.0;
	};

	std::printf( "\n  THE ROTATING BORDER COLUMNS IN CLOSED FORM\n" );
	std::printf( "    %6s %8s %16s %16s %12s %16s %16s %12s\n",
	             "R", "psi", "dF/dpsi_ax", "difference", "rel",
	             "dF/dpsi_bnd", "difference", "rel" );

	for ( double const radius : testRadii )
		for ( double const psi : { 0.1, 0.4, 0.7 } )
		{
			source.setNormalisation( 0.9, -0.2 );
			double analyticAxis = 0.0, analyticBoundary = 0.0;
			BOOST_TEST_REQUIRE( source.normalisationDerivatives(
				radius, z, psi, analyticAxis, analyticBoundary ),
				"NormalisedRotatingSource must supply its normalisation "
				"derivatives" );

			double const numericAxis = difference( radius, psi, true );
			double const numericBoundary = difference( radius, psi, false );

			double const relA = std::abs( analyticAxis - numericAxis )
			                    /std::max( std::abs( numericAxis ), 1.0 );
			double const relB = std::abs( analyticBoundary - numericBoundary )
			                    /std::max( std::abs( numericBoundary ), 1.0 );

			std::printf( "    %6.2f %8.2f %16.8e %16.8e %12.2e %16.8e %16.8e "
			             "%12.2e\n", radius, psi, analyticAxis, numericAxis, relA,
			             analyticBoundary, numericBoundary, relB );

			BOOST_TEST( relA < 1.0e-8,
			            "dF/dpsi_ax at R = " << radius << ", psi = " << psi
			            << " reads " << analyticAxis << " against a differenced "
			            << numericAxis );
			BOOST_TEST( relB < 1.0e-8,
			            "dF/dpsi_bnd at R = " << radius << ", psi = " << psi
			            << " reads " << analyticBoundary
			            << " against a differenced " << numericBoundary );
		}
	std::fflush( stdout );

	/*
	 * AND OUTSIDE THE PLASMA BOTH ARE EXACTLY ZERO, which is the half a
	 * difference structurally cannot reproduce: perturbing the normalisation
	 * moves the edge, so a point just outside is INSIDE for one of the two
	 * evaluations and the difference reports a spurious F/( 2 step ) -- which
	 * diverges as the step is refined rather than converging to anything.
	 */
	source.setPlasmaSupport( true );
	source.setNormalisation( 0.9, -0.2 );
	double outsideAxis = 1.0, outsideBoundary = 1.0;
	BOOST_TEST_REQUIRE( source.normalisationDerivatives( testRadii[ 0 ], z, -0.5,
	                                                     outsideAxis,
	                                                     outsideBoundary ) );
	BOOST_TEST( outsideAxis == 0.0 );
	BOOST_TEST( outsideBoundary == 0.0 );
}

/*
 * A MIRROR'S SUPPORT IS A BOX IN ( Psi, z ), AND THE SOURCE MUST VANISH ON ALL
 * THREE OF ITS SIDES. MIRROR-PLAN.md MR-1.
 *
 * A mirror has no magnetic axis. Its plasma is an ANNULUS labelled by two
 * prescribed flux surfaces -- psi_ax the INNER one and psi_bnd the OUTER -- so
 * the span psi_ax - psi_bnd is NEGATIVE, and the one-sided `Psi > 0` a tokamak
 * needs excludes only the radially-outer side. Two more are owed:
 *
 *   inside the inner surface   Psi > 1, toward R = 0, where there is no plasma
 *                              at all -- and the profiles CLAMP, so ungated the
 *                              source sits at its inner-edge value down to the
 *                              axis;
 *   beyond the throat          |z| > z_edge, where the tube flares into the
 *                              expander and the same Psi recurs. M-168
 *                              measured geq's flux-only mask putting the
 *                              end-plane current five orders above the
 *                              midplane's.
 *
 * ALL THREE QUANTITIES THE SOLVER ASSEMBLES ARE CHECKED AT EVERY PROBE -- f(),
 * dFdPsi() and normalisationDerivatives() -- because CLAUDE_FLOW.md records
 * ConfineToPlasma reaching one of a pair of evaluations and not the other, and
 * half a confinement is harder to see in an answer than none.
 *
 * EVERY PROBE IS FIRST SHOWN TO CARRY A NON-ZERO UNCONFINED SOURCE, so each
 * zero asserted below is one the gate produced and not one the profiles
 * happened to give. Without that this case could not fail.
 *
 * The edges are CLOSED: Psi = 1 (the inner label itself) and |z| = z_edge are
 * inside, which is what lets GradShafranovSolver::checkAxialCutoff() evaluate
 * the source AT the cutoff and read the interior limit of the jump.
 */
BOOST_AUTO_TEST_CASE( aMirrorSourceVanishesOnAllThreeSidesOfItsBox )
{
	double const innerLabel = 0.2;   // psi_ax: the inner surface
	double const outerLabel = 1.0;   // psi_bnd: the outer surface
	double const zEdge = 0.8;

	auto mirror = [ & ]()
	{
		meq::NormalisedRotatingSource source( hydrogenicSpecies(), rotationProfile(),
			ggPrimeProfile(), referenceRadius, 1.0 );
		source.setNormalisation( innerLabel, outerLabel );
		return source;
	};

	meq::NormalisedRotatingSource boxed = mirror();
	meq::NormalisedRotatingSource const unconfined = mirror();

	// THE SPAN IS NEGATIVE, which is the orientation the bound has to survive.
	BOOST_TEST_REQUIRE( innerLabel - outerLabel < 0.0 );

	boxed.setPlasmaSupportShape( meq::PlasmaSupportShape::Box, zEdge );
	BOOST_TEST_REQUIRE( ( boxed.plasmaSupportShape() == meq::PlasmaSupportShape::Box ) );
	BOOST_TEST_REQUIRE( boxed.plasmaSupportAxialEdge() == zEdge );
	BOOST_TEST_REQUIRE( boxed.plasmaSupport(),
	                    "asking for a box must switch the support on" );

	// The three sides, as ( psi, z ) probes. Psi = ( psi - 1 )/( -0.8 ).
	struct Probe
	{
		double psi;
		double z;
	};
	std::vector<Probe> const innerSide{    // Psi > 1
		{ 0.19, 0.0 }, { 0.1, 0.0 }, { 0.0, 0.3 }, { -0.1, -0.5 } };
	std::vector<Probe> const outerSide{    // Psi < 0
		{ 1.01, 0.0 }, { 1.2, 0.0 }, { 1.5, 0.4 }, { 2.0, -0.7 } };
	std::vector<Probe> const beyondThroat{ // 0 < Psi <= 1, |z| > z_edge
		{ 0.6, 0.81 }, { 0.6, -0.81 }, { 0.3, 1.5 }, { 0.9, -3.0 } };
	std::vector<Probe> const insideBox{    // 0 < Psi <= 1, |z| <= z_edge
		{ 0.3, 0.0 }, { 0.6, 0.5 }, { 0.9, -0.5 },
		{ innerLabel, 0.0 },               // Psi = 1 exactly: closed
		{ 0.6, zEdge }, { 0.6, -zEdge } }; // |z| = z_edge: closed

	auto vanishes = [ & ]( std::vector<Probe> const &side, char const *name )
	{
		for ( double const radius : testRadii )
			for ( Probe const &p : side )
			{
				BOOST_TEST_REQUIRE( unconfined.f( radius, p.z, p.psi ) != 0.0,
				                    "the unconfined source is already zero at R = "
				                    << radius << ", psi = " << p.psi << ", z = "
				                    << p.z << ", so a zero there proves nothing" );

				BOOST_TEST( boxed.f( radius, p.z, p.psi ) == 0.0,
				            "F = " << boxed.f( radius, p.z, p.psi ) << " " << name
				            << " at R = " << radius << ", psi = " << p.psi
				            << ", z = " << p.z );
				BOOST_TEST( boxed.dFdPsi( radius, p.z, p.psi ) == 0.0,
				            "dF/dpsi = " << boxed.dFdPsi( radius, p.z, p.psi ) << " "
				            << name << " at R = " << radius << ", psi = " << p.psi
				            << ", z = " << p.z );

				double againstAxis = 1.0, againstBoundary = 1.0;
				BOOST_TEST_REQUIRE( boxed.normalisationDerivatives(
					radius, p.z, p.psi, againstAxis, againstBoundary ) );
				BOOST_TEST( againstAxis == 0.0,
				            "dF/dpsi_ax = " << againstAxis << " " << name
				            << " at R = " << radius << ", psi = " << p.psi
				            << ", z = " << p.z );
				BOOST_TEST( againstBoundary == 0.0,
				            "dF/dpsi_bnd = " << againstBoundary << " " << name
				            << " at R = " << radius << ", psi = " << p.psi
				            << ", z = " << p.z );
			}
	};

	vanishes( innerSide, "inside the inner surface (Psi > 1)" );
	vanishes( outerSide, "outside the outer surface (Psi < 0)" );
	vanishes( beyondThroat, "beyond the throat (|z| > z_edge)" );

	// AND NOTHING INSIDE IS SWITCHED OFF OR CHANGED, which is the half a bare
	// "return 0" would break: every quantity is the unconfined one, to the bit.
	for ( double const radius : testRadii )
		for ( Probe const &p : insideBox )
		{
			double const want = unconfined.f( radius, p.z, p.psi );
			BOOST_TEST_REQUIRE( want != 0.0 );
			BOOST_TEST( boxed.f( radius, p.z, p.psi ) == want,
			            "F inside the box at R = " << radius << ", psi = " << p.psi
			            << ", z = " << p.z << " reads "
			            << boxed.f( radius, p.z, p.psi ) << " against " << want );
			BOOST_TEST( boxed.dFdPsi( radius, p.z, p.psi )
			            == unconfined.dFdPsi( radius, p.z, p.psi ) );

			double boxAxis = 0.0, boxBoundary = 0.0, freeAxis = 0.0, freeBoundary = 0.0;
			BOOST_TEST_REQUIRE( boxed.normalisationDerivatives(
				radius, p.z, p.psi, boxAxis, boxBoundary ) );
			BOOST_TEST_REQUIRE( unconfined.normalisationDerivatives(
				radius, p.z, p.psi, freeAxis, freeBoundary ) );
			BOOST_TEST( boxAxis == freeAxis );
			BOOST_TEST( boxBoundary == freeBoundary );
			BOOST_TEST( boxAxis != 0.0 );
		}

	// setPlasmaSupport( false ) still switches the whole gate off, which is how
	// a caller runs the same source unconfined.
	boxed.setPlasmaSupport( false );
	BOOST_TEST( boxed.f( testRadii[ 0 ], 0.0, 0.1 )
	            == unconfined.f( testRadii[ 0 ], 0.0, 0.1 ) );
}

/*
 * THE TOKAMAK'S SUPPORT IS UNCHANGED BY THE BOX EXISTING, AND IN PARTICULAR IT
 * DOES NOT IMPOSE Psi <= 1.
 *
 * In a tokamak psi_ax is the EXTREMUM of psi, so Psi <= 1 holds by
 * construction -- except that psi_h can overshoot the LOCATED psi_ax by a
 * discretisation error, and a Psi <= 1 bound there would cut a hole round the
 * magnetic axis. Nor is there any |z| bound. So under the default shape a point
 * at Psi = 1.05, or at any z, must carry exactly the unconfined source, and the
 * default must BE the one-sided shape without anybody asking for it.
 */
BOOST_AUTO_TEST_CASE( theTokamakSupportIsOneSidedAndUnchanged )
{
	meq::NormalisedRotatingSource confined( hydrogenicSpecies(), rotationProfile(),
		ggPrimeProfile(), referenceRadius, 1.0 );
	meq::NormalisedRotatingSource const unconfined( hydrogenicSpecies(),
		rotationProfile(), ggPrimeProfile(), referenceRadius, 1.0 );
	confined.setPlasmaSupport( true );

	BOOST_TEST_REQUIRE( ( confined.plasmaSupportShape()
	                      == meq::PlasmaSupportShape::OneSided ),
	                    "the default support shape must be the tokamak's" );
	BOOST_TEST( std::isinf( confined.plasmaSupportAxialEdge() ) );

	for ( double const radius : testRadii )
		for ( double const z : { 0.0, 2.0, -50.0 } )
			for ( double const psi : { 0.3, 1.0, 1.05, 1.5 } )
			{
				BOOST_TEST_REQUIRE( unconfined.f( radius, z, psi ) != 0.0 );
				BOOST_TEST( confined.f( radius, z, psi ) == unconfined.f( radius, z, psi ),
				            "a one-sided support changed F at R = " << radius
				            << ", z = " << z << ", Psi = " << psi );
				BOOST_TEST( confined.dFdPsi( radius, z, psi )
				            == unconfined.dFdPsi( radius, z, psi ) );
			}

	// Outside, Psi <= 0, it is zero -- the one bound a tokamak has.
	BOOST_TEST( confined.f( testRadii[ 0 ], 0.0, -0.1 ) == 0.0 );

	// AND THE BOX IS REVERSIBLE: asking for the one-sided shape back puts the
	// source exactly where it started.
	confined.setPlasmaSupportShape( meq::PlasmaSupportShape::Box, 0.5 );
	BOOST_TEST( confined.f( testRadii[ 0 ], 0.0, 1.05 ) == 0.0 );
	confined.setPlasmaSupportShape( meq::PlasmaSupportShape::OneSided );
	BOOST_TEST( ( confined.plasmaSupportShape() == meq::PlasmaSupportShape::OneSided ) );
	BOOST_TEST( confined.f( testRadii[ 0 ], 2.0, 1.05 )
	            == unconfined.f( testRadii[ 0 ], 2.0, 1.05 ) );
}

/*
 * THE AXIAL EDGE IS PRESCRIBED, AND A BOX WITHOUT ONE IS REFUSED.
 *
 * An infinite edge would be a mirror support with no axial bound -- geq's
 * flux-only mask, which M-168 measured enrolling the expander as plasma -- so
 * it takes more than a default to ask for it. A refused call must leave the
 * source as it was.
 */
BOOST_AUTO_TEST_CASE( aBoxWithoutAFinitePositiveEdgeIsRefused )
{
	meq::NormalisedRotatingSource source( hydrogenicSpecies(), rotationProfile(),
		ggPrimeProfile(), referenceRadius, 1.0 );

	for ( double const bad : { 0.0, -0.5, std::numeric_limits<double>::infinity(),
	                           std::numeric_limits<double>::quiet_NaN() } )
	{
		BOOST_CHECK_THROW( source.setPlasmaSupportShape( meq::PlasmaSupportShape::Box,
		                                                 bad ),
		                   std::invalid_argument );
		BOOST_TEST( ( source.plasmaSupportShape() == meq::PlasmaSupportShape::OneSided ) );
		BOOST_TEST( source.plasmaSupport() == false,
		            "a refused box must not have switched the support on" );
	}
	BOOST_CHECK_THROW( source.setPlasmaSupportShape( meq::PlasmaSupportShape::Box ),
	                   std::invalid_argument );
}

// ---------------------------------------------------------------------------
// THE MIDPLANE GAUGE. MIRROR-PLAN.md section 5: phi_0 pinned where each flux
// surface crosses z = 0, R_0( psi ), rather than at a constant R_ref.
// ---------------------------------------------------------------------------

namespace
{
	/*
	 * A CROSSING IN CLOSED FORM: R_0^2 = a + b psi + c psi^2, falling from 1.4
	 * at psi = 0 to 0.6 at psi = 1 across testRadii, as a mirror's does from
	 * its outer label to its inner. QUADRATIC IN R_0^2 ON PURPOSE: the midplane
	 * gauge's half-delta is ( R^2 - R_0^2 )/2, and with R_0^2 linear in psi its
	 * SECOND derivative would vanish identically and the h'' term of the
	 * Jacobian would be checked against nothing.
	 */
	double const crossingA = 1.96;
	double const crossingB = -1.2;
	double const crossingC = -0.4;

	double crossingSquared( double psi ) { return crossingA + crossingB*psi + crossingC*psi*psi; }
	double crossingSquaredPrime( double psi ) { return crossingB + 2.0*crossingC*psi; }

	double crossing( double psi ) { return std::sqrt( crossingSquared( psi ) ); }
	double crossingPrime( double psi ) { return crossingSquaredPrime( psi )/( 2.0*crossing( psi ) ); }
	double crossingDoublePrime( double psi )
	{
		double const r = crossing( psi );
		double const q = crossingSquaredPrime( psi );
		return crossingC/r - q*q/( 4.0*r*r*r );
	}

	std::shared_ptr<meq::Profile const> crossingProfile()
	{
		return analytic( crossing, crossingPrime, crossingDoublePrime );
	}

	/// The control: R_0 right, its slope and curvature dropped -- the
	/// midplane gauge as it would be with the new Jacobian terms missing.
	std::shared_ptr<meq::Profile const> crossingWithoutSlope()
	{
		return analytic( crossing, []( double ) { return 0.0; }, []( double ) { return 0.0; } );
	}

	meq::RotatingSource midplaneSource( std::vector<meq::Species> species, double omegaScale,
	                                    std::shared_ptr<meq::Profile const> table )
	{
		meq::RotatingSource source( std::move( species ),
			omegaScale == 0.0 ? nullptr : rotationProfile( omegaScale ),
			ggPrimeProfile(), 0.0, meq::vacuumPermeability,
			meq::RotatingSource::Closure::Automatic, meq::RotatingSource::Gauge::Midplane );
		source.setMidplaneCrossing( std::move( table ) );
		return source;
	}

	/// The midplane of a mirror in closed form: psi( R, 0 ) = B R^2 ( 1 + k R )/2,
	/// which rises monotonically from the axis, normalised between an inner
	/// label at 0.1 and an outer at 0.35 so that Psi = 1 inside and 0 outside,
	/// as MR-4's.
	double const mirrorField = 4.0;
	double const mirrorBend = 0.3;
	double const mirrorInner = 0.1;
	double const mirrorOuter = 0.35;

	double mirrorPsi( double radius )
	{
		return mirrorField*radius*radius*( 1.0 + mirrorBend*radius )/2.0;
	}

	double mirrorPsiPrime( double radius )
	{
		return mirrorField*radius*( 1.0 + 1.5*mirrorBend*radius );
	}

	meq::MidplaneField mirrorNormalised()
	{
		return []( double radius, double &x, double &dXdR )
		{
			double const span = mirrorPsi( mirrorInner ) - mirrorPsi( mirrorOuter );
			x = ( mirrorPsi( radius ) - mirrorPsi( mirrorOuter ) )/span;
			dXdR = mirrorPsiPrime( radius )/span;
			return true;
		};
	}

	/// R with Psi( R ) = x, by bisection on the closed form -- the
	/// independent inverse a tabulated crossing is held to.
	double mirrorCrossing( double x )
	{
		double lo = 0.05, hi = 0.5;
		double const span = mirrorPsi( mirrorInner ) - mirrorPsi( mirrorOuter );
		for ( int i = 0; i < 200; ++i )
		{
			double const mid = 0.5*( lo + hi );
			double const value = ( mirrorPsi( mid ) - mirrorPsi( mirrorOuter ) )/span;
			( value > x ? lo : hi ) = mid;
		}
		return 0.5*( lo + hi );
	}

	/// The message a refusal carries, or "" when there was none.
	template<typename Call>
	std::string refusal( Call call )
	{
		try
		{
			call();
		}
		catch ( std::exception const &e )
		{
			return e.what();
		}
		return "";
	}
}

BOOST_AUTO_TEST_SUITE( midplane_gauge )

/*
 * THE GAUGE IS WHAT IT SAYS: phi_0 vanishes and every density is its table's
 * value WHERE THE SURFACE CROSSES THE MIDPLANE, and nowhere else in general.
 */
BOOST_AUTO_TEST_CASE( theMidplaneGaugePinsThePotentialAtEachCrossing )
{
	for ( bool three : { false, true } )
	{
		std::vector<meq::Species> const species = three ? impuritySpecies() : hydrogenicSpecies();
		meq::RotatingSource const source = midplaneSource( species, 1.0, crossingProfile() );
		BOOST_TEST( ( source.gauge() == meq::RotatingSource::Gauge::Midplane ) );
		BOOST_TEST( source.referenceRadius() == 0.0 );

		for ( double psi : testFluxes )
		{
			double const r0 = crossing( psi );
			checkClose( source.potential( r0, psi )/keV, 0.0, 1.0e-12,
				three ? "e phi_0 at the crossing, three species" : "e phi_0 at the crossing", psi );
			for ( std::size_t s = 0; s < species.size(); ++s )
				checkRelative( source.density( s, r0, psi ), ( *species[ s ].density )( psi ), 1.0e-13,
					"n_s at the crossing against n_s0", r0, psi );
			// Off the crossing it is not pinned: a source that ignored R_0 and
			// pinned nothing would pass the two checks above at phi_0 = 0
			// everywhere, and fails this.
			BOOST_TEST( std::fabs( source.potential( r0 + 0.2, psi ) ) > 1.0e-6*keV,
				"phi_0 vanishes off the crossing at psi = " << psi );
		}
	}
}

/*
 * ACCEPTANCE 1: THE NEW JACOBIAN TERM AGAINST A CENTRAL DIFFERENCE.
 *
 * Under the midplane gauge dA_s/dpsi gains -m_s omega^2 R_0 R_0'/T_s, and
 * d2A_s/dpsi2 the term below it, -( R_0'^2 + R_0 R_0'' ) times the same. F is
 * mu0 R^2 dp/dpsi, so the first enters the RESIDUAL and the second the
 * Jacobian. Both closures, over a Mach sweep, against a central difference of
 * pressure(), f() and potential() -- exactly as the fixed-radius gauge's FD
 * checks above, and at their tolerances.
 *
 * AND THE CHECK IS SHOWN TO SEE THE TERM: the same source handed R_0 with its
 * slope and curvature zeroed -- the gauge with the new terms missing -- is
 * asserted to FAIL the pressure comparison by orders of magnitude.
 */
BOOST_AUTO_TEST_CASE( theMidplaneJacobianAgreesWithACentralDifference )
{
	double worstControl = 0.0;
	std::printf( "\n  THE MIDPLANE GAUGE'S DERIVATIVES AGAINST CENTRAL DIFFERENCES\n" );

	for ( bool three : { false, true } )
		for ( double omegaScale : { 0.5, 1.0, 2.0 } )
		{
			std::vector<meq::Species> const species = three ? impuritySpecies() : hydrogenicSpecies();
			meq::RotatingSource const source = midplaneSource( species, omegaScale, crossingProfile() );
			meq::RotatingSource const control = midplaneSource( species, omegaScale, crossingWithoutSlope() );
			BOOST_TEST_REQUIRE( ( source.closure() == ( three ? meq::RotatingSource::Closure::RootFind
			                                                  : meq::RotatingSource::Closure::ClosedForm ) ) );

			double worstHere = 0.0;
			for ( double radius : testRadii )
				for ( double psi : testFluxes )
				{
					double const hp = 1.0e-6;
					double const dp = ( source.pressure( radius, psi + hp ) - source.pressure( radius, psi - hp ) )/( 2.0*hp );
					checkRelative( source.dPressureDPsi( radius, psi ), dp, 1.0e-7,
						"midplane gauge, dp/dpsi against a central difference", radius, psi );

					double const dphi = ( source.potential( radius, psi + hp ) - source.potential( radius, psi - hp ) )/( 2.0*hp );
					checkRelative( source.dPotentialDPsi( radius, psi ), dphi, 1.0e-6,
						"midplane gauge, dphi_0/dpsi against a central difference", radius, psi, keV );

					double const dpControl = ( control.pressure( radius, psi + hp ) - control.pressure( radius, psi - hp ) )/( 2.0*hp );
					worstControl = std::max( worstControl, std::fabs( control.dPressureDPsi( radius, psi ) - dpControl )
						/std::max( std::fabs( dpControl ), 1.0e-300 ) );

					for ( double h : { 1.0e-5, 1.0e-4 } )
					{
						double const difference = ( source.f( radius, 0.0, psi + h ) - source.f( radius, 0.0, psi - h ) )/( 2.0*h );
						double const analyticValue = source.dFdPsi( radius, 0.0, psi );
						double const error = std::fabs( analyticValue - difference )/std::max( 1.0, std::fabs( difference ) );
						worstHere = std::max( worstHere, error );
						BOOST_CHECK_MESSAGE( error <= 1.0e-6,
							( three ? "three species" : "two species" ) << ", omega scale " << omegaScale
							<< ", h = " << h << ": the midplane gauge's dFdPsi disagrees with a central "
							"difference at R = " << radius << ", psi = " << psi << ": analytic "
							<< analyticValue << ", difference " << difference );
					}
				}
			std::printf( "    %-13s omega x %.1f: worst | dF/dpsi - difference |/max( 1, |difference| ) %.3e\n",
			             three ? "three species" : "two species", omegaScale, worstHere );
		}

	std::printf( "    control, R_0' and R_0'' dropped: worst relative dp/dpsi error %.3e\n", worstControl );
	std::fflush( stdout );
	BOOST_TEST( worstControl > 1.0e-2,
		"dropping R_0' and R_0'' moves dp/dpsi by only " << worstControl << " relative, so the "
		"central-difference check above cannot see the midplane gauge's Jacobian term" );
}

/*
 * THE TWO GAUGES DESCRIBE ONE PLASMA. Two species with constant T and omega,
 * so the shared exponent is C( R^2 - R_*^2 )/2 with C constant and the gauge
 * transfer is closed form:
 *
 *     n_s0^mid( psi ) = n_s0^fixed( psi ) exp( C ( R_0( psi )^2 - R_ref^2 )/2 ),
 *
 * which with R_0^2 quadratic in psi has both derivatives in closed form too. So
 * every PHYSICAL quantity -- n_s, p, F and dF/dpsi -- must agree between the two
 * at every ( R, psi ), and phi_0 must differ by a function of psi alone, which
 * is what a gauge is.
 *
 * It tests the R_0 chain rule against a second route that does not go through
 * it: the fixed-radius source has no R_0 anywhere.
 */
BOOST_AUTO_TEST_CASE( theTwoGaugesDescribeOnePlasma )
{
	double const tIon = 1.8*keV;
	double const tElectron = 1.3*keV;
	double const omega = 1.3*sonicOmega;
	double const c = omega*omega*( electronMass + deuteronMass )/( tElectron + tIon );

	auto transfer = [ c ]( double psi )
	{
		return std::exp( c*( crossingSquared( psi ) - referenceRadius*referenceRadius )/2.0 );
	};
	auto midplaneDensity = analytic(
		[ transfer ]( double psi ) { return ionDensity( psi )*transfer( psi ); },
		[ transfer, c ]( double psi )
		{
			double const g = c*crossingSquaredPrime( psi )/2.0;
			return ( ionDensityPrime( psi ) + ionDensity( psi )*g )*transfer( psi );
		},
		[ transfer, c ]( double psi )
		{
			double const g = c*crossingSquaredPrime( psi )/2.0;
			double const gPrime = c*crossingC;
			return ( ionDensityDoublePrime( psi ) + 2.0*ionDensityPrime( psi )*g
			         + ionDensity( psi )*( g*g + gPrime ) )*transfer( psi );
		} );

	auto pair = [ & ]( std::shared_ptr<meq::Profile const> ionDensityTable )
	{
		std::vector<meq::Species> species( 2 );
		species[ 0 ].mass = deuteronMass;
		species[ 0 ].charge = 1.0;
		species[ 0 ].temperature = std::make_shared<meq::ConstantProfile const>( tIon );
		species[ 0 ].density = std::move( ionDensityTable );
		species[ 1 ].mass = electronMass;
		species[ 1 ].charge = -1.0;
		species[ 1 ].temperature = std::make_shared<meq::ConstantProfile const>( tElectron );
		species[ 1 ].density = meq::neutralisingDensity( species, 1 );
		return species;
	};

	auto rotationConstant = std::make_shared<meq::ConstantProfile const>( omega );
	meq::RotatingSource const fixed( pair( ionDensityProfile() ), rotationConstant, ggPrimeProfile(),
	                                 referenceRadius );
	meq::RotatingSource midplane( pair( midplaneDensity ), rotationConstant, ggPrimeProfile(), 0.0,
	                              meq::vacuumPermeability, meq::RotatingSource::Closure::Automatic,
	                              meq::RotatingSource::Gauge::Midplane );
	midplane.setMidplaneCrossing( crossingProfile() );

	double worst = 0.0;
	for ( double psi : testFluxes )
	{
		double const offset = fixed.potential( testRadii[ 0 ], psi ) - midplane.potential( testRadii[ 0 ], psi );
		for ( double radius : testRadii )
		{
			for ( std::size_t s = 0; s < 2; ++s )
				checkRelative( midplane.density( s, radius, psi ), fixed.density( s, radius, psi ), 1.0e-13,
					"n_s, midplane against fixed-radius gauge", radius, psi );
			checkRelative( midplane.pressure( radius, psi ), fixed.pressure( radius, psi ), 1.0e-13,
				"p, midplane against fixed-radius gauge", radius, psi );
			checkRelative( midplane.f( radius, 0.0, psi ), fixed.f( radius, 0.0, psi ), 1.0e-12,
				"F, midplane against fixed-radius gauge", radius, psi );
			checkRelative( midplane.dFdPsi( radius, 0.0, psi ), fixed.dFdPsi( radius, 0.0, psi ), 1.0e-11,
				"dF/dpsi, midplane against fixed-radius gauge", radius, psi );
			checkRelative( fixed.potential( radius, psi ) - midplane.potential( radius, psi ), offset, 1.0e-12,
				"phi_0's gauge offset, which must not depend on R", radius, psi, keV );
			worst = std::max( worst, std::fabs( midplane.dFdPsi( radius, 0.0, psi ) - fixed.dFdPsi( radius, 0.0, psi ) )
				/std::max( std::fabs( fixed.dFdPsi( radius, 0.0, psi ) ), 1.0e-300 ) );
		}
	}
	std::printf( "\n  ONE PLASMA IN TWO GAUGES: worst relative dF/dpsi difference %.3e\n", worst );
}

/*
 * tabulateMidplaneCrossing() ON A MIRROR'S MIDPLANE IN CLOSED FORM: every knot
 * on the root to round-off, its slope the inverse function theorem's, and the
 * Hermite interpolant between knots against an independent bisection.
 */
BOOST_AUTO_TEST_CASE( theCrossingIsTabulatedOnAMirrorsMidplane )
{
	std::shared_ptr<meq::SplineProfile const> const table =
		meq::tabulateMidplaneCrossing( mirrorNormalised(), 0.0, 3.0, 0.0, 1.0 );
	BOOST_TEST_REQUIRE( table->knots().size() == 1025u );
	BOOST_TEST( table->uniformKnots() );

	double const span = mirrorPsi( mirrorInner ) - mirrorPsi( mirrorOuter );
	double worstKnot = 0.0, worstSlope = 0.0, worstBetween = 0.0, worstBetweenSlope = 0.0;
	for ( meq::Knot const &k : table->knots() )
	{
		double const x = ( mirrorPsi( k.value ) - mirrorPsi( mirrorOuter ) )/span;
		worstKnot = std::max( worstKnot, std::fabs( x - k.psi ) );
		double const slope = span/mirrorPsiPrime( k.value );
		worstSlope = std::max( worstSlope, std::fabs( k.derivative - slope )/std::fabs( slope ) );
	}
	for ( std::size_t i = 0; i + 1 < table->knots().size(); ++i )
	{
		double const x = 0.5*( table->knots()[ i ].psi + table->knots()[ i + 1 ].psi );
		double const r = mirrorCrossing( x );
		worstBetween = std::max( worstBetween, std::fabs( ( *table )( x ) - r ) );
		double const slope = span/mirrorPsiPrime( r );
		worstBetweenSlope = std::max( worstBetweenSlope, std::fabs( table->prime( x ) - slope )/std::fabs( slope ) );
	}
	std::printf( "\n  THE MIDPLANE CROSSING, 1024 INTERVALS: knot residual in Psi %.3e, knot slope %.3e "
	             "relative; between knots R_0 %.3e m and its slope %.3e relative\n",
	             worstKnot, worstSlope, worstBetween, worstBetweenSlope );

	BOOST_TEST( worstKnot < 1.0e-14 );
	BOOST_TEST( worstSlope < 1.0e-13 );
	BOOST_TEST( worstBetween < 1.0e-11 );
	BOOST_TEST( worstBetweenSlope < 1.0e-7 );
	BOOST_TEST( ( *table )( 1.0 ) == mirrorInner, boost::test_tools::tolerance( 1.0e-14 ) );
	BOOST_TEST( ( *table )( 0.0 ) == mirrorOuter, boost::test_tools::tolerance( 1.0e-14 ) );
}

/*
 * THE PRECONDITION IS CHECKED, AND EACH WAY IT FAILS IS REFUSED BY NAME.
 *
 * A tokamak's midplane, Psi = 1 - ( ( R - 1 )/0.4 )^2: every surface crosses
 * z = 0 twice, inboard and outboard, so R_0 is not a function. A slope of the
 * wrong sign -- what a q read with the wrong orientation gives. A midplane that
 * stops short of a label, and one with a hole in it.
 */
BOOST_AUTO_TEST_CASE( aMidplaneThatIsNotMonotoneIsRefused )
{
	meq::MidplaneField const tokamak = []( double radius, double &x, double &dXdR )
	{
		x = 1.0 - ( radius - 1.0 )*( radius - 1.0 )/0.16;
		dXdR = -2.0*( radius - 1.0 )/0.16;
		return true;
	};
	std::string const folded = refusal( [ & ]
		{ meq::tabulateMidplaneCrossing( tokamak, 0.3, 1.8, 0.0, 1.0 ); } );
	BOOST_TEST( folded.find( "NOT MONOTONE" ) != std::string::npos, folded );

	meq::MidplaneField const mirror = mirrorNormalised();
	meq::MidplaneField const wrongSign = [ & ]( double radius, double &x, double &dXdR )
	{
		mirror( radius, x, dXdR );
		dXdR = -dXdR;
		return true;
	};
	std::string const sense = refusal( [ & ]
		{ meq::tabulateMidplaneCrossing( wrongSign, 0.0, 3.0, 0.0, 1.0 ); } );
	BOOST_TEST( sense.find( "against the sense" ) != std::string::npos, sense );

	std::string const shortOf = refusal( [ & ]
		{ meq::tabulateMidplaneCrossing( mirror, 0.0, 0.3, 0.0, 1.0 ); } );
	BOOST_TEST( shortOf.find( "not the whole of" ) != std::string::npos, shortOf );

	meq::MidplaneField const holed = [ & ]( double radius, double &x, double &dXdR )
	{
		mirror( radius, x, dXdR );
		return !( radius > 0.2 && radius < 0.21 );
	};
	std::string const hole = refusal( [ & ]
		{ meq::tabulateMidplaneCrossing( holed, 0.0, 3.0, 0.0, 1.0 ); } );
	BOOST_TEST( hole.find( "leaves the domain" ) != std::string::npos, hole );

	// And the one that is right, so the four above are refusals of their own
	// defects and not of the mirror.
	BOOST_TEST( refusal( [ & ] { meq::tabulateMidplaneCrossing( mirror, 0.0, 3.0, 0.0, 1.0 ); } ).empty() );
}

/*
 * WHAT THE GAUGE REFUSES, AND THAT THE FIXED RADIUS IS UNTOUCHED.
 */
BOOST_AUTO_TEST_CASE( theMidplaneGaugeRefusesWhatItCannotMean )
{
	// A radius beside the midplane gauge would be accepted and not read.
	BOOST_CHECK_THROW( meq::RotatingSource( hydrogenicSpecies(), rotationProfile(), ggPrimeProfile(), 1.0,
		meq::vacuumPermeability, meq::RotatingSource::Closure::Automatic,
		meq::RotatingSource::Gauge::Midplane ), std::invalid_argument );

	// No crossing yet: it says so rather than evaluating at some radius.
	meq::RotatingSource bare( hydrogenicSpecies(), rotationProfile(), ggPrimeProfile(), 0.0,
		meq::vacuumPermeability, meq::RotatingSource::Closure::Automatic,
		meq::RotatingSource::Gauge::Midplane );
	BOOST_TEST( bare.midplaneCrossing() == nullptr );
	BOOST_CHECK_THROW( bare.f( 1.0, 0.0, 0.5 ), std::logic_error );
	BOOST_CHECK_THROW( bare.setMidplaneCrossing( nullptr ), std::invalid_argument );

	// The fixed radius has no crossing to be given, and is the default.
	meq::RotatingSource fixed = makeSource();
	BOOST_TEST( ( fixed.gauge() == meq::RotatingSource::Gauge::FixedRadius ) );
	BOOST_CHECK_THROW( fixed.setMidplaneCrossing( crossingProfile() ), std::logic_error );

	meq::NormalisedRotatingSource normalised( hydrogenicSpecies(), rotationProfile(), ggPrimeProfile(),
		referenceRadius, 1.0 );
	BOOST_TEST( !normalised.readsMidplane() );
	BOOST_CHECK_THROW( normalised.setMidplaneField( mirrorNormalised(), 0.0, 3.0 ), std::logic_error );
}

/*
 * THE NORMALISED SOURCE TABULATES IN Psi, AT THE NORMALISATION IN FORCE, from a
 * field in PHYSICAL psi -- and its normalisation derivatives stay analytic for
 * the frozen table, which is the property the table being in Psi buys.
 */
BOOST_AUTO_TEST_CASE( theNormalisedSourceTabulatesInPsiAndKeepsItsColumnsAnalytic )
{
	double const psiInner = mirrorPsi( mirrorInner );
	double const psiOuter = mirrorPsi( mirrorOuter );

	meq::NormalisedRotatingSource source( hydrogenicSpecies(), rotationProfile(), ggPrimeProfile(),
		0.0, psiInner, meq::vacuumPermeability, meq::RotatingSource::Closure::Automatic,
		meq::RotatingSource::Gauge::Midplane );
	BOOST_TEST_REQUIRE( source.readsMidplane() );
	source.setNormalisation( psiInner, psiOuter );

	meq::MidplaneField const physical = []( double radius, double &psi, double &dPsiDR )
	{
		psi = mirrorPsi( radius );
		dPsiDR = mirrorPsiPrime( radius );
		return true;
	};
	source.setMidplaneField( physical, 0.0, 3.0 );
	meq::Profile const *table = source.unnormalised().midplaneCrossing();
	BOOST_TEST_REQUIRE( table != nullptr );
	for ( double x : { 0.0, 0.25, 0.5, 0.75, 1.0 } )
		BOOST_TEST( ( *table )( x ) == mirrorCrossing( x ), boost::test_tools::tolerance( 1.0e-11 ) );

	// phi_0 vanishes at each surface's crossing, read in physical psi.
	for ( double x : { 0.1, 0.5, 0.9 } )
	{
		double const psi = psiOuter + x*( psiInner - psiOuter );
		checkClose( source.potential( mirrorCrossing( x ), psi )/keV, 0.0, 1.0e-10,
			"normalised midplane source, e phi_0 at the crossing", x );
	}

	// The border columns, against a Richardson difference over the
	// normalisation with the table held -- the freeze.
	double const z = 0.0;
	double worst = 0.0;
	for ( double radius : { 0.15, 0.22, 0.3 } )
		for ( double x : { 0.2, 0.5, 0.8 } )
		{
			double const psi = psiOuter + x*( psiInner - psiOuter );
			source.setNormalisation( psiInner, psiOuter );
			double dAxis = 0.0, dBoundary = 0.0;
			BOOST_TEST_REQUIRE( source.normalisationDerivatives( radius, z, psi, dAxis, dBoundary ) );
			auto at = [ & ]( double step, bool axis )
			{
				double const scale = psiInner - psiOuter;
				source.setNormalisation( psiInner + ( axis ? step*scale : 0.0 ), psiOuter + ( axis ? 0.0 : step*scale ) );
				double const plus = source.f( radius, z, psi );
				source.setNormalisation( psiInner - ( axis ? step*scale : 0.0 ), psiOuter - ( axis ? 0.0 : step*scale ) );
				double const minus = source.f( radius, z, psi );
				return ( plus - minus )/( 2.0*step*scale );
			};
			double const h = 1.0e-4;
			double const numericAxis = ( 4.0*at( h/2.0, true ) - at( h, true ) )/3.0;
			double const numericBoundary = ( 4.0*at( h/2.0, false ) - at( h, false ) )/3.0;
			worst = std::max( { worst, std::fabs( dAxis - numericAxis )/std::fabs( numericAxis ),
			                    std::fabs( dBoundary - numericBoundary )/std::fabs( numericBoundary ) } );
		}
	source.setNormalisation( psiInner, psiOuter );
	std::printf( "\n  MIDPLANE GAUGE, THE BORDER COLUMNS FOR THE FROZEN TABLE: worst relative %.3e\n", worst );
	BOOST_TEST( worst < 1.0e-7 );
}

BOOST_AUTO_TEST_SUITE_END()

/*
 * THE FIXED-CURVE GAUGE, meq::RotatingSource::Gauge::FixedCurve: phi_0 pinned on
 * a curve R_g( Psi ) given once and never replaced -- the vacuum's crossing, in
 * the solver. MIRROR-PLAN.md section 5, [M-VACGAUGE](MEASUREMENTS.md#m-vacgauge).
 */
namespace
{
	/// Densities that vanish with their slope at both labels, so F does too
	/// whatever omega is: n = 16 n_peak Psi^2 ( 1 - Psi )^2.
	std::shared_ptr<meq::Profile const> vanishingDensity( double peak )
	{
		return analytic(
			[ peak ]( double x ) { return 16.0*peak*x*x*( 1.0 - x )*( 1.0 - x ); },
			[ peak ]( double x ) { return 32.0*peak*x*( 1.0 - x )*( 1.0 - 2.0*x ); },
			[ peak ]( double x ) { return 32.0*peak*( ( 1.0 - 2.0*x )*( 1.0 - 2.0*x ) - 2.0*x*( 1.0 - x ) ); } );
	}

	/// Densities FLAT at both labels and not vanishing there: n' = 0 at
	/// Psi = 0 and 1, n = n_peak there and 1.5 n_peak between.
	std::shared_ptr<meq::Profile const> flatDensity( double peak )
	{
		return analytic(
			[ peak ]( double x ) { return peak*( 1.0 + 0.5*x*x*( 3.0 - 2.0*x ) ); },
			[ peak ]( double x ) { return peak*3.0*x*( 1.0 - x ); },
			[ peak ]( double x ) { return peak*3.0*( 1.0 - 2.0*x ); } );
	}

	/// Two species at constant temperatures on @a ionDensity.
	std::vector<meq::Species> labelledPair( std::shared_ptr<meq::Profile const> ionDensity )
	{
		std::vector<meq::Species> species( 2 );
		species[ 0 ].mass = deuteronMass;
		species[ 0 ].charge = 1.0;
		species[ 0 ].temperature = std::make_shared<meq::ConstantProfile const>( 1.8*keV );
		species[ 0 ].density = std::move( ionDensity );
		species[ 1 ].mass = electronMass;
		species[ 1 ].charge = -1.0;
		species[ 1 ].temperature = std::make_shared<meq::ConstantProfile const>( 1.3*keV );
		species[ 1 ].density = meq::neutralisingDensity( species, 1 );
		return species;
	}

	/// The closed-form mirror's midplane in PHYSICAL psi.
	meq::MidplaneField mirrorPhysical()
	{
		return []( double radius, double &psi, double &dPsiDR )
		{
			psi = mirrorPsi( radius );
			dPsiDR = mirrorPsiPrime( radius );
			return true;
		};
	}

	/// A normalised fixed-curve source, g g' = 0 as on a mirror, at a STARTING
	/// normalisation that is not the vacuum's -- what the solve may start from.
	meq::NormalisedRotatingSource fixedCurveSource( std::vector<meq::Species> species,
	                                                std::shared_ptr<meq::Profile const> omega )
	{
		double const psiInner = mirrorPsi( mirrorInner );
		double const psiOuter = mirrorPsi( mirrorOuter );
		meq::NormalisedRotatingSource source( std::move( species ), std::move( omega ),
			std::make_shared<meq::ConstantProfile const>( 0.0 ), 0.0, 1.3*psiInner,
			meq::vacuumPermeability, meq::RotatingSource::Closure::Automatic,
			meq::RotatingSource::Gauge::FixedCurve );
		source.setNormalisation( 1.3*psiInner, 0.7*psiOuter );
		return source;
	}
}

BOOST_AUTO_TEST_SUITE( fixed_curve_gauge )

/*
 * ON THE SAME TABLE THE FIXED CURVE IS THE MIDPLANE GAUGE'S ALGEBRA, VALUE FOR
 * VALUE -- so theMidplaneJacobianAgreesWithACentralDifference covers its
 * derivatives -- and what differs is who may replace the table: nobody.
 */
BOOST_AUTO_TEST_CASE( theFixedCurveIsTheMidplaneAlgebraOnATableNobodyReplaces )
{
	for ( bool three : { false, true } )
	{
		std::vector<meq::Species> const species = three ? impuritySpecies() : hydrogenicSpecies();
		meq::RotatingSource const midplane = midplaneSource( species, 1.0, crossingProfile() );
		meq::RotatingSource curve( species, rotationProfile(), ggPrimeProfile(), 0.0,
			meq::vacuumPermeability, meq::RotatingSource::Closure::Automatic,
			meq::RotatingSource::Gauge::FixedCurve );
		BOOST_CHECK_THROW( curve.f( 1.0, 0.0, 0.5 ), std::logic_error );
		curve.setMidplaneCrossing( crossingProfile() );
		for ( double radius : testRadii )
			for ( double psi : testFluxes )
			{
				BOOST_TEST( curve.f( radius, 0.0, psi ) == midplane.f( radius, 0.0, psi ) );
				BOOST_TEST( curve.dFdPsi( radius, 0.0, psi ) == midplane.dFdPsi( radius, 0.0, psi ) );
			}
		for ( double psi : testFluxes )
			checkClose( curve.potential( crossing( psi ), psi )/keV, 0.0, 1.0e-12,
				"e phi_0 on the fixed curve", psi );
	}

	// A radius beside a curve would be accepted and not read.
	BOOST_CHECK_THROW( meq::RotatingSource( hydrogenicSpecies(), rotationProfile(), ggPrimeProfile(), 1.0,
		meq::vacuumPermeability, meq::RotatingSource::Closure::Automatic,
		meq::RotatingSource::Gauge::FixedCurve ), std::invalid_argument );
	// The fixed radius has no curve to check on.
	BOOST_CHECK_THROW( makeSource().labelCheck( crossingProfile() ), std::logic_error );
}

/*
 * THE VACUUM'S CURVE IS TABULATED AGAINST THE VACUUM'S LABELS, NOT THE
 * NORMALISATION THE SOLVE STARTS FROM, AND IS NEVER REPLACED. And with the curve
 * fixed, the border columns -- where moving psi_ax or psi_bnd moves Psi and so
 * R_g( Psi ) -- are EXACT, against Richardson differences over the labels at
 * two normalisations, with nothing frozen; so is dF/dpsi in physical psi.
 */
BOOST_AUTO_TEST_CASE( theVacuumCurveIsTheVacuumsAndItsColumnsAreExact )
{
	double const psiInner = mirrorPsi( mirrorInner );
	double const psiOuter = mirrorPsi( mirrorOuter );
	meq::NormalisedRotatingSource source = fixedCurveSource( labelledPair( vanishingDensity( 1.0e19 ) ),
	                                                         rotationProfile( 2.0 ) );
	BOOST_TEST( source.readsVacuumMidplane() );
	BOOST_TEST( !source.readsMidplane() );
	BOOST_CHECK_THROW( source.setMidplaneField( mirrorPhysical(), 0.0, 3.0 ), std::logic_error );

	source.setVacuumMidplaneField( mirrorPhysical(), 0.0, 3.0, psiInner, psiOuter );
	BOOST_TEST( !source.readsVacuumMidplane() );
	meq::Profile const *table = source.unnormalised().midplaneCrossing();
	BOOST_TEST_REQUIRE( table != nullptr );
	for ( double x : { 0.0, 0.25, 0.5, 0.75, 1.0 } )
		BOOST_TEST( ( *table )( x ) == mirrorCrossing( x ), boost::test_tools::tolerance( 1.0e-11 ) );
	BOOST_TEST( ( *table )( 1.0 ) == mirrorInner, boost::test_tools::tolerance( 1.0e-14 ) );
	BOOST_TEST( ( *table )( 0.0 ) == mirrorOuter, boost::test_tools::tolerance( 1.0e-14 ) );
	BOOST_CHECK_THROW( source.setVacuumMidplaneField( mirrorPhysical(), 0.0, 3.0, psiInner, psiOuter ),
	                   std::logic_error );
	BOOST_TEST( table == source.unnormalised().midplaneCrossing() );

	double const z = 0.0;
	double worstColumn = 0.0, worstField = 0.0;
	for ( double scale : { 1.0, 1.1 } )
	{
		double const axis = scale*psiInner;
		double const boundary = psiOuter/scale;
		double const span = axis - boundary;
		for ( double radius : { 0.15, 0.22, 0.3 } )
			for ( double x : { 0.2, 0.5, 0.8 } )
			{
				double const psi = boundary + x*span;
				source.setNormalisation( axis, boundary );
				double dAxis = 0.0, dBoundary = 0.0;
				BOOST_TEST_REQUIRE( source.normalisationDerivatives( radius, z, psi, dAxis, dBoundary ) );
				double const dField = source.dFdPsi( radius, z, psi );
				auto labels = [ & ]( double step, bool onAxis )
				{
					source.setNormalisation( axis + ( onAxis ? step*span : 0.0 ), boundary + ( onAxis ? 0.0 : step*span ) );
					double const plus = source.f( radius, z, psi );
					source.setNormalisation( axis - ( onAxis ? step*span : 0.0 ), boundary - ( onAxis ? 0.0 : step*span ) );
					double const minus = source.f( radius, z, psi );
					source.setNormalisation( axis, boundary );
					return ( plus - minus )/( 2.0*step*span );
				};
				auto field = [ & ]( double step )
				{
					return ( source.f( radius, z, psi + step*span ) - source.f( radius, z, psi - step*span ) )/( 2.0*step*span );
				};
				double const h = 1.0e-4;
				double const numericAxis = ( 4.0*labels( h/2.0, true ) - labels( h, true ) )/3.0;
				double const numericBoundary = ( 4.0*labels( h/2.0, false ) - labels( h, false ) )/3.0;
				double const numericField = ( 4.0*field( h/2.0 ) - field( h ) )/3.0;
				worstColumn = std::max( { worstColumn, std::fabs( dAxis - numericAxis )/std::fabs( numericAxis ),
				                          std::fabs( dBoundary - numericBoundary )/std::fabs( numericBoundary ) } );
				worstField = std::max( worstField, std::fabs( dField - numericField )/std::fabs( numericField ) );
			}
	}
	std::printf( "\n  FIXED CURVE, THE BORDER COLUMNS AND dF/dpsi WITH NOTHING FROZEN: worst relative %.3e and %.3e\n",
	             worstColumn, worstField );
	BOOST_TEST( worstColumn < 1.0e-7 );
	BOOST_TEST( worstField < 1.0e-7 );
}

/*
 * THE BOX'S PRECONDITION, CHECKED BEFORE THE SOLVE. On the vacuum's curve a
 * density FLAT at a label still carries -n m omega^2 R_g R_g'/T in dn/dpsi where
 * omega does not vanish, so the flat profile set is refused at rotation and
 * accepted without it -- which is exactly the claim -- and refused WITHOUT the
 * curve going in. Densities vanishing with their slope pass at any rotation.
 */
BOOST_AUTO_TEST_CASE( aDensityFlatAtALabelIsRefusedOnTheVacuumCurve )
{
	double const psiInner = mirrorPsi( mirrorInner );
	double const psiOuter = mirrorPsi( mirrorOuter );
	auto rotationConstant = std::make_shared<meq::ConstantProfile const>( sonicOmega );

	meq::NormalisedRotatingSource flat = fixedCurveSource( labelledPair( flatDensity( 1.0e19 ) ), rotationConstant );
	std::string const refused = refusal( [ & ]
		{ flat.setVacuumMidplaneField( mirrorPhysical(), 0.0, 3.0, psiInner, psiOuter ); } );
	std::printf( "\n  A DENSITY FLAT AT A LABEL, ROTATING: %s\n", refused.c_str() );
	BOOST_TEST( refused.find( "label" ) != std::string::npos, refused );
	BOOST_TEST( flat.readsVacuumMidplane(), "a refused curve must not be installed" );

	meq::NormalisedRotatingSource still = fixedCurveSource( labelledPair( flatDensity( 1.0e19 ) ), nullptr );
	BOOST_TEST( refusal( [ & ] { still.setVacuumMidplaneField( mirrorPhysical(), 0.0, 3.0, psiInner, psiOuter ); } ).empty() );

	meq::NormalisedRotatingSource vanishing = fixedCurveSource( labelledPair( vanishingDensity( 1.0e19 ) ),
	                                                            rotationProfile( 2.0 ) );
	BOOST_TEST( refusal( [ & ] { vanishing.setVacuumMidplaneField( mirrorPhysical(), 0.0, 3.0, psiInner, psiOuter ); } ).empty() );

	// What the check measured, on the curve the solver would have installed.
	meq::NormalisedRotatingSource probe = fixedCurveSource( labelledPair( flatDensity( 1.0e19 ) ), rotationConstant );
	meq::RotatingSource::LabelCheck const c = probe.unnormalised().labelCheck(
		meq::tabulateMidplaneCrossing( mirrorNormalised(), 0.0, 3.0, 0.0, 1.0 ) );
	std::printf( "    flat and rotating: F at the label %.3e of the peak, at Psi = %g, R = %.4f\n",
	             c.relative, c.label, c.radius );
	BOOST_TEST( c.relative > 1.0e3*meq::labelTolerance );
	BOOST_TEST( c.radius <= ( c.label == 1.0 ? mirrorInner : mirrorOuter )*( 1.0 + 1.0e-12 ),
		"a label must be swept only over the radii its surface reaches" );
}

/*
 * dn_s/dR AND dn_s/dpsi AGAINST CENTRAL DIFFERENCES, both closures and the two
 * gauges with a radius dependence of their own: the partials densityOnCurve()
 * builds a table's slopes from.
 */
BOOST_AUTO_TEST_CASE( theDensityDerivativesAgreeWithCentralDifferences )
{
	std::vector<meq::RotatingSource> sources;
	sources.push_back( makeSource() );
	sources.push_back( makeImpuritySource() );
	sources.push_back( midplaneSource( hydrogenicSpecies(), 1.0, crossingProfile() ) );
	sources.push_back( midplaneSource( impuritySpecies(), 1.0, crossingProfile() ) );
	double worst = 0.0;
	for ( meq::RotatingSource const &source : sources )
		for ( std::size_t s = 0; s < source.species().size(); ++s )
			for ( double radius : testRadii )
				for ( double psi : { 0.15, 0.4, 0.62, 0.85 } )
				{
					double value = 0.0, dRadius = 0.0, dPsi = 0.0;
					source.densityDerivatives( s, radius, psi, value, dRadius, dPsi );
					checkRelative( value, source.density( s, radius, psi ), 1.0e-14, "n_s", radius, psi );
					double const h = 1.0e-5;
					double const byR = ( source.density( s, radius + h, psi ) - source.density( s, radius - h, psi ) )/( 2.0*h );
					double const byPsi = ( source.density( s, radius, psi + h ) - source.density( s, radius, psi - h ) )/( 2.0*h );
					checkRelative( dRadius, byR, 1.0e-7, "dn_s/dR against a central difference", radius, psi, value );
					checkRelative( dPsi, byPsi, 1.0e-7, "dn_s/dpsi against a central difference", radius, psi, value );
					worst = std::max( { worst, std::fabs( dRadius - byR )/std::max( std::fabs( byR ), value ),
					                    std::fabs( dPsi - byPsi )/std::max( std::fabs( byPsi ), value ) } );
				}
	std::printf( "\n  dn_s/dR AND dn_s/dpsi AGAINST CENTRAL DIFFERENCES: worst relative %.3e\n", worst );
}

/*
 * A DENSITY RELABELLED ONTO ITS OWN GAUGE'S CURVE IS ITS TABLE: the exponent is
 * zero along the curve, so the value at each knot is n_s0 and the slope -- the
 * chain rule dn/dpsi + ( dn/dR ) R_c' -- is n_s0', the two partials cancelling
 * against each other to what their own accuracy allows. And the fixed-radius
 * plasma relabelled onto the same curve is the closed-form transfer.
 */
BOOST_AUTO_TEST_CASE( aDensityRelabelledOntoItsOwnCurveIsItsTable )
{
	std::vector<meq::Knot> knots;
	for ( int i = 0; i <= 64; ++i )
	{
		double const x = i/64.0;
		knots.push_back( meq::Knot{ x, crossing( x ), crossingPrime( x ) } );
	}
	meq::SplineProfile const curve( knots );

	double worstValue = 0.0, worstSlope = 0.0;
	for ( bool three : { false, true } )
	{
		std::vector<meq::Species> const species = three ? impuritySpecies() : hydrogenicSpecies();
		meq::RotatingSource const source = midplaneSource( species, 1.0, crossingProfile() );
		for ( std::size_t s = 0; s < species.size(); ++s )
		{
			std::shared_ptr<meq::SplineProfile const> const table = meq::densityOnCurve( source, s, curve );
			BOOST_TEST_REQUIRE( table->knots().size() == knots.size() );
			for ( meq::Knot const &k : table->knots() )
			{
				double const n = ( *species[ s ].density )( k.psi );
				double const nPrime = species[ s ].density->prime( k.psi );
				worstValue = std::max( worstValue, std::fabs( k.value - n )/n );
				worstSlope = std::max( worstSlope, std::fabs( k.derivative - nPrime )/std::max( std::fabs( nPrime ), n ) );
			}
		}
	}

	// The fixed-radius plasma at constant T and omega, onto the curve: the
	// transfer exp( C ( R_0^2 - R_ref^2 )/2 ) of theTwoGaugesDescribeOnePlasma.
	double const tIon = 1.8*keV, tElectron = 1.3*keV, omega = 1.3*sonicOmega;
	double const c = omega*omega*( electronMass + deuteronMass )/( tElectron + tIon );
	std::vector<meq::Species> pair = labelledPair( ionDensityProfile() );
	meq::RotatingSource const fixed( pair, std::make_shared<meq::ConstantProfile const>( omega ),
	                                 ggPrimeProfile(), referenceRadius );
	std::shared_ptr<meq::SplineProfile const> const transferred = meq::densityOnCurve( fixed, 0, curve );
	double worstTransfer = 0.0;
	for ( meq::Knot const &k : transferred->knots() )
	{
		double const factor = std::exp( c*( crossingSquared( k.psi ) - referenceRadius*referenceRadius )/2.0 );
		double const g = c*crossingSquaredPrime( k.psi )/2.0;
		worstTransfer = std::max( { worstTransfer,
			std::fabs( k.value - ionDensity( k.psi )*factor )/( ionDensity( k.psi )*factor ),
			std::fabs( k.derivative - ( ionDensityPrime( k.psi ) + ionDensity( k.psi )*g )*factor )
				/( ionDensity( k.psi )*factor ) } );
	}
	std::printf( "\n  RELABELLED ONTO ITS OWN CURVE: value %.3e, slope %.3e relative; the fixed radius onto "
	             "it against the closed-form transfer %.3e\n", worstValue, worstSlope, worstTransfer );
	BOOST_TEST( worstValue < 1.0e-13 );
	BOOST_TEST( worstSlope < 1.0e-12 );
	BOOST_TEST( worstTransfer < 1.0e-13 );
}

BOOST_AUTO_TEST_SUITE_END()

/*
 * THE ROOT FIND ENDS WHEN ITS BRACKET DOES. Far from a plasma's crossings at
 * high Mach the potential that balances (97) is so large against T that one
 * ulp of it moves Sum_s Z_s n_s by more than the relative tolerance the
 * residual test asks for, so the residual can never pass and the bracket
 * collapses onto two adjacent doubles first. MEASURED on
 * examples/mirror-three-species.toml, whose output pass evaluates the
 * densities at every grid node and threw at the first node far enough out.
 * The root is then known to the last bit, which is the answer.
 */
BOOST_AUTO_TEST_CASE( theRootFindEndsWhenItsBracketDoes )
{
	meq::RotatingSource const source = makeImpuritySource( 12.0 );
	int balanced = 0;
	for ( double radius : { 0.6, 1.4, 2.0, 2.6 } )
		for ( double psi : testFluxes )
		{
			double phi = 0.0;
			BOOST_CHECK_NO_THROW( phi = source.potential( radius, psi ) );
			BOOST_TEST( std::isfinite( phi ), "R = " << radius << ", psi = " << psi );
			// Where the densities themselves still fit a double, they balance.
			// Out at R = 2 at this Mach they do not -- exp( A_s ) passes
			// exp( 709.8 ) -- and only the potential is asked for.
			double sum = 0.0, scale = 0.0;
			for ( std::size_t s = 0; s < 3; ++s )
			{
				double const n = source.density( s, radius, psi );
				sum += source.species()[ s ].charge*n;
				scale += std::fabs( source.species()[ s ].charge*n );
			}
			if ( std::isfinite( scale ) )
			{
				++balanced;
				BOOST_TEST( std::fabs( sum ) <= 1.0e-10*scale,
					"quasineutrality at R = " << radius << ", psi = " << psi << ": " << sum/scale );
			}
		}
	BOOST_TEST( balanced >= 6 );
}
