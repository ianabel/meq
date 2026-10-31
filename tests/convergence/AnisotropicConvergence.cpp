#define BOOST_TEST_MODULE AnisotropicConvergence
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <cmath>

#include "analytic/ManufacturedAnisotropic.hpp"
#include "analytic/ManufacturedNonlinear.hpp"

/*
 * ANISOTROPIC PRESSURE: the acceptance of ANISOTROPY-PLAN.md's AN-0 to AN-2.
 *
 * AN-0 is everything in this file so far: THE BENCHMARK BEFORE THE SOLVER.
 * tests/analytic/ManufacturedAnisotropic.hpp is checked here against its own
 * operator, against the symbolic closed form tools/anisotropic_equation.py
 * derives, and against Example 5, which it must become at Delta = 0. Nothing
 * here solves anything; AN-1 adds the solves of the position model and AN-2
 * those of the flux model, on NewtonConvergence.cpp's box.
 *
 * The equation and its reference are refs/AlphaGK.pdf eq (6.18), in MEQ's
 * convention
 *
 *     -div_bar( ( 1 - Delta )/R grad_bar psi ) = F/R.
 */

namespace
{
	/// NewtonConvergence.cpp's box, on which the fixture is posed.
	double const minRadius = 0.6;
	double const maxRadius = 1.4;
	double const zMin = -0.6;
	double const zMax = 0.6;

	using Fixture = meq::analytic::ManufacturedAnisotropic;

	struct Reference
	{
		double radius;
		double z;
		double f;
	};

	/// tools/anisotropic_equation.py's output, the symbolic F at the fixture's
	/// default parameters to 17 digits. Re-run it, and re-paste, if those
	/// parameters change.
	Reference const positionReference[] = {
		{ 0.7, -0.5, 11.777095845763922 },
		{ 1.0, 0.0, 17.074809983704749 },
		{ 1.25, 0.3, 2.6717596599129035 },
		{ 1.3, -0.2, -0.58176170764865970 },
	};

	Reference const fluxReference[] = {
		{ 0.7, -0.5, 9.9901315861702581 },
		{ 1.0, 0.0, 18.652143338077114 },
		{ 1.25, 0.3, -0.32517430295125946 },
		{ 1.3, -0.2, -1.9269188437797595 },
	};
}

/// The benchmark before the solver, part one: -R div( ( 1 - Delta )/R grad psi )
/// by central differences of psi() and delta() must equal f(), under both
/// models. A closed form that is consistent with itself and wrong -- a sign on
/// the grad Delta term, a missing R -- converges at the full rate to the wrong
/// function, and only this sees it.
BOOST_AUTO_TEST_CASE( theSourceMatchesTheAnisotropicOperator )
{
	for ( Fixture const &eq : { Fixture::position(), Fixture::flux() } )
	{
		for ( double radius = minRadius; radius <= maxRadius + 1.0e-12; radius += 0.1 )
		{
			for ( double z = zMin; z <= zMax + 1.0e-12; z += 0.15 )
			{
				double const operatorFD = eq.anisotropicOperatorFD( radius, z );
				double const minusF = -eq.f( radius, z );
				BOOST_TEST( std::abs( operatorFD - minusF ) < 1.0e-5*( 1.0 + std::abs( minusF ) ),
				            "model " << static_cast<int>( eq.model() ) << " at ( " << radius
				            << ", " << z << " ): R div( (1 - Delta)/R grad psi ) = " << operatorFD
				            << " but -F = " << minusF );
			}
		}
	}
}

/// Part two: f() is the closed form sympy derives, to round-off. The finite
/// difference above is good to about 1e-7 and cannot tell a transcription slip
/// in the last digits from the truth; this can.
BOOST_AUTO_TEST_CASE( theSourceIsTheSymbolicClosedForm )
{
	auto compare = []( Fixture const &eq, Reference const ( &table )[ 4 ] )
	{
		for ( Reference const &row : table )
		{
			double const value = eq.f( row.radius, row.z );
			BOOST_TEST( std::abs( value - row.f ) < 1.0e-13*std::max( 1.0, std::abs( row.f ) ),
			            "model " << static_cast<int>( eq.model() ) << " at ( " << row.radius
			            << ", " << row.z << " ): f = " << value << " but sympy gives " << row.f );
		}
	};
	compare( Fixture::position(), positionReference );
	compare( Fixture::flux(), fluxReference );
}

/// Part three: at Delta = 0 the fixture is Example 5 -- its F at its own
/// solution, and its flux the isotropic q = grad psi/R. Both models reach it,
/// the flux model through a zero amplitude.
BOOST_AUTO_TEST_CASE( zeroAnisotropyIsExample5 )
{
	meq::analytic::ManufacturedNonlinear const example5
		= meq::analytic::ManufacturedNonlinear::example5();

	for ( Fixture const &eq : { Fixture::position( 0.0 ), Fixture::flux( 0.0 ) } )
	{
		for ( double radius = minRadius; radius <= maxRadius + 1.0e-12; radius += 0.2 )
		{
			for ( double z = zMin; z <= zMax + 1.0e-12; z += 0.3 )
			{
				BOOST_TEST( eq.psi( radius, z ) == example5.psi( radius, z ) );

				double const isotropic = example5.f( radius, z, example5.psi( radius, z ) );
				BOOST_TEST( std::abs( eq.f( radius, z ) - isotropic )
				            < 1.0e-13*std::max( 1.0, std::abs( isotropic ) ),
				            "at ( " << radius << ", " << z << " ): F = " << eq.f( radius, z )
				            << " but Example 5's is " << isotropic );

				double qR = 0.0;
				double qZ = 0.0;
				eq.flux( radius, z, qR, qZ );
				double qR5 = 0.0;
				double qZ5 = 0.0;
				example5.flux( radius, z, qR5, qZ5 );
				// To one rounding and not bitwise: ( 1 - 0 )/R times grad psi
				// and grad psi divided by R are different operations.
				BOOST_TEST( std::abs( qR - qR5 ) < 1.0e-15*std::max( 1.0, std::abs( qR5 ) ) );
				BOOST_TEST( std::abs( qZ - qZ5 ) < 1.0e-15*std::max( 1.0, std::abs( qZ5 ) ) );
			}
		}
	}
}

/// Part four: grad Delta is the gradient of delta(), and D'( psi ) the
/// derivative of D( psi ). The first enters F; the second is what AN-2's Newton
/// carries in the flux block, where a wrong derivative costs quadratic
/// convergence and not the answer.
BOOST_AUTO_TEST_CASE( theAnisotropyDerivativesMatchFiniteDifferences )
{
	double const h = 1.0e-6;

	for ( Fixture const &eq : { Fixture::position(), Fixture::flux() } )
	{
		for ( double radius = minRadius; radius <= maxRadius + 1.0e-12; radius += 0.2 )
		{
			for ( double z = zMin; z <= zMax + 1.0e-12; z += 0.3 )
			{
				double dDeltaDr = 0.0;
				double dDeltaDz = 0.0;
				eq.gradDelta( radius, z, dDeltaDr, dDeltaDz );
				double const fdR = ( eq.delta( radius + h, z ) - eq.delta( radius - h, z ) )/( 2.0*h );
				double const fdZ = ( eq.delta( radius, z + h ) - eq.delta( radius, z - h ) )/( 2.0*h );
				BOOST_TEST( std::abs( dDeltaDr - fdR ) < 1.0e-8, "d Delta/dR at ( " << radius
				            << ", " << z << " ): " << dDeltaDr << " against " << fdR );
				BOOST_TEST( std::abs( dDeltaDz - fdZ ) < 1.0e-8, "d Delta/dz at ( " << radius
				            << ", " << z << " ): " << dDeltaDz << " against " << fdZ );
			}
		}
	}

	Fixture const eq = Fixture::flux();
	for ( double psi : { -1.0, -0.3, 0.0, 0.4, 1.0 } )
	{
		double const fd = ( eq.deltaOfPsi( psi + h ) - eq.deltaOfPsi( psi - h ) )/( 2.0*h );
		BOOST_TEST( std::abs( eq.dDeltaDPsi( psi ) - fd ) < 1.0e-8,
		            "D'( " << psi << " ) = " << eq.dDeltaDPsi( psi ) << " against " << fd );
	}
}

/// Part five: the fixture is elliptic on its box, 1 - Delta > 0 everywhere, by
/// a margin. Below the firehose limit is the precondition of the operator, and
/// AN-1's refusal at 1 - Delta <= 0 must never fire on this fixture.
BOOST_AUTO_TEST_CASE( theFixtureIsEllipticOnItsBox )
{
	for ( Fixture const &eq : { Fixture::position(), Fixture::flux() } )
	{
		double smallest = 1.0e300;
		for ( double radius = minRadius; radius <= maxRadius + 1.0e-12; radius += 0.01 )
			for ( double z = zMin; z <= zMax + 1.0e-12; z += 0.01 )
				smallest = std::min( smallest, 1.0 - eq.delta( radius, z ) );

		BOOST_TEST( smallest > 0.5, "model " << static_cast<int>( eq.model() )
		            << ": min( 1 - Delta ) = " << smallest );
	}
}
