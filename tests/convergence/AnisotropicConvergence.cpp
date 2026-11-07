#define BOOST_TEST_MODULE AnisotropicConvergence
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "mfem.hpp"

#if defined( MFEM_USE_OPENMP ) && defined( MFEM_THREAD_SAFE )
	// For omp_get_max_threads() alone, so that the threaded flux-law case can
	// say whether the run it reports could have exhibited a race.
	#include <omp.h>
#endif

#include "RoundOffFloor.hpp"

#include "meq/Anisotropy.hpp"
#include "meq/Estimator.hpp"
#include "meq/GradShafranov.hpp"

#include "analytic/ManufacturedAnisotropic.hpp"
#include "analytic/ManufacturedNonlinear.hpp"
#include "convergence/ConvergenceHarness.hpp"

/*
 * ANISOTROPIC PRESSURE: the acceptance of ANISOTROPY-PLAN.md's AN-0 to AN-2.
 *
 * AN-0, THE BENCHMARK BEFORE THE SOLVER: tests/analytic/ManufacturedAnisotropic.hpp
 * checked against its own operator, against the symbolic closed form
 * tools/anisotropic_equation.py derives, and against Example 5, which it must
 * become at Delta = 0.
 *
 * AN-1, THE OPERATOR WITH Delta A FUNCTION OF POSITION, through
 * GradShafranovSolver::setAnisotropy(), on NewtonConvergence.cpp's box:
 *
 *  - Delta = 0 reproduces the isotropic solve BIT FOR BIT;
 *  - a constant Delta IS the isotropic problem with F/( 1 - Delta ) and
 *    tau/( 1 - Delta ), to round-off, with q_Delta = ( 1 - Delta ) q_0. Not
 *    psi_0/( 1 - Delta ): the HDG stabilisation does not scale with the
 *    operator, so the exact discrete statement carries tau with it;
 *  - the position model converges at k+1 in psi and q and k+2 in psi*;
 *  - the refusals.
 *
 * AN-2a, Delta A FUNCTION OF psi, the fixture's flux model Delta = d psi^2,
 * through the non-linear flux law and Newton:
 *
 *  - the assembled Jacobian against a central difference of the residual;
 *  - the flux model converges at k+1 in psi and q and k+2 in psi*, Newton
 *    quadratic;
 *  - a zero-amplitude flux model, which takes the flux-law path, agrees with
 *    the isotropic Newton solve to the solver's tolerance, and the cost of
 *    that path per Newton step against the isotropic one is printed;
 *  - the refusals.
 *
 * AN-2b, Delta A FUNCTION OF THE NORMALISED FLUX Psi, on a normalised source
 * and the bordered Newton, the axis label read at a prescribed point:
 *
 *  - psi and q at k+1, psi_ax tending to the exact 1, Newton quadratic -- on a
 *    source that CLAIMS analytic label derivatives of zero, so that a solver
 *    taking the source's column rather than differencing would miss the flux
 *    law's dependence on psi_ax (measured as a mutation: order 1.000);
 *  - psi* NOT asserted at k+2: the label is psi_h at a point, so Delta carries
 *    psi_h's pointwise error, k+1 with a constant that moves with where the
 *    point sits in its element, and psi* cannot beat its own operator;
 *  - the refusals.
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

namespace
{
	/// The fixture's Delta, through the solver's interface.
	class FixtureAnisotropy : public meq::PositionAnisotropy
	{
		public:
			explicit FixtureAnisotropy( Fixture const &eqIn )
				: eq( eqIn )
			{
			}

			double delta( double radius, double z ) const override
			{
				return eq.delta( radius, z );
			}

		private:
			Fixture eq;
	};

	class ConstantAnisotropy : public meq::PositionAnisotropy
	{
		public:
			explicit ConstantAnisotropy( double valueIn )
				: value( valueIn )
			{
			}

			double delta( double, double ) const override
			{
				return value;
			}

		private:
			double value;
	};

	/// One solve of the position model on n cells a side, and its errors.
	struct AnisotropicMeasurement
	{
		double h = 0.0;
		int traceDofs = 0;
		double errorPsi = 0.0;
		double errorFlux = 0.0;
		double errorPsiStar = 0.0;
		double normPsi = 0.0;
		double normFlux = 0.0;
	};

	AnisotropicMeasurement measurePosition( Fixture const &eq, int order, int n )
	{
		meq::tests::Rectangle const box = meq::tests::standardBox();
		mfem::Mesh mesh = meq::tests::makeMesh( box, n );

		mfem::FunctionCoefficient sourceCoeff( [ &eq ]( mfem::Vector const &x )
		{
			return eq.f( x( 0 ), x( 1 ) );
		} );
		mfem::FunctionCoefficient psiCoeff( [ &eq ]( mfem::Vector const &x )
		{
			return eq.psi( x( 0 ), x( 1 ) );
		} );
		// The ANISOTROPIC flux, ( 1 - Delta )/R grad psi: the q the solve
		// carries. See the fixture's flux().
		mfem::VectorFunctionCoefficient fluxCoeff( 2, [ &eq ]( mfem::Vector const &x,
		                                                       mfem::Vector &value )
		{
			double fluxR = 0.0, fluxZ = 0.0;
			eq.flux( x( 0 ), x( 1 ), fluxR, fluxZ );
			value( 0 ) = fluxR;
			value( 1 ) = fluxZ;
		} );

		FixtureAnisotropy const anisotropy( eq );
		meq::GradShafranovSolver solver( mesh, order );
		solver.setAnisotropy( anisotropy );
		solver.setSource( sourceCoeff );
		solver.setBoundaryData( psiCoeff );
		solver.solve();
		solver.postProcess();

		AnisotropicMeasurement point;
		point.h = box.width()/static_cast<double>( n );
		point.traceDofs = solver.numTraceDofs();
		point.errorPsi = solver.potentialError( psiCoeff );
		point.errorFlux = solver.fluxError( fluxCoeff );
		point.errorPsiStar = solver.postProcessedPotentialError( psiCoeff );
		point.normPsi = solver.exactPotentialNorm( psiCoeff );
		point.normFlux = solver.exactFluxNorm( fluxCoeff );
		return point;
	}

	/// Solve a ManufacturedNonlinear-sourced linear problem, F evaluated at
	/// Example 5's exact psi and scaled by @a sourceScale, under @a anisotropy
	/// (or none) at constant tau, and hand back the solver for inspection.
	struct Example5Solve
	{
		mfem::Mesh mesh;
		mfem::FunctionCoefficient sourceCoeff;
		mfem::FunctionCoefficient psiCoeff;
		meq::GradShafranovSolver solver;

		Example5Solve( int order, int n, double tau, double sourceScale,
		               meq::Anisotropy const *anisotropy )
			: mesh( meq::tests::makeMesh( meq::tests::standardBox(), n ) ),
			  sourceCoeff( [ sourceScale ]( mfem::Vector const &x )
			  {
				  meq::analytic::ManufacturedNonlinear const eq
					  = meq::analytic::ManufacturedNonlinear::example5();
				  return sourceScale*eq.f( x( 0 ), x( 1 ), eq.psi( x( 0 ), x( 1 ) ) );
			  } ),
			  psiCoeff( []( mfem::Vector const &x )
			  {
				  return meq::analytic::ManufacturedNonlinear::example5().psi( x( 0 ), x( 1 ) );
			  } ),
			  solver( mesh, order, tau )
		{
			if ( anisotropy != nullptr )
				solver.setAnisotropy( *anisotropy );
			solver.setSource( sourceCoeff );
			solver.setBoundaryData( psiCoeff );
			solver.solve();
		}
	};

	double maxDifference( mfem::GridFunction const &a, mfem::GridFunction const &b,
	                      double scale = 1.0 )
	{
		double worst = 0.0;
		for ( int i = 0; i < a.Size(); ++i )
			worst = std::max( worst, std::abs( a( i ) - scale*b( i ) ) );
		return worst;
	}
}

/// Delta = 0 through setAnisotropy() must be the isotropic solve to the bit: the
/// weight is R times 1/( 1 - 0 ) = R exactly, on the same rule. Anything less
/// means the anisotropic path assembles something other than the operator it
/// claims to at the point where the two coincide.
BOOST_AUTO_TEST_CASE( zeroAnisotropyIsTheIsotropicSolveBitForBit )
{
	ConstantAnisotropy const zero( 0.0 );
	for ( int order : { 1, 3 } )
	{
		Example5Solve const isotropic( order, 8, 1.0, 1.0, nullptr );
		Example5Solve const anisotropic( order, 8, 1.0, 1.0, &zero );

		BOOST_TEST( maxDifference( anisotropic.solver.potential(), isotropic.solver.potential() ) == 0.0,
		            "k = " << order << ": psi differs by "
		            << maxDifference( anisotropic.solver.potential(), isotropic.solver.potential() ) );
		BOOST_TEST( maxDifference( anisotropic.solver.flux(), isotropic.solver.flux() ) == 0.0,
		            "k = " << order << ": q differs by "
		            << maxDifference( anisotropic.solver.flux(), isotropic.solver.flux() ) );
	}
}

/// A constant Delta is the isotropic problem with F/( 1 - Delta ) and
/// tau/( 1 - Delta ): substituting q = ( 1 - Delta ) q' into the anisotropic
/// HDG system divides its conservation and constraint rows by 1 - Delta, which
/// scales the source and the stabilisation together. So the anisotropic solve
/// at tau = 1 - Delta must equal the isotropic one at tau = 1 and source
/// F/( 1 - Delta ) to round-off, with q_Delta = ( 1 - Delta ) q_0. Both signs
/// of Delta, since a sign error in how the weight is formed would show only on
/// one side of zero for some mistakes and not others.
BOOST_AUTO_TEST_CASE( aConstantAnisotropyIsTheIsotropicProblemRescaled )
{
	for ( double delta : { -0.4, 0.3 } )
	{
		ConstantAnisotropy const constant( delta );
		double const tension = 1.0 - delta;
		for ( int order : { 1, 3 } )
		{
			Example5Solve const isotropic( order, 8, 1.0, 1.0/tension, nullptr );
			Example5Solve const anisotropic( order, 8, tension, 1.0, &constant );

			double const psiScale = isotropic.solver.potential().Normlinf();
			double const fluxScale = isotropic.solver.flux().Normlinf();
			double const psiDifference
				= maxDifference( anisotropic.solver.potential(), isotropic.solver.potential() );
			double const fluxDifference
				= maxDifference( anisotropic.solver.flux(), isotropic.solver.flux(), tension );

			// TWO SOLVES OF ONE PROBLEM POSED TWO WAYS, so they meet at
			// round-off and not at zero. 1e-12 relative is 4500 of double's
			// eps; single reads 160 to 1250 of its own (Delta = 0.3, k = 3,
			// psi), so exactAgreementFloor()'s 4096 eps is the same statement
			// in that format, and exactly zero in double. A weight formed as
			// R( 1 - Delta ) is O( Delta ) out. MEASUREMENTS.md M-ANISO-1.
			double const psiBound = std::max( 1.0e-12*psiScale,
			                                  meq::tests::exactAgreementFloor( psiScale ) );
			double const fluxBound = std::max( 1.0e-12*fluxScale,
			                                   meq::tests::exactAgreementFloor( fluxScale ) );
			BOOST_TEST( psiDifference < psiBound,
			            "Delta = " << delta << ", k = " << order << ": psi differs by "
			            << psiDifference << " on a scale of " << psiScale );
			BOOST_TEST( fluxDifference < fluxBound,
			            "Delta = " << delta << ", k = " << order << ": q - ( 1 - Delta ) q_0 is "
			            << fluxDifference << " on a scale of " << fluxScale );
		}
	}
}

namespace
{
	std::vector<int> const meshSizes = { 4, 8, 16, 32 };
	double const rateSlack = 0.15;

	void printTable( char const *model, int order,
	                 std::vector<AnisotropicMeasurement> const &points,
	                 char const *box = "standard box" )
	{
		std::printf( "\n  manufactured anisotropic, %s, k = %d, %s\n", model, order, box );
		std::printf( "  %8s %9s %13s %6s %13s %6s %13s %6s\n",
		             "h", "trace", "L2(psi)", "rate", "L2(q)", "rate", "L2(psi*)", "rate" );
		for ( std::size_t i = 0; i < points.size(); ++i )
		{
			AnisotropicMeasurement const &p = points[ i ];
			if ( i == 0 )
			{
				std::printf( "  %8.5f %9d %13.6e %6s %13.6e %6s %13.6e %6s\n", p.h, p.traceDofs,
				             p.errorPsi, "-", p.errorFlux, "-", p.errorPsiStar, "-" );
				continue;
			}
			AnisotropicMeasurement const &c = points[ i - 1 ];
			double const ratio = c.h/p.h;
			std::printf( "  %8.5f %9d %13.6e %6.3f %13.6e %6.3f %13.6e %6.3f\n", p.h, p.traceDofs,
			             p.errorPsi, meq::tests::rate( c.errorPsi, p.errorPsi, ratio ),
			             p.errorFlux, meq::tests::rate( c.errorFlux, p.errorFlux, ratio ),
			             p.errorPsiStar, meq::tests::rate( c.errorPsiStar, p.errorPsiStar, ratio ) );
		}
		std::fflush( stdout );
	}

	/// The rates, and the finest errors against ceilings at about three times
	/// their measured values: the rate is blind to a solution wrong by a factor,
	/// and a weight of R( 1 - Delta ) in place of R/( 1 - Delta ) converges at
	/// the full rate to a different function.
	void checkPositionOrder( int order, double psiCeiling, double fluxCeiling )
	{
		Fixture const eq = Fixture::position();
		std::vector<AnisotropicMeasurement> points;
		for ( int n : meshSizes )
			points.push_back( measurePosition( eq, order, n ) );
		printTable( "position Delta", order, points );

		if ( !meq::tests::rateIsMeasurable( order ) )
			BOOST_TEST_MESSAGE( meq::tests::rateNotMeasurableNotice( order ) );
		else
		{
			double const expected = order + 1.0 - rateSlack;
			double const expectedStar = order + 2.0 - rateSlack;
			for ( std::size_t i = 1; i < points.size(); ++i )
			{
				AnisotropicMeasurement const &c = points[ i - 1 ];
				AnisotropicMeasurement const &p = points[ i ];
				double const ratio = c.h/p.h;
				double const ratePsi = meq::tests::rate( c.errorPsi, p.errorPsi, ratio );
				double const rateFlux = meq::tests::rate( c.errorFlux, p.errorFlux, ratio );
				double const rateStar = meq::tests::rate( c.errorPsiStar, p.errorPsiStar, ratio );
				BOOST_TEST( ratePsi >= expected, "k = " << order << ", h = " << p.h
				            << ": psi converged at " << ratePsi << ", wanted " << expected );
				BOOST_TEST( rateFlux >= expected, "k = " << order << ", h = " << p.h
				            << ": q converged at " << rateFlux << ", wanted " << expected );
				if ( p.errorPsiStar < meq::tests::roundOffFloor( p.normPsi ) )
					BOOST_TEST_MESSAGE( "  k = " << order << ", h = " << p.h
					                    << ": psi* under the round-off floor, rate not asserted" );
				else
					BOOST_TEST( rateStar >= expectedStar, "k = " << order << ", h = " << p.h
					            << ": psi* converged at " << rateStar << ", wanted " << expectedStar );
			}
		}

		meq::tests::checkCeilings( order, points.back().errorPsi, points.back().normPsi,
		                           psiCeiling, points.back().errorFlux, points.back().normFlux,
		                           fluxCeiling );
	}
}

BOOST_AUTO_TEST_CASE( thePositionModelConvergesAtK1 )
{
	// Measured at h = 0.025: psi 5.55e-04, q 1.02e-03.
	checkPositionOrder( 1, 1.7e-03, 3.1e-03 );
}

BOOST_AUTO_TEST_CASE( thePositionModelConvergesAtK2 )
{
	// Measured at h = 0.025: psi 5.62e-06, q 1.09e-05.
	checkPositionOrder( 2, 1.7e-05, 3.3e-05 );
}

BOOST_AUTO_TEST_CASE( thePositionModelConvergesAtK3 )
{
	// Measured at h = 0.025: psi 4.44e-08, q 9.45e-08.
	checkPositionOrder( 3, 1.4e-07, 2.9e-07 );
}

BOOST_AUTO_TEST_CASE( thePositionModelConvergesAtK4 )
{
	// Measured at h = 0.025: psi 3.79e-10, q 8.77e-10.
	checkPositionOrder( 4, 1.2e-09, 2.7e-09 );
}

/// 1 - Delta <= 0 anywhere the solve reaches is refused at the first solve, and
/// named, rather than clipped or assembled into an operator that is not
/// elliptic. A Delta reaching 1 only near one corner of the box is the case a
/// check at element centres alone would miss.
BOOST_AUTO_TEST_CASE( theFirehoseLimitIsRefused )
{
	class CornerFirehose : public meq::PositionAnisotropy
	{
		public:
			double delta( double radius, double z ) const override
			{
				return ( radius > 1.39 && z > 0.59 ) ? 1.0 : 0.0;
			}
	};
	CornerFirehose const firehose;

	mfem::Mesh mesh = meq::tests::makeMesh( meq::tests::standardBox(), 4 );
	mfem::ConstantCoefficient zero( 0.0 );
	meq::GradShafranovSolver solver( mesh, 2 );
	solver.setAnisotropy( firehose );
	solver.setSource( zero );
	solver.setBoundaryData( zero );
	BOOST_CHECK_THROW( solver.solve(), std::invalid_argument );
}

/// One solver holds one anisotropy, set before the forms are built.
BOOST_AUTO_TEST_CASE( theAnisotropyIsSetOnceAndBeforeTheFirstSolve )
{
	ConstantAnisotropy const a( -0.1 );
	ConstantAnisotropy const b( -0.2 );

	mfem::Mesh mesh = meq::tests::makeMesh( meq::tests::standardBox(), 4 );
	mfem::ConstantCoefficient zero( 0.0 );
	meq::GradShafranovSolver solver( mesh, 1 );
	solver.setAnisotropy( a );
	BOOST_CHECK_THROW( solver.setAnisotropy( b ), std::logic_error );
	BOOST_TEST( solver.anisotropy() == &a );

	meq::GradShafranovSolver built( mesh, 1 );
	built.setSource( zero );
	built.setBoundaryData( zero );
	built.solve();
	BOOST_CHECK_THROW( built.setAnisotropy( a ), std::logic_error );
}

/// The residual estimator's terms are the isotropic operator's, so it refuses a
/// solver carrying an anisotropy rather than reporting an indicator of a
/// different equation.
BOOST_AUTO_TEST_CASE( theEstimatorRefusesAnAnisotropicSolve )
{
	ConstantAnisotropy const constant( -0.2 );
	Example5Solve solve( 1, 4, 1.0, 1.0, &constant );
	meq::ResidualEstimator estimator( solve.solver, solve.sourceCoeff );
	BOOST_CHECK_THROW( estimator.GetLocalErrors(), std::logic_error );
}

/*
 * ---------------------------------------------------------------------------
 * AN-2a: Delta a function of psi
 * ---------------------------------------------------------------------------
 */
namespace
{
	/// The fixture's flux model, Delta = d psi^2, through the solver's interface.
	class FluxModelAnisotropy : public meq::Anisotropy
	{
		public:
			explicit FluxModelAnisotropy( Fixture const &eqIn )
				: eq( eqIn )
			{
			}

			double delta( double, double, double psi ) const override
			{
				return eq.deltaOfPsi( psi );
			}

			double dDeltaDPsi( double, double, double psi ) const override
			{
				return eq.dDeltaDPsi( psi );
			}

			Argument argument() const override
			{
				return Argument::Flux;
			}

		private:
			Fixture eq;
	};

	/// The fixture's F, which is a function of position, as a meq::Source: the
	/// flux law is non-linear, so the solve must take the Newton path.
	class FixtureSource : public meq::Source
	{
		public:
			explicit FixtureSource( Fixture const &eqIn )
				: eq( eqIn )
			{
			}

			double f( double radius, double z, double ) const override
			{
				return eq.f( radius, z );
			}

			double dFdPsi( double, double, double ) const override
			{
				return 0.0;
			}

		private:
			Fixture eq;
	};

	struct FluxMeasurement
	{
		AnisotropicMeasurement point;
		int newtonIterations = 0;
		std::vector<double> history;
	};

	FluxMeasurement measureFlux( Fixture const &eq, int order, int n )
	{
		meq::tests::Rectangle const box = meq::tests::standardBox();
		mfem::Mesh mesh = meq::tests::makeMesh( box, n );

		mfem::FunctionCoefficient psiCoeff( [ &eq ]( mfem::Vector const &x )
		{
			return eq.psi( x( 0 ), x( 1 ) );
		} );
		mfem::VectorFunctionCoefficient fluxCoeff( 2, [ &eq ]( mfem::Vector const &x,
		                                                       mfem::Vector &value )
		{
			double fluxR = 0.0, fluxZ = 0.0;
			eq.flux( x( 0 ), x( 1 ), fluxR, fluxZ );
			value( 0 ) = fluxR;
			value( 1 ) = fluxZ;
		} );

		FluxModelAnisotropy const anisotropy( eq );
		FixtureSource const source( eq );
		meq::GradShafranovSolver solver( mesh, order );
		solver.setAnisotropy( anisotropy );
		solver.setSource( source );
		solver.setBoundaryData( psiCoeff );
		solver.solve();
		solver.postProcess();

		FluxMeasurement result;
		AnisotropicMeasurement &point = result.point;
		point.h = box.width()/static_cast<double>( n );
		point.traceDofs = solver.numTraceDofs();
		point.errorPsi = solver.potentialError( psiCoeff );
		point.errorFlux = solver.fluxError( fluxCoeff );
		point.errorPsiStar = solver.postProcessedPotentialError( psiCoeff );
		point.normPsi = solver.exactPotentialNorm( psiCoeff );
		point.normFlux = solver.exactFluxNorm( fluxCoeff );
		result.newtonIterations = solver.newtonIterations();
		result.history = solver.newtonResiduals();
		return result;
	}

	void checkFluxOrder( int order, double psiCeiling, double fluxCeiling )
	{
		Fixture const eq = Fixture::flux();
		std::vector<AnisotropicMeasurement> points;
		std::vector<int> iterations;
		double best = 0.0;
		for ( int n : meshSizes )
		{
			FluxMeasurement const m = measureFlux( eq, order, n );
			if ( n == 16 )
			{
				std::printf( "\n  [flux model] k = %d, n = 16, Newton residuals:", order );
				for ( double r : m.history )
					std::printf( " %.3e", r );
				std::printf( "\n" );
			}
			points.push_back( m.point );
			iterations.push_back( m.newtonIterations );
			best = std::max( best, meq::tests::bestNewtonOrder( m.history ) );
		}
		std::printf( "\n  [flux model, Delta = -0.5 psi^2] Newton iterations:" );
		for ( int it : iterations )
			std::printf( " %d", it );
		std::printf( "; best observed order %.3f\n", best );
		printTable( "flux Delta = -0.5 psi^2", order, points );

		if ( !meq::tests::rateIsMeasurable( order ) )
			BOOST_TEST_MESSAGE( meq::tests::rateNotMeasurableNotice( order ) );
		else
		{
			double const expected = order + 1.0 - rateSlack;
			double const expectedStar = order + 2.0 - rateSlack;
			for ( std::size_t i = 1; i < points.size(); ++i )
			{
				AnisotropicMeasurement const &c = points[ i - 1 ];
				AnisotropicMeasurement const &p = points[ i ];
				double const ratio = c.h/p.h;
				double const ratePsi = meq::tests::rate( c.errorPsi, p.errorPsi, ratio );
				double const rateFlux = meq::tests::rate( c.errorFlux, p.errorFlux, ratio );
				double const rateStar = meq::tests::rate( c.errorPsiStar, p.errorPsiStar, ratio );
				BOOST_TEST( ratePsi >= expected, "k = " << order << ", h = " << p.h
				            << ": psi converged at " << ratePsi << ", wanted " << expected );
				BOOST_TEST( rateFlux >= expected, "k = " << order << ", h = " << p.h
				            << ": q converged at " << rateFlux << ", wanted " << expected );
				// psi* FROM THE SECOND PAIR: at k = 4 the coarsest pair reads 5.78
				// and the next two 5.98 and 6.00, pre-asymptotic on 280 trace
				// dofs rather than short of the rate. A precondition relaxed, not
				// the tolerance.
				if ( i == 1 )
					continue;
				if ( p.errorPsiStar < meq::tests::roundOffFloor( p.normPsi ) )
					BOOST_TEST_MESSAGE( "  k = " << order << ", h = " << p.h
					                    << ": psi* under the round-off floor, rate not asserted" );
				else
					BOOST_TEST( rateStar >= expectedStar, "k = " << order << ", h = " << p.h
					            << ": psi* converged at " << rateStar << ", wanted " << expectedStar );
			}
		}

		// The flux law is non-linear, so Newton must take more than one step --
		// one would mean the law's psi-dependence never reached the operator --
		// and a Jacobian consistent with its residual makes it quadratic.
		for ( int it : iterations )
			BOOST_TEST( it >= 2, "k = " << order << ": Newton converged in " << it
			            << " step, so the flux law behaved linearly" );
		//
		// IN SINGLE THE TRIPLE THAT READS 1.8 IS UNDER THE FLOOR. Double reads
		// it on ( 1.06e-02, 4.88e-05, 2.67e-09 ), and 2.67e-09 is 170 times
		// under single's floor at 4.6e-07; what single keeps is the triple
		// before it, which reads 1.47 to 1.52 in BOTH formats -- the
		// trajectory's, before the quadratic regime. On that triple, measured
		// in single, dDelta/dpsi 10% high reads 1.28 to 1.31, 10% low 1.22 to
		// 1.24, and dropped 1.16, against 1.49 to 1.52 exact: 1.4 sits between
		// with 0.09 either side, thin and recorded. MEASUREMENTS.md M-ANISO-1.
#ifdef MFEM_USE_SINGLE
		double const orderBound = 1.4;
#else
		double const orderBound = 1.8;
#endif
		BOOST_TEST( best >= orderBound, "k = " << order << ": the best observed Newton order is "
		            << best << " against " << orderBound
		            << "; the flux law's Jacobian disagrees with its residual" );

		meq::tests::checkCeilings( order, points.back().errorPsi, points.back().normPsi,
		                           psiCeiling, points.back().errorFlux, points.back().normFlux,
		                           fluxCeiling );
	}
}

/// THE FLUX LAW'S JACOBIAN AGAINST A CENTRAL DIFFERENCE OF THE RESIDUAL, as
/// NewtonConvergence.cpp's jacobianMatchesAFiniteDifferenceOfTheResidual does
/// it and for its reasons: on the CONDENSED operator, at a state away from the
/// solution, the linearisation held at that state, the essential trace dofs
/// left out of every direction. Under a flux law the hybridization's local
/// operator is FullNL, so this reaches the law's df through the element-local
/// elimination. Beside it, checkFluxOrder()'s observed Newton order is the
/// instrument on the NPC path the solve drives.
BOOST_AUTO_TEST_CASE( theFluxLawsJacobianMatchesAFiniteDifferenceOfTheResidual )
{
	Fixture const eq = Fixture::flux();
	FluxModelAnisotropy const anisotropy( eq );
	FixtureSource const source( eq );
	mfem::Mesh mesh = meq::tests::makeMesh( meq::tests::standardBox(), 3 );
	mfem::FunctionCoefficient psiCoeff( [ &eq ]( mfem::Vector const &x )
	{
		return eq.psi( x( 0 ), x( 1 ) );
	} );

	meq::GradShafranovSolver solver( mesh, 2 );
	solver.setAnisotropy( anisotropy );
	solver.setSource( source );
	solver.setBoundaryData( psiCoeff );
	solver.prepare();

	mfem::Operator &residual = solver.reducedOperator();
	int const size = solver.reducedSolution().Size();
	std::vector<bool> isEssential( static_cast<std::size_t>( size ), false );
	mfem::Array<int> const &essential = solver.essentialTraceDofs();
	for ( int i = 0; i < essential.Size(); ++i )
		isEssential[ static_cast<std::size_t>( essential[ i ] ) ] = true;
	BOOST_TEST_REQUIRE( essential.Size() > 0 );

	mfem::Vector state( solver.reducedSolution() );
	for ( int i = 0; i < size; ++i )
		if ( !isEssential[ static_cast<std::size_t>( i ) ] )
			state( i ) += 0.3*std::sin( 1.7*i + 0.4 );

	residual.GetGradient( state );

	double const step = meq::tests::jacobianDifferenceStep( 1.0e-5 );
	int const numDirections = 4;
	std::vector<mfem::Vector> directions( numDirections );
	std::vector<mfem::Vector> differences( numDirections );
	for ( int trial = 0; trial < numDirections; ++trial )
	{
		mfem::Vector &direction = directions[ static_cast<std::size_t>( trial ) ];
		direction.SetSize( size );
		for ( int i = 0; i < size; ++i )
			direction( i ) = isEssential[ static_cast<std::size_t>( i ) ]
			                 ? 0.0
			                 : std::cos( 2.3*i + 1.1*trial ) + 0.5*std::sin( 0.7*i - 0.3*trial );
		direction /= direction.Norml2();

		mfem::Vector forward( size );
		mfem::Vector backward( size );
		mfem::Vector shifted( state );
		shifted.Add( step, direction );
		residual.Mult( shifted, forward );
		shifted = state;
		shifted.Add( -step, direction );
		residual.Mult( shifted, backward );

		mfem::Vector &difference = differences[ static_cast<std::size_t>( trial ) ];
		difference.SetSize( size );
		difference = forward;
		difference -= backward;
		difference /= 2.0*step;
	}

	mfem::Operator &jacobian = residual.GetGradient( state );
	double worst = 0.0;
	for ( int trial = 0; trial < numDirections; ++trial )
	{
		mfem::Vector applied( size );
		jacobian.Mult( directions[ static_cast<std::size_t>( trial ) ], applied );
		mfem::Vector error( differences[ static_cast<std::size_t>( trial ) ] );
		error -= applied;
		double const relative = error.Norml2()/applied.Norml2();
		std::printf( "  flux law: Jacobian vs central difference, direction %d: "
		             "||J d|| = %.6e, relative error %.3e\n",
		             trial, applied.Norml2(), relative );
		worst = std::max( worst, relative );
	}
	std::fflush( stdout );

	BOOST_TEST( worst < meq::tests::jacobianDifferenceTolerance( 1.0e-8 ),
	            "the flux law's condensed Jacobian differs from a central difference "
	            "of the residual by " << worst << " relative" );
}

/// THREADED ASSEMBLY BESIDE A FLUX LAW REPRODUCES SERIAL, bit for bit, as
/// SolverContract.cpp's threadedAssemblyReproducesSerialAssemblyOnANonlinearSource
/// requires of the source integrator and for its reasons: psi, q, psi* -- the
/// threaded reconstruction lifts the law through one frozen coefficient every
/// thread evaluates -- and the Newton count, which is what a racing gradient
/// moves and nothing else does. omp_get_max_threads() is printed: one thread
/// cannot exhibit a race.
BOOST_AUTO_TEST_CASE( aFluxLawUnderThreadedAssemblyIsTheSerialSolve )
{
#if !defined( MFEM_USE_OPENMP ) || !defined( MFEM_THREAD_SAFE )
	std::printf( "\n  flux law, threaded assembly: skipped, this MFEM has no OpenMP "
	             "or no thread safety\n" );
#else
	using AM = meq::GradShafranovSolver::AssemblyMode;
	std::printf( "\n  flux law, threaded against serial: omp_get_max_threads() = %d%s\n",
	             omp_get_max_threads(), omp_get_max_threads() > 1 ? ""
	             : "  <- ONE THREAD: this run cannot exhibit a race" );

	Fixture const eq = Fixture::flux();
	FluxModelAnisotropy const anisotropy( eq );
	FixtureSource const source( eq );
	mfem::FunctionCoefficient psiCoeff( [ &eq ]( mfem::Vector const &x )
	{
		return eq.psi( x( 0 ), x( 1 ) );
	} );

	for ( int order : { 1, 3 } )
	{
		mfem::Mesh mesh = meq::tests::makeMesh( meq::tests::standardBox(), 16 );
		auto solveIn = [ & ]( AM mode )
		{
			auto solver = std::make_unique<meq::GradShafranovSolver>( mesh, order );
			solver->setAssemblyMode( mode );
			solver->setAnisotropy( anisotropy );
			solver->setSource( source );
			solver->setBoundaryData( psiCoeff );
			solver->solve();
			solver->postProcess();
			return solver;
		};
		auto const serial = solveIn( AM::Serial );
		auto const threaded = solveIn( AM::Threaded );

		double const psi = maxDifference( threaded->potential(), serial->potential() );
		double const flux = maxDifference( threaded->flux(), serial->flux() );
		double const star = maxDifference( threaded->postProcessedPotential(),
		                                   serial->postProcessedPotential() );
		std::printf( "    k = %d, n = 16: psi %.3e, q %.3e, psi* %.3e; Newton %d and %d\n",
		             order, psi, flux, star, serial->newtonIterations(),
		             threaded->newtonIterations() );
		std::fflush( stdout );
		BOOST_TEST( psi == 0.0, "k = " << order << ": psi differs by " << psi );
		BOOST_TEST( flux == 0.0, "k = " << order << ": q differs by " << flux );
		BOOST_TEST( star == 0.0, "k = " << order << ": psi* differs by " << star );
		BOOST_TEST( serial->newtonIterations() == threaded->newtonIterations() );
	}
#endif
}

BOOST_AUTO_TEST_CASE( theFluxModelConvergesAtK1 )
{
	// Measured at h = 0.025: psi 6.60e-04, q 1.26e-03.
	checkFluxOrder( 1, 2.0e-03, 3.8e-03 );
}

BOOST_AUTO_TEST_CASE( theFluxModelConvergesAtK2 )
{
	// Measured at h = 0.025: psi 1.06e-05, q 2.02e-05.
	checkFluxOrder( 2, 3.2e-05, 6.1e-05 );
}

BOOST_AUTO_TEST_CASE( theFluxModelConvergesAtK3 )
{
	// Measured at h = 0.025: psi 1.75e-07, q 3.43e-07.
	checkFluxOrder( 3, 5.3e-07, 1.0e-06 );
}

BOOST_AUTO_TEST_CASE( theFluxModelConvergesAtK4 )
{
	// Measured at h = 0.025: psi 2.44e-09, q 4.76e-09.
	checkFluxOrder( 4, 7.3e-09, 1.4e-08 );
}

/// A zero-amplitude flux model reads psi -- so it takes the flux-law path,
/// FullNL and no linear flux mass -- and must solve the same problem as the
/// isotropic Newton path, to the solver's tolerance rather than to the bit:
/// the two assemble the flux mass through different integrators. And what the
/// flux-law path costs per Newton step, against PotNL's, printed and not
/// asserted, since it is a timing.
BOOST_AUTO_TEST_CASE( aZeroFluxModelIsTheIsotropicNewtonSolve )
{
	Fixture const isotropicEq = Fixture::position( 0.0 );
	Fixture const zeroFlux = Fixture::flux( 0.0 );
	FluxModelAnisotropy const anisotropy( zeroFlux );
	FixtureSource const source( isotropicEq );

	for ( int order : { 1, 3 } )
	{
		mfem::Mesh mesh = meq::tests::makeMesh( meq::tests::standardBox(), 16 );
		mfem::FunctionCoefficient psiCoeff( [ &isotropicEq ]( mfem::Vector const &x )
		{
			return isotropicEq.psi( x( 0 ), x( 1 ) );
		} );

		mfem::StopWatch isotropicClock;
		meq::GradShafranovSolver isotropic( mesh, order );
		isotropic.setSource( source );
		isotropic.setBoundaryData( psiCoeff );
		isotropicClock.Start();
		isotropic.solve();
		isotropicClock.Stop();

		mfem::StopWatch lawClock;
		meq::GradShafranovSolver law( mesh, order );
		law.setAnisotropy( anisotropy );
		law.setSource( source );
		law.setBoundaryData( psiCoeff );
		lawClock.Start();
		law.solve();
		lawClock.Stop();

		double const psiScale = isotropic.potential().Normlinf();
		double const fluxScale = isotropic.flux().Normlinf();
		double const psiDifference = maxDifference( law.potential(), isotropic.potential() );
		double const fluxDifference = maxDifference( law.flux(), isotropic.flux() );

		std::printf( "  zero flux model, k = %d, n = 16: psi differs by %.3e, q by %.3e; "
		             "Newton %d and %d steps; %.3f s and %.3f s, %.2fx per step\n",
		             order, psiDifference/psiScale, fluxDifference/fluxScale,
		             isotropic.newtonIterations(), law.newtonIterations(),
		             isotropicClock.RealTime(), lawClock.RealTime(),
		             ( lawClock.RealTime()/std::max( 1, law.newtonIterations() ) )
		             /( isotropicClock.RealTime()/std::max( 1, isotropic.newtonIterations() ) ) );
		std::fflush( stdout );

		// Two converged solves of one problem: in single each stops at its
		// round-off floor, and they meet to 5 to 174 eps (q, k = 3) where
		// double meets the solver's tolerance. exactAgreementFloor() is zero
		// in double. MEASUREMENTS.md M-ANISO-1.
		double const psiBound = std::max( 1.0e-10*psiScale,
		                                  meq::tests::exactAgreementFloor( psiScale ) );
		double const fluxBound = std::max( 1.0e-10*fluxScale,
		                                   meq::tests::exactAgreementFloor( fluxScale ) );
		BOOST_TEST( psiDifference < psiBound, "k = " << order << ": psi differs by "
		            << psiDifference << " on a scale of " << psiScale );
		BOOST_TEST( fluxDifference < fluxBound, "k = " << order << ": q differs by "
		            << fluxDifference << " on a scale of " << fluxScale );
	}
}

/// What a Delta that reads psi is refused beside, at the first solve: a
/// psi-independent source, which takes the linear path. And an ANSWER past
/// the firehose limit, which cannot be checked before the solve.
BOOST_AUTO_TEST_CASE( aFluxModelsRefusals )
{
	Fixture const eq = Fixture::flux();
	FluxModelAnisotropy const anisotropy( eq );
	FixtureSource const source( eq );
	mfem::FunctionCoefficient psiCoeff( [ &eq ]( mfem::Vector const &x )
	{
		return eq.psi( x( 0 ), x( 1 ) );
	} );
	mfem::FunctionCoefficient sourceCoeff( [ &eq ]( mfem::Vector const &x )
	{
		return eq.f( x( 0 ), x( 1 ) );
	} );

	{
		mfem::Mesh mesh = meq::tests::makeMesh( meq::tests::standardBox(), 4 );
		meq::GradShafranovSolver solver( mesh, 1 );
		solver.setAnisotropy( anisotropy );
		solver.setSource( sourceCoeff );
		solver.setBoundaryData( psiCoeff );
		BOOST_CHECK_THROW( solver.solve(), std::logic_error );
	}
	{
		// Delta = +3 psi^2 reaches 1 wherever |psi| > 0.577, which Example 5's
		// boundary data already has. Newton does not converge on the
		// indefinite operator, and the refusal must name the limit rather than
		// report the iteration.
		Fixture const steep = Fixture::flux( 3.0 );
		FluxModelAnisotropy const firehose( steep );
		mfem::Mesh mesh = meq::tests::makeMesh( meq::tests::standardBox(), 4 );
		meq::GradShafranovSolver solver( mesh, 1 );
		solver.setAnisotropy( firehose );
		solver.setSource( source );
		solver.setBoundaryData( psiCoeff );
		// By its message, so that a Newton failing on the way is not mistaken
		// for the refusal.
		BOOST_CHECK_EXCEPTION( solver.solve(), std::runtime_error,
		                       []( std::runtime_error const &error )
		                       {
			                       return std::string( error.what() ).find( "firehose" )
			                              != std::string::npos;
		                       } );
	}
}

/*
 * ---------------------------------------------------------------------------
 * AN-2b: Delta a function of Psi, through the labels
 * ---------------------------------------------------------------------------
 */
namespace
{
	/// The flux model in the NORMALISED flux, Delta = d Psi^2. On Example 5's box
	/// psi_bnd = 0 and the exact axis value is 1, so at the exact solution
	/// Psi = psi and the manufactured F is the flux model's.
	class NormalisedFluxModelAnisotropy : public meq::Anisotropy
	{
		public:
			explicit NormalisedFluxModelAnisotropy( Fixture const &eqIn )
				: eq( eqIn )
			{
			}

			double delta( double, double, double normalisedPsi ) const override
			{
				return eq.deltaOfPsi( normalisedPsi );
			}

			double dDeltaDPsi( double, double, double normalisedPsi ) const override
			{
				return eq.dDeltaDPsi( normalisedPsi );
			}

			Argument argument() const override
			{
				return Argument::NormalisedFlux;
			}

		private:
			Fixture eq;
	};

	/// The fixture's F as a normalised source. F does not depend on the labels,
	/// and this says so ANALYTICALLY: normalisationDerivatives() answers zero.
	/// That is the trap the solver must not fall into -- the flux law does
	/// depend on psi_ax, and only a differenced column carries it.
	class NormalisedFixtureSource : public meq::NormalisedSource
	{
		public:
			explicit NormalisedFixtureSource( Fixture const &eqIn )
				: eq( eqIn )
			{
			}

			double f( double radius, double z, double ) const override
			{
				return eq.f( radius, z );
			}

			double dFdPsi( double, double, double ) const override
			{
				return 0.0;
			}

			void setNormalisation( double psiAxis, double psiBoundary ) override
			{
				axis = psiAxis;
				boundary = psiBoundary;
			}
			using meq::NormalisedSource::setNormalisation;

			double normalisation() const override
			{
				return axis;
			}

			double boundaryNormalisation() const override
			{
				return boundary;
			}

			bool normalisationDerivatives( double, double, double, double &dFdAxis,
			                               double &dFdBoundary ) const override
			{
				dFdAxis = 0.0;
				dFdBoundary = 0.0;
				return true;
			}

		private:
			Fixture eq;
			double axis = 1.0;
			double boundary = 0.0;
	};

	struct NormalisedFluxMeasurement
	{
		AnisotropicMeasurement point;
		int newtonIterations = 0;
		std::vector<double> history;
		double psiAxis = 0.0;
		int passes = 0;
	};

	/// The best Newton order over triples whose MIDDLE residual is already
	/// below a tenth of the first. meq::tests::bestNewtonOrder() takes every
	/// triple above the round-off floor, and on a bordered solve whose first
	/// steps mostly move psi_ax -- 0.99, then 0.88 -- the ratio of two logs of
	/// numbers near 1 reads 47 to 62, which says nothing.
	double asymptoticNewtonOrder( std::vector<double> const &history )
	{
		double best = 0.0;
		double const floor = meq::newtonTolerance( 1.0e-12 )*history.front();
		for ( std::size_t i = 2; i < history.size(); ++i )
		{
			if ( history[ i - 1 ] > 0.1*history.front() )
				continue;
			if ( history[ i ] < floor || history[ i - 1 ] < floor )
				continue;
			best = std::max( best, meq::tests::newtonOrder( history[ i - 2 ], history[ i - 1 ],
			                                                history[ i ] ) );
		}
		return best;
	}

	/**
	 * In single, the order read through the triple that LANDS on the floor
	 * too, as a lower bound; in double, @a asymptotic unchanged.
	 *
	 * **A RESIDUAL HELD UP BY THE FLOOR CAN ONLY UNDER-STATE THE ORDER**: the
	 * step that landed there would have gone at least as low, so a triple
	 * ending on it reads less than the method's order, never more. That makes
	 * it useless for refuting quadratic convergence -- M-306's 0.36 on a
	 * quadratic tail -- and sound for certifying it. The bordered solve here
	 * has no other triple in single: the one asymptoticNewtonOrder() admits
	 * is the first, which mostly moves psi_ax and reads 0.92. Measured, the
	 * floor-landing triple reads 1.93 to 2.00 at n = 4 and 1.34 to 1.46 at
	 * n = 32, exact; on the source's analytic label column -- the mutation the
	 * case names -- every triple reads 1.00 at best. MEASUREMENTS.md M-ANISO-1.
	 */
	double singleOrderLowerBound( double asymptotic, std::vector<double> const &history )
	{
#ifdef MFEM_USE_SINGLE
		double best = asymptotic;
		for ( std::size_t i = 2; i < history.size(); ++i )
			if ( history[ i - 1 ] <= 0.1*history.front() )
				best = std::max( best, meq::tests::newtonOrder( history[ i - 2 ], history[ i - 1 ],
				                                                history[ i ] ) );
		return best;
#else
		(void)history;
		return asymptotic;
#endif
	}

	NormalisedFluxMeasurement measureNormalisedFlux( Fixture const &eq, int order, int n,
	                                                 bool postProcessedLabels = false )
	{
		meq::tests::Rectangle const box = meq::tests::standardBox();
		mfem::Mesh mesh = meq::tests::makeMesh( box, n );

		mfem::FunctionCoefficient psiCoeff( [ &eq ]( mfem::Vector const &x )
		{
			return eq.psi( x( 0 ), x( 1 ) );
		} );
		mfem::VectorFunctionCoefficient fluxCoeff( 2, [ &eq ]( mfem::Vector const &x,
		                                                       mfem::Vector &value )
		{
			double fluxR = 0.0, fluxZ = 0.0;
			eq.flux( x( 0 ), x( 1 ), fluxR, fluxZ );
			value( 0 ) = fluxR;
			value( 1 ) = fluxZ;
		} );

		NormalisedFluxModelAnisotropy const anisotropy( eq );
		NormalisedFixtureSource source( eq );
		meq::GradShafranovSolver solver( mesh, order );
		solver.setAnisotropy( anisotropy );
		// 0.9 and not the exact 1, so that psi_ax is a genuine unknown the
		// border has to move.
		solver.setSource( source, 0.9 );
		solver.setBoundaryData( psiCoeff );
		// THE AXIS LABEL AT A PRESCRIBED POINT, as a mirror reads its labels,
		// and not the located axis. Delta reads psi_ax, so every error in it is
		// an error in the operator: under the located axis |psi_ax - 1| at k = 1
		// read 1.05e-02, 1.65e-03, 2.49e-03, 3.28e-04 and psi* went 2.88,
		// 0.63, 2.94. The point is Example 5's exact maximum,
		// kr( R + R_0 ) = pi/2, z = 0.
		solver.setAxisFluxPoint( 0.5 + 1.0/( 2.0*1.15 ), 0.0 );
		if ( postProcessedLabels )
			solver.setLabelReading( meq::GradShafranovSolver::LabelReading::PostProcessed );
		solver.solve();
		solver.postProcess();

		NormalisedFluxMeasurement result;
		result.passes = solver.labelReadingPasses();
		AnisotropicMeasurement &point = result.point;
		point.h = box.width()/static_cast<double>( n );
		point.traceDofs = solver.numTraceDofs();
		point.errorPsi = solver.potentialError( psiCoeff );
		point.errorFlux = solver.fluxError( fluxCoeff );
		point.errorPsiStar = solver.postProcessedPotentialError( psiCoeff );
		point.normPsi = solver.exactPotentialNorm( psiCoeff );
		point.normFlux = solver.exactFluxNorm( fluxCoeff );
		result.newtonIterations = solver.newtonIterations();
		result.history = solver.newtonResiduals();
		result.psiAxis = solver.psiAxis();
		return result;
	}
}

/// Delta in Psi on the bordered Newton. The rates as AN-2a's; psi_ax tending to
/// the exact 1; and Newton quadratic, which it is NOT if the label column is
/// the source's analytic zero -- the flux law moves with psi_ax, and only a
/// differenced column sees it.
BOOST_AUTO_TEST_CASE( theNormalisedFluxModelConvergesOnTheBorderedNewton )
{
	Fixture const eq = Fixture::flux();
	for ( int order : { 1, 2, 3 } )
	{
		std::vector<AnisotropicMeasurement> points;
		std::vector<double> axisErrors;
		double best = 0.0;
		for ( int n : meshSizes )
		{
			NormalisedFluxMeasurement const m = measureNormalisedFlux( eq, order, n );
			if ( n == 16 )
			{
				std::printf( "\n  [Psi model] k = %d, n = 16, Newton residuals:", order );
				for ( double r : m.history )
					std::printf( " %.3e", r );
				std::printf( "\n" );
			}
			points.push_back( m.point );
			axisErrors.push_back( std::abs( m.psiAxis - 1.0 ) );
			best = std::max( best, singleOrderLowerBound( asymptoticNewtonOrder( m.history ),
			                                              m.history ) );
			BOOST_TEST( m.newtonIterations >= 2, "k = " << order << ", n = " << n
			            << ": Newton converged in " << m.newtonIterations << " step" );
		}
		printTable( "Psi Delta = -0.5 Psi^2, bordered", order, points );
		std::printf( "  |psi_ax - 1|:" );
		for ( double e : axisErrors )
			std::printf( " %.3e", e );
		std::printf( "; best observed Newton order %.3f\n", best );
		std::fflush( stdout );

		// psi AND q ONLY. psi* is printed and not asserted: Delta reads psi_ax,
		// which is psi_h at a point, so the operator itself carries psi_h's
		// POINTWISE error -- k+1, its constant moving with where the point sits
		// in its element (|psi_ax - 1| at k = 1 reads 9.2e-03, 1.4e-03,
		// 2.4e-03, 3.2e-04) -- and psi* at k+2 cannot beat its own operator.
		// At k = 1 it reads 2.92, 0.63, 2.91.
		// The rates by rateIsMeasurable(), as AN-2a's: in single q at k = 3
		// reads 2.14e-05 and 4.01e-05 on the last two meshes against double's
		// 5.49e-06 and 3.43e-07, the floor the unbordered flux model meets at
		// the same 2.06e-05 and 3.96e-05. MEASUREMENTS.md M-ANISO-1.
		double const expected = order + 1.0 - rateSlack;
		if ( !meq::tests::rateIsMeasurable( order ) )
			BOOST_TEST_MESSAGE( meq::tests::rateNotMeasurableNotice( order ) );
		else
			for ( std::size_t i = 1; i < points.size(); ++i )
			{
				AnisotropicMeasurement const &c = points[ i - 1 ];
				AnisotropicMeasurement const &p = points[ i ];
				double const ratio = c.h/p.h;
				double const ratePsi = meq::tests::rate( c.errorPsi, p.errorPsi, ratio );
				double const rateFlux = meq::tests::rate( c.errorFlux, p.errorFlux, ratio );
				BOOST_TEST( ratePsi >= expected, "k = " << order << ", h = " << p.h
				            << ": psi converged at " << ratePsi << ", wanted " << expected );
				BOOST_TEST( rateFlux >= expected, "k = " << order << ", h = " << p.h
				            << ": q converged at " << rateFlux << ", wanted " << expected );
			}
		BOOST_TEST( axisErrors.back() < axisErrors.front(),
		            "k = " << order << ": psi_ax does not tend to the exact 1" );
		BOOST_TEST( best >= meq::tests::quadraticOrderThreshold( 1.8 ), "k = " << order
		            << ": the best observed Newton order is " << best
		            << "; the label column misses the flux law's dependence on psi_ax" );
	}
}

/// What a Delta in Psi, and any anisotropy beside a prescribed current, is
/// refused beside.
BOOST_AUTO_TEST_CASE( aNormalisedFluxModelsRefusals )
{
	Fixture const eq = Fixture::flux();
	mfem::FunctionCoefficient psiCoeff( [ &eq ]( mfem::Vector const &x )
	{
		return eq.psi( x( 0 ), x( 1 ) );
	} );

	// Psi has no labels to read without a normalised source.
	{
		NormalisedFluxModelAnisotropy const anisotropy( eq );
		FixtureSource const source( eq );
		mfem::Mesh mesh = meq::tests::makeMesh( meq::tests::standardBox(), 4 );
		meq::GradShafranovSolver solver( mesh, 1 );
		solver.setAssemblyMode( meq::GradShafranovSolver::AssemblyMode::Serial );
		solver.setAnisotropy( anisotropy );
		solver.setSource( source );
		solver.setBoundaryData( psiCoeff );
		BOOST_CHECK_THROW( solver.solve(), std::logic_error );
	}
	// The physical psi beside a normalised source is not measured.
	{
		FluxModelAnisotropy const anisotropy( eq );
		NormalisedFixtureSource source( eq );
		mfem::Mesh mesh = meq::tests::makeMesh( meq::tests::standardBox(), 4 );
		meq::GradShafranovSolver solver( mesh, 1 );
		solver.setAnisotropy( anisotropy );
		solver.setSource( source, 0.9 );
		solver.setBoundaryData( psiCoeff );
		BOOST_CHECK_THROW( solver.solve(), std::logic_error );
	}
	// A prescribed current constrains F/R, which is not the toroidal current
	// under any Delta -- a position one included.
	{
		ConstantAnisotropy const anisotropy( -0.2 );
		NormalisedFixtureSource source( eq );
		mfem::Mesh mesh = meq::tests::makeMesh( meq::tests::standardBox(), 4 );
		meq::GradShafranovSolver solver( mesh, 1 );
		solver.setAnisotropy( anisotropy );
		solver.setSource( source, 0.9 );
		solver.setBoundaryData( psiCoeff );
		solver.setPlasmaCurrent( 1.0 );
		BOOST_CHECK_THROW( solver.solve(), std::logic_error );
	}
}

/// THE LABELS FROM psi*. The same problem, the axis label read from psi* at its
/// point by LabelReading::PostProcessed: psi* back at k+2, which the label read
/// from psi_h cannot give it, and the passes settling in a handful.
BOOST_AUTO_TEST_CASE( labelsReadFromPsiStarGivePsiStarItsOrder )
{
	Fixture const eq = Fixture::flux();
	for ( int order : { 1, 2, 3 } )
	{
		std::vector<AnisotropicMeasurement> points;
		std::vector<double> axisErrors;
		std::vector<int> passes;
		for ( int n : meshSizes )
		{
			NormalisedFluxMeasurement const m = measureNormalisedFlux( eq, order, n, true );
			points.push_back( m.point );
			axisErrors.push_back( std::abs( m.psiAxis - 1.0 ) );
			passes.push_back( m.passes );
		}
		printTable( "Psi Delta = -0.5 Psi^2, labels from psi*", order, points );
		std::printf( "  |psi_ax - 1|:" );
		for ( double e : axisErrors )
			std::printf( " %.3e", e );
		std::printf( "; passes:" );
		for ( int p : passes )
			std::printf( " %d", p );
		std::printf( "\n" );
		std::fflush( stdout );

		// The rates by rateIsMeasurable(), as above. In single at k >= 2 the
		// later passes start under the round-off gate and take no step, so
		// the answer is the first pass's to every printed digit: the label's
		// correction is under the format's floor. MEASUREMENTS.md M-ANISO-1.
		double const expected = order + 1.0 - rateSlack;
		double const expectedStar = order + 2.0 - rateSlack;
		if ( !meq::tests::rateIsMeasurable( order ) )
			BOOST_TEST_MESSAGE( meq::tests::rateNotMeasurableNotice( order ) );
		else
			for ( std::size_t i = 1; i < points.size(); ++i )
			{
				AnisotropicMeasurement const &c = points[ i - 1 ];
				AnisotropicMeasurement const &p = points[ i ];
				double const ratio = c.h/p.h;
				BOOST_TEST( meq::tests::rate( c.errorPsi, p.errorPsi, ratio ) >= expected,
				            "k = " << order << ", h = " << p.h << ": psi below k+1" );
				BOOST_TEST( meq::tests::rate( c.errorFlux, p.errorFlux, ratio ) >= expected,
				            "k = " << order << ", h = " << p.h << ": q below k+1" );
				double const rateStar = meq::tests::rate( c.errorPsiStar, p.errorPsiStar, ratio );
				if ( p.errorPsiStar < meq::tests::roundOffFloor( p.normPsi ) )
					BOOST_TEST_MESSAGE( "  k = " << order << ", h = " << p.h
					                    << ": psi* under the round-off floor, rate not asserted" );
				else
					BOOST_TEST( rateStar >= expectedStar, "k = " << order << ", h = " << p.h
					            << ": psi* converged at " << rateStar << ", wanted " << expectedStar
					            << " -- the label is not being read from psi*" );
			}
		for ( int p : passes )
			BOOST_TEST( ( p >= 2 && p <= 6 ), "k = " << order << ": " << p
			            << " passes; the first carries no correction, and more than six"
			               " means the correction is not contracting" );
	}
}

/*
 * ---------------------------------------------------------------------------
 * AN-6a: Delta that reads B
 * ---------------------------------------------------------------------------
 *
 * THE POSITION MODEL POSED AS Delta( B ). A meq::FieldAnisotropy whose
 * numerator is mu0 Delta_p = Delta( R, z ) |grad psi|^2/R^2 at the EXACT
 * solution has Delta = mu0 Delta_p/B^2 = Delta( R, z ) there, so its exact
 * solution and source are the position model's -- and its solve forms Delta
 * from the discrete flux, through meq::fieldDelta(), with d Delta/d u in the
 * law's flux-flux Jacobian. psi* reads the reconstruction's Delta frozen at
 * the solved flux.
 *
 * ON A BOX WITHOUT ITS O-POINT. The standard box holds psi's O-point at
 * ( 0.935, 0 ), where B = 0 and the numerator vanishes as B^2 with it, so
 * Delta there is a ratio of two discretely small numbers. A mirror has no
 * such point inside its plasma. The shifted box R in [ 1.0, 1.75 ] keeps
 * |grad psi| above 0.4.
 */
namespace
{
	meq::tests::Rectangle shiftedBox()
	{
		return meq::tests::Rectangle{ 1.0, 1.75, -0.6, 0.6 };
	}

	/// The position model's Delta, through its numerator at the exact B.
	class FieldModelAnisotropy : public meq::FieldAnisotropy
	{
		public:
			explicit FieldModelAnisotropy( Fixture const &eqIn )
				: eq( eqIn )
			{
			}

			double pressure( double radius, double z, double ) const override
			{
				double dPsiDr = 0.0;
				double dPsiDz = 0.0;
				eq.gradPsi( radius, z, dPsiDr, dPsiDz );
				return eq.delta( radius, z )*( dPsiDr*dPsiDr + dPsiDz*dPsiDz )/( radius*radius );
			}

			double dPressureDPsi( double, double, double ) const override
			{
				return 0.0;
			}

			Argument argument() const override
			{
				return Argument::Position;
			}

		private:
			Fixture eq;
	};

	FluxMeasurement measureField( Fixture const &eq, int order, int n )
	{
		meq::tests::Rectangle const box = shiftedBox();
		mfem::Mesh mesh = meq::tests::makeMesh( box, n );

		mfem::FunctionCoefficient psiCoeff( [ &eq ]( mfem::Vector const &x )
		{
			return eq.psi( x( 0 ), x( 1 ) );
		} );
		mfem::VectorFunctionCoefficient fluxCoeff( 2, [ &eq ]( mfem::Vector const &x,
		                                                       mfem::Vector &value )
		{
			double fluxR = 0.0, fluxZ = 0.0;
			eq.flux( x( 0 ), x( 1 ), fluxR, fluxZ );
			value( 0 ) = fluxR;
			value( 1 ) = fluxZ;
		} );

		FieldModelAnisotropy const anisotropy( eq );
		FixtureSource const source( eq );
		meq::GradShafranovSolver solver( mesh, order );
		solver.setAnisotropy( anisotropy );
		solver.setSource( source );
		solver.setBoundaryData( psiCoeff );
		solver.solve();
		solver.postProcess();

		FluxMeasurement result;
		AnisotropicMeasurement &point = result.point;
		point.h = box.width()/static_cast<double>( n );
		point.traceDofs = solver.numTraceDofs();
		point.errorPsi = solver.potentialError( psiCoeff );
		point.errorFlux = solver.fluxError( fluxCoeff );
		point.errorPsiStar = solver.postProcessedPotentialError( psiCoeff );
		point.normPsi = solver.exactPotentialNorm( psiCoeff );
		point.normFlux = solver.exactFluxNorm( fluxCoeff );
		result.newtonIterations = solver.newtonIterations();
		result.history = solver.newtonResiduals();
		return result;
	}

	/// @param coldOrder whether the best Newton order read off the cold
	///        start's residuals is asserted; the instrument for Newton's
	///        order itself is theFieldLawsNewtonIsQuadraticAtBothSigns.
	void checkFieldOrder( Fixture const &eq, char const *label, int order, double psiCeiling,
	                      double fluxCeiling, bool coldOrder = true )
	{
		std::vector<AnisotropicMeasurement> points;
		std::vector<int> iterations;
		double best = 0.0;
		for ( int n : meshSizes )
		{
			FluxMeasurement const m = measureField( eq, order, n );
			if ( n == 16 )
			{
				std::printf( "\n  [%s] k = %d, n = 16, Newton residuals:", label, order );
				for ( double r : m.history )
					std::printf( " %.3e", r );
				std::printf( "\n" );
			}
			points.push_back( m.point );
			iterations.push_back( m.newtonIterations );
			best = std::max( best, meq::tests::bestNewtonOrder( m.history ) );
		}
		std::printf( "\n  [%s] Newton iterations:", label );
		for ( int it : iterations )
			std::printf( " %d", it );
		std::printf( "; best observed order %.3f\n", best );
		printTable( label, order, points, "R in [ 1.0, 1.75 ]" );

		if ( !meq::tests::rateIsMeasurable( order ) )
			BOOST_TEST_MESSAGE( meq::tests::rateNotMeasurableNotice( order ) );
		else
		{
			double const expected = order + 1.0 - rateSlack;
			double const expectedStar = order + 2.0 - rateSlack;
			for ( std::size_t i = 1; i < points.size(); ++i )
			{
				AnisotropicMeasurement const &c = points[ i - 1 ];
				AnisotropicMeasurement const &p = points[ i ];
				double const ratio = c.h/p.h;
				double const ratePsi = meq::tests::rate( c.errorPsi, p.errorPsi, ratio );
				double const rateFlux = meq::tests::rate( c.errorFlux, p.errorFlux, ratio );
				double const rateStar = meq::tests::rate( c.errorPsiStar, p.errorPsiStar, ratio );
				BOOST_TEST( ratePsi >= expected, label << ", k = " << order << ", h = " << p.h
				            << ": psi converged at " << ratePsi << ", wanted " << expected );
				BOOST_TEST( rateFlux >= expected, label << ", k = " << order << ", h = " << p.h
				            << ": q converged at " << rateFlux << ", wanted " << expected );
				if ( p.errorPsiStar < meq::tests::roundOffFloor( p.normPsi ) )
					BOOST_TEST_MESSAGE( "  k = " << order << ", h = " << p.h
					                    << ": psi* under the round-off floor, rate not asserted" );
				else
					BOOST_TEST( rateStar >= expectedStar, label << ", k = " << order << ", h = "
					            << p.h << ": psi* converged at " << rateStar << ", wanted "
					            << expectedStar );
			}
		}

		for ( int it : iterations )
			BOOST_TEST( it >= 2, label << ", k = " << order << ": Newton converged in " << it
			            << " step, so the law's dependence on B never reached the operator" );
		// IN SINGLE THE COLD ORDER IS NOT READ. Double reads 2.1 on
		// ( 1.95e-03, 3.36e-05, 6.70e-09 ), whose last residual is 50 times
		// under single's floor, and the triples single keeps read 1.62 to 1.72
		// in either format -- and do not discriminate: d Delta/d |q|^2 10% low
		// reads 1.52 to 1.61 there and 10% HIGH 1.75 to 1.86, above the exact.
		// theFieldLawsNewtonIsQuadraticAtBothSigns is the order's instrument,
		// in single as in double. MEASUREMENTS.md M-ANISO-1.
#ifdef MFEM_USE_SINGLE
		if ( coldOrder )
			BOOST_TEST_MESSAGE( "  *** " << label << ", k = " << order << ": SKIPPED the cold"
			                    " Newton order, " << best << ": in single its quadratic triple"
			                    " is under the floor and the triples left read the trajectory."
			                    " The restarts assert the order. M-ANISO-1" );
		coldOrder = false;
#endif
		if ( coldOrder )
			BOOST_TEST( best >= 1.8, label << ", k = " << order << ": the best observed Newton"
			            " order is " << best << "; the law's Jacobian disagrees with its residual" );

		meq::tests::checkCeilings( order, points.back().errorPsi, points.back().normPsi,
		                           psiCeiling, points.back().errorFlux, points.back().normFlux,
		                           fluxCeiling );
	}
}

/// meq::fieldDelta() is the root of Delta w^2 = a( 1 - Delta )^2 that
/// vanishes with a, with its two slopes, on either side of a = 0, at the fold
/// and past it. The slopes against central differences; J = w sqrt( w^2 + 4a ).
BOOST_AUTO_TEST_CASE( fieldDeltaIsTheRootThatVanishesWithThePressure )
{
	for ( double pressure : { -0.24, -0.1, -1.0e-9, 0.0, 1.0e-9, 0.3, 5.0, 1.0e4 } )
	{
		double const fluxSquared = 1.0;
		meq::FieldDelta const d = meq::fieldDelta( pressure, fluxSquared );
		BOOST_TEST_REQUIRE( d.exists, "a = " << pressure );
		double const tension = 1.0 - d.delta;
		BOOST_TEST( std::abs( d.delta*fluxSquared - pressure*tension*tension )
		            < 1.0e-14*( 1.0 + std::abs( pressure ) ), "a = " << pressure
		            << ": Delta = " << d.delta << " is not a root" );
		BOOST_TEST( tension > 0.0, "a = " << pressure );
		if ( pressure > 0.0 )
			BOOST_TEST( ( d.delta > 0.0 && d.delta < 1.0 ), "a = " << pressure );
		if ( pressure < 0.0 )
			BOOST_TEST( ( d.delta < 0.0 && d.delta > -1.0 ), "a = " << pressure );
		// The small-a root is a/w^2 to first order, with nothing cancelled.
		if ( std::abs( pressure ) < 1.0e-6 )
			BOOST_TEST( std::abs( d.delta - pressure ) <= 3.0e-9*std::abs( pressure ) );

		double const h = 1.0e-6*std::max( 1.0, std::abs( pressure ) );
		double const dPressure = ( meq::fieldDelta( pressure + h, fluxSquared ).delta
		                           - meq::fieldDelta( pressure - h, fluxSquared ).delta )/( 2.0*h );
		double const g = 1.0e-6;
		double const dFlux = ( meq::fieldDelta( pressure, fluxSquared + g ).delta
		                       - meq::fieldDelta( pressure, fluxSquared - g ).delta )/( 2.0*g );
		BOOST_TEST( std::abs( dPressure - d.dDeltaDPressure ) < 1.0e-6*( 1.0 + std::abs( dPressure ) ),
		            "a = " << pressure << ": dDelta/da " << d.dDeltaDPressure << " against "
		            << dPressure );
		BOOST_TEST( std::abs( dFlux - d.dDeltaDFluxSquared ) < 1.0e-6*( 1.0 + std::abs( dFlux ) ),
		            "a = " << pressure << ": dDelta/dw^2 " << d.dDeltaDFluxSquared << " against "
		            << dFlux );
		double const jacobian = fluxSquared + 2.0*pressure*tension;
		BOOST_TEST( std::abs( jacobian - std::sqrt( fluxSquared*( fluxSquared + 4.0*pressure ) ) )
		            < 1.0e-12*( 1.0 + std::abs( pressure ) ), "a = " << pressure );
	}

	// At the fold Delta = -1, and past it there is no root: held at -1.
	meq::FieldDelta const fold = meq::fieldDelta( -0.25 + 1.0e-12, 1.0 );
	BOOST_TEST( std::abs( fold.delta + 1.0 ) < 1.0e-5 );
	meq::FieldDelta const past = meq::fieldDelta( -0.3, 1.0 );
	BOOST_TEST( !past.exists );
	BOOST_TEST( past.delta == -1.0 );
	BOOST_TEST( past.dDeltaDPressure == 0.0 );
	BOOST_TEST( past.dDeltaDFluxSquared == 0.0 );
	// No flux and a positive numerator: the firehose limit, held at 0.
	meq::FieldDelta const null = meq::fieldDelta( 0.1, 0.0 );
	BOOST_TEST( !null.exists );
	BOOST_TEST( null.delta == 0.0 );
}

/// The law's Jacobian against a central difference of the residual, as
/// theFluxLawsJacobianMatchesAFiniteDifferenceOfTheResidual has it, for the
/// field model: what it adds is d Delta/d u in the flux-flux block, which the
/// element-local elimination reaches.
BOOST_AUTO_TEST_CASE( theFieldLawsJacobianMatchesAFiniteDifferenceOfTheResidual )
{
	for ( double amplitude : { -0.3, 0.3 } )
	{
	Fixture const eq = Fixture::position( amplitude );
	FieldModelAnisotropy const anisotropy( eq );
	FixtureSource const source( eq );
	mfem::Mesh mesh = meq::tests::makeMesh( shiftedBox(), 3 );
	mfem::FunctionCoefficient psiCoeff( [ &eq ]( mfem::Vector const &x )
	{
		return eq.psi( x( 0 ), x( 1 ) );
	} );

	// AT THE SOLVED ANSWER, PERTURBED BY 0.01, and not at prepare()'s zero
	// trace and zero flux, perturbed by the psi law's 0.3. The condensed
	// residual solves each element's flux row by a local Newton started from
	// the stored flux, and under Delta( B ) that row has no solution where the
	// trace asks for B^2 < -mu0 Delta_p: the root that vanishes with Delta_p
	// ends at the fold, Delta = -1, and the flux cannot carry less. Zero flux
	// is past it, and so is a perturbation of 0.03 on one element of this
	// mesh: there the local solve stalls at 5e-2 and the difference reads
	// 9e-2. At 0, 1e-3 and 1e-2 it reads 1e-10 (measured). So the solver
	// starts from a solved twin's psi, through setInitialGuess(), which seeds
	// the flux, and the state is the twin's trace. The global Newton is
	// NPC's and makes no local solve.
	meq::GradShafranovSolver solved( mesh, 2 );
	solved.setAnisotropy( anisotropy );
	solved.setSource( source );
	solved.setBoundaryData( psiCoeff );
	solved.solve();

	meq::GradShafranovSolver solver( mesh, 2 );
	solver.setAnisotropy( anisotropy );
	solver.setSource( source );
	solver.setBoundaryData( psiCoeff );
	solver.setInitialGuess( solved.potential() );
	solver.prepare();

	mfem::Operator &residual = solver.reducedOperator();
	int const size = solver.reducedSolution().Size();
	std::vector<bool> isEssential( static_cast<std::size_t>( size ), false );
	mfem::Array<int> const &essential = solver.essentialTraceDofs();
	for ( int i = 0; i < essential.Size(); ++i )
		isEssential[ static_cast<std::size_t>( essential[ i ] ) ] = true;
	BOOST_TEST_REQUIRE( essential.Size() > 0 );

	BOOST_TEST_REQUIRE( solved.reducedSolution().Size() == size );
	mfem::Vector state( solved.reducedSolution() );
	for ( int i = 0; i < size; ++i )
		if ( !isEssential[ static_cast<std::size_t>( i ) ] )
			state( i ) += 0.01*std::sin( 1.7*i + 0.4 );

	residual.GetGradient( state );

	double const step = meq::tests::jacobianDifferenceStep( 1.0e-5 );
	int const numDirections = 4;
	double worst = 0.0;
	mfem::Operator &jacobian = residual.GetGradient( state );
	for ( int trial = 0; trial < numDirections; ++trial )
	{
		mfem::Vector direction( size );
		for ( int i = 0; i < size; ++i )
			direction( i ) = isEssential[ static_cast<std::size_t>( i ) ]
			                 ? 0.0
			                 : std::cos( 2.3*i + 1.1*trial ) + 0.5*std::sin( 0.7*i - 0.3*trial );
		direction /= direction.Norml2();

		mfem::Vector forward( size );
		mfem::Vector backward( size );
		mfem::Vector shifted( state );
		shifted.Add( step, direction );
		residual.Mult( shifted, forward );
		shifted = state;
		shifted.Add( -step, direction );
		residual.Mult( shifted, backward );
		mfem::Vector difference( forward );
		difference -= backward;
		difference /= 2.0*step;

		residual.GetGradient( state );
		mfem::Vector applied( size );
		jacobian.Mult( direction, applied );
		mfem::Vector error( difference );
		error -= applied;
		double const relative = error.Norml2()/applied.Norml2();
		std::printf( "  field law, Delta %+.1f: Jacobian vs central difference, direction %d: "
		             "||J d|| = %.6e, relative error %.3e\n",
		             amplitude, trial, applied.Norml2(), relative );
		worst = std::max( worst, relative );
	}
	std::fflush( stdout );

	BOOST_TEST( worst < meq::tests::jacobianDifferenceTolerance( 1.0e-8 ),
	            "Delta " << amplitude << ": the field law's condensed Jacobian differs from a"
	            " central difference of the residual by " << worst << " relative" );
	}
}

BOOST_AUTO_TEST_CASE( theFieldModelConvergesAtK1 )
{
	// Measured at h = 0.0234: psi 4.39e-04, q 7.63e-04.
	checkFieldOrder( Fixture::position(), "field, Delta = -0.3 cos cos", 1, 1.3e-03, 2.3e-03 );
}

BOOST_AUTO_TEST_CASE( theFieldModelConvergesAtK2 )
{
	// Measured at h = 0.0234: psi 3.32e-06, q 6.30e-06.
	checkFieldOrder( Fixture::position(), "field, Delta = -0.3 cos cos", 2, 1.0e-05, 1.9e-05 );
}

BOOST_AUTO_TEST_CASE( theFieldModelConvergesAtK3 )
{
	// Measured at h = 0.0234: psi 3.02e-08, q 6.28e-08.
	checkFieldOrder( Fixture::position(), "field, Delta = -0.3 cos cos", 3, 9.1e-08, 1.9e-07 );
}

BOOST_AUTO_TEST_CASE( theFieldModelConvergesAtK4 )
{
	// Measured at h = 0.0234: psi 2.36e-10, q 5.34e-10.
	checkFieldOrder( Fixture::position(), "field, Delta = -0.3 cos cos", 4, 7.1e-10, 1.6e-09 );
}

/// The firehose side: Delta up to +0.3 through its numerator, where the
/// root's slope in |q|^2 has the other sign. The rates, and not the order read
/// off the cold start: its residuals read 1.71 at best, the quadratic
/// constant growing from 6 to 62 over the last two steps, on PARDISO and
/// UMFPack alike, while the restarts below find Newton quadratic to 0.1% over
/// three decades. So the cold reading is the trajectory's and not the
/// Jacobian's, and the restarts are what assert the order.
BOOST_AUTO_TEST_CASE( thePositiveFieldModelConvergesAtK2 )
{
	// Measured at h = 0.0234: psi 1.91e-06, q 3.14e-06.
	checkFieldOrder( Fixture::position( 0.3 ), "field, Delta = +0.3 cos cos", 2, 5.7e-06, 9.4e-06,
	                 false );
}

/// NEWTON IS QUADRATIC, AT BOTH SIGNS OF Delta: restarted from the solved
/// state with q and psi perturbed by delta = 1e-2, 1e-3, 1e-4, the first
/// step's residual is C r_0^2 with one C to within 1% across the three -- an
/// exact Jacobian's signature, on the NPC path the solve drives, which the
/// condensed difference above does not reach.
BOOST_AUTO_TEST_CASE( theFieldLawsNewtonIsQuadraticAtBothSigns )
{
	for ( double amplitude : { -0.3, 0.3 } )
	{
		Fixture const eq = Fixture::position( amplitude );
		FieldModelAnisotropy const anisotropy( eq );
		FixtureSource const source( eq );
		mfem::Mesh mesh = meq::tests::makeMesh( shiftedBox(), 16 );
		mfem::FunctionCoefficient psiCoeff( [ &eq ]( mfem::Vector const &x )
		{
			return eq.psi( x( 0 ), x( 1 ) );
		} );
		meq::GradShafranovSolver solved( mesh, 2 );
		solved.setAnisotropy( anisotropy );
		solved.setSource( source );
		solved.setBoundaryData( psiCoeff );
		solved.solve();
		int const fieldSize = solved.fluxSpace().GetVSize() + solved.potentialSpace().GetVSize();

		std::vector<double> constants;
		std::vector<double> firstSteps;
		for ( double delta : { 1.0e-2, 1.0e-3, 1.0e-4 } )
		{
			mfem::Vector state( solved.state() );
			for ( int i = 0; i < fieldSize; ++i )
				state( i ) += delta*std::sin( 1.7*i + 0.4 );
			meq::GradShafranovSolver solver( mesh, 2 );
			solver.setAnisotropy( anisotropy );
			solver.setSource( source );
			solver.setBoundaryData( psiCoeff );
			solver.setInitialState( state, {} );
			solver.solve();
			std::vector<double> const &history = solver.newtonResiduals();
			BOOST_TEST_REQUIRE( history.size() >= 2 );
			constants.push_back( history[ 1 ]/( history[ 0 ]*history[ 0 ] ) );
			firstSteps.push_back( history[ 1 ] );
			std::printf( "  field law, Delta %+.1f, restart at delta = %.0e: residuals", amplitude,
			             delta );
			for ( double r : history )
				std::printf( " %.3e", r );
			std::printf( "; r_1/r_0^2 = %.4e\n", constants.back() );
			std::fflush( stdout );
		}
		/*
		 * IN SINGLE THE CONSTANT IS UNDER THE FLOOR AT EVERY RESTART. Double's
		 * C r_0^2 is 1.2e-07 at the LARGEST, under single's floor of about
		 * 3e-07, so all three first steps land on the floor and r_1/r_0^2 is
		 * the floor over r_0^2. What single can read is whether the first step
		 * from 1e-2 reaches the same floor as the one from 1e-4, which an
		 * exact Jacobian does: measured, 1.56 and 1.87 times that floor at the
		 * two signs, and d Delta/d |q|^2 10% high or low 9.6 to 11.5 times, since
		 * its linear part, 1.5e-04 of r_0, is above the floor from 1e-2 and
		 * not from 1e-4. 4 sits between. Double asserts the constant to 1%.
		 * MEASUREMENTS.md M-ANISO-1.
		 */
#ifdef MFEM_USE_SINGLE
		for ( double r : firstSteps )
			BOOST_TEST( r <= 4.0*firstSteps.back(), "Delta " << amplitude << ": a first step"
			            " left " << r << " against " << firstSteps.back() << " from the smallest"
			            " restart, so it did not reach the floor; the law's Jacobian disagrees"
			            " with its residual" );
#else
		for ( double c : constants )
			BOOST_TEST( std::abs( c/constants.back() - 1.0 ) < 1.0e-2, "Delta " << amplitude
			            << ": r_1/r_0^2 = " << c << " against " << constants.back()
			            << " at the smallest restart; Newton is not quadratic" );
#endif
	}
}

/// Threaded assembly and reconstruction beside the field law are the serial
/// solve, bit for bit: the reconstruction reads the solved flux through a host
/// pointer from every thread.
BOOST_AUTO_TEST_CASE( aFieldLawUnderThreadedAssemblyIsTheSerialSolve )
{
#if !defined( MFEM_USE_OPENMP ) || !defined( MFEM_THREAD_SAFE )
	std::printf( "\n  field law, threaded assembly: skipped, this MFEM has no OpenMP "
	             "or no thread safety\n" );
#else
	using AM = meq::GradShafranovSolver::AssemblyMode;
	std::printf( "\n  field law, threaded against serial: omp_get_max_threads() = %d%s\n",
	             omp_get_max_threads(), omp_get_max_threads() > 1 ? ""
	             : "  <- ONE THREAD: this run cannot exhibit a race" );

	Fixture const eq = Fixture::position();
	FieldModelAnisotropy const anisotropy( eq );
	FixtureSource const source( eq );
	mfem::FunctionCoefficient psiCoeff( [ &eq ]( mfem::Vector const &x )
	{
		return eq.psi( x( 0 ), x( 1 ) );
	} );

	for ( int order : { 1, 3 } )
	{
		mfem::Mesh mesh = meq::tests::makeMesh( shiftedBox(), 16 );
		auto solveIn = [ & ]( AM mode )
		{
			auto solver = std::make_unique<meq::GradShafranovSolver>( mesh, order );
			solver->setAssemblyMode( mode );
			solver->setAnisotropy( anisotropy );
			solver->setSource( source );
			solver->setBoundaryData( psiCoeff );
			solver->solve();
			solver->postProcess();
			return solver;
		};
		auto const serial = solveIn( AM::Serial );
		auto const threaded = solveIn( AM::Threaded );

		double const psi = maxDifference( threaded->potential(), serial->potential() );
		double const flux = maxDifference( threaded->flux(), serial->flux() );
		double const star = maxDifference( threaded->postProcessedPotential(),
		                                   serial->postProcessedPotential() );
		std::printf( "    k = %d, n = 16: psi %.3e, q %.3e, psi* %.3e; Newton %d and %d\n",
		             order, psi, flux, star, serial->newtonIterations(),
		             threaded->newtonIterations() );
		std::fflush( stdout );
		BOOST_TEST( psi == 0.0, "k = " << order << ": psi differs by " << psi );
		BOOST_TEST( flux == 0.0, "k = " << order << ": q differs by " << flux );
		BOOST_TEST( star == 0.0, "k = " << order << ": psi* differs by " << star );
		BOOST_TEST( serial->newtonIterations() == threaded->newtonIterations() );
	}
#endif
}

/// What a Delta that reads B is refused: an answer past fieldDelta()'s fold,
/// on the solved flux and by name; beside a psi-independent source; and beside
/// setStabilisationFollowsAnisotropy(), whatever its argument.
BOOST_AUTO_TEST_CASE( aFieldModelsRefusals )
{
	/// A numerator far past the fold everywhere: |q|^2 < -4 mu0 Delta_p.
	class PastTheFold : public meq::FieldAnisotropy
	{
		public:
			double pressure( double, double, double ) const override
			{
				return -1.0e3;
			}

			double dPressureDPsi( double, double, double ) const override
			{
				return 0.0;
			}

			Argument argument() const override
			{
				return Argument::Position;
			}
	};
	PastTheFold const past;
	Fixture const eq = Fixture::position();
	FixtureSource const source( eq );
	mfem::FunctionCoefficient psiCoeff( [ &eq ]( mfem::Vector const &x )
	{
		return eq.psi( x( 0 ), x( 1 ) );
	} );

	{
		mfem::Mesh mesh = meq::tests::makeMesh( shiftedBox(), 4 );
		meq::GradShafranovSolver solver( mesh, 1 );
		solver.setAnisotropy( past );
		solver.setSource( source );
		solver.setBoundaryData( psiCoeff );
		std::string what;
		try
		{
			solver.solve();
		}
		catch ( std::runtime_error const &e )
		{
			what = e.what();
		}
		BOOST_TEST( what.find( "fold" ) != std::string::npos, what );
	}
	{
		mfem::Mesh mesh = meq::tests::makeMesh( shiftedBox(), 4 );
		mfem::ConstantCoefficient zero( 0.0 );
		meq::GradShafranovSolver solver( mesh, 1 );
		solver.setAnisotropy( past );
		solver.setSource( zero );
		solver.setBoundaryData( zero );
		BOOST_CHECK_THROW( solver.solve(), std::logic_error );
	}
	{
		mfem::Mesh mesh = meq::tests::makeMesh( shiftedBox(), 4 );
		meq::GradShafranovSolver solver( mesh, 1 );
		solver.setAnisotropy( past );
		solver.setStabilisationFollowsAnisotropy( true );
		solver.setSource( source );
		solver.setBoundaryData( psiCoeff );
		BOOST_CHECK_THROW( solver.solve(), std::logic_error );
	}
	BOOST_CHECK_THROW( past.delta( 1.0, 0.0, 0.0 ), std::logic_error );
}

