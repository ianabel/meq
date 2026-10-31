#define BOOST_TEST_MODULE AnisotropicConvergence
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <vector>

#include "mfem.hpp"

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
 * AN-2 adds the solves of the flux model.
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
	class FixtureAnisotropy : public meq::Anisotropy
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

	class ConstantAnisotropy : public meq::Anisotropy
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

			BOOST_TEST( psiDifference < 1.0e-12*psiScale,
			            "Delta = " << delta << ", k = " << order << ": psi differs by "
			            << psiDifference << " on a scale of " << psiScale );
			BOOST_TEST( fluxDifference < 1.0e-12*fluxScale,
			            "Delta = " << delta << ", k = " << order << ": q - ( 1 - Delta ) q_0 is "
			            << fluxDifference << " on a scale of " << fluxScale );
		}
	}
}

namespace
{
	std::vector<int> const meshSizes = { 4, 8, 16, 32 };
	double const rateSlack = 0.15;

	void printTable( int order, std::vector<AnisotropicMeasurement> const &points )
	{
		std::printf( "\n  manufactured anisotropic, position Delta, k = %d, standard box\n", order );
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
		printTable( order, points );

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
	class CornerFirehose : public meq::Anisotropy
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
