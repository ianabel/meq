#define BOOST_TEST_MODULE RotatingSourceConvergence

#include <boost/test/unit_test.hpp>

#include <cmath>
#include <cstdio>
#include <memory>
#include <vector>

#include "mfem.hpp"

#include "meq/GradShafranov.hpp"
#include "meq/RotatingSource.hpp"
#include "meq/Source.hpp"

#include "analytic/RotatingSoloviev.hpp"
#include "analytic/Soloviev.hpp"

#include "convergence/ConvergenceHarness.hpp"

/*
 * FL-2 OF FLOW-PLAN.md, THROUGH THE SOLVER.
 *
 * RotatingSourceTests.cpp asserts that meq::RotatingSource at omega = 0 is
 * meq::MHDSource pointwise, which is the stronger statement of the two and is
 * where a wrong 4 pi or a wrong sign of Delta* would show. What it cannot show
 * is that a meq::RotatingSource is USABLE: nothing there ever hands one to
 * GradShafranovSolver, and the FL-4 study drives the solver from the analytic
 * fixture rather than from the source class. So a rotating source that computed
 * F perfectly and could not be passed to setSource() would pass every other
 * test in this stage.
 *
 * This file closes that. It builds a rotating source whose profiles are chosen
 * so that at omega = 0 its F is EXACTLY the Solov'ev source -- mu0 p' = -( 1 - A )
 * and g g' = -A, which is refs/HDG-GradShafranov.pdf eq (10) -- and then runs
 * the whole Solov'ev convergence study through it. The rates and the absolute
 * error ceilings are the ones SolovievConvergence.cpp already measured, because
 * it is the same discrete problem reached by a different route.
 *
 * NORMALISED UNITS, mu0 = 1, which meq::RotatingSource takes a constructor
 * argument for. The Solov'ev benchmark is dimensionless, so putting SI constants
 * in and dividing them out again would only add round-off to a comparison whose
 * point is that there is none.
 */

namespace
{

	using meq::tests::standardBox;

	// The NSTX Solov'ev case the rest of the suite uses. F = -( ( 1 - A ) r^2 + A ),
	// so mu0 p' = -( 1 - A ) and g g' = -A.
	double const solovievA = -0.52;

	double const pressureSlope = -( 1.0 - solovievA );
	double const ggPrimeConstant = -solovievA;

	// T_i = T_e = 1/2 so that T_i + T_e = 1 and the reference pressure P0 is the
	// ion density itself. The offset keeps the density positive over the range of
	// psi this benchmark visits, which is O( 0.1 ).
	double const speciesTemperature = 0.5;
	double const densityOffset = 10.0;

	double const referenceRadius = 1.0;

	class AffineProfile : public meq::Profile
	{
		public:
			AffineProfile( double value, double slope ) : v( value ), s( slope ) {}

			double operator()( double psi ) const override { return v + s*psi; }
			double prime( double ) const override { return s; }
			double doublePrime( double ) const override { return 0.0; }

		private:
			double v, s;
	};

	class ConstantMassProfile : public meq::Profile
	{
		public:
			explicit ConstantMassProfile( double value ) : v( value ) {}

			double operator()( double ) const override { return v; }
			double prime( double ) const override { return 0.0; }
			double doublePrime( double ) const override { return 0.0; }

		private:
			double v;
	};

	/// A hydrogenic pair whose total pressure on r = rRef is
	/// P0( psi ) = densityOffset + pressureSlope*psi, so that P0' is the
	/// Solov'ev p'. The masses are equal and opposite in the combination
	/// Z_1 m_2 - Z_2 m_1, which is all the exponent depends on.
	std::vector<meq::Species> solovievSpecies( double mass )
	{
		std::vector<meq::Species> species( 2 );

		species[ 0 ].mass = mass;
		species[ 0 ].charge = 1.0;
		species[ 0 ].temperature = std::make_shared<ConstantMassProfile const>( speciesTemperature );
		species[ 0 ].density = std::make_shared<AffineProfile const>( densityOffset, pressureSlope );

		species[ 1 ].mass = mass*1.0e-4;
		species[ 1 ].charge = -1.0;
		species[ 1 ].temperature = std::make_shared<ConstantMassProfile const>( speciesTemperature );
		species[ 1 ].density = meq::neutralisingDensity( species, 1 );

		return species;
	}

	meq::RotatingSource makeRotatingSource( double omega )
	{
		// mass = 1 and omega chosen by the caller: the exponent coefficient is
		// C = omega^2 ( m_1 + m_2 )/( T_1 + T_2 ), so with these constants
		// C = omega^2 ( 1 + 1e-4 ).
		return meq::RotatingSource( solovievSpecies( 1.0 ),
			omega == 0.0 ? nullptr : std::make_shared<ConstantMassProfile const>( omega ),
			std::make_shared<ConstantMassProfile const>( ggPrimeConstant ),
			referenceRadius,
			1.0 );
	}

	/// A rotating source built to reproduce a RotatingSolovievEquilibrium exactly.
	/// The mapping is forced: with T_1 + T_2 = 1 and constant profiles the
	/// exponent coefficient is C = omega^2 ( m_1 + m_2 ), and the fixture's
	/// exponent is machSquared ( r^2/R0^2 - 1 ) = C ( r^2 - R0^2 )/2, so
	/// C = 2 machSquared/R0^2. P0' is the fixture's p1 because mu0 is 1 here,
	/// and g g' is its F0.
	meq::RotatingSource matching( meq::analytic::RotatingSolovievEquilibrium const & eq )
	{
		double const r0 = eq.getMajorRadius();
		double const massSum = 1.0 + 1.0e-4;
		double const omega = std::sqrt( 2.0*eq.getMachSquared()/( r0*r0*massSum ) );

		std::vector<meq::Species> species( 2 );

		species[ 0 ].mass = 1.0;
		species[ 0 ].charge = 1.0;
		species[ 0 ].temperature = std::make_shared<ConstantMassProfile const>( speciesTemperature );
		species[ 0 ].density = std::make_shared<AffineProfile const>( densityOffset, eq.getP1() );

		species[ 1 ].mass = 1.0e-4;
		species[ 1 ].charge = -1.0;
		species[ 1 ].temperature = std::make_shared<ConstantMassProfile const>( speciesTemperature );
		species[ 1 ].density = meq::neutralisingDensity( species, 1 );

		return meq::RotatingSource( species,
			eq.getMachSquared() == 0.0 ? nullptr : std::make_shared<ConstantMassProfile const>( omega ),
			std::make_shared<ConstantMassProfile const>( eq.getF0() ),
			r0, 1.0 );
	}

	struct Measurement
	{
		double h;
		int traceDofs;
		double errorPsi;
		double errorFlux;
		int newtonIterations;
	};

	Measurement measureThroughSource( meq::Source const & source,
		meq::analytic::SolovievEquilibrium const & eq, int order, int n )
	{
		mfem::Mesh mesh = meq::tests::makeMesh( standardBox(), n );

		mfem::FunctionCoefficient psiCoeff( [ &eq ]( mfem::Vector const & x )
		{
			return eq.psi( x( 0 ), x( 1 ) );
		} );
		mfem::VectorFunctionCoefficient fluxCoeff( 2, [ &eq ]( mfem::Vector const & x, mfem::Vector & value )
		{
			eq.flux( x( 0 ), x( 1 ), value( 0 ), value( 1 ) );
		} );

		meq::GradShafranovSolver solver( mesh, order );
		solver.setSource( source );
		solver.setBoundaryData( psiCoeff );
		solver.solve();

		Measurement point;
		point.h = standardBox().width()/static_cast<double>( n );
		point.traceDofs = solver.numTraceDofs();
		point.errorPsi = solver.potentialError( psiCoeff );
		point.errorFlux = solver.fluxError( fluxCoeff );
		point.newtonIterations = solver.newtonIterations();
		return point;
	}

	double rate( double coarse, double fine, double ratio )
	{
		return std::log( coarse/fine )/std::log( ratio );
	}

	void printTable( int order, std::vector<Measurement> const & points )
	{
		std::printf( "\n  Solov'ev NSTX driven by meq::RotatingSource at omega = 0, k = %d\n", order );
		std::printf( "  %8s %9s %14s %7s %14s %7s %7s\n",
		             "h", "trace", "L2(psi)", "rate", "L2(q)", "rate", "Newton" );
		for ( std::size_t i = 0; i < points.size(); ++i )
		{
			Measurement const & p = points[ i ];
			if ( i == 0 )
			{
				std::printf( "  %8.5f %9d %14.6e %7s %14.6e %7s %7d\n",
				             p.h, p.traceDofs, p.errorPsi, "-", p.errorFlux, "-", p.newtonIterations );
			}
			else
			{
				double const ratio = points[ i - 1 ].h/p.h;
				std::printf( "  %8.5f %9d %14.6e %7.3f %14.6e %7.3f %7d\n",
				             p.h, p.traceDofs,
				             p.errorPsi, rate( points[ i - 1 ].errorPsi, p.errorPsi, ratio ),
				             p.errorFlux, rate( points[ i - 1 ].errorFlux, p.errorFlux, ratio ),
				             p.newtonIterations );
			}
		}
		std::fflush( stdout );
	}

	double const rateSlack = 0.15;

	void checkOrder( int order, double psiCeiling, double fluxCeiling )
	{
		meq::analytic::SolovievEquilibrium const eq = meq::analytic::SolovievEquilibrium::nstx();
		meq::RotatingSource const source = makeRotatingSource( 0.0 );

		std::vector<Measurement> points;
		for ( int n : meq::tests::dyadicMeshes() )
			points.push_back( measureThroughSource( source, eq, order, n ) );

		printTable( order, points );

		double const expected = order + 1.0 - rateSlack;

		for ( std::size_t i = 1; i < points.size(); ++i )
		{
			double const ratio = points[ i - 1 ].h/points[ i ].h;
			double const ratePsi = rate( points[ i - 1 ].errorPsi, points[ i ].errorPsi, ratio );
			double const rateFlux = rate( points[ i - 1 ].errorFlux, points[ i ].errorFlux, ratio );

			BOOST_TEST( ratePsi >= expected,
			            "k = " << order << ", h = " << points[ i ].h
			            << ": psi converged at " << ratePsi << ", wanted " << expected );
			BOOST_TEST( rateFlux >= expected,
			            "k = " << order << ", h = " << points[ i ].h
			            << ": q converged at " << rateFlux << ", wanted " << expected );
		}

		BOOST_TEST( points.back().errorPsi < psiCeiling,
		            "k = " << order << ": L2 error in psi is " << points.back().errorPsi
		            << ", above the ceiling " << psiCeiling );
		BOOST_TEST( points.back().errorFlux < fluxCeiling,
		            "k = " << order << ": L2 error in q is " << points.back().errorFlux
		            << ", above the ceiling " << fluxCeiling );
	}

}

/// The source before the solver. If these two disagree pointwise then every
/// number below is measuring the wrong equation, and the rates would not say so
/// -- a wrong sign convention converges, at the right rate, to the wrong
/// function.
BOOST_AUTO_TEST_CASE( atRestTheRotatingSourceIsTheSolovievSource )
{
	meq::RotatingSource const rotating = makeRotatingSource( 0.0 );
	meq::SolovievSource const soloviev( solovievA );

	meq::tests::Rectangle const box = standardBox();

	for ( double r = box.rMin; r <= box.rMax + 1.0e-12; r += 0.05 )
	{
		for ( double psi : { -0.3, -0.1, 0.0, 0.1, 0.25 } )
		{
			double const expected = soloviev.f( r, 0.0, psi );
			double const actual = rotating.f( r, 0.0, psi );

			BOOST_TEST( std::fabs( actual - expected ) <= 1.0e-13*std::max( 1.0, std::fabs( expected ) ),
			            "at r = " << r << ", psi = " << psi << ": the rotating source gives F = "
			            << actual << " where the Solov'ev source gives " << expected );
			BOOST_TEST( std::fabs( rotating.dFdPsi( r, 0.0, psi ) ) <= 1.0e-13,
			            "at r = " << r << ", psi = " << psi << ": dF/dpsi is "
			            << rotating.dFdPsi( r, 0.0, psi ) << " where the Solov'ev problem is affine in psi" );
		}
	}
}

/// Rotation must actually reach the source, or the omega = 0 agreement above is
/// the only thing this class ever does. The exponent is chosen to be O( 1 ) at
/// the outboard edge, and the check is on the SHAPE -- outboard against inboard
/// at the same psi -- because that is what centrifugal confinement is.
BOOST_AUTO_TEST_CASE( rotationChangesTheSourceOutboardAgainstInboard )
{
	meq::RotatingSource const rotating = makeRotatingSource( 1.0 );
	meq::tests::Rectangle const box = standardBox();

	double const outboard = rotating.pressure( box.rMax, 0.0 );
	double const inboard = rotating.pressure( box.rMin, 0.0 );

	BOOST_TEST( outboard > 2.0*inboard,
	            "the pressure is " << outboard << " outboard against " << inboard
	            << " inboard, so the centrifugal term is not reaching the source" );
	BOOST_TEST( std::fabs( rotating.potential( referenceRadius, 0.0 ) ) == 0.0,
	            "phi_0 does not vanish on the reference curve, so the gauge is not the one documented" );
}

/*
 * THE RATES. Same ceilings as SolovievConvergence.cpp, because this is the same
 * discrete problem: if driving the solver from meq::RotatingSource rather than
 * from a FunctionCoefficient changed the answer at all, these would move.
 */

BOOST_AUTO_TEST_CASE( orderOneConvergesAtTwo )
{
	// SolovievConvergence.cpp measures psi 3.597e-05, q 6.084e-05 at h = 0.025.
	checkOrder( 1, 1.1e-4, 1.9e-4 );
}

BOOST_AUTO_TEST_CASE( orderTwoConvergesAtThree )
{
	// Measured there: psi 1.916e-07, q 2.510e-07.
	checkOrder( 2, 5.8e-7, 7.6e-7 );
}

BOOST_AUTO_TEST_CASE( orderThreeConvergesAtFour )
{
	// Measured there: psi 5.510e-10, q 8.816e-10.
	checkOrder( 3, 1.7e-9, 2.7e-9 );
}

/*
 * THE INTEGRATION CHECK, AND IT IS THE ONE THAT CLOSES THE STAGE.
 *
 * RotatingSoloviev.hpp is a closed form verified against Delta* by central
 * differences; meq::RotatingSource is an independent implementation of the same
 * physics from RoPP (96), (97) and (136). Nothing so far has required them to
 * agree -- the FL-4 study drives the solver from the fixture's own f(), and the
 * study above drives it from the source at omega = 0. If the two disagree under
 * rotation then one of them is wrong and every other test in this stage would
 * still pass.
 */
BOOST_AUTO_TEST_CASE( theSourceAndTheFixtureAgreeUnderRotation )
{
	meq::tests::Rectangle const box = standardBox();

	for ( auto const & eq : { meq::analytic::RotatingSolovievEquilibrium::stationary(),
	                          meq::analytic::RotatingSolovievEquilibrium::rotating(),
	                          meq::analytic::RotatingSolovievEquilibrium::fastRotating() } )
	{
		meq::RotatingSource const source = matching( eq );

		for ( double r = box.rMin; r <= box.rMax + 1.0e-12; r += 0.05 )
		{
			for ( double z : { -0.4, 0.0, 0.4 } )
			{
				double const expected = eq.f( r, z, eq.psi( r, z ) );
				double const actual = source.f( r, z, eq.psi( r, z ) );

				BOOST_TEST( std::fabs( actual - expected ) <= 1.0e-12*std::max( 1.0, std::fabs( expected ) ),
				            "machSquared = " << eq.getMachSquared() << " at r = " << r
				            << ": meq::RotatingSource gives F = " << actual
				            << " where the fixture gives " << expected );
			}
		}
	}
}

/// And the same, all the way through the solver: drive the rotating benchmark
/// from meq::RotatingSource rather than from the fixture, and require the same
/// k+1. A source that agreed pointwise but could not be assembled would pass
/// the check above and fail this one.
BOOST_AUTO_TEST_CASE( theSourceDrivesTheRotatingBenchmarkAtDesignOrder )
{
	meq::analytic::RotatingSolovievEquilibrium const eq
		= meq::analytic::RotatingSolovievEquilibrium::rotating();
	meq::RotatingSource const source = matching( eq );

	mfem::FunctionCoefficient psiCoeff( [ &eq ]( mfem::Vector const & x )
	{
		return eq.psi( x( 0 ), x( 1 ) );
	} );
	mfem::VectorFunctionCoefficient fluxCoeff( 2, [ &eq ]( mfem::Vector const & x, mfem::Vector & value )
	{
		eq.flux( x( 0 ), x( 1 ), value( 0 ), value( 1 ) );
	} );

	std::vector<double> errors;
	std::vector<double> sizes;

	std::printf( "\n  Rotating Solov'ev driven by meq::RotatingSource, k = 2\n" );
	std::printf( "  %8s %14s %7s %7s\n", "h", "L2(psi)", "rate", "Newton" );

	for ( int n : meq::tests::dyadicMeshes() )
	{
		mfem::Mesh mesh = meq::tests::makeMesh( standardBox(), n );
		meq::GradShafranovSolver solver( mesh, 2 );
		solver.setSource( source );
		solver.setBoundaryData( psiCoeff );
		solver.solve();

		double const h = standardBox().width()/static_cast<double>( n );
		double const e = solver.potentialError( psiCoeff );

		if ( errors.empty() )
			std::printf( "  %8.5f %14.6e %7s %7d\n", h, e, "-", solver.newtonIterations() );
		else
			std::printf( "  %8.5f %14.6e %7.3f %7d\n", h, e,
			             rate( errors.back(), e, sizes.back()/h ), solver.newtonIterations() );

		errors.push_back( e );
		sizes.push_back( h );
	}
	std::fflush( stdout );

	for ( std::size_t i = 1; i < errors.size(); ++i )
	{
		double const r = rate( errors[ i - 1 ], errors[ i ], sizes[ i - 1 ]/sizes[ i ] );
		BOOST_TEST( r >= 3.0 - rateSlack,
		            "h = " << sizes[ i ] << ": psi converged at " << r << " driven from the source, "
		            "wanted " << 3.0 - rateSlack );
	}
}

/*
 * THE GENERAL CLOSURE, THROUGH A REAL SOLVE.
 *
 * RotatingSourceTests exercises Closure::RootFind thoroughly, but only ever at
 * points a test chose. A solve is different in one way that matters: the root
 * find runs at every quadrature point of every element on every Newton
 * iteration, at whatever psi the iterate has wandered to, and it THROWS if it
 * fails. So a closure that were merely fragile -- rather than wrong -- would
 * pass every unit test in this stage and take a production run down from inside
 * a quadrature loop.
 *
 * The impurity here is a trace one, 1e-6 of the electron density, so the answer
 * must be the two-species one to within that. Two things are therefore asserted
 * at once: that the general path survives the solve, and that it agrees with the
 * closed form in the limit where they must agree.
 *
 * THE NUMBER COMPARED IS A YARDSTICK AND NOT AN ERROR, WHICH IS WORTH SAYING
 * PLAINLY. The source rotates, so it does not solve the static Solov'ev problem
 * and its L2 distance from that solution is large and means nothing on its own.
 * What it is good for is being computed IDENTICALLY for both closures, so that
 * the difference between the two numbers is a difference between the two
 * solutions and nothing else. Reading it as an error would be reading 6.6e-03
 * as a convergence failure; it is the rotation.
 */
BOOST_AUTO_TEST_CASE( theGeneralClosureSurvivesASolveAndAgreesInTheTraceLimit )
{
	double const impurityFraction = 1.0e-6;

	// Deuterium, a trace of fully stripped carbon, and electrons closed by
	// neutrality. Same total pressure slope as the two-species case, so the
	// equilibrium is the same one to O( impurityFraction ).
	std::vector<meq::Species> three( 3 );

	three[ 0 ].mass = 1.0;
	three[ 0 ].charge = 1.0;
	three[ 0 ].temperature = std::make_shared<ConstantMassProfile const>( speciesTemperature );
	three[ 0 ].density = std::make_shared<AffineProfile const>( densityOffset, pressureSlope );

	three[ 1 ].mass = 12.0;
	three[ 1 ].charge = 6.0;
	three[ 1 ].temperature = std::make_shared<ConstantMassProfile const>( speciesTemperature );
	three[ 1 ].density = std::make_shared<AffineProfile const>( impurityFraction*densityOffset, 0.0 );

	three[ 2 ].mass = 1.0e-4;
	three[ 2 ].charge = -1.0;
	three[ 2 ].temperature = std::make_shared<ConstantMassProfile const>( speciesTemperature );
	three[ 2 ].density = meq::neutralisingDensity( three, 2 );

	meq::RotatingSource const general( three, std::make_shared<ConstantMassProfile const>( 1.0 ),
		std::make_shared<ConstantMassProfile const>( ggPrimeConstant ), referenceRadius, 1.0 );
	BOOST_CHECK_MESSAGE( general.closure() == meq::RotatingSource::Closure::RootFind,
	                     "three species should not have taken the closed-form path" );

	meq::RotatingSource const closed = makeRotatingSource( 1.0 );

	meq::analytic::SolovievEquilibrium const eq = meq::analytic::SolovievEquilibrium::nstx();

	// The solve itself. A root-find failure anywhere in it throws, so reaching
	// the assertions at all is half of what this test is for.
	Measurement generalPoint{};
	BOOST_REQUIRE_NO_THROW( generalPoint = measureThroughSource( general, eq, 2, 16 ) );
	Measurement const closedPoint = measureThroughSource( closed, eq, 2, 16 );

	std::printf( "\n  three species through the solver, k = 2, n = 16\n" );
	std::printf( "    (the L2 column is a yardstick against the STATIC Solov'ev solution,\n" );
	std::printf( "     not an error: these sources rotate. Only the difference is meaningful.)\n" );
	std::printf( "    root find  L2 %14.6e, Newton %d\n", generalPoint.errorPsi, generalPoint.newtonIterations );
	std::printf( "    closed     L2 %14.6e, Newton %d\n", closedPoint.errorPsi, closedPoint.newtonIterations );
	std::fflush( stdout );

	// The impurity is trace, so the two must agree to about its fraction. Not to
	// round-off: the carbon really is there, and a test asserting it made no
	// difference at all would be asserting that the third species is ignored.
	double const relative = std::fabs( generalPoint.errorPsi - closedPoint.errorPsi )
		/std::max( closedPoint.errorPsi, 1.0e-300 );

	BOOST_TEST( relative < 1.0e-3,
	            "the trace-impurity solve differs from the two-species one by " << relative
	            << ", far more than the impurity fraction " << impurityFraction
	            << ", so the general closure is not reducing to the closed form" );
	BOOST_TEST( relative > 0.0,
	            "the trace-impurity solve is bit-identical to the two-species one, so the "
	            "third species is being ignored rather than solved for" );
}
