#define BOOST_TEST_MODULE SafetyFactorSolver
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <functional>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include "mfem.hpp"

#include "meq/FluxExtraction.hpp"
#include "meq/GradShafranov.hpp"
#include "meq/Profiles.hpp"
#include "meq/SafetyFactor.hpp"
#include "meq/SafetyFactorSolve.hpp"
#include "meq/SafetyFactorTangent.hpp"
#include "meq/Source.hpp"

#include "convergence/ConvergenceHarness.hpp"

/*
 * DRIVING AN EQUILIBRIUM BY `q( psi )` -- ROADMAP.md item 10, the solver half,
 * through the library: meq::SafetyFactorEquilibrium solves and differentiates,
 * meq::SafetyFactorSystem poses q at the knots and the pin, and
 * meq::solveLeastSquares iterates. SafetyFactorTests covers the MFEM-free
 * half -- the spline, the rows and the least squares on synthetic problems.
 *
 * THE ACCEPTANCE IS A ROUND TRIP WITH A KNOWN ANSWER. The reference
 * equilibrium is built from `g = g0 + g1 Psi`, whose `g^2` is quadratic and so
 * reproduced exactly by the knots' not-a-knot cubic; its `q` at the knots, on
 * the same mesh and by the same extraction, is the target. The loop is then
 * started from a `g` 40% larger -- a different equilibrium, not a
 * perturbation -- pinned by the reference's own `g` on the boundary, and has
 * to come back. A loop that merely reached A fixed point would satisfy a
 * self-consistency check while sitting on the wrong one.
 */
namespace
{
	using meq::GradShafranovSolver;

	double const g0 = 2.20;
	double const g1 = 0.55;

	double referenceGSquared( double psi )
	{
		return ( g0 + g1*psi )*( g0 + g1*psi );
	}

	/// One solver on the standard box, normalised MHD source, the cubic's
	/// gg' of @a field installed by the caller.
	struct Machine
	{
		std::unique_ptr<mfem::Mesh> mesh;
		std::unique_ptr<meq::NormalisedMHDSource> source;
		std::unique_ptr<GradShafranovSolver> solver;
		std::unique_ptr<mfem::FunctionCoefficient> guess;
		mfem::ConstantCoefficient zero{ 0.0 };
		std::unique_ptr<meq::SafetyFactorEquilibrium> equilibrium;

		Machine( int cells, int order )
		{
			meq::tests::Rectangle const box = meq::tests::standardBox();
			mesh = std::make_unique<mfem::Mesh>( meq::tests::makeMesh( box, cells ) );
			auto pPrime = std::make_shared<meq::ConstantProfile const>( 0.45 );
			auto flat = std::make_shared<meq::ConstantProfile const>( 0.0 );
			source = std::make_unique<meq::NormalisedMHDSource>( pPrime, flat, 0.30, 1.0 );
			solver = std::make_unique<GradShafranovSolver>( *mesh, order );
			guess = std::make_unique<mfem::FunctionCoefficient>( [ box ]( mfem::Vector const &x )
			{
				return 0.30*std::sin( M_PI*( x( 0 ) - box.minRadius )/box.width() )
				       *std::sin( M_PI*( x( 1 ) - box.zMin )/box.height() );
			} );
			solver->setSource( *source, 0.30 );
			solver->setBoundaryData( zero );
			solver->setInitialGuess( *guess );
			solver->setNewtonControl( 1.0e-12, 1.0e-14, 60 );

			std::size_t const knots = meq::toroidalKnotCount( meq::plasmaCells( *mesh ), order );
			meq::ToroidalKnotCut const cut = meq::toroidalKnotCut( knots );
			meq::SafetyFactorEquilibrium::Options options;
			options.family.surfaces = cut.surfaces;
			options.family.innerCut = cut.innerCut;
			options.family.outerCut = cut.outerCut;
			options.pin = meq::ScalePin::EdgeField;
			equilibrium = std::make_unique<meq::SafetyFactorEquilibrium>( *solver, *source, options );
		}

		meq::ToroidalFieldSpline splineOf( std::function<double( double )> const &gSquared ) const
		{
			meq::ToroidalFieldSpline spline( equilibrium->labels() );
			std::vector<double> values;
			for ( double label : spline.labels() )
				values.push_back( gSquared( 1.0 - label ) );
			spline.setValues( values );
			return spline;
		}
	};

	/// The target: q at the knots of the closed form's own equilibrium.
	std::vector<double> referenceTarget( Machine &machine )
	{
		meq::ToroidalFieldSpline const exact = machine.splineOf( referenceGSquared );
		meq::SafetyFactorEvaluation evaluation;
		BOOST_TEST_REQUIRE( machine.equilibrium->evaluate( exact, evaluation ) );
		std::vector<double> q;
		for ( std::size_t i = 0; i < evaluation.qOverG.size(); ++i )
			q.push_back( evaluation.qOverG[ i ]*std::sqrt( exact.values()[ i ] ) );
		return q;
	}

	meq::SafetyFactorSystem systemFor( Machine &machine, std::vector<double> const &knotTarget,
	                                   double tolerance )
	{
		std::vector<double> const labels = machine.equilibrium->labels();
		auto target = [ labels, knotTarget ]( double psi )
		{
			std::size_t best = 0;
			for ( std::size_t i = 1; i < labels.size(); ++i )
				if ( std::abs( 1.0 - labels[ i ] - psi ) < std::abs( 1.0 - labels[ best ] - psi ) )
					best = i;
			return knotTarget[ best ];
		};
		meq::SafetyFactorSystem::Pin pin;
		pin.kind = meq::ScalePin::EdgeField;
		pin.value = g0;
		meq::SafetyFactorEquilibrium *equilibrium = machine.equilibrium.get();
		return meq::SafetyFactorSystem(
			labels, target, pin,
			[ equilibrium ]( meq::ToroidalFieldSpline const &s, meq::SafetyFactorEvaluation &out )
			{
				return equilibrium->evaluate( s, out );
			},
			[ equilibrium ]( meq::ToroidalFieldSpline const &s, meq::SafetyFactorDerivative &out )
			{
				return equilibrium->differentiate( s, out );
			},
			tolerance );
	}
}

/*
 * THE LOOP RECOVERS THE FIELD ITS TARGET WAS BUILT FROM, from a g 40% too
 * large, and quadratically: the Jacobian is the linearised solve's exact
 * tangent (M-QR-1), so the removable residual falls by orders a step once
 * near, not by a constant.
 */
BOOST_AUTO_TEST_CASE( theLoopRecoversTheFieldItsTargetWasBuiltFrom )
{
	int const cells = 16, order = 2;
	Machine machine( cells, order );
	std::vector<double> const target = referenceTarget( machine );
	meq::SafetyFactorSystem system = systemFor( machine, target, meq::newtonTolerance( 1.0e-12 ) );

	std::vector<double> start;
	for ( double label : machine.equilibrium->labels() )
		start.push_back( 1.96*referenceGSquared( 1.0 - label ) );

	meq::LeastSquaresOptions options;
	options.tolerance = 1.0e-10;
	options.floor = [ &system ]() { return system.floor(); };
	options.constraints = 1;
	options.level = meq::SafetyFactorSystem::scaleGradient;
	options.epsilon = static_cast<double>( std::numeric_limits<mfem::real_t>::epsilon() );
	meq::LeastSquaresResult const result = meq::solveLeastSquares(
		start,
		[ &system ]( std::vector<double> const &u, std::vector<double> &r ) { return system.residual( u, r ); },
		[ &system ]( std::vector<double> const &u, std::vector<double> &j ) { return system.jacobian( u, j ); },
		options );

	std::printf( "\n  THE LOOP AGAINST A KNOWN g: %zu knots, %s\n", start.size(), result.status.c_str() );
	std::printf( "    %d iterations, %d equilibria, misfit %.3e, conditioning %.3e; removable by iteration:",
	             result.iterations, system.solves(), result.misfit, result.conditioning );
	for ( double h : result.history )
		std::printf( " %.1e", h );
	std::printf( "\n" );
	BOOST_TEST_REQUIRE( result.converged, result.status );

	double worst = 0.0;
	for ( std::size_t i = 0; i < result.x.size(); ++i )
	{
		double const exact = referenceGSquared( 1.0 - machine.equilibrium->labels()[ i ] );
		worst = std::max( worst, std::abs( result.x[ i ] - exact )/exact );
	}
	std::printf( "    worst |g^2/g^2_exact - 1| at the knots %.3e\n", worst );
	BOOST_TEST( worst < 1.0e-7,
		"the loop did not recover the closed form its target was measured from: "
		<< worst << " at worst. The target is exact at every knot and the spline "
		"represents the closed form exactly, so a converged loop that lands away "
		"from it converged to the wrong fixed point" );
	BOOST_TEST( result.iterations <= 10,
		"the loop took " << result.iterations << " iterations: a crawl that "
		"converges is a defect report, and the Jacobian is meant to be exact" );
}

/*
 * THE TANGENT JACOBIAN IS THE DIFFERENCED ONE (M-QR-1): every column of
 * d( q/g )/du from SafetyFactorEquilibrium::differentiate() against central
 * differences of re-solved equilibria, Richardson-extrapolated, at the
 * closed form.
 */
BOOST_AUTO_TEST_CASE( theTangentJacobianIsTheDifferencedOne )
{
	int const cells = 8, order = 2;
	Machine machine( cells, order );
	meq::ToroidalFieldSpline const exact = machine.splineOf( referenceGSquared );
	meq::SafetyFactorEvaluation base;
	BOOST_TEST_REQUIRE( machine.equilibrium->evaluate( exact, base ) );
	meq::SafetyFactorDerivative tangent;
	BOOST_TEST_REQUIRE( machine.equilibrium->differentiate( exact, tangent ) );
	std::size_t const m = exact.size();

	auto qOverGAt = [ & ]( std::size_t j, double h )
	{
		meq::ToroidalFieldSpline moved( exact );
		std::vector<double> values = exact.values();
		values[ j ] += h;
		moved.setValues( values );
		meq::SafetyFactorEvaluation out;
		BOOST_TEST_REQUIRE( machine.equilibrium->evaluate( moved, out ) );
		return out.qOverG;
	};
	double worst = 0.0, largest = 0.0;
	for ( std::size_t j = 0; j < m; ++j )
	{
		double const h = 1.0e-6*exact.values()[ j ];
		std::vector<double> const p2 = qOverGAt( j, 2.0*h ), m2 = qOverGAt( j, -2.0*h );
		std::vector<double> const p1 = qOverGAt( j, h ), m1 = qOverGAt( j, -h );
		for ( std::size_t i = 0; i < m; ++i )
		{
			double const d2 = ( p2[ i ] - m2[ i ] )/( 4.0*h );
			double const d1 = ( p1[ i ] - m1[ i ] )/( 2.0*h );
			double const differenced = ( 4.0*d1 - d2 )/3.0;
			largest = std::max( largest, std::abs( differenced ) );
			worst = std::max( worst, std::abs( tangent.qOverG[ i*m + j ] - differenced ) );
		}
	}
	std::printf( "\n  THE TANGENT JACOBIAN: %zu knots, max | J_tangent - J_resolved | %.3e of max | J | %.3e, "
	             "relative %.3e\n", m, worst, largest, worst/largest );
	BOOST_TEST( worst < 1.0e-6*largest,
		"the tangent Jacobian is " << worst/largest << " of its size from re-solved "
		"differences: a term of the map is missing from it, as the ray origin's "
		"motion was (M-QR-1)" );

	/*
	 * AND THE PIN ROWS: the delivered flux and current linearised as the
	 * quadratures they are read with, against the same re-solved differences
	 * of the delivered quantity.
	 */
	for ( meq::ScalePin kind : { meq::ScalePin::ToroidalFlux, meq::ScalePin::PlasmaCurrent } )
	{
		meq::SafetyFactorEquilibrium::Options options;
		std::size_t const knots = meq::toroidalKnotCount( meq::plasmaCells( *machine.mesh ), order );
		meq::ToroidalKnotCut const cut = meq::toroidalKnotCut( knots );
		options.family.surfaces = cut.surfaces;
		options.family.innerCut = cut.innerCut;
		options.family.outerCut = cut.outerCut;
		options.pin = kind;
		meq::SafetyFactorEquilibrium pinned( *machine.solver, *machine.source, options );
		meq::SafetyFactorEvaluation here;
		BOOST_TEST_REQUIRE( pinned.evaluate( exact, here ) );
		meq::SafetyFactorDerivative derivative;
		BOOST_TEST_REQUIRE( pinned.differentiate( exact, derivative ) );
		auto measuredAt = [ & ]( std::size_t j, double h )
		{
			meq::ToroidalFieldSpline moved( exact );
			std::vector<double> values = exact.values();
			values[ j ] += h;
			moved.setValues( values );
			meq::SafetyFactorEvaluation out;
			BOOST_TEST_REQUIRE( pinned.evaluate( moved, out ) );
			return out.measured;
		};
		double rowWorst = 0.0, rowLargest = 0.0;
		for ( std::size_t j = 0; j < m; ++j )
		{
			double const h = 1.0e-6*exact.values()[ j ];
			double const d2 = ( measuredAt( j, 2.0*h ) - measuredAt( j, -2.0*h ) )/( 4.0*h );
			double const d1 = ( measuredAt( j, h ) - measuredAt( j, -h ) )/( 2.0*h );
			double const differenced = ( 4.0*d1 - d2 )/3.0;
			rowLargest = std::max( rowLargest, std::abs( differenced ) );
			rowWorst = std::max( rowWorst, std::abs( derivative.measured[ j ] - differenced ) );
		}
		char const *name = kind == meq::ScalePin::ToroidalFlux ? "ToroidalFlux" : "PlasmaCurrent";
		std::printf( "  the %s row: max | d_tangent - d_resolved | %.3e of %.3e, relative %.3e\n",
		             name, rowWorst, rowLargest, rowWorst/rowLargest );
		BOOST_TEST( rowWorst < 1.0e-6*rowLargest,
			"the " << name << " pin row's linearisation is " << rowWorst/rowLargest
			<< " of its size from re-solved differences of the delivered quantity" );
	}
}
