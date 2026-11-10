/*
 * THE q-DRIVEN SOLVE ON A FREE PLASMA EDGE: LIMITED AND DIVERTED.
 *
 * Where psi_bnd is an unknown -- pinned at a limiter point, or at the X-point
 * with ( R_X, z_X ) two more -- the retained linearisation carries them as
 * borders (GradShafranovSolver::linearise()) and the q-driven tangent moves
 * the boundary label with them (meq::SafetyFactorEquilibrium). Two machines
 * on the standard box, nine cells a side at k = 2, the q-driven source of
 * SafetyFactorSolver.cpp:
 *
 *   * LIMITED: psi_bnd = psi_h( 1.3, 0 ), a limiter point inside the box;
 *   * DIVERTED: the box's datum a divertor filament at ( 1.0, -0.68 ), below
 *     the box, less its harmonic Taylor part about the box's centre (1, R^2,
 *     z, R^2 z: each Delta*-harmonic), so the plasma keeps its size and the
 *     filament makes a saddle between it and the bottom wall; psi_bnd the
 *     saddle's, the X-point border bootstrapped by XP-2's pinned solves.
 *     Nine cells and not eight: on eight the saddle lies within 0.007 of the
 *     mesh line R = 1 and the X-point's Newton cycles between the two
 *     elements' roots.
 *
 * Both are UNCONFINED -- the source is smooth across Psi = 0 -- so the
 * discrete map is differentiable and the tangent can be held to re-solved
 * differences: a confined source's edge is a jump in the quadrature, which a
 * difference of solves resolves to the quadrature points, not to 1e-6.
 *
 * WHAT IT ASSERTS, per machine: the tangent Jacobian -- the q rows and the
 * current pin -- against Richardson differences of re-solved equilibria, to
 * 1e-6 of its size (SafetyFactorSolver's bound), where the reference is
 * differentiable; the flux pin refused; and the loop, started from a g 40%
 * off (20% diverted), recovering the g^2 its target was measured from to
 * 1e-7 in at most 10 iterations.
 */
#define BOOST_TEST_MODULE QDrivenFreeEdge
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include "mfem.hpp"
#include "meq/Coils.hpp"
#include "meq/ConductorField.hpp"
#include "meq/CriticalPoints.hpp"
#include "meq/GradShafranov.hpp"
#include "meq/Profiles.hpp"
#include "meq/SafetyFactor.hpp"
#include "meq/SafetyFactorSolve.hpp"
#include "meq/SafetyFactorTangent.hpp"
#include "meq/Source.hpp"
#include "convergence/ConvergenceHarness.hpp"

/*
 * UNDER MEQ_TEST_DEVICE, an mfem::Device configured before any case runs: the
 * q-driven tangent's device work -- the retained solve's funnels, the
 * tangent fields' vector operations, the pin rows' doubleDotKernel -- on the
 * device arms registered in tests/CMakeLists.txt.
 */
namespace
{
	struct DeviceFromEnvironment
	{
		DeviceFromEnvironment()
		{
			char const *const requested = std::getenv( "MEQ_TEST_DEVICE" );
			if ( requested != nullptr && *requested != '\0' )
				static mfem::Device device( requested );
		}
	};
}

BOOST_TEST_GLOBAL_FIXTURE( DeviceFromEnvironment );

namespace
{
	using meq::GradShafranovSolver;

	double const g0 = 2.20;
	double const g1 = 0.55;

	double referenceGSquared( double psi )
	{
		return ( g0 + g1*psi )*( g0 + g1*psi );
	}

	enum class Edge { Limited, Diverted };

	char const *nameOf( Edge edge )
	{
		return edge == Edge::Limited ? "limited" : "diverted";
	}

	/// The divertor filament's flux less its harmonic Taylor part at ( 1, 0 ):
	/// zero, with zero gradient and zero cross derivative, at the box's centre.
	struct Divertor
	{
		double radius = 1.0, z = -0.68, current = 1.5;
		double c0 = 0.0, c1 = 0.0, c2 = 0.0, c3 = 0.0;

		double filament( double r, double zz ) const
		{
			return meq::filamentPsi( r, zz, radius, z, 1.0, 1.0 );
		}

		Divertor()
		{
			double const e = 1.0e-4;
			double const value = filament( 1.0, 0.0 );
			double const dR = ( filament( 1.0 + e, 0.0 ) - filament( 1.0 - e, 0.0 ) )/( 2.0*e );
			double const dZ = ( filament( 1.0, e ) - filament( 1.0, -e ) )/( 2.0*e );
			double const dRZ = ( filament( 1.0 + e, e ) - filament( 1.0 + e, -e )
			                     - filament( 1.0 - e, e ) + filament( 1.0 - e, -e ) )/( 4.0*e*e );
			// T = c0 + c1 R^2 + c2 z + c3 R^2 z matched in value, d_R, d_z, d_Rz.
			c1 = dR/2.0;
			c0 = value - c1;
			c3 = dRZ/2.0;
			c2 = dZ - c3;
		}

		double operator()( double r, double zz ) const
		{
			return current*( filament( r, zz ) - ( c0 + c1*r*r + c2*zz + c3*r*r*zz ) );
		}
	};

	struct Machine
	{
		std::unique_ptr<mfem::Mesh> mesh;
		std::unique_ptr<meq::NormalisedMHDSource> source;
		std::unique_ptr<GradShafranovSolver> solver;
		std::unique_ptr<mfem::FunctionCoefficient> guess;
		std::unique_ptr<mfem::FunctionCoefficient> datum;
		std::unique_ptr<meq::SafetyFactorEquilibrium> equilibrium;
		Edge edge;
		int order;

		std::unique_ptr<meq::ConductorField> conductors;

		/// @a withConductor: a filament below the box, taken out of the
		/// mesh by the conductor split, so that the extraction reads
		/// psi_p + psi_c and q_p + q_c.
		Machine( Edge edgeIn, int cells, int orderIn, meq::ScalePin pin = meq::ScalePin::EdgeField,
		         bool withConductor = false )
			: edge( edgeIn ), order( orderIn )
		{
			meq::tests::Rectangle const box = meq::tests::standardBox();
			mesh = std::make_unique<mfem::Mesh>( meq::tests::makeMesh( box, cells ) );
			auto pPrime = std::make_shared<meq::ConstantProfile const>( 0.45 );
			auto flat = std::make_shared<meq::ConstantProfile const>( 0.0 );
			source = std::make_unique<meq::NormalisedMHDSource>( pPrime, flat, 0.30, 1.0 );
			Divertor const divertor;
			bool const diverted = edge == Edge::Diverted;
			auto vacuum = [ divertor, diverted ]( double r, double z )
			{
				return diverted ? divertor( r, z ) : 0.0;
			};
			guess = std::make_unique<mfem::FunctionCoefficient>( [ box, vacuum ]( mfem::Vector const &x )
			{
				return 0.30*std::sin( M_PI*( x( 0 ) - box.minRadius )/box.width() )
				       *std::sin( M_PI*( x( 1 ) - box.zMin )/box.height() ) + vacuum( x( 0 ), x( 1 ) );
			} );
			datum = std::make_unique<mfem::FunctionCoefficient>( [ vacuum ]( mfem::Vector const &x )
			{
				return vacuum( x( 0 ), x( 1 ) );
			} );

			auto configure = [ & ]( GradShafranovSolver &s )
			{
				s.setSource( *source, 0.30 );
				s.setBoundaryData( *datum );
				s.setInitialGuess( *guess );
				s.setNewtonControl( 1.0e-12, 1.0e-14, 60 );
				meq::tests::namedRoutesFromEnvironment( s );
			};
			solver = std::make_unique<GradShafranovSolver>( *mesh, order );
			if ( withConductor )
			{
				conductors = std::make_unique<meq::ConductorField>( 1.0 );
				conductors->add( meq::CurrentFilament( 1.0, -0.9, 0.2 ) );
				solver->setConductorField( *conductors );
			}
			double const seedR = 1.014, seedZ = -0.504;
			if ( edge == Edge::Limited )
				solver->setBoundaryFluxPoint( 1.3, 0.0 );
			configure( *solver );
			std::size_t const knots = meq::toroidalKnotCount( meq::plasmaCells( *mesh ), order );
			meq::ToroidalKnotCut const cut = meq::toroidalKnotCut( knots );
			meq::SafetyFactorEquilibrium::Options options;
			options.family.surfaces = cut.surfaces;
			options.family.innerCut = cut.innerCut;
			options.family.outerCut = cut.outerCut;
			options.pin = pin;
			equilibrium = std::make_unique<meq::SafetyFactorEquilibrium>( *solver, *source, options );
			if ( diverted )
			{
				/*
				 * XP-2'S START: psi_bnd pinned at the seed, solved, and the
				 * X-point border continued from there, at the reference g -- from the guess itself
				 * the bordered Newton's first corner is singular (TODO, *The
				 * X-point border from a cold start*). The loop's first
				 * evaluate() consumes the state.
				 */
				source->setGGPrime( splineOf( referenceGSquared ).ggPrime() );
				double radius = seedR, z = seedZ;
				std::unique_ptr<GradShafranovSolver> pinned;
				for ( int sweep = 0; sweep < 12; ++sweep )
				{
					double const previousR = radius, previousZ = z;
					// XP-2: pin at the estimate, solve, re-locate the saddle.
					pinned = std::make_unique<GradShafranovSolver>( *mesh, order );
					pinned->setBoundaryFluxPoint( radius, z );
					configure( *pinned );
					pinned->solve();
					meq::CriticalPoint saddle;
					meq::CriticalPointFinder const finder( *pinned );
					BOOST_TEST_REQUIRE( finder.tryFindCriticalPointFrom( radius, z, meq::AxisSense::Saddle, saddle ),
					                    "no saddle near ( " << radius << ", " << z << " )" );
					radius = saddle.radius;
					z = saddle.z;
					std::printf( "    XP-2 sweep %d: psi_ax %.6f, psi_bnd %.6f, saddle ( %.6f, %.6f ) psi %.6f\n",
					             sweep, pinned->psiAxis(), pinned->psiBoundary(), radius, z, saddle.psi );
					if ( std::hypot( radius - previousR, z - previousZ ) < 1.0e-6 )
						break;
				}
				solver->setXPointBoundary( radius, z );
				solver->setInitialState( pinned->state(), {} );
				solver->setBoundaryFluxInitial( pinned->psiBoundary() );
				solver->setAxisFluxInitial( pinned->psiAxis() );
			}

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

		void describe() const
		{
			std::printf( "    %s: psi_ax %.6f, psi_bnd %.6f", nameOf( edge ), solver->psiAxis(),
			             solver->psiBoundary() );
			if ( edge == Edge::Diverted )
				std::printf( ", X-point ( %.4f, %.4f )", solver->xPointR(), solver->xPointZ() );
			std::printf( ", axis ( %.4f, %.4f )\n", solver->axisR(), solver->axisZ() );
		}
	};

	/*
	 * THE TANGENT JACOBIAN AGAINST RE-SOLVED DIFFERENCES, the q rows and both
	 * quadrature pins, at the closed form's g^2.
	 */
	void tangentAgainstSolves( Edge edge )
	{
		int const cells = 9, order = 2;
		Machine machine( edge, cells, order );
		meq::ToroidalFieldSpline const exact = machine.splineOf( referenceGSquared );
		meq::SafetyFactorEvaluation base;
		bool const evaluated = machine.equilibrium->evaluate( exact, base );
		if ( !evaluated )
		{
			std::printf( "    residuals:" );
			for ( double r : machine.solver->newtonResiduals() )
				std::printf( " %.2e", r );
			std::printf( "\n    X ( %.5f, %.5f )\n", machine.solver->xPointR(), machine.solver->xPointZ() );
		}
		BOOST_TEST_REQUIRE( evaluated, machine.equilibrium->refusal() );
		std::printf( "\n  THE TANGENT JACOBIAN, %s, %zu knots\n", nameOf( edge ), exact.size() );
		machine.describe();
		meq::SafetyFactorDerivative tangent;
		BOOST_TEST_REQUIRE( machine.equilibrium->differentiate( exact, tangent ) );
		std::size_t const m = exact.size();
		double boundaryMoves = 0.0;
		for ( double d : machine.equilibrium->boundaryDerivative() )
			boundaryMoves = std::max( boundaryMoves, std::abs( d ) );
		std::printf( "    max | d psi_bnd/d u_j | %.3e\n", boundaryMoves );
		BOOST_TEST( boundaryMoves > 0.0, "psi_bnd is free and its tangent is zero" );

		auto qOverGAt = [ & ]( std::size_t j, double h )
		{
			meq::ToroidalFieldSpline moved( exact );
			std::vector<double> values = exact.values();
			values[ j ] += h;
			moved.setValues( values );
			meq::SafetyFactorEvaluation out;
			BOOST_TEST_REQUIRE( machine.equilibrium->evaluate( moved, out ),
			                    machine.equilibrium->refusal() );
			return out.qOverG;
		};
		/*
		 * THE REFERENCE MUST BE DIFFERENTIABLE WHERE IT JUDGES, and that is a
		 * precondition checked, not assumed: the re-solved map re-locates the
		 * axis on every field, and where the axis sits within a few 1e-4 of a
		 * face the root can change element inside 2e-5 of the state -- a jump
		 * of O( h^( k+1 ) ) in the reference, which the tangent, moving the
		 * origin by the root's derivative, does not have (M-391, M-401).
		 * So each entry is differenced twice, Richardson at h and at h/4, and
		 * judged only where the two agree to 1e-7 of the Jacobian; at least
		 * 90% must.
		 */
		std::vector<std::vector<double>> richardson( 2, std::vector<double>( m*m, 0.0 ) );
		for ( std::size_t j = 0; j < m; ++j )
			for ( int pass = 0; pass < 2; ++pass )
			{
				double const h = ( pass == 0 ? 1.0e-6 : 0.25e-6 )*exact.values()[ j ];
				std::vector<double> const p2 = qOverGAt( j, 2.0*h ), m2 = qOverGAt( j, -2.0*h );
				std::vector<double> const p1 = qOverGAt( j, h ), m1 = qOverGAt( j, -h );
				for ( std::size_t i = 0; i < m; ++i )
				{
					double const d2 = ( p2[ i ] - m2[ i ] )/( 4.0*h );
					double const d1 = ( p1[ i ] - m1[ i ] )/( 2.0*h );
					richardson[ static_cast<std::size_t>( pass ) ][ i*m + j ] = ( 4.0*d1 - d2 )/3.0;
				}
			}
		double largest = 0.0;
		for ( double d : richardson[ 0 ] )
			largest = std::max( largest, std::abs( d ) );
		double worst = 0.0, worstSkipped = 0.0;
		std::size_t judged = 0;
		for ( std::size_t k = 0; k < m*m; ++k )
		{
			double const apart = std::abs( tangent.qOverG[ k ] - richardson[ 0 ][ k ] );
			if ( std::abs( richardson[ 0 ][ k ] - richardson[ 1 ][ k ] ) <= 1.0e-7*largest )
			{
				++judged;
				worst = std::max( worst, apart );
			}
			else
				worstSkipped = std::max( worstSkipped, apart );
		}
		double const fraction = static_cast<double>( judged )/static_cast<double>( m*m );
		std::printf( "    q rows: max | J_tangent - J_resolved | %.3e of max | J | %.3e, relative %.3e, over "
		             "%zu of %zu entries; the rest, where the reference is not differentiable, %.3e\n",
		             worst, largest, worst/largest, judged, m*m, worstSkipped/largest );
		bool const resolvable = std::numeric_limits<mfem::real_t>::epsilon() < 1.0e-10;
		if ( !resolvable )
		{
			std::printf( "    (single: the re-solved reference cannot resolve the derivative; not asserted)\n" );
			return;
		}
		BOOST_TEST( fraction >= 0.9, nameOf( edge ) << ": the re-solved reference is differentiable at only "
		            << judged << " of " << m*m << " entries -- a fixture too near an element switch to judge by" );
		BOOST_TEST( worst < 1.0e-6*largest,
			nameOf( edge ) << ": the tangent Jacobian is " << worst/largest << " of its size from "
			"re-solved differences -- a border the tangent does not carry, or a label it does not move" );

		/*
		 * THE FLUX PIN IS REFUSED BESIDE A FREE EDGE: g( 0 ) = R0B0 is not
		 * zero, so the delivered Phi = int_{Psi > 0} g/R jumps as each
		 * quadrature point crosses the moving edge, and is not differentiable
		 * in the state (TODO, *The flux pin beside a free plasma edge*).
		 */
		{
			meq::SafetyFactorEquilibrium::Options options;
			options.pin = meq::ScalePin::ToroidalFlux;
			BOOST_CHECK_THROW( meq::SafetyFactorEquilibrium( *machine.solver, *machine.source, options ),
			                   std::invalid_argument );
		}

		for ( meq::ScalePin kind : { meq::ScalePin::PlasmaCurrent } )
		{
			meq::SafetyFactorEquilibrium::Options options;
			std::size_t const knots = meq::toroidalKnotCount( meq::plasmaCells( *machine.mesh ), order );
			meq::ToroidalKnotCut const cut = meq::toroidalKnotCut( knots );
			options.family.surfaces = cut.surfaces;
			options.family.innerCut = cut.innerCut;
			options.family.outerCut = cut.outerCut;
			options.pin = kind;
			meq::SafetyFactorEquilibrium pinned( *machine.solver, *machine.source, options );
			pinned.warmFromSolver();
			meq::SafetyFactorEvaluation here;
			BOOST_TEST_REQUIRE( pinned.evaluate( exact, here ), pinned.refusal() );
			meq::SafetyFactorDerivative derivative;
			BOOST_TEST_REQUIRE( pinned.differentiate( exact, derivative ) );
			auto measuredAt = [ & ]( std::size_t j, double h )
			{
				meq::ToroidalFieldSpline moved( exact );
				std::vector<double> values = exact.values();
				values[ j ] += h;
				moved.setValues( values );
				meq::SafetyFactorEvaluation out;
				BOOST_TEST_REQUIRE( pinned.evaluate( moved, out ), pinned.refusal() );
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
			char const *name = "PlasmaCurrent";
			std::printf( "    the %s row: max | d_tangent - d_resolved | %.3e of %.3e, relative %.3e\n",
			             name, rowWorst, rowLargest, rowWorst/rowLargest );
			BOOST_TEST( rowWorst < 1.0e-6*rowLargest,
				nameOf( edge ) << ": the " << name << " pin row is " << rowWorst/rowLargest
				<< " of its size from re-solved differences of the delivered quantity" );
		}
		std::fflush( stdout );
	}

	/*
	 * THE LOOP RECOVERS THE FIELD ITS TARGET WAS BUILT FROM: the target q at
	 * the knots of the closed form's own equilibrium, the start g 40% large
	 * (20% diverted, below), pinned by R0B0.
	 */
	void roundTrip( Edge edge )
	{
		int const cells = 9, order = 2;
		Machine machine( edge, cells, order );
		meq::ToroidalFieldSpline const exact = machine.splineOf( referenceGSquared );
		meq::SafetyFactorEvaluation evaluation;
		BOOST_TEST_REQUIRE( machine.equilibrium->evaluate( exact, evaluation ),
		                    machine.equilibrium->refusal() );
		std::vector<double> target;
		for ( std::size_t i = 0; i < evaluation.qOverG.size(); ++i )
			target.push_back( evaluation.qOverG[ i ]*std::sqrt( exact.values()[ i ] ) );
		std::vector<double> const labels = machine.equilibrium->labels();
		auto targetAt = [ labels, target ]( double psi )
		{
			std::size_t best = 0;
			for ( std::size_t i = 1; i < labels.size(); ++i )
				if ( std::abs( 1.0 - labels[ i ] - psi ) < std::abs( 1.0 - labels[ best ] - psi ) )
					best = i;
			return target[ best ];
		};
		meq::SafetyFactorSystem::Pin pin;
		pin.kind = meq::ScalePin::EdgeField;
		pin.value = g0;
		meq::SafetyFactorEquilibrium *equilibrium = machine.equilibrium.get();
		meq::SafetyFactorSystem system(
			labels, targetAt, pin,
			[ equilibrium ]( meq::ToroidalFieldSpline const &s, meq::SafetyFactorEvaluation &out )
			{
				return equilibrium->evaluate( s, out );
			},
			[ equilibrium ]( meq::ToroidalFieldSpline const &s, meq::SafetyFactorDerivative &out )
			{
				return equilibrium->differentiate( s, out );
			},
			meq::newtonTolerance( 1.0e-12 ) );

		/*
		 * g 40% LARGE ON THE LIMITED BOX, 20% ON THE DIVERTED ONE -- A
		 * PRECONDITION, NOT A TOLERANCE. From 40% the diverted start's
		 * equilibrium is bounded by another saddle: the first three steps are
		 * refused as the X-point border leaves the mesh, and the loop settles
		 * on that branch with a misfit of 5e-5, which the exact target cannot
		 * have. From 20% and 10% it converges quadratically (M-401).
		 */
		double const startFactor = edge == Edge::Diverted ? 1.44 : 1.96;
		std::vector<double> start;
		for ( double label : labels )
			start.push_back( startFactor*referenceGSquared( 1.0 - label ) );
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

		std::printf( "\n  THE LOOP AGAINST A KNOWN g, %s: %zu knots, %s\n", nameOf( edge ), start.size(),
		             result.status.c_str() );
		machine.describe();
		std::printf( "    %d iterations, %d equilibria (%d refused), misfit %.3e, conditioning %.3e; "
		             "removable by iteration:", result.iterations, system.solves(), system.refusals(),
		             result.misfit, result.conditioning );
		for ( double h : result.history )
			std::printf( " %.1e", h );
		std::printf( "\n" );
		if ( system.refusals() > 0 )
			std::printf( "    the last refusal: %s\n", machine.equilibrium->refusal().c_str() );
		BOOST_TEST_REQUIRE( result.converged, nameOf( edge ) << ": " << result.status
			<< " (a Jacobian that misleads the step has its trials refused: the differenced one did in "
			"single, MEASUREMENTS.md M-406)" );

		double worst = 0.0;
		for ( std::size_t i = 0; i < result.x.size(); ++i )
		{
			double const exactValue = referenceGSquared( 1.0 - labels[ i ] );
			worst = std::max( worst, std::abs( result.x[ i ] - exactValue )/exactValue );
		}
		std::printf( "    worst |g^2/g^2_exact - 1| at the knots %.3e\n", worst );
		std::fflush( stdout );
		bool const resolvable = std::numeric_limits<mfem::real_t>::epsilon() < 1.0e-10;
		if ( resolvable )
			BOOST_TEST( worst < 1.0e-7, nameOf( edge ) << ": the loop landed " << worst
			            << " from the closed form its target was measured from" );
		else
			std::printf( "    (single: printed, not asserted -- TODO, *The q-driven solve in single "
			             "precision*)\n" );
		BOOST_TEST( result.iterations <= 10,
			nameOf( edge ) << ": the loop took " << result.iterations << " iterations: a crawl "
			"that converges is a defect report" );
	}
}

namespace
{
	/*
	 * THE SHAPE DERIVATIVE ON A FREE EDGE (M-408): the extraction's
	 * analytic derivative against central differences of the extraction, the
	 * labels psi_ax and psi_bnd moving with each tangent, at a ladder of
	 * steps. In double the best step must come within 1e-8 of the Jacobian's
	 * size; in single the differences are noise (M-406), printed only.
	 */
	void shapeAgainstDifferences( Edge edge, bool withConductor = false )
	{
		Machine machine( edge, 9, 2, meq::ScalePin::EdgeField, withConductor );
		meq::ToroidalFieldSpline const exact = machine.splineOf( referenceGSquared );
		meq::SafetyFactorEvaluation base;
		BOOST_TEST_REQUIRE( machine.equilibrium->evaluate( exact, base ), machine.equilibrium->refusal() );
		meq::SafetyFactorDerivative analytic;
		BOOST_TEST_REQUIRE( machine.equilibrium->differentiate( exact, analytic ) );
		double largest = 0.0;
		for ( double v : analytic.qOverG )
			largest = std::max( largest, std::abs( v ) );
		std::printf( "\n  %s%s: THE SHAPE DERIVATIVE AGAINST THE DIFFERENCED EXTRACTION, %d nodes, max | J | %.6e\n",
		             nameOf( edge ), withConductor ? ", a conductor out of the mesh" : "",
		             machine.equilibrium->shapeNodes(), largest );
		double best = std::numeric_limits<double>::infinity();
		for ( double step : { 1.0e-3, 1.0e-4, 1.0e-5, 1.0e-6, 1.0e-7 } )
		{
			meq::SafetyFactorEquilibrium::Options options;
			std::size_t const knots = meq::toroidalKnotCount( meq::plasmaCells( *machine.mesh ), 2 );
			meq::ToroidalKnotCut const cut = meq::toroidalKnotCut( knots );
			options.family.surfaces = cut.surfaces;
			options.family.innerCut = cut.innerCut;
			options.family.outerCut = cut.outerCut;
			options.differenced = true;
			options.step = step;
			meq::SafetyFactorEquilibrium differenced( *machine.solver, *machine.source, options );
			differenced.warmFromSolver();
			meq::SafetyFactorEvaluation again;
			BOOST_TEST_REQUIRE( differenced.evaluate( exact, again ), differenced.refusal() );
			meq::SafetyFactorDerivative reference;
			if ( !differenced.differentiate( exact, reference ) )
			{
				std::printf( "    step %.0e: a perturbed family could not be extracted\n", step );
				continue;
			}
			double apart = 0.0;
			for ( std::size_t k = 0; k < reference.qOverG.size(); ++k )
				apart = std::max( apart, std::abs( analytic.qOverG[ k ] - reference.qOverG[ k ] ) );
			apart /= largest;
			best = std::min( best, apart );
			std::printf( "    step %.0e: max | J_shape - J_differenced |/max | J | %.3e\n", step, apart );
		}
		std::printf( "    at the best step %.3e\n", best );
		if ( std::numeric_limits<mfem::real_t>::epsilon() > 1.0e-10 )
		{
			std::printf( "  (single: the differences are noise-dominated, M-406; not asserted)\n" );
			return;
		}
		BOOST_TEST( best < 1.0e-8,
			nameOf( edge ) << ": the shape derivative is " << best << " of the Jacobian's size from the "
			"differenced extraction at its best step" );
	}
}

BOOST_AUTO_TEST_CASE( theLimitedShapeDerivativeIsTheDifferencedOne ) { shapeAgainstDifferences( Edge::Limited ); }
BOOST_AUTO_TEST_CASE( theDivertedShapeDerivativeIsTheDifferencedOne ) { shapeAgainstDifferences( Edge::Diverted ); }
BOOST_AUTO_TEST_CASE( theShapeDerivativeCarriesTheConductors ) { shapeAgainstDifferences( Edge::Limited, true ); }
BOOST_AUTO_TEST_CASE( theLimitedTangentIsTheDifferencedOne ) { tangentAgainstSolves( Edge::Limited ); }
BOOST_AUTO_TEST_CASE( theDivertedTangentIsTheDifferencedOne ) { tangentAgainstSolves( Edge::Diverted ); }
BOOST_AUTO_TEST_CASE( theLimitedLoopRecoversItsField ) { roundTrip( Edge::Limited ); }
BOOST_AUTO_TEST_CASE( theDivertedLoopRecoversItsField ) { roundTrip( Edge::Diverted ); }

/*
 * THE EVALUATION'S REPEATABILITY, A MEASUREMENT (disabled; run by name): q at
 * the closed form, then at g^2 1.96 times it, then at the closed form again from
 * that warm start -- how far the same spline's q moves with the path to it,
 * which is the floor the loop's rows carry.
 */
BOOST_AUTO_TEST_CASE( theEvaluationIsRepeatable, * boost::unit_test::disabled() )
{
	for ( Edge edge : { Edge::Limited, Edge::Diverted } )
	{
		Machine machine( edge, 9, 2 );
		meq::ToroidalFieldSpline const exact = machine.splineOf( referenceGSquared );
		double const factor = 1.96;
		meq::ToroidalFieldSpline const moved = machine.splineOf(
			[ factor ]( double psi ) { return factor*referenceGSquared( psi ); } );
		meq::SafetyFactorEvaluation first, away, again;
		BOOST_TEST_REQUIRE( machine.equilibrium->evaluate( exact, first ), machine.equilibrium->refusal() );
		double const axis0 = machine.solver->psiAxis(), boundary0 = machine.solver->psiBoundary();
		mfem::Vector const state0( machine.solver->state() );
		std::vector<double> const residuals0 = machine.solver->newtonResiduals();
		BOOST_TEST_REQUIRE( machine.equilibrium->evaluate( moved, away ), machine.equilibrium->refusal() );
		BOOST_TEST_REQUIRE( machine.equilibrium->evaluate( exact, again ), machine.equilibrium->refusal() );
		double worst = 0.0;
		for ( std::size_t i = 0; i < first.qOverG.size(); ++i )
			worst = std::max( worst, std::abs( again.qOverG[ i ]/first.qOverG[ i ] - 1.0 ) );
		mfem::Vector d( machine.solver->state() );
		d -= state0;
		std::printf( "  %s: q/g moves %.3e with the path; psi_ax %.3e, psi_bnd %.3e, state %.3e of %.3e\n",
		             nameOf( edge ), worst, machine.solver->psiAxis() - axis0,
		             machine.solver->psiBoundary() - boundary0, d.Normlinf(), state0.Normlinf() );
		std::printf( "    first solve's residuals:" );
		for ( double r : residuals0 )
			std::printf( " %.2e", r );
		std::printf( "\n    again's:" );
		for ( double r : machine.solver->newtonResiduals() )
			std::printf( " %.2e", r );
		std::printf( "\n" );
	}
}

/*
 * THE TANGENT JACOBIAN WRITTEN OUT, A MEASUREMENT (disabled; run by name in
 * each precision): the limited and diverted boxes at the closed form,
 * jacobian-<edge>-<precision>.txt, to hold a single build's to a double's.
 */
BOOST_AUTO_TEST_CASE( theTangentJacobianWrittenOut, * boost::unit_test::disabled() )
{
	bool const single = std::numeric_limits<mfem::real_t>::epsilon() > 1.0e-10;
	for ( Edge edge : { Edge::Limited, Edge::Diverted } )
	{
		Machine machine( edge, 9, 2 );
		meq::ToroidalFieldSpline const exact = machine.splineOf( referenceGSquared );
		meq::SafetyFactorEvaluation base;
		BOOST_TEST_REQUIRE( machine.equilibrium->evaluate( exact, base ), machine.equilibrium->refusal() );
		meq::SafetyFactorDerivative tangent;
		BOOST_TEST_REQUIRE( machine.equilibrium->differentiate( exact, tangent ) );
		std::string const name = std::string( "jacobian-" ) + nameOf( edge ) + ( single ? "-single" : "-double" ) + ".txt";
		std::FILE *out = std::fopen( name.c_str(), "w" );
		for ( double v : tangent.qOverG )
			std::fprintf( out, "%.17g\n", v );
		for ( double v : base.qOverG )
			std::fprintf( out, "# %.17g\n", v );
		for ( double v : machine.equilibrium->axisDerivative() )
			std::fprintf( out, "a %.17g\n", v );
		for ( double v : machine.equilibrium->boundaryDerivative() )
			std::fprintf( out, "b %.17g\n", v );
		std::fclose( out );
		std::printf( "  wrote %s\n", name.c_str() );
	}
}
