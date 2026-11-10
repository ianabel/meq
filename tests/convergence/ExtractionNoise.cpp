#define BOOST_TEST_MODULE ExtractionNoise
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "mfem.hpp"

#include "meq/CriticalPoints.hpp"
#include "meq/Profiles.hpp"
#include "meq/Source.hpp"

#include "ExtractionCrossing.hpp"
#include "RecoveryInstrument.hpp"

/*
 * WHAT CARRIES THE FAMILY'S NOISE, LINK BY LINK -- the instrument behind
 * MEASUREMENTS.md's M-371. BUILT, NOT REGISTERED: it asserts nothing.
 *
 * ExtractionCrossing.hpp's linear box, A = 0.30 + j * 1e-5 for j = 0 ... 20,
 * each state its own solve. Over 21 states 1e-5 apart a smooth quantity's
 * second difference is ~1e-10 of it, so the second differences are the noise:
 * for independent draws of size sigma their root mean square is sqrt( 6 )
 * sigma, which is what is printed as sigma. Every link of the chain is read:
 * psi_h, psi*, psi~ (dofs, of their largest value), the located axis (its
 * position, its psi of the span, its element), psi~ and psi_h at the located
 * axis and at a fixed point, and the family -- six averages per surface, of
 * their value -- under three axes:
 *
 *   located  the finder's axis as it comes: the point, psi_h there;
 *   psi~     the finder's point, psi~ there;
 *   fixed    state 10's point for every state, psi~ there -- no axis noise
 *            but what psi~ carries at a fixed point.
 *
 * and on psi* and q* with the located axis, for comparison. "psi~" is what
 * meq::labelledAxis() gives a solve with no psi_ax unknown, as this one is.
 *
 * The normalised case is the same reading on a bordered solve, whose psi_ax
 * is an unknown and labels the family: SafetyFactorSolver's source, gg'
 * linear in Psi with its slope scaled by 1 + j * 1e-5, each state a cold
 * solve. It reads the solve's psi_ax beside the finder's psi_h and psi~ there.
 *
 * MEQ_TEST_DEVICE, MEQ_TEST_BATCHED and MEQ_LOCAL_REFINEMENT choose the
 * device, the batched routes and the local refinement
 * (RecoveryInstrument.hpp); each case prints the routes it took.
 */
using meq::tests::recovery::DeviceFromEnvironment;
BOOST_TEST_GLOBAL_FIXTURE( DeviceFromEnvironment );
namespace
{
	double const centre = 0.30;
	double const spacing = 1.0e-5;
	int const states = 21;

	double valueAt( mfem::GridFunction const &field, mfem::Mesh &mesh, double r, double z )
	{
		mfem::DenseMatrix point( 2, 1 );
		point( 0, 0 ) = r;
		point( 1, 0 ) = z;
		mfem::Array<int> elements;
		mfem::Array<mfem::IntegrationPoint> ips;
		mesh.FindPoints( point, elements, ips, false );
		if ( elements[ 0 ] < 0 )
			return std::nan( "" );
		return field.GetValue( elements[ 0 ], ips[ 0 ] );
	}

	std::vector<double> averages( meq::FluxSurfaceFamily const &family )
	{
		std::vector<double> out;
		double const span = family.psiBoundary - family.psiAxis;
		for ( meq::FluxSurface const &s : family.surfaces )
		{
			double const rho = std::sqrt( s.normalisedFlux );
			double const dRhoDPsi = 1.0/( 2.0*rho*span );
			out.insert( out.end(), { s.vPrime, s.gradPsiSquared, s.inverseRSquared,
			                         s.gradPsiSquaredOverRSquared, s.vPrime/std::abs( dRhoDPsi ),
			                         s.gradPsiSquared*dRhoDPsi*dRhoDPsi } );
		}
		return out;
	}

	char const *averageName( std::size_t q )
	{
		static char const *const names[] = { "V'", "<|grad psi|^2>", "<R^-2>", "<|grad psi|^2/R^2>",
		                                     "V'_rho", "<|grad rho|^2>" };
		return names[ q ];
	}

	/// sigma of each component over the states, from second differences,
	/// relative to @a scale ( component-wise when empty ).
	std::vector<double> sigma( std::vector<std::vector<double>> const &series, std::vector<double> scale )
	{
		std::size_t const n = series.front().size();
		if ( scale.empty() )
			for ( std::size_t i = 0; i < n; ++i )
				scale.push_back( std::abs( series[ series.size()/2 ][ i ] ) );
		std::vector<double> out( n, 0.0 );
		for ( std::size_t i = 0; i < n; ++i )
		{
			double sum = 0.0;
			for ( std::size_t j = 1; j + 1 < series.size(); ++j )
			{
				double const d2 = series[ j + 1 ][ i ] - 2.0*series[ j ][ i ] + series[ j - 1 ][ i ];
				sum += d2*d2;
			}
			out[ i ] = std::sqrt( sum/static_cast<double>( series.size() - 2 ) )/std::sqrt( 6.0 )
			           /std::max( scale.size() == 1 ? scale[ 0 ] : scale[ i ], 1.0e-300 );
		}
		return out;
	}

	double worst( std::vector<double> const &v )
	{
		double w = 0.0;
		for ( double x : v )
			w = std::max( w, x );
		return w;
	}

	std::vector<double> dofs( mfem::GridFunction const &field )
	{
		mfem::Vector copy( field );
		copy.HostRead();
		return std::vector<double>( copy.begin(), copy.end() );
	}

	double largest( std::vector<double> const &v )
	{
		double m = 0.0;
		for ( double x : v )
			m = std::max( m, std::abs( x ) );
		return m;
	}

	void reportFamily( char const *name, std::vector<std::vector<double>> const &series )
	{
		std::vector<double> const s = sigma( series, {} );
		std::size_t const surfaces = s.size()/6;
		std::printf( "    %-24s worst %.2e;", name, worst( s ) );
		for ( std::size_t q = 0; q < 6; ++q )
		{
			double w = 0.0;
			for ( std::size_t i = 0; i < surfaces; ++i )
				w = std::max( w, s[ 6*i + q ] );
			std::printf( " %s %.1e", averageName( q ), w );
		}
		std::printf( "\n      by surface (worst average):" );
		for ( std::size_t i = 0; i < surfaces; ++i )
		{
			double w = 0.0;
			for ( std::size_t q = 0; q < 6; ++q )
				w = std::max( w, s[ 6*i + q ] );
			std::printf( " %.1e", w );
		}
		std::printf( "\n" );
	}

	void measure( int n, int order )
	{
		std::vector<std::vector<double>> psiH, psiStar, psiTilde, trace, fluxH, psiHAtFixed;
		std::vector<std::vector<double>> axisPoint, axisPsi, nodes;
		std::vector<std::vector<double>> located, tildeAxis, fixedAxis, broken;
		std::vector<int> axisElements;
		std::vector<double> overshoots;
		double fixedR = 0.0, fixedZ = 0.0, span = 0.0;

		// The fixed point first: state 10's located axis.
		{
			meq::tests::crossing::Solved s( centre + 10*spacing, n, order, meq::tests::recovery::configure );
			meq::tests::recovery::reportRoutes( s.theSolver() );
			meq::CriticalPoint const axis = meq::CriticalPointFinder( s.theSolver() ).findAxis();
			fixedR = axis.radius;
			fixedZ = axis.z;
		}

		meq::FluxFamilyOptions options;
		options.surfaces = 12;
		options.angles = 96;
		options.innerCut = 0.05;
		options.outerCut = 0.90;

		for ( int j = 0; j < states; ++j )
		{
			meq::tests::crossing::Solved s( centre + j*spacing, n, order, meq::tests::recovery::configure );
			meq::GradShafranovSolver const &solver = s.theSolver();
			psiH.push_back( dofs( solver.potential() ) );
			psiStar.push_back( dofs( solver.postProcessedPotential() ) );
			psiTilde.push_back( dofs( solver.continuousPotential() ) );
			trace.push_back( dofs( solver.trace() ) );
			fluxH.push_back( dofs( solver.flux() ) );

			meq::CriticalPoint const axis = meq::CriticalPointFinder( solver ).findAxis();
			double const tildeAtAxis = valueAt( solver.continuousPotential(), s.theMesh(), axis.radius, axis.z );
			double const tildeAtFixed = valueAt( solver.continuousPotential(), s.theMesh(), fixedR, fixedZ );
			span = axis.psi - solver.psiBoundary();
			axisPoint.push_back( { axis.radius, axis.z } );
			axisPsi.push_back( { axis.psi, tildeAtAxis, tildeAtFixed,
			                     valueAt( solver.potential(), s.theMesh(), fixedR, fixedZ ),
			                     valueAt( solver.postProcessedPotential(), s.theMesh(), fixedR, fixedZ ) } );
			axisElements.push_back( axis.element );
			overshoots.push_back( axis.overshoot );

			meq::ContourTracer const continuous( solver, meq::Potential::Continuous );
			meq::ContourTracer const postProcessed( solver, meq::Potential::PostProcessed );
			meq::FluxSurfaceFamily const a = meq::extractFluxSurfaces( continuous, axis, solver.psiBoundary(), options );
			located.push_back( averages( a ) );
			std::vector<double> positions;
			for ( meq::FluxSurface const &surface : a.surfaces )
				for ( std::size_t node = 0; node < surface.count(); ++node )
				{
					positions.push_back( surface.radius[ node ] );
					positions.push_back( surface.z[ node ] );
				}
			nodes.push_back( positions );

			meq::CriticalPoint tilde = axis;
			tilde.psi = tildeAtAxis;
			tildeAxis.push_back( averages( meq::extractFluxSurfaces( continuous, tilde, solver.psiBoundary(), options ) ) );

			meq::CriticalPoint fixed = axis;
			fixed.radius = fixedR;
			fixed.z = fixedZ;
			fixed.psi = tildeAtFixed;
			fixedAxis.push_back( averages( meq::extractFluxSurfaces( continuous, fixed, solver.psiBoundary(), options ) ) );

			broken.push_back( averages( meq::extractFluxSurfaces( postProcessed, axis, solver.psiBoundary(), options ) ) );
		}

		std::size_t elementChanges = 0;
		for ( std::size_t j = 1; j < axisElements.size(); ++j )
			elementChanges += axisElements[ j ] != axisElements[ j - 1 ];
		std::vector<int> distinct( axisElements );
		std::sort( distinct.begin(), distinct.end() );
		distinct.erase( std::unique( distinct.begin(), distinct.end() ), distinct.end() );

		std::printf( "\n  box, n = %d, k = %d, %s: sigma over %d states %.0e apart in A\n", n, order,
		             sizeof( mfem::real_t ) == 4 ? "single" : "double", states, spacing );
		std::printf( "    psi_h %.2e, psi* %.2e, psi~ %.2e of their largest value (worst dof)\n",
		             worst( sigma( psiH, { largest( psiH[ 10 ] ) } ) ),
		             worst( sigma( psiStar, { largest( psiStar[ 10 ] ) } ) ),
		             worst( sigma( psiTilde, { largest( psiTilde[ 10 ] ) } ) ) );
		auto rms = []( std::vector<double> const &v )
		{
			double sum = 0.0;
			for ( double x : v )
				sum += x*x;
			return std::sqrt( sum/static_cast<double>( v.size() ) );
		};
		std::printf( "    root mean square over dofs: psi_h %.2e, q_h %.2e, trace %.2e, psi* %.2e, psi~ %.2e; "
		             "worst q_h %.2e, trace %.2e\n",
		             rms( sigma( psiH, { largest( psiH[ 10 ] ) } ) ), rms( sigma( fluxH, { largest( fluxH[ 10 ] ) } ) ),
		             rms( sigma( trace, { largest( trace[ 10 ] ) } ) ),
		             rms( sigma( psiStar, { largest( psiStar[ 10 ] ) } ) ),
		             rms( sigma( psiTilde, { largest( psiTilde[ 10 ] ) } ) ),
		             worst( sigma( fluxH, { largest( fluxH[ 10 ] ) } ) ),
		             worst( sigma( trace, { largest( trace[ 10 ] ) } ) ) );
		std::vector<double> const point = sigma( axisPoint, { 1.0 } );
		std::vector<double> const psis = sigma( axisPsi, { std::abs( span ) } );
		std::printf( "    axis: R %.2e m, z %.2e m; of the span: psi_h there %.2e, psi~ there %.2e; "
		             "at a fixed point psi~ %.2e, psi_h %.2e, psi* %.2e\n", point[ 0 ], point[ 1 ], psis[ 0 ],
		             psis[ 1 ], psis[ 2 ], psis[ 3 ], psis[ 4 ] );
		std::printf( "    axis element changes %zu over %zu elements (", elementChanges, distinct.size() );
		for ( int e : distinct )
			std::printf( " %d", e );
		std::printf( " ); z from %.3e to %.3e; overshoot from %.2e to %.2e\n",
		             std::min_element( axisPoint.begin(), axisPoint.end(),
		                               []( auto const &x, auto const &y ) { return x[ 1 ] < y[ 1 ]; } )->at( 1 ),
		             std::max_element( axisPoint.begin(), axisPoint.end(),
		                               []( auto const &x, auto const &y ) { return x[ 1 ] < y[ 1 ]; } )->at( 1 ),
		             *std::min_element( overshoots.begin(), overshoots.end() ),
		             *std::max_element( overshoots.begin(), overshoots.end() ) );
		std::printf( "    family nodes, worst position: %.2e m\n", worst( sigma( nodes, { 1.0 } ) ) );
		{
			// Not noise: what taking psi_ax from psi~ MOVES, at state 10.
			double move = 0.0;
			for ( std::size_t i = 0; i < located[ 10 ].size(); ++i )
				move = std::max( move, std::abs( tildeAxis[ 10 ][ i ] - located[ 10 ][ i ] )
				                       /std::abs( located[ 10 ][ i ] ) );
			std::printf( "    at state 10, psi~ - psi_h at the located axis is %.2e of the span, and the "
			             "family moves %.2e (worst average) taking psi_ax from psi~\n",
			             ( axisPsi[ 10 ][ 1 ] - axisPsi[ 10 ][ 0 ] )/span, move );
		}
		reportFamily( "psi~, located axis", located );
		reportFamily( "psi~, psi~ at the axis", tildeAxis );
		reportFamily( "psi~, fixed axis point", fixedAxis );
		reportFamily( "psi*, located axis", broken );
	}
}

/*
 * ACCURACY, NOT NOISE: the two axis labels' families at A = 0.30 against a
 * fine reference ( n = 48, k = 4, psi~ at the axis; its two labels agree to
 * the printed figure ). meq::extractFluxSurfaces() on a tracer takes the
 * axis's psi as given, which is what lets both labels be read here.
 */
BOOST_AUTO_TEST_CASE( accuracy )
{
	meq::FluxFamilyOptions options;
	options.surfaces = 12;
	options.angles = 96;
	options.innerCut = 0.05;
	options.outerCut = 0.90;
	auto families = [ & ]( int n, int order )
	{
		meq::tests::crossing::Solved s( 0.30, n, order, meq::tests::recovery::configure );
		meq::GradShafranovSolver const &solver = s.theSolver();
		meq::CriticalPoint const axis = meq::CriticalPointFinder( solver ).findAxis();
		meq::ContourTracer const continuous( solver, meq::Potential::Continuous );
		meq::CriticalPoint tilde = axis;
		tilde.psi = valueAt( solver.continuousPotential(), s.theMesh(), axis.radius, axis.z );
		return std::make_pair( averages( meq::extractFluxSurfaces( continuous, axis, solver.psiBoundary(), options ) ),
		                       averages( meq::extractFluxSurfaces( continuous, tilde, solver.psiBoundary(), options ) ) );
	};
	auto const reference = families( 48, 4 );
	auto difference = []( std::vector<double> const &a, std::vector<double> const &b, std::size_t firstSurface )
	{
		double w = 0.0;
		for ( std::size_t i = 6*firstSurface; i < a.size(); ++i )
			w = std::max( w, std::abs( a[ i ] - b[ i ] )/std::abs( b[ i ] ) );
		return w;
	};
	std::printf( "\n  accuracy against n = 48, k = 4 (its labels differ by %.1e): worst average, all surfaces | "
	             "surfaces 3 on\n", difference( reference.first, reference.second, 0 ) );
	for ( auto [ n, order ] : { std::pair{ 8, 2 }, { 16, 2 }, { 32, 2 }, { 8, 3 }, { 16, 3 } } )
	{
		auto const f = families( n, order );
		std::printf( "    n = %2d, k = %d: psi_h at the axis %.2e | %.2e; psi~ at the axis %.2e | %.2e\n", n, order,
		             difference( f.first, reference.second, 0 ), difference( f.first, reference.second, 3 ),
		             difference( f.second, reference.second, 0 ), difference( f.second, reference.second, 3 ) );
	}
}

BOOST_AUTO_TEST_CASE( normalisedN16k2 )
{
	int const n = 16;
	int const order = 2;
	meq::tests::Rectangle const box = meq::tests::standardBox();
	mfem::Mesh mesh = meq::tests::makeMesh( box, n );
	meq::FluxFamilyOptions options;
	options.surfaces = 12;
	options.angles = 96;
	options.innerCut = 0.05;
	options.outerCut = 0.90;

	std::vector<std::vector<double>> psiH, psiTilde, psiStarN, traceN, fluxN, axisPsi, bySolve, byTilde, byFinder;
	double span = 0.0;
	for ( int j = 0; j < states; ++j )
	{
		double const slope = 0.55*( 1.0 + j*spacing );
		std::vector<meq::Knot> knots;
		for ( double x : { 0.0, 1.0 } )
			knots.push_back( meq::Knot{ x, slope*( 2.20 + slope*x ), slope*slope } );
		auto pPrime = std::make_shared<meq::ConstantProfile const>( 0.45 );
		meq::NormalisedMHDSource source( pPrime, std::make_shared<meq::SplineProfile const>( std::move( knots ) ),
		                                 1.0, 1.0 );
		meq::GradShafranovSolver solver( mesh, order );
		mfem::FunctionCoefficient guess( [ box ]( mfem::Vector const &x )
		{
			return 0.30*std::sin( M_PI*( x( 0 ) - box.minRadius )/box.width() )
			       *std::sin( M_PI*( x( 1 ) - box.zMin )/box.height() );
		} );
		mfem::ConstantCoefficient zero( 0.0 );
		solver.setSource( source, 0.30 );
		solver.setBoundaryData( zero );
		solver.setInitialGuess( guess );
		solver.setNewtonControl( 1.0e-12, 1.0e-14, 60 );
		meq::tests::recovery::configure( solver );
		solver.solve();
		solver.postProcess();
		if ( j == 0 )
			meq::tests::recovery::reportRoutes( solver );

		psiH.push_back( dofs( solver.potential() ) );
		psiTilde.push_back( dofs( solver.continuousPotential() ) );
		psiStarN.push_back( dofs( solver.postProcessedPotential() ) );
		traceN.push_back( dofs( solver.trace() ) );
		fluxN.push_back( dofs( solver.flux() ) );
		meq::CriticalPoint const axis = meq::CriticalPointFinder( solver ).findAxis();
		double const tilde = valueAt( solver.continuousPotential(), mesh, axis.radius, axis.z );
		span = solver.psiAxis() - solver.psiBoundary();
		axisPsi.push_back( { solver.psiAxis(), axis.psi, tilde } );

		meq::ContourTracer const tracer( solver, meq::Potential::Continuous );
		meq::CriticalPoint a = axis;
		a.psi = solver.psiAxis();
		bySolve.push_back( averages( meq::extractFluxSurfaces( tracer, a, solver.psiBoundary(), options ) ) );
		a.psi = tilde;
		byTilde.push_back( averages( meq::extractFluxSurfaces( tracer, a, solver.psiBoundary(), options ) ) );
		byFinder.push_back( averages( meq::extractFluxSurfaces( tracer, axis, solver.psiBoundary(), options ) ) );
	}
	std::vector<double> const psis = sigma( axisPsi, { std::abs( span ) } );
	std::printf( "\n  normalised box, n = %d, k = %d, %s: sigma over %d states, gg' slope 1e-5 apart\n", n, order,
	             sizeof( mfem::real_t ) == 4 ? "single" : "double", states );
	auto both = [ & ]( std::vector<std::vector<double>> const &x )
	{
		std::vector<double> const v = sigma( x, { largest( x[ 10 ] ) } );
		double sum = 0.0;
		for ( double y : v )
			sum += y*y;
		return std::make_pair( worst( v ), std::sqrt( sum/static_cast<double>( v.size() ) ) );
	};
	std::printf( "    of their largest value, worst dof / root mean square: trace %.2e / %.2e, psi_h %.2e / %.2e, "
	             "q_h %.2e / %.2e, psi* %.2e / %.2e, psi~ %.2e / %.2e\n",
	             both( traceN ).first, both( traceN ).second, both( psiH ).first, both( psiH ).second,
	             both( fluxN ).first, both( fluxN ).second, both( psiStarN ).first, both( psiStarN ).second,
	             both( psiTilde ).first, both( psiTilde ).second );
	std::printf( "    of the span: the solve's psi_ax %.2e, the finder's psi_h %.2e, psi~ at the axis %.2e\n",
	             psis[ 0 ], psis[ 1 ], psis[ 2 ] );
	reportFamily( "labelled by the solve", bySolve );
	reportFamily( "labelled by psi~", byTilde );
	reportFamily( "labelled by the finder", byFinder );
}

BOOST_AUTO_TEST_CASE( n16k2 ) { measure( 16, 2 ); }
BOOST_AUTO_TEST_CASE( n12k2 ) { measure( 12, 2 ); }
BOOST_AUTO_TEST_CASE( n16k3 ) { measure( 16, 3 ); }
BOOST_AUTO_TEST_CASE( n32k2 ) { measure( 32, 2 ); }

