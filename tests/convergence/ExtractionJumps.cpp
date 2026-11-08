#define BOOST_TEST_MODULE ExtractionJumps
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <memory>
#include <string>
#include <vector>

#include "ExtractionCrossing.hpp"

/*
 * WHAT MAKES THE GEOMETRY JUMP, AND HOW THE JUMP SCALES WITH THE MESH -- the
 * instrument behind MEASUREMENTS.md's M-364 and M-365, and the by-hand study
 * that would catch a jump coming back. BUILT, NOT REGISTERED: it asserts
 * nothing, and each case is a measurement (`--run_test=<case>`).
 * ExtractionContinuity is the suite's assertion; this is the picture, on
 * either pair the extraction can root: psi~ and q~, its default, or psi* and
 * q*, broken across every face.
 *
 * One parameter, the source mix A of ExtractionCrossing.hpp, stepped
 * A = A0 + j * spacing for j = -10 ... 10, every state its own linear solve,
 * so the family is a function of A alone. On EVERY surface of the family, the
 * four averages a transport code reads (V'_rho, < |grad rho|^2 >, < R^-2 >,
 * < |grad psi|^2/R^2 >). A jump J between scan states j and j + 1 is the
 * increment there less the mean of its neighbours' increments -- smooth, that
 * is spacing^3 G''' of the increment, plus the noise; a jump puts J at j and
 * -J/2 either side. Then, per surface and step:
 *
 *   * which nodes changed element (FluxSurface::elements), band face, ray
 *     route or ray iteration count, and whether the axis element moved;
 *   * every step with a jump above the threshold, and whether a node of that
 *     surface changed element across it, or across a neighbouring step whose
 *     echo it is (attributed), or not -- and for one that is not, the node
 *     that moved and its decisions over four states;
 *   * every element change, and the jump it made.
 *
 * And the two numbers a differencing caller pays for: how often a central
 * difference at 1.5e-8 straddles a crossing, and what it is wrong by when it
 * does -- the straddle cases bisect up to four crossings to 1e-11 and take
 * that difference centred on each.
 */
namespace
{
	using Family = std::shared_ptr<meq::FluxSurfaceFamily const>;

	struct Quantities
	{
		// [ surface ][ quantity ]
		std::vector<std::vector<double>> values;
	};

	Quantities quantities( meq::FluxSurfaceFamily const &family )
	{
		Quantities out;
		double const span = family.psiBoundary - family.psiAxis;
		for ( meq::FluxSurface const &s : family.surfaces )
		{
			double const rho = std::sqrt( s.normalisedFlux );
			double const dRhoDPsi = 1.0/( 2.0*rho*span );
			out.values.push_back( { s.vPrime/std::abs( dRhoDPsi ),
			                        s.gradPsiSquared*dRhoDPsi*dRhoDPsi,
			                        s.inverseRSquared, s.gradPsiSquaredOverRSquared } );
		}
		return out;
	}

	char const *quantityName( std::size_t q )
	{
		static char const *const names[] = { "V'_rho", "<|grad rho|^2>", "<R^-2>",
		                                     "<|grad psi|^2/R^2>" };
		return names[ q ];
	}

	/// Every jump of one scan's families, attributed: what M-364's table
	/// and M-365's continuity column print. @a families are the states
	/// at A = A0 + j * spacing, j = -10 ... 10.
	void analyse( char const *name, std::vector<Family> const &families, double spacing )
	{
		std::vector<Quantities> g;
		for ( Family const &f : families )
			g.push_back( quantities( *f ) );
		std::size_t const surfaces = g.front().values.size();
		std::size_t const steps = families.size() - 1;
		std::size_t const angles = families.front()->surfaces.front().count();
		double const threshold = 1.0e-11;

		std::size_t axisChanges = 0;
		for ( std::size_t j = 0; j < steps; ++j )
			axisChanges += families[ j ]->axisElement != families[ j + 1 ]->axisElement;

		// The jump estimate per surface, quantity and step: the increment less
		// the mean of its neighbours', relative to the value, less its median
		// along the scan (the smooth third-derivative floor).
		std::vector<std::vector<std::vector<double>>> estimates( surfaces,
			std::vector<std::vector<double>>( 4, std::vector<double>( steps, 0.0 ) ) );
		for ( std::size_t i = 0; i < surfaces; ++i )
			for ( std::size_t q = 0; q < 4; ++q )
			{
				auto increment = [ & ]( std::size_t k )
				{
					return g[ k + 1 ].values[ i ][ q ] - g[ k ].values[ i ][ q ];
				};
				double const scale = std::abs( g[ 10 ].values[ i ][ q ] );
				std::vector<double> &e = estimates[ i ][ q ];
				for ( std::size_t j = 1; j + 1 < steps; ++j )
					e[ j ] = ( increment( j ) - 0.5*( increment( j - 1 ) + increment( j + 1 ) ) )/scale;
				std::vector<double> sorted( e.begin() + 1, e.end() - 1 );
				std::nth_element( sorted.begin(), sorted.begin() + sorted.size()/2, sorted.end() );
				double const median = sorted[ sorted.size()/2 ];
				for ( std::size_t j = 1; j + 1 < steps; ++j )
					e[ j ] -= median;
			}

		// Which ( surface, step ) pairs have a node changing element, so that
		// a jump can be attributed to a crossing at its own step or, when two
		// crossings one step either side sum their echoes into it, at both
		// neighbours.
		std::vector<std::vector<char>> crossed( surfaces, std::vector<char>( steps, 0 ) );
		for ( std::size_t i = 0; i < surfaces; ++i )
			for ( std::size_t j = 0; j < steps; ++j )
			{
				meq::FluxSurface const &a = families[ j ]->surfaces[ i ];
				meq::FluxSurface const &b = families[ j + 1 ]->surfaces[ i ];
				for ( std::size_t node = 0; node < a.count(); ++node )
					crossed[ i ][ j ] |= a.elements[ node ] != b.elements[ node ];
			}
		std::size_t neighbours = 0, edges = 0;

		std::vector<double> moves;
		std::size_t crossings = 0, crossingsMeasured = 0, crossingsBelow = 0;
		std::size_t jumps = 0, attributed = 0;
		std::size_t faceChanges = 0, routeChanges = 0, iterationChanges = 0, bandNodes = 0;
		std::vector<double> sizes;
		std::vector<double> errors;
		std::vector<std::string> unattributed;

		for ( std::size_t i = 0; i < surfaces; ++i )
		{
			for ( std::size_t j = 0; j < steps; ++j )
			{
				meq::FluxSurface const &a = families[ j ]->surfaces[ i ];
				meq::FluxSurface const &b = families[ j + 1 ]->surfaces[ i ];
				std::size_t changed = 0, routes = 0, iterations = 0;
				double moved = 0.0;
				for ( std::size_t node = 0; node < a.count(); ++node )
				{
					if ( a.elements[ node ] != b.elements[ node ] )
					{
						++changed;
						moved = std::max( moved, std::hypot( b.radius[ node ] - a.radius[ node ],
						                                     b.z[ node ] - a.z[ node ] ) );
					}
					faceChanges += a.bandFaces[ node ] != b.bandFaces[ node ];
					routes += a.rayRoutes[ node ] != b.rayRoutes[ node ];
					iterations += a.rayIterations[ node ] != b.rayIterations[ node ];
				}
				if ( j == 0 )
					bandNodes += static_cast<std::size_t>( a.extendedNodes );
				crossings += changed;
				routeChanges += routes;
				iterationChanges += iterations;

				// The jump at this step, where both neighbouring increments exist:
				// a step counts where its estimate is the local maximum, so that
				// a jump's -J/2 echoes either side are not counted again.
				if ( j == 0 || j + 1 == steps )
					continue;
				double worst = 0.0, error = 0.0;
				std::size_t worstQuantity = 0;
				bool peak = false;
				for ( std::size_t q = 0; q < 4; ++q )
				{
					std::vector<double> const &e = estimates[ i ][ q ];
					double const jump = std::abs( e[ j ] );
					if ( jump > worst )
					{
						worst = jump;
						worstQuantity = q;
						peak = jump >= 1.5*std::abs( e[ j - 1 ] ) && jump >= 1.5*std::abs( e[ j + 1 ] );
						// What a central difference at 1.5e-8 straddling it is
						// wrong by, relative to the derivative: J / ( 2 delta ) / G'.
						double const derivative = std::abs( g[ j + 1 ].values[ i ][ q ]
						                                    - g[ j ].values[ i ][ q ] )/spacing;
						double const scale = std::abs( g[ 10 ].values[ i ][ q ] );
						error = jump*scale/( 2.0*1.5e-8 )/std::max( derivative, 1.0e-300 );
					}
				}
				if ( changed > 0 )
					moves.push_back( moved );
				if ( changed > 0 )
				{
					++crossingsMeasured;
					if ( worst > threshold )
					{
						sizes.push_back( worst );
						errors.push_back( error );
					}
					else
						++crossingsBelow;
				}
				if ( worst > threshold && peak )
				{
					++jumps;
					if ( changed > 0 )
						++attributed;
					else if ( ( j == 1 && crossed[ i ][ 0 ] ) || ( j + 2 == steps && crossed[ i ][ steps - 1 ] ) )
						++edges;
					else if ( crossed[ i ][ j - 1 ] || crossed[ i ][ j + 1 ] )
						++neighbours;
					else
					{
						char line[ 256 ];
						std::snprintf( line, sizeof line, "surface %zu, step %+d -> %+d: %.2e of %s, "
						               "route changes %zu, iteration changes %zu", i,
						               static_cast<int>( j ) - 10, static_cast<int>( j ) - 9, worst,
						               quantityName( worstQuantity ), routes, iterations );
						unattributed.push_back( line );

						// The node whose move departs most from smooth at this step,
						// and the axis's own departure.
						auto departure = [ & ]( auto const &field, std::size_t node )
						{
							auto at = [ & ]( std::size_t k )
							{
								return field( *families[ k ], node );
							};
							return std::abs( ( at( j + 1 ) - at( j ) )
							                 - 0.5*( ( at( j ) - at( j - 1 ) ) + ( at( j + 2 ) - at( j + 1 ) ) ) );
						};
						auto nodeR = []( meq::FluxSurfaceFamily const &f, std::size_t n, std::size_t i0 )
						{
							return f.surfaces[ i0 ].radius[ n ];
						};
						std::size_t culprit = 0;
						double worstMove = 0.0;
						for ( std::size_t node = 0; node < a.count(); ++node )
						{
							double const dR = departure( [ & ]( meq::FluxSurfaceFamily const &f, std::size_t n )
							                             { return nodeR( f, n, i ); }, node );
							double const dZ = departure( [ & ]( meq::FluxSurfaceFamily const &f, std::size_t n )
							                             { return f.surfaces[ i ].z[ n ]; }, node );
							if ( std::hypot( dR, dZ ) > worstMove )
							{
								worstMove = std::hypot( dR, dZ );
								culprit = node;
							}
						}
						double const axisMove = std::hypot(
							departure( []( meq::FluxSurfaceFamily const &f, std::size_t ) { return f.axisR; }, 0 ),
							departure( []( meq::FluxSurfaceFamily const &f, std::size_t ) { return f.axisZ; }, 0 ) );
						double const axisPsi = departure( []( meq::FluxSurfaceFamily const &f, std::size_t )
						                                  { return f.psiAxis; }, 0 );
						std::snprintf( line, sizeof line, "  node %zu departs %.2e m (elements %d -> %d, "
						               "routes %d -> %d, iterations %d -> %d); the axis departs %.2e m, "
						               "psi_ax %.2e; traced points %zu -> %zu", culprit, worstMove,
						               a.elements[ culprit ], b.elements[ culprit ],
						               static_cast<int>( a.rayRoutes[ culprit ] ),
						               static_cast<int>( b.rayRoutes[ culprit ] ),
						               a.rayIterations[ culprit ], b.rayIterations[ culprit ], axisMove,
						               axisPsi, a.tracedPoints, b.tracedPoints );
						unattributed.push_back( line );
						for ( std::size_t k = j - 1; k <= j + 2; ++k )
						{
							meq::FluxSurface const &c = families[ k ]->surfaces[ i ];
							std::snprintf( line, sizeof line, "    state %+d: element %d, band face %d, "
							               "route %d, iterations %d, ( %.15f, %.15f ), stalled rays %d",
							               static_cast<int>( k ) - 10, c.elements[ culprit ],
							               c.bandFaces[ culprit ], static_cast<int>( c.rayRoutes[ culprit ] ),
							               c.rayIterations[ culprit ], c.radius[ culprit ], c.z[ culprit ],
							               c.stalledRays );
							unattributed.push_back( line );
						}
					}
				}
			}
		}

		auto quantile = []( std::vector<double> v, double p )
		{
			if ( v.empty() )
				return 0.0;
			std::sort( v.begin(), v.end() );
			return v[ static_cast<std::size_t>( p*static_cast<double>( v.size() - 1 ) ) ];
		};

		double const window = spacing*static_cast<double>( steps );
		double const rate = static_cast<double>( crossings )/window;
		std::printf( "\n  %s: %zu surfaces x %zu angles, %zu band nodes at the base\n",
		             name, surfaces, angles, bandNodes );
		std::printf( "    over a span of %.0e in A: %zu node element changes "
		             "(%.3g per unit of A), axis element changes %zu, band "
		             "face changes %zu, ray route changes %zu, ray iteration changes %zu\n",
		             window, crossings, rate, axisChanges, faceChanges, routeChanges,
		             iterationChanges );
		std::printf( "    steps with a jump above %.0e of the value: %zu, of which %zu coincide "
		             "with a node of that surface changing element\n", threshold, jumps,
		             attributed );
		std::printf( "    and %zu more next to a crossing on the scan's first or last step, whose "
		             "own jump is not estimated and whose echo this is, and %zu next to a crossing "
		             "on a neighbouring step whose own estimate another crossing on that surface "
		             "cancels -- the echo is then the peak\n", edges, neighbours );
		for ( std::string const &line : unattributed )
			std::printf( "      NOT attributed: %s\n", line.c_str() );
		std::printf( "    element changes with both neighbouring steps in the scan: %zu; their "
		             "jump, relative to the value: %zu below %.0e, median %.2e, 90%% %.2e, max %.2e\n",
		             crossingsMeasured, crossingsBelow, threshold, quantile( sizes, 0.5 ),
		             quantile( sizes, 0.9 ), quantile( sizes, 1.0 ) );
		std::printf( "    the crossing node's own move across the step: median %.2e m, max %.2e m\n",
		             quantile( moves, 0.5 ), quantile( moves, 1.0 ) );
		std::printf( "    a central difference at 1.5e-8 straddling one is wrong by, "
		             "relative to the derivative: median %.2e, max %.2e; it straddles one with "
		             "probability %.2e per column (all %zu surfaces read)\n",
		             quantile( errors, 0.5 ), quantile( errors, 1.0 ), rate*3.0e-8, surfaces );
	}

	double const centre = 0.30;
	double const spacing = 5.0e-4;

	std::vector<Family> scanAround( int n, int order, meq::Potential which )
	{
		std::vector<Family> families;
		for ( int j = -10; j <= 10; ++j )
			families.push_back( std::make_shared<meq::FluxSurfaceFamily const>(
				meq::tests::crossing::familyAt( centre + j*spacing, n, order, which ) ) );
		return families;
	}

	char const *pairName( meq::Potential which )
	{
		return which == meq::Potential::Continuous ? "psi~ and q~" : "psi* and q*";
	}

	/// M-364: the extraction along the scan.
	void scan( int n, int order, meq::Potential which )
	{
		char name[ 128 ];
		std::snprintf( name, sizeof name, "box, n = %d, k = %d, on %s", n, order, pairName( which ) );
		analyse( name, scanAround( n, order, which ), spacing );
	}
}

BOOST_AUTO_TEST_CASE( continuousN8k2 ) { scan( 8, 2, meq::Potential::Continuous ); }
BOOST_AUTO_TEST_CASE( continuousN16k2 ) { scan( 16, 2, meq::Potential::Continuous ); }
BOOST_AUTO_TEST_CASE( continuousN32k2 ) { scan( 32, 2, meq::Potential::Continuous ); }
BOOST_AUTO_TEST_CASE( continuousN8k3 ) { scan( 8, 3, meq::Potential::Continuous ); }
BOOST_AUTO_TEST_CASE( continuousN16k3 ) { scan( 16, 3, meq::Potential::Continuous ); }
BOOST_AUTO_TEST_CASE( brokenN8k2 ) { scan( 8, 2, meq::Potential::PostProcessed ); }
BOOST_AUTO_TEST_CASE( brokenN16k2 ) { scan( 16, 2, meq::Potential::PostProcessed ); }
BOOST_AUTO_TEST_CASE( brokenN32k2 ) { scan( 32, 2, meq::Potential::PostProcessed ); }
BOOST_AUTO_TEST_CASE( brokenN8k3 ) { scan( 8, 3, meq::Potential::PostProcessed ); }
BOOST_AUTO_TEST_CASE( brokenN16k3 ) { scan( 16, 3, meq::Potential::PostProcessed ); }

/*
 * THE STRADDLES -- M-365's measurement, kept as the by-hand study. Along the
 * scan of `analyse`, the first node to change element on each of up to four
 * surfaces; the crossing bisected to 1e-11 in A; the central difference at
 * 1.5e-8 centred on it, against the one-sided derivatives either side. On the
 * continuous pair it is off by the noise; on the broken pair, or after
 * anything that made the geometry jump again, by J/( 2 delta ).
 */
namespace
{
	/// Where along A a node of a surface changes element, located by
	/// bisection to @a resolution; then the central difference at a
	/// 1.5e-8 centred there, against the one-sided derivatives either side
	/// (differenced at 1e-6, 2e-6 off it). For a jump J the straddled
	/// difference is off by J/( 2 delta ); for a kink by at most half the
	/// change in G', whatever the step.
	struct Straddle
	{
		double error = 0.0;   ///< | D_straddled - ( D_L + D_R )/2 | / | G' |, worst quantity
		double kink = 0.0;    ///< | D_R - D_L | / | G' |, worst quantity
		double floor = 0.0;   ///< the noise the straddled difference carries, ~1e-14 | G |/( 2 delta | G' | )
	};

	Straddle straddle( int n, int order, meq::Potential which, std::size_t surface, std::size_t node,
	                   double lo, double hi )
	{
		auto at = [ & ]( double mix )
		{
			return std::make_shared<meq::FluxSurfaceFamily const>(
				meq::tests::crossing::familyAt( mix, n, order, which ) );
		};
		int const left = at( lo )->surfaces[ surface ].elements[ node ];
		while ( hi - lo > 1.0e-11 )
		{
			double const mid = 0.5*( lo + hi );
			if ( at( mid )->surfaces[ surface ].elements[ node ] == left )
				lo = mid;
			else
				hi = mid;
		}
		double const centre = 0.5*( lo + hi );
		double const delta = 1.5e-8;
		auto values = [ & ]( double h ) { return quantities( *at( h ) ).values[ surface ]; };
		std::vector<double> const minus = values( centre - delta ), plus = values( centre + delta );
		std::vector<double> const l3 = values( centre - 3.0e-6 ), l1 = values( centre - 1.0e-6 );
		std::vector<double> const r1 = values( centre + 1.0e-6 ), r3 = values( centre + 3.0e-6 );
		Straddle out;
		for ( std::size_t q = 0; q < 4; ++q )
		{
			double const dl = ( l1[ q ] - l3[ q ] )/2.0e-6;
			double const dr = ( r3[ q ] - r1[ q ] )/2.0e-6;
			double const reference = 0.5*( dl + dr );
			double const straddled = ( plus[ q ] - minus[ q ] )/( 2.0*delta );
			double const scale = std::max( std::abs( reference ), 1.0e-300 );
			if ( std::abs( straddled - reference )/scale > out.error )
			{
				out.error = std::abs( straddled - reference )/scale;
				out.floor = 1.0e-14*std::abs( plus[ q ] )/( 2.0*delta )/scale;
			}
			out.kink = std::max( out.kink, std::abs( dr - dl )/scale );
		}
		return out;
	}

	/// Up to four crossings along a scan, one per surface: ( surface, node, step ).
	std::vector<std::array<std::size_t, 3>> crossings( std::vector<Family> const &scan )
	{
		std::vector<std::array<std::size_t, 3>> found;
		for ( std::size_t i = 0; i < scan.front()->surfaces.size() && found.size() < 4; ++i )
			for ( std::size_t j = 0; j + 1 < scan.size() && ( found.empty() || found.back()[ 0 ] != i ); ++j )
				for ( std::size_t node = 0; node < scan[ j ]->surfaces[ i ].count(); ++node )
					if ( scan[ j ]->surfaces[ i ].elements[ node ] != scan[ j + 1 ]->surfaces[ i ].elements[ node ] )
					{
						found.push_back( { i, node, j } );
						break;
					}
		return found;
	}

	void straddles( int n, int order, meq::Potential which )
	{
		char name[ 128 ];
		std::snprintf( name, sizeof name, "box, n = %d, k = %d, on %s", n, order, pairName( which ) );
		std::vector<Family> const scan = scanAround( n, order, which );
		analyse( name, scan, spacing );
		for ( auto const &c : crossings( scan ) )
		{
			double const lo = centre + ( static_cast<double>( c[ 2 ] ) - 10.0 )*spacing;
			Straddle const st = straddle( n, order, which, c[ 0 ], c[ 1 ], lo, lo + spacing );
			std::printf( "    surface %zu node %zu: a central difference at 1.5e-8 straddling the crossing "
			             "is off by %.2e of the derivative (noise there ~%.1e); the kink, "
			             "| G'_R - G'_L |/| G' |, is %.2e\n", c[ 0 ], c[ 1 ], st.error, st.floor, st.kink );
		}
	}
}

BOOST_AUTO_TEST_CASE( straddleContinuousN8k2 ) { straddles( 8, 2, meq::Potential::Continuous ); }
BOOST_AUTO_TEST_CASE( straddleContinuousN16k2 ) { straddles( 16, 2, meq::Potential::Continuous ); }
BOOST_AUTO_TEST_CASE( straddleContinuousN16k3 ) { straddles( 16, 3, meq::Potential::Continuous ); }
BOOST_AUTO_TEST_CASE( straddleBrokenN8k2 ) { straddles( 8, 2, meq::Potential::PostProcessed ); }
BOOST_AUTO_TEST_CASE( straddleBrokenN16k2 ) { straddles( 16, 2, meq::Potential::PostProcessed ); }
BOOST_AUTO_TEST_CASE( straddleBrokenN16k3 ) { straddles( 16, 3, meq::Potential::PostProcessed ); }
