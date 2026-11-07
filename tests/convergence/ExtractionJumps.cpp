#define BOOST_TEST_MODULE ExtractionJumps
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <memory>
#include <regex>
#include <string>
#include <vector>

#include "meq/EquilibriumSnapshot.hpp"
#include "meq/Run.hpp"

#include "SessionCases.hpp"
#include "SourceTree.hpp"

/*
 * WHAT MAKES THE GEOMETRY JUMP, AND HOW THE JUMP SCALES WITH THE MESH -- the
 * instrument behind MEASUREMENTS.md's M-SHAPE-1. BUILT, NOT REGISTERED: it
 * asserts nothing, and each case is a measurement (`--run_test=<case>`).
 *
 * One input, Omega scaled by 1 + h, h = j * 1e-4 for j = -10 ... 10, every
 * state a resolve() warm from the base, so G is a function of the table
 * alone. On EVERY surface of the family, the four quantities MaNTA reads
 * (V'_rho, < |grad rho|^2 >, < R^-2 >, < |grad psi|^2/R^2 >). A jump J
 * between scan states j and j + 1 is the increment there less the mean of
 * its neighbours' increments -- smooth, that is h^3 G''' ~ 1e-12 of the
 * increment, plus the noise; a jump puts J at j and -J/2 either side. Then,
 * per surface and step:
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
 * difference at MaNTA's step straddles a crossing, and what it is wrong by
 * when it does.
 */
namespace
{
	using namespace meq::test::session;

	struct Quantities
	{
		// [ surface ][ quantity ]
		std::vector<std::vector<double>> values;
	};

	Quantities quantities( meq::EquilibriumSnapshot const &snapshot )
	{
		Quantities out;
		meq::FluxSurfaceFamily const &family = *snapshot.family;
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

	void scan( char const *name, std::string const &path, std::string const &text )
	{
		auto run = solvedRun( path, text );
		auto tables = []( double h )
		{
			return meq::ProfileTables{ { "Omega", omega( 4.0e5*( 1.0 + h ) ) },
			                           { "species.D.Density", density( 0.0 ) } };
		};
		auto at = [ & ]( double h, std::shared_ptr<meq::EquilibriumSnapshot const> const &warm )
		{
			BOOST_TEST_REQUIRE( run->resolve( tables( h ), warm ) == meq::Solved );
			return run->snapshot();
		};

		std::shared_ptr<meq::EquilibriumSnapshot const> const base = at( 0.0, run->snapshot() );
		double const spacing = 1.0e-4;
		std::vector<std::shared_ptr<meq::EquilibriumSnapshot const>> states;
		std::vector<Quantities> g;
		for ( int j = -10; j <= 10; ++j )
		{
			states.push_back( at( j*spacing, base ) );
			g.push_back( quantities( *states.back() ) );
		}

		std::size_t const surfaces = g.front().values.size();
		std::size_t const steps = states.size() - 1;
		std::size_t const angles = states.front()->family->surfaces.front().count();
		double const threshold = 1.0e-11;

		std::size_t axisChanges = 0;
		for ( std::size_t j = 0; j < steps; ++j )
			axisChanges += states[ j ]->family->axisElement != states[ j + 1 ]->family->axisElement;

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
				meq::FluxSurface const &a = states[ j ]->family->surfaces[ i ];
				meq::FluxSurface const &b = states[ j + 1 ]->family->surfaces[ i ];
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
				meq::FluxSurface const &a = states[ j ]->family->surfaces[ i ];
				meq::FluxSurface const &b = states[ j + 1 ]->family->surfaces[ i ];
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
						// What a central difference at MaNTA's 1.5e-8 straddling it is
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
								return field( *states[ k ]->family, node );
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
							meq::FluxSurface const &c = states[ k ]->family->surfaces[ i ];
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
		std::printf( "    over a relative span of %.0e in Omega: %zu node element changes "
		             "(%.3g per unit of relative change), axis element changes %zu, band "
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
		std::printf( "    a central difference at MaNTA's 1.5e-8 straddling one is wrong by, "
		             "relative to the derivative: median %.2e, max %.2e; it straddles one with "
		             "probability %.2e per column (all %zu surfaces read)\n",
		             quantile( errors, 0.5 ), quantile( errors, 1.0 ), rate*3.0e-8, surfaces );
	}

	std::string box( int levels, int degree )
	{
		std::string text = adapted( meq::test::shippedConfig( "examples/rotating-normalised.toml" ),
		                            1.0e-16 );
		text = std::regex_replace( text, std::regex( R"(RefinementLevels = 1)" ),
		                           "RefinementLevels = " + std::to_string( levels ) );
		return std::regex_replace( text, std::regex( R"(PolynomialDegree = 2)" ),
		                           "PolynomialDegree = " + std::to_string( degree ) );
	}

	std::string tcv( int levels, int degree )
	{
		std::string text = adapted( tcvRotating(), 1.0e-16 );
		text = std::regex_replace( text, std::regex( R"(RefinementLevels = 1)" ),
		                           "RefinementLevels = " + std::to_string( levels ) );
		return std::regex_replace( text, std::regex( R"(PolynomialDegree = 2)" ),
		                           "PolynomialDegree = " + std::to_string( degree ) );
	}

	std::string const boxPath = meq::test::sourcePath( "examples/rotating-normalised.toml" );
	std::string const tcvPath = meq::test::sourcePath( "examples/fixed-d-tcv.toml" );
}

BOOST_AUTO_TEST_CASE( boxL1k2 ) { scan( "box, level 1, k = 2", boxPath, box( 1, 2 ) ); }
BOOST_AUTO_TEST_CASE( boxL2k2 ) { scan( "box, level 2, k = 2", boxPath, box( 2, 2 ) ); }
BOOST_AUTO_TEST_CASE( boxL3k2 ) { scan( "box, level 3, k = 2", boxPath, box( 3, 2 ) ); }
BOOST_AUTO_TEST_CASE( boxL1k3 ) { scan( "box, level 1, k = 3", boxPath, box( 1, 3 ) ); }
BOOST_AUTO_TEST_CASE( boxL2k3 ) { scan( "box, level 2, k = 3", boxPath, box( 2, 3 ) ); }
BOOST_AUTO_TEST_CASE( tcvL1k2 ) { scan( "TCV, level 1, k = 2", tcvPath, tcv( 1, 2 ) ); }
BOOST_AUTO_TEST_CASE( tcvL2k2 ) { scan( "TCV, level 2, k = 2", tcvPath, tcv( 2, 2 ) ); }
BOOST_AUTO_TEST_CASE( tcvL1k3 ) { scan( "TCV, level 1, k = 3", tcvPath, tcv( 1, 3 ) ); }
