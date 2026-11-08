#define BOOST_TEST_MODULE ExtractionJumps
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <functional>
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
 * instrument behind MEASUREMENTS.md's M-SHAPE-1 and M-SHAPE-2, and the
 * by-hand study that would catch a jump coming back. BUILT, NOT REGISTERED:
 * it asserts nothing, and each case is a measurement (`--run_test=<case>`).
 * SessionSmoothness is the suite's assertion; this is the picture, on the
 * run's own extraction -- psi~ and q~, the continuous pair. M-SHAPE-1's
 * numbers were taken on psi* and q*, at aa5adcb.
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

	using Family = std::shared_ptr<meq::FluxSurfaceFamily const>;

	/// Every jump of one scan's families, attributed: what M-SHAPE-1's table
	/// and M-SHAPE-2's continuity column print. @a families are the states
	/// at h = j * spacing, j = -10 ... 10.
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

	meq::ProfileTables tables( double h )
	{
		return meq::ProfileTables{ { "Omega", omega( 4.0e5*( 1.0 + h ) ) },
		                           { "species.D.Density", density( 0.0 ) } };
	}

	/// M-SHAPE-1: the shipped extraction along the scan.
	void scan( char const *name, std::string const &path, std::string const &text )
	{
		auto run = solvedRun( path, text );
		BOOST_TEST_REQUIRE( run->resolve( tables( 0.0 ), run->snapshot() ) == meq::Solved );
		std::shared_ptr<meq::EquilibriumSnapshot const> const base = run->snapshot();
		double const spacing = 1.0e-4;
		std::vector<Family> families;
		for ( int j = -10; j <= 10; ++j )
		{
			BOOST_TEST_REQUIRE( run->resolve( tables( j*spacing ), base ) == meq::Solved );
			families.push_back( run->snapshot()->family );
		}
		analyse( name, families, spacing );
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

/*
 * THE STRADDLES, ON THE RUN'S OWN EXTRACTION -- M-SHAPE-2's measurement, kept
 * as the by-hand study (MANTA-EQUILIBRIUM-PLAN.md section 4.5). Along the
 * scan of `analyse`, the first node to change element on each of up to four
 * surfaces; the crossing bisected to 1e-11; the central difference at
 * MaNTA's 1.5e-8 centred on it, against the one-sided derivatives either side.
 * On the continuous pair it is off by the noise; a band face switching, an
 * axis element switching, or anything else that made the geometry jump again
 * would read J/( 2 delta ) here. The prototype that compared psi~, q~ against
 * psi*, q* through a seam on meq::Run is at aa5adcb (M-SHAPE-2), and the axis
 * on q~'s root at 9a7dc4f (M-SHAPE-3).
 */
namespace
{
	using Extract = std::function<Family( meq::Run & )>;

	Extract const runFamily = []( meq::Run &run ) { return run.snapshot()->family; };

	/// Where along the table a node of a surface changes element, located by
	/// bisection to @a resolution; then the central difference at MaNTA's
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

	Straddle straddle( meq::Run &run, std::shared_ptr<meq::EquilibriumSnapshot const> const &base,
	                   Extract const &extract, std::size_t surface, std::size_t node, double lo, double hi )
	{
		auto at = [ & ]( double h )
		{
			BOOST_TEST_REQUIRE( run.resolve( tables( h ), base ) == meq::Solved );
			return extract( run );
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

	void straddles( char const *name, std::string const &path, std::string const &text )
	{
		auto run = solvedRun( path, text );
		BOOST_TEST_REQUIRE( run->resolve( tables( 0.0 ), run->snapshot() ) == meq::Solved );
		std::shared_ptr<meq::EquilibriumSnapshot const> const base = run->snapshot();

		// The noise: the base tables from the base and from a state 1e-4 away.
		Quantities const fromBase = quantities( *runFamily( *run ) );
		BOOST_TEST_REQUIRE( run->resolve( tables( 1.0e-4 ), base ) == meq::Solved );
		std::shared_ptr<meq::EquilibriumSnapshot const> const away = run->snapshot();
		BOOST_TEST_REQUIRE( run->resolve( tables( 0.0 ), away ) == meq::Solved );
		Quantities const fromAway = quantities( *runFamily( *run ) );
		double noise = 0.0;
		for ( std::size_t i = 0; i < fromBase.values.size(); ++i )
			for ( std::size_t q = 0; q < 4; ++q )
				noise = std::max( noise, std::abs( fromAway.values[ i ][ q ] - fromBase.values[ i ][ q ] )
				                         /std::abs( fromBase.values[ i ][ q ] ) );

		double const spacing = 1.0e-4;
		std::vector<Family> scan;
		for ( int j = -10; j <= 10; ++j )
		{
			BOOST_TEST_REQUIRE( run->resolve( tables( j*spacing ), base ) == meq::Solved );
			scan.push_back( runFamily( *run ) );
		}
		std::printf( "\n  %s: the same tables from two warm starts differ by %.2e (relative, worst)\n",
		             name, noise );
		analyse( name, scan, spacing );
		for ( auto const &c : crossings( scan ) )
		{
			double const lo = ( static_cast<double>( c[ 2 ] ) - 10.0 )*spacing;
			Straddle const st = straddle( *run, base, runFamily, c[ 0 ], c[ 1 ], lo, lo + spacing );
			std::printf( "    surface %zu node %zu: a central difference at 1.5e-8 straddling the crossing "
			             "is off by %.2e of the derivative (noise there ~%.1e); the kink, "
			             "| G'_R - G'_L |/| G' |, is %.2e\n", c[ 0 ], c[ 1 ], st.error, st.floor, st.kink );
		}
	}
}

BOOST_AUTO_TEST_CASE( straddleBoxL1k2 ) { straddles( "box, level 1, k = 2", boxPath, box( 1, 2 ) ); }
BOOST_AUTO_TEST_CASE( straddleBoxL2k2 ) { straddles( "box, level 2, k = 2", boxPath, box( 2, 2 ) ); }
BOOST_AUTO_TEST_CASE( straddleBoxL3k2 ) { straddles( "box, level 3, k = 2", boxPath, box( 3, 2 ) ); }
BOOST_AUTO_TEST_CASE( straddleBoxL1k3 ) { straddles( "box, level 1, k = 3", boxPath, box( 1, 3 ) ); }
BOOST_AUTO_TEST_CASE( straddleBoxL2k3 ) { straddles( "box, level 2, k = 3", boxPath, box( 2, 3 ) ); }
BOOST_AUTO_TEST_CASE( straddleTcvL1k2 ) { straddles( "TCV, level 1, k = 2", tcvPath, tcv( 1, 2 ) ); }
BOOST_AUTO_TEST_CASE( straddleTcvL2k2 ) { straddles( "TCV, level 2, k = 2", tcvPath, tcv( 2, 2 ) ); }
BOOST_AUTO_TEST_CASE( straddleTcvL1k3 ) { straddles( "TCV, level 1, k = 3", tcvPath, tcv( 1, 3 ) ); }

/*
 * M-SESSION-4's SWEEP, kept as a measurement: the central difference of the
 * geometry MaNTA reads on the innermost, middle and outermost surfaces and
 * psi_ax, h from 1e-2 to 1e-10, Omega and one density knot, each state a
 * resolve() warm from one base; the fine scan over +-1e-3 and its decision
 * report. The suite's SessionSmoothness asserts the property on one located
 * crossing in about a minute; this prints the whole picture, box and TCV,
 * in about ten.
 */
namespace sweep
{
	using namespace meq::test::session;

	struct Output
	{
		std::vector<std::string> names;
		std::vector<double> values;
	};

	Output geometry( meq::EquilibriumSnapshot const &snapshot )
	{
		Output out;
		BOOST_TEST_REQUIRE( snapshot.family );
		meq::FluxSurfaceFamily const &family = *snapshot.family;
		std::size_t const n = family.surfaces.size();
		double const span = family.psiBoundary - family.psiAxis;
		for ( std::size_t i : { std::size_t( 0 ), n/2, n - 1 } )
		{
			meq::FluxSurface const &s = family.surfaces[ i ];
			double const rho = std::sqrt( s.normalisedFlux );
			double const dRhoDPsi = 1.0/( 2.0*rho*span );
			char label[ 16 ];
			std::snprintf( label, sizeof label, "%.3f", rho );
			out.names.push_back( std::string( "V'_rho@" ) + label );
			out.values.push_back( s.vPrime/std::abs( dRhoDPsi ) );
			out.names.push_back( std::string( "<|grad rho|^2>@" ) + label );
			out.values.push_back( s.gradPsiSquared*dRhoDPsi*dRhoDPsi );
			out.names.push_back( std::string( "<R^-2>@" ) + label );
			out.values.push_back( s.inverseRSquared );
			out.names.push_back( std::string( "<|grad psi|^2/R^2>@" ) + label );
			out.values.push_back( s.gradPsiSquaredOverRSquared );
		}
		out.names.push_back( "psi_ax" );
		out.values.push_back( snapshot.psiAxis );
		return out;
	}

	using Perturbed = std::function<meq::ProfileTables( double )>;

	/*
	 * WHICH DECISION MAKES THE JUMP -- the instrument behind M-SHAPE-1. Along
	 * the scan, every node of every surface keeps the decisions that placed it
	 * (FluxSurface::elements, bandFaces, rayRoutes, rayIterations) and the
	 * axis its element. A jump in an output is a jump in some node's position
	 * or flux; the node is found as the one whose position has the largest
	 * second difference across the spike, and its decisions are printed on
	 * the three states either side. Along the whole scan, how many nodes ever
	 * change each decision.
	 */
	void reportDecisions( std::vector<std::shared_ptr<meq::EquilibriumSnapshot const>> const &states,
	                      std::size_t output, int at )
	{
		std::size_t const n = states.front()->family->surfaces.size();
		std::size_t const surfaces[] = { 0, n/2, n - 1 };
		std::size_t const which = output/4 < 3 ? surfaces[ output/4 ] : 0;

		std::size_t axisChanges = 0;
		for ( std::size_t j = 1; j < states.size(); ++j )
			if ( states[ j ]->family->axisElement != states[ j - 1 ]->family->axisElement )
				++axisChanges;
		std::printf( "    decisions along the scan: the axis element changes %zu times "
		             "(element %d at the base)\n", axisChanges,
		             states[ states.size()/2 ]->family->axisElement );

		for ( std::size_t i = 0; i < n; ++i )
		{
			std::size_t elementNodes = 0, faceNodes = 0, routeNodes = 0, iterationNodes = 0;
			std::size_t traced = 0;
			meq::FluxSurface const &first = states.front()->family->surfaces[ i ];
			for ( std::size_t node = 0; node < first.count(); ++node )
			{
				bool element = false, face = false, route = false, iterations = false;
				for ( std::size_t j = 1; j < states.size(); ++j )
				{
					meq::FluxSurface const &a = states[ j - 1 ]->family->surfaces[ i ];
					meq::FluxSurface const &b = states[ j ]->family->surfaces[ i ];
					element |= a.elements[ node ] != b.elements[ node ];
					face |= a.bandFaces[ node ] != b.bandFaces[ node ];
					route |= a.rayRoutes[ node ] != b.rayRoutes[ node ];
					iterations |= a.rayIterations[ node ] != b.rayIterations[ node ];
				}
				elementNodes += element;
				faceNodes += face;
				routeNodes += route;
				iterationNodes += iterations;
			}
			for ( std::size_t j = 1; j < states.size(); ++j )
				traced += states[ j ]->family->surfaces[ i ].tracedPoints
				          != states[ j - 1 ]->family->surfaces[ i ].tracedPoints;
			if ( elementNodes + faceNodes + routeNodes + iterationNodes + traced > 0 )
				std::printf( "      surface %2zu (rho %.3f): nodes changing element %zu, band face %zu, "
				             "route %zu, ray iterations %zu; traced point count changes %zu times\n",
				             i, first.radial, elementNodes, faceNodes, routeNodes, iterationNodes,
				             traced );
		}

		if ( at < 1 || at + 1 >= static_cast<int>( states.size() ) )
			return;
		auto surfaceAt = [ & ]( int j ) -> meq::FluxSurface const &
		{
			return states[ static_cast<std::size_t>( j ) ]->family->surfaces[ which ];
		};
		meq::FluxSurface const &mid = surfaceAt( at );
		double worst = 0.0;
		std::size_t culprit = 0;
		for ( std::size_t node = 0; node < mid.count(); ++node )
		{
			double const dR = surfaceAt( at + 1 ).radius[ node ] - 2.0*mid.radius[ node ]
			                  + surfaceAt( at - 1 ).radius[ node ];
			double const dZ = surfaceAt( at + 1 ).z[ node ] - 2.0*mid.z[ node ]
			                  + surfaceAt( at - 1 ).z[ node ];
			double const size = std::hypot( dR, dZ );
			if ( size > worst )
			{
				worst = size;
				culprit = node;
			}
		}
		std::printf( "      at the spike: surface %zu's node %zu moves by a second difference of %.2e m\n",
		             which, culprit, worst );
		for ( int j = at - 1; j <= at + 1; ++j )
		{
			meq::FluxSurface const &s = surfaceAt( j );
			std::printf( "        state %+d: element %d, band face %d, route %d, iterations %d, "
			             "( %.15f, %.15f ), axis element %d\n", j - 10, s.elements[ culprit ],
			             s.bandFaces[ culprit ], static_cast<int>( s.rayRoutes[ culprit ] ),
			             s.rayIterations[ culprit ], s.radius[ culprit ], s.z[ culprit ],
			             states[ static_cast<std::size_t>( j ) ]->family->axisElement );
		}
	}

	struct Result
	{
		std::string name;
		double floor;
		double spike;
		std::string spikeOutput;
		std::vector<std::vector<double>> differences;
	};
	std::vector<Result> results;

	/// The table of differences for one input, and the floor: what it prints
	/// is MEASUREMENTS.md's M-SESSION-4.
	void measure( char const *name, meq::Run &run, Perturbed const &tables )
	{
		auto at = [ & ]( double h, std::shared_ptr<meq::EquilibriumSnapshot const> const &warm )
		{
			BOOST_TEST_REQUIRE( run.resolve( tables( h ), warm ) == meq::Solved );
			return run.snapshot();
		};

		// The base is the tables at h = 0, solved as every perturbed state is.
		std::shared_ptr<meq::EquilibriumSnapshot const> const base = at( 0.0, run.snapshot() );
		Output const g0 = geometry( *base );

		// THE FLOOR: the base tables, solved warm from the base and from a
		// state one 1e-4 step away.
		std::shared_ptr<meq::EquilibriumSnapshot const> const away = at( 1.0e-4, base );
		Output const fromBase = geometry( *at( 0.0, base ) );
		Output const fromAway = geometry( *at( 0.0, away ) );
		double floor = 0.0;
		for ( std::size_t k = 0; k < g0.values.size(); ++k )
		{
			BOOST_TEST( fromBase.values[ k ] == g0.values[ k ] );
			floor = std::max( floor, std::abs( fromAway.values[ k ] - fromBase.values[ k ] )
			                         /std::abs( g0.values[ k ] ) );
		}

		std::printf( "\n  %s\n", name );
		std::printf( "    the same tables from two warm starts differ by %.2e "
		             "(relative, worst output): the map's noise\n", floor );
		std::printf( "    %8s %12s %12s %12s  %s\n", "h", "max |dD|/|D|",
		             "worst", "Newton +/-", "|D| at the worst" );

		std::vector<double> const steps = { 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10 };
		std::vector<std::vector<double>> differences;
		for ( double h : steps )
		{
			std::shared_ptr<meq::EquilibriumSnapshot const> const plus = at( h, base );
			std::shared_ptr<meq::EquilibriumSnapshot const> const minus = at( -h, base );
			Output const up = geometry( *plus );
			Output const down = geometry( *minus );
			std::vector<double> d( g0.values.size() );
			for ( std::size_t k = 0; k < d.size(); ++k )
				d[ k ] = ( up.values[ k ] - down.values[ k ] )/( 2.0*h );
			double worst = 0.0;
			std::size_t which = 0;
			if ( !differences.empty() )
				for ( std::size_t k = 0; k < d.size(); ++k )
				{
					double const scale = std::max( std::abs( d[ k ] ), 1.0e-300 );
					double const change = std::abs( d[ k ] - differences.back()[ k ] )/scale;
					if ( change > worst )
					{
						worst = change;
						which = k;
					}
				}
			std::printf( "    %8.0e %12.3e %12s %6d %5d  %.6e\n", h, worst,
			             differences.empty() ? "-" : g0.names[ which ].c_str(),
			             plus->newtonIterations, minus->newtonIterations, d[ which ] );
			differences.push_back( d );
		}
		std::printf( "    the derivative at h = 1e-4, by output:\n" );
		for ( std::size_t k = 0; k < g0.values.size(); ++k )
			std::printf( "      %-26s %14.6e  (value %.6e)\n", g0.names[ k ].c_str(),
			             differences[ 2 ][ k ], g0.values[ k ] );

		/*
		 * THE FINE SCAN: G at h = j * 1e-4, j = -10 ... 10, and its second
		 * differences. Smooth, those are all h^2 G'' -- 1e-8 of G'', the same
		 * along the scan to h^3 G''' -- plus the noise; a jump J between two
		 * samples puts +J and -J into two neighbouring ones. So the spike is
		 * each second difference's departure from their median.
		 */
		double const spacing = 1.0e-4;
		std::vector<Output> scan;
		std::vector<std::shared_ptr<meq::EquilibriumSnapshot const>> states;
		for ( int j = -10; j <= 10; ++j )
		{
			states.push_back( at( j*spacing, base ) );
			scan.push_back( geometry( *states.back() ) );
		}
		double spike = 0.0;
		std::size_t spikeOutput = 0;
		int spikeAt = 0;
		for ( std::size_t k = 0; k < g0.values.size(); ++k )
		{
			std::vector<double> second;
			for ( std::size_t j = 1; j + 1 < scan.size(); ++j )
				second.push_back( ( scan[ j + 1 ].values[ k ] - 2.0*scan[ j ].values[ k ]
				                    + scan[ j - 1 ].values[ k ] )/std::abs( g0.values[ k ] ) );
			std::vector<double> sorted = second;
			std::nth_element( sorted.begin(), sorted.begin() + sorted.size()/2, sorted.end() );
			double const median = sorted[ sorted.size()/2 ];
			for ( std::size_t j = 0; j < second.size(); ++j )
				if ( std::abs( second[ j ] - median ) > spike )
				{
					spike = std::abs( second[ j ] - median );
					spikeOutput = k;
					spikeAt = static_cast<int>( j + 1 ) - 10;
				}
		}
		std::printf( "    fine scan over +-1e-3: the largest departure of a second difference "
		             "from their median is %.2e of the value, %s at h = %+.0e\n", spike,
		             g0.names[ spikeOutput ].c_str(), spikeAt*spacing );
		reportDecisions( states, spikeOutput, spikeAt + 10 );
		results.push_back( Result{ name, floor, spike, g0.names[ spikeOutput ], differences } );

		// Back to the base, so the next input starts where this one did.
		BOOST_TEST_REQUIRE( run.resolve( tables( 0.0 ), base ) == meq::Solved );
	}

	void measureCase( char const *name, meq::Run &run )
	{
		measure( ( std::string( name ) + ": Omega, a constant table, scaled by 1 + h" ).c_str(),
		         run, []( double h )
		         {
			         return meq::ProfileTables{ { "Omega", omega( 4.0e5*( 1.0 + h ) ) },
			                                    { "species.D.Density", density( 0.0 ) } };
		         } );
		measure( ( std::string( name ) + ": n_D's knot at Psi = 0.5 moved by h" ).c_str(),
		         run, []( double h )
		         {
			         return meq::ProfileTables{ { "Omega", omega( 4.0e5 ) },
			                                    { "species.D.Density", density( h ) } };
		         } );
	}
}

BOOST_AUTO_TEST_CASE( sweepBox )
{
	auto run = solvedRun( meq::test::sourcePath( "examples/rotating-normalised.toml" ),
	                      adapted( meq::test::shippedConfig( "examples/rotating-normalised.toml" ), 1.0e-16 ) );
	sweep::measureCase( "rotating-normalised (rectangle), NewtonRelativeTolerance 1e-16", *run );
}

BOOST_AUTO_TEST_CASE( sweepTcv )
{
	auto run = solvedRun( meq::test::sourcePath( "examples/fixed-d-tcv.toml" ),
	                      adapted( tcvRotating(), 1.0e-16 ) );
	sweep::measureCase( "fixed-d-tcv's MXH boundary, rotating source, NewtonRelativeTolerance 1e-16", *run );
}
