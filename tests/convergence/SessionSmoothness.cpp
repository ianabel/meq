#define BOOST_TEST_MODULE SessionSmoothness
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
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
 * IS THE EQUILIBRIUM A SMOOTH FUNCTION OF ITS TABLES? -- MANTA-EQUILIBRIUM-
 * PLAN.md section 4.6, measured through the session before any derivative is
 * built on it.
 *
 * MaNTA differences its equilibrium model at relative steps of 1e-8, and its
 * Newton is pinned to the model being smooth and deterministic there. So: one
 * table entry perturbed by +-h, h from 1e-2 to 1e-10, each a resolve() warm
 * from the SAME snapshot (the answer is then a function of the table alone),
 * and the central difference of the geometry MaNTA would read at fixed
 * rho = sqrt( Psi_N ):
 *
 *     V'_rho = dV/drho,  < |grad rho|^2 >,  < R^-2 >,  < |grad psi|^2/R^2 >
 *
 * on the family's innermost, middle and outermost surfaces, and psi_ax.
 *
 * A smooth map gives the V: the difference settles as h^2 and then rises as
 * noise/h. A piecewise one plateaus or jumps. The noise is read separately,
 * as the spread of the same tables solved from two warm starts -- the floor
 * below which no difference means anything.
 */
namespace
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

namespace
{
	/*
	 * WHAT IS WANTED, which is what MaNTA's contract says: at its relative step
	 * of 1e-8, G differenced to better than 1e-4 of the derivative. That needs
	 * the map's noise below 1e-12 of its value -- a resolve converged to
	 * round-off, not to a tolerance, whatever warm start it had -- and no jump
	 * anywhere a step could straddle one: the scan's spike at the noise.
	 */
	void requireSmooth()
	{
		for ( Result const &r : results )
		{
			BOOST_TEST( r.floor < 1.0e-12,
				r.name << ": the same tables solved from two warm starts differ "
				"by " << r.floor << " of the value, so a difference at MaNTA's 1e-8 step "
				"carries " << r.floor/1.0e-8 << " of noise. A session must converge to "
				"round-off, not to a relative tolerance (NewtonRelativeTolerance below ten "
				"eps stops at the round-off floor, M-306)." );
			BOOST_TEST( r.spike < 1.0e-11,
				r.name << ": G jumps by " << r.spike << " of " << r.spikeOutput
				<< " between two tables 1e-4 apart -- a discrete decision in the solve or "
				"the extraction changing with the input. A step straddling it is wrong by "
				"the jump over the step. What removes it is an extraction continuous in the "
				"solved state, chosen self-consistently and not by freezing the choice: "
				"MANTA-EQUILIBRIUM-PLAN.md section 4.5 and M-SHAPE-1 weigh the candidates "
				"(manta/shape-derivative). The decision report printed above names the "
				"choice that changed." );
		}
		results.clear();
	}
}

BOOST_AUTO_TEST_CASE( theRectangleIsSmoothInItsTables )
{
	std::string const path = meq::test::sourcePath( "examples/rotating-normalised.toml" );
	// TO ROUND-OFF, as a session serving MaNTA solves: at the shipped 1e-10
	// the noise across warm starts is 2.8e-11 on the machine below, and a
	// perturbation under ~1e-9 takes no Newton step at all, so G is flat
	// there (M-SESSION-4).
	for ( double tolerance : { 1.0e-16 } )
	{
		auto run = solvedRun( path, adapted( meq::test::shippedConfig(
			"examples/rotating-normalised.toml" ), tolerance ) );
		char name[ 96 ];
		std::snprintf( name, sizeof name, "rotating-normalised (rectangle), "
		               "NewtonRelativeTolerance %.0e", tolerance );
		measureCase( name, *run );
		if ( tolerance < 1.0e-12 )
			requireSmooth();
		results.clear();
	}
}

BOOST_AUTO_TEST_CASE( theShapedMachineIsSmoothInItsTables )
{
	std::string const path = meq::test::sourcePath( "examples/fixed-d-tcv.toml" );
	for ( double tolerance : { 1.0e-16 } )
	{
		auto run = solvedRun( path, adapted( tcvRotating(), tolerance ) );
		char name[ 96 ];
		std::snprintf( name, sizeof name, "fixed-d-tcv's MXH boundary, rotating source, "
		               "NewtonRelativeTolerance %.0e", tolerance );
		measureCase( name, *run );
		if ( tolerance < 1.0e-12 )
			requireSmooth();
		results.clear();
	}
}
