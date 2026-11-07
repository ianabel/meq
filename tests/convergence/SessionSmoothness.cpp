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
	std::string tcvRotating()
	{
		// examples/fixed-d-tcv.toml's machine -- an MXH boundary, curved, so
		// the outer surfaces cross the band -- with rotating-normalised's
		// source: RotatingSource is the one MaNTA's tables feed.
		std::string text = meq::test::shippedConfig( "examples/fixed-d-tcv.toml" );
		std::string const density = meq::test::sourcePath( "examples/rotating-density-normalised.dat" );
		std::string const source =
			"[source]\n"
			"Type = \"rotating\"\n"
			"Normalised = true\n"
			"PsiAxis = 0.043\n"
			"ReferenceRadius = 0.9\n"
			"Omega = 4.0e5\n"
			"GGPrime = 0.02\n\n"
			"[[source.species]]\n"
			"Name = \"D\"\n"
			"Mass = 3.3435837768e-27\n"
			"Charge = 1.0\n"
			"Temperature = 1.0\n"
			"TemperatureScale = 1.602176634e-16\n"
			"DensityFile = \"" + density + "\"\n"
			"DensityScale = 1.0e20\n\n"
			"[[source.species]]\n"
			"Name = \"e\"\n"
			"Mass = 9.1093837015e-31\n"
			"Charge = -1.0\n"
			"Temperature = 0.8\n"
			"TemperatureScale = 1.602176634e-16\n"
			"Neutralising = true\n\n";
		text = std::regex_replace( text, std::regex( R"(\n\[source\]\n[\s\S]*?\n(?=\[solver\]))" ),
		                           "\n" + source );
		text = std::regex_replace( text, std::regex( R"(RefinementLevels = 2)" ),
		                           "RefinementLevels = 1" );
		text = std::regex_replace( text, std::regex( R"(\nAmplitude = [0-9.eE+-]+)" ),
		                           "\nAmplitude = 0.043" );
		return text;
	}

	/// The configuration with the D density's scale taken into its table, so
	/// the table can be handed in, and the Newton's relative tolerance set.
	std::string adapted( std::string text, double tolerance )
	{
		text = std::regex_replace( text, std::regex( R"(\nDensityScale = 1.0e20\n)" ), "\n" );
		char line[ 64 ];
		std::snprintf( line, sizeof line, "\nNewtonRelativeTolerance = %.1e\n", tolerance );
		return std::regex_replace( text, std::regex( R"(\nNewtonRelativeTolerance = [0-9.eE+-]+\n)" ),
		                           line );
	}

	meq::ProfileTable density( double h );

	std::unique_ptr<meq::Run> solvedRun( std::string const &path, std::string const &text )
	{
		meq::RunOptions options;
		options.configPath = path;
		options.configText = text;
		options.profiles = { { "species.D.Density", density( 0.0 ) } };
		options.version = "SessionSmoothness";
		auto run = std::make_unique<meq::Run>( options );
		BOOST_TEST_REQUIRE( run->prepare() == meq::Solved );
		BOOST_TEST_REQUIRE( run->solve() == meq::Solved );
		return run;
	}

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
		for ( int j = -10; j <= 10; ++j )
			scan.push_back( geometry( *at( j*spacing, base ) ) );
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
		results.push_back( Result{ name, floor, spike, g0.names[ spikeOutput ], differences } );

		// Back to the base, so the next input starts where this one did.
		BOOST_TEST_REQUIRE( run.resolve( tables( 0.0 ), base ) == meq::Solved );
	}

	meq::ProfileTable omega( double value )
	{
		meq::ProfileTable table;
		table.knots = { meq::Knot{ -0.5, value, 0.0 }, meq::Knot{ 1.5, value, 0.0 } };
		return table;
	}

	/// rotating-density-normalised.dat, with the knot at Psi = 0.5 moved by
	/// the relative step h.
	meq::ProfileTable density( double h )
	{
		double const scale = 1.0e20;
		meq::ProfileTable table;
		table.knots = { meq::Knot{ -0.2, 0.1*scale, 2.0*scale }, meq::Knot{ 0.0, 0.5*scale, 2.0*scale },
		                meq::Knot{ 0.5, 1.5*scale*( 1.0 + h ), 2.0*scale },
		                meq::Knot{ 1.0, 2.5*scale, 2.0*scale }, meq::Knot{ 1.5, 3.5*scale, 2.0*scale } };
		return table;
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
				"the jump over the step. MANTA-EQUILIBRIUM-PLAN.md section 4.5's frozen "
				"extraction is what removes it (manta/shape-derivative)." );
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
