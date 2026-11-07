#define BOOST_TEST_MODULE RunSession
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <memory>
#include <stdexcept>
#include <string>
#include <variant>
#include <vector>

#include "meq/EquilibriumSnapshot.hpp"
#include "meq/Run.hpp"

#include "SourceTree.hpp"

/*
 * THE SESSION -- MANTA-EQUILIBRIUM-PLAN.md section 4.1: meq::Run::resolve()
 * and meq::EquilibriumSnapshot.
 *
 * A transport code calls an equilibrium hundreds of times a transient with
 * tables that move a little each time, keeps several answers at once, and
 * pins its own Newton to the equilibrium being a deterministic function of
 * the input and the warm start. These are its contract, case by case, on
 * examples/rotating-normalised.toml -- RotatingSource, the source MaNTA's
 * tables feed -- and each is what PYTHON-PLAN.md's PY-3 asserts of the
 * file-restart route, made of the in-memory one:
 *
 *   - an unchanged resolve is the answer in no Newton step, to the bit;
 *   - a changed table resolves warm to the cold solve of the same tables;
 *   - a resolve IS the .nc restart of the same tables, to the bit;
 *   - a failed resolve leaves the run as it was, and it resolves after;
 *   - the answer depends on the tables and the warm start and on nothing
 *     else the run holds, and a snapshot is a value;
 *   - what a session cannot restart on is refused, with the run untouched.
 */
namespace
{
	std::string const rotating = meq::test::sourcePath( "examples/rotating-normalised.toml" );

	/// rotating-normalised's constant Omega as a table, in Psi, across the
	/// range its own density table covers.
	meq::ProfileTable omega( double value )
	{
		meq::ProfileTable table;
		table.abscissa = meq::ProfileAbscissa::Psi;
		table.knots = { meq::Knot{ -0.5, value, 0.0 }, meq::Knot{ 1.5, value, 0.0 } };
		return table;
	}

	std::unique_ptr<meq::Run> prepared( std::string const &config,
	                                    meq::ProfileTables const &profiles = {},
	                                    std::string const &restart = "" )
	{
		meq::RunOptions options;
		options.configPath = config;
		options.profiles = profiles;
		options.restartFrom = restart;
		options.version = "RunSession";
		auto run = std::make_unique<meq::Run>( options );
		BOOST_TEST_REQUIRE( run->prepare() == meq::Solved );
		return run;
	}

	std::unique_ptr<meq::Run> solved( std::string const &config,
	                                  meq::ProfileTables const &profiles = {},
	                                  std::string const &restart = "" )
	{
		auto run = prepared( config, profiles, restart );
		BOOST_TEST_REQUIRE( run->solve() == meq::Solved );
		return run;
	}

	bool sameBits( std::vector<double> const &a, std::vector<double> const &b )
	{
		return a.size() == b.size()
		       && std::memcmp( a.data(), b.data(), a.size()*sizeof( double ) ) == 0;
	}

	/// max |a - b| over max |b|.
	double relativeGap( std::vector<double> const &a, std::vector<double> const &b )
	{
		BOOST_TEST_REQUIRE( a.size() == b.size() );
		double gap = 0.0, scale = 0.0;
		for ( std::size_t i = 0; i < a.size(); ++i )
		{
			gap = std::max( gap, std::abs( a[ i ] - b[ i ] ) );
			scale = std::max( scale, std::abs( b[ i ] ) );
		}
		return gap/scale;
	}

	/// Every variable of two datasets, byte for byte, and every global
	/// attribute but those saying how long a phase took or which binary ran.
	void requireSameDataset( meq::Dataset const &a, meq::Dataset const &b )
	{
		BOOST_TEST_REQUIRE( a.variables().size() == b.variables().size() );
		for ( std::size_t i = 0; i < a.variables().size(); ++i )
		{
			meq::DatasetVariable const &x = a.variables()[ i ];
			meq::DatasetVariable const &y = b.variables()[ i ];
			BOOST_TEST_REQUIRE( x.name == y.name );
			bool same = x.values.index() == y.values.index();
			if ( same && std::holds_alternative<std::vector<double>>( x.values ) )
				same = sameBits( std::get<std::vector<double>>( x.values ),
				                 std::get<std::vector<double>>( y.values ) );
			else if ( same )
				same = x.values == y.values;
			BOOST_TEST( same, "variable " << x.name << " differs" );
		}
		auto timing = []( std::string const &name )
		{
			return name.find( "seconds" ) != std::string::npos
			       || name.find( "time" ) != std::string::npos
			       || name == "meq_version" || name == "history";
		};
		BOOST_TEST_REQUIRE( a.attributes().size() == b.attributes().size() );
		for ( std::size_t i = 0; i < a.attributes().size(); ++i )
		{
			meq::DatasetAttribute const &x = a.attributes()[ i ];
			meq::DatasetAttribute const &y = b.attributes()[ i ];
			BOOST_TEST_REQUIRE( x.name == y.name );
			if ( !timing( x.name ) )
				BOOST_TEST( ( x.value == y.value ), "attribute " << x.name << " differs" );
		}
	}

	std::string stem( std::string const &name )
	{
		std::filesystem::create_directories( "RunSession-files" );
		return std::filesystem::absolute( "RunSession-files/" + name ).string();
	}
}

BOOST_AUTO_TEST_CASE( anUnchangedResolveIsTheAnswerInNoNewtonStep )
{
	auto run = solved( rotating );
	std::shared_ptr<meq::EquilibriumSnapshot const> const first = run->snapshot();
	BOOST_TEST_REQUIRE( first->family );
	BOOST_TEST( first->sequence == 1u );

	BOOST_TEST_REQUIRE( run->resolve( {} ) == meq::Solved );
	std::shared_ptr<meq::EquilibriumSnapshot const> const again = run->snapshot();

	BOOST_TEST( again->sequence == 2u );
	BOOST_TEST( again->newtonIterations == 0,
		"an unchanged resolve took " << again->newtonIterations << " Newton steps; "
		"it starts AT the answer, so its first residual is the converged one" );
	BOOST_TEST( sameBits( again->state, first->state ),
		"an unchanged resolve moved the state by " << relativeGap( again->state, first->state ) );
	BOOST_TEST( again->psiAxis == first->psiAxis );
	BOOST_TEST_REQUIRE( again->family );
	BOOST_TEST( again->family->surfaces.size() == first->family->surfaces.size() );
	for ( std::size_t i = 0; i < first->family->surfaces.size(); ++i )
		BOOST_TEST( again->family->surfaces[ i ].vPrime
		            == first->family->surfaces[ i ].vPrime );
}

BOOST_AUTO_TEST_CASE( aChangedTableResolvesWarmToTheColdSolve )
{
	meq::ProfileTables const faster = { { "Omega", omega( 4.4e5 ) } };
	auto run = solved( rotating );
	BOOST_TEST_REQUIRE( run->resolve( faster ) == meq::Solved );
	std::shared_ptr<meq::EquilibriumSnapshot const> const warm = run->snapshot();

	auto cold = solved( rotating, faster );
	std::shared_ptr<meq::EquilibriumSnapshot const> const reference = cold->snapshot();

	double const gap = relativeGap( warm->state, reference->state );
	std::printf( "\n  Omega 10%% up, warm against cold: %.3e of the state's max, "
	             "psi_ax %.3e apart; Newton %d against %d\n", gap,
	             std::abs( warm->psiAxis - reference->psiAxis )/std::abs( reference->psiAxis ),
	             warm->newtonIterations, reference->newtonIterations );
	BOOST_TEST( gap < 1.0e-8 );
	BOOST_TEST( std::abs( warm->psiAxis - reference->psiAxis ) < 1.0e-9*std::abs( reference->psiAxis ) );
	BOOST_TEST( warm->newtonIterations < reference->newtonIterations,
		"the warm resolve took no fewer Newton steps than the cold solve" );
	BOOST_TEST( ( warm->profiles.count( "Omega" ) == 1 ) );
}

BOOST_AUTO_TEST_CASE( aResolveIsTheFileRestartToTheBit )
{
	meq::ProfileTables const faster = { { "Omega", omega( 4.4e5 ) } };

	auto session = solved( rotating );
	BOOST_TEST_REQUIRE( session->writeRestart( stem( "restart-from" ) ) == meq::Solved );
	BOOST_TEST_REQUIRE( session->resolve( faster ) == meq::Solved );

	// PY-3's route: a new run of the same configuration, restarted from the
	// file the last one wrote, with the new tables.
	auto restarted = solved( rotating, faster, stem( "restart-from" ) + ".nc" );

	std::shared_ptr<meq::EquilibriumSnapshot const> const a = session->snapshot();
	std::shared_ptr<meq::EquilibriumSnapshot const> const b = restarted->snapshot();
	BOOST_TEST( sameBits( a->state, b->state ),
		"the in-memory resolve and the file restart differ by "
		<< relativeGap( a->state, b->state ) << " of the state" );
	BOOST_TEST( a->psiAxis == b->psiAxis );
	BOOST_TEST( a->newtonIterations == b->newtonIterations );
}

BOOST_AUTO_TEST_CASE( aFailedResolveLeavesTheRunAsItWas )
{
	auto run = solved( rotating );
	std::shared_ptr<meq::EquilibriumSnapshot const> const before = run->snapshot();
	BOOST_TEST_REQUIRE( run->writeRestart( stem( "before" ) ) == meq::Solved );
	std::shared_ptr<meq::Dataset const> const written = run->dataset();

	// Omega at 3e6 does not solve (PythonResolve's failing case).
	int const code = run->resolve( { { "Omega", omega( 3.0e6 ) } } );
	BOOST_TEST_REQUIRE( code != meq::Solved );
	BOOST_TEST( code == meq::SolveFailed );

	BOOST_TEST( run->dataset() == written, "a failed resolve replaced dataset()" );
	BOOST_TEST( run->snapshot() == before,
		"a failed resolve replaced the snapshot the run holds" );
	BOOST_TEST_REQUIRE( run->writeRestart( stem( "after" ) ) == meq::Solved );
	requireSameDataset( *written, *run->dataset() );

	// And it resolves again from there.
	BOOST_TEST_REQUIRE( run->resolve( { { "Omega", omega( 4.4e5 ) } } ) == meq::Solved );
	BOOST_TEST( run->snapshot()->sequence == before->sequence + 1 );
}

BOOST_AUTO_TEST_CASE( theAnswerDependsOnTheTablesAndTheWarmStartAlone )
{
	meq::ProfileTables const one = { { "Omega", omega( 4.4e5 ) } };
	meq::ProfileTables const two = { { "Omega", omega( 3.6e5 ) } };

	auto run = solved( rotating );
	std::shared_ptr<meq::EquilibriumSnapshot const> const start = run->snapshot();
	BOOST_TEST_REQUIRE( run->resolve( one ) == meq::Solved );
	std::shared_ptr<meq::EquilibriumSnapshot const> const first = run->snapshot();
	std::vector<double> const firstState = first->state;

	BOOST_TEST_REQUIRE( run->resolve( two ) == meq::Solved );
	// From the start, with the run now holding another state.
	BOOST_TEST_REQUIRE( run->resolve( one, start ) == meq::Solved );
	std::shared_ptr<meq::EquilibriumSnapshot const> const repeat = run->snapshot();

	BOOST_TEST( sameBits( repeat->state, firstState ),
		"the same tables from the same warm start differ by "
		<< relativeGap( repeat->state, firstState )
		<< ": the answer depends on the state the run held when it was called" );
	BOOST_TEST( sameBits( first->state, firstState ), "a later resolve changed a snapshot" );
	BOOST_TEST( repeat->sequence == 4u );
}

BOOST_AUTO_TEST_CASE( whatASessionCannotRestartOnIsRefused )
{
	auto unsolved = prepared( rotating );
	BOOST_CHECK_THROW( unsolved->resolve( {} ), std::logic_error );
	BOOST_CHECK_THROW( unsolved->snapshot(), std::logic_error );

	auto one = solved( rotating );
	auto other = solved( rotating );
	BOOST_CHECK_THROW( one->resolve( {}, other->snapshot() ), std::invalid_argument );

	// A table no slot takes is prepare()'s refusal, and the run is untouched.
	std::shared_ptr<meq::EquilibriumSnapshot const> const held = one->snapshot();
	meq::ProfileTable nothing = omega( 1.0 );
	BOOST_TEST( one->resolve( { { "NoSuchSlot", nothing } } ) == meq::ConfigurationError );
	BOOST_TEST( one->snapshot() == held );
	BOOST_TEST( one->resolve( {} ) == meq::Solved );

	// [adaptivity] rebuilds the mesh a session would restart on.
	auto adaptive = solved( meq::test::sourcePath( "examples/miller-adaptive.toml" ) );
	std::shared_ptr<meq::EquilibriumSnapshot const> const adapted = adaptive->snapshot();
	BOOST_TEST( adaptive->resolve( {} ) == meq::ConfigurationError );
	BOOST_TEST( adaptive->snapshot() == adapted );
	BOOST_TEST( adaptive->writeRestart( stem( "adaptive" ) ) == meq::Solved );
}

/*
 * WHAT A RESOLVE COSTS AGAINST PY-3's FILE RESTART, which pays a whole
 * prepare() -- configuration, mesh, spaces -- and a file each way. Printed,
 * not asserted: no wall time on a shared machine is a measurement
 * (MEASUREMENTS.md M-SESSION-3 reads ratios).
 */
BOOST_AUTO_TEST_CASE( theCostOfAResolveAgainstTheFileRestart )
{
	using Clock = std::chrono::steady_clock;
	auto seconds = []( Clock::time_point from )
	{
		return std::chrono::duration<double>( Clock::now() - from ).count();
	};
	int const repeats = 4;
	auto run = solved( rotating );
	double inMemory = 0.0, viaFile = 0.0;
	for ( int i = 0; i < repeats; ++i )
	{
		meq::ProfileTables const step = { { "Omega", omega( 4.0e5*( 1.0 + 1.0e-3*( i + 1 ) ) ) } };
		// Both from the same state, to the same tables.
		Clock::time_point t0 = Clock::now();
		BOOST_TEST_REQUIRE( run->writeRestart( stem( "cost" ) ) == meq::Solved );
		auto next = solved( rotating, step, stem( "cost" ) + ".nc" );
		viaFile += seconds( t0 );

		t0 = Clock::now();
		BOOST_TEST_REQUIRE( run->resolve( step ) == meq::Solved );
		inMemory += seconds( t0 );
	}
	std::printf( "\n  a 0.1%% Omega step, mean of %d: resolve %.3f s, file restart %.3f s, "
	             "ratio %.2f\n", repeats, inMemory/repeats, viaFile/repeats,
	             viaFile/inMemory );
}
