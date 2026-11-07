#define BOOST_TEST_MODULE QDrivenScalePin
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <memory>
#include <sstream>
#include <string>
#include <variant>
#include <vector>

#include "meq/Dataset.hpp"
#include "meq/Profiles.hpp"
#include "meq/Run.hpp"

#include "SourceTree.hpp"

/*
 * THE q-DRIVEN SOLVE'S SCALE PINS COMPOSE -- MEASUREMENTS.md M-QSCALE-3.
 *
 * Equilibria sharing p( Psi_N ) and q( Psi_N ) form a one-parameter family, and
 * a q-driven run names exactly one of R0B0, ToroidalFlux and PlasmaCurrent to
 * select the member (meq::PinnedToroidalFieldMap). Each pin is a different
 * square system. THE TEST THAT MATTERS is that they pose ONE problem: solve
 * once, read the answer's own g on the boundary, toroidal flux and current,
 * and re-solve from cold under each pin set to what was read. All three must
 * return the same equilibrium.
 *
 * THE TARGET IS NOT REPRESENTABLE, ON PURPOSE. q-driven-q.dat was measured
 * from an exactly quadratic g^2, so at degree 2 every q condition holds at the
 * answer and any choice of the one that gives way returns it. Here the target
 * is that table times ( 1 + 0.2 Psi^3 ), whose g^2 no quadratic follows: the
 * fit leaves a misfit, it lies along whichever condition gave way, and two
 * pins that dropped different conditions return equilibria that differ by it.
 * Mutation-tested on this fixture (M-QSCALE-3): R0B0 giving way on its
 * axis-limit condition instead of y_2 puts its answer 0.13 of psi's peak from
 * the others, and PlasmaCurrent giving way on y_1 puts every answer 0.05 to
 * 0.13 apart; this case asks for 1e-6.
 *
 * Coarse on purpose (the q-driven box at RefinementLevels = 0, k = 2): the
 * agreement asked for is the LOOP's, between solves on one mesh, and does not
 * improve with resolution.
 */
namespace
{
	std::string configuration( std::string const &prefix, std::string const &pin )
	{
		return
			"[mesh]\nRMin = 0.6\nRMax = 1.4\nZMin = -0.6\nZMax = 0.6\nNR = 8\nNZ = 8\n"
			"RefinementLevels = 0\n\n"
			"[discretisation]\nPolynomialDegree = 2\nTau = 1.0\n\n"
			"[source]\nType = \"mhd\"\nNormalised = true\nPsiAxis = 0.3\nMu0 = 1.0\n"
			"PPrimeFile = \"" + meq::test::sourcePath( "examples/q-driven-pprime.dat" ) + "\"\n"
			"SafetyFactorFile = \"pin-q.dat\"\nSafetyFactorDegree = 2\n" + pin + "\n"
			"[boundary]\nType = \"zero\"\n\n"
			"[initialguess]\nType = \"bump\"\nAmplitude = 0.3\nCentreR = 1.0\nCentreZ = 0.0\n"
			"RadiusR = 0.4\nRadiusZ = 0.6\n\n"
			"[solver]\nNewtonMaxIterations = 60\nNewtonRelativeTolerance = 1.0e-12\n"
			"NewtonAbsoluteTolerance = 1.0e-14\n\n"
			"[output]\nDirectory = \".\"\nPrefix = \"" + prefix + "\"\nGridNR = 17\nGridNZ = 17\n";
	}

	struct Answer
	{
		std::vector<double> psi;
		std::vector<double> gSquared;
		double edgeField = 0.0;
		double toroidalFlux = 0.0;
		double plasmaCurrent = 0.0;
		int iterations = 0;
		int equilibria = 0;
	};

	double number( meq::Dataset const &data, std::string const &name )
	{
		meq::DatasetAttribute const *a = data.attribute( name );
		BOOST_TEST_REQUIRE( a != nullptr, "no attribute " << name );
		if ( auto const *d = std::get_if<std::vector<double>>( &a->value ) )
			return d->front();
		return static_cast<double>( std::get<std::vector<int>>( a->value ).front() );
	}

	Answer solve( std::string const &directory, std::string const &prefix,
	              std::string const &pin )
	{
		meq::RunOptions options;
		options.configPath = directory + "/" + prefix + ".toml";
		options.configText = configuration( prefix, pin );
		options.baseDirectory = directory;
		options.version = "QDrivenScalePin";
		meq::Run run( options );
		BOOST_TEST_REQUIRE( run.prepare() == meq::Solved );
		int const code = run.solve();
		BOOST_TEST_REQUIRE( code == meq::Solved,
			"the q-driven run under " << pin << " exits " << code );
		BOOST_TEST_REQUIRE( run.writeRestart( directory + "/" + prefix ) == meq::Solved );
		meq::Dataset const &data = *run.dataset();

		Answer out;
		meq::DatasetVariable const *psi = data.variable( "psi_coefficients" );
		BOOST_TEST_REQUIRE( psi != nullptr );
		out.psi = std::get<std::vector<double>>( psi->values );
		std::istringstream in( std::get<std::string>(
			data.attribute( "g_squared_coefficients" )->value ) );
		for ( double value = 0.0; in >> value; )
			out.gSquared.push_back( value );
		out.edgeField = number( data, "edge_toroidal_field" );
		out.toroidalFlux = number( data, "toroidal_flux" );
		out.plasmaCurrent = number( data, "plasma_current" );
		out.iterations = static_cast<int>( number( data, "toroidal_loop_iterations" ) );
		out.equilibria = static_cast<int>( number( data, "toroidal_loop_equilibria" ) );
		std::printf( "    %-34s %2d KINSOL iterations, %3d equilibria: R0B0 %.9f, "
		             "Phi %.9e, I_p %.9e\n", pin.c_str(), out.iterations, out.equilibria,
		             out.edgeField, out.toroidalFlux, out.plasmaCurrent );
		return out;
	}

	/// max | a - b |/max | a |, over the coefficients of psi_h.
	double apart( std::vector<double> const &a, std::vector<double> const &b )
	{
		BOOST_TEST_REQUIRE( a.size() == b.size() );
		double peak = 0.0, worst = 0.0;
		for ( std::size_t i = 0; i < a.size(); ++i )
		{
			peak = std::max( peak, std::abs( a[ i ] ) );
			worst = std::max( worst, std::abs( a[ i ] - b[ i ] ) );
		}
		return worst/peak;
	}

	std::string fixed( double value )
	{
		std::ostringstream out;
		out.precision( 17 );
		out << value;
		return out.str();
	}
}

BOOST_AUTO_TEST_CASE( theThreePinsReturnOneEquilibrium )
{
	std::string const directory = std::filesystem::absolute( "QDrivenScalePin-files" ).string();
	std::filesystem::create_directories( directory );
	{
		// q-driven-q.dat times ( 1 + 0.2 Psi^3 ), at its own abscissae, with
		// the product rule on the slopes.
		std::ifstream in( meq::test::sourcePath( "examples/q-driven-q.dat" ) );
		std::ofstream out( directory + "/pin-q.dat" );
		out.precision( 17 );
		std::string line;
		while ( std::getline( in, line ) )
		{
			std::istringstream fields( line );
			double psi = 0.0, q = 0.0, slope = 0.0;
			if ( line.empty() || line[ 0 ] == '#' || !( fields >> psi >> q >> slope ) )
				continue;
			double const factor = 1.0 + 0.2*psi*psi*psi;
			out << psi << " " << q*factor << " " << slope*factor + q*0.6*psi*psi << "\n";
		}
	}

	std::printf( "\n  THE SCALE PINS COMPOSE (M-QSCALE-3)\n" );
	// THE BASE UNDER A FOURTH POSING, so that each re-solve below starts
	// somewhere else: the current pin from g = 1.8, where its re-solve opens
	// at g = 2.2 and R0B0's at the base's own boundary g.
	Answer const base = solve( directory, "pin-base",
	                           "ToroidalFieldGuess = 1.8\nPlasmaCurrent = 5.0\n" );

	std::vector<std::pair<std::string, std::string>> const pins = {
		{ "pin-edge", "R0B0 = " + fixed( base.edgeField ) + "\n" },
		{ "pin-flux", "ToroidalFieldGuess = 2.2\nToroidalFlux = " + fixed( base.toroidalFlux ) + "\n" },
		{ "pin-current", "ToroidalFieldGuess = 2.2\nPlasmaCurrent = " + fixed( base.plasmaCurrent ) + "\n" } };

	std::vector<Answer> answers = { base };
	for ( auto const &[ prefix, pin ] : pins )
	{
		Answer const again = solve( directory, prefix, pin );
		double const psiApart = apart( base.psi, again.psi );
		double const gApart = apart( base.gSquared, again.gSquared );
		std::printf( "      psi %.3e and g^2 %.3e from the base, relative to their peaks\n",
		             psiApart, gApart );
		BOOST_TEST( psiApart < 1.0e-6,
			"re-solved under " << pin << " -- the value the base equilibrium reads -- "
			"psi_h is " << psiApart << " of its peak from the base. The pins pose "
			"different problems: each must drop the SAME q condition "
			"(meq::PinnedToroidalFieldMap), or the fit's misfit lands in different "
			"directions and so do the answers" );
		BOOST_TEST( gApart < 1.0e-6, "g^2's coefficients are " << gApart << " apart under " << pin );
		answers.push_back( again );
	}

	/*
	 * AND NONE OF THEM CRAWLED. The loop is a differenced Newton started at
	 * gg' = 0 and finishes in 4 or 5 iterations under every pin here
	 * (M-QSCALE-3); a crawl that converged would still pass the agreement
	 * above, and is a defect.
	 */
	for ( Answer const &a : answers )
		BOOST_TEST( a.iterations <= 8,
			"a pinned q-driven loop took " << a.iterations << " KINSOL iterations "
			"and " << a.equilibria << " equilibria: a crawl that converges is a defect "
			"report, not a slow case" );
}
