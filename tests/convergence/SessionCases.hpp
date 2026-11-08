#ifndef MEQ_TESTS_SESSION_CASES_HPP
#define MEQ_TESTS_SESSION_CASES_HPP

/*
 * THE SESSION'S SMOOTHNESS CASES, shared by SessionSmoothness (the assertion)
 * and SessionJumps (the instrument behind M-364), so that both pose the
 * same equilibria: rotating-normalised's box and fixed-d-tcv's MXH boundary
 * with the same RotatingSource, the D density as a table, and the Newton run
 * to the tolerance asked for.
 */

#include <boost/test/unit_test.hpp>

#include <cstdio>
#include <memory>
#include <regex>
#include <string>

#include "meq/Run.hpp"

#include "SourceTree.hpp"

namespace meq::test::session
{
	inline meq::ProfileTable omega( double value )
	{
		meq::ProfileTable table;
		table.knots = { meq::Knot{ -0.5, value, 0.0 }, meq::Knot{ 1.5, value, 0.0 } };
		return table;
	}

	/// rotating-density-normalised.dat, with the knot at Psi = 0.5 moved by
	/// the relative step h.
	inline meq::ProfileTable density( double h )
	{
		double const scale = 1.0e20;
		meq::ProfileTable table;
		table.knots = { meq::Knot{ -0.2, 0.1*scale, 2.0*scale }, meq::Knot{ 0.0, 0.5*scale, 2.0*scale },
		                meq::Knot{ 0.5, 1.5*scale*( 1.0 + h ), 2.0*scale },
		                meq::Knot{ 1.0, 2.5*scale, 2.0*scale }, meq::Knot{ 1.5, 3.5*scale, 2.0*scale } };
		return table;
	}

	inline std::string tcvRotating()
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
	inline std::string adapted( std::string text, double tolerance )
	{
		text = std::regex_replace( text, std::regex( R"(\nDensityScale = 1.0e20\n)" ), "\n" );
		char line[ 64 ];
		std::snprintf( line, sizeof line, "\nNewtonRelativeTolerance = %.1e\n", tolerance );
		return std::regex_replace( text, std::regex( R"(\nNewtonRelativeTolerance = [0-9.eE+-]+\n)" ),
		                           line );
	}

	inline std::unique_ptr<meq::Run> solvedRun( std::string const &path, std::string const &text )
	{
		meq::RunOptions options;
		options.configPath = path;
		options.configText = text;
		options.profiles = { { "species.D.Density", density( 0.0 ) } };
		options.version = "session";
		auto run = std::make_unique<meq::Run>( options );
		BOOST_TEST_REQUIRE( run->prepare() == meq::Solved );
		BOOST_TEST_REQUIRE( run->solve() == meq::Solved );
		return run;
	}
}

#endif // MEQ_TESTS_SESSION_CASES_HPP
