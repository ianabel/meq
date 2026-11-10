#define BOOST_TEST_MODULE QDrivenCurvedBoundary
#include <boost/test/unit_test.hpp>

#include <chrono>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <variant>
#include <vector>

#include "mfem.hpp"

#include "meq/Dataset.hpp"
#include "meq/Run.hpp"

#include "convergence/QDrivenTarget.hpp"

/*
 * THE q-DRIVEN SOLVE ON A CURVED FIXED BOUNDARY -- the band.
 *
 * On a curved Gamma, Omega_h is inscribed in it and the outer flux surfaces
 * lie partly in the band outside the mesh. The write phase traces them
 * through the band by the transfer lift; the q-driven loop's map must too, or
 * every family fails at its outer surfaces, every map evaluation is refused,
 * and the outer Newton stops on an exactly zero Jacobian (MEASUREMENTS.md
 * M-353) -- so a q-driven run could not be posed on any machine shape.
 *
 * The round trip q-driven.toml makes on its box, made on a Miller boundary:
 * an equilibrium from a closed-form g = 2.20 + 0.55 Psi, its q measured on
 * the same mesh through q = V' g < R^-2 >/4 pi^2 at the q-driven run's own
 * knots (QDrivenTarget.hpp), and the q-driven run of that target recovering
 * g^2 there. The discretisation error is common to both, so what is read is
 * the loop's accuracy.
 *
 * UNDER THE R0B0 PIN, g = 2.20 on the boundary. The target is exact at every
 * knot and g^2 is quadratic, which the knots' cubic reproduces, so every row
 * holds at the answer and the recovered values are the closed form's.
 */
namespace
{
	double const g0 = 2.20;
	double const g1 = 0.55;

	std::string configuration( std::string const &prefix, std::string const &toroidal,
	                           std::string const &surfaces = "FluxSurfaces = true\nFluxSurfaceCount = 24\n" )
	{
		return
			"[mesh]\nRMin = 0.7\nRMax = 2.3\nZMin = -1.9\nZMax = 1.9\nNR = 8\nNZ = 10\n"
			// Unrefined: 29 equilibria, and the band is the same band.
			"RefinementLevels = 0\n\n"
			"[discretisation]\nPolynomialDegree = 2\nTau = 1.0\n\n"
			"[source]\nType = \"mhd\"\nNormalised = true\nPsiAxis = 0.3\nMu0 = 1.0\n"
			"PPrimeFile = \"curved-pprime.dat\"\n" + toroidal + "\n"
			"[boundary]\nType = \"zero\"\n\n"
			"[boundary.shape]\nType = \"miller\"\nR0 = 1.5\nZ0 = 0.0\nMinorRadius = 0.5\n"
			"Elongation = 1.6\nTriangularity = 0.35\nSquareness = 0.0\n\n"
			"[initialguess]\nType = \"bump\"\nAmplitude = 0.3\nCentreR = 1.5\nCentreZ = 0.0\n"
			"RadiusR = 0.5\nRadiusZ = 0.8\n\n"
			"[solver]\nNewtonMaxIterations = 60\nNewtonRelativeTolerance = 1.0e-12\n"
			"NewtonAbsoluteTolerance = 1.0e-14\n\n"
			"[output]\nDirectory = \".\"\nPrefix = \"" + prefix + "\"\nGridNR = 33\nGridNZ = 33\n"
			+ surfaces;
	}

	std::unique_ptr<meq::Run> run( std::string const &text, std::string const &directory )
	{
		meq::RunOptions options;
		options.configPath = directory + "/curved.toml";
		options.configText = text;
		options.baseDirectory = directory;
		options.version = "QDrivenCurvedBoundary";
		return std::make_unique<meq::Run>( options );
	}

	std::vector<double> variable( meq::Dataset const &data, std::string const &name )
	{
		meq::DatasetVariable const *v = data.variable( name );
		BOOST_TEST_REQUIRE( v != nullptr, "no variable " << name );
		return std::get<std::vector<double>>( v->values );
	}
}

BOOST_AUTO_TEST_CASE( theLoopRecoversItsFieldOnACurvedBoundary )
{
	std::string const directory = std::filesystem::absolute( "QDrivenCurvedBoundary-files" ).string();
	std::filesystem::create_directories( directory );
	{
		std::ofstream pPrime( directory + "/curved-pprime.dat" );
		pPrime << "0.0 0.45 0.0\n1.0 0.45 0.0\n";
		// gg' = g1( g0 + g1 Psi ), exact in one Hermite pair.
		std::ofstream ggPrime( directory + "/curved-ggprime.dat" );
		ggPrime.precision( 17 );
		ggPrime << 0.0 << " " << g1*g0 << " " << g1*g1 << "\n"
		        << 1.0 << " " << g1*( g0 + g1 ) << " " << g1*g1 << "\n";
	}

	// ---- the knots, the reference and its q there, on this mesh -------------
	using Clock = std::chrono::steady_clock;
	Clock::time_point const started = Clock::now();
	{
		std::ofstream provisional( directory + "/curved-q.dat" );
		provisional << "0.0 1.0 0.0\n1.0 1.0 0.0\n";
	}
	std::string const driven = "SafetyFactorFile = \"curved-q.dat\"\nR0B0 = 2.2\n";
	std::vector<double> const knots = meq::test::knotsOf( directory, configuration( "curved-driven", driven ) );
	std::shared_ptr<meq::Dataset const> const surfaces = meq::test::forward(
		directory, "curved-reference",
		configuration( "curved-reference", "GGPrimeFile = \"curved-ggprime.dat\"\n",
		               meq::test::surfacesAt( knots, 4 ) ) );
	meq::test::writeTarget( *surfaces, []( double psi ) { return g0 + g1*psi; },
	                        directory + "/curved-q.dat" );

	// ---- the q-driven run of that target ------------------------------------
	auto drivenRun = run( configuration( "curved-driven", driven ), directory );
	BOOST_TEST_REQUIRE( drivenRun->prepare() == meq::Solved );
	Clock::time_point const loop = Clock::now();
	int const code = drivenRun->solve();
	std::printf( "\n  wall: the reference and its target %.1f s, the q-driven solve %.1f s\n",
	             std::chrono::duration<double>( loop - started ).count(),
	             std::chrono::duration<double>( Clock::now() - loop ).count() );
	BOOST_TEST_REQUIRE( code == meq::Solved,
		"the q-driven run on a curved boundary exits " << code << ": if every "
		"evaluation was refused, its families are traced without the band "
		"(meq::SafetyFactorEquilibrium::Options::gammaHMarker)" );
	BOOST_TEST_REQUIRE( drivenRun->writeRestart( directory + "/curved-driven" ) == meq::Solved );
	auto numbers = [ & ]( char const *name )
	{
		meq::DatasetAttribute const *a = drivenRun->dataset()->attribute( name );
		BOOST_TEST_REQUIRE( a != nullptr, "no attribute " << name );
		std::istringstream in( std::get<std::string>( a->value ) );
		std::vector<double> out;
		for ( double value = 0.0; in >> value; )
			out.push_back( value );
		return out;
	};
	std::vector<double> const label = numbers( "g_squared_knots" );
	std::vector<double> const value = numbers( "g_squared" );
	BOOST_TEST_REQUIRE( label.size() == knots.size() );
	BOOST_TEST_REQUIRE( value.size() == knots.size() );

	double worst = 0.0;
	for ( std::size_t i = 0; i < label.size(); ++i )
	{
		double const psi = 1.0 - label[ i ];
		double const exact = ( g0 + g1*psi )*( g0 + g1*psi );
		worst = std::max( worst, std::abs( value[ i ] - exact )/exact );
	}
	std::printf( "  q-DRIVEN ON A MILLER BOUNDARY: %zu knots, worst |g^2/g^2_exact - 1| %.3e\n",
	             label.size(), worst );
	// WHAT THE LOOP OWES, NOT A FIGURE IT HAPPENS TO REACH: the target is exact
	// at every knot and the spline represents the closed form exactly, so
	// nothing but the loop's own convergence stands between the two. 1e-7 is a
	// tenth of what MEQ owes, and a map biased by 1e-8 (M-372) is caught at it.
	// In single the inner solve is good to ten eps, 1.2e-6, and the bound is
	// 20 eps, the rule QDrivenScalePin states.
	double const bound = std::max( 1.0e-7,
		20.0*static_cast<double>( std::numeric_limits<mfem::real_t>::epsilon() ) );
	BOOST_TEST( worst < bound,
		"g^2 at the knots is " << worst << " from the closed form at worst, against " << bound
		<< ": the answer is not the one its target was measured from" );
}
