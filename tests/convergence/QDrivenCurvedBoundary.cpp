#define BOOST_TEST_MODULE QDrivenCurvedBoundary
#include <boost/test/unit_test.hpp>

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
#include "meq/Run.hpp"

/*
 * THE q-DRIVEN SOLVE ON A CURVED FIXED BOUNDARY -- the band.
 *
 * On a curved Gamma, Omega_h is inscribed in it and the outer flux surfaces
 * lie partly in the band outside the mesh. The write phase traces them
 * through the band by the transfer lift; the q-driven loop's map must too, or
 * every family fails at its outer surfaces, every map evaluation is refused,
 * and the outer Newton stops on an exactly zero Jacobian (MEASUREMENTS.md
 * M-FIX-2) -- so a q-driven run could not be posed on any machine shape.
 *
 * The round trip q-driven.toml makes on its box, made on a Miller boundary:
 * an equilibrium from a closed-form g = 2.20 + 0.55 Psi, its q measured on
 * the same mesh through q = V' g < R^-2 >/4 pi^2, and the q-driven run of
 * that target recovering g^2's three coefficients. The discretisation error
 * is common to both, so what is read is the loop's accuracy.
 */
namespace
{
	double const g0 = 2.20;
	double const g1 = 0.55;

	std::string configuration( std::string const &prefix, std::string const &toroidal )
	{
		return
			"[mesh]\nRMin = 0.7\nRMax = 2.3\nZMin = -1.9\nZMax = 1.9\nNR = 8\nNZ = 10\n"
			"RefinementLevels = 1\n\n"
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
			"FluxSurfaces = true\nFluxSurfaceCount = 24\n";
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

	// ---- the reference and its q, on this mesh ------------------------------
	auto reference = run( configuration( "curved-reference", "GGPrimeFile = \"curved-ggprime.dat\"\n" ),
	                      directory );
	BOOST_TEST_REQUIRE( reference->prepare() == meq::Solved );
	BOOST_TEST_REQUIRE( reference->solve() == meq::Solved );
	std::shared_ptr<meq::Dataset const> const surfaces =
		reference->fluxSurfaces( directory + "/curved-reference_surfaces.nc" );
	std::vector<double> const psiN = variable( *surfaces, "normalised_flux" );
	std::vector<double> const vPrime = variable( *surfaces, "V_prime" );
	std::vector<double> const inverseR2 = variable( *surfaces, "inverse_R_squared" );

	// The target in the SOURCE's Psi, one on the axis, ascending, with
	// centred slopes: q-driven-q.dat's construction.
	std::vector<double> psi, q;
	for ( std::size_t i = psiN.size(); i-- > 0; )
	{
		double const sourcePsi = 1.0 - psiN[ i ];
		psi.push_back( sourcePsi );
		q.push_back( vPrime[ i ]*( g0 + g1*sourcePsi )*inverseR2[ i ]/( 4.0*M_PI*M_PI ) );
	}
	{
		std::ofstream target( directory + "/curved-q.dat" );
		target.precision( 17 );
		std::size_t const n = psi.size();
		for ( std::size_t i = 0; i < n; ++i )
		{
			std::size_t const a = i == 0 ? 0 : i - 1;
			std::size_t const b = i + 1 == n ? n - 1 : i + 1;
			target << psi[ i ] << " " << q[ i ] << " "
			       << ( q[ b ] - q[ a ] )/( psi[ b ] - psi[ a ] ) << "\n";
		}
	}

	// ---- the q-driven run of that target ------------------------------------
	auto driven = run( configuration( "curved-driven",
		"SafetyFactorFile = \"curved-q.dat\"\nSafetyFactorDegree = 2\n"
		"ToroidalFieldGuess = 2.2\n" ), directory );
	BOOST_TEST_REQUIRE( driven->prepare() == meq::Solved );
	int const code = driven->solve();
	BOOST_TEST_REQUIRE( code == meq::Solved,
		"the q-driven run on a curved boundary exits " << code << ": if every map "
		"evaluation was refused, its families are traced without the band "
		"(Run::Impl::extractFamily())" );
	BOOST_TEST_REQUIRE( driven->writeRestart( directory + "/curved-driven" ) == meq::Solved );
	meq::DatasetAttribute const *attribute = driven->dataset()->attribute( "g_squared_coefficients" );
	BOOST_TEST_REQUIRE( attribute != nullptr );
	std::istringstream in( std::get<std::string>( attribute->value ) );
	std::vector<double> c;
	for ( double value = 0.0; in >> value; )
		c.push_back( value );
	BOOST_TEST_REQUIRE( c.size() == 3u );

	double const exact[ 3 ] = { g0*g0, 2.0*g0*g1, g1*g1 };
	std::printf( "\n  q-DRIVEN ON A MILLER BOUNDARY: g^2 = %.8f + %.8f Psi + %.8f Psi^2 "
	             "against %.4f + %.4f Psi + %.4f Psi^2\n", c[ 0 ], c[ 1 ], c[ 2 ],
	             exact[ 0 ], exact[ 1 ], exact[ 2 ] );
	for ( int k = 0; k < 3; ++k )
		BOOST_TEST( std::abs( c[ k ] - exact[ k ] ) < 1.0e-4*exact[ 0 ],
			"coefficient " << k << " is " << c[ k ] << " against " << exact[ k ] );
}
