#define BOOST_TEST_MODULE MeqSourceFactoryTests

#include <boost/test/unit_test.hpp>

#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>

#include "meq/Config.hpp"
#include "meq/SourceFactory.hpp"

#include "analytic/ManufacturedNonlinear.hpp"

namespace
{
	/// A configuration built from a TOML string, so the fixtures are visible
	/// beside the assertions rather than living in a file somewhere.
	meq::Configuration configure( std::string const &body )
	{
		std::string const preamble =
			"[mesh]\n"
			"RMin = 0.6\nRMax = 1.4\nZMin = -0.6\nZMax = 0.6\n"
			"NR = 4\nNZ = 6\n"
			"\n[discretisation]\nPolynomialDegree = 2\n\n";
		return meq::Configuration::fromString( preamble + body, "<SourceFactoryTests>" );
	}

	/// A scratch file removed when the test leaves scope, whatever happens.
	class ScratchFile
	{
		public:
			explicit ScratchFile( std::string const &nameIn, std::string const &contents )
				: name( nameIn )
			{
				std::ofstream out( name );
				out << contents;
			}

			~ScratchFile()
			{
				std::remove( name.c_str() );
			}

			std::string const &fileName() const
			{
				return name;
			}

		private:
			std::string name;
	};

	/// A two-knot spline table in the format SplineProfile::fromStream reads.
	std::string const constantProfileTable =
		"# psi\tf\tf'\n"
		"0.0\t2.5\t0.0\n"
		"1.0\t2.5\t0.0\n"
		"\n";
}

BOOST_AUTO_TEST_SUITE( source_factory_tests )

BOOST_AUTO_TEST_CASE( soloviev_matches_a_directly_constructed_source )
{
	meq::Configuration const config = configure(
		"[source]\nType = \"soloviev\"\nA = -0.52\n" );

	auto const built = meq::makeSource( config.getSource(), config.getFileName() );
	meq::SolovievSource const direct( -0.52 );

	for ( double r = 0.7; r < 1.35; r += 0.13 )
	{
		for ( double z = -0.5; z < 0.55; z += 0.21 )
		{
			BOOST_TEST( built->f( r, z, 0.3 ) == direct.f( r, z, 0.3 ),
			            boost::test_tools::tolerance( 1.0e-15 ) );
			BOOST_TEST( built->dFdPsi( r, z, 0.3 ) == direct.dFdPsi( r, z, 0.3 ),
			            boost::test_tools::tolerance( 1.0e-15 ) );
		}
	}
}

/*
 * THE DRIFT TEST, and the reason SourceFactory.cpp is allowed to restate the
 * Example 5 formula that tests/analytic/ManufacturedNonlinear.hpp also carries.
 *
 * The library must not depend on a test fixture, and the fixture carries the
 * exact solution and its gradient, which the library has no use for. So the
 * source term exists twice. That is only safe while something checks the two
 * agree -- if it did not, one could be corrected and the other left, and a
 * convergence study would quietly measure a different problem from the one the
 * driver runs.
 */
BOOST_AUTO_TEST_CASE( manufactured_agrees_with_the_analytic_fixture )
{
	double const r0 = -0.5;
	double const kr = 1.15*3.14159265358979323846;
	double const kz = 1.15;

	meq::Configuration const config = configure(
		"[source]\nType = \"manufactured\"\n"
		"R0 = -0.5\nKr = 3.6128315516282616\nKz = 1.15\n"
		"\n[boundary]\nType = \"exact\"\n" );

	auto const built = meq::makeSource( config.getSource(), config.getFileName() );
	meq::analytic::ManufacturedNonlinear const fixture( r0, kr, kz );

	double worstF = 0.0;
	double worstDerivative = 0.0;

	for ( double r = 0.65; r < 1.36; r += 0.07 )
	{
		for ( double z = -0.55; z < 0.56; z += 0.11 )
		{
			// Probed away from the exact solution as well as on it: the two
			// implementations must agree as functions of psi, not merely at
			// the value that solves the problem.
			for ( double psi : { -0.4, 0.0, fixture.psi( r, z ), 0.9 } )
			{
				double const scale = 1.0 + std::fabs( fixture.f( r, z, psi ) );
				worstF = std::max( worstF,
					std::fabs( built->f( r, z, psi ) - fixture.f( r, z, psi ) )/scale );
				worstDerivative = std::max( worstDerivative,
					std::fabs( built->dFdPsi( r, z, psi ) - fixture.dFdPsi( r, z, psi ) )
					/( 1.0 + std::fabs( fixture.dFdPsi( r, z, psi ) ) ) );
			}
		}
	}

	BOOST_TEST_MESSAGE( "  factory vs analytic fixture: worst relative difference "
	                    << "in F " << worstF << ", in dF/dpsi " << worstDerivative );
	BOOST_TEST( worstF < 1.0e-14 );
	BOOST_TEST( worstDerivative < 1.0e-14 );
}

BOOST_AUTO_TEST_CASE( manufactured_derivative_matches_a_finite_difference )
{
	meq::Configuration const config = configure(
		"[source]\nType = \"manufactured\"\n"
		"R0 = -0.5\nKr = 3.6128315516282616\nKz = 1.15\n"
		"\n[boundary]\nType = \"exact\"\n" );

	auto const built = meq::makeSource( config.getSource(), config.getFileName() );

	double const step = 1.0e-6;
	for ( double r = 0.7; r < 1.35; r += 0.17 )
	{
		for ( double psi = -0.5; psi < 0.6; psi += 0.25 )
		{
			double const difference =
				( built->f( r, 0.13, psi + step ) - built->f( r, 0.13, psi - step ) )
				/( 2.0*step );
			BOOST_TEST( built->dFdPsi( r, 0.13, psi ) == difference,
			            boost::test_tools::tolerance( 1.0e-6 ) );
		}
	}
}

BOOST_AUTO_TEST_CASE( mhd_reads_its_profile_files )
{
	ScratchFile const pPrime( "meq_test_pprime.dat", constantProfileTable );
	ScratchFile const ggPrime( "meq_test_ggprime.dat", constantProfileTable );

	meq::Configuration const config = configure(
		"[source]\nType = \"mhd\"\n"
		"PPrimeFile = \"" + pPrime.fileName() + "\"\n"
		"GGPrimeFile = \"" + ggPrime.fileName() + "\"\n"
		"Mu0 = 1.0\n" );

	auto const built = meq::makeSource( config.getSource(), config.getFileName() );

	// F = mu0 r^2 p' + gg' = r^2 * 2.5 + 2.5 with mu0 = 1 and both profiles
	// constant at 2.5. Constant in psi, so the derivative vanishes.
	double const r = 1.2;
	BOOST_TEST( built->f( r, 0.0, 0.4 ) == r*r*2.5 + 2.5,
	            boost::test_tools::tolerance( 1.0e-12 ) );
	BOOST_TEST( built->dFdPsi( r, 0.0, 0.4 ) == 0.0,
	            boost::test_tools::tolerance( 1.0e-12 ) );
}

BOOST_AUTO_TEST_CASE( a_missing_profile_file_names_itself )
{
	meq::Configuration const config = configure(
		"[source]\nType = \"mhd\"\n"
		"PPrimeFile = \"meq_test_no_such_profile.dat\"\n"
		"GGPrimeFile = \"meq_test_no_such_profile.dat\"\n" );

	BOOST_CHECK_EXCEPTION( meq::makeSource( config.getSource(), config.getFileName() ),
		meq::ConfigError,
		[]( meq::ConfigError const &error )
		{
			// The failure is a configuration error, reported as one, naming
			// both the key at fault and the file it could not open -- rather
			// than an I/O error surfacing from somewhere inside the numerics.
			std::string const message( error.what() );
			return error.getKey() == "source.PPrimeFile"
			    && message.find( "meq_test_no_such_profile.dat" ) != std::string::npos;
		} );
}

BOOST_AUTO_TEST_SUITE_END()
