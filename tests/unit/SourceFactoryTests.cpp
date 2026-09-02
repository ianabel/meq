#define BOOST_TEST_MODULE MeqSourceFactoryTests

#include <boost/test/unit_test.hpp>

#include <cmath>
#include <cstdio>
#include <fstream>
#include <memory>
#include <string>

#include "meq/Config.hpp"
#include "meq/RotatingSource.hpp"
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

	/// n( psi ) = 1 + 2 psi. A Hermite cubic through two knots carrying the
	/// exact end slopes IS the linear function -- its cubic and quadratic
	/// coefficients cancel to zero in floating point, not merely to round-off
	/// -- so n'' is exactly zero and the rotating fixtures below have a
	/// derivative that can be written down.
	std::string const linearDensityTable =
		"# psi\tn\tn'\n"
		"0.0\t1.0\t2.0\n"
		"1.0\t3.0\t2.0\n"
		"\n";

	/// T( psi ) = 4, as a file rather than a constant, so that the same profile
	/// can be read both ways and the two required to agree.
	std::string const constantTemperatureTable =
		"# psi\tT\tT'\n"
		"0.0\t4.0\t0.0\n"
		"1.0\t4.0\t0.0\n"
		"\n";

	/// g g'( psi ) = 0.5 + psi, which varies with its argument -- which is what
	/// makes a change of normalisation visible as more than a scale factor.
	std::string const linearGGPrimeTable =
		"# psi\tgg'\td(gg')/dpsi\n"
		"0.0\t0.5\t1.0\n"
		"1.0\t1.5\t1.0\n"
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


/*
 * ROTATING SOURCES -- FLOW-PLAN.md's FL-8, from the file to the object.
 *
 * The fixtures below are not a plasma. The masses are kilogrammes of order one
 * and the temperatures Joules of order one, chosen so that the exponent of (96)
 * comes out at exactly 1.25 and the whole answer can be written down; a
 * physical set would put e^(m omega^2 r^2 / 2T) at the mercy of round-off in
 * the twentieth significant figure of the input, which tests the arithmetic of
 * the fixture rather than the arithmetic of the factory. Every physical claim
 * about meq::RotatingSource is made in RotatingSourceTests and in the
 * Rotating*Convergence studies; what is tested here is that the TOML reaches
 * the constructor intact.
 */

BOOST_AUTO_TEST_CASE( rotating_matches_a_hand_computed_value )
{
	ScratchFile const density( "meq_test_rot_density.dat", linearDensityTable );

	meq::Configuration const config = configure(
		"[source]\nType = \"rotating\"\n"
		"Omega = 2.0\n"
		"GGPrime = 0.5\n"
		"ReferenceRadius = 1.0\n"
		"Mu0 = 1.0\n"
		"\n[[source.species]]\n"
		"Name = \"ion\"\n"
		"Mass = 3.0\nCharge = 1.0\nTemperature = 4.0\n"
		"DensityFile = \"" + density.fileName() + "\"\n"
		"\n[[source.species]]\n"
		"Name = \"electron\"\n"
		"Mass = 1.0\nCharge = -1.0\nTemperature = 4.0\n"
		"Neutralising = true\n" );

	auto const built = meq::makeSource( config.getSource(), config.getFileName() );

	/*
	 * THE ARITHMETIC, from FLOW-PLAN.md section 4.1. With both temperatures
	 * constant and equal, omega constant, and Z = +1 / -1, the two species
	 * share one exponent and the pressure is P0( psi ) exp[ C ( r^2 - rRef^2 )
	 * / 2 ] with
	 *
	 *     C  = omega^2 ( Z_1 m_2 - Z_2 m_1 )/( Z_1 T_2 - Z_2 T_1 )
	 *        = omega^2 ( m_1 + m_2 )/( T_1 + T_2 )
	 *        = 4 * ( 3 + 1 )/( 4 + 4 )                              = 2
	 *     P0 = n_i0 ( T_i + T_e )                                   = 8 n_i0
	 *
	 * the second line using n_e0 = n_i0, which is what Neutralising asked the
	 * factory to derive. The table gives n_i0 = 1 + 2 psi, so P0' = 16, exactly
	 * and at every psi. At r = 1.5 and rRef = 1,
	 *
	 *     ( r^2 - rRef^2 )/2 = ( 2.25 - 1 )/2                       = 0.625
	 *     exponent           = C * 0.625                            = 1.25
	 *     F = mu0 r^2 P0' e^1.25 + g g'
	 *       = 1 * 2.25 * 16 * 3.4903429574597902 + 0.5
	 *       = 125.65234646862629 + 0.5                    = 126.15234646862629
	 *
	 * and dF/dpsi vanishes: P0 is linear, so P0'' = 0, and C is constant, so
	 * the exponent contributes nothing either. That the derivative of a
	 * deliberately affine case is exactly zero is worth asserting -- it is the
	 * Solov'ev rung of the ladder in CLAUDE.md's Testing stance, one step down
	 * from RotatingSourceTests' finite-difference sweep.
	 */
	double const r = 1.5;
	BOOST_TEST( built->f( r, 0.1, 0.4 ) == 126.15234646862629,
	            boost::test_tools::tolerance( 1.0e-12 ) );
	BOOST_TEST( built->dFdPsi( r, 0.1, 0.4 ) == 0.0,
	            boost::test_tools::tolerance( 1.0e-12 ) );

	// P0' is constant here, so F is too -- which says the density table was
	// read as the linear profile it is rather than clamped or defaulted.
	BOOST_TEST( built->f( r, -0.2, 0.9 ) == 126.15234646862629,
	            boost::test_tools::tolerance( 1.0e-12 ) );

	// At r = rRef the exponent vanishes by the gauge, so F is the static
	// answer: mu0 rRef^2 P0' + g g' = 16 + 0.5.
	BOOST_TEST( built->f( 1.0, 0.0, 0.4 ) == 16.5,
	            boost::test_tools::tolerance( 1.0e-12 ) );

	// The Neutralising species got the density charge neutrality implies,
	// rather than none: n_e0 = -( Z_i/Z_e ) n_i0 = n_i0.
	auto const *rotating = dynamic_cast<meq::RotatingSource const *>( built.get() );
	BOOST_TEST_REQUIRE( rotating != nullptr );
	BOOST_TEST_REQUIRE( rotating->species().size() == 2u );
	BOOST_TEST( rotating->referenceRadius() == 1.0 );
	BOOST_TEST( ( *rotating->species()[ 1 ].density )( 0.4 ) == 1.8,
	            boost::test_tools::tolerance( 1.0e-14 ) );
	BOOST_TEST( ( *rotating->species()[ 0 ].density )( 0.4 ) == 1.8,
	            boost::test_tools::tolerance( 1.0e-14 ) );
}

/*
 * A SCALE ACTUALLY SCALES, measured through the answer rather than through the
 * profile alone.
 *
 * The same file read twice, once at TemperatureScale = 1 and once at 2. A table
 * arrives in whatever units its author wrote it in -- keV is the usual one --
 * and editing the file to suit meq would make the file a function of which code
 * reads it, so the conversion belongs in the configuration. If it went missing
 * the run would still converge, to an equilibrium of the wrong plasma.
 *
 * Doubling T moves the answer twice over, in opposite directions, which is what
 * makes this a check on the profile rather than on one multiplication:
 *
 *     C  = omega^2 ( m_1 + m_2 )/( T_1 + T_2 )  HALVES,  2    -> 1
 *     P0'= n_i0' ( T_i + T_e )                  DOUBLES, 16   -> 32
 *     F  = 2.25 * 32 * e^0.625 + 0.5
 *        = 134.51370893512 + 0.5                       = 135.01370893512
 */
BOOST_AUTO_TEST_CASE( a_temperature_scale_scales_the_temperature )
{
	ScratchFile const density( "meq_test_scale_density.dat", linearDensityTable );
	ScratchFile const temperature( "meq_test_scale_temperature.dat", constantTemperatureTable );

	auto const source = [ & ]( std::string const &scale )
	{
		return
			"[source]\nType = \"rotating\"\n"
			"Omega = 2.0\nGGPrime = 0.5\nReferenceRadius = 1.0\nMu0 = 1.0\n"
			"\n[[source.species]]\n"
			"Mass = 3.0\nCharge = 1.0\n"
			"TemperatureFile = \"" + temperature.fileName() + "\"\n"
			"TemperatureScale = " + scale + "\n"
			"DensityFile = \"" + density.fileName() + "\"\n"
			"\n[[source.species]]\n"
			"Mass = 1.0\nCharge = -1.0\n"
			"TemperatureFile = \"" + temperature.fileName() + "\"\n"
			"TemperatureScale = " + scale + "\n"
			"Neutralising = true\n";
	};

	meq::Configuration const plain = configure( source( "1.0" ) );
	meq::Configuration const scaled = configure( source( "2.0" ) );

	auto const plainSource = meq::makeSource( plain.getSource(), plain.getFileName() );
	auto const scaledSource = meq::makeSource( scaled.getSource(), scaled.getFileName() );

	// At scale 1 the file reproduces the constant-key fixture above to the last
	// bit, which is what says the two spellings of a profile are one profile.
	BOOST_TEST( plainSource->f( 1.5, 0.1, 0.4 ) == 126.15234646862629,
	            boost::test_tools::tolerance( 1.0e-12 ) );
	BOOST_TEST( scaledSource->f( 1.5, 0.1, 0.4 ) == 135.01370893512,
	            boost::test_tools::tolerance( 1.0e-12 ) );

	// And directly: the scale reached the profile, at the value and not merely
	// somewhere in the answer.
	auto const *rotating = dynamic_cast<meq::RotatingSource const *>( scaledSource.get() );
	BOOST_TEST_REQUIRE( rotating != nullptr );
	BOOST_TEST( ( *rotating->species()[ 0 ].temperature )( 0.4 ) == 8.0,
	            boost::test_tools::tolerance( 1.0e-14 ) );
	BOOST_TEST( rotating->species()[ 0 ].temperature->prime( 0.4 ) == 0.0,
	            boost::test_tools::tolerance( 1.0e-14 ) );
}

/*
 * makeNormalisedSource, and what setNormalisation does to the answer.
 *
 * A normalised source has the form F( r, z, psi ) = H( r, z, psi/psi_ax )/psi_ax,
 * so psi_ax enters TWICE -- once as the profile argument and once as the
 * chain-rule factor. The g g' table below is deliberately linear rather than
 * constant so that both are visible: were only the factor moving, halving
 * psi_ax would exactly double the answer, and 632.76/2 = 316.38 is not the
 * 315.88 measured.
 *
 *     psi = 0.08, psi_ax = 0.2  ->  Psi = 0.4, g g' = 0.9
 *         F = ( 125.65234646862629 + 0.9 )/0.2      = 632.7617323431315
 *     psi = 0.08, psi_ax = 0.4  ->  Psi = 0.2, g g' = 0.7
 *         F = ( 125.65234646862629 + 0.7 )/0.4      = 315.88086617156574
 *
 * with 125.652... the mu0 r^2 P0' e^1.25 of the case above, unchanged because
 * P0' is constant in Psi.
 */
BOOST_AUTO_TEST_CASE( a_normalised_rotating_source_answers_for_its_normalisation )
{
	ScratchFile const density( "meq_test_norm_density.dat", linearDensityTable );
	ScratchFile const ggPrime( "meq_test_norm_ggprime.dat", linearGGPrimeTable );

	meq::Configuration const config = configure(
		"[source]\nType = \"rotating\"\n"
		"Omega = 2.0\n"
		"GGPrimeFile = \"" + ggPrime.fileName() + "\"\n"
		"ReferenceRadius = 1.0\nMu0 = 1.0\n"
		"Normalised = true\nPsiAxis = 0.2\n"
		"\n[[source.species]]\n"
		"Mass = 3.0\nCharge = 1.0\nTemperature = 4.0\n"
		"DensityFile = \"" + density.fileName() + "\"\n"
		"\n[[source.species]]\n"
		"Mass = 1.0\nCharge = -1.0\nTemperature = 4.0\n"
		"Neutralising = true\n" );

	BOOST_TEST( config.getSource().isNormalised() == true );
	BOOST_TEST( config.getSource().psiAxisGuess() == 0.2 );

	std::shared_ptr<meq::NormalisedSource> built =
		meq::makeNormalisedSource( config.getSource(), config.getFileName() );
	BOOST_TEST_REQUIRE( built != nullptr );

	// The guess out of the file is where it starts.
	BOOST_TEST( built->normalisation() == 0.2, boost::test_tools::tolerance( 1.0e-14 ) );
	BOOST_TEST( built->f( 1.5, 0.1, 0.08 ) == 632.7617323431315,
	            boost::test_tools::tolerance( 1.0e-12 ) );

	// And the solver moves it, which is what the non-const return is for: this
	// is the call GradShafranovSolver makes before every residual evaluation.
	built->setNormalisation( 0.4 );
	BOOST_TEST( built->normalisation() == 0.4, boost::test_tools::tolerance( 1.0e-14 ) );
	BOOST_TEST( built->f( 1.5, 0.1, 0.08 ) == 315.88086617156574,
	            boost::test_tools::tolerance( 1.0e-12 ) );
}

/// The other normalised source out of the same factory, which is the half that
/// says the normalised plumbing is not rotation's alone. With both tables
/// constant at 2.5 and mu0 = 1, H = r^2 2.5 + 2.5 whatever Psi is, so the whole
/// psi_ax dependence is the chain-rule factor and can be read off.
BOOST_AUTO_TEST_CASE( a_normalised_mhd_source_comes_out_of_the_same_factory )
{
	ScratchFile const pPrime( "meq_test_norm_pprime.dat", constantProfileTable );
	ScratchFile const ggPrime( "meq_test_norm_mhd_ggprime.dat", constantProfileTable );

	meq::Configuration const config = configure(
		"[source]\nType = \"mhd\"\n"
		"PPrimeFile = \"" + pPrime.fileName() + "\"\n"
		"GGPrimeFile = \"" + ggPrime.fileName() + "\"\n"
		"Mu0 = 1.0\n"
		"Normalised = true\nPsiAxis = 0.5\n" );

	std::shared_ptr<meq::NormalisedSource> built =
		meq::makeNormalisedSource( config.getSource(), config.getFileName() );
	BOOST_TEST_REQUIRE( built != nullptr );

	double const r = 1.2;
	BOOST_TEST( built->f( r, 0.0, 0.1 ) == ( r*r*2.5 + 2.5 )/0.5,
	            boost::test_tools::tolerance( 1.0e-12 ) );

	built->setNormalisation( 0.25 );
	BOOST_TEST( built->f( r, 0.0, 0.1 ) == ( r*r*2.5 + 2.5 )/0.25,
	            boost::test_tools::tolerance( 1.0e-12 ) );
}

/*
 * THE TWO FACTORIES REFUSE EACH OTHER'S CONFIGURATIONS, and the first of the
 * two refusals is the one that matters.
 *
 * meq::NormalisedRotatingSource IS-A meq::Source, so returning one from
 * makeSource would compile and would solve -- with psi_ax frozen at the guess
 * for ever, because nothing would ever call setNormalisation() on it. The run
 * would converge, print rates, and answer a problem the file did not describe.
 * That is the failure mode this project keeps meeting and it is why the guard
 * is a throw rather than a cast.
 *
 * Both messages name the function to use instead, because the reader of this
 * error is holding a file that is right and calling the wrong entry point.
 */
BOOST_AUTO_TEST_CASE( a_normalised_source_must_not_come_out_of_makeSource )
{
	ScratchFile const density( "meq_test_refuse_density.dat", linearDensityTable );

	meq::Configuration const config = configure(
		"[source]\nType = \"rotating\"\n"
		"Omega = 2.0\nGGPrime = 0.5\nReferenceRadius = 1.0\nMu0 = 1.0\n"
		"Normalised = true\nPsiAxis = 0.2\n"
		"\n[[source.species]]\n"
		"Mass = 3.0\nCharge = 1.0\nTemperature = 4.0\n"
		"DensityFile = \"" + density.fileName() + "\"\n"
		"\n[[source.species]]\n"
		"Mass = 1.0\nCharge = -1.0\nTemperature = 4.0\n"
		"Neutralising = true\n" );

	BOOST_CHECK_EXCEPTION( meq::makeSource( config.getSource(), config.getFileName() ),
		meq::ConfigError,
		[]( meq::ConfigError const &error )
		{
			std::string const message( error.what() );
			return error.getKey() == "source.Normalised"
			    && message.find( "makeNormalisedSource" ) != std::string::npos
			    && message.find( "setSource( NormalisedSource &, double )" ) != std::string::npos;
		} );
}

BOOST_AUTO_TEST_CASE( makeNormalisedSource_refuses_an_unnormalised_configuration )
{
	meq::Configuration const config = configure(
		"[source]\nType = \"soloviev\"\nA = -0.52\n" );

	BOOST_CHECK_EXCEPTION( meq::makeNormalisedSource( config.getSource(), config.getFileName() ),
		meq::ConfigError,
		[]( meq::ConfigError const &error )
		{
			std::string const message( error.what() );
			return error.getKey() == "source.Normalised"
			    && message.find( "makeSource" ) != std::string::npos;
		} );
}

/*
 * THE SPECIES CAP, AND WHERE IT ACTUALLY FIRES TODAY.
 *
 * meq::RotatingSource refuses more than meq::maxSpecies = 8, so that the
 * per-quadrature-point work allocates nothing; the factory translates that
 * std::invalid_argument into a ConfigError naming source.species, which is
 * what makes a library rule read as the configuration error it is.
 *
 * ConfigTests asserts the same refusal at PARSE time, which is what Config.hpp
 * documents ("Between two and meq::maxSpecies"), and that assertion is RED --
 * the parser checks the floor and not the ceiling. This case is deliberately
 * kept beside it rather than instead of it: it pins the behaviour that exists,
 * and the red one names the behaviour that is wanted. When the parser gains the
 * check, this case goes on passing, because the factory is still entitled to
 * refuse a set nobody parsed.
 */
BOOST_AUTO_TEST_CASE( the_factory_refuses_more_species_than_maxSpecies )
{
	// BUILT BY HAND RATHER THAN PARSED, DELIBERATELY. The parser refuses an
	// over-long species list itself now, so this configuration cannot be
	// reached through configure() at all -- ConfigTests owns that half. What is
	// left to test here is the factory's OWN guard, which is not redundant: a
	// SourceConfig is a plain struct with public members, so a caller that
	// builds one directly (a driver assembling a scan, a future GUI) bypasses
	// the parser entirely and would otherwise reach meq::RotatingSource's
	// std::invalid_argument instead of a ConfigError naming the key.
	meq::SourceConfig config;
	config.type = meq::SourceType::Rotating;

	meq::RotatingParameters parameters;
	parameters.referenceRadius = 1.0;
	parameters.mu0 = 1.0;
	parameters.ggPrime = 0.5;

	meq::SpeciesParameters electrons;
	electrons.mass = 1.0;
	electrons.charge = -8.0;
	electrons.temperature = 4.0;
	electrons.neutralising = true;
	parameters.species.push_back( electrons );

	for ( int i = 0; i < 8; ++i )
	{
		meq::SpeciesParameters ion;
		ion.mass = 3.0;
		ion.charge = 1.0;
		ion.temperature = 4.0;
		ion.density = 1.0;
		parameters.species.push_back( ion );
	}

	config.parameters = parameters;
	BOOST_TEST_REQUIRE( config.getRotating().species.size() == 9u );

	BOOST_CHECK_EXCEPTION( meq::makeSource( config, "hand-built.toml" ),
		meq::ConfigError,
		[]( meq::ConfigError const &error )
		{
			std::string const message( error.what() );
			return error.getKey() == "source.species"
			    && message.find( "maxSpecies" ) != std::string::npos;
		} );
}

BOOST_AUTO_TEST_SUITE_END()
