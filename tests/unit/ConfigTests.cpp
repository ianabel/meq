/*
 * Tests for meq::Configuration -- the TOML schema in src/meq/Config.{hpp,cpp}.
 *
 * Three of these guard defects that the code this replaces actually had, and
 * they are the reason the file is worth its length:
 *
 *   * soloviev_example_parses, manufactured_example_parses. The old
 *     examples/helmholtz.conf could not be read by the parser shipped beside
 *     it: it wrote OutputMeshFile and PsiResolution where the parser demanded
 *     FinalMeshFile and CellSize, and omitted the NetCDFFile the parser
 *     required. Nothing in the tree noticed, because nothing ever loaded the
 *     example. These two cases load both shipped examples.
 *
 *   * integer_valued_keys_are_accepted_where_a_number_belongs. Every value used
 *     to be read with .as_floating(), so `RMin = 0` -- a perfectly ordinary
 *     thing to write -- threw where `RMin = 0.0` worked. toml11's own
 *     find<double> has the same behaviour, and find_or<double> is worse: on an
 *     integer node it silently returns the *default*.
 *
 *   * unknown_key_is_rejected. The cause of the first defect: a key that no
 *     reader looks at used to be ignored in silence.
 *
 * Note while reading: TOML key names are UpperCamelCase and the C++ members
 * they land in are lowerCamelCase, deliberately (Config.hpp says why). So a
 * test asserts on the key "mesh.RMin" and reads the value out of
 * getMesh().rMin, and neither spelling is a typo for the other.
 *
 * Boost.Test entry point: this translation unit supplies it, by defining
 * BOOST_TEST_MODULE below. That works both for the static library and for the
 * shared one (whose imported target, Boost::unit_test_framework, defines
 * BOOST_UNIT_TEST_FRAMEWORK_DYN_LINK itself -- which is how the CMake build
 * links). A build that instead provides its own main should compile this file
 * with -DMEQ_TEST_HAS_MAIN, which suppresses the definition here.
 */

#ifndef MEQ_TEST_HAS_MAIN
#ifndef BOOST_TEST_MODULE
#define BOOST_TEST_MODULE MeqConfigTests
#endif
#endif

#include <boost/test/unit_test.hpp>

#include "Config.hpp"

#include <filesystem>
#include <fstream>
#include <ostream>
#include <random>
#include <string>

namespace meq
{

	// Boost.Test insists on being able to print anything it compares, and a
	// scoped enumeration has no operator<< of its own. These live here rather
	// than in Config.hpp because they exist for the diagnostics of a failing
	// assertion and for nothing else.
	std::ostream & operator<<( std::ostream & os, SourceType type )
	{
		switch ( type )
		{
			case SourceType::Soloviev:     return os << "soloviev";
			case SourceType::MHD:          return os << "mhd";
			case SourceType::Manufactured: return os << "manufactured";
		}
		return os << "<invalid SourceType>";
	}

	std::ostream & operator<<( std::ostream & os, BoundaryDataType type )
	{
		switch ( type )
		{
			case BoundaryDataType::Zero:  return os << "zero";
			case BoundaryDataType::Exact: return os << "exact";
		}
		return os << "<invalid BoundaryDataType>";
	}

}

using meq::BoundaryDataType;
using meq::ConfigError;
using meq::Configuration;
using meq::SourceType;

namespace
{

	// A configuration with every required key and no optional one, so that a
	// test can perturb exactly the thing it is about.
	std::string minimal()
	{
		return
			"[mesh]\n"
			"RMin = 0.2\n"
			"RMax = 1.8\n"
			"ZMin = -1.6\n"
			"ZMax = 1.6\n"
			"\n"
			"[discretisation]\n"
			"PolynomialDegree = 2\n"
			"\n"
			"[source]\n"
			"Type = \"soloviev\"\n"
			"A = -0.52\n";
	}

	Configuration parse( std::string const & text )
	{
		return Configuration::fromString( text, "test.toml" );
	}

	// The directory holding the shipped examples. A build system that knows
	// where the source tree is can say so with -DMEQ_EXAMPLES_DIR="..."; failing
	// that, walk up from the working directory looking for it, so the tests run
	// from the source root, from build/, or from build/tests alike.
	std::filesystem::path examplesDirectory()
	{
#ifdef MEQ_EXAMPLES_DIR
		return std::filesystem::path( MEQ_EXAMPLES_DIR );
#else
		std::filesystem::path here = std::filesystem::current_path();
		for ( int level = 0; level < 6; ++level )
		{
			std::filesystem::path candidate = here / "examples" / "soloviev-nstx.toml";
			if ( std::filesystem::exists( candidate ) )
				return here / "examples";
			if ( !here.has_parent_path() || here.parent_path() == here )
				break;
			here = here.parent_path();
		}
		return std::filesystem::path( "examples" );
#endif
	}

	// Writes TOML to a real file and removes it again, for the paths that only
	// the file-taking constructor exercises.
	class TemporaryConfigFile
	{
		public:
			explicit TemporaryConfigFile( std::string const & contents )
				: path( std::filesystem::temp_directory_path() / uniqueName() )
			{
				std::ofstream out( path );
				out << contents;
			};

			~TemporaryConfigFile()
			{
				std::error_code ignored;
				std::filesystem::remove( path, ignored );
			};

			TemporaryConfigFile( TemporaryConfigFile const & ) = delete;
			TemporaryConfigFile & operator=( TemporaryConfigFile const & ) = delete;

			std::string name() const { return path.string(); };

		private:
			// Unique per process and per instance, so two test binaries
			// running at once cannot collide over the fixture.
			static std::string uniqueName()
			{
				static std::random_device entropy;
				static unsigned const run = entropy();
				static int counter = 0;
				return "meq-config-test-" + std::to_string( run ) + "-" + std::to_string( counter++ ) + ".toml";
			};

			std::filesystem::path path;
	};

	bool mentions( std::exception const & error, std::string const & needle )
	{
		return std::string( error.what() ).find( needle ) != std::string::npos;
	}

}

BOOST_AUTO_TEST_SUITE( config_tests, *boost::unit_test::tolerance( 1.0e-12 ) )

// ---------------------------------------------------------------------------
// The shipped examples
// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( soloviev_example_parses )
{
	std::string const file = ( examplesDirectory() / "soloviev-nstx.toml" ).string();
	BOOST_REQUIRE_MESSAGE( std::filesystem::exists( file ), "cannot find " << file );

	Configuration config( file );

	BOOST_TEST( config.getSource().type == SourceType::Soloviev );
	BOOST_TEST( config.getSource().getSoloviev().a == -0.52 );
	BOOST_TEST( config.getBoundary().type == BoundaryDataType::Zero );
	BOOST_TEST( config.getDiscretisation().tau == 1.0 );

	// The box has to contain the NSTX plasma boundary, 0.22 <= r <= 1.78,
	// |z| <= 1.56 in units of the major radius, or the benchmark is not the
	// benchmark.
	BOOST_TEST( config.getMesh().rMin < 0.22 );
	BOOST_TEST( config.getMesh().rMax > 1.78 );
	BOOST_TEST( config.getMesh().zMin < -1.56 );
	BOOST_TEST( config.getMesh().zMax > 1.56 );
	BOOST_TEST( config.getMesh().rMin > 0.0 );
}

BOOST_AUTO_TEST_CASE( manufactured_example_parses )
{
	std::string const file = ( examplesDirectory() / "manufactured.toml" ).string();
	BOOST_REQUIRE_MESSAGE( std::filesystem::exists( file ), "cannot find " << file );

	Configuration config( file );

	BOOST_TEST( config.getSource().type == SourceType::Manufactured );

	// Example 5 of refs/HDG-GradShafranov.pdf: r0 = -0.5, kr = 1.15 pi,
	// kz = 1.15.
	meq::ManufacturedParameters const & parameters = config.getSource().getManufactured();
	BOOST_TEST( parameters.r0 == -0.5 );
	BOOST_TEST( parameters.kr == 1.15*3.14159265358979323846 );
	BOOST_TEST( parameters.kz == 1.15 );

	// The manufactured psi is not zero on the boundary of the box.
	BOOST_TEST( config.getBoundary().type == BoundaryDataType::Exact );
}

// ---------------------------------------------------------------------------
// Integers where a number is expected
// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( integer_valued_keys_are_accepted_where_a_number_belongs )
{
	Configuration config = parse(
		"[mesh]\n"
		"RMin = 0\n"            // not 0.0
		"RMax = 2\n"
		"ZMin = -2\n"
		"ZMax = 2\n"
		"\n"
		"[discretisation]\n"
		"PolynomialDegree = 1\n"
		"Tau = 2\n"             // not 2.0
		"\n"
		"[source]\n"
		"Type = \"soloviev\"\n"
		"A = 1\n"               // not 1.0
		"\n"
		"[solver]\n"
		"NewtonRelativeTolerance = 1\n" );

	BOOST_TEST( config.getMesh().rMin == 0.0 );
	BOOST_TEST( config.getMesh().rMax == 2.0 );
	BOOST_TEST( config.getMesh().zMin == -2.0 );
	BOOST_TEST( config.getMesh().zMax == 2.0 );
	BOOST_TEST( config.getDiscretisation().tau == 2.0 );
	BOOST_TEST( config.getSource().getSoloviev().a == 1.0 );
	BOOST_TEST( config.getSolver().newtonRelativeTolerance == 1.0 );
}

BOOST_AUTO_TEST_CASE( floating_values_are_rejected_where_a_count_belongs )
{
	// The converse of the above: a count is a count. Truncating 4.7 cells
	// silently would be worse than saying so.
	std::string const fractionalCellCount =
		"[mesh]\n"
		"RMin = 0.2\n"
		"RMax = 1.8\n"
		"ZMin = -1.6\n"
		"ZMax = 1.6\n"
		"NR = 4.0\n"
		"\n"
		"[discretisation]\n"
		"PolynomialDegree = 2\n"
		"\n"
		"[source]\n"
		"Type = \"soloviev\"\n"
		"A = -0.52\n";

	BOOST_CHECK_EXCEPTION( parse( fractionalCellCount ), ConfigError,
		[]( ConfigError const & e ) { return e.getKey() == "mesh.NR" && mentions( e, "integer" ); } );
}

// ---------------------------------------------------------------------------
// Missing and mistyped keys
// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( missing_required_key_throws_and_names_the_key_and_file )
{
	std::string const text =
		"[mesh]\n"
		"RMin = 0.2\n"
		"ZMin = -1.6\n"
		"ZMax = 1.6\n"
		"\n"
		"[discretisation]\n"
		"PolynomialDegree = 2\n"
		"\n"
		"[source]\n"
		"Type = \"soloviev\"\n"
		"A = -0.52\n";

	BOOST_CHECK_EXCEPTION( parse( text ), ConfigError,
		[]( ConfigError const & e )
		{
			return e.getKey() == "mesh.RMax"
			    && e.getFile() == "test.toml"
			    && mentions( e, "RMax" )
			    && mentions( e, "test.toml" );
		} );
}

BOOST_AUTO_TEST_CASE( missing_required_table_throws_and_names_it )
{
	std::string const noSource =
		"[mesh]\n"
		"RMin = 0.2\n"
		"RMax = 1.8\n"
		"ZMin = -1.6\n"
		"ZMax = 1.6\n"
		"\n"
		"[discretisation]\n"
		"PolynomialDegree = 2\n";

	BOOST_CHECK_EXCEPTION( parse( noSource ), ConfigError,
		[]( ConfigError const & e ) { return e.getKey() == "source" && mentions( e, "[source]" ); } );
}

BOOST_AUTO_TEST_CASE( a_key_of_the_wrong_type_throws )
{
	BOOST_CHECK_EXCEPTION( parse( minimal() + "\n[output]\nPrefix = 3\n" ), ConfigError,
		[]( ConfigError const & e ) { return e.getKey() == "output.Prefix" && mentions( e, "string" ); } );

	std::string const stringForNumber =
		"[mesh]\n"
		"RMin = \"inner wall\"\n"
		"RMax = 1.8\n"
		"ZMin = -1.6\n"
		"ZMax = 1.6\n"
		"\n"
		"[discretisation]\n"
		"PolynomialDegree = 2\n"
		"\n"
		"[source]\n"
		"Type = \"soloviev\"\n"
		"A = -0.52\n";

	BOOST_CHECK_EXCEPTION( parse( stringForNumber ), ConfigError,
		[]( ConfigError const & e ) { return e.getKey() == "mesh.RMin" && mentions( e, "number" ); } );

	// A table where a scalar belongs, and a scalar where a table belongs.
	BOOST_CHECK_EXCEPTION( parse( minimal() + "\n[solver]\n[solver.NewtonMaxIterations]\nx = 1\n" ), ConfigError,
		[]( ConfigError const & e ) { return e.getKey() == "solver.NewtonMaxIterations"; } );
	BOOST_CHECK_EXCEPTION( parse( "mesh = 3\n" ), ConfigError,
		[]( ConfigError const & e ) { return e.getKey() == "mesh" && mentions( e, "must be a table" ); } );
}

BOOST_AUTO_TEST_CASE( unknown_key_is_rejected )
{
	// This is the failure the old example shipped with: a key nothing reads.
	BOOST_CHECK_EXCEPTION( parse( minimal() + "\n[output]\nPsiResolution = 0.01\n" ), ConfigError,
		[]( ConfigError const & e ) { return e.getKey() == "output.PsiResolution" && mentions( e, "accepted keys" ); } );

	BOOST_CHECK_EXCEPTION( parse( minimal() + "\n[coils]\nR = 1.0\n" ), ConfigError,
		[]( ConfigError const & e ) { return e.getKey() == "coils" && mentions( e, "schema" ); } );

	// A key belonging to a different source than the one selected.
	BOOST_CHECK_EXCEPTION( parse( minimal() + "Kr = 3.6\n" ), ConfigError,
		[]( ConfigError const & e ) { return e.getKey() == "source.Kr"; } );
}

BOOST_AUTO_TEST_CASE( a_misspelt_key_is_told_what_it_probably_meant )
{
	// A near miss gets a suggestion...
	BOOST_CHECK_EXCEPTION( parse( minimal() + "\n[output]\nPrefx = \"run\"\n" ), ConfigError,
		[]( ConfigError const & e ) { return mentions( e, "did you mean 'Prefix'?" ); } );

	BOOST_CHECK_EXCEPTION( parse( minimal() + "\n[misc]\nx = 1\n" ), ConfigError,
		[]( ConfigError const & e ) { return mentions( e, "did you mean [mesh]?" ); } );

	// ... and something unrelated does not get a misleading one.
	BOOST_CHECK_EXCEPTION( parse( minimal() + "\n[output]\nNetCDFFile = \"run.nc\"\n" ), ConfigError,
		[]( ConfigError const & e )
		{
			return e.getKey() == "output.NetCDFFile" && !mentions( e, "did you mean" );
		} );
}

// ---------------------------------------------------------------------------
// Defaults
// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( optional_keys_fall_back_to_their_documented_defaults )
{
	Configuration config = parse( minimal() );

	// tau = 1 is what both papers use; the old code used 5.0 with no
	// justification recorded anywhere.
	BOOST_TEST( config.getDiscretisation().tau == 1.0 );

	BOOST_TEST( config.getMesh().nR == 1 );
	BOOST_TEST( config.getMesh().nZ == 1 );
	BOOST_TEST( config.getMesh().refinementLevels == 0 );
	BOOST_TEST( config.getMesh().file == "" );
	BOOST_TEST( config.getMesh().fromFile() == false );

	BOOST_TEST( config.getBoundary().type == BoundaryDataType::Zero );

	BOOST_TEST( config.getSolver().newtonMaxIterations == 20 );
	BOOST_TEST( config.getSolver().newtonRelativeTolerance == 1.0e-8 );
	BOOST_TEST( config.getSolver().newtonAbsoluteTolerance == 1.0e-12 );
	BOOST_TEST( config.getSolver().linearMaxIterations == 1000 );
	BOOST_TEST( config.getSolver().linearTolerance == 1.0e-12 );

	BOOST_TEST( config.getOutput().directory == "." );
	BOOST_TEST( config.getOutput().prefix == "meq" );
	BOOST_TEST( config.getOutput().getMeshFile() == "./meq.mesh" );
	BOOST_TEST( config.getOutput().getPsiFile() == "./meq_psi.gf" );
	BOOST_TEST( config.getOutput().getGradPsiFile() == "./meq_grad_psi.gf" );
}

BOOST_AUTO_TEST_CASE( an_empty_optional_table_is_the_same_as_an_absent_one )
{
	Configuration config = parse( minimal() + "\n[solver]\n\n[output]\n\n[boundary]\n" );

	BOOST_TEST( config.getSolver().newtonMaxIterations == 20 );
	BOOST_TEST( config.getOutput().prefix == "meq" );
	BOOST_TEST( config.getBoundary().type == BoundaryDataType::Zero );
}

// ---------------------------------------------------------------------------
// Round trip
// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( values_round_trip_from_the_file_to_the_accessors )
{
	Configuration config = parse(
		"[mesh]\n"
		"RMin = 0.15\n"
		"RMax = 1.85\n"
		"ZMin = -1.45\n"
		"ZMax = 1.55\n"
		"NR = 7\n"
		"NZ = 11\n"
		"RefinementLevels = 4\n"
		"\n"
		"[discretisation]\n"
		"PolynomialDegree = 5\n"
		"Tau = 0.75\n"
		"\n"
		"[source]\n"
		"Type = \"manufactured\"\n"
		"R0 = -0.5\n"
		"Kr = 3.6128315516282616\n"
		"Kz = 1.15\n"
		"\n"
		"[boundary]\n"
		"Type = \"exact\"\n"
		"\n"
		"[solver]\n"
		"NewtonMaxIterations = 37\n"
		"NewtonRelativeTolerance = 1.0e-9\n"
		"NewtonAbsoluteTolerance = 2.5e-13\n"
		"LinearMaxIterations = 512\n"
		"LinearTolerance = 3.0e-11\n"
		"\n"
		"[output]\n"
		"Directory = \"/tmp/meq-run\"\n"
		"Prefix = \"case17\"\n" );

	BOOST_TEST( config.getFileName() == "test.toml" );

	BOOST_TEST( config.getMesh().rMin == 0.15 );
	BOOST_TEST( config.getMesh().rMax == 1.85 );
	BOOST_TEST( config.getMesh().zMin == -1.45 );
	BOOST_TEST( config.getMesh().zMax == 1.55 );
	BOOST_TEST( config.getMesh().nR == 7 );
	BOOST_TEST( config.getMesh().nZ == 11 );
	BOOST_TEST( config.getMesh().refinementLevels == 4 );

	BOOST_TEST( config.getDiscretisation().polynomialDegree == 5 );
	BOOST_TEST( config.getDiscretisation().tau == 0.75 );

	BOOST_TEST( config.getSource().type == SourceType::Manufactured );
	BOOST_TEST( config.getSource().getManufactured().r0 == -0.5 );
	BOOST_TEST( config.getSource().getManufactured().kr == 3.6128315516282616 );
	BOOST_TEST( config.getSource().getManufactured().kz == 1.15 );

	BOOST_TEST( config.getBoundary().type == BoundaryDataType::Exact );

	BOOST_TEST( config.getSolver().newtonMaxIterations == 37 );
	BOOST_TEST( config.getSolver().newtonRelativeTolerance == 1.0e-9 );
	BOOST_TEST( config.getSolver().newtonAbsoluteTolerance == 2.5e-13 );
	BOOST_TEST( config.getSolver().linearMaxIterations == 512 );
	BOOST_TEST( config.getSolver().linearTolerance == 3.0e-11 );

	BOOST_TEST( config.getOutput().directory == "/tmp/meq-run" );
	BOOST_TEST( config.getOutput().prefix == "case17" );
	BOOST_TEST( config.getOutput().getMeshFile() == "/tmp/meq-run/case17.mesh" );
	BOOST_TEST( config.getOutput().getPsiFile() == "/tmp/meq-run/case17_psi.gf" );
	BOOST_TEST( config.getOutput().getGradPsiFile() == "/tmp/meq-run/case17_grad_psi.gf" );
}

// ---------------------------------------------------------------------------
// Sources
// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( the_mhd_source_carries_its_profile_paths )
{
	Configuration config = parse(
		"[mesh]\n"
		"RMin = 4.0\n"
		"RMax = 8.4\n"
		"ZMin = -4.6\n"
		"ZMax = 4.6\n"
		"\n"
		"[discretisation]\n"
		"PolynomialDegree = 3\n"
		"\n"
		"[source]\n"
		"Type = \"mhd\"\n"
		"PPrimeFile = \"profiles/pprime.dat\"\n"
		"GGPrimeFile = \"profiles/ggprime.dat\"\n" );

	BOOST_TEST( config.getSource().type == SourceType::MHD );
	BOOST_TEST( config.getSource().getMHD().pPrimeFile == "profiles/pprime.dat" );
	BOOST_TEST( config.getSource().getMHD().ggPrimeFile == "profiles/ggprime.dat" );
	// SI by default; a normalised run sets Mu0 = 1.
	BOOST_TEST( config.getSource().getMHD().mu0 == 4.0e-7*3.14159265358979323846 );

	BOOST_CHECK_THROW( config.getSource().getSoloviev(), ConfigError );
	BOOST_CHECK_THROW( config.getSource().getManufactured(), ConfigError );
}

BOOST_AUTO_TEST_CASE( an_unknown_source_type_names_the_alternatives )
{
	std::string const text =
		"[mesh]\n"
		"RMin = 0.2\n"
		"RMax = 1.8\n"
		"ZMin = -1.6\n"
		"ZMax = 1.6\n"
		"\n"
		"[discretisation]\n"
		"PolynomialDegree = 2\n"
		"\n"
		"[source]\n"
		"Type = \"solovev\"\n";

	BOOST_CHECK_EXCEPTION( parse( text ), ConfigError,
		[]( ConfigError const & e )
		{
			return e.getKey() == "source.Type" && mentions( e, "soloviev" ) && mentions( e, "manufactured" );
		} );
}

BOOST_AUTO_TEST_CASE( exact_boundary_data_needs_a_source_that_has_an_exact_solution )
{
	std::string const text =
		"[mesh]\n"
		"RMin = 4.0\n"
		"RMax = 8.4\n"
		"ZMin = -4.6\n"
		"ZMax = 4.6\n"
		"\n"
		"[discretisation]\n"
		"PolynomialDegree = 3\n"
		"\n"
		"[source]\n"
		"Type = \"mhd\"\n"
		"PPrimeFile = \"p.dat\"\n"
		"GGPrimeFile = \"gg.dat\"\n"
		"\n"
		"[boundary]\n"
		"Type = \"exact\"\n";

	BOOST_CHECK_EXCEPTION( parse( text ), ConfigError,
		[]( ConfigError const & e ) { return e.getKey() == "boundary.Type" && mentions( e, "mhd" ); } );
}

// ---------------------------------------------------------------------------
// Values that parse but cannot describe a run
// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( an_inverted_or_negative_box_is_rejected )
{
	std::string const inverted =
		"[mesh]\n"
		"RMin = 1.8\n"
		"RMax = 0.2\n"
		"ZMin = -1.6\n"
		"ZMax = 1.6\n"
		"\n"
		"[discretisation]\n"
		"PolynomialDegree = 2\n"
		"\n"
		"[source]\n"
		"Type = \"soloviev\"\n"
		"A = -0.52\n";

	BOOST_CHECK_EXCEPTION( parse( inverted ), ConfigError,
		[]( ConfigError const & e ) { return e.getKey() == "mesh.RMax" && mentions( e, "RMin" ); } );

	std::string const negativeRadius =
		"[mesh]\n"
		"RMin = -0.2\n"
		"RMax = 1.8\n"
		"ZMin = -1.6\n"
		"ZMax = 1.6\n"
		"\n"
		"[discretisation]\n"
		"PolynomialDegree = 2\n"
		"\n"
		"[source]\n"
		"Type = \"soloviev\"\n"
		"A = -0.52\n";

	BOOST_CHECK_EXCEPTION( parse( negativeRadius ), ConfigError,
		[]( ConfigError const & e ) { return e.getKey() == "mesh.RMin"; } );
}

BOOST_AUTO_TEST_CASE( out_of_range_scalars_are_rejected )
{
	std::string const zeroTau =
		"[mesh]\n"
		"RMin = 0.2\n"
		"RMax = 1.8\n"
		"ZMin = -1.6\n"
		"ZMax = 1.6\n"
		"\n"
		"[discretisation]\n"
		"PolynomialDegree = 2\n"
		"Tau = 0.0\n"
		"\n"
		"[source]\n"
		"Type = \"soloviev\"\n"
		"A = -0.52\n";

	BOOST_CHECK_EXCEPTION( parse( zeroTau ), ConfigError,
		[]( ConfigError const & e ) { return e.getKey() == "discretisation.Tau"; } );

	BOOST_CHECK_EXCEPTION( parse( minimal() + "\n[solver]\nNewtonMaxIterations = 0\n" ), ConfigError,
		[]( ConfigError const & e ) { return e.getKey() == "solver.NewtonMaxIterations"; } );

	BOOST_CHECK_EXCEPTION( parse( minimal() + "\n[solver]\nLinearTolerance = -1.0e-8\n" ), ConfigError,
		[]( ConfigError const & e ) { return e.getKey() == "solver.LinearTolerance"; } );

	BOOST_CHECK_EXCEPTION( parse( minimal() + "\n[output]\nPrefix = \"\"\n" ), ConfigError,
		[]( ConfigError const & e ) { return e.getKey() == "output.Prefix"; } );
}

// ---------------------------------------------------------------------------
// A mesh read from a file instead of generated
// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( a_mesh_file_replaces_the_box )
{
	Configuration config = parse(
		"[mesh]\n"
		"File = \"meshes/iter.msh\"\n"
		"RefinementLevels = 2\n"
		"\n"
		"[discretisation]\n"
		"PolynomialDegree = 3\n"
		"\n"
		"[source]\n"
		"Type = \"soloviev\"\n"
		"A = -0.52\n" );

	BOOST_TEST( config.getMesh().fromFile() == true );
	BOOST_TEST( config.getMesh().file == "meshes/iter.msh" );
	BOOST_TEST( config.getMesh().refinementLevels == 2 );
}

// ---------------------------------------------------------------------------
// Broken input
// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( a_syntax_error_is_reported_as_one )
{
	BOOST_CHECK_EXCEPTION( parse( "[mesh\nRMin = 0.2\n" ), ConfigError,
		[]( ConfigError const & e )
		{
			return e.getKey().empty() && mentions( e, "syntax" ) && mentions( e, "test.toml" );
		} );

	BOOST_CHECK_EXCEPTION( parse( "[mesh]\nRMin = = 0.2\n" ), ConfigError,
		[]( ConfigError const & e ) { return mentions( e, "syntax" ); } );

	// A table defined twice is a TOML error, not a meq one, and arrives here
	// through the same path.
	BOOST_CHECK_EXCEPTION( parse( minimal() + "\n[discretisation]\n" ), ConfigError,
		[]( ConfigError const & e ) { return mentions( e, "syntax" ); } );

	// Also from a file, which is a different code path.
	TemporaryConfigFile broken( "[mesh]\nRMin 0.2\n" );
	BOOST_CHECK_EXCEPTION( Configuration( broken.name() ), ConfigError,
		[ &broken ]( ConfigError const & e ) { return mentions( e, "syntax" ) && e.getFile() == broken.name(); } );
}

BOOST_AUTO_TEST_CASE( an_unreadable_file_is_reported_with_its_name )
{
	std::string const missing = ( std::filesystem::temp_directory_path() / "meq-no-such-config.toml" ).string();
	std::error_code ignored;
	std::filesystem::remove( missing, ignored );

	BOOST_CHECK_EXCEPTION( Configuration{ missing }, ConfigError,
		[ &missing ]( ConfigError const & e ) { return e.getFile() == missing && mentions( e, missing ); } );
}

BOOST_AUTO_TEST_CASE( a_good_file_on_disk_reads_the_same_as_the_same_text_in_memory )
{
	TemporaryConfigFile good( minimal() );

	Configuration fromFile( good.name() );
	Configuration fromString = parse( minimal() );

	BOOST_TEST( fromFile.getFileName() == good.name() );
	BOOST_TEST( fromFile.getMesh().rMin == fromString.getMesh().rMin );
	BOOST_TEST( fromFile.getMesh().rMax == fromString.getMesh().rMax );
	BOOST_TEST( fromFile.getDiscretisation().polynomialDegree == fromString.getDiscretisation().polynomialDegree );
	BOOST_TEST( fromFile.getSource().getSoloviev().a == fromString.getSource().getSoloviev().a );
}

BOOST_AUTO_TEST_SUITE_END()
