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
			case SourceType::Rotating:     return os << "rotating";
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

	// The two performance keys. Spelled as the TOML spells them, so a failing
	// assertion prints the value a reader would put in the file.
	std::ostream & operator<<( std::ostream & os, AssemblyModeType mode )
	{
		switch ( mode )
		{
			case AssemblyModeType::Serial:   return os << "serial";
			case AssemblyModeType::Threaded: return os << "threaded";
		}
		return os << "<invalid AssemblyModeType>";
	}

	std::ostream & operator<<( std::ostream & os, TraceSolverType solver )
	{
		switch ( solver )
		{
			case TraceSolverType::UMFPack: return os << "umfpack";
			case TraceSolverType::Pardiso: return os << "pardiso";
			case TraceSolverType::cuDSS:   return os << "cudss";
		}
		return os << "<invalid TraceSolverType>";
	}

}

using meq::AssemblyModeType;
using meq::BoundaryDataType;
using meq::ConfigError;
using meq::Configuration;
using meq::SourceType;
using meq::TraceSolverType;

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

	// The mesh and discretisation of minimal(), with the [source] table left to
	// the caller. A rotating source cannot be appended to minimal(), whose
	// [source] is already a Solov'ev one -- and [[source.species]] has to follow
	// the [source] scalars it belongs to, so the whole table is built at once.
	std::string withSource( std::string const & source )
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
			"\n" + source;
	}

	// A well-formed two-species rotating source, so that a refusal test can
	// perturb exactly the one thing it is about. `head` is appended to the
	// [source] scalars, `ion` and `electron` to their species tables; a key
	// given in one of those overrides nothing, so a test that wants a key
	// *changed* rather than added spells the whole source out itself.
	std::string rotating( std::string const & head = "",
	                      std::string const & ion = "",
	                      std::string const & electron = "" )
	{
		return withSource(
			"[source]\n"
			"Type = \"rotating\"\n"
			"GGPrime = 0.5\n"
			+ head +
			"\n[[source.species]]\n"
			"Name = \"deuterium\"\n"
			"Mass = 3.344e-27\n"
			"Charge = 1.0\n"
			"Temperature = 1.6e-16\n"
			"Density = 5.0e19\n"
			+ ion +
			"\n[[source.species]]\n"
			"Name = \"electrons\"\n"
			"Mass = 9.109e-31\n"
			"Charge = -1.0\n"
			"Temperature = 1.6e-16\n"
			"Neutralising = true\n"
			+ electron );
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

	// A TABLE THAT IS NOT PART OF THE SCHEMA. This used to be spelled
	// [coils], and [[coils]] is part of the schema as of FB-6 -- so the test
	// started failing the day it was added, which is the suite doing its job.
	// The example is now a name nothing is likely to claim.
	BOOST_CHECK_EXCEPTION( parse( minimal() + "\n[windings]\nR = 1.0\n" ), ConfigError,
		[]( ConfigError const & e ) { return e.getKey() == "windings" && mentions( e, "schema" ); } );

	// And [coils] written as a plain table rather than an array of tables is
	// caught too -- it IS in the schema now, so it reaches the array check
	// rather than the table-name one.
	BOOST_CHECK_EXCEPTION( parse( minimal() + "\n[coils]\nR = 1.0\n" ), ConfigError,
		[]( ConfigError const & e ) { return e.getKey() == "coils"; } );

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

	BOOST_TEST( config.getOutput().directory == "." );
	BOOST_TEST( config.getOutput().prefix == "meq" );
	BOOST_TEST( config.getOutput().getMeshFile() == "./meq.mesh" );
	BOOST_TEST( config.getOutput().getPsiFile() == "./meq_psi.gf" );
	BOOST_TEST( config.getOutput().getGradPsiFile() == "./meq_grad_psi.gf" );
	BOOST_TEST( config.getOutput().getPsiStarFile() == "./meq_psistar.gf" );
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

	BOOST_TEST( config.getOutput().directory == "/tmp/meq-run" );
	BOOST_TEST( config.getOutput().prefix == "case17" );
	BOOST_TEST( config.getOutput().getMeshFile() == "/tmp/meq-run/case17.mesh" );
	BOOST_TEST( config.getOutput().getPsiFile() == "/tmp/meq-run/case17_psi.gf" );
	BOOST_TEST( config.getOutput().getGradPsiFile() == "/tmp/meq-run/case17_grad_psi.gf" );
	BOOST_TEST( config.getOutput().getPsiStarFile() == "/tmp/meq-run/case17_psistar.gf" );
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
			return e.getKey() == "source.Type" && mentions( e, "soloviev" )
			    && mentions( e, "manufactured" ) && mentions( e, "rotating" );
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

	// REFUSED, NOT VALIDATED, AND THE DIFFERENCE IS THE POINT. These once
	// parsed, validated and were read by nothing -- so a file could set them,
	// pass every check, and change no behaviour whatever. They describe an
	// iterative inner solve and MEQ's trace solve is direct. A VALID value is
	// therefore refused just as an invalid one is, which is what these two
	// assert: 512 is a perfectly reasonable iteration cap and is still an error.
	BOOST_CHECK_EXCEPTION( parse( minimal() + "\n[solver]\nLinearMaxIterations = 512\n" ), ConfigError,
		[]( ConfigError const & e ) { return e.getKey() == "solver.LinearMaxIterations"; } );
	BOOST_CHECK_EXCEPTION( parse( minimal() + "\n[solver]\nLinearTolerance = 3.0e-11\n" ), ConfigError,
		[]( ConfigError const & e ) { return e.getKey() == "solver.LinearTolerance"; } );

	BOOST_CHECK_EXCEPTION( parse( minimal() + "\n[output]\nPrefix = \"\"\n" ), ConfigError,
		[]( ConfigError const & e ) { return e.getKey() == "output.Prefix"; } );

	// The two performance keys. A wrong SPELLING is a parse fault; an
	// unavailable CHOICE is not, and is refused by apps/meq.cpp instead -- this
	// header is one of the four CI compiles without MFEM, so it cannot know
	// whether this build has OpenMP or PARDISO. Two faults, two messages.
	//
	// The capitalised spelling is the one worth asserting: there is no
	// case-folding anywhere in Config, so "Threaded" must fail rather than be
	// quietly accepted, and a reader who assumes TOML VALUES follow the
	// UpperCamelCase rule that TOML KEYS follow would write exactly that.
	BOOST_CHECK_EXCEPTION( parse( minimal() + "\n[solver]\nAssemblyMode = \"parallel\"\n" ), ConfigError,
		[]( ConfigError const & e ) { return e.getKey() == "solver.AssemblyMode"; } );
	BOOST_CHECK_EXCEPTION( parse( minimal() + "\n[solver]\nAssemblyMode = \"Threaded\"\n" ), ConfigError,
		[]( ConfigError const & e ) { return e.getKey() == "solver.AssemblyMode"; } );
	BOOST_CHECK_EXCEPTION( parse( minimal() + "\n[solver]\nTraceSolver = \"mumps\"\n" ), ConfigError,
		[]( ConfigError const & e ) { return e.getKey() == "solver.TraceSolver"; } );
	BOOST_CHECK_EXCEPTION( parse( minimal() + "\n[solver]\nTraceSolver = \"cuDSS\"\n" ), ConfigError,
		[]( ConfigError const & e ) { return e.getKey() == "solver.TraceSolver"; } );
}

/*
 * The two performance keys parse, and their defaults are the measured ones.
 *
 * Threaded and PARDISO are not arbitrary defaults and the test says so, because
 * the temptation to "improve" them in a config file is exactly what this pair of
 * keys creates. Threaded became the default when MFEM extended the flag to
 * cover the residual and the Jacobian as well as the assembly, and the test
 * that had most argued against it -- HighBetaConvergence, 1.8x SLOWER under an
 * automatic gate -- came back 2.4x faster. See setAssemblyMode()'s own
 * documentation for the table.
 *
 * PARDISO IS THE DEFAULT BECAUSE OF THE OTHER KEY, WHICH IS WHY THE PAIR IS
 * TESTED TOGETHER. It is faster than UMFPack single-threaded and it is the only
 * one of the two whose MKL threads are spendable under threaded assembly: MKL
 * suppresses its own threading inside an active OpenMP region, so the
 * element-local dense work is nested and free while the trace solve runs on the
 * master thread outside every region and takes them all. UMFPack's BLAS sits
 * outside any parallel region and pays full MKL threading per frontal matrix,
 * so the two respond to MKL_NUM_THREADS in OPPOSITE directions.
 *
 * THE DEFAULT HERE IS UNCONDITIONAL AND THE LIBRARY'S IS NOT, and that
 * asymmetry is the MFEM-free rule rather than an inconsistency. Config cannot
 * ask whether this build has oneMKL -- it is one of the four translation units
 * CI compiles without MFEM at all -- so it names PARDISO always and
 * apps/meq.cpp downgrades an INHERITED one to UMFPack where the build has no
 * PARDISO, refusing only one the file asked for. traceSolverWasGiven is what
 * carries that distinction, and it is the reason the flag exists.
 *
 * Neither key may change the ANSWER, which is what makes exposing them safe at
 * all: the two assembly modes are bit-identical at MKL_NUM_THREADS=1 and the
 * three trace solvers agree to 1e-14. That property is asserted in
 * tests/convergence/, on both a linear and a nonlinear source; here it is only
 * the parsing.
 */
BOOST_AUTO_TEST_CASE( the_two_performance_keys_parse_and_default_to_the_measured_choice )
{
	Configuration const bare = parse( minimal() );
	BOOST_TEST( bare.getSolver().assemblyMode == AssemblyModeType::Threaded );
	BOOST_TEST( bare.getSolver().traceSolver == TraceSolverType::Pardiso );

	// AND NEITHER WAS ASKED FOR, which the driver has to know: an inherited
	// "threaded" is downgraded to Serial on a build without OpenMP and an
	// inherited "pardiso" to UMFPack on a build without oneMKL, while one the
	// file states is refused. Without that distinction the defaults would make
	// every example in this tree fail on a stock MFEM.
	BOOST_TEST( bare.getSolver().assemblyModeWasGiven == false );
	BOOST_TEST( bare.getSolver().traceSolverWasGiven == false );

	Configuration const serial = parse( minimal()
		+ "\n[solver]\nAssemblyMode = \"serial\"\nTraceSolver = \"umfpack\"\n" );
	BOOST_TEST( serial.getSolver().assemblyMode == AssemblyModeType::Serial );
	BOOST_TEST( serial.getSolver().traceSolver == TraceSolverType::UMFPack );
	BOOST_TEST( serial.getSolver().assemblyModeWasGiven == true );
	BOOST_TEST( serial.getSolver().traceSolverWasGiven == true );

	Configuration const device = parse( minimal()
		+ "\n[solver]\nTraceSolver = \"cudss\"\n" );
	BOOST_TEST( device.getSolver().traceSolver == TraceSolverType::cuDSS );
	BOOST_TEST( device.getSolver().traceSolverWasGiven == true );

	// And the explicit spellings of the defaults, so that writing them down
	// cannot mean something different from leaving them out -- except in the
	// two WasGiven flags, which is the whole point of them.
	Configuration const explicitDefaults = parse( minimal()
		+ "\n[solver]\nAssemblyMode = \"threaded\"\nTraceSolver = \"pardiso\"\n" );
	BOOST_TEST( explicitDefaults.getSolver().assemblyMode == AssemblyModeType::Threaded );
	BOOST_TEST( explicitDefaults.getSolver().traceSolver == TraceSolverType::Pardiso );
	BOOST_TEST( explicitDefaults.getSolver().assemblyModeWasGiven == true );
	BOOST_TEST( explicitDefaults.getSolver().traceSolverWasGiven == true );
}

// ---------------------------------------------------------------------------
// The keys reserved for the unwritten NetCDF profile facility
// ---------------------------------------------------------------------------

/// A RESERVED KEY MUST BE REFUSED ON EVERY PATH THAT ACCEPTS IT, AND THIS IS
/// THE TEST THAT WAS MISSING WHILE ONE OF THEM WAS NOT.
///
/// ProfileFile, and <Profile>Variable / <Profile>Fit beside it, name a facility
/// that is not implemented: several profiles read out of one NetCDF file. They
/// are listed as accepted keys so that a reader who tries them is told that,
/// rather than told the key does not exist -- which means the refusal is the
/// only thing standing between the key and being silently ignored.
///
/// It was not standing on the "mhd" path. That source requires PPrimeFile and
/// GGPrimeFile and offers no constant form, so it reads them directly instead
/// of through readEitherProfile, and never reached the refusal the rotating
/// path got for free. GGPrimeVariable was therefore a hard error under
/// Type = "rotating" and inert under Type = "mhd" -- the same key name, two
/// behaviours, decided by a neighbouring key. Nothing tested either path, which
/// is how it survived.
///
/// So the assertions below are paired deliberately: the interesting property is
/// not that a key is refused but that BOTH source types refuse it, which is
/// what a single-path test would have passed without establishing.
BOOST_AUTO_TEST_CASE( the_reserved_profile_keys_are_refused_on_every_source )
{
	auto const mhd = []( std::string const & extra )
	{
		return withSource(
			"[source]\n"
			"Type = \"mhd\"\n"
			"PPrimeFile = \"profiles/pprime.dat\"\n"
			"GGPrimeFile = \"profiles/ggprime.dat\"\n" + extra );
	};

	auto const refuses = []( std::string const & text, std::string const & key )
	{
		BOOST_CHECK_EXCEPTION( parse( text ), ConfigError,
			[&]( ConfigError const & e ) { return e.getKey() == key; } );
	};

	// The four on the "mhd" path, which is the one that used to ignore them.
	refuses( mhd( "PPrimeVariable = \"Var0\"\n" ), "source.PPrimeVariable" );
	refuses( mhd( "PPrimeFit = \"pchip\"\n" ), "source.PPrimeFit" );
	refuses( mhd( "GGPrimeVariable = \"Var1\"\n" ), "source.GGPrimeVariable" );
	refuses( mhd( "GGPrimeFit = \"pchip\"\n" ), "source.GGPrimeFit" );
	refuses( mhd( "ProfileFile = \"equilibrium.nc\"\n" ), "source.ProfileFile" );

	// The same names on the rotating path, plus the two it has of its own.
	refuses( rotating( "GGPrimeVariable = \"Var1\"\n" ), "source.GGPrimeVariable" );
	refuses( rotating( "GGPrimeFit = \"pchip\"\n" ), "source.GGPrimeFit" );
	refuses( rotating( "OmegaVariable = \"Var2\"\n" ), "source.OmegaVariable" );
	refuses( rotating( "OmegaFit = \"pchip\"\n" ), "source.OmegaFit" );
	refuses( rotating( "ProfileFile = \"equilibrium.nc\"\n" ), "source.ProfileFile" );

	// And per species, where the same helper is reached through Temperature and
	// Density rather than through a [source] scalar.
	refuses( rotating( "", "TemperatureVariable = \"Var3\"\n" ), "source.species[0].TemperatureVariable" );
	refuses( rotating( "", "DensityFit = \"pchip\"\n" ), "source.species[0].DensityFit" );

	// THE MESSAGE IS THE POINT OF RESERVING THEM, so it is asserted rather than
	// left to the key name: a reader who meets one has to learn that the
	// facility is unwritten, not merely that this file will not start.
	try
	{
		parse( mhd( "PPrimeVariable = \"Var0\"\n" ) );
		BOOST_FAIL( "PPrimeVariable was accepted" );
	}
	catch ( ConfigError const & error )
	{
		std::string const text = error.what();
		BOOST_TEST( text.find( "not implemented yet" ) != std::string::npos );
		// The advice names PPrimeFile and NOT a bare PPrime, because the "mhd"
		// source has no constant form. The rotating path has both and says so.
		BOOST_TEST( text.find( "PPrimeFile" ) != std::string::npos );
		BOOST_TEST( text.find( "give PPrime as a constant" ) == std::string::npos );
	}

	try
	{
		parse( rotating( "GGPrimeVariable = \"Var1\"\n" ) );
		BOOST_FAIL( "GGPrimeVariable was accepted" );
	}
	catch ( ConfigError const & error )
	{
		std::string const text = error.what();
		BOOST_TEST( text.find( "not implemented yet" ) != std::string::npos );
		BOOST_TEST( text.find( "give GGPrime as a constant" ) != std::string::npos );
	}
}

/// THE BANNER IS THE CALLER'S, NOT THE MESSAGE'S. apps/meq.cpp prints every
/// line it writes as "MEQ: <text>" and catches std::exception, so a banner
/// inside what() gave "MEQ: MEQ: configuration error in ...". It went unread
/// while the project spelled itself in lower case and became obvious when it
/// did not. Asserted here rather than left to a reader's eye, since nothing
/// else in the suite looks at the whole string.
BOOST_AUTO_TEST_CASE( a_configuration_error_carries_no_banner_of_its_own )
{
	// An unknown key, chosen because its message is the one place in the file
	// where MEQ does not name itself in the body -- the refusals that explain a
	// design decision generally do, so they cannot tell a banner from prose.
	try
	{
		parse( minimal() + "\n[output]\nNoSuchKey = 1\n" );
		BOOST_FAIL( "NoSuchKey was accepted" );
	}
	catch ( ConfigError const & error )
	{
		std::string const text = error.what();
		BOOST_TEST( text.rfind( "configuration error", 0 ) == 0u );
		BOOST_TEST( text.find( "MEQ" ) == std::string::npos );
	}
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

/*
 * [boundary.shape] -- the curved Gamma of stage 5, as configuration.
 *
 * Config only reports what the TOML said; meq::BoundaryShape is what decides
 * whether the curve is usable, because star-shapedness is a property of the
 * whole curve and not of any one key. So these check the PARSE, and
 * BoundaryShapeTests.cpp checks the geometry.
 */
BOOST_AUTO_TEST_CASE( a_miller_shape_is_read )
{
	meq::Configuration const config = parse( minimal() +
		"\n[boundary]\nType = \"zero\"\n"
		"[boundary.shape]\nType = \"miller\"\n"
		"R0 = 1.0\nMinorRadius = 0.32\nElongation = 1.7\nTriangularity = 0.33\n" );

	meq::ShapeConfig const &shape = config.getBoundary().shape;
	BOOST_TEST( ( shape.type == meq::ShapeType::Miller ) );
	BOOST_TEST( shape.majorRadius == 1.0, boost::test_tools::tolerance( 1.0e-12 ) );
	BOOST_TEST( shape.minorRadius == 0.32, boost::test_tools::tolerance( 1.0e-12 ) );
	BOOST_TEST( shape.elongation == 1.7, boost::test_tools::tolerance( 1.0e-12 ) );
	BOOST_TEST( shape.triangularity == 0.33, boost::test_tools::tolerance( 1.0e-12 ) );
	BOOST_TEST( shape.centreHeight == 0.0, boost::test_tools::tolerance( 1.0e-12 ) );
}

BOOST_AUTO_TEST_CASE( an_mxh_shape_reads_its_harmonics )
{
	meq::Configuration const config = parse( minimal() +
		"\n[boundary.shape]\nType = \"mxh\"\n"
		"R0 = 1.7\nZ0 = 0.1\nMinorRadius = 0.6\nElongation = 1.8\n"
		"CosCoefficients = [ -0.1, 0.02 ]\nSinCoefficients = [ 0.35, 0.08 ]\n" );

	meq::ShapeConfig const &shape = config.getBoundary().shape;
	BOOST_TEST( ( shape.type == meq::ShapeType::Mxh ) );
	BOOST_TEST( shape.centreHeight == 0.1, boost::test_tools::tolerance( 1.0e-12 ) );
	BOOST_TEST_REQUIRE( shape.cosCoefficients.size() == 2u );
	BOOST_TEST_REQUIRE( shape.sinCoefficients.size() == 2u );
	BOOST_TEST( shape.cosCoefficients[ 0 ] == -0.1, boost::test_tools::tolerance( 1.0e-12 ) );
	BOOST_TEST( shape.sinCoefficients[ 1 ] == 0.08, boost::test_tools::tolerance( 1.0e-12 ) );
}

/// Integers are numbers. The same case asFloat() exists for, now reached
/// through an array: Elongation = 2 and CosCoefficients = [ 0, 1 ] are both
/// legitimate TOML and both used to be a silent default.
BOOST_AUTO_TEST_CASE( integer_valued_shape_keys_are_numbers )
{
	meq::Configuration const config = parse( minimal() +
		"\n[boundary.shape]\nType = \"mxh\"\n"
		"R0 = 2\nMinorRadius = 1\nElongation = 2\n"
		"SinCoefficients = [ 0, 1 ]\n" );

	meq::ShapeConfig const &shape = config.getBoundary().shape;
	BOOST_TEST( shape.majorRadius == 2.0, boost::test_tools::tolerance( 1.0e-12 ) );
	BOOST_TEST( shape.elongation == 2.0, boost::test_tools::tolerance( 1.0e-12 ) );
	BOOST_TEST_REQUIRE( shape.sinCoefficients.size() == 2u );
	BOOST_TEST( shape.sinCoefficients[ 1 ] == 1.0, boost::test_tools::tolerance( 1.0e-12 ) );
}

/// No shape is the default and is not an error: that is the fitted path, which
/// is what every convergence test in the suite uses.
BOOST_AUTO_TEST_CASE( no_shape_section_means_the_fitted_path )
{
	meq::Configuration const config = parse( minimal() + "\n[boundary]\nType = \"zero\"\n" );
	BOOST_TEST( ( config.getBoundary().shape.type == meq::ShapeType::None ) );
}

/// The keys of the other shape are REFUSED rather than ignored. A Triangularity
/// sitting unread under Type = "mxh" is a file that says one thing and does
/// another.
BOOST_AUTO_TEST_CASE( the_wrong_shapes_keys_are_refused )
{
	BOOST_CHECK_THROW( parse( minimal() +
		"\n[boundary.shape]\nType = \"mxh\"\nR0 = 1.0\nMinorRadius = 0.32\n"
		"Triangularity = 0.33\nSinCoefficients = [ 0.3 ]\n" ), ConfigError );

	BOOST_CHECK_THROW( parse( minimal() +
		"\n[boundary.shape]\nType = \"miller\"\nR0 = 1.0\nMinorRadius = 0.32\n"
		"SinCoefficients = [ 0.3 ]\n" ), ConfigError );
}

BOOST_AUTO_TEST_CASE( a_shape_missing_its_geometry_is_refused )
{
	// MinorRadius has no sensible default, so it is required rather than
	// silently zero -- which would be a degenerate curve the geometry layer
	// would then reject with a less helpful message.
	BOOST_CHECK_THROW( parse( minimal() +
		"\n[boundary.shape]\nType = \"miller\"\nR0 = 1.0\n" ), ConfigError );

	BOOST_CHECK_THROW( parse( minimal() +
		"\n[boundary.shape]\nType = \"lozenge\"\n" ), ConfigError );
}


/*
 * [source] Type = "rotating" -- sonic toroidal rotation, as configuration.
 * docs/rotation.rst is the physics.
 *
 * Two things here are new to meq's schema and are what these cases are mostly
 * about. [[source.species]] is the first ARRAY OF TABLES the parser reads, so
 * an element's diagnostics have to stay qualified all the way down to
 * source.species[2].Mass rather than collapsing to source.species; and a
 * profile is now given EITHER as a constant in Temperature OR as a path in
 * TemperatureFile, with TemperatureScale multiplying whichever it was.
 *
 * As everywhere else in this file these check the PARSE. Whether the species
 * set describes a plasma -- charge neutrality on the reference curve, in
 * particular -- is meq::RotatingSource's to decide, and SourceFactoryTests is
 * where that is asserted.
 */

BOOST_AUTO_TEST_CASE( a_rotating_source_carries_its_species )
{
	Configuration config = parse( withSource(
		"[source]\n"
		"Type = \"rotating\"\n"
		"Omega = 1.2e5\n"
		"OmegaScale = 2.0\n"
		"GGPrimeFile = \"profiles/ggprime.dat\"\n"
		"GGPrimeScale = 0.5\n"
		"ReferenceRadius = 1.65\n"
		"Mu0 = 1.0\n"
		"\n[[source.species]]\n"
		"Name = \"deuterium\"\n"
		"Mass = 3.344e-27\n"
		"Charge = 1.0\n"
		"TemperatureFile = \"profiles/ti.dat\"\n"
		"TemperatureScale = 1.602176634e-16\n"
		"Density = 5.0e19\n"
		"DensityScale = 2.0\n"
		"\n[[source.species]]\n"
		"Name = \"electrons\"\n"
		"Mass = 9.1093837015e-31\n"
		"Charge = -1.0\n"
		"Temperature = 3.2e-16\n"
		"Neutralising = true\n" ) );

	BOOST_TEST( config.getSource().type == SourceType::Rotating );

	meq::RotatingParameters const & rotating = config.getSource().getRotating();

	// The [source] scalars. Omega given as a constant, so omegaFile is empty
	// and omegaGiven says the rotation was asked for at all.
	BOOST_TEST( rotating.omega == 1.2e5 );
	BOOST_TEST( rotating.omegaFile == "" );
	BOOST_TEST( rotating.omegaScale == 2.0 );
	BOOST_TEST( rotating.omegaGiven == true );
	// ... and g g' the other way round: a path, and no constant.
	BOOST_TEST( rotating.ggPrime == 0.0 );
	BOOST_TEST( rotating.ggPrimeFile == "profiles/ggprime.dat" );
	BOOST_TEST( rotating.ggPrimeScale == 0.5 );
	BOOST_TEST( rotating.referenceRadius == 1.65 );
	BOOST_TEST( rotating.mu0 == 1.0 );

	// The species, in file order.
	BOOST_TEST_REQUIRE( rotating.species.size() == 2u );

	meq::SpeciesParameters const & ion = rotating.species[ 0 ];
	BOOST_TEST( ion.name == "deuterium" );
	BOOST_TEST( ion.mass == 3.344e-27 );
	BOOST_TEST( ion.charge == 1.0 );
	BOOST_TEST( ion.temperature == 0.0 );
	BOOST_TEST( ion.temperatureFile == "profiles/ti.dat" );
	BOOST_TEST( ion.temperatureScale == 1.602176634e-16 );
	BOOST_TEST( ion.density == 5.0e19 );
	BOOST_TEST( ion.densityFile == "" );
	BOOST_TEST( ion.densityScale == 2.0 );
	BOOST_TEST( ion.neutralising == false );

	meq::SpeciesParameters const & electrons = rotating.species[ 1 ];
	BOOST_TEST( electrons.name == "electrons" );
	BOOST_TEST( electrons.mass == 9.1093837015e-31 );
	BOOST_TEST( electrons.charge == -1.0 );
	BOOST_TEST( electrons.temperature == 3.2e-16 );
	BOOST_TEST( electrons.temperatureFile == "" );
	BOOST_TEST( electrons.temperatureScale == 1.0 );
	BOOST_TEST( electrons.neutralising == true );
	// Neither form of the density, which is what Neutralising means: it is
	// derived from the others rather than given.
	BOOST_TEST( electrons.density == 0.0 );
	BOOST_TEST( electrons.densityFile == "" );

	// Not normalised unless it says so, so psi_ax is not an unknown.
	BOOST_TEST( config.getSource().isNormalised() == false );
	BOOST_TEST( config.getSource().psiAxisGuess() == 0.0 );

	BOOST_CHECK_THROW( config.getSource().getMHD(), ConfigError );
}

/// BOTH OMEGA KEYS ABSENT MEANS NO ROTATION, which is a documented state and
/// not an omission: the source then reduces to the static equation
/// meq::MHDSource solves, and being able to ask for that from the same file is
/// how the omega -> 0 collapse gets run from a configuration file.
BOOST_AUTO_TEST_CASE( a_rotating_source_without_omega_is_not_rotating )
{
	Configuration config = parse( withSource(
		"[source]\n"
		"Type = \"rotating\"\n"
		"GGPrime = 0.25\n"
		"\n[[source.species]]\n"
		"Mass = 3.344e-27\n"
		"Charge = 1.0\n"
		"Temperature = 1.6e-16\n"
		"Density = 5.0e19\n"
		"\n[[source.species]]\n"
		"Mass = 9.109e-31\n"
		"Charge = -1.0\n"
		"Temperature = 1.6e-16\n"
		"Neutralising = true\n" ) );

	meq::RotatingParameters const & rotating = config.getSource().getRotating();

	BOOST_TEST( rotating.omegaGiven == false );
	BOOST_TEST( rotating.omega == 0.0 );
	BOOST_TEST( rotating.omegaFile == "" );
	BOOST_TEST( rotating.omegaScale == 1.0 );

	// A constant g g' this time, and the defaults for everything optional.
	BOOST_TEST( rotating.ggPrime == 0.25 );
	BOOST_TEST( rotating.ggPrimeFile == "" );
	BOOST_TEST( rotating.ggPrimeScale == 1.0 );
	BOOST_TEST( rotating.referenceRadius == 1.0 );
	BOOST_TEST( rotating.mu0 == 4.0e-7*3.14159265358979323846 );

	// An unnamed species is named by its index, so the output variables of a
	// file that did not bother still differ from one another.
	BOOST_TEST_REQUIRE( rotating.species.size() == 2u );
	BOOST_TEST( rotating.species[ 0 ].name == "species0" );
	BOOST_TEST( rotating.species[ 1 ].name == "species1" );
	BOOST_TEST( rotating.species[ 0 ].temperatureScale == 1.0 );
	BOOST_TEST( rotating.species[ 0 ].densityScale == 1.0 );
}

/*
 * Normalised flux, on both sources that take it.
 *
 * Normalised = true says the profiles are functions of Psi = psi/psi_ax rather
 * than of psi, which makes psi_ax a functional of the solution and therefore an
 * UNKNOWN of the non-linear system. isNormalised() is what lets the driver
 * branch to setSource( NormalisedSource &, double ) without switching over the
 * source type itself, and psiAxisGuess() is where the bordered Newton starts.
 */
BOOST_AUTO_TEST_CASE( normalised_profiles_make_psi_axis_an_unknown )
{
	Configuration rotatingConfig = parse( rotating(
		"Normalised = true\nPsiAxis = 0.115\n" ) );

	BOOST_TEST( rotatingConfig.getSource().getRotating().normalised == true );
	BOOST_TEST( rotatingConfig.getSource().getRotating().psiAxis == 0.115 );
	BOOST_TEST( rotatingConfig.getSource().isNormalised() == true );
	BOOST_TEST( rotatingConfig.getSource().psiAxisGuess() == 0.115 );

	Configuration mhdConfig = parse( withSource(
		"[source]\n"
		"Type = \"mhd\"\n"
		"PPrimeFile = \"profiles/pprime.dat\"\n"
		"GGPrimeFile = \"profiles/ggprime.dat\"\n"
		"PPrimeScale = 1.0e6\n"
		"GGPrimeScale = -2.0\n"
		"Normalised = true\n"
		"PsiAxis = -0.42\n" ) );

	meq::MHDParameters const & mhd = mhdConfig.getSource().getMHD();
	BOOST_TEST( mhd.normalised == true );
	BOOST_TEST( mhd.psiAxis == -0.42 );
	BOOST_TEST( mhd.pPrimeScale == 1.0e6 );
	// A negative scale is a sign convention like any other, and is accepted.
	BOOST_TEST( mhd.ggPrimeScale == -2.0 );
	BOOST_TEST( mhdConfig.getSource().isNormalised() == true );
	BOOST_TEST( mhdConfig.getSource().psiAxisGuess() == -0.42 );

	// A source with no normalisation at all answers both questions rather than
	// throwing: the driver asks before it knows which type it has.
	Configuration soloviev = parse( minimal() );
	BOOST_TEST( soloviev.getSource().isNormalised() == false );
	BOOST_TEST( soloviev.getSource().psiAxisGuess() == 0.0 );
}

/// PsiAxis is required exactly when Normalised is true, and refused otherwise.
/// The refusal is the interesting half: a file carrying a PsiAxis that nothing
/// reads is a file whose author believed something false about what the run was
/// doing, and it would produce a plausible answer to a different problem.
BOOST_AUTO_TEST_CASE( psi_axis_is_required_exactly_when_normalised )
{
	BOOST_CHECK_EXCEPTION( parse( rotating( "PsiAxis = 0.1\n" ) ), ConfigError,
		[]( ConfigError const & e )
		{
			return e.getKey() == "source.PsiAxis" && mentions( e, "Normalised = true" );
		} );

	BOOST_CHECK_EXCEPTION( parse( rotating( "Normalised = true\n" ) ), ConfigError,
		[]( ConfigError const & e )
		{
			return e.getKey() == "source.PsiAxis" && mentions( e, "required" );
		} );

	// Psi = psi/psi_ax is undefined at zero, so the guess cannot be it.
	BOOST_CHECK_EXCEPTION( parse( rotating( "Normalised = true\nPsiAxis = 0.0\n" ) ), ConfigError,
		[]( ConfigError const & e )
		{
			return e.getKey() == "source.PsiAxis" && mentions( e, "non-zero" );
		} );

	// The same three, on the other source that takes them.
	std::string const mhd =
		"[source]\n"
		"Type = \"mhd\"\n"
		"PPrimeFile = \"p.dat\"\n"
		"GGPrimeFile = \"gg.dat\"\n";

	BOOST_CHECK_EXCEPTION( parse( withSource( mhd + "PsiAxis = 0.3\n" ) ), ConfigError,
		[]( ConfigError const & e )
		{
			return e.getKey() == "source.PsiAxis" && mentions( e, "Normalised = true" );
		} );

	BOOST_CHECK_EXCEPTION( parse( withSource( mhd + "Normalised = true\n" ) ), ConfigError,
		[]( ConfigError const & e ) { return e.getKey() == "source.PsiAxis"; } );
}

/// Two species is the floor, and it is quasineutrality's requirement rather
/// than a convenience: (97) determines phi_0 by Sum_s Z_s n_s = 0, which needs
/// the sum to be able to change sign.
BOOST_AUTO_TEST_CASE( a_rotating_source_needs_two_species_of_both_signs )
{
	std::string const oneSpecies = withSource(
		"[source]\n"
		"Type = \"rotating\"\n"
		"GGPrime = 0.5\n"
		"\n[[source.species]]\n"
		"Mass = 3.344e-27\n"
		"Charge = 1.0\n"
		"Temperature = 1.6e-16\n"
		"Neutralising = true\n" );

	BOOST_CHECK_EXCEPTION( parse( oneSpecies ), ConfigError,
		[]( ConfigError const & e )
		{
			return e.getKey() == "source.species"
			    && mentions( e, "at least two species" )
			    && mentions( e, "[[source.species]]" );
		} );

	// No species at all reaches the same check rather than a null dereference.
	std::string const noSpecies = withSource(
		"[source]\n"
		"Type = \"rotating\"\n"
		"GGPrime = 0.5\n" );

	BOOST_CHECK_EXCEPTION( parse( noSpecies ), ConfigError,
		[]( ConfigError const & e )
		{
			return e.getKey() == "source.species" && mentions( e, "at least two species" );
		} );

	// Two species, but both positive: nothing for the potential to hold in
	// balance, and (97) has no root at all.
	std::string const sameSign = withSource(
		"[source]\n"
		"Type = \"rotating\"\n"
		"GGPrime = 0.5\n"
		"\n[[source.species]]\n"
		"Mass = 3.344e-27\n"
		"Charge = 1.0\n"
		"Temperature = 1.6e-16\n"
		"Density = 5.0e19\n"
		"\n[[source.species]]\n"
		"Mass = 2.0e-26\n"
		"Charge = 6.0\n"
		"Temperature = 1.6e-16\n"
		"Neutralising = true\n" );

	BOOST_CHECK_EXCEPTION( parse( sameSign ), ConfigError,
		[]( ConfigError const & e )
		{
			return e.getKey() == "source.species"
			    && mentions( e, "both signs" )
			    && mentions( e, "quasineutrality" );
		} );
}

/// EXACTLY ONE species is Neutralising, and the message says why it is one
/// rather than none or all: fixing the gauge removes exactly one function's
/// worth of freedom from the densities, so n species carry n - 1 independent
/// density profiles and the missing one has to come from somewhere.
BOOST_AUTO_TEST_CASE( exactly_one_species_is_neutralising )
{
	std::string const none = withSource(
		"[source]\n"
		"Type = \"rotating\"\n"
		"GGPrime = 0.5\n"
		"\n[[source.species]]\n"
		"Mass = 3.344e-27\nCharge = 1.0\nTemperature = 1.6e-16\nDensity = 5.0e19\n"
		"\n[[source.species]]\n"
		"Mass = 9.109e-31\nCharge = -1.0\nTemperature = 1.6e-16\nDensity = 5.0e19\n" );

	BOOST_CHECK_EXCEPTION( parse( none ), ConfigError,
		[]( ConfigError const & e )
		{
			return e.getKey() == "source.species"
			    && mentions( e, "exactly one" )
			    && mentions( e, "n - 1" );
		} );

	std::string const both = withSource(
		"[source]\n"
		"Type = \"rotating\"\n"
		"GGPrime = 0.5\n"
		"\n[[source.species]]\n"
		"Mass = 3.344e-27\nCharge = 1.0\nTemperature = 1.6e-16\nNeutralising = true\n"
		"\n[[source.species]]\n"
		"Mass = 9.109e-31\nCharge = -1.0\nTemperature = 1.6e-16\nNeutralising = true\n" );

	BOOST_CHECK_EXCEPTION( parse( both ), ConfigError,
		[]( ConfigError const & e )
		{
			return e.getKey() == "source.species"
			    && mentions( e, "exactly one" )
			    && mentions( e, "n - 1" );
		} );
}

/// The density is given for every species but the Neutralising one, and given
/// for that one is as much an error as missing from another. A Density that
/// charge neutrality then overwrites is a number the file states and the run
/// ignores.
BOOST_AUTO_TEST_CASE( the_neutralising_species_is_the_one_without_a_density )
{
	BOOST_CHECK_EXCEPTION( parse( rotating( "", "", "Density = 5.0e19\n" ) ), ConfigError,
		[]( ConfigError const & e )
		{
			return e.getKey() == "source.species[1].Density"
			    && mentions( e, "Neutralising" );
		} );

	// The same thing said as a file rather than a constant.
	BOOST_CHECK_EXCEPTION( parse( rotating( "", "", "DensityFile = \"ne.dat\"\n" ) ), ConfigError,
		[]( ConfigError const & e ) { return e.getKey() == "source.species[1].Density"; } );

	std::string const noDensity = withSource(
		"[source]\n"
		"Type = \"rotating\"\n"
		"GGPrime = 0.5\n"
		"\n[[source.species]]\n"
		"Mass = 3.344e-27\nCharge = 1.0\nTemperature = 1.6e-16\n"
		"\n[[source.species]]\n"
		"Mass = 9.109e-31\nCharge = -1.0\nTemperature = 1.6e-16\nNeutralising = true\n" );

	BOOST_CHECK_EXCEPTION( parse( noDensity ), ConfigError,
		[]( ConfigError const & e )
		{
			return e.getKey() == "source.species[0].Density"
			    && mentions( e, "required" )
			    && mentions( e, "Neutralising" );
		} );
}

/// Constants that cannot describe a species. Charge is Z_s, signed and
/// dimensionless, so zero is not a plasma species; and a non-positive mass or
/// reference radius is a number the closure divides by.
BOOST_AUTO_TEST_CASE( unphysical_species_constants_are_refused )
{
	std::string const zeroCharge = withSource(
		"[source]\n"
		"Type = \"rotating\"\n"
		"GGPrime = 0.5\n"
		"\n[[source.species]]\n"
		"Mass = 3.344e-27\nCharge = 1.0\nTemperature = 1.6e-16\nDensity = 5.0e19\n"
		"\n[[source.species]]\n"
		"Mass = 9.109e-31\nCharge = 0.0\nTemperature = 1.6e-16\nNeutralising = true\n" );

	BOOST_CHECK_EXCEPTION( parse( zeroCharge ), ConfigError,
		[]( ConfigError const & e )
		{
			return e.getKey() == "source.species[1].Charge" && mentions( e, "Z_s" );
		} );

	std::string const zeroMass = withSource(
		"[source]\n"
		"Type = \"rotating\"\n"
		"GGPrime = 0.5\n"
		"\n[[source.species]]\n"
		"Mass = 3.344e-27\nCharge = 1.0\nTemperature = 1.6e-16\nDensity = 5.0e19\n"
		"\n[[source.species]]\n"
		"Mass = 0.0\nCharge = -1.0\nTemperature = 1.6e-16\nNeutralising = true\n" );

	BOOST_CHECK_EXCEPTION( parse( zeroMass ), ConfigError,
		[]( ConfigError const & e )
		{
			return e.getKey() == "source.species[1].Mass" && mentions( e, "positive" );
		} );

	std::string const negativeMass = withSource(
		"[source]\n"
		"Type = \"rotating\"\n"
		"GGPrime = 0.5\n"
		"\n[[source.species]]\n"
		"Mass = -3.344e-27\nCharge = 1.0\nTemperature = 1.6e-16\nDensity = 5.0e19\n"
		"\n[[source.species]]\n"
		"Mass = 9.109e-31\nCharge = -1.0\nTemperature = 1.6e-16\nNeutralising = true\n" );

	BOOST_CHECK_EXCEPTION( parse( negativeMass ), ConfigError,
		[]( ConfigError const & e ) { return e.getKey() == "source.species[0].Mass"; } );

	// The gauge: phi_0 vanishes on r = ReferenceRadius, so a non-positive one
	// is not a curve in the domain at all.
	BOOST_CHECK_EXCEPTION( parse( rotating( "ReferenceRadius = 0.0\n" ) ), ConfigError,
		[]( ConfigError const & e )
		{
			return e.getKey() == "source.ReferenceRadius" && mentions( e, "phi_0" );
		} );

	BOOST_CHECK_EXCEPTION( parse( rotating( "ReferenceRadius = -1.65\n" ) ), ConfigError,
		[]( ConfigError const & e ) { return e.getKey() == "source.ReferenceRadius"; } );

	BOOST_CHECK_EXCEPTION( parse( rotating( "Mu0 = 0.0\n" ) ), ConfigError,
		[]( ConfigError const & e ) { return e.getKey() == "source.Mu0"; } );
}

/*
 * A profile is a constant OR a table, never both and never a scale for
 * neither.
 *
 * Two keys rather than one that dispatches on node type, for the reason
 * recorded at the top of Config.cpp: TOML distinguishes 1 from 1.0, so a single
 * key meaning "a number or a filename" would fail by silently defaulting rather
 * than by refusing. These are what makes that pay.
 */
BOOST_AUTO_TEST_CASE( a_profile_is_a_constant_or_a_file_but_not_both )
{
	BOOST_CHECK_EXCEPTION(
		parse( rotating( "", "TemperatureFile = \"ti.dat\"\n" ) ), ConfigError,
		[]( ConfigError const & e )
		{
			return e.getKey() == "source.species[0].Temperature"
			    && mentions( e, "TemperatureFile" )
			    && mentions( e, "must not both be given" );
		} );

	// A scale with nothing to scale is the misspelling that would otherwise be
	// silent: the profile is missing and the file says how to convert it.
	std::string const scaleAlone = withSource(
		"[source]\n"
		"Type = \"rotating\"\n"
		"GGPrime = 0.5\n"
		"\n[[source.species]]\n"
		"Mass = 3.344e-27\nCharge = 1.0\nTemperature = 1.6e-16\nDensity = 5.0e19\n"
		"\n[[source.species]]\n"
		"Mass = 9.109e-31\nCharge = -1.0\nTemperature = 1.6e-16\n"
		"DensityScale = 2.0\nNeutralising = true\n" );

	BOOST_CHECK_EXCEPTION( parse( scaleAlone ), ConfigError,
		[]( ConfigError const & e )
		{
			return e.getKey() == "source.species[1].DensityScale"
			    && mentions( e, "means nothing without" )
			    && mentions( e, "DensityFile" );
		} );

	// Temperature is required of every species, Neutralising or not: only the
	// density is what neutrality determines.
	std::string const noTemperature = withSource(
		"[source]\n"
		"Type = \"rotating\"\n"
		"GGPrime = 0.5\n"
		"\n[[source.species]]\n"
		"Mass = 3.344e-27\nCharge = 1.0\nDensity = 5.0e19\n"
		"\n[[source.species]]\n"
		"Mass = 9.109e-31\nCharge = -1.0\nTemperature = 1.6e-16\nNeutralising = true\n" );

	BOOST_CHECK_EXCEPTION( parse( noTemperature ), ConfigError,
		[]( ConfigError const & e )
		{
			return e.getKey() == "source.species[0].Temperature"
			    && mentions( e, "TemperatureFile" );
		} );

	// The same rule on [source]'s own profiles.
	BOOST_CHECK_EXCEPTION(
		parse( rotating( "Omega = 1.0e5\nOmegaFile = \"omega.dat\"\n" ) ), ConfigError,
		[]( ConfigError const & e )
		{
			return e.getKey() == "source.Omega" && mentions( e, "OmegaFile" );
		} );

	BOOST_CHECK_EXCEPTION( parse( rotating( "OmegaScale = 2.0\n" ) ), ConfigError,
		[]( ConfigError const & e )
		{
			return e.getKey() == "source.OmegaScale" && mentions( e, "means nothing without" );
		} );

	// g g' has no "absent means zero" reading, unlike omega: a source with no
	// toroidal field function is not a tokamak.
	std::string const noGGPrime = withSource(
		"[source]\n"
		"Type = \"rotating\"\n"
		"\n[[source.species]]\n"
		"Mass = 3.344e-27\nCharge = 1.0\nTemperature = 1.6e-16\nDensity = 5.0e19\n"
		"\n[[source.species]]\n"
		"Mass = 9.109e-31\nCharge = -1.0\nTemperature = 1.6e-16\nNeutralising = true\n" );

	BOOST_CHECK_EXCEPTION( parse( noGGPrime ), ConfigError,
		[]( ConfigError const & e )
		{
			return e.getKey() == "source.GGPrime"
			    && mentions( e, "required" )
			    && mentions( e, "GGPrimeFile" );
		} );
}

/// A MISSPELT KEY INSIDE A SPECIES IS QUALIFIED ALL THE WAY DOWN, index and
/// all. That is what [[source.species]] cost: an unqualified "Mas is not a key"
/// on a five-species file sends the reader to read all five.
BOOST_AUTO_TEST_CASE( a_misspelt_species_key_names_its_element )
{
	BOOST_CHECK_EXCEPTION( parse( rotating( "", "Mas = 1.0\n" ) ), ConfigError,
		[]( ConfigError const & e )
		{
			return e.getKey() == "source.species[0].Mas"
			    && mentions( e, "did you mean 'Mass'?" )
			    && mentions( e, "accepted keys" );
		} );

	// The index is the element's own, not always zero.
	BOOST_CHECK_EXCEPTION( parse( rotating( "", "", "Charg = -1.0\n" ) ), ConfigError,
		[]( ConfigError const & e )
		{
			return e.getKey() == "source.species[1].Charg"
			    && mentions( e, "did you mean 'Charge'?" );
		} );

	// And the table it names is the species, not [source]: ReferenceRadius is
	// a real key of the parent and still has no business here.
	BOOST_CHECK_EXCEPTION( parse( rotating( "", "ReferenceRadius = 1.0\n" ) ), ConfigError,
		[]( ConfigError const & e )
		{
			return e.getKey() == "source.species[0].ReferenceRadius"
			    && mentions( e, "[source.species[0]]" );
		} );
}

/// species IS AN ARRAY OF TABLES, and the ways of writing something else are
/// worth naming because [[a.b]] and [a.b] differ by two characters and nest
/// identically. This is the first array of tables in meq's schema, so the
/// diagnostic is new code rather than a path anything else exercises.
BOOST_AUTO_TEST_CASE( species_must_be_an_array_of_tables )
{
	// A single [source.species] table -- one bracket, not two.
	std::string const singleTable = withSource(
		"[source]\n"
		"Type = \"rotating\"\n"
		"GGPrime = 0.5\n"
		"\n[source.species]\n"
		"Mass = 3.344e-27\nCharge = 1.0\nTemperature = 1.6e-16\nDensity = 5.0e19\n" );

	BOOST_CHECK_EXCEPTION( parse( singleTable ), ConfigError,
		[]( ConfigError const & e )
		{
			return e.getKey() == "source.species"
			    && mentions( e, "[[source.species]]" )
			    && mentions( e, "table" );
		} );

	// A scalar.
	std::string const scalar = withSource(
		"[source]\n"
		"Type = \"rotating\"\n"
		"GGPrime = 0.5\n"
		"species = 3\n" );

	BOOST_CHECK_EXCEPTION( parse( scalar ), ConfigError,
		[]( ConfigError const & e )
		{
			return e.getKey() == "source.species" && mentions( e, "array of tables" );
		} );

	// An array of something that is not a table, which is refused BY ELEMENT
	// so the index says which one.
	std::string const arrayOfNumbers = withSource(
		"[source]\n"
		"Type = \"rotating\"\n"
		"GGPrime = 0.5\n"
		"species = [ 1.0, 2.0 ]\n" );

	BOOST_CHECK_EXCEPTION( parse( arrayOfNumbers ), ConfigError,
		[]( ConfigError const & e )
		{
			return e.getKey() == "source.species[0]" && mentions( e, "must be a table" );
		} );
}

/*
 * MORE SPECIES THAN THE SOURCE CAN CARRY.
 *
 * Config.hpp says of RotatingParameters::species: "Between two and
 * meq::maxSpecies, of both charge signs, exactly one of them Neutralising."
 * The parser enforces the second and third of those. This asserts the first,
 * which is what the header claims and is the behaviour that is wanted -- a cap
 * the file is refused against, rather than one it is built against and thrown
 * out of a layer down.
 *
 * meq::RotatingSource caps the count at maxSpecies = 8 so that the
 * per-quadrature-point work allocates nothing and holds no mutable scratch,
 * which a source evaluated from a threaded assembly must not have. It is not a
 * number a configuration file can talk its way past, so the file is the place
 * to say so -- with the two keys around it, in one message, rather than after
 * every profile file named by the nine species has been opened and read.
 */
BOOST_AUTO_TEST_CASE( more_species_than_the_source_can_carry_are_refused )
{
	std::string source =
		"[source]\n"
		"Type = \"rotating\"\n"
		"GGPrime = 0.5\n";

	// Nine: one electron species to be the Neutralising one, and eight ions,
	// which is one past meq::maxSpecies.
	source +=
		"\n[[source.species]]\n"
		"Mass = 9.109e-31\nCharge = -1.0\nTemperature = 1.6e-16\nNeutralising = true\n";

	for ( int i = 0; i < 8; ++i )
		source +=
			"\n[[source.species]]\n"
			"Mass = 3.344e-27\nCharge = 1.0\nTemperature = 1.6e-16\nDensity = 5.0e18\n";

	BOOST_TEST_CONTEXT( "RED, AND THE DEFECT IS IN Config.cpp RATHER THAN HERE: the "
	                    "rotating branch checks species.size() < 2 and not "
	                    "species.size() > meq::maxSpecies, so nine species parse. The "
	                    "fix is one check beside the existing one, failing on "
	                    "source.species with meq::maxSpecies named in the message; "
	                    "SourceFactoryTests::the_factory_refuses_more_species_than_maxSpecies "
	                    "pins where it fires today." )
	{
		BOOST_CHECK_EXCEPTION( parse( withSource( source ) ), ConfigError,
			[]( ConfigError const & e )
			{
				return e.getKey() == "source.species" && mentions( e, "maxSpecies" );
			} );
	}
}

BOOST_AUTO_TEST_SUITE_END()

/*
 * [initialguess] and [adaptivity], added with the driver.
 *
 * The Amplitude check is the one that matters. Every GS-2 section 4.2-4.5
 * source vanishes at psi = 0, so with homogeneous data psi = 0 SOLVES the
 * problem and Newton stops on it in zero iterations -- see CLAUDE.md under
 * Traps. A ramp of amplitude zero IS homogeneous data, so accepting it would
 * hand the user the exact failure the option exists to avoid, and it would look
 * like a converged run rather than an error.
 */
BOOST_AUTO_TEST_CASE( initialGuessAndAdaptivityParse )
{
	auto const write = []( std::string const &body )
	{
		std::string const path = "config-test-driver-sections.toml";
		std::ofstream file( path );
		file << "[mesh]\nRMin = 0.1\nRMax = 1.9\nZMin = -1.7\nZMax = 1.7\n"
		     << "NR = 3\nNZ = 4\n\n[discretisation]\nPolynomialDegree = 2\n\n"
		     << "[source]\nType = \"soloviev\"\nA = -0.52\n\n"
		     << body;
		return path;
	};

	// Defaults: no guess, no adaptivity, and a 129x129 output grid.
	{
		meq::Configuration const config( write( "" ) );
		BOOST_TEST( ( config.getInitialGuess().type == meq::InitialGuessType::None ) );
		BOOST_TEST( config.getAdaptivity().enabled == false );
		BOOST_TEST( config.getOutput().gridNR == 129 );
		BOOST_TEST( config.getOutput().gridNZ == 129 );
	}

	{
		meq::Configuration const config( write(
			"[initialguess]\nType = \"ramp\"\nAmplitude = 0.25\n\n"
			"[adaptivity]\nEnabled = true\nStrategy = \"maximum\"\n"
			"MaxIterations = 4\nTheta = 0.3\nTargetError = 1.0e-4\n\n"
			"[output]\nGridNR = 65\nGridNZ = 33\n" ) );

		BOOST_TEST( ( config.getInitialGuess().type == meq::InitialGuessType::Ramp ) );
		BOOST_TEST( config.getInitialGuess().amplitude == 0.25 );
		BOOST_TEST( config.getAdaptivity().enabled == true );
		BOOST_TEST( ( config.getAdaptivity().strategy == meq::MarkingStrategy::Maximum ) );
		BOOST_TEST( config.getAdaptivity().maxIterations == 4 );
		BOOST_TEST( config.getAdaptivity().theta == 0.3 );
		BOOST_TEST( config.getOutput().gridNR == 65 );
		BOOST_TEST( config.getOutput().gridNZ == 33 );
	}

	// A ramp of zero amplitude is homogeneous data by another name.
	BOOST_CHECK_THROW(
		meq::Configuration( write( "[initialguess]\nType = \"ramp\"\nAmplitude = 0.0\n" ) ),
		meq::ConfigError );

	// gridfunction without the mesh it lives on cannot be read at all.
	BOOST_CHECK_THROW(
		meq::Configuration( write( "[initialguess]\nType = \"gridfunction\"\nFile = \"a.gf\"\n" ) ),
		meq::ConfigError );

	BOOST_CHECK_THROW(
		meq::Configuration( write( "[initialguess]\nType = \"warm\"\n" ) ),
		meq::ConfigError );
	BOOST_CHECK_THROW(
		meq::Configuration( write( "[adaptivity]\nStrategy = \"greedy\"\n" ) ),
		meq::ConfigError );
	BOOST_CHECK_THROW(
		meq::Configuration( write( "[adaptivity]\nTheta = 1.5\n" ) ),
		meq::ConfigError );
	BOOST_CHECK_THROW(
		meq::Configuration( write( "[output]\nGridNR = 1\n" ) ),
		meq::ConfigError );

	// Enabled is a boolean, and a string that looks like one is not one. This
	// is the find_or<double> trap in the other direction, recorded at the top
	// of Config.cpp: a silently defaulted false would disable adaptivity for a
	// user who asked for it.
	BOOST_CHECK_THROW(
		meq::Configuration( write( "[adaptivity]\nEnabled = \"true\"\n" ) ),
		meq::ConfigError );

	std::remove( "config-test-driver-sections.toml" );
}

/*
 * ============================================================================
 * [[coils]] -- FB-6's configuration half
 * ============================================================================
 *
 * FREE-BOUNDARY-PLAN.md section 5.4 calls the coils "ordinary": coil currents
 * are data, and the source adds F_coil = mu0 r I_k / |Omega_ck| on each coil
 * subdomain. meq::Coil has been the library half since FB-2 and is measured;
 * this is the file half.
 *
 * MEQ's SECOND array of tables, after [[source.species]], and the first at the
 * TOP LEVEL -- which is why Table gained a root() factory. Its elements are
 * named "coils[i]", so a fault in the third block says so.
 */
BOOST_AUTO_TEST_CASE( coils_are_given_a_total_current_or_a_uniform_density )
{
	auto const withCoil = []( std::string const & body )
	{
		return withSource( "[source]\nType = \"soloviev\"\nA = 0.5\n\n"
		                   "[boundary]\nType = \"zero\"\n\n"
		                   "[[coils]]\n" + body );
	};

	// The TOTAL current, in amperes.
	{
		Configuration const c = parse( withCoil( "Name = \"PF1\"\n"
		                                  "CentreR = 1.6\nCentreZ = 0.8\n"
		                                  "HalfWidth = 0.1\nHalfHeight = 0.05\n"
		                                  "Current = 1.25e6\n" ) );
		BOOST_TEST_REQUIRE( c.getCoils().coils.size() == 1u );
		meq::CoilParameters const & coil = c.getCoils().coils.front();
		BOOST_TEST( coil.name == "PF1" );
		BOOST_TEST( coil.centreR == 1.6 );
		BOOST_TEST( coil.current == 1.25e6 );
		BOOST_TEST( coil.densityGiven == false );
	}

	// THE DENSITY, RESOLVED THROUGH THE AREA. 4 * 0.1 * 0.05 = 0.02 m^2, so
	// 5e7 A/m^2 is 1e6 A.
	//
	// TO A TIGHT RELATIVE TOLERANCE AND NOT EXACTLY, and the first version of
	// this asserted equality on the grounds that "it is one multiplication".
	// It is not: floating-point multiplication does not associate, and the
	// parse computes area = 4*hw*hh and then density*area where the assertion
	// wrote density*4*hw*hh left to right. Those differ in the last bit --
	// 1000000.0000000002 against 1000000 -- so an exact test would have been a
	// test of the code's ASSOCIATION ORDER, which is not the property wanted.
	// The property wanted is the area, and a wrong area is wrong by a factor,
	// not by an ulp.
	{
		Configuration const c = parse( withCoil( "CentreR = 1.6\nCentreZ = 0.8\n"
		                                  "HalfWidth = 0.1\nHalfHeight = 0.05\n"
		                                  "CurrentDensity = 5.0e7\n" ) );
		meq::CoilParameters const & coil = c.getCoils().coils.front();
		double const expected = 5.0e7*4.0*0.1*0.05;
		BOOST_TEST( std::fabs( coil.current - expected )
		            < 1.0e-14*std::fabs( expected ) );
		BOOST_TEST( coil.densityGiven == true );
		// Unnamed blocks get a positional name, so a diagnostic can still say
		// which one.
		BOOST_TEST( coil.name == "coil0" );
	}

	// NEITHER, and BOTH, are refused. Both is the interesting one: an author
	// who writes both has two numbers in mind, and honouring one by precedence
	// is how a coil set ends up carrying a current nobody chose.
	auto const refuses = []( std::string const & text, std::string const & key )
	{
		BOOST_CHECK_EXCEPTION( parse( text ), ConfigError,
			[&]( ConfigError const & e ) { return e.getKey() == key; } );
	};

	refuses( withCoil( "CentreR = 1.6\nCentreZ = 0.8\n"
	                   "HalfWidth = 0.1\nHalfHeight = 0.05\n" ), "coils[0].Current" );
	refuses( withCoil( "CentreR = 1.6\nCentreZ = 0.8\n"
	                   "HalfWidth = 0.1\nHalfHeight = 0.05\n"
	                   "Current = 1.0e6\nCurrentDensity = 5.0e7\n" ), "coils[0].Current" );
}

/// A coil that reaches the axis is refused HERE as well as in meq::Coil, and
/// the reason is the operator rather than the class: Grad-Shafranov carries a
/// 1/r that is not integrable through r = 0. meq::BoundaryShape makes the same
/// refusal. Catching it at the parse means the diagnostic names the coil and
/// the key instead of arriving from a constructor three layers down.
BOOST_AUTO_TEST_CASE( a_coil_may_not_reach_the_axis_or_be_degenerate )
{
	auto const withCoil = []( std::string const & body )
	{
		return withSource( "[source]\nType = \"soloviev\"\nA = 0.5\n\n"
		                   "[boundary]\nType = \"zero\"\n\n"
		                   "[[coils]]\n" + body + "Current = 1.0e6\n" );
	};
	auto const refuses = []( std::string const & text, std::string const & key )
	{
		BOOST_CHECK_EXCEPTION( parse( text ), ConfigError,
			[&]( ConfigError const & e ) { return e.getKey() == key; } );
	};

	// Straddling the axis, and touching it exactly.
	refuses( withCoil( "CentreR = 0.05\nCentreZ = 0.0\n"
	                   "HalfWidth = 0.1\nHalfHeight = 0.05\n" ), "coils[0].CentreR" );
	refuses( withCoil( "CentreR = 0.1\nCentreZ = 0.0\n"
	                   "HalfWidth = 0.1\nHalfHeight = 0.05\n" ), "coils[0].CentreR" );

	// A cross-section with no area cannot carry a current density.
	refuses( withCoil( "CentreR = 1.6\nCentreZ = 0.0\n"
	                   "HalfWidth = 0.0\nHalfHeight = 0.05\n" ), "coils[0].HalfWidth" );
	refuses( withCoil( "CentreR = 1.6\nCentreZ = 0.0\n"
	                   "HalfWidth = 0.1\nHalfHeight = -0.05\n" ), "coils[0].HalfHeight" );

	// A misspelt key names its own block, which is the whole reason the
	// elements are named rather than numbered in the message.
	refuses( withCoil( "CentreR = 1.6\nCentreZ = 0.0\n"
	                   "HalfWidth = 0.1\nHalfHeight = 0.05\nCurrentDensty = 1.0\n" ),
	         "coils[0].CurrentDensty" );
}

/// NO COILS IS THE COMMON CASE AND IS NOT AN ERROR. Every fixed-boundary
/// configuration in examples/ has none, so absence has to parse; what must not
/// is a coil that cannot be built.
BOOST_AUTO_TEST_CASE( a_configuration_without_coils_has_an_empty_coil_set )
{
	Configuration const c = parse( withSource( "[source]\nType = \"soloviev\"\nA = 0.5\n\n"
	                                    "[boundary]\nType = \"zero\"\n" ) );
	BOOST_TEST( c.getCoils().coils.empty() );
}

/// Several blocks, in file order, which is what a machine is.
BOOST_AUTO_TEST_CASE( the_coil_set_keeps_the_order_the_file_gives )
{
	Configuration const c = parse( withSource(
		"[source]\nType = \"soloviev\"\nA = 0.5\n\n"
		"[boundary]\nType = \"zero\"\n\n"
		"[[coils]]\nName = \"upper\"\nCentreR = 1.6\nCentreZ = 0.8\n"
		"HalfWidth = 0.1\nHalfHeight = 0.05\nCurrent = 1.0e6\n\n"
		"[[coils]]\nName = \"lower\"\nCentreR = 1.6\nCentreZ = -0.8\n"
		"HalfWidth = 0.1\nHalfHeight = 0.05\nCurrent = -1.0e6\n" ) );

	BOOST_TEST_REQUIRE( c.getCoils().coils.size() == 2u );
	BOOST_TEST( c.getCoils().coils[ 0 ].name == "upper" );
	BOOST_TEST( c.getCoils().coils[ 1 ].name == "lower" );
	BOOST_TEST( c.getCoils().coils[ 1 ].centreZ == -0.8 );
	// Signed, and zero is allowed: a coil set with a coil switched off is a
	// perfectly ordinary machine state.
	BOOST_TEST( c.getCoils().coils[ 1 ].current == -1.0e6 );
}

/*
 * ConfineToPlasma -- the moving support, and it is normalised-only.
 *
 * The plasma is where the NORMALISED flux is positive, so the key means nothing
 * without a normalisation and is refused rather than ignored there, on exactly
 * the principle PsiAxis is: a file that asks for a moving plasma boundary and
 * silently gets a fixed one is a file whose author believed something false
 * about what the run was doing.
 *
 * IT IS OFF BY DEFAULT AND MUST STAY SO. Every fixed-boundary configuration in
 * examples/ has Omega = the plasma, Psi > 0 throughout by construction, and
 * switching this on would put a branch in the inner loop for nothing --
 * whereas switching it on by ACCIDENT would switch the source off wherever psi
 * changed sign, which on a curved boundary it does.
 */
BOOST_AUTO_TEST_CASE( the_moving_plasma_support_is_normalised_only )
{
	Configuration on = parse( rotating(
		"Normalised = true\nPsiAxis = 0.115\nConfineToPlasma = true\n" ) );
	BOOST_TEST( on.getSource().getRotating().confineToPlasma == true );
	BOOST_TEST( on.getSource().confinesToPlasma() == true );

	Configuration off = parse( rotating( "Normalised = true\nPsiAxis = 0.115\n" ) );
	BOOST_TEST( off.getSource().getRotating().confineToPlasma == false );
	BOOST_TEST( off.getSource().confinesToPlasma() == false );

	std::string const mhdBody =
		"[source]\n"
		"Type = \"mhd\"\n"
		"PPrimeFile = \"profiles/pprime.dat\"\n"
		"GGPrimeFile = \"profiles/ggprime.dat\"\n";

	Configuration mhdOn = parse( withSource(
		mhdBody + "Normalised = true\nPsiAxis = 0.3\nConfineToPlasma = true\n" ) );
	BOOST_TEST( mhdOn.getSource().getMHD().confineToPlasma == true );
	BOOST_TEST( mhdOn.getSource().confinesToPlasma() == true );

	// Refused without a normalisation, on both sources that take it.
	for ( std::string const &text : { rotating( "ConfineToPlasma = true\n" ),
	                                  withSource( mhdBody + "ConfineToPlasma = true\n" ) } )
		BOOST_CHECK_EXCEPTION( parse( text ), ConfigError,
			[]( ConfigError const & e )
			{
				return e.getKey() == "source.ConfineToPlasma"
				       && mentions( e, "Normalised = true" );
			} );

	// And a source that has no normalised form at all answers the question
	// rather than throwing, because the driver asks before it knows the type.
	BOOST_TEST( parse( minimal() ).getSource().confinesToPlasma() == false );
}

/*
 * mu0 IS THE SOURCE'S AND THE COILS SHARE IT, which is why [[coils]] has no
 * Mu0 key of its own -- two keys that must agree are a way of writing down a
 * disagreement. A run in normalised units setting [source] Mu0 = 1 while the
 * coils kept the SI value would be summing two terms scaled a million-fold
 * apart, and would converge, at full order, to a machine nobody described.
 */
BOOST_AUTO_TEST_CASE( the_coils_take_the_sources_permeability )
{
	std::string const body =
		"[source]\n"
		"Type = \"mhd\"\n"
		"PPrimeFile = \"profiles/pprime.dat\"\n"
		"GGPrimeFile = \"profiles/ggprime.dat\"\n";

	BOOST_TEST( parse( withSource( body ) ).getSource().permeability()
	            == 4.0e-7*3.14159265358979323846 );
	BOOST_TEST( parse( withSource( body + "Mu0 = 1.0\n" ) ).getSource().permeability()
	            == 1.0 );

	// A benchmark source folds mu0 into its own coefficients and has no key
	// for it; SI is the only answer that makes an SI coil current mean
	// anything, and that is what it gives.
	BOOST_TEST( parse( minimal() ).getSource().permeability()
	            == 4.0e-7*3.14159265358979323846 );

	// There is deliberately no Mu0 on a coil block.
	BOOST_CHECK_EXCEPTION(
		parse( withSource( body + "\n[[coils]]\nCentreR = 1.6\nCentreZ = 0.0\n"
		                   "HalfWidth = 0.1\nHalfHeight = 0.05\nCurrent = 1.0e6\n"
		                   "Mu0 = 1.0\n" ) ),
		ConfigError,
		[]( ConfigError const & e ) { return e.getKey() == "coils[0].Mu0"; } );
}


/*
 * THE TWO FREE-BOUNDARY BORDERS, AND WHAT THE SCHEMA REFUSES RATHER THAN
 * RESOLVES.
 *
 * setBoundaryFluxPoint() and setExteriorCoupling() have been library capability
 * since FB-3 and FB-5 and had no route from a configuration file. Wiring them
 * creates three ways for a file to describe something incoherent, and each is
 * refused with a message naming what to change -- because all three would
 * otherwise PARSE, SOLVE and CONVERGE, at full order, to an equilibrium the
 * file did not describe. That is this project's standing failure mode and the
 * reason the refusals are tested rather than the acceptances.
 */
BOOST_AUTO_TEST_CASE( the_free_boundary_borders_refuse_what_they_cannot_honour )
{
	auto const refuses = []( std::string const & text, std::string const & key )
	{
		BOOST_CHECK_EXCEPTION( parse( text ), ConfigError,
			[&]( ConfigError const & e ) { return e.getKey() == key; } );
	};

	std::string const halfDisc =
		"[mesh]\n"
		"RMin = 0.0\n"
		"RMax = 1.7\n"
		"ZMin = -1.7\n"
		"ZMax = 1.7\n"
		"\n"
		"[discretisation]\n"
		"PolynomialDegree = 2\n"
		"\n";

	std::string const normalisedSource =
		"[source]\n"
		"Type = \"mhd\"\n"
		"Normalised = true\n"
		"PsiAxis = 0.1\n"
		"PPrimeFile = \"examples/fb-pprime.dat\"\n"
		"GGPrimeFile = \"examples/fb-ggprime.dat\"\n";

	std::string const plainSource =
		"[source]\n"
		"Type = \"mhd\"\n"
		"PPrimeFile = \"examples/fb-pprime.dat\"\n"
		"GGPrimeFile = \"examples/fb-ggprime.dat\"\n";

	// THE COUPLING ALONE IS ACCEPTED, which is the control: without it the
	// refusals below would pass on a schema that refuses everything.
	BOOST_CHECK_NO_THROW( parse( halfDisc + normalisedSource
		+ "\n[boundary.exterior]\nRadius = 1.5\nModes = 4\n" ) );

	Configuration const good = parse( halfDisc + normalisedSource
		+ "\n[boundary.exterior]\nRadius = 1.5\nCentreZ = 0.25\nModes = 6\n"
		+ "\n[boundary.limiter]\nR = 1.2\nZ = 0.0\n" );
	BOOST_TEST( good.getBoundary().exterior.given );
	BOOST_TEST( good.getBoundary().exterior.radius == 1.5 );
	BOOST_TEST( good.getBoundary().exterior.centreZ == 0.25 );
	BOOST_TEST( good.getBoundary().exterior.modes == 6 );
	BOOST_TEST( good.getBoundary().limiter.given );
	BOOST_TEST( good.getBoundary().limiter.r == 1.2 );

	// A LIMITER WITHOUT A NORMALISATION. psi_bnd is read ONLY through
	// Psi = ( psi - psi_bnd )/( psi_ax - psi_bnd ), so on a plain source it is
	// an unknown nothing consumes: the border would be solved and the answer
	// discarded.
	refuses( halfDisc + plainSource + "\n[boundary.limiter]\nR = 1.2\nZ = 0.0\n",
	         "boundary.limiter.R" );

	// A LIMITER ON THE AXIS. r = 0 is not a limiter, and the nearest potential
	// dof to it is in an element whose flux mass ( r q, v ) degenerates.
	refuses( halfDisc + normalisedSource
	         + "\n[boundary.limiter]\nR = 0.0\nZ = 0.0\n",
	         "boundary.limiter.R" );

	// BOTH COORDINATES OR NEITHER. A limiter given only its height sits at
	// r = 0, which is the case above wearing a different spelling.
	refuses( halfDisc + normalisedSource + "\n[boundary.limiter]\nZ = 0.0\n",
	         "boundary.limiter.R" );

	// TWO DESCRIPTIONS OF Gamma. [boundary.exterior] makes it a semicircle about
	// the axis; [boundary.shape] makes it a closed surface that may not reach
	// r = 0. They are alternatives, and silently taking one is how a run ends up
	// solving on a domain nobody asked for.
	refuses( halfDisc + normalisedSource
	         + "\n[boundary.exterior]\nRadius = 1.5\nModes = 4\n"
	         + "\n[boundary.shape]\nType = \"miller\"\nR0 = 1.0\n"
	           "MinorRadius = 0.3\nElongation = 1.4\n",
	         "boundary.exterior.Radius" );

	// A DATUM ON Gamma BESIDE THE TRANSMISSION CONDITION. Gamma carries the
	// exterior map, not prescribed data; imposing both is two conditions on one
	// curve.
	refuses( "[mesh]\nRMin = 0.0\nRMax = 1.7\nZMin = -1.7\nZMax = 1.7\n\n"
	         "[discretisation]\nPolynomialDegree = 2\n\n"
	         "[source]\nType = \"soloviev\"\nA = -0.52\n\n"
	         "[boundary]\nType = \"exact\"\n\n"
	         "[boundary.exterior]\nRadius = 1.5\nModes = 4\n",
	         "boundary.exterior.Radius" );

	// A TRUNCATION OF NOTHING. The exterior map is exact mode by mode, so the
	// count is the only approximation in it; zero modes is a coupling that
	// transmits nothing.
	refuses( halfDisc + normalisedSource
	         + "\n[boundary.exterior]\nRadius = 1.5\nModes = 0\n",
	         "boundary.exterior.Modes" );
	refuses( halfDisc + normalisedSource
	         + "\n[boundary.exterior]\nRadius = -1.5\nModes = 4\n",
	         "boundary.exterior.Radius" );
}

/*
 * THE BUMP GUESS, AND WHY A RAMP IS NOT ENOUGH ONCE THE BOUNDARY IS FREE.
 *
 * A ramp is antisymmetric in z: it exists so that psi = 0 falls in the INTERIOR
 * and the trivial branch is not a fixed point. It describes no plasma. A
 * free-boundary problem has more than one converged solution -- the freegs4e
 * rehearsal measured three, with the physical one BETWEEN the two a ramp sweep
 * reaches -- so the guess is part of the problem statement rather than an
 * optimisation, and a core is what selects the branch that is one.
 */
BOOST_AUTO_TEST_CASE( the_bump_guess_describes_a_core_and_refuses_one_that_is_not )
{
	auto const guess = []( std::string const & body )
	{
		return minimal() + "\n[initialguess]\nType = \"bump\"\n" + body;
	};
	auto const refuses = []( std::string const & text, std::string const & key )
	{
		BOOST_CHECK_EXCEPTION( parse( text ), ConfigError,
			[&]( ConfigError const & e ) { return e.getKey() == key; } );
	};

	Configuration const good = parse( guess(
		"CentreR = 1.0\nCentreZ = 0.1\nRadiusR = 0.4\nRadiusZ = 0.6\n"
		"Amplitude = 0.2\n" ) );
	BOOST_TEST( ( good.getInitialGuess().type == meq::InitialGuessType::Bump ) );
	BOOST_TEST( good.getInitialGuess().centreR == 1.0 );
	BOOST_TEST( good.getInitialGuess().radiusZ == 0.6 );

	// RadiusZ defaults to RadiusR, which is the circular case.
	Configuration const circular = parse( guess(
		"CentreR = 1.0\nRadiusR = 0.4\nAmplitude = 0.2\n" ) );
	BOOST_TEST( circular.getInitialGuess().radiusZ == 0.4 );

	// A BUMP REACHING THE AXIS. psi is pinned there and the operator's 1/r is
	// not integrable through it, so a guess straddling r = 0 describes no
	// plasma however plausible its peak.
	refuses( guess( "CentreR = 0.3\nRadiusR = 0.4\nAmplitude = 0.2\n" ),
	         "initialguess.RadiusR" );
	refuses( guess( "CentreR = 0.0\nRadiusR = 0.4\nAmplitude = 0.2\n" ),
	         "initialguess.CentreR" );

	// A FLAT BUMP IS THE TRIVIAL BRANCH, which is exactly what a guess exists
	// to keep the iteration off.
	refuses( guess( "CentreR = 1.0\nRadiusR = 0.4\nAmplitude = 0.0\n" ),
	         "initialguess.Amplitude" );
}


/*
 * [source] PlasmaCurrent: the third border, and the unit it is spelled in.
 *
 * meq::GradShafranovSolver::setPlasmaCurrent takes `mu0 I_p` deliberately --
 * everything inside the solver already speaks in it, and a solver taking amperes
 * would need a mu0 of its own that could disagree with the source's. THE
 * CONFIGURATION LAYER HAS NO SUCH PROBLEM: the file names exactly one mu0, under
 * [source], and it is the same one [[coils]] uses for its Current. So the key is
 * in AMPERES and the driver multiplies once, which is the only spelling a user
 * can write without knowing which of two mu0 values is meant.
 */
BOOST_AUTO_TEST_CASE( the_prescribed_plasma_current_is_amperes_and_needs_a_normalisation )
{
	auto const refuses = []( std::string const & text, std::string const & key )
	{
		BOOST_CHECK_EXCEPTION( parse( text ), ConfigError,
			[&]( ConfigError const & e ) { return e.getKey() == key; } );
	};

	auto const mhd = []( std::string const & extra )
	{
		return withSource(
			"[source]\n"
			"Type = \"mhd\"\n"
			"PPrimeFile = \"examples/fb-pprime.dat\"\n"
			"GGPrimeFile = \"examples/fb-ggprime.dat\"\n" + extra );
	};

	Configuration const good = parse( mhd(
		"Normalised = true\nPsiAxis = 0.1\nPlasmaCurrent = 3.0e5\n" ) );
	BOOST_TEST( good.getSource().plasmaCurrent() == 3.0e5 );

	// ABSENT MEANS "DO NOT CONSTRAIN", and reads as zero.
	Configuration const none = parse( mhd( "Normalised = true\nPsiAxis = 0.1\n" ) );
	BOOST_TEST( none.getSource().plasmaCurrent() == 0.0 );

	// A NEGATIVE CURRENT IS ORDINARY. The sign is the direction of the toroidal
	// current and refusing it would refuse half the machines in the world.
	Configuration const reversed = parse( mhd(
		"Normalised = true\nPsiAxis = 0.1\nPlasmaCurrent = -3.0e5\n" ) );
	BOOST_TEST( reversed.getSource().plasmaCurrent() == -3.0e5 );

	// WITHOUT A NORMALISATION there is no profile scale to make an unknown, so
	// the border would be solved and its answer would multiply nothing.
	refuses( mhd( "PlasmaCurrent = 3.0e5\n" ), "source.PlasmaCurrent" );

	// AN EXPLICIT ZERO IS NOT THE DEFAULT WEARING A DIFFERENT SPELLING. It is a
	// prescribed current of zero, which the border satisfies by driving the
	// profile scale to zero -- a converged solve describing no plasma.
	refuses( mhd( "Normalised = true\nPsiAxis = 0.1\nPlasmaCurrent = 0.0\n" ),
	         "source.PlasmaCurrent" );
}

/// [output] FluxSurfaces and the four sizes beside it: the ( Psi, theta ) file
/// of INVERSION-PLAN.md stage IN-6.
BOOST_AUTO_TEST_CASE( the_flux_surface_output_is_off_by_default_and_its_cut_is_validated )
{
	auto const refuses = []( std::string const & text, std::string const & key )
	{
		BOOST_CHECK_EXCEPTION( parse( text ), ConfigError,
			[&]( ConfigError const & e ) { return e.getKey() == key; } );
	};

	// OFF BY DEFAULT, and the defaults are the ones src/meq/FluxFamily.hpp
	// records the decision for. Asserted rather than assumed, because a
	// reduction that started writing itself on every run would cost every
	// existing configuration a contour trace per surface -- and could FAIL on a
	// run that solved perfectly well.
	Configuration const plain = parse( minimal() );
	BOOST_TEST( !plain.getOutput().fluxSurfaces );
	BOOST_TEST( plain.getOutput().fluxSurfaceCount == 24 );
	BOOST_TEST( plain.getOutput().fluxAngleCount == 128 );
	BOOST_TEST( plain.getOutput().fluxInnerCut == 0.05 );
	BOOST_TEST( plain.getOutput().fluxOuterCut == 0.95 );

	Configuration const asked = parse( minimal() +
		"\n[output]\n"
		"Prefix = \"run\"\n"
		"FluxSurfaces = true\n"
		"FluxSurfaceCount = 40\n"
		"FluxAngleCount = 512\n"
		"FluxInnerCut = 0.02\n"
		"FluxOuterCut = 0.99\n" );
	BOOST_TEST( asked.getOutput().fluxSurfaces );
	BOOST_TEST( asked.getOutput().fluxSurfaceCount == 40 );
	BOOST_TEST( asked.getOutput().fluxAngleCount == 512 );
	BOOST_TEST( asked.getOutput().fluxInnerCut == 0.02 );
	BOOST_TEST( asked.getOutput().fluxOuterCut == 0.99 );
	BOOST_TEST( asked.getOutput().getFluxSurfaceFile() == "./run_surfaces.nc" );

	// THE CUT IS REFUSED AT ITS ENDS AND NOT CLAMPED TO THEM. Psi_N = 0 is the
	// magnetic axis, where the surface is a point and drho/dpsi is unbounded;
	// Psi_N = 1 is the plasma boundary. Silently moving a caller's number to the
	// nearest one it could honour is how a run comes to describe a family
	// nobody asked for.
	refuses( minimal() + "\n[output]\nFluxInnerCut = 0.0\n",
	         "output.FluxInnerCut" );
	refuses( minimal() + "\n[output]\nFluxInnerCut = -0.1\n",
	         "output.FluxInnerCut" );
	refuses( minimal() + "\n[output]\nFluxOuterCut = 1.0\n",
	         "output.FluxOuterCut" );
	refuses( minimal() + "\n[output]\nFluxInnerCut = 0.9\nFluxOuterCut = 0.5\n",
	         "output.FluxInnerCut" );
	refuses( minimal() + "\n[output]\nFluxSurfaceCount = 1\n",
	         "output.FluxSurfaceCount" );
	refuses( minimal() + "\n[output]\nFluxAngleCount = 2\n",
	         "output.FluxAngleCount" );

	// VALIDATED WHETHER OR NOT THE FILE IS BEING WRITTEN. A key that is present
	// is a key the author meant, and every refusal above is taken with
	// FluxSurfaces absent -- so a nonsense cut costs milliseconds at parse
	// rather than the solve it would otherwise be discovered after.

	// And an integer where a float is expected still reads: toml11's
	// find_or<double> returns the DEFAULT on an integer node rather than
	// converting, which is the silent-wrong-answer trap Config.cpp's asFloat()
	// exists to close. A cut written as 1 rather than 1.0 must therefore be
	// refused rather than quietly becoming 0.95.
	refuses( minimal() + "\n[output]\nFluxOuterCut = 1\n",
	         "output.FluxOuterCut" );
}

/*
 * [ source ] SafetyFactorFile: DRIVING BY q RATHER THAN BY g g'.
 *
 * ROADMAP.md item 10 inverts the usual direction, so the key that carries the
 * toroidal field becomes an ALTERNATIVE to GGPrimeFile rather than a companion
 * to it -- and the refusal matters more than the acceptance, because both
 * spellings converge. A run given both would do whichever physics key order
 * happened to select, at full order, with nothing in its output to say so.
 */
BOOST_AUTO_TEST_CASE( the_safety_factor_target_replaces_the_toroidal_field_table )
{
	auto const refuses = []( std::string const & text, std::string const & key )
	{
		BOOST_CHECK_EXCEPTION( parse( text ), ConfigError,
			[&]( ConfigError const & e ) { return e.getKey() == key; } );
	};

	auto const mhd = []( std::string const & extra )
	{
		return withSource(
			"[source]\n"
			"Type = \"mhd\"\n"
			"PPrimeFile = \"examples/fb-pprime.dat\"\n" + extra );
	};

	std::string const driven =
		"Normalised = true\nPsiAxis = 0.3\n"
		"SafetyFactorFile = \"examples/q-driven-q.dat\"\n"
		"ToroidalFieldGuess = 2.2\n";

	Configuration const good = parse( mhd( driven ) );
	BOOST_TEST( good.getSource().getMHD().safetyFactorFile
	            == "examples/q-driven-q.dat" );
	BOOST_TEST( good.getSource().getMHD().toroidalFieldGuess == 2.2 );
	BOOST_TEST( good.getSource().getMHD().safetyFactorDegree == 2u );
	BOOST_TEST( good.getSource().getMHD().ggPrimeFile.empty(),
		"a q-driven run has no gg' table, and leaving one here would let the "
		"factory read a profile the loop is about to replace" );

	// THE ORDINARY ROUTE IS UNTOUCHED, which is what says the alternative was
	// added rather than substituted.
	Configuration const prescribed = parse( mhd(
		"GGPrimeFile = \"examples/fb-ggprime.dat\"\n" ) );
	BOOST_TEST( prescribed.getSource().getMHD().safetyFactorFile.empty() );
	BOOST_TEST( !prescribed.getSource().getMHD().ggPrimeFile.empty() );

	// BOTH IS REFUSED. One prescribes the toroidal field and the other asks for
	// whichever field delivers a given q; there is no safe precedence.
	refuses( mhd( driven + "GGPrimeFile = \"examples/fb-ggprime.dat\"\n" ),
	         "source.SafetyFactorFile" );

	// NEITHER IS REFUSED TOO, and against GGPrimeFile, because that is the key
	// an ordinary run is missing.
	refuses( mhd( "Normalised = true\nPsiAxis = 0.3\n" ), "source.GGPrimeFile" );

	// THE TARGET IS q( Psi ), so without the normalisation there is no Psi.
	refuses( mhd( "SafetyFactorFile = \"examples/q-driven-q.dat\"\n"
	              "ToroidalFieldGuess = 2.2\n" ),
	         "source.SafetyFactorFile" );

	// AND THE OPENING g HAS NO DEFAULT. q determines g through the geometry,
	// but only once there is a geometry; a machine's vacuum R0 B0 is the number
	// a user has and inventing one here would start the loop somewhere nobody
	// chose.
	refuses( mhd( "Normalised = true\nPsiAxis = 0.3\n"
	              "SafetyFactorFile = \"examples/q-driven-q.dat\"\n" ),
	         "source.ToroidalFieldGuess" );
	// A duplicate key would be refused by the TOML parser rather than by MEQ,
	// so the negative case is built WITHOUT the good guess rather than beside it.
	refuses( mhd( "Normalised = true\nPsiAxis = 0.3\n"
	              "SafetyFactorFile = \"examples/q-driven-q.dat\"\n"
	              "ToroidalFieldGuess = -1.0\n" ),
	         "source.ToroidalFieldGuess" );

	// DEGREE ZERO IS A CONSTANT g^2, so gg' is identically zero and no
	// coefficient the loop moves can change the equilibrium. A degree the
	// surface family cannot determine is refused at the other end.
	refuses( mhd( driven + "SafetyFactorDegree = 0\n" ),
	         "source.SafetyFactorDegree" );
	refuses( mhd( driven + "SafetyFactorDegree = 9\n" ),
	         "source.SafetyFactorDegree" );

	Configuration const cubic = parse( mhd( driven + "SafetyFactorDegree = 3\n" ) );
	BOOST_TEST( cubic.getSource().getMHD().safetyFactorDegree == 3u );
}
