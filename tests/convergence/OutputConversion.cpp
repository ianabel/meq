/*
 * THE `.nc` IS THE OUTPUT, AND THE OTHER FORMATS ARE MADE FROM IT.
 *
 * `meq` writes one file by default -- `<stem>.nc` -- because it carries the
 * whole discrete answer: the solve's mesh as MFEM's own text, every
 * coefficient of psi_h, q and psi*, and the boundary curve. `--gf`/`--vtk` (or
 * `[output] GridFunctions`/`VTK`) ask the run for the GLVis and ParaView files
 * as well, and `meq-convert` makes the same files from the `.nc` afterwards.
 * MEASUREMENTS.md M-297 is what that buys.
 *
 * WHAT HAS TO HOLD FOR THAT TO BE A REPRESENTATION AND NOT A LOSS, and each
 * case below is one clause:
 *
 *   * a default run writes the `.nc` and nothing else, and the `.nc` carries
 *     what the formats need;
 *   * the file's keys ask for the formats exactly as the flags do;
 *   * `meq-convert` on a run's own `.nc` writes BYTE FOR BYTE the files that
 *     run wrote -- the `.mesh`, the three `.gf`s and the VTK, bent onto Gamma
 *     -- on a curved case and on one with conductors;
 *   * the `.nc` restarts a run at its answer, from the solver's whole iterate.
 *
 * BYTES AND NOT A TOLERANCE, because nothing in the conversion should touch a
 * number: the mesh text is stored verbatim, the coefficients as doubles, and
 * the VTK stage is one library call both programs make. A tolerance would hide
 * the first defect that does touch one. THE COMPARISON IS OF ONE RUN WITH ITS
 * OWN `.nc` -- two runs are two solves, and threaded assembly makes those
 * differ in the last bits, which is the solve's business and not this file's.
 */

#define BOOST_TEST_MODULE OutputConversion
#include <boost/test/unit_test.hpp>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <regex>
#include <string>
#include <vector>

#include <sys/wait.h>

#include "SourceTree.hpp"

namespace
{
	char const *driver() { return MEQ_DRIVER_PATH; }
	char const *converter() { return MEQ_CONVERT_PATH; }

	std::string slurp( std::string const &path )
	{
		std::ifstream stream( path, std::ios::binary );
		return std::string( ( std::istreambuf_iterator<char>( stream ) ),
		                    std::istreambuf_iterator<char>() );
	}

	bool exists( std::string const &path ) { return std::filesystem::exists( path ); }

	/// A fresh directory under the build tree's working directory.
	std::string fresh( std::string const &name )
	{
		std::string const directory = "output-conversion/" + name;
		std::filesystem::remove_all( directory );
		std::filesystem::create_directories( directory );
		return directory;
	}

	/// Run @a command in @a directory, its output into @a directory/log.
	int runIn( std::string const &directory, std::string const &command )
	{
		std::string const line = "cd " + directory + " && " + command + " > log 2>&1";
		int const status = std::system( line.c_str() );
		return WIFEXITED( status ) ? WEXITSTATUS( status ) : -1;
	}

	std::string header( std::string const &nc )
	{
		std::string const scratch = nc + ".header";
		if ( std::system( ( "ncdump -h " + nc + " > " + scratch + " 2>&1" ).c_str() ) != 0 )
			return std::string();
		std::string const text = slurp( scratch );
		std::remove( scratch.c_str() );
		return text;
	}

	/// Every file the formats write for @a stem, relative to the run directory.
	std::vector<std::string> formatFiles( std::string const &stem )
	{
		return { stem + ".mesh", stem + "_psi.gf", stem + "_grad_psi.gf",
		         stem + "_psistar.gf", stem + "/" + stem + ".pvd",
		         stem + "/Cycle000000/data.pvtu", stem + "/Cycle000000/proc000000.vtu" };
	}

	/// The shipped @a example with its profile tables made absolute and
	/// @a edit applied, written into @a directory; the path of the copy.
	std::string variant( std::string const &directory, std::string const &example,
	                     std::string const &name,
	                     std::string ( *edit )( std::string const & ) )
	{
		std::string text = slurp( meq::test::sourcePath( "examples/" + example + ".toml" ) );
		std::string const examples = meq::test::sourcePath( "examples" ) + "/";
		text = std::regex_replace( text,
		                           std::regex( "(File = \")([^\"/]+\\.(dat|msh|mesh|gf)\")" ),
		                           "$1" + examples + "$2" );
		text = edit( text );
		std::string const path = directory + "/" + name + ".toml";
		std::ofstream( path ) << text;
		return name + ".toml";
	}

	/// The integer after "converged in " in a driver log, or -1.
	int newtonSteps( std::string const &log )
	{
		std::smatch match;
		if ( std::regex_search( log, match, std::regex( "converged in ([0-9]+) Newton" ) ) )
			return std::stoi( match[ 1 ] );
		return -1;
	}

	/// Whether the driver under test is a single-precision build, from its own
	/// --version. A re-solve from another start lands within the format's
	/// floor of the first, and in single that floor is measured and not
	/// assumed to be double's (CLAUDE.md, the single-precision traps).
	bool singlePrecision()
	{
		std::string const scratch = "output-conversion-version.txt";
		if ( std::system( ( std::string( driver() ) + " --version > " + scratch ).c_str() ) != 0 )
			return false;
		bool const single = slurp( scratch ).find( ", single," ) != std::string::npos;
		std::remove( scratch.c_str() );
		return single;
	}

	/// The same psi_ax: the same printed digits in double, and in single within
	/// 5e-5 relative -- measured, a re-solve there lands 2.7e-6 (fixed-h-circular
	/// from its .gf) and 1.1e-5 (limited-tokamak from its .nc) from the first,
	/// where single's own answer for the latter is 1.2e-4 from double's.
	bool sameAxis( std::string const &a, std::string const &b )
	{
		if ( a.empty() || b.empty() )
			return a == b;
		if ( !singlePrecision() )
			return a == b;
		double const x = std::stod( a ), y = std::stod( b );
		return std::abs( x - y ) <= 5.0e-5*std::abs( y );
	}

	/// The converged psi_ax as the driver prints it, without the constraint
	/// residual beside it: that is round-off, and differs with the start.
	std::string psiAxisValue( std::string const &log )
	{
		std::smatch match;
		if ( std::regex_search( log, match, std::regex( "psi_ax = ([^ ]+) Wb/rad" ) ) )
			return match[ 1 ];
		return std::string();
	}
}

// A DEFAULT RUN WRITES THE `.nc` AND NOTHING ELSE, AND THE `.nc` CARRIES WHAT
// THE OTHER FORMATS ARE MADE OF: the mesh, psi* beside psi_h and q, and the
// curve the VTK is bent onto. And it names no `mesh_file`, there being none.
BOOST_AUTO_TEST_CASE( theDefaultRunWritesTheNetCDFAndNothingElse )
{
	std::string const directory = fresh( "default" );
	BOOST_TEST_REQUIRE( runIn( directory, std::string( driver() ) + " "
	                           + meq::test::sourcePath( "examples/miller-curved.toml" ) ) == 0 );

	BOOST_TEST_REQUIRE( exists( directory + "/miller-curved.nc" ) );
	for ( std::string const &file : formatFiles( "miller-curved" ) )
		BOOST_TEST( !exists( directory + "/" + file ),
		            file << " was written by a run that asked for no format beside the .nc" );

	std::string const text = header( directory + "/miller-curved.nc" );
	for ( char const *needed : { "string mesh ;", "double psi_coefficients(",
	                             "double flux_coefficients(", "double psistar_coefficients(",
	                             ":psistar_collection", ":shape_minor_radius" } )
		BOOST_TEST( text.find( needed ) != std::string::npos,
		            "the default .nc lacks " << needed
		            << ", which meq-convert and the restart need" );
	BOOST_TEST( text.find( ":mesh_file" ) == std::string::npos,
	            "the .nc names a mesh_file that no run wrote" );
}

// THE FILE ASKS FOR THE FORMATS EXACTLY AS THE COMMAND LINE DOES.
BOOST_AUTO_TEST_CASE( theConfigurationKeysAskForTheFormats )
{
	std::string const directory = fresh( "keys" );
	std::string const config = variant( directory, "miller-curved", "keys",
		[]( std::string const &text )
		{
			return std::regex_replace( text, std::regex( "\n\\[output\\]\n" ),
			                           "\n[output]\nGridFunctions = true\nVTK = true\n" );
		} );
	BOOST_TEST_REQUIRE( runIn( directory, std::string( driver() ) + " " + config ) == 0,
	                    slurp( directory + "/log" ) );
	for ( std::string const &file : formatFiles( "miller-curved" ) )
		BOOST_TEST( exists( directory + "/" + file ),
		            file << " was not written under [output] GridFunctions and VTK" );
	BOOST_TEST( header( directory + "/miller-curved.nc" ).find( ":mesh_file" )
	            != std::string::npos,
	            "a run that wrote its .mesh does not name it in the .nc" );
}

// THE CONVERTER WRITES THE RUN'S OWN FILES, BYTE FOR BYTE. A curved case, whose
// VTK is bent onto Gamma from the curve the .nc records, and a case with
// conductors, whose VTK carries the element record's cells.
BOOST_AUTO_TEST_CASE( theConverterWritesTheRunsOwnFiles )
{
	for ( char const *example : { "miller-curved", "coils-rectangle" } )
	{
		std::string const stem = example;
		std::string const directory = fresh( "convert-" + stem );
		std::filesystem::create_directories( directory + "/converted" );
		BOOST_TEST_REQUIRE( runIn( directory, std::string( driver() ) + " --gf --vtk "
		                           + meq::test::sourcePath( "examples/" + stem + ".toml" ) )
		                    == 0 );
		BOOST_TEST_REQUIRE( runIn( directory, std::string( converter() ) + " " + stem
		                           + ".nc --output converted/" + stem ) == 0,
		                    slurp( directory + "/log" ) );

		for ( std::string const &file : formatFiles( stem ) )
		{
			std::string const ran = directory + "/" + file;
			std::string const made = directory + "/converted/" + file;
			BOOST_TEST_REQUIRE( exists( ran ), ran << " was not written by --gf --vtk" );
			BOOST_TEST_REQUIRE( exists( made ), "meq-convert did not write " << file );
			BOOST_TEST( slurp( ran ) == slurp( made ),
			            "meq-convert's " << file << " is not the run's: the .nc does "
			            "not carry what that file holds, or the conversion touched a "
			            "number" );
		}
	}
}

// THE `.nc` RESTARTS A RUN AT ITS ANSWER. It carries the solver's whole
// iterate ( q, psi, psihat ), so a restart on the mesh it was solved on starts
// there: no Newton step, and the run that results writes the same psi_h and
// psi* to the byte. The `.gf` pair carries psi_h alone -- q is re-derived and
// psihat projected -- and still takes a step or two to the same answer, which
// is the control. On a non-linear case, so that a restart that is not exact
// costs iterations and shows. MEASUREMENTS.md M-298.
BOOST_AUTO_TEST_CASE( theNetCDFRestartsAtTheAnswer )
{
	std::string const directory = fresh( "restart" );
	BOOST_TEST_REQUIRE( runIn( directory, std::string( driver() ) + " --gf "
	                           + meq::test::sourcePath( "examples/fixed-h-circular.toml" ) )
	                    == 0 );
	std::string const original = slurp( directory + "/log" );
	int const cold = newtonSteps( original );
	BOOST_TEST_REQUIRE( cold > 0 );
	std::string const here = std::filesystem::absolute( directory ).string();

	static std::string guessBlock, prefix;
	auto restart = [ & ]( std::string const &name, std::string const &block )
	{
		guessBlock = block;
		prefix = name;
		std::string const config = variant( directory, "fixed-h-circular", name,
			[]( std::string const &text )
			{
				std::string out = std::regex_replace(
					text, std::regex( "\n\\[initialguess\\]\n(?:[A-Za-z]+ = [^\n]*\n)+" ), "\n" );
				out = std::regex_replace( out, std::regex( "\n\\[source\\]\n" ),
				                          "\n" + guessBlock + "\n[source]\n" );
				return std::regex_replace( out, std::regex( "Prefix = \"[^\"]*\"" ),
				                           "Prefix = \"" + prefix + "\"" );
			} );
		BOOST_TEST_REQUIRE( runIn( directory, std::string( driver() ) + " --gf " + config )
		                    == 0, slurp( directory + "/log" ) );
		return slurp( directory + "/log" );
	};

	std::string const fromNc = restart( "from-nc",
		"[initialguess]\nType = \"gridfunction\"\nFile = \"" + here
		+ "/fixed-h-circular.nc\"\n" );
	std::string const fromGf = restart( "from-gf",
		"[initialguess]\nType = \"gridfunction\"\nFile = \"" + here
		+ "/fixed-h-circular_psi.gf\"\nMeshFile = \"" + here
		+ "/fixed-h-circular.mesh\"\nContent = \"total\"\n" );

	BOOST_TEST( fromNc.find( "psi_ax starts at the file's" ) != std::string::npos,
	            "the .nc restart did not start psi_ax from the file's psi_axis" );
	BOOST_TEST( fromNc.find( "an exact restart" ) != std::string::npos,
	            "the .nc restart did not take the file's whole iterate:\n" << fromNc );
	BOOST_TEST( newtonSteps( fromNc ) == 0,
	            "a restart from the file's own iterate took " << newtonSteps( fromNc )
	            << " Newton steps; at the answer it takes none" );
	for ( char const *suffix : { "_psi.gf", "_grad_psi.gf", "_psistar.gf" } )
		BOOST_TEST( slurp( directory + "/from-nc" + suffix )
		            == slurp( directory + "/fixed-h-circular" + suffix ),
		            "the run restarted at the answer did not write the answer's "
		            << suffix );

	// The control: psi_h alone is a restart, and not an exact one.
	BOOST_TEST( newtonSteps( fromGf ) > 0 );
	BOOST_TEST( newtonSteps( fromGf ) < cold );
	BOOST_TEST( !psiAxisValue( fromNc ).empty() );
	BOOST_TEST( psiAxisValue( fromNc ) == psiAxisValue( original ) );
	BOOST_TEST( sameAxis( psiAxisValue( fromGf ), psiAxisValue( original ) ),
	            "the .gf restart reached psi_ax = " << psiAxisValue( fromGf )
	            << " where the run it came from reached " << psiAxisValue( original ) );
}

// AND THE REST OF THE STATE: a run whose only unknowns are the field, on the
// plain Newton path, and a free-boundary run whose border carries psi_bnd and
// the I_p profile scale as well as psi_ax. Each restarted from its own .nc
// starts at its answer, or a step from it -- the support's edge is frozen at
// the answer's own labels where the cold run froze it at the sweep before's --
// and reaches the same equilibrium. MEASUREMENTS.md M-298.
BOOST_AUTO_TEST_CASE( theRestartCarriesTheWholeState )
{
	struct Case { char const *example; int most; };
	for ( Case const &c : { Case{ "coils-rectangle", 0 }, Case{ "limited-tokamak", 1 } } )
	{
		std::string const stem = c.example;
		std::string const directory = fresh( "state-" + stem );
		std::string const original = variant( directory, stem, "original",
			[]( std::string const &text ) { return text; } );
		BOOST_TEST_REQUIRE( runIn( directory, std::string( driver() ) + " " + original ) == 0,
		                    slurp( directory + "/log" ) );
		std::string const cold = slurp( directory + "/log" );
		std::string const here = std::filesystem::absolute( directory ).string();

		static std::string guessBlock;
		guessBlock = "[initialguess]\nType = \"gridfunction\"\nFile = \"" + here + "/"
		             + stem + ".nc\"\n";
		std::string const restarted = variant( directory, stem, "restarted",
			[]( std::string const &text )
			{
				std::string out = std::regex_replace(
					text, std::regex( "\n\\[initialguess\\]\n(?:[A-Za-z]+ = [^\n]*\n)+" ), "\n" );
				out = std::regex_replace( out, std::regex( "\n\\[source\\]\n" ),
				                          "\n" + guessBlock + "\n[source]\n" );
				return std::regex_replace( out, std::regex( "Prefix = \"[^\"]*\"" ),
				                           "Prefix = \"restarted\"" );
			} );
		BOOST_TEST_REQUIRE( runIn( directory, std::string( driver() ) + " " + restarted ) == 0,
		                    slurp( directory + "/log" ) );
		std::string const warm = slurp( directory + "/log" );

		BOOST_TEST( warm.find( "an exact restart" ) != std::string::npos,
		            stem << ": the restart did not take the whole iterate:\n" << warm );
		// One step more in single, measured on limited-tokamak: its stored state
		// is at single's floor rather than at round-off.
		int const most = c.most + ( singlePrecision() ? 1 : 0 );
		BOOST_TEST( newtonSteps( warm ) <= most,
		            stem << ": restarted from its own answer it took " << newtonSteps( warm )
		            << " Newton steps, against " << newtonSteps( cold ) << " cold" );
		BOOST_TEST( newtonSteps( warm ) < newtonSteps( cold ) );
		BOOST_TEST( sameAxis( psiAxisValue( warm ), psiAxisValue( cold ) ),
		            stem << ": the restart reached psi_ax = " << psiAxisValue( warm )
		            << " where the run it came from reached " << psiAxisValue( cold ) );
	}
}

// AND THE CONVERTER REFUSES WHAT IT CANNOT READ, with meq's own exit codes.
BOOST_AUTO_TEST_CASE( theConverterRefusesWhatItCannotRead )
{
	std::string const directory = fresh( "refusals" );
	BOOST_TEST( runIn( directory, converter() ) == 1, "no file named is a usage error" );
	BOOST_TEST( runIn( directory, std::string( converter() ) + " absent.nc" ) == 3,
	            "an unreadable .nc is meq's exit 3" );
	BOOST_TEST( runIn( directory, std::string( converter() ) + " --nonsense x.nc" ) == 1 );
}
