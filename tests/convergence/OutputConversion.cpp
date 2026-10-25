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
#include <sstream>
#include <string>
#include <vector>

#include <sys/wait.h>

#include "SourceTree.hpp"
#include "convergence/RoundOffFloor.hpp"

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

	/// Newton steps over a whole run: the sum of the support sweeps' `its`
	/// where the run swept, which is what a restart is meant to save, and the
	/// single solve's count where it did not.
	int totalSteps( std::string const &log )
	{
		std::size_t at = log.find( "sweep     its" );
		if ( at == std::string::npos )
			return newtonSteps( log );
		int total = 0;
		std::regex const row( "^ {5}(\\d+) +(\\d+) " );
		std::istringstream in( log.substr( at ) );
		std::string line;
		std::getline( in, line );
		while ( std::getline( in, line ) )
		{
			std::smatch m;
			if ( !std::regex_search( line, m, row ) )
				break;
			total += std::stoi( m[ 2 ] );
		}
		return total;
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

	/// The same psi_ax from a re-solve that does NOT start at the answer: the
	/// same printed digits in double, and in single within 1e-5 relative.
	/// Measured on MKL 2026.2 and 2026.1, the one such re-solve that moves is
	/// fixed-h-circular restarted from psi_h alone, its .gf, which takes a
	/// step and lands 3.5e-06 and 3.0e-06 from the first (30 and 25 eps); the
	/// axis on a mesh line moves by nothing. A restart AT the answer -- the
	/// .nc's whole state -- takes no step and is asserted to the digit, not
	/// through this.
	bool sameAxis( std::string const &a, std::string const &b )
	{
		if ( a.empty() || b.empty() )
			return a == b;
		if ( !singlePrecision() )
			return a == b;
		double const x = std::stod( a ), y = std::stod( b );
		return std::abs( x - y ) <= 1.0e-5*std::abs( y );
	}

	/// The plasma current as the driver prints it, or empty without one.
	std::string plasmaCurrentValue( std::string const &log )
	{
		std::smatch match;
		if ( std::regex_search( log, match, std::regex( "I_p = ([^ ]+) A against" ) ) )
			return match[ 1 ];
		return std::string();
	}

	/// How far the support loop's LAST sweep moved psi_ax, from the table
	/// the driver prints, or -1 for a run that did not sweep. The digits are
	/// the table's, which are psi_ax's own.
	double lastSweepMove( std::string const &log )
	{
		std::size_t const at = log.find( "sweep     its" );
		if ( at == std::string::npos )
			return -1.0;
		std::regex const row( "^ {5}\\d+ +\\d+ +([^ ]+) " );
		std::istringstream in( log.substr( at ) );
		std::string line;
		std::getline( in, line );
		std::vector<double> axes;
		while ( std::getline( in, line ) )
		{
			std::smatch m;
			if ( !std::regex_search( line, m, row ) )
				break;
			axes.push_back( std::stod( m[ 1 ] ) );
		}
		return axes.size() < 2 ? -1.0
		                       : std::abs( axes.back() - axes[ axes.size() - 2 ] );
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
// plain Newton path; a free-boundary run whose border carries psi_bnd and the
// I_p profile scale as well as psi_ax; and a diverted one, with an X-point and
// coils, whose first support sweep is seeded at the file's axis. MEASUREMENTS.md
// M-298.
//
// THE TWO THAT DO NOT SWEEP RESTART AT THEIR ANSWER, IN EITHER PRECISION: no
// Newton step, and the psi_ax and I_p the run that wrote the file printed, to
// the digit. That is what a whole state is. The stored iterate reads back at
// exactly the field residual the run ended on, 1.433e-07 on limited-tokamak
// in single, and every border row at the same value; what a restart can still
// get wrong is a border unknown the file holds and the solver is not given.
// The profile scale is the one measured: started at one rather than at the
// file's 1.0009 (GradShafranovSolver::setPlasmaCurrentScaleInitial()), the
// first field residual is 2.6e-05, the restart takes two steps, and psi_ax
// moves by 3.1e-05 and 5.6e-05 relative (MKL 2026.2, 2026.1) -- inside a
// 5e-05 band on one MKL and not the other, so no band is asserted here.
//
// THE DIVERTED ONE IS ONE MORE SUPPORT SWEEP, and is bounded by the loop's own
// convergence rather than by a tolerance. Its first sweep starts at the answer
// and takes no step; the second re-decides the support at the answer's labels,
// where the cold run's last sweep froze it at the sweep before's, so it is a
// different frozen problem and takes a step or two. It may move psi_ax by no
// more than the cold run's last sweep did, or than the 64 eps of |psi_ax| a
// solve stopped at its round-off floor is worth on this machine, whichever is
// larger: each sweep is one draw of that floor's noise (M-306), so the last
// sweep's move alone is a coin toss -- measured in single, the restart moved
// 3.8e-07 against a last sweep of 8.0e-08 and a 64 eps of 6.3e-07. In double
// the last sweep moved no printed digit, and nor does the restart.
BOOST_AUTO_TEST_CASE( theRestartCarriesTheWholeState )
{
	struct Case { char const *example; int most; };
	for ( Case const &c : { Case{ "coils-rectangle", 0 }, Case{ "limited-tokamak", 0 },
	                        Case{ "diverted-tokamak-xpoint", 2 } } )
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
		auto restartAs = [ & ]( std::string const &name,
		                        std::string ( *edit )( std::string const & ) )
		{
			static std::string prefix;
			prefix = name;
			static std::string ( *more )( std::string const & );
			more = edit;
			std::string const config = variant( directory, stem, name,
				[]( std::string const &text )
				{
					std::string out = std::regex_replace(
						text, std::regex( "\n\\[initialguess\\]\n(?:[A-Za-z]+ = [^\n]*\n)+" ), "\n" );
					out = std::regex_replace( out, std::regex( "\n\\[source\\]\n" ),
					                          "\n" + guessBlock + "\n[source]\n" );
					out = std::regex_replace( out, std::regex( "Prefix = \"[^\"]*\"" ),
					                          "Prefix = \"" + prefix + "\"" );
					return more( out );
				} );
			BOOST_TEST_REQUIRE( runIn( directory, std::string( driver() ) + " " + config ) == 0,
			                    slurp( directory + "/log" ) );
			return slurp( directory + "/log" );
		};
		std::string const warm = restartAs( "restarted",
			[]( std::string const &text ) { return text; } );

		BOOST_TEST( warm.find( "an exact restart" ) != std::string::npos,
		            stem << ": the restart did not take the whole iterate:\n" << warm );
		BOOST_TEST( totalSteps( warm ) <= c.most,
		            stem << ": restarted from its own answer it took " << totalSteps( warm )
		            << " Newton steps, against " << totalSteps( cold ) << " cold" );
		BOOST_TEST( totalSteps( warm ) < totalSteps( cold ) );
		// A restart at the answer that needs the retry rung has left the branch
		// on its first sweep: diverted-tokamak-xpoint did, when its first plasma
		// fill was seeded at a coil rather than at the file's axis.
		BOOST_TEST( warm.find( "Retrying" ) == std::string::npos,
		            stem << ": restarted at its own answer, the run needed the retry "
		            "rung:\n" << warm );

		double const sweepMove = lastSweepMove( cold );
		if ( sweepMove < 0.0 )
		{
			BOOST_TEST( psiAxisValue( warm ) == psiAxisValue( cold ),
			            stem << ": the restart reached psi_ax = " << psiAxisValue( warm )
			            << " where the run it came from reached " << psiAxisValue( cold ) );
			BOOST_TEST( plasmaCurrentValue( warm ) == plasmaCurrentValue( cold ),
			            stem << ": the restart reports I_p = " << plasmaCurrentValue( warm )
			            << " where the run it came from reported "
			            << plasmaCurrentValue( cold ) );
		}
		else
		{
			BOOST_TEST_REQUIRE( !psiAxisValue( warm ).empty() );
			double const moved = std::abs( std::stod( psiAxisValue( warm ) )
			                               - std::stod( psiAxisValue( cold ) ) );
			double const noise = meq::tests::scalarAgreementFloor(
				0.0, std::stod( psiAxisValue( cold ) ) );
			BOOST_TEST( moved <= std::max( sweepMove, noise ),
			            stem << ": the restart moved psi_ax by " << moved << " ( "
			            << psiAxisValue( cold ) << " -> " << psiAxisValue( warm )
			            << " ), more than the cold run's last support sweep did, "
			            << sweepMove << ", and more than 64 eps of it, " << noise );
		}

		/*
		 * AND AT ANOTHER CURRENT THE FILE IS NOT THE ANSWER. The field rows do
		 * not read the target, so a restart at the answer meets every one of
		 * them and only the current row can say the problem has changed --
		 * which it did not, at the first iterate, until it was evaluated there:
		 * limited-tokamak restarted at 3.3e5 A "converged in 0 Newton
		 * iterations" to the 3.0e5 equilibrium and reported I_p = 0. Measured
		 * since, 6 steps to 3.300000e+05 in single.
		 */
		if ( stem == "limited-tokamak" )
		{
			std::string const other = restartAs( "other-current",
				[]( std::string const &text )
				{
					return std::regex_replace( text, std::regex( "\nPlasmaCurrent = [^\n]*\n" ),
					                           "\nPlasmaCurrent = 3.3e5\n" );
				} );
			BOOST_TEST( other.find( "an exact restart" ) != std::string::npos, other );
			BOOST_TEST( newtonSteps( other ) > 0,
			            "restarted at another plasma current, the solve took no step:\n"
			            << other );
			BOOST_TEST_REQUIRE( !plasmaCurrentValue( other ).empty(), other );
			BOOST_TEST( std::abs( std::stod( plasmaCurrentValue( other ) ) - 3.3e5 )
			            <= 1.0e-5*3.3e5,
			            "restarted at 3.3e5 A, the solve reports I_p = "
			            << plasmaCurrentValue( other ) );
		}
	}
}

// AN AXIS ON A MESH LINE IS ONE POINT ON EVERY ROUTE. examples/rotating-
// normalised.toml is up-down symmetric on a box whose mesh has the line z = 0,
// so its axis sits on a face and q_h's jump gives it a twin in each neighbour,
// O( h^(k+1) ) apart in psi. The bordered solve's warm search, the sweep that
// reports the axis afterwards, and a restart's cold search must all name the
// same twin, or the equilibrium reported depends on where the solve started:
// it did, 2.7e-05 apart in psi_ax. MEASUREMENTS.md M-299.
BOOST_AUTO_TEST_CASE( anAxisOnAMeshLineIsOnePointOnEveryRoute )
{
	std::string const directory = fresh( "axis-on-a-line" );
	std::string const original = variant( directory, "rotating-normalised", "original",
		[]( std::string const &text ) { return text; } );
	BOOST_TEST_REQUIRE( runIn( directory, std::string( driver() ) + " " + original ) == 0,
	                    slurp( directory + "/log" ) );
	std::string const cold = slurp( directory + "/log" );

	// The constrained psi_ax against the sweep's axis, as the run reports both.
	std::smatch reported;
	BOOST_TEST_REQUIRE( std::regex_search( cold, reported,
		std::regex( "the axis, as a zero of q_h: psi = ([^ ]+) at \\([^)]*\\), normalised flux ([0-9.]+)" ) ),
		cold );
	BOOST_TEST( sameAxis( reported[ 1 ].str(), psiAxisValue( cold ) ),
	            "the solve constrained psi_ax = " << psiAxisValue( cold ) << " and the "
	            "sweep finds the axis at psi = " << reported[ 1 ].str() << ": the two "
	            "routes took different twins of an axis on a mesh line" );
	BOOST_TEST( reported[ 2 ].str() == "1.0000" );
	BOOST_TEST( cold.find( "NO O-point was reachable" ) == std::string::npos );

	std::string const here = std::filesystem::absolute( directory ).string();
	static std::string guessBlock;
	guessBlock = "[initialguess]\nType = \"gridfunction\"\nFile = \"" + here
	             + "/rotating-normalised.nc\"\n";
	std::string const restarted = variant( directory, "rotating-normalised", "restarted",
		[]( std::string const &text )
		{
			std::string out = std::regex_replace( text, std::regex( "\n\\[source\\]\n" ),
			                                      "\n" + guessBlock + "\n[source]\n" );
			return std::regex_replace( out, std::regex( "Prefix = \"[^\"]*\"" ),
			                           "Prefix = \"restarted\"" );
		} );
	BOOST_TEST_REQUIRE( runIn( directory, std::string( driver() ) + " " + restarted ) == 0,
	                    slurp( directory + "/log" ) );
	std::string const warm = slurp( directory + "/log" );
	BOOST_TEST( newtonSteps( warm ) <= ( singlePrecision() ? 1 : 0 ),
	            "restarted at its own answer, the run took " << newtonSteps( warm )
	            << " Newton steps: its axis constraint is not the cold run's" );
	BOOST_TEST( sameAxis( psiAxisValue( warm ), psiAxisValue( cold ) ),
	            "restarted, psi_ax = " << psiAxisValue( warm ) << " against the cold "
	            "run's " << psiAxisValue( cold ) );
	// And a solve that converged at once still says where its axis is.
	BOOST_TEST( warm.find( "NO O-point was reachable" ) == std::string::npos,
	            "a solve converged on its first residual reports the nodal fallback:\n"
	            << warm );
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
