/*
 * THE SHIPPED MACHINE-SHAPED FIXED-BOUNDARY CASES, AND WHAT THEY EXIST TO SAY.
 *
 * `examples/fixed-*.toml` are six real machines posed as fixed-boundary
 * problems: freegs4e solves each one FREE boundary at 257^2, its psi_n = 0.95
 * surface is fitted to MXH, and MEQ is handed that curve with psi = 0 on it.
 * psi is defined only up to an additive constant, so that is the SAME
 * equilibrium in a different gauge -- which is what makes freegs4e's answer a
 * REFERENCE here rather than a second opinion.
 *
 * THE TREE HAD NO SUCH CASE AND THE GAP WAS STRUCTURAL. Every machine-geometry
 * case in it is FREE boundary -- the half-disc reaching the axis, the DtN
 * coupling, the two borders, the moving plasma edge -- and the fixed-boundary
 * examples are analytic fixtures or one 768-element box. So there was no rung
 * between the two, and three separate questions had nowhere to be asked. The
 * TODO entry these were written for has the list; the one this file is the
 * acceptance for is the second:
 *
 *     THE PLASMA EDGE CAPS THE ORDER ON EVERY MACHINE CASE WE HAVE.
 *     `ConstrainPaxisIp` at alpha_n = 1.2 makes the edge a fractional power of
 *     1 - Psi, so the rate is capped at about 1.2 whatever the polynomial
 *     degree. A SMOOTH-EDGE machine-shaped problem is where a high-order claim
 *     can actually be made.
 *
 * At psi_n = 0.95 the edge singularity is OUTSIDE the domain and the profiles
 * are analytic across Gamma -- p' at Gamma is 0.6% of its axis value on machine
 * A and gg' is 4.1e-03, both comfortably nonzero. MEASUREMENTS.md M-147 is the
 * rate, and it reads k+2 in psi*.
 *
 * WHAT THIS FILE ASSERTS AND WHY IT IS THE DRIVER RATHER THAN THE LIBRARY.
 * These are EXAMPLES, and what rots about an example is the file: a key
 * renamed, a default moved, a profile table's path. A library-level fixture
 * would re-implement the curved path and then not notice any of that. So this
 * runs the shipped TOML through the shipped binary, exactly as a user would.
 *
 * THE REFERENCE IS AN INDEPENDENT CODE AND THE GATE IS ITS OWN ACCURACY.
 * freegs4e is finite differences and Picard where MEQ is HDG and Newton; they
 * share the equation and essentially no code. What limits the agreement is the
 * MXH fit of Gamma -- 2.7e-05 m on a minor radius of 0.3193 m for the circular
 * case, 1.9e-04 relative for machine A -- and freegs4e's own 257^2 grid.
 * Neither refines with MEQ's mesh, so the assertion below is that refinement
 * moves TOWARD the reference and then STOPS at the reference's floor. A gate
 * that kept tightening with h would be asserting something false.
 */

#define BOOST_TEST_MODULE MachineFixedBoundary
#include <boost/test/unit_test.hpp>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iterator>
#include <string>
#include <vector>

namespace
{
	/// Where the driver binary is. CMake passes it, for the reason
	/// DriverAcceptance records: guessing "../meq" is how a test starts
	/// silently not running the thing it claims to.
	char const *driver()
	{
		return MEQ_DRIVER_PATH;
	}

	std::string slurp( std::string const &path )
	{
		std::ifstream stream( path );
		return std::string( ( std::istreambuf_iterator<char>( stream ) ),
		                      std::istreambuf_iterator<char>() );
	}

	/// The scratch directory every run below writes into, so that a test run
	/// does not litter the repository root the way a bare `meq examples/x.toml`
	/// would. Relative, because ctest sets the working directory to the source
	/// tree and the shipped files name their profile tables `examples/...`.
	char const *scratch()
	{
		return "machine-fixed-boundary-scratch";
	}

	/// A shipped example with three keys overridden: the output directory, the
	/// prefix, and RefinementLevels. Returns the path of the written file.
	///
	/// TEXT SUBSTITUTION ON THE SHIPPED FILE RATHER THAN A FILE OF ITS OWN,
	/// because the whole point is to exercise what is shipped. A copy would
	/// drift, and a drifted copy passing is worse than no test.
	std::string configure( std::string const &stem, int refinement,
	                       std::string const &prefix )
	{
		std::string text = slurp( "examples/" + stem + ".toml" );
		BOOST_TEST_REQUIRE( !text.empty(),
		                    "examples/" + stem + ".toml is missing or empty" );

		auto replaceLine = [ &text ]( std::string const &key,
		                              std::string const &line )
		{
			std::size_t at = text.find( "\n" + key );
			BOOST_TEST_REQUIRE( at != std::string::npos,
			                    "the shipped example has no " + key + " line" );
			++at;
			std::size_t const end = text.find( '\n', at );
			text.replace( at, end - at, line );
		};

		replaceLine( "Directory = ", "Directory = \"" + std::string( scratch() )
		                             + "\"" );
		replaceLine( "Prefix = ", "Prefix = \"" + prefix + "\"" );
		replaceLine( "RefinementLevels = ",
		             "RefinementLevels = " + std::to_string( refinement ) );

		std::string const path = std::string( scratch() ) + "/" + prefix
		                         + ".toml";
		std::ofstream out( path );
		out << text;
		return path;
	}

	/// What one run reports. `status` is the driver's exit code; the rest are
	/// parsed out of its standard output.
	struct Run
	{
		int status = -1;
		double psiAxis = 0.0;
		double normalisedFlux = 0.0;
		int elements = 0;
		int newton = 0;
		bool parsed = false;
	};

	/// The number after `marker` in `text`, or false if the marker is absent.
	bool after( std::string const &text, std::string const &marker,
	            double &value )
	{
		std::size_t const at = text.find( marker );
		if ( at == std::string::npos )
			return false;
		value = std::strtod( text.c_str() + at + marker.size(), nullptr );
		return true;
	}

	Run solve( std::string const &stem, int refinement,
	           std::string const &prefix )
	{
		std::string const config = configure( stem, refinement, prefix );
		std::string const log = std::string( scratch() ) + "/" + prefix + ".log";

		std::string const command = std::string( driver() ) + " " + config
		                            + " > " + log + " 2>&1";
		int const raw = std::system( command.c_str() );

		Run run;
		run.status = WIFEXITED( raw ) ? WEXITSTATUS( raw ) : -1;

		std::string const text = slurp( log );
		double value = 0.0;
		run.parsed = after( text, "psi_ax = ", run.psiAxis );
		if ( after( text, "normalised flux ", value ) )
			run.normalisedFlux = value;
		if ( after( text, "converged in ", value ) )
			run.newton = static_cast<int>( value );
		if ( after( text, " iterations on ", value ) )
			run.elements = static_cast<int>( value );
		return run;
	}

	/// The reference's own axis flux in MEQ's gauge -- freegs4e's psi_axis less
	/// the flux of the psi_n = 0.95 surface -- read out of the `-meta.json` the
	/// generator wrote beside each example. Quoted here rather than parsed
	/// because a test that reads its expectation out of a file the same script
	/// wrote is asserting that the script is self-consistent.
	struct Machine
	{
		char const *stem;
		double referencePsiAxis;
		char const *what;
	};

	Machine const circular
		{ "fixed-h-circular", 6.354433603318e-02,
		  "limited, kappa 1.01: the nearly-circular control" };
	Machine const testTokamak
		{ "fixed-a-testtokamak", 4.778846259077e-02,
		  "diverted, kappa 1.19, delta 0.32" };
}

// EVERY SHIPPED CASE STILL PARSES, STILL SOLVES, AND STILL REPORTS A MAGNETIC
// AXIS -- which is the whole of what "keeps an example green" means, and is the
// question the TODO entry left open about whether these should be ctests at
// all. They are cheap: no exterior coupling, no borders, no support sweep, one
// Dirichlet problem each.
BOOST_AUTO_TEST_CASE( theShippedMachineCasesSolveAndFindTheirAxis )
{
	std::system( ( "mkdir -p " + std::string( scratch() ) ).c_str() );

	Machine const machines[] = { circular, testTokamak };

	std::printf( "\n  THE SHIPPED MACHINE-SHAPED FIXED-BOUNDARY CASES\n" );
	std::printf( "    %-22s %8s %7s %14s %14s %10s\n",
	             "case", "elements", "Newton", "psi_ax", "freegs4e", "relative" );

	for ( Machine const &machine : machines )
	{
		Run const run = solve( machine.stem, 0, std::string( machine.stem ) );

		BOOST_TEST_REQUIRE( run.status == 0,
		                    std::string( "examples/" ) + machine.stem
		                    + ".toml did not solve: the driver exited "
		                    + std::to_string( run.status )
		                    + ". " + machine.what );
		BOOST_TEST_REQUIRE( run.parsed,
		                    "the driver printed no psi_ax, so [source] "
		                    "Normalised = true is no longer reaching the "
		                    "bordered Newton on this file" );

		double const relative = std::fabs( run.psiAxis
		                                   - machine.referencePsiAxis )
		                        / std::fabs( machine.referencePsiAxis );

		std::printf( "    %-22s %8d %7d %14.6e %14.6e %10.2e\n",
		             machine.stem, run.elements, run.newton, run.psiAxis,
		             machine.referencePsiAxis, relative );

		// psi_ax IS AN UNKNOWN OF THIS SOLVE, SO IT HAVING A VALUE AT ALL IS
		// A RESULT. The bordered Newton closes psi_ax - max psi_h = 0, and the
		// driver additionally reports the normalised flux at the LOCATED
		// O-point of q_h, which must be one. CLAUDE.md records that asserting
		// it under the located-axis constraint would be a tautology; here the
		// constraint is against the nodal peak, so this has teeth.
		BOOST_TEST( std::fabs( run.normalisedFlux - 1.0 ) < 1.0e-3,
		            "the axis is not where psi_ax says it is, so the profiles "
		            "were normalised by something that is not the axis flux" );

		// AGAINST AN INDEPENDENT CODE, at the coarsest mesh each file ships
		// with. The gate is loose on purpose -- it is a regression on the
		// EXAMPLE, and the convergence claim is the case below.
		BOOST_TEST( relative < 1.0e-3,
		            std::string( machine.stem )
		            + " no longer reproduces the equilibrium freegs4e solved, "
		              "which means the file, the profile tables or the curved "
		              "path have moved rather than that the mesh is coarse" );
	}
}

// AND REFINING MOVES IT TOWARD THE REFERENCE AND THEN STOPS THERE.
//
// THE SECOND HALF IS THE INTERESTING ONE AND IT IS WHY THIS IS NOT A RATE
// ASSERTION. What separates MEQ from freegs4e here is not MEQ's mesh: it is the
// MXH fit of Gamma, 2.7e-05 m on a = 0.3193, and freegs4e's own 257^2 grid.
// Neither moves when MEQ refines, so the agreement must improve and then FLOOR
// -- and a test asserting that it keeps improving would be asserting something
// false about the instrument rather than something true about the solver.
//
// RefinementLevels IS THE RIGHT KNOB and NR x NZ is not: bisection of triangles
// is exact, so each level is a mesh NESTED in the last and Gamma_h approaches
// Gamma through a nested sequence. Raising NR instead changes WHICH background
// elements fall inside Gamma, which is a different problem rather than a finer
// one -- see make_case.py's choose_mesh().
BOOST_AUTO_TEST_CASE( refiningTheMachineCaseMovesTowardTheIndependentAnswer )
{
	std::system( ( "mkdir -p " + std::string( scratch() ) ).c_str() );

	std::vector<double> gap;
	std::printf( "\n  REFINING %s TOWARD freegs4e\n", circular.stem );
	std::printf( "    %6s %8s %7s %14s %12s\n",
	             "levels", "elements", "Newton", "psi_ax", "rel to ref" );

	for ( int level = 0; level <= 2; ++level )
	{
		Run const run = solve( circular.stem, level,
		                       "refine" + std::to_string( level ) );
		BOOST_TEST_REQUIRE( run.status == 0 );
		BOOST_TEST_REQUIRE( run.parsed );

		double const relative = std::fabs( run.psiAxis
		                                   - circular.referencePsiAxis )
		                        / std::fabs( circular.referencePsiAxis );
		gap.push_back( relative );

		std::printf( "    %6d %8d %7d %14.6e %12.2e\n",
		             level, run.elements, run.newton, run.psiAxis, relative );
	}

	BOOST_TEST_REQUIRE( gap.size() == 3u );

	std::printf( "    first refinement closes the gap by %.1fx\n",
	             gap[ 0 ]/gap[ 1 ] );

	BOOST_TEST( gap[ 1 ] < 0.5*gap[ 0 ],
	            "one refinement did not halve the distance to freegs4e's own "
	            "answer, so either the curved path is not converging on this "
	            "geometry or the coarsest mesh was already at the reference's "
	            "floor -- print the sequence and look before relaxing this" );

	// AND THE FLOOR. The second refinement may not improve on the first,
	// because by then what is left is the MXH fit and freegs4e's grid; what it
	// must not do is get WORSE, which would say the extension is losing
	// something as Gamma_h approaches Gamma.
	BOOST_TEST( gap[ 2 ] < 2.0*gap[ 1 ],
	            "the second refinement moved AWAY from the reference by more "
	            "than the floor can explain" );

	BOOST_TEST( gap[ 2 ] < 1.0e-4,
	            "the refined answer does not reach the reference's own "
	            "accuracy, which is where this comparison bottoms out" );
}
