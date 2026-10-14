/*
 * A DEVICE MUST NOT CHANGE THE ANSWER OR THE WORK, ON EVERY NONLINEAR EXAMPLE
 * AND NOT ON ONE.
 *
 * `DEVICE-PLAN.md` section 4's stage 4 states its acceptance as *"psi_ax,
 * psi_bnd, the X-point and the Newton count identical to the host ... on
 * machine-f-diiid.toml"*, and that row has carried a note since
 * MEASUREMENTS.md M-205 saying it should name every nonlinear example rather
 * than one machine. This is that test.
 *
 * WHY IT IS WORTH A FILE, given the suite already solves these cases. Three
 * defects reached a shipped example while every host test stayed green:
 *
 *   M-187  the answer right to every digit and the work FOUR TIMES the host's
 *   M-203  nine Newton steps against five, with a no-op interleaved
 *   M-206  psi_ax wrong by a FACTOR OF TWO, reported as converged
 *
 * None of them is visible to a host run, and the first is not visible to a
 * comparison of answers either -- which is why the iteration count is asserted
 * beside the numbers and not instead of them. **A device that gets the right
 * answer by doing twice the work has a defect**, and M-187 sat in this
 * project's measurements for weeks described as a property of the device path.
 *
 * `debug` IS THE PORTABLE ARM. `mfem::Device( "debug" )` has device MEMORY
 * semantics with host arithmetic, so it exercises the whole validity-flag path
 * -- which is where every one of the three defects above lived -- and needs no
 * GPU, so this runs anywhere MFEM does. It is also the stricter arm for one
 * class: it `mprotect`s the host page, so a raw read of a device-valid buffer
 * faults here and is merely a wrong number under CUDA.
 *
 * **AND `cuda` IS A THIRD ARM, BECAUSE `debug` IS BLIND TO THE CLASS THAT COST
 * M-206 A FACTOR OF TWO AND THE ITERATION COUNT DOES NOT COVER IT EITHER.**
 * This header used to say the count was what the stale-copy class needed. It
 * is what **M-187**'s needs -- an answer right to every digit on four times
 * the work. **M-206 defeats both**: a 43-solve outer loop ran to completion
 * under `debug` without faulting once, `psi_ax` came back wrong by a factor of
 * two, and the inner bordered solve took the SAME four Newton steps on both
 * arms. Honest flags over a stale copy fault nowhere and change no count. What
 * found it was `cpu` against `cuda`, compared digit for digit.
 *
 * So the cuda arm runs where the build and the machine have one, and is
 * **skipped loudly rather than failed** where they do not -- a build without
 * CUDA is an absence, not a disagreement. It is nearly free: `debug` is the
 * expensive arm at some hundreds of seconds a case and a card run is seconds.
 * **It is added BEFORE item 6 on purpose**: `DEVICE-PLAN.md` section 3.5 ends
 * the digit-for-digit instrument when cuDSS becomes the device default, so
 * this is the last point at which a CUDA arm can be pinned exactly rather than
 * to a tolerance.
 *
 * THE COMPARISON IS ON THE PRINTED TEXT, not on a parsed double within a
 * tolerance. "Identical to every printed digit" is what the acceptance says
 * and it is the only form with no number in it for a later session to widen.
 * Where a device legitimately cannot match -- it does not, on any case here --
 * the right repair is to say so per case, not to introduce a tolerance.
 */

#define BOOST_TEST_MODULE DeviceParity
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "mfem.hpp"

#include "RoundOffFloor.hpp"

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iterator>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "SourceTree.hpp"

namespace
{
	/// Where the driver binary is; CMake passes it, for the reason
	/// DriverAcceptance records.
	char const *driver()
	{
		return MEQ_DRIVER_PATH;
	}

	char const *scratch()
	{
		return "device-parity-scratch";
	}

	std::string slurp( std::string const &path )
	{
		std::ifstream stream( path );
		return std::string( ( std::istreambuf_iterator<char>( stream ) ),
		                      std::istreambuf_iterator<char>() );
	}

	/// A shipped example to run on both devices.
	struct Case
	{
		char const *stem;      ///< examples/<stem>.toml
		bool meshReady;        ///< carries [mesh.generate]: made here, then --mesh-ready
		char const *what;
		/// Whether to pay for the `debug` arm. It is the expensive one -- it
		/// `mprotect`s every host page -- and a case added for a defect the
		/// card arm catches need not buy it. Defaults on.
		bool debugArm = true;
		/// A key to REMOVE from the shipped file before running it, or null to
		/// run the file as shipped. The copy is written into the scratch
		/// directory; see the moving-support case for why one exists.
		char const *strip = nullptr;
	};

	/// The file to hand the driver: the shipped one, or a copy in the scratch
	/// directory with `which.strip`'s line removed and, where the file
	/// generates its mesh, that mesh redirected into the scratch directory and
	/// made there -- once per process, since every arm reads the same one.
	///
	/// THE MESH IS MADE HERE BECAUSE A GENERATED [mesh] File IS THE RUN'S
	/// PRODUCT, written relative to the working directory, so nothing ships it;
	/// and into this test's own scratch because the working directory is shared
	/// with every other test, InverseSolve's meq-run among them, which makes the
	/// same file under the same name. Made the way meq-run makes it: the
	/// command `meq --mesh-command` prints, run by tools/mesh/halfdisc.py.
	std::string configFor( Case const &which )
	{
		std::string const relative = "examples/" + std::string( which.stem ) + ".toml";
		std::string const shipped = meq::test::sourcePath( relative );
		if ( which.strip == nullptr && !which.meshReady )
			return shipped;

		static std::map<std::string, std::string> prepared;
		auto const found = prepared.find( which.stem );
		if ( found != prepared.end() )
			return found->second;

		std::string const copy = std::string( scratch() ) + "/" + which.stem
		                         + ( which.strip ? std::string( "-without-" ) + which.strip
		                                         : std::string() )
		                         + ".toml";
		std::istringstream in( meq::test::shippedConfig( relative ) );
		std::ostringstream text;
		std::string line;
		std::string table;
		bool removed = false;
		bool redirected = false;
		while ( std::getline( in, line ) )
		{
			std::size_t const first = line.find_first_not_of( " \t" );
			if ( first != std::string::npos && line[ first ] == '[' )
				table = line.substr( first );
			if ( which.strip != nullptr && first != std::string::npos
			     && line.compare( first, std::string( which.strip ).size(),
			                      which.strip ) == 0 )
			{
				removed = true;
				continue;
			}
			if ( which.meshReady && table == "[mesh]" && first != std::string::npos
			     && line.compare( first, 7, "File = " ) == 0 )
			{
				line = "File = \"" + std::string( scratch() ) + "/" + which.stem + ".msh\"";
				redirected = true;
			}
			text << line << '\n';
		}
		if ( which.strip != nullptr && !removed )
			throw std::runtime_error(
				shipped + " has no " + which.strip + " line to remove, so this "
				"case would run the shipped configuration and test nothing it "
				"was added for" );
		if ( which.meshReady && !redirected )
			throw std::runtime_error( shipped + " has no [mesh] File line to redirect" );
		{
			std::ofstream out( copy );
			out << text.str();
		}

		if ( which.meshReady )
		{
			std::string const commandFile = copy + ".mesh-command";
			std::string const commandLog = copy + ".mesh-command-log";
			if ( std::system( ( std::string( driver() ) + " --mesh-command " + copy
			                    + " > " + commandFile + " 2> " + commandLog ).c_str() ) != 0 )
				throw std::runtime_error( "meq --mesh-command " + copy + " failed: "
				                          + slurp( commandLog ) );
			std::istringstream commandText( slurp( commandFile ) );
			std::string command;
			std::getline( commandText, command );
			if ( command.rfind( "halfdisc ", 0 ) != 0 )
				throw std::runtime_error( "meq --mesh-command " + copy + " printed '"
				                          + command + "', not a halfdisc command" );
			std::string const meshLog = copy + ".mesh-log";
			if ( std::system( ( "python3 " + meq::test::sourcePath( "tools/mesh/halfdisc.py" )
			                    + " " + command.substr( 9 ) + " > " + meshLog
			                    + " 2>&1" ).c_str() ) != 0 )
				throw std::runtime_error( "the mesh generator failed for " + copy + ": "
				                          + slurp( meshLog ) );
		}

		prepared.emplace( which.stem, copy );
		return copy;
	}

	/// What one run reported, as TEXT -- see the header for why.
	struct Run
	{
		int status = -1;
		double seconds = 0.0;
		std::string psiAxis;
		std::string psiBoundary;
		std::string xPoint;
		std::string newton;

		/// Whether MFEM printed its own device configuration naming this
		/// backend. **A build carrying CUDA on a machine with no card does not
		/// get that line**, and the difference between *the device was never
		/// configured* and *the device gave a different answer* is the whole
		/// of what makes a skip honest. Asked of MFEM's own output rather than
		/// inferred from an exit code, because an abort and a refusal share one.
		bool configured = false;
	};

	/// Whether this BUILD can reach CUDA at all.
	///
	/// **THIS TRANSLATION UNIT INCLUDES `mfem.hpp`, WHICH IS WHAT MAKES THE
	/// TEST LEGITIMATE.** `MachineFixedBoundary.cpp` drives the shipped binary,
	/// needs nothing of the library, included no MFEM header -- and its
	/// `#ifdef MFEM_USE_SINGLE` compiled, ran, and silently took the other
	/// branch, because a precision macro in a TU that has not included MFEM is
	/// not false but UNDEFINED. `mfem/config/config.hpp` is what defines these.
	bool cudaIsBuilt()
	{
#ifdef MFEM_USE_CUDA
		return true;
#else
		return false;
#endif
	}

	/// The characters following @a marker up to the end of the number, or an
	/// empty string when the marker is absent. Absence is itself comparable:
	/// a field one device prints and the other does not is a difference.
	std::string after( std::string const &text, std::string const &marker )
	{
		std::size_t const at = text.find( marker );
		if ( at == std::string::npos )
			return {};
		std::size_t i = at + marker.size();
		std::size_t const begin = i;
		while ( i < text.size()
		        && ( std::isdigit( static_cast<unsigned char>( text[ i ] ) )
		             || text[ i ] == '.' || text[ i ] == '-'
		             || text[ i ] == '+' || text[ i ] == 'e' ) )
			++i;
		return text.substr( begin, i - begin );
	}

	/// @param forced run with `MEQ_AXIS_DEVICE=force`, so every host reader of
	///        the iterate takes its device route whatever the iterate's flags
	///        say. See the fourth arm below.
	Run solve( Case const &which, char const *device, bool forced = false )
	{
		std::string const log = std::string( scratch() ) + "/"
		                        + which.stem + "-" + device
		                        + ( forced ? "-forced" : "" ) + ".log";
		/*
		 * `OMP_NUM_THREADS=1` ON BOTH ARMS, AND IT IS NOT A CONVENIENCE.
		 * A Device and a host OpenMP region that touches MFEM memory are
		 * mutually exclusive -- MEASUREMENTS.md M-188 -- so the driver REFUSES
		 * `AssemblyMode = "threaded"` above one thread under `--device`, and
		 * several shipped examples ask for exactly that. Without this the
		 * device arm exits 1 and the test reports a parity failure when what
		 * it has actually met is a refusal working correctly. Pinned on the
		 * HOST arm too, so the two sides differ in the device and in nothing
		 * else.
		 */
		std::string command = std::string( forced ? "MEQ_AXIS_DEVICE=force " : "" )
		                      + "OMP_NUM_THREADS=1 " + std::string( driver() )
		                      + " --device " + device;
		if ( which.meshReady )
			command += " --mesh-ready";
		command += " " + configFor( which ) + " > " + log + " 2>&1";

		auto const started = std::chrono::steady_clock::now();
		int const raw = std::system( command.c_str() );
		auto const ended = std::chrono::steady_clock::now();

		Run run;
		run.status = WIFEXITED( raw ) ? WEXITSTATUS( raw ) : -1;
		run.seconds = std::chrono::duration<double>( ended - started ).count();

		std::string const text = slurp( log );
		run.psiAxis = after( text, "psi_ax = " );
		run.psiBoundary = after( text, "psi_bnd = " );
		run.xPoint = after( text, "at the X-point ( " );
		run.newton = after( text, "converged in " );
		run.configured =
			text.find( std::string( "Device configuration: " ) + device )
			!= std::string::npos;
		return run;
	}

	Case const rectangle
		{ "mhd-rectangle", false,
		  "nonlinear, fixed boundary, no border -- M-203's case" };
	Case const diiid
		{ "machine-f-diiid", true,
		  "bordered free boundary, four border rows -- M-187's case" };
	/*
	 * A MOVING PLASMA SUPPORT WITH XP-1's FILL, AND ITS ABSENCE FROM THIS LIST
	 * IS HOW A WRONG ANSWER SURVIVED. Under any Device the integrand's element
	 * flags were an upload refreshed only when the operator was rebuilt, while
	 * the host's followed refreshPlasmaComponent() at every residual -- so
	 * `--device cuda` parted from the host at the first Newton step and
	 * settled with its axis inside a coil, at `HEAD` as well as after.
	 * MEASUREMENTS.md M-239. The two cases above freeze the support within a
	 * solve, which is why neither could see it.
	 *
	 * **RUN WITH `PlasmaSupportSweeps` STRIPPED, BECAUSE THE SHIPPED FILE NOW
	 * FREEZES IT** -- that is what puts it on the benchmarked branch (M-239) --
	 * and no other shipped example moves a FILLED support: `limited-tokamak`
	 * moves one without the fill and agreed under CUDA with the stale flags.
	 * So the copy is the only configuration in the tree that reaches the path,
	 * and the strip throws if the key it removes ever disappears.
	 */
	Case const moving
		{ "diverted-tokamak", false,
		  "moving, filled plasma support -- M-239's case", false,
		  "PlasmaSupportSweeps" };
}

// EVERY NONLINEAR SHIPPED EXAMPLE GIVES THE HOST'S ANSWER AND THE HOST'S WORK
// UNDER A CONFIGURED mfem::Device.
BOOST_AUTO_TEST_CASE( aDeviceChangesNeitherTheAnswerNorTheWork )
{
	std::system( ( "mkdir -p " + std::string( scratch() ) ).c_str() );

	Case const cases[] = { rectangle, diiid, moving };

	std::printf( "\n  A DEVICE CHANGES NEITHER THE ANSWER NOR THE WORK\n" );
	std::printf( "    %-18s %-8s %6s %8s %16s %10s\n",
	             "case", "device", "exit", "s", "psi_ax", "Newton" );

	for ( Case const &which : cases )
	{
		Run const host = solve( which, "cpu" );
		Run const device = which.debugArm ? solve( which, "debug" ) : host;

		/*
		 * THE CARD ARM. Run only where the build has CUDA; `configured` then
		 * says whether the MACHINE did too, read off MFEM's own output rather
		 * than guessed from an exit code. See the header for why `debug` alone
		 * is not enough and why this is added before item 6 rather than after.
		 */
		Run const card = cudaIsBuilt() ? solve( which, "cuda" ) : Run{};

		for ( auto const &pair : { std::make_pair( "cpu", &host ),
		                           std::make_pair( "debug", &device ),
		                           std::make_pair( "cuda", &card ) } )
		{
			if ( std::string( pair.first ) == "cuda" && !card.configured )
				continue;
			if ( std::string( pair.first ) == "debug" && !which.debugArm )
			{
				std::printf( "    %-18s %-8s %s\n", which.stem, "debug",
				             "SKIPPED -- not paid for on this case; the cuda "
				             "arm is the one it was added for" );
				continue;
			}
			std::printf( "    %-18s %-8s %6d %8.1f %16s %10s\n",
			             which.stem, pair.first, pair.second->status,
			             pair.second->seconds,
			             pair.second->psiAxis.empty()
			                 ? "-" : pair.second->psiAxis.c_str(),
			             pair.second->newton.empty()
			                 ? "-" : pair.second->newton.c_str() );
		}
		if ( !card.configured )
			std::printf( "    %-18s %-8s %s\n", which.stem, "cuda",
			             cudaIsBuilt()
			                 ? "SKIPPED -- built with CUDA, no device configured"
			                 : "SKIPPED -- this build has no CUDA" );

		BOOST_TEST_REQUIRE( host.status == 0,
		                    std::string( which.stem )
		                        + " does not solve on the host at all, so "
		                          "there is nothing to compare a device "
		                          "against" );
		/*
		 * A RED HERE MUST SAY WHAT FIXES IT. `CLAUDE.md`'s testing stance:
		 * *"a failing test must say what fixes it, because a red suite nobody
		 * can action is one people learn to ignore"* -- and the message below
		 * is written for a reader who has just been handed a red they did not
		 * cause, because the defect this case has caught twice was upstream's
		 * both times.
		 *
		 * THE CASE IS GREEN. M-214's ReconstructTotalFlux() defect is fixed at
		 * upstream's `d720f02d62` and M-207's `mhd-rectangle` at `ce744b84de`,
		 * so both arms exit 0 and the whole binary passes. The named diagnosis
		 * stays in the message anyway: it is the cheapest thing to check first
		 * if the exit status ever moves again, and LocalFactorMode is still the
		 * one key that separates a device defect in the SOLVE from one in the
		 * post-processing -- `"serial"` there leaves the answer untouched and
		 * moves only which route reads the element blocks.
		 */
		BOOST_TEST( device.status == 0,
		            std::string( which.stem ) + " exits "
		                + std::to_string( device.status )
		                + " under --device debug and 0 on the host. CHECK "
		                  "LocalFactorMode FIRST: a fault of exactly this shape "
		                  "was DarcyHybridization::ReconstructTotalFlux() "
		                  "reading the element-local blocks on the host after "
		                  "LocalFactorMode::Batched left them on the device, and "
		                  "setting LocalFactorMode = \"serial\" in the file made "
		                  "the same run exit 0 with the host's answer to every "
		                  "digit -- which is what separates a device fault in "
		                  "the SOLVE from one in the post-processing. Fixed "
		                  "upstream at d720f02d62; MEASUREMENTS.md M-214. If "
		                  "that key does not move it, this is new -- walk it "
		                  "with MEQ_STAGE_CHECKSUM=1" );
		/*
		 * THE NEWTON COUNT IS THE ONE FIELD EVERY CASE PRINTS, and on a
		 * fixed-boundary case it is the ONLY one: `psi_ax` is a bordered
		 * constraint's report and `mhd-rectangle` has no border. That is not
		 * a weakness here -- M-203's defect on exactly that case WAS the
		 * count, 9 against 5, with the answer agreeing to 1.05e-14.
		 */
		BOOST_TEST_REQUIRE( !host.newton.empty(),
		                    std::string( which.stem )
		                        + " printed no Newton count on the host, so "
		                          "there is no field to compare" );

		/*
		 * THE ANSWER, TO EVERY PRINTED DIGIT -- AND THE NUMBER OF DIGITS A
		 * FORMAT HAS IS A PROPERTY OF THE FORMAT.
		 *
		 * The strings are compared first and that is the whole test in double:
		 * *"--device cuda reproduces --device cpu to every printed digit"* is
		 * the instrument M-79, M-130 and M-187 were root-caused with, and it
		 * has to stay exact there. Under MFEM_PRECISION = single the driver
		 * prints seven significant figures of a `real_t`, which IS `eps`, so
		 * the device's different summation order moves the last one and nothing
		 * else: measured on machine-f-diiid, psi_ax 3.759859e-01 against
		 * 3.759861e-01, psi_bnd 7.081396e-02 against 7.081404e-02 and the
		 * X-point 1.200931 against 1.200930 -- five parts in ten million, one
		 * ulp of the printed field.
		 *
		 * exactAgreementFloor() is EXACTLY ZERO in double, so the fallback
		 * below cannot admit anything there that the strings did not already
		 * admit. What it costs in single is the last printed digit; what it
		 * keeps is everything a wrong kernel does, which M-206 measured at a
		 * FACTOR OF TWO and M-79 at an identically zero potential. M-221.
		 */
		auto printedDigitsAgree = []( std::string const &a, std::string const &b )
		{
			if ( a == b )
				return true;
			try
			{
				double const x = std::stod( a );
				double const y = std::stod( b );
				return std::abs( x - y )
				       <= meq::tests::exactAgreementFloor(
					       std::max( std::abs( x ), std::abs( y ) ) );
			}
			catch ( std::exception const & )
			{
				return false;
			}
		};

		BOOST_TEST( printedDigitsAgree( device.psiAxis, host.psiAxis ),
		            std::string( which.stem ) + ": psi_ax is "
		                + device.psiAxis + " under a device against "
		                + host.psiAxis + " on the host" );
		BOOST_TEST( printedDigitsAgree( device.psiBoundary, host.psiBoundary ),
		            std::string( which.stem ) + ": psi_bnd is "
		                + device.psiBoundary + " against " + host.psiBoundary );
		BOOST_TEST( printedDigitsAgree( device.xPoint, host.xPoint ),
		            std::string( which.stem ) + ": the X-point is "
		                + device.xPoint + " against " + host.xPoint );

		// AND THE WORK, which is the half a comparison of answers cannot see:
		// M-187 agreed on every printed number and did four times as much.
		BOOST_TEST( device.newton == host.newton,
		            std::string( which.stem ) + ": " + device.newton
		                + " Newton iterations under a device against "
		                + host.newton + " on the host -- the answer may be "
		                  "right and the work is not (MEASUREMENTS.md M-187)" );

		/*
		 * AND THE SAME FOUR AGAINST THE CARD, WHICH IS A DIFFERENT QUESTION
		 * AND NOT A REPEAT.
		 *
		 * `debug` and `cuda` disagree about exactly one thing and it is the
		 * one that matters here: `debug` protects the host page, so a stale
		 * copy read raw is a FAULT, while under CUDA nothing is protected and
		 * the same read returns the wrong number quietly. M-206 is the worked
		 * example -- a whole outer loop clean under `debug`, `psi_ax` wrong by
		 * a factor of two under CUDA, the same Newton count on both. So a
		 * green `debug` arm is not evidence about this one.
		 *
		 * Skipped and not failed where no device was configured: *an absent
		 * GPU is an absence, not a disagreement*, and `configured` is read off
		 * MFEM's own line so the two cannot be confused.
		 */
		if ( !card.configured )
			continue;

		BOOST_TEST( card.status == 0,
		            std::string( which.stem ) + " exits "
		                + std::to_string( card.status )
		                + " under --device cuda and 0 on the host" );
		BOOST_TEST( printedDigitsAgree( card.psiAxis, host.psiAxis ),
		            std::string( which.stem ) + ": psi_ax is "
		                + card.psiAxis + " under CUDA against "
		                + host.psiAxis + " on the host. THE debug ARM BEING "
		                  "GREEN SAYS NOTHING HERE -- debug faults on a raw "
		                  "read of a device-valid buffer and is blind to "
		                  "honest flags over a STALE copy, which is what cost "
		                  "q-driven a factor of two (MEASUREMENTS.md M-206). "
		                  "Walk it with MEQ_STAGE_CHECKSUM=1 on both arms and "
		                  "diff: the first stage that disagrees names it" );
		BOOST_TEST( printedDigitsAgree( card.psiBoundary, host.psiBoundary ),
		            std::string( which.stem ) + ": psi_bnd is "
		                + card.psiBoundary + " under CUDA against "
		                + host.psiBoundary );
		BOOST_TEST( printedDigitsAgree( card.xPoint, host.xPoint ),
		            std::string( which.stem ) + ": the X-point is "
		                + card.xPoint + " under CUDA against " + host.xPoint );
		BOOST_TEST( card.newton == host.newton,
		            std::string( which.stem ) + ": " + card.newton
		                + " Newton iterations under CUDA against "
		                + host.newton + " on the host" );

		/*
		 * AND A FOURTH ARM, WITH EVERY DEVICE ROUTE FORCED OPEN.
		 *
		 * The arms above run the routes ResidentIterate.hpp and AxisSearch.hpp
		 * built, and NOT ONE OF THEM IS TAKEN: each is gated on the iterate being
		 * device-valid and host-invalid, which nothing yet arranges, so every
		 * shipped case reads `0 taken`. A route nobody takes is a route nobody
		 * tests, and the day residency arrives they all switch on at once.
		 *
		 * `MEQ_AXIS_DEVICE=force` takes them now, on the same file, and the same
		 * four fields are compared with the host's as printed text. It is the
		 * acceptance for the ten host readers of the iterate being ported --
		 * the plasma loops, the transmission and psi_bnd rows, the cold state --
		 * and it costs a card run, which is seconds. The `debug` arm is not
		 * repeated forced: it is the expensive one, and a forced run under it is
		 * recorded in MEASUREMENTS.md rather than paid for on every ctest.
		 */
		Run const forced = solve( which, "cuda", true );
		std::printf( "    %-18s %-8s %6d %8.1f %16s %10s\n", which.stem,
		             "cuda+f", forced.status, forced.seconds,
		             forced.psiAxis.empty() ? "-" : forced.psiAxis.c_str(),
		             forced.newton.empty() ? "-" : forced.newton.c_str() );

		BOOST_TEST( forced.status == 0,
		            std::string( which.stem ) + " exits "
		                + std::to_string( forced.status )
		                + " under --device cuda with MEQ_AXIS_DEVICE=force" );
		BOOST_TEST( printedDigitsAgree( forced.psiAxis, host.psiAxis ),
		            std::string( which.stem ) + ": psi_ax is "
		                + forced.psiAxis + " with the device routes forced "
		                  "against " + host.psiAxis + " on the host. Every "
		                  "route reproduces a host loop operation for operation "
		                  "on host-computed radii and weights; a digit here is a "
		                  "kernel that does not" );
		BOOST_TEST( printedDigitsAgree( forced.psiBoundary, host.psiBoundary ),
		            std::string( which.stem ) + ": psi_bnd is "
		                + forced.psiBoundary + " with the device routes forced "
		                  "against " + host.psiBoundary );
		BOOST_TEST( printedDigitsAgree( forced.xPoint, host.xPoint ),
		            std::string( which.stem ) + ": the X-point is "
		                + forced.xPoint + " with the device routes forced "
		                  "against " + host.xPoint );
		BOOST_TEST( forced.newton == host.newton,
		            std::string( which.stem ) + ": " + forced.newton
		                + " Newton iterations with the device routes forced "
		                  "against " + host.newton + " on the host" );
	}
}
