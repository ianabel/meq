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
 * `debug` AND NOT `cuda`, DELIBERATELY. `mfem::Device( "debug" )` has device
 * MEMORY semantics with host arithmetic, so it exercises the whole validity-
 * flag path -- which is where every one of the three defects above lived --
 * and needs no GPU, so this runs anywhere MFEM does. It is also the stricter
 * arm: it `mprotect`s the host page, so a raw read of a device-valid buffer
 * faults here and is merely a wrong number under CUDA. What it CANNOT see is
 * honest flags over a stale copy (M-187), and that is exactly what the
 * iteration-count assertion is for.
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
#include <string>
#include <vector>

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
		bool meshReady;        ///< carries [mesh.generate], so --mesh-ready
		char const *what;
	};

	/// What one run reported, as TEXT -- see the header for why.
	struct Run
	{
		int status = -1;
		double seconds = 0.0;
		std::string psiAxis;
		std::string psiBoundary;
		std::string xPoint;
		std::string newton;
	};

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

	Run solve( Case const &which, char const *device )
	{
		std::string const log = std::string( scratch() ) + "/"
		                        + which.stem + "-" + device + ".log";
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
		std::string command = "OMP_NUM_THREADS=1 " + std::string( driver() )
		                      + " --device " + device;
		if ( which.meshReady )
			command += " --mesh-ready";
		command += " examples/" + std::string( which.stem ) + ".toml > " + log
		           + " 2>&1";

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
		return run;
	}

	Case const rectangle
		{ "mhd-rectangle", false,
		  "nonlinear, fixed boundary, no border -- M-203's case" };
	Case const diiid
		{ "machine-f-diiid", true,
		  "bordered free boundary, four border rows -- M-187's case" };
}

// EVERY NONLINEAR SHIPPED EXAMPLE GIVES THE HOST'S ANSWER AND THE HOST'S WORK
// UNDER A CONFIGURED mfem::Device.
BOOST_AUTO_TEST_CASE( aDeviceChangesNeitherTheAnswerNorTheWork )
{
	std::system( ( "mkdir -p " + std::string( scratch() ) ).c_str() );

	Case const cases[] = { rectangle, diiid };

	std::printf( "\n  A DEVICE CHANGES NEITHER THE ANSWER NOR THE WORK\n" );
	std::printf( "    %-18s %-8s %6s %8s %16s %10s\n",
	             "case", "device", "exit", "s", "psi_ax", "Newton" );

	for ( Case const &which : cases )
	{
		Run const host = solve( which, "cpu" );
		Run const device = solve( which, "debug" );

		for ( auto const &pair : { std::make_pair( "cpu", &host ),
		                           std::make_pair( "debug", &device ) } )
			std::printf( "    %-18s %-8s %6d %8.1f %16s %10s\n",
			             which.stem, pair.first, pair.second->status,
			             pair.second->seconds,
			             pair.second->psiAxis.empty()
			                 ? "-" : pair.second->psiAxis.c_str(),
			             pair.second->newton.empty()
			                 ? "-" : pair.second->newton.c_str() );

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
	}
}
