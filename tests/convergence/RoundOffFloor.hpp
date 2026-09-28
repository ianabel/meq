#ifndef MEQ_TESTS_ROUNDOFFFLOOR_HPP
#define MEQ_TESTS_ROUNDOFFFLOOR_HPP

/*
 * THE ROUND-OFF FLOOR, AND WHY IT IS ITS OWN HEADER.
 *
 * An absolute L2 error ceiling is a relative one with the scale left out: every
 * ceiling in this directory was set at about three times a measured
 * double-precision error, so each is `3 C h^(k+1) |psi|` with the `|psi|` folded
 * into the digits. That is a good assertion while discretisation limits the
 * error and a meaningless one once the FORMAT does -- `SolovievConvergence`
 * reads 1.08e-05 at `k = 2` and 7.60e-06 at `k = 3` against ceilings that fall
 * by 340 between them.
 *
 * **THIS IS A HEADER RATHER THAN PART OF `ConvergenceHarness.hpp` BECAUSE FOUR
 * OF THE FIVE STUDIES THAT NEED IT DO NOT INCLUDE THE HARNESS.**
 * `SolovievConvergence`, `ExtensionConvergence`, `EstimatorConvergence` and
 * `RotatingSourceConvergence` each carry their own `Measurement`, their own
 * `rate()`, their own `printTable()` and -- until M-218 -- their own copy of the
 * same six-line ceiling assertion. Putting the floor where only the harness
 * could reach it would have floored one of four call sites and left three
 * reporting a format as a method, which is M-165 and M-189 repeating for a third
 * time. Small, dependency-free, and included by all five.
 *
 * M-218.
 */

#include <algorithm>
#include <limits>
#include <sstream>
#include <string>

#include <boost/test/unit_test.hpp>

#include "mfem.hpp"

namespace meq
{
namespace tests
{

	/**
	 * The smallest L2 error a solve in this build's `real_t` can be expected to
	 * reach on a field of size @a norm.
	 *
	 * **AN ABSOLUTE L2 CEILING IS A RELATIVE ONE WITH THE SCALE LEFT OUT.**
	 * Every ceiling in this directory was set at about three times a measured
	 * double-precision error, so each is `3 C h^(k+1) |psi|` with the `|psi|`
	 * folded into the digits -- a good assertion while discretisation limits the
	 * error and a meaningless one once the FORMAT does. `SolovievConvergence`
	 * reads 1.04e-05 at `k = 2` and 7.75e-06 at `k = 3` against ceilings that
	 * fall by 340 between them.
	 *
	 * **THE CONSTANT IS 2048 AND IT IS MEASURED AGAINST THE DOUBLE BUILD**, on
	 * the same working tree, which is the only control that separates a format's
	 * floor from a method's error. The worst `|| error || / ( eps || exact || )`
	 * over the sequences in this directory is **1412** -- `Solov'ev NSTX`, the
	 * flux, `k = 3` on the finest mesh -- so the margin is 1.45. Thin for a
	 * discriminator and ample for a ceiling: what a ceiling exists to catch is a
	 * solution wrong by a constant factor, a sign or a convention, and those are
	 * `O( || exact || )` errors, seven orders above this.
	 *
	 * **AND IT CANNOT REACH THE DOUBLE BUILD**, which was checked rather than
	 * argued: `2048 eps_double || q ||` is `1.1e-12` on the largest `|| q ||` in
	 * the directory (2.398), and the tightest historical ceiling anywhere in it
	 * is `6.6e-10`. There is no sequence in which the `std::max` selects the
	 * floor in double. M-218.
	 */
	inline double roundOffFloor( double norm )
	{
		return 2048.0*norm*static_cast<double>(
			std::numeric_limits<mfem::real_t>::epsilon() );
	}

	/**
	 * Whether this build's `real_t` can measure a convergence rate of
	 * @a order + 1 at all.
	 *
	 * **IN SINGLE PRECISION IT CANNOT, ABOVE `k = 1`, AND THAT IS A
	 * MEASUREMENT.** Against the double build on the same working tree -- the
	 * only control that separates a format's floor from a method's error --
	 * `Solov'ev NSTX` reads, in the flux, single divided by double:
	 *
	 *     k = 1   n=4 1.00x   n=8 1.00x   n=16 1.00x   n=32   1.09x
	 *     k = 2   n=4 1.00x   n=8 1.05x   n=16  6.9x   n=32    203x
	 *     k = 3   n=4 1.50x   n=8  36x    n=16 1793x   n=32  61000x
	 *
	 * At `k = 1` every mesh is clean and every pair reports 1.87 to 1.99. At
	 * `k = 2` exactly ONE pair of three survives. At `k = 3` the COARSEST mesh
	 * of the sequence is already 1.5x contaminated, so there is no window at
	 * all: coarsening leaves the asymptotic regime and refining goes deeper into
	 * round-off.
	 *
	 * **AND THE TAIL BENDS BEFORE IT FLATTENS, WHICH IS WHY NO PER-PAIR TEST
	 * WORKS.** Three were tried and measured. A floor of `C eps || q ||` times
	 * the dof count puts `k = 1, n = 32` -- 9% contaminated and reporting a good
	 * 1.873 -- at **0.277** and `k = 3, n = 8` -- 36x contaminated -- at
	 * **0.249**, and no constant separates them. A constant-free "the error must
	 * at least halve" is a weak necessary condition and correctly drops the
	 * flattened rows, but `MaschkePerrinConvergence` at `k = 3` has a pair that
	 * halves comfortably and still reports **3.03** against a wanted 3.85 --
	 * partially contaminated, genuinely converging, and not at `k+1`. Asserting
	 * on the best pair does not save it either, because that IS the best pair.
	 *
	 * **SO THE SKIP IS BY DEGREE AND IT IS HONEST RATHER THAN TUNED.** `k+1` at
	 * `k >= 2` is a double-build acceptance; a single build verifies it at
	 * `k = 1` and, above that, verifies the ABSOLUTE error against a floored
	 * ceiling instead. Nothing is left unable to fail: a solution wrong by a
	 * constant factor, a sign or a convention is `O( || exact || )` and is
	 * exactly what the ceiling catches, with seven orders of margin.
	 *
	 * **THE ALTERNATIVE WOULD HAVE BEEN TO WIDEN `rateSlack` UNTIL THE TABLES
	 * PASSED**, which is the vacuity `CLAUDE.md`'s testing stance forbids: it
	 * would weaken the assertion on every mesh in every build, including the
	 * ones still measuring the method, to accommodate meshes that are not.
	 * M-218.
	 */
	inline bool rateIsMeasurable( int order )
	{
#ifdef MFEM_USE_SINGLE
		return order <= 1;
#else
		(void)order;
		return true;
#endif
	}

	/// The message a caller prints when [[rateIsMeasurable]] refuses, so that the
	/// weakening is visible in the output rather than inferred from a gap.
	inline std::string rateNotMeasurableNotice( int order )
	{
		std::ostringstream text;
		text << "  *** k = " << order << ": SKIPPED the rate assertion. This build "
		        "is MFEM_PRECISION = single, in which the discretisation error at "
		        "k >= 2 falls under the format's own floor before the sequence ends "
		        "-- measured against the double build at 6.9x contamination by the "
		        "third mesh at k = 2 and 1.5x on the FIRST at k = 3. The absolute "
		        "ceilings below are still asserted. M-218";
		return text.str();
	}

	inline void checkCeilings( int order, double errorPsi, double normPsi,
	                           double psiCeiling, double errorFlux,
	                           double normFlux, double fluxCeiling )
	{
		double const psiBound = std::max( psiCeiling, roundOffFloor( normPsi ) );
		double const fluxBound = std::max( fluxCeiling, roundOffFloor( normFlux ) );

		/*
		 * SAY SO WHEN THE FLOOR IS WHAT CARRIED THE CASE. A ceiling lifted by
		 * the format is a weaker assertion than the one the caller wrote, and a
		 * weakening nobody can see in the output is the kind that outlives its
		 * reason. This prints on a PASS, which is the only run in which anybody
		 * would otherwise not look.
		 */
		if ( errorPsi >= psiCeiling && errorPsi < psiBound )
			BOOST_TEST_MESSAGE(
				"  k = " << order << ": psi is inside the ROUND-OFF FLOOR " << psiBound
				<< " and outside the ceiling " << psiCeiling << " this caller set"
				<< " -- format-limited, not method-limited" );
		if ( errorFlux >= fluxCeiling && errorFlux < fluxBound )
			BOOST_TEST_MESSAGE(
				"  k = " << order << ": q is inside the ROUND-OFF FLOOR " << fluxBound
				<< " and outside the ceiling " << fluxCeiling << " this caller set"
				<< " -- format-limited, not method-limited" );

		BOOST_TEST( errorPsi < psiBound,
		            "k = " << order << ": L2 error in psi is " << errorPsi
		            << ", above the ceiling " << psiBound << " ( set at " << psiCeiling
		            << ", floored at round-off on a || psi || of " << normPsi << " )" );
		BOOST_TEST( errorFlux < fluxBound,
		            "k = " << order << ": L2 error in q is " << errorFlux
		            << ", above the ceiling " << fluxBound << " ( set at " << fluxCeiling
		            << ", floored at round-off on a || q || of " << normFlux << " )" );
	}

}
}

#endif
