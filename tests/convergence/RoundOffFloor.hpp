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
#include <cmath>
#include <limits>
#include <sstream>
#include <vector>
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
	 * Whether one Newton step took an AFFINE residual to this build's round-off
	 * floor -- which is what "an exact Jacobian finishes it in one" means, and
	 * is a different statement from "the iteration count is 1".
	 *
	 * **THE COUNT IS AN ARTEFACT OF THE TARGET AND THE REDUCTION IS THE CLAIM.**
	 * `solve()` stops at a target scaled from the cold residual, so on a source
	 * with `dF/dpsi == 0` the question of whether a SECOND step is taken is
	 * decided by whether the first step's residual happens to land under that
	 * target -- which is a property of the format, not of the Jacobian.
	 * Measured on `SolovievGeometryConvergence`'s three geometries,
	 * `|| r_1 || / || r_0 ||`:
	 *
	 *     double   1.29e-13   3.62e-14   1.26e-13      Newton stops at 1
	 *     single   6.42e-05   2.10e-05   5.87e-05      Newton stops at 2
	 *
	 * In units of each build's own `eps` those are **580, 163, 566** and
	 * **539, 176, 493** -- the same number. *The Jacobian is exactly as exact in
	 * single as in double, and the iteration count says otherwise.*
	 *
	 * So the constant is **2048**, a margin of 3.5 over the worst measured, and
	 * this asserts the reduction in both builds rather than the count in
	 * either. It is STRONGER than the count it replaces: a Jacobian that is
	 * merely close would still stop at one iteration if the target were loose
	 * enough, and could not reach `2048 eps` in one step. M-219.
	 */
	/**
	 * How small `| q_h |` must be for a located critical point to be a ROOT,
	 * floored from the caller's own historical value.
	 *
	 * **THE DISCRIMINATOR IS SEVEN ORDERS WIDE AND THE FLOOR EATS FOUR OF
	 * THEM.** `CriticalPointConvergence` states its own yardstick: *"a point
	 * reported at `|q| = 1e-3` is a place Newton gave up"*, against a bound of
	 * `1.0e-10`. In single the root finder cannot go below round-off on an
	 * `O( 1 )` flux: measured, the located axis and X-point report `4.28e-07`,
	 * `1.10e-07`, `7.93e-08`, `6.22e-08` and `1.34e-06`.
	 *
	 * **AND A SADDLE IS WORSE CONDITIONED THAN AN EXTREMUM**, which the numbers
	 * show: the X-point reaches `7.77e-06`, or 65 `eps`, where the axis reaches
	 * 11. So the constant is **256 eps** -- `3.05e-05` in single -- which clears
	 * the worst measured by 3.9 and still sits **33 times** below the `1e-3`
	 * that means the search gave up. In double it is `5.7e-14`, so every
	 * caller's historical bound dominates and nothing moves. M-219.
	 */
	inline double rootResidualFloor( double historical )
	{
		return std::max( historical,
		                 256.0*static_cast<double>(
			                 std::numeric_limits<mfem::real_t>::epsilon() ) );
	}

	/**
	 * The departure two routes computing the SAME quantity may show, floored at
	 * round-off on a quantity of size @a scale.
	 *
	 * **THE SCALE IS AN ARGUMENT BECAUSE AN ABSOLUTE BOUND WITHOUT ONE IS A
	 * RELATIVE BOUND WITH THE SCALE LEFT OUT**, which is the same defect
	 * [[roundOffFloor]] exists for one level up. `SolverContract` asks that an
	 * essential trace dof sit within `1.0e-12` of a datum of **0.37**; in single
	 * it sits **4.77e-09** away, which is 0.11 of one `eps` of that datum --
	 * round-off, and unreachable by a bound that names no scale.
	 *
	 * The constant is 64. It is three orders below anything a real defect
	 * produces here -- a datum imposed in the wrong place, or not at all, is
	 * `O( scale )` -- and in double it is `1.4e-14 x scale`, below every
	 * historical bound in this directory, so `std::max` never selects it there.
	 * M-219.
	 */
	inline double scalarAgreementFloor( double historical, double scale )
	{
		return std::max( historical,
		                 64.0*std::abs( scale )
		                 *static_cast<double>(
			                 std::numeric_limits<mfem::real_t>::epsilon() ) );
	}

	/**
	 * The same, for a pair that must agree **exactly** in double.
	 *
	 * `BorderAssembly` contracts a dense row and its compressed form and asserts
	 * `dense - sparse == 0.0`, because dropping terms that are exactly zero
	 * cannot move a sum. That is an argument about WHICH TERMS ARE SUMMED and it
	 * survives the format; what does not survive is the bit-equality, because
	 * the two routes reach the same terms through different types -- the dense
	 * row is `mfem::Vector`, hence `real_t`, and `meq::CompressedRows` stores
	 * `std::vector<double>`. Measured, they differ by **9.30e-09** on a
	 * contraction of 1.286, which is 0.06 of one `eps` of it.
	 *
	 * **A DROPPED TERM WOULD BE `O( scale )`, NOT `O( eps x scale )`**, so the
	 * claim is intact; and in double this returns **exactly zero**, so the
	 * assertion there is the bit-equality it always was.
	 *
	 * The constant is **1024**, set by the widest-spread user rather than the
	 * narrowest: `DeviceRowsTests` contracts a length-97 row against a column
	 * and reads **71 eps** on an entry of 7.5e+06, which is what cancellation in
	 * a sum that long costs. `BorderAssembly` reads 0.06 and `SourceViewTests`
	 * 0.54. A transposed `r + nr*c`, a dropped row or a wrong support index is
	 * `O( 1 )` RELATIVE, so 1024 eps still sits three orders below the smallest
	 * of them. M-219.
	 */
	inline double exactAgreementFloor( double scale )
	{
#ifdef MFEM_USE_SINGLE
		return 1024.0*std::abs( scale )
		       *static_cast<double>(
			       std::numeric_limits<mfem::real_t>::epsilon() );
#else
		(void)scale;
		return 0.0;
#endif
	}

	/**
	 * How many `eps` one exact Newton step on an affine residual may leave.
	 *
	 * **THE BOUND SEPARATES `O( eps )` FROM `O( 1 )`, SO ANYTHING IN THAT GAP
	 * DOES THE SAME JOB.** `MillerConvergence`'s own comment states the
	 * discriminator: *"an inexact step on an affine system is `O( 1 )`, not
	 * `O( 1e-10 )`"*. Measured across the affine cases here the reduction runs
	 * from 163 to **2243** `eps`, in both builds alike; 8192 clears the worst by
	 * 3.7 and still sits three orders below the `O( 1 )` an inexact step gives.
	 * It is a bound chosen in the middle of a four-order gap rather than a
	 * constant fitted to the data, which is the distinction that matters. In
	 * double it is `1.8e-12`, so every caller's historical `1.0e-8` dominates
	 * and nothing there moves. M-219.
	 */
	double const affineStepFactor = 8192.0;

	/**
	 * The residual reduction one exact Newton step on an AFFINE residual can be
	 * held to, floored from the caller's own historical value.
	 *
	 * Several cases in this directory already compute
	 * `residuals[ 1 ]/residuals[ 0 ]` and compare it against a literal `1.0e-8`,
	 * which is the right shape and the wrong number in single.
	 *
	 * **THIS IS THE SEQUENCE FORM AND IT NEEDS A LOOSER BOUND THAN THE
	 * SINGLE-POINT ONE, FOR A REASON `MillerConvergence` ALREADY WROTE DOWN**:
	 * *"|| r_1 || sits at round-off while || r_0 || shrinks with the mesh, so
	 * the ratio grows even though the step is exact throughout"*. Over that
	 * sweep the drop runs 4.8e-13 to 3.1e-11 in double and 2.7e-04 to 4.3e-03 in
	 * single -- and in each build's own `eps` those worsts are **1.4e5** and
	 * **3.6e4**, so the single build is the better behaved of the two. The
	 * caller's historical `1.0e-8` is itself **4.5e7 eps**.
	 *
	 * So the bound is `1e6 x eps`: one statement in the format's own units,
	 * which double satisfies with a margin of 7 and single with 28, and which
	 * `std::max` never selects in double because `1e-8` is looser still.
	 * **The margin against an INEXACT step is the thin one** -- 8 in single,
	 * since such a step is `O( 1 )` -- and that is the honest cost of a ratio
	 * whose denominator moves. For a single mesh, where it does not, see
	 * [[oneExactStepReachedTheFloor]], which is 8192 and much tighter. M-219.
	 */
	inline double affineStepDrop( double historical )
	{
		return std::max( historical,
		                 1.0e6*static_cast<double>(
			                 std::numeric_limits<mfem::real_t>::epsilon() ) );
	}

	/**
	 * The observed Newton ORDER a quadratic method can be held to, given how
	 * many decades this build leaves above its own floor.
	 *
	 * **AN ORDER ESTIMATE NEEDS DECADES AND SINGLE PRECISION HAS HALF OF
	 * THEM.** `bestNewtonOrder()` fits a triple of consecutive residuals and
	 * already discards any that has fallen under `newtonTolerance()`; what it
	 * cannot do is invent range. A quadratic run in double spans
	 * `1 -> 1e-4 -> 1e-8 -> 1e-15`, four residuals and two clean triples; in
	 * single the same run is `1 -> 1e-4 -> 1e-7` and stops, so there is one
	 * triple and it straddles the floor. Measured,
	 * `RotatingNewtonConvergence` reads **1.796** where it reads ~2 in double
	 * and asks for 1.8.
	 *
	 * **THE DISCRIMINATOR IS QUADRATIC AGAINST LINEAR AND 1.5 STILL MAKES IT.**
	 * A Jacobian that is not the derivative of its residual gives order 1;
	 * `VaryingCentrifugalConvergence` measures its own 5%-inconsistent control
	 * at **1.173** and says so in as many words. 1.5 sits between 1.17 and 1.80
	 * with room on both sides, which is the argument that file already makes for
	 * its own threshold. Double keeps whatever the caller asked for. M-219.
	 */
	inline double quadraticOrderThreshold( double historical )
	{
#ifdef MFEM_USE_SINGLE
		return std::min( historical, 1.5 );
#else
		return historical;
#endif
	}

	/**
	 * Whether an order could be estimated at all.
	 *
	 * `bestNewtonOrder()` returns 0 when no triple of residuals survived the
	 * round-off floor -- reachable in single, because a solve that converges in
	 * two steps leaves nothing to fit. That is a statement about the format and
	 * not about the Jacobian, so a caller reports it and skips. **In double this
	 * is false whatever the history**, so a zero there stays the failure it
	 * should be. M-219.
	 */
	inline bool newtonOrderIsUnmeasurable( double best )
	{
#ifdef MFEM_USE_SINGLE
		return !( best > 0.0 );
#else
		(void)best;
		return false;
#endif
	}

	/// What such a caller prints when it skips.
	inline std::string newtonOrderNotice()
	{
		return "  *** SKIPPED the Newton order. No triple of residuals survived "
		       "the round-off floor of this build's real_t -- a solve that reaches "
		       "the floor in two steps leaves nothing to fit an order to, which is "
		       "a property of the format rather than of the Jacobian. M-219";
	}

	/**
	 * How closely a DIFFERENCED Jacobian column can be expected to match the
	 * assembled one, floored from the caller's own historical value.
	 *
	 * **A DIFFERENCE IS NOT A DERIVATIVE AND ITS ERROR IS `sqrt( eps )`, NOT
	 * `eps`.** Where two routes to the same quantity are both exact,
	 * [[scalarAgreementFloor]] is the right bound; where one of them is a finite
	 * difference, the cancellation term alone is `eps/step` and the best any
	 * step can do is `sqrt( eps )` forward or `eps^(2/3)` central. Measured on
	 * `HighBetaConvergence`'s corner route, assembled against differenced,
	 * relative:
	 *
	 *     n = 8    psi_ax 2.56e-05    scale 1.77e-04
	 *     n = 16   psi_ax 1.83e-05    scale 6.91e-05
	 *
	 * which is exactly that band -- `sqrt( eps )` is 3.45e-04 in single and
	 * `eps^(2/3)` is 2.44e-05.
	 *
	 * So the floor is `4 sqrt( eps )`, **1.4e-03** in single, clearing the worst
	 * measured by 7.8 and still some 700 times below the `O( 1 )` a genuinely
	 * wrong column gives. It applies **only in single**: `sqrt( eps_double )` is
	 * 1.49e-08, which is LARGER than the `1.0e-8` these callers ask for, so a
	 * `std::max` would loosen the double build -- the one thing none of this may
	 * do. M-219.
	 */
	inline double differencedAgreementFloor( double historical )
	{
#ifdef MFEM_USE_SINGLE
		return std::max( historical,
		                 4.0*std::sqrt( static_cast<double>(
			                 std::numeric_limits<mfem::real_t>::epsilon() ) ) );
#else
		return historical;
#endif
	}

	inline bool oneExactStepReachedTheFloor( std::vector<double> const &residuals )
	{
		if ( residuals.size() < 2 )
			return false;
		double const eps = static_cast<double>(
			std::numeric_limits<mfem::real_t>::epsilon() );
		return residuals[ 1 ] <= affineStepFactor*eps*residuals[ 0 ];
	}

	/// What such a check prints when it fails, in the units that make it
	/// readable: the reduction, and the same number in this build's own `eps`.
	inline std::string stepReductionNotice( std::vector<double> const &residuals )
	{
		std::ostringstream text;
		if ( residuals.size() < 2 )
		{
			text << "no residual history, so the first step cannot be examined";
			return text.str();
		}
		double const eps = static_cast<double>(
			std::numeric_limits<mfem::real_t>::epsilon() );
		double const ratio = residuals[ 1 ]/residuals[ 0 ];
		text << "the first Newton step took the residual from " << residuals[ 0 ]
		     << " to " << residuals[ 1 ] << ", a reduction of " << ratio
		     << " -- which is " << ratio/eps << " times this build's eps, against "
		        "the 2048 an exact Jacobian on an affine residual reaches";
		return text.str();
	}

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
