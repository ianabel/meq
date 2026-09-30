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
	 * How far past the tolerance a routine was GIVEN its achieved residual may
	 * sit, floored from the caller's own historical bound.
	 *
	 * **A TOLERANCE THAT FOLLOWS `real_t` IS NOT A CONSTANT, AND AN ASSERTION
	 * ON WHAT IT ACHIEVES CANNOT BE ONE EITHER.** `meq::newtonTolerance()`
	 * floors the tracer's `1e-12` at ten `eps`, so the corrector and the ray
	 * Newton are asked for `1.19e-06` in single and reach it: `OpenSurfaces`
	 * reads a target of `1.335e-06` and a worst residual of exactly
	 * `1.335e-06`, and `SurfaceFitConvergence` reads `8.87e-07` against the
	 * same target. Both were asserted against `1.0e-08` and `1.0e-09`, numbers
	 * chosen when the target was `1e-12` and three orders below what the
	 * routine is now allowed to stop at.
	 *
	 * **THE TARGET IS AN ARGUMENT AND NOT AN `eps` EXPRESSION**, because it is
	 * not `eps` alone: `Contour::correctorTarget` is the tolerance times the
	 * scale of `psi_h` over the mesh, and only the caller has that. Passing it
	 * is also what makes this a statement about the routine rather than about
	 * the format -- a corrector that stalls two orders above its target still
	 * fails.
	 *
	 * @a slack is 4 by default. It is not the double build's own ratio, which
	 * is 274 on `OpenSurfaces`: a corrector keeping its best iterate after four
	 * non-improving steps does not always meet its target, and in double the
	 * historical bound absorbs that. In single every measured residual sits
	 * AT the target, so 4 is slack over a routine that converged and is three
	 * orders below a trace that has wandered off the level set. M-220.
	 */
	inline double achievedToleranceFloor( double historical, double target,
	                                      double slack = 4.0 )
	{
		return std::max( historical, slack*target );
	}

	/**
	 * The residual a CONVERGED Newton may leave, floored from the caller's own
	 * historical bound.
	 *
	 * **A CONVERGENCE PRECONDITION IS A STATEMENT ABOUT THE TOLERANCE THE
	 * SOLVER WAS GIVEN, AND SEVEN SITES IN THIS DIRECTORY NAME A LITERAL
	 * INSTEAD.** Each asks `newtonResiduals().back() < 1.0e-08` while handing
	 * `solve()` a tolerance built by `meq::newtonTolerance()`, which floors at
	 * ten `eps` -- `1.19e-06` in single. So the solve reaches exactly what it
	 * was asked for and the test rejects it at a hundred times tighter. M-217's
	 * lesson one level out: *a default is only a default on the paths that do
	 * not override it*, and here the override is in the caller.
	 *
	 * The floor is `80 eps` times `max( 1, || r_0 || )`, the reference being
	 * the relative target's own multiplier read back out of the history the
	 * caller already holds -- `9.5e-06` in single on an `O( 1 )` problem. A
	 * solve that did NOT converge leaves `O( || r_0 || )`, which on the machine
	 * cases here is `1e-02` to `1e+00`, so the bound sits three to five orders
	 * inside that gap. In double it is `1.8e-14` and every historical value
	 * dominates, so the double path is unmoved. M-220.
	 */
	inline double convergedResidualFloor( double historical,
	                                      std::vector<double> const &residuals )
	{
		double const reference =
			residuals.empty() ? 1.0 : std::max( 1.0, residuals.front() );
		return std::max( historical,
		                 80.0*reference
		                 *static_cast<double>(
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
	 * The constant is **4096**, set by the widest-spread user rather than the
	 * narrowest. Within one contraction it is `DeviceRowsTests`, which reads
	 * **71 eps** on a length-97 row against an entry of 7.5e+06 -- what
	 * cancellation in a sum that long costs -- against `BorderAssembly`'s 0.06
	 * and `SourceViewTests`' 0.54. **The widest is a whole SOLVE compared
	 * through two assembly routes**: `NewtonConvergence` reads **1270 eps** on
	 * the flux at `k = 4, n = 16`, where the difference between the two paths
	 * is `eps` times the trace system's conditioning and that conditioning
	 * grows with the degree. A transposed `r + nr*c`, a dropped row, a wrong
	 * support index, a sign in `SourceIntegrator` -- every failure this bound
	 * guards is `O( 1 )` RELATIVE, so 4096 eps is `4.9e-04` and still sits
	 * three orders below the smallest of them. M-219, M-220.
	 */
	inline double exactAgreementFloor( double scale )
	{
#ifdef MFEM_USE_SINGLE
		return 4096.0*std::abs( scale )
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
	 * from 163 `eps` to **19,087**, the largest being
	 * `PedestalConvergence`'s corner-singular control at `k = 3` on the finest
	 * mesh -- where one exact step takes the residual from 6.637e-03 to
	 * 1.510e-05, which is `eps` times a trace system whose conditioning the
	 * corner has made about 1.9e+04. 65,536 clears that by 3.4 and still sits
	 * **two orders** below the `O( 1 )` an inexact step gives, an inexact
	 * Jacobian leaving the residual within a factor of about two of where it
	 * started. That is a narrower gap than this header's other constants enjoy
	 * and the number is recorded for that reason.
	 *
	 * In double it is `1.5e-11`, and it is consulted only where the iteration
	 * COUNT has already failed -- `newtonIterations <= 1 || reachedTheFloor` --
	 * which on every affine case in this tree does not happen there, so no
	 * double-build assertion reaches it at all. M-219, M-220.
	 */
	double const affineStepFactor = 65536.0;

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

	/**
	 * Whether the CEILING on an observed Newton order can be read.
	 *
	 * **THE FLOOR AND THE CEILING ON ONE NUMBER FAIL IN OPPOSITE DIRECTIONS AND
	 * ONLY ONE OF THEM SURVIVES A NARROW FORMAT.** The floor -- *this is at
	 * least superlinear* -- is a statement about the Jacobian and holds at any
	 * precision, because an inexact Jacobian contracts linearly whatever the
	 * arithmetic. The ceiling -- *this is a rate at all, rather than a short or
	 * non-monotone history read as one* -- is a statement about the FIT, and
	 * an order fitted through a last residual that is sitting on the floor is
	 * the floor's arithmetic rather than Newton's. It is unbounded above:
	 * measured **6.15** on `XPointBorder`'s bootstrap in single, where the
	 * double build reads inside its own band.
	 *
	 * So the ceiling is skipped and the floor is not, which leaves the
	 * assertion that has teeth in place. `#ifdef`'d, so the double build asks
	 * both questions exactly as it did. M-222.
	 */
	inline bool newtonOrderCeilingIsMeasurable()
	{
#ifdef MFEM_USE_SINGLE
		return false;
#else
		return true;
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
		        "the " << affineStepFactor
		     << " an exact Jacobian on an affine residual is allowed";
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

	/**
	 * The same question asked of the RATE rather than of the degree, for the
	 * studies whose observed order is not `k+1`.
	 *
	 * **THE FORMAT DOES NOT CARE WHICH POLYNOMIAL DEGREE PRODUCED AN ORDER, IT
	 * CARES HOW MANY DECADES THE SEQUENCE HAS TO CROSS.** `psi*` converges at
	 * `k+2`, so a `k = 1` post-processing study is asking for THREE orders over
	 * four dyadic meshes -- 3.6 decades -- where a `k = 1` potential study asks
	 * for two. Measured, `ExtensionConvergence`'s `psi*` at `k = 1` reads 2.354
	 * across the sequence against a design 3, and 1.480 on the finest pair,
	 * while its potential at the same degree is clean. So the predicate is the
	 * EXPECTED RATE and the cut is the same one
	 * [[rateIsMeasurable]] makes: two orders is measurable in single, three is
	 * not. `rateIsMeasurable( k )` is exactly this asked of `k + 1`. M-220.
	 */
	inline bool rateOfOrderIsMeasurable( double expected )
	{
#ifdef MFEM_USE_SINGLE
		return expected <= 2.0;
#else
		(void)expected;
		return true;
#endif
	}

	/**
	 * Whether an error is still measuring the DISCRETISATION rather than the
	 * format, on a quantity of size @a norm.
	 *
	 * **THE DEGREE-BASED SKIP IS THE RIGHT INSTRUMENT FOR A RATE ACROSS A
	 * SEQUENCE AND THE WRONG ONE FOR A SINGLE READING.** M-218 settled that no
	 * per-pair criterion can separate a good rate from a contaminated one --
	 * the tail bends before it flattens -- so `rateIsMeasurable` cuts by
	 * degree. A COMPARISON is a different question: *"psi* beats psi_h by ten
	 * times"* is about two numbers, and it is answerable exactly while both of
	 * them are above the floor. Measured on `NewtonConvergence`'s Example 5 in
	 * single, `psi*` reads 4.10e-03, 4.80e-04, 5.81e-05, 7.23e-06 at `k = 1` --
	 * the design third order, every pair -- and floors at 3.73e-07 by `k = 2`.
	 * A degree cut would throw the `k = 1` column away with the rest.
	 *
	 * The constant is **16 eps**. The floors measured here are 3.5 `eps` for
	 * `psi*`, 2.5 for a level-set position and 6.5 for a Chebyshev fit over
	 * 2886 samples, so 16 clears the worst by 2.5 and still admits everything
	 * an order above it.
	 *
	 * **IT IS `#ifdef`'d RATHER THAN SCALED, AND THAT IS THE SAFETY PROPERTY.**
	 * In double it returns true unconditionally, so a guard written with it
	 * cannot skip an assertion the double build makes today, however wrong the
	 * norm handed to it. Every helper in this header is either a `std::max`
	 * against the caller's own historical value or is `#ifdef`'d for the same
	 * reason. M-220.
	 */
	inline bool aboveTheFormatsFloor( double error, double norm,
	                                  double constant = 16.0 )
	{
#ifdef MFEM_USE_SINGLE
		return error > constant*std::abs( norm )
		       *static_cast<double>(
			       std::numeric_limits<mfem::real_t>::epsilon() );
#else
		(void)error;
		(void)norm;
		(void)constant;
		return true;
#endif
	}

	/**
	 * How well a STATIONARY point's position can be known, on a feature of
	 * size @a scale.
	 *
	 * **A MAXIMUM IS FLAT, SO LOCATING ONE COSTS HALF THE DIGITS.** Near a
	 * stationary point `psi = psi_max - c dz^2`, so a flux known to `eps` pins
	 * the position only to `sqrt( eps/c )` -- `3.45e-04` in single against
	 * `1.49e-08` in double, four and a half orders apart where the format is
	 * one. It is not a tolerance anybody chose and no amount of refinement
	 * moves it: it is the square root in the relation.
	 *
	 * Measured, and the shape is exactly that. `LimiterCurve`'s located
	 * contact stands 1.177e-03, 2.408e-04, 1.536e-04 from the midplane over
	 * `n = 12, 24, 48` -- converging at 2.289 and then FLOORING at 1.5e-04,
	 * which is `sqrt( eps )` on a half-width of 0.5. `AxisAgreement`'s located
	 * axis reads `R = 0.99999583` against an exact 1 and `z = -3.16e-06`
	 * against 0. Both were asserted against `1.0e-08`.
	 *
	 * **THE VALUE AT THE POINT IS UNAFFECTED, WHICH IS THE OTHER HALF AND IS
	 * WHY THIS IS NOT A LOSS.** The same flatness that hides the position
	 * makes `psi` there accurate to `eps`: an `O( sqrt( eps ) )` error in
	 * WHERE costs `O( eps )` in the VALUE. That is the envelope argument
	 * `LimiterCurve` already rests on, arriving from the other side.
	 *
	 * `#ifdef`'d true in double, so a guard written with it cannot skip an
	 * assertion that build makes today. M-220.
	 */
	inline bool stationaryPositionIsMeasurable( double offset, double scale )
	{
#ifdef MFEM_USE_SINGLE
		return offset > 4.0*std::abs( scale )
		       *std::sqrt( static_cast<double>(
			       std::numeric_limits<mfem::real_t>::epsilon() ) );
#else
		(void)offset;
		(void)scale;
		return true;
#endif
	}

	/// The same floor as a NUMBER, to `std::max` against a caller's own
	/// historical bound on how far a located stationary point may sit from a
	/// known one: `1.38e-03` on a unit feature in single and **exactly zero**
	/// in double, so the caller's own bound is what runs there, to the bit.
	/// See [[stationaryPositionIsMeasurable]] for why the exponent is a half,
	/// and [[exactAgreementFloor]] for why it is `#ifdef`'d rather than
	/// evaluated in both formats -- `sqrt( eps )` in double is 1.49e-08, which
	/// is LARGER than several of the bounds this would be maxed against, so a
	/// scaled version would loosen the double build rather than leave it
	/// alone.
	inline double stationaryPositionFloor( double scale )
	{
#ifdef MFEM_USE_SINGLE
		return 4.0*std::abs( scale )
		       *std::sqrt( static_cast<double>(
			       std::numeric_limits<mfem::real_t>::epsilon() ) );
#else
		(void)scale;
		return 0.0;
#endif
	}

	/**
	 * Whether @a value stands far enough above a KNOWN contaminant for a rate
	 * through it to be a rate of @a value rather than of the two together.
	 *
	 * **THIS IS NOT A ROUND-OFF FLOOR, IT IS A SECOND ERROR THAT DOES NOT
	 * MOVE.** `WarmStartConvergence` measures a structured grid's bilinear
	 * interpolation against the coarse solve it interpolates, and that solve's
	 * own discretisation error is a constant the grid sequence cannot go below.
	 * The file already says so: *"by 257 nodes the grid error is within a
	 * factor of two of the coarse solve's own, which the grid cannot represent
	 * better than"*. In single the solve's error carries the format in it as
	 * well, the gap narrows further, and the last rate reads 1.689 against a
	 * bound of 1.7 set from the double build's 1.735.
	 *
	 * `#ifdef`'d true in double for the same reason as
	 * [[aboveTheFormatsFloor]]: the published rates are that build's, and a
	 * guard must not be able to drop one of them. M-220.
	 */
	inline bool clearOfContaminant( double value, double contaminant,
	                                double factor = 2.0 )
	{
#ifdef MFEM_USE_SINGLE
		return value > factor*std::abs( contaminant );
#else
		(void)value;
		(void)contaminant;
		(void)factor;
		return true;
#endif
	}

	/**
	 * Whether a SECOND difference of measured errors is readable in this
	 * format -- a drift, a ratio of one measured ratio to the next.
	 *
	 * **A RATIO OF RATIOS HAS THE ROUND-OFF IN IT TWICE AND THE SIGNAL ONCE.**
	 * `EstimatorConvergence`'s effectivity index is `eta/L2( q )` and it is
	 * asserted to settle to within 1.5 per cent between consecutive meshes,
	 * which is the sharpest statement in this directory -- it is what separates
	 * an estimator with the right power of `h` in it from one that merely
	 * shares a rate. Reading it needs both columns clean to well under a per
	 * cent. Measured in single at `k = 1`, the cleanest case in that file, the
	 * index runs 10.71, 10.60, 10.53, **9.66** against a double build's 10.72,
	 * 10.61, 10.55, 10.52: the finest `L2( q )` is 9 per cent contaminated and
	 * the drift reads 0.918 against a band of 0.95 to 1.05.
	 *
	 * **THE BOUNDS ON THE INDEX ITSELF ARE A DIFFERENT QUESTION AND STILL RUN**
	 * -- reliable above and efficient below, 1 to 60 -- and they hold at every
	 * mesh and every degree in single. What is skipped is the drift alone.
	 *
	 * `#ifdef`'d, so the double build asserts everything it always did. M-220.
	 */
	inline bool driftIsMeasurable()
	{
#ifdef MFEM_USE_SINGLE
		return false;
#else
		return true;
#endif
	}

	/**
	 * Whether an ITERATION COUNT may be held to a value measured in another
	 * format.
	 *
	 * **A COUNT IS AN ARTEFACT OF THE ARITHMETIC AS MUCH AS OF THE METHOD**, and
	 * it moves in the direction nobody expects: a narrower format makes the
	 * element-local non-linear solves noisier, so the OUTER Newton takes more
	 * steps to reach a target that is itself looser. Measured on
	 * `PedestalConvergence`'s `k = 2, n = 16` -- the case
	 * `pedestalConvergenceIsAResolutionThreshold` exists to document as a knife
	 * edge -- single takes **31** iterations where double takes 8 to 11.
	 *
	 * **WHAT THE CASE ACTUALLY CLAIMS IS UNAFFECTED AND IS STILL ASSERTED**:
	 * that `p`-refinement cures the coarse case, which is `converged` being
	 * true. M-219's lesson is to assert the claim rather than the count, and
	 * here the claim already has its own assertion beside the count.
	 *
	 * `#ifdef`'d, so the double build holds every count it always did. M-220.
	 */
	inline bool iterationCountIsComparable()
	{
#ifdef MFEM_USE_SINGLE
		return false;
#else
		return true;
#endif
	}

	/**
	 * The floor on a RESIDUAL ESTIMATOR, on a potential of size @a norm over a
	 * mesh of spacing @a h, to `std::max` against a caller's own ceiling.
	 *
	 * **AN ESTIMATOR DIFFERENTIATES, SO ITS FLOOR RISES AS THE MESH FALLS --
	 * THE OPPOSITE DIRECTION TO THE ERROR IT ESTIMATES.** `eta_1` is
	 * `h_K || f + Delta* psi_h ||_K` and that residual is a cancellation of
	 * second derivatives, so its round-off part carries a negative power of
	 * `h` where every term it is compared against carries a positive one. The
	 * two therefore CROSS, and past the crossing refining makes `eta` larger.
	 *
	 * Measured on `EstimatorConvergence` at `k = 3` in single, over
	 * `h = 0.2, 0.1, 0.05, 0.025`: `eta` reads 8.93e-05, 9.73e-05, 1.97e-04,
	 * 3.93e-04 -- rising, and doubling with each halving of `h`, so the floor
	 * goes as `1/h` on this fixture. The double build's own column falls to
	 * 1.93e-08 on the same meshes and never approaches it. The constant is
	 * **512**, which puts the floor at 1.22e-03 on the finest mesh here --
	 * above the 3.93e-04 measured by 3.1, and 2.27e-12 in double, below every
	 * ceiling in that file by eight orders, so nothing there moves.
	 *
	 * **AND IT IS WHY AN ADAPTIVE LOOP IN A NARROW FORMAT STOPS BEING A
	 * DESCENT.** M-220's `AdaptiveRefinement` entry is the same crossing seen
	 * from the loop's side: past it the marker is steering on round-off.
	 */
	inline double residualEstimatorFloor( double norm, double h )
	{
		if ( !( h > 0.0 ) )
			return 0.0;
		return 512.0*std::abs( norm )
		       *static_cast<double>(
			       std::numeric_limits<mfem::real_t>::epsilon() )/h;
	}

	/**
	 * The RELATIVE round-off floor on a flux-surface average, on a mesh of
	 * spacing @a h.
	 *
	 * **THE SAME LAW AS [[residualEstimatorFloor]] ONE DERIVATIVE FURTHER ON,
	 * AND IT IS WHY REFINING MAKES AN AVERAGE WORSE IN A NARROW FORMAT.** `V'`
	 * and `< R^-2 >` are quadratures over a traced level set with weight
	 * `1/( R | grad psi | )`, and the arc-length metric that parametrises them
	 * comes off `q` -- a first derivative of `psi` -- whose own variation along
	 * the contour is a second. So the round-off part carries `h^-2` where every
	 * term it is compared against carries `h^(k+1)`, the two cross, and past
	 * the crossing the average is worse on the finer mesh.
	 *
	 * Measured on `SurfaceAverageConvergence` in single at `k = 2`, over
	 * `n = 12, 24, 48, 96`: the rms relative error in `V'` reads 6.64e-06,
	 * 8.80e-06, 1.54e-04, 3.31e-04 -- RISING, and by a factor of 50 over an
	 * eight-fold refinement, which is `h^-1.9`. The same column at `k = 3`
	 * behaves the same way. Against `eps/h^2` those four readings are 0.39,
	 * 0.13, 0.36 and 0.30 of it, so the law fits to within a factor of three
	 * across the sweep and the constant here is **4**, an order above the worst
	 * of them.
	 *
	 * In double the same expression at `n = 96` is `8.2e-12`, eight orders
	 * below the `1e-04` ceilings in that file, so nothing there moves. M-220.
	 */
	inline double averagedMetricFloor( double h )
	{
		if ( !( h > 0.0 ) )
			return 0.0;
		return 4.0*static_cast<double>(
			       std::numeric_limits<mfem::real_t>::epsilon() )/( h*h );
	}

	/**
	 * How far a TRACED point may sit from the level set it was corrected onto,
	 * on a mesh of spacing @a h and a field of size @a norm.
	 *
	 * **THE THIRD MEMBER OF THE SAME FAMILY**, between
	 * [[residualEstimatorFloor]]'s `1/h` and [[averagedMetricFloor]]'s `1/h^2`.
	 * The corrector is Newton along a line on `psi_h`, so its step is
	 * `( psi_h - c )/| grad psi_h |` and the round-off in the numerator is
	 * divided by a gradient the mesh resolves to `eps/h`; the distance it
	 * leaves therefore carries `1/h` as well.
	 *
	 * Measured on `FluxSurfaceConvergence`'s field-error study in single, the
	 * traced contour's distance from the exact one over `n = 16, 32, 64`:
	 *
	 *     k = 1   9.120e-03  1.076e-03  5.042e-04     rates 3.083, 1.093
	 *     k = 2   6.365e-05  8.641e-05  1.991e-04     rates -0.441, -1.205
	 *
	 * -- the `k = 2` column RISING with refinement, which is the crossing. Each
	 * of the floored readings is 21 to 53 times `eps/h`, so the constant is
	 * **64**: it admits `k = 1`'s first pair, whose 3.083 is a real measurement,
	 * and turns down every reading past the crossing. In double it is
	 * `1.1e-12` at the finest `h` here, ten orders below the coarsest distance
	 * in that table. M-220.
	 */
	inline double tracedPositionFloor( double norm, double h )
	{
		if ( !( h > 0.0 ) )
			return 0.0;
		return 64.0*std::abs( norm )
		       *static_cast<double>(
			       std::numeric_limits<mfem::real_t>::epsilon() )/h;
	}

	/**
	 * Whether an adaptive loop's estimator can still be asked to DESCEND.
	 *
	 * **PAST THE CROSSING, REFINING RAISES THE ERROR, AND THAT IS THE
	 * ESTIMATOR REPORTING WHAT THE ARITHMETIC GIVES IT RATHER THAN THE
	 * ESTIMATOR FAILING.** [[residualEstimatorFloor]] is the quantity: `eta`
	 * differentiates, so its floor goes as `eps/h` and RISES as the mesh falls,
	 * while the truncation falls as `h^(k+1)`. The two cross, and an adaptive
	 * loop past that point is steering on round-off -- neither `eta` nor the
	 * true error is monotone, and the marker spends its budget wherever the
	 * conditioning is worst.
	 *
	 * Measured on `DriverAcceptance`'s coupled half-disc, four cycles from
	 * `n_R = 32` on a box of 1.7, in single:
	 *
	 *     4.46e-03   1.08e-02   2.68e-02   7.29e-02
	 *
	 * -- a factor of about 2.4 a cycle, against the 2 that halving `h` adds to
	 * `eps/h`. The same loop in double comes down over the same cycles.
	 *
	 * **`#ifdef`'d RATHER THAN COMPARED AGAINST A FLOOR, AND THAT IS THE
	 * HONEST FORM HERE.** A caller cannot evaluate `residualEstimatorFloor()`
	 * at the right `h`: the loop refines ADAPTIVELY, so the mesh that sets the
	 * floor at cycle 3 is a handful of elements several bisections below
	 * anything the configuration names, and the coarsest `h` gives a floor four
	 * times too small. Pretending to a threshold that is computed from the
	 * wrong `h` is worse than a skip that says what it is.
	 *
	 * **WHAT IS SKIPPED IS THE DESCENT ALONE.** That the loop refined, that it
	 * refined adaptively rather than uniformly, and that whatever the case
	 * couples to survived it are all separate assertions and all still run --
	 * and the defect this particular descent was written against, a spurious
	 * `eta_5` term on the coupled path, is one this format cannot tell from its
	 * own floor. M-219's third distinction: a controlled experiment whose
	 * signal is below the floor is a double-build acceptance, not a tolerance
	 * to widen. M-220, M-222.
	 */
	inline bool estimatorDescentIsMeasurable()
	{
#ifdef MFEM_USE_SINGLE
		return false;
#else
		return true;
#endif
	}

	/**
	 * How far an answer may sit from a known root, given the residual the
	 * solver actually STOPPED at.
	 *
	 * **THE SIBLING OF [[achievedToleranceFloor]] WITH THE REQUESTED NUMBER
	 * REPLACED BY THE ACHIEVED ONE, AND THE DIFFERENCE MATTERS WHENEVER A
	 * SOLVER CAN STOP FOR A SECOND REASON.** That one says a solver cannot beat
	 * the tolerance it was given. This one is needed where it did not MEET that
	 * tolerance and returned anyway: KINSOL stops on a step that has become too
	 * small exactly as readily as on a function norm that has become small
	 * enough, reports `converged` for both, and the two are told apart only by
	 * the norm it came to rest at.
	 *
	 * Measured on `SafetyFactorSolver`'s closed-form case, whose map is one
	 * whole equilibrium solve and a flux-surface extraction. Asked for `1e-09`
	 * in both builds:
	 *
	 *     double   final norm 3.68e-06   recovered g within 1e-04
	 *     single   final norm 3.71e-03   recovered g out by 7.0e-03
	 *
	 * -- and `7.0e-03` is `1.9` times that residual, which is the outer
	 * Jacobian's own conditioning working on it. **The single figure is not a
	 * tolerance to be relaxed and is not KINSOL's**: a safety factor is a
	 * flux-surface average, whose floor rises as `eps/h^2` (M-220), so at
	 * `h = 1/24` the map's own evaluation noise is `6.9e-05` before the fit
	 * amplifies it. Nothing downstream of that map can be asked for more.
	 *
	 * The constant is **8**, which clears the measured 1.9 by four. It
	 * `std::max`es against the caller's own value as everything here does, and
	 * in double that is what it picks: `8 x 3.68e-06` is `2.9e-05` against a
	 * bound of `1e-04`, so the double assertion is unchanged with a margin of
	 * 3.4 -- stated rather than assumed, because unlike an `eps` expression
	 * this one is inert in double by MEASUREMENT and not by construction. If a
	 * double run ever comes to rest above `1.3e-05` this bound starts to give,
	 * and that would itself be the finding. M-222.
	 */
	inline double achievedResidualFloor( double historical, double finalResidual )
	{
		return std::max( historical, 8.0*std::abs( finalResidual ) );
	}

	/**
	 * How far a bordered CONSTRAINT may sit from the value it prescribes, given
	 * the residual the bordered Newton came to rest at.
	 *
	 * **A BORDER ROW IS ONE ENTRY OF A RESIDUAL AND THE CONSTRAINT IT CARRIES
	 * IS THAT ENTRY WORKED ON BY THE BORDERED JACOBIAN'S CONDITIONING**, so
	 * neither the tolerance ([[achievedToleranceFloor]]) nor the residual
	 * ([[achievedResidualFloor]]) is the whole of it. Measured on
	 * `FreeBoundaryCoupling`'s two-border machine in single, seven rows of one
	 * table, prescribing `mu0 I_p = 0.12`:
	 *
	 *     R_lim  coils   delivered - 0.12    ||r||      ratio
	 *     1.08   yes         1.271e-05     3.258e-07      39
	 *     1.10   yes         9.434e-05     6.330e-07     149
	 *     1.10   NO          1.097e-05     1.635e-07      67
	 *     1.12   yes         8.672e-05     1.778e-07     488
	 *     1.12   NO          1.021e-06     8.338e-08      12
	 *     1.15   yes         1.690e-06     1.356e-08     125
	 *     1.18   yes         2.405e-06     1.118e-08     215
	 *
	 * -- a ratio between 12 and 488, which is the conditioning and is a
	 * property of the geometry rather than of the format. The constant is
	 * **4096**, clearing the worst row by 8.4.
	 *
	 * **`#ifdef`'d RATHER THAN `std::max`'d, AND THAT IS NOT STYLE.** The
	 * double build reaches `||r|| <= 7.7e-11` on these rows -- `1e-09` relative
	 * against a reference of 0.077 -- and `4096` of that is `3.2e-07` against a
	 * historical bound of `1.2e-06`. It would not loosen today, and it would be
	 * within a factor of four of doing so, which is not a margin to leave a
	 * bound standing on. A floor is owed only to a format that has one. M-222.
	 */
	inline double borderedConstraintFloor( double historical,
	                                       double finalResidual )
	{
#ifdef MFEM_USE_SINGLE
		return std::max( historical, 4096.0*std::abs( finalResidual ) );
#else
		( void ) finalResidual;
		return historical;
#endif
	}

	/**
	 * How far two positions ROOTED on one level set may sit apart, on a field
	 * whose coordinates are of size @a scale.
	 *
	 * **THE SIBLING OF [[tracedPositionFloor]] WITH THE MESH TAKEN OUT, AND
	 * THE DIFFERENCE IS WHICH QUESTION WAS ASKED OF THE FIELD.** That one is a
	 * position obtained by DIFFERENCING -- a corrector step is
	 * `( psi_h - c )/| grad psi_h |`, whose denominator the mesh resolves to
	 * `eps/h`, so the distance carries `1/h`. This one is the gap between two
	 * points that are each already ON the level set: the implicit function
	 * theorem moves a root by `dpsi/| grad psi |`, and BOTH factors carry the
	 * field's own size, so what survives is `eps` times the length over which
	 * the field varies and `h` does not appear at all.
	 *
	 * Measured on `thePointErrorDoesNotAccumulateOverManyCircuits` in single,
	 * where the closure gap's normal component is exactly such a pair: **1.00e-07,
	 * 7.39e-07 and 6.42e-07** over 1, 5 and 10 circuits of a contour whose
	 * coordinates are `R = 1.05`, which is 0.8 to 6.2 `eps` of that scale --
	 * and FLAT in the path length, which is the property the case exists to
	 * assert and which survives the format intact. The bound it is added to is
	 * the curve's own bending over the tangential gap, `kappa g_t^2`, which is
	 * 7.5e-11 there: a geometric bound with no arithmetic in it, three orders
	 * under a floor the format cannot go below.
	 *
	 * **AND THE SAME NUMBER FALLS OUT OF A SECOND TABLE**, which is why the
	 * constant is 64 rather than the 6 the closure gap alone would set.
	 * `theOutermostSurfaceIsTracedEntirelyInTheBand` traces Gamma itself and
	 * measures its distance from the exact Gamma; the transfer lift's column
	 * settles at **4.73e-06 at k = 2 and 5.88e-06 at k = 3** -- the same number
	 * at two degrees, which is what says it is the format and not the method --
	 * and that is 40 to 50 `eps`. 64 admits every reading either table still
	 * calls a measurement and turns down both floored ones.
	 *
	 * **EXACTLY ZERO IN DOUBLE**, as [[stationaryPositionFloor]] is, rather
	 * than a small number: one caller adds it to an additive term of `1e-14`,
	 * and `64 eps` of a double is `1.4e-14`, so it would loosen by half an
	 * order the one bound the measurement is real under. A floor is owed only
	 * to a format that has one. M-222.
	 */
	inline double rootedPositionFloor( double scale )
	{
#ifdef MFEM_USE_SINGLE
		return 64.0*std::abs( scale )
		       *static_cast<double>(
			       std::numeric_limits<mfem::real_t>::epsilon() );
#else
		( void ) scale;
		return 0.0;
#endif
	}

	/// The message a caller prints when [[rateIsMeasurable]] refuses, so that the
	/// weakening is visible in the output rather than inferred from a gap.
	/// The same notice for a quantity whose order is not `k+1` -- `psi*` at
	/// `k+2`, an estimator component, a differenced metric.
	inline std::string rateNotMeasurableNotice( int order, char const *quantity,
	                                            double expected )
	{
		std::ostringstream text;
		text << "  *** k = " << order << ": SKIPPED the rate assertion on "
		     << quantity << ", whose design order here is " << expected
		     << ". This build is MFEM_PRECISION = single, in which a sequence "
		        "has to cross that many decades before it ends and the error "
		        "reaches the format's floor first. The absolute ceilings and "
		        "every comparison that is not a rate are still asserted. M-220";
		return text.str();
	}

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
