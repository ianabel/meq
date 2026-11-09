#ifndef MEQ_SAFETY_FACTOR_SOLVE_HPP
#define MEQ_SAFETY_FACTOR_SOLVE_HPP

#include <functional>
#include <string>
#include <vector>

/**
 * @file SafetyFactorSolve.hpp
 * THE GAUSS-NEWTON THAT DRIVES `q( psi )`. MFEM-free.
 *
 * meq::SafetyFactorSystem poses `q` at every knot and one scale pin -- `m + 1`
 * rows in the `m` knot values of `g^2` -- and this solves it in least squares.
 * The rows are a whole equilibrium and a surface extraction each, so
 *
 *
 * 1. IT HAS TO BE A NEWTON, AND THAT IS A THEOREM RATHER THAN A PREFERENCE
 *
 * A relaxed Picard iteration on the inversion `g = 4 pi^2 q/( V' < R^-2 > )`
 * has derivative `1 + w( G' - 1 )` at its fixed point, above one for EVERY
 * `w > 0` when `G' > 1`. Measured, it runs away: 0.362 down to 0.019 and back
 * to 0.50, the step never shrinking.
 *
 *
 * 2. GAUSS-NEWTON BY QR, NEVER THE NORMAL EQUATIONS
 *
 * The step solves `min || J dx + r ||_2` by a column-pivoting Householder QR
 * of `J`. The `q` rows are nearly singular along the scale family, which only
 * the pin row resolves, so `J^T J` would square a conditioning that is
 * already the problem's weakest number.
 *
 *
 * 2a. AND AN EXACT CONSTRAINT, ON A CURVE THE CONSTRAINT DOES NOT CHOOSE
 *
 * A scale pin is an equality: the delivered R0B0, flux or current is the one
 * asked for. Weighted into the least squares it trades against the misfit,
 * and differently for each pin -- the answers then differ by the misfit's
 * share, 1.5e-6 of psi on QDrivenScalePin's inconsistent target (M-QR-6).
 * So the rows are met in least squares within each level set of a fixed
 * functional (the q-driven loop's: the sum of ln g^2 at the knots, the
 * scale's own coordinate), which traces one curve of answers, and the
 * constraint picks the point on it. Every pin picks from the same curve.
 *
 *
 * 3. AN ARMIJO SEARCH, AND A REFUSED TRIAL IS A REJECTED STEP
 *
 * The inner bordered Newton has no globalisation while `psi_ax` is a border
 * unknown, so the outer step length is the only control there is. A trial is
 * accepted when the merit -- `|| r ||^2`, the constraint rows weighted by
 * 1e6 -- falls by at least `c alpha` of what the linear model predicts for the
 * full step (`c` = 1e-4, `alpha` the step fraction): a sufficient decrease,
 * so that the search cannot creep (`HDG-NPC-GLOBALISATION-FROM-MEQ.md`'s
 * recorded mistake). Against the merit itself it would ask for what only a
 * zero misfit can give, and stall at the misfit.
 * A trial the residual refuses -- a solve that does not converge, a family
 * that cannot be traced, `g^2 <= 0` -- is a rejected step: there is no value
 * to compare, so nothing is compared.
 *
 *
 * 4. CONVERGED MEANS THE PART OF THE RESIDUAL A STEP CAN REMOVE IS AT THE
 *    FLOOR
 *
 * In least squares the residual does not go to zero: what remains at the
 * answer is the discrete problem's inconsistency, the misfit. What a step can
 * remove is `J dx = -Q_1 Q_1^T r`, the residual's projection on the range of
 * `J`, and the loop is converged when each row's is at or below the larger
 * of the tolerance and that row's floor -- its own round-off noise
 * (LeastSquaresOptions::floor) plus `|| J_i ||_1 eps || x ||`, what holding
 * `x` in `real_t` costs. Per row, because the rows are not equally noisy: a
 * weighted pin row carries its weight times its quantity's noise.
 *
 * **AND THE FLOOR IS MEASURED WHERE IT GOVERNS** (MEASUREMENTS.md M-363): the
 * first time the removable part comes within ten floors, and at a stall above
 * the floor before it is called one, the point is evaluated again after its
 * last unknown is moved each way -- two other warm starts -- and four
 * standard deviations of the spread become the floor if larger. A
 * measurement that is refused raises nothing.
 */

namespace meq
{

	struct LeastSquaresOptions
	{
		/// The removable residual the answer must reach, max norm, in the
		/// rows' units (relative misfits for meq::SafetyFactorSystem).
		double tolerance = 1.0e-10;

		int maxIterations = 30;

		/// Halvings of the step before the search gives up.
		int backtracks = 12;

		/// The Armijo constant.
		double sufficientDecrease = 1.0e-4;

		/// The residual's modelled round-off floor, one per row; empty is a
		/// residual exact to working precision.
		std::function<std::vector<double>()> floor;

		/// The noise measurement's move of the last unknown, relative to
		/// `max( |x_last|, 1 )`.
		double noiseProbe = 1.0e-2;

		/**
		 * CONSTRAINT ROWS: the last @ref constraints rows are met exactly, as
		 * equalities, and the rest in least squares WITHIN a level set of the
		 * functional whose gradient @ref level gives -- so the set of answers
		 * the constraints choose among is the same whichever constraint it is.
		 * Zero or one; one needs @ref level.
		 */
		std::size_t constraints = 0;
		std::function<std::vector<double>( std::vector<double> const & )> level;

		/// The epsilon of the format the iterate is held in -- mfem::real_t's,
		/// which this MFEM-free file cannot name: `|| J || epsilon || x ||` is
		/// the floor an exact residual has there.
		double epsilon = 2.220446049250313e-16;
	};

	struct LeastSquaresResult
	{
		std::vector<double> x;
		std::vector<double> residual;
		int iterations = 0;
		int evaluations = 0;
		int rejected = 0;
		bool converged = false;

		/// max |r| at @ref x: the misfit, what least squares leaves.
		double misfit = 0.0;

		/// max |J dx| at @ref x: the part a full step would remove, which
		/// @ref converged is judged on, row by row, against the larger of the
		/// tolerance and the row's floor.
		double removable = 0.0;

		/// The largest row floor at @ref x: modelled or measured, plus the
		/// iterate's own.
		double floor = 0.0;

		/// sigma_min/sigma_max of J, worst over the iteration.
		double conditioning = 1.0;

		/// The removable residual at each iterate: the convergence history.
		std::vector<double> history;

		std::string status;
	};

	/**
	 * `min || r( x ) ||_2` by Gauss-Newton with an Armijo search.
	 *
	 * @param residual `r( x )`; false refuses @a x. It may throw: a throw ends
	 *        the loop and is reported, not rethrown.
	 * @param jacobian `dr/dx` row-major at the point @a residual last
	 *        admitted -- the loop only ever asks there; false ends the loop.
	 *
	 * A rank-deficient Jacobian (`sigma_min/sigma_max` under `n eps`) is
	 * reported and not stepped through: the rows do not determine the
	 * unknowns. A zero removable residual at the start is an answer, returned
	 * at once.
	 *
	 * @throws std::invalid_argument for an empty start or a missing callable.
	 */
	LeastSquaresResult solveLeastSquares(
		std::vector<double> start,
		std::function<bool( std::vector<double> const &, std::vector<double> & )> const &residual,
		std::function<bool( std::vector<double> const &, std::vector<double> & )> const &jacobian,
		LeastSquaresOptions const &options = {} );

}

#endif
