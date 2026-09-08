#ifndef MEQ_SAFETY_FACTOR_SOLVE_HPP
#define MEQ_SAFETY_FACTOR_SOLVE_HPP

#include <functional>
#include <string>
#include <vector>

/**
 * @file SafetyFactorSolve.hpp
 * THE OUTER NEWTON THAT DRIVES `q( psi )`, ON KINSOL.
 *
 * `meq::SafetyFactor` is the inversion -- `g = 4 pi^2 q/( V' < R^-2 > )`, one
 * division per surface, MFEM-free and gated by CI. This is the loop around it:
 * `V'` and `< R^-2 >` are functionals of the solution, so the fixed point of
 *
 *     c  ->  solve, extract, invert, fit  ->  c
 *
 * is the equilibrium whose own `q` is the one asked for, with `c` the
 * coefficients `g^2` is fitted in.
 *
 *
 * 1. IT HAS TO BE A NEWTON, AND THAT IS A THEOREM RATHER THAN A PREFERENCE
 *
 * A relaxed Picard iteration `c <- c + w( G( c ) - c )` has derivative
 * `1 + w( G' - 1 )` at the fixed point. For `G' > 1` that is above one for
 * EVERY `w > 0`: under-relaxation stabilises a map that oscillates and can do
 * nothing at all for one that runs away. Measured, this map runs away -- the
 * damped loop walks the error from 0.362 down to 0.019 and then back up to
 * 0.50, with the step never shrinking, including where the error passes through
 * zero. It does not stall at the fixed point, it crosses it.
 *
 *
 * 2. AND IT IS AFFORDABLE BECAUSE THE PROFILE IS FITTED
 *
 * `c` is a handful of coefficients rather than `nFieldDOF`, so a differenced
 * Jacobian costs a few map evaluations per step -- each one a solve and an
 * extraction -- against `INVERSION-PLAN.md` section 11.1's **5.7 hours** for
 * `dGeometry_dpsi`, which differences against every field degree of freedom.
 * The fit was put there for conditioning and it pays for the Newton as well.
 *
 *
 * 3. KINSOL RATHER THAN A HAND-ROLLED ITERATION, AND THE REASON IS THE LINE
 *    SEARCH
 *
 * An undamped Newton on this map takes a first step the INNER bordered Newton
 * cannot solve at -- and that inner solve has no globalisation of its own,
 * because `GradShafranovSolver` refuses every `Globalisation` but `None` while
 * `psi_ax` is a border unknown. So the outer step length is the only control
 * there is, and it needs a real line search.
 *
 * **WRITING ONE HERE WOULD BE REPEATING A MISTAKE THIS TREE HAS ALREADY PAID
 * FOR.** `HDG-NPC-GLOBALISATION-FROM-MEQ.md` records a hand-rolled backtracking
 * search accepting on a MONOTONE test with no sufficient-decrease constant: for
 * Newton on an `l2` merit the direction is always a descent direction, so any
 * small enough step "succeeds", and the iteration creeps by about 1% a step
 * instead of failing honestly. KINSOL's `KIN_LINESEARCH` applies an Armijo
 * condition and does not have that failure mode. It is also already linked,
 * already used by `GradShafranovSolver`, and somebody else's to maintain --
 * which is this tree's standing preference for a well-known algorithm.
 *
 *
 * 4. AND THE PRICE OF THAT LINE SEARCH IS THAT A SINGULAR JACOBIAN MUST NEVER
 *    REACH IT
 *
 * `KIN_LINESEARCH` interpolates its step length on a quotient whose numerator
 * and denominator both carry the directional derivative `< F, J p >`. A
 * singular `J` makes that zero whatever the step is, so the quotient is `0/0`,
 * the iterate goes to NaN, and **KINSOL then never returns** -- every one of
 * its convergence and failure tests is a comparison against NaN and every
 * comparison against NaN is false. Measured on a two-variable map with no root
 * at all: still calling the map, at NaN, after **two million evaluations**,
 * each of which on the real loop is an equilibrium solve.
 *
 * **NEITHER OBVIOUS REPAIR REACHES IT**, and both were tried. Bounding the
 * linear solver's iteration count does not, because `mfem::GMRESSolver` divides
 * by a zero pivot in its own back-substitution on the FIRST iteration and comes
 * back infinite. Replacing it with a rank-revealing dense solve does not
 * either, even though that correctly returns a ZERO step -- because it is
 * `< F, J p >` and not the step that has gone to zero. **The line search cannot
 * be rescued from outside it**, so the degeneracy is detected before KINSOL is
 * entered at all.
 *
 * The dense truncated-SVD solve is kept regardless: it is the right solver for
 * a small dense system, and it bounds a Jacobian that goes singular partway
 * through, where the pre-flight check cannot see it.
 */

namespace meq
{

	struct OuterNewtonOptions
	{
		/// Passed to KINSOL as its function-norm and step tolerances.
		double functionTolerance = 1.0e-9;
		double stepTolerance = 0.0;      ///< zero leaves KINSOL's own default
		int maxIterations = 30;

		/// Column step for the differenced Jacobian, relative to each
		/// coefficient with this as the floor. Central differences: two map
		/// evaluations a column, for the reason this tree gives everywhere else
		/// it differences a residual -- a forward difference leaves an
		/// `O( step )` truncation on top of the map's own noise, and here that
		/// noise is a surface extraction's.
		double difference = 1.0e-4;

		/// KIN_LINESEARCH by default. KIN_NONE is kept as the control: it is
		/// what says the line search is doing the work rather than the Newton
		/// direction being good enough on its own.
		bool lineSearch = true;
	};

	struct OuterNewtonResult
	{
		std::vector<double> coefficients;
		int iterations = 0;
		int mapEvaluations = 0;
		bool converged = false;

		/// `sigma_min/sigma_max` of the outer Jacobian, worst over the
		/// iteration -- 1 for an orthogonal one and 0 for a singular one. It is
		/// what says whether the coefficients determine their own fixed point,
		/// and it is reported rather than asserted on because the answer is a
		/// property of the degree asked for and of how many surfaces the family
		/// has, both of which are the caller's.
		double jacobianConditioning = 1.0;

		/// KINSOL's own account of how it stopped, for a caller that has to
		/// report rather than only branch.
		std::string status;
	};

	/**
	 * Solve `G( c ) = c` for the coefficients of `g^2`.
	 *
	 * @param map one whole Picard step: given the coefficients, solve the
	 *        equilibrium, extract its surfaces, invert the target safety factor
	 *        against them, and return the coefficients that implies. Every
	 *        solve lives inside this callable.
	 *
	 * @throws std::invalid_argument if @a start is empty, if the map returns a
	 *         different number of coefficients than it was given, or if the
	 *         options do not describe an iteration.
	 * @throws whatever @a map throws -- in particular a failed inner solve,
	 *         which is a state the caller knows how to describe and this does
	 *         not. **A map that cannot be evaluated at a trial point should
	 *         return a large finite residual rather than throw**: KINSOL is C,
	 *         an exception unwinds through its frames, and the line search
	 *         needs a value there to reject the step with.
	 *
	 * A target the coefficients do not determine -- a degree higher than the
	 * surface family supports -- leaves the outer Jacobian rank deficient.
	 * That is REPORTED and not thrown: `converged` is false,
	 * `jacobianConditioning` is the measurement, and `status` says which
	 * directions were undetermined and what to do about it. Reporting is the
	 * right answer rather than a soft one, since the caller knows what its
	 * degree means and this does not.
	 *
	 * **A ZERO JACOBIAN IS ALSO WHAT AN ALREADY-SOLVED PROBLEM HAS**, and the
	 * two are opposite answers: `G = identity` makes every point a fixed point,
	 * so returning at once is right there where refusing would be wrong. The
	 * residual separates them and both are reported as such.
	 */
	OuterNewtonResult solveForToroidalField(
		std::vector<double> start,
		std::function<std::vector<double>( std::vector<double> const & )> const
			&map,
		OuterNewtonOptions const &options = {} );

}

#endif
