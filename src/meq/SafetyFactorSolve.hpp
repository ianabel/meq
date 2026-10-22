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
 *
 *
 * 5. CONVERGED MEANS THE RESIDUAL REACHED max( TOLERANCE, FLOOR ), AND THE
 *    FLOOR IS THE MAP'S
 *
 * Every map evaluation is a whole equilibrium solved to a relative accuracy,
 * and a surface extraction of it, so the map is NOISY: asked twice at the same
 * coefficients from two warm starts it answers differently, by up to 1.3e-04
 * in single on `examples/q-driven.toml`. No residual below that means
 * anything, so the loop is judged against the larger of the caller's tolerance
 * and the map's own floor -- `OuterNewtonOptions::mapFloor`, which for the
 * q-driven map is meq::toroidalMapFloor(), derived there from the inner
 * solve's accuracy -- plus the floor of holding the iterate in `real_t` at all.
 * At or below that the loop is converged, whatever KINSOL's flag; above it, it
 * is not, whatever KINSOL's flag. `KIN_STEP_LT_STPTOL` above the floor is a
 * STALL, and `mfem::KINSolver::GetConverged()` counts it as success.
 *
 * Three settings follow from taking the floor seriously, and each was a
 * measured defect on `examples/q-driven.toml` in single before it was set:
 *
 *   - THE TOLERANCE IS ABSOLUTE AND IT IS KINSOL'S. MFEM's relative tolerance
 *     gave way to its default absolute one, `eps^( 1/3 )` -- 4.9e-03 in
 *     single -- so the loop reported success there, at g^2 1-2% from the
 *     closed form.
 *   - THE DIFFERENCE STEP IS SET BY THE FLOOR -- the least at which the
 *     noise costs a column 1%, `100 floor` -- not a fixed `1e-4`, which
 *     beside a noise of 1e-4 made Jacobian columns wrong at O( 1 ). In double
 *     the `1e-4` is the larger and is what is used.
 *   - A FRESH JACOBIAN EVERY ITERATION, and its product given to KINSOL. Its
 *     default keeps a Jacobian ten iterations -- here the vacuum start's, so
 *     the loop crawled linearly -- and forms `J p` for the line search by a
 *     difference of its own: an equilibrium per iteration spent on a slope
 *     made of noise.
 *
 * Together, over both MKLs and one, two and four MKL threads in single: 4
 * KINSOL iterations and 30 equilibria in every draw, to a residual of 1.6e-05
 * to 1.4e-04 against a floor of 3.8e-04, the coefficients of g^2 1.3e-05 to
 * 2.4e-04 from the closed form. The same six draws had stopped in 8 to 13
 * iterations and 33 to 87 equilibria, at 2.9e-04 once and 2.2e-03 to 4.6e-03
 * otherwise, the coefficients 6.3e-04 to 4.4e-02 out -- and three of the six
 * then failed the q the run reports. In double: 5 iterations and 36
 * equilibria to 3.5e-10, where it had stopped in 12 and 43 at 1.8e-06.
 */

namespace meq
{

	struct OuterNewtonOptions
	{
		/// The residual the answer must reach: `|| G( c ) - c ||_inf`, in the
		/// units of the coefficients, ABSOLUTE. The loop is converged when the
		/// residual at the point it returns is at or below the larger of this
		/// and mapFloor()'s value there, and not otherwise.
		double functionTolerance = 1.0e-9;

		/// KINSOL's scaled step tolerance; zero leaves KINSOL's own,
		/// `eps^( 2/3 )`. A step below it STOPS the iteration, and that stop
		/// is judged by the residual like any other.
		double stepTolerance = 0.0;
		int maxIterations = 30;

		/// Column step for the differenced Jacobian, relative to each
		/// coefficient: the LEAST step, raised where the map's floor needs it
		/// to `100 floor` -- see section 5. Central
		/// differences: two map evaluations a column, for the reason this tree
		/// gives everywhere else it differences a residual -- a forward
		/// difference leaves an `O( step )` truncation on top of the map's own
		/// noise, and here that noise is a surface extraction's.
		double difference = 1.0e-4;

		/// THE MAP'S OWN ROUND-OFF FLOOR, asked after each evaluation: the
		/// level, in the same norm and units as @ref functionTolerance, below
		/// which the residual just returned is noise. For the q-driven map
		/// that is meq::ToroidalFieldMap::floor(). Empty is a map exact to
		/// working precision -- a floor of zero.
		std::function<double()> mapFloor;

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

		/// `|| G( c ) - c ||_inf` at the coefficients returned -- the number
		/// @ref converged is judged on, against max( tolerance, @ref floor ).
		///
		/// **THE ANSWER IS WORTH THE RESIDUAL IT STOPPED AT, DIVIDED BY THE
		/// CONDITIONING BESIDE IT**, and a caller that wants to know whether
		/// its recovered `g` means anything has no other way to ask.
		double finalResidual = 0.0;

		/// The round-off floor at the coefficients returned: the map's own,
		/// from OuterNewtonOptions::mapFloor, plus that of holding the
		/// coefficients in `real_t` at all -- `|| J || eps || c ||`, which is
		/// what an EXACT map has in single. @ref converged is
		/// @ref finalResidual against max( tolerance, this ).
		double floor = 0.0;

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
