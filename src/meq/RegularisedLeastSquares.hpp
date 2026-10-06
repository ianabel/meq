#ifndef MEQ_REGULARISED_LEAST_SQUARES_HPP
#define MEQ_REGULARISED_LEAST_SQUARES_HPP

#include <vector>

/**
 * @file RegularisedLeastSquares.hpp
 * THE INVERSE SOLVE'S DENSE CORNER, AS A FREE FUNCTION WITH NO SOLVER ATTACHED.
 *
 * `INVERSE-AND-RECONSTRUCTION-PLAN.md` section 4.3 and stage IV-2. An inverse
 * solve makes some coil currents unknowns of the bordered Newton and adds
 * shape targets as rows, so the dense corner the bordered step inverts stops
 * being square: `nBorders + N_targets` rows against `nBorders + N_free`
 * columns, 21 x 33 on DIII-D. Everything else in the step -- the Jacobian, its
 * factorisation, the backsolves -- is unchanged, so the whole of the new
 * numerical content is this one small solve:
 *
 *     minimise  || M x - f ||^2 + || diag( gamma ) x ||^2
 *
 * which is algebraically `freegs4e`'s `x = ( A^T A + gamma^2 I )^-1 A^T b`
 * applied to the Newton step in the corner.
 *
 * **ONE CODE PATH FOR BOTH REGIMES, AND THE MACHINES NEED BOTH.** With targets
 * the coils can reach, the residual goes to zero and this is a Newton step on an
 * underdetermined system -- `gamma` selects WHICH root, not whether one exists.
 * With targets they cannot reach, it is the best-fit Gauss-Newton step. Machine A
 * has four coils against five target rows and is overdetermined; DIII-D has
 * eighteen against six and is under. A square solve fails the first and an
 * unregularised one is arbitrary on the second.
 *
 * **gamma IS PER COLUMN AND ZERO MEANS UNREGULARISED**, so the border unknowns
 * that are not currents -- `psi_ax`, `psi_bnd`, the exterior modes -- carry a
 * zero and only the current block is damped. Section 4.4: the natural weight on
 * coil c is `1/M_c`, an inverse mutual inductance, and `gamma -> infinity`
 * damps the current step to nothing, which is the forward problem.
 *
 * **SOLVED AS THE AUGMENTED SYSTEM `[ M ; diag( gamma ) ] x ~ [ f ; 0 ]` BY A
 * SINGULAR VALUE DECOMPOSITION, NEVER BY THE NORMAL EQUATIONS**, which square the
 * condition number of a matrix whose columns are coil responses of very
 * different sizes. Among minimisers the minimum-norm one is returned, so a
 * rank-deficient corner -- targets that leave some currents undetermined -- has
 * a defined step rather than an arbitrary one.
 *
 * **THE RANK IS REPORTED BECAUSE IT IS THE HONEST ANSWER TO "HOW MANY CURRENTS
 * DID THE TARGETS DETERMINE"**: section 1.4, Lackner's warning that a shape can
 * be matched to plotting accuracy by a current set that is nothing like the
 * machine's. `rank` and `nullity` are of `M` alone, before the regularisation
 * makes the system full rank.
 *
 * **EIGEN AND NOT `mfem::DenseMatrixSVD`, AND IT WAS MEASURED.** At the corner's
 * size both are accurate to 3e-11 against a closed form at a condition number of
 * 1e6, and both cost under a tenth of a millisecond against a Newton step's
 * sparse factorisation. What decided it: Eigen's `JacobiSVD` solves for the
 * minimum-norm minimiser itself where MFEM's returns factors and leaves the
 * pseudo-inverse to be written and maintained here; this file stays MFEM-free
 * and so CI-tested, as `Profiles` and `Config` are; and `DenseMatrixSVD`
 * OVERFLOWS ITS RIGHT SINGULAR VECTORS for a wide matrix with `jobvt = 'S'` --
 * it allocates `min( m, n ) x n` and hands LAPACK a leading dimension of `n` --
 * which is the shape this corner has in the normal, underdetermined case.
 * Reported upstream as `DENSEMATRIXSVD-LDVT-FROM-MEQ.md`. Eigen stays inside the
 * `.cpp`, as `SurfaceFit`'s does.
 */
namespace meq
{
	/// What solveRegularisedLeastSquares() returns.
	struct RegularisedStep
	{
		/// The minimum-norm minimiser, one entry per column of `M`.
		std::vector<double> step;

		/// The numerical rank of `M` alone, and `columns - rank`: how many
		/// directions of the unknowns the rows leave undetermined, before the
		/// regularisation fixes them.
		int rank = 0;
		int nullity = 0;

		/// The singular values of `M`, largest first. Their ratio is the
		/// conditioning a report of the rank should be read beside.
		std::vector<double> singularValues;

		/// `|| M step - f ||`: zero to round-off where the targets are
		/// reachable, and the size of the best fit where they are not.
		double residualNorm = 0.0;

		/// `|| diag( gamma ) step ||`, the regularisation's share.
		double penaltyNorm = 0.0;
	};

	/// The default relative threshold below which a singular value counts as
	/// zero: `max( rows, columns ) * eps` of the largest, which is LAPACK's and
	/// Eigen's own default for the same decision.
	double defaultRankThreshold( int rows, int columns );

	/**
	 * @brief The minimum-norm minimiser of
	 *        `|| M x - f ||^2 + || diag( gamma ) x ||^2`.
	 *
	 * @param matrix   `M`, `rows x columns`, ROW-MAJOR: entry ( i, j ) is
	 *                 `matrix[ i*columns + j ]`.
	 * @param rows     at least one.
	 * @param columns  at least one.
	 * @param rhs      `f`, `rows` entries.
	 * @param gamma    one non-negative, finite weight per column; zero leaves
	 *                 that column unregularised. Empty means all zero.
	 * @param relativeThreshold  singular values at or below this fraction of
	 *                 the largest count as zero, both for the rank of `M` and
	 *                 in the solve. Non-positive selects
	 *                 defaultRankThreshold().
	 *
	 * @throws std::invalid_argument naming the inconsistency: a size that does
	 *         not match, a negative or non-finite weight, a non-finite entry.
	 */
	RegularisedStep solveRegularisedLeastSquares( std::vector<double> const &matrix,
	                                              int rows, int columns,
	                                              std::vector<double> const &rhs,
	                                              std::vector<double> const &gamma = {},
	                                              double relativeThreshold = 0.0 );
}

#endif
