#ifndef MEQ_ANISOTROPY_HPP
#define MEQ_ANISOTROPY_HPP

namespace meq
{
	/**
	 * THE PRESSURE ANISOTROPY Delta = mu0 sum_s ( p_par,s - p_perp,s )/B^2, as
	 * the operator sees it. refs/AlphaGK.pdf eq (6.18), ANISOTROPY-PLAN.md.
	 *
	 * The anisotropic Grad-Shafranov equation, in MEQ's convention, is
	 *
	 *     -div_bar( ( 1 - Delta )/R grad_bar psi ) = F/R,
	 *
	 * so Delta enters the OPERATOR and not the source: the diffusion coefficient
	 * 1/R becomes ( 1 - Delta )/R, and meq::GradShafranovSolver::setAnisotropy()
	 * carries it in the flux mass. The sum is over every species: a Maxwellian
	 * contributes Delta_s = 0, and nothing in the equation tells a thermal species
	 * from a fast one.
	 *
	 * AN-1: Delta IS A FUNCTION OF POSITION ONLY, so the flux mass stays linear.
	 * Its dependence on psi, through Psi and the labels, is AN-2's, and on B,
	 * AN-6's.
	 *
	 * MFEM-free on purpose, as meq::Source is, so a closed form is unit-testable
	 * without the library.
	 */
	class Anisotropy
	{
		public:
			virtual ~Anisotropy() = default;

			/// Delta at ( R, z ), dimensionless. 1 - Delta must be positive
			/// wherever the solve reaches: Delta = 1 is the firehose limit, where
			/// the operator stops being elliptic, and the solver refuses a Delta
			/// at or beyond it rather than clipping.
			///
			/// MUST BE REENTRANT: the solver evaluates it from threaded element
			/// loops.
			virtual double delta( double radius, double z ) const = 0;
	};
}

#endif
