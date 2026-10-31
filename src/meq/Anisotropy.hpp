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
	 * WHAT delta() IS HANDED is argument()'s choice:
	 *
	 *   Position        nothing; Delta( R, z ). Derive from meq::PositionAnisotropy,
	 *                   which keeps the flux mass linear (AN-1)
	 *   Flux            the physical flux psi the solver carries (AN-2a)
	 *   NormalisedFlux  Psi = ( psi - psi_bnd )/( psi_ax - psi_bnd ), from the
	 *                   normalised source's labels as the bordered Newton moves
	 *                   them (AN-2b) -- the label a mirror's profiles are
	 *                   written in
	 *
	 * Under the last two the flux mass is a non-linear flux law with the
	 * derivative below in its Jacobian, and dDeltaDPsi() is with respect to the
	 * argument handed: the solver applies 1/( psi_ax - psi_bnd ) itself. The
	 * dependence on B is AN-6's.
	 *
	 * MFEM-free on purpose, as meq::Source is, so a closed form is unit-testable
	 * without the library. EVERY METHOD MUST BE REENTRANT: the solver evaluates
	 * them from threaded element loops.
	 */
	class Anisotropy
	{
		public:
			enum class Argument
			{
				Position,
				Flux,
				NormalisedFlux
			};

			virtual ~Anisotropy() = default;

			/// Delta at ( R, z ) and the flux argument() names, dimensionless.
			/// 1 - Delta must be positive wherever the solve reaches: Delta = 1 is
			/// the firehose limit, where the operator stops being elliptic, and
			/// the solver refuses a Delta at or beyond it rather than clipping.
			virtual double delta( double radius, double z, double psi ) const = 0;

			/// d Delta/d( argument ) at fixed ( R, z ). What the Newton Jacobian's
			/// flux-potential block carries; a wrong one costs quadratic
			/// convergence and not the answer, so it is checked by its own test.
			virtual double dDeltaDPsi( double radius, double z, double psi ) const = 0;

			/// What delta() reads. Position selects the linear flux mass, and a
			/// Position answer from a Delta that does read its argument is a
			/// different equation solved without complaint.
			virtual Argument argument() const = 0;
	};

	/// A Delta of position only, AN-1's: the flux mass stays a bilinear form.
	class PositionAnisotropy : public Anisotropy
	{
		public:
			/// Delta at ( R, z ).
			virtual double delta( double radius, double z ) const = 0;

			double delta( double radius, double z, double ) const final
			{
				return delta( radius, z );
			}

			double dDeltaDPsi( double, double, double ) const final
			{
				return 0.0;
			}

			Argument argument() const final
			{
				return Argument::Position;
			}
	};
}

#endif
