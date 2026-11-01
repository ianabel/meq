#ifndef MEQ_ANISOTROPY_HPP
#define MEQ_ANISOTROPY_HPP

#include <cmath>
#include <stdexcept>

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
	 * argument handed: the solver applies 1/( psi_ax - psi_bnd ) itself. A
	 * Delta that reads B is meq::FieldAnisotropy, below.
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

	/**
	 * A Delta THAT READS B, ANISOTROPY-PLAN.md AN-6a: Delta = mu0 Delta_p/B^2
	 * with B the iterate's |B|, the form a kinetic closure hands an
	 * equilibrium (mirrorfp's, through geq's KineticClosure). The class hands
	 * the solver only the numerator,
	 *
	 *     pressure() = mu0 sum_s ( p_par,s - p_perp,s ),
	 *
	 * in tesla squared, of position and of the flux argument() names, as an
	 * Anisotropy's Delta is. The solver forms Delta.
	 *
	 * WITH g = 0, B = |q|/( 1 - Delta ), q the total flux the law carries, so
	 * at each law point Delta solves Delta |q|^2 = mu0 Delta_p ( 1 - Delta )^2,
	 * through fieldDelta() below. Every such Delta is a flux law, whatever
	 * argument() says. The mirror has g = 0; a toroidal field would make B
	 * implicit in Delta through Ibar/( 1 - Delta ) too, and is not built.
	 *
	 * delta() and dDeltaDPsi() need B, which only the solver has, so they
	 * throw: nothing should reach them.
	 */
	class FieldAnisotropy : public Anisotropy
	{
		public:
			/// mu0 sum_s ( p_par,s - p_perp,s ) at ( R, z ) and the argument, in
			/// tesla squared: Delta times B^2.
			virtual double pressure( double radius, double z, double psi ) const = 0;

			/// d pressure()/d( argument ) at fixed ( R, z ), with the same
			/// convention as Anisotropy::dDeltaDPsi().
			virtual double dPressureDPsi( double radius, double z, double psi ) const = 0;

			double delta( double, double, double ) const final
			{
				throw std::logic_error( "meq::FieldAnisotropy::delta: Delta reads B, which only the solver has; it is formed by meq::fieldDelta()" );
			}

			double dDeltaDPsi( double, double, double ) const final
			{
				throw std::logic_error( "meq::FieldAnisotropy::dDeltaDPsi: Delta reads B, which only the solver has; it is formed by meq::fieldDelta()" );
			}
	};

	/// fieldDelta()'s answer.
	struct FieldDelta
	{
		/// Delta, or where no root exists the value it is held at.
		double delta = 0.0;

		/// d Delta/d pressure at fixed |q|^2.
		double dDeltaDPressure = 0.0;

		/// d Delta/d |q|^2 at fixed pressure.
		double dDeltaDFluxSquared = 0.0;

		/// False where no root exists; delta is then held and the slopes are 0.
		bool exists = true;
	};

	/**
	 * Delta FROM ITS NUMERATOR AND THE FLUX, meq::FieldAnisotropy's: the root
	 * of
	 *
	 *     Delta w^2 = a ( 1 - Delta )^2,   a = mu0 Delta_p,  w = |q|,
	 *
	 * that vanishes with a,
	 *
	 *     Delta = 2a/( 2a + w^2 + w sqrt( w^2 + 4a ) ),
	 *
	 * written so that nothing cancels as a goes to zero. Implicitly
	 * differentiated, with J = w^2 + 2a( 1 - Delta ) = w sqrt( w^2 + 4a ),
	 *
	 *     dDelta/da = ( 1 - Delta )^2/J,   dDelta/dw^2 = -Delta/J.
	 *
	 * FOR a > 0 THE ROOT LIES IN ( 0, 1 ) AT EVERY w > 0: at fixed flux the
	 * firehose limit is approached, as B falls to sqrt( a ), and never reached.
	 * FOR a < 0 IT EXISTS ONLY WHILE w^2 >= -4a. There Delta = -1 and J = 0:
	 * the flux w = B - |a|/B has a minimum over B, at B^2 = -a, and no B
	 * carries less flux. Past that fold, and at w = 0 with a > 0, exists is
	 * false, Delta is held at the fold's -1 or at 0, and the slopes are zero.
	 * The solver refuses an ANSWER with either; an iterate may pass through.
	 */
	inline FieldDelta fieldDelta( double pressure, double fluxSquared )
	{
		FieldDelta result;
		if ( pressure == 0.0 )
		{
			if ( fluxSquared > 0.0 )
				result.dDeltaDPressure = 1.0/fluxSquared;
			return result;
		}
		double const discriminant = fluxSquared + 4.0*pressure;
		if ( !( discriminant > 0.0 ) || !( fluxSquared > 0.0 ) )
		{
			result.delta = pressure < 0.0 ? -1.0 : 0.0;
			result.exists = false;
			return result;
		}
		double const flux = std::sqrt( fluxSquared );
		double const fold = flux*std::sqrt( discriminant );
		result.delta = 2.0*pressure/( 2.0*pressure + fluxSquared + fold );
		double const tension = 1.0 - result.delta;
		double const jacobian = fluxSquared + 2.0*pressure*tension;
		result.dDeltaDPressure = tension*tension/jacobian;
		result.dDeltaDFluxSquared = -result.delta/jacobian;
		return result;
	}
}

#endif
