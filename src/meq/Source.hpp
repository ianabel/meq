#ifndef MEQ_SOURCE_HPP
#define MEQ_SOURCE_HPP

#include <memory>

#include "Profiles.hpp"

/*
 * The right hand side of the Grad-Shafranov equation.
 *
 * meq solves the fixed-boundary problem
 *
 *     -div_bar( ( 1/r ) grad_bar( psi ) ) = F( r, z, psi ) / r    in Omega
 *                                    psi  = 0                    on Gamma
 *
 *     F( r, z, psi ) := mu0 r^2 dp/dpsi + g dg/dpsi
 *
 * (Sanchez-Vizuet & Solano, CPC 235 (2019) 120-132, eqs (1)-(4)), where p( psi )
 * is the plasma pressure and g( psi ) the toroidal field function.
 *
 * A Source is exactly that F, and -- unlike the papers, which use an
 * Anderson-accelerated Picard iteration -- meq closes the nonlinearity with
 * Newton, so a Source must also supply dF/dpsi. The 1/r on the right hand side
 * belongs to the weak form, not here: F is F as written above.
 *
 * A note on the method names, since one of them is a single letter. The value is
 * f() and the derivative is dFdPsi(): F is the symbol both papers, the manual and
 * the comments in this tree use for this exact quantity, and the house rule that
 * methods start lower case demotes it to f. Spelling it value()/derivative()
 * instead would read more smoothly in isolation but would cut the one thread that
 * ties the code back to eq (2) -- and a source term that quietly disagrees with
 * eq (2) by a factor or a sign is the failure this whole file is arranged to
 * prevent. dFdPsi keeps its capital F for the same reason: it is dF/dpsi, spelled
 * the way it is written in the algorithm notes.
 *
 * As with Profiles, MFEM is deliberately absent; a mfem::Coefficient adapter
 * wraps a Source elsewhere.
 */

namespace meq
{

	/// Vacuum permeability mu0 in SI units, H/m. The conventional 4 pi x 10^-7,
	/// which since the 2019 SI redefinition is a measured rather than an exact
	/// value; the difference is 1e-10 relative and irrelevant here. Sources take
	/// it as a constructor argument so that a run in normalised units can set it
	/// to 1.
	inline constexpr double vacuumPermeability = 4.0e-7*3.14159265358979323846;

	/**
	 * The right hand side F( r, z, psi ) of the Grad-Shafranov equation, and its
	 * derivative with respect to psi.
	 *
	 * Coordinates are cylindrical ( r, z ) in metres, r > 0. z is part of the
	 * interface because a source is allowed to depend on position however it
	 * likes (a manufactured solution does); the physical MHD source does not use
	 * it.
	 */
	class Source
	{
		public:
			virtual ~Source() = default;

			/// F at ( r, z ) for the flux value psi, in the units of eq (2): the
			/// full right hand side numerator, no 1/r applied.
			virtual double f( double r, double z, double psi ) const = 0;

			/// dF/dpsi at fixed ( r, z ). Required by the Newton solve: it is the
			/// only term the source contributes to the Jacobian, and an error here
			/// does not change the converged answer, it only wrecks (or silently
			/// slows) the convergence to it. Must be the exact derivative of f().
			virtual double dFdPsi( double r, double z, double psi ) const = 0;

		protected:
			Source() = default;
			Source( Source const & ) = default;
			Source( Source && ) = default;
			Source & operator=( Source const & ) = default;
			Source & operator=( Source && ) = default;
	};

	/**
	 * The physical MHD source built from a pair of user-supplied profiles:
	 *
	 *     F( r, z, psi ) = mu0 r^2 p'( psi ) + ( g g' )( psi )
	 *
	 * Convention -- read this before wiring up a configuration file, because a
	 * factor or sign error here produces a converged but wrong equilibrium:
	 *
	 *   * The profiles handed in are the *derivative* quantities appearing in F,
	 *     not p and g themselves.
	 *       - pPrime  is dp/dpsi, the pressure gradient with respect to flux, in
	 *         Pa per Wb/rad. Multiplied by mu0 r^2 inside f().
	 *       - ggPrime is the single product g dg/dpsi, in T^2 m^2 per Wb/rad.
	 *         This is what equilibrium files tabulate ("FF'" in EQDSK, where that
	 *         F is this g); it is *not* multiplied by mu0.
	 *   * Storing the products rather than p and g is what keeps the Newton
	 *     derivative honest: no chain rule is needed anywhere, and
	 *
	 *         dF/dpsi = mu0 r^2 p''( psi ) + ( g g' )'( psi )
	 *
	 *     is just prime() of each profile. Had this class stored p and g it would
	 *     have to differentiate a product of interpolants and their derivatives,
	 *     and the Jacobian would stop matching the residual the moment either
	 *     profile's prime() disagreed with a difference of its own values.
	 *   * No sign is applied here. F is the right hand side exactly as written in
	 *     eq (2); with the Solov'ev profiles mu0 p' = -C and g g' = -A, both
	 *     negative, F comes out negative, and that is correct.
	 *   * The profiles are evaluated at the psi passed to f(), unaltered. Profiles
	 *     are tabulated against normalised flux on [ 0, 1 ] (see meq::Profile), so
	 *     whoever builds the Coefficient that feeds this class is responsible for
	 *     normalising the solver's psi the same way the table was built. This
	 *     class does not, because the normalisation moves between Newton iterates
	 *     -- psi on the magnetic axis is part of the solution -- and an MHDSource
	 *     is a fixed function of its arguments. **NormalisedMHDSource below is
	 *     the one to use when the profiles really are in normalised flux**, and it
	 *     is a different object rather than a flag on this one because psi_ax
	 *     becomes an unknown of the non-linear system rather than an input to it.
	 *
	 * Ownership: shared_ptr, so one profile can back several sources and a Source
	 * can outlive the Configuration that parsed it. Neither profile may be null.
	 */
	class MHDSource : public Source
	{
		public:
			/// Build from dp/dpsi and ( g dg/dpsi ). mu0 defaults to the SI value;
			/// pass 1 to work in normalised units. Throws std::invalid_argument if
			/// either profile is null or mu0 is not finite.
			MHDSource( std::shared_ptr<Profile const> pPrime, std::shared_ptr<Profile const> ggPrime, double mu0 = vacuumPermeability );

			/// mu0 r^2 p'( psi ) + ( g g' )( psi ). Independent of z.
			double f( double r, double z, double psi ) const override;

			/// mu0 r^2 p''( psi ) + ( g g' )'( psi ), i.e. mu0 r^2 times the
			/// pressure profile's prime() plus the g g' profile's prime().
			double dFdPsi( double r, double z, double psi ) const override;

			/// The dp/dpsi profile.
			Profile const & pPrime() const;

			/// The g dg/dpsi profile.
			Profile const & ggPrime() const;

			/// The permeability this source multiplies r^2 p' by.
			double mu0() const;

		private:
			std::shared_ptr<Profile const> pPrimeProfile;
			std::shared_ptr<Profile const> ggPrimeProfile;
			double permeability;
	};

	/**
	 * A source whose profiles are functions of NORMALISED flux,
	 *
	 *     Psi = ( psi - psi_bnd ) / ( psi_ax - psi_bnd )
	 *
	 * which is how refs/GourdainContour.pdf section V eq (39) poses them, how
	 * meq::Profile is tabulated, and how every equilibrium code specifies an
	 * equilibrium. meq solves the fixed-boundary problem with psi = 0 on Gamma,
	 * so psi_bnd is zero and Psi = psi / psi_ax throughout; free boundary makes
	 * psi_bnd an unknown too, and this interface is where that will go.
	 *
	 * WHY THIS IS NOT JUST AN MHDSource WITH A SCALED ARGUMENT, which is the
	 * whole reason it needs a class of its own. psi_ax is psi on the magnetic
	 * axis, which is to say max psi over the domain -- a GLOBAL FUNCTIONAL OF THE
	 * SOLUTION, not data. Three consequences, each of which was measured before
	 * it was believed:
	 *
	 *   * FIXING psi_ax DOES NOT APPROXIMATE THE PROBLEM, IT REPLACES IT. Hand
	 *     the solver a psi_ax the solution does not reach and the profile is
	 *     never sampled: with psi_ax = 1 on the standard box a peaked pressure
	 *     drove solutions that agreed to every digit at amplitudes 1 and 512,
	 *     because Psi never exceeded 0.0013 and a Psi^(nu-1) gradient is then
	 *     1e-9 of itself.
	 *
	 *   * THE SELF-CONSISTENT PROBLEM IS NOT THE psi_ax-PARAMETERISED ONE. With
	 *     psi_ax held fixed the equation has a small solution that Newton finds
	 *     from zero and a large one that is the equilibrium; only the large one
	 *     satisfies max psi = psi_ax. Closing the loop with an outer iteration on
	 *     psi_ax does not fix that -- the outer map has a pole beside its own
	 *     fixed point, and it falls off the branch.
	 *
	 *   * SO psi_ax BELONGS INSIDE THE RESIDUAL, as an unknown of the non-linear
	 *     system, where the Jacobian can see the non-local terms it contributes.
	 *     GradShafranovSolver::setSource( NormalisedSource &, double ) is what
	 *     does that, and the solver -- not the caller -- owns the value from then
	 *     on: it calls setNormalisation() before every residual evaluation.
	 *
	 * A source of this kind necessarily has the form F( r, z, psi ) =
	 * H( r, z, psi/psi_ax )/psi_ax, and the solver relies on nothing beyond
	 * f() and dFdPsi() answering for whatever normalisation was last set.
	 */
	class NormalisedSource : public Source
	{
		public:
			/// Set psi on the magnetic axis. The next calls to f() and dFdPsi()
			/// must answer for this value. Called by the solver once per residual
			/// evaluation, so it has to be cheap and must not allocate.
			///
			/// @throws std::invalid_argument if @a psiAxis is not finite or is
			///         zero: Psi = psi/psi_ax is undefined there, and a solver
			///         that has wandered onto psi_ax = 0 should say so rather than
			///         return infinities.
			virtual void setNormalisation( double psiAxis ) = 0;

			/// The value the next f() and dFdPsi() will use.
			virtual double normalisation() const = 0;

		protected:
			NormalisedSource() = default;
			NormalisedSource( NormalisedSource const & ) = default;
			NormalisedSource( NormalisedSource && ) = default;
			NormalisedSource & operator=( NormalisedSource const & ) = default;
			NormalisedSource & operator=( NormalisedSource && ) = default;
	};

	/**
	 * The physical MHD source with both profiles tabulated against normalised
	 * flux, which is the form meq::Profile documents and the form an equilibrium
	 * file carries:
	 *
	 *     F( r, z, psi ) = [ mu0 r^2 ( dp/dPsi )( Psi ) + ( g dg/dPsi )( Psi ) ]
	 *                      / psi_ax,          Psi = psi / psi_ax
	 *
	 * The single factor of 1/psi_ax is the chain rule and it is the whole
	 * difference between this class and MHDSource: dp/dpsi = ( dp/dPsi )/psi_ax.
	 * dF/dpsi picks up a second factor for the same reason, which is exactly
	 * where a normalisation goes missing, and SourceTests checks it against a
	 * finite difference.
	 *
	 * The profiles handed in are the DERIVATIVE quantities with respect to Psi --
	 * dp/dPsi and ( g dg/dPsi ) -- for the reason MHDSource records: storing the
	 * products rather than p and g is what keeps the Newton derivative free of a
	 * chain rule through an interpolant.
	 */
	class NormalisedMHDSource : public NormalisedSource
	{
		public:
			/// @param psiAxis  the initial normalisation. It is a starting value
			///                 and nothing more: the solver overwrites it at every
			///                 residual evaluation.
			/// @throws std::invalid_argument if either profile is null, if mu0 is
			///         not finite, or if psiAxis is not a usable normalisation.
			NormalisedMHDSource( std::shared_ptr<Profile const> pPrime,
			                     std::shared_ptr<Profile const> ggPrime,
			                     double psiAxis,
			                     double mu0 = vacuumPermeability );

			double f( double r, double z, double psi ) const override;
			double dFdPsi( double r, double z, double psi ) const override;

			void setNormalisation( double psiAxis ) override;
			double normalisation() const override;

			/// The dp/dPsi profile.
			Profile const & pPrime() const;

			/// The g dg/dPsi profile.
			Profile const & ggPrime() const;

			double mu0() const;

		private:
			std::shared_ptr<Profile const> pPrimeProfile;
			std::shared_ptr<Profile const> ggPrimeProfile;
			double psiAxisValue;
			double permeability;
	};

	/**
	 * The Solov'ev source, HDG-GS-1 eq (10):
	 *
	 *     mu0 dp/dpsi = -C,   g dg/dpsi = -A,   A + C = 1
	 *  => F( r, z, psi ) = -( ( 1 - A ) r^2 + A )
	 *
	 * The flux normalisation A + C = 1 is baked in: only A is a parameter, and
	 * C = 1 - A. Being independent of psi, F is linear in the unknown and
	 * dFdPsi() is identically zero -- Newton converges in a single step on this
	 * problem, which is precisely what makes it the first test to run.
	 *
	 * Equivalent to an MHDSource with ConstantProfile( -( 1 - A )/mu0 ) and
	 * ConstantProfile( -A ); it exists as its own class because it is the
	 * benchmark the analytic solutions in tests/analytic are built on, and
	 * because writing it out leaves nothing to get wrong.
	 */
	class SolovievSource : public Source
	{
		public:
			/// `a` is the paper's A: the coefficient of the g g' term. The pressure
			/// term then carries C = 1 - A. Throws std::invalid_argument if it is
			/// not finite.
			explicit SolovievSource( double a );

			/// -( ( 1 - A ) r^2 + A ). Independent of z and of psi.
			double f( double r, double z, double psi ) const override;

			/// Zero, exactly, for every argument: this source is linear in psi.
			double dFdPsi( double r, double z, double psi ) const override;

			/// The A this source was built with.
			double a() const;

			/// C = 1 - A, the pressure coefficient implied by the normalisation.
			double c() const;

		private:
			double aValue;
	};

}

#endif // MEQ_SOURCE_HPP
