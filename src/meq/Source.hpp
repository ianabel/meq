#ifndef MEQ_SOURCE_HPP
#define MEQ_SOURCE_HPP

#include <memory>

#include "Profiles.hpp"

/*
 * The right hand side of the Grad-Shafranov equation.
 *
 * MEQ solves the fixed-boundary problem
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
 * Anderson-accelerated Picard iteration -- MEQ closes the nonlinearity with
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
	 * equilibrium. MEQ solves the fixed-boundary problem with psi = 0 on Gamma,
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
			/// **THE BOUNDARY FLUX IS A SECOND NORMALISATION AND IT USED TO BE
			/// ASSUMED ZERO.** The profiles are functions of
			///
			///     Psi = ( psi - psi_bnd )/( psi_ax - psi_bnd )
			///
			/// and MEQ's fixed-boundary problem has psi = 0 on Gamma, so
			/// psi_bnd vanished and this class was written as Psi = psi/psi_ax.
			/// Free boundary makes it an unknown as well -- the flux at the
			/// limiter contact or the X-point -- and FREE-BOUNDARY-PLAN.md's
			/// FB-3 is a SECOND border row of the same shape as psi_ax's.
			///
			/// The one-argument form below keeps psi_bnd = 0 and is what every
			/// fixed-boundary caller wants.
			///
			/// @throws std::invalid_argument if either is not finite, or if the
			///         SPAN psi_ax - psi_bnd is zero: Psi is undefined there,
			///         and a solver that has wandered onto it should say so
			///         rather than return infinities.
			virtual void setNormalisation( double psiAxis,
			                               double psiBoundary ) = 0;

			/// psi_bnd = 0, the fixed-boundary case.
			void setNormalisation( double psiAxis )
			{
				setNormalisation( psiAxis, 0.0 );
			}

			/// The axis value the next f() and dFdPsi() will use.
			virtual double normalisation() const = 0;

			/**
			 * CONFINE THE SOURCE TO THE PLASMA, which is what makes the support
			 * MOVE with the solution and is the whole of FB-4's
			 * `chi_{Omega_p( psi )}`.
			 *
			 * The plasma is where the normalised flux is positive,
			 * `Psi = ( psi - psi_bnd )/( psi_ax - psi_bnd ) > 0`, and outside it
			 * `F` and `dF/dpsi` are ZERO. That is all the moving support needs
			 * from a source: nothing here knows where the boundary is, and the
			 * boundary is not an input -- it is wherever `psi` currently puts
			 * it, so it moves as Newton moves and converges as `psi` converges.
			 *
			 * **OFF BY DEFAULT**, so every fixed-boundary caller is untouched.
			 * There the domain IS the plasma, `Psi > 0` throughout by
			 * construction, and switching this on would change nothing except
			 * to put a branch in the inner loop.
			 *
			 * **AND IT IS THE POINTWISE TEST, NOT A CONNECTIVITY ONE.**
			 * `{ Psi > 0 }` can pick up private-flux regions beyond an X-point
			 * and pockets near the coils, which are not the plasma;
			 * FREE-BOUNDARY-PLAN.md section 5.3 records that CEDRES++ needs a
			 * connectivity test for exactly this and that
			 * meq::CriticalPointFinder is what would make one cheap. Until that
			 * exists this is a limiter plasma's support, and a diverted one's
			 * only while the search region excludes the private flux.
			 *
			 * FREE-BOUNDARY-PLAN.md section 10.3 is that note promoted to a
			 * plan, and it says which DIRECTION this is wrong in: the private
			 * flux region is the sector OPPOSITE the plasma across the saddle,
			 * so it carries Psi > 0 by value and this test switches the source
			 * ON there. A diverted run would converge and would describe a
			 * machine with a second current channel under the divertor.
			 *
			 * **A PRECONDITION, MEASURED RATHER THAN ASSUMED.** The profiles
			 * must vanish at the plasma edge -- `p'( 0 ) = 0` -- or Newton does
			 * not converge at all. With `p'( 0 ) != 0` the source JUMPS across
			 * the edge, so the assembled residual is discontinuous in the
			 * unknowns and there is no Jacobian to iterate with: measured, it
			 * fails at every degree and every mesh, from the exact solution, and
			 * under PicardThenNewton. See CLAUDE.md's *At j = 0 the question
			 * does not arise*.
			 */
			/// **VIRTUAL, AND THAT IS NOT DECORATION.** A source that WRAPS
			/// another -- meq::CoilAugmentedNormalisedSource is the one in
			/// this tree -- evaluates the plasma term through the source it
			/// holds, so it is the HELD source's flag that insidePlasma()
			/// consults. Non-virtual, a caller holding a NormalisedSource &
			/// would set the wrapper's flag, the wrapper's f() would delegate
			/// to a plasma source still unconfined, and the moving support
			/// would silently do nothing. An override forwards; this base is
			/// what an override calls to keep plasmaSupport() honest.
			virtual void setPlasmaSupport( bool confined )
			{
				confinedToPlasma = confined;
			}

			/// Whether setPlasmaSupport() is on.
			bool plasmaSupport() const
			{
				return confinedToPlasma;
			}

			/**
			 * True where the plasma is: `Psi > 0`, with `Psi` built from the
			 * normalisation the source currently carries.
			 *
			 * Always true when setPlasmaSupport() is off, so a subclass may call
			 * it unconditionally and a fixed-boundary solve pays one comparison.
			 *
			 * THE SIGN OF THE SPAN IS NOT ASSUMED. `psi_ax - psi_bnd` is
			 * negative wherever `F` is single-signed negative -- every Solov'ev
			 * fixture in this tree has its magnetic axis at an interior MINIMUM,
			 * which CriticalPoints.hpp records as the reason findAxis() seeds
			 * from both nodal extremes. So the test is on the PRODUCT rather
			 * than on the difference, and it is the same test either way round.
			 */
			bool insidePlasma( double psi ) const
			{
				if ( !confinedToPlasma )
					return true;
				double const span = normalisation() - boundaryNormalisation();
				return ( psi - boundaryNormalisation() )*span > 0.0;
			}

			/**
			 * dF/d(psi_ax) and dF/d(psi_bnd) AT FIXED psi, analytically.
			 *
			 * WHY THIS EXISTS. The bordered Newton's COLUMN is dR/ds, the
			 * derivative of the residual with respect to a normalisation, and
			 * solveWithNormalisation() obtains it by a CENTRAL DIFFERENCE of two
			 * full residual evaluations. The row is exact under NPC -- psi_ax's
			 * is -e_j and psi_bnd's is too -- so the column is the only
			 * differenced thing left, and it is what floors the iteration:
			 * measured on the half-disc, a coupled solve descends to about 3e-09
			 * and then sits there for as many iterations as it is given, which
			 * is the difference's own accuracy and not the discretisation's.
			 *
			 * AND IT IS WORSE THAN A FLOOR WHEN THE SUPPORT MOVES.
			 * setPlasmaSupport() makes the edge a function of the normalisation,
			 * so perturbing s by a step moves the edge across quadrature points,
			 * and a central difference then STRADDLES a kink rather than
			 * measuring a derivative. An analytic value straddles nothing: it is
			 * evaluated pointwise at the current state, and outside the plasma
			 * it is zero for the same reason f() is.
			 *
			 * @return false if the source cannot supply these, in which case the
			 *         solver differences the column as before. Deliberately NOT
			 *         pure: this is an optimisation of the Jacobian and never of
			 *         the answer, so a source that has not implemented it must
			 *         keep working rather than fail to compile.
			 *
			 * PRECONDITION, AND IT IS THE SAME ONE dFdPsi() CARRIES. With a
			 * moving support the true derivative of the assembled residual picks
			 * up a SURFACE term where the edge sweeps, and that term vanishes
			 * exactly when the profiles vanish at the edge. So this is the
			 * derivative of the assembled residual when p'( 0 ) = 0 and is
			 * missing a term when it is not -- which is the condition
			 * setPlasmaSupport() already documents as its precondition.
			 */
			virtual bool normalisationDerivatives( double r, double z, double psi,
			                                       double &dFdAxis,
			                                       double &dFdBoundary ) const
			{
				(void)r; (void)z; (void)psi;
				(void)dFdAxis; (void)dFdBoundary;
				return false;
			}

			/**
			 * SCALE THE PLASMA TERM, so that the TOTAL PLASMA CURRENT can be
			 * prescribed instead of the profile amplitude.
			 *
			 * WHY THIS EXISTS, AND IT IS NOT A CONVENIENCE. With the amplitude
			 * fixed, a confined equilibrium is a non-linear EIGENVALUE problem:
			 * writing `p' ~ A Psi^j`, the equation reduces to
			 * `Delta* u = -Lambda u^j (...)` with `Lambda = A/span^2`, so
			 * `Lambda` has to be an eigenvalue of the linearised operator ON THE
			 * PLASMA REGION -- and the region is itself unknown. Scaling `A`
			 * therefore changes nothing, because `span` moves with `sqrt( A )`
			 * and the reaction ratio is amplitude-independent; measured, and the
			 * same statement CLAUDE.md records for the high-beta source.
			 *
			 * Prescribing `I_p` and solving for this scale turns that balance
			 * into an ordinary unknown of the same bordered Newton. It is what
			 * CEDRES++ and FreeGS both do -- FREE-BOUNDARY-PLAN.md section 7.11
			 * names `ConstrainBetapIp` as FB-6's reference model -- and it is
			 * why every production free-boundary code asks for a current rather
			 * than for an amplitude.
			 *
			 * VIRTUAL FOR THE REASON setPlasmaSupport() IS: a wrapper must
			 * forward it to the source that evaluates the profiles, and
			 * meq::CoilAugmentedNormalisedSource does. **A coil is not scaled**
			 * -- its current is amperes and is prescribed input, so the scale
			 * multiplies the plasma term alone.
			 */
			virtual void setCurrentScale( double scale )
			{
				currentScaleValue = scale;
			}

			/// The scale in force. One unless it has been set.
			virtual double currentScale() const
			{
				return currentScaleValue;
			}

			/**
			 * The part of `f()` that carries setCurrentScale()'s factor, and its
			 * `psi`-derivative.
			 *
			 * These are `f()` and `dFdPsi()` for an ordinary plasma source and
			 * are NOT for a wrapped one: a coil-augmented source's `f()` is the
			 * sum, and the current constraint is about the plasma alone. So a
			 * wrapper overrides these to forward to what it holds, exactly as it
			 * forwards setCurrentScale().
			 */
			virtual double scaledF( double r, double z, double psi ) const
			{
				return f( r, z, psi );
			}

			/// @see scaledF
			virtual double scaledDFdPsi( double r, double z, double psi ) const
			{
				return dFdPsi( r, z, psi );
			}

			/// And the boundary value; zero unless it has been set.
			virtual double boundaryNormalisation() const = 0;

		private:
			/// setPlasmaSupport(). Off by default; see it for why.
			bool confinedToPlasma = false;

			/// setCurrentScale(). One unless a current is prescribed.
			double currentScaleValue = 1.0;

		public:

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

			/// Analytic, so the bordered Newton need not difference its columns.
			bool normalisationDerivatives( double r, double z, double psi,
			                               double &dFdAxis,
			                               double &dFdBoundary ) const override;

			void setNormalisation( double psiAxis, double psiBoundary ) override;
			using NormalisedSource::setNormalisation;
			double normalisation() const override;
			double boundaryNormalisation() const override;

			/// The dp/dPsi profile.
			Profile const & pPrime() const;

			/// The g dg/dPsi profile.
			Profile const & ggPrime() const;

			double mu0() const;

		private:
			std::shared_ptr<Profile const> pPrimeProfile;
			std::shared_ptr<Profile const> ggPrimeProfile;
			double psiAxisValue;
			double psiBoundaryValue = 0.0;
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
