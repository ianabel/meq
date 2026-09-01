#ifndef MEQ_ROTATINGSOURCE_HPP
#define MEQ_ROTATINGSOURCE_HPP

#include "Profiles.hpp"
#include "Source.hpp"

#include <array>
#include <cstddef>
#include <memory>
#include <vector>

/*
 * The Grad-Shafranov source for a plasma in sonic toroidal rotation.
 *
 * The equation is refs/RotatingGK.pdf eq (136) -- Abel, Plunk, Wang, Barnes,
 * Cowley, Dorland & Schekochihin, Rep. Prog. Phys. 76 (2013) 116201 -- closed by
 * its (96) for the poloidal density variation and (97) for the electrostatic
 * potential phi_0 that holds quasineutrality against it. FLOW-PLAN.md is the
 * design; this file implements FL-0 to FL-3 of it.
 *
 * WHAT ROTATION CHANGES, AND IT IS ONE THING: the density is no longer a flux
 * function. Centrifugal force sweeps heavy species to the outboard side, and an
 * electrostatic potential arises to stop that separating the charges. So
 *
 *     n_s( r, psi ) = n_s0( psi ) exp[ m_s omega^2 ( r^2 - rRef^2 )/2T_s
 *                                      - Z_s e phi_0/T_s ]          (96)
 *     Sum_s Z_s n_s( r, psi ) = 0                                   (97)
 *
 * and (97) is what determines phi_0 at each ( r, psi ). Everything else about
 * meq is untouched: the operator, the discretisation, tau, the hybridization,
 * the estimator, the adaptive loop and the curved boundary all stay as they are.
 * meq::Source's signature already carries r, which is the whole reason a
 * non-flux-function density needs no interface change anywhere in the solver.
 *
 * THE GAUGE IS FREE AND THIS CLASS TAKES THE LOCAL ONE. RoPP fixes phi_0 by
 * <phi_0>_psi = 0 and says at its (59) that this is a convention rather than a
 * closure -- "we can add any function of psi to it". That choice would make
 * phi_0 proportional to r^2 - <R^2>_psi, a flux-surface average of the unknown
 * on every surface, which meq has no machinery for and which buys nothing
 * physical. This class instead pins
 *
 *     phi_0( rRef, psi ) = 0
 *
 * with rRef a CONSTANT -- the geometric axis, given to the constructor, not the
 * magnetic axis and not an average. Two consequences, both good: F is a
 * pointwise function of ( r, z, psi ), and n_s0 is then the PHYSICAL density of
 * species s on the curve r = rRef, which is a quantity a user can state and
 * another code can be compared against. Li & Zhu (Comput. Phys. Commun. 260
 * (2021) 107264) make the same choice independently, referencing their exponent
 * to the magnetic axis.
 *
 * THE PRICE OF THE GAUGE IS THAT TWO SETS OF n_s0 DIFFERING BY IT DESCRIBE THE
 * SAME PLASMA. Anyone comparing profiles against RoPP's own notation, against a
 * gyrokinetic code or against a transport code has to know which convention each
 * is in. That is the one part of the gauge freedom that is not a simplification.
 *
 * WHAT F IS. Differentiating p = Sum_s n_s T_s at fixed r, the dphi_0/dpsi terms
 * collect into -e ( dphi_0/dpsi ) Sum_s Z_s n_s, which vanishes IDENTICALLY by
 * (97). What is left is exactly (136)'s brace plus its omega omega' term, so
 *
 *     F( r, z, psi ) = mu0 r^2 dp/dpsi|_r + g g'
 *
 * which is meq::MHDSource's shape with an r-dependent p. The residual therefore
 * needs phi_0 but never its psi-derivative; only the Jacobian does. Checked two
 * ways: it is RoPP's own force balance (128), and at omega -> 0 it gives RoPP
 * (243), which is the static equation meq already solves.
 *
 * UNITS AND SIGN, pinned deliberately because two sign errors have already been
 * found in the papers meq does follow. RoPP is Gaussian, with the 4 pi of (136)
 * and the c of (135); meq is SI, so 4 pi -> mu0 and RoPP's I is meq's g = r B_phi.
 * RoPP (135) is j.grad(phi) = -( c/4 pi r^2 ) Delta* psi, so (136) reads
 * Delta* psi = -4 pi r^2 { ... }, so F = -Delta* psi = +mu0 r^2 { ... } -- the
 * same positive convention meq::MHDSource uses. RoPP's psi is poloidal flux per
 * radian, from its (31) B = I grad(phi) + grad(psi) x grad(phi), which is what
 * meq's g g' in T^2 m^2 per Wb/rad already assumes and what EQDSK tabulates.
 */

namespace meq
{

	/// The most species meq::RotatingSource will carry. A cap rather than a
	/// vector so that the per-quadrature-point work allocates nothing and needs
	/// no mutable scratch, which a source evaluated from a threaded assembly
	/// must not have. Eight covers a main ion, a bulk impurity, a seeded
	/// impurity, an alpha population and electrons with room to spare.
	inline constexpr std::size_t maxSpecies = 8;

	/**
	 * One plasma species: its constants, its temperature and its density on the
	 * reference curve.
	 *
	 * A plain value type with no MFEM in it, for the same reason meq::Profile and
	 * meq::Source have none -- so that the physics can be unit tested without the
	 * library.
	 */
	struct Species
	{
		/// Particle mass in kg. Must be finite and positive.
		double mass = 0.0;

		/// Charge in units of e, signed and dimensionless: +1 for a proton, -1
		/// for an electron, +6 for fully stripped carbon. Must be non-zero.
		double charge = 0.0;

		/// T_s( psi ) in JOULES, not eV and not keV. Must not be null, and must
		/// be positive everywhere the solve reaches -- a species with a
		/// non-positive temperature has no Maxwellian.
		std::shared_ptr<Profile const> temperature;

		/// n_s0( psi ) in m^-3: the density of this species ON THE CURVE
		/// r = rRef, which in the local gauge of this file is what N_s of RoPP
		/// (96) reduces to. Must not be null.
		std::shared_ptr<Profile const> density;
	};

	/**
	 * Sum_s Z_s n_s0( psi ), which charge neutrality on r = rRef requires to
	 * vanish.
	 *
	 * Exposed because it is the quantity meq::RotatingSource's constructor
	 * checks, so a caller can see why a set was refused, and because it is what a
	 * configuration layer should drive to zero rather than asking the user to.
	 */
	double chargeNeutralityResidual( std::vector<Species> const & species, double psi );

	/**
	 * The density profile that charge neutrality implies for one species, given
	 * all the others.
	 *
	 * WHY THIS EXISTS RATHER THAN A CHECK ALONE. Fixing the gauge removes exactly
	 * one function's worth of freedom from the set of densities -- the
	 * transformation phi_0 -> phi_0 + delta( psi ), n_s0 -> n_s0 exp( Z_s e
	 * delta/T_s ) leaves every physical quantity alone -- so for n species there
	 * are n - 1 independent density flux functions. For two species that is one,
	 * which is the familiar statement. Asking a user to supply n profiles that
	 * happen to satisfy a constraint invites them to supply n that do not.
	 *
	 * The returned profile is -( 1/Z_index ) Sum_{s != index} Z_s n_s0, and its
	 * prime() and doublePrime() are the matching combinations, so it is exact at
	 * every level rather than differenced.
	 *
	 * The nominated species' own density is never read, so it may be null when
	 * this is called -- which is what lets a caller fill the set in one pass.
	 *
	 * @throws std::out_of_range if @a index is not a species.
	 * @throws std::invalid_argument if the nominated species has zero charge, if
	 *         fewer than two species are given, or if any density OTHER than the
	 *         nominated one is null.
	 */
	std::shared_ptr<Profile const> neutralisingDensity( std::vector<Species> const & species, std::size_t index );

	/**
	 * The rotating Grad-Shafranov source.
	 *
	 * WHY THIS IS NOT A GENERALISED meq::MHDSource. That class stores the
	 * pre-multiplied PRODUCTS p' and g g' precisely so that no chain rule exists
	 * to get wrong: F is one profile evaluation and dF/dpsi is one prime(). The
	 * trick is unavailable here, because p depends on the density, the
	 * temperatures and omega separately and on r as well, so there is no product
	 * to pre-store. Two classes with different invariants, not one class with a
	 * flag -- and meq::MHDSource remains the omega = 0 path.
	 *
	 * IT COSTS A THIRD DERIVATIVE LEVEL, which is why meq::Profile grew
	 * doublePrime(). F is already dp/dpsi, so the Jacobian is d2p/dpsi2 and every
	 * input profile is asked for two derivatives. There is no reparametrisation
	 * that avoids it; see Profile::doublePrime.
	 *
	 * TWO CLOSURES, AND THE CHOICE IS NOT NUMERICAL. With exactly two species
	 * (97) is linear in phi_0 after taking logarithms, so the whole closure is a
	 * closed form: no root find, no inner tolerance, no implicit
	 * differentiation. Three or more needs a safeguarded scalar Newton per
	 * evaluation point and phi_0's two psi-derivatives by implicit
	 * differentiation. Both are implemented and Closure::Automatic picks the
	 * cheap one when it applies; Closure::RootFind can be forced at two species
	 * so that a test can drive the same problem down both paths and require the
	 * same answer, which is what says the general path is right.
	 */
	class RotatingSource : public Source
	{
		public:
			/// How (97) is solved for phi_0.
			enum class Closure
			{
				Automatic,  ///< ClosedForm at two species, RootFind above.
				ClosedForm, ///< Two species only; refused for more.
				RootFind    ///< Any number, and the cross-check at two.
			};

			/**
			 * @param species  exactly two, of opposite charge sign, satisfying
			 *                 charge neutrality on r = rRef.
			 * @param omega    the rigid rotation frequency omega( psi ) in rad/s.
			 *                 May be null, which means omega = 0 and reduces this
			 *                 source to meq::MHDSource's equation.
			 * @param ggPrime  the single product g dg/dpsi, in T^2 m^2 per Wb/rad,
			 *                 exactly as meq::MHDSource takes it. Must not be null.
			 * @param referenceRadius  rRef in metres, where phi_0 vanishes and
			 *                 where each n_s0 is the physical density. Must be
			 *                 finite and positive.
			 * @param mu0      so that a run in normalised units can set it to 1.
			 *
			 * @param closure  which route through (97); see Closure.
			 *
			 * @throws std::invalid_argument for a null profile, a non-finite or
			 *         non-positive mass, radius or mu0, a zero charge, fewer than
			 *         two species or more than maxSpecies, a species set whose
			 *         charges are all one sign, or Closure::ClosedForm asked for
			 *         with more than two species.
			 * @throws std::invalid_argument if charge neutrality on r = rRef is
			 *         violated. It is checked by sampling n_s0 over [ 0, 1 ], the
			 *         range meq::Profile documents itself on, because the closed
			 *         form below is derived from Z_1 n_10 = -Z_2 n_20 and is
			 *         simply wrong without it. Use neutralisingDensity() to build
			 *         a conforming set rather than hand-balancing one.
			 */
			RotatingSource( std::vector<Species> species,
				std::shared_ptr<Profile const> omega,
				std::shared_ptr<Profile const> ggPrime,
				double referenceRadius,
				double mu0 = vacuumPermeability,
				Closure closure = Closure::Automatic );

			/// F at ( r, z ) for the flux value psi, in the units of
			/// meq::Source::f: mu0 r^2 dp/dpsi|_r + g g', no 1/r applied.
			double f( double r, double z, double psi ) const override;

			/// dF/dpsi at fixed ( r, z ), which is mu0 r^2 d2p/dpsi2|_r plus the
			/// derivative of g g'. Exact, and checked against a central difference
			/// of f() over a Mach sweep in RotatingSourceTests.
			double dFdPsi( double r, double z, double psi ) const override;

			// ---- the closure, exposed because it is what FL-1 asserts on ----

			/// e phi_0( r, psi ) in JOULES -- the potential of (97) times the
			/// elementary charge, which is the combination that appears in every
			/// exponent and the one that carries no factor of e. Zero at r = rRef
			/// by construction, exactly.
			double potential( double r, double psi ) const;

			/// n_s( r, psi ) in m^-3, from (96). At r = rRef this is n_s0( psi ).
			/// @throws std::out_of_range if @a index is not a species.
			double density( std::size_t index, double r, double psi ) const;

			/// d( e phi_0 )/dpsi at fixed r, by implicit differentiation of (97)
			/// rather than by differencing it. Differencing an inner root find
			/// from outside gives a derivative whose accuracy is the inner
			/// tolerance, which is how a Newton iteration silently degrades from
			/// quadratic to linear with no wrong answer anywhere to point at.
			double dPotentialDPsi( double r, double psi ) const;

			/// The total pressure Sum_s n_s T_s at ( r, psi ), in Pa. This is the
			/// p whose psi-derivative F is built from.
			double pressure( double r, double psi ) const;

			/// dp/dpsi at fixed r, in Pa per Wb/rad. F is mu0 r^2 times this plus
			/// g g'; exposed separately so a test can see which half is wrong.
			double dPressureDPsi( double r, double psi ) const;

			/// The exponent of (96) for one species. AT TWO SPECIES ALL SPECIES
			/// SHARE IT, which is what makes Sum_s Z_s n_s vanish at every r once
			/// it vanishes at rRef; above two they differ, and the cancellation
			/// is the root find's job instead. Dimensionless, and equal to M^2/2
			/// at the outboard edge in the usual Mach-number sense.
			/// @throws std::out_of_range if @a index is not a species.
			double densityExponent( std::size_t index, double r, double psi ) const;

			// ---- accessors ----

			std::vector<Species> const & species() const;
			Closure closure() const;
			Profile const * omega() const;
			Profile const & ggPrime() const;
			double referenceRadius() const;
			double mu0() const;

		private:
			// Everything the source needs at one ( r, psi ), in one pass, on the
			// stack. f() and dFdPsi() are called at every quadrature point of
			// every element of every residual evaluation, so this allocates
			// nothing and holds no state between calls.
			struct State
			{
				double potential;             // e phi_0
				double potentialPrime;        // and its two psi-derivatives
				double potentialDoublePrime;
				std::array<double, maxSpecies> exponent;          // A_s of (96)
				std::array<double, maxSpecies> exponentPrime;     // dA_s/dpsi
				std::array<double, maxSpecies> exponentDoublePrime;
			};

			State stateAt( double r, double psi ) const;
			State closedFormState( double r, double psi ) const;
			State rootFindState( double r, double psi ) const;

			// p, dp/dpsi and d2p/dpsi2 at fixed r, from a State.
			void pressureFrom( State const & state, double psi,
				double & p, double & pPrime, double & pDoublePrime ) const;

			std::vector<Species> speciesData;
			std::shared_ptr<Profile const> omegaProfile;
			std::shared_ptr<Profile const> ggPrimeProfile;
			double rRef;
			double permeability;
			Closure closureChoice;
	};

	/**
	 * The rotating source with its profiles given against NORMALISED flux.
	 *
	 * Psi = ( psi - psi_bnd )/( psi_ax - psi_bnd ), and meq is fixed boundary
	 * with psi = 0 on Gamma, so Psi = psi/psi_ax -- which makes psi_ax a
	 * functional of the solution and therefore an unknown of the non-linear
	 * system, not an input. meq::NormalisedSource is the interface the solver
	 * closes by a bordered Newton; see CLAUDE.md, "Newton, and the obligation it
	 * creates", before using this.
	 *
	 * IT IS A WRAPPER AND NOT A REIMPLEMENTATION, which is the whole design. A
	 * normalised source has the form F( r, z, psi ) = H( r, z, psi/psi_ax )/psi_ax,
	 * so with the profiles read at Psi the answer is the unnormalised source's,
	 * divided by psi_ax once for f() and TWICE for dFdPsi() -- once for the
	 * profile argument and once for the derivative itself. Every property
	 * RotatingSourceTests establishes for meq::RotatingSource therefore carries
	 * over unchanged, and there is no second copy of the closure to keep in step.
	 *
	 * The profiles mean what they always meant, read at Psi: n_s0( Psi ) is the
	 * density on r = rRef, T_s( Psi ) is in Joules, omega( Psi ) is in rad/s, and
	 * ggPrime is d( g^2/2 )/dPsi -- note dPsi, not dpsi, which is what the second
	 * factor of 1/psi_ax accounts for.
	 */
	class NormalisedRotatingSource : public NormalisedSource
	{
		public:
			/// Arguments as meq::RotatingSource's, plus the starting psi_ax.
			/// @throws std::invalid_argument for anything RotatingSource refuses,
			///         or a psi_ax that is not finite or is zero.
			NormalisedRotatingSource( std::vector<Species> species,
				std::shared_ptr<Profile const> omega,
				std::shared_ptr<Profile const> ggPrime,
				double referenceRadius,
				double psiAxis,
				double mu0 = vacuumPermeability,
				RotatingSource::Closure closure = RotatingSource::Closure::Automatic );

			double f( double r, double z, double psi ) const override;
			double dFdPsi( double r, double z, double psi ) const override;

			/// @throws std::invalid_argument if psi_ax is not finite or is zero:
			///         Psi = psi/psi_ax is undefined there, and a solver that has
			///         wandered onto psi_ax = 0 should say so rather than return
			///         infinities.
			void setNormalisation( double psiAxis ) override;
			double normalisation() const override;

			// Diagnostics, taking PHYSICAL psi and converting internally, so that
			// a caller never has to remember which of the two it is holding.

			double potential( double r, double psi ) const;
			double density( std::size_t index, double r, double psi ) const;
			double pressure( double r, double psi ) const;

			/// The source underneath, whose arguments are Psi rather than psi.
			RotatingSource const & unnormalised() const;

		private:
			RotatingSource inner;
			double psiAxisValue;
	};

}

#endif // MEQ_ROTATINGSOURCE_HPP
