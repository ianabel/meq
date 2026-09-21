#ifndef MEQ_CONDUCTORFIELD_HPP
#define MEQ_CONDUCTORFIELD_HPP

#include <cstddef>
#include <vector>

#include "Coils.hpp"

namespace meq
{
	/**
	 * `psi_c`: the field of the conductors ALONE, evaluated exactly, for the
	 * split `psi = psi_c + psi_p` of `COIL-SUBTRACTION-PLAN.md`.
	 *
	 * `Delta*` is linear and a forward solve's conductor currents are
	 * prescribed inputs, so the coil term can be taken off the right-hand side
	 * entirely and what is solved for is the remainder
	 *
	 *     Delta* psi_p = -mu0 r J_plasma( psi_c + psi_p )
	 *
	 * with `psi_c` supplied here. **The conductors then need not be in the mesh
	 * at all**, which is the point of the plan.
	 *
	 * **A FILAMENT IS THE CONDUCTOR THIS EXISTS FOR.** A filament
	 * is the one conductor MEQ structurally CANNOT carry any other way: a point
	 * source has no finite-element representation as a current density, `psi`
	 * near it is logarithmic, and `meq::CoilSet::f()` has nothing to add for it
	 * because the current density is infinite on a set of measure zero. So
	 * subtraction is not one route to filament support, it is the only one --
	 * and the remainder is well behaved for exactly the reason the subtraction
	 * works: `psi_c` carries the whole logarithm, so `psi_p` sees a bounded
	 * right-hand side supported on the plasma alone.
	 *
	 * **RECTANGLES ARE HERE TOO (CS-1b), AND THEY ARE A DIFFERENT KIND OF
	 * CONDUCTOR RATHER THAN A BIGGER FILAMENT.** A `meq::Coil` carries a uniform
	 * current density over an area, so its field is a quadrature of the filament
	 * kernel over the cross-section and is **finite everywhere, including
	 * inside the coil** -- `meq::coilPsi()` says so in those words. Two
	 * consequences, and both are contracts rather than conveniences:
	 *
	 *   * **coincides() is about FILAMENTS ONLY.** A rectangle has no line
	 *     singularity for a mesh point to land on, so there is nothing to
	 *     refuse, and refusing a node inside a coil would reject the ordinary
	 *     configuration this plan exists to make cheap.
	 *   * a rectangle CAN be carried as a domain source instead, through
	 *     `meq::CoilAugmentedSource`, and for a rectangle the split is a
	 *     performance and accuracy choice rather than the only option. For a
	 *     filament there is no alternative at all.
	 *
	 * The quadrature order is forwarded to the underlying `meq::CoilSet` and is
	 * raisable, which is what `COIL-SUBTRACTION-PLAN.md` §7.2's replacement for
	 * CS-5 needs: a reference field built far beyond what a solve would use.
	 *
	 * **THIS CLASS IS MFEM-FREE AND THAT DECIDES ITS INTERFACE.** It lives
	 * beside `meq::Coils` in the half of `src/meq` that CI can build, so it
	 * cannot take an `mfem::Mesh` -- and the mesh check below is therefore a
	 * per-POINT predicate that the caller, which owns the mesh, loops. The
	 * physics stays unit-testable and the mesh walking stays where meshes live.
	 *
	 * `mu0` is a constructor argument for the reason `meq::Source` gives: a run
	 * in normalised units sets it to 1, and it should be visible in one place
	 * rather than compiled in.
	 */
	class ConductorField
	{
		public:
			/// The default coincidence tolerance, RELATIVE to a filament's own
			/// radius. See coincides() for what it is and is not for.
			static constexpr double defaultCoincidenceTolerance = 1.0e-12;

			/// @param mu0In  the permeability. Must be finite and positive; a
			///               zero would make every conductor silently inert,
			///               which is worse than an error.
			///
			/// @throws std::invalid_argument if mu0In is not finite or not
			///         positive.
			explicit ConductorField( double mu0In = vacuumPermeability );

			/// Append a filament. It is copied; a CurrentFilament is three
			/// doubles.
			void add( CurrentFilament const &filament );

			/// Append a rectangle carrying a uniform current density. Copied,
			/// as meq::CoilSet::add() copies.
			void add( Coil const &coil );

			/// Every conductor, filaments and rectangles together. `empty()` is
			/// true only when there are none of either.
			std::size_t size() const;
			bool empty() const;

			std::size_t filamentCount() const;
			std::size_t coilCount() const;

			/// @throws std::out_of_range naming the index and the count.
			CurrentFilament const &filament( std::size_t index ) const;

			/// The filaments, in the order they were added.
			std::vector<CurrentFilament> const &filaments() const;

			/// The rectangles, as the set that evaluates them.
			CoilSet const &coils() const;

			/// The cross-section quadrature the RECTANGLES are integrated at;
			/// filaments have no quadrature and are unaffected. Forwarded to
			/// meq::CoilSet, whose default is what a solve uses -- raise it to
			/// build the reference field §7.2 of the plan asks CS-5 for.
			void setQuadratureOrder( int order );
			int quadratureOrder() const;

			/// The signed sum of EVERY conductor's current, filaments and
			/// rectangles alike, in amperes. The check a boundary integral of
			/// this field is made against, exactly as
			/// meq::CoilSet::totalCurrent() is.
			double totalCurrent() const;

			double mu0() const;

			/// psi_c at a point: the sum over the filaments of
			/// meq::filamentPsi() and over the rectangles of meq::coilPsi(),
			/// which is the superposition `Delta*`'s linearity licenses.
			///
			/// **Exactly zero on the axis**, bit for bit, because `k^2` is an
			/// exact factor of the kernel -- which is the boundary condition
			/// the free-boundary problem imposes at `r = 0` and is the reason
			/// the split does not disturb it.
			///
			/// **AND EVEN IN `r`, SO A POINT PAST THE AXIS IS ANSWERED RATHER
			/// THAN REFUSED.** `psi = r A_phi` and both factors change sign
			/// under `r -> -r`, so `psi( -r, z ) == psi( r, z )` exactly; this
			/// is the analytic continuation of the flux and not a clamp. Two
			/// of MEQ's evaluations extrapolate off the half-plane by design
			/// -- the exterior datum, whose transfer paths target a `Gamma`
			/// that MEETS the axis, and the critical-point Newton, which is
			/// allowed to leave its element -- and this is what lets them.
			/// gradPsi(), flux() and poloidalField() still refuse, because
			/// `d_r psi` is ODD where `psi` is even and a vector cannot take
			/// one rule for both entries; see ConductorField::psi()'s body.
			///
			/// @throws std::invalid_argument if the point is ON a filament,
			///         where psi is genuinely infinite. That refusal is
			///         meq::filamentPsi()'s and it is the LAST line of defence
			///         rather than the intended one -- see coincides(), which
			///         is how a caller finds out at setup instead of at the
			///         first evaluation.
			double psi( double r, double z ) const;

			/// grad_bar( psi_c ) of the whole set.
			/// @throws std::invalid_argument as psi() does.
			void gradPsi( double r, double z,
			              double &dPsiDr, double &dPsiDz ) const;

			/// q_c = ( 1/r ) grad_bar( psi_c ). NaN on the axis, as
			/// meq::filamentFlux() is.
			/// @throws std::invalid_argument as psi() does.
			void flux( double r, double z, double &qR, double &qZ ) const;

			/**
			 * B_pol of the conductors: `B_R = -q_z`, `B_Z = +q_r`, the same
			 * relabelling meq::poloidalField() applies to a solved flux.
			 *
			 * **AND IT IS FINITE ON THE AXIS WHERE flux() IS NaN, WHICH IS THE
			 * WHOLE REASON IT EXISTS.** `B_R` there is exactly zero and `B_Z`
			 * is the closed-form limit meq::filamentAxisFlux() and
			 * meq::coilAxisFlux() supply -- see those for the derivation and
			 * for why only one component is a limit. An output grid on a
			 * half-disc machine has its whole first column on `r = 0`, so a
			 * caller adding the conductors' field to a sampled one meets this
			 * at every such node rather than occasionally.
			 *
			 * This is the ONLY entry point here that special-cases the axis.
			 * psi() needs none -- it is exactly zero there -- and flux() keeps
			 * its NaN deliberately, so that a caller who has not thought about
			 * the axis is told rather than handed a plausible number.
			 *
			 * @throws std::invalid_argument as psi() does.
			 */
			void poloidalField( double r, double z,
			                    double &bR, double &bZ ) const;

			/**
			 * Does this point lie ON a filament?
			 *
			 * **THIS IS A COINCIDENCE TEST AND NOT A CLEARANCE TEST, AND THE
			 * DIFFERENCE IS MEASURED RATHER THAN ASSUMED.** A point NEAR a
			 * filament is not a problem: `Coils.hpp` records `filamentPsi()`
			 * agreeing with an independent transcription to 1.4e-12 and still
			 * correct at `eps = 1e-13` from the conductor -- sitting a constant
			 * 0.0397 above `-(1/2) ln eps` -- where the textbook form has been
			 * NaN since `1e-09`. And a large value at one evaluation point is
			 * not a large error, because nothing here approximates `psi_c` by a
			 * polynomial: `psi_c` is EVALUATED, and what is approximated is
			 * `psi_p`, which is smooth precisely BECAUSE `psi_c` carries the
			 * whole logarithm.
			 *
			 * So the only thing that fails is a point exactly on the ring, and
			 * the tolerance exists to catch "the mesher put a node here" -- a
			 * coordinate that came from the same double, possibly through a
			 * text round trip, which is `1e-16` relative at worst -- rather
			 * than to enforce any distance. **Do not widen it into a clearance
			 * rule**: nothing measured supports needing one, and a clearance is
			 * a much stronger claim about the discretisation than this plan
			 * makes.
			 *
			 * The distance is `hypot( r - radius, z - height )` in the
			 * half-plane and the tolerance is relative to the filament's own
			 * radius, which is strictly positive by CurrentFilament's own
			 * refusal.
			 *
			 * **RECTANGLES ARE NOT CONSIDERED AND THAT IS NOT AN OVERSIGHT.**
			 * A coil's field is a quadrature over its cross-section and is
			 * finite everywhere, inside it included, so a mesh point in a
			 * rectangle is an ordinary point -- and refusing one would reject
			 * the configuration this plan exists to make cheap.
			 */
			bool coincides( double r, double z ) const;

			/// The index of the FIRST filament this point coincides with, in
			/// insertion order, or `-1`. The return type is signed so that -1
			/// can be the sentinel; a non-negative answer may be cast to
			/// std::size_t and passed to filament(). Two filaments at one
			/// location is a configuration error rather than an ambiguity, so
			/// the first is the only one a caller needs to be told about.
			int indexAt( double r, double z ) const;

			double coincidenceTolerance() const;

			/// @throws std::invalid_argument if the tolerance is not finite or
			///         is negative. Zero is allowed and means exact equality,
			///         which is a defensible choice for a caller that knows
			///         its coordinates never round-trip.
			void setCoincidenceTolerance( double toleranceIn );

			/**
			 * How far the FARTHEST conductor lies inside a semicircular `Gamma`
			 * of radius @a rhoGamma centred at `( 0, centreZ )`, in metres.
			 * Positive when every conductor is strictly inside it; an empty
			 * field gives infinity.
			 *
			 * **THIS IS THE MIRROR IMAGE OF meq::ExteriorCoilSet::clearance()
			 * AND THE TWO PRECONDITIONS ARE OPPOSITE ONES.** A conductor handed
			 * to setExteriorConductors() must be OUTSIDE `Gamma`, so that
			 * `psi_coil` is `Delta*`-harmonic in `Omega` and the conductor can
			 * enter through the boundary alone. A conductor handed to
			 * setConductorField() must be INSIDE it, and for the dual reason:
			 * the exterior field is represented by a Gegenbauer series in
			 * `rho^( 1 - n )`, which converges outside `Gamma` only when every
			 * source it stands for is within it. A subtracted conductor beyond
			 * `Gamma` puts a singularity in the region that series describes,
			 * and the run converges to a machine nobody described.
			 *
			 * So a conductor is required to be strictly on one side or the
			 * other depending on WHICH ROUTE carries it, and one straddling
			 * `Gamma` belongs to neither -- which is exactly what
			 * `ExteriorCoilSet::clearance()`'s own documentation says from its
			 * side.
			 *
			 * **AND IT CONSTRAINS NOTHING WITHOUT AN EXTERIOR COUPLING.** On a
			 * fixed-boundary problem there is no series and no `Gamma`, so a
			 * subtracted conductor may sit anywhere at all -- inside the mesh,
			 * outside it, or straddling its edge. The caller is
			 * GradShafranovSolver, which knows whether a DtN is installed; this
			 * function only measures.
			 *
			 * The distance is to the FARTHEST point of each member -- the most
			 * distant corner of a rectangle, the ring itself for a filament --
			 * so a conductor straddling `Gamma` reports a negative containment
			 * rather than being judged by its centre.
			 *
			 * @throws std::invalid_argument on a non-finite argument or a
			 *         non-positive radius.
			 */
			double containment( double centreZ, double rhoGamma ) const;

		private:
			double mu0Value;
			double toleranceValue;
			std::vector<CurrentFilament> filamentList;
			CoilSet coilList;
	};
}

#endif
