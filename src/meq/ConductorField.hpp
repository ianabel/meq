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
	 * **FILAMENTS ONLY, AND THAT IS CS-1 RATHER THAN A LIMITATION.** A filament
	 * is the one conductor MEQ structurally CANNOT carry any other way: a point
	 * source has no finite-element representation as a current density, `psi`
	 * near it is logarithmic, and `meq::CoilSet::f()` has nothing to add for it
	 * because the current density is infinite on a set of measure zero. So
	 * subtraction is not one route to filament support, it is the only one --
	 * and the remainder is well behaved for exactly the reason the subtraction
	 * works: `psi_c` carries the whole logarithm, so `psi_p` sees a bounded
	 * right-hand side supported on the plasma alone. A rectangle is the same
	 * code with a quadrature rule around it and is CS-1b.
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

			std::size_t size() const;
			bool empty() const;

			/// @throws std::out_of_range naming the index and the size.
			CurrentFilament const &filament( std::size_t index ) const;

			/// The filaments, in the order they were added.
			std::vector<CurrentFilament> const &filaments() const;

			/// The signed sum of the filament currents, in amperes. The check
			/// a boundary integral of this field is made against, exactly as
			/// meq::CoilSet::totalCurrent() is.
			double totalCurrent() const;

			double mu0() const;

			/// psi_c at a point: the sum over filaments of meq::filamentPsi().
			///
			/// **Exactly zero on the axis**, bit for bit, because `k^2` is an
			/// exact factor of the kernel -- which is the boundary condition
			/// the free-boundary problem imposes at `r = 0` and is the reason
			/// the split does not disturb it.
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

		private:
			double mu0Value;
			double toleranceValue;
			std::vector<CurrentFilament> filamentList;
	};
}

#endif
