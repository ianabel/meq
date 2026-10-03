#ifndef MEQ_RESIDENTITERATE_HPP
#define MEQ_RESIDENTITERATE_HPP

#include "DeviceNPC.hpp"
#include "SourceView.hpp"

#include "mfem.hpp"

#include <vector>

/**
 * @file
 *
 * THE LAST HOST READERS OF THE BORDERED ITERATE, GIVEN DEVICE ROUTES.
 *
 * **WHAT THIS IS FOR, AND IT IS A COUNT RATHER THAN A LEG.**
 * [M-232](../MEASUREMENTS.md) forced the axis and X-point routes open on
 * `machine-f-diiid` under `--device cuda` and watched the same 12.5 MB move
 * from the axis constraint to `refreshXPoint()` to `limiterValue()`: *a leg
 * ported in isolation does not remove a transfer, it RENAMES one*, and the
 * count falls only when the last host reader of the iterate goes. Watching the
 * iterate's own buffer through every funnel of that run named exactly ten:
 *
 *   * the five plasma element loops -- `assemblePlasmaCurrent()`,
 *     `assembleCurrentColumn()`, `assembleCurrentRow()`,
 *     `assembleCurrentNormalisationCorner()` and
 *     `assembleNormalisationColumn()` -- which are meq::DevicePlasmaLoops;
 *   * the transmission constraint's gather, which is meq::DeviceRows and was
 *     already built;
 *   * the `psi_bnd` row, `limiterValue()`, which is orderedDotOnDevice();
 *   * and three COPY CONSTRUCTIONS of the iterate, whose funnels were the
 *     read. A copy of a device-valid vector is born device-valid
 *     ([M-205](../MEASUREMENTS.md)), so those need no kernel at all -- only
 *     the funnel left out and every later reader of the copy device-aware,
 *     which is copyEntriesOnDevice() for the one that is masked afterwards.
 *
 * **NONE OF IT IS A SEARCH**, which is what `DEVICE-PLAN.md` said and what
 * makes this file short beside AxisSearch.hpp: every one is a reduction, a
 * scatter or a gather over indices the host already holds.
 *
 * **EVERY ROUTE IS GATED BY THE SOLVER ON THE SAME PREDICATE AS THE AXIS
 * ROUTE**, `DeviceIsValid() && !HostIsValid()` of the vector it is handed,
 * with `MEQ_AXIS_DEVICE=force` overriding it. That gate is a predicate about
 * the DATA and is shut on every shipped case today, so every host number is
 * unchanged by construction; what the routes buy shows only once something
 * keeps the iterate on the device -- and the force switch is how that is
 * measured before it does.
 *
 * WHAT THE PLASMA LOOPS REPRODUCE, AND HOW.
 *
 *   **THE HOST'S OWN RADII AND WEIGHTS, NOT A RECOMPUTATION.** The host loops
 *   take `scratch.Transform( ip )` for the radius and `ip.weight*Weight()`
 *   for the weight at every quadrature point. Both are constants of the mesh,
 *   so they are evaluated ONCE, on the host, by those same calls, and uploaded
 *   -- phase B. `mfem::Mesh::GetGeometricFactors()` would have computed them
 *   on the device, by a different expression, AND it calls `EnsureNodes()` on
 *   the solver's own mesh, which `CLAUDE.md` records costing `GridSampler` two
 *   orders of magnitude through a mesh acquiring nodes without becoming
 *   curved. Neither is acceptable for a route whose acceptance is the host's
 *   digits.
 *
 *   **EVERY EXPRESSION IN THE HOST'S ORDER**, and the device file is compiled
 *   without FMA contraction for meq::DeviceRows' recorded reason: the host
 *   build carries no `-march`, so it rounds a multiply and an add twice.
 */

namespace meq
{

	/// Which of the five plasma element loops a launch is.
	enum class PlasmaQuantity : int
	{
		/// Per element: `sum_q ( w W ) scaledF / R`. `assemblePlasmaCurrent()`.
		Current = 0,
		/// Per element, TWO sums: `sum_q ( w W / R ) dF/dpsi_ax` and the same
		/// against `psi_bnd`. `assembleCurrentNormalisationCorner()`.
		NormalisationCorner = 1,
		/// Scattered: `-( w W )( scaledF/scale )/R phi_j`. `assembleCurrentColumn()`.
		CurrentColumn = 2,
		/// Scattered: `( w W ) scaledDFdPsi/R phi_j`. `assembleCurrentRow()`.
		CurrentRow = 3,
		/// Scattered: `-( w W ) dF/dpsi_ax/R phi_j`. `assembleNormalisationColumn()`.
		AxisColumn = 4,
		/// @copydoc AxisColumn, against `psi_bnd`.
		BoundaryColumn = 5
	};

	/**
	 * Everything one element's thread needs, as plain types and pointers.
	 *
	 * No default member initialisers, for the reason meq::AxisCandidate gives.
	 */
	struct PlasmaLoopView
	{
		SourceView source;
		KnotView const *knots;
		/// XP-1's per-element flag; zero skips the element whole.
		int const *flags;
		/// `psi_c` per quadrature point, `shift[ shiftOffsets[ e ] + q ]`.
		int const *shiftOffsets;
		mfem::real_t const *shift;
		int hasShift;
		/// `shape[ q*dofs + j ]`: the reference element's shape functions at
		/// the rule's points, exactly `CalcShape()`'s.
		mfem::real_t const *shape;
		/// `radius[ e*quadPoints + q ]`, likewise `z` and `weight`.
		mfem::real_t const *radius;
		mfem::real_t const *z;
		mfem::real_t const *weight;
		/// `elementDofs[ e*dofs + j ]`, the potential-space dof of local `j`.
		int const *elementDofs;
		mfem::real_t const *state;
		int potentialStart;
		int elements;
		int quadPoints;
		int dofs;
		/// The current scale, which CurrentColumn divides by.
		mfem::real_t scale;
		int quantity;
	};

	/**
	 * The five plasma element loops of the bordered Newton, on the device.
	 *
	 * **PHASE B IS prepare() AND IT RUNS ONCE PER MESH.** It refuses -- and
	 * the caller keeps the host loop -- on anything the host loop handles and
	 * these cannot: a vector or non-L2 potential space (the scatter is
	 * race-free only because no two elements share a dof), or elements that do
	 * not share ONE finite element and ONE rule, since the host builds its rule
	 * per element from `2k + OrderW() + extra` and a curved mesh can disagree
	 * with itself.
	 *
	 * It is keyed on the space, the mesh and the mesh's sequence number, so an
	 * adaptive cycle that refines the mesh under a surviving solver rebuilds it
	 * rather than indexing a stale table -- which would not throw, it would
	 * name other elements.
	 */
	class DevicePlasmaLoops
	{
		public:
			/// Build or confirm phase B for @a potential at @a extraOrder.
			/// @return whether the routes are available.
			bool prepare( mfem::FiniteElementSpace &potential, int extraOrder );

			/**
			 * One or two per-element partial sums, left DEVICE-VALID.
			 *
			 * Per element and not a total, for `assemblePlasmaCurrent()`'s
			 * reason: the host sums an array indexed by element, in element
			 * order, so that the answer does not depend on a thread count --
			 * and the caller does the same with these.
			 *
			 * @param second required for NormalisationCorner and ignored
			 *        otherwise.
			 */
			bool integrate( PlasmaQuantity quantity, SourceView const &source,
			                DeviceConstants const &constants,
			                mfem::Vector const &state, int potentialStart,
			                mfem::Vector &first, mfem::Vector *second ) const;

			/**
			 * A column or a row the size of @a state, zero outside the potential
			 * dofs of plasma elements, left DEVICE-VALID.
			 */
			bool scatter( PlasmaQuantity quantity, SourceView const &source,
			              DeviceConstants const &constants,
			              mfem::Vector const &state, int potentialStart,
			              mfem::real_t scale, mfem::Vector &out ) const;

		private:
			PlasmaLoopView view( PlasmaQuantity quantity, SourceView const &source,
			                     DeviceConstants const &constants,
			                     mfem::Vector const &state, int potentialStart,
			                     mfem::real_t scale ) const;

			mfem::FiniteElementSpace const *space = nullptr;
			mfem::Mesh const *mesh = nullptr;
			long sequence = -1;
			int extra = 0;
			/// 0 not asked for this signature, 1 available, -1 refused.
			int status = 0;

			int elements = 0;
			int quadPoints = 0;
			int dofs = 0;
			mfem::Vector shape;
			mfem::Vector radius;
			mfem::Vector z;
			mfem::Vector weight;
			mfem::Array<int> elementDofs;
	};

	/**
	 * `out[ k ] = source[ indices[ k ] ]`, on the device; @a out is left
	 * device-valid. A support's worth of a device-valid vector, so that
	 * meq::CompressedRows can compress a row without the row coming over.
	 */
	bool gatherOnDevice( mfem::Vector const &source,
	                     mfem::Array<int> const &indices, mfem::Vector &out );

	/**
	 * `out[ 0 ] = sum_k weights[ k ] source[ indices[ k ] ]`, summed in @a k
	 * order from zero by ONE thread, which is the host's loop exactly. Used for
	 * the `psi_bnd` row, a handful of entries, where a parallel reduction would
	 * buy nothing and cost the order.
	 */
	bool orderedDotOnDevice( mfem::Vector const &source,
	                         mfem::Array<int> const &indices,
	                         mfem::Vector const &weights, mfem::Vector &out );

	/**
	 * `destination[ indices[ k ] ] = source[ indices[ k ] ]`, on the device.
	 * The cold state's mask, which keeps the essential trace entries of the
	 * iterate and nothing else.
	 */
	bool copyEntriesOnDevice( mfem::Vector const &source,
	                          mfem::Array<int> const &indices,
	                          mfem::Vector &destination );

	/**
	 * The kernel halves, in the device translation unit; no-ops returning
	 * false without a device language. See AxisSearch.hpp for why the split.
	 */
	bool plasmaElementKernel( PlasmaLoopView const &inputs,
	                          mfem::real_t *first, mfem::real_t *second );
	/// @copydoc plasmaElementKernel
	bool plasmaScatterKernel( PlasmaLoopView const &inputs, mfem::real_t *out );
	/// @copydoc plasmaElementKernel
	bool gatherKernel( mfem::Vector const &source, mfem::Array<int> const &indices,
	                   mfem::Vector &out );
	/// @copydoc plasmaElementKernel
	bool orderedDotKernel( mfem::Vector const &source,
	                       mfem::Array<int> const &indices,
	                       mfem::Vector const &weights, mfem::Vector &out );
	/// @copydoc plasmaElementKernel
	bool copyEntriesKernel( mfem::Vector const &source,
	                        mfem::Array<int> const &indices,
	                        mfem::Vector &destination );

}

#endif // MEQ_RESIDENTITERATE_HPP
