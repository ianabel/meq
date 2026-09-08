#ifndef MEQ_FLUXEXTRACTION_HPP
#define MEQ_FLUXEXTRACTION_HPP

#include <cstddef>
#include <functional>

#include "CriticalPoints.hpp"
#include "FluxFamily.hpp"
#include "FluxSurfaces.hpp"

/*
 * Turning a solved psi_h and q_h into a meq::FluxSurfaceFamily: the MFEM half of
 * INVERSION-PLAN.md stage IN-6.
 *
 * The container, the flux label, the interpolation and the per-psi cache are in
 * FluxFamily.hpp and are MFEM-free so that continuous integration can gate the
 * contract. This file is the part that needs a mesh: trace, fit, average, once
 * per level.
 *
 * IT ADDS NOTHING TO THE PHYSICS AND THAT IS DELIBERATE. Every quantity here is
 * IN-0's tracer, IN-1's pointwise metric and IN-2's averages, called in a loop
 * and copied into a value. The two things it computes for itself are the
 * enclosed volume and the enclosed cross-section area, and they are computed
 * from the SAME pointwise rho' the metric uses -- see below.
 *
 * A LEVEL THAT CANNOT BE TRACED KILLS THE WHOLE FAMILY, WHICH IS THE POINT.
 * extractFluxSurfaces() throws on the first level it cannot reach rather than
 * returning a family with a hole in it. A family with a hole is worse than no
 * family: FluxSurfaceFamily::at() interpolates between the surfaces it has, so
 * a missing surface is silently bridged by its neighbours and the consumer is
 * handed a plausible number for a place meq could not look at. That is
 * MANTA-COUPLING.md section 8's "fail by throwing, never by fudging", and the
 * message names the level and its Psi_N so that the cut can be moved rather
 * than guessed at.
 *
 * THE VOLUME IS A SECOND, INDEPENDENT ROUTE TO V' AND IT IS WHY IT IS HERE.
 * Green's theorem gives the enclosed volume as a contour integral with no
 * gradient in it at all,
 *
 *     V = closed-integral pi R^2 dz,        A = closed-integral R dz
 *
 * while V' = closed-integral 2 pi R dl / | grad psi | is a weighted line
 * integral that divides by the flux everywhere. The coarea formula says
 * V'( psi ) = -dV/dpsi, so the two must agree -- with nothing in common but the
 * node positions, and needing no reference value. That is the third leg
 * SurfaceAverage.hpp says it is short of, on one quantity rather than on all of
 * them, and tests/convergence/FluxGridConvergence.cpp asserts it.
 *
 * dz/dtheta IS TAKEN POINTWISE FROM rho' AND NEVER BY DIFFERENCING THE z_j.
 * With z = z_ax + rho( theta ) sin( theta ),
 *
 *     dz/dtheta = rho'( theta ) sin( theta ) + rho( theta ) cos( theta )
 *
 * and rho' comes from the solved flux by IN-1's identity. Differencing
 * neighbouring z_j instead is INVERSION-PLAN.md section 3.2's metric trap, which
 * IN-1 measured at 7.03 against 1.97 on an arc length and IN-2 measured again on
 * an average. This is the same trap on a third integrand, and the differenced
 * column is kept live as a control in the test rather than argued away.
 */

namespace meq
{

	/// How the family's levels are laid out between the cuts.
	enum class FluxLevelSpacing
	{
		/// Equispaced in rho = sqrt( Psi_N ). THE DEFAULT, because rho is what
		/// the geometry is smooth in -- FluxFamily.hpp section 1 -- so equal
		/// steps in it are equal steps in the thing being interpolated.
		Radial,

		/// Equispaced in Psi_N. Kept because a consumer whose own grid is in
		/// normalised flux wants its surfaces where its nodes are, and because
		/// it is the control that says the choice of layout was measured.
		NormalisedFlux
	};

	/// What to extract.
	struct FluxFamilyOptions
	{
		/// Surfaces in the family, and nodes on each.
		std::size_t surfaces = 24;
		std::size_t angles = 128;

		/// The cut, in Psi_N. See FluxFamily.hpp section 2 for the decision and
		/// CLAUDE.md's IN-6 section for the measurement behind the defaults.
		double innerCut = 0.05;
		double outerCut = 0.95;

		FluxLevelSpacing spacing = FluxLevelSpacing::Radial;

		/// g( psi ) = R B_toroidal, for RoPP (142)'s safety factor. Empty by
		/// default, in which case FluxSurface::safetyFactor is zero and
		/// FluxSurfaceFamily::safetyFactorAvailable is false -- a meq::Source
		/// carries g g' and not g, so this is the caller's to supply from its
		/// own profile or not at all. Absent and zero are told apart by the
		/// flag rather than by the value.
		std::function<double( double psi )> toroidalField;
	};

	/**
	 * Extract the family.
	 *
	 * @param tracer      already pointed at whichever potential the caller wants.
	 *                    Potential::PostProcessed is the tracer's own default and
	 *                    is the right one: IN-0 measured the traced curve 60x,
	 *                    54x and 83x closer to the truth at k = 1, 2, 3 on the
	 *                    same mesh, at FEWER corrector iterations per point.
	 * @param axis        IN-A's magnetic axis, a root of q_h. Its psi is the
	 *                    psi_ax the label is normalised by.
	 * @param psiBoundary psi on the plasma boundary -- zero on every
	 *                    fixed-boundary run in this tree, and
	 *                    GradShafranovSolver::psiBoundary() once FB-3's limiter
	 *                    border makes it an unknown.
	 *
	 * @throws std::invalid_argument if the options do not describe a family: a
	 *         cut outside ( 0, 1 ), an inner cut at or past the outer one, fewer
	 *         than two surfaces, fewer than three angles.
	 * @throws std::runtime_error if psi_ax equals psi_bnd, or if any level
	 *         cannot be traced, fitted or closed -- see the header on why one
	 *         bad level takes the family with it.
	 */
	FluxSurfaceFamily extractFluxSurfaces( ContourTracer const &tracer,
	                                       CriticalPoint const &axis,
	                                       double psiBoundary,
	                                       FluxFamilyOptions const &options );

	/**
	 * The same from a solved solver, which is what a driver has in hand.
	 *
	 * Builds the tracer, locates the axis with meq::CriticalPointFinder and
	 * takes psi_bnd from the solver. Nothing more, and it is here so that the
	 * four lines are written once.
	 *
	 * @throws std::invalid_argument if @a which is Potential::PostProcessed and
	 *         GradShafranovSolver::postProcess() has not been called, which is
	 *         ContourTracer's own refusal.
	 * @throws whatever findAxis() throws when there is no unambiguous magnetic
	 *         axis -- which is a state meq cannot express a flux label against
	 *         at all, and is exactly the refusal MANTA-COUPLING.md section 8
	 *         asks for rather than a plausible answer.
	 */
	FluxSurfaceFamily extractFluxSurfaces( GradShafranovSolver const &solver,
	                                       FluxFamilyOptions const &options,
	                                       Potential which
	                                           = Potential::PostProcessed );

	/**
	 * The volume enclosed by one fitted surface, by Green's theorem on its own
	 * nodes: V = closed-integral pi R^2 dz, with dz/dtheta pointwise from rho'.
	 *
	 * meq::FluxSurface::volume is this, and it is exposed separately because
	 * the coarea check needs the volume of surfaces that are NOT in the family
	 * -- the four the Richardson difference straddles the level with. Calling
	 * the shipping routine rather than re-implementing it in a test is what
	 * makes that check a check of what ships.
	 *
	 * @throws std::invalid_argument if the fit carries fewer than three nodes.
	 */
	double enclosedVolume( AngleParametrisation const &fit );

	/**
	 * The enclosed volume of one fitted surface with dz/dtheta obtained by
	 * CENTRAL DIFFERENCING the neighbouring rho_j instead of pointwise from the
	 * solved flux.
	 *
	 * THE CONTROL, AND IT IS NOT TEST SCAFFOLDING. INVERSION-PLAN.md section
	 * 3.2's metric trap is that a spectrally accurate rule fed a second-order
	 * Jacobian is a second-order scheme with nothing in its output to say so.
	 * IN-1 measured it on an arc length -- 7.03 against 1.97 -- and IN-2 on an
	 * average; this is the same trap on a third integrand, and it lives beside
	 * the answer for the same reason SurfaceAverages::averageDifferenced() does,
	 * so that the claim stays a live column rather than a remark.
	 */
	double enclosedVolumeByDifferencedSlope( AngleParametrisation const &fit );

}

#endif // MEQ_FLUXEXTRACTION_HPP
