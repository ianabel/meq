#ifndef MEQ_DEVICEQUADRATURE_HPP
#define MEQ_DEVICEQUADRATURE_HPP

#include "SourceView.hpp"

#include "mfem.hpp"

/**
 * @file
 *
 * THE PLASMA-CURRENT QUADRATURE AS A DEVICE ELEMENT LOOP.
 *
 * `GradShafranovSolver::assemblePlasmaCurrent()` integrates
 * `scaledF( R, z, psi )/R` over the plasma support, once per constraint
 * evaluation and 28 times in a twelve-step DIII-D run. It is the third item of
 * DEVICE-PLAN.md section 3.4a, and that section originally called it a
 * reduction. **It is not**: its output is a scalar but its body is an element
 * loop that builds shape functions, sums `psi = shape . dofs` per quadrature
 * point and evaluates a VIRTUAL on the polymorphic source. *A scalar output
 * does not make a reduction; the body does.*
 *
 * SO THIS IS ITEM 1's SHAPE ON MEQ'S OWN SIDE OF THE FENCE. Upstream owns the
 * batched local residual inside `DarcyHybridization`; this loop is MEQ's, over
 * MEQ's own space, and needs nothing from them. What it needed was the POD --
 * meq::SourceView, and specifically meq::makePlasmaSourceView(), since
 * `scaledF()` is the plasma term alone and every production case is
 * coil-augmented (M-191).
 *
 * **IT IS BUILT OUT OF MFEM'S ABSTRACTIONS AND NOT OUT OF A HAND-WRITTEN
 * KERNEL, WHICH IS THE STANDING PREFERENCE AND ALSO THE PORTABLE CHOICE.**
 * `mfem::FiniteElementSpace::GetElementRestriction()` gathers the L-vector into
 * an E-vector on the device; `mfem::Mesh::GetGeometricFactors()` supplies the
 * physical coordinates and the Jacobian determinants at every quadrature point;
 * `mfem::FiniteElement::GetDofToQuad()` supplies the shape functions there; and
 * `mfem::forall` runs the loop. MEQ writes no CUDA and no HIP, so the same
 * source serves both -- see CLAUDE.md, *Prefer a maintained library*.
 *
 * WHAT IT REFUSES, AND EACH REFUSAL IS A CASE THE HOST LOOP HANDLES AND THIS
 * ONE CANNOT YET. A refusal returns false and the caller keeps the host path,
 * which is meq::makeSourceView()'s contract and for the same reason: there is
 * no approximate mode worth having in a constraint functional.
 *
 *   * **A quadrature rule that is not the same on every element.** The host
 *     builds its rule per element from `2k + tr.OrderW() + extra`, and
 *     `GetGeometricFactors()` takes ONE rule for the mesh. On a straight-sided
 *     mesh `OrderW()` is 0 everywhere and the rules agree; on a curved one they
 *     need not. Checked rather than assumed.
 *   * **A vector or non-L2 space.** The potential space is scalar L2 and the
 *     restriction, the shape layout and the dof gather all assume it.
 *   * **A conductor cache whose offsets do not match the rule.** `psi_c` is
 *     carried now -- see below -- as meq::SourceIntegrator holds it, a flat
 *     value array indexed `values[ offsets[ e ] + q ]`. What is refused is a
 *     cache built against a DIFFERENT rule, detected by the offsets not being
 *     `e*nq`. Refused rather than silently dropped, because dropping `psi_c`
 *     would put the coil flux out of `psi` and change the equilibrium, which
 *     is the class of silent error COIL-SUBTRACTION-PLAN.md exists to avoid.
 *
 * **AND THE CONDUCTOR SHIFT IS WHY UPSTREAM'S SEAM PASSES `(element, point)`
 * RATHER THAN A FLAG ARRAY.** MEQ proposed an `int const *elementFlags`
 * parameter for the plasma-connectivity test; upstream passed the two loop
 * counters instead and gave this cache as the reason, *"a per-point cached
 * field ... which an elementFlags parameter would not have covered"*. They
 * were right, and this loop indexes exactly that way, so it transfers to their
 * `f( element, point, x, y, u )` without reshaping.
 *
 * AND IT RETURNS PER-ELEMENT PARTIALS RATHER THAN A TOTAL, DELIBERATELY.
 * `assemblePlasmaCurrent()` carries a long comment refusing `reduction( + : )`
 * because the association would then depend on the thread count and take
 * `psi_ax` and every published digit with it; it sums an array indexed by
 * element, IN ELEMENT ORDER, for that reason. This returns the same array so
 * the caller can sum it the same way. **An `mfem::reduce` with `SumReducer`
 * would be the natural primitive for a total and is the wrong one here**, for
 * that documented reason and no other.
 */

namespace meq
{

	/**
	 * `int_K scaledF( R, z, psi )/R` per element, on the device.
	 *
	 * @param fes        the scalar L2 potential space.
	 * @param potential  the potential block of the state, as an L-vector.
	 * @param view       the PLASMA source as POD -- meq::makePlasmaSourceView().
	 * @param knots      the flat knot array @a view indexes into.
	 * @param inPlasma   one entry per element, non-zero where XP-1's fill
	 *                   reached. The pointwise `insidePlasma( psi )` test is
	 *                   inside `sourceF` and is independent of this.
	 * @param conductorValues `psi_c` per quadrature point, as
	 *                   meq::SourceIntegrator::conductorCacheValues() holds
	 *                   it. **Empty is not a special case**: the shift is then
	 *                   exactly zero and every path without the
	 *                   COIL-SUBTRACTION-PLAN.md split is bit-unchanged, which
	 *                   is the same contract `conductorShiftAt()` has.
	 * @param conductorOffsets the matching
	 *                   meq::SourceIntegrator::conductorCacheOffsets().
	 * @param extraOrder `sourceQuadratureExtra`, the solver's own.
	 * @param perElement resized to the element count and left **DEVICE-VALID**;
	 *                   the caller funnels. Bringing it over is `NE` doubles --
	 *                   39 KB on the shipped machine case against the state's
	 *                   698 KB -- and whether that is wanted is the caller's
	 *                   question, not this kernel's.
	 *
	 * @return false, having touched nothing, on any refusal above or with no
	 *         mfem::Device configured.
	 */
	bool plasmaCurrentPerElementOnDevice( mfem::FiniteElementSpace &fes,
	                                      mfem::Vector const &potential,
	                                      SourceView const &view,
	                                      std::vector<KnotView> const &knots,
	                                      mfem::Array<int> const &inPlasma,
	                                      std::vector<double> const &conductorValues,
	                                      std::vector<int> const &conductorOffsets,
	                                      int extraOrder,
	                                      mfem::Vector &perElement );

	/// The launch, in the one translation unit a device compiler owns. See
	/// DeviceRows.hpp for why that split exists.
	bool plasmaCurrentKernel( mfem::Vector const &elementPotential,
	                          mfem::Vector const &quadWeights,
	                          mfem::Vector const &detJ,
	                          mfem::Vector const &coordinates,
	                          mfem::Array<mfem::real_t> const &shape,
	                          mfem::Array<int> const &inPlasma,
	                          mfem::Vector const &conductorShift,
	                          bool hasConductorShift,
	                          SourceView const &view,
	                          mfem::Array<KnotView> const &knots,
	                          int elements, int quadPoints, int dofs,
	                          int spaceDimension,
	                          mfem::Vector &perElement );

}

#endif // MEQ_DEVICEQUADRATURE_HPP
