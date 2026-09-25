#ifndef MEQ_DEVICENPC_HPP
#define MEQ_DEVICENPC_HPP

#include "SourceView.hpp"

#include "mfem.hpp"

#include <memory>
#include <vector>

/**
 * @file
 *
 * MEQ'S NPC OPERATOR, WITH THE ELEMENT-LOCAL NONLINEARITY SUPPLIED AS A DEVICE
 * CALLABLE.
 *
 * **THIS IS THE CALL SITE UPSTREAM'S SEAM IS FOR, AND MEQ HAD TO HOLD IT.**
 * `DarcyHybridization::NPCResidual( b, x, x_tr, r, r_tr, f )` and
 * `NPCGradient( x, x_tr, g )` take the caller's integrand and fall back to the
 * element loop when `CanBatchNonlinearResidual()` declines. MEQ's Newton holds
 * an `mfem::Operator &`, not a hybridization, and `mfem::DarcyNPCOperator` has
 * no way to take an integrand -- so MEQ subclasses it and overrides the two
 * methods that would otherwise call the untemplated pair.
 *
 * **MEQ FILED AN OBJECTION TO DOING IT THIS WAY AND THE OBJECTION WAS WRONG.**
 * It argued that a consumer-side operator would restate *"the prolongation, the
 * state, the skeleton load and the essential-row convention"*. It does not:
 * those live in `NPCResidualStart()` and `NPCResidualFinish()`, INSIDE
 * `NPCResidual()`. What the operator restates is block splitting and two
 * synchronisation hops, about twenty lines, and the hops are the part to copy
 * exactly rather than to reason about -- `DarcyNPCOperator::Mult`'s own comment
 * records that dropping the first turns M-79's symptom back on.
 *
 * WHAT THE INTEGRAND IS, AND IT IS ONE EXPRESSION RATHER THAN TWO BRANCHES.
 * `meq::SourceIntegrator` accumulates
 * `-weight * sourceValue( R, z, psi + psi_c, e )/R` against the basis, and
 * `sourceValue` is `elementCarriesPlasma( e ) ? source->f( ... )
 *                                             : fOutsidePlasma( R, z )`.
 * For a coil-augmented source those are `plasma + coil` and `coil`, so the POD
 * form is **( plasma, gated by the element flag ) + coil** and the gate is the
 * only difference between the branches. `SourceViewTests` pins that identity to
 * the bit, including across a coil overlap and on a closed rectangle edge.
 *
 * The weight and the contraction with the basis are upstream's, so `f` returns
 * the integrand and nothing else.
 *
 * **THE VIEW IS REBUILT PER CALL AND THAT IS REQUIRED RATHER THAN TIDY.**
 * `psi_ax`, `psi_bnd` and the current scale are unknowns of the bordered Newton
 * and move between residual evaluations within one step. Upstream enforces it
 * from their side -- `f` is taken BY VALUE and captured by value, so a fresh
 * lambda per evaluation is the ordinary usage. What does NOT move within a
 * solve is the knots, the coils and the plasma fill (XP-2 freezes the support),
 * so those are uploaded once when the operator is built and the operator is
 * rebuilt whenever `prepare()` re-finalises the hybridization.
 */

namespace meq
{

	class Source;

	/// What the integrand needs besides the live source, all of it fixed for
	/// the operator's lifetime. Pointers are BORROWED and must outlive it.
	struct DeviceIntegrandInputs
	{
		/// The live source, re-read per residual for the moving scalars.
		Source const *source = nullptr;

		/// XP-1's per-element connectivity flag, one per element, non-zero
		/// where the fill reached. Null means every element carries plasma,
		/// which is meq::SourceIntegrator::elementCarriesPlasma()'s own
		/// constant-true case and keeps an unconfined solve bit-identical.
		int const *elementInPlasma = nullptr;
		int elementCount = 0;

		/// meq::SourceIntegrator's conductor cache, `values[ offsets[ e ] + q ]`.
		/// Empty is the ordinary no-split case and shifts by exactly zero.
		double const *conductorValues = nullptr;
		int const *conductorOffsets = nullptr;
		int conductorValueCount = 0;
	};

	/**
	 * An `mfem::DarcyNPCOperator` that routes through the caller-integrand
	 * seam, or **null** when this configuration cannot use it.
	 *
	 * @return null when no device is configured, when the POD cannot be built
	 *         for @a inputs.source -- `meq::makePlasmaSourceView()`'s refusal,
	 *         which is a rotating source or an unrepresentable profile -- or
	 *         when MEQ was not built with a device language. **The caller then
	 *         constructs the stock operator and nothing else changes.** That is
	 *         the same convention every batched route here follows, and it
	 *         means the fallback is a construction-time decision rather than a
	 *         per-call branch.
	 *
	 * Even when this returns an operator, `CanBatchNonlinearResidual()` may
	 * still decline per call -- a mixed mesh, a rule that is not explicit --
	 * and upstream's templated entry point then runs the element loop. So there
	 * are two independent fallbacks and neither is a property of the build.
	 */
	std::unique_ptr<mfem::DarcyNPCOperator> makeDeviceNPCOperator(
		mfem::DarcyHybridization &hybridization,
		mfem::Array<int> const &offsets,
		mfem::BlockVector const &load,
		DeviceIntegrandInputs const &inputs );

}

#endif // MEQ_DEVICENPC_HPP
