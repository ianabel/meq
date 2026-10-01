#ifndef MEQ_DEVICENPC_HPP
#define MEQ_DEVICENPC_HPP

#include "Crossings.hpp"
#include "SourceView.hpp"

#include "mfem.hpp"

#include <cstring>
#include <memory>
#include <type_traits>
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
	 * **PHASE B'S PAYLOAD, OWNED BY THE CALLER SO THAT IT OUTLIVES THE
	 * OPERATOR THAT READS IT.**
	 *
	 * DEVICE-PLAN.md section 3.3 draws phase B as a prologue -- *configure the
	 * device problem*, one upload, then run -- and MEQ was running it **sixteen
	 * times inside phase C** -> [M-226](../MEASUREMENTS.md). The five buffers
	 * below are per-mesh and per-source constants; the operator that reads them
	 * is rebuilt whenever the hybridization is re-finalised, which is once per
	 * support sweep and again after every right-hand-side re-assembly, and
	 * every rebuild re-uploaded all five.
	 *
	 * **THE MEASUREMENT THAT FOUND IT COULD NOT HAVE BEEN A CLOCK, AND WAS NOT
	 * A DOWNLOAD COUNT EITHER.** 1.22 MB of h2d over a run whose whole
	 * download traffic is 0.888 MB, in a phase the published acceptance reports
	 * as `0.2 transfers a step` -- a figure that counts the other direction and
	 * is in any case three transfers in phase A. *An instrument that can see
	 * one direction of the bus cannot measure a plan whose first phase is an
	 * upload.*
	 *
	 * **WHAT IT IS NOT IS A CACHE OF THE STATE.** Nothing here is written by
	 * phase C: `knots` and `coils` are the source's own description, `flags` is
	 * XP-1's per-element connectivity -- frozen within a solve by XP-2 and
	 * moving only between support sweeps -- and `shift`/`offsetsOfShift` are
	 * the conductor cache, which is per mesh. A buffer phase C wrote would have
	 * to come back down, and none of these does, which is exactly what makes
	 * them phase B's and not phase C's.
	 *
	 * The caller declares one beside the flag vector it already declares for
	 * the same reason, and hands it in; `refresh()` re-uploads a buffer only
	 * when its contents have actually changed, so the support sweep that
	 * genuinely moves `flags` still pays for `flags` and for nothing else.
	 */
	struct DeviceConstants
	{
		mfem::Array<KnotView> knots;
		mfem::Array<CoilView> coils;
		mfem::Array<int> flags;
		mfem::Array<int> offsetsOfShift;
		mfem::Vector shift;

		int coilCount = 0;
		int elementCount = 0;
		int shiftCount = 0;

		/*
		 * HOST SHADOWS, KEPT ONLY SO refresh() CAN TELL A CHANGED BUFFER FROM
		 * AN UNCHANGED ONE.
		 *
		 * The alternative is to read the device copy back and compare, which
		 * is the transfer this whole object exists to avoid -- so the
		 * comparison has to be against something that never left the host.
		 * They cost one host copy of each constant, which is the smaller half
		 * of what one avoided upload saves.
		 */
		std::vector<KnotView> knotShadow;
		std::vector<CoilView> coilShadow;
		std::vector<int> flagShadow;
		std::vector<int> offsetShadow;
		std::vector<double> shiftShadow;
		bool populated = false;

		/**
		 * Bring every buffer up to date with @a inputs, uploading only what
		 * changed.
		 *
		 * Call it inside a `meq::PhaseScope( meq::Phase::Upload )` so the
		 * uploads land in the phase they belong to -- which is what lets the
		 * instrument say whether "one upload" is what is happening.
		 */
		void refresh( DeviceIntegrandInputs const &inputs,
		              std::vector<KnotView> const &knotList,
		              std::vector<CoilView> const &coilList );
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
		DeviceIntegrandInputs const &inputs,
		DeviceConstants &constants );

}

#endif // MEQ_DEVICENPC_HPP
