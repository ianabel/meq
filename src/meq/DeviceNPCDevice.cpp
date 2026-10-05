#include "DeviceNPC.hpp"

#include "Crossings.hpp"
#include "Source.hpp"

/*
 * THE DEVICE TRANSLATION UNIT FOR MEQ'S NPC OPERATOR.
 *
 * It is here and not in GradShafranov.cpp because instantiating upstream's
 * templated NPCResidual()/NPCGradient() pulls in HDGNLElementQuadrature, which
 * contains mfem::forall -- and that static_asserts under g++. One file per
 * feature is compiled by the device compiler; PointLocatorDevice.cpp carries
 * the full account.
 */

namespace meq
{

#if defined( MFEM_USE_CUDA_OR_HIP_LANG )

	namespace
	{

		class DeviceNPCOperator : public mfem::DarcyNPCOperator
		{
			public:
				DeviceNPCOperator( mfem::DarcyHybridization &hybridization,
				                   mfem::Array<int> const &offsetsIn,
				                   mfem::BlockVector const &loadIn,
				                   DeviceIntegrandInputs const &inputs,
				                   SourceView const &view,
				                   mfem::real_t permeabilityIn,
				                   DeviceConstants const &constantsIn )
					: mfem::DarcyNPCOperator( hybridization, offsetsIn, loadIn ),
					  hyb( &hybridization ), load( &loadIn ),
					  source( inputs.source ), permeability( permeabilityIn ),
					  cached( view ), constants( &constantsIn )
				{
					offsets = offsetsIn;
					locOffsets = LocalOffsets();

					/*
					 * THE FIVE CONSTANTS ARE NOT BUILT HERE ANY MORE, AND
					 * THAT IS PHASE B LEAVING PHASE C.
					 *
					 * They were members, filled and uploaded in this
					 * constructor -- which is correct once and was happening
					 * SIXTEEN times, because this operator is rebuilt
					 * whenever the hybridization is re-finalised
					 * -> M-226. meq::DeviceConstants now owns them, the
					 * caller declares it beside the flag vector it already
					 * declares for the same reason, and it outlives every
					 * rebuild.
					 */
				}

				/*
				 * THE TWO INTEGRANDS ARE FUNCTOR STRUCTS AND NOT LAMBDAS, AND
				 * THAT IS nvcc's DOING RATHER THAN A STYLE CHOICE.
				 *
				 * An extended __host__ __device__ lambda may not sit in a
				 * private member function, may not sit in one with a deduced
				 * return type, and may not first-capture inside a constexpr-if
				 * -- the last of which is what stopped upstream's tagged kernel
				 * compiling. A lambda cannot be RETURNED from a function at all
				 * without `auto`, so "build the callable and hand it over" is
				 * closed. A plain struct with a MFEM_HOST_DEVICE operator() is
				 * not an extended lambda, carries none of those restrictions,
				 * and is nameable -- so it can be built once and passed twice.
				 *
				 * **A HOST-ONLY BUILD SAYS NOTHING ABOUT WHETHER ANY OF THIS IS
				 * LEGAL**, which is the transferable part: every one of those
				 * rules is enforced only by the device compiler.
				 */
				struct ResidualIntegrand
				{
					SourceView view;
					KnotView const *knots;
					CoilView const *coils;
					int const *flags;
					int const *shiftOffsets;
					mfem::real_t const *shift;
					mfem::real_t permeability;
					int coilCount;
					bool hasShift;

					MFEM_HOST_DEVICE mfem::real_t operator()(
						int element, int point, mfem::real_t radius,
						mfem::real_t z, mfem::real_t u ) const
					{
						mfem::real_t const psi = hasShift
							? u + shift[ shiftOffsets[ element ] + point ] : u;

						// ( plasma, gated by the element flag ) + coil. One
						// expression for sourceValue()'s two branches; the gate
						// is the only difference between them.
						mfem::real_t total =
							coilF( coils, coilCount, permeability, radius, z );
						if ( flags[ element ] != 0 )
							total += sourceF( view, knots, radius, z, psi );

						// The weight and the basis contraction are upstream's,
						// so this returns the integrand and nothing else.
						return -total/radius;
					}
				};

				struct GradientIntegrand
				{
					SourceView view;
					KnotView const *knots;
					int const *flags;
					int const *shiftOffsets;
					mfem::real_t const *shift;
					bool hasShift;

					MFEM_HOST_DEVICE mfem::real_t operator()(
						int element, int point, mfem::real_t radius,
						mfem::real_t z, mfem::real_t u ) const
					{
						// **THE COILS CONTRIBUTE EXACTLY ZERO HERE**, which is
						// CoilAugmentedNormalisedSource::dFdPsi()'s own
						// definition: a coil current is amperes and depends on
						// no flux. So the gradient integrand is the plasma term
						// alone and the element gate is the whole of it.
						if ( flags[ element ] == 0 )
							return mfem::real_t( 0 );

						mfem::real_t const psi = hasShift
							? u + shift[ shiftOffsets[ element ] + point ] : u;
						return -sourceDFdPsi( view, knots, radius, z, psi )
						       /radius;
					}
				};

				/*
				 * THE VIEW IS REFRESHED ON EVERY CALL, because psi_ax, psi_bnd
				 * and the current scale are unknowns of the bordered Newton and
				 * move between residual evaluations inside one step. The knots
				 * do not move within a solve, so only the scalars are taken.
				 */
				void refresh() const
				{
					SourceView fresh{};
					std::vector<KnotView> scratch;
					if ( source != nullptr
					     && makePlasmaSourceView( *source, fresh, scratch ) )
					{
						cached.psiAxis = fresh.psiAxis;
						cached.psiBoundary = fresh.psiBoundary;
						cached.currentScale = fresh.currentScale;
						cached.supportAxis = fresh.supportAxis;
						cached.supportBoundary = fresh.supportBoundary;
						cached.confineToPlasma = fresh.confineToPlasma;
					}
				}

				ResidualIntegrand residualIntegrand() const
				{
					refresh();
					// Counted, not bare: these are phase B's payload being
					// re-read once per residual, and whether that MOVES
					// anything is the measurement, not an assumption.
					return ResidualIntegrand{
						cached,
						countedDeviceRead( constants->knots, __LINE__ ),
						countedDeviceRead( constants->coils, __LINE__ ),
						countedDeviceRead( constants->flags, __LINE__ ),
						countedDeviceRead( constants->offsetsOfShift, __LINE__ ),
						countedDeviceRead( constants->shift, __LINE__ ),
						permeability, constants->coilCount,
						constants->shiftCount > 0 };
				}

				GradientIntegrand gradientIntegrand() const
				{
					refresh();
					return GradientIntegrand{
						cached,
						countedDeviceRead( constants->knots, __LINE__ ),
						countedDeviceRead( constants->flags, __LINE__ ),
						countedDeviceRead( constants->offsetsOfShift, __LINE__ ),
						countedDeviceRead( constants->shift, __LINE__ ),
						constants->shiftCount > 0 };
				}

				void Mult( mfem::Vector const &x, mfem::Vector &y ) const override
				{
					mfem::BlockVector const xb(
						const_cast<mfem::Vector &>( x ), offsets );
					/*
					 * THE OUTPUT AND THE LOCAL STATE ARE MADE DEVICE-VALID
					 * BEFORE THEIR BLOCKS ARE WRITTEN, OR THE TWO
					 * SyncFromBlocks() BELOW ARE DOWNLOADS.
					 *
					 * `BlockVector::SyncFromBlocks()` is `block.SyncAlias(
					 * base )` per block: it makes each BLOCK match its BASE,
					 * and `SyncAlias_` copies the block's device range down
					 * whenever the base says host. `y` arrives host-valid --
					 * every residual vector MEQ declares is a fresh
					 * `mfem::Vector v( n )` -- and `xLoc` is host-valid after
					 * `Update()`, so every residual evaluation downloaded five
					 * blocks the kernels had just written: 147 copies and
					 * 45 MB over `machine-f-diiid`, upstream's
					 * `DarcyNPCSolver::ArrayMult()` shape in MEQ's own copy
					 * (M-237). `Write()` marks the device copy live without
					 * moving anything, and a base that says device leaves the
					 * blocks nothing to catch up to -- while a block written on
					 * the HOST is then carried UP, so the flags stay honest on
					 * either route.
					 *
					 * Before `yb` is taken: an alias copies its base's flags
					 * when it is made.
					 */
					y.UseDevice( true );
					y.Write();
					mfem::BlockVector yb( y, offsets );

					/*
					 * AND EVERY BLOCK IS MADE DEVICE-CAPABLE BY NAME, because
					 * the copies below are `Vector::operator=`, which runs
					 * where either side's UseDevice() says -- and a block keeps
					 * the flags it was MADE with. `Update()` at an unchanged
					 * size returns early and keeps blocks made before the base
					 * was device-capable, so the copy ran on the host into a
					 * page the base had just been told was the device's:
					 * `mfem::Device( "debug" )` faulted in its memcpy, and under
					 * CUDA it was a silent host write repaired by an upload.
					 */
					for ( int b = 0; b < yb.NumBlocks(); ++b )
						yb.GetBlock( b ).UseDevice( true );

					xLoc.Update( locOffsets );
					xLoc.UseDevice( true );
					xLoc.Write();
					for ( int b = 0; b < xLoc.NumBlocks(); ++b )
						xLoc.GetBlock( b ).UseDevice( true );
					xLoc.GetBlock( 0 ) = xb.GetBlock( 0 );
					xLoc.GetBlock( 1 ) = xb.GetBlock( 1 );
					/*
					 * HYGIENE, AND NOT LOAD BEARING -- WHICH THIS COMMENT
					 * GOT WRONG AND UPSTREAM CAUGHT.
					 *
					 * This routine writes `xLoc`'s two blocks through aliases
					 * with flags of their own, so `xLoc` itself is left
					 * believing a host copy is current that is not. Carrying
					 * them up is the right thing to do with a vector whose
					 * flags would otherwise be a lie, and upstream keeps the
					 * same call in their own copy on the same grounds.
					 *
					 * **WHAT IT DOES NOT DO IS FIX A READ, BECAUSE NOTHING
					 * READS THE BASE.** `NPCResidual()` and both
					 * `NPCGradient()` overloads take `x.GetBlock( 0 )` and
					 * `x.GetBlock( 1 )` -- read out of
					 * `fem/darcy/darcyhybridization.{hpp,cpp}` rather than
					 * assumed -- so a sync that carries the blocks UP to the
					 * base is inert for every consumer. The version of this
					 * comment that stood here said `NPCGradient()` reads the
					 * base and faults on it. It does not, and the fault MEQ
					 * actually observed was `Df_data += Df_lin_data` several
					 * frames further in: a different defect with a different
					 * cause, and upstream's to fix.
					 *
					 * **AND THE STEP-COUNT DEFECT THIS WAS CREDITED WITH
					 * CURING IS `b_tr`.** `DarcyNPCSolver` reuses it across
					 * calls; `NPCReduce()`'s host route wrote it through an
					 * alias and never carried that back, so from the second
					 * Newton step on the reduced right-hand side was wrong.
					 * Upstream reproduced it, added THIS sync, and the step
					 * count did not move. Fixed at `ce744b84de`.
					 * MEASUREMENTS.md M-207 and M-210.
					 */
					xLoc.SyncFromBlocks();
					/*
					 * AND THE RESIDUAL'S LOCAL BLOCKS GO TO THE DEVICE FOR
					 * THE SAME REASON AS `xLoc`.
					 *
					 * Upstream's batched residual writes `rLoc`'s blocks on
					 * the device and then syncs each block to its BASE -- so a
					 * host-valid `rLoc` had both blocks downloaded inside
					 * `NonlinearResidualFinish()`, and the copies into `yb`
					 * below uploaded them again: 28 of each over
					 * `machine-f-diiid`, the potential block's re-upload
					 * alone 6.5 MB. Upstream named the cause in its reply to
					 * HDG-DEVICE-BUS-FROM-MEQ.md -- *MEQ's `r` has
					 * `UseDevice()` false* -- and it is the caller's to fix.
					 */
					rLoc.Update( locOffsets );
					rLoc.UseDevice( true );
					rLoc.Write();
					for ( int b = 0; b < rLoc.NumBlocks(); ++b )
						rLoc.GetBlock( b ).UseDevice( true );
					/*
					 * AND THE TRACE ROW, WHICH IS THE SAME DEFECT ONE ARGUMENT
					 * OVER. With no prolongation `NPCResidualStart()` makes
					 * its local trace row an ALIAS of `rTr`, the kernel writes
					 * the alias, and `NPCResidualFinish()` syncs the alias to
					 * this base -- so a host-valid `rTr` downloaded the row on
					 * every residual and the copy into `yb` below uploaded it
					 * again. Sized here because `NPCResidualStart()` resizes
					 * only on a mismatch, and a resize would reallocate past
					 * the `Write()`; with a prolongation the true-dof size is
					 * still the iterate's trace block.
					 */
					rTr.SetSize( xb.GetBlock( 2 ).Size() );
					rTr.UseDevice( true );
					rTr.Write();

					hyb->NPCResidual( *load, xLoc, xb.GetBlock( 2 ), rLoc, rTr,
					                  residualIntegrand() );

					yb.GetBlock( 0 ) = rLoc.GetBlock( 0 );
					yb.GetBlock( 1 ) = rLoc.GetBlock( 1 );
					yb.GetBlock( 2 ) = rTr;

					/*
					 * BOTH HOPS, COPIED FROM DarcyNPCOperator::Mult RATHER THAN
					 * REASONED ABOUT. yb is an ALIAS of y and each block is an
					 * alias of yb -- three Memory objects with three sets of
					 * flags. Upstream's own comment records that without the
					 * first, a fresh Newton work vector reads zero and M-79's
					 * symptom comes back: a converged Newton that takes no
					 * steps.
					 *
					 * **WHAT `SyncFromBlocks()` DOES IS MAKE EACH BLOCK MATCH
					 * ITS BASE, NOT CARRY THE BLOCKS UP** -- it is
					 * `block.SyncAlias( base )`, and `SyncAlias_` reconciles the
					 * alias to the base's validity, copying the block's other
					 * copy across where the base says a location is live and the
					 * block does not. This comment said *up*. With `y` arriving
					 * host-valid, what happened was the device-written blocks
					 * being copied DOWN to match it -- the right values, paid
					 * for with a download per block per residual. See the top
					 * of this function for why `y` is now device-valid going
					 * in. M-238.
					 */
					yb.SyncFromBlocks();
					yb.SyncAliasMemory( y );
				}

				mfem::Operator &GetGradient( mfem::Vector const &x ) const override
				{
					mfem::BlockVector const xb(
						const_cast<mfem::Vector &>( x ), offsets );

					// Device-valid before the blocks are written, and the blocks
					// device-capable by name; see Mult().
					xLoc.Update( locOffsets );
					xLoc.UseDevice( true );
					xLoc.Write();
					for ( int b = 0; b < xLoc.NumBlocks(); ++b )
						xLoc.GetBlock( b ).UseDevice( true );
					xLoc.GetBlock( 0 ) = xb.GetBlock( 0 );
					xLoc.GetBlock( 1 ) = xb.GetBlock( 1 );
					/*
					 * HYGIENE, AND NOT LOAD BEARING -- WHICH THIS COMMENT
					 * GOT WRONG AND UPSTREAM CAUGHT.
					 *
					 * This routine writes `xLoc`'s two blocks through aliases
					 * with flags of their own, so `xLoc` itself is left
					 * believing a host copy is current that is not. Carrying
					 * them up is the right thing to do with a vector whose
					 * flags would otherwise be a lie, and upstream keeps the
					 * same call in their own copy on the same grounds.
					 *
					 * **WHAT IT DOES NOT DO IS FIX A READ, BECAUSE NOTHING
					 * READS THE BASE.** `NPCResidual()` and both
					 * `NPCGradient()` overloads take `x.GetBlock( 0 )` and
					 * `x.GetBlock( 1 )` -- read out of
					 * `fem/darcy/darcyhybridization.{hpp,cpp}` rather than
					 * assumed -- so a sync that carries the blocks UP to the
					 * base is inert for every consumer. The version of this
					 * comment that stood here said `NPCGradient()` reads the
					 * base and faults on it. It does not, and the fault MEQ
					 * actually observed was `Df_data += Df_lin_data` several
					 * frames further in: a different defect with a different
					 * cause, and upstream's to fix.
					 *
					 * **AND THE STEP-COUNT DEFECT THIS WAS CREDITED WITH
					 * CURING IS `b_tr`.** `DarcyNPCSolver` reuses it across
					 * calls; `NPCReduce()`'s host route wrote it through an
					 * alias and never carried that back, so from the second
					 * Newton step on the reduced right-hand side was wrong.
					 * Upstream reproduced it, added THIS sync, and the step
					 * count did not move. Fixed at `ce744b84de`.
					 * MEASUREMENTS.md M-207 and M-210.
					 */
					xLoc.SyncFromBlocks();

					mfem::Operator &s = hyb->NPCGradient(
						xLoc, xb.GetBlock( 2 ), gradientIntegrand() );
					jacobian.reset( new Jacobian( *hyb, s, offsets, locOffsets ) );
					return *jacobian;
				}

			private:



				mfem::DarcyHybridization *hyb;
				mfem::BlockVector const *load;
				Source const *source;
				mfem::real_t permeability;

				mfem::Array<int> offsets, locOffsets;

				/// Phase B's payload, BORROWED. The caller owns it so that it
				/// outlives this operator's rebuilds -- see DeviceConstants.
				DeviceConstants const *constants;

				mutable SourceView cached;
				mutable mfem::BlockVector xLoc, rLoc;
				mutable mfem::Vector rTr;
				mutable std::unique_ptr<Jacobian> jacobian;
		};

	}

	std::unique_ptr<mfem::DarcyNPCOperator> makeDeviceNPCOperator(
		mfem::DarcyHybridization &hybridization,
		mfem::Array<int> const &offsets,
		mfem::BlockVector const &load,
		DeviceIntegrandInputs const &inputs,
		DeviceConstants &constants )
	{
		if ( !mfem::Device::Allows( mfem::Backend::DEVICE_MASK ) )
			return nullptr;
		if ( inputs.source == nullptr )
			return nullptr;

		SourceView view{};
		std::vector<KnotView> knotList;
		if ( !makePlasmaSourceView( *inputs.source, view, knotList ) )
			return nullptr;

		std::vector<CoilView> coilList;
		mfem::real_t permeability = 0;
		if ( !makeCoilViews( *inputs.source, coilList, permeability ) )
			return nullptr;

		// PHASE B. Uploads only what has actually changed since the last
		// rebuild, which on every rebuild but the first and the support
		// sweeps is nothing at all.
		constants.refresh( inputs, knotList, coilList );

		return std::make_unique<DeviceNPCOperator>(
			hybridization, offsets, load, inputs, view, permeability,
			constants );
	}

#else

	std::unique_ptr<mfem::DarcyNPCOperator> makeDeviceNPCOperator(
		mfem::DarcyHybridization &, mfem::Array<int> const &,
		mfem::BlockVector const &, DeviceIntegrandInputs const &,
		DeviceConstants & )
	{
		return nullptr;
	}

#endif

	/*
	 * PHASE B, AND IT IS OUTSIDE THE #if ON PURPOSE.
	 *
	 * A host build has no kernels to feed, so every comparison below finds
	 * nothing to do and every Read() is a no-op on an unregistered buffer --
	 * but the function still has to EXIST, because the caller declares a
	 * DeviceConstants whether or not a device compiler was present. The five
	 * device translation units each compile to a host stub for the same
	 * reason; this is the one piece of them that is not a stub because it is
	 * not device code.
	 */
	void DeviceConstants::refresh( DeviceIntegrandInputs const &inputs,
	                               std::vector<KnotView> const &knotList,
	                               std::vector<CoilView> const &coilList )
	{
		static_assert( std::is_trivially_copyable<KnotView>::value,
		               "KnotView is compared by memcmp below" );
		static_assert( std::is_trivially_copyable<CoilView>::value,
		               "CoilView is compared by memcmp below" );

		/*
		 * A ZERO-SIZED mfem::Array GIVES Read() NOTHING TO HAND A KERNEL, so
		 * an empty set still allocates one entry the kernel will not look at
		 * and the COUNTS are what gate it. Same convention the operator used
		 * when it owned these, and the kernels already rely on it.
		 */
		auto const changed = []( auto const &shadow, auto const *source,
		                         std::size_t count ) -> bool
		{
			if ( shadow.size() != count )
				return true;
			if ( count == 0 )
				return false;
			return std::memcmp( shadow.data(), source,
			                    count*sizeof( *source ) ) != 0;
		};

		/*
		 * EVERY FILL BELOW OPENS WITH HostWrite(), AND M-211 IS WHY.
		 *
		 * These buffers are now REFRESHED IN PLACE rather than freshly
		 * constructed, and `SetSize()` at an unchanged size reallocates
		 * nothing -- so `Array::operator[]` is a raw host write into a buffer
		 * whose live copy is on the device, and the `Read()` that follows sees
		 * a device copy MFEM still believes is current and copies nothing.
		 * The support changes between sweeps, the device keeps the old flags,
		 * and the run converges to a neighbouring answer: measured at
		 * `psi_ax = 3.759845e-01` against the host's `3.759851e-01`.
		 *
		 * **IT IS MEQ'S OWN meq::DeviceRows::upload() SHAPE, MET AGAIN** -- *a
		 * buffer refreshed in place is device-valid from its own previous use*
		 * -- and the first version of this function reproduced it exactly,
		 * because the code it was hoisted out of was CONSTRUCTING these arrays
		 * each time and a fresh mfem::Array has UseDevice() false. *A fill that
		 * is safe in a constructor is not safe in a refresh.*
		 */
		bool const first = !populated;

		if ( first || changed( knotShadow, knotList.data(), knotList.size() ) )
		{
			knotShadow = knotList;
			knots.SetSize( knotList.empty()
			               ? 1 : static_cast<int>( knotList.size() ) );
			knots.HostWrite();
			for ( std::size_t k = 0; k < knotList.size(); ++k )
				knots[ static_cast<int>( k ) ] = knotList[ k ];
			if ( knotList.empty() )
				knots[ 0 ] = KnotView{ 0, 0, 0 };
			countedDeviceRead( knots, __LINE__ );
		}

		if ( first || changed( coilShadow, coilList.data(), coilList.size() ) )
		{
			coilShadow = coilList;
			coilCount = static_cast<int>( coilList.size() );
			coils.SetSize( coilCount > 0 ? coilCount : 1 );
			coils.HostWrite();
			for ( int k = 0; k < coilCount; ++k )
				coils[ k ] = coilList[ static_cast<std::size_t>( k ) ];
			if ( coilCount == 0 )
				coils[ 0 ] = CoilView{ 0, 0, 0, 0, 0 };
			countedDeviceRead( coils, __LINE__ );
		}

		/*
		 * THE ONE THAT GENUINELY MOVES, AND IT MOVES ONCE PER SUPPORT SWEEP.
		 * XP-2 freezes the plasma support for the duration of a solve, so
		 * within a sweep this comparison finds nothing; between sweeps the
		 * fill has re-decided and the upload is real work rather than waste.
		 * That is the whole point of comparing instead of re-uploading: the
		 * buffer that changes still pays, and the four that do not stop
		 * paying.
		 */
		elementCount = inputs.elementCount;
		std::vector<int> freshFlags(
			static_cast<std::size_t>( elementCount > 0 ? elementCount : 0 ) );
		for ( int e = 0; e < elementCount; ++e )
			freshFlags[ static_cast<std::size_t>( e ) ] =
				inputs.elementInPlasma ? inputs.elementInPlasma[ e ] : 1;

		if ( first || changed( flagShadow, freshFlags.data(),
		                       freshFlags.size() ) )
		{
			flagShadow = freshFlags;
			flags.SetSize( elementCount > 0 ? elementCount : 1 );
			flags.HostWrite();
			for ( int e = 0; e < elementCount; ++e )
				flags[ e ] = freshFlags[ static_cast<std::size_t>( e ) ];
			if ( elementCount == 0 )
				flags[ 0 ] = 1;
			countedDeviceRead( flags, __LINE__ );
		}

		std::vector<int> freshOffsets(
			static_cast<std::size_t>( elementCount > 0 ? elementCount : 0 ) );
		for ( int e = 0; e < elementCount; ++e )
			freshOffsets[ static_cast<std::size_t>( e ) ] =
				inputs.conductorOffsets ? inputs.conductorOffsets[ e ] : 0;

		if ( first || changed( offsetShadow, freshOffsets.data(),
		                       freshOffsets.size() ) )
		{
			offsetShadow = freshOffsets;
			offsetsOfShift.SetSize( elementCount > 0 ? elementCount : 1 );
			offsetsOfShift.HostWrite();
			for ( int e = 0; e < elementCount; ++e )
				offsetsOfShift[ e ] =
					freshOffsets[ static_cast<std::size_t>( e ) ];
			if ( elementCount == 0 )
				offsetsOfShift[ 0 ] = 0;
			countedDeviceRead( offsetsOfShift, __LINE__ );
		}

		shiftCount = inputs.conductorValueCount;
		if ( first || changed( shiftShadow, inputs.conductorValues,
		                       static_cast<std::size_t>(
		                           shiftCount > 0 ? shiftCount : 0 ) ) )
		{
			shiftShadow.clear();
			if ( inputs.conductorValues != nullptr && shiftCount > 0 )
				shiftShadow.assign( inputs.conductorValues,
				                    inputs.conductorValues + shiftCount );
			shift.SetSize( shiftCount > 0 ? shiftCount : 1 );
			shift.HostWrite();
			for ( int k = 0; k < shiftCount; ++k )
				shift( k ) = static_cast<mfem::real_t>(
					inputs.conductorValues[ k ] );
			if ( shiftCount == 0 )
				shift( 0 ) = mfem::real_t( 0 );
			shift.UseDevice( true );
			countedDeviceRead( shift, __LINE__ );
		}

		populated = true;
	}

}
