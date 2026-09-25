#include "DeviceNPC.hpp"

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
				                   std::vector<KnotView> const &knotList,
				                   std::vector<CoilView> const &coilList,
				                   mfem::real_t permeabilityIn )
					: mfem::DarcyNPCOperator( hybridization, offsetsIn, loadIn ),
					  hyb( &hybridization ), load( &loadIn ),
					  source( inputs.source ), permeability( permeabilityIn ),
					  cached( view )
				{
					offsets = offsetsIn;
					locOffsets = LocalOffsets();

					// The three fixed arrays, uploaded once. A zero-sized
					// mfem::Array gives Read() nothing to hand a kernel, so an
					// empty set still allocates one entry the kernel will not
					// look at -- the counts are what gate it.
					knots.SetSize( knotList.empty() ? 1
					               : static_cast<int>( knotList.size() ) );
					for ( std::size_t k = 0; k < knotList.size(); ++k )
						knots[ static_cast<int>( k ) ] = knotList[ k ];
					if ( knotList.empty() )
						knots[ 0 ] = KnotView{ 0, 0, 0 };

					coilCount = static_cast<int>( coilList.size() );
					coils.SetSize( coilCount > 0 ? coilCount : 1 );
					for ( int k = 0; k < coilCount; ++k )
						coils[ k ] = coilList[ static_cast<std::size_t>( k ) ];
					if ( coilCount == 0 )
						coils[ 0 ] = CoilView{ 0, 0, 0, 0, 0 };

					elementCount = inputs.elementCount;
					flags.SetSize( elementCount > 0 ? elementCount : 1 );
					for ( int e = 0; e < elementCount; ++e )
						flags[ e ] = inputs.elementInPlasma
							? inputs.elementInPlasma[ e ] : 1;
					if ( elementCount == 0 )
						flags[ 0 ] = 1;

					shiftCount = inputs.conductorValueCount;
					shift.SetSize( shiftCount > 0 ? shiftCount : 1 );
					shift.HostWrite();
					for ( int k = 0; k < shiftCount; ++k )
						shift( k ) = static_cast<mfem::real_t>(
							inputs.conductorValues[ k ] );
					if ( shiftCount == 0 )
						shift( 0 ) = mfem::real_t( 0 );
					shift.UseDevice( true );

					// The cache's own offsets, so the kernel indexes it exactly
					// as conductorShiftAt() does rather than assuming e*nq.
					offsetsOfShift.SetSize( elementCount > 0 ? elementCount : 1 );
					for ( int e = 0; e < elementCount; ++e )
						offsetsOfShift[ e ] = inputs.conductorOffsets
							? inputs.conductorOffsets[ e ] : 0;
					if ( elementCount == 0 )
						offsetsOfShift[ 0 ] = 0;

					knots.Read();
					coils.Read();
					flags.Read();
					offsetsOfShift.Read();
					shift.Read();
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
					return ResidualIntegrand{
						cached, knots.Read(), coils.Read(), flags.Read(),
						offsetsOfShift.Read(), shift.Read(), permeability,
						coilCount, shiftCount > 0 };
				}

				GradientIntegrand gradientIntegrand() const
				{
					refresh();
					return GradientIntegrand{
						cached, knots.Read(), flags.Read(),
						offsetsOfShift.Read(), shift.Read(), shiftCount > 0 };
				}

				void Mult( mfem::Vector const &x, mfem::Vector &y ) const override
				{
					mfem::BlockVector const xb(
						const_cast<mfem::Vector &>( x ), offsets );
					mfem::BlockVector yb( y, offsets );

					xLoc.Update( locOffsets );
					xLoc.GetBlock( 0 ) = xb.GetBlock( 0 );
					xLoc.GetBlock( 1 ) = xb.GetBlock( 1 );
					rLoc.Update( locOffsets );

					hyb->NPCResidual( *load, xLoc, xb.GetBlock( 2 ), rLoc, rTr,
					                  residualIntegrand() );

					yb.GetBlock( 0 ) = rLoc.GetBlock( 0 );
					yb.GetBlock( 1 ) = rLoc.GetBlock( 1 );
					yb.GetBlock( 2 ) = rTr;

					/*
					 * BOTH HOPS, COPIED FROM DarcyNPCOperator::Mult RATHER THAN
					 * REASONED ABOUT. yb is an ALIAS of y and each block is an
					 * alias of yb -- three Memory objects with three sets of
					 * flags. SyncFromBlocks() carries the blocks up to yb and
					 * SyncAliasMemory() carries yb up to y. Upstream's own
					 * comment records that without the first, a fresh Newton
					 * work vector reads zero and M-79's symptom comes back: a
					 * converged Newton that takes no steps.
					 */
					yb.SyncFromBlocks();
					yb.SyncAliasMemory( y );
				}

				mfem::Operator &GetGradient( mfem::Vector const &x ) const override
				{
					mfem::BlockVector const xb(
						const_cast<mfem::Vector &>( x ), offsets );

					xLoc.Update( locOffsets );
					xLoc.GetBlock( 0 ) = xb.GetBlock( 0 );
					xLoc.GetBlock( 1 ) = xb.GetBlock( 1 );

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
				mfem::Array<KnotView> knots;
				mfem::Array<CoilView> coils;
				mfem::Array<int> flags;
				mfem::Array<int> offsetsOfShift;
				mfem::Vector shift;
				int coilCount = 0;
				int elementCount = 0;
				int shiftCount = 0;

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
		DeviceIntegrandInputs const &inputs )
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

		return std::make_unique<DeviceNPCOperator>(
			hybridization, offsets, load, inputs, view, knotList, coilList,
			permeability );
	}

#else

	std::unique_ptr<mfem::DarcyNPCOperator> makeDeviceNPCOperator(
		mfem::DarcyHybridization &, mfem::Array<int> const &,
		mfem::BlockVector const &, DeviceIntegrandInputs const & )
	{
		return nullptr;
	}

#endif

}
