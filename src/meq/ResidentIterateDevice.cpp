#include "ResidentIterate.hpp"

/*
 * THE LAUNCHES, AND NOTHING ELSE. Everything a reader looks up -- the
 * refusals, the phase-B tables, the gate -- is in ResidentIterate.cpp, which
 * stays under the naming check; this file is handed to nvcc and cannot be.
 *
 * COMPILED WITHOUT FMA CONTRACTION, for meq::DeviceRows' recorded reason and
 * not for speed: these reproduce host loops whose build has no -march and
 * therefore rounds a multiply and an add twice. CMakeLists.txt carries the
 * flag beside DeviceRowsDevice.cpp's.
 */

namespace meq
{

#if defined( MFEM_USE_CUDA_OR_HIP_LANG )

	namespace
	{
		/// The integrand's pointwise factor, which is every expression the five
		/// host loops write between `psi` and the accumulation -- in their
		/// order, operation for operation.
		MFEM_HOST_DEVICE inline void pointFactors( PlasmaLoopView const &in,
		                                           mfem::real_t pointWeight,
		                                           mfem::real_t radius,
		                                           mfem::real_t z,
		                                           mfem::real_t psi,
		                                           mfem::real_t &first,
		                                           mfem::real_t &second )
		{
			first = mfem::real_t( 0 );
			second = mfem::real_t( 0 );
			switch ( in.quantity )
			{
				case static_cast<int>( PlasmaQuantity::Current ):
					// ip.weight*Weight()*scaledF( R, z, psi )/R
					first = pointWeight*sourceF( in.source, in.knots, radius, z, psi )
					        /radius;
					break;
				case static_cast<int>( PlasmaQuantity::NormalisationCorner ):
				{
					// w = ip.weight*Weight()/R; w*dAxis and w*dBoundary.
					mfem::real_t dAxis = 0, dBoundary = 0;
					sourceNormalisationDerivatives( in.source, in.knots, radius, z,
					                                psi, dAxis, dBoundary );
					mfem::real_t const w = pointWeight/radius;
					first = w*dAxis;
					second = w*dBoundary;
					break;
				}
				case static_cast<int>( PlasmaQuantity::CurrentColumn ):
				{
					// derivative = scaledF/scale; -ip.weight*Weight()*derivative/R.
					// The unary minus binds to ip.weight on the host, and negating
					// a product's factor is exact, so -( w W ) is the same double.
					mfem::real_t const derivative =
						sourceF( in.source, in.knots, radius, z, psi )/in.scale;
					first = -pointWeight*derivative/radius;
					break;
				}
				case static_cast<int>( PlasmaQuantity::CurrentRow ):
					// ip.weight*Weight()*scaledDFdPsi/R, NOT negated.
					first = pointWeight
					        *sourceDFdPsi( in.source, in.knots, radius, z, psi )
					        /radius;
					break;
				case static_cast<int>( PlasmaQuantity::AxisColumn ):
				case static_cast<int>( PlasmaQuantity::BoundaryColumn ):
				{
					// -weight*derivative/radius, weight = ip.weight*Weight().
					mfem::real_t dAxis = 0, dBoundary = 0;
					sourceNormalisationDerivatives( in.source, in.knots, radius, z,
					                                psi, dAxis, dBoundary );
					mfem::real_t const derivative =
						in.quantity == static_cast<int>( PlasmaQuantity::AxisColumn )
						? dAxis : dBoundary;
					first = -pointWeight*derivative/radius;
					break;
				}
				default:
					break;
			}
		}

		/// `psi` at one quadrature point, as the host forms it: the shape
		/// contraction from zero in local-dof order, then `psi_c` added.
		MFEM_HOST_DEVICE inline mfem::real_t psiAt( PlasmaLoopView const &in,
		                                            int e, int q )
		{
			mfem::real_t psi = mfem::real_t( 0 );
			for ( int j = 0; j < in.dofs; ++j )
				psi += in.shape[ q*in.dofs + j ]
				       *in.state[ in.potentialStart + in.elementDofs[ e*in.dofs + j ] ];
			if ( in.hasShift )
				psi += in.shift[ in.shiftOffsets[ e ] + q ];
			return psi;
		}
	}

	bool plasmaElementKernel( PlasmaLoopView const &inputs,
	                          mfem::real_t *first, mfem::real_t *second )
	{
		if ( !mfem::Device::Allows( mfem::Backend::DEVICE_MASK ) )
			return false;
		if ( inputs.elements <= 0 )
			return false;

		// By value: every member is a scalar or a device pointer.
		PlasmaLoopView const in = inputs;
		mfem::forall( in.elements, [ = ] MFEM_HOST_DEVICE ( int e )
		{
			mfem::real_t a = mfem::real_t( 0 );
			mfem::real_t b = mfem::real_t( 0 );
			if ( in.flags[ e ] != 0 )
			{
				for ( int q = 0; q < in.quadPoints; ++q )
				{
					int const at = e*in.quadPoints + q;
					mfem::real_t const radius = in.radius[ at ];
					mfem::real_t const psi = psiAt( in, e, q );
					mfem::real_t pa = 0, pb = 0;
					pointFactors( in, in.weight[ at ], radius, in.z[ at ], psi,
					              pa, pb );
					a += pa;
					b += pb;
				}
			}
			first[ e ] = a;
			if ( second != nullptr )
				second[ e ] = b;
		} );
		return true;
	}

	bool plasmaScatterKernel( PlasmaLoopView const &inputs, mfem::real_t *out )
	{
		if ( !mfem::Device::Allows( mfem::Backend::DEVICE_MASK ) )
			return false;
		if ( inputs.elements <= 0 )
			return false;

		PlasmaLoopView const in = inputs;
		/*
		 * ONE THREAD PER ELEMENT, ACCUMULATING STRAIGHT INTO ITS OWN DOFS. The
		 * potential space is L2, so no other thread touches them, and the host
		 * loop adds to each dof once per quadrature point in point order -- so
		 * accumulating in place, point by point, is the host's sum in the
		 * host's order rather than a local array summed and written once.
		 */
		mfem::forall( in.elements, [ = ] MFEM_HOST_DEVICE ( int e )
		{
			if ( in.flags[ e ] == 0 )
				return;
			for ( int q = 0; q < in.quadPoints; ++q )
			{
				int const at = e*in.quadPoints + q;
				mfem::real_t const radius = in.radius[ at ];
				mfem::real_t const psi = psiAt( in, e, q );
				mfem::real_t factor = 0, unused = 0;
				pointFactors( in, in.weight[ at ], radius, in.z[ at ], psi,
				              factor, unused );
				for ( int j = 0; j < in.dofs; ++j )
					out[ in.potentialStart + in.elementDofs[ e*in.dofs + j ] ] +=
						factor*in.shape[ q*in.dofs + j ];
			}
		} );
		return true;
	}

	bool gatherKernel( mfem::Vector const &source, mfem::Array<int> const &indices,
	                   mfem::Vector &out )
	{
		auto const *d_source = source.Read();
		auto const *d_indices = indices.Read();
		auto *d_out = out.Write();
		mfem::forall( indices.Size(), [ = ] MFEM_HOST_DEVICE ( int k )
		{
			d_out[ k ] = d_source[ d_indices[ k ] ];
		} );
		return true;
	}

	bool orderedDotKernel( mfem::Vector const &source,
	                       mfem::Array<int> const &indices,
	                       mfem::Vector const &weights, mfem::Vector &out )
	{
		auto const *d_source = source.Read();
		auto const *d_indices = indices.Read();
		auto const *d_weights = weights.Read();
		auto *d_out = out.Write();
		int const count = indices.Size();
		// ONE THREAD, deliberately: a handful of entries, summed from zero in
		// the host's order. See orderedDotOnDevice().
		mfem::forall( 1, [ = ] MFEM_HOST_DEVICE ( int )
		{
			mfem::real_t total = mfem::real_t( 0 );
			for ( int k = 0; k < count; ++k )
				total += d_weights[ k ]*d_source[ d_indices[ k ] ];
			d_out[ 0 ] = total;
		} );
		return true;
	}

	bool copyEntriesKernel( mfem::Vector const &source,
	                        mfem::Array<int> const &indices,
	                        mfem::Vector &destination )
	{
		auto const *d_source = source.Read();
		auto const *d_indices = indices.Read();
		auto *d_destination = destination.ReadWrite();
		mfem::forall( indices.Size(), [ = ] MFEM_HOST_DEVICE ( int k )
		{
			d_destination[ d_indices[ k ] ] = d_source[ d_indices[ k ] ];
		} );
		return true;
	}

#else

	bool plasmaElementKernel( PlasmaLoopView const &, mfem::real_t *,
	                          mfem::real_t * )
	{
		return false;
	}

	bool plasmaScatterKernel( PlasmaLoopView const &, mfem::real_t * )
	{
		return false;
	}

	bool gatherKernel( mfem::Vector const &, mfem::Array<int> const &,
	                   mfem::Vector & )
	{
		return false;
	}

	bool orderedDotKernel( mfem::Vector const &, mfem::Array<int> const &,
	                       mfem::Vector const &, mfem::Vector & )
	{
		return false;
	}

	bool copyEntriesKernel( mfem::Vector const &, mfem::Array<int> const &,
	                        mfem::Vector & )
	{
		return false;
	}

#endif

}
