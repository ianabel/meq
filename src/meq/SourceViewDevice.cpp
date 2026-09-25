#include "SourceView.hpp"

/*
 * THE DEVICE TRANSLATION UNIT FOR meq::SourceView, and its first job is to
 * PROVE THE HEADER COMPILES AS DEVICE CODE.
 *
 * SourceView.hpp is all `MFEM_HOST_DEVICE` inline functions, and every one of
 * them is compiled by g++ on MEQ's ordinary path -- which says nothing whatever
 * about whether nvcc will take them. A `MFEM_HOST_DEVICE` function that calls
 * something host-only, or that uses a `std::` facility with no device form,
 * fails only in a translation unit compiled as CUDA, and MEQ has exactly one
 * of those per feature by design (see PointLocatorDevice.cpp for why:
 * mfem::forall static_asserts under g++). Without this file the header would be
 * device-annotated and untested, which is the same class of claim as a comment
 * asserting somebody else's code is guarded.
 *
 * AND IT IS A REAL ENTRY POINT RATHER THAN ONLY A COMPILE PROBE.
 * evaluateSourceOnDevice() is what a caller with the flux already resident
 * wants: F at N points without a host round trip. MEQ does not call it yet --
 * the batched local residual that would is upstream's
 * `CanBatchLinearResidual()`, which refuses on meq::SourceIntegrator being
 * nonlinear, and closing that is the interface this header exists to propose.
 * It is here so that the day the other side is built, the MEQ half is already
 * compiled, tested and device-clean.
 */

namespace meq
{

#if defined( MFEM_USE_CUDA_OR_HIP_LANG )

	bool evaluateSourceOnDevice( SourceView const &view,
	                             KnotView const *knots, int knotCount,
	                             mfem::real_t const *radius,
	                             mfem::real_t const *z,
	                             mfem::real_t const *psi,
	                             mfem::real_t *value,
	                             mfem::real_t *derivative,
	                             int count )
	{
		if ( !mfem::Device::Allows( mfem::Backend::DEVICE_MASK ) || count == 0 )
			return false;

		// Registered here rather than held, for the reason M-79 and M-130 both
		// are: a long-lived view onto a buffer whose device copy moves is the
		// defect this project has already written twice.
		mfem::Array<KnotView> knotBuf( const_cast<KnotView *>( knots ),
		                               knotCount > 0 ? knotCount : 1 );
		mfem::Vector radiusBuf( const_cast<mfem::real_t *>( radius ), count );
		mfem::Vector zBuf( const_cast<mfem::real_t *>( z ), count );
		mfem::Vector psiBuf( const_cast<mfem::real_t *>( psi ), count );
		mfem::Vector valueBuf( value, count );
		mfem::Vector derivativeBuf( derivative, count );

		auto const *d_knots  = knotBuf.Read();
		auto const *d_radius = radiusBuf.Read();
		auto const *d_z      = zBuf.Read();
		auto const *d_psi    = psiBuf.Read();
		auto *d_value        = valueBuf.Write();
		auto *d_derivative   = derivativeBuf.Write();

		// By value into the lambda: SourceView holds no pointers, deliberately,
		// so it captures whole and needs no device allocation of its own.
		SourceView const local = view;

		mfem::forall( count, [ = ] MFEM_HOST_DEVICE ( int i )
		{
			d_value[ i ] = sourceF( local, d_knots, d_radius[ i ], d_z[ i ],
			                        d_psi[ i ] );
			d_derivative[ i ] = sourceDFdPsi( local, d_knots, d_radius[ i ],
			                                  d_z[ i ], d_psi[ i ] );
		} );

		valueBuf.HostRead();
		derivativeBuf.HostRead();
		return true;
	}

#else

	bool evaluateSourceOnDevice( SourceView const &, KnotView const *, int,
	                             mfem::real_t const *, mfem::real_t const *,
	                             mfem::real_t const *, mfem::real_t *,
	                             mfem::real_t *, int )
	{
		return false;
	}

#endif

}
