#include "DeviceQuadrature.hpp"

/*
 * THE DEVICE TRANSLATION UNIT FOR THE PLASMA-CURRENT QUADRATURE. One file per
 * feature is compiled by the device compiler; PointLocatorDevice.cpp carries
 * the full account of why.
 */

namespace meq
{

#if defined( MFEM_USE_CUDA_OR_HIP_LANG )

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
	                          mfem::Vector &perElement )
	{
		if ( !mfem::Device::Allows( mfem::Backend::DEVICE_MASK ) )
			return false;
		if ( elements <= 0 || quadPoints <= 0 || dofs <= 0 )
			return false;

		auto const *d_u       = elementPotential.Read();
		auto const *d_w       = quadWeights.Read();
		auto const *d_detJ    = detJ.Read();
		auto const *d_x       = coordinates.Read();
		auto const *d_B       = shape.Read();
		auto const *d_flag    = inPlasma.Read();
		auto const *d_shift   = conductorShift.Read();
		auto const *d_knots   = knots.Read();
		auto *d_out           = perElement.Write();

		// By value: SourceView holds no pointers, so it captures whole.
		SourceView const local = view;
		int const nq = quadPoints;
		int const nd = dofs;
		int const sdim = spaceDimension;
		bool const shift = hasConductorShift;

		/*
		 * ONE THREAD PER ELEMENT, AND THE SUM INSIDE IT IS SEQUENTIAL OVER THE
		 * QUADRATURE POINTS IN RULE ORDER -- which is the order the host loop
		 * takes them in. The per-element partials are then summed by the
		 * CALLER in element order, exactly as assemblePlasmaCurrent() does and
		 * for the reason its own comment gives: an association that depends on
		 * the thread count takes psi_ax with it.
		 *
		 * LAYOUTS, read from MFEM's own headers rather than assumed:
		 *   B   is ( nq x nd ), column major   -> B[ q + nq*j ]
		 *   X   is ( nq x sdim x ne )          -> X[ q + nq*( d + sdim*e ) ]
		 *   detJ is ( nq x ne )                -> detJ[ q + nq*e ]
		 *   the E-vector is ( nd x ne )        -> u[ j + nd*e ]
		 */
		mfem::forall( elements, [ = ] MFEM_HOST_DEVICE ( int e )
		{
			if ( d_flag[ e ] == 0 )
			{
				d_out[ e ] = mfem::real_t( 0 );
				return;
			}

			mfem::real_t elementSum = mfem::real_t( 0 );
			for ( int q = 0; q < nq; ++q )
			{
				mfem::real_t psi = mfem::real_t( 0 );
				for ( int j = 0; j < nd; ++j )
					psi += d_B[ q + nq*j ]*d_u[ j + nd*e ];

				/*
				 * THE SOURCE SEES THE TOTAL FLUX AND THE SOLVER HOLDS THE
				 * REMAINDER. psi_c is added BEFORE F is evaluated and never
				 * afterwards -- COIL-SUBTRACTION-PLAN.md CS-2 -- and the
				 * offsets were checked host-side to be e*nq, so this index is
				 * the cache's own `values[ offsets[ e ] + q ]`.
				 */
				if ( shift )
					psi += d_shift[ q + nq*e ];

				mfem::real_t const radius = d_x[ q + nq*( 0 + sdim*e ) ];
				mfem::real_t const z      = d_x[ q + nq*( 1 + sdim*e ) ];

				elementSum += d_w[ q ]*d_detJ[ q + nq*e ]
				              *sourceF( local, d_knots, radius, z, psi )
				              /radius;
			}
			d_out[ e ] = elementSum;
		} );

		/*
		 * LEFT DEVICE-VALID, AND THE CALLER FUNNELS -- see
		 * dotCompressedRowsOnDevice() for the argument. The in-element-order
		 * sum is host arithmetic today and is 39 KB on the shipped machine
		 * case against the state's 698 KB, but a kernel that forced the read
		 * would fix that choice for a resident step that should not want it.
		 */
		return true;
	}

#else

	bool plasmaCurrentKernel( mfem::Vector const &, mfem::Vector const &,
	                          mfem::Vector const &, mfem::Vector const &,
	                          mfem::Array<mfem::real_t> const &,
	                          mfem::Array<int> const &, mfem::Vector const &,
	                          bool, SourceView const &,
	                          mfem::Array<KnotView> const &,
	                          int, int, int, int, mfem::Vector & )
	{
		return false;
	}

#endif

}
