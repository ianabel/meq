#include "ResidentIterate.hpp"

#include "mfem/general/reducers.hpp"

#include <limits>

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

	bool nodalPeakKernel( mfem::Vector const &source, int begin, int end,
	                      mfem::real_t &best, int &bestIndex )
	{
		best = -std::numeric_limits<mfem::real_t>::infinity();
		bestIndex = -1;
		int const count = end - begin;
		if ( count <= 0 )
			return true;

		auto const *d_source = source.Read();

		// THE VALUE FIRST, which is exact in any order: fmax is a selection
		// and not arithmetic, and it skips NaN exactly as the host's `>` does.
		mfem::real_t peak = best;
		mfem::Array<mfem::real_t> valueWork;
		mfem::reduce( count, peak,
		              [ = ] MFEM_HOST_DEVICE ( int i, mfem::real_t &r )
		              {
		                      r = fmax( r, d_source[ begin + i ] );
		              },
		              mfem::MaxReducer<mfem::real_t>{}, true, valueWork );
		if ( !( peak > best ) )
			return true;

		// THEN THE EARLIEST INDEX HOLDING IT, which is the walk's tie-break.
		int first = end;
		mfem::Array<int> indexWork;
		mfem::reduce( count, first,
		              [ = ] MFEM_HOST_DEVICE ( int i, int &r )
		              {
		                      if ( d_source[ begin + i ] == peak && begin + i < r )
		                              r = begin + i;
		              },
		              mfem::MinReducer<int>{}, true, indexWork );
		best = peak;
		bestIndex = first;
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

	bool plasmaLoopGeometryKernel( mfem::Vector const &coordinates,
	                               mfem::Vector const &determinants,
	                               mfem::Vector const &quadWeights,
	                               int elements, int quadPoints, int dofs,
	                               mfem::Vector &radius, mfem::Vector &z,
	                               mfem::Vector &weight,
	                               mfem::Array<int> &elementDofs )
	{
		// X is ( nq x sdim x ne ) and detJ ( nq x ne ), as
		// DeviceQuadratureDevice.cpp reads them from MFEM's headers.
		auto const *d_x = coordinates.Read();
		auto const *d_detJ = determinants.Read();
		auto const *d_w = quadWeights.Read();
		auto *d_radius = radius.Write();
		auto *d_z = z.Write();
		auto *d_weight = weight.Write();
		auto *d_dofs = elementDofs.Write();
		int const nq = quadPoints;
		mfem::forall( elements*nq, [ = ] MFEM_HOST_DEVICE ( int k )
		{
			int const e = k/nq;
			int const q = k - e*nq;
			d_radius[ k ] = d_x[ q + nq*( 0 + 2*e ) ];
			d_z[ k ] = d_x[ q + nq*( 1 + 2*e ) ];
			d_weight[ k ] = d_w[ q ]*d_detJ[ q + nq*e ];
		} );
		mfem::forall( elements*dofs, [ = ] MFEM_HOST_DEVICE ( int k )
		{
			d_dofs[ k ] = k;
		} );
		return true;
	}

	bool setEntryKernel( mfem::Vector &v, int index, mfem::real_t value,
	                     bool add )
	{
		auto *d_v = v.ReadWrite();
		mfem::forall( 1, [ = ] MFEM_HOST_DEVICE ( int )
		{
			d_v[ index ] = add ? d_v[ index ] + value : value;
		} );
		return true;
	}

	bool pointPotentialKernel( mfem::Vector const &shape,
	                           mfem::Array<int> const &elementDofs,
	                           mfem::Vector const &potential, int potentialStart,
	                           int elements, int quadPoints, int dofs,
	                           mfem::Vector &psi )
	{
		auto const *d_shape = shape.Read();
		auto const *d_dofs = elementDofs.Read();
		auto const *d_potential = potential.Read();
		auto *d_psi = psi.Write();
		int const nq = quadPoints;
		int const nd = dofs;
		int const start = potentialStart;
		mfem::forall( elements*nq, [ = ] MFEM_HOST_DEVICE ( int k )
		{
			int const e = k/nq;
			int const q = k - e*nq;
			mfem::real_t sum = mfem::real_t( 0 );
			for ( int j = 0; j < nd; ++j )
				sum += d_shape[ q*nd + j ]*d_potential[ start + d_dofs[ e*nd + j ] ];
			d_psi[ k ] = sum;
		} );
		return true;
	}

	bool axisSourceScanKernel( SourceView const &view,
	                           mfem::Array<KnotView> const &knots,
	                           mfem::Array<CoilView> const &coils, int coilCount,
	                           mfem::real_t permeability,
	                           mfem::Array<AffineTriangle> const &maps,
	                           mfem::Vector const &nodeXi,
	                           mfem::Vector const &nodeEta,
	                           mfem::Vector const &potential,
	                           int elements, int nodes,
	                           mfem::real_t &scale, mfem::Array<int> &onAxis )
	{
		auto const *d_knots = knots.Read();
		auto const *d_coils = coils.Read();
		auto const *d_maps = maps.Read();
		auto const *d_xi = nodeXi.Read();
		auto const *d_eta = nodeEta.Read();
		auto const *d_psi = potential.Read();
		auto *d_axis = onAxis.Write();
		SourceView const source = view;
		int const count = coilCount;
		mfem::real_t const mu = permeability;
		int const nn = nodes;

		mfem::forall( elements, [ = ] MFEM_HOST_DEVICE ( int e )
		{
			AffineTriangle const map = d_maps[ e ];
			int any = 0;
			for ( int i = 0; i < nn; ++i )
			{
				mfem::real_t radius, z;
				affineTransform( map, d_xi[ i ], d_eta[ i ], radius, z );
				if ( radius == mfem::real_t( 0 ) )
					any = 1;
			}
			d_axis[ e ] = any;
		} );

		mfem::real_t peak = mfem::real_t( 0 );
		mfem::Array<mfem::real_t> work;
		mfem::reduce( elements*nn, peak,
		              [ = ] MFEM_HOST_DEVICE ( int k, mfem::real_t &r )
		              {
		                      int const e = k/nn;
		                      int const i = k - e*nn;
		                      mfem::real_t radius, z;
		                      affineTransform( d_maps[ e ], d_xi[ i ], d_eta[ i ],
		                                       radius, z );
		                      mfem::real_t const value =
		                              sourceF( source, d_knots, radius, z, d_psi[ k ] )
		                              + coilF( d_coils, count, mu, radius, z );
		                      r = fmax( r, fabs( value ) );
		              },
		              mfem::MaxReducer<mfem::real_t>{}, true, work );
		scale = peak;
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

	bool nodalPeakKernel( mfem::Vector const &, int, int, mfem::real_t &,
	                      int & )
	{
		return false;
	}

	bool copyEntriesKernel( mfem::Vector const &, mfem::Array<int> const &,
	                        mfem::Vector & )
	{
		return false;
	}

	bool plasmaLoopGeometryKernel( mfem::Vector const &, mfem::Vector const &,
	                               mfem::Vector const &, int, int, int,
	                               mfem::Vector &, mfem::Vector &,
	                               mfem::Vector &, mfem::Array<int> & )
	{
		return false;
	}

	bool setEntryKernel( mfem::Vector &, int, mfem::real_t, bool )
	{
		return false;
	}

	bool pointPotentialKernel( mfem::Vector const &, mfem::Array<int> const &,
	                           mfem::Vector const &, int, int, int, int,
	                           mfem::Vector & )
	{
		return false;
	}

	bool axisSourceScanKernel( SourceView const &, mfem::Array<KnotView> const &,
	                           mfem::Array<CoilView> const &, int, mfem::real_t,
	                           mfem::Array<AffineTriangle> const &,
	                           mfem::Vector const &, mfem::Vector const &,
	                           mfem::Vector const &, int, int, mfem::real_t &,
	                           mfem::Array<int> & )
	{
		return false;
	}

#endif

}
