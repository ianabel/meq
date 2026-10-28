#include "ResidentIterate.hpp"

#include "Crossings.hpp"

namespace meq
{

	bool DevicePlasmaLoops::prepare( mfem::FiniteElementSpace &potential,
	                                 int extraOrder )
	{
		if ( !mfem::Device::Allows( mfem::Backend::DEVICE_MASK ) )
			return false;

		mfem::Mesh *const meshIn = potential.GetMesh();
		if ( meshIn == nullptr )
			return false;

		if ( status != 0 && space == &potential && mesh == meshIn
		     && sequence == meshIn->GetSequence() && extra == extraOrder
		     && elements == meshIn->GetNE() )
			return status > 0;

		// Refused until proved otherwise, so every `return false` below is also
		// a decision not to ask again for this signature.
		space = &potential;
		mesh = meshIn;
		sequence = meshIn->GetSequence();
		extra = extraOrder;
		elements = meshIn->GetNE();
		status = -1;

		if ( elements == 0 || potential.GetVDim() != 1 )
			return false;
		// THE SCATTER IS RACE-FREE BECAUSE NO TWO ELEMENTS SHARE A DOF, which is
		// a property of L2 and of nothing else MEQ builds. One thread per
		// element writes its own dofs, in the host's order, with no atomics.
		if ( dynamic_cast<mfem::L2_FECollection const *>( potential.FEColl() )
		     == nullptr )
			return false;

		mfem::FiniteElement const *const first = potential.GetFE( 0 );
		if ( first == nullptr )
			return false;

		mfem::IsoparametricTransformation scratch;
		meshIn->GetElementTransformation( 0, &scratch );
		// The rule the host loops take from sourceRule(), and therefore the one
		// meq::SourceIntegrator assembles on: the derivative of the assembled
		// residual has to be taken on the points it was assembled on.
		int const order = 2*first->GetOrder() + scratch.OrderW() + extraOrder;
		mfem::IntegrationRule const &rule =
			mfem::IntRules.Get( first->GetGeomType(), order );
		quadPoints = rule.GetNPoints();
		dofs = first->GetDof();
		if ( quadPoints <= 0 || dofs <= 0 )
			return false;

		/*
		 * PHASE B, AND ATTRIBUTED THERE, which is M-226's point: these are
		 * constants of the discretisation, and an upload counted into phase C
		 * would be read as the iterate crossing.
		 */
		PhaseScope const residency( Phase::Upload );

		shape.SetSize( quadPoints*dofs );
		shape.UseDevice( true );
		{
			mfem::real_t *const host = shape.HostWrite();
			mfem::Vector values( dofs );
			for ( int q = 0; q < quadPoints; ++q )
			{
				first->CalcShape( rule.IntPoint( q ), values );
				for ( int j = 0; j < dofs; ++j )
					host[ q*dofs + j ] = values( j );
			}
		}

		int const points = elements*quadPoints;
		radius.SetSize( points );
		z.SetSize( points );
		weight.SetSize( points );
		elementDofs.SetSize( elements*dofs );
		radius.UseDevice( true );
		z.UseDevice( true );
		weight.UseDevice( true );
		elementDofs.UseDevice( true );

		/*
		 * ON AN AFFINE MESH OF ONE GEOMETRY THE TABLES ARE BUILT ON THE
		 * DEVICE, out of MFEM's geometric factors -- the coordinates and
		 * determinants meq::DeviceQuadrature already takes for the source
		 * integrand. The host loop below took 97% of the border leg on
		 * fixed-f-diiid at 102,840 triangles, on one core, since a Device
		 * switches MEQ's host threading off.
		 *
		 * Every element then has the same rule, because OrderW() is a constant
		 * of an affine map, and the same finite element, because an
		 * L2_FECollection hands out one per geometry and order. An L2 space
		 * numbers its dofs element by element, which the kernel's table
		 * `e*dofs + j` assumes and three elements are checked against.
		 *
		 * NOT BIT-IDENTICAL to the host loop: MFEM's factors interpolate the
		 * vertices where Transform() applies the affine map, so a radius or a
		 * weight can differ in the last bit.
		 */
		mfem::FiniteElementSpace const *nodal = meshIn->GetNodalFESpace();
		bool const affine = nodal == nullptr || nodal->GetMaxElementOrder() == 1;
		bool const uniform =
			meshIn->GetNumGeometries( meshIn->Dimension() ) == 1
			&& !potential.IsVariableOrder()
			&& meshIn->SpaceDimension() == 2;
		if ( affine && uniform )
		{
			bool numbered = true;
			mfem::Array<int> local;
			int const probes[ 3 ] = { 0, elements/2, elements - 1 };
			for ( int const e : probes )
			{
				potential.GetElementDofs( e, local );
				if ( local.Size() != dofs )
					numbered = false;
				for ( int j = 0; numbered && j < dofs; ++j )
					numbered = local[ j ] == e*dofs + j;
			}

			mfem::GeometricFactors const *factors = numbered
				? meshIn->GetGeometricFactors(
					rule, mfem::GeometricFactors::COORDINATES
					      | mfem::GeometricFactors::DETERMINANTS )
				: nullptr;
			if ( factors != nullptr )
			{
				mfem::Vector quadWeights( quadPoints );
				{
					mfem::real_t *const host = quadWeights.HostWrite();
					for ( int q = 0; q < quadPoints; ++q )
						host[ q ] = rule.IntPoint( q ).weight;
				}
				quadWeights.UseDevice( true );
				if ( plasmaLoopGeometryKernel( factors->X, factors->detJ,
				                               quadWeights, elements, quadPoints,
				                               dofs, radius, z, weight,
				                               elementDofs ) )
				{
					// The mesh's nodes back on the host, for every element
					// transformation after this; see affineTrianglesOf().
					if ( mfem::GridFunction const *meshNodes = meshIn->GetNodes() )
						meshNodes->HostRead();
					countedDeviceRead( shape, __LINE__ );
					status = 1;
					return true;
				}
			}
		}

		{
			mfem::real_t *const hostRadius = radius.HostWrite();
			mfem::real_t *const hostZ = z.HostWrite();
			mfem::real_t *const hostWeight = weight.HostWrite();
			int *const hostDofs = elementDofs.HostWrite();

			mfem::Array<int> local;
			mfem::DofTransformation doftrans;
			mfem::Vector point;
			for ( int e = 0; e < elements; ++e )
			{
				// ONE FINITE ELEMENT FOR THE WHOLE MESH, and checked by
				// identity: L2_FECollection hands out one object per geometry
				// and order, so a different pointer is a different element and
				// `shape` above would be the wrong table for it.
				if ( potential.GetFE( e ) != first )
					return false;

				meshIn->GetElementTransformation( e, &scratch );
				if ( 2*first->GetOrder() + scratch.OrderW() + extraOrder != order )
					return false;

				potential.GetElementDofs( e, local, doftrans );
				if ( local.Size() != dofs || !doftrans.IsIdentity() )
					return false;
				for ( int j = 0; j < dofs; ++j )
				{
					if ( local[ j ] < 0 )
						return false;
					hostDofs[ e*dofs + j ] = local[ j ];
				}

				/*
				 * THE HOST LOOPS' OWN CALLS, IN THEIR ORDER: SetIntPoint, then
				 * Transform for the point and Weight for the measure. The
				 * product `ip.weight*Weight()` is formed here exactly as every
				 * host loop forms it, so the device multiplies the same double.
				 */
				for ( int q = 0; q < quadPoints; ++q )
				{
					mfem::IntegrationPoint const &ip = rule.IntPoint( q );
					scratch.SetIntPoint( &ip );
					scratch.Transform( ip, point );
					hostRadius[ e*quadPoints + q ] = point( 0 );
					hostZ[ e*quadPoints + q ] = point( 1 );
					hostWeight[ e*quadPoints + q ] = ip.weight*scratch.Weight();
				}
			}
		}

		countedDeviceRead( shape, __LINE__ );
		countedDeviceRead( radius, __LINE__ );
		countedDeviceRead( z, __LINE__ );
		countedDeviceRead( weight, __LINE__ );
		countedDeviceRead( elementDofs, __LINE__ );

		status = 1;
		return true;
	}

	PlasmaLoopView DevicePlasmaLoops::view( PlasmaQuantity quantity,
	                                        SourceView const &source,
	                                        DeviceConstants const &constants,
	                                        mfem::Vector const &state,
	                                        int potentialStart,
	                                        mfem::real_t scale ) const
	{
		PlasmaLoopView inputs{};
		inputs.source = source;
		inputs.knots = constants.knots.Read();
		inputs.flags = constants.flags.Read();
		inputs.shiftOffsets = constants.offsetsOfShift.Read();
		inputs.shift = constants.shift.Read();
		inputs.hasShift = constants.shiftCount > 0 ? 1 : 0;
		inputs.shape = shape.Read();
		inputs.radius = radius.Read();
		inputs.z = z.Read();
		inputs.weight = weight.Read();
		inputs.elementDofs = elementDofs.Read();
		// The one buffer here that is not a constant: counted, because under a
		// forced route on a host-valid iterate this is an upload, and that is
		// exactly what the measurement is for.
		inputs.state = countedDeviceRead( state, __LINE__ );
		inputs.potentialStart = potentialStart;
		inputs.elements = elements;
		inputs.quadPoints = quadPoints;
		inputs.dofs = dofs;
		inputs.scale = scale;
		inputs.quantity = static_cast<int>( quantity );
		return inputs;
	}

	bool DevicePlasmaLoops::integrate( PlasmaQuantity quantity,
	                                   SourceView const &source,
	                                   DeviceConstants const &constants,
	                                   mfem::Vector const &state,
	                                   int potentialStart,
	                                   mfem::Vector &first,
	                                   mfem::Vector *second ) const
	{
		if ( status <= 0 )
			return false;
		bool const two = quantity == PlasmaQuantity::NormalisationCorner;
		if ( quantity != PlasmaQuantity::Current && !two )
			return false;
		if ( two && second == nullptr )
			return false;
		// A table built for a different element count is a table of other
		// elements; DeviceConstants refreshes by comparison, so this is the
		// check that the two were built for one mesh.
		if ( constants.elementCount != elements )
			return false;

		first.SetSize( elements );
		first.UseDevice( true );
		mfem::real_t *const a = first.Write();
		mfem::real_t *b = nullptr;
		if ( two )
		{
			second->SetSize( elements );
			second->UseDevice( true );
			b = second->Write();
		}
		return plasmaElementKernel(
			view( quantity, source, constants, state, potentialStart, 1.0 ), a, b );
	}

	bool DevicePlasmaLoops::scatter( PlasmaQuantity quantity,
	                                 SourceView const &source,
	                                 DeviceConstants const &constants,
	                                 mfem::Vector const &state,
	                                 int potentialStart, mfem::real_t scale,
	                                 mfem::Vector &out ) const
	{
		if ( status <= 0 )
			return false;
		if ( quantity == PlasmaQuantity::Current
		     || quantity == PlasmaQuantity::NormalisationCorner )
			return false;
		if ( constants.elementCount != elements )
			return false;

		out.SetSize( state.Size() );
		out.UseDevice( true );
		// Device-aware, because UseDevice() is now true: a forall rather than a
		// host loop, so the zero is written where the scatter will add to it.
		out = 0.0;
		return plasmaScatterKernel(
			view( quantity, source, constants, state, potentialStart, scale ),
			out.ReadWrite() );
	}

	bool DevicePlasmaLoops::pointValues( mfem::Vector const &potential,
	                                     int potentialStart, mfem::Vector &psi,
	                                     mfem::Vector const *&radiusOut,
	                                     mfem::Vector const *&zOut,
	                                     mfem::Vector const *&weightOut,
	                                     int &pointsPerElement ) const
	{
		if ( status <= 0 )
			return false;
		psi.SetSize( elements*quadPoints );
		psi.UseDevice( true );
		if ( !pointPotentialKernel( shape, elementDofs, potential, potentialStart,
		                            elements, quadPoints, dofs, psi ) )
			return false;
		radiusOut = &radius;
		zOut = &z;
		weightOut = &weight;
		pointsPerElement = quadPoints;
		return true;
	}

	bool axisSourceScanOnDevice( Source const &source,
	                             mfem::FiniteElementSpace &space,
	                             mfem::Vector const &potential,
	                             mfem::Array<AffineTriangle> const &maps,
	                             mfem::real_t &scale, mfem::Array<int> &onAxis )
	{
		if ( !mfem::Device::Allows( mfem::Backend::DEVICE_MASK ) )
			return false;

		SourceView view{};
		std::vector<KnotView> knotList;
		if ( !makePlasmaSourceView( source, view, knotList ) )
			return false;
		std::vector<CoilView> coilList;
		mfem::real_t permeability = 0.0;
		if ( !makeCoilViews( source, coilList, permeability ) )
			return false;

		int const elements = space.GetNE();
		if ( elements == 0 || space.GetVDim() != 1 || space.IsVariableOrder()
		     || maps.Size() != elements
		     || dynamic_cast<mfem::L2_FECollection const *>( space.FEColl() ) == nullptr )
			return false;
		mfem::FiniteElement const *const first = space.GetFE( 0 );
		mfem::Mesh *const mesh = space.GetMesh();
		if ( first == nullptr || mesh == nullptr
		     || mesh->GetNumGeometries( mesh->Dimension() ) != 1 )
			return false;
		mfem::IntegrationRule const &nodeRule = first->GetNodes();
		int const nodes = first->GetDof();
		if ( nodeRule.GetNPoints() != nodes )
			return false;

		// Element by element, which the kernel's `e*nodes + i` assumes; three
		// elements are checked against it, as DevicePlasmaLoops::prepare() does.
		mfem::Array<int> local;
		int const probes[ 3 ] = { 0, elements/2, elements - 1 };
		for ( int const e : probes )
		{
			space.GetElementDofs( e, local );
			if ( local.Size() != nodes )
				return false;
			for ( int i = 0; i < nodes; ++i )
				if ( local[ i ] != e*nodes + i )
					return false;
		}

		mfem::Vector nodeXi( nodes ), nodeEta( nodes );
		for ( int i = 0; i < nodes; ++i )
		{
			nodeXi( i ) = nodeRule.IntPoint( i ).x;
			nodeEta( i ) = nodeRule.IntPoint( i ).y;
		}
		nodeXi.UseDevice( true );
		nodeEta.UseDevice( true );

		// Never zero-sized: Read() on an empty Array hands a kernel nothing to
		// capture, and the counts say how much of each is real.
		mfem::Array<KnotView> knots( std::max<int>( 1, static_cast<int>( knotList.size() ) ) );
		for ( std::size_t k = 0; k < knotList.size(); ++k )
			knots[ static_cast<int>( k ) ] = knotList[ k ];
		knots.UseDevice( true );
		mfem::Array<CoilView> coils( std::max<int>( 1, static_cast<int>( coilList.size() ) ) );
		for ( std::size_t k = 0; k < coilList.size(); ++k )
			coils[ static_cast<int>( k ) ] = coilList[ k ];
		coils.UseDevice( true );

		onAxis.SetSize( elements );
		onAxis.UseDevice( true );
		if ( !axisSourceScanKernel( view, knots, coils, static_cast<int>( coilList.size() ),
		                            permeability, maps, nodeXi, nodeEta, potential,
		                            elements, nodes, scale, onAxis ) )
			return false;
		onAxis.HostRead();
		return true;
	}

	bool setEntryOnDevice( mfem::Vector &v, int index, mfem::real_t value,
	                       bool add )
	{
		if ( !mfem::Device::Allows( mfem::Backend::DEVICE_MASK ) )
			return false;
		if ( index < 0 || index >= v.Size() )
			return false;
		v.UseDevice( true );
		return setEntryKernel( v, index, value, add );
	}

	bool gatherOnDevice( mfem::Vector const &source,
	                     mfem::Array<int> const &indices, mfem::Vector &out )
	{
		if ( !mfem::Device::Allows( mfem::Backend::DEVICE_MASK ) )
			return false;
		out.SetSize( indices.Size() );
		out.UseDevice( true );
		if ( indices.Size() == 0 )
			return true;
		// Counted, so an index list re-sent every call shows as the upload it
		// is; the kernel's own Read()s then move nothing.
		countedDeviceRead( source, __LINE__ );
		countedDeviceRead( indices, __LINE__ );
		return gatherKernel( source, indices, out );
	}

	bool orderedDotOnDevice( mfem::Vector const &source,
	                         mfem::Array<int> const &indices,
	                         mfem::Vector const &weights, mfem::Vector &out )
	{
		if ( !mfem::Device::Allows( mfem::Backend::DEVICE_MASK ) )
			return false;
		if ( weights.Size() != indices.Size() )
			return false;
		out.SetSize( 1 );
		out.UseDevice( true );
		countedDeviceRead( source, __LINE__ );
		countedDeviceRead( indices, __LINE__ );
		countedDeviceRead( weights, __LINE__ );
		return orderedDotKernel( source, indices, weights, out );
	}

	bool nodalPeakOnDevice( mfem::Vector const &source, int begin, int end,
	                        mfem::real_t &best, int &bestIndex )
	{
		if ( !mfem::Device::Allows( mfem::Backend::DEVICE_MASK ) )
			return false;
		if ( begin < 0 || end > source.Size() || begin > end )
			return false;
		countedDeviceRead( source, __LINE__ );
		return nodalPeakKernel( source, begin, end, best, bestIndex );
	}

	bool copyEntriesOnDevice( mfem::Vector const &source,
	                          mfem::Array<int> const &indices,
	                          mfem::Vector &destination )
	{
		if ( !mfem::Device::Allows( mfem::Backend::DEVICE_MASK ) )
			return false;
		if ( source.Size() != destination.Size() )
			return false;
		if ( indices.Size() == 0 )
			return true;
		countedDeviceRead( source, __LINE__ );
		countedDeviceRead( indices, __LINE__ );
		return copyEntriesKernel( source, indices, destination );
	}

}
