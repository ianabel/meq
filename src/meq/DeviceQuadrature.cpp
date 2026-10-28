#include "DeviceQuadrature.hpp"

namespace meq
{

	bool plasmaCurrentPerElementOnDevice( mfem::FiniteElementSpace &fes,
	                                      mfem::Vector const &potential,
	                                      SourceView const &view,
	                                      std::vector<KnotView> const &knots,
	                                      mfem::Array<int> const &inPlasma,
	                                      std::vector<double> const &conductorValues,
	                                      std::vector<int> const &conductorOffsets,
	                                      int extraOrder,
	                                      mfem::Vector &perElement )
	{
		if ( !mfem::Device::Allows( mfem::Backend::DEVICE_MASK ) )
			return false;

		mfem::Mesh *mesh = fes.GetMesh();
		if ( mesh == nullptr )
			return false;

		int const elements = mesh->GetNE();
		if ( elements == 0 )
			return false;
		if ( fes.GetVDim() != 1 )
			return false;
		if ( inPlasma.Size() != elements )
			return false;

		/*
		 * ONE RULE FOR THE MESH, AND IT HAS TO BE THE HOST'S RULE ON EVERY
		 * ELEMENT OR THIS REFUSES.
		 *
		 * sourceRule() is `2*el.GetOrder() + tr.OrderW() + extraOrder` per
		 * element; GetGeometricFactors() takes a single rule. On a
		 * straight-sided mesh OrderW() is 0 everywhere and they agree, which is
		 * every case this is used on today -- but "today" is not a check, so
		 * the order is computed for every element and compared. A mesh that
		 * disagrees keeps the host path rather than being integrated on a rule
		 * it did not ask for.
		 */
		mfem::FiniteElement const *first = fes.GetFE( 0 );
		if ( first == nullptr )
			return false;
		mfem::Geometry::Type const geometry = first->GetGeomType();

		int wantedOrder = -1;
		for ( int e = 0; e < elements; ++e )
		{
			mfem::FiniteElement const *el = fes.GetFE( e );
			if ( el == nullptr || el->GetGeomType() != geometry )
				return false;

			mfem::IsoparametricTransformation transformation;
			mesh->GetElementTransformation( e, &transformation );
			int const order =
				2*el->GetOrder() + transformation.OrderW() + extraOrder;
			if ( e == 0 )
				wantedOrder = order;
			else if ( order != wantedOrder )
				return false;
		}

		mfem::IntegrationRule const &rule =
			mfem::IntRules.Get( geometry, wantedOrder );
		int const quadPoints = rule.GetNPoints();
		int const dofs = first->GetDof();
		int const spaceDimension = mesh->SpaceDimension();

		/*
		 * THE THREE MFEM ABSTRACTIONS THAT MAKE THIS A LOOP RATHER THAN A
		 * KERNEL. Each is device-resident and each is somebody else's to
		 * maintain, which is the whole of why they are used -- CLAUDE.md,
		 * *Prefer a maintained library*.
		 */
		mfem::Operator const *restriction =
			fes.GetElementRestriction( mfem::ElementDofOrdering::NATIVE );
		if ( restriction == nullptr )
			return false;

		// Not Mesh::GetGeometricFactors(), which would give a mesh without
		// nodes a nodal field; see meq::CoordinateField.
		CoordinateField coordinates( *mesh );
		mfem::GeometricFactors const *geometryFactors = &coordinates.factors(
			rule,
			mfem::GeometricFactors::COORDINATES
			| mfem::GeometricFactors::DETERMINANTS );

		mfem::DofToQuad const &maps =
			first->GetDofToQuad( rule, mfem::DofToQuad::FULL );

		mfem::Vector elementPotential( restriction->Height() );
		elementPotential.UseDevice( true );
		restriction->Mult( potential, elementPotential );

		// The rule's weights, as a device vector of its own: IntegrationRule
		// holds IntegrationPoints and a kernel wants the weights alone.
		mfem::Vector quadWeights( quadPoints );
		quadWeights.HostWrite();
		for ( int q = 0; q < quadPoints; ++q )
			quadWeights( q ) = rule.IntPoint( q ).weight;
		quadWeights.UseDevice( true );

		/*
		 * THE KNOTS AS AN mfem::Array SO THE MEMORY MANAGER OWNS THE COPY.
		 * meq::SourceView holds no pointers by design, so the view captures
		 * whole into the lambda; the knots are a separate flat array and have
		 * to be moved deliberately. An empty profile set still needs a
		 * non-empty allocation, because Read() on a zero-sized Array gives
		 * nothing for the kernel to hold -- and sourceF() will not dereference
		 * it for a Soloviev or constant profile.
		 */
		mfem::Array<KnotView> knotArray(
			static_cast<int>( knots.empty() ? 1 : knots.size() ) );
		for ( std::size_t k = 0; k < knots.size(); ++k )
			knotArray[ static_cast<int>( k ) ] = knots[ k ];
		if ( knots.empty() )
			knotArray[ 0 ] = KnotView{ 0, 0, 0 };

		/*
		 * psi_c PER QUADRATURE POINT, FLATTENED AND CHECKED AGAINST THE RULE.
		 *
		 * meq::SourceIntegrator holds it as `values[ offsets[ e ] + q ]`, with
		 * per-element offsets because the rule it was built against could in
		 * principle differ per element. This loop has already refused a
		 * non-uniform rule, so the offsets must be exactly `e*nq` -- and if
		 * they are not, the cache was built against a DIFFERENT rule and its
		 * points are not these points. Refused rather than indexed anyway:
		 * a psi_c read at the wrong quadrature point is a plausible number and
		 * a different equilibrium.
		 *
		 * An EMPTY cache is not a refusal. It is the ordinary case with no
		 * conductor split, the shift is exactly zero, and
		 * `conductorShiftAt()` makes the same promise from the host side.
		 */
		bool const hasShift = !conductorValues.empty();
		mfem::Vector conductorShift( hasShift
			? static_cast<int>( conductorValues.size() ) : 1 );
		if ( hasShift )
		{
			if ( static_cast<int>( conductorOffsets.size() ) != elements )
				return false;
			for ( int e = 0; e < elements; ++e )
				if ( conductorOffsets[ static_cast<std::size_t>( e ) ]
				     != e*quadPoints )
					return false;
			if ( static_cast<int>( conductorValues.size() )
			     != elements*quadPoints )
				return false;

			conductorShift.HostWrite();
			for ( std::size_t k = 0; k < conductorValues.size(); ++k )
				conductorShift( static_cast<int>( k ) ) = conductorValues[ k ];
		}
		else
		{
			conductorShift.HostWrite();
			conductorShift( 0 ) = 0.0;
		}
		conductorShift.UseDevice( true );

		perElement.SetSize( elements );
		perElement.UseDevice( true );

		return plasmaCurrentKernel( elementPotential, quadWeights,
		                            geometryFactors->detJ,
		                            geometryFactors->X,
		                            maps.B, inPlasma,
		                            conductorShift, hasShift,
		                            view, knotArray,
		                            elements, quadPoints, dofs,
		                            spaceDimension, perElement );
	}

}
