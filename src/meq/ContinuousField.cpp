#include "ContinuousField.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

#include "Threading.hpp"

/*
 * The header has the account: what the two fields are, why their degrees,
 * why psihat on the interior edges, and why the build is the serial one bit
 * for bit. This file is the arithmetic.
 *
 * WHAT THE THREADED LOOPS RELY ON FROM MFEM, and it was read rather than
 * assumed: FiniteElementSpace::GetElementDofs(), GetEdgeInteriorDofs() and
 * GetElementInteriorDofs() read the mesh's own tables and build nothing
 * lazily on a fixed-order space; GridFunction::GetValue() and
 * GetVectorValue() touch the space's shared element transformation only for
 * a map type other than VALUE, which buildContinuousFields() refuses; and the
 * face and element transformations are the caller-allocated overloads into
 * thread_local objects (CLAUDE.md, MFEM's API).
 */

namespace meq
{

	namespace
	{
		/// P_m on [-1, 1], by the three-term recurrence.
		double legendre( int m, double x )
		{
			double p0 = 1.0;
			double p1 = x;
			if ( m == 0 )
				return p0;
			for ( int n = 1; n < m; ++n )
			{
				double const p2 = ( ( 2.0*n + 1.0 )*x*p1 - n*p0 )/( n + 1.0 );
				p0 = p1;
				p1 = p2;
			}
			return p1;
		}

		/// The cell moments' test functions at a reference point: monomials
		/// of total degree at most p - 3 on a triangle, of degree at most
		/// p - 2 in each variable on a quadrilateral -- as many as the
		/// element's interior dofs either way.
		void cellBasis( mfem::Geometry::Type geometry, int p, mfem::IntegrationPoint const &ip,
		                std::vector<double> &basis )
		{
			basis.clear();
			if ( geometry == mfem::Geometry::TRIANGLE )
			{
				for ( int total = 0; total <= p - 3; ++total )
					for ( int ax = total; ax >= 0; --ax )
						basis.push_back( std::pow( ip.x, ax )*std::pow( ip.y, total - ax ) );
				return;
			}
			for ( int ay = 0; ay <= p - 2; ++ay )
				for ( int ax = 0; ax <= p - 2; ++ax )
					basis.push_back( std::pow( ip.x, ax )*std::pow( ip.y, ay ) );
		}

		/// Solve the small moment system for the unknown dofs and write them.
		void writeSolution( mfem::DenseMatrix &a, mfem::Vector const &b, mfem::Array<int> const &dofs,
		                    mfem::real_t *data )
		{
			mfem::DenseMatrixInverse const inverse( a );
			mfem::Vector x( b.Size() );
			inverse.Mult( b, x );
			for ( int j = 0; j < dofs.Size(); ++j )
				data[ dofs[ j ] ] = x( j );
		}
	}

	void momentInterpolant( mfem::GridFunction &out, CellValue const &cell, FaceValue const &face )
	{
		mfem::FiniteElementSpace const &fes = *out.FESpace();
		mfem::Mesh const &mesh = *fes.GetMesh();
		if ( fes.GetVDim() != 1 )
			throw std::invalid_argument( "meq::momentInterpolant: the space must be scalar" );
		for ( int e = 0; e < mesh.GetNE(); ++e )
		{
			mfem::Geometry::Type const geometry = mesh.GetElementBaseGeometry( e );
			if ( geometry != mfem::Geometry::TRIANGLE && geometry != mfem::Geometry::SQUARE )
				throw std::invalid_argument( "meq::momentInterpolant: element " + std::to_string( e )
				                             + " is neither a triangle nor a quadrilateral" );
		}
		int const p = fes.GetElementOrder( 0 );
		int const elements = mesh.GetNE();

		// ZEROED THROUGH THE HOST POINTER. `out = 0.0` is device-aware, and a
		// GridFunction is born UseDevice( true ): under a Device it zeroes the
		// device copy and invalidates the host one, so every raw write below
		// lands in a page MFEM believes stale, and the next HostRead() copies
		// the device's zeros over it (and the debug backend faults on the
		// first raw read).
		mfem::real_t *const data = out.HostWrite();
		std::fill_n( data, out.Size(), mfem::real_t( 0 ) );

		// VERTICES: each element's value at its corners, threaded, then the
		// means summed serially in element order.
		std::vector<std::vector<double>> corners( static_cast<std::size_t>( elements ) );
		MEQ_OMP_ELEMENT_LOOP
		for ( int e = 0; e < elements; ++e )
		{
			mfem::IntegrationRule const *vertices = mfem::Geometries.GetVertices( mesh.GetElementBaseGeometry( e ) );
			std::vector<double> &values = corners[ static_cast<std::size_t>( e ) ];
			values.resize( static_cast<std::size_t>( vertices->GetNPoints() ) );
			for ( int i = 0; i < vertices->GetNPoints(); ++i )
				values[ static_cast<std::size_t>( i ) ] = cell( e, vertices->IntPoint( i ) );
		}
		std::vector<double> sum( static_cast<std::size_t>( mesh.GetNV() ), 0.0 );
		std::vector<int> count( sum.size(), 0 );
		mfem::Array<int> vertices;
		for ( int e = 0; e < elements; ++e )
		{
			mesh.GetElementVertices( e, vertices );
			for ( int i = 0; i < vertices.Size(); ++i )
			{
				std::size_t const v = static_cast<std::size_t>( vertices[ i ] );
				sum[ v ] += corners[ static_cast<std::size_t>( e ) ][ static_cast<std::size_t>( i ) ];
				++count[ v ];
			}
		}
		mfem::Array<int> vertexDofs;
		for ( int v = 0; v < mesh.GetNV(); ++v )
		{
			fes.GetVertexDofs( v, vertexDofs );
			if ( vertexDofs.Size() == 1 && count[ static_cast<std::size_t>( v ) ] > 0 )
				data[ vertexDofs[ 0 ] ] = sum[ static_cast<std::size_t>( v ) ]/count[ static_cast<std::size_t>( v ) ];
		}

		// EDGES: each edge's interior dofs from @a face's moments against
		// P_(p-2), the vertex dofs at either end already written. In two
		// dimensions face i is edge i.
		//
		// ONLY VERTEX DOFS ARE SUBTRACTED AS KNOWN. The element's other edges'
		// shapes vanish on this edge only to round-off, and their dofs are
		// being written by other threads -- so counting them made the answer
		// depend on the order the edges were done in, threaded or not. H1
		// numbers the vertex dofs first.
		int const vertexDofCount = fes.GetNVDofs();
		int const faces = mesh.GetNumFaces();
		MEQ_OMP_ELEMENT_LOOP
		for ( int f = 0; f < faces; ++f )
		{
			thread_local mfem::FaceElementTransformations faceTransform;
			thread_local mfem::IsoparametricTransformation one;
			thread_local mfem::IsoparametricTransformation two;
			mfem::Array<int> interior;
			fes.GetEdgeInteriorDofs( f, interior );
			int const n = interior.Size();
			if ( n == 0 )
				continue;
			mesh.GetFaceElementTransformations( f, faceTransform, one, two );
			int const element = faceTransform.Elem1No;
			mfem::FiniteElement const *fe = fes.GetFE( element );
			mfem::Array<int> local;
			fes.GetElementDofs( element, local );
			std::vector<int> unknown( static_cast<std::size_t>( local.Size() ) );
			for ( int i = 0; i < local.Size(); ++i )
				unknown[ static_cast<std::size_t>( i ) ] = interior.Find( local[ i ] );
			mfem::Vector shape( fe->GetDof() );
			mfem::DenseMatrix a( n, n );
			mfem::Vector b( n );
			a = 0.0;
			b = 0.0;
			mfem::IntegrationRule const &rule = mfem::IntRules.Get( mfem::Geometry::SEGMENT, 2*p + 2 );
			for ( int q = 0; q < rule.GetNPoints(); ++q )
			{
				mfem::IntegrationPoint const &ip = rule.IntPoint( q );
				faceTransform.SetAllIntPoints( &ip );
				fe->CalcShape( faceTransform.GetElement1IntPoint(), shape );
				double const weight = ip.weight*faceTransform.Weight();
				double known = 0.0;
				for ( int i = 0; i < local.Size(); ++i )
					if ( unknown[ static_cast<std::size_t>( i ) ] < 0 && local[ i ] < vertexDofCount )
						known += shape( i )*data[ local[ i ] ];
				double const target = face( f, faceTransform, ip ) - known;
				for ( int m = 0; m < n; ++m )
				{
					double const basis = legendre( m, 2.0*ip.x - 1.0 )*weight;
					b( m ) += basis*target;
					for ( int i = 0; i < local.Size(); ++i )
						if ( unknown[ static_cast<std::size_t>( i ) ] >= 0 )
							a( m, unknown[ static_cast<std::size_t>( i ) ] ) += basis*shape( i );
				}
			}
			writeSolution( a, b, interior, data );
		}

		// CELLS: each element's interior dofs from @a cell's moments, every
		// boundary dof already written.
		MEQ_OMP_ELEMENT_LOOP
		for ( int e = 0; e < elements; ++e )
		{
			thread_local mfem::IsoparametricTransformation transform;
			mfem::Array<int> interior;
			fes.GetElementInteriorDofs( e, interior );
			int const n = interior.Size();
			if ( n == 0 )
				continue;
			mesh.GetElementTransformation( e, &transform );
			mfem::FiniteElement const *fe = fes.GetFE( e );
			mfem::Array<int> local;
			fes.GetElementDofs( e, local );
			std::vector<int> unknown( static_cast<std::size_t>( local.Size() ) );
			for ( int i = 0; i < local.Size(); ++i )
				unknown[ static_cast<std::size_t>( i ) ] = interior.Find( local[ i ] );
			mfem::Vector shape( fe->GetDof() );
			mfem::DenseMatrix a( n, n );
			mfem::Vector b( n );
			a = 0.0;
			b = 0.0;
			std::vector<double> basis;
			mfem::IntegrationRule const &rule = mfem::IntRules.Get( fe->GetGeomType(), 2*p + 2 );
			for ( int q = 0; q < rule.GetNPoints(); ++q )
			{
				mfem::IntegrationPoint const &ip = rule.IntPoint( q );
				transform.SetIntPoint( &ip );
				fe->CalcShape( ip, shape );
				double const weight = ip.weight*transform.Weight();
				double known = 0.0;
				for ( int i = 0; i < local.Size(); ++i )
					if ( unknown[ static_cast<std::size_t>( i ) ] < 0 )
						known += shape( i )*data[ local[ i ] ];
				double const target = cell( e, ip ) - known;
				cellBasis( fe->GetGeomType(), p, ip, basis );
				for ( int m = 0; m < n; ++m )
				{
					double const test = basis[ static_cast<std::size_t>( m ) ]*weight;
					b( m ) += test*target;
					for ( int i = 0; i < local.Size(); ++i )
						if ( unknown[ static_cast<std::size_t>( i ) ] >= 0 )
							a( m, unknown[ static_cast<std::size_t>( i ) ] ) += test*shape( i );
				}
			}
			writeSolution( a, b, interior, data );
		}
	}

	ContinuousFields buildContinuousFields( mfem::GridFunction const &potential,
	                                        mfem::GridFunction const &flux,
	                                        mfem::GridFunction const &trace )
	{
		mfem::FiniteElementSpace const &potentialSpace = *potential.FESpace();
		mfem::FiniteElementSpace const &fluxSpace = *flux.FESpace();
		mfem::FiniteElementSpace const &traceSpace = *trace.FESpace();
		mfem::Mesh *mesh = potentialSpace.GetMesh();
		if ( mesh->Dimension() != 2 || fluxSpace.GetVDim() != 2 || fluxSpace.GetMesh() != mesh )
			throw std::invalid_argument( "meq::buildContinuousFields: a two-dimensional psi* and a vdim-2 "
			                             "q* on one mesh are needed" );
		// VALUE maps only: GridFunction::GetValue() reaches the space's shared
		// element transformation for any other, and these are evaluated from
		// several threads.
		if ( potentialSpace.GetFE( 0 )->GetMapType() != mfem::FiniteElement::MapType::VALUE
		     || fluxSpace.GetFE( 0 )->GetMapType() != mfem::FiniteElement::MapType::VALUE )
			throw std::invalid_argument( "meq::buildContinuousFields: psi* and q* must map by VALUE" );

		// One host sync each, before any thread reads them.
		static_cast<void>( potential.HostRead() );
		static_cast<void>( flux.HostRead() );
		static_cast<void>( trace.HostRead() );

		int const potentialOrder = potentialSpace.GetElementOrder( 0 );
		int const fluxOrder = fluxSpace.GetElementOrder( 0 );

		ContinuousFields out;
		auto *potentialCollection = new mfem::H1_FECollection( potentialOrder + 1, 2 );
		out.potential = std::make_unique<mfem::GridFunction>( new mfem::FiniteElementSpace( mesh, potentialCollection ) );
		out.potential->MakeOwner( potentialCollection );

		momentInterpolant( *out.potential,
			[ & ]( int e, mfem::IntegrationPoint const &ip ) -> double { return potential.GetValue( e, ip ); },
			[ & ]( int f, mfem::FaceElementTransformations &transform, mfem::IntegrationPoint const &ip ) -> double
			{
				if ( transform.Elem2No < 0 )
					return potential.GetValue( transform.Elem1No, transform.GetElement1IntPoint() );
				// psihat: DG_Interface_FECollection, nodal with VALUE map type,
				// so shape . dofs is the value, in the face's own frame, which
				// is the frame the face transformation's point is in.
				mfem::FiniteElement const *fe = traceSpace.GetFaceElement( f );
				mfem::Array<int> dofs;
				traceSpace.GetFaceVDofs( f, dofs );
				mfem::Vector local;
				trace.GetSubVector( dofs, local );
				mfem::Vector shape( fe->GetDof() );
				fe->CalcShape( ip, shape );
				return shape*local;
			} );

		auto *fluxCollection = new mfem::H1_FECollection( fluxOrder, 2 );
		auto *vectorSpace = new mfem::FiniteElementSpace( mesh, fluxCollection, 2, mfem::Ordering::byNODES );
		out.flux = std::make_unique<mfem::GridFunction>( vectorSpace );
		out.flux->MakeOwner( fluxCollection );
		mfem::FiniteElementSpace scalarSpace( mesh, fluxCollection );
		mfem::real_t *const fluxData = out.flux->HostWrite();
		int const scalarDofs = scalarSpace.GetNDofs();
		for ( int component = 0; component < 2; ++component )
		{
			mfem::GridFunction part( &scalarSpace );
			auto at = [ & ]( int e, mfem::IntegrationPoint const &ip ) -> double
			{
				mfem::Vector value( 2 );
				flux.GetVectorValue( e, ip, value );
				return value( component );
			};
			momentInterpolant( part, at,
				[ & ]( int, mfem::FaceElementTransformations &transform, mfem::IntegrationPoint const & ) -> double
				{
					if ( transform.Elem2No < 0 )
						return at( transform.Elem1No, transform.GetElement1IntPoint() );
					return 0.5*( at( transform.Elem1No, transform.GetElement1IntPoint() )
					             + at( transform.Elem2No, transform.GetElement2IntPoint() ) );
				} );
			mfem::real_t const *const partData = part.HostRead();
			for ( int i = 0; i < scalarDofs; ++i )
				fluxData[ i + component*scalarDofs ] = partData[ i ];
		}
		return out;
	}

}
