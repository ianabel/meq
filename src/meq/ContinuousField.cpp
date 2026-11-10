#include "ContinuousField.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "SourceView.hpp"
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
		// first raw read). MEASUREMENTS.md M-382.
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

	/*
	 * ---------------------------------------------------------------------
	 * THE INTERPOLANT AS BLOCKS. ContinuousField.hpp, ContinuousInterpolant,
	 * has the account; this is the arithmetic.
	 * ---------------------------------------------------------------------
	 */
	namespace
	{
		/// One block before it is packed: the rows it writes, the inputs and
		/// earlier outputs it reads, and its matrix, row-major, in double.
		struct Block
		{
			std::vector<int> out, in, known;
			std::vector<double> matrix;
		};

		LinearBlocks pack( std::vector<Block> const &blocks )
		{
			LinearBlocks packed;
			packed.count = static_cast<int>( blocks.size() );
			std::size_t outs = 0, ins = 0, knowns = 0, entries = 0;
			for ( Block const &b : blocks )
			{
				outs += b.out.size();
				ins += b.in.size();
				knowns += b.known.size();
				entries += b.matrix.size();
			}
			packed.outOffset.SetSize( packed.count + 1 );
			packed.inOffset.SetSize( packed.count + 1 );
			packed.knownOffset.SetSize( packed.count + 1 );
			packed.matrixOffset.SetSize( packed.count + 1 );
			packed.outIndex.SetSize( static_cast<int>( outs ) );
			packed.inIndex.SetSize( static_cast<int>( ins ) );
			packed.knownIndex.SetSize( static_cast<int>( knowns ) );
			packed.matrix.SetSize( static_cast<int>( entries ) );
			// Fresh arrays, so these raw fills are the host copies MFEM vouches
			// for; the kernel's Read() is their one upload.
			int o = 0, i = 0, k = 0, m = 0;
			for ( int b = 0; b < packed.count; ++b )
			{
				Block const &block = blocks[ static_cast<std::size_t>( b ) ];
				packed.outOffset[ b ] = o;
				packed.inOffset[ b ] = i;
				packed.knownOffset[ b ] = k;
				packed.matrixOffset[ b ] = m;
				for ( int v : block.out )
					packed.outIndex[ o++ ] = v;
				for ( int v : block.in )
					packed.inIndex[ i++ ] = v;
				for ( int v : block.known )
					packed.knownIndex[ k++ ] = v;
				for ( double v : block.matrix )
					packed.matrix( m++ ) = static_cast<mfem::real_t>( v );
			}
			packed.outOffset[ packed.count ] = o;
			packed.inOffset[ packed.count ] = i;
			packed.knownOffset[ packed.count ] = k;
			packed.matrixOffset[ packed.count ] = m;
			return packed;
		}

		/// Solve A X = B in place of B, A n x n and B n x m, both row-major,
		/// by Gaussian elimination with partial pivoting in double. The moment
		/// systems are at most ( p - 1 )^2 square.
		void solveInPlace( std::vector<double> a, std::vector<double> &b, int n, int m )
		{
			for ( int c = 0; c < n; ++c )
			{
				int pivot = c;
				for ( int r = c + 1; r < n; ++r )
					if ( std::abs( a[ r*n + c ] ) > std::abs( a[ pivot*n + c ] ) )
						pivot = r;
				if ( a[ pivot*n + c ] == 0.0 )
					throw std::runtime_error( "meq::ContinuousInterpolant: a moment system is singular" );
				if ( pivot != c )
				{
					for ( int j = 0; j < n; ++j )
						std::swap( a[ c*n + j ], a[ pivot*n + j ] );
					for ( int j = 0; j < m; ++j )
						std::swap( b[ c*m + j ], b[ pivot*m + j ] );
				}
				for ( int r = c + 1; r < n; ++r )
				{
					double const factor = a[ r*n + c ]/a[ c*n + c ];
					if ( factor == 0.0 )
						continue;
					for ( int j = c; j < n; ++j )
						a[ r*n + j ] -= factor*a[ c*n + j ];
					for ( int j = 0; j < m; ++j )
						b[ r*m + j ] -= factor*b[ c*m + j ];
				}
			}
			for ( int c = n - 1; c >= 0; --c )
				for ( int j = 0; j < m; ++j )
				{
					double sum = b[ c*m + j ];
					for ( int k = c + 1; k < n; ++k )
						sum -= a[ c*n + k ]*b[ k*m + j ];
					b[ c*m + j ] = sum/a[ c*n + c ];
				}
		}

		/// An input's value at a point, as a list of coefficients on input
		/// entries: VALUE-mapped L2 shapes against one component's vdofs.
		void addElementInput( mfem::FiniteElementSpace const &space, int component, int element,
		                      mfem::IntegrationPoint const &ip, double scale,
		                      std::vector<std::pair<int, double>> &terms )
		{
			mfem::FiniteElement const *fe = space.GetFE( element );
			mfem::Array<int> vdofs;
			space.GetElementVDofs( element, vdofs );
			int const dof = fe->GetDof();
			mfem::Vector shape( dof );
			fe->CalcShape( ip, shape );
			for ( int j = 0; j < dof; ++j )
				terms.emplace_back( vdofs[ j + component*dof ], scale*shape( j ) );
		}

		/// The column of each distinct input entry, in first-seen order.
		int columnOf( std::vector<int> &columns, int index )
		{
			for ( std::size_t c = 0; c < columns.size(); ++c )
				if ( columns[ c ] == index )
					return static_cast<int>( c );
			columns.push_back( index );
			return static_cast<int>( columns.size() ) - 1;
		}

		/*
		 * THE THREE STAGES FOR ONE SCALAR OUTPUT. `cellInput` gives the field
		 * to interpolate inside element e at ip; `faceInput` its value on face
		 * f at the face's point, the face transformation set there. Both
		 * append coefficients on input entries. The quadratures are
		 * momentInterpolant()'s, rule for rule; the moments go into the
		 * matrix instead of into one solution.
		 */
		using CellInput = std::function<void( int, mfem::IntegrationPoint const &,
		                                      std::vector<std::pair<int, double>> & )>;
		using FaceInput = std::function<void( int, mfem::FaceElementTransformations &,
		                                      mfem::IntegrationPoint const &,
		                                      std::vector<std::pair<int, double>> & )>;

		void buildStages( mfem::FiniteElementSpace const &out, CellInput const &cellInput,
		                  FaceInput const &faceInput, LinearBlocks stages[ 3 ] )
		{
			mfem::Mesh const &mesh = *out.GetMesh();
			int const p = out.GetElementOrder( 0 );
			int const elements = mesh.GetNE();

			// VERTICES: one row per vertex, the mean of each meeting element's
			// value there, its terms in element order.
			{
				std::vector<std::vector<std::pair<int, double>>> corner( static_cast<std::size_t>( elements ) );
				std::vector<std::vector<int>> cornerStart( static_cast<std::size_t>( elements ) );
				MEQ_OMP_ELEMENT_LOOP
				for ( int e = 0; e < elements; ++e )
				{
					mfem::IntegrationRule const *vertices = mfem::Geometries.GetVertices( mesh.GetElementBaseGeometry( e ) );
					auto &terms = corner[ static_cast<std::size_t>( e ) ];
					auto &start = cornerStart[ static_cast<std::size_t>( e ) ];
					for ( int i = 0; i < vertices->GetNPoints(); ++i )
					{
						start.push_back( static_cast<int>( terms.size() ) );
						cellInput( e, vertices->IntPoint( i ), terms );
					}
					start.push_back( static_cast<int>( terms.size() ) );
				}
				std::vector<std::vector<std::pair<int, int>>> meeting( static_cast<std::size_t>( mesh.GetNV() ) );
				mfem::Array<int> vertices;
				for ( int e = 0; e < elements; ++e )
				{
					mesh.GetElementVertices( e, vertices );
					for ( int i = 0; i < vertices.Size(); ++i )
						meeting[ static_cast<std::size_t>( vertices[ i ] ) ].emplace_back( e, i );
				}
				std::vector<Block> blocks;
				mfem::Array<int> vertexDofs;
				for ( int v = 0; v < mesh.GetNV(); ++v )
				{
					auto const &list = meeting[ static_cast<std::size_t>( v ) ];
					out.GetVertexDofs( v, vertexDofs );
					if ( vertexDofs.Size() != 1 || list.empty() )
						continue;
					Block block;
					block.out.push_back( vertexDofs[ 0 ] );
					double const share = 1.0/static_cast<double>( list.size() );
					for ( auto const &[ e, i ] : list )
					{
						auto const &terms = corner[ static_cast<std::size_t>( e ) ];
						auto const &start = cornerStart[ static_cast<std::size_t>( e ) ];
						for ( int t = start[ static_cast<std::size_t>( i ) ]; t < start[ static_cast<std::size_t>( i ) + 1 ]; ++t )
						{
							block.in.push_back( terms[ static_cast<std::size_t>( t ) ].first );
							block.matrix.push_back( share*terms[ static_cast<std::size_t>( t ) ].second );
						}
					}
					blocks.push_back( std::move( block ) );
				}
				stages[ 0 ] = pack( blocks );
			}

			// EDGES: the interior dofs from the face values' moments against
			// P_(p-2), the vertex dofs known. In two dimensions face i is edge i.
			int const vertexDofCount = out.GetNVDofs();
			int const faces = mesh.GetNumFaces();
			std::vector<Block> edgeBlocks( static_cast<std::size_t>( faces ) );
			MEQ_OMP_ELEMENT_LOOP
			for ( int f = 0; f < faces; ++f )
			{
				thread_local mfem::FaceElementTransformations faceTransform;
				thread_local mfem::IsoparametricTransformation one;
				thread_local mfem::IsoparametricTransformation two;
				mfem::Array<int> interior;
				out.GetEdgeInteriorDofs( f, interior );
				int const n = interior.Size();
				if ( n == 0 )
					continue;
				mesh.GetFaceElementTransformations( f, faceTransform, one, two );
				int const element = faceTransform.Elem1No;
				mfem::FiniteElement const *fe = out.GetFE( element );
				mfem::Array<int> local;
				out.GetElementDofs( element, local );
				std::vector<int> unknown( static_cast<std::size_t>( local.Size() ) );
				for ( int i = 0; i < local.Size(); ++i )
					unknown[ static_cast<std::size_t>( i ) ] = interior.Find( local[ i ] );
				Block &block = edgeBlocks[ static_cast<std::size_t>( f ) ];
				for ( int i = 0; i < local.Size(); ++i )
					if ( unknown[ static_cast<std::size_t>( i ) ] < 0 && local[ i ] < vertexDofCount )
						block.known.push_back( local[ i ] );
				std::vector<double> a( static_cast<std::size_t>( n*n ), 0.0 );
				std::vector<std::vector<double>> inRows( static_cast<std::size_t>( n ) );
				std::vector<std::vector<double>> knownRows( static_cast<std::size_t>( n ),
					std::vector<double>( block.known.size(), 0.0 ) );
				mfem::Vector shape( fe->GetDof() );
				std::vector<std::pair<int, double>> terms;
				mfem::IntegrationRule const &rule = mfem::IntRules.Get( mfem::Geometry::SEGMENT, 2*p + 2 );
				for ( int q = 0; q < rule.GetNPoints(); ++q )
				{
					mfem::IntegrationPoint const &ip = rule.IntPoint( q );
					faceTransform.SetAllIntPoints( &ip );
					fe->CalcShape( faceTransform.GetElement1IntPoint(), shape );
					double const weight = ip.weight*faceTransform.Weight();
					terms.clear();
					faceInput( f, faceTransform, ip, terms );
					for ( int m = 0; m < n; ++m )
					{
						double const basis = legendre( m, 2.0*ip.x - 1.0 )*weight;
						for ( auto const &[ index, value ] : terms )
						{
							int const c = columnOf( block.in, index );
							auto &row = inRows[ static_cast<std::size_t>( m ) ];
							if ( row.size() < block.in.size() )
								row.resize( block.in.size(), 0.0 );
							row[ static_cast<std::size_t>( c ) ] += basis*value;
						}
						int kc = 0;
						for ( int i = 0; i < local.Size(); ++i )
						{
							int const u = unknown[ static_cast<std::size_t>( i ) ];
							if ( u >= 0 )
								a[ static_cast<std::size_t>( m*n + u ) ] += basis*shape( i );
							else if ( local[ i ] < vertexDofCount )
								knownRows[ static_cast<std::size_t>( m ) ][ static_cast<std::size_t>( kc++ ) ] -= basis*shape( i );
						}
					}
				}
				int const columns = static_cast<int>( block.in.size() + block.known.size() );
				std::vector<double> rhs( static_cast<std::size_t>( n*columns ), 0.0 );
				for ( int m = 0; m < n; ++m )
				{
					auto &row = inRows[ static_cast<std::size_t>( m ) ];
					row.resize( block.in.size(), 0.0 );
					for ( std::size_t c = 0; c < block.in.size(); ++c )
						rhs[ static_cast<std::size_t>( m*columns ) + c ] = row[ c ];
					for ( std::size_t c = 0; c < block.known.size(); ++c )
						rhs[ static_cast<std::size_t>( m*columns ) + block.in.size() + c ] = knownRows[ static_cast<std::size_t>( m ) ][ c ];
				}
				solveInPlace( a, rhs, n, columns );
				block.matrix = std::move( rhs );
				for ( int j = 0; j < n; ++j )
					block.out.push_back( interior[ j ] );
			}
			{
				std::vector<Block> kept;
				for ( Block &b : edgeBlocks )
					if ( !b.out.empty() )
						kept.push_back( std::move( b ) );
				stages[ 1 ] = pack( kept );
			}

			// CELLS: the interior dofs from the field's moments against
			// cellBasis(), every boundary dof known.
			std::vector<Block> cellBlocks( static_cast<std::size_t>( elements ) );
			MEQ_OMP_ELEMENT_LOOP
			for ( int e = 0; e < elements; ++e )
			{
				thread_local mfem::IsoparametricTransformation transform;
				mfem::Array<int> interior;
				out.GetElementInteriorDofs( e, interior );
				int const n = interior.Size();
				if ( n == 0 )
					continue;
				mesh.GetElementTransformation( e, &transform );
				mfem::FiniteElement const *fe = out.GetFE( e );
				mfem::Array<int> local;
				out.GetElementDofs( e, local );
				std::vector<int> unknown( static_cast<std::size_t>( local.Size() ) );
				for ( int i = 0; i < local.Size(); ++i )
					unknown[ static_cast<std::size_t>( i ) ] = interior.Find( local[ i ] );
				Block &block = cellBlocks[ static_cast<std::size_t>( e ) ];
				for ( int i = 0; i < local.Size(); ++i )
					if ( unknown[ static_cast<std::size_t>( i ) ] < 0 )
						block.known.push_back( local[ i ] );
				std::vector<double> a( static_cast<std::size_t>( n*n ), 0.0 );
				std::vector<std::vector<double>> inRows( static_cast<std::size_t>( n ) );
				std::vector<std::vector<double>> knownRows( static_cast<std::size_t>( n ),
					std::vector<double>( block.known.size(), 0.0 ) );
				mfem::Vector shape( fe->GetDof() );
				std::vector<double> basis;
				std::vector<std::pair<int, double>> terms;
				mfem::IntegrationRule const &rule = mfem::IntRules.Get( fe->GetGeomType(), 2*p + 2 );
				for ( int q = 0; q < rule.GetNPoints(); ++q )
				{
					mfem::IntegrationPoint const &ip = rule.IntPoint( q );
					transform.SetIntPoint( &ip );
					fe->CalcShape( ip, shape );
					double const weight = ip.weight*transform.Weight();
					cellBasis( fe->GetGeomType(), p, ip, basis );
					terms.clear();
					cellInput( e, ip, terms );
					for ( int m = 0; m < n; ++m )
					{
						double const test = basis[ static_cast<std::size_t>( m ) ]*weight;
						for ( auto const &[ index, value ] : terms )
						{
							int const c = columnOf( block.in, index );
							auto &row = inRows[ static_cast<std::size_t>( m ) ];
							if ( row.size() < block.in.size() )
								row.resize( block.in.size(), 0.0 );
							row[ static_cast<std::size_t>( c ) ] += test*value;
						}
						int kc = 0;
						for ( int i = 0; i < local.Size(); ++i )
						{
							int const u = unknown[ static_cast<std::size_t>( i ) ];
							if ( u >= 0 )
								a[ static_cast<std::size_t>( m*n + u ) ] += test*shape( i );
							else
								knownRows[ static_cast<std::size_t>( m ) ][ static_cast<std::size_t>( kc++ ) ] -= test*shape( i );
						}
					}
				}
				int const columns = static_cast<int>( block.in.size() + block.known.size() );
				std::vector<double> rhs( static_cast<std::size_t>( n*columns ), 0.0 );
				for ( int m = 0; m < n; ++m )
				{
					auto &row = inRows[ static_cast<std::size_t>( m ) ];
					row.resize( block.in.size(), 0.0 );
					for ( std::size_t c = 0; c < block.in.size(); ++c )
						rhs[ static_cast<std::size_t>( m*columns ) + c ] = row[ c ];
					for ( std::size_t c = 0; c < block.known.size(); ++c )
						rhs[ static_cast<std::size_t>( m*columns ) + block.in.size() + c ] = knownRows[ static_cast<std::size_t>( m ) ][ c ];
				}
				solveInPlace( a, rhs, n, columns );
				block.matrix = std::move( rhs );
				for ( int j = 0; j < n; ++j )
					block.out.push_back( interior[ j ] );
			}
			{
				std::vector<Block> kept;
				for ( Block &b : cellBlocks )
					if ( !b.out.empty() )
						kept.push_back( std::move( b ) );
				stages[ 2 ] = pack( kept );
			}
		}

		void checkSpaces( mfem::FiniteElementSpace const &potentialSpace,
		                  mfem::FiniteElementSpace const &fluxSpace )
		{
			mfem::Mesh const *mesh = potentialSpace.GetMesh();
			if ( mesh->Dimension() != 2 || fluxSpace.GetVDim() != 2 || fluxSpace.GetMesh() != mesh )
				throw std::invalid_argument( "meq::buildContinuousFields: a two-dimensional psi* and a vdim-2 "
				                             "q* on one mesh are needed" );
			if ( potentialSpace.GetFE( 0 )->GetMapType() != mfem::FiniteElement::MapType::VALUE
			     || fluxSpace.GetFE( 0 )->GetMapType() != mfem::FiniteElement::MapType::VALUE )
				throw std::invalid_argument( "meq::buildContinuousFields: psi* and q* must map by VALUE" );
			for ( int e = 0; e < mesh->GetNE(); ++e )
			{
				mfem::Geometry::Type const geometry = mesh->GetElementBaseGeometry( e );
				if ( geometry != mfem::Geometry::TRIANGLE && geometry != mfem::Geometry::SQUARE )
					throw std::invalid_argument( "meq::momentInterpolant: element " + std::to_string( e )
					                             + " is neither a triangle nor a quadrilateral" );
			}
		}
	}

	void applyLinearBlocks( LinearBlocks const &blocks, mfem::Vector const &first,
	                        mfem::Vector const *second, int split, mfem::Vector &out, int outShift )
	{
		int const *outOffset = blocks.outOffset.HostRead();
		int const *inOffset = blocks.inOffset.HostRead();
		int const *knownOffset = blocks.knownOffset.HostRead();
		int const *matrixOffset = blocks.matrixOffset.HostRead();
		int const *outIndex = blocks.outIndex.HostRead();
		int const *inIndex = blocks.inIndex.HostRead();
		int const *knownIndex = blocks.knownIndex.HostRead();
		mfem::real_t const *matrix = blocks.matrix.HostRead();
		mfem::real_t const *firstData = first.HostRead();
		mfem::real_t const *secondData = second != nullptr ? second->HostRead() : firstData;
		mfem::real_t *data = out.HostReadWrite();
		int const count = blocks.count;
		// THE KERNEL'S LOOP, OPERATION FOR OPERATION (ContinuousFieldDevice.cpp).
		MEQ_OMP_ELEMENT_LOOP
		for ( int b = 0; b < count; ++b )
		{
			int const outBegin = outOffset[ b ];
			int const rows = outOffset[ b + 1 ] - outBegin;
			int const inBegin = inOffset[ b ];
			int const inputs = inOffset[ b + 1 ] - inBegin;
			int const knownBegin = knownOffset[ b ];
			int const knowns = knownOffset[ b + 1 ] - knownBegin;
			int const columns = inputs + knowns;
			mfem::real_t const *m = matrix + matrixOffset[ b ];
			for ( int r = 0; r < rows; ++r )
			{
				mfem::real_t sum = mfem::real_t( 0 );
				for ( int c = 0; c < inputs; ++c )
				{
					int const at = inIndex[ inBegin + c ];
					mfem::real_t const value = at < split ? firstData[ at ] : secondData[ at - split ];
					sum += m[ r*columns + c ]*value;
				}
				for ( int c = 0; c < knowns; ++c )
					sum += m[ r*columns + inputs + c ]*data[ knownIndex[ knownBegin + c ] + outShift ];
				data[ outIndex[ outBegin + r ] + outShift ] = sum;
			}
		}
	}

	struct ContinuousInterpolant::OutputSpaces
	{
		mfem::H1_FECollection potentialCollection, fluxCollection;
		mfem::FiniteElementSpace potential, flux;
		OutputSpaces( mfem::Mesh *mesh, int potentialOrder, int fluxOrder )
			: potentialCollection( potentialOrder, 2 ), fluxCollection( fluxOrder, 2 ),
			  potential( mesh, &potentialCollection ),
			  flux( mesh, &fluxCollection, 2, mfem::Ordering::byNODES ) {}
	};

	ContinuousInterpolant::ContinuousInterpolant( mfem::FiniteElementSpace const &potential,
	                                              mfem::FiniteElementSpace const &flux,
	                                              mfem::FiniteElementSpace const &trace )
		: potentialSpace( &potential ), fluxSpace( &flux ), traceSpace( &trace ),
		  potentialSequence( potential.GetSequence() ), fluxSequence( flux.GetSequence() ),
		  traceSequence( trace.GetSequence() ), meshSequence( potential.GetMesh()->GetSequence() ),
		  potentialOrder( potential.GetElementOrder( 0 ) ), fluxOrder( flux.GetElementOrder( 0 ) ),
		  potentialSize( potential.GetVSize() ), fluxSize( flux.GetVSize() ), traceSize( trace.GetVSize() ),
		  elementCount( potential.GetMesh()->GetNE() )
	{
		checkSpaces( potential, flux );
		mfem::Mesh *mesh = potential.GetMesh();
		int const split = potential.GetVSize();
		output = std::make_shared<OutputSpaces>( mesh, potentialOrder + 1, fluxOrder );

		// psi~: psi* in the cells and on the boundary edges, psihat's moments on
		// the interior ones. psihat's entries sit after psi*'s.
		{
			buildStages( output->potential,
				[ & ]( int e, mfem::IntegrationPoint const &ip, std::vector<std::pair<int, double>> &terms )
				{
					addElementInput( potential, 0, e, ip, 1.0, terms );
				},
				[ & ]( int f, mfem::FaceElementTransformations &transform, mfem::IntegrationPoint const &ip,
				       std::vector<std::pair<int, double>> &terms )
				{
					if ( transform.Elem2No < 0 )
					{
						addElementInput( potential, 0, transform.Elem1No, transform.GetElement1IntPoint(), 1.0, terms );
						return;
					}
					mfem::FiniteElement const *fe = trace.GetFaceElement( f );
					mfem::Array<int> dofs;
					trace.GetFaceVDofs( f, dofs );
					mfem::Vector shape( fe->GetDof() );
					fe->CalcShape( ip, shape );
					for ( int k = 0; k < dofs.Size(); ++k )
						terms.emplace_back( split + dofs[ k ], shape( k ) );
				},
				potentialBlocks );
		}

		// q~, a component at a time: q* in the cells and the mean of the two
		// sides' traces on the edges.
		{
			mfem::H1_FECollection collection( fluxOrder, 2 );
			mfem::FiniteElementSpace out( mesh, &collection );
			for ( int component = 0; component < 2; ++component )
				buildStages( out,
					[ & ]( int e, mfem::IntegrationPoint const &ip, std::vector<std::pair<int, double>> &terms )
					{
						addElementInput( flux, component, e, ip, 1.0, terms );
					},
					[ & ]( int, mfem::FaceElementTransformations &transform, mfem::IntegrationPoint const &,
					       std::vector<std::pair<int, double>> &terms )
					{
						if ( transform.Elem2No < 0 )
						{
							addElementInput( flux, component, transform.Elem1No, transform.GetElement1IntPoint(),
							                 1.0, terms );
							return;
						}
						addElementInput( flux, component, transform.Elem1No, transform.GetElement1IntPoint(),
						                 0.5, terms );
						addElementInput( flux, component, transform.Elem2No, transform.GetElement2IntPoint(),
						                 0.5, terms );
					},
					fluxBlocks[ component ] );
		}
	}

	bool ContinuousInterpolant::builtFor( mfem::FiniteElementSpace const &potential,
	                                      mfem::FiniteElementSpace const &flux,
	                                      mfem::FiniteElementSpace const &trace ) const
	{
		return &potential == potentialSpace && &flux == fluxSpace && &trace == traceSpace
		       && potential.GetSequence() == potentialSequence && flux.GetSequence() == fluxSequence
		       && trace.GetSequence() == traceSequence
		       && potential.GetMesh()->GetSequence() == meshSequence
		       // and the sizes, which a space rebuilt at the same address with
		       // a fresh sequence would not otherwise give away
		       && potential.GetVSize() == potentialSize && flux.GetVSize() == fluxSize
		       && trace.GetVSize() == traceSize && potential.GetMesh()->GetNE() == elementCount;
	}

	ContinuousFields ContinuousInterpolant::apply( mfem::GridFunction const &potential,
	                                               mfem::GridFunction const &flux,
	                                               mfem::GridFunction const &trace,
	                                               Route route ) const
	{
		if ( !builtFor( *potential.FESpace(), *flux.FESpace(), *trace.FESpace() ) )
			throw std::logic_error( "meq::ContinuousInterpolant::apply: the fields are not on the spaces "
			                        "this interpolant was built for, or those have changed since" );

		ContinuousFields out;
		out.spaces = output;
		out.potential = std::make_unique<mfem::GridFunction>( &output->potential );
		out.flux = std::make_unique<mfem::GridFunction>( &output->flux );
		int const scalarDofs = output->flux.GetNDofs();
		int const split = potential.Size();

		bool const wantDevice = route == Route::Device
		                        || ( route == Route::Automatic
		                             && mfem::Device::Allows( mfem::Backend::DEVICE_MASK ) );
		ranOnDevice = false;
		if ( wantDevice && deviceKernelsBuilt() )
		{
			// Every entry is written by some stage, so nothing is uploaded.
			out.potential->Write();
			out.flux->Write();
			bool ok = true;
			for ( int stage = 0; stage < 3 && ok; ++stage )
				ok = applyLinearBlocksOnDevice( potentialBlocks[ stage ], potential, &trace, split,
				                                *out.potential, 0 );
			for ( int component = 0; component < 2 && ok; ++component )
				for ( int stage = 0; stage < 3 && ok; ++stage )
					ok = applyLinearBlocksOnDevice( fluxBlocks[ component ][ stage ], flux, nullptr,
					                                flux.Size(), *out.flux, component*scalarDofs );
			ranOnDevice = ok;
			// HANDED BACK: the tracer reads both fields raw, and a device read of
			// the inputs leaves their host pages protected under `debug`.
			static_cast<void>( out.potential->HostRead() );
			static_cast<void>( out.flux->HostRead() );
			static_cast<void>( potential.HostRead() );
			static_cast<void>( flux.HostRead() );
			static_cast<void>( trace.HostRead() );
			if ( ok )
				return out;
		}

		// THE HOST ROUTE, the same loop. Every entry is written by some stage.
		out.potential->HostWrite();
		out.flux->HostWrite();
		for ( int stage = 0; stage < 3; ++stage )
			applyLinearBlocks( potentialBlocks[ stage ], potential, &trace, split, *out.potential, 0 );
		for ( int component = 0; component < 2; ++component )
			for ( int stage = 0; stage < 3; ++stage )
				applyLinearBlocks( fluxBlocks[ component ][ stage ], flux, nullptr, flux.Size(), *out.flux,
				                   component*scalarDofs );
		return out;
	}

	ContinuousFields buildContinuousFields( mfem::GridFunction const &potential,
	                                        mfem::GridFunction const &flux,
	                                        mfem::GridFunction const &trace )
	{
		// One set of quadratures, used once. A caller applying the map again
		// on the same spaces keeps a ContinuousInterpolant instead -- the
		// solver does (GradShafranovSolver::continuousInterpolant()).
		ContinuousInterpolant const interpolant( *potential.FESpace(), *flux.FESpace(), *trace.FESpace() );
		return interpolant.apply( potential, flux, trace );
	}

}
