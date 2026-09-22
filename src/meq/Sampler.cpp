#include "Sampler.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace meq
{
	namespace
	{
		/// Whether a field on this element can be evaluated as shape . dofs.
		///
		/// THE GROUPED PASSES BELOW DO WHAT GridFunction::GetValue() DOES, with
		/// the per-element half hoisted out of the node loop, and they are only
		/// entitled to do that for the case it handles that way: a scalar-valued
		/// element whose basis needs no transformation to evaluate. MEQ's own
		/// spaces are all L2, so this is always true here -- but sample() is a
		/// public entry point and an H(div) field handed to it must still get
		/// the right answer rather than a fast wrong one.
		bool plainNodal( mfem::FiniteElement const &fe )
		{
			return fe.GetMapType() == mfem::FiniteElement::VALUE
			       && fe.GetRangeType() == mfem::FiniteElement::SCALAR;
		}

		/// Whether the mesh's geometry map is affine on simplices -- which is
		/// what the fast inverse below needs, and is NOT the same question as
		/// whether a nodal field exists.
		///
		/// A MESH ACQUIRES NODES WITHOUT BECOMING CURVED, AND THE DIFFERENCE IS
		/// WORTH TWO ORDERS OF MAGNITUDE. Mesh::GetNodes() is null for every
		/// mesh MEQ builds with MakeCartesian2D and every order-1 .msh gmsh
		/// writes -- but any call to Mesh::EnsureNodes() installs one, and
		/// meq::FieldTransfer's constructor makes exactly that call, because
		/// FindPointsGSLIB refuses a mesh without nodes. Testing the POINTER
		/// therefore drops this sampler onto the Newton route, permanently and
		/// for the whole process, on any mesh a warm start has touched:
		/// measured at 159x at 129^2 and 83x at 513^2, with the answer
		/// unchanged, so nothing but a clock can see it.
		///
		/// An order-1 nodal field on a triangle describes the same affine map
		/// the vertices do, so the test is on the geometry's DEGREE. The values
		/// are compared with the vertex array rather than assumed equal to it:
		/// MFEM lets a caller deform a mesh through its nodes while leaving the
		/// vertices behind, and the loop below reads vertices. One pass over
		/// the vertices, against an element loop that is the whole cost here.
		bool affineGeometry( mfem::Mesh &mesh )
		{
			mfem::GridFunction const *nodes = mesh.GetNodes();
			if ( nodes == nullptr )
				return true;

			mfem::FiniteElementSpace const *space = nodes->FESpace();
			if ( space == nullptr || space->IsVariableOrder()
			     || space->GetMaxElementOrder() != 1 )
				return false;

			// Nodes and vertices must be the same points, EXACTLY: this is a
			// question about which array describes the geometry, not a
			// tolerance on how far two descriptions of it have drifted apart.
			//
			// The loop is over the nodal field's own component count rather
			// than the mesh dimension, the two differing for a surface mesh,
			// because what is being read is the field.
			int const components = space->GetVDim();
			for ( int v = 0; v < mesh.GetNV(); ++v )
			{
				double const *vertex = mesh.GetVertex( v );
				for ( int d = 0; d < components; ++d )
					if ( ( *nodes )( space->DofToVDof( v, d ) ) != vertex[ d ] )
						return false;
			}
			return true;
		}
	}

	GridSampler::GridSampler( mfem::Mesh &meshIn,
	                          double rMinIn, double rMaxIn, int nRIn,
	                          double zMinIn, double zMaxIn, int nZIn )
		: mesh( meshIn ), minRadius( rMinIn ), maxRadius( rMaxIn ),
		  zMin( zMinIn ), zMax( zMaxIn ), nR( nRIn ), nZ( nZIn ), found( 0 )
	{
		if ( nR < 2 || nZ < 2 )
			throw std::logic_error( "meq::GridSampler: a grid needs at least two nodes in each direction" );
		if ( !( maxRadius > minRadius ) || !( zMax > zMin ) )
			throw std::logic_error( "meq::GridSampler: the grid extent must be positive in both directions" );

		element.assign( static_cast<std::size_t>( nR )*nZ, -1 );
		point.resize( static_cast<std::size_t>( nR )*nZ );
		blend.assign( static_cast<std::size_t>( nR )*nZ, 0.0 );
		offsetR.assign( static_cast<std::size_t>( nR )*nZ, 0.0 );
		offsetZ.assign( static_cast<std::size_t>( nR )*nZ, 0.0 );

		double const dR = ( maxRadius - minRadius )/( nR - 1 );
		double const dZ = ( zMax - zMin )/( nZ - 1 );

		mfem::Vector physical( 2 );
		mfem::Vector lower( 2 ), upper( 2 );

		/*
		 * A STRAIGHT-SIDED ELEMENT'S MAP IS AFFINE, SO ITS INVERSE IS A 2x2
		 * SOLVE AND NOT A NEWTON ITERATION.
		 *
		 * ElementTransformation::TransformBack() constructs an
		 * InverseElementTransformation, searches a point set for a starting
		 * guess and iterates -- measured at 1.45 us per call on the DIII-D
		 * case, which is the whole cost of this constructor. For a triangle
		 * whose geometry is its three vertices the answer is exact in four
		 * multiplications, and MOST OF THOSE CALLS ARE REJECTIONS: the index
		 * range below is padded by an element in each direction, so at a grid
		 * spacing near the mesh's own only about one candidate in five is
		 * inside.
		 *
		 * A curved mesh keeps the Newton route, and so does anything that is
		 * not a triangle. What decides is the geometry's DEGREE and not whether
		 * a nodal field happens to exist -- see affineGeometry(), which is the
		 * difference between this path and a 159x slower one on any mesh a
		 * meq::FieldTransfer has been constructed on.
		 */
		bool const straightSided = affineGeometry( mesh );

		// InverseElementTransformation's own default, applied to the same test
		// Geometry::CheckPoint() applies -- so the affine path accepts exactly
		// the candidates TransformBack() accepts, rather than nearly them.
		double const insideTolerance = 1.0e-8;

		// The reentrant overload, into a local. Mesh::GetElementTransformation(
		// int ) hands out ONE shared member of the mesh, which CLAUDE.md
		// records as a silent wrong answer under threading; this loop is the
		// obvious thing to parallelise and should not have to be rewritten
		// first.
		mfem::IsoparametricTransformation transformation;

		for ( int e = 0; e < mesh.GetNE(); ++e )
		{
			// The element's bounding box. GetBoundingBox() is a mesh-wide call,
			// so the box is taken from the element's own vertices -- exact for
			// the straight-sided meshes MEQ builds, and an inner bound for a
			// curved one, which is why the containment test below is what
			// decides rather than the box.
			mfem::Array<int> vertices;
			mesh.GetElementVertices( e, vertices );
			lower( 0 ) = lower( 1 ) = 1.0e300;
			upper( 0 ) = upper( 1 ) = -1.0e300;
			for ( int v = 0; v < vertices.Size(); ++v )
			{
				double const *coordinates = mesh.GetVertex( vertices[ v ] );
				for ( int d = 0; d < 2; ++d )
				{
					lower( d ) = std::min( lower( d ), coordinates[ d ] );
					upper( d ) = std::max( upper( d ), coordinates[ d ] );
				}
			}

			// THE INVERSION: box to index range, by arithmetic. One element pad
			// so that a node sitting exactly on a boundary is not missed to
			// round-off.
			int const i0 = std::max( 0, static_cast<int>( std::floor( ( lower( 0 ) - minRadius )/dR ) ) - 1 );
			int const i1 = std::min( nR - 1, static_cast<int>( std::ceil( ( upper( 0 ) - minRadius )/dR ) ) + 1 );
			int const j0 = std::max( 0, static_cast<int>( std::floor( ( lower( 1 ) - zMin )/dZ ) ) - 1 );
			int const j1 = std::min( nZ - 1, static_cast<int>( std::ceil( ( upper( 1 ) - zMin )/dZ ) ) + 1 );
			if ( i1 < i0 || j1 < j0 )
				continue;

			// x = v0 + J ( xi, eta ), inverted once for the whole candidate
			// box. A degenerate element falls back rather than dividing by
			// zero -- Newton will report it Outside, which is the right answer.
			bool affine = straightSided && vertices.Size() == 3
			              && mesh.GetElementBaseGeometry( e ) == mfem::Geometry::TRIANGLE;
			double originR = 0.0, originZ = 0.0;
			double inverse[ 4 ] = { 0.0, 0.0, 0.0, 0.0 };
			if ( affine )
			{
				double const *p0 = mesh.GetVertex( vertices[ 0 ] );
				double const *p1 = mesh.GetVertex( vertices[ 1 ] );
				double const *p2 = mesh.GetVertex( vertices[ 2 ] );
				double const j00 = p1[ 0 ] - p0[ 0 ], j01 = p2[ 0 ] - p0[ 0 ];
				double const j10 = p1[ 1 ] - p0[ 1 ], j11 = p2[ 1 ] - p0[ 1 ];
				double const determinant = j00*j11 - j01*j10;
				if ( determinant == 0.0 )
					affine = false;
				else
				{
					originR = p0[ 0 ];
					originZ = p0[ 1 ];
					inverse[ 0 ] =  j11/determinant;
					inverse[ 1 ] = -j01/determinant;
					inverse[ 2 ] = -j10/determinant;
					inverse[ 3 ] =  j00/determinant;
				}
			}

			if ( !affine )
				mesh.GetElementTransformation( e, &transformation );

			for ( int j = j0; j <= j1; ++j )
				for ( int i = i0; i <= i1; ++i )
				{
					std::size_t const at = static_cast<std::size_t>( index( i, j ) );
					if ( element[ at ] >= 0 )
						continue;      // first element to claim it keeps it

					physical( 0 ) = rAt( i );
					physical( 1 ) = zAt( j );

					if ( affine )
					{
						double const dr = physical( 0 ) - originR;
						double const dz = physical( 1 ) - originZ;
						double const xi  = inverse[ 0 ]*dr + inverse[ 1 ]*dz;
						double const eta = inverse[ 2 ]*dr + inverse[ 3 ]*dz;
						if ( xi < -insideTolerance || eta < -insideTolerance
						     || xi + eta > 1.0 + insideTolerance )
							continue;
						point[ at ].Set2( xi, eta );
						element[ at ] = e;
						++found;
						continue;
					}

					mfem::IntegrationPoint reference;
					if ( transformation.TransformBack( physical, reference )
					     == mfem::InverseElementTransformation::Inside )
					{
						element[ at ] = e;
						point[ at ] = reference;
						++found;
					}
				}
		}
	}

	void GridSampler::buildGroups() const
	{
		if ( groupsValid )
			return;

		int const elements = mesh.GetNE();
		groupStart.assign( static_cast<std::size_t>( elements ) + 1, 0 );

		// A counting sort, which is what this is: the key is the element index
		// and it is already an integer in [ 0, NE ).
		for ( std::size_t at = 0; at < element.size(); ++at )
			if ( element[ at ] >= 0 )
				++groupStart[ static_cast<std::size_t>( element[ at ] ) + 1 ];
		for ( int e = 0; e < elements; ++e )
			groupStart[ static_cast<std::size_t>( e ) + 1 ] +=
				groupStart[ static_cast<std::size_t>( e ) ];

		groupNode.resize( static_cast<std::size_t>( groupStart.back() ) );
		std::vector<int> cursor( groupStart.begin(), groupStart.end() - 1 );
		for ( std::size_t at = 0; at < element.size(); ++at )
			if ( element[ at ] >= 0 )
				groupNode[ static_cast<std::size_t>(
					cursor[ static_cast<std::size_t>( element[ at ] ) ]++ ) ] =
					static_cast<int>( at );

		groupsValid = true;
	}

	void GridSampler::samplePotentialWithFlux( mfem::GridFunction const &potential,
	                                           mfem::GridFunction const &flux,
	                                           std::vector<double> &values,
	                                           double fill ) const
	{
		values.assign( static_cast<std::size_t>( nR )*nZ, fill );
		buildGroups();

		mfem::FiniteElementSpace const &potentialSpace = *potential.FESpace();
		mfem::FiniteElementSpace const &fluxSpace = *flux.FESpace();
		mfem::Array<int> dofs, vdofs;
		mfem::DofTransformation potentialTransform, fluxTransform;
		mfem::Vector potentialLocal, fluxLocal, shape, fluxShape;

		for ( int e = 0; e + 1 < static_cast<int>( groupStart.size() ); ++e )
		{
			int const from = groupStart[ static_cast<std::size_t>( e ) ];
			int const to = groupStart[ static_cast<std::size_t>( e ) + 1 ];
			if ( from == to )
				continue;

			mfem::FiniteElement const &fe = *potentialSpace.GetFE( e );
			if ( !plainNodal( fe ) )
			{
				for ( int k = from; k < to; ++k )
				{
					std::size_t const at =
						static_cast<std::size_t>( groupNode[ static_cast<std::size_t>( k ) ] );
					values[ at ] = potential.GetValue( e, point[ at ] );
				}
			}
			else
			{
				potentialSpace.GetElementDofs( e, dofs, potentialTransform );
				potential.GetSubVector( dofs, potentialLocal );
				if ( !potentialTransform.IsIdentity() )
					potentialTransform.InvTransformPrimal( potentialLocal );
				shape.SetSize( fe.GetDof() );

				for ( int k = from; k < to; ++k )
				{
					std::size_t const at =
						static_cast<std::size_t>( groupNode[ static_cast<std::size_t>( k ) ] );
					fe.CalcShape( point[ at ], shape );
					values[ at ] = shape*potentialLocal;
				}
			}

			// THE BAND, AND ONLY IF THIS ELEMENT HAS ANY. Interior nodes have a
			// zero offset, so on every element inside Gamma_h the flux is never
			// fetched at all -- which on the ordinary run is all of them.
			bool banded = false;
			for ( int k = from; k < to && !banded; ++k )
			{
				std::size_t const at =
					static_cast<std::size_t>( groupNode[ static_cast<std::size_t>( k ) ] );
				banded = ( offsetR[ at ] != 0.0 || offsetZ[ at ] != 0.0 );
			}
			if ( !banded )
				continue;

			mfem::FiniteElement const &fluxFe = *fluxSpace.GetFE( e );
			bool const fluxNodal = plainNodal( fluxFe );
			int const fluxDof = fluxFe.GetDof();
			if ( fluxNodal )
			{
				fluxSpace.GetElementVDofs( e, vdofs, fluxTransform );
				flux.GetSubVector( vdofs, fluxLocal );
				if ( !fluxTransform.IsIdentity() )
					fluxTransform.InvTransformPrimal( fluxLocal );
				fluxShape.SetSize( fluxDof );
			}

			mfem::Vector q( 2 );
			for ( int k = from; k < to; ++k )
			{
				std::size_t const at =
					static_cast<std::size_t>( groupNode[ static_cast<std::size_t>( k ) ] );
				if ( offsetR[ at ] == 0.0 && offsetZ[ at ] == 0.0 )
					continue;

				if ( fluxNodal )
				{
					fluxFe.CalcShape( point[ at ], fluxShape );
					q( 0 ) = fluxShape*( &fluxLocal[ 0 ] );
					q( 1 ) = fluxShape*( &fluxLocal[ fluxDof ] );
				}
				else
					flux.GetVectorValue( e, point[ at ], q );

				// grad psi = R q, with R taken at the foot -- the point the
				// flux was actually read at.
				int const i = static_cast<int>( at ) % nR;
				double const footR = rAt( i ) - offsetR[ at ];
				values[ at ] += footR*( q( 0 )*offsetR[ at ] + q( 1 )*offsetZ[ at ] );
			}
		}
	}

	void GridSampler::sampleComponentWithGradient( mfem::GridFunction const &field,
	                                               int component,
	                                               std::vector<double> &values,
	                                               double fill ) const
	{
		values.assign( static_cast<std::size_t>( nR )*nZ, fill );
		buildGroups();

		mfem::FiniteElementSpace const &space = *field.FESpace();
		mfem::Array<int> vdofs;
		mfem::DofTransformation doftrans;
		mfem::Vector local, shape, vector;
		mfem::DenseMatrix gradient;
		mfem::IsoparametricTransformation transformation;

		for ( int e = 0; e + 1 < static_cast<int>( groupStart.size() ); ++e )
		{
			int const from = groupStart[ static_cast<std::size_t>( e ) ];
			int const to = groupStart[ static_cast<std::size_t>( e ) + 1 ];
			if ( from == to )
				continue;

			mfem::FiniteElement const &fe = *space.GetFE( e );
			bool const nodal = plainNodal( fe );
			int const dof = fe.GetDof();
			if ( nodal )
			{
				space.GetElementVDofs( e, vdofs, doftrans );
				field.GetSubVector( vdofs, local );
				if ( !doftrans.IsIdentity() )
					doftrans.InvTransformPrimal( local );
				shape.SetSize( dof );
			}

			bool transformed = false;
			for ( int k = from; k < to; ++k )
			{
				std::size_t const at =
					static_cast<std::size_t>( groupNode[ static_cast<std::size_t>( k ) ] );

				double value;
				if ( nodal )
				{
					fe.CalcShape( point[ at ], shape );
					value = shape*( &local[ dof*component ] );
				}
				else
				{
					field.GetVectorValue( e, point[ at ], vector );
					value = vector( component );
				}

				// Interior nodes have a zero offset, so this costs them a branch
				// and nothing else -- the same shape as samplePotentialWithFlux().
				if ( offsetR[ at ] != 0.0 || offsetZ[ at ] != 0.0 )
				{
					if ( !transformed )
					{
						mesh.GetElementTransformation( e, &transformation );
						transformed = true;
					}
					transformation.SetIntPoint( &point[ at ] );
					field.GetVectorGradient( transformation, gradient );

					// grad( i, j ) = d u_i / d x_j: COMPONENT first, DIRECTION
					// second. Transposing it is wrong at every node whose gradient
					// is not symmetric and silent everywhere else, so the test field
					// in theBandVectorContinuesAtItsGradientsOrder is deliberately
					// one with an asymmetric gradient.
					value += gradient( component, 0 )*offsetR[ at ]
					         + gradient( component, 1 )*offsetZ[ at ];
				}
				values[ at ] = value;
			}
		}
	}

	bool GridSampler::wasExtended( int i, int j ) const
	{
		std::size_t const at = static_cast<std::size_t>( index( i, j ) );
		return element[ at ] >= 0
		       && ( offsetR[ at ] != 0.0 || offsetZ[ at ] != 0.0 );
	}

	double GridSampler::blendWeight( int i, int j ) const
	{
		return blend[ static_cast<std::size_t>( index( i, j ) ) ];
	}

	int GridSampler::extendOutward( double reach,
	                                std::function<bool( double, double )> const &accept,
	                                std::function<double( double, double )> const
	                                    &gapToBoundary )
	{
		if ( !( reach > 0.0 ) )
			return 0;

		// The boundary faces, with the element each one belongs to. These are
		// what the band is measured from: a node in the sliver is outside the
		// mesh across some face of Gamma_h, and that face's element is the one
		// whose polynomial continues into it.
		struct Face { double radius0, z0, radius1, z1; int element; double length; };
		std::vector<Face> faces;
		faces.reserve( mesh.GetNBE() );
		for ( int b = 0; b < mesh.GetNBE(); ++b )
		{
			mfem::Array<int> vertices;
			mesh.GetBdrElementVertices( b, vertices );
			if ( vertices.Size() != 2 )
				continue;               // not a 2D mesh; nothing to do here
			int element = -1, info = 0;
			mesh.GetBdrElementAdjacentElement( b, element, info );
			if ( element < 0 )
				continue;

			double const *a = mesh.GetVertex( vertices[ 0 ] );
			double const *c = mesh.GetVertex( vertices[ 1 ] );
			double const length = std::hypot( c[ 0 ] - a[ 0 ], c[ 1 ] - a[ 1 ] );
			faces.push_back( Face{ a[0], a[1], c[0], c[1], element, length } );
		}
		if ( faces.empty() )
			return 0;

		int filled = 0;
		for ( int j = 0; j < nZ; ++j )
			for ( int i = 0; i < nR; ++i )
			{
				std::size_t const at = static_cast<std::size_t>( index( i, j ) );
				if ( element[ at ] >= 0 )
					continue;

				double const radius = rAt( i ), z = zAt( j );
				if ( accept && !accept( radius, z ) )
					continue;

				// Nearest boundary face, and the FOOT of the node on it.
				int best = -1;
				double bestDistance = 0.0, footR = radius, footZ = z;
				for ( std::size_t f = 0; f < faces.size(); ++f )
				{
					Face const &face = faces[ f ];
					double const dr = face.radius1 - face.radius0, dz = face.z1 - face.z0;
					double const square = dr*dr + dz*dz;
					double parameter = 0.0;
					if ( square > 0.0 )
						parameter = ( ( radius - face.radius0 )*dr + ( z - face.z0 )*dz )/square;
					parameter = std::min( 1.0, std::max( 0.0, parameter ) );
					double const onFaceR = face.radius0 + parameter*dr;
					double const onFaceZ = face.z0 + parameter*dz;
					double const distance = std::hypot( radius - onFaceR, z - onFaceZ );
					if ( best < 0 || distance < bestDistance )
					{
						best = static_cast<int>( f );
						bestDistance = distance;
						footR = onFaceR;
						footZ = onFaceZ;
					}
				}

				// The limit is per face rather than mesh-wide, so a graded mesh
				// extrapolates further where its elements are larger, which is
				// where the discretisation error is larger anyway.
				if ( best < 0 || bestDistance > reach*faces[ best ].length )
					continue;

				// THE FOOT, NOT THE NODE. x0 lies on the element's own
				// boundary, so the inverse map converges and every field
				// evaluation is inside the element that owns it. Inverting at
				// the node instead is what made this an extrapolation.
				mfem::Vector physical( 2 );
				physical( 0 ) = footR;
				physical( 1 ) = footZ;
				mfem::IntegrationPoint reference;
				mfem::IsoparametricTransformation transformation;
				mesh.GetElementTransformation( faces[ best ].element,
				                               &transformation );

				// TransformBack reports Outside for exactly the nodes this
				// function is for, and still returns reference coordinates --
				// evaluating the basis there is the extrapolation. So the
				// RESULT CODE is not the test.
				//
				// THE REFERENCE POINT IS, THOUGH, and skipping this check was a
				// measured mistake. InverseElementTransformation runs a Newton
				// iteration, and for a point outside the element it can fail to
				// converge and leave `reference` anywhere at all. The basis is
				// then evaluated a long way outside its element and returns
				// whatever a degree-k polynomial does out there -- which is not
				// an approximation of anything. Measured on the Miller case
				// before this guard: psi overshot past zero to +1.4e-02 against
				// a peak of 2.5e-01, about a hundred times the O( h^(k+1) ) an
				// honest one-element extrapolation gives.
				//
				// The reference triangle is 0 <= x, y and x + y <= 1, so a point
				// within `slack` of it is at most that far outside in reference
				// units. One half is generous for a band one face deep.
				transformation.TransformBack( physical, reference );
				// The foot is ON the element, so this should always pass; it is
				// kept because a degenerate face could still defeat the inverse
				// map, and a silent wild reference point is the failure this
				// whole function had before.
				double const slack = 0.5;
				double const x = reference.x, y = reference.y;
				if ( x < -slack || y < -slack || x + y > 1.0 + slack )
					continue;

				element[ at ] = faces[ best ].element;
				point[ at ] = reference;
				offsetR[ at ] = radius - footR;
				offsetZ[ at ] = z - footZ;

				// How far through the band: 0 where it meets Gamma_h, 1 on
				// Gamma. bestDistance is the distance to Gamma_h and the caller
				// supplies the distance to Gamma, so the two bracket the node
				// and their ratio needs no geometry this class does not have.
				if ( gapToBoundary )
				{
					double const toGamma = std::max( 0.0, gapToBoundary( radius, z ) );
					double const total = bestDistance + toGamma;
					blend[ at ] = total > 0.0
						? std::min( 1.0, bestDistance/total )
						: 1.0;
				}
				++found;
				++extended;
				++filled;
			}

		// The band nodes are new members of their elements' groups, so whatever
		// a previous sampling pass built is now short of them.
		if ( filled > 0 )
			groupsValid = false;
		return filled;
	}

	double GridSampler::rAt( int i ) const
	{
		return minRadius + ( maxRadius - minRadius )*i/( nR - 1 );
	}

	double GridSampler::zAt( int j ) const
	{
		return zMin + ( zMax - zMin )*j/( nZ - 1 );
	}

	bool GridSampler::located( int i, int j ) const
	{
		return element[ static_cast<std::size_t>( index( i, j ) ) ] >= 0;
	}

	void GridSampler::sample( mfem::GridFunction const &field,
	                          std::vector<double> &values, double fill ) const
	{
		values.assign( static_cast<std::size_t>( nR )*nZ, fill );
		buildGroups();

		mfem::FiniteElementSpace const &space = *field.FESpace();
		mfem::Array<int> dofs;
		mfem::DofTransformation doftrans;
		mfem::Vector local, shape;

		for ( int e = 0; e + 1 < static_cast<int>( groupStart.size() ); ++e )
		{
			int const from = groupStart[ static_cast<std::size_t>( e ) ];
			int const to = groupStart[ static_cast<std::size_t>( e ) + 1 ];
			if ( from == to )
				continue;

			mfem::FiniteElement const &fe = *space.GetFE( e );
			if ( !plainNodal( fe ) )
			{
				for ( int k = from; k < to; ++k )
				{
					std::size_t const at =
						static_cast<std::size_t>( groupNode[ static_cast<std::size_t>( k ) ] );
					values[ at ] = field.GetValue( e, point[ at ] );
				}
				continue;
			}

			space.GetElementDofs( e, dofs, doftrans );
			field.GetSubVector( dofs, local );
			if ( !doftrans.IsIdentity() )
				doftrans.InvTransformPrimal( local );
			shape.SetSize( fe.GetDof() );

			for ( int k = from; k < to; ++k )
			{
				std::size_t const at =
					static_cast<std::size_t>( groupNode[ static_cast<std::size_t>( k ) ] );
				fe.CalcShape( point[ at ], shape );
				values[ at ] = shape*local;
			}
		}
	}

	void GridSampler::sampleComponent( mfem::GridFunction const &field, int component,
	                                   std::vector<double> &values, double fill ) const
	{
		values.assign( static_cast<std::size_t>( nR )*nZ, fill );
		buildGroups();

		mfem::FiniteElementSpace const &space = *field.FESpace();
		mfem::Array<int> vdofs;
		mfem::DofTransformation doftrans;
		mfem::Vector local, shape, vector;

		for ( int e = 0; e + 1 < static_cast<int>( groupStart.size() ); ++e )
		{
			int const from = groupStart[ static_cast<std::size_t>( e ) ];
			int const to = groupStart[ static_cast<std::size_t>( e ) + 1 ];
			if ( from == to )
				continue;

			mfem::FiniteElement const &fe = *space.GetFE( e );
			if ( !plainNodal( fe ) )
			{
				for ( int k = from; k < to; ++k )
				{
					std::size_t const at =
						static_cast<std::size_t>( groupNode[ static_cast<std::size_t>( k ) ] );
					field.GetVectorValue( e, point[ at ], vector );
					values[ at ] = vector( component );
				}
				continue;
			}

			int const dof = fe.GetDof();
			space.GetElementVDofs( e, vdofs, doftrans );
			field.GetSubVector( vdofs, local );
			if ( !doftrans.IsIdentity() )
				doftrans.InvTransformPrimal( local );
			shape.SetSize( dof );

			for ( int k = from; k < to; ++k )
			{
				std::size_t const at =
					static_cast<std::size_t>( groupNode[ static_cast<std::size_t>( k ) ] );
				fe.CalcShape( point[ at ], shape );
				values[ at ] = shape*( &local[ dof*component ] );
			}
		}
	}

	void GridSampler::sampleCoefficient( mfem::Coefficient &coefficient,
	                                     std::vector<double> &values,
	                                     double fill ) const
	{
		values.assign( static_cast<std::size_t>( nR )*nZ, fill );
		buildGroups();

		// One transformation per ELEMENT rather than per node, and a local one
		// rather than the mesh's shared member -- see the constructor.
		mfem::IsoparametricTransformation transformation;

		for ( int e = 0; e + 1 < static_cast<int>( groupStart.size() ); ++e )
		{
			int const from = groupStart[ static_cast<std::size_t>( e ) ];
			int const to = groupStart[ static_cast<std::size_t>( e ) + 1 ];
			if ( from == to )
				continue;

			mesh.GetElementTransformation( e, &transformation );
			for ( int k = from; k < to; ++k )
			{
				std::size_t const at =
					static_cast<std::size_t>( groupNode[ static_cast<std::size_t>( k ) ] );
				transformation.SetIntPoint( &point[ at ] );
				values[ at ] = coefficient.Eval( transformation, point[ at ] );
			}
		}
	}
}
