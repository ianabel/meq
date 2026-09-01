#include "Sampler.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace meq
{
	GridSampler::GridSampler( mfem::Mesh &meshIn,
	                          double rMinIn, double rMaxIn, int nRIn,
	                          double zMinIn, double zMaxIn, int nZIn )
		: mesh( meshIn ), rMin( rMinIn ), rMax( rMaxIn ),
		  zMin( zMinIn ), zMax( zMaxIn ), nR( nRIn ), nZ( nZIn ), found( 0 )
	{
		if ( nR < 2 || nZ < 2 )
			throw std::logic_error( "meq::GridSampler: a grid needs at least two nodes in each direction" );
		if ( !( rMax > rMin ) || !( zMax > zMin ) )
			throw std::logic_error( "meq::GridSampler: the grid extent must be positive in both directions" );

		element.assign( static_cast<std::size_t>( nR )*nZ, -1 );
		point.resize( static_cast<std::size_t>( nR )*nZ );

		double const dR = ( rMax - rMin )/( nR - 1 );
		double const dZ = ( zMax - zMin )/( nZ - 1 );

		mfem::Vector physical( 2 );
		mfem::Vector lower( 2 ), upper( 2 );

		for ( int e = 0; e < mesh.GetNE(); ++e )
		{
			// The element's bounding box. GetBoundingBox() is a mesh-wide call,
			// so the box is taken from the element's own vertices -- exact for
			// the straight-sided meshes meq builds, and an inner bound for a
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
			int const i0 = std::max( 0, static_cast<int>( std::floor( ( lower( 0 ) - rMin )/dR ) ) - 1 );
			int const i1 = std::min( nR - 1, static_cast<int>( std::ceil( ( upper( 0 ) - rMin )/dR ) ) + 1 );
			int const j0 = std::max( 0, static_cast<int>( std::floor( ( lower( 1 ) - zMin )/dZ ) ) - 1 );
			int const j1 = std::min( nZ - 1, static_cast<int>( std::ceil( ( upper( 1 ) - zMin )/dZ ) ) + 1 );
			if ( i1 < i0 || j1 < j0 )
				continue;

			mfem::ElementTransformation *transformation = mesh.GetElementTransformation( e );

			for ( int j = j0; j <= j1; ++j )
				for ( int i = i0; i <= i1; ++i )
				{
					std::size_t const at = static_cast<std::size_t>( index( i, j ) );
					if ( element[ at ] >= 0 )
						continue;      // first element to claim it keeps it

					physical( 0 ) = rAt( i );
					physical( 1 ) = zAt( j );

					mfem::IntegrationPoint reference;
					if ( transformation->TransformBack( physical, reference )
					     == mfem::InverseElementTransformation::Inside )
					{
						element[ at ] = e;
						point[ at ] = reference;
						++found;
					}
				}
		}
	}

	int GridSampler::extendOutward( double reach,
	                                std::function<bool( double, double )> const &accept )
	{
		if ( !( reach > 0.0 ) )
			return 0;

		// The boundary faces, with the element each one belongs to. These are
		// what the band is measured from: a node in the sliver is outside the
		// mesh across some face of Gamma_h, and that face's element is the one
		// whose polynomial continues into it.
		struct Face { double r0, z0, r1, z1; int element; double length; };
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

				double const r = rAt( i ), z = zAt( j );
				if ( accept && !accept( r, z ) )
					continue;

				// Nearest boundary face, by point-to-segment distance.
				int best = -1;
				double bestDistance = 0.0;
				for ( std::size_t f = 0; f < faces.size(); ++f )
				{
					Face const &face = faces[ f ];
					double const dr = face.r1 - face.r0, dz = face.z1 - face.z0;
					double const square = dr*dr + dz*dz;
					double parameter = 0.0;
					if ( square > 0.0 )
						parameter = ( ( r - face.r0 )*dr + ( z - face.z0 )*dz )/square;
					parameter = std::min( 1.0, std::max( 0.0, parameter ) );
					double const distance =
						std::hypot( r - ( face.r0 + parameter*dr ),
						            z - ( face.z0 + parameter*dz ) );
					if ( best < 0 || distance < bestDistance )
					{
						best = static_cast<int>( f );
						bestDistance = distance;
					}
				}

				// The limit is per face rather than mesh-wide, so a graded mesh
				// extrapolates further where its elements are larger, which is
				// where the discretisation error is larger anyway.
				if ( best < 0 || bestDistance > reach*faces[ best ].length )
					continue;

				mfem::Vector physical( 2 );
				physical( 0 ) = r;
				physical( 1 ) = z;
				mfem::IntegrationPoint reference;
				mfem::ElementTransformation *transformation =
					mesh.GetElementTransformation( faces[ best ].element );

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
				transformation->TransformBack( physical, reference );
				double const slack = 0.5;
				double const x = reference.x, y = reference.y;
				if ( x < -slack || y < -slack || x + y > 1.0 + slack )
					continue;

				element[ at ] = faces[ best ].element;
				point[ at ] = reference;
				++found;
				++extended;
				++filled;
			}
		return filled;
	}

	double GridSampler::rAt( int i ) const
	{
		return rMin + ( rMax - rMin )*i/( nR - 1 );
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
		for ( std::size_t at = 0; at < element.size(); ++at )
			if ( element[ at ] >= 0 )
				values[ at ] = field.GetValue( element[ at ], point[ at ] );
	}

	void GridSampler::sampleComponent( mfem::GridFunction const &field, int component,
	                                   std::vector<double> &values, double fill ) const
	{
		values.assign( static_cast<std::size_t>( nR )*nZ, fill );
		mfem::Vector vector;
		for ( std::size_t at = 0; at < element.size(); ++at )
			if ( element[ at ] >= 0 )
			{
				field.GetVectorValue( element[ at ], point[ at ], vector );
				values[ at ] = vector( component );
			}
	}

	void GridSampler::sampleCoefficient( mfem::Coefficient &coefficient,
	                                     std::vector<double> &values,
	                                     double fill ) const
	{
		values.assign( static_cast<std::size_t>( nR )*nZ, fill );
		for ( std::size_t at = 0; at < element.size(); ++at )
			if ( element[ at ] >= 0 )
			{
				mfem::ElementTransformation *transformation =
					mesh.GetElementTransformation( element[ at ] );
				transformation->SetIntPoint( &point[ at ] );
				values[ at ] = coefficient.Eval( *transformation, point[ at ] );
			}
	}
}
