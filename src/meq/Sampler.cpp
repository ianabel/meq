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
