#include "BoundaryShape.hpp"

#include <cmath>
#include <sstream>

namespace meq
{
	namespace
	{
		double const twoPi = 2.0*M_PI;

		/// The polar angle of ( dr, dz ) about the origin, in [ 0, 2 pi ).
		double wrappedAngle( double dr, double dz )
		{
			double angle = std::atan2( dz, dr );
			if ( angle < 0.0 )
				angle += twoPi;
			return angle;
		}
	}

	BoundaryShape::BoundaryShape( double r0In, double z0In, double minorIn,
	                              double elongationIn,
	                              std::vector<double> cosIn,
	                              std::vector<double> sinIn )
		: r0( r0In ), z0( z0In ), minor( minorIn ), elongationValue( elongationIn ),
		  cosCoefficients( std::move( cosIn ) ),
		  sinCoefficients( std::move( sinIn ) )
	{
		if ( !( minor > 0.0 ) )
			throw ShapeError( "meq::BoundaryShape: the minor radius must be positive" );
		if ( !( elongationValue > 0.0 ) )
			throw ShapeError( "meq::BoundaryShape: the elongation must be positive" );
		if ( !( r0 > 0.0 ) )
			throw ShapeError( "meq::BoundaryShape: the major radius must be positive" );

		// The operator carries a 1/r, which is not integrable through the axis, so
		// a surface reaching r <= 0 is not merely unusual -- it is unsolvable. The
		// bound is on the bounding box rather than on r0 - minor, because the
		// harmonics move the innermost point.
		double rMin = 0.0, rMax = 0.0, zMin = 0.0, zMax = 0.0;
		boundingBox( rMin, rMax, zMin, zMax );
		if ( !( rMin > 0.0 ) )
		{
			std::ostringstream message;
			message << "meq::BoundaryShape: the surface reaches r = " << rMin
			        << ", which is on or beyond the axis; the Grad-Shafranov "
			           "operator's 1/r is not integrable there";
			throw ShapeError( message.str() );
		}

		requireStarShaped();
	}

	BoundaryShape BoundaryShape::miller( double r0In, double z0In, double minorIn,
	                                     double elongationIn, double deltaIn,
	                                     double squarenessIn )
	{
		if ( !( std::abs( deltaIn ) < 1.0 ) )
			throw ShapeError( "meq::BoundaryShape::miller: the triangularity must satisfy |delta| < 1, since it enters as arcsin( delta )" );

		// s_1 = arcsin( delta ), s_2 = -zeta: MXH eq (4). arcsin, NOT delta --
		// the two differ by 1.1 per cent at delta = 0.35 and the mistake is
		// invisible in every convergence rate.
		std::vector<double> sines{ std::asin( deltaIn ) };
		if ( squarenessIn != 0.0 )
			sines.push_back( -squarenessIn );

		return BoundaryShape( r0In, z0In, minorIn, elongationIn, {}, sines );
	}

	double BoundaryShape::shiftedAngle( double theta ) const
	{
		double shifted = theta;
		if ( !cosCoefficients.empty() )
			shifted += cosCoefficients[ 0 ];                 // c_0, the tilt

		for ( std::size_t n = 1; n < cosCoefficients.size(); ++n )
			shifted += cosCoefficients[ n ]*std::cos( static_cast<double>( n )*theta );

		for ( std::size_t i = 0; i < sinCoefficients.size(); ++i )
		{
			double const n = static_cast<double>( i + 1 );   // s_1 is index 0
			shifted += sinCoefficients[ i ]*std::sin( n*theta );
		}
		return shifted;
	}

	void BoundaryShape::point( double theta, double &r, double &z ) const
	{
		r = r0 + minor*std::cos( shiftedAngle( theta ) );
		z = z0 + elongationValue*minor*std::sin( theta );
	}

	double BoundaryShape::polarAngle( double theta ) const
	{
		double r = 0.0, z = 0.0;
		point( theta, r, z );
		return wrappedAngle( r - r0, z - z0 );
	}

	void BoundaryShape::boundingBox( double &rMin, double &rMax,
	                                 double &zMin, double &zMax ) const
	{
		// Sampled rather than solved. The extrema of R( theta ) satisfy
		// sin( theta_R ) theta_R' = 0, which for a general harmonic sum has no
		// closed form; 2000 points put the box within O( 1e-6 ) of the true one,
		// which is far inside the margin any caller checks it against.
		int const samples = 2000;
		rMin = rMax = r0 + minor*std::cos( shiftedAngle( 0.0 ) );
		zMin = zMax = z0;

		for ( int i = 0; i <= samples; ++i )
		{
			double const theta = twoPi*static_cast<double>( i )
			                     /static_cast<double>( samples );
			double r = 0.0, z = 0.0;
			point( theta, r, z );
			rMin = std::min( rMin, r );
			rMax = std::max( rMax, r );
			zMin = std::min( zMin, z );
			zMax = std::max( zMax, z );
		}
	}

	void BoundaryShape::requireStarShaped() const
	{
		// The polar angle must increase strictly with theta. It is checked on a
		// fine sample rather than proved, which is the same standard
		// tests/analytic/MillerDShape.hpp holds itself to -- but here it runs at
		// construction, because these parameters come from a user rather than
		// from a paper.
		//
		// This is not a technicality. levelSet() inverts the polar angle by
		// bisection, and on a curve that folds back there are several roots; it
		// would return one of them and the domain would be quietly wrong.
		int const samples = 4000;
		double previous = 0.0;

		for ( int i = 1; i <= samples; ++i )
		{
			double const theta = twoPi*static_cast<double>( i )
			                     /static_cast<double>( samples );
			double const angle = ( i == samples ) ? twoPi : polarAngle( theta );

			if ( !( angle > previous ) )
			{
				std::ostringstream message;
				message << "meq::BoundaryShape: the surface is not star shaped about "
				           "( " << r0 << ", " << z0 << " ): the polar angle stops "
				           "increasing at theta = " << theta << ", where it is "
				        << angle << " against " << previous << " just before. "
				           "A level set cannot be built by radial bisection on such "
				           "a curve. Reduce the harmonic coefficients, or the tilt.";
				throw ShapeError( message.str() );
			}
			previous = angle;
		}
	}

	double BoundaryShape::parameterAtPolarAngle( double target ) const
	{
		// Bisection on [ 0, 2 pi ], where the polar angle runs 0 to 2 pi
		// monotonically -- guaranteed by requireStarShaped() at construction.
		// Fifty halvings take the bracket to 2 pi * 2^-50, far below anything the
		// caller can resolve.
		double low = 0.0;
		double high = twoPi;

		for ( int i = 0; i < 50; ++i )
		{
			double const middle = 0.5*( low + high );
			double const angle = polarAngle( middle );
			if ( angle < target )
				low = middle;
			else
				high = middle;
		}
		return 0.5*( low + high );
	}

	double BoundaryShape::levelSet( double r, double z ) const
	{
		double const dr = r - r0;
		double const dz = z - z0;
		double const radius = std::hypot( dr, dz );

		// The centre itself. Any angle serves; the gap is the whole radius to the
		// curve, and it is negative because the centre is inside.
		if ( radius < 1.0e-14 )
		{
			double br = 0.0, bz = 0.0;
			point( 0.0, br, bz );
			return -std::hypot( br - r0, bz - z0 );
		}

		double const theta = parameterAtPolarAngle( wrappedAngle( dr, dz ) );
		double br = 0.0, bz = 0.0;
		point( theta, br, bz );
		return radius - std::hypot( br - r0, bz - z0 );
	}
}
