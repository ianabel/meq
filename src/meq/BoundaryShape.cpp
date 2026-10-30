#include "BoundaryShape.hpp"

#include <algorithm>
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

		/// @a angle less a multiple of 2 pi, into [ -pi, pi ): a difference of
		/// polar angles that is continuous through the cut at 0 = 2 pi.
		double wrappedDifference( double angle )
		{
			double difference = std::fmod( angle + M_PI, twoPi );
			if ( difference < 0.0 )
				difference += twoPi;
			return difference - M_PI;
		}
	}

	BoundaryShape::BoundaryShape( double radius0In, double z0In, double minorIn,
	                              double elongationIn,
	                              std::vector<double> cosIn,
	                              std::vector<double> sinIn )
		: radius0( radius0In ), z0( z0In ), minor( minorIn ), elongationValue( elongationIn ),
		  cosCoefficients( std::move( cosIn ) ),
		  sinCoefficients( std::move( sinIn ) )
	{
		if ( !( minor > 0.0 ) )
			throw ShapeError( "meq::BoundaryShape: the minor radius must be positive" );
		if ( !( elongationValue > 0.0 ) )
			throw ShapeError( "meq::BoundaryShape: the elongation must be positive" );
		if ( !( radius0 > 0.0 ) )
			throw ShapeError( "meq::BoundaryShape: the major radius must be positive" );

		// The operator carries a 1/R, which is not integrable through the axis, so
		// a surface reaching R <= 0 is not merely unusual -- it is unsolvable. The
		// bound is on the bounding box rather than on R_0 - minor, because the
		// harmonics move the innermost point.
		double minRadius = 0.0, maxRadius = 0.0, zMin = 0.0, zMax = 0.0;
		boundingBox( minRadius, maxRadius, zMin, zMax );
		if ( !( minRadius > 0.0 ) )
		{
			std::ostringstream message;
			message << "meq::BoundaryShape: the surface reaches R = " << minRadius
			        << ", which is on or beyond the axis; the Grad-Shafranov "
			           "operator's 1/R is not integrable there";
			throw ShapeError( message.str() );
		}

		requireStarShaped();
		buildPolarTable();
	}

	BoundaryShape BoundaryShape::miller( double radius0In, double z0In, double minorIn,
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

		return BoundaryShape( radius0In, z0In, minorIn, elongationIn, {}, sines );
	}

	double BoundaryShape::shiftedAngle( double theta, double &slope ) const
	{
		double shifted = theta;
		slope = 1.0;
		if ( !cosCoefficients.empty() )
			shifted += cosCoefficients[ 0 ];                 // c_0, the tilt

		std::size_t const harmonics = std::max(
			cosCoefficients.empty() ? std::size_t( 0 ) : cosCoefficients.size() - 1,
			sinCoefficients.size() );
		double const cosine = std::cos( theta ), sine = std::sin( theta );
		double cosN = 1.0, sinN = 0.0;
		for ( std::size_t n = 1; n <= harmonics; ++n )
		{
			double const cosNext = cosN*cosine - sinN*sine;
			sinN = sinN*cosine + cosN*sine;
			cosN = cosNext;
			double const order = static_cast<double>( n );
			if ( n < cosCoefficients.size() )
			{
				shifted += cosCoefficients[ n ]*cosN;
				slope -= order*cosCoefficients[ n ]*sinN;
			}
			if ( n <= sinCoefficients.size() )             // s_1 is index 0
			{
				shifted += sinCoefficients[ n - 1 ]*sinN;
				slope += order*sinCoefficients[ n - 1 ]*cosN;
			}
		}
		return shifted;
	}

	void BoundaryShape::point( double theta, double &radius, double &z ) const
	{
		double slope = 0.0;
		radius = radius0 + minor*std::cos( shiftedAngle( theta, slope ) );
		z = z0 + elongationValue*minor*std::sin( theta );
	}

	double BoundaryShape::polarAngle( double theta ) const
	{
		double radius = 0.0, z = 0.0;
		point( theta, radius, z );
		return wrappedAngle( radius - radius0, z - z0 );
	}

	void BoundaryShape::boundingBox( double &minRadius, double &maxRadius,
	                                 double &zMin, double &zMax ) const
	{
		// Sampled rather than solved. The extrema of R( theta ) satisfy
		// sin( theta_R ) theta_R' = 0, which for a general harmonic sum has no
		// closed form; 2000 points put the box within O( 1e-6 ) of the true one,
		// which is far inside the margin any caller checks it against.
		int const samples = 2000;
		double slope = 0.0;
		minRadius = maxRadius = radius0 + minor*std::cos( shiftedAngle( 0.0, slope ) );
		zMin = zMax = z0;

		for ( int i = 0; i <= samples; ++i )
		{
			double const theta = twoPi*static_cast<double>( i )
			                     /static_cast<double>( samples );
			double radius = 0.0, z = 0.0;
			point( theta, radius, z );
			minRadius = std::min( minRadius, radius );
			maxRadius = std::max( maxRadius, radius );
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
		// This is not a technicality. levelSet() inverts the polar angle by a
		// bracketed root find, and on a curve that folds back there are several
		// roots; it would return one of them and the domain would be quietly
		// wrong.
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
				           "( " << radius0 << ", " << z0 << " ): the polar angle stops "
				           "increasing at theta = " << theta << ", where it is "
				        << angle << " against " << previous << " just before. "
				           "A level set cannot be built by radial bisection on such "
				           "a curve. Reduce the harmonic coefficients, or the tilt.";
				throw ShapeError( message.str() );
			}
			previous = angle;
		}
	}

	void BoundaryShape::buildPolarTable()
	{
		polarTable.resize( tableSize + 1 );
		polarTable[ 0 ] = 0.0;
		polarTable[ tableSize ] = twoPi;
		for ( int i = 1; i < tableSize; ++i )
			polarTable[ i ] = polarAngle( twoPi*static_cast<double>( i )
			                              /static_cast<double>( tableSize ) );
	}

	double BoundaryShape::parameterAtPolarAngle( double target ) const
	{
		// The bracket is the table interval holding the target, and the start
		// is the linear interpolant across it. The polar angle runs 0 to 2 pi
		// monotonically in theta -- requireStarShaped() at construction -- so
		// the root is unique and stays bracketed.
		int low = 0, high = tableSize;
		while ( high - low > 1 )
		{
			int const middle = ( low + high )/2;
			if ( polarTable[ middle ] < target )
				low = middle;
			else
				high = middle;
		}
		double const spacing = twoPi/static_cast<double>( tableSize );
		double lower = spacing*static_cast<double>( low );
		double upper = spacing*static_cast<double>( high );
		double theta = lower + ( upper - lower )*( target - polarTable[ low ] )
		                       /( polarTable[ high ] - polarTable[ low ] );

		// Newton on the polar angle, falling back to bisection whenever a step
		// would leave the bracket. Newton is quadratic from the interpolant and
		// stops in three or four evaluations; the fallback alone would take the
		// bracket below the tolerance in under 45, so 64 is never reached by a
		// root that exists. 1e-15 in theta is below anything the caller resolves
		// and at the round-off of the curve evaluation itself.
		double const tolerance = 1.0e-15;
		for ( int iteration = 0; iteration < 64; ++iteration )
		{
			double slope = 0.0;
			double const shifted = shiftedAngle( theta, slope );
			double const dr = minor*std::cos( shifted );
			double const dz = elongationValue*minor*std::sin( theta );
			double const drdTheta = -minor*std::sin( shifted )*slope;
			double const dzdTheta = elongationValue*minor*std::cos( theta );

			double const mismatch = wrappedDifference( std::atan2( dz, dr ) - target );
			if ( mismatch == 0.0 )
				return theta;
			if ( mismatch < 0.0 )
				lower = theta;
			else
				upper = theta;

			// d/dtheta of atan2( dz, dr ), positive on a star shaped curve.
			double const derivative = ( dr*dzdTheta - dz*drdTheta )/( dr*dr + dz*dz );
			double const step = mismatch/derivative;
			if ( std::abs( step ) <= tolerance )
				return theta - step;

			double next = theta - step;
			if ( !( next > lower && next < upper ) )
				next = 0.5*( lower + upper );
			if ( upper - lower <= tolerance )
				return next;
			theta = next;
		}
		return theta;
	}

	double BoundaryShape::levelSet( double radius, double z ) const
	{
		double const dr = radius - radius0;
		double const dz = z - z0;
		double const centreDistance = std::hypot( dr, dz );

		// The centre itself. Any angle serves; the gap is the whole distance to the
		// curve, and it is negative because the centre is inside.
		if ( centreDistance < 1.0e-14 )
		{
			double br = 0.0, bz = 0.0;
			point( 0.0, br, bz );
			return -std::hypot( br - radius0, bz - z0 );
		}

		double const theta = parameterAtPolarAngle( wrappedAngle( dr, dz ) );
		double br = 0.0, bz = 0.0;
		point( theta, br, bz );
		return centreDistance - std::hypot( br - radius0, bz - z0 );
	}

	void BoundaryShape::levelSetGradient( double radius, double z,
	                                      double &dRadius, double &dZ ) const
	{
		double const dr = radius - radius0;
		double const dz = z - z0;
		double const centreDistance = std::hypot( dr, dz );
		if ( centreDistance < 1.0e-14 )
		{
			dRadius = 0.0;
			dZ = 0.0;
			return;
		}

		// The curve point q at the parameter whose polar angle is x's, and its
		// tangent dq/dtheta, both from the centre.
		double const theta = parameterAtPolarAngle( wrappedAngle( dr, dz ) );
		double slope = 0.0;
		double const shifted = shiftedAngle( theta, slope );
		double const qR = minor*std::cos( shifted );
		double const qZ = elongationValue*minor*std::sin( theta );
		double const tangentR = -minor*std::sin( shifted )*slope;
		double const tangentZ = elongationValue*minor*std::cos( theta );

		// d rho/d alpha = ( d rho/d theta )/( d alpha/d theta )
		//               = rho ( q . q' )/( q x q' ),
		// where q x q' > 0 is the star-shaped condition the constructor checks.
		double const rho = std::hypot( qR, qZ );
		double const cross = qR*tangentZ - qZ*tangentR;
		double const dRhoDAlpha = rho*( qR*tangentR + qZ*tangentZ )/cross;

		double const eR = dr/centreDistance;
		double const eZ = dz/centreDistance;
		double const angular = dRhoDAlpha/centreDistance;
		dRadius = eR + angular*eZ;    // e_alpha = ( -eZ, eR )
		dZ = eZ - angular*eR;
	}
}
