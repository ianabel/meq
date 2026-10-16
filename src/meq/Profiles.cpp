#include "Profiles.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <istream>
#include <limits>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace
{

	// The Hermite basis on the unit interval, and its derivatives with respect to
	// the local coordinate t = ( x - x_l )/delta.
	//
	//   h00 = 2t^3 - 3t^2 + 1   h10 = t^3 - 2t^2 + t
	//   h01 = -2t^3 + 3t^2      h11 = t^3 - t^2

	inline double h00( double t )
	{
		return ( 1.0 + 2.0*t )*( 1.0 - t )*( 1.0 - t );
	}

	inline double h10( double t )
	{
		return t*( 1.0 - t )*( 1.0 - t );
	}

	inline double h01( double t )
	{
		return t*t*( 3.0 - 2.0*t );
	}

	inline double h11( double t )
	{
		return t*t*( t - 1.0 );
	}

	inline double h00Prime( double t )
	{
		return 6.0*t*( t - 1.0 );
	}

	inline double h10Prime( double t )
	{
		return 3.0*t*t - 4.0*t + 1.0;
	}

	inline double h01Prime( double t )
	{
		return 6.0*t*( 1.0 - t );
	}

	inline double h11Prime( double t )
	{
		return t*( 3.0*t - 2.0 );
	}

	// ... and the second derivatives, which a source whose F is itself a
	// psi-derivative needs. Linear in t, hence the jump at a shared knot.

	inline double h00DoublePrime( double t )
	{
		return 12.0*t - 6.0;
	}

	inline double h10DoublePrime( double t )
	{
		return 6.0*t - 4.0;
	}

	inline double h01DoublePrime( double t )
	{
		return 6.0 - 12.0*t;
	}

	inline double h11DoublePrime( double t )
	{
		return 6.0*t - 2.0;
	}

	// Reject the values that would turn into a silently wrong interpolant later.
	void checkTable( std::vector<meq::Knot> const & data )
	{
		if ( data.size() < 2 )
			throw std::invalid_argument( "meq::SplineProfile: a spline needs at least two knots" );

		for ( std::size_t i = 0; i < data.size(); ++i )
		{
			if ( !std::isfinite( data[ i ].psi ) || !std::isfinite( data[ i ].value ) || !std::isfinite( data[ i ].derivative ) )
				throw std::invalid_argument( "meq::SplineProfile: knot data must be finite" );

			if ( i > 0 && !( data[ i - 1 ].psi < data[ i ].psi ) )
				throw std::invalid_argument( "meq::SplineProfile: knots must be strictly increasing in psi" );
		}
	}

}

namespace meq
{

	ConstantProfile::ConstantProfile( double value )
		: constantValue( value )
	{
	}

	double ConstantProfile::operator()( double ) const
	{
		return constantValue;
	}

	double ConstantProfile::prime( double ) const
	{
		return 0.0;
	}

	double ConstantProfile::doublePrime( double ) const
	{
		return 0.0;
	}

	double ConstantProfile::value() const
	{
		return constantValue;
	}

	ScaledProfile::ScaledProfile( std::shared_ptr<Profile const> inner, double scale )
		: innerProfile( std::move( inner ) ), scaleFactor( scale )
	{
		if ( !innerProfile )
			throw std::invalid_argument( "meq::ScaledProfile: the profile to scale must not be null" );
		if ( !std::isfinite( scaleFactor ) )
			throw std::invalid_argument( "meq::ScaledProfile: the scale must be finite" );
	}

	double ScaledProfile::operator()( double psi ) const
	{
		return scaleFactor*( *innerProfile )( psi );
	}

	double ScaledProfile::prime( double psi ) const
	{
		return scaleFactor*innerProfile->prime( psi );
	}

	double ScaledProfile::doublePrime( double psi ) const
	{
		return scaleFactor*innerProfile->doublePrime( psi );
	}

	Profile const & ScaledProfile::unscaled() const
	{
		return *innerProfile;
	}

	double ScaledProfile::scale() const
	{
		return scaleFactor;
	}

	HermiteCubicSpline::HermiteCubicSpline( double lower, double upper, double fLower, double fUpper, double fPrimeLower, double fPrimeUpper )
		: xLower( lower ), xUpper( upper ), delta( upper - lower ),
		  fLower( fLower ), fUpper( fUpper ),
		  fPrimeLower( fPrimeLower ), fPrimeUpper( fPrimeUpper )
	{
		if ( !( lower < upper ) )
			throw std::invalid_argument( "meq::HermiteCubicSpline: the interval must have positive width" );
	}

	HermiteCubicSpline::HermiteCubicSpline( Knot const & lower, Knot const & upper )
		: HermiteCubicSpline( lower.psi, upper.psi, lower.value, upper.value, lower.derivative, upper.derivative )
	{
	}

	double HermiteCubicSpline::operator()( double x ) const
	{
		if ( x <= xLower )
			return fLower;
		if ( x >= xUpper )
			return fUpper;

		double const t = ( x - xLower )/delta;

		return fLower*h00( t ) + fPrimeLower*delta*h10( t ) + fUpper*h01( t ) + fPrimeUpper*delta*h11( t );
	}

	double HermiteCubicSpline::prime( double x ) const
	{
		// Strictly outside the interval the interpolant is extended by a constant,
		// whose derivative is zero. At the endpoints themselves the tabulated
		// one-sided derivative is the right answer, so the comparisons here are
		// strict where operator()'s are not.
		if ( x < xLower || x > xUpper )
			return 0.0;

		double const t = ( x - xLower )/delta;

		return ( fLower*h00Prime( t ) + fPrimeLower*delta*h10Prime( t ) + fUpper*h01Prime( t ) + fPrimeUpper*delta*h11Prime( t ) )/delta;
	}

	double HermiteCubicSpline::doublePrime( double x ) const
	{
		// Same convention as prime(): the constant extension outside the interval
		// has a zero second derivative, and the comparisons are strict so that an
		// endpoint gets this interval's answer rather than the extension's.
		if ( x < xLower || x > xUpper )
			return 0.0;

		double const t = ( x - xLower )/delta;

		return ( fLower*h00DoublePrime( t ) + fPrimeLower*delta*h10DoublePrime( t )
		         + fUpper*h01DoublePrime( t ) + fPrimeUpper*delta*h11DoublePrime( t ) )/( delta*delta );
	}

	std::pair<double,double> HermiteCubicSpline::interval() const
	{
		return std::make_pair( xLower, xUpper );
	}

	std::pair<double,double> HermiteCubicSpline::values() const
	{
		return std::make_pair( fLower, fUpper );
	}

	std::pair<double,double> HermiteCubicSpline::derivatives() const
	{
		return std::make_pair( fPrimeLower, fPrimeUpper );
	}

	SplineProfile::SplineProfile( std::vector<Knot> data )
		: knotData( std::move( data ) )
	{
		checkTable( knotData );
	}

	SplineProfile::SplineProfile( RealFunction f, RealFunction fPrime, unsigned int intervals )
	{
		if ( intervals == 0 )
			throw std::invalid_argument( "meq::SplineProfile: need at least one interval" );
		if ( !f || !fPrime )
			throw std::invalid_argument( "meq::SplineProfile: both the function and its derivative must be supplied" );

		double const deltaPsi = 1.0/intervals;

		knotData.reserve( intervals + 1 );
		for ( unsigned int i = 0; i <= intervals; ++i )
		{
			// The last knot is set to exactly 1.0 rather than intervals*deltaPsi so
			// that domain() is exactly [ 0, 1 ].
			double const psi = ( i == intervals ) ? 1.0 : i*deltaPsi;
			knotData.push_back( Knot{ psi, f( psi ), fPrime( psi ) } );
		}

		checkTable( knotData );
	}

	SplineProfile SplineProfile::fromStream( std::istream & is )
	{
		std::vector<Knot> data;
		std::string line;

		while ( std::getline( is, line ) )
		{
			std::size_t const firstNonBlank = line.find_first_not_of( " \t\r\f\v" );

			// A blank line marks the end of the table, so that several profiles can
			// live in one stream.
			if ( firstNonBlank == std::string::npos )
				break;

			if ( line[ firstNonBlank ] == '#' )
				continue;

			std::istringstream lineStream( line );
			Knot knot{ 0.0, 0.0, 0.0 };
			if ( !( lineStream >> knot.psi >> knot.value >> knot.derivative ) )
				throw std::runtime_error( "meq::SplineProfile: cannot parse 'psi f(psi) f'(psi)' from line: " + line );

			data.push_back( knot );
		}

		if ( data.size() < 2 )
			throw std::runtime_error( "meq::SplineProfile: fewer than two knots read from stream" );

		return SplineProfile( std::move( data ) );
	}

	SplineProfile SplineProfile::fromFile( std::string const & fileName )
	{
		std::ifstream file( fileName );

		if ( !file )
			throw std::runtime_error( "meq::SplineProfile::fromFile: cannot open '" + fileName + "'" );

		try
		{
			return fromStream( file );
		}
		catch ( std::exception const & e )
		{
			throw std::runtime_error( "meq::SplineProfile::fromFile: error reading '" + fileName + "': " + e.what() );
		}
	}

	std::size_t SplineProfile::findInterval( double psi ) const
	{
		// First knot strictly beyond psi; the interval we want starts at the knot
		// before it. A psi sitting exactly on an interior knot lands on the
		// interval starting there, and since the interpolant matches both value and
		// derivative at every knot, the neighbouring interval would have answered
		// identically.
		auto const above = std::upper_bound( knotData.begin(), knotData.end(), psi,
			[]( double x, Knot const & knot ) { return x < knot.psi; } );

		if ( above == knotData.begin() )
			return 0;

		std::size_t const index = static_cast<std::size_t>( above - knotData.begin() ) - 1;

		return std::min( index, numIntervals() - 1 );
	}

	double SplineProfile::operator()( double psi ) const
	{
		if ( psi <= knotData.front().psi )
			return knotData.front().value;
		if ( psi >= knotData.back().psi )
			return knotData.back().value;

		return intervalAt( findInterval( psi ) )( psi );
	}

	double SplineProfile::prime( double psi ) const
	{
		if ( psi < knotData.front().psi || psi > knotData.back().psi )
			return 0.0;

		return intervalAt( findInterval( psi ) ).prime( psi );
	}

	double SplineProfile::doublePrime( double psi ) const
	{
		if ( psi < knotData.front().psi || psi > knotData.back().psi )
			return 0.0;

		return intervalAt( findInterval( psi ) ).doublePrime( psi );
	}

	bool SplineProfile::uniformKnots() const
	{
		if ( knotData.size() < 3 )
			return true;

		double const h = knotData[ 1 ].psi - knotData[ 0 ].psi;
		if ( !( h > 0.0 ) )
			return false;

		for ( std::size_t i = 1; i + 1 < knotData.size(); ++i )
			if ( std::abs( ( knotData[ i + 1 ].psi - knotData[ i ].psi ) - h )
			     > 1.0e-9*std::abs( h ) )
				return false;

		return true;
	}

	SplineProfile SplineProfile::resampledUniform( unsigned int intervals ) const
	{
		if ( intervals == 0 )
			throw std::invalid_argument(
				"meq::SplineProfile::resampledUniform: need at least one interval" );

		double const lower = knotData.front().psi;
		double const upper = knotData.back().psi;
		double const delta = ( upper - lower )/intervals;

		std::vector<Knot> data;
		data.reserve( intervals + 1 );
		for ( unsigned int i = 0; i <= intervals; ++i )
		{
			// The last abscissa is written as the endpoint rather than
			// accumulated, so the resampled range ENDS where this one does --
			// lower + intervals*delta is the same number only in exact
			// arithmetic, and a last knot a few ulps short moves the constant
			// extension that operator() applies beyond it.
			double const psi = ( i == intervals ) ? upper : lower + i*delta;
			data.push_back( Knot{ psi, ( *this )( psi ), prime( psi ) } );
		}

		return SplineProfile( std::move( data ) );
	}

	void SplineProfile::write( std::ostream & os ) const
	{
		// max_digits10 so that write() followed by fromStream() is the identity.
		std::streamsize const oldPrecision = os.precision( std::numeric_limits<double>::max_digits10 );

		os << "# MEQ spline profile: " << knotData.size() << " knots\n";
		os << "# psi\tf(psi)\tf'(psi)\n";

		for ( auto const & knot : knotData )
			os << knot.psi << "\t" << knot.value << "\t" << knot.derivative << "\n";

		// Blank line marks end-of-data.
		os << std::endl;

		os.precision( oldPrecision );
	}

	std::vector<Knot> const & SplineProfile::knots() const
	{
		return knotData;
	}

	std::pair<double,double> SplineProfile::domain() const
	{
		return std::make_pair( knotData.front().psi, knotData.back().psi );
	}

	std::size_t SplineProfile::numIntervals() const
	{
		return knotData.size() - 1;
	}

	HermiteCubicSpline SplineProfile::intervalAt( std::size_t i ) const
	{
		if ( i >= numIntervals() )
			throw std::out_of_range( "meq::SplineProfile::intervalAt: no such interval" );

		return HermiteCubicSpline( knotData[ i ], knotData[ i + 1 ] );
	}

	std::ostream & operator<<( std::ostream & os, SplineProfile const & profile )
	{
		profile.write( os );
		return os;
	}


	CardinalSpline::CardinalSpline( double lower, double upper, std::size_t knots )
		: xLower( lower ), xUpper( upper ), spacing( 0.0 ), count( knots )
	{
		if ( knots < 4 )
			throw std::invalid_argument( "meq::CardinalSpline: not-a-knot needs at least four knots" );
		if ( !std::isfinite( lower ) || !std::isfinite( upper ) || !( lower < upper ) )
			throw std::invalid_argument( "meq::CardinalSpline: the range must be finite and increasing" );
		spacing = ( upper - lower )/static_cast<double>( knots - 1 );

		/*
		 * THE CURVATURE MAP, BY ONE DENSE ELIMINATION ON THE SPLINE'S OWN SYSTEM.
		 *
		 * Interior rows are C^2 continuity on uniform knots,
		 *
		 *     M_{i-1} + 4 M_i + M_{i+1} = 6( f_{i-1} - 2 f_i + f_{i+1} )/h^2,
		 *
		 * and the end rows not-a-knot -- f''' continuous across knots 1 and
		 * K - 2, so M_0 - 2 M_1 + M_2 = 0 and its mirror. Solved once for the
		 * identity right hand side, which is every phi_j's curvatures at once.
		 * Partial pivoting: the end rows are not diagonally dominant.
		 */
		std::size_t const n = count;
		std::vector<double> system( n*n, 0.0 ), right( n*n, 0.0 );
		auto at = [ n ]( std::vector<double> &m, std::size_t i, std::size_t j ) -> double &
		{
			return m[ i*n + j ];
		};
		double const scale = 6.0/( spacing*spacing );
		at( system, 0, 0 ) = 1.0;
		at( system, 0, 1 ) = -2.0;
		at( system, 0, 2 ) = 1.0;
		at( system, n - 1, n - 3 ) = 1.0;
		at( system, n - 1, n - 2 ) = -2.0;
		at( system, n - 1, n - 1 ) = 1.0;
		for ( std::size_t i = 1; i + 1 < n; ++i )
		{
			at( system, i, i - 1 ) = 1.0;
			at( system, i, i ) = 4.0;
			at( system, i, i + 1 ) = 1.0;
			at( right, i, i - 1 ) = scale;
			at( right, i, i ) = -2.0*scale;
			at( right, i, i + 1 ) = scale;
		}
		for ( std::size_t c = 0; c < n; ++c )
		{
			std::size_t pivot = c;
			for ( std::size_t r = c + 1; r < n; ++r )
				if ( std::fabs( at( system, r, c ) ) > std::fabs( at( system, pivot, c ) ) )
					pivot = r;
			if ( pivot != c )
				for ( std::size_t j = 0; j < n; ++j )
				{
					std::swap( at( system, c, j ), at( system, pivot, j ) );
					std::swap( at( right, c, j ), at( right, pivot, j ) );
				}
			double const d = at( system, c, c );
			for ( std::size_t r = c + 1; r < n; ++r )
			{
				double const factor = at( system, r, c )/d;
				if ( factor == 0.0 )
					continue;
				for ( std::size_t j = c; j < n; ++j )
					at( system, r, j ) -= factor*at( system, c, j );
				for ( std::size_t j = 0; j < n; ++j )
					at( right, r, j ) -= factor*at( right, c, j );
			}
		}
		for ( std::size_t c = n; c-- > 0; )
			for ( std::size_t j = 0; j < n; ++j )
			{
				double value = at( right, c, j );
				for ( std::size_t k = c + 1; k < n; ++k )
					value -= at( system, c, k )*at( right, k, j );
				at( right, c, j ) = value/at( system, c, c );
			}
		curvatureMap = std::move( right );
	}

	std::size_t CardinalSpline::size() const
	{
		return count;
	}

	double CardinalSpline::lower() const
	{
		return xLower;
	}

	double CardinalSpline::upper() const
	{
		return xUpper;
	}

	double CardinalSpline::knot( std::size_t k ) const
	{
		return k + 1 == count ? xUpper : xLower + spacing*static_cast<double>( k );
	}

	std::size_t CardinalSpline::interval( double x ) const
	{
		double const t = ( x - xLower )/spacing;
		if ( !( t > 0.0 ) )
			return 0;
		return std::min( static_cast<std::size_t>( t ), count - 2 );
	}

	void CardinalSpline::basis( double x, double *value, double *slope, double *curvature ) const
	{
		// Clamped outside, with the derivatives of the constant extension.
		bool const outside = x < xLower || x > xUpper;
		double const xc = std::min( xUpper, std::max( xLower, x ) );
		std::size_t const i = interval( xc );
		double const h = spacing;
		double const b = ( xc - knot( i ) )/h;
		double const a = 1.0 - b;

		/*
		 * S( x ) = a f_i + b f_{i+1} + ( a^3 - a ) h^2 M_i/6 + ( b^3 - b ) h^2 M_{i+1}/6
		 *
		 * on [ x_i, x_{i+1} ], with M = A f, so phi_k reads rows i and i + 1 of
		 * A and the two unit values.
		 */
		double const cA = ( a*a*a - a )*h*h/6.0, cB = ( b*b*b - b )*h*h/6.0;
		double const sA = -( 3.0*a*a - 1.0 )*h/6.0, sB = ( 3.0*b*b - 1.0 )*h/6.0;
		double const *rowA = &curvatureMap[ i*count ];
		double const *rowB = &curvatureMap[ ( i + 1 )*count ];
		for ( std::size_t k = 0; k < count; ++k )
		{
			if ( value )
				value[ k ] = cA*rowA[ k ] + cB*rowB[ k ];
			if ( slope )
				slope[ k ] = outside ? 0.0 : sA*rowA[ k ] + sB*rowB[ k ];
			if ( curvature )
				curvature[ k ] = outside ? 0.0 : a*rowA[ k ] + b*rowB[ k ];
		}
		if ( value )
		{
			value[ i ] += a;
			value[ i + 1 ] += b;
		}
		if ( slope && !outside )
		{
			slope[ i ] -= 1.0/h;
			slope[ i + 1 ] += 1.0/h;
		}
	}

	std::vector<double> CardinalSpline::curvatures( std::vector<double> const &values ) const
	{
		if ( values.size() != count )
			throw std::invalid_argument( "meq::CardinalSpline::curvatures: one value per knot" );
		std::vector<double> out( count, 0.0 );
		for ( std::size_t k = 0; k < count; ++k )
		{
			double total = 0.0;
			for ( std::size_t j = 0; j < count; ++j )
				total += curvatureMap[ k*count + j ]*values[ j ];
			out[ k ] = total;
		}
		return out;
	}

	void CardinalSpline::evaluate( std::vector<double> const &values,
	                               std::vector<double> const &curvatures, double x,
	                               double &value, double &slope, double &curvature ) const
	{
		bool const outside = x < xLower || x > xUpper;
		double const xc = std::min( xUpper, std::max( xLower, x ) );
		std::size_t const i = interval( xc );
		double const h = spacing;
		double const b = ( xc - knot( i ) )/h;
		double const a = 1.0 - b;
		double const mA = curvatures[ i ], mB = curvatures[ i + 1 ];
		value = a*values[ i ] + b*values[ i + 1 ]
		        + ( ( a*a*a - a )*mA + ( b*b*b - b )*mB )*h*h/6.0;
		if ( outside )
		{
			slope = 0.0;
			curvature = 0.0;
			return;
		}
		slope = ( values[ i + 1 ] - values[ i ] )/h
		        + ( -( 3.0*a*a - 1.0 )*mA + ( 3.0*b*b - 1.0 )*mB )*h/6.0;
		curvature = a*mA + b*mB;
	}

	CardinalSplineProfile::CardinalSplineProfile( std::shared_ptr<CardinalSpline const> basis,
	                                              std::vector<double> values )
		: splineBasis( std::move( basis ) ), valueData( std::move( values ) )
	{
		if ( !splineBasis )
			throw std::invalid_argument( "meq::CardinalSplineProfile: the basis must not be null" );
		if ( valueData.size() != splineBasis->size() )
			throw std::invalid_argument( "meq::CardinalSplineProfile: one value per knot" );
		curvatureData = splineBasis->curvatures( valueData );
	}

	double CardinalSplineProfile::operator()( double psi ) const
	{
		double v = 0.0, s = 0.0, c = 0.0;
		splineBasis->evaluate( valueData, curvatureData, psi, v, s, c );
		return v;
	}

	double CardinalSplineProfile::prime( double psi ) const
	{
		double v = 0.0, s = 0.0, c = 0.0;
		splineBasis->evaluate( valueData, curvatureData, psi, v, s, c );
		return s;
	}

	double CardinalSplineProfile::doublePrime( double psi ) const
	{
		double v = 0.0, s = 0.0, c = 0.0;
		splineBasis->evaluate( valueData, curvatureData, psi, v, s, c );
		return c;
	}

	CardinalSpline const & CardinalSplineProfile::basis() const
	{
		return *splineBasis;
	}

	std::vector<double> const & CardinalSplineProfile::values() const
	{
		return valueData;
	}

}
