#include "Profiles.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <istream>
#include <limits>
#include <ostream>
#include <sstream>
#include <stdexcept>

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

	void SplineProfile::write( std::ostream & os ) const
	{
		// max_digits10 so that write() followed by fromStream() is the identity.
		std::streamsize const oldPrecision = os.precision( std::numeric_limits<double>::max_digits10 );

		os << "# meq spline profile: " << knotData.size() << " knots\n";
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

}
