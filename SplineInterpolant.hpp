#ifndef SPLINEINTERPOLANT_HPP
#define SPLINEINTERPOLANT_HPP

#include <functional>
#include <vector>
#include <string>
#include <iostream>


/*
 * Utility class for interpolating functions on [ 0, 1 ] 
 * accurately and efficiently, using cubic splines
 */

class HermiteCubicSpline {
	private:
		double x_l,x_u;
		double delta;
		double f_x_l,f_x_u;
		double delta_f;
		double fprime_l,fprime_u;

		// Static copies for the basis splines and their derivatives
		static double h00( double t )
		{
			return ( 1. + 2*t )*( 1 - t )*( 1 - t );
		};
		static double h10( double t )
		{
			return ( t )*( 1 - t )*( 1 - t );
		};
		static double h01( double t )
		{
			return t*t*( 3.0 - 2.0*t );
		};
		static double h11( double t )
		{
			return t*t*( t - 1.0 );
		};
		
		static double h00_p( double t )
		{
			return 6.0*t*( t - 1.0 );
		};
		static double h10_p( double t )
		{
			return 3.0*t*t - 4.0*t + 1.0;
		};
		static double h01_p( double t )
		{
			return 6.0*t*( 1.0 - t );
		};
		static double h11_p( double t )
		{
			return t * ( 3*t - 2 );
		};
	public:
		HermiteCubicSpline( double l, double u, double fl, double fu, double fprime_l_in, double fprime_u_in )
		{
			if ( !( l < u ) )
				throw std::logic_error( "Ftagn" );

			x_l = l;
			x_u = u;
			delta = x_u - x_l;

			f_x_l = fl;
			f_x_u = fu;
			delta_f = f_x_u - f_x_l;

			fprime_l = fprime_l_in;
			fprime_u = fprime_u_in;
		};

		double operator()( double x ) const
		{
			if ( x < x_l )
				return f_x_l;
			if ( x > x_u )
				return f_x_u;

			double t = ( x - x_l )/delta;

			return f_x_l * h00( t ) + fprime_l * delta * h10( t ) + f_x_u * h01( t ) + fprime_u * delta * h11( t );
		};

		double prime( double x ) const
		{
			if ( x < x_l )
				return f_x_l;
			if ( x > x_u )
				return f_x_u;

			double t = ( x - x_l )/delta;

			return ( f_x_l * h00_p( t ) + fprime_l * delta * h10_p( t ) + f_x_u * h01_p( t ) + fprime_u * delta * h11_p( t ) ) / delta;
		};

		std::pair<double,double> getInterval() const
		{
			return std::make_pair( x_l, x_u );
		};

		std::pair<double,double> getValues() const
		{
			return std::make_pair( f_x_l, f_x_u );
		};

		std::pair<double,double> getDerivatives() const
		{
			return std::make_pair( fprime_l, fprime_u );
		};
};

class SplineInterpolant {
	public:
		using RealFunc = std::function<double( double )>;
		using Spline = HermiteCubicSpline;
		// SplineInterpolant( RealFunc F, unsigned int N ); // Split [0,1] into N intervals and spline away
		// If we know F' then we can do everything directly
		SplineInterpolant( RealFunc F, RealFunc Fprime, unsigned int N ) {
			piecewise_data.clear();
			double delta_x = 1.0/N;
			for ( unsigned int i=0; i < N; ++i )
			{
				double x_l = i*delta_x;
				double x_u = ( i + 1 )*delta_x;
				piecewise_data.emplace_back( x_l, x_u, F( x_l ), F( x_u ), Fprime( x_l ), Fprime( x_u ) );
			}
		};
		// SplineInterpolant( std::vector< std::pair< double, double > > ); // Construct from function values
		// Construct from function values and derivatives
		SplineInterpolant( std::vector< std::tuple< double, double, double > > data )
		{
			for ( size_t i=0; i < data.size()-1; i++ )
			{
				auto [ x_l, f_x_l, fprim_l ] = data[ i ];
				auto [ x_u, f_x_u, fprim_u ] = data[ i + 1 ];
				piecewise_data.emplace_back( x_l, x_u, f_x_l, f_x_u, fprim_l, fprim_u );
			}
		};

		// Read from stream.
		SplineInterpolant( std::istream &is )
		{
			// Loop over lines
			char c;
			std::vector< std::tuple<double,double,double> > data;
			data.clear();
			while ( is.good() ) {
				c = is.peek();
				if ( c == '\n' )
					break;
				if ( c == std::char_traits<char>::eof() )
					break;
				if ( c == '#' )
				{
					std::string line;
					std::getline( is, line );
					continue;
				}
				else
				{
					double x,f_x,f_prime_x;
					is >> x >> f_x >> f_prime_x;
					std::cerr << "Read x=" << x << " f(x)="<<f_x << " f'(x)=" << f_prime_x << std::endl;
					// Eat the rest of the line.
					std::string line;
					std::getline( is, line );

					data.emplace_back( x, f_x, f_prime_x );
				}
			}
			if ( data.size() <= 1 )
				throw std::logic_error( "Unable to construct interpolant from input data" );
			else
			{
				piecewise_data.clear();
				for ( size_t i=0; i < data.size()-1; i++ )
				{
					auto [ x_l, f_x_l, fprim_l ] = data[ i ];
					auto [ x_u, f_x_u, fprim_u ] = data[ i + 1 ];
					piecewise_data.emplace_back( x_l, x_u, f_x_l, f_x_u, fprim_l, fprim_u );
				}
			}
		}

		SplineInterpolant( SplineInterpolant const& other )
			: piecewise_data( other.piecewise_data )
		{
		}; // copy
		
		SplineInterpolant( SplineInterpolant && other )
			: piecewise_data( other.piecewise_data )
		{
		}; // move

		SplineInterpolant & operator=( SplineInterpolant const& other )
		{
			piecewise_data = other.piecewise_data;
			return *this;
		}; // copy

		void Write( std::ostream &os ) const
		{
			os << "# Spline Interpolation data" << std::endl;
			os << "# x\t f(x)\t f'(x)" << std::endl;
			for ( auto const& spl : piecewise_data )
			{
				double x = spl.getInterval().first;
				double f_x = spl.getValues().first;
				double fprime_x = spl.getDerivatives().first;
				os << x << "\t" << f_x << "\t" << fprime_x << std::endl;
			}

			{
				auto const & spl = piecewise_data.back();
				double x = spl.getInterval().second;
				double f_x = spl.getValues().second;
				double fprime_x = spl.getDerivatives().second;
				os << x << "\t" << f_x << "\t" << fprime_x << std::endl;
			}
			// Blank line marks end-of-data
			os << std::endl;
		};

		friend std::ostream & operator<<( std::ostream &os, SplineInterpolant const& spl )
		{
			spl.Write( os );
			return os;
		};

		double operator()( double s ) {
			if ( s < 0.0 || s > 1.0 )
				throw std::invalid_argument( "This class only interpolates on [0,1]" );
			Spline const& spl = FindInterval( piecewise_data, s );
			return spl( s );
		};
		double prime( double s ) {
			if ( s < 0.0 || s > 1.0 )
				throw std::invalid_argument( "This class only interpolates on [0,1]" );
			Spline const& spl = FindInterval( piecewise_data, s );
			return spl.prime( s );
		};
	private:
		static Spline const& FindInterval( std::vector<Spline> partition, double x )
		{
			unsigned int lower = 0;
			unsigned int upper = partition.size();

			if ( ( x >= partition[ lower ].getInterval().first ) && ( x <= partition[ lower ].getInterval().second ) )
				return partition[ lower ];


			do {
				unsigned int mid = ( lower + upper )/2;
				if ( x < partition[ mid ].getInterval().first )
					upper = mid;
				else if ( x > partition[ mid ].getInterval().second )
					lower = mid;
				else if ( ( x >= partition[ mid ].getInterval().first ) && ( x <= partition[ mid ].getInterval().second ) )
					return partition[ mid ];

			} while ( upper - lower > 1 );

			throw std::logic_error( "Needle not in partition" );
			return partition[ 0 ];
		};
		std::vector<Spline> piecewise_data;

};
#endif // SPLINEINTERPOLANT_HPP
