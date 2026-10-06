#include "meq/RegularisedLeastSquares.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

#include <Eigen/Dense>
#include <Eigen/SVD>

namespace meq
{
	double defaultRankThreshold( int rows, int columns )
	{
		return static_cast<double>( std::max( rows, columns ) )
		       *std::numeric_limits<double>::epsilon();
	}

	RegularisedStep solveRegularisedLeastSquares( std::vector<double> const &matrix,
	                                              int rows, int columns,
	                                              std::vector<double> const &rhs,
	                                              std::vector<double> const &gamma,
	                                              double relativeThreshold )
	{
		auto refuse = []( std::string const &why )
		{
			throw std::invalid_argument( "meq::solveRegularisedLeastSquares: " + why );
		};

		if ( rows < 1 || columns < 1 )
			refuse( "the matrix must have at least one row and one column, got "
			        + std::to_string( rows ) + " x " + std::to_string( columns ) );
		std::size_t const r = static_cast<std::size_t>( rows );
		std::size_t const c = static_cast<std::size_t>( columns );
		if ( matrix.size() != r*c )
			refuse( "the matrix holds " + std::to_string( matrix.size() )
			        + " entries for a " + std::to_string( rows ) + " x "
			        + std::to_string( columns ) + " matrix" );
		if ( rhs.size() != r )
			refuse( "the right-hand side has " + std::to_string( rhs.size() )
			        + " entries against " + std::to_string( rows ) + " rows" );
		if ( !gamma.empty() && gamma.size() != c )
			refuse( "gamma has " + std::to_string( gamma.size() )
			        + " entries against " + std::to_string( columns )
			        + " columns; give one per column, zero where a column is "
			          "not regularised, or none at all" );
		for ( double v : matrix )
			if ( !std::isfinite( v ) )
				refuse( "the matrix has a non-finite entry" );
		for ( double v : rhs )
			if ( !std::isfinite( v ) )
				refuse( "the right-hand side has a non-finite entry" );
		for ( double v : gamma )
			if ( !( v >= 0.0 ) || !std::isfinite( v ) )
				refuse( "every gamma must be finite and non-negative, got "
				        + std::to_string( v ) );

		double const threshold = relativeThreshold > 0.0
		                         ? relativeThreshold
		                         : defaultRankThreshold( rows, columns );

		Eigen::MatrixXd system( rows, columns );
		for ( int i = 0; i < rows; ++i )
			for ( int j = 0; j < columns; ++j )
				system( i, j ) = matrix[ static_cast<std::size_t>( i )*c
				                    + static_cast<std::size_t>( j ) ];
		Eigen::VectorXd f( rows );
		for ( int i = 0; i < rows; ++i )
			f( i ) = rhs[ static_cast<std::size_t>( i ) ];

		RegularisedStep out;

		// THE RANK OF M ALONE, before the regularisation makes every system
		// full rank and hides it.
		{
			Eigen::JacobiSVD<Eigen::MatrixXd> bare( system );
			Eigen::VectorXd const &s = bare.singularValues();
			out.singularValues.assign( s.data(), s.data() + s.size() );
			double const cut = s.size() > 0 ? threshold*s( 0 ) : 0.0;
			for ( Eigen::Index k = 0; k < s.size(); ++k )
				if ( s( k ) > cut )
					++out.rank;
			out.nullity = columns - out.rank;
		}

		// THE AUGMENTED SYSTEM: one extra row per regularised column, so a
		// column with gamma = 0 adds nothing and stays exactly as M has it.
		int extra = 0;
		for ( double v : gamma )
			if ( v > 0.0 )
				++extra;
		Eigen::MatrixXd augmented = Eigen::MatrixXd::Zero( rows + extra, columns );
		Eigen::VectorXd b = Eigen::VectorXd::Zero( rows + extra );
		augmented.topRows( rows ) = system;
		b.head( rows ) = f;
		int next = rows;
		for ( int j = 0; j < static_cast<int>( gamma.size() ); ++j )
			if ( gamma[ static_cast<std::size_t>( j ) ] > 0.0 )
				augmented( next++, j ) = gamma[ static_cast<std::size_t>( j ) ];

		Eigen::JacobiSVD<Eigen::MatrixXd> svd( augmented, Eigen::ComputeThinU | Eigen::ComputeThinV );
		svd.setThreshold( threshold );
		Eigen::VectorXd const x = svd.solve( b );

		out.step.assign( x.data(), x.data() + x.size() );
		out.residualNorm = ( system*x - f ).norm();
		double penalty = 0.0;
		for ( int j = 0; j < static_cast<int>( gamma.size() ); ++j )
		{
			double const term = gamma[ static_cast<std::size_t>( j ) ]*x( j );
			penalty += term*term;
		}
		out.penaltyNorm = std::sqrt( penalty );
		return out;
	}
}
