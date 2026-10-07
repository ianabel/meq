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

	ConstrainedStep solveConstrainedLeastSquares( std::vector<double> const &constraints,
	                                              int constraintRows,
	                                              std::vector<double> const &right,
	                                              std::vector<double> const &matrix,
	                                              int rows, int columns,
	                                              std::vector<double> const &rhs,
	                                              std::vector<double> const &gamma,
	                                              double relativeThreshold )
	{
		auto refuse = []( std::string const &why )
		{
			throw std::invalid_argument( "meq::solveConstrainedLeastSquares: " + why );
		};
		if ( columns < 1 || rows < 0 || constraintRows < 0 )
			refuse( "the shape " + std::to_string( constraintRows ) + " + "
			        + std::to_string( rows ) + " rows by " + std::to_string( columns )
			        + " columns is not a problem" );
		std::size_t const c = static_cast<std::size_t>( columns );
		if ( constraints.size() != static_cast<std::size_t>( constraintRows )*c )
			refuse( "the constraint matrix holds " + std::to_string( constraints.size() )
			        + " entries for " + std::to_string( constraintRows ) + " x "
			        + std::to_string( columns ) );
		if ( right.size() != static_cast<std::size_t>( constraintRows ) )
			refuse( "the constraints' right-hand side has " + std::to_string( right.size() )
			        + " entries against " + std::to_string( constraintRows ) + " rows" );
		if ( matrix.size() != static_cast<std::size_t>( rows )*c )
			refuse( "the matrix holds " + std::to_string( matrix.size() )
			        + " entries for " + std::to_string( rows ) + " x "
			        + std::to_string( columns ) );
		if ( rhs.size() != static_cast<std::size_t>( rows ) )
			refuse( "the right-hand side has " + std::to_string( rhs.size() )
			        + " entries against " + std::to_string( rows ) + " rows" );
		if ( !gamma.empty() && gamma.size() != c )
			refuse( "gamma has " + std::to_string( gamma.size() ) + " entries against "
			        + std::to_string( columns ) + " columns" );
		for ( std::vector<double> const *list : { &constraints, &right, &matrix, &rhs } )
			for ( double v : *list )
				if ( !std::isfinite( v ) )
					refuse( "an entry is not finite" );
		for ( double v : gamma )
			if ( !( v >= 0.0 ) || !std::isfinite( v ) )
				refuse( "every gamma must be finite and non-negative, got "
				        + std::to_string( v ) );

		auto load = [ columns ]( std::vector<double> const &values, int count )
		{
			Eigen::MatrixXd out( count, columns );
			for ( int i = 0; i < count; ++i )
				for ( int j = 0; j < columns; ++j )
					out( i, j ) = values[ static_cast<std::size_t>( i )
					                      *static_cast<std::size_t>( columns )
					                      + static_cast<std::size_t>( j ) ];
			return out;
		};
		Eigen::MatrixXd const bound = load( constraints, constraintRows );
		Eigen::MatrixXd const fitted = load( matrix, rows );
		Eigen::VectorXd d = Eigen::VectorXd::Zero( constraintRows );
		for ( int i = 0; i < constraintRows; ++i )
			d( i ) = right[ static_cast<std::size_t>( i ) ];
		Eigen::VectorXd b = Eigen::VectorXd::Zero( rows );
		for ( int i = 0; i < rows; ++i )
			b( i ) = rhs[ static_cast<std::size_t>( i ) ];
		Eigen::VectorXd g = Eigen::VectorXd::Zero( columns );
		for ( int j = 0; j < static_cast<int>( gamma.size() ); ++j )
			g( j ) = gamma[ static_cast<std::size_t>( j ) ];

		// THE CONSTRAINTS: a particular solution and the null space, from one
		// decomposition with the full right singular vectors.
		Eigen::VectorXd particular = Eigen::VectorXd::Zero( columns );
		Eigen::MatrixXd nullBasis = Eigen::MatrixXd::Identity( columns, columns );
		if ( constraintRows > 0 )
		{
			double const threshold = relativeThreshold > 0.0
			                         ? relativeThreshold
			                         : defaultRankThreshold( constraintRows, columns );
			Eigen::JacobiSVD<Eigen::MatrixXd> svd( bound, Eigen::ComputeThinU
			                                              | Eigen::ComputeFullV );
			svd.setThreshold( threshold );
			particular = svd.solve( d );
			int const boundRank = static_cast<int>( svd.rank() );
			nullBasis = svd.matrixV().rightCols( columns - boundRank );
		}
		int const openDirections = static_cast<int>( nullBasis.cols() );

		ConstrainedStep out;
		Eigen::VectorXd x = particular;
		if ( openDirections > 0 && rows > 0 )
		{
			Eigen::MatrixXd const reduced = fitted*nullBasis;
			double const threshold = relativeThreshold > 0.0
			                         ? relativeThreshold
			                         : defaultRankThreshold( rows, openDirections );
			{
				Eigen::JacobiSVD<Eigen::MatrixXd> bare( reduced );
				Eigen::VectorXd const &s = bare.singularValues();
				out.singularValues.assign( s.data(), s.data() + s.size() );
				double const cut = s.size() > 0 ? threshold*s( 0 ) : 0.0;
				for ( Eigen::Index k = 0; k < s.size(); ++k )
					if ( s( k ) > cut )
						++out.rank;
			}
			out.nullity = openDirections - out.rank;

			// [ A Z ; diag( gamma ) Z ] y ~ [ b - A x_p ; -diag( gamma ) x_p ],
			// so the penalty is on the whole step x and not on y.
			Eigen::MatrixXd augmented( rows + columns, openDirections );
			augmented.topRows( rows ) = reduced;
			augmented.bottomRows( columns ) = g.asDiagonal()*nullBasis;
			Eigen::VectorXd target( rows + columns );
			target.head( rows ) = b - fitted*particular;
			target.tail( columns ) = -( g.asDiagonal()*particular );
			Eigen::JacobiSVD<Eigen::MatrixXd> svd( augmented, Eigen::ComputeThinU
			                                                  | Eigen::ComputeThinV );
			svd.setThreshold( threshold );
			x += nullBasis*svd.solve( target );
		}
		else
			out.nullity = openDirections;

		out.step.assign( x.data(), x.data() + x.size() );
		out.residualNorm = rows > 0 ? ( fitted*x - b ).norm() : 0.0;
		out.constraintResidualNorm =
			constraintRows > 0 ? ( bound*x - d ).norm() : 0.0;
		out.penaltyNorm = ( g.asDiagonal()*x ).norm();
		return out;
	}
}
