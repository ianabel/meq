#include "meq/SafetyFactorSolve.hpp"

#include <algorithm>
#include <numeric>
#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>

// Eigen for the step: a column-pivoting QR of a dense ( m + 1 ) x m, and the
// singular values the conditioning is reported in. Header only, as in
// SafetyFactor.cpp.
#include <Eigen/Dense>
#include <Eigen/SVD>

namespace meq
{

	namespace
	{
		double maxAbs( std::vector<double> const &v )
		{
			double out = 0.0;
			for ( double x : v )
				out = std::max( out, std::abs( x ) );
			return out;
		}

		/// The constraint rows' weight in the line search's merit.
		double const constraintWeight = 1.0e6;
	}

	LeastSquaresResult solveLeastSquares(
		std::vector<double> start,
		std::function<bool( std::vector<double> const &, std::vector<double> & )> const &residual,
		std::function<bool( std::vector<double> const &, std::vector<double> & )> const &jacobian,
		LeastSquaresOptions const &options )
	{
		if ( start.empty() )
			throw std::invalid_argument( "meq::solveLeastSquares: there are no unknowns" );
		if ( !residual || !jacobian )
			throw std::invalid_argument( "meq::solveLeastSquares: a residual and its Jacobian are both needed" );
		if ( options.maxIterations < 1 || options.backtracks < 0
		     || !( options.tolerance >= 0.0 ) || !( options.sufficientDecrease > 0.0 ) )
			throw std::invalid_argument( "meq::solveLeastSquares: the options do not describe an iteration" );

		LeastSquaresResult out;
		std::size_t const n = start.size();
		std::vector<double> x = std::move( start );
		std::vector<double> r;

		// Whether the last evaluation was at x: the Jacobian is only ever
		// formed at the state the residual last left behind.
		bool atX = false;
		auto evaluate = [ & ]( std::vector<double> const &at, std::vector<double> &into )
		{
			++out.evaluations;
			bool const admitted = residual( at, into );
			atX = admitted && at == x;
			return admitted;
		};

		auto finish = [ & ]( bool converged, std::string status )
		{
			out.x = x;
			out.residual = r;
			out.misfit = maxAbs( r );
			out.converged = converged;
			out.status = std::move( status );
			return out;
		};

		std::vector<double> measuredFloor;
		bool measured = false;
		/*
		 * M-363's measurement: x evaluated again after its last unknown is
		 * moved each way, so from two other warm starts; four standard
		 * deviations of one draw. Ends with x evaluated, or with atX false if
		 * x was refused on the way back.
		 */
		auto measureNoise = [ & ]()
		{
			measured = true;
			std::vector<double> spread( r.size(), 0.0 );
			for ( double sign : { 1.0, -1.0 } )
			{
				std::vector<double> away( x ), discard, again;
				away[ n - 1 ] += sign*options.noiseProbe*std::max( std::abs( x[ n - 1 ] ), 1.0 );
				evaluate( away, discard );
				if ( !evaluate( x, again ) || again.size() != r.size() )
					continue;
				for ( std::size_t i = 0; i < r.size(); ++i )
					spread[ i ] = std::max( spread[ i ], std::abs( again[ i ] - r[ i ] ) );
			}
			measuredFloor.resize( r.size() );
			for ( std::size_t i = 0; i < r.size(); ++i )
				measuredFloor[ i ] = 4.0*spread[ i ]/std::sqrt( 2.0 );
		};

		try
		{
			if ( !evaluate( x, r ) )
				return finish( false, "the starting point was refused, so there is nothing to step from" );

			for ( int it = 0; it < options.maxIterations; ++it )
			{
				out.iterations = it;
				if ( !atX && !evaluate( x, r ) )
					return finish( false, "the accepted point no longer evaluates" );
				std::vector<double> j;
				if ( !jacobian( x, j ) || j.size() != r.size()*n )
					return finish( false, "the Jacobian could not be formed" );

				Eigen::Index const rows = static_cast<Eigen::Index>( r.size() );
				Eigen::Index const cols = static_cast<Eigen::Index>( n );
				Eigen::MatrixXd a( rows, cols );
				for ( Eigen::Index i = 0; i < rows; ++i )
					for ( Eigen::Index k = 0; k < cols; ++k )
						a( i, k ) = j[ static_cast<std::size_t>( i*cols + k ) ];
				Eigen::VectorXd b( rows );
				for ( Eigen::Index i = 0; i < rows; ++i )
					b( i ) = -r[ static_cast<std::size_t>( i ) ];

				Eigen::JacobiSVD<Eigen::MatrixXd> const svd( a );
				double const largest = svd.singularValues()( 0 );
				double const ratio = largest > 0.0
				                     ? svd.singularValues()( svd.singularValues().size() - 1 )/largest : 0.0;
				out.conditioning = std::min( out.conditioning, ratio );
				if ( ratio < static_cast<double>( n )*std::numeric_limits<double>::epsilon() )
				{
					std::ostringstream message;
					message << "the Jacobian is RANK DEFICIENT (sigma_min/sigma_max " << ratio
					        << "): the rows do not determine the unknowns";
					return finish( false, message.str() );
				}

				/*
				 * THE STEP. Without constraints, min || J dx + r || by QR.
				 *
				 * With the last `c` rows equality constraints and a level
				 * functional l: the least squares of the other rows WITHIN a
				 * level set of l, the constraints choosing the level. Written
				 * dx = Z y + t l^ (Z an orthonormal basis of l's null space,
				 * l^ = l/|l|), the rows' least squares in y is y0 + t y1 by
				 * one QR of A Z, and t solves the linearised constraint. So the
				 * curve the answers lie on -- least squares at each level of
				 * l -- is the same whichever constraint picks the point on it.
				 */
				Eigen::VectorXd step;
				// An orthonormal basis of what the step can remove from the
				// least-squares rows: range( J ), or range( J Z ) with a
				// constraint. The search's merit is measured in it, below.
				Eigen::MatrixXd removableBasis;
				std::size_t const constraints = options.constraints;
				if ( constraints == 0 )
				{
					Eigen::ColPivHouseholderQR<Eigen::MatrixXd> const qr( a );
					step = qr.solve( b );
					removableBasis = qr.householderQ()*Eigen::MatrixXd::Identity( rows, qr.rank() );
				}
				else
				{
					if ( constraints != 1 || !options.level || r.size() <= 1 )
						throw std::invalid_argument( "meq::solveLeastSquares: one constraint row and a "
						                             "level functional, or none" );
					std::vector<double> const level = options.level( x );
					if ( level.size() != n )
						throw std::invalid_argument( "meq::solveLeastSquares: the level functional's "
						                             "gradient has one entry per unknown" );
					Eigen::VectorXd l( cols );
					for ( Eigen::Index k = 0; k < cols; ++k )
						l( k ) = level[ static_cast<std::size_t>( k ) ];
					l.normalize();
					// Z: the last n - 1 columns of a Householder reflection
					// taking l to a multiple of e_1.
					Eigen::HouseholderQR<Eigen::MatrixXd> const reflect( l );
					Eigen::MatrixXd const q = reflect.householderQ()*Eigen::MatrixXd::Identity( cols, cols );
					Eigen::MatrixXd const z = q.rightCols( cols - 1 );
					Eigen::Index const rowsFree = rows - 1;
					Eigen::MatrixXd const az = a.topRows( rowsFree )*z;
					Eigen::ColPivHouseholderQR<Eigen::MatrixXd> const qr( az );
					removableBasis = qr.householderQ()*Eigen::MatrixXd::Identity( rowsFree, qr.rank() );
					Eigen::VectorXd const y0 = qr.solve( b.head( rowsFree ) );
					Eigen::VectorXd const y1 = qr.solve( -( a.topRows( rowsFree )*l ) );
					Eigen::RowVectorXd const c = a.row( rows - 1 );
					double const denominator = c.dot( z*y1 ) + c.dot( l );
					if ( !( std::abs( denominator ) > 0.0 ) )
						return finish( false, "the constraint does not move the level: it cannot fix it" );
					double const t = ( b( rows - 1 ) - c.dot( z*y0 ) )/denominator;
					step = z*( y0 + t*y1 ) + t*l;
				}
				Eigen::VectorXd const change = a*step;
				out.removable = change.lpNorm<Eigen::Infinity>();
				out.history.push_back( out.removable );

				/*
				 * THE FLOOR PER ROW: the row's own noise, modelled or measured,
				 * and what holding x in real_t costs it. A row is at its
				 * threshold when its removable part is under the larger of the
				 * tolerance and its floor; converged is every row at it.
				 */
				std::vector<double> model( r.size(), 0.0 );
				if ( options.floor )
				{
					std::vector<double> const given = options.floor();
					for ( std::size_t i = 0; i < std::min( given.size(), model.size() ); ++i )
						model[ i ] = given[ i ];
				}
				double const size = maxAbs( x );
				auto floorOf = [ & ]( std::size_t i )
				{
					double const ofRow = a.row( static_cast<Eigen::Index>( i ) ).lpNorm<1>();
					double const noise = measuredFloor.empty() ? model[ i ]
					                     : std::max( model[ i ], measuredFloor[ i ] );
					return noise + ofRow*options.epsilon*size;
				};
				// The worst row's removable over its threshold; at or under one
				// is converged. Within `within` is near enough to measure.
				auto worst = [ & ]( double within )
				{
					double ratio = 0.0;
					out.floor = 0.0;
					for ( std::size_t i = 0; i < r.size(); ++i )
					{
						double const floorRow = floorOf( i );
						out.floor = std::max( out.floor, floorRow );
						ratio = std::max( ratio, std::abs( change( static_cast<Eigen::Index>( i ) ) )
						                         /std::max( options.tolerance, within*floorRow ) );
					}
					return ratio;
				};
				bool const floorJudges =
					std::any_of( model.begin(), model.end(),
					             [ & ]( double f ) { return f > options.tolerance; } );

				if ( worst( 1.0 ) <= 1.0 )
					return finish( true, "converged: the removable residual is at the threshold" );

				// WHERE THE FLOOR JUDGES, MEASURE IT ON THE WAY IN (M-363).
				if ( !measured && floorJudges && worst( 10.0 ) <= 1.0 )
				{
					measureNoise();
					if ( worst( 1.0 ) <= 1.0 )
					{
						if ( !atX && !evaluate( x, r ) )
							return finish( false, "the accepted point no longer evaluates" );
						return finish( true, "converged at the measured noise" );
					}
					if ( !atX && !evaluate( x, r ) )
						return finish( false, "the accepted point no longer evaluates" );
				}

				/*
				 * THE SEARCH: Armijo, a refusal a rejection, on the part of the
				 * residual the step can remove -- the rows projected on
				 * range( J Z ) at x, and the constraint rows weighted:
				 *
				 *     phi( u ) = || Q^T r_free( u ) ||^2 + w || r_con( u ) ||^2.
				 *
				 * NOT ON || r ||^2, AND THAT IS MEASURED (M-QS-2). With a
				 * constraint the answer is least squares within a level set
				 * of `level`, the constraint picking the level (M-396). It is
				 * not a minimum of || r ||^2 on the constraint's surface, so
				 * near it a right step can raise || r ||^2 to first order.
				 * On the limited box at level 1 every trial of such a step
				 * rose linearly in its length, and the loop stalled at a
				 * removable 8e-7, on one gg' in three perturbed by 1e-14.
				 * phi is zero exactly at the answer, and the step takes its
				 * model to zero, so the step is a descent direction of phi with
				 * slope -2 phi. Without a constraint Z is the identity, and
				 * the misfit the step cannot remove leaves the comparison
				 * along with its noise.
				 */
				std::size_t const free = r.size() - options.constraints;
				auto merit = [ & ]( std::vector<double> const &v )
				{
					Eigen::VectorXd freeRows( static_cast<Eigen::Index>( free ) );
					for ( std::size_t i = 0; i < free; ++i )
						freeRows( static_cast<Eigen::Index>( i ) ) = v[ i ];
					double sum = ( removableBasis.transpose()*freeRows ).squaredNorm();
					for ( std::size_t i = free; i < v.size(); ++i )
						sum += constraintWeight*v[ i ]*v[ i ];
					return sum;
				};
				double const f0 = merit( r );
				double predicted = 0.0;
				{
					std::vector<double> modelled( r.size() );
					for ( std::size_t i = 0; i < r.size(); ++i )
						modelled[ i ] = r[ i ] + change( static_cast<Eigen::Index>( i ) );
					predicted = std::max( f0 - merit( modelled ), 0.0 );
				}
				double alpha = 1.0;
				bool accepted = false;
				std::vector<double> trial( n ), rt;
				for ( int back = 0; back <= options.backtracks; ++back )
				{
					for ( std::size_t k = 0; k < n; ++k )
						trial[ k ] = x[ k ] + alpha*step( static_cast<Eigen::Index>( k ) );
					if ( evaluate( trial, rt ) && rt.size() == r.size()
					     && merit( rt ) <= f0 - options.sufficientDecrease*alpha*predicted )
					{
						accepted = true;
						break;
					}
					++out.rejected;
					alpha *= 0.5;
				}
				if ( !accepted )
				{
					// A STALL IS MEASURED BEFORE IT IS JUDGED.
					if ( !atX && !evaluate( x, r ) )
						return finish( false, "the accepted point no longer evaluates" );
					if ( !measured )
						measureNoise();
					if ( !atX )
						evaluate( x, r );
					if ( worst( 1.0 ) <= 1.0 )
						return finish( true, "converged: the search stopped at the noise" );
					std::ostringstream message;
					message << "did not converge: the line search found no sufficient decrease, "
					        << "the removable residual " << out.removable << " above the threshold "
					        << std::max( options.tolerance, out.floor );
					return finish( false, message.str() );
				}
				x = trial;
				r = rt;
				atX = true;
			}
			out.iterations = options.maxIterations;
			return finish( false, "did not converge within the iteration cap" );
		}
		catch ( std::exception const &error )
		{
			return finish( false, std::string( "the residual or its Jacobian threw: " ) + error.what() );
		}
	}

}
