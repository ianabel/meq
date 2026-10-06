/*
 * meq::solveRegularisedLeastSquares -- INVERSE-AND-RECONSTRUCTION-PLAN.md IV-2.
 *
 * The inverse solve's dense corner, tested in isolation and against closed
 * forms, before any solver is attached to it. Every case builds its matrix as
 * U diag( s ) V^T from orthonormal U and V the test makes itself, so the
 * minimum-norm answer is V_r diag( 1/s_r ) U_r^T f WITHOUT reference to the code
 * under test, and the rank is the number of non-zero s by construction.
 *
 * WHAT EACH CASE WOULD CATCH, since the plan asks for that to be checked rather
 * than assumed: a TRANSPOSED read of the row-major matrix fails the closed form
 * (the matrices are not square, and not symmetric where they are); a gamma
 * applied to the WRONG COLUMN fails the normal-equations identity and the
 * damping case; and a solve that is not the minimum-norm one fails the
 * rank-deficient case, which has a whole null space to wander in.
 */
#define BOOST_TEST_MODULE MeqRegularisedLeastSquaresTests
#include <boost/test/unit_test.hpp>

#include "RegularisedLeastSquares.hpp"

#include <cmath>
#include <cstdio>
#include <limits>
#include <random>
#include <stdexcept>
#include <vector>

namespace
{
	/// A dense matrix, row-major, for the fixtures alone.
	struct Dense
	{
		int rows = 0;
		int columns = 0;
		std::vector<double> a;

		Dense( int r, int c ) : rows( r ), columns( c ),
		                        a( static_cast<std::size_t>( r*c ), 0.0 ) {}
		double &operator()( int i, int j ) { return a[ static_cast<std::size_t>( i*columns + j ) ]; }
		double operator()( int i, int j ) const { return a[ static_cast<std::size_t>( i*columns + j ) ]; }
	};

	/// n orthonormal columns of length n, by modified Gram-Schmidt run twice
	/// over random vectors -- twice, because once loses orthogonality at the
	/// 1e-13 level this file asserts below.
	Dense orthonormal( int n, unsigned seed )
	{
		std::mt19937 generator( seed );
		std::normal_distribution<double> normal;
		Dense q( n, n );
		for ( int i = 0; i < n; ++i )
			for ( int j = 0; j < n; ++j )
				q( i, j ) = normal( generator );
		for ( int pass = 0; pass < 2; ++pass )
			for ( int j = 0; j < n; ++j )
			{
				for ( int k = 0; k < j; ++k )
				{
					double dot = 0.0;
					for ( int i = 0; i < n; ++i )
						dot += q( i, k )*q( i, j );
					for ( int i = 0; i < n; ++i )
						q( i, j ) -= dot*q( i, k );
				}
				double norm = 0.0;
				for ( int i = 0; i < n; ++i )
					norm += q( i, j )*q( i, j );
				norm = std::sqrt( norm );
				for ( int i = 0; i < n; ++i )
					q( i, j ) /= norm;
			}
		return q;
	}

	/// The fixture: M = U_r diag( s ) V_r^T, f random, and the closed-form
	/// minimum-norm least-squares solution.
	struct Known
	{
		Dense m;
		std::vector<double> f;
		std::vector<double> exact;
		int rank;
	};

	Known known( int rows, int columns, std::vector<double> const &s, unsigned seed )
	{
		int const r = static_cast<int>( s.size() );
		Dense const u = orthonormal( rows, seed );
		Dense const v = orthonormal( columns, seed + 1 );
		Known k{ Dense( rows, columns ), std::vector<double>( static_cast<std::size_t>( rows ) ),
		         std::vector<double>( static_cast<std::size_t>( columns ), 0.0 ), r };
		for ( int i = 0; i < rows; ++i )
			for ( int j = 0; j < columns; ++j )
				for ( int q = 0; q < r; ++q )
					k.m( i, j ) += u( i, q )*s[ static_cast<std::size_t>( q ) ]*v( j, q );

		std::mt19937 generator( seed + 2 );
		std::normal_distribution<double> normal;
		for ( double &x : k.f )
			x = normal( generator );

		for ( int q = 0; q < r; ++q )
		{
			double c = 0.0;
			for ( int i = 0; i < rows; ++i )
				c += u( i, q )*k.f[ static_cast<std::size_t>( i ) ];
			c /= s[ static_cast<std::size_t>( q ) ];
			for ( int j = 0; j < columns; ++j )
				k.exact[ static_cast<std::size_t>( j ) ] += c*v( j, q );
		}
		return k;
	}

	double relativeError( std::vector<double> const &got, std::vector<double> const &want )
	{
		double d = 0.0, w = 0.0;
		for ( std::size_t i = 0; i < want.size(); ++i )
		{
			d += ( got[ i ] - want[ i ] )*( got[ i ] - want[ i ] );
			w += want[ i ]*want[ i ];
		}
		return std::sqrt( d/w );
	}
}

/*
 * THE UNDERDETERMINED, RANK-DEFICIENT CORNER, WHICH IS THE NORMAL CASE.
 *
 * DIII-D's shape: 21 rows against 33 columns, and of rank 17 -- so the targets
 * leave sixteen directions of the unknowns free, and only the minimum-norm
 * minimiser is a defined answer. Singular values from 1 to 1/4, so the closed
 * form is good to a few eps and 1e-14 is a real bound.
 */
BOOST_AUTO_TEST_CASE( theStepIsTheMinimumNormSolutionAndTheRankIsExact )
{
	std::vector<double> s;
	for ( int q = 0; q < 17; ++q )
		s.push_back( 1.0 - 0.75*q/16.0 );
	Known const k = known( 21, 33, s, 11 );

	meq::RegularisedStep const got =
		meq::solveRegularisedLeastSquares( k.m.a, 21, 33, k.f );

	BOOST_TEST( got.rank == 17 );
	BOOST_TEST( got.nullity == 16 );
	BOOST_TEST_REQUIRE( got.step.size() == 33u );
	double const error = relativeError( got.step, k.exact );
	std::printf( "  21 x 33 of rank 17: minimum-norm step to %.2e, rank %d\n",
	             error, got.rank );
	BOOST_TEST( error < 1.0e-14 );

	// The singular values are reported largest first and the zero ones are
	// zero to round-off.
	BOOST_TEST_REQUIRE( got.singularValues.size() == 21u );
	BOOST_TEST( std::abs( got.singularValues.front() - 1.0 ) < 1.0e-14 );
	BOOST_TEST( std::abs( got.singularValues[ 16 ] - 0.25 ) < 1.0e-14 );
	BOOST_TEST( got.singularValues[ 17 ] < 1.0e-14 );

	// f has a component outside M's range, so the residual is that component
	// and not zero: this is the best fit, and its size is reported.
	BOOST_TEST( got.residualNorm > 0.1 );
	BOOST_TEST( got.penaltyNorm == 0.0 );
}

/*
 * THE OVERDETERMINED CORNER: MACHINE A, FOUR COILS AGAINST FIVE TARGET ROWS.
 * Full column rank, so the minimiser is unique and the normal equations hold;
 * a square solve could not have answered it at all.
 */
BOOST_AUTO_TEST_CASE( anOverdeterminedCornerGetsItsBestFit )
{
	Known const k = known( 5, 4, { 2.0, 1.5, 1.0, 0.5 }, 23 );
	meq::RegularisedStep const got =
		meq::solveRegularisedLeastSquares( k.m.a, 5, 4, k.f );
	BOOST_TEST( got.rank == 4 );
	BOOST_TEST( got.nullity == 0 );
	BOOST_TEST( relativeError( got.step, k.exact ) < 1.0e-14 );

	// M^T ( M x - f ) = 0 at the minimiser.
	for ( int j = 0; j < 4; ++j )
	{
		double g = 0.0;
		for ( int i = 0; i < 5; ++i )
		{
			double mx = 0.0;
			for ( int q = 0; q < 4; ++q )
				mx += k.m( i, q )*got.step[ static_cast<std::size_t>( q ) ];
			g += k.m( i, j )*( mx - k.f[ static_cast<std::size_t>( i ) ] );
		}
		BOOST_TEST( std::abs( g ) < 1.0e-14 );
	}
}

/*
 * THE REGULARISED MINIMISER SATISFIES ( M^T M + diag( gamma )^2 ) x = M^T f.
 *
 * That is freegs4e's x = ( A^T A + gamma^2 I )^-1 A^T b with gamma per column,
 * and checking it -- rather than computing it -- is how a gamma applied to the
 * wrong column, or a dropped transpose, fails: both satisfy SOME least-squares
 * identity, and not this one. Gamma on three of eight columns only, as a
 * current block beside border unknowns.
 */
BOOST_AUTO_TEST_CASE( theRegularisedStepSatisfiesTheNormalEquations )
{
	Known const k = known( 6, 8, { 1.0, 0.8, 0.6, 0.4, 0.2, 0.1 }, 37 );
	std::vector<double> gamma( 8, 0.0 );
	gamma[ 2 ] = 0.3;
	gamma[ 5 ] = 1.7;
	gamma[ 7 ] = 0.05;

	meq::RegularisedStep const got =
		meq::solveRegularisedLeastSquares( k.m.a, 6, 8, k.f, gamma );

	double worst = 0.0;
	for ( int j = 0; j < 8; ++j )
	{
		double lhs = gamma[ static_cast<std::size_t>( j ) ]*gamma[ static_cast<std::size_t>( j ) ]
		             *got.step[ static_cast<std::size_t>( j ) ];
		double rhs = 0.0;
		for ( int i = 0; i < 6; ++i )
		{
			double mx = 0.0;
			for ( int q = 0; q < 8; ++q )
				mx += k.m( i, q )*got.step[ static_cast<std::size_t>( q ) ];
			lhs += k.m( i, j )*mx;
			rhs += k.m( i, j )*k.f[ static_cast<std::size_t>( i ) ];
		}
		worst = std::max( worst, std::abs( lhs - rhs ) );
	}
	std::printf( "  normal equations with gamma on three of eight columns: %.2e\n", worst );
	BOOST_TEST( worst < 1.0e-14 );

	double penalty = 0.0;
	for ( int j = 0; j < 8; ++j )
	{
		double const t = gamma[ static_cast<std::size_t>( j ) ]*got.step[ static_cast<std::size_t>( j ) ];
		penalty += t*t;
	}
	BOOST_TEST( std::abs( got.penaltyNorm - std::sqrt( penalty ) ) < 1.0e-15 );
}

/*
 * gamma -> infinity DAMPS THE CURRENT BLOCK TO NOTHING, MONOTONICALLY, AND
 * LEAVES THE REST SOLVING WHAT IS LEFT.
 *
 * Section 4.4's continuation knob: large gamma is the forward problem, since
 * the currents do not move. The damped block has to fall monotonically to zero;
 * the undamped columns have to converge to the least-squares solution of M with
 * the damped columns deleted, which is what "the rest solve what is left" means.
 */
BOOST_AUTO_TEST_CASE( aLargeGammaDampsTheCurrentBlockToNothing )
{
	Known const k = known( 7, 6, { 1.5, 1.2, 1.0, 0.7, 0.5, 0.3 }, 53 );
	std::vector<int> const current = { 1, 4 };

	// The limit: the least squares over the four undamped columns alone.
	Dense reduced( 7, 4 );
	std::vector<int> kept;
	for ( int j = 0; j < 6; ++j )
		if ( j != 1 && j != 4 )
			kept.push_back( j );
	for ( int i = 0; i < 7; ++i )
		for ( int q = 0; q < 4; ++q )
			reduced( i, q ) = k.m( i, kept[ static_cast<std::size_t>( q ) ] );
	meq::RegularisedStep const limit =
		meq::solveRegularisedLeastSquares( reduced.a, 7, 4, k.f );

	double previous = std::numeric_limits<double>::infinity();
	double lastBlock = 0.0;
	double lastRest = 0.0;
	for ( double g : { 1.0e-3, 1.0e-2, 1.0e-1, 1.0, 1.0e1, 1.0e2, 1.0e3, 1.0e6 } )
	{
		std::vector<double> gamma( 6, 0.0 );
		for ( int j : current )
			gamma[ static_cast<std::size_t>( j ) ] = g;
		meq::RegularisedStep const got =
			meq::solveRegularisedLeastSquares( k.m.a, 7, 6, k.f, gamma );

		double block = 0.0;
		for ( int j : current )
			block += got.step[ static_cast<std::size_t>( j ) ]*got.step[ static_cast<std::size_t>( j ) ];
		block = std::sqrt( block );
		BOOST_TEST( block < previous, "the current block must shrink as gamma grows, at gamma = " << g );
		previous = block;
		lastBlock = block;

		std::vector<double> rest;
		for ( int j : kept )
			rest.push_back( got.step[ static_cast<std::size_t>( j ) ] );
		lastRest = relativeError( rest, limit.step );
	}
	std::printf( "  gamma = 1e6: current block %.2e, the rest %.2e from the reduced solve\n",
	             lastBlock, lastRest );
	BOOST_TEST( lastBlock < 1.0e-5 );
	BOOST_TEST( lastRest < 1.0e-10 );
}

/*
 * EVERY INCONSISTENT CALL IS REFUSED, BECAUSE A SILENT ONE WOULD BE A STEP.
 */
BOOST_AUTO_TEST_CASE( anInconsistentCallIsRefused )
{
	std::vector<double> const m( 6, 1.0 );
	std::vector<double> const f( 2, 1.0 );
	auto solve = []( std::vector<double> const &a, int r, int c, std::vector<double> const &b,
	                 std::vector<double> const &g = {} )
	{
		return meq::solveRegularisedLeastSquares( a, r, c, b, g );
	};
	BOOST_CHECK_THROW( solve( m, 2, 4, f ), std::invalid_argument );
	BOOST_CHECK_THROW( solve( m, 2, 3, { 1.0 } ), std::invalid_argument );
	BOOST_CHECK_THROW( solve( m, 2, 3, f, { 1.0, 1.0 } ), std::invalid_argument );
	BOOST_CHECK_THROW( solve( m, 2, 3, f, { 1.0, -1.0, 0.0 } ), std::invalid_argument );
	BOOST_CHECK_THROW( solve( m, 2, 3, f, { 1.0, std::nan( "" ), 0.0 } ), std::invalid_argument );
	BOOST_CHECK_THROW( solve( {}, 0, 3, {} ), std::invalid_argument );
	std::vector<double> bad = m;
	bad[ 3 ] = std::numeric_limits<double>::infinity();
	BOOST_CHECK_THROW( solve( bad, 2, 3, f ), std::invalid_argument );
	BOOST_CHECK_NO_THROW( solve( m, 2, 3, f, { 0.0, 1.0, 0.0 } ) );
}
