// Unit tests for the Zernike disc basis in src/meq/Zernike.
//
// INVERSION-PLAN.md section 4.1 builds the flux-surface representation
// R( rho, theta ), z( rho, theta ) on this basis, with rho = sqrt( Psi_N ). Two
// properties are the reason it was chosen over a naive "polynomial in rho times
// Fourier in theta", and they are what this file exists to MEASURE rather than
// to assert from the literature:
//
//   * SPECTRAL DECAY. A function analytic on the disc has Zernike coefficients
//     that fall GEOMETRICALLY with degree. That is the property IN-3 depends on
//     -- it is what makes a fit of a handful of modes a usable representation
//     rather than an interpolation table. Measured below against a function
//     that is merely continuous at the centre, whose coefficients fall
//     ALGEBRAICALLY; a decay measurement with no such control cannot tell the
//     two apart.
//
//   * AXIS REGULARITY. Every admissible mode -- l >= |m| with l - |m| even -- is
//     a polynomial in the Cartesian coordinates ( x, y ), so the centre of the
//     disc, which for meq is the magnetic axis, is an ordinary point. The
//     contrast is measured directly: a parity-violating mode rho^2 cos( theta )
//     has a second Cartesian derivative that jumps across the origin by a fixed
//     amount that does not shrink as the step does, while every admissible mode's
//     jump falls at O( h ).
//
// The other measurements here are the ones that say the evaluation and the
// normalisation are right at all: closed forms at low degree checked exactly,
// and orthogonality on the disc to the quadrature's own accuracy. The radial
// polynomial is evaluated through Boost.Math's Jacobi polynomial under a change
// of variable, so what those tests are really checking is the IDENTITY and its
// two sign conventions, which is exactly where such a thing goes wrong and is
// the reason to keep them rather than trusting the library.
//
// NOTHING HERE INCLUDES MFEM, and neither does the unit under test. That is what
// lets continuous integration build and run this file; see INSTALL.md.

#define BOOST_TEST_MODULE ZernikeTests
#ifndef BOOST_TEST_DYN_LINK
#define BOOST_TEST_DYN_LINK
#endif

#include <boost/test/unit_test.hpp>

#include "meq/Zernike.hpp"

// Only so that the two RIVAL conventions of the Zernike-to-Jacobi identity can
// be evaluated alongside the right one below. meq's own header pulls in nothing
// but <cstddef> and <vector>; this is the test reaching past it deliberately.
#include <boost/math/special_functions/jacobi.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{

	double const pi = M_PI;

	// -----------------------------------------------------------------------
	// Small numerical helpers.
	// -----------------------------------------------------------------------

	// Absolute comparison with a floor of one, so that a quantity passing through
	// zero -- which a Zernike polynomial does repeatedly -- does not make a
	// relative tolerance meaningless. Same helper, same reasoning, as
	// ProfilesTests.cpp.
	void checkClose( double actual, double expected, double tolerance, std::string const & what )
	{
		double const scale = std::max( 1.0, std::fabs( expected ) );

		BOOST_CHECK_MESSAGE( std::fabs( actual - expected ) <= tolerance*scale,
			what << ": got " << actual << ", expected " << expected
			<< " (difference " << actual - expected << ", allowed " << tolerance*scale << ")" );
	}

	// Gauss-Legendre nodes and weights on [ -1, 1 ], by Newton on the Legendre
	// polynomial. Used because the orthogonality and projection integrands are
	// polynomials in rho of known degree, so an n-point rule is EXACT for them
	// and the measured departure from orthogonality is then the basis evaluation's
	// error alone rather than the quadrature's.
	void gaussLegendre( int n, std::vector<double> & nodes, std::vector<double> & weights )
	{
		nodes.assign( static_cast<std::size_t>( n ), 0.0 );
		weights.assign( static_cast<std::size_t>( n ), 0.0 );

		for ( int i = 0; i < n; ++i )
		{
			double x = std::cos( pi*( i + 0.75 )/( n + 0.5 ) );
			double value = 0.0;
			double derivative = 1.0;

			for ( int iteration = 0; iteration < 200; ++iteration )
			{
				double previous = 0.0;
				value = 1.0;

				for ( int k = 0; k < n; ++k )
				{
					double const older = previous;
					previous = value;
					value = ( ( 2.0*k + 1.0 )*x*previous - k*older )/( k + 1.0 );
				}

				derivative = n*( x*value - previous )/( x*x - 1.0 );

				double const step = -value/derivative;
				x += step;

				if ( std::fabs( step ) < 1.0e-15 )
					break;
			}

			nodes[ static_cast<std::size_t>( i ) ] = x;
			weights[ static_cast<std::size_t>( i ) ] = 2.0/( ( 1.0 - x*x )*derivative*derivative );
		}
	}

	// The same rule mapped onto [ 0, 1 ], which is where rho lives.
	void gaussLegendreUnitInterval( int n, std::vector<double> & nodes, std::vector<double> & weights )
	{
		gaussLegendre( n, nodes, weights );

		for ( std::size_t i = 0; i < nodes.size(); ++i )
		{
			nodes[ i ] = 0.5*( nodes[ i ] + 1.0 );
			weights[ i ] *= 0.5;
		}
	}

	// A least-squares straight line y = slope x + intercept, with the coefficient
	// of determination. R^2 is what makes "geometric, not algebraic" a
	// measurement: the same data is fitted against log( c ) = a + b l and against
	// log( c ) = a + b log( l ), and the model that fits is the answer.
	struct LineFit
	{
		double slope;
		double intercept;
		double rSquared;
	};

	LineFit fitLine( std::vector<double> const & x, std::vector<double> const & y )
	{
		std::size_t const n = x.size();
		double sumX = 0.0, sumY = 0.0, sumXX = 0.0, sumXY = 0.0;

		for ( std::size_t i = 0; i < n; ++i )
		{
			sumX += x[ i ];
			sumY += y[ i ];
			sumXX += x[ i ]*x[ i ];
			sumXY += x[ i ]*y[ i ];
		}

		double const count = static_cast<double>( n );
		double const slope = ( count*sumXY - sumX*sumY )/( count*sumXX - sumX*sumX );
		double const intercept = ( sumY - slope*sumX )/count;

		double totalSumSquares = 0.0, residualSumSquares = 0.0;
		double const mean = sumY/count;

		for ( std::size_t i = 0; i < n; ++i )
		{
			double const predicted = slope*x[ i ] + intercept;
			totalSumSquares += ( y[ i ] - mean )*( y[ i ] - mean );
			residualSumSquares += ( y[ i ] - predicted )*( y[ i ] - predicted );
		}

		return LineFit{ slope, intercept, 1.0 - residualSumSquares/totalSumSquares };
	}

	// -----------------------------------------------------------------------
	// The explicit factorial sum, kept as a CONTROL and never as an
	// implementation.
	//
	// This is the textbook definition of R_l^m. It is here to do two jobs: to
	// check the Jacobi route against an independent formula at low degree, where
	// both are accurate; and to MEASURE the cancellation that makes it unusable
	// at moderate degree, which is why Zernike.cpp goes through Boost.Math's
	// Jacobi polynomial instead. If this ever agrees with that route to
	// round-off at l = 40, the comparison below is empty and should be deleted
	// rather than left as decoration.
	// -----------------------------------------------------------------------

	double factorial( int n )
	{
		double result = 1.0;

		for ( int i = 2; i <= n; ++i )
			result *= i;

		return result;
	}

	double naiveRadial( int l, int m, double rho )
	{
		int const absM = std::abs( m );
		int const half = ( l - absM )/2;
		double sum = 0.0;

		for ( int s = 0; s <= half; ++s )
		{
			double const term = factorial( l - s )
				/( factorial( s )*factorial( ( l + absM )/2 - s )*factorial( ( l - absM )/2 - s ) )
				*std::pow( rho, l - 2*s );

			sum += ( s % 2 == 0 ) ? term : -term;
		}

		return sum;
	}

	// The largest term of that sum, which is what sets how many digits it loses.
	double naiveLargestTerm( int l, int m, double rho )
	{
		int const absM = std::abs( m );
		int const half = ( l - absM )/2;
		double largest = 0.0;

		for ( int s = 0; s <= half; ++s )
		{
			double const term = factorial( l - s )
				/( factorial( s )*factorial( ( l + absM )/2 - s )*factorial( ( l - absM )/2 - s ) )
				*std::pow( rho, l - 2*s );

			largest = std::max( largest, std::fabs( term ) );
		}

		return largest;
	}

	// -----------------------------------------------------------------------
	// Projection onto the basis, by quadrature. Deliberately NOT part of the
	// library: see the note in Zernike.hpp about IN-3 fitting a point cloud
	// rather than a function. It is the sharpest instrument available for
	// measuring coefficient decay, which is a different job.
	// -----------------------------------------------------------------------

	using DiscFunction = std::function<double( double, double )>;   // ( x, y ) -> f

	meq::ZernikeExpansion project( DiscFunction const & f, int maxDegree, int radialPoints, int angularPoints )
	{
		std::vector<double> radialNodes, radialWeights;
		gaussLegendreUnitInterval( radialPoints, radialNodes, radialWeights );

		std::vector<meq::ZernikeMode> const modes = meq::zernikeModes( maxDegree );
		std::vector<double> coefficients( modes.size(), 0.0 );

		// f at every quadrature point, once.
		std::vector<double> sampled( static_cast<std::size_t>( radialPoints*angularPoints ), 0.0 );

		for ( int i = 0; i < radialPoints; ++i )
			for ( int j = 0; j < angularPoints; ++j )
			{
				double const rho = radialNodes[ static_cast<std::size_t>( i ) ];
				double const theta = 2.0*pi*j/angularPoints;
				sampled[ static_cast<std::size_t>( i*angularPoints + j ) ] = f( rho*std::cos( theta ), rho*std::sin( theta ) );
			}

		for ( std::size_t mode = 0; mode < modes.size(); ++mode )
		{
			int const l = modes[ mode ].l;
			int const m = modes[ mode ].m;
			double inner = 0.0;

			for ( int i = 0; i < radialPoints; ++i )
			{
				double const rho = radialNodes[ static_cast<std::size_t>( i ) ];
				double const radial = meq::zernikeRadial( l, m, rho );
				double angularSum = 0.0;

				for ( int j = 0; j < angularPoints; ++j )
				{
					double const theta = 2.0*pi*j/angularPoints;
					double const angular = ( m >= 0 ) ? std::cos( m*theta ) : std::sin( -m*theta );
					angularSum += sampled[ static_cast<std::size_t>( i*angularPoints + j ) ]*angular;
				}

				// The equispaced trapezoidal rule in theta: weight 2 pi / N, and
				// no endpoint halving because the integrand is periodic.
				inner += radialWeights[ static_cast<std::size_t>( i ) ]*rho*radial*angularSum*( 2.0*pi/angularPoints );
			}

			coefficients[ mode ] = inner/meq::zernikeNormSquared( l, m );
		}

		return meq::ZernikeExpansion( maxDegree, std::move( coefficients ) );
	}

	// The largest coefficient magnitude at each degree -- the decay envelope,
	// which is what a statement about "mode number" is actually about.
	std::vector<double> degreeEnvelope( meq::ZernikeExpansion const & expansion )
	{
		std::vector<double> envelope( static_cast<std::size_t>( expansion.maxDegree() + 1 ), 0.0 );
		std::vector<meq::ZernikeMode> const & modes = expansion.modes();

		for ( std::size_t i = 0; i < modes.size(); ++i )
		{
			std::size_t const degree = static_cast<std::size_t>( modes[ i ].l );
			envelope[ degree ] = std::max( envelope[ degree ], std::fabs( expansion.coefficients()[ i ] ) );
		}

		return envelope;
	}

	// The largest error of the expansion against f over a sample of the disc.
	double maximumError( meq::ZernikeExpansion const & expansion, DiscFunction const & f )
	{
		double worst = 0.0;

		for ( int i = 0; i <= 40; ++i )
			for ( int j = 0; j < 37; ++j )
			{
				double const rho = i/40.0;
				double const theta = 2.0*pi*j/37.0;
				double const exact = f( rho*std::cos( theta ), rho*std::sin( theta ) );
				worst = std::max( worst, std::fabs( expansion( rho, theta ) - exact ) );
			}

		return worst;
	}

	// -----------------------------------------------------------------------
	// The smoothness probe used by the axis tests.
	//
	// Restrict a function on the disc to a straight line through the origin at
	// angle alpha and take a ONE-SIDED second difference from each side. For a
	// polynomial in ( x, y ) the restriction is a polynomial in the line
	// parameter, so the two agree to O( h ). For a mode violating the parity
	// constraint the restriction is a function such as t|t|, whose two one-sided
	// second derivatives differ by a fixed amount at every h.
	//
	// A CENTRED second difference is useless here and that is worth recording:
	// t|t| is ODD, so ( g(h) - 2g(0) + g(-h) )/h^2 is identically zero and the
	// defect cancels itself out of the measurement.
	// -----------------------------------------------------------------------

	using PolarFunction = std::function<double( double, double )>;   // ( rho, theta ) -> value

	double alongLine( PolarFunction const & g, double alpha, double t )
	{
		return ( t >= 0.0 ) ? g( t, alpha ) : g( -t, alpha + pi );
	}

	double oneSidedSecondDifferenceJump( PolarFunction const & g, double alpha, double h )
	{
		double const centre = alongLine( g, alpha, 0.0 );
		double const right = ( alongLine( g, alpha, 2.0*h ) - 2.0*alongLine( g, alpha, h ) + centre )/( h*h );
		double const left = ( alongLine( g, alpha, -2.0*h ) - 2.0*alongLine( g, alpha, -h ) + centre )/( h*h );

		return std::fabs( right - left );
	}

}

// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE( mode_enumeration_tests )

BOOST_AUTO_TEST_CASE( admits_exactly_the_modes_with_l_minus_m_even )
{
	// The constraint, stated as a predicate and checked against its own
	// definition over a block of index pairs -- including the ones a naive
	// tensor-product basis would have admitted.
	for ( int l = 0; l <= 8; ++l )
		for ( int m = -10; m <= 10; ++m )
		{
			bool const expected = std::abs( m ) <= l && ( l - std::abs( m ) ) % 2 == 0;
			BOOST_CHECK_MESSAGE( meq::isValidZernikeMode( l, m ) == expected,
				"mode ( " << l << ", " << m << " ) admissibility" );
		}

	BOOST_CHECK( !meq::isValidZernikeMode( -1, 0 ) );

	// The three that matter: rho^2 cos( theta ) and rho cos( 2 theta ) are the
	// modes the parity constraint exists to exclude, and rho^3 cos( theta ) is
	// admissible even though it looks like it should not be.
	BOOST_CHECK( !meq::isValidZernikeMode( 2, 1 ) );
	BOOST_CHECK( !meq::isValidZernikeMode( 1, 2 ) );
	BOOST_CHECK( meq::isValidZernikeMode( 3, 1 ) );
}

BOOST_AUTO_TEST_CASE( the_mode_count_is_triangular_and_the_order_is_prefix_stable )
{
	for ( int degree = 0; degree <= 12; ++degree )
	{
		std::vector<meq::ZernikeMode> const modes = meq::zernikeModes( degree );

		BOOST_CHECK_EQUAL( modes.size(), meq::zernikeModeCount( degree ) );
		BOOST_CHECK_EQUAL( modes.size(), static_cast<std::size_t>( ( degree + 1 )*( degree + 2 )/2 ) );

		// Every enumerated mode is admissible, and its index is where it sits.
		for ( std::size_t i = 0; i < modes.size(); ++i )
		{
			BOOST_CHECK( meq::isValidZernikeMode( modes[ i ].l, modes[ i ].m ) );
			BOOST_CHECK_EQUAL( meq::zernikeModeIndex( modes[ i ].l, modes[ i ].m ), i );
		}
	}

	// Prefix stability: truncating the mode list IS truncating the coefficient
	// vector, which is what a convergence study in mode number relies on.
	std::vector<meq::ZernikeMode> const large = meq::zernikeModes( 12 );

	for ( int degree = 0; degree <= 12; ++degree )
	{
		std::vector<meq::ZernikeMode> const small = meq::zernikeModes( degree );

		for ( std::size_t i = 0; i < small.size(); ++i )
			BOOST_CHECK( small[ i ] == large[ i ] );
	}
}

BOOST_AUTO_TEST_CASE( refuses_an_inadmissible_mode_loudly )
{
	// The whole point of the basis is that these modes are not in it, so an
	// entry point that quietly returned something for them would be worse than
	// useless. std::invalid_argument, as the Profiles constructors throw.
	BOOST_CHECK_THROW( meq::zernikeRadial( 2, 1, 0.5 ), std::invalid_argument );
	BOOST_CHECK_THROW( meq::zernikeRadialPrime( 2, 1, 0.5 ), std::invalid_argument );
	BOOST_CHECK_THROW( meq::zernike( 2, 1, 0.5, 0.3 ), std::invalid_argument );
	BOOST_CHECK_THROW( meq::zernikeRadialDerivative( 2, 1, 0.5, 0.3 ), std::invalid_argument );
	BOOST_CHECK_THROW( meq::zernikeAngularDerivative( 2, 1, 0.5, 0.3 ), std::invalid_argument );
	BOOST_CHECK_THROW( meq::zernikeNormSquared( 2, 1 ), std::invalid_argument );
	BOOST_CHECK_THROW( meq::zernikeModeIndex( 2, 1 ), std::invalid_argument );

	BOOST_CHECK_THROW( meq::zernike( 1, 3, 0.5, 0.3 ), std::invalid_argument );
	BOOST_CHECK_THROW( meq::zernike( -1, -1, 0.5, 0.3 ), std::invalid_argument );

	BOOST_CHECK_THROW( meq::zernikeModes( -1 ), std::invalid_argument );
	BOOST_CHECK_THROW( meq::zernikeModes( meq::maxZernikeDegree + 1 ), std::invalid_argument );

	BOOST_CHECK_THROW( meq::ZernikeExpansion( 3, std::vector<double>( 9, 0.0 ) ), std::invalid_argument );

	// And the message names the parity condition rather than saying "invalid",
	// because a caller who has just looped over every |m| <= l has made a
	// modelling error and this is the only place they are told which.
	try
	{
		meq::zernikeRadial( 2, 1, 0.5 );
		BOOST_FAIL( "expected an exception" );
	}
	catch ( std::invalid_argument const & error )
	{
		BOOST_CHECK( std::string( error.what() ).find( "even" ) != std::string::npos );
	}
}

BOOST_AUTO_TEST_SUITE_END()

// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE( closed_form_tests )

BOOST_AUTO_TEST_CASE( the_jacobi_identity_uses_the_right_argument_and_the_right_sign )
{
	// THE CHEAPEST DISCRIMINATOR AMONG THE THREE CANDIDATE CONVENTIONS, written
	// down because a future reader will wonder which way round it goes and
	// because two of the three are silent for the first several modes.
	//
	// src/meq/Zernike.cpp evaluates
	//
	//     R_l^m( rho ) = ( -1 )^k rho^|m| P_k^( |m|, 0 )( 1 - 2 rho^2 ),
	//                    k = ( l - |m| )/2.
	//
	// The two ways to get that wrong are to take the argument as 2 rho^2 - 1 and
	// to drop the leading ( -1 )^k. BOTH REPRODUCE EVERY PURE POWER R_l^l
	// EXACTLY, because there k = 0, the Jacobi factor is the constant one, and
	// neither the argument nor the sign can be seen -- so a test built only on
	// R_0^0, R_1^1 and R_2^2 passes for all three conventions. l = 3, m = 1 is
	// the smallest case that separates them, and R_l^m( 1 ) = 1 is the other.
	//
	// The rivals are EVALUATED here rather than described, so "the other
	// convention passes neither check" is a measurement.
	auto const candidate = []( int l, int m, double rho, bool flipArgument, bool dropSign )
	{
		unsigned const order = static_cast<unsigned>( ( l - m )/2 );
		double const argument = flipArgument ? ( 2.0*rho*rho - 1.0 ) : ( 1.0 - 2.0*rho*rho );
		double const parity = ( dropSign || order % 2 == 0 ) ? 1.0 : -1.0;

		return parity*std::pow( rho, m )*boost::math::jacobi( order, static_cast<double>( m ), 0.0, argument );
	};

	// First: all three agree on the pure powers, which is what makes the other
	// two dangerous rather than obviously wrong.
	for ( int l : { 0, 1, 2, 3, 4, 5 } )
		for ( double rho : { 0.3, 0.7, 1.0 } )
		{
			checkClose( candidate( l, l, rho, true, false ), meq::zernikeRadial( l, l, rho ), 1.0e-15, "flipped argument on a pure power R_l^l" );
			checkClose( candidate( l, l, rho, false, true ), meq::zernikeRadial( l, l, rho ), 1.0e-15, "dropped sign on a pure power R_l^l" );
		}

	// Then: R_3^1 = 3 rho^3 - 2 rho, which only the right convention gives. The
	// flipped argument gives rho - 3 rho^3 and the dropped sign gives
	// 2 rho - 3 rho^3, both of which are checked against here so that a silent
	// swap of the convention would fail this test in a way that names the cause.
	std::cout << "\n  R_3^1 under the three candidate conventions, at rho = 0.7\n"
	          << "    correct  ( -1 )^k, 1 - 2 rho^2 : " << std::fixed << std::setprecision( 12 ) << meq::zernikeRadial( 3, 1, 0.7 ) << "\n"
	          << "    rival    ( -1 )^k, 2 rho^2 - 1 : " << candidate( 3, 1, 0.7, true, false ) << "\n"
	          << "    rival    no sign, 1 - 2 rho^2  : " << candidate( 3, 1, 0.7, false, true ) << "\n"
	          << "    the closed form 3 rho^3 - 2 rho: " << 3.0*0.343 - 1.4 << std::defaultfloat << std::endl;

	for ( double rho : { 0.0, 0.2, 0.45, 0.7, 0.9, 1.0 } )
	{
		double const r3 = rho*rho*rho;

		checkClose( meq::zernikeRadial( 3, 1, rho ), 3.0*r3 - 2.0*rho, 1.0e-15, "R_3^1 against its closed form" );
		checkClose( candidate( 3, 1, rho, true, false ), rho - 3.0*r3, 1.0e-15, "the flipped-argument rival's own closed form" );
		checkClose( candidate( 3, 1, rho, false, true ), 2.0*rho - 3.0*r3, 1.0e-15, "the dropped-sign rival's own closed form" );
	}

	// And the normalisation, over every admissible mode to degree 30 -- the
	// second discriminator, and the one that is cheap enough to run over the
	// whole basis rather than over one mode.
	int flippedFailures = 0;
	int unsignedFailures = 0;
	int modes = 0;

	for ( int l = 0; l <= 30; ++l )
		for ( int m = 0; m <= l; m += 2 )
		{
			if ( ( l - m ) % 2 != 0 )
				continue;

			++modes;

			checkClose( meq::zernikeRadial( l, m, 1.0 ), 1.0, 2.0e-14,
				"R_" + std::to_string( l ) + "^" + std::to_string( m ) + "( 1 ) = 1" );

			if ( std::fabs( candidate( l, m, 1.0, true, false ) - 1.0 ) > 1.0e-10 )
				++flippedFailures;

			if ( std::fabs( candidate( l, m, 1.0, false, true ) - 1.0 ) > 1.0e-10 )
				++unsignedFailures;
		}

	std::cout << "  R_l^m( 1 ) = 1 over " << modes << " modes to degree 30: the right convention holds everywhere; "
	          << "the flipped argument fails on " << flippedFailures << " of them and the dropped sign on "
	          << unsignedFailures << std::endl;

	// If either rival ever passed, this whole test would be measuring nothing.
	BOOST_CHECK_MESSAGE( flippedFailures > 0, "the 2 rho^2 - 1 convention was expected to break the normalisation and did not" );
	BOOST_CHECK_MESSAGE( unsignedFailures > 0, "dropping the ( -1 )^k was expected to break the normalisation and did not" );
}

BOOST_AUTO_TEST_CASE( reproduces_the_low_order_radial_polynomials_exactly )
{
	// Every one of these is a published closed form. Checked at a spread of rho
	// rather than at a point, so a polynomial of the right degree with a wrong
	// coefficient cannot pass.
	for ( double rho : { 0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0 } )
	{
		double const r2 = rho*rho;
		double const r3 = r2*rho;
		double const r4 = r3*rho;
		double const r5 = r4*rho;
		double const r6 = r5*rho;

		checkClose( meq::zernikeRadial( 0, 0, rho ), 1.0, 1.0e-15, "R_0^0" );
		checkClose( meq::zernikeRadial( 1, 1, rho ), rho, 1.0e-15, "R_1^1" );
		checkClose( meq::zernikeRadial( 2, 0, rho ), 2.0*r2 - 1.0, 1.0e-15, "R_2^0" );
		checkClose( meq::zernikeRadial( 2, 2, rho ), r2, 1.0e-15, "R_2^2" );
		checkClose( meq::zernikeRadial( 3, 1, rho ), 3.0*r3 - 2.0*rho, 1.0e-15, "R_3^1" );
		checkClose( meq::zernikeRadial( 3, 3, rho ), r3, 1.0e-15, "R_3^3" );
		checkClose( meq::zernikeRadial( 4, 0, rho ), 6.0*r4 - 6.0*r2 + 1.0, 1.0e-15, "R_4^0" );
		checkClose( meq::zernikeRadial( 4, 2, rho ), 4.0*r4 - 3.0*r2, 1.0e-15, "R_4^2" );
		checkClose( meq::zernikeRadial( 4, 4, rho ), r4, 1.0e-15, "R_4^4" );
		checkClose( meq::zernikeRadial( 5, 1, rho ), 10.0*r5 - 12.0*r3 + 3.0*rho, 1.0e-15, "R_5^1" );
		checkClose( meq::zernikeRadial( 5, 3, rho ), 5.0*r5 - 4.0*r3, 1.0e-15, "R_5^3" );
		checkClose( meq::zernikeRadial( 6, 0, rho ), 20.0*r6 - 30.0*r4 + 12.0*r2 - 1.0, 1.0e-15, "R_6^0" );
		checkClose( meq::zernikeRadial( 6, 2, rho ), 15.0*r6 - 20.0*r4 + 6.0*r2, 1.0e-15, "R_6^2" );

		// The sign of m selects the angular function and nothing else, so the
		// radial polynomial must not notice it.
		checkClose( meq::zernikeRadial( 3, -1, rho ), meq::zernikeRadial( 3, 1, rho ), 1.0e-15, "R_3^-1 against R_3^1" );
	}
}

BOOST_AUTO_TEST_CASE( the_low_order_modes_are_the_expected_cartesian_polynomials )
{
	// THIS IS THE PARITY CONSTRAINT MADE CONCRETE. Z_1^1 is x, Z_2^-2 is 2xy,
	// Z_2^2 is x^2 - y^2: polynomials in the Cartesian coordinates, with no
	// square root and no branch anywhere, which is exactly what an admissible
	// mode always is and an inadmissible one never is.
	for ( double rho : { 0.0, 0.3, 0.7, 1.0 } )
		for ( int j = 0; j < 11; ++j )
		{
			double const theta = 2.0*pi*j/11.0;
			double const x = rho*std::cos( theta );
			double const y = rho*std::sin( theta );

			checkClose( meq::zernike( 0, 0, rho, theta ), 1.0, 1.0e-15, "Z_0^0 = 1" );
			checkClose( meq::zernike( 1, 1, rho, theta ), x, 1.0e-15, "Z_1^1 = x" );
			checkClose( meq::zernike( 1, -1, rho, theta ), y, 1.0e-15, "Z_1^-1 = y" );
			checkClose( meq::zernike( 2, 2, rho, theta ), x*x - y*y, 1.0e-15, "Z_2^2 = x^2 - y^2" );
			checkClose( meq::zernike( 2, -2, rho, theta ), 2.0*x*y, 1.0e-15, "Z_2^-2 = 2xy" );
			checkClose( meq::zernike( 2, 0, rho, theta ), 2.0*( x*x + y*y ) - 1.0, 1.0e-15, "Z_2^0" );
			checkClose( meq::zernike( 3, 3, rho, theta ), x*x*x - 3.0*x*y*y, 1.0e-15, "Z_3^3" );
			checkClose( meq::zernike( 3, 1, rho, theta ), 3.0*x*( x*x + y*y ) - 2.0*x, 1.0e-15, "Z_3^1" );
		}
}

BOOST_AUTO_TEST_CASE( is_normalised_to_one_at_the_rim_and_correct_on_the_axis )
{
	// R_l^m( 1 ) = 1 is the normalisation the whole basis is stated in, and it
	// is a sharp test at high degree: the naive sum below loses eight digits on
	// it at l = 30. R_l^m( 0 ) is zero for m /= 0 -- which is why an expansion
	// is single valued at the centre whatever the coefficients are -- and
	// ( -1 )^( l/2 ) for m = 0.
	for ( int l = 0; l <= 40; ++l )
		for ( int m = -l; m <= l; m += 2 )
		{
			checkClose( meq::zernikeRadial( l, m, 1.0 ), 1.0, 2.0e-14, "R_" + std::to_string( l ) + "^" + std::to_string( m ) + "( 1 )" );

			double const expectedAtZero = ( m == 0 ) ? ( ( ( l/2 ) % 2 == 0 ) ? 1.0 : -1.0 ) : 0.0;
			checkClose( meq::zernikeRadial( l, m, 0.0 ), expectedAtZero, 1.0e-15,
				"R_" + std::to_string( l ) + "^" + std::to_string( m ) + "( 0 )" );
		}
}

BOOST_AUTO_TEST_CASE( the_jacobi_route_agrees_with_the_factorial_sum_and_then_beats_it )
{
	// Two jobs. At low degree both formulae are accurate, so agreement is an
	// independent check that the Jacobi identity produces the right polynomial.
	for ( int l = 0; l <= 10; ++l )
		for ( int m = -l; m <= l; m += 2 )
			for ( double rho : { 0.0, 0.2, 0.45, 0.8, 1.0 } )
				checkClose( meq::zernikeRadial( l, m, rho ), naiveRadial( l, m, rho ), 1.0e-12,
					"R_" + std::to_string( l ) + "^" + std::to_string( m ) + " against the factorial sum at rho = " + std::to_string( rho ) );

	// At moderate degree they part company, and this is the measurement that says
	// the sum is not an acceptable implementation. R_l^m( 1 ) = 1
	// exactly, so the error of each route is readable directly.
	//
	// A CAVEAT ON THIS COLUMN, because it makes the table UNDERSTATE the problem
	// and was surprising when first run: at rho = 1 every term of the factorial
	// sum is an INTEGER, and a sum of integers is exact in binary floating point
	// until one of them passes 2^53. So the naive sum is not merely accurate at
	// rho = 1 for l <= 30, it is exact, and the cancellation only shows once the
	// terms clear about 9e15 -- which happens between l = 40 and l = 50 here.
	// Away from rho = 1 there is no such reprieve, which is what the second
	// table below measures.
	std::cout << "\n  the factorial sum against the Jacobi route, at R_l^0( 1 ) = 1 exactly\n"
	          << "     l   largest term   naive error     Jacobi error\n";

	bool naiveIsWorse = false;

	for ( int l : { 10, 20, 30, 40, 50 } )
	{
		double const naiveError = std::fabs( naiveRadial( l, 0, 1.0 ) - 1.0 );
		double const jacobiError = std::fabs( meq::zernikeRadial( l, 0, 1.0 ) - 1.0 );

		std::cout << "  " << std::setw( 4 ) << l
		          << std::scientific << std::setprecision( 3 )
		          << "   " << std::setw( 12 ) << naiveLargestTerm( l, 0, 1.0 )
		          << "   " << std::setw( 11 ) << naiveError
		          << "   " << std::setw( 16 ) << jacobiError << "\n";

		if ( naiveError > 1.0e-8 )
			naiveIsWorse = true;
	}

	std::cout << std::defaultfloat << std::endl;

	// Off the rim, where the terms are not integers, the two routes are compared
	// against each other. THE LOGIC OF ATTRIBUTING THE DIFFERENCE TO THE NAIVE
	// SUM IS NOT CIRCULAR, and it is worth spelling out: the Jacobi route has
	// already been shown orthogonal to 1e-14 out to degree 30 by
	// the_radial_polynomials_are_orthogonal_at_high_degree, which the sum played
	// no part in. A polynomial that is orthogonal to all its neighbours to
	// round-off is the right polynomial to round-off, so a disagreement at 1e-7
	// is the sum's.
	std::cout << "  the same comparison at rho = 0.83, where the terms are not integers\n"
	          << "     l   largest term      |naive - Jacobi|\n";

	double worstDisagreement = 0.0;

	for ( int l : { 10, 20, 30, 40 } )
	{
		double const difference = std::fabs( naiveRadial( l, 0, 0.83 ) - meq::zernikeRadial( l, 0, 0.83 ) );
		worstDisagreement = std::max( worstDisagreement, difference );

		std::cout << "  " << std::setw( 4 ) << l
		          << std::scientific << std::setprecision( 3 )
		          << "   " << std::setw( 12 ) << naiveLargestTerm( l, 0, 0.83 )
		          << "   " << std::setw( 20 ) << difference << "\n";
	}

	std::cout << std::defaultfloat << std::endl;

	BOOST_CHECK_MESSAGE( worstDisagreement > 1.0e-10,
		"the factorial sum was expected to have lost several digits by degree 40; it disagreed with the Jacobi route by only "
		<< worstDisagreement << ", so this comparison measures nothing" );

	// If this ever fails, the factorial sum has stopped being a control and the
	// comparison above is empty -- delete it rather than relaxing it.
	BOOST_CHECK_MESSAGE( naiveIsWorse,
		"the factorial sum was expected to lose several digits to cancellation by l = 30; if it no longer does, this comparison measures nothing" );

	// And the Jacobi route itself is still at round-off there, which is the claim
	// that matters for a surface fit wanting twenty or thirty modes.
	for ( int l : { 20, 30, 40 } )
		for ( int m : { 0, 2, l/2 - ( l/2 ) % 2, l } )
			if ( meq::isValidZernikeMode( l, m ) )
				checkClose( meq::zernikeRadial( l, m, 1.0 ), 1.0, 2.0e-14,
					"R_" + std::to_string( l ) + "^" + std::to_string( m ) + "( 1 ) at high degree" );
}

BOOST_AUTO_TEST_SUITE_END()

// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE( orthogonality_tests )

BOOST_AUTO_TEST_CASE( the_radial_polynomials_are_orthogonal_at_high_degree )
{
	// The strongest single statement about the evaluation, and the reason to keep
	// this test rather than trusting Boost: it is where a coordinate-change
	// convention error in the Zernike-to-Jacobi identity would hide. The integrand
	// R_l^m R_l'^m rho is a polynomial of degree at most 2*30 + 1 = 61, so a
	// 50-point Gauss-Legendre rule -- exact to degree 99 -- contributes nothing
	// of its own and what is measured is the identity and its normalisation.
	int const maxDegree = 30;
	std::vector<double> nodes, weights;
	gaussLegendreUnitInterval( 50, nodes, weights );

	double worstOffDiagonal = 0.0;
	double worstDiagonal = 0.0;

	for ( int m = 0; m <= 6; ++m )
		for ( int l = m; l <= maxDegree; l += 2 )
			for ( int lPrime = m; lPrime <= maxDegree; lPrime += 2 )
			{
				double inner = 0.0;

				for ( std::size_t i = 0; i < nodes.size(); ++i )
					inner += weights[ i ]*nodes[ i ]*meq::zernikeRadial( l, m, nodes[ i ] )*meq::zernikeRadial( lPrime, m, nodes[ i ] );

				if ( l == lPrime )
					worstDiagonal = std::max( worstDiagonal, std::fabs( inner - meq::zernikeRadialNormSquared( l, m ) )*2.0*( l + 1.0 ) );
				else
					worstOffDiagonal = std::max( worstOffDiagonal, std::fabs( inner )*2.0*( l + 1.0 ) );
			}

	std::cout << "\n  radial orthogonality to degree " << maxDegree
	          << ": worst off-diagonal " << std::scientific << std::setprecision( 3 ) << worstOffDiagonal
	          << ", worst diagonal error " << worstDiagonal << std::defaultfloat << std::endl;

	BOOST_CHECK_MESSAGE( worstOffDiagonal < 1.0e-12, "worst relative off-diagonal entry " << worstOffDiagonal );
	BOOST_CHECK_MESSAGE( worstDiagonal < 1.0e-12, "worst relative diagonal error " << worstDiagonal );
}

BOOST_AUTO_TEST_CASE( the_modes_are_orthogonal_on_the_disc_with_the_stated_norms )
{
	// The full two-dimensional statement, including the angular factor -- and so
	// including the factor of two on m = 0 that zernikeNormSquared() carries. Get
	// that wrong and every m = 0 coefficient of a projection is out by two, which
	// looks like a fit that is merely poor.
	int const maxDegree = 8;
	int const angularPoints = 64;

	std::vector<double> nodes, weights;
	gaussLegendreUnitInterval( 24, nodes, weights );

	std::vector<meq::ZernikeMode> const modes = meq::zernikeModes( maxDegree );
	double worstOffDiagonal = 0.0;
	double worstDiagonal = 0.0;

	for ( std::size_t a = 0; a < modes.size(); ++a )
		for ( std::size_t b = a; b < modes.size(); ++b )
		{
			double inner = 0.0;

			for ( std::size_t i = 0; i < nodes.size(); ++i )
				for ( int j = 0; j < angularPoints; ++j )
				{
					double const theta = 2.0*pi*j/angularPoints;
					inner += weights[ i ]*nodes[ i ]*( 2.0*pi/angularPoints )
						*meq::zernike( modes[ a ].l, modes[ a ].m, nodes[ i ], theta )
						*meq::zernike( modes[ b ].l, modes[ b ].m, nodes[ i ], theta );
				}

			double const scale = std::sqrt( meq::zernikeNormSquared( modes[ a ].l, modes[ a ].m )*meq::zernikeNormSquared( modes[ b ].l, modes[ b ].m ) );

			if ( a == b )
				worstDiagonal = std::max( worstDiagonal, std::fabs( inner - meq::zernikeNormSquared( modes[ a ].l, modes[ a ].m ) )/scale );
			else
				worstOffDiagonal = std::max( worstOffDiagonal, std::fabs( inner )/scale );
		}

	std::cout << "  disc orthogonality to degree " << maxDegree
	          << ": worst off-diagonal " << std::scientific << std::setprecision( 3 ) << worstOffDiagonal
	          << ", worst diagonal error " << worstDiagonal << std::defaultfloat << std::endl;

	BOOST_CHECK_MESSAGE( worstOffDiagonal < 1.0e-13, "worst relative off-diagonal entry " << worstOffDiagonal );
	BOOST_CHECK_MESSAGE( worstDiagonal < 1.0e-13, "worst relative diagonal error " << worstDiagonal );
}

BOOST_AUTO_TEST_SUITE_END()

// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE( derivative_tests )

BOOST_AUTO_TEST_CASE( the_radial_derivative_matches_a_central_difference )
{
	// The cases deliberately included: m = 0, where the polynomial is even in rho
	// and the m rho^( m - 1 ) term of the chain rule must drop out exactly rather
	// than being computed as zero times an infinity; and l = |m|, where the
	// polynomial is the pure power rho^|m| and the Jacobi factor is the constant
	// one, so any error in the seeding is naked.
	struct Case { int l; int m; };
	std::vector<Case> const cases = {
		{ 0, 0 }, { 1, 1 }, { 1, -1 }, { 2, 0 }, { 2, 2 }, { 3, 1 }, { 3, 3 },
		{ 4, 0 }, { 4, 2 }, { 4, 4 }, { 6, 0 }, { 7, 1 }, { 8, 4 }, { 10, 0 },
		{ 12, 2 }, { 15, 5 }, { 20, 0 }, { 20, 20 } };

	double const step = 1.0e-5;
	double worst = 0.0;
	double worstRichardson = 0.0;

	for ( Case const & c : cases )
		for ( double rho : { 0.05, 0.2, 0.4, 0.6, 0.8, 0.95 } )
		{
			double const difference = ( meq::zernikeRadial( c.l, c.m, rho + step ) - meq::zernikeRadial( c.l, c.m, rho - step ) )/( 2.0*step );
			double const exact = meq::zernikeRadialPrime( c.l, c.m, rho );
			double const scale = std::max( 1.0, std::fabs( exact ) );

			worst = std::max( worst, std::fabs( exact - difference )/scale );

			checkClose( difference, exact, 1.0e-6,
				"dR_" + std::to_string( c.l ) + "^" + std::to_string( c.m ) + "/drho at rho = " + std::to_string( rho ) );

			// AND THE SAME COMPARISON WITH THE DIFFERENCE'S OWN ERROR TAKEN OUT,
			// which is the only way to say anything about how good the ANALYTIC
			// derivative is.
			//
			// The number in the plain column above is not a statement about
			// zernikeRadialPrime() at all: a central difference carries its own
			// O( h^2 ) truncation, so the comparison is floored by the
			// instrument whatever the derivative does. Improving the derivative
			// cannot move it. Richardson extrapolation -- ( 4 D( h/2 ) - D( h ) )/3
			// -- cancels that h^2 term and leaves O( h^4 ) plus round-off, and
			// what is measured then really is the derivative.
			double const coarse = ( meq::zernikeRadial( c.l, c.m, rho + 1.0e-4 ) - meq::zernikeRadial( c.l, c.m, rho - 1.0e-4 ) )/2.0e-4;
			double const fine = ( meq::zernikeRadial( c.l, c.m, rho + 5.0e-5 ) - meq::zernikeRadial( c.l, c.m, rho - 5.0e-5 ) )/1.0e-4;
			double const richardson = ( 4.0*fine - coarse )/3.0;

			worstRichardson = std::max( worstRichardson, std::fabs( exact - richardson )/scale );
		}

	std::cout << "\n  radial derivative against a difference, over " << cases.size() << " modes to degree 20 and six radii\n"
	          << "    plain central difference, h = 1e-5      : " << std::scientific << std::setprecision( 3 ) << worst
	          << "   <- the difference's own O( h^2 ), not the derivative's error\n"
	          << "    Richardson extrapolated, h = 1e-4, 5e-5 : " << worstRichardson
	          << "   <- what the analytic derivative actually reaches\n" << std::defaultfloat << std::endl;

	// With the instrument's own second-order error removed, the analytic
	// derivative and the difference agree at what is left: the O( h^4 ) tail and
	// the round-off floor of a difference taken at h = 5e-4. That is four to five
	// orders tighter than the plain column, and it is the assertion with content.
	BOOST_CHECK_MESSAGE( worstRichardson < 1.0e-9,
		"worst Richardson-extrapolated mismatch was " << worstRichardson << ", which is not the round-off floor an analytic derivative should reach" );
}

BOOST_AUTO_TEST_CASE( the_radial_derivative_error_is_at_the_second_order_floor )
{
	// Not merely "close": the disagreement must fall as h^2, which is what says
	// the exact derivative is the limit of the difference rather than something
	// near it. Below about h = 1e-5 the round-off in the difference takes over,
	// so the sweep stops there.
	std::cout << "\n  central-difference convergence of dR_9^3/drho at rho = 0.37\n"
	          << "         h          error       ratio\n";

	int const l = 9;
	int const m = 3;
	double const rho = 0.37;
	double const exact = meq::zernikeRadialPrime( l, m, rho );

	double previous = 0.0;
	double worstRatio = 0.0;

	for ( int i = 0; i < 4; ++i )
	{
		double const h = 1.0e-2/std::pow( 2.0, i );
		double const difference = ( meq::zernikeRadial( l, m, rho + h ) - meq::zernikeRadial( l, m, rho - h ) )/( 2.0*h );
		double const error = std::fabs( difference - exact );
		double const ratio = ( i == 0 ) ? 0.0 : previous/error;

		std::cout << "   " << std::scientific << std::setprecision( 3 ) << std::setw( 10 ) << h
		          << "   " << std::setw( 12 ) << error;

		if ( i > 0 )
		{
			std::cout << "   " << std::fixed << std::setprecision( 3 ) << std::setw( 8 ) << ratio;
			worstRatio = ( worstRatio == 0.0 ) ? ratio : std::min( worstRatio, ratio );
		}

		std::cout << std::defaultfloat << "\n";
		previous = error;
	}

	std::cout << std::endl;

	// Four is the second-order value; anything from 3.5 up says the difference is
	// converging at h^2 and not at h.
	BOOST_CHECK_MESSAGE( worstRatio > 3.5, "worst error ratio between successive halvings was " << worstRatio << ", expected about 4" );
}

BOOST_AUTO_TEST_CASE( the_angular_derivative_is_exact_and_vanishes_identically_for_m_zero )
{
	double const step = 1.0e-6;

	for ( int l = 0; l <= 8; ++l )
		for ( int m = -l; m <= l; m += 2 )
			for ( double rho : { 0.15, 0.5, 0.9 } )
				for ( int j = 0; j < 7; ++j )
				{
					double const theta = 2.0*pi*j/7.0;
					double const difference = ( meq::zernike( l, m, rho, theta + step ) - meq::zernike( l, m, rho, theta - step ) )/( 2.0*step );

					checkClose( meq::zernikeAngularDerivative( l, m, rho, theta ), difference, 1.0e-8,
						"dZ_" + std::to_string( l ) + "^" + std::to_string( m ) + "/dtheta" );

					// Exactly zero, not merely small: an m = 0 mode has no
					// angular dependence at all and a consumer differentiating a
					// flux-surface shape is entitled to rely on it.
					if ( m == 0 )
						BOOST_CHECK_EQUAL( meq::zernikeAngularDerivative( l, m, rho, theta ), 0.0 );
				}
}

BOOST_AUTO_TEST_SUITE_END()

// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE( flux_coordinate_tests )

BOOST_AUTO_TEST_CASE( the_radial_coordinate_is_the_square_root_of_the_normalised_flux )
{
	for ( double flux : { 0.0, 1.0e-6, 0.01, 0.25, 0.5, 0.81, 1.0 } )
	{
		checkClose( meq::radiusFromNormalisedFlux( flux ), std::sqrt( flux ), 1.0e-15, "rho = sqrt( Psi_N )" );
		checkClose( meq::normalisedFluxFromRadius( meq::radiusFromNormalisedFlux( flux ) ), flux, 1.0e-15, "round trip" );
	}

	// Rejected rather than clamped: a negative normalised flux means the
	// caller's normalisation is wrong, and a silent zero would put the surface
	// on the magnetic axis.
	BOOST_CHECK_THROW( meq::radiusFromNormalisedFlux( -1.0e-12 ), std::invalid_argument );
	BOOST_CHECK_THROW( meq::normalisedFluxFromRadius( -1.0 ), std::invalid_argument );
}

BOOST_AUTO_TEST_CASE( the_chain_rule_factor_is_one_over_two_rho )
{
	// The factor that gets dropped. Measured on a function whose Psi_N
	// derivative is known in closed form: R_4^0, R_2^0 and R_0^0 combine to
	// give rho^4 = Psi_N^2 exactly, whose derivative with respect to Psi_N is
	// 2 Psi_N. Matching powers of rho in
	// R_4^0 = 6 rho^4 - 6 rho^2 + 1 and R_2^0 = 2 rho^2 - 1,
	//
	//     rho^4 = R_4^0/6 + R_2^0/2 + R_0^0/3
	//
	// -- checked first, so that a failure below is about the chain rule and not
	// about the identity. (It caught a wrong identity on the first run, which is
	// the only reason it is a separate assertion rather than a comment.)
	meq::ZernikeExpansion quartic( 4 );
	quartic.setCoefficient( 4, 0, 1.0/6.0 );
	quartic.setCoefficient( 2, 0, 1.0/2.0 );
	quartic.setCoefficient( 0, 0, 1.0/3.0 );

	for ( double rho : { 0.0, 0.2, 0.5, 0.77, 1.0 } )
		checkClose( quartic( rho, 0.4 ), rho*rho*rho*rho, 1.0e-15, "the rho^4 identity" );

	std::cout << "\n  d/dPsi_N of Psi_N^2, against the exact 2 Psi_N\n"
	          << "     Psi_N        exact     from d/drho\n";

	for ( double flux : { 0.04, 0.16, 0.36, 0.64, 1.0 } )
	{
		double const rho = meq::radiusFromNormalisedFlux( flux );
		double const converted = meq::fluxDerivativeFromRadial( quartic.radialDerivative( rho, 0.4 ), rho );

		std::cout << "   " << std::fixed << std::setprecision( 6 ) << std::setw( 8 ) << flux
		          << "   " << std::setw( 10 ) << 2.0*flux
		          << "   " << std::setw( 12 ) << converted << std::defaultfloat << "\n";

		checkClose( converted, 2.0*flux, 1.0e-13, "d/dPsi_N of Psi_N^2" );
		checkClose( quartic.fluxDerivative( flux, 0.4 ), 2.0*flux, 1.0e-13, "ZernikeExpansion::fluxDerivative" );

		// The two mistakes that produce a plausible answer: dropping the factor
		// altogether, and inverting it. Both are wrong by a factor of 2 rho or
		// its square here, and neither is small.
		BOOST_CHECK_MESSAGE( std::fabs( quartic.radialDerivative( rho, 0.4 ) - 2.0*flux ) > 1.0e-3*std::max( 1.0e-3, 2.0*flux ) || flux > 0.99,
			"at Psi_N = " << flux << " a dropped chain-rule factor would have been indistinguishable, so this row measures nothing" );
	}

	std::cout << std::endl;
}

BOOST_AUTO_TEST_SUITE_END()

// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE( expansion_tests )

BOOST_AUTO_TEST_CASE( the_shared_factor_evaluation_agrees_with_a_mode_by_mode_sum )
{
	// ZernikeExpansion walks the modes one azimuthal order at a time and shares
	// the trigonometric factors down a family, rather than calling zernike() per
	// mode. That is a different code path from everything else tested in this
	// file, so it is pinned against the mode-by-mode sum it replaces -- for the
	// value AND for both derivatives.
	int const maxDegree = 11;
	std::vector<meq::ZernikeMode> const modes = meq::zernikeModes( maxDegree );
	std::vector<double> coefficients( modes.size(), 0.0 );

	for ( std::size_t i = 0; i < coefficients.size(); ++i )
		coefficients[ i ] = std::sin( 1.7*i + 0.4 )/( 1.0 + i );

	meq::ZernikeExpansion const expansion( maxDegree, coefficients );

	double worst = 0.0;

	for ( double rho : { 0.0, 0.05, 0.3, 0.62, 0.9, 1.0 } )
		for ( int j = 0; j < 13; ++j )
		{
			double const theta = 2.0*pi*j/13.0;
			double value = 0.0, dRho = 0.0, dTheta = 0.0;

			for ( std::size_t i = 0; i < modes.size(); ++i )
			{
				value += coefficients[ i ]*meq::zernike( modes[ i ].l, modes[ i ].m, rho, theta );
				dRho += coefficients[ i ]*meq::zernikeRadialDerivative( modes[ i ].l, modes[ i ].m, rho, theta );
				dTheta += coefficients[ i ]*meq::zernikeAngularDerivative( modes[ i ].l, modes[ i ].m, rho, theta );
			}

			double fastValue = 0.0, fastRho = 0.0, fastTheta = 0.0;
			expansion.evaluate( rho, theta, fastValue, fastRho, fastTheta );

			worst = std::max( { worst, std::fabs( fastValue - value ), std::fabs( fastRho - dRho ), std::fabs( fastTheta - dTheta ) } );

			checkClose( fastValue, value, 1.0e-13, "expansion value" );
			checkClose( fastRho, dRho, 1.0e-13, "expansion d/drho" );
			checkClose( fastTheta, dTheta, 1.0e-13, "expansion d/dtheta" );
			checkClose( expansion( rho, theta ), value, 1.0e-13, "operator()" );
			checkClose( expansion.radialDerivative( rho, theta ), dRho, 1.0e-13, "radialDerivative()" );
			checkClose( expansion.angularDerivative( rho, theta ), dTheta, 1.0e-13, "angularDerivative()" );
		}

	std::cout << "\n  shared-factor evaluation against the mode-by-mode sum: worst absolute difference "
	          << std::scientific << std::setprecision( 3 ) << worst << std::defaultfloat << std::endl;
}

BOOST_AUTO_TEST_CASE( truncation_is_a_prefix_of_the_coefficients )
{
	int const maxDegree = 9;
	std::vector<double> coefficients( meq::zernikeModeCount( maxDegree ), 0.0 );

	for ( std::size_t i = 0; i < coefficients.size(); ++i )
		coefficients[ i ] = 1.0 + 0.5*i;

	meq::ZernikeExpansion const expansion( maxDegree, coefficients );

	for ( int degree = 0; degree <= maxDegree; ++degree )
	{
		meq::ZernikeExpansion const truncated = expansion.truncated( degree );

		BOOST_CHECK_EQUAL( truncated.maxDegree(), degree );
		BOOST_CHECK_EQUAL( truncated.coefficients().size(), meq::zernikeModeCount( degree ) );

		for ( std::size_t i = 0; i < truncated.coefficients().size(); ++i )
			BOOST_CHECK_EQUAL( truncated.coefficients()[ i ], coefficients[ i ] );

		BOOST_CHECK_EQUAL( truncated.coefficient( degree, degree ), expansion.coefficient( degree, degree ) );
	}

	BOOST_CHECK_THROW( expansion.truncated( maxDegree + 1 ), std::invalid_argument );
	BOOST_CHECK_THROW( expansion.coefficient( maxDegree + 1, 1 ), std::invalid_argument );

	// ( 4, 1 ) has l - |m| odd, so it is not a mode of the basis and the setter
	// must say so rather than writing into whatever index the arithmetic gives.
	meq::ZernikeExpansion mutableExpansion( maxDegree );
	BOOST_CHECK_THROW( mutableExpansion.setCoefficient( 4, 1, 1.0 ), std::invalid_argument );
	BOOST_CHECK_THROW( mutableExpansion.setCoefficient( maxDegree + 1, 1, 1.0 ), std::invalid_argument );
}

BOOST_AUTO_TEST_SUITE_END()

// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE( spectral_decay_tests )

BOOST_AUTO_TEST_CASE( an_analytic_function_has_geometrically_decaying_coefficients )
{
	// THE MEASUREMENT IN-3 RESTS ON.
	//
	// f = 1/( a - u ) with u = ( 3x + 4y )/5, so u runs over [ -1, 1 ] on the
	// disc and the nearest singularity sits at distance a - 1 = 0.4 outside it.
	// Chosen over an entire function such as exp deliberately: an entire
	// function's coefficients fall FASTER than geometrically, so log( c ) is
	// concave in the degree and a straight-line fit is a poor description of
	// something even better than the claim. A pole just outside the disc gives
	// clean geometric decay, which is the honest thing to fit a line to. It also
	// exercises both the cosine and the sine modes, u being neither aligned with
	// x nor with y.
	double const a = 1.4;
	DiscFunction const f = []( double x, double y ) { return 1.0/( 1.4 - ( 3.0*x + 4.0*y )/5.0 ); };

	int const maxDegree = 24;
	meq::ZernikeExpansion const expansion = project( f, maxDegree, 64, 128 );
	std::vector<double> const envelope = degreeEnvelope( expansion );

	std::cout << "\n  Zernike coefficients of 1/( " << a << " - u ), u = ( 3x + 4y )/5\n"
	          << "     l   max |c_lm|      ratio     truncation error\n";

	std::vector<double> degrees, logEnvelope, logDegrees, logEnvelopeForAlgebraic;

	for ( int l = 0; l <= maxDegree; ++l )
	{
		double const current = envelope[ static_cast<std::size_t>( l ) ];
		double const ratio = ( l == 0 ) ? 0.0 : current/envelope[ static_cast<std::size_t>( l - 1 ) ];

		std::cout << "  " << std::setw( 4 ) << l
		          << std::scientific << std::setprecision( 3 ) << "   " << std::setw( 11 ) << current;

		if ( l > 0 )
			std::cout << "   " << std::fixed << std::setprecision( 4 ) << std::setw( 8 ) << ratio;
		else
			std::cout << "           ";

		if ( l % 4 == 0 )
			std::cout << "   " << std::scientific << std::setprecision( 3 ) << std::setw( 12 ) << maximumError( expansion.truncated( l ), f );

		std::cout << std::defaultfloat << "\n";

		// The fit range stops where the coefficients reach the quadrature's own
		// noise floor; below about 1e-13 what is being fitted is round-off.
		if ( current > 1.0e-13 )
		{
			degrees.push_back( l );
			logEnvelope.push_back( std::log10( current ) );

			if ( l > 0 )
			{
				logDegrees.push_back( std::log10( static_cast<double>( l ) ) );
				logEnvelopeForAlgebraic.push_back( std::log10( current ) );
			}
		}
	}

	LineFit const geometric = fitLine( degrees, logEnvelope );
	LineFit const algebraic = fitLine( logDegrees, logEnvelopeForAlgebraic );

	std::cout << std::fixed << std::setprecision( 4 )
	          << "\n  geometric model  log10 c = " << geometric.intercept << " + " << geometric.slope << " l"
	          << "        R^2 = " << geometric.rSquared << "\n"
	          << "  algebraic model  log10 c = " << algebraic.intercept << " + " << algebraic.slope << " log10 l"
	          << "  R^2 = " << algebraic.rSquared << "\n"
	          << "  measured decay per degree: " << std::pow( 10.0, geometric.slope ) << std::defaultfloat << std::endl;

	// Geometric: a fixed factor per degree, and a good straight line in the
	// degree itself.
	BOOST_CHECK_MESSAGE( geometric.slope < -0.15, "decay per degree was " << std::pow( 10.0, geometric.slope ) << ", expected well below one" );
	BOOST_CHECK_MESSAGE( geometric.rSquared > 0.99, "geometric model R^2 was " << geometric.rSquared );

	// And NOT algebraic -- which is the half of the statement that makes it a
	// measurement rather than a description. The algebraic model must fit worse.
	BOOST_CHECK_MESSAGE( geometric.rSquared > algebraic.rSquared,
		"the geometric model must describe this decay better than an algebraic one: R^2 " << geometric.rSquared << " against " << algebraic.rSquared );

	// The consequence, stated where a consumer will read it: a fit of a few
	// dozen modes IS the function, to round-off.
	BOOST_CHECK_MESSAGE( maximumError( expansion, f ) < 1.0e-7,
		"the degree-" << maxDegree << " expansion was " << maximumError( expansion, f ) << " from the function" );
}

BOOST_AUTO_TEST_CASE( a_function_that_is_not_smooth_at_the_centre_decays_only_algebraically )
{
	// THE CONTROL, and it is the same statement as the axis one from the other
	// side. f = sqrt( x^2 + y^2 ) is a perfectly nice function of rho -- it IS
	// rho -- and it is continuous everywhere on the disc. What it is not is
	// smooth at the CENTRE, where it has a cone point.
	//
	// The Zernike radial polynomials at m = 0 are Legendre polynomials in
	// rho^2, so expanding rho in them is expanding sqrt of the variable, and the
	// coefficients fall like a power of the degree rather than like a power to
	// the degree. Nothing about the quadrature changes between this test and the
	// one above; only the smoothness of the function does.
	DiscFunction const f = []( double x, double y ) { return std::sqrt( x*x + y*y ); };

	int const maxDegree = 24;
	meq::ZernikeExpansion const expansion = project( f, maxDegree, 64, 128 );
	std::vector<double> const envelope = degreeEnvelope( expansion );

	std::cout << "\n  Zernike coefficients of sqrt( x^2 + y^2 ), which is C^0 but not C^1 at the centre\n"
	          << "     l   max |c_lm|      ratio   truncation error\n";

	std::vector<double> degrees, logEnvelope, logDegrees;

	for ( int l = 2; l <= maxDegree; l += 2 )
	{
		double const current = envelope[ static_cast<std::size_t>( l ) ];
		double const ratio = current/envelope[ static_cast<std::size_t>( l - 2 ) ];

		std::cout << "  " << std::setw( 4 ) << l
		          << std::scientific << std::setprecision( 3 ) << "   " << std::setw( 11 ) << current
		          << "   " << std::fixed << std::setprecision( 4 ) << std::setw( 8 ) << ratio;

		if ( l % 8 == 0 )
			std::cout << "   " << std::scientific << std::setprecision( 3 ) << std::setw( 12 ) << maximumError( expansion.truncated( l ), f );

		std::cout << std::defaultfloat << "\n";

		degrees.push_back( l );
		logEnvelope.push_back( std::log10( current ) );
		logDegrees.push_back( std::log10( static_cast<double>( l ) ) );
	}

	// The odd degrees, and every m /= 0, must be identically absent: the function
	// is radially symmetric. A projection that put content there would say the
	// quadrature was leaking between modes.
	double worstSpurious = 0.0;

	for ( std::size_t i = 0; i < expansion.modes().size(); ++i )
		if ( expansion.modes()[ i ].m != 0 )
			worstSpurious = std::max( worstSpurious, std::fabs( expansion.coefficients()[ i ] ) );

	LineFit const geometric = fitLine( degrees, logEnvelope );
	LineFit const algebraic = fitLine( logDegrees, logEnvelope );

	std::cout << std::fixed << std::setprecision( 4 )
	          << "\n  geometric model  R^2 = " << geometric.rSquared
	          << "        algebraic model  log10 c = " << algebraic.intercept << " + " << algebraic.slope << " log10 l"
	          << "  R^2 = " << algebraic.rSquared << "\n"
	          << "  largest coefficient with m /= 0: " << std::scientific << std::setprecision( 3 ) << worstSpurious
	          << std::defaultfloat << std::endl;

	BOOST_CHECK_MESSAGE( worstSpurious < 1.0e-14, "a radially symmetric function put " << worstSpurious << " into a mode with m /= 0" );

	// Algebraic, at a measurable power, and better described that way than
	// geometrically. This is the comparison that gives the word "spectral" in the
	// previous test its content.
	BOOST_CHECK_MESSAGE( algebraic.rSquared > geometric.rSquared,
		"this decay must be better described as algebraic than as geometric: R^2 " << algebraic.rSquared << " against " << geometric.rSquared );
	BOOST_CHECK_MESSAGE( algebraic.slope < -1.0 && algebraic.slope > -5.0,
		"the algebraic exponent was " << algebraic.slope << ", which is not the modest power an endpoint singularity gives" );

	// And the consequence: twenty-four degrees buys three or four digits here
	// where it bought fourteen for the analytic function above.
	double const finalError = maximumError( expansion, f );
	std::cout << "  degree-" << maxDegree << " truncation error: " << std::scientific << std::setprecision( 3 ) << finalError << std::defaultfloat << std::endl;

	BOOST_CHECK_MESSAGE( finalError > 1.0e-6,
		"this function was expected to be hard for the basis; if a degree-" << maxDegree << " fit now reaches " << finalError << " the contrast above is empty" );
}

BOOST_AUTO_TEST_SUITE_END()

// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE( axis_regularity_tests )

BOOST_AUTO_TEST_CASE( an_expansion_is_single_valued_at_the_centre_whatever_its_coefficients )
{
	// STRUCTURAL, NOT A CONSEQUENCE OF THE FIT. R_l^m( 0 ) = 0 for every m /= 0,
	// so at rho = 0 only the m = 0 modes survive and the value cannot depend on
	// theta. This is the property a naive "polynomial in rho times Fourier in
	// theta" basis does not have, and it is why the magnetic axis needs no
	// special case downstream.
	int const maxDegree = 10;
	std::vector<double> coefficients( meq::zernikeModeCount( maxDegree ), 0.0 );

	for ( std::size_t i = 0; i < coefficients.size(); ++i )
		coefficients[ i ] = std::cos( 2.3*i );

	meq::ZernikeExpansion const expansion( maxDegree, coefficients );

	double const reference = expansion( 0.0, 0.0 );
	double worst = 0.0;

	for ( int j = 0; j < 64; ++j )
		worst = std::max( worst, std::fabs( expansion( 0.0, 2.0*pi*j/64.0 ) - reference ) );

	std::cout << "\n  spread of a random degree-" << maxDegree << " expansion over theta at rho = 0: "
	          << std::scientific << std::setprecision( 3 ) << worst << std::defaultfloat << std::endl;

	BOOST_CHECK_MESSAGE( worst == 0.0, "the value at the centre depended on theta by " << worst );
}

BOOST_AUTO_TEST_CASE( a_smooth_function_is_represented_smoothly_through_the_centre )
{
	// The control the axis claim needs: something smooth in CARTESIAN
	// coordinates on the disc -- exp( 0.7x - 0.4y ), which is entire -- fitted
	// and then interrogated at rho = 0, where a basis with a coordinate
	// singularity would show one.
	DiscFunction const f = []( double x, double y ) { return std::exp( 0.7*x - 0.4*y ); };

	int const maxDegree = 16;
	meq::ZernikeExpansion const expansion = project( f, maxDegree, 64, 128 );

	double const fitError = maximumError( expansion, f );
	std::cout << "\n  exp( 0.7x - 0.4y ) fitted to degree " << maxDegree
	          << ": worst error on the disc " << std::scientific << std::setprecision( 3 ) << fitError << std::defaultfloat << std::endl;
	BOOST_CHECK_MESSAGE( fitError < 1.0e-12, "the fit itself was only good to " << fitError );

	// The value at the centre.
	checkClose( expansion( 0.0, 0.0 ), 1.0, 1.0e-12, "the fitted value on the axis" );

	// And the Cartesian gradient there, by a central difference of the
	// representation. For a polynomial in ( x, y ) this converges; for a
	// representation with a branch at the centre it would not.
	double const h = 1.0e-4;
	auto const cartesian = [ &expansion ]( double x, double y )
	{
		return expansion( std::hypot( x, y ), std::atan2( y, x ) );
	};

	double const dx = ( cartesian( h, 0.0 ) - cartesian( -h, 0.0 ) )/( 2.0*h );
	double const dy = ( cartesian( 0.0, h ) - cartesian( 0.0, -h ) )/( 2.0*h );

	std::cout << "  gradient at the axis: ( " << std::fixed << std::setprecision( 10 ) << dx << ", " << dy
	          << " ) against the exact ( 0.7, -0.4 )" << std::defaultfloat << std::endl;

	checkClose( dx, 0.7, 1.0e-7, "d/dx at the axis" );
	checkClose( dy, -0.4, 1.0e-7, "d/dy at the axis" );

	// The sharper statement: the representation restricted to a line through the
	// centre is smooth THROUGH it, so the one-sided second derivatives from
	// either side agree in the limit. See the note beside
	// oneSidedSecondDifferenceJump().
	PolarFunction const g = [ &expansion ]( double rho, double theta ) { return expansion( rho, theta ); };

	std::cout << "  second-derivative jump across the centre, along the line at 0.3 rad:\n"
	          << "         h           jump      ratio\n";

	double previous = 0.0;
	double worstRatio = 0.0;

	for ( int i = 0; i < 3; ++i )
	{
		double const step = 1.0e-2/std::pow( 2.0, i );
		double const jump = oneSidedSecondDifferenceJump( g, 0.3, step );

		std::cout << "   " << std::scientific << std::setprecision( 3 ) << std::setw( 10 ) << step
		          << "   " << std::setw( 12 ) << jump;

		if ( i > 0 )
		{
			double const ratio = previous/jump;
			std::cout << "   " << std::fixed << std::setprecision( 3 ) << std::setw( 8 ) << ratio;
			worstRatio = ( worstRatio == 0.0 ) ? ratio : std::min( worstRatio, ratio );
		}

		std::cout << std::defaultfloat << "\n";
		previous = jump;
	}

	std::cout << std::endl;

	// O( h ): halving the step halves the jump. That is what "smooth through the
	// centre" reduces to numerically.
	BOOST_CHECK_MESSAGE( worstRatio > 1.6, "the jump was expected to halve with the step; worst ratio was " << worstRatio );
}

BOOST_AUTO_TEST_CASE( violating_the_parity_constraint_breaks_smoothness_at_the_centre )
{
	// THE CONTRAST THAT MAKES THE WHOLE CHOICE OF BASIS A MEASUREMENT.
	//
	// rho^2 cos( theta ) is what a naive tensor product of "polynomial in rho"
	// and "Fourier in theta" would happily admit: l = 2, m = 1, l - |m| odd. In
	// Cartesian coordinates it is x sqrt( x^2 + y^2 ), which is C^1 at the origin
	// and not C^2 -- its second derivative reads 3 cos( theta ) - cos^3( theta ),
	// a different number down every ray.
	//
	// The library refuses to build it, so the test builds it by hand and measures
	// what the refusal is worth: the second-derivative jump across the centre
	// does not shrink with the step at all, where every admissible mode's does.
	BOOST_CHECK_THROW( meq::zernike( 2, 1, 0.5, 0.3 ), std::invalid_argument );

	PolarFunction const inadmissible = []( double rho, double theta ) { return rho*rho*std::cos( theta ); };

	// Three admissible neighbours of it, chosen to bracket the case: the same
	// angular order at the admissible degrees either side, and one m = 0 mode.
	struct Mode { int l; int m; char const * name; };
	std::vector<Mode> const admissible = { { 1, 1, "Z_1^1" }, { 3, 1, "Z_3^1" }, { 2, 0, "Z_2^0" }, { 4, 2, "Z_4^2" }, { 5, 3, "Z_5^3" } };

	double const alpha = 0.3;

	std::cout << "\n  one-sided second-derivative jump across the disc centre, along the line at "
	          << alpha << " rad\n"
	          << "  mode          h = 1e-2      h = 5e-3     h = 2.5e-3      verdict\n";

	for ( Mode const & mode : admissible )
	{
		PolarFunction const g = [ &mode ]( double rho, double theta ) { return meq::zernike( mode.l, mode.m, rho, theta ); };

		double const coarse = oneSidedSecondDifferenceJump( g, alpha, 1.0e-2 );
		double const middle = oneSidedSecondDifferenceJump( g, alpha, 5.0e-3 );
		double const fine = oneSidedSecondDifferenceJump( g, alpha, 2.5e-3 );

		std::cout << "  " << std::setw( 8 ) << mode.name
		          << std::scientific << std::setprecision( 3 )
		          << "   " << std::setw( 11 ) << coarse
		          << "   " << std::setw( 11 ) << middle
		          << "   " << std::setw( 11 ) << fine
		          << "     shrinks" << std::defaultfloat << "\n";

		// Either the jump is already at round-off -- the mode is a polynomial of
		// degree two or less along the line, so the second difference is exact --
		// or it falls with the step.
		BOOST_CHECK_MESSAGE( fine < 1.0e-9 || fine < 0.6*coarse,
			mode.name << " jump did not shrink with the step: " << coarse << " -> " << fine );
	}

	double const coarse = oneSidedSecondDifferenceJump( inadmissible, alpha, 1.0e-2 );
	double const middle = oneSidedSecondDifferenceJump( inadmissible, alpha, 5.0e-3 );
	double const fine = oneSidedSecondDifferenceJump( inadmissible, alpha, 2.5e-3 );

	std::cout << "  " << std::setw( 8 ) << "r^2cosT"
	          << std::scientific << std::setprecision( 3 )
	          << "   " << std::setw( 11 ) << coarse
	          << "   " << std::setw( 11 ) << middle
	          << "   " << std::setw( 11 ) << fine
	          << "     DOES NOT" << std::defaultfloat << std::endl;

	// The exact value is 4 |cos alpha|, at every step. Asserted as such, because
	// "does not converge" is a weaker and vaguer claim than "sits at a computable
	// constant".
	double const expected = 4.0*std::fabs( std::cos( alpha ) );

	checkClose( coarse, expected, 1.0e-8, "the inadmissible mode's jump at h = 1e-2" );
	checkClose( middle, expected, 1.0e-8, "the inadmissible mode's jump at h = 5e-3" );
	checkClose( fine, expected, 1.0e-8, "the inadmissible mode's jump at h = 2.5e-3" );

	std::cout << "  the inadmissible mode's jump is 4 |cos( " << alpha << " )| = "
	          << std::fixed << std::setprecision( 6 ) << expected
	          << " at every step, which is what a coordinate singularity at the magnetic axis looks like"
	          << std::defaultfloat << std::endl;
}

BOOST_AUTO_TEST_SUITE_END()
