#include "Zernike.hpp"

#include <boost/math/special_functions/jacobi.hpp>

#include <cmath>
#include <stdexcept>
#include <string>
#include <utility>

/*
 * HOW THE RADIAL POLYNOMIAL IS EVALUATED: BOOST.MATH'S JACOBI POLYNOMIAL, AND
 * NOT A HAND-ROLLED ANYTHING.
 *
 * R_l^m IS a Jacobi polynomial under a change of variable. The identity is
 *
 *     R_l^m( rho ) = ( -1 )^k rho^|m| P_k^( |m|, 0 )( 1 - 2 rho^2 ),
 *                    k = ( l - |m| )/2
 *
 * so writing a recurrence here would be reimplementing a well-tested special
 * function and taking on its accuracy characteristics and its maintenance for
 * nothing. boost::math::jacobi() and boost::math::jacobi_prime() are what this
 * file calls.
 *
 * TWO CONVENTIONS IN THAT IDENTITY ARE EASY TO GET WRONG AND BOTH ARE SILENT
 * FOR THE FIRST FEW MODES.
 *
 *   * THE ARGUMENT IS 1 - 2 rho^2, NOT 2 rho^2 - 1.
 *   * THERE IS A LEADING ( -1 )^k.
 *
 * Either mistake reproduces R_0^0, R_1^1, R_2^2 and every R_l^l -- the pure
 * powers, where k = 0 and the Jacobi factor is the constant one -- so the
 * cheapest case that discriminates is l = 3, m = 1. The correct form gives
 * 3 rho^3 - 2 rho; swapping the argument gives rho - 3 rho^3 and dropping the
 * sign gives 2 rho - 3 rho^3. Neither of those is 1 at rho = 1, so the
 * normalisation R_l^m( 1 ) = 1 is the other cheap discriminator.
 * ZernikeTests.cpp asserts both, and evaluates the two rival conventions
 * alongside so the discrimination is measured rather than asserted.
 *
 * WHAT MUST NOT BE USED IS THE EXPLICIT FACTORIAL SUM,
 *
 *     R_l^m( rho ) = sum_{s=0}^{k} ( -1 )^s ( l - s )!
 *                    / [ s! ( (l+m)/2 - s )! ( (l-m)/2 - s )! ] rho^( l - 2s ),
 *
 * which is the definition most references print. Its terms are binomial-sized
 * while its answer is O( 1 ) -- at l = 30 the largest is about 1e10 -- so it
 * costs roughly log10( largest term ) digits to cancellation. Measured in
 * ZernikeTests.cpp at rho = 0.83, it disagrees with the Jacobi route by 4.5e-8
 * at l = 30 and 1.3e-4 at l = 40, and that sum is kept there as a control
 * rather than as an implementation for exactly that reason. A flux-surface fit
 * wanting twenty or thirty modes would be reading noise.
 *
 * THE DERIVATIVE IS EXACT, WITH NO DIFFERENCING ANYWHERE. jacobi_prime() is the
 * derivative with respect to the Jacobi argument, so with x = 1 - 2 rho^2 and
 * dx/drho = -4 rho the chain rule gives
 *
 *     dR_l^m/drho = ( -1 )^k [ m rho^( m - 1 ) P_k( x ) - 4 rho^( m + 1 ) P_k'( x ) ].
 *
 * There is no singularity in that at rho = 0: for m >= 1 the power rho^( m - 1 )
 * is perfectly finite, and for m = 0 the factor m kills the term. It is written
 * as an explicit branch on m = 0 below rather than relying on the product,
 * because at rho = 0 and m = 0 the factor rho^( -1 ) is an infinity and
 * 0 * infinity is a NaN -- so the "obviously zero" term would poison the answer
 * at exactly the point this whole basis exists to make ordinary.
 *
 * Boost.Math is a header-only dependency and it is confined to this .cpp:
 * Zernike.hpp includes nothing but <cstddef> and <vector>, so a consumer of the
 * basis takes on no dependency of its own. Nothing here includes MFEM; see the
 * header for why that matters.
 */

namespace
{

	// Exact for the small non-negative exponents a Zernike degree produces, and
	// unlike std::pow it has no special cases to think about at rho = 0.
	double integerPower( double base, int exponent )
	{
		double result = 1.0;

		for ( int i = 0; i < exponent; ++i )
			result *= base;

		return result;
	}

	// Reject an index pair that would otherwise produce a plausible-looking
	// answer for a mode that is not in the basis. The caller's name is passed in
	// so the message says which entry point was handed the bad pair, in the style
	// of src/meq/Profiles.cpp.
	void checkMode( int l, int m, char const * where )
	{
		if ( l < 0 )
			throw std::invalid_argument( std::string( where ) + ": the Zernike degree l must be non-negative, got l = " + std::to_string( l ) );

		int const absM = ( m < 0 ) ? -m : m;

		if ( absM > l )
			throw std::invalid_argument( std::string( where ) + ": a Zernike mode needs |m| <= l, got l = " + std::to_string( l ) + ", m = " + std::to_string( m ) );

		// The parity condition, and the one worth an explanation in the message
		// rather than a bare "invalid": a caller who has just written a loop over
		// every ( l, m ) with |m| <= l has not made a typographical error, they
		// have made the modelling error this basis exists to prevent, and the
		// message is the only place they will be told so.
		if ( ( l - absM ) % 2 != 0 )
			throw std::invalid_argument( std::string( where ) + ": a Zernike mode needs l - |m| even, got l = " + std::to_string( l ) + ", m = " + std::to_string( m )
				+ " -- an odd difference is a mode that is not smooth at the disc centre, which for meq is the magnetic axis" );
	}

	void checkDegree( int maxDegree, char const * where )
	{
		if ( maxDegree < 0 )
			throw std::invalid_argument( std::string( where ) + ": the maximum degree must be non-negative, got " + std::to_string( maxDegree ) );

		if ( maxDegree > meq::maxZernikeDegree )
			throw std::invalid_argument( std::string( where ) + ": the maximum degree must not exceed meq::maxZernikeDegree = "
				+ std::to_string( meq::maxZernikeDegree ) + ", got " + std::to_string( maxDegree ) );
	}

	/// R_l^m( rho ) together with its rho-derivative, which are wanted together
	/// often enough -- every expansion evaluation wants both -- to be worth
	/// returning as a pair rather than computing the Jacobi polynomial twice.
	struct RadialValue
	{
		double value;
		double derivative;
	};

	/// The identity at the top of this file, evaluated. @a absM is |m|, which is
	/// all the radial polynomial depends on; the sign of m selects the angular
	/// function and is applied by the caller.
	RadialValue radialValue( int absM, int l, double rho )
	{
		unsigned const order = static_cast<unsigned>( ( l - absM )/2 );
		double const alpha = static_cast<double>( absM );
		double const argument = 1.0 - 2.0*rho*rho;

		double const polynomial = boost::math::jacobi( order, alpha, 0.0, argument );
		double const polynomialPrime = boost::math::jacobi_prime( order, alpha, 0.0, argument );

		double const parity = ( order % 2 == 0 ) ? 1.0 : -1.0;
		double const rhoPower = integerPower( rho, absM );

		// The m = 0 branch is explicit rather than arithmetic; see the header of
		// this file for why 0 * rho^( -1 ) is not an acceptable way to write zero.
		double const fromPower = ( absM == 0 ) ? 0.0 : absM*integerPower( rho, absM - 1 )*polynomial;

		return RadialValue{ parity*rhoPower*polynomial,
		                    parity*( fromPower - 4.0*rho*rhoPower*polynomialPrime ) };
	}

	// The angular factor and its theta-derivative. m >= 0 is the cosine branch,
	// m < 0 the sine branch; see the note on ZernikeMode about the sign.
	double angularFactor( int m, double theta )
	{
		return ( m >= 0 ) ? std::cos( m*theta ) : std::sin( -m*theta );
	}

	double angularFactorPrime( int m, double theta )
	{
		return ( m >= 0 ) ? -m*std::sin( m*theta ) : -m*std::cos( -m*theta );
	}

}

namespace meq
{

	bool operator==( ZernikeMode left, ZernikeMode right )
	{
		return left.l == right.l && left.m == right.m;
	}

	bool operator!=( ZernikeMode left, ZernikeMode right )
	{
		return !( left == right );
	}

	bool isValidZernikeMode( int l, int m )
	{
		int const absM = ( m < 0 ) ? -m : m;

		return l >= 0 && absM <= l && ( l - absM ) % 2 == 0;
	}

	std::size_t zernikeModeCount( int maxDegree )
	{
		checkDegree( maxDegree, "meq::zernikeModeCount" );

		return static_cast<std::size_t>( maxDegree + 1 )*static_cast<std::size_t>( maxDegree + 2 )/2;
	}

	std::vector<ZernikeMode> zernikeModes( int maxDegree )
	{
		checkDegree( maxDegree, "meq::zernikeModes" );

		std::vector<ZernikeMode> result;
		result.reserve( zernikeModeCount( maxDegree ) );

		// Degree ascending, and within a degree m from -l to +l in steps of two.
		// That is the ANSI/OSA order, and it is what makes the first
		// zernikeModeCount( d ) entries be exactly zernikeModes( d ).
		for ( int l = 0; l <= maxDegree; ++l )
			for ( int m = -l; m <= l; m += 2 )
				result.push_back( ZernikeMode{ l, m } );

		return result;
	}

	std::size_t zernikeModeIndex( int l, int m )
	{
		checkMode( l, m, "meq::zernikeModeIndex" );

		// ( l( l + 2 ) + m )/2. The numerator is always even: l( l + 2 ) has the
		// parity of l, and l - |m| even makes m have the parity of l too.
		return static_cast<std::size_t>( ( l*( l + 2 ) + m )/2 );
	}

	double zernikeRadial( int l, int m, double rho )
	{
		checkMode( l, m, "meq::zernikeRadial" );

		return radialValue( ( m < 0 ) ? -m : m, l, rho ).value;
	}

	double zernikeRadialPrime( int l, int m, double rho )
	{
		checkMode( l, m, "meq::zernikeRadialPrime" );

		return radialValue( ( m < 0 ) ? -m : m, l, rho ).derivative;
	}

	double zernike( int l, int m, double rho, double theta )
	{
		checkMode( l, m, "meq::zernike" );

		return zernikeRadial( l, m, rho )*angularFactor( m, theta );
	}

	double zernikeRadialDerivative( int l, int m, double rho, double theta )
	{
		checkMode( l, m, "meq::zernikeRadialDerivative" );

		return zernikeRadialPrime( l, m, rho )*angularFactor( m, theta );
	}

	double zernikeAngularDerivative( int l, int m, double rho, double theta )
	{
		checkMode( l, m, "meq::zernikeAngularDerivative" );

		return zernikeRadial( l, m, rho )*angularFactorPrime( m, theta );
	}

	double zernikeRadialNormSquared( int l, int m )
	{
		checkMode( l, m, "meq::zernikeRadialNormSquared" );

		return 1.0/( 2.0*( l + 1.0 ) );
	}

	double zernikeNormSquared( int l, int m )
	{
		checkMode( l, m, "meq::zernikeNormSquared" );

		// The angular integral is 2 pi at m = 0 and pi otherwise. That factor of
		// two is the one worth being careful about; see the header.
		double const angular = ( m == 0 ) ? 2.0*M_PI : M_PI;

		return angular*zernikeRadialNormSquared( l, m );
	}

	double radiusFromNormalisedFlux( double normalisedFlux )
	{
		if ( !std::isfinite( normalisedFlux ) || normalisedFlux < 0.0 )
			throw std::invalid_argument( "meq::radiusFromNormalisedFlux: the normalised flux must be finite and non-negative, got "
				+ std::to_string( normalisedFlux ) );

		return std::sqrt( normalisedFlux );
	}

	double normalisedFluxFromRadius( double rho )
	{
		if ( !std::isfinite( rho ) || rho < 0.0 )
			throw std::invalid_argument( "meq::normalisedFluxFromRadius: the radius must be finite and non-negative, got " + std::to_string( rho ) );

		return rho*rho;
	}

	double fluxDerivativeFromRadial( double radialDerivative, double rho )
	{
		// The factor. Not 1, not 2 rho. See the header.
		return radialDerivative/( 2.0*rho );
	}

	ZernikeExpansion::ZernikeExpansion( int maxDegree )
		: degree( maxDegree ),
		  modeList( zernikeModes( maxDegree ) ),
		  coefficientList( zernikeModeCount( maxDegree ), 0.0 )
	{
	}

	ZernikeExpansion::ZernikeExpansion( int maxDegree, std::vector<double> coefficients )
		: degree( maxDegree ),
		  modeList( zernikeModes( maxDegree ) ),
		  coefficientList( std::move( coefficients ) )
	{
		if ( coefficientList.size() != modeList.size() )
			throw std::invalid_argument( "meq::ZernikeExpansion: degree " + std::to_string( maxDegree ) + " needs "
				+ std::to_string( modeList.size() ) + " coefficients, got " + std::to_string( coefficientList.size() ) );
	}

	void ZernikeExpansion::evaluate( double rho, double theta, double & value, double & dRho, double & dTheta ) const
	{
		value = 0.0;
		dRho = 0.0;
		dTheta = 0.0;

		// One pass per |m|, so that the two trigonometric functions are evaluated
		// once for a whole family of degrees rather than once per mode, and the
		// cosine and sine coefficients of a degree are taken together. It
		// allocates nothing and holds no state, so an expansion may be evaluated
		// from several threads at once.
		for ( int absM = 0; absM <= degree; ++absM )
		{
			double const cosine = std::cos( absM*theta );
			double const sine = std::sin( absM*theta );

			for ( int l = absM; l <= degree; l += 2 )
			{
				RadialValue const radialTerm = radialValue( absM, l, rho );
				double const radial = radialTerm.value;
				double const radialPrime = radialTerm.derivative;

				double const cosineCoefficient = coefficientList[ static_cast<std::size_t>( ( l*( l + 2 ) + absM )/2 ) ];
				double const sineCoefficient = ( absM == 0 ) ? 0.0 : coefficientList[ static_cast<std::size_t>( ( l*( l + 2 ) - absM )/2 ) ];

				double const angular = cosineCoefficient*cosine + sineCoefficient*sine;
				double const angularPrime = absM*( sineCoefficient*cosine - cosineCoefficient*sine );

				value += radial*angular;
				dRho += radialPrime*angular;
				dTheta += radial*angularPrime;
			}
		}
	}

	double ZernikeExpansion::operator()( double rho, double theta ) const
	{
		double value = 0.0;
		double dRho = 0.0;
		double dTheta = 0.0;

		evaluate( rho, theta, value, dRho, dTheta );

		return value;
	}

	double ZernikeExpansion::radialDerivative( double rho, double theta ) const
	{
		double value = 0.0;
		double dRho = 0.0;
		double dTheta = 0.0;

		evaluate( rho, theta, value, dRho, dTheta );

		return dRho;
	}

	double ZernikeExpansion::angularDerivative( double rho, double theta ) const
	{
		double value = 0.0;
		double dRho = 0.0;
		double dTheta = 0.0;

		evaluate( rho, theta, value, dRho, dTheta );

		return dTheta;
	}

	double ZernikeExpansion::fluxDerivative( double normalisedFlux, double theta ) const
	{
		double const rho = radiusFromNormalisedFlux( normalisedFlux );

		return fluxDerivativeFromRadial( radialDerivative( rho, theta ), rho );
	}

	int ZernikeExpansion::maxDegree() const
	{
		return degree;
	}

	std::vector<ZernikeMode> const & ZernikeExpansion::modes() const
	{
		return modeList;
	}

	std::vector<double> const & ZernikeExpansion::coefficients() const
	{
		return coefficientList;
	}

	double ZernikeExpansion::coefficient( int l, int m ) const
	{
		if ( l > degree )
			throw std::invalid_argument( "meq::ZernikeExpansion::coefficient: degree " + std::to_string( l ) + " exceeds the expansion's "
				+ std::to_string( degree ) );

		return coefficientList[ zernikeModeIndex( l, m ) ];
	}

	void ZernikeExpansion::setCoefficient( int l, int m, double value )
	{
		if ( l > degree )
			throw std::invalid_argument( "meq::ZernikeExpansion::setCoefficient: degree " + std::to_string( l ) + " exceeds the expansion's "
				+ std::to_string( degree ) );

		coefficientList[ zernikeModeIndex( l, m ) ] = value;
	}

	ZernikeExpansion ZernikeExpansion::truncated( int newDegree ) const
	{
		if ( newDegree < 0 || newDegree > degree )
			throw std::invalid_argument( "meq::ZernikeExpansion::truncated: degree " + std::to_string( newDegree )
				+ " is not in [ 0, " + std::to_string( degree ) + " ]" );

		std::size_t const wanted = zernikeModeCount( newDegree );

		return ZernikeExpansion( newDegree, std::vector<double>( coefficientList.begin(), coefficientList.begin() + static_cast<std::ptrdiff_t>( wanted ) ) );
	}

}
