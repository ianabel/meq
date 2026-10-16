#include "RotatingSource.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

namespace
{

	// Charge neutrality on R = R_ref is checked by sampling, because the densities
	// are profiles and the constraint has to hold at every psi. meq::Profile
	// documents itself on [ 0, 1 ], so that is where the samples go.
	int const neutralitySamples = 21;
	double const neutralityTolerance = 1.0e-10;

	// The root find for (97). The residual is compared against the sum of the
	// magnitudes it cancels, so the test is relative and does not care what units
	// the densities are in.
	int const potentialIterations = 300;
	double const potentialTolerance = 1.0e-14;

	void requireProfile( std::shared_ptr<meq::Profile const> const & profile, char const * what )
	{
		if ( !profile )
			throw std::invalid_argument( std::string( "meq::RotatingSource: " ) + what + " must not be null" );
	}

	// A profile that is a fixed linear combination of others, exact at all three
	// derivative levels rather than differenced. Private to this file: nothing
	// outside neutralisingDensity() should be able to build one.
	class CombinationProfile : public meq::Profile
	{
		public:
			CombinationProfile( std::vector<double> weights, std::vector<std::shared_ptr<meq::Profile const>> parts )
				: weightData( std::move( weights ) ), partData( std::move( parts ) )
			{
			}

			double operator()( double psi ) const override
			{
				double sum = 0.0;
				for ( std::size_t i = 0; i < partData.size(); ++i )
					sum += weightData[ i ]*( *partData[ i ] )( psi );

				return sum;
			}

			double prime( double psi ) const override
			{
				double sum = 0.0;
				for ( std::size_t i = 0; i < partData.size(); ++i )
					sum += weightData[ i ]*partData[ i ]->prime( psi );

				return sum;
			}

			double doublePrime( double psi ) const override
			{
				double sum = 0.0;
				for ( std::size_t i = 0; i < partData.size(); ++i )
					sum += weightData[ i ]*partData[ i ]->doublePrime( psi );

				return sum;
			}

		private:
			std::vector<double> weightData;
			std::vector<std::shared_ptr<meq::Profile const>> partData;
	};

	// R_0 = sqrt( 2W ) over a CardinalSplineProfile of W; see
	// meq::crossingThroughKnots().
	class SquaredCrossingProfile : public meq::Profile
	{
		public:
			explicit SquaredCrossingProfile( meq::CardinalSplineProfile half )
				: halfSquare( std::move( half ) )
			{
			}

			double operator()( double psi ) const override
			{
				return std::sqrt( 2.0*halfSquare( psi ) );
			}

			double prime( double psi ) const override
			{
				return halfSquare.prime( psi )/( *this )( psi );
			}

			double doublePrime( double psi ) const override
			{
				double const r = ( *this )( psi );
				double const slope = halfSquare.prime( psi )/r;
				return ( halfSquare.doublePrime( psi ) - slope*slope )/r;
			}

		private:
			meq::CardinalSplineProfile halfSquare;
	};

}

namespace meq
{

	double chargeNeutralityResidual( std::vector<Species> const & species, double psi )
	{
		double sum = 0.0;
		for ( Species const & s : species )
		{
			if ( !s.density )
				throw std::invalid_argument( "meq::chargeNeutralityResidual: a density profile is null" );

			sum += s.charge*( *s.density )( psi );
		}

		return sum;
	}

	std::shared_ptr<Profile const> neutralisingDensity( std::vector<Species> const & species, std::size_t index )
	{
		if ( index >= species.size() )
			throw std::out_of_range( "meq::neutralisingDensity: no such species" );
		if ( species.size() < 2 )
			throw std::invalid_argument( "meq::neutralisingDensity: charge neutrality needs at least two species" );
		if ( species[ index ].charge == 0.0 )
			throw std::invalid_argument( "meq::neutralisingDensity: the nominated species has zero charge, so its density is not determined by neutrality" );

		std::vector<double> weights;
		std::vector<std::shared_ptr<Profile const>> parts;
		for ( std::size_t i = 0; i < species.size(); ++i )
		{
			if ( i == index )
				continue;

			requireProfile( species[ i ].density, "a density profile" );
			weights.push_back( -species[ i ].charge/species[ index ].charge );
			parts.push_back( species[ i ].density );
		}

		return std::make_shared<CombinationProfile const>( std::move( weights ), std::move( parts ) );
	}

	std::shared_ptr<Profile const> crossingThroughKnots(
		std::shared_ptr<CardinalSpline const> basis, std::vector<double> const &crossings )
	{
		if ( !basis )
			throw std::invalid_argument( "meq::crossingThroughKnots: the basis must not be null" );
		if ( crossings.size() != basis->size() )
			throw std::invalid_argument( "meq::crossingThroughKnots: one crossing per knot" );
		std::vector<double> half( crossings.size() );
		for ( std::size_t k = 0; k < crossings.size(); ++k )
		{
			if ( !std::isfinite( crossings[ k ] ) || !( crossings[ k ] > 0.0 ) )
				throw std::invalid_argument( "meq::crossingThroughKnots: a crossing radius must be finite and positive, and knot "
					+ std::to_string( k ) + "'s is " + std::to_string( crossings[ k ] ) );
			half[ k ] = 0.5*crossings[ k ]*crossings[ k ];
		}
		return std::make_shared<SquaredCrossingProfile const>(
			CardinalSplineProfile( std::move( basis ), std::move( half ) ) );
	}

	std::shared_ptr<SplineProfile const> tabulateMidplaneCrossing(
		MidplaneField const &field, double rMin, double rMax,
		double lower, double upper, unsigned int intervals, unsigned int samples )
	{
		if ( !field )
			throw std::invalid_argument( "meq::tabulateMidplaneCrossing: the field is empty" );
		if ( !std::isfinite( rMin ) || !std::isfinite( rMax ) || !( rMin < rMax ) )
			throw std::invalid_argument( "meq::tabulateMidplaneCrossing: [ rMin, rMax ] must be finite and non-empty" );
		if ( !std::isfinite( lower ) || !std::isfinite( upper ) || !( lower < upper ) )
			throw std::invalid_argument( "meq::tabulateMidplaneCrossing: [ lower, upper ] must be finite and non-empty" );
		if ( intervals < 2 || samples < 2 )
			throw std::invalid_argument( "meq::tabulateMidplaneCrossing: at least two intervals and two samples" );

		// THE SAMPLES, which bracket every root and are what monotonicity is
		// judged on.
		std::size_t const count = static_cast<std::size_t>( samples ) + 1;
		std::vector<double> radius( count ), x( count ), slope( count );
		std::vector<char> found( count, 0 );
		for ( std::size_t j = 0; j < count; ++j )
		{
			radius[ j ] = rMin + ( rMax - rMin )*static_cast<double>( j )/samples;
			found[ j ] = field( radius[ j ], x[ j ], slope[ j ] )
			             && std::isfinite( x[ j ] ) && std::isfinite( slope[ j ] );
		}

		// The run carrying x in [ lower, upper ], first to last.
		std::size_t first = count, last = 0;
		for ( std::size_t j = 0; j < count; ++j )
			if ( found[ j ] && x[ j ] >= lower && x[ j ] <= upper )
			{
				first = std::min( first, j );
				last = std::max( last, j );
			}
		if ( first == count )
			throw std::runtime_error( "meq::tabulateMidplaneCrossing: no point of the midplane in [ "
				+ std::to_string( rMin ) + ", " + std::to_string( rMax ) + " ] carries a flux in [ "
				+ std::to_string( lower ) + ", " + std::to_string( upper ) + " ], so no flux surface of the support crosses z = 0" );

		// One sample beyond each end, where there is one in the domain, so the
		// end surfaces are bracketed.
		std::size_t const a = ( first > 0 && found[ first - 1 ] ) ? first - 1 : first;
		std::size_t const b = ( last + 1 < count && found[ last + 1 ] ) ? last + 1 : last;

		/*
		 * THE PRECONDITION. Strictly monotone across the run, without a gap, and
		 * with the field's own slope agreeing in sense at every sample inside --
		 * the second test is what catches a fold the samples straddle, and a
		 * slope of the wrong SIGN, which is what a q read with the wrong
		 * orientation would give.
		 */
		double const sense = x[ b ] > x[ a ] ? 1.0 : -1.0;
		for ( std::size_t j = a; j < b; ++j )
		{
			if ( !found[ j + 1 ] )
				throw std::runtime_error( "meq::tabulateMidplaneCrossing: the midplane leaves the domain at R = "
					+ std::to_string( radius[ j + 1 ] ) + " inside the support, so a flux surface's crossing is not defined there" );
			if ( !( sense*( x[ j + 1 ] - x[ j ] ) > 0.0 ) )
				throw std::runtime_error( "meq::tabulateMidplaneCrossing: psi is NOT MONOTONE along z = 0 across the support -- it turns between R = "
					+ std::to_string( radius[ j ] ) + " and " + std::to_string( radius[ j + 1 ] )
					+ " -- so z = 0 crosses a flux surface more than once and the midplane gauge is not one condition per surface. That is a tokamak's midplane, or a mirror whose midplane flux has folded; use the fixed-radius gauge" );
		}
		for ( std::size_t j = first; j <= last; ++j )
			if ( !( sense*slope[ j ] > 0.0 ) )
				throw std::runtime_error( "meq::tabulateMidplaneCrossing: dpsi/dR read off the field at R = "
					+ std::to_string( radius[ j ] ) + " is " + std::to_string( slope[ j ] )
					+ ", against the sense of psi along the midplane, so the crossing's slope 1/( dpsi/dR ) is not defined there" );

		double const xLo = std::min( x[ a ], x[ b ] );
		double const xHi = std::max( x[ a ], x[ b ] );
		if ( xLo > lower || xHi < upper )
			throw std::runtime_error( "meq::tabulateMidplaneCrossing: the midplane inside the domain reaches fluxes [ "
				+ std::to_string( xLo ) + ", " + std::to_string( xHi ) + " ] and not the whole of [ "
				+ std::to_string( lower ) + ", " + std::to_string( upper ) + " ], so some surface of the support does not cross z = 0 in the domain" );

		/*
		 * EACH KNOT BY A SAFEGUARDED NEWTON on its bracketing samples. Newton
		 * takes the field's slope, which is q's and not exactly the derivative
		 * of psi_h, so it converges fast rather than quadratically; the bracket
		 * is what makes it certain. psi_h is discontinuous across element faces
		 * at its own discretisation error, so the end condition is the bracket's
		 * width as well as the residual's size.
		 */
		std::vector<Knot> knots;
		knots.reserve( intervals + 1 );
		for ( unsigned int i = 0; i <= intervals; ++i )
		{
			double const target = lower + ( upper - lower )*static_cast<double>( i )/intervals;

			// The bracket: the samples either side of target, by bisection over
			// a run that is monotone in the sense found above.
			std::size_t lo = a, hi = b;
			while ( hi - lo > 1 )
			{
				std::size_t const mid = lo + ( hi - lo )/2;
				if ( sense*( x[ mid ] - target ) <= 0.0 )
					lo = mid;
				else
					hi = mid;
			}

			double rLo = radius[ lo ], rHi = radius[ hi ];
			double r = 0.0, value = 0.0, d = 0.0;
			bool done = false;
			if ( x[ lo ] == target )
			{
				r = radius[ lo ];
				d = slope[ lo ];
				done = true;
			}
			else if ( x[ hi ] == target )
			{
				r = radius[ hi ];
				d = slope[ hi ];
				done = true;
			}
			else
			{
				r = rLo + ( rHi - rLo )*( target - x[ lo ] )/( x[ hi ] - x[ lo ] );
			}

			for ( int k = 0; k < 100 && !done; ++k )
			{
				if ( !field( r, value, d ) )
					throw std::runtime_error( "meq::tabulateMidplaneCrossing: the field was not found at R = "
						+ std::to_string( r ) + ", between two samples it was found at" );

				double const residual = value - target;
				if ( std::fabs( residual ) <= 4.0*std::numeric_limits<double>::epsilon()*std::max( 1.0, std::fabs( target ) ) )
					break;

				if ( sense*residual < 0.0 )
					rLo = r;
				else
					rHi = r;
				if ( rHi - rLo <= 4.0*std::numeric_limits<double>::epsilon()*rHi )
					break;

				double next = r - residual/d;
				if ( !( next > rLo && next < rHi ) )
					next = 0.5*( rLo + rHi );
				r = next;
			}

			if ( !done && !field( r, value, d ) )
				throw std::runtime_error( "meq::tabulateMidplaneCrossing: the field was not found at R = " + std::to_string( r ) );
			if ( !( sense*d > 0.0 ) )
				throw std::runtime_error( "meq::tabulateMidplaneCrossing: dpsi/dR at the crossing R = " + std::to_string( r )
					+ " has the wrong sense, so 1/( dpsi/dR ) is not the crossing's slope" );

			knots.push_back( Knot{ target, r, 1.0/d } );
		}

		return std::make_shared<SplineProfile const>( std::move( knots ) );
	}

	RotatingSource::RotatingSource( std::vector<Species> species,
		std::shared_ptr<Profile const> omega,
		std::shared_ptr<Profile const> ggPrime,
		double referenceRadius,
		double mu0,
		Closure closure,
		Gauge gauge )
		: speciesData( std::move( species ) ), omegaProfile( std::move( omega ) ),
		  ggPrimeProfile( std::move( ggPrime ) ), refRadius( referenceRadius ), permeability( mu0 ),
		  closureChoice( closure ), gaugeChoice( gauge )
	{
		if ( speciesData.size() < 2 )
			throw std::invalid_argument( "meq::RotatingSource: quasineutrality needs at least two species" );
		if ( speciesData.size() > maxSpecies )
			throw std::invalid_argument( "meq::RotatingSource: more species than meq::maxSpecies; raise that constant rather than allocating per evaluation" );

		if ( closureChoice == Closure::Automatic )
			closureChoice = speciesData.size() == 2 ? Closure::ClosedForm : Closure::RootFind;
		if ( closureChoice == Closure::ClosedForm && speciesData.size() != 2 )
			throw std::invalid_argument( "meq::RotatingSource: the closed-form closure exists only for two species; use Closure::RootFind" );

		requireProfile( ggPrimeProfile, "the g dg/dpsi profile" );

		if ( gaugeChoice == Gauge::Midplane )
		{
			if ( refRadius != 0.0 )
				throw std::invalid_argument( "meq::RotatingSource: Gauge::Midplane pins phi_0 at each flux surface's own midplane crossing and has no reference radius; pass 0 rather than a radius it would not read" );
		}
		else if ( gaugeChoice == Gauge::FixedCurve )
		{
			if ( refRadius != 0.0 )
				throw std::invalid_argument( "meq::RotatingSource: Gauge::FixedCurve pins phi_0 on a curve R_g( psi ) and has no reference radius; pass 0 rather than a radius it would not read" );
		}
		else if ( !std::isfinite( refRadius ) || refRadius <= 0.0 )
			throw std::invalid_argument( "meq::RotatingSource: the reference radius must be finite and positive" );
		if ( !std::isfinite( permeability ) )
			throw std::invalid_argument( "meq::RotatingSource: mu0 must be finite" );

		bool anyPositive = false;
		bool anyNegative = false;
		for ( Species const & s : speciesData )
		{
			requireProfile( s.temperature, "a temperature profile" );
			requireProfile( s.density, "a density profile" );

			if ( !std::isfinite( s.mass ) || s.mass <= 0.0 )
				throw std::invalid_argument( "meq::RotatingSource: every species mass must be finite and positive" );
			if ( !std::isfinite( s.charge ) || s.charge == 0.0 )
				throw std::invalid_argument( "meq::RotatingSource: every species charge must be finite and non-zero" );

			anyPositive = anyPositive || s.charge > 0.0;
			anyNegative = anyNegative || s.charge < 0.0;
		}

		// Charges of both signs. At two species this is what makes the closed
		// form's denominator Z_1 T_2 - Z_2 T_1 positive; at any number it is what
		// gives (97) a root at all, since the left hand side has to change sign.
		if ( !anyPositive || !anyNegative )
			throw std::invalid_argument( "meq::RotatingSource: the species must carry charges of both signs, or quasineutrality has no solution" );

		// The closures are derived from Sum_s Z_s n_s0 = 0 and are simply wrong
		// without it -- not inaccurate, wrong -- so this is a throw and not a
		// warning. neutralisingDensity() is how a caller builds a conforming set.
		for ( int i = 0; i <= neutralitySamples; ++i )
		{
			double const psi = static_cast<double>( i )/neutralitySamples;

			double scale = 0.0;
			for ( Species const & s : speciesData )
				scale += std::fabs( s.charge*( *s.density )( psi ) );

			double const residual = chargeNeutralityResidual( speciesData, psi );

			if ( std::fabs( residual ) > neutralityTolerance*std::max( 1.0, scale ) )
				throw std::invalid_argument( "meq::RotatingSource: the species violate charge neutrality on the reference curve at psi = "
					+ std::to_string( psi ) + "; sum of Z_s n_s0 is " + std::to_string( residual ) );
		}
	}

	RotatingSource::State RotatingSource::closedFormState( double radius, double psi ) const
	{
		// Two species. Taking logarithms of (97) makes it linear in phi_0, and
		// the two exponents come out EQUAL:
		//
		//     C( psi ) = omega^2 ( Z_1 m_2 - Z_2 m_1 )/( Z_1 T_2 - Z_2 T_1 )
		//     e phi_0  = omega^2 ( R^2 - R_ref^2 )( m_1 T_2 - m_2 T_1 )/2( Z_1 T_2 - Z_2 T_1 )
		//
		// which for ions and electrons ( Z_2 = -1 ) is the familiar
		// omega^2 ( m_i + Z_i m_e )/( T_i + Z_i T_e ). The electron mass is kept
		// rather than dropped: it costs one term and removes a question. That the
		// exponents are equal is what makes Sum_s Z_s n_s vanish at every R once
		// it vanishes at R_ref.
		State state{};

		Species const & one = speciesData[ 0 ];
		Species const & two = speciesData[ 1 ];

		if ( !omegaProfile )
		{
			state.potential = 0.0;
			state.potentialPrime = 0.0;
			state.potentialDoublePrime = 0.0;
			for ( std::size_t s = 0; s < 2; ++s )
			{
				state.exponent[ s ] = 0.0;
				state.exponentPrime[ s ] = 0.0;
				state.exponentDoublePrime[ s ] = 0.0;
			}
			return state;
		}

		double const k = one.charge*two.mass - two.charge*one.mass;

		double const t1 = ( *one.temperature )( psi );
		double const t2 = ( *two.temperature )( psi );
		double const d = one.charge*t2 - two.charge*t1;

		if ( d == 0.0 )
			throw std::invalid_argument( "meq::RotatingSource: Z_1 T_2 - Z_2 T_1 vanishes, so the quasineutrality closure is singular; with opposite charges this needs a temperature to have gone non-positive" );

		double const t1Prime = one.temperature->prime( psi );
		double const t2Prime = two.temperature->prime( psi );
		double const dPrime = one.charge*t2Prime - two.charge*t1Prime;
		double const dDoublePrime = one.charge*two.temperature->doublePrime( psi )
			- two.charge*one.temperature->doublePrime( psi );

		// W = omega^2, so the chain rule below is in one variable rather than two.
		double const w0 = ( *omegaProfile )( psi );
		double const wPrime0 = omegaProfile->prime( psi );
		double const wDoublePrime0 = omegaProfile->doublePrime( psi );

		double const w = w0*w0;
		double const wPrime = 2.0*w0*wPrime0;
		double const wDoublePrime = 2.0*( wPrime0*wPrime0 + w0*wDoublePrime0 );

		double const c = k*w/d;
		double const cPrime = k*( wPrime/d - w*dPrime/( d*d ) );
		double const cDoublePrime = k*( wDoublePrime/d - 2.0*wPrime*dPrime/( d*d ) - w*dDoublePrime/( d*d )
			+ 2.0*w*dPrime*dPrime/( d*d*d ) );

		// ( R^2 - R_*^2 )/2 and its derivatives: constant in psi under the
		// fixed-radius gauge, so the h' and h'' terms below are exactly zero
		// there and every fixed-radius answer is the product it always was.
		double h = 0.0, hPrime = 0.0, hDoublePrime = 0.0;
		halfDelta( radius, psi, h, hPrime, hDoublePrime );

		for ( std::size_t s = 0; s < 2; ++s )
		{
			state.exponent[ s ] = c*h;
			state.exponentPrime[ s ] = cPrime*h + c*hPrime;
			state.exponentDoublePrime[ s ] = cDoublePrime*h + 2.0*cPrime*hPrime + c*hDoublePrime;
		}

		// e phi_0 and its derivatives, from the same denominator. The numerator is
		// a mass-weighted temperature difference: two species with equal m/T leave
		// nothing for the field to separate and phi_0 vanishes identically.
		double const num = one.mass*t2 - two.mass*t1;
		double const numPrime = one.mass*t2Prime - two.mass*t1Prime;
		double const numDoublePrime = one.mass*two.temperature->doublePrime( psi )
			- two.mass*one.temperature->doublePrime( psi );

		double const ratio = num/d;
		double const ratioPrime = numPrime/d - num*dPrime/( d*d );
		double const ratioDoublePrime = numDoublePrime/d - 2.0*numPrime*dPrime/( d*d )
			- num*dDoublePrime/( d*d ) + 2.0*num*dPrime*dPrime/( d*d*d );

		// e phi_0 = u h with u = omega^2 num/d, differentiated as a product.
		double const u = w*ratio;
		double const uPrime = wPrime*ratio + w*ratioPrime;
		double const uDoublePrime = wDoublePrime*ratio + 2.0*wPrime*ratioPrime + w*ratioDoublePrime;

		state.potential = u*h;
		state.potentialPrime = uPrime*h + u*hPrime;
		state.potentialDoublePrime = uDoublePrime*h + 2.0*uPrime*hPrime + u*hDoublePrime;

		return state;
	}

	RotatingSource::State RotatingSource::rootFindState( double radius, double psi ) const
	{
		// The general closure. (97) is transcendental for three or more species,
		// but as well behaved as such a thing gets: its derivative in phi_0 is
		//
		//     -e Sum_s ( Z_s^2/T_s ) n_s  <  0   strictly,
		//
		// every term carrying Z_s^2, so the left hand side is strictly decreasing
		// for positive densities and temperatures. With charges of both signs it
		// runs from +infinity to -infinity, so the root exists, is unique, and can
		// be bracketed -- a safeguarded Newton cannot fail.
		std::size_t const n = speciesData.size();
		double h = 0.0, hPrime = 0.0, hDoublePrime = 0.0;
		halfDelta( radius, psi, h, hPrime, hDoublePrime );

		std::array<double, maxSpecies> t{}, tPrime{}, tDoublePrime{};
		std::array<double, maxSpecies> nRef{}, nRefPrime{}, nRefDoublePrime{};
		double scale = 0.0;

		for ( std::size_t s = 0; s < n; ++s )
		{
			t[ s ] = ( *speciesData[ s ].temperature )( psi );
			tPrime[ s ] = speciesData[ s ].temperature->prime( psi );
			tDoublePrime[ s ] = speciesData[ s ].temperature->doublePrime( psi );

			if ( !( t[ s ] > 0.0 ) )
				throw std::invalid_argument( "meq::RotatingSource: a species temperature is not positive at psi = "
					+ std::to_string( psi ) + ", so (97) has no Maxwellian to solve for" );

			nRef[ s ] = ( *speciesData[ s ].density )( psi );
			nRefPrime[ s ] = speciesData[ s ].density->prime( psi );
			nRefDoublePrime[ s ] = speciesData[ s ].density->doublePrime( psi );

			scale = std::max( scale, t[ s ] );
		}

		double const w0 = omegaProfile ? ( *omegaProfile )( psi ) : 0.0;
		double const wPrime0 = omegaProfile ? omegaProfile->prime( psi ) : 0.0;
		double const wDoublePrime0 = omegaProfile ? omegaProfile->doublePrime( psi ) : 0.0;

		double const w = w0*w0;
		double const wPrime = 2.0*w0*wPrime0;
		double const wDoublePrime = 2.0*( wPrime0*wPrime0 + w0*wDoublePrime0 );

		// The centrifugal half of the exponent, which does not move during the
		// solve, and the value of Sum_s Z_s n_s at a trial potential. Evaluated
		// with the largest exponent factored out, so that a high Mach number
		// cannot overflow the residual: every shifted exponent is at most one, and
		// the common factor cancels out of both the sign test and the Newton step.
		auto shifted = [ & ]( double y, std::array<double, maxSpecies> & value )
		{
			double top = -std::numeric_limits<double>::infinity();
			for ( std::size_t s = 0; s < n; ++s )
			{
				value[ s ] = speciesData[ s ].mass*w*h/t[ s ] - speciesData[ s ].charge*y/t[ s ];
				top = std::max( top, value[ s ] );
			}
			for ( std::size_t s = 0; s < n; ++s )
				value[ s ] = nRef[ s ]*std::exp( value[ s ] - top );
		};

		std::array<double, maxSpecies> v{};

		auto residual = [ & ]( double y )
		{
			shifted( y, v );
			double sum = 0.0;
			for ( std::size_t s = 0; s < n; ++s )
				sum += speciesData[ s ].charge*v[ s ];

			return sum;
		};

		// Bracket. G is decreasing, so the low end is where it is positive.
		double step = std::max( scale, 1.0e-300 );
		double lo = 0.0, hi = 0.0;
		double const atZero = residual( 0.0 );

		if ( atZero > 0.0 )
		{
			lo = 0.0;
			hi = step;
			for ( int i = 0; i < 200 && residual( hi ) > 0.0; ++i )
			{
				lo = hi;
				hi *= 2.0;
			}
		}
		else
		{
			hi = 0.0;
			lo = -step;
			for ( int i = 0; i < 200 && residual( lo ) < 0.0; ++i )
			{
				hi = lo;
				lo *= 2.0;
			}
		}

		double y = 0.5*( lo + hi );
		bool converged = false;
		double previousStep = hi - lo;

		for ( int i = 0; i < potentialIterations; ++i )
		{
			shifted( y, v );

			double numer = 0.0;
			double denom = 0.0;
			double magnitude = 0.0;
			for ( std::size_t s = 0; s < n; ++s )
			{
				double const z = speciesData[ s ].charge;
				numer += z*v[ s ];
				denom += z*z*v[ s ]/t[ s ];
				magnitude += std::fabs( z*v[ s ] );
			}

			if ( std::fabs( numer ) <= potentialTolerance*magnitude )
			{
				converged = true;
				break;
			}

			if ( numer > 0.0 )
				lo = y;
			else
				hi = y;

			// THE BRACKET HAS CLOSED: the root lies between two adjacent
			// doubles. Far from a plasma's reference curve at high Mach, y/T_s is
			// so large that one ulp of y moves the residual by more than the
			// relative tolerance above, which then cannot pass -- MEASURED on
			// examples/mirror-three-species.toml's output grid, with
			// RotatingSourceTests' theRootFindEndsWhenItsBracketDoes the
			// regression. The root is known to the last bit, which is the answer.
			if ( hi - lo <= 4.0*std::numeric_limits<double>::epsilon()*std::max( std::fabs( lo ), std::fabs( hi ) ) )
			{
				converged = true;
				break;
			}

			// Newton on G, which is y + numer/denom because G' = -denom. Bisect
			// whenever that leaves the bracket, which is what makes it safeguarded
			// rather than merely fast -- AND whenever it has not halved the step
			// before it. Where one species' exponential dominates the sum, G'/G
			// is that species' Z/T and Newton walks by T a step across a bracket
			// of thousands of T: inside the bracket every time, and never there.
			double next = y + numer/denom;
			if ( !( next > lo && next < hi ) || 2.0*std::fabs( next - y ) > std::fabs( previousStep ) )
				next = 0.5*( lo + hi );
			previousStep = next - y;

			if ( next == y )
			{
				converged = true;
				break;
			}
			y = next;
		}

		if ( !converged )
			throw std::runtime_error( "meq::RotatingSource: the quasineutrality root find did not converge at psi = "
				+ std::to_string( psi ) + "; the closure is strictly monotone, so this means a profile has gone non-physical" );

		// phi_0's psi-derivatives, by IMPLICIT differentiation of (97). Sum_s Z_s
		// n_s vanishes identically in psi, so its first and second derivatives do
		// too, and each gives one linear equation for one derivative of y against
		// the same denominator the Newton step used. Differencing the root find
		// from outside would instead give a derivative whose accuracy is the inner
		// tolerance -- see the class documentation.
		State state{};
		state.potential = y;

		std::array<double, maxSpecies> b{}, bPrime{};
		double denom = 0.0;
		for ( std::size_t s = 0; s < n; ++s )
			denom += speciesData[ s ].charge*speciesData[ s ].charge*v[ s ]/t[ s ];

		// v[ s ] is n_s0 exp( A_s - max A ), so dividing it by n_s0 recovers the
		// shifted exponential alone -- which is what the numerator wants, since
		// n_s0 appears there explicitly. Guarded because a density may
		// legitimately pass through zero, where the whole term does too.
		double firstNumer = 0.0;
		for ( std::size_t s = 0; s < n; ++s )
		{
			double const z = speciesData[ s ].charge;
			double const m = speciesData[ s ].mass;
			double const ti = t[ s ], tp = tPrime[ s ];

			// B_s = dA_s/dpsi at FIXED y. The centrifugal half is m u_s h with
			// u_s = omega^2/T_s, so its derivative is m( u_s' h + u_s h' ), and
			// h' is the midplane gauge's -R_0 R_0' -- zero under a fixed radius.
			b[ s ] = m*( h*( wPrime/ti - w*tp/( ti*ti ) ) + hPrime*w/ti ) + z*y*tp/( ti*ti );

			double const e = nRef[ s ] == 0.0 ? 0.0 : v[ s ]/nRef[ s ];
			firstNumer += z*( nRefPrime[ s ] + nRef[ s ]*b[ s ] )*e;
		}

		double const yPrime = firstNumer/denom;
		state.potentialPrime = yPrime;

		double secondNumer = 0.0;
		std::array<double, maxSpecies> l{};
		for ( std::size_t s = 0; s < n; ++s )
		{
			double const z = speciesData[ s ].charge;
			double const m = speciesData[ s ].mass;
			double const ti = t[ s ], tp = tPrime[ s ], tpp = tDoublePrime[ s ];

			l[ s ] = b[ s ] - z*yPrime/ti;

			bPrime[ s ] = m*( h*( wDoublePrime/ti - 2.0*wPrime*tp/( ti*ti )
					- w*tpp/( ti*ti ) + 2.0*w*tp*tp/( ti*ti*ti ) )
					+ 2.0*hPrime*( wPrime/ti - w*tp/( ti*ti ) ) + hDoublePrime*w/ti )
				+ z*yPrime*tp/( ti*ti ) + z*y*( tpp/( ti*ti ) - 2.0*tp*tp/( ti*ti*ti ) );

			double const e = nRef[ s ] == 0.0 ? 0.0 : v[ s ]/nRef[ s ];
			secondNumer += z*( nRefDoublePrime[ s ] + 2.0*nRefPrime[ s ]*l[ s ]
				+ nRef[ s ]*( l[ s ]*l[ s ] + bPrime[ s ] ) + nRef[ s ]*z*yPrime*tp/( ti*ti ) )*e;
		}

		double const yDoublePrime = secondNumer/denom;
		state.potentialDoublePrime = yDoublePrime;

		for ( std::size_t s = 0; s < n; ++s )
		{
			double const z = speciesData[ s ].charge;
			double const ti = t[ s ], tp = tPrime[ s ];

			state.exponent[ s ] = speciesData[ s ].mass*w*h/ti - z*y/ti;
			state.exponentPrime[ s ] = l[ s ];
			state.exponentDoublePrime[ s ] = bPrime[ s ] - z*yDoublePrime/ti + z*yPrime*tp/( ti*ti );
		}

		return state;
	}

	void RotatingSource::halfDelta( double radius, double psi, double & value,
		double & prime, double & doublePrime ) const
	{
		if ( gaugeChoice == Gauge::FixedRadius )
		{
			value = ( radius*radius - refRadius*refRadius )/2.0;
			prime = 0.0;
			doublePrime = 0.0;
			return;
		}

		// THE MIDPLANE GAUGE, AND THE FIXED CURVE, WHICH IS THE SAME
		// ALGEBRA ON A TABLE NOBODY REPLACES. h = ( R^2 - R_0^2 )/2 at fixed R, so
		//
		//     h'  = -R_0 R_0'
		//     h'' = -( R_0'^2 + R_0 R_0'' )
		//
		// with R_0', R_0'' the table's own derivatives -- the knots' C^2
		// spline, or a Hermite table whose R_0'' is piecewise linear and jumps
		// at its knots -- exactly the derivatives of the R_0 that h uses, so
		// dFdPsi() is the derivative of f() for the table in hand.
		if ( !midplaneTable && gaugeChoice == Gauge::FixedCurve )
			throw std::logic_error( "meq::RotatingSource: Gauge::FixedCurve has no curve R_g( psi ) yet. GradShafranovSolver supplies the vacuum's midplane crossing through NormalisedSource::setVacuumMidplaneField() before the first residual of a normalised solve; elsewhere call setMidplaneCrossing() first" );
		if ( !midplaneTable )
			throw std::logic_error( "meq::RotatingSource: Gauge::Midplane has no midplane crossing R_0( psi ) yet. GradShafranovSolver supplies it at every residual of the bordered solve, as the interpolant through its crossings' knots, so a midplane-gauged source is solved as a meq::NormalisedRotatingSource; elsewhere call setMidplaneCrossing() first" );

		double const r0 = ( *midplaneTable )( psi );
		double const r0Prime = midplaneTable->prime( psi );
		double const r0DoublePrime = midplaneTable->doublePrime( psi );

		value = ( radius*radius - r0*r0 )/2.0;
		prime = -r0*r0Prime;
		doublePrime = -( r0Prime*r0Prime + r0*r0DoublePrime );
	}

	void RotatingSource::setMidplaneCrossing( std::shared_ptr<Profile const> table )
	{
		if ( gaugeChoice == Gauge::FixedRadius )
			throw std::logic_error( "meq::RotatingSource::setMidplaneCrossing: this source is in Gauge::FixedRadius, which pins phi_0 at a constant radius and reads no crossing" );
		if ( !table )
			throw std::invalid_argument( "meq::RotatingSource::setMidplaneCrossing: the table must not be null" );

		midplaneTable = std::move( table );
	}

	Profile const * RotatingSource::midplaneCrossing() const
	{
		return midplaneTable.get();
	}

	RotatingSource::LabelCheck RotatingSource::labelCheck( std::shared_ptr<Profile const> curve ) const
	{
		if ( gaugeChoice == Gauge::FixedRadius )
			throw std::logic_error( "meq::RotatingSource::labelCheck: this source is in Gauge::FixedRadius, which has no curve to measure on" );
		if ( !curve )
			throw std::invalid_argument( "meq::RotatingSource::labelCheck: the curve must not be null" );

		/*
		 * ON A COPY CARRYING THE CURVE, so this source's own table is neither
		 * read nor replaced. The peak is taken ALONG the curve, where on a
		 * mirror every surface is widest and its centrifugal density largest.
		 * Each label is swept over every radius it can reach, ( 0, R_g( label ) ]
		 * -- a mirror's surface is widest where it crosses z = 0, and the
		 * solved label surface passes through the label point -- because the
		 * solved surface does not lie on the vacuum's off the midplane, and a
		 * profile set meeting the precondition makes F vanish there at EVERY
		 * radius: at two species F's bracket is linear in ( R^2 - R_g^2 )/2,
		 * and both its coefficients must vanish. NOT beyond R_g( label ),
		 * which no label surface reaches and where the centrifugal factor of
		 * a heavy species is astronomical -- MEASURED, 3e+68 of the peak at the
		 * inner label read at the outer one's radius, for profiles whose
		 * densities fall to 1e-6 there.
		 */
		RotatingSource probe( *this );
		probe.midplaneTable = std::move( curve );

		int const along = 256;
		int const across = 128;
		LabelCheck out;
		for ( int i = 0; i <= along; ++i )
		{
			double const x = static_cast<double>( i )/along;
			double const radius = ( *probe.midplaneTable )( x );
			if ( !( radius > 0.0 ) )
				continue;
			out.peak = std::max( out.peak, std::fabs( probe.f( radius, 0.0, x ) )/radius );
		}
		for ( double label : { 0.0, 1.0 } )
			for ( int k = 1; k <= across; ++k )
			{
				double const radius = ( *probe.midplaneTable )( label )*static_cast<double>( k )/across;
				double const value = std::fabs( probe.f( radius, 0.0, label ) )/radius;
				if ( value > out.worst )
				{
					out.worst = value;
					out.label = label;
					out.radius = radius;
				}
			}
		out.relative = out.peak > 0.0 ? out.worst/out.peak
		                              : ( out.worst > 0.0 ? std::numeric_limits<double>::infinity() : 0.0 );
		return out;
	}

	void RotatingSource::crossingDerivatives( double radius, double psi, double &dFdW,
	                                          double &dFdWPrime ) const
	{
		if ( gaugeChoice != Gauge::Midplane )
			throw std::logic_error( "meq::RotatingSource::crossingDerivatives: this source is in Gauge::FixedRadius, whose F reads no crossing" );

		dFdW = 0.0;
		dFdWPrime = 0.0;
		if ( !omegaProfile )
			return;

		/*
		 * dp'/dh AND dp'/dh' AT FIXED ( R, psi ), h = R^2/2 - W.
		 *
		 * p' = Sum_s e^A_s [ ( n_s0' + n_s0 A_s' ) T_s + n_s0 T_s' ], so with
		 * alpha_s = dA_s/dh, and A_s independent of h' -- (97) holds at each
		 * point and h' is not in it --
		 *
		 *     dp'/dh  = Sum_s e^A_s [ alpha_s ( ( n_s0' + n_s0 A_s' ) T_s + n_s0 T_s' )
		 *                             + n_s0 T_s dA_s'/dh ]
		 *     dp'/dh' = Sum_s e^A_s n_s0 T_s dA_s'/dh'
		 *
		 * and the closures differ only in alpha_s, dA_s'/dh and dA_s'/dh'.
		 * W = R_0^2/2 enters as -h, so both derivatives change sign at the end.
		 */
		State const state = stateAt( radius, psi );
		std::size_t const n = speciesData.size();
		std::array<double, maxSpecies> alpha{}, lH{}, lHPrime{};

		double const w0 = ( *omegaProfile )( psi );
		double const wPrime0 = omegaProfile->prime( psi );
		double const w = w0*w0;
		double const wPrime = 2.0*w0*wPrime0;

		if ( closureChoice == Closure::ClosedForm )
		{
			// A = C h and A' = C' h + C h', so alpha = C, dA'/dh = C' and
			// dA'/dh' = C -- closedFormState()'s C, rebuilt here.
			Species const & one = speciesData[ 0 ];
			Species const & two = speciesData[ 1 ];
			double const k = one.charge*two.mass - two.charge*one.mass;
			double const d = one.charge*( *two.temperature )( psi ) - two.charge*( *one.temperature )( psi );
			double const dPrime = one.charge*two.temperature->prime( psi ) - two.charge*one.temperature->prime( psi );
			double const c = k*w/d;
			double const cPrime = k*( wPrime/d - w*dPrime/( d*d ) );
			for ( std::size_t s = 0; s < 2; ++s )
			{
				alpha[ s ] = c;
				lH[ s ] = cPrime;
				lHPrime[ s ] = c;
			}
		}
		else
		{
			/*
			 * THE ROOT FIND, DIFFERENTIATED IMPLICITLY AS rootFindState() IS.
			 *
			 * y = e phi_0 solves G = Sum_s Z_s n_s0 e^A_s = 0 with
			 * A_s = m_s W_s h - Z_s y/T_s, W_s = omega^2/T_s, so at fixed psi
			 *
			 *     dy/dh = Sum_s Z_s n_s m_s W_s / D,   D = Sum_s Z_s^2 n_s/T_s
			 *
			 * and alpha_s = m_s W_s - Z_s y_h/T_s. The exponent's slope is
			 * A_s' = b_s - Z_s y'/T_s with y' = N/D, N = Sum_s Z_s ( n_s0' + n_s0 b_s ) e^A_s
			 * and b_s rootFindState()'s B_s, whose h- and y-dependence is
			 * m_s beta_s h + m_s W_s h' + Z_s y T_s'/T_s^2. Every sum is over
			 * n_s = n_s0 e^A_s, so a common factor cancels from every ratio.
			 */
			std::array<double, maxSpecies> e{}, b{}, nRef{}, nRefPrime{}, t{}, tPrime{}, beta{};
			double denom = 0.0, yH = 0.0;
			for ( std::size_t s = 0; s < n; ++s )
			{
				double const z = speciesData[ s ].charge;
				double const m = speciesData[ s ].mass;
				t[ s ] = ( *speciesData[ s ].temperature )( psi );
				tPrime[ s ] = speciesData[ s ].temperature->prime( psi );
				nRef[ s ] = ( *speciesData[ s ].density )( psi );
				nRefPrime[ s ] = speciesData[ s ].density->prime( psi );
				e[ s ] = std::exp( state.exponent[ s ] );
				beta[ s ] = wPrime/t[ s ] - w*tPrime[ s ]/( t[ s ]*t[ s ] );
				b[ s ] = state.exponentPrime[ s ] + z*state.potentialPrime/t[ s ];
				double const ns = nRef[ s ]*e[ s ];
				denom += z*z*ns/t[ s ];
				yH += z*ns*m*w/t[ s ];
			}
			yH /= denom;

			double dNdH = 0.0, dDdH = 0.0, dNdHPrime = 0.0;
			for ( std::size_t s = 0; s < n; ++s )
			{
				double const z = speciesData[ s ].charge;
				double const m = speciesData[ s ].mass;
				alpha[ s ] = m*w/t[ s ] - z*yH/t[ s ];
				double const bH = m*beta[ s ] + z*tPrime[ s ]*yH/( t[ s ]*t[ s ] );
				dNdH += z*( nRef[ s ]*bH + ( nRefPrime[ s ] + nRef[ s ]*b[ s ] )*alpha[ s ] )*e[ s ];
				dDdH += z*z*nRef[ s ]*e[ s ]*alpha[ s ]/t[ s ];
				dNdHPrime += z*nRef[ s ]*e[ s ]*m*w/t[ s ];
			}
			double const yPrimeH = ( dNdH - state.potentialPrime*dDdH )/denom;
			double const yPrimeHPrime = dNdHPrime/denom;
			for ( std::size_t s = 0; s < n; ++s )
			{
				double const z = speciesData[ s ].charge;
				double const m = speciesData[ s ].mass;
				lH[ s ] = m*beta[ s ] + z*tPrime[ s ]*yH/( t[ s ]*t[ s ] ) - z*yPrimeH/t[ s ];
				lHPrime[ s ] = m*w/t[ s ] - z*yPrimeHPrime/t[ s ];
			}
		}

		double pH = 0.0, pHPrime = 0.0;
		for ( std::size_t s = 0; s < n; ++s )
		{
			double const nRef = ( *speciesData[ s ].density )( psi );
			double const nRefPrime = speciesData[ s ].density->prime( psi );
			double const t = ( *speciesData[ s ].temperature )( psi );
			double const tPrime = speciesData[ s ].temperature->prime( psi );
			double const e = std::exp( state.exponent[ s ] );
			double const l = state.exponentPrime[ s ];

			pH += e*( alpha[ s ]*( ( nRefPrime + nRef*l )*t + nRef*tPrime ) + nRef*t*lH[ s ] );
			pHPrime += e*nRef*t*lHPrime[ s ];
		}

		dFdW = -permeability*radius*radius*pH;
		dFdWPrime = -permeability*radius*radius*pHPrime;
	}

	RotatingSource::State RotatingSource::stateAt( double radius, double psi ) const
	{
		return closureChoice == Closure::ClosedForm ? closedFormState( radius, psi ) : rootFindState( radius, psi );
	}

	void RotatingSource::pressureFrom( State const & state, double psi,
		double & p, double & pPrime, double & pDoublePrime ) const
	{
		// p = Sum_s n_s0 exp( A_s ) T_s, differentiated twice. Written through
		// n_s0 and its derivatives rather than through ( ln n_s0 )', so that a
		// density passing through zero is not a division.
		p = 0.0;
		pPrime = 0.0;
		pDoublePrime = 0.0;

		for ( std::size_t s = 0; s < speciesData.size(); ++s )
		{
			double const nRef = ( *speciesData[ s ].density )( psi );
			double const nRefPrime = speciesData[ s ].density->prime( psi );
			double const nRefDoublePrime = speciesData[ s ].density->doublePrime( psi );

			double const t = ( *speciesData[ s ].temperature )( psi );
			double const tPrime = speciesData[ s ].temperature->prime( psi );
			double const tDoublePrime = speciesData[ s ].temperature->doublePrime( psi );

			double const e = std::exp( state.exponent[ s ] );
			double const l = state.exponentPrime[ s ];
			double const lPrime = state.exponentDoublePrime[ s ];

			double const nS = nRef*e;
			double const nSPrime = ( nRefPrime + nRef*l )*e;
			double const nSDoublePrime = ( nRefDoublePrime + 2.0*nRefPrime*l + nRef*( l*l + lPrime ) )*e;

			p += nS*t;
			pPrime += nSPrime*t + nS*tPrime;
			pDoublePrime += nSDoublePrime*t + 2.0*nSPrime*tPrime + nS*tDoublePrime;
		}
	}

	double RotatingSource::densityExponent( std::size_t index, double radius, double psi ) const
	{
		if ( index >= speciesData.size() )
			throw std::out_of_range( "meq::RotatingSource::densityExponent: no such species" );

		return stateAt( radius, psi ).exponent[ index ];
	}

	double RotatingSource::potential( double radius, double psi ) const
	{
		return stateAt( radius, psi ).potential;
	}

	double RotatingSource::dPotentialDPsi( double radius, double psi ) const
	{
		return stateAt( radius, psi ).potentialPrime;
	}

	double RotatingSource::density( std::size_t index, double radius, double psi ) const
	{
		if ( index >= speciesData.size() )
			throw std::out_of_range( "meq::RotatingSource::density: no such species" );

		return ( *speciesData[ index ].density )( psi )*std::exp( stateAt( radius, psi ).exponent[ index ] );
	}

	void RotatingSource::densityDerivatives( std::size_t index, double radius, double psi,
	                                         double &value, double &dRadius, double &dPsi ) const
	{
		if ( index >= speciesData.size() )
			throw std::out_of_range( "meq::RotatingSource::densityDerivatives: no such species" );

		State const state = stateAt( radius, psi );
		Species const &species = speciesData[ index ];
		double const nRef = ( *species.density )( psi );
		double const e = std::exp( state.exponent[ index ] );
		value = nRef*e;
		dPsi = ( species.density->prime( psi ) + nRef*state.exponentPrime[ index ] )*e;

		/*
		 * d/dR AT FIXED psi. The centrifugal half of A_s is m_s omega^2 h/T_s
		 * with dh/dR = R in every gauge, and the potential's half is
		 * -Z_s y/T_s, so dA_s/dR = ( m_s omega^2 R - Z_s y_R )/T_s with y_R
		 * from (97) held along R:
		 *
		 *     Sum_s Z_s n_s ( m_s omega^2 R - Z_s y_R )/T_s = 0,
		 *
		 * the same denominator Sum_s Z_s^2 n_s/T_s the root find's Newton
		 * step uses. At two species it is the closed form's omega^2 R num/d.
		 * The n_s enter only as ratios, so they are taken with the largest
		 * exponent factored out, which a high Mach number cannot overflow.
		 */
		double const w0 = omegaProfile ? ( *omegaProfile )( psi ) : 0.0;
		double const w = w0*w0;
		double top = -std::numeric_limits<double>::infinity();
		for ( std::size_t s = 0; s < speciesData.size(); ++s )
			top = std::max( top, state.exponent[ s ] );
		double numer = 0.0, denom = 0.0;
		for ( std::size_t s = 0; s < speciesData.size(); ++s )
		{
			Species const &other = speciesData[ s ];
			double const t = ( *other.temperature )( psi );
			double const v = ( *other.density )( psi )*std::exp( state.exponent[ s ] - top );
			numer += other.charge*v*other.mass*w*radius/t;
			denom += other.charge*other.charge*v/t;
		}
		double const yR = denom > 0.0 ? numer/denom : 0.0;
		double const t = ( *species.temperature )( psi );
		dRadius = value*( species.mass*w*radius - species.charge*yR )/t;
	}

	double RotatingSource::pressure( double radius, double psi ) const
	{
		double p = 0.0, pPrime = 0.0, pDoublePrime = 0.0;
		pressureFrom( stateAt( radius, psi ), psi, p, pPrime, pDoublePrime );

		return p;
	}

	double RotatingSource::dPressureDPsi( double radius, double psi ) const
	{
		double p = 0.0, pPrime = 0.0, pDoublePrime = 0.0;
		pressureFrom( stateAt( radius, psi ), psi, p, pPrime, pDoublePrime );

		return pPrime;
	}

	double RotatingSource::f( double radius, double, double psi ) const
	{
		double p = 0.0, pPrime = 0.0, pDoublePrime = 0.0;
		pressureFrom( stateAt( radius, psi ), psi, p, pPrime, pDoublePrime );

		return permeability*radius*radius*pPrime + ( *ggPrimeProfile )( psi );
	}

	double RotatingSource::dFdPsi( double radius, double, double psi ) const
	{
		double p = 0.0, pPrime = 0.0, pDoublePrime = 0.0;
		pressureFrom( stateAt( radius, psi ), psi, p, pPrime, pDoublePrime );

		return permeability*radius*radius*pDoublePrime + ggPrimeProfile->prime( psi );
	}

	std::vector<Species> const & RotatingSource::species() const
	{
		return speciesData;
	}

	RotatingSource::Closure RotatingSource::closure() const
	{
		return closureChoice;
	}

	Profile const * RotatingSource::omega() const
	{
		return omegaProfile.get();
	}

	Profile const & RotatingSource::ggPrime() const
	{
		return *ggPrimeProfile;
	}

	RotatingSource::Gauge RotatingSource::gauge() const
	{
		return gaugeChoice;
	}

	double RotatingSource::referenceRadius() const
	{
		return refRadius;
	}

	double RotatingSource::mu0() const
	{
		return permeability;
	}

	NormalisedRotatingSource::NormalisedRotatingSource( std::vector<Species> species,
		std::shared_ptr<Profile const> omega,
		std::shared_ptr<Profile const> ggPrime,
		double referenceRadius,
		double psiAxis,
		double mu0,
		RotatingSource::Closure closure,
		RotatingSource::Gauge gauge )
		: inner( std::move( species ), std::move( omega ), std::move( ggPrime ), referenceRadius, mu0,
		         closure, gauge ),
		  psiAxisValue( 1.0 )
	{
		// Validation lives in setNormalisation, in one place, exactly as
		// meq::NormalisedMHDSource does it.
		setNormalisation( psiAxis );
	}

	void NormalisedRotatingSource::setNormalisation( double psiAxis,
	                                                 double psiBoundary )
	{
		// The SPAN is what must not vanish, not psi_ax; see
		// meq::NormalisedSource for why that generalisation is FB-3's.
		if ( !std::isfinite( psiAxis ) || !std::isfinite( psiBoundary )
		     || psiAxis == psiBoundary )
			throw std::invalid_argument( "meq::NormalisedRotatingSource::setNormalisation: psi_ax and psi_bnd must be finite and must differ" );

		psiAxisValue = psiAxis;
		psiBoundaryValue = psiBoundary;
	}

	double NormalisedRotatingSource::boundaryNormalisation() const
	{
		return psiBoundaryValue;
	}

	double NormalisedRotatingSource::normalisation() const
	{
		return psiAxisValue;
	}

	double NormalisedRotatingSource::f( double radius, double z, double psi ) const
	{
		// OUTSIDE THE PLASMA THERE IS NO SOURCE, exactly as for
		// meq::NormalisedMHDSource, and the flag reaching here is the whole of
		// what setPlasmaSupport() means. Off unless it was asked for, and then
		// this is one comparison.
		//
		// AND IT MATTERS MORE HERE THAN A MISSING ZERO USUALLY WOULD, because
		// the profiles CLAMP: meq::SplineProfile holds its endpoint value
		// outside the knot range, so an ungated source does not tail off out
		// here -- it sits at its plasma-edge value all the way to the wall,
		// which is a current density filling the vacuum region.
		if ( !insidePlasma( radius, z, psi ) )
			return 0.0;

		// One factor of 1/psi_ax, because the profiles are functions of Psi and
		// F is a psi-derivative of what they build.
		//
		// AND currentScale() MULTIPLIES IT, exactly as it multiplies
		// meq::NormalisedMHDSource::f(). Leaving it off is not a missing
		// convenience: `[source] Type = "rotating"` accepts `PlasmaCurrent`
		// (Config.cpp's rotating key list) and apps/meq.cpp calls
		// setPlasmaCurrent() with no source-type guard, so an unscaled f() gives
		// a current border whose constraint cannot respond to the unknown it is
		// solving for -- while assembleCurrentColumn() and cornerEntry( I, I )
		// both compute their entries as scaledF/lambda and so report a
		// sensitivity that does not exist. The Newton is then driven by a
		// Jacobian describing a different problem from the residual.
		double const span = psiAxisValue - psiBoundaryValue;
		return currentScale()
		       *inner.f( radius, z, ( psi - psiBoundaryValue )/span )/span;
	}

	double NormalisedRotatingSource::dFdPsi( double radius, double z, double psi ) const
	{
		// The derivative of a source that is identically zero out here is zero.
		// NOT the one-sided limit from inside: at the edge itself dF/dpsi picks
		// up a surface term F delta( Psi ), which this interface structurally
		// cannot carry -- and which vanishes exactly when the profiles vanish at
		// the edge, the condition setPlasmaSupport() documents as its
		// precondition. Word for word meq::NormalisedMHDSource's position, and
		// for the same reason.
		if ( !insidePlasma( radius, z, psi ) )
			return 0.0;

		// TWO factors, not one: the chain rule supplies a second whenever another
		// psi-derivative is taken. meq::NormalisedMHDSource carries the same
		// asymmetry, and RotatingSourceTests checks it against a difference.
		// The scale multiplies F, so it multiplies every psi-derivative of F.
		double const span = psiAxisValue - psiBoundaryValue;
		return currentScale()
		       *inner.dFdPsi( radius, z, ( psi - psiBoundaryValue )/span )/( span*span );
	}

	bool NormalisedRotatingSource::normalisationDerivatives( double radius, double z,
	                                                         double psi,
	                                                         double &dFdAxis,
	                                                         double &dFdBoundary ) const
	{
		// Outside the plasma F is identically zero however the normalisation
		// moves, so both derivatives are too -- and this branch is the whole
		// reason an analytic column beats a differenced one rather than merely
		// costing less: a difference perturbs the normalisation, which MOVES THE
		// EDGE, and then straddles it. Word for word
		// meq::NormalisedMHDSource's position.
		if ( !insidePlasma( radius, z, psi ) )
		{
			dFdAxis = 0.0;
			dFdBoundary = 0.0;
			return true;
		}

		/*
		 * THE CLOSURE IS NOT DIFFERENTIATED AGAIN, WHICH IS WHY THIS IS SIX
		 * LINES AND NOT A DERIVATION.
		 *
		 * This class is a pure wrapper: F = S H( Psi )/span with
		 * H( Psi ) := inner.f( R, z, Psi ), and inner.dFdPsi( R, z, Psi ) is
		 * exactly dH/dPsi -- RotatingSource::f() and ::dFdPsi() differ by one
		 * level of pressureFrom()'s output and nothing else. So the
		 * normalisation enters ONLY through the argument Psi and the overall
		 * 1/span, and whatever the species root find cost to differentiate once
		 * has already been paid inside inner.dFdPsi(). The two-species closed
		 * form and the n-species root find are both invisible here.
		 *
		 * That leaves meq::NormalisedMHDSource's algebra verbatim under
		 * g -> H, g' -> H':
		 *
		 *   dPsi/dpsi_ax  = -Psi/span,        dspan/dpsi_ax  = +1
		 *   dPsi/dpsi_bnd = ( Psi - 1 )/span, dspan/dpsi_bnd = -1
		 */
		double const span = psiAxisValue - psiBoundaryValue;
		double const psiN = ( psi - psiBoundaryValue )/span;
		double const h = inner.f( radius, z, psiN );
		double const hPrime = inner.dFdPsi( radius, z, psiN );

		dFdAxis = -currentScale()*( hPrime*psiN + h )/( span*span );
		dFdBoundary = currentScale()*( hPrime*( psiN - 1.0 ) + h )/( span*span );
		return true;
	}

	double NormalisedRotatingSource::potential( double radius, double psi ) const
	{
		return inner.potential( radius, ( psi - psiBoundaryValue )/( psiAxisValue - psiBoundaryValue ) );
	}

	double NormalisedRotatingSource::density( std::size_t index, double radius, double psi ) const
	{
		return inner.density( index, radius,
		                      ( psi - psiBoundaryValue )/( psiAxisValue - psiBoundaryValue ) );
	}

	double NormalisedRotatingSource::pressure( double radius, double psi ) const
	{
		return inner.pressure( radius, ( psi - psiBoundaryValue )/( psiAxisValue - psiBoundaryValue ) );
	}

	RotatingSource const & NormalisedRotatingSource::unnormalised() const
	{
		return inner;
	}

	RotatingSource & NormalisedRotatingSource::unnormalised()
	{
		return inner;
	}

	bool NormalisedRotatingSource::readsMidplane() const
	{
		return inner.gauge() == RotatingSource::Gauge::Midplane;
	}

	void NormalisedRotatingSource::setMidplaneField( MidplaneField const &field,
	                                                 double rMin, double rMax )
	{
		if ( !readsMidplane() )
			throw std::logic_error( "meq::NormalisedRotatingSource::setMidplaneField: this source is not in Gauge::Midplane and reads no midplane" );

		// The field is in PHYSICAL psi; the inner source's argument is Psi, so
		// the table is taken in Psi at the normalisation in force -- see the
		// declaration for why that keeps this class a pure wrapper.
		double const base = psiBoundaryValue;
		double const span = psiAxisValue - psiBoundaryValue;
		MidplaneField const normalised = [ & ]( double radius, double &x, double &dXdR )
		{
			double psi = 0.0, dPsiDR = 0.0;
			if ( !field( radius, psi, dPsiDR ) )
				return false;
			x = ( psi - base )/span;
			dXdR = dPsiDR/span;
			return true;
		};

		inner.setMidplaneCrossing( tabulateMidplaneCrossing( normalised, rMin, rMax, 0.0, 1.0 ) );
	}

	bool NormalisedRotatingSource::readsVacuumMidplane() const
	{
		return inner.gauge() == RotatingSource::Gauge::FixedCurve && !inner.midplaneCrossing();
	}

	void NormalisedRotatingSource::setVacuumMidplaneField( MidplaneField const &field,
	                                                       double rMin, double rMax,
	                                                       double psiAxis, double psiBoundary )
	{
		if ( inner.gauge() != RotatingSource::Gauge::FixedCurve )
			throw std::logic_error( "meq::NormalisedRotatingSource::setVacuumMidplaneField: this source is not in Gauge::FixedCurve and has no fixed curve to be given" );
		if ( inner.midplaneCrossing() )
			throw std::logic_error( "meq::NormalisedRotatingSource::setVacuumMidplaneField: the fixed curve is already in, and a fixed curve is never replaced -- replacing it would change what every Density means" );
		if ( !std::isfinite( psiAxis ) || !std::isfinite( psiBoundary ) || psiAxis == psiBoundary )
			throw std::invalid_argument( "meq::NormalisedRotatingSource::setVacuumMidplaneField: the vacuum's labels must be finite and must differ" );

		// In Psi, against the VACUUM's labels and not the normalisation in
		// force: the solve's starting labels are a guess, and the curve says
		// what the profiles mean.
		double const span = psiAxis - psiBoundary;
		MidplaneField const normalised = [ & ]( double radius, double &x, double &dXdR )
		{
			double psi = 0.0, dPsiDR = 0.0;
			if ( !field( radius, psi, dPsiDR ) )
				return false;
			x = ( psi - psiBoundary )/span;
			dXdR = dPsiDR/span;
			return true;
		};
		std::shared_ptr<SplineProfile const> const curve =
			tabulateMidplaneCrossing( normalised, rMin, rMax, 0.0, 1.0 );

		/*
		 * THE BOX'S PRECONDITION, CHECKED BECAUSE IT NOW CAN BE. F must vanish
		 * at both labels -- a jump across a label's contour caps the order,
		 * and the Jacobian misses its edge term -- and on a fixed curve
		 * whether it does is a property of the profiles and R_g alone, so a
		 * set that fails it is a configuration error and is refused as one.
		 * See meq::labelTolerance for what the bound guards, MEASURED --
		 * [M-VACGAUGE](MEASUREMENTS.md#m-vacgauge).
		 */
		RotatingSource::LabelCheck const check = inner.labelCheck( curve );
		if ( !( check.relative <= labelTolerance ) )
		{
			char buffer[ 768 ];
			std::snprintf( buffer, sizeof buffer,
				"F at the %s label, Psi = %g, is %.3e of its peak at R = %.4g m on the "
				"vacuum's midplane crossing, above the %.0e the box allows: the support "
				"needs F to vanish at both labels, since a jump across a label's contour "
				"caps the convergence order and the Newton Jacobian misses its edge "
				"term. On this gauge's curve R_g a Density flat at a label still "
				"carries -n m omega^2 R_g R_g'/T in its psi-derivative wherever omega "
				"does not vanish there, so let the densities or the rotation fall away "
				"at the labels, or give each Density the slope that cancels it",
				check.label == 1.0 ? "inner" : "outer", check.label, check.relative,
				check.radius, labelTolerance );
			throw std::invalid_argument( buffer );
		}
		inner.setMidplaneCrossing( curve );
	}

	std::shared_ptr<SplineProfile const> densityOnCurve( RotatingSource const &source,
		std::size_t index, SplineProfile const &curve )
	{
		std::vector<Knot> knots;
		knots.reserve( curve.knots().size() );
		for ( Knot const &k : curve.knots() )
		{
			double value = 0.0, dRadius = 0.0, dPsi = 0.0;
			source.densityDerivatives( index, k.value, k.psi, value, dRadius, dPsi );
			knots.push_back( Knot{ k.psi, value, dPsi + dRadius*k.derivative } );
		}
		return std::make_shared<SplineProfile const>( std::move( knots ) );
	}

	void NormalisedRotatingSource::setMidplaneCrossing( std::shared_ptr<Profile const> crossing )
	{
		inner.setMidplaneCrossing( std::move( crossing ) );
	}

	bool NormalisedRotatingSource::crossingDerivatives( double radius, double z, double psi,
	                                                    double &dFdW, double &dFdWPrime ) const
	{
		if ( !readsMidplane() )
			throw std::logic_error( "meq::NormalisedRotatingSource::crossingDerivatives: this source is in Gauge::FixedRadius, whose F reads no crossing" );
		dFdW = 0.0;
		dFdWPrime = 0.0;
		if ( !insidePlasma( radius, z, psi ) )
			return true;

		// F = S H( Psi )/span with W a function of Psi, so the span enters
		// once, as it does in f(): W and W' are the table's at Psi and are
		// held while they move.
		double const span = psiAxisValue - psiBoundaryValue;
		inner.crossingDerivatives( radius, ( psi - psiBoundaryValue )/span, dFdW, dFdWPrime );
		dFdW *= currentScale()/span;
		dFdWPrime *= currentScale()/span;
		return true;
	}

}
