#include "RotatingSource.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

namespace
{

	// Charge neutrality on r = rRef is checked by sampling, because the densities
	// are profiles and the constraint has to hold at every psi. meq::Profile
	// documents itself on [ 0, 1 ], so that is where the samples go.
	int const neutralitySamples = 21;
	double const neutralityTolerance = 1.0e-10;

	// The root find for (97). The residual is compared against the sum of the
	// magnitudes it cancels, so the test is relative and does not care what units
	// the densities are in.
	int const potentialIterations = 100;
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

	RotatingSource::RotatingSource( std::vector<Species> species,
		std::shared_ptr<Profile const> omega,
		std::shared_ptr<Profile const> ggPrime,
		double referenceRadius,
		double mu0,
		Closure closure )
		: speciesData( std::move( species ) ), omegaProfile( std::move( omega ) ),
		  ggPrimeProfile( std::move( ggPrime ) ), rRef( referenceRadius ), permeability( mu0 ),
		  closureChoice( closure )
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

		if ( !std::isfinite( rRef ) || rRef <= 0.0 )
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

	RotatingSource::State RotatingSource::closedFormState( double r, double psi ) const
	{
		// Two species. Taking logarithms of (97) makes it linear in phi_0, and
		// the two exponents come out EQUAL:
		//
		//     C( psi ) = omega^2 ( Z_1 m_2 - Z_2 m_1 )/( Z_1 T_2 - Z_2 T_1 )
		//     e phi_0  = omega^2 ( r^2 - rRef^2 )( m_1 T_2 - m_2 T_1 )/2( Z_1 T_2 - Z_2 T_1 )
		//
		// which for ions and electrons ( Z_2 = -1 ) is the familiar
		// omega^2 ( m_i + Z_i m_e )/( T_i + Z_i T_e ). The electron mass is kept
		// rather than dropped: it costs one term and removes a question. That the
		// exponents are equal is what makes Sum_s Z_s n_s vanish at every r once
		// it vanishes at rRef.
		State state{};

		Species const & one = speciesData[ 0 ];
		Species const & two = speciesData[ 1 ];

		double const halfDelta = ( r*r - rRef*rRef )/2.0;

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

		for ( std::size_t s = 0; s < 2; ++s )
		{
			state.exponent[ s ] = c*halfDelta;
			state.exponentPrime[ s ] = cPrime*halfDelta;
			state.exponentDoublePrime[ s ] = cDoublePrime*halfDelta;
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

		state.potential = w*halfDelta*ratio;
		state.potentialPrime = halfDelta*( wPrime*ratio + w*ratioPrime );
		state.potentialDoublePrime = halfDelta*( wDoublePrime*ratio + 2.0*wPrime*ratioPrime + w*ratioDoublePrime );

		return state;
	}

	RotatingSource::State RotatingSource::rootFindState( double r, double psi ) const
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
		double const halfDelta = ( r*r - rRef*rRef )/2.0;

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
				value[ s ] = speciesData[ s ].mass*w*halfDelta/t[ s ] - speciesData[ s ].charge*y/t[ s ];
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

			// Newton on G, which is y + numer/denom because G' = -denom. Bisect
			// whenever that leaves the bracket, which is what makes it safeguarded
			// rather than merely fast.
			double next = y + numer/denom;
			if ( !( next > lo && next < hi ) )
				next = 0.5*( lo + hi );

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

			// B_s = dA_s/dpsi at FIXED y.
			b[ s ] = m*halfDelta*( wPrime/ti - w*tp/( ti*ti ) ) + z*y*tp/( ti*ti );

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

			bPrime[ s ] = m*halfDelta*( wDoublePrime/ti - 2.0*wPrime*tp/( ti*ti )
					- w*tpp/( ti*ti ) + 2.0*w*tp*tp/( ti*ti*ti ) )
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

			state.exponent[ s ] = speciesData[ s ].mass*w*halfDelta/ti - z*y/ti;
			state.exponentPrime[ s ] = l[ s ];
			state.exponentDoublePrime[ s ] = bPrime[ s ] - z*yDoublePrime/ti + z*yPrime*tp/( ti*ti );
		}

		return state;
	}

	RotatingSource::State RotatingSource::stateAt( double r, double psi ) const
	{
		return closureChoice == Closure::ClosedForm ? closedFormState( r, psi ) : rootFindState( r, psi );
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

	double RotatingSource::densityExponent( std::size_t index, double r, double psi ) const
	{
		if ( index >= speciesData.size() )
			throw std::out_of_range( "meq::RotatingSource::densityExponent: no such species" );

		return stateAt( r, psi ).exponent[ index ];
	}

	double RotatingSource::potential( double r, double psi ) const
	{
		return stateAt( r, psi ).potential;
	}

	double RotatingSource::dPotentialDPsi( double r, double psi ) const
	{
		return stateAt( r, psi ).potentialPrime;
	}

	double RotatingSource::density( std::size_t index, double r, double psi ) const
	{
		if ( index >= speciesData.size() )
			throw std::out_of_range( "meq::RotatingSource::density: no such species" );

		return ( *speciesData[ index ].density )( psi )*std::exp( stateAt( r, psi ).exponent[ index ] );
	}

	double RotatingSource::pressure( double r, double psi ) const
	{
		double p = 0.0, pPrime = 0.0, pDoublePrime = 0.0;
		pressureFrom( stateAt( r, psi ), psi, p, pPrime, pDoublePrime );

		return p;
	}

	double RotatingSource::dPressureDPsi( double r, double psi ) const
	{
		double p = 0.0, pPrime = 0.0, pDoublePrime = 0.0;
		pressureFrom( stateAt( r, psi ), psi, p, pPrime, pDoublePrime );

		return pPrime;
	}

	double RotatingSource::f( double r, double, double psi ) const
	{
		double p = 0.0, pPrime = 0.0, pDoublePrime = 0.0;
		pressureFrom( stateAt( r, psi ), psi, p, pPrime, pDoublePrime );

		return permeability*r*r*pPrime + ( *ggPrimeProfile )( psi );
	}

	double RotatingSource::dFdPsi( double r, double, double psi ) const
	{
		double p = 0.0, pPrime = 0.0, pDoublePrime = 0.0;
		pressureFrom( stateAt( r, psi ), psi, p, pPrime, pDoublePrime );

		return permeability*r*r*pDoublePrime + ggPrimeProfile->prime( psi );
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

	double RotatingSource::referenceRadius() const
	{
		return rRef;
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
		RotatingSource::Closure closure )
		: inner( std::move( species ), std::move( omega ), std::move( ggPrime ), referenceRadius, mu0, closure ),
		  psiAxisValue( 1.0 )
	{
		// Validation lives in setNormalisation, in one place, exactly as
		// meq::NormalisedMHDSource does it.
		setNormalisation( psiAxis );
	}

	void NormalisedRotatingSource::setNormalisation( double psiAxis )
	{
		if ( !std::isfinite( psiAxis ) || psiAxis == 0.0 )
			throw std::invalid_argument( "meq::NormalisedRotatingSource::setNormalisation: psi_ax must be finite and non-zero" );

		psiAxisValue = psiAxis;
	}

	double NormalisedRotatingSource::normalisation() const
	{
		return psiAxisValue;
	}

	double NormalisedRotatingSource::f( double r, double z, double psi ) const
	{
		// One factor of 1/psi_ax, because the profiles are functions of Psi and
		// F is a psi-derivative of what they build.
		return inner.f( r, z, psi/psiAxisValue )/psiAxisValue;
	}

	double NormalisedRotatingSource::dFdPsi( double r, double z, double psi ) const
	{
		// TWO factors, not one: the chain rule supplies a second whenever another
		// psi-derivative is taken. meq::NormalisedMHDSource carries the same
		// asymmetry, and RotatingSourceTests checks it against a difference.
		return inner.dFdPsi( r, z, psi/psiAxisValue )/( psiAxisValue*psiAxisValue );
	}

	double NormalisedRotatingSource::potential( double r, double psi ) const
	{
		return inner.potential( r, psi/psiAxisValue );
	}

	double NormalisedRotatingSource::density( std::size_t index, double r, double psi ) const
	{
		return inner.density( index, r, psi/psiAxisValue );
	}

	double NormalisedRotatingSource::pressure( double r, double psi ) const
	{
		return inner.pressure( r, psi/psiAxisValue );
	}

	RotatingSource const & NormalisedRotatingSource::unnormalised() const
	{
		return inner;
	}

}
