#include "RotatingSource.hpp"

#include <cmath>
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
		double mu0 )
		: speciesData( std::move( species ) ), omegaProfile( std::move( omega ) ),
		  ggPrimeProfile( std::move( ggPrime ) ), rRef( referenceRadius ), permeability( mu0 )
	{
		// Two species only. With exactly two, (97) is linear in phi_0 after taking
		// logarithms and there is no root find; three or more needs the
		// safeguarded scalar Newton of FLOW-PLAN.md's FL-6, which is not written.
		// Refusing is the house behaviour: approximating here would give a
		// converged answer to a problem nobody posed.
		if ( speciesData.size() != 2 )
			throw std::invalid_argument( "meq::RotatingSource: exactly two species are supported; three or more needs the general phi_0 root find, which is not implemented" );

		requireProfile( ggPrimeProfile, "the g dg/dpsi profile" );

		if ( !std::isfinite( rRef ) || rRef <= 0.0 )
			throw std::invalid_argument( "meq::RotatingSource: the reference radius must be finite and positive" );
		if ( !std::isfinite( permeability ) )
			throw std::invalid_argument( "meq::RotatingSource: mu0 must be finite" );

		for ( Species const & s : speciesData )
		{
			requireProfile( s.temperature, "a temperature profile" );
			requireProfile( s.density, "a density profile" );

			if ( !std::isfinite( s.mass ) || s.mass <= 0.0 )
				throw std::invalid_argument( "meq::RotatingSource: every species mass must be finite and positive" );
			if ( !std::isfinite( s.charge ) || s.charge == 0.0 )
				throw std::invalid_argument( "meq::RotatingSource: every species charge must be finite and non-zero" );
		}

		// Opposite signs, which is what makes the denominator Z_1 T_2 - Z_2 T_1
		// positive for positive temperatures and hence makes the closed form below
		// well defined. It is also what gives (97) a root at all: the left hand
		// side has to change sign.
		if ( speciesData[ 0 ].charge*speciesData[ 1 ].charge > 0.0 )
			throw std::invalid_argument( "meq::RotatingSource: the two species must carry opposite charge, or quasineutrality has no solution" );

		// The closed form is derived from Z_1 n_10 = -Z_2 n_20 and is simply wrong
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

	void RotatingSource::referencePressure( double psi, double & p0, double & p0Prime, double & p0DoublePrime ) const
	{
		p0 = 0.0;
		p0Prime = 0.0;
		p0DoublePrime = 0.0;

		for ( Species const & s : speciesData )
		{
			double const n = ( *s.density )( psi );
			double const nPrime = s.density->prime( psi );
			double const nDoublePrime = s.density->doublePrime( psi );

			double const t = ( *s.temperature )( psi );
			double const tPrime = s.temperature->prime( psi );
			double const tDoublePrime = s.temperature->doublePrime( psi );

			p0 += n*t;
			p0Prime += nPrime*t + n*tPrime;
			p0DoublePrime += nDoublePrime*t + 2.0*nPrime*tPrime + n*tDoublePrime;
		}
	}

	void RotatingSource::exponentCoefficient( double psi, double & c, double & cPrime, double & cDoublePrime ) const
	{
		// Both species carry the SAME exponent, which is what makes Sum_s Z_s n_s
		// vanish at every r once it vanishes at rRef. Solving (97) for two species
		// gives that shared coefficient as
		//
		//     C( psi ) = omega^2 ( Z_1 m_2 - Z_2 m_1 )/( Z_1 T_2 - Z_2 T_1 )
		//
		// which for ions and electrons ( Z_2 = -1 ) is the familiar
		// omega^2 ( m_i + Z_i m_e )/( T_i + Z_i T_e ). The electron mass is kept
		// rather than dropped: it costs one term and removes a question.
		if ( !omegaProfile )
		{
			c = 0.0;
			cPrime = 0.0;
			cDoublePrime = 0.0;
			return;
		}

		Species const & one = speciesData[ 0 ];
		Species const & two = speciesData[ 1 ];

		double const k = one.charge*two.mass - two.charge*one.mass;

		double const t1 = ( *one.temperature )( psi );
		double const t2 = ( *two.temperature )( psi );
		double const d = one.charge*t2 - two.charge*t1;

		if ( d == 0.0 )
			throw std::invalid_argument( "meq::RotatingSource: Z_1 T_2 - Z_2 T_1 vanishes, so the quasineutrality closure is singular; with opposite charges this needs a temperature to have gone non-positive" );

		double const dPrime = one.charge*two.temperature->prime( psi ) - two.charge*one.temperature->prime( psi );
		double const dDoublePrime = one.charge*two.temperature->doublePrime( psi ) - two.charge*one.temperature->doublePrime( psi );

		// W = omega^2, so that the chain rule below is in one variable rather than
		// two. omega itself never appears alone in the answer.
		double const w0 = ( *omegaProfile )( psi );
		double const wPrime0 = omegaProfile->prime( psi );
		double const wDoublePrime0 = omegaProfile->doublePrime( psi );

		double const w = w0*w0;
		double const wPrime = 2.0*w0*wPrime0;
		double const wDoublePrime = 2.0*( wPrime0*wPrime0 + w0*wDoublePrime0 );

		c = k*w/d;
		cPrime = k*( wPrime/d - w*dPrime/( d*d ) );
		cDoublePrime = k*( wDoublePrime/d - 2.0*wPrime*dPrime/( d*d ) - w*dDoublePrime/( d*d )
			+ 2.0*w*dPrime*dPrime/( d*d*d ) );
	}

	double RotatingSource::densityExponent( double r, double psi ) const
	{
		double c = 0.0, cPrime = 0.0, cDoublePrime = 0.0;
		exponentCoefficient( psi, c, cPrime, cDoublePrime );

		return c*( r*r - rRef*rRef )/2.0;
	}

	double RotatingSource::potential( double r, double psi ) const
	{
		// e phi_0 = omega^2 ( r^2 - rRef^2 )( m_1 T_2 - m_2 T_1 )/( 2 ( Z_1 T_2 - Z_2 T_1 ) ),
		// which is (97) solved for the potential in the gauge phi_0( rRef ) = 0.
		// The numerator is a mass-weighted temperature difference: with equal
		// mass-to-temperature ratios there is nothing for the field to separate
		// and phi_0 vanishes identically.
		if ( !omegaProfile )
			return 0.0;

		Species const & one = speciesData[ 0 ];
		Species const & two = speciesData[ 1 ];

		double const t1 = ( *one.temperature )( psi );
		double const t2 = ( *two.temperature )( psi );
		double const d = one.charge*t2 - two.charge*t1;

		if ( d == 0.0 )
			throw std::invalid_argument( "meq::RotatingSource::potential: Z_1 T_2 - Z_2 T_1 vanishes" );

		double const w0 = ( *omegaProfile )( psi );

		return w0*w0*( r*r - rRef*rRef )*( one.mass*t2 - two.mass*t1 )/( 2.0*d );
	}

	double RotatingSource::density( std::size_t index, double r, double psi ) const
	{
		if ( index >= speciesData.size() )
			throw std::out_of_range( "meq::RotatingSource::density: no such species" );

		return ( *speciesData[ index ].density )( psi )*std::exp( densityExponent( r, psi ) );
	}

	double RotatingSource::pressure( double r, double psi ) const
	{
		double p0 = 0.0, p0Prime = 0.0, p0DoublePrime = 0.0;
		referencePressure( psi, p0, p0Prime, p0DoublePrime );

		return p0*std::exp( densityExponent( r, psi ) );
	}

	double RotatingSource::dPressureDPsi( double r, double psi ) const
	{
		double p0 = 0.0, p0Prime = 0.0, p0DoublePrime = 0.0;
		referencePressure( psi, p0, p0Prime, p0DoublePrime );

		double c = 0.0, cPrime = 0.0, cDoublePrime = 0.0;
		exponentCoefficient( psi, c, cPrime, cDoublePrime );

		double const halfDelta = ( r*r - rRef*rRef )/2.0;
		double const q = c*halfDelta;
		double const qPrime = cPrime*halfDelta;

		return std::exp( q )*( p0Prime + p0*qPrime );
	}

	double RotatingSource::f( double r, double, double psi ) const
	{
		return permeability*r*r*dPressureDPsi( r, psi ) + ( *ggPrimeProfile )( psi );
	}

	double RotatingSource::dFdPsi( double r, double, double psi ) const
	{
		double p0 = 0.0, p0Prime = 0.0, p0DoublePrime = 0.0;
		referencePressure( psi, p0, p0Prime, p0DoublePrime );

		double c = 0.0, cPrime = 0.0, cDoublePrime = 0.0;
		exponentCoefficient( psi, c, cPrime, cDoublePrime );

		double const halfDelta = ( r*r - rRef*rRef )/2.0;
		double const q = c*halfDelta;
		double const qPrime = cPrime*halfDelta;
		double const qDoublePrime = cDoublePrime*halfDelta;

		// d2p/dpsi2 of p = P0 exp( q ), both P0 and q being functions of psi. This
		// is where the second derivative of every input profile is spent, and why
		// meq::Profile had to grow doublePrime().
		double const pDoublePrime = std::exp( q )
			*( p0DoublePrime + 2.0*p0Prime*qPrime + p0*( qPrime*qPrime + qDoublePrime ) );

		return permeability*r*r*pDoublePrime + ggPrimeProfile->prime( psi );
	}

	std::vector<Species> const & RotatingSource::species() const
	{
		return speciesData;
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

}
