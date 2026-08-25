#include "Source.hpp"

#include <cmath>
#include <stdexcept>
#include <utility>

namespace meq
{

	MHDSource::MHDSource( std::shared_ptr<Profile const> pPrime, std::shared_ptr<Profile const> ggPrime, double mu0 )
		: pPrimeProfile( std::move( pPrime ) ), ggPrimeProfile( std::move( ggPrime ) ), permeability( mu0 )
	{
		if ( !pPrimeProfile )
			throw std::invalid_argument( "meq::MHDSource: the dp/dpsi profile must not be null" );
		if ( !ggPrimeProfile )
			throw std::invalid_argument( "meq::MHDSource: the g dg/dpsi profile must not be null" );
		if ( !std::isfinite( permeability ) )
			throw std::invalid_argument( "meq::MHDSource: mu0 must be finite" );
	}

	double MHDSource::f( double r, double, double psi ) const
	{
		return permeability*r*r*( *pPrimeProfile )( psi ) + ( *ggPrimeProfile )( psi );
	}

	double MHDSource::dFdPsi( double r, double, double psi ) const
	{
		return permeability*r*r*pPrimeProfile->prime( psi ) + ggPrimeProfile->prime( psi );
	}

	Profile const & MHDSource::pPrime() const
	{
		return *pPrimeProfile;
	}

	Profile const & MHDSource::ggPrime() const
	{
		return *ggPrimeProfile;
	}

	double MHDSource::mu0() const
	{
		return permeability;
	}

	SolovievSource::SolovievSource( double a )
		: aValue( a )
	{
		if ( !std::isfinite( a ) )
			throw std::invalid_argument( "meq::SolovievSource: A must be finite" );
	}

	double SolovievSource::f( double r, double, double ) const
	{
		return -( ( 1.0 - aValue )*r*r + aValue );
	}

	double SolovievSource::dFdPsi( double, double, double ) const
	{
		return 0.0;
	}

	double SolovievSource::a() const
	{
		return aValue;
	}

	double SolovievSource::c() const
	{
		return 1.0 - aValue;
	}

}
