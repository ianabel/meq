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

	NormalisedMHDSource::NormalisedMHDSource( std::shared_ptr<Profile const> pPrime,
	                                          std::shared_ptr<Profile const> ggPrime,
	                                          double psiAxis, double mu0 )
		: pPrimeProfile( std::move( pPrime ) ), ggPrimeProfile( std::move( ggPrime ) ),
		  psiAxisValue( 1.0 ), permeability( mu0 )
	{
		if ( !pPrimeProfile )
			throw std::invalid_argument( "meq::NormalisedMHDSource: the dp/dPsi profile must not be null" );
		if ( !ggPrimeProfile )
			throw std::invalid_argument( "meq::NormalisedMHDSource: the g dg/dPsi profile must not be null" );
		if ( !std::isfinite( permeability ) )
			throw std::invalid_argument( "meq::NormalisedMHDSource: mu0 must be finite" );
		setNormalisation( psiAxis );
	}

	void NormalisedMHDSource::setNormalisation( double psiAxis )
	{
		// Loudly rather than by returning infinities. A solver whose iterate has
		// reached psi_ax = 0 has left the branch, and the degenerate fixed point
		// where psi and psi_ax shrink together is exactly the failure this class
		// exists to make impossible -- so it has to be a throw and not a floor.
		if ( !std::isfinite( psiAxis ) || psiAxis == 0.0 )
			throw std::invalid_argument( "meq::NormalisedMHDSource::setNormalisation: psi_ax must be finite and non-zero" );
		psiAxisValue = psiAxis;
	}

	double NormalisedMHDSource::normalisation() const
	{
		return psiAxisValue;
	}

	double NormalisedMHDSource::f( double r, double, double psi ) const
	{
		double const psiN = psi/psiAxisValue;
		return ( permeability*r*r*( *pPrimeProfile )( psiN ) + ( *ggPrimeProfile )( psiN ) )
		       /psiAxisValue;
	}

	double NormalisedMHDSource::dFdPsi( double r, double, double psi ) const
	{
		// Two factors of 1/psi_ax, not one: the profiles are differentiated with
		// respect to Psi and the argument carries a further 1/psi_ax. Dropping the
		// second is the classic error here, and it does not move the converged
		// answer -- only the convergence to it.
		double const psiN = psi/psiAxisValue;
		return ( permeability*r*r*pPrimeProfile->prime( psiN ) + ggPrimeProfile->prime( psiN ) )
		       /( psiAxisValue*psiAxisValue );
	}

	Profile const & NormalisedMHDSource::pPrime() const
	{
		return *pPrimeProfile;
	}

	Profile const & NormalisedMHDSource::ggPrime() const
	{
		return *ggPrimeProfile;
	}

	double NormalisedMHDSource::mu0() const
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
