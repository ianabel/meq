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

	void NormalisedMHDSource::setNormalisation( double psiAxis, double psiBoundary )
	{
		// Loudly rather than by returning infinities. A solver whose iterate has
		// collapsed the SPAN has left the branch, and the degenerate fixed point
		// where psi and the span shrink together is exactly the failure this
		// class exists to make impossible -- so it has to be a throw and not a
		// floor.
		//
		// THE TEST IS ON THE SPAN AND NOT ON psi_ax, which is the generalisation
		// FB-3 needs: with psi_bnd = 0 the two coincide, and with psi_bnd set
		// they do not. psi_ax = 0 is perfectly admissible once the boundary flux
		// is non-zero -- it is psi_ax = psi_bnd that is not.
		if ( !std::isfinite( psiAxis ) || !std::isfinite( psiBoundary )
		     || psiAxis == psiBoundary )
			throw std::invalid_argument( "meq::NormalisedMHDSource::setNormalisation: psi_ax and psi_bnd must be finite and must differ" );
		psiAxisValue = psiAxis;
		psiBoundaryValue = psiBoundary;
	}

	double NormalisedMHDSource::boundaryNormalisation() const
	{
		return psiBoundaryValue;
	}

	double NormalisedMHDSource::normalisation() const
	{
		return psiAxisValue;
	}

	double NormalisedMHDSource::f( double r, double, double psi ) const
	{
		// Psi = ( psi - psi_bnd )/span, and the profiles are differentiated with
		// respect to Psi, so F carries one factor of 1/span. With psi_bnd = 0
		// the span IS psi_ax and this is the expression it always was.
		// OUTSIDE THE PLASMA THERE IS NO SOURCE, and that is the whole of the
		// moving support: nothing here knows where the boundary is, and the
		// boundary is not an input -- it is wherever psi currently puts it. Off
		// unless setPlasmaSupport() asked for it, and then this is one
		// comparison. See NormalisedSource::setPlasmaSupport.
		if ( !insidePlasma( psi ) )
			return 0.0;

		double const span = psiAxisValue - psiBoundaryValue;
		double const psiN = ( psi - psiBoundaryValue )/span;
		return ( permeability*r*r*( *pPrimeProfile )( psiN ) + ( *ggPrimeProfile )( psiN ) )
		       /span;
	}

	double NormalisedMHDSource::dFdPsi( double r, double, double psi ) const
	{
		// Two factors of 1/span, not one: the profiles are differentiated with
		// respect to Psi and the argument carries a further 1/span. Dropping the
		// second is the classic error here, and it does not move the converged
		// answer -- only the convergence to it.
		// The derivative of a source that is identically zero out here is zero.
		// NOT the one-sided limit from inside: at the edge itself dF/dpsi picks
		// up a surface term F delta( Psi ), which this interface structurally
		// cannot carry -- and which VANISHES exactly when the profiles vanish
		// at the edge, which setPlasmaSupport() documents as its precondition
		// and CLAUDE.md measures the cost of violating.
		if ( !insidePlasma( psi ) )
			return 0.0;

		double const span = psiAxisValue - psiBoundaryValue;
		double const psiN = ( psi - psiBoundaryValue )/span;
		return ( permeability*r*r*pPrimeProfile->prime( psiN ) + ggPrimeProfile->prime( psiN ) )
		       /( span*span );
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
