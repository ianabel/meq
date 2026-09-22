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

	double MHDSource::f( double radius, double, double psi ) const
	{
		return permeability*radius*radius*( *pPrimeProfile )( psi ) + ( *ggPrimeProfile )( psi );
	}

	double MHDSource::dFdPsi( double radius, double, double psi ) const
	{
		return permeability*radius*radius*pPrimeProfile->prime( psi ) + ggPrimeProfile->prime( psi );
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

	double NormalisedMHDSource::f( double radius, double, double psi ) const
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
		return currentScale()
		       *( permeability*radius*radius*( *pPrimeProfile )( psiN ) + ( *ggPrimeProfile )( psiN ) )
		       /span;
	}

	double NormalisedMHDSource::dFdPsi( double radius, double, double psi ) const
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
		return currentScale()
		       *( permeability*radius*radius*pPrimeProfile->prime( psiN ) + ggPrimeProfile->prime( psiN ) )
		       /( span*span );
	}

	bool NormalisedMHDSource::normalisationDerivatives( double radius, double /*z*/,
	                                                    double psi,
	                                                    double &dFdAxis,
	                                                    double &dFdBoundary ) const
	{
		// Outside the plasma F is identically zero however the normalisation
		// moves, so both derivatives are too. This is the branch that makes an
		// analytic column better than a differenced one rather than merely
		// cheaper: a difference perturbs the normalisation, which MOVES THE
		// EDGE, and then straddles it.
		if ( !insidePlasma( psi ) )
		{
			dFdAxis = 0.0;
			dFdBoundary = 0.0;
			return true;
		}

		double const span = psiAxisValue - psiBoundaryValue;
		double const psiN = ( psi - psiBoundaryValue )/span;

		// g( Psi ) is exactly what f() assembles before dividing by the span,
		// and g'( Psi ) is one Profile::prime() of each stored profile -- the
		// SECOND derivative of p and of gg, which is the level meq::Profile
		// already carries because the rotating source needed it.
		double const g = permeability*radius*radius*( *pPrimeProfile )( psiN )
		                 + ( *ggPrimeProfile )( psiN );
		double const gPrime = permeability*radius*radius*pPrimeProfile->prime( psiN )
		                      + ggPrimeProfile->prime( psiN );

		// F = g( Psi )/span with dPsi/dpsi_ax = -Psi/span, dspan/dpsi_ax = +1,
		// dPsi/dpsi_bnd = ( Psi - 1 )/span, dspan/dpsi_bnd = -1.
		// The scale multiplies F, so it multiplies both derivatives too.
		dFdAxis = -currentScale()*( gPrime*psiN + g )/( span*span );
		dFdBoundary = currentScale()*( gPrime*( psiN - 1.0 ) + g )/( span*span );
		return true;
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
	
	void NormalisedMHDSource::setGGPrime( std::shared_ptr<Profile const> ggPrime )
	{
		if ( !ggPrime )
			throw std::invalid_argument(
				"meq::NormalisedMHDSource::setGGPrime: no profile. A null g dg/"
				"dPsi is not a vacuum -- it is a source that cannot be "
				"evaluated" );
		ggPrimeProfile = std::move( ggPrime );
	}

	SolovievSource::SolovievSource( double a )
		: aValue( a )
	{
		if ( !std::isfinite( a ) )
			throw std::invalid_argument( "meq::SolovievSource: A must be finite" );
	}

	double SolovievSource::f( double radius, double, double ) const
	{
		return -( ( 1.0 - aValue )*radius*radius + aValue );
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
