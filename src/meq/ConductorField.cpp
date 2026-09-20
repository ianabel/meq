#include "ConductorField.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <sstream>
#include <stdexcept>

namespace meq
{
	namespace
	{
		// A local copy rather than a shared one: Coils.cpp's requireFinite() is
		// in an anonymous namespace there, which is where a helper this small
		// belongs. Exporting it to share four lines would put a name in the
		// public header that no caller wants.
		void requireFinite( double value, char const *what, char const *where )
		{
			if ( !std::isfinite( value ) )
			{
				std::ostringstream message;
				message << where << ": " << what << " must be finite, but is "
				        << value;
				throw std::invalid_argument( message.str() );
			}
		}
	}

	ConductorField::ConductorField( double mu0In )
		: mu0Value( mu0In ),
		  toleranceValue( defaultCoincidenceTolerance ),
		  filamentList(),
		  // ONE mu0 FOR BOTH KINDS. meq::CoilSet carries its own, and letting
		  // the two disagree would make psi_c a sum of fields in different
		  // unit systems -- which converges, at the full rate, to a machine
		  // nobody described. The CoilSet's constructor refuses a non-positive
		  // mu0 as this one does, so the check below is not duplicated here.
		  coilList( mu0In )
	{
		requireFinite( mu0In, "mu0", "meq::ConductorField" );
		if ( !( mu0In > 0.0 ) )
			throw std::invalid_argument(
				"meq::ConductorField: mu0 must be positive; a zero would make "
				"every conductor silently inert" );
	}

	void ConductorField::add( CurrentFilament const &filament )
	{
		filamentList.push_back( filament );
	}

	void ConductorField::add( Coil const &coil )
	{
		coilList.add( coil );
	}

	std::size_t ConductorField::size() const
	{
		return filamentList.size() + coilList.size();
	}

	bool ConductorField::empty() const
	{
		return filamentList.empty() && coilList.empty();
	}

	std::size_t ConductorField::filamentCount() const
	{
		return filamentList.size();
	}

	std::size_t ConductorField::coilCount() const
	{
		return coilList.size();
	}

	CoilSet const &ConductorField::coils() const
	{
		return coilList;
	}

	void ConductorField::setQuadratureOrder( int order )
	{
		coilList.setQuadratureOrder( order );
	}

	int ConductorField::quadratureOrder() const
	{
		return coilList.quadratureOrder();
	}

	CurrentFilament const &ConductorField::filament( std::size_t index ) const
	{
		if ( index >= filamentList.size() )
		{
			std::ostringstream message;
			message << "meq::ConductorField::filament: index " << index
			        << " is out of range; the set holds "
			        << filamentList.size() << " filaments";
			throw std::out_of_range( message.str() );
		}
		return filamentList[ index ];
	}

	std::vector<CurrentFilament> const &ConductorField::filaments() const
	{
		return filamentList;
	}

	double ConductorField::totalCurrent() const
	{
		double total = 0.0;
		for ( CurrentFilament const &f : filamentList )
			total += f.current();
		return total + coilList.totalCurrent();
	}

	double ConductorField::mu0() const
	{
		return mu0Value;
	}

	double ConductorField::psi( double r, double z ) const
	{
		double total = 0.0;
		for ( CurrentFilament const &f : filamentList )
			total += filamentPsi( f, r, z, mu0Value );
		// The rectangles, through their own set, which carries the quadrature
		// and the same mu0. Empty is exactly zero rather than a special case.
		return total + coilList.psi( r, z );
	}

	void ConductorField::gradPsi( double r, double z,
	                              double &dPsiDr, double &dPsiDz ) const
	{
		dPsiDr = 0.0;
		dPsiDz = 0.0;
		for ( CurrentFilament const &f : filamentList )
		{
			double dr = 0.0;
			double dz = 0.0;
			filamentGradPsi( f, r, z, dr, dz, mu0Value );
			dPsiDr += dr;
			dPsiDz += dz;
		}

		double coilR = 0.0;
		double coilZ = 0.0;
		coilList.gradPsi( r, z, coilR, coilZ );
		dPsiDr += coilR;
		dPsiDz += coilZ;
	}

	void ConductorField::flux( double r, double z,
	                           double &qR, double &qZ ) const
	{
		// NOT a sum of filamentFlux() calls, deliberately. q = ( 1/r ) grad_bar
		// psi, so summing the per-filament fluxes divides by r once per
		// filament and then adds -- which is the same number in exact
		// arithmetic and a different one in floating point, and is NaN for
		// every filament on the axis rather than once. Take the gradient of the
		// sum and divide once.
		double dPsiDr = 0.0;
		double dPsiDz = 0.0;
		gradPsi( r, z, dPsiDr, dPsiDz );
		qR = dPsiDr/r;
		qZ = dPsiDz/r;
	}

	double ConductorField::containment( double centreZ, double rhoGamma ) const
	{
		requireFinite( centreZ, "the centre height",
		               "meq::ConductorField::containment" );
		requireFinite( rhoGamma, "the radius of Gamma",
		               "meq::ConductorField::containment" );
		if ( !( rhoGamma > 0.0 ) )
			throw std::invalid_argument(
				"meq::ConductorField::containment: the radius of Gamma must be "
				"positive; got " + std::to_string( rhoGamma ) );

		// An empty field is contained by everything, and says so with an
		// infinity rather than with a large number a caller might read as a
		// measurement. The same convention ExteriorCoilSet::clearance() uses.
		double farthest = 0.0;

		for ( CurrentFilament const &one : filamentList )
			farthest = std::max( farthest,
			                     std::hypot( one.radius(),
			                                 one.height() - centreZ ) );

		for ( Coil const &one : coilList.coils() )
		{
			// The FARTHEST corner, which is rMax paired with whichever of the
			// two heights is further from the centre. Taking the farthest point
			// rather than the centre is what makes a conductor straddling Gamma
			// report a negative containment.
			double const reach = std::max( std::abs( one.zMin() - centreZ ),
			                               std::abs( one.zMax() - centreZ ) );
			farthest = std::max( farthest, std::hypot( one.rMax(), reach ) );
		}

		if ( empty() )
			return std::numeric_limits<double>::infinity();

		return rhoGamma - farthest;
	}

	void ConductorField::poloidalField( double r, double z,
	                                    double &bR, double &bZ ) const
	{
		if ( r > 0.0 )
		{
			double qR = 0.0;
			double qZ = 0.0;
			flux( r, z, qR, qZ );
			bR = -qZ;
			bZ = qR;
			return;
		}

		// ON THE AXIS, WHERE flux() IS 0/0 AND THE LIMIT IS CLOSED FORM.
		//
		// B_R = -q_z is EXACTLY zero rather than approximately: psi ~ c( z ) r^2
		// for a field regular on the axis, so d_z psi ~ c'( z ) r^2 and
		// q_z ~ c'( z ) r. B_Z = q_r tends to c( z ), which is the on-axis
		// field and is what the two AxisFlux kernels return.
		//
		// SUMMED PER CONDUCTOR AND NOT TAKEN OFF A SUMMED GRADIENT, which is
		// the opposite of what flux() does eight lines above and is right for
		// the opposite reason: there the division by r is what must happen once,
		// here there is no division at all and each conductor's limit is its own
		// closed form.
		bR = 0.0;
		bZ = 0.0;

		for ( CurrentFilament const &one : filamentList )
			bZ += filamentAxisFlux( one, z, mu0Value );

		for ( Coil const &one : coilList.coils() )
			bZ += coilAxisFlux( one, z, coilList.quadratureOrder(), mu0Value );
	}

	bool ConductorField::coincides( double r, double z ) const
	{
		return indexAt( r, z ) >= 0;
	}

	int ConductorField::indexAt( double r, double z ) const
	{
		requireFinite( r, "the field point radius",
		               "meq::ConductorField::indexAt" );
		requireFinite( z, "the field point height",
		               "meq::ConductorField::indexAt" );

		for ( std::size_t i = 0; i < filamentList.size(); ++i )
		{
			CurrentFilament const &f = filamentList[ i ];
			double const dr = r - f.radius();
			double const dz = z - f.height();
			// std::hypot rather than sqrt( dr*dr + dz*dz ): the squares can
			// underflow to zero for the very separations this test exists to
			// resolve, which would report a coincidence that is not there.
			if ( std::hypot( dr, dz ) <= toleranceValue*f.radius() )
				return static_cast<int>( i );
		}
		return -1;
	}

	double ConductorField::coincidenceTolerance() const
	{
		return toleranceValue;
	}

	void ConductorField::setCoincidenceTolerance( double toleranceIn )
	{
		requireFinite( toleranceIn, "the coincidence tolerance",
		               "meq::ConductorField::setCoincidenceTolerance" );
		if ( toleranceIn < 0.0 )
			throw std::invalid_argument(
				"meq::ConductorField::setCoincidenceTolerance: the tolerance "
				"must not be negative; zero is allowed and means exact "
				"equality" );
		toleranceValue = toleranceIn;
	}
}
