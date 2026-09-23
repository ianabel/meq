#include "ConductorField.hpp"

#include "ConductorStore.hpp"

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

	std::uint64_t ConductorField::digest() const
	{
		/*
		 * THE COUNTS GO IN AS WELL AS THE VALUES, which is not redundant: a
		 * filament of zero current and no filament at all give the same psi_c
		 * and must still give different digests, because the quadrature cache
		 * is indexed by element and the SECOND of those can legitimately be
		 * built on a different mesh. Mixing the count first is also what stops
		 * a set of three filaments and a set of one with three times the data
		 * colliding.
		 */
		std::uint64_t digest = digestSeed();

		digest = digestAppend(
			digest, static_cast< std::int64_t >( filaments().size() ) );
		for ( CurrentFilament const &filament : filaments() )
		{
			digest = digestAppend( digest, filament.radius() );
			digest = digestAppend( digest, filament.height() );
			digest = digestAppend( digest, filament.current() );
		}

		digest = digestAppend(
			digest, static_cast< std::int64_t >( coils().size() ) );
		for ( std::size_t i = 0; i < coils().size(); ++i )
		{
			Coil const &conductor = coils().coil( i );
			digest = digestAppend( digest, conductor.centreR() );
			digest = digestAppend( digest, conductor.centreZ() );
			digest = digestAppend( digest, conductor.halfWidth() );
			digest = digestAppend( digest, conductor.halfHeight() );
			digest = digestAppend( digest, conductor.current() );
		}

		digest = digestAppend( digest, mu0Value );

		// The rectangles' cross-section rule changes psi_c by more than round
		// off -- it is the model, not the arithmetic -- so it belongs here.
		digest = digestAppend(
			digest, static_cast< std::int64_t >( quadratureOrder() ) );

		return digest;
	}

	double ConductorField::psi( double radius, double z ) const
	{
		/*
		 * EVEN IN R, AND THAT IS AN ANALYTIC CONTINUATION RATHER THAN A CLAMP.
		 *
		 * psi = R A_phi. Under R -> -R at fixed z the point is the same
		 * physical point rotated by pi in phi, so phi-hat reverses and A_phi
		 * changes sign with it; the product does not. So psi_c( -R, z ) is
		 * psi_c( R, z ) EXACTLY, and reflecting is the unique smooth extension
		 * of the flux across the axis rather than a convenient substitute for
		 * one. Near the axis psi_c ~ c( z ) R^2, so the value this returns for
		 * a point a hair past R = 0 is a hair above zero, which is what the
		 * physics says it should be.
		 *
		 * IT IS HERE AND NOT AT A CALLER BECAUSE IT IS A PROPERTY OF THE
		 * FUNCTION, and because MEQ has two evaluations of psi_c that
		 * legitimately extrapolate off the half-plane and would each otherwise
		 * need their own rule: the exterior datum, whose transfer paths target
		 * a Gamma that MEETS the axis, and CriticalPointFinder, whose
		 * element-local Newton is allowed to leave its element on purpose.
		 *
		 * gradPsi(), flux() and poloidalField() still REFUSE, and the
		 * asymmetry is the parity rather than an oversight: d_r psi is ODD
		 * where psi is even, so the three of them continue with a sign that
		 * differs per component, and a caller holding a vector cannot apply
		 * one rule to both of its entries. meq::CriticalPointFinder::totalFlux
		 * is the seam that meets this for q and it abandons the evaluation,
		 * for the reason recorded against it.
		 */
		double const rho = std::abs( radius );

		double total = 0.0;
		for ( CurrentFilament const &f : filamentList )
			total += filamentPsi( f, rho, z, mu0Value );
		// The rectangles, through their own set, which carries the quadrature
		// and the same mu0. Empty is exactly zero rather than a special case.
		return total + coilList.psi( rho, z );
	}

	void ConductorField::gradPsi( double radius, double z,
	                              double &dPsiDr, double &dPsiDz ) const
	{
		dPsiDr = 0.0;
		dPsiDz = 0.0;
		for ( CurrentFilament const &f : filamentList )
		{
			double dr = 0.0;
			double dz = 0.0;
			filamentGradPsi( f, radius, z, dr, dz, mu0Value );
			dPsiDr += dr;
			dPsiDz += dz;
		}

		double coilR = 0.0;
		double coilZ = 0.0;
		coilList.gradPsi( radius, z, coilR, coilZ );
		dPsiDr += coilR;
		dPsiDz += coilZ;
	}

	void ConductorField::flux( double radius, double z,
	                           double &qR, double &qZ ) const
	{
		// NOT a sum of filamentFlux() calls, deliberately. q = ( 1/R ) grad_bar
		// psi, so summing the per-filament fluxes divides by R once per
		// filament and then adds -- which is the same number in exact
		// arithmetic and a different one in floating point, and is NaN for
		// every filament on the axis rather than once. Take the gradient of the
		// sum and divide once.
		double dPsiDr = 0.0;
		double dPsiDz = 0.0;
		gradPsi( radius, z, dPsiDr, dPsiDz );
		qR = dPsiDr/radius;
		qZ = dPsiDz/radius;
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
			// The FARTHEST corner, which is R_max paired with whichever of the
			// two heights is further from the centre. Taking the farthest point
			// rather than the centre is what makes a conductor straddling Gamma
			// report a negative containment.
			double const reach = std::max( std::abs( one.zMin() - centreZ ),
			                               std::abs( one.zMax() - centreZ ) );
			farthest = std::max( farthest, std::hypot( one.maxRadius(), reach ) );
		}

		if ( empty() )
			return std::numeric_limits<double>::infinity();

		return rhoGamma - farthest;
	}

	void ConductorField::poloidalField( double radius, double z,
	                                    double &bR, double &bZ ) const
	{
		if ( radius > 0.0 )
		{
			double qR = 0.0;
			double qZ = 0.0;
			flux( radius, z, qR, qZ );
			bR = -qZ;
			bZ = qR;
			return;
		}

		// AND R < 0 IS NOT "ON THE AXIS", WHICH THE R > 0 TEST ABOVE LEAVES IT
		// INDISTINGUISHABLE FROM. B is a VECTOR -- B_R = -q_z is odd under
		// R -> -R where B_Z = q_r is even -- so it continues past the axis with
		// a sign that differs between its two entries, exactly as gradPsi() and
		// flux() do, and this refuses for their reason. psi() is the one entry
		// point here that continues, because it is a scalar and even.
		if ( radius < 0.0 )
			throw std::invalid_argument(
				"meq::ConductorField::poloidalField: the field point radius "
				"must not be negative; R = 0 is the axis and gives the closed-"
				"form limit. psi() continues past the axis by reflection and B "
				"cannot, because B_R and B_Z carry opposite parities" );

		// ON THE AXIS, WHERE flux() IS 0/0 AND THE LIMIT IS CLOSED FORM.
		//
		// B_R = -q_z is EXACTLY zero rather than approximately: psi ~ c( z ) R^2
		// for a field regular on the axis, so d_z psi ~ c'( z ) R^2 and
		// q_z ~ c'( z ) R. B_Z = q_r tends to c( z ), which is the on-axis
		// field and is what the two AxisFlux kernels return.
		//
		// SUMMED PER CONDUCTOR AND NOT TAKEN OFF A SUMMED GRADIENT, which is
		// the opposite of what flux() does eight lines above and is right for
		// the opposite reason: there the division by R is what must happen once,
		// here there is no division at all and each conductor's limit is its own
		// closed form.
		bR = 0.0;
		bZ = 0.0;

		for ( CurrentFilament const &one : filamentList )
			bZ += filamentAxisFlux( one, z, mu0Value );

		for ( Coil const &one : coilList.coils() )
			bZ += coilAxisFlux( one, z, coilList.quadratureOrder(), mu0Value );
	}

	bool ConductorField::coincides( double radius, double z ) const
	{
		return indexAt( radius, z ) >= 0;
	}

	int ConductorField::indexAt( double radius, double z ) const
	{
		requireFinite( radius, "the field point radius",
		               "meq::ConductorField::indexAt" );
		requireFinite( z, "the field point height",
		               "meq::ConductorField::indexAt" );

		for ( std::size_t i = 0; i < filamentList.size(); ++i )
		{
			CurrentFilament const &f = filamentList[ i ];
			double const dr = radius - f.radius();
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
