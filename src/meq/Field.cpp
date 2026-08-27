#include "Field.hpp"

#include <stdexcept>

namespace meq
{
	void poloidalField( mfem::GridFunction const &q, mfem::GridFunction &field )
	{
		mfem::FiniteElementSpace *space =
			const_cast<mfem::FiniteElementSpace *>( q.FESpace() );
		if ( space == nullptr )
			throw std::logic_error( "meq::poloidalField: the flux has no finite element space" );
		if ( space->GetVDim() != 2 )
			throw std::logic_error( "meq::poloidalField: the flux must have vdim 2" );

		if ( field.FESpace() == nullptr )
			field.SetSpace( space );
		else if ( field.FESpace() != space )
			throw std::logic_error( "meq::poloidalField: the field and the flux must share a space" );

		// B_R = -q_z, B_Z = +q_r: a component swap and one sign. Done by dof
		// rather than by interpolation, so the field is exact wherever q is and
		// carries q's own k+1 convergence with no quadrature in between.
		int const perComponent = space->GetNDofs();
		mfem::Ordering::Type const ordering = space->GetOrdering();

		for ( int i = 0; i < perComponent; ++i )
		{
			int const rIndex = ( ordering == mfem::Ordering::byNODES )
			                   ? i : 2*i;
			int const zIndex = ( ordering == mfem::Ordering::byNODES )
			                   ? i + perComponent : 2*i + 1;

			double const qR = q( rIndex );
			double const qZ = q( zIndex );

			field( rIndex ) = -qZ;   // B_R
			field( zIndex ) = qR;    // B_Z
		}
	}

	PoloidalFieldCoefficient::PoloidalFieldCoefficient( mfem::GridFunction const &qIn,
	                                                    int componentIn )
		: q( qIn ), component( componentIn ), value( 2 )
	{
		if ( componentIn != 0 && componentIn != 1 )
			throw std::logic_error( "meq::PoloidalFieldCoefficient: the component must be 0 for B_R or 1 for B_Z" );
	}

	double PoloidalFieldCoefficient::Eval( mfem::ElementTransformation &tr,
	                                       mfem::IntegrationPoint const &ip )
	{
		q.GetVectorValue( tr, ip, value );
		// value( 0 ) is q_r and value( 1 ) is q_z.
		return ( component == 0 ) ? -value( 1 ) : value( 0 );
	}
}
