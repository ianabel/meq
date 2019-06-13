
#include "mfem.hpp"

double ElementDiameter( Element* e )
{
	unsigned int nVerts = e.GetNVertices();
	mfem::Array<int> VertIndices( nVerts );
	e.GetVertices( VertIndices );
	switch ( e.GetType() )
	{
		case mfem::Element::Type::SEGMENT:
			Vector a,b;
			return <<>>;
		case mfem::Element::Type::Triangle:
			return <<>>;
		default:
			throw new std::logic_error( "Unsupported." );
	}
	return -1.0;
}

class CockburnZhangEstimator : public mfem::ErrorEstimator
{
	protected:
		mfem::Mesh *mesh;
		mfem::FiniteElementSpace *q_space;
		mfem::FiniteElementSpace *u_space;

		mfem::Vector localErrors;
		GridFunction &q_sol,&u_sol;
		Coefficient &kappa;
		Coefficient &rhs;
		bool valid;
	public:
	CockburnZhangEstimator( GridFunction &q, GridFunction &u, Coefficient &kappa_ref, Coefficient &rhs_ref )
		: q_sol( q ), u_sol( u ), kappa( kappa_ref ), rhs( rhs_ref ),valid( false )
	{
		q_space = q.FESpace();
		u_space = u.FESpace();
		mesh = u_space->GetMesh();
	};
	

	double ComputeElementError( unsigned int i ) 
	{
		mfem::IntegrationRule CellIntegrator;
		mfem::ElementTransformation trans = mesh->GetElementTransformation( i );
		mfem::FiniteElement u_fe = u_space->GetFE( i );
		mfem::FiniteElement q_fe = q_space->GetFE( i );

		double eta_1 = 0;
		double eta_2 = 0;

		double h_K = 

		for ()
		{
			// Integrate over K.
		}

		double eta_3 = 0;
		double eta_4 = 0;
		double eta_5 = 0;

		for ()
		{
			// Sum over all faces
			for ()
			{
				// Integrate over e
			}
		}
	}

   /// Get a Vector with all element errors.
   virtual const Vector &GetLocalErrors() 
	{
		if ( valid )
			return localErrors;

		unsigned int nElems = mesh->GetNE();
		localErrors.SetSize( nElems );
		for ( unsigned int j=0; j <= nElems; j++ )
			localErrors( j ) = ComputeElementError( j );

		valid = true;
		return localErrors;
	}

   /// Force recomputation of the estimates on the next call to GetLocalErrors.
   virtual void Reset() 
	{
		valid = false;
	}
}
