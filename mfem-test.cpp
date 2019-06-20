#include "mfem.hpp"
#include <fstream>
#include <iostream>

using namespace mfem;

/*
double ComputeElementError( unsigned int i ) 
	{
		mfem::ElementTransformation *trans = mesh->GetElementTransformation( i );
		mfem::Element const *K = mesh->GetElement( i );
		mfem::FiniteElement const *u_fe = u_space->GetFE( i );
		mfem::FiniteElement const *q_fe = q_space->GetFE( i );

		double eta_1 = 0;
		double eta_2 = 0;

		unsigned int order = 2 * q_fe->GetOrder() + 3;
		mfem::IntegrationRule const& CellIntegrator = mfem::IntRules.Get( K->GetType(), order );
		

		double h_K = ElementDiameter( mesh, K );

		// Sum over integration points
		for (int j=0; j < CellIntegrator.GetNPoints(); j++)
		{
			const IntegrationPoint &ip = CellIntegrator.IntPoint( j );
			trans->SetIntPoint( &ip );

			double eta_1_tmp;
			// eta_1 = h_K^2 * || F_RHS + div q ||^2
			eta_1_tmp = ( rhs.Eval( *trans, ip ) + q_sol.GetDivergence( *trans ) );
			eta_1 += h_K * h_K * eta_1_tmp * eta_1_tmp;

			Vector eta_2_tmp;
			// eta_2 = || q - kappa*(grad u) ||^2
			Vector q_val;
			Vector GradU;
			q_sol.GetVectorValue( i, ip, q_val );
			u_sol.GetGradient( *trans, GradU );
			eta_2_tmp = q_val - GradU;	
			eta_2 += eta_2_tmp * eta_2_tmp;
		}

		double eta_3 = 0;
		double eta_4 = 0;
		double eta_5 = 0;

		mfem::Array<int> faces,orientations;
		mesh->GetElementEdges( i, faces, orientations );

		// Sum over all faces
		for (int j=0; j < faces.Size(); j++)
		{
			mfem::Element const *e = mesh->GetFace( faces[ j ] );
			double h_e = ElementDiameter( mesh, e );

			mfem::FaceElementTransformations *feTrans = mesh->GetFaceElementTransformations( faces[ j ] );

			// Integrate over e
			mfem::IntegrationRule const& EdgeIntegrator = mfem::IntRules.Get( e->GetType(), order );
			for (int k=0; k < EdgeIntegrator.GetNPoints(); k++)
			{
				const IntegrationPoint& ip = EdgeIntegrator.IntPoint( k );
				IntegrationPoint eip1,eip2; 
				feTrans->Loc1.Transform( ip, eip1 );
				feTrans->Loc2.Transform( ip, eip2 );
				feTrans->Elem1->SetIntPoint( &eip1 );
				feTrans->Elem2->SetIntPoint( &eip2 );
				feTrans->Face->SetIntPoint( &ip );

				// eta_3 = (1/2) * h_e * || (q_plus - q_minus) . n ||^2
				Vector q_plus,q_minus;
				q_sol.GetVectorValue( feTrans->Elem1No, eip1, q_plus );
				q_sol.GetVectorValue( feTrans->Elem2No, eip2, q_minus );
				
				Vector normal( mesh->SpaceDimension() );	
				mfem::CalcOrtho( feTrans->Face->Jacobian(), normal );
				normal /= normal.Norml2();

				double q_jump = q_plus * normal - q_minus * normal;
				eta_3 += 0.5 * h_e * q_jump * q_jump;

				// eta_4 = .5 * h_e^{-1} || u_plus - u_minus ||^2
				double u_plus,u_minus;

				u_plus  = u_sol.GetValue( feTrans->Elem1No, eip1 );
				u_minus = u_sol.GetValue( feTrans->Elem2No, eip2 );
				
				eta_4 += ( .5 / h_e ) * ( u_plus - u_minus ) * ( u_plus - u_minus );

				// eta_5 = (.5/h_e) * || lambda - u ||^2
				double lambda_val,u_val;
				lambda_val = lambda.GetValue( faces[ j ], ip );
				if ( feTrans->Elem1No == i )
					u_val = u_sol.GetValue( feTrans->Elem1No, eip1 );
				else if ( feTrans->Elem2No == i )
					u_val = u_sol.GetValue( feTrans->Elem2No, eip2 );
				else 
					throw new std::logic_error( "Element is neither of the ones attached to the face. Wat." );

				eta_5 += ( .5/h_e )*( lambda_val - u_val )*( lambda_val - u_val );

			}
		}
		return eta_1 + eta_2 + eta_3 + eta_4 + eta_5;
	}
*/

double sinX( mfem::Vector const& pt )
{
	double pi = 3.14159265358979323844;
	return ::sin( pi * pt( 0 ) / 2 );
}

int main(int argc, char *argv[])
{
	const char mesh_file[] = "trivial.msh";
	int order;
	if ( argc > 1 )
		order = ::atoi( argv[ 1 ] );
	else
		order = 1;
   Mesh *mesh = new Mesh(mesh_file, 1, 1);
	int dim = mesh->Dimension();

	FiniteElementCollection *fec_cells = new DG_FECollection( order, dim );
   FiniteElementSpace *fe_cell_space = new FiniteElementSpace( mesh, fec_cells );
	FiniteElementSpace *fe_cell_vector_space = new FiniteElementSpace( mesh, fec_cells, dim );
	
	FiniteElementCollection *fec_edges = new DG_Interface_FECollection( order, dim );
	FiniteElementSpace *fe_edge_space = new FiniteElementSpace( mesh, fec_edges );

	GridFunction lambda( fe_edge_space );
	FunctionCoefficient fc( sinX );
	lambda.ProjectCoefficientSkeletonDG( fc );


	std::cout << "Mesh is dimension " << dim << std::endl;
	std::cout << "Number of finite element unknowns in cells for scalar: " << fe_cell_space->GetNDofs() << std::endl;
	std::cout << "Number of finite element unknowns in cells for vector: " << fe_cell_vector_space->GetNDofs() << std::endl;
	std::cout << "Number of boundary elements: " << fe_cell_space->GetNBE() << std::endl;

	std::cout << std::endl;

	std::cout << "Number of scalar elements: " << fe_cell_space->GetNE() << std::endl;
	std::cout << "Number of vector elements: " << fe_cell_vector_space->GetNE() << std::endl;
	std::cout << "Number of edge elements: " << fe_edge_space->GetNBE() << std::endl;

	int elemId = 0;
	
	mfem::Array<int> faces,orientations;
	mesh->GetElementEdges( elemId, faces, orientations );

	// Sum over all faces
	for (int j=0; j < faces.Size(); j++)
	{
		mfem::Element const *e = mesh->GetFace( faces[ j ] );
		mfem::FaceElementTransformations *feTrans = mesh->GetFaceElementTransformations( faces[ j ] );

		std::cout << "The face we are looking at is ";
		IntegrationPoint testPoint;
		testPoint.x = 0.0; testPoint.y = 0.0; testPoint.z = 0.0; testPoint.weight = 0.0;
		Vector testActPt( 2 );
		feTrans->Face->Transform( testPoint, testActPt );
		std::cout << " (" << testActPt( 0 ) << ", " << testActPt( 1 ) << ") to (";
		testPoint.x = 1.0; testPoint.y = 0.0; testPoint.z = 0.0; testPoint.weight = 0.0;
		feTrans->Face->Transform( testPoint, testActPt );
		std::cout << testActPt( 0 ) << ", " << testActPt( 1 ) << ")" << std::endl;




		mfem::IntegrationRule const& EdgeIntegrator = mfem::IntRules.Get( e->GetType(), order );
		for (int k=0; k < EdgeIntegrator.GetNPoints(); k++)
		{
			const IntegrationPoint& ip = EdgeIntegrator.IntPoint( k );
			IntegrationPoint eip1,eip2; 
			feTrans->Loc1.Transform( ip, eip1 );
			feTrans->Loc2.Transform( ip, eip2 );
			feTrans->Elem1->SetIntPoint( &eip1 );
			feTrans->Elem2->SetIntPoint( &eip2 );
			feTrans->Face->SetIntPoint( &ip );

			Vector pt( 2 );
			feTrans->Face->Transform( ip, pt );
			std::cout << "Ref Point is (" << ip.x << ", " << ip.y << ") Target Point is  (" << pt( 0 ) << ", " << pt( 1 ) << ")" << std::endl;

			double lambda_val,u_val;
			lambda_val = lambda.GetValue( feTrans->Elem1No, eip1 );
			std::cout << " The value of lambda at (" << pt( 0 ) << ", " << pt( 1 ) << ") is " << lambda_val << std::endl;

		}
	}

	return 0;
}


