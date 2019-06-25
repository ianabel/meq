#include "mfem.hpp"
#include <fstream>
#include <iostream>

using namespace mfem;

double GetTraceValue( GridFunction const& lambda, int i, const IntegrationPoint& ip )
{
	mfem::Array<int> dofs;
	mfem::FiniteElementSpace const *fes = lambda.FESpace();
   fes->GetFaceDofs(i, dofs);
	int vdim = fes->GetVDim();
   fes->DofsToVDofs(vdim-1, dofs);
	mfem::Vector DofVal(dofs.Size()), LocVec;
   const FiniteElement *fe = fes->GetFaceElement(i);
   MFEM_ASSERT(fe->GetMapType() == FiniteElement::VALUE, "invalid FE map type");
   fe->CalcShape(ip, DofVal);
   lambda.GetSubVector(dofs, LocVec);

   return (DofVal * LocVec);
}

double sinXY( mfem::Vector const& pt )
{
	double pi = 3.14159265358979323844;
	return ::sin( pi * pt( 0 ) * pt( 1 )  );
}

int main(int argc, char *argv[])
{
	int order = 5;
   Mesh *mesh = new Mesh(10, 10, Element::Type::TRIANGLE );
	int dim = mesh->Dimension();

	FiniteElementCollection *fec_cells = new DG_FECollection( order, dim );
	FiniteElementCollection *fec_edges = new DG_Interface_FECollection( order, dim );
	FiniteElementSpace *fe_edge_space = new FiniteElementSpace( mesh, fec_edges );

	GridFunction lambda( fe_edge_space );
	FunctionCoefficient fc( sinXY );
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
			lambda_val = GetTraceValue( lambda, feTrans->Face->ElementNo, ip );
			std::cout << " The value of lambda at (" << pt( 0 ) << ", " << pt( 1 ) << ") is " << lambda_val << std::endl;

		}
		break;
	}

	return 0;
}


