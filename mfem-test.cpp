#include "mfem.hpp"
#include <fstream>
#include <iostream>

using namespace mfem;


double sinX( mfem::Vector const& pt )
{
	double pi = 3.14159265358979323844;
	return ::sin( pi * pt( 0 ) / 2 );
}

void qFun( mfem::Vector const & in , mfem::Vector & out )
{
	double x = in( 0 );
	double y = in( 1 );

	out( 0 ) = 2.0 * ( x - 0.5 );
	out( 1 ) = 2.0 * ( y - 0.5 );

	return;
}





double u_fun( mfem::Vector const& in )
{
	double y = in( 1 );
	return 4.0*( y - 0.5 )*( y - 0.5 );
}

int main(int argc, char *argv[])
{
	int order = 4;
   Mesh *mesh = new Mesh(4, 4, mfem::Element::Type::TRIANGLE );
	int dim = mesh->Dimension();

	FiniteElementCollection *fec_cells = new DG_FECollection( order, dim );
   FiniteElementSpace *fe_cell_space = new FiniteElementSpace( mesh, fec_cells );
	FiniteElementSpace *fe_cell_vector_space = new FiniteElementSpace( mesh, fec_cells, dim );
	
	FiniteElementCollection *fec_edges = new DG_Interface_FECollection( order, dim );
	FiniteElementSpace *fe_edge_space = new FiniteElementSpace( mesh, fec_edges );

	GridFunction lambda( fe_edge_space );
	Array<int> attr( mesh->bdr_attributes.Max() );
	attr = 1;
	ConstantCoefficient one( 1.0 );
	std::cout << attr.Size() << std::endl;
	lambda.ProjectBdrCoefficient( one, attr );

	GridFunction u( fe_cell_space );
	FunctionCoefficient u_c( u_fun );
	u.ProjectCoefficient( u_c );

	GridFunction q( fe_cell_vector_space );
	VectorFunctionCoefficient q_c( dim, qFun );
	q.ProjectCoefficient( q_c );
	

	std::cout << "Mesh is dimension " << dim << std::endl;
	std::cout << "Number of finite element unknowns in cells for scalar: " << fe_cell_space->GetNDofs() << std::endl;
	std::cout << "Number of finite element unknowns in cells for vector: " << fe_cell_vector_space->GetNDofs() << std::endl;
	std::cout << "Number of boundary elements: " << fe_cell_space->GetNBE() << std::endl;

	std::cout << std::endl;

	std::cout << "Number of scalar elements: " << fe_cell_space->GetNE() << std::endl;
	std::cout << "Number of vector elements: " << fe_cell_vector_space->GetNE() << std::endl;
	std::cout << "Number of edge elements: " << fe_edge_space->GetNBE() << std::endl;

	LinearForm biForm( fe_cell_space );
	biForm.AddBdrFaceIntegrator( new HDGBoundaryTraceIntegrator( one ) );
	biForm.Assemble();

	std::cout << "Integral of 1 * 1 around the boundary is " << biForm( u ) << std::endl;

	LinearForm biVForm( fe_cell_vector_space );
	biVForm.AddBdrFaceIntegrator( new HDGBoundaryNormalTraceIntegrator( one ) );
	biVForm.Assemble();	

	std::cout << "Integral of n . q around the boundary is " << biVForm( q ) << std::endl;

	return 0;
}


