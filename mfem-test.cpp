#include "mfem.hpp"
#include "HDGGSIntegrator.hpp"
#include "FreeBoundary.hpp"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <boost/format.hpp>

using namespace mfem;

void qFun( mfem::Vector const & in , mfem::Vector & out )
{
	double x = in( 0 );
	double y = in( 1 );

	out( 0 ) = ( x - 0.55 )/.45;
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
	auto xform = []( const Vector& in, Vector& out ) { 
		constexpr double R_min = 0.1;
		out( 1 ) = in( 1 );
		out( 0 ) = R_min + in( 0 )*( 1 - R_min );
	};
	mesh->Transform( xform );


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
	
	mfem::Vector test_pt( 2 );
	double bpsi;

	boost::format output_form(" BPsi[%1%,%2%] = %3$.12d \n" );
	
	test_pt( 0 ) = .1;
	test_pt( 1 ) = .5;

	bpsi = BoundaryPsi( fe_cell_vector_space, q, test_pt );

	std::cout << output_form % test_pt( 0 ) % test_pt( 1 ) % bpsi;

	test_pt( 0 ) = .1;
	test_pt( 1 ) = 0;

	bpsi = BoundaryPsi( fe_cell_vector_space, q, test_pt );

	std::cout << output_form % test_pt( 0 ) % test_pt( 1 ) % bpsi;

	test_pt( 0 ) = .1;
	test_pt( 1 ) = 1;

	bpsi = BoundaryPsi( fe_cell_vector_space, q, test_pt );

	std::cout << output_form % test_pt( 0 ) % test_pt( 1 ) % bpsi;

	test_pt( 0 ) = .5;
	test_pt( 1 ) = 1;

	bpsi = BoundaryPsi( fe_cell_vector_space, q, test_pt );

	std::cout << output_form % test_pt( 0 ) % test_pt( 1 ) % bpsi;
	
	test_pt( 0 ) = .5;
	test_pt( 1 ) = 0;

	bpsi = BoundaryPsi( fe_cell_vector_space, q, test_pt );

	std::cout << output_form % test_pt( 0 ) % test_pt( 1 ) % bpsi;

	test_pt( 0 ) = 1;
	test_pt( 1 ) = .5;

	bpsi = BoundaryPsi( fe_cell_vector_space, q, test_pt );

	std::cout << output_form % test_pt( 0 ) % test_pt( 1 ) % bpsi;

	test_pt( 0 ) = 1;
	test_pt( 1 ) = 0;

	bpsi = BoundaryPsi( fe_cell_vector_space, q, test_pt );

	std::cout << output_form % test_pt( 0 ) % test_pt( 1 ) % bpsi;

	test_pt( 0 ) = 1;
	test_pt( 1 ) = 1;

	bpsi = BoundaryPsi( fe_cell_vector_space, q, test_pt );

	std::cout << output_form % test_pt( 0 ) % test_pt( 1 ) % bpsi;










	return 0;
}


