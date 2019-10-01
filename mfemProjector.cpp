#include "mfem.hpp"
#include <fstream>
#include <iostream>
#include <vector>
#include <utility>
#include <list>

using namespace mfem;

std::ostream& operator<<( std::ostream& os, Vector const& vec )
{
	os << "(" << vec( 0 ) << "," << vec( 1 ) << ")";
	return os;
}

int main(int argc, char *argv[])
{
	if ( argc != 3 )
	{
		std::cerr << "./mfemProjector <mesh> <gridfunction>" << std::endl;
		return -2;
	}

	std::string mesh_file = argv[ 1 ];
	std::string grid_fn_file = argv[ 2 ];



	std::list< std::pair< double, double> > outputPoints{
		{0.01,0.75},{0.01,0.0},{0.01,-0.75}, {0.5,0.5},{0.5,0.0},{0.5,-0.5}, {0.75,0.75},{0.75,0.0},{0.75,-0.75}
	};



   //    Read the mesh from the given mesh file. We can handle triangular,
   //    quadrilateral, tetrahedral, hexahedral, surface and volume meshes with
   //    the same code.
   Mesh *mesh = new Mesh(mesh_file.c_str(), 1, 1);

	//    Read the grid function
	std::fstream grid_fs( grid_fn_file.c_str() );
	GridFunction grid_func( mesh, grid_fs );
	grid_fs.close();

	FiniteElementSpace *fes = grid_func.FESpace();

	unsigned int nElems = fes->GetNE();

	for ( unsigned int iElem = 0; iElem < nElems; iElem++ )
	{
		ElementTransformation *Tr = mesh->GetElementTransformation( iElem );
		IntegrationRule ir( 0 );
		std::vector<std::pair<double,double>> pointsInElem;


		Vector TestVec( 3 );
		for ( std::pair<double, double> const &x : outputPoints )
		{
			IntegrationPoint ip;
			TestVec( 0 ) = x.first; TestVec( 1 ) = x.second; TestVec( 2 ) = 0.0;
			int res = Tr->TransformBack( TestVec, ip );
			if ( res == 0 )
			{
				ir.Append( ip );
				pointsInElem.push_back( x );
			}
		}


		unsigned int nPts = ir.GetNPoints();
		if ( nPts != pointsInElem.size() )
		{
			std::cerr << "Bork" << std::endl; return -1;
		}

		Vector values( nPts );
		grid_func.GetValues( iElem, ir, values );
		for ( unsigned int j = 0 ; j < nPts; j++ )
		{
			std::cout << pointsInElem[ j ].first << "\t" << pointsInElem[ j ].second << "\t" << values( j ) << std::endl;
		}

	}

	delete mesh;
	return 0;
}
