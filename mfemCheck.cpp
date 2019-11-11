#include "mfem.hpp"
#include <fstream>
#include <iostream>
#include <vector>
#include <utility>
#include <list>

using namespace mfem;

int main(int argc, char *argv[])
{
	if ( argc != 4 )
	{
		std::cerr << "./mfemCheck <prefix_1> <prefix_2>" << std::endl;
		return -2;
	}

	std::string 

	std::string mesh_file = argv[ 1 ];
	std::string grid_fn_1_file = argv[ 2 ];
	std::string grid_fn_2_file = argv[ 3 ];

   Mesh *mesh = new Mesh(mesh_file.c_str(), 1, 1);

	//    Read the grid functions
	std::fstream grid_fs( grid_fn_1_file.c_str() );
	GridFunction grid_func_1( mesh, grid_fs );
	grid_fs.close();
	grid_fs.open( grid_fn_2_file.c_str() );
	GridFunction grid_func_2( mesh, grid_fs );
	grid_fs.close();


	GridFunctionCoefficient gf_2( &grid_func_2 );

	double l1_err = grid_func_1.ComputeL1Error( gf_2 );
	double l2_err = grid_func_1.ComputeL2Error( gf_2 );
	double l_inf_err = grid_func_1.ComputeMaxError( gf_2 );

	std::cout << "\t|| u_1 - u_2 ||_1 = " << l1_err << "\n";
	std::cout << "\t|| u_1 - u_2 ||_2 = " << l2_err << "\n";
	std::cout << "\t|| u_1 - u_2 ||_2 = " << l_inf_err << "\n";


	delete mesh;
	return 0;
}
