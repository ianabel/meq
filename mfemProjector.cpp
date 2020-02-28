#include "mfem.hpp"
#include <fstream>
#include <iostream>
#include <vector>
#include <utility>
#include <list>
#include <boost/multi_array.hpp>
#include <netcdf>

#include "config.hpp"

using namespace mfem;

std::ostream& operator<<( std::ostream& os, Vector const& vec )
{
	os << "(" << vec( 0 ) << "," << vec( 1 ) << ")";
	return os;
}

// Wrap data in boost::multi_array
using BoostArray2D = boost::multi_array<double, 2>;

// Projects grid_func onto a tensor-product grid
int Project( Mesh * mesh, GridFunction const &grid_func, std::vector<double> const &X_coords, std::vector<double> const &Y_coords, BoostArray2D &data )
{
	size_t N_X = X_coords.size();
	size_t N_Y = Y_coords.size();

	FiniteElementSpace const *fes = grid_func.FESpace();
	unsigned int nElems = fes->GetNE();

	// Set all points to NaN -- if a point isn't inside one of the elements this will persist
	for ( unsigned int i = 0; i < N_X; i++ )
		for ( unsigned int j = 0; j < N_Y; j++ )
			data[ i ][ j ] = std::nan( "" );

	unsigned int nPts = 0;
	for ( unsigned int iElem = 0; iElem < nElems; iElem++ )
	{
		ElementTransformation *Tr = mesh->GetElementTransformation( iElem );

		Vector TestVec( 3 );
		TestVec( 2 ) = 0.0; // MFEM is implicitly 3D in places
		for ( unsigned int i = 0; i < N_X; i++ )
		{
			TestVec( 0 ) = X_coords[ i ];
			for ( unsigned int j = 0; j < N_Y; j++ )
			{
				IntegrationPoint ip;
				TestVec( 1 ) = Y_coords[ j ];

				int res = Tr->TransformBack( TestVec, ip );
				if ( res == 0 )
				{
					data[ i ][ j ] = grid_func.GetValue( iElem, ip );
					nPts++;
				}
			}
		}
	}

	return nPts;

}

int main(int argc, char *argv[])
{
	if ( argc != 2 )
	{
		std::cerr << "./mfemProjector <config-file>" << std::endl;
		return -2;
	}

	meq::Configuration meq_conf( argv[ 1 ] );

	std::string const & mesh_file = meq_conf.GetFinalMeshFile();
	std::string const & grid_fn_file = meq_conf.GetPsiFile();

	meq::Domain const & domain = *( meq_conf.GetDomain() );

	const int N_R( 151 );
	const int N_Z( 151 );
	double R_min = domain.RMin;
	double R_max = domain.RMax;
	double Z_min = domain.ZMin;
	double Z_max = domain.ZMax;
	std::vector<double> R_pts( N_R );
	std::vector<double> Z_pts( N_Z );
	for ( unsigned int i=0; i < N_R; i++ )
		R_pts[ i ] = R_min + i*( R_max - R_min )/( N_R - 1 );
	for ( unsigned int i=0; i < N_Z; i++ )
		Z_pts[ i ] = Z_min + i*( Z_max - Z_min )/( N_Z - 1 );


   //    Read the mesh from the given mesh file. We can handle triangular,
   //    quadrilateral, tetrahedral, hexahedral, surface and volume meshes with
   //    the same code.
   Mesh *mesh = new Mesh(mesh_file.c_str(), 1, 1);

	//    Read the grid function
	std::fstream grid_fs( grid_fn_file.c_str() );
	GridFunction grid_func( mesh, grid_fs );
	grid_fs.close();

	BoostArray2D data( boost::extents[ R_pts.size() ][ Z_pts.size() ], boost::c_storage_order() );

	int nFound = Project( mesh, grid_func, R_pts, Z_pts, data );
	std::cout << " Of " << R_pts.size() * Z_pts.size() << " possible points, " << nFound << " were found in the mesh" << std::endl;

	{
		netCDF::NcFile data_file( meq_conf.GetNetCDFFile(), netCDF::NcFile::FileMode::replace );
		netCDF::NcDim R_dim = data_file.addDim( "R", R_pts.size() );
		netCDF::NcDim Z_dim = data_file.addDim( "Z", Z_pts.size() );
		netCDF::NcVar R_var = data_file.addVar( "R", netCDF::NcDouble(), R_dim );
		netCDF::NcVar Z_var = data_file.addVar( "Z", netCDF::NcDouble(), Z_dim );
		R_var.putVar( R_pts.data() );
		Z_var.putVar( Z_pts.data() );

		std::vector<netCDF::NcDim> psi_dims = { R_dim, Z_dim };
		netCDF::NcVar psi = data_file.addVar( "psi", netCDF::NcDouble(), psi_dims );

		psi.putVar( data.data() );
		data_file.close();
	}


	delete mesh;
	return 0;
}
