
#include <netcdf>
#include <boost/multi_array.hpp>
using BoostArray2D = boost::multi_array<double, 2>;

#include "meq.hpp"

// Code for NetCDF interface
//

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

int ProjectBField( Mesh * mesh, GridFunction const &q_field, std::vector<double> const &X_coords, std::vector<double> const &Y_coords, BoostArray2D &B_R_data, BoostArray2D &B_Z_data )
{
	size_t N_R = R_coords.size();
	size_t N_Z = Z_coords.size();

	FiniteElementSpace const *fes = q_field.FESpace();
	unsigned int nElems = fes->GetNE();

	// Set all points to NaN -- if a point isn't inside one of the elements this will persist
	for ( unsigned int i = 0; i < N_R; i++ )
		for ( unsigned int j = 0; j < N_Z; j++ )
		{
			B_R_data[ i ][ j ] = std::nan( "" );
			B_Z_data[ i ][ j ] = std::nan( "" );
		}

	unsigned int nPts = 0;
	for ( unsigned int iElem = 0; iElem < nElems; iElem++ )
	{
		ElementTransformation *Tr = mesh->GetElementTransformation( iElem );

		Vector TestVec( 3 );
		TestVec( 2 ) = 0.0; // MFEM is implicitly 3D in places
		for ( unsigned int i = 0; i < N_R; i++ )
		{
			TestVec( 0 ) = R_coords[ i ];
			for ( unsigned int j = 0; j < N_Z; j++ )
			{
				IntegrationPoint ip;
				TestVec( 1 ) = Z_coords[ j ];

				int res = Tr->TransformBack( TestVec, ip );
				if ( res == 0 )
				{
					mfem::Vector q_val( 2 );
					q_field.GetVectorValue( iElem, ip, q_val );
					double R = R_coords[ i ];

					// B_pol = Grad psi x Grad phi
					// q = Grad Psi / R
					B_R_data[ i ][ j ] = -R*q_val( 1 );
					B_R_data[ i ][ j ] =  R*q_val( 0 );
					nPts++;
				}
			}
		}
	}

	return nPts;

}

void WriteNetCDF( const meq::NetCDFConfig &netcdf_conf, const meq::Configuration &meq_config, const meq::Solution &soln )
{
	std::vector<double> const &R_pts = netcdf_conf.R_dim;
	std::vector<double> const &Z_pts = netcdf_conf.Z_dim;



	netCDF::NcFile data_file( netcdf_file, netCDF::NcFile::FileMode::replace );
	netCDF::NcDim R_dim = data_file.addDim( "R", R_pts.size() );
	netCDF::NcDim Z_dim = data_file.addDim( "Z", Z_pts.size() );
	netCDF::NcVar R_var = data_file.addVar( "R", netCDF::NcDouble(), R_dim );
	netCDF::NcVar Z_var = data_file.addVar( "Z", netCDF::NcDouble(), Z_dim );
	R_var.putVar( R_pts.data() );
	Z_var.putVar( Z_pts.data() );

	std::vector<netCDF::NcDim> psi_dims = { R_dim, Z_dim };
	netCDF::NcVar psi = data_file.addVar( "psi", netCDF::NcDouble(), psi_dims );

	BoostArray2D data( boost::extents[ R_pts.size() ][ Z_pts.size() ], boost::c_storage_order() );
	mfem::Mesh *mesh = soln.u_variable.FESpace()->GetMesh();
	int nFound;
	if ( soln.hasUStar )
		nFound = Project( mesh, soln.u_star_variable, R_pts, Z_pts, data );
	else
		nFound = Project( mesh, soln.u_variable, R_pts, Z_pts, data );
	psi.putVar( data.data() );
	data_file.close();
}

