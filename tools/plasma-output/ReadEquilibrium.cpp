
#include <netcdf>
#include "Equilbrium.hpp"


Equilbrium::Equilibrium( std::string const& netcdf_file )
{
	netCDF::NcFile data_file;
	data_file.open( netcdf_file, netCDF::NcFile::FileMode::read );
	netCDF::NcDim R_dim,Z_dim;
	R_dim = data_file.getDim( "R" );
	Z_dim = data_file.getDim( "Z" );

	size_t N_R = R_dim.getSize(), N_Z = Z_dim.getSize();

	R_data.resize( N_R );
	Z_data.resize( N_Z );

	netCDF::NcVar R_var,Z_var;

	R_var = data_file.getVar( "R" );
	Z_var = data_file.getVar( "Z" );

	R_var.getVar( R_data.data() );
	Z_var.getVar( Z_data.data() );

	psiData = new double[ N_R ][ N_Z ];

	netCDF::NcVar psi = data_file.getVar( "Psi" );
	psi.getVar( psiData );
}

double Equilibrium::Psi( double R, double Z )
{
	if ( R < R_data.front || R > R_data.back ) {
		throw std::runtime_error( "R not in range" );
	} else if ( Z < Z_data.front || Z > Z_data.back ) {
		throw std::runtime_error( "Z not in range" );
	}


	auto find_R = [ = ]( double x_l, double x_u ) {
		return ( x_l <= R ) && ( R <= x_u );
	};

	auto find_Z = [ = ]( double x_l, double x_u ) {
		return ( x_l <= Z ) && ( Z <= x_u );
	};

	std::vector<double>::iterator R_interval = std::adjacent_find( R_data.begin(), R_data.end(), find_R );
	std::vector<double>::iterator Z_interval =  std::adjacent_find( Z_data.begin(), Z_data.end(), find_Z );

	ptrdiff_t R_index = std::difference( R_data.begin(), R_interval );
	ptrdiff_t Z_index = std::difference( Z_data.begin(), Z_interval );

	// Target point is in [ *R_interval, *(R_interval + 1) ] x [ *Z_interval, *(Z_interval + 1) ]
	
	// Interpolate in R at Z = *Z_interval;
	double R_u = *( R_interval + 1 );
	double R_l = *( R_interval );
	double psiZ_l = ( ( R_u - R )/( R_u - R_l ) ) * psiData[ R_index ][ Z_index ]
	                 + ( ( R - R_l )/( R_u - R_l ) ) * psiData[ R_index + 1 ][ Z_index ];
	double psiZ_u = ( ( R_u - R )/( R_u - R_l ) ) * psiData[ R_index ][ Z_index + 1 ]
	                 + ( ( R - R_l )/( R_u - R_l ) ) * psiData[ R_index + 1 ][ Z_index + 1 ];

	double Z_u = *( Z_interval + 1 );
	double Z_l = *( Z_interval );
	double result = ( ( Z_u - Z ) / ( Z_u - Z_l ) ) * psiZ_l + ( ( Z - Z_l ) / ( Z_u - Z_l ) ) * psiZ_u;

	return result;
}
