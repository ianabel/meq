
#include <netcdf>
#include "Equilbrium.hpp"


Equilbrium::Equilibrium( std::string const& netcdf_file )
{
	netCDF::NcFile data_file;
	data_file.open( netcdf_file, netCDF::NcFile::FileMode::read );
	netCDF::NcDim R_dim,Z_dim;
	R_dim = data_file.getDim( "R" );
	Z_dim = data_file.getDim( "Z" );
}
