/*
 * Computes the exact Vacuum solution by integrating the greens function.
 * takes coilset from meq config files.
 */
#include "mfem.hpp"
#include <fstream>
#include <iostream>
#include <vector>
#include <utility>
#include <list>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <boost/multi_array.hpp>
#include <netcdf>

using namespace mfem;

#include "meq.hpp"

double GreensFunctionSolution( std::vector<meq::Coil> const & CoilSet, std::pair<double,double> pt )
{
	mfem::Vector ptVec( 2 );
	ptVec( 0 ) = pt.first;
	ptVec( 1 ) = pt.second;

	using Integrator = boost::math::quadrature::gauss_kronrod<double, 31>;

	double result = 0;
	const static double max_h = .05;
	for ( auto const & coil : CoilSet )
	{
		double partial_integral = 0;
		/*
		 * We have to integrate over the rectangle
		 *  [ R_min, R_max ] x [ Z_min, Z_max ]
		 */
	
		/*
		 * First split into multiple rectangles:
		 */	
		int N_z = ( coil.h / max_h + .5 );
		int N_R = ( coil.w / max_h + .5 );

		double delta_z = coil.h / N_z;
		double delta_R = coil.w / N_R;

		double R_min = coil.R - coil.w / 2;
		double Z_min = coil.z - coil.h / 2;

		auto GFintegrand = [ & ]( double R, double z ){
			mfem::Vector ptStarVec( 2 );
			ptStarVec( 0 ) = R;
			ptStarVec( 1 ) = z;
			return meq::GreensFunction( ptVec, ptStarVec )/R;
		};

		auto IntegrateZ = [&]( double R, double Z_min, double Z_max ){
			auto integrand = std::bind( GFintegrand, R, std::placeholders::_1 );
			return Integrator::integrate( integrand, Z_min, Z_max );
		};

		for ( int j=0; j < N_z; j++ )
			for ( int i=0; i < N_R; i++ )
			{
				// Integrate over rectangle [ R_min + i*delta_R, R_min + (i+1)*delta_R ] x [ Z_min + j*delta_z, Z_min + (j+1)*delta_z ]
				auto IntegratedInZ = std::bind( IntegrateZ, std::placeholders::_1, Z_min + j*delta_z, Z_min + ( j+1 )*delta_z );
				partial_integral += Integrator::integrate( IntegratedInZ, R_min + i*delta_R, R_min + ( i+1 )*delta_R );
			}

		result += partial_integral * coil.J;
	}
	return result;
}


int main(int argc, char *argv[])
{
	std::string config_file;

	if ( argc == 1 )
		config_file = "meq.conf";
	else if ( argc == 2 )
		config_file = argv[ 1 ];
	else
	{
		std::cerr << "./vacuum-test <meq config file>" << std::endl;
		return -1;
	}

	std::shared_ptr<meq::Configuration> config = nullptr;

	try {
		config = std::make_shared<meq::Configuration>( config_file );
	} catch ( toml::syntax_error &tomlErr ) {
		
		std::cerr << "Error parsing the configuration in " << config_file << std::endl;
		std::cerr << tomlErr.what() << std::endl;
		return 2;
	} catch ( std::exception &other ) {
		std::cerr << "Unknown Error in loading configuration: " << other.what() << std::endl;
		return 3;
	}

	std::cout << "Using configuration in " << config_file << std::endl;

	meq::Domain const & domain = *( config->GetDomain() );

	const int N_R( 51 );
	const int N_Z( 51 );
	double R_min = 0.1;
	double R_max = 0.9;
	double Z_min = -1.0;
	double Z_max = 1.0;
	std::vector<double> R_pts( N_R );
	std::vector<double> Z_pts( N_Z );
	for ( unsigned int i=0; i < N_R; i++ )
		R_pts[ i ] = R_min + i*( R_max - R_min )/( N_R - 1 );
	for ( unsigned int i=0; i < N_Z; i++ )
		Z_pts[ i ] = Z_min + i*( Z_max - Z_min )/( N_Z - 1 );


using BoostArray2D = boost::multi_array<double, 2>;

	BoostArray2D data( boost::extents[ R_pts.size() ][ Z_pts.size() ], boost::c_storage_order() );

	for ( unsigned int i = 0; i < R_pts.size(); i++ )
		for ( unsigned int j = 0; j < Z_pts.size(); j++ )
			data[ i ][ j ] = GreensFunctionSolution( config->GetCoils(), {R_pts[ i ],Z_pts[ j ]} );

	{
		netCDF::NcFile data_file( "exact-vacuum.nc", netCDF::NcFile::FileMode::replace );
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

	return 0;
}
