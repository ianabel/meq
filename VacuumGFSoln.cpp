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

using namespace mfem;

#include "meq.hpp"

double GreensFunctionSolution( std::vector<meq::Coil> const & CoilSet, std::pair<double,double> pt )
{
	mfem::Vector ptVec( 2 );
	ptVec( 0 ) = pt.first;
	ptVec( 1 ) = pt.second;

	using Integrator = boost::math::quadrature::gauss_kronrod<double, 15>;

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

	std::list< std::pair< double, double> > outputPoints{
		{0.01,0.75},{0.01,0.0},{0.01,-0.75}, {0.5,0.5},{0.5,0.0},{0.5,-0.5}, {0.75,0.75},{0.75,0.0},{0.75,-0.75}
	};


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

	for ( auto const &p : outputPoints )
	{
		std::cout << std::setprecision( 11 );
		std::cout << p.first << "\t" << p.second << "\t" << GreensFunctionSolution( config->GetCoils(), p ) << std::endl;
	}

	return 0;
}
