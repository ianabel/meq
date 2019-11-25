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

	std::shared_ptr<Mesh> mesh = nullptr;
	try {
		mesh = std::make_shared<Mesh>( config->GetFinalMeshFile().c_str() );
	} catch ( std::exception &e ) {
		std::cerr << "Unable to load the mesh from " << config->GetMeshFile() << "." << std::endl;
		std::cerr << " Error: " << e.what() << std::endl << std::endl;
		return -1;
	}


	std::cout << "Using configuration in " << config_file << std::endl;
	meq::Jtor FieldCoils;

	for ( const auto coil : config->GetCoils() )
		FieldCoils.AddCoil( coil );


	for ( auto const &p : outputPoints )
	{
		std::cout << std::setprecision( 11 );

		mfem::Vector ptVec( 2 );
		ptVec( 0 ) = p.first;
		ptVec( 1 ) = p.second;

		std::cout << p.first << "\t" << p.second << "\t" << meq::GreensFunctionPsi( mesh.get(), ptVec, FieldCoils ) << std::endl;
	}

	return 0;
}
