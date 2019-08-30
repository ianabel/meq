
#include <exception>
#include <iostream>

#include <gmsh.h>

#include "meq.hpp"
#include "MEQConf.hpp"


using namespace meq;


int main( int argc, char** argv )
{
	unsigned int CellsPerCoilEdge = 4;
	std::string configFile;
	if ( argc == 1 )
	{
		configFile = "meq.conf";
	}
	else if ( argc == 2 )
	{
		configFile = argv[ 1 ];
	}
	else
	{
		std::cerr << "Usage is: meshgen <configfile>" << std::endl;
		return 1;
	}

	meq::Configuration *config = nullptr;

	try {
		config = new meq::Configuration( configFile );
	} catch ( toml::syntax_error &tomlErr ) {
		std::cerr << "Error in Configuration File" << std::endl;
		std::cerr << tomlErr.what() << std::endl;
		return 2;
	} catch ( std::exception &other ) {
		std::cerr << "Unknown Error: " << other.what() << std::endl;
		return 3;
	}

	meq::GenerateMesh( config );

	return 0;

}
