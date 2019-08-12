#include "toml11/toml.hpp"
#include <vector>
#include <iostream>

int main( int argc, char** argv )
{
	const auto config = toml::parse( "meq.conf" );
	const auto coils = toml::find< std::vector< toml::table > >( config, "coils" );
	const auto options = toml::find< toml::table >( config, "options" );
	


	for ( const auto coil : coils )
	{
		double R = coil.at( "R" ).as_floating();
		double Z = coil.at( "Z" ).as_floating();
		double w = coil.at( "Width" ).as_floating();
		double h = coil.at( "Height" ).as_floating();
		std::cout << "Coil at (" << R << ", " << Z << ") with width " << w << " and height " << h << std::endl;
	}


}
