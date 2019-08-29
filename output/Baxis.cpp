#include "mfem.hpp"
#include "toml11/toml.hpp"
#include <fstream>
#include <iostream>
#include <vector>
#include <utility>
#include <list>

using namespace mfem;

bool inline inElement( IntegrationPoint const& ip ) { return ( ip.x > 0 ) && ( ip.x < 1 ) && ( ip.y > 0 ) && ( ip.y < 1 ) && (  ip.x + ip.y < 1 ); };

std::ostream& operator<<( std::ostream& os, Vector const& vec )
{
	os << "(" << vec( 0 ) << "," << vec( 1 ) << ")";
	return os;
}

int main(int argc, char *argv[])
{

	const char* config_file = "meq.conf";

	if ( argc > 1 )
		config_file = argv[ 1 ];

	const auto config = toml::parse( config_file );
	const auto options = toml::find< toml::table >( config, "options" );
	const auto domain = toml::find< toml::table >( config, "domain" );

	std::string MeshFile;
	MeshFile = options.at( "MeshFile" ).as_string();

	double R_axis = domain.at( "RMin" ).as_floating();
	double ZMin = domain.at( "ZMin" ).as_floating();
	double ZMax = domain.at( "ZMax" ).as_floating();


	std::list< std::pair< double, double> > outputPoints;

	unsigned int NPoints = 250;
	double dz = ( ZMax - ZMin )/( NPoints + 1 );

	for ( unsigned int i = 0; i < NPoints; i++ )
			outputPoints.emplace_back( R_axis, ( i + 1 ) * dz + ZMin );

   //    Read the mesh from the given mesh file. We can handle triangular,
   //    quadrilateral, tetrahedral, hexahedral, surface and volume meshes with
   //    the same code.
   Mesh *mesh = new Mesh(MeshFile.c_str(), 1, 1);

	std::string grid_fn_file;
	grid_fn_file = options.at( "GradPsiFile" ).as_string();
	//    Read the grid function
	std::fstream grid_fs( grid_fn_file.c_str() );
	GridFunction grid_func( mesh, grid_fs );
	grid_fs.close();

	FiniteElementSpace *fes = grid_func.FESpace();

	unsigned int nElems = fes->GetNE();
	unsigned int vDim = 2;

	for ( unsigned int iElem = 0; iElem < nElems; iElem++ )
	{
		ElementTransformation *Tr = mesh->GetElementTransformation( iElem );
		IntegrationRule ir( 0 );
		std::vector<std::pair<double,double>> pointsInElem;


		Vector TestVec( 2 );
		for ( std::pair<double, double> const &x : outputPoints )
		{
			IntegrationPoint ip;
			TestVec( 0 ) = x.first; TestVec( 1 ) = x.second;
			int res = Tr->TransformBack( TestVec, ip );
			if ( res == 0 )
			{
				ir.Append( ip );
				pointsInElem.push_back( x );
			}
		}


		unsigned int nPts = ir.GetNPoints();
		if ( nPts != pointsInElem.size() )
		{
			std::cerr << "Bork" << std::endl; return -1;
		}

		DenseMatrix values;
		grid_func.GetVectorValues( *Tr, ir, values );
		for ( unsigned int j = 0 ; j < nPts; j++ )
		{
			std::cout << pointsInElem[ j ].first << "\t" << pointsInElem[ j ].second << "\t";
			for ( unsigned int i = 0; i < vDim; i++ )
				std::cout << values( i, j ) << "\t";
			std::cout << std::endl;
		}

	}

	delete mesh;
	return 0;
}
