
#include <exception>
#include <iostream>

#include <gmsh.h>

#include "meq.hpp"
#include "MEQConf.hpp"


namespace meq {

void GenerateMesh( std::shared_ptr<meq::Configuration> config )
{
	unsigned int CellsPerCoilEdge = 4;
	std::shared_ptr<meq::Domain> domain = config->GetDomain();

	char opts[] = "-format vtk";
	char *gmsh_args[ 1 ];
	gmsh_args[ 0 ] = opts;
	gmsh::initialize( 1, gmsh_args );

	namespace geo = gmsh::model::geo;

	double RMin = domain->RMin;
	double RMax = domain->RMax;
	double ZMin = domain->ZMin;
	double ZMax = domain->ZMax;
	double h = domain->CellSize;

	int LowerLeft,LowerRight,UpperLeft,UpperRight;

	LowerLeft  = geo::addPoint( RMin, ZMin, 0, h );
	LowerRight = geo::addPoint( RMax, ZMin, 0, h );
	UpperLeft  = geo::addPoint( RMin, ZMax, 0, h );
	UpperRight = geo::addPoint( RMax, ZMax, 0, h );

	std::vector<int> domainCurve;
	domainCurve.emplace_back( geo::addLine( LowerLeft,  LowerRight ) );
	domainCurve.emplace_back( geo::addLine( LowerRight, UpperRight ) );
	domainCurve.emplace_back( geo::addLine( UpperRight, UpperLeft  ) );
	domainCurve.emplace_back( geo::addLine( UpperLeft,  LowerLeft  ) );
	int domainCurveTag = geo::addCurveLoop( domainCurve );

	std::vector<int> coilCurveTags = {domainCurveTag}; // Initialise with the Domain Curve so appending coil curves carves out holes
	std::vector<int> coilSurfaceTags;


	for ( const auto &coil : config->GetCoils() )
	{
		double R = coil.R;
		double Z = coil.z;
		double wCoil = coil.w;
		double hCoil = coil.h;
		RMin = R - wCoil/2.0;
		RMax = R + wCoil/2.0;
		ZMin = Z - hCoil/2.0;
		ZMax = Z + hCoil/2.0;
		h = std::min( hCoil, wCoil )/CellsPerCoilEdge;
		LowerLeft  = geo::addPoint( RMin, ZMin, 0, h );
		LowerRight = geo::addPoint( RMax, ZMin, 0, h );
		UpperLeft  = geo::addPoint( RMin, ZMax, 0, h );
		UpperRight = geo::addPoint( RMax, ZMax, 0, h );

		std::vector<int> coilCurve;
		coilCurve.clear();
		coilCurve.emplace_back( geo::addLine( LowerLeft,  LowerRight ) );
		coilCurve.emplace_back( geo::addLine( LowerRight, UpperRight ) );
		coilCurve.emplace_back( geo::addLine( UpperRight, UpperLeft  ) );
		coilCurve.emplace_back( geo::addLine( UpperLeft,  LowerLeft  ) );
		int coilCurveTag = geo::addCurveLoop( coilCurve );

		coilCurveTags.emplace_back( coilCurveTag );

		std::vector<int> ccT = {coilCurveTag};
		coilSurfaceTags.emplace_back( geo::addPlaneSurface( ccT ) );
		geo::synchronize();
	};

	int domainSurfaceTag = geo::addPlaneSurface( coilCurveTags );

	coilSurfaceTags.emplace_back( domainSurfaceTag );

	gmsh::model::addPhysicalGroup( 2, coilSurfaceTags, 1 );
	gmsh::model::mesh::generate( 2 );

	gmsh::write( config->GetMeshFile() );
	gmsh::finalize();
}

}
