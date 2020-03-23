
/*
 * Traces curves of a function psi, given gridfunctions for 
 * psi and psi'
 */


#include <functional>
#include <vector>
#include <fstream>
#include <iostream>
#include <netcdf>

#include "mfem.hpp"

using namespace mfem;

struct Curve {
	std::vector<double> x_pts,y_pts,arclen_data;
};

Curve Trace( double x, double y, mfem::GridFunction const& psi, mfem::GridFunction const& gradPsi, double s_max, double delta_s, double tolerance = 1e-6 )
{
	Curve curve;
	Mesh* mesh = psi.FESpace()->GetMesh();

	mfem::Array<int> elemId( 1 );
	mfem::Array<IntegrationPoint> startPointRef( 1 );
	mfem::Vector current( 2 );
	mfem::DenseMatrix Coords( 2, 1 );

	current( 0 ) = x;
	current( 1 ) = y;
	Coords.UseExternalData( current, 2, 1 );

	int nPts = 0;
	nPts = mesh->FindPoints( Coords, elemId, startPointRef );
	if ( nPts != 1 )
	{
		throw std::logic_error( "Starting point wasn't found" );
	}
	mfem::Vector start( current );

	int currentElement = elemId[ 0 ];
	mfem::IntegrationPoint currentRefPoint = startPointRef[ 0 ];
	mfem::Vector gradPsiVal;

	double PsiVal = psi.GetValue( currentElement, currentRefPoint );
	gradPsi.GetVectorValue( currentElement, currentRefPoint, gradPsiVal );

	curve.x_pts.emplace_back( x );
	curve.y_pts.emplace_back( y );
	curve.arclen_data.emplace_back( 0.0 );

	double arcLength = 0;

	// Step along the curve
	mfem::Vector delta( 2 );
	while ( arcLength < s_max )
	{

		Coords.UseExternalData( current, 2, 1 );
		mesh->FindPoints( Coords, elemId, startPointRef );

		currentElement = elemId[ 0 ];
		currentRefPoint = startPointRef[ 0 ];

		gradPsi.GetVectorValue( currentElement, currentRefPoint, gradPsiVal );

		delta( 0 ) = gradPsiVal( 1 );
		delta( 1 ) = -gradPsiVal( 0 );

		double scaling = delta_s/::sqrt( delta( 0 )*delta( 0 ) + delta( 1 )*delta( 1 ) );

		delta *= scaling;

		mfem::Vector guess = current;
		guess += delta;
		// Check if we're on the level curve, refine if not
		Coords.UseExternalData( guess, 2, 1 );
		nPts = mesh->FindPoints( Coords, elemId, startPointRef );
		if ( nPts == 0 || elemId[ 0 ] == -1 )
		{
			std::cerr << "Terminating because curve left domain" << std::endl;
			return curve;
		}
		double psiGuess = psi.GetValue( elemId[ 0 ],startPointRef[ 0 ] );

		if ( ::fabs( psiGuess - PsiVal ) > tolerance )
		{
			do {
				// Refine with a Newton step.
				gradPsiVal = 0.0;
				gradPsi.GetVectorValue( elemId[ 0 ], startPointRef[ 0 ], gradPsiVal );
				double t = ( PsiVal - psiGuess )/( gradPsiVal( 0 )*gradPsiVal( 0 ) + gradPsiVal( 1 )*gradPsiVal( 1 ) );
				gradPsiVal *= t;
				guess += gradPsiVal;
				Coords.UseExternalData( guess, 2, 1 );
				nPts = mesh->FindPoints( Coords, elemId, startPointRef );
				if ( ( nPts == 0 ) || ( elemId[ 0 ] == -1 ) )
				{
					std::cerr << "Terminating because curve left domain" << std::endl;
					return curve;
				}
				psiGuess = psi.GetValue( elemId[ 0 ],startPointRef[ 0 ] );
			} while ( ::fabs( psiGuess - PsiVal ) > tolerance );
		}


		mfem::Vector tmp = guess;
		tmp -= current;
		arcLength += ::sqrt( tmp( 0 )*tmp( 0 ) + tmp( 1 )*tmp( 1 ) );
		current = guess;

		curve.x_pts.emplace_back( current( 0 ) );
		curve.y_pts.emplace_back( current( 1 ) );
		curve.arclen_data.emplace_back( arcLength );
		tmp = current; 
		tmp -= start;
		if ( ( tmp( 0 )*tmp( 0 ) + tmp( 1 )*tmp( 1 ) ) < .01*delta_s*delta_s )
		{
			std::cerr << "Stopping because curve is closed" << std::endl;
			break;
		}
	}
	return curve;
}

double HimmelblauFn( Vector const &r ) {
	double x=r( 0 ),y=r( 1 );
	return ( x*x + y - 11. )*( x*x + y - 11. ) + ( x + y*y - 7. )*( x + y*y - 7. );
};

	// {4 x (-11 + x^2 + y) + 2 (-7 + x + y^2), 2 (-11 + x^2 + y) + 4 y (-7 + x + y^2)}
void HimmelblauGradFn( Vector const &r, Vector& g ) {
	double x=r( 0 ),y=r( 1 );
	g( 0 ) = 4*x*( -11.0 + x*x + y ) + 2.0*( -7.0 + x + y*y );
	g( 1 ) = 2*(-11.0 + x*x + y) + 4*y*(-7.0 + x + y*y);
};

int main( int argc, char** argv )
{

	// std::string grid_fn_file = argv[ 2 ];
	// std::string vec_grid_fn_file = argv[ 3 ];


	//    Read the mesh from the given mesh file. We can handle triangular,
   //    quadrilateral, tetrahedral, hexahedral, surface and volume meshes with
   //    the same code.
   Mesh *mesh = new Mesh(5, 5, Element::Type::TRIANGLE, false, 1.0, 1.0, true );
	auto xform = []( const Vector& in, Vector& out ) { 
		constexpr double R_min = -6.0;
		constexpr double R_max = 6.0;
		constexpr double Z_min = -6.0;
		constexpr double Z_max = 6.0;
		out( 0 ) = R_min + in( 0 )*( R_max - R_min );
		out( 1 ) = Z_min + in( 1 )*( Z_max - Z_min );
	};
	mesh->Transform( xform );
	mesh->Finalize(true, true);

	/*
	//    Read the grid function
	std::fstream grid_fs( grid_fn_file.c_str() );
	GridFunction psi( mesh, grid_fs );
	grid_fs.close();
	std::fstream vgrid_fs( vec_grid_fn_file.c_str() );
	GridFunction gradPsi( mesh, vgrid_fs );
	vgrid_fs.close();
	*/

	int ord = 4;
	FiniteElementCollection *dg_coll   = new DG_FECollection(ord, 2);
	FiniteElementSpace *fes = new FiniteElementSpace( mesh, dg_coll );
	FiniteElementSpace *v_fes = new FiniteElementSpace( mesh, dg_coll, 2 );

	GridFunction HimmelblauGrad( v_fes );
	GridFunction Himmelblau( fes );

	FunctionCoefficient HimmelblauCf( HimmelblauFn );
	VectorFunctionCoefficient HimmelblauGradCf( 2, HimmelblauGradFn );

	Himmelblau.ProjectCoefficient( HimmelblauCf );
	HimmelblauGrad.ProjectCoefficient( HimmelblauGradCf );

	double y_0 = 0;
	double y_1 = 5.5;
	double h = 0.15;

	std::vector<double> initial_y_values;
	for ( double y=y_0; y < y_1; y += h )
		initial_y_values.emplace_back( y );

	unsigned int N_surfaces = 
	
	netCDF::NcFile data_file( "contours.nc", netCDF::NcFile::FileMode::replace );

	netCDF::NcDim psi_dim = data_file.addDim( "Psi", psi_values.size() );
	netCDF::NcVar psi_var = data_file.addVar( "Psi", netCDF::ncDouble, psi_dim );

	// R(psi_val) , Z(psi_val), s(psi_val) are all VLENs of NC_DOUBLEs,
	// define that type
	
	netCDF::NcVlenType double_vector = data_file.addVlenType( "DoubleVector", netCDF::ncDouble );

	netCDF::NcVar R_var = data_file.addVar( "R", double_vector, psi_dim );
	netCDF::NcVar Z_var = data_file.addVar( "Z", double_vector, psi_dim );
	netCDF::NcVar s_var = data_file.addVar( "s", double_vector, psi_dim );

	std::vector< Curve > Curves( psi_values.size() );
	for ( size_t j = 0; j < psi_values.size(); j++ )
	{
		Curves[ j ] = Trace( 0.0, psi_values[ j ], Himmelblau, HimmelblauGrad, 50, 0.05 );

		psi_var.putVar( psi_values.data() );
	}
	
	std::ofstream DataFile( "contours.data" );
	for ( size_t j = 0; j < Curves.size(); j++ )
	{
		Curve& curve = Curves[ j ];

		size_t n_pts = curve.x_pts.size();
		if ( curve.y_pts.size() != n_pts )
			throw std::logic_error( "VLENs should be the same length" );
		nc_vlen_t r_data{n_pts,curve.x_pts.data()},z_data{n_pts,curve.y_pts.data()},s_data{n_pts,curve.arclen_data.data()};
		std::vector<size_t> index{j};

		R_var.putVar( index, &r_data );
		Z_var.putVar( index, &z_data );
		s_var.putVar( index, &s_data );

		for ( size_t i=0; i < n_pts; i++ )
		{
			DataFile << curve.x_pts[ i ] << "\t" << curve.y_pts[ i ] << std::endl;
		}
		DataFile << std::endl;
		DataFile << "# " << std::endl;
		DataFile << std::endl;
	}

	data_file.close();

	DataFile.close();
	return 0;

}

