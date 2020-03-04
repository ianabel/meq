
/*
 * Traces curves of a function psi, given gridfunctions for 
 * psi and psi'
 */


#include <eigen3/Eigen/Dense>
#include <functional>
#include <vector>
#include <fstream>
#include <iostream>

#include "mfem.hpp"

using namespace mfem;

class ControlPoint
{
	public:
		double x,y,dx,dy;
		ControlPoint( double a, double b, double da, double db ) : x( a ), y( b ), dx( da ), dy( db ) {};
};

std::vector<ControlPoint> Trace( double x, double y, mfem::GridFunction const& psi, mfem::GridFunction const& gradPsi, double s_max, double delta_s, double tolerance = 1e-6 )
{
	std::vector<ControlPoint> curve;
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

	curve.emplace_back( x, y, gradPsiVal( 0 ), gradPsiVal( 1 ) );

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
		gradPsi.GetVectorValue( elemId[ 0 ], startPointRef[ 0 ], gradPsiVal );
		scaling = ::sqrt( gradPsiVal( 0 )*gradPsiVal( 0 ) + gradPsiVal( 1 )*gradPsiVal( 1 ) );

		curve.emplace_back( current( 0 ), current( 1 ), 0, 0 );
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

	std::string mesh_file = argv[ 1 ];
	// std::string grid_fn_file = argv[ 2 ];
	// std::string vec_grid_fn_file = argv[ 3 ];


	//    Read the mesh from the given mesh file. We can handle triangular,
   //    quadrilateral, tetrahedral, hexahedral, surface and volume meshes with
   //    the same code.
   Mesh *mesh = new Mesh(5, 5, Element::Type::Triangle, false, 1.0, 1.0, true );
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

	FiniteElementCollection *dg_coll   = new DG_FECollection(5, 2);
	FiniteElementSpace *fes = new FiniteElementSpace( mesh, dg_coll );
	FiniteElementSpace *v_fes = new FiniteElementSpace( mesh, dg_coll, 2 );

	GridFunction HimmelblauGrad( v_fes );
	GridFunction Himmelblau( fes );

	FunctionCoefficient HimmelblauCf( HimmelblauFn );
	VectorFunctionCoefficient HimmelblauGradCf( 2, HimmelblauGradFn );

	Himmelblau.ProjectCoefficient( HimmelblauCf );
	HimmelblauGrad.ProjectCoefficient( HimmelblauGradCf );

	std::vector< std::vector<ControlPoint> > Curves;
	double y_0 = 0;
	double y_1 = 5.5;
	double h = 0.15;

	for ( double y=y_0; y < y_1; y += h )
	{
		Curves.emplace_back( Trace( 0.0, y, Himmelblau, HimmelblauGrad, 50, 0.05 ) );
	}
	
	std::ofstream DataFile( "contours.data" );
	for ( auto const &curve : Curves )
	{
		for ( auto const &p : curve )
			DataFile << p.x << "\t" << p.y << std::endl;
		DataFile << std::endl;
		DataFile << "# " << std::endl;
		DataFile << std::endl;
	}

	DataFile.close();
	return 0;

}

