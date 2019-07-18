#include "mfem.hpp"
#include "HDGGSIntegrator.hpp"
#include "FreeBoundary.hpp"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <boost/format.hpp>

using namespace mfem;

class Coil {
	private:
		double R_l,R_u,Z_l,Z_u;
	public:
		double R,z;
		double J;
		double h,w;
		Coil( double R_0, double z_0, double height, double width, double Current )
			: R( R_0 ), z( z_0 ), h( height ), w( width ), J( Current )
		{
			R_l = R - w/2;
			R_u = R + w/2;
			Z_l = z - h/2;
			Z_u = z + h/2;
		};

		bool inline Contains( mfem::Vector const& pt ) const
		{
			return ( ( pt( 0 ) >= R_l ) && ( pt( 0 ) <= R_u ) && ( pt( 1 ) >= Z_l ) && ( pt( 1 ) <= Z_u ) );
		}

		~Coil() {};
};

class Jtor {
	protected:
		std::vector<Coil> Coils;
	public:
		Jtor() {Coils.clear();};
		~Jtor() {};
		void AddCoil( double R, double z, double h, double w, double J ) { Coils.emplace_back( R, z, h, w, J );};

		// Returns R*j_tor
		double operator()( mfem::Vector const& pt ) {
			double JtorPt = 0.0;
			for ( auto const &coil : Coils )
			{
				if ( coil.Contains( pt ) )
					JtorPt += coil.J;
			}
			return JtorPt;
		};
};

int main(int argc, char *argv[])
{
	int order = 6;
	Mesh *mesh = new Mesh(58, 120, Element::Type::TRIANGLE);
	auto xlate = []( const mfem::Vector& in, mfem::Vector & out ) {
		double R_min = 0.05;
		double R_max = 1.5;
		double Z_min = -1.5;
		double Z_max = +1.5;
		out( 0 ) = in( 0 )*( R_max - R_min ) + R_min;
		out( 1 ) = in( 1 )*( Z_max - Z_min ) + Z_min;
		return;
	};
	mesh->Transform( xlate );

	Jtor FieldCoils;
	FieldCoils.AddCoil( 1.,  1., .05, .05, 300. );
	FieldCoils.AddCoil( 1., -1., .05, .05, 300. );
	

	mfem::Array<mfem::Vector> points( 6 );
	points[ 0 ].SetSize( 2 );
	points[ 0 ]( 0 ) = 0.5;
	points[ 0 ]( 1 ) = 0.0;
	points[ 1 ].SetSize( 2 );
	points[ 1 ]( 0 ) = 0.5;
	points[ 1 ]( 1 ) = 0.5;
	points[ 2 ].SetSize( 2 );
	points[ 2 ]( 0 ) = 0.5;
	points[ 2 ]( 1 ) = -0.5;
	points[ 3 ].SetSize( 2 );
	points[ 3 ]( 0 ) = 0.75;
	points[ 3 ]( 1 ) = 0.0;
	points[ 4 ].SetSize( 2 );
	points[ 4 ]( 0 ) = 0.75;
	points[ 4 ]( 1 ) = 0.75;
	points[ 5 ].SetSize( 2 );
	points[ 5 ]( 0 ) = 0.75;
	points[ 5 ]( 1 ) = -0.75;



	std::cout << std::setprecision( 12 );
	/*
	mfem::Vector r( 2 );
	r( 0 ) = 1; r( 1 ) = 1;
	for ( int j=0; j<points.Size(); j++ )
		std::cout << " G at (" << points[ j ]( 0 ) << ", " << points[ j ]( 1 )  << ",1,1) = " << GreensFunction( points[ j ], r ) << std::endl;
		*/

	std::cout << std::endl;

	for ( int j=0; j<points.Size(); j++ )
	{
		std::cout << " Psi at (" << points[ j ]( 0 ) << ", " << points[ j ]( 1 )  << ") = " << GreensFunctionPsi( mesh, points[ j ], FieldCoils ) << std::endl;
	}



	
	return 0;
}


