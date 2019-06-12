
/*
 * An MFEM-based nonlinear poisson solver, derived from the HDG Example Poisson Solver
 * of T. Horvath, S. Rhebergen, and A. Sivas.
 *
 * The discretization is based on the paper:
 * N.C. Nguyen, J. Peraire, B. Cockburn, An implicit high-order hybridizable
 * discontinuous Galerkin method for linear convection–diffusion equations,
 * J. Comput. Phys., 2009, 228:9, 3232--3254.
 *
 * Modifications are by I. G. Abel, University of Maryland
 */

#include "mfem.hpp"
#include "GSInverter.hpp"

#include <fstream>
#include <iostream>
#include <algorithm>
#include <functional>


using namespace std;
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

	// Parse command-line options.
	const char *mesh_file = "grid.mesh";
	int order = 3;
	int initial_ref_levels = 0;
	bool visualization = true;
	bool post = true;
	bool save = true;

	// TODO -- replace with boost option parser
	OptionsParser args(argc, argv);
	args.AddOption(&mesh_file, "-m", "--mesh",
			"Mesh file to use.");
	args.AddOption(&order, "-o", "--order",
			"Finite element order (polynomial degree).");
	args.AddOption(&visualization, "-vis", "--visualization", "-no-vis",
			"--no-visualization",
			"Enable or disable GLVis visualization.");
	args.AddOption(&post, "-post", "--postprocessing",
			"-no-post", "--no-postprocessing",
			"Enable or disable postprocessing.");
	args.AddOption(&save, "-save", "--save-files", "-no-save",
			"--no-save-files",
			"Enable or disable file saving.");
	args.AddOption(&initial_ref_levels, "-mr", "--mesh-refinement-levels",
			"The number of levels of uniform refinement to apply to the grid.");
	args.Parse();
	if (!args.Good())
	{
		args.PrintUsage(cout);
		return 1;
	}
	args.PrintOptions(cout);

	
	
	// Do the vacuum problem
	// One coil at R=1,Z=0, h=w=0.05, J = 400 kA/m^2 => j_tot = 1kA
	
	Jtor FieldCoils;
	FieldCoils.AddCoil( 1., 0., .05, .05, 400. );
	
	// Mesh goes from R_min to R_max, and Z_min to Z_max
	
	Mesh *mesh = new Mesh(10, 10, Element::Type::TRIANGLE);
	auto xlate = []( const mfem::Vector& in, mfem::Vector & out ) {
		double R_min = 0.05;
		double R_max = 1.5;
		double Z_min = -1;
		double Z_max = +1;
		out( 0 ) = in( 0 )*( R_max - R_min ) + R_min;
		out( 1 ) = in( 1 )*( Z_max - Z_min ) + Z_min;
		return;
	};
	mesh->Transform( xlate );

	int dim = mesh->Dimension();

	for (int ii=0; ii<initial_ref_levels; ii++)
	{
		mesh->UniformRefinement();
	}

	// Thus the RHS function is FieldCoils::operator()
	
	ConstantCoefficient zero( 0.0 );
	GSInverter solver( mesh, order, FieldCoils );
	solver.SetBCs( zero );
	
	mfem::Vector qu_zero_bc;

	solver.Mult( qu_zero_bc, qu_zero_bc );


	GridFunction q_solution,u_solution;

	q_solution.MakeRef( solver.GetQSpace(), qu_zero_bc, 0 );
	u_solution.MakeRef( solver.GetUSpace(), qu_zero_bc, solver.GetQSpace()->GetVSize() );



	if (save)
	{
		ofstream mesh_ofs("refined_vacuum.mesh");
		mesh_ofs.precision(8);
		mesh->Print(mesh_ofs);

		ofstream q_solution_ofs("sol_q.gf");
		q_solution_ofs.precision(8);
		q_solution.Save(q_solution_ofs);

		ofstream u_solution_ofs("sol_u.gf");
		u_solution_ofs.precision(8);
		u_solution.Save(u_solution_ofs);

	}

	// 14. Send the solution by socket to a GLVis server.
	if (visualization)
	{
		char vishost[] = "localhost";
		int  visport   = 19916;
		socketstream u_sock(vishost, visport);
		u_sock.precision(8);
		u_sock << "solution\n" << *mesh << u_solution << "window_title 'Solution u'" <<
			endl;

		socketstream q_sock(vishost, visport);
		q_sock.precision(8);
		q_sock << "solution\n" << *mesh << q_solution << "window_title 'Solution q'" <<
			endl;
	}

	// 18. Free the used memory.
	delete mesh;

	return 0;

}

