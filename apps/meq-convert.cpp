/*
 * meq-convert: the GLVis and ParaView files of an equilibrium, from its `.nc`.
 *
 *     meq-convert [--gf] [--vtk] [--output <stem>] <run>.nc
 *
 * MEQ writes one file by default, the `.nc`, because it is the whole answer:
 * the solve's mesh as MFEM's own text, every coefficient of psi_h, q and psi*
 * with the space each lives in, and the boundary curve. The other formats are
 * VIEWS of that, and on a large mesh they were most of the run's output time
 * -- the VTK bends the mesh onto Gamma before it draws it -- for files most runs
 * are never opened in. MEASUREMENTS.md M-297.
 *
 * THE FILES ARE THE ONES THE RUN WOULD HAVE WRITTEN. The same library calls
 * make them -- meq::writeMfem() with the embedded mesh text written verbatim,
 * meq::writePostProcessed(), meq::writeEquilibriumVtu() -- from fields read
 * back exactly (meq::StoredEquilibrium), so a run's own `--gf --vtk` files and
 * this program's from that run's `.nc` agree byte for byte, the bent VTK
 * included (OutputConversion asserts it).
 *
 * WHAT IT CANNOT MAKE is a subtracted run's `_psi_total.gf`, which needs the
 * conductor field evaluated on the mesh; `meq --gf` writes that one.
 */

#include <cstdio>
#include <exception>
#include <memory>
#include <string>

#include "mfem.hpp"

#include "meq/Field.hpp"
#include "meq/Output.hpp"
#include "meq/StoredEquilibrium.hpp"

namespace
{
	void usage()
	{
		std::printf(
			"meq-convert -- the GLVis and ParaView files of a MEQ equilibrium\n"
			"\n"
			"  meq-convert [--gf] [--vtk] [--output <stem>] <run>.nc\n"
			"\n"
			"  --gf             <stem>.mesh, _psi.gf, _grad_psi.gf, _psistar.gf\n"
			"  --vtk            <stem>/<name>.pvd, bent onto Gamma as the run\n"
			"                   would have drawn it\n"
			"                   (neither: both)\n"
			"  --output <stem>  the path stem to write at; default the .nc's\n"
			"                   own, beside it\n"
			"\n"
			"The files are those `meq --gf --vtk` writes. A subtracted run's\n"
			"_psi_total.gf needs the conductors and only `meq --gf` writes it.\n"
			"\n"
			"Exit codes: 0 written, 1 usage, 3 the .nc could not be read or a\n"
			"file could not be written, as meq's own exit 3.\n" );
	}
}

int main( int argc, char **argv )
{
	std::string input, stem;
	bool wantGridFunctions = false, wantVtk = false;

	for ( int i = 1; i < argc; ++i )
	{
		std::string const arg = argv[ i ];
		if ( arg == "--help" || arg == "-h" )
		{
			usage();
			return 0;
		}
		else if ( arg == "--gf" )
			wantGridFunctions = true;
		else if ( arg == "--vtk" )
			wantVtk = true;
		else if ( arg == "--output" && i + 1 < argc )
			stem = argv[ ++i ];
		else if ( !arg.empty() && arg[ 0 ] == '-' )
		{
			std::fprintf( stderr, "meq-convert: unknown option \"%s\"\n", arg.c_str() );
			usage();
			return 1;
		}
		else if ( input.empty() )
			input = arg;
		else
		{
			std::fprintf( stderr, "meq-convert: more than one file named (\"%s\" and \"%s\")\n",
			              input.c_str(), arg.c_str() );
			return 1;
		}
	}
	if ( input.empty() )
	{
		usage();
		return 1;
	}
	if ( !wantGridFunctions && !wantVtk )
		wantGridFunctions = wantVtk = true;
	if ( stem.empty() )
		stem = input.size() > 3 && input.compare( input.size() - 3, 3, ".nc" ) == 0
		       ? input.substr( 0, input.size() - 3 ) : input;

	try
	{
		meq::StoredEquilibrium stored( input );
		int const degree = stored.polynomialDegree();

		if ( wantGridFunctions )
		{
			meq::writeMfem( stem, stored.meshText(), stored.potential(), stored.flux() );
			meq::writePostProcessed( stem, stored.postProcessed() );
			std::printf( "meq-convert: wrote %s.mesh, %s_psi.gf (psi_h, degree %d),\n"
			             "             %s_grad_psi.gf and %s_psistar.gf (psi*, degree %d)\n",
			             stem.c_str(), stem.c_str(), degree, stem.c_str(), stem.c_str(),
			             degree + 1 );
			if ( !stored.holdsTotal() )
				std::printf( "meq-convert: these hold psi_p, the REMAINDER, as the run's own\n"
				             "             would; its _psi_total.gf needs `meq --gf`\n" );
		}

		if ( wantVtk )
		{
			// B = ( -q_z, q_r ), the relabelling the run draws; see Field.hpp.
			mfem::GridFunction field( &stored.fluxSpaceOf() );
			meq::poloidalField( stored.flux(), field );

			std::unique_ptr<meq::ElementCellData> cells;
			if ( stored.hasElementColumns() )
				cells = std::make_unique<meq::ElementCellData>(
					stored.mesh(), stored.regionColumn(), stored.plasmaFractionColumn(),
					stored.conductorMaterialColumn(), stored.conductorFractionColumn() );

			meq::writeEquilibriumVtu( stem, stored.mesh(), stored.shape(), degree,
			                          stored.postProcessed(), field, degree + 1,
			                          cells.get() );
			std::string const name = stem.substr( stem.find_last_of( '/' ) + 1 );
			std::printf( "meq-convert: wrote %s/%s.pvd (VTK at degree %d)\n",
			             stem.c_str(), name.c_str(), degree + 1 );
		}
	}
	catch ( std::exception const &error )
	{
		std::fprintf( stderr, "meq-convert: %s\n", error.what() );
		return 3;
	}
	return 0;
}
