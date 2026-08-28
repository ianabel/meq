/*
 * meq -- the Maryland Equilibrium Solver.
 *
 *     meq config.toml
 *
 * One binary, no subcommands, and deliberately NOT mfem::OptionsParser, which
 * the program this replaces used and which wants to own argument parsing for
 * the whole executable. See DRIVER-PLAN.md section 5.
 *
 * EXIT CODES, because a shell script driving a parameter scan is a first-class
 * caller and "it printed something" is not an interface:
 *
 *     0  solved, and everything asked for was written
 *     1  the configuration is wrong, or the run cannot be set up from it
 *     2  the non-linear solve did not converge
 *     3  output could not be written
 *
 * Exit 2 is only reportable because ../mfem/install is built with
 * MFEM_USE_EXCEPTIONS: without it a failed solve takes the process down with
 * SIGABRT before anything here runs. See CLAUDE.md under "On SUNDIALS".
 */

#include "mfem.hpp"

#include "meq_version.hpp"

#include "meq/BoundaryShape.hpp"
#include "meq/Config.hpp"
#include "meq/Field.hpp"
#include "meq/GradShafranov.hpp"
#include "meq/Output.hpp"
#include "meq/Sampler.hpp"
#include "meq/SourceFactory.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <limits>
#include <memory>
#include <string>
#include <vector>

namespace
{
	enum ExitCode
	{
		Solved = 0,
		ConfigurationError = 1,
		SolveFailed = 2,
		OutputFailed = 3
	};

	void usage()
	{
		std::printf(
			"meq -- the Maryland Equilibrium Solver\n"
			"\n"
			"  meq <config.toml>    solve the equilibrium the file describes\n"
			"  meq --help           this\n"
			"  meq --version        the build this is\n"
			"\n"
			"Exit codes: 0 solved, 1 configuration, 2 solve did not converge,\n"
			"3 output could not be written.\n" );
	}

	/*
	 * The one trap in CLAUDE.md a released binary can protect its user from.
	 *
	 * libblas.so.3 on the development machine resolves to libmkl_rt.so, which
	 * without MKL_THREADING_LAYER=GNU silently corrupts UMFPACK's BLAS-3 -- you
	 * get numbers, and they are wrong. CMake sets it on every registered ctest;
	 * a user running this binary by hand gets no such help.
	 *
	 * Refusing outright would be wrong, since most builds do not link MKL at
	 * all and there is no portable way to ask the loader what libblas resolved
	 * to. So this warns, loudly, on stderr, and says what to do. A wrong answer
	 * the user was warned about beats a wrong answer in silence.
	 */
	void warnAboutThreadingLayer()
	{
		if ( std::getenv( "MKL_THREADING_LAYER" ) != nullptr )
			return;

		std::fprintf( stderr,
			"meq: warning: MKL_THREADING_LAYER is not set.\n"
			"     If this build's BLAS resolves to MKL, UMFPACK's BLAS-3 is\n"
			"     silently wrong without it -- you get numbers, and they are not\n"
			"     the answer. Re-run as:\n"
			"\n"
			"         MKL_THREADING_LAYER=GNU meq <config.toml>\n"
			"\n" );
	}

	/// The background mesh: the box from [mesh], or a file, then refinement.
	mfem::Mesh buildMesh( meq::MeshConfig const &config )
	{
		mfem::Mesh mesh = config.fromFile()
			? mfem::Mesh( config.file.c_str(), 1, 1 )
			: mfem::Mesh::MakeCartesian2D( config.nR, config.nZ,
			                               mfem::Element::TRIANGLE, false,
			                               config.rMax - config.rMin,
			                               config.zMax - config.zMin );

		// MakeCartesian2D puts the box at the origin; move it to [RMin,RMax] x
		// [ZMin,ZMax]. A mesh read from a file is already where it is.
		if ( !config.fromFile() )
		{
			mfem::Vector shift( 2 );
			shift( 0 ) = config.rMin;
			shift( 1 ) = config.zMin;
			mesh.Transform( [ shift ]( mfem::Vector const &in, mfem::Vector &out )
			{
				out = in;
				out( 0 ) += shift( 0 );
				out( 1 ) += shift( 1 );
			} );
		}

		for ( int i = 0; i < config.refinementLevels; ++i )
			mesh.UniformRefinement();

		return mesh;
	}

	/*
	 * The computational subdomain D_h and the transfer path, for a curved Gamma.
	 *
	 * This is GS-2's technique: Gamma is a level set of the shape rather than a
	 * union of mesh faces, so the mesh is NOT fitted to it. D_h is the union of
	 * background elements lying inside Gamma, and the Dirichlet datum is
	 * transferred from Gamma onto Gamma_h = boundary( D_h ) along short paths.
	 * The point of it is that a polygonal approximation to a curved boundary
	 * caps the convergence rate whatever the polynomial degree -- CLAUDE.md
	 * measures 2.12 at k = 3 on a fixed 40-gon against 4.00 for the same
	 * discretisation on a fitted domain.
	 *
	 * The sequence follows tests/convergence/ExtensionConvergence.cpp, which is
	 * where it is measured, and miniapps/hdg/extension.cpp in the MFEM tree.
	 */
	struct Subdomain
	{
		std::unique_ptr<mfem::SubMesh> mesh;
		std::unique_ptr<mfem::VertexConePath> path;
		mfem::Array<int> gammaHMarker;
		int gammaH = 0;
		int widened = 0;
	};

	/// The characteristic background cell size, which sets how far the transfer
	/// paths are allowed to search.
	double backgroundCellSize( meq::MeshConfig const &config )
	{
		double const levels = static_cast<double>( 1 << config.refinementLevels );
		double const hR = ( config.rMax - config.rMin )/( config.nR*levels );
		double const hZ = ( config.zMax - config.zMin )/( config.nZ*levels );
		return std::max( hR, hZ );
	}

	Subdomain buildSubdomain( mfem::Mesh &background,
	                          meq::BoundaryShape const &shape,
	                          double h )
	{
		meq::BoundaryShape const *shapePointer = &shape;
		mfem::PositionFunction const levelSet =
			[ shapePointer ]( mfem::Vector const &x )
			{
				return shapePointer->levelSet( x( 0 ), x( 1 ) );
			};

		// extra_refine = 1: the vertex test alone is exact only where Omega is
		// convex, and a flux surface with triangularity is not obviously so.
		mfem::Array<int> marker;
		int const inside = mfem::MarkLevelSetSubdomain( background, levelSet, 0.0,
		                                               marker, 1 );
		if ( inside == 0 )
			throw std::runtime_error(
				"[boundary.shape] encloses no background element. The mesh is too "
				"coarse for the surface, or the surface lies outside the [mesh] box" );

		for ( int e = 0; e < background.GetNE(); ++e )
			background.SetAttribute( e, marker[ e ] ? 1 : 2 );
		background.SetAttributes();

		mfem::Array<int> domainAttribute( 1 );
		domainAttribute[ 0 ] = 1;

		Subdomain subdomain;
		subdomain.mesh = std::make_unique<mfem::SubMesh>(
			mfem::SubMesh::CreateFromDomain( background, domainAttribute ) );

		// SubMesh gives the boundary it had to generate ONE new attribute, and
		// leaves anything inherited from the background with the attribute it
		// already had. So more than one attribute means D_h reaches the edge of
		// the box: part of Gamma_h would then be a fitted mesh boundary carrying
		// no transferred datum, and the solve would silently impose zero there.
		if ( subdomain.mesh->bdr_attributes.Size() != 1 )
			throw std::runtime_error(
				"[boundary.shape] touches the edge of the [mesh] box: the subdomain "
				"has boundary inherited from it, so part of Gamma_h is fitted and "
				"carries no transferred datum. Enlarge the box" );

		subdomain.gammaH = subdomain.mesh->bdr_attributes.Max();

		// Six h of search. The paths are about 1.3 h long, so this is a factor of
		// four of slack. VertexConePath is the family Cockburn and Solano analyse,
		// and it widens rather than failing when a ray finds no root -- which is
		// why the count is reported rather than ignored.
		subdomain.path = std::make_unique<mfem::VertexConePath>(
			*subdomain.mesh, subdomain.gammaH, levelSet, 6.0*h );
		subdomain.widened = subdomain.path->NumWidened();

		subdomain.gammaHMarker.SetSize( subdomain.gammaH );
		subdomain.gammaHMarker = 0;
		subdomain.gammaHMarker[ subdomain.gammaH - 1 ] = 1;
		return subdomain;
	}

	/// The residual history, in the shape CLAUDE.md records. It is the
	/// diagnostic that separates a wrong Jacobian from a hard problem -- a run
	/// that grinds down linearly means the two disagree -- and it costs nothing.
	void reportResiduals( std::vector<double> const &residuals )
	{
		if ( residuals.empty() )
			return;

		std::printf( "\n   it            ||r||    ||r||/||r_0||    order\n" );
		double const first = residuals.front();
		for ( std::size_t i = 0; i < residuals.size(); ++i )
		{
			std::printf( "  %3zu     %12.6e     %12.6e", i, residuals[ i ],
			             first > 0.0 ? residuals[ i ]/first : 0.0 );

			// Observed order needs three points, and is meaningless once the
			// residual is at the round-off floor.
			if ( i >= 2 && residuals[ i ] > 0.0 && residuals[ i-1 ] > 0.0
			     && residuals[ i-2 ] > 0.0 && residuals[ i ] > 1.0e-14 )
			{
				double const a = std::log( residuals[ i ]/residuals[ i-1 ] );
				double const b = std::log( residuals[ i-1 ]/residuals[ i-2 ] );
				if ( std::abs( b ) > 1.0e-14 )
					std::printf( "    %5.3f", a/b );
				else
					std::printf( "        -" );
			}
			else
			{
				std::printf( "        -" );
			}
			std::printf( "\n" );
		}
		std::printf( "\n" );
	}
}

int main( int argc, char **argv )
{
	std::string const argument = argc > 1 ? argv[ 1 ] : std::string();

	if ( argc != 2 || argument == "--help" || argument == "-h" )
	{
		usage();
		return argument == "--help" || argument == "-h"
			? Solved : ConfigurationError;
	}

	if ( argument == "--version" )
	{
		std::printf( "meq %s (MFEM %s)\n", MEQ_VERSION, MFEM_VERSION_STRING );
		return Solved;
	}

	warnAboutThreadingLayer();

	// ---- configuration -------------------------------------------------
	std::unique_ptr<meq::Configuration> config;
	std::shared_ptr<meq::Source const> source;
	try
	{
		config = std::make_unique<meq::Configuration>( argument );
		source = meq::makeSource( config->getSource(), argument );

		if ( config->getBoundary().type == meq::BoundaryDataType::Exact )
		{
			std::fprintf( stderr,
				"meq: [boundary] Type = \"exact\" needs the source's closed-form\n"
				"     solution, which meq::Source does not carry -- it is a\n"
				"     convergence-study device and lives in the test fixtures.\n"
				"     Use Type = \"zero\" for the fixed-boundary problem.\n" );
			return ConfigurationError;
		}

		if ( config->getAdaptivity().enabled )
		{
			std::fprintf( stderr,
				"meq: [adaptivity] Enabled = true, and the driver does not yet run\n"
				"     the loop. Refusing rather than quietly solving once.\n" );
			return ConfigurationError;
		}
	}
	catch ( std::exception const &error )
	{
		std::fprintf( stderr, "meq: %s\n", error.what() );
		return ConfigurationError;
	}

	// ---- set the run up ------------------------------------------------
	std::unique_ptr<meq::GradShafranovSolver> solver;
	mfem::Mesh background;
	std::unique_ptr<meq::BoundaryShape> shape;
	Subdomain subdomain;
	mfem::ConstantCoefficient zero( 0.0 );
	std::unique_ptr<mfem::FunctionCoefficient> ramp;
	std::unique_ptr<mfem::Mesh> guessMesh;
	std::unique_ptr<mfem::GridFunction> guess;
	// The mesh actually solved on: D_h on the curved path, the background mesh
	// on the fitted one. Everything downstream uses it.
	mfem::Mesh *solveMesh = nullptr;
	try
	{
		background = buildMesh( config->getMesh() );

		// The curved path solves on D_h, a SUBSET of the background mesh; the
		// fitted path solves on the background mesh itself. Everything
		// downstream -- the solver, the sampler, the files -- follows this.
		meq::ShapeConfig const &shapeConfig = config->getBoundary().shape;
		if ( shapeConfig.type != meq::ShapeType::None )
		{
			shape = std::make_unique<meq::BoundaryShape>(
				shapeConfig.type == meq::ShapeType::Miller
					? meq::BoundaryShape::miller( shapeConfig.majorRadius,
					                              shapeConfig.centreHeight,
					                              shapeConfig.minorRadius,
					                              shapeConfig.elongation,
					                              shapeConfig.triangularity,
					                              shapeConfig.squareness )
					: meq::BoundaryShape( shapeConfig.majorRadius,
					                      shapeConfig.centreHeight,
					                      shapeConfig.minorRadius,
					                      shapeConfig.elongation,
					                      shapeConfig.cosCoefficients,
					                      shapeConfig.sinCoefficients ) );

			subdomain = buildSubdomain( background, *shape,
			                            backgroundCellSize( config->getMesh() ) );
		}

		solveMesh = subdomain.mesh ? subdomain.mesh.get() : &background;
		mfem::Mesh &mesh = *solveMesh;

		solver = std::make_unique<meq::GradShafranovSolver>(
			mesh, config->getDiscretisation().polynomialDegree,
			config->getDiscretisation().tau );

		if ( subdomain.mesh )
			solver->setExtension( *subdomain.path, subdomain.gammaHMarker );
		solver->setSource( *source );
		solver->setBoundaryData( zero );
		solver->setNewtonControl( config->getSolver().newtonRelativeTolerance,
		                          config->getSolver().newtonAbsoluteTolerance,
		                          config->getSolver().newtonMaxIterations );

		switch ( config->getInitialGuess().type )
		{
			case meq::InitialGuessType::None:
				break;

			case meq::InitialGuessType::Ramp:
			{
				// psi = 0 in the INTERIOR, not on the boundary. See Config.hpp:
				// every GS-2 section 4.2-4.5 source vanishes at psi = 0, so a
				// homogeneous start is a fixed point of the iteration.
				double const amplitude = config->getInitialGuess().amplitude;
				double const zMin = config->getMesh().zMin;
				double const zMax = config->getMesh().zMax;
				ramp = std::make_unique<mfem::FunctionCoefficient>(
					[ amplitude, zMin, zMax ]( mfem::Vector const &x )
					{
						double const half = 0.5*( zMax + zMin );
						double const span = 0.5*( zMax - zMin );
						return span > 0.0 ? amplitude*( x( 1 ) - half )/span : 0.0;
					} );
				solver->setInitialGuess( *ramp );
				break;
			}

			case meq::InitialGuessType::GridFunction:
			{
				guessMesh = std::make_unique<mfem::Mesh>(
					config->getInitialGuess().meshFile.c_str(), 1, 1 );
				std::ifstream stream( config->getInitialGuess().file );
				if ( !stream )
					throw std::runtime_error( "cannot read [initialguess] File \""
					                          + config->getInitialGuess().file + "\"" );
				guess = std::make_unique<mfem::GridFunction>( guessMesh.get(), stream );

				// The EXACT restart only: same mesh, same degree. The
				// interpolating restart needs FindPointsGSLIB and is
				// DRIVER-PLAN.md section 4's second route, not written yet.
				// Refuse rather than interpolate badly and call it warm.
				if ( guessMesh->GetNE() != mesh.GetNE() )
					throw std::runtime_error(
						"[initialguess] MeshFile has " + std::to_string( guessMesh->GetNE() )
						+ " elements where the run's mesh has " + std::to_string( mesh.GetNE() )
						+ ". Only the exact restart -- same mesh, same degree -- is "
						"implemented; the interpolating one needs GSLIB" );

				solver->setInitialGuess( *guess );
				break;
			}
		}
	}
	catch ( std::exception const &error )
	{
		std::fprintf( stderr, "meq: %s\n", error.what() );
		return ConfigurationError;
	}

	// ---- solve ---------------------------------------------------------
	try
	{
		solver->solve();
	}
	catch ( std::exception const &error )
	{
		reportResiduals( solver->newtonResiduals() );
		std::fprintf( stderr, "meq: the solve did not converge: %s\n", error.what() );
		return SolveFailed;
	}

	if ( shape )
		std::printf( "meq: curved Gamma: %d of %d background elements inside, "
		             "%d transfer paths widened\n",
		             solveMesh->GetNE(), background.GetNE(), subdomain.widened );

	std::printf( "meq: converged in %d Newton iterations on %d elements, "
	             "degree %d\n",
	             solver->newtonIterations(), solveMesh->GetNE(),
	             config->getDiscretisation().polynomialDegree );
	reportResiduals( solver->newtonResiduals() );

	// ---- write ---------------------------------------------------------
	try
	{
		meq::OutputConfig const &output = config->getOutput();
		std::string const stem = output.directory + "/" + output.prefix;

		meq::writeMfem( stem, *solveMesh, solver->potential(), solver->flux() );

		// B_poloidal, a RELABELLING of the solved flux and not a derivative of
		// psi -- B_R = -q_z, B_Z = +q_r. That is the payoff for the mixed
		// method: the field comes out at the same order as the potential
		// rather than one order down. See Field.hpp.
		mfem::GridFunction field( solver->flux().FESpace() );
		meq::poloidalField( solver->flux(), field );

		// The grid is [output] GridNR x GridNZ, NOT the mesh's NR/NZ -- those are
		// pre-refinement cell counts of the SOLVE and have nothing to do with
		// how finely the answer should be sampled for a plot.
		//
		// It spans GAMMA's bounding box where there is one, not the background
		// box: on the curved path most of the background lies outside the plasma
		// and sampling it would spend the grid on NaN.
		double gridRMin = config->getMesh().rMin, gridRMax = config->getMesh().rMax;
		double gridZMin = config->getMesh().zMin, gridZMax = config->getMesh().zMax;
		if ( shape )
			shape->boundingBox( gridRMin, gridRMax, gridZMin, gridZMax );

		meq::GridSampler sampler( *solveMesh,
			gridRMin, gridRMax, output.gridNR,
			gridZMin, gridZMax, output.gridNZ );

		// NaN outside the domain, which is what the file format documents and
		// what a plotting library will mask for free. A zero there would be a
		// physical claim, and a wrong one.
		double const outside = std::numeric_limits<double>::quiet_NaN();
		std::vector<double> psi, bR, bZ;
		sampler.sample( solver->potential(), psi, outside );
		sampler.sampleComponent( field, 0, bR, outside );
		sampler.sampleComponent( field, 1, bZ, outside );

		meq::NetCDFWriter writer( stem + ".nc", sampler );
		writer.attribute( "title", "meq equilibrium" );
		writer.attribute( "meq_version", MEQ_VERSION );
		// The file says which input produced it and how well. A directory of
		// scan output is unreadable otherwise, and "which commit was this?" is
		// the first question asked of any result that looks wrong.
		writer.attribute( "config_file", argument );
		writer.attribute( "polynomial_degree",
		                  config->getDiscretisation().polynomialDegree );
		writer.attribute( "refinement_levels", config->getMesh().refinementLevels );
		writer.attribute( "elements", solveMesh->GetNE() );
		if ( shape )
		{
			// So a reader can tell a curved run from a fitted one without
			// re-deriving Gamma from the configuration.
			writer.attribute( "boundary", "curved (transferred datum on Gamma_h)" );
			writer.attribute( "paths_widened", subdomain.widened );
		}
		else
		{
			writer.attribute( "boundary", "fitted (psi = 0 on the mesh boundary)" );
		}
		writer.attribute( "newton_iterations", solver->newtonIterations() );
		writer.attribute( "final_residual",
		                  solver->newtonResiduals().empty()
		                      ? 0.0 : solver->newtonResiduals().back() );
		writer.field( "psi", psi, "poloidal flux function", "Wb/rad" );
		writer.field( "B_R", bR, "poloidal field, R component", "T" );
		writer.field( "B_Z", bZ, "poloidal field, Z component", "T" );
		writer.close();

		std::printf( "meq: wrote %s.mesh, %s_psi.gf, %s_grad_psi.gf, %s.nc "
		             "(%d of %d grid nodes inside the domain)\n",
		             stem.c_str(), stem.c_str(), stem.c_str(), stem.c_str(),
		             sampler.locatedCount(),
		             sampler.nodesR()*sampler.nodesZ() );
	}
	catch ( std::exception const &error )
	{
		std::fprintf( stderr, "meq: could not write output: %s\n", error.what() );
		return OutputFailed;
	}

	return Solved;
}
