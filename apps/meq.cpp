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
#include "meq/Estimator.hpp"
#include "meq/Field.hpp"
#include "meq/GradShafranov.hpp"
#include "meq/Output.hpp"
#include "meq/RotatingSource.hpp"
#include "meq/Sampler.hpp"
#include "meq/SourceFactory.hpp"

#include <algorithm>
#include <cctype>
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
	 * A WARNING ABOUT MKL_THREADING_LAYER STOOD HERE AND IS DELETED.
	 *
	 * It fired on every run where that variable was unset, telling the user
	 * that UMFPACK's BLAS-3 might be silently wrong. That was true while meq
	 * resolved BLAS through Debian's libblas.so.3, an alternatives symlink to
	 * the libmkl_rt dispatcher. meq now builds against its own SuiteSparse,
	 * which links oneAPI's threading layer directly, so there is no dispatcher
	 * to misconfigure and the variable is inert.
	 *
	 * Keeping it would have been worse than useless: most builds link no MKL
	 * at all, so the warning told the majority of users to set a variable
	 * naming a library they do not have, about a failure they cannot suffer.
	 * A warning that is usually wrong is one people learn to ignore, which
	 * spends the credibility of every other message this driver prints.
	 */

	/// What the file said [source] Type was, for the output's provenance. A
	/// directory of scan output is unreadable without it: every other attribute
	/// describes the discretisation, and two runs differing only in the physics
	/// would otherwise be indistinguishable from their files alone.
	char const *sourceTypeName( meq::SourceType type )
	{
		switch ( type )
		{
			case meq::SourceType::Soloviev:     return "soloviev";
			case meq::SourceType::MHD:          return "mhd";
			case meq::SourceType::Manufactured: return "manufactured";
			case meq::SourceType::Rotating:     return "rotating";
		}
		return "unknown";
	}

	/// A species name made safe to use as a NetCDF variable name. Names come
	/// from the TOML, where nothing stops "C 6+", and netCDF-4 would refuse the
	/// space -- AFTER the solve, which is an expensive place to learn about a
	/// typo. Anything that is not alphanumeric or '_' becomes '_', and a name
	/// that does not start with a letter gains one.
	std::string variableName( std::string const &name )
	{
		std::string safe;
		safe.reserve( name.size() + 1 );
		for ( char const c : name )
			safe += ( std::isalnum( static_cast<unsigned char>( c ) ) || c == '_' )
				? c : '_';
		if ( safe.empty() || !std::isalpha( static_cast<unsigned char>( safe[ 0 ] ) ) )
			safe = "s" + safe;
		return safe;
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

	/// Negative inside Omega, which is the sign convention every piece of the
	/// extension machinery uses. Built once and shared: buildSubdomain() below
	/// and meq::AdaptiveDomain on the adaptive path must be given the SAME
	/// function, or D_h differs between a one-shot run and cycle 0 of an adaptive
	/// one for no reason a user could see.
	mfem::PositionFunction levelSetOf( meq::BoundaryShape const &shape )
	{
		meq::BoundaryShape const *shapePointer = &shape;
		return [ shapePointer ]( mfem::Vector const &x )
		{
			return shapePointer->levelSet( x( 0 ), x( 1 ) );
		};
	}

	Subdomain buildSubdomain( mfem::Mesh &background,
	                          mfem::PositionFunction const &levelSet,
	                          double h )
	{
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

	/// The marking step, GS-2 section 3.2. Doerfler is what the convergence
	/// analysis of Cockburn, Nochetto and Zhang assumes and is the default;
	/// maximum is what GS-2's own experiments used, at gamma = 0.3, and its
	/// analysis for HDG is still open. The two respond to gamma in OPPOSITE
	/// directions -- see meq::markMaximum -- which is why the configuration names
	/// the strategy rather than inferring it from the value.
	void markElements( meq::MarkingStrategy strategy, double theta,
	                   mfem::Vector const &local, mfem::Array<int> &marked )
	{
		if ( strategy == meq::MarkingStrategy::Doerfler )
			meq::markDoerfler( local, theta, marked );
		else
			meq::markMaximum( local, theta, marked );
	}

	/// One turn of the loop, kept for the report at the end.
	struct Cycle
	{
		int elements;
		int traceDofs;
		int marked;
		int widened;
		double eta;
		int iterations;
		bool globalised;
	};

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


	// ---- configuration -------------------------------------------------
	std::unique_ptr<meq::Configuration> config;
	std::shared_ptr<meq::Source const> source;
	/*
	 * A SECOND HANDLE ON THE SAME OBJECT, AND BOTH ARE NEEDED.
	 *
	 * A source whose profiles are functions of NORMALISED flux
	 * Psi = psi/psi_ax makes psi_ax a functional of the solution, so it is an
	 * UNKNOWN of the non-linear system rather than data, and the solver closes
	 * the pair by a bordered Newton. That path wants a non-const
	 * meq::NormalisedSource, because setNormalisation() is called on the source
	 * before every residual evaluation.
	 *
	 * Everything else here -- the estimator above all -- wants a
	 * meq::Source const &, and a normalised source IS one. So `source` aliases
	 * `normalised` when there is one, and nothing downstream has to know which
	 * kind it was given.
	 */
	std::shared_ptr<meq::NormalisedSource> normalised;
	/*
	 * psi_ax, CARRIED FORWARD BETWEEN ADAPTIVE CYCLES.
	 *
	 * makeSolver() builds a FRESH solver every cycle -- it must, since the mesh
	 * changes -- and re-arms the source with this value. Re-arming with the
	 * TOML's guess every time would throw away the one thing the previous cycle
	 * established: psi_ax is a physical quantity, converged on the coarser mesh
	 * to within that mesh's discretisation error, and it is the best guess
	 * available for the finer one. It is the interpolating warm start the FIELD
	 * cannot have -- DRIVER-PLAN section 4 needs GSLIB for that -- arriving free
	 * for the one unknown that happens to be a scalar.
	 *
	 * WHAT IT IS MEASURED TO BUY IS NOTHING YET, AND THAT IS WORTH SAYING RATHER
	 * THAN LEAVING FOR SOMEBODY TO REDISCOVER. On examples/rotating-normalised.toml
	 * run adaptively, every cycle takes the same 6 Newton steps whether the
	 * border starts from the file's 0.3 or from the converged 0.354775928 -- and
	 * a one-shot run handed the exact answer as its PsiAxis takes 6 as well. The
	 * iteration is dominated by the FIELD, which starts cold on every cycle
	 * because carrying it forward is the interpolation that is not written. So
	 * this is the right thing to do rather than a fast one, and it will start
	 * paying only when the field is warm too.
	 */
	double psiAxisGuess = 0.0;
	try
	{
		config = std::make_unique<meq::Configuration>( argument );
		if ( config->getSource().isNormalised() )
		{
			normalised = meq::makeNormalisedSource( config->getSource(), argument );
			psiAxisGuess = config->getSource().psiAxisGuess();
			source = normalised;
		}
		else
		{
			source = meq::makeSource( config->getSource(), argument );
		}

		if ( config->getBoundary().type == meq::BoundaryDataType::Exact )
		{
			std::fprintf( stderr,
				"meq: [boundary] Type = \"exact\" needs the source's closed-form\n"
				"     solution, which meq::Source does not carry -- it is a\n"
				"     convergence-study device and lives in the test fixtures.\n"
				"     Use Type = \"zero\" for the fixed-boundary problem.\n" );
			return ConfigurationError;
		}
	}
	catch ( std::exception const &error )
	{
		std::fprintf( stderr, "meq: %s\n", error.what() );
		return ConfigurationError;
	}

	// ---- set the run up ------------------------------------------------
	mfem::Mesh background;
	std::unique_ptr<meq::BoundaryShape> shape;
	mfem::PositionFunction levelSet;
	mfem::ConstantCoefficient zero( 0.0 );
	std::unique_ptr<mfem::FunctionCoefficient> ramp;
	std::unique_ptr<mfem::Mesh> guessMesh;
	std::unique_ptr<mfem::GridFunction> guess;

	/*
	 * REBUILT EVERY ADAPTIVE CYCLE, AND THE LAST ONE IS WHAT THE WRITE PHASE
	 * READS, so these are declared out here and must outlive the loop.
	 *
	 * The ordering inside the loop is the part that is easy to get wrong.
	 * AdaptiveDomain::refine() builds a NEW SubMesh -- element numbering does not
	 * survive a refinement and nothing there pretends it does -- so the transfer
	 * path, the solver and every finite element space built on the old one are
	 * dangling the moment it is called. They are therefore destroyed BEFORE
	 * refine(), and rebuilt after it.
	 * tests/convergence/AdaptiveRefinement.cpp scopes them for the same reason.
	 */
	std::unique_ptr<meq::GradShafranovSolver> solver;
	/// Curved AND adaptive. See the branch below for why this is not the only
	/// construction of D_h.
	std::unique_ptr<meq::AdaptiveDomain> domain;
	/// Curved, adaptive or not; rebuilt with the domain.
	std::unique_ptr<mfem::VertexConePath> path;
	/// Curved and NOT adaptive: the one-shot construction, kept as it was.
	Subdomain subdomain;
	/// Null on the fitted path, which is how everything downstream tells the two
	/// apart -- including whether the estimator has to exclude Gamma_h.
	mfem::Array<int> const *gammaHMarker = nullptr;
	int widened = 0;
	// The mesh actually solved on: D_h on the curved path, the background mesh
	// on the fitted one. Everything downstream uses it.
	mfem::Mesh *solveMesh = nullptr;

	meq::AdaptivityConfig const &adapt = config->getAdaptivity();
	// MaxIterations counts SOLVES, not refinements: a run with MaxIterations = 1
	// is a plain single solve with an estimate printed, and MaxIterations = 10
	// refines at most nine times.
	int const maxCycles = adapt.enabled ? adapt.maxIterations : 1;

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
			levelSet = levelSetOf( *shape );

			/*
			 * TWO CONSTRUCTIONS OF D_h, AND THE DIFFERENCE IS NOT COSMETIC.
			 *
			 * A one-shot curved run wants nothing but T_h, and buildSubdomain()
			 * gives it, with the two validations worded in terms of the TOML keys
			 * a user would have to change.
			 *
			 * An adaptive one needs the COMPANION mesh of GS-2 section 3.3 as
			 * well, and that is the whole content of the difference. Refine an
			 * element of T_h and its children are still inside Omega, so Gamma_h
			 * does not move: the gap to Gamma stays where it was while h_loc
			 * halves, and dist/h_loc DOUBLES every cycle -- measured at 0.98,
			 * 1.97, 3.94, 7.88 with the domain held fixed against 0.98, 1.02,
			 * 1.39, 1.36 with the companion update in place. The transfer is only
			 * optimal while that ratio is O(1), so without the companion mesh an
			 * adaptive curved run silently leaves the regime the method is
			 * analysed in. meq::AdaptiveDomain is that update.
			 *
			 * So the non-adaptive path is left exactly as it was rather than
			 * routed through AdaptiveDomain for tidiness:
			 * theDriverSolvesOnACurvedBoundary pins it against the library at
			 * 1.6e-16, and the two differ in the transfer path's search length --
			 * six h against twelve times the largest element, which is the right
			 * number on a graded mesh and a needlessly long one on a uniform.
			 */
			if ( adapt.enabled )
			{
				try
				{
					domain = std::make_unique<meq::AdaptiveDomain>( background, levelSet );
				}
				catch ( std::exception const &error )
				{
					throw std::runtime_error(
						std::string( "[boundary.shape] " ) + error.what()
						+ " -- either the [mesh] box is too coarse for the surface, "
						"or the surface is not strictly inside it" );
				}
			}
			else
			{
				subdomain = buildSubdomain( background, levelSet,
				                            backgroundCellSize( config->getMesh() ) );
				solveMesh = subdomain.mesh.get();
				gammaHMarker = &subdomain.gammaHMarker;
				path = std::move( subdomain.path );
				widened = subdomain.widened;
			}
		}

		if ( !shape )
			solveMesh = &background;

		// The guess objects are built once. The RAMP is a coefficient and so is
		// valid on every mesh; the GRID FUNCTION is not, and is applied on the
		// first cycle only -- see makeSolver below.
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
				break;
			}
		}
	}
	catch ( std::exception const &error )
	{
		std::fprintf( stderr, "meq: %s\n", error.what() );
		return ConfigurationError;
	}

	/*
	 * A FRESH SOLVER, FROM SCRATCH, AND BOTH CALLERS NEED IT TO BE.
	 *
	 * Every adaptive cycle solves on a different mesh, and GradShafranovSolver
	 * builds its three spaces in its constructor -- DarcyForm's hybridization
	 * takes what it finds when EnableHybridization() runs -- so a refined mesh
	 * needs a new solver rather than an Update(). That is what CLAUDE.md records
	 * the deleted Solution::Prolong() and Update() as having been for.
	 *
	 * And the reactive ladder needs it for a different reason: with
	 * MFEM_USE_EXCEPTIONS a failed solve throws from the middle of MFEM -- a NaN
	 * detected inside NewtonSolver::Mult, or deeper still from an element-local
	 * solve -- and leaves its objects as the throw found them. CLAUDE.md is
	 * explicit that a GradShafranovSolver must not be assumed reusable
	 * afterwards, so the retry builds another one.
	 */
	auto makeSolver = [ & ]( mfem::Mesh &mesh, bool firstCycle )
		-> std::unique_ptr<meq::GradShafranovSolver>
	{
		auto fresh = std::make_unique<meq::GradShafranovSolver>(
			mesh, config->getDiscretisation().polynomialDegree,
			config->getDiscretisation().tau );

		if ( gammaHMarker )
			fresh->setExtension( *path, *gammaHMarker );

		// The bordered Newton, or the plain one. psiAxisGuess is the TOML's
		// value on the first cycle and the previous cycle's answer afterwards;
		// see its declaration.
		if ( normalised )
			fresh->setSource( *normalised, psiAxisGuess );
		else
			fresh->setSource( *source );

		fresh->setBoundaryData( zero );
		fresh->setNewtonControl( config->getSolver().newtonRelativeTolerance,
		                         config->getSolver().newtonAbsoluteTolerance,
		                         config->getSolver().newtonMaxIterations );

		if ( config->getInitialGuess().type == meq::InitialGuessType::Ramp )
		{
			fresh->setInitialGuess( *ramp );
		}
		else if ( config->getInitialGuess().type == meq::InitialGuessType::GridFunction
		          && firstCycle )
		{
			// The EXACT restart only: same mesh, same degree. The interpolating
			// restart needs FindPointsGSLIB and is DRIVER-PLAN.md section 4's
			// second route, not written yet. Refuse rather than interpolate badly
			// and call it warm.
			//
			// It is a FIRST-CYCLE guess for the same reason: after a refinement
			// the mesh no longer matches, and carrying the previous cycle's answer
			// forward is that same missing interpolation. So cycles after the
			// first start cold -- which costs less than it sounds, since they
			// start on FINER meshes, where Newton is the more reliable and not the
			// less. The coarse first solve is the one at risk, and it is the one
			// the ladder below is for.
			if ( guessMesh->GetNE() != mesh.GetNE() )
				throw std::runtime_error(
					"[initialguess] MeshFile has " + std::to_string( guessMesh->GetNE() )
					+ " elements where the run's mesh has " + std::to_string( mesh.GetNE() )
					+ ". Only the exact restart -- same mesh, same degree -- is "
					"implemented; the interpolating one needs GSLIB" );

			fresh->setInitialGuess( *guess );
		}

		return fresh;
	};

	// ---- solve, and refine if that is what was asked for ----------------
	std::vector<Cycle> history;

	// One VTK frame per adaptive cycle, so the refinement can be watched rather
	// than inferred from the table of element counts at the end. Only for an
	// adaptive run: without the loop there is one state and "<stem>" already
	// holds it. See Output.hpp for why this is a separate collection from the
	// answer, and why its frames are not bent onto Gamma.
	std::unique_ptr<meq::VtuSeries> series;
	if ( adapt.enabled )
		series = std::make_unique<meq::VtuSeries>(
			config->getOutput().directory + "/" + config->getOutput().prefix,
			config->getDiscretisation().polynomialDegree );

	for ( int cycle = 0; cycle < maxCycles; ++cycle )
	{
		// The curved adaptive path rebuilds D_h's dependants every cycle. The
		// other two set these once, in the setup block above.
		if ( domain )
		{
			solveMesh = &domain->computational();
			// TWELVE times the LARGEST element, not the mesh parameter. On a
			// graded mesh the coarse part needs the long search and the fine part
			// is not harmed by being given one; ExtensionConvergence.cpp uses six
			// h on a uniform mesh, where the two are the same number.
			path = std::make_unique<mfem::VertexConePath>(
				domain->computational(), domain->gammaHAttribute(), levelSet,
				12.0*domain->largestElement() );
			widened = path->NumWidened();
			gammaHMarker = &domain->gammaHMarker();
		}

		/*
		 * THE NON-LINEAR PATH THE DRIVER SHIPS: Newton, and on OBSERVED failure
		 * Anderson-Picard into Newton. A REACTIVE ladder, never a predictive one,
		 * and the distinction is load bearing.
		 *
		 * Nothing may be inferred from F about which solver to run, because
		 * nothing can be. The ratio max|dF/dpsi| / lambda_1 is computable from the
		 * black-box interface, dFdPsi being mandatory -- but the GS-2 pressure
		 * pedestal converges at a ratio of 7 where the current hole fails at 26,
		 * which is two points and not a threshold, and the ratio needs the range
		 * of psi, which is not known before solving. A detector calibrated on that
		 * would be fitting noise. Failure is therefore OBSERVED, which needs
		 * nothing from F beyond the existing interface.
		 *
		 * Why this pairing and not a line search: globalising the outer trace
		 * iteration does not globalise the element-local ones, and measured,
		 * KIN_LINESEARCH is WORSE than the undamped iteration on the case that
		 * motivated it -- failing at 18 where plain Newton takes 42, having spent
		 * 1.4M element-local iterations. Anderson-Picard freezes F at the previous
		 * iterate, which leaves every local problem LINEAR, and walks the iterate
		 * into Newton's basin; Newton then supplies the quadratic endgame Picard
		 * structurally cannot. Measured on three cases that plain Newton cannot
		 * reach at coarse resolution, it finishes in 3 to 28 Newton steps.
		 *
		 * It is not cheap -- stage 1 spends 122 to 290 full linear solves -- which
		 * is exactly why it is the fallback and not the default.
		 */
		bool globalised = false;
		try
		{
			// Built here rather than inside the retry so that a configuration
			// which cannot produce a solver at all stays exit 1. That has nothing
			// to do with whether Newton converges.
			solver = makeSolver( *solveMesh, cycle == 0 );
		}
		catch ( std::exception const &error )
		{
			std::fprintf( stderr, "meq: %s\n", error.what() );
			return ConfigurationError;
		}

		try
		{
			solver->solve();
		}
		catch ( std::exception const &firstAttempt )
		{
			reportResiduals( solver->newtonResiduals() );

			/*
			 * THE LADDER HAS NO SECOND RUNG WHEN psi_ax IS AN UNKNOWN, and
			 * saying so beats letting the retry throw.
			 *
			 * Every globalisation meq has drives a residual of its own: the
			 * KINSOL paths solve oper(x) = 0 through their own adapter, and the
			 * Picard ones build no Jacobian at all. The border is a row and a
			 * column ON THE NEWTON JACOBIAN, so neither has anywhere to put it,
			 * and GradShafranovSolver refuses the combination -- correctly, and
			 * with a std::logic_error that would surface here as
			 * "Picard-then-Newton did not converge either", which is a
			 * diagnosis of the wrong thing entirely.
			 */
			if ( normalised )
			{
				std::fprintf( stderr,
					"meq: Newton did not converge: %s\n"
					"     [source] Normalised = true makes psi_ax an unknown, closed by a\n"
					"     BORDERED Newton, and no globalisation meq has can carry that\n"
					"     border -- so there is no fallback to try. The levers are a\n"
					"     better [source] PsiAxis guess, an [initialguess] that puts psi\n"
					"     near the right size, and RESOLUTION.\n",
					firstAttempt.what() );
				return SolveFailed;
			}

			std::fprintf( stderr,
				"meq: Newton did not converge: %s\n"
				"     Retrying with Anderson-Picard to reach Newton's basin, then\n"
				"     Newton for the endgame. This is the observed-failure fallback,\n"
				"     and it costs hundreds of linear solves.\n",
				firstAttempt.what() );

			try
			{
				solver = makeSolver( *solveMesh, cycle == 0 );
				solver->setGlobalisation(
					meq::GradShafranovSolver::Globalisation::PicardThenNewton );
				solver->solve();
				globalised = true;
			}
			catch ( std::exception const &secondAttempt )
			{
				if ( solver )
					reportResiduals( solver->newtonResiduals() );
				std::fprintf( stderr,
					"meq: Picard-then-Newton did not converge either: %s\n"
					"     The remedy for a hard source is RESOLUTION -- try a finer\n"
					"     [mesh], a higher [discretisation] Degree, or [adaptivity].\n",
					secondAttempt.what() );
				return SolveFailed;
			}
		}

		// THIS CYCLE'S ANSWER IS THE NEXT CYCLE'S GUESS. See psiAxisGuess's
		// declaration: a fresh solver per cycle would otherwise re-arm the
		// border with the file's original number and discard what the coarser
		// mesh had already established.
		if ( normalised )
			psiAxisGuess = solver->psiAxis();

		Cycle record{ solveMesh->GetNE(), solver->numTraceDofs(), 0, widened,
		              -1.0, solver->newtonIterations(), globalised };

		if ( !adapt.enabled )
		{
			history.push_back( record );
			break;
		}

		/*
		 * THE ESTIMATE, ON psi* -- WHICH IT COULD NOT USE UNTIL MFEM WAS FIXED.
		 *
		 * Eq (20) builds four of its five terms on the post-processed potential.
		 * ReconstructFluxAndPot() used to skip its mean-value regularisation
		 * whenever a non-linear potential integrator was merely PRESENT rather
		 * than when it had contributed, so on any element where dF/dpsi vanished
		 * the local problem was singular and was factored anyway -- and psi* there
		 * was a different function, 20x to 64x out. meq's Newton path puts every
		 * source on the non-linear form, so that reached every configuration the
		 * driver could be given, and this loop ran on Potential::Raw instead: one
		 * order down at every k, correct but blunt.
		 *
		 * The fix landed as "The postprocessing closes on the element average,
		 * always" -- the close is unconditional now, because the local problem is
		 * a pure Neumann one by construction and there was never anything to
		 * decide. Measured from meq's side, the same case that read 20.3, 64.1 and
		 * 61.6 now reads 1.0069 against 1.0048 for elements where dF/dpsi does not
		 * vanish at all: psi* is a post-processing everywhere.
		 * NewtonConvergence.cpp's
		 * thePostProcessedPotentialIsCorrectWhereTheJacobianVanishes is what says
		 * so, and it is what would say if it came back.
		 */
		solver->postProcess();
		meq::ResidualEstimator estimator( *solver, *source );

		/*
		 * ON THE EXTENSION PATH eta_5 HAS TO LEAVE Gamma_h OUT, AND THE DRIVER
		 * SAYS SO RATHER THAN LEAVING THE USER TO FIND OUT.
		 *
		 * On such a face psihat_h is not the condition that was imposed: phi_h
		 * is, and the trace dofs there are pinned to zero because nothing
		 * references them. So eta_5 compares the potential against zero and the
		 * difference is O( dist( Gamma_h, Gamma ) ) = O( h ), not O( h^(k+2) ).
		 * Unmitigated that gives eta = 4.09e-1 where eta_1 = 2.12e-3, converging at
		 * about 0.5 -- the loop would run, produce plausible pictures and refine
		 * the WRONG ELEMENTS, which is the quietest possible way for this to fail.
		 *
		 * THIS IS AN OMISSION, NOT A REPAIR. The proper fix is to evaluate phi_h,
		 * which MFEM now makes reachable through TransferredDatumCoefficient; eta_5
		 * has not yet been rebuilt on it, and that wants its own convergence
		 * measurement rather than a switch.
		 */
		if ( gammaHMarker )
			estimator.setTransferredBoundary( *gammaHMarker );

		mfem::Vector const &local = estimator.GetLocalErrors();
		record.eta = estimator.GetTotalError();

		mfem::Array<int> marked;
		bool const lastCycle = cycle + 1 == maxCycles;
		bool const reachedTarget = record.eta <= adapt.targetError;

		if ( !lastCycle && !reachedTarget )
		{
			markElements( adapt.strategy, adapt.theta, local, marked );
			record.marked = marked.Size();
		}
		history.push_back( record );

		// This cycle's frame, written while the solver still owns its fields --
		// a few lines below they are destroyed so the mesh can be refined.
		if ( series )
		{
			mfem::GridFunction cycleField( solver->flux().FESpace() );
			meq::poloidalField( solver->flux(), cycleField );
			series->append( *solveMesh, solver->potential(), cycleField,
			                cycle, static_cast<double>( cycle ) );
		}

		if ( reachedTarget )
		{
			std::printf( "meq: eta = %.4e is at or below TargetError = %.4e "
			             "after %d cycle%s\n", record.eta, adapt.targetError,
			             cycle + 1, cycle == 0 ? "" : "s" );
			break;
		}
		if ( lastCycle )
		{
			std::printf( "meq: reached MaxIterations = %d with eta = %.4e, above "
			             "TargetError = %.4e. The answer is the finest one "
			             "computed, not a converged one.\n",
			             maxCycles, record.eta, adapt.targetError );
			break;
		}
		if ( marked.Size() == 0 )
		{
			std::printf( "meq: the marking strategy selected no elements at "
			             "eta = %.4e, so there is nothing to refine\n", record.eta );
			break;
		}

		// EVERYTHING BUILT ON THIS MESH DIES FIRST. See the declarations above.
		solver.reset();
		path.reset();
		gammaHMarker = nullptr;

		try
		{
			if ( domain )
			{
				domain->refine( marked );
			}
			else
			{
				// Conforming refinement, so it propagates beyond the marked set.
				solveMesh->GeneralRefinement( marked );
			}
		}
		catch ( std::exception const &error )
		{
			/*
			 * A refinement that cannot be carried out ends the run rather than the
			 * loop, and that is a deliberate choice rather than an easy one. The
			 * cycle just completed had converged and been estimated, so there was
			 * a correct answer on a coarser mesh than was asked for -- but the
			 * solver holding it was released above, before the mesh changed,
			 * because it had to be. Nothing is left to write, so exit 2 is the
			 * honest report.
			 */
			std::fprintf( stderr,
				"meq: cycle %d could not be refined: %s\n"
				"     The previous cycle had converged, but its solver was released\n"
				"     before the mesh changed -- as it must be -- so there is\n"
				"     nothing left to write. Lower [adaptivity] MaxIterations.\n",
				cycle, error.what() );
			return SolveFailed;
		}
	}

	{
		// The background element count is the DOMAIN's after any refinement, not
		// this scope's copy: AdaptiveDomain owns and refines its own.
		int const backgroundElements = domain ? domain->numBackground()
		                                      : background.GetNE();
		if ( shape )
			std::printf( "meq: curved Gamma: %d of %d background elements inside, "
			             "%d transfer paths widened\n",
			             solveMesh->GetNE(), backgroundElements, widened );

		Cycle const &last = history.back();
		std::printf( "meq: converged in %d Newton iterations on %d elements, "
		             "degree %d%s\n",
		             last.iterations, last.elements,
		             config->getDiscretisation().polynomialDegree,
		             last.globalised ? " (via Picard-then-Newton)" : "" );

		/*
		 * BOTH NUMBERS, BECAUSE THEY HAVE DIFFERENT UNITS AND EITHER CAN BE
		 * SATISFIED WITHOUT THE OTHER.
		 *
		 * The bordered system is R( lambda, s ) = 0 together with
		 * G( lambda, s ) = s - max psi_h = 0. R is a trace residual and G is a
		 * flux, so the printed ||r|| above is a weighted combination and not a
		 * statement about either half on its own: a run can drive R to round-off
		 * while psi_ax still disagrees with the peak it is supposed to BE, which
		 * is a solved equation for a plasma nobody asked for. So the constraint
		 * is reported separately, in its own units.
		 */
		if ( normalised )
			std::printf( "     psi_ax = %.6e Wb/rad, constraint psi_ax - max psi_h "
			             "= %.3e\n",
			             solver->psiAxis(), solver->normalisationResidual() );

		if ( adapt.enabled )
		{
			std::printf( "\n  the adaptive loop: %s marking at %.2f\n",
			             adapt.strategy == meq::MarkingStrategy::Doerfler
			                 ? "Doerfler" : "maximum", adapt.theta );
			std::printf( "  %6s %8s %8s %8s %6s %12s %5s\n",
			             "cycle", "elem", "trace", "marked", "wide", "eta", "it" );
			for ( std::size_t c = 0; c < history.size(); ++c )
			{
				Cycle const &e = history[ c ];
				std::printf( "  %6zu %8d %8d %8d %6d %12.4e %5d%s\n",
				             c, e.elements, e.traceDofs, e.marked, e.widened,
				             e.eta, e.iterations, e.globalised ? "  P->N" : "" );
			}
			std::fflush( stdout );
		}

		reportResiduals( solver->newtonResiduals() );
	}

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

		// THE SLIVER BETWEEN Gamma_h AND Gamma, which only the curved path has.
		// Omega_h is the union of background elements lying INSIDE Gamma, so
		// Gamma_h is inscribed and there is a band O(h) wide that is inside the
		// plasma and outside the mesh. Left alone it rasterises as NaN and the
		// picture gets a ragged polygonal edge where the boundary is smooth.
		//
		// Filled by extrapolating the element across each boundary face, and
		// ONLY where the node is inside Gamma -- the level set is what stops
		// this painting a band outside the plasma, where the solve claims
		// nothing. Reach 1.0: one face length, which is the O(h) regime the
		// transfer technique is analysed in and no further.
		int extended = 0;
		if ( shape )
		{
			meq::BoundaryShape const *curve = shape.get();
			extended = sampler.extendOutward( 1.0,
				[ curve ]( double r, double z )
				{
					return curve->levelSet( r, z ) <= 0.0;
				},
				[ curve ]( double r, double z )
				{
					// levelSet is the RADIAL gap and is negative inside, so its
					// magnitude is the distance to Gamma along the ray. Not the
					// perpendicular distance, which is smaller -- but the ratio
					// below only needs the two gaps measured the same way as
					// each other, and both are along the ray.
					return -curve->levelSet( r, z );
				} );
		}

		// NaN outside the domain, which is what the file format documents and
		// what a plotting library will mask for free. A zero there would be a
		// physical claim, and a wrong one.
		double const outside = std::numeric_limits<double>::quiet_NaN();
		std::vector<double> psi, bR, bZ;

		// THE FLUX CARRIES THE BAND, which is the payoff for the mixed method
		// showing up somewhere unexpected. q is computed at the SAME order as
		// psi, so continuing psi across the Gamma_h-to-Gamma band as
		// psi( x0 ) + r q( x0 ) . ( p - x0 ) evaluates both fields at a point on
		// the element's own boundary and never outside it. The alternative --
		// evaluating psi_h's polynomial outside its element -- is bounded by
		// nothing, and was measured crossing psi = 0, a value this
		// configuration knows exactly.
		sampler.samplePotentialWithFlux( solver->potential(), solver->flux(),
		                                 psi, outside );
		sampler.sampleComponent( field, 0, bR, outside );
		sampler.sampleComponent( field, 1, bZ, outside );

		/*
		 * THE ROTATING FIELDS: n_s( R, Z ) AND e phi_0( R, Z ).
		 *
		 * With rotation the density is NOT a flux function -- centrifugal force
		 * sweeps the heavy species outboard and an electrostatic potential
		 * arises to stop that separating the charges -- so n_s and phi_0 are
		 * genuine two-dimensional fields and a rotating equilibrium is not
		 * interpretable without them. Both are ALGEBRAIC in ( r, psi ), so this
		 * costs one pass over the grid and no solve. See FLOW-PLAN.md section
		 * 5.6 and meq::RotatingSource.
		 *
		 * SAMPLED BY HAND, AND NOT WITH sampleCoefficient(), WHICH IS THE TRAP
		 * THIS BLOCK EXISTS TO AVOID. sample(), sampleComponent() and
		 * sampleCoefficient() all evaluate a field at the node's location IN AN
		 * ELEMENT, and for a node in the Gamma_h-to-Gamma band that location is
		 * the node's FOOT on Gamma_h, not the node itself: only
		 * samplePotentialWithFlux() applies the Taylor step outward. For psi the
		 * difference is O( h ) and visible; for n_s it is worse, because the
		 * whole physics of RoPP (96) is the r^2 in the exponent, so being handed
		 * the foot's radius rather than the node's misplaces the centrifugal
		 * enrichment by a band width. So the loop below uses the ALREADY
		 * SAMPLED psi -- which carries the flux continuation -- together with
		 * sampler.rAt( i ), the node's own R.
		 *
		 * The guard keeps the fields consistent with `inside`: a node the
		 * sampler did not locate, or one whose psi came back NaN, keeps the fill
		 * value here too. OutputConvergence's theMaskAgreesWithTheData is what
		 * says that matters.
		 */
		std::vector<std::vector<double>> densities;
		std::vector<double> ePhi;
		if ( config->getSource().type == meq::SourceType::Rotating )
		{
			meq::RotatingParameters const &rotating = config->getSource().getRotating();

			// One of the two is non-null, and both take PHYSICAL psi and convert
			// internally, so the loop below does not care which path it is on.
			meq::RotatingSource const *plain =
				dynamic_cast<meq::RotatingSource const *>( source.get() );
			meq::NormalisedRotatingSource const *scaled =
				dynamic_cast<meq::NormalisedRotatingSource const *>( source.get() );

			if ( !plain && !scaled )
				throw std::runtime_error(
					"[source] Type = \"rotating\" did not produce a rotating source, so "
					"the density and potential fields cannot be written" );

			std::size_t const count =
				static_cast<std::size_t>( sampler.nodesR() )*sampler.nodesZ();
			densities.assign( rotating.species.size(),
			                  std::vector<double>( count, outside ) );
			ePhi.assign( count, outside );

			for ( int j = 0; j < sampler.nodesZ(); ++j )
				for ( int i = 0; i < sampler.nodesR(); ++i )
				{
					std::size_t const at =
						static_cast<std::size_t>( j )*sampler.nodesR() + i;
					if ( !sampler.located( i, j ) || !std::isfinite( psi[ at ] ) )
						continue;

					double const r = sampler.rAt( i );
					ePhi[ at ] = plain ? plain->potential( r, psi[ at ] )
					                   : scaled->potential( r, psi[ at ] );
					for ( std::size_t sp = 0; sp < rotating.species.size(); ++sp )
						densities[ sp ][ at ] =
							plain ? plain->density( sp, r, psi[ at ] )
							      : scaled->density( sp, r, psi[ at ] );
				}
		}

		meq::NetCDFWriter writer( stem + ".nc", sampler );
		writer.attribute( "title", "meq equilibrium" );
		writer.attribute( "meq_version", MEQ_VERSION );
		// The file says which input produced it and how well. A directory of
		// scan output is unreadable otherwise, and "which commit was this?" is
		// the first question asked of any result that looks wrong.
		writer.attribute( "config_file", argument );
		// WHICH PHYSICS, not just which discretisation. DRIVER-PLAN section 3
		// lists this and it was not being written: two runs differing only in
		// [source] Type produced files that were identical in every attribute.
		writer.attribute( "source_type", sourceTypeName( config->getSource().type ) );
		writer.attribute( "polynomial_degree",
		                  config->getDiscretisation().polynomialDegree );
		writer.attribute( "refinement_levels", config->getMesh().refinementLevels );
		writer.attribute( "elements", solveMesh->GetNE() );
		if ( shape )
		{
			// So a reader can tell a curved run from a fitted one without
			// re-deriving Gamma from the configuration.
			writer.attribute( "boundary", "curved (transferred datum on Gamma_h)" );
			writer.attribute( "paths_widened", widened );
		}
		else
		{
			writer.attribute( "boundary", "fitted (psi = 0 on the mesh boundary)" );
		}
		writer.attribute( "newton_iterations", solver->newtonIterations() );
		// psi_ax IS AN ANSWER HERE, NOT AN INPUT, so the file has to carry it:
		// a reader given only the [source] PsiAxis guess from the TOML would be
		// reading the starting point of a Newton iteration and calling it the
		// axis flux. normalisation_residual is the half of the augmented
		// residual that says whether it is self consistent -- final_residual
		// below mixes it with the trace residual, which is a different quantity
		// in different units.
		if ( normalised )
		{
			writer.attribute( "psi_axis", solver->psiAxis() );
			writer.attribute( "normalisation_residual",
			                  solver->normalisationResidual() );
		}
		// A reader is entitled to know which nodes are the solution and which
		// are a continuation of it past Gamma_h. Zero on the fitted path.
		writer.attribute( "extrapolated_nodes", extended );
		writer.attribute( "final_residual",
		                  solver->newtonResiduals().empty()
		                      ? 0.0 : solver->newtonResiduals().back() );

		// AN ADAPTIVE RUN'S FILE MUST SAY THAT IT WAS ONE. "elements" above is
		// then the count the loop arrived at rather than anything derivable from
		// [mesh], so a directory of scan output is otherwise unreadable -- and eta
		// is the number that says whether it was worth stopping there.
		//
		// estimator_potential is recorded beside it because the two etas are NOT
		// comparable: the degraded one is an order lower and 124x to 407x larger,
		// so a scan that mixes them and plots eta against dofs is plotting two
		// different quantities.
		if ( adapt.enabled )
		{
			Cycle const &last = history.back();
			writer.attribute( "adaptive_cycles", static_cast<int>( history.size() ) );
			writer.attribute( "adaptive_eta", last.eta );
			writer.attribute( "adaptive_target_error", adapt.targetError );
			writer.attribute( "marking_strategy",
			                  adapt.strategy == meq::MarkingStrategy::Doerfler
			                      ? "doerfler" : "maximum" );
			writer.attribute( "marking_theta", adapt.theta );
			writer.attribute( "estimator_potential", "post-processed" );
		}
		// THE BOUNDARY meq WAS GIVEN, so a plot can draw the answer against
		// what was asked for rather than against the mesh's own edge. On the
		// curved path that is the smooth Gamma, sampled from the shape itself;
		// on the fitted path Gamma IS the mesh boundary, so it is walked out of
		// the mesh. Either way it is what psi = 0 was imposed on.
		{
			std::vector<double> boundaryR, boundaryZ;
			if ( shape )
			{
				// 512 points: enough that a Miller shape with triangularity
				// draws smoothly at any figure size, and trivial beside the
				// grid itself.
				int const samples = 512;
				boundaryR.resize( samples );
				boundaryZ.resize( samples );
				for ( int s = 0; s < samples; ++s )
					shape->point( 2.0*M_PI*s/samples,
					              boundaryR[ s ], boundaryZ[ s ] );
				writer.attribute( "boundary_source", "shape (smooth Gamma)" );
			}
			else
			{
				int unreached = 0;
				meq::boundaryPolyline( *solveMesh, boundaryR, boundaryZ, unreached );
				writer.attribute( "boundary_source", "mesh boundary" );
				if ( unreached > 0 )
					std::fprintf( stderr,
						"meq: warning: the mesh boundary is not a single loop; "
						"%d boundary vertices are not in boundary_R/Z\n",
						unreached );
			}
			writer.boundary( boundaryR, boundaryZ );
		}

		writer.field( "psi", psi, "poloidal flux function", "Wb/rad" );
		writer.field( "B_R", bR, "poloidal field, R component", "T" );
		writer.field( "B_Z", bZ, "poloidal field, Z component", "T" );

		if ( !ePhi.empty() )
		{
			meq::RotatingParameters const &rotating = config->getSource().getRotating();
			for ( std::size_t sp = 0; sp < densities.size(); ++sp )
				writer.field( "n_" + variableName( rotating.species[ sp ].name ),
				              densities[ sp ],
				              "number density of " + rotating.species[ sp ].name
				                  + ", RoPP (96); equals its Density profile on "
				                    "r = ReferenceRadius",
				              "m^-3" );

			// e phi_0 AND NOT phi_0, in JOULES AND NOT VOLTS. It is the
			// elementary charge times the potential, which is the combination
			// every exponent of (96) actually contains, and carrying it that way
			// is what keeps the elementary charge in exactly one place. Divide by
			// 1.602176634e-19 for volts.
			writer.field( "e_phi_0", ePhi,
			              "elementary charge times the electrostatic potential "
			              "phi_0 of RoPP (97), in JOULES not volts; zero on "
			              "r = ReferenceRadius by the gauge",
			              "J" );
		}

		writer.close();

		// VTK LAST, AND FOR A REASON. Curving the mesh changes the map from
		// reference to physical space, so everything that reads geometry --
		// writeMfem(), the sampler -- has to have finished. Nothing below
		// touches the mesh again.
		//
		// On the curved path the boundary is bent out onto the true Gamma, so
		// the drawn domain is Omega rather than the inscribed polygon
		// Omega_h. writeVtu() already emits VTK Lagrange cells, so a
		// curvilinear mesh needs nothing further from the format.
		double curved = 0.0;
		int curvedNodes = 0;
		if ( shape )
		{
			meq::BoundaryShape const *curve = shape.get();
			curvedNodes = meq::curveBoundaryOnto( *solveMesh,
				config->getDiscretisation().polynomialDegree,
				[ curve ]( double r, double z, double &outR, double &outZ )
				{
					// Radial projection onto Gamma. levelSet() IS the radial
					// gap -- |p - centre| minus the curve's radius at the same
					// polar angle -- so subtracting it along the ray lands
					// exactly on the curve. Exact for a star shaped boundary,
					// which BoundaryShape's constructor already requires.
					double const centreR = curve->majorRadius();
					double const centreZ = curve->centreHeight();
					double const vR = r - centreR, vZ = z - centreZ;
					double const rho = std::hypot( vR, vZ );
					outR = r;
					outZ = z;
					if ( rho <= 0.0 )
						return;
					double const gap = curve->levelSet( r, z );
					outR = r - vR*gap/rho;
					outZ = z - vZ*gap/rho;
				}, curved );

			// Not a warning below 100%: a face that cannot reach Gamma without
			// folding its element is a statement about the mesh being coarse
			// there, and the rest of the boundary still reaches. It is reported
			// because the VTK boundary is then not uniformly Gamma, which a
			// reader comparing it against the .nc boundary_R/Z would notice.
			if ( curvedNodes > 0 && curved < 1.0 )
				std::printf( "meq: %.0f%% of the boundary reached Gamma; the "
				             "rest was held back to keep its elements from "
				             "folding\n", 100.0*curved );
			else if ( curvedNodes == 0 )
				std::fprintf( stderr,
					"meq: warning: could not bend the boundary onto Gamma "
					"without tangling an element. The VTK boundary is the "
					"polygon Gamma_h.\n" );
		}

		meq::writeVtu( stem, *solveMesh, solver->potential(), field,
		               config->getDiscretisation().polynomialDegree );

		// The description goes FIRST on each line and the path last, because
		// the path can be long and absolute and any attempt to align a
		// trailing comment against it collapses. Named individually rather
		// than as "<stem>.*" because the three formats are for three different
		// tools; tools/README.md says which is which.
		//
		// The ParaView line names "<stem>/<name>.pvd" and not "<stem>.pvd":
		// that collection writes a directory with the index inside it. See
		// Output.hpp.
		std::string const name =
			stem.substr( stem.find_last_of( '/' ) + 1 );
		std::printf(
			"meq: wrote\n"
			"  exact, for GLVis and restart:  %s.mesh\n"
			"                                 %s_psi.gf\n"
			"                                 %s_grad_psi.gf\n"
			"  VTK at degree %d, ParaView:     %s/%s.pvd\n"
			"  (R, Z) grid, %d/%d inside:  %s.nc\n",
			stem.c_str(), stem.c_str(), stem.c_str(),
			config->getDiscretisation().polynomialDegree,
			stem.c_str(), name.c_str(),
			sampler.locatedCount(), sampler.nodesR()*sampler.nodesZ(),
			stem.c_str() );
		if ( series && series->frames() > 0 )
			std::printf( "  %d refinement frames, scrubbable:  %s_cycles/%s_cycles.pvd\n",
			             series->frames(), stem.c_str(), name.c_str() );
	}
	catch ( std::exception const &error )
	{
		std::fprintf( stderr, "meq: could not write output: %s\n", error.what() );
		return OutputFailed;
	}

	return Solved;
}
