/*
 * MEQ -- the Maryland Equilibrium Solver.
 *
 *     MEQ config.toml
 *
 * One binary, no subcommands, and deliberately NOT mfem::OptionsParser, which
 * the program this replaces used and which wants to own argument parsing for
 * the whole executable. docs/running.rst is what the driver owes the user.
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
#include "meq/Coils.hpp"
#include "meq/Config.hpp"
#include "meq/CriticalPoints.hpp"
#include "meq/Estimator.hpp"
#include "meq/ExteriorDtN.hpp"
#include "meq/Field.hpp"
#include "meq/FluxExtraction.hpp"
#include "meq/FluxFamily.hpp"
#include "meq/GradShafranov.hpp"
#include "meq/Output.hpp"
#include "meq/RotatingSource.hpp"
#include "meq/Sampler.hpp"
#include "meq/SourceFactory.hpp"
#include "meq/WarmStart.hpp"

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
			"MEQ -- the Maryland Equilibrium Solver\n"
			"\n"
			"  MEQ <config.toml>    solve the equilibrium the file describes\n"
			"  MEQ --help           this\n"
			"  MEQ --version        the build this is\n"
			"\n"
			"Exit codes: 0 solved, 1 configuration, 2 solve did not converge,\n"
			"3 output could not be written.\n" );
	}

	/*
	 * A WARNING ABOUT MKL_THREADING_LAYER STOOD HERE AND IS DELETED.
	 *
	 * It fired on every run where that variable was unset, telling the user
	 * that UMFPACK's BLAS-3 might be silently wrong. That was true while MEQ
	 * resolved BLAS through Debian's libblas.so.3, an alternatives symlink to
	 * the libmkl_rt dispatcher. MEQ now builds against its own SuiteSparse,
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
	/// The background mesh parameter, which sets the transfer path's search
	/// length in buildSubdomain().
	///
	/// A MESH READ FROM A FILE HAS TO BE MEASURED RATHER THAN COMPUTED, AND
	/// THIS USED TO TAKE ONLY THE CONFIG. `RMin` .. `ZMax`, `NR`, `NZ` and
	/// `RefinementLevels` describe a box the driver BUILDS; with `[mesh] File`
	/// they are absent, so the formula below reads `( 0 - 0 )/( 0*1 )` and the
	/// search length came out as zero -- and a zero search length is not a
	/// tolerance that merely fails on the odd vertex, it is
	/// mfem::VertexConePath aborting on the FIRST vertex of Gamma_h, since no
	/// ray reaches Gamma within nothing. `[mesh] File` beside a curved boundary
	/// is exactly what FB-6 wants (a gmsh half-disc with the conductors meshed
	/// to), and it was the one path this quantity had never been asked for on.
	///
	/// The measured branch is the LARGEST element diameter, which is what the
	/// adaptive path already uses and is the right number on a graded mesh --
	/// the coarse part needs the long search and the fine part is not harmed by
	/// being given one. The built branch is left computing exactly what it
	/// always computed, bit for bit, so no existing configuration moves: a
	/// search length is a bound on where a ray may look, and widening it can
	/// change which rays succeed and therefore NumWidened().
	double backgroundCellSize( meq::MeshConfig const &config, mfem::Mesh &mesh )
	{
		if ( config.fromFile() )
		{
			double largest = 0.0;
			for ( int e = 0; e < mesh.GetNE(); ++e )
				largest = std::max( largest, meq::elementDiameter( mesh, e ) );
			return largest;
		}

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

	/// The level set of the SEMICIRCLE `[boundary.exterior]` describes: negative
	/// inside, zero on `Gamma`, positive outside.
	///
	/// **THIS IS NOT A meq::BoundaryShape AND CANNOT BE ONE.** That class refuses
	/// a surface reaching `r <= 0`, rightly, since a closed plasma surface
	/// through the axis carries a non-integrable `1/r`. `Gamma` here is an
	/// ARTIFICIAL boundary whose flat side IS the axis, and the axis half of it
	/// is ordinary fitted boundary that SubMesh leaves with its inherited
	/// attribute -- so `Gamma_h` is the arc alone and nothing is ever
	/// transferred across `r = 0`. meq::AdaptiveDomain was relaxed for exactly
	/// this geometry.
	mfem::PositionFunction semicircleLevelSet( double radius, double centreZ )
	{
		return [ radius, centreZ ]( mfem::Vector const &x )
		{
			return std::hypot( x( 0 ), x( 1 ) - centreZ ) - radius;
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

		// The parent's largest boundary attribute, taken BEFORE the cut, because
		// it is what tells a generated attribute from an inherited one below.
		int const parentBoundaryMax = background.bdr_attributes.Size() > 0
			? background.bdr_attributes.Max() : 0;

		mfem::Array<int> domainAttribute( 1 );
		domainAttribute[ 0 ] = 1;

		Subdomain subdomain;
		subdomain.mesh = std::make_unique<mfem::SubMesh>(
			mfem::SubMesh::CreateFromDomain( background, domainAttribute ) );

		// SubMesh gives the boundary it had to generate ONE new attribute, one
		// past whatever the parent already used, and leaves inherited boundary
		// with the attribute it already had. So Gamma_h is the generated one and
		// it is the largest.
		//
		// INHERITED BOUNDARY IS LEGITIMATE AND THIS USED TO REFUSE IT, WHICH IS
		// THE SAME RELAXATION meq::AdaptiveDomain NEEDED AND FOR THE SAME
		// GEOMETRY. The guard required EXACTLY ONE attribute -- Omega strictly
		// inside the box -- which every [boundary.shape] run satisfies and
		// [boundary.exterior] never can: the half-disc's flat side IS the box's
		// r = 0 edge, deliberately, because the exterior expansion is valid only
		// on a semicircle centred on the axis. That edge is the AXIS. It is
		// ordinary fitted boundary, it wants no transfer, and every Gegenbauer
		// mode vanishes on it identically.
		//
		// What has to hold is that SOME boundary was generated -- that there is a
		// Gamma_h to transfer to at all -- and that is what is checked now. A
		// domain aligned with the box on every side generates nothing, leaving
		// the largest attribute inherited and gammaH naming a fitted edge, which
		// is the failure this guard exists to prevent.
		subdomain.gammaH = subdomain.mesh->bdr_attributes.Max();
		if ( subdomain.gammaH <= parentBoundaryMax )
			throw std::runtime_error(
				"no boundary was generated when cutting the subdomain, so Omega is "
				"aligned with the [mesh] box everywhere and there is no Gamma_h to "
				"transfer to" );

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
		std::printf( "MEQ %s (MFEM %s)\n", MEQ_VERSION, MFEM_VERSION_STRING );
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
	 * THE PLASMA SOURCE BEFORE THE COILS WERE ADDED TO IT, KEPT SEPARATELY.
	 *
	 * With `[[coils]]` blocks present, `source` above is a
	 * meq::CoilAugmentedSource wrapping this one, because the coil current is
	 * part of F and everything that reads F -- the solve and the residual
	 * ESTIMATOR above all -- must see the sum. But the rotating output fields
	 * are recovered by dynamic_cast, and a cast to meq::RotatingSource through
	 * a wrapper fails: without this handle a rotating run with coils would
	 * report "did not produce a rotating source" and write no densities.
	 *
	 * It is null when there are no coils, in which case `source` IS the plasma
	 * source and the cast sites below fall back to it.
	 */
	std::shared_ptr<meq::Source const> plasmaSource;
	/// The coils of `[[coils]]`, or null if the file described none.
	std::shared_ptr<meq::CoilSet const> coils;
	/*
	 * psi_ax, CARRIED FORWARD BETWEEN ADAPTIVE CYCLES.
	 *
	 * makeSolver() builds a FRESH solver every cycle -- it must, since the mesh
	 * changes -- and re-arms the source with this value. Re-arming with the
	 * TOML's guess every time would throw away the one thing the previous cycle
	 * established: psi_ax is a physical quantity, converged on the coarser mesh
	 * to within that mesh's discretisation error, and it is the best guess
	 * available for the finer one. It used to be described here as the warm start
	 * the FIELD could not have; the field has it now, by meq::FieldTransfer a few
	 * hundred lines below, so this is the same idea arriving free for the one
	 * unknown that happens to be a scalar rather than a consolation for its
	 * absence.
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
		coils = meq::makeCoilSet( config->getCoils(),
		                          config->getSource().permeability(), argument );

		if ( config->getSource().isNormalised() )
		{
			auto plasma = meq::makeNormalisedSource( config->getSource(), argument );
			psiAxisGuess = config->getSource().psiAxisGuess();

			/*
			 * THE MOVING SUPPORT IS SET ON THE PLASMA SOURCE AND NOT ON THE
			 * WRAPPER, WHICH IS WHY THE WRAPPER OVERRIDES IT. The coils are
			 * OUTSIDE the plasma by construction -- confining them to it would
			 * switch off every coil in the machine -- so the confinement
			 * applies to the plasma term alone and the sum is taken after it.
			 * meq::CoilAugmentedNormalisedSource::setPlasmaSupport forwards,
			 * and meq::NormalisedSource::setPlasmaSupport is virtual so that
			 * this call reaches the forwarding one whichever handle it goes
			 * through.
			 */
			if ( config->getSource().confinesToPlasma() )
				plasma->setPlasmaSupport( true );

			if ( coils )
			{
				plasmaSource = plasma;
				normalised = std::make_shared<meq::CoilAugmentedNormalisedSource>(
					std::move( plasma ), coils );
			}
			else
			{
				normalised = std::move( plasma );
			}

			source = normalised;
		}
		else
		{
			auto plasma = meq::makeSource( config->getSource(), argument );

			if ( coils )
			{
				plasmaSource = plasma;
				source = std::make_shared<meq::CoilAugmentedSource const>(
					std::move( plasma ), coils );
			}
			else
			{
				source = std::move( plasma );
			}
		}

		if ( config->getBoundary().type == meq::BoundaryDataType::Exact )
		{
			std::fprintf( stderr,
				"MEQ: [boundary] Type = \"exact\" needs the source's closed-form\n"
				"     solution, which meq::Source does not carry -- it is a\n"
				"     convergence-study device and lives in the test fixtures.\n"
				"     Use Type = \"zero\" for the fixed-boundary problem.\n" );
			return ConfigurationError;
		}
	}
	catch ( std::exception const &error )
	{
		std::fprintf( stderr, "MEQ: %s\n", error.what() );
		return ConfigurationError;
	}

	/*
	 * ---- the two performance keys, resolved before any work is done ----
	 *
	 * [solver] AssemblyMode and TraceSolver are the only keys in the file whose
	 * validity depends on how MFEM was BUILT rather than on what the file says.
	 * Config parses the spelling and stops there, deliberately: it is one of the
	 * four translation units CI compiles without MFEM at all, so it cannot ask
	 * whether this build has OpenMP, or PARDISO, or cuDSS.
	 *
	 * So the question is asked here, and asked ONCE, before a mesh is built or a
	 * source is assembled -- a run that is going to be refused for its build
	 * should be refused in milliseconds and not after the first cycle. That also
	 * means setAssemblyMode() and setTraceSolver() inside makeSolver() cannot
	 * throw, which matters because makeSolver() runs once per adaptive cycle and
	 * an exception there would be indistinguishable from a failed solve.
	 *
	 * Exit code 1, because an unavailable solver is a configuration fault: the
	 * file asks for something this binary cannot do. The message names the CMake
	 * option rather than the symptom, since the fix is a rebuild.
	 */
	using AM = meq::GradShafranovSolver::AssemblyMode;
	using TS = meq::GradShafranovSolver::TraceSolver;

	AM assemblyMode =
		( config->getSolver().assemblyMode == meq::AssemblyModeType::Threaded )
		? AM::Threaded : AM::Serial;

	TS traceSolver = TS::UMFPack;
	switch ( config->getSolver().traceSolver )
	{
		case meq::TraceSolverType::UMFPack: traceSolver = TS::UMFPack; break;
		case meq::TraceSolverType::Pardiso: traceSolver = TS::Pardiso; break;
		case meq::TraceSolverType::cuDSS:   traceSolver = TS::cuDSS;   break;
	}

	if ( !meq::GradShafranovSolver::assemblyModeAvailable( assemblyMode ) )
	{
		/*
		 * ASKED FOR IS REFUSED; INHERITED IS DOWNGRADED. The default is
		 * "threaded" and most builds of MFEM have neither MFEM_USE_OPENMP nor
		 * MFEM_THREAD_SAFE, so refusing unconditionally would make every
		 * example in this repository fail on a stock build -- including the one
		 * CI compiles. A file that says nothing has expressed no preference and
		 * gets the mode that works.
		 *
		 * A file that SAYS "threaded" is a different matter and is refused, for
		 * the reason the trace solver is refused two paragraphs down: a caller
		 * naming a mode has a reason, and silently doing something else would
		 * be invisible in the answer, both modes reaching the same equilibrium.
		 */
		if ( config->getSolver().assemblyModeWasGiven )
		{
			std::fprintf( stderr,
				"MEQ: [solver] AssemblyMode = \"threaded\" needs an MFEM built with\n"
				"     both MFEM_USE_OPENMP and MFEM_THREAD_SAFE, and this one has at\n"
				"     least one of them off. MFEM aborts rather than falling back to\n"
				"     the serial loop, so this is refused here instead.\n"
				"     Use AssemblyMode = \"serial\", or rebuild MFEM.\n" );
			return ConfigurationError;
		}
		assemblyMode = AM::Serial;
	}

	/*
	 * cuDSS IS REFUSED BY THE DRIVER EVEN WHERE THE BUILD HAS IT, AND THE REASON
	 * IS NOT THAT THE DRIVER CONFIGURES NO mfem::Device. IT IS THAT CONFIGURING
	 * ONE WOULD BE THE WRONG TRADE UNTIL THE REST OF THE SOLVE GOES WITH IT.
	 *
	 * A device solver is only worth having if the data STAYS on the device.
	 * MFEM's own plan for this -- doc/HDG-DEVICE-OFFLOAD.md on the branch MEQ
	 * builds from, which is explicitly under construction -- divides the
	 * element-local work into four groups and measures their shares:
	 *
	 *     group 1  local dense linear algebra     7-10% of an NPC step
	 *     group 2  the integrators               46-53%
	 *     group 3  the scatter into the matrix   12-17%
	 *     group 4  THE TRACE SOLVE               26-31%   <- this is cuDSS
	 *
	 * and it gates the whole thing on group 2, in its own words: "Two of the
	 * four groups are nearly free and doing only those is worse than doing
	 * nothing. Groups 1 and 4 leave the integrators on the host, so every
	 * iteration would copy the local blocks host-device around host-side
	 * integrator work -- plausibly slower than staying on the host throughout."
	 *
	 * Group 2 needs a partial-assembly rewrite and is not built. So exposing
	 * cuDSS from a config file today is precisely the group-4-alone case: MEQ
	 * would configure a Device, every Vector in the process would allocate
	 * through it, the integrators and the scatter -- 58% to 70% of an NPC step
	 * -- would still run on the host, and each Newton iteration would pay a
	 * round trip for the one part that moved.
	 *
	 * THERE IS ALSO AN IMMEDIATE FAILURE, and it is what made the refusal
	 * urgent rather than merely principled. Without a Device, CuDSSSolver does
	 * not fall back: it reads its matrix through SparseMatrix::ReadI/ReadJ/
	 * ReadData and its vectors through Read()/Write(), which hand back HOST
	 * pointers, and it aborts inside CUDA with
	 * `cudaMemcpyDeviceToDevice ... invalid argument` -- a message with nothing
	 * in it about the key that caused it.
	 *
	 * The LIBRARY still offers cuDSS and should: it is how correctness on the
	 * device path is checked at all. tests/performance/TraceSolverScaling.cpp
	 * constructs the Device first, and the `cuDSSTraceSolver` ctest pins cuDSS
	 * against UMFPACK. What is refused here is only the config-file route, and
	 * only until the field data has a reason to be on the device.
	 */
	if ( traceSolver == TS::cuDSS )
	{
		std::fprintf( stderr,
			"MEQ: [solver] TraceSolver = \"cudss\" is not available through this\n"
			"     driver, and it is withheld rather than merely unimplemented.\n"
			"\n"
			"     A device solver is only worth having if the data stays on the\n"
			"     device. MEQ's element-local integrators and its scatter into\n"
			"     the trace matrix -- together most of a Newton step -- have no\n"
			"     device kernels yet, so a device trace solve would copy the\n"
			"     system across the bus once per iteration to accelerate one\n"
			"     part of it. That is slower than staying on the host, which is\n"
			"     what MFEM's own HDG device-offload plan concludes about doing\n"
			"     exactly this group on its own.\n"
			"\n"
			"     Use \"umfpack\" or \"pardiso\". cuDSS stays reachable from the\n"
			"     library, which is how its agreement with UMFPACK is checked.\n" );
		return ConfigurationError;
	}

	if ( !meq::GradShafranovSolver::traceSolverAvailable( traceSolver ) )
	{
		char const *needs =
			( traceSolver == TS::Pardiso ) ? "MFEM_USE_MKL_PARDISO"
			: ( traceSolver == TS::cuDSS ) ? "MFEM_USE_CUDSS" : "MFEM_USE_SUITESPARSE";
		std::fprintf( stderr,
			"MEQ: [solver] TraceSolver names a solver this build does not have;\n"
			"     it needs %s. It is refused rather than\n"
			"     silently replaced, because a caller naming a solver has a\n"
			"     reason -- and all three reach the same equilibrium, so the\n"
			"     substitution would be invisible in the answer.\n", needs );
		return ConfigurationError;
	}

	/*
	 * AND A WARNING ABOUT THE THREAD COUNTS, WHICH IS THE ONE COMBINATION A USER
	 * CAN REACH BY ACCIDENT AND WHICH IS CATASTROPHIC RATHER THAN MERELY SLOW.
	 *
	 * MKL suppresses its own threading inside an active OpenMP parallel region,
	 * so under AssemblyMode::Threaded the element-local dgetrs and dgemm cost
	 * the same whatever MKL_NUM_THREADS says -- which is what makes PARDISO's
	 * threads spendable at all. But a team of ONE thread does not get that
	 * suppression: measured on the isolated kernels, OMP_NUM_THREADS=1 with
	 * MKL_NUM_THREADS=8 cost 12.4 s against 0.069 s serial at k = 3.
	 *
	 * This is a warning and not a refusal, because it is a run-time environment
	 * question rather than a fault in the file, and because a user who has
	 * genuinely set both deliberately should not be stopped. It is printed only
	 * when the combination is actually present, which is the standing rule this
	 * file learned from the MKL_THREADING_LAYER warning it used to carry.
	 */
	if ( assemblyMode == AM::Threaded )
	{
		char const *ompEnv = std::getenv( "OMP_NUM_THREADS" );
		char const *mklEnv = std::getenv( "MKL_NUM_THREADS" );
		int const ompCount = ompEnv ? std::atoi( ompEnv ) : 0;
		int const mklCount = mklEnv ? std::atoi( mklEnv ) : 0;
		if ( ompCount == 1 && mklCount > 1 )
			std::fprintf( stderr,
				"MEQ: warning -- AssemblyMode = \"threaded\" with OMP_NUM_THREADS=1\n"
				"     and MKL_NUM_THREADS=%d. A one-thread OpenMP team does not get\n"
				"     MKL's nested-region suppression, and the element-local dense\n"
				"     work then pays full MKL threading per call: measured 12.4 s\n"
				"     against 0.069 s at k = 3. Either raise OMP_NUM_THREADS or set\n"
				"     MKL_NUM_THREADS=1.\n", mklCount );
	}

	// ---- set the run up ------------------------------------------------
	mfem::Mesh background;
	std::unique_ptr<meq::BoundaryShape> shape;
	/// The exterior coupling of [boundary.exterior]. BORROWED by the solver, so
	/// it is declared out here and outlives every cycle -- and it is the SAME
	/// object at every cycle deliberately: Gamma is fixed while only Gamma_h
	/// climbs toward it, which is what makes the exterior coefficients
	/// comparable across a refinement.
	std::unique_ptr<meq::ExteriorDtN> exterior;
	mfem::PositionFunction levelSet;
	/// Curved: D_h is a strict subset of the background mesh, cut by levelSet.
	/// True for [boundary.shape] AND for [boundary.exterior], which are two ways
	/// of describing a Gamma that is not the mesh boundary.
	bool curved = false;
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

	/*
	 * THE PREVIOUS CYCLE'S ANSWER, KEPT ALIVE ACROSS THE TEARDOWN.
	 *
	 * An adaptive cycle destroys its solver before refining -- it must, since
	 * the mesh moves under it -- so the converged potential dies with the space
	 * it lives on unless it is copied out first. These four are that copy, and
	 * they are four rather than one because a GridFunction is only meaningful
	 * beside its space, its collection and its mesh, all of which are about to
	 * be replaced.
	 *
	 * The collection is rebuilt BY NAME from the solver's own, rather than
	 * assumed: MEQ's volume spaces are on the closed Gauss-Lobatto basis, and a
	 * dof-for-dof copy into a Legendre space of the same degree would be a
	 * different function while looking like a successful copy.
	 */
	std::unique_ptr<mfem::Mesh> previousMesh;
	std::unique_ptr<mfem::FiniteElementCollection> previousCollection;
	std::unique_ptr<mfem::FiniteElementSpace> previousSpace;
	std::unique_ptr<mfem::GridFunction> previousPotential;

	/// The guess actually handed to the next solver, on ITS space.
	/// setInitialGuess() borrows, so this has to outlive the solve.
	std::unique_ptr<mfem::GridFunction> carried;

	meq::AdaptivityConfig const &adapt = config->getAdaptivity();
	// MaxIterations counts SOLVES, not refinements: a run with MaxIterations = 1
	// is a plain single solve with an estimate printed, and MaxIterations = 10
	// refines at most nine times.
	int const maxCycles = adapt.enabled ? adapt.maxIterations : 1;

	try
	{
		background = buildMesh( config->getMesh() );

		/*
		 * A COIL OUTSIDE THE MESH CONTRIBUTES EXACTLY NOTHING, SILENTLY.
		 *
		 * F is assembled by quadrature over the elements, so a coil the mesh
		 * does not reach is never sampled: the run converges, writes its files,
		 * and describes a machine with that coil switched off. Nothing else in
		 * the output can show it -- the current appears in `coil_current`
		 * whether or not it did any work.
		 *
		 * A WARNING AND NOT A REFUSAL, because the configuration is not wrong
		 * in principle: FREE-BOUNDARY-PLAN.md section 5.4 offers exactly this
		 * -- coils outside Omega, entering through the exterior coupling rather
		 * than through F. That route is not wired yet, so today the coil is
		 * inert, and saying so is the honest thing. The test is against the
		 * BACKGROUND box, which is the coarsest true statement on every path:
		 * on the curved path D_h is smaller still, so a coil this check passes
		 * may yet be only partly covered.
		 */
		if ( coils )
		{
			mfem::Vector low, high;
			background.GetBoundingBox( low, high );

			for ( std::size_t i = 0; i < coils->size(); ++i )
			{
				meq::Coil const &one = coils->coil( i );
				bool const overlaps = one.rMax() > low( 0 ) && one.rMin() < high( 0 )
				                   && one.zMax() > low( 1 ) && one.zMin() < high( 1 );
				if ( !overlaps )
					std::fprintf( stderr,
						"MEQ: warning: coil %d spans r [%g, %g], z [%g, %g], which is\n"
						"     entirely outside the mesh box r [%g, %g], z [%g, %g]. Its\n"
						"     current enters F nowhere, so it contributes NOTHING to this\n"
						"     solve.\n",
						static_cast<int>( i ), one.rMin(), one.rMax(),
						one.zMin(), one.zMax(),
						low( 0 ), high( 0 ), low( 1 ), high( 1 ) );
			}
		}

		// The curved path solves on D_h, a SUBSET of the background mesh; the
		// fitted path solves on the background mesh itself. Everything
		// downstream -- the solver, the sampler, the files -- follows this.
		meq::ShapeConfig const &shapeConfig = config->getBoundary().shape;
		meq::ExteriorConfig const &exteriorConfig = config->getBoundary().exterior;
		curved = shapeConfig.type != meq::ShapeType::None || exteriorConfig.given;

		if ( exteriorConfig.given )
		{
			/*
			 * FREE BOUNDARY, AND Gamma IS THE SEMICIRCLE RATHER THAN A FLUX
			 * SURFACE. Config has already refused this alongside
			 * [boundary.shape], so exactly one of the two branches runs.
			 *
			 * THE DOMAIN MUST REACH THE AXIS EXACTLY, and this is where that
			 * is checked rather than discovered. meq::ExteriorDtN is diagonal
			 * because the Gegenbauer separation holds on a semicircle CENTRED
			 * ON THE AXIS; a domain starting at r = 0.05 has a flat side that
			 * is an arbitrary vertical line, the modes do not span its
			 * exterior, and the run would converge at full order to a machine
			 * nobody described. That is FB-A's requirement arriving through
			 * the configuration layer.
			 *
			 * BOTH TESTS ARE AGAINST THE MESH AND NOT AGAINST [mesh] RMin ..
			 * ZMax, AND THAT IS THE WHOLE POINT OF THIS PARAGRAPH. Those keys
			 * describe a box the driver BUILDS; with [mesh] File they are
			 * absent and default to zero, so the axis test passed vacuously on
			 * a mesh nobody had looked at and the radius test compared Gamma
			 * against an rMax of 0 and refused every file outright, with a
			 * message about a box the run does not have. That combination --
			 * a gmsh half-disc with the conductors meshed to, plus the
			 * exterior coupling -- is exactly the one FB-6 needs, and it was
			 * the one path on which neither precondition was enforced.
			 * mfem::Mesh::GetBoundingBox is the coarsest true statement on
			 * every path and reproduces the box keys exactly where they apply.
			 */
			mfem::Vector meshLow, meshHigh;
			background.GetBoundingBox( meshLow, meshHigh );

			if ( meshLow( 0 ) != 0.0 )
				throw std::runtime_error(
					"[boundary.exterior] needs a mesh reaching r = 0 exactly, "
					"and this one starts at r = "
					+ std::to_string( meshLow( 0 ) ) + ": the exterior "
					"expansion is valid only on a semicircle centred on the "
					"axis, and a domain stopping short of r = 0 is not a "
					"slightly worse one -- the Gegenbauer modes do not span "
					"its exterior at all" );

			double const zLow = exteriorConfig.centreZ - exteriorConfig.radius;
			double const zHigh = exteriorConfig.centreZ + exteriorConfig.radius;
			if ( exteriorConfig.radius >= meshHigh( 0 )
			     || zLow <= meshLow( 1 ) || zHigh >= meshHigh( 1 ) )
				throw std::runtime_error(
					"[boundary.exterior] Radius puts Gamma outside or on the "
					"mesh, which spans r [0, "
					+ std::to_string( meshHigh( 0 ) ) + "], z ["
					+ std::to_string( meshLow( 1 ) ) + ", "
					+ std::to_string( meshHigh( 1 ) ) + "]; D_h is cut FROM "
					"that mesh, so Gamma must fit strictly inside it" );

			exterior = std::make_unique<meq::ExteriorDtN>(
				exteriorConfig.centreZ, exteriorConfig.radius,
				exteriorConfig.modes );
			levelSet = semicircleLevelSet( exteriorConfig.radius,
			                               exteriorConfig.centreZ );
		}

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
		}

		if ( curved )
		{
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
						std::string( curved && exterior ? "[boundary.exterior] "
						                                : "[boundary.shape] " )
						+ error.what()
						+ " -- either the [mesh] box is too coarse for the surface, "
						"or the surface is not strictly inside it" );
				}
			}
			else
			{
				subdomain = buildSubdomain( background, levelSet,
				                            backgroundCellSize( config->getMesh(),
				                                                background ) );
				solveMesh = subdomain.mesh.get();
				gammaHMarker = &subdomain.gammaHMarker;
				path = std::move( subdomain.path );
				widened = subdomain.widened;
			}
		}

		if ( !curved )
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

			case meq::InitialGuessType::Bump:
			{
				// A CORE, WHICH IS WHAT SELECTS THE BRANCH. See Config.hpp: a
				// free-boundary problem has more than one converged solution and
				// the guess is what says which of them is wanted. A ramp is
				// antisymmetric in z and describes no plasma; this is a
				// paraboloid, positive inside its ellipse and zero outside.
				meq::InitialGuessConfig const &g = config->getInitialGuess();
				double const amplitude = g.amplitude;
				double const centreR = g.centreR;
				double const centreZ = g.centreZ;
				double const radiusR = g.radiusR;
				double const radiusZ = g.radiusZ;
				ramp = std::make_unique<mfem::FunctionCoefficient>(
					[ amplitude, centreR, centreZ, radiusR, radiusZ ]
					( mfem::Vector const &x )
					{
						double const dr = ( x( 0 ) - centreR )/radiusR;
						double const dz = ( x( 1 ) - centreZ )/radiusZ;
						double const t = 1.0 - ( dr*dr + dz*dz );
						return t > 0.0 ? amplitude*t : 0.0;
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
		std::fprintf( stderr, "MEQ: %s\n", error.what() );
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

		/*
		 * THE OTHER TWO BORDERS. Both are unknowns of the SAME bordered Newton
		 * psi_ax already lives in -- N + 2 borders against one factorisation --
		 * so the order they are set in does not matter and neither costs a
		 * second solve.
		 *
		 * ORDER AGAINST setSource() DOES MATTER, though, and it is why these sit
		 * here rather than beside setExtension() above: setBoundaryFluxPoint()
		 * makes psi_bnd an unknown of a NORMALISATION that has to exist first.
		 */
		meq::LimiterConfig const &limiterConfig = config->getBoundary().limiter;
		if ( limiterConfig.surfaceAttribute > 0 )
			// THE CURVE. psi_bnd = max psi_h over the meshed limiter surface,
			// with the contact found rather than prescribed. Config refuses
			// this beside R/Z, so the two branches are exclusive here by
			// construction rather than by precedence.
			fresh->setLimiterSurface( limiterConfig.surfaceAttribute );
		else if ( limiterConfig.given )
			fresh->setBoundaryFluxPoint( limiterConfig.r, limiterConfig.z );

		/*
		 * THE THIRD BORDER, AND THE FILE SPEAKS AMPERES WHERE THE SOLVER SPEAKS
		 * mu0 I_p. setPlasmaCurrent() takes mu0 I_p deliberately -- everything
		 * inside the solver already does, Ampere's law reads the flux integral
		 * as -mu0 I_p and the constraint is assembled as int F/r which IS
		 * mu0 I_p -- and a solver that took amperes would need a mu0 of its own,
		 * which could disagree with the source's and scale two terms of one
		 * equation differently.
		 *
		 * THE CONFIGURATION LAYER HAS NO SUCH PROBLEM AND SO TAKES THE UNIT A
		 * USER HAS. The file names exactly one mu0, under [source], and it is
		 * the same one [[coils]] uses for its Current -- SourceConfig
		 * ::permeability() exists for precisely that sharing. So the conversion
		 * happens here, once, against the only mu0 in the file.
		 */
		if ( config->getSource().plasmaCurrent() != 0.0 )
			fresh->setPlasmaCurrent( config->getSource().permeability()
			                         *config->getSource().plasmaCurrent() );

		if ( exterior )
			fresh->setExteriorCoupling( *exterior );

		fresh->setBoundaryData( zero );
		fresh->setNewtonControl( config->getSolver().newtonRelativeTolerance,
		                         config->getSolver().newtonAbsoluteTolerance,
		                         config->getSolver().newtonMaxIterations );

		// The two performance keys. Availability was checked once at startup,
		// before any mesh was built -- see checkSolverCapabilities() -- so
		// these cannot throw here, and a run that is going to be refused for
		// its build is refused before it does any work.
		fresh->setAssemblyMode( assemblyMode );
		fresh->setTraceSolver( traceSolver );

		if ( !firstCycle && previousPotential )
		{
			/*
			 * THE INTERPOLATING WARM START, AND THE REASON CYCLES AFTER THE FIRST
			 * NO LONGER BEGIN COLD. The FULL-ORDER route of docs/running.rst, which
			 * this file used to refuse on the grounds that it "needs GSLIB" --
			 * and GSLIB has been on, with meq::FieldTransfer written and pinned
			 * by WarmStartConvergence, for longer than that comment survived.
			 *
			 * The previous cycle converged on a coarser mesh, so its answer is
			 * the best starting point available for this one: same equilibrium,
			 * one refinement apart. The fallback is the Dirichlet datum, which is
			 * what a cold start would have had at a node the old mesh does not
			 * cover -- on the curved path the computational domain GROWS as it
			 * refines, so there genuinely are such nodes.
			 */
			meq::FieldTransfer transfer( *previousMesh );
			carried = std::make_unique<mfem::GridFunction>(
				fresh->potential().FESpace() );
			int const missed = transfer.transfer( *previousPotential, zero, *carried );
			fresh->setInitialGuess( *carried );

			std::printf( "MEQ: warm start from the previous cycle: %d of %d nodes "
			             "interpolated, %d fell outside it\n",
			             transfer.queried() - missed, transfer.queried(), missed );
		}
		else if ( config->getInitialGuess().type == meq::InitialGuessType::Ramp
		          || config->getInitialGuess().type == meq::InitialGuessType::Bump )
		{
			fresh->setInitialGuess( *ramp );
		}
		else if ( config->getInitialGuess().type == meq::InitialGuessType::GridFunction
		          && firstCycle )
		{
			// TWO ROUTES, AND THE MESH DECIDES WHICH. Same mesh, same degree:
			// the EXACT restart, every coefficient. A different mesh: the
			// INTERPOLATING restart of docs/running.rst, through
			// meq::FieldTransfer. This branch used to refuse the second and say
			// it "needs FindPointsGSLIB and is not written yet", which had been
			// untrue for months -- see the else below.
			//
			// It is a FIRST-CYCLE guess because cycles after the first have a
			// better one available: the previous cycle's converged answer,
			// interpolated onto the refined mesh, which is the branch at the top
			// of this lambda. A stored file is what the FIRST solve starts from,
			// and that coarse first solve is the one at risk -- it is what the
			// ladder below is for.
			if ( guessMesh->GetNE() == mesh.GetNE() )
			{
				// The EXACT restart: same mesh, same degree, every coefficient.
				fresh->setInitialGuess( *guess );
			}
			else
			{
				// AND A DIFFERENT MESH IS NO LONGER REFUSED. This threw until
				// 2026-09-02, saying the interpolating restart "needs GSLIB" --
				// which was written before meq::FieldTransfer existed and was
				// never revisited. Restarting from a run at another resolution is
				// the ordinary way to use a stored answer, and it works.
				meq::FieldTransfer transfer( *guessMesh );
				carried = std::make_unique<mfem::GridFunction>(
					fresh->potential().FESpace() );
				int const missed = transfer.transfer( *guess, zero, *carried );
				fresh->setInitialGuess( *carried );

				std::printf( "MEQ: interpolating warm start from %s: %d elements "
				             "to this run's %d, %d of %d nodes fell outside the "
				             "stored mesh\n",
				             config->getInitialGuess().meshFile.c_str(),
				             guessMesh->GetNE(), mesh.GetNE(), missed,
				             transfer.queried() );
			}
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
	//
	// The frames carry psi* at degree k+1, exactly as "<stem>" does, because
	// which FIELD is drawn and which GEOMETRY it is drawn on are separate
	// questions. The geometry stays as solved -- Gamma_h, faceted, unbent --
	// since bending mid-loop would hand the next refinement a domain the
	// estimator never saw. The field is the better one either way, and a series
	// whose last frame disagreed with the answer for no stated reason would be
	// a small trap of its own.
	std::unique_ptr<meq::VtuSeries> series;
	if ( adapt.enabled )
		series = std::make_unique<meq::VtuSeries>(
			config->getOutput().directory + "/" + config->getOutput().prefix,
			config->getDiscretisation().polynomialDegree + 1 );

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
			std::fprintf( stderr, "MEQ: %s\n", error.what() );
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
			 * Every globalisation MEQ has drives a residual of its own: the
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
					"MEQ: Newton did not converge: %s\n"
					"     [source] Normalised = true makes psi_ax an unknown, closed by a\n"
					"     BORDERED Newton, and no globalisation MEQ has can carry that\n"
					"     border -- so there is no fallback to try. The levers are a\n"
					"     better [source] PsiAxis guess, an [initialguess] that puts psi\n"
					"     near the right size, and RESOLUTION.\n",
					firstAttempt.what() );
				return SolveFailed;
			}

			std::fprintf( stderr,
				"MEQ: Newton did not converge: %s\n"
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
					"MEQ: Picard-then-Newton did not converge either: %s\n"
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
		 * was a different function, 20x to 64x out. MEQ's Newton path puts every
		 * source on the non-linear form, so that reached every configuration the
		 * driver could be given, and this loop ran on Potential::Raw instead: one
		 * order down at every k, correct but blunt.
		 *
		 * The fix landed as "The postprocessing closes on the element average,
		 * always" -- the close is unconditional now, because the local problem is
		 * a pure Neumann one by construction and there was never anything to
		 * decide. Measured from MEQ's side, the same case that read 20.3, 64.1 and
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
		 * THAT WAS AN OMISSION AND IT IS NOW A REPAIR, so the paragraph above is
		 * history rather than current behaviour. eta_5 is rebuilt on the datum
		 * actually imposed -- GradShafranovSolver::transferredDatum(), which is
		 * mfem::TransferredDatumCoefficient -- so those faces are IN, compared
		 * against phi_h instead of against a zero standing in for it. Measured on
		 * the extension benchmark at k = 2, eta_5 drops from 4.07e-1 to 9.6e-5 on
		 * the coarsest mesh and converges at 2.78 rather than 0.40, and eta keeps
		 * its rate. theTransferredDatumRestoresEtaFive has the table.
		 *
		 * IT IS REBUILT EVERY CYCLE ON PURPOSE. The datum is a lifting of the
		 * SOLVED FLUX along the transfer paths, so it goes stale the moment the
		 * solver is solved again -- and in this loop it is, on a new mesh with a
		 * new path family. Hoisting it out of the loop would compare psi* against
		 * the previous cycle's boundary condition, which is exactly the class of
		 * quiet error this whole block is about.
		 */
		std::unique_ptr<mfem::Coefficient> datum;
		if ( gammaHMarker )
		{
			/*
			 * THE DATUM ON THE TRUE GAMMA, AND ON THE COUPLED PATH IT IS NOT
			 * ZERO. transferredDatum()'s default `g` is the zero function, which
			 * is right for every fixed-boundary case in this tree -- Gamma
			 * carries psi = 0 there -- and WRONG the moment [boundary.exterior]
			 * is present, where Gamma carries the Gegenbauer trace
			 * `sum a_n C_n` that setExteriorDatum() deposits as a load.
			 *
			 * LEFT AT ZERO, eta_5 COMPARES psi* AGAINST A CONDITION NOBODY
			 * IMPOSED, and the term does not merely become inaccurate -- it
			 * DIVERGES. eta_5^2 carries an h_e^-1 weight, so an O(1) per-face
			 * mismatch contributes one copy of its square per face and eta_5
			 * grows as sqrt( faces ): measured under near-uniform refinement,
			 * 2.864e-01, 4.262e-01, 6.151e-01, 8.787e-01 over 71, 142, 282, 570
			 * faces, a ratio settling on 1.429 against sqrt( 2 ) = 1.414. Passing
			 * the datum takes eta_5 from 2.864e-01 to 1.554e-04 -- a factor of
			 * 1844 -- and eta from RISING ( 3.100e-01 -> 5.519e-01 ) to falling
			 * ( 1.186e-01 -> 7.530e-02 ) over the same three cycles.
			 *
			 * AND IT WAS CORRUPTING THE MARKING, NOT ONLY THE NUMBER. The loop
			 * spent its budget crowding Gamma_h against an indicator measuring
			 * nothing: 33 and 68 elements marked against 13 and 29 once repaired.
			 * That is the "refines the wrong elements" failure Estimator.hpp
			 * warns about, in production.
			 *
			 * Rebuilt every cycle for the reason above: `a` moves with the
			 * Newton, so a hoisted lambda would carry the previous cycle's trace.
			 */
			if ( exterior )
			{
				meq::ExteriorDtN const *dtn = exterior.get();
				std::vector<double> const *a = &solver->exteriorCoefficients();
				datum = solver->transferredDatum(
					[ dtn, a ]( mfem::Vector const &x )
				{
					double total = 0.0;
					for ( std::size_t m = 0; m < a->size(); ++m )
						total += ( *a )[ m ]*dtn->basis(
							meq::ExteriorDtN::firstMode() + static_cast<int>( m ),
							x( 0 ), x( 1 ) );
					return total;
				} );
			}
			else
			{
				datum = solver->transferredDatum();
			}
			estimator.setTransferredBoundary( *gammaHMarker, datum.get() );
		}
		if ( exterior )
			estimator.setExteriorCoupling( exterior.get() );

		mfem::Vector const &local = estimator.GetLocalErrors();
		record.eta = estimator.GetTotalError();

		mfem::Array<int> marked;
		bool const lastCycle = cycle + 1 == maxCycles;
		bool const reachedTarget = record.eta <= adapt.targetError;

		if ( !lastCycle && !reachedTarget )
		{
			markElements( adapt.strategy, adapt.theta, local, marked );

			/*
			 * A SECOND MARKING PASS ON THE BOUNDARY TERM, AND SUMMING IT INTO
			 * eta IS NOT ENOUGH WITHOUT ONE. eta_6 is the transmission residual
			 * of the exterior coupling and it is the only error the paper's five
			 * terms structurally cannot see -- but it is SMALL against them:
			 * measured on FB-5's loop, 8.58e-04 against an eta of 2.51e-01, so
			 * its share of eta^2 is about 1e-5 and a Doerfler competition never
			 * reaches it. Gamma_h kept all 34 of its faces for four cycles with
			 * the term summed in.
			 *
			 * THAT IS NOT A THRESHOLD TO LOWER. The two are different quantities
			 * in different units -- an interior discretisation error and a
			 * boundary functional -- so one sum over both is a comparison with no
			 * meaning however it is weighted. The boundary term marks on its OWN
			 * distribution and the sets are unioned, which makes the loop drive
			 * both errors down rather than whichever happens to be larger.
			 *
			 * Measured with this in place: Gamma_h refines 34 -> 46 -> 57 -> 64
			 * faces and the exterior coefficients' error falls 1.32e-03 ->
			 * 1.53e-04, where the same loop without it leaves them at 1.32e-03
			 * to five digits.
			 *
			 * eta_6 stays IN eta as well, because the STOPPING rule does have to
			 * see it: a loop that halted on the interior error alone would report
			 * success with the boundary unresolved.
			 */
			if ( exterior )
			{
				mfem::Vector boundary( estimator.localSquares(
					meq::ResidualEstimator::Term::Transmission ) );
				for ( int e = 0; e < boundary.Size(); ++e )
					boundary( e ) = std::sqrt( boundary( e ) );

				mfem::Array<int> boundaryMarked;
				markElements( adapt.strategy, adapt.theta, boundary, boundaryMarked );

				// local.Size() and not solveMesh->GetNE(): the estimator's own
				// element count is what both marking passes indexed into.
				std::vector<char> already( static_cast<std::size_t>(
					local.Size() ), 0 );
				for ( int i = 0; i < marked.Size(); ++i )
					already[ static_cast<std::size_t>( marked[ i ] ) ] = 1;
				for ( int i = 0; i < boundaryMarked.Size(); ++i )
				{
					std::size_t const e =
						static_cast<std::size_t>( boundaryMarked[ i ] );
					if ( !already[ e ] )
					{
						already[ e ] = 1;
						marked.Append( boundaryMarked[ i ] );
					}
				}
				marked.Sort();
			}

			record.marked = marked.Size();
		}
		history.push_back( record );

		// This cycle's frame, written while the solver still owns its fields --
		// a few lines below they are destroyed so the mesh can be refined.
		// postProcess() ran above, for the estimator, so psi* is current here
		// and costs this nothing.
		if ( series )
		{
			mfem::GridFunction cycleField( solver->flux().FESpace() );
			meq::poloidalField( solver->flux(), cycleField );
			series->append( *solveMesh, solver->postProcessedPotential(),
			                cycleField, cycle, static_cast<double>( cycle ) );
		}

		if ( reachedTarget )
		{
			std::printf( "MEQ: eta = %.4e is at or below TargetError = %.4e "
			             "after %d cycle%s\n", record.eta, adapt.targetError,
			             cycle + 1, cycle == 0 ? "" : "s" );
			break;
		}
		if ( lastCycle )
		{
			std::printf( "MEQ: reached MaxIterations = %d with eta = %.4e, above "
			             "TargetError = %.4e. The answer is the finest one "
			             "computed, not a converged one.\n",
			             maxCycles, record.eta, adapt.targetError );
			break;
		}
		if ( marked.Size() == 0 )
		{
			std::printf( "MEQ: the marking strategy selected no elements at "
			             "eta = %.4e, so there is nothing to refine\n", record.eta );
			break;
		}

		/*
		 * THIS CYCLE'S ANSWER IS COPIED OUT BEFORE ANYTHING IS DESTROYED, so the
		 * next cycle can start from it. On the fitted path solveMesh is refined
		 * IN PLACE, so the coarse mesh ceases to exist a few lines below and a
		 * borrowed pointer would dangle rather than merely go stale -- hence a
		 * deep copy of the mesh and not a reference to it.
		 *
		 * Reset in dependency order. previousSpace holds raw pointers to the mesh
		 * and the collection, so replacing either first would leave it pointing
		 * at freed memory even though nothing reads it in between.
		 */
		previousPotential.reset();
		previousSpace.reset();
		previousMesh.reset();
		previousCollection.reset( mfem::FiniteElementCollection::New(
			solver->potential().FESpace()->FEColl()->Name() ) );
		previousMesh = std::make_unique<mfem::Mesh>( *solveMesh );
		previousSpace = std::make_unique<mfem::FiniteElementSpace>(
			previousMesh.get(), previousCollection.get() );
		previousPotential = std::make_unique<mfem::GridFunction>( previousSpace.get() );
		*previousPotential = solver->potential();

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
				"MEQ: cycle %d could not be refined: %s\n"
				"     The previous cycle had converged, but its solver was released\n"
				"     before the mesh changed -- as it must be -- so there is\n"
				"     nothing left to write. Lower [adaptivity] MaxIterations.\n",
				cycle, error.what() );
			return SolveFailed;
		}
	}

	/*
	 * WHETHER psi_ax IS THE FLUX AT A MAGNETIC AXIS, computed in the report block
	 * below and read again by the writer, which is why it is declared out here.
	 * The `.nc` carries the ratio for the same reason it carries `coil_current`:
	 * a consumer differencing two runs cannot otherwise tell a good psi_axis from
	 * a spurious one.
	 */
	meq::AxisAgreement axisCheck;
	bool axisChecked = false;
	// Deferred rather than acted on where it is found, so that the SOURCE
	// check below can speak first: a psi_ax that is not an axis is the
	// CONSEQUENCE of an unbounded axis current and its advice -- look at the
	// guess, look at the mesh -- is wrong when that is the cause.
	bool axisRefused = false;

	{
		// The background element count is the DOMAIN's after any refinement, not
		// this scope's copy: AdaptiveDomain owns and refines its own.
		int const backgroundElements = domain ? domain->numBackground()
		                                      : background.GetNE();
		if ( shape )
			std::printf( "MEQ: curved Gamma: %d of %d background elements inside, "
			             "%d transfer paths widened\n",
			             solveMesh->GetNE(), backgroundElements, widened );

		/*
		 * THE TOTAL CURRENT IS PRINTED BECAUSE IT IS THE ONE COIL NUMBER A
		 * READER CAN CHECK. FB-2's acceptance identity is
		 * `oint ( 1/r ) dpsi/dn dl = -mu0 I`, a property of the trace alone --
		 * so a sign error in a `[[coils]]` block, which is the mistake
		 * FREE-BOUNDARY-PLAN.md section 7 predicts will be made at least once,
		 * shows up here before anything is plotted.
		 */
		if ( coils )
			std::printf( "MEQ: %d coil%s, total current %+.6e A\n",
			             static_cast<int>( coils->size() ),
			             coils->size() == 1 ? "" : "s",
			             coils->totalCurrent() );

		if ( config->getSource().confinesToPlasma() )
			std::printf( "MEQ: the plasma support MOVES: F = 0 wherever the "
			             "normalised flux is non-positive\n" );

		Cycle const &last = history.back();
		std::printf( "MEQ: converged in %d Newton iterations on %d elements, "
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

		/*
		 * THE EXTERIOR COEFFICIENTS ARE AN ANSWER AND NOT A DIAGNOSTIC, so they
		 * are reported rather than left in the solver. On a free-boundary run
		 * they ARE the boundary condition: psi on Gamma is their sum against the
		 * Gegenbauer basis, and everything outside the mesh is determined by
		 * them. A run that printed only psi_ax would be reporting the interior
		 * of a problem whose whole point is what happens outside it.
		 *
		 * AND THEY ARE WHAT AN ADAPTIVE RUN CANNOT SEE. eta estimates the
		 * INTERIOR error; these are a boundary functional, and FB-5 measured
		 * them frozen to five digits across four refinement cycles while eta
		 * fell by a factor of eight. Printing them per run is the cheapest way
		 * for somebody to notice that before trusting a refinement.
		 */
		/*
		 * THE CURRENT IS REPORTED IN AMPERES, WHICH IS WHAT WAS ASKED FOR, and
		 * beside it the SCALE the border solved for -- because a scale far from
		 * one is the statement that the profiles as written carry a very
		 * different current from the one prescribed, which is a modelling fact a
		 * user wants rather than an internal.
		 */
		if ( config->getSource().plasmaCurrent() != 0.0 )
		{
			double const mu0 = config->getSource().permeability();
			std::printf( "     I_p = %.6e A against the %.6e A asked for, "
			             "profile scale %.6e\n",
			             mu0 != 0.0 ? solver->plasmaCurrent()/mu0 : 0.0,
			             config->getSource().plasmaCurrent(),
			             solver->plasmaCurrentScale() );
		}

		if ( exterior )
		{
			std::vector<double> const &a = solver->exteriorCoefficients();
			std::printf( "     Gamma at rho = %g, %d Gegenbauer modes:",
			             config->getBoundary().exterior.radius,
			             static_cast<int>( a.size() ) );
			for ( std::size_t m = 0; m < a.size(); ++m )
				std::printf( " a%d=%.4e", static_cast<int>( m ) + 2, a[ m ] );
			std::printf( "\n" );

			/*
			 * IS `Modes` ENOUGH? THE RAW COEFFICIENTS CANNOT SAY, AND THAT IS
			 * WHY THIS IS PRINTED RATHER THAN LEFT TO THE READER.
			 *
			 * The exterior block is diagonal in the weight `dGamma/r` with a
			 * mass `2/( n( n-1 )( 2n-1 ) )`, so `a_n` carries `n^{3/2}` of its
			 * own before any physics: a spectrum that is genuinely DECAYING
			 * reads FLAT in the printed row above, and one that is flat is
			 * growing. The energy-normalised amplitude `|a_n| sqrt( mass( n ) )`
			 * is what decays, and the ratio of the tail to the largest is the
			 * one number that says whether the truncation is converged. Measured
			 * over ten modes the two views differ by a factor of 19.6.
			 *
			 * meq::ExteriorDtN::truncationRatio() takes the larger of the LAST
			 * TWO amplitudes, for a parity reason recorded there: an up-down
			 * symmetric trace has identically zero odd modes, so reading the
			 * last one alone returns an exact zero -- perfect convergence -- for
			 * half of all mode counts.
			 */
			/*
			 * THE THRESHOLD IS CALIBRATED AND NOT CHOSEN. Swept over
			 * `Modes` on both shipped couplings, against the relative move in
			 * `psi_ax` from the most-resolved run of each:
			 *
			 *   tail     6.8e-01  3.7e-01  1.6e-01  1.0e-01  1.5e-02  1.1e-02
			 *   psi_ax   3.1e-02  9.5e-03  2.7e-03  2.0e-03  4.2e-05  2.9e-04
			 *
			 * so a tail of 1e-1 is about a per cent in `psi_ax` and 1e-2 is
			 * a few parts in ten thousand. A first cut advised above 1e-2 and
			 * fired on examples/limited-tokamak.toml at ten modes -- which
			 * FREE-BOUNDARY-PLAN.md section 11.6 records as CONVERGED, and
			 * which the sweep confirms at 2.9e-04. Advice that fires on the
			 * converged production case is noise.
			 */
			double const tail = exterior->truncationRatio( a );
			std::printf( "     the retained spectrum's tail is %.3e of its "
			             "largest mode%s\n", tail,
			             tail > 1.0e-1
			                 ? "; raise [boundary.exterior] Modes" : "" );
			if ( config->getBoundary().limiter.given )
			{
				// THE CONTACT IS AN OUTPUT ON THE CURVE ROUTE AND AN INPUT ON
				// THE POINT ONE, so it is read from the solver where it was
				// found and echoed from the file where it was given. Printing
				// the configured r and z under SurfaceAttribute would report
				// ( 0, 0 ), which is on the axis and is not a limiter.
				bool const located = solver->limiterContactWasLocated();
				std::printf( "     psi_bnd = %.6e Wb/rad at the limiter "
				             "( %g, %g )%s\n", solver->psiBoundary(),
				             located ? solver->limiterContactR()
				                     : config->getBoundary().limiter.r,
				             located ? solver->limiterContactZ()
				                     : config->getBoundary().limiter.z,
				             located ? ", found on the limiter surface"
				                     : ", as prescribed" );
			}
		}

		/*
		 * AND THAT CONSTRAINT CANNOT TELL AN AXIS FROM A SPIKE, WHICH IS WHY THE
		 * LINE ABOVE IS NOT ENOUGH.
		 *
		 * psi_ax is the largest NODAL value of psi_h -- deliberately, since one
		 * nodal value is one entry of the discrete unknown and the bordered
		 * Newton's row is then exactly -e_j. But nothing in that definition says
		 * the largest nodal value is a MAGNETIC AXIS, and G = psi_ax - max psi_h
		 * is satisfied at machine zero by a spurious nodal spike exactly as it is
		 * by an axis. Measured on a free-boundary machine case at k = 2: 17 Newton
		 * steps, the constraint at 0.000e+00, the prescribed current delivered to
		 * seven figures, and psi_ax reported twenty-nine times too large from a
		 * single dof of one element. Every number the run printed was green.
		 *
		 * meq::CriticalPointFinder locates the axis properly, as a zero of q_h --
		 * a SOLVED field carrying the potential's own order rather than a
		 * derivative of one -- and the normalised flux there must be 1, because
		 * that is what Psi = ( psi - psi_bnd )/( psi_ax - psi_bnd ) means at the
		 * axis. So the comparison is in the units the profiles actually consume.
		 *
		 * A WARNING AND NOT A REFUSAL, on the same footing as the coil-outside-
		 * the-mesh one above: the run converged, the files are worth writing, and
		 * the honest thing is to say what was found. It is also opt-in by
		 * accident of scope -- psi_ax is an answer only where the source is
		 * normalised, so on every other path there is nothing to check.
		 *
		 * AND IT IS UNCONDITIONAL BECAUSE IT IS CHEAP, which was measured rather
		 * than assumed: on examples/rotating-normalised.toml it costs 0.041 s over
		 * 768 elements and 0.214 s over 12,288, against solves of 1.0 s and 36 s.
		 * The sweep is two Newtons on a 2x2 system per element and is linear in
		 * the mesh.
		 */
		if ( normalised )
		{
			try
			{
				meq::CriticalPointFinder finder( *solver );
				axisCheck = finder.checkAxis( solver->psiAxis(),
				                              solver->psiBoundary() );
				axisChecked = true;
			}
			catch ( std::exception const &error )
			{
				// A diagnostic that cannot run is not a run that failed. Say so
				// and carry on to the files, which are the answer.
				std::fflush( stdout );
				std::fprintf( stderr,
					"MEQ: warning: the magnetic axis could not be located, so\n"
					"     psi_ax was not checked against one: %s\n", error.what() );
			}
		}

		// WHICH DEFINITION psi_ax WAS CONSTRAINED BY, because the two differ by
		// O( h ) in position and O( h^2 ) in value and a reader comparing runs
		// needs to know which they have. The fallback is not an error -- it is
		// the annulus branch, where there is no closed surface to be an axis of.
		if ( normalised && solver
		     && solver->axisConstraint()
		        == meq::GradShafranovSolver::AxisConstraint::LocatedAxis )
			std::printf( "     psi_ax is constrained at %s\n",
			             solver->axisWasLocated()
			                 ? "the located magnetic axis, a zero of q_h"
			                 : "the largest NODAL value: NO O-point was reachable" );

		if ( axisChecked && axisCheck.located )
		{
			std::printf( "     the axis, as a zero of q_h: psi = %.6e at "
			             "( %.4f, %.4f ), normalised flux %.4f\n",
			             axisCheck.axis.psi, axisCheck.axis.r, axisCheck.axis.z,
			             axisCheck.normalisedFlux );

			/*
			 * ABOVE 1 IS INFORMATIVE UNDER OPTION 3, WHICH IT IS NOT UNDER THE
			 * NODAL MAXIMUM, AND THAT IS WHY THIS IS HERE RATHER THAN IN THE
			 * LIBRARY GUARD.
			 *
			 * checkAxis() is one sided BY DESIGN: psi_ax is the largest NODAL
			 * value under AxisConstraint::NodalMaximum, and the peak of a
			 * polynomial over a closed element is at least that, so a healthy
			 * field approaches 1 FROM ABOVE and only the low side can indicate a
			 * defect.
			 *
			 * Under AxisConstraint::LocatedAxis psi_ax IS the polynomial value at
			 * an O-point, so a healthy reading is 1.0000 and a HIGH one says
			 * something specific: checkAxis() found a DIFFERENT, higher O-point
			 * than the one the constraint followed. That is the branch-selection
			 * risk this constraint carries -- the search is warm started from the
			 * previous iterate's axis, so it follows one extremum rather than
			 * re-picking the largest each step, and on a field with several it
			 * can settle on one that is not the core. Measured on a half-disc
			 * with a limiter: the constraint held psi_ax = 1.176e-01 while an
			 * O-point at 1.351e-01 sat elsewhere, reading Psi = 1.2043.
			 *
			 * A WARNING AND NOT A REFUSAL. Which O-point is the core is not
			 * something MEQ can settle -- an equilibrium with several is a real
			 * thing, and the initial guess is what chooses the branch. What can
			 * be said is that the choice was not the obvious one, and said with
			 * the number that shows it.
			 */
			if ( axisCheck.agrees && solver
			     && solver->axisConstraint()
			        == meq::GradShafranovSolver::AxisConstraint::LocatedAxis
			     && solver->axisWasLocated()
			     && axisCheck.normalisedFlux > 1.10 )
			{
				std::fflush( stdout );
				std::fprintf( stderr,
					"MEQ: warning: psi_ax was constrained at the magnetic axis, so the\n"
					"     normalised flux there should read 1 -- and it reads %.4f. A\n"
					"     DIFFERENT O-point of q_h, at ( %.4f, %.4f ), carries psi =\n"
					"     %.6e, which is above the one the solve followed. The axis\n"
					"     search is warm started from the previous iterate, so it\n"
					"     follows one extremum rather than re-picking the largest each\n"
					"     step: on a flux with several, the branch is chosen by the\n"
					"     initial guess. If the core is meant to be the other one, that\n"
					"     is what to change.\n",
					axisCheck.normalisedFlux, axisCheck.axis.r, axisCheck.axis.z,
					axisCheck.axis.psi );
			}

			/*
			 * AND AN "AXIS" INSIDE A CONDUCTOR IS NOT ONE, WHICH IS THE ONE
			 * NON-CIRCULAR THING THAT CAN BE SAID HERE.
			 *
			 * Found 2026-09-07 while pushing FB-6. On a coarse mesh with a
			 * consistent set of coil currents the solve CONVERGED -- every
			 * border satisfied, the prescribed current delivered -- to
			 * psi_ax = -1.2327e-01 at ( 1.7516, 0.9000 ), which is a coil
			 * centre, with a profile scale of -1.03. The coils carry negative
			 * current, so each has an O-point of its own field there; the span
			 * came out negative, so AxisConstraint::LocatedAxis went looking
			 * for a MINIMUM and found the coil's.
			 *
			 * NOTHING ELSE IN THIS BLOCK CAN SEE IT. checkAxis() reports the
			 * normalised flux at the located axis, which under that constraint
			 * is 1 BY CONSTRUCTION -- it read 1.0000 on that run -- and
			 * checkAxisSource() asks about r = 0 and was clean. Every test that
			 * refers to the plasma is circular at that point, because Psi = 1
			 * there is what the constraint imposes. The conductor GEOMETRY is
			 * outside all of it, and the driver is where it is known.
			 *
			 * A WARNING, on the precedent of the paragraph above: which O-point
			 * is the core is the guess's choice, and refinement cures this one
			 * (both finer meshes found the plasma). What can be said is that
			 * the answer is inside a piece of copper, and said with the number.
			 */
			if ( solver && coils && solver->axisWasLocated() )
			{
				for ( std::size_t i = 0; i < coils->size(); ++i )
				{
					meq::Coil const &one = coils->coil( i );
					if ( solver->axisR() < one.rMin()
					     || solver->axisR() > one.rMax()
					     || solver->axisZ() < one.zMin()
					     || solver->axisZ() > one.zMax() )
						continue;

					std::fflush( stdout );
					std::fprintf( stderr,
						"MEQ: warning: the located magnetic axis ( %.4f, %.4f ) is\n"
						"     INSIDE conductor %zu, which spans r [ %.4f, %.4f ]\n"
						"     z [ %.4f, %.4f ].  A plasma has no magnetic axis inside a\n"
						"     coil, so psi_ax is a coil's own O-point and not this\n"
						"     equilibrium's -- and every profile is normalised by it.\n"
						"     The normalised-flux check CANNOT see this: psi_ax is\n"
						"     constrained at the located axis, so Psi there reads 1\n"
						"     whatever the axis is.  Look at the sign of the profile\n"
						"     scale and at the initial guess; refinement has been\n"
						"     measured to cure it.\n",
						solver->axisR(), solver->axisZ(), i,
						one.rMin(), one.rMax(), one.zMin(), one.zMax() );
					break;
				}
			}

			if ( !axisCheck.agrees )
			{
				/*
				 * A REFUSAL AND NOT A WARNING, AND THIS WAS ARGUED RATHER THAN
				 * INHERITED -- FREE-BOUNDARY-PLAN.md section 11.4.
				 *
				 * It used to warn, on the coil-outside-the-mesh precedent. That
				 * precedent is about a configuration which is NOT WRONG IN
				 * PRINCIPLE -- section 5.4 of that plan offers exactly it -- and
				 * this is not that. psi_ax is what the profiles are NORMALISED
				 * by, so a psi_ax that is not the flux at a magnetic axis is not
				 * a bad number in one field: it is a different equilibrium, with
				 * the current, the geometry and every profile-derived quantity
				 * downstream of it. Writing three files describing a machine
				 * nobody asked for is worse than writing none.
				 *
				 * AND IT IS ONLY THE POSITIVE DETECTION THAT REFUSES, which is
				 * what the guard's own one-sidedness entitles us to. checkAxis()
				 * is largest-Psi-wins and its sweep is seeded rather than
				 * exhaustive, so it MISSES rather than false-alarms: reaching
				 * here means an O-point was located AND its normalised flux is
				 * far from 1, which is evidence and not an absence of it. The
				 * `located == false` branch below stays a warning for the
				 * mirror-image reason -- a wall-hugging annulus is a real
				 * equilibrium somebody may want to look at, and "no extremum was
				 * FOUND" is not "no extremum EXISTS".
				 */
				std::fflush( stdout );
				std::fprintf( stderr,
					"MEQ: psi_ax = %.6e is the largest NODAL value of psi_h,\n"
					"     at ( %.4f, %.4f ), and that point is NOT a magnetic axis. The\n"
					"     nearest zero of q_h carries psi = %.6e, which is a normalised\n"
					"     flux of %.4f where the axis must read 1, and it sits %.3e away\n"
					"     -- %.1f diameters of its own element. The profiles have been\n"
					"     evaluated over a normalised flux the plasma never reaches, so\n"
					"     this equilibrium is NOT the one [source] describes. Look first\n"
					"     at the initial guess, which chooses the branch, and at whether\n"
					"     [mesh] resolves the plasma. Nothing has been written.\n",
					axisCheck.psiAxis, axisCheck.nodeR, axisCheck.nodeZ,
					axisCheck.axis.psi, axisCheck.normalisedFlux,
					axisCheck.separation, axisCheck.separationInElements );
				axisRefused = true;
			}
		}
		else if ( axisChecked )
		{
			std::fflush( stdout );
			std::fprintf( stderr,
				"MEQ: warning: no interior extremum of psi_h was found anywhere on\n"
				"     the mesh, so psi_ax = %.6e is not the flux at a magnetic axis:\n"
				"     this solve carries no closed flux surface around one. A search\n"
				"     for zeros of q_h is seeded rather than exhaustive, so this is\n"
				"     evidence and not proof -- but the branch it usually means is a\n"
				"     wall-hugging annulus, where psi rises monotonically to the\n"
				"     boundary and every constraint the solve imposes is satisfied by\n"
				"     a plasma nobody asked for.\n",
				axisCheck.psiAxis );
		}

		/*
		 * AND IS THE TOROIDAL CURRENT DENSITY BOUNDED ON THE SYMMETRY AXIS?
		 *
		 * The load MEQ assembles is -( F/r, w ), and F/r IS mu_0 j_phi:
		 *
		 *     j_phi  =  r p'( Psi )  +  g g'( Psi ) / ( mu_0 r )
		 *
		 * so a finite current density on the axis requires F( 0, z ) = 0, and
		 * F = mu_0 r^2 p' + g g' leaves only g g' there. Which Psi the axis sits
		 * at decides whether that vanishes: psi( 0, z ) = 0 exactly, so
		 * Psi_axis = -psi_bnd/span, and a run with [boundary.limiter] has
		 * psi_bnd > 0 and therefore evaluates the profiles at NEGATIVE Psi -- in
		 * the vacuum, where the physics is g = const and g g' = 0, and where an
		 * unconfined profile extrapolates and returns a current instead.
		 *
		 * A REFUSAL, for the same reason as the one above and one more: the
		 * DISCRETE load functional is then unbounded -- L2 polynomials are free
		 * to be nonzero at r = 0 where the energy space's members are not -- so
		 * the quadrature is the only thing making the assembly finite and the
		 * answer depends on the RULE rather than on the mesh. It is not a bad
		 * number in one place either: measured, psi_h develops an O( 1 ) layer
		 * along the WHOLE axis, at a value that can exceed the true peak.
		 *
		 * IT CANNOT BE CHECKED AT STARTUP, which is why this costs a solve.
		 * psi_bnd is an UNKNOWN of the bordered Newton, so the Psi the axis sits
		 * at is not known until the solve has closed. A startup heuristic on the
		 * keys -- limiter present, mesh reaching r = 0, ConfineToPlasma absent --
		 * would refuse a profile that vanishes below Psi = 0 by construction,
		 * which is a legitimate configuration. This asks F itself, so it cannot
		 * false-refuse.
		 */
		if ( normalised && solver )
		{
			meq::GradShafranovSolver::AxisSourceCheck const axisSource =
				solver->checkAxisSource();

			/*
			 * THE PLASMA CONTAINS THE SYMMETRY AXIS, WHICH IS NOT A LARGE ERROR
			 * BUT THE WRONG TOPOLOGY -- AND IT IS REPORTED FIRST BECAUSE IT IS
			 * THE WORSE OF THE TWO.
			 *
			 * psi( 0, z ) = 0 exactly, so Psi on the axis is -psi_bnd/span, and
			 * a POSITIVE reading needs psi_bnd and the span to carry opposite
			 * signs -- an ordinary positive span with a NEGATIVE psi_bnd. The
			 * plasma then contains r = 0: current threading the machine's own
			 * centre line. A tokamak's axis is in the vacuum by construction, so
			 * no refinement turns this into the equilibrium that was asked for,
			 * and reporting it as a converged answer would be reporting a
			 * different device.
			 *
			 * THE ONE ESCAPE IS A g THAT VANISHES EXACTLY. F( 0, z, . ) is g g'
			 * and nothing else, p' being killed by its own r^2, so g g' == 0
			 * identically leaves j_phi = r p' vanishing with r whatever Psi
			 * reads on the axis. That is the only configuration in which this is
			 * a curiosity rather than a defect, and it is tested over a spread
			 * of Psi rather than at one value.
			 */
			if ( axisSource.reachesAxis && axisSource.axisInsidePlasma
			     && !axisSource.sourceVanishesOnAxis )
			{
				std::fflush( stdout );
				std::fprintf( stderr,
					"MEQ: THE PLASMA CONTAINS THE SYMMETRY AXIS. This is not a\n"
					"     large error, it is the wrong topology, and the run is\n"
					"     refused rather than reported.\n"
					"\n"
					"     psi = 0 on r = 0 exactly -- it is the poloidal flux through a\n"
					"     circle of vanishing area -- so the normalised flux there is\n"
					"     -psi_bnd/( psi_ax - psi_bnd ) = %+.4e, and the plasma is\n"
					"     wherever that is POSITIVE. It is. psi_bnd came out %+.6e\n"
					"     against a span of %+.6e: the two carry opposite signs, which\n"
					"     is what puts r = 0 inside the plasma.\n"
					"\n"
					"     A TOKAMAK is a torus about R_0 > 0 and its symmetry axis is\n"
					"     in the vacuum, so what this describes is toroidal current\n"
					"     threading the machine's own centre line: | F | reads %.6e\n"
					"     there, and F/r is mu_0 j_phi.\n"
					"\n"
					"     A LEVITATED DIPOLE OR A MAGNETIC MIRROR GENUINELY REACHES THE\n"
					"     AXIS, and this is not a refusal of those -- but neither has a\n"
					"     toroidal field, so g vanishes identically in both. That is the\n"
					"     same fact twice: B_phi = g/r must be finite on the axis, so a\n"
					"     plasma that reaches r = 0 cannot carry a toroidal field there.\n"
					"     F( 0, z, . ) is g g' and nothing else, so g g' == 0 leaves\n"
					"     j_phi = r p' going to zero with r whatever Psi reads. THIS\n"
					"     source's does not: checked across Psi = -0.5 .. 1.5, F on the\n"
					"     axis is non-zero. If a dipole or a mirror is what you meant,\n"
					"     set GGPrime to zero and this passes.\n"
					"\n"
					"     WHAT TO LOOK AT, in order: psi_bnd is what went negative, so\n"
					"     [boundary.limiter] and the initial guess, which together\n"
					"     choose the branch; then whether [source] confines to the\n"
					"     plasma at all. Nothing has been written.\n",
					axisSource.normalisedFluxOnAxis, solver->psiBoundary(),
					solver->psiAxis() - solver->psiBoundary(),
					axisSource.worstOnAxis );
				return ConfigurationError;
			}

			if ( axisSource.reachesAxis && !axisSource.bounded )
			{
				std::fflush( stdout );
				std::fprintf( stderr,
					"MEQ: the source does not vanish on the symmetry axis: | F | is\n"
					"     %.6e at ( %.4f, %.4f ), which is %.3e of | F |'s own scale\n"
					"     over the mesh. F/r is mu_0 j_phi, so this is an UNBOUNDED\n"
					"     toroidal current density on r = 0, and the load ( F/r, w ) is\n"
					"     not integrable against a discrete w that does not vanish there\n"
					"     -- psi_h picks up a layer along the whole axis whose size is\n"
					"     set by the quadrature rule rather than by the mesh.\n"
					"     WHAT PUTS IT THERE: psi = 0 on the axis exactly, so the\n"
					"     profiles are being evaluated at Psi = %.4e -- OUTSIDE the\n"
					"     plasma, where the vacuum carries g = const and g g' must be\n"
					"     zero. Set [source] ConfineToPlasma = true, which is what says\n"
					"     the vacuum carries no current, or give a GGPrime that vanishes\n"
					"     for Psi <= 0. Nothing has been written.\n",
					axisSource.worstOnAxis, axisSource.worstR, axisSource.worstZ,
					axisSource.relative,
					solver->psiAxis() > solver->psiBoundary()
					|| solver->psiAxis() < solver->psiBoundary()
						? -solver->psiBoundary()
						  /( solver->psiAxis() - solver->psiBoundary() )
						: 0.0 );
				return ConfigurationError;
			}
		}

		// The deferred one, now that the cause has had its chance.
		if ( axisRefused )
			return ConfigurationError;

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

		/*
		 * psi* IS WHAT MEQ REPORTS, SO IT IS POST-PROCESSED ON EVERY RUN AND NOT
		 * ONLY ON AN ADAPTIVE ONE.
		 *
		 * psi_h lives in P_k; psi* is the element-local post-processing of it in
		 * P_(k+1), and it converges at k+2 rather than k+1 -- measured 47x, 113x
		 * and 125x smaller than psi_h on the finest mesh at k = 1, 2, 3. It costs
		 * one small dense solve per element, once, against a Newton iteration
		 * that has already done far more, so there is no configuration in which
		 * reporting the worse field is the right trade.
		 *
		 * THE GUARD IS NOT DEFENSIVE, IT IS ARITHMETIC. The adaptive loop has
		 * already called postProcess() on this very solver -- eq (20) builds four
		 * of its five terms on psi* -- and every exit from that loop is after the
		 * call, so on that path this would otherwise be the second reconstruction
		 * of the same state.
		 *
		 * AND THE HAZARD, WHICH IS NEW IN KIND RATHER THAN IN DEGREE. psi* comes
		 * from MFEM's DarcyForm::Reconstruct(). The local problem is a pure
		 * Neumann one and its mean-value constraint is what regularises it; that
		 * regularisation was once SKIPPED wherever a non-linear potential
		 * integrator was merely present, so on any element where dF/dpsi vanishes
		 * -- which a tabulated profile with a flat segment produces routinely --
		 * the local matrix was factored singular and psi* there was a different
		 * function, 20x to 64x out. It is PER ELEMENT, so no whole-domain norm
		 * can see it. The fix is "The postprocessing closes on the element
		 * average, always", and it is on the gf-hdg-dev branch and on NO OTHER:
		 * an MFEM built without it loses this silently. See INSTALL.md for which
		 * branches MEQ-integration must contain, and CLAUDE.md's "Post-processing
		 * is back" for the measurement.
		 *
		 * Until today that defect could only reach the error estimator, and a
		 * blunt estimator refines the wrong elements. It now reaches the .nc and
		 * the .vtu -- the primary outputs -- so the regression that catches it,
		 * NewtonConvergence.cpp's
		 * thePostProcessedPotentialIsCorrectWhereTheJacobianVanishes, is guarding
		 * the answer MEQ reports rather than the mesh it refines.
		 */
		if ( !solver->isPostProcessed() )
			solver->postProcess();

		// psi_h AND NOT psi*, AND THAT ASYMMETRY IS DELIBERATE. These three files
		// are the exact restart format: "<stem>_psi.gf" is read back into the
		// degree-k potential space of a solver built from the same configuration,
		// and psi* is degree k+1 -- it does not fit, and a "consistency" edit
		// putting it here would break restart rather than improve it. psi* is
		// written beside them, in its own file, on the next line.
		meq::writeMfem( stem, *solveMesh, solver->potential(), solver->flux() );
		meq::writePostProcessed( stem, solver->postProcessedPotential() );

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
		else if ( config->getMesh().fromFile() )
		{
			/*
			 * A MESH FROM A FILE CARRIES NO RMin..ZMax, SO THE GRID TAKES THE
			 * MESH'S OWN BOUNDING BOX.
			 *
			 * Those four keys describe the box `[mesh]` would have BUILT, and a
			 * file supplies none of them -- they default to zero, the extent
			 * comes out empty, and meq::GridSampler refuses with "the grid
			 * extent must be positive in both directions". The run had already
			 * solved by then, so the whole answer was lost to the output stage.
			 *
			 * This is the ordinary case for free boundary rather than a corner:
			 * the half-disc reaching the axis is not a box and cannot come from
			 * `MakeCartesian2D`, so a free-boundary run ALWAYS reads its mesh
			 * from a file -- see tools/mesh/README.md. The mesh's own extent is
			 * also the right answer, being exactly the region the solve claims
			 * anything about.
			 */
			mfem::Vector low, high;
			solveMesh->GetBoundingBox( low, high );
			gridRMin = low( 0 ); gridRMax = high( 0 );
			gridZMin = low( 1 ); gridZMax = high( 1 );
		}

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
		//
		// THE POTENTIAL SAMPLED HERE IS psi*, NOT psi_h: the file reports the
		// better field. The FLUX ARGUMENT STAYS solver->flux(), which is q at
		// degree k, and that pairing is consistent rather than a leftover --
		// psi*'s defining local problem matches its gradient to r q, so r q IS
		// grad psi* to the order psi* is claimed at. Substituting the enriched
		// flux would change nothing the band can resolve and would pair the
		// Taylor step with a field the reconstruction produced rather than with
		// the one it was built from.
		sampler.samplePotentialWithFlux( solver->postProcessedPotential(),
		                                 solver->flux(), psi, outside );

		// AND B GETS THE BAND TOO, WHICH IT DID NOT UNTIL 2026-09-02. It used to
		// go through sampleComponent(), which reads the FOOT on Gamma_h and so
		// left about one node in ten piecewise constant -- O( h ), in the field
		// most readers open this file for, behind a mask saying inside = 1.
		// sampleComponentWithGradient() takes the same Taylor step psi gets,
		// using B's own gradient inside the element.
		//
		// IT IS NOT AS GOOD AS psi's AND THE DIFFERENCE IS STRUCTURAL: psi is
		// continued with a SOLVED variable at full order, B with a
		// DIFFERENTIATED one at k-1, so this is O( h^2 ) at every k rather than
		// the flux's own order. A full order better than what it replaces, and
		// not the same thing. Sampler.hpp says what it would take to beat it.
		sampler.sampleComponentWithGradient( field, 0, bR, outside );
		sampler.sampleComponentWithGradient( field, 1, bZ, outside );

		/*
		 * THE ROTATING FIELDS: n_s( R, Z ) AND e phi_0( R, Z ).
		 *
		 * With rotation the density is NOT a flux function -- centrifugal force
		 * sweeps the heavy species outboard and an electrostatic potential
		 * arises to stop that separating the charges -- so n_s and phi_0 are
		 * genuine two-dimensional fields and a rotating equilibrium is not
		 * interpretable without them. Both are ALGEBRAIC in ( r, psi ), so this
		 * costs one pass over the grid and no solve. See docs/rotation.rst under
		 * Output, and meq::RotatingSource.
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
			// Through the plasma handle when the coils wrapped the source; see
			// where plasmaSource is declared for why the cast cannot go
			// through `source` in that case.
			meq::Source const *rotatingHandle =
				plasmaSource ? plasmaSource.get() : source.get();
			meq::RotatingSource const *plain =
				dynamic_cast<meq::RotatingSource const *>( rotatingHandle );
			meq::NormalisedRotatingSource const *scaled =
				dynamic_cast<meq::NormalisedRotatingSource const *>( rotatingHandle );

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
		writer.attribute( "title", "MEQ equilibrium" );
		writer.attribute( "meq_version", MEQ_VERSION );
		// The file says which input produced it and how well. A directory of
		// scan output is unreadable otherwise, and "which commit was this?" is
		// the first question asked of any result that looks wrong.
		writer.attribute( "config_file", argument );
		// WHICH PHYSICS, not just which discretisation. This was specified and
		// was not being written: two runs differing only in [source] Type
		// produced files that were identical in every attribute. The variables
		// and attributes are documented in docs/output.rst.
		writer.attribute( "source_type", sourceTypeName( config->getSource().type ) );
		writer.attribute( "polynomial_degree",
		                  config->getDiscretisation().polynomialDegree );
		/*
		 * WHICH POTENTIAL THE `psi` VARIABLE ACTUALLY IS, because there are two
		 * and a reader cannot tell them apart from the numbers.
		 *
		 * psi_h is the solved potential in P_k. psi* is the element-local
		 * post-processing of it in P_(k+1), converging at k+2 -- a different
		 * field, agreeing with psi_h only to psi_h's own accuracy. Files written
		 * before 2026-09-02 carry psi_h and no attribute at all, so the ABSENCE
		 * of this is itself informative; differencing such a file against a new
		 * one at the same resolution measures the post-processing and not the
		 * physics.
		 *
		 * polynomial_degree above stays k, the degree of the SOLVE, which is what
		 * every other attribute here is about. The degree is repeated in this
		 * value so the two cannot be read as contradicting each other.
		 */
		writer.attribute( "potential",
		                  "post-processed (psi*, degree "
		                      + std::to_string(
		                            config->getDiscretisation().polynomialDegree + 1 )
		                      + ")" );
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
		/*
		 * AND WHETHER psi_axis IS A MAGNETIC AXIS, which normalisation_residual
		 * structurally cannot say: it is satisfied at machine zero by a spurious
		 * nodal spike, since psi_ax IS the largest nodal value by definition.
		 *
		 * `axis_normalised_flux` is Psi at the located O-point and MUST BE 1. It
		 * is written for the reason `extrapolated_nodes` had to become a mask:
		 * this is the interchange format, a consumer differencing two runs reads
		 * `psi_axis` and has nothing else to judge it by, and a file that carries
		 * a number without carrying whether it means anything is the shape of
		 * defect this tree keeps finding. `axis_r` and `axis_z` come with it
		 * because "where" is the next question and they cost nothing.
		 *
		 * Their ABSENCE is informative too: a run whose axis could not be located
		 * writes none of the three, which is the annulus branch and the warning
		 * on stderr above.
		 */
		if ( axisChecked && axisCheck.located )
		{
			writer.attribute( "axis_normalised_flux", axisCheck.normalisedFlux );
			writer.attribute( "axis_r", axisCheck.axis.r );
			writer.attribute( "axis_z", axisCheck.axis.z );
		}
		/*
		 * THE COILS ARE PART OF F, SO THE FILE HAS TO SAY SO.
		 *
		 * Without this a reader differencing two `.nc` files -- which is what
		 * the freegs4e benchmark does -- has no way to tell a run with coils
		 * from one without, and the difference between them is the whole of a
		 * free-boundary comparison. `coil_current` is the signed total, which
		 * is the number Ampere's law over the domain checks the outward flux
		 * against; see meq::CoilSet::totalCurrent.
		 */
		if ( coils )
		{
			writer.attribute( "coils", static_cast<int>( coils->size() ) );
			writer.attribute( "coil_current", coils->totalCurrent() );
		}
		if ( config->getSource().confinesToPlasma() )
			writer.attribute( "plasma_support", "moving (F = 0 where Psi <= 0)" );

		// FREE BOUNDARY, and the file has to say so: without these a reader
		// differencing two .nc files cannot tell a truncated vacuum from a
		// prescribed datum on the same Gamma, and the coefficients are the
		// exterior solution in full -- psi outside the mesh is their sum against
		// the Gegenbauer basis and nothing else.
		if ( exterior )
		{
			writer.attribute( "exterior_coupling", "Gegenbauer DtN on a semicircle about the axis" );
			writer.attribute( "exterior_radius", config->getBoundary().exterior.radius );
			writer.attribute( "exterior_centre_z", config->getBoundary().exterior.centreZ );
			writer.attribute( "exterior_modes",
			                  static_cast<int>( solver->exteriorCoefficients().size() ) );
			std::vector<double> const &a = solver->exteriorCoefficients();
			for ( std::size_t m = 0; m < a.size(); ++m )
				writer.attribute( "exterior_a" + std::to_string( m + 2 ), a[ m ] );
		}
		if ( config->getSource().plasmaCurrent() != 0.0 )
		{
			double const mu0 = config->getSource().permeability();
			writer.attribute( "plasma_current_target", config->getSource().plasmaCurrent() );
			writer.attribute( "plasma_current",
			                  mu0 != 0.0 ? solver->plasmaCurrent()/mu0 : 0.0 );
			writer.attribute( "profile_scale", solver->plasmaCurrentScale() );
		}
		if ( config->getBoundary().limiter.given )
		{
			// Where the contact ACTUALLY was, which on the curve route is a
			// solved quantity and not the configuration echoed back. A reader
			// differencing two runs needs the contact that produced psi_bnd.
			bool const located = solver->limiterContactWasLocated();
			writer.attribute( "limiter_r",
			                  located ? solver->limiterContactR()
			                          : config->getBoundary().limiter.r );
			writer.attribute( "limiter_z",
			                  located ? solver->limiterContactZ()
			                          : config->getBoundary().limiter.z );
			writer.attribute( "limiter_contact_located", located ? 1 : 0 );
			writer.attribute( "psi_boundary", solver->psiBoundary() );
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
		// THE BOUNDARY MEQ WAS GIVEN, so a plot can draw the answer against
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
						"MEQ: warning: the mesh boundary is not a single loop; "
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

		/*
		 * ---- the ( Psi, theta ) flux-surface file, INVERSION-PLAN.md IN-6 ----
		 *
		 * The fourth output format, and the one a 1-D transport code reads. The
		 * other three are a field on a mesh, a picture, and a field on a
		 * rasterised ( R, Z ) lattice; none of them is a set of scalar functions
		 * of a flux label, and reconstructing those at the far end means redoing
		 * the whole of the inversion item in the consumer. tools/README.md says
		 * which reader goes with which.
		 *
		 * BEFORE THE VTK, AND THAT ORDER IS LOAD BEARING. The block below bends
		 * the boundary out onto the true Gamma, which changes the map from
		 * reference to physical space -- and the tracer reads geometry at every
		 * corrector step. Extracting after the bend would trace contours of a
		 * field on a mesh that is no longer the mesh it was solved on. This is
		 * the same reason the comment below gives for the sampler, one consumer
		 * further along.
		 *
		 * A WARNING RATHER THAN AN EXIT WHEN IT FAILS, and the reason is what
		 * the file is. The equilibrium has already been solved, checked and
		 * written; this is a derived product, and there are configurations that
		 * solve perfectly well and have no closed flux surfaces at the levels
		 * asked for -- an annulus, a level outside the plasma, a surface that
		 * leaves the mesh. meq refuses to invent one, per MANTA-COUPLING.md
		 * section 8, and refusing produces no file rather than a plausible one.
		 * Losing the whole run to it would be the wrong trade, so this says so
		 * on stderr and the exit code stays as the solve left it.
		 */
		bool wroteFluxSurfaces = false;
		if ( output.fluxSurfaces )
		{
			try
			{
				meq::ContourTracer tracer( *solver );

				// THE BAND, ON THE CURVED PATH ONLY. Omega_h is inscribed in
				// Gamma, so the outermost surfaces of the family are partly
				// outside the mesh; BandExtension::TransferLift continues the
				// field along the same transfer paths the solve imposed its
				// boundary condition through. Measured at k+2 against the flux
				// Taylor step's second order -- see FluxSurfaces.hpp -- and it
				// is what lets the outer cut sit anywhere near Gamma at all.
				// Every node it answers for is marked, per node, in the file.
				if ( gammaHMarker && path )
					tracer.setBandExtension( meq::BandExtension::TransferLift,
					                         *gammaHMarker, path.get() );

				meq::CriticalPointFinder const finder( *solver );
				meq::CriticalPoint const axis = finder.findAxis();

				meq::FluxFamilyOptions options;
				options.surfaces =
					static_cast<std::size_t>( output.fluxSurfaceCount );
				options.angles =
					static_cast<std::size_t>( output.fluxAngleCount );
				options.innerCut = output.fluxInnerCut;
				options.outerCut = output.fluxOuterCut;

				// NO g( psi ), SO NO SAFETY FACTOR COLUMN. g = R B_toroidal is
				// not something a meq::Source carries -- it carries g g' -- and
				// a column of zeroes would be indistinguishable from a machine
				// with no toroidal field. The file says which by not having the
				// variable at all.
				meq::FluxSurfaceFamily const family = meq::extractFluxSurfaces(
					tracer, axis, solver->psiBoundary(), options );

				meq::FluxGridWriter surfaces( output.getFluxSurfaceFile(),
				                              family );
				surfaces.attribute( "title", "MEQ flux-surface geometry" );
				surfaces.attribute( "meq_version", MEQ_VERSION );
				surfaces.attribute( "config_file", argument );
				surfaces.attribute( "source_type",
				                    sourceTypeName( config->getSource().type ) );
				surfaces.attribute( "polynomial_degree",
				                    config->getDiscretisation().polynomialDegree );
				surfaces.attribute( "potential", "post-processed" );
				surfaces.attribute( "band_extension",
				                    gammaHMarker && path ? "transfer lift"
				                                         : "none (fitted)" );
				surfaces.close();

				wroteFluxSurfaces = true;
				std::printf( "MEQ: %zu flux surfaces over Psi_N in [ %.3f, "
				             "%.3f ], %zu nodes each; %d of %zu nodes are band "
				             "data, worst | psi_h - level | %.3e\n",
				             family.size(), family.innerCut, family.outerCut,
				             family.angles, family.extendedNodes(),
				             family.size()*family.angles,
				             family.worstResidual() );
			}
			catch ( std::exception const &error )
			{
				std::fprintf( stderr,
					"MEQ: warning: no flux-surface file was written: %s\n"
					"MEQ:   the equilibrium itself is unaffected and the other "
					"outputs are complete.\n"
					"MEQ:   Narrow [output] FluxInnerCut/FluxOuterCut, or look "
					"at the mesh where it gave out.\n",
					error.what() );
			}
		}

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
				std::printf( "MEQ: %.0f%% of the boundary reached Gamma; the "
				             "rest was held back to keep its elements from "
				             "folding\n", 100.0*curved );
			else if ( curvedNodes == 0 )
				std::fprintf( stderr,
					"MEQ: warning: could not bend the boundary onto Gamma "
					"without tangling an element. The VTK boundary is the "
					"polygon Gamma_h.\n" );
		}

		// psi* AGAIN, AND THE SUBDIVISION GOES UP WITH IT. writeVtu()'s last
		// argument is the Lagrange-cell subdivision, and it must be the degree of
		// the field being drawn: psi* is P_(k+1), so passing k would sample a
		// degree-(k+1) field at degree-k nodes and throw away exactly the order
		// the post-processing was done to gain. The picture would look like a
		// slightly coarse mesh, which is the failure mode SetHighOrderOutput()
		// exists to prevent, one degree further along.
		meq::writeVtu( stem, *solveMesh, solver->postProcessedPotential(), field,
		               config->getDiscretisation().polynomialDegree + 1 );

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
			"MEQ: wrote\n"
			"  exact, for GLVis and restart:  %s.mesh\n"
			"                                 %s_psi.gf        (psi_h, degree %d)\n"
			"                                 %s_grad_psi.gf\n"
			"  post-processed, GLVis:         %s_psistar.gf    (psi*, degree %d)\n"
			"  VTK at degree %d, ParaView:     %s/%s.pvd\n"
			"  (R, Z) grid, %d/%d inside:  %s.nc\n"
			"  psi in the last two is psi*; _psi.gf keeps psi_h, which is what a\n"
			"  restart reads back.\n",
			stem.c_str(),
			stem.c_str(), config->getDiscretisation().polynomialDegree,
			stem.c_str(),
			stem.c_str(), config->getDiscretisation().polynomialDegree + 1,
			config->getDiscretisation().polynomialDegree + 1,
			stem.c_str(), name.c_str(),
			sampler.locatedCount(), sampler.nodesR()*sampler.nodesZ(),
			stem.c_str() );
		if ( series && series->frames() > 0 )
			std::printf( "  %d refinement frames, scrubbable:  %s_cycles/%s_cycles.pvd\n",
			             series->frames(), stem.c_str(), name.c_str() );
		// Named on its own line rather than folded into the block above,
		// because it is the one output that a run may legitimately not have:
		// see the extraction block for why a failure there is a warning.
		if ( wroteFluxSurfaces )
			std::printf( "  (Psi, theta) surfaces, 1-D transport:  %s_surfaces.nc\n",
			             stem.c_str() );
	}
	catch ( std::exception const &error )
	{
		std::fprintf( stderr, "MEQ: could not write output: %s\n", error.what() );
		return OutputFailed;
	}

	return Solved;
}
