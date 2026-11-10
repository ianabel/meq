/*
 * MEQ -- the Maryland Equilibrium Solver.
 *
 *     MEQ config.toml
 *
 * One binary, no subcommands, and deliberately NOT mfem::OptionsParser, which
 * the program this replaces used and which wants to own argument parsing for
 * the whole executable. docs/running.rst is what the driver owes the user.
 *
 * THIS FILE IS THE COMMAND LINE AND NOTHING ELSE. The run -- the
 * configuration read, the run set up, solved and written -- is meq::Run
 * (src/meq/Run.hpp), so that every front end runs the same one: this binary,
 * and the Python module PYTHON-PLAN.md builds. What stays here is what is
 * about the PROCESS: argv, --help, --version, --mesh-command, and the exit
 * code.
 *
 * EXIT CODES, because a shell script driving a parameter scan is a first-class
 * caller and "it printed something" is not an interface. meq::ExitCode has
 * them:
 *
 *     0  solved, and everything asked for was written
 *     1  the configuration is wrong, or the run cannot be set up from it
 *     2  the non-linear solve did not converge -- including a q-driven
 *        run's OUTER loop, whatever its last inner equilibrium did
 *     3  output could not be written
 *     4  an inverse solve was written at the best fit its free currents
 *        allow, and the targets were not met
 *
 * Exit 2 is only reportable because ../mfem/install is built with
 * MFEM_USE_EXCEPTIONS: without it a failed solve takes the process down with
 * SIGABRT before anything here runs. See CLAUDE.md under "On SUNDIALS".
 */

#include "mfem.hpp"

#include "meq_version.hpp"

#include "meq/BuildInfo.hpp"
#include "meq/Coils.hpp"
#include "meq/Config.hpp"
#include "meq/Run.hpp"

#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <optional>
#include <string>
#include <vector>

namespace
{
	using meq::ConfigurationError;
	using meq::Solved;

	void usage()
	{
		std::printf(
			"MEQ -- the Maryland Equilibrium Solver\n"
			"\n"
			"  MEQ <config.toml>    solve the equilibrium the file describes\n"
			"  MEQ --help           this\n"
			"  MEQ --version        the build this is\n"
			"\n"
			"  --mesh-command       print the generator command [mesh.generate]\n"
			"                       describes, shell-quoted, and exit without\n"
			"                       solving. MEQ links MFEM and not gmsh, so it\n"
			"                       says what to run rather than running it.\n"
			"                       Nothing on stdout, and exit 0, when the file\n"
			"                       has no mesh to make.\n"
			"  --mesh-ready         the caller has already made that mesh. A\n"
			"                       configuration carrying [mesh.generate] is\n"
			"                       REFUSED without this, so an edited geometry\n"
			"                       cannot be answered from the previous one's\n"
			"                       mesh. meq-run passes it; you should not.\n"
			"  --profile            print the last solve\'s leg split -- residual,\n"
			"                       gradient, trace factorisation, backsolve,\n"
			"                       plasma fill, and an explicit remainder. The\n"
			"                       timers are always on; this only prints them.\n"
			"                       A leg share is only a measurement on an idle\n"
			"                       machine.\n"
			"  --trace-solver-stats let PARDISO print its own statistics for\n"
			"                       every phase, perturbed pivots and refinement\n"
			"                       steps among them -- what decides the\n"
			"                       backsolve leg. Several lines a phase.\n"
			"  --trace-refinement-cap N\n"
			"                       PARDISO's iterative refinement cap, iparm[7]:\n"
			"                       0 by default, which refines only after a\n"
			"                       perturbed pivot; MFEM's own default is 2. A\n"
			"                       diagnostic: --profile reports the steps each\n"
			"                       solve took and the pivots perturbed.\n"
			"  --q-iteration-cap N  the q-driven loop's Gauss-Newton iteration\n"
			"                       cap, 30 by default. Each iteration costs a\n"
			"                       whole equilibrium or more, m tangents and\n"
			"                       two flux-surface extractions, m the knots,\n"
			"                       so this bounds the run. It cannot\n"
			"                       change the g a converging loop reports, only\n"
			"                       whether it reports one: a loop that has not\n"
			"                       converged within it exits 2. Refused on a run\n"
			"                       without [source] SafetyFactorFile.\n"
			"\n"
			"  --gf                 also write the GLVis files -- <stem>.mesh,\n"
			"                       _psi.gf, _grad_psi.gf, _psistar.gf -- as\n"
			"                       [output] GridFunctions = true does.\n"
			"  --vtk                also write the ParaView files, as [output]\n"
			"                       VTK = true does. Without either only the\n"
			"                       .nc is written: it carries the mesh and\n"
			"                       every coefficient, it is the restart, and\n"
			"                       meq-convert makes both formats from it.\n"
			"\n"
			"  --device <name>      run under an mfem::Device: cpu (the\n"
			"                       default), cuda, or debug. Answers agree with\n"
			"                       the host to every printed digit, and the\n"
			"                       SOLVE leg is at parity with the host-threaded\n"
			"                       default; what a device run still costs is the\n"
			"                       output stage, where the reconstruction gives\n"
			"                       up its host threading. \"debug\" has device\n"
			"                       memory semantics with host arithmetic and\n"
			"                       turns a stale host read into a named fault.\n"
			"\n"
			"Exit codes: 0 solved, 1 configuration, 2 solve did not converge,\n"
			"3 output could not be written, 4 an inverse solve's targets were not\n"
			"met and the best fit was written.\n" );
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

	/**
	 * A double spelled so that reading it back gives the same double, and no
	 * longer than that. `%.17g` always round-trips and always looks wrong --
	 * it writes 0.35 as 0.34999999999999998 -- and this command line is meant
	 * to be read and pasted as well as executed, so the shortest exact
	 * spelling is worth the two extra formats it costs.
	 */
	std::string exactly( double value )
	{
		char buffer[ 64 ];
		for ( int digits = 15; digits <= 17; ++digits )
		{
			std::snprintf( buffer, sizeof buffer, "%.*g", digits, value );
			if ( std::strtod( buffer, nullptr ) == value )
				break;
		}
		return buffer;
	}

	/// One argument of a command line, quoted so a POSIX shell -- and Python's
	/// shlex.split, which is what actually consumes this -- sees exactly the
	/// bytes given. Single quotes, because they quote everything; the only
	/// character that cannot appear inside them is a single quote, which is
	/// spliced.
	std::string shellQuoted( std::string const &word )
	{
		bool safe = !word.empty();
		for ( char const c : word )
			if ( !( std::isalnum( static_cast<unsigned char>( c ) ) || c == '_' || c == '-'
			        || c == '.' || c == '/' || c == '=' || c == '+' || c == ',' || c == ':' ) )
				safe = false;
		if ( safe )
			return word;

		std::string quoted = "'";
		for ( char const c : word )
		{
			if ( c == '\'' )
				quoted += "'\\''";
			else
				quoted += c;
		}
		return quoted + "'";
	}

	/**
	 * THE GENERATOR'S ARGUMENT LIST, BUILT FROM THE CONFIGURATION.
	 *
	 * This is the whole of what the driver knows about meshing, and it is
	 * deliberately a STRING rather than a subprocess: `meq` links MFEM and not
	 * gmsh -- tools/mesh/README.md records why, and the short form is that the
	 * coupling between MEQ and a mesher is a FILE -- so the driver's honest
	 * contribution is to say what to run. `meq-run` runs it.
	 *
	 * TWO CONVERSIONS HAPPEN HERE AND THEY ARE THE REASON THE SCHEMA IS NOT
	 * JUST halfdisc.py's FLAGS SPELLED DIFFERENTLY:
	 *
	 *   * the refined box is four BOUNDS in the file, as [mesh]'s own box is,
	 *     and a corner plus two EXTENTS on the command line;
	 *   * the conductors are CENTRE and HALF-EXTENTS in [[coils]], because
	 *     that is what meq::Coil takes, and corner plus extents here.
	 *
	 * Either convention is fine and having both in one file is not, so the
	 * file gets MEQ's and this function pays the difference. The coils are
	 * emitted in [[coils]] order, which is the order halfdisc.py assigns its
	 * `10 + i` element attributes in -- that agreement is what lets a machine's
	 * conductors be written once.
	 *
	 * **AND A SUBTRACTED CONDUCTOR IS NOT EMITTED AT ALL, WHICH IS THE ONE
	 * PLACE `[conductors] Model` REACHES THE MESH.** The whole claim of the
	 * conductor split is that the mesh need not carry the conductors; a
	 * `--coil` here would make it carry them anyway, and the
	 * generator would fragment and grade a region whose current the solve
	 * evaluates analytically. Every other consequence of the key is inside the
	 * solver, so this loop is the whole of its geometry half -- and
	 * `[mesh.generate] CoilSize`, which exists only to grade around a `--coil`,
	 * is refused at parse under a subtracting model rather than emitted here
	 * and ignored there.
	 *
	 * **THE MESH IS STILL REBUILT WHEN THE KEY CHANGES**, and nothing extra
	 * was needed for that: `meq-run` stamps the mesh with this COMMAND, so
	 * flipping `Model` changes the argv, the stamp misses, and the coil-free
	 * mesh is generated. That is the same mechanism that re-meshes on an
	 * edited coil, and it is why the stamp holds the command rather than the
	 * file's mtime.
	 */
	std::string meshCommand( meq::Configuration const &config )
	{
		meq::MeshGeneratorConfig const &g = config.getMesh().generate;

		std::vector<std::string> words;
		auto word = [ &words ]( std::string const &text ) { words.push_back( text ); };
		auto number = [ &words ]( double value ) { words.push_back( exactly( value ) ); };

		word( g.tool );
		word( "--rho" );    number( g.radius );
		word( "--size" );   number( g.size );
		if ( g.order != 1 )
		{
			word( "--order" );
			word( std::to_string( g.order ) );
		}
		if ( g.coilSize > 0.0 )
		{
			word( "--coil-size" ); number( g.coilSize );
		}
		if ( g.plasmaGiven )
		{
			word( "--plasma" );
			number( g.plasmaRMin );
			number( g.plasmaZMin );
			number( g.plasmaRMax - g.plasmaRMin );
			number( g.plasmaZMax - g.plasmaZMin );
			word( "--plasma-size" ); number( g.plasmaSize );
			// THE BOX AS GEOMETRY: its edges become element faces, which is
			// what keeps a mirror's axial cutoff from capping the order.
			if ( g.alignPlasma )
				word( "--align-plasma" );
		}
		if ( !g.vessel.empty() )
		{
			// EVERY NUMBER AS THE PARSER READ IT, which is the same rule the
			// conductors follow one block down: the mesher fragments this
			// polygon into the geometry, so a rounded vertex is a mesh that does
			// not describe the vessel the file names.
			word( "--vessel" );
			for ( double const value : g.vessel )
				number( value );
		}

		if ( g.limiterGiven )
		{
			word( "--limiter" );
			number( g.limiterR );
			number( g.limiterZ );
			number( g.limiterRadius );
		}
		if ( g.transition > 0.0 )
		{
			word( "--transition" ); number( g.transition );
		}
		// NOT UNDER A SPLIT -- see this function's header. halfdisc.py takes
		// --coil as action="append", so emitting none is already the coil-free
		// mesh and no generator change was needed.
		if ( !config.getConductors().subtracts() )
		{
			for ( meq::CoilParameters const &coil : config.getCoils().coils )
			{
				word( "--coil" );
				number( coil.centreR - coil.halfWidth );
				number( coil.centreZ - coil.halfHeight );
				number( 2.0 * coil.halfWidth );
				number( 2.0 * coil.halfHeight );
			}
		}
		// [mesh.generate] AlignConductors: under a split, the same rectangles
		// as GEOMETRY ONLY -- fragmented in so their edges are element faces,
		// with no attribute and no grading, which is what --align-coil means to
		// the generator. docs/output.rst, What each element is made of; the
		// parse refuses the key under "meshed", where the --coil above
		// already does this.
		else if ( g.alignConductors )
		{
			for ( meq::CoilParameters const &coil : config.getCoils().coils )
			{
				word( "--align-coil" );
				number( coil.centreR - coil.halfWidth );
				number( coil.centreZ - coil.halfHeight );
				number( 2.0 * coil.halfWidth );
				number( 2.0 * coil.halfHeight );
			}
		}
		// AFTER the conductors, because --symmetric checks them: the generator
		// pairs every --coil about z = 0 and refuses an unpaired one, so the
		// flag is a statement about the argv it follows. Under a split there
		// are none to pair, and the flag still means the same thing about the
		// DISC -- halfdisc.py meshes z >= 0 and reflects it whatever the
		// conductors are.
		if ( g.symmetric )
			word( "--symmetric" );
		word( "-o" );
		word( config.getMesh().file );
		if ( g.check )
			word( "--check" );

		std::string line;
		for ( std::string const &one : words )
		{
			if ( !line.empty() )
				line += ' ';
			line += shellQuoted( one );
		}
		return line;
	}
}

int main( int argc, char **argv )
{
	/*
	 * ONE POSITIONAL ARGUMENT -- the configuration file -- AND THREE OPTIONS,
	 * NONE OF WHICH DESCRIBES THE EQUILIBRIUM.
	 *
	 * This was `argc != 2` and a single `argv[1]`, which is the right shape for
	 * a driver whose whole interface is a TOML file and is why the file carries
	 * everything else. The exceptions are exceptions on purpose, and all three
	 * are about the PROCESS rather than about the physics:
	 *
	 *   --device        mfem::Device is PROCESS-WIDE state that has to be
	 *                   configured before the first Vector is allocated, so it
	 *                   cannot be a key in a file that is parsed after MFEM is
	 *                   already running. It is also not a property of the
	 *                   equilibrium -- the same file must describe the same
	 *                   problem on a machine with a GPU and on one without --
	 *                   so putting it in [solver] beside AssemblyMode would
	 *                   make a configuration unportable in a way none of the
	 *                   other keys are.
	 *   --mesh-command  a QUERY about the file, like --version is a query about
	 *                   the build. It does not solve.
	 *   --mesh-ready    a claim by the CALLER about the state of the working
	 *                   directory, which is exactly what a file describing an
	 *                   equilibrium cannot carry.
	 *   --gf, --vtk     which FILES to write beside the .nc, and so about what
	 *                   this caller wants to look at rather than about the
	 *                   physics. They only ever ADD a format: [output]
	 *                   GridFunctions and VTK say the same in the file, and
	 *                   either asking is enough.
	 */
	std::string configPath;
	std::string deviceName;
	bool deviceGiven = false;
	bool wantHelp = false;
	bool wantVersion = false;
	bool wantMeshCommand = false;
	bool meshReady = false;
	bool wantProfile = false;
	bool wantTraceSolverStats = false;
	// Unset leaves the solver's own default.
	std::optional<int> traceRefinementCap;
	// Unset leaves meq::LeastSquaresOptions' own default.
	std::optional<int> qIterationCap;
	bool wantGridFunctions = false;
	bool wantVtk = false;
	bool badUsage = false;

	for ( int i = 1; i < argc; ++i )
	{
		std::string const arg = argv[ i ];
		if ( arg == "--help" || arg == "-h" )
			wantHelp = true;
		else if ( arg == "--version" )
			wantVersion = true;
		else if ( arg == "--mesh-command" )
			wantMeshCommand = true;
		else if ( arg == "--mesh-ready" )
			meshReady = true;
		else if ( arg == "--profile" )
			wantProfile = true;
		else if ( arg == "--trace-solver-stats" )
			wantTraceSolverStats = true;
		else if ( arg == "--trace-refinement-cap" )
		{
			// Any integer is MKL's to read, a negative one included, so only
			// a missing or non-numeric value is refused.
			char *end = nullptr;
			long const value = ( i + 1 < argc )
				? std::strtol( argv[ i + 1 ], &end, 10 ) : 0;
			if ( i + 1 >= argc || end == argv[ i + 1 ] || *end != '\0' )
			{
				std::fprintf( stderr,
					"MEQ: --trace-refinement-cap needs an integer, for example "
					"--trace-refinement-cap 0.\n" );
				return ConfigurationError;
			}
			traceRefinementCap = static_cast<int>( value );
			++i;
		}
		else if ( arg == "--q-iteration-cap" )
		{
			// A count of KINSOL iterations, so it has to be one at least:
			// KINSOL reads zero as "its own default" and a negative as an
			// error, and neither is what a user writing 0 meant.
			char *end = nullptr;
			long const value = ( i + 1 < argc )
				? std::strtol( argv[ i + 1 ], &end, 10 ) : 0;
			if ( i + 1 >= argc || end == argv[ i + 1 ] || *end != '\0'
			     || value < 1 || value > 100000 )
			{
				std::fprintf( stderr,
					"MEQ: --q-iteration-cap needs a positive integer, for example "
					"--q-iteration-cap 10.\n" );
				return ConfigurationError;
			}
			qIterationCap = static_cast<int>( value );
			++i;
		}
		else if ( arg == "--gf" )
			wantGridFunctions = true;
		else if ( arg == "--vtk" )
			wantVtk = true;
		else if ( arg == "--device" )
		{
			if ( i + 1 >= argc )
			{
				std::fprintf( stderr,
					"MEQ: --device needs a name, for example --device cuda.\n"
					"     Names are MFEM's: cpu, cuda, debug, or any backend\n"
					"     string mfem::Device accepts.\n" );
				return ConfigurationError;
			}
			deviceName = argv[ ++i ];
			deviceGiven = true;
		}
		else if ( !arg.empty() && arg[ 0 ] == '-' )
		{
			std::fprintf( stderr, "MEQ: unknown option \"%s\".\n", arg.c_str() );
			badUsage = true;
		}
		else if ( configPath.empty() )
			configPath = arg;
		else
		{
			// REFUSED RATHER THAN IGNORED. Two files named is a mistake with
			// two readings -- which one did they mean? -- and silently solving
			// the first would put the answer to a different problem in the
			// output directory under a name taken from the file that ran.
			std::fprintf( stderr,
				"MEQ: more than one configuration file named (\"%s\" and "
				"\"%s\").\n", configPath.c_str(), arg.c_str() );
			badUsage = true;
		}
	}

	if ( wantHelp )
	{
		usage();
		return Solved;
	}

	if ( wantVersion )
	{
		/*
		 * THE PRECISION AND THE DEVICE BACKENDS BELONG HERE BECAUSE THEY ARE
		 * NOT MEQ'S TO CHOOSE AND ARE NOT IN THE CONFIGURATION FILE.
		 * `MFEM_PRECISION` decides `real_t` throughout the library, so it is a
		 * property of the MFEM this binary linked and is inherited through
		 * `MFEM_DIR`; there is no TOML key for it and there cannot be. So two
		 * installs produce two `meq` binaries that are indistinguishable from
		 * the outside and answer to different floors -- DEVICE-PLAN.md section
		 * 7.5 measures one equilibrium at 2.2 eps and another at 1,260 -- and
		 * this line is the only place a user can ask which one they have.
		 */
		std::printf( "MEQ %s (MFEM %s, %s)\n", MEQ_VERSION, MFEM_VERSION_STRING,
		             meq::buildDescription().c_str() );
		return Solved;
	}

	if ( badUsage || configPath.empty() )
	{
		usage();
		return ConfigurationError;
	}

	std::string const argument = configPath;

	/*
	 * --mesh-command: A QUERY, NOT A SOLVE, AND IT RETURNS BEFORE THE DEVICE.
	 *
	 * It answers "what makes this file's mesh?" out of the configuration
	 * alone, which is MFEM-free -- so it costs a parse and nothing else, and
	 * in particular it does NOT need the mesh to exist. That is the whole
	 * point: it is what you run to bring the mesh into existence.
	 *
	 * ONE LINE ON STDOUT AND EVERY DIAGNOSTIC ON STDERR, because the consumer
	 * is a program. `meq-run` reads it with shlex.split, so the quoting here
	 * is the contract between the two.
	 */
	if ( wantMeshCommand )
	{
		try
		{
			meq::Configuration const query( argument );
			// NOTHING TO MAKE IS NOT AN ERROR, AND THE SILENCE IS THE
			// ANSWER. Most configurations name a mesh that already exists or
			// a box this driver builds itself, and a wrapper asking "what
			// makes this file's mesh?" has to be able to tell that apart from
			// a file it could not read -- so empty output with exit 0 means
			// "nothing", and a non-zero exit means the question was bad.
			if ( !query.getMesh().generate.given )
			{
				std::fprintf( stderr,
					"MEQ: %s has no [mesh.generate], so there is no mesh to\n"
					"     make -- [mesh] names one that already exists, or a\n"
					"     box this driver builds itself.\n", argument.c_str() );
				return Solved;
			}
			std::printf( "%s\n", meshCommand( query ).c_str() );
			return Solved;
		}
		catch ( meq::ConfigError const &error )
		{
			std::fprintf( stderr, "MEQ: %s\n", error.what() );
			return ConfigurationError;
		}
		catch ( std::exception const &error )
		{
			std::fprintf( stderr, "MEQ: %s\n", error.what() );
			return ConfigurationError;
		}
	}

	meq::RunOptions options;
	options.configPath = configPath;
	options.version = MEQ_VERSION;
	options.deviceName = deviceName;
	options.deviceGiven = deviceGiven;
	options.meshReady = meshReady;
	options.wantProfile = wantProfile;
	options.wantTraceSolverStats = wantTraceSolverStats;
	options.traceRefinementCap = traceRefinementCap;
	options.qIterationCap = qIterationCap;
	options.wantGridFunctions = wantGridFunctions;
	options.wantVtk = wantVtk;
	return meq::Run( options ).execute();
}
