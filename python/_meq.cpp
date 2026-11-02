/*
 * _meq -- the `meq` Python package's extension: meq::Run, from Python.
 *
 * A THIRD FRONT END ON THE ONE RUN, beside `meq` and `meq-convert`, and not a
 * second driver (PYTHON-PLAN.md). Equilibrium( path ) is `meq path`, a phase
 * at a time: construction is meq::Run::prepare(), solve() and write() its
 * namesakes, and a run made through it writes what `meq` writes, coefficient
 * for coefficient -- PythonDriverParity holds that. python/meq/__init__.py is
 * what a user imports; this module is private to it.
 *
 * THE GIL IS HELD FOR THE WHOLE RUN, and the module does not declare
 * py::mod_gil_not_used(). That is MaNTA's reasoning
 * (../../manta/src/Python.cpp, above PYBIND11_MODULE) and it applies here:
 * meq_core has statics nobody has audited for a free-threaded caller, and
 * OpenMP's threads are not Python's, so holding the GIL costs a sequential
 * caller nothing.
 *
 * THE REPORT GOES TO THE PROCESS'S STDOUT AND STDERR, file descriptors 1 and
 * 2, as `meq`'s does: printed in a terminal, not captured by a notebook.
 * Python's own buffered stdout is flushed before each phase so the two
 * interleave in the order they were written.
 */

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <toml.hpp>

#include "meq_version.hpp"

#include "meq/BuildInfo.hpp"
#include "meq/Run.hpp"

#include <cstdint>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <utility>

namespace py = pybind11;

namespace
{
	/// meq.Error and its three kinds, one per failing exit code. The run has
	/// already said why on stderr; the exception says which kind.
	struct Failure : std::runtime_error
	{
		using std::runtime_error::runtime_error;
	};
	struct ConfigFailure : Failure
	{
		using Failure::Failure;
	};
	struct SolveFailure : Failure
	{
		using Failure::Failure;
	};
	struct OutputFailure : Failure
	{
		using Failure::Failure;
	};

	/// A phase's exit code as Python sees it: nothing for success and for an
	/// inverse best fit, which is a status and not a failure, and an exception
	/// for the rest.
	void raiseFor( int code, char const *phase )
	{
		std::fflush( stdout );
		std::fflush( stderr );
		std::string const where = std::string( "meq: " ) + phase + "() ";
		switch ( code )
		{
			case meq::Solved:
			case meq::TargetsNotMet:
				return;
			case meq::ConfigurationError:
				throw ConfigFailure( where + "refused the configuration, or could "
				                     "not set the run up from it (exit code 1); "
				                     "the reason is on stderr" );
			case meq::SolveFailed:
				throw SolveFailure( where + "did not converge (exit code 2); "
				                    "the reason is on stderr" );
			case meq::OutputFailed:
				throw OutputFailure( where + "could not write the output "
				                     "(exit code 3); the reason is on stderr" );
			default:
				throw Failure( where + "ended with exit code "
				               + std::to_string( code ) );
		}
	}

	/// Before each phase: Python's stdout before the run's, so a caller's
	/// print() and the report come out in the order they were made; and
	/// meq/__init__.py's check for a library another module has loaded since
	/// meq was imported, which warns.
	void flushPython()
	{
		py::module_::import( "sys" ).attr( "stdout" ).attr( "flush" )();
		py::module_::import( "meq" ).attr( "_warn_if_foreign_libraries" )();
	}

	/*
	 * Python -> toml, recursively: ../../manta/src/PyToml.hpp's tomlFromPy(),
	 * so a dict means here what it means to a MaNTA Runner.
	 *
	 * BOOL BEFORE INT, because Python's bool is an int, and the other order
	 * turns True into 1. AN INTEGRAL FLOAT STAYS A FLOAT, because TOML tells 1
	 * from 1.0 and so does the configuration's reader. A NUMPY SCALAR is
	 * unwrapped through item() and tested again: np.float64 subclasses float
	 * and would pass anyway, which is the trap, since np.float32 and every
	 * np.int* subclass nothing.
	 */
	toml::value tomlFromPy( py::handle obj )
	{
		if ( py::isinstance<py::bool_>( obj ) )
			return toml::value( py::cast<bool>( obj ) );
		if ( py::isinstance<py::int_>( obj ) )
			return toml::value( py::cast<std::int64_t>( obj ) );
		if ( py::isinstance<py::float_>( obj ) )
			return toml::value( py::cast<double>( obj ) );
		if ( py::isinstance<py::str>( obj ) )
			return toml::value( py::cast<std::string>( obj ) );

		if ( py::isinstance<py::dict>( obj ) )
		{
			toml::table table;
			for ( auto const &item : py::reinterpret_borrow<py::dict>( obj ) )
				table.emplace( py::cast<std::string>( py::str( item.first ) ),
				               tomlFromPy( item.second ) );
			return toml::value( std::move( table ) );
		}

		if ( py::isinstance<py::list>( obj ) || py::isinstance<py::tuple>( obj ) )
		{
			toml::array array;
			for ( auto const &item : py::reinterpret_borrow<py::sequence>( obj ) )
				array.push_back( tomlFromPy( item ) );
			return toml::value( std::move( array ) );
		}

		if ( py::hasattr( obj, "item" ) )
			return tomlFromPy( obj.attr( "item" )() );

		throw std::invalid_argument(
			"meq: cannot express a Python "
			+ py::cast<std::string>( py::str( py::type::of( obj ).attr( "__name__" ) ) )
			+ " as a TOML value. A configuration holds booleans, integers, "
			  "floats, strings, and lists or dicts of those." );
	}

	/// The options a Python caller may set: the process-level ones meq's
	/// command line sets, less the device. mfem::Device is process-wide and
	/// configured once, so it is not an argument of one equilibrium among
	/// several in a process.
	meq::RunOptions optionsFrom( bool meshReady, bool gridFunctions, bool vtk,
	                             bool profile )
	{
		meq::RunOptions options;
		options.version = MEQ_VERSION;
		options.meshReady = meshReady;
		options.wantGridFunctions = gridFunctions;
		options.wantVtk = vtk;
		options.wantProfile = profile;
		return options;
	}

	/// One run, a phase at a time.
	class Equilibrium
	{
		public:
			explicit Equilibrium( meq::RunOptions const &options )
				: run( options )
			{
				flushPython();
				raiseFor( run.prepare(), "Equilibrium" );
			}

			void solve()
			{
				flushPython();
				raiseFor( run.solve(), "solve" );
			}

			/// True unless an inverse solve wrote its best fit with the
			/// targets unmet, which `meq` reports as exit code 4.
			bool write()
			{
				flushPython();
				int const code = run.write();
				raiseFor( code, "write" );
				return code == meq::Solved;
			}

		private:
			meq::Run run;
	};
}

PYBIND11_MODULE( _meq, m )
{
	m.doc() = "meq::Run from Python. Import meq, not this.";

	auto const &error = py::register_exception<Failure>( m, "Error" );
	py::register_exception<ConfigFailure>( m, "ConfigError", error.ptr() );
	py::register_exception<SolveFailure>( m, "SolveError", error.ptr() );
	py::register_exception<OutputFailure>( m, "OutputError", error.ptr() );

	m.attr( "__version__" ) = MEQ_VERSION;
	m.def( "build_info", []()
	{
		py::dict info;
		info[ "version" ] = MEQ_VERSION;
		info[ "precision" ] = meq::realPrecision();
		info[ "device_backends" ] = meq::deviceBackends();
		info[ "description" ] = meq::buildDescription();
		return info;
	}, "The build: version and commit, precision, and device backends, as "
	   "`meq --version` and the .nc report them." );

	py::class_<Equilibrium>( m, "Equilibrium" )
		.def( py::init( []( std::string const &path, bool meshReady,
		                    bool gridFunctions, bool vtk, bool profile )
		      {
		          meq::RunOptions options
		              = optionsFrom( meshReady, gridFunctions, vtk, profile );
		          options.configPath = path;
		          return std::make_unique<Equilibrium>( options );
		      } ),
		      py::arg( "path" ), py::kw_only(),
		      py::arg( "mesh_ready" ) = false, py::arg( "grid_functions" ) = false,
		      py::arg( "vtk" ) = false, py::arg( "profile" ) = false,
		      "Read a configuration file and set the run up, as `meq path` "
		      "does before it solves. Raises meq.ConfigError where meq exits 1." )
		.def_static( "from_dict",
		      []( py::dict const &config, std::string const &baseDir,
		          std::string const &name, bool meshReady, bool gridFunctions,
		          bool vtk, bool profile )
		      {
		          meq::RunOptions options
		              = optionsFrom( meshReady, gridFunctions, vtk, profile );
		          options.configPath = name;
		          options.configText = toml::format( tomlFromPy( config ) );
		          options.baseDirectory = baseDir;
		          return std::make_unique<Equilibrium>( options );
		      },
		      py::arg( "config" ), py::kw_only(), py::arg( "base_dir" ) = "",
		      py::arg( "name" ) = "<dict>",
		      py::arg( "mesh_ready" ) = false, py::arg( "grid_functions" ) = false,
		      py::arg( "vtk" ) = false, py::arg( "profile" ) = false,
		      "The same from a dict of the TOML's tables. Relative input paths "
		      "resolve against base_dir, as a file's do against its own "
		      "directory; empty leaves them relative to the working directory. "
		      "`name` stands for the file in messages and in the output." )
		.def( "solve", &Equilibrium::solve,
		      "Solve, and refine if the configuration asks. Raises "
		      "meq.SolveError where meq exits 2. Once per Equilibrium." )
		.def( "write", &Equilibrium::write,
		      "Write what [output] and the options ask for. False where an "
		      "inverse solve met its targets only at a best fit (meq's exit "
		      "code 4); raises meq.OutputError where meq exits 3." );
}
