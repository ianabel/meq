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
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <toml.hpp>

#include "meq_version.hpp"

#include "meq/BuildInfo.hpp"
#include "meq/Run.hpp"

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <memory>
#include <stdexcept>
#include <system_error>
#include <type_traits>
#include <variant>

#include <unistd.h>
#include <string>
#include <utility>
#include <vector>

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

	/// Python's stdout before the run's, so a caller's print() and the
	/// report come out in the order they were made.
	void flushPython()
	{
		py::module_::import( "sys" ).attr( "stdout" ).attr( "flush" )();
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
	                             bool profile, meq::ProfileTables const &profiles )
	{
		meq::RunOptions options;
		options.version = MEQ_VERSION;
		options.profiles = profiles;
		options.meshReady = meshReady;
		options.wantGridFunctions = gridFunctions;
		options.wantVtk = vtk;
		options.wantProfile = profile;
		return options;
	}

	/**
	 * ONE EQUILIBRIUM, SOLVED AS OFTEN AS ITS PROFILES CHANGE.
	 *
	 * The first solve() is the run construction prepared. EACH LATER ONE IS A
	 * NEW meq::Run, RESTARTED EXACTLY FROM THE LAST: the current run writes
	 * its state to a private .nc (Run::writeRestart()), and the new one is
	 * prepared from the same configuration with that file as its
	 * [initialguess] and the profiles as they now stand -- the whole iterate
	 * ( q, psi, psihat ), labels and exterior coefficients, on the same mesh,
	 * through the restart path `meq` itself takes. The new run replaces the
	 * current one only if it solves, so a solve that fails leaves this object
	 * as it was: its state, its profiles' effect, what write() writes.
	 */
	class Equilibrium
	{
		public:
			explicit Equilibrium( meq::RunOptions const &optionsIn )
				: options( optionsIn ), run( std::make_unique<meq::Run>( optionsIn ) )
			{
				flushPython();
				raiseFor( run->prepare(), "Equilibrium" );
			}

			~Equilibrium()
			{
				if ( !scratch.empty() )
				{
					std::error_code ignored;
					std::filesystem::remove_all( scratch, ignored );
				}
			}

			Equilibrium( Equilibrium const & ) = delete;
			Equilibrium &operator=( Equilibrium const & ) = delete;

			/// Replace one profile from the next solve() on. Checked there,
			/// where the source is built, as at construction.
			void setProfile( std::string const &slot, meq::ProfileTable const &table )
			{
				options.profiles[ slot ] = table;
				profilesChanged = true;
			}

			void solve()
			{
				flushPython();
				if ( !solved )
				{
					// Profiles set before the first solve: the prepared run
					// built its source without them, so it is prepared again,
					// cold, as construction would have been.
					if ( profilesChanged )
					{
						auto fresh = std::make_unique<meq::Run>( options );
						raiseFor( fresh->prepare(), "solve" );
						run = std::move( fresh );
					}
					raiseFor( run->solve(), "solve" );
					solved = true;
					profilesChanged = false;
					return;
				}

				if ( !stateWritten )
					writeState( "solve" );
				meq::RunOptions again = options;
				again.restartFrom = scratchDirectory() + "/state.nc";
				// AN R_midplane TABLE IS HELD AT THE LAST ANSWER'S CROSSING, so
				// solving again with the tables unchanged moves the map to the
				// new answer's: the held map's outer loop. A support that is not
				// a box has no crossing to hold, and the source factory refuses
				// the table by its key, source.SupportShape.
				for ( auto const &entry : options.profiles )
					if ( entry.second.abscissa == meq::ProfileAbscissa::RMidplane )
					{
						if ( run->supportIsBox() )
							again.heldCrossing = crossingOrRaise();
						break;
					}
				auto next = std::make_unique<meq::Run>( again );
				raiseFor( next->prepare(), "solve" );
				raiseFor( next->solve(), "solve" );
				run = std::move( next );
				profilesChanged = false;
				stateWritten = false;
				written = false;
			}

			/// True unless an inverse solve wrote its best fit with the
			/// targets unmet, which `meq` reports as exit code 4.
			bool write()
			{
				flushPython();
				int const code = run->write();
				raiseFor( code, "write" );
				written = true;
				return code == meq::Solved;
			}

			/**
			 * THE .nc OF THE SOLVED STATE, IN MEMORY (PYTHON-PLAN.md, PY-4):
			 * what write() wrote, or, before any write(), the private state
			 * file a re-solve restarts from, written now and kept for it. The
			 * same write phase either way, so the same numbers.
			 */
			std::shared_ptr<meq::Dataset const> dataset()
			{
				if ( !solved )
					throw Failure( "meq: there is no solved state to read: solve() first" );
				if ( !written && !stateWritten )
					writeState( "results" );
				return run->dataset();
			}

			/// The solved state's _surfaces.nc, written privately if write()
			/// wrote none.
			std::shared_ptr<meq::Dataset const> surfaces()
			{
				if ( !solved )
					throw Failure( "meq: there is no solved state to read: solve() first" );
				return translated( [ & ]() { return run->fluxSurfaces( surfacesFile() ); } );
			}

			std::vector<std::pair<std::string, std::vector<double>>>
			surfacesAt( std::vector<double> const &labels )
			{
				if ( !solved )
					throw Failure( "meq: there is no solved state to read: solve() first" );
				return translated( [ & ]() { return run->fluxSurfacesAt( labels, surfacesFile() ); } );
			}

			/// The solved state's fields at points, as the .nc's grid has them.
			meq::Run::Samples at( std::vector<double> const &radius,
			                      std::vector<double> const &z )
			{
				if ( !solved )
					throw Failure( "meq: there is no solved state to read: solve() first" );
				if ( radius.size() != z.size() )
					throw Failure( "meq: at() takes as many R as Z, and has "
					               + std::to_string( radius.size() ) + " and "
					               + std::to_string( z.size() ) );
				if ( radius.empty() )
					return meq::Run::Samples();
				return run->sampleAt( radius, z );
			}

			/// The solved state's midplane crossing, R_0( Psi ).
			std::vector<meq::Knot> crossing()
			{
				return crossingOrRaise();
			}

		private:
			std::vector<meq::Knot> crossingOrRaise()
			{
				try
				{
					return run->midplaneCrossing();
				}
				catch ( std::logic_error const &error )
				{
					throw ConfigFailure( std::string( "meq: " ) + error.what() );
				}
				catch ( std::runtime_error const &error )
				{
					throw SolveFailure( std::string( "meq: the answer has no single midplane "
					                                 "crossing: " ) + error.what() );
				}
			}

			std::string surfacesFile()
			{
				return scratchDirectory() + "/state_surfaces.nc";
			}

			/// A library refusal as the exception it is in Python: a question
			/// this state cannot answer is meq.ConfigError, a label outside
			/// the cut or an extraction that fails meq.Error.
			template <typename Call>
			static auto translated( Call const &call ) -> decltype( call() )
			{
				try
				{
					return call();
				}
				catch ( std::logic_error const &error )
				{
					if ( dynamic_cast<std::domain_error const *>( &error ) )
						throw Failure( std::string( "meq: " ) + error.what() );
					throw ConfigFailure( std::string( "meq: " ) + error.what() );
				}
				catch ( std::runtime_error const &error )
				{
					throw Failure( std::string( "meq: no flux surfaces: " ) + error.what() );
				}
			}

			/// The state, to the private file a re-solve restarts from.
			void writeState( char const *phase )
			{
				raiseFor( run->writeRestart( scratchDirectory() + "/state" ), phase );
				stateWritten = true;
			}

			meq::RunOptions options;
			std::unique_ptr<meq::Run> run;
			bool solved = false;
			bool profilesChanged = false;
			/// Whether this run's state file, and whether its write(), are written.
			bool stateWritten = false;
			bool written = false;
			std::string scratch;

			/// A directory of this object's own for its restart file, made on
			/// the first re-solve and removed with the object.
			std::string const &scratchDirectory()
			{
				if ( scratch.empty() )
				{
					static int made = 0;
					std::filesystem::path const path = std::filesystem::temp_directory_path()
						/ ( "meq-" + std::to_string( ::getpid() ) + "-" + std::to_string( made++ ) );
					std::filesystem::create_directories( path );
					scratch = path.string();
				}
				return scratch;
			}
	};
}

namespace
{
	/// Numbers as a numpy array of @a shape; text as a str, or a list of them.
	template <typename T>
	py::object numbers( std::vector<T> const &values, std::vector<py::ssize_t> const &shape )
	{
		py::array_t<T> out( shape );
		std::copy( values.begin(), values.end(), out.mutable_data() );
		return std::move( out );
	}

	py::object valuesToPython( meq::DatasetValues const &values,
	                           std::vector<py::ssize_t> const &shape )
	{
		return std::visit( [ & ]( auto const &held ) -> py::object
		{
			using Element = typename std::decay_t<decltype( held )>::value_type;
			if constexpr ( std::is_same_v<Element, std::string> )
			{
				if ( shape.empty() )
					return py::str( held.empty() ? std::string() : held.front() );
				py::list out;
				for ( std::string const &one : held )
					out.append( py::str( one ) );
				return std::move( out );
			}
			else if ( shape.empty() )
				return py::cast( held.empty() ? Element() : held.front() );
			else
				return numbers( held, shape );
		}, values );
	}

	/// A scalar attribute as a Python scalar, an array of several as numpy.
	py::object attributeToPython( meq::DatasetAttribute const &attribute )
	{
		return std::visit( [ & ]( auto const &held ) -> py::object
		{
			using Held = std::decay_t<decltype( held )>;
			if constexpr ( std::is_same_v<Held, std::string> )
				return py::str( held );
			else if ( held.size() == 1 )
				return py::cast( held.front() );
			else
				return numbers( held, { static_cast<py::ssize_t>( held.size() ) } );
		}, attribute.value );
	}

	py::dict attributesToPython( std::vector<meq::DatasetAttribute> const &attributes )
	{
		py::dict out;
		for ( meq::DatasetAttribute const &one : attributes )
			out[ py::str( one.name ) ] = attributeToPython( one );
		return out;
	}

	std::vector<py::ssize_t> shapeOf( meq::Dataset const &data,
	                                  meq::DatasetVariable const &variable )
	{
		std::vector<py::ssize_t> shape;
		for ( std::string const &one : variable.dimensions )
			shape.push_back( static_cast<py::ssize_t>( data.dimension( one ) ) );
		return shape;
	}

	/// The whole root group: {"dimensions", "variables", "attributes"}, each
	/// variable {"dimensions", "values", "attributes"}.
	py::dict datasetToPython( meq::Dataset const &data )
	{
		py::dict dimensions;
		for ( auto const &one : data.dimensions() )
			dimensions[ py::str( one.first ) ] = one.second;
		py::dict variables;
		for ( meq::DatasetVariable const &one : data.variables() )
		{
			py::dict variable;
			variable[ "dimensions" ] = py::tuple( py::cast( one.dimensions ) );
			variable[ "values" ] = valuesToPython( one.values, shapeOf( data, one ) );
			variable[ "attributes" ] = attributesToPython( one.attributes );
			variables[ py::str( one.name ) ] = variable;
		}
		py::dict out;
		out[ "dimensions" ] = dimensions;
		out[ "variables" ] = variables;
		out[ "attributes" ] = attributesToPython( data.attributes() );
		return out;
	}

	/// The variables on the ( Z, R ) lattice and its two coordinates.
	py::dict gridToPython( meq::Dataset const &data )
	{
		py::dict out;
		for ( meq::DatasetVariable const &one : data.variables() )
		{
			std::vector<std::string> const &on = one.dimensions;
			bool const lattice = on == std::vector<std::string>{ "Z", "R" }
			                     || on == std::vector<std::string>{ "R" }
			                     || on == std::vector<std::string>{ "Z" };
			if ( lattice )
				out[ py::str( one.name ) ] = valuesToPython( one.values, shapeOf( data, one ) );
		}
		return out;
	}
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

	/*
	 * A PROFILE AS ARRAYS: knots in `abscissa`, the derivative in that variable
	 * too. THE ABSCISSA IS A REQUIRED KEYWORD AND HAS NO DEFAULT, because
	 * reading psi_n as Psi converges, at full order, to the wrong equilibrium
	 * (meq::ProfileAbscissa). Checked where the source is built, slot by slot.
	 */
	py::class_<meq::ProfileTable>( m, "Table",
		"A profile as arrays: x, the value and its derivative in x, in the "
		"abscissa named -- \"Psi\" (MEQ's normalised flux, 1 at the axis or a "
		"mirror's inner label, what a profile file is in), \"psi_n\" (1 - Psi), "
		"\"sqrt_psi_n\", \"psi\" (the poloidal flux itself, Wb/rad, which a "
		"normalised source takes only with every profile in psi), or "
		"\"R_midplane\" (m, a mirror's, held at the last answer's crossing). "
		"Units are the configuration key's." )
		.def( py::init( []( std::vector<double> const &x,
		                    std::vector<double> const &value,
		                    std::vector<double> const &derivative,
		                    std::string const &abscissa )
		      {
		          if ( value.size() != x.size() || derivative.size() != x.size() )
		              throw py::value_error(
		                  "meq.Table: x, value and derivative have lengths "
		                  + std::to_string( x.size() ) + ", "
		                  + std::to_string( value.size() ) + " and "
		                  + std::to_string( derivative.size() ) );
		          meq::ProfileTable table;
		          try
		          {
		              table.abscissa = meq::profileAbscissaNamed( abscissa );
		          }
		          catch ( std::invalid_argument const &error )
		          {
		              throw py::value_error( error.what() );
		          }
		          for ( std::size_t i = 0; i < x.size(); ++i )
		              table.knots.push_back( meq::Knot{ x[ i ], value[ i ], derivative[ i ] } );
		          return table;
		      } ),
		      py::arg( "x" ), py::arg( "value" ), py::arg( "derivative" ),
		      py::kw_only(), py::arg( "abscissa" ) )
		.def_property_readonly( "abscissa", []( meq::ProfileTable const &t )
		      { return std::string( meq::profileAbscissaName( t.abscissa ) ); } )
		.def_property_readonly( "x", []( meq::ProfileTable const &t )
		      {
		          std::vector<double> out;
		          for ( meq::Knot const &k : t.knots )
		              out.push_back( k.psi );
		          return out;
		      } )
		.def_property_readonly( "value", []( meq::ProfileTable const &t )
		      {
		          std::vector<double> out;
		          for ( meq::Knot const &k : t.knots )
		              out.push_back( k.value );
		          return out;
		      } )
		.def_property_readonly( "derivative", []( meq::ProfileTable const &t )
		      {
		          std::vector<double> out;
		          for ( meq::Knot const &k : t.knots )
		              out.push_back( k.derivative );
		          return out;
		      } )
		.def( "__repr__", []( meq::ProfileTable const &t )
		      {
		          return "meq.Table(" + std::to_string( t.knots.size() ) + " knots, abscissa=\""
		                 + meq::profileAbscissaName( t.abscissa ) + "\")";
		      } );

	py::class_<Equilibrium>( m, "Equilibrium" )
		.def( py::init( []( std::string const &path, bool meshReady,
		                    bool gridFunctions, bool vtk, bool profile,
		                    meq::ProfileTables const &profiles )
		      {
		          meq::RunOptions options
		              = optionsFrom( meshReady, gridFunctions, vtk, profile, profiles );
		          options.configPath = path;
		          return std::make_unique<Equilibrium>( options );
		      } ),
		      py::arg( "path" ), py::kw_only(),
		      py::arg( "mesh_ready" ) = false, py::arg( "grid_functions" ) = false,
		      py::arg( "vtk" ) = false, py::arg( "profile" ) = false,
		      py::arg( "profiles" ) = meq::ProfileTables(),
		      "Read a configuration file and set the run up, as `meq path` "
		      "does before it solves. `profiles` maps a slot -- \"Omega\", "
		      "\"GGPrime\", \"species.<Name>.Density\" or "
		      "\"species.<Name>.Temperature\", a rotating source's -- to a "
		      "meq.Table that replaces the profile the configuration gives it. "
		      "Raises meq.ConfigError where meq exits 1." )
		.def_static( "from_dict",
		      []( py::dict const &config, std::string const &baseDir,
		          std::string const &name, bool meshReady, bool gridFunctions,
		          bool vtk, bool profile, meq::ProfileTables const &profiles )
		      {
		          meq::RunOptions options
		              = optionsFrom( meshReady, gridFunctions, vtk, profile, profiles );
		          options.configPath = name;
		          options.configText = toml::format( tomlFromPy( config ) );
		          options.baseDirectory = baseDir;
		          return std::make_unique<Equilibrium>( options );
		      },
		      py::arg( "config" ), py::kw_only(), py::arg( "base_dir" ) = "",
		      py::arg( "name" ) = "<dict>",
		      py::arg( "mesh_ready" ) = false, py::arg( "grid_functions" ) = false,
		      py::arg( "vtk" ) = false, py::arg( "profile" ) = false,
		      py::arg( "profiles" ) = meq::ProfileTables(),
		      "The same from a dict of the TOML's tables. Relative input paths "
		      "resolve against base_dir, as a file's do against its own "
		      "directory; empty leaves them relative to the working directory. "
		      "`name` stands for the file in messages and in the output." )
		.def( "midplane_crossing", []( Equilibrium &self )
		      {
		          std::vector<double> psi, radius, slope;
		          for ( meq::Knot const &k : self.crossing() )
		          {
		              psi.push_back( k.psi );
		              radius.push_back( k.value );
		              slope.push_back( k.derivative );
		          }
		          py::dict out;
		          out[ "Psi" ] = psi;
		          out[ "R" ] = radius;
		          out[ "dR_dPsi" ] = slope;
		          return out;
		      },
		      "The solved state's midplane crossing: R where each surface crosses "
		      "z = 0, against Psi, with its slope, as a dict of lists. A mirror's; "
		      "an R_midplane table of the next solve is held at it." )
		.def( "attributes", []( Equilibrium &self )
		      { return attributesToPython( self.dataset()->attributes() ); },
		      "The .nc's global attributes for the solved state, under its names: "
		      "psi_axis, psi_boundary, newton_iterations, the gauge and the rest. "
		      "Raises meq.Error before a solve." )
		.def( "grid", []( Equilibrium &self ) { return gridToPython( *self.dataset() ); },
		      "The .nc's ( R, Z ) lattice: R and Z, inside and extrapolated, and "
		      "every field on it -- psi, B_R, B_Z and the rest -- as numpy arrays "
		      "shaped ( Z, R ), under the .nc's names and with its NaN outside." )
		.def( "dataset", []( Equilibrium &self ) { return datasetToPython( *self.dataset() ); },
		      "The whole .nc root group in memory: {'dimensions', 'variables', "
		      "'attributes'}, each variable {'dimensions', 'values', "
		      "'attributes'}. The numbers the file holds, from the calls that "
		      "write it." )
		.def( "surfaces", []( Equilibrium &self, py::object rho ) -> py::object
		      {
		          if ( rho.is_none() )
		              return datasetToPython( *self.surfaces() );
		          py::array_t<double, py::array::c_style | py::array::forcecast> labels( rho );
		          std::vector<py::ssize_t> const shape( labels.shape(), labels.shape() + labels.ndim() );
		          std::vector<double> const at( labels.data(), labels.data() + labels.size() );
		          py::dict out;
		          for ( auto const &quantity : self.surfacesAt( at ) )
		              out[ py::str( quantity.first ) ] = numbers( quantity.second, shape );
		          return std::move( out );
		      },
		      py::arg( "rho" ) = py::none(),
		      "The flux surfaces of the solved state. With no argument, the "
		      "_surfaces.nc root group in memory, as dataset() is the .nc's: the "
		      "family at [output]'s FluxSurfaceCount, angles and cuts. With rho, "
		      "an array of sqrt(Psi_N), the family's 1-D quantities there under "
		      "the file's names, by its monotone cubic; a label outside the cut "
		      "is refused, not extrapolated. A box support has no closed "
		      "surfaces and raises meq.ConfigError." )
		.def( "at", []( Equilibrium &self, py::array_t<double, py::array::c_style | py::array::forcecast> radius,
		                py::array_t<double, py::array::c_style | py::array::forcecast> z )
		      {
		          if ( radius.ndim() != z.ndim() || !std::equal( radius.shape(), radius.shape() + radius.ndim(), z.shape() ) )
		              throw Failure( "meq: at() takes R and Z of one shape" );
		          std::vector<py::ssize_t> const shape( radius.shape(), radius.shape() + radius.ndim() );
		          std::vector<double> const r( radius.data(), radius.data() + radius.size() );
		          std::vector<double> const height( z.data(), z.data() + z.size() );
		          meq::Run::Samples const found = self.at( r, height );
		          py::dict out;
		          out[ "psi" ] = numbers( found.psi, shape );
		          out[ "B_R" ] = numbers( found.bR, shape );
		          out[ "B_Z" ] = numbers( found.bZ, shape );
		          out[ "inside" ] = numbers( found.inside, shape );
		          out[ "extrapolated" ] = numbers( found.extrapolated, shape );
		          return out;
		      },
		      py::arg( "R" ), py::arg( "Z" ),
		      "psi, B_R and B_Z at points (R, Z), arrays of one shape, with the "
		      "inside and extrapolated masks, under the .nc's names: the write "
		      "phase's own sampling, so at a node of grid() the values are its. "
		      "NaN where inside is 0, outside the domain and its band." )
		.def( "set_profile", &Equilibrium::setProfile, py::arg( "slot" ), py::arg( "table" ),
		      "Replace one profile, as the constructors' `profiles` does, from "
		      "the next solve() on." )
		.def( "solve", &Equilibrium::solve,
		      "Solve, and refine if the configuration asks. Raises "
		      "meq.SolveError where meq exits 2. AGAIN after set_profile(): a "
		      "warm re-solve, restarted exactly from the last state on the same "
		      "mesh; one that fails leaves the Equilibrium as it was. A "
		      "configuration with [adaptivity] solves once." )
		.def( "write", &Equilibrium::write,
		      "Write what [output] and the options ask for. False where an "
		      "inverse solve met its targets only at a best fit (meq's exit "
		      "code 4); raises meq.OutputError where meq exits 3." );
}
