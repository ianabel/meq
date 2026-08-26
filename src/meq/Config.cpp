#include "Config.hpp"

#include <algorithm>
#include <cstdint>
#include <initializer_list>
#include <limits>
#include <map>
#include <string>
#include <vector>

/*
 * Parsing and validation of a meq configuration file. See Config.hpp for the
 * shape of the schema; this file is the only place that knows the spelling of a
 * key or the value of a default.
 *
 * The key names below are UpperCamelCase and the C++ members they are read into
 * are lowerCamelCase. That is deliberate; Config.hpp says why.
 *
 * Two things here are worth the space they take.
 *
 * First, every number is read through asFloat()/asInteger() rather than through
 * toml::find<double>. TOML distinguishes 1 from 1.0, and toml11 does too:
 * find<double> on the node `RMin = 0` throws toml::type_error, and -- worse --
 * find_or<double>( v, "RMin", 1.0 ) silently returns the *default* 1.0 for that
 * same node, because a failed conversion is indistinguishable from a missing
 * key. An author writing `RMin = 0` means the number zero, so both spellings
 * are accepted below. (The converse is not true: a float where a count belongs,
 * `NR = 4.0`, is rejected rather than truncated.)
 *
 * Second, an unknown key is an error, and the message names the nearest
 * accepted key. The configuration this replaces ignored unknown keys, and its
 * one shipped example paid for it: the example set OutputMeshFile and
 * PsiResolution where the parser read FinalMeshFile and CellSize, so half of it
 * did nothing -- and the file could not be loaded at all, since the key the
 * parser did require was missing. A misspelt key that changes nothing is the
 * worst failure mode a configuration format has, because the run still produces
 * a plausible answer.
 */

namespace meq
{

	namespace
	{

		std::string formatMessage( std::string const & file, std::string const & key, std::string const & message )
		{
			std::string text = "meq: configuration error";
			if ( !file.empty() )
				text += " in '" + file + "'";
			if ( !key.empty() )
				text += ", key '" + key + "'";
			text += ": " + message;
			return text;
		}

		// Levenshtein distance, for the did-you-mean on an unknown key.
		std::size_t editDistance( std::string const & a, std::string const & b )
		{
			std::vector< std::size_t > previous( b.size() + 1 );
			std::vector< std::size_t > current( b.size() + 1 );

			for ( std::size_t j = 0; j <= b.size(); ++j )
				previous[ j ] = j;

			for ( std::size_t i = 1; i <= a.size(); ++i )
			{
				current[ 0 ] = i;
				for ( std::size_t j = 1; j <= b.size(); ++j )
				{
					std::size_t const substitution = previous[ j - 1 ] + ( a[ i - 1 ] == b[ j - 1 ] ? 0 : 1 );
					current[ j ] = std::min( { previous[ j ] + 1, current[ j - 1 ] + 1, substitution } );
				}
				previous = current;
			}

			return previous[ b.size() ];
		}

		// The closest candidate to `key`, or "" when nothing is close enough to
		// be worth printing. A suggestion that looks nothing like what was
		// written sends the reader off to check a key they never typed.
		std::string nearestKey( std::string const & key, std::initializer_list< char const * > candidates )
		{
			std::string best;
			std::size_t bestDistance = std::numeric_limits< std::size_t >::max();

			for ( char const * candidate : candidates )
			{
				std::size_t const distance = editDistance( key, candidate );
				if ( distance < bestDistance )
				{
					bestDistance = distance;
					best = candidate;
				}
			}

			std::size_t const tolerance = std::max< std::size_t >( 2, key.size()/3 );
			return ( bestDistance <= tolerance ) ? best : std::string();
		}

		std::string listOf( std::initializer_list< char const * > names )
		{
			std::string text;
			for ( char const * name : names )
				text += ( text.empty() ? "" : ", " ) + std::string( name );
			return text;
		}

		// Accessors for one TOML table, carrying enough context to name the
		// offending key and file in any error they throw.
		//
		// An optional table that is absent is represented by a null `values`,
		// and every getFooOr() then returns its default -- so "[solver] omitted
		// entirely" and "[solver] present but empty" behave identically.
		class Table
		{
			public:
				Table( toml::value const & document, std::string const & tableName, std::string const & configFile, bool required )
					: values( nullptr ), name( tableName ), file( configFile )
				{
					if ( !document.contains( tableName ) )
					{
						if ( required )
							throw ConfigError( file, name, "required table [" + name + "] is missing" );
						return;
					}

					toml::value const & table = document.at( tableName );
					if ( !table.is_table() )
						throw ConfigError( file, name, "[" + name + "] must be a table, but is a " + toml::to_string( table.type() ) );

					values = &table;
				};

				/// A table nested inside another, as [boundary.shape] is inside
				/// [boundary]. Needed because the constructor above resolves a
				/// FLAT key -- document.contains( "boundary.shape" ) looks for a
				/// key spelled with a dot, which TOML does not create; the dotted
				/// header nests instead. An absent parent gives an absent child,
				/// which is what lets [boundary.shape] be optional without the
				/// caller checking twice.
				Table( Table const & parent, std::string const & childName,
				       std::string const & configFile, bool required )
					: values( nullptr ), name( parent.name + "." + childName ),
					  file( configFile )
				{
					if ( !parent.isPresent() || !parent.values->contains( childName ) )
					{
						if ( required )
							throw ConfigError( file, name, "required table [" + name + "] is missing" );
						return;
					}

					toml::value const & table = parent.values->at( childName );
					if ( !table.is_table() )
						throw ConfigError( file, name, "[" + name + "] must be a table, but is a " + toml::to_string( table.type() ) );

					values = &table;
				};

				bool isPresent() const { return values != nullptr; };

				bool has( std::string const & key ) const { return find( key ) != nullptr; };

				double getFloat( std::string const & key ) const { return asFloat( key, require( key ) ); };
				int getInteger( std::string const & key ) const { return asInteger( key, require( key ) ); };
				std::string getString( std::string const & key ) const { return asString( key, require( key ) ); };

				double getFloatOr( std::string const & key, double fallback ) const
				{
					toml::value const * value = find( key );
					return ( value == nullptr ) ? fallback : asFloat( key, *value );
				};

				int getIntegerOr( std::string const & key, int fallback ) const
				{
					toml::value const * value = find( key );
					return ( value == nullptr ) ? fallback : asInteger( key, *value );
				};

				/// An array of numbers, empty if the key is absent. Every element
				/// goes through asFloat(), so [ 1, 0.5 ] is accepted and a
				/// mistyped [ 1, "0.5" ] is refused by element rather than
				/// silently becoming a default -- the same reasoning as asFloat()
				/// itself, recorded at the top of this file.
				std::vector< double > getFloatArrayOr( std::string const & key ) const
				{
					toml::value const * value = find( key );
					if ( value == nullptr )
						return {};
					if ( !value->is_array() )
						fail( key, "must be an array of numbers, but is a " + toml::to_string( value->type() ) );

					std::vector< double > numbers;
					auto const & elements = value->as_array();
					numbers.reserve( elements.size() );
					for ( std::size_t i = 0; i < elements.size(); ++i )
						numbers.push_back( asFloat( key + "[" + std::to_string( i ) + "]",
						                            elements[ i ] ) );
					return numbers;
				};

				std::string getStringOr( std::string const & key, std::string const & fallback ) const
				{
					toml::value const * value = find( key );
					return ( value == nullptr ) ? fallback : asString( key, *value );
				};

				// Every key present in the table must appear in `accepted`.
				//
				// TODO (schema evolution, see Config.hpp): when a key is
				// renamed, its old spelling belongs here as a deprecated alias
				// that warns and still resolves, rather than being deleted.
				void rejectUnknownKeys( std::initializer_list< char const * > accepted ) const
				{
					if ( values == nullptr )
						return;

					for ( auto const & entry : values->as_table() )
					{
						auto matches = [ &entry ]( char const * candidate ) { return entry.first == candidate; };
						if ( std::any_of( accepted.begin(), accepted.end(), matches ) )
							continue;

						std::string message = "is not a key of [" + name + "]";
						std::string const suggestion = nearestKey( entry.first, accepted );
						if ( !suggestion.empty() )
							message += "; did you mean '" + suggestion + "'?";
						message += " accepted keys are: " + listOf( accepted );

						throw ConfigError( file, qualify( entry.first ), message );
					}
				};

				[[noreturn]] void fail( std::string const & key, std::string const & message ) const
				{
					throw ConfigError( file, qualify( key ), message );
				};

			private:
				std::string qualify( std::string const & key ) const { return name + "." + key; };

				toml::value const * find( std::string const & key ) const
				{
					if ( values == nullptr || !values->contains( key ) )
						return nullptr;
					return &values->at( key );
				};

				toml::value const & require( std::string const & key ) const
				{
					toml::value const * value = find( key );
					if ( value == nullptr )
						fail( key, "is required but was not specified" );
					return *value;
				};

				double asFloat( std::string const & key, toml::value const & value ) const
				{
					// Both spellings of a number are a number.
					if ( value.is_floating() )
						return value.as_floating();
					if ( value.is_integer() )
						return static_cast< double >( value.as_integer() );

					fail( key, "must be a number, but is a " + toml::to_string( value.type() ) );
				};

				int asInteger( std::string const & key, toml::value const & value ) const
				{
					if ( !value.is_integer() )
						fail( key, "must be an integer, but is a " + toml::to_string( value.type() ) );

					std::int64_t const wide = value.as_integer();
					if ( wide < std::numeric_limits< int >::min() || wide > std::numeric_limits< int >::max() )
						fail( key, "is out of range for an integer" );

					return static_cast< int >( wide );
				};

				std::string asString( std::string const & key, toml::value const & value ) const
				{
					if ( !value.is_string() )
						fail( key, "must be a string, but is a " + toml::to_string( value.type() ) );

					return value.as_string();
				};

				toml::value const * values;
				std::string name;
				std::string file;
		};

		toml::value parseFile( std::string const & fileName )
		{
			try
			{
				return toml::parse( fileName );
			}
			catch ( toml::file_io_error const & error )
			{
				throw ConfigError( fileName, "", std::string( "cannot be read: " ) + error.what() );
			}
			catch ( toml::syntax_error const & error )
			{
				throw ConfigError( fileName, "", std::string( "TOML syntax error:\n" ) + error.what() );
			}
		}

		toml::value parseString( std::string const & text, std::string const & source )
		{
			try
			{
				return toml::parse_str( text );
			}
			catch ( toml::syntax_error const & error )
			{
				throw ConfigError( source, "", std::string( "TOML syntax error:\n" ) + error.what() );
			}
		}

		SourceType toSourceType( Table const & source, std::string const & spelling )
		{
			static std::map< std::string, SourceType > const types =
			{
				{ "soloviev",     SourceType::Soloviev },
				{ "mhd",          SourceType::MHD },
				{ "manufactured", SourceType::Manufactured }
			};

			auto found = types.find( spelling );
			if ( found == types.end() )
				source.fail( "Type", "'" + spelling + "' is not a known source; accepted values are: \"soloviev\", \"mhd\", \"manufactured\"" );

			return found->second;
		}

		ShapeType toShapeType( Table const & shape, std::string const & spelling )
		{
			if ( spelling == "none" )
				return ShapeType::None;
			if ( spelling == "miller" )
				return ShapeType::Miller;
			if ( spelling == "mxh" )
				return ShapeType::Mxh;

			shape.fail( "Type", "'" + spelling + "' is not a known boundary shape; accepted values are: \"none\", \"miller\", \"mxh\"" );
		}

		BoundaryDataType toBoundaryDataType( Table const & boundary, std::string const & spelling )
		{
			static std::map< std::string, BoundaryDataType > const types =
			{
				{ "zero",  BoundaryDataType::Zero },
				{ "exact", BoundaryDataType::Exact }
			};

			auto found = types.find( spelling );
			if ( found == types.end() )
				boundary.fail( "Type", "'" + spelling + "' is not a known boundary condition; accepted values are: \"zero\", \"exact\"" );

			return found->second;
		}

	}

	ConfigError::ConfigError( std::string const & file, std::string const & key, std::string const & message )
		: std::runtime_error( formatMessage( file, key, message ) ), fileName( file ), keyName( key )
	{
	}

	SolovievParameters const & SourceConfig::getSoloviev() const
	{
		if ( type != SourceType::Soloviev )
			throw ConfigError( "", "source.Type", "the configured source is not a Solov'ev source" );
		return std::get< SolovievParameters >( parameters );
	}

	MHDParameters const & SourceConfig::getMHD() const
	{
		if ( type != SourceType::MHD )
			throw ConfigError( "", "source.Type", "the configured source is not an MHD source" );
		return std::get< MHDParameters >( parameters );
	}

	ManufacturedParameters const & SourceConfig::getManufactured() const
	{
		if ( type != SourceType::Manufactured )
			throw ConfigError( "", "source.Type", "the configured source is not a manufactured source" );
		return std::get< ManufacturedParameters >( parameters );
	}

	std::string OutputConfig::getMeshFile() const
	{
		return directory.empty() ? prefix + ".mesh" : directory + "/" + prefix + ".mesh";
	}

	std::string OutputConfig::getPsiFile() const
	{
		return directory.empty() ? prefix + "_psi.gf" : directory + "/" + prefix + "_psi.gf";
	}

	std::string OutputConfig::getGradPsiFile() const
	{
		return directory.empty() ? prefix + "_grad_psi.gf" : directory + "/" + prefix + "_grad_psi.gf";
	}

	Configuration::Configuration( std::string const & fileName )
		: Configuration( parseFile( fileName ), fileName )
	{
	}

	Configuration::Configuration( toml::value const & document, std::string const & source )
		: sourceName( source )
	{
		parse( document );
	}

	Configuration Configuration::fromString( std::string const & text, std::string const & source )
	{
		return Configuration( parseString( text, source ), source );
	}

	void Configuration::parse( toml::value const & document )
	{
		if ( !document.is_table() )
			throw ConfigError( sourceName, "", "the configuration must be a table of tables" );

		// Catch a misspelt or misplaced table before anything reports a key
		// missing from a table that is not the one the author meant to write.
		{
			std::initializer_list< char const * > const tables = { "mesh", "discretisation", "source", "boundary", "solver", "output" };
			for ( auto const & entry : document.as_table() )
			{
				auto matches = [ &entry ]( char const * candidate ) { return entry.first == candidate; };
				if ( std::any_of( tables.begin(), tables.end(), matches ) )
					continue;

				std::string message = "is not part of the meq schema";
				std::string const suggestion = nearestKey( entry.first, tables );
				if ( !suggestion.empty() )
					message += "; did you mean [" + suggestion + "]?";
				message += " the configuration consists of the tables [" + listOf( tables ) + "]";

				throw ConfigError( sourceName, entry.first, message );
			}
		}

		// [mesh]
		{
			Table mesh( document, "mesh", sourceName, true );
			mesh.rejectUnknownKeys( { "RMin", "RMax", "ZMin", "ZMax", "NR", "NZ", "RefinementLevels", "File" } );

			meshOptions.file = mesh.getStringOr( "File", "" );
			meshOptions.refinementLevels = mesh.getIntegerOr( "RefinementLevels", 0 );

			if ( meshOptions.refinementLevels < 0 )
				mesh.fail( "RefinementLevels", "must not be negative" );

			if ( !meshOptions.fromFile() )
			{
				meshOptions.rMin = mesh.getFloat( "RMin" );
				meshOptions.rMax = mesh.getFloat( "RMax" );
				meshOptions.zMin = mesh.getFloat( "ZMin" );
				meshOptions.zMax = mesh.getFloat( "ZMax" );
				meshOptions.nR = mesh.getIntegerOr( "NR", 1 );
				meshOptions.nZ = mesh.getIntegerOr( "NZ", 1 );

				if ( meshOptions.rMin < 0.0 )
					mesh.fail( "RMin", "must not be negative: r is a cylindrical radius" );
				if ( meshOptions.rMax <= meshOptions.rMin )
					mesh.fail( "RMax", "must be greater than RMin (RMin = " + std::to_string( meshOptions.rMin ) + ")" );
				if ( meshOptions.zMax <= meshOptions.zMin )
					mesh.fail( "ZMax", "must be greater than ZMin (ZMin = " + std::to_string( meshOptions.zMin ) + ")" );
				if ( meshOptions.nR < 1 )
					mesh.fail( "NR", "must be at least 1" );
				if ( meshOptions.nZ < 1 )
					mesh.fail( "NZ", "must be at least 1" );
			}
		}

		// [discretisation]
		{
			Table discretisation( document, "discretisation", sourceName, true );
			discretisation.rejectUnknownKeys( { "PolynomialDegree", "Tau" } );

			discretisationOptions.polynomialDegree = discretisation.getInteger( "PolynomialDegree" );
			discretisationOptions.tau = discretisation.getFloatOr( "Tau", 1.0 );

			if ( discretisationOptions.polynomialDegree < 0 )
				discretisation.fail( "PolynomialDegree", "must not be negative" );
			if ( !( discretisationOptions.tau > 0.0 ) )
				discretisation.fail( "Tau", "must be positive; the HDG stabilisation is tau = O(1) and defaults to 1.0" );
		}

		// [source]
		{
			Table source( document, "source", sourceName, true );

			sourceOptions.type = toSourceType( source, source.getString( "Type" ) );

			switch ( sourceOptions.type )
			{
				case SourceType::Soloviev:
				{
					source.rejectUnknownKeys( { "Type", "A" } );
					SolovievParameters parameters;
					parameters.a = source.getFloat( "A" );
					sourceOptions.parameters = parameters;
					break;
				}
				case SourceType::MHD:
				{
					source.rejectUnknownKeys( { "Type", "PPrimeFile", "GGPrimeFile", "Mu0" } );
					MHDParameters parameters;
					parameters.pPrimeFile = source.getString( "PPrimeFile" );
					parameters.ggPrimeFile = source.getString( "GGPrimeFile" );
					parameters.mu0 = source.getFloatOr( "Mu0", parameters.mu0 );
					if ( parameters.pPrimeFile.empty() )
						source.fail( "PPrimeFile", "must not be empty" );
					if ( parameters.ggPrimeFile.empty() )
						source.fail( "GGPrimeFile", "must not be empty" );
					if ( !( parameters.mu0 > 0.0 ) )
						source.fail( "Mu0", "must be positive" );
					sourceOptions.parameters = parameters;
					break;
				}
				case SourceType::Manufactured:
				{
					source.rejectUnknownKeys( { "Type", "R0", "Kr", "Kz" } );
					ManufacturedParameters parameters;
					parameters.r0 = source.getFloat( "R0" );
					parameters.kr = source.getFloat( "Kr" );
					parameters.kz = source.getFloat( "Kz" );
					sourceOptions.parameters = parameters;
					break;
				}
			}
		}

		// [boundary]
		{
			Table boundary( document, "boundary", sourceName, false );
			boundary.rejectUnknownKeys( { "Type", "shape" } );

			boundaryOptions.type = toBoundaryDataType( boundary, boundary.getStringOr( "Type", "zero" ) );

			if ( boundaryOptions.type == BoundaryDataType::Exact && sourceOptions.type == SourceType::MHD )
				boundary.fail( "Type", "\"exact\" needs a source with a known exact solution; the \"mhd\" source has none" );

			// [boundary.shape]
			{
				Table shape( boundary, "shape", sourceName, false );
				shape.rejectUnknownKeys( { "Type", "R0", "Z0", "MinorRadius", "Elongation",
				                           "Triangularity", "Squareness",
				                           "CosCoefficients", "SinCoefficients" } );

				ShapeConfig & s = boundaryOptions.shape;
				s.type = toShapeType( shape, shape.getStringOr( "Type", "none" ) );

				if ( s.type != ShapeType::None )
				{
					s.majorRadius = shape.getFloat( "R0" );
					s.centreHeight = shape.getFloatOr( "Z0", 0.0 );
					s.minorRadius = shape.getFloat( "MinorRadius" );
					s.elongation = shape.getFloatOr( "Elongation", 1.0 );
				}

				// Refuse the keys of the other shape rather than ignore them. A
				// Triangularity sitting unread under Type = "mxh" is a config that
				// says one thing and does another, which is exactly the class of
				// silent wrong answer this file's asFloat() note is about.
				if ( s.type == ShapeType::Miller )
				{
					s.triangularity = shape.getFloatOr( "Triangularity", 0.0 );
					s.squareness = shape.getFloatOr( "Squareness", 0.0 );
					if ( shape.has( "CosCoefficients" ) || shape.has( "SinCoefficients" ) )
						shape.fail( "Type", "\"miller\" takes Triangularity and Squareness, not harmonic coefficients; use \"mxh\" for those" );
				}
				else if ( s.type == ShapeType::Mxh )
				{
					s.cosCoefficients = shape.getFloatArrayOr( "CosCoefficients" );
					s.sinCoefficients = shape.getFloatArrayOr( "SinCoefficients" );
					if ( shape.has( "Triangularity" ) || shape.has( "Squareness" ) )
						shape.fail( "Type", "\"mxh\" takes CosCoefficients and SinCoefficients, not Triangularity or Squareness; those are \"miller\" spellings" );
					if ( s.cosCoefficients.empty() && s.sinCoefficients.empty() )
						shape.fail( "Type", "\"mxh\" was given no harmonics at all, which describes an ellipse; say Type = \"miller\" with Triangularity = 0 if that is what is meant" );
				}
			}
		}

		// [solver]
		{
			Table solver( document, "solver", sourceName, false );
			solver.rejectUnknownKeys( { "NewtonMaxIterations", "NewtonRelativeTolerance", "NewtonAbsoluteTolerance",
			                            "LinearMaxIterations", "LinearTolerance" } );

			solverOptions.newtonMaxIterations = solver.getIntegerOr( "NewtonMaxIterations", solverOptions.newtonMaxIterations );
			solverOptions.newtonRelativeTolerance = solver.getFloatOr( "NewtonRelativeTolerance", solverOptions.newtonRelativeTolerance );
			solverOptions.newtonAbsoluteTolerance = solver.getFloatOr( "NewtonAbsoluteTolerance", solverOptions.newtonAbsoluteTolerance );
			solverOptions.linearMaxIterations = solver.getIntegerOr( "LinearMaxIterations", solverOptions.linearMaxIterations );
			solverOptions.linearTolerance = solver.getFloatOr( "LinearTolerance", solverOptions.linearTolerance );

			if ( solverOptions.newtonMaxIterations < 1 )
				solver.fail( "NewtonMaxIterations", "must be at least 1" );
			if ( !( solverOptions.newtonRelativeTolerance > 0.0 ) )
				solver.fail( "NewtonRelativeTolerance", "must be positive" );
			if ( !( solverOptions.newtonAbsoluteTolerance > 0.0 ) )
				solver.fail( "NewtonAbsoluteTolerance", "must be positive" );
			if ( solverOptions.linearMaxIterations < 1 )
				solver.fail( "LinearMaxIterations", "must be at least 1" );
			if ( !( solverOptions.linearTolerance > 0.0 ) )
				solver.fail( "LinearTolerance", "must be positive" );
		}

		// [output]
		{
			Table output( document, "output", sourceName, false );
			output.rejectUnknownKeys( { "Directory", "Prefix" } );

			outputOptions.directory = output.getStringOr( "Directory", outputOptions.directory );
			outputOptions.prefix = output.getStringOr( "Prefix", outputOptions.prefix );

			if ( outputOptions.prefix.empty() )
				output.fail( "Prefix", "must not be empty; it is the stem of every output file name" );
		}
	}

}
