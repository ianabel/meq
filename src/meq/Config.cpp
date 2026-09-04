#include "Config.hpp"

#include "RotatingSource.hpp"

#include <algorithm>
#include <cstdint>
#include <initializer_list>
#include <limits>
#include <map>
#include <string>
#include <vector>

/*
 * Parsing and validation of a MEQ configuration file. See Config.hpp for the
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
			std::string text = "MEQ: configuration error";
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

				/// The elements of an ARRAY OF TABLES, as [[source.species]] is
				/// inside [source]. Empty if the key is absent.
				///
				/// The nested-table constructor above cannot do this: an array of
				/// tables *is* an array, so its is_table() check refuses one --
				/// even though TOML nests [[a.b]] under [a] exactly as [a.b]
				/// does. This is the first array of tables in MEQ's schema.
				///
				/// Each element is named "parent.key[i]", so fail() and
				/// rejectUnknownKeys() on a returned Table qualify a diagnostic
				/// all the way down to source.species[2].Mass without the caller
				/// doing anything. That is the same convention getFloatArrayOr()
				/// uses for an element of a scalar array.
				std::vector< Table > getTableArrayOr( std::string const & key ) const
				{
					std::vector< Table > tables;

					toml::value const * value = find( key );
					if ( value == nullptr )
						return tables;
					if ( !value->is_array() )
						fail( key, "must be an array of tables, written [[" + qualify( key )
						           + "]], but is a " + toml::to_string( value->type() ) );

					auto const & elements = value->as_array();
					tables.reserve( elements.size() );
					for ( std::size_t i = 0; i < elements.size(); ++i )
					{
						std::string const elementName = qualify( key ) + "[" + std::to_string( i ) + "]";
						if ( !elements[ i ].is_table() )
							throw ConfigError( file, elementName, "must be a table, but is a "
							                   + toml::to_string( elements[ i ].type() ) );

						tables.push_back( Table( &elements[ i ], elementName, file ) );
					}

					return tables;
				};
				bool getBooleanOr( std::string const & key, bool fallback ) const
				{
					toml::value const * value = find( key );
					return ( value == nullptr ) ? fallback : asBoolean( key, *value );
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
				/// Wrap a node that has already been resolved, with a name given
				/// outright rather than built from a parent's. Used only by
				/// getTableArrayOr(), whose elements are named "parent.key[i]".
				Table( toml::value const * node, std::string const & elementName, std::string const & configFile )
					: values( node ), name( elementName ), file( configFile )
				{
				};

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

				// Explicit, for the reason recorded at the top of this file: a
				// find_or<bool> would turn Enabled = "true" into a silent false.
				bool asBoolean( std::string const & key, toml::value const & value ) const
				{
					if ( !value.is_boolean() )
						fail( key, "must be true or false, but is a " + toml::to_string( value.type() ) );

					return value.as_boolean();
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

		/// Read a profile given EITHER as a constant in `key` OR as a path in
		/// `key + "File"`. Exactly one, or -- when not required -- neither.
		///
		/// Two keys rather than one that dispatches on node type. TOML
		/// distinguishes 1 from 1.0 and toml11 does too, so a single key meaning
		/// "a number or a filename" would inherit exactly the trap recorded at
		/// the top of this file: the failure mode is a silent default rather
		/// than a refusal. Two keys make "which did you mean" a question the
		/// parser answers rather than one it guesses at.
		void readEitherProfile( Table const & table, std::string const & key,
		                        double & value, std::string & fileName, double & scale,
		                        bool & given, bool required )
		{
			std::string const fileKey = key + "File";
			std::string const scaleKey = key + "Scale";

			// RESERVED, AND REFUSED WITH THEIR OWN MESSAGE RATHER THAN AS
			// UNKNOWN KEYS. A NetCDF file holding several profiles at once --
			// MaNTA writes one, as groups Var<i> carrying u and q -- is wanted
			// and is not written. Naming the keys now costs nothing and buys
			// two things: a reader who tries it is told it is not implemented
			// rather than told the key does not exist, and adding it later is
			// purely additive rather than a schema change.
			//
			// When it is written, the rule is settled: a variable carrying
			// VALUES BUT NO DERIVATIVE is refused unless the configuration opts
			// in per profile, with e.g. TemperatureFit = "pchip". Profile's
			// contract is that prime() is the exact derivative of operator(),
			// and a fitted one is a modelling choice -- a monotone fit and a
			// natural spline disagree about whether an edge pedestal overshoots
			// into a negative density. That choice belongs to the person who
			// knows what the file means, in the file that records the run.
			for ( char const * reserved : { "Variable", "Fit" } )
			{
				if ( table.has( key + reserved ) )
					table.fail( key + reserved, std::string( "is reserved for reading profiles out of a NetCDF file holding several, which is not implemented yet; give " )
					            + key + " as a constant or " + fileKey + " as a path to a text table" );
			}

			bool const hasValue = table.has( key );
			bool const hasFile = table.has( fileKey );

			if ( hasValue && hasFile )
				table.fail( key, "and " + fileKey + " must not both be given: one is a constant profile and the other a tabulated one" );

			scale = table.getFloatOr( scaleKey, 1.0 );
			if ( !std::isfinite( scale ) )
				table.fail( scaleKey, "must be finite" );
			if ( !hasValue && !hasFile && table.has( scaleKey ) )
				table.fail( scaleKey, "means nothing without " + key + " or " + fileKey + " to scale" );

			if ( hasValue )
			{
				value = table.getFloat( key );
				given = true;
				return;
			}
			if ( hasFile )
			{
				fileName = table.getString( fileKey );
				if ( fileName.empty() )
					table.fail( fileKey, "must not be empty" );
				given = true;
				return;
			}

			given = false;
			if ( required )
				table.fail( key, "is required: give it as a constant, or give " + fileKey + " as a path to a table" );
		}

		/// [source] ProfileFile names one NetCDF holding several profiles, which
		/// the per-profile Variable keys would then select from. Reserved, for
		/// the reason recorded in readEitherProfile.
		void refuseReservedProfileFile( Table const & source )
		{
			if ( source.has( "ProfileFile" ) )
				source.fail( "ProfileFile", "is reserved for a NetCDF file holding several profiles at once, which is not implemented yet; give each profile as a constant or name a text table with its own <Profile>File key" );
		}

		/// [solver] LinearMaxIterations and LinearTolerance describe an
		/// ITERATIVE inner solve, and MEQ's trace solve is DIRECT -- UMFPACK,
		/// PARDISO or cuDSS -- so neither has anything to control.
		///
		/// THEY USED TO BE PARSED AND VALIDATED AND READ BY NOTHING, which is
		/// worse than either accepting or rejecting them: a key that validates
		/// is a key its author believes is doing something. The only Krylov
		/// solver in the tree is the fallback for a build with no direct solver
		/// at all, and its three sites in GradShafranov.cpp set their own
		/// numbers -- which do not even agree with the defaults these keys had,
		/// 5000 against 1000 and 1e-14 against 1e-12, so passing them through
		/// would have quietly changed that path rather than configured it.
		///
		/// Refused rather than ignored, on the same principle as
		/// refuseReservedProfileFile below and the PsiAxis case after it. They
		/// become meaningful the day the trace solve stops being direct, which
		/// 3D or a parallel build would force; reintroducing a key then is
		/// cheap, and a dead key in the meantime is a standing invitation to
		/// believe it works.
		void refuseIterativeSolverKeys( Table const & solver )
		{
			for ( char const * key : { "LinearMaxIterations", "LinearTolerance" } )
				if ( solver.has( key ) )
					solver.fail( key, "configures an iterative linear solve, and MEQ solves the hybridized trace system with a DIRECT solver -- UMFPACK, PARDISO or cuDSS -- which has no iteration count and no tolerance to set. Remove the key" );
		}

		/// The psi_ax guess, which is required exactly when the profiles are in
		/// normalised flux and meaningless otherwise. Refusing it in the
		/// un-normalised case is deliberate: a file carrying a PsiAxis that
		/// nothing reads is a file whose author believed something false about
		/// what the run was doing.
		void readNormalisation( Table const & source, bool normalised, double & psiAxis )
		{
			if ( !normalised )
			{
				if ( source.has( "PsiAxis" ) )
					source.fail( "PsiAxis", "means nothing unless Normalised = true, when psi_ax becomes an unknown of the system; without it the profiles are functions of psi itself" );
				return;
			}

			psiAxis = source.getFloat( "PsiAxis" );
			if ( !std::isfinite( psiAxis ) || psiAxis == 0.0 )
				source.fail( "PsiAxis", "must be finite and non-zero: Psi = psi/psi_ax is undefined at zero" );
		}

		SourceType toSourceType( Table const & source, std::string const & spelling )
		{
			static std::map< std::string, SourceType > const types =
			{
				{ "soloviev",     SourceType::Soloviev },
				{ "mhd",          SourceType::MHD },
				{ "manufactured", SourceType::Manufactured },
				{ "rotating",     SourceType::Rotating }
			};

			auto found = types.find( spelling );
			if ( found == types.end() )
				source.fail( "Type", "'" + spelling + "' is not a known source; accepted values are: \"soloviev\", \"mhd\", \"manufactured\", \"rotating\"" );

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

	RotatingParameters const & SourceConfig::getRotating() const
	{
		if ( type != SourceType::Rotating )
			throw ConfigError( "", "source.Type", "the configured source is not a rotating source" );
		return std::get< RotatingParameters >( parameters );
	}

	bool SourceConfig::isNormalised() const
	{
		switch ( type )
		{
			case SourceType::MHD:
				return std::get< MHDParameters >( parameters ).normalised;
			case SourceType::Rotating:
				return std::get< RotatingParameters >( parameters ).normalised;
			case SourceType::Soloviev:
			case SourceType::Manufactured:
				return false;
		}

		return false;
	}

	double SourceConfig::psiAxisGuess() const
	{
		switch ( type )
		{
			case SourceType::MHD:
				return std::get< MHDParameters >( parameters ).psiAxis;
			case SourceType::Rotating:
				return std::get< RotatingParameters >( parameters ).psiAxis;
			case SourceType::Soloviev:
			case SourceType::Manufactured:
				return 0.0;
		}

		return 0.0;
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

	std::string OutputConfig::getPsiStarFile() const
	{
		return directory.empty() ? prefix + "_psistar.gf" : directory + "/" + prefix + "_psistar.gf";
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
			std::initializer_list< char const * > const tables = { "mesh", "discretisation", "source", "boundary", "solver", "output", "initialguess", "adaptivity" };
			for ( auto const & entry : document.as_table() )
			{
				auto matches = [ &entry ]( char const * candidate ) { return entry.first == candidate; };
				if ( std::any_of( tables.begin(), tables.end(), matches ) )
					continue;

				std::string message = "is not part of the MEQ schema";
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
					source.rejectUnknownKeys( { "Type", "PPrimeFile", "GGPrimeFile", "Mu0",
					                            "PPrimeScale", "GGPrimeScale",
					                            "PPrimeVariable", "PPrimeFit",
					                            "GGPrimeVariable", "GGPrimeFit",
					                            "Normalised", "PsiAxis", "ProfileFile" } );
					refuseReservedProfileFile( source );
					MHDParameters parameters;
					parameters.pPrimeFile = source.getString( "PPrimeFile" );
					parameters.ggPrimeFile = source.getString( "GGPrimeFile" );
					parameters.mu0 = source.getFloatOr( "Mu0", parameters.mu0 );
					parameters.pPrimeScale = source.getFloatOr( "PPrimeScale", 1.0 );
					parameters.ggPrimeScale = source.getFloatOr( "GGPrimeScale", 1.0 );
					parameters.normalised = source.getBooleanOr( "Normalised", false );
					if ( !std::isfinite( parameters.pPrimeScale ) )
						source.fail( "PPrimeScale", "must be finite" );
					if ( !std::isfinite( parameters.ggPrimeScale ) )
						source.fail( "GGPrimeScale", "must be finite" );
					if ( parameters.pPrimeFile.empty() )
						source.fail( "PPrimeFile", "must not be empty" );
					if ( parameters.ggPrimeFile.empty() )
						source.fail( "GGPrimeFile", "must not be empty" );
					if ( !( parameters.mu0 > 0.0 ) )
						source.fail( "Mu0", "must be positive" );
					readNormalisation( source, parameters.normalised, parameters.psiAxis );
					sourceOptions.parameters = parameters;
					break;
				}
				case SourceType::Rotating:
				{
					source.rejectUnknownKeys( { "Type", "species", "Omega", "OmegaFile", "OmegaScale",
					                            "OmegaVariable", "OmegaFit",
					                            "GGPrime", "GGPrimeFile", "GGPrimeScale",
					                            "GGPrimeVariable", "GGPrimeFit",
					                            "ReferenceRadius", "Mu0", "Normalised", "PsiAxis",
					                            "ProfileFile" } );
					refuseReservedProfileFile( source );
					RotatingParameters parameters;

					parameters.referenceRadius = source.getFloatOr( "ReferenceRadius", parameters.referenceRadius );
					parameters.mu0 = source.getFloatOr( "Mu0", parameters.mu0 );
					parameters.normalised = source.getBooleanOr( "Normalised", false );

					if ( !( parameters.referenceRadius > 0.0 ) )
						source.fail( "ReferenceRadius", "must be positive: it is the radius at which phi_0 vanishes and at which each Density is the physical density" );
					if ( !( parameters.mu0 > 0.0 ) )
						source.fail( "Mu0", "must be positive" );

					// Rotation is OPTIONAL. Both keys absent means omega = 0, and
					// the source then reduces to the static equation, which is a
					// useful thing to be able to ask for from the same file.
					readEitherProfile( source, "Omega", parameters.omega, parameters.omegaFile,
					                   parameters.omegaScale, parameters.omegaGiven, false );
					bool ggGiven = false;
					readEitherProfile( source, "GGPrime", parameters.ggPrime, parameters.ggPrimeFile,
					                   parameters.ggPrimeScale, ggGiven, true );

					readNormalisation( source, parameters.normalised, parameters.psiAxis );

					std::vector< Table > const species = source.getTableArrayOr( "species" );
					if ( species.size() < 2 )
						source.fail( "species", "a rotating source needs at least two species, written as [[source.species]] blocks; quasineutrality is what determines phi_0 and it needs charges of both signs" );

					// The ceiling belongs here as well as in the source, and not
					// only for symmetry: meq::RotatingSource does refuse an
					// over-long set, but it does so after the factory has opened
					// and parsed every profile file the species named. Failing at
					// the parse means the diagnostic is about the configuration
					// rather than about the ninth file.
					if ( species.size() > maxSpecies )
						source.fail( "species", "a rotating source carries at most " + std::to_string( maxSpecies )
						             + " species and this one has " + std::to_string( species.size() )
						             + "; the cap is meq::maxSpecies, which exists so that the per-quadrature-point work allocates nothing" );

					int neutralising = 0;
					bool anyPositive = false;
					bool anyNegative = false;
					for ( std::size_t i = 0; i < species.size(); ++i )
					{
						Table const & one = species[ i ];
						one.rejectUnknownKeys( { "Name", "Mass", "Charge", "Temperature",
						                         "TemperatureFile", "TemperatureScale",
						                         "TemperatureVariable", "TemperatureFit",
						                         "Density", "DensityFile", "DensityScale",
						                         "DensityVariable", "DensityFit",
						                         "Neutralising" } );

						SpeciesParameters entry;
						entry.name = one.getStringOr( "Name", "species" + std::to_string( i ) );
						entry.mass = one.getFloat( "Mass" );
						entry.charge = one.getFloat( "Charge" );
						entry.neutralising = one.getBooleanOr( "Neutralising", false );

						if ( !( entry.mass > 0.0 ) )
							one.fail( "Mass", "must be positive" );
						if ( entry.charge == 0.0 )
							one.fail( "Charge", "must not be zero: it is Z_s, signed and dimensionless, so +1 for a proton and -1 for an electron" );

						bool temperatureGiven = false;
						readEitherProfile( one, "Temperature", entry.temperature,
						                   entry.temperatureFile, entry.temperatureScale,
						                   temperatureGiven, true );

						bool densityGiven = false;
						readEitherProfile( one, "Density", entry.density, entry.densityFile,
						                   entry.densityScale, densityGiven, false );

						if ( entry.neutralising )
						{
							++neutralising;
							if ( densityGiven )
								one.fail( "Density", "must not be given for the Neutralising species: its density is what charge neutrality determines, which is the whole point of marking it" );
						}
						else if ( !densityGiven )
						{
							one.fail( "Density", "is required unless this species is Neutralising; give it, or a DensityFile, or set Neutralising = true on exactly one species" );
						}

						anyPositive = anyPositive || entry.charge > 0.0;
						anyNegative = anyNegative || entry.charge < 0.0;

						parameters.species.push_back( entry );
					}

					if ( neutralising != 1 )
						source.fail( "species", "exactly one species must set Neutralising = true, and " + std::to_string( neutralising ) + " did. Fixing the gauge removes one function's worth of freedom from the densities, so for n species there are n - 1 independent ones" );
					if ( !anyPositive || !anyNegative )
						source.fail( "species", "the species must carry charges of both signs, or quasineutrality has no solution" );

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
			refuseIterativeSolverKeys( solver );

			if ( solverOptions.newtonMaxIterations < 1 )
				solver.fail( "NewtonMaxIterations", "must be at least 1" );
			if ( !( solverOptions.newtonRelativeTolerance > 0.0 ) )
				solver.fail( "NewtonRelativeTolerance", "must be positive" );
			if ( !( solverOptions.newtonAbsoluteTolerance > 0.0 ) )
				solver.fail( "NewtonAbsoluteTolerance", "must be positive" );
		}

		// [output]
		{
			Table output( document, "output", sourceName, false );
			output.rejectUnknownKeys( { "Directory", "Prefix", "GridNR", "GridNZ" } );

			outputOptions.directory = output.getStringOr( "Directory", outputOptions.directory );
			outputOptions.prefix = output.getStringOr( "Prefix", outputOptions.prefix );
			outputOptions.gridNR = output.getIntegerOr( "GridNR", outputOptions.gridNR );
			outputOptions.gridNZ = output.getIntegerOr( "GridNZ", outputOptions.gridNZ );

			if ( outputOptions.prefix.empty() )
				output.fail( "Prefix", "must not be empty; it is the stem of every output file name" );
			if ( outputOptions.gridNR < 2 )
				output.fail( "GridNR", "must be at least 2: these are grid NODES, so the "
				             "spacing is ( RMax - RMin )/( GridNR - 1 )" );
			if ( outputOptions.gridNZ < 2 )
				output.fail( "GridNZ", "must be at least 2: these are grid NODES, so the "
				             "spacing is ( ZMax - ZMin )/( GridNZ - 1 )" );
		}

		// [initialguess]
		{
			Table guess( document, "initialguess", sourceName, false );
			guess.rejectUnknownKeys( { "Type", "File", "MeshFile", "Amplitude" } );

			std::string const type = guess.getStringOr( "Type", "none" );
			if ( type == "none" )
				initialGuessOptions.type = InitialGuessType::None;
			else if ( type == "ramp" )
				initialGuessOptions.type = InitialGuessType::Ramp;
			else if ( type == "gridfunction" )
				initialGuessOptions.type = InitialGuessType::GridFunction;
			else
				guess.fail( "Type", "must be one of none, ramp, gridfunction, but is \"" + type + "\"" );

			initialGuessOptions.file = guess.getStringOr( "File", "" );
			initialGuessOptions.meshFile = guess.getStringOr( "MeshFile", "" );
			initialGuessOptions.amplitude =
				guess.getFloatOr( "Amplitude", initialGuessOptions.amplitude );

			if ( initialGuessOptions.type == InitialGuessType::GridFunction )
			{
				if ( initialGuessOptions.file.empty() )
					guess.fail( "File", "is required when Type = \"gridfunction\"" );
				if ( initialGuessOptions.meshFile.empty() )
					guess.fail( "MeshFile", "is required when Type = \"gridfunction\": a "
					            "GridFunction cannot be read without the mesh it lives on" );
			}

			if ( initialGuessOptions.type == InitialGuessType::Ramp
			     && !( initialGuessOptions.amplitude > 0.0 ) )
				guess.fail( "Amplitude", "must be positive: the point of the ramp is that psi "
				            "crosses zero in the INTERIOR, and an amplitude of zero puts the "
				            "iteration straight onto the trivial branch it exists to avoid" );
		}

		// [adaptivity]
		{
			Table adaptivity( document, "adaptivity", sourceName, false );
			adaptivity.rejectUnknownKeys( { "Enabled", "MaxIterations", "Strategy",
			                                "Theta", "TargetError" } );

			adaptivityOptions.enabled =
				adaptivity.getBooleanOr( "Enabled", adaptivityOptions.enabled );
			adaptivityOptions.maxIterations =
				adaptivity.getIntegerOr( "MaxIterations", adaptivityOptions.maxIterations );
			adaptivityOptions.theta =
				adaptivity.getFloatOr( "Theta", adaptivityOptions.theta );
			adaptivityOptions.targetError =
				adaptivity.getFloatOr( "TargetError", adaptivityOptions.targetError );

			std::string const strategy = adaptivity.getStringOr( "Strategy", "doerfler" );
			if ( strategy == "doerfler" )
				adaptivityOptions.strategy = MarkingStrategy::Doerfler;
			else if ( strategy == "maximum" )
				adaptivityOptions.strategy = MarkingStrategy::Maximum;
			else
				adaptivity.fail( "Strategy", "must be doerfler or maximum, but is \"" + strategy + "\"" );

			if ( adaptivityOptions.maxIterations < 1 )
				adaptivity.fail( "MaxIterations", "must be at least 1" );
			if ( !( adaptivityOptions.theta > 0.0 ) || adaptivityOptions.theta > 1.0 )
				adaptivity.fail( "Theta", "must be in ( 0, 1 ]: it is the fraction of the total "
				                 "estimated error the marked elements must carry" );
			if ( !( adaptivityOptions.targetError > 0.0 ) )
				adaptivity.fail( "TargetError", "must be positive" );
		}
	}

}
