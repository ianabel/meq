#include "Config.hpp"

// maximumFilamentsPerCoil, so that a stack too fine to build is refused at
// parse with the block's name in the message rather than out of
// meq::filamentStack three layers down. MFEM-free, like everything Config
// reaches: Coils.hpp pulls Source.hpp and Profiles.hpp and nothing else.
#include "Coils.hpp"
// defaultCircuitName, the one circuit name [[circuits]] may not take.
#include "Materials.hpp"
#include "RotatingSource.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <initializer_list>
#include <limits>
#include <map>
#include <sstream>
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

		// THE BANNER BELONGS TO WHOEVER IS PRINTING, NOT TO THE MESSAGE.
		// apps/meq.cpp prefixes every line it writes with "MEQ: ", and its
		// catch blocks take std::exception, so a banner here as well produced
		// "MEQ: MEQ: configuration error in ...". That went unread for as long
		// as the project spelled itself in lower case, and the rename to MEQ
		// made it obvious. Every other exception MEQ throws identifies itself
		// by a namespace-qualified function name -- meq::SplineProfile::fromFile
		// and its kind -- rather than by a banner, so this one now says what
		// went wrong and where, and nothing about who it came from.
		std::string formatMessage( std::string const & file, std::string const & key, std::string const & message )
		{
			std::string text = "configuration error";
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

				/// True for an absent table AND for a present but empty one.
				/// Both mean "the author wrote nothing here", which is what a
				/// caller asking the question wants to know.
				bool isEmpty() const
				{
					return values == nullptr || values->as_table().empty();
				};

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

				/// As getFloatArrayOr(), for an array of INTEGERS -- a list of
				/// mesh element attributes, which is the only thing MEQ reads
				/// this way. An entry that is not an integer fails rather than
				/// truncating: `30.5` is not an attribute and silently becoming
				/// 30 is the class of error asFloat() exists to stop.
				std::vector< int > getIntegerArrayOr( std::string const & key ) const
				{
					toml::value const * value = find( key );
					if ( value == nullptr )
						return {};
					if ( !value->is_array() )
						fail( key, "must be an array of integers, but is a " + toml::to_string( value->type() ) );

					std::vector< int > numbers;
					auto const & elements = value->as_array();
					numbers.reserve( elements.size() );
					for ( std::size_t i = 0; i < elements.size(); ++i )
						numbers.push_back( asInteger( key + "[" + std::to_string( i ) + "]",
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
				/// The document itself, so that a TOP-LEVEL array of tables
				/// can be read the same way a nested one is. [[coils]] is the
				/// first of those; [[source.species]] is nested inside [source]
				/// and needed no such thing.
				static Table root( toml::value const & document,
				                   std::string const & configFile )
				{
					return Table( &document, "", configFile );
				};

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

				/// An empty name is the ROOT document, whose keys qualify to
				/// themselves -- [[coils]] rather than [[.coils]].
				std::string qualify( std::string const & key ) const
				{ return name.empty() ? key : name + "." + key; };

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

		/// Refuse `<Profile>Variable` and `<Profile>Fit`, the two keys that
		/// would select a profile out of a NetCDF file holding several.
		///
		/// RESERVED, AND REFUSED WITH THEIR OWN MESSAGE RATHER THAN AS UNKNOWN
		/// KEYS. A NetCDF file holding several profiles at once -- MaNTA writes
		/// one, as groups Var<i> carrying u and q -- is wanted and is not
		/// written. Naming the keys now costs nothing and buys two things: a
		/// reader who tries it is told it is not implemented rather than told
		/// the key does not exist, and adding it later is purely additive
		/// rather than a schema change. That is why they appear in the
		/// accepted-key lists: rejectUnknownKeys runs first, so a key has to be
		/// accepted before it can reach a message of its own.
		///
		/// When it is written, the rule is settled: a variable carrying VALUES
		/// BUT NO DERIVATIVE is refused unless the configuration opts in per
		/// profile, with e.g. TemperatureFit = "pchip". Profile's contract is
		/// that prime() is the exact derivative of operator(), and a fitted one
		/// is a modelling choice -- a monotone fit and a natural spline
		/// disagree about whether an edge pedestal overshoots into a negative
		/// density. That choice belongs to the person who knows what the file
		/// means, in the file that records the run.
		///
		/// A FREE FUNCTION RATHER THAN PART OF readEitherProfile, BECAUSE ONE
		/// PATH DOES NOT CALL THAT. The "mhd" source requires PPrimeFile and
		/// GGPrimeFile and offers no constant form, so it reads them directly
		/// and never reached the loop this replaces -- while still listing all
		/// four keys as accepted. The result was that GGPrimeVariable was a
		/// hard error under Type = "rotating" and silently inert under
		/// Type = "mhd": the same key name, two behaviours, decided by a
		/// neighbouring key. Accepted-and-ignored is the failure mode this file
		/// opens by warning about.
		///
		/// `constantAccepted` is what the advice at the end of the message
		/// turns on, since half the callers have no constant form to offer.
		void refuseReservedVariableKeys( Table const & table, std::string const & key,
		                                 bool constantAccepted )
		{
			std::string const fileKey = key + "File";

			for ( char const * reserved : { "Variable", "Fit" } )
			{
				if ( !table.has( key + reserved ) )
					continue;

				std::string const advice = constantAccepted
					? "give " + key + " as a constant or " + fileKey + " as a path to a text table"
					: "give " + fileKey + " as a path to a text table";

				table.fail( key + reserved, "is reserved for reading profiles out of a NetCDF file holding several, which is not implemented yet; " + advice );
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
		/**
		 * GGPrimeFile and SafetyFactorFile, which are alternatives.
		 *
		 * Prescribing g g' and asking for the g that delivers a target q are
		 * opposite statements about the same quantity, and BOTH SPELLINGS
		 * CONVERGE -- to different equilibria. So naming both is refused rather
		 * than resolved by precedence, which is the rule
		 * [ boundary.limiter ] SurfaceAttribute already follows against R / Z
		 * and for the same reason: a key order deciding which physics a run did
		 * is a silent wrong answer.
		 */
		void readToroidalFieldTarget( Table const &source,
		                              MHDParameters &parameters )
		{
			parameters.ggPrimeFile = source.getStringOr( "GGPrimeFile", "" );
			parameters.safetyFactorFile =
				source.getStringOr( "SafetyFactorFile", "" );

			bool const prescribed = !parameters.ggPrimeFile.empty();
			bool const driven = !parameters.safetyFactorFile.empty();

			if ( prescribed && driven )
				source.fail( "SafetyFactorFile",
					"cannot be given beside GGPrimeFile. GGPrimeFile PRESCRIBES "
					"the toroidal field and SafetyFactorFile asks for whichever "
					"field delivers a target q, so the two are alternatives "
					"rather than a pair -- and both converge, to different "
					"equilibria, so there is no safe precedence between them. "
					"Give one" );

			if ( !prescribed && !driven )
				source.fail( "GGPrimeFile",
					"must be given, or SafetyFactorFile in its place. The "
					"source needs g dg/dPsi either as input or as the output of "
					"an outer solve for it" );

			if ( !driven )
				return;

			// ROADMAP.md item 10 drives the loop through the bordered Newton,
			// whose unknowns are psi_ax and psi_bnd, and the target is a
			// function of NORMALISED flux. An un-normalised source has no Psi
			// to write a q table against.
			if ( !parameters.normalised )
				source.fail( "SafetyFactorFile",
					"requires Normalised = true. The target is q( Psi ) against "
					"the normalised flux, and without the normalisation there is "
					"no Psi for it to be a function of" );

			parameters.safetyFactorDegree = static_cast<unsigned int>(
				source.getIntegerOr( "SafetyFactorDegree", 2 ) );
			if ( parameters.safetyFactorDegree < 1
			     || parameters.safetyFactorDegree > 8 )
				source.fail( "SafetyFactorDegree",
					"must be between 1 and 8. Degree zero makes g^2 a constant, "
					"so gg' is identically zero and no coefficient the loop "
					"moves can change the equilibrium; and a degree the surface "
					"family does not determine leaves the outer Jacobian rank "
					"deficient, which is reported rather than solved" );

			parameters.toroidalFieldGuess =
				source.getFloatOr( "ToroidalFieldGuess", 0.0 );
			if ( !( parameters.toroidalFieldGuess > 0.0 )
			     || !std::isfinite( parameters.toroidalFieldGuess ) )
				source.fail( "ToroidalFieldGuess",
					"must be given and positive beside SafetyFactorFile. It is "
					"g = R B_phi to open the loop at, as a constant, and a "
					"machine's vacuum R0 B0 is the number to use. There is no "
					"default: q determines g through the geometry, but only "
					"once there is a geometry" );
		}

		void readEitherProfile( Table const & table, std::string const & key,
		                        double & value, std::string & fileName, double & scale,
		                        bool & given, bool required )
		{
			std::string const fileKey = key + "File";
			std::string const scaleKey = key + "Scale";

			refuseReservedVariableKeys( table, key, true );

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
		/// the reason recorded on refuseReservedVariableKeys.
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

		/// The moving plasma support, which is normalised-only for the same
		/// reason PsiAxis is: the test is Psi > 0, and without a normalisation
		/// there is no Psi. Refused rather than ignored in the un-normalised
		/// case, on the same principle -- a file that asks for a moving plasma
		/// boundary and silently gets a fixed one is a file whose author
		/// believed something false about what the run was doing.
		void readPlasmaSupport( Table const & source, bool normalised, bool & confine )
		{
			if ( !normalised )
			{
				if ( source.has( "ConfineToPlasma" ) )
					source.fail( "ConfineToPlasma", "means nothing unless Normalised = true: the plasma is where the NORMALISED flux Psi = ( psi - psi_bnd )/( psi_ax - psi_bnd ) is positive, and without a normalisation there is no Psi to test" );
				return;
			}

			confine = source.getBooleanOr( "ConfineToPlasma", false );
		}

		/*
		 * `[source] ExcludeAttributes` -- MESH ELEMENT ATTRIBUTES THAT CAN NEVER
		 * BE PLASMA, WHATEVER THE FLUX SAYS THERE.
		 *
		 * A statement about the DEVICE and not about the solution: the far side
		 * of a vessel wall, a port, a pocket the mesh carries for the coils'
		 * sake. `psi_bnd` already confines the plasma and the connectivity fill
		 * already separates the lobes of `{ Psi > 0 }`; what neither can do is
		 * know that a lobe is behind a wall. It earns its place where several
		 * O-points sit across a saddle and connectivity alone cannot say which
		 * of them is the plasma.
		 *
		 * AN ATTRIBUTE AND NOT A POLYGON, because the support is re-decided on
		 * every sweep and a geometric test would be a point-in-polygon per
		 * element per residual evaluation for an answer that cannot move.
		 *
		 * REFUSED WITHOUT A NORMALISATION for ConfineToPlasma's reason: there is
		 * no plasma support to exclude from.
		 */
		void readPlasmaExclusion( Table const & source, bool normalised,
		                          std::vector< int > & attributes )
		{
			if ( !normalised )
			{
				if ( source.has( "ExcludeAttributes" ) )
					source.fail( "ExcludeAttributes", "means nothing unless Normalised = true: it removes elements from the PLASMA SUPPORT, and a source that is not normalised has none -- F is whatever the profiles say everywhere" );
				return;
			}

			attributes = source.getIntegerArrayOr( "ExcludeAttributes" );
			for ( std::size_t i = 0; i < attributes.size(); ++i )
				if ( attributes[ i ] <= 0 )
					source.fail( "ExcludeAttributes[" + std::to_string( i ) + "]",
					             "must be a positive mesh element attribute" );
			for ( std::size_t i = 0; i < attributes.size(); ++i )
				for ( std::size_t j = i + 1; j < attributes.size(); ++j )
					if ( attributes[ i ] == attributes[ j ] )
						source.fail( "ExcludeAttributes", "names the same attribute twice" );
		}

		/// `[source] PlasmaCurrent`, in amperes.
		///
		/// REFUSED WITHOUT A NORMALISATION, for the reason ConfineToPlasma is:
		/// what this makes an unknown is the profile SCALE, and a scale is only
		/// meaningful against profiles read at a normalised flux. On a plain
		/// source the border would be solved and its answer would multiply
		/// nothing.
		void readPlasmaCurrent( Table const & source, bool normalised, double & current )
		{
			if ( !normalised )
			{
				if ( source.has( "PlasmaCurrent" ) )
					source.fail( "PlasmaCurrent", "means nothing unless Normalised = true: prescribing the current makes the profile SCALE an unknown of the bordered Newton, and without a normalisation there is no bordered Newton for it to join" );
				return;
			}

			current = source.getFloatOr( "PlasmaCurrent", 0.0 );

			// ZERO IS THE DEFAULT AND MEANS "DO NOT CONSTRAIN". An explicit zero
			// is refused rather than silently read as the default, because a
			// prescribed current of zero describes no plasma -- the scale would
			// be driven to make int F/R vanish, which it does at scale zero.
			if ( source.has( "PlasmaCurrent" ) && current == 0.0 )
				source.fail( "PlasmaCurrent", "is zero, which is not a plasma: the border would drive the profile scale to zero to satisfy it. Remove the key to leave the amplitude fixed" );
			if ( source.has( "PlasmaCurrent" ) && !std::isfinite( current ) )
				source.fail( "PlasmaCurrent", "must be finite" );
		}

		/*
		 * `[source.axis] R, Z` -- psi_ax AT A PRESCRIBED POINT. See
		 * AxisPointConfig.
		 *
		 * ONE HELPER FOR BOTH NORMALISED SOURCE TYPES, for the reason
		 * refuseReservedVariableKeys() is one: a key that one path reads and
		 * another lists without reading is accepted and ignored on the second.
		 *
		 * A SUB-TABLE, SO ITS KEYS ARE ITS OWN. TOML puts every key after
		 * `[source.axis]` into it until the next header, so a `[source]` key
		 * written below it by mistake lands here -- and rejectUnknownKeys() says
		 * so rather than letting it vanish.
		 */
		void readAxisPoint( Table const & source, std::string const & configFile,
		                    bool normalised, AxisPointConfig & point )
		{
			Table axis( source, "axis", configFile, false );
			axis.rejectUnknownKeys( { "R", "Z" } );
			if ( !axis.isPresent() )
				return;

			if ( !normalised )
				axis.fail( "R", "pins psi_ax at a point, and psi_ax is an unknown only when Normalised = true; without it the profiles are functions of psi itself and there is no psi_ax to pin" );

			// BOTH OR NEITHER, on [boundary.limiter]'s precedent: a point given
			// only its height sits on R = 0, where psi = 0 by construction.
			point.radius = axis.getFloat( "R" );
			point.z = axis.getFloat( "Z" );
			point.given = true;

			if ( !( point.radius > 0.0 ) || !std::isfinite( point.radius ) )
				axis.fail( "R", "the axis label must be at a finite R > 0: psi( 0, z ) = 0 exactly, so a label on the symmetry axis would pin psi_ax to zero whatever the plasma does" );
			if ( !std::isfinite( point.z ) )
				axis.fail( "Z", "must be finite" );
		}

		/*
		 * `[source] SupportShape` and `AxialEdge` -- THE SHAPE OF THE PLASMA
		 * SUPPORT. See SourceConfig::supportShape.
		 *
		 * Every refusal here is a combination that would otherwise be accepted
		 * and ignored, or that poses a problem nothing has measured, and each
		 * is asserted on both source types that can carry the keys.
		 */
		void readSupportShape( Table const & source, bool normalised, bool confine,
		                       bool qDriven, double plasmaCurrent,
		                       SourceConfig & options )
		{
			bool const shapeGiven = source.has( "SupportShape" );
			bool const edgeGiven = source.has( "AxialEdge" );

			if ( !normalised )
			{
				if ( shapeGiven )
					source.fail( "SupportShape", "shapes the PLASMA SUPPORT, which is decided on the normalised flux Psi; without Normalised = true there is no Psi and no support" );
				if ( edgeGiven )
					source.fail( "AxialEdge", "is the axial end of a SupportShape = \"box\" support, and a source that is not Normalised has no support at all" );
				return;
			}

			std::string const shape = source.getStringOr( "SupportShape", "one-sided" );
			if ( shape == "one-sided" )
				options.supportShape = PlasmaSupportShapeType::OneSided;
			else if ( shape == "box" )
				options.supportShape = PlasmaSupportShapeType::Box;
			else
				source.fail( "SupportShape", "must be \"one-sided\" or \"box\", but is \"" + shape + "\"" );

			if ( !options.supportIsBox() )
			{
				// ACCEPTED AND IGNORED IS THE OUTCOME THIS FILE REFUSES: a
				// one-sided support has no axial end, so an edge beside it is a
				// number nothing reads.
				if ( edgeGiven )
					source.fail( "AxialEdge", "is the axial end of a box support and means nothing without SupportShape = \"box\"; the one-sided support is Psi > 0 alone" );
				return;
			}

			// A SHAPE OF THE CONFINEMENT, so it does not switch confinement on
			// by implication: every consumer of confinesToPlasma() reads one
			// flag, and a support with a shape but no switch is a contradiction
			// to refuse rather than a default to guess.
			if ( !confine )
				source.fail( "SupportShape", "\"box\" is a shape of the plasma support, and the support is off: set ConfineToPlasma = true as well" );

			if ( !edgeGiven )
				source.fail( "AxialEdge", "is required with SupportShape = \"box\": the box stops the plasma at the throat, |z| <= AxialEdge, and a mirror support without it is the flux-only mask M-168 measured putting the end-plane current five orders above the midplane's" );
			options.axialEdge = source.getFloat( "AxialEdge" );
			if ( !( options.axialEdge > 0.0 ) || !std::isfinite( options.axialEdge ) )
				source.fail( "AxialEdge", "must be finite and positive: it is the throat, in metres from z = 0" );

			// THE BOX'S UPPER BOUND, Psi <= 1, IS A RING CUT OFF A TOKAMAK'S
			// CORE: psi_h overshoots a located psi_ax by a discretisation error.
			// So the box is posed only with the inner label at a point.
			if ( !options.axisPoint.given )
				source.fail( "SupportShape", "\"box\" bounds the support by Psi <= 1, which is right between a mirror's two prescribed labels and cuts a ring off the core around a LOCATED axis; give the inner label as [source.axis] R and Z" );

			// INERT UNDER THE BOX, SO REFUSED. Exclusion acts through the
			// element mask the connectivity fill builds, and the fill is off
			// under the box by a stated branch -- so the attributes would be
			// read and would remove nothing.
			if ( !options.excludeAttributes.empty() )
				source.fail( "ExcludeAttributes", "acts through the connectivity fill's element mask, and SupportShape = \"box\" switches the fill off: the box is pointwise. Remove one of the two" );

			// The q-driven loop extracts flux surfaces about an O-point every
			// map evaluation, and a mirror has none.
			if ( qDriven )
				source.fail( "SafetyFactorFile", "drives g g' through q on closed flux surfaces about a magnetic axis, and a SupportShape = \"box\" plasma has neither" );

			// NOT MEASURED ON THE BOX, AND THE CURRENT BORDER'S ROW IS
			// ASSEMBLED ON THE CONNECTIVITY FILL'S SUPPORT, which the box
			// switches off. Refused until a measured case lifts it.
			if ( plasmaCurrent != 0.0 )
				source.fail( "PlasmaCurrent", "is not measured beside SupportShape = \"box\": the current border has not been solved on a mirror. Remove it, and the profiles' own amplitude sets the current" );
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

	double SourceConfig::plasmaCurrent() const
	{
		switch ( type )
		{
			case SourceType::MHD:
				return std::get< MHDParameters >( parameters ).plasmaCurrent;
			case SourceType::Rotating:
				return std::get< RotatingParameters >( parameters ).plasmaCurrent;
			case SourceType::Soloviev:
			case SourceType::Manufactured:
				return 0.0;
		}

		return 0.0;
	}

	bool SourceConfig::confinesToPlasma() const
	{
		switch ( type )
		{
			case SourceType::MHD:
				return std::get< MHDParameters >( parameters ).confineToPlasma;
			case SourceType::Rotating:
				return std::get< RotatingParameters >( parameters ).confineToPlasma;
			case SourceType::Soloviev:
			case SourceType::Manufactured:
				return false;
		}

		return false;
	}

	double SourceConfig::permeability() const
	{
		switch ( type )
		{
			case SourceType::MHD:
				return std::get< MHDParameters >( parameters ).mu0;
			case SourceType::Rotating:
				return std::get< RotatingParameters >( parameters ).mu0;
			case SourceType::Soloviev:
			case SourceType::Manufactured:
				break;
		}

		return 4.0e-7*3.14159265358979323846;
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

	std::string OutputConfig::getFluxSurfaceFile() const
	{
		return directory.empty() ? prefix + "_surfaces.nc" : directory + "/" + prefix + "_surfaces.nc";
	}

	Configuration::Configuration( std::string const & fileName )
		: Configuration( parseFile( fileName ), fileName )
	{
		resolveInputPaths( std::filesystem::path( fileName ).parent_path().string() );
	}

	void Configuration::resolveInputPaths( std::string const & directory )
	{
		if ( directory.empty() )
			return;

		auto resolve = [ &directory ]( std::string & path )
		{
			if ( path.empty() || std::filesystem::path( path ).is_absolute() )
				return;
			path = ( std::filesystem::path( directory ) / path ).lexically_normal().string();
		};

		// A generated mesh is the run's product, not the file's input.
		if ( !meshOptions.generate.given )
			resolve( meshOptions.file );

		if ( auto * mhd = std::get_if<MHDParameters>( &sourceOptions.parameters ) )
		{
			resolve( mhd->pPrimeFile );
			resolve( mhd->ggPrimeFile );
			resolve( mhd->safetyFactorFile );
		}
		if ( auto * rotating = std::get_if<RotatingParameters>( &sourceOptions.parameters ) )
		{
			resolve( rotating->omegaFile );
			resolve( rotating->ggPrimeFile );
			for ( SpeciesParameters & one : rotating->species )
			{
				resolve( one.temperatureFile );
				resolve( one.densityFile );
			}
		}

		resolve( initialGuessOptions.file );
		resolve( initialGuessOptions.meshFile );
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
			std::initializer_list< char const * > const tables = { "mesh", "discretisation", "source", "boundary", "solver", "output", "initialguess", "adaptivity", "coils", "circuits", "conductors", "inverse" };
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

		// [[coils]]
		//
		// An ARRAY OF TABLES, MEQ's second after [[source.species]], and read
		// the same way: getTableArrayOr() names its elements "coils[i]" so a
		// fault in the third one says so rather than reporting a key missing
		// from a table the author never wrote.
		//
		// EMPTY IS LEGAL AND IS THE COMMON CASE. Every fixed-boundary
		// configuration in examples/ has no coils, so absence is not an error;
		// what is an error is a coil that cannot be built.
		{
			Table const rootTable = Table::root( document, sourceName );
			std::vector< Table > const blocks = rootTable.getTableArrayOr( "coils" );
			std::vector< bool > namesCurrent;
			for ( std::size_t i = 0; i < blocks.size(); ++i )
			{
				Table const & one = blocks[ i ];
				one.rejectUnknownKeys( { "Name", "CentreR", "CentreZ", "HalfWidth",
				                         "HalfHeight", "Current", "CurrentDensity",
				                         "FilamentsR", "FilamentsZ", "Free" } );

				CoilParameters coil;
				coil.name = one.getStringOr( "Name", "coil" + std::to_string( i ) );
				coil.centreR = one.getFloat( "CentreR" );
				coil.centreZ = one.getFloat( "CentreZ" );
				coil.halfWidth = one.getFloat( "HalfWidth" );
				coil.halfHeight = one.getFloat( "HalfHeight" );

				if ( !( coil.halfWidth > 0.0 ) )
					one.fail( "HalfWidth", "must be strictly positive: a coil of zero width has no cross-section to carry a current density over" );
				if ( !( coil.halfHeight > 0.0 ) )
					one.fail( "HalfHeight", "must be strictly positive: a coil of zero height has no cross-section to carry a current density over" );

				// THE SAME REFUSAL meq::Coil AND meq::BoundaryShape MAKE, and
				// for the same reason: the Grad-Shafranov operator carries a
				// 1/R that is not integrable through R = 0, so a conductor
				// reaching the axis is not a modelling choice this code can
				// honour. Caught here as well as there so that the diagnostic
				// names the coil and the key rather than arriving from a
				// constructor three layers down.
				if ( !( coil.centreR - coil.halfWidth > 0.0 ) )
					one.fail( "CentreR", "the coil reaches or crosses the axis: CentreR - HalfWidth = "
					          + std::to_string( coil.centreR - coil.halfWidth )
					          + " and must be strictly positive, because the operator's 1/R is not integrable through R = 0" );

				bool const hasCurrent = one.has( "Current" );
				bool const hasDensity = one.has( "CurrentDensity" );

				// EXACTLY ONE. Both is refused rather than resolved by
				// precedence -- an author who writes both has two numbers in
				// mind, and silently honouring one is how a coil set ends up
				// carrying a current nobody chose.
				if ( hasCurrent && hasDensity )
					one.fail( "Current", "names both Current and CurrentDensity; give exactly one. They are related by the area 4 * HalfWidth * HalfHeight = "
					          + std::to_string( 4.0*coil.halfWidth*coil.halfHeight )
					          + " m^2, so writing both says the same thing twice or contradicts itself, and this cannot tell which" );
				// NEITHER IS REFUSED BELOW, once [[circuits]] is read: a coil
				// wired into a circuit takes its current from the circuit and
				// must name none, so whether "no current" is an error depends
				// on a table this loop has not reached.
				namesCurrent.push_back( hasCurrent || hasDensity );

				double const area = 4.0*coil.halfWidth*coil.halfHeight;
				if ( hasDensity )
				{
					coil.densityGiven = true;
					coil.current = one.getFloat( "CurrentDensity" )*area;
				}
				else if ( hasCurrent )
				{
					coil.current = one.getFloat( "Current" );
				}

				// AN UNKNOWN, NOT DATA: `Current` is where the inverse solve
				// starts, which is why it is still required beside it.
				coil.free = one.getBooleanOr( "Free", false );

				/*
				 * THIS BLOCK'S OWN FILAMENT STACK, overriding [conductors]
				 * FilamentSize.
				 *
				 * BOTH OR NEITHER. An author who names one has a number in
				 * mind for that direction and none for the other, and taking
				 * the missing one from the size would leave a block's stack
				 * half explicit and half derived with nothing saying which.
				 * Refused rather than filled in, exactly as Current beside
				 * CurrentDensity is -- and 1 is a legal count, so the
				 * instruction the message gives is always writable.
				 *
				 * WHETHER THE MODEL CAN HONOUR THEM IS NOT ASKED HERE, because
				 * [conductors] has not been read yet: that refusal is below,
				 * where both halves are in hand.
				 */
				bool const hasFilamentsR = one.has( "FilamentsR" );
				bool const hasFilamentsZ = one.has( "FilamentsZ" );
				if ( hasFilamentsR != hasFilamentsZ )
					one.fail( hasFilamentsR ? "FilamentsR" : "FilamentsZ",
					          "names one of FilamentsR and FilamentsZ and not the other; give both. They are a stack in two directions and there is no default for the one left out -- taking it from [conductors] FilamentSize would make this block's stack half explicit and half derived, with the file saying neither. 1 is a legal count and is the single filament in that direction" );

				if ( hasFilamentsR )
				{
					coil.filamentsR = one.getInteger( "FilamentsR" );
					coil.filamentsZ = one.getInteger( "FilamentsZ" );

					if ( coil.filamentsR < 1 || coil.filamentsZ < 1 )
						one.fail( coil.filamentsR < 1 ? "FilamentsR" : "FilamentsZ",
						          "must be at least 1; 1 x 1 is the single filament at the rectangle's centre, which is what a block without these keys gets, and there is no reading of a count of zero that is a conductor" );

					// The library's own cap, caught here so the diagnostic
					// names the block and the keys rather than arriving from
					// meq::filamentStack three layers down.
					long long const total = static_cast<long long>( coil.filamentsR )
					                        *static_cast<long long>( coil.filamentsZ );
					if ( total > maximumFilamentsPerCoil )
						one.fail( "FilamentsR", "FilamentsR x FilamentsZ = "
						          + std::to_string( total ) + " filaments for one rectangle, against a limit of "
						          + std::to_string( maximumFilamentsPerCoil )
						          + ". A stack this fine is a cell size given in the wrong units far more often than a modelling choice; if it is the choice, Model = \"subtracted\" integrates the same rectangle exactly at a fixed cost" );
				}

				coilOptions.coils.push_back( coil );
			}

			/*
			 * [[circuits]] -- CONDUCTORS ON ONE SUPPLY.
			 * INVERSE-AND-RECONSTRUCTION-PLAN.md section 7.2, and see
			 * CircuitParameters for why a flat coil list is not enough.
			 *
			 * Every member is a [[coils]] block BY NAME, so a member must exist
			 * and be the only block of that name; it belongs to one circuit;
			 * and it names no current or Free of its own, since the circuit's
			 * current is the one number and the circuit's Free the one
			 * unknown. Each of those is refused rather than resolved by
			 * precedence, for the reason Current beside CurrentDensity is.
			 *
			 * THE MULTIPLIER IS REQUIRED. It is signed -- an anti-series pair
			 * is a -1 -- and a default of 1 would wire a pair in series that
			 * the author meant the other way, silently and with a plausible
			 * equilibrium at the end.
			 */
			std::vector< Table > const circuitBlocks = rootTable.getTableArrayOr( "circuits" );
			std::vector< std::string > circuitNames;
			for ( std::size_t c = 0; c < circuitBlocks.size(); ++c )
			{
				Table const & one = circuitBlocks[ c ];
				one.rejectUnknownKeys( { "Name", "Current", "Free", "Coils" } );

				CircuitParameters circuit;
				circuit.name = one.getString( "Name" );
				// "default" IS THE OUTPUT'S NAME FOR THE UNWIRED COILS, a
				// grouping with no supply (meq::defaultCircuitName). A supplied
				// circuit of that name would merge with them in the file.
				if ( circuit.name == defaultCircuitName )
					one.fail( "Name", "\"" + circuit.name + "\" is reserved: the output groups every coil outside a [[circuits]] block into a circuit of that name, which has no supply, and a supplied circuit called the same would merge with it. Choose another name" );
				if ( std::find( circuitNames.begin(), circuitNames.end(), circuit.name )
				     != circuitNames.end() )
					one.fail( "Name", "\"" + circuit.name + "\" names a second circuit; circuit names must be distinct, because they are how a result reports which current it chose" );
				circuitNames.push_back( circuit.name );

				if ( !one.has( "Current" ) )
					one.fail( "Current", "a circuit needs a current: the TOTAL its supply drives, in amperes, which each member carries times its Multiplier -- and the starting value when Free = true" );
				circuit.current = one.getFloat( "Current" );
				circuit.free = one.getBooleanOr( "Free", false );

				std::vector< Table > const members = one.getTableArrayOr( "Coils" );
				if ( members.empty() )
					one.fail( "Coils", "a circuit needs at least one member: Coils = [ { Name = \"...\", Multiplier = 1.0 }, ... ], each Name a [[coils]] block" );

				for ( Table const & member : members )
				{
					member.rejectUnknownKeys( { "Name", "Multiplier" } );
					CircuitMember entry;
					entry.coil = member.getString( "Name" );
					if ( !member.has( "Multiplier" ) )
						member.fail( "Multiplier", "is required: it is this coil's signed share of the circuit's current, and a default of 1 would wire an anti-series pair in series without saying so" );
					entry.multiplier = member.getFloat( "Multiplier" );
					if ( !std::isfinite( entry.multiplier ) || entry.multiplier == 0.0 )
						member.fail( "Multiplier", "must be finite and non-zero; a coil that carries none of the circuit's current is not in the circuit" );

					std::size_t found = coilOptions.coils.size();
					for ( std::size_t k = 0; k < coilOptions.coils.size(); ++k )
					{
						if ( coilOptions.coils[ k ].name != entry.coil )
							continue;
						if ( found != coilOptions.coils.size() )
							member.fail( "Name", "\"" + entry.coil + "\" names more than one [[coils]] block, so this cannot tell which is in the circuit; give the coils distinct names" );
						found = k;
					}
					if ( found == coilOptions.coils.size() )
						member.fail( "Name", "\"" + entry.coil + "\" is not the Name of any [[coils]] block" );

					CoilParameters & coil = coilOptions.coils[ found ];
					Table const & coilTable = blocks[ found ];
					if ( !coil.circuit.empty() )
						member.fail( "Name", "\"" + entry.coil + "\" is already in circuit \"" + coil.circuit + "\"; a conductor is on one supply" );
					if ( namesCurrent[ found ] )
						coilTable.fail( coilTable.has( "Current" ) ? "Current" : "CurrentDensity",
						                "this coil is in circuit \"" + circuit.name + "\", whose Current times this coil's Multiplier IS its current; naming one here as well says the same thing twice or contradicts it, and this cannot tell which" );
					if ( coilTable.has( "Free" ) )
						coilTable.fail( "Free", "this coil is in circuit \"" + circuit.name + "\", and a circuit has one current: write Free on the [[circuits]] block, where the unknown is" );

					coil.circuit = circuit.name;
					coil.multiplier = entry.multiplier;
					coil.current = entry.multiplier*circuit.current;
					entry.index = found;
					circuit.members.push_back( entry );
				}

				coilOptions.circuits.push_back( circuit );
			}

			// AND NOW "NO CURRENT" CAN BE DECIDED: a coil outside every
			// circuit still needs exactly one, with the message it always had.
			for ( std::size_t i = 0; i < blocks.size(); ++i )
				if ( !namesCurrent[ i ] && coilOptions.coils[ i ].circuit.empty() )
					blocks[ i ].fail( "Current", "a coil needs a current: give either Current, the TOTAL through the cross-section in amperes, or CurrentDensity, the uniform j_phi in A/m^2 -- or wire it into a [[circuits]] block, which then supplies it" );
		}

		// [conductors] -- HOW the blocks above enter the equation.
		//
		// READ AFTER [[coils]] BECAUSE IT IS ABOUT THEM, and the refusal below
		// needs to know whether there are any. A table with no [[coils]] to
		// act on is refused rather than accepted and ignored: an
		// accepted-and-ignored key is exactly the failure CLAUDE.md records
		// twice already -- ConfineToPlasma on the rotating source, and the
		// reserved profile keys on "mhd" -- and this one would be worse, since
		// a user who asked for the split and got the meshed route would read a
		// run that is slower than it should be rather than one that is wrong,
		// and would have nothing to look at.
		{
			Table conductors( document, "conductors", sourceName, false );
			conductors.rejectUnknownKeys( { "Model", "QuadratureOrder",
			                                "FilamentSize", "CacheFile",
			                                "InterpolatedFlux" } );

			std::string const model = conductors.getStringOr( "Model", "meshed" );
			if ( model == "meshed" )
				conductorOptions.model = ConductorModel::Meshed;
			else if ( model == "subtracted" )
				conductorOptions.model = ConductorModel::Subtracted;
			else if ( model == "filament" )
				conductorOptions.model = ConductorModel::Filament;
			else
				conductors.fail( "Model", "must be meshed, subtracted or filament, "
				                 "but is \"" + model + "\". \"meshed\" carries each "
				                 "[[coils]] rectangle as a domain source and needs the "
				                 "mesh to resolve it; \"subtracted\" integrates the same "
				                 "rectangle analytically and solves for psi - psi_c, so "
				                 "the mesh need not carry it; \"filament\" does the same "
				                 "with a point filament at each rectangle's centre, "
				                 "which is a different and cheaper conductor model "
				                 "rather than the same one computed differently" );

			conductorOptions.quadratureOrder =
				conductors.getIntegerOr( "QuadratureOrder",
				                         conductorOptions.quadratureOrder );

			if ( conductorOptions.subtracts() && coilOptions.coils.empty() )
				conductors.fail( "Model", "= \"" + model + "\" takes the conductors out "
				                 "of the mesh, and this file describes none: there are "
				                 "no [[coils]] blocks for it to act on. Written as it "
				                 "stands the key would be accepted and do nothing" );

			if ( conductorOptions.quadratureOrder != 0 )
			{
				if ( conductorOptions.model == ConductorModel::Filament )
					conductors.fail( "QuadratureOrder", "a filament has no "
					                 "cross-section to integrate over -- its field is "
					                 "one closed-form evaluation -- so this key would "
					                 "be accepted and ignored under Model = "
					                 "\"filament\"" );
				if ( conductorOptions.model == ConductorModel::Meshed )
					conductors.fail( "QuadratureOrder", "this is the quadrature "
					                 "meq::ConductorField integrates a SUBTRACTED "
					                 "rectangle's cross-section with, and Model = "
					                 "\"meshed\" subtracts nothing: the mesh's own "
					                 "quadrature is what integrates a meshed coil, and "
					                 "[discretisation] is where that is set" );
				if ( conductorOptions.quadratureOrder < 2 )
					conductors.fail( "QuadratureOrder", "must be at least 2 Gauss "
					                 "points per direction; omit the key for "
					                 "meq::ConductorField's own default" );
			}

			/*
			 * THE FILAMENT STACK, AND BOTH HALVES OF IT ARE CHECKED HERE.
			 *
			 * [conductors] FilamentSize sets every block's stack by a target
			 * CELL size; a [[coils]] block naming FilamentsR and FilamentsZ
			 * overrides it for that block alone. Both are meaningless under
			 * the other two models and are refused rather than accepted and
			 * ignored -- CLAUDE.md records that failure three times now, and
			 * this one would be the quiet kind: a file asking for a subdivided
			 * solenoid and getting one filament at its centre reads as a
			 * converged machine that is not the machine asked for.
			 *
			 * THE PER-BLOCK KEYS ARE CHECKED HERE AND NOT IN THE LOOP ABOVE
			 * because the loop runs before this table is read, so it cannot
			 * know which model is in force.
			 */
			conductorOptions.filamentSize =
				conductors.getFloatOr( "FilamentSize",
				                       conductorOptions.filamentSize );

			std::size_t stackedBlock = coilOptions.coils.size();
			for ( std::size_t i = 0; i < coilOptions.coils.size(); ++i )
				if ( coilOptions.coils[ i ].stackGiven() )
				{
					stackedBlock = i;
					break;
				}
			bool const anyBlockStacked = stackedBlock < coilOptions.coils.size();

			if ( conductorOptions.model != ConductorModel::Filament )
			{
				if ( conductorOptions.filamentSize != 0.0 )
					conductors.fail( "FilamentSize", "divides each [[coils]] rectangle into a stack of point filaments, which only Model = \"filament\" carries. \"meshed\" puts the rectangle in the mesh and \"subtracted\" integrates it exactly -- both of which already resolve the cross-section, so a stack would be a worse answer at more cost" );

				if ( anyBlockStacked )
					conductors.fail( "Model", "= \"" + model + "\" does not divide a rectangle into filaments, and [[coils]] block "
					                 + std::to_string( stackedBlock ) + " ("
					                 + coilOptions.coils[ stackedBlock ].name
					                 + ") names FilamentsR and FilamentsZ. Set Model = \"filament\" or remove the keys; written as it stands they would be accepted and do nothing" );
			}
			/*
			 * REFUSED UNDER "meshed" RATHER THAN ACCEPTED AND IGNORED. A
			 * meshed run subtracts nothing, so there is no psi_c to cache and
			 * the key would name a file that is never written and never read
			 * -- the accepted-and-ignored failure this file records under the
			 * reserved profile keys, and which recurred once already with
			 * ConfineToPlasma on the rotating source.
			 */
			conductorOptions.cacheFile =
				conductors.getStringOr( "CacheFile",
				                        conductorOptions.cacheFile );

			conductorOptions.interpolatedFlux =
				conductors.getBooleanOr( "InterpolatedFlux",
				                      conductorOptions.interpolatedFlux );

			if ( conductorOptions.interpolatedFlux
			     && !conductorOptions.subtracts() )
				conductors.fail( "InterpolatedFlux", "interpolates q_c, the field the SUBTRACTED conductors contribute, and Model = \"meshed\" puts the conductors in the mesh instead -- so there is no q_c for this run to interpolate. Set Model = \"subtracted\" or \"filament\", or remove the key" );

			if ( !conductorOptions.cacheFile.empty()
			     && !conductorOptions.subtracts() )
				conductors.fail( "CacheFile", "caches psi_c, the field the SUBTRACTED conductors contribute, and Model = \"meshed\" puts the conductors in the mesh instead -- so there is no psi_c for this run to cache. Set Model = \"subtracted\" or \"filament\", or remove the key" );

			if ( conductorOptions.model == ConductorModel::Filament
			     && conductorOptions.filamentSize != 0.0
			     && !( conductorOptions.filamentSize > 0.0 ) )
				conductors.fail( "FilamentSize", "must be strictly positive, in metres: it is an upper bound on a filament cell's extent, so each block gets ceil( 2*HalfWidth/FilamentSize ) x ceil( 2*HalfHeight/FilamentSize ) filaments. Omit the key for one filament per block, which is what a file without it means" );
		}

		// [mesh]
		{
			Table mesh( document, "mesh", sourceName, true );
			mesh.rejectUnknownKeys( { "RMin", "RMax", "ZMin", "ZMax", "NR", "NZ", "RefinementLevels", "File", "generate" } );

			meshOptions.file = mesh.getStringOr( "File", "" );
			meshOptions.refinementLevels = mesh.getIntegerOr( "RefinementLevels", 0 );

			if ( meshOptions.refinementLevels < 0 )
				mesh.fail( "RefinementLevels", "must not be negative" );

			if ( !meshOptions.fromFile() )
			{
				meshOptions.minRadius = mesh.getFloat( "RMin" );
				meshOptions.maxRadius = mesh.getFloat( "RMax" );
				meshOptions.zMin = mesh.getFloat( "ZMin" );
				meshOptions.zMax = mesh.getFloat( "ZMax" );
				meshOptions.nR = mesh.getIntegerOr( "NR", 1 );
				meshOptions.nZ = mesh.getIntegerOr( "NZ", 1 );

				if ( meshOptions.minRadius < 0.0 )
					mesh.fail( "RMin", "must not be negative: R is a cylindrical radius" );
				if ( meshOptions.maxRadius <= meshOptions.minRadius )
					mesh.fail( "RMax", "must be greater than RMin (RMin = " + std::to_string( meshOptions.minRadius ) + ")" );
				if ( meshOptions.zMax <= meshOptions.zMin )
					mesh.fail( "ZMax", "must be greater than ZMin (ZMin = " + std::to_string( meshOptions.zMin ) + ")" );
				if ( meshOptions.nR < 1 )
					mesh.fail( "NR", "must be at least 1" );
				if ( meshOptions.nZ < 1 )
					mesh.fail( "NZ", "must be at least 1" );
			}

			// [mesh.generate] -- the mesh as a BUILD PRODUCT of this file.
			//
			// MEQ does not run the generator: `meq` links MFEM and not gmsh,
			// for the reasons tools/mesh/README.md records. What it does with
			// this block is print the generator's argument list on request
			// (`meq --mesh-command`) and REFUSE to solve without
			// `--mesh-ready`, so that an edited geometry cannot be answered
			// from the mesh the previous geometry made. `meq-run` is the
			// caller that does both, and it is the single executable the user
			// runs.
			//
			// THE COILS ARE NOT HERE. They come from the [[coils]] blocks, in
			// file order, which is the order halfdisc.py assigns its 10 + i
			// element attributes in. Writing a machine's conductors twice --
			// once on a command line and once in the file the solve reads --
			// is what this block exists to stop.
			{
				Table generate( mesh, "generate", sourceName, false );
				generate.rejectUnknownKeys( { "Tool", "Radius", "Size", "Order", "CoilSize",
				                              "PlasmaRMin", "PlasmaRMax", "PlasmaZMin", "PlasmaZMax",
				                              "PlasmaSize", "LimiterR", "LimiterZ", "LimiterRadius",
				                              "Vessel", "Transition", "Check",
				                              "Symmetric", "AlignPlasma", "AlignConductors" } );

				MeshGeneratorConfig & g = meshOptions.generate;
				g.given = generate.has( "Tool" );

				// TOOL IS THE BLOCK, AND THE TEST IS "ANY KEY AT ALL" RATHER
				// THAN A LIST. A geometry with no generator named is a set of
				// numbers nothing will read -- which is the accepted-and-
				// ignored failure this schema refuses everywhere else -- and
				// reporting it against the key that is MISSING is more useful
				// than reporting it against whichever one happens to be first.
				if ( !g.given && !generate.isEmpty() )
					generate.fail( "Tool", "[mesh.generate] needs a Tool to name which generator makes the mesh; \"halfdisc\" is the one there is" );

				if ( g.given )
				{
					g.tool = generate.getString( "Tool" );
					if ( g.tool != "halfdisc" )
						generate.fail( "Tool", "unknown mesh generator \"" + g.tool + "\"; the generators are [halfdisc], which is tools/mesh/halfdisc.py -- a semicircle reaching R = 0 exactly, with the conductors fragmented in" );

					// A generator with nowhere to write is not a run. File is
					// also what the solve then READS, so the two are the same
					// path by construction rather than by the author keeping
					// two lines in step.
					if ( !meshOptions.fromFile() )
						generate.fail( "Tool", "[mesh.generate] makes the mesh [mesh] File names, so File is required with it: it is where the generator writes and where the solve reads" );

					g.radius = generate.getFloat( "Radius" );
					g.size = generate.getFloat( "Size" );
					g.order = generate.getIntegerOr( "Order", 1 );
					g.coilSize = generate.getFloatOr( "CoilSize", 0.0 );
					g.transition = generate.getFloatOr( "Transition", 0.0 );
					g.check = generate.getBooleanOr( "Check", true );
					g.symmetric = generate.getBooleanOr( "Symmetric", false );
					g.alignPlasma = generate.getBooleanOr( "AlignPlasma", false );
					g.alignConductors = generate.getBooleanOr( "AlignConductors", false );

					if ( !( g.radius > 0.0 ) )
						generate.fail( "Radius", "the disc's radius must be positive. NOTE it is the BACKGROUND's radius and not Gamma's: with an exterior coupling, D_h is cut from this mesh at [boundary.exterior] Radius, which has to fit strictly inside it" );
					if ( !( g.size > 0.0 ) )
						generate.fail( "Size", "the background element size must be positive" );
					if ( g.order < 1 )
						generate.fail( "Order", "the geometric order must be at least 1; above 1 the arc's mid-edge nodes are placed on the true circle rather than on the chord" );
					if ( g.coilSize < 0.0 )
						generate.fail( "CoilSize", "the element size inside the conductors must not be negative; omit it for the background Size" );

					// COILSIZE GRADES THE MESH AROUND A --coil, AND A
					// SUBTRACTED CONDUCTOR IS NEVER EMITTED AS ONE.
					//
					// The driver's mesh command omits every --coil under a
					// subtracting [conductors] Model -- that is the whole of
					// how the key reaches the geometry -- so a CoilSize beside
					// it names an element size for a region the generator is
					// never told about. Accepted and ignored is the one
					// outcome this schema refuses everywhere, and here it is
					// the outcome a user would MOST easily misread: the point
					// of the split is that the mesh gets cheaper, and a key
					// that looks like it is still paying for the conductors is
					// exactly the wrong thing to leave sitting in the file.
					if ( generate.has( "CoilSize" ) && conductorOptions.subtracts() )
						generate.fail( "CoilSize", "[conductors] Model takes the conductors OUT of the mesh, so the generator is never told where they are and there is nothing for this element size to grade. Remove the key -- the background Size is what meshes that region now -- or set Model = \"meshed\" to carry the conductors as a domain source again" );
					// ALIGNCONDUCTORS IS THE SUBTRACTING MODELS' --coil, AS
					// GEOMETRY ONLY. Under "meshed" every rectangle is already
					// fragmented in with its attribute, so the key would be
					// accepted and do nothing. A subtracting model with no
					// [[coils]] is refused above, so there is always something
					// to align when this passes.
					if ( g.alignConductors && !conductorOptions.subtracts() )
						generate.fail( "AlignConductors", "aligns the mesh to the [[coils]] rectangles under a subtracting [conductors] Model, and Model = \"meshed\" already fragments every one of them in as a --coil with its own attribute. Remove the key, or set Model = \"subtracted\" or \"filament\"" );
					if ( g.transition < 0.0 )
						generate.fail( "Transition", "the graded transition's width must not be negative; omit it for four background sizes" );

					// THE REFINED BOX, IN MEQ's OWN CONVENTION. halfdisc.py
					// takes a corner and two extents and [mesh] takes four
					// bounds; the driver converts when it prints the command,
					// so one file never carries two meanings of four numbers.
					bool const boxGiven = generate.has( "PlasmaRMin" ) || generate.has( "PlasmaRMax" )
					                      || generate.has( "PlasmaZMin" ) || generate.has( "PlasmaZMax" );
					g.plasmaGiven = boxGiven || generate.has( "PlasmaSize" );
					if ( g.plasmaGiven )
					{
						// BOTH OR NEITHER, for halfdisc.py's own reason: a
						// region with no size refines nothing and a size with
						// no region has nowhere to act.
						g.plasmaRMin = generate.getFloat( "PlasmaRMin" );
						g.plasmaRMax = generate.getFloat( "PlasmaRMax" );
						g.plasmaZMin = generate.getFloat( "PlasmaZMin" );
						g.plasmaZMax = generate.getFloat( "PlasmaZMax" );
						g.plasmaSize = generate.getFloat( "PlasmaSize" );

						if ( g.plasmaRMin < 0.0 )
							generate.fail( "PlasmaRMin", "must not be negative: R is a cylindrical radius" );
						if ( g.plasmaRMax <= g.plasmaRMin )
							generate.fail( "PlasmaRMax", "must be greater than PlasmaRMin" );
						if ( g.plasmaZMax <= g.plasmaZMin )
							generate.fail( "PlasmaZMax", "must be greater than PlasmaZMin" );
						if ( !( g.plasmaSize > 0.0 ) )
							generate.fail( "PlasmaSize", "the element size in the refined box must be positive" );
					}

					/*
					 * `AlignPlasma` -- THE BOX AS GEOMETRY, and halfdisc.py's own
					 * refusals made again at parse time, where they cost nothing:
					 * there is a box to align, it stays inside the arc, and
					 * under Symmetric it is its own reflection, since aligned it
					 * is geometry rather than a size field.
					 */
					if ( g.alignPlasma )
					{
						if ( !g.plasmaGiven )
							generate.fail( "AlignPlasma", "aligns the mesh to the Plasma* box, and there is none: give PlasmaRMin, PlasmaRMax, PlasmaZMin, PlasmaZMax and PlasmaSize" );
						if ( std::hypot( g.plasmaRMax, std::max( std::abs( g.plasmaZMin ), std::abs( g.plasmaZMax ) ) ) >= g.radius )
							generate.fail( "AlignPlasma", "the Plasma* box reaches the disc's arc at Radius = " + std::to_string( g.radius ) + ", and an aligned box is geometry that has to lie inside it" );
						if ( g.symmetric && g.plasmaZMin != -g.plasmaZMax )
							generate.fail( "AlignPlasma", "makes the Plasma* box geometry, so Symmetric = true needs it symmetric about z = 0, PlasmaZMin = -PlasmaZMax; it is z in [ " + std::to_string( g.plasmaZMin ) + ", " + std::to_string( g.plasmaZMax ) + " ]" );
					}

					// THE LIMITER, FRAGMENTED IN RATHER THAN CUT. Written as
					// element attribute 20, which is what [boundary.limiter]
					// SurfaceAttribute reads; the two are checked against each
					// other where that block is parsed.
					g.limiterGiven = generate.has( "LimiterR" ) || generate.has( "LimiterZ" )
					                 || generate.has( "LimiterRadius" );
					if ( g.limiterGiven )
					{
						g.limiterR = generate.getFloat( "LimiterR" );
						g.limiterZ = generate.getFloat( "LimiterZ" );
						g.limiterRadius = generate.getFloat( "LimiterRadius" );

						if ( !( g.limiterRadius > 0.0 ) )
							generate.fail( "LimiterRadius", "a limiter of zero radius has no interior to give an attribute to" );
						if ( !( g.limiterR - g.limiterRadius > 0.0 ) )
							generate.fail( "LimiterR", "the limiter circle reaches or crosses the axis: LimiterR - LimiterRadius = "
							               + std::to_string( g.limiterR - g.limiterRadius )
							               + " and must be strictly positive, because a closed plasma surface through R = 0 carries a non-integrable 1/R" );
					}

					/*
					 * `Vessel` -- THE REGION THAT CAN NEVER BE PLASMA, CHECKED
					 * HERE BECAUSE A SOLVE CANNOT CHECK IT AT ALL.
					 *
					 * The mesher writes attribute 30 outside this polygon and
					 * `[source] ExcludeAttributes` names it; a polygon with
					 * fewer than three points, an odd count of numbers, a
					 * negative radius or zero area produces a mesh that is
					 * wrong in a way the solve reads as an ordinary geometry.
					 * Putting the geometry in the file is what turns those into
					 * parse errors, which is this block's whole argument.
					 */
					g.vessel = generate.getFloatArrayOr( "Vessel" );
					if ( !g.vessel.empty() )
					{
						if ( g.vessel.size() % 2 != 0 )
							generate.fail( "Vessel", "is alternating R and Z, so it takes an EVEN count of numbers; this has "
							               + std::to_string( g.vessel.size() ) );
						if ( g.vessel.size() < 6 )
							generate.fail( "Vessel", "is a closed polygon and needs at least three points, so at least six numbers; this has "
							               + std::to_string( g.vessel.size() ) );
						double area = 0.0;
						for ( std::size_t k = 0; k < g.vessel.size(); k += 2 )
						{
							if ( g.vessel[ k ] < 0.0 )
								generate.fail( "Vessel", "reaches R = "
								               + std::to_string( g.vessel[ k ] )
								               + ", and the domain is R >= 0" );
							std::size_t const n = ( k + 2 ) % g.vessel.size();
							area += g.vessel[ k ]*g.vessel[ n + 1 ]
							        - g.vessel[ n ]*g.vessel[ k + 1 ];
						}
						if ( std::abs( area )/2.0 <= 0.0 )
							generate.fail( "Vessel", "has zero area: its points are collinear or repeated" );
					}
				}
			}
		}

		// [discretisation] -- the degree here, the stabilisation after
		// [source], since a box support changes its default.
		{
			Table discretisation( document, "discretisation", sourceName, true );
			discretisation.rejectUnknownKeys( { "PolynomialDegree", "Tau", "TauScaling", "TauFloorRadius" } );

			discretisationOptions.polynomialDegree = discretisation.getInteger( "PolynomialDegree" );
			if ( discretisationOptions.polynomialDegree < 0 )
				discretisation.fail( "PolynomialDegree", "must not be negative" );
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
					                            "SafetyFactorFile", "SafetyFactorDegree",
					                            "ToroidalFieldGuess",
					                            "Normalised", "PsiAxis", "ConfineToPlasma", "PlasmaCurrent",
					                            "ExcludeAttributes", "axis", "SupportShape", "AxialEdge",
					                            "ProfileFile" } );
					refuseReservedProfileFile( source );
					// The "mhd" source has no constant form for either profile,
					// so it reads both files directly and does not go through
					// readEitherProfile. These two calls are what keeps the
					// reserved keys refused here as well; without them all four
					// were accepted and ignored, which is the one failure mode
					// this file exists to prevent.
					refuseReservedVariableKeys( source, "PPrime", false );
					refuseReservedVariableKeys( source, "GGPrime", false );
					MHDParameters parameters;
					parameters.pPrimeFile = source.getString( "PPrimeFile" );
					parameters.mu0 = source.getFloatOr( "Mu0", parameters.mu0 );
					parameters.pPrimeScale = source.getFloatOr( "PPrimeScale", 1.0 );
					parameters.ggPrimeScale = source.getFloatOr( "GGPrimeScale", 1.0 );
					parameters.normalised = source.getBooleanOr( "Normalised", false );
					readToroidalFieldTarget( source, parameters );
					if ( !std::isfinite( parameters.pPrimeScale ) )
						source.fail( "PPrimeScale", "must be finite" );
					if ( !std::isfinite( parameters.ggPrimeScale ) )
						source.fail( "GGPrimeScale", "must be finite" );
					if ( parameters.pPrimeFile.empty() )
						source.fail( "PPrimeFile", "must not be empty" );
					if ( !( parameters.mu0 > 0.0 ) )
						source.fail( "Mu0", "must be positive" );
					readNormalisation( source, parameters.normalised, parameters.psiAxis );
					readPlasmaSupport( source, parameters.normalised, parameters.confineToPlasma );
					readPlasmaExclusion( source, parameters.normalised, sourceOptions.excludeAttributes );
					readPlasmaCurrent( source, parameters.normalised, parameters.plasmaCurrent );
					readAxisPoint( source, sourceName, parameters.normalised, sourceOptions.axisPoint );
					readSupportShape( source, parameters.normalised, parameters.confineToPlasma,
					                  !parameters.safetyFactorFile.empty(),
					                  parameters.plasmaCurrent, sourceOptions );
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
					                            "ConfineToPlasma", "ExcludeAttributes", "PlasmaCurrent", "ProfileFile",
					                            "axis", "SupportShape", "AxialEdge", "Gauge" } );
					refuseReservedProfileFile( source );
					RotatingParameters parameters;

					std::string const gauge = source.getStringOr( "Gauge", "fixed-radius" );
					if ( gauge == "fixed-radius" )
						parameters.gauge = RotatingGaugeType::FixedRadius;
					else if ( gauge == "midplane" )
						parameters.gauge = RotatingGaugeType::Midplane;
					else if ( gauge == "vacuum-midplane" )
						parameters.gauge = RotatingGaugeType::VacuumMidplane;
					else
						source.fail( "Gauge", "must be \"fixed-radius\", \"midplane\" or \"vacuum-midplane\", but is \"" + gauge + "\"" );
					bool const vacuum = parameters.gauge == RotatingGaugeType::VacuumMidplane;
					// Both crossing gauges, which share every refusal below but
					// one: neither has a radius, and both are a mirror's.
					bool const midplane = parameters.gauge == RotatingGaugeType::Midplane || vacuum;

					// REFUSED RATHER THAN IGNORED: the midplane gauges pin phi_0
					// on a crossing of z = 0 and read no radius, so a
					// ReferenceRadius beside one would parse and do nothing.
					if ( midplane && source.has( "ReferenceRadius" ) )
						source.fail( "ReferenceRadius", "is the fixed-radius gauge's constant R_ref, and Gauge = \"" + gauge + "\" pins phi_0 where flux surfaces cross z = 0 instead, so there is no radius to give. Remove one of the two" );
					parameters.referenceRadius = midplane ? 0.0
						: source.getFloatOr( "ReferenceRadius", parameters.referenceRadius );
					parameters.mu0 = source.getFloatOr( "Mu0", parameters.mu0 );
					parameters.normalised = source.getBooleanOr( "Normalised", false );

					if ( !midplane && !( parameters.referenceRadius > 0.0 ) )
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
					readPlasmaSupport( source, parameters.normalised, parameters.confineToPlasma );
					readPlasmaExclusion( source, parameters.normalised, sourceOptions.excludeAttributes );
					readPlasmaCurrent( source, parameters.normalised, parameters.plasmaCurrent );
					readAxisPoint( source, sourceName, parameters.normalised, sourceOptions.axisPoint );
					readSupportShape( source, parameters.normalised, parameters.confineToPlasma,
					                  false, parameters.plasmaCurrent, sourceOptions );

					// THE MIDPLANE GAUGE IS A MIRROR'S, and both refusals say why.
					// Its crossings are re-tabulated once per Newton iteration of
					// the BORDERED solve, which a source without Normalised does
					// not take; and it needs psi monotone along z = 0, which a box
					// between two midplane labels has and a tokamak -- whose
					// midplane crosses every surface twice -- does not.
					if ( midplane && !parameters.normalised )
						source.fail( "Gauge", "\"" + gauge + "\" needs Normalised = true: the solver tabulates the surfaces' midplane crossings in Psi, at the labels of the bordered Newton, which is the normalised path" );
					if ( midplane && sourceOptions.supportShape != PlasmaSupportShapeType::Box )
						source.fail( "Gauge", "\"" + gauge + "\" is a mirror's gauge and needs SupportShape = \"box\": it pins phi_0 where each flux surface crosses z = 0, which is one condition per surface only where psi is monotone along the midplane. A tokamak's midplane crosses every surface twice; use the fixed-radius gauge" );

					// THE VACUUM IS THE SUBTRACTED CONDUCTORS' OWN FIELD, psi_c,
					// and only under the split is it known before the solve:
					// meshed coils make the vacuum a solve of its own, and a
					// curve read off anything else -- the starting guess, a
					// previous answer -- would make what each Density means
					// depend on it. MIRROR-PLAN.md section 5.
					if ( vacuum && !conductorOptions.subtracts() )
						source.fail( "Gauge", "\"vacuum-midplane\" states each Density where the VACUUM's flux surface crosses z = 0, and reads the vacuum as the conductors' own field, which is known before the solve only when they are subtracted; set [conductors] Model = \"filament\" or \"subtracted\", or use Gauge = \"midplane\"" );

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
					parameters.radius0 = source.getFloat( "R0" );
					parameters.kr = source.getFloat( "Kr" );
					parameters.kz = source.getFloat( "Kz" );
					sourceOptions.parameters = parameters;
					break;
				}
			}

			/*
			 * THE BOX'S CUTOFF AND THE MESH, WHERE THE MESH IS THIS FILE'S OWN
			 * BUILD PRODUCT. The source JUMPS at |z| = AxialEdge, and a jump
			 * across a line that cuts elements costs the order silently --
			 * theAlignedCutoffKeepsTheOrder, M-262. Where the file also says how
			 * the mesh is made, MEQ can see whether the line will be a set of
			 * element faces, so it is checked here rather than measured later;
			 * SurfaceAttribute against LimiterR is the precedent. A mesh made
			 * somewhere else is left alone, since this parser cannot read it.
			 */
			if ( sourceOptions.supportIsBox() && meshOptions.generate.given )
			{
				MeshGeneratorConfig const & g = meshOptions.generate;
				double const edge = sourceOptions.axialEdge;
				if ( !g.alignPlasma )
					throw ConfigError( sourceName, "mesh.generate.AlignPlasma", "SupportShape = \"box\" cuts the source off at |z| = AxialEdge, where it jumps, and [mesh.generate] would cut elements across that line, which costs the order silently. Set [mesh.generate] AlignPlasma = true with PlasmaZMin = -AxialEdge and PlasmaZMax = AxialEdge" );
				if ( g.plasmaZMin != -edge )
					throw ConfigError( sourceName, "mesh.generate.PlasmaZMin", "is " + std::to_string( g.plasmaZMin ) + ", and the box support's cutoff is at z = -AxialEdge = " + std::to_string( -edge ) + ": the aligned box has to BE the cutoff for its faces to carry the jump" );
				if ( g.plasmaZMax != edge )
					throw ConfigError( sourceName, "mesh.generate.PlasmaZMax", "is " + std::to_string( g.plasmaZMax ) + ", and the box support's cutoff is at z = AxialEdge = " + std::to_string( edge ) + ": the aligned box has to BE the cutoff for its faces to carry the jump" );
			}
		}

		/*
		 * [discretisation]'s STABILISATION, READ AFTER [source] BECAUSE A BOX
		 * SUPPORT CHANGES ITS DEFAULT. A mirror's plasma sits where nu = 1/R is
		 * 3 to 10 and a constant tau = 1 folds inside one element there, so a
		 * SupportShape = "box" run defaults to MirrorStabilisation, the
		 * diffusion-scaled tau the mirror suites were swept to choose
		 * (MEASUREMENTS.md M-MIRRORTAU); every other run keeps constant
		 * tau = 1, which its fixtures were measured at. Beside the box, Tau
		 * without TauScaling is AMBIGUOUS -- tau0 in m^-1 by the box's
		 * default, a constant tau in m^-2 by every other run's, a factor 1/R
		 * between them -- so it is refused rather than read either way.
		 */
		{
			Table discretisation( document, "discretisation", sourceName, true );
			bool const box = sourceOptions.supportIsBox();
			bool const scalingGiven = discretisation.has( "TauScaling" );
			std::string const scaling = discretisation.getStringOr( "TauScaling",
				box ? "diffusion" : "constant" );
			if ( scaling == "constant" )
				discretisationOptions.tauScaling = TauScalingType::Constant;
			else if ( scaling == "diffusion" )
				discretisationOptions.tauScaling = TauScalingType::Diffusion;
			else
				discretisation.fail( "TauScaling", "must be one of constant, diffusion, but is \""
				                     + scaling + "\"" );
			bool const diffusion = discretisationOptions.tauScaling == TauScalingType::Diffusion;

			if ( box && !scalingGiven && discretisation.has( "Tau" ) )
			{
				std::ostringstream mirror;
				mirror << "tau = " << MirrorStabilisation::tau0 << "/max( R, "
				       << MirrorStabilisation::floorRadius << " m )";
				discretisation.fail( "Tau", "is ambiguous beside [source] SupportShape = \"box\" without TauScaling: "
				                     "a box support's default stabilisation reads it as tau0 in "
				                     "tau = Tau/max( R, R_f ), in m^-1, and every other run as a constant tau in "
				                     "m^-2. Say which with TauScaling = \"diffusion\" or \"constant\"; with "
				                     "neither key the mirror's default is " + mirror.str() );
			}
			discretisationOptions.tau = discretisation.getFloatOr( "Tau",
				box && diffusion ? MirrorStabilisation::tau0 : 1.0 );

			if ( discretisation.has( "TauFloorRadius" ) )
			{
				if ( !diffusion )
					discretisation.fail( "TauFloorRadius", "means nothing under a constant tau: it is R_f in "
					                     "tau = Tau/max( R, R_f ), which TauScaling = \"diffusion\" -- or a box "
					                     "support's default -- applies, and a constant tau never reads R" );
				discretisationOptions.tauFloorRadius = discretisation.getFloat( "TauFloorRadius" );
				if ( !( std::isfinite( discretisationOptions.tauFloorRadius )
				        && discretisationOptions.tauFloorRadius > 0.0 ) )
					discretisation.fail( "TauFloorRadius", "must be a positive length in m; R_f = 0 is the unbounded "
					                     "tau = Tau/R, which costs the flux half an order on a mesh reaching the axis" );
			}

			if ( !( discretisationOptions.tau > 0.0 ) )
				discretisation.fail( "Tau", "must be positive; the HDG stabilisation is tau = O(1) and defaults to 1.0" );
		}

		// [boundary]
		{
			Table boundary( document, "boundary", sourceName, false );
			boundary.rejectUnknownKeys( { "Type", "shape", "limiter", "xpoint", "exterior" } );

			boundaryOptions.type = toBoundaryDataType( boundary, boundary.getStringOr( "Type", "zero" ) );

			if ( boundaryOptions.type == BoundaryDataType::Exact && sourceOptions.type == SourceType::MHD )
				boundary.fail( "Type", "\"exact\" needs a source with a known exact solution; the \"mhd\" source has none" );
			// THE VACUUM GAUGE READS THE VACUUM AS psi_c WITH A ZERO REMAINDER,
			// which non-zero data on Gamma would make false.
			if ( boundaryOptions.type != BoundaryDataType::Zero && sourceOptions.type == SourceType::Rotating
			     && sourceOptions.getRotating().gauge == RotatingGaugeType::VacuumMidplane )
				boundary.fail( "Type", "is not \"zero\", and [source] Gauge = \"vacuum-midplane\" reads the vacuum as the subtracted conductors' field with a ZERO remainder, which data on Gamma would contradict" );

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

			// [boundary.limiter] -- what pins psi_bnd, FB-3: a prescribed
			// contact, or the meshed limiter surface to find one on.
			{
				Table limiter( boundary, "limiter", sourceName, false );
				limiter.rejectUnknownKeys( { "R", "Z", "SurfaceAttribute" } );

				LimiterConfig & l = boundaryOptions.limiter;
				bool const pointGiven = limiter.has( "R" ) || limiter.has( "Z" );
				bool const surfaceGiven = limiter.has( "SurfaceAttribute" );

				// THE POINT AND THE CURVE ARE ALTERNATIVES. Both converge, and
				// to equilibria that differ by the O( h ) a prescribed contact
				// costs -- so a precedence rule here would decide which answer
				// the run reports on the strength of key order. Refuse instead.
				if ( pointGiven && surfaceGiven )
					limiter.fail( "SurfaceAttribute", "a limiter is a point OR a curve, not both: [boundary.limiter] SurfaceAttribute finds the contact on the meshed limiter surface, while R and Z prescribe it. Remove one" );

				l.given = pointGiven || surfaceGiven;
				if ( surfaceGiven )
				{
					l.surfaceAttribute = limiter.getInteger( "SurfaceAttribute" );
					if ( l.surfaceAttribute <= 0 )
						limiter.fail( "SurfaceAttribute", "the limiter region's element attribute must be positive; MFEM numbers attributes from 1, and tools/mesh/halfdisc.py --limiter writes the enclosed region as 20" );

					// WHERE THE MESH IS THIS FILE'S OWN BUILD PRODUCT, THE
					// ATTRIBUTE IS KNOWN AND CAN BE CHECKED. Without this the
					// mistake is found at the solve, where the symptom is
					// "psi_bnd = max psi_h over the empty set" -- an attribute
					// no element carries -- rather than a key that is wrong.
					if ( meshOptions.generate.given )
					{
						if ( !meshOptions.generate.limiterGiven )
							limiter.fail( "SurfaceAttribute", "[mesh.generate] makes this mesh and was not asked for a limiter, so no element will carry this attribute; give [mesh.generate] LimiterR, LimiterZ and LimiterRadius, or pin the contact with [boundary.limiter] R and Z" );
						if ( l.surfaceAttribute != generatedLimiterAttribute )
							limiter.fail( "SurfaceAttribute", "[mesh.generate] writes the limiter's interior as attribute "
							              + std::to_string( generatedLimiterAttribute ) + ", so that is what this must be" );
					}
				}
				if ( pointGiven )
				{
					// BOTH OR NEITHER. A limiter given only its height sits at
					// R = 0, which is the axis; honouring that would pin psi_bnd
					// on the one part of the boundary that is not a limiter.
					l.radius = limiter.getFloat( "R" );
					l.z = limiter.getFloat( "Z" );

					if ( !( l.radius > 0.0 ) )
						limiter.fail( "R", "the limiter contact must be at R > 0; the axis is not a limiter, and the operator's 1/R is not integrable there" );
				}
				// psi_bnd is read ONLY through Psi, so without a normalisation
				// this is an unknown nothing consumes -- a file that would
				// converge, at full order, to an equilibrium it did not
				// describe. Reported against whichever key the file actually
				// used: naming "R" on a file that wrote SurfaceAttribute would
				// point at a key that is not there.
				if ( l.given && !sourceOptions.isNormalised() )
					limiter.fail( surfaceGiven ? "SurfaceAttribute" : "R",
					              "a limiter pins psi_bnd, which only enters through the normalised flux; set [source] Normalised = true or remove [boundary.limiter]" );

				// A MIRROR'S OUTER LABEL IS THIS POINT. The box is posed between
				// two labels read at prescribed points, and a contact FOUND on a
				// limiter surface is a maximum over a curve -- a tokamak's
				// question, with no meaning on a negative span.
				if ( sourceOptions.supportIsBox() )
				{
					if ( surfaceGiven )
						limiter.fail( "SurfaceAttribute", "[source] SupportShape = \"box\" reads the outer label at a PRESCRIBED point; give [boundary.limiter] R and Z instead of a surface" );
					if ( !pointGiven )
						limiter.fail( "R", "is required with [source] SupportShape = \"box\": the box lies between two labels at prescribed points, [source.axis] the inner and [boundary.limiter] R, Z the outer" );
				}

				// THE TWO LABELS AT ONE POINT make the span identically zero, and
				// Psi = ( psi - psi_bnd )/( psi_ax - psi_bnd ) is then undefined
				// everywhere.
				if ( pointGiven && sourceOptions.axisPoint.given
				     && l.radius == sourceOptions.axisPoint.radius
				     && l.z == sourceOptions.axisPoint.z )
					limiter.fail( "R", "is the same point as [source.axis], so psi_ax and psi_bnd are the same flux and the normalised flux is undefined; the two labels must be different points" );
			}

			// [boundary.xpoint] -- psi_bnd at an X-POINT THAT IS ITSELF TWO
			// UNKNOWNS, XP-3. The block above prescribes the point; this one
			// seeds it and lets the Newton move it.
			{
				Table xpoint( boundary, "xpoint", sourceName, false );
				xpoint.rejectUnknownKeys( { "R", "Z" } );

				XPointConfig & x = boundaryOptions.xpoint;
				x.given = xpoint.has( "R" ) || xpoint.has( "Z" );
				if ( x.given )
				{
					// BOTH OR NEITHER, for [boundary.limiter]'s reason: a seed
					// given only its height starts the search on the symmetry
					// axis, where q has near-zeros that sweep() reports as
					// saddles and that no divertor put there.
					x.radius = xpoint.getFloat( "R" );
					x.z = xpoint.getFloat( "Z" );

					if ( !( x.radius > 0.0 ) )
						xpoint.fail( "R", "the X-point seed must be at R > 0; the symmetry axis carries near-zeros of q that are not X-points, and the flux mass ( R q, v ) degenerates there" );

					// BEFORE the limiter pairing below, so that a mirror file
					// is told about its box rather than about the limiter it
					// was required to give.
					if ( sourceOptions.supportIsBox() )
						xpoint.fail( "R", "[source] SupportShape = \"box\" reads psi_bnd at the prescribed point [boundary.limiter] R, Z, and a mirror's outer label is not a separatrix; remove [boundary.xpoint]" );

					// ALL THREE ROUTES PIN ONE UNKNOWN. GradShafranovSolver
					// refuses the combination too, but arriving there would
					// name a method rather than the two blocks a file wrote.
					if ( boundaryOptions.limiter.given )
						xpoint.fail( "R", "[boundary.limiter] and [boundary.xpoint] both pin psi_bnd and are alternatives: the limiter PRESCRIBES the bounding point, which is right for a material limiter, and the X-point is an UNKNOWN the Newton moves. Remove one" );

					// psi_bnd is read ONLY through Psi, exactly as for a
					// limiter: without a normalisation this solves two extra
					// unknowns and a constraint nothing consumes.
					if ( !sourceOptions.isNormalised() )
						xpoint.fail( "R", "an X-point pins psi_bnd, which only enters through the normalised flux; set [source] Normalised = true or remove [boundary.xpoint]" );
				}
			}

			// [boundary.exterior] -- the exact exterior DtN, FB-5, and what
			// makes a run FREE boundary.
			{
				Table exterior( boundary, "exterior", sourceName, false );
				exterior.rejectUnknownKeys( { "Radius", "CentreZ", "Modes" } );

				ExteriorConfig & e = boundaryOptions.exterior;
				e.given = exterior.has( "Radius" ) || exterior.has( "Modes" )
				          || exterior.has( "CentreZ" );
				if ( e.given )
				{
					e.radius = exterior.getFloat( "Radius" );
					e.centreZ = exterior.getFloatOr( "CentreZ", 0.0 );
					e.modes = exterior.getIntegerOr( "Modes", 4 );

					if ( !( e.radius > 0.0 ) )
						exterior.fail( "Radius", "the artificial boundary's radius must be positive" );
					if ( e.modes < 1 )
						exterior.fail( "Modes", "at least one Gegenbauer mode is needed; the exterior map is diagonal, so this is a truncation and not a discretisation" );

					// ALTERNATIVES, NOT LAYERS. This block defines Gamma as a
					// semicircle centred on the axis; [boundary.shape] defines
					// it as a closed MXH surface that may not reach R = 0. A
					// file naming both has described two different curves and
					// silently taking one is how a run ends up solving on a
					// domain nobody asked for.
					if ( boundaryOptions.shape.type != ShapeType::None )
						exterior.fail( "Radius", "[boundary.exterior] defines Gamma as a semicircle about the axis and [boundary.shape] defines it as a closed surface; they are alternatives, so remove one" );

					// The exterior is a VACUUM: psi -> 0 at infinity is what the
					// decaying modes represent. A datum of "exact" on Gamma
					// would impose a second, contradictory condition there.
					if ( boundaryOptions.type != BoundaryDataType::Zero )
						exterior.fail( "Radius", "[boundary] Type must be \"zero\" with an exterior coupling: Gamma carries the transmission condition, not a prescribed datum" );

					// TWO SEMICIRCLES, AND GAMMA IS THE INNER ONE. It is
					// tempting to read [mesh.generate] Radius and this one as
					// the same number and they are not: the generated arc is
					// the BACKGROUND mesh's outer edge, and D_h is cut FROM
					// that mesh as the elements inside Gamma, with the band
					// between the resulting staircase and Gamma bridged by the
					// Cockburn-Solano transfer. So Gamma has to fit STRICTLY
					// inside the disc, which is what the driver checks against
					// the loaded mesh's bounding box -- and what can be
					// checked here, before gmsh has run at all, whenever the
					// disc is this file's own build product.
					if ( meshOptions.generate.given
					     && std::abs( e.centreZ ) + e.radius >= meshOptions.generate.radius )
						exterior.fail( "Radius", "Gamma must fit strictly inside the generated disc, and [mesh.generate] Radius = "
						               + std::to_string( meshOptions.generate.radius )
						               + " does not leave room for it: D_h is CUT FROM that mesh, so the arc gmsh draws is the background's outer edge and not Gamma. Give the generator the larger radius" );
				}
			}
		}

		// [inverse] -- the targets that determine the Free currents.
		//
		// READ AFTER [[coils]] AND [boundary], because both refusals that
		// matter are pairings with them: a Free current and a target come
		// together or not at all, and Track = true pins the X-point that
		// [boundary.xpoint] makes an unknown.
		{
			Table inverse( document, "inverse", sourceName, false );
			inverse.rejectUnknownKeys( { "Gamma", "isoflux", "xpoint" } );
			InverseConfig & v = inverseOptions;
			v.given = inverse.isPresent() && !inverse.isEmpty();

			v.gamma = inverse.getFloatOr( "Gamma", 0.0 );
			if ( !( v.gamma >= 0.0 ) || !std::isfinite( v.gamma ) )
				inverse.fail( "Gamma", "must be finite and non-negative: it is the weight on the step in each free current, in flux per ampere, and zero is no regularisation" );

			for ( Table const & one : inverse.getTableArrayOr( "isoflux" ) )
			{
				one.rejectUnknownKeys( { "R1", "Z1", "R2", "Z2" } );
				IsofluxTargetConfig t;
				t.r1 = one.getFloat( "R1" );
				t.z1 = one.getFloat( "Z1" );
				t.r2 = one.getFloat( "R2" );
				t.z2 = one.getFloat( "Z2" );
				if ( !std::isfinite( t.r1 ) || !std::isfinite( t.z1 )
				     || !std::isfinite( t.r2 ) || !std::isfinite( t.z2 ) )
					one.fail( "R1", "every coordinate must be finite" );
				if ( t.r1 == t.r2 && t.z1 == t.z2 )
					one.fail( "R2", "the two points coincide, so the row is identically zero and constrains nothing" );
				v.isoflux.push_back( t );
			}

			int tracked = 0;
			for ( Table const & one : inverse.getTableArrayOr( "xpoint" ) )
			{
				one.rejectUnknownKeys( { "R", "Z", "Track" } );
				XPointTargetConfig t;
				t.radius = one.getFloat( "R" );
				t.z = one.getFloat( "Z" );
				t.track = one.getBooleanOr( "Track", true );
				if ( !( t.radius > 0.0 ) || !std::isfinite( t.radius ) || !std::isfinite( t.z ) )
					one.fail( "R", "the target must be finite and at R > 0: the symmetry axis carries near-zeros of q that are not X-points" );
				if ( t.track && !boundaryOptions.xpoint.given )
					one.fail( "Track", "Track = true pins the X-point the solve TRACKS, and there is none without [boundary.xpoint]; give one, or Track = false to ask for a critical point at this place without tracking it" );
				if ( t.track && ++tracked > 1 )
					one.fail( "Track", "a second tracked X-point target pins the same two unknowns somewhere else; a solve tracks one X-point" );
				v.xpoint.push_back( t );
			}

			if ( v.given && v.rowCount() == 0 )
				inverse.fail( "isoflux", "[inverse] names no target: give [[inverse.isoflux]] or [[inverse.xpoint]] blocks, which are what determine the Free currents" );
			if ( v.given && !sourceOptions.isNormalised() )
				inverse.fail( v.isoflux.empty() ? "xpoint" : "isoflux", "an inverse solve's currents are border unknowns of the bordered Newton, which runs only on a normalised source; set [source] Normalised = true" );
			if ( v.given && sourceOptions.supportIsBox() )
				inverse.fail( v.isoflux.empty() ? "xpoint" : "isoflux", "an inverse solve has not been measured on a [source] SupportShape = \"box\" plasma -- a mirror -- and its targets are posed on a tokamak's separatrix; remove [inverse]" );
			if ( v.given && conductorOptions.subtracts() )
				inverse.fail( v.isoflux.empty() ? "xpoint" : "isoflux", "free currents are solved for MESHED conductors only: under [conductors] Model = \"subtracted\" or \"filament\" a current also enters the analytic conductor field, which the inverse border's column does not measure" );
			if ( v.given && !coilOptions.anyFree() )
				inverse.fail( v.isoflux.empty() ? "xpoint" : "isoflux", "there are targets and no Free current to meet them: mark a [[coils]] or [[circuits]] block Free = true, or remove [inverse]" );
			if ( !v.given && coilOptions.anyFree() )
			{
				std::string where;
				for ( std::size_t i = 0; i < coilOptions.coils.size() && where.empty(); ++i )
					if ( coilOptions.coils[ i ].free )
						where = "coils[" + std::to_string( i ) + "].Free";
				for ( std::size_t c = 0; c < coilOptions.circuits.size() && where.empty(); ++c )
					if ( coilOptions.circuits[ c ].free )
						where = "circuits[" + std::to_string( c ) + "].Free";
				throw ConfigError( sourceName, where, "makes a current an unknown and nothing determines it: an inverse solve needs targets, so add an [inverse] table with [[inverse.isoflux]] or [[inverse.xpoint]] blocks, or remove Free to solve with the given Current as data" );
			}
		}

		// [solver]
		{
			Table solver( document, "solver", sourceName, false );
			solver.rejectUnknownKeys( { "PicardSweeps", "PicardBlend", "BorderRegularisation",
			                            "BorderCollinearityRegularisation", "TopologyRetry",
			                            "UpDownSymmetry",
			                            "NewtonMaxIterations", "NewtonRelativeTolerance", "NewtonAbsoluteTolerance",
			                            "PlasmaSupportSweeps", "XPointMeritWeight", "LineSearchMerit",
			                            			                            "AssemblyMode",
			                            "LocalFactorMode", "TraceAssemblyMode", "TraceSolver",
			                            "LinearMaxIterations", "LinearTolerance" } );

			solverOptions.newtonMaxIterations = solver.getIntegerOr( "NewtonMaxIterations", solverOptions.newtonMaxIterations );
			solverOptions.newtonRelativeTolerance = solver.getFloatOr( "NewtonRelativeTolerance", solverOptions.newtonRelativeTolerance );
			solverOptions.newtonAbsoluteTolerance = solver.getFloatOr( "NewtonAbsoluteTolerance", solverOptions.newtonAbsoluteTolerance );
			refuseIterativeSolverKeys( solver );

			// PlasmaSupportSweeps -- the support's own outer loop. See
			// SolverConfig for what it buys and M-82 for what its absence
			// costs on a diverted machine.
			solverOptions.plasmaSupportSweeps = solver.getIntegerOr( "PlasmaSupportSweeps", solverOptions.plasmaSupportSweeps );

			/*
			 * `XPointMeritWeight` -- HOW HEAVILY XP-3's TWO ROWS COUNT IN THE
			 * LINE SEARCH, and NOT in the equation. The border still solves
			 * q_r = q_z = 0 either way, so this may change how many iterations
			 * a solve costs and must not change what it converges to.
			 *
			 * Refused without an X-point border, for the reason every other key
			 * here is: a weight on rows that do not exist is a number nothing
			 * will read, which is the accepted-and-ignored failure this schema
			 * refuses everywhere.
			 */
			solverOptions.borderRegularisation = solver.getFloatOr( "BorderRegularisation", solverOptions.borderRegularisation );
			if ( !( solverOptions.borderRegularisation >= 0.0 ) )
				solver.fail( "BorderRegularisation", "must be finite and non-negative; zero is off and leaves the border solve bit-identical" );
			solverOptions.borderCollinearityRegularisation = solver.getFloatOr( "BorderCollinearityRegularisation", solverOptions.borderCollinearityRegularisation );
			if ( !( solverOptions.borderCollinearityRegularisation >= 0.0 ) )
				solver.fail( "BorderCollinearityRegularisation", "must be finite and non-negative; zero leaves the plain Levenberg form" );
			if ( solverOptions.borderCollinearityRegularisation > 0.0
			     && !( solverOptions.borderRegularisation > 0.0 ) )
				solver.fail( "BorderCollinearityRegularisation", "weights a penalty that multiplies BorderRegularisation's lambda, so it does nothing on its own: set BorderRegularisation as well, or remove this key" );
			solverOptions.topologyRetry = solver.getIntegerOr( "TopologyRetry", solverOptions.topologyRetry );
			solverOptions.upDownSymmetry = solver.getBooleanOr( "UpDownSymmetry", solverOptions.upDownSymmetry );
			/*
			 * AND IT IS REFUSED ON A GENERATED MESH THAT IS NOT SYMMETRIC,
			 * BECAUSE THE ALTERNATIVE IS A THROW A LONG WAY IN.
			 *
			 * buildMirrorMaps() runs inside the bordered Newton driver and not
			 * in prepare(), so the refusal arrives after the mesh, the spaces,
			 * the assembly and the guess -- measured, past 150 s on MAST-U at
			 * k = 2. Where the file also says HOW the mesh is made, MEQ knows
			 * at parse time that the answer will be no.
			 *
			 * ONLY WHERE [mesh.generate] IS PRESENT. A mesh somebody else made
			 * may be symmetric and this parser cannot read it, so a bare
			 * [mesh] File is left alone and finds out the old way.
			 */
			if ( solverOptions.upDownSymmetry && meshOptions.generate.given
			     && !meshOptions.generate.symmetric )
				solver.fail( "UpDownSymmetry", "averages every dof with the dof at its own reflection, so it needs a mirror-symmetric MESH -- and gmsh's triangulation of a symmetric geometry is not symmetric, it picks a diagonal and picks freely. This file generates its mesh, so say [mesh.generate] Symmetric = true, which meshes z >= 0 and reflects it" );
			if ( solverOptions.topologyRetry < 0 )
				solver.fail( "TopologyRetry", "cannot be negative; zero is off" );
			solverOptions.xPointMeritWeight = solver.getFloatOr( "XPointMeritWeight", solverOptions.xPointMeritWeight );
			if ( !( solverOptions.xPointMeritWeight > 0.0 ) )
				solver.fail( "XPointMeritWeight", "must be positive: it multiplies the length R h that converts q into a flux, and zero or negative would make the X-point rows count for nothing or against themselves in the merit" );
			/*
			 * `LineSearchMerit` -- WHAT THE ARMIJO BACKTRACKING COMPARES, and
			 * not what is solved or when it stops. `augmented` is the field
			 * residual plus every weighted border constraint, which is what MEQ
			 * has always used; `field` is `|| R ||` alone, which is the only
			 * part of it the Newton step LINEARISES.
			 */
			{
				std::string const merit =
					solver.getStringOr( "LineSearchMerit", "augmented" );
				if ( merit == "augmented" )
					solverOptions.lineSearchMerit = LineSearchMeritChoice::Augmented;
				else if ( merit == "field" )
					solverOptions.lineSearchMerit = LineSearchMeritChoice::Field;
				else
					solver.fail( "LineSearchMerit", "must be \"augmented\" or \"field\"" );
			}

			/*
			 * THE X-POINT BORDER'S OWN CONTROLS, WHICH THE DRIVER SETS ONLY
			 * BESIDE THAT BORDER -- so on a mirror, which has no X-point and
			 * refuses [boundary.xpoint], each would be accepted and ignored.
			 */
			if ( sourceOptions.supportIsBox() )
			{
				if ( solverOptions.upDownSymmetry )
					solver.fail( "UpDownSymmetry", "is applied with the X-point border, and a [source] SupportShape = \"box\" run has none, so the key would be read and do nothing. Remove it" );
				if ( solverOptions.topologyRetry > 0 )
					solver.fail( "TopologyRetry", "is applied with the X-point border, and a [source] SupportShape = \"box\" run has none, so the key would be read and do nothing. Remove it" );
				if ( solverOptions.borderRegularisation > 0.0 )
					solver.fail( "BorderRegularisation", "is applied with the X-point border, and a [source] SupportShape = \"box\" run has none, so the key would be read and do nothing. Remove it" );
			}

			if ( solver.has( "XPointMeritWeight" ) && !boundaryOptions.xpoint.given )
				solver.fail( "XPointMeritWeight", "weights the X-point rows of the bordered Newton, and this file has no [boundary.xpoint] for them to weight" );
			if ( solverOptions.plasmaSupportSweeps < 0 )
				solver.fail( "PlasmaSupportSweeps", "cannot be negative; 0 leaves the support moving inside Newton, which is what a run without this key does" );

			// INERT WITHOUT A CONFINED SOURCE, and refused rather than ignored.
			// The freeze fixes the threshold insidePlasma() tests against and
			// the component the fill reaches -- and NormalisedSource::
			// freezePlasmaEdge is documented inert unless setPlasmaSupport() is
			// on, so on an unconfined source this key buys a loop that re-solves
			// the identical problem until the iteration cap.
			if ( solverOptions.plasmaSupportSweeps > 0 && !sourceOptions.confinesToPlasma() )
				solver.fail( "PlasmaSupportSweeps", "freezes the plasma SUPPORT between solves, and there is no support to freeze unless the source is confined to one: set [source] ConfineToPlasma = true or remove this key" );

			/*
			 * AND NOT UNDER THE BOX, WHICH HAS NO DISCRETE STATE TO SETTLE.
			 *
			 * The loop freezes an ELEMENT-GRANULAR moving support -- XP-1's
			 * connectivity component -- and judges a sweep settled by that
			 * component repeating. The box switches the fill off by a stated
			 * branch, so the component is identically empty and the loop would
			 * "settle" at its second sweep whatever the labels did; and the
			 * labels themselves are doubles, which never repeat exactly, so
			 * judging by them would need a tolerance -- a new solver knob. The
			 * box needs neither: F vanishes at both labels, so a live edge
			 * leaves the Jacobian nothing to miss, and MR-4b converged that way.
			 */
			if ( solverOptions.plasmaSupportSweeps > 0 && sourceOptions.supportIsBox() )
				solver.fail( "PlasmaSupportSweeps", "freezes and re-decides the connectivity fill's element support, and [source] SupportShape = \"box\" has none -- the box is pointwise and its edge stays live within Newton. Remove this key" );

			// PicardSweeps -- the basin, before the border. See SolverConfig.
			solverOptions.picardSweeps = solver.getIntegerOr( "PicardSweeps", solverOptions.picardSweeps );
			solverOptions.picardBlend = solver.getFloatOr( "PicardBlend", solverOptions.picardBlend );

			if ( solverOptions.picardSweeps < 0 )
				solver.fail( "PicardSweeps", "cannot be negative; 0 drives the bordered Newton from the initial guess, which is what a run without this key does" );
			if ( !( solverOptions.picardBlend > 0.0 && solverOptions.picardBlend <= 1.0 ) )
				solver.fail( "PicardBlend", "must lie in ( 0, 1 ]: 1 is the plain fixed point and anything below it under-relaxes. Zero would take none of each sweep's answer, which is not an iteration" );

			// NOTHING TO ITERATE WITHOUT A NORMALISATION, and refused rather
			// than ignored. The whole of what these sweeps do is hold
			// ( psi_ax, psi_bnd ) fixed so the source becomes an ordinary
			// meq::Source, and on an unnormalised source it already is one --
			// so the loop would re-solve the identical problem N times.
			if ( solverOptions.picardSweeps > 0 && !sourceOptions.isNormalised() )
				solver.fail( "PicardSweeps", "iterates the NORMALISATION, and an unnormalised source has none: set [source] Normalised = true or remove this key" );

			/*
			 * NOT BESIDE A PRESCRIBED AXIS, WHICH IS MEASURED RATHER THAN
			 * ARGUED. Each sweep solves the unbordered problem at frozen labels
			 * under Picard-then-Newton, and on both configurations tried with
			 * [source.axis] that solve does not do its job: on the shipped mirror
			 * every sweep's solve diverges -- a non-finite first Newton residual
			 * of 5e+47 -- while the bordered Newton needs no initialiser there,
			 * converging from the cold start in four steps; on
			 * examples/limited-tokamak.toml with its axis pinned the run ends at
			 * exit 2. An initialiser that cannot be shown to help a constraint
			 * is refused beside it until a case measures one.
			 */
			if ( solverOptions.picardSweeps > 0 && sourceOptions.axisPoint.given )
				solver.fail( "PicardSweeps", "is not measured to work beside [source.axis]: on a mirror every sweep's frozen-label solve diverges, and on a tokamak with its axis pinned the run fails. The bordered Newton reads psi_ax at the point from the first step; remove this key" );

			// AssemblyMode and TraceSolver. Both are performance keys and
			// neither may change the answer, which is why they can be exposed at
			// all: the two assembly modes are bit-identical and the three trace
			// solvers agree to 1e-14.
			//
			// The spellings are lower case, as every other Type key in this file
			// is, and the values are compared literally -- there is no
			// case-folding anywhere in Config, so "Threaded" is a fault and says
			// so rather than being quietly accepted.
			//
			// AVAILABILITY IS NOT CHECKED HERE, and that is the MFEM-free rule
			// rather than an omission: whether this build has OpenMP, or
			// PARDISO, or cuDSS is a question about the linked library, and this
			// translation unit is one CI compiles without it. apps/meq.cpp asks
			// GradShafranovSolver::assemblyModeAvailable() and
			// traceSolverAvailable() and refuses there, with a message naming
			// the build option. So a wrong SPELLING fails at parse and an
			// unavailable CHOICE fails at startup -- two different faults with
			// two different messages.
			{
				solverOptions.assemblyModeWasGiven = solver.has( "AssemblyMode" );
				std::string const mode = solver.getStringOr( "AssemblyMode", "threaded" );
				if ( mode == "serial" )
					solverOptions.assemblyMode = AssemblyModeType::Serial;
				else if ( mode == "threaded" )
					solverOptions.assemblyMode = AssemblyModeType::Threaded;
				else if ( mode == "batched" )
					solverOptions.assemblyMode = AssemblyModeType::Batched;
				else
					solver.fail( "AssemblyMode", "must be one of serial, threaded, batched, but is \""
					             + mode + "\"" );

				// THE OTHER TWO BATCHED AXES, AS SEPARATE KEYS. They have
				// different preconditions, different measured host costs, and
				// -- for the trace one -- different bit-exactness, so one key
				// covering all three would make a regression unattributable.
				solverOptions.localFactorModeWasGiven = solver.has( "LocalFactorMode" );
				std::string const local = solver.getStringOr( "LocalFactorMode", "serial" );
				if ( local == "serial" )
					solverOptions.localFactorMode = LocalFactorModeType::Serial;
				else if ( local == "batched" )
					solverOptions.localFactorMode = LocalFactorModeType::Batched;
				else
					solver.fail( "LocalFactorMode", "must be one of serial, batched, but is \""
					             + local + "\"" );

				std::string const traceAsm = solver.getStringOr( "TraceAssemblyMode", "batched" );
				if ( traceAsm == "serial" )
					solverOptions.traceAssemblyMode = TraceAssemblyModeType::Serial;
				else if ( traceAsm == "batched" )
					solverOptions.traceAssemblyMode = TraceAssemblyModeType::Batched;
				else
					solver.fail( "TraceAssemblyMode", "must be one of serial, batched, but is \""
					             + traceAsm + "\"" );

				solverOptions.traceSolverWasGiven = solver.has( "TraceSolver" );
				std::string const trace = solver.getStringOr( "TraceSolver", "pardiso" );
				if ( trace == "umfpack" )
					solverOptions.traceSolver = TraceSolverType::UMFPack;
				else if ( trace == "pardiso" )
					solverOptions.traceSolver = TraceSolverType::Pardiso;
				else if ( trace == "cudss" )
					solverOptions.traceSolver = TraceSolverType::cuDSS;
				else if ( trace == "gmres" )
					solverOptions.traceSolver = TraceSolverType::GMRES;
				else
					solver.fail( "TraceSolver", "must be one of umfpack, pardiso, cudss, gmres, "
					             "but is \"" + trace + "\"" );
			}

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
			output.rejectUnknownKeys( { "Directory", "Prefix", "GridNR", "GridNZ",
			                            "GridRMin", "GridRMax", "GridZMin", "GridZMax",
			                            "FluxSurfaces", "FluxSurfaceCount",
			                            "FluxAngleCount", "FluxInnerCut",
			                            "FluxOuterCut" } );

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

			// THE EXTENT, IN PAIRS. One end without the other would leave the
			// grid half the default box and half the author's, which is a
			// box nobody wrote down.
			auto extent = [ &output ]( char const *lowKey, char const *highKey,
			                           double &low, double &high )
			{
				if ( output.has( lowKey ) != output.has( highKey ) )
					output.fail( output.has( lowKey ) ? highKey : lowKey,
					             std::string( "must be given with " )
					             + ( output.has( lowKey ) ? lowKey : highKey )
					             + ": the grid's extent is set as a pair or not at all" );
				if ( !output.has( lowKey ) )
					return;
				low = output.getFloat( lowKey );
				high = output.getFloat( highKey );
				if ( !std::isfinite( low ) )
					output.fail( lowKey, "must be finite" );
				if ( !std::isfinite( high ) )
					output.fail( highKey, "must be finite" );
				if ( !( high > low ) )
					output.fail( highKey, std::string( "must exceed " ) + lowKey
					             + ": the grid's extent must be positive" );
			};
			extent( "GridRMin", "GridRMax", outputOptions.gridRMin, outputOptions.gridRMax );
			extent( "GridZMin", "GridZMax", outputOptions.gridZMin, outputOptions.gridZMax );
			if ( outputOptions.hasGridRadii() && outputOptions.gridRMin < 0.0 )
				output.fail( "GridRMin", "must not be negative: R is the major radius, "
				             "and the solve is the half-plane R >= 0" );

			// The ( Psi, theta ) file: INVERSION-PLAN.md stage IN-6. See
			// OutputConfig for why it is off by default and why both ends of
			// the cut are cut.
			outputOptions.fluxSurfaces =
				output.getBooleanOr( "FluxSurfaces", outputOptions.fluxSurfaces );
			outputOptions.fluxSurfaceCount =
				output.getIntegerOr( "FluxSurfaceCount", outputOptions.fluxSurfaceCount );
			outputOptions.fluxAngleCount =
				output.getIntegerOr( "FluxAngleCount", outputOptions.fluxAngleCount );
			outputOptions.fluxInnerCut =
				output.getFloatOr( "FluxInnerCut", outputOptions.fluxInnerCut );
			outputOptions.fluxOuterCut =
				output.getFloatOr( "FluxOuterCut", outputOptions.fluxOuterCut );

			// VALIDATED WHETHER OR NOT THE FILE IS BEING WRITTEN, which is the
			// same stance the grid sizes take: a key that is present is a key
			// the author meant, and refusing it at parse costs milliseconds
			// where refusing it after the solve costs the solve.
			if ( outputOptions.fluxSurfaceCount < 2 )
				output.fail( "FluxSurfaceCount", "must be at least 2: a family with one "
				             "surface has nothing to interpolate between" );
			if ( outputOptions.fluxAngleCount < 3 )
				output.fail( "FluxAngleCount", "must be at least 3: fewer nodes than that "
				             "enclose nothing" );
			if ( !( outputOptions.fluxInnerCut > 0.0 ) )
				output.fail( "FluxInnerCut", "must be greater than zero: Psi_N = 0 is the "
				             "magnetic axis, where the surface is a point and the geometry "
				             "derivative is unbounded" );
			if ( !( outputOptions.fluxOuterCut < 1.0 ) )
				output.fail( "FluxOuterCut", "must be less than one: Psi_N = 1 is the "
				             "plasma boundary, where the surface is the boundary itself and "
				             "1/| grad psi | diverges if it is a separatrix" );
			if ( !( outputOptions.fluxInnerCut < outputOptions.fluxOuterCut ) )
				output.fail( "FluxInnerCut", "must be less than FluxOuterCut; the two are "
				             "the ends of the range of normalised flux the family covers" );

			// The ( Psi, theta ) file traces CLOSED surfaces about a magnetic
			// axis, and a mirror's are open field lines through the throat with
			// no axis among them; refused here rather than discovered after
			// the solve.
			if ( outputOptions.fluxSurfaces && sourceOptions.supportIsBox() )
				output.fail( "FluxSurfaces", "traces closed flux surfaces about a magnetic axis, and a [source] SupportShape = \"box\" plasma -- a mirror -- has neither" );
		}

		// [initialguess]
		{
			Table guess( document, "initialguess", sourceName, false );
			guess.rejectUnknownKeys( { "Type", "File", "MeshFile", "Content",
			                           "Amplitude",
			                           "CentreR", "CentreZ", "RadiusR", "RadiusZ" } );

			std::string const type = guess.getStringOr( "Type", "none" );
			if ( type == "none" )
				initialGuessOptions.type = InitialGuessType::None;
			else if ( type == "ramp" )
				initialGuessOptions.type = InitialGuessType::Ramp;
			else if ( type == "bump" )
				initialGuessOptions.type = InitialGuessType::Bump;
			else if ( type == "conductors" )
				initialGuessOptions.type = InitialGuessType::Conductors;
			else if ( type == "gridfunction" )
				initialGuessOptions.type = InitialGuessType::GridFunction;
			else
				guess.fail( "Type", "must be one of none, ramp, bump, conductors, gridfunction, but is \"" + type + "\"" );

			// CONDUCTORS NEEDS CONDUCTORS, and saying so at parse time beats a
			// solve that starts from an identically zero guess and calls it a
			// cold start. An empty [[coils]] with this Type is a file that asked
			// for the vacuum field of nothing.
			if ( initialGuessOptions.type == InitialGuessType::Conductors
			     && coilOptions.coils.empty() )
				guess.fail( "Type", "= \"conductors\" builds the guess from this file's own [[coils]], and there are none" );

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

			// WHAT THE STORED FILE HOLDS. See Config.hpp for why a .gf cannot
			// say this for itself and why the default is "remainder".
			std::string const content = guess.getStringOr( "Content", "remainder" );
			if ( content == "remainder" )
				initialGuessOptions.content = GuessContent::Remainder;
			else if ( content == "total" )
				initialGuessOptions.content = GuessContent::Total;
			else
				guess.fail( "Content", "must be \"total\" or \"remainder\", but is \"" + content
				            + "\". A .gf carries the coefficients and the space and nothing else, so "
				              "under [conductors] Model it cannot say whether it holds psi or "
				              "psi - psi_c; this is where the file's author says. \"remainder\" is "
				              "what a subtracted run writes and what every file written before this "
				              "key existed means" );

			if ( guess.has( "Content" )
			     && initialGuessOptions.type != InitialGuessType::GridFunction )
				guess.fail( "Content", "describes a STORED guess and there is none: it is read "
				            "only with Type = \"gridfunction\". Every other Type builds the guess "
				            "from this file, in whichever field the run solves for" );

			if ( initialGuessOptions.type == InitialGuessType::Conductors )
			{
				/*
				 * THE PLASMA COLUMN IS OPTIONAL AND COSTS TWO NUMBERS.
				 *
				 * With no CentreR this is the conductors' own vacuum field and
				 * nothing else, which is a complete and useful guess: it
				 * carries the machine's scale and topology and is
				 * Delta*-harmonic off the conductors. MEASURED on freegsnke's
				 * MAST-U, it converges in 9 iterations -- to an X-point 1.574 m
				 * from its seed, which is M-26's hazard reached by a guess that
				 * describes no plasma.
				 *
				 * CentreR is therefore the guessed MAGNETIC AXIS, and
				 * `[source] PlasmaCurrent` is spread over an ELLIPSE about it.
				 * That is closer to the physics than a paraboloid -- a plasma
				 * IS a current distribution, and a column of one has a
				 * separatrix and an X-point where a paraboloid has no null at
				 * all -- and it asks the user for something they know.
				 *
				 * **AND IT IS AN ELLIPSE RATHER THAN A FILAMENT BECAUSE A
				 * FILAMENT IS SINGULAR AT THE ONE POINT THAT MATTERS.**
				 * meq::filamentPsi() diverges logarithmically AT the filament,
				 * so carrying I_p on one at the guessed axis puts an unbounded
				 * spike exactly where the axis search has to look: on MAST-U it
				 * reads 6.667e-01 at 3 mm from the guess against a reference
				 * psi_axis of 9.187e-02. A rectangle would be bounded, and the
				 * ellipse is preferred over one because a uniform current
				 * density on a rectangle carries a logarithm in its second
				 * derivatives at each of four corners, which is an artefact of
				 * the shape and not of the equilibrium.
				 *
				 * RadiusR and RadiusZ are its semi-axes, and default to half
				 * the guessed major radius, circular. That default cannot
				 * reach the axis whatever CentreR is, which is the refusal it
				 * would otherwise meet; an explicit RadiusR can, and does.
				 * CentreZ defaults to 0 because an up-down symmetric machine
				 * is the common case. Amplitude is NOT read here and is a
				 * bump's business.
				 */
				initialGuessOptions.centreR = guess.getFloatOr( "CentreR", 0.0 );
				initialGuessOptions.centreZ = guess.getFloatOr( "CentreZ", 0.0 );
				if ( initialGuessOptions.centreR < 0.0 )
					guess.fail( "CentreR", "is the guessed magnetic axis and must be at R > 0. "
					            "Leave it out altogether for the conductors' own field with no "
					            "plasma column in it" );

				initialGuessOptions.radiusR =
					guess.getFloatOr( "RadiusR", 0.5*initialGuessOptions.centreR );
				initialGuessOptions.radiusZ =
					guess.getFloatOr( "RadiusZ", initialGuessOptions.radiusR );

				if ( initialGuessOptions.centreR > 0.0 )
				{
					if ( !( initialGuessOptions.radiusR > 0.0 )
					     || !( initialGuessOptions.radiusZ > 0.0 ) )
						guess.fail( "RadiusR", "the semi-axes of the plasma column must both be "
						            "positive" );

					if ( !( initialGuessOptions.centreR
					        - initialGuessOptions.radiusR > 0.0 ) )
						guess.fail( "RadiusR", "the plasma column reaches the axis: "
						            "CentreR - RadiusR must be strictly positive, because psi is "
						            "identically zero on R = 0 and a conductor there would sit in "
						            "its own zero" );
				}
			}

			if ( initialGuessOptions.type == InitialGuessType::Bump )
			{
				initialGuessOptions.centreR = guess.getFloat( "CentreR" );
				initialGuessOptions.centreZ = guess.getFloatOr( "CentreZ", 0.0 );
				initialGuessOptions.radiusR = guess.getFloat( "RadiusR" );
				initialGuessOptions.radiusZ =
					guess.getFloatOr( "RadiusZ", initialGuessOptions.radiusR );

				if ( !( initialGuessOptions.centreR > 0.0 ) )
					guess.fail( "CentreR", "the magnetic axis is at R > 0; a bump centred on "
					            "or across the axis describes no plasma" );
				if ( !( initialGuessOptions.radiusR > 0.0 )
				     || !( initialGuessOptions.radiusZ > 0.0 ) )
					guess.fail( "RadiusR", "must be positive" );
				if ( !( initialGuessOptions.amplitude > 0.0 ) )
					guess.fail( "Amplitude", "must be positive: it is the peak of the bump, and "
					            "a flat one is the trivial branch this guess exists to avoid" );
				// The bump must not straddle the axis, where psi is pinned and the
				// operator's 1/R is not integrable.
				if ( initialGuessOptions.centreR - initialGuessOptions.radiusR <= 0.0 )
					guess.fail( "RadiusR", "the bump reaches R <= 0; CentreR - RadiusR must be "
					            "strictly positive" );
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
