#include "Config.hpp"

#include "RotatingSource.hpp"

#include <algorithm>
#include <cmath>
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
			// be driven to make int F/r vanish, which it does at scale zero.
			if ( source.has( "PlasmaCurrent" ) && current == 0.0 )
				source.fail( "PlasmaCurrent", "is zero, which is not a plasma: the border would drive the profile scale to zero to satisfy it. Remove the key to leave the amplitude fixed" );
			if ( source.has( "PlasmaCurrent" ) && !std::isfinite( current ) )
				source.fail( "PlasmaCurrent", "must be finite" );
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
			std::initializer_list< char const * > const tables = { "mesh", "discretisation", "source", "boundary", "solver", "output", "initialguess", "adaptivity", "coils" };
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
			for ( std::size_t i = 0; i < blocks.size(); ++i )
			{
				Table const & one = blocks[ i ];
				one.rejectUnknownKeys( { "Name", "CentreR", "CentreZ", "HalfWidth",
				                         "HalfHeight", "Current", "CurrentDensity" } );

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
				// 1/r that is not integrable through r = 0, so a conductor
				// reaching the axis is not a modelling choice this code can
				// honour. Caught here as well as there so that the diagnostic
				// names the coil and the key rather than arriving from a
				// constructor three layers down.
				if ( !( coil.centreR - coil.halfWidth > 0.0 ) )
					one.fail( "CentreR", "the coil reaches or crosses the axis: CentreR - HalfWidth = "
					          + std::to_string( coil.centreR - coil.halfWidth )
					          + " and must be strictly positive, because the operator's 1/r is not integrable through r = 0" );

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
				if ( !hasCurrent && !hasDensity )
					one.fail( "Current", "a coil needs a current: give either Current, the TOTAL through the cross-section in amperes, or CurrentDensity, the uniform j_phi in A/m^2" );

				double const area = 4.0*coil.halfWidth*coil.halfHeight;
				if ( hasDensity )
				{
					coil.densityGiven = true;
					coil.current = one.getFloat( "CurrentDensity" )*area;
				}
				else
				{
					coil.current = one.getFloat( "Current" );
				}

				coilOptions.coils.push_back( coil );
			}
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
				                              "Vessel", "Transition", "Check" } );

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
						generate.fail( "Tool", "unknown mesh generator \"" + g.tool + "\"; the generators are [halfdisc], which is tools/mesh/halfdisc.py -- a semicircle reaching r = 0 exactly, with the conductors fragmented in" );

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

					if ( !( g.radius > 0.0 ) )
						generate.fail( "Radius", "the disc's radius must be positive. NOTE it is the BACKGROUND's radius and not Gamma's: with an exterior coupling, D_h is cut from this mesh at [boundary.exterior] Radius, which has to fit strictly inside it" );
					if ( !( g.size > 0.0 ) )
						generate.fail( "Size", "the background element size must be positive" );
					if ( g.order < 1 )
						generate.fail( "Order", "the geometric order must be at least 1; above 1 the arc's mid-edge nodes are placed on the true circle rather than on the chord" );
					if ( g.coilSize < 0.0 )
						generate.fail( "CoilSize", "the element size inside the conductors must not be negative; omit it for the background Size" );
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
							generate.fail( "PlasmaRMin", "must not be negative: r is a cylindrical radius" );
						if ( g.plasmaRMax <= g.plasmaRMin )
							generate.fail( "PlasmaRMax", "must be greater than PlasmaRMin" );
						if ( g.plasmaZMax <= g.plasmaZMin )
							generate.fail( "PlasmaZMax", "must be greater than PlasmaZMin" );
						if ( !( g.plasmaSize > 0.0 ) )
							generate.fail( "PlasmaSize", "the element size in the refined box must be positive" );
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
							               + " and must be strictly positive, because a closed plasma surface through r = 0 carries a non-integrable 1/r" );
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
								generate.fail( "Vessel", "reaches r = "
								               + std::to_string( g.vessel[ k ] )
								               + ", and the domain is r >= 0" );
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
					                            "SafetyFactorFile", "SafetyFactorDegree",
					                            "ToroidalFieldGuess",
					                            "Normalised", "PsiAxis", "ConfineToPlasma", "PlasmaCurrent",
					                            "ExcludeAttributes",
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
					                            "ConfineToPlasma", "ExcludeAttributes", "PlasmaCurrent", "ProfileFile" } );
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
					readPlasmaSupport( source, parameters.normalised, parameters.confineToPlasma );
					readPlasmaExclusion( source, parameters.normalised, sourceOptions.excludeAttributes );
					readPlasmaCurrent( source, parameters.normalised, parameters.plasmaCurrent );

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
			boundary.rejectUnknownKeys( { "Type", "shape", "limiter", "xpoint", "exterior" } );

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
					// r = 0, which is the axis; honouring that would pin psi_bnd
					// on the one part of the boundary that is not a limiter.
					l.r = limiter.getFloat( "R" );
					l.z = limiter.getFloat( "Z" );

					if ( !( l.r > 0.0 ) )
						limiter.fail( "R", "the limiter contact must be at R > 0; the axis is not a limiter, and the operator's 1/r is not integrable there" );
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
					x.r = xpoint.getFloat( "R" );
					x.z = xpoint.getFloat( "Z" );

					if ( !( x.r > 0.0 ) )
						xpoint.fail( "R", "the X-point seed must be at R > 0; the symmetry axis carries near-zeros of q that are not X-points, and the flux mass ( r q, v ) degenerates there" );

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
					// it as a closed MXH surface that may not reach r = 0. A
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

		// [solver]
		{
			Table solver( document, "solver", sourceName, false );
			solver.rejectUnknownKeys( { "PicardSweeps", "PicardBlend",
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
			solverOptions.xPointMeritWeight = solver.getFloatOr( "XPointMeritWeight", solverOptions.xPointMeritWeight );
			if ( !( solverOptions.xPointMeritWeight > 0.0 ) )
				solver.fail( "XPointMeritWeight", "must be positive: it multiplies the length r h that converts q into a flux, and zero or negative would make the X-point rows count for nothing or against themselves in the merit" );
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
				else
					solver.fail( "TraceSolver", "must be one of umfpack, pardiso, cudss, but is \""
					             + trace + "\"" );
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
		}

		// [initialguess]
		{
			Table guess( document, "initialguess", sourceName, false );
			guess.rejectUnknownKeys( { "Type", "File", "MeshFile", "Amplitude",
			                           "CentreR", "CentreZ", "RadiusR", "RadiusZ" } );

			std::string const type = guess.getStringOr( "Type", "none" );
			if ( type == "none" )
				initialGuessOptions.type = InitialGuessType::None;
			else if ( type == "ramp" )
				initialGuessOptions.type = InitialGuessType::Ramp;
			else if ( type == "bump" )
				initialGuessOptions.type = InitialGuessType::Bump;
			else if ( type == "gridfunction" )
				initialGuessOptions.type = InitialGuessType::GridFunction;
			else
				guess.fail( "Type", "must be one of none, ramp, bump, gridfunction, but is \"" + type + "\"" );

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

			if ( initialGuessOptions.type == InitialGuessType::Bump )
			{
				initialGuessOptions.centreR = guess.getFloat( "CentreR" );
				initialGuessOptions.centreZ = guess.getFloatOr( "CentreZ", 0.0 );
				initialGuessOptions.radiusR = guess.getFloat( "RadiusR" );
				initialGuessOptions.radiusZ =
					guess.getFloatOr( "RadiusZ", initialGuessOptions.radiusR );

				if ( !( initialGuessOptions.centreR > 0.0 ) )
					guess.fail( "CentreR", "the magnetic axis is at r > 0; a bump centred on "
					            "or across the axis describes no plasma" );
				if ( !( initialGuessOptions.radiusR > 0.0 )
				     || !( initialGuessOptions.radiusZ > 0.0 ) )
					guess.fail( "RadiusR", "must be positive" );
				if ( !( initialGuessOptions.amplitude > 0.0 ) )
					guess.fail( "Amplitude", "must be positive: it is the peak of the bump, and "
					            "a flat one is the trivial branch this guess exists to avoid" );
				// The bump must not straddle the axis, where psi is pinned and the
				// operator's 1/r is not integrable.
				if ( initialGuessOptions.centreR - initialGuessOptions.radiusR <= 0.0 )
					guess.fail( "RadiusR", "the bump reaches r <= 0; CentreR - RadiusR must be "
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
