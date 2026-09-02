#include "SourceFactory.hpp"

#include "RotatingSource.hpp"

#include <cmath>
#include <stdexcept>
#include <utility>
#include <string>
#include <vector>

namespace meq
{
	namespace
	{
		/*
		 * The manufactured source of refs/HDG-GradShafranov.pdf Example 5,
		 * for the exact solution psi = sin( kr ( r + r0 ) ) cos( kz z ):
		 *
		 *     F( r, z, psi ) = ( kr^2 + kz^2 ) psi
		 *                    + ( kr / r ) cos( kr ( r + r0 ) ) cos( kz z )
		 *                    + r ( sc^2 - psi^2 + exp( -sc ) - exp( -psi ) )
		 *
		 * with sc the exact flux at ( r, z ). Note that F depends on position
		 * both through r and through sc, and that it returns F rather than F/r,
		 * as meq::Source::f() requires.
		 *
		 * It is file-local because nothing outside the factory should reach it:
		 * a manufactured source is a thing to run, not a thing to build against.
		 *
		 * THE SAME FORMULA IS IN tests/analytic/ManufacturedNonlinear.hpp, which
		 * additionally carries the exact solution and its gradient. That is a
		 * duplication, and deliberately not resolved by having the library
		 * depend on a test fixture or the reverse. It is closed instead by
		 * tests/unit/SourceFactoryTests.cpp, which asserts the two agree
		 * pointwise -- so if either drifts, that test fails rather than a
		 * convergence study quietly measuring the wrong problem.
		 */
		class ManufacturedSource : public Source
		{
			public:
				ManufacturedSource( double r0In, double krIn, double kzIn )
					: r0Value( r0In ), krValue( krIn ), kzValue( kzIn )
				{
				}

				double f( double r, double z, double psiValue ) const override
				{
					double const sc = exactPsi( r, z );
					double const linear = ( krValue*krValue + kzValue*kzValue )*psiValue;
					double const geometric = ( krValue/r )
					                         *std::cos( krValue*( r + r0Value ) )
					                         *std::cos( kzValue*z );
					double const nonlinear = sc*sc - psiValue*psiValue
					                         + std::exp( -sc ) - std::exp( -psiValue );
					return linear + geometric + r*nonlinear;
				}

				double dFdPsi( double r, double, double psiValue ) const override
				{
					return krValue*krValue + kzValue*kzValue
					       + r*( -2.0*psiValue + std::exp( -psiValue ) );
				}

			private:
				double exactPsi( double r, double z ) const
				{
					return std::sin( krValue*( r + r0Value ) )*std::cos( kzValue*z );
				}

				double r0Value, krValue, kzValue;
		};

		/// Load one profile, turning any failure into a ConfigError that names
		/// the file and the key it came from.
		std::shared_ptr<Profile const> loadProfile( std::string const &fileName,
		                                            std::string const &key,
		                                            std::string const &configFileName )
		{
			try
			{
				return std::make_shared<SplineProfile const>( SplineProfile::fromFile( fileName ) );
			}
			catch ( std::exception const &error )
			{
				throw ConfigError( configFileName, "source." + key,
				                   "could not read the profile file '" + fileName
				                   + "': " + error.what() );
			}
		}

		/// A profile given EITHER as a constant OR as a table, with a scale
		/// applied to whichever it was. Null when neither was given, which is
		/// how an optional profile -- omega, meaning no rotation -- is spelled.
		///
		/// The scale is folded into a constant rather than wrapped, and a table
		/// at scale 1 is returned unwrapped, so the common cases cost nothing:
		/// this is called once per profile at construction, but the result is
		/// evaluated at every quadrature point of every Newton iteration.
		std::shared_ptr<Profile const> loadEitherProfile( double value, std::string const &fileName,
		                                                  double scale, bool given,
		                                                  std::string const &key,
		                                                  std::string const &configFileName )
		{
			if ( !given )
				return nullptr;

			if ( fileName.empty() )
				return std::make_shared<ConstantProfile const>( value*scale );

			std::shared_ptr<Profile const> table = loadProfile( fileName, key + "File", configFileName );
			if ( scale == 1.0 )
				return table;

			return std::make_shared<ScaledProfile const>( std::move( table ), scale );
		}

		/// The species of a rotating source, with the one marked Neutralising
		/// given the density charge neutrality implies.
		///
		/// Profiles are loaded in index order into named locals, for the reason
		/// the MHD case records below: C++17 does not fix argument evaluation
		/// order, and a diagnostic that blames an arbitrary one of several bad
		/// files sends the reader to fix a file that may be perfectly fine.
		std::vector<Species> buildSpecies( RotatingParameters const &parameters,
		                                   std::string const &configFileName )
		{
			std::vector<Species> species( parameters.species.size() );
			std::size_t neutralising = parameters.species.size();

			for ( std::size_t i = 0; i < parameters.species.size(); ++i )
			{
				SpeciesParameters const &one = parameters.species[ i ];
				std::string const where = "species[" + std::to_string( i ) + "].";

				species[ i ].mass = one.mass;
				species[ i ].charge = one.charge;
				species[ i ].temperature = loadEitherProfile( one.temperature, one.temperatureFile,
					one.temperatureScale, true, where + "Temperature", configFileName );

				if ( one.neutralising )
				{
					neutralising = i;
					continue;
				}

				species[ i ].density = loadEitherProfile( one.density, one.densityFile,
					one.densityScale, true, where + "Density", configFileName );
			}

			// The parser has already refused a set without exactly one of these,
			// so this is a consistency check on that rather than a second policy.
			if ( neutralising == parameters.species.size() )
				throw ConfigError( configFileName, "source.species",
				                   "no species is marked Neutralising, so one density is undetermined" );

			species[ neutralising ].density = neutralisingDensity( species, neutralising );
			return species;
		}

		/// Everything a rotating source needs except its normalisation, which
		/// is what the two factories below differ in.
		struct RotatingPieces
		{
			std::vector<Species> species;
			std::shared_ptr<Profile const> omega;
			std::shared_ptr<Profile const> ggPrime;
		};

		RotatingPieces loadRotating( RotatingParameters const &parameters,
		                             std::string const &configFileName )
		{
			RotatingPieces pieces;
			pieces.species = buildSpecies( parameters, configFileName );
			pieces.omega = loadEitherProfile( parameters.omega, parameters.omegaFile,
				parameters.omegaScale, parameters.omegaGiven, "Omega", configFileName );
			pieces.ggPrime = loadEitherProfile( parameters.ggPrime, parameters.ggPrimeFile,
				parameters.ggPrimeScale, true, "GGPrime", configFileName );
			return pieces;
		}

		/// meq::RotatingSource and friends validate their own arguments and throw
		/// std::invalid_argument. That is right for a library call and wrong for
		/// a configuration file, whose author wants to know which key to edit --
		/// so every construction below goes through this.
		template<typename Build>
		auto guarded( Build build, std::string const &configFileName, char const *key )
		{
			try
			{
				return build();
			}
			catch ( ConfigError const & )
			{
				throw;
			}
			catch ( std::exception const &error )
			{
				throw ConfigError( configFileName, key, error.what() );
			}
		}

	}

	std::shared_ptr<Source const> makeSource( SourceConfig const &config,
	                                          std::string const &configFileName )
	{
		// A NORMALISED SOURCE MUST NOT COME OUT OF HERE, and the reason is that
		// it would work. meq::NormalisedRotatingSource is-a meq::Source, so
		// returning one would compile and solve -- with psi_ax frozen at the
		// guess for ever, because nothing would call setNormalisation() on it.
		// The answer would be a converged solution to a problem the file did not
		// describe. Hence a throw rather than a static_cast.
		if ( config.isNormalised() )
			throw ConfigError( configFileName, "source.Normalised",
			                   "this source is specified in normalised flux, so psi_ax is an "
			                   "unknown of the system; it must be built with "
			                   "meq::makeNormalisedSource and handed to the solver through "
			                   "setSource( NormalisedSource &, double )" );

		switch ( config.type )
		{
			case SourceType::Soloviev:
			{
				return std::make_shared<SolovievSource const>( config.getSoloviev().a );
			}

			case SourceType::MHD:
			{
				MHDParameters const &parameters = config.getMHD();

				// Loaded into named locals, in order, rather than as arguments
				// to make_shared. C++17 does not fix the evaluation order of
				// function arguments, so with both files missing the compiler
				// chose which one the error message blamed -- gcc reported the
				// second. A diagnostic that names an arbitrary one of two bad
				// keys is worse than useless, because it sends the reader to
				// fix a file that may be perfectly fine.
				auto pPrime = loadEitherProfile( 0.0, parameters.pPrimeFile, parameters.pPrimeScale,
				                                 true, "PPrime", configFileName );
				auto ggPrime = loadEitherProfile( 0.0, parameters.ggPrimeFile, parameters.ggPrimeScale,
				                                  true, "GGPrime", configFileName );

				return std::make_shared<MHDSource const>( std::move( pPrime ),
				                                          std::move( ggPrime ),
				                                          parameters.mu0 );
			}

			case SourceType::Manufactured:
			{
				ManufacturedParameters const &parameters = config.getManufactured();
				return std::make_shared<ManufacturedSource const>(
					parameters.r0, parameters.kr, parameters.kz );
			}

			case SourceType::Rotating:
			{
				RotatingParameters const &parameters = config.getRotating();
				RotatingPieces pieces = loadRotating( parameters, configFileName );

				return guarded( [ & ]() -> std::shared_ptr<Source const>
				{
					return std::make_shared<RotatingSource const>( std::move( pieces.species ),
						std::move( pieces.omega ), std::move( pieces.ggPrime ),
						parameters.referenceRadius, parameters.mu0 );
				}, configFileName, "source.species" );
			}
		}

		// Unreachable for a Configuration, which validates the discriminator
		// before it is stored. Present so that adding a SourceType without
		// adding a case here fails loudly rather than returning null.
		throw ConfigError( configFileName, "source.Type",
		                   "no factory case for this source type" );
	}

	std::shared_ptr<NormalisedSource> makeNormalisedSource( SourceConfig const &config,
	                                                        std::string const &configFileName )
	{
		if ( !config.isNormalised() )
			throw ConfigError( configFileName, "source.Normalised",
			                   "this source is not specified in normalised flux, so psi_ax is "
			                   "data rather than an unknown; build it with meq::makeSource" );

		double const psiAxis = config.psiAxisGuess();

		switch ( config.type )
		{
			case SourceType::MHD:
			{
				MHDParameters const &parameters = config.getMHD();

				auto pPrime = loadEitherProfile( 0.0, parameters.pPrimeFile, parameters.pPrimeScale,
				                                 true, "PPrime", configFileName );
				auto ggPrime = loadEitherProfile( 0.0, parameters.ggPrimeFile, parameters.ggPrimeScale,
				                                  true, "GGPrime", configFileName );

				return guarded( [ & ]() -> std::shared_ptr<NormalisedSource>
				{
					return std::make_shared<NormalisedMHDSource>( std::move( pPrime ),
						std::move( ggPrime ), psiAxis, parameters.mu0 );
				}, configFileName, "source.PsiAxis" );
			}

			case SourceType::Rotating:
			{
				RotatingParameters const &parameters = config.getRotating();
				RotatingPieces pieces = loadRotating( parameters, configFileName );

				return guarded( [ & ]() -> std::shared_ptr<NormalisedSource>
				{
					return std::make_shared<NormalisedRotatingSource>( std::move( pieces.species ),
						std::move( pieces.omega ), std::move( pieces.ggPrime ),
						parameters.referenceRadius, psiAxis, parameters.mu0 );
				}, configFileName, "source.species" );
			}

			case SourceType::Soloviev:
			case SourceType::Manufactured:
				break;
		}

		// Unreachable: isNormalised() is false for both of the types above, so
		// the guard at the top has already thrown. Present for the same reason
		// the one in makeSource is.
		throw ConfigError( configFileName, "source.Type",
		                   "this source type has no normalised form" );
	}
}
