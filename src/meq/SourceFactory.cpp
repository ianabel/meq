#include "SourceFactory.hpp"

#include <cmath>
#include <stdexcept>
#include <utility>
#include <string>

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
	}

	std::shared_ptr<Source const> makeSource( SourceConfig const &config,
	                                          std::string const &configFileName )
	{
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
				auto pPrime = loadProfile( parameters.pPrimeFile, "PPrimeFile", configFileName );
				auto ggPrime = loadProfile( parameters.ggPrimeFile, "GGPrimeFile", configFileName );

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
		}

		// Unreachable for a Configuration, which validates the discriminator
		// before it is stored. Present so that adding a SourceType without
		// adding a case here fails loudly rather than returning null.
		throw ConfigError( configFileName, "source.Type",
		                   "no factory case for this source type" );
	}
}
