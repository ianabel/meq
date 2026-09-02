#ifndef MEQ_SOURCEFACTORY_HPP
#define MEQ_SOURCEFACTORY_HPP

/*
 * Builds a meq::Source from a parsed meq::SourceConfig.
 *
 * This is deliberately not a method on Configuration. Config.hpp reports what
 * the TOML said and nothing else -- it does not know what a Source is, does not
 * include Source.hpp, and does not read the profile files that an "mhd" source
 * names. Keeping the two apart is what lets the configuration layer be tested
 * without the numerics and the numerics without a config file.
 */

#include <memory>
#include <string>

#include "Config.hpp"
#include "Source.hpp"

namespace meq
{
	/// Construct the source the configuration describes.
	///
	/// For SourceType::MHD this reads the two profile files named in the
	/// configuration, so it touches the filesystem and can fail for reasons the
	/// configuration parse could not have caught -- a missing file, or one whose
	/// contents are not a spline table.
	///
	/// @throws ConfigError if a profile file cannot be read or parsed. The
	///         message names the file and the key it came from, so a failure
	///         here reads like the configuration error it usually is rather
	///         than like an I/O error from somewhere in the numerics.
	/// @param configFileName  where the configuration came from, used only to
	///                        make a thrown ConfigError name its own source.
	std::shared_ptr<Source const> makeSource( SourceConfig const &config,
	                                          std::string const &configFileName = std::string() );

	/// Construct a source whose profiles are functions of NORMALISED flux, so
	/// that psi_ax is an unknown of the non-linear system rather than data.
	///
	/// NON-CONST, AND THAT IS THE POINT. The solver calls setNormalisation() on
	/// such a source before every residual evaluation, so it cannot be held by
	/// const reference and cannot come back from makeSource(). Hand the result
	/// to GradShafranovSolver::setSource( NormalisedSource &, double ) together
	/// with SourceConfig::psiAxisGuess(); the object must outlive the solve.
	///
	/// @throws ConfigError if the configuration is not a normalised one, if a
	///         profile file cannot be read, or if the source itself refuses its
	///         arguments -- the last translated so that a failure reads as the
	///         configuration error it is rather than as a library exception.
	std::shared_ptr<NormalisedSource> makeNormalisedSource( SourceConfig const &config,
	                                                        std::string const &configFileName = std::string() );
}

#endif // MEQ_SOURCEFACTORY_HPP
