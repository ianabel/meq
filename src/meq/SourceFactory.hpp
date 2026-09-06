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

#include "Coils.hpp"
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

	/// Construct the coil set the `[[coils]]` blocks describe, or a null
	/// pointer if there are none.
	///
	/// **A NULL RETURN IS THE ORDINARY CASE AND IS NOT AN ERROR.** Every
	/// fixed-boundary configuration in `examples/` has no coils at all, and a
	/// caller distinguishes "no coils" from "an empty CoilSet" only by which
	/// of the two it gets -- so the caller can skip the augmentation entirely
	/// rather than wrapping its source around a set that adds zero.
	///
	/// The parse has already resolved `Current` / `CurrentDensity` to a total
	/// current and refused a coil that names both or neither. What is left for
	/// meq::Coil to refuse is the geometry: a non-positive half-extent, and a
	/// coil reaching the axis. Those refusals are translated here so that they
	/// read as the configuration errors they are and NAME THE COIL, which the
	/// library exception cannot do -- `meq::Coil` does not know it came from a
	/// file or which block it was.
	///
	/// @throws ConfigError naming the offending `[[coils]]` block by the name
	///         it was given, or by `coils[i]` if it was not named.
	std::shared_ptr<CoilSet const> makeCoilSet( CoilConfig const &config,
	                                            double mu0 = vacuumPermeability,
	                                            std::string const &configFileName = std::string() );
}

#endif // MEQ_SOURCEFACTORY_HPP
