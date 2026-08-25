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
}

#endif // MEQ_SOURCEFACTORY_HPP
