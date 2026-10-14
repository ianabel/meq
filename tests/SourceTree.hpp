#ifndef MEQ_TESTS_SOURCETREE_HPP
#define MEQ_TESTS_SOURCETREE_HPP

/*
 * WHERE A TEST FINDS WHAT IT READS.
 *
 * Every test runs in its build tree (tests/CMakeLists.txt) and writes only
 * there, so nothing it reads from the source tree can be named relative to the
 * working directory. MEQ_SOURCE_DIR and MEQ_RUNS_DIR are defined for every
 * test target by tests/CMakeLists.txt.
 *
 * A SHIPPED CONFIGURATION RUN AS IT IS needs nothing more than sourcePath():
 * meq::Configuration resolves the paths a file reads against the file's own
 * directory. A VARIANT of one -- its text, edited, written into the build tree
 * -- has moved away from its data, and shippedConfig() is what moves the data
 * paths with it.
 */

#include <filesystem>
#include <fstream>
#include <regex>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>

#ifndef MEQ_SOURCE_DIR
#error "MEQ_SOURCE_DIR is defined for every test by tests/CMakeLists.txt"
#endif
#ifndef MEQ_RUNS_DIR
#error "MEQ_RUNS_DIR is defined for every test by tests/CMakeLists.txt"
#endif

namespace meq::test
{
	/// `relative` from the root of the source tree: "examples/mhd-rectangle.toml".
	inline std::string sourcePath( std::string const & relative )
	{
		return ( std::filesystem::path( MEQ_SOURCE_DIR ) / relative ).lexically_normal().string();
	}

	/// `relative` from the runs directory, where tools/ write -- an independent
	/// code's answer that a test compares against.
	inline std::string runsPath( std::string const & relative )
	{
		return ( std::filesystem::path( MEQ_RUNS_DIR ) / relative ).lexically_normal().string();
	}

	/**
	 * The text of `text`, a configuration that lives in `directory`, with every
	 * relative path it READS made absolute against that directory -- exactly
	 * the resolution meq::Configuration applies to a file, done to the text so
	 * that the text can be written anywhere.
	 *
	 * The keys are meq::Configuration::resolveInputPaths()'s. CacheFile and
	 * [output] Directory are written by the run and stay relative to the
	 * working directory, and so does [mesh] File when a [mesh.generate] block
	 * makes it.
	 */
	inline std::string relocateInputs( std::string const & text, std::string const & directory )
	{
		static std::set<std::string> const inputs = {
			"PPrimeFile", "GGPrimeFile", "SafetyFactorFile", "OmegaFile",
			"TemperatureFile", "DensityFile", "MeshFile", "File" };
		static std::regex const header( R"(^\s*\[+\s*([^\]\s]+)\s*\]+)" );
		static std::regex const assignment( R"(^(\s*)([A-Za-z]+)(\s*=\s*")([^"]*)(".*)$)" );

		bool const meshGenerated = text.find( "[mesh.generate]" ) != std::string::npos;

		std::istringstream in( text );
		std::ostringstream out;
		std::string table;
		std::string line;
		while ( std::getline( in, line ) )
		{
			std::smatch match;
			if ( std::regex_search( line, match, header ) )
				table = match[ 1 ];
			else if ( std::regex_match( line, match, assignment )
			          && inputs.count( match[ 2 ] ) != 0
			          && !( table == "mesh" && match[ 2 ] == "File" && meshGenerated )
			          && !match[ 4 ].str().empty()
			          && !std::filesystem::path( match[ 4 ].str() ).is_absolute() )
			{
				line = match[ 1 ].str() + match[ 2 ].str() + match[ 3 ].str()
				     + ( std::filesystem::path( directory ) / match[ 4 ].str() ).lexically_normal().string()
				     + match[ 5 ].str();
			}
			out << line << '\n';
		}
		return out.str();
	}

	/// The text of the shipped configuration at `relative` from the source
	/// root, with its input paths made absolute by relocateInputs(), ready to
	/// be edited and written into the build tree.
	inline std::string shippedConfig( std::string const & relative )
	{
		std::string const path = sourcePath( relative );
		std::ifstream file( path );
		if ( !file )
			throw std::runtime_error( "cannot open " + path );
		std::stringstream buffer;
		buffer << file.rdbuf();
		return relocateInputs( buffer.str(), std::filesystem::path( path ).parent_path().string() );
	}
}

#endif // MEQ_TESTS_SOURCETREE_HPP
