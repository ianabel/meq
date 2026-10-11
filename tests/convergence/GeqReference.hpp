#ifndef MEQ_TESTS_GEQ_REFERENCE_HPP
#define MEQ_TESTS_GEQ_REFERENCE_HPP

/*
 * ONE READER FOR AN EXPORTED geq CASE: meta.txt, and the species and rotation
 * tables it names. tools/geq-benchmark/README.md documents the format.
 *
 * GeqSourceComparison evaluates the source at a prescribed geq state and
 * RotatingMirror solves geq's free-boundary mirror; both are handed the same
 * kind of case, and a second reader would let the two disagree about what a
 * case says -- which column is the slope, which density file a species names,
 * whether a raw table is taken. So the tables are read here and only here.
 *
 * MFEM-FREE AND BOOST-FREE ON PURPOSE. GeqSourceComparison builds without MFEM,
 * so this reports a missing file by throwing, with the path, and leaves the
 * test framework to its caller.
 */

#include <fstream>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "meq/Profiles.hpp"
#include "meq/RotatingSource.hpp"

namespace meqtest
{
namespace geq
{
	/**
	 * THE COMMAND THAT REGENERATES A CASE'S OUTPUT, by the case's directory
	 * name, as tools/geq-benchmark/README.md gives it.
	 *
	 * WHAT IS TRACKED IS WHAT MEQ IS GIVEN, AND WHAT IS NOT IS geq'S ANSWER.
	 * meta.txt and the tables are committed; the point sets, the psi files and
	 * solution.txt are a geq run's output and are gitignored, so a fresh clone
	 * has none of them and a reader of one must say how to make it rather than
	 * fail to open a file.
	 */
	inline std::string regenerationCommand( std::string directory )
	{
		while ( !directory.empty() && directory.back() == '/' )
			directory.pop_back();
		std::string const name = directory.substr( directory.find_last_of( '/' ) + 1 );
		std::string const run = "PYTHONPATH=../geq/src tools/freegs4e-benchmark/venv/bin/python "
		                        "tools/geq-benchmark/";
		std::string const out = " --out tools/geq-benchmark/reference/" + name;
		if ( name == "two" )
			return run + "export_mirror.py --case two --knots 2049 --omega-axis 4.758e6" + out;
		if ( name == "three" )
			return run + "export_mirror.py --case three --knots 1025" + out;
		if ( name == "vacuum" )
			return run + "export_vacuum.py" + out;
		if ( name == "rotating" )
			return run + "export_rotating.py" + out + " ( about eleven minutes )";
		return "the export script tools/geq-benchmark/README.md names for " + directory;
	}

	/**
	 * Opens one of a case's OUTPUT files -- geq's answer, untracked -- from the
	 * source tree, which is ctest's working directory.
	 *
	 * @throws std::runtime_error naming the file and the command that
	 *         regenerates it, when it cannot be opened.
	 */
	inline std::ifstream openOutput( std::string const &directory, std::string const &name )
	{
		std::string const path = directory
			+ ( !directory.empty() && directory.back() == '/' ? "" : "/" ) + name;
		std::ifstream in( path );
		if ( !in )
			throw std::runtime_error( "cannot open " + path + ". It is geq's ANSWER, an "
				"untracked output of the export rather than part of the committed "
				"reference, so a fresh clone does not have it. Regenerate it, from the "
				"repository root, with: " + regenerationCommand( directory ) );
		return in;
	}

	/// What one case's meta.txt says. A key a case does not write keeps the
	/// value below: a prescribed state has no points, a solved one no labels.
	struct GeqCase
	{
		std::string directory;  ///< with its trailing '/'
		double mu0 = 0.0;
		double referenceRadius = 0.0;
		/// The labels of a PRESCRIBED state ( reference/two, reference/three ).
		double psiAxis = 0.0;
		double psiBoundary = 0.0;
		/// Where a SOLVED state's labels are read ( reference/rotating ).
		double axisPointR = std::numeric_limits<double>::quiet_NaN();
		double axisPointZ = std::numeric_limits<double>::quiet_NaN();
		double boundaryPointR = std::numeric_limits<double>::quiet_NaN();
		double boundaryPointZ = std::numeric_limits<double>::quiet_NaN();
		/// The box end, |z| <= axialEdge; infinite where the case has none.
		double axialEdge = std::numeric_limits<double>::infinity();
		/// geq's grid box, and the count of nodes a side that every resolution
		/// it was solved at shares -- the nodes its psi files are written on.
		/// Zero where the case is not a solve.
		double gridRMin = 0.0;
		double gridRMax = 0.0;
		double gridZMin = 0.0;
		double gridZMax = 0.0;
		int gridCommon = 0;

		struct SpeciesEntry
		{
			std::string name;
			double mass = 0.0;
			double charge = 0.0;
			std::string temperatureFile;
			std::string densityFile;
		};

		std::vector<SpeciesEntry> species;
	};

	/// @param directory the case's directory, with or without a trailing '/'.
	/// @throws std::runtime_error naming meta.txt when it cannot be opened.
	inline GeqCase readCase( std::string directory )
	{
		if ( directory.empty() || directory.back() != '/' )
			directory += '/';
		GeqCase out;
		out.directory = directory;

		std::string const metaPath = directory + "meta.txt";
		std::ifstream meta( metaPath );
		if ( !meta )
			throw std::runtime_error(
				"cannot open " + metaPath + " -- the exported geq reference is missing. "
				"It is committed, so a build tree that lacks it is running from the "
				"wrong working directory; ctest pins CMAKE_SOURCE_DIR. To regenerate "
				"it see tools/geq-benchmark/README.md." );

		std::string line;
		while ( std::getline( meta, line ) )
		{
			std::istringstream in( line );
			std::string key;
			in >> key;

			if ( key == "R_ref" )
				in >> out.referenceRadius;
			else if ( key == "mu0" )
				in >> out.mu0;
			else if ( key == "psi_axis" )
				in >> out.psiAxis;
			else if ( key == "psi_boundary" )
				in >> out.psiBoundary;
			else if ( key == "axis_point" )
				in >> out.axisPointR >> out.axisPointZ;
			else if ( key == "boundary_point" )
				in >> out.boundaryPointR >> out.boundaryPointZ;
			else if ( key == "z_edge" )
				in >> out.axialEdge;
			else if ( key == "grid" )
			{
				// grid R <min> <max> z <min> <max> common <n>
				std::string label;
				in >> label >> out.gridRMin >> out.gridRMax
				   >> label >> out.gridZMin >> out.gridZMax
				   >> label >> out.gridCommon;
			}
			else if ( key == "name" )
			{
				GeqCase::SpeciesEntry entry;
				std::string label;
				in >> entry.name
				   >> label >> entry.mass
				   >> label >> entry.charge
				   >> label >> entry.temperatureFile
				   >> label >> entry.densityFile;
				out.species.push_back( entry );
			}
		}
		return out;
	}

	/// The species set in MEQ's gauge, or -- with @a raw -- in geq's, which is
	/// the mistake each export writes a second set of tables to make checkable.
	inline std::vector<meq::Species> speciesFrom( GeqCase const &reference, bool raw )
	{
		std::vector<meq::Species> species;

		for ( auto const &entry : reference.species )
		{
			std::string const density = raw
				? "Nraw_" + entry.name + ".dat"
				: entry.densityFile;

			meq::Species one;
			one.mass = entry.mass;
			one.charge = entry.charge;
			one.temperature = std::make_shared<meq::SplineProfile const>(
				meq::SplineProfile::fromFile( reference.directory + entry.temperatureFile ) );
			one.density = std::make_shared<meq::SplineProfile const>(
				meq::SplineProfile::fromFile( reference.directory + density ) );
			species.push_back( one );
		}

		return species;
	}

	inline std::shared_ptr<meq::Profile const> omegaFrom( GeqCase const &reference )
	{
		return std::make_shared<meq::SplineProfile const>(
			meq::SplineProfile::fromFile( reference.directory + "omega.dat" ) );
	}
}
}

#endif
