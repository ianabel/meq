#define BOOST_TEST_MODULE DriverAcceptance
#include <boost/test/unit_test.hpp>

#include "mfem.hpp"

#include "meq/BoundaryShape.hpp"
#include "meq/Estimator.hpp"
#include "meq/GradShafranov.hpp"
#include "meq/Output.hpp"
#include "meq/Profiles.hpp"
#include "meq/RotatingSource.hpp"
#include "meq/Source.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iterator>
#include <utility>
#include <memory>
#include <string>
#include <vector>

/*
 * THE DRIVER, END TO END.
 *
 * Everything else in this suite calls GradShafranovSolver directly. Until stage
 * 7e nothing exercised the path a USER takes -- config file in, files out --
 * and CLAUDE.md recorded for six stages that "meq is reachable only through its
 * test suite". This is the test that stops that being true.
 *
 * What it checks that the library tests cannot:
 *
 *   - the configuration actually reaches the solver. A driver that parsed
 *     PolynomialDegree and then built degree 1 anyway would pass every
 *     convergence test in the tree, because none of them goes through Config.
 *   - psi is written, and is the psi that was SOLVED FOR. The stored
 *     GridFunction is read back and compared against the same problem driven
 *     through the library directly, which is what catches a driver that wrote
 *     the flux into the potential's file, or wrote the iterate before
 *     RecoverFEMSolution, or flipped a sign on the way out.
 *   - the exit codes mean what they say. A shell script driving a parameter
 *     scan is a first-class caller, and "it printed something" is not an
 *     interface.
 *
 * WHY THE COMPARISON IS AGAINST THE LIBRARY AND NOT AGAINST THE CLOSED FORM.
 * The obvious test -- run the Solov'ev example, compare to Soloviev.hpp -- is
 * wrong, and measuring it is what showed that: it differs by 7.2e-1 in L2.
 * examples/soloviev-nstx.toml poses psi = 0 on the RECTANGLE, and the NSTX
 * Solov'ev's psi = 0 is the separatrix, a curved contour through an X-point.
 * Those are different problems and the driver is solving the one it was asked
 * for. Reproducing the closed form needs the extension path, which the driver
 * refuses to configure today rather than silently solving on the mesh boundary.
 *
 * So this asserts what a driver test can honestly assert: that the same
 * configuration reaches the same answer whether it goes through a TOML file or
 * through the API. SolovievConvergence is where the discretisation is measured
 * against a closed form, and it imposes the exact trace to do it.
 *
 */

namespace
{
	/// Where the driver binary is. CMake passes it, because the test cannot
	/// know the build directory's layout and guessing "../meq" is how a test
	/// starts silently not running the thing it claims to.
	char const *driver()
	{
		return MEQ_DRIVER_PATH;
	}

	int run( std::string const &configFile )
	{
		// MKL_NUM_THREADS=1 is inherited from the ctest environment, which
		// CMake sets for exactly the reason the other convergence tests record.
		// MKL_THREADING_LAYER is deliberately not set anywhere any more; see
		// tests/CMakeLists.txt.
		std::string const command = std::string( driver() ) + " " + configFile
		                            + " > /dev/null 2>&1";
		int const status = std::system( command.c_str() );
		return WIFEXITED( status ) ? WEXITSTATUS( status ) : -1;
	}

	bool exists( std::string const &path )
	{
		std::ifstream stream( path );
		return stream.good();
	}

	/// A file's whole text, for the ncdump reads below.
	std::string slurp( std::string const &path )
	{
		std::ifstream stream( path );
		return std::string( ( std::istreambuf_iterator<char>( stream ) ),
		                      std::istreambuf_iterator<char>() );
	}

	/// The header of a NetCDF file, through ncdump rather than through meq's own
	/// writer -- so the check does not share code with the thing it checks.
	/// Empty if ncdump could not read it.
	std::string ncdumpHeader( std::string const &path )
	{
		std::string const scratch = "driver-acceptance-header.txt";
		std::string const command = "ncdump -h " + path + " > " + scratch + " 2>&1";
		if ( std::system( command.c_str() ) != 0 )
			return std::string();

		std::string const text = slurp( scratch );
		std::remove( scratch.c_str() );
		return text;
	}

	/// One global attribute out of an ncdump header, as a double. NaN if it is
	/// not there, which every caller asserts against rather than ignores.
	double headerAttribute( std::string const &header, std::string const &name )
	{
		std::string const needle = ":" + name + " = ";
		std::size_t const at = header.find( needle );
		if ( at == std::string::npos )
			return std::nan( "" );
		return std::strtod( header.c_str() + at + needle.size(), nullptr );
	}

	/*
	 * THE ROTATING SPECIES OF examples/rotating-*.toml, BY HAND.
	 *
	 * Every number is repeated from the file rather than read from it, for the
	 * reason the Solov'ev case records: a driver that parsed TemperatureScale
	 * and then ignored it would agree with a test that also read the file, and
	 * disagree with this one. The two examples share their species exactly, so
	 * they share this.
	 */
	double const rotatingKeV = 1.602176634e-16;      // J per keV
	double const rotatingOmega = 4.0e5;              // rad/s
	double const rotatingGGPrime = 0.8;              // T^2 m^2 per (Wb/rad)
	double const rotatingNormalisedGGPrime = 0.08;   // per unit Psi
	double const rotatingDensityScale = 1.0e20;      // m^-3, the table's unit
	double const rotatingReferenceRadius = 1.0;      // m

	std::shared_ptr<meq::Profile const> constantProfile( double value )
	{
		return std::make_shared<meq::ConstantProfile const>( value );
	}

	/// A constant written as `value` with `scale` beside it, exactly as
	/// loadEitherProfile() builds it from the TOML: the scale is a wrapper and
	/// not a multiplication done at parse time, and this reproduces that rather
	/// than assuming the two are the same.
	std::shared_ptr<meq::Profile const> scaledConstant( double value, double scale )
	{
		return std::make_shared<meq::ScaledProfile const>( constantProfile( value ),
		                                                   scale );
	}

	/*
	 * THE DENSITY IS A TABLE AND SO IT IS READ, WHICH IS NOT THE SAME
	 * CONCESSION AS READING THE TOML.
	 *
	 * The rule this file works to is that the CONFIGURATION's numbers are
	 * repeated here rather than parsed, so that a driver which mis-plumbs a key
	 * disagrees with the test. The profile table is not configuration: it is the
	 * run's data, the driver and this test are both entitled to it, and the
	 * alternative -- transcribing five knots and their slopes -- would fail for
	 * a typo rather than for a defect. What stays repeated is everything the
	 * TOML says ABOUT the table, which is the file name and DensityScale; an
	 * ignored scale still fails here by a factor of 1e20.
	 */
	std::shared_ptr<meq::Profile const> densityTable( std::string const &file )
	{
		return std::make_shared<meq::ScaledProfile const>(
			std::make_shared<meq::SplineProfile const>(
				meq::SplineProfile::fromFile( file ) ),
			rotatingDensityScale );
	}

	/// @param densityFile  the table the matching TOML names. The two examples
	///                     differ in it and in nothing else about their species:
	///                     one is a function of psi and the other of Psi, which
	///                     is the whole difference between the two files.
	std::vector<meq::Species> rotatingSpecies( std::string const &densityFile )
	{
		std::vector<meq::Species> species( 2 );

		species[ 0 ].mass = 3.3435837768e-27;                 // deuterium, kg
		species[ 0 ].charge = 1.0;
		species[ 0 ].temperature = scaledConstant( 1.0, rotatingKeV );
		species[ 0 ].density = densityTable( densityFile );    // m^-3 on rRef

		species[ 1 ].mass = 9.1093837015e-31;                 // electron, kg
		species[ 1 ].charge = -1.0;
		species[ 1 ].temperature = scaledConstant( 0.8, rotatingKeV );
		// Neutralising = true in the file: the density is what charge neutrality
		// determines, and the driver must derive the same one.
		species[ 1 ].density = meq::neutralisingDensity( species, 1 );

		return species;
	}

	/// The mesh both rotating examples ask for: the standard rectangle, 8 x 12
	/// cells, one uniform refinement.
	mfem::Mesh rotatingMesh()
	{
		double const rMin = 0.6, rMax = 1.4, zMin = -0.6, zMax = 0.6;
		mfem::Mesh mesh = mfem::Mesh::MakeCartesian2D(
			8, 12, mfem::Element::TRIANGLE, false, rMax - rMin, zMax - zMin );
		mesh.Transform( [ rMin, zMin ]( mfem::Vector const &in, mfem::Vector &out )
		{
			out( 0 ) = in( 0 ) + rMin;
			out( 1 ) = in( 1 ) + zMin;
		} );
		mesh.UniformRefinement();
		return mesh;
	}

	/// The first and last values of a NetCDF coordinate variable, through
	/// ncdump. The `.nc` carries the grid extent in `R` and `Z` themselves
	/// rather than as attributes, so this is where a caller reads it from --
	/// and going through ncdump keeps the check off MEQ's own writer.
	std::pair<double, double> coordinateRange( std::string const &path,
	                                           std::string const &name )
	{
		std::string const scratch = "driver-acceptance-coord.txt";
		std::string const command = "ncdump -v " + name + " " + path
		                            + " > " + scratch + " 2>&1";
		if ( std::system( command.c_str() ) != 0 )
			return { std::nan( "" ), std::nan( "" ) };

		std::string text = slurp( scratch );
		std::remove( scratch.c_str() );

		// The data section, which is where the values are -- the header above
		// it also contains the variable's name.
		std::size_t at = text.find( "data:" );
		if ( at == std::string::npos ) { return { std::nan( "" ), std::nan( "" ) }; }
		at = text.find( "\n " + name + " = ", at );
		if ( at == std::string::npos ) { return { std::nan( "" ), std::nan( "" ) }; }
		at += name.size() + 5;
		std::size_t const end = text.find( ';', at );
		if ( end == std::string::npos ) { return { std::nan( "" ), std::nan( "" ) }; }

		std::string const values = text.substr( at, end - at );
		double const first = std::strtod( values.c_str(), nullptr );
		std::size_t const last = values.find_last_of( ',' );
		double const back = last == std::string::npos
			? first : std::strtod( values.c_str() + last + 1, nullptr );
		return { first, back };
	}

	/// A stored GridFunction, on a mesh the caller keeps alive.
	mfem::GridFunction readGridFunction( std::string const &path, mfem::Mesh &mesh )
	{
		std::ifstream stream( path );
		BOOST_TEST_REQUIRE( stream.good(), "cannot open " << path );
		return mfem::GridFunction( &mesh, stream );
	}

	/// Every occurrence, because a coil pair has two of each key and replacing
	/// one of them would leave a configuration nobody wrote.
	std::string replaceAll( std::string text, std::string const &from,
	                        std::string const &to )
	{
		for ( std::size_t at = text.find( from ); at != std::string::npos;
		      at = text.find( from, at + to.size() ) )
			text.replace( at, from.size(), to );
		return text;
	}

	/// ||stored - computed|| / ||computed||, the figure every case here prints.
	double relativeDifference( mfem::GridFunction const &stored,
	                           mfem::GridFunction const &computed )
	{
		mfem::Vector difference( stored );
		difference -= computed;
		return difference.Norml2()/std::max( 1.0e-300, computed.Norml2() );
	}

}

BOOST_AUTO_TEST_CASE( theDriverSolvesTheSolovievBenchmarkAndWritesIt )
{
	BOOST_TEST_REQUIRE( run( "examples/soloviev-nstx.toml" ) == 0,
	                    "the driver did not exit 0 on its own shipped example" );

	BOOST_TEST_REQUIRE( exists( "soloviev-nstx.mesh" ) );
	BOOST_TEST_REQUIRE( exists( "soloviev-nstx_psi.gf" ) );
	BOOST_TEST_REQUIRE( exists( "soloviev-nstx_grad_psi.gf" ) );

	mfem::Mesh storedMesh( "soloviev-nstx.mesh", 1, 1 );
	std::ifstream stream( "soloviev-nstx_psi.gf" );
	BOOST_TEST_REQUIRE( stream.good() );
	mfem::GridFunction stored( &storedMesh, stream );

	// The same problem, set up by hand. The numbers below are the example's,
	// repeated deliberately rather than read from it: a driver that parsed
	// PolynomialDegree and then built degree 1 anyway would agree with a test
	// that also read the file, and disagree with this one.
	mfem::Mesh mesh = mfem::Mesh::MakeCartesian2D(
		3, 4, mfem::Element::TRIANGLE, false, 1.9 - 0.1, 1.7 - ( -1.7 ) );
	mesh.Transform( []( mfem::Vector const &in, mfem::Vector &out )
	{
		out = in;
		out( 0 ) += 0.1;
		out( 1 ) += -1.7;
	} );
	for ( int i = 0; i < 3; ++i )
		mesh.UniformRefinement();

	meq::SolovievSource const source( -0.52 );
	mfem::ConstantCoefficient zero( 0.0 );

	meq::GradShafranovSolver solver( mesh, 3, 1.0 );
	solver.setSource( source );
	solver.setBoundaryData( zero );
	solver.solve();

	BOOST_TEST_REQUIRE( stored.Size() == solver.potential().Size(),
	                    "the driver wrote " << stored.Size() << " potential dofs "
	                    "where the same configuration gives "
	                    << solver.potential().Size() << ": the driver is not "
	                    "building the discretisation the file asks for" );

	mfem::Vector difference( stored );
	difference -= solver.potential();
	double const scale = std::max( 1.0e-300, solver.potential().Norml2() );
	double const relative = difference.Norml2()/scale;

	std::printf( "\n  driver vs library: ||psi_driver - psi_library|| / ||psi|| "
	             "= %.3e over %d dofs\n", relative, stored.Size() );
	std::fflush( stdout );

	// Both ran the same arithmetic in the same order, so this is round-off and
	// not a tolerance to be tuned. 1e-10 leaves room for the file round-trip,
	// which writes at precision 16.
	BOOST_TEST( relative < 1.0e-10,
	            "the driver's stored psi differs from the library's by "
	            << relative << " relative. Same mesh, same degree, same source, "
	            "so this is the driver mis-plumbing the configuration or writing "
	            "the wrong field" );
}

/*
 * The NON-LINEAR path through the driver.
 *
 * The Solov'ev case above is linear in psi -- dF/dpsi is identically zero and
 * Newton finishes in one step -- so it cannot tell a working Jacobian from an
 * absent one, and cannot tell that the driver reaches the Newton path at all.
 * Example 5's source depends on psi linearly, quadratically and exponentially.
 *
 * This asserts only that it runs and writes. The RATE work is
 * NewtonConvergence's, which imposes the exact trace to measure against the
 * closed form; what is being checked here is that a configuration file reaches
 * the non-linear solver rather than silently taking the linear branch.
 */
BOOST_AUTO_TEST_CASE( theDriverRunsTheNonlinearPath )
{
	BOOST_TEST_REQUIRE( run( "examples/manufactured-driver.toml" ) == 0,
	                    "the driver did not exit 0 on the non-linear example" );

	BOOST_TEST( exists( "manufactured-driver.mesh" ) );
	BOOST_TEST( exists( "manufactured-driver_psi.gf" ) );
	BOOST_TEST( exists( "manufactured-driver_grad_psi.gf" ) );

	// A non-trivial solution. A driver that took the linear branch on a
	// non-linear source, or that wrote before recovering the solution, would
	// leave this at zero.
	mfem::Mesh mesh( "manufactured-driver.mesh", 1, 1 );
	std::ifstream stream( "manufactured-driver_psi.gf" );
	BOOST_TEST_REQUIRE( stream.good() );
	mfem::GridFunction psi( &mesh, stream );
	BOOST_TEST( psi.Normlinf() > 1.0e-3,
	            "the stored psi is essentially zero, max |psi| = "
	            << psi.Normlinf() );
}

/*
 * THE CURVED BOUNDARY, THROUGH THE DRIVER.
 *
 * This is the configuration meq exists to run, and the one the driver refused
 * until it was wired. Gamma is a level set of the Miller shape; the background
 * mesh knows nothing about it; D_h is the union of background elements inside
 * Gamma; and the datum is transferred onto Gamma_h along short paths.
 *
 * As with the fitted case, the comparison is against the LIBRARY on the same
 * configuration, and the subdomain is rebuilt here by hand rather than read
 * from the driver. That is the point: it checks the driver marks the same
 * elements, finds the same Gamma_h attribute, builds the same paths and passes
 * the same marker to setExtension(). Getting any of those wrong gives a
 * plausible-looking answer -- CLAUDE.md records that a Gamma_h with no
 * transferred datum silently imposes zero, and that the one configuration which
 * really fails reaches 5e13 rather than failing loudly.
 *
 * The rates for the technique itself are ExtensionConvergence's, against a
 * closed form. This is a plumbing test.
 */
BOOST_AUTO_TEST_CASE( theDriverSolvesOnACurvedBoundary )
{
	BOOST_TEST_REQUIRE( run( "examples/miller-curved.toml" ) == 0,
	                    "the driver did not exit 0 on the curved example" );

	BOOST_TEST_REQUIRE( exists( "miller-curved.mesh" ) );
	BOOST_TEST_REQUIRE( exists( "miller-curved_psi.gf" ) );

	mfem::Mesh storedMesh( "miller-curved.mesh", 1, 1 );
	std::ifstream stream( "miller-curved_psi.gf" );
	BOOST_TEST_REQUIRE( stream.good() );
	mfem::GridFunction stored( &storedMesh, stream );

	// The same run, by hand. Numbers repeated from the example deliberately.
	mfem::Mesh background = mfem::Mesh::MakeCartesian2D(
		8, 10, mfem::Element::TRIANGLE, false, 2.3 - 0.7, 1.9 - ( -1.9 ) );
	background.Transform( []( mfem::Vector const &in, mfem::Vector &out )
	{
		out( 0 ) = in( 0 ) + 0.7;
		out( 1 ) = in( 1 ) + ( -1.9 );
	} );
	for ( int i = 0; i < 2; ++i )
		background.UniformRefinement();

	meq::BoundaryShape const shape =
		meq::BoundaryShape::miller( 1.5, 0.0, 0.5, 1.6, 0.35, 0.0 );
	mfem::PositionFunction const levelSet = [ &shape ]( mfem::Vector const &x )
	{
		return shape.levelSet( x( 0 ), x( 1 ) );
	};

	mfem::Array<int> marker;
	int const insideCount =
		mfem::MarkLevelSetSubdomain( background, levelSet, 0.0, marker, 1 );
	BOOST_TEST_REQUIRE( insideCount > 0 );

	for ( int e = 0; e < background.GetNE(); ++e )
		background.SetAttribute( e, marker[ e ] ? 1 : 2 );
	background.SetAttributes();

	mfem::Array<int> domainAttribute( 1 );
	domainAttribute[ 0 ] = 1;
	mfem::SubMesh sub = mfem::SubMesh::CreateFromDomain( background, domainAttribute );

	BOOST_TEST_REQUIRE( sub.GetNE() == storedMesh.GetNE(),
	                    "the driver solved on " << storedMesh.GetNE()
	                    << " elements where the same shape gives " << sub.GetNE()
	                    << ": it is not marking the same subdomain" );

	int const gammaH = sub.bdr_attributes.Max();
	double const h = std::max( ( 2.3 - 0.7 )/( 8*4.0 ), ( 1.9 + 1.9 )/( 10*4.0 ) );
	mfem::VertexConePath path( sub, gammaH, levelSet, 6.0*h );

	mfem::Array<int> gammaHMarker( gammaH );
	gammaHMarker = 0;
	gammaHMarker[ gammaH - 1 ] = 1;

	meq::SolovievSource const source( -0.52 );
	mfem::ConstantCoefficient zero( 0.0 );

	meq::GradShafranovSolver solver( sub, 2, 1.0 );
	solver.setSource( source );
	solver.setBoundaryData( zero );
	solver.setExtension( path, gammaHMarker );
	solver.solve();

	BOOST_TEST_REQUIRE( stored.Size() == solver.potential().Size() );

	mfem::Vector difference( stored );
	difference -= solver.potential();
	double const scale = std::max( 1.0e-300, solver.potential().Norml2() );
	double const relative = difference.Norml2()/scale;

	std::printf( "\n  curved driver vs library: %.3e relative over %d dofs, "
	             "%d of %d elements inside\n",
	             relative, stored.Size(), sub.GetNE(), background.GetNE() );
	std::fflush( stdout );

	BOOST_TEST( relative < 1.0e-10,
	            "the driver's curved solve differs from the library's by "
	            << relative << " relative" );

	// And it must not be the fitted answer wearing a curved hat. If
	// setExtension() were never called the solve would impose zero on Gamma_h
	// and still converge, so this is the assertion that says the transfer
	// happened at all.
	meq::GradShafranovSolver fitted( sub, 2, 1.0 );
	fitted.setSource( source );
	fitted.setBoundaryData( zero );
	fitted.solve();

	mfem::Vector fittedDifference( solver.potential() );
	fittedDifference -= fitted.potential();
	double const fittedRelative = fittedDifference.Norml2()/scale;
	std::printf( "  extension vs zero-on-Gamma_h: %.3e relative\n", fittedRelative );
	std::fflush( stdout );

	BOOST_TEST( fittedRelative > 1.0e-6,
	            "the curved solve agrees with imposing zero on Gamma_h to "
	            << fittedRelative << ", so the transferred datum is not being "
	            "applied and the extension is inert" );
}

/*
 * Exit codes, which are the driver's actual interface to a scan script.
 *
 * 1 is checked and not 2 or 3: a solve that fails to converge needs a source
 * that does not converge, and the only one meq has is GS-2 section 4.4, which
 * takes minutes. That belongs in a slower test than this one if it is ever
 * worth pinning.
 */
/*
 * THE ADAPTIVE LOOP, THROUGH THE DRIVER, ON A CURVED BOUNDARY.
 *
 * This is the last thing the driver refused to do, and it is the configuration
 * meq exists to run: Gamma is a level set, the mesh is not fitted to it, and
 * which elements to refine is decided by the residual estimator rather than by
 * the person writing the TOML.
 *
 * WHAT IS PINNED, AND WHY IT IS THE LIBRARY AGAIN. The same reasoning as
 * theDriverSolvesOnACurvedBoundary: examples/miller-adaptive.toml poses a
 * problem with no closed form -- the Solov'ev source on a Miller boundary is not
 * the Solov'ev equilibrium -- so what can honestly be asserted is that the loop
 * reaches the same place whether it is driven from a file or from the API. The
 * discretisation is measured against closed forms elsewhere; the ADAPTIVE
 * MACHINERY is measured in AdaptiveRefinement.cpp, which asserts that eta and
 * the true error both come down and that the proximity condition holds.
 *
 * THREE THINGS HERE THAT THE LIBRARY TESTS DO NOT COVER, each of which is a way
 * for a driver to pass everything else while doing the wrong thing:
 *
 *   - it must REFINE. A driver that parsed [adaptivity], solved once and wrote
 *     the answer would produce a correct file and satisfy every assertion about
 *     agreement, because agreement with one library cycle is agreement. The
 *     element count is what catches that.
 *   - it must refine ADAPTIVELY. Marking every element is uniform refinement
 *     wearing an estimator, and it is what maximum marking degenerates to at
 *     small gamma -- so the marked count must stay under the element count.
 *   - it must use the COMPANION mesh on the curved path. Refining D_h alone
 *     leaves Gamma_h where it is while h_loc halves, and the transfer stops
 *     being optimal. The hand-rolled loop below uses meq::AdaptiveDomain, so a
 *     driver that reached for a plain SubMesh would disagree at once.
 */
BOOST_AUTO_TEST_CASE( theDriverRunsTheAdaptiveLoop )
{
	BOOST_TEST_REQUIRE( run( "examples/miller-adaptive.toml" ) == 0,
	                    "the driver did not exit 0 on the adaptive example" );

	BOOST_TEST_REQUIRE( exists( "miller-adaptive.mesh" ) );
	BOOST_TEST_REQUIRE( exists( "miller-adaptive_psi.gf" ) );

	mfem::Mesh storedMesh( "miller-adaptive.mesh", 1, 1 );
	std::ifstream stream( "miller-adaptive_psi.gf" );
	BOOST_TEST_REQUIRE( stream.good() );
	mfem::GridFunction stored( &storedMesh, stream );

	// The same run, by hand. Numbers repeated from the example deliberately, so
	// that editing the example without editing this test is a failure and not a
	// silent divergence.
	int const order = 2;
	int const cycles = 4;
	double const theta = 0.6;
	double const rMin = 0.7, rMax = 2.3, zMin = -1.9, zMax = 1.9;

	mfem::Mesh background = mfem::Mesh::MakeCartesian2D(
		8, 10, mfem::Element::TRIANGLE, false, rMax - rMin, zMax - zMin );
	background.Transform( [ rMin, zMin ]( mfem::Vector const &in, mfem::Vector &out )
	{
		out( 0 ) = in( 0 ) + rMin;
		out( 1 ) = in( 1 ) + zMin;
	} );
	background.UniformRefinement();      // RefinementLevels = 1

	meq::BoundaryShape const shape =
		meq::BoundaryShape::miller( 1.5, 0.0, 0.5, 1.6, 0.35, 0.0 );
	mfem::PositionFunction const levelSet = [ &shape ]( mfem::Vector const &x )
	{
		return shape.levelSet( x( 0 ), x( 1 ) );
	};

	meq::SolovievSource const source( -0.52 );
	mfem::ConstantCoefficient zero( 0.0 );

	meq::AdaptiveDomain domain( background, levelSet );

	struct Turn { int elements; int marked; int widened; double eta; };
	std::vector<Turn> history;
	std::unique_ptr<meq::GradShafranovSolver> solver;
	std::unique_ptr<mfem::VertexConePath> path;

	for ( int cycle = 0; cycle < cycles; ++cycle )
	{
		mfem::Array<int> marked;
		Turn turn{ 0, 0, 0, -1.0 };

		// Twelve times the LARGEST element, as the driver uses: on a graded mesh
		// the coarse part needs the long search.
		path = std::make_unique<mfem::VertexConePath>(
			domain.computational(), domain.gammaHAttribute(), levelSet,
			12.0*domain.largestElement() );

		solver = std::make_unique<meq::GradShafranovSolver>(
			domain.computational(), order, 1.0 );
		solver->setSource( source );
		solver->setBoundaryData( zero );
		solver->setExtension( *path, domain.gammaHMarker() );
		solver->solve();

		// The driver's own choice of potential, reproduced rather than assumed:
		// psi*, the published estimator, which it could not use until MFEM's
		// reconstruction stopped skipping its mean-value close. See apps/meq.cpp.
		solver->postProcess();
		meq::ResidualEstimator estimator( *solver, source );

		// And the driver's treatment of Gamma_h, which is no longer to leave
		// those faces out: eta_5 compares psi* against the datum actually
		// imposed. Rebuilt per cycle for the reason apps/meq.cpp gives -- the
		// datum lifts the SOLVED flux, so last cycle's is the wrong boundary
		// condition. If this drifts from the driver the two etas part company
		// and this test says so, which is what it is for.
		std::unique_ptr<mfem::Coefficient> const datum =
			solver->transferredDatum();
		estimator.setTransferredBoundary( domain.gammaHMarker(), datum.get() );

		mfem::Vector const &local = estimator.GetLocalErrors();
		turn.elements = domain.numComputational();
		turn.widened = path->NumWidened();
		turn.eta = estimator.GetTotalError();

		if ( cycle + 1 < cycles )
		{
			meq::markDoerfler( local, theta, marked );
			turn.marked = marked.Size();
		}
		history.push_back( turn );

		if ( cycle + 1 == cycles )
			break;

		BOOST_TEST_REQUIRE( marked.Size() > 0,
		                    "nothing was marked at cycle " << cycle );
		solver.reset();
		path.reset();
		domain.refine( marked );
	}

	std::printf( "\n  the driver's adaptive loop, reproduced by hand\n" );
	std::printf( "  %6s %8s %8s %6s %12s\n", "cycle", "elem", "marked", "wide", "eta" );
	for ( std::size_t c = 0; c < history.size(); ++c )
		std::printf( "  %6zu %8d %8d %6d %12.4e\n", c, history[ c ].elements,
		             history[ c ].marked, history[ c ].widened, history[ c ].eta );
	std::fflush( stdout );

	// IT REFINED, AND IT REFINED ADAPTIVELY.
	BOOST_TEST( history.back().elements > 2*history.front().elements,
	            "the loop went from " << history.front().elements << " elements to "
	            << history.back().elements << ", which is not a refinement worth "
	            "four cycles -- check that the driver is not solving once and "
	            "calling it a loop" );
	for ( std::size_t c = 0; c + 1 < history.size(); ++c )
		BOOST_TEST( history[ c ].marked < history[ c ].elements,
		            "cycle " << c << " marked all " << history[ c ].elements
		            << " elements, which is uniform refinement" );

	// ETA CAME DOWN. Monotone here, and it is worth knowing that it is: the loop
	// would still be a loop without this, and it would be refining on an
	// indicator that is telling it nothing.
	for ( std::size_t c = 1; c < history.size(); ++c )
		BOOST_TEST( history[ c ].eta < history[ c - 1 ].eta,
		            "eta went from " << history[ c - 1 ].eta << " to "
		            << history[ c ].eta << " at cycle " << c );

	// ASSUMPTION P.1, on the graded Gamma_h the loop produces. VertexConePath
	// widens its fan rather than failing when a ray finds no root, and the method
	// still runs -- but the Cockburn-Solano estimate no longer covers it, so the
	// count is asserted rather than reported.
	for ( std::size_t c = 0; c < history.size(); ++c )
		BOOST_TEST( history[ c ].widened == 0,
		            "cycle " << c << ": " << history[ c ].widened
		            << " vertices of Gamma_h needed a widened fan" );

	// AND THE DRIVER LANDED IN THE SAME PLACE.
	BOOST_TEST_REQUIRE( storedMesh.GetNE() == history.back().elements,
	                    "the driver wrote " << storedMesh.GetNE()
	                    << " elements where the same loop by hand ends at "
	                    << history.back().elements
	                    << ": the refinement sequences have diverged" );
	BOOST_TEST_REQUIRE( stored.Size() == solver->potential().Size() );

	mfem::Vector difference( stored );
	difference -= solver->potential();
	double const scale = std::max( 1.0e-300, solver->potential().Norml2() );
	double const relative = difference.Norml2()/scale;

	std::printf( "  adaptive driver vs library: %.3e relative over %d dofs\n",
	             relative, stored.Size() );
	std::fflush( stdout );

	/*
	 * AND THE TWO NO LONGER START THEIR CYCLES FROM THE SAME ITERATE, WHICH MAKES
	 * THIS A STRONGER STATEMENT THAN IT WAS RATHER THAN A WEAKER ONE.
	 *
	 * Since 2026-09-02 the driver interpolates each cycle's answer onto the
	 * refined mesh and starts the next one from it; the loop above deliberately
	 * still starts cold. So this no longer reads 1.7e-16 -- it reads about
	 * 4e-14 -- and what it now asserts is the property a warm start has to have:
	 * that it changes the WORK and not the ANSWER. The library's own
	 * aWarmStartCutsTheWorkAndNotTheAnswer says the same thing one level down.
	 *
	 * Mirroring the warm start here would restore bitwise agreement and assert
	 * less, since a difference in the starting iterate is exactly what this is
	 * now able to see.
	 */
	BOOST_TEST( relative < 1.0e-10,
	            "the driver's adaptive solve differs from the library's by "
	            << relative << " relative. The driver warm-starts each cycle and "
	            "this loop does not, so a difference here is the warm start "
	            "changing the answer rather than only the work" );
}

/*
 * A ROTATING PLASMA, THROUGH THE DRIVER.
 *
 * The rotating source is a bigger piece of configuration
 * than anything else meq takes -- an array of tables, two profiles per species,
 * a scale on each, a derived density, a reference radius and a rotation
 * frequency -- and every one of those is a way for the file to reach the solver
 * as something other than what it says.
 *
 * THREE THINGS ARE PINNED HERE THAT NOTHING ELSE CAN SEE.
 *
 * psi ITSELF, against the same problem built by hand. That catches the whole
 * configuration path: a TemperatureScale or a DensityScale ignored, a
 * Neutralising species given the wrong sign, ReferenceRadius defaulted, Omega
 * dropped because toml11's find_or<double> returns the default for an integer
 * node, the wrong density table read. All of those reach psi, because the
 * example's n_D0 has a slope in psi and so mu0 r^2 dp/dpsi does not vanish --
 * which was not true of the first version of this example and is the whole
 * reason it was changed.
 *
 * THAT ROTATION REACHES psi AT ALL, by solving the same problem at omega = 0 and
 * requiring the two to differ. Every other assertion in this case passes whether
 * or not it does; see the block below.
 *
 * THE OUTPUT FIELDS, which are the new thing FL-8 adds. n_s and phi_0 are the
 * whole content of (96) and (97), and a rotating equilibrium is not
 * interpretable without them.
 */
BOOST_AUTO_TEST_CASE( theDriverSolvesARotatingEquilibrium )
{
	BOOST_TEST_REQUIRE( run( "examples/rotating-rectangle.toml" ) == 0,
	                    "the driver did not exit 0 on the rotating example" );

	BOOST_TEST_REQUIRE( exists( "rotating-rectangle.mesh" ) );
	BOOST_TEST_REQUIRE( exists( "rotating-rectangle_psi.gf" ) );

	mfem::Mesh storedMesh( "rotating-rectangle.mesh", 1, 1 );
	std::ifstream stream( "rotating-rectangle_psi.gf" );
	BOOST_TEST_REQUIRE( stream.good() );
	mfem::GridFunction stored( &storedMesh, stream );

	// The same run, by hand, with the numbers repeated from the file.
	mfem::Mesh mesh = rotatingMesh();
	meq::RotatingSource const source(
		rotatingSpecies( "examples/rotating-density.dat" ),
		constantProfile( rotatingOmega ),
		constantProfile( rotatingGGPrime ),
		rotatingReferenceRadius );

	mfem::ConstantCoefficient zero( 0.0 );
	meq::GradShafranovSolver solver( mesh, 2, 1.0 );
	solver.setSource( source );
	solver.setBoundaryData( zero );
	solver.setNewtonControl( 1.0e-10, 1.0e-14, 20 );
	solver.solve();

	BOOST_TEST_REQUIRE( stored.Size() == solver.potential().Size(),
	                    "the driver wrote " << stored.Size() << " potential dofs "
	                    "where the same configuration gives "
	                    << solver.potential().Size() );

	double const relative = relativeDifference( stored, solver.potential() );
	std::printf( "\n  rotating driver vs library: %.3e relative over %d dofs, "
	             "max psi = %.6e Wb/rad\n",
	             relative, stored.Size(), solver.potential().Max() );
	std::fflush( stdout );

	BOOST_TEST( relative < 1.0e-10,
	            "the driver's rotating solve differs from the library's by "
	            << relative << " relative" );

	/*
	 * AND ROTATION REACHED psi, WHICH IS THE ONE THING THE EXAMPLE EXISTS TO
	 * SHOW AND THE ONE THING EVERY ASSERTION ABOVE WOULD PASS WITHOUT.
	 *
	 * The driver-against-library check is a plumbing check and is satisfied by
	 * any source at all. Nothing in it notices if the equilibrium is the
	 * NON-ROTATING one -- which is exactly what happens when every profile is a
	 * constant: p is then a function of r alone, dp/dpsi vanishes, F reduces to
	 * g g', and psi comes out BIT-IDENTICAL to the plasma at rest while n_s and
	 * phi_0 still look convincingly centrifugal. The example was written that way
	 * first. Its header quotes the number below, so the number is measured here
	 * rather than asserted there, and a later edit that shrinks the pressure term
	 * -- raising GGPrime, flattening the density -- fails instead of quietly
	 * making the example vacuous again.
	 *
	 * Same species, same mesh, same everything, and no rotation.
	 */
	meq::RotatingSource const still(
		rotatingSpecies( "examples/rotating-density.dat" ),
		nullptr,                              // omega = 0
		constantProfile( rotatingGGPrime ),
		rotatingReferenceRadius );

	meq::GradShafranovSolver rest( mesh, 2, 1.0 );
	rest.setSource( still );
	rest.setBoundaryData( zero );
	rest.setNewtonControl( 1.0e-10, 1.0e-14, 20 );
	rest.solve();

	double const shift = relativeDifference( solver.potential(), rest.potential() );
	std::printf( "  rotation moves psi by %.4e relative in L2 over %d dofs; "
	             "max psi %.6e spinning against %.6e at rest\n",
	             shift, stored.Size(), solver.potential().Max(),
	             rest.potential().Max() );
	std::fflush( stdout );

	BOOST_TEST( shift > 5.0e-2,
	            "psi differs from the non-rotating answer by only " << shift
	            << " relative. The example is then documenting its own vacuity: "
	            "a reader who copies it and edits Omega will see nothing move. "
	            "The cause is the pressure term being small beside GGPrime -- "
	            "check that the density profile still has a slope in psi" );

	if ( !meq::hasNetCDF() )
	{
		BOOST_TEST_MESSAGE( "  built without netcdf-cxx4, skipping the fields" );
		return;
	}

	// THE FIELDS. Declared with the names the species were given, and with
	// source_type saying which physics produced the file -- without which two
	// runs differing only in [source] Type have identical headers.
	std::string const header = ncdumpHeader( "rotating-rectangle.nc" );
	BOOST_TEST_REQUIRE( !header.empty(),
	                    "ncdump could not read rotating-rectangle.nc" );

	for ( std::string const &needle :
	      { std::string( "double n_D(Z, R)" ), std::string( "double n_e(Z, R)" ),
	        std::string( "double e_phi_0(Z, R)" ),
	        std::string( "n_D:units = \"m^-3\"" ),
	        std::string( "e_phi_0:units = \"J\"" ),
	        std::string( ":source_type = \"rotating\"" ) } )
		BOOST_TEST( header.find( needle ) != std::string::npos,
		            "the file does not carry '" << needle << "'" );
}

/*
 * NORMALISED FLUX THROUGH THE DRIVER, WHICH IS A DIFFERENT SOLVE AND NOT A
 * DIFFERENT SET OF NUMBERS.
 *
 * Psi = psi/psi_ax makes psi_ax a functional of the solution, so the driver has
 * to reach setSource( NormalisedSource &, double ) rather than
 * setSource( Source const & ), and the solver closes the pair by a bordered
 * Newton. A driver that took the ordinary branch would still converge -- to the
 * equilibrium at whatever psi_ax the file happened to guess, which is a solved
 * equation for a plasma nobody asked for.
 *
 * SO psi_ax IS PINNED TWICE: against the library's own bordered Newton, and
 * against the value the driver WROTE to the file. The second is not redundant.
 * psi_ax is an answer here rather than an input, so a reader has nowhere else to
 * get it -- and a driver that solved correctly and then wrote the TOML's guess
 * into the attribute would pass every other assertion in this file.
 */
BOOST_AUTO_TEST_CASE( theDriverSolvesForPsiAxisAsAnUnknown )
{
	BOOST_TEST_REQUIRE( run( "examples/rotating-normalised.toml" ) == 0,
	                    "the driver did not exit 0 on the normalised example" );

	BOOST_TEST_REQUIRE( exists( "rotating-normalised.mesh" ) );
	BOOST_TEST_REQUIRE( exists( "rotating-normalised_psi.gf" ) );

	mfem::Mesh storedMesh( "rotating-normalised.mesh", 1, 1 );
	std::ifstream stream( "rotating-normalised_psi.gf" );
	BOOST_TEST_REQUIRE( stream.good() );
	mfem::GridFunction stored( &storedMesh, stream );

	// The same run, by hand. PsiAxis = 0.3 is the file's GUESS, and the point of
	// the case is that the answer is not it.
	double const psiAxisGuess = 0.08;
	mfem::Mesh mesh = rotatingMesh();
	meq::NormalisedRotatingSource source(
		rotatingSpecies( "examples/rotating-density-normalised.dat" ),
		constantProfile( rotatingOmega ),
		constantProfile( rotatingNormalisedGGPrime ),
		rotatingReferenceRadius, psiAxisGuess );

	mfem::ConstantCoefficient zero( 0.0 );
	meq::GradShafranovSolver solver( mesh, 2, 1.0 );
	solver.setSource( source, psiAxisGuess );
	solver.setBoundaryData( zero );
	solver.setNewtonControl( 1.0e-10, 1.0e-14, 20 );
	solver.solve();

	BOOST_TEST_REQUIRE( stored.Size() == solver.potential().Size() );

	double const relative = relativeDifference( stored, solver.potential() );
	std::printf( "\n  normalised driver vs library: %.3e relative over %d dofs\n"
	             "  psi_ax = %.9e, constraint psi_ax - max psi_h = %.3e, "
	             "%d Newton steps\n",
	             relative, stored.Size(), solver.psiAxis(),
	             solver.normalisationResidual(), solver.newtonIterations() );
	std::fflush( stdout );

	BOOST_TEST( relative < 1.0e-10,
	            "the driver's normalised solve differs from the library's by "
	            << relative << " relative" );

	// It solved FOR psi_ax rather than accepting the guess. Without this the
	// case above would pass with the border doing nothing at all.
	BOOST_TEST( std::fabs( solver.psiAxis() - psiAxisGuess ) > 1.0e-2*psiAxisGuess,
	            "psi_ax came back at the initial guess " << psiAxisGuess
	            << ", so nothing solved for it" );

	if ( !meq::hasNetCDF() )
	{
		BOOST_TEST_MESSAGE( "  built without netcdf-cxx4, skipping the attributes" );
		return;
	}

	std::string const header = ncdumpHeader( "rotating-normalised.nc" );
	BOOST_TEST_REQUIRE( !header.empty(),
	                    "ncdump could not read rotating-normalised.nc" );

	double const written = headerAttribute( header, "psi_axis" );
	double const constraint = headerAttribute( header, "normalisation_residual" );

	std::printf( "  the file says psi_axis = %.9e, "
	             "normalisation_residual = %.3e\n", written, constraint );
	std::fflush( stdout );

	BOOST_TEST( std::isfinite( written ),
	            "the file carries no psi_axis attribute, so a reader has no way "
	            "to know the axis flux of an equilibrium whose psi_ax was an "
	            "unknown" );
	BOOST_TEST( std::isfinite( constraint ),
	            "the file carries no normalisation_residual attribute" );

	// The attribute is written at full precision, so this is the file round trip
	// and not a tolerance.
	BOOST_TEST( std::fabs( written - solver.psiAxis() )
	                <= 1.0e-12*std::fabs( solver.psiAxis() ),
	            "the file says psi_axis = " << written << " where the same run "
	            "gives " << solver.psiAxis() );
	BOOST_TEST( std::fabs( constraint ) <= 1.0e-8*std::fabs( solver.psiAxis() ),
	            "the driver reports a normalisation residual of " << constraint
	            << " against psi_ax = " << solver.psiAxis()
	            << ", so the constraint psi_ax = max psi_h is not satisfied and "
	            "the border is not closing the system" );
}

/*
 * THE COILS, THROUGH THE DRIVER.
 *
 * `[[coils]]` puts a current that is DATA into F. Nothing else in this suite
 * goes through that path from a file, and the two things it can get wrong are
 * both quiet:
 *
 *   * THE COIL TERM NEVER REACHING THE ASSEMBLY. F is built by quadrature over
 *     the elements, so a coil set parsed and then dropped changes nothing at
 *     all -- the run converges, writes its files, and reports `coil_current` in
 *     the `.nc` whether or not that current did any work. The first half below
 *     is the check that it did.
 *   * THE WRAPPER PERTURBING SOMETHING ELSE. Carrying the coils means the
 *     driver hands the solver a meq::CoilAugmentedSource around the plasma
 *     source rather than the plasma source itself, and that is a change to
 *     EVERY run with a coil block in it. The second half pins zero-current
 *     coils against no coils at all and requires BIT IDENTITY -- 0.000e+00 and
 *     not "agrees to round-off", because an empty contribution is exactly
 *     empty and anything else is the wrapper doing arithmetic it should not.
 *
 * The equilibrium itself is not compared against anything: examples/
 * coils-rectangle.toml is a fixed-boundary rectangle with the conductors inside
 * the plasma, so there is no closed form to hold it to. What the coils are
 * worth on a SOLVE is measured in FreeBoundaryCoupling.cpp, against Ampere's
 * law, on the half-disc geometry the technique is for.
 */
BOOST_AUTO_TEST_CASE( theDriverAddsTheCoilsToF )
{
	BOOST_TEST_REQUIRE( run( "examples/coils-rectangle.toml" ) == 0,
	                    "the driver did not exit 0 on its own shipped coil example" );

	// The provenance, so a reader differencing two files can tell a run with
	// conductors from one without.
	std::string const header = ncdumpHeader( "coils-rectangle.nc" );
	BOOST_TEST_REQUIRE( !header.empty(), "ncdump could not read coils-rectangle.nc" );
	BOOST_TEST( headerAttribute( header, "coils" ) == 2.0 );
	// 1.5e5 given as Current plus 1.5e7 A/m^2 over 0.01 m^2 as CurrentDensity,
	// so this also says the two spellings resolve to the same number.
	BOOST_TEST( headerAttribute( header, "coil_current" ) == 3.0e5,
	            boost::test_tools::tolerance( 1.0e-12 ) );

	mfem::Mesh mesh( "coils-rectangle.mesh", 1, 1 );
	mfem::GridFunction driven = readGridFunction( "coils-rectangle_psi.gf", mesh );

	/*
	 * The same file with the two currents zeroed, and the same file with the
	 * coil blocks deleted outright. Written out here rather than shipped,
	 * because they are controls rather than examples -- and edited by text
	 * substitution on the shipped file so that they cannot drift away from it.
	 */
	std::string const shipped = slurp( "examples/coils-rectangle.toml" );
	BOOST_TEST_REQUIRE( !shipped.empty() );

	auto write = []( std::string const &path, std::string text,
	                 std::string const &prefix )
	{
		text = replaceAll( text, "Prefix = \"coils-rectangle\"",
		                   "Prefix = \"" + prefix + "\"" );
		std::ofstream file( path );
		file << text;
		return file.good();
	};

	std::string zeroed = replaceAll( shipped, "\nCurrent = 1.5e5", "\nCurrent = 0.0" );
	zeroed = replaceAll( zeroed, "\nCurrentDensity = 1.5e7", "\nCurrentDensity = 0.0" );
	BOOST_TEST_REQUIRE( zeroed != shipped, "the current substitution matched nothing" );

	// LINE-ANCHORED, because the file's own prose mentions [[coils]] and
	// [boundary] several times before either table appears -- the first
	// version of this cut the file at a COMMENT and produced a configuration
	// nobody wrote, which the driver then refused with exit 1.
	std::string stripped = shipped.substr( 0, shipped.find( "\n[[coils]]\n" ) + 1 )
	                       + shipped.substr( shipped.find( "\n[boundary]\n" ) + 1 );
	BOOST_TEST_REQUIRE( stripped.find( "\n[[coils]]\n" ) == std::string::npos );
	BOOST_TEST_REQUIRE( stripped.find( "\n[boundary]\n" ) != std::string::npos );

	BOOST_TEST_REQUIRE( write( "driver-acceptance-zerocoil.toml", zeroed, "zerocoil" ) );
	BOOST_TEST_REQUIRE( write( "driver-acceptance-nocoil.toml", stripped, "nocoil" ) );

	BOOST_TEST_REQUIRE( run( "driver-acceptance-zerocoil.toml" ) == 0 );
	BOOST_TEST_REQUIRE( run( "driver-acceptance-nocoil.toml" ) == 0 );

	mfem::Mesh zeroMesh( "zerocoil.mesh", 1, 1 );
	mfem::Mesh noMesh( "nocoil.mesh", 1, 1 );
	mfem::GridFunction zeroPsi = readGridFunction( "zerocoil_psi.gf", zeroMesh );
	mfem::GridFunction noPsi = readGridFunction( "nocoil_psi.gf", noMesh );

	// AN EMPTY CONTRIBUTION IS EXACTLY EMPTY.
	mfem::Vector control( zeroPsi );
	control -= noPsi;
	BOOST_TEST( control.Normlinf() == 0.0,
	            "coils carrying no current changed the answer by "
	            << control.Normlinf() << ", so the wrapper is not inert" );

	// AND A REAL CURRENT DOES REAL WORK. 3.0e5 A moves psi by 22% in L2
	// on this configuration; the gate is well below that and well above the
	// round-off the control sits at.
	double const moved = relativeDifference( driven, noPsi );
	BOOST_TEST( moved > 1.0e-2,
	            "300 kA of coil current moved psi by only " << moved
	            << " relative, so the coil term is not reaching F" );
	std::printf( "\n  the coils through the driver\n"
	             "    zero current vs no coils   %.3e   (must be exactly zero)\n"
	             "    300 kA vs no coils         %.3e   relative in L2\n",
	             control.Normlinf(), moved );

	for ( char const *path : { "driver-acceptance-zerocoil.toml",
	                           "driver-acceptance-nocoil.toml" } )
		std::remove( path );
}

/*
 * A MESH FROM A FILE, AND THE GRID EXTENT THAT USED TO BE LOST WITH IT.
 *
 * `[mesh] File` supplies no `RMin`..`ZMax`, and those four keys are what the
 * gridded output samples over -- so the extent came out empty, `meq::GridSampler`
 * refused, and a run that had ALREADY SOLVED was lost at the output stage. The
 * driver takes the mesh's own bounding box instead.
 *
 * NOT A CORNER CASE. The half-disc reaching the axis that free boundary needs
 * cannot come from `MakeCartesian2D` at all -- it is a semicircle centred on
 * r = 0, and the exterior expansion is a statement about exactly that geometry
 * -- so a free-boundary run ALWAYS reads its mesh from a file. See
 * tools/mesh/README.md and tools/mesh/halfdisc.py, which generates it.
 *
 * The mesh here is written by MFEM rather than by gmsh, deliberately: what is
 * being tested is the driver's handling of a mesh whose extent it did not
 * choose, and generating one needs gmsh, which the test suite does not require.
 * The gmsh path is exercised by `halfdisc.py --check`, which re-reads its own
 * output.
 */
BOOST_AUTO_TEST_CASE( theDriverTakesItsGridFromAMeshItDidNotBuild )
{
	// A box that is NOT at the origin and is NOT square, so a driver that fell
	// back to a default extent, or transposed the two directions, disagrees.
	double const rMin = 1.3, rMax = 2.1, zMin = -0.55, zMax = 0.75;
	{
		mfem::Mesh mesh = mfem::Mesh::MakeCartesian2D(
			6, 8, mfem::Element::TRIANGLE, false, rMax - rMin, zMax - zMin );
		mesh.Transform( [ & ]( mfem::Vector const &in, mfem::Vector &out )
		{
			out = in;
			out( 0 ) += rMin;
			out( 1 ) += zMin;
		} );
		std::ofstream file( "driver-acceptance-frommesh.mesh" );
		mesh.Print( file );
	}

	{
		std::ofstream file( "driver-acceptance-frommesh.toml" );
		file << "[mesh]\n"
		        "File = \"driver-acceptance-frommesh.mesh\"\n"
		        "\n[discretisation]\nPolynomialDegree = 2\n"
		        "\n[source]\nType = \"soloviev\"\nA = -0.52\n"
		        "\n[boundary]\nType = \"zero\"\n"
		        "\n[output]\nDirectory = \".\"\nPrefix = \"frommesh\"\n"
		        "GridNR = 33\nGridNZ = 33\n";
	}

	BOOST_TEST_REQUIRE( run( "driver-acceptance-frommesh.toml" ) == 0,
	                    "the driver did not exit 0 on a mesh read from a file" );
	BOOST_TEST_REQUIRE( exists( "frommesh.nc" ),
	                    "the gridded output was not written, which is the "
	                    "failure this case exists for -- the solve succeeded "
	                    "and the answer was lost at the output stage" );

	// The extent is the MESH's, read back through ncdump rather than through
	// MEQ's own writer.
	std::pair<double, double> const r = coordinateRange( "frommesh.nc", "R" );
	std::pair<double, double> const z = coordinateRange( "frommesh.nc", "Z" );
	BOOST_TEST( r.first == rMin, boost::test_tools::tolerance( 1.0e-12 ) );
	BOOST_TEST( r.second == rMax, boost::test_tools::tolerance( 1.0e-12 ) );
	BOOST_TEST( z.first == zMin, boost::test_tools::tolerance( 1.0e-12 ) );
	BOOST_TEST( z.second == zMax, boost::test_tools::tolerance( 1.0e-12 ) );
	std::printf( "\n  a mesh MEQ did not build\n"
	             "    grid R [%.4f, %.4f]  Z [%.4f, %.4f]  from the mesh itself\n",
	             r.first, r.second, z.first, z.second );

	for ( char const *path : { "driver-acceptance-frommesh.toml",
	                           "driver-acceptance-frommesh.mesh" } )
		std::remove( path );
}

BOOST_AUTO_TEST_CASE( theDriverReportsConfigurationErrorsAsExitOne )
{
	BOOST_TEST( run( "examples/does-not-exist.toml" ) == 1,
	            "a missing configuration file must exit 1" );

	// A file that parses as TOML but is not a valid meq configuration.
	{
		std::ofstream bad( "driver-acceptance-bad.toml" );
		// Complete but for the one thing under test, so that exit 1 is earned by
		// the unknown source type and not by an unrelated missing table.
		bad << "[mesh]\nRMin = 0.1\nRMax = 1.9\nZMin = -1.7\nZMax = 1.7\n"
		    << "NR = 3\nNZ = 4\n\n[discretisation]\nPolynomialDegree = 2\n"
		    << "\n[source]\nType = \"not-a-source\"\n";
	}
	BOOST_TEST( run( "driver-acceptance-bad.toml" ) == 1,
	            "an unknown [source] Type must exit 1" );

	// An unknown key, which the parser rejects rather than ignores -- the
	// reason rejectUnknownKeys() exists is that a typo in a key name otherwise
	// becomes a silently defaulted value.
	{
		std::ofstream bad( "driver-acceptance-typo.toml" );
		bad << "[mesh]\nRMin = 0.1\nRMax = 1.9\nZMin = -1.7\nZMax = 1.7\n"
		    << "NR = 3\nNZ = 4\n\n[discretisation]\nPolynomialDegree = 2\n"
		    << "\n[source]\nType = \"soloviev\"\nA = -0.52\n"
		    << "\n[output]\nPrefix = \"typo\"\nGirdNR = 65\n";
	}
	BOOST_TEST( run( "driver-acceptance-typo.toml" ) == 1,
	            "a misspelled key must exit 1 rather than be silently ignored" );

	std::remove( "driver-acceptance-bad.toml" );
	std::remove( "driver-acceptance-typo.toml" );
}

/*
 * THE TWO PERFORMANCE KEYS, AND WHAT THE DRIVER REFUSES RATHER THAN
 * APPROXIMATES.
 *
 * `[solver] AssemblyMode` and `TraceSolver` are the only solver knobs exposed to
 * TOML, and they are exposed because neither can change the answer -- the two
 * assembly modes are asserted bit for bit and the three trace solvers agree to
 * 1e-14. That is also exactly why an unavailable one must be REFUSED rather than
 * substituted: a silent fallback would be invisible in the result.
 *
 * THIS CASE EXISTS BECAUSE EXERCISING THE KEYS FOUND A REAL BUG THAT REASONING
 * ABOUT THEM DID NOT. `TraceSolver = "cudss"` passed every check -- the spelling
 * is valid and this build genuinely has cuDSS -- and then aborted inside CUDA
 * with `cudaMemcpyDeviceToDevice ... invalid argument`, because cuDSS reads its
 * data through MFEM's device-aware accessors and the driver configures no
 * mfem::Device. The library had documented that requirement all along; putting
 * the choice in a config file is what made it reachable by somebody who had not
 * read it. The driver now refuses it with a message naming the reason, and this
 * is what stops that regressing.
 */
BOOST_AUTO_TEST_CASE( theDriverRefusesASolverItCannotHonour )
{
	auto write = []( char const *path, char const *solverBody )
	{
		std::ofstream file( path );
		file << "[mesh]\nRMin = 0.1\nRMax = 1.9\nZMin = -1.7\nZMax = 1.7\n"
		     << "NR = 4\nNZ = 4\n\n[discretisation]\nPolynomialDegree = 1\n"
		     << "\n[source]\nType = \"soloviev\"\nA = -0.52\n"
		     << "\n[solver]\n" << solverBody
		     << "\n[output]\nPrefix = \"driver-acceptance-solver\"\n";
	};

	// A misspelling, refused at parse -- Config can check the string and does.
	write( "driver-acceptance-mode.toml", "AssemblyMode = \"parallel\"\n" );
	BOOST_TEST( run( "driver-acceptance-mode.toml" ) == 1,
	            "an unknown AssemblyMode spelling must exit 1" );

	// The capitalised spelling, which is the plausible mistake: TOML KEYS in
	// meq are UpperCamelCase, so a reader may expect the VALUES to be too.
	// There is no case folding anywhere in Config, deliberately.
	write( "driver-acceptance-case.toml", "AssemblyMode = \"Threaded\"\n" );
	BOOST_TEST( run( "driver-acceptance-case.toml" ) == 1,
	            "AssemblyMode values are lower case and are compared literally, "
	            "so \"Threaded\" must be refused rather than quietly accepted" );

	// cuDSS: a valid spelling, possibly present in the build, and refused by the
	// DRIVER because it needs an mfem::Device that the driver does not create.
	// Asserted unconditionally, because the refusal does not depend on whether
	// this build has cuDSS -- that is the point of it.
	write( "driver-acceptance-cudss.toml", "TraceSolver = \"cudss\"\n" );
	BOOST_TEST( run( "driver-acceptance-cudss.toml" ) == 1,
	            "TraceSolver = \"cudss\" must exit 1 from the driver. If this "
	            "passes, the driver has started accepting it -- and without an "
	            "mfem::Device that is not a solve, it is a raw CUDA abort with "
	            "nothing in the message about the key that caused it" );

	// And the defaults, plus an explicit umfpack, must still SOLVE. A test that
	// only checked refusals would pass with every value refused.
	write( "driver-acceptance-ok.toml", "TraceSolver = \"umfpack\"\n" );
	BOOST_TEST( run( "driver-acceptance-ok.toml" ) == 0,
	            "an explicit TraceSolver = \"umfpack\" must solve, or the "
	            "refusals above are testing nothing" );

	std::remove( "driver-acceptance-mode.toml" );
	std::remove( "driver-acceptance-case.toml" );
	std::remove( "driver-acceptance-cudss.toml" );
	std::remove( "driver-acceptance-ok.toml" );
}

/*
 * --help and --version are the two things a user tries first, and a binary that
 * exits non-zero on --help is a binary that looks broken.
 */
BOOST_AUTO_TEST_CASE( theDriverHasAUsableCommandLine )
{
	BOOST_TEST( run( "--help" ) == 0 );
	BOOST_TEST( run( "--version" ) == 0 );

	// No arguments is a usage error, not success: a scan script that lost its
	// argument must not be told everything went fine.
	std::string const command = std::string( driver() ) + " > /dev/null 2>&1";
	int const status = std::system( command.c_str() );
	BOOST_TEST( WEXITSTATUS( status ) == 1,
	            "running with no arguments must exit 1" );
}
