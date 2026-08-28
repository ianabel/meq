#define BOOST_TEST_MODULE DriverAcceptance
#include <boost/test/unit_test.hpp>

#include "mfem.hpp"

#include "meq/BoundaryShape.hpp"
#include "meq/GradShafranov.hpp"
#include "meq/Source.hpp"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <string>

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
		// MKL_THREADING_LAYER is inherited from the ctest environment, which
		// CMake sets for exactly the reason the other convergence tests record.
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
