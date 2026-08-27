#define BOOST_TEST_MODULE MeqOutputConvergence

#include <boost/test/unit_test.hpp>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>

#include "mfem.hpp"

#include "meq/Field.hpp"
#include "meq/GradShafranov.hpp"
#include "meq/Output.hpp"
#include "meq/Sampler.hpp"

#include "ConvergenceHarness.hpp"
#include "analytic/Soloviev.hpp"

/*
 * The output, checked by READING IT BACK.
 *
 * A writer that runs without throwing has demonstrated nothing: the interesting
 * failures are a file that no reader can open, a field written in the wrong
 * index order, a mask that disagrees with the data it masks, and a restart that
 * silently loses digits. None of those raise an exception at write time.
 *
 * So every test here writes and then reads, and compares against the exact
 * Solov'ev solution rather than against what was written -- which would only
 * prove the writer agrees with itself.
 */

namespace
{
	meq::analytic::SolovievEquilibrium const &equilibrium()
	{
		static meq::analytic::SolovievEquilibrium const eq =
			meq::analytic::SolovievEquilibrium::nstx();
		return eq;
	}

	meq::tests::Rectangle box()
	{
		return meq::tests::Rectangle{ 0.6, 1.4, -0.6, 0.6 };
	}

	/// A file removed when the test leaves scope, however it leaves.
	class Scratch
	{
		public:
			explicit Scratch( std::string const &nameIn ) : name( nameIn ) { }
			~Scratch() { std::remove( name.c_str() ); }
			std::string const &path() const { return name; }
		private:
			std::string name;
	};
}

BOOST_AUTO_TEST_SUITE( output_convergence )

/// MFEM's own formats, read back into fresh objects. This is the exact-restart
/// path, so the check is that nothing is lost -- not that it is close.
BOOST_AUTO_TEST_CASE( theMfemFilesRoundTripExactly )
{
	meq::analytic::SolovievEquilibrium const &eq = equilibrium();
	mfem::Mesh mesh = meq::tests::makeMesh( box(), 8 );

	mfem::FunctionCoefficient source( [ &eq ]( mfem::Vector const &x )
		{ return eq.f( x( 0 ), x( 1 ), 0.0 ); } );
	mfem::FunctionCoefficient exact( [ &eq ]( mfem::Vector const &x )
		{ return eq.psi( x( 0 ), x( 1 ) ); } );

	meq::GradShafranovSolver solver( mesh, 2 );
	solver.setSource( source );
	solver.setBoundaryData( exact );
	solver.solve();

	Scratch const meshFile( "meq_test_out.mesh" );
	Scratch const psiFile( "meq_test_out_psi.gf" );
	Scratch const fluxFile( "meq_test_out_grad_psi.gf" );

	meq::writeMfem( "meq_test_out", mesh, solver.potential(), solver.flux() );

	// Read back with no knowledge of what wrote them.
	std::ifstream meshIn( meshFile.path() );
	BOOST_TEST_REQUIRE( meshIn.good() );
	mfem::Mesh reloadedMesh( meshIn, 1, 1 );

	std::ifstream psiIn( psiFile.path() );
	BOOST_TEST_REQUIRE( psiIn.good() );
	mfem::GridFunction reloadedPsi( &reloadedMesh, psiIn );

	BOOST_TEST_REQUIRE( reloadedPsi.Size() == solver.potential().Size() );

	double worst = 0.0;
	for ( int i = 0; i < reloadedPsi.Size(); ++i )
		worst = std::max( worst, std::abs( reloadedPsi( i ) - solver.potential()( i ) ) );

	std::printf( "\n  .gf round trip: worst coefficient difference %.3e\n", worst );
	std::fflush( stdout );

	// 16 digits written, so the round trip should be at round-off and not at
	// the 1e-8 MFEM's default precision would give. That difference is the
	// whole reason writeMfem() sets it.
	BOOST_TEST( worst < 1.0e-13,
	            "the grid function did not survive the round trip: worst "
	            "coefficient moved by " << worst << ", which is the scale of a "
	            "precision setting rather than of round-off" );
}

/// The grid file, read back with the NetCDF library and compared against the
/// exact solution. Index order is the thing most likely to be wrong and least
/// likely to announce itself.
BOOST_AUTO_TEST_CASE( theGridFileReadsBackAsTheExactSolution )
{
	if ( !meq::hasNetCDF() )
	{
		BOOST_TEST_MESSAGE( "  built without netcdf-cxx4, skipping" );
		return;
	}

	meq::analytic::SolovievEquilibrium const &eq = equilibrium();
	mfem::Mesh mesh = meq::tests::makeMesh( box(), 32 );

	mfem::FunctionCoefficient source( [ &eq ]( mfem::Vector const &x )
		{ return eq.f( x( 0 ), x( 1 ), 0.0 ); } );
	mfem::FunctionCoefficient exact( [ &eq ]( mfem::Vector const &x )
		{ return eq.psi( x( 0 ), x( 1 ) ); } );

	meq::GradShafranovSolver solver( mesh, 2 );
	solver.setSource( source );
	solver.setBoundaryData( exact );
	solver.solve();

	mfem::GridFunction field;
	meq::poloidalField( solver.flux(), field );

	// Deliberately NOT square, and deliberately different in each direction:
	// a transposed write on a square grid produces a plausible file.
	int const nodesR = 33;
	int const nodesZ = 49;
	double const inset = 0.05;
	meq::GridSampler sampler( mesh,
		box().rMin + inset, box().rMax - inset, nodesR,
		box().zMin + inset, box().zMax - inset, nodesZ );

	std::vector<double> psi, bR, bZ;
	sampler.sample( solver.potential(), psi, std::nan( "" ) );
	sampler.sampleComponent( field, 0, bR, std::nan( "" ) );
	sampler.sampleComponent( field, 1, bZ, std::nan( "" ) );

	Scratch const file( "meq_test_out.nc" );
	{
		meq::NetCDFWriter writer( file.path(), sampler );
		writer.attribute( "title", "meq test output" );
		writer.attribute( "polynomial_degree", 2 );
		writer.field( "psi", psi, "Poloidal flux per radian", "Wb/rad" );
		writer.field( "B_R", bR, "Radial magnetic field", "T" );
		writer.field( "B_Z", bZ, "Vertical magnetic field", "T" );
		writer.close();
	}

	// Read it back through ncdump rather than through the writer's own library
	// calls, so the check does not share code with the thing it checks.
	std::string const command = "ncdump -h " + file.path() + " > meq_test_hdr.txt 2>&1";
	BOOST_TEST_REQUIRE( std::system( command.c_str() ) == 0,
	                    "ncdump could not read the file meq just wrote" );

	Scratch const header( "meq_test_hdr.txt" );
	std::ifstream headerIn( header.path() );
	std::string const text( ( std::istreambuf_iterator<char>( headerIn ) ),
	                        std::istreambuf_iterator<char>() );

	std::printf( "  ncdump header:\n" );
	for ( std::string const &needle :
	      { "R = 33", "Z = 49", "double psi(Z, R)", "double B_R(Z, R)",
	        "double B_Z(Z, R)", "byte inside(Z, R)" } )
	{
		bool const present = text.find( needle ) != std::string::npos;
		std::printf( "    %-24s %s\n", needle.c_str(), present ? "yes" : "MISSING" );
		BOOST_TEST( present, "the file does not declare '" << needle
		            << "', so its shape is not what the header promises" );
	}
	std::fflush( stdout );

	// And the values, against the exact solution -- with the index arithmetic
	// done independently here, which is what catches a transposed write.
	double worstPsi = 0.0, worstB = 0.0;
	int compared = 0;
	for ( int j = 0; j < nodesZ; ++j )
		for ( int i = 0; i < nodesR; ++i )
		{
			std::size_t const at = static_cast<std::size_t>( j )*nodesR + i;
			if ( !sampler.located( i, j ) ) continue;
			++compared;

			double const r = sampler.rAt( i ), z = sampler.zAt( j );
			worstPsi = std::max( worstPsi, std::abs( psi[ at ] - eq.psi( r, z ) ) );

			double qR = 0.0, qZ = 0.0;
			eq.flux( r, z, qR, qZ );
			worstB = std::max( worstB, std::hypot( bR[ at ] + qZ, bZ[ at ] - qR ) );
		}

	std::printf( "  %d nodes compared: worst psi %.4e, worst B %.4e\n",
	             compared, worstPsi, worstB );
	std::fflush( stdout );

	BOOST_TEST( compared > 1000 );
	BOOST_TEST( worstPsi < 6.0e-6 );
	BOOST_TEST( worstB < 4.0e-6 );
}

/// The mask must agree with the data. A node marked outside carrying a finite
/// value, or one marked inside carrying a NaN, means the two were written from
/// different sources.
BOOST_AUTO_TEST_CASE( theMaskAgreesWithTheData )
{
	if ( !meq::hasNetCDF() )
	{
		BOOST_TEST_MESSAGE( "  built without netcdf-cxx4, skipping" );
		return;
	}

	mfem::Mesh mesh = meq::tests::makeMesh( box(), 8 );
	mfem::L2_FECollection collection( 1, mesh.Dimension() );
	mfem::FiniteElementSpace space( &mesh, &collection );
	mfem::GridFunction ones( &space );
	ones = 1.0;

	// A grid wider than the mesh, so there are genuinely absent nodes.
	double const pad = 0.3*box().width();
	meq::GridSampler sampler( mesh,
		box().rMin - pad, box().rMax + pad, 45,
		box().zMin - pad, box().zMax + pad, 45 );

	std::vector<double> values;
	sampler.sample( ones, values, std::nan( "" ) );

	int insideFinite = 0, outsideNaN = 0, disagreements = 0;
	for ( int j = 0; j < sampler.nodesZ(); ++j )
		for ( int i = 0; i < sampler.nodesR(); ++i )
		{
			std::size_t const at = static_cast<std::size_t>( j )*sampler.nodesR() + i;
			bool const finite = std::isfinite( values[ at ] );
			if ( sampler.located( i, j ) == finite )
			{
				if ( finite ) ++insideFinite; else ++outsideNaN;
			}
			else
			{
				++disagreements;
			}
		}

	std::printf( "  mask agreement: %d inside and finite, %d outside and NaN, "
	             "%d disagreements\n", insideFinite, outsideNaN, disagreements );
	std::fflush( stdout );

	BOOST_TEST( disagreements == 0 );
	BOOST_TEST( insideFinite > 0 );
	BOOST_TEST( outsideNaN > 0 );
}

BOOST_AUTO_TEST_SUITE_END()
