#define BOOST_TEST_MODULE MeqOutputConvergence

#include <boost/test/unit_test.hpp>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iterator>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include "mfem.hpp"

#include "meq/Field.hpp"
#include "meq/GradShafranov.hpp"
#include "meq/Output.hpp"
#include "meq/Profiles.hpp"
#include "meq/RotatingSource.hpp"
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

	/// n( psi ) = a + b psi, for the rotating case below. It needs a genuine
	/// psi-dependence: with a constant density the field depends on r alone, and
	/// the test would say nothing about pairing the node's R with the node's psi
	/// -- only about the R.
	class LinearProfile : public meq::Profile
	{
		public:
			LinearProfile( double aIn, double bIn ) : a( aIn ), b( bIn ) { }
			double operator()( double psi ) const override { return a + b*psi; }
			double prime( double ) const override { return b; }
			double doublePrime( double ) const override { return 0.0; }
		private:
			double a, b;
	};

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

/// CONTINUING THE POTENTIAL ACROSS THE Gamma_h-TO-Gamma BAND, and the rate
/// that says it is a Taylor extension rather than an extrapolation.
///
/// THE OBVIOUS TEST DOES NOT WORK, AND THE REASON IS WORTH RECORDING. For a
/// LINEAR psi the extension psi( x0 ) + grad psi( x0 ) . ( p - x0 ) is exactly
/// psi, so the error ought to be round-off -- except that meq's flux is
/// q = ( 1/r ) grad psi, and for a linear psi that is ( alpha/r, beta/r ),
/// which is NOT a polynomial. Projecting it into the flux space costs
/// O( h^(k+1) ) and the band inherits it: measured 1.2e-06 where the interior
/// nodes were exact at 2.2e-15. No psi makes both psi and q exactly
/// representable at once -- q polynomial needs grad psi proportional to r, and
/// the extension is only exact for psi linear.
///
/// So the assertion is on the RATE. The band error must fall at the flux's own
/// order, which is what distinguishes a Taylor extension carried by an
/// accurately computed derivative from a polynomial continued outside its
/// element -- the latter is bounded by nothing and does not converge in the
/// band at all.
BOOST_AUTO_TEST_CASE( theBandIsContinuedByTheFluxAtTheFluxesOwnOrder )
{
	double const alpha = 0.7, beta = -1.3;
	int const degree = 2;
	meq::tests::Rectangle const inner{ 0.8, 1.2, -0.2, 0.2 };

	auto measure = [ & ]( int cells, double &interiorWorst )
	{
		mfem::Mesh mesh = meq::tests::makeMesh( inner, cells );
		mfem::L2_FECollection collection( degree, mesh.Dimension() );
		mfem::FiniteElementSpace scalar( &mesh, &collection );
		mfem::FiniteElementSpace vector( &mesh, &collection, 2 );

		mfem::FunctionCoefficient exact( [ = ]( mfem::Vector const &x )
			{ return alpha*x( 0 ) + beta*x( 1 ); } );
		mfem::VectorFunctionCoefficient fluxOf( 2,
			[ = ]( mfem::Vector const &x, mfem::Vector &q )
			{
				q( 0 ) = alpha/x( 0 );
				q( 1 ) = beta/x( 0 );
			} );

		mfem::GridFunction potential( &scalar );
		mfem::GridFunction flux( &vector );
		potential.ProjectCoefficient( exact );
		flux.ProjectCoefficient( fluxOf );

		meq::GridSampler sampler( mesh, 0.7, 1.3, 41, -0.3, 0.3, 41 );
		int const filled = sampler.extendOutward( 1.0 );
		BOOST_TEST_REQUIRE( filled > 0, "no node was continued at " << cells );

		std::vector<double> values;
		sampler.samplePotentialWithFlux( potential, flux, values,
			std::numeric_limits<double>::quiet_NaN() );

		double band = 0.0;
		interiorWorst = 0.0;
		for ( int j = 0; j < sampler.nodesZ(); ++j )
			for ( int i = 0; i < sampler.nodesR(); ++i )
			{
				if ( !sampler.located( i, j ) )
					continue;
				double const r = sampler.rAt( i ), z = sampler.zAt( j );
				double const got =
					values[ static_cast<std::size_t>( j )*sampler.nodesR() + i ];
				double const error = std::abs( got - ( alpha*r + beta*z ) );
				bool const outside = r < inner.rMin || r > inner.rMax
				                     || z < inner.zMin || z > inner.zMax;
				( outside ? band : interiorWorst ) =
					std::max( outside ? band : interiorWorst, error );
			}
		return band;
	};

	double coarseInterior = 0.0, fineInterior = 0.0;
	double const coarse = measure( 8, coarseInterior );
	double const fine = measure( 16, fineInterior );
	double const rate = std::log2( coarse/std::max( 1.0e-300, fine ) );

	std::printf( "\n  band continuation: %.3e at n = 8, %.3e at n = 16, "
	             "rate %.2f (flux order is k+1 = %d)\n"
	             "                     interior nodes exact to %.3e\n",
	             coarse, fine, rate, degree + 1,
	             std::max( coarseInterior, fineInterior ) );
	std::fflush( stdout );

	// Interior nodes involve no continuation at all, so they ARE exact.
	BOOST_TEST( std::max( coarseInterior, fineInterior ) < 1.0e-12,
	            "an interior node is out by "
	            << std::max( coarseInterior, fineInterior )
	            << " on a linear potential, which needs no continuation" );

	// The floor is 2.0 rather than k+1 = 3 because the band error mixes the
	// flux's projection error with the O( |p - x0|^2 ) of the extension itself,
	// and the second term does not care how good q is. What it rules out is the
	// thing that matters: a continuation that does not converge.
	BOOST_TEST( rate > 2.0,
	            "the band error converges at " << rate << ", which is too slow "
	            "to be a Taylor extension carried by the flux. A polynomial "
	            "continued outside its own element behaves like this" );
}

/// `B` USED TO BE ZEROTH ORDER IN THE BAND. THIS IS THE TEST THAT SAYS IT IS
/// NOT ANY MORE, AND IT MEASURES BOTH SIDES SO THE CLAIM IS A COMPARISON.
///
/// Until 2026-09-02 the driver sent `B` through sampleComponent(), which reads
/// the FOOT on Gamma_h -- so about one node in ten of the interchange file
/// carried a piecewise-constant field, `O( h )`, while `psi` beside it got the
/// flux's Taylor step and `q` itself resolves at `k+1`. Nothing could see it:
/// the `inside` mask says 1 in the band, and the only band test in this file
/// tested `psi`.
///
/// THE FIELD IS A QUADRATIC AND P2 HOLDS IT EXACTLY, which is the design of the
/// test rather than a convenience. With no projection error anywhere, an
/// interior node must be exact to round-off and the band error is the Taylor
/// truncation ALONE -- so the rate measures the continuation and nothing else,
/// unlike the `psi` test above whose 3.75 mixes in the flux's own projection
/// error.
///
/// AND ITS GRADIENT IS DELIBERATELY ASYMMETRIC -- `d u0/dz = 2*gamma*z` against
/// `d u1/dr = beta*z` -- because MFEM's `grad( i, j )` is component-first and
/// reading it transposed is otherwise a silent error that a symmetric test
/// field would pass.
BOOST_AUTO_TEST_CASE( theBandVectorContinuesAtItsGradientsOrder )
{
	double const alpha = 0.7, gamma = 0.5, beta = -1.3;
	int const degree = 2;
	meq::tests::Rectangle const inner{ 0.8, 1.2, -0.2, 0.2 };

	auto exactly = [ = ]( int component, double r, double z )
	{
		return component == 0 ? alpha*r*r + gamma*z*z : beta*r*z;
	};

	auto measure = [ & ]( int cells, double &viaFoot, double &interiorWorst )
	{
		mfem::Mesh mesh = meq::tests::makeMesh( inner, cells );
		mfem::L2_FECollection collection( degree, mesh.Dimension() );
		mfem::FiniteElementSpace vector( &mesh, &collection, 2 );

		mfem::VectorFunctionCoefficient fieldOf( 2,
			[ = ]( mfem::Vector const &x, mfem::Vector &u )
			{
				u( 0 ) = alpha*x( 0 )*x( 0 ) + gamma*x( 1 )*x( 1 );
				u( 1 ) = beta*x( 0 )*x( 1 );
			} );

		mfem::GridFunction field( &vector );
		field.ProjectCoefficient( fieldOf );

		meq::GridSampler sampler( mesh, 0.7, 1.3, 41, -0.3, 0.3, 41 );
		int const filled = sampler.extendOutward( 1.0 );
		BOOST_TEST_REQUIRE( filled > 0, "no node was continued at " << cells );

		double const nan = std::numeric_limits<double>::quiet_NaN();
		std::vector<double> withGradient, atFoot;

		double band = 0.0;
		viaFoot = 0.0;
		interiorWorst = 0.0;
		for ( int component = 0; component < 2; ++component )
		{
			sampler.sampleComponentWithGradient( field, component,
			                                     withGradient, nan );
			sampler.sampleComponent( field, component, atFoot, nan );

			for ( int j = 0; j < sampler.nodesZ(); ++j )
				for ( int i = 0; i < sampler.nodesR(); ++i )
				{
					if ( !sampler.located( i, j ) )
						continue;

					std::size_t const at =
						static_cast<std::size_t>( j )*sampler.nodesR() + i;
					double const want =
						exactly( component, sampler.rAt( i ), sampler.zAt( j ) );

					if ( sampler.wasExtended( i, j ) )
					{
						band = std::max( band,
						                 std::abs( withGradient[ at ] - want ) );
						viaFoot = std::max( viaFoot,
						                    std::abs( atFoot[ at ] - want ) );
					}
					else
						interiorWorst =
							std::max( interiorWorst,
							          std::abs( withGradient[ at ] - want ) );
				}
		}
		return band;
	};

	double coarseFoot = 0.0, fineFoot = 0.0;
	double coarseInterior = 0.0, fineInterior = 0.0;
	double const coarse = measure( 8, coarseFoot, coarseInterior );
	double const fine = measure( 16, fineFoot, fineInterior );

	double const tiny = 1.0e-300;
	double const rate = std::log2( coarse/std::max( tiny, fine ) );
	double const footRate = std::log2( coarseFoot/std::max( tiny, fineFoot ) );
	double const gain = fineFoot/std::max( tiny, fine );

	std::printf( "\n  band vector continuation, k = %d:\n"
	             "    with grad u:   %.3e at n = 8, %.3e at n = 16, rate %.2f\n"
	             "    at the foot:   %.3e at n = 8, %.3e at n = 16, rate %.2f\n"
	             "    gain at n = 16: %.1fx;  interior nodes exact to %.3e\n",
	             degree, coarse, fine, rate,
	             coarseFoot, fineFoot, footRate, gain,
	             std::max( coarseInterior, fineInterior ) );
	std::fflush( stdout );

	// The field is in the space, so an interior node has no discretisation error
	// to show. Anything above round-off here is the sampler, not the extension.
	BOOST_TEST( std::max( coarseInterior, fineInterior ) < 1.0e-11,
	            "an interior node is out by "
	            << std::max( coarseInterior, fineInterior )
	            << " on a field its own space represents exactly" );

	// O( h^2 ), and 2 is the whole truncation error here rather than a floor
	// mixed with something else.
	BOOST_TEST( rate > 1.7,
	            "the band converges at " << rate << ", not the 2 a first-order "
	            "Taylor step gives. Is grad u being applied at all?" );

	// THE CONTROL, and it is the assertion that would have caught the original
	// defect: reading the foot is a zeroth-order continuation and must converge
	// at about 1. If this ever reads 2 as well, the two paths have silently
	// become the same path.
	BOOST_TEST( footRate < 1.4,
	            "reading the foot converges at " << footRate << ", which is too "
	            "good for a zeroth-order continuation -- the two sampling paths "
	            "may have become the same one" );

	BOOST_TEST( gain > 5.0,
	            "the gradient continuation is only " << gain << " times better "
	            "than reading the foot, which does not justify its cost" );
}

/// THE ADAPTIVE SERIES, AND THE ONE ASSERTION THAT CATCHES THE FAILURE THAT
/// ACTUALLY HAPPENED.
///
/// Writing a frame per cycle is easy to get almost right: rebuild the
/// collection each time and every Cycle directory appears, holding the correct
/// refined mesh, while the .pvd index lists only the LAST of them. ParaView
/// then opens the file and shows a single frame. Nothing errors, the data is
/// all on disk, and the animation is silently missing.
///
/// So the assertion is on the INDEX -- one DataSet entry per append -- rather
/// than on the pieces, which were never the thing that broke.
BOOST_AUTO_TEST_CASE( theAdaptiveSeriesIndexesEveryFrame )
{
	mfem::Mesh mesh = meq::tests::makeMesh( box(), 4 );

	std::string const stem = "meq_test_series";
	int const degree = 1;
	int const frames = 3;
	{
		meq::VtuSeries series( stem, degree );
		for ( int cycle = 0; cycle < frames; ++cycle )
		{
			// A different mesh every frame, which is the situation the loop
			// creates and the reason the collection cannot simply be reused.
			if ( cycle > 0 )
				mesh.UniformRefinement();

			mfem::L2_FECollection collection( degree, mesh.Dimension() );
			mfem::FiniteElementSpace scalar( &mesh, &collection );
			mfem::FiniteElementSpace vector( &mesh, &collection, 2 );
			mfem::GridFunction potential( &scalar );
			mfem::GridFunction field( &vector );
			potential = static_cast<double>( cycle );
			field = 0.0;

			series.append( mesh, potential, field, cycle,
			               static_cast<double>( cycle ) );
		}
		BOOST_TEST( series.frames() == frames );
	}

	std::string const pvd = stem + "_cycles/" + stem + "_cycles.pvd";
	std::ifstream index( pvd );
	BOOST_TEST_REQUIRE( index.good(), "no series index at " << pvd );
	std::string const content( ( std::istreambuf_iterator<char>( index ) ),
	                             std::istreambuf_iterator<char>() );

	int entries = 0;
	for ( std::size_t at = content.find( "<DataSet" );
	      at != std::string::npos;
	      at = content.find( "<DataSet", at + 1 ) )
		++entries;

	std::printf( "\n  series index: %d DataSet entries for %d appends\n",
	             entries, frames );
	std::fflush( stdout );

	BOOST_TEST( entries == frames,
	            "the .pvd indexes " << entries << " frames for " << frames
	            << " appends. Every Cycle directory is probably on disk and "
	            "correct -- it is the INDEX that is short, and ParaView shows "
	            "only what the index lists" );

	for ( int cycle = 0; cycle < frames; ++cycle )
	{
		char directory[ 64 ];
		std::snprintf( directory, sizeof directory, "%s_cycles/Cycle%06d",
		               stem.c_str(), cycle );
		std::string const piece = std::string( directory ) + "/proc000000.vtu";
		std::ifstream vtu( piece );
		BOOST_TEST( vtu.good(), "no piece for frame " << cycle << " at " << piece );
		vtu.close();
		std::remove( piece.c_str() );
		std::remove( ( std::string( directory ) + "/data.pvtu" ).c_str() );
		std::remove( directory );
	}
	index.close();
	std::remove( pvd.c_str() );
	std::remove( ( stem + "_cycles" ).c_str() );
}

/// BENDING Gamma_h OUT ONTO Gamma, which is what makes a curved-path picture
/// show the boundary that was asked for rather than the polygon inscribed in
/// it.
///
/// The property is exact and so is the assertion: after curveBoundaryOnto with
/// a circle projector, EVERY boundary node must lie on the circle. Not close to
/// it -- on it, because the projector puts it there and nothing afterwards
/// moves it. A test at a tolerance would pass while a fraction of the
/// displacement was quietly being backed off for tangling, which is the failure
/// this is for.
BOOST_AUTO_TEST_CASE( theBoundaryBendsOntoTheTrueGamma )
{
	// A background box, and the subdomain of it inside a circle -- the same
	// construction the driver's curved path uses, so Gamma_h here is a real
	// inscribed polygon and not a contrivance.
	double const centreR = 1.0, centreZ = 0.0, radius = 0.35;
	auto const circle = [ = ]( mfem::Vector const &x )
	{
		return std::hypot( x( 0 ) - centreR, x( 1 ) - centreZ ) - radius;
	};

	mfem::Mesh background = meq::tests::makeMesh(
		meq::tests::Rectangle{ 0.5, 1.5, -0.5, 0.5 }, 16 );

	mfem::Array<int> marker;
	int const inside = mfem::MarkLevelSetSubdomain( background, circle, 0.0,
	                                               marker, 1 );
	BOOST_TEST_REQUIRE( inside > 0 );
	for ( int e = 0; e < background.GetNE(); ++e )
		background.SetAttribute( e, marker[ e ] ? 1 : 2 );
	background.SetAttributes();

	mfem::Array<int> domainAttribute( 1 );
	domainAttribute[ 0 ] = 1;
	mfem::SubMesh mesh = mfem::SubMesh::CreateFromDomain( background,
	                                                      domainAttribute );

	// Gamma_h is inscribed, so before bending every BOUNDARY vertex is strictly
	// inside the circle. That gap is the thing being fixed, so it is worth
	// pinning -- and it has to be measured over the boundary vertices alone.
	// Taken over every vertex it would report the one nearest the centre, which
	// is the radius itself and says nothing about Gamma_h.
	double worstBefore = 0.0;
	for ( int b = 0; b < mesh.GetNBE(); ++b )
	{
		mfem::Array<int> vertices;
		mesh.GetBdrElementVertices( b, vertices );
		for ( int n = 0; n < vertices.Size(); ++n )
		{
			double const *p = mesh.GetVertex( vertices[ n ] );
			worstBefore = std::max( worstBefore,
				radius - std::hypot( p[ 0 ] - centreR, p[ 1 ] - centreZ ) );
		}
	}
	BOOST_TEST_REQUIRE( worstBefore > 0.0,
	                    "the subdomain already reaches the circle, so this "
	                    "measures nothing" );

	double applied = 0.0;
	int const moved = meq::curveBoundaryOnto( mesh, 2,
		[ = ]( double r, double z, double &outR, double &outZ )
		{
			double const vR = r - centreR, vZ = z - centreZ;
			double const rho = std::hypot( vR, vZ );
			outR = r;
			outZ = z;
			if ( rho > 0.0 )
			{
				outR = centreR + vR*radius/rho;
				outZ = centreZ + vZ*radius/rho;
			}
		}, applied );

	BOOST_TEST_REQUIRE( moved > 0, "no boundary nodes were moved at all" );
	BOOST_TEST( applied == 1.0,
	            "only " << applied << " of the displacement was applied, so an "
	            "element tangled. The boundary is then somewhere between "
	            "Gamma_h and Gamma, which is worse than either" );

	// Every node of every boundary face, in the CURVED nodal space -- vertices
	// and the high-order nodes between them alike. A curvature that moved only
	// the vertices would leave the edges straight and fail here.
	mfem::GridFunction *nodes = mesh.GetNodes();
	BOOST_TEST_REQUIRE( nodes != nullptr );
	mfem::FiniteElementSpace const *space = nodes->FESpace();

	double worst = 0.0;
	int checked = 0;
	for ( int b = 0; b < mesh.GetNBE(); ++b )
	{
		mfem::Array<int> vdofs;
		space->GetBdrElementVDofs( b, vdofs );
		int const perComponent = vdofs.Size()/2;
		for ( int n = 0; n < perComponent; ++n )
		{
			double const r = ( *nodes )( vdofs[ n ] );
			double const z = ( *nodes )( vdofs[ perComponent + n ] );
			worst = std::max( worst,
				std::abs( std::hypot( r - centreR, z - centreZ ) - radius ) );
			++checked;
		}
	}

	std::printf( "\n  boundary bent onto Gamma: %d nodes over %d faces, "
	             "was up to %.3e inside, now off the circle by %.3e\n",
	             checked, mesh.GetNBE(), worstBefore, worst );
	std::fflush( stdout );

	BOOST_TEST( worst < 1.0e-12,
	            "a boundary node sits " << worst << " off the circle after "
	            "projection onto it. Either a node was missed -- the high-order "
	            "ones between the vertices are the likely ones -- or the "
	            "displacement was scaled back" );
}

/// THE VTK FILES, AND THE ONE PROPERTY THAT CAN SILENTLY BE WRONG.
///
/// Writing VTK is hard to get catastrophically wrong -- ParaView either opens
/// the file or it does not, and a human notices. What CAN go wrong silently is
/// the resolution: VTK's native cells are linear, so the default path samples a
/// P_k field at the element vertices, and a k = 3 solution is then drawn as if
/// it were k = 1. The picture still looks like a plausible equilibrium. It
/// looks like a coarser mesh, which is exactly what nobody investigates.
///
/// So the assertion is on the POINT COUNT, and it is EXACT rather than an
/// inequality, because the number is derivable and a derived number catches a
/// half-regression that a bound does not. A Lagrange triangle of order d carries
/// ( d + 1 )( d + 2 ) / 2 points, one cell per element, so the piece must hold
/// exactly numElements * ( d + 1 )( d + 2 ) / 2 of them. A regression to linear
/// output gives d = 1 and three points per element; a regression to the OLD
/// subdivision gives d = k, which is a bound away from nothing and is the case
/// worth catching.
///
/// AND THE FIELD IS psi*, WHICH IS WHAT THE DRIVER NOW DRAWS. psi* is P_(k+1),
/// so the subdivision handed to writeVtu() is k+1 and not k -- drawing a
/// degree-(k+1) field at degree-k nodes discards exactly the order the
/// post-processing was done to gain, and it looks like a slightly coarse mesh.
/// This is the same failure SetHighOrderOutput() exists to prevent, one degree
/// further along, and it is why the number below moved from 320 to 480.
BOOST_AUTO_TEST_CASE( theVtkFilesCarryTheHighOrderSolution )
{
	meq::analytic::SolovievEquilibrium const &eq = equilibrium();
	mfem::Mesh mesh = meq::tests::makeMesh( box(), 4 );

	mfem::FunctionCoefficient source( [ &eq ]( mfem::Vector const &x )
		{ return eq.f( x( 0 ), x( 1 ), 0.0 ); } );
	mfem::FunctionCoefficient exact( [ &eq ]( mfem::Vector const &x )
		{ return eq.psi( x( 0 ), x( 1 ) ); } );

	int const degree = 3;
	meq::GradShafranovSolver solver( mesh, degree );
	solver.setSource( source );
	solver.setBoundaryData( exact );
	solver.solve();

	// What the driver writes: psi*, at psi*'s own degree.
	solver.postProcess();
	int const drawn = degree + 1;

	mfem::GridFunction field( solver.flux().FESpace() );
	meq::poloidalField( solver.flux(), field );

	std::string const stem = "meq_test_vtk";
	meq::writeVtu( stem, mesh, solver.postProcessedPotential(), field, drawn );

	// The .pvd is INSIDE the collection directory, which is the thing about
	// this format most likely to be got wrong by a caller. Asserting the
	// documented path is what would catch Output.hpp drifting from MFEM.
	std::string const pvd = stem + "/" + stem + ".pvd";
	std::ifstream index( pvd );
	BOOST_TEST_REQUIRE( index.good(),
	                    "no .pvd at " << pvd << ". ParaViewDataCollection "
	                    "writes the index inside the collection directory, not "
	                    "beside it -- if that changed, Output.hpp and the "
	                    "driver's message are both now wrong" );

	std::string const piece = stem + "/Cycle000000/proc000000.vtu";
	std::ifstream vtu( piece, std::ios::binary );
	BOOST_TEST_REQUIRE( vtu.good(), "no .vtu piece at " << piece );
	std::string const content( ( std::istreambuf_iterator<char>( vtu ) ),
	                            std::istreambuf_iterator<char>() );

	// Base64 payload, plain-XML header: the field names are greppable.
	BOOST_TEST( content.find( "Name=\"psi\"" ) != std::string::npos,
	            "the .vtu declares no psi array" );
	BOOST_TEST( content.find( "Name=\"B\"" ) != std::string::npos,
	            "the .vtu declares no B array -- the poloidal field is what "
	            "this file is written to look at" );

	std::size_t const points = content.find( "NumberOfPoints=\"" );
	BOOST_TEST_REQUIRE( points != std::string::npos );
	long const written = std::strtol(
		content.c_str() + points + std::string( "NumberOfPoints=\"" ).size(),
		nullptr, 10 );

	// The derivation, in code rather than as a literal: one Lagrange cell per
	// element, ( d + 1 )( d + 2 ) / 2 points on a triangle of order d.
	// makeMesh() is Cartesian triangles, so 4 x 4 x 2 = 32 elements and 25
	// vertices; at d = k + 1 = 4 that is 32 x 15 = 480 points, where the old
	// d = k = 3 gave 32 x 10 = 320.
	BOOST_TEST_REQUIRE( mesh.GetElement( 0 )->GetGeometryType()
	                        == mfem::Geometry::TRIANGLE,
	                    "this derivation is the triangular Lagrange cell's; the "
	                    "harness mesh is no longer triangles" );
	long const expected = static_cast<long>( mesh.GetNE() )
	                      *( drawn + 1 )*( drawn + 2 )/2;

	std::printf( "\n  .vtu at k = %d drawing psi* at degree %d: %ld points "
	             "against %d mesh vertices, %d elements (expected %ld)\n",
	             degree, drawn, written, mesh.GetNV(), mesh.GetNE(), expected );
	std::fflush( stdout );

	BOOST_TEST( written == expected,
	            "the .vtu carries " << written << " points where "
	            << mesh.GetNE() << " Lagrange triangles of order " << drawn
	            << " give " << expected << ". At "
	            << static_cast<long>( mesh.GetNE() )*3 << " it is linear output "
	            "-- SetHighOrderOutput() off, or SetLevelsOfDetail() at 1. At "
	            << static_cast<long>( mesh.GetNE() )*( degree + 1 )*( degree + 2 )/2
	            << " it is subdividing at the SOLVE's degree k = " << degree
	            << " while drawing psi*, which is P_" << drawn << ": every "
	            "picture then understates the field by the order the "
	            "post-processing was done to gain" );

	// Written by MFEM, so removed by hand: Scratch takes one file.
	std::remove( piece.c_str() );
	std::remove( ( stem + "/Cycle000000/data.pvtu" ).c_str() );
	std::remove( pvd.c_str() );
	std::remove( ( stem + "/Cycle000000" ).c_str() );
	std::remove( stem.c_str() );
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
	        "double B_Z(Z, R)", "byte inside(Z, R)",
	        "byte extrapolated(Z, R)" } )
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

/// THE ROTATING FIELDS -- n_s AND e phi_0 -- AND THE ONE WAY OF WRITING THEM
/// THAT LOOKS RIGHT AND IS NOT.
///
/// With rotation the density is not a flux function: RoPP (96) makes it
/// n_s( r, psi ), and the whole physics of it is the r^2 in the exponent. So a
/// node's own R is not incidental to the answer, it IS the answer, and that is
/// what makes the obvious route wrong.
///
/// sample(), sampleComponent() and sampleCoefficient() all evaluate at the
/// node's location INSIDE AN ELEMENT, and for a node in the Gamma_h-to-Gamma
/// band that location is its FOOT on Gamma_h rather than the node itself. Only
/// samplePotentialWithFlux() applies the Taylor step outward. So writing n_s
/// through an mfem::Coefficient hands the exponent the foot's radius, and the
/// centrifugal enrichment comes out misplaced by a band width -- silently, in a
/// picture that still looks like a rotating plasma.
///
/// The driver therefore loops the grid itself, pairing the already-sampled psi
/// with sampler.rAt( i ). This measures both routes against the closed form and
/// asserts the gap, so the comment above is a measurement rather than a claim.
BOOST_AUTO_TEST_CASE( theRotatingFieldsUseTheNodesOwnRadius )
{
	// A linear psi, which the flux continuation reproduces exactly in the band --
	// so anything left over is the radius, which is the variable under test.
	double const alpha = 0.7, beta = -1.3;
	int const degree = 2;
	meq::tests::Rectangle const inner{ 0.8, 1.2, -0.2, 0.2 };

	mfem::Mesh mesh = meq::tests::makeMesh( inner, 8 );
	mfem::L2_FECollection collection( degree, mesh.Dimension() );
	mfem::FiniteElementSpace scalar( &mesh, &collection );
	mfem::FiniteElementSpace vector( &mesh, &collection, 2 );

	mfem::FunctionCoefficient exact( [ = ]( mfem::Vector const &x )
		{ return alpha*x( 0 ) + beta*x( 1 ); } );
	mfem::VectorFunctionCoefficient fluxOf( 2,
		[ = ]( mfem::Vector const &x, mfem::Vector &q )
		{
			q( 0 ) = alpha/x( 0 );
			q( 1 ) = beta/x( 0 );
		} );

	mfem::GridFunction potential( &scalar );
	mfem::GridFunction flux( &vector );
	potential.ProjectCoefficient( exact );
	flux.ProjectCoefficient( fluxOf );

	// Two species in normalised units, unequal temperatures so that phi_0 is not
	// the symmetric special case, and a density with a real slope in psi -- so
	// that a node's value depends on BOTH the coordinates the routes disagree
	// about. The psi-dependence of F is measured in the Rotating*Convergence
	// files; what is under test here is where the fields are EVALUATED.
	double const referenceRadius = 1.0;
	double const ionDensity = 0.1;
	std::vector<meq::Species> species( 2 );
	species[ 0 ].mass = 1.0;
	species[ 0 ].charge = 1.0;
	species[ 0 ].temperature = std::make_shared<meq::ConstantProfile const>( 0.6 );
	species[ 0 ].density = std::make_shared<LinearProfile const>( ionDensity,
	                                                              0.5*ionDensity );
	species[ 1 ].mass = 1.0e-4;
	species[ 1 ].charge = -1.0;
	species[ 1 ].temperature = std::make_shared<meq::ConstantProfile const>( 0.4 );
	species[ 1 ].density = meq::neutralisingDensity( species, 1 );

	meq::RotatingSource const rotating( species,
		std::make_shared<meq::ConstantProfile const>( 1.2 ),
		std::make_shared<meq::ConstantProfile const>( 0.05 ),
		referenceRadius, 1.0 );

	// A grid wider than the mesh, and a band around it: this is the curved
	// path's sliver, reproduced without needing a curved boundary.
	meq::GridSampler sampler( mesh, 0.7, 1.3, 41, -0.3, 0.3, 41 );
	int const filled = sampler.extendOutward( 1.0 );
	BOOST_TEST_REQUIRE( filled > 0, "no node was continued, so there is no band "
	                                "and this test measures nothing" );

	double const outside = std::numeric_limits<double>::quiet_NaN();
	std::vector<double> psi;
	sampler.samplePotentialWithFlux( potential, flux, psi, outside );

	// THE ROUTE THE DRIVER TAKES: the node's own R, and the psi that carries the
	// flux continuation.
	std::vector<double> density( psi.size(), outside );
	std::vector<double> ePhi( psi.size(), outside );
	for ( int j = 0; j < sampler.nodesZ(); ++j )
		for ( int i = 0; i < sampler.nodesR(); ++i )
		{
			std::size_t const at = static_cast<std::size_t>( j )*sampler.nodesR() + i;
			if ( !sampler.located( i, j ) || !std::isfinite( psi[ at ] ) )
				continue;

			double const r = sampler.rAt( i );
			density[ at ] = rotating.density( 0, r, psi[ at ] );
			ePhi[ at ] = rotating.potential( r, psi[ at ] );
		}

	// THE ROUTE THAT LOOKS RIGHT: the same closed form as a Coefficient, which
	// the sampler evaluates at the foot.
	mfem::FunctionCoefficient densityOf( [ & ]( mfem::Vector const &x )
		{ return rotating.density( 0, x( 0 ), alpha*x( 0 ) + beta*x( 1 ) ); } );
	std::vector<double> viaCoefficient;
	sampler.sampleCoefficient( densityOf, viaCoefficient, outside );

	double worstInterior = 0.0, worstBand = 0.0, worstBandCoefficient = 0.0;
	int bandNodes = 0;
	for ( int j = 0; j < sampler.nodesZ(); ++j )
		for ( int i = 0; i < sampler.nodesR(); ++i )
		{
			std::size_t const at = static_cast<std::size_t>( j )*sampler.nodesR() + i;
			if ( !sampler.located( i, j ) )
				continue;

			double const r = sampler.rAt( i ), z = sampler.zAt( j );
			// The closed form at the NODE, which is what both routes are trying
			// to compute.
			double const want = rotating.density( 0, r, alpha*r + beta*z );
			double const mine = std::abs( density[ at ] - want )/want;
			double const theirs = std::abs( viaCoefficient[ at ] - want )/want;

			bool const band = r < inner.rMin || r > inner.rMax
			                  || z < inner.zMin || z > inner.zMax;
			if ( band )
			{
				++bandNodes;
				worstBand = std::max( worstBand, mine );
				worstBandCoefficient = std::max( worstBandCoefficient, theirs );
			}
			else
			{
				worstInterior = std::max( worstInterior, mine );
			}
		}

	std::printf( "\n  n_s over %d band nodes: node's own R is %.3e wrong, the "
	             "foot's R is %.3e wrong\n"
	             "                          interior nodes %.3e\n",
	             bandNodes, worstBand, worstBandCoefficient, worstInterior );
	std::fflush( stdout );

	BOOST_TEST( worstInterior < 1.0e-12,
	            "an interior node's density is out by " << worstInterior
	            << ", where no continuation happens at all and psi is linear" );
	BOOST_TEST( worstBand < 1.0e-5,
	            "a band node's density is out by " << worstBand
	            << " using the node's own radius, which is more than the flux "
	            "continuation's own error in psi can account for" );
	BOOST_TEST( worstBandCoefficient > 1.0e-3,
	            "sampling n_s through a Coefficient agrees with the closed form "
	            "at the node to " << worstBandCoefficient << ", so this test is "
	            "no longer measuring the trap it was written for -- either the "
	            "band vanished or the sampler now applies the offset" );
	BOOST_TEST( worstBandCoefficient > 100.0*worstBand,
	            "the two routes differ by only a factor of "
	            << worstBandCoefficient/std::max( 1.0e-16, worstBand )
	            << ", which is not the separation the driver's choice rests on" );

	// AND THEY GO INTO THE FILE, under the names and units a reader is promised.
	if ( !meq::hasNetCDF() )
	{
		BOOST_TEST_MESSAGE( "  built without netcdf-cxx4, skipping the header" );
		return;
	}

	Scratch const file( "meq_test_rotating.nc" );
	{
		meq::NetCDFWriter writer( file.path(), sampler );
		writer.attribute( "source_type", "rotating" );
		writer.attribute( "psi_axis", 0.25 );
		writer.attribute( "normalisation_residual", 0.0 );
		writer.field( "psi", psi, "poloidal flux function", "Wb/rad" );
		writer.field( "n_D", density, "number density of D, RoPP (96)", "m^-3" );
		writer.field( "e_phi_0", ePhi,
		              "elementary charge times phi_0 of RoPP (97), in JOULES",
		              "J" );
		writer.close();
	}

	std::string const command =
		"ncdump -h " + file.path() + " > meq_test_rot_hdr.txt 2>&1";
	BOOST_TEST_REQUIRE( std::system( command.c_str() ) == 0,
	                    "ncdump could not read the rotating file" );

	Scratch const header( "meq_test_rot_hdr.txt" );
	std::ifstream headerIn( header.path() );
	std::string const text( ( std::istreambuf_iterator<char>( headerIn ) ),
	                        std::istreambuf_iterator<char>() );

	std::printf( "  ncdump header:\n" );
	for ( std::string const &needle :
	      { "double n_D(Z, R)", "n_D:units = \"m^-3\"",
	        "double e_phi_0(Z, R)", "e_phi_0:units = \"J\"",
	        ":source_type = \"rotating\"", ":psi_axis = 0.25",
	        ":normalisation_residual = 0." } )
	{
		bool const present = text.find( needle ) != std::string::npos;
		std::printf( "    %-34s %s\n", needle.c_str(), present ? "yes" : "MISSING" );
		BOOST_TEST( present, "the file does not declare '" << needle << "'" );
	}
	std::fflush( stdout );
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
