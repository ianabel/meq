#define BOOST_TEST_MODULE ExtensionConvergence
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <memory>
#include <queue>
#include <vector>

#include "mfem.hpp"

#include "meq/GradShafranov.hpp"
#include "meq/Source.hpp"

#include "analytic/Soloviev.hpp"

/*
 * The stage-5 acceptance test: the NSTX Solov'ev equilibrium on a CURVED
 * boundary that the mesh does not follow, solved on a polygonal subdomain of a
 * uniform background triangulation with the Dirichlet datum carried in from the
 * true boundary along transferring paths.
 *
 * This is the technique of Cockburn and Solano, as
 * refs/HDG-GradShafranov-Adaptive.pdf sections 2.1-2.2 uses it. The claim being
 * tested is a claim about RATES, and a strong one: the design order survives
 * even though dist( Gamma_h, Gamma ) is only O( h ), where every earlier
 * technique needed O( h^(k+1) ). The dist column of the table below is there to
 * show that the gap really is O( h ) -- it sits at about 1.3 h throughout -- so
 * that the rates are not being obtained by accidentally fitting the boundary.
 *
 * The trace column is the size of the trace space, as in the other convergence
 * tables, but here it is not quite the size of the system that is solved: the
 * trace dofs of Gamma_h -- 21 per cent of them on the coarsest mesh, 4 per cent
 * on the finest -- carry no unknown at all and are pinned. See
 * meq::GradShafranovSolver::setExtension.
 *
 * The control is the point of the "without transfer" tables. Solving the same
 * problem on the same D_h with psi = 0 imposed on Gamma_h itself, ignoring the
 * gap, gives 0.84-0.91 for psi and 0.54-0.58 for q at EVERY k. That is what
 * this test would measure if the extension were decorative.
 *
 * WHICH LEVEL SET, AND WHY NOT THE PUBLISHED SEPARATRIX.
 *
 * The plasma boundary is a level set of psi, and psi is defined only up to an
 * additive constant -- psi_1 = 1 is one of the twelve Delta*-harmonic functions
 * of the expansion, so displacing psi by a constant gives the same equilibrium,
 * with the same source F and, to the last bit, the same flux q. Which surface is
 * called psi = 0 is therefore a choice, and it is fixed here by requiring that
 * Gamma be a closed, smooth flux surface.
 *
 * That choice is forced, and the reason is worth recording because it looks at
 * first like a transcription bug. With the twelve coefficients as
 * tests/analytic/Soloviev.hpp carries them, the level set psi = 0 IS NOT CLOSED.
 * Measured: psi has a saddle at ( 0.6960, -1.8070 ) where psi = -8.718e-3, not
 * zero, so { psi < 0 } runs through the saddle into the region below it and on
 * to infinity -- psi ~ -3.5e-4 z^6 for large |z|. The published separatrix
 * cannot be recovered from those coefficients:
 *
 *   - refs/HDG-GradShafranov-Adaptive.pdf prints c_7 and c_10 as the same
 *     number, which Soloviev.hpp already flags as suspicious. Solving Cerfon and
 *     Freidberg's conditions for the NSTX single null with the other eleven
 *     magnitudes fixed and c_10 free gives c_10 = -2.87e-3, sixty-five times the
 *     printed value, and brings the residual of their ten geometric conditions
 *     down from 4.6e-2 to 8.9e-5. c_10 in the paper is very likely a
 *     typesetting duplicate of c_7.
 *   - the paper's minus signs do not survive text extraction from the PDF at
 *     all, so the printed coefficients cannot be checked against it. NOTHING IN
 *     A CONVERGENCE TEST CAN CATCH A WRONG ONE: psi_1 ... psi_12 are
 *     Delta*-harmonic, so every c_i leaves Delta*(psi) = -F untouched and every
 *     rate with it.
 *
 *     That last sentence used to end "the coefficients are the one part of the
 *     benchmark that is asserted nowhere", and that is no longer true:
 *     SolovievGeometryConvergence.cpp evaluates all twelve of Cerfon and
 *     Freidberg's constraints on every set in the fixture, which is the only
 *     kind of check that can see them. It caught a second error in the
 *     replacement coefficients after it caught the first in the published
 *     ones -- see the note at the head of tests/analytic/Soloviev.hpp.
 *
 * And even with the right coefficients the NSTX separatrix would be the wrong
 * geometry for this technique: it passes through an X-point, which is a CORNER
 * of Gamma, and both path families give out there. Measured, with psi as
 * Soloviev.hpp has it: LevelSetPath aborts because the outward normal below the
 * plasma tip never meets the level set -- it runs straight through the X-point
 * into the private-flux region, where psi is negative again -- and VertexConePath
 * aborts for the same reason, its whole admissible fan missing Gamma. A curved
 * reentrant corner on Gamma is a separate problem from the one stage 5 is about.
 *
 * So Gamma here is the closed flux surface psi_nstx = -0.03, written as the zero
 * set of psi := psi_nstx + 0.03. The offset has to clear the saddle at
 * -8.718e-3; 0.03 clears it comfortably and encloses r in [0.467, 1.762],
 * z in [-1.465, 1.272], which is most of the plasma and a proper D shape. The
 * source, the flux, and every convention under test are exactly those of
 * tests/convergence/SolovievConvergence.cpp, which measures the same
 * equilibrium on a fitted rectangle -- so the rates here are comparable with
 * that table even though the errors are not.
 *
 * WHAT THE CONVENTIONS COST WHEN THEY ARE WRONG.
 *
 * Every choice in the extension branch of src/meq/GradShafranov.cpp was settled
 * by running this test with the alternative and reading the table, as stage 2's
 * were. At k = 1 and k = 2 over the first three meshes:
 *
 *   HDGExtensionIntegrator's sign = -1     psi 4.0e-1, rates 1.63 then 0.58
 *   instead of +1                          q   4.1e+0, rates 0.13 then 0.19;
 *                                          at k = 2 the q rate goes negative
 *   its coefficient C = 1/r instead of r   psi flat at 2-5e-2, rates 0.56-0.79
 *                                          q   flat at 0.9-1.8e-1, rates 0.4-0.7
 *   the HDG stabilisation left on          psi 2.5e-2, rates 1.62, 1.67
 *   Gamma_h as well as the interior        q   4.9e-2, rates 1.41, 1.52
 *                                          -- an order lost at k=1, two at k=2
 *   the Gamma_h trace dofs left unpinned   no change in any digit, and the
 *   (SetEssentialBC on the fitted          reduced matrix has no zero row: they
 *   attributes only)                       are coupled to nothing but themselves
 *   B keeps its boundary face integrator   no change in any digit
 *   on Gamma_h (so the flux constraint
 *   is registered there)
 *   both of the last two at once           psi 4.8e+13. With the constraint
 *                                          registered and nothing pinning the
 *                                          dofs it constrains, psihat on Gamma_h
 *                                          becomes a free unknown answering
 *                                          < qhat.n, mu > = 0 -- a natural
 *                                          condition, and the wrong problem.
 *
 * The path family, by contrast, hardly matters here. VertexConePath with 16,
 * 64 or 256 rays and LevelSetPath all give the same errors to three significant
 * figures. LevelSetPath needs a search length of 12 h rather than 6 h to find an
 * endpoint at all, because the staircase boundary of D_h has faces whose normal
 * is nearly tangent to Gamma, and its longest path is then about 3 h rather than
 * 1.3 h. VertexConePath is used because it is the family Cockburn and Solano
 * build for a boundary with no closed-form closest point, and because it is the
 * one whose swept regions tile.
 *
 * WHY THE PER-PAIR SLACK IS LARGER HERE THAN ON THE FITTED RECTANGLE.
 *
 * D_h changes shape with h in a way a fitted mesh does not: which background
 * elements happen to fall entirely inside Omega is not a smooth function of h,
 * and at k = 3 the error wobbles by a factor of two or three between
 * neighbouring meshes. Measured over n = 6 ... 64 at k = 3, the L2 error in psi
 * is not even monotone -- 2.02e-5, 3.03e-5, 3.88e-5 at n = 6, 7, 8. It is not
 * the path family: the four families above give the same numbers. So a rate
 * estimated from one pair of meshes is noisier than on the fitted problem, and
 * the assertions below are made twice: on every pair with a slack of 0.30, and
 * on the rate across the whole sequence with the 0.15 that
 * SolovievConvergence.cpp uses. The tightest pair measured is q at k = 3 on the
 * coarsest step, 3.824.
 */

namespace
{

	meq::analytic::SolovievEquilibrium const &equilibrium()
	{
		static meq::analytic::SolovievEquilibrium const eq
			= meq::analytic::SolovievEquilibrium::nstx();
		return eq;
	}

	/// The additive constant of psi, refixed so that Gamma is closed. See the
	/// header comment for how it was chosen and why it has to be.
	double const psiOffset = 0.03;

	/// The exact solution. Displacing psi by a constant changes neither the
	/// source nor the flux, so this is nstx() throughout.
	double psiExact( double r, double z )
	{
		return equilibrium().psi( r, z ) + psiOffset;
	}

	/// The level set the subdomain and the paths are built from: the exact
	/// solution itself, negative inside Omega. Which sign that is was measured,
	/// not assumed -- see theInteriorOfOmegaIsWherePsiIsNegative below.
	double levelSet( mfem::Vector const &x )
	{
		return psiExact( x( 0 ), x( 1 ) );
	}

	// The background box. It contains Omega with room to spare, keeps r well
	// away from zero -- the operator carries a 1/r and psi carries a log r -- and
	// has sides in the ratio 1:2, so n by 2n cells are square.
	double const rMin = 0.25;
	double const rMax = 1.95;
	double const zMin = -1.75;
	double const zMax = 1.65;

	/// A point inside Omega: the magnetic axis, found by scanning psi.
	double const axisR = 1.322;
	double const axisZ = 0.008;

	struct Measurement
	{
		double h;
		int elements;
		int traceDofs;
		double distance;
		int widened;
		double errorPsi;
		double errorFlux;
	};

	/// Whether every marked element is reachable from every other through shared
	/// faces. { psi < 0 } has components far from the plasma -- psi tends to
	/// minus infinity with |z| -- and the box is placed to exclude them, so this
	/// is a check on the box rather than a filter on the marking.
	bool markedSetIsConnected( mfem::Mesh &mesh, mfem::Array<int> const &marker )
	{
		int seed = -1;
		int marked = 0;
		for ( int e = 0; e < mesh.GetNE(); ++e )
		{
			if ( !marker[ e ] )
				continue;
			marked++;
			if ( seed < 0 )
				seed = e;
		}
		if ( seed < 0 )
			return false;

		mfem::Table const &elementToElement = mesh.ElementToElementTable();
		std::vector<int> seen( mesh.GetNE(), 0 );
		std::queue<int> pending;
		pending.push( seed );
		seen[ seed ] = 1;
		int reached = 1;
		while ( !pending.empty() )
		{
			int const e = pending.front();
			pending.pop();
			int const *row = elementToElement.GetRow( e );
			for ( int i = 0; i < elementToElement.RowSize( e ); ++i )
			{
				int const neighbour = row[ i ];
				if ( neighbour >= 0 && marker[ neighbour ] && !seen[ neighbour ] )
				{
					seen[ neighbour ] = 1;
					reached++;
					pending.push( neighbour );
				}
			}
		}
		return reached == marked;
	}

	/// D_h: the elements of a uniform background triangulation of the box that
	/// lie entirely inside Omega.
	std::unique_ptr<mfem::SubMesh> makeSubdomain( int n, int &gammaH, double &h )
	{
		mfem::Mesh background = mfem::Mesh::MakeCartesian2D(
			n, 2*n, mfem::Element::TRIANGLE, false, rMax - rMin, zMax - zMin );
		background.Transform( []( mfem::Vector const &in, mfem::Vector &out )
		{
			out( 0 ) = in( 0 ) + rMin;
			out( 1 ) = in( 1 ) + zMin;
		} );
		h = ( rMax - rMin )/static_cast<double>( n );

		// extra_refine = 1. The vertex test alone is exact only where Omega is
		// convex, and a flux surface with triangularity is not obviously so.
		// Measured, it makes no difference here: 0, 1, 2 and 3 select exactly the
		// same elements at every mesh in the study, so this surface is convex
		// enough at every resolution used. Kept as the cheap insurance it is.
		mfem::Array<int> marker;
		int const inside = mfem::MarkLevelSetSubdomain( background, levelSet, 0.0,
		                                                marker, 1 );
		BOOST_TEST_REQUIRE( inside > 0, "the subdomain is empty at n = " << n );
		BOOST_TEST_REQUIRE( markedSetIsConnected( background, marker ),
		                    "the elements inside Omega form more than one piece at n = "
		                    << n << ", so the background box has caught a component of "
		                    "{ psi < 0 } that is not the plasma" );

		for ( int e = 0; e < background.GetNE(); ++e )
			background.SetAttribute( e, marker[ e ] ? 1 : 2 );
		background.SetAttributes();

		mfem::Array<int> domainAttr( 1 );
		domainAttr[ 0 ] = 1;
		auto sub = std::make_unique<mfem::SubMesh>(
			mfem::SubMesh::CreateFromDomain( background, domainAttr ) );

		// SubMesh gives the boundary it had to generate one new attribute and
		// leaves any inherited from the parent with the attributes it already had.
		// Here there is nothing inherited: Omega does not reach the box, so the
		// whole of Gamma_h is transferred and none of it is fitted.
		gammaH = sub->bdr_attributes.Max();
		BOOST_TEST_REQUIRE( sub->bdr_attributes.Size() == 1,
		                    "D_h has boundary inherited from the background box at n = "
		                    << n << ", so part of Gamma_h is fitted and the box is too "
		                    "small" );
		BOOST_TEST_REQUIRE( gammaH == background.bdr_attributes.Max() + 1,
		                    "D_h has no boundary of its own at n = " << n );
		return sub;
	}

	/// The largest distance from Gamma_h to Gamma along the paths. This is the
	/// quantity the whole construction is about, and it is O( h ) here -- which
	/// is exactly what earlier techniques could not afford.
	double transferDistance( mfem::SubMesh &sub, int gammaH,
	                         mfem::TransferPath const &path )
	{
		double largest = 0.0;
		mfem::Vector x, xbar;
		for ( int be = 0; be < sub.GetNBE(); ++be )
		{
			if ( sub.GetBdrAttribute( be ) != gammaH )
				continue;
			mfem::FaceElementTransformations *ftr = sub.GetBdrFaceTransformations( be );
			if ( !ftr )
				continue;
			mfem::IntegrationRule const &ir = mfem::IntRules.Get( ftr->GetGeometryType(), 4 );
			for ( int q = 0; q < ir.GetNPoints(); ++q )
			{
				path.Endpoint( *ftr, ir.IntPoint( q ), xbar );
				ftr->Transform( ir.IntPoint( q ), x );
				xbar -= x;
				largest = std::max( largest, xbar.Norml2() );
			}
		}
		return largest;
	}

	/// Solve once and measure. With @a extend the homogeneous datum is carried in
	/// from Gamma along the paths; without it, psi = 0 is imposed on Gamma_h
	/// itself, which ignores the O( h ) gap and is a different problem.
	Measurement measure( int order, int n, bool extend )
	{
		meq::analytic::SolovievEquilibrium const &eq = equilibrium();

		int gammaH = 0;
		double h = 0.0;
		std::unique_ptr<mfem::SubMesh> sub = makeSubdomain( n, gammaH, h );

		mfem::FunctionCoefficient sourceCoeff( [ &eq ]( mfem::Vector const &x )
		{
			return eq.f( x( 0 ), x( 1 ), 0.0 );
		} );
		mfem::FunctionCoefficient psiCoeff( []( mfem::Vector const &x )
		{
			return psiExact( x( 0 ), x( 1 ) );
		} );
		mfem::VectorFunctionCoefficient fluxCoeff( 2, [ &eq ]( mfem::Vector const &x,
		                                                       mfem::Vector &value )
		{
			eq.flux( x( 0 ), x( 1 ), value( 0 ), value( 1 ) );
		} );
		mfem::ConstantCoefficient zero( 0.0 );

		// Six h of search length. The paths are about 1.3 h long, so this is a
		// factor of four of slack; LevelSetPath would need twelve, see the header.
		mfem::VertexConePath path( *sub, gammaH, levelSet, 6.0*h );

		mfem::Array<int> gammaHMarker( gammaH );
		gammaHMarker = 0;
		gammaHMarker[ gammaH - 1 ] = 1;

		meq::GradShafranovSolver solver( *sub, order );
		solver.setSource( sourceCoeff );
		solver.setBoundaryData( zero );
		if ( extend )
			solver.setExtension( path, gammaHMarker );
		solver.solve();

		Measurement point;
		point.h = h;
		point.elements = sub->GetNE();
		point.traceDofs = solver.numTraceDofs();
		point.distance = transferDistance( *sub, gammaH, path );
		point.widened = path.NumWidened();
		point.errorPsi = solver.potentialError( psiCoeff );
		point.errorFlux = solver.fluxError( fluxCoeff );
		return point;
	}

	double rate( double coarseError, double fineError, double refinementRatio )
	{
		return std::log( coarseError/fineError )/std::log( refinementRatio );
	}

	void printTable( char const *what, int order, std::vector<Measurement> const &points )
	{
		std::printf( "\n  Solov'ev NSTX on a curved boundary, %s, k = %d\n", what, order );
		std::printf( "  %8s %7s %8s %9s %5s %14s %6s %14s %6s\n",
		             "h", "elem", "trace", "dist", "wide", "L2(psi)", "rate",
		             "L2(q)", "rate" );
		for ( std::size_t i = 0; i < points.size(); ++i )
		{
			Measurement const &p = points[ i ];
			if ( i == 0 )
			{
				std::printf( "  %8.5f %7d %8d %9.3e %5d %14.6e %6s %14.6e %6s\n",
				             p.h, p.elements, p.traceDofs, p.distance, p.widened,
				             p.errorPsi, "-", p.errorFlux, "-" );
			}
			else
			{
				double const ratio = points[ i - 1 ].h/p.h;
				std::printf( "  %8.5f %7d %8d %9.3e %5d %14.6e %6.3f %14.6e %6.3f\n",
				             p.h, p.elements, p.traceDofs, p.distance, p.widened,
				             p.errorPsi, rate( points[ i - 1 ].errorPsi, p.errorPsi, ratio ),
				             p.errorFlux, rate( points[ i - 1 ].errorFlux, p.errorFlux, ratio ) );
			}
		}
		std::fflush( stdout );
	}

	/// Four dyadic background meshes, so three measured rates per quantity per
	/// order. The finest carries 7334 elements and, at k = 3, 44516 trace dofs;
	/// the whole file runs in about ten seconds.
	std::vector<int> const meshSizes = { 8, 16, 32, 64 };

	/// k+1, less the slack allowed for a rate estimated across the whole
	/// sequence. The same number SolovievConvergence.cpp uses.
	double const rateSlack = 0.15;

	/// And the slack allowed for a rate estimated from a single pair of meshes,
	/// which on an unfitted boundary is a noisier thing than on a fitted one. See
	/// the header comment; the tightest measured is 3.824 against 4.
	double const pairSlack = 0.30;

	std::vector<Measurement> study( int order, bool extend )
	{
		std::vector<Measurement> points;
		points.reserve( meshSizes.size() );
		for ( int n : meshSizes )
			points.push_back( measure( order, n, extend ) );
		return points;
	}

	void checkOrder( int order, double psiCeiling, double fluxCeiling )
	{
		std::vector<Measurement> points = study( order, true );
		printTable( "with transfer", order, points );

		double const expectedPair = order + 1.0 - pairSlack;
		for ( std::size_t i = 1; i < points.size(); ++i )
		{
			double const ratio = points[ i - 1 ].h/points[ i ].h;
			double const ratePsi = rate( points[ i - 1 ].errorPsi, points[ i ].errorPsi, ratio );
			double const rateFlux = rate( points[ i - 1 ].errorFlux, points[ i ].errorFlux, ratio );

			BOOST_TEST( ratePsi >= expectedPair,
			            "k = " << order << ", h = " << points[ i ].h
			            << ": psi converged at " << ratePsi << ", wanted " << expectedPair );
			BOOST_TEST( rateFlux >= expectedPair,
			            "k = " << order << ", h = " << points[ i ].h
			            << ": q converged at " << rateFlux << ", wanted " << expectedPair );
		}

		// The rate across the whole sequence, which is the robust statement: the
		// per-pair rates wobble with the shape of D_h, this one does not.
		double const span = points.front().h/points.back().h;
		double const overallPsi = rate( points.front().errorPsi, points.back().errorPsi, span );
		double const overallFlux = rate( points.front().errorFlux, points.back().errorFlux, span );
		double const expected = order + 1.0 - rateSlack;

		std::printf( "    across the whole sequence: psi %.3f, q %.3f, wanted %.2f\n",
		             overallPsi, overallFlux, expected );
		std::fflush( stdout );

		BOOST_TEST( overallPsi >= expected,
		            "k = " << order << ": psi converged at " << overallPsi
		            << " across the sequence, wanted " << expected );
		BOOST_TEST( overallFlux >= expected,
		            "k = " << order << ": q converged at " << overallFlux
		            << " across the sequence, wanted " << expected );

		// A rate is blind to a solution wrong by a constant factor or a sign, so
		// the absolute error is checked too. The ceilings sit at about three times
		// the measured values.
		BOOST_TEST( points.back().errorPsi < psiCeiling,
		            "k = " << order << ": L2 error in psi is " << points.back().errorPsi
		            << ", above the ceiling " << psiCeiling );
		BOOST_TEST( points.back().errorFlux < fluxCeiling,
		            "k = " << order << ": L2 error in q is " << points.back().errorFlux
		            << ", above the ceiling " << fluxCeiling );

		// Assumption P.1 of the analysis: every vertex of Gamma_h must have a path
		// leaving D_h through both of the faces meeting there. VertexConePath
		// widens the admissible fan when it cannot, and says so; the method may
		// still run, but the estimate no longer covers it.
		for ( Measurement const &p : points )
			BOOST_TEST( p.widened == 0,
			            "k = " << order << ", h = " << p.h << ": " << p.widened
			            << " vertices of Gamma_h needed a widened fan" );
	}

}

/// Which side of the level set is the plasma. Measured rather than assumed: the
/// subdomain marking, the transfer paths and the sign of the datum all read it
/// the same way, so getting it backwards would select the complement of Omega
/// and nothing downstream would say so until the rates did.
BOOST_AUTO_TEST_CASE( theInteriorOfOmegaIsWherePsiIsNegative )
{
	BOOST_TEST( psiExact( axisR, axisZ ) < 0.0,
	            "psi is " << psiExact( axisR, axisZ ) << " at the magnetic axis, so the "
	            "interior is not the negative side of the level set" );

	// And the box's own boundary is entirely outside, which is what makes the
	// subdomain a subdomain.
	double largestOnBox = -1.0e300;
	for ( int i = 0; i <= 200; ++i )
	{
		double const s = static_cast<double>( i )/200.0;
		double const r = rMin + s*( rMax - rMin );
		double const z = zMin + s*( zMax - zMin );
		for ( double value : { psiExact( r, zMin ), psiExact( r, zMax ),
		                       psiExact( rMin, z ), psiExact( rMax, z ) } )
			largestOnBox = std::max( largestOnBox, -value );
	}
	BOOST_TEST( largestOnBox < 0.0,
	            "psi reaches " << -largestOnBox << " on the boundary of the background "
	            "box, so Omega is not strictly inside it" );
}

/// The benchmark before the solver, on the domain that is actually used.
/// SolovievConvergence.cpp makes the same check on its rectangle; this repeats
/// it over Omega, and on the DISPLACED psi, which is what pins the claim that
/// the displacement leaves Delta*(psi) alone.
BOOST_AUTO_TEST_CASE( theDisplacedFluxStillSolvesTheEquation )
{
	meq::analytic::SolovievEquilibrium const &eq = equilibrium();

	int sampled = 0;
	for ( double r = 0.5; r <= 1.75; r += 0.15 )
	{
		for ( double z = -1.4; z <= 1.21; z += 0.2 )
		{
			if ( psiExact( r, z ) > 0.0 )
				continue;
			sampled++;
			double const deltaStar = eq.deltaStarFD( r, z );
			double const minusF = -eq.f( r, z, 0.0 );
			BOOST_TEST( std::abs( deltaStar - minusF ) < 1.0e-5,
			            "at ( " << r << ", " << z << " ): Delta*(psi) = " << deltaStar
			            << " but -F = " << minusF );
		}
	}
	BOOST_TEST( sampled > 20, "only " << sampled << " sample points fell inside Omega" );
}

/// The gap is O( h ), which is the whole claim. If Gamma_h were quietly fitting
/// Gamma the rates below would be measuring the fitted method, so this asserts
/// that it is not: the gap stays a good fraction of h and shrinks no faster.
BOOST_AUTO_TEST_CASE( theGapBetweenGammaHAndGammaIsOrderH )
{
	std::vector<Measurement> points;
	for ( int n : meshSizes )
		points.push_back( measure( 1, n, true ) );

	std::printf( "\n  the transfer distance against h\n" );
	for ( Measurement const &p : points )
		std::printf( "    h = %8.5f   dist = %9.3e   dist/h = %6.3f\n",
		             p.h, p.distance, p.distance/p.h );
	std::fflush( stdout );

	for ( Measurement const &p : points )
	{
		BOOST_TEST( p.distance/p.h > 0.5,
		            "at h = " << p.h << " the gap is only " << p.distance/p.h
		            << " h, so Gamma_h is closer to Gamma than an unfitted mesh should "
		            "put it" );
		BOOST_TEST( p.distance/p.h < 3.0,
		            "at h = " << p.h << " the gap is " << p.distance/p.h
		            << " h, which is further than the paths should have to cross" );
	}
}

BOOST_AUTO_TEST_CASE( orderOneConvergesAtTwo )
{
	// Measured with the coefficient set re-solved with alpha = arcsin(delta):
	// psi 7.403e-05, q 1.149e-04 at h = 0.0266. Kept at ~3x rather than
	// tightened further: which background elements fall inside Omega is not a
	// smooth function of h, so the finest error moves when the domain does.
	checkOrder( 1, 2.3e-4, 3.5e-4 );
}

BOOST_AUTO_TEST_CASE( orderTwoConvergesAtThree )
{
	// Measured: psi 3.589e-07, q 6.219e-07.
	checkOrder( 2, 1.1e-6, 1.9e-6 );
}

BOOST_AUTO_TEST_CASE( orderThreeConvergesAtFour )
{
	// Measured: psi 2.797e-09, q 1.695e-08. Both ceilings were too tight for
	// the corrected coefficients and are raised here, which is the recalibration
	// working as intended rather than a tolerance being relaxed to pass.
	checkOrder( 3, 8.6e-9, 5.2e-8 );
}

/// The cheap, strong check that the transfer is doing something. Solve exactly
/// the same problem on exactly the same D_h with psi = 0 imposed on Gamma_h
/// itself -- which is what a naive treatment of an unfitted boundary does -- and
/// the order collapses to about one for psi and about a half for q, at every k.
/// The absolute error is then four to seven orders of magnitude larger, and does
/// not improve with k at all, because it is the O( h ) gap and not the
/// approximation that limits it.
BOOST_AUTO_TEST_CASE( withoutTheTransferTheRateCollapses )
{
	for ( int order = 1; order <= 3; ++order )
	{
		std::vector<Measurement> const points = study( order, false );
		printTable( "without transfer", order, points );

		std::vector<Measurement> const extended = study( order, true );

		for ( std::size_t i = 1; i < points.size(); ++i )
		{
			double const ratio = points[ i - 1 ].h/points[ i ].h;
			double const ratePsi = rate( points[ i - 1 ].errorPsi, points[ i ].errorPsi, ratio );
			double const rateFlux = rate( points[ i - 1 ].errorFlux, points[ i ].errorFlux, ratio );

			BOOST_TEST( ratePsi < 1.5,
			            "k = " << order << ", h = " << points[ i ].h
			            << ": psi converged at " << ratePsi << " with the datum imposed on "
			            "Gamma_h, which is better than ignoring an O(h) gap should manage" );
			BOOST_TEST( rateFlux < 1.5,
			            "k = " << order << ", h = " << points[ i ].h
			            << ": q converged at " << rateFlux << " with the datum imposed on "
			            "Gamma_h" );
		}

		BOOST_TEST( points.back().errorPsi > 30.0*extended.back().errorPsi,
		            "k = " << order << ": imposing the datum on Gamma_h gives "
		            << points.back().errorPsi << " against " << extended.back().errorPsi
		            << " with the transfer, which is not the difference the technique "
		            "is supposed to make" );
	}
}

/// The extension has to survive the Newton path as well as the linear one, and
/// the Solov'ev source is where that can be checked exactly: F does not depend
/// on psi, so the same problem can be handed over either as a coefficient on the
/// right hand side or as a meq::Source on the non-linear potential mass form,
/// and the two must produce the same numbers on the same D_h.
///
/// What it catches is the one thing the extension branch does differently on the
/// two paths: the HDG stabilisation is registered on the non-linear form there
/// and on the linear one here, and it has to come off Gamma_h in both. Leave it
/// on in one of them and this disagrees -- and the rate tables above, which
/// exercise only the linear path, would not have noticed.
BOOST_AUTO_TEST_CASE( theNewtonPathReproducesTheLinearPathOnTheCurvedBoundary )
{
	meq::analytic::SolovievEquilibrium const &eq = equilibrium();
	meq::SolovievSource const source( eq.getA() );

	int const order = 2;
	int const cells = 16;

	mfem::FunctionCoefficient sourceCoeff( [ &eq ]( mfem::Vector const &x )
	{
		return eq.f( x( 0 ), x( 1 ), 0.0 );
	} );
	mfem::FunctionCoefficient psiCoeff( []( mfem::Vector const &x )
	{
		return psiExact( x( 0 ), x( 1 ) );
	} );
	mfem::ConstantCoefficient zero( 0.0 );

	int gammaH = 0;
	double h = 0.0;
	std::unique_ptr<mfem::SubMesh> linearMesh = makeSubdomain( cells, gammaH, h );
	std::unique_ptr<mfem::SubMesh> newtonMesh = makeSubdomain( cells, gammaH, h );

	mfem::Array<int> gammaHMarker( gammaH );
	gammaHMarker = 0;
	gammaHMarker[ gammaH - 1 ] = 1;

	mfem::VertexConePath linearPath( *linearMesh, gammaH, levelSet, 6.0*h );
	mfem::VertexConePath newtonPath( *newtonMesh, gammaH, levelSet, 6.0*h );

	meq::GradShafranovSolver linearSolver( *linearMesh, order );
	linearSolver.setSource( sourceCoeff );
	linearSolver.setBoundaryData( zero );
	linearSolver.setExtension( linearPath, gammaHMarker );
	linearSolver.solve();

	meq::GradShafranovSolver newtonSolver( *newtonMesh, order );
	newtonSolver.setSource( source );
	newtonSolver.setBoundaryData( zero );
	newtonSolver.setExtension( newtonPath, gammaHMarker );
	newtonSolver.solve();

	BOOST_TEST( ( linearSolver.isExtended() && newtonSolver.isExtended() ),
	            "setExtension() did not take on both paths" );
	BOOST_TEST( newtonSolver.isNonlinear() == true,
	            "the meq::Source overload should have selected the Newton path" );

	mfem::Vector difference( newtonSolver.potential() );
	difference -= linearSolver.potential();
	double const relative = difference.Norml2()/linearSolver.potential().Norml2();

	std::printf( "\n  the curved boundary through both paths, k = %d, n = %d:\n"
	             "    L2(psi) linear %.6e, Newton %.6e, coefficients differ by %.3e "
	             "relative, %d Newton iteration(s)\n",
	             order, cells,
	             linearSolver.potentialError( psiCoeff ),
	             newtonSolver.potentialError( psiCoeff ),
	             relative, newtonSolver.newtonIterations() );
	std::fflush( stdout );

	BOOST_TEST( relative < 1.0e-10,
	            "the linear and Newton paths disagree by " << relative
	            << " relative on a curved boundary, so the extension is not being "
	            "assembled the same way on both" );
	BOOST_TEST( newtonSolver.newtonIterations() <= 2,
	            "Newton took " << newtonSolver.newtonIterations() << " iterations on a "
	            "problem whose dF/dpsi is identically zero" );
}
