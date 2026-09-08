// The ( Psi, theta ) flux-surface output and the per-psi geometry cache:
// INVERSION-PLAN.md stage IN-6.
//
// tests/unit/FluxFamilyTests.cpp is the MFEM-free half -- the cache's
// invalidation contract, the interpolation, and the refusals -- driven by a
// synthetic extractor so that "was this recomputed?" is an exact question. This
// file is the other half: the same machinery over a real solve, plus the
// geometry itself and the file it is written to.
//
// WHAT IS NEW HERE AND WHAT IS INHERITED. The averages are IN-2's and their
// rates are asserted in SurfaceAverageConvergence.cpp; nothing is re-measured.
// What this file measures for the first time is
//
//   * THE ENCLOSED VOLUME AGAINST V', WHICH NEEDS NO REFERENCE VALUE.
//     Green's theorem gives V = closed-integral pi R^2 dz -- a contour integral
//     with no gradient in it at all -- and the coarea formula says
//     V' = | dV/dpsi |. Two completely different integrals over the same nodes,
//     which must agree. SurfaceAverage.hpp records that the third leg of
//     section 3.3's cross-check is missing; this is that leg on one quantity.
//   * THE SEPARATRIX CUT, swept rather than inherited from LIUQE.
//   * THE BAND MASK PER NODE, on the curved path where there really is a band,
//     through the file rather than only in memory.
//   * THE CACHE OVER THE REAL EXTRACTION, bit for bit.
//
// THE d/dpsi IS RICHARDSON-EXTRAPOLATED AND THE PLAIN DIFFERENCE IS KEPT AS THE
// CONTROL. This tree has repeatedly found that a plain central difference floors
// a check at the INSTRUMENT rather than at the quantity -- Zernike's derivative,
// IN-2's averaged Grad-Shafranov identity, IN-3's fit derivative. The enclosed
// volume is a new integrand and is entitled to its own reading.

#define BOOST_TEST_MODULE FluxGridConvergence
#include <boost/test/unit_test.hpp>

#include <cctype>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iterator>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "mfem.hpp"

#include "meq/CriticalPoints.hpp"
#include "meq/FluxExtraction.hpp"
#include "meq/FluxFamily.hpp"
#include "meq/FluxSurfaces.hpp"
#include "meq/GradShafranov.hpp"
#include "meq/Output.hpp"
#include "meq/SurfaceAverage.hpp"

#include "analytic/FluxSurfaceReference.hpp"
#include "analytic/Soloviev.hpp"
#include "convergence/ConvergenceHarness.hpp"

namespace
{
	using meq::tests::Rectangle;
	using Equilibrium = meq::analytic::SolovievEquilibrium;

	double const twoPi = 6.283185307179586476925286766559;

	/// The box. Lifted from SurfaceAverageConvergence.cpp, and for the reason
	/// recorded there: standardBox() cannot hold nstx()'s outer surfaces.
	Rectangle nstxBox()
	{
		return Rectangle{ 0.60, 1.90, -1.10, 1.10 };
	}

	struct ExactAxis
	{
		double r;
		double z;
		double psi;
	};

	/// The magnetic axis of the CLOSED FORM, to round-off. Lifted from
	/// SurfaceAverageConvergence.cpp, which lifted it from
	/// FluxSurfaceConvergence.cpp, where the argument for it is.
	ExactAxis exactAxis( Equilibrium const &eq, double rGuess, double zGuess )
	{
		double r = rGuess;
		double z = zGuess;
		double const step = 1.0e-5;

		for ( int iteration = 0; iteration < 200; ++iteration )
		{
			double gr = 0.0;
			double gz = 0.0;
			eq.gradPsi( r, z, gr, gz );

			double a0 = 0.0;
			double a1 = 0.0;
			double b0 = 0.0;
			double b1 = 0.0;
			double hessian[ 2 ][ 2 ];
			eq.gradPsi( r + step, z, a0, b0 );
			eq.gradPsi( r - step, z, a1, b1 );
			hessian[ 0 ][ 0 ] = ( a0 - a1 )/( 2.0*step );
			hessian[ 1 ][ 0 ] = ( b0 - b1 )/( 2.0*step );
			eq.gradPsi( r, z + step, a0, b0 );
			eq.gradPsi( r, z - step, a1, b1 );
			hessian[ 0 ][ 1 ] = ( a0 - a1 )/( 2.0*step );
			hessian[ 1 ][ 1 ] = ( b0 - b1 )/( 2.0*step );

			double const det = hessian[ 0 ][ 0 ]*hessian[ 1 ][ 1 ]
			                   - hessian[ 0 ][ 1 ]*hessian[ 1 ][ 0 ];
			r += -(  hessian[ 1 ][ 1 ]*gr - hessian[ 0 ][ 1 ]*gz )/det;
			z += -( -hessian[ 1 ][ 0 ]*gr + hessian[ 0 ][ 0 ]*gz )/det;
		}

		ExactAxis axis;
		axis.r = r;
		axis.z = z;
		axis.psi = eq.psi( r, z );
		return axis;
	}

	Equilibrium const &equilibrium()
	{
		static Equilibrium const eq = Equilibrium::nstx();
		return eq;
	}

	ExactAxis const &axis()
	{
		static ExactAxis const found = exactAxis( equilibrium(), 1.05, 0.0 );
		return found;
	}

	/// One solve, kept alive so that a tracer can borrow the fields. Lifted from
	/// SurfaceAverageConvergence.cpp.
	class SolvedEquilibrium
	{
		public:
			SolvedEquilibrium( int orderIn, int n )
				: mesh( meq::tests::makeMesh( nstxBox(), n ) ),
				  sourceCoeff( []( mfem::Vector const &x )
				  {
					  return equilibrium().f( x( 0 ), x( 1 ), 0.0 );
				  } ),
				  psiCoeff( []( mfem::Vector const &x )
				  {
					  return equilibrium().psi( x( 0 ), x( 1 ) );
				  } ),
				  solver( mesh, orderIn )
			{
				solver.setSource( sourceCoeff );
				solver.setBoundaryData( psiCoeff );
				solver.solve();
				solver.postProcess();
			}

			SolvedEquilibrium( SolvedEquilibrium const & ) = delete;
			SolvedEquilibrium &operator=( SolvedEquilibrium const & ) = delete;

			meq::GradShafranovSolver &theSolver()
			{
				return solver;
			}

			meq::CriticalPoint magneticAxis() const
			{
				meq::CriticalPointFinder const finder( solver );
				return finder.findAxis();
			}

		private:
			mfem::Mesh mesh;
			mfem::FunctionCoefficient sourceCoeff;
			mfem::FunctionCoefficient psiCoeff;
			meq::GradShafranovSolver solver;
	};

	/// THE OUTER CUT ON THE FITTED FIXTURE IS SET BY THE BOX AND NOT BY THE
	/// PHYSICS, and it is worth saying so where the number appears.
	/// SurfaceAverageConvergence.cpp records that Psi_N = 0.75 on nstx() is not
	/// measurable on a fitted rectangle at all -- enclosing it leaves under one
	/// cell of margin -- and Psi_N = 1 is the separatrix, which passes through
	/// an X-point no rectangle in this tree contains. Measured here: the family
	/// stops closing at about Psi_N = 0.70 on nstxBox().
	///
	/// So the fitted cases below run to 0.50 and the OUTER end of the cut is
	/// asked on the curved fixture, where Gamma is a closed flux surface and
	/// Psi_N = 1 is a real place on the mesh.
	meq::FluxFamilyOptions defaultOptions()
	{
		meq::FluxFamilyOptions options;
		options.surfaces = 12;
		options.angles = 256;
		options.innerCut = 0.05;
		options.outerCut = 0.50;
		return options;
	}

	/// A temporary file that removes itself. Lifted from OutputConvergence.cpp.
	class Scratch
	{
		public:
			explicit Scratch( std::string nameIn ) : name( std::move( nameIn ) ) { }
			~Scratch() { std::remove( name.c_str() ); }
			Scratch( Scratch const & ) = delete;
			Scratch &operator=( Scratch const & ) = delete;
			std::string const &path() const { return name; }

		private:
			std::string name;
	};

	double relative( double value, double target )
	{
		return std::abs( value - target )/std::abs( target );
	}

	// -----------------------------------------------------------------------
	// The curved fixture, lifted from FluxSurfaceConvergence.cpp.
	//
	// WHY THE OUTER END OF THE CUT IS ASKED HERE AND NOT ON THE FITTED BOX.
	// Gamma is the closed flux surface psi_nstx = -0.03, so Psi_N = 1 is a real
	// place on the mesh rather than a separatrix through an X-point no
	// rectangle contains -- and Omega_h is INSCRIBED in it, so the outermost
	// surfaces cross the band between Gamma_h and Gamma. That is the
	// configuration a production run is in, and it is the only one in which the
	// per-node band mask is anything but zero.
	// -----------------------------------------------------------------------

	double const curvedOffset = 0.03;

	Equilibrium const &curvedEquilibrium()
	{
		static Equilibrium const eq = Equilibrium::nstx();
		return eq;
	}

	double curvedLevelSet( mfem::Vector const &x )
	{
		return curvedEquilibrium().psi( x( 0 ), x( 1 ) ) + curvedOffset;
	}

	Rectangle curvedBox()
	{
		return Rectangle{ 0.25, 1.95, -1.75, 1.65 };
	}

	std::unique_ptr<mfem::SubMesh> makeCurvedSubdomain( int n, int &gammaH,
	                                                    double &h )
	{
		Rectangle const box = curvedBox();

		// THE BACKGROUND MUST OUTLIVE THE SubMesh CUT FROM IT. mfem::SubMesh
		// keeps a POINTER to its parent, mfem::VertexConePath's cone reads the
		// parent's edges, and this function returns the SubMesh -- so a local
		// background is walked as freed memory. CLAUDE.md records this costing a
		// segfault three MFEM frames deep with no meq frame in the trace.
		static std::vector<std::unique_ptr<mfem::Mesh>> backgrounds;
		backgrounds.push_back( std::make_unique<mfem::Mesh>(
			mfem::Mesh::MakeCartesian2D( n, 2*n, mfem::Element::TRIANGLE, false,
			                             box.width(), box.height() ) ) );
		mfem::Mesh &background = *backgrounds.back();
		background.Transform( [ box ]( mfem::Vector const &in, mfem::Vector &out )
		{
			out( 0 ) = in( 0 ) + box.rMin;
			out( 1 ) = in( 1 ) + box.zMin;
		} );
		h = box.width()/static_cast<double>( n );

		mfem::Array<int> marker;
		int const inside = mfem::MarkLevelSetSubdomain( background,
		                                                curvedLevelSet, 0.0,
		                                                marker, 1 );
		BOOST_TEST_REQUIRE( inside > 0, "the subdomain is empty at n = " << n );

		for ( int e = 0; e < background.GetNE(); ++e )
			background.SetAttribute( e, marker[ e ] ? 1 : 2 );
		background.SetAttributes();

		mfem::Array<int> domainAttr( 1 );
		domainAttr[ 0 ] = 1;
		auto sub = std::make_unique<mfem::SubMesh>(
			mfem::SubMesh::CreateFromDomain( background, domainAttr ) );

		gammaH = sub->bdr_attributes.Max();
		BOOST_TEST_REQUIRE( sub->bdr_attributes.Size() == 1,
		                    "D_h has boundary inherited from the background box "
		                    "at n = " << n );
		return sub;
	}

	/// One curved solve, kept alive with its subdomain and its path family.
	/// Member order is load bearing: the solver borrows the path, the path
	/// borrows the mesh, and destruction runs the other way.
	class CurvedSolve
	{
		public:
			CurvedSolve( int orderIn, int n )
				: sub( makeCurvedSubdomain( n, gammaH, meshSize ) ),
				  sourceCoeff( []( mfem::Vector const &x )
				  {
					  return curvedEquilibrium().f( x( 0 ), x( 1 ), 0.0 );
				  } ),
				  zero( 0.0 ),
				  path( *sub, gammaH, curvedLevelSet, 6.0*meshSize ),
				  solver( *sub, orderIn )
			{
				marker.SetSize( gammaH );
				marker = 0;
				marker[ gammaH - 1 ] = 1;

				solver.setSource( sourceCoeff );
				solver.setBoundaryData( zero );
				solver.setExtension( path, marker );
				solver.solve();
				solver.postProcess();
			}

			CurvedSolve( CurvedSolve const & ) = delete;
			CurvedSolve &operator=( CurvedSolve const & ) = delete;

			meq::GradShafranovSolver &theSolver() { return solver; }
			mfem::VertexConePath const &paths() const { return path; }
			mfem::Array<int> const &gammaHMarker() const { return marker; }
			double h() const { return meshSize; }

			meq::CriticalPoint magneticAxis() const
			{
				meq::CriticalPointFinder const finder( solver );
				return finder.findAxis();
			}

		private:
			int gammaH = 0;
			double meshSize = 0.0;
			std::unique_ptr<mfem::SubMesh> sub;
			mfem::Array<int> marker;
			mfem::FunctionCoefficient sourceCoeff;
			mfem::ConstantCoefficient zero;
			mfem::VertexConePath path;
			meq::GradShafranovSolver solver;
	};
}

// ---------------------------------------------------------------------------
// The geometry
// ---------------------------------------------------------------------------

/// The family reproduces the exact surfaces and the exact averages.
///
/// NOT A RATE STUDY: IN-2 measured the averages at k+1 in
/// SurfaceAverageConvergence.cpp and nothing here would add to it. What this
/// asserts is that assembling a whole family through meq::extractFluxSurfaces()
/// -- levels laid out in rho, one trace and one fit each -- delivers what the
/// one-surface route already delivers, against the converged reference on the
/// EXACT field.
BOOST_AUTO_TEST_CASE( theFamilyReproducesTheExactSurfaces )
{
	SolvedEquilibrium solved( 2, 48 );
	meq::FluxSurfaceFamily const family =
		meq::extractFluxSurfaces( solved.theSolver(), defaultOptions() );

	BOOST_TEST_REQUIRE( family.size() == 12u );
	BOOST_TEST( family.angles == 256u );
	BOOST_TEST( !family.safetyFactorAvailable );

	std::printf( "  %8s %10s %14s %14s %10s %14s %14s %10s\n",
	             "Psi_N", "rho", "V'", "V' exact", "rel", "<R^-2>",
	             "<R^-2> exact", "rel" );

	double worstVPrime = 0.0;
	double worstInverse = 0.0;
	double worstPosition = 0.0;

	for ( meq::FluxSurface const &surface : family.surfaces )
	{
		meq::analytic::SurfaceQuadrature const exact =
			meq::analytic::surfaceQuadrature( equilibrium(), axis().r, axis().z,
			                                  surface.level, 2048, 2.0 );

		double const exactInverse = exact.average(
			[]( meq::analytic::SurfacePoint const &p )
			{ return 1.0/( p.r*p.r ); } );

		double const vRel = relative( surface.vPrime, exact.vPrime );
		double const iRel = relative( surface.inverseRSquared, exactInverse );
		worstVPrime = std::max( worstVPrime, vRel );
		worstInverse = std::max( worstInverse, iRel );

		// The position error: the first-order distance from each node to the
		// exact contour, | psi( x ) - c | / | grad psi |.
		for ( std::size_t j = 0; j < surface.count(); ++j )
		{
			double gr = 0.0;
			double gz = 0.0;
			equilibrium().gradPsi( surface.r[ j ], surface.z[ j ], gr, gz );
			double const magnitude = std::hypot( gr, gz );
			double const distance =
				std::abs( equilibrium().psi( surface.r[ j ], surface.z[ j ] )
				          - surface.level )/magnitude;
			worstPosition = std::max( worstPosition, distance );
		}

		std::printf( "  %8.4f %10.4f %14.6e %14.6e %10.2e %14.6e %14.6e %10.2e\n",
		             surface.normalisedFlux, surface.radial, surface.vPrime,
		             exact.vPrime, vRel, surface.inverseRSquared, exactInverse,
		             iRel );
	}

	std::printf( "  worst: V' %.3e, <R^-2> %.3e, position %.3e m\n",
	             worstVPrime, worstInverse, worstPosition );
	std::fflush( stdout );

	BOOST_TEST( worstVPrime < 1.0e-6 );
	BOOST_TEST( worstInverse < 1.0e-6 );
	BOOST_TEST( worstPosition < 1.0e-6 );

	// The fitted path has no band, so every mask is zero. Asserted rather than
	// assumed: a zero is the correct answer here and not an absent one.
	BOOST_TEST( family.extendedNodes() == 0 );
	BOOST_TEST( !family.crossesBand() );
	for ( meq::FluxSurface const &surface : family.surfaces )
	{
		BOOST_TEST_REQUIRE( surface.transverse );
		BOOST_TEST( surface.stalledRays == 0 );
	}
}

/// V' against the enclosed volume, through the coarea formula, with nothing
/// else in it.
///
/// V' = closed-integral 2 pi R dl / | grad psi | divides by the flux at every
/// node. V = closed-integral pi R^2 dz has no gradient in it at all. The coarea
/// formula says | dV/dpsi | is the first, so the two are a check on each other
/// with no reference value anywhere -- which is the property this tree keeps
/// looking for and SurfaceAverage.hpp records as missing a third leg.
///
/// MEQ'S V' IS THE MAGNITUDE, AND THAT IS WORTH SAYING ONCE. It is a sum of
/// positive weights, so it is positive by construction, while dV/dpsi carries
/// the sign of the direction the enclosed volume grows in -- which on nstx(),
/// where the axis is an interior MINIMUM of psi, is positive, and on a fixture
/// whose axis is a maximum is negative. The identity is | dV/dpsi | = V'.
BOOST_AUTO_TEST_CASE( theEnclosedVolumeAndVPrimeAgreeThroughTheCoareaFormula )
{
	SolvedEquilibrium solved( 3, 48 );
	meq::ContourTracer const tracer( solved.theSolver() );
	meq::CriticalPoint const found = solved.magneticAxis();

	double const span = std::abs( found.psi );
	double const level = meq::fluxAtNormalised( 0.5, found.psi, 0.0 );

	auto fitAt = [ & ]( double at )
	{
		return tracer.fitByAngle( tracer.traceFromAxis( at, found ), found, 512 );
	};

	meq::AngleParametrisation const here = fitAt( level );
	meq::SurfaceAverages const averages = meq::surfaceAverages( tracer, here );

	std::printf( "  V' from the average: %.10e\n", averages.vPrime );
	std::printf( "  %10s %18s %10s %18s %10s %10s\n",
	             "step/|psi|", "|dV/dpsi| plain", "rel", "Richardson", "rel",
	             "rho' diffd" );

	double bestPlain = 1.0;
	double bestRichardson = 1.0;
	double bestDifferenced = 1.0;

	for ( double fraction : { 0.04, 0.02, 0.01, 0.005 } )
	{
		double const step = fraction*span;

		// meq::enclosedVolume() is the SHIPPING routine, not a copy of it here:
		// meq::FluxSurface::volume is exactly this, so what the coarea formula
		// checks is what the family carries.
		double const dPlain =
			( meq::enclosedVolume( fitAt( level + step ) )
			  - meq::enclosedVolume( fitAt( level - step ) ) )/( 2.0*step );
		double const dHalf =
			( meq::enclosedVolume( fitAt( level + 0.5*step ) )
			  - meq::enclosedVolume( fitAt( level - 0.5*step ) ) )/step;
		double const dRichardson = ( 4.0*dHalf - dPlain )/3.0;

		// THE METRIC TRAP, ON A THIRD INTEGRAND. rho' by central difference of
		// the neighbouring rho_j instead of pointwise from the solved flux. IN-1
		// measured this at 7.03 against 1.97 on an arc length and IN-2 measured
		// it again on an average; the Green's-theorem volume is a new integrand
		// and is entitled to its own reading.
		double const dDifferenced =
			( meq::enclosedVolumeByDifferencedSlope( fitAt( level + 0.5*step ) )
			  - meq::enclosedVolumeByDifferencedSlope(
			      fitAt( level - 0.5*step ) ) )/step;
		double const dDifferencedPlain =
			( meq::enclosedVolumeByDifferencedSlope( fitAt( level + step ) )
			  - meq::enclosedVolumeByDifferencedSlope(
			      fitAt( level - step ) ) )/( 2.0*step );
		double const dDifferencedRich =
			( 4.0*dDifferenced - dDifferencedPlain )/3.0;

		double const relPlain = relative( std::abs( dPlain ), averages.vPrime );
		double const relRich = relative( std::abs( dRichardson ),
		                                 averages.vPrime );
		double const relDiff = relative( std::abs( dDifferencedRich ),
		                                 averages.vPrime );
		bestPlain = std::min( bestPlain, relPlain );
		bestRichardson = std::min( bestRichardson, relRich );
		bestDifferenced = std::min( bestDifferenced, relDiff );

		std::printf( "  %10.4f %18.10e %10.2e %18.10e %10.2e %10.2e\n",
		             fraction, std::abs( dPlain ), relPlain,
		             std::abs( dRichardson ), relRich, relDiff );
	}

	std::printf( "  best plain %.3e, best Richardson %.3e, ratio %.1f\n",
	             bestPlain, bestRichardson, bestPlain/bestRichardson );
	std::printf( "  best with rho' DIFFERENCED %.3e, which is %.1fx worse\n",
	             bestDifferenced, bestDifferenced/bestRichardson );
	std::fflush( stdout );

	BOOST_TEST( bestRichardson < bestPlain,
	            "the Richardson extrapolation did not beat the plain central "
	            "difference, so the check is floored by the instrument" );
	BOOST_TEST( bestRichardson < 1.0e-5,
	            "the enclosed volume and V' disagree, and they are two "
	            "completely different integrals over the same nodes" );
	BOOST_TEST( bestDifferenced > 10.0*bestRichardson,
	            "the differenced metric did as well as the pointwise one, so "
	            "either this fixture cannot see the trap or the pointwise route "
	            "is not being used" );
}

// ---------------------------------------------------------------------------
// The cut
// ---------------------------------------------------------------------------

/// Where the outer end of the cut has to go, measured rather than inherited.
///
/// INVERSION-PLAN.md section 8's second risk says to decide meq's cut
/// deliberately and record it, and names what production codes do: LIUQE cuts at
/// Psi_N = 0.95, FreeGS extrapolates outside [ 0.01, 0.99 ]. Repeating LIUQE's
/// number with LIUQE's reason would be inheriting a decision rather than making
/// one, so this is the sweep the default is set from.
///
/// AND THE ANSWER IS NOT "IT FAILS AT X", WHICH IS THE FINDING. Nothing fails.
/// The trace closes, the fit converges, no ray stalls, transversality barely
/// moves, and | psi_h - c | sits at 1e-13 at EVERY level right up to
/// Psi_N = 0.995. What changes is what the surface is MADE OF: on the curved
/// path Omega_h is inscribed in Gamma, so the outer surfaces cross the band and
/// the extension answers for them -- as confidently as an element does, and with
/// nothing in the residual to say which answered. The cut is therefore a
/// decision about data provenance and not about convergence, and the per-node
/// mask is the only signal there is.
BOOST_AUTO_TEST_CASE( theOuterCutIsWhereTheSurfaceStopsBeingSolvedData )
{
	std::printf( "  %4s %8s %8s %8s %10s %12s %10s %10s %14s\n",
	             "n", "h", "Psi_N", "band", "band/N", "deepest", "transv",
	             "stalled", "worst resid" );

	double fractionCoarse = -1.0;
	double fractionFine = -1.0;
	double residualInside = 0.0;
	double residualOutside = 0.0;

	for ( int n : { 40, 80 } )
	{
		CurvedSolve solved( 2, n );
		meq::ContourTracer tracer( solved.theSolver() );
		tracer.setBandExtension( meq::BandExtension::TransferLift,
		                         solved.gammaHMarker(), &solved.paths() );

		meq::CriticalPoint const found = solved.magneticAxis();

		for ( double target : { 0.50, 0.80, 0.90, 0.95, 0.98, 0.99, 0.995 } )
		{
			double const level = meq::fluxAtNormalised( target, found.psi, 0.0 );
			meq::Contour const contour = tracer.traceFromAxis( level, found );
			BOOST_TEST_REQUIRE( contour.closed(),
			                    "the contour at Psi_N = " << target
			                    << " did not close at n = " << n );
			meq::AngleParametrisation const fit =
				tracer.fitByAngle( contour, found, 256 );

			double const fraction = static_cast<double>( fit.extendedNodes )
				/static_cast<double>( fit.count() );

			std::printf( "  %4d %8.4f %8.4f %8d %10.4f %12.4e %10.4f %10d "
			             "%14.4e\n",
			             n, solved.h(), target, fit.extendedNodes, fraction,
			             fit.deepestBandNode, fit.transversality,
			             fit.stalledRays, fit.worstResidual );

			if ( target == 0.95 )
			{
				if ( n == 40 )
					fractionCoarse = fraction;
				else
					fractionFine = fraction;
			}

			if ( n == 80 && target == 0.50 )
				residualInside = fit.worstResidual;
			if ( n == 80 && target == 0.995 )
				residualOutside = fit.worstResidual;

			// Nothing gives out. Asserted rather than printed, because the whole
			// point of this case is that the extraction stays healthy while the
			// data behind it changes.
			BOOST_TEST( fit.stalledRays == 0 );
			BOOST_TEST( fit.transverse );
		}
	}

	std::printf( "  band fraction at Psi_N = 0.95: %.4f at n = 40, %.4f at "
	             "n = 80, ratio %.2f\n",
	             fractionCoarse, fractionFine,
	             fractionCoarse/std::max( fractionFine, 1.0e-12 ) );
	std::printf( "  worst residual at n = 80: %.4e at Psi_N = 0.50, %.4e at "
	             "0.995\n", residualInside, residualOutside );
	std::fflush( stdout );

	// THE BAND FRACTION FALLS WITH THE MESH, so the outer cut is a resolution
	// statement rather than a geometric one: the band is O( h ) wide and the
	// surface at a fixed Psi_N is a fixed distance inside Gamma.
	BOOST_TEST( fractionCoarse > fractionFine,
	            "refining the mesh did not reduce the share of the surface that "
	            "is band data, so the band is not O( h ) here and the reasoning "
	            "behind the default cut does not apply" );

	// AND THE RESIDUAL SAYS NOTHING, which is why the mask is load bearing. If
	// this ever fails there IS another signal and the cut could be discovered
	// rather than decided -- which would be worth knowing.
	BOOST_TEST( residualOutside < 100.0*residualInside,
	            "the fit residual grew as the surface entered the band, so it is "
	            "a signal after all" );

	// The cut is a decision, and it is refused when it is not one.
	CurvedSolve small( 2, 24 );
	meq::FluxFamilyOptions bad;
	bad.outerCut = 1.0;
	BOOST_CHECK_THROW( meq::extractFluxSurfaces( small.theSolver(), bad ),
	                   std::invalid_argument );
	bad.outerCut = 0.95;
	bad.innerCut = 0.0;
	BOOST_CHECK_THROW( meq::extractFluxSurfaces( small.theSolver(), bad ),
	                   std::invalid_argument );
	bad.innerCut = 0.96;
	BOOST_CHECK_THROW( meq::extractFluxSurfaces( small.theSolver(), bad ),
	                   std::invalid_argument );
	bad.innerCut = 0.05;
	bad.surfaces = 1;
	BOOST_CHECK_THROW( meq::extractFluxSurfaces( small.theSolver(), bad ),
	                   std::invalid_argument );
}

// ---------------------------------------------------------------------------
// The file
// ---------------------------------------------------------------------------

namespace
{
	/// Read a whole file into a string.
	std::string slurp( std::string const &path )
	{
		std::ifstream in( path );
		return std::string( ( std::istreambuf_iterator<char>( in ) ),
		                    std::istreambuf_iterator<char>() );
	}

	/// The values of one variable out of an `ncdump -v <name>` dump, as
	/// doubles.
	///
	/// PARSED FROM ncdump RATHER THAN READ BACK THROUGH THE WRITER'S OWN
	/// LIBRARY CALLS, so that the check does not share code with the thing it
	/// checks -- which is what OutputConvergence.cpp does for the ( R, Z ) grid
	/// and for the same reason.
	std::vector<double> dumpedValues( std::string const &text,
	                                  std::string const &name )
	{
		std::size_t const data = text.find( "\ndata:" );
		BOOST_TEST_REQUIRE( data != std::string::npos,
		                    "ncdump printed no data section" );

		std::size_t const at = text.find( " " + name + " =", data );
		BOOST_TEST_REQUIRE( at != std::string::npos,
		                    "ncdump printed no values for " << name );

		std::size_t const begin = text.find( '=', at ) + 1;
		std::size_t const end = text.find( ';', begin );
		BOOST_TEST_REQUIRE( end != std::string::npos );

		std::vector<double> values;
		std::string const body = text.substr( begin, end - begin );
		std::size_t cursor = 0;
		while ( cursor < body.size() )
		{
			while ( cursor < body.size()
			        && ( body[ cursor ] == ',' || std::isspace(
			                 static_cast<unsigned char>( body[ cursor ] ) ) ) )
				++cursor;
			if ( cursor >= body.size() )
				break;

			std::size_t used = 0;
			values.push_back( std::stod( body.substr( cursor ), &used ) );
			cursor += used;
		}
		return values;
	}
}

/// The ( Psi, theta ) file, with the band mask read back out of it.
///
/// A COUNT IS NOT A MASK, and this is that obligation one dimension up.
/// CLAUDE.md records that the ( R, Z ) grid file carried `extrapolated_nodes`
/// as an attribute and nothing else, so nothing downstream could tell WHICH of
/// 1667 nodes had been continued outward from Gamma_h. A flux surface is worse
/// placed still: it can be inside Omega_h at one theta and outside it at the
/// next, so a per-SURFACE flag under-reports in the middle of a band excursion,
/// which is exactly where a q( psi ) profile is being read.
///
/// So the file carries `extrapolated( flux, theta )` per node, `band( flux )`
/// as the per-surface summary, and `extrapolated_nodes` as the total -- and the
/// three are checked against each other here, through ncdump rather than
/// through the writer's own library calls.
BOOST_AUTO_TEST_CASE( theFluxGridFileCarriesTheBandMaskPerNode )
{
	if ( !meq::hasNetCDF() )
	{
		BOOST_TEST_MESSAGE( "  built without netcdf-cxx4, skipping" );
		return;
	}

	CurvedSolve solved( 2, 40 );
	meq::ContourTracer tracer( solved.theSolver() );
	tracer.setBandExtension( meq::BandExtension::TransferLift,
	                         solved.gammaHMarker(), &solved.paths() );

	meq::FluxFamilyOptions options;
	options.surfaces = 10;
	options.angles = 64;
	options.innerCut = 0.05;
	options.outerCut = 0.95;
	// A g( psi ) so that the safety factor column is exercised. Constant, which
	// is what a vacuum toroidal field is, and enough to pin the algebra.
	options.toroidalField = []( double ) { return 2.5; };

	meq::FluxSurfaceFamily const family =
		meq::extractFluxSurfaces( tracer, solved.magneticAxis(), 0.0, options );

	BOOST_TEST_REQUIRE( family.size() == 10u );
	BOOST_TEST( family.safetyFactorAvailable );

	// The mask has to discriminate, or the check below is vacuous.
	BOOST_TEST_REQUIRE( family.extendedNodes() > 0,
	                    "no node of this family is band data, so the mask "
	                    "cannot be tested here" );
	BOOST_TEST_REQUIRE( family.extendedNodes()
	                    < static_cast<int>( family.size()*family.angles ),
	                    "every node is band data, so a mask of all ones would "
	                    "pass this case" );

	Scratch const file( "meq_test_fluxgrid.nc" );
	{
		meq::FluxGridWriter writer( file.path(), family );
		writer.attribute( "title", "meq flux-surface grid test" );
		writer.attribute( "polynomial_degree", 2 );
		writer.close();
	}

	Scratch const header( "meq_test_fluxgrid_hdr.txt" );
	BOOST_TEST_REQUIRE( std::system( ( "ncdump -h " + file.path() + " > "
	                                   + header.path() + " 2>&1" ).c_str() )
	                    == 0,
	                    "ncdump could not read the file meq just wrote" );

	std::string const text = slurp( header.path() );
	std::printf( "  ncdump header:\n" );
	for ( std::string const &needle :
	      { "flux = 10", "theta = 64", "double R(flux, theta)",
	        "double Z(flux, theta)", "byte extrapolated(flux, theta)",
	        "double V_prime(flux)", "double volume(flux)",
	        "double inverse_R_squared(flux)", "double safety_factor(flux)",
	        "byte band(flux)", "double rho(flux)",
	        "double normalised_flux(flux)", ":flux_label" } )
	{
		bool const present = text.find( needle ) != std::string::npos;
		std::printf( "    %-38s %s\n", needle.c_str(),
		             present ? "yes" : "MISSING" );
		BOOST_TEST( present, "the file does not declare '" << needle
		            << "', so its shape is not what the header promises" );
	}
	std::fflush( stdout );

	// And the mask itself, value by value.
	Scratch const dump( "meq_test_fluxgrid_mask.txt" );
	BOOST_TEST_REQUIRE( std::system(
		( "ncdump -v extrapolated,band " + file.path() + " > " + dump.path()
		  + " 2>&1" ).c_str() ) == 0 );

	std::string const dumped = slurp( dump.path() );
	std::vector<double> const mask = dumpedValues( dumped, "extrapolated" );
	std::vector<double> const perSurface = dumpedValues( dumped, "band" );

	BOOST_TEST_REQUIRE( mask.size() == family.size()*family.angles );
	BOOST_TEST_REQUIRE( perSurface.size() == family.size() );

	// THE ROW SUMS ARE THE INDEPENDENT CHECK. extendedNodes comes from the fit
	// and never passes through the writer, so a transposed write -- the failure
	// a rectangular array is most likely to have and least likely to announce --
	// would scramble them.
	int total = 0;
	std::printf( "  %8s %8s %10s %8s %8s\n",
	             "Psi_N", "band", "from fit", "band()", "expect" );
	for ( std::size_t i = 0; i < family.size(); ++i )
	{
		int row = 0;
		for ( std::size_t j = 0; j < family.angles; ++j )
			row += mask[ i*family.angles + j ] != 0.0 ? 1 : 0;
		total += row;

		std::printf( "  %8.4f %8d %10d %8.0f %8d\n",
		             family.surfaces[ i ].normalisedFlux, row,
		             family.surfaces[ i ].extendedNodes, perSurface[ i ],
		             family.surfaces[ i ].crossesBand ? 1 : 0 );

		BOOST_TEST( row == family.surfaces[ i ].extendedNodes,
		            "the mask written for surface " << i << " does not agree "
		            "with the count the fit reported" );
		BOOST_TEST( ( perSurface[ i ] != 0.0 )
		            == family.surfaces[ i ].crossesBand );
		BOOST_TEST( ( perSurface[ i ] != 0.0 ) == ( row > 0 ),
		            "band( flux ) disagrees with its own per-node mask" );
	}
	std::fflush( stdout );

	BOOST_TEST( total == family.extendedNodes() );

	// The attribute is the total and never a substitute for the mask.
	std::string const needle = ":extrapolated_nodes = " + std::to_string( total );
	BOOST_TEST( text.find( needle ) != std::string::npos,
	            "the file's extrapolated_nodes attribute is not the sum of its "
	            "own mask" );

	// The safety factor is algebra over two measured quantities, and a constant
	// is exactly as easy to get wrong here as anywhere.
	for ( meq::FluxSurface const &surface : family.surfaces )
	{
		double const expected = surface.vPrime*2.5*surface.inverseRSquared
			/( 4.0*3.141592653589793*3.141592653589793 );
		BOOST_TEST( surface.safetyFactor == expected,
		            boost::test_tools::tolerance( 1.0e-14 ) );
	}

	// An empty family has no grid, and a ragged one cannot be a rectangle.
	meq::FluxSurfaceFamily const nothing;
	BOOST_CHECK_THROW( meq::FluxGridWriter( file.path(), nothing ),
	                   std::invalid_argument );

	meq::FluxSurfaceFamily ragged = family;
	ragged.surfaces[ 3 ].r.pop_back();
	BOOST_CHECK_THROW( meq::FluxGridWriter( file.path(), ragged ),
	                   std::invalid_argument );
}

// ---------------------------------------------------------------------------
// The cache, over the real extraction
// ---------------------------------------------------------------------------

/// The per-psi cache driven by the extraction it exists for.
///
/// tests/unit/FluxFamilyTests.cpp pins the contract against a synthetic
/// extractor, where "was this recomputed?" is exact. What this adds is that the
/// thing being cached is REPRODUCIBLE: a family extracted twice from the same
/// field has to be bit-identical, or "the cache is transparent" is not a
/// statement anyone can check. MaNTA pins the analogous property with a test
/// asserting that a reused coupled solver matches a fresh one bit for bit --
/// MANTA-COUPLING.md section 8 -- and the reason it is zero-tolerance there is
/// that the last defect of this kind left the second run completing, plausible,
/// and wrong in the eleventh digit.
///
/// THE PERTURBATION IS A UNIFORM SHIFT OF psi, WHICH IS A REAL FIELD AND NOT A
/// SCRIBBLE. Adding a constant leaves grad psi -- and so q -- exactly as it was,
/// so ( psi + delta, q ) is a consistent pair; what it moves is psi_ax and with
/// it every level the family is cut at. So a miss really does produce a
/// different family, which is what makes the extraction counts below mean
/// something.
BOOST_AUTO_TEST_CASE( theCacheServesTheFamilyTheExtractionWouldRebuild )
{
	SolvedEquilibrium solved( 2, 40 );
	meq::CriticalPoint const found = solved.magneticAxis();

	mfem::GridFunction &reference = solved.theSolver().postProcessedPotential();
	mfem::GridFunction const &flux = solved.theSolver().postProcessedFlux();

	std::vector<double> const psi( reference.GetData(),
	                               reference.GetData() + reference.Size() );

	int extractions = 0;
	auto extract = [ & ]( std::vector<double> const &state )
	{
		++extractions;

		mfem::GridFunction potential( reference.FESpace() );
		BOOST_TEST_REQUIRE( static_cast<int>( state.size() )
		                    == potential.Size() );
		for ( std::size_t i = 0; i < state.size(); ++i )
			potential( static_cast<int>( i ) ) = state[ i ];

		meq::ContourTracer const tracer( potential, flux );

		// psi_ax from THIS field: a uniform shift moves it, which is what makes
		// the perturbed family a different family.
		meq::CriticalPoint axisHere = found;
		double qR = 0.0;
		double qZ = 0.0;
		BOOST_TEST_REQUIRE( tracer.sampleAt( found.r, found.z, axisHere.psi,
		                                     qR, qZ ) );

		return meq::extractFluxSurfaces( tracer, axisHere, 0.0,
		                                 defaultOptions() );
	};

	// Twice from the same field, with no cache anywhere.
	meq::FluxSurfaceFamily const once = extract( psi );
	meq::FluxSurfaceFamily const twice = extract( psi );

	double worstRepeat = 0.0;
	BOOST_TEST_REQUIRE( once.size() == twice.size() );
	for ( std::size_t i = 0; i < once.size(); ++i )
	{
		worstRepeat = std::max( worstRepeat,
			std::abs( once.surfaces[ i ].vPrime - twice.surfaces[ i ].vPrime ) );
		for ( std::size_t j = 0; j < once.angles; ++j )
		{
			worstRepeat = std::max( worstRepeat,
				std::abs( once.surfaces[ i ].r[ j ]
				          - twice.surfaces[ i ].r[ j ] ) );
			worstRepeat = std::max( worstRepeat,
				std::abs( once.surfaces[ i ].z[ j ]
				          - twice.surfaces[ i ].z[ j ] ) );
		}
	}

	std::printf( "  the same field extracted twice: worst difference %.3e\n",
	             worstRepeat );
	BOOST_TEST( worstRepeat == 0.0,
	            "the extraction is not reproducible, so nothing downstream can "
	            "assert that a served answer equals a cold one" );

	// Now through the cache, in the consumer's own call pattern.
	meq::GeometryCache cache( extract );
	int const before = extractions;

	std::size_t const nodes = 40;
	std::size_t const residuals = 3;
	double worstServed = 0.0;

	for ( std::size_t residual = 0; residual < residuals; ++residual )
		for ( std::size_t node = 0; node < nodes; ++node )
		{
			double const label = once.innerLabel()
				+ ( once.outerLabel() - once.innerLabel() )
				  *static_cast<double>( node )/static_cast<double>( nodes );

			double const served = cache.at( psi, label,
				[]( meq::FluxSurface const &s ) { return s.vPrime; } );
			double const cold = once.vPrimeAt( label );
			worstServed = std::max( worstServed, std::abs( served - cold ) );
		}

	std::printf( "  %zu queries: %d extraction(s), %zu hits, served minus cold "
	             "%.3e\n",
	             cache.queries(), extractions - before, cache.hits(),
	             worstServed );

	BOOST_TEST( extractions - before == 1,
	            "the family was rebuilt inside the pointwise loop" );
	BOOST_TEST( cache.queries() == residuals*nodes );
	BOOST_TEST( cache.hits() == residuals*nodes - 1 );
	BOOST_TEST( worstServed == 0.0,
	            "a served answer differs from a cold one" );

	// resetForRun() rebuilds, and the rebuild is bit-identical because nothing
	// behind the extractor moved. That is the property MaNTA asserts at zero
	// tolerance.
	cache.resetForRun();
	BOOST_TEST( !cache.held() );
	double const afterReset = cache.at( psi, 0.5,
		[]( meq::FluxSurface const &s ) { return s.vPrime; } );
	BOOST_TEST( extractions - before == 2 );
	BOOST_TEST( afterReset == once.vPrimeAt( 0.5 ) );

	// A DIFFERENT STATE IS A DIFFERENT FAMILY, or the counts above would be
	// measuring nothing. A uniform shift of one per cent of the axis flux.
	std::vector<double> shifted = psi;
	double const delta = 0.01*std::abs( found.psi );
	for ( double &value : shifted )
		value += delta;

	double const moved = cache.at( shifted, 0.5,
		[]( meq::FluxSurface const &s ) { return s.vPrime; } );
	std::printf( "  V' at rho = 0.5: %.10e, shifted %.10e\n",
	             afterReset, moved );
	BOOST_TEST( extractions - before == 3 );
	BOOST_TEST( moved != afterReset,
	            "shifting psi did not move the family, so the extraction does "
	            "not depend on the state it is handed" );

	// And back again: the old state is NOT served from what the new one left.
	BOOST_TEST( cache.at( psi, 0.5,
		[]( meq::FluxSurface const &s ) { return s.vPrime; } ) == afterReset );
	BOOST_TEST( extractions - before == 4 );

	// WHAT IT IS WORTH, PRINTED AND NOT ASSERTED ON. IN-P measured the cache at
	// nodes / surfaces -- 5.1x on a sixty-node case -- and the count above is
	// the property; this is what the count costs in seconds on this machine,
	// which is a measurement about the machine and belongs in
	// tests/performance if it is ever to be tracked.
	meq::GeometryCache warm( extract );
	auto const startWarm = std::chrono::steady_clock::now();
	for ( std::size_t node = 0; node < nodes; ++node )
		warm.at( psi, 0.5, []( meq::FluxSurface const &s )
			{ return s.vPrime; } );
	auto const stopWarm = std::chrono::steady_clock::now();

	auto const startCold = std::chrono::steady_clock::now();
	for ( std::size_t node = 0; node < nodes; ++node )
	{
		meq::GeometryCache once2( extract );
		once2.at( psi, 0.5, []( meq::FluxSurface const &s )
			{ return s.vPrime; } );
	}
	auto const stopCold = std::chrono::steady_clock::now();

	// AND THE SHARPER BASELINE, WHICH IS THE ONE IN-P MEASURED. A consumer
	// written without a cache would not rebuild the whole family per node -- it
	// would locate the ONE surface through that node and integrate on it. That
	// is the honest naive, and against it the saving is nodes / surfaces rather
	// than nodes, which is why IN-P reads 5.1x on a sixty-node case and the
	// column above reads a much larger number against a much worse baseline.
	meq::ContourTracer const single( reference, flux );
	auto const startOne = std::chrono::steady_clock::now();
	for ( std::size_t node = 0; node < nodes; ++node )
	{
		double const label = once.innerLabel()
			+ ( once.outerLabel() - once.innerLabel() )
			  *static_cast<double>( node )/static_cast<double>( nodes );
		double const level = meq::fluxAtNormalised( label*label, found.psi, 0.0 );
		meq::AngleParametrisation const fit =
			single.fitByAngle( single.traceFromAxis( level, found ), found,
			                   defaultOptions().angles );
		// Computed and discarded: what is being timed is the extraction,
		// and the value is what the cached route above already asserted on.
		( void )meq::surfaceAverages( single, fit ).vPrime;
	}
	auto const stopOne = std::chrono::steady_clock::now();

	double const warmSeconds =
		std::chrono::duration<double>( stopWarm - startWarm ).count();
	double const coldSeconds =
		std::chrono::duration<double>( stopCold - startCold ).count();
	double const oneSeconds =
		std::chrono::duration<double>( stopOne - startOne ).count();
	std::printf( "  %zu nodes, %zu surfaces: cached %.4f s; naive family per "
	             "node %.4f s ( %.1fx ); naive SURFACE per node %.4f s "
	             "( %.1fx )\n",
	             nodes, once.size(), warmSeconds, coldSeconds,
	             coldSeconds/warmSeconds, oneSeconds,
	             oneSeconds/warmSeconds );
	std::fflush( stdout );
}
