#define BOOST_TEST_MODULE SolovievConvergence
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

#include "mfem.hpp"

#include "RoundOffFloor.hpp"

#include "meq/GradShafranov.hpp"

#include "analytic/Soloviev.hpp"

/*
 * The stage-2 acceptance test: the linear HDG Grad-Shafranov operator measured
 * against an exact Solov'ev equilibrium.
 *
 * A convergence rate on its own cannot pass this test, and that is the point.
 * A wrong sign or a missing R weight converges at the right rate to the wrong
 * function, so the comparison is against a closed form -- the NSTX-like
 * equilibrium of refs/HDG-GradShafranov-Adaptive.pdf section 4.1, whose twelve
 * coefficients are published to fifteen digits -- and the absolute error is
 * asserted alongside the rate.
 *
 * refs/HDG-GradShafranov.pdf Tables 1-5 report k+1 for both psi and q. That is
 * what is asserted, with 0.15 of slack for the fact that a rate estimated from
 * two meshes is not the asymptotic one.
 *
 * The domain is a rectangle in ( R, z ) with R bounded well away from zero: the
 * Solov'ev expansion contains log R, and the operator carries 1/R. It is not the
 * plasma boundary, so psi is not zero on it -- which makes this a
 * non-homogeneous Dirichlet problem, and exercises a boundary-data path that a
 * homogeneous one would leave untested. Being a rectangle it is also fitted, so
 * Gamma_h == Gamma and the extension-from-subdomains machinery of stage 5 is not
 * involved.
 *
 * Rates stop improving around k ~ 5, or after enough refinements, when round-off
 * takes over; both papers report this. The table below stays well short: at k=3
 * on the finest mesh the error is still 5e-10.
 *
 * WHAT THE CONVENTIONS COST WHEN THEY ARE WRONG.
 *
 * Every choice in src/meq/GradShafranov.cpp was settled by running this test
 * with the alternative and reading the table, rather than by argument. The
 * numbers, so that they do not have to be rediscovered:
 *
 *   flux mass form holds 1/R instead of R    psi flat at 1.9e-2, rate 0.00
 *   potential R.h.s. is +F/R instead of -F/R psi flat at 7.3e-2, rate 0.00
 *   flux reported without the sign flip      q   flat at 6.5e-1, rate 0.00
 *   no boundary face integrator on B         psi flat at 1.5e-1, rate 0.01
 *   B's face integrator coefficients changed no change in any digit
 *   DarcyForm built with bsymmetrize = false psi flat at 3.1e-1 for all four
 *                                            sign combinations
 *
 * and the one that converges, at the right rate, to the wrong place:
 *
 *   HDGDiffusionIntegrator's built-in        psi at k+1, q at k only:
 *   h^-1 Q stabilisation, tau not overridden 1.03, 2.01, 3.00 for k = 1, 2, 3
 *
 * That last row is the argument for ConstantStabilization. It is also the one a
 * study of psi alone would have passed.
 */

namespace
{

	double const minRadius = 0.6;
	double const maxRadius = 1.4;
	double const zMin = -0.6;
	double const zMax = 0.6;

	/// One point on the convergence curve.
	struct Measurement
	{
		double h;
		int traceDofs;
		double errorPsi;
		double errorFlux;
		/// || psi || and || q || of the EXACT solution, in the same measure and
		/// on the same rule the two errors above are taken in. Zero means "no
		/// floor", which is what an unset one degrades to. M-218.
		double normPsi = 0.0;
		double normFlux = 0.0;
	};

	/// A triangulated rectangle [R_min,R_max] x [zMin,zMax] with n cells a side.
	/// Triangles rather than quadrilaterals because that is what both papers use
	/// and because it is measurably better here: on quadrilaterals the flux rate
	/// at k = 2 comes out at 2.78 rather than 3.00.
	mfem::Mesh makeMesh( int n )
	{
		mfem::Mesh mesh = mfem::Mesh::MakeCartesian2D( n, n, mfem::Element::TRIANGLE, false,
		                                               maxRadius - minRadius, zMax - zMin );
		mesh.Transform( []( mfem::Vector const &in, mfem::Vector &out )
		{
			out( 0 ) = in( 0 ) + minRadius;
			out( 1 ) = in( 1 ) + zMin;
		} );
		return mesh;
	}

	/// The stabilisation a sweep runs under: the solver's default, constant
	/// tau = 1, unless a case says otherwise.
	struct Stabilisation
	{
		double tau = 1.0;
		bool diffusionScaled = false;
	};

	/// Solve once and measure. The source is F, not F/R: the solver applies the
	/// 1/R itself, and the Dirichlet datum is the exact psi on all four sides.
	Measurement measure( meq::analytic::SolovievEquilibrium const &eq, int order, int n,
	                     Stabilisation const &stabilisation )
	{
		mfem::Mesh mesh = makeMesh( n );

		mfem::FunctionCoefficient sourceCoeff( [ &eq ]( mfem::Vector const &x )
		{
			return eq.f( x( 0 ), x( 1 ), 0.0 );
		} );
		mfem::FunctionCoefficient psiCoeff( [ &eq ]( mfem::Vector const &x )
		{
			return eq.psi( x( 0 ), x( 1 ) );
		} );
		mfem::VectorFunctionCoefficient fluxCoeff( 2, [ &eq ]( mfem::Vector const &x,
		                                                       mfem::Vector &value )
		{
			/*
		 * THE MFEM-FREE SEAM, WRITTEN OUT. `meq::analytic`'s `flux()` takes
		 * `double &` out-parameters because CLAUDE.md's *Layout* keeps the
		 * Source and Profile hierarchies free of MFEM -- so they stay `double`
		 * whatever `MFEM_PRECISION` says -- while `mfem::Vector::operator()`
		 * yields `real_t &`. In the double build those bind; in single they do
		 * not, and no cast can fix a REFERENCE. MEASUREMENTS.md M-217.
		 */
		double fluxR = 0.0, fluxZ = 0.0;
		eq.flux( x( 0 ), x( 1 ), fluxR, fluxZ );
		value( 0 ) = fluxR;
		value( 1 ) = fluxZ;
		} );

		meq::GradShafranovSolver solver( mesh, order, stabilisation.tau );
		if ( stabilisation.diffusionScaled )
			solver.setTauScaling( meq::GradShafranovSolver::TauScaling::Diffusion );
		solver.setSource( sourceCoeff );
		solver.setBoundaryData( psiCoeff );
		solver.solve();

		Measurement point;
		point.h = ( maxRadius - minRadius )/static_cast<double>( n );
		point.traceDofs = solver.numTraceDofs();
		point.errorPsi = solver.potentialError( psiCoeff );
		point.errorFlux = solver.fluxError( fluxCoeff );
		point.normPsi = solver.exactPotentialNorm( psiCoeff );
		point.normFlux = solver.exactFluxNorm( fluxCoeff );
		return point;
	}

	double rate( double coarseError, double fineError, double refinementRatio )
	{
		return std::log( coarseError/fineError )/std::log( refinementRatio );
	}

	/// The whole diagnostic. An assertion message says which number was wrong; a
	/// rate table says why, so it is printed unconditionally and ctest's
	/// --output-on-failure decides whether anyone sees it.
	void printTable( int order, std::vector<Measurement> const &points,
	                 Stabilisation const &stabilisation )
	{
		std::printf( "\n  Solov'ev NSTX, k = %d, triangles on [%.1f,%.1f]x[%.1f,%.1f], tau = %s%g\n",
		             order, minRadius, maxRadius, zMin, zMax,
		             stabilisation.diffusionScaled ? "tau0/R, tau0 = " : "", stabilisation.tau );
		std::printf( "  %8s %9s %14s %7s %14s %7s\n",
		             "h", "trace", "L2(psi)", "rate", "L2(q)", "rate" );
		for ( std::size_t i = 0; i < points.size(); ++i )
		{
			Measurement const &p = points[ i ];
			if ( i == 0 )
			{
				std::printf( "  %8.5f %9d %14.6e %7s %14.6e %7s\n",
				             p.h, p.traceDofs, p.errorPsi, "-", p.errorFlux, "-" );
			}
			else
			{
				double ratio = points[ i - 1 ].h/p.h;
				std::printf( "  %8.5f %9d %14.6e %7.3f %14.6e %7.3f\n",
				             p.h, p.traceDofs,
				             p.errorPsi, rate( points[ i - 1 ].errorPsi, p.errorPsi, ratio ),
				             p.errorFlux, rate( points[ i - 1 ].errorFlux, p.errorFlux, ratio ) );
			}
		}
		std::fflush( stdout );
	}

	/// The meshes: four dyadic refinements starting at 4 cells a side, so three
	/// measured rates per quantity per order. The whole suite runs in under a
	/// second, so this is not where to economise; it stops at 32 because that is
	/// already three orders of magnitude of error reduction at k = 3.
	std::vector<int> const meshSizes = { 4, 8, 16, 32 };

	/// k+1, less the slack allowed for a two-mesh rate estimate.
	double const rateSlack = 0.15;

	void checkOrder( int order, double psiCeiling, double fluxCeiling,
	                 Stabilisation const &stabilisation = Stabilisation() )
	{
		meq::analytic::SolovievEquilibrium const eq
			= meq::analytic::SolovievEquilibrium::nstx();

		std::vector<Measurement> points;
		points.reserve( meshSizes.size() );
		for ( int n : meshSizes )
			points.push_back( measure( eq, order, n, stabilisation ) );

		printTable( order, points, stabilisation );

		double const expected = order + 1.0 - rateSlack;

		// THE RATE IS NOT MEASURABLE AT THIS DEGREE IN THIS FORMAT. Skipped
		// loudly rather than loosened; the ceilings below still run. M-218.
		if ( !meq::tests::rateIsMeasurable( order ) )
		{
			BOOST_TEST_MESSAGE( meq::tests::rateNotMeasurableNotice( order ) );
		}
		else
		{
			// A pair whose finer error is under this build's round-off floor
			// measures the format and not the method; skip it, say so, and keep
			// a quorum so that a sequence entirely at the floor FAILS rather
			// than passing having asserted nothing. The floor cannot reach the
			// double build (roundOffFloor()), so this is inert there. Measured
			// in single at k = 1, the finest q error is 6.3e-05 or 6.7e-05
			// against a floor of 7.8e-05 depending on an ulp in R, and the pair
			// reads 1.933 or 1.846. M-322.
			int assertedPsi = 0;
			int assertedFlux = 0;
			for ( std::size_t i = 1; i < points.size(); ++i )
			{
				double ratio = points[ i - 1 ].h/points[ i ].h;
				double ratePsi = rate( points[ i - 1 ].errorPsi, points[ i ].errorPsi, ratio );
				double rateFlux = rate( points[ i - 1 ].errorFlux, points[ i ].errorFlux, ratio );

				if ( points[ i ].errorPsi < meq::tests::roundOffFloor( points[ i ].normPsi ) )
					BOOST_TEST_MESSAGE( "  *** k = " << order << ", h = " << points[ i ].h
					                    << ": SKIPPED the psi rate, the error is under the "
					                       "round-off floor. M-322" );
				else
				{
					++assertedPsi;
					BOOST_TEST( ratePsi >= expected,
					            "k = " << order << ", h = " << points[ i ].h
					            << ": psi converged at " << ratePsi << ", wanted " << expected );
				}
				if ( points[ i ].errorFlux < meq::tests::roundOffFloor( points[ i ].normFlux ) )
					BOOST_TEST_MESSAGE( "  *** k = " << order << ", h = " << points[ i ].h
					                    << ": SKIPPED the q rate, the error is under the "
					                       "round-off floor. M-322" );
				else
				{
					++assertedFlux;
					BOOST_TEST( rateFlux >= expected,
					            "k = " << order << ", h = " << points[ i ].h
					            << ": q converged at " << rateFlux << ", wanted " << expected );
				}
			}
			// A majority of the size() - 1 pairs.
			std::size_t const quorum = points.size()/2;
			BOOST_TEST( static_cast<std::size_t>( assertedPsi ) >= quorum,
			            "k = " << order << ": only " << assertedPsi
			            << " psi rates were above the round-off floor" );
			BOOST_TEST( static_cast<std::size_t>( assertedFlux ) >= quorum,
			            "k = " << order << ": only " << assertedFlux
			            << " q rates were above the round-off floor" );
		}

		// The rate is blind to a solution that is wrong by a constant factor or a
		// sign, so the absolute error is checked too. The ceilings sit at roughly
		// three times the measured values, which is loose enough not to be a
		// running cost and tight enough to catch a changed convention.
		//
		// They are also the ONLY thing in the suite that can see a wrong
		// coefficient c_i, because psi_1 ... psi_12 are Delta*-harmonic: any
		// c_i leaves F, Delta*(psi) and every convergence rate exact. So the
		// measured value each ceiling was set from is recorded beside it, and a
		// change to Soloviev.hpp that moves them should move these too,
		// deliberately, rather than passing because the ceiling had slack.
		// ONE IMPLEMENTATION, in ConvergenceHarness.hpp: these six lines stood
		// here and in three other files, and a floor added to one of four
		// copies is M-165 and M-189 repeating. M-218.
		meq::tests::checkCeilings( order,
		                           points.back().errorPsi, points.back().normPsi,
		                           psiCeiling,
		                           points.back().errorFlux, points.back().normFlux,
		                           fluxCeiling );
	}

}

/// The benchmark before the solver: Delta*( psi ) must equal -F, or everything
/// measured against it is measured against the wrong thing. Soloviev.hpp's own
/// comment argues this out from the two papers' disagreeing signs; this is the
/// arithmetic, on the mesh's own coordinate range.
BOOST_AUTO_TEST_CASE( solovievSourceMatchesTheOperator )
{
	meq::analytic::SolovievEquilibrium const eq
		= meq::analytic::SolovievEquilibrium::nstx();

	for ( double radius = minRadius; radius <= maxRadius + 1.0e-12; radius += 0.2 )
	{
		for ( double z = zMin; z <= zMax + 1.0e-12; z += 0.3 )
		{
			double const deltaStar = eq.deltaStarFD( radius, z );
			double const minusF = -eq.f( radius, z, 0.0 );
			BOOST_TEST( std::abs( deltaStar - minusF ) < 1.0e-5,
			            "at ( " << radius << ", " << z << " ): Delta*(psi) = " << deltaStar
			            << " but -F = " << minusF );
		}
	}
}

/// The Dirichlet data is not zero on this rectangle, which is what makes the
/// test exercise the non-homogeneous boundary path. Asserted rather than assumed,
/// because a benchmark that had quietly become homogeneous would still pass every
/// rate check below while testing less.
BOOST_AUTO_TEST_CASE( boundaryDataIsNonHomogeneous )
{
	meq::analytic::SolovievEquilibrium const eq
		= meq::analytic::SolovievEquilibrium::nstx();

	double largest = 0.0;
	for ( int i = 0; i <= 20; ++i )
	{
		double const s = static_cast<double>( i )/20.0;
		double const radius = minRadius + s*( maxRadius - minRadius );
		double const z = zMin + s*( zMax - zMin );
		for ( double value : { eq.psi( radius, zMin ), eq.psi( radius, zMax ),
		                       eq.psi( minRadius, z ), eq.psi( maxRadius, z ) } )
			largest = std::max( largest, std::abs( value ) );
	}

	BOOST_TEST( largest > 1.0e-2,
	            "the exact psi is only " << largest << " on the boundary of the test "
	            "rectangle, so the Dirichlet path is barely being exercised" );
}

BOOST_AUTO_TEST_CASE( orderOneConvergesAtTwo )
{
	// Ceilings at 3x the finest-mesh error, measured with the coefficient set
	// re-solved from Cerfon & Freidberg with alpha = arcsin(delta):
	// psi 3.597e-05, q 6.084e-05 at h = 0.025.
	checkOrder( 1, 1.1e-4, 1.9e-4 );
}

BOOST_AUTO_TEST_CASE( orderTwoConvergesAtThree )
{
	// Measured: psi 1.916e-07, q 2.510e-07.
	checkOrder( 2, 5.8e-7, 7.6e-7 );
}

BOOST_AUTO_TEST_CASE( orderThreeConvergesAtFour )
{
	// Measured: psi 5.510e-10, q 8.816e-10.
	checkOrder( 3, 1.7e-9, 2.7e-9 );
}

/*
 * THE DIFFUSION-SCALED tau KEEPS k+1 IN BOTH psi AND q, at tau0 = 3 m^-1 --
 * the value MEASUREMENTS.md M-289 recommends for it. Not a claim that
 * it is the better choice here: on this box nu = 1/R stays within [0.7, 1.7],
 * so tau = 3/R is two to five times the default everywhere, and the flux
 * pays for that in its constant (2.6x at k = 2), not in its order. What is
 * pinned is the order and the size of the answer, so that the option cannot
 * drift into MFEM's h^-1 LDG scaling -- which keeps psi at k+1 and loses an
 * order in q -- without this going red.
 */
BOOST_AUTO_TEST_CASE( theDiffusionScaledTauKeepsTheOrder )
{
	Stabilisation scaled;
	scaled.tau = 3.0;
	scaled.diffusionScaled = true;
	// Ceilings at 3x the measured finest-mesh errors, as above. Measured:
	// k = 1 psi 2.524e-05, q 8.858e-05; k = 2 psi 1.859e-07, q 6.554e-07;
	// k = 3 psi 4.607e-10, q 1.999e-09.
	checkOrder( 1, 7.6e-5, 2.7e-4, scaled );
	checkOrder( 2, 5.6e-7, 2.0e-6, scaled );
	checkOrder( 3, 1.4e-9, 6.0e-9, scaled );
}
