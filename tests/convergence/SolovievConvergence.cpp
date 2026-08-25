#define BOOST_TEST_MODULE SolovievConvergence
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

#include "mfem.hpp"

#include "meq/GradShafranov.hpp"

#include "analytic/Soloviev.hpp"

/*
 * The stage-2 acceptance test: the linear HDG Grad-Shafranov operator measured
 * against an exact Solov'ev equilibrium.
 *
 * A convergence rate on its own cannot pass this test, and that is the point.
 * A wrong sign or a missing r weight converges at the right rate to the wrong
 * function, so the comparison is against a closed form -- the NSTX-like
 * equilibrium of refs/HDG-GradShafranov-Adaptive.pdf section 4.1, whose twelve
 * coefficients are published to fifteen digits -- and the absolute error is
 * asserted alongside the rate.
 *
 * refs/HDG-GradShafranov.pdf Tables 1-5 report k+1 for both psi and q. That is
 * what is asserted, with 0.15 of slack for the fact that a rate estimated from
 * two meshes is not the asymptotic one.
 *
 * The domain is a rectangle in ( r, z ) with r bounded well away from zero: the
 * Solov'ev expansion contains log r, and the operator carries 1/r. It is not the
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
 *   flux mass form holds 1/r instead of r    psi flat at 1.9e-2, rate 0.00
 *   potential r.h.s. is +F/r instead of -F/r psi flat at 7.3e-2, rate 0.00
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

	double const rMin = 0.6;
	double const rMax = 1.4;
	double const zMin = -0.6;
	double const zMax = 0.6;

	/// One point on the convergence curve.
	struct Measurement
	{
		double h;
		int traceDofs;
		double errorPsi;
		double errorFlux;
	};

	/// A triangulated rectangle [rMin,rMax] x [zMin,zMax] with n cells a side.
	/// Triangles rather than quadrilaterals because that is what both papers use
	/// and because it is measurably better here: on quadrilaterals the flux rate
	/// at k = 2 comes out at 2.78 rather than 3.00.
	mfem::Mesh makeMesh( int n )
	{
		mfem::Mesh mesh = mfem::Mesh::MakeCartesian2D( n, n, mfem::Element::TRIANGLE, false,
		                                               rMax - rMin, zMax - zMin );
		mesh.Transform( []( mfem::Vector const &in, mfem::Vector &out )
		{
			out( 0 ) = in( 0 ) + rMin;
			out( 1 ) = in( 1 ) + zMin;
		} );
		return mesh;
	}

	/// Solve once and measure. The source is F, not F/r: the solver applies the
	/// 1/r itself, and the Dirichlet datum is the exact psi on all four sides.
	Measurement measure( meq::analytic::SolovievEquilibrium const &eq, int order, int n )
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
			eq.flux( x( 0 ), x( 1 ), value( 0 ), value( 1 ) );
		} );

		meq::GradShafranovSolver solver( mesh, order );
		solver.setSource( sourceCoeff );
		solver.setBoundaryData( psiCoeff );
		solver.solve();

		Measurement point;
		point.h = ( rMax - rMin )/static_cast<double>( n );
		point.traceDofs = solver.numTraceDofs();
		point.errorPsi = solver.potentialError( psiCoeff );
		point.errorFlux = solver.fluxError( fluxCoeff );
		return point;
	}

	double rate( double coarseError, double fineError, double refinementRatio )
	{
		return std::log( coarseError/fineError )/std::log( refinementRatio );
	}

	/// The whole diagnostic. An assertion message says which number was wrong; a
	/// rate table says why, so it is printed unconditionally and ctest's
	/// --output-on-failure decides whether anyone sees it.
	void printTable( int order, std::vector<Measurement> const &points )
	{
		std::printf( "\n  Solov'ev NSTX, k = %d, triangles on [%.1f,%.1f]x[%.1f,%.1f]\n",
		             order, rMin, rMax, zMin, zMax );
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

	void checkOrder( int order, double psiCeiling, double fluxCeiling )
	{
		meq::analytic::SolovievEquilibrium const eq
			= meq::analytic::SolovievEquilibrium::nstx();

		std::vector<Measurement> points;
		points.reserve( meshSizes.size() );
		for ( int n : meshSizes )
			points.push_back( measure( eq, order, n ) );

		printTable( order, points );

		double const expected = order + 1.0 - rateSlack;

		for ( std::size_t i = 1; i < points.size(); ++i )
		{
			double ratio = points[ i - 1 ].h/points[ i ].h;
			double ratePsi = rate( points[ i - 1 ].errorPsi, points[ i ].errorPsi, ratio );
			double rateFlux = rate( points[ i - 1 ].errorFlux, points[ i ].errorFlux, ratio );

			BOOST_TEST( ratePsi >= expected,
			            "k = " << order << ", h = " << points[ i ].h
			            << ": psi converged at " << ratePsi << ", wanted " << expected );
			BOOST_TEST( rateFlux >= expected,
			            "k = " << order << ", h = " << points[ i ].h
			            << ": q converged at " << rateFlux << ", wanted " << expected );
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
		BOOST_TEST( points.back().errorPsi < psiCeiling,
		            "k = " << order << ": L2 error in psi is " << points.back().errorPsi
		            << ", above the ceiling " << psiCeiling );
		BOOST_TEST( points.back().errorFlux < fluxCeiling,
		            "k = " << order << ": L2 error in q is " << points.back().errorFlux
		            << ", above the ceiling " << fluxCeiling );
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

	for ( double r = rMin; r <= rMax + 1.0e-12; r += 0.2 )
	{
		for ( double z = zMin; z <= zMax + 1.0e-12; z += 0.3 )
		{
			double const deltaStar = eq.deltaStarFD( r, z );
			double const minusF = -eq.f( r, z, 0.0 );
			BOOST_TEST( std::abs( deltaStar - minusF ) < 1.0e-5,
			            "at ( " << r << ", " << z << " ): Delta*(psi) = " << deltaStar
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
		double const r = rMin + s*( rMax - rMin );
		double const z = zMin + s*( zMax - zMin );
		for ( double value : { eq.psi( r, zMin ), eq.psi( r, zMax ),
		                       eq.psi( rMin, z ), eq.psi( rMax, z ) } )
			largest = std::max( largest, std::abs( value ) );
	}

	BOOST_TEST( largest > 1.0e-2,
	            "the exact psi is only " << largest << " on the boundary of the test "
	            "rectangle, so the Dirichlet path is barely being exercised" );
}

BOOST_AUTO_TEST_CASE( orderOneConvergesAtTwo )
{
	// Ceilings at 3x the finest-mesh error measured with the corrected c_10:
	// psi 3.585e-05, q 5.971e-05 at h = 0.025.
	checkOrder( 1, 1.1e-4, 1.8e-4 );
}

BOOST_AUTO_TEST_CASE( orderTwoConvergesAtThree )
{
	// Measured: psi 1.895e-07, q 2.483e-07.
	checkOrder( 2, 5.7e-7, 7.5e-7 );
}

BOOST_AUTO_TEST_CASE( orderThreeConvergesAtFour )
{
	// Measured: psi 5.510e-10, q 8.349e-10.
	checkOrder( 3, 1.7e-9, 2.6e-9 );
}
