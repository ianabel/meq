#define BOOST_TEST_MODULE RotatingSolovievConvergence
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <vector>

#include "mfem.hpp"

#include "meq/GradShafranov.hpp"

#include "analytic/RotatingSoloviev.hpp"
#include "convergence/ConvergenceHarness.hpp"

/*
 * FLOW-PLAN.md stage FL-4: the HDG Grad-Shafranov operator measured against an
 * exact ROTATING Solov'ev equilibrium, refs/SpectralElementGSRotation.pdf
 * section 3.1 (Li & Zhu, CPC 260 (2021) 107264).
 *
 * WHAT IS NEW HERE, AND IT IS EXACTLY ONE THING. The source
 *
 *     F( r, z, psi ) = p1 r^2 exp[ M2 ( r^2/R0^2 - 1 ) ] + F0
 *
 * is still constant in psi -- so the problem is linear and dF/dpsi is
 * identically zero -- but it is EXPONENTIAL IN r^2 at fixed psi, which is the
 * whole structural consequence of sonic rotation: the pressure is no longer a
 * flux function, so mu0 dp/dpsi picks up an r. FLOW-PLAN.md section 1 records
 * that meq::Source's signature already carries r and so needs no change for
 * this; this file is where that claim stops being a claim.
 *
 * WHAT THIS TEST CANNOT SEE. dF/dpsi = 0, so it is the same rung as
 * SolovievConvergence.cpp and has the same blind spot: a Newton Jacobian is
 * neither exercised nor checked, because every Jacobian converges in one step
 * on an affine system. FLOW-PLAN.md section 6.4 requires a manufactured
 * nonlinear rotating case for that, and FL-4 comes before FL-5 deliberately --
 * a failure here is the discretisation or the source and CANNOT be the
 * Jacobian.
 *
 * ORDER OF THE TEST CASES IS LOAD BEARING. The Delta* scan comes first,
 * because everything below it is measured against a closed form that would
 * otherwise be taken on trust; a wrong transcription of eq (15), or a closure
 * that is not meq's, converges at the right rate to the wrong function and no
 * rate table in this file could tell. That is the standing hazard CLAUDE.md's
 * "Testing stance" is arranged around, and it has already fired twice in this
 * tree on published coefficients.
 *
 * The domain is meq::tests::standardBox(), the same rectangle
 * [0.6,1.4]x[-0.6,0.6] SolovievConvergence.cpp uses, so the two studies differ
 * in the source and in nothing else. It is not the plasma boundary -- the
 * fixture's psi = 0 contour is a closed curve strictly inside it -- so the
 * exact trace is imposed on all four sides and this is a non-homogeneous
 * Dirichlet problem, as there.
 */

namespace
{

	using Equilibrium = meq::analytic::RotatingSolovievEquilibrium;

	/// The scan SolovievConvergence.cpp uses, at the same step: 0.2 in r and
	/// 0.3 in z over the benchmark box.
	template<typename Check>
	void overTheBox( Check check )
	{
		meq::tests::Rectangle const box = meq::tests::standardBox();
		for ( double r = box.rMin; r <= box.rMax + 1.0e-12; r += 0.2 )
		{
			for ( double z = box.zMin; z <= box.zMax + 1.0e-12; z += 0.3 )
			{
				check( r, z );
			}
		}
	}

	/// The four Delta*-harmonic terms of eq (15) with the fixture's own
	/// coefficients -- the part of psi that the rotation does not touch. Used
	/// by the small-Mach test, which needs to compare the particular solution
	/// against eq (16) with the geometry divided out.
	double harmonicPart( Equilibrium const &eq, double r, double z )
	{
		std::array<double, 4> const c = eq.getCoefficients();
		double const r2 = r*r;
		double const z2 = z*z;
		return c[ 0 ] + c[ 1 ]*r2 + c[ 2 ]*( r2*r2 - 4.0*r2*z2 )
		       + c[ 3 ]*( r2*std::log( r ) - z2 );
	}

	/// eq (15)'s rotating term EXACTLY AS THE PAPER PRINTS IT, prefactor
	/// 1/M2^2 and all. Not what the fixture evaluates, and that is the point:
	/// it is the control for theSmallMachLimitIsContinuous, which asserts that
	/// this form loses everything for small machSquared while the fixture does
	/// not.
	double naiveRotatingTerm( double r, double machSquared, double p1,
	                          double majorRadius )
	{
		double const v = r*r/( majorRadius*majorRadius ) - 1.0;
		double const prefactor = majorRadius*majorRadius/( 2.0*machSquared );
		return -p1*prefactor*prefactor
		       *( std::exp( machSquared*v ) - machSquared*v - 1.0 );
	}

}

/// The benchmark before the solver: Delta*( psi ) must equal -F, or everything
/// measured against it is measured against the wrong thing.
///
/// THIS IS THE GUARD THE WHOLE FILE RESTS ON. It catches a mistyped term in
/// eq (15), a sign slip between Li & Zhu's Delta* and meq's, and -- the reason
/// FLOW-PLAN.md section 6.2 insists on it -- a closure that is not (136)'s
/// isothermal one at all. Without it a wrong fixture would produce a perfect
/// k+1 table for somebody else's equilibrium.
///
/// Run at machSquared = 0 and 1 only. fastRotating() is excluded because the
/// CENTRAL DIFFERENCE, not the fixture, gives out there: see the note on that
/// factory in RotatingSoloviev.hpp.
BOOST_AUTO_TEST_CASE( theRotatingSolovievSourceMatchesTheOperator )
{
	for ( Equilibrium const &eq : { Equilibrium::stationary(), Equilibrium::rotating() } )
	{
		double worst = 0.0;
		overTheBox( [ &eq, &worst ]( double r, double z )
		{
			double const deltaStar = eq.deltaStarFD( r, z );
			double const minusF = -eq.f( r, z, 0.0 );
			worst = std::max( worst, std::abs( deltaStar - minusF ) );

			BOOST_TEST( std::abs( deltaStar - minusF ) < 1.0e-5,
			            "at machSquared = " << eq.getMachSquared() << ", ( " << r
			            << ", " << z << " ): Delta*(psi) = " << deltaStar
			            << " but -F = " << minusF );
		} );

		// Measured: 1.615e-08 at machSquared = 0 and 2.768e-07 at machSquared
		// = 1, both of which are the h^2 truncation floor of the difference and
		// not a property of the fixture.
		std::printf( "  Delta* check at machSquared = %.1f: worst %10.3e\n",
		             eq.getMachSquared(), worst );
	}
	std::fflush( stdout );
}

/// gradPsi() is differentiated by hand and the solver's flux error is measured
/// against it, so a slip there would move every L2( q ) in this file without
/// moving a single rate. deltaStarFD() cannot see it: that is built on psi()
/// alone, deliberately.
BOOST_AUTO_TEST_CASE( theGradientsMatchFiniteDifferences )
{
	double const h = 1.0e-6;

	for ( Equilibrium const &eq : { Equilibrium::stationary(), Equilibrium::rotating(),
	                                Equilibrium::fastRotating() } )
	{
		overTheBox( [ &eq, h ]( double r, double z )
		{
			double analyticR, analyticZ;
			eq.gradPsi( r, z, analyticR, analyticZ );

			double const differencedR = ( eq.psi( r + h, z ) - eq.psi( r - h, z ) )
			                            /( 2.0*h );
			double const differencedZ = ( eq.psi( r, z + h ) - eq.psi( r, z - h ) )
			                            /( 2.0*h );

			// Measured worst over the three configurations: 4.0e-10, which is
			// the O( h^2 ) truncation of the difference.
			BOOST_TEST( std::abs( analyticR - differencedR ) < 1.0e-8,
			            "at machSquared = " << eq.getMachSquared() << ", ( " << r
			            << ", " << z << " ): d_r psi = " << analyticR
			            << " but the difference is " << differencedR );
			BOOST_TEST( std::abs( analyticZ - differencedZ ) < 1.0e-8,
			            "at machSquared = " << eq.getMachSquared() << ", ( " << r
			            << ", " << z << " ): d_z psi = " << analyticZ
			            << " but the difference is " << differencedZ );
		} );

		// And the flux really is grad_bar( psi )/r, since that is what the
		// solver's fluxError() is handed.
		double qR, qZ, gR, gZ;
		eq.flux( 1.1, 0.2, qR, qZ );
		eq.gradPsi( 1.1, 0.2, gR, gZ );
		BOOST_TEST( std::abs( qR - gR/1.1 ) < 1.0e-15 );
		BOOST_TEST( std::abs( qZ - gZ/1.1 ) < 1.0e-15 );
	}
}

/// eq (15) as printed is 0/0 as machSquared -> 0 -- a 1/M2^2 prefactor on a
/// brace that vanishes like M2^2 -- and the fixture evaluates an algebraically
/// equivalent form in which M2 has cancelled. This is the test of that, in
/// three parts, and the middle one is the control that makes the other two
/// mean something.
BOOST_AUTO_TEST_CASE( theSmallMachLimitIsContinuous )
{
	Equilibrium const stationary = Equilibrium::stationary();
	std::array<double, 4> const c = stationary.getCoefficients();
	double const p1 = stationary.getP1();
	double const f0 = stationary.getF0();
	double const r0 = stationary.getMajorRadius();

	// PART 1: at machSquared = 0 the fixture must BE the static Solov'ev
	// particular solution of eq (16), psi_h - p1( r^2 - R0^2 )^2/8 - ( F0/2 )z^2,
	// which is a different expression and not merely a limit of one.
	// Measured: 0.0e+00 over the scan this test runs, and 5.6e-17 over a
	// 0.05-stepped one.
	{
		double worst = 0.0;
		overTheBox( [ & ]( double r, double z )
		{
			double const reference = harmonicPart( stationary, r, z )
			                         - p1*std::pow( r*r - r0*r0, 2 )/8.0
			                         - 0.5*f0*z*z;
			worst = std::max( worst, std::abs( stationary.psi( r, z ) - reference ) );
		} );
		BOOST_TEST( worst < 1.0e-15,
		            "at machSquared = 0 the fixture differs from eq (16) by "
		            << worst );
		std::printf( "  eq (15) at machSquared = 0 against eq (16): %10.3e\n", worst );
	}

	// PART 2, THE CONTROL. The naive form is not merely less accurate for small
	// machSquared, it is wrong -- and a test that only checked the fixture
	// against itself would pass with the naive form in place. MEASURED, at
	// r = 1.4 where v = 0.96, comparing the naive rotating term against the
	// fixture's AT THE SAME machSquared:
	//
	//     machSquared    naive, relative error against the fixture
	//        1e-1               1.1e-14
	//        1e-2               2.7e-12
	//        1e-4               8.8e-09
	//        1e-6               3.6e-04
	//        1e-8               3.8e+00      <-- 380 per cent
	//        1e-10              1.0e+00      <-- returns exactly zero
	//
	// Asserted at 1e-4 and 1e-8 together, so the test states both halves: the
	// naive form is fine where there is no cancellation and useless where there
	// is. The first half matters as much as the second -- it is what says the
	// two expressions really are the same algebra, so that the second half is
	// measuring precision loss and not a transcription error in either.
	//
	// The comparison is at z = 0, where psi minus its harmonic part IS the
	// rotating term, the -( F0/2 )z^2 piece having dropped out.
	{
		double const r = 1.4;
		double const v = r*r/( r0*r0 ) - 1.0;
		double const limit = -p1*r0*r0*r0*r0*v*v/8.0;

		auto rotatingTermOf = [ & ]( Equilibrium const &eq )
		{
			return eq.psi( r, 0.0 ) - harmonicPart( eq, r, 0.0 );
		};

		double const mildMach = 1.0e-4;
		Equilibrium const mild( r0, mildMach, p1, f0, c );
		double const stableMild = rotatingTermOf( mild );
		double const naiveMild = naiveRotatingTerm( r, mildMach, p1, r0 );
		BOOST_TEST( std::abs( naiveMild - stableMild )/std::abs( stableMild ) < 1.0e-6,
		            "the two forms disagree by "
		            << std::abs( naiveMild - stableMild )/std::abs( stableMild )
		            << " at machSquared = " << mildMach
		            << ", where there is no cancellation to blame -- so they are "
		               "not the same algebra and one of them is mistyped" );

		double const smallMach = 1.0e-8;
		Equilibrium const small( r0, smallMach, p1, f0, c );
		double const stableSmall = rotatingTermOf( small );
		double const naiveSmall = naiveRotatingTerm( r, smallMach, p1, r0 );
		BOOST_TEST( std::abs( naiveSmall - stableSmall )/std::abs( stableSmall ) > 1.0e-2,
		            "eq (15) as printed is accurate to "
		            << std::abs( naiveSmall - stableSmall )/std::abs( stableSmall )
		            << " at machSquared = " << smallMach
		            << ", so the stable form is not buying anything and this "
		               "control has stopped being a control" );

		// And it is the naive one that is wrong, not the fixture: at
		// machSquared = 1e-8 the true rotating term is within 1e-8 of its
		// machSquared = 0 limit, and the fixture is. Measured: 3.2e-09.
		BOOST_TEST( std::abs( stableSmall - limit )/std::abs( limit ) < 1.0e-7,
		            "the fixture is out by "
		            << std::abs( stableSmall - limit )/std::abs( limit )
		            << " at machSquared = " << smallMach );
	}

	// PART 3: the fixture converges to its own machSquared = 0 value linearly
	// in machSquared -- which is what d/dmachSquared being finite means -- and
	// does so all the way down to 1e-12, where the naive form has no digits
	// left at all. Both gradient components go the same way.
	{
		for ( double mach : { 1.0e-2, 1.0e-4, 1.0e-6, 1.0e-8, 1.0e-10, 1.0e-12 } )
		{
			Equilibrium const eq( r0, mach, p1, f0, c );
			overTheBox( [ & ]( double r, double z )
			{
				double referenceR, referenceZ, valueR, valueZ;
				stationary.gradPsi( r, z, referenceR, referenceZ );
				eq.gradPsi( r, z, valueR, valueZ );

				// The slope is -( p1 R0^4/4 ) v^3/6 in psi, so at most 0.04 on
				// this box; 1.0 is a bound on it and not a fitted constant.
				BOOST_TEST( std::abs( eq.psi( r, z ) - stationary.psi( r, z ) )
				            < 1.0*mach + 1.0e-15,
				            "psi at machSquared = " << mach << " is "
				            << eq.psi( r, z ) << " against " << stationary.psi( r, z )
				            << " at ( " << r << ", " << z << " )" );
				BOOST_TEST( std::abs( valueR - referenceR ) < 1.0*mach + 1.0e-15 );
				BOOST_TEST( std::abs( valueZ - referenceZ ) < 1.0e-15 );
			} );
		}
	}

	// PART 4: the crossover itself. g1() and g2() change branch at
	// |u| = seriesThreshold(), and a series that disagreed with the closed form
	// there would put a step into psi and into the flux -- small enough to
	// survive every other assertion in this file and large enough to spoil a
	// convergence rate. Straddle it in machSquared at fixed r and require the
	// jump to scale like the perturbation, which is what "no step" means.
	//
	// Measured at r = 1.4, z = 0.3, where the crossover is machSquared =
	// 0.5/0.96 = 0.520833:
	//
	//     delta      jump in psi      jump in d_r psi
	//      1e-6         9.5e-08           9.1e-07
	//      1e-9         9.5e-11           9.1e-10
	//      1e-12        9.5e-14           9.1e-13
	//
	// i.e. exactly the local slope times 2 delta at every scale, with no floor
	// underneath it -- there is no step.
	{
		double const r = 1.4;
		double const z = 0.3;
		double const v = r*r/( r0*r0 ) - 1.0;
		double const crossover = Equilibrium::seriesThreshold()/v;

		for ( double delta : { 1.0e-6, 1.0e-9, 1.0e-12 } )
		{
			Equilibrium const below( r0, crossover - delta, p1, f0, c );
			Equilibrium const above( r0, crossover + delta, p1, f0, c );

			double belowR, belowZ, aboveR, aboveZ;
			below.gradPsi( r, z, belowR, belowZ );
			above.gradPsi( r, z, aboveR, aboveZ );

			double const jumpPsi = std::abs( above.psi( r, z ) - below.psi( r, z ) );
			double const jumpFlux = std::abs( aboveR - belowR );

			// 10 delta is two orders above the measured 0.1 and 0.9 delta, and
			// is a bound on the slope rather than a fitted tolerance. A genuine
			// branch mismatch does not scale with delta at all, so it fails this
			// at the smallest one however loose the constant is.
			BOOST_TEST( jumpPsi < 10.0*delta,
			            "psi jumps by " << jumpPsi << " across the series "
			            "crossover at machSquared = " << crossover
			            << " for a perturbation of " << delta );
			BOOST_TEST( jumpFlux < 10.0*delta,
			            "d_r psi jumps by " << jumpFlux << " across the series "
			            "crossover at machSquared = " << crossover
			            << " for a perturbation of " << delta );
		}
	}

	// And the crossover is reached the other way too, at fixed machSquared with
	// r sweeping through R0, where u -> 0 for any Mach number whatever. psi is
	// even in v to leading order there, so what this catches is a series that
	// is wrong rather than merely truncated.
	{
		Equilibrium const eq = Equilibrium::rotating();
		double previous = 0.0;
		bool first = true;
		for ( int i = -200; i <= 200; ++i )
		{
			double const r = r0 + i*1.0e-3;
			double const value = eq.psi( r, 0.0 );
			if ( !first )
			{
				BOOST_TEST( std::abs( value - previous ) < 1.0e-3,
				            "psi steps by " << std::abs( value - previous )
				            << " between r = " << r - 1.0e-3 << " and " << r );
			}
			previous = value;
			first = false;
		}
	}
}

/// The fixture's own coefficients, checked against the four geometric
/// conditions they were solved from.
///
/// This is the ONLY assertion in the file that can see a mistyped c_i, for the
/// same reason Soloviev.hpp records: all four terms are Delta*-harmonic, so a
/// wrong coefficient leaves F, Delta*( psi ) and every convergence rate exact
/// and changes only which equilibrium is being solved. The absolute-error
/// ceilings below would eventually notice, but they would notice as an
/// unexplained ceiling failure rather than as this.
BOOST_AUTO_TEST_CASE( theCoefficientsPutTheZeroContourWhereItWasDesigned )
{
	for ( Equilibrium const &eq : { Equilibrium::stationary(), Equilibrium::rotating(),
	                                Equilibrium::fastRotating() } )
	{
		double gradR, gradZ;
		eq.gradPsi( 1.0, 0.45, gradR, gradZ );

		// Measured: every one of these is at or below 7e-17.
		BOOST_TEST( std::abs( eq.psi( 1.3, 0.0 ) ) < 1.0e-15,
		            "machSquared = " << eq.getMachSquared()
		            << ": psi at the outer equatorial point is "
		            << eq.psi( 1.3, 0.0 ) );
		BOOST_TEST( std::abs( eq.psi( 0.7, 0.0 ) ) < 1.0e-15,
		            "machSquared = " << eq.getMachSquared()
		            << ": psi at the inner equatorial point is "
		            << eq.psi( 0.7, 0.0 ) );
		BOOST_TEST( std::abs( eq.psi( 1.0, 0.45 ) ) < 1.0e-15,
		            "machSquared = " << eq.getMachSquared()
		            << ": psi at the high point is " << eq.psi( 1.0, 0.45 ) );
		BOOST_TEST( std::abs( gradR ) < 1.0e-14,
		            "machSquared = " << eq.getMachSquared()
		            << ": d_r psi at the high point is " << gradR );

		// psi is strictly negative on the whole boundary of the benchmark box
		// and positive at the axis, so psi = 0 is a closed curve inside it.
		// That is what makes this a physically sensible equilibrium rather than
		// merely a solution of the equation, and it is also what makes the
		// Dirichlet data non-homogeneous -- the boundary path this study
		// exercises and a homogeneous one would not.
		meq::tests::Rectangle const box = meq::tests::standardBox();
		double worstOnBoundary = -1.0e300;
		for ( int i = 0; i <= 200; ++i )
		{
			double const s = static_cast<double>( i )/200.0;
			double const r = box.rMin + s*box.width();
			double const z = box.zMin + s*box.height();
			for ( double value : { eq.psi( r, box.zMin ), eq.psi( r, box.zMax ),
			                       eq.psi( box.rMin, z ), eq.psi( box.rMax, z ) } )
				worstOnBoundary = std::max( worstOnBoundary, value );
		}

		// Measured: -0.0346, -0.0327 and -0.0368 for the three configurations.
		BOOST_TEST( worstOnBoundary < -1.0e-2,
		            "machSquared = " << eq.getMachSquared() << ": psi reaches "
		            << worstOnBoundary << " on the boundary of the benchmark box, "
		            "so the psi = 0 contour is not closed inside it" );
		BOOST_TEST( eq.psi( 1.0, 0.0 ) > 1.0e-2,
		            "machSquared = " << eq.getMachSquared()
		            << ": psi on the axis is only " << eq.psi( 1.0, 0.0 ) );
	}
}

/*
 * THE RATES. meq::tests::checkOrder asserts k+1 in psi AND in q over the four
 * dyadic meshes, plus an absolute ceiling on the finest -- which is the only
 * assertion here that can see a solution converging beautifully to the wrong
 * function.
 *
 * The slack is 0.15 rather than the harness default of 0.2, matching
 * SolovievConvergence.cpp: this problem is linear and has an exact solution, so
 * there is no reason to allow the width the non-linear studies need.
 *
 * The source reaches the solver as a meq::Source rather than as a bare
 * mfem::Coefficient, so these run down the NEWTON path even though dF/dpsi is
 * zero. That is deliberate and free: the Newton column of the printed table
 * must read 1 at every mesh and every order, which is the McCarthy ladder's
 * bottom rung asserted in passing -- an affine system whose first step is not
 * exact means the Jacobian's mass term has appeared from somewhere it should
 * not have.
 */

BOOST_AUTO_TEST_CASE( orderOneConvergesAtTwo )
{
	// Ceilings at roughly 3x the measured finest-mesh error:
	// psi 8.308524e-05, q 1.369991e-04 at h = 0.025.
	meq::tests::checkOrder( Equilibrium::rotating(), "rotating Solov'ev, M2 = 1", 1,
	                        2.5e-4, 4.1e-4, meq::tests::standardBox(),
	                        meq::tests::dyadicMeshes(), 0.15 );
}

BOOST_AUTO_TEST_CASE( orderTwoConvergesAtThree )
{
	// Measured: psi 6.798255e-07, q 1.056370e-06.
	meq::tests::checkOrder( Equilibrium::rotating(), "rotating Solov'ev, M2 = 1", 2,
	                        2.0e-6, 3.2e-6, meq::tests::standardBox(),
	                        meq::tests::dyadicMeshes(), 0.15 );
}

BOOST_AUTO_TEST_CASE( orderThreeConvergesAtFour )
{
	// Measured: psi 3.981334e-09, q 6.436949e-09.
	meq::tests::checkOrder( Equilibrium::rotating(), "rotating Solov'ev, M2 = 1", 3,
	                        1.2e-8, 1.9e-8, meq::tests::standardBox(),
	                        meq::tests::dyadicMeshes(), 0.15 );
}
