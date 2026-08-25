#define BOOST_TEST_MODULE SolovievGeometryConvergence
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "mfem.hpp"

#include "meq/GradShafranov.hpp"

#include "analytic/Soloviev.hpp"
#include "convergence/ConvergenceHarness.hpp"

/*
 * The three Solov'ev geometries of refs/HDG-GradShafranov.pdf section 4.2 --
 * its Examples 1, 2 and 3 -- as exact-solution convergence benchmarks.
 *
 *   Example 1   FRC          eps = 0.99, delta = 0.7,   kappa = 10,  A = 0
 *   Example 2   ITER-like    eps = 0.32, delta = 0.33,  kappa = 2,   A = -0.115
 *   Example 3   NSTX-like    eps = 0.78, delta = 0.335, kappa = 1.7, A = -0.115
 *
 * Example 1 is up-down symmetric and smooth, so it is Cerfon & Freidberg's
 * seven-coefficient expansion closed by their seven smooth-surface constraints.
 * Examples 2 and 3 are "up-down asymmetric with a downwards oriented x-point",
 * in the paper's own words, so they are the twelve-coefficient expansion closed
 * by the twelve single-null constraints. tests/analytic/Soloviev.hpp carries the
 * derivations and the verification record.
 *
 * WHY THESE ARE WORTH HAVING, given that SolovievConvergence.cpp already
 * measures k+1 against a Solov'ev solution. Two reasons, and only the second is
 * about the solver.
 *
 * The first is that this file is where the COEFFICIENTS are finally asserted.
 * ExtensionConvergence.cpp records the gap in as many words: "the coefficients
 * are the one part of the benchmark that is asserted nowhere", because every
 * psi_i is Delta*-harmonic and so any coefficients whatever leave the source,
 * the operator and every convergence rate exact. What they do decide is the
 * geometry, and the geometry IS checkable -- against Cerfon & Freidberg's own
 * boundary constraints. coefficientsSatisfyTheGeometricConstraints below does
 * that for all four sets in the fixture, including the two that were there
 * before, and it is the only test in the tree that can see a mistyped c_i.
 *
 * The second is range. The existing study runs one configuration, at
 * A = -0.52 on a mesh whose flux varies by a factor of a few. These three span
 * A = 0 to A = -0.115, aspect ratios from 0.32 to 0.99, and elongations from
 * 1.7 to 10, so the balance between the r^4, the r^6 and the log r terms of the
 * expansion is quite different in each -- and the log r terms are where a weight
 * convention on the 1/r in the operator would show up.
 *
 * THE DOMAIN. All three are posed on the same rectangle as
 * SolovievConvergence.cpp and NewtonConvergence.cpp, so the error levels are
 * directly comparable with those tables. It is not any of the three plasma
 * boundaries -- fitting a curved separatrix is the extension path, measured in
 * ExtensionConvergence.cpp -- which means the Dirichlet data is non-homogeneous,
 * and that is the point: it is the exact solution restricted to the box, so the
 * solution of the discrete problem converges to a function that is known in
 * closed form. Note that the box pokes slightly outside Example 2's plasma in r;
 * the expansion is analytic for every r > 0, so that is harmless, and keeping
 * one box for all three is worth more than staying inside each separatrix.
 *
 * These sources do not depend on psi, so dF/dpsi is identically zero and the
 * Newton path must finish in a single step. That is asserted: it is the same
 * statement NewtonConvergence.cpp makes about the A = -0.52 case, and it is what
 * makes these tables a check on the discretisation alone.
 */

namespace
{
	using meq::analytic::SolovievEquilibrium;
	using meq::tests::standardBox;

	/// The geometric parameters of one configuration, as the paper states them.
	struct Geometry
	{
		char const *name;
		double eps, delta, kappa;
		bool asymmetric;
	};

	Geometry const frc = { "Example 1 FRC", 0.99, 0.7, 10.0, false };
	Geometry const iter = { "Example 2 ITER-like", 0.32, 0.33, 2.0, true };
	Geometry const nstxThree = { "Example 3 NSTX-like", 0.78, 0.335, 1.7, true };
	/// refs/HDG-GradShafranov-Adaptive.pdf section 4.1, for comparison.
	Geometry const nstxAdaptive = { "Adaptive 4.1 NSTX", 0.78, 0.35, 2.0, true };

	/// Cerfon & Freidberg eq (11), with sin( alpha ) = delta as their section III
	/// defines it. Verified independently by differentiating their model surface
	/// eq (9): see the note in tests/analytic/Soloviev.hpp.
	struct Curvatures
	{
		double n1, n2, n3;
	};

	Curvatures curvatures( Geometry const &g, bool alphaIsDelta = false )
	{
		double const alpha = alphaIsDelta ? g.delta : std::asin( g.delta );
		Curvatures c;
		c.n1 = -( 1.0 + alpha )*( 1.0 + alpha )/( g.eps*g.kappa*g.kappa );
		c.n2 = ( 1.0 - alpha )*( 1.0 - alpha )/( g.eps*g.kappa*g.kappa );
		c.n3 = -g.kappa/( g.eps*std::cos( alpha )*std::cos( alpha ) );
		return c;
	}

	/// d_zz psi, by a central difference of the analytic d_z psi. Second order,
	/// so good to about 1e-9 here -- six orders of magnitude tighter than the
	/// discrepancy this is used to detect.
	double psiZZ( SolovievEquilibrium const &eq, double r, double z,
	              double h = 1.0e-4 )
	{
		double gr, up, down;
		eq.gradPsi( r, z + h, gr, up );
		eq.gradPsi( r, z - h, gr, down );
		return ( up - down )/( 2.0*h );
	}

	/// d_rr psi, likewise.
	double psiRR( SolovievEquilibrium const &eq, double r, double z,
	              double h = 1.0e-4 )
	{
		double gz, up, down;
		eq.gradPsi( r + h, z, up, gz );
		eq.gradPsi( r - h, z, down, gz );
		return ( up - down )/( 2.0*h );
	}

	double psiR( SolovievEquilibrium const &eq, double r, double z )
	{
		double gr, gz;
		eq.gradPsi( r, z, gr, gz );
		return gr;
	}

	double psiZ( SolovievEquilibrium const &eq, double r, double z )
	{
		double gr, gz;
		eq.gradPsi( r, z, gr, gz );
		return gz;
	}
}

/*
 * -------------------------------------------------------------------------
 * The coefficients, asserted at last
 * -------------------------------------------------------------------------
 */

/// Cerfon & Freidberg's boundary constraints, evaluated on every coefficient set
/// in the fixture. This is the check that no convergence rate can make: the
/// twelve psi_i are Delta*-harmonic, so the coefficients are invisible to the
/// equation and visible only in the geometry.
///
/// The conditions, from refs/CerfonFreidberg.pdf:
///
///   eq (10), up-down symmetric and smooth -- seven of them:
///     psi = 0 at ( 1+eps, 0 ), ( 1-eps, 0 ) and ( 1-delta eps, kappa eps )
///     psi_r = 0 at the high point
///     psi_zz = -N1 psi_r at the outer equatorial point
///     psi_zz = -N2 psi_r at the inner equatorial point
///     psi_rr = -N3 psi_z at the high point
///
///   eq (28), up-down asymmetric single null -- twelve: the above, plus
///     psi_z = 0 at both equatorial points, and
///     psi = psi_r = psi_z = 0 at the X-point ( 1-1.1 delta eps, -1.1 kappa eps )
BOOST_AUTO_TEST_CASE( coefficientsSatisfyTheGeometricConstraints )
{
	struct Case
	{
		Geometry geometry;
		SolovievEquilibrium equilibrium;
	};

	std::vector<Case> const cases = {
		{ frc, SolovievEquilibrium::frcExample1() },
		{ iter, SolovievEquilibrium::iterExample2() },
		{ nstxThree, SolovievEquilibrium::nstxExample3() },
		{ nstxAdaptive, SolovievEquilibrium::nstx() }
	};

	for ( Case const &one : cases )
	{
		Geometry const &g = one.geometry;
		SolovievEquilibrium const &eq = one.equilibrium;
		Curvatures const n = curvatures( g );

		double const outer = 1.0 + g.eps;
		double const inner = 1.0 - g.eps;
		double const highR = 1.0 - g.delta*g.eps;
		double const highZ = g.kappa*g.eps;
		double const sepR = 1.0 - 1.1*g.delta*g.eps;
		double const sepZ = -1.1*g.kappa*g.eps;

		// The flux at the magnetic axis, as the scale every residual is measured
		// against: a condition satisfied to 1e-17 absolute means nothing without
		// knowing that psi itself is O( 0.1 ) here.
		double const scale = std::abs( eq.psi( 1.0, 0.0 ) );

		struct Check
		{
			char const *what;
			double value;
		};

		std::vector<Check> checks = {
			{ "psi( outer equatorial )", eq.psi( outer, 0.0 ) },
			{ "psi( inner equatorial )", eq.psi( inner, 0.0 ) },
			{ "psi( high point )", eq.psi( highR, highZ ) },
			{ "psi_r( high point )", psiR( eq, highR, highZ ) },
			{ "curvature, outer", psiZZ( eq, outer, 0.0 ) + n.n1*psiR( eq, outer, 0.0 ) },
			{ "curvature, inner", psiZZ( eq, inner, 0.0 ) + n.n2*psiR( eq, inner, 0.0 ) },
			{ "curvature, high", psiRR( eq, highR, highZ )
			                     + n.n3*psiZ( eq, highR, highZ ) }
		};

		if ( g.asymmetric )
		{
			checks.push_back( { "psi_z( outer equatorial )", psiZ( eq, outer, 0.0 ) } );
			checks.push_back( { "psi_z( inner equatorial )", psiZ( eq, inner, 0.0 ) } );
			checks.push_back( { "psi( X-point )", eq.psi( sepR, sepZ ) } );
			checks.push_back( { "psi_r( X-point )", psiR( eq, sepR, sepZ ) } );
			checks.push_back( { "psi_z( X-point )", psiZ( eq, sepR, sepZ ) } );
		}

		std::printf( "\n  %s: %zu Cerfon-Freidberg conditions, psi( 1, 0 ) = %.4e\n",
		             g.name, checks.size(), eq.psi( 1.0, 0.0 ) );
		for ( Check const &check : checks )
		{
			std::printf( "    %-26s %12.4e  (%.2e of the axis flux)\n",
			             check.what, check.value, std::abs( check.value )/scale );
		}
		std::fflush( stdout );

		// The curvature conditions get 1e-6 relative and the rest 1e-13. That is
		// not slack, it is the finite-difference floor: psi_zz and psi_rr come
		// from a second-order central difference of the analytic gradient at
		// h = 1e-4, so their truncation error is O( h^2 ) times a fourth
		// derivative, which measures 1e-8 to 1e-7 of the axis flux here. The
		// defect this has to be able to see -- the wrong alpha in nstx() -- is
		// 4.4e-2 of the axis flux, four orders above the tolerance.
		for ( Check const &check : checks )
		{
			bool const isCurvature = std::string( check.what ).substr( 0, 9 ) == "curvature";
			double const tolerance = isCurvature ? 1.0e-6 : 1.0e-13;

			BOOST_TEST( std::abs( check.value )/scale < tolerance,
			            g.name << ": " << check.what << " = " << check.value
			            << ", which is " << std::abs( check.value )/scale
			            << " of the axis flux. A coefficient is wrong, and nothing "
			            "else in this suite can see it" );
		}
	}
}

/// nstx() must satisfy Cerfon & Freidberg's two equatorial curvature conditions
/// with alpha = arcsin( delta ), which is what their eq (11) means.
///
/// This test exists because it did not, twice. The published coefficients were
/// wrong throughout; the replacement solved here was then wrong again, more
/// subtly, because eq (11)'s alpha was read as sin( alpha ) = delta in N_1 and
/// N_2. N_3 is immune to that confusion, so the error showed up in two of three
/// conditions and in nothing else at all -- every psi_i is Delta*-harmonic, so
/// no convergence rate moved.
///
/// The reason this is a separate case rather than folded into
/// coefficientsSatisfyTheGeometricConstraints is that it compares the two
/// readings side by side and prints both, so a future failure says which
/// convention the coefficients were solved under rather than merely that they
/// are wrong.
BOOST_AUTO_TEST_CASE( nstxUsesTheCorrectAlpha )
{
	SolovievEquilibrium const eq = SolovievEquilibrium::nstx();
	Geometry const &g = nstxAdaptive;

	Curvatures const correct = curvatures( g, false );
	Curvatures const asDelta = curvatures( g, true );

	double const outer = 1.0 + g.eps;
	double const inner = 1.0 - g.eps;

	double const outerCorrect = psiZZ( eq, outer, 0.0 )
	                            + correct.n1*psiR( eq, outer, 0.0 );
	double const outerAsDelta = psiZZ( eq, outer, 0.0 )
	                            + asDelta.n1*psiR( eq, outer, 0.0 );
	double const innerCorrect = psiZZ( eq, inner, 0.0 )
	                            + correct.n2*psiR( eq, inner, 0.0 );
	double const innerAsDelta = psiZZ( eq, inner, 0.0 )
	                           + asDelta.n2*psiR( eq, inner, 0.0 );

	std::printf( "\n  nstx(), Cerfon-Freidberg equatorial curvature conditions:\n"
	             "    outer   alpha = arcsin(delta): %11.4e    alpha = delta: %11.4e\n"
	             "    inner   alpha = arcsin(delta): %11.4e    alpha = delta: %11.4e\n"
	             "    N1 = %.6f (correct) vs %.6f (the misreading)\n",
	             outerCorrect, outerAsDelta, innerCorrect, innerAsDelta,
	             correct.n1, asDelta.n1 );
	std::fflush( stdout );

	BOOST_TEST( std::abs( outerCorrect ) < 1.0e-6,
	            "nstx() does not satisfy the outer equatorial curvature condition "
	            "with alpha = arcsin( delta ): residual " << outerCorrect );
	BOOST_TEST( std::abs( innerCorrect ) < 1.0e-6,
	            "nstx() does not satisfy the inner equatorial curvature condition "
	            "with alpha = arcsin( delta ): residual " << innerCorrect );

	// And the misreading must now be visibly wrong, so that a silent revert to
	// it fails here rather than passing both ways.
	BOOST_TEST( std::abs( outerAsDelta ) > 1.0e-4,
	            "nstx() satisfies the outer condition under BOTH readings of "
	            "alpha, which it cannot -- the test has lost its discrimination" );
}

/// -Delta*( psi ) must equal F at every point, or the whole expansion is
/// mistyped. Checked by central differences, per fixture, as
/// SolovievConvergence.cpp does for the one set it uses.
BOOST_AUTO_TEST_CASE( everyGeometrySourceMatchesTheOperator )
{
	meq::tests::Rectangle const box = standardBox();
	std::vector<SolovievEquilibrium> const cases = {
		SolovievEquilibrium::frcExample1(),
		SolovievEquilibrium::iterExample2(),
		SolovievEquilibrium::nstxExample3()
	};

	for ( SolovievEquilibrium const &eq : cases )
	{
		double worst = 0.0;
		for ( double r = box.rMin; r <= box.rMax + 1.0e-12; r += 0.1 )
		{
			for ( double z = box.zMin; z <= box.zMax + 1.0e-12; z += 0.15 )
			{
				double const deltaStar = eq.deltaStarFD( r, z );
				double const minusF = -eq.f( r, z, 0.0 );
				worst = std::max( worst, std::abs( deltaStar - minusF ) );
			}
		}
		BOOST_TEST( worst < 1.0e-5,
		            "Delta*( psi ) differs from -F by " << worst );
	}
}

/*
 * -------------------------------------------------------------------------
 * The convergence tables
 * -------------------------------------------------------------------------
 */

/// One Newton step, exactly, for every one of these: dF/dpsi is identically
/// zero, so the discrete residual is affine and an exact Jacobian finishes it in
/// one. More than one means the Jacobian is not the derivative of the residual.
namespace
{
	void checkOneStep( SolovievEquilibrium const &eq, char const *label )
	{
		meq::tests::Measurement const point
			= meq::tests::measure( eq, standardBox(), 2, 8 );
		BOOST_TEST( point.newtonIterations <= 1,
		            label << ": Newton took " << point.newtonIterations
		            << " iterations on a source that does not depend on psi" );
	}
}

BOOST_AUTO_TEST_CASE( everyGeometryNeedsOneNewtonStep )
{
	checkOneStep( SolovievEquilibrium::frcExample1(), "Example 1 FRC" );
	checkOneStep( SolovievEquilibrium::iterExample2(), "Example 2 ITER" );
	checkOneStep( SolovievEquilibrium::nstxExample3(), "Example 3 NSTX" );
}

// The ceilings below sit at about three times the measured L2 error on the
// finest mesh, h = 0.025, and the measurement each was set from is recorded
// beside it. They are the only assertion here that can see a solution wrong by a
// constant factor or a sign; the rates cannot. A change to Soloviev.hpp that
// moves them should move them deliberately.
//
// MEASURED RATES, for the record. Every one of the eighteen sequences below
// converges at k+1 in both psi and q, and every solve finishes in one Newton
// step:
//
//   Example 1 FRC        k=1  psi 1.975 1.990 1.996   q 1.960 1.980 1.990
//                        k=2  psi 2.999 3.000 3.000   q 2.984 2.992 2.996
//                        k=3  psi 4.000 4.000 4.000   q 3.984 3.992 3.996
//   Example 2 ITER-like  k=1  psi 1.977 1.991 1.996   q 1.874 1.951 1.979
//                        k=2  psi 2.957 2.989 2.998   q 2.899 2.966 2.987
//                        k=3  psi 3.950 3.988 3.997   q 3.942 3.980 3.992
//   Example 3 NSTX-like  k=1  psi 1.980 1.992 1.996   q 1.978 1.990 1.995
//                        k=2  psi 2.998 3.000 3.000   q 2.952 2.980 2.991
//                        k=3  psi 3.994 3.998 4.000   q 3.950 3.981 3.993
//
// For comparison, refs/HDG-GradShafranov.pdf's own Tables 1, 2 and 3 report, for
// the finest pair and the L2 norm: Example 1 psi 2.13 / 3.14 / 4.44 at
// k = 1 / 2 / 3, Example 2 psi 1.98 / 3.04 / 4.49, Example 3 psi 1.99 / 2.99 /
// 4.23. Those are measured on their own curved domains rather than this
// rectangle, so the error levels are not comparable, but the orders are.

BOOST_AUTO_TEST_CASE( frcConvergesAtDesignOrder )
{
	// Measured at h = 0.025: psi 2.0474e-05, 1.1452e-07, 1.7833e-10 and
	// q 2.9798e-05, 1.3541e-07, 2.1855e-10 for k = 1, 2, 3.
	double const psiCeiling[] = { 0.0, 6.2e-5, 3.5e-7, 5.4e-10 };
	double const fluxCeiling[] = { 0.0, 9.0e-5, 4.1e-7, 6.6e-10 };
	for ( int order = 1; order <= 3; ++order )
	{
		meq::tests::checkOrder( SolovievEquilibrium::frcExample1(),
		                        "Solov'ev Example 1, FRC", order,
		                        psiCeiling[ order ], fluxCeiling[ order ] );
	}
}

BOOST_AUTO_TEST_CASE( iterConvergesAtDesignOrder )
{
	// Measured at h = 0.025: psi 3.0568e-05, 2.6035e-07, 2.0908e-09 and
	// q 1.0634e-04, 1.6946e-06, 1.7670e-08. The flux errors are the largest of
	// the three geometries by a factor of four, which is the price of putting an
	// aspect ratio of 0.32 on a box that reaches past its own separatrix.
	double const psiCeiling[] = { 0.0, 9.2e-5, 7.9e-7, 6.3e-9 };
	double const fluxCeiling[] = { 0.0, 3.2e-4, 5.1e-6, 5.3e-8 };
	for ( int order = 1; order <= 3; ++order )
	{
		meq::tests::checkOrder( SolovievEquilibrium::iterExample2(),
		                        "Solov'ev Example 2, ITER-like", order,
		                        psiCeiling[ order ], fluxCeiling[ order ] );
	}
}

BOOST_AUTO_TEST_CASE( nstxConvergesAtDesignOrder )
{
	// Measured at h = 0.025: psi 2.7726e-05, 1.6382e-07, 4.0421e-10 and
	// q 4.9015e-05, 2.3657e-07, 1.1414e-09.
	double const psiCeiling[] = { 0.0, 8.4e-5, 5.0e-7, 1.3e-9 };
	double const fluxCeiling[] = { 0.0, 1.5e-4, 7.2e-7, 3.5e-9 };
	for ( int order = 1; order <= 3; ++order )
	{
		meq::tests::checkOrder( SolovievEquilibrium::nstxExample3(),
		                        "Solov'ev Example 3, NSTX-like", order,
		                        psiCeiling[ order ], fluxCeiling[ order ] );
	}
}
