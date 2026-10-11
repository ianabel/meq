#define BOOST_TEST_MODULE PrescribedAxisPoint
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <memory>
#include <vector>

#include "mfem.hpp"

#include "meq/GradShafranov.hpp"
#include "meq/Profiles.hpp"
#include "meq/Source.hpp"

#include "analytic/HighBetaPoloidal.hpp"
#include "convergence/ConvergenceHarness.hpp"

/*
 * MIRROR-PLAN.md MR-2: AxisConstraint::PrescribedPoint, MEASURED ON A TOKAMAK.
 *
 * A mirror has no magnetic axis -- psi rises monotonically outward from R = 0 --
 * so neither NodalMaximum nor LocatedAxis has anything to constrain psi_ax to.
 * PrescribedPoint pins it at a point the caller names, with setBoundaryFluxPoint()'s
 * row reused for the other label: the containing element's potential shape
 * functions at the point, located once per solve, exact and undifferenced.
 *
 * **IT NEEDS NO MIRROR TO BE TESTED, WHICH IS WHAT MAKES IT A SAFE STAGE.** On a
 * tokamak, pin psi_ax at the point LocatedAxis converged to and the two
 * constraints define the SAME discrete equation at that solution: LocatedAxis
 * says psi_ax = psi_h( x* ) with x* a zero of q_h, PrescribedPoint says
 * psi_ax = psi_h( x_p ), and x_p = x* there. So LocatedAxis's answer solves the
 * prescribed system, and a Newton started from the same guess must reach it --
 * to the root finder's accuracy in x*, not to a discretisation error.
 *
 * WHAT THAT AGREEMENT IS WORTH DEPENDS ON A CONTROL, and the control is the same
 * pin moved a quarter of an element outboard. A constraint that ignored its point, or a row that
 * ignored its constraint, would agree with LocatedAxis on BOTH -- so the
 * displaced pin must DISAGREE, by orders of magnitude more than the agreement
 * tolerance, or the agreement says nothing.
 *
 * TWO FIXTURES, BOTH BORDERED NEWTONS CLOSED THROUGH psi_ax:
 *
 *   * HighBetaConvergence's nu = 4 high-beta source, psi_bnd = 0 -- the stiff
 *     case, |dF/dpsi|/lambda_1 about 13, where the nodal maximum sits nine per
 *     cent above the field's only O-point. Pinning at the wrong one of those two
 *     is a different equilibrium, which is what makes it a fixture with teeth.
 *   * FB-3's two-border case -- constant profiles, so F depends on psi only
 *     through the span, and psi_bnd an unknown pinned at its own point. Both
 *     labels are then prescribed-point rows, which is the configuration a mirror
 *     will run in.
 *
 * THE REFUSALS ARE IN SolverContract.cpp, beside every other refusal the class
 * makes, as aPrescribedAxisPointIsRefusedRatherThanDefaulted: they need no
 * solve.
 */

namespace
{
	using meq::analytic::HighBetaPoloidal;
	using meq::tests::NormalisedEquilibriumSource;
	using meq::tests::standardBox;

	int const order = 2;

	/// HighBetaConvergence's guess, for its reason: a separable sine bump of the
	/// expected height keeps Newton off the small-psi branch.
	mfem::FunctionCoefficient bump( double height )
	{
		meq::tests::Rectangle const box = standardBox();
		double const minRadius = box.minRadius;
		double const zMin = box.zMin;
		double const width = box.width();
		double const depth = box.height();
		return mfem::FunctionCoefficient(
			[ height, minRadius, zMin, width, depth ]( mfem::Vector const &x )
			{
				return height*std::sin( M_PI*( x( 0 ) - minRadius )/width )
				       *std::sin( M_PI*( x( 1 ) - zMin )/depth );
			} );
	}

	/// sqrt( nu A/lambda_1 ), HighBetaConvergence's dimensional estimate of psi_ax.
	double axisEstimate( int nu, double amplitude )
	{
		meq::tests::Rectangle const box = standardBox();
		double const w = box.width();
		double const h = box.height();
		double const lambda = M_PI*M_PI*( 1.0/( w*w ) + 1.0/( h*h ) );
		return std::sqrt( nu*amplitude/lambda );
	}

	/// HighBetaConvergence's bestOrder(), for its reasons: every member of the
	/// triple above the floor, and the first step already contracting by two.
	double bestOrder( std::vector<double> const &history, double floorValue )
	{
		double best = 0.0;
		for ( std::size_t i = 2; i < history.size(); ++i )
		{
			if ( history[ i ] < floorValue || history[ i - 1 ] < floorValue
			     || history[ i - 2 ] < floorValue )
				continue;
			if ( !( history[ i ] < history[ i - 1 ] && history[ i - 1 ] < history[ i - 2 ] ) )
				continue;
			if ( history[ i - 1 ] > 0.5*history[ i - 2 ] )
				continue;
			best = std::max( best, meq::tests::newtonOrder( history[ i - 2 ],
			                                                history[ i - 1 ],
			                                                history[ i ] ) );
		}
		return best;
	}

	/// What one solve leaves behind, read off the solver before it goes.
	struct Solved
	{
		bool converged = false;
		double psiAxis = 0.0;
		double psiBoundary = 0.0;
		double axisR = 0.0;
		double axisZ = 0.0;
		bool located = false;
		double constraint = 0.0;
		int iterations = 0;
		std::vector<double> residuals;
		/// psi_h's coefficients, so that the field can be differenced against
		/// another solve on an identically built mesh.
		mfem::Vector psi;
		double psiNorm = 0.0;
	};

	/// The axis pin, or none: `set` false is LocatedAxis, the default.
	struct Pin
	{
		bool set = false;
		double radius = 0.0;
		double z = 0.0;
	};

	void readBack( meq::GradShafranovSolver const &solver, Solved &out )
	{
		out.psiAxis = solver.psiAxis();
		out.psiBoundary = solver.psiBoundary();
		out.axisR = solver.axisR();
		out.axisZ = solver.axisZ();
		out.located = solver.axisWasLocated();
		out.constraint = solver.normalisationResidual();
		out.iterations = solver.newtonIterations();
		out.residuals = solver.newtonResiduals();
		out.psi = solver.potential();
		mfem::ConstantCoefficient zero( 0.0 );
		out.psiNorm = solver.potential().ComputeL2Error( zero );
	}

	/// The L2 norm of the difference of two psi_h on identically built meshes,
	/// relative to the first. The space is the second solve's, which is the
	/// same space by construction.
	double fieldGap( meq::GradShafranovSolver const &solver, Solved const &reference )
	{
		mfem::GridFunction difference( solver.potential() );
		difference -= reference.psi;
		mfem::ConstantCoefficient zero( 0.0 );
		return difference.ComputeL2Error( zero )
		       /std::max( reference.psiNorm, 1.0e-300 );
	}

	/// The stiff high-beta fixture, psi_bnd = 0.
	Solved solveHighBeta( int n, Pin const &pin, double *gap = nullptr,
	                      Solved const *reference = nullptr )
	{
		int const nu = 4;
		double const amplitude = 1.0;
		double const estimate = axisEstimate( nu, amplitude );

		mfem::Mesh mesh = meq::tests::makeMesh( standardBox(), n );
		NormalisedEquilibriumSource<HighBetaPoloidal> source(
			HighBetaPoloidal::peaked( nu, amplitude, estimate ) );
		mfem::ConstantCoefficient zero( 0.0 );
		mfem::FunctionCoefficient guess = bump( estimate );

		meq::GradShafranovSolver solver( mesh, order );
		if ( pin.set )
			solver.setAxisFluxPoint( pin.radius, pin.z );
		solver.setSource( source, estimate );
		solver.setBoundaryData( zero );
		solver.setInitialGuess( guess );
		solver.setNewtonControl( 1.0e-12, 1.0e-14, 40 );

		Solved out;
		out.converged = true;
		try
		{
			solver.solve();
		}
		catch ( std::exception const & )
		{
			out.converged = false;
		}
		readBack( solver, out );
		if ( gap && reference )
			*gap = fieldGap( solver, *reference );
		return out;
	}

	/// FB-3's two-border fixture: psi_bnd pinned at a limiter point too.
	Solved solveTwoBorder( int n, Pin const &pin, double *gap = nullptr,
	                       Solved const *reference = nullptr )
	{
		meq::tests::Rectangle const box = standardBox();
		double const limiterR = box.minRadius + 0.68*box.width();
		double const limiterZ = box.zMin + 0.31*box.height();

		mfem::Mesh mesh = meq::tests::makeMesh( box, n );
		auto pPrime = std::make_shared<meq::ConstantProfile const>( 0.45 );
		auto ggPrime = std::make_shared<meq::ConstantProfile const>( 0.30 );
		meq::NormalisedMHDSource source( pPrime, ggPrime, 1.0, 1.0 );
		mfem::ConstantCoefficient zero( 0.0 );
		mfem::FunctionCoefficient guess = bump( 0.30 );

		meq::GradShafranovSolver solver( mesh, order );
		solver.setBoundaryFluxPoint( limiterR, limiterZ );
		if ( pin.set )
			solver.setAxisFluxPoint( pin.radius, pin.z );
		solver.setSource( source, 0.30 );
		solver.setBoundaryData( zero );
		solver.setInitialGuess( guess );
		solver.setNewtonControl( 1.0e-12, 1.0e-14, 40 );

		Solved out;
		out.converged = true;
		try
		{
			solver.solve();
		}
		catch ( std::exception const & )
		{
			out.converged = false;
		}
		readBack( solver, out );
		if ( gap && reference )
			*gap = fieldGap( solver, *reference );
		return out;
	}

	void printRow( char const *label, int n, Solved const &one, double gap )
	{
		std::printf( "    %-10s %4d %17.10e %17.10e %9.4f %9.4f %10.2e %10.2e %6d%s\n",
		             label, n, one.psiAxis, one.psiBoundary, one.axisR, one.axisZ,
		             gap, std::abs( one.constraint ), one.iterations,
		             one.converged ? "" : "   NOT CONVERGED" );
	}

	void printHistory( std::vector<double> const &history )
	{
		std::printf( "        ||( R, gamma G )||:" );
		for ( double value : history )
			std::printf( " %.2e", value );
		std::printf( "\n" );
	}

	/*
	 * ONE FIXTURE, THREE SOLVES: LocatedAxis, the pin at its axis, and the pin
	 * a quarter of an element outboard. Returns nothing; asserts as it goes.
	 */
	template<typename SolveFunction>
	void reproduceLocatedAxis( char const *name, int n, SolveFunction solveOne )
	{
		meq::tests::Rectangle const box = standardBox();

		std::printf( "\n  %s, k = %d, n = %d\n", name, order, n );
		std::printf( "    %-10s %4s %17s %17s %9s %9s %10s %10s %6s\n",
		             "axis", "n", "psi_ax", "psi_bnd", "R", "z", "L2 gap",
		             "|G|", "newton" );

		Solved const located = solveOne( n, Pin{}, nullptr, nullptr );
		printRow( "located", n, located, 0.0 );
		printHistory( located.residuals );

		BOOST_TEST_REQUIRE( located.converged,
		                    name << ", n = " << n << ": the LocatedAxis reference "
		                    "did not converge, so there is nothing to reproduce" );
		BOOST_TEST_REQUIRE( located.located,
		                    name << ", n = " << n << ": the reference converged "
		                    "without LOCATING an axis -- it fell back to the nodal "
		                    "maximum -- so its axis point is not a zero of q_h and "
		                    "pinning there reproduces nothing in particular" );

		double gap = 0.0;
		Pin const atAxis{ true, located.axisR, located.axisZ };
		Solved const pinned = solveOne( n, atAxis, &gap, &located );
		printRow( "prescribed", n, pinned, gap );
		printHistory( pinned.residuals );

		/*
		 * THE CONTROL: THE SAME PIN, A QUARTER OF AN ELEMENT OUTBOARD. The
		 * O-point is a maximum, so psi there is lower by about half the
		 * curvature times the offset squared, and the profiles are normalised
		 * by it -- a different equilibrium, by far more than any agreement
		 * tolerance below.
		 *
		 * A QUARTER AND NOT A WHOLE ELEMENT, MEASURED. A whole element at
		 * n = 8 on the stiff nu = 4 source is 0.1 m off the axis, and that
		 * pin does not converge in 40 steps; a control that fails to solve
		 * disagrees trivially and so tests nothing. It must CONVERGE to a
		 * different answer, which is required below.
		 */
		double controlGap = 0.0;
		Pin const displaced{ true, located.axisR + 0.25*box.width()/n,
		                     located.axisZ };
		Solved const control = solveOne( n, displaced, &controlGap, &located );
		printRow( "displaced", n, control, controlGap );
		std::fflush( stdout );

		BOOST_TEST_REQUIRE( pinned.converged,
		                    name << ", n = " << n << ": the prescribed-point solve "
		                    "did not converge from the guess LocatedAxis converged "
		                    "from" );

		double const scale = std::abs( located.psiAxis );
		double const axisGap = std::abs( pinned.psiAxis - located.psiAxis )/scale;
		double const boundaryGap =
			std::abs( pinned.psiBoundary - located.psiBoundary )/scale;

		/*
		 * TO THE ROOT FINDER'S ACCURACY, NOT TO A DISCRETISATION ERROR.
		 *
		 * x_p is LocatedAxis's x*, which a Newton on q_h = 0 put there to its
		 * own tolerance; psi_h( x_p ) - psi_h( x* ) is grad psi_h . dx, and
		 * grad psi_h at a zero of q_h is O( h^k ) rather than zero (M-120). So
		 * the two equilibria differ by that and by the two Newtons' own
		 * stopping points, and nothing else. 1e-9 is HighBetaConvergence's own
		 * bound for a LocatedAxis constraint residual, for the same reason.
		 */
		BOOST_TEST( axisGap < 1.0e-9,
		            name << ", n = " << n << ": psi_ax came out "
		            << pinned.psiAxis << " pinned at the located axis against "
		            << located.psiAxis << " located -- " << axisGap
		            << " relative. At the located point the two constraints are "
		            "the same equation, so this is a wrong row, a wrong point or "
		            "psi_c added on one path and not the other" );
		BOOST_TEST( boundaryGap < 1.0e-9,
		            name << ", n = " << n << ": psi_bnd came out "
		            << pinned.psiBoundary << " against " << located.psiBoundary );
		BOOST_TEST( gap < 1.0e-9,
		            name << ", n = " << n << ": psi_h differs from LocatedAxis's "
		            "by " << gap << " in relative L2" );

		// The constraint closes as an exact row should: no root search stands
		// between it and machine zero, so it is held to Newton's own floor.
		BOOST_TEST( std::abs( pinned.constraint ) < 1.0e-12*scale,
		            name << ", n = " << n << ": the prescribed constraint reads "
		            << pinned.constraint << ". The row is exact, so it closes to "
		            "round-off or it is not the row" );

		// What it reports: the point, and NOT a located axis.
		BOOST_TEST( !pinned.located,
		            name << ": axisWasLocated() is true under PrescribedPoint, "
		            "which looked for nothing" );
		BOOST_TEST( pinned.axisR == located.axisR );
		BOOST_TEST( pinned.axisZ == located.axisZ );

		/*
		 * AND THE ORDER, which is what the row's Jacobian entries buy and what no
		 * agreement of answers can see -- a row that was dropped or wrong
		 * converges to the same equilibrium, linearly. Floor well above where
		 * the history bottoms out, as HighBetaConvergence sets it.
		 */
		double const observed = bestOrder( pinned.residuals, 1.0e-10 );
		std::printf( "    best observed Newton order, prescribed: %.3f\n", observed );
		BOOST_TEST( observed > 1.5,
		            name << ", n = " << n << ": the best observed Newton order "
		            "under PrescribedPoint is " << observed << ". The row is the "
		            "shape functions at a fixed point and the corner is 1, both "
		            "exact, so anything short of quadratic is a Jacobian defect" );

		/*
		 * AND WHERE THE TAIL LANDS, WHICH IS THE SHARPER OF THE TWO.
		 *
		 * The stopping threshold is max( 1e-12 r_0, 1e-14 ), about 1e-13 here.
		 * A quadratic step from just above it lands at round-off, a few times
		 * 1e-16 on these meshes; a LINEAR tail crosses it by its contraction
		 * factor and stops just under it. Measured with the row scaled by 0.99
		 * -- a one per cent Jacobian error -- all four runs still converge,
		 * the best-triple order above still reads 1.54 to 2.16 and passes, and
		 * the histories finish at 7.9e-14 to 4.3e-13 on a contraction of 0.02
		 * to 0.04 a step. The exact row finishes at 2e-16 to 5e-16.
		 *
		 * LocatedAxis's own row FAILS this, and that is a measurement rather
		 * than a defect of this test: its tail on the n = 8 high-beta case is a
		 * clean 0.073 a step, because AxisRow::PositionDropped omits the
		 * envelope term M-120 measures. A prescribed point has no position term
		 * to omit.
		 */
		double const last = pinned.residuals.empty()
		                    ? std::numeric_limits<double>::infinity()
		                    : pinned.residuals.back();
		BOOST_TEST( last < 1.0e-14,
		            name << ", n = " << n << ": the prescribed-point Newton "
		            "stopped at " << last << ", just under its threshold rather "
		            "than at round-off. That is a linear tail, and with a row and "
		            "a corner that are both exact it means one of them is not" );

		// THE CONTROL HAS TO DISAGREE, or the agreement above is not evidence.
		BOOST_TEST( control.converged,
		            name << ", n = " << n << ": the displaced control did not "
		            "converge, so its disagreement below would be a failure to "
		            "solve rather than a different equilibrium" );
		double const controlAxisGap =
			std::abs( control.psiAxis - located.psiAxis )/scale;
		BOOST_TEST( controlAxisGap > 1.0e-6,
		            name << ", n = " << n << ": pinning psi_ax a quarter element "
		            "from the axis moved it by only " << controlAxisGap
		            << " relative. Then the agreement above would hold whatever "
		            "the point, and it is not evidence that the point was used" );
	}
}

BOOST_AUTO_TEST_CASE( aPinAtTheLocatedAxisReproducesTheStiffHighBetaEquilibrium )
{
	reproduceLocatedAxis( "high beta, nu = 4, psi_bnd = 0", 8, solveHighBeta );
	reproduceLocatedAxis( "high beta, nu = 4, psi_bnd = 0", 16, solveHighBeta );
}

BOOST_AUTO_TEST_CASE( aPinAtTheLocatedAxisReproducesTheTwoBorderEquilibrium )
{
	reproduceLocatedAxis( "FB-3 two borders", 8, solveTwoBorder );
	reproduceLocatedAxis( "FB-3 two borders", 16, solveTwoBorder );
}
