#define BOOST_TEST_MODULE HighBetaConvergence
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

#include "mfem.hpp"

#include "meq/GradShafranov.hpp"
#include "meq/Source.hpp"

#include "analytic/HighBetaPoloidal.hpp"
#include "convergence/ConvergenceHarness.hpp"

/*
 * PROFILES IN NORMALISED FLUX, WHICH IS HOW AN EQUILIBRIUM IS ACTUALLY POSED --
 * and psi on the magnetic axis is therefore an unknown of the problem, not an
 * input to it.
 *
 * refs/GourdainContour.pdf section V eq (39) writes p and F^2 as polynomials in
 * Psi = ( psi - psi_bnd )/( psi_ax - psi_bnd ). meq::Profile is tabulated the
 * same way and free boundary has no other option, so this is the shape meq has
 * to live with. What it costs is that psi_ax is a global functional of the
 * solution, which is what CEDRES++ (refs/CEDRES.pdf) means when it warns that a
 * normalised profile "leads to non-local entries in the stiffness matrix".
 *
 * THIS FILE WAS RED FOR A LONG TIME, and what it was red about is worth keeping
 * in view, because two cheaper answers were measured and killed before the
 * expensive one was written.
 *
 *   * psi_ax FIXED. The profile is then inert: at psi_ax = 1 the solution
 *     reaches Psi = 0.0013, so a Psi^(nu-1) pressure gradient is 1e-9 of itself
 *     and amplitudes of 1 and 512 gave answers identical to every digit.
 *     theFixedNormalisationIsADifferentProblem still records that.
 *
 *   * psi_ax ITERATED OUTSIDE THE SOLVER. The outer map psi_ax -> max psi has a
 *     POLE beside its own fixed point -- at nu = 2, amplitude 1 the fixed point
 *     is 0.3059 and the pole is at 0.2996, six parts in a thousand away -- so
 *     the iteration falls off the branch and settles at psi_ax ~ 1e-12, where
 *     psi and psi_ax shrink together and the pressure gradient nu A/psi_ax runs
 *     to 1e12 while the solution it drives sits at 1e-12. Relaxing it harder is
 *     not a fix for a pole.
 *
 * So psi_ax went INSIDE the residual, as a genuine unknown, in
 * GradShafranovSolver::setSource( NormalisedSource &, double ) -- a bordered
 * Newton whose extra row and column ARE the non-local terms. That is what these
 * tests measure.
 */

namespace
{
	using meq::analytic::HighBetaPoloidal;
	using meq::tests::NormalisedEquilibriumSource;
	using meq::tests::NormalisedMeasurement;
	using meq::tests::standardBox;

	double zeroDatum( double, double )
	{
		return 0.0;
	}

	/// The first Dirichlet eigenvalue of the benchmark box, pi^2( 1/w^2 + 1/h^2 ).
	/// CLAUDE.md uses this as the scale dF/dpsi has to be compared against: the
	/// linearised operator is -div_bar( ( 1/r ) grad_bar( . ) ) - ( dF/dpsi )/r,
	/// so a reaction term past lambda_1 has pushed the operator indefinite and
	/// the continuous problem multi-valued.
	double firstEigenvalue()
	{
		meq::tests::Rectangle const box = standardBox();
		double const w = box.rMax - box.rMin;
		double const h = box.zMax - box.zMin;
		return M_PI*M_PI*( 1.0/( w*w ) + 1.0/( h*h ) );
	}

	/*
	 * THE DIMENSIONAL ESTIMATE OF psi_ax, which is what makes "did it find the
	 * physical answer" an assertion rather than a hope.
	 *
	 * Write psi = psi_ax u with max u = 1, which is the self-consistency
	 * condition. The equation -div_bar( ( 1/r ) grad_bar psi ) = F/r with
	 * p = A Psi^nu and psi_bnd = 0 becomes
	 *
	 *     psi_ax^2 L( u ) = r nu A u^(nu-1) + ( g g' term )
	 *
	 * and at the peak L( u ) is about lambda_1 while r is about 1 on this box, so
	 * psi_ax ~ sqrt( nu A / lambda_1 ). It is a scaling and not a solution -- the
	 * measured ratio to it runs from 0.64 at nu = 4 to 1.02 at nu = 2 -- so what
	 * is asserted on is a factor of three either side, which a degenerate fixed
	 * point at 1e-12 misses by ten orders of magnitude.
	 */
	double axisEstimate( int nu, double amplitude )
	{
		return std::sqrt( nu*amplitude/firstEigenvalue() );
	}

	/*
	 * THE STARTING POINT, AND IT IS PART OF THE PROBLEM STATEMENT.
	 *
	 * At a fixed normalisation this equation has a small positive solution and a
	 * large one; measured at nu = 4, amplitude 1, psi_ax = 0.42, they are 3.0e-3
	 * and 6.2e-1. Only the large one can satisfy max psi = psi_ax, so adding the
	 * constraint takes the small branch out of the SOLUTION SET -- but not out of
	 * the iteration's reach, and Newton from the Dirichlet datum walks straight
	 * onto it.
	 *
	 * A separable sine bump of the height psi_ax is expected to have is the
	 * obvious guess and it is a robust one: measured at nu = 4, amplitude 1,
	 * heights of 0.5, 0.8, 1.0 and 1.5 times the estimate all reach
	 * psi_ax = 2.86077085e-01, to every digit printed, in 5 to 9 iterations.
	 * Twice the estimate does not.
	 */
	mfem::FunctionCoefficient bump( double height )
	{
		meq::tests::Rectangle const box = standardBox();
		double const rMin = box.rMin;
		double const zMin = box.zMin;
		double const width = box.width();
		double const depth = box.height();
		return mfem::FunctionCoefficient(
			[ height, rMin, zMin, width, depth ]( mfem::Vector const &x )
			{
				return height*std::sin( M_PI*( x( 0 ) - rMin )/width )
				       *std::sin( M_PI*( x( 1 ) - zMin )/depth );
			} );
	}

	/// max | dF/dpsi | over the box and over the range the SOLUTION occupies,
	/// which is the only range worth sampling and which is why this needs a
	/// converged normalisation to mean anything. Sampling [ 0, psi_ax ] with a
	/// psi_ax the solve never reaches is what once reported ratios of 2523 on
	/// problems that took one Newton step.
	double worstReaction( HighBetaPoloidal const &eq, double psiMin, double psiMax )
	{
		meq::tests::Rectangle const box = standardBox();
		double worst = 0.0;
		for ( int i = 0; i <= 20; ++i )
		{
			double const r = box.rMin + ( box.rMax - box.rMin )*i/20.0;
			for ( int j = 0; j <= 20; ++j )
			{
				double const psi = psiMin + ( psiMax - psiMin )*j/20.0;
				worst = std::max( worst, std::abs( eq.dFdPsi( r, 0.0, psi ) ) );
			}
		}
		return worst;
	}

	/*
	 * THE BEST OBSERVED NEWTON ORDER OVER A TRIPLE THAT IS ENTITLED TO SUPPORT
	 * ONE, which takes two guards and not the one CLAUDE.md's rule names.
	 *
	 * The rule is "assert on the best triple above the round-off floor". Applied
	 * to the histories in this file it reported orders of 5.9, 8.3 and 12.2, and
	 * the 12.2 is manufactured exactly the way that rule warns about: the
	 * amplitude-100 run has a nearly flat pair, 1.12e-1 to 1.07e-1, and a ratio
	 * of 0.955 in the denominator of log( r2/r1 )/log( r1/r0 ) turns an ordinary
	 * following step into a large number. It is a statement about a damped step,
	 * not about the order.
	 *
	 * So: EVERY member of the triple above the floor -- not just the last, or the
	 * round-off floor itself lands in the denominator -- and the first step must
	 * already be contracting, taken as a factor of two. With both, the same six
	 * histories report 5.96, 5.88, 5.87, 5.37, 2.38 and 2.59.
	 *
	 * A value above two is a real super-quadratic step and not an artefact:
	 * these iterations are short, and a damped step followed by an undamped
	 * quadratic one genuinely reduces the residual by more than the square. What
	 * the assertion is for is the OTHER end -- a border that is wrong or missing
	 * gives a constant contraction and an order of one, which is what
	 * theNonLocalTermsAreWhatMakeItConverge measures directly.
	 */
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
			best = std::max( best, meq::tests::newtonOrder( history[ i - 2 ], history[ i - 1 ],
			                                                history[ i ] ) );
		}
		return best;
	}

	/// The whole augmented residual history on one line. Printed rather than
	/// summarised because an order is a claim about a history and the reader has
	/// to be able to check it.
	void printHistory( std::vector<double> const &history )
	{
		std::printf( "        ||( R, gamma G )||:" );
		for ( double value : history )
			std::printf( " %.2e", value );
		std::printf( "\n" );
	}

	/// The six ( nu, amplitude ) pairs the table is run over.
	struct Case
	{
		int nu;
		double amplitude;
	};

	std::vector<Case> const &cases()
	{
		static std::vector<Case> const list = {
			{ 2, 1.0 }, { 2, 10.0 }, { 2, 100.0 },
			{ 4, 1.0 }, { 4, 10.0 }, { 4, 100.0 } };
		return list;
	}
}

/*
 * THE FIXTURE BEFORE THE SOLVER. dF/dpsi is what Newton puts in the Jacobian and
 * an error in it does not change the converged answer -- it only wrecks, or
 * silently slows, the convergence to it. The polynomial chain rule through
 * Psi = psi/psi_axis is exactly where a factor of psi_axis goes missing.
 */
BOOST_AUTO_TEST_CASE( theDerivativeMatchesAFiniteDifference )
{
	HighBetaPoloidal const eq = HighBetaPoloidal::moderate();
	meq::tests::Rectangle const box = standardBox();

	double worst = 0.0;
	for ( int i = 0; i <= 8; ++i )
	{
		double const r = box.rMin + ( box.rMax - box.rMin )*i/8.0;
		for ( int j = 1; j < 20; ++j )
		{
			double const psi = eq.psiAxis()*j/20.0;
			double const step = 1.0e-6*eq.psiAxis();
			double const difference = ( eq.f( r, 0.0, psi + step )
			                            - eq.f( r, 0.0, psi - step ) )/( 2.0*step );
			double const analytic = eq.dFdPsi( r, 0.0, psi );
			double const scale = std::max( 1.0, std::abs( analytic ) );
			worst = std::max( worst, std::abs( difference - analytic )/scale );
		}
	}

	std::printf( "\n  dF/dpsi against a central difference: worst %.3e relative\n", worst );
	std::fflush( stdout );

	BOOST_TEST( worst < 1.0e-8,
	            "dF/dpsi disagrees with a central difference of f by " << worst
	            << " relative, so the chain rule through the normalisation is wrong" );
}

/*
 * AND THERE IS NO TRIVIAL BRANCH HERE, WHICH IS WHY THIS SOURCE IS POSED WITH
 * THE PAPER'S OWN HOMOGENEOUS DATA and the GS-2 sources cannot be.
 *
 * Every one of GS-2 sections 4.2 to 4.5 has F( r, 0 ) = 0, so psi == 0 solves the
 * homogeneous problem and Newton stops on it in zero iterations -- CLAUDE.md's
 * standing trap, and the reason those four are run with an artificial ramp. A
 * profile equilibrium does not have that property: p'( 0 ) = a_1/psi_axis is
 * whatever the first pressure coefficient says, and it is not zero.
 *
 * The guess these tests do supply is a different thing entirely, and the
 * distinction matters. It is not there to keep the iterate off psi == 0, which
 * is not a solution; it is there to keep it off the SMALL branch, which is.
 */
BOOST_AUTO_TEST_CASE( theSourceDoesNotVanishOnTheTrivialBranch )
{
	HighBetaPoloidal const eq = HighBetaPoloidal::moderate();
	meq::tests::Rectangle const box = standardBox();

	double smallest = 1.0e300;
	for ( int i = 0; i <= 8; ++i )
	{
		double const r = box.rMin + ( box.rMax - box.rMin )*i/8.0;
		smallest = std::min( smallest, std::abs( eq.f( r, 0.0, 0.0 ) ) );
	}

	std::printf( "  the smallest | F( r, z, 0 ) | over the box is %.4e\n", smallest );
	std::fflush( stdout );

	BOOST_TEST( smallest > 1.0e-3,
	            "F vanishes at psi = 0 somewhere on the box, so psi == 0 solves "
	            "the homogeneous problem and this source needs a non-trivial "
	            "guess like the GS-2 ones do" );
}

/*
 * THE CONTROL: A FIXED NORMALISATION IS NOT AN APPROXIMATION OF THIS PROBLEM,
 * IT IS A DIFFERENT ONE.
 *
 * Kept because it is the measurement that justifies every line of the bordered
 * Newton. With psi_ax pinned at 1 and a peaked pressure, the solve reaches
 * Psi ~ 1e-3, a Psi^(nu-1) gradient is ~1e-9 of itself, and the pressure profile
 * might as well not be there -- so a five-hundred-fold change in the pressure
 * amplitude changes nothing.
 *
 * What is ASSERTED is the wanted behaviour and not the defect: that the two
 * solves agree, which is what says the profile was never sampled. The next test
 * is where the same amplitudes are required to matter.
 */
BOOST_AUTO_TEST_CASE( theFixedNormalisationIsADifferentProblem )
{
	std::printf( "\n  psi_ax pinned at 1, k = 2, n = 16 -- the control\n" );
	std::printf( "    %6s %10s %14s %10s %10s\n",
	             "nu", "amplitude", "psi_max", "Psi_max", "iterations" );

	std::vector<double> peaks;
	for ( double amplitude : { 1.0, 512.0 } )
	{
		HighBetaPoloidal const eq = HighBetaPoloidal::peaked( 4, amplitude, 1.0 );
		meq::tests::SampleCloud const cloud( standardBox(), 400 );
		meq::tests::SelfMeasurement const point = meq::tests::measureSelf(
			eq, standardBox(), 2, 16, cloud, zeroDatum, 60, 1.0e-10 );

		std::printf( "    %6d %10.1f %14.4e %10.5f %10d\n",
		             4, amplitude, point.psiMax, point.psiMax/eq.psiAxis(),
		             point.newtonIterations );
		std::fflush( stdout );
		peaks.push_back( point.psiMax );
	}

	BOOST_TEST( std::abs( peaks[ 1 ] - peaks[ 0 ] ) < 0.01*peaks[ 0 ],
	            "with psi_ax pinned, raising the pressure amplitude from 1 to 512 "
	            "changed psi_max by "
	            << 100.0*std::abs( peaks[ 1 ] - peaks[ 0 ] )/peaks[ 0 ]
	            << "%. It is supposed to change it by nothing: the point of this "
	            "control is that a pinned normalisation never samples the profile. "
	            "If this now fails the normalisation has leaked into the pinned "
	            "path and theSelfConsistentNormalisation's comparison against it "
	            "no longer isolates anything" );
}

/*
 * THE NORMALISATION SOLVED FOR, WHICH IS THE POINT OF THE FILE.
 *
 * psi_ax is an unknown; the system Newton closes is the trace residual together
 * with psi_ax - max psi_h = 0, and the Jacobian is bordered by the two non-local
 * terms. Four things are asserted on every case, and each of them fails for a
 * different reason:
 *
 *   converged            the bordered Newton reached its tolerance at all
 *   self consistent      psi_ax IS max psi_h, to round-off -- which is the
 *                        equation that was missing, so this is the one that
 *                        distinguishes a solved normalisation from a guessed one
 *   physical             psi_ax is within a factor of three of the dimensional
 *                        estimate. The degenerate fixed point the outer
 *                        iteration used to find misses this by ten orders
 *   superlinear          an observed order above 1.5 on a triple entitled to
 *                        support one, which is the only thing here that can see
 *                        a WRONG border: CLAUDE.md's standing result is that a
 *                        wrong Jacobian leaves every error and every rate
 *                        unchanged and shows up only in the path. Measured, the
 *                        six read 5.96, 5.88, 5.87, 5.37, 2.38 and 2.59 -- these
 *                        iterations are short enough that a damped step followed
 *                        by an undamped quadratic one beats the square
 *
 * and one thing is asserted across them: that the pressure amplitude moves the
 * answer, which the control above shows a pinned normalisation does not.
 *
 * The stiffness ratio is reported here rather than in a test of its own, and now
 * over the range the solution actually occupies. That question -- is a high-beta
 * profile a stiff case for Newton? -- has been unanswerable until now, because
 * with the profile inert there was no non-linearity switched on to be stiff.
 */
BOOST_AUTO_TEST_CASE( theSelfConsistentNormalisation )
{
	double const lambda = firstEigenvalue();
	std::printf( "\n  psi_ax as an unknown, k = 2, n = 8, lambda_1 = %.3f\n", lambda );
	std::printf( "    %4s %10s %13s %13s %11s %9s %9s %7s %6s %s\n",
	             "nu", "amplitude", "psi_ax", "max psi_h", "psi_ax-max", "estimate",
	             "ratio", "|dF|/l1", "Newton", "order" );

	struct Converged
	{
		int nu;
		double amplitude;
		double psiAxis;
	};
	std::vector<Converged> converged;

	for ( Case const &one : cases() )
	{
		double const estimate = axisEstimate( one.nu, one.amplitude );
		mfem::FunctionCoefficient guess = bump( estimate );

		NormalisedMeasurement const point = meq::tests::measureNormalised(
			HighBetaPoloidal::peaked( one.nu, one.amplitude, estimate ),
			standardBox(), 2, 8, zeroDatum, estimate, guess, 40, 1.0e-10 );

		double reaction = 0.0;
		if ( point.converged )
		{
			HighBetaPoloidal converged =
				HighBetaPoloidal::peaked( one.nu, one.amplitude, point.psiAxis );
			reaction = worstReaction( converged, point.psiMin, point.psiMax );
		}

		// The floor is set well above the 1e-12 the augmented residual bottoms
		// out at on these cases, so that no triple used for an order claim has a
		// member sitting on round-off.
		double const order = bestOrder( point.residuals, 1.0e-10 );

		std::printf( "    %4d %10.1f %13.6e %13.6e %11.2e %9.4f %9.4f %7.2f %6d %6.3f%s\n",
		             one.nu, one.amplitude, point.psiAxis, point.psiMax,
		             point.psiAxis - point.psiMax, estimate,
		             point.psiAxis/estimate, reaction/lambda,
		             point.newtonIterations, order,
		             point.converged ? "" : "   FAILED" );
		printHistory( point.residuals );
		std::fflush( stdout );

		BOOST_TEST( point.converged,
		            "nu = " << one.nu << ", amplitude = " << one.amplitude
		            << ": the bordered Newton did not converge" );

		if ( !point.converged )
			continue;

		BOOST_TEST( std::abs( point.psiAxis - point.psiMax ) < 1.0e-8*point.psiAxis,
		            "nu = " << one.nu << ", amplitude = " << one.amplitude
		            << ": psi_ax = " << point.psiAxis << " but max psi_h = "
		            << point.psiMax << ". The normalisation constraint is the extra "
		            "equation the bordered Newton exists to satisfy, so this is what "
		            "says psi_ax was solved for rather than assumed" );

		BOOST_TEST( point.psiAxis > estimate/3.0,
		            "nu = " << one.nu << ", amplitude = " << one.amplitude
		            << ": psi_ax = " << point.psiAxis << " against a dimensional "
		            "estimate of " << estimate << ". An answer this far below the "
		            "scaling is the degenerate branch where psi and psi_ax shrink "
		            "together -- which is what psi_ax OUTSIDE the residual finds" );

		BOOST_TEST( point.psiAxis < 3.0*estimate,
		            "nu = " << one.nu << ", amplitude = " << one.amplitude
		            << ": psi_ax = " << point.psiAxis << " against a dimensional "
		            "estimate of " << estimate );

		BOOST_TEST( order > 1.5,
		            "nu = " << one.nu << ", amplitude = " << one.amplitude
		            << ": the best observed Newton order over the tail is " << order
		            << ", not 2. The bordered Jacobian's two non-local terms are "
		            "the only part of it that a rate table cannot see -- an error "
		            "there converges to the same answer and only costs the order. "
		            "See CLAUDE.md, 'A wrong Jacobian is invisible to a convergence "
		            "table'" );

		converged.push_back( Converged{ one.nu, one.amplitude, point.psiAxis } );
	}

	/*
	 * AND THE AMPLITUDE HAS TO MATTER, which is the whole difference from the
	 * control. Scaling the pressure at fixed F is what raises beta_p, and the
	 * scaling above says psi_ax should go as its square root: a hundredfold
	 * amplitude is a tenfold psi_ax. Asserted loosely, as a factor of three,
	 * because what is being ruled out is a profile that is inert.
	 */
	auto lookup = [ &converged ]( int nu, double amplitude ) -> double
	{
		for ( Converged const &one : converged )
			if ( one.nu == nu && one.amplitude == amplitude )
				return one.psiAxis;
		return 0.0;
	};

	// Paired by ( nu, amplitude ) rather than by position, so that a case that
	// failed above takes its own assertion down and not somebody else's.
	for ( int nu : { 2, 4 } )
	{
		double const weak = lookup( nu, 1.0 );
		double const strong = lookup( nu, 100.0 );
		if ( weak == 0.0 || strong == 0.0 )
			continue;

		BOOST_TEST( strong > 3.0*weak,
		            "nu = " << nu << ": psi_ax went from " << weak << " to "
		            << strong << " when the pressure amplitude was raised a "
		            "hundredfold. It should go as the square root of the "
		            "amplitude, so tenfold. A profile that does not respond to "
		            "its own amplitude has not been sampled" );
	}
}

/*
 * WHAT THE NON-LOCAL TERMS BUY, MEASURED AGAINST THE SAME SOLVER WITHOUT THEM.
 *
 * This is the test ROADMAP.md asked for when it said the work "needs a test that
 * can SEE the missing terms", and the reason it is needed is CLAUDE.md's
 * standing result: a wrong Jacobian leaves every error and every convergence
 * rate untouched and shows up only in the PATH. SourceTests' finite-difference
 * check on dFdPsi structurally cannot see this one, because f() and dFdPsi()
 * would both be evaluated at the same wrong normalisation and would agree with
 * each other.
 *
 * Normalisation::Decoupled is that control. It is the same solver, the same
 * mesh, the same guess, the same line search and the same stopping rule, with
 * exactly three quantities zeroed: dR/dpsi_ax, d( max psi_h )/d lambda and
 * d( max psi_h )/d psi_ax. The step in psi_ax then reduces to
 * psi_ax <- max psi_h, which is precisely psi_ax outside the residual, and the
 * trace step is a Newton step that does not know psi_ax is about to move.
 *
 * Measured at nu = 2, amplitude 1, k = 2, n = 8: the coupled path reaches
 * 4.45e-15 in four iterations and the decoupled one is at 8.24e-2 after fifteen,
 * having started at 8.31e-2 -- and at 8.30e-2 after forty, when it is given
 * them. It does not converge slowly. It does not move. Fifteen is the budget
 * here only because forty of them cost time to prove the same thing.
 */
BOOST_AUTO_TEST_CASE( theNonLocalTermsAreWhatMakeItConverge )
{
	int const order = 2;
	int const n = 8;
	int const nu = 2;
	double const amplitude = 1.0;
	double const estimate = axisEstimate( nu, amplitude );
	int const budget = 15;

	std::printf( "\n  the same solve with and without the border, k = %d, n = %d\n",
	             order, n );

	std::vector<double> reduction;
	for ( int decoupled = 0; decoupled < 2; ++decoupled )
	{
		mfem::Mesh mesh = meq::tests::makeMesh( standardBox(), n );
		NormalisedEquilibriumSource<HighBetaPoloidal> source(
			HighBetaPoloidal::peaked( nu, amplitude, estimate ) );
		mfem::ConstantCoefficient zero( 0.0 );
		mfem::FunctionCoefficient guess = bump( estimate );

		meq::GradShafranovSolver solver( mesh, order );
		solver.setSource( source, estimate );
		solver.setBoundaryData( zero );
		solver.setInitialGuess( guess );
		solver.setNewtonControl( 1.0e-10, 1.0e-14, budget );
		if ( decoupled )
			solver.setNormalisationCoupling(
				meq::GradShafranovSolver::Normalisation::Decoupled );

		bool converged = true;
		try
		{
			solver.solve();
		}
		catch ( std::exception const & )
		{
			converged = false;
		}

		std::vector<double> const &history = solver.newtonResiduals();
		double const first = history.empty() ? 0.0 : history.front();
		double const last = history.empty() ? 0.0 : history.back();

		std::printf( "    %-10s %11s in %2d  psi_ax %13.6e  residual %.3e -> %.3e\n",
		             decoupled ? "decoupled" : "coupled",
		             converged ? "converged" : "NOT converged",
		             solver.newtonIterations(), solver.psiAxis(), first, last );
		printHistory( history );
		std::fflush( stdout );

		reduction.push_back( first > 0.0 && last > 0.0 ? first/last : 1.0 );

		if ( !decoupled )
			BOOST_TEST( converged,
			            "the coupled path did not converge within " << budget
			            << " iterations, so this comparison has no control to be "
			            "a control of" );
	}

	BOOST_TEST( reduction[ 0 ] > 1.0e8*reduction[ 1 ],
	            "the bordered Newton reduced the augmented residual by a factor of "
	            << reduction[ 0 ] << " and the same iteration WITHOUT the two "
	            "non-local terms reduced it by " << reduction[ 1 ]
	            << ". The terms are supposed to be what makes this problem "
	            "solvable at all -- if the gap has closed, either the border has "
	            "stopped being used or the control has stopped being a control, "
	            "and the second is the more likely" );
}

/*
 * THE BORDER'S SPARSITY, MEASURED RATHER THAN ARGUED.
 *
 * The bordered Jacobian's extra ROW is b = -d( max psi_h )/d lambda, and
 * solveWithNormalisation() computes it on the trace dofs of one element only --
 * the element that attained the maximum -- and takes the rest to be exactly
 * zero. That is a claim about the hybridization: static condensation expresses
 * an element's flux and potential in terms of the trace on its own faces and
 * nothing else, so a trace dof elsewhere cannot reach it.
 *
 * If it were false the border would be wrong by whatever was left out, and the
 * symptom would be a Newton order below two rather than a wrong answer -- which
 * is precisely the failure this project has learned a convergence table cannot
 * see. So it is checked directly: perturb a trace dof off the element and the
 * maximum must not move AT ALL, not merely move a little.
 */
BOOST_AUTO_TEST_CASE( theAxisSensitivityIsLocalToItsElement )
{
	int const order = 2;
	int const n = 6;
	double const estimate = axisEstimate( 2, 1.0 );

	mfem::Mesh mesh = meq::tests::makeMesh( standardBox(), n );
	NormalisedEquilibriumSource<HighBetaPoloidal> source(
		HighBetaPoloidal::peaked( 2, 1.0, estimate ) );
	mfem::ConstantCoefficient zero( 0.0 );
	mfem::FunctionCoefficient guess = bump( estimate );

	meq::GradShafranovSolver solver( mesh, order );
	solver.setSource( source, estimate );
	solver.setBoundaryData( zero );
	solver.setInitialGuess( guess );
	solver.setNewtonControl( 1.0e-10, 1.0e-14, 40 );
	solver.solve();

	mfem::Vector trace( solver.trace() );
	int element = -1;
	double const base = solver.axisFlux( trace, &element );

	BOOST_TEST_REQUIRE( element >= 0 );

	// The trace dofs of that element's faces -- in two dimensions, its edges.
	mfem::Array<int> faces, orientations, faceDofs;
	mesh.GetElementEdges( element, faces, orientations );
	std::vector<bool> onElement( trace.Size(), false );
	int onCount = 0;
	for ( int i = 0; i < faces.Size(); ++i )
	{
		solver.traceSpace().GetFaceVDofs( faces[ i ], faceDofs );
		for ( int j = 0; j < faceDofs.Size(); ++j )
		{
			int const d = faceDofs[ j ] >= 0 ? faceDofs[ j ] : -1 - faceDofs[ j ];
			if ( !onElement[ d ] )
				++onCount;
			onElement[ d ] = true;
		}
	}

	double const step = 1.0e-4*std::max( trace.Normlinf(), 1.0 );
	double worstOn = 0.0;
	double worstOff = 0.0;
	int offSampled = 0;

	for ( int d = 0; d < trace.Size(); ++d )
	{
		// Every dof of the element, and a spread sample of the rest: a recovery
		// apiece, and there are far more of the rest.
		bool const sampled = onElement[ d ] || ( d % 23 == 0 && offSampled < 20 );
		if ( !sampled )
			continue;
		if ( !onElement[ d ] )
			++offSampled;

		double const saved = trace( d );
		trace( d ) = saved + step;
		double const moved = std::abs( solver.axisFlux( trace ) - base );
		trace( d ) = saved;

		if ( onElement[ d ] )
			worstOn = std::max( worstOn, moved );
		else
			worstOff = std::max( worstOff, moved );
	}

	std::printf( "\n  the axis sits on element %d of %d, %d trace dofs on its faces\n",
	             element, mesh.GetNE(), onCount );
	std::printf( "    perturbing a trace dof by %.2e moves max psi_h by\n", step );
	std::printf( "      on  the element:  %.4e   (%d dofs)\n", worstOn, onCount );
	std::printf( "      off the element:  %.4e   (%d dofs sampled)\n",
	             worstOff, offSampled );
	std::fflush( stdout );

	BOOST_TEST_REQUIRE( offSampled > 0 );

	BOOST_TEST( worstOff < 1.0e-14*base,
	            "perturbing a trace dof that is NOT on the axis element's faces "
	            "moved max psi_h by " << worstOff << ", against " << base
	            << " for the value itself. The border of the bordered Jacobian "
	            "takes that to be exactly zero, so it is missing a term and the "
	            "Newton order is what will pay for it" );

	BOOST_TEST( worstOn > 1.0e-6*base,
	            "perturbing the axis element's own trace dofs moved max psi_h by "
	            "only " << worstOn << ", against " << base << " for the value "
	            "itself -- so the border is being built where there is nothing to "
	            "find, and this test is not measuring what it claims to" );
}
