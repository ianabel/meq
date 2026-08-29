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

namespace
{
	using meq::tests::EquilibriumSource;
	using meq::tests::SampleCloud;
	using meq::tests::SelfMeasurement;
	using meq::tests::standardBox;

	SampleCloud const &cloud()
	{
		static SampleCloud const c( standardBox() );
		return c;
	}

	double zeroDatum( double, double )
	{
		return 0.0;
	}

	/// The first Dirichlet eigenvalue of the benchmark box, pi^2( 1/w^2 + 1/h^2 ).
	/// CLAUDE.md uses this as the scale dF/dpsi has to be compared against: the
	/// linearised operator is -div_bar( ( 1/r ) grad_bar ) - ( dF/dpsi )/r, so a
	/// reaction term past lambda_1 has pushed the operator indefinite and the
	/// continuous problem multi-valued.
	double firstEigenvalue()
	{
		meq::tests::Rectangle const box = standardBox();
		double const w = box.rMax - box.rMin;
		double const h = box.zMax - box.zMin;
		return M_PI*M_PI*( 1.0/( w*w ) + 1.0/( h*h ) );
	}

	/// max | dF/dpsi | over the box and over a range of psi, which is the half of
	/// the ratio that needs a solution to be honest about. Sampled over
	/// [ 0, psi_axis ] because that is the range a peaked equilibrium occupies.
	double worstReaction( meq::analytic::HighBetaPoloidal const &eq )
	{
		meq::tests::Rectangle const box = standardBox();
		double worst = 0.0;
		for ( int i = 0; i <= 20; ++i )
		{
			double const r = box.rMin + ( box.rMax - box.rMin )*i/20.0;
			for ( int j = 0; j <= 20; ++j )
			{
				double const psi = eq.psiAxis()*j/20.0;
				worst = std::max( worst, std::abs( eq.dFdPsi( r, 0.0, psi ) ) );
			}
		}
		return worst;
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
	meq::analytic::HighBetaPoloidal const eq = meq::analytic::HighBetaPoloidal::moderate();
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
 */
BOOST_AUTO_TEST_CASE( theSourceDoesNotVanishOnTheTrivialBranch )
{
	meq::analytic::HighBetaPoloidal const eq = meq::analytic::HighBetaPoloidal::moderate();
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
 * IS IT ACTUALLY STIFF? MEASURED, NOT ASSUMED.
 *
 * This project has been wrong about stiffness before: GS-2 sections 4.2, 4.3 and
 * 4.5 were called stiff for months, an MFEM work item was requested to fix them,
 * and they were under-resolved. The diagnostic that settled that is the ratio
 * max| dF/dpsi | / lambda_1 -- how far the Jacobian's reaction term has pushed
 * past the first eigenvalue of the operator it is added to. The pedestal sits at
 * 7 and converges; the current hole sits at 26 and is multi-valued.
 *
 * The sweep raises the pressure amplitude at fixed F, which is what raises
 * beta_p, and reports the ratio and the iteration count together.
 */
BOOST_AUTO_TEST_CASE( aFixedNormalisationNeverSamplesTheProfile )
{
	// psi_axis = 1, and the solution does not come close to it. Recorded as a
	// test rather than a comment because it is the argument for doing the
	// normalisation properly, and because a reader looking at the sweep below
	// needs to know why psi_axis is iterated on rather than set.
	std::printf( "\n  a fixed psi_axis = 1, with a peaked pressure\n" );
	std::printf( "    %6s %10s %14s %10s %10s\n",
	             "nu", "amplitude", "psi_max", "Psi_max", "iterations" );

	std::vector<double> peaks;
	for ( double amplitude : { 1.0, 512.0 } )
	{
		meq::analytic::HighBetaPoloidal const eq =
			meq::analytic::HighBetaPoloidal::peaked( 4, amplitude, 1.0 );
		SelfMeasurement const point = meq::tests::measureSelf(
			eq, standardBox(), 2, 16, cloud(), zeroDatum, 60, 1.0e-10 );

		std::printf( "    %6d %10.1f %14.4e %10.5f %10d\n",
		             4, amplitude, point.psiMax, point.psiMax/eq.psiAxis(),
		             point.newtonIterations );
		std::fflush( stdout );
		peaks.push_back( point.psiMax );
	}

	// A five-hundred-fold change in the pressure amplitude that changes NOTHING
	// is the sharpest statement of the problem: at Psi ~ 1e-3 a Psi^3 pressure
	// gradient is 1e-9 of itself, so the solve never sees it.
	// FIVE HUNDRED TIMES THE PRESSURE MUST CHANGE THE EQUILIBRIUM. Ten per cent
	// is a very weak demand on that; measured, it changes it by 0.008%.
	BOOST_TEST( std::abs( peaks[ 1 ] - peaks[ 0 ] ) > 0.1*peaks[ 0 ],
	            "raising the pressure amplitude from 1 to 512 changed psi_max by "
	            << 100.0*std::abs( peaks[ 1 ] - peaks[ 0 ] )/peaks[ 0 ]
	            << "%, so the pressure profile is inert: with psi_axis fixed at 1 "
	            "the solution reaches Psi ~ 1e-3, where a Psi^(nu-1) gradient is "
	            "1e-9 of itself. A normalised profile needs a psi_axis consistent "
	            "with the solution it produces, and meq has no way to impose one -- "
	            "that is the gap this file exists to record. See "
	            "theSelfConsistentNormalisation" );
}

/*
 * THE NORMALISATION CLOSED OUTSIDE THE SOLVER, which is the honest version of
 * "psi_axis fixed": pick one, solve, take the peak, use it as the next psi_axis,
 * repeat. A Picard iteration on the normalisation.
 *
 * That is not what a production solver should do -- psi_axis belongs inside the
 * residual, where the Jacobian can see it -- but it is what makes a normalised
 * profile mean something with the solver meq has today, and it isolates the two
 * questions. Whether the PROFILE SHAPE is hard to solve is answered here, with
 * an ordinary local Jacobian. Whether the NORMALISATION is hard is a separate
 * question about non-local terms, and CLAUDE.md records that the existing
 * finite-difference test cannot even see them.
 */
BOOST_AUTO_TEST_CASE( theSelfConsistentNormalisation )
{
	double const lambda = firstEigenvalue();
	std::printf( "\n  psi_axis iterated to consistency, lambda_1 = %.3f\n", lambda );
	std::printf( "    %4s %9s %4s %11s %11s %9s %8s %s\n",
	             "nu", "amplitude", "it", "psi_axis", "psi_max", "ratio", "Newton",
	             "outcome" );

	for ( int nu : { 2, 4 } )
	{
		for ( double amplitude : { 1.0, 10.0, 100.0 } )
		{
			double psiAxis = 1.0e-3;
			bool consistent = false;
			bool solved = true;
			int newton = 0;
			double peak = 0.0;
			double ratio = 0.0;

			for ( int sweep = 0; sweep < 40 && solved; ++sweep )
			{
				meq::analytic::HighBetaPoloidal const eq =
					meq::analytic::HighBetaPoloidal::peaked( nu, amplitude, psiAxis );

				SelfMeasurement const point = meq::tests::measureSelf(
					eq, standardBox(), 2, 16, cloud(), zeroDatum, 60, 1.0e-10 );
				solved = point.converged;
				newton = point.newtonIterations;
				peak = point.psiMax;

				// The reaction term over the range the SOLUTION occupies, not over
				// [ 0, psi_axis ]: sampling where the solve never goes is what made
				// the first version of this sweep report ratios of 2500 on a
				// problem that took one Newton step.
				ratio = 0.0;
				meq::tests::Rectangle const box = standardBox();
				for ( int i = 0; i <= 20; ++i )
				{
					double const r = box.rMin + ( box.rMax - box.rMin )*i/20.0;
					for ( int j = 0; j <= 20; ++j )
						ratio = std::max( ratio,
						                  std::abs( eq.dFdPsi( r, 0.0, peak*j/20.0 ) ) );
				}
				ratio /= lambda;

				if ( std::abs( peak - psiAxis ) < 1.0e-6*std::abs( peak ) )
				{
					consistent = true;
					break;
				}
				// Relaxed, because the undamped map oscillates: the pressure
				// responds to psi_axis faster than linearly.
				psiAxis += 0.5*( peak - psiAxis );
			}

			std::printf( "    %4d %9.1f %4s %11.4e %11.4e %9.2e %8d %s\n",
			             nu, amplitude, consistent ? "ok" : "--", psiAxis, peak,
			             ratio, newton,
			             !solved ? "SOLVE FAILED"
			                     : ( consistent ? "consistent" : "not consistent" ) );
			std::fflush( stdout );

			/*
			 * THE FIXED POINT HAS TO BE A PHYSICAL ONE, and it is not.
			 *
			 * Dimensionally, psi_axis ~ sqrt( nu A L^2 / 8 ) for this box and this
			 * profile -- of order 0.5, not of order zero. What the outer Picard
			 * actually finds is psi_axis ~ 1e-12, a degenerate fixed point where
			 * psi and psi_axis shrink together, Psi stays O( 1 ), and the pressure
			 * gradient nu A / psi_axis runs away to 1e12 while the solution it is
			 * supposed to drive stays at 1e-12. That is not an equilibrium, it is
			 * the iteration falling off the branch.
			 *
			 * So this asserts the physical answer and FAILS. Relaxing the outer
			 * map harder is not the fix and should not be tried: psi_axis is a
			 * functional of the solution and belongs INSIDE the residual, where
			 * the Jacobian can see the non-local terms it contributes.
			 * CLAUDE.md, under "Newton, and the obligation it creates", records
			 * that those terms are missing and that SourceTests' finite-difference
			 * check cannot detect their absence, because f() and dFdPsi() would be
			 * missing the same ones.
			 */
			BOOST_TEST( solved,
			            "nu = " << nu << ", amplitude = " << amplitude
			            << ": the solve failed while iterating on psi_axis" );
			BOOST_TEST( peak > 1.0e-3,
			            "nu = " << nu << ", amplitude = " << amplitude
			            << ": the self-consistent normalisation settled at psi_axis = "
			            << psiAxis << ", psi_max = " << peak
			            << ". Dimensionally it should be of order "
			            << std::sqrt( nu*amplitude*0.09 )
			            << ". The outer iteration has found the degenerate fixed "
			            "point where psi and psi_axis shrink together, which is "
			            "what happens when the normalisation is left outside the "
			            "residual. See the comment above" );
		}
	}
}

BOOST_AUTO_TEST_CASE( theStiffnessRatioAgainstBeta )
{
	meq::analytic::HighBetaPoloidal const base = meq::analytic::HighBetaPoloidal::moderate();
	double const lambda = firstEigenvalue();

	std::printf( "\n  high beta poloidal: lambda_1 = %.3f on the benchmark box\n", lambda );
	std::printf( "    %14s %14s %10s %8s %14s\n",
	             "profile", "max|dF/dpsi|", "ratio", "n", "outcome" );

	(void)base;
	for ( int nu : { 2, 4, 8 } )
	{
		for ( double amplitude : { 1.0, 8.0, 64.0, 512.0 } )
		{
			meq::analytic::HighBetaPoloidal const eq =
				meq::analytic::HighBetaPoloidal::peaked( nu, amplitude );
			double const reaction = worstReaction( eq );

			for ( int n : { 16, 32 } )
			{
				SelfMeasurement const point = meq::tests::measureSelf(
					eq, standardBox(), 2, n, cloud(), zeroDatum, 60, 1.0e-10 );

				std::printf( "    nu %d  A %7.1f %12.3e %9.2f %6d %11s (%d)  psi in [%.3e, %.3e]  Psi_max %.4f\n",
				             nu, amplitude, reaction, reaction/lambda, n,
				             point.converged ? "converged" : "FAILED",
				             point.newtonIterations, point.psiMin, point.psiMax,
				             point.psiMax/eq.psiAxis() );
				std::fflush( stdout );
			}
		}
	}
}
