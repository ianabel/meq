#define BOOST_TEST_MODULE RotatingNormalisedConvergence

#include <boost/test/unit_test.hpp>

#include <cmath>
#include <cstdio>
#include <limits>
#include <memory>
#include <vector>

#include "mfem.hpp"

#include "meq/GradShafranov.hpp"
#include "meq/RotatingSource.hpp"

#include "convergence/ConvergenceHarness.hpp"

/*
 * FL-7 OF FLOW-PLAN.md: A ROTATING SOURCE IN NORMALISED FLUX, THROUGH THE
 * BORDERED NEWTON.
 *
 * Profiles specified against Psi = psi/psi_ax make psi_ax a functional of the
 * solution, so it is an unknown of the non-linear system rather than an input,
 * and GradShafranovSolver::setSource( NormalisedSource &, double ) closes the
 * pair by a bordered Newton. CLAUDE.md's "Newton, and the obligation it creates"
 * is the account of that machinery; HighBetaConvergence is its acceptance
 * criterion. This file asks a narrower question: does a ROTATING source work
 * through it?
 *
 * WHY THAT IS NOT OBVIOUS FROM EITHER PIECE ALONE. meq::NormalisedRotatingSource
 * is a wrapper, and RotatingSourceTests checks its scaling pointwise -- one
 * factor of 1/psi_ax in f() and two in dFdPsi(). But the border differences the
 * ASSEMBLED residual with respect to psi_ax, so it sees the source through a
 * whole hybridized solve, and a source whose setNormalisation() allocated, or
 * cached, or answered for the previous value would pass every pointwise test and
 * corrupt the border. The constraint residual below is what sees that.
 *
 * THE PROFILES ARE CHOSEN SO THAT THERE IS NO TRIVIAL BRANCH. Every GS-2
 * section 4.2-4.5 source vanishes at psi = 0, so psi == 0 solves the homogeneous
 * problem and Newton lands on it in zero iterations; see CLAUDE.md, Traps. Here
 * n_s0( Psi ) has a non-zero slope at Psi = 0, so dp/dPsi and hence F are
 * non-zero on the trivial branch and the iteration cannot rest there. That is a
 * property of the profiles and is asserted below rather than assumed.
 */

namespace
{

	using meq::tests::standardBox;

	// Normalised units, mu0 = 1, as the rest of the rotating tests use.
	double const referenceRadius = 1.0;
	double const speciesTemperature = 0.5;     // T_i = T_e, so T_i + T_e = 1
	double const densityScale = 0.1;
	double const densitySlope = 1.0;           // a: sets F on the trivial branch
	double const densityCurvature = 0.2;       // b: sets dF/dpsi, and must stay
	                                           //    well under lambda_1 = 22.3
	double const ggPrimeConstant = 0.05;
	double const rotationRate = 1.2;

	class PolynomialProfile : public meq::Profile
	{
		public:
			PolynomialProfile( double c0, double c1, double c2 ) : a( c0 ), b( c1 ), c( c2 ) {}

			double operator()( double psi ) const override { return a + b*psi + c*psi*psi; }
			double prime( double psi ) const override { return b + 2.0*c*psi; }
			double doublePrime( double ) const override { return 2.0*c; }

		private:
			double a, b, c;
	};

	std::shared_ptr<meq::Profile const> constantProfile( double v )
	{
		return std::make_shared<meq::ConstantProfile const>( v );
	}

	std::vector<meq::Species> rotatingSpecies()
	{
		std::vector<meq::Species> species( 2 );

		species[ 0 ].mass = 1.0;
		species[ 0 ].charge = 1.0;
		species[ 0 ].temperature = constantProfile( speciesTemperature );
		species[ 0 ].density = std::make_shared<PolynomialProfile const>(
			densityScale, densityScale*densitySlope, densityScale*densityCurvature );

		species[ 1 ].mass = 1.0e-4;
		species[ 1 ].charge = -1.0;
		species[ 1 ].temperature = constantProfile( speciesTemperature );
		species[ 1 ].density = meq::neutralisingDensity( species, 1 );

		return species;
	}

	std::shared_ptr<meq::NormalisedRotatingSource> makeSource( double psiAxis, double omega )
	{
		return std::make_shared<meq::NormalisedRotatingSource>( rotatingSpecies(),
			omega == 0.0 ? nullptr : constantProfile( omega ),
			constantProfile( ggPrimeConstant ),
			referenceRadius, psiAxis, 1.0 );
	}

	struct Solved
	{
		bool converged;
		double psiAxis;
		double psiMax;
		double constraint;
		int newtonIterations;
		std::vector<double> residuals;
		int traceDofs;
	};

	Solved solve( double omega, int order, int n, double psiAxisGuess = 0.1 )
	{
		mfem::Mesh mesh = meq::tests::makeMesh( standardBox(), n );
		std::shared_ptr<meq::NormalisedRotatingSource> source = makeSource( psiAxisGuess, omega );

		mfem::ConstantCoefficient zero( 0.0 );

		meq::GradShafranovSolver solver( mesh, order );
		solver.setSource( *source, psiAxisGuess );
		solver.setBoundaryData( zero );
		solver.setNewtonControl( 1.0e-11, 1.0e-14, 40 );

		Solved out{};
		out.converged = true;
		try
		{
			solver.solve();
		}
		catch ( std::exception const & )
		{
			out.converged = false;
		}

		out.psiAxis = solver.psiAxis();
		out.constraint = solver.normalisationResidual();
		out.newtonIterations = solver.newtonIterations();
		out.residuals = solver.newtonResiduals();
		out.traceDofs = solver.numTraceDofs();
		out.psiMax = out.converged ? solver.potential().Max() : 0.0;
		return out;
	}

}

/// The trivial branch first, because if F vanished at psi = 0 every number below
/// would be a solve that never started. See CLAUDE.md, Traps.
BOOST_AUTO_TEST_CASE( theSourceDoesNotVanishOnTheTrivialBranch )
{
	std::shared_ptr<meq::NormalisedRotatingSource> const source = makeSource( 0.1, rotationRate );

	double smallest = std::numeric_limits<double>::infinity();
	meq::tests::Rectangle const box = standardBox();
	for ( double r = box.rMin; r <= box.rMax + 1.0e-12; r += 0.05 )
		smallest = std::min( smallest, std::fabs( source->f( r, 0.0, 0.0 ) ) );

	BOOST_TEST( smallest > 0.1,
	            "the smallest |F( r, z, 0 )| over the box is " << smallest
	            << ", so psi = 0 nearly solves the homogeneous problem and Newton will "
	            "land on the trivial branch rather than on an equilibrium" );
}

/// The property the border exists to enforce. psi_ax is not an input here: it is
/// an unknown, and the constraint psi_ax = max psi_h is the extra equation the
/// bordered row supplies. If the border were doing nothing this would come back
/// at the initial guess instead.
BOOST_AUTO_TEST_CASE( theNormalisationIsSelfConsistent )
{
	std::printf( "\n  Rotating source in normalised flux, k = 2\n" );
	std::printf( "  %8s %9s %14s %14s %12s %7s\n",
	             "h", "trace", "psi_ax", "max psi_h", "constraint", "Newton" );

	for ( int n : { 8, 16, 24 } )
	{
		Solved const s = solve( rotationRate, 2, n );

		std::printf( "  %8.5f %9d %14.6e %14.6e %12.3e %7d\n",
		             standardBox().width()/n, s.traceDofs, s.psiAxis, s.psiMax,
		             s.psiAxis - s.psiMax, s.newtonIterations );
		std::fflush( stdout );

		BOOST_TEST( s.converged, "n = " << n << ": the bordered Newton did not converge" );
		BOOST_TEST( std::fabs( s.psiAxis - s.psiMax ) < 1.0e-8*std::fabs( s.psiAxis ),
		            "n = " << n << ": psi_ax = " << s.psiAxis << " but max psi_h = " << s.psiMax
		            << ", so the normalisation constraint is not satisfied and the border "
		            "is not closing the system" );
		BOOST_TEST( s.psiAxis > 0.0,
		            "n = " << n << ": psi_ax = " << s.psiAxis << ", but F > 0 with a zero "
		            "datum must give a positive psi by the maximum principle" );

		// And it moved off the guess, or the constraint is satisfied trivially.
		BOOST_TEST( std::fabs( s.psiAxis - 0.1 ) > 1.0e-3*0.1,
		            "n = " << n << ": psi_ax came back at the initial guess, so nothing solved it" );
	}
}

/// Rotation must reach the answer, or this file is HighBetaConvergence with more
/// steps. The comparison is against the same profiles at rest.
BOOST_AUTO_TEST_CASE( rotationMovesTheEquilibrium )
{
	Solved const still = solve( 0.0, 2, 16 );
	Solved const spinning = solve( rotationRate, 2, 16 );

	BOOST_TEST( still.converged, "the solve at rest failed" );
	BOOST_TEST( spinning.converged, "the rotating solve failed" );

	double const relative = std::fabs( spinning.psiAxis - still.psiAxis )/std::fabs( still.psiAxis );
	std::printf( "\n  psi_ax at rest %.6e, rotating %.6e, relative difference %.3e\n",
	             still.psiAxis, spinning.psiAxis, relative );
	std::fflush( stdout );

	BOOST_TEST( relative > 1.0e-2,
	            "psi_ax moves by only " << relative << " between omega = 0 and omega = " << rotationRate
	            << ", so the centrifugal term is not reaching the equilibrium" );
}

/// Newton's order, which is the only thing that can see a wrong Jacobian: a
/// wrong dF/dpsi leaves the converged answer and every rate untouched and drops
/// the observed order to one.
///
/// THE ASSERTION IS ON THE TAIL, NOT ON THE BEST TRIPLE, AND THAT MATTERS HERE.
/// This history opens 6.24e-02 -> 4.01e-02 -> 7.44e-03, which is the iterate
/// walking into the basin rather than converging in it, and the "order" of that
/// first triple reads 3.81. Taking the best over all triples would therefore
/// report 3.81 and pass -- for entirely the wrong reason, and it would go on
/// passing with a Jacobian degraded enough to destroy the tail. CLAUDE.md's
/// testing stance says only a monotone tail supports an order claim; this takes
/// the last triple whose three residuals are all above the round-off floor, and
/// bounds it on BOTH sides, because 1 is a broken Jacobian and 3.8 is an
/// artefact.
BOOST_AUTO_TEST_CASE( theBorderedNewtonKeepsQuadraticOrderInItsTail )
{
	Solved const s = solve( rotationRate, 2, 16 );
	BOOST_TEST( s.converged, "the solve did not converge, so there is no order to measure" );
	BOOST_TEST( s.residuals.size() >= 4u, "too short a history to measure an order from" );

	// The floor: below this the residual is round-off on the initial one, and a
	// ratio of two round-offs is not a convergence rate.
	double const floor = 1.0e-13*s.residuals.front();

	std::printf( "\n  bordered Newton residual history, rotating, k = 2, n = 16\n" );
	std::printf( "    %2s %14s %8s\n", "it", "residual", "order" );
	for ( std::size_t i = 0; i < s.residuals.size(); ++i )
	{
		if ( i >= 2 )
			std::printf( "    %2zu %14.6e %8.3f\n", i, s.residuals[ i ],
			             meq::tests::newtonOrder( s.residuals[ i - 2 ], s.residuals[ i - 1 ], s.residuals[ i ] ) );
		else
			std::printf( "    %2zu %14.6e %8s\n", i, s.residuals[ i ], "-" );
	}
	std::fflush( stdout );

	bool found = false;
	double tailOrder = 0.0;
	std::size_t tailAt = 0;
	for ( std::size_t i = 2; i < s.residuals.size(); ++i )
	{
		if ( s.residuals[ i ] < floor || s.residuals[ i - 1 ] < floor || s.residuals[ i - 2 ] < floor )
			continue;

		tailOrder = meq::tests::newtonOrder( s.residuals[ i - 2 ], s.residuals[ i - 1 ], s.residuals[ i ] );
		tailAt = i;
		found = true;
	}

	BOOST_TEST( found, "no triple of residuals sits above the round-off floor " << floor );
	std::printf( "  tail order %.3f, at iteration %zu, floor %.3e\n", tailOrder, tailAt, floor );
	std::fflush( stdout );

	BOOST_TEST( tailOrder > 1.7,
	            "the observed Newton order in the tail is " << tailOrder << ", not the 2 a correct "
	            "Jacobian gives. A wrong dF/dpsi presents exactly like this, with the converged "
	            "answer and every convergence rate unchanged" );
	BOOST_TEST( tailOrder < 2.5,
	            "the observed Newton order in the tail is " << tailOrder << ", above 2. That is not "
	            "a better Jacobian; it means the triple straddles the round-off floor and the "
	            "number is an artefact rather than a rate" );
}
