#define BOOST_TEST_MODULE SafetyFactorSolver
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <memory>
#include <vector>

#include "mfem.hpp"

#include "meq/FluxExtraction.hpp"
#include "meq/GradShafranov.hpp"
#include "meq/Profiles.hpp"
#include "meq/SafetyFactor.hpp"
#include "meq/SafetyFactorSolve.hpp"
#include "meq/Source.hpp"

#include "convergence/ConvergenceHarness.hpp"

/*
 * DRIVING AN EQUILIBRIUM BY `q( psi )` -- ROADMAP.md item 10, the solver half.
 *
 * `SafetyFactorTests` covers the inversion, which is one division per surface
 * and needs no mesh. What is here is the loop around it: `V'` and `< R^-2 >`
 * are functionals of the solution and the solution depends on `g`, so
 *
 *     solve -> extract the surfaces -> invert -> rebuild gg' -> solve
 *
 * and the fixed point is an equilibrium whose own `q` is the one asked for.
 *
 * THIS IS A PICARD ITERATION ON THE GEOMETRY AND IT IS DELIBERATELY THE HALFWAY
 * HOUSE. Making it quadratic needs `d( V', < R^-2 > )/dpsi`, which
 * INVERSION-PLAN.md section 11.1 measures at about **5.7 hours per Jacobian**
 * by differencing -- `nFieldDOF` complete re-extractions of the whole family --
 * so the Newton wants the shape derivative first and that is not built. The
 * same staging XP-2 takes against XP-3, and for the same reason: an outer fixed
 * point that converges is worth having, and it is what says a later Newton
 * would be a change of ALGORITHM rather than a change of PROBLEM.
 *
 * THE ACCEPTANCE IS A ROUND TRIP WITH A KNOWN ANSWER, which is available here
 * and is worth more than a self-consistency check. The reference equilibrium is
 * built from a `g( Psi )` given in CLOSED FORM, so the `q` it produces is known
 * exactly rather than measured; the loop is then started from a DIFFERENT `g`
 * and has to come back. A loop that merely reached a fixed point would satisfy
 * a self-consistency test while sitting on the wrong one.
 */
namespace
{
	using meq::GradShafranovSolver;

	/*
	 * THE REFERENCE FIELD, LINEAR IN THE SOURCE'S Psi AND SO EXACT IN A KNOT
	 * PAIR. g = g0 + g1 Psi gives gg' = g1( g0 + g1 Psi ), which is linear, and
	 * a Hermite cubic reproduces a cubic exactly -- so the profile the solver
	 * sees IS the closed form and not an approximation of it. That is what lets
	 * the target q be computed rather than measured.
	 */
	double const g0 = 2.20;
	double const g1 = 0.55;

	double referenceG( double psiNormalisedSource )
	{
		return g0 + g1*psiNormalisedSource;
	}

	std::shared_ptr<meq::Profile const> linearGGPrime( double a, double b )
	{
		// gg'( Psi ) = b( a + b Psi ), derivative b^2.
		std::vector<meq::Knot> knots;
		for ( double x : { 0.0, 1.0 } )
			knots.push_back( meq::Knot{ x, b*( a + b*x ), b*b } );
		return std::make_shared<meq::SplineProfile const>( std::move( knots ) );
	}

	meq::FluxFamilyOptions familyOptions()
	{
		meq::FluxFamilyOptions options;
		options.surfaces = 12;
		options.angles = 96;
		return options;
	}

	/// One solve on the standard box with the given gg', to a converged
	/// bordered Newton. Everything but the toroidal profile is held.
	struct Solved
	{
		std::unique_ptr<mfem::Mesh> mesh;
		std::unique_ptr<meq::NormalisedMHDSource> source;
		std::unique_ptr<GradShafranovSolver> solver;
		double psiAxis = 0.0;
		bool converged = true;
	};

	/*
	 * ONE SOLVE, WARM-STARTED WHERE THERE IS SOMETHING TO START FROM.
	 *
	 * THE MESH IS THE CALLER'S AND OUTLIVES EVERY SOLVE, which matters: a
	 * warm start hands the next solver a GridFunction owned by the previous
	 * one, so the previous solver has to still be alive when the next is
	 * built. CLAUDE.md records a SubMesh dangling on its parent for months
	 * before anything asked; this is the same contract one level down.
	 *
	 * AND THE COLD START IS WHAT FAILS. Restarting from the bump at every
	 * iteration does not converge past the second: the profile has moved away
	 * from the closed form the bump was chosen for, and a cold Newton on a
	 * moving normalised source is exactly the case this tree records
	 * PicardThenNewton existing for. Carrying the answer forward is not an
	 * optimisation here, it is what makes the loop run at all.
	 */
	Solved solveWith( mfem::Mesh &mesh,
	                  std::shared_ptr<meq::Profile const> const &ggPrime,
	                  int order, mfem::GridFunction const *previous,
	                  double psiAxisGuess )
	{
		Solved out;
		meq::tests::Rectangle const box = meq::tests::standardBox();

		auto pPrime = std::make_shared<meq::ConstantProfile const>( 0.45 );
		out.source = std::make_unique<meq::NormalisedMHDSource>(
			pPrime, ggPrime, 1.0, 1.0 );

		out.solver = std::make_unique<GradShafranovSolver>( mesh, order );

		mfem::FunctionCoefficient guess(
			[ box ]( mfem::Vector const &x )
			{
				return 0.30*std::sin( M_PI*( x( 0 ) - box.minRadius )/box.width() )
				       *std::sin( M_PI*( x( 1 ) - box.zMin )/box.height() );
			} );
		mfem::ConstantCoefficient zero( 0.0 );

		out.solver->setSource( *out.source, psiAxisGuess );
		out.solver->setBoundaryData( zero );
		if ( previous != nullptr )
			out.solver->setInitialGuess( *previous );
		else
			out.solver->setInitialGuess( guess );
		out.solver->setNewtonControl( 1.0e-12, 1.0e-14, 60 );

		// THE REACTIVE LADDER IS NOT AVAILABLE HERE, which is a real constraint
		// on this loop rather than an omission: psi_ax is a BORDER UNKNOWN, and
		// GradShafranovSolver refuses every Globalisation but None on that path --
		// the KINSOL routes drive a residual of their own and the Picard ones
		// build no Jacobian to border. So a hard intermediate profile has no
		// fallback, and the outer step is the only thing that can keep the inner
		// solve in its basin. That is why the relaxation below is adaptive.
		try
		{
			out.solver->solve();
		}
		catch ( std::exception const & )
		{
			out.converged = false;
			return out;
		}
		out.solver->postProcess();
		out.psiAxis = out.solver->psiAxis();
		return out;
	}
}

/*
 * THE LOOP RECOVERS THE FIELD ITS TARGET WAS BUILT FROM.
 *
 * The reference equilibrium is built from a g( Psi ) given in CLOSED FORM, so
 * the q it produces is known rather than measured; the loop is then started
 * from a g 40% larger everywhere -- a different equilibrium, not a
 * perturbation -- and has to come back. A loop that merely reached A fixed
 * point would satisfy a self-consistency check while sitting on the wrong one.
 *
 * THREE THINGS HAD TO BE TRUE AT ONCE AND EACH FAILED FIRST.
 *
 * A DAMPED PICARD CANNOT DO IT, AND THAT IS A THEOREM. The relaxed iteration
 * has derivative 1 + w( G' - 1 ) at the fixed point, above one for every w > 0
 * when G' > 1. Measured, the damped loop walks the error from 0.362 down to
 * 0.019 and straight back up to 0.50, the step never shrinking -- including
 * where the error passes through zero. It does not stall at the fixed point, it
 * crosses it.
 *
 * THE MAP MUST BE TOTAL, because KINSOL is C. An exception raised inside the
 * residual unwinds through its frames and denies the line search the one thing
 * it needs: a finite value at the trial point. With the map throwing, the
 * Newton's first full step lands on g^2 < 0 across the whole profile and the
 * run dies without ever backtracking.
 *
 * AND THE KNOTS MUST SPAN THE WHOLE OF Psi, WHICH IS THE ONE THAT LOOKED LIKE
 * PHYSICS AND WAS ARITHMETIC. Laying them only over the family's own range and
 * letting SplineProfile clamp gives a profile that differs from the closed form
 * exactly where it clamps -- so the outer residual AT THE ANSWER was 1.3e-01
 * rather than zero, and both outer methods converged, correctly, to a fixed
 * point that was not it. theOuterResidualsConditioningAtTheAnswer is the
 * measurement that found it and now guards it: 4.9e-13, on a Jacobian whose
 * degeneracy measure is 0.12.
 */
BOOST_AUTO_TEST_CASE( theLoopRecoversTheFieldItsTargetWasBuiltFrom )
{
	int const n = 24;
	int const order = 2;

	// ---- the reference, and the target read off it -------------------------
	meq::tests::Rectangle const box = meq::tests::standardBox();
	mfem::Mesh mesh = meq::tests::makeMesh( box, n );

	Solved const reference =
		solveWith( mesh, linearGGPrime( g0, g1 ), order, nullptr, 0.30 );

	meq::FluxFamilyOptions options = familyOptions();
	double const referenceAxis = reference.psiAxis;
	double const referenceSpan = referenceAxis - reference.solver->psiBoundary();
	options.toroidalField = [ referenceAxis, referenceSpan ]( double psi )
	{
		return referenceG( ( psi - ( referenceAxis - referenceSpan ) )
		                   /referenceSpan );
	};

	meq::FluxSurfaceFamily const referenceFamily =
		meq::extractFluxSurfaces( *reference.solver, options );
	BOOST_TEST_REQUIRE( referenceFamily.safetyFactorAvailable,
		"the reference family carries no safety factor, so there is no target" );

	std::vector<double> const label = [ & ]
	{
		std::vector<double> out;
		for ( meq::FluxSurface const &s : referenceFamily.surfaces )
			out.push_back( s.normalisedFlux );
		return out;
	}();
	std::vector<double> const targetQ = [ & ]
	{
		std::vector<double> out;
		for ( meq::FluxSurface const &s : referenceFamily.surfaces )
			out.push_back( s.safetyFactor );
		return out;
	}();

	// The target is a function of the family's OWN label, and every family in
	// the loop uses the same cut, count and spacing -- so the grid does not
	// move and this is an interpolation only against round-off.
	auto target = [ & ]( double psiN )
	{
		std::size_t best = 0;
		for ( std::size_t i = 1; i < label.size(); ++i )
			if ( std::abs( label[ i ] - psiN ) < std::abs( label[ best ] - psiN ) )
				best = i;
		return targetQ[ best ];
	};

	std::printf( "\n  DRIVING BY q( psi ): the loop against a known g\n" );
	std::printf( "    reference g = %.3f + %.3f Psi, psi_ax %.6e\n",
	             g0, g1, referenceAxis );
	/*
	 * ---- THE OUTER NEWTON ON g^2's COEFFICIENTS ---------------------------
	 *
	 * The unknown is the handful of coefficients the profile is fitted in, and
	 * the residual is G( c ) - c where G is one whole Picard step: build gg'
	 * from c, solve, extract the surfaces, invert the target against them, and
	 * fit what that implies. Its root is the fixed point
	 * theReferenceIsAFixedPointOfItsOwnTarget shows is the answer.
	 */
	unsigned int const degree = 2;

	auto fieldFromCoefficients = [ & ]( std::vector<double> const &c )
	{
		meq::ToroidalField field;
		field.normalisedFlux = label;
		for ( double x : label )
		{
			double const psi = 1.0 - x;
			double value = 0.0, power = 1.0;
			for ( double coefficient : c )
			{
				value += coefficient*power;
				power *= psi;
			}
			field.gSquared.push_back( value );
			field.g.push_back( value > 0.0 ? std::sqrt( value ) : 0.0 );
			field.safetyFactor.push_back( 0.0 );
		}
		return field;
	};

	std::unique_ptr<Solved> current;
	int solves = 0;

	/*
	 * WHAT THE MAP RETURNS WHERE IT CANNOT BE EVALUATED, AND WHY IT MUST NOT
	 * THROW.
	 *
	 * KINSOL is C. An exception raised inside mfem::Operator::Mult unwinds
	 * through its frames, which is undefined behaviour and in any case denies
	 * the line search the one thing it needs -- a finite residual at the trial
	 * point, so that Armijo can reject the step and halve it. Measured: with a
	 * throw here the Newton's first full step lands on g^2 < 0 across the whole
	 * profile and the run dies, never having backtracked once.
	 *
	 * So an inadmissible c returns a residual that points back to the last
	 * state that WAS admissible. It is large there and vanishes nowhere, which
	 * is what a merit function needs; and it can introduce no spurious fixed
	 * point, because g^2 < 0 is not an equilibrium and the normal branch runs
	 * everywhere it is.
	 */
	std::vector<double> lastGood;

	auto unreachable = [ & ]( std::vector<double> const &c )
	{
		std::vector<double> out( c.size() );
		for ( std::size_t i = 0; i < c.size(); ++i )
			out[ i ] = lastGood.empty() ? c[ i ] - 1.0 : lastGood[ i ];
		return out;
	};

	auto picardStep = [ & ]( std::vector<double> const &c )
	{
		++solves;
		{
			meq::ToroidalField const trial = fieldFromCoefficients( c );
			for ( double value : trial.gSquared )
				if ( !( value > 0.0 ) || !std::isfinite( value ) )
					return unreachable( c );
		}

		Solved attempt = solveWith(
			mesh, std::make_shared<meq::SplineProfile const>(
				meq::ggPrimeKnotsFromFit( fieldFromCoefficients( c ), degree ) ),
			order, current ? &current->solver->potential() : nullptr,
			current ? current->psiAxis : 0.30 );

		if ( !attempt.converged )
			return unreachable( c );

		meq::FluxFamilyOptions loopOptions = familyOptions();
		meq::FluxSurfaceFamily const family =
			meq::extractFluxSurfaces( *attempt.solver, loopOptions );
		meq::ToroidalField const inverted =
			meq::invertSafetyFactor( family, target );

		// The BASE iterate's solve is what every column is warm-started from,
		// so it is kept and the perturbed ones are not.
		lastGood = c;
		std::vector<double> const out =
			meq::fitToroidalFieldSquared( inverted, degree );
		return out;
	};

	// Start 40% away, as the damped loop did, so the two are comparable.
	double const startScale = 1.40;
	std::vector<double> startCoefficients;
	{
		meq::ToroidalField warm;
		warm.normalisedFlux = label;
		for ( double x : label )
		{
			double const g = startScale*referenceG( 1.0 - x );
			warm.g.push_back( g );
			warm.gSquared.push_back( g*g );
			warm.safetyFactor.push_back( 0.0 );
		}
		startCoefficients = meq::fitToroidalFieldSquared( warm, degree );
	}

	// One solve at the start, kept as the warm start for everything after.
	current = std::make_unique<Solved>( solveWith(
		mesh, std::make_shared<meq::SplineProfile const>(
			meq::ggPrimeKnotsFromFit( fieldFromCoefficients( startCoefficients ),
			                          degree ) ),
		order, nullptr, 0.30 ) );
	BOOST_TEST_REQUIRE( current->converged,
		"the loop's own starting profile does not solve" );

	meq::OuterNewtonOptions newtonOptions;
	newtonOptions.functionTolerance = 1.0e-8;
	newtonOptions.maxIterations = 20;
	meq::OuterNewtonResult const result =
		meq::solveForToroidalField( startCoefficients, picardStep,
		                            newtonOptions );

	std::printf( "    outer Newton: %s\n", result.status.c_str() );
	std::printf( "    %d inner solves in all\n", solves );

	BOOST_TEST_REQUIRE( result.converged,
		"the outer Newton did not converge: " << result.status );

	// The recovered g against the closed form its target was built from.
	meq::ToroidalField const recovered =
		fieldFromCoefficients( result.coefficients );
	double worstError = 0.0;
	for ( std::size_t i = 0; i < label.size(); ++i )
	{
		double const exact = referenceG( 1.0 - label[ i ] );
		worstError = std::max( worstError,
		                       std::abs( recovered.g[ i ] - exact )
		                       /std::abs( exact ) );
	}

	// And the equilibrium the answer implies.
	Solved const final = solveWith(
		mesh, std::make_shared<meq::SplineProfile const>(
			meq::ggPrimeKnotsFromFit( recovered, degree ) ),
		order, &current->solver->potential(), current->psiAxis );
	BOOST_TEST_REQUIRE( final.converged );
	double const lastAxis = final.psiAxis;

	std::printf( "    worst |g - g_exact|/g  %.6e\n", worstError );
	std::printf( "    psi_ax %.6e against the reference's %.6e (%.2e)\n",
	             lastAxis, referenceAxis,
	             std::abs( lastAxis - referenceAxis )/referenceAxis );

	BOOST_TEST( worstError < 1.0e-4,
		"the Newton settled with g still " << worstError << " from the closed "
		"form its own target was built from, having started 0.40 away" );
	BOOST_TEST( std::abs( lastAxis - referenceAxis ) < 1.0e-4*referenceAxis,
		"the recovered equilibrium's psi_ax is " << lastAxis
		<< " against the reference's " << referenceAxis );
}

/*
 * THE REFERENCE IS A FIXED POINT OF ITS OWN TARGET, WHICH IS THE ONE THING THE
 * LOOP ABOVE CANNOT ESTABLISH FOR ITSELF.
 *
 * Extract the reference family, read its q, and invert with that q as the
 * target: the g that comes back must be the g the family was built with, to the
 * extraction's own round-off. If it is not, the loop is chasing a fixed point
 * that is not the answer, and no amount of relaxation or step control would
 * ever find one -- which is a DIFFERENT failure from an unstable map and wants a
 * different repair.
 *
 * This is the diagnostic that separates the two, and it is cheap: one solve.
 */
BOOST_AUTO_TEST_CASE( theReferenceIsAFixedPointOfItsOwnTarget )
{
	int const n = 24;
	int const order = 2;

	meq::tests::Rectangle const box = meq::tests::standardBox();
	mfem::Mesh mesh = meq::tests::makeMesh( box, n );
	Solved const reference =
		solveWith( mesh, linearGGPrime( g0, g1 ), order, nullptr, 0.30 );
	BOOST_TEST_REQUIRE( reference.converged );

	double const axis = reference.psiAxis;
	double const span = axis - reference.solver->psiBoundary();

	meq::FluxFamilyOptions options = familyOptions();
	options.toroidalField = [ axis, span ]( double psi )
	{
		return referenceG( ( psi - ( axis - span ) )/span );
	};
	meq::FluxSurfaceFamily const family =
		meq::extractFluxSurfaces( *reference.solver, options );
	BOOST_TEST_REQUIRE( family.safetyFactorAvailable );

	std::vector<double> label, targetQ;
	for ( meq::FluxSurface const &s : family.surfaces )
	{
		label.push_back( s.normalisedFlux );
		targetQ.push_back( s.safetyFactor );
	}

	meq::ToroidalField const inverted = meq::invertSafetyFactor(
		family, [ & ]( double psiN )
		{
			std::size_t best = 0;
			for ( std::size_t i = 1; i < label.size(); ++i )
				if ( std::abs( label[ i ] - psiN )
				     < std::abs( label[ best ] - psiN ) )
					best = i;
			return targetQ[ best ];
		} );

	std::printf( "\n  IS THE REFERENCE A FIXED POINT OF ITS OWN TARGET?\n" );
	std::printf( "    %10s %14s %14s %12s\n",
	             "Psi_N", "g inverted", "g exact", "relative" );
	double worst = 0.0;
	for ( std::size_t i = 0; i < inverted.size(); ++i )
	{
		double const exact = referenceG( 1.0 - label[ i ] );
		double const relative = std::abs( inverted.g[ i ] - exact )
		                        /std::abs( exact );
		worst = std::max( worst, relative );
		std::printf( "    %10.4f %14.6f %14.6f %12.3e\n",
		             label[ i ], inverted.g[ i ], exact, relative );
	}
	std::printf( "    worst %.3e\n", worst );

	BOOST_TEST( worst < 1.0e-10,
		"inverting the reference's OWN q against the reference's OWN geometry "
		"returns a g that is " << worst << " from the one it was built with. "
		"The two are the same division, so this is not the geometry's accuracy "
		"-- it is an inconsistency between how the target was read and how it "
		"is inverted, and it means the loop's fixed point is not the answer." );
}

/*
 * NEWTON FINISHES AN AFFINE MAP IN ONE STEP -- INCLUDING AN UNSTABLE ONE.
 *
 * The acceptance FB-5 uses for the same reason: if the residual is affine in the
 * unknown then an EXACT Jacobian must land on the root in a single step, so one
 * step is a statement about the Jacobian rather than about the problem being
 * easy. Here the Jacobian is differenced rather than assembled, so it is exact
 * only to the difference's own truncation -- which on an affine map is zero,
 * because a central difference of a linear function is that function's slope
 * with no remainder at all.
 *
 * THE MAP IS DELIBERATELY UNSTABLE. Its multiplier is 1.8, so the fixed point
 * REPELS and the Picard iteration this replaces runs away from it. That is the
 * case the outer Newton exists for, and the control below is what says damping
 * could not have done it.
 */
BOOST_AUTO_TEST_CASE( newton_finishes_an_unstable_affine_map_in_one_step )
{
	std::vector<double> const fixedPoint = { 4.84, 2.42, 0.3025 };
	double const multiplier = 1.8;

	auto map = [ & ]( std::vector<double> const &c )
	{
		std::vector<double> out( c.size() );
		for ( std::size_t i = 0; i < c.size(); ++i )
			out[ i ] = fixedPoint[ i ] + multiplier*( c[ i ] - fixedPoint[ i ] );
		return out;
	};

	std::vector<double> start = { 6.0, 3.4, 0.5 };
	meq::OuterNewtonResult const result =
		meq::solveForToroidalField( start, map );

	BOOST_TEST( result.converged,
		"Newton did not converge on an AFFINE map, where an exact Jacobian must "
		"land on the root in one step" );
	BOOST_TEST( result.iterations == 1,
		"Newton took " << result.iterations << " steps on an affine map. One is "
		"the whole assertion: more means the differenced Jacobian is not the "
		"map's own slope." );
	for ( std::size_t i = 0; i < fixedPoint.size(); ++i )
		BOOST_TEST( std::abs( result.coefficients[ i ] - fixedPoint[ i ] )
		            < 1.0e-10,
			"coefficient " << i << " came back " << result.coefficients[ i ]
			<< " against the map's own fixed point " << fixedPoint[ i ] );

	/*
	 * THE CONTROL, AND IT IS THE REASON THE NEWTON IS HERE. The same map under
	 * a RELAXED Picard iteration, at three dampings spanning two decades. Every
	 * one of them moves AWAY: the relaxed derivative is 1 + w( G' - 1 ), which
	 * for G' = 1.8 is 1 + 0.8 w, above one for every w > 0. Under-relaxation
	 * stabilises a map that oscillates and can do nothing for one that runs
	 * away, and that is a theorem rather than a measurement -- this column is
	 * the measurement that it applies here.
	 */
	for ( double omega : { 1.0, 0.5, 0.01 } )
	{
		std::vector<double> c = start;
		double const before = std::abs( c[ 0 ] - fixedPoint[ 0 ] );
		for ( int step = 0; step < 200; ++step )
		{
			std::vector<double> const image = map( c );
			for ( std::size_t i = 0; i < c.size(); ++i )
				c[ i ] += omega*( image[ i ] - c[ i ] );
		}
		double const after = std::abs( c[ 0 ] - fixedPoint[ 0 ] );
		BOOST_TEST( after > before,
			"damped Picard at omega = " << omega << " moved the error from "
			<< before << " to " << after << ", i.e. TOWARD the fixed point. If "
			"a relaxation can reach this map then it is not the unstable one "
			"this control is meant to be, and the case above proves nothing "
			"the Picard loop could not have done." );
	}
}

/*
 * A MAP WHOSE FIXED POINT IS NOT DETERMINED IS REFUSED RATHER THAN SOLVED, AND
 * PROMPTLY -- WHICH IS THE HALF THAT WAS BROKEN.
 *
 * A rank-deficient outer Jacobian is a real state: it means the degree `g^2` is
 * fitted in is higher than the surface family determines, so some combination
 * of coefficients does not change the equilibrium the loop comes back with.
 * KINSOL's contract is to REPORT that rather than throw, and this asserts the
 * report.
 *
 * **WHAT IT REALLY GUARDS IS THAT THE ANSWER ARRIVES AT ALL.** `KIN_LINESEARCH`
 * interpolates its step length on a quotient whose numerator and denominator
 * both carry the directional derivative `< F, J p >`; a singular `J` makes that
 * zero whatever the step is, so the quotient is `0/0`, the iterate goes to NaN,
 * and **KINSOL never returns** -- every one of its stopping tests is a
 * comparison against NaN and every comparison against NaN is false. This case
 * ran for **36 minutes** without finishing, and the map was still being called,
 * at NaN, after two million evaluations. On the real loop each of those is an
 * equilibrium solve.
 *
 * So the assertion below is on the EVALUATION COUNT, not on a wall clock: a
 * timing here would be a measurement about the machine where the count is a
 * measurement about the code. The same currency, and the same reason, as the
 * symbolic-factorisation reuse case.
 *
 * **AND THE IDENTITY MAP IS THE CONTROL RATHER THAN THE EXAMPLE.** Its Jacobian
 * is exactly zero too, and it is the opposite answer -- every point is a fixed
 * point of it, so returning at once is right. The two are indistinguishable by
 * their Jacobians, so a refusal keyed on rank alone refuses both; the residual
 * is what separates them. They are asserted side by side for that reason.
 */
BOOST_AUTO_TEST_CASE( an_undetermined_outer_jacobian_is_refused )
{
	std::vector<double> const start = { 1.0, 2.0 };

	// THE IDENTITY IS NOT THE DEGENERATE CASE, WHICH IS WORTH SAYING BECAUSE IT
	// LOOKS LIKE IT. Every point is a fixed point of it, so the residual at the
	// start is already zero and returning immediately is the RIGHT answer --
	// there is nothing to solve, and KINSOL is never entered.
	auto identity = []( std::vector<double> const &c ) { return c; };
	meq::OuterNewtonResult const trivial =
		meq::solveForToroidalField( start, identity );
	BOOST_TEST( trivial.converged );
	BOOST_TEST( trivial.iterations == 0 );

	// ITS JACOBIAN IS SINGULAR TOO, AND THAT IS WHY THE RESIDUAL HAS TO BE THE
	// DISCRIMINATOR. The identity and the no-root map below are
	// INDISTINGUISHABLE by their Jacobians -- both exactly zero -- and they are
	// opposite answers. A refusal keyed on rank alone would refuse this one.
	BOOST_TEST( trivial.jacobianConditioning == 0.0 );

	// THE DEGENERATE CASE IS A NONZERO RESIDUAL WITH NO SLOPE UNDER IT: this
	// map has no fixed point at all, so R is a nonzero constant and its
	// Jacobian is the zero matrix. That is what has to be refused rather than
	// solved, because a singular solve returns a step and the caller would
	// never know.
	auto noFixedPoint = []( std::vector<double> const &c )
	{
		std::vector<double> out( c );
		for ( double &value : out )
			value += 0.1;
		return out;
	};
	meq::OuterNewtonResult const hopeless =
		meq::solveForToroidalField( start, noFixedPoint );
	BOOST_TEST( !hopeless.converged,
		"a map with no fixed point was reported as converged: "
		<< hopeless.status );
	BOOST_TEST( hopeless.jacobianConditioning == 0.0,
		"the Jacobian of a constant residual is the zero matrix, so its "
		"smallest singular value must be exactly zero, not "
		<< hopeless.jacobianConditioning );

	// **AND IT MUST BE REFUSED PROMPTLY, WHICH IS THE ASSERTION WITH TEETH.**
	// Before the degeneracy was caught ahead of KINSOL this case did not fail,
	// it HUNG: KIN_LINESEARCH interpolates on a quotient carrying < F, J p >,
	// which a singular J makes 0/0, and once the iterate is NaN every one of
	// KINSOL's stopping tests is a comparison against NaN and so false.
	// Measured, the map was still being called after two million evaluations.
	//
	// A COUNT AND NOT A TIMING, for the reason this tree gives for the
	// symbolic-factorisation reuse: a wall clock here would be a measurement
	// about the machine, where the evaluation count is a measurement about the
	// code -- and on the real loop every one of those evaluations is a whole
	// equilibrium solve. One residual and 2n differencing calls is 5 here.
	BOOST_TEST( hopeless.mapEvaluations <= 4*static_cast<int>( start.size() ),
		"a map with no root cost " << hopeless.mapEvaluations
		<< " evaluations, which is more than forming one Jacobian and its "
		"residual -- the degeneracy is being discovered inside the iteration "
		"rather than before it" );
	std::printf( "\n  a map with no root: %d evaluations, %s\n",
	             hopeless.mapEvaluations, hopeless.status.c_str() );

	// THE CONTROL. A map that DOES determine its fixed point must be accepted,
	// or the refusal above is compatible with a Newton that refuses everything.
	auto determined = []( std::vector<double> const &c )
	{
		std::vector<double> out( c.size() );
		for ( std::size_t i = 0; i < c.size(); ++i )
			out[ i ] = 0.5*c[ i ] + 1.0;
		return out;
	};
	BOOST_CHECK_NO_THROW( meq::solveForToroidalField( start, determined ) );

	// And a map that changes the length of its argument is a different space.
	auto ragged = []( std::vector<double> const &c )
	{
		return std::vector<double>( c.size() + 1, 0.0 );
	};
	BOOST_CHECK_THROW( meq::solveForToroidalField( start, ragged ),
	                   std::invalid_argument );
}

/*
 * HOW WELL CONDITIONED IS THE ROOT? -- the measurement that decides what an
 * outer solver can do here at all.
 *
 * theReferenceIsAFixedPointOfItsOwnTarget shows the answer IS a root of
 * R( c ) = G( c ) - c, to 2.1e-14 -- but in the SURFACE values, and the outer
 * solver does not iterate in those. This measures the residual in the fitted
 * COEFFICIENTS, which is the space it does iterate in, and differences R'
 * there by the same central differences the Newton uses.
 *
 * **THE ASSERTION IS THAT R( c_exact ) IS ZERO, AND IT IS THE SHARPEST THING
 * IN THIS FILE, BECAUSE IT IS WHAT FAILED.** With ggPrimeKnotsFromFit()
 * sampling only the family's own range -- `[ 0.05, 0.95 ]` rather than the
 * whole of `Psi` -- the residual at the exact answer read **1.3e-01**, so the
 * closed-form `g` was NOT a fixed point of the loop built around it and
 * nothing the outer solver found could have been right. Both outer methods
 * then converged, from 40% away and from 5% away, to the same wrong root, and
 * every symptom looked like a conditioning problem: a damped Picard drifting
 * through at constant speed and a Newton landing somewhere else twice over are
 * exactly what a near-singular R' would produce. **The residual at a known
 * answer is what separates the two, and it is available here only because the
 * reference `g` is a closed form.** Spanning the knots over `[ 0, 1 ]` takes
 * it to 4.9e-13.
 *
 * The conditioning numbers are printed rather than asserted on: the degeneracy
 * measure is a property of the fixture's degree and surface count, and what
 * this case is for is the root.
 */
BOOST_AUTO_TEST_CASE( theOuterResidualsConditioningAtTheAnswer )
{
	int const n = 24;
	int const order = 2;
	unsigned int const degree = 2;

	meq::tests::Rectangle const box = meq::tests::standardBox();
	mfem::Mesh mesh = meq::tests::makeMesh( box, n );
	Solved reference =
		solveWith( mesh, linearGGPrime( g0, g1 ), order, nullptr, 0.30 );
	BOOST_TEST_REQUIRE( reference.converged );

	double const axis = reference.psiAxis;
	double const span = axis - reference.solver->psiBoundary();
	meq::FluxFamilyOptions options = familyOptions();
	options.toroidalField = [ axis, span ]( double psi )
	{
		return referenceG( ( psi - ( axis - span ) )/span );
	};
	meq::FluxSurfaceFamily const referenceFamily =
		meq::extractFluxSurfaces( *reference.solver, options );

	std::vector<double> label, targetQ;
	for ( meq::FluxSurface const &s : referenceFamily.surfaces )
	{
		label.push_back( s.normalisedFlux );
		targetQ.push_back( s.safetyFactor );
	}
	auto target = [ & ]( double psiN )
	{
		std::size_t best = 0;
		for ( std::size_t i = 1; i < label.size(); ++i )
			if ( std::abs( label[ i ] - psiN ) < std::abs( label[ best ] - psiN ) )
				best = i;
		return targetQ[ best ];
	};

	auto fieldFromCoefficients = [ & ]( std::vector<double> const &c )
	{
		meq::ToroidalField field;
		field.normalisedFlux = label;
		for ( double x : label )
		{
			double const psi = 1.0 - x;
			double value = 0.0, power = 1.0;
			for ( double coefficient : c )
			{
				value += coefficient*power;
				power *= psi;
			}
			field.gSquared.push_back( value );
			field.g.push_back( std::sqrt( std::max( value, 0.0 ) ) );
			field.safetyFactor.push_back( 0.0 );
		}
		return field;
	};

	auto residual = [ & ]( std::vector<double> const &c )
	{
		Solved const step = solveWith(
			mesh, std::make_shared<meq::SplineProfile const>(
				meq::ggPrimeKnotsFromFit( fieldFromCoefficients( c ), degree ) ),
			order, &reference.solver->potential(), reference.psiAxis );
		BOOST_TEST_REQUIRE( step.converged,
			"the solve failed at a point one difference step from the answer" );
		meq::FluxFamilyOptions loopOptions = familyOptions();
		meq::ToroidalField const inverted = meq::invertSafetyFactor(
			meq::extractFluxSurfaces( *step.solver, loopOptions ), target );
		std::vector<double> out =
			meq::fitToroidalFieldSquared( inverted, degree );
		for ( std::size_t i = 0; i < out.size(); ++i )
			out[ i ] -= c[ i ];
		return out;
	};

	// g^2 = ( g0 + g1 Psi )^2, exactly quadratic, so these ARE the answer.
	std::vector<double> const exact = { g0*g0, 2.0*g0*g1, g1*g1 };
	std::vector<double> const atAnswer = residual( exact );

	std::printf( "\n  THE OUTER RESIDUAL AT THE ANSWER\n" );
	std::printf( "    c_exact  [ %.6f %.6f %.6f ]\n",
	             exact[ 0 ], exact[ 1 ], exact[ 2 ] );
	std::printf( "    R        [ %.3e %.3e %.3e ]\n",
	             atAnswer[ 0 ], atAnswer[ 1 ], atAnswer[ 2 ] );

	std::size_t const m = exact.size();
	std::vector<std::vector<double>> columns( m );
	for ( std::size_t j = 0; j < m; ++j )
	{
		double const h = 1.0e-4*std::max( std::abs( exact[ j ] ), 1.0 );
		std::vector<double> plus = exact, minus = exact;
		plus[ j ] += h;
		minus[ j ] -= h;
		std::vector<double> const rPlus = residual( plus );
		std::vector<double> const rMinus = residual( minus );
		columns[ j ].resize( m );
		for ( std::size_t i = 0; i < m; ++i )
			columns[ j ][ i ] = ( rPlus[ i ] - rMinus[ i ] )/( 2.0*h );
	}

	// A 3x3 by hand rather than a matrix library: Eigen is PRIVATE to meq_core
	// and a test does not see it. The degeneracy measure is |det| against the
	// product of the row norms, which is 1 for orthogonal rows and goes to zero
	// as they become dependent -- enough to say whether R' has slope under it,
	// which is the whole question.
	double j[ 3 ][ 3 ];
	for ( std::size_t i = 0; i < m; ++i )
		for ( std::size_t k = 0; k < m; ++k )
			j[ i ][ k ] = columns[ k ][ i ];

	std::printf( "    R' =\n" );
	for ( std::size_t i = 0; i < m; ++i )
		std::printf( "        [ %12.5e %12.5e %12.5e ]\n",
		             j[ i ][ 0 ], j[ i ][ 1 ], j[ i ][ 2 ] );

	double const determinant =
		j[0][0]*( j[1][1]*j[2][2] - j[1][2]*j[2][1] )
		- j[0][1]*( j[1][0]*j[2][2] - j[1][2]*j[2][0] )
		+ j[0][2]*( j[1][0]*j[2][1] - j[1][1]*j[2][0] );

	double rowNorms = 1.0;
	for ( std::size_t i = 0; i < m; ++i )
		rowNorms *= std::sqrt( j[i][0]*j[i][0] + j[i][1]*j[i][1]
		                       + j[i][2]*j[i][2] );

	std::printf( "    det %.5e, |det|/(product of row norms) %.5e\n",
	             determinant, std::abs( determinant )/rowNorms );

	// THE ROOT IS A ROOT. This restates theReferenceIsAFixedPointOfItsOwnTarget
	// in the coefficients the outer solver actually iterates in, which is the
	// space the conditioning below is measured in.
	double worst = 0.0;
	for ( double value : atAnswer )
		worst = std::max( worst, std::abs( value ) );
	BOOST_TEST( worst < 1.0e-8,
		"the residual at the closed-form answer is " << worst
		<< ", so the coefficients the outer solver iterates in do not have the "
		"answer as a root and nothing it finds could be right" );
}
