#define BOOST_TEST_MODULE SafetyFactorSolver
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <functional>
#include <memory>
#include <stdexcept>
#include <string>
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
		out.solver->setNewtonControl( 1.0e-12,
		                  1.0e-14, 60 );

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
 * THE MAP MUST BE TOTAL. An exception raised inside the residual ends the loop
 * -- KINSOL is C, so the callback catches it and the loop reports failure --
 * and denies the line search the one thing it needs: a finite value at the
 * trial point. With the map throwing, the Newton's first full step lands on
 * g^2 < 0 across the whole profile and the run ends without ever
 * backtracking.
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
	 * KINSOL is C, so an exception raised inside mfem::Operator::Mult is
	 * caught in the callback and ends the loop as a reported failure
	 * (a_map_that_throws_ends_the_loop_and_is_reported). That denies the line
	 * search the one thing it needs -- a finite residual at the trial point,
	 * so that Armijo can reject the step and halve it. Measured: with a throw
	 * here the Newton's first full step lands on g^2 < 0 across the whole
	 * profile and the run ends, never having backtracked once.
	 *
	 * So an inadmissible c returns a residual that points back to the last
	 * state that WAS admissible. It is large there and vanishes nowhere, which
	 * is what a merit function needs; and it can introduce no spurious fixed
	 * point, because g^2 < 0 is not an equilibrium and the normal branch runs
	 * everywhere it is.
	 */
	std::vector<double> lastGood;
	double mapFloor = 0.0;

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
		// THE MAP'S OWN FLOOR, from the accuracy every solve here is converged
		// to at worst: solveWith()'s relative tolerance floored at ten eps,
		// the gate under which a single solve stops at its round-off floor
		// (M-306).
		mapFloor = meq::toroidalMapFloor( inverted, degree,
		                                  meq::newtonTolerance( 1.0e-12 ) );
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
	newtonOptions.mapFloor = [ &mapFloor ]() { return mapFloor; };
	meq::OuterNewtonResult const result =
		meq::solveForToroidalField( startCoefficients, picardStep,
		                            newtonOptions );

	std::printf( "    outer Newton: %s\n", result.status.c_str() );
	std::printf( "    %d inner solves in all\n", solves );

	BOOST_TEST_REQUIRE( result.converged,
		"the outer Newton did not converge: " << result.status );

	/*
	 * **A NEWTON AND NOT A CRAWL.** 4 KINSOL iterations and 29 solves in
	 * double, 2 or 3 in single, the last step quadratic. The map is piecewise
	 * smooth -- a surface node crossing an element face steps it by the DG
	 * jump there, 2.6e-05 in the coefficients on this fixture -- and a
	 * Jacobian column whose difference straddles one is off by that over the
	 * step: at a step of 4.4e-4 one sits 5.3e-4 below the root's c_1, every
	 * Jacobian near the root straddles it, and the loop contracts linearly at
	 * 0.19 for 8 iterations and 57 solves. MEASUREMENTS.md M-308b. Red here
	 * means a column straddles a face crossing; it is the step, not the map.
	 */
	BOOST_TEST( result.iterations <= 5,
		"the outer Newton took " << result.iterations << " iterations, where "
		"a Newton takes 4: a differenced Jacobian column straddles a step in "
		"the map, a surface node crossing an element face (M-308b). "
		<< result.status );

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

	/*
	 * **AND THE ANSWER IS WORTH THE RESIDUAL THE LOOP CAME TO REST AT**, which
	 * is now the residual it was asked for or the map's floor, whichever is
	 * larger: the loop is converged only there (meq::solveForToroidalField).
	 * Asked for 1e-08, it rests at 1.6e-09 in double -- 4 iterations, 29
	 * solves -- and at 4.5e-05 to 4.6e-04 in single against a floor of
	 * 4.9e-04, in 2 or 3 iterations and 15 or 22 solves on both MKLs and at
	 * one, two and four MKL threads, with g 3.1e-05 to 5.8e-04 out.
	 *
	 * **THE FIGURES THIS USED TO QUOTE WERE KINSOL'S AND NOT THE MAP'S.**
	 * 3.68e-06 in double and 3.71e-03 in single are under MFEM's default
	 * absolute tolerance, eps^( 1/3 ) -- 6.1e-06 and 4.9e-03 -- which is
	 * where the loop stopped, reporting success, before it was given a
	 * tolerance of its own. The single one was read here as the map's eps/h^2
	 * noise; measured, that noise is 1.5 to 4.6 times the inner solve's 10
	 * eps per surface and flat in h, and the single loop's three failures in
	 * six MKL draws ("no sufficiently distinct iterate") were its line search
	 * working on a Jacobian differenced at 1e-4 against it.
	 *
	 * achievedResidualFloor() `std::max`es, so in double it picks the 1e-4 that
	 * has always been asserted. M-222.
	 */
	double const settled =
		meq::tests::achievedResidualFloor( 1.0e-4, result.finalResidual );
	BOOST_TEST( worstError < settled,
		"the Newton settled with g still " << worstError << " from the closed "
		"form its own target was built from, having started 0.40 away, against "
		"a floor of " << settled << " set by the " << result.finalResidual
		<< " its own residual came to rest at" );
	BOOST_TEST( std::abs( lastAxis - referenceAxis )
	            < meq::tests::achievedResidualFloor( 1.0e-4, result.finalResidual )
	              *referenceAxis,
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

	// THE TWO ARE THE SAME DIVISION, so in double this is a bit-level identity
	// and 1e-10 is what it has always been asserted at. In single the same
	// division is carried out on a q the extraction delivers to its own
	// round-off, and the measured disagreement is 1.68e-08 -- 0.14 of ONE eps
	// of the relative scale, which is what an identity looks like in a narrower
	// format. scalarAgreementFloor() is 64 eps and picks 1e-10 in double.
	BOOST_TEST( worst < meq::tests::scalarAgreementFloor( 1.0e-10, 1.0 ),
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

	/*
	 * ONE STEP, AND WHERE IT LANDS.
	 *
	 * **WHERE IT LANDS IS THE DIFFERENCE'S ACCURACY TIMES THE DISTANCE IT
	 * TRAVELLED, NOT THE MAP'S EXACTNESS.** A central difference of a linear
	 * function is that function's slope with no truncation remainder -- which
	 * is the whole reason this fixture is affine -- but the CANCELLATION term
	 * survives, and it is eps/step, so the best a differenced Jacobian can know
	 * its slope to is sqrt( eps ). One Newton step with a Jacobian carrying a
	 * relative error d lands within d times the distance it moved. Measured at
	 * a step of 1e-4: coefficient 0 travels 1.16 and comes back 2.25e-04 out,
	 * coefficient 1 travels 0.98 and comes back 2.08e-05 out, against
	 * sqrt( eps ) = 3.45e-04. M-219, M-222.
	 *
	 * AND IT COSTS ONE RESIDUAL, ONE JACOBIAN AND ONE TRIAL: 2n + 2
	 * evaluations. KINSOL's own difference-quotient `J p` for the line search,
	 * which is what it forms unless it is handed the Jacobian's product, is a
	 * ninth -- an equilibrium per outer iteration on the real loop.
	 */
	meq::OuterNewtonOptions oneStep;
	oneStep.maxIterations = 1;
	meq::OuterNewtonResult const first =
		meq::solveForToroidalField( start, map, oneStep );
	BOOST_TEST( first.iterations == 1 );
	BOOST_TEST( first.mapEvaluations == 2*static_cast<int>( start.size() ) + 2,
		"one Newton step on a three-coefficient map cost "
		<< first.mapEvaluations << " evaluations, against the 8 of a residual, "
		"a Jacobian and a trial point. One more is KINSOL forming J p by a "
		"difference of its own instead of from the Jacobian" );
	for ( std::size_t i = 0; i < fixedPoint.size(); ++i )
	{
		double const travelled = std::abs( start[ i ] - fixedPoint[ i ] );
		double const landing =
			std::max( 1.0e-10,
			          meq::tests::differencedAgreementFloor( 0.0 )*travelled );
		BOOST_TEST( std::abs( first.coefficients[ i ] - fixedPoint[ i ] )
		            < landing,
			"coefficient " << i << " came back " << first.coefficients[ i ]
			<< " against the map's own fixed point " << fixedPoint[ i ]
			<< ", outside the " << landing << " that a differenced Jacobian "
			"leaves after one step over a distance of " << travelled );
	}

	/*
	 * AND THE WHOLE SOLVE, JUDGED AGAINST ITS TOLERANCE.
	 *
	 * **ONE STEP IN DOUBLE, TWO IN SINGLE, AND THE SECOND IS THE FORMAT AND
	 * NOT THE JACOBIAN.** The first step lands within sqrt( eps ) of the
	 * distance it travelled, as above; the loop is converged only once the
	 * residual is at max( tolerance, floor ), and the floor of an EXACT map is
	 * the iterate's own rounding, `|| J || eps || c ||` -- 6e-07 in single and
	 * 1e-15 in double, against a tolerance of 1e-09. Double's first landing is
	 * below its tolerance; single's 2e-04 is not, and one more step with the
	 * same Jacobian squares the relative error into the iterate's floor. A
	 * Jacobian that was NOT the map's slope would converge linearly, taking
	 * many more than two.
	 */
	meq::OuterNewtonResult const result =
		meq::solveForToroidalField( start, map );
	BOOST_TEST( result.converged,
		"Newton did not converge on an AFFINE map: " << result.status );
#ifdef MFEM_USE_SINGLE
	int const steps = 2;
#else
	int const steps = 1;
#endif
	BOOST_TEST( result.iterations <= steps,
		"Newton took " << result.iterations << " steps on an affine map, "
		"against " << steps << " in this format. More means the differenced "
		"Jacobian is not the map's own slope: " << result.status );

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

	// And a map that changes the length of its argument is a different space:
	// the loop ends there, reported as a map failure (see
	// a_map_that_throws_ends_the_loop_and_is_reported).
	auto ragged = []( std::vector<double> const &c )
	{
		return std::vector<double>( c.size() + 1, 0.0 );
	};
	meq::OuterNewtonResult const wrongSpace =
		meq::solveForToroidalField( start, ragged );
	BOOST_TEST( !wrongSpace.converged );
	BOOST_TEST( wrongSpace.status.find( "different number of coefficients" )
	            != std::string::npos, wrongSpace.status );
}

namespace
{
	/*
	 * A DETERMINISTIC NOISE IN [ -1, 1 ], a hash of the bits of the point and
	 * the component. It is what a map built on a whole equilibrium solve looks
	 * like from outside: smooth on the scale of a Newton step, rough on the
	 * scale of round-off, and the same number every time it is asked at the
	 * same point.
	 */
	double hashedNoise( std::vector<double> const &c, std::size_t component )
	{
		std::uint64_t h = 1469598103934665603ull;
		for ( double value : c )
		{
			std::uint64_t bits = 0;
			std::memcpy( &bits, &value, sizeof bits );
			h = ( h ^ bits )*1099511628211ull;
		}
		h ^= ( component + 1 )*0x9E3779B97F4A7C15ull;
		h ^= h >> 33;
		h *= 0xff51afd7ed558ccdull;
		h ^= h >> 33;
		return static_cast<double>( h >> 11 )
		       /static_cast<double>( 1ull << 53 )*2.0 - 1.0;
	}

	/// The unstable affine map of newton_finishes_an_unstable_affine_map_in_
	/// one_step, made non-linear so that Newton takes several steps, and noisy
	/// at @a noise.
	std::function<std::vector<double>( std::vector<double> const & )>
	noisyMap( double noise )
	{
		return [ noise ]( std::vector<double> const &c )
		{
			std::vector<double> const root = { 4.84, 2.42, 0.3025 };
			std::vector<double> out( c.size() );
			for ( std::size_t i = 0; i < c.size(); ++i )
			{
				double const away = c[ i ] - root[ i ];
				out[ i ] = root[ i ] + 1.8*away + 0.3*away*away
				           + noise*hashedNoise( c, i );
			}
			return out;
		};
	}
}

/*
 * A MAP THAT THROWS ENDS THE LOOP, AND THE LOOP SAYS SO -- WHEREVER IN THE
 * ITERATION IT THROWS.
 *
 * Most map evaluations happen inside KINSOL, which is C, so an exception there
 * must not unwind through its frames: each of OuterKINSolver's callbacks
 * catches, hands KINSOL its unrecoverable flag, and the result is a loop that
 * did not converge, with the map's message in its status. This is the path the
 * driver turns into ToroidalLoopFailed and exit 2, and no valid run reaches
 * it, so it is posed here with a map that throws on a chosen evaluation.
 *
 * On the noisy map's smooth form with three coefficients the evaluations go:
 * 1 the start, 2 to 7 the Jacobian there (both ahead of KINSOL), 8 the line
 * search's first trial -- KINSOL's RESIDUAL callback -- and from 9 the
 * Jacobian at the first accepted iterate, KINSOL's SETUP callback. Each of
 * the three places is asserted, and that the map is not called again after
 * it threw: KINSOL returns at once on the flag.
 */
BOOST_AUTO_TEST_CASE( a_map_that_throws_ends_the_loop_and_is_reported )
{
	std::vector<double> const start = { 6.0, 3.4, 0.5 };
	auto smooth = noisyMap( 0.0 );
	std::vector<double> const startImage = smooth( start );
	double startResidual = 0.0;
	for ( std::size_t i = 0; i < start.size(); ++i )
		startResidual = std::max( startResidual,
		                          std::abs( startImage[ i ] - start[ i ] ) );

	struct Where
	{
		int throwAt;
		char const *callback;
		bool ragged;
	};
	for ( Where const where : { Where{ 1, "starting coefficients", false },
	                            Where{ 8, "the residual", false },
	                            Where{ 10, "the differenced outer Jacobian",
	                                   false },
	                            Where{ 8, "the residual", true } } )
	{
		int calls = 0;
		std::string const message =
			"the inner solve failed at evaluation "
			+ std::to_string( where.throwAt );
		auto failing = [ & ]( std::vector<double> const &c )
		{
			if ( ++calls == where.throwAt )
			{
				if ( where.ragged )
					return std::vector<double>( c.size() + 1, 0.0 );
				throw std::runtime_error( message );
			}
			return smooth( c );
		};

		meq::OuterNewtonResult result;
		BOOST_REQUIRE_NO_THROW( result =
			meq::solveForToroidalField( start, failing ) );
		std::printf( "\n  A MAP FAILING AT EVALUATION %d%s: %s\n",
		             where.throwAt, where.ragged ? " (wrong length)" : "",
		             result.status.c_str() );

		BOOST_TEST( !result.converged,
			"a loop whose map threw reads as converged: " << result.status );
		BOOST_TEST( result.status.find( where.ragged
			? std::string( "different number of coefficients" ) : message )
			!= std::string::npos,
			"the status does not carry the map's own message: "
			<< result.status );
		BOOST_TEST( result.status.find( where.callback ) != std::string::npos,
			"the status does not say where the map failed (" << where.callback
			<< "): " << result.status );
		BOOST_TEST( calls == where.throwAt,
			"the map was called " << calls << " times, past the evaluation "
			<< where.throwAt << " that failed: KINSOL went on after the "
			"failure flag" );
		BOOST_TEST( result.mapEvaluations == where.throwAt );

		// UP TO THE FIRST ACCEPTED STEP THE LOOP IS STILL AT ITS START, and
		// reports the start's residual -- or NaN, where even that was never
		// evaluated. After it, the accepted iterate and its smaller residual.
		// The start is held in real_t, so the residual is compared to 1e-5.
		// The count is KINSOL's, which includes the iteration it failed in.
		int const failedIn = where.throwAt == 1 ? 0 : where.throwAt <= 8 ? 1 : 2;
		BOOST_TEST( result.iterations == failedIn );
		if ( where.throwAt <= 8 )
		{
			for ( std::size_t i = 0; i < start.size(); ++i )
				BOOST_TEST( result.coefficients[ i ]
				            == static_cast<double>(
				                   static_cast<mfem::real_t>( start[ i ] ) ) );
		}
		if ( where.throwAt == 1 )
			BOOST_TEST( std::isnan( result.finalResidual ) );
		else if ( where.throwAt == 8 )
			BOOST_TEST( std::abs( result.finalResidual - startResidual )
			            <= 1.0e-5*startResidual,
				"a loop that failed before its first step reports a residual "
				"of " << result.finalResidual << " against the start's "
				<< startResidual );
		else
		{
			BOOST_TEST( result.finalResidual < startResidual,
				"a loop that failed after one accepted step reports a residual "
				"of " << result.finalResidual << ", not below the start's "
				<< startResidual );
		}
	}
}

/*
 * A LOOP THAT STALLS ABOVE ITS MAP'S FLOOR IS NOT CONVERGED, AND ONE THAT
 * REACHES THE FLOOR IS -- WHATEVER KINSOL'S OWN FLAG SAYS ABOUT EITHER.
 *
 * **THE RULE: CONVERGED MEANS `|| G( c ) - c ||_inf` AT THE ANSWER IS AT OR
 * BELOW max( TOLERANCE, FLOOR ).** `mfem::KINSolver::GetConverged()` is
 * `flag >= 0`, which counts `KIN_STEP_LT_STPTOL` -- KINSOL's "the step got too
 * small to move", its report that it can make no more progress -- as
 * converged at any residual at all. And MFEM's default absolute tolerance,
 * `eps^( 1/3 )`, is 4.9e-03 in single, so the q-driven loop on
 * `examples/q-driven.toml` stopped there on `KIN_SUCCESS` and read "converged
 * ... final norm 0.00502" beside a requested 1e-9, with g^2 1-2% out.
 *
 * Three maps, each asserted both ways round, because a rule that only ever
 * says one thing passes half of this:
 *
 *   A STEP-TOLERANCE STALL, made deterministic: a smooth map, a step tolerance
 *   of 1e-2, so KINSOL stops on `KIN_STEP_LT_STPTOL` at a residual well above
 *   the 1e-12 asked for. Not converged -- and converged once the map declares
 *   a floor above where it stopped, which is the control.
 *
 *   A NOISY MAP WITH ITS FLOOR STATED: noise 1e-4 and a floor of 1e-3. It
 *   reaches the floor, in a handful of Newton steps.
 *
 *   THE SAME MAP WITH ITS FLOOR UNDERSTATED a hundredfold: no residual it can
 *   reach is below the threshold, so whatever KINSOL stops on, it is not
 *   converged -- and it must still stop, in a bounded number of evaluations.
 */
BOOST_AUTO_TEST_CASE( a_stall_above_the_floor_is_not_converged )
{
	std::vector<double> const start = { 6.0, 3.4, 0.5 };

	// ---- a deterministic step-tolerance stop -----------------------------
	meq::OuterNewtonOptions stalling;
	stalling.functionTolerance = 1.0e-12;
	stalling.stepTolerance = 1.0e-2;
	meq::OuterNewtonResult const stalled =
		meq::solveForToroidalField( start, noisyMap( 0.0 ), stalling );
	std::printf( "\n  A STEP-TOLERANCE STOP: %s\n", stalled.status.c_str() );
	BOOST_TEST( !stalled.converged,
		"a loop that stopped on its step tolerance at a residual of "
		<< stalled.finalResidual << ", above the 1e-12 asked for and with no "
		"floor declared, reads as converged: " << stalled.status );
	BOOST_TEST( stalled.finalResidual > stalling.functionTolerance,
		"the precondition: KINSOL must have stopped ABOVE the tolerance, or "
		"this case is not a stall" );

	meq::OuterNewtonOptions excused = stalling;
	double const declared = 2.0*stalled.finalResidual;
	excused.mapFloor = [ declared ]() { return declared; };
	meq::OuterNewtonResult const atFloor =
		meq::solveForToroidalField( start, noisyMap( 0.0 ), excused );
	BOOST_TEST( atFloor.converged,
		"the same stop with a floor declared at " << declared << " above it "
		"is not converged, so the floor is not part of the rule: "
		<< atFloor.status );
	BOOST_TEST( atFloor.floor >= declared );

	// ---- a noisy map, its floor stated -----------------------------------
	double const noise = 1.0e-4;
	meq::OuterNewtonOptions honest;
	honest.mapFloor = [ noise ]() { return 10.0*noise; };
	meq::OuterNewtonResult const reached =
		meq::solveForToroidalField( start, noisyMap( noise ), honest );
	std::printf( "  A NOISY MAP, FLOOR STATED: %s\n", reached.status.c_str() );
	BOOST_TEST( reached.converged, reached.status );
	BOOST_TEST( reached.finalResidual <= 10.0*noise );
	// A NEWTON, NOT A CRAWL: a fresh Jacobian at every iterate takes this
	// quadratic map from 1.2 away to its floor in a handful of steps, where
	// KINSOL's default -- a Jacobian kept for ten iterations -- contracts
	// linearly.
	BOOST_TEST( reached.iterations <= 6,
		"the loop took " << reached.iterations << " iterations to reach a "
		"floor of " << 10.0*noise << " on a quadratic map" );

	// ---- the same map, its floor understated -----------------------------
	meq::OuterNewtonOptions understated;
	understated.functionTolerance = 1.0e-12;
	understated.mapFloor = [ noise ]() { return noise/100.0; };
	meq::OuterNewtonResult const lost =
		meq::solveForToroidalField( start, noisyMap( noise ), understated );
	std::printf( "  A NOISY MAP, FLOOR UNDERSTATED: %s\n", lost.status.c_str() );
	BOOST_TEST( !lost.converged,
		"a loop whose residual " << lost.finalResidual << " sits above the "
		"threshold " << noise/100.0 << " reads as converged: " << lost.status );
	BOOST_TEST( lost.finalResidual > noise/100.0 );
	BOOST_TEST( lost.mapEvaluations
	            <= ( understated.maxIterations + 1 )
	               *( 2*static_cast<int>( start.size() ) + 6 ) );
}

/*
 * THE OUTER NEWTON IS A NEWTON: A FRESH JACOBIAN AT EVERY ITERATE.
 *
 * KINSOL's default keeps a Jacobian for ten iterations, which on the q-driven
 * loop meant the vacuum start's for the whole run: measured on
 * `examples/q-driven.toml` in single, the max norm went 0.244, 0.191, 0.149,
 * ..., 0.008, 0.0042 over eight iterations -- a linear crawl, which this tree
 * treats as a defect rather than as convergence. This map makes the same
 * mistake visible: cubic, so the slope at the start is eleven times the
 * slope at the root, and a Jacobian kept from the start takes steps eleven
 * times too short. Newton's own count on it is asserted, with the cost.
 */
BOOST_AUTO_TEST_CASE( the_outer_newton_forms_a_fresh_jacobian_every_iteration )
{
	std::vector<double> const root = { 4.84, 2.42, 0.3025 };
	auto cubic = [ &root ]( std::vector<double> const &c )
	{
		std::vector<double> out( c.size() );
		for ( std::size_t i = 0; i < c.size(); ++i )
		{
			double const away = c[ i ] - root[ i ];
			out[ i ] = root[ i ] + 1.8*away + 2.0*away*away*away;
		}
		return out;
	};

	std::vector<double> const start = { 6.0, 3.4, 0.5 };
	meq::OuterNewtonOptions options;
	options.functionTolerance = 1.0e-6;
	meq::OuterNewtonResult const result =
		meq::solveForToroidalField( start, cubic, options );
	std::printf( "\n  A CUBIC MAP: %s\n", result.status.c_str() );
	BOOST_TEST( result.converged, result.status );
	// Measured: 6 iterations and 43 evaluations, one residual and one fresh
	// Jacobian each. With KINSOL's default the same run takes 9 and 64 -- MORE
	// evaluations, not fewer, because a kept Jacobian's short steps are what
	// its line search then spends equilibria on. One iteration of margin.
	BOOST_TEST( result.iterations <= 7,
		"the outer Newton took " << result.iterations << " iterations on a "
		"cubic map that Newton itself finishes in 6: a Jacobian is being kept "
		"from an earlier iterate. " << result.status );
	BOOST_TEST( result.mapEvaluations
	            <= 1 + 7*( 2*static_cast<int>( start.size() ) + 1 ) );
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
	// A SAFETY FACTOR IS A FLUX-SURFACE AVERAGE, so its floor is not eps but
	// eps/h^2 -- M-220's fourth law, and the sharpest case of it in the tree.
	// At n = 24 that is 6.9e-05 and the residual reads 5.78e-05 in single,
	// which is the average's own round-off and not a root that has moved. In
	// double averagedMetricFloor() is 1.3e-13 and the 1e-8 stands.
	double const averageFloor =
		std::max( 1.0e-8, meq::tests::averagedMetricFloor( 1.0/n ) );
	BOOST_TEST( worst < averageFloor,
		"the residual at the closed-form answer is " << worst
		<< ", so the coefficients the outer solver iterates in do not have the "
		"answer as a root and nothing it finds could be right" );
}
