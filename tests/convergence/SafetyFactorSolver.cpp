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
				return 0.30*std::sin( M_PI*( x( 0 ) - box.rMin )/box.width() )
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
 * ==> THIS CASE IS RED ON PURPOSE AND NAMES WHAT WOULD FIX IT. <==
 *
 * Per *Testing stance*: a defect gets a FAILING test naming it, never a passing
 * test recording it, and the suite is expected to be red while a known defect
 * stands. What stands here is that the damped Picard outer loop does not
 * converge, and the diagnosis is below rather than in a commit message because
 * the next person to touch this needs it.
 *
 * WHAT IS ESTABLISHED. theReferenceIsAFixedPointOfItsOwnTarget shows the
 * answer IS a fixed point of this map, to **2.1e-14** -- inverting the
 * reference's own q against the reference's own geometry returns the g it was
 * built with. So the target is consistent with the inversion and the loop is
 * not chasing the wrong point.
 *
 * WHAT FAILS. Started 40% away and stepped adaptively, the loop walks the error
 * down 0.362 -> 0.019 at iteration 7 -- essentially the answer -- and then
 * straight back up, 0.048, 0.094, 0.140, to 0.50 by iteration 29, with psi_ax
 * marching monotonically past the reference's 3.457e-01 and out the other side.
 * The step never shrinks: |dg/g| sits at 0.04 to 0.05 throughout, including
 * where the error passes through zero. **It does not stall at the fixed point,
 * it crosses it.**
 *
 * AND DAMPING CANNOT BE THE ANSWER, WHICH IS THE PART WORTH KEEPING. A relaxed
 * iteration x <- x + w( G( x ) - x ) has derivative 1 + w( G' - 1 ) at the fixed
 * point. For G' > 1 that is above one for EVERY w > 0: under-relaxation
 * stabilises a map that oscillates ( G' < -1 ) and can do nothing at all for one
 * that runs away. So no step control, no relaxation and no continuation in w
 * will make this converge, and the measurements above are what say so rather
 * than an argument from the shape of the problem.
 *
 * WHAT WOULD FIX IT is an outer NEWTON or an Anderson acceleration on the g
 * table -- ROADMAP.md item 10's own "the non-local dependence inside dF/dpsi if
 * Newton is to stay quadratic". The Jacobian is d( V', < R^-2 > )/dg, which is
 * small: g is a handful of coefficients once fitted, not nFieldDOF, so
 * differencing it costs one extraction per COEFFICIENT rather than
 * INVERSION-PLAN.md section 11.1's 5.7 hours per field degree of freedom. At
 * degree two that is three extra solves an outer step. **That is the next
 * thing to build, and it is affordable.**
 *
 *
 * THE LOOP CLOSES ON THE FIELD IT WAS BUILT FROM.
 *
 * Reference: g = g0 + g1 Psi in closed form. Its q on the family's own grid is
 * V' g < R^-2 >/4 pi^2 with g EVALUATED rather than extracted, so the target is
 * exact to the geometry's own accuracy and no further.
 *
 * The loop then starts from a g that is 40% too large everywhere -- a different
 * equilibrium, not a perturbation -- and has to walk back. What is asserted is
 * that it converges, that it converges TO THE REFERENCE, and that the residual
 * falls monotonically enough to be a contraction rather than a wander.
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
	std::printf( "    %4s %14s %14s %14s   %s\n",
	             "it", "max |dg/g|", "worst g err", "psi_ax", "gg' range" );

	/*
	 * ---- the loop, started 40% away, UNDER-RELAXED ------------------------
	 *
	 * UNDAMPED THIS DIVERGES, AND THAT IS THE FINDING RATHER THAN AN
	 * INCONVENIENCE. Taking the inverted g whole gives 0.324 -> 0.277 -> 0.364
	 * in the distance from the reference and then a solve that does not
	 * converge at all: the map g -> ( V', < R^-2 > ) -> g has gain above one
	 * here, because a larger g raises gg', which raises F, which moves the
	 * surfaces the next inversion divides by. That is the same shape as this
	 * tree's other Picard result -- undamped stalls, omega = 0.5 reaches
	 * 2.8e-08 -- and it is why the loop STATE is the g table rather than the
	 * profile: the relaxation has to happen on g, where it is a physical
	 * quantity, and not on the knots.
	 */
	double const relaxation = 0.5;
	double const startScale = 1.40;

	std::vector<double> gTable;
	for ( double x : label )
		gTable.push_back( startScale*referenceG( 1.0 - x ) );

	auto fieldOf = [ & ]( std::vector<double> const &g )
	{
		meq::ToroidalField field;
		field.normalisedFlux = label;
		field.g = g;
		for ( double value : g )
		{
			field.gSquared.push_back( value*value );
			field.safetyFactor.push_back( 0.0 );
		}
		return field;
	};

	std::vector<double> change, worstError;
	double lastAxis = 0.0;
	double step = relaxation;
	int retreats = 0;

	std::unique_ptr<Solved> current = std::make_unique<Solved>( solveWith(
		mesh, std::make_shared<meq::SplineProfile const>(
			meq::ggPrimeKnotsFromFit( fieldOf( gTable ), 2 ) ),
		order, nullptr, 0.30 ) );
	BOOST_TEST_REQUIRE( current->converged,
		"the loop's own starting profile does not solve, so there is nothing to "
		"iterate from" );
	lastAxis = current->psiAxis;

	for ( int iteration = 0; iteration < 30; ++iteration )
	{
		meq::FluxFamilyOptions loopOptions = familyOptions();
		meq::FluxSurfaceFamily const family =
			meq::extractFluxSurfaces( *current->solver, loopOptions );
		meq::ToroidalField const inverted =
			meq::invertSafetyFactor( family, target );
		BOOST_TEST_REQUIRE( inverted.size() == gTable.size(),
			"the family's grid moved between iterations" );

		/*
		 * ADAPTIVE RELAXATION: halve on a failed solve, grow by 1.3 on a
		 * successful one. The same step control section 4.4's continuation
		 * uses, and for the same reason -- there, uniform steps stall at
		 * c3 = -10.8 while the adaptive one walks the whole way in nine solves
		 * with two retreats. Here the thing being protected is the INNER
		 * bordered Newton, which has no globalisation of its own.
		 */
		std::vector<double> trial( gTable.size() );
		std::unique_ptr<Solved> attempt;
		double moved = 0.0;
		bool accepted = false;

		for ( int retry = 0; retry < 8 && !accepted; ++retry )
		{
			moved = 0.0;
			for ( std::size_t i = 0; i < gTable.size(); ++i )
			{
				trial[ i ] = gTable[ i ]
				             + step*( inverted.g[ i ] - gTable[ i ] );
				moved = std::max( moved, std::abs( trial[ i ] - gTable[ i ] )
				                         /std::abs( trial[ i ] ) );
			}

			attempt = std::make_unique<Solved>( solveWith(
				mesh, std::make_shared<meq::SplineProfile const>(
					meq::ggPrimeKnotsFromFit( fieldOf( trial ), 2 ) ),
				order, &current->solver->potential(), lastAxis ) );

			if ( attempt->converged )
				accepted = true;
			else
			{
				step *= 0.5;
				++retreats;
			}
		}

		BOOST_TEST_REQUIRE( accepted,
			"the inner solve failed at every relaxation down to " << step
			<< ". The bordered Newton has no globalisation available, so the "
			"outer step is the only control there is." );

		gTable = trial;
		current = std::move( attempt );
		lastAxis = current->psiAxis;
		change.push_back( moved );

		double worst = 0.0;
		for ( std::size_t i = 0; i < gTable.size(); ++i )
		{
			double const exact = referenceG( 1.0 - label[ i ] );
			worst = std::max( worst,
			                  std::abs( gTable[ i ] - exact )/std::abs( exact ) );
		}
		worstError.push_back( worst );

		std::printf( "    %4d %14.6e %14.6e %14.6e   step %.4f\n",
		             iteration, moved, worst, lastAxis, step );
		std::fflush( stdout );

		step = std::min( 1.0, 1.3*step );
		if ( moved < 1.0e-9 )
			break;
	}
	std::printf( "    %d retreats\n", retreats );

	BOOST_TEST_REQUIRE( change.size() >= 3u,
		"the loop did not run long enough to say anything about contraction" );

	// IT CONTRACTS. Not a rate -- a Picard iteration on the geometry has no
	// order to claim -- but the step must be shrinking by a real factor, or
	// this is a wander that happened to stop.
	double const contraction = change.back()/change.front();
	std::printf( "    step fell by %.3e over %zu iterations\n",
	             contraction, change.size() );
	BOOST_TEST( contraction < 1.0e-2,
		"the loop's step fell only by " << contraction << " over "
		<< change.size() << " iterations, which is a wander rather than a "
		"contraction" );

	// AND IT CONVERGES TO THE REFERENCE, WHICH IS THE ASSERTION WITH TEETH. A
	// loop that reached any fixed point at all would pass the line above.
	std::printf( "    worst |g - g_exact|/g: %.6e -> %.6e\n",
	             worstError.front(), worstError.back() );
	BOOST_TEST( worstError.back() < 0.02,
		"the loop settled with g still " << worstError.back()
		<< " from the closed form its own target was built from, having started "
		<< worstError.front() << " away. It reached A fixed point and not THE "
		"one." );
	BOOST_TEST( worstError.back() < 0.2*worstError.front(),
		"the loop did not move g materially closer to the reference: "
		<< worstError.front() << " -> " << worstError.back() );

	// The equilibrium comes back too, which is the statement a consumer cares
	// about -- q is an input now, and psi_ax is one of the things it decides.
	std::printf( "    psi_ax %.6e against the reference's %.6e (%.2e)\n",
	             lastAxis, referenceAxis,
	             std::abs( lastAxis - referenceAxis )/referenceAxis );
	BOOST_TEST( std::abs( lastAxis - referenceAxis ) < 0.02*referenceAxis,
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
